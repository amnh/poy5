(* POY 5.1.1. A phylogenetic analysis program using Dynamic Homologies.       *)
(* Copyright (C) 2014 Andrés Varón, Lin Hong, Nicholas Lucaroni, Ward Wheeler,*)
(* and the American Museum of Natural History.                                *)
(*                                                                            *)
(* This program is free software; you can redistribute it and/or modify       *)
(* it under the terms of the GNU General Public License as published by       *)
(* the Free Software Foundation; either version 2 of the License, or          *)
(* (at your option) any later version.                                        *)
(*                                                                            *)
(* This program is distributed in the hope that it will be useful,            *)
(* but WITHOUT ANY WARRANTY; without even the implied warranty of             *)
(* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the              *)
(* GNU General Public License for more details.                               *)
(*                                                                            *)
(* You should have received a copy of the GNU General Public License          *)
(* along with this program; if not, write to the Free Software                *)
(* Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301   *)
(* USA                                                                        *)

let () = SadmanOutput.register "MlModel" "$Revision: 3653 $"

open Numerical.FPInfix

let (-->) b a = a b

let failwithf format = Printf.ksprintf failwith format

let debug = false

exception LikelihoodModelError of string

let likelihood_not_enabled =
    "Likelihood not enabled: download different binary or contact mailing list" 

let lfailwith () =
    raise (LikelihoodModelError likelihood_not_enabled)

let dyno_likelihood_warning = 
    "Gap as an additional character is required for the dynamic "^
    "likelihood criteria. Please add argument gap:coupled, or gap:character."

let dyno_gamma_warning = 
    "Gamma classes for dynamic MPL are un-necessary, and are being removed."

let debug_printf msg format =
    Printf.ksprintf (fun x -> if debug then print_string x; flush stdout) msg format

let ba_of_array1 x = Bigarray.Array1.of_array Bigarray.float64 Bigarray.c_layout x
and ba_of_array2 x = Bigarray.Array2.of_array Bigarray.float64 Bigarray.c_layout x

let create_ba1 x   = Bigarray.Array1.create Bigarray.float64 Bigarray.c_layout x
and create_ba2 x y = Bigarray.Array2.create Bigarray.float64 Bigarray.c_layout x y

let barray_to_array2 bray =
    let a = Bigarray.Array2.dim1 bray and b = Bigarray.Array2.dim2 bray in
    let r = Array.make_matrix a b 0.0 in
    for i = 0 to a-1 do for j = 0 to b-1 do
        r.(i).(j) <- bray.{i,j};
    done; done; r

let print_barray1 a =
    for i = 0 to (Bigarray.Array1.dim a)-1 do
        Printf.printf "%2.10f\t" a.{i};
    done; Printf.printf "\n%!"; ()

and print_barray2 a =
    for i = 0 to (Bigarray.Array2.dim1 a)-1 do
        for j = 0 to (Bigarray.Array2.dim2 a)-1 do
            Printf.printf "%2.10f\t" a.{i,j};
        done; Printf.printf "\n";
    done; Printf.printf "\n%!"; ()

(* type to help the parsing of specification data *)
type string_spec = string * (string * string * string * string) * float list
                 * [ `Given of (string * float) list | `Estimate of float array option
                   | `Consistent of float array option | `Equal ]
                 * (string * float option) * string * string option
let empty_str_spec : string_spec = ("",("","","",""),[],`Consistent None,("",None),"",None)

(* the default when commands are left off from the interactive console / scripts *)
let default_command = (`Max, `MAL,`GTR [],None,`Consistent,`Missing)

(**  [default_tstv] default ratio for models wth tstv rates; same as phyml *)
let default_tstv  = 4.0

(** [default_gtr] the default parameter rates for GTR; all 1's. False for
    coupled. otherwise, true *)
let default_gtr a =
    function | false -> Array.make (((a-2)*(a+1))/2) 1.0
             | true  -> Array.make (((a)*(a-3))/2)   1.0 (* Coupled *)

(** [default_alpha] default value for alpha parameter; these are the same as
    phyml. p is if theta is being used (as in is a percentage) *)
let default_alpha p = if p then 0.2 else 1.0

(** [default_invar] default percent for invariant sites *)
let default_invar   = 0.2

(** [default_gap_r] default gap rate, for coupled gap parameters *)
let default_gap_r   = 0.15

(** Define how characters can be sent to the classification function *)
type chars = [ `List of int list | `Packed of int ]

(** Define a distribution that branch lengths can be part of *)
type site_var = 
    | Gamma of int * float
    | Theta of int * float * float
    (* | Given of (float * float) array *)
    | Constant 

(** Define the substitution rate model, between characters *)
type subst_model =
    | JC69
    | F81
    | K2P    of float
    | F84    of float
    | HKY85  of float
    | TN93   of (float * float)
    | GTR    of float array
    | File   of float array array * string
    | Custom of (int All_sets.IntegerMap.t * float array * string)

(** Define how base frequencies are calculated and used *)
type priors = 
    | Estimated of float array
    | Given     of float array
    | Equal

(** This is a specification of a model; and provides a one-to-one mapping to a
    full model defined below **)
type spec = {
    substitution : subst_model;
    site_variation : site_var;
    base_priors : priors;
    cost_fn : Methods.ml_costfn;
    use_gap : Methods.ml_gap;
    alphabet : Alphabet.a * int;
}

(** Define the Model; here we store a diagonalized matrix, priors in C-format
    big-array, along with some other minor details *)
type model = {
    spec  : spec;
    pi_0  : (float, Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array1.t;
    invar : float option;
    rate  : (float, Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array1.t;
    prob  : (float, Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array1.t;
    s     : (float, Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array2.t;
    u     : (float, Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array2.t;
    d     : (float, Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array2.t;
    ui    : (float, Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array2.t option; 
}

(** This is a mapping to the codes on the C-side for determining median
    calculation functions from the cost-style *)
let get_costfn_code a = match a.spec.cost_fn with 
    | `MPL -> 1 
    | `MAL -> 0 

(** Return the alphabet for the characters *)
let get_alphabet m = fst m.spec.alphabet

(** Return the size of the alphabet for this model; includes gaps if necessary *)
let get_alphabet_size m = snd m.spec.alphabet

(** Compare two models; not to be used for a total ordering (ret neg or zero) *)
let compare a b = 
    let compare_array x y =
        let results = ref true in
        for i = 0 to (snd a.spec.alphabet) - 1 do
            results := !results && (x.(i) =. y.(i));
        done;
        !results
    in
    let compare_priors a b = match a.spec.base_priors ,a.spec.base_priors with
            | _ when (snd a.spec.alphabet) != (snd b.spec.alphabet) -> false
            | Equal , Equal            -> true
            | Estimated x, Estimated y
            | Given     x, Given     y -> compare_array x y
            | (Equal | Estimated _ | Given _), _ -> false
    in
    let m_compare = match a.spec.substitution,b.spec.substitution with
        | JC69 , JC69 | F81 , F81 -> 0
        | K2P x, K2P y | F84 x, F84 y | HKY85 x, HKY85 y when x = y -> 0
        | Custom _, Custom _ -> 0
        | TN93 (x1,x2),  TN93 (y1,y2) when x1 = y1 && x2 = y2 -> 0
        | GTR xs, GTR ys when compare_array xs ys -> 0
        | File (_,x), File (_,y) when x = y -> 0
        | _,_ -> ~-1
    and c_compare = if a.spec.cost_fn = b.spec.cost_fn then 0 else ~-1
    and v_compare = match a.spec.site_variation,b.spec.site_variation with
        | Gamma (ix,ax), Gamma (iy,ay) ->
            if ix=iy && ax=ay then 0 else ~-1
        | Theta (ix,ax,bx), Theta (iy,ay,by) ->
            if ix=iy && ax=ay && bx=by then 0 else ~-1
        | Constant, Constant            -> 0
        | (Gamma _|Theta _|Constant), _ -> ~-1
    and g_compare = match a.spec.use_gap,b.spec.use_gap with
        | `Missing, `Missing
        | `Independent, `Independent -> 0
        | `Coupled x, `Coupled y when x=y-> 0
        | (`Missing | `Independent | `Coupled _), _ -> ~-1
    and p_compare = if compare_priors a b then 0 else ~-1 in
    (* just knowing that they are different is enough *)
    if debug then begin
        Printf.printf "Models : %b\n%!" (m_compare = 0);
        Printf.printf "Rate   : %b\n%!" (v_compare = 0);
        Printf.printf "Use Gap: %b\n%!" (g_compare = 0);
        Printf.printf "Priors : %b\n%!" (p_compare = 0);
        Printf.printf "Cost   : %b\n%!" (c_compare = 0);
    end;
    m_compare + v_compare + g_compare + p_compare + c_compare
    

module OrderedML = struct
    (* we could choose model or spec, but spec can use Pervasives.compare *)
    type t = spec 
    let compare a b = Pervasives.compare a b
end
module MlModelMap = Map.Make (OrderedML)

(** Categorize a list of codes by model *)
let categorize_by_model get_fn codes =
    let set_codes,non_lk_codes =
        List.fold_left
            (fun (acc,oth) code ->
                try let spec = get_fn code in
                    try let old = MlModelMap.find spec acc in
                        (MlModelMap.add spec (code::old) acc,oth)
                    with | Not_found ->
                        (MlModelMap.add spec ([code]) acc,oth)
                with | _ -> (acc,code::oth))
            (MlModelMap.empty,[])
            codes
    in
    let init = match non_lk_codes with | [] -> [] | xs -> [xs] in
    MlModelMap.fold (fun _ e a -> e :: a) set_codes init

(** Count the number of parameters in the model; used for xIC functions **)
let count_parameters model : int =
    let num_subst = match model.spec.substitution with
        | Custom (_,x,_)          -> (Array.length x) - 1
        | JC69  | F81 | File _    -> 0
        | K2P _ | F84 _ | HKY85 _ -> 1
        | TN93 _                  -> 2
        | GTR _  ->
            let a = snd model.spec.alphabet in
            ((a*(a-1))/2)-1
    and num_rates = match model.spec.site_variation with
        | Constant        -> 0
        | Gamma _         -> 1
        | Theta (x,_,_) when x = 1 -> 1
        | Theta _         -> 2
    and num_priors = match model.spec.base_priors with
        | Estimated f -> (Array.length f) - 1
        | Given     _ -> 0
        | Equal       -> 0
    in
    num_subst + num_priors + num_rates

(** Return a list of all the models with default parameters; requires alphabet
    to determine if tstv models are valid and length of GTR parameters *)
let get_all_models alph gap est_prior =
    [ (JC69,Equal);
      (F81,est_prior);
      (K2P default_tstv,Equal);
      (F84 default_tstv,est_prior);
      (HKY85 default_tstv,est_prior);
      (TN93 (default_tstv,default_tstv),est_prior);
      (GTR [||],est_prior)]


IFDEF USE_LIKELIHOOD THEN

(* ------------------------------------------------ *)
(* EXTERNAL FUNCTIONS -- maintained in likelihood.c *)

(* diagonalize a symmetric or gtr matrix, WARNING: modifies passed matrices *)
external diagonalize_gtr: (* U D Ui *) FMatrix.m ->
    (float,Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array2.t ->
    (float,Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array2.t ->
    (float,Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array2.t ->
    unit = "likelihood_CAML_diagonalize_gtr"
external diagonalize_sym: (* U D *) FMatrix.m ->
    (float,Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array2.t ->
    (float,Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array2.t ->
    unit = "likelihood_CAML_diagonalize_sym"

(* compose matrices -- for testing purposes, as this composition is
 * usually done on the C side exclusively. If the time is less then zero, the
 * instantaneious rate matrix will be returned instead (which is just UDUi). *)
external compose_gtr: (* U D Ui t *) FMatrix.m ->
    (float,Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array2.t -> 
    (float,Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array2.t ->
    (float,Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array2.t -> float
    -> (float,Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array2.t =
        "likelihood_CAML_compose_gtr"
external compose_sym: (* U D *) FMatrix.m ->
    (float,Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array2.t ->
    (float,Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array2.t -> float
    -> (float,Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array2.t =
        "likelihood_CAML_compose_sym"

(* ------------------------------------------------ *)
(* MODEL CALCULATION FUNCTIONS                      *)

(* divide a matrix by the mean rate so it will equal 1 *)
let m_meanrate srm pi_ =
    let mr = ref 0.0 and a_size = Bigarray.Array2.dim1 srm in
    for i = 0 to (a_size-1) do
        mr := !mr +. (~-.(srm.{i,i}) *. pi_.{i});
    done;
    for i = 0 to (a_size-1) do
        for j = 0 to (a_size-1) do
            srm.{i,j} <- srm.{i,j} /. !mr;
        done;
    done

(* val jc69 :: ANY ALPHABET size *)
let m_jc69 pi_ mu a_size gap_r =
    let srm = create_ba2 a_size a_size in
    Bigarray.Array2.fill srm mu;
    let () = match gap_r with
        | None -> 
            let diag = -. mu *. float (a_size-1) in
            for i = 0 to (a_size-1) do
                srm.{i,i} <- diag
            done;
        | Some (i,r) ->
            let diag = -. mu *. (r +. float (a_size-2)) in
            for j = 0 to (a_size-1) do
                srm.{i,j} <- mu *. r;
                srm.{j,i} <- mu *. r;
                srm.{j,j} <- diag;
            done;
            srm.{i,i} <- -. mu *. r *. float (a_size-1)
    in
    (* normalize by mean-rate *)
    let mr = ref 0.0 and wght = 1.0 /. float a_size in
    for i = 0 to (a_size-1) do
        mr := !mr +. (~-.(srm.{i,i}) *. wght );
    done;
    for i = 0 to (a_size-1) do
        for j = 0 to (a_size-1) do
            srm.{i,j} <- srm.{i,j} /. !mr;
        done;
    done;
    srm

(* val k2p :: only 4 or 5 characters *)
let m_k2p pi_ alpha beta a_size gap_r =
    if not ((a_size = 4) || (a_size = 5)) then begin
        raise (LikelihoodModelError "Alphabet does not support this model")
    end;
    let srm = create_ba2 a_size a_size in
    let beta = if beta >. 0.0 then beta else Numerical.minimum in
    Bigarray.Array2.fill srm beta;
    (* modify transition elements to alpha *)
    srm.{1, 3} <- alpha; srm.{3, 1} <- alpha;
    srm.{2, 0} <- alpha; srm.{0, 2} <- alpha;
    (* set up the diagonal elements *)
    let () = match gap_r with
        | None -> 
            let diag =
                if a_size = 4 then
                    -. alpha -. beta -. beta
                else begin
                    assert ( a_size = 5 );
                    -. alpha -. beta *. 3.0
                end
            in
            for i = 0 to (a_size-1) do
                srm.{i,i} <- diag
            done;
        | Some (i,r) ->
            assert( a_size = 5 );
            let diag = -. alpha -. beta *. 3.0 -. beta *. r in
            for j = 0 to (a_size-1) do
                srm.{i,j} <- beta *. r;
                srm.{j,i} <- beta *. r;
                srm.{j,j} <- diag;
            done;
            srm.{i,i} <- -. beta *. r *. float (a_size-1)
    in
    (* normalize by mean-rate *)
    let mr = ref 0.0 and wght = 1.0 /. float a_size in
    for i = 0 to (a_size-1) do
        mr := !mr +. (~-.(srm.{i,i}) *. wght );
    done;
    for i = 0 to (a_size-1) do
        for j = 0 to (a_size-1) do
            srm.{i,j} <- srm.{i,j} /. !mr;
        done;
    done;
    srm

(* val tn93 :: only 4 or 5 characters *)
let m_tn93 pi_ alpha beta gamma a_size gap_r =
    if not ((a_size = 4) || (a_size = 5)) then begin
        raise (LikelihoodModelError "Alphabet size does not support this model")
    end;
    let srm = create_ba2 a_size a_size in
    let gamma = if gamma >. 0.0 then gamma else Numerical.minimum in
    Bigarray.Array2.fill srm gamma;
    srm.{0,2} <- alpha; srm.{1,3} <- beta; (* ACGT -- R=AG -- Y=CT *)
    srm.{2,0} <- alpha; srm.{3,1} <- beta; (* 0123 -- R=02 -- Y=13 *)
    let () = match gap_r with
        | None -> 
            for i = 0 to (a_size-1) do
                for j = 0 to (a_size-1) do
                    srm.{i,j} <- srm.{i,j} *. pi_.{j};
                done;
            done;
        | Some (k,r) -> 
            for i = 0 to (a_size-1) do
                for j = 0 to (a_size-1) do
                    srm.{i,j} <- srm.{i,j} *. pi_.{j};
                done;
                srm.{i,k} <- gamma *. r *. pi_.{k};
                srm.{k,i} <- gamma *. r *. pi_.{i};
            done;
    in
    (* normalize diagonal so row sums to 0 *)
    for i = 0 to (a_size-1) do
        let diag = ref 0.0 in
        for j = 0 to (a_size-1) do
            if (i <> j) then diag := !diag +. srm.{i,j};
        done;
        srm.{i,i} <- -. !diag;
    done;
    m_meanrate srm pi_;
    srm

let m_tn93_ratio pi_ kappa1 kappa2 a_size gap_r =
    let beta = (pi_.{0} *. pi_.{2} *. kappa1) +. (pi_.{1} *. pi_.{3} *. kappa2) +.
                ((pi_.{0} +. pi_.{2}) *. (pi_.{1}+.pi_.{3})) in
    let beta = 1.0 /. (2.0 *. beta) in
    let alpha1 = kappa1 *. beta and alpha2 = kappa2 *. beta in
    m_tn93 pi_ alpha1 alpha2 beta a_size gap_r

(* val f81 :: ANY ALPHABET size *)
let m_f81 pi_ lambda a_size gap_r =
    let srm = create_ba2 a_size a_size in
    let lambda = if lambda >. 0.0 then lambda else Numerical.minimum in
    let () = match gap_r with
        | None -> 
            for i = 0 to (a_size-1) do
                for j = 0 to (a_size-1) do
                    if i = j then ()
                             else srm.{i,j} <- pi_.{j} *. lambda;
                done;
            done;
        | Some (k,r) ->
            for i = 0 to (a_size-1) do
                for j = 0 to (a_size-1) do
                         if i = j then ()
                    else if k = j then srm.{i,j} <- pi_.{j} *. r
                                  else srm.{i,j} <- pi_.{j} *. lambda;
                done;
            done;
    in
    for i = 0 to (a_size-1) do
        let diag = ref 0.0 in
        for j = 0 to (a_size-1) do
            if (i <> j) then diag := !diag +. srm.{i,j};
        done;
        srm.{i,i} <- -. !diag;
    done;
    m_meanrate srm pi_;
    srm

(* val hky85 :: only 4 or 5 characters *)
let m_hky85 pi_ kappa a_size gap_r =
    let beta = (pi_.{0} *. pi_.{2} *. kappa) +. (pi_.{1} *. pi_.{3} *. kappa) +.
                ((pi_.{0} +. pi_.{2}) *. (pi_.{1}+.pi_.{3}) ) in
    let beta = 1.0 /. (2.0 *. beta) in
    let alpha = kappa *. beta in
    m_tn93 pi_ alpha alpha beta a_size gap_r

(* val f84 :: only 4 or 5 characters *)
let m_f84 pi_ gamma kappa a_size gap_r =
    let gamma = if gamma >. 0.0 then gamma else Numerical.minimum in
    let y = pi_.{1} +. pi_.{3} in (* Y = C + T *)
    let r = pi_.{0} +. pi_.{2} in (* R = A + G *)
    let alpha = (1.0+.kappa/.r) *. gamma in
    let beta = (1.0+.kappa/.y) *. gamma in
    m_tn93 pi_ alpha beta gamma a_size gap_r

(* normalize against two characters that are not gaps; unless its the only choice *)
let normalize ?(m=Numerical.minimum) alph gap_state vec = match gap_state with
    | `Independent ->
        let normalize_factor = max m vec.((Array.length vec) - 1) in
        let vec = Array.map (fun i -> (max m i) /. normalize_factor) vec in
        vec, `Independent
    | `Coupled n ->
        let normalize_factor = max m vec.((Array.length vec) - 1) in
        let vec = Array.map (fun i -> (max m i) /. normalize_factor) vec in
        vec, `Coupled (n/.normalize_factor)
    | `Missing ->
        let normalize_factor = max m vec.((Array.length vec) - 1) in
        let vec = Array.map (fun i -> (max m i) /. normalize_factor) vec in
        vec,`Missing

(* val gtr :: ANY ALPHABET size
   pi_ == n; co_ == ((n-1)*n)/2
   form of lower packed storage mode, excluding diagonal, *)
let m_gtr_independent pi_ co_ a_size =
    if (((a_size+1)*(a_size-2))/2) <> Array.length co_ then begin
        failwithf ("Length of GTR parameters (I) is incorrect for the alphabet."
                   ^^"They should be %d, but are %d.")
                  (((a_size+1)*(a_size-2))/2) (Array.length co_);
    end;
    (* last element of GTR = 1.0 *)
    let co_ =
        let size = (((a_size-1)*a_size)/2) in
        Array.init (size)
                   (fun i -> if i = (size-1) then 1.0 else co_.(i))
    in
    (* create matrix *)
    let n = ref 0 in (* array index *)
    let srm = create_ba2 a_size a_size in
    for i = 0 to (a_size-1) do
        for j = (i+1) to (a_size-1) do
            srm.{i,j} <- co_.(!n) *. pi_.{j};
            srm.{j,i} <- co_.(!n) *. pi_.{i};
            incr n;
        done;
    done;
    (* set diagonal so row sums to 0 *)
    for i = 0 to (a_size-1) do begin
        let diag = ref 0.0 in
        for j = 0 to (a_size-1) do
            if (i <> j) then diag := !diag +. srm.{i,j};
        done;
        srm.{i,i} <- -. !diag;
    end; done;
    (* divide through by mean-rate *)
    m_meanrate srm pi_;
    srm

let m_gtr_coupled pi_ co_ a_size i_gap r_gap =
    if (((a_size)*(a_size-3))/2) <> Array.length co_ then begin
        failwithf ("Length of GTR parameters (%d) is incorrect for alphabet "^^
                   "(%d). It should be %d.")
                  (Array.length co_) a_size (((a_size)*(a_size-3))/2)
    end;
     (* last element of GTR = 1.0; add back *)
    let co_ =
        let size = Array.length co_ in
        Array.init (size+1)
                   (fun i -> if i = (size) then 1.0 else co_.(i))
    in
    (* create matrix *)
    let n = ref 0 in (* array index *)
    let srm = create_ba2 a_size a_size in
    for i = 0 to (a_size-1) do
        (* set the gap and gap coefficient col/row *)
        if i = i_gap then begin
            for j = 0 to (a_size-1) do
                srm.{j,i} <- r_gap *. pi_.{i};
                srm.{i,j} <- r_gap *. pi_.{j};
            done;
        end else begin
            for j = i+1 to (a_size-1) do
                if j = i_gap then ()
                else begin
                    srm.{i,j} <- co_.(!n) *. pi_.{j};
                    srm.{j,i} <- co_.(!n) *. pi_.{i};
                    incr n;
                end;
            done;
        end;
    done;
    (* set diagonal so row sums to 0 *)
    for i = 0 to (a_size-1) do begin
        let diag = ref 0.0 in
        for j = 0 to (a_size-1) do
            if (i <> j) then diag := !diag +. srm.{i,j};
        done;
        srm.{i,i} <- -. !diag;
    end; done;
    (* divide through by mean-rate *)
    m_meanrate srm pi_;
    srm

let m_gtr pi_ co_ a_size gap_r =
    let srm = match gap_r with
        | Some (i,r) -> m_gtr_coupled pi_ co_ a_size i r
        | None       -> m_gtr_independent pi_ co_ a_size
    in
    srm

(* val m_file :: any alphabet size -- recomputes diagonal and divides by meanrate *)
let m_file pi_ f_rr a_size =
    assert(a_size = Array.length f_rr);
    let srm = create_ba2 a_size a_size in
    for r = 0 to (a_size-1) do
        assert(a_size = Array.length f_rr.(r));
        let diag = ref 0.0 in
        for c = 0 to (a_size-1) do
            if (c <> r) then begin
                diag := !diag +. f_rr.(r).(c);
                srm.{r,c} <- f_rr.(r).(c);
            end
        done;
        srm.{r,r} <- ~-. !diag;
    done;
    m_meanrate srm pi_;
    srm

let m_custom pi_ idxs ray a_size =
    let idx = ref 0 in
    let srm = create_ba2 a_size a_size in
    for i = 0 to (a_size-1) do
        for j=i+1 to (a_size-1) do
            let value = ray.( All_sets.IntegerMap.find !idx idxs ) in
            srm.{i,j} <- value *. pi_.{j};
            srm.{j,i} <- value *. pi_.{i};
            incr idx;
        done;
    done;
    for i = 0 to (a_size-1) do begin
        let diag = ref 0.0 in
        for j = 0 to (a_size-1) do
            if (i <> j) then diag := !diag +. srm.{i,j};
        done;
        srm.{i,i} <- -. !diag;
    end; done;
    m_meanrate srm pi_;
    srm


(* create a custom model by two Maps. One that maps an index to it's linked index
   and another that maps that to a value. The default value if the element does
   not exist is set to 1.0. The mean rate is normalized and the priors are
   multiplied through. The indexes (see below) ensure that we are dealing with a
   symmetric matrix for the parameters,

   ex,  - a a b  The index association would be, 1->1, 2->1, 3->2, 4->2,
        a - b a                                  5->1, 6->1 (upper triangular).
        a b - a  Diagonal elements are always ignored, a dash is recommended.
        b a a -  We ensure that the parameters are coupled symmetrically.  *)
let process_custom_model alph_size (f_aa: char array array) =
    let found = ref All_sets.Integers.empty in
    let assoc = ref All_sets.IntegerMap.empty in
    let idx = ref 0 in
    let a_size = Array.length f_aa in
    if not (a_size = alph_size) then begin
        raise (LikelihoodModelError "Alphabet size does not match model")
    end;
    for i = 0 to a_size -1 do
        assert( Array.length f_aa.(i) = a_size );
        for j = i+1 to a_size-1 do
            let letter = Char.code f_aa.(i).(j) in
            if not ( letter = Char.code f_aa.(j).(i) ) then begin
                raise (LikelihoodModelError "Custom model should be symmetric")
            end;
            assoc := All_sets.IntegerMap.add !idx letter !assoc;
            found := All_sets.Integers.add letter !found;
            incr idx;
        done;
    done;
    let length,map =
        All_sets.Integers.fold
            (fun v1 (i,map) ->
                let map =
                    All_sets.IntegerMap.map
                        (fun v2 -> if v1 = v2 then i else v2) map
                in
                (i+1, map))
            (!found)
            (0,!assoc)
    in
    (map, Array.create length 1.0)



(* ------------------------------------------------ *)
(* MATRIX CALCULATION FUNCTIONS                     *)

(* A function to check for nan values in matrix; catch before a diagonalization,
   since the lapack routines reeturn DGEBAL parameter illegal value instead of a
   backtrace *)
let check_for_nan mat =
    try for i = 0 to (Bigarray.Array2.dim1 mat)-1 do
            for j = 0 to (Bigarray.Array2.dim1 mat)-1 do
                assert( not (Numerical.is_nan mat.{i,j}));
            done;
        done;
        true
    with _ ->
        false

(* functions to diagonalize the two types of substitution matrices *)
let diagonalize (sym : bool) mat =
    assert( check_for_nan mat );
    let alph = Bigarray.Array2.dim1 mat in
    let n_u  = create_ba2 alph alph in 
    let () = Bigarray.Array2.blit mat n_u in
    let n_d  = create_ba2 alph alph in 
    let () = Bigarray.Array2.fill n_d 0.0 in
    assert( alph = Bigarray.Array2.dim2 mat);
    let apply_sym sub_mat =
        let () = diagonalize_sym FMatrix.scratch_space n_u n_d in
        n_u, n_d, None
    and apply_gtr sub_mat = 
        let n_ui = create_ba2 alph alph in
        let () = Bigarray.Array2.fill n_ui 0.0 in
        let () = diagonalize_gtr FMatrix.scratch_space n_u n_d n_ui in
        n_u, n_d, Some n_ui
    in 
    match sym with
    | true  -> apply_sym mat
    | false -> apply_gtr mat

let compose model t = match model.ui with
    | Some ui -> compose_gtr FMatrix.scratch_space model.u model.d ui t
    | None    -> compose_sym FMatrix.scratch_space model.u model.d t

let subst_matrix model topt =
    let _gapr = match model.spec.use_gap with
        | `Coupled x   -> Some (Alphabet.get_gap (fst model.spec.alphabet), x)
        | `Independent -> None
        | `Missing     -> None
    and a_size = snd model.spec.alphabet
    and priors = model.pi_0 in
    let m = match model.spec.substitution with
        | JC69          -> m_jc69 priors 1.0 a_size _gapr
        | F81           -> m_f81 priors 1.0 a_size _gapr
        | K2P t         -> m_k2p priors t 1.0 a_size _gapr
        | F84 t         -> m_f84 priors t 1.0 a_size _gapr
        | HKY85 t       -> m_hky85 priors t a_size _gapr
        | TN93 (ts,tv)  -> m_tn93 priors ts tv 1.0 a_size _gapr
        | GTR c         -> m_gtr priors c a_size _gapr
        | File (m,s)    -> m_file priors m a_size
        | Custom (assoc,ray,_) -> m_custom priors assoc ray a_size
    in
    match topt with
    | Some t ->
        for i = 0 to (Bigarray.Array2.dim1 m)-1 do
            for j = 0 to (Bigarray.Array2.dim2 m)-1 do
                m.{i,j} <- m.{i,j} *. t;
            done;
        done;
        m
    | None ->
        m

let get_optimization_method m =
  let meth = match m.spec.cost_fn with
  | `MAL -> [Numerical.BFGS None]
  | `MPL -> [Numerical.Brent_Multi None; Numerical.BFGS None]
  in
  Numerical.default_numerical_optimization_strategy
      meth !Methods.opt_mode (count_parameters m)


(* print output in our nexus format or Phyml output *)
let output_model output output_table nexus model set = 
    let printf format = Printf.ksprintf output format in
    let gtr_mod = ref false in
    if nexus = `Nexus then begin
        printf "@[Likelihood@.";
        let () = match model.spec.substitution with
            | JC69   -> printf "@[Model = JC69;@]@\n";
            | F81    -> printf "@[Model = F81;@]";
            | K2P x  -> printf "@[Model = K2P;@]";
                        printf "@[Parameters = %f;@]" x
            | F84 x  -> printf "@[Model = F84;@]";
                        printf "@[Parameters = %f;@]" x
            | HKY85 x-> printf "@[Model = HKY85;@]";
                        printf "@[Parameters = %f;@]" x
            | TN93 (a,b) -> printf "@[Model = TN93;@]";
                        printf "@[Parameters = %f %f;@]" a b
            | GTR xs -> printf "@[Model = GTR;@]";
                        printf "@[Parameters = ";
                        Array.iter (printf "%f ") xs;
                        printf ";@]";
                        gtr_mod := true
            | File (_,str) -> 
                        printf "@[Model = File:%s;@]" str
            | Custom (_,xs,str) ->
                        printf "@[Model = Custom:%s;@]" str;
                        printf "@[Parameters = ";
                        Array.iter (printf "%f ") xs;
                        printf ";@]";
        in
        let () = match model.spec.base_priors with
            | Equal ->
                printf "@[Priors = Equal;@]"
            | Estimated x | Given x ->
                let first = ref true in
                printf "@[Priors =";
                List.iter
                    (fun (s,i) -> 
                        try if !first then begin
                            printf "@[%s %.5f" s x.(i); first := false
                        end else
                            printf ",@]@[%s %.5f" s x.(i) 
                        with _ -> ())
                    (Alphabet.to_list (fst model.spec.alphabet));
                printf ";@]@]"
        in
        let () = match model.spec.cost_fn with
            | `MPL -> printf "@[Cost = mpl;@]";
            | `MAL -> printf "@[Cost = mal;@]"; 
        in
        let () = match model.spec.site_variation with
            | Constant        -> ()
            | Gamma (c, p)    ->
                printf "@[Variation = gamma;@]@[alpha = %.5f@]@[sites = %d;@]" p c
            | Theta (c, p, i) ->
                printf ("@[Variation = theta;@]@[alpha = %.5f@]@[sites = %d;@]"
                        ^^"@[percent = %.5f;@]") p c i
        in
        let () = match model.spec.use_gap with
            (* the meaning of independent/coupled switch under gtr *)
            | `Independent -> printf "@[gap = independent;@]"
            | `Coupled x   -> printf "@[gap = coupled:%f;@]" x
            | `Missing     -> printf "@[gap = missing;@]"
        in
        let () = match set with
            | Some namelist ->
                let first = ref true in
                printf "@[CharSet = ";
                List.iter
                    (fun s -> 
                        if !first then begin printf "@[%s" s; first := false
                        end else printf ",@]@[%s" s)
                    namelist;
                printf ";@]@]";
            | None -> ()
        in
        printf ";@]@."
    (* ---------------------------- *)
    end else (* phylip *) begin
        printf "@[<hov 0>Discrete gamma model: ";
        let () = match model.spec.site_variation with
            | Constant -> printf "No@]\n";
            | Gamma (cats,param) ->
                printf ("Yes@]@\n@[<hov 1>- Number of categories: %d@]\n"^^
                        "@[<hov 1>- Gamma Shape Parameter: %.4f@]\n") cats param
            | Theta (cats,param,inv) ->
                printf ("Yes@]@\n@[<hov 1>- Number of categories: %d@]\n"^^
                        "@[<hov 1>- Gamma Shape Parameter: %.4f@]\n") cats param;
                printf ("@[<hov 1>- Proportion of invariant: %.4f@]\n") inv
        in
        let () = match model.spec.cost_fn with
            | `MPL -> printf "@[<hov 0>Cost mode: mpl;@]@\n";
            | `MAL -> printf "@[<hov 0>Cost mode: mal;@]@\n";
        in
        printf "@[<hov 0>Priors / Base frequencies:@\n";
        let () = match model.spec.base_priors with
            | Equal -> printf "@[Equal@]@]@\n"
            | Estimated x | Given x ->
                List.iter
                    (fun (s,i) ->
                        (* this expection handling avoids gaps when they are not
                         * enabled in the alphabet as an additional character *)
                        try printf "@[<hov 1>- f(%s)= %.5f@]@\n" s x.(i) with _ -> ())
                    (Alphabet.to_list (fst model.spec.alphabet));
        in
        printf "@[<hov 0>Model Parameters: ";
        let a = (snd model.spec.alphabet) in
        let () = match model.spec.substitution with
            | JC69  -> printf "JC69@]@\n"
            | F81   -> printf "F81@]@\n"
            | K2P x ->
                printf "K2P@]@\n@[<hov 1>- Transition/transversion ratio: %.5f@]@\n" x
            | F84 x ->
                printf "F84@]@\n@[<hov 1>- Transition/transversion ratio: %.5f@]@\n" x
            | HKY85 x    ->
                printf "HKY85@]@\n@[<hov 1>- Transition/transversion ratio:%.5f@]@\n" x
            | TN93 (a,b) -> 
                printf "tn93@]@\n@[<hov 1>- transition/transversion ratio:%.5f/%.5f@]@\n" a b
            | GTR ray -> gtr_mod := true;
                printf "GTR@]@\n@[<hov 1>- Rate Parameters: @]@\n";
                let get_str i = Alphabet.match_code i (fst model.spec.alphabet)
                and convert s r c = (c + (r*(s-1)) - ((r*(r+1))/2)) - 1 in
                begin match model.spec.use_gap with
                    | `Coupled x ->
                        let ray =
                            let size = (((a-3)*a)/2) in
                            Array.init (size+1)
                                       (fun i -> if i = (size) then 1.0 else ray.(i))
                        in
                        for i = 0 to a - 2 do
                            for j = i+1 to a - 2 do
                                printf "@[<hov 1>%s <-> %s - %.5f@]@\n"
                                    (get_str i) (get_str j) ray.(convert (a-1) i j)
                            done;
                        done;
                        printf "@[<hov 1>%s <-> N - %.5f@]@\n"
                            (get_str (Alphabet.get_gap (fst model.spec.alphabet))) x
                    | `Missing | `Independent ->
                        let ray =
                            let size = (((a-1)*a)/2) in
                            Array.init (size)
                                       (fun i -> if i = (size-1) then 1.0 else ray.(i))
                        in
                        for i = 0 to a - 1 do
                            for j = i+1 to a - 1 do
                                printf "@[<hov 1>%s <-> %s - %.5f@]@\n" (get_str i)
                                        (get_str j) ray.(convert a i j)
                            done;
                        done
                end
            | File (ray,name) ->
                printf "File:%s@]@\n" name;
                let mat = compose model 0.0 in
                printf "@[<hov 1>[";
                for i = 0 to a - 1 do
                    printf "%s ------- " (Alphabet.match_code i (fst model.spec.alphabet))
                done;
                printf "]@]@\n";
                for i = 0 to a - 1 do
                    for j = 0 to a - 1 do
                        printf "%.5f\t" mat.{i,j}
                    done;
                    printf "@\n";
                done;
            | Custom (_,xs,str) ->
                printf "@[Custom:%s;@]@\n" str;
                printf "@[Parameters : ";
                Array.iter (printf "%f ") xs;
                printf ";@]";
        in
        let () = match model.spec.use_gap with
            | `Independent -> printf "@[<hov 0>Gap property: independent;@]@\n"
            | `Coupled x   -> printf "@[<hov 0>Gap property: coupled, Ratio: %f;@]@\n" x
            | `Missing     -> printf "@[<hov 0>Gap property: missing;@]@\n"
        in
        printf "@]";
        printf "@[@[<hov 0>Instantaneous rate matrix:@]@\n";
        match output_table with
        | Some output_table ->
            let table = Array.make_matrix (a+1) a "" in
            let () =
                let mat = compose model ~-.1.0 in
                for i = 0 to a - 1 do
                    table.(0).(i) <- Alphabet.match_code i (fst model.spec.alphabet);
                done;
                for i = 0 to a - 1 do 
                    for j = 0 to a - 1 do
                        table.(i+1).(j) <- string_of_float mat.{i,j}
                    done;
                done;
            in
            let () = output_table table in
            printf "@\n@]"
        | None ->
            let () =
                let mat = compose model ~-.1.0 in
                printf "@[<hov 1>[";
                for i = 0 to a - 1 do
                    printf "%s ------- " (Alphabet.match_code i (fst model.spec.alphabet))
                done;
                printf "]";
                for i = 0 to a - 1 do
                    printf "@]@\n@[<hov 1>";
                    for j = 0 to a - 1 do
                        printf "%8.5f  " mat.{i,j}
                    done;
                done;
            in
            printf "@]@\n@]"
    end

(** [compose_model] compose a substitution probability matrix from branch length
    and substitution rate matrix. *)
let compose_model sub_mat t = 
    let a_size = Bigarray.Array2.dim1 sub_mat in
    let (u_,d_,ui_) = 
        let n_d = Bigarray.Array2.create Bigarray.float64 Bigarray.c_layout a_size a_size
        and n_ui = Bigarray.Array2.create Bigarray.float64 Bigarray.c_layout a_size a_size in
        Bigarray.Array2.fill n_d 0.0;
        Bigarray.Array2.fill n_ui 0.0;
        diagonalize_gtr FMatrix.scratch_space sub_mat n_d n_ui;
        (sub_mat, n_d, n_ui)
    in
    compose_gtr FMatrix.scratch_space u_ d_ ui_ t

(** [integerized_model] convert a model and branch length to a probability rate
    matrix and then scale by a factor of accuracy and convert to integers. *)
let integerized_model ?(sigma=4) model t =
    let sigma = 10.0 ** (float_of_int sigma)
    and create = match model.spec.substitution with
        | JC69 | K2P _ ->
            (fun s m i j -> ~-(int_of_float (s *.(log m.(i).(j)))))
        | _ ->
            (fun s m i j -> ~-(int_of_float (s *.(log (model.pi_0.{i} *. m.(i).(j))))))
    and matrix =
        barray_to_array2 (match model.ui with
        | Some ui -> compose_gtr FMatrix.scratch_space model.u model.d ui t
        | None    -> compose_sym FMatrix.scratch_space model.u model.d t)
    and imatrix =
        Array.make_matrix (snd model.spec.alphabet) (snd model.spec.alphabet) 0
    in
    assert( (snd model.spec.alphabet) = Array.length matrix);
    assert( (snd model.spec.alphabet) = Array.length matrix.(0));
    for i = 0 to (Array.length matrix) - 1 do
        for j = 0 to (Array.length matrix.(0)) - 1 do
            imatrix.(i).(j) <- create sigma matrix i j
        done;
    done;
    imatrix

(** create a cost matrix from a model *)
let model_to_cm model t =
    let input = let t = max Numerical.minimum t in integerized_model model t in
    let llst = Array.to_list (Array.map Array.to_list input) in
    let res = Cost_matrix.Two_D.of_list ~suppress:true llst (snd model.spec.alphabet) in
    res

ELSE

    let output_model _ _ _ _ = lfailwith ()
    let compose _ _  = lfailwith ()
    let spec_from_classification _ _ _ _ _ _ = lfailwith ()
    let compare _ _  = lfailwith ()
    let classify_seq_pairs _ _ _ _ _ = lfailwith ()
    let subst_matrix _ _  = lfailwith ()
    let process_custom_model _ = lfailwith ()
    let model_to_cm _ _   = lfailwith ()
    let compose_model _ _ = lfailwith ()
    let m_gtr  _ _ _ _ = lfailwith ()
    let m_file _ _ _   = lfailwith ()
    let m_jc69 _ _ _   = lfailwith ()

END


(* ------------------------------------------------ *)
(* CONVERSION/MODEL CREATION FUNCTIONS             *)

(* convert a string spec (from nexus, for example) to a specification *)
let convert_string_spec alph ((name,(var,site,alpha,invar),param,priors,gap,cost,file):string_spec) =
  IFDEF USE_LIKELIHOOD THEN
    let gap_info = match String.uppercase (fst gap),snd gap with
        | "COUPLED", None   -> `Coupled default_gap_r
        | "COUPLED", Some x -> `Coupled x
        | "INDEPENDENT",_   -> `Independent
        | "MISSING",_       -> `Missing
        | "",_              -> `Missing
        | x,_        ->
            raise (LikelihoodModelError ("Invalid gap property, "^x))
    in
    let submatrix = match String.uppercase name with
        | "JC69" -> begin match param with
            | [] -> JC69
            | _ -> failwith "Parameters don't match model" end
        | "F81" -> begin match param with
            | [] -> F81
            | _ -> failwith "Parameters don't match model" end
        | "K80" | "K2P" -> begin match param with
            | ratio::[] -> K2P ratio
            | []        -> K2P default_tstv
            | _ -> failwith "Parameters don't match model" end
        | "F84" -> begin match param with
            | ratio::[] -> F84 ratio
            | []        -> F84 default_tstv
            | _ -> failwith "Parameters don't match model" end
        | "HKY" | "HKY85" -> begin match param with
            | ratio::[] -> HKY85 ratio
            | []        -> HKY85 default_tstv
            | _ -> failwith "Parameters don't match model" end
        | "TN93" -> begin match param with
            | ts::tv::[] -> TN93 (ts,tv)
            | []         -> TN93 (default_tstv,default_tstv)
            | _ -> failwith "Parameters don't match model" end
        | "GTR" -> begin match param with
            | [] -> GTR [||]
            | ls -> GTR (Array.of_list ls) end
        | "GIVEN"->
            begin match file with
                | Some name ->
                    name --> FileStream.read_floatmatrix
                         --> List.map (Array.of_list)
                         --> Array.of_list
                         --> (fun x -> File (x,name))
                | None ->
                    raise (LikelihoodModelError "No File specified for model")
            end
        | "CUSTOM" ->
            let alph_size = snd alph in
            begin match file with
                | Some name ->
                    let convert str = assert( String.length str = 1 ); String.get str 0 in
                    let matrix = Cost_matrix.Two_D.matrix_of_file convert (`Local name) in
                    let matrix = Array.of_list (List.map (Array.of_list) matrix) in
                    let assoc,ray = process_custom_model alph_size matrix in
                    Custom (assoc,ray,name)
                | None ->
                    raise (LikelihoodModelError "No File specified for model")
            end
        (* ERROR *)
        | "" -> raise (LikelihoodModelError "No model specified")
        | xx -> raise (LikelihoodModelError ("Unknown likelihood model "^xx))
    in
    let cost_fn : Methods.ml_costfn = match String.uppercase cost with
        | "MPL" -> `MPL
        | "MAL" -> `MAL
        | ""    -> `MAL (* unmentioned default *)
        | x     -> raise (LikelihoodModelError ("Unknown cost mode "^x))
    and variation = match String.uppercase var with
        | "GAMMA" ->
            let alpha = 
                try float_of_string alpha
                with _ -> default_alpha false
            in
            Gamma (int_of_string site, alpha)
        | "THETA" -> 
            let alpha =
                try float_of_string alpha
                with _ -> default_alpha true
            in
            let invar =
                try float_of_string invar
                with _ -> default_invar
            in
            Theta (int_of_string site, alpha, invar)
        | "NONE" | "CONSTANT" | "" ->
            Constant
        | x -> raise (LikelihoodModelError ("Unknown rate variation mode "^x))
    and priors = match priors with
        | `Given priors        -> Given (Array.of_list (List.map snd priors))
        | `Equal               -> Equal
        | `Estimate (Some x)   -> Estimated x
        | `Consistent pre_calc ->
            begin match submatrix, pre_calc with
                | JC69, _
                | K2P _, _   -> Equal
                | _, Some pi -> Estimated pi
                | _ , None   -> assert false
            end
        | `Estimate None -> assert false
    in
    {   substitution = submatrix;
        site_variation = variation;
        base_priors = priors;
        use_gap = gap_info;
        alphabet = alph;
        cost_fn = cost_fn; }
  ELSE
    lfailwith ()
  END


(** Convert Methods.ml specification to that of an MlModel spec *)
let convert_methods_spec (alph,alph_size) (compute_priors) 
        ((_,talph,cst,subst,site_variation,base_priors,use_gap):Methods.ml_spec) =
    let u_gap = match use_gap with 
        | `Independent | `Coupled _ -> true | `Missing -> false in
    let alph_size = 
        let w_gap = if u_gap then alph_size else alph_size - 1 in
        match talph with | `Min | `Max  -> w_gap | `Int x -> x
    in
    let base_priors = match base_priors with
        | `Estimate  -> Estimated (compute_priors ())
        | `Equal     -> Equal
        | `Given arr -> Given (Array.of_list arr)
        | `Consistent ->
            begin match subst with
                | `JC69 | `K2P _ -> Equal
                | _ -> Estimated (compute_priors ())
            end
    and site_variation = match site_variation with
        | None -> Constant 
        | Some (`Gamma (w,y)) -> 
            let y = match y with
                | Some x -> x 
                | None   -> default_alpha false
            in
            Gamma (w,y)
        | Some (`Theta (w,y)) ->
            let y,p = match y with 
                | Some x -> x 
                | None   -> default_alpha true, default_invar
            in
            Theta (w,y,p)
    and substitution = match subst with
        | `AIC _ | `BIC _ | `AICC _ | `NCM -> assert false
        | `JC69    -> JC69
        | `F81     -> F81
        | `K2P [x] -> K2P x
        | `K2P []  -> K2P default_tstv
        | `K2P _   ->
                raise (LikelihoodModelError
                    "Likelihood model K2P requires 1 or 0 parameters")
        | `HKY85 [x] -> HKY85 x
        | `HKY85 []  -> HKY85 default_tstv
        | `HKY85 _   ->
                raise (LikelihoodModelError
                    "Likelihood model HKY85 requires 1 or 0 parameters")
        | `F84 [x]   -> F84 x
        | `F84 []    -> F84 default_tstv
        | `F84 _     ->
                raise (LikelihoodModelError
                    "Likelihood model F84 requires 1 or 0 parameters")
        | `TN93 [x;y] -> TN93 (x,y)
        | `TN93 []    -> TN93 (default_tstv,default_tstv)
        | `TN93 _     ->
                raise (LikelihoodModelError
                    "Likelihood model TN93 requires 2 or 0 parameters")
        | `GTR xs -> GTR (Array.of_list xs)
        | `File str ->
            (* this needs to be changed to allow remote files as well *)
            let matrix = Cost_matrix.Two_D.matrix_of_file float_of_string (`Local str) in
            let matrix = Array.of_list (List.map (Array.of_list) matrix) in
            Array.iter
                (fun x ->
                    if Array.length x = alph_size then ()
                    else begin
                        raise (LikelihoodModelError "Likelihood model TN93 requires 2 or 0 parameters")
                    end)
                (matrix);
            if Array.length matrix = alph_size then
                File (matrix,str)
            else begin
                raise (LikelihoodModelError "Likelihood model TN93 requires 2 or 0 parameters")
            end;
        | `Custom str ->
            let convert str = assert( String.length str = 1 ); String.get str 0 in
            let matrix = Cost_matrix.Two_D.matrix_of_file convert (`Local str) in
            let matrix = Array.of_list (List.map (Array.of_list) matrix) in
            let assoc,ray = process_custom_model alph_size matrix in
            Custom (assoc,ray,str)
    in
    { substitution = substitution;
    site_variation = site_variation;
       base_priors = base_priors;
          alphabet = (alph,alph_size);
           cost_fn = cst;
           use_gap = use_gap; }


(** check the rates so SUM(r_k*p_k) == 1 and SUM(p_k) == 1 |p| == |r| *)
let verify_rates probs rates =
    let p1 = (Bigarray.Array1.dim probs) = (Bigarray.Array1.dim rates)
    and p2 =
        let rsps = ref 0.0 and ps = ref 0.0 in
        for i = 0 to (Bigarray.Array1.dim probs) - 1 do
            rsps := !rsps +. (probs.{i} *.  rates.{i});
            ps := !ps +. probs.{i};
        done;
        !rsps =. 1.0 && !ps =. 1.0
    in
    p1 && p2

(** create a model based on a specification and an alphabet *)
let create ?(min_prior=Numerical.minimum) lk_spec = 
  IFDEF USE_LIKELIHOOD THEN
    let (alph,a_size) = lk_spec.alphabet in
    assert( a_size > 1 );
    (* set up all the probability and rates *)
    let variation,probabilities,invar =
        match lk_spec.site_variation with
            | Constant ->
                ba_of_array1 [| 1.0 |], ba_of_array1 [| 1.0 |], None
            | Gamma (x,y) -> (* SITES,ALPHA *)
                let p = create_ba1 x in
                Bigarray.Array1.fill p (1.0 /. (float_of_int x));
                Numerical.gamma_rates y y x,p,None
            | Theta (w,x,z) -> (* GAMMA_SITES,ALPHA,PERCENT_INVAR *)
                if w < 1 then
                    failwith "Number of rate categories must be >= 1 (if invar w/out gamma)"
                else if w = 1 then begin
                    ba_of_array1 [| 1.0 |], ba_of_array1 [| 1.0 |], Some z
                end else begin
                    let p = create_ba1 w in
                    Bigarray.Array1.fill p (1.0 /. (float_of_int w));
                    let r = Numerical.gamma_rates x x w in
                    r,p,Some z
                end
    in
    assert( verify_rates probabilities variation );
    (* extract the prior probability *)
    let priors =
        let p = match lk_spec.base_priors with 
            | Equal ->
                Array.make a_size (1.0 /. (float_of_int a_size))
            | Estimated p | Given p -> 
                let p = Array.map (fun i -> max min_prior i) p in
                let sum = Array.fold_left (fun a b -> a +. b) 0.0 p in 
                if sum =. 1.0 
                    then p
                    else Array.map (fun i -> i /. sum) p
        in
        if not (a_size = Array.length p) then begin
            debug_printf "Alphabet = %d; Priors = %d\n%!" a_size (Array.length p);
            raise (LikelihoodModelError "Priors don't match alphabet");
        end;
        ba_of_array1 p
    in
    let _gapr = match lk_spec.use_gap with
        | `Coupled x   ->
            let gap_i =
              if Alphabet.zero_indexed alph
                then Alphabet.get_gap alph
                else (Alphabet.get_gap alph)-1
            in
            Some (gap_i, x)
        | `Independent -> None
        | `Missing     -> None
    in
    (*  get the substitution rate matrix and set sym variable and to_formatter vars *)
    let sym, sub_mat, subst_model = match lk_spec.substitution with
        | JC69  -> true,  m_jc69 priors 1.0 a_size _gapr, JC69
        | F81   -> false, m_f81 priors 1.0 a_size _gapr, F81
        | K2P t -> true,  m_k2p priors t 1.0 a_size _gapr, lk_spec.substitution
        | F84 t -> false, m_f84 priors t 1.0 a_size _gapr, lk_spec.substitution
        | HKY85 t -> false, m_hky85 priors t a_size _gapr, lk_spec.substitution
        | TN93 (ts,tv) ->
            false, m_tn93 priors ts tv 1.0 a_size _gapr, lk_spec.substitution
        | GTR c ->
            let c = match (Array.length c), lk_spec.use_gap with
                | 0, `Coupled _ -> default_gtr a_size true
                | 0, _          -> default_gtr a_size false
                | _, _          -> c
            in
            false, m_gtr priors c a_size _gapr, (GTR c)
        | File (m,s)->
            false, m_file priors m a_size, lk_spec.substitution
        | Custom (assoc,xs,_) ->
            false, m_custom priors assoc xs a_size, lk_spec.substitution
    in
    (* ensure that when priors are not =, we use a model that asserts that *)
    let lk_spec = match lk_spec.base_priors with
        | Estimated _  when sym ->
            Status.user_message Status.Warning 
                "I@ am@ using@ equal@ priors@ as@ required@ in@ the@ selected@ likelihood@ model.";
            { lk_spec with base_priors = Equal }
        | Given _ when sym ->
            Status.user_message Status.Warning
                "I@ am@ using@ equal@ priors@ as@ required@ in@ the@ selected@ likelihood@ model.";
            { lk_spec with base_priors = Equal }
        | (Estimated _ | Equal | Given _) -> lk_spec
    in
    let (u_,d_,ui_) = diagonalize sym sub_mat in
    {
      rate = variation;
      prob = probabilities;
      spec = {lk_spec with substitution = subst_model; };
     invar = invar;
      pi_0 = priors;
         s = sub_mat;
         u = u_;
         d = d_;
        ui = ui_;
    }
  ELSE
    lfailwith ()
  END


(** Print a substitution probability matrix *)
let debug_model model t =
    let subst = subst_matrix model t in
    print_barray2 subst; print_newline ();
    match t with | Some t -> print_barray2 (compose model t) | None -> ()


(** Replace the priors in a model with that of an array *)
let replace_priors model array = 
    if debug then begin
        Printf.printf "Replacing Priors\n\t%!";
        Array.iter (fun x -> Printf.printf "%f, " x) array;
        print_newline ();
    end;
    create {model.spec with base_priors = Estimated array}


(** Compute the priors of a dataset by frequency and gap-counts *)
let compute_priors (alph,u_gap) freq_ (count,gcount) lengths : float array = 
    let size = if u_gap then (Alphabet.size alph) else (Alphabet.size alph)-1 in
    let gap_contribution = (float_of_int gcount) /. (float_of_int size) in
    let gap_char = Alphabet.get_gap alph in
    if debug then begin
        Printf.printf "Computed Priors of %d char + %d gaps: " count gcount;
        Array.iter (Printf.printf "|%f") freq_;
        Printf.printf "|]\n%!";
    end;
    let final_priors = 
        if u_gap then begin
            let total_added_gaps =
                let longest : int = List.fold_left (fun a x-> max a x) 0 lengths in
                let add_gap : int = List.fold_left (fun acc x -> (longest - x) + acc) 0 lengths in
                float_of_int add_gap
            in
            Printf.printf "Total added gaps = %f\n%!" total_added_gaps;
            freq_.(gap_char) <- freq_.(gap_char) +. total_added_gaps;
            let count = (float_of_int (count - gcount)) +. total_added_gaps;
            and weight  = (float_of_int gcount) /. (float_of_int size) in
            Array.map (fun x -> (x -. weight) /. count) freq_
        end else begin
            Array.map (fun x -> (x -. gap_contribution) /. (float_of_int count)) freq_
        end
    in
    let sum = Array.fold_left (fun a x -> a +. x) 0.0 final_priors in
    if debug then begin 
        Printf.printf "Final Priors (%f): [" sum;
        Array.iter (Printf.printf "|%f") final_priors;
        Printf.printf "|]\n%!";
    end;
    final_priors

(** Add Independent Gap to a model *)
let add_gap_to_model compute_priors model = 
    raise (LikelihoodModelError dyno_likelihood_warning)

let add_gap_to_model compute_priors model = 
    match model.spec.use_gap with
    | `Missing -> add_gap_to_model compute_priors model
    | `Independent | `Coupled _ -> model

let remove_gamma_from_spec spec =
    match spec.site_variation with
    | Constant          -> spec
    | Gamma _ | Theta _ ->
        Status.user_message Status.Warning
            (Str.global_replace (Str.regexp " ") "@ " dyno_gamma_warning);
        { spec with site_variation = Constant; }

IFDEF USE_LIKELIHOOD THEN

(* ------------------------------------------------ *)
(* MODEL ESTIMATION FUNCTIONS                       *)
let chars2str = function
    | `Packed s -> string_of_int s
    | `List s   -> "[ "^(String.concat " | " (List.map string_of_int s))^" ]"

(** estimate the model based on two sequences with attached weights *)
let classify_seq_pairs leaf1 leaf2 seq1 seq2 acc =
    let chars_to_list = function
        | `Packed s -> BitSet.Int.list_of_packed s
        | `List   s -> s
    and incr_map k v map =
        let v = match All_sets.IntegerMap.mem k map with
            | true  -> v +. (All_sets.IntegerMap.find k map)
            | false -> v 
        in
        All_sets.IntegerMap.add k v map 
    and incr_tuple ((a,b) as k) v map =
        let v = match All_sets.FullTupleMap.mem k map with
            | true -> v +. (All_sets.FullTupleMap.find k map)
            | false -> v 
        in
        All_sets.FullTupleMap.add k v map
    in
    (* classify the mutation from a1 to a2 in a map *)
    let mk_transitions (tmap,fmap) (w1,c1,s1) (w2,c2,s2) =
        assert( w1 = w2 );
        let s1 = chars_to_list s1 and s2 = chars_to_list s2 in 
        let n1 = List.length s1   and n2 = List.length s2   in
        (* Count all the transitions A->B. Polymorphisms are counted as
         * 1/(na*nb) each, where na and nb are the number of polymorphisms
         * in and b respectively. *)
        let v = w1 /. (float_of_int (n1 * n2) ) in
        let tmap = (* cross product of states a and b *)
            List.fold_left 
                (fun map1 a ->
                    List.fold_left
                        (fun map2 b -> incr_tuple (a,b) v map2)
                        map1 s2)
                tmap s1
        in
        (* Counts the base frequencies. Polymorphisms are counted as 1/n in
         * each base, where, as above, n is the number of polymorphisms. *)
        let n1 = w1 /. (float_of_int n1) and n2 = w2 /. (float_of_int n2) in
        let fmap = 
            if leaf1 then 
                List.fold_left (fun map a -> incr_map a n1 map) fmap s1 
            else fmap in
        let fmap = 
            if leaf2 then
                List.fold_left (fun map a -> incr_map a n2 map) fmap s2
            else fmap in
        tmap,fmap
    in
    List.fold_left2 mk_transitions acc seq1 seq2

let get_priors f_prior alph code = match f_prior with
    | Estimated x -> x.(Alphabet.match_base code alph)
    | Given     x -> x.(Alphabet.match_base code alph)
    | Equal       -> 1.0 /. (float_of_int (Alphabet.size alph))

(* Develop a model from a classification --created above *)
let spec_from_classification alph gap kind rates (priors:Methods.ml_priors) costfn (comp_map,pis) =
    let tuple_sum = 
        All_sets.FullTupleMap.fold (fun k v a -> a +. v) comp_map 0.0
    and ugap = match gap with
        | `Missing      -> false
        | `Independent  -> true
        | `Coupled _    -> true
    in
    let f_priors,a_size = match priors,kind with
        | `Equal,_ | (`Consistent | `Estimate), (`JC69 | `K2P _) ->
            let size =
                if ugap then (Alphabet.size alph) else (Alphabet.size alph)-1
            in
            Equal,size
        | `Consistent,_ | `Estimate,_ ->
            let sum = All_sets.IntegerMap.fold (fun k v x -> v +. x) pis 0.0
            and gap_size =
                try All_sets.IntegerMap.find (Alphabet.get_gap alph) pis
                with | Not_found -> 0.0
            in
            let l =
                let sum = if ugap then sum else sum -. gap_size in
                List.fold_left
                    (fun acc (r,b) -> match b with
                        | b when (not ugap) && (b = Alphabet.get_gap alph) -> acc
                        | b ->
                            let c =
                                try (All_sets.IntegerMap.find b pis) /. sum
                                with Not_found -> 0.0 in
                            c :: acc)
                    []
                    (Alphabet.to_list alph)
            in
            let ray = Array.of_list (List.rev l) in
            Estimated ray, Array.length ray
        | `Given xs, _ ->
            let ray = Array.of_list xs in
            Given ray, Array.length ray
    and is_comp a b =
        (* these models assume nucleotides only; T<->C=1, A<->G=2, this is
            because this is used in DNA/nucleotide models only (k2p,hky85...) *)
             if (Alphabet.match_base "T" alph) = a then
                 if (Alphabet.match_base "C" alph) = b then 1 else 0
        else if (Alphabet.match_base "C" alph) = a then
                 if (Alphabet.match_base "T" alph) = b then 1 else 0

        else if (Alphabet.match_base "A" alph) = a then
                 if (Alphabet.match_base "G" alph) = b then 2 else 0
        else if (Alphabet.match_base "G" alph) = a then
                 if (Alphabet.match_base "A" alph) = b then 2 else 0
        else
            0
    in
    let get_sva comp_map =
        let s1,s2,v,a =
            All_sets.FullTupleMap.fold
                (fun k v (sc1,sc2,vc,all) -> match k with
                    | k1,k2 when 1 == is_comp k1 k2 -> (sc1+.v,sc2,vc,all+.v)
                    | k1,k2 when 2 == is_comp k1 k2 -> (sc1,sc2+.v,vc,all+.v)
                    | k1,k2 when k1 != k2           -> (sc1,sc2,vc+.v,all+.v)
                    | _                             -> (sc1,sc2,vc,all+.v) )
                comp_map
                (0.0,0.0,0.0,0.0)
        in
        s1 /. a, s2 /. a, v /. a, a
    in
    (* build the model *)
    let m,gap = try match kind with
        | `JC69 -> JC69,gap
        | `F81  -> F81,gap
        | `K2P _ ->
            (* YANG: 1.12                                       *)
            (*  alpha*t = -0.5 * log(1-2S-V) + 0.25 * log(1-2V) *)
            (* 2*beta*t = -0.5 * log(1-2V)                      *)   
            (*        k = alpha / beta                          *)  
            (* k = [log(1-2S-V)-0.5*log(1-2V)] / log(1-2V)      *)
            let s1,s2,v,a = get_sva comp_map in
            let s = s1 +. s2 in
            let numer = 2.0 *. (log (1.0-.(2.0*.s)-.v))
            and denom = log (1.0-.(2.0*.v)) in
            K2P ((numer /. denom)-.1.0),gap
        | `GTR _ ->
            let tuple_sum a1 a2 map =
                let one = try All_sets.FullTupleMap.find (a1,a2) map 
                          with | Not_found -> 0.0
                and two = try All_sets.FullTupleMap.find (a2,a1) map
                          with | Not_found -> 0.0
                in
                one +. two
            in
            (* create list of transitions for GTR model creation *)
            (* 1 -> 2, 1 -> 3, 1 -> 4 ... 2 -> 3 ... *)
            begin match gap with
                | `Independent | `Missing ->
                    let cgap = Alphabet.get_gap alph in
                    let lst =
                        List.fold_right
                            (fun (s1,alph1) acc1 ->
                                if (not ugap) && (cgap = alph1) then
                                    acc1
                                else
                                    List.fold_right
                                        (fun (s2,alph2) acc2 ->
                                            if alph2 <= alph1 then acc2
                                            else if (not ugap) && (cgap = alph2) then
                                                acc2
                                            else begin
                                                let sum = tuple_sum alph1 alph2 comp_map in
                                                sum :: acc2 
                                            end)
                                        (Alphabet.to_list alph) acc1)
                            (Alphabet.to_list alph) []
                    in
                    let sum = List.fold_left (fun a x -> x +. a) 0.0 lst in
                    let lst = List.map (fun x -> x /. sum) lst in
                    let arr,gap = normalize alph gap (Array.of_list lst) in
                    let arr = Array.init ((List.length lst)-1) (fun i -> arr.(i)) in
                    GTR arr,gap
                | `Coupled _ ->
                    let cgap = Alphabet.get_gap alph in
                    let lst,gap_trans =
                        List.fold_right
                            (fun (s1,alph1) acc1 ->
                                if cgap = alph1 then acc1
                                else
                                    List.fold_right
                                        (fun (s2,alph2) ((chrt,gapt) as acc2) ->
                                            if alph2 <= alph1    then acc2
                                            else if cgap = alph2 then
                                                let sum = tuple_sum alph1 alph2 comp_map in
                                                (chrt, sum +. gapt)
                                            else begin
                                                let sum = tuple_sum alph1 alph2 comp_map in
                                                (sum :: chrt,gapt)
                                            end)
                                        (Alphabet.to_list alph) acc1)
                            (Alphabet.to_list alph)
                            ([],0.0)
                    in
                    let sum = List.fold_left (fun a x -> x +. a) gap_trans lst in
                    let lst = List.map (fun x -> x /. sum) lst in
                    let gap = `Coupled (gap_trans /. (sum *. (float_of_int a_size))) in
                    let arr,gap = normalize alph gap (Array.of_list lst) in
                    let arr = Array.init ((List.length lst)-1) (fun i -> arr.(i)) in
                    GTR arr,gap
            end
        | `F84 _ ->
            let pi_a = get_priors f_priors alph "A"
            and pi_c = get_priors f_priors alph "C"
            and pi_g = get_priors f_priors alph "G"
            and pi_t = get_priors f_priors alph "T" in
            let pi_y = pi_c +. pi_t and pi_r = pi_a +. pi_g in
            let s1,s2,v,a = get_sva comp_map in let s = s1 +. s2 in
            let a = ~-. (log (1.0 -. (s /. (2.0 *. (pi_t*.pi_c/.pi_y +. pi_a*.pi_g/.pi_r)))
                                 -. (v *.(pi_t*.pi_c *.pi_r/.pi_y +.
                                         (pi_a*.pi_g*.pi_y/.pi_r))) /.
                                    (2.0 *. (pi_t*.pi_c *. pi_r +. pi_a *. pi_g *. pi_y))))
            and b = ~-. (log (1.0 -. (v/.(2.0 *. pi_y*.pi_r)))) in
            F84 (a/.b -. 1.0),gap
        | `HKY85 _ ->
            let pi_a = get_priors f_priors alph "A"
            and pi_c = get_priors f_priors alph "C"
            and pi_g = get_priors f_priors alph "G"
            and pi_t = get_priors f_priors alph "T" in
            let pi_y = pi_c +. pi_t and pi_r = pi_a +. pi_g in
            let s1,s2,v,a = get_sva comp_map in let s = s1 +. s2 in
            let a = ~-. (log (1.0 -. (s /. (2.0 *. (pi_t*.pi_c/.pi_y +. pi_a*.pi_g/.pi_r)))
                                 -. (v *.(pi_t*.pi_c *.pi_r/.pi_y +.
                                         (pi_a*.pi_g*.pi_y/.pi_r))) /.
                                    (2.0 *. (pi_t*.pi_c *. pi_r +. pi_a *. pi_g *. pi_y))))
            and b = ~-. (log (1.0 -. (v/.(2.0 *. pi_y*.pi_r)))) in
            HKY85 (a/.b -. 1.0),gap
        | `TN93 _ ->
            let pi_a = get_priors f_priors alph "A"
            and pi_c = get_priors f_priors alph "C"
            and pi_g = get_priors f_priors alph "G"
            and pi_t = get_priors f_priors alph "T" in
            let pi_y = pi_c +. pi_t and pi_r = pi_a +. pi_g
            and s1,s2,v,a = get_sva comp_map in
            let a1 = ~-. (log (1.0 -. (pi_y*.s1/.(2.0*.pi_t*.pi_c)) -.  (v/.(2.0*.pi_y))))
            and a2 = ~-. (log (1.0 -. (pi_r*.s2/.(2.0*.pi_a*.pi_g)) -.  (v/.(2.0*.pi_r))))
            and b  = ~-. (log (1.0 -. (v /. (2.0 *. pi_y *. pi_r)))) in
            (* finally compute the ratios *)
            let k1 = (a1 -. (pi_r *. b)) /. (pi_y *. b)
            and k2 = (a2 -. (pi_y *. b)) /. (pi_r *. b) in 
            TN93 (k1,k2),gap
        | `Custom _
        | `File _ -> failwith "I cannot estimate this type of model"
        with | Not_found -> failwith "Cannot find something for model"
    and calc_invar all comp_map =
        let same = 
            List.fold_left
                (fun acc (_,ac) ->
                    acc +. (try (All_sets.FullTupleMap.find (ac,ac) comp_map)
                            with | Not_found -> 0.0))
                0.0
                (Alphabet.to_list alph)
        in
        same /. all
    in
    let v = match rates with
        | None                -> Constant
        | Some (`Gamma (i,_)) -> Gamma (i,default_alpha false)
        | Some (`Theta (i,_)) -> Theta (i,default_alpha true,calc_invar tuple_sum comp_map)
    in
    {
        substitution = m;
        site_variation = v;
        base_priors = f_priors;
        cost_fn = costfn;
        use_gap = gap;
        alphabet = (alph, a_size);
    }

let convert_gapr m = function
    | `Coupled x   -> Some (Alphabet.get_gap (fst m.spec.alphabet),x)
    | `Independent -> None
    | `Missing     -> None

let update_jc69 old_model gap_r = 
    let subst_spec = { old_model.spec with use_gap = gap_r; }
    and gap_r = convert_gapr old_model gap_r in
    let subst_model =
        m_jc69 old_model.pi_0 1.0 (snd old_model.spec.alphabet) gap_r in
    let u,d,ui = diagonalize true subst_model in
    { old_model with spec = subst_spec; s  = subst_model; u  = u; d  = d; ui = ui; }

and update_f81 old_model gap_r = 
    let subst_spec = { old_model.spec with use_gap = gap_r; }
    and gap_r = convert_gapr old_model gap_r in
    let subst_model =
        m_f81 old_model.pi_0 1.0 (snd old_model.spec.alphabet) gap_r in
    let u,d,ui = diagonalize false subst_model in
    { old_model with spec = subst_spec; s  = subst_model; u  = u; d  = d; ui = ui; }

and update_k2p old_model new_value gap_r =
    let subst_spec = { old_model.spec with substitution = K2P new_value;
                                           use_gap = gap_r; }
    and gap_r = convert_gapr old_model gap_r in
    let subst_model =
        m_k2p old_model.pi_0 new_value 1.0 (snd old_model.spec.alphabet) gap_r in
    let u,d,ui = diagonalize true subst_model in
    { old_model with spec = subst_spec; s  = subst_model; u  = u; d  = d; ui = ui; }

and update_hky old_model new_value gap_r =
    let subst_spec = { old_model.spec with substitution = HKY85 new_value;
                                           use_gap = gap_r; }
    and gap_r = convert_gapr old_model gap_r in
    let subst_model =
        m_hky85 old_model.pi_0 new_value (snd old_model.spec.alphabet) gap_r in
    let u,d,ui = diagonalize false subst_model in
    { old_model with spec = subst_spec; s  = subst_model; u  = u; d  = d; ui = ui; }

and update_tn93 old_model ((x,y) as new_value) gap_r =
    let subst_spec = { old_model.spec with substitution = TN93 new_value;
                                           use_gap = gap_r; }
    and gap_r = convert_gapr old_model gap_r in
    let subst_model =
        m_tn93 old_model.pi_0 x y 1.0 (snd old_model.spec.alphabet) gap_r in
    let u,d,ui = diagonalize false subst_model in
    { old_model with spec = subst_spec; s  = subst_model; u  = u; d  = d; ui = ui; }

and update_f84 old_model new_value gap_r =
    let subst_spec = { old_model.spec with substitution = F84 new_value;
                                           use_gap = gap_r; }
    and gap_r = convert_gapr old_model gap_r in
    let subst_model =
        m_f84 old_model.pi_0 new_value 1.0 (snd old_model.spec.alphabet) gap_r in
    let u,d,ui = diagonalize false subst_model in
    { old_model with spec = subst_spec; s  = subst_model; u  = u; d  = d; ui = ui; }

and update_gtr old_model new_values gap_r =
    let subst_spec = { old_model.spec with substitution = GTR new_values;
                                           use_gap = gap_r; }
    and gap_r = convert_gapr old_model gap_r in
    let subst_model =
        m_gtr old_model.pi_0 new_values (snd old_model.spec.alphabet) gap_r in
    let u,d,ui = diagonalize false subst_model in
    { old_model with spec = subst_spec; s  = subst_model; u  = u; d  = d; ui = ui; }

and update_custom old_model new_values =
    let subst_spec,assoc = match old_model.spec.substitution with
        | Custom (assoc,_,s) -> 
            { old_model.spec with substitution = Custom (assoc,new_values,s); },assoc
        | _ -> assert false
    in
    let subst_model =
        m_custom old_model.pi_0 assoc new_values (snd old_model.spec.alphabet) in
    let u,d,ui = diagonalize false subst_model in
    { old_model with spec = subst_spec; s = subst_model; u = u; d = d; ui = ui; }

and update_alpha old_model new_value =
    let new_spec_var,new_array = match old_model.spec.site_variation with
        | Gamma (i,_)   -> Gamma (i,new_value),
                           Numerical.gamma_rates new_value new_value i
        | Theta (i,_,p) -> Theta (i,new_value,p),
                           Numerical.gamma_rates new_value new_value i
        | Constant      -> Constant,old_model.rate
    in
    { old_model with rate = new_array;
                     spec = { old_model.spec with site_variation = new_spec_var;} }

and update_alpha_invar old_model new_alpha new_invar =
    let new_spec_var,new_array = match old_model.spec.site_variation with
        | Gamma (i,_)   -> Gamma (i,new_alpha),
                           Numerical.gamma_rates new_alpha new_alpha i
        | Theta (i,_,_) -> Theta (i,new_alpha,new_invar),
                           Numerical.gamma_rates new_alpha new_alpha i
        | Constant      -> Constant, old_model.rate
    in
    { old_model with rate = new_array;
                     spec = { old_model.spec with site_variation = new_spec_var; };
                     invar = Some new_invar; }
END

(* [to_formatter m]
 * builds a formatted output of the model [m] spec. Since the spec can be
 * transformed to the model, this is a much more readable form then other model
 * data --the decomposed matrix for example. *)
let to_formatter (model: model) : Xml.xml Sexpr.t list =
    let alph,alph_s = model.spec.alphabet in
    let priors = 
        let inner = 
            Array.mapi
                (fun i v ->
                    (PXML -[Xml.Characters.vector]
                        ([Xml.Alphabet.value] = [`Float v])
                        { `String (Alphabet.match_code i alph)} --))
                (match model.spec.base_priors with
                    | Equal ->
                        Array.make alph_s (1.0 /. (float_of_int alph_s))
                    | Given x | Estimated x -> x)
        in
        (PXML -[Xml.Characters.priors] { `Set (Array.to_list inner) } --)
    and get_model model = match model.spec.substitution with
            | JC69    -> "jc69"     | F81     -> "f81"
            | K2P _   -> "k2p"      | F84 _   -> "f84"
            | HKY85 _ -> "hky85"    | TN93 _  -> "tn93"
            | GTR _   -> "gtr"      | File _  -> "file"
            | Custom _ -> "custom"
    and get_alpha m = match m.spec.site_variation with
        | Constant -> `Float 0.0
        | Gamma (_,x)
        | Theta (_,x,_) -> `Float x
    and get_invar m = match m.spec.site_variation with
        | Gamma _ | Constant -> `Float 0.0
        | Theta (_,_,p) -> `Float p
    and get_gap m = match m.spec.use_gap with
        | `Missing     -> `String "missing"
        | `Independent -> `String "independent"
        | `Coupled x   -> `String ("coupled:"^string_of_float x)
    and get_cats m = match m.spec.site_variation with
        | Constant -> `Int 1
        | Gamma (c,_)
        | Theta (c,_,_) -> `Int c
    and parameters m = match m.spec.substitution with
        | JC69  | F81 -> []
        | K2P x | F84 x | HKY85 x ->
            [(PXML -[Xml.Data.param 1]
                ([Xml.Alphabet.value] = [`Float x])
                { `String "" } --)]
        | TN93 (a,b) -> 
            (PXML -[Xml.Data.param 1]
                ([Xml.Alphabet.value] = [`Float a])
                { `String "" } --) ::
            [(PXML -[Xml.Data.param 2]
                ([Xml.Alphabet.value] = [`Float b])
                { `String "" } --)]
        | GTR ray ->
            let r,_ =
                Array.fold_left
                    (fun (acc,i) x ->
                        (PXML -[Xml.Data.param i]
                            ([Xml.Alphabet.value] = [`Float x])
                            { `String "" } --) :: acc,i+1)
                    ([],1)
                    ray
            in
            List.rev r
        | Custom (_,_,str) ->
            [(PXML -[Xml.Data.param 1]
                ([Xml.Alphabet.value] = [`String str])
                { `String "" } --)]
        | File (_,str) ->
            [(PXML -[Xml.Data.param 1]
                ([Xml.Alphabet.value] = [`String str])
                { `String "" } --)]
    in
    (PXML
        -[Xml.Characters.model]
            ([Xml.Characters.name] = [`String (get_model model)])
            ([Xml.Characters.categories] = [get_cats model])
            ([Xml.Characters.alpha] = [get_alpha model])
            ([Xml.Characters.invar] = [get_invar model]) 
            ([Xml.Characters.gapascharacter] = [get_gap model]) 
        { `Set (priors :: (parameters model)) }
    --) :: []
(* -> Xml.xml Sexpr.t list *)

(** produce a short readable name for the model: like JC69+G+I. We ignore data
    dependent and optimality paramters (gap as missing/independent) and mpl/mal *)
let short_name model =
    let model_name = match model.spec.substitution with
        | JC69    -> "JC69" | F81     -> "F81"
        | K2P _   -> "K81"  | F84 _   -> "F84"
        | HKY85 _ -> "HKY"  | TN93 _  -> "TN93"
        | GTR _   -> "GTR"  | File _  -> "FILE"
        | Custom _-> "CUSTOM"
    and variation_name = match model.spec.site_variation with
        | Theta (i,_,_) when i > 1 -> "+G+I"
        | Gamma _ -> "+G"
        | Theta _ -> "+I"
        | Constant -> ""
    in
    model_name ^ variation_name


(* this assumes that the spec is consistent with the model itself, the returned
 * update function ensures this consistency *)
and get_update_function_for_model model =
  IFDEF USE_LIKELIHOOD THEN
    let split_array ray = 
        ray.((Array.length ray)-1),
        Array.init ((Array.length ray)-1) (fun i -> ray.(i)) in
    match model.spec.substitution,model.spec.use_gap with
        | JC69,`Coupled _   -> Some (fun y x -> update_jc69 y (`Coupled x.(0)))
        | F81,`Coupled _    -> Some (fun y x -> update_f81 y (`Coupled x.(0)))
        | K2P _,`Coupled _  -> Some (fun y x -> update_k2p y x.(0) (`Coupled x.(1)))
        | TN93 _,`Coupled _ -> Some (fun y x -> update_tn93 y (x.(0),x.(1)) (`Coupled x.(2)))
        | F84 _,`Coupled _  -> Some (fun y x -> update_f84 y x.(0) (`Coupled x.(1)))
        | HKY85 _,`Coupled _-> Some (fun y x -> update_hky y x.(0) (`Coupled x.(1)))
        | GTR _,`Coupled _  ->
            Some (fun y x -> 
                    let h,t = split_array x in update_gtr y t (`Coupled h))
        (* no fifth state *)
        | JC69,_ | F81,_ | File _,_ -> None
        | TN93 _,z  -> Some (fun y x -> update_tn93 y (x.(0),x.(1)) z)
        | F84 _,z   -> Some (fun y x -> update_f84 y x.(0) z)
        | GTR _,z   -> Some (fun y x -> update_gtr y x z)
        | K2P _,z   -> Some (fun y x -> update_k2p y x.(0) z)
        | HKY85 _,z -> Some (fun y x -> update_hky y x.(0) z)
        | Custom _,_-> Some (fun y x -> update_custom y x)
  ELSE
    None
  END

and get_update_function_for_alpha model = 
  IFDEF USE_LIKELIHOOD THEN
    match model.spec.site_variation with
        | Gamma _  -> Some update_alpha
        | Theta _  -> Some update_alpha
        | Constant -> None
  ELSE
    None
  END

(* for model parameters *)
let get_current_parameters_for_model model =
  IFDEF USE_LIKELIHOOD THEN
    match model.spec.substitution, model.spec.use_gap with
        | JC69,`Coupled gapr | F81,`Coupled gapr ->
            (Array.make 1 gapr)
        | F84 x,`Coupled gapr | K2P x,`Coupled gapr | HKY85 x,`Coupled gapr ->
            let x = Array.make 2 x in
            x.(1) <- gapr;
            x
        | TN93 (x1,x2),`Coupled gapr ->
            let y = Array.make 3 x2 in
            y.(0) <- x1;
            y.(2) <- gapr;
            y
        | GTR y,`Coupled gapr ->
            let y = Array.init ((Array.length y)+1)
                       (fun i -> if i = Array.length y then gapr else y.(i))
            in
            y
        (* no gap below *)
        | JC69,_ | F81,_ | File _,_ -> [||]
        | F84 x,_ | K2P x,_ | HKY85 x,_ ->
            Array.make 1 x
        | TN93 (x1,x2),_  ->
            let y = Array.make 2 x2 in
            y.(0) <- x1;
            y
        | GTR x,_ -> x
        | Custom (_,xs,_),_ -> xs
  ELSE
    [||]
  END

(* for alpha parameter *)
and get_current_parameters_for_alpha model =
  IFDEF USE_LIKELIHOOD THEN
    match model.spec.site_variation with
        | Gamma (cat,alpha) -> Some alpha
        | Theta (cat,alpha,invar) -> Some alpha
        | Constant -> None
  ELSE
    None
  END
