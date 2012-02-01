(* POY 5.0 Alpha. A phylogenetic analysis program using Dynamic Homologies.   *)
(* Copyright (C) 2011 Andrés Varón, Lin Hong, Nicholas Lucaroni, Ward Wheeler,*)
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

let () = SadmanOutput.register "MlMatrix" "$Revision"

let (=.) a b = abs_float (a-.b) < Numerical.tolerance
let (-->) b a = a b
let failwithf format = Printf.ksprintf failwith format

let likelihood_not_enabled =
    "Likelihood@ not@ enabled:@ download@ different@ binary@ or@ contact@ mailing@ list" 
let dyno_likelihood_warning = 
    "Gap@ as@ an@ additional@ character@ is@ required@ for@ the@ dynamic@ "^
    "likelihood@ criteria.@ I@ am@ enabling@ this@ setting@ for@ the@ transformation."
let dyno_gamma_warning = 
    "Gamma@ classes@ for@ dynamic@ MPL/FLK@ are@ un-necessary,@ and@ are@ being@ removed."

let debug = false
let debug_printf msg format =
    Printf.ksprintf (fun x -> if debug then print_string x; flush stdout) msg format

let ba_of_array1 x = Bigarray.Array1.of_array Bigarray.float64 Bigarray.c_layout x
and ba_of_array2 x = Bigarray.Array2.of_array Bigarray.float64 Bigarray.c_layout x

let create_ba1 x     = Bigarray.Array1.create Bigarray.float64 Bigarray.c_layout x
and create_ba2 x y   = Bigarray.Array2.create Bigarray.float64 Bigarray.c_layout x y

let barray_to_array2 bray =
    let a = Bigarray.Array2.dim1 bray and b = Bigarray.Array2.dim2 bray in
    let r = Array.make_matrix a b 0.0 in
    for i = 0 to a-1 do for j = 0 to b-1 do
        r.(i).(j) <- bray.{i,j};
    done; done; r

let print_barray1 a =
    for i = 0 to (Bigarray.Array1.dim a)-1 do
        Printf.printf "%2.10f\t" a.{i};
    done; Printf.printf "\n"; ()

and print_barray2 a =
    for i = 0 to (Bigarray.Array2.dim1 a)-1 do
        for j = 0 to (Bigarray.Array2.dim2 a)-1 do
            Printf.printf "%2.10f\t" a.{i,j};
        done; Printf.printf "\n";
    done; Printf.printf "\n"; ()

(* type to help the parsing of specification data *)
type string_spec = string * (string * string * string * string) * float list
                 * [ `Given of (string * float) list | `Estimate of float array option
                   | `Consistent of float array option | `Equal ]
                 * (string * float option) * string * string option
let empty_str_spec : string_spec = ("",("","","",""),[],`Consistent None,("",None),"",None)

(* the default when commands are left off from the interactive console / scripts *)
let default_command = (`MAL,`GTR None,None,`Consistent,`Missing)

(* --- DEFAULTS FOR MODELS FROM PHYML --- *)
(*  These are used for DNA sequences, and
 *  might not be applicable for all types
 *  of datasets.                          *)
(* See, Set_Model_Parameters,
 *      Init_Model
 *      Translate_Custom_Mod_String       *)
let default_tstv    = 4.0
let default_gtr a   = Array.make (((a-1)*a)/2) 1.0
let default_alpha p = if p then 0.2 else 1.0
let default_invar   = 0.2
let default_gap_r   = 0.15

(* list of set bits, or packed integer of set bits *)
type chars = [ `List of int list | `Packed of int ]

type site_var = 
    | Gamma of int * float
    | Theta of int * float * float
    (* | Given of (float * float) array *)
    | Constant 

type subst_model =
    | JC69
    | F81
    | K2P    of float option
    | F84    of float option
    | HKY85  of float option
    | TN93   of (float * float) option
    | GTR    of float array option
    | File   of float array array * string
    | Custom of (int All_sets.IntegerMap.t * float array * string)

type priors = 
    | Estimated of float array
    | Given     of float array
    | Equal

type spec = {
    substitution : subst_model;
    site_variation : site_var option;
    base_priors : priors;
    cost_fn : Methods.ml_costfn;
    use_gap : Methods.ml_gap;
    iterate_model : bool;
    iterate_alpha : bool;
}

type model = {
    spec  : spec;
    alph  : Alphabet.a;
    alph_s: int;
    pi_0  : (float, Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array1.t;
    invar : float option;
    rate  : (float, Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array1.t;
    prob  : (float, Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array1.t;
    s     : (float, Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array2.t;
    u     : (float, Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array2.t;
    d     : (float, Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array2.t;
    ui    : (float, Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array2.t option; 
}

let jc69_5 = { substitution = JC69; site_variation= None;cost_fn=`MAL;
               iterate_alpha=true;iterate_model=true;use_gap=`Independent;
               base_priors = Equal;  }

and jc69_4 = { substitution = JC69;site_variation= None;cost_fn=`MAL;
               iterate_alpha=true;iterate_model=true;use_gap=`Missing;
               base_priors=Equal; }

let jc69_5_gap flo = { jc69_5 with use_gap = `Coupled flo; }

module OrderedML = struct
    (* we could choose model or spec, but spec can use Pervasives.compare *)
    type t = spec 
    let compare a b = Pervasives.compare a b
end
module MlModelMap = Map.Make (OrderedML)

let get_costfn_code a = match a.spec.cost_fn with 
    | `MPL -> 1 
    | `MAL -> 0 
    | `FLK -> ~-1 (* should not call C functions; yet *)

let categorize_by_model get_fn codes =
    let set_codes =
        List.fold_left
            (fun acc code ->
                try let spec = get_fn code in
                    try let old = MlModelMap.find spec acc in
                        MlModelMap.add spec (code::old) acc
                    with | Not_found ->
                        MlModelMap.add spec ([code]) acc
                with | _ -> acc)
            MlModelMap.empty
            codes
    in
    MlModelMap.fold (fun _ e a -> e :: a) set_codes []

let compare_priors a b =
    let compare_array x y = 
        let results = ref true in
        for i = 0 to a.alph_s - 1 do
            results := !results && (x.(i) =. y.(i));
        done;
        !results
    in
    match a.spec.base_priors ,a.spec.base_priors with
        | _ when a.alph_s != b.alph_s -> false
        | Equal , Equal            -> true
        | Estimated x, Estimated y
        | Given     x, Given     y -> compare_array x y
        | (Equal | Estimated _ | Given _), _ -> false

IFDEF USE_LIKELIHOOD THEN

(* a gentler compare that excludes the parameters of the model itself
    Not to be used for a total ordering *)
let compare a b = 
    let m_compare = match a.spec.substitution,b.spec.substitution with
        | JC69 , JC69 | F81 , F81 -> 0
        | K2P m, K2P x | F84 m, F84 x | HKY85 m, HKY85 x ->
                begin match m,x with
                | None,None -> 0
                | Some v,Some x when v =. x-> 0
                | _,_ -> ~-1
                end
        | Custom _, Custom _
        | TN93 _,  TN93 _ | GTR _, GTR _ -> 0
        | File (_,x), File (_,y) when x = y -> 0
        | _,_ -> ~-1
    and c_compare = if a.spec.cost_fn = b.spec.cost_fn then 0 else ~-1
    and v_compare = match a.spec.site_variation,b.spec.site_variation with
        | Some (Gamma (i,a)), Some (Gamma (ix,ax)) ->
                if i = ix && a = ax then 0 else ~-1
        | Some (Theta (i,a,p)), Some (Theta (ix,ax,px)) ->
                if i = ix && a = ax && p = px then 0 else ~-1
        | None, None | Some Constant, Some Constant 
        | Some Constant, None | None, Some Constant -> 0
        | _,_ -> ~-1
    and g_compare = match a.spec.use_gap,b.spec.use_gap with
        | `Missing, `Missing
        | `Independent, `Independent -> 0
        | `Coupled x, `Coupled y when x = y -> 0
        | _, _ -> ~-1
    and p_compare = if compare_priors a b then 0 else ~-1 in
    (* just knowing that they are different is enough *)
    if debug then begin
        Printf.printf "Models : %b\n%!" (m_compare >= 0);
        Printf.printf "Rate   : %b\n%!" (v_compare >= 0);
        Printf.printf "Use Gap: %b\n%!" (g_compare >= 0);
        Printf.printf "Priors : %b\n%!" (p_compare >= 0);
    end;
    m_compare + v_compare + g_compare + p_compare + c_compare
    
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
    let srm = create_ba2 a_size a_size in
    Bigarray.Array2.fill srm beta;
    (* modify transition elements to alpha *)
    srm.{1, 3} <- alpha; srm.{3, 1} <- alpha;
    srm.{2, 0} <- alpha; srm.{0, 2} <- alpha;
    (* set up the diagonal elements *)
    let () = match gap_r with
        | None -> 
            let diag = if a_size = 4 then
                            -. alpha -. beta -. beta
                       else -. alpha -. beta *. 3.0  in
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
    let srm = create_ba2 a_size a_size in
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
    let () = match gap_r with
        | None -> 
            for i = 0 to (a_size-1) do
                for j = 0 to (a_size-1) do
                    srm.{i,j} <- 
                        if i = j then -. lambda *. (1.0 -. pi_.{i})
                                 else pi_.{j} *. lambda;
                done;
            done;
        | Some (k,r) ->
            for i = 0 to (a_size-1) do
                for j = 0 to (a_size-1) do
                    srm.{i,j} <- 
                        if i = j then -. lambda *. (1.0 -. pi_.{i})
                                 else pi_.{j} *. lambda;
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
    let y = pi_.{1} +. pi_.{3} in (* Y = C + T *)
    let r = pi_.{0} +. pi_.{2} in (* R = A + G *)
    let alpha = (1.0+.kappa/.r) *. gamma in
    let beta = (1.0+.kappa/.y) *. gamma in
    m_tn93 pi_ alpha beta gamma a_size gap_r

(* normalize against two characters that are not gaps; unless its the only choice *)
let normalize ?(m=Numerical.minimum) spec ray = match spec.use_gap with
    | `Independent -> 
        let normalize_factor = max m ray.((Array.length ray) - 3) in
        let vec = Array.map (fun i -> (max m i) /. normalize_factor) ray in
        vec, `Independent
    | `Coupled n ->
        let normalize_factor = max m ray.((Array.length ray) - 1) in
        let vec = Array.map (fun i -> (max m i) /. normalize_factor) ray in
        vec, `Coupled ( n/.normalize_factor )
    | `Missing ->
        let normalize_factor = max m ray.((Array.length ray) - 1) in
        let vec = Array.map (fun i -> (max m i) /. normalize_factor) ray in
        vec,`Missing

(* val gtr :: ANY ALPHABET size
   pi_ == n; co_ == ((n-1)*n)/2
   form of lower packed storage mode, excluding diagonal, *)
let m_gtr_independent pi_ co_ a_size =
    if ((a_size*(a_size-1))/2) <> Array.length co_ then begin
        failwithf "Length of GTR parameters is insufficient for alphabet";
    end;
    (* last element of GTR = 1.0 *)
    let last = co_.( (Array.length co_) - 1) in
    let co_ = Array.map (fun x -> x /. last) co_ in
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
    (* normalize diagonal so row sums to 0 *)
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
    if (((a_size-1)*(a_size-2))/2) <> Array.length co_ then begin
        failwithf "Length of GTR parameters (%d) is insufficient for alphabet with coupled gap parameter (%d)"
            (Array.length co_) (((a_size-1)*(a_size-2))/2);
    end;
     (* last element of GTR = 1.0 *)
    let last = co_.( (Array.length co_) - 1) in
    let co_ = Array.map (fun x -> x /. last) co_ in
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
    (* normalize diagonal so row sums to 0 *)
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

let m_gtr pi_ co_ a_size gap_r = match gap_r with
    | Some (i,r) -> m_gtr_coupled pi_ co_ a_size i r
    | None       -> m_gtr_independent pi_ co_ a_size

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
    if not (a_size = alph_size) then
        failwithf "Alphabet size for priors and model inconsistent: %d =/= %d" a_size alph_size;
    for i = 0 to a_size -1 do
        assert( Array.length f_aa.(i) = a_size );
        for j = i+1 to a_size-1 do
            let letter = Char.code f_aa.(i).(j) in
            if not ( letter = Char.code f_aa.(j).(i) ) then
                failwithf "Custom model should be symmetric: %d,%d =/= %d,%d" i j j i;
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

(* functions to diagonalize the two types of substitution matrices *)
let diagonalize sym mat =
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

(* for pure ocaml median functions *)
let compose model t = match model.ui with
    | Some ui -> compose_gtr FMatrix.scratch_space model.u model.d ui t
    | None    -> compose_sym FMatrix.scratch_space model.u model.d t

let subst_matrix model topt =
    let _gapr = match model.spec.use_gap with
        | `Coupled x   -> Some (Alphabet.get_gap model.alph, x)
        | `Independent -> None
        | `Missing     -> None
    and a_size = model.alph_s 
    and priors = model.pi_0 in
    let m = match model.spec.substitution with
        | JC69  -> m_jc69 priors 1.0 a_size _gapr
        | F81   -> m_f81 priors 1.0 a_size _gapr
        | K2P t -> let t = match t with | Some x -> x | None -> default_tstv in
                   m_k2p priors t 1.0 a_size _gapr
        | F84 t -> let t = match t with | Some x -> x | None -> default_tstv in
                   m_f84 priors t 1.0 a_size _gapr
        | HKY85 t -> let t = match t with | Some x -> x | None -> default_tstv in
                   m_hky85 priors t a_size _gapr
        | TN93 tstv -> let ts,tv = match tstv with | Some (x,y) -> x,y | None -> default_tstv,default_tstv in
                   m_tn93 priors ts tv 1.0 a_size _gapr
        | GTR c ->
                let c,gapr = match c,_gapr with
                    | Some xs, _   -> normalize model.spec xs
                    | None, Some _ -> default_gtr (a_size-1), model.spec.use_gap
                    | None, _      -> default_gtr a_size, model.spec.use_gap
                in
                m_gtr priors c a_size _gapr
        | File (m,s)-> m_file priors m a_size
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
    | None -> m

(* print output in our nexus format or Phyml output *)
let output_model output nexus model set = 
    let printf format = Printf.ksprintf output format in
    let gtr_mod = ref false in
    if nexus = `Nexus then begin
        printf "@[Likelihood@.";
        let () = match model.spec.substitution with
            | JC69   -> printf "@[Model = JC69;@]@\n";
            | F81    -> printf "@[Model = F81;@]";
            | K2P None -> printf "@[Model = K2P;@]";
            | K2P (Some x) -> printf "@[Model = K2P;@]";
                        printf "@[Parameters = %f;@]" x
            | F84 None -> printf "@[Model = F84;@]";
            | F84 (Some x) -> printf "@[Model = F84;@]";
                        printf "@[Parameters = %f;@]" x
            | HKY85 None -> printf "@[Model = HKY85;@]";
            | HKY85 (Some x) -> printf "@[Model = HKY85;@]";
                        printf "@[Parameters = %f;@]" x
            | TN93 None -> printf "@[Model = TN93;@]";
            | TN93 (Some (a,b)) -> printf "@[Model = TN93;@]";
                        printf "@[Parameters = %f %f;@]" a b
            | GTR None -> printf "@[Model = GTR;@]"
            | GTR (Some xs) -> printf "@[Model = GTR;@]";
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
                    (Alphabet.to_list model.alph);
                printf ";@]@]"
        in
        let () = match model.spec.cost_fn with
            | `MPL -> printf "@[Cost = mpl;@]";
            | `MAL -> printf "@[Cost = mal;@]"; 
            | `FLK -> printf "@[Cost = flk;@]";
        in
        let () = match model.spec.site_variation with
            | Some Constant | None -> ()
            | Some (Gamma (c, p)) ->
                printf "@[Variation = gamma;@]@[alpha = %.5f@]@[sites = %d;@]" p c
            | Some (Theta (c, p, i)) ->
                printf ("@[Variation = gamma;@]@[alpha = %.5f@]@[sites = %d;@]"
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
            | Some Constant
            | None -> printf "No@]@\n"
            | Some (Gamma (cats,param)) ->
                printf ("Yes@]@\n@[<hov 1>- Number of categories: %d@]"^^
                        "@[<hov 1>- Gamma Shape Parameter: %.4f@]") cats param
            | Some (Theta (cats,param,inv)) ->
                printf ("Yes@]@\n[@<hov 1>- Number of categories: %d@]"^^
                        "@[<hov 1>- Gamma Shape Parameter: %.4f@]") cats param;
                printf ("@[<hov 1>- Proportion of invariant: %.4f@]") inv
        in
        printf "@]";
        let () = match model.spec.cost_fn with
            | `MPL -> printf "@[<hov 0>Cost mode: mpl;@]@\n";
            | `MAL -> printf "@[<hov 0>Cost mode: mal;@]@\n"; 
            | `FLK -> printf "@[<hov 0>Cost mode: flk;@]@\n"; 
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
                    (Alphabet.to_list model.alph);
        in
        printf "@[<hov 0>Model Parameters: ";
        let () = match model.spec.substitution with
            | JC69  -> printf "JC69@]@\n"
            | F81   -> printf "F81@]@\n"
            | K2P x -> printf "K2P@]@\n@[<hov 1>- Transition/transversion ratio: %.5f@]@\n"
                        (match x with Some x -> x | None -> default_tstv)
            | F84 x -> printf "F84@]@\n@[<hov 1>- Transition/transversion ratio: %.5f@]@\n"
                        (match x with Some x -> x | None -> default_tstv)
            | HKY85 x->printf "HKY85@]@\n@[<hov 1>- Transition/transversion ratio:%.5f@]@\n"
                        (match x with Some x -> x | None -> default_tstv)
            | TN93 None->printf "tn93@]@\n@[<hov 1>- transition/transversion ratio:%.5f/%.5f@]@\n"
                            default_tstv default_tstv
            | TN93 (Some (a,b)) -> 
                printf "tn93@]@\n@[<hov 1>- transition/transversion ratio:%.5f/%.5f@]@\n" a b 
            | GTR x -> gtr_mod := true;
                printf "GTR@]@\n@[<hov 1>- Rate Parameters: @]@\n";
                let get_str i = Alphabet.match_code i model.alph
                and convert s r c = (c + (r * (s-1)) - ((r*(r+1))/2)) - 1 in
                let ray = match x with | Some x -> x | None -> default_gtr model.alph_s in
                begin match model.spec.use_gap with
                    | `Coupled x ->
                        for i = 0 to model.alph_s - 2 do for j = i+1 to model.alph_s - 2 do
                            printf "@[<hov 1>%s <-> %s - %.5f@]@\n" (get_str i) 
                                    (get_str j) ray.(convert (model.alph_s-1) i j)
                        done; done;
                        printf "@[<hov 1>%s <-> N - %.5f@]@\n" 
                            (get_str (Alphabet.get_gap model.alph)) x
                    | `Missing | `Independent ->
                        for i = 0 to model.alph_s - 1 do for j = i+1 to model.alph_s - 1 do
                            printf "@[<hov 1>%s <-> %s - %.5f@]@\n" (get_str i) 
                                    (get_str j) ray.(convert model.alph_s i j)
                        done; done
                end
            | File (ray,name) ->
                printf "File:%s@]@\n" name;
                let mat = compose model 0.0 in
                printf "@[<hov 1>[";
                for i = 0 to model.alph_s - 1 do
                    printf "%s ------- " (Alphabet.match_code i model.alph)
                done;
                printf "]@]@\n";
                for i = 0 to model.alph_s - 1 do 
                    for j = 0 to model.alph_s - 1 do
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
        let () = 
            let mat = compose model ~-.1.0 in
            printf "@[<hov 1>[";
            for i = 0 to model.alph_s - 1 do
                printf "%s ------- " (Alphabet.match_code i model.alph)
            done;
            printf "]";
            for i = 0 to model.alph_s - 1 do 
                printf "@]@\n@[<hov 1>";
                for j = 0 to model.alph_s - 1 do
                    printf "%8.5f  " mat.{i,j}
                done; 
            done;
        in 
        printf "@]@\n@]"
    end

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

(* check the metricity of the composed matrix by issuing a passive warning *)
let check_metricity model t1 t2 log_mat1 log_mat2 : unit =
    let flag = ref false in
    for i = 0 to model.alph_s-1 do
        for j = 0 to model.alph_s-1 do
            let a_state = ref [] and a_cost = ref infinity in
            for n = 0 to model.alph_s-1 do
                let n_cost = log_mat1.{i,n} +. log_mat2.{j,n} in
                if !a_cost =. n_cost then begin
                    a_state := n :: !a_state;
                    a_cost  := min !a_cost n_cost;
                end else if !a_cost > n_cost then begin
                    a_state := [n];
                    a_cost  := n_cost;
                end;
            done;
            if not ((List.mem i !a_state) || (List.mem j !a_state)) then
                flag := true;
        done;
    done;
    let gap_ratio = match model.spec.use_gap with
        | `Coupled f -> string_of_float f
        | `Independent | `Missing -> ""
    in
    if !flag then begin
        Status.user_message Status.Warning 
            (Printf.sprintf "Matrix has metric issues: %f:%f/%s" t1 t2 gap_ratio);
(*        Status.user_message Status.Warning (Printf.sprintf "Branch Length: %f/%f\n" t1 t2);*)
(*        output_model (Status.user_message Status.Warning) `Nexus model None;*)
    end;
    ()

(* integerized converstion and composition of a model and branch lenghth *)
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
    and imatrix = Array.make_matrix model.alph_s model.alph_s 0 in
    assert( model.alph_s = Array.length matrix);
    assert( model.alph_s = Array.length matrix.(0));
    for i = 0 to (Array.length matrix) - 1 do
        for j = 0 to (Array.length matrix.(0)) - 1 do
            imatrix.(i).(j) <- create sigma matrix i j
        done;
    done;
    imatrix

ELSE

    let output_model _ _ _ _ = failwith likelihood_not_enabled
    let compose _ _  = failwith likelihood_not_enabled
    let spec_from_classification _ _ _ _ _ _ = failwith likelihood_not_enabled
    let compare _ _ = failwith likelihood_not_enabled
    let classify_seq_pairs _ _ _ _ _ = failwith likelihood_not_enabled
    let subst_matrix _ _ = failwith likelihood_not_enabled
    let check_metricity _ _ _ _ _ = failwith likelihood_not_enabled
    let process_custom_model _ = failwith likelihood_not_enabled

    (** durability functions for outside toplevel access **)
    let compose_model _ _ = failwith likelihood_not_enabled
    let m_file _ _ _ = failwith likelihood_not_enabled
    let m_jc69 _ _ _ = failwith likelihood_not_enabled
    let m_gtr  _ _ _ = failwith likelihood_not_enabled

END

(* create a cost matrix from a model *)
let model_to_cm model t =
    let input =
        IFDEF USE_LIKELIHOOD THEN 
            let t = max Numerical.minimum t in
            integerized_model model t
        ELSE
            failwith likelihood_not_enabled
        END
    in
    let llst = Array.to_list (Array.map Array.to_list input) in
    let res = Cost_matrix.Two_D.of_list ~suppress:true llst model.alph_s in
    res

(* ------------------------------------------------ *)
(* CONVERSION/MODEL CREATION FUNCTIONS             *)

(* convert a string spec to a specification, used in Parser for nexus *)
let convert_string_spec ((name,(var,site,alpha,invar),param,priors,gap,cost,file):string_spec) =
  IFDEF USE_LIKELIHOOD THEN
    let iterate_model = ref true in
    let iterate_alpha = ref true in
    let submatrix = match String.uppercase name with
        | "JC69" -> 
            iterate_model := false;
            (match param with
            | [] -> JC69
            | _ -> failwith "Parameters don't match model")
        | "F81" ->
            iterate_model := false;
            (match param with
            | [] -> F81
            | _ -> failwith "Parameters don't match model")
        | "K80" | "K2P" -> (match param with
            | ratio::[] ->
                iterate_model := false;
                K2P (Some ratio)
            | []        -> K2P None
            | _ -> failwith "Parameters don't match model")
        | "F84" -> (match param with
            | ratio::[] ->
                iterate_model := false;
                F84 (Some ratio)
            | []        -> F84 None
            | _ -> failwith "Parameters don't match model")
        | "HKY" | "HKY85" -> (match param with
            | ratio::[] ->
                iterate_model := false;
                HKY85 (Some ratio)
            | []        -> HKY85 None
            | _ -> failwith "Parameters don't match model")
        | "TN93" -> (match param with
            | ts::tv::[] ->
                iterate_model := false;
                TN93 (Some (ts,tv))
            | []         -> TN93 None
            | _ -> failwith "Parameters don't match model")
        | "GTR" -> (match param with
            | [] -> GTR None
            | ls -> 
                iterate_model := false;    
                GTR (Some (Array.of_list ls)) )
        | "GIVEN"->
            begin match file with
            | Some name ->
                iterate_model := false;    
                name
                    --> FileStream.read_floatmatrix
                    --> List.map (Array.of_list)
                    --> Array.of_list
                    --> (fun x -> File (x,name))
            | None -> failwith "File not specified for Likelihood Model."
            end
        | "CUSTOM" ->
            begin match file with
            | Some name ->
                iterate_model := true;
                failwith "not done"
            | None -> failwith "File not specified for Likelihood Model."
            end
        (* ERROR *)
        | "" -> failwith "No Model Specified"
        | _  -> failwith "Incorrect Model"
    in
    let cost_fn = match String.uppercase cost with
        | "MPL" -> `MPL
        | "MAL" -> `MAL
        | "FLK" -> `FLK
        | ""    -> `MAL (* unmentioned default *)
        | x     -> 
            failwith ("I don't know "^x^" as a cost mode for likelihood.")
    and variation = match String.uppercase var with
        | "GAMMA" ->
            let alpha = try let res = float_of_string alpha in
                            iterate_alpha := false;
                            res
                        with | _ -> default_alpha false in
            Some (Gamma (int_of_string site, alpha))
        | "THETA" -> 
            let alpha = try let res = float_of_string alpha in
                            iterate_alpha := false;
                            res
                        with | _ -> default_alpha true in
            let invar = try let res = float_of_string invar in
                            iterate_alpha := !iterate_alpha || false;
                            res
                        with | _ -> default_invar in
            Some (Theta (int_of_string site, alpha, invar))
        | "NONE" | "CONSTANT" | "" ->
            iterate_alpha := false;
            None
        | _ -> failwith "Unrecognized variation method"
    and priors = match priors with
        | `Given priors  -> Given (Array.of_list (List.map snd priors))
        | `Equal -> Equal
        | `Estimate (Some x) -> Estimated x
        | `Consistent pre_calc ->
            begin match submatrix, pre_calc with
                | JC69, _
                | K2P _, _   -> Equal
                | _, Some pi -> Estimated pi
                | _ , None   -> assert false
            end
        | `Estimate None -> assert false
    and gap_info =
        match String.uppercase (fst gap),snd gap with
        | "COUPLED", None   -> `Coupled default_gap_r
        | "COUPLED", Some x -> `Coupled x
        | "INDEPENDENT",_   -> `Independent
        | "MISSING",_       -> `Missing
        | "",_              -> `Missing
        | x,_               -> failwith ("Invalid gap property: "^x)
    in
    {   substitution = submatrix;
        site_variation = variation;
        base_priors = priors;
        use_gap = gap_info;
        cost_fn = cost_fn;
        iterate_model = !iterate_model;
        iterate_alpha = !iterate_alpha; }
  ELSE
    failwith likelihood_not_enabled
  END


let convert_methods_spec alph_size compute_priors (_,cst,subst,site_variation,base_priors,use_gap) =
    let u_gap = match use_gap with 
            | `Independent | `Coupled _ -> true | `Missing -> false in
    let alph_size = if u_gap then alph_size else alph_size - 1 in
    let i_alpha = ref true and i_model = ref true in
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
        | None   -> Some Constant 
        | Some (`Gamma (w,y)) -> 
            let y = match y with
                | Some x -> i_alpha := false; x 
                | None   -> default_alpha false
            in
            Some (Gamma (w,y))
        | Some (`Theta (w,y)) ->
            let y,p = match y with 
                | Some x -> i_alpha := false; x 
                | None   -> default_alpha true, default_invar
            in
            Some (Theta (w,y,p))
    and substitution = match subst with
        | `JC69 -> i_model := false;
                   JC69
        | `F81  -> i_model := false;
                   F81
        | `K2P (Some x) ->
            let aray = Array.of_list x in
            if Array.length aray = 1 then begin
                i_model := false;
                K2P (Some aray.(0))
            end else if Array.length aray = 0 then begin
                K2P None
            end else
                let () =
                    Status.user_message Status.Error
                        "Likelihood@ model@ K80@ requires@ 1@ or@ 0@ parameters"
                in
                failwith("Incorrect Parameters");
        | `K2P None -> K2P None
        | `HKY85 (Some x) ->
            let aray = Array.of_list x in
            if Array.length aray = 1 then begin
                i_model := false;
                HKY85 (Some aray.(0))
            end else if Array.length aray = 0 then begin
                HKY85 None
            end else
                let () =
                    Status.user_message Status.Error
                        "Likelihood@ model@ HKY85@ requires@ 1@ or@ 0@ parameters"
                in
                failwith("Incorrect Parameters");
        | `HKY85 None -> HKY85 None
        | `F84 (Some x) ->
            let aray = Array.of_list x in
            if Array.length aray = 1 then begin
                i_model := false;
                F84 (Some aray.(0))
            end else if Array.length aray = 0 then begin
                F84 None
            end else
                let () =
                    Status.user_message Status.Error
                        "Likelihood@ model@ F84@ requires@ 1@ parameters"
                in
                failwith("Incorrect Parameters");
        | `F84 None -> F84 None
        | `TN93 (Some x) ->
            let aray = Array.of_list x in
            if Array.length aray <> 2 then
                let () =
                    Status.user_message Status.Error
                        "Likelihood@ model@ TN93@ requires@ 2@ or@ 0@ parameters"
                in
                failwith("Incorrect Parameters");
            else if Array.length aray = 0 then begin
                TN93 None
            end else begin
                i_model := false;
                TN93 (Some (aray.(0),aray.(1)))
            end
        | `TN93 None -> TN93 None
        | `GTR (Some x) ->
            let aray = Array.of_list x in
            let n_a = (alph_size * (alph_size-1)) / 2 in
            if (Array.length aray) <> n_a then
                let () = Status.user_message Status.Error
                    ("Likelihood@ model@ GTR@ requires@ (a-1)*(a/2)@ "^
                     "parameters@ with@ alphabet@ size@ a. In@ this@ case@ "^
                     (string_of_int n_a) ^",@ with@ a@ =@ "^ (string_of_int alph_size) ^".")
                in
                failwith "Incorrect Parameters";
            else if Array.length aray = 0 then begin
                GTR None
            end else begin
                i_model := false;
                GTR (Some aray)
            end
        | `GTR None -> GTR None
        | `File str ->
                (* this needs to be changed to allow remote files as well *)
            let matrix = Cost_matrix.Two_D.matrix_of_file float_of_string (`Local str) in
            let matrix = Array.of_list (List.map (Array.of_list) matrix) in
            (* check the array size == a_size *)
            (* check the array array size == a_size *)
            Array.iter 
                (fun x ->
                    if Array.length x = alph_size then ()
                    else failwith "Input@ file@ is@ inconsistent@ with@ alphabet.")
                (matrix);
            if Array.length matrix = alph_size then
                File (matrix,str)
            else
                failwith "Input@ file@ is@ inconsistent@ with@ alphabet."
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
     iterate_model = !i_model;
     iterate_alpha = !i_alpha;
           cost_fn = cst;
           use_gap = use_gap; }


(* check the rates so SUM(r_k*p_k) == 1 and SUM(p_k) == 1 |p| == |r| *)
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

(* create a model based on a specification and an alphabet *)
let create ?(min_prior=Numerical.minimum) alph lk_spec = 
  IFDEF USE_LIKELIHOOD THEN
    let alph = Alphabet.to_sequential alph in
    let (a_size) = match lk_spec.use_gap with
        | `Missing -> (Alphabet.size alph) - 1
        | `Independent 
        | `Coupled _ -> Alphabet.size alph
    in
    (* set up all the probability and rates *)
    let variation,probabilities,invar =
        match lk_spec.site_variation with
        | None ->
            ba_of_array1 [| 1.0 |], ba_of_array1 [| 1.0 |], None
        | Some a -> begin match a with
            | Constant ->  (* same as NONE *)
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
            failwith "MlModel.create: Priors don't match alphabet"
        end;
        ba_of_array1 p
    in
    let _gapr = match lk_spec.use_gap with
        | `Coupled x   -> Some (Alphabet.get_gap alph, x)
        | `Independent -> None
        | `Missing     -> None
    in
    (*  get the substitution rate matrix and set sym variable and to_formatter vars *)
    let sym, sub_mat, subst_model, rep_gap = match lk_spec.substitution with
        | JC69  -> true,  m_jc69 priors 1.0 a_size _gapr, JC69, lk_spec.use_gap
        | F81   -> false, m_f81 priors 1.0 a_size _gapr, F81, lk_spec.use_gap
        | K2P t -> 
            let t = match t with | Some x -> x | None -> default_tstv in
            true,  m_k2p priors t 1.0 a_size _gapr, lk_spec.substitution, lk_spec.use_gap
        | F84 t -> 
            let t = match t with | Some x -> x | None -> default_tstv in
            false, m_f84 priors t 1.0 a_size _gapr, lk_spec.substitution, lk_spec.use_gap  
        | HKY85 t ->
            let t = match t with | Some x -> x | None -> default_tstv in
            false, m_hky85 priors t a_size _gapr, lk_spec.substitution, lk_spec.use_gap
        (* more complex models *)
        | TN93 tstv ->
            let ts,tv = match tstv with 
                | Some (x,y) -> x,y 
                | None       -> default_tstv,default_tstv
            in
            false, m_tn93 priors ts tv 1.0 a_size _gapr, lk_spec.substitution, lk_spec.use_gap
        | GTR c ->
            let c,gapr = match c,_gapr with 
                | Some xs, _   -> normalize lk_spec xs
                | None, Some _ -> default_gtr (a_size-1), lk_spec.use_gap 
                | None, _      -> default_gtr a_size, lk_spec.use_gap
            in
            false, m_gtr priors c a_size _gapr, (GTR (Some c)),gapr
        | File (m,s)->
            false, m_file priors m a_size, lk_spec.substitution, lk_spec.use_gap
        | Custom (assoc,xs,_) ->
            false, m_custom priors assoc xs a_size, lk_spec.substitution, lk_spec.use_gap
    in
    (* ensure that when priors are not =, we use a model that asserts that *)
    let lk_spec = match lk_spec.base_priors with
        | Estimated _  when sym ->
            Status.user_message Status.Warning 
                "I am using equal priors as required in the selected likelihood model";
            { lk_spec with base_priors = Equal }
        | Given _ when sym ->
            Status.user_message Status.Warning
                "I am using equal priors as required in the selected likelihood model";
            { lk_spec with base_priors = Equal }
        | (Estimated _ | Equal | Given _) -> lk_spec
    in
    let (u_,d_,ui_) = diagonalize sym sub_mat in
    {
      rate = variation;
      prob = probabilities;
      spec = {lk_spec with substitution = subst_model; 
                           use_gap = rep_gap; };
     invar = invar;
      pi_0 = priors;
      alph = alph;
    alph_s = a_size;
         s = sub_mat;
         u = u_;
         d = d_;
        ui = ui_;
    }
  ELSE
    failwith likelihood_not_enabled
  END

let debug_model model t =
    let subst = subst_matrix model t in
    print_barray2 subst; print_newline ();
    match t with | Some t -> print_barray2 (compose model t) | None -> ()

let replace_priors model array = 
    if debug then begin
        Printf.printf "Replacing Priors\n\t%!";
        Array.iter (fun x -> Printf.printf "%f, " x) array;
        print_newline ();
    end;
    create model.alph {model.spec with base_priors = Estimated array}

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

let add_gap_to_model compute_priors model = 
    Status.user_message Status.Warning dyno_likelihood_warning;
    let size = model.alph_s + 1 in
    let priors = match model.spec.base_priors with
        | Estimated _  -> Estimated (compute_priors ())
        | Equal        -> Equal
        | Given    _   ->
            failwith ("I cannot transform the specified model to add gap as a"^
                      " character. The given priors requires a prior for the "^
                      "gap character.")
    and rates = match model.spec.substitution with
        (* these models require no changes *)
        | JC69    | F81 | K2P _   | F84 _
        | HKY85 _ | TN93 _ | GTR None -> model.spec.substitution
        (* user defined rate matrices cannot on transfered *)
        | Custom _
        | File _ -> 
            failwith ("I cannot transform the specified characters to"^
                      " dynamic likelihood characters; the given rate"^
                      " matrix requires gap transformation rates.")
        | GTR (Some xs) ->
            let convert_i_to_rc l i = 
                let rec convert_ r v =
                    let sub = l - r - 1 in
                    if v < sub then (r,l+(v-sub))
                    else convert_ (r+1) (v-sub)
                in
                convert_ 0 i
            and convert_rc_to_i l (r,c) = 
                let r,c = (min r c)+1,(max r c)+1 in
                c + (l*r) - l - 1 - (((r+1)*r)/2)
            in
            let ngtr =
                Array.init 
                    (((size-1)*size)/2)
                    (fun i -> 
                        let r,c = convert_i_to_rc size i in
                        if c = size - 1 then 0.01
                        else xs.(convert_rc_to_i (size-1) (r,c)))
            in
            GTR (Some ngtr)
    in
    let new_spec = {model.spec with base_priors = priors; 
                                        use_gap = `Independent;
                                   substitution = rates; }
    in
    create model.alph new_spec

let add_gap_to_model compute_priors model = 
    match model.spec.use_gap with
    | `Missing -> add_gap_to_model compute_priors model
    | `Independent | `Coupled _ -> model

let remove_gamma_from_spec spec = 
    match spec.site_variation with
    | Some Constant | None -> spec
    | _ -> 
        Status.user_message Status.Warning dyno_gamma_warning;
        { spec with site_variation = None; }

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

(* Develop a model from a classification --created above *)
let spec_from_classification alph gap kind rates costfn (comp_map,pis) =
    let tuple_sum = 
        All_sets.FullTupleMap.fold (fun k v a -> a +. v) comp_map 0.0
    and ugap = match gap with
        | `Missing      -> false
        | `Independent  -> true
        | `Coupled _    -> true
    in
    let f_priors = 
        let sum = All_sets.IntegerMap.fold (fun k v x -> v +. x) pis 0.0
        and gap_size =
            try All_sets.IntegerMap.find (Alphabet.get_gap alph) pis
            with | Not_found -> 0.0
        in
        let l =
            (* subtract gaps from totals *)
            if not ugap then begin
                let sum = sum -. gap_size in
                List.fold_left
                    (fun acc (r,b) -> match b with
                        | b when b = Alphabet.get_gap alph -> acc
                        | b ->
                            let c = 
                                try (All_sets.IntegerMap.find b pis) /. sum
                                with | Not_found -> 0.0 in
                            c :: acc)
                    []
                    (Alphabet.to_list alph)
            (* gaps are charcters *)
            end else begin
                List.fold_left
                    (fun acc (r,b) -> 
                        let c = try (All_sets.IntegerMap.find b pis) /. sum
                                with | Not_found -> 0.0 in
                        c :: acc)
                    []
                    (Alphabet.to_list alph)
            end
        in
        (* TODO: what to do if we don't have any of a state? *)
        List.rev l
    and is_comp a b = match Alphabet.complement a alph with
        | Some x -> x = b
        | None   -> false
    in
    (* build the model *)
    let m = try match kind with
        | `JC69 -> JC69
        | `F81  -> F81
        | `K2P _ -> 
            (* YANG: 1.12                                       *)
            (*  alpha*t = -0.5 * log(1-2S-V) + 0.25 * log(1-2V) *)
            (* 2*beta*t = -0.5 * log(1-2V)                      *)   
            (*        k = alpha / beta                          *)  
            (* k = [log(1-2S-V)-0.5*log(1-2V)] / log(1-2V)      *)
            let s,v,a =
                let s,v,a = All_sets.FullTupleMap.fold
                    (fun k v (sc,vc,all) -> match k with
                            | k1,k2 when is_comp k1 k2 -> (sc+.v,vc,all+.v)
                            | k1,k2 when k1 != k2 -> (sc,vc+.v,all+.v)
                            | _ -> (sc,vc,all+.v) )
                    comp_map
                    (0.0,0.0,0.0)
                in
                s /. a, v /. a, a
            in
            let kappa = (log (1.0-.(2.0*.s)-.v))-.(log (1.0-.(2.0*.v)))
            and denom = log (1.0-.(2.0*.v)) in
            K2P (Some (kappa /. denom) )
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
            let lst = 
                List.fold_right
                    (fun (s1,alph1) acc1 ->
                        if (not ugap) && ((Alphabet.get_gap alph) = alph1) then
                            acc1
                        else
                            List.fold_right
                                (fun (s2,alph2) acc2 ->
                                    if alph2 <= alph1 then acc2
                                    else if (not ugap) && ((Alphabet.get_gap alph) = alph2) then
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
            GTR (Some (Array.of_list lst))
        | `F84 _
        | `Custom _
        | `HKY85 _ 
        | `TN93 _ 
        | `File _ -> failwith "I need you to specify a model first"
        with | Not_found -> failwith "Cannot find something for model"
    and calc_invar all comp_map =
        let same = 
            List.fold_left
                (fun acc (_,ac) ->
                    acc +. (All_sets.FullTupleMap.find (ac,ac) comp_map) )
                0.0
                (Alphabet.to_list alph)
        in
        same /. all
    in
    let v = match rates with
        | None -> None
        | Some (`Gamma (i,_)) -> Some (Gamma (i,1.0))
        | Some (`Theta (i,_)) -> Some (Theta (i,0.2,calc_invar tuple_sum comp_map))
    in
    {
        substitution = m;
        site_variation = v;
        base_priors = Estimated (Array.of_list f_priors);
        cost_fn = costfn;
        use_gap = gap;
        iterate_model = true;
        iterate_alpha = true;
    }

let convert_gapr m = function
    | `Coupled x   -> Some (Alphabet.get_gap m.alph,x)
    | `Independent -> None
    | `Missing     -> None

let update_jc69 old_model gap_r = 
    let subst_spec = { old_model.spec with use_gap = gap_r; }
    and gap_r = convert_gapr old_model gap_r in
    let subst_model = m_jc69 old_model.pi_0 1.0 old_model.alph_s gap_r in
    let u,d,ui = diagonalize true subst_model in
    { old_model with spec = subst_spec; s  = subst_model; u  = u; d  = d; ui = ui; }

and update_f81 old_model gap_r = 
    let subst_spec = { old_model.spec with use_gap = gap_r; }
    and gap_r = convert_gapr old_model gap_r in
    let subst_model = m_f81 old_model.pi_0 1.0 old_model.alph_s gap_r in
    let u,d,ui = diagonalize true subst_model in
    { old_model with spec = subst_spec; s  = subst_model; u  = u; d  = d; ui = ui; }

and update_k2p old_model new_value gap_r =
    let subst_spec = { old_model.spec with substitution = K2P (Some new_value);
                                           use_gap = gap_r; }
    and gap_r = convert_gapr old_model gap_r in
    let subst_model = m_k2p old_model.pi_0 new_value 1.0 old_model.alph_s gap_r in
    let u,d,ui = diagonalize true subst_model in
    { old_model with spec = subst_spec; s  = subst_model; u  = u; d  = d; ui = ui; }

and update_hky old_model new_value gap_r =
    let subst_spec = { old_model.spec with substitution = HKY85 (Some new_value);
                                           use_gap = gap_r; }
    and gap_r = convert_gapr old_model gap_r in
    let subst_model = m_hky85 old_model.pi_0 new_value old_model.alph_s gap_r in
    let u,d,ui = diagonalize false subst_model in
    { old_model with spec = subst_spec; s  = subst_model; u  = u; d  = d; ui = ui; }

and update_tn93 old_model ((x,y) as new_value) gap_r =
    let subst_spec = { old_model.spec with substitution = TN93 (Some new_value);
                                           use_gap = gap_r; }
    and gap_r = convert_gapr old_model gap_r in
    let subst_model = m_tn93 old_model.pi_0 x y 1.0 old_model.alph_s gap_r in
    let u,d,ui = diagonalize false subst_model in
    { old_model with spec = subst_spec; s  = subst_model; u  = u; d  = d; ui = ui; }

and update_f84 old_model new_value gap_r =
    let subst_spec = { old_model.spec with substitution = F84 (Some new_value);
                                           use_gap = gap_r; }
    and gap_r = convert_gapr old_model gap_r in
    let subst_model = m_f84 old_model.pi_0 new_value 1.0 old_model.alph_s gap_r in
    let u,d,ui = diagonalize false subst_model in
    { old_model with spec = subst_spec; s  = subst_model; u  = u; d  = d; ui = ui; }

and update_gtr old_model new_values gap_r =
    let normalized,gap_r = normalize {old_model.spec with use_gap = gap_r} new_values in
    let subst_spec = { old_model.spec with substitution = GTR (Some normalized);
                                           use_gap = gap_r; }
    and gap_r = convert_gapr old_model gap_r in
    let subst_model = m_gtr old_model.pi_0 normalized old_model.alph_s gap_r in
    let u,d,ui = diagonalize false subst_model in
    { old_model with spec = subst_spec; s  = subst_model; u  = u; d  = d; ui = ui; }

and update_custom old_model new_values =
    let subst_spec,assoc = match old_model.spec.substitution with
        | Custom (assoc,_,s) -> 
            { old_model.spec with substitution = Custom (assoc,new_values,s); },assoc
        | _ -> assert false
    in
    let subst_model = m_custom old_model.pi_0 assoc new_values old_model.alph_s in
    let u,d,ui = diagonalize false subst_model in
    { old_model with spec = subst_spec; s = subst_model; u = u; d = d; ui = ui; }

and update_alpha old_model new_value =
    let new_spec_var,new_array = match old_model.spec.site_variation with
        | Some (Gamma (i,_))   -> Some (Gamma (i,new_value)),
                                  Numerical.gamma_rates new_value new_value i
        | Some (Theta (i,_,p)) -> Some (Theta (i,new_value,p)),
                                  Numerical.gamma_rates new_value new_value i
        | _ -> failwith "Model doesn't specify site variation"
    in
    { old_model with rate = new_array;
                     spec = { old_model.spec with site_variation = new_spec_var;} }

and update_alpha_invar old_model new_alpha new_invar =
    let new_spec_var,new_array = match old_model.spec.site_variation with
        | Some (Gamma (i,_)) -> Some (Gamma (i,new_alpha)),
                                Numerical.gamma_rates new_alpha new_alpha i
        | Some (Theta (i,_,_)) -> Some (Theta (i,new_alpha,new_invar)),
                                  Numerical.gamma_rates new_alpha new_alpha i
        | _ -> failwith "Model doesn't specify site variation"
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
    let alph = model.alph in
    let priors = 
        let inner = 
            Array.mapi
                (fun i v ->
                    (PXML -[Xml.Characters.vector]
                        ([Xml.Alphabet.value] = [`Float v])
                        { `String (Alphabet.match_code i alph)} --))
                (match model.spec.base_priors with
                    | Equal ->
                        Array.make (model.alph_s) (1.0 /. (float_of_int model.alph_s))
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
        | None | Some Constant -> `Float 0.0
        | Some (Gamma (_,x)) 
        | Some (Theta (_,x,_))  -> `Float x
    and get_invar m = match m.spec.site_variation with
        | Some (Gamma _) | None | Some Constant -> `Float 0.0
        | Some (Theta (_,_,p))  -> `Float p
    and get_gap m = match m.spec.use_gap with
        | `Missing     -> `String "missing"
        | `Independent -> `String "independent"
        | `Coupled x   -> `String ("coupled:"^string_of_float x)
    and get_cats m = match m.spec.site_variation with
        | None | Some Constant -> `Int 1
        | Some (Gamma (c,_)) 
        | Some (Theta (c,_,_))  -> `Int c
    and parameters m = match m.spec.substitution with
        | JC69  | F81 -> []
        | K2P x | F84 x | HKY85 x ->
                let x = match x with 
                        | Some x -> x 
                        | None -> default_tstv
                in
                [(PXML -[Xml.Data.param 1]
                    ([Xml.Alphabet.value] = [`Float x])
                    { `String "" } --)]
        | TN93 x -> 
                let a,b = match x with
                        | Some (a,b) -> a,b
                        | None -> default_tstv,default_tstv
                in
                (PXML -[Xml.Data.param 1]
                    ([Xml.Alphabet.value] = [`Float a])
                    { `String "" } --) ::
                [(PXML -[Xml.Data.param 2]
                    ([Xml.Alphabet.value] = [`Float b])
                    { `String "" } --)]
        | GTR ray ->
                let ray = match ray with 
                    | None -> default_gtr m.alph_s
                    | Some x -> x
                in
                let r,_ =
                    Array.fold_right
                        (fun x (acc,i) ->
                            (PXML -[Xml.Data.param i]
                                ([Xml.Alphabet.value] = [`Float x])
                                { `String "" } --) :: acc,i+1)
                        ray
                        ([],1)
                in
                r
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


(* this assumes that the spec is consistent with the model itself, the returned
 * update function ensures this consistency *)
and get_update_function_for_model model =
  IFDEF USE_LIKELIHOOD THEN
    let split_array ray = 
        ray.((Array.length ray)-1),
        Array.init ((Array.length ray)-1) (fun i -> ray.(i)) in
    if model.spec.iterate_model then begin
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
    end else begin
        None
    end
  ELSE
    None
  END

and get_update_function_for_alpha model = 
  IFDEF USE_LIKELIHOOD THEN
    if model.spec.iterate_alpha then begin
        match model.spec.site_variation with
        | Some (Gamma _) -> Some update_alpha
        | Some (Theta _) -> Some update_alpha
        | None | Some Constant -> None
    end else begin
        None
    end
  ELSE
    None
  END

(* for model parameters *)
let get_current_parameters_for_model model =
  IFDEF USE_LIKELIHOOD THEN
    match model.spec.substitution, model.spec.use_gap with
        | JC69,`Coupled gapr | F81,`Coupled gapr ->
            Some (Array.make 1 gapr)
        | F84 x,`Coupled gapr | K2P x,`Coupled gapr | HKY85 x,`Coupled gapr ->
            let x = Array.make 2 (match x with | None -> default_tstv | Some x -> x) in
            x.(1) <- gapr;
            Some x
        | TN93 x,`Coupled gapr ->
            let x1,x2 = match x with 
                | Some (x,y) -> x,y
                | None -> default_tstv,default_tstv
            in
            let y = Array.make 3 x2 in
            y.(0) <- x1;
            y.(2) <- gapr;
            Some y
        | GTR (Some y),`Coupled gapr ->
            let y = Array.init ((Array.length y)+1)
                       (fun i -> if i = Array.length y then gapr else y.(i))
            in
            Some y
        | GTR None,`Coupled gapr ->
            let a = Array.make ((((model.alph_s-1)*model.alph_s)/2)+1) 1.0 in 
            a.((Array.length a)-1) <- gapr;
            Some a
        (* no gap below *)
        | JC69,_ | F81,_ | File _,_ -> None
        | F84 x,_ | K2P x,_ | HKY85 x,_ ->
            let x = match x with | None -> default_tstv | Some x -> x in
            Some (Array.make 1 x)
        | TN93 x,_  ->
            let x1,x2 = match x with 
                | Some (x,y) -> x,y
                | None -> default_tstv,default_tstv
            in
            let y = Array.make 2 x2 in
            y.(0) <- x1;
            Some y
        | GTR ((Some _) as x),_ -> x
        | GTR None,_ -> Some (default_gtr model.alph_s)
        | Custom (_,xs,_),_ -> Some xs
  ELSE
    None
  END

(* for alpha parameter *)
and get_current_parameters_for_alpha model =
  IFDEF USE_LIKELIHOOD THEN
    match model.spec.site_variation with
        | Some (Gamma (cat,alpha)) -> Some alpha
        | Some (Theta (cat,alpha,invar)) -> Some alpha
        | None | Some Constant -> None
  ELSE
    None
  END
