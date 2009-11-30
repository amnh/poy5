(* POY 4.0 Beta. A phylogenetic analysis program using Dynamic Homologies.    *)
(* Copyright (C) 2007  Andrés Varón, Le Sy Vinh, Illya Bomash, Ward Wheeler,  *)
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

let epsilon = 0.000001
let (=.) a b = abs_float (a-.b) < epsilon
let (>=.) a b = abs_float (a-.b) > ~-.epsilon
let (-->) b a = a b
let failwithf format = Printf.ksprintf failwith format

let likelihood_not_enabled =
    "Likelihood not enabled: download different binary or contact mailing list" 

let debug = false
let debug_printf msg format = 
    Printf.ksprintf (fun x -> if debug then print_string x; flush stdout) msg format

let ba_of_array1 x = Bigarray.Array1.of_array Bigarray.float64 Bigarray.c_layout x
and ba_of_array2 x = Bigarray.Array2.of_array Bigarray.float64 Bigarray.c_layout x

let create_ba1 x     = Bigarray.Array1.create Bigarray.float64 Bigarray.c_layout x
and create_ba2 x y   = Bigarray.Array2.create Bigarray.float64 Bigarray.c_layout x y

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
and pp_farray xs =
    (Array.fold_left (fun acc x -> acc^"| "^(string_of_float x)^" ") "[" xs)^" |]"

(* type to help the parsing of specification data *)
type string_spec = string * (string * string * string * string)
                          * float list * ( string * float ) list * string * string option

let empty_str_spec = ("",("","","",""),[],[],"",None)

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
(* --- ------------------------------ --- *)

(* list of set bits, or packed integer of set bits *)
type chars = [ `List of int list | `Packed of int ]

let list_of_packed d =
    let rec loop_ c i d = match d land 1 with
        | 0 when d = 0 -> c
        | 0  -> loop_ c (i+1) (d lsr 1)
        | 1  -> loop_ (i::c) (i+1) (d lsr 1)
        | _  -> failwith "MlModel.classify_seq_pairs.build_lst"
    in loop_ [] 0 d

and packed_of_list lst =
    List.fold_left (fun acc x -> (1 lsl x) + acc) 0 lst

type site_var = 
    | Gamma of int * float
    | Theta of int * float * float
    (* | Given of (float * float) array *)
    | Constant 

type subst_model =
    | JC69
    | F81
    | K2P   of float option
    | F84   of float option
    | HKY85 of float option
    | TN93  of (float * float) option
    | GTR   of float array option
    | File  of float array array 

type priors = 
    | Estimated of float array
    | Given     of float array

type spec = {
    substitution : subst_model;
    site_variation : site_var option;
    base_priors : priors;
    use_gap : bool;
}

type model = {
    spec  : spec;
    alph  : int;
    pi_0  : (float, Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array1.t;
    invar : float option;
    rate  : (float, Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array1.t;
    prob  : (float, Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array1.t;
    s     : (float, Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array2.t;
    u     : (float, Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array2.t;
    d     : (float, Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array2.t;
    ui    : (float, Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array2.t option; 
}


IFDEF USE_LIKELIHOOD THEN

(* a gentler compare that excludes the parameters of the model itself *)
let compare a b = 
    let m_compare = match a.spec.substitution,b.spec.substitution with
        | JC69 , JC69 | F81 , F81 | K2P _, K2P _ | F84 _, F84 _ 
        | HKY85 _, HKY85 _ | TN93 _,  TN93 _ | GTR _, GTR _ 
        | File _, File _ -> 0
        | _,_ -> ~-1
    and v_compare = match a.spec.site_variation,b.spec.site_variation with
        | Some (Gamma _), Some (Gamma _) | Some (Theta _), Some (Theta _) -> 0
        | None, None | Some Constant, Some Constant 
        | Some Constant, None | None, Some Constant -> 0
        | _,_ -> ~-1
    in
    (* just knowing that they are different is enough *)
    m_compare + v_compare
    
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
 * usually done on the C side exclusively *)
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

(* calculates the gamma rates for specific alpha, beta and #classes *)
external gamma_rates: float -> float -> int ->
    (float,Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array1.t = 
        "gamma_CAML_rates"

(* ------------------------------------------------ *)
(* MODEL CALCULATION FUNCTIONS                      *)

let read_file file =
    (* explode a string around a character;filtering empty results *)
    let explode str ch =
        let rec expl s i l =
            if String.contains_from s i ch then
                let spac = String.index_from s i ch in
                let word = String.sub s i (spac-i) in
                expl s (spac+1) (word::l)
            else
                let final = String.sub s i ((String.length s)-i) in
                final::l
        in
        List.filter (fun x-> if x = "" then false else true)
                    (List.rev (expl str 0 []))
    in
    (* read a channel line by line and applying f into a list *)
    let rec read_loop f chan =
        try let line = FileStream.Pervasives.input_line chan in
            (List.map (float_of_string) (f line ' ') ) :: read_loop f chan
        with e -> []
    in
    let f = FileStream.Pervasives.open_in file in
    let mat = read_loop (explode) f in
    let _ = FileStream.Pervasives.close_in f in 
    mat

(* divide a matrix by the mean rate so it will equal 1 *)
let m_meanrate srm pi_ =
    let mr = ref 0.0 and a_size = Bigarray.Array2.dim1 srm in
    for i = 0 to (a_size-1) do
        mr := !mr +. (~-.(srm.{i,i}) *. pi_.(i));
    done;
    for i = 0 to (a_size-1) do
        for j = 0 to (a_size-1) do
            srm.{i,j} <- srm.{i,j} /. !mr;
        done;
    done

(* val jc69 :: ANY ALPHABET size *)
let m_jc69 pi_ mu a_size =
    let srm = create_ba2 a_size a_size in
    Bigarray.Array2.fill srm mu;
    let diag = -. mu *. float (a_size-1) in
    for i = 0 to (a_size-1) do
        srm.{i,i} <- diag
    done;
    m_meanrate srm pi_;
    srm

(* val k2p :: only 4 or 5 characters *)
let m_k2p pi_ alpha beta a_size =
    let srm = create_ba2 a_size a_size in
    Bigarray.Array2.fill srm beta;
    (* modify transition elements to alpha *)
    srm.{1, 3} <- alpha; srm.{3, 1} <- alpha;
    srm.{2, 0} <- alpha; srm.{0, 2} <- alpha;
    (* set up the diagonal elements *)
    let diag = if a_size = 4 then
                    -. alpha -. beta -. beta
               else -. alpha -. beta *. 3.0  in
    for i = 0 to (a_size-1) do
        srm.{i,i} <- diag
    done;
    m_meanrate srm pi_;
    srm

(* val tn93 :: only 4 or 5 characters *)
let m_tn93 pi_ alpha beta gamma a_size =
    let srm = create_ba2 a_size a_size in
    Bigarray.Array2.fill srm gamma;
    srm.{0,2} <- alpha; srm.{1,3} <- beta; (* ACGT -- R=AG -- Y=CT *)
    srm.{2,0} <- alpha; srm.{3,1} <- beta; (* 0123 -- R=02 -- Y=13 *)
    for i = 0 to (a_size-1) do
        for j = 0 to (a_size-1) do
            srm.{i,j} <- srm.{i,j} *. pi_.(j);
        done;
    done;
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

let m_tn93_ratio pi_ kappa1 kappa2 a_size = 
    let beta = (pi_.(0) *. pi_.(2) *. kappa1) +. (pi_.(1) *. pi_.(3) *. kappa2) +.
                ((pi_.(0) +. pi_.(2)) *. (pi_.(1)+.pi_.(3)) ) in
    let beta = 1.0 /. (2.0 *. beta) in
    let alpha1 = kappa1 *. beta and alpha2 = kappa2 *. beta in
    m_tn93 pi_ alpha1 alpha2 beta a_size

(* val f81 :: ANY ALPHABET size *)
let m_f81 pi_ lambda a_size =
    let srm = create_ba2 a_size a_size in
    for i = 0 to (a_size-1) do
        for j = 0 to (a_size-1) do
            srm.{i,j} <- 
                if i = j then
                     -. lambda *. (1.0 -. pi_.(i))
                else 
                    pi_.(j)*.lambda;
        done;
    done;
    m_meanrate srm pi_;
    srm

(* val hky85 :: only 4 or 5 characters *)
let m_hky85 pi_ kappa a_size =
    let beta = (pi_.(0) *. pi_.(2) *. kappa) +. (pi_.(1) *. pi_.(3) *. kappa) +.
                ((pi_.(0) +. pi_.(2)) *. (pi_.(1)+.pi_.(3)) ) in
    let beta = 1.0 /. (2.0 *. beta) in
    let alpha = kappa *. beta in
    m_tn93 pi_ alpha alpha beta a_size

(* val f84 :: only 4 or 5 characters *)
let m_f84 pi_ gamma kappa a_size =
    let y = pi_.(1) +. pi_.(3) in (* Y = C + T *)
    let r = pi_.(0) +. pi_.(2) in (* R = A + G *)
    let alpha = (1.0+.kappa/.r) *. gamma in
    let beta = (1.0+.kappa/.y) *. gamma in
    m_tn93 pi_ alpha beta gamma a_size 

(* val gtr :: ANY ALPHABET size
   pi_ == n; co_ == ((n-1)*n)/2
   form of lower packed storage mode, excluding diagonal, *)
let m_gtr pi_ co_ a_size =
    if ((a_size*(a_size-1))/2) <> Array.length co_ then begin
        failwithf "Length of GTR insufficient against alphabet: %d != %d\n" 
                    (Array.length co_) a_size;
    end;
    (* last element of GTR = 1.0 *)
    let last = co_.( (Array.length co_) - 1) in
    let co_ = Array.map (fun x -> x /. last) co_ in
    (* create matrix *)
    let n = ref 0 in (* array index *)
    let srm = create_ba2 a_size a_size in
    for i = 0 to (a_size-1) do
        for j = (i+1) to (a_size-1) do
            srm.{i,j} <- co_.(!n) *. pi_.(j);
            srm.{j,i} <- co_.(!n) *. pi_.(i);
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

(* val m_file :: any alphabet size -- recomputes diagonal and divides by meanrate *)
let m_file pi_ f_rr a_size =
    assert(a_size = Array.length f_rr);
    let srm = create_ba2 a_size a_size in
    for r = 0 to (a_size-1) do
        assert(a_size = Array.length f_rr.(r));
        let diag = ref 0.0 in
        for c = 0 to (a_size-1) do
            if (c <> r) then 
                (diag := !diag +. f_rr.(r).(c); srm.{r,c} <- f_rr.(r).(c);)
        done;
        srm.{r,r} <- ~-. !diag;
    done;
    m_meanrate srm pi_; (* TODO: should this be done if given? *)
    srm


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

END

(* ------------------------------------------------ *)
(* CONVERSTION/MODEL CREATION FUNCTIONS             *)

(* convert a string spec to a specification, used in Parser for nexus *)
let convert_string_spec ((name,(var,site,alpha,invar),param,priors,gap,file):string_spec) =
  IFDEF USE_LIKELIHOOD THEN
    let submatrix = match String.uppercase name with
        | "JC69" -> (match param with
            | [] -> JC69
            | _ -> failwith "Parameters don't match model")
        | "F81" -> (match param with
            | [] -> F81
            | _ -> failwith "Parameters don't match model")
        | "K80" | "K2P" -> (match param with
            | ratio::[] -> K2P (Some ratio)
            | []        -> K2P None
            | _ -> failwith "Parameters don't match model")
        | "F84" -> (match param with
            | ratio::[] -> F84 (Some ratio)
            | []        -> F84 None
            | _ -> failwith "Parameters don't match model")
        | "HKY" | "HKY85" -> (match param with
            | ratio::[] -> HKY85 (Some ratio)
            | []        -> HKY85 None
            | _ -> failwith "Parameters don't match model")
        | "TN93" -> (match param with
            | ts::tv::[] -> TN93 (Some (ts,tv))
            | []         -> TN93 None
            | _ -> failwith "Parameters don't match model")
        | "GTR" -> (match param with
            | [] -> GTR None
            | ls -> GTR (Some (Array.of_list ls)) )
        | "GIVEN"-> (match file with
            | Some name ->
                (`Local name)
                    --> read_file
                    --> List.map (Array.of_list)
                    --> Array.of_list
                    --> (fun x -> File x)
            | _ -> failwith "File not specified; Erasing hard drive. done.")
        (* ERROR *)
        | "" -> failwith "No Model Specified"
        | _  -> failwith "Incorrect Model"
    and variation = match String.uppercase var with
        | "GAMMA" ->
            Gamma (int_of_string site, float_of_string alpha)
        | "THETA" -> 
            Theta (int_of_string site,
                   float_of_string alpha,
                   float_of_string invar)
        | "NONE" | "CONSTANT" | "" | _ -> Constant
    and gap = match String.uppercase gap with
        | "TRUE" -> true
        | "FALSE" | _ -> false
    and priors = List.map snd priors
    in
    {   substitution = submatrix;
        site_variation = Some variation;
        base_priors = Given (Array.of_list priors);
        use_gap = gap;
    }
  ELSE
    failwith likelihood_not_enabled
  END

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
let create alph lk_spec = 
  IFDEF USE_LIKELIHOOD THEN
    let alph = Alphabet.to_sequential alph in
    let (a_size,a_gap) = match lk_spec.use_gap with
        | true -> Alphabet.size alph, ~-1
        | false -> (Alphabet.size alph) - 1, Alphabet.get_gap alph
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
                gamma_rates y y x,p,None
            | Theta (w,x,z) -> (* GAMMA_SITES,ALPHA,PERCENT_INVAR *)
                if w < 1 then
                    failwith "Number of rate categories must be >= 1 (if invar w/out gamma)"
                else if w = 1 then begin
                    ba_of_array1 [| 1.0 |], ba_of_array1 [| 1.0 |], Some z
                end else begin
                    let p = create_ba1 w in
                    Bigarray.Array1.fill p (1.0 /. (float_of_int w));
                    let r = gamma_rates x x w in
                    r,p,Some z
                end
        end
    in
    assert( verify_rates probabilities variation );
    (* extract the prior probability *)
    let priors = match lk_spec.base_priors with
        | Estimated p
        | Given p -> assert(a_size = Array.length p); p
    in
    (*  get the substitution rate matrix and set sym variable and to_formatter vars *)
    let sym, sub_mat, subst_model = match lk_spec.substitution with
        | JC69  -> true,  m_jc69 priors 1.0 a_size, JC69
        | F81   -> false, m_f81 priors 1.0 a_size, F81
        | K2P t -> 
            let t = match t with | Some x -> x | None -> default_tstv in
            true,  m_k2p priors t 1.0 a_size, lk_spec.substitution
        | F84 t -> 
            let t = match t with | Some x -> x | None -> default_tstv in
            false, m_f84 priors t 1.0 a_size, lk_spec.substitution  
        | HKY85 t ->
            let t = match t with | Some x -> x | None -> default_tstv in
            false, m_hky85 priors t a_size, lk_spec.substitution
        (* more complex models *)
        | TN93 tstv ->
            let ts,tv = match tstv with 
                | Some (x,y) -> x,y 
                | None       -> default_tstv,default_tstv
            in
            false, m_tn93 priors ts tv 1.0 a_size, lk_spec.substitution
        | GTR c ->
            let c = match c with 
                | Some xs ->
                    (* normalize so last element is = 1.0 *)
                    let last = xs.( (Array.length xs) - 1) in
                    Array.map (fun x -> x /. last) xs
                | None -> default_gtr a_size
            in
            false, m_gtr priors c a_size, (GTR (Some c))
        | File m-> false, m_file priors m a_size, lk_spec.substitution
    in
    let (u_,d_,ui_) = diagonalize sym sub_mat in
    {
      rate = variation;
      prob = probabilities;
      spec = {lk_spec with substitution = subst_model; };
     invar = invar;
      pi_0 = ba_of_array1 priors;
      alph = a_size;
         s = sub_mat;
         u = u_;
         d = d_;
        ui = ui_;
    }
  ELSE
    failwith likelihood_not_enabled
  END

IFDEF USE_LIKELIHOOD THEN

(* ------------------------------------------------ *)
(* MODEL ESTIMATION FUNCTIONS                       *)
let chars2str = function
    | `Packed s -> string_of_int s
    | `List s   -> "[ "^(String.concat " | " (List.map string_of_int s))^" ]"

(** estimate the model based on two sequences with attached weights *)
let classify_seq_pairs leaf1 leaf2 seq1 seq2 acc =
    let chars_to_list = function
        | `Packed s -> list_of_packed s
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
let spec_from_classification alph gap (kind:Methods.ml_substitution) rates (comp_map,pis) =
    let alph_size = 
        float_of_int (if gap then Alphabet.size alph else (Alphabet.size alph) - 1)
    and tuple_sum = 
        All_sets.FullTupleMap.fold (fun k v a -> a +. v) comp_map 0.0 in
    let f_priors = 
        let sum = All_sets.IntegerMap.fold (fun k v x -> v +. x) pis 0.0
        and gap_size =
            let s = try All_sets.IntegerMap.find (Alphabet.get_gap alph) pis
                    with | Not_found -> 0.0 in
            s /. alph_size
        in
        let l = List.fold_left
            (fun acc (r,b) ->
                (* if gaps aren't characters then incorporate them *)
                if not gap then begin
                    if b = Alphabet.get_gap alph then acc
                    else begin
                        let c = try (gap_size +. (All_sets.IntegerMap.find b pis)) /. sum
                                with | Not_found -> gap_size /. sum in
                        c :: acc
                    end
                (* gaps are characters, seperate them -- minimum is 1/sum.
                 * essentially saying one character must exist, else a different
                 * alphabet should have been chosen *)
                end else begin
                    let c = try (All_sets.IntegerMap.find b pis) /. sum
                            with | Not_found -> 0.0 in
                    c :: acc
                end)
            []
            (Alphabet.to_list alph)
        in 
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
                        if (not gap) && ((Alphabet.get_gap alph) = alph1) then
                            acc1
                        else
                            List.fold_right
                                (fun (s2,alph2) acc2 ->
                                    if alph2 <= alph1 then acc2
                                    else if (not gap) && ((Alphabet.get_gap alph) = alph2) then
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
        | Some (`Theta (i,_,_)) -> Some (Theta (i,0.2,calc_invar tuple_sum comp_map))
    in
    {
        substitution = m;
        site_variation = v;
        base_priors = Estimated (Array.of_list f_priors);
        use_gap = gap;
    }

let update_k2p old_model new_value =
    let subst_spec = { old_model.spec with substitution = K2P (Some new_value) }
    and priors = match old_model.spec.base_priors with | Estimated x | Given x -> x in
    let subst_model = m_k2p priors new_value 1.0 old_model.alph in
    let u,d,ui = diagonalize true subst_model in
    { old_model with spec = subst_spec; s  = subst_model; u  = u; d  = d; ui = ui; }

and update_hky old_model new_value =
    let subst_spec = { old_model.spec with substitution = HKY85 (Some new_value) }
    and priors = match old_model.spec.base_priors with | Estimated x | Given x -> x in
    let subst_model = m_hky85 priors new_value old_model.alph in
    let u,d,ui = diagonalize false subst_model in
    { old_model with spec = subst_spec; s  = subst_model; u  = u; d  = d; ui = ui; }

and update_tn93 old_model ((x,y) as new_value) =
    let subst_spec = { old_model.spec with substitution = TN93 (Some new_value) }
    and priors = match old_model.spec.base_priors with | Estimated x | Given x -> x in
    let subst_model = m_tn93 priors x y 1.0 old_model.alph in
    let u,d,ui = diagonalize false subst_model in
    { old_model with spec = subst_spec; s  = subst_model; u  = u; d  = d; ui = ui; }

and update_f84 old_model new_value =  
    let subst_spec = { old_model.spec with substitution = F84 (Some new_value) }
    and priors = match old_model.spec.base_priors with | Estimated x | Given x -> x in
    let subst_model = m_f84 priors new_value 1.0 old_model.alph in
    let u,d,ui = diagonalize false subst_model in
    { old_model with spec = subst_spec; s  = subst_model; u  = u; d  = d; ui = ui; }

and update_gtr old_model new_values =  
    let subst_spec = { old_model.spec with substitution = GTR (Some new_values) }
    and priors = match old_model.spec.base_priors with | Estimated x | Given x -> x in
    let subst_model = m_gtr priors new_values old_model.alph in
    let u,d,ui = diagonalize false subst_model in
    { old_model with spec = subst_spec; s  = subst_model; u  = u; d  = d; ui = ui; }

and update_alpha old_model new_value = 
    let new_spec_var,new_array = match old_model.spec.site_variation with
        | Some (Gamma (i,_))   -> Some (Gamma (i,new_value)),
                                  gamma_rates new_value new_value i
        | Some (Theta (i,_,p)) -> Some (Theta (i,new_value,p)),
                                  gamma_rates new_value new_value i
        | _ -> failwith "Model doesn't specify site variation"
    in
    { old_model with rate = new_array;
                     spec = { old_model.spec with site_variation = new_spec_var;} }

and update_alpha_invar old_model new_alpha new_invar =
    let new_spec_var,new_array = match old_model.spec.site_variation with
        | Some (Gamma (i,_)) -> Some (Gamma (i,new_alpha)),
                                gamma_rates new_alpha new_alpha i
        | Some (Theta (i,_,_)) -> Some (Theta (i,new_alpha,new_invar)),
                                  gamma_rates new_alpha new_alpha i
        | _ -> failwith "Model doesn't specify site variation"
    in
    { old_model with rate = new_array;
                     spec = { old_model.spec with site_variation = new_spec_var; };
                     invar = Some new_invar; }

let update_all fun_model old_model new_values =
    let n = Array.length new_values in
    let alpha = new_values.(n - 1) in
    let model = Array.init (n-1) (fun x -> new_values.(x)) in
    update_alpha (fun_model old_model model) alpha

let update_alli fun_model old_model new_values = 
    let n = Array.length new_values in
    let alpha = new_values.(n - 2) 
    and invar = new_values.(n - 1) in
    let model = Array.init (n - 2) (fun x -> new_values.(x)) in
    update_alpha_invar (fun_model old_model model) alpha invar

END



(** GENERAL BRENTS METHOD **)
exception Colinear
(* function and braketed area and error *)
let brents_method ?(iter_max=1000) ?(epsilon=epsilon) ((orig_val,orig_ll) as orig) f =
    (* equality and greater-than functions with epsilon error *)
    let (=.) a b = abs_float (a-.b) < epsilon
    and (>=.) a b = abs_float (a-.b) > ~-.epsilon
    (* get and decrement, because it makes more sense *)
    and decr x = decr x; !x 
    (* a point between a and b such that (a-x)/(b-x) = phi *)
    and golden_middle a b = 
        let a,b = if a < b then a,b else b,a in
        a +. ((b -. a) *. 2.0 /. (1.0 +. sqrt 5.0))
    (* a point outside of a and b such that (a-b)/(x-b) = phi *)
    and golden_exterior a b =
        a +. ((b -. a) *. 2.0 /. ((sqrt 5.0) -. 1.0))
    (* the minimum of a parabola based on three points *)
    and abscissa (a,fa) (b,fb) (c,fc) =
        let numer = ((b-.a)*.(b-.a)*.(fb-.fc)) -. ((b-.c)*.(b-.c)*.(fb-.fa))
        and denom = ((b-.a)*.(fb-.fc)) -. ((b-.c)*.(fb-.fa)) in
        if denom = 0.0 then raise Colinear
        else b -. (0.5 *. (numer /. denom))
    (* set the total number of iterations to take place *)
    and iter = ref iter_max in
    (* order the three types *)
    let order_triples ((a1,_) as a) ((b1,_) as b) ((c1,_) as c) =
        if a1 < b1 then
           if b1 < c1 then (a,b,c)
           else if a1 < c1 then (a,c,b)
           else (c,a,b)
        else begin
           if a1 < c1 then (b,a,c)
           else if c1 < b1 then (c,b,a)
           else (b,a,c)
        end
    (* choose best likelihood of two values *)
    and best_of ((_,x) as x') ((_,y) as y') = if x <= y then x' else y' in
    (* parabolic interpolation *)
    let rec parabolic_interp ((best_t,best_l) as best) ((av,fa) as a) ((bv,fb) as b) ((cv,fc) as c) = 
        assert( av >=. 0.0 && bv >=. 0.0 && cv >=. 0.0 );
        assert( fa >=. 0.0 && fb >=. 0.0 && fc >=. 0.0 );
        if ((abs_float (fb -. fa)) <= epsilon) or 
           ((abs_float (fc -. fa)) <= epsilon) or (decr iter) = 0 then best
        else
            try
                let xv = abscissa a b c in let fx = f xv in let x = xv,fx in
                if fx =. fb || xv =. bv || xv <= epsilon then best
                else begin
                    let best = best_of best x in
                    if av < xv && xv < bv      then parabolic_interp best a x b
                    else if xv < av            then parabolic_interp best x a b
                    else if bv < xv && xv < cv then parabolic_interp best b x c
                    else if cv < xv            then parabolic_interp best b c x
                    else failwith "shouldn't happen"
                end
            with | Colinear -> brent_decision best a b c
    (* golden section search, when function is crappy --this isn't needed when
     * there is only one minima, as in branch lengths *)
    and golden_ratio ((best_t,best_l) as best) ((av,fa) as a) ((bv,fb) as b) ((cv,fc) as c)  =
        assert( av >=. 0.0 && bv >=. 0.0 && cv >=. 0.0 );
        let best,a,nb,c = 
            if (abs_float (cv-.bv)) > (abs_float (bv-.av)) then
                let other = golden_middle bv cv in let other = (other,f other) in
                let best = best_of best other in
                best,b,other,c
            else 
                let other = golden_middle av bv in let other = (other,f other) in
                let best = best_of best other in
                best,a,other,b
        in
        if fb =. (snd nb) or (decr iter) = 0 then best
        else brent_decision best a nb c
    (* brent exponential search, when points are colinear or monotonic
     * does not return a result since we are widening the search area. *)
    and brent_exp ((best_t,best_l) as best) a b c : float * float =
        assert (snd c < snd a); (* since we estimate "past" c, always *)
        let n = match golden_exterior (fst a) (fst c) with
            | x when x <= 0.0 -> epsilon
            | x -> x 
        in
        let other = n,f n in
        brent_decision (best_of best other) b c other
    (* What to do, what to do? Well, if the three points are monotonic, then
     * call brent_exp until something better comes up, if there is a minimum
     * between we can do a parabolic interpolation, else we use a shitty golden
     * bisect search method each iteration to find a better spot. *) 
    and brent_decision ((bv,fb) as b) l m u : float * float =
        let ((lv,fl) as l), ((mv,fm) as m), ((uv,fu) as u) = order_triples l m u in
        debug_printf "Iteration %d: B:(%f,%f) L:(%f,%f) M:(%f,%f) U:(%f,%f)\n%!"
                            (iter_max - !iter) bv fb lv fl mv fm uv fu;
        if lv =. uv then                  (* ending condition(s) *)
            b
        else if fl <= fm && fm <= fu then (* monotonically increasing *)
            brent_exp b u m l
        else if fl >= fm && fm >= fu then (* monotonically decreasing *)
            brent_exp b l m u
        else if fm <= fl && fm <= fu then (* minimum between *)
            parabolic_interp b l m u
        else                              (* shouldn't really happen *)
            golden_ratio b l m u
    in
    (* set up variables.. order arguments,find golden middle and evaluate *)
    let lower = 0.20 *. orig_val and upper = 2.0 *. orig_val in
    let middle = golden_middle lower upper in
    let l = lower,f lower and m = middle,f middle and u = upper,f upper in
    let best = best_of orig (best_of l (best_of m u)) in
    let ((new_val,new_ll) as new_) = brent_decision best l m u in
    debug_printf "Iterated: %d\tImprovement: %f\tVariable: %f -> %f\n%!"
                  (iter_max - !iter) (orig_ll -. new_ll) orig_val new_val;
    new_

(* find the derivative of a single variable function *)
let derivative_at_x ?(epsilon=3.0e-8) f x fx =
    let _,f_new = f (x +. epsilon) in
    (f_new -. fx) /. epsilon
(* find the magnitude of a vector x_array *)
let magnitude x_array = sqrt (Array.fold_left (fun acc x -> acc +. (x *. x)) 0.00 x_array)
(* find the gradient of a multi-variant function at a point x_array *)
let gradient_at_x ?(epsilon=3.0e-8) f_ x_array f_array : float array = 
    let i_replace i x v = let y = Array.copy x in Array.set y i v; y in
    Array.mapi 
        (fun i i_val ->
            derivative_at_x ~epsilon 
                            (fun x -> 
                                let newvec = i_replace i x_array x in
                                let newlk = f_ newvec in
                                debug_printf "\t[%s] -- %f\n" (pp_farray newvec) (snd newlk);
                                newlk)
                            i_val
                            f_array)
        x_array
(* dot product of two arrays *)
let dot_product x_array y_array = 
    let n = Array.length x_array and r = ref 0.0 in
    assert (n = Array.length y_array);
    for i = 0 to n-1 do
        r := !r +. (x_array.(i) *. y_array.(i));
    done;
    !r
(* map a matrix with a function *)
let matrix_map f mat = 
    let n1 = Array.length mat and n2 = Array.length mat.(0) in
    let output_matrix = Array.create_matrix n1 n2 0.0 in
    for i = 0 to n1 - 1 do for j = 0 to n2 -1 do 
        output_matrix.(i).(j) <- f i j mat.(i).(j);
    done; done;
    output_matrix
(* line search along a specified direction *)
(* Numerical Recipes in C : 9.7            *)
let line_search ?(epsilon=1.0e-7) ?(alf=1.0e-4) f point fpoint gradient maxstep direction =
    let (=.) a b = epsilon > (abs_float (a -. b)) and get_score x = snd x in
    (* set up some globals for the function to avoid tons of arguments *)
    let n = Array.length point and origfpoint = get_score fpoint in
    (* scale direction so, |pstep| <= maxstep *)
    let setup_function point direction gradient = 
        let direction = (* ||dir|| <= maxstep *)
            let magstep = magnitude direction in
            if magstep > maxstep then Array.map (fun x -> x *. maxstep /. magstep) direction 
            else direction in
        let slope = dot_product direction gradient
        and minstep = 
            let r = ref 0.0 in
            Array.iteri 
                (fun i x -> 
                    let tmp = (abs_float direction.(i)) /. (max (abs_float x) 1.0) in
                    if tmp > !r then r := tmp)
                point;
            epsilon /. (!r)
        and step = 1.0 in
        direction, slope, minstep, step
    (* find a new point and *)
    and next_step step prevstep slope newfpoint prevfpoint = 
        let newstep = 
            if step =. 1.0 then 
                ~-. slope /. (2.0 *. (newfpoint -. origfpoint -.  slope))
            else begin
                let tstep = 
                    let rhs1 = newfpoint -. origfpoint -. (step *. slope)
                    and rhs2 = prevfpoint -. origfpoint -. (prevstep *. slope) in
                    let rhs1divstepstep = rhs1 /. (step *. step)
                    and rhs2divpsteppstep = rhs2 /. (prevstep *. prevstep) in
                    let a =  rhs1divstepstep -. rhs2divpsteppstep
                    and b = (step *. rhs2divpsteppstep) -. (prevstep *. rhs1divstepstep) in
                    let a = a /. (step -. prevstep) and b = b /. (step -. prevstep) in

                    if a =. 0.0 then ~-. slope /. (2.0 *. b)
                    else begin
                        let disc = (b *. b) -. (3.0 *. a *. slope) in
                        if disc < 0.0 then 0.5 *. step
                        else if b <= 0.0 then (~-. b +. (sqrt disc)) /. (3.0 *. a)
                        else ~-. slope /. (b +. (sqrt disc))
                    end
                in
                min tstep (0.5 *. step)
            end
        in
        max newstep (0.1 *. step)
    in
    (* main algorithm -- first instance sets up some variables *)
    let rec main_ prevfpoint slope direction step prevstep minstep = 
        if step < minstep then
            (point,fpoint,true) (* -~verify convergence~- *)
        else begin
            let newpoint = Array.init n (fun i -> point.(i) +. (step *. direction.(i))) in
            let newfpoint = f newpoint in
            if (get_score newfpoint) <= origfpoint then begin
                debug_printf "\t\t%f--Accepting %f\n" prevstep (get_score newfpoint);
                (newpoint,newfpoint,false)
            end else begin
                let newstep = next_step step prevstep slope (get_score newfpoint) prevfpoint in
                debug_printf "\t\t%f--Rejecting %f\n" step (get_score newfpoint);
                main_ (get_score newfpoint) slope direction newstep step minstep
            end
        end
    in
    (* initialize and run... *)
    let direction, slope, minstep, step = setup_function point direction gradient in
    debug_printf "\tInitial LineSearch: %f, slope: %f\n" origfpoint slope;
    main_ origfpoint slope direction step step minstep 

(** BFGS Algorithm                   **)
(** Numerical Recipes in C : 10.7    **)
let bfgs_method ?(max_iter=200) ?(epsilon=3.0e-8) ?(mx_step=10.0) ?(g_tol=1.0e-5) f p fp =
   let n = Array.length p and get_score x = snd x in
    (* test convergence of a point --that it equals the direction, essentially *)
    let converged_l direction test_array =
        let test = ref 0.0 in
        Array.iteri
            (fun i x->
                let temp = max (abs_float x) 1.0 in
                test := max ((abs_float direction.(i)) /. temp) !test)
            test_array;
        (!test < (epsilon *. 4.0))
    (* Test tolerance for zeroing the gradient *)
    and converged_g fp gradient test_array =
        let test = ref 0.0 in
        Array.iteri
            (fun i x ->
                let temp = (max (abs_float x) 1.0) /. (max fp 1.0) in
                let temp = (abs_float gradient.(i)) *. temp in
                if temp > !test then test := temp)
            test_array;
        (!test < g_tol)
    (* Setup initial hessian (identity), initial gradiant vector, maximum step and direction *)
    and setup_function f_array x_array fx_array =
        let hessian =
            let h = Array.make_matrix n n 0.0 in
            for i = 0 to n-1 do h.(i).(i)  <- 1.0 done;
            h
        and x_grad = gradient_at_x f_array x_array fx_array in
        let dir = Array.init n (fun i -> ~-. (x_grad.(i)) )
        and mxstep = mx_step *. (max (magnitude x_array) (float_of_int n)) in
        hessian, x_grad, mxstep, dir in
    (* Do a line search step --return new p, new fp, new dir, if converged *)
    let line_searcher f p fp gradient step dir =
        let np,nfp,_ = line_search f p fp gradient step dir in
        let dir = Array.init n (fun i -> np.(i) -. p.(i) ) in
        np, nfp, dir, (converged_l dir np)
    (* update gradient --ret new gradient, difference of gradients,
     * difference of gradient times hessian matrix, if converged *)
    and gradient_update hessian ograd f p fp = 
        let ngrad = gradient_at_x f p fp in
        let dgrad = Array.init n (fun i -> ngrad.(i) -. ograd.(i)) in
        let hgrad = 
            Array.init n
                (fun i -> 
                    let res = ref 0.0 in
                    for j = 0 to n-1 do
                        res := !res +. (hessian.(i).(j) *. dgrad.(j));
                    done; !res)
        in
        ngrad, dgrad, hgrad, (converged_g fp ngrad p)
    (* Update the hessian matrix --skips the update if fac not sufficiently positive *)
    and bfgs_update_matrix dgrad hgrad direc hessian =
        let fac = dot_product dgrad direc
        and fae = dot_product dgrad hgrad
        and sumdgr = Array.fold_left (fun a x -> (x *. x) +. a) 0.0 dgrad
        and sumdir = Array.fold_left (fun a x -> (x *. x) +. a) 0.0 direc in
        if (fac *. fac) <= (epsilon *. sumdgr *. sumdir) then
            hessian
        else begin
            let fac = 1.0 /. fac and fad = 1.0 /. fae in
            let dgrad = 
                Array.init n
                    (fun i -> (fac *. direc.(i)) -. (fad *. hgrad.(i)))
            in
            matrix_map
                (fun i j x -> x +. (fac *. direc.(i) *. direc.(j))
                                -. (fad *. hgrad.(i) *. hgrad.(j))
                                +. (fae *. dgrad.(i) *. dgrad.(j)) )
                hessian
        end
    (* calculate the new direction for the line search by the hessian matrix *)
    and calculate_direction (hessian:float array array) (gradient:float array) : float array =
        Array.init n
            (fun i -> 
                let acc = ref 0.0 in
                for j = 0 to n-1 do
                    acc := !acc -. (hessian.(i).(j) *. gradient.(j))
                done;
                !acc) in
    (* main loop of algorithm *)
    let iter = ref 0 in
    let rec main_loop hessian f p fp step direction grad =
        incr iter;
        let np, nfp, direction, c = line_searcher f p fp grad step direction in
        if c || (!iter > max_iter) then begin
            np,nfp
        end else begin
            let grad, dgrad, hgrad, c = gradient_update hessian grad f np (get_score nfp) in
            debug_printf "\tNext Gradient: %s\n%!" (pp_farray grad);
            if c then begin
                np,nfp
            end else begin
                let hessian = bfgs_update_matrix dgrad hgrad direction hessian in
                let direction = calculate_direction hessian grad in
                debug_printf "\tNew Direction: %s\n%!" (pp_farray direction);
                main_loop hessian f np nfp step direction grad
            end
        end in
    (* initiate algorithm *)
    let hessian, pgrad, mxstep, dir = setup_function f p (get_score fp) in
    debug_printf "Initial Gradient: %s\n%!" (pp_farray pgrad);
    debug_printf "Initial Direction: %s\n%!" (pp_farray dir);
    let pf,fpf = main_loop hessian f p fp mxstep dir pgrad in
    debug_printf "Performed BFGS:\n\t(%s,%f)\n\t\t--[%d]-->\n\t(%s,%f)\n%!" (pp_farray p)
                    (get_score fp) !iter (pp_farray pf) (get_score fpf);
    (pf,fpf)


(* this assumes that the spec is consistent with the model itself, the returned
 * update function ensures this consistency *)
let get_update_function_for_model model =
  IFDEF USE_LIKELIHOOD THEN
    let add_alpha model_fun = match model.spec.site_variation with
        | Some (Gamma _) -> Some (update_all model_fun)
        | Some (Theta _) -> Some (update_alli model_fun)
        | None | Some Constant -> Some model_fun
    in
    match model.spec.substitution with
        | JC69 | F81 | File _ ->
            begin match model.spec.site_variation with
                | Some (Gamma _) -> Some (fun y x -> update_alpha y x.(0))
                | Some (Theta _) -> Some (fun y x -> update_alpha_invar y x.(0) x.(1))
                | None | Some Constant -> None
            end
        | TN93 _  -> add_alpha (fun y x -> update_tn93 y (x.(0),x.(1)) )
        | F84 _   -> add_alpha (fun y x -> update_f84 y x.(0))
        | GTR _   -> add_alpha update_gtr
        | K2P _   -> add_alpha (fun y x -> update_k2p y x.(0))
        | HKY85 _ -> add_alpha (fun y x -> update_hky y x.(0))
  ELSE
    None
  END

let get_current_parameters_for_model model =
  IFDEF USE_LIKELIHOOD THEN
    let make_single_with x = match model.spec.site_variation with
        | Some (Gamma (cat,alpha)) -> 
            begin match x with 
                | None   -> Some (Array.make 1 alpha)
                | Some x -> let y = Array.make 2 x in y.(1) <- alpha; Some y
            end
        | Some (Theta (cat,alpha,invar)) ->
            begin match x with 
                | None   -> Some (Array.make 1 alpha)
                | Some x -> let y = Array.make 2 x in y.(1) <- alpha; Some y
            end
        | None | Some Constant -> 
            begin match x with
                | None   -> None
                | Some x -> Some (Array.make 1 x)
            end
    in
    match model.spec.substitution with
        | JC69 | F81 | File _ -> make_single_with None
        | F84 x | K2P x | HKY85 x ->
            make_single_with (match x with | None -> Some default_tstv | _ -> x)
        | TN93 x  ->
            let x1,x2 = match x with 
                | Some (x,y) -> x,y
                | None -> default_tstv,default_tstv
            in
            begin match model.spec.site_variation with
                | Some (Gamma (cat,alpha)) ->
                    let y = Array.make 3 alpha in
                    y.(0) <- x1;
                    y.(1) <- x2;
                    Some y
                | Some (Theta (cat,alpha,invar)) ->
                    let y = Array.make 4 invar in
                    y.(0) <- x1;
                    y.(1) <- x1;
                    y.(2) <- alpha;
                    Some y
                | None | Some Constant -> 
                    let y = Array.make 2 x2 in
                    y.(0) <- x1;
                    Some y
            end
        | GTR x ->
            let x = match x with | Some x -> x | None -> default_gtr model.alph in
            begin match model.spec.site_variation with
                | Some (Gamma (cat,alpha)) ->
                    let y = Array.make (1 + Array.length x) alpha in
                    Array.blit x 0 y 0 (Array.length x);
                    Some y
                | Some (Theta (cat,alpha,invar)) ->
                    let y = Array.make (2 + Array.length x) alpha in
                    Array.blit x 0 y 0 (Array.length x);
                    y.(1 + Array.length x) <- invar;
                    Some y
                | None | Some Constant -> Some x
            end
  ELSE
    None
  END
