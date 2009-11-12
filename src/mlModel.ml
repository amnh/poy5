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
let debug = true

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
    let default_gtr a   = Array.make a 1.0
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

exception Inconsistent_model
type update_model = 
    | Model  of subst_model
    | Rates  of site_var
    | Percent of float
    | Priors of priors

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

    assert( ((a_size*(a_size-1))/2) = Array.length co_ );
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


(* ------------------------------------------------ *)
(* CONVERSTION/MODEL CREATION FUNCTIONS             *)

(* convert a string spec to a specification, used in Parser for nexus *)
let convert_string_spec ((name,(var,site,alpha,invar),param,priors,gap,file):string_spec) =
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
    and priors = 
        let priors = List.map snd priors in 
        assert( 1.0 =. (List.fold_left (+.) 0.0 priors) );
        priors
    in
    {   substitution = submatrix;
        site_variation = Some variation;
        base_priors = Given (Array.of_list priors);
        use_gap = gap;
    }

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
    let sym, sub_mat = match lk_spec.substitution with
        | JC69  -> true,  m_jc69 priors 1.0 a_size
        | F81   -> false, m_f81 priors 1.0 a_size
        | K2P t -> 
            let t = match t with | Some x -> x | None -> default_tstv in
            true,  m_k2p priors t 1.0 a_size
        | F84 t -> 
            let t = match t with | Some x -> x | None -> default_tstv in
            false, m_f84 priors t 1.0 a_size
        | HKY85 t ->
            let t = match t with | Some x -> x | None -> default_tstv in
            false, m_hky85 priors t a_size
        (* more complex models *)
        | TN93 tstv ->
            let ts,tv = match tstv with 
                | Some (x,y) -> x,y 
                | None       -> default_tstv,default_tstv
            in
            false, m_tn93 priors ts tv 1.0 a_size
        | GTR c ->
            let c = match c with | Some xs -> xs | None -> default_gtr a_size in
            false, m_gtr priors c a_size
        | File m-> false, m_file priors m a_size
    in
    let (u_,d_,ui_) = diagonalize sym sub_mat in
    {
      rate = variation;
      prob = probabilities;
      spec = lk_spec;
     invar = invar;
      pi_0 = ba_of_array1 priors;
      alph = a_size;
         s = sub_mat;
         u = u_;
         d = d_;
        ui = ui_;
    }


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
    (* print alphabet and complements *)
    let () = 
        Printf.printf "Alphabet:\n";
        List.iter
            (fun (s,i) -> 
                let cs,ci = match Alphabet.complement i alph with
                    | Some x -> (Alphabet.match_code x alph),x
                    | None   -> "none",~-1
                in
                Printf.printf "\t%s -- %d -- %s -- %d\n" s i cs ci)
            (Alphabet.to_list alph)
    in

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
                        Printf.printf "Priors %s{%d} = %f\n%!" r b c;
                        c :: acc
                    end
                (* gaps are characters, seperate them -- minimum is 1/sum.
                 * essentially saying one character must exist, else a different
                 * alphabet should have been chosen *)
                end else begin
                    let c = try (All_sets.IntegerMap.find b pis) /. sum
                            with | Not_found -> 0.0 in
                    Printf.printf "Priors %s{%d} = %f\n%!" r b c;
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
            Printf.printf "S:%f\tV:%f\tA:%f\n" s v a;
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
                                        Printf.printf "Transition %s(%d)->%s(%d) = %f\n%!"
                                                        s1 alph1 s2 alph2 sum;
                                        sum :: acc2 
                                    end)
                                (Alphabet.to_list alph) acc1)
                        (Alphabet.to_list alph) []
            in
            (* we don't normalize the list, since in model creation the values
             * will be adjusted so the mean rate = 1.0 *)
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
        if debug then
            Printf.printf "Iteration %d: B:(%f,%f) L:(%f,%f) M:(%f,%f) U:(%f,%f)\n"
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
    if debug then
        Printf.printf "Iterated: %d\tImprovement: %f\tVariable: %f -> %f\n%!"
                  (iter_max - !iter) (orig_ll -. new_ll) orig_val new_val;
    new_

(* functions to update the models of evolution of a single (independent) parameter *)
let update_k2p old_model new_value =
    let subst_spec = { old_model.spec with substitution = K2P (Some new_value) }
    and priors = match old_model.spec.base_priors with | Estimated x | Given x -> x in
    let subst_model = m_k2p priors new_value 1.0 (Array.length priors) in
    let u,d,ui = diagonalize true subst_model in
    { old_model with spec = subst_spec; s  = subst_model; u  = u; d  = d; ui = ui; }

and update_hky old_model new_value =
    let subst_spec = { old_model.spec with substitution = HKY85 (Some new_value) }
    and priors = match old_model.spec.base_priors with | Estimated x | Given x -> x in
    let subst_model = m_hky85 priors new_value (Array.length priors) in
    let u,d,ui = diagonalize false subst_model in
    { old_model with spec = subst_spec; s  = subst_model; u  = u; d  = d; ui = ui; }

and update_f84 old_model new_value =  
    let subst_spec = { old_model.spec with substitution = F84 (Some new_value) }
    and priors = match old_model.spec.base_priors with | Estimated x | Given x -> x in
    let subst_model = m_f84 priors new_value 1.0 (Array.length priors) in
    let u,d,ui = diagonalize false subst_model in
    { old_model with spec = subst_spec; s  = subst_model; u  = u; d  = d; ui = ui; }

and update_alpha old_model new_value = 
    let new_spec_var,new_array = match old_model.spec.site_variation with
        | Some (Gamma (i,_)) -> Some (Gamma (i,new_value)), gamma_rates new_value new_value i
        | Some (Theta (i,_,p)) -> Some (Theta (i,new_value,p)), gamma_rates new_value new_value i
        | _ -> failwith "Model doesn't specify site variation"
    in
    { old_model with rate = new_array;
                     spec = { old_model.spec with site_variation = new_spec_var; } }
