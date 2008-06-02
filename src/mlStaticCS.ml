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
let () = SadmanOutput.register "MlStaticCS" "$Revision %r $"

let debug = false 

(** caml links to garbage collection for deserialization **)
external register : unit -> unit = "likelihood_CAML_register"
let _ = register

type s
type cm = { (* character model *)
    name: string; (* JC69, F81 etcetera : for to_formatter output *)
    param: float array;     (* for to_formatter *)
    pi_0: (float, Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array1.t;
    (* for discrete distributions; gamma, invarient, custom... *)
    rate: (float, Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array1.t;
    prob: (float, Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array1.t;

    (* diagonalization of probability matrix   *)
    (* if ui == None then symmetric matrix..   *)
    u: (float, Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array2.t;
    d: (float, Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array2.t;
    ui:(float, Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array2.t option; }
type t = {
    code: int;
    mle: float;          (* -log likelihood of char set *)
    model: cm;
    codes: int array;
    chars: s; }

external diagonalize_gtr: (* U D Ui *)
    (float,Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array2.t ->
    (float,Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array2.t ->
    (float,Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array2.t ->
    unit = "likelihood_CAML_diagonalize_gtr"
external diagonalize_sym: (* U D *)
    (float,Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array2.t ->
    (float,Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array2.t ->
    unit = "likelihood_CAML_diagonalize_sym"
external compose_gtr: (* U D Ui *)
    (float,Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array2.t -> 
    (float,Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array2.t ->
    (float,Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array2.t -> float
    -> (float,Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array2.t =
        "likelihood_CAML_compose_gtr"
external compose_sym: (* U D *)
    (float,Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array2.t ->
    (float,Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array2.t -> float
    -> (float,Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array2.t =
        "likelihood_CAML_compose_sym"
external median_gtr: (* median_gtr U D Ui ta tb a b r p -> output_c *)
    (float,Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array2.t ->
    (float,Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array2.t ->
    (float,Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array2.t ->
     float -> float -> s -> s -> 
    (float,Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array1.t ->
    (float,Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array1.t -> s =
         "likelihood_CAML_median_gtr" "likelihood_CAML_median_wrapped_gtr" 
external median_sym: (* median sym U D ta tb a b r p -> output_c *)
    (float,Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array2.t ->
    (float,Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array2.t ->
     float -> float -> s-> s ->
    (float,Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array1.t ->
    (float,Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array1.t -> s = 
          "likelihood_CAML_median_sym" "likelihood_CAML_median_wrapped_sym"
external readjust_sym: (* readjust_sym U D a b c ta tb r p pi ll -> ll*ta*tb *)
    (float,Bigarray.float64_elt,Bigarray.c_layout) Bigarray.Array2.t ->
    (float,Bigarray.float64_elt,Bigarray.c_layout) Bigarray.Array2.t ->
    s -> s -> s -> float -> float -> 
    (float,Bigarray.float64_elt,Bigarray.c_layout) Bigarray.Array1.t ->
    (float,Bigarray.float64_elt,Bigarray.c_layout) Bigarray.Array1.t ->
    (float,Bigarray.float64_elt,Bigarray.c_layout) Bigarray.Array1.t ->
    float -> float*float*float = 
        "likelihood_CAML_readjust_sym" "likelihood_CAML_readjust_wrapped_sym"
external readjust_gtr:(* readjust_sym U D Ui a b c ta tb r p pi ll -> ll*ta*tb *)
    (float,Bigarray.float64_elt,Bigarray.c_layout) Bigarray.Array2.t ->
    (float,Bigarray.float64_elt,Bigarray.c_layout) Bigarray.Array2.t ->
    (float,Bigarray.float64_elt,Bigarray.c_layout) Bigarray.Array2.t ->
    s -> s -> s -> float -> float -> 
    (float,Bigarray.float64_elt,Bigarray.c_layout) Bigarray.Array1.t ->
    (float,Bigarray.float64_elt,Bigarray.c_layout) Bigarray.Array1.t ->
    (float,Bigarray.float64_elt,Bigarray.c_layout) Bigarray.Array1.t ->
    float -> float*float*float = 
        "likelihood_CAML_readjust_gtr" "likelihood_CAML_readjust_wrapped_gtr"

external s_bigarray: 
    s -> (float,Bigarray.float64_elt,Bigarray.c_layout) Bigarray.Array2.t =
    "likelihood_CAML_StoBigarray"
external bigarray_s: 
    (float,Bigarray.float64_elt,Bigarray.c_layout) Bigarray.Array2.t -> s =
    "likelihood_CAML_BigarraytoS"
external loglikelihood:
    s -> (float,Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array1.t ->
    float = "likelihood_CAML_loglikelihood"
external filter: s -> int array -> s = "likelihood_CAML_filter"
external compare_chars: s -> s -> int = "likelihood_CAML_compare"
external gamma_rates: float -> float -> int ->
    (float,Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array1.t = 
        "gamma_CAML_rates"

(* ------------------------------------------------------------------------- *)
(* printing functions *)
let print_float x = Printf.printf "%2.10f\t" x
let print_array xs = Array.iter print_float xs; print_newline ()
let print_array2 xxs = Array.iter print_array xxs; print_newline ()
let print_list x = List.map (fun y -> Printf.printf "%d\t" y) x
let print_barray a = 
    for i = 0 to (Bigarray.Array1.dim a)-1 do 
        Printf.printf "%2.10f\t" a.{i};
    done; Printf.printf "\n"; ()
let print_barray2 a = 
    for i = 0 to (Bigarray.Array2.dim1 a)-1 do 
        for j = 0 to (Bigarray.Array2.dim2 a)-1 do
            Printf.printf "%2.10f\t" a.{i,j};
        done; Printf.printf "\n"; 
    done; Printf.printf "\n"; ()

(* ------------------------------------------------------------------------- *)
(* conversion/utility functions *)

(* bigarray --> float array array *)
let barray_matrix bray =  
    let a = Bigarray.Array2.dim1 bray and b = Bigarray.Array2.dim2 bray in
    let r = Array.make_matrix a b 0.0 in
    for i = 0 to a-1 do for j = 0 to b-1 do
        r.(i).(j) <- bray.{i,j};
    done; done; r

(* bigarray1 to float array *)
let ba2array a = 
    let asize = Bigarray.Array1.dim a in
    let afray = Array.make asize a.{0} in
    for i = 1 to (asize-1) do
        afray.(i) <- a.{i};
    done; afray

(* copy elements of a into b *)
let barray_copy src des =
    assert( (Bigarray.Array1.dim src) <= (Bigarray.Array1.dim des));
    for i = 0 to (Bigarray.Array1.dim src)-1 do
        des.{i} <- src.{i}
    done; ()

(* ------------------------------------------------------------------------- *)
(* model creation functions *)

(* val jc69 :: ANY ALPHABET size *)
let m_jc69 mu a_size =
    let srm = Bigarray.Array2.create Bigarray.float64 Bigarray.c_layout a_size a_size in
    Bigarray.Array2.fill srm mu;
    let diag = -. mu *. float (a_size-1) in
    for i = 0 to (a_size-1) do
        srm.{i,i} <- diag
    done; srm

(* val k2p :: only 4 or 5 characters *)
let m_k2p alpha beta a_size =
    let srm = Bigarray.Array2.create Bigarray.float64 Bigarray.c_layout a_size a_size in
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
    done; srm

(* val tn93 :: only 4 or 5 characters *)
let m_tn93 pi_ alpha beta gamma a_size=
    let srm = Bigarray.Array2.create Bigarray.float64 Bigarray.c_layout a_size a_size in
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
    done; srm

(* val f81 :: ANY ALPHABET size *)
let m_f81 pi_ lambda a_size =
    let srm = Bigarray.Array2.create Bigarray.float64 Bigarray.c_layout a_size a_size in
    for i = 0 to (a_size-1) do
        for j = 0 to (a_size-1) do
            srm.{i,j} <- 
                if i = j then
                     -. lambda *. (1.0 -. pi_.(i))
                else 
                    pi_.(j)*.lambda;
        done;
    done; srm 

(* val hky85 :: only 4 or 5 characters *)
let m_hky85 pi_ alpha beta a_size = m_tn93 pi_ alpha alpha beta a_size

(* val f84 :: only 4 or 5 characters *)
let m_f84 pi_ gamma kappa a_size =
    let y = pi_.(1) +. pi_.(3) in (* Y = C + T *)
    let r = pi_.(0) +. pi_.(2) in (* R = A + G *)
    let alpha = (1.0+.kappa/.r) *. gamma in
    let beta = (1.0+.kappa/.y) *. gamma in
    m_tn93 pi_ alpha beta gamma a_size 

(* val gtr :: ANY ALPHABET size *)
(* pi_ == n; co_ == (1+n)*(n/2) *)
let m_gtr pi_ co_ a_size =
    let srm = Bigarray.Array2.create Bigarray.float64 Bigarray.c_layout a_size a_size in
    for i = 0 to (a_size-1) do
        for j = (i+1) to (a_size-1) do
            srm.{i,j} <- pi_.(j) *. co_.(i+j-1);
            srm.{j,i} <- pi_.(i) *. co_.(i+j-1);
        done;
    done;
    (* normalize diagonal so row sums to 0 *)
    for i = 0 to (a_size-1) do
        let diag = ref 0.0 in
        for j = 0 to (a_size-1) do
            if (i <> j) then diag := !diag +. srm.{i,j};
        done;
        srm.{i,i} <- -. !diag;
    done; srm

(* val m_file :: any alphabet size *)
(* re-computes diagonal --no verification against alphabet size*)
let m_file f_rr a_size =
    let srm = Bigarray.Array2.create Bigarray.float64 Bigarray.c_layout a_size a_size in
    for r = 0 to (a_size-1) do
        let diag = ref 0.0 in
        for c = 0 to (a_size-1) do
            if (c <> r) then 
                (diag := !diag +. f_rr.(r).(c); srm.{r,c} <- f_rr.(r).(c);)
        done;
        srm.{r,r} <- -. !diag;
    done; srm

(* ------------------------------------------------------------------------- *)
(* required functions *)

(* [median] calculate the new new node between [an] and [bn] with
 * distance [t1] + [t2], being applied to[an], [bn] respectively    *)
let median an bn t1 t2 =
    let am = an.model in
    let n_chars = match am.ui with
        | None -> 
            median_sym am.u am.d t1 t2 an.chars bn.chars am.rate am.prob
        | Some ui -> 
            median_gtr am.u am.d ui t1 t2 an.chars bn.chars am.rate am.prob in
    { an with
        chars = n_chars;
        mle = loglikelihood n_chars an.model.pi_0; 
    }

(* of_parser stuff *)
let rec list_of n x =
    match n with
    | 0 -> []
    | e -> x :: (list_of (n-1) x)
let rec set_in num lst =
    match lst with
    | hd :: tl ->
       (match num with
        | 0 -> 1.0::tl
        | n -> hd:: set_in (n-1) tl)
    | [] ->
        match num with
        | 0 -> [1.0]
        | e -> (list_of (e) 0.0) @ [1.0]
let rec sublist l a b =
    match l with
    | hd::tl ->
       (match a,b with
        | 0,0 -> []
        | 0,_ -> hd::sublist tl 0 (b-1)
        | _,_ -> sublist tl (a-1) b)
    | _ -> []
(* Parser.SC.static_spec -> ((int list option * int) array) -> t *)
let of_parser spec characters =
    let model = spec.Parser.SC.st_type in
    let model =
        match model with
        | Parser.SC.STLikelihood m -> m 
        | _ -> failwith "not a likelihood model" in 

    let (a_size,a_gap) = 
        match model.Parser.SC.use_gap with
        | true -> 
            ( Alphabet.size(Alphabet.to_sequential(spec.Parser.SC.st_alph)), -1 ) 
        | false ->
            ( Alphabet.size(Alphabet.to_sequential(spec.Parser.SC.st_alph))-1,
              Alphabet.get_gap( spec.Parser.SC.st_alph ) ) in

    let variation,probabilities =
        (* set up all the probability and rates *)
        match model.Parser.SC.site_variation with
        | None ->
           (Bigarray.Array1.of_array Bigarray.float64 Bigarray.c_layout ([| 1.0 |]),
            Bigarray.Array1.of_array Bigarray.float64 Bigarray.c_layout ([| 1.0 |]))
        | Some a -> ( match a with
            | Parser.SC.Invariant ->  (* same as NONE *)
                (Bigarray.Array1.of_array Bigarray.float64 Bigarray.c_layout ([| 1.0 |]),
                 Bigarray.Array1.of_array Bigarray.float64 Bigarray.c_layout ([| 1.0 |]))
            | Parser.SC.Gamma (x,y,z) -> (* SITES,ALPHA,BETA *)
                let p = Bigarray.Array1.create Bigarray.float64 Bigarray.c_layout x in
                Bigarray.Array1.fill p (1.0 /. (float_of_int x));
                (gamma_rates y z x, p )
            | Parser.SC.Theta (w,x,y,z) -> (* SITES,ALPHA,BETA,PERCENT_INVAR *)
                let p = Bigarray.Array1.create Bigarray.float64 Bigarray.c_layout (w+1) in
                let r = Bigarray.Array1.create Bigarray.float64 Bigarray.c_layout (w+1) in
                let rs = gamma_rates x y w in
                    (* [GAMMA@(1-z)]+[INVAR@z] *)
                    Bigarray.Array1.fill p ((1.0 -. z) /. (float_of_int w)); 
                    (* copy array to master array since one extra site @ 1.0 & z *)
                    barray_copy rs r; p.{w} <- z; r.{w} <- 1.0;
                (r, p)
            )
    in

    (* return the prior probability *)
    let priors =
        match model.Parser.SC.base_priors with
        | Parser.SC.Estimated p
        | Parser.SC.Given p -> 
            assert(a_size = Array.length p); p in

    (*  get the substitution rate matrix and set sym variable and to_formatter vars *)
    let sym = ref false in
    let mname = ref "none" in
    let (sub_mat,m_params) =
        match model.Parser.SC.substitution with
        | Parser.SC.JC69 lambda -> 
                sym := true;
                mname := "JC69";
                (m_jc69 lambda a_size, [| lambda |])
        | Parser.SC.K2P (alpha,beta) ->
                sym := true;
                mname := "K2P";
                (m_k2p alpha beta a_size, [|alpha;beta|])
        | Parser.SC.F81 lambda ->
                mname := "F81";
                (m_f81 priors lambda a_size, [| lambda |])
        | Parser.SC.F84 (kappa,beta) ->
                mname := "F84";
                (m_f84 priors kappa beta a_size, [| kappa;beta |])
        | Parser.SC.HKY85 (alpha,beta) ->
                mname := "HKY85";
                (m_hky85 priors alpha beta a_size, [| alpha;beta |])
        | Parser.SC.TN93 (alpha1,alpha2,beta) ->
                mname := "TN93";
                (m_tn93 priors alpha1 alpha2 beta a_size,[|alpha1;alpha2;beta|])
        | Parser.SC.GTR coeff -> 
                mname := "GTR"; (m_gtr priors coeff a_size, coeff)
        | Parser.SC.File matrix ->
                mname := "File"; (m_file matrix a_size, [| |])
    in
    (* diagonalize to get factored probability matrix --submat is destroyed
     * but can be reconstructed through the model parameters *)
    let (u_,d_,ui_) = 
        let n_d = Bigarray.Array2.create Bigarray.float64 Bigarray.c_layout a_size a_size in
        Bigarray.Array2.fill n_d 0.0;
        if !sym then (
            diagonalize_sym sub_mat n_d; (sub_mat, n_d, None) )
        else (
            let n_ui = Bigarray.Array2.create
                            Bigarray.float64 Bigarray.c_layout a_size a_size in
            Bigarray.Array2.fill n_ui 0.0;
            diagonalize_gtr sub_mat n_d n_ui; (sub_mat, n_d, Some n_ui) ) in

    (* loop to create array for each character *)
    let loop_ (states,code) =
        match states with 
        | None -> Array.make a_size 1.0
        | Some s -> 
            let lst = match s with
                | `List s -> s
                | `Bits s -> BitSet.to_list s
            in

            if List.mem a_gap lst then
                Array.make a_size 1.0 
            else
                let pl = List.fold_right set_in lst (list_of a_size 0.0) in
                Array.of_list (sublist pl 0 a_size)
    in

    (* convert character array to abstract type --redo *)
    let ba_chars = (Array.map loop_ characters) in (* create initial arrays *)
    let ba_chars = Bigarray.Array2.of_array Bigarray.float64 Bigarray.c_layout ba_chars in
    {
        code = model.Parser.SC.set_code;
         mle = 0.0;
       model = {
            rate = variation;
            prob = probabilities;
            name = !mname;
           param = m_params;
            pi_0 = Bigarray.Array1.of_array Bigarray.float64 Bigarray.c_layout priors;
           u = u_; d = d_;ui = ui_ };
       codes = Array.map (fun (x,y) -> y) characters; 
       chars = bigarray_s ba_chars;
    }

let f_abs x = if x < 0.0 then -.x else x
let pack lst = `Structured (`Set (List.map (fun x -> `Single x) lst))
(* Tags.attributes ->t ->t option ->Data.d *)
let to_formatter attr mine minet _ data :Tags.output list =
    (** get the alphabet **)
    let alpha = Data.get_alphabet data (Array.get mine.codes 0) in
    (** map functions over a set, into a single, then pack: with index, and with 2 lists **)
    let _ray f v = pack (Array.to_list (Array.map f v)) in
    let _ray2 f v c = pack (List.map2 f v c) in
    let _rayi f v = pack (Array.to_list (Array.mapi f v)) in
    (** pack individual elements into an ID **)
    let f_element c p :Tags.output =
        (Alphabet.find_code c alpha, [], `String (string_of_float p)) in
    (** pack a single character into an ID **)
    let f_char c cc:Tags.output = 
        let a_char = (Tags.Data.code, string_of_int cc) in
        (Tags.Characters.character, [a_char], (_rayi f_element c)) in
    (** pack the priors of the model **)
    let f_prior priors :Tags.output = let pi = ba2array priors in
        (Tags.Characters.prior, [], (_rayi f_element pi)) in
    (** function to pack a whole character set **)
    let f_chars ch_s co_s =
        let r = Array.to_list (barray_matrix ch_s) in
        _ray2 f_char r (Array.to_list co_s) in

    let attrib :Tags.attribute list = 
        (Tags.Data.code,string_of_int mine.code) :: 
        (Tags.Nodes.time, string_of_float minet) ::
        (Tags.Characters.mle, string_of_float mine.mle) :: attr in

    let con = pack (
                (Tags.Characters.model, (Tags.Data.modeltype, mine.model.name) :: [], 
                    pack ((f_prior mine.model.pi_0) :: [])
                ) :: 
                (Tags.Data.characters,[],
                    (f_chars (s_bigarray mine.chars) mine.codes)) :: []) in
    [(Tags.Characters.likelihood, attrib, con)]
(* -> Tags.output list *)

(* readjust the branch lengths to create better mle score *)
let readjust xopt x c1 c2 mine c_t1 c_t2 mine_t =
    let model = c1.model and mcpy = mine in
    let (nta,ntb,nl) = match model.ui with
        | None ->
            readjust_sym model.u model.d c1.chars c2.chars mcpy.chars c_t1
                    c_t2 model.rate model.prob model.pi_0 mine.mle
         | Some ui ->
            readjust_gtr model.u model.d ui c1.chars c2.chars mcpy.chars c_t1
                    c_t2 model.rate model.prob model.pi_0 mine.mle in
    if (nta = c_t1 && ntb = c_t2) then
        (x,mine.mle,mine.mle,(c_t1,c_t2),mine)
    else
        (* let _ = Printf.printf "Updating cost: %f --> %f\n" mine.mle nl in *)
        let x = Array.fold_right (* bottle neck? *)
                (fun c s -> All_sets.Integers.add c s) mine.codes x in
        (x,mine.mle,nl,(nta,ntb), {mine with mle = nl; chars = mcpy.chars} )

let distance a_node b_node t1 t2 =
    let t = median a_node b_node t1 t2 in t.mle

(* insert a node between two *)
(* (c, b) -> (c,(a),b)  *)
let dist_2 n a b nt at bt xt= 
    let x = median a b at bt in
    let tt = 0.0 in 
    let y = median x n tt nt in
    y.mle

let median_3 p x c1 c2  =  x
let reroot_median a b at bt = median a b at bt 
let median_cost ta = loglikelihood ta.chars ta.model.pi_0
let root_cost t = t.mle
let to_string a = "MLStaticCS"
let cardinal ta = Array.length ta.codes
let union prev ch1 ch2 = prev

let f_codes t codes = 
    let loopi_ i c = match (All_sets.Integers.exists (fun x -> x = c) codes) with
        | true -> Some i
        | false -> None in
    let loopOpt_ a = match a with
        | None -> false
        | _ -> true in
    let loopStrip_ a = match a with | Some i -> i | _ -> 0 in
    let opt_idx = Array.mapi (loopi_) t.codes in
    let opt_idx = Array.of_list (List.filter loopOpt_ (Array.to_list opt_idx)) in
    let opt_idx = Array.map loopStrip_ opt_idx in
    { t with chars = filter t.chars opt_idx }

let f_codes_comp t codes = 
    let loopi_ i c = match (All_sets.Integers.exists (fun x -> x = c) codes) with
        | true -> None
        | false -> Some i in
    let loopOpt_ a = match a with
        | None -> false
        | _ -> true in
    let loopStrip_ a = match a with | Some i -> i | _ -> 0 in
    let opt_idx = Array.mapi (loopi_) t.codes in
    let opt_idx = Array.of_list (List.filter loopOpt_ (Array.to_list opt_idx)) in
    let opt_idx = Array.map loopStrip_ opt_idx in
    { t with chars = filter t.chars opt_idx }

let compare_data a b = compare_chars a.chars b.chars

(* ------------------------------------------------------------------------- *)
(* Unused functions

let explode nd =
    let rec ba_lst ba i =
        if i = Bigarray.Array2.dim1 ba then []
        else (Bigarray.Array2.sub_left ba i 1) :: ba_lst ba (i+1)
    in
    (* convert chars to bigarray *) 
    let ba_s = s_bigarray nd.chars in
    (* combine codes and exploded bigarray, and convert bigarray back to s *)
    let stuff = List.map2  (fun x y -> (x,bigarray_s y))
                           (Array.to_list nd.codes)
                           (ba_lst ba_s 0) in
    List.map (fun (x,y) ->  {nd with codes = [|x|]; chars = y;code=nd.code;
                             mle = loglikelihood y nd.model.pi_0;})
             stuff *)
