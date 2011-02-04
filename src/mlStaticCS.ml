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

let compress = true

IFDEF USE_LIKELIHOOD THEN
let failwithf format = Printf.ksprintf failwith format

(** caml links to garbage collection for deserialization **)
external register : unit -> unit = "likelihood_CAML_register"
let () = register ()

let (-->) a b = b a 
let (=.) a b = abs_float (a-.b) < Numerical.epsilon 

type s

type t = {
    mle: float;
    model: MlModel.model;
    weights : (float,Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array1.t;
    codes   : int array;
    chars   : s;
}

external median1_gtr: (* median_gtr U D Ui t a b r p -> output_c *)
    FMatrix.m ->
    (float,Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array2.t ->
    (float,Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array2.t ->
    (float,Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array2.t ->
     float -> s -> s -> 
    (float,Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array1.t -> int -> s =
         "likelihood_CAML_median1_gtr" "likelihood_CAML_median1_wrapped_gtr" 
external median1_sym: (* median_sym U D t a b r p -> output_c *)
    FMatrix.m ->
    (float,Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array2.t ->
    (float,Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array2.t ->
     float -> s-> s ->
    (float,Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array1.t -> int -> s = 
          "likelihood_CAML_median1_sym" "likelihood_CAML_median1_wrapped_sym"

external median2_gtr: (* median_gtr U D Ui ta tb a b r p -> output_c *)
    FMatrix.m ->
    (float,Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array2.t ->
    (float,Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array2.t ->
    (float,Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array2.t ->
     float -> float -> s -> s -> 
    (float,Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array1.t -> int -> s =
         "likelihood_CAML_median2_gtr" "likelihood_CAML_median2_wrapped_gtr" 
external median2_sym: (* median_sym U D ta tb a b r p -> output_c *)
    FMatrix.m ->
    (float,Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array2.t ->
    (float,Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array2.t ->
     float -> float -> s-> s ->
    (float,Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array1.t -> int -> s = 
          "likelihood_CAML_median2_sym" "likelihood_CAML_median2_wrapped_sym"

external median3_gtr: (* median_gtr U D Ui ta tb tc a b c r p -> output_d *)
    FMatrix.m ->
    (float,Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array2.t ->
    (float,Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array2.t ->
    (float,Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array2.t ->
     float -> float -> float -> s -> s -> s -> 
    (float,Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array1.t -> int -> s =
         "likelihood_CAML_median3_gtr" "likelihood_CAML_median3_wrapped_gtr" 
external median3_sym: (* median_sym U D ta tb tc a b c r p -> output_d *)
    FMatrix.m ->
    (float,Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array2.t ->
    (float,Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array2.t ->
     float -> float -> float -> s-> s -> s ->
    (float,Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array1.t -> int -> s = 
          "likelihood_CAML_median3_sym" "likelihood_CAML_median3_wrapped_sym"

external readjust_sym: (* readjust_sym U D a b c ta tb %i r p pi ll -> ll*branch *)
    FMatrix.m ->
    (float,Bigarray.float64_elt,Bigarray.c_layout) Bigarray.Array2.t ->
    (float,Bigarray.float64_elt,Bigarray.c_layout) Bigarray.Array2.t ->
    s -> s -> s -> float -> float -> float ->
    (float,Bigarray.float64_elt,Bigarray.c_layout) Bigarray.Array1.t ->
    (float,Bigarray.float64_elt,Bigarray.c_layout) Bigarray.Array1.t ->
    (float,Bigarray.float64_elt,Bigarray.c_layout) Bigarray.Array1.t ->
    (float,Bigarray.float64_elt,Bigarray.c_layout) Bigarray.Array1.t ->
    float -> int -> float*float = 
        "likelihood_CAML_readjust_sym" "likelihood_CAML_readjust_sym_wrapped"
external readjust_gtr:(* readjust_sym U D Ui a b c ta tb %i r p pi ll -> ll*branch *)
    FMatrix.m ->
    (float,Bigarray.float64_elt,Bigarray.c_layout) Bigarray.Array2.t ->
    (float,Bigarray.float64_elt,Bigarray.c_layout) Bigarray.Array2.t ->
    (float,Bigarray.float64_elt,Bigarray.c_layout) Bigarray.Array2.t ->
    s -> s -> s -> float -> float -> float ->
    (float,Bigarray.float64_elt,Bigarray.c_layout) Bigarray.Array1.t ->
    (float,Bigarray.float64_elt,Bigarray.c_layout) Bigarray.Array1.t ->
    (float,Bigarray.float64_elt,Bigarray.c_layout) Bigarray.Array1.t ->
    (float,Bigarray.float64_elt,Bigarray.c_layout) Bigarray.Array1.t ->
    float -> int -> float*float = 
        "likelihood_CAML_readjust_gtr" "likelihood_CAML_readjust_gtr_wrapped"

external proportion: s -> s -> float = "likelihood_CAML_proportion"
external minimum_bl: unit -> float = "likelihood_CAML_minimum_bl"
external gc_alloc_max : int -> unit = "likelihood_GC_custom_max"
external copy : s -> s = "likelihood_CAML_copy"
external loglikelihood: (* vector, weight, priors, probabilities, and %invar -> loglk *)
    s -> (float,Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array1.t
      -> (float,Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array1.t
      -> (float,Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array1.t
      -> float -> int -> float =
          "likelihood_CAML_loglikelihood" "likelihood_CAML_loglikelihood_wrapped"
external filter: s -> int array -> s = "likelihood_CAML_filter"
external compare_chars: s -> s -> int = "likelihood_CAML_compare"

(* ------------------------------------------------------------------------- *)
(* printing functions *)
let print_barray1 a = 
    for i = 0 to (Bigarray.Array1.dim a)-1 do 
        Printf.printf "%2.10f\t" a.{i};
    done; Printf.printf "\n"; ()
let print_barray2 a = 
    for i = 0 to (Bigarray.Array2.dim1 a)-1 do 
        for j = 0 to (Bigarray.Array2.dim2 a)-1 do
            Printf.printf "%2.10f\t" a.{i,j};
        done; Printf.printf "\n"; 
    done; Printf.printf "\n"; ()
let print_barray3 a = 
    for i = 0 to (Bigarray.Array3.dim1 a)-1 do 
        for j = 0 to (Bigarray.Array3.dim2 a)-1 do
            for k = 0 to (Bigarray.Array3.dim3 a)-1 do
                Printf.printf "%2.10f\t" a.{i,j,k};
            done; Printf.printf "\n"; 
        done; Printf.printf "\n";
    done; Printf.printf "\n"; ()
let print_barray3with1 a b =
    for i = 0 to (Bigarray.Array3.dim1 a)-1 do 
        for j = 0 to (Bigarray.Array3.dim2 a)-1 do
            Printf.printf "%f:\t" b.{j};
            for k = 0 to (Bigarray.Array3.dim3 a)-1 do
                Printf.printf "%2.10f\t" a.{i,j,k};
            done; Printf.printf "\n"; 
        done; Printf.printf "\n";
    done; Printf.printf "\n%!"; ()

(* ------------------------------------------------------------------------- *)
(* conversion/utility functions *)

(* bigarray2 --> float array array *)
let barray_matrix bray =  
    let a = Bigarray.Array2.dim1 bray and b = Bigarray.Array2.dim2 bray in
    let r = Array.make_matrix a b 0.0 in
    for i = 0 to a-1 do for j = 0 to b-1 do
        r.(i).(j) <- bray.{i,j};
    done; done; r

(* bigarray3 --> float array array array *)
let barray_3matrix bray =  
    let b = Bigarray.Array3.dim2 bray and c = Bigarray.Array3.dim3 bray in
    Array.init
        (Bigarray.Array3.dim1 bray)
        (fun r -> 
            let x = Array.make_matrix b c 0.0 in
            for i = 0 to b-1 do for j = 0 to c-1 do
                x.(i).(j) <- bray.{r,i,j};
            done; done; x)

(* bigarray1 to float array *)
let ba2array a = 
    let asize = Bigarray.Array1.dim a in
    let afray = Array.make asize a.{0} in
    for i = 1 to (asize-1) do
        afray.(i) <- a.{i};
    done; afray

external s_bigarray: 
    s -> (float,Bigarray.float64_elt,Bigarray.c_layout) Bigarray.Array3.t *
    ((int32,Bigarray.int32_elt,Bigarray.c_layout) Bigarray.Array1.t) option =
    "likelihood_CAML_StoBigarray"
external bigarray_s: 
    (float,Bigarray.float64_elt,Bigarray.c_layout) Bigarray.Array3.t ->
    ((int32,Bigarray.int32_elt,Bigarray.c_layout) Bigarray.Array1.t) option ->
        int -> s = "likelihood_CAML_BigarraytoS"

(* ------------------------------------------------------------------------- *)
(* basic functions for interfacing with mlStatic module *)
let root_cost t = t.mle
let to_string _ = "MLStaticCS"
let print a = print_barray3 (fst (s_bigarray a.chars));
              MlModel.output_model print_string `Nexus a.model None;
              Printf.printf "\nLikelihood: %f\n%!" (a.mle)

let cardinal ta = Array.length ta.codes
let union prev ch1 ch2 = prev
let get_codes a = a.codes
let get_model a = a.model

(* ------------------------------------------------------------------------- *)
(* initial estimation functions --jc69 *)
let min_bl = minimum_bl ()
let estimate_time a b = 
    let p = match (1.0 -. (proportion a.chars b.chars)) with
        | x when x < 0.75 -> x
        | x -> 0.70
    in
    let nt2 = ~-. 0.75 *. (log (1.0 -. (p *. 4.0 /. 3.0))) in
    let nt = if nt2 <= min_bl then min_bl else nt2 /. 2.0 in
    (nt,nt)

(* ------------------------------------------------------------------------- *)
(* median functions *)
let median2 an bn t1 t2 acode bcode =
    let am = an.model in
    let n_chars = match am.MlModel.ui with
        | None -> 
            median2_sym FMatrix.scratch_space 
                am.MlModel.u am.MlModel.d t1 t2 an.chars bn.chars 
                am.MlModel.rate (MlModel.get_costfn_code am)
        | Some ui -> 
            median2_gtr FMatrix.scratch_space 
                am.MlModel.u am.MlModel.d ui t1 t2 an.chars bn.chars
                am.MlModel.rate (MlModel.get_costfn_code am)
    in
    let pinvar = match an.model.MlModel.invar with | Some x -> x | None -> ~-.1.0 in
    let loglike = 
        loglikelihood n_chars an.weights an.model.MlModel.pi_0 
                      an.model.MlModel.prob pinvar
                      (MlModel.get_costfn_code an.model)
    in
    assert( loglike >= 0.0 );
    { an with
        chars = n_chars;
        mle = loglike; 
    }

let median1 an bn t1 = 
    let am = an.model in
    let n_chars = match am.MlModel.ui with
        | None -> 
            median1_sym FMatrix.scratch_space 
                am.MlModel.u am.MlModel.d t1 an.chars bn.chars 
                am.MlModel.rate (MlModel.get_costfn_code am)
        | Some ui -> 
            median1_gtr FMatrix.scratch_space 
                am.MlModel.u am.MlModel.d ui t1 an.chars bn.chars
                am.MlModel.rate (MlModel.get_costfn_code am)
    in
    let pinvar = match am.MlModel.invar with | Some x -> x | None -> ~-.1.0 in
    let loglike = 
        loglikelihood n_chars an.weights am.MlModel.pi_0 
                      am.MlModel.prob pinvar (MlModel.get_costfn_code am)
    in
    assert( loglike >= 0.0 );
    { an with chars = n_chars; mle = loglike; }

let median3 an bn cn t1 t2 t3 =
    let am = an.model in
    let n_chars = match am.MlModel.ui with
        | None -> 
            median3_sym FMatrix.scratch_space 
                am.MlModel.u am.MlModel.d t1 t2 t3 an.chars bn.chars 
                cn.chars am.MlModel.rate (MlModel.get_costfn_code am)
        | Some ui -> 
            median3_gtr FMatrix.scratch_space 
                am.MlModel.u am.MlModel.d ui t1 t2 t3 an.chars bn.chars
                cn.chars am.MlModel.rate (MlModel.get_costfn_code am)
    in
    let pinvar = match an.model.MlModel.invar with | Some x -> x | None -> ~-.1.0 in
    let loglike = 
        loglikelihood n_chars an.weights an.model.MlModel.pi_0 
                      an.model.MlModel.prob pinvar
                      (MlModel.get_costfn_code an.model)
    in
    assert( loglike >= 0.0 );
    { an with
        chars = n_chars;
        mle = loglike; 
    }


(* ------------------------------------------------------------------------- *)
(* parser/formatter stuff *)
let rec list_of n x =
    match n with
    | 0 -> []
    | e -> x :: (list_of (e-1) x)
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
let farray_to_int32 x =
    let ipow =
        let rec ipow acc a =
            function | 0 -> acc | n -> ipow (a*acc) a (n-1)
        in ipow 1
    in
    Int32.of_int
        (Array_ops.fold_righti 
            (fun idx elm v ->
                if elm = 1.0
                    then (ipow 2 idx) + v
                    else v)
            0
            x)
let list_of_packed d =
    let rec loop_ c i d = match d land 1 with
        | 0 when d = 0 -> c
        | 0  -> loop_ c (i+1) (d lsr 1)
        | 1  -> loop_ (i::c) (i+1) (d lsr 1)
        | _  -> assert( false );
    in loop_ [] 0 d

let of_parser_simple seq model =
    let to_str_list string = 
        let rec to_ i len acc = 
            if i = len then List.rev acc
                       else to_ (i+1) len ((String.get string i)::acc)
        in
        List.map (String.make 1) 
                 (to_ 0 (String.length string) [])
    in
    let load_characters alph ss = 
        let alph = Alphabet.to_sequential alph in
        let gap  = Alphabet.get_gap alph
        and size,ugap = match model.MlModel.spec.MlModel.use_gap with
            | `Missing -> (Alphabet.size alph) - 1,false
            | `Independent | `Coupled _ -> Alphabet.size alph,true
        in
        let loop_ x =
            let elm = Alphabet.match_base x alph in
            if gap = elm && not ugap 
                then Array.make size 1.0 
                else list_of size 0.0 --> set_in elm --> Array.of_list
        in
        ss  --> List.map loop_
            --> Array.of_list
            --> Array.create 1
            --> Bigarray.Array3.of_array Bigarray.float64 Bigarray.c_layout
    in
    let seq   = to_str_list seq in
    let nchar = List.length seq in
    let chars = load_characters model.MlModel.alph seq in
    let wghts = Bigarray.Array1.of_array Bigarray.float64 Bigarray.c_layout 
                                         (Array.create nchar 1.0) in
    let schar = bigarray_s chars None  (MlModel.get_costfn_code model) in
    let loglk = 
        loglikelihood schar wghts model.MlModel.pi_0
                      model.MlModel.prob (~-.1.0)
                      (MlModel.get_costfn_code model)
    in
    {   mle = loglk;
      model = model;
      codes = Array.init nchar (fun i -> i);
    weights = wghts;
      chars = schar; }

(* Parser.SC.static_spec -> ((int list option * int) array) -> t *)
let of_parser spec weights characters =
    assert( (compress) ||
            (Array.fold_left (fun acc x -> acc && (x = 1.0)) true weights) );
    let computed_model = match spec.Nexus.File.st_type with
        | Nexus.File.STLikelihood x -> x
        | _ -> failwith "Not a likelihood model" in
    let () = (* ensure cost mode is acceptable: MAL/MPL *)
        match computed_model.MlModel.spec.MlModel.cost_fn with
        | `MPL | `MAL -> ()
        | `FLK | `ILK -> failwith "Cannot apply cost mode to static characters"
    in
    let (a_size,a_gap,u_gap) = 
        let alph = Alphabet.to_sequential computed_model.MlModel.alph in
        match computed_model.MlModel.spec.MlModel.use_gap with
        | `Missing-> (Alphabet.size alph)-1,Alphabet.get_gap alph,false
        | `Independent | `Coupled _ -> Alphabet.size alph,Alphabet.get_gap alph,true
    in
    (* loop to create array for each character *)
    let loop_ (states,code) = match states with 
        | None -> Array.make a_size 1.0
        | Some s -> 
            let lst = Nexus.File.static_state_to_list s in
            if (not u_gap) && (List.mem a_gap lst) then
                Array.make a_size 1.0 
            else begin
                list_of a_size 0.0
                    --> List.fold_right set_in lst
                    --> Array.of_list
            end
    in
    (* convert character array to abstract type --redo *)
    let aa_chars = Array.map loop_ characters in (* create initial arrays *)
    let ba_chars = (* construct Array3 dim1=rates,dim2=chars,dim3=states *)
        Array.init (Bigarray.Array1.dim computed_model.MlModel.rate)
                   (fun _ -> Array.copy aa_chars)
            --> Bigarray.Array3.of_array Bigarray.float64 Bigarray.c_layout
    and aa_chars =
        Array.map farray_to_int32 aa_chars -->
            Bigarray.Array1.of_array Bigarray.int32 Bigarray.c_layout
    in
    let lk_chars = match computed_model.MlModel.invar with
        | Some _ -> bigarray_s ba_chars (Some aa_chars) 
                        (MlModel.get_costfn_code computed_model)
        | None   -> bigarray_s ba_chars None
                        (MlModel.get_costfn_code computed_model)
    in
    let pinvar  = match computed_model.MlModel.invar with | Some x -> x | None -> ~-.1.0
    and weights = Bigarray.Array1.of_array Bigarray.float64 Bigarray.c_layout weights 
    and codes   = Array.map (fun (x,y) -> y) characters in
    assert( (Bigarray.Array1.dim weights) = (Bigarray.Array3.dim2 ba_chars));
    assert( (Array.length codes) = (Bigarray.Array3.dim2 ba_chars));
    let loglike = 
        loglikelihood lk_chars weights computed_model.MlModel.pi_0
                      computed_model.MlModel.prob pinvar
                      (MlModel.get_costfn_code computed_model)
    in
    (* print_barray3with1 ba_chars weights; *)
    assert( loglike >= 0.0 );
    {    mle  = loglike;
       model  = computed_model;
       codes  = codes;
       weights= weights;
       chars  = lk_chars; }

let to_formatter attr mine (t1,t2) data : Xml.xml Sexpr.t list =
    let str_time = function | Some x -> `Float x | None -> `String "None"
    and alphabet = mine.model.MlModel.alph in
    let rec make_single_vec char_code single_ray =
        (Array.to_list 
            (Array.mapi 
                (fun state_code value ->
                    let alph = Alphabet.match_code state_code alphabet in
                    (PXML -[Xml.Characters.vector]
                        ([Xml.Alphabet.value] = [`Float value])
                        {`String alph} --))
                single_ray))

    and make_single char_code single_ray = 
        (PXML
            -[Xml.Characters.likelihood]
                ([Xml.Data.code] = [`Int char_code])
                { `Set (make_single_vec char_code single_ray) }
        --)
    in
    let sequence =
        let likelihood_vec,invariant_vec = s_bigarray mine.chars in
        (PXML
            -[Xml.Characters.characters]
            {  
                let r = Array.to_list (barray_3matrix (likelihood_vec)).(0) in
                `Set (List.map2 (make_single) (Array.to_list mine.codes) r)
            }
        --)
    and model_data = MlModel.to_formatter mine.model in
    (PXML
        (* tag *)
        -[Xml.Characters.likelihood]
            (* attributes *)
            ([Xml.Characters.llike] = [`Float mine.mle])
            ([Xml.Nodes.min_time] = [str_time t1])
            ([Xml.Nodes.oth_time] = [str_time t2])
            ([attr])
            (* data *)
            { `Set (model_data @ [sequence]) }
        --) :: []
(* -> Xml.xml Sexpr.t list *)


(* ------------------------------------------------------------------------- *)
(* readjust the branch lengths to create better mle score *)
let readjust xopt x c1 c2 mine c_t1 c_t2 =
    (* copy characters to a new set *)
    let new_mine = {mine with chars = copy mine.chars} in
    let model = c1.model in
    let pinv  = match model.MlModel.invar with | Some x -> x | None -> ~-.1.0 in
    (* Printf.printf "S: %f\t%f\t%f\n%!" c_t1 c_t2 new_mine.mle; *)
    let (nta,nl) = match model.MlModel.ui with
        | None ->
            readjust_sym FMatrix.scratch_space model.MlModel.u model.MlModel.d 
                         c1.chars c2.chars new_mine.chars c_t1 c_t2 pinv
                         c1.weights model.MlModel.rate model.MlModel.prob
                         model.MlModel.pi_0 new_mine.mle (MlModel.get_costfn_code model)

        | Some ui ->
            readjust_gtr FMatrix.scratch_space model.MlModel.u
                         model.MlModel.d ui c1.chars c2.chars new_mine.chars
                         c_t1 c_t2 pinv c1.weights model.MlModel.rate
                         model.MlModel.prob model.MlModel.pi_0 new_mine.mle
                         (MlModel.get_costfn_code model)
    and ntb = c_t2 in
    (* Printf.printf "E: %f\t%f\t%f\n%!" nta ntb nl; *)
    if nta =. c_t1 then
        (x,new_mine.mle,new_mine.mle,(c_t1,c_t2),new_mine)
    else
        let x = Array.fold_right (* bottle neck? *)
                (fun c s -> All_sets.Integers.add c s) new_mine.codes x in
        (x,new_mine.mle,nl,(nta,ntb), {new_mine with mle = nl;} )

(* ------------------------------------------------------------------------- *)
(* extract maximum state from all the characters *)
let extract_states a_node =
    let ray, _ = s_bigarray a_node.chars in (* ignore invar *)
    let nchars = Bigarray.Array3.dim2 ray
    and nrates = Bigarray.Array3.dim1 ray
    and nalpha = Bigarray.Array3.dim3 ray
    and priors = a_node.model.MlModel.pi_0 in
    let result = ref [] in
    for i = 0 to nchars - 1 do
        let highest = ref 0.0 and state_ids = ref [] in
        for j = 0 to nrates - 1 do
            for k = 0 to nalpha - 1 do
                let v = (ray.{j,i,k}) *. (priors.{k}) in
                if v =. !highest then begin
                    state_ids := k::(!state_ids);
                    highest := max !highest v;
                end else if v > !highest then begin
                    state_ids := [k];
                    highest := v;
                end
            done;
        done;
        result := (a_node.weights.{i},i,`List !state_ids)::(!result);
    done;
    !result

let distance a_node b_node t1 t2 = (* codes don't matter here *)
    let t = median2 a_node b_node t1 t2 0 0 in t.mle

(* insert a node between two *)
(* (c, b) -> (c,(a),b)  *)
let dist_2 n a b nt at bt xt= 
    let x = median2 a b at bt 0 0 in
    let tt = 0.5 in (* estimate time here *)
    let y = median2 x n tt nt 0 0 in
    y.mle

let median_3 p x c1 c2  =  x
let reroot_median a b at bt = median2 a b at bt 0 0
let median_cost ta = ta.mle (*
    let pinvar = match ta.model.invar with | Some x -> x | None -> ~-.1.0 in
    loglikelihood ta.chars ta.model.pi_0 ta.model.prob pinvar *)

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

let compare a b = MlModel.compare a.model b.model

ELSE

let likelihood_error = 
    "Likelihood not enabled: download different binary or contact mailing list" 

type t = unit
END
