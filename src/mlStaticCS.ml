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
let () = SadmanOutput.register "MlStaticCS" "$Revision: 3670 $"

let compress = true

let use_ocaml_readjust = false

IFDEF USE_LIKELIHOOD THEN
let failwithf format = Printf.ksprintf failwith format

(** caml links to garbage collection for deserialization **)
external register : unit -> unit = "likelihood_CAML_register"
let () = register ()

let (-->) a b = b a 
open Numerical.FPInfix

type s

type t = {
    mle     : float;
    model   : MlModel.model;
    weights : (float,Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array1.t;
    codes   : (int * int list) array;
    chars   : s;
}

let mem cs t = match cs with
    | None    -> true
    | Some [] -> false
    | Some xs -> 
        let rec found i = 
            if i = Array.length t.codes 
                then false
                else begin
                    if List.mem (fst t.codes.(i)) xs
                        then true
                        else found (i+1)
                end
        in
        found 0

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

external loglikelihood: (* vector, weight, priors, probabilities, %invar, mpl -> loglk *)
    s -> (float,Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array1.t
      -> (float,Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array1.t
      -> (float,Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array1.t
      -> float -> int -> float =
          "likelihood_CAML_loglikelihood" "likelihood_CAML_loglikelihood_wrapped"

external loglikelihood_site: (* vector, priors, probabilities, %invar, mpl, i -> loglk *)
    s -> (float,Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array1.t
      -> (float,Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array1.t
      -> float -> int -> int -> float =
          "likelihood_CAML_loglikelihood_site" "likelihood_CAML_loglikelihood_site_wrapped"

external rell_bootstrap: (* vector, weight, priors, probabilities, %invar, mpl, rel -> mean * var *)
    s -> (float,Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array1.t
      -> (float,Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array1.t
      -> (float,Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array1.t
      -> float -> int -> int -> int option -> float * float =
          "likelihood_CAML_rell_bootstrap" "likelihood_CAML_rell_bootstrap_wrapped"

external lk_variance_ratio :
    s -> s -> (float,Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array1.t
           -> (float,Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array1.t
           -> (float,Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array1.t
        -> float -> int -> float * float =
          "likelihood_CAML_variance_ratio" "likelihood_CAML_variance_ratio_wrapped"

external lk_variance :
    s -> (float,Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array1.t
      -> (float,Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array1.t
      -> (float,Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array1.t
      -> float -> int -> float * float =
          "likelihood_CAML_variance" "likelihood_CAML_variance_wrapped"

external filter: s -> int array -> s = "likelihood_CAML_filter"

external compare_chars: s -> s -> int = "likelihood_CAML_compare"

external debug: s -> unit = "likelihood_CAML_debug"

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
              MlModel.output_model print_string None `Nexus a.model None;
              Printf.printf "\nLikelihood: %f\n%!" (a.mle);
              ()

let cardinal ta = Array.length ta.codes

let union prev ch1 ch2 = prev

let get_codes a = Array.map fst a.codes

let get_model a = a.model

let set_model m a = {a with model = m; }

(* ------------------------------------------------------------------------- *)
(* initial estimation functions --jc69 *)
let min_bl = minimum_bl ()
let estimate_time a b =
    let r = float_of_int (MlModel.get_alphabet_size a.model) in
    let p = match (1.0 -. (proportion a.chars b.chars)) with
        | x when x < ((r-.1.0)/.r) -> x
        | x                      -> 0.50
    in
    let nt = ~-. ((r-.1.0)/.r) *. (log (1.0 -. (p *. (r/.(r-.1.0))))) in
    let nt = if nt <= min_bl then min_bl else nt /. 2.0 in
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
        loglikelihood n_chars an.weights am.MlModel.pi_0 
                      am.MlModel.prob pinvar (MlModel.get_costfn_code am)
    in
    let med = { an with chars = n_chars; mle = loglike; } in
    assert(loglike >= -0.0);
    med

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
    assert( loglike >= -0.0 );
    { an with
        chars = n_chars;
        mle = loglike;
    }

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
        loglikelihood n_chars an.weights am.MlModel.pi_0 am.MlModel.prob pinvar
                      (MlModel.get_costfn_code am)
    in
    assert( loglike >= -0.0 );
    { an with
        chars = n_chars;
        mle = loglike; 
    }

(* ------------------------------------------------------------------------- *)
(* Bootstrap / statistical confidence tests *)

let rell_bootstrap ?chars an num_replicates =
    let pinvar = match an.model.MlModel.invar with | Some x -> x | None -> ~-.1.0 in
    rell_bootstrap an.chars an.weights an.model.MlModel.pi_0 an.model.MlModel.prob
                pinvar (MlModel.get_costfn_code an.model) num_replicates chars

let variance_ratio an bn =
    let pinvar = match an.model.MlModel.invar with | Some x -> x | None -> ~-.1.0 in
    lk_variance_ratio an.chars bn.chars an.weights an.model.MlModel.pi_0
                an.model.MlModel.prob pinvar (MlModel.get_costfn_code an.model)

let variance an =
    let pinvar = match an.model.MlModel.invar with | Some x -> x | None -> ~-.1.0 in
    lk_variance an.chars an.weights an.model.MlModel.pi_0 an.model.MlModel.prob
                pinvar (MlModel.get_costfn_code an.model)

let site_likelihood an : (int * float * float) array =
    let m = an.model in
    let pinvar = match m.MlModel.invar with | Some x -> x | None -> ~-.1.0 in
    Array.init
        (Bigarray.Array1.dim an.weights)
        (fun i ->
            let s = 
                loglikelihood_site an.chars m.MlModel.pi_0 m.MlModel.prob
                                    pinvar (MlModel.get_costfn_code m) i in
            (an.codes.(i),an.weights.{i},~-.s))
    --> Array.to_list
    --> List.map
            (fun ((code,codes),weight,lk) ->
                assert (List.mem code codes);
                List.map (fun c -> (c,1.0,lk)) codes)
    --> List.flatten
    --> Array.of_list

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
            --> Array.make 1
            --> Bigarray.Array3.of_array Bigarray.float64 Bigarray.c_layout
    in
    let seq   = to_str_list seq in
    let nchar = List.length seq in
    let chars = load_characters (MlModel.get_alphabet model) seq in
    let wghts = Bigarray.Array1.of_array Bigarray.float64 Bigarray.c_layout 
                                         (Array.make nchar 1.0) in
    let schar = bigarray_s chars None  (MlModel.get_costfn_code model) in
    let loglk = 
        loglikelihood schar wghts model.MlModel.pi_0
                      model.MlModel.prob (-1.0) (MlModel.get_costfn_code model)
    in
    {   mle = loglk;
      model = model;
      codes = Array.init nchar (fun i -> (i,[i]));
    weights = wghts;
      chars = schar; }

(* Parser.SC.static_spec -> ((int list option * int) array) -> t *)
let of_parser spec weights characters =
    assert( (compress) ||
            (Array.fold_left (fun acc x -> acc && (x = 1.0)) true weights) );
    let computed_model = match spec.Nexus.File.st_type with
        | Nexus.File.STLikelihood x -> x
        | _ -> assert false
    in
    let (a_size,a_gap,u_gap) = 
        let alph = Alphabet.to_sequential (MlModel.get_alphabet computed_model) in
        match computed_model.MlModel.spec.MlModel.use_gap with
        | `Missing-> (Alphabet.size alph)-1,Alphabet.get_gap alph,false
        | `Independent | `Coupled _ -> Alphabet.size alph,Alphabet.get_gap alph,true
    in
    assert( a_size = (Bigarray.Array1.dim computed_model.MlModel.pi_0) );
    (* loop to create array for each character *)
    let loop_ (states,_) = match states with
        | None -> Array.make a_size 1.0
        | Some s ->
            let lst = Nexus.File.static_state_to_list s in
            let lst_minus_gap = List.filter (fun x -> not (x = a_gap)) lst in
            (*  Under a Missing gap model, if the gap character is set along
                with another character (ie, A/-) then we use the other chars
                minus the gap, else we set the state as we would a gap. *)
            if (not u_gap) && (List.mem a_gap lst) && ((List.length lst_minus_gap) = 0) then
                Array.make a_size 1.0
            else begin
                let lst = if (not u_gap) then lst_minus_gap else lst in
                Array.make a_size 0.0
                  --> List.fold_right (fun x acc -> acc.(x) <- 1.0;acc) lst
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
    let costfn = MlModel.get_costfn_code computed_model in
    let lk_chars = match computed_model.MlModel.invar with
        | Some _ -> bigarray_s ba_chars (Some aa_chars) costfn
        | None   -> bigarray_s ba_chars None costfn
    in
    let pinvar  = match computed_model.MlModel.invar with | Some x -> x | None -> ~-.1.0
    and weights = Bigarray.Array1.of_array Bigarray.float64 Bigarray.c_layout weights 
    and codes   = Array.map (snd) characters in
    assert( (Bigarray.Array1.dim weights) = (Bigarray.Array3.dim2 ba_chars));
    assert( (Array.length codes) = (Bigarray.Array3.dim2 ba_chars));
    let loglike =
        loglikelihood lk_chars weights computed_model.MlModel.pi_0
                      computed_model.MlModel.prob pinvar
                      (MlModel.get_costfn_code computed_model)
    in
    assert( loglike >= -0.0 );
    { mle  = loglike;
    model  = computed_model;
    codes  = codes;
    weights= weights;
    chars  = lk_chars; }

let to_formatter attr mine (t1,t2) data : Xml.xml Sexpr.t list =
    let str_time = function | Some x -> `Float x | None -> `String "None"
    and alphabet = MlModel.get_alphabet mine.model in
    let rec make_single_vec char_code single_ray =
        (Array.to_list
            (Array.mapi
                (fun state_code value ->
                    let alph = Alphabet.match_code state_code alphabet in
                    (PXML -[Xml.Characters.state]
                        ([Xml.Alphabet.element] = [`String alph])
                        ([Xml.Alphabet.value] = [`Float value])
                        { `Empty } --))
                single_ray))
    and make_single char_code single_ray =
        let name = try Data.code_character char_code data
                   with | _ -> failwithf "Cannot find code : %d" char_code in
        (PXML
            -[Xml.Characters.vector]
                ([Xml.Data.code] = [`Int char_code])
                ([Xml.Characters.name] = [`String name])
                { `Set (make_single_vec char_code single_ray) }
        --)
    in
    let sequence =
        let likelihood_vec,invariant_vec = s_bigarray mine.chars in
        (PXML
            -[Xml.Characters.characters]
            {
                let r = (barray_3matrix (likelihood_vec)).(0) in
                assert( (Array.length r) = (Array.length mine.codes) );
                `Set (List.map2 (make_single)
                                (List.map fst (Array.to_list mine.codes))
                                (Array.to_list r))
            }
        --)
    and model_data = MlModel.to_formatter mine.model in
    (PXML
        (* tag *)
        -[Xml.Characters.likelihood]
            (* attributes *)
            ([Xml.Characters.name] = [`String "-"])
            ([Xml.Characters.cost] = [`Float mine.mle])
            ([Xml.Nodes.min_time] = [str_time t1])
            ([Xml.Nodes.oth_time] = [str_time t2])
            ([attr])
            (* data *)
            { `Set (model_data @ [sequence]) }
        --) :: []
(* -> Xml.xml Sexpr.t list *)


(* ------------------------------------------------------------------------- *)
(* readjust the branch lengths to create better mle score *)
let readjust_c xopt x c1 c2 mine c_t1 c_t2 =
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
                (fun (c,_) s -> All_sets.Integers.add c s) new_mine.codes x in
        (x,new_mine.mle,nl,(nta,ntb), {new_mine with mle = nl;} )

let readjust_ocaml xopt x c1 c2 mine c_t1 c_t2 =
    (* Printf.printf "S: %f\t%f\t%f\n%!" c_t1 c_t2 new_mine.mle; *)
    let (nta,(new_mine,nl)) =
        let f nt1 = let x = median2 c1 c2 nt1 c_t2 0 0 in x, x.mle in
        Numerical.analyzer f (c_t1,f c_t1)
    and ntb = c_t2 in
    (* Printf.printf "E: %f\t%f\t%f\n%!" nta ntb nl; *)
    if nta =. c_t1 then
        (x,new_mine.mle,new_mine.mle,(c_t1,c_t2),new_mine)
    else
        let x = Array.fold_right (* bottle neck? *)
                (fun (c,_) s -> All_sets.Integers.add c s) new_mine.codes x in
        (x,new_mine.mle,nl,(nta,ntb), {new_mine with mle = nl;} )
    
let readjust = 
  if use_ocaml_readjust then readjust_ocaml else readjust_c

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

(** Resolve the assignment at a node to the be expected value from the median *)
let resolve ?(single=false) t =
    let comp,init = match t.model.MlModel.spec.MlModel.cost_fn with
        | `MPL -> max,(log 0.0)
        | `MAL -> max,(log 0.0)
    in
    let ray, _ = s_bigarray t.chars in (* ignore invar *)
    let nchars = Bigarray.Array3.dim2 ray
    and nrates = Bigarray.Array3.dim1 ray
    and nalpha = Bigarray.Array3.dim3 ray
    and priors = t.model.MlModel.pi_0
    and result = ref []
    and applyf = if single then List.hd else BitSet.Int.packed_of_list in
    for i = 0 to nchars - 1 do
        let best = ref init and state_ids = ref [] in
        for j = 0 to nrates - 1 do
            for k = 0 to nalpha - 1 do
                let v = (ray.{j,i,k}) *. (priors.{k}) in
                if v =. !best then begin
                    state_ids := k::(!state_ids);
                    best := comp !best v;
                end else if v > !best then begin
                    state_ids := [k];
                    best := v;
                end
            done;
        done;
        result := (applyf !state_ids)::(!result);
    done;
    Sequence.of_array (Array.of_list (List.rev !result))

(** Return the distance (the mle score) in joining two nodes *)
let distance a_node b_node t1 t2 =
    let t = median2 a_node b_node t1 t2 0 0 in t.mle

(** Return the cost of inserting a node between two: (c, b) -> (a,(n,b))  *)
let dist_2 n a b nt at bt xt= 
    let x = median2 a b at bt 0 0 in
    let tt = 0.5 in (* estimate time here *)
    let y = median2 x n tt nt 0 0 in
    y.mle

(** Median_3 is not necessary in static characters *)
let median_3 p x c1 c2 = x

let reroot_median a b at bt = median2 a b at bt 0 0

(** Median cost is the MLE score of the node *)
let median_cost ta = ta.mle

(** Process codes for striping out characters in the data-set *)
and process_codes comp node_codes codes =
    let loopi_ i (c,x) = 
        let isin = (All_sets.Integers.exists (fun x -> x = c) codes) in
        if isin = comp
            then Some (i,x)
            else None
    and loopOpt_ a = match a with
        | None   -> false
        | Some _ -> true
    and loopStrip_ a = match a with | Some i -> i | None -> assert false in
    node_codes 
        --> Array.mapi (loopi_)
        --> Array.to_list --> List.filter loopOpt_ --> Array.of_list
        --> Array.map loopStrip_


(** Filter out codes in the codes data-set *)
let f_codes_comp t codes =
    let filter chars opt_idx = filter chars (Array.map fst opt_idx) in
    let opt_idx = process_codes true t.codes codes in
    { t with chars = filter t.chars opt_idx;
             codes = process_codes false t.codes codes; }

(** Filter out codes NOT in the codes data-set *)
let f_codes t codes =
    let filter chars opt_idx = filter chars (Array.map fst opt_idx) in
    let opt_idx = process_codes false t.codes codes in
    { t with chars = filter t.chars opt_idx;
             codes = process_codes true t.codes codes; }

(** Compare to characters data *)
let compare_data a b = compare_chars a.chars b.chars

(** Compare the model of two characters; to determine if a median between them
    can be done; for consistency purposes *)
let compare a b = MlModel.compare a.model b.model

(* Does the tree in yang for testing; for decimal approximation. *)
let yang () =
    let model =
        let spec =
            let a = Alphabet.to_sequential Alphabet.nucleotides in
            {
                MlModel.alphabet = (a,4); MlModel.use_gap = `Missing;
                MlModel.site_variation = MlModel.Constant;
                MlModel.base_priors = MlModel.Equal; MlModel.cost_fn = `MAL;
                MlModel.substitution = MlModel.K2P 2.0;
            }
        in
        MlModel.create spec
    in
    let a = of_parser_simple "A" model and c = of_parser_simple "C" model
    (*and g = of_parser_simple "G" model*) and t = of_parser_simple "T" model in
    let tc  = median2 t   c  0.2 0.2 0 0 in
    let tca = median2 tc  a  0.1 0.2 0 0 in
    let cc  = median2 c   c  0.2 0.2 0 0 in
    let root= median2 tca cc 0.1 0.1 0 0 in
    Printf.printf "YANG TREE COST: %f\n" (root_cost root)

ELSE

let likelihood_error = 
    "Likelihood not enabled: download different binary or contact mailing list" 

let minimum_bl () = failwith likelihood_error

type t = unit

END

(** Return the NCM priors for the codes specified; here we assume the data has
    integrity, and that any subset of codes be sent here to obtain a prior cost.
    This is for ease of use; we don't categorize characters in node, since we
    need a cost for that function anyways. *)
let ncm_priors data codes =
    (* very analogous to the Data.assign_ncm_weight *)
    let rec get_prior_of_spec t = match t with
        | Nexus.File.STNCM (a,w,_) when a > 1 -> w *. (ncm_prior a)
        | Nexus.File.STSankoff _ | Nexus.File.STUnordered
        | Nexus.File.STOrdered   | Nexus.File.STLikelihood _ 
        | Nexus.File.STNCM _ -> 0.0
    (* the prior is, log (1 - 1/s)) *)
    and ncm_prior s =
        ~-. (log ((float_of_int (s-1)) /. (float_of_int s)))
    in
    (* calculate the overall-cost for static characters from codes *)
    let codes = match codes with
        | Some x -> x
        | None   -> Data.get_chars_codes_comp data `All
    in
    List.fold_left
        (fun acc code ->
            match Hashtbl.find data.Data.character_specs code with
            | Data.Static (Data.NexusFile spec) ->
                acc +. (get_prior_of_spec spec.Nexus.File.st_type)
            | Data.Dynamic _ | Data.Kolmogorov _ | Data.Set
            | Data.Static (Data.FixedStates _) -> acc)
        0.0
        codes

