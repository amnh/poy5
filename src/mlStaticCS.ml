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

IFDEF USE_LIKELIHOOD THEN
let pure_ocaml   = false    (* ONLY use pure ocaml version *)
let graph_output = false    (* graph all medians in %d--%d format *)
let phyml_mode   = true     (* divide by the mean rate *)

(** caml links to garbage collection for deserialization **)
external register : unit -> unit = "likelihood_CAML_register"
let () = register ()

let (-->) a b = b a 

let epsilon = 0.000001
let (=.) a b = abs_float (a-.b) < epsilon (*
        match classify_float ( a -. b ) with
        | FP_subnormal | FP_zero -> true
        | FP_infinite | FP_nan | FP_normal -> false *)

let scratch_space = FMatrix.create 20

type s
type cm = { (* character model *)
    name: string; (* JC69, F81 etcetera : for to_formatter output *)
    pi_0: (float, Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array1.t;
    (* for discrete distributions; gamma, invarient, custom... *)
    alpha : float option; invar : float option; sites: int;
    rate: (float, Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array1.t;
    prob: (float, Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array1.t;
    (* diagonalization of probability matrix   *)
    (* if ui == None then symmetric matrix..   *)
    u: (float, Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array2.t;
    d: (float, Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array2.t;
    ui:(float, Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array2.t option; }

type t = {
    mle: float;
    model: cm;
    codes: int array;
    chars: s;
}

external diagonalize_gtr: (* U D Ui *) FMatrix.m ->
    (float,Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array2.t ->
    (float,Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array2.t ->
    (float,Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array2.t ->
    unit = "likelihood_CAML_diagonalize_gtr"
external diagonalize_sym: (* U D *) FMatrix.m ->
    (float,Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array2.t ->
    (float,Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array2.t ->
    unit = "likelihood_CAML_diagonalize_sym"
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
external median_gtr: (* median_gtr U D Ui ta tb a b r p -> output_c *)
    FMatrix.m ->
    (float,Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array2.t ->
    (float,Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array2.t ->
    (float,Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array2.t ->
     float -> float -> s -> s -> 
    (float,Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array1.t ->
    (float,Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array1.t -> s =
         "likelihood_CAML_median_gtr" "likelihood_CAML_median_wrapped_gtr" 
external median_sym: (* median sym U D ta tb a b r p -> output_c *)
    FMatrix.m ->
    (float,Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array2.t ->
    (float,Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array2.t ->
     float -> float -> s-> s ->
    (float,Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array1.t ->
    (float,Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array1.t -> s = 
          "likelihood_CAML_median_sym" "likelihood_CAML_median_wrapped_sym"
external readjust_sym: (* readjust_sym U D a b c ta tb %i r p pi ll -> ll*branch *)
    FMatrix.m ->
    (float,Bigarray.float64_elt,Bigarray.c_layout) Bigarray.Array2.t ->
    (float,Bigarray.float64_elt,Bigarray.c_layout) Bigarray.Array2.t ->
    s -> s -> s -> float -> float -> float ->
    (float,Bigarray.float64_elt,Bigarray.c_layout) Bigarray.Array1.t ->
    (float,Bigarray.float64_elt,Bigarray.c_layout) Bigarray.Array1.t ->
    (float,Bigarray.float64_elt,Bigarray.c_layout) Bigarray.Array1.t ->
    float -> float*float = 
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
    float -> float*float = 
        "likelihood_CAML_readjust_gtr" "likelihood_CAML_readjust_gtr_wrapped"

external proportion: s -> s -> float = "likelihood_CAML_proportion"
external minimum_bl: unit -> float = "likelihood_CAML_minimum_bl"

external loglikelihood: (* vector, priors, probabilities, and %invar -> loglk *)
    s -> (float,Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array1.t ->
    (float,Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array1.t ->
    float -> float = "likelihood_CAML_loglikelihood"
external filter: s -> int array -> s = "likelihood_CAML_filter"
external compare_chars: s -> s -> int = "likelihood_CAML_compare"
external gamma_rates: float -> float -> int ->
    (float,Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array1.t = 
        "gamma_CAML_rates"

(* ------------------------------------------------------------------------- *)
(* printing functions *)
let print_float x = Printf.printf "%2.10f\t" x
let print_int x = Printf.printf "%d\t" x
let print_array xs = Array.iter print_float xs; print_newline ()
let print_arrayi xs = Array.iter print_int xs; print_newline ()
let print_array2 xxs = Array.iter print_array xxs; print_newline ()
let print_array3 xxxs = Array.iter print_array2 xxxs; print_newline ()
let print_list x = List.map (fun y -> Printf.printf "%d\t" y) x
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

let pp_fopt chan v = 
    output_string chan 
                  (match v with
                    | Some v -> string_of_float v
                    | None -> "none")

let get_codes a = a.codes

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
    ((int32,Bigarray.int32_elt,Bigarray.c_layout) Bigarray.Array1.t) option -> s =
    "likelihood_CAML_BigarraytoS"

(* ------------------------------------------------------
 *  A pure ocaml implementation of the median functions
*)
let array_map2 f a1 a2 =
    let x = Array.make (Array.length a1) (f a1.(0) a2.(0)) in
        for i = 1 to ((Array.length a1)-1) do
            x.(i) <- f a1.(i) a2.(i)
        done; x

(* sum of product of two vectors *)
let dot_product v1 v2 = 
    let r = array_map2 ( *. ) v1 v2 in
    let res = Array.fold_right (+.) r 0.0 in
    res

(* negative log liklihood *)
let mle a priors probs =
    let total = ref 0.0 in
    Array.iteri
        (fun i _ ->
            let res x = -. log (dot_product x priors) in
            let curr = Array.fold_right (+.) (Array.map res a.(i)) 0.0 in
            total := !total +. (probs.(i) *. curr)
        )
        a;
    !total

(* calculate a new nodes mle and prob_vectors *)
let median_char p_1 p_2 a b = 
    (* calculates one element of the probability vector *)
    let median_element a b p1 p2 x = 
        let x1 = Array.get p1 x and x2 = Array.get p2 x in
        let right_sum = dot_product a x1 and lefts_sum = dot_product b x2 in
        lefts_sum *. right_sum
    in
    let curried_c = median_element a b p_1 p_2 in
    let npv = Array.init (Array.length a) curried_c in
    npv

(* empty argument is 'previous' *)
let median_fmat a_vec b_vec a_mat b_mat =
    array_map2 (median_char a_mat b_mat) a_vec b_vec

(** --- SEARCH METHODS --- **)

(* search by brent's method *)
exception Colinear
let decr x = decr x; !x (* get and decrement, because it makes more sense *)
let golden_middle a b = 
    let a,b = if a < b then a,b else b,a in
    a +. ((b -. a) *. 2.0 /. (1.0 +. sqrt 5.0))
let golden_exterior a b = (* when fb < fa *)
    a +. ((b -. a) *. 2.0 /. ((sqrt 5.0) -. 1.0))
let abscissa (a,fa) (b,fb) (c,fc) =
    let numer = ((b-.a)*.(b-.a)*.(fb-.fc)) -. ((b-.c)*.(b-.c)*.(fb-.fa))
    and denom = ((b-.a)*.(fb-.fc)) -. ((b-.c)*.(fb-.fa)) in
    if denom = 0.0 then raise Colinear
    else b -. (0.5 *. (numer /. denom))

(* function and braketed area and error *)
let brents_method ((orig_bl,orig_ll) as orig) f lower upper epsilon = 
    let iter = ref 1000 in
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

    and best_of ((_,x) as x') ((_,y) as y') = if x <= y then x' else y' in

    (* parabolic interpolation *)
    let rec parabolic_interp ((best_t,best_l) as best) a fa b fb c fc : float * float = 
        assert( a > 0.0 && b > 0.0 && c > 0.0 );
        if ((abs_float (fb -. fa)) <= epsilon) or 
           ((abs_float (fc -. fa)) <= epsilon) or (decr iter) = 0 then best
        else
            try
                let x = abscissa (a,fa) (b,fb) (c,fc) in let fx = f x in
                if fx =. fb || x =. b || x <= epsilon then best
                else begin
                    let best = best_of best (x,fx) in
                    if a < x && x < b then parabolic_interp best a fa x fx b fb
                    else if x < a then parabolic_interp best x fx a fa b fb
                    else if b < x && x < c then parabolic_interp best b fb x fx c fc
                    else if c < x then parabolic_interp best b fb c fc x fx
                    else failwith "shouldn't happen"
                end
            with | Colinear -> brent_decision best a fa b fb c fc

    (* golden section search, when function is crappy --this isn't needed when
     * there is only one minima, as in branch lengths *)
    and golden_ratio ((best_t,best_l) as best) a af b bf c cf  =
        assert( a > 0.0 && b > 0.0 && c > 0.0 );
        let best,(a,fa),(nb,nfb),(c,fc) = 
            if (abs_float (c-.b)) > (abs_float (b-.a)) then
                let other = golden_middle b c in let other = (other,f other) in
                let best = best_of best other in
                best,(b,bf),other,(c,cf)
            else 
                let other = golden_middle a b in let other = (other,f other) in
                let best = best_of best other in
                best,(a,af),other,(b,bf)
        in
        if bf =. nfb or (decr iter) = 0 then best
        else brent_decision best a fa nb nfb c fc

    (* brent exponential search, when points are colinear or monotonic
     *  does not return a result since we are widening the search area. *)
    and brent_exp ((best_t,best_l) as best) a fa b fb c fc : float * float =
        assert (fc < fa ); (* since we estimate "past" c *)
        let n = match golden_exterior a c with
            | x when x <= 0.0 -> epsilon
            | x -> x 
        in
        let other = n,f n in
        brent_decision (best_of best other) b fb c fc n (snd other)
    (**
     * What to do, what to do? Well, if the three points are monotonic, then
     * call brent_exp until something better comes up, if there is a minimum
     * between we can do a parabolic interpolation, else we use a shitty golden
     * bisect search method each iteration to find a better spot.
    *) 
    and brent_decision best lower fl middle fm upper fu : float * float =
        let (lower,fl),(middle,fm),(upper,fu) = 
            order_triples (lower,fl) (middle,fm) (upper,fu) in
        if lower =. upper then best
        else if fl <= fm && fm <= fu then (* monotonically increasing *)
            brent_exp best upper fu middle fm lower fl
        else if fl >= fm && fm >= fu then (* monotonically decreasing *)
            brent_exp best lower fl middle fm upper fu
        else if fm <= fl && fm <= fu then (* minimum between *)
            parabolic_interp best lower fl middle fm upper fu
        else begin
            (* let step = 0.0001 in
            let stepval xref = xref := !xref +. step;!xref
            and time = ref (lower -. step) and out_chan = open_out "curve.tsv" in
            Printf.printf "Making file 'curve.tsv': %f -- %f\n" lower upper;
            while (stepval time) < upper do
                Printf.fprintf out_chan "%f\t%f\n" !time (f !time);
            done;
            close_out out_chan; *)
            golden_ratio best lower fl middle fm upper fu
        end
    in
    (* set up variables.. order arguments,find golden middle and evaluate *)
    assert( lower != upper );
    let middle = golden_middle lower upper in
    let fl = f lower and fm = f middle and fu = f upper in
    let best = best_of orig (best_of (lower,fl) (best_of (middle,fm) (upper,fu))) in
    let ((t,ll) as x) = brent_decision best lower fl middle fm upper fu in
    Printf.printf "Iterated: %d\tImprovement: %f\tBranch: %f -> %f\n%!"
                  (1000 - !iter) (orig_ll -. ll) orig_bl t;
    x

(**
 * converts the stored variables into float array/matrices and computes mle
*)
let ocaml_median (a:t) (b:t) (acode:int) (bcode:int) (t1:float) (t2:float) =
    let make_matrix model t = 
        barray_matrix
            (match model.ui with
            | Some ui -> compose_gtr scratch_space model.u model.d ui t
            | None -> compose_sym scratch_space model.u model.d t)
    in
    (* convert abstract types to C types *)
    let ach,ain = s_bigarray a.chars and bch,bin = s_bigarray b.chars in
    (* convert to ocaml primatives *)
    let ach = barray_3matrix ach and bch = barray_3matrix bch 
    (** TODO :: invar **)
    and pi_ = ba2array (a.model.pi_0)
    and prob= ba2array (a.model.prob) in
    let root =
        Array.mapi
            (fun i _ -> 
                median_fmat ach.(i) bch.(i) 
                            (make_matrix a.model (t1 *. a.model.rate.{i}) )
                            (make_matrix b.model (t2 *. b.model.rate.{i}) ) )
            ach (* arbitrary, as long as it's the same length *)
    in
    root, mle root pi_ prob

let ocaml_graph (a:t) (b:t) (min:float) (max:float) (step:float) (f:string) = 
    let stepval xref = xref := !xref +. step;!xref
    and time = ref (min -. step) and out_chan = open_out f in
    while (stepval time) < max do
        let c_time = !time /. 2.0 in
        let _,ll = ocaml_median a b 0 0 c_time c_time in
        Printf.fprintf out_chan "%f\t%f\n" !time ll;
    done;
    close_out out_chan;
    ()

(**
 * converts the stored variables into float array/matrices and adjusts branches *)
let ocaml_readjust (a:t) (b:t) (t1:float) (t2:float) (b_ll:float) : float * float * float =
    let dist = t1 +. t2 in
    let median_2 t1 t2 = let _,ll = ocaml_median a b 0 0 t1 t2 in ll in
    let t,ll = brents_method (dist,b_ll)
                             (fun x -> let half = x /. 2.0 in median_2 half half)
                             (dist /. 10.0) (dist *. 1.5) epsilon
    in
    let new_halves = t /. 2.0 in
    new_halves,new_halves,ll

(* ------------------------------------------------------------------------- *)
(* model creation functions *)

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
    let srm = Bigarray.Array2.create Bigarray.float64 Bigarray.c_layout a_size a_size in
    Bigarray.Array2.fill srm mu;
    let diag = -. mu *. float (a_size-1) in
    for i = 0 to (a_size-1) do
        srm.{i,i} <- diag
    done;
    if phyml_mode then m_meanrate srm pi_;
    srm

(* val k2p :: only 4 or 5 characters *)
let m_k2p pi_ alpha beta a_size =
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
    done;
    if phyml_mode then m_meanrate srm pi_;
    srm

(* val tn93 :: only 4 or 5 characters *)
let m_tn93 pi_ alpha beta gamma a_size =
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
    done;
    if phyml_mode then m_meanrate srm pi_;
    srm

let m_tn93_ratio pi_ kappa1 kappa2 a_size = 
    let beta = (pi_.(0) *. pi_.(2) *. kappa1) +. (pi_.(1) *. pi_.(3) *. kappa2) +.
                ((pi_.(0) +. pi_.(2)) *. (pi_.(1)+.pi_.(3)) ) in
    let beta = 1.0 /. (2.0 *. beta) in
    let alpha1 = kappa1 *. beta and alpha2 = kappa2 *. beta in
    m_tn93 pi_ alpha1 alpha2 beta a_size

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
    done;
    if phyml_mode then m_meanrate srm pi_;
    srm

(* val hky85_ratio :: only 4 or 5 characters *)
let m_hky85_ratio pi_ kappa a_size =
    let beta = (pi_.(0) *. pi_.(2) *. kappa) +. (pi_.(1) *. pi_.(3) *. kappa) +.
                ((pi_.(0) +. pi_.(2)) *. (pi_.(1)+.pi_.(3)) ) in
    let beta = 1.0 /. (2.0 *. beta) in
    let alpha = kappa *. beta in
    m_tn93 pi_ alpha alpha beta a_size

(* val hky85 :: only 4 or 5 characters *)
let m_hky85 pi_ alpha beta a_size = m_tn93 pi_ alpha alpha beta a_size

(* val f84 :: only 4 or 5 characters *)
let m_f84 pi_ gamma kappa a_size =
    let y = pi_.(1) +. pi_.(3) in (* Y = C + T *)
    let r = pi_.(0) +. pi_.(2) in (* R = A + G *)
    let alpha = (1.0+.kappa/.r) *. gamma in
    let beta = (1.0+.kappa/.y) *. gamma in
    m_tn93 pi_ alpha beta gamma a_size 

(* val gtr :: ANY ALPHABET size
   pi_ == n; co_ == ((1+n)*n)/2
   form of lower packed storage mode, excluding diagonal, *)
let m_gtr pi_ co_ a_size =
    let n = ref 0 in (* array index *)
    let srm = Bigarray.Array2.create Bigarray.float64 Bigarray.c_layout a_size a_size in
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
    if phyml_mode then m_meanrate srm pi_;
    srm

(* val m_file :: any alphabet size -- recomputes diagonal and divides by meanrate *)
let m_file pi_ f_rr a_size =
    assert(a_size = Array.length f_rr);
    let srm = Bigarray.Array2.create Bigarray.float64 Bigarray.c_layout a_size a_size in
    for r = 0 to (a_size-1) do
        assert(a_size = Array.length f_rr.(r));
        let diag = ref 0.0 in
        for c = 0 to (a_size-1) do
            if (c <> r) then 
                (diag := !diag +. f_rr.(r).(c); srm.{r,c} <- f_rr.(r).(c);)
        done;
        srm.{r,r} <- ~-. !diag;
    done;
    if phyml_mode then m_meanrate srm pi_;
    srm

let test_model sub_mat t = 
    let a_size = Bigarray.Array2.dim1 sub_mat in
    let (u_,d_,ui_) = 
        let n_d = Bigarray.Array2.create Bigarray.float64 Bigarray.c_layout a_size a_size
        and n_ui = Bigarray.Array2.create Bigarray.float64 Bigarray.c_layout a_size a_size in
        Bigarray.Array2.fill n_d 0.0;
        Bigarray.Array2.fill n_ui 0.0;
        diagonalize_gtr scratch_space sub_mat n_d n_ui;
        (sub_mat, n_d, n_ui)
    in
    print_barray2 (compose_gtr scratch_space u_ d_ ui_ t)

let create_lk_model (spec : Parser.SC.static_spec) : cm =
    let model = match spec.Parser.SC.st_type with
        | Parser.SC.STLikelihood m -> m 
        | _ -> failwith "not a likelihood model"
    in 
    let (a_size,a_gap) = 
        let alph = spec.Parser.SC.st_alph in
        match model.Parser.SC.use_gap with
        | true -> Alphabet.size alph, (-1)
        | false -> (Alphabet.size alph) - 1, Alphabet.get_gap alph
    in
    (* set up all the probability and rates *)
    let variation,probabilities,alpha,invar,sites =
        match model.Parser.SC.site_variation with
        | None ->
            Bigarray.Array1.of_array Bigarray.float64 Bigarray.c_layout ([| 1.0 |]),
            Bigarray.Array1.of_array Bigarray.float64 Bigarray.c_layout ([| 1.0 |]),
            None,None,1
        | Some a -> begin match a with
            | Parser.SC.Invariant ->  (* same as NONE *)
                 Bigarray.Array1.of_array Bigarray.float64 Bigarray.c_layout ([| 1.0 |]),
                 Bigarray.Array1.of_array Bigarray.float64 Bigarray.c_layout ([| 1.0 |]),
                 None,None,1
            | Parser.SC.Gamma (x,y,z) -> (* SITES,ALPHA,BETA *)
                let p = Bigarray.Array1.create Bigarray.float64 Bigarray.c_layout x in
                Bigarray.Array1.fill p (1.0 /. (float_of_int x));
                gamma_rates y z x,p,Some y,None,x
            | Parser.SC.Theta (w,x,y,z) -> (* GAMMA_SITES,ALPHA,BETA,PERCENT_INVAR *)
                if w < 1 then
                    failwith "Number of rate categories must be >= 1 (if invar w/out gamma)"
                else if w = 1 then begin
                    Bigarray.Array1.of_array Bigarray.float64 Bigarray.c_layout ([| 1.0 |]),
                    Bigarray.Array1.of_array Bigarray.float64 Bigarray.c_layout ([| 1.0 |]),
                    None,Some z,1
                end else begin
                    let p = Bigarray.Array1.create Bigarray.float64 Bigarray.c_layout w in
                    Bigarray.Array1.fill p (1.0 /. (float_of_int w));
                    let r = gamma_rates x y w in
                    r,p,Some x,Some z,w
                end
        end
    in
    (* check the rates so SUM(r_k*p_k) == 1 and SUM(p_k) == 1 |p| == |r| *)
    assert( (Bigarray.Array1.dim probabilities) = (Bigarray.Array1.dim variation) );
    let rsps = ref 0.0 and ps = ref 0.0 in
    for i = 0 to (Bigarray.Array1.dim probabilities) - 1 do
        rsps := !rsps +. ((Bigarray.Array1.get probabilities i) *. (Bigarray.Array1.get variation i));
        ps := !ps +. Bigarray.Array1.get probabilities i;
    done;
    assert( !rsps =. 1.0 && !ps =. 1.0 );
    (* extract the prior probability *)
    let priors = match model.Parser.SC.base_priors with
        | Parser.SC.Estimated p
        | Parser.SC.Given p -> assert(a_size = Array.length p); p
    in
    (*  get the substitution rate matrix and set sym variable and to_formatter vars *)
    let sym,mname,sub_mat = match model.Parser.SC.substitution with
        | Parser.SC.JC69 lambda -> 
                true,"JC69", m_jc69 priors lambda a_size
        | Parser.SC.K2P (alpha,beta) ->
                true, "K2P", m_k2p priors alpha beta a_size
        | Parser.SC.F81 lambda ->
                false, "F81", m_f81 priors lambda a_size
        | Parser.SC.F84 (kappa,beta) ->
                false, "F84", m_f84 priors kappa beta a_size
        | Parser.SC.HKY85 (alpha,beta) ->
                let model = 
                    if phyml_mode then m_hky85_ratio priors (alpha /. beta) a_size
                                  else m_hky85 priors alpha beta a_size
                in
                false, "HKY85", model
        | Parser.SC.TN93 (alpha1,alpha2,beta) ->
                false, "TN93", m_tn93 priors alpha1 alpha2 beta a_size
        | Parser.SC.GTR coeff -> 
                false, "GTR", m_gtr priors coeff a_size
        | Parser.SC.File matrix ->
                false, "File", m_file priors matrix a_size
    in
    (* diagonalize to get factored probability matrix --submat is destroyed
     * but can be reconstructed through the model parameters *)
    let (u_,d_,ui_) = 
        let n_d = Bigarray.Array2.create Bigarray.float64 Bigarray.c_layout a_size a_size in
        Bigarray.Array2.fill n_d 0.0;
        if sym then (
            diagonalize_sym scratch_space sub_mat n_d; (sub_mat, n_d, None) )
        else (
            let n_ui = Bigarray.Array2.create
                            Bigarray.float64 Bigarray.c_layout a_size a_size in
            Bigarray.Array2.fill n_ui 0.0;
            diagonalize_gtr scratch_space sub_mat n_d n_ui;
            (sub_mat, n_d, Some n_ui)
        )
    in
    { rate = variation;
      prob = probabilities;
      name = mname;
     alpha = alpha;
     invar = invar;
     sites = sites;
      pi_0 = Bigarray.Array1.of_array Bigarray.float64 Bigarray.c_layout priors;
         u = u_; d = d_;ui = ui_ }

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
(* required functions *)

(* [median] calculate the new new node between [an] and [bn] with
 * distance [t1] + [t2], being applied to[an], [bn] respectively    *)
let median an bn t1 t2 acode bcode =
    if pure_ocaml then begin
        let faa,loglike = ocaml_median an bn acode bcode t1 t2 in
        { an with
            chars = 
                bigarray_s 
                    (Bigarray.Array3.of_array Bigarray.float64
                                              Bigarray.c_layout
                                              faa)
                    None; (* TODO *)
            mle = loglike;
        }
    end else
        let am = an.model in
        let n_chars = match am.ui with
            | None -> 
                median_sym scratch_space am.u am.d t1 t2 an.chars bn.chars am.rate am.prob
            | Some ui -> 
                median_gtr scratch_space am.u am.d ui t1 t2 an.chars bn.chars am.rate am.prob
        in
        let pinvar = match an.model.invar with | Some x -> x | None -> ~-.1.0 in
        let loglike = loglikelihood n_chars an.model.pi_0 an.model.prob pinvar in
        { an with
            chars = n_chars;
            mle = loglike; 
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

(* Parser.SC.static_spec -> ((int list option * int) array) -> t *)
let of_parser spec characters =
    let computed_model = create_lk_model spec in
    let (a_size,a_gap) = 
        let model = match spec.Parser.SC.st_type with
            | Parser.SC.STLikelihood m -> m 
            | _ -> failwith "not a likelihood model" in 
        let alph = spec.Parser.SC.st_alph in
        match model.Parser.SC.use_gap with
        | true -> Alphabet.size alph, (-1)
        | false -> (Alphabet.size alph) - 1, Alphabet.get_gap alph
    in
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
                assert( a_size = List.length pl);
                Array.of_list (sublist pl 0 a_size)
    in
    (* convert character array to abstract type --redo *)
    let aa_chars = Array.map loop_ characters in (* create initial arrays *)
    let ba_chars = (* construct Array3 dim1=rates,dim2=chars,dim3=states *)
        Array.init (Bigarray.Array1.dim computed_model.rate)
                   (fun _ -> Array.copy aa_chars)
            --> Bigarray.Array3.of_array Bigarray.float64 Bigarray.c_layout
    and aa_chars =
        Array.map farray_to_int32 aa_chars -->
            Bigarray.Array1.of_array Bigarray.int32 Bigarray.c_layout
    in
    let lk_chars = match computed_model.invar with
                | Some _ -> bigarray_s ba_chars (Some aa_chars)
                | None   -> bigarray_s ba_chars None
    in
    let pinvar = match computed_model.invar with | Some x -> x | None -> ~-.1.0 in
    let loglike = loglikelihood lk_chars
                                computed_model.pi_0
                                computed_model.prob
                                pinvar
    in
    {    mle = loglike;
       model = computed_model;
       codes = Array.map (fun (x,y) -> y) characters; 
       chars = lk_chars; }

let to_formatter attr mine (t1,t2) data : Xml.xml Sexpr.t list =
    let str_time = function | Some x -> `Float x | None -> `String "None"
    and some_or x = function | Some x -> x | None -> x in

    let rec make_single_vec char_code single_ray = 
        (Array.to_list 
            (Array.mapi 
                (fun state_code value ->
                    let alph = Data.to_human_readable data char_code state_code in
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
    let priors = 
        [(PXML
            -[Xml.Characters.priors]
                { `Set (make_single_vec (Array.get mine.codes 0) 
                                        (ba2array mine.model.pi_0)) }
        --)]
    in
    let model =
        (PXML
            -[Xml.Characters.model]
                ([Xml.Characters.name] = [`String mine.model.name])
                ([Xml.Characters.sites] = [`Int mine.model.sites])
                ([Xml.Characters.alpha] = [str_time mine.model.alpha])
                ([Xml.Characters.invar] = [`Float (some_or 0.0 mine.model.invar)])
                { `Set priors }
        --)
    and sequence =
        let likelihood_vec,invariant_vec = s_bigarray mine.chars in
        (** TODO :: invar vector, and group rates in output(?) **)
        (PXML
            -[Xml.Characters.characters]
            {  
                let r = Array.to_list (barray_3matrix (likelihood_vec)).(0) in
                `Set (List.map2 (make_single) (Array.to_list mine.codes) r)
            }
        --)
    in

    (PXML
        (* tag *)
        -[Xml.Characters.likelihood]
            (* attributes *)
            ([Xml.Characters.llike] = [`Float mine.mle])
            ([Xml.Nodes.min_time] = [str_time t1])
            ([Xml.Nodes.oth_time] = [str_time t2])
            ([attr])
            (* data *)
            { `Set [model;sequence] }
        --) :: []
(* -> Xml.xml Sexpr.t list *)

(* readjust the branch lengths to create better mle score *)
let readjust xopt x c1 c2 _ c_t1 c_t2 =
    if pure_ocaml then begin
        let mine = median c1 c2 c_t1 c_t2 0 0 in
        let t1,t2,ll = ocaml_readjust c1 c2 c_t1 c_t2 mine.mle in
        let x = Array.fold_right 
                (fun c s -> All_sets.Integers.add c s) mine.codes x in
        (x,mine.mle,ll,(t1,t2), {mine with mle = ll; } )
    end else begin
        let model = c1.model
        and pinv  = match c1.model.invar with | Some x -> x | None -> ~-.1.0 in
        let new_mine = median c1 c2 c_t1 c_t2 0 0 in
        (*Printf.printf "S: %f\t%f\t%f\n%!" c_t1 c_t2 new_mine.mle;*)
        let (nta,nl) = match model.ui with
            | None ->
                readjust_sym scratch_space model.u model.d 
                             c1.chars c2.chars new_mine.chars c_t1 c_t2 pinv
                             model.rate model.prob model.pi_0 new_mine.mle
            | Some ui ->
                readjust_gtr scratch_space model.u model.d ui
                             c1.chars c2.chars new_mine.chars c_t1 c_t2 pinv
                             model.rate model.prob model.pi_0 new_mine.mle
        and ntb = c_t2
        in
        (* Printf.printf "E: %f\t%f\t%f\n%!" nta ntb nl; *)
        if nta =. c_t1 then
            (x,new_mine.mle,new_mine.mle,(c_t1,c_t2),new_mine)
        else
            let x = Array.fold_right (* bottle neck? *)
                    (fun c s -> All_sets.Integers.add c s) new_mine.codes x in
            (x,new_mine.mle,nl,(nta,ntb), {new_mine with mle = nl;} )
    end

let distance a_node b_node t1 t2 = (* codes don't matter here *)
    let t = median a_node b_node t1 t2 0 0 in t.mle

(* insert a node between two *)
(* (c, b) -> (c,(a),b)  *)
let dist_2 n a b nt at bt xt= 
    let x = median a b at bt 0 0 in
    let tt = 0.5 in (* estimate time here *)
    let y = median x n tt nt 0 0 in
    y.mle

let median_3 p x c1 c2  =  x
let reroot_median a b at bt = median a b at bt 0 0
let median_cost ta = ta.mle (*
    let pinvar = match ta.model.invar with | Some x -> x | None -> ~-.1.0 in
    loglikelihood ta.chars ta.model.pi_0 ta.model.prob pinvar *)
let root_cost t = t.mle
let to_string _ = "MLStaticCS"
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

ELSE
type t = unit
END
