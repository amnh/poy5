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

(** metric.ml:  *)
let () = SadmanOutput.register "Metric" "$Revision: 1644 $"


module type Norm = sig
    type elt
    type space
    val sub : space -> elt -> elt -> elt
    val norm : space -> elt -> float
end

module type Metric = sig
    type elt
    type space
    val dist : space -> elt -> elt -> float
end

module MetricOfNorm(N : Norm) = struct
    type elt = N.elt
    type space = N.space
    let dist s a b =
        let c = N.sub s a b in
        N.norm s c
end

module DiscreteMetric = struct
    type elt = int
    type space = int * float array array
    let dist (n, d) a b =
        assert (0 <= a); assert (a < n);
        assert (0 <= b); assert (b < n);
        d.(a).(b)

    let make_space n ar =
        if n = Array.length ar
        then (n, ar)
        else invalid_arg "make_space"

    let make_elts n =
        let rec make_upto i acc =
            if i < 0 then acc
            else make_upto (i - 1) (i :: acc)
        in make_upto (n - 1) []
end

module LpdNorm = struct
    type space = float * int
    type elt = float array

    let add (_, d) a b =
        Array.init d (fun i -> a.(i) +. b.(i))
    let sub (_, d) a b =
        Array.init d (fun i -> a.(i) -. b.(i))
    let zero (_, d) = Array.make d 0.

    let norm (p, d) a =
        let rec sum i acc =
            if i = d
            then acc ** (1. /. p)
            else sum (succ i) (acc +. (a.(i) ** p))
        in sum 0 0.

    let make_space p d = (p, d)
    let make_elt d q = Array.init d q
end
module LpdMetric = MetricOfNorm(LpdNorm)

module MetricUtil(M : Metric) = struct
    let min_dist space a bs =
        List.fold_right (fun b prev ->
                             min (M.dist space a b) prev) bs infinity
    
    let embed_Bourgain space elts c rand =
(*         let log2 a = (log a) /. (log 2.) in *)

        let n = float_of_int (List.length elts) in
        let m = int_of_float (10. *. (log n)) in (* i = 1..m *)
        let q = int_of_float (c *. (log n)) in (* j = 1..q *)

        let fromindex i = (i / q, i mod q) in
        let rand j =
            let rec r j =
                if j = 0
                then true
                else if rand ()
                then false
                else r (j - 1)
(*             (\* let outof = int_of_float (2. ** (float_of_int j)) in *\) *)
(*             if  *)
(*             0 = rand outof in *)
            in r j
        in

        let a = Array.init m
            (fun i -> 
             let i = i + 1 in
             ignore(i + i);
             Array.init q
                 (fun j ->
                      let j = j + 1 in
                      List.filter (fun _ -> rand j) elts)) in

        let new_elts = List.map
            (fun x ->
                 LpdNorm.make_elt (m * q)
                     (fun i ->
                          let i, j = fromindex i in
                          min_dist space x a.(i).(j))) elts in
        (LpdNorm.make_space 1. (m * q), new_elts)
end

module DiscreteUtil = MetricUtil(DiscreteMetric)

type 'a distfn = 'a -> 'a -> float
type 'a eltset = 'a list * 'a distfn

let distortion (elts1, dist1) (elts2, dist2) =
    assert (List.length elts1 = List.length elts2);
    let minc, maxc =
        List.fold_right2
            (fun a1 b1 (mincoeff, maxcoeff) ->
                 List.fold_right2
                     (fun a2 b2 (mincoeff, maxcoeff) ->
                          if (a1 = a2) && (b1 = b2)
                          then (mincoeff, maxcoeff)
                          else let da = dist1 a1 a2 in
                          let db = dist2 b1 b2 in
                          let coeff = da /. db in
                          print_endline (string_of_float da ^ " "
                                         ^ string_of_float db ^ " "
                                         ^ string_of_float coeff);
                          (min mincoeff coeff,
                           max maxcoeff coeff))
                     elts1 elts2 (mincoeff, maxcoeff))
            elts1 elts2 (infinity, neg_infinity)
    in
(*     print_endline (string_of_float maxc ^ " " ^ string_of_float minc *)
(*                    ^ " " ^ string_of_float (maxc /. minc)); *)
    maxc /. minc

let convert_lpd_Bourgain a c =
    let n = (Array.length a) in
    let space = DiscreteMetric.make_space n a in
    let elts = DiscreteMetric.make_elts n in
    let bspace, belts = DiscreteUtil.embed_Bourgain space elts c Random.bool in
    (bspace, belts)

let compare_Bourgain a c =
    let n = Array.length a in
    let dspace, delts = (DiscreteMetric.make_space n a,
                         DiscreteMetric.make_elts n) in
    let bspace, belts = convert_lpd_Bourgain a c in
    let ddist = DiscreteMetric.dist dspace in
    let bdist = LpdMetric.dist bspace in
    distortion (delts, ddist) (belts, bdist)

let make_rand_a n max =
    Array.init n
        (fun i -> Array.init n
             (fun j ->
                  if i = j
                  then 0.
                  else Random.float max))

let make_rand_metric_a n max =
    let max = max ** (1. /. float_of_int n) in
    let space = LpdNorm.make_space 2. n in
    let elts = Array.init n
        (fun i -> LpdNorm.make_elt n
             (fun _ -> Random.float max)) in
    Array.init n
        (fun i -> Array.init n
             (fun j ->
                  LpdMetric.dist space elts.(i) elts.(j)))

let pairwise_dist (elts, dist) =
    let n = List.length elts in
    Array.init n
        (fun i -> Array.init n
             (fun j ->
                  dist (List.nth elts i) (List.nth elts j)))

let make_discrete a =
    let n = Array.length a in
    let space = DiscreteMetric.make_space n a in
    let elts = DiscreteMetric.make_elts n in
    (elts, DiscreteMetric.dist space)

let make_lpd (space, elts) = (elts, LpdMetric.dist space)

let avg_n_trials f n arg =
    let acc = ref 0. in
    for i = 0 to n - 1 do
        acc := !acc +. f arg
    done;
    !acc /. (float_of_int n)


let verify_triangle matr =
    let check a b c =
        if (a +. b) < c
        then failwith "No triangle ineq.!"
    in
    let n = Array.length matr in
    for i = 0 to n - 1 do
        print_endline ("Checking elt. " ^ string_of_int i ^ "..");
        let matri = matr.(i) in
        for j = 0 to n - 1 do
            print_string ("..with elt. " ^ string_of_int j);
            let d1 = matri.(j) in
            flush stdout;
            for k = 0 to n - 1 do
                let d2 = matr.(j).(k) in
                let d3 = matri.(k) in
                check d1 d2 d3;
                check d1 d3 d2;
                check d2 d3 d1;
                print_string "."
            done;
            print_newline ();
        done;
    done;
    ()

(* let a = make_rand_metric_a 500 50000. *)

let make_float_list from fto step =
    let rec f x acc =
        if x < from
        then acc
        else f (x -. step) (x :: acc)
    in f fto []

let input = (make_float_list 0.4 1.0 0.025)
let input = [0.95]


(* Convert Sankoff characters to additive characters *)
let get_sankoff_cm r =
    let sank = r.Node.preliminary in
    Array.map (fun x -> Array.map float_of_int x) (SankCS.tcm sank)

let sankoff_convert_fn ?(c=0.95) tcm =
    let a = Array.map (fun x -> Array.map float_of_int x) tcm in
    let bspace, belts = convert_lpd_Bourgain a c in (* TODO: regen if bad *)

(*     let ndims = Array.length (List.hd belts) in *)
(*     let addcs_init_array = Array.init ndims (fun i -> i) in *)
    let elt_to_add elt =
        let minstates = SankCS.get_minstates elt in
        let belt = List.fold_right
            (fun state elt ->
                 LpdNorm.add
                     bspace (List.nth belts state) elt) minstates
            (LpdNorm.zero bspace)
        in
        
        (* now take the element and make each dimension an additive character
           ... *)
        let scale = 1. in
        let belt_array = Array.mapi
            (fun i d -> (int_of_float (scale *. d),
                         int_of_float (scale *. d),
                         i)) belt in
        let add = AddCS.of_array belt_array (SankCS.ecode elt) in
        add
    in

    fun cs acc ->
        let sank = cs.Node.preliminary in
        let elts = SankCS.get_elt_array sank in
        let elts = Array.to_list elts in
        let additives = List.map elt_to_add elts in
        let additives = List.map (fun a -> { Node.preliminary = a;
                                             Node.final = a;
                                             Node.cost = 0.;
                                             Node.sum_cost = 0.;
                                             Node.weight = 1.}) additives in
        let additives = List.map (fun a -> Node.Add a) additives in
        additives @ acc


let sankoff_null_convert cs acc =
    failwith "unknown sankoff character code in conversion"
let remember_sankoff_convert ?(c=0.95) sank fn =
    let thisfn = sankoff_convert_fn ~c:c (SankCS.tcm sank) in
    let thiscode = SankCS.code sank in
    fun cs acc ->
        let sank = cs.Node.preliminary in
        let code = SankCS.code sank in
        if code = thiscode
        then thisfn cs acc
        else fn cs acc


(* let distortions = List.map *)
(*     (fun x -> (x, *)
(*                avg_n_trials (compare_Bourgain a) 5 x)) *)
(*     input *)

(*
let _ = List.iter (fun (x,y) -> print_float x; print_string " ";
(*                        print_float y; print_newline ()) *)
(*     distortions *)
*)
