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

let () = SadmanOutput.register "Char_nonadd" "$Revision: 1644 $"

(* $Id: char_nonadd.ml 1644 2007-02-14 19:05:47Z andres $ *)

module IntSet = All_sets.Integers

type t = int * IntSet.t
type gen = int * int * bool

let rand_gen () =
    let justone = true in
    let outofmin = 50 in
    let outofmax = 100 in
    let code = Random.int 1073741823 in
    (code, (Random.int (outofmax - outofmin)) + outofmin, justone)

let make_rand (code, n, justone) =
    let ntomake =
        if justone
        then 1
        else (Random.int (n - 1)) + 1 in
    let myset = ref IntSet.empty in
    
    for i = 1 to ntomake do
        myset := IntSet.add (Random.int n) !myset
    done;
    (code, !myset)

let make_char code i =
    (code, IntSet.singleton i)


let color = Character.White
let code (c, _) = c

let median _ (ci, di) (cj, dj) =
    assert (ci = cj);
    let intersect =
        IntSet.inter di dj in
    if IntSet.is_empty intersect
    then
        (ci, IntSet.union di dj)
    else (ci, intersect)

(* let median_3 prev i j k = *)
(*     (\* This should not actually get called.  (We don't take 3-medians of *)
(*        non-additive characters.) *\) *)
(*     let intersect = *)
(*         (\* Parallelizable? *\) *)
(*         IntSet.inter (IntSet.inter i j) k in *)
(*     if IntSet.is_empty intersect *)
(*     then *)
(*         (\* Parallelizable? *\) *)
(*         IntSet.union (IntSet.union i j) k *)
(*     else intersect *)
        
(* Simplistic set distance through counting (undefined for two empty sets) *)
let distance (c1, i) (c2, j) =
    assert (c1 = c2);
    let intersect =
        (* Parallelizable? *)
        IntSet.inter i j in
    if IntSet.is_empty intersect
    then 1.
    else 0.

let (|>) a b =
    if a = 0
    then b
    else a

let compare_codes (c1, _) (c2, _) = compare c1 c2

let compare_data (c1, s1) (c2, s2) =
    (compare c1 c2) |> (IntSet.compare s1 s2)

let parse i = []
let to_string (c, i) =
    "{" ^ string_of_int c ^ ": " ^
        (IntSet.fold
             (fun i s -> (string_of_int i) ^
                  if s = ""
                  then "}"
                  else " " ^ s)
             i "")

(* let update old n1 n2 n3 = *)
(*     let median = median_3 old n1 n2 n3 in *)
(*     match old with *)
(*     | None -> (median, Character.rContinue) *)
(*     | Some old -> *)
(*           if *)
(*               (\* Parallelizable? *\) *)
(*               0 = compare_data old median *)
(*           then (median, Character.rStop) *)
(*           else (median, Character.rContinue) *)
