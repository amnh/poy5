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

let () = SadmanOutput.register "AlphSpec" "$Revision: 1644 $"

open StdLabels

exception Illegal_Element of string 

(* An element in the alphabet, and the length of the coded element in nats *)
type t = (string * float) list
type p = Equal | Enum of (string * float) list

let create ~elem:e =
    let size = List.length e in
    let prob = -. (log (1.0 /. (float_of_int size))) in
    List.map (fun x -> x, prob) e

let probs ~alph:a ~probs:p =
    match p with
    | Equal ->
            let size = List.length a in
            let prob = -. (log (1.0 /. (float_of_int size))) in
            List.map (fun (x, _) -> x, prob) a
    | Enum p ->
            let replace  lst (it, newp) =
                let rec do_replace = function
                    | (a, b) :: tl when a = it -> (a, -. (log newp)) :: tl
                    | hd :: tl -> hd :: (do_replace tl)
                    | [] -> raise (Illegal_Element it)
                in
                do_replace lst
            in
            List.fold_left ~f:replace ~init:a p

let rec length ~alph:a ~elem:e =
    match a with
    | (x, y) :: tl when x = e -> y
    | _ :: tl -> length ~alph:tl ~elem:e
    | [] -> raise (Illegal_Element e)

let to_list x = x
(* this is just the upper bound of the decoder complexity *)

let decoder x = 
    let e = exp 1.0 in
    let len = float_of_int (List.length x) in
    ((e *. len) -. 1.0) /. (e -. 1.0)

let to_formatter t = 
    let mapper (item, prob) =
        `Single 
            (Tags.KolSpecs.alph_element, 
            [(Tags.KolSpecs.value, item); (Tags.KolSpecs.prob, string_of_float
            prob)], `Structured `Empty)
    in
    `Set (List.map mapper t)
