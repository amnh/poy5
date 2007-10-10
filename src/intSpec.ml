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

let () = SadmanOutput.register "IntSpec" "$Revision: 1644 $"

open StdLabels

exception Illegal_Element of int

type func = Constant of float array | Variable of float array

type t = { 
    min : int; 
    max : int; 
    codes: func;
}

type p = Equal | Function of (int -> float)

let create ~min:i ~max:a =
    let dif = (log (float_of_int (a - i))) in
    let f = Array.make (a - i - 1) dif in
    { min = i; max = a; codes = Constant f }

let codes ~ints:i ~codes:c =
    match c with
    | Equal -> 
            let dif = (log (float_of_int (i.max - i.min))) in
            let f = Array.make (i.max - i.min - 1) dif in
            { i with codes = Constant f }
    | Function f -> 
            let f = Array.init (i.max - i.min - 1) (fun x -> f x) in
            let f = 
                Array.map (fun x -> (-. (log x))) f
            in
            { i with codes = Variable f }

let length ~ints:i ~elem:x = 
    match i.codes with
    | Constant y
    | Variable y -> y.(x - i.min + 1)

let to_list x =
    let res = ref [] in
    for i = x.max downto x.min do
        let len = length x i in
        res := (i, len) :: !res;
    done;
    !res

let ( --> ) a b = b a

let decoder c = 
    match c.codes with
    | Constant _ ->
        let dif = float_of_int (c.max - c.min + 1) in
        log dif
    | Variable probs ->
            let rec sorted_insert item lst = 
                match lst with
                | f :: ((s :: t) as tl) ->
                        if f < item && item < s then f :: item :: tl
                        else f :: (sorted_insert item tl)
                | f :: [] -> [f; item]
                | [] -> [item]
            in
            let rec total_length len lst =
                match lst with
                | x :: y :: t -> 
                        let new_prob = x +. y in
                        let t = sorted_insert new_prob t in
                        total_length (len + 1) t
                | _ -> len
            in
            let len = total_length 0 (Array.to_list probs) in
            (float_of_int len) /. (log 2.)
            

let to_formatter t =
    `Single (Tags.KolSpecs.int_set, 
    [(Tags.KolSpecs.min, string_of_int t.min); 
    (Tags.KolSpecs.max, string_of_int t.max);
    (Tags.KolSpecs.prob, string_of_float (decoder t))], `Structured `Empty)

let bounds x = x.min, x.max

let get_class t = t.codes
