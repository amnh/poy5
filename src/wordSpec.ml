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

let () = SadmanOutput.register "WordSpec" "$Revision: 1644 $"

open StdLabels

exception Illegal_Element of int
type t = { 
    alph:AlphSpec.t;
    minlen:int;
    maxlen:int;
    codes_spec:int -> float 
}

type p = Equal | Function of (int -> float)

let codes ~word:w ~codes:c =
    match c with
    | Equal -> 
            let dif = (log (float_of_int (w.maxlen - w.minlen))) in
            { w with codes_spec = fun _ -> dif }
    | Function f -> { w with codes_spec = f }


let create ~min:i ~max:a ~alph:c = 
    let init = 
        { alph = c; minlen = i; maxlen = a; codes_spec = (fun _ -> 1.0) }
    in
    codes init Equal

let length ~word:w ~length:l = w.codes_spec l

let wordlength ~word:w ~elem:e = 
    let len = String.length e in
    let res = ref (length w len) in
    for i = 0 to len - 1 do
        let tmp = String.sub e ~pos:i ~len:1 in
        let tmp = AlphSpec.length ~alph:w.alph ~elem:tmp in
        res := !res +. tmp;
    done;
    !res

let to_list w =
    let res = ref [] in
    for i = w.maxlen downto w.minlen do
        res := (i, w.codes_spec i) :: !res;
    done;
    !res

let to_formatter t =
    `Single 
    (Tags.KolSpecs.word_set, 
    [(Tags.KolSpecs.min, string_of_int t.minlen);
    (Tags.KolSpecs.max, string_of_int t.maxlen)], 
    `Structured (AlphSpec.to_formatter t.alph))
