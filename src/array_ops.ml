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

(* $Id: array_ops.ml 1644 2007-02-14 19:05:47Z andres $ *)
let () = SadmanOutput.register "Array_ops" "$Revision: 1644 $"


exception Empty

let rec _calculate_size a len lst sum =
    match len with 
    | (-1) -> sum, lst
    | _ ->
            let s = Array.length (a.(len)) in
            _calculate_size a (len - 1) (s :: lst) (sum + s)

let rec _transfer_cont a tgt lst cura curtgt =
    match lst with
    | [] -> tgt
    | 0 :: t -> _transfer_cont a tgt t (cura + 1) curtgt
    | h :: t ->
            let new_curtgt = curtgt + h in
            for i = curtgt to new_curtgt - 1 do
                tgt.(i) <- a.(cura).(i - curtgt);
            done;
            _transfer_cont a tgt t (cura + 1) new_curtgt

(* Returns the first element in any of the non empty arrays in the array array
* a. *)
let rec _get_first a lst count =
    match lst with
    | 0 :: t -> _get_first a t (count + 1)
    | _ :: t -> a.(count).(0)
    | [] -> raise (Empty)

(* Note that tgt should have the same size as the sum of all the
 arrays pointed by a. The lst, cura and curtgt are variables for tail
 recursion. *)
let flatten_array a =
    let len = Array.length a in
    let size, lst = _calculate_size a (len - 1) [] 0 in
    match size with
    | 0 -> [||]
    | _ ->
            let sample = _get_first a lst 0 in
            let tmp = Array.make size sample in
            _transfer_cont a tmp lst 0 0

let rec _rec_fold_min f v a l it max =
    if max = it then (Array.of_list l)
    else begin
        let res, cost = f a.(it) in
        if cost > v then _rec_fold_min f v a l (it + 1) max
        else if cost = v then _rec_fold_min f v a (res :: l) (it + 1) max
        else _rec_fold_min f cost a [res] (it + 1) max
    end

let rec _rec_fold_max f v a l it max =
    if max = it then (Array.of_list l)
    else begin
        let res, cost = f a.(it) in
        if cost < v then _rec_fold_max f v a l (it + 1) max
        else if cost = v then _rec_fold_max f v a (res :: l) (it + 1) max
        else _rec_fold_max f cost a [res] (it + 1) max
    end

let int_f func c elem = 
    let r = func elem in
    r, (c r)

let fold_min (f : 'a -> 'b) (c : 'b -> int) a =
    _rec_fold_min (int_f f c) max_int a [] 0 (Array.length a)

let fold_min_p f c a min =
    _rec_fold_min (int_f f c) min a [] 0 (Array.length a)


let fold_max f c a =
    _rec_fold_max (int_f f c) min_int a [] 0 (Array.length a)

let fold_max_p f c a max =
    _rec_fold_max (int_f f c) max a [] 0 (Array.length a)

let map_ip f a =
    let len = Array.length a in
    if 0 < len then begin
        for i = 0 to len - 1 do
            a.(i) <- f (a.(i));
        done;
    end else ()

let mapi_ip f a =
    let len = Array.length a in
    if 0 < len then begin
        for i = 0 to len - 1 do
            a.(i) <- f i (a.(i));
        done;
    end else ()

let randomize ar =
    let l = Array.length ar - 1 in
    for i = 0 to l do
        let rnd = i + Random.int (l - i + 1) in
        let tmp = ar.(i) in
        ar.(i) <- ar.(rnd);
        ar.(rnd) <- tmp;
    done;;

let filter f arr =
    let res = Array.fold_right (fun x acc -> 
        if (f x) then x :: acc else acc) arr [] in
    Array.of_list res

(* vim: set sw=4 ts=4 et tw=80 : *)
