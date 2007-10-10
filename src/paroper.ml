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

(* $Id: paroper.ml 1644 2007-02-14 19:05:47Z andres $ *)
let () = SadmanOutput.register "Paroper" "$Revision: 1644 $"


exception Empty_List;;

let _comm_context = Mpi.comm_world;;

(* [generic_extract l] returns the starting point and length of the elements
* that should be processed in a hypothetical array fo size l by the calling
* process.  Converts the start, end format returned by cal_range to start, end *)
let generic_extract l =
    let s, l = Parparam.cal_range l in
    match l with
    | 0 -> s, s;
    | _ -> s, s + l - 1;;

(* [_all_min_max_generic f a] is a common starting operation for all the
* parallel processing. After dividing the array into subarrays for each process,
* applies the function f to each member in a that is responsability of the
* calling task. *)
let _all_min_max_generic f a =
    let s, e = generic_extract (Array.length !a) in
    let sample = f (!a).(s) 
    and sub_arr_size = e - s + 1 in
    let arr = Array.make sub_arr_size sample in
    for i = 1 to sub_arr_size - 1 do
        arr.(i) <- f (!a).(s + i);
    done;
    ref arr;;

let all_minimize f c a =
    let arr = _all_min_max_generic f a 
    and minima = ref Pervasives.max_float in
    for i = 0 to (Array.length !arr) - 1 do
        let cost = c (!arr).(i) in
        if !minima > cost then minima := cost;
    done;
    Mpi.allreduce_float !minima Mpi.Float_min _comm_context;;

let all_maximize f c a =
    let arr = _all_min_max_generic f a 
    and maxima = ref Pervasives.min_float in 
    for i = 0 to (Array.length !arr) - 1 do
        let cost = c (!arr).(i) in
        if !maxima < cost then maxima := cost;
    done;
    Mpi.allreduce_float !maxima Mpi.Float_max _comm_context;;

let gather_minimum f c a th =
    let arr = _all_min_max_generic f a 
    and minima = ref Pervasives.max_float in
    let arr_minima = Array.make (Array.length !arr) !minima in 
    for i = 0 to (Array.length !arr) - 1 do
        arr_minima.(i) <- c (!arr).(i);
        if !minima > arr_minima.(i) then minima := arr_minima.(i);
    done;
    let minimum = Mpi.allreduce_float !minima  Mpi.Float_min _comm_context 
    and best_list = ref [] in
    for i = 0 to (Array.length !arr) - 1 do
        if arr_minima.(i) < minimum +. th then best_list := (!arr).(i) ::
            !best_list;
    done;
    let arr = ref (Mpi.allgather (Array.of_list !best_list) _comm_context) in
    (Array_ops.flatten_array arr);;

let gather_maximum f c a th =
    let arr = _all_min_max_generic f a 
    and maxima = ref Pervasives.min_float in
    let len = Array.length !arr in
    let arr_maxima = Array.make len !maxima in
    for i = 0 to len - 1 do
        arr_maxima.(i) <- c (!arr).(i);
        if !maxima < arr_maxima.(i) then maxima := arr_maxima.(i);
    done;
    let maximum = Mpi.allreduce_float !maxima Mpi.Float_max _comm_context 
    and best_list = ref [] in
    for i = 0 to len - 1 do
        if arr_maxima.(i) > maximum -. th then best_list := (!arr).(i) ::
            !best_list;
    done;
    let arr = ref (Mpi.allgather (Array.of_list !best_list) _comm_context) in
    Array_ops.flatten_array arr;;

let map f a =
    let arr = _all_min_max_generic f a in
    let arr = ref (Mpi.allgather !arr _comm_context) in
    Array_ops.flatten_array arr;;

let scan f a =
    let s, e = generic_extract (Array.length !a) 
    and lst = ref [] in
    for i = s to e - 1 do
        if f !a.(i) then lst := !a.(i) :: !lst;
    done;
    let arr = ref (Array.of_list !lst) in
    Array_ops.flatten_array (ref (Mpi.allgather !arr _comm_context));;

(* vim:sw=4 et tw=80
 *)
