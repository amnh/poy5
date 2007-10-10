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

(* $Id: paroperas.ml 1644 2007-02-14 19:05:47Z andres $ *)
let () = SadmanOutput.register "Paroperas" "$Revision: 1644 $"


(* The only difference between paroper and paroperas is the usage of Mpi.scatter
* to split the work among the slaves. Note that we will depend on the
* implementation of Mpi.scatter for performing this operation, as the data has
* to be to each of the processes running concurrently from the master. *)

exception Empty_List;;

let _comm_context = Mpi.comm_world
and master = 0;;

let inc a = a := !a + 1;;
let dec a = a := !a - 1;;

(* [_break_arr_for_slaves a] partition the array reference a in same size sets 
* among the processors available in the context communicator. The rest r of the 
* division is spread among the first r processes: these will have one more task to 
* perform. *)
let _break_arr_for_slaves a =
    let _comm_context_size = Mpi.comm_size _comm_context in
    let len = Array.length !a in
    let rest = ref (len mod _comm_context_size)
    and each =
        (try len / _comm_context_size
         with
         | e -> prerr_string "Divide by zero in paroperas.ml"; raise e)
    and counter = ref 0
    and result = Array.make _comm_context_size [||] in
    for i = 0 to _comm_context_size - 1 do
        if !rest > 0 then begin
            result.(i) <- Array.make (each + 1) (!a.(0));
            for j = !counter to !counter + each do
                let k = j - !counter in
                (result.(i)).(k) <- (!a.(j)); 
            done;
            counter := !counter + each + 1;
            dec rest;
        end else begin
            result.(i) <- Array.make each !a.(0);
            for j = !counter to !counter + each - 1 do
                let k  = j - !counter in
                (result.(i)).(k) <- !a.(j);                 
            done;
            counter := !counter + each;
        end;
    done;
    ref result;;

(* generic operation for the min and max functions. Sends the arrays referenced
* by a to each of the slave processors. *)
let _all_min_max_generic f a =
    let to_process = Mpi.scatter !a master _comm_context in
    let result = Array.map f to_process in
    (Array.length result), ref result;;



let all_minimize f c a =
    let a = _break_arr_for_slaves a in
    let len, result = _all_min_max_generic f a 
    and minima = ref Pervasives.max_float in
    for i = 0 to len - 1 do
        let cost = c (!result.(i)) in
        if  cost < !minima then minima := cost;
    done;
    Mpi.reduce_float !minima Mpi.Float_min master _comm_context;;

let all_maximize f c a =
    let a = _break_arr_for_slaves a in
    let len, result = _all_min_max_generic f a 
    and maxima = ref Pervasives.min_float in
    for i = 0 to len - 1 do
        let cost = c (!result).(i) in
        if cost > !maxima then maxima := cost;
    done;
    Mpi.reduce_float !maxima Mpi.Float_max master _comm_context;;

let gather_minimum f c a th =
    let a = _break_arr_for_slaves a in
    let len, arr = _all_min_max_generic f a 
    and minima = ref Pervasives.max_float in
    let arr_minima = Array.make len Pervasives.max_float in
    for i = 0 to len - 1 do
        arr_minima.(i) <- c (!arr).(i);
        if !minima > arr_minima.(i) then minima := arr_minima.(i);
    done;
    let minimum = Mpi.allreduce_float !minima Mpi.Float_min _comm_context
    and best_list = ref [] in
    for i = 0 to len - 1 do
        if arr_minima.(i) < minimum +. th then best_list := (!arr).(i) ::
            !best_list;
    done;
    let best_list = Array.of_list !best_list in
    let best_array = ref (Mpi.gather best_list master _comm_context) in
    Array_ops.flatten_array best_array;;

let gather_maximum f c a th =
    let a = _break_arr_for_slaves a in
    let len, arr = _all_min_max_generic f a 
    and maxima = ref Pervasives.min_float in
    let arr_maxima = Array.make len Pervasives.min_float in
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
    let best_list = Array.of_list !best_list in
    let best_array = ref (Mpi.gather best_list master _comm_context) in
    Array_ops.flatten_array best_array;;

let map f a =
    let a = _break_arr_for_slaves a in
    let len, result = _all_min_max_generic f a in
    let arr = ref (Mpi.gather !result 0 _comm_context) in
    Array_ops.flatten_array arr;;

let _all_scan f a =
    let res = Mpi.scatter !a master _comm_context 
    and lst = ref [] in
    for i = 0 to (Array.length res) - 1 do
        if (f res.(i)) then lst := res.(i) :: !lst;
    done;
    Array.of_list !lst;;

let scan f a =
    let a = _break_arr_for_slaves a in
    let result = _all_scan f a in
    let arr = ref (Mpi.gather result 0 _comm_context) in
    Array_ops.flatten_array arr;;

(* vim:sw=4 et tw=80
 *)
