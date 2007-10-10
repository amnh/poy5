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

(* $Id: parparam.ml 1644 2007-02-14 19:05:47Z andres $ *)
let () = SadmanOutput.register "Parparam" "$Revision: 1644 $"


exception Impossible of string;;

let rec _cal_spread l n s p =
    if n < p then _cal_spread l (n + p) (n :: s) p
    else s;;

let cal_spread l =
    _cal_spread l (Mpi.comm_rank Mpi.comm_world) [] 
    (Mpi.comm_size Mpi.comm_world);;

let cal_range l = 
    let n_proc = Mpi.comm_size Mpi.comm_world in
    if l > n_proc then begin
        let size =
            try l / n_proc
            with
            | e -> prerr_string "Div by zero in parparam.ml"; raise e
        in
        (size * (Mpi.comm_rank Mpi.comm_world), size)
    end else
        raise (Impossible "The total number of processes is less than the
        number of jobs to be partitioned.");;

(* vim:sw=4 et tw=80
 *)
