(* POY 5.1.1. A phylogenetic analysis program using Dynamic Homologies.       *)
(* Copyright (C) 2011  Andrés Varón, Lin Hong, Nicholas Lucaroni, Ward Wheeler*)
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

(** [get_matcharr_and_costmatrix] is the main function of this module. it take
* two sequences, set of parameters (min lcb ratio, min lcb length, etc) and a
* costmatrix between charactors as input, output two code array of matching
* blocks, a cost_mat between lcb blocks and non-lcb blocks, ali_mat is the cost
* and alignment sequence between lcb blocks.
* 
* NOTE: edit_cost is the total editing cost between lcb blocks (which is not
* included is cost_mat), full_cost_listlst is a little bit redundant here, for
* it contains the two code array, but it also has the begin and end point of
* each block. *)
val get_matcharr_and_costmatrix :
    Sequence.s -> Sequence.s -> float -> float -> float -> float
        -> int*int -> Cost_matrix.Two_D.m ->  bool -> 
            int array * int array * int array array * 
            (int*Sequence.s*Sequence.s) array array * 
            int * int * int * (int * (int * int) * int list) list list 

val output2mauvefile :
    string -> int -> (int option) -> int array -> int array
        -> (int * (int * int) * int list) list list
            -> (int*Sequence.s*Sequence.s) array array ->int -> int -> int ->
    unit 
