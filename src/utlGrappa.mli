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

(** utlGrappa module provides functions for rearrangement 
* operations such as computing inversion, breakpoint distances 
* between two gene orders arrays *)

val fprintf : out_channel -> ('a, out_channel, unit) format -> 'a


(** [standardize genomeX genomeY] standardizes  
 * gene order [genomeX=(x1, x2,..., xk)] and  one of its 
 * permutations [genomeY=(y1, y2,...,yk)] 
 * into [sta_X] and [sta_Y] such that [genomeX=(1,...,k)]) 
 * For example: [genomeX] = (5, -3, 2) and [genomeY] = (-2, 3, 5)
 * [sta_X] = (1, 2, 3) and [sta_Y] = (-3, -2, 1 *)
val standardize : int array -> int array -> int array * int array

(** [cmp_inversion_dis genomeX genomeY circular] computes
 * the inversion distance between two given gene orders 
 * [genomeX=(x1, x2, ... xk)] and [genomeY=(y1, y2,..yk)]. 
 * Note [genomeY] is a permutation of [genomeX].
 * Compute the inversion distance between [genomeX] and [genomeY] using GRAPPA functions
 * For example: [genomeX] = (-6, 1, 5), [genomeY] = (-5, 1, 6) *)
val cmp_inversion_dis : int array -> int array -> int -> int

(** [cmp_breakpoint_dis genomeX genomeY circular] computes
 * the breakpoint distance between two given gene orders 
 * [genomeX=(x1, x2, ... xk)] and [genomeY=(y1, y2,..yk)]. 
 * Note that orientations are ignored. *)
val cmp_breakpoint_dis : int array -> int array -> int -> int

(** [cmp_oriented_breakpoint_dis genomeX genomeY circular] computes
 * the breakpoint distance between two given gene orders 
 * [genomeX=(x1, x2, ... xk)] and [genomeY=(y1, y2,..yk)]. 
 * Note that orientations are taken into account. 
 * For example: [genomeX] = (-6, 1, 5), [genomeY] = (-5, 1, 6) *)
val cmp_oriented_breakpoint_dis : int array -> int array -> int -> int

(** [get_ordered_permutation genomeX] returns 
 * the ordered permutation [genomeY]=(y1,..,yk | y1 < ... < yk) 
 * of [genomeX] =(x1, x2, ... xk), 
 * For example. [genomeX] = (-6, 1, 5), [genomeY] = (1, 5, 6) *)
val get_ordered_permutation : int array -> int array

(** [cmp_self_inversion_dis genome circular] 
 * computes the inversion distance between a gene orders [genome=(x1, x2, ... xk)]
 * and [Y=(y1,..,yk | y1 < ... < yk) using GRAPPA functions
 * where [Y] is an ordered permutation of |genome|
 * For example. [genome] = (-6, 1, 5), Y = (1, 5, 6) *)
val cmp_self_inversion_dis : int array -> int -> int

(** [cmp_self_breakpoint_dis genome circular] 
 * computes the breakpoint distance between a gene orders [genome=(x1, x2, ... xk)]
 * and [Y=(y1,..,yk | y1 < ... < yk) using GRAPPA functions
 * where [Y] is an ordered permutation of |genome|
 * For example. [genome] = (-6, 1, 5), Y = (1, 5, 6) *)
val cmp_self_breakpoint_dis : int array -> int -> int

(** [cmp_self_oriented_breakpoint_dis genome circular] 
 * computes the breakpoint distance between a gene orders [genome=(x1, x2, ... xk)]
 * and [Y=(y1,..,yk | y1 < ... < yk) using GRAPPA functions
 * where [Y] is an ordered permutation of |genome|
 * For example. [genome] = (-6, 1, 5), Y = (1, 5, 6) *)
val cmp_self_oriented_breakpoint_dis : int array -> int -> int
