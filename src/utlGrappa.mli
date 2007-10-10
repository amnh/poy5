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

val fprintf : out_channel -> ('a, out_channel, unit) format -> 'a
val write_2genome_file : string -> int array -> int array -> unit
val standardize : int array -> int array -> int array * int array
val cmp_inversion_dis : int array -> int array -> int -> int
val cmp_breakpoint_dis : int array -> int array -> int -> int
val cmp_oriented_breakpoint_dis : int array -> int array -> int -> int
val get_ordered_permutation : int array -> int array
val cmp_self_inversion_dis : int array -> int -> int
val cmp_self_breakpoint_dis : int array -> int -> int
val cmp_self_oriented_breakpoint_dis : int array -> int -> int
