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

type type_t = Alied | Deleted | Both
val fprintf : out_channel -> ('a, out_channel, unit) format -> 'a
type subseq_t = {
  mutable id : int;
  mutable sta : int;
  mutable en : int;
  mutable block_id_ls : int list;
}
val print : ?channel:out_channel -> subseq_t -> unit
val get_sum : subseq_t list -> int
val delete_block : subseq_t -> int -> unit
val is_equal : subseq_t -> subseq_t -> bool
val get_intersection : 'a -> 'b -> 'a -> 'b -> 'a * 'b
val cmp_del_cost : subseq_t -> Sequence.s -> Cost_matrix.Two_D.m -> int
val is_free : subseq_t -> bool
val get_subseq : Sequence.s -> subseq_t -> Sequence.s
val get_len : subseq_t -> int
