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

val deref : 'a option -> 'a
type rule = Rule1 | Rule21 | Rule22 | Rule31 | Rule32
val max_located_node : int
val cur_phase_ptr : int ref
type sufTree_t = {
  mutable code_arr : int array;
  mutable len : int;
  mutable root : SufNode.sufNode_t;
  mutable num_located_node : int;
  mutable num_used_node : int;
  mutable located_node_arr : SufNode.sufNode_t ref array;
  mutable gre_node_id : int;
}
val create_fresh_suf_tree : unit -> sufTree_t
val check_tree : sufTree_t -> unit
val do_first_phase : sufTree_t -> unit
val locate_mem : sufTree_t -> unit
val do_rule22 : sufTree_t -> int -> SufNode.sufNode_t -> int -> unit
val do_rule21 :
  sufTree_t ->
  int ->
  SufNode.sufNode_t ->
  int ->
  rule * SufNode.sufNode_t option * SufNode.sufNode_t option *
  SufNode.sufNode_t
val do_explicit_extension :
  sufTree_t ->
  int ->
  SufNode.sufNode_t ->
  int ->
  rule * SufNode.sufNode_t option * SufNode.sufNode_t option *
  SufNode.sufNode_t
val travel : sufTree_t -> SufNode.sufNode_t -> int -> SufNode.sufNode_t * int
val do_extension_rec :
  sufTree_t ->
  int ->
  SufNode.sufNode_t ->
  int -> SufNode.sufNode_t option -> SufNode.sufNode_t * int * int
val create_suf_tree : int array -> sufTree_t
val create_k_mer_tree : int array -> int -> sufTree_t
val find_stop_node : sufTree_t -> int list -> SufNode.sufNode_t * int list
