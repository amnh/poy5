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
val deref : 'a option -> 'a
type sufNode_t = {
  mutable id : int;
  mutable father_ptr : sufNode_t ref option;
  mutable suf_link_ptr : sufNode_t ref option;
  mutable first_child_ptr : sufNode_t ref option;
  mutable right_sibling_ptr : sufNode_t ref option;
  mutable depth : int;
  mutable leaf_ls : int list;
  mutable start_pos : int;
  mutable end_pos_ptr : int ref;
}
val create_fresh_node : unit -> sufNode_t
val get_label_len : sufNode_t -> int
val set_leaf_label_pos : sufNode_t -> int -> int ref -> unit
val set_inter_label_pos : sufNode_t -> int -> int -> unit
val is_root : sufNode_t -> bool
val is_leaf : sufNode_t -> bool
val create_depth : sufNode_t -> sufNode_t -> unit
val cmp_num_child : sufNode_t -> int
val create_leaf_ls : sufNode_t -> int -> unit
val find_child :
  ?exact:bool -> sufNode_t -> int -> int array -> sufNode_t option
val get_left_sibling : sufNode_t -> sufNode_t -> sufNode_t
val travel : sufNode_t -> int array -> int -> int -> bool
val print_node : sufNode_t -> int array -> unit
val check_node :
  ?exact:bool -> sufNode_t -> int array -> int list -> bool * int list
val find_stop_node :
  ?exact:bool -> sufNode_t -> int array -> int list -> sufNode_t * int list
val find_leaves :
  ?exact:bool -> sufNode_t -> int array -> int list -> int list
