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

(** rtree.mli - This file defines the interface for the rooted tree/forest. The
 data structures are very similar to the ones in the Tree module, with the
 exception that there are no explicitly directed edges. The edges are implicit. *)

(** Thrown when the node id is invalid. *)
exception Invalid_Node_Id

(** the node type of the rooted tree. *)
type node =
    Root of int * int * int (** (id, left_id, right_id) *)
  | Interior of int * int * int * int (** (id, parent id, left id, right id) *)
  | Leaf of int * int (** (id, parent id) *)

(** The rooteed tree data structure. *)
type r_tree = {
  
    r_topo : node All_sets.IntegerMap.t;
    (** The topology of the rooted tree represented as a map from node_id ->
      node. *)
  
    roots : All_sets.Integers.t;
    (** The set of roots of the rooted tree. *)
}

val empty : r_tree
(** The empty rooted tree. *)

val get_id : node -> int
(** returns the id of the corresponding node. *)

val get_node : int -> r_tree -> node
(** returns the node corresponding to the id. *)

val get_roots : r_tree -> All_sets.Integers.t
(** gets the set of roots of the tree. *)

val is_root : int -> r_tree -> bool
(** checks whether the passed in id belongs to a root. *)

val is_leaf : int -> r_tree -> bool
(** checks whether the id belongs to a leaf node or not. *)

val pre_order_node_visitor :
  (int -> 'a -> 'a) ->
  int -> r_tree -> 'a -> 'a
val post_order_node_visitor :
  (int -> 'a -> 'a) ->
  int -> r_tree -> 'a -> 'a
(** functions to visit the nodes in pre, post order. *)
  
val get_num_nodes : int -> r_tree -> int
(** gets the number of nodes in the tree. O(n) *)

val index_nodes :
  int * int ->
  int ->
  int ->
  r_tree -> (int, int) Hashtbl.t * (int, int) Hashtbl.t
(** indexes the nodes in a sequential order. The returned hash tables are the
* mappings from node_id (need not be sequential) to node_index (is sequential)
* and the mapping from node_index to node_id. *)

val of_tree :
  int ->
  Tree.break_jxn ->
  Tree.u_tree -> r_tree
(** roots the passed in unrooted tree. The passed in rt_id is used to add a root
* node that has as left and right child the start and end nodes of the break
* jxn. *)

val root_at : int -> int -> Tree.u_tree -> r_tree
(** function to root the given unrooted tree at the passed in outgroup using the
* given root_id. *)

val print_tree : All_sets.Integers.elt -> r_tree -> unit
val print_forest : r_tree -> unit
(** very simple functions to print the trees. *)
