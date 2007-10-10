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

(** ctree.mli - this file defines an interface to a linear time randomized 
* algorithm to create a consensus tree for a given set of trees.
* Ref: A Linear Time Majority Algorithm for Majority Tree, Amenta,
* Clarke and St. John, LNCS ... Some of the terminology used here reflects
* the terminology used in that paper. *)

(** Indicates an invalid node id. *)
exception Invalid_Node_Id

(** The nodes of the consensus tree. The nodes are flexible enough to allow
 * for polytomy. *)
type node = Root of int * int list (** (id, child list) *)
          | Interior of int * int * int list (** (id, parent id, child list) *)
          | Leaf of int * int (** (id, parent id) *)

(** The consensus tree/forest. The following representation is flexible enough
 * to allow for the representation of both trees and forests. *)
type c_tree = {
  c_topo : node All_sets.IntegerMap.t ; 
  (** the topology of the consensus tree stored as a map of tree nodes. *)
  roots : All_sets.Integers.t ;         
  (** the set of roots of the forest. *)
}

val empty : c_tree
(** the empty consensus tree. *)

val get_id : node -> int
(** function to return the id of any node. *)

val get_node : int -> c_tree -> node
(** function to return the node, given its id. *)

val get_parent : int -> c_tree -> int
(** function to return the parent id of a node. *)

val is_node : int -> c_tree -> bool
(** checks whether the passed in int corresponds to a node id. *)

val is_root : int -> c_tree -> bool
(** checks whether the passed in id, belongs to a root node in the forest. *)

val is_leaf : int -> c_tree -> bool
(** checks whether the passed in int belongs to a leaf node. *)

val pre_order_node_visitor :
  (int -> 'a -> 'a) ->
  int -> c_tree -> 'a -> 'a
val post_order_node_visitor :
  (int -> 'a -> 'a) ->
  int -> c_tree -> 'a -> 'a
(** functions to perform pre and post order traversals on the consensus tree. *)
  
val get_num_nodes : int -> c_tree -> int
(** function to get the number of nodes in the consensus tree. *)

val print_node : node -> unit
(** prints a node's type and the ids of its parent and its children. *)

val build_ctree :
  int ->
  int ->
  (int * Rtree.r_tree) list ->
  int * c_tree
(** builds a consensus tree given 1. maj 2. num_of_nodes and 
    3. (the root, rooted_tree) list. The return value is the consensus 
    tree root and consensus tree pair. *)
