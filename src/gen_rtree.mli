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

(** gen_rtree.mli - this file defines an interface to a generalized rooted tree.
* The nodes of this type of tree can have any number of children. Nodes with no
* children could be additionally distinguished by a separate type. 
* Two consensus methods are implemented 1) a linear time randomized 
* algorithm to create a consensus tree for a given set of trees.
* Ref: A Linear Time Majority Algorithm for Majority Tree, Amenta,
* Clarke and St. John, LNCS ... Some of the terminology used here reflects
* the terminology used in that paper. 
* 2) A straight-forward implementation of the semi-strict consensus trees. The
* algorithm replicates the mathematical definition of semi-strict consensus
* using a O((n * t)^2) algorithm. *)

(* Indicates an invalid node id. *)
exception Invalid_Node_Id

(** The nodes of the consensus tree. The nodes are flexible enough to allow
* for polytomy. *)
type node = Root of int * int list (* (id, child list) *)
          | Interior of int * int * int list (* (id, parent id, child list) *)
          | Leaf of int * int (* (id, parent id) *)

(** The consensus tree/forest. The following representation is flexible enough
* to allow for the representation of both trees and forests. *)
type g_rtree = {
  g_rtopo : node All_sets.IntegerMap.t ; 
  (** the topology of the consensus tree stored as a map of tree nodes. *)

  roots : All_sets.Integers.t ;         
  (** the set of roots of the forest. *)
}

val empty : g_rtree
(** the empty consensus tree. *)

val get_id : node -> int
(** function to return the id of any node. *)

val get_node : int -> g_rtree -> node
(** function to return the node, given its id. *)

val get_parent : int -> g_rtree -> int
(** function to return the parent id of a node. *)

val is_node : int -> g_rtree -> bool
(** checks whether the passed in int corresponds to a node id. *)

val is_root : int -> g_rtree -> bool
(** checks whether the passed in id, belongs to a root node in the forest. *)

val is_leaf : int -> g_rtree -> bool
(** checks whether the passed in int belongs to a leaf node. *)

val pre_order_node_visitor :
  (int -> 'a -> 'a) ->
  int -> g_rtree -> 'a -> 'a
val post_order_node_visitor :
  (int -> 'a -> 'a) ->
  int -> g_rtree -> 'a -> 'a
(** functions to perform pre and post order traversals on the consensus tree. *)
  
val get_num_nodes : int -> g_rtree -> int * int
(** function to get the number of nodes in the consensus tree. The returned
* values are the (num_leaves, num_nodes) *)

val print_node : node -> unit
(** prints a node's type and the ids of its parent and its children. *)

val of_tree : Rtree.r_tree -> g_rtree
(** function to convert a rooted tree into a generalized rooted tree. *)

val build_ctree :
    int ->
    int ->
    (int * Rtree.r_tree) list ->
    All_sets.IntegerMap.key * g_rtree * (int, int) Hashtbl.t *
    (int, int) Hashtbl.t
(** builds a consensus tree given 1. maj 2. num_of_nodes and 
   3. (the root, rooted_tree) list. The return value is the consensus 
   tree root and consensus tree pair. *)

val get_cluster_lst :
    (int, int) Hashtbl.t ->
    int -> int * g_rtree -> (int * int * BitSet.t) list
(** Function to get the list of all clusters of the tree. *)

val build_semi_strict_ctree :
  int ->
  (int * g_rtree) list -> g_rtree
(** builds the semi-strict consensus tree of the passed in rooted generalized
* trees. *)

val to_parser :
    g_rtree -> (int -> string) -> string Parser.Tree.t
