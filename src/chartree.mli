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

(** [Chartree] works with trees of characters (things of type {!Node.node_data})
*)

(** This module works with trees with Node.node_data as the node data *)
type 'a p_tree = (Node.node_data, 'a) Ptree.p_tree

type 'a id = 
    | Same 
    | Cost of (Node.node_data, 'a) Ptree.p_tree 
    | All of (Node.node_data, 'a) Ptree.p_tree 

type tdp = MSame | MCost | MAll

type incremental = [
    | `Children of int
    | `No_Children of int
    | `HandleC of (int * int)
    | `HandleNC of (int * int)
]

(** {2 Printing} *)

val features : 
    Methods.local_optimum -> 
        (string * string) list -> (string * string) list

val print_nodes : 'a p_tree -> unit
(** [print_nodes tree] simply prints all the nodes in the tree (outputs a lot of
    information) *)
val print_node_data_indent :
  ?indent:int -> 'a p_tree -> unit
(** [print_node_data_indent ?(indent=3) tree] prints a simple, indented tree
    representation, including the {!Node.to_string} value of each node *)

(** {2 Downpass and Uppass} *)

val downpass_handle :
  int ->
  'a p_tree -> 'a p_tree
(** [downpass_handle handle tree] performs a downpass operation on the tree in
    [tree] with handle [handle].  This downpass uses existing nodes as previous
    alignments.  This operation is {i not} sufficient for having a correct cost
    associated with the tree. *)
val downpass :
  'a p_tree -> 'a p_tree
(** [downpass tree] performs a downpass on each handle of [tree] *)
val force_downpass_handle :
  int ->
  'a p_tree -> 'a p_tree
(** As above, but this downpass discards existing HTU nodes.  Useful for
    replacing the nodes in the OTU leaves. *)
val force_downpass :
  'a p_tree -> 'a p_tree
(** Performs {!Chartree.force_downpass} on all handles. *)
val uppass_handle :
  int ->
  'a p_tree -> 'a p_tree
(** [uppass_handle handle tree] performs an uppass on the tree in [tree] with
    handle [handle].  This establishes the correct cost for this [handle]. *)
val uppass :
  'a p_tree -> 'a p_tree
(** [uppass tree] performs an uppass for all handles in the tree.  This
    establishes the correct cost for the tree. *)

(** {2 Search interface} *)
type 'a break_fn_t = 
  Tree.id * Tree.id ->
  'a p_tree ->
  'a p_tree * Tree.break_delta * float * int * Node.node_data * 
  incremental list

(* Clears whatever node and root information should be for the given join_delta,
* and returns the vertex from which a downpass should be started to update the
* internal vertices. *)
val prepare_tree_for_downpass :
  ('a, 'b) Ptree.p_tree ->
  Tree.join_delta -> int * ('a, 'b) Ptree.p_tree

val break_fn : 'a break_fn_t
(** [break_fn] is a break function usable for Ptree's searches.  See {!Ptree}
    for more information. *)
val join_fn :
    incremental list ->
  Tree.join_jxn ->
  Tree.join_jxn ->
  'a p_tree ->
  'a p_tree * Tree.join_delta
(** [join_fn] is a break function usable for Ptree's searches.  See {!Ptree}
    for more information. *)
val cost_fn :
  Tree.join_jxn ->
  Tree.join_jxn ->
  float ->
  Node.node_data ->
    'a p_tree -> Ptree.clade_cost
(** [cost_fn] is a break function usable for Ptree's searches.  See {!Ptree}
    for more information. *)

val reroot_fn :
    bool -> 
  Tree.edge ->
  'a p_tree -> 'a p_tree * incremental list
(** [reroot_fn] is a break function usable for Ptree's searches.  See {!Ptree}
    for more information. *)


(** Module for performing these operations without skipping based on the
    min-cost-of-clade criterion (see above) *)
module TreeOps : 
    Ptree.Tree_Operations with type a = Node.node_data with type b =
        Node.node_data

val incremental_downpass : int -> Ptree.phylogeny -> Ptree.phylogeny *
incremental list

val incremental_uppass : Ptree.phylogeny -> incremental list -> Ptree.phylogeny

val aux_incremental_downpass :   
    tdp ->
    ([ `Children of int
    | `HandleC of int * int
    | `HandleNC of int * int
    | `No_Children of int ] as 'a) list ->
        int option ->
            (Node.node_data, 'b) Ptree.p_tree ->
                int -> (Node.node_data, 'b) Ptree.p_tree * 'a list

val incremental_downpass :   
    int ->
        (Node.node_data, 'a) Ptree.p_tree ->
            (Node.node_data, 'a) Ptree.p_tree *
            [> `Children of int
            | `HandleC of int * int
            | `HandleNC of int * int
            | `No_Children of int ] list

(** [subtree_to_formatter a b c d e f g] creates Tags.output of the contents of
* the subtree rooted by the node with code [e] in the phylogenetic tree [d],
* using as attributes of the root of the subtree [b], and as parent contents,
* and single assignment to the parent contents the tuple [g]. *)
val subtree_to_formatter : 
    ChromCS.IntSet.t * ChromCS.IntSet.t ->
    Tags.attributes -> Data.d -> Ptree.phylogeny -> int ->
    (Node.node_data * Node.node_data) option -> Tags.output

val handle_to_formatter : 
    All_sets.Integers.t * All_sets.Integers.t -> 
        Tags.attributes -> Data.d -> Ptree.phylogeny -> 
            int -> Tags.output

val to_formatter :
        Tags.attributes ->
        Data.d -> (Node.node_data, 'a) Ptree.p_tree -> Tags.output

val get_active_ref_code :(Node.node_data, 'a) Ptree.p_tree -> All_sets.Integers.t * All_sets.Integers.t
