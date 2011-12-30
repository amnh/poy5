(* POY 5.0 Alpha. A phylogenetic analysis program using Dynamic Homologies.   *)
(* Copyright (C) 2011 Andrés Varón, Lin Hong, Nicholas Lucaroni, Ward Wheeler,*)
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

(* $Id: supports.mli 2871 2008-05-23 17:48:34Z andres $ *)
(* Created Thu Feb  2 16:04:01 2006 (Illya Bomash) *)

(** This module implements computing a support diagnosis of a tree. *)
module type S = sig
    type a 
    type b

    val support : 
        (a, b) Ptree.p_tree Sexpr.t ->
            a list -> 
        Methods.support_method -> Data.d -> 
         Sampler.ft_queue ->
             int Tree.CladeFPMap.t

    val bremer_support : 
        (a, b) Ptree.p_tree Sexpr.t ->
        int -> int ->
        a list -> (a, b) Ptree.p_tree Sexpr.t -> Methods.local_optimum
        -> Methods.build -> Data.d -> Sampler.ft_queue ->
            Methods.support_tree Sexpr.t

val support_to_string_tree : Data.d -> Methods.support_tree ->
                             Tree.Parse.tree_types

(** [join_support_trees trees] takes a list of [(iterations, support_tree)]
    pairs and combines them into a single support tree *)
val join_support_trees : (int * Methods.support_tree) list -> Methods.support_tree


(** [bremer_of_input_file cg r ts d f t] returns an sexpr which is a map of the
* sexpr of the phylogenetic trees [t], where the input file [f] holds a list of
* phylogenetic trees for the loaded data [d], where each tree holds some annotated
* information in square brackets. The root of the trees is [r], and the, and the
* taxon name generation function (for each node code) is [ts]. The resulting
* sexpr is a set of printable trees with their associated bremer support values,
* as inferred from the input trees in [f]. *) 
val bremer_of_input_file :
    (Tree.u_tree -> string -> int) -> int ->
        (int -> string) -> Data.d -> Methods.filename list -> 
            (a, b) Ptree.p_tree Sexpr.t -> Tree.Parse.tree_types Sexpr.t

(** Like [bremer_of_input_file] but trust whatever input cost is provided with
* each tree in it's annotated information.*)
val bremer_of_input_file_but_trust_input_cost : int ->
    (int -> string) -> Data.d -> Methods.filename list -> 
        (a, b) Ptree.p_tree Sexpr.t -> Tree.Parse.tree_types Sexpr.t

end

module Make (Node : NodeSig.S with type other_n = Node.Standard.n) (Edge : Edge.EdgeSig with type n = Node.n) 
    (TreeOps : 
        Ptree.Tree_Operations 
        with type a = Node.n with type b = Edge.e) : S 
    with type a = Node.n
    with type b = Edge.e
