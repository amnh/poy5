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

(** Building trees from initial data.
*
* This module provides several methods to build starting tree positions for
* phylogenetic analysis. *)

module type S = sig
      type a 
      type b

      val report_mst : Data.d -> a list -> string option -> unit

        (** [prebuilt a b] generates a list of trees using the list of trees [a] and the
    * node data [b]. It is required that the names of the taxa in [a] have data
    * associated in [b]. [b] is generated by the {!Node.load_data} function. *)
      val prebuilt :
          (string option * Tree.Parse.tree_types list) list ->
        Data.d * a list ->
        [> `Set of [> `Single of (a, b) Ptree.p_tree ] list ]

    val build_initial_trees :
      (a, b) Ptree.p_tree Sexpr.t -> Data.d ->
      a list -> Methods.build -> (a, b) Ptree.p_tree Sexpr.t

    end

module MakeNormal :
  functor
      (Node : NodeSig.S) -> 
          functor (Edge : Edge.EdgeSig with type n = Node.n) -> 
              functor
                  (TreeOps : Ptree.Tree_Operations with type a = Node.n with type b = Edge.e) ->
                      S with type a = Node.n with type b = Edge.e
module Make :
  functor
      (Node : NodeSig.S with type other_n = Node.Standard.n) -> 
          functor (Edge : Edge.EdgeSig with type n = Node.n) -> 
              functor
                  (TreeOps : Ptree.Tree_Operations with type a = Node.n with type b = Edge.e) ->
                      S with type a = Node.n with type b = Edge.e

val get_transformations : Methods.build -> Methods.transform list
