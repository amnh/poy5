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

val sets : 
    Methods.tabu_join_strategy -> 
        Data.d ->
        ('a, 'b) Ptree.p_tree Sexpr.t -> All_sets.IntSet.t Lazy.t

val search_time_and_trees_considered : float -> int -> (string * string) list

val get_transformations : Methods.local_optimum -> Methods.transform list

module type S = 
    sig
        type a 
        type b

        val process_trees :
            Methods.information_contained list -> (a, b) Ptree.p_tree Sexpr.t
                -> (string option * Tree.Parse.tree_types) list *
                    (string * string, (int array * float option) list) Hashtbl.t

        val report_trees :
            Methods.information_contained list -> string option -> (a, b) Ptree.p_tree Sexpr.t -> unit

        val forest_break_search_tree :
            (a, b) Ptree.nodes_manager option -> float -> (a, b) Ptree.p_tree -> (a, b) Ptree.p_tree

        val diagnose :
            (a, b) Ptree.p_tree -> (a, b) Ptree.p_tree

        val find_local_optimum :
            ?base_sampler:(a, b) Sampler.search_manager_sampler  ->
          ?queue : (float array * (int * int) list * int * Status.status *
          int ref * float) ->
          Data.d ->
          Sampler.ft_queue ->
          (a, b) Ptree.p_tree Sexpr.t ->
          All_sets.IntSet.t Lazy.t ->
          Methods.local_optimum -> (a, b) Ptree.p_tree Sexpr.t

        val forest_search :
          Data.d ->
          Sampler.ft_queue ->
          float ->
          Methods.local_optimum ->
          (a, b) Ptree.p_tree
          Sexpr.t ->
          (a, b) Ptree.p_tree Sexpr.t

        val fusing :
          Data.d ->
          Sampler.ft_queue ->
          (a, b) Ptree.p_tree Sexpr.t ->
          int option * int option * Methods.tree_weights * 'a * Methods.local_optimum *
          (int * int) -> (a, b) Ptree.p_tree Sexpr.t

        val output_consensus : Data.d -> (a, b) Ptree.p_tree Sexpr.t -> string
        option -> float option -> bool -> unit
    end

module MakeNormal :
  functor (Node : NodeSig.S) ->
    functor (Edge : Edge.EdgeSig with type n = Node.n) ->
      functor
        (TreeOps : 
            Ptree.Tree_Operations with type a = Node.n with type b =
            Edge.e) -> S with type a = Node.n with type b = Edge.e

module Make :
  functor (NodeH : NodeSig.S with type other_n = Node.Standard.n) ->
    functor (EdgeH : Edge.EdgeSig with type n = NodeH.n) ->
      functor
        (TreeOpsH :
            Ptree.Tree_Operations with type a = NodeH.n with type b =
            EdgeH.e) -> 
                S with type a = NodeH.n with type b = EdgeH.e
