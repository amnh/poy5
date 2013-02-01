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

module type S = sig

    (** {6 Types *)

    (** abstract type for node in functorization *)
    type a 

    (** abstract type for edge in functorization *)
    type b

    (** tree composed of a and b *)
    type tree = (a,b) Ptree.p_tree

    (** wrapped tree; this is so we don't call the site_likelihood command for
        higher-efficency. *)
    type wtree = { t : tree; slk : float array; root : MlStaticCS.t; }

    (** Define a type to determine the way to optimize the tree to obtain
        likelihood scores for replicates of the tree. *)
    type replicate = { m:bool; b:bool }

    (** The RELL - relative estimated log-likelihood - no optimizations *)
    val rell : replicate

    (** Full optimizations for each replicate --very long time to optimize *)
    val full : replicate


    (** {6 Helper Functions *)

    (** [generate_wrapped_tree] create the wtree type from a tree. *)
    val create_wrapped_tree : Data.bool_characters -> tree -> wtree

    (** [can_perform_stat_tests t] test if we are using mlstatic data and that
        the tree is not disjoint. *)
    val can_perform_stat_tests : Data.bool_characters -> tree -> bool


    (** {6 Replicate Procedures *)

    (** [bootstrap_rep] generate a bootstrap replicate and return the weights
    * for each of the characters defined in n or the size of the cdf. *)
    val bootstrap_weights : ?n:int -> float array -> float array
    
    (** [replicate_cost] generate the cost of a replicate *)
    val replicate_cost : replicate -> wtree -> float array -> float

    (** [return a cdf from a tree that represent weights *)
    val get_cdf : wtree -> float array

    
    (** {6 Test statistics on trees *)

    (** Return the two-tailed P-Value for the KH test. It should not be used
        against the ML tree, the intention is a priori trees. *)
    val kh : ?n:int -> ?p:float -> ?rep:replicate -> ?chars:Data.bool_characters -> tree -> tree -> unit
 
    (** Return the P-Values and Trees for a candidate set of trees *) 
    val sh : ?n:int -> ?p:float -> ?rep:replicate -> ?chars:Data.bool_characters -> tree list -> unit

    (** Return the P-Value to support the best tree passed *)
    val au : ?n:int -> ?rep:replicate -> ?k:int -> ?chars:Data.bool_characters -> tree list -> unit

end

module Make :
    functor (NodeF : NodeSig.S with type other_n = Node.Standard.n) ->
        functor (Edge : Edge.EdgeSig with type n = NodeF.n) ->
            functor (TreeOps : Ptree.Tree_Operations with type a = NodeF.n with type b = Edge.e) ->
                S with type a = NodeF.n with type b = Edge.e

