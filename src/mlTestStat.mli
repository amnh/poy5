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

    (** Define a type to determine the way to optimize the tree to obtain
        likelihood scores for replicates of the tree. *)
    type replicate = { m:bool; b:bool }

    (** The RELL - relative estimated log-likelihood - no optimizations *)
    val rell : replicate

    (** Full optimizations for each replicate --very long time to optimize *)
    val full : replicate


    (** {6 Helper Functions *)

    (** [can_perform_stat_tests t] test if we are using mlstatic data and that
        the tree is not disjoint *)
    val can_perform_stat_tests : tree -> bool

    (** [analyze_tree] return RELL mean/variance and Site mean/variance *)
    val analyze_tree : tree -> (float * float) * (float * float)

    (** [analyze_pair] return the sitewise likelihood ratio variance and mean *)
    val analyze_pair : tree -> tree -> float * float


    (** {6 Replicate Procedures *)

    (** [bootstrap_data] generate a site-wise set of weights for characters
        based on a generated CDF, see get_cdf *)
    val bootstrap_data : ?n:int -> replicate -> float array -> float array
        
    (** [return a cdf from an array of floats that represent weights *)
    val get_cdf : float array -> float array

    
    (** {6 Test statistics on trees *)

    (** Return the two-tailed P-Value for the KH test. It should not be used
        against the ML tree, the intention is a priori trees. *)
    val kh : ?n:int -> ?p_star:float -> ?rep:replicate -> tree -> tree -> ()
 
    (** Return the P-Values and Trees for a candidate set of trees *) 
    val sh : ?n:int -> ?p_star:float -> ?rep:replicate -> tree list -> ()

    (** Return the P-Value to support the best tree passed *)
    val au : ?n:int -> ?p_star:float -> ?rep:replicate -> tree list -> ()


    (** {6 Testing Functions *)

    val analyze : tree list -> unit

end

module Make :
    functor (Node : NodeSig.S with type other_n = Node.Standard.n) ->
        functor (Edge : Edge.EdgeSig with type n = Node.n) ->
            functor (TreeOps : Ptree.Tree_Operations with type a = Node.n with type b = Edge.e) ->
                S with type a = Node.n with type b = Edge.e
