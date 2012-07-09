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


(** {6 Define the functor of nodes and edges and tree operations *)

module type S = sig

    (** {6 Types *)

    (** abstract type for node in functorization *)
    type a 

    (** abstract type for edge in functorization *)
    type b

    (** tree composed of a and b *)
    type tree = (a,b) Ptree.p_tree

    (** Define a list/subset of models to use in the particular IC *)
    type models = MlModel.subst_model list

    (** Types of Information Criterias used in the tree-stats *)
    type ic = | AIC | AICC | BIC

    (** Define stats for a tree; likelihood, aic, and bic *)
    type tree_stats =
        {   tree : tree;
              lk : float;
                 (* IC, delta_IC, IC_weight *)
              ic : float * float * float; }

    (** Define a collection of models for information testing *)
    type stats =
        { tree_stats : tree_stats array;
              min_ic : int * float; }


    (** {6 Information metrics for Model selection *)

    val aic  : int -> int -> float -> float
    (** [aic n k l_max] Calculate the Akaike Information Criterion *)

    val aicc : int -> int -> float -> float
    (** [aic n k l_max] Calculate the AICc; as n->inf this is equal to aic, but has
        a higher parameter penalty for smaller n *)

    val bic  : int -> int -> float -> float
    (** [bic n k l_max] Calculate the Bayesian Information Criterion *)


    (** {6 Return stats of trees for information criterias *)

    (** [best_model] return the optimial tree from stats *)
    val best_model : stats -> tree 

    (** [stats_of_models] fills and calculates the *IC stats for a set of
        models, used as a basis for all *IC methods. *)
    val stats_of_models : models -> ic -> Data.bool_characters -> tree -> stats


end
