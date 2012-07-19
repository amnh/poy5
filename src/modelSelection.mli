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

    (** Defines the properties of a model to update **)
    type model_diff =
        { subst : MlModel.subst_model;
          rates : MlModel.site_var;
          prior : MlModel.priors;}

    (** Define a list/subset of models to use in the particular IC *)
    type model_diffs = model_diff list

    (** Types of information criteria used in the tree-stats *)
    type ic = [ `AIC | `AICC | `BIC ]

    (** Define stats for a tree; likelihood, aic, and bic *)
    type tree_stats =
        {   tree : tree;
              lk : float;
                 (* IC, delta_IC, IC_weight *)
              ic : float * float * float; }

    (** Define a collection of models for information testing *)
    type stats =
        { tree_stats : tree_stats array;
             type_ic : ic;
              min_ic : int * float; }

    (** Helper to create a list of all combinations of subsitutions and rate
        variations from a candidate spec. *)
    val select_models : tree -> Methods.ml_spec -> MlModel.model * model_diffs

    (** {6 Information metrics for Model selection *)

    (** [aic] is defined as, AIC = 2 * k - 2 (ln L_max), where n is
        the number of observations, k are the number of parameters, and L_max is the
        log-likelihood of the model under the data. *)
    val aic  : int -> int -> float -> float

    (** [aicc] is defined as AICc = AIC + (2 k (k -1)) / (n-k-1); this criteria is the
        same as AIC as n -> inf (since n is in the denominator), as is a stricter
        parameter penalized version AIC. *)
    val aicc : int -> int -> float -> float

    (** [bic is defined as, BIC = ln(n) * k - 2 * ln(L_max), where n is the number
        of observations, k are the number of parameters, and L_max is the
        log-likelihood of the model under the data. *)
    val bic  : int -> int -> float -> float


    (** {6 Return stats of trees for information criteria *)

    (** [best_model] return the optimial tree from stats *)
    val best_model : stats -> tree

    (** [stats_of_models] fills and calculates the *IC stats for a set of
        models, used as a basis for all *IC methods. *)
    val stats_of_models :
        model_diffs -> ic -> Data.bool_characters -> MlModel.model -> tree -> stats

    (** [merge_stats] merges a list of tree_stats into and updates the info *)
    val merge_stats : stats list -> stats
   
    (** [report_stats] report the stats in a chart; ic, weight, lk, ... *)
    val report_stats : stats -> Data.bool_characters -> string array array

    (** [model_average_parameter] return the model average parameter for a
        function that access the given parameter.  The function should return
        "None" if the model should be excluded from the model averaging. *)
    val model_averaged_parameter :
        (MlModel.model -> float option) -> stats -> Data.bool_characters
            -> float * float

    (** [generate_stats] built from the natural Methods.spec from the parser to
        create a stats object; most convenient entry point. *)
    val generate_stats : tree -> Methods.ml_spec -> stats
        
end

module Make :
    functor (Node : NodeSig.S with type other_n = Node.Standard.n) ->
        functor (Edge : Edge.EdgeSig with type n = Node.n) ->
            functor (TreeOps : Ptree.Tree_Operations with type a = Node.n with type b = Edge.e) ->
                S with type a = Node.n with type b = Edge.e

(** Load a run script that should load data, build trees and transform to
    likelihood. run method does not need to be iterative, adjusted in method. *)
(*val test : string -> unit*)
