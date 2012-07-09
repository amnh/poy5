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

let () = SadmanOutput.register "ModelSelection" "$Revision"

let (-->) b a = a b
let failwithf format = Printf.ksprintf failwith format

module type S = sig

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

    val aic  : int -> int -> float -> float
    (** [aic n k l_max] Calculate the Akaike Information Criterion *)

    val aicc : int -> int -> float -> float
    (** [aic n k l_max] Calculate the AICc; as n->inf this is equal to aic, but has
        a higher parameter penalty for smaller n *)

    val bic  : int -> int -> float -> float
    (** [bic n k l_max] Calculate the Bayesian Information Criterion *)

    (** [best_model] return the optimial tree from stats *)
    val best_model : stats -> tree 

    (** [stats_of_models] fills and calculates the *IC stats for a set of
        models, used as a basis for all *IC methods. *)
    val stats_of_models : models -> ic -> Data.bool_characters -> tree -> stats

end

module Make (Node : NodeSig.S with type other_n = Node.Standard.n)
            (Edge : Edge.EdgeSig with type n = Node.n) 
            (TreeOps : Ptree.Tree_Operations with type a = Node.n with type b = Edge.e) : S =
struct
    type a = Node.n
    type b = Edge.e
    type tree = (a, b) Ptree.p_tree 

    (** {6 Types for collection of data and defining what type of IC to use *)

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

    (** {6 Helper functions for determining model properties *)

    (** list of all models for optimizing against; default *)
    let all_models =
        [ MlModel.JC69; MlModel.F81; MlModel.K2P None; MlModel.F84 None;
          MlModel.HKY85 None; MlModel.TN93 None; MlModel.GTR None; ]

    (** [parameter_cardinality] Helper function to determine the total number of
        parameters in the model; this includes branches, model rates, gamma, and
        for each character set if there are more than one. *)
    let parameter_cardinality tree chars : int =
        let sets =
            Data.categorize_likelihood_chars_by_model chars tree.Ptree.data in
        let model_params = 
            List.fold_left
                (fun acc xs ->
                    let model = Data.get_likelihood_model tree.Ptree.data xs in
                    acc + (MlModel.count_parameters model))
                0 sets
        and branch_params =
            (List.length sets) * (Tree.EdgeSet.cardinal tree.Ptree.tree.Tree.d_edges)
        in
        model_params + branch_params

    (** [sample_size] return the size of [n] in methods; this is to abstract the
        methodology used to determine [n] *)
    let sample_size tree chars : int = assert false

    (** [negative_loglikelihood] return the negative log-likelihood of a tree *)
    let negative_loglikelihood tree chars : float = assert false

    (** [diagnose_tree_with_model] determine the maximum likelihood value for
        the tree and model. Return the ML parameters and negative log-likelihood *)
    let diagnose_tree_with_model tree chars model : tree = assert false


    (** {6 General Information Metrics for Model Selection *)

    (** [aic] is defined as, AIC = 2 * k - 2 (ln L_max), where n is
        the number of observations, k are the number of parameters, and L_max is the
        log-likelihood of the model under the data. *)
    let aic k n nl_max =
        let k = float_of_int k and n = float_of_int n in
        (2.0 *. k) +. 2.0 *. nl_max

    (** [aicc] is defined as AICc = AIC + (2 k (k -1)) / (n-k-1); this criteria is the
        same as AIC as n -> inf (since n is in the denominator), as is a stricter
        parameter penalized version AIC. *)
    let aicc k n nl_max =
        let aic = aic k n nl_max and n = float_of_int n and k = float_of_int k in
        let c   = (2.0 *. k *. ( k +. 1.0)) /. (n -. k -. 1.0) in
        aic +. c

    (** [bic is defined as, BIC = ln(n) * k - 2 * ln(L_max), where n is the number
        of observations, k are the number of parameters, and L_max is the
        log-likelihood of the model under the data. *)
    let bic k n nl_max =
        let k = float_of_int k and n = float_of_int n in
        ((log n) *. k) +. 2.0 *. nl_max


    (** {6 Building the stats for Model Selection *)

    (** [best_model] return the optimial tree from stats *)
    let best_model stats =
        stats.tree_stats.( fst stats.min_ic ).tree

    (** [stats_of_tree] fill stats, minus decision theory and hlrt, of
        information criteria tests from a ML tree *)
    let stats_of_tree ic chars tree =
        let k = parameter_cardinality tree chars in
        let n = sample_size tree chars in
        let l = negative_loglikelihood tree chars in
        let ic = match ic with
            | AIC  -> aic k n l
            | AICC -> aicc k n l
            | BIC  -> bic k n l
        in
        tree, ic

    (** [stats_of_models] fills and calculates the *IC stats for a set of
        models, used as a basis for all *IC methods. *)
    let stats_of_models models ic chars tree =
        let tree_stats = 
            models --> Array.of_list
                   --> Array.map (diagnose_tree_with_model tree chars)
                   --> Array.map (stats_of_tree ic chars)
        in
        let _,min_ic =
            Array.fold_left
                (fun (i,((_,ma) as aa)) (_,ca) ->
                    let acc = if ca < ma then aa else (i,ca) in
                    (i+1,acc))
                (0,(~-1,max_float))
                (tree_stats)
        in
        let tree_stats, sum_ic =
            Array.fold_left
                (fun (data,sum_ic) (t,ic) ->
                    let w_ic  = exp (~-. (ic -. (snd min_ic)) /. 2.0) in
                    ((t,ic,w_ic)::data), sum_ic +. w_ic )
                ([],0.0)
                (tree_stats)
        in
        let tree_stats : tree_stats list =
            List.map
                (fun (t,ic,w_ic) ->
                    let d_ic = ic  -. (snd min_ic) in
                    let w_ic = w_ic /. sum_ic in
                    { tree = t;
                        lk = negative_loglikelihood t chars;
                        ic = (ic, d_ic, w_ic); })
                tree_stats
        in
        { tree_stats = Array.of_list tree_stats; min_ic = min_ic; }

end
