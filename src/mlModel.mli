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

type string_spec = string * (string * string * string * string) * float list
                 * [ `Given of (string * float) list | `Equal | `Estimate of float array option ]
                 * (string * float option) * string * string option
(** [string_spec] primative representation of the spec used for parsing as an
    intermediary step
    the type follows: 
        modelname,(variation,#sites,alpha,%invar),params,priors,gap?,file *)

val empty_str_spec : string_spec
(** an empty string spec, contains no default values for options *)

val likelihood_not_enabled : string
(** string for erroring when likelihood is not enabled. *)

type site_var = 
    (** [Gamma] #categories, alpha, beta *)
    | Gamma of int * float
    (** [Theta] #categories, alpha, beta, %invar *) 
    | Theta of int * float * float
    (* | Given of (float * float) array *)
    (** [Constant] constant rate parameter *)
    | Constant 
(** [site_var] the site variation rate parameters. *)

type subst_model =
    | JC69
    | F81
    | K2P   of float option
    | F84   of float option
    | HKY85 of float option
    | TN93  of (float * float) option
    (** [GTR] alphabetical order describtion of the transition rates *)
    | GTR   of (float array) option
    (** [File] matrix read from a file, diagonal is readjusted so row = 0 *)
    | File  of float array array * string
(** [subst_model] the model defining the substitution rate matrix. **)

type priors = 
    | Estimated  of float array
    | Given      of float array
    | Equal
(** [priors] the prior probabilities, or base frequencies. **)

type gap_properties = [ `Missing | `Independent | `Coupled of float ]
(** How to deal with gaps; are they missing? (the traditional implemenatation),
    are their rates consistent with the model, or do are their rates coupled
    with an additional parameter. *)

type spec = {
    substitution : subst_model;
    site_variation : site_var option;
    base_priors : priors;
    cost_fn : Methods.ml_costfn;
    use_gap : gap_properties;
    iterate_model : bool;
    iterate_alpha : bool;
}
(** [spec] the specification of the model. *)

type model = {
    (** [spec] specification this model was created from *)
    spec  : spec;
    (** [alph] the alphabet. *)
    alph : Alphabet.a;
    (** [alph_s] size of the alphabet for the model. *)
    alph_s : int;
    (** [pi_0] priors for the model *)
    pi_0  : (float, Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array1.t;
    (** [invar] the percent used for the proportion of invariant sites *)
    invar : float option;
    (** [rate] rate classes *)
    rate  : (float, Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array1.t;
    (** [prob] probabilities for the rate classes *)
    prob  : (float, Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array1.t;
    (** [s] the instantaneous rate matrix of the model *)
    s     : (float, Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array2.t;
    (** [u] the left eigenvectors of the substitution rate matrix *)
    u     : (float, Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array2.t;
    (** [d] the eigenvalues as a diagonal matrix *)
    d     : (float, Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array2.t;
    (** [ui] the inverse of u. In symmetric cases this isn't needed, ui=ut *)
    ui    : (float, Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array2.t option; 
}
(** [model] type to define the model used in the calculation of the likelihood
   vector, and the loglikelihood value associated with the result. We also hold
   the specification of the model for output and possible iteration. *) 

val jc69_5_gap : float -> spec
(** sample spec for testing of 5 state JC69 with different gap rate *)

val jc69_5 : spec
(** sample spec for testing of 5 state JC69 *)

val jc69_4 : spec
(** sample spec for testing of 4 state JC69 *)

type chars = [ `List of int list | `Packed of int ]
(** list of set bits, or packed integer of set bits *)

val get_costfn_code : model -> int
(** code for costfn for c-side *)

val categorize_by_model : 'a list -> ('a -> spec) -> 'a list list
(** categorize a list of values into a list list of values; usually codes *)

val replace_priors : model -> float array -> model 
(** replace the priors in the model; used in dynamic likelihood to reestimate
    when an implied alignment is performed *)

val compute_priors :
    (Alphabet.a * bool) -> float array -> (int * int) -> int list -> float array
(** compute priors from some basic information; the alphabet, if we are
    counting gaps as an additional character, the frequency of the bases in an
    array, the total number of bases, number of gaps counted, and the lengths of
    all the sequences to estimate additional gaps in dynamic likelihood. The
    frequencies of the alphabet is calculated earlier by functions that
    understand the current representation of the data. *)

val compare_priors : model -> model -> bool
(** compare two sets of priors *)

val convert_string_spec : string_spec -> spec
(** [convert_string_spec] convert a string spec from nexus and other formats to
   the basic specification in for a likelihood model *)

val convert_methods_spec : int -> (unit -> float array) -> Methods.ml_spec -> spec
(** [convert_methods_spec] convert the specification from Methods into the
    proper MlModel.spec; this then can be converted to MlModel.model *)

val create : ?min_prior:float -> Alphabet.a -> spec -> model
(** [create_lk_model s] create the model for likelihood from parser *)

val add_gap_to_model : (unit -> float array) -> model -> model
(** [enable_gaps f m] a function to add gaps to the model; this is used to
   ensure that dynamic characters have gaps as a character during transforms.
   The compute priors function is only used if the previous model had estimated
   priors, and they need to be computed again *)

val remove_gamma_from_spec : spec -> spec
(** [remove_gamma_from_spec] set the gamma property to None/Constant. This is
    for un-partitioned dynamic likelihood characters *)

val compare : model -> model -> int
(** compare two models; not metric *)

val classify_seq_pairs :
    bool -> bool -> (float * int * chars) list -> (float * int * chars) list ->
        (float All_sets.FullTupleMap.t) * (float All_sets.IntegerMap.t) ->
            (float All_sets.FullTupleMap.t) * (float All_sets.IntegerMap.t)
(** [classify_seq_pairs l1 l2 seq1 seq2] does initial classification of the
    transitions from [seq1] to [seq2] Returns a map of complements, the number of
    of exact matches, and total number of transitions. Only base frequency in the
    leaves are added to the IntegerMap, indicated by [l1] and [l2]. *)

val spec_from_classification :
    Alphabet.a -> Methods.ml_gap -> Methods.ml_substitution -> Methods.ml_site_variation option -> 
        Methods.ml_costfn -> (float All_sets.FullTupleMap.t) * (float All_sets.IntegerMap.t) -> spec
(** Create the specification from the classification of the columns; above *)

val compose : model -> float -> 
    (float, Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array2.t
(** [compse m t] compose a matrix with the time *)

val subst_matrix : model -> float option -> 
    (float, Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array2.t
(** [subst_matrix] return the substitution rate matrix; Q matrix if t = None,
    else Qt. *)

val debug_model : model -> float option -> unit

val check_metricity : model -> float -> float ->
    (float, Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array2.t ->
    (float, Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array2.t -> unit

val model_to_cm : model -> float -> Cost_matrix.Two_D.m
(** return an integerized cost matrix of the model *)

val output_model : (string -> unit) -> [`Nexus | `Phylip | `Hennig] -> model -> string list option -> unit
(** print the model in a format similar to Phylip's output, or a formatted nexus
    output. Hennig output is the same as Phylip, and used to match treeoutput
    functions. *)

val to_formatter : model -> Xml.xml Sexpr.t list
(** return XML reprentation of the model *)

val get_update_function_for_model    : model -> (model -> float array -> model) option
(** [get_update_function_for_model] based on the model provided return a function
    that will update the model based on an input value *)

val get_current_parameters_for_model : model -> float array option
(** [get_current_parameters_for_model] will return a properly formatted array
    based on the model to use for further updates. *)

val get_update_function_for_alpha    : model -> (model -> float -> model) option
(** [get_update_function_for_alpha] based on the alpha parameters in model *)

val get_current_parameters_for_alpha : model -> float option
(** [get_current_parameters_for_alpha] return the alpha parameter if it exists *) 
