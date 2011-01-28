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

(** [string_spec] primative representation of the spec used to fold over a list
 * of properties in nexus format and build a complete model. 
 *
 * the type follows: 
 *      modelname,(variation,#sites,alpha,%invar),params,priors,gap?,file *)
type string_spec = string * (string * string * string * string) * float list
                 * [ `Given of (string * float) list | `Equal | `Estimate of float array option ]
                 * (string * float option) * string * string option

val empty_str_spec : string_spec

val likelihood_not_enabled : string

(** [site_var] the site variation rate parameters. *)
type site_var = 
    (** [Gamma] #categories, alpha, beta *)
    | Gamma of int * float
    (** [Theta] #categories, alpha, beta, %invar *) 
    | Theta of int * float * float
    (* | Given of (float * float) array *)
    (** [Constant] constant rate parameter *)
    | Constant 

(* [subst_model] the model defining the substitution rate matrix. An extra
 * parameter is present (as in JC69) and represents a scale factor for the branch
 * lengths of the resultant tree. *)
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

(** [priors] the prior probabilities. This is only used in specification, in
 * favor of a C-type big array for compatibility.  *)
type priors = 
    | Estimated  of float array
    | Given      of float array
    | Equal

type gap_properties = [ `Missing | `Independent | `Coupled of float ]

(** [spec] the specification of the model. *)
type spec = {
    substitution : subst_model;
    site_variation : site_var option;
    base_priors : priors;
    cost_fn : Methods.ml_costfn;
    use_gap : gap_properties;
    iterate_model : bool;
    iterate_alpha : bool;
}

type model = {
    (** [spec] specification this model was created from *)
    spec  : spec;
    (** [alph] size of the alphabet for the model *)
    alph : Alphabet.a;
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


val jc69_5 : spec
(** sample spec for testing of 5 state JC69 *)

val jc69_4 : spec
(** sample spec for testing of 4 state JC69 *)

type chars = [ `List of int list | `Packed of int ]
(** list of set bits, or packed integer of set bits *)

val list_of_packed : ?zerobase:bool -> int -> int list
(** convert packed integer to a list of indexes *)

val packed_of_list : int list -> int
(** convert a list of indexes to a packed integer *)

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
    all the sequences to estimate additional gaps in dynamic likelihood *)

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

(** [enable_gaps f m] a function to add gaps to the model; this is used to
 * ensure that dynamic characters have gaps as a character during transforms.
 * The compute priors function is only used if the previous model had estimated
 * priors, and they need to be computed again *)
val add_gap_to_model : (unit -> float array) -> model -> model

val remove_gamma_from_spec : spec -> spec

val compare : model -> model -> int

(** [classify_seq_pairs l1 l2 seq1 seq2] does initial classification of the
* transitions from [seq1] to [seq2] Returns a map of complements, the number of
* of exact matches, and total number of transitions. Only base frequency in the
* leaves are added to the IntegerMap, indicated by [l1] and [l2]. *)
val classify_seq_pairs :
    bool -> bool -> (float * int * chars) list -> (float * int * chars) list ->
        (float All_sets.FullTupleMap.t) * (float All_sets.IntegerMap.t) ->
            (float All_sets.FullTupleMap.t) * (float All_sets.IntegerMap.t)

val spec_from_classification :
    Alphabet.a -> Methods.ml_gap -> Methods.ml_substitution -> Methods.ml_site_variation option -> 
        Methods.ml_costfn -> (float All_sets.FullTupleMap.t) * (float All_sets.IntegerMap.t) -> spec

(** [compse m t] compose a matrix with the time *)
val compose : model -> float -> 
    (float, Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array2.t

val model_to_cm : model -> float -> Cost_matrix.Two_D.m

val output_model : (string -> unit) -> [`Nexus | `Phylip | `Hennig] -> model -> string list option -> unit

val to_formatter : model -> Xml.xml Sexpr.t list

(* [get_update_function_for_model] based on the model provided return a function
 * that will update the model based on an input value *)
val get_update_function_for_model    : model -> (model -> float array -> model) option
val get_current_parameters_for_model : model -> float array option

(* [get_update_function_for_alpha] based on the alpha parameters in model *)
val get_update_function_for_alpha    : model -> (model -> float -> model) option
val get_current_parameters_for_alpha : model -> float option
