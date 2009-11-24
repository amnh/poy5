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
type string_spec = string * (string * string * string * string)
                          * float list * ( string * float ) list * string * string option

val empty_str_spec : string_spec

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
    | File  of float array array 

(** [priors] the prior probabilities. This is only used in specification, in
 * favor of a C-type big array for compatibility.  *)
type priors = 
    | Estimated of float array
    | Given     of float array

(** [spec] the specification of the model. *)
type spec = {
    substitution : subst_model;
    site_variation : site_var option;
    base_priors : priors;
    use_gap : bool;
}

(** [model] type to define the model used in the calculation of the likelihood
 * vector, and the loglikelihood value associated with the result. We also hold
 * the specification of the model for output and possible iteration. *) 
type model = {
    (** [spec] specification this model was created from *)
    spec  : spec;
    (** [alph] size of the alphabet for the model *)
    alph : int;
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

IFDEF USE_LIKELIHOOD THEN

val compare : model -> model -> int

(** [diagonalize_*** Q D [Ui] ] 
 * Diagonalize [Q], and places the eigenvalues along the diagonal of [D],
 * thus [D] and [Q] must be nxn --where n is the size of the alphabet. The
 * eigenvectors are entered along the rows of [Q], this function destroyes
 * the substution rate matrix, but that can be copied or the primative
 * arguments saved and this matrix can be reconstructed later --it should
 * be fairly cheap to recreate *)
external diagonalize_gtr: (* U D Ui *) FMatrix.m ->
    (float,Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array2.t -> 
    (float,Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array2.t -> 
    (float,Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array2.t ->
    unit = "likelihood_CAML_diagonalize_gtr"
external diagonalize_sym: (* U D *) FMatrix.m ->
    (float,Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array2.t -> 
    (float,Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array2.t -> 
    unit = "likelihood_CAML_diagonalize_sym"
(** [compose_*** U D [Ui] t]
 * Composes the probability matrix from it's parts P = U*exp(t*D)*Ui
 * Function is used for testing and to_formatter function usage (output) *)     
external compose_gtr: (* U D Ui t -> P *) FMatrix.m ->
    (float,Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array2.t -> 
    (float,Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array2.t ->
    (float,Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array2.t -> float
    -> (float,Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array2.t = 
        "likelihood_CAML_compose_gtr"
external compose_sym: (* U D t -> P *) FMatrix.m ->     
    (float,Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array2.t -> 
    (float,Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array2.t -> float
    -> (float,Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array2.t = 
        "likelihood_CAML_compose_sym" 

(** [gamma_rates alpha beta cats] -> rates
 * takes alpha, beta gamma parameters and number of categories to cut the gamma
 * function into, and returns the mean rates in those cuts of 1/cats parts. *)
external gamma_rates: float -> float -> int -> 
    (float,Bigarray.float64_elt,Bigarray.c_layout) Bigarray.Array1.t =
        "gamma_CAML_rates"

(** models to be used outside likelihood if necessary **)
val m_gtr   : float array -> float array -> int ->
    (float,Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array2.t
val m_f84   : float array -> float -> float -> int ->
    (float,Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array2.t
val m_hky85 : float array -> float -> int ->
    (float,Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array2.t
val m_f81   : float array -> float -> int -> 
    (float,Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array2.t
val m_tn93  : float array -> float -> float -> float -> int -> 
    (float,Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array2.t
val m_k2p   : float array -> float -> float -> int -> 
    (float,Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array2.t
val m_jc69  : float array -> float -> int -> 
    (float,Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array2.t
val m_file  : float array -> float array array -> int -> 
    (float,Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array2.t

(** [convert_string_spec] convert a string spec from nexus and other formats to
 * the basic specification in for a likelihood model *)
val convert_string_spec : string_spec -> spec

(** [create_lk_model s] create the model for likelihood from parser *)
val create  : Alphabet.a -> spec -> model

(** [estimate_model_pair l1 l2 seq1 seq2] does initial classification of the
* transitions from [seq1] to [seq2] Returns a map of complements, the number of
* of exact matches, and total number of transitions. Only base frequency in the
* leaves are added to the IntegerMap, indicated by [l1] and [l2]. *)
type chars = [ `List of int list | `Packed of int ]
val classify_seq_pairs :
    bool -> bool -> (float * int * chars) list -> (float * int * chars) list ->
        (float All_sets.FullTupleMap.t) * (float All_sets.IntegerMap.t) ->
            (float All_sets.FullTupleMap.t) * (float All_sets.IntegerMap.t)

val spec_from_classification :
    Alphabet.a -> bool -> Methods.ml_substitution -> Methods.ml_site_variation option ->
        (float All_sets.FullTupleMap.t) * (float All_sets.IntegerMap.t) -> spec

(** [compse m t] compose a matrix with the time *)
val compose : model -> float -> 
    (float, Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array2.t

(** [diagonalize s m] diagonalize [m] *)
val diagonalize :
    bool ->
    (float, Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array2.t ->
    (float, Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array2.t *
    (float, Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array2.t *
    (float, Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array2.t option

(** [compose_model s t] compose a matrix from substitution rate matrix *)
val compose_model : (float, Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array2.t 
        -> float -> (float, Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array2.t

(** [brents_method ?i ?e o f] uses brents method of parabolic interpolation to
 * find the local minimum near [o] of the function [f]. [i] and [e] are used to
 * control the number of iterations and tolerance, respectively. [o] is a pair
 * of floating point numbers, essentially, representing a point *)
val brents_method :
    ?iter_max:int -> ?epsilon:float 
        -> float * float -> (float -> float) -> float * float

(* [line_search ?e ?a ?i ?min f p fp g s d] does a line search along the
 * gradient [g] and direction [d] of function [f] by point [p], attempting the
 * longest step with maximum [s] *)
val line_search : 
    ?epsilon:float -> ?alf:float -> (float array -> 'a * float) -> float array 
        -> 'a * float -> float array -> float -> float array -> float array * ('a * float) * bool

(* [bfgs_method ?i ?e ?s f init] uses bfgs method to approximate the hessian
 * matrix of a function f with starting point init. *)
val bfgs_method :
    ?max_iter:int -> ?epsilon:float -> ?mx_step:float -> ?g_tol:float ->
        (float array -> 'a * float) -> float array -> 'a * float -> float array * ('a * float)

(** [update_k2p m v] update model [m] with new value [v] *)
val update_k2p : model -> float -> model

(** [update_hky m v] update model [m] with new value [v] *)
val update_hky : model -> float -> model

(** [update_f84 m v] update model [m] with new value [v] *)
val update_f84 : model -> float -> model

(** [update_gtr m v] update model [m] with new value [v] *)
val update_gtr : model -> float array -> model

(** [update_alpha m v] update model [m] with rate parameter [v] *)
val update_alpha : model -> float -> model

(* [get_update_function_for_model] based on the model provided return a function
 * that will update the model based on an input value *)
val get_update_function_for_model : model -> (model -> float array -> model) option
val get_current_parameters_for_model : model -> float array option

END
