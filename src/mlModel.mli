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

type site_var = 
    (* #categories, alpha, beta, %invar *)
    | Gamma of int * float * float
    | Theta of int * float * float * float
    | Constant 
type subst_model =
    | JC69  of float
    | F81   of float
    | K2P   of float * float
    | F84   of float * float
    | HKY85 of float * float
    | TN93  of float * float * float
    | GTR   of float array
    | File  of float array array 
type priors = 
    | Estimated of float array
    | Given     of float array

(* specification of a model *)
type spec = {
    substitution : subst_model;
    site_variation : site_var option;
    base_priors : priors;
    use_gap : bool;
}

(* type of a model *)
type model = {
    name  : string;
    pi_0  : (float, Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array1.t;
    alpha : float option;
    invar : float option;
    sites : int;
    rate  : (float, Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array1.t;
    prob  : (float, Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array1.t;
    u     : (float, Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array2.t;
    d     : (float, Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array2.t;
    ui    : (float, Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array2.t option; 
}

(** models to be used outside likelihood if necessary **)
val m_gtr   : float array -> float array -> int ->
    (float,Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array2.t
val m_f84   : float array -> float -> float -> int ->
    (float,Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array2.t
val m_hky85 : float array -> float -> float -> int ->
    (float,Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array2.t
val m_hky85_ratio : float array -> float -> int ->
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

(** [create_lk_model s] create the model for likelihood from parser *)
val create  : Alphabet.a -> spec -> model

(** [compse m t] compose a matrix with the time *)
val compose : model -> float -> (float, Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array2.t

(** [compose_model s t] compose a matrix from substitution rate matrix *)
val compose_model : (float, Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array2.t 
        -> float -> (float, Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array2.t
