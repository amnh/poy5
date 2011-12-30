(* POY 5.0 Alpha. A phylogenetic analysis program using Dynamic Homologies.   *)
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

(** Numerical Functions *)

(** {6 Constants} *)

val tolerance : float
(** [Tolerance] Tolerance of floating point calculations for convergence. *)

val epsilon   : float
(** [Epsilon] Machine floating point precision. *)

val minimum   : float
(** [Minimum] Minimum for numerical calculations (2.0 * tolerance). *)

(** {6 Floating Point Functions} *)

val is_nan : float -> bool
(** is the floating point number NAN *)

val is_inf : float -> bool
(** is the floating point number +/- infinity *)

val is_zero : float -> bool
(** is the floating point number zero or sub-normal (effectively zero). *)

(** {6 Special Functions} *)

external gamma : float -> float = "gamma_CAML_gamma"
(** [gamma x] calculates the gamma of x. Using Lanczos Approximation, with
    error around 1e-15. Max value for overflow is around 142. *)

external lngamma : float -> float = "gamma_CAML_lngamma"
(** [lngamma x] calculates the ln gamma of x using Lanczos Approximation. *)

external gamma_rates: float -> float -> int -> 
    (float,Bigarray.float64_elt,Bigarray.c_layout) Bigarray.Array1.t = "gamma_CAML_rates"
(** [gamma_rates alpha beta cats] -> rates takes alpha, beta gamma parameters
    and number of categories to cut the gamma function into, and returns the mean
    rates in those cuts of 1/cats parts. *)

(** {6 Numerical Optimization Functions} *)

val brents_method :
    ?max_iter:int -> ?v_min:float -> ?v_max:float -> ?tol:float -> ?epsilon:float 
        -> float * ('a * float) -> (float -> 'a * float) -> float * ('a * float)
(** [brents_method ?i ?min ?max ?tol ?e o f] -> n Uses brents method, a
    combination of parabolic interpolation and golden section search, to find
    the local minimum near [o] of the function [f], bounded by [min] and [max].
    [i] and [tol] are used to control the number of iterations, tolerance,
    respectively. [o] is a pair of floating point numbers, with additional data,
    representing a point.  (See, Numerical Recipes in C; 10.2) *)

val line_search : 
    ?epsilon:float -> (float array -> 'a * float) -> float array -> 
        'a * float -> float array -> float -> float array -> float array * ('a * float) * bool
(** [line_search ?e ?a ?i ?min f p fp g s d] does a line search along the
    gradient [g] and direction [d] of function [f] by point [p], attempting the
    longest step, of maximum distance [s] *)

val bfgs_method :
    ?max_iter:int -> ?epsilon:float -> ?mx_step:float -> ?g_tol:float ->
        (float array -> 'a * float) -> float array -> 'a * float -> float array * ('a * float)
(** [bfgs_method ?i ?e ?s f init] uses bfgs method to approximate the hessian
    matrix of a function f with starting point init. *)

(** {6 Infix Module} *)

module type I =
    sig
        val set_eps : float -> unit
        val set_ops : int -> unit
        val reset   : unit -> unit

        val (=.) : float -> float -> bool
        val (>.) : float -> float -> bool
        val (<.) : float -> float -> bool
    end
module FPInfix : I
(** This module creates fuzzy infix operators for floating point numbers. This
    still must be used with care since they will not hold transitivity. Epsilon
    can be set directly or calculated from the machine precision and the number
    of operations where error has accumulated. Ensure that this modification is
    reverted, or the change will persist. These functions are probably only
    useful on floating point alignment --when many operations are going to
    accumulate a large number of ulps. *)
