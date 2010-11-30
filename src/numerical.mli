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

(** These constants define the ending conditions for numerical functions, and
 * provide a common baseline throughout the program for comparing floating point
 * numbers.
 *  Tolerance   - Tolerance of floating point calculations for convergence.
 *  Epsilon     - Machine floating point precision.
 *  Minimum     - minimum for numerical calculations (2.0 * tolerance). *)
val tolerance : float
val epsilon   : float
val minimum   : float

(** [gamma_rates alpha beta cats] -> rates
 * takes alpha, beta gamma parameters and number of categories to cut the gamma
 * function into, and returns the mean rates in those cuts of 1/cats parts. *)
external gamma_rates: float -> float -> int -> 
    (float,Bigarray.float64_elt,Bigarray.c_layout) Bigarray.Array1.t = "gamma_CAML_rates"

(** [brents_method ?i ?min ?max ?tol ?e o f] -> n 
 * Uses brents method, a combination of parabolic interpolation and golden
 * section search, to find the local minimum near [o] of the function [f],
 * bounded by [min] and [max]. [i] and [tol] are used to control the number
 * of iterations, tolerance, respectively. [o] is a pair of floating point
 * numbers, with additional data, representing a point.
 *
 * Numerical Recipes in C; 10.2 *)
val brents_method :
    ?max_iter:int -> ?v_min:float -> ?v_max:float -> ?tol:float -> ?epsilon:float 
        -> float * ('a * float) -> (float -> 'a * float) -> float * ('a * float)

(* [line_search ?e ?a ?i ?min f p fp g s d] does a line search along the
 * gradient [g] and direction [d] of function [f] by point [p], attempting the
 * longest step, of maximum distance [s] *)
val line_search : 
    ?epsilon:float -> (float array -> 'a * float) -> float array -> 
        'a * float -> float array -> float -> float array -> float array * ('a * float) * bool

(* [bfgs_method ?i ?e ?s f init] uses bfgs method to approximate the hessian
 * matrix of a function f with starting point init. *)
val bfgs_method :
    ?max_iter:int -> ?epsilon:float -> ?mx_step:float -> ?g_tol:float ->
        (float array -> 'a * float) -> float array -> 'a * float -> float array * ('a * float)

