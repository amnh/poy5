(* POY 5.0 Alpha. A phylogenetic analysis program using Dynamic Homologies.   *)
(* Copyright (C) 2011  Andrés Varón, Lin Hong, Nicholas Lucaroni, Ward Wheeler*)
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

val pi : float
(** [pi] IEEE float closest to the true value of pi; (defined, 4.0 * atan 1.0) *)

(** {6 Types} *)

type 'a simplex
(** How we internally represent a simplex, keeping it abstract, of course *)

type simplex_strategy = 
    {   alpha : float;  (** The Reflection factor *)
         beta : float;  (** The Contraction factor *)
        gamma : float;  (** The Expansion factor *)
        delta : float;  (** The Shrinkage (Massive Contraction) factor *)
    } 
(** Simplex strategy defines how to expand, contract, and reflect a simplex *)

type subplex_strategy = 
    { simplex : simplex_strategy;   (** How to perform the simplex *)
          psi : float;              (** The Simplex Reduction coefficient *)
        omega : float;              (** The Step Reduction coefficient *)
        nsmin : int;                (** Minimum subspace dimension; or 2 *)
        nsmax : int;                (** Maximum subspace dimension; or 5 *)
    }
(** The Subplex strategy contains a simplex strategy and added features *)


(** {6 Floating Point Functions} *)

val is_nan : float -> bool
(** is the floating point number the IEEE representation of NAN *)

val is_inf : float -> bool
(** is the floating point number +/- infinity *)

val is_zero : float -> bool
(** is the floating point number zero or sub-normal (effectively zero). *)


(** {6 Special Functions} *)

val gamma : float -> float
(** [gamma x] calculates the gamma of x. Using Lanczos Approximation, with
    error around 1e-15. Max value for overflow is around 142. *)

val lngamma : float -> float
(** [lngamma x] calculates the ln gamma of x using Lanczos Approximation. *)

val gamma_rates:
    float -> float -> int -> 
        (float,Bigarray.float64_elt,Bigarray.c_layout) Bigarray.Array1.t
(** [gamma_rates alpha beta cats] -> rates takes alpha and beta parameters of
    gamma and number of categories to cut the gamma function into, and returns
    the mean rates in those cuts of 1/cats parts. *)

val rand_normal : float -> float -> float
(** [rand_normal m stdev] generate a random variable with a normal distribution
    with the mean [m] and standard deviation [stdev]. **)

val rand_exp : float -> float
(** [rand_exp m] Generates a random variable with the exponential distribution
    of mean [m]. *)

val rand_gamma : float -> float -> float
(** [rand_gamma shape scale] Generates a random variable with a gamma
    distribuation of the shape and scale parameters passed. *)

val pnorm : ?mean:float -> ?sd:float -> float -> float
(** Cumulative distribution function of the normal distribution; \Phi *)

val dnorm : ?mean:float -> ?sd:float -> float -> float
(** Density function for normal distribution function; \phi *)

val qnorm : ?mean:float -> ?sd:float -> float -> float
(** Quantile function for the normal distribution function; \Phi^{-1} *)


(** {6 Statisical Objects *)

class running_stats : object
(** [running_stats] is an object that can calculate the variance,
    standard-deviation, min, max, and mean from a set of data of unknown size,
    and will update for each additional number. Use [push] to add numbers, then
    other functions will update accordingly. *)

    method push : float -> unit
    (** Push a new value to the running stats object *)

    method iter : unit -> int
    (** the current number of iterations *)

    method min  : unit -> float
    (** the current minimum found in the data-set *)

    method max  : unit -> float
    (** the current maximum found in the data-set *)

    method variance : unit -> float
    (** calculate the current variance in the data-set *)

    method mean     : unit -> float
    (** calculate the current mean of the data-set *)

    method std_dev  : unit -> float
    (** calculate the current standard-deviation of the data-set *)
end


(** {6 Numerical Optimization Functions} *)

val brents_method :
    ?max_iter:int -> ?v_min:float -> ?v_max:float -> ?tol:float -> ?epsilon:float 
        -> (float -> 'a * float) -> float * ('a * float) -> float * ('a * float)
(** [brents_method ?i ?min ?max ?tol ?e o f] -> n Uses brents method, a
    combination of parabolic interpolation and golden section search, to find
    the local minimum near [o] of the function [f], bounded by [min] and [max].
    [i] and [tol] are used to control the number of iterations, tolerance,
    respectively. [o] is a pair of floating point numbers, with additional data,
    representing a point.  (See, Numerical Recipes in C; 10.2) *)

val brents_method_multi :
    ?max_iter:int -> ?v_min:float -> ?v_max:float -> ?tol:float -> ?epsilon:float 
        -> (float array -> 'a * float) -> (float array * ('a * float))
            -> (float array * ('a * float))
(** brents_method_multi ...] Generalization of the above function for brents
    method. We optimize each value in the array one by one and only ONCE. This
    was seen in RAxML 7.04; and possibility has some utility since we do not
    need to calculate the derivative of a vector which can be costly and may not
    work on routines with a number of discontinuities. *)

(* Not necessary to be exposed; here in case a situation requires it.
val line_search :
    ?epsilon:float -> (float array -> 'a * float) -> float array -> 
        'a * float -> float array -> float -> float array -> float array * ('a * float) * bool
(** [line_search ?e ?a ?i ?min f p fp g s d] does a line search along the
    gradient [g] and direction [d] of function [f] by point [p], attempting the
    longest step, of maximum distance [s] *) 
*)

val bfgs_method :
    ?max_iter:int -> ?epsilon:float -> ?max_step:float -> ?tol:float
        -> (float array -> 'a * float) -> (float array * ('a * float))
            -> (float array * ('a * float))
(** [bfgs_method ?i ?e ?s f init] uses bfgs method to approximate the hessian
    matrix of a function f with starting point init. *)

val simplex_method : 
    ?termination_test : (float -> 'a simplex -> bool) -> ?tol : float ->
    ?simplex_strategy : simplex_strategy -> ?max_iter:int -> ?step:float array option
        -> (float array -> 'a * float) -> (float array * ('a * float))
            -> (float array * ('a * float))
(** The simplex uses an n+1 dimensional figure to move, like an amoeba, around
    the surface. The extermities of the object are evaluated, where the maximal
    values are modified, while the others are progressed with improvements. The
    four moves that can be done on the simplex are, reflection, expansion,
    contraction, and shrinkage. The degree to which these are done is modified
    by the strategy used. *)

val subplex_method : 
    ?subplex_strategy : subplex_strategy -> ?tol:float -> ?max_iter:int
        -> (float array -> 'a * float) -> (float array * ('a * float))
            -> (float array * ('a * float))
(** The Subplex method is a generalization to a number of algorithms, paramount
    the Nelder-Mead simplex method, with alternating variables, and Nelder-Mead
    Simplex with restart. The advantages are outlined in the previously
    mentioned paper, in section 5.3.6. *)

(** {6 Unification of Optimization Strategies *)

type optimization_strategy = 
    {   routine : [ `Simplex of simplex_strategy option
                  | `Subplex of subplex_strategy option
                  | `BFGS    of float option
                  | `Brent_Multi ];
        (* Below are optional to use the algorithm default *)
        max_iter : int option;
        tol      : float option; 
    }
(** Define an optimization strategy. This will call the appropriate optimization
    routine with the specified parameters. A wrapper around the four methods we
    have to optimize multi-dimensional functions. *)

type opt_modes =
    [ `None
    | `Coarse of int option
    | `Exhaustive of int option
    | `Custom of optimization_strategy list ]

val compare_opt_mode : opt_modes -> opt_modes -> int 
(** Provides an ordering of optimization modes from least intense to more intense *)

val default_numerical_optimization_strategy :
    opt_modes -> int -> optimization_strategy list
(** Take the cost model of the Methods module and determine an opt strategy *)

val default_number_of_passes : opt_modes -> int
(** Determine the default number of passes, in an optimization routine to do
    under the specified optimization level *)

val set_tol_c_brents_method : float -> unit
(** Set the C tolerance in brents method; for optimization of branch lenghths *)

val get_tol : opt_modes -> float
(** Get the tolerance for that optimization level; used to communicate with
    external routines --for example, the C-brents method routine. *)

val run_method :
    optimization_strategy list -> (float array -> 'a * float)  
        -> (float array * ('a * float)) -> (float array * ('a * float))
(** Run an optimization strategy based on the optimization_strategy. We accept a
    list, thus multiple tests can be run to obtain a robust strategy. *)


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
