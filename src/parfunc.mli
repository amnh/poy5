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

(* $Id: parfunc.mli 1644 2007-02-14 19:05:47Z andres $ *)

(** Parallel operations *)

(** [all_minimize f c a] applies the function f on each element in a and
* calculates the cost of the resulting value using c. Minimum of all the
* calculated values is returned.
* The function can raise a Process.All_Children_Dead exception if there are no
* more available processes to continue the computation. *)
val all_minimize : Register.all_parallel -> Register.all_parallel -> 'a array ref -> float

(** [all_maximize f c a] provides the same functionality of all_minimize but
* finds the maximum instead of the minimum. 
* The function can raise a Process.All_Children_Dead exception if there are no
* more available processes to continue the computation. *)
val all_maximize : Register.all_parallel -> Register.all_parallel -> 'a array ref -> float

(** [gather_minimum f c a th] applies the function f to the elements in a,
* calculates the cost each element in the resulting array using the function c
* and returns the subarray of elements at a distance less than the minimum + th
* (a threshold). 
* The function can raise a Process.All_Children_Dead exception if there are no
* more available processes to continue the computation. *)
val gather_minimum :
 Register.all_parallel -> Register.all_parallel -> 'a array ref -> float -> 'b array ref

(** [gather_maximum f c a th] provides the same functionality of gather_maximum
* but finds the maximum instead of the minimum, that is the set of elements with
* cost bigger than maximum - th. 
* The function can raise a Process.All_Children_Dead exception if there are no
* more available processes to continue the computation. *)
val gather_maximum :
  Register.all_parallel -> Register.all_parallel -> 'a array ref -> float -> 'b array ref

(** [map f a] applies the function f to each element in a and returns the
* resulting array. *)
val map : Register.all_parallel -> 'a array ref -> 'b array ref

(** [scan f a] returns all the members in a that have true value as a result of
* applying f to them. *)
val scan : Register.all_parallel -> 'a array ref -> 'a array ref

(* vim:sw=4 et tw=80
 *)
