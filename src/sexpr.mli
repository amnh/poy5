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

(** Structured Expressions
*
* This module specifies a lightweight version of a tree that can be used to
* express sets of data in tree like form. The module uses polymorphic variant
* types to eliminate the possibility of circular dependencies, and simplify the
* use of more complex structured types, like those defined in the [Tags] module.
* *)


(** {2 Types} *)

(** The structured type, holding contents of type ['a]. *)
type 'a t = [ `Empty | `Set of 'a t list | `Single of 'a ]


(** {2 Functions} *)

(** [fold_left f init s] applies the function [f] on the accumulator and each 
* leaf of [s], outputting the accumulator for the next leaf. The initial value
* of the accumulator is [init], that is, [fold_left f init `Empty] outputs
* [init]. This function is applied in a right to left order, and so the stack
* size is given by the maximum depth of [s]. *)
val fold_left :
  ('a -> 'b -> 'a) -> 'a -> 'b t -> 'a

(** [fold_right f init s] is similar to [fold_left] but the order of aplication
* is rightwise. The maximum stack size is given by the maximum between the depth
* of [s] and the longest list inside a set in [s]. *)
val fold_right :
  ('a -> 'b -> 'a) -> 'a -> 'b t -> 'a

(** [rev s] outputs [s] in reverse order, that is, each list inside [s] is
* reversed. *)
val rev: 'a t -> 'a t

(** [map f s] applies the function [f] on each leaf of [s] producing the
* corresponding structured expression. *)
val map :
  ('a -> 'b) -> 'a t -> 'b t


(** [cardinal s] outputs the number of `Single in [s]. *)
val cardinal : 'a t -> int

(** [map_feedback feed f s] is equivalent to [map] but applies the function
* [feed] with an argument which is a counter of the leaf number where [f] is
* being applied. [feed] is typically used to generate user interface feedback to
* represent the advancement the map application. *)
val map_feedback :
  (int -> unit) -> ('a -> 'b) -> 'a t -> 'b t

(** [map_status message eta f s] is equivalent to [map f s] but outpus a user
* [message] in the POY status box, which would include the time estimation to
* finish the map if [eta] is [true]. *)
val map_status :
    string -> ?eta:bool -> ('a -> 'b) -> 'a t -> 'b t

(** [fold_status message eta f init s]  is equivalent to [fold_left f init s]
 * but outputs a user [message] in the POY status bos, which would include the
 * time estimation to finish the fold oepration if [eta] is [true]. *)
val fold_status :
    string -> ?eta:bool -> ('a -> 'b -> 'a) -> 'a -> 'b t -> 'a

val compose_status :
    string -> ?eta:bool -> (int -> 'a -> 'a) -> int -> 'a -> 'a

(** [singleton x] produce [`Single x] *)
val singleton : 'a -> 'a t

(** [of_list lst] produces a set containing as leaves all the elements in [lst].
* The funtion is not tail recursive. *)
val of_list : 'a list -> 'a t

(** [to_list s] produces a list containing the elements in [s] in the same order
* of appearance in [s]. *)
val to_list : 'a t -> 'a list

(** [first s] outputs the first leaf that would be visited in a left to right
* traversals of [s]. *)
val first : 'a t -> 'a

(** [all_to_all f s t] applies the function [f] on every pair of leaves in [s]
 * and [t], to produce a new structured type holding them all. The overall size
 * of the output is [(cardinal s) * (cardinal t)]. *)
val all_to_all :
  ('a -> 'b -> 'c) -> 'a t -> 'b t -> 'c t

(** [flatten s] takes a structured type [s] holding structured types, and merges
* them in a single structured type. Its equivalent to [List.flatten]. *)
val flatten : 'a t t -> 'a t

(** [iter f s] iteratively applies the function [f] on each element of [s],
* including interior vertices. [leaf_iter] applies the function on each leaf
* only. *)
val iter : ('a t -> unit) -> 'a t -> unit

val shallow_all_to_all : ('a t -> 'a t -> 'b) -> 'a t -> 'a t -> 'b t

(** [leaf_iter f s] applies the function [f] on each leaf of [s]. *)
val leaf_iter : ('a -> unit) -> 'a t -> unit

val to_array : 'a t -> 'a array t

(** Equivalent to [cardinal]. *)
val length : 'a t -> int

(** [nth it s] retrieves the [it]'th leaf that would be visited using
* [fold_left] in [s], starting with [0]. If [s] has less elements than [it],
* then the function raises [Failure "Sexpr.nth"]. *)
val nth : int -> 'a t -> 'a

(** [filter f s] produces a structured type holding all those leaves [l] in [s]
 * where [f l] is [true]. *)
val filter : ('b -> bool) -> 'b t -> 'b t

(** [split f s] returns a tuple of structured types, where the first element in
 * the tuple holds those leaves * [l] in [s] where [f l = true], and the second
 * contains those where [f l = false]. *)
val split : ('b -> bool) -> 'b t -> ('b t * 'b t)

val combine : 'a t * 'b t -> ('a * 'b t) t

val full_combine : 'a t * 'b t -> ('a * 'b) t

val map_insexpr : ('a -> ('b t) list) -> 'a t -> 'b t

val choose_random : 'a t -> 'a option

(** [union s t] simply produces the union of [s] and [t]. *)
val union : 'a t -> 'a t -> 'a t
