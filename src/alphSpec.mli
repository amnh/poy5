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

(** Alphabet specification in the U language. 
*
* During the specification of a set of operations over sequential characters,
* the alphabet of those sequences and the elements probabilities must be
* assigned. This library handles such specifications once parsed by
* {!PoyParser.of_channel} and {!PoyParser.process_tree}. *)

(** {2 Exceptions} *)

exception Illegal_Element of string

(** {2 Types} *)

(** An alphabet specification *)
type t 

(** A probability specification of an element in an alphabet. *)
type p =
    | Equal     (** All the elements in a set have the same probability *)
    | Enum of (string * float) list   (** For each pair [(a, b)], the element
    [a] of the alphabet has probability [0.0 <= b <= 1.0] *)

(** {2 Conversion} *)

(** [create lst] creates a fresh alphabet using as elements [lst], and assigns
* the same probability to every element in [lst]. *)
val create : elem:string list -> t

(** [to_list a] converts the alphabet specification [a] to an association list
* containing the elements in the alphabet and their respective probabilities. *)
val to_list : t -> (string * float) list

(** {2 Complexity Related Functions} *)

(** [probs t p] assigns the probability specified by [p] to the elements in [t].
* The sum of the probabilities in [p] must be less than or equal to 1.0. If not
* all the elements in the alphabet are included in [p], the missing elements are
* assigned equal probability so that the sum of [p] and theirs is [1.0]. 
* If one of the elements in [p] does not belong to the alphabet [t] an
* [Illegal_Element] exception is raised. *)
val probs : alph:t -> probs:p -> t

(** [length a i] returns the length in nats of the element [i] from the alphabet
* [t]. If [i] is not part of the alphabet, an [Illegal_Element] exception is
* raised. *)
val length : alph:t -> elem:string -> float

(** [decoder a] calculates the [a] decoder length in nats. *)
val decoder : t -> float

val to_formatter : t -> Tags.output Sexpr.t
