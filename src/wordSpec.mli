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

(** Word or Sequence specification for the U language. *)

(** {2 Exceptions } *)

exception Illegal_Element of int

(** {2 Types } *)

(** The word spec type *)
type t

(** The probability of each word length *)
type p = Equal | Function of (int -> float)

(** {2 Functions } *)

(* [create min max alph] creates a word specification with length between [min]
 * and [max] over the alphaabet [alph]. *)
val create : min:int -> max:int -> alph:AlphSpec.t -> t

(** [codes word codes] returns a fresh copy of [word] with probability [codes].
* *)
val codes : word:t -> codes:p -> t

(** [length word length] returns the length in nats of a word of length [length]
 * according to the specs [word]. *)
val length : word:t -> length:int -> float

(** [wordlength word elem] returns the length of the word [elem] according to
* the word specs [word]. *)
val wordlength : word:t -> elem:string -> float

(** [to_list t] returns an association list with the probability of each word
* length in the word specifications [t]. *)
val to_list : t -> (int * float) list

val to_formatter : t -> Tags.output Sexpr.t
