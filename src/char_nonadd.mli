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

(** type of a single nonadditive character *)
type t
(** type to generate new nonadditive characters *)
type gen
(** make a character generator *)
val rand_gen : unit -> gen
(** make random nonadditive characters given a generator *)
val make_rand : gen -> t
(** [make_char code val] makes a nonadd char with given code and initial value *)
val make_char : int -> int -> t
(** get the code of an element *)
val code : t -> int
(** take the median between two characters using an optional previous median *)
val median : t option -> t -> t -> t
(** distance between two characters *)
val distance : t -> t -> float
(** convert an element to a string *)
val to_string : t -> string
(** compare function for codes *)
val compare_codes : t -> t -> int
(** compare function for codes and data *)
val compare_data : t -> t -> int
