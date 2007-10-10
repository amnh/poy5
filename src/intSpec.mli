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

(** Integer specifications for the U language. *)

(** {2 Exceptions} *)

exception Illegal_Element of int

(** {2 Types} *)

type t
type p = Equal | Function of (int -> float)
type func = Constant of float array | Variable of float array

(** {2 Creation} *)

val create : min:int -> max:int -> t
val codes : ints:t -> codes:p -> t
val length : ints:t -> elem:int -> float
val to_list : t -> (int * float) list
val decoder : t -> float
val to_formatter : t -> Tags.output Sexpr.t
val bounds : t -> (int * int)
val get_class : t -> func

