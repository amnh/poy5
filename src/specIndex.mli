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

(** Index of Specifications for the U language *)

(** {2 Exceptions} *)

exception Not_Defined of string
exception Illegal_Element of string
exception Name_Used
exception Illegal_Combination

(** {2 Types} *)

(** The specification of alphabets, integers and words *)
type s = Alph of AlphSpec.t | Int of IntSpec.t | Word of WordSpec.t

type e = InInt of int | InAlph of string | InWord of int

type t = (string * s) list

val empty : unit -> 'a list
val add : index:t -> name:string -> spec:s -> t
val remove : index:t -> name:string -> t
val length : index:t -> spec:s -> elem:e -> float
val find : index:t -> name:string -> s
val name : index:t -> item:s -> string
val to_list : t -> (string * s) list
val k : t -> float
val to_formatter : t -> Tags.output
