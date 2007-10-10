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

(** Convenient sets and mappings.
*
* Implementations of various instances of modules generated using the standard
* library {!Set.Make} and {!Map.Make} functors, widely used within the Poy
* source code. For further information on each check their respective
* documentation. *)

module OrderedString : Set.OrderedType with type t = string
module OrderedFloat : Set.OrderedType with type t = float
module OrderedInt : Set.OrderedType with type t = int
module OrderedTuple : Set.OrderedType with type t = (int * int)
module FullOrderedTuple : Set.OrderedType with type t = (int * int)

module OrderedIntList : Set.OrderedType with type t = int list

module Integers : Set.S with type elt = int

module Strings : Set.S with type elt = string

module IntegerMap : Map.S with type key = int

module IntegersList : Set.S with type elt = int list

module IntegerListMap : Map.S with type key = int list

module Floats : Set.S with type elt = float

module StringMap : Map.S with type key = string

module TupleMap : Map.S with type key = (int * int)

module FullTuples : Set.S with type elt = (int * int)

module OrderedIntSet : Set.OrderedType with type t = Integers.t

module IntSet : Set.S with type elt = Integers.t

module FullTupleMap :Map.S with type key = (int * int)
