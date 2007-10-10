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

let () = SadmanOutput.register "All_sets" "$Revision: 1728 $"

module OrderedString = struct
    type t = string 
    let compare a b = Pervasives.compare a b
end

module OrderedFloat = struct
    type t = float
    let compare a b = Pervasives.compare a b
end

module OrderedInt = struct
    type t = int
    let compare a b = a - b
end

module OrderedTuple = struct
    type t = (int * int)
    let compare (a, _) (c, _) = a - c
end

module FullOrderedTuple = struct
    type t = (int * int) 
    let compare (a, b) (c, d) =
        match a - c with
        | 0 -> b - d
        | x -> x
end

module OrderedIntList = struct
    type t = int list 
    let rec compare a b =
        match a, b with
        | ha :: ta, hb :: tb -> 
                begin match ha - hb with
                | 0 -> compare ta tb
                | v -> v
                end
        | [], [] -> 0
        | [], _ -> (-1)
        | _, [] -> 1
end

module Integers = Set.Make (OrderedInt)
module IntegerMap = Map.Make (OrderedInt)

module IntegersList = Set.Make (OrderedIntList)
module IntegerListMap = Map.Make (OrderedIntList)

module Floats = Set.Make (OrderedFloat)
module FloatMap = Map.Make (OrderedFloat)

module Strings = Set.Make (OrderedString)
module StringMap = Map.Make (OrderedString)

module TupleMap = Map.Make (OrderedTuple)

module FullTuples = Set.Make (FullOrderedTuple)

module FullTupleMap = Map.Make (FullOrderedTuple)

module OrderedIntSet = struct
    type t = Integers.t
    let compare = Integers.compare
end

module IntSet = Set.Make (OrderedIntSet)

