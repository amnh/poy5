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

type characters = 
    | Range of (int * int)
    | Single of int
    | All

type cost_change = (bool * int * (string * string) * int)

type char_change =
    | Additive of characters list
    | NonAdditive of characters list
    | Active of characters list
    | Inactive of characters list
    | Sankoff of characters list
    | Weight of (int * characters list)

type char_name = string
type command = 
    | Nstates of [ `Dna | `Rna | `Proteins | `Number of int ] option
    | Ccode of char_change list
    | Cost of cost_change list
    | Tread of string
    | Xread of string
    | Ignore
    | Charname of char_name list

