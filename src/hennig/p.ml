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
type gappy = [ `Nogap | `Gap ] option
type command = 
    | Nstates of [ `Dna of gappy | `Protein of gappy | `Number of int ] option
    | Ccode of char_change list
    | Cost of cost_change list
    | Tread of string
    | Xread of string
    | Ignore
    | Charname of char_name list


let is_hennig file =
    let ch = FileStream.Pervasives.open_in file in
    try 
        while true do
            let line = FileStream.Pervasives.input_line ch in
            if String.length line > 4 then
                let lst = [ line.[0]; line.[1]; line.[2]; line.[3]; line.[4] ] in
                let lst = List.map Char.uppercase lst in
                (* Does the line start with xread or tread? *)
                if lst = ['X';'R';'E';'A';'D'] || lst = ['T';'R';'E';'A';'D'] then 
                    raise Exit;
        done;
        assert false
    with
    | Exit -> true
    | End_of_file -> false

