(* POY 5.0 Alpha. A phylogenetic analysis program using Dynamic Homologies.   *)
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

module File : sig
    val of_channel : in_channel -> string -> Nexus.File.nexus
end

module P : sig
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


    val is_hennig : FileStream.f -> bool
end

module Grammar : sig
    type token
    val xread : (Lexing.lexbuf -> token) -> Lexing.lexbuf -> (int * int * string)
    val command : (Lexing.lexbuf -> token) -> Lexing.lexbuf -> P.command
end

module Lexer : sig
    exception Eof
    val token : Lexing.lexbuf -> Grammar.token
    val xread : Lexing.lexbuf -> Grammar.token
end
