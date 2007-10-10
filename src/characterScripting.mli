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

type characters = [
    | `All
    | `Some of (bool * int list)
    | `Names of (bool * string list)
    | `Random of float
    | `AllDynamic
    | `AllStatic
    | `Missing of (bool * int)
]

type taxa = characters

type std_cs =
    | Nonadd8 of NonaddCS8.t           (** A set of non additive characters with
                                           at most 8 states *)
    | Nonadd16 of NonaddCS16.t         (** A set of non additive characters with
                                           at most 16 states *)
    | Nonadd32 of NonaddCS32.t         (** A set of non additive characters with
                                           at most 32 states *)
    | Add of AddCS.t                    (** A set of additive characters *)
    | Sank of SankCS.t                  (** A set of sankoff characters *)
    | Dynamic of DynamicCS.t                    (** A set of sequences *)

    | Set of std_cs list                    (** A set of characters *)

module type S = sig
    type cs
    type n
    (* This type must be consistant with (cs, float) Methods.character_input_output
    * *)
    type character_input_output = [ 
        | `Characters of cs Sexpr.t
        | `Floats of float Sexpr.t
    ]

    val distance : cs -> cs -> float

    val median : cs -> cs -> cs

    val scriptchar_operations :
        n list ->
        Data.d ->
        [ `Distance of (characters * characters) * characters
        | `Median of (characters * characters) * characters ] -> character_input_output

    val filter_char_operations : n list -> Data.d -> (taxa * taxa) ->
        characters -> (cs Sexpr.t * cs Sexpr.t) list

    val extract_character : int -> n -> cs 

    val character_operations : [ `Distance of (cs Sexpr.t * cs Sexpr.t) | `Median of (cs
        Sexpr.t * cs Sexpr.t) ] -> [ `Characters of cs Sexpr.t | `Floats of float Sexpr.t ]
end

module Standard : S with type cs = std_cs with type n = Node.node_data
