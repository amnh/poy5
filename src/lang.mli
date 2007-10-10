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

module Character : sig
    exception Empty_function
    exception Undefined_function of string
    exception Inconsistent_function_parameter
    exception Inconsistent_function_type
    exception Illegal_argument
    exception Illegal_parameter
    (* The only primitive types supported by the U language *)
    type utype = | Base | Seq | Int | Bool | Unkn
    (* A character definition specification statistics, this is for type *)
    type chardef = {
        mutable functions : string list;
        mutable function_param : (string, utype list) Hashtbl.t;
        mutable function_types : (string, utype) Hashtbl.t;
        mutable function_count : (string, int) Hashtbl.t;
    }
    (* Reinitialize the internal state of the parser *)
    val clean : unit -> unit
    (* Initialize the contents of the default functions *)
    val init : unit -> unit
    (* Returns the complete internal state *)
    val get : unit -> chardef
    (* Checks if a function has already been defined *)
    val exists : string -> bool
    (* Adds a function definition to the current state of the parser *)
    val addf : string list -> unit
    (* Gets a function definition elements in the current state of the parser.
    * Raises Undefined_function if the searched function has not been defined. *)
    val getf : string -> utype list * utype * int
    (* Sets a function type in the internal state. Raises Undefined_function if
    * the searched function has not been defined. *)
    val setft : string -> utype -> unit
    (* Sets the type of a parameter of a function. The first parameter is 0.
    * Raises Undefined_function if the searched function has not been defined. *)
    val setpt : string -> int -> utype -> unit
    (* Gets the type of a parameter of a function. The first parameter is 0.
    * Raises Illegal_argument if the parameter index is invalid and
    * Undefined_function if the function searched has not been defined. *)
    val getpt : string -> int -> utype
    (* Increments the internal counter for a particular function. Raises
    * Undefined_function if the function searcher has not been defined yet . *)
    val incr_count : string -> unit
    (* Gets the current counter for a particular function. Raises
    * Undefined_function if the function searched has not been defined yet. *)
    val get_count : string -> int
    val elements : chardef -> string list * (string, utype list) Hashtbl.t *
        (string, utype) Hashtbl.t * (string, int) Hashtbl.t
end
