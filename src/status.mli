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

exception Illegal_update

type formatter_output = StatusCommon.formatter_output
type c = SearchReport | Status | Warning | Error | Information 
         | Output of (string option * bool * formatter_output list)


type status

val main_loop : (string -> unit) -> unit

val type_io : string -> int -> c -> [ 
    | `Status of (string * int)
    | `Warning of (string * int)
    | `Error of (string * int)
    | `Information of (string * int)
    | `Output of (string * int)
]

val set_verbosity : [ `Low | `Medium | `High ] -> unit

val get_verbosity : unit -> [`Low | `Medium | `High ]

val send_output : Pervasives.out_channel list -> unit

val create : string -> int option -> string -> status

val achieved : status -> int -> unit

val get_achieved : status -> int

val message : status -> string -> unit

val report : status -> unit

val finished : status -> unit

val user_message : c -> string -> unit

val full_report : ?msg:string -> ?adv:int -> status -> unit

val is_parallel : (c -> string -> unit) option -> unit

val rank : int -> unit

val init : unit -> unit

val error_location : int -> int -> unit

val output_table : c -> string array array -> unit

val redraw_screen : unit -> unit

val resize_history : int -> unit

val clear_status_subwindows : unit -> unit

val is_interactive : unit -> bool
