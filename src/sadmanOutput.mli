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

(* $Id: sadmanOutput.mli 1800 2007-05-08 01:02:29Z andres $ *)
(* Created Mon Jan 23 11:54:28 2006 (Illya Bomash) *)

(** sadmanOutput.mli *)

val establish_connection : ?port:int -> unit -> unit
val output : string -> unit
val do_sadman : bool ref

val start : string -> (string * string) list -> unit
val finish : (string * string) list -> unit


val register : string -> string -> unit
val get_module_versions : unit -> (string * string) list
val runtime_properties : (string * string) list
