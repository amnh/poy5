(* POY 5.1.1. A phylogenetic analysis program using Dynamic Homologies.       *)
(* Copyright (C) 2014 Andrés Varón, Lin Hong, Nicholas Lucaroni, Ward Wheeler,*)
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

let () = SadmanOutput.register "Matrix" "$Revision: 3649 $"

type m
external create : unit -> m = "mat_CAML_create_general"
external print_2d : m -> int -> int -> unit = "mat_CAML_print_algn_2d"
external print_3d : m -> int -> int -> int -> unit = "mat_CAML_print_algn_3d"
external get_dir : m -> int -> int -> int -> int -> int = "mat_CAML_get_value"
external register : unit -> unit = "mat_CAML_initialize"
external c_flush : m -> unit = "mat_CAML_flush_memory"
external clear_direction : m -> unit = "mat_CAML_clear_direction"

let () = register ()

let default = create ()

let flush () = c_flush default

