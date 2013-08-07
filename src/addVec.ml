(* POY 5.0. A phylogenetic analysis program using Dynamic Homologies.         *)
(* Copyright (C) 2013 Andrés Varón, Lin Hong, Nicholas Lucaroni, Ward Wheeler,*)
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

type ct

external register : unit -> unit = "add_CAML_register"
external register_globals : unit -> unit = "add_CAML_register_mem"

let () = register (); register_globals ()

external create : int array -> int array -> ct = "add_CAML_create"

external median : ct -> ct -> ct = "add_CAML_median"

external distance : ct -> ct -> float = "add_CAML_distance"

external distance_2 : ct -> ct -> ct -> float = "add_CAML_distance_2"

external distance_median : ct -> ct -> float * ct = "add_CAML_distance_and_median"

external median_cost : ct -> float = "add_CAML_total"

external compare_data : ct -> ct -> int = "add_CAML_compare_data"

external copy : ct -> ct -> unit = "add_CAML_copy"

external clone : ct -> ct = "add_CAML_dup"

external pos_set_state : ct -> int -> int -> int -> unit = "add_CAML_set_state"

external pos_get_max : ct -> int -> int = "add_CAML_get_max"

external pos_get_min : ct -> int -> int = "add_CAML_get_min"

external pos_get_cost : ct -> int -> float = "add_CAML_get_cost"

external median_3 : ct -> ct -> ct -> ct -> ct = "add_CAML_median_3"

external full_unioni : ct -> ct -> ct -> unit = "add_CAML_full_union"
    
let full_union a b = 
    let r = clone a in
    full_unioni a b r;
    r

external mediani : ct -> ct -> ct -> unit = "add_CAML_median_imp"

external to_string : ct -> string = "add_CAML_to_string"

external cardinal : ct -> int = "add_CAML_cardinal"

let print ct = Printf.printf "%s" (to_string ct)
