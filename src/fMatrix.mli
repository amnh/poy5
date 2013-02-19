(* POY 5.0 Beta. A phylogenetic analysis program using Dynamic Homologies.    *)
(* Copyright (C) 2011  Andrés Varón, Lin Hong, Nicholas Lucaroni, Ward Wheeler*)
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

(** FMatrix; 64-bit Floating Point Temporary Storage Module *)

type m
(** Abstract type to hold scratch space for 64bit floating point values.
    Auto-expands; does not shrink. *)

val scratch_space : m
(** The default instantiation; initial size is '20'. This is very small, but
    will instantly grow to a reasonable size when used first. *)

external print : m -> unit = "floatmatrix_CAML_print"
(** [print m] Print the contents of the matrix; used for debugging. **)

external freeall : m -> unit = "floatmatrix_CAML_freeall"
(** [freeall m i j] Clear the matrix. **)

