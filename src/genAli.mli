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

val fprintf : out_channel -> ('a, out_channel, unit) format -> 'a
type dyna_state_t = Data.dyna_state_t
val cmp_recost :
  [> `Breakinv ] ->
  int array ->
  int array ->
  int array ->
  [< `Breakpoint of int | `Inversion of int ] -> int -> int * int
val cmp_cost :
  [> `Breakinv ] ->
  int array ->
  int array ->
  int array ->
  int array array ->
  int ->
  [< `Breakpoint of int | `Inversion of int ] ->
  int -> int * (int * int) * int array * int array
val find_wagner_ali :
  [> `Breakinv ] ->
  int array ->
  int array ->
  int array array ->
  int -> [< `Breakpoint of int | `Inversion of int ] -> int -> int array
val multi_swap_locus :
  [> `Breakinv ] ->
  int array ->
  int array ->
  int array ->
  int ->
  int array array ->
  int ->
  [< `Breakpoint of int | `Inversion of int ] ->
  int -> int -> int -> int * int array
val create_gen_ali :
  [> `Breakinv ] ->
  Sequence.s ->
  Sequence.s ->
  Cost_matrix.Two_D.m ->
  int array array ->
  Alphabet.a ->
  [< `Breakpoint of int | `Inversion of int ] ->
  int -> int -> int * (int * int) * Sequence.s * Sequence.s
val create_gen_ali_code :
  [> `Breakinv ] ->
  int array ->
  int array ->
  int array array ->
  int ->
  [< `Breakpoint of int | `Inversion of int ] ->
  int -> int -> int * (int * int) * int array * int array

val create_gen_ali3 :
  Sequence.s ->
  Sequence.s ->
  Sequence.s ->
  int array array ->
  Alphabet.a ->
  [< `Breakpoint of int | `Inversion of int ] ->
  'a -> int -> Sequence.s * int
