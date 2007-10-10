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



(** The implementation of funtions to calculate the cost, alignments and medians
    between general sequences where both point mutations and rearrangement operations
    are considered *)


val fprintf : out_channel -> ('a, out_channel, unit) format -> 'a
type breakinv_t = {
    seq : Sequence.s;
    alied_med : Sequence.s;
    alied_seq1 : Sequence.s;
    alied_seq2 : Sequence.s;
    ref_code : int;
    ref_code1 : int;
    ref_code2 : int;
    cost1 : int;
    cost2 : int;
    recost1 : int;
    recost2 : int;
  }

type breakinvPam_t = {
  re_meth : Data.re_meth_t;
  keep_median : int;
  circular : int;
  swap_med : int;
  symmetric : bool;
}


val breakinvPam_default : breakinvPam_t
val init : Sequence.s -> breakinv_t
val get_breakinv_pam : Data.dyna_pam_t -> breakinvPam_t
val cmp_cost :
  breakinv_t ->
  breakinv_t ->
  Cost_matrix.Two_D.m ->
  int array array -> Alphabet.a -> Data.dyna_pam_t -> int * (int * int)
val find_med2_ls :
  breakinv_t ->
  breakinv_t ->
  Cost_matrix.Two_D.m ->
  int array array ->
  Alphabet.a -> Data.dyna_pam_t -> int * (int * int) * breakinv_t list
val get_costs : breakinv_t -> int -> int * int
