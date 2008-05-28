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

(** The implementation of funtions to calculate the cost, 
* alignments and medians between general sequences 
* where both point mutations and rearrangement operations
* are considered *)
val fprintf : out_channel -> ('a, out_channel, unit) format -> 'a

(** Data structure to contain a breakinv sequence *)
type breakinv_t = {
    seq : Sequence.s; (** the breakinv sequence *)

    alied_med : Sequence.s; (** the aligned median sequence  *)
    alied_seq1 : Sequence.s; (** the aligned sequence of its first child *)
    alied_seq2 : Sequence.s; (** the aligned sequence of its second child *)

    ref_code : int; (** the reference code of this breakinv *)
    ref_code1 : int; (** the reference code of its first child *)
    ref_code2 : int; (** the reference code of its second child *)
    cost1 : int; (** the cost between this breakinv and its first child *)
    cost2 : int; (** the cost between this breakinv and its second child *)
    recost1 : int; (** the recost between this breakinv and its first child *)
    recost2 : int; (** the recost between this breakinv and its second child *)
}

(** Data structure to contain parameters 
* used to align two breakinv sequences *)
type breakinvPam_t = {
  re_meth : Data.re_meth_t;
  keep_median : int;
  circular : int;
  swap_med : int;
  symmetric : bool;
}


val breakinvPam_default : breakinvPam_t

(** [init seq] creates a new breakinv sequence from [seq] *)
val init : Sequence.s -> breakinv_t

(** [get_breakinv_pam user_breakinv_pam] returns 
* user defined parameters used to align two breakinv sequences *)
val get_breakinv_pam : Data.dyna_pam_t -> breakinvPam_t

(** [get_recost pams] returns the rearrangement cost in [pams] *)
val get_recost : Data.dyna_pam_t -> int

(** [cmp_cost med1 med2 gen_cost_mat pure_gen_cost_mat alpha breakinv_pam]
* computes total cost between two breakinv sequences [med1] and [med2].
* the total cost = editing cost + rearrangement cost *)
val cmp_cost :
  breakinv_t ->
  breakinv_t ->
  Cost_matrix.Two_D.m ->
  int array array -> Alphabet.a -> Data.dyna_pam_t -> int * (int * int)

(** [find_med2_ls med1 med2 gen_cost_mat pure_gen_cost_mat alpha breakinv_pam]
* finds all median sequences between [med1] and [med2] and vice versa.
* Rearrangements are allowed *) 
val find_med2_ls :
  breakinv_t ->
  breakinv_t ->
  Cost_matrix.Two_D.m ->
  int array array ->
  Alphabet.a -> Data.dyna_pam_t -> int * (int * int) * breakinv_t list

(** [get_costs med child_ref] returns the cost
* from this breakinv median to its child [child_ref] *)
val get_costs : breakinv_t -> int -> int * int
