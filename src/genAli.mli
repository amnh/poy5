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

(** This module implements methods to align two general
* characters allowing rearrangements *)

val fprintf : out_channel -> ('a, out_channel, unit) format -> 'a

val to_ori_code : int -> int

val from_ori_code : int -> int

val from_ori_arr : int array -> int array

val to_ori_arr : int array -> int array

(* [cmp_recost_simple seq1 seq2 ali_pam.re_meth ali_pam.circular ]
 * will return the breakpoint or inversion distance of seq1 and seq2.
 * It's like function [cmp_recost state seq1 seq2 reseq2 ...] when reseq2=seq2 
* *)
val cmp_recost_simple : int array -> int array ->
  [< `Locus_Breakpoint of int | `Locus_Inversion of int ] -> int -> int


(** [cmp_recost state seq1 seq2 reseq2 re_meth circular] returns
* the rearrangement distance between two sequence [seq1] and [seq2] *)

val cmp_recost :
  [> `Annotated | `Breakinv | `Chromosome | `Genome ] ->
  int array ->
  int array ->
  int array ->
  [< `Locus_Breakpoint of int | `Locus_Inversion of int ] ->
  int -> bool -> int * int



(** [cmp_cost state code1_arr code2_arr recode2_arr 
*              cost_mat gap re_meth circular] returns
* the total cost between [seq1] and [reseq2]. Precisely,
* total cost = editing cost ([seq1], [reseq2]) + rearrangement cost ([seq2], [reseq2]) *)
val cmp_cost :
  [> `Annotated | `Breakinv | `Chromosome | `Genome ] ->
  int array ->
  int array ->
  int array ->
    Cost_matrix.Two_D.m ->
  int ->
  [< `Locus_Breakpoint of int | `Locus_Inversion of int ] ->
  int -> bool -> bool -> int * (int * int) * int array * int array

(**[cmp_cost3 seq1 seq2 seq3 med cost_mat gap re_meth cir sym] returns
* the total cost between [med] and three sequences [seq1], [seq2], [seq3] *)
val cmp_cost3 :
  int ->   
  int array ->
  int array ->
  int array ->
  int array ->
  int array array ->
    int -> [< `Locus_Breakpoint of int | `Locus_Inversion of int ] -> int ->
        bool -> bool -> (int * int * int)


(** [find_wagner_ali state seq1 seq2 gen_cost_mat gap re_meth circular]
 * returns rearranged sequence [reseq2] of sequence [seq2] using stepwise addition method 
 * such that the total cost is minimum where 
 * total cost = editing cost ([seq1], [reseq2]) + rearrangement cost ([seq2], [reseq2]) *)
val find_wagner_ali :
  int ->
  [> `Annotated | `Breakinv | `Chromosome | `Genome ] ->
  int array ->
  int array ->
    Cost_matrix.Two_D.m ->
  int -> [< `Locus_Breakpoint of int | `Locus_Inversion of int ] -> int -> bool
  -> bool -> int array

(** [multi_swap_locus state seq1 seq2 best_seq2 best_cost 
*                     gen_cost_mat gap re_meth max_swap_med circular num_done_swap] 
* swaps [reseq2] in order to minimize the total cost between [seq1] and [reseq2] where 
* total cost = editing cost ([seq1], [reseq2]) + rearrangement cost ([seq2], [reseq2]) *)
val multi_swap_locus :
  [> `Annotated | `Breakinv | `Chromosome | `Genome ] ->
  int array ->
  int array ->
  int array ->
  int ->
    Cost_matrix.Two_D.m ->
  int ->
  [< `Locus_Breakpoint of int | `Locus_Inversion of int ] ->
  int -> int -> bool -> int -> bool -> int * int array


(** [create_gen_ali state seq1 seq1 gen_cost_mat alpha re_meth max_swap_med circular]
* creates the general alignment between [seq1] and [seq2] with minimum total cost 
* where total cost = editing cost + rearrangement cost *)
val create_gen_ali :
  int ->
  [> `Annotated | `Breakinv | `Chromosome | `Genome ] ->
  Sequence.s ->
  Sequence.s ->
  int array array ->
  Alphabet.a ->
  [< `Locus_Breakpoint of int | `Locus_Inversion of int ] ->
  int -> int -> bool -> int * (int * int) * Sequence.s * Sequence.s

val create_gen_ali_new :  
    int array ->
    int array ->
    int array array -> 
    int ->
    [< `Locus_Breakpoint of int | `Locus_Inversion of int ] ->
    int -> 
    bool ->
    int * int * int array * int array


(** [create_gen_ali_code state seq1 seq2 gen_cost_mat gen_gap_code 
*        re_meth max_swap_med circular] creates the general 
* alignment between [seq1] and [seq2] with minimum total cost
* where total cost = editing cost + rearrangement cost *)
val create_gen_ali_code :
  int ->
  [> `Annotated | `Breakinv | `Chromosome | `Genome ] ->
  int array ->
  int array ->
  int array array ->
  int ->
  [< `Locus_Breakpoint of int | `Locus_Inversion of int ] ->
  int -> int -> bool -> bool -> int * (int * int) * int array * int array


(** [create_gen_ali3 seq1 seq2 seq3 med gen_cost_mat 
*     alpha re_meth  max_swap_med circular sym] create
* the general alignment among [seq1], [seq2], and [seq3] 
* such that total cost = editing cost + rearrangement cost is minimized *)
val create_gen_ali3 :
    int ->
    Sequence.s ->
    Sequence.s ->
    Sequence.s ->
    Sequence.s ->
    int array array ->
    Alphabet.a ->
    [< `Locus_Breakpoint of int | `Locus_Inversion of int ] ->
    'a -> int -> bool -> bool -> Sequence.s * int * int * int

val create_gen_ali3_by_medsov_seqarr :
    int ->
    int ->
    Sequence.s ->
    Sequence.s ->
    Sequence.s ->
    int list list ->
    Cost_matrix.Two_D.m ->
    Alphabet.a ->
    [< `Locus_Breakpoint of int | `Locus_Inversion of int ] ->
    'a -> int -> bool -> bool -> 
        Sequence.s * int list * int * int * int * int * int * int 
        * int array * int array * int array * int array * int array * int array

(*
val create_gen_ali3_by_medsov_codearr :
    Data.median_solver_t ->
    int ->
    int array ->
    int array ->
    int array ->
    int list list ->
    Cost_matrix.Two_D.m ->
    Alphabet.a ->
    [< `Locus_Breakpoint of int | `Locus_Inversion of int ] ->
    'a -> int -> bool -> bool -> 
        Sequence.s * int list * int * int * int * int * int * int 
        * int array * int array * int array * int array * int array * int array

*)
