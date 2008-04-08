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
    between chromosomes where both point mutations and rearrangement operations
    are considered *)

type seed_t = Seed.seed_t
type pairChromPam_t = ChromPam.chromPairAliPam_t
type block_t = Block.block_t
type subseq_t = Subseq.subseq_t
val fprintf : out_channel -> ('a, out_channel, unit) format -> 'a
val deref : 'a option -> 'a
val debug : bool
type direction_t = ChromPam.direction_t
type seg_t = {
  sta : int;
  en : int;
  cost : int;
  alied_med : Sequence.s;
  sta1 : int;
  en1 : int;
  alied_seq1 : Sequence.s;
  dir1 : direction_t;
  sta2 : int;
  en2 : int;
  alied_seq2 : Sequence.s;
  dir2 : direction_t;
}
type med_t = {
  seq : Sequence.s;

  ref_code : int;
  ref_code1 : int;
  ref_code2 : int;
  cost1 : int;
  recost1 : int;
  cost2 : int;
  recost2 : int;
  chrom_map : seg_t list;
}
val create_med : Sequence.s -> med_t
val to_single : med_t -> int -> Cost_matrix.Two_D.m -> Data.dyna_pam_t -> Sequence.s

val init_med : Sequence.s -> med_t

val print_map : seg_t list -> unit
val print_median : med_t list -> string -> unit

val create_map : med_t -> int -> int * int * Tags.output
val create_single_map : med_t -> Tags.output

val create_global_map :
  Sequence.s ->
  Sequence.s ->
  Cost_matrix.Two_D.m ->
  Seed.pairChromPam_t ->
  Block.block_t list * Block.subseq_t list * Block.subseq_t list


val create_median :
  Subseq.subseq_t list ->
  Subseq.subseq_t list ->
    int ->
  Sequence.s * int ->
  Sequence.s * int ->
  Block.block_t list ->
  (Sequence.s * Sequence.s) array array ->
  int array ->
  int array ->
  int array * int * int * int ->
  Cost_matrix.Two_D.m -> ChromPam.chromPairAliPam_t -> med_t


val cmp_cost :
  med_t ->
  med_t ->
  Cost_matrix.Two_D.m ->
  Data.dyna_pam_t -> [< `Chromosome | `Genome ] -> int * int

val find_med2_ls :
  med_t ->
  med_t -> Cost_matrix.Two_D.m -> Data.dyna_pam_t -> int * int * med_t list
val test : unit -> 'a



val find_approx_med2 : med_t -> med_t -> med_t -> med_t 
val print : med_t -> unit
val copy_chrom_map : med_t -> med_t -> med_t
val to_single_root : med_t -> int -> Cost_matrix.Two_D.m -> Sequence.s

val change_to_single : med_t -> Sequence.s -> Cost_matrix.Two_D.m -> med_t


val find_med3 :
  med_t ->
  med_t ->
  med_t ->
  med_t ->
  Cost_matrix.Two_D.m ->
  Cost_matrix.Three_D.m -> Data.dyna_pam_t -> int * med_t
