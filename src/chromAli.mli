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

(** A pair of segements on both chromosomes and the median sequence *)
type seg_t = {
    sta : int; (** start position of this segment on the median sequence *)
    en : int; (** end position of this segment on the median sequence *)
    cost : int; (** the cost of aligning the segment *)
    alied_med : Sequence.s; (** the median sequence of the segment *)

    sta1 : int; (** start position of the segment on the first chromosome *)
    en1 : int; (** end position of the segment on the first chromosome *)
    alied_seq1 : Sequence.s; (** aligned sequence of this segment on the first chromosome *)
    dir1 : direction_t; (** the orientation of this segment on the first chromosome *) 

    sta2 : int; (** start position of the segment on the second chromosome *)
    en2 : int;  (** end position of the segment on the second chromosome *)
    alied_seq2 : Sequence.s; (** alied_seq1 <-> reversed alied_seq2 
                                if dir2 is Negative *)
    dir2 : direction_t; (** the orientation of this segment on the first chromosome *) 
}

(** the median data structure of two chromosomes *)
type med_t = {
    seq : Sequence.s; (** the median sequence *)
    ref_code : int; (** reference code of this median *)
    ref_code1 : int;    (** Child's code *)    
    ref_code2 : int;  (** Child's code *)
    cost1 : int;  
    recost1 : int;
    cost2 : int;
    recost2 : int;
    chrom_map : seg_t list;    
}

(** [create_med seq] return a new median from a seq *)
val create_med : Sequence.s -> med_t


(** [to_single single_parent child_ref c2 pam] returns the
* single state sequence for the node [child_ref] *)
val to_single : med_t -> int -> Cost_matrix.Two_D.m -> Data.dyna_pam_t -> Sequence.s

val init_med : Sequence.s -> med_t

val print_map : seg_t list -> unit
val print_median : med_t list -> string -> unit

(** [create_map anc_med des_ref] creates a map from the 
* ancestor sequence  [anc_med] to descendant sequence whose ref_code is [des_ref]*)
val create_map : med_t -> int -> int * int * Tags.output

val create_single_map : med_t -> Tags.output

(** [create_global_map seq1 seq2 cost_mat ali_pam] creates the 
*  global map between two chromosomes [seq1] and [seq2] 
* Rearrangements are allowed in the global map *)
val create_global_map :
  Sequence.s ->
  Sequence.s ->
  Cost_matrix.Two_D.m ->
  Seed.pairChromPam_t ->
  Block.block_t list * Block.subseq_t list * Block.subseq_t list


(** [create_median subseq1_ls subseq2_ls gen_gap_code (seq1, chrom1_id) (seq2, chrom2_id) global_map 
*                  ali_mat alied_gen_seq1 alied_gen_seq2 (order2_arr, total_cost, recost1, recost2) 
*                   cost_mat ali_pam] creates the median between 
* two chromosome [seq1] and [seq2] given a [global_map] *)
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


(** [cmp_simple_cost med1 med2 cost_mat chrom_pams state] computes
* the cost wthich is the min of cost between [med1] and [med2],  
* and cost between [med2] and [med1]. Rearrangement operations are allowed *)
val cmp_cost :
  med_t ->
  med_t ->
  Cost_matrix.Two_D.m ->
  Data.dyna_pam_t -> [< `Chromosome | `Genome ] -> int * int

(** [find_med2_ls med1 med2 cost_mat user_chrom_pam] find the
* median list whose cost is minimum of cost between [med1] and [med2]
* and between [med2] and [med1] *) 
val find_med2_ls :
  med_t ->
  med_t -> Cost_matrix.Two_D.m -> Data.dyna_pam_t -> int * int * med_t list




(** [find_approx_med2 med1 med2 med12] returns the median 
* between [med1] and [med2] based on all information 
* from [med12] which is a computed median of [med1] and [med2] *)
val find_approx_med2 : med_t -> med_t -> med_t -> med_t 

(** [find_med3 ch1 ch2 ch3 mine c2 c3 pam] returns
* the median of [ch1], [ch2] and [ch3] where [mine]
* is the current median of [ch1], [ch2] and [ch3] *) 
val find_med3 :
  med_t ->
  med_t ->
  med_t ->
  med_t ->
  Cost_matrix.Two_D.m ->
  Cost_matrix.Three_D.m -> Data.dyna_pam_t -> int * med_t

val print : med_t -> unit

(** [copy_chrom_map s d] copy the chromosome map
* of chromosome [s] into chromosome map of chromosome [d] *)
val copy_chrom_map : med_t -> med_t -> med_t

(** [to_single_root root other_code c2] returns a 
* the map of the root for single state *)
val to_single_root : med_t -> int -> Cost_matrix.Two_D.m -> Sequence.s

(** [change_to_single med single_seq c2] 
* changes the single state of the [med] by [single_seq] *)
val change_to_single : med_t -> Sequence.s -> Cost_matrix.Two_D.m -> med_t

