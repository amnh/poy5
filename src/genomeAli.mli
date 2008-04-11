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
val gap : int
type direction_t = ChromPam.direction_t
(** [seg_t] is a data structure to contain a pair of segments 
* on a chromosome of the genome *)
type seg_t = {
    sta : int; (** start position on the median sequence *)
    en : int; (** end position on the median sequence *)
    cost : int; (** the cost segment *)
    med_chrom_id : int; (** median chromosome id *)
    alied_med : Sequence.s; (** aligned median sequence *)


    ref_code1 : int;    (** first child's code *)
    sta1 : int; (** start position on the first child *)
    en1 : int; (** end position on the first child *)
    chi1_chrom_id : int; (** first child chromosome id *)
    alied_seq1 : Sequence.s; (** aligned sequence of the first child *)
    dir1 : direction_t; (** sequence orientation of the first child *)

    ref_code2 : int;    (** second child's code *)
    sta2 : int; (** start position on the second child *)
    en2 : int; (** end position on the second child *)
    chi2_chrom_id : int; (** second child chromosome id *)
    alied_seq2 : Sequence.s; (** aligned sequence of the second child *)
    dir2 : direction_t; (** sequence orientation of the second child *)
}

(** [chrom_t] is a data structure to contain a chromosome median *)
type chrom_t = {
    chrom_id : int ref; (** chromosome id *)
    main_chrom1_id : int; (** main first child's chromosome's id *)
    main_chrom2_id : int; (** main second child's chromosome's id *)

    chrom_ref_code : int; (** chromosome  reference code *)
    map : seg_t list;    
    seq : Sequence.s;   (** the chromosome sequence *)
}

(** [med_t] is data structure to contain a genome*)
type med_t = {
    chrom_arr : chrom_t array; (** chromosome array *)
    genome_ref_code : int; (** genome reference code *)
    genome_ref_code1 : int; (** the first genome child's code *)
    genome_ref_code2 : int; (** the second genome child's code *)
    cost1 :int; (** the cost from this genome to its first child *)
    recost1 : int; (** the recost from this genome to its first child *)
    cost2 : int; (** the cost from this genome to its second child *)
    recost2 : int; (** the recost from this genome to its second child *)

}
type genome_block_t = {
  block_id : int;
  chrom1_id : int;
  chrom2_id : int;
  block : Block.block_t;
}
(** [to_string med] converts all information in 
* median [med] into a string *) 
val to_string : med_t -> string

val print_genome : med_t -> unit
val get_chroms : med_t -> Sequence.s array

(** [chrom_map_to_string chrom_map] converts the information
* in the chromosome map [chrom_ma] into a string *)
val chrom_map_to_string : seg_t list -> string

(** [genome_map_to_string genome] converts all information
* in the genome [genome] into a string *)
val genome_map_to_string : med_t -> string

val ref_genome : med_t ref

(** [assign_hom_chrom med cost_mat user_chrom_pams]
* assigns homologous Id chrom all chromosomes in
* genome [med] *)
val assign_hom_chrom :
  med_t -> Cost_matrix.Two_D.m -> Data.dyna_pam_t -> unit

(** [init genome] returns a new genome created
* from an array of chromosomes [genome] passed from
* Data module *)
val init : Sequence.s Data.dyna_data -> med_t

(** [create_med_from_seq chrom_arr] returns a genome
* created from an array chromosomes [chrom_arr] *)
val create_med_from_seq : Sequence.s array ->med_t

(** [find_conserved_areas seq1 seq2 cost_mat chrom_pams]
* returns a global map containing conserved areas between
* two sequences [seq1] and [seq2] *)
val find_conserved_areas :
  Sequence.s ->
  Sequence.s -> Cost_matrix.Two_D.m -> Data.dyna_pam_t -> Block.block_t list

(** [find_chrom chrom_id genome] returns the chromosome [chrom_id]
* in the genome [genome] *)
val find_chrom : int -> med_t -> chrom_t option

val create_loci : Sequence.s -> Subseq.subseq_t list -> Subseq.subseq_t list
val create_fast_general_ali :
  int ->
  int ->
  Sequence.s ->
  Subseq.subseq_t list ->
  int ->
  Sequence.s ->
  Subseq.subseq_t list ->
  genome_block_t list ->
  Cost_matrix.Two_D.m ->
  ChromPam.chromPairAliPam_t -> Sequence.s * seg_t list * int * int * int
val create_chrom_med :
  int * chrom_t ->
  int * chrom_t ->
  genome_block_t list ->
  Cost_matrix.Two_D.m ->
  ChromPam.chromPairAliPam_t -> chrom_t * int * int * int
val create_genome_blocks :
  med_t ->
  med_t -> Cost_matrix.Two_D.m -> Data.dyna_pam_t -> genome_block_t list
val create_med :
  med_t ->
  med_t -> Cost_matrix.Two_D.m -> Data.dyna_pam_t -> med_t * int * (int * int)
val cmp_cost :
  med_t -> med_t -> Cost_matrix.Two_D.m -> Data.dyna_pam_t -> int * (int * int)
val find_med2_ls :
  med_t ->
  med_t ->
  Cost_matrix.Two_D.m -> Data.dyna_pam_t -> int * (int * int) * med_t list
val compare : med_t -> med_t -> int
val create_map : med_t -> int -> int * int * Tags.output
val to_single :
  med_t ->
  med_t ->
  Cost_matrix.Two_D.m -> Data.dyna_pam_t -> Sequence.s array
val change_to_single : med_t -> Sequence.s array -> med_t
val copy_chrom_map : med_t -> med_t -> med_t
val to_single_root : med_t -> Cost_matrix.Two_D.m -> Sequence.s array
val create_single_map : med_t -> Tags.output

val find_med3 :
  med_t ->
  med_t ->
  med_t ->
  med_t ->
  Cost_matrix.Two_D.m ->
  Cost_matrix.Three_D.m -> Data.dyna_pam_t -> int * med_t



