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
* alignments and medians between annotated chromosomes 
* where both point mutations and rearrangement operations 
* are considered *)

val fprintf : out_channel -> ('a, out_channel, unit) format -> 'a

(** [seq_t] is data structure to contain a segment of an annotated chromosome *)
type seq_t = {
    seq : Sequence.s; (** the segment sequence *)
    seq_ref_code : int; (** the reference code of this segment *)
    alied_med : Sequence.s; (** the aligned sequence this segement *)

    seq_ord1 : int; (** the segment order of its first child *)
    alied_seq1 : Sequence.s; (** the aligned sequence of its first child *)

    seq_ord2 : int; (** the segment order of its second child *)
    alied_seq2 : Sequence.s; (** the aligned sequence of its second child *)
}

(** [annchrom_t] is a data structure to contain an annoated chromosome *)
type annchrom_t = {
    seq_arr : seq_t array;  (** annotated chromosome *)
    ref_code : int;  (** reference code of this chromosome *)
    ref_code1 : int;    (** first child's code *)    
    ref_code2 : int;  (** second child's code *)
    cost1 : int; (** cost from this chromosome to its first child *)
    recost1 : int; (** recost from this chromosome to its first child *)
    cost2 : int; (** cost from this chromosome to its second child *)
    recost2 : int; (** recost from this chromosome to its second child *) 
}

(** [annchromPam_t] is data structure to 
* contains parameters used to align two annotated chromosomes *)
type annchromPam_t = {
  re_meth : Data.re_meth_t;
  keep_median : int;
  circular : int;
  swap_med : int;
  approx : ChromPam.order_t;
  symmetric : bool;
  locus_indel_cost : int * int;
}

(** [clone_seq] returns a fresh clone of segment [s] *)
val clone_seq : seq_t -> seq_t

(** [clone_med] returns a fresh clone of annotated chromosome [m] *)
val clone_med : annchrom_t -> annchrom_t

val annchromPam_default : annchromPam_t

(** [init_seq_t (seq, code)] returns a segment from [seq] and [code]*)
val init_seq_t : Sequence.s * int -> seq_t

(** [init seq_arr] returns an annotated chromosome 
* from sequence array [seq_arr] *)
val init : (Sequence.s * int) array -> annchrom_t

val printMap : seq_t array -> unit
val get_seq_arr : annchrom_t -> Sequence.s array

(** [get_annchrom_pam] returns user defined parameters
* to align two annotated chromosomes *)
val get_annchrom_pam : Data.dyna_pam_t -> annchromPam_t

val print : annchrom_t -> Alphabet.a -> unit

(** [split chrom] returns an array code sequences and an array
* of codes coresponding to the sequence array*)
val split : annchrom_t -> Sequence.s array * int array

(** [create_pure_gen_cost_mat seq1_arr seq2_arr cost_mat ali_pam]
* return a cost matrix between segments of [seq1_arr] and segments of [seq2_arr] *)
val create_pure_gen_cost_mat :
  Sequence.s array ->
  Sequence.s array ->
  Cost_matrix.Two_D.m ->
  annchromPam_t -> int array array * int array * int array * int

(** [cmp_cost chrom1 chrom2 cost_mat alpha annchrom_pam] 
* compute the minimum cost from annotated chromosome [chrom1]
* to annotated chromosome [chrom2], and from [chrom2] to [chrom1] *)
val cmp_cost :
  annchrom_t ->
  annchrom_t -> Cost_matrix.Two_D.m -> 'a -> Data.dyna_pam_t -> int * int

(** [find_med2_ls chrom1 chrom2 cost_mat alpha annchrom_pam]
* finds medians between [chrom1] and [chrom2] and vice versa
* allowing rearrangements *)
val find_med2_ls :
  annchrom_t ->
  annchrom_t ->
  Cost_matrix.Two_D.m -> 'a -> Data.dyna_pam_t -> int * int * annchrom_t list

(** [find_med3 ch1 ch2 ch3 mine c2 c3 alpha annchrom_pam]
* finds the median of annotated chromosome [ch1], [ch2],
* and [ch3] based on the current median [mine] *) 
val find_med3 :
  annchrom_t ->
  annchrom_t ->
  annchrom_t ->
  annchrom_t ->
  Cost_matrix.Two_D.m ->
  Cost_matrix.Three_D.m -> 'a -> Data.dyna_pam_t -> int * annchrom_t

(** [compare annchrom1 annchrom2] compares
* annotated chromosome [annchrom1] and annotated chromosome [annchrom2] *)
val compare : annchrom_t -> annchrom_t -> int

(** [find_approx_med2 med1 med2 med12] returns an approximate
* median between annotated chromosomes [med1] and [med2] based 
* on computed median [med12] *)
val find_approx_med2 : annchrom_t -> annchrom_t -> annchrom_t -> annchrom_t

(** [assign_seq_ref annchrom seq_ref_code] assigns
* ID for all segments of annotated chromosome [annchrom] *)
val assign_seq_ref : annchrom_t -> int -> annchrom_t * int

(** [create_map med child_ref] returns the map
* between chromosome [med] and its child [child_ref] *)
val create_map : annchrom_t -> int -> int * int * Tags.output

val create_single_map : annchrom_t -> Tags.output


(** [to_single single_parent child_ref c2] returns
* thesingle state sequence for annotated chromosome [child_ref] *)
val to_single : annchrom_t -> int -> Cost_matrix.Two_D.m -> Sequence.s array

(** [to_single_root root other_code c2] returns the single
* state sequence for the [root] *)
val to_single_root :
  annchrom_t -> int -> Cost_matrix.Two_D.m -> Sequence.s array

(** [change_to_single med single_seq_arr c2] changes the
* single state sequence of [med] by [single_seq] *)
val change_to_single : annchrom_t -> Sequence.s array -> Cost_matrix.Two_D.m -> annchrom_t
val to_formater : annchrom_t -> Alphabet.a -> string
val copy_chrom_map : annchrom_t -> annchrom_t -> annchrom_t

(** [convert_map med] converts the map of the 
* median [med] into the format suitable for
* creating implied alignments *)
val convert_map : annchrom_t -> (int * int * int * int * int * int) list