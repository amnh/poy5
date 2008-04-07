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

(** Blocks are conserved areas between two chromosomes
* which do not require identical nucleotide segments but 
* highly similar. Blocks are considered as homologus segments and used
* as milestones to divide chromosomes into sequences of loci*)

type pairChromPam_t = ChromPam.chromPairAliPam_t
type seed_t = Seed.seed_t
type direction_t = ChromPam.direction_t
type order_t = ChromPam.order_t
type subseq_t = Subseq.subseq_t
val deref : 'a option -> 'a
val fprintf : out_channel -> ('a, out_channel, unit) format -> 'a
module IntSet :
  sig
    type elt = int
    type t = All_sets.Integers.t
    val empty : t
    val is_empty : t -> bool
    val mem : elt -> t -> bool
    val add : elt -> t -> t
    val singleton : elt -> t
    val remove : elt -> t -> t
    val union : t -> t -> t
    val inter : t -> t -> t
    val diff : t -> t -> t
    val compare : t -> t -> int
    val equal : t -> t -> bool
    val subset : t -> t -> bool
    val iter : (elt -> unit) -> t -> unit
    val fold : (elt -> 'a -> 'a) -> t -> 'a -> 'a
    val for_all : (elt -> bool) -> t -> bool
    val exists : (elt -> bool) -> t -> bool
    val filter : (elt -> bool) -> t -> t
    val partition : (elt -> bool) -> t -> t * t
    val cardinal : t -> int
    val elements : t -> elt list
    val min_elt : t -> elt
    val max_elt : t -> elt
    val choose : t -> elt
    val split : elt -> t -> t * bool * t
  end

(** Parameters are used to create blocks between two chromosomes*)
type blockPam_t = {
    (** Two consecutive blocks are connected if their distance and shift are
        smaller than thresholds *)
    max_connecting_dis : int;
    max_connecting_shift: int;
}   
 
(** A block is created by connecting a list of seeds together*)
type block_t = {
    mutable id : int; (** The block id *)
    mutable is_dum : bool; (* Dummy blocks are used as boundary *)

    mutable sta1 : int; (* sta1 is the start of block in the first chromosome *)
    mutable sta2 : int; (* sta2 is the start of block in the second chromosome *)
    mutable en1 : int; (* end1 is the end of block in the first chromosome *)
    mutable en2 : int; (* end2 is the end of block in the second chromosome *)
    mutable direction : direction_t; (* The direction of this block, either postive or negative *)
    mutable cost : int; (* The alignment cost of this block *)
    mutable alied_seq1 : Sequence.s option; (* alied_seq1 and alied_seq2 are aligned sequences of this block *)
    mutable alied_seq2 : Sequence.s option;

    mutable seed_ls : seed_t list; (* The list of seeds constituted this block *)

    (** A chromosome is divided into consecutive sub-sequences *)
    mutable subseq1_id : int; (* The identification of this block in the first chromosome *)
    mutable subseq2_id : int; (* The identification of this block in the second chromosome *)
}
        
val blockPam_default : blockPam_t
val cloneBlockPam : blockPam_t -> blockPam_t

(** [create_from_seed] returns a block 
* with ID [block_id], and contains only one seed [seed] *)
val create_from_seed : int -> seed_t -> block_t

val create_simple_block : int -> int -> int -> int -> int -> block_t

(** [get_dum_first_block ali_pam]  returns  a dummy block which is used 
* as a start point for dynamic programming to connect blocks together *)
val get_dum_first_block : ChromPam.chromPairAliPam_t -> block_t

(** [get_dum_last_block ali_pam] returns  a dummy block which is used 
* as an end point for dynamic programming to connect blocks together *)
val get_dum_last_block : ChromPam.chromPairAliPam_t -> block_t

val max_len : block_t -> int

(** [invert block min_pos2  max_pos2] returns [block']
* which is inverted from [block] due to the inversion 
* from [min_pos2] to [max_pos2] in the second chromosome *)
val invert : block_t -> int -> int -> unit

val get_pos : block_t -> order_t -> int * int
val cmp_dia_dis : block_t -> block_t -> int
val print : block_t -> unit
val add_seed : block_t -> seed_t -> unit

(** [cmp_cost_based_seed block ali_pam] returns an integer
* number as the cost of this [block]. The cost is calculated from 
* its seed list and chromosome parameters [ali_pam] *)
val cmp_cost_based_seed : block_t -> pairChromPam_t -> int

(** [create_from_seed_ls block_id seed_ls ali_pam] returns
* a new block whose ID is [block_id]. The new block is
* created from the [seed_ls] *)
val create_from_seed_ls : int -> seed_t list -> pairChromPam_t -> block_t

(** [cmp_ali_cost alied_seq1 alied_seq2 direction ali_pam] returns
* an integer number as the alignment cost between [alied_seq1]
* and [alied_seq2] based on chromosome parameters [ali_pam] *)
val cmp_ali_cost :
  Sequence.s -> Sequence.s -> direction_t -> pairChromPam_t -> int

(** [find_local_block seed_ls ali_pam] returns a list of blocks 
* created by connecting seeds which are near each other *)
val find_local_block : seed_t list -> pairChromPam_t -> block_t list

(** [is_free b sig1_arr sig2_arr] where
*    - b: a block 
*    - sig1_arr: sig1_arr[i] = -1 if position i 
*                 in the first chromosome is not yet
*                 occupied by any block, otherwise occupied
*    - sig2_arr: sig2_arr[i] = -1 if position i 
*                 in the second chromosome is not yet
*                 occupied by any block, otherwise occupied
*    returns true if whole block [b] is not occupied, otherwise false *)
val is_free : block_t -> int array -> int array -> bool

(** [assign block sig1_arr sig2_ar] marks all positions
* of [block] as occupied position in both genomes *)
val assign : block_t -> int array -> int array -> unit

(** [create_sig_arr block_ls ali_pam] returns
* two signiture arrays [sig1_arr] and [sig2_arr] where
* - sig1_arr[i] is the block ID covering position i in the first chromosome
* - sig2_arr[i] is the block ID covering position i in the second chromosome *)
val create_sig_arr :
  block_t list -> ChromPam.chromPairAliPam_t -> int array * int array

(** [select_separated_block b_ls ali_pam] returns a list
* of blocks which are not overlaped each other. Blocks
* are selected according to their scores. Thus, higher score
* blocks are given higher selection priority *)
val select_separated_block :
  block_t list -> ChromPam.chromPairAliPam_t -> block_t list

(** [create_pos_alied_block block seq1 seq2 cost_mat ali_pam]
* creates the alignment for [block] based on its seed_ls *)
val create_pos_alied_block :
  block_t ->
  Sequence.s -> Sequence.s -> Cost_matrix.Two_D.m -> pairChromPam_t -> unit

(** [is_inide seed block] returns true if [seed] 
* is in side [block], otherwise false *)    
val is_inside : seed_t -> block_t -> bool

(** [determine_separated_subseq block_ls order max_pos subseq_type] returns
* a list of separated subsequences which are created by
* using blocks as milestones. If order = First,the 
* separated subseqs is for first chromosome, otherwise the second chromosome *)
val determine_separated_subseq :
  block_t list ->
  order_t -> int -> [> `Alied | `Both | `Deleted ] -> subseq_t list

(** [create_alied_block_ls block_ls ali_pam seq1 seq2 cost_mat]
* creates the alignment for all blocks of [block_ls] based on their seed_ls *)
val create_alied_block_ls :
  block_t list ->
  pairChromPam_t -> Sequence.s -> Sequence.s -> Cost_matrix.Two_D.m -> unit
val check_sep_block : block_t list -> unit

(** [prepen b1 b2 block_pam seq1 seq2 cost_mat ali_pam] appends
block [b2] to block [b1] *)
val prepen :
  block_t ->
  block_t ->
  blockPam_t ->
  Sequence.s -> Sequence.s -> Cost_matrix.Two_D.m -> pairChromPam_t -> unit

(** [connect_pos_consecutive_block block_ls block_pam seq1 seq2 cost_mat ali_pam]
* connect consecutive positive blocks together to create large blocks. This functions
* returns a concatenated blocks list *)
val connect_pos_consecutive_block :
  block_t list ->
  blockPam_t ->
  Sequence.s ->
  Sequence.s -> Cost_matrix.Two_D.m -> pairChromPam_t -> block_t list

(** [connect_consecutive_block block_ls block_pam seq1 seq2 cost_mat ali_pam]
* connect consecutive  blocks together to create large blocks. This functions
* returns a concatenated blocks list *)
val connect_consecutive_block :
  block_t list ->
  blockPam_t ->
  Sequence.s ->
  Sequence.s -> Cost_matrix.Two_D.m -> pairChromPam_t -> block_t list

(** [create_subseq_id subseq_type sep_block_ls ali_pam] returns
* two lists  of separated subsequences for two chromosomes 
* which are created by using [sep_block_ls] as milestones *)
val create_subseq_id :
  [> `Alied | `Both | `Deleted ] ->
  block_t list ->
  pairChromPam_t -> block_t list * subseq_t list * subseq_t list


(** [create_median] approx block cost_mat] returns
* the median sequence and the cost of [block] *)
val create_median : ?approx:ChromPam.order_t -> block_t -> Cost_matrix.Two_D.m -> Sequence.s * int

(** [find_block block_ls subseq1_id subseq2_id] returns
* the blocks whose subseq ids are [subseq1_id] and [subseq2_id] *)
val find_block : block_t list -> int -> int -> block_t option

(** [find_subseq1 block_ls subseq1_id] return the subseq1
* of the block whose subseq1 is [subseq1_id] *)
val find_subseq1 : block_t list -> int -> block_t option
