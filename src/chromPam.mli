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

(** Chromosome parameters
 *
 * This module contains all parameters used to create alignments
 * as well as medians of chromosomes, basically two or there sequences.
 * Two kinds of operations are allowed between two chromosomes: 
 * point mutations and rearrangements. Point mutations comprise
 * deletion, insertion and substitution of nucleotides. 
 * Rearrangements comprise translocation, inversion and
 * translocation-inversion of sub-sequences in the chromosome.
 * 
 *) 

(** order_t is a data type with three values
* - `First: Consider only the first chromsome
* - `Second: Consider only the second chromosome
* - `BothSeq: Consider both chromosomes *)
type order_t = [ `BothSeq | `First | `Second ]


(** order_t is a data type with three values
* - `Positive: The orientation is positive, or consider only positive segments and characters
* - `Negative: The orientation is negative, or consider only negative segments and characters
* - `BothDir:  Consider only negative segments and characters *)
type direction_t = [ `BothDir | `Negative | `Positive ]

type re_meth_t = Data.re_meth_t

(** Parameters used to align two chromosomes *)
type chromPairAliPam_t = {
    k                       : int;  (** min basic seed length *)
    sig_k                   : int;  (** min significant seed length *)
    sig_block_len           : int;  (** min signification block length *)
    max_gap                 : int; (** max distance between two basic seeds *)

    (** It's believed that no rearrangments or reversions happened 
        within a segment whose length < unbreaked_len *)
    rearranged_len           : int;

    (** Align the subsequence (min_pos1->max_pos1) of chromosome 1 to the
        subsequence (min_pos2->max_pos2) of chromosome 2. This information is
        reserved for devide-conquer techniques *)
    min_pos1                : int;
    max_pos1                : int;
    min_pos2                : int;
    max_pos2                : int;      

    (** self-meaningful parameters. 
        These costs are used to compute the cost of seeds and blocks.
        Note that these cost pams are different from cost_matrix used to align
        two sequences in Sequence.Align.align_2*)
    gap_opening_cost        : int;
    gap_ext_cost            : int;
    mat_cost                : int;
    mismat_cost             : int;

    
    (** Cost parameters of rearrangement function which is either
    * breakpoint distance or inversion distance *)
    re_meth : re_meth_t;

    (** The cost of a breakpoint happing between two chromosome *)
    chrom_breakpoint : int;
    
    (** The maximum number of medians at one node kept during the search*)
    keep_median : int;
    (** The maximum number of swap iterations emplying
        * to improve the pairwise alignment with rearrangements*)
    swap_med : int;

    approx : order_t;
    symmetric : bool;
    negative : bool;

    (** Indel cost of a locus in a chromosome *)
    locus_indel_cost : (int * int);

    (** Indel cost of a chromosome in a genome *)
    chrom_indel_cost : (int * int);

    chrom_hom : int; (* The maximum cost between two chromosomes
                      * at which they are considered as homologous chromosomes *)
    (** Circular = 0 means linear chromosome, otherwise circular chromosome *)
    circular : int;

    (** maximum length of sequences aligned by 3D-alignment *)
    max_3d_len : int;
    detected_3d_len : int;
    kept_wag : int;
}

val chromPairAliPam_default : chromPairAliPam_t
val get_chrom_pam : Data.dyna_pam_t -> chromPairAliPam_t
val cloneChromPairPam : chromPairAliPam_t -> chromPairAliPam_t
val print_pair_ali_pam : chromPairAliPam_t -> unit
val locus_indel_cost_default : int * int
