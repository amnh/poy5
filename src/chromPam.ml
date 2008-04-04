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
let () = SadmanOutput.register "ChromPam" "$Revision: 2665 $"

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

type order_t = [`First | `Second | `BothSeq]
type direction_t = [`Positive | `Negative | `BothDir]

type re_meth_t = Data.re_meth_t



(** Parameters used to align two chromosomes *)
type chromPairAliPam_t = {
    (** The minimum length of a basic seed *)
    k                       : int; 
    sig_k                   : int;
    sig_block_len           : int;
    (** The maximum distance between two basic seeds *)
    max_gap                 : int;

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

    
    (** The cost of an inversion operation. At the moment, the total cost of
        rearrangements is computed in the term of inversion distance  *)
    re_meth : re_meth_t;
    chrom_breakpoint : int;

    keep_median : int;
    swap_med : int;

    approx : order_t;
    symmetric : bool;
    negative : bool;

    (** The cost to delete or insert a segment (loci) in a chromosome *)
    locus_indel_cost : (int * int);
    chrom_indel_cost : (int * int);
    chrom_hom : int; (* if cost > threshold * min_cost,  then not homologous *)
    circular : int;


    max_3d_len : int;
    detected_3d_len : int;
}


let locus_indel_cost_default = (10, 100)

(** TODO: Get rid of this set of defaults!!!
* WARNING: The true defaults are in Data!!! this defaults need to be cleaned up.
* *)
let chromPairAliPam_default = {
    min_pos1 = -1; max_pos1 = -1; min_pos2 = -1; max_pos2 = -1; 


    k = 9;  (* Seed length *)
    max_gap = 10; 
    sig_k = 12;
    sig_block_len = 100;
    rearranged_len = 100;
    
    gap_opening_cost = 200;
    gap_ext_cost = 10;
    mat_cost = -100;    
    mismat_cost = 50;

    re_meth  = `Locus_Breakpoint 10;
    chrom_breakpoint = 100;

    keep_median = 1;
    swap_med = 1;

    approx = `BothSeq;
    symmetric = false;
    negative = true;
    
    locus_indel_cost = (10, 100);
    chrom_indel_cost = (10, 100);
    chrom_hom = 200;
    
    circular = 0;

    max_3d_len = max_int;
    detected_3d_len = 200;
}


let get_chrom_pam user_chrom_pam =     
    let chrom_pam = chromPairAliPam_default in 
    let chrom_pam = 
        match user_chrom_pam.Data.seed_len with 
        | None -> chrom_pam
        | Some k -> {chrom_pam with k = k}
    in 

    let chrom_pam = 
        match user_chrom_pam.Data.re_meth with
        | None -> chrom_pam
        | Some re_meth -> {chrom_pam with re_meth = re_meth}
    in 

    let chrom_pam = 
        match user_chrom_pam.Data.chrom_breakpoint with
        | None -> chrom_pam
        | Some chrom_breakpoint -> {chrom_pam with chrom_breakpoint = chrom_breakpoint}
    in 

    let chrom_pam = 
        match user_chrom_pam.Data.circular with
        | None -> chrom_pam
        | Some circular -> {chrom_pam with circular = circular}
    in 

    let chrom_pam = 
        match user_chrom_pam.Data.keep_median with
        | None -> chrom_pam
        | Some keep_median -> {chrom_pam with keep_median = keep_median}
    in 

    let chrom_pam = 
        match user_chrom_pam.Data.swap_med with
        | None -> chrom_pam
        | Some swap_med -> {chrom_pam with swap_med = swap_med}
    in 

    let chrom_pam = 
        match user_chrom_pam.Data.sig_block_len with
        | None -> chrom_pam
        | Some sig_block_len -> {chrom_pam with sig_block_len = sig_block_len}
    in 


    let chrom_pam = 
        match user_chrom_pam.Data.rearranged_len with
        | None -> chrom_pam
        | Some rearranged_len -> {chrom_pam with rearranged_len = rearranged_len}
    in 


    let chrom_pam = 
        match user_chrom_pam.Data.approx with
        | None -> chrom_pam
        | Some approx -> 
              if approx then {chrom_pam with approx = `First}
              else {chrom_pam with approx = `BothSeq}
    in 


    let chrom_pam = 
        match user_chrom_pam.Data.symmetric with
        | None -> chrom_pam
        | Some sym -> {chrom_pam with symmetric = sym}
    in 

    let chrom_pam =
        match user_chrom_pam.Data.locus_indel_cost with  
        | None -> chrom_pam
        | Some locus_indel_cost -> {chrom_pam with locus_indel_cost = locus_indel_cost}
    in 

    let chrom_pam =
        match user_chrom_pam.Data.chrom_indel_cost with  
        | None -> chrom_pam
        | Some chrom_indel_cost -> {chrom_pam with chrom_indel_cost = chrom_indel_cost}
    in 


    let chrom_pam =
        match user_chrom_pam.Data.max_3d_len with  
        | None -> chrom_pam
        | Some l -> {chrom_pam with max_3d_len = l}
    in 

    chrom_pam
    


let cloneChromPairPam (donor : chromPairAliPam_t) = {
    k = donor.k;
    max_gap = donor.max_gap;
    sig_block_len = donor.sig_block_len;
    sig_k = donor.sig_k;
    rearranged_len = donor.rearranged_len;

    
    min_pos1 = donor.min_pos1;
    min_pos2 = donor.min_pos2;
    max_pos1 = donor.max_pos1;
    max_pos2 = donor.max_pos2;
    
    gap_opening_cost = donor.gap_opening_cost;
    gap_ext_cost = donor.gap_ext_cost;
    mat_cost = donor.mat_cost;
    mismat_cost = donor.mismat_cost;

    re_meth  = donor.re_meth;
    keep_median  = donor.keep_median;
    swap_med  = donor.swap_med;
    approx  = donor.approx;
    negative = donor.negative;

    locus_indel_cost = donor.locus_indel_cost;
    chrom_indel_cost = donor.chrom_indel_cost;
    chrom_breakpoint = donor.chrom_breakpoint;
    chrom_hom = donor.chrom_hom;

    circular = donor.circular;
    symmetric = donor.symmetric;

    max_3d_len = donor.max_3d_len;
    detected_3d_len = donor.detected_3d_len;
}



let print_pair_ali_pam (pair_ali_pam : chromPairAliPam_t) = 
    print_endline "==================================";

    let print = Printf.fprintf in
    print stdout "(%i, %i) <-> (%i, %i) \n" pair_ali_pam.min_pos1
        pair_ali_pam.max_pos1 pair_ali_pam.min_pos2 pair_ali_pam.max_pos2;

    print stdout "k: %i\n" pair_ali_pam.k;

    print stdout "max_gap: %i\n" pair_ali_pam.max_gap;


    print stdout "seq1: (%i, %i) -> %i\n" pair_ali_pam.min_pos1

    pair_ali_pam.max_pos1 (pair_ali_pam.max_pos1 - pair_ali_pam.min_pos1 + 1);

	print stdout "seq2: (%i, %i) -> %i \n" pair_ali_pam.min_pos2
        pair_ali_pam.max_pos2 (pair_ali_pam.max_pos2 - pair_ali_pam.min_pos2 + 1);

	print_endline "=================================="
