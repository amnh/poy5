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
let () = SadmanOutput.register "ChromPam" "$Revision: 2678 $"

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
type order_t = [`First | `Second | `BothSeq]

(** order_t is a data type with three values
* - `Positive: The orientation is positive, or consider only positive segments and characters
* - `Negative: The orientation is negative, or consider only negative segments and characters
* - `BothDir:  Consider only negative segments and characters *)
type direction_t = [`Positive | `Negative | `BothDir]

type re_meth_t = Data.re_meth_t

type annotate_tool_t = Data.annotate_tool_t

(** Parameters used to align two chromosomes *)
type chromPairAliPam_t = {
    max_gap                 : int; (** max distance between two basic seeds *)
    sig_k                   : int;  (** min significant seed length 
    what is this?*)
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

    annotate_tool : annotate_tool_t;
}

let locus_indel_cost_default = (10, 100)

let chromPairAliPam_default = {
    min_pos1 = -1; max_pos1 = -1; 
    min_pos2 = -1; max_pos2 = -1; 
    max_gap = 10; 
    sig_k = 12;
    (*k = 9;  
    sig_block_len = 100;
    rearranged_len = 100; *)
    
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
    kept_wag = 1;

    annotate_tool = `Default (100,100,9)
}


let get_chrom_pam user_chrom_pam =   
    let chrom_pam = chromPairAliPam_default  in

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

    let chrom_pam = 
        match user_chrom_pam.Data.max_kept_wag with
        | None -> chrom_pam
        | Some max_kept_wag -> {chrom_pam with kept_wag = max_kept_wag}
    in 

    let chrom_pam =
        match user_chrom_pam.Data.annotate_tool with
        | None -> chrom_pam
        | Some annotate_tool -> {chrom_pam with annotate_tool = annotate_tool}
    in

    chrom_pam
    


let cloneChromPairPam (donor : chromPairAliPam_t) = {
    max_gap = donor.max_gap;
    sig_k = donor.sig_k;
    (*sig_block_len = donor.sig_block_len;
    k = donor.k;
    rearranged_len = donor.rearranged_len; *)
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
    kept_wag = donor.kept_wag;

    annotate_tool = donor.annotate_tool;

}



let print_pair_ali_pam (pair_ali_pam : chromPairAliPam_t) = 
    print_endline "==================================";

    let print = Printf.fprintf in
    print stdout "(%i, %i) <-> (%i, %i) \n" pair_ali_pam.min_pos1
        pair_ali_pam.max_pos1 pair_ali_pam.min_pos2 pair_ali_pam.max_pos2;

    (*print stdout "k: %i\n" pair_ali_pam.k;*)

    print stdout "max_gap: %i\n" pair_ali_pam.max_gap;


    print stdout "seq1: (%i, %i) -> %i\n" pair_ali_pam.min_pos1

    pair_ali_pam.max_pos1 (pair_ali_pam.max_pos1 - pair_ali_pam.min_pos1 + 1);

	print stdout "seq2: (%i, %i) -> %i \n" pair_ali_pam.min_pos2
        pair_ali_pam.max_pos2 (pair_ali_pam.max_pos2 - pair_ali_pam.min_pos2 + 1);

	print_endline "=================================="


let get_min_seed_len pam =
    match pam.annotate_tool with
    | `Mauve (_,_,_)  -> failwith "Mauve annotator does not have this" 
    | `Default (min_loci_len,min_rearrangement_len,min_seed_length) ->
            min_seed_length

let get_min_loci_len pam =
    match pam.annotate_tool with
    | `Mauve (_,_,_)  -> failwith "Mauve annotator does not have this" 
    | `Default (min_loci_len,min_rearrangement_len,min_seed_length) ->
            min_loci_len

let get_min_rearrangement_len pam =
    match pam.annotate_tool with
    | `Mauve (_,_,_)  -> failwith "Mauve annotator does not have this" 
    | `Default (min_loci_len,min_rearrangement_len,min_seed_length) ->
            min_rearrangement_len

let get_min_cover_ratio pam = 
    match pam.annotate_tool with
    | `Mauve (min_lcb_ratio,min_cover_ratio,bk_penalty) -> min_cover_ratio
    | `Default (_,_,_) ->
            failwith "Default annotator does not have this"

let get_min_lcb_ratio pam =
    match pam.annotate_tool with
    | `Mauve (min_lcb_ratio,min_seed_num,bk_penalty) ->  min_lcb_ratio
    | `Default (_,_,_) ->
            failwith "Default annotator does not have this"

let get_min_bk_penalty pam =
    match pam.annotate_tool with
    | `Mauve (min_lcb_ratio,min_seed_num,min_bk_penalty) -> min_bk_penalty
    | `Default (_,_,_) ->
            failwith "Default annotator does not have this"



