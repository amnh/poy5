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

(** Genome module implements functions to create medians
    between two lists of genomes *)

val fprintf : out_channel -> ('a, out_channel, unit) format -> 'a
type med_t = GenomeAli.med_t

(** [meds_t] is a data structure for a list of  medians 
* between two genomes. Rearrangements are allowed *)
type meds_t = {
    med_ls : med_t list; (** genome list *)
    total_cost : int;   (** the total cost to create this genome list *)
    total_recost : int;  (** total recost to create this genome list *)
    c2 : Cost_matrix.Two_D.m; 

    approx_cost_arr : int array; (* An array of genome medians 
                                  * created from this median to other medians.
                                  * This is used for approximation *)
    approx_recost_arr : int array;
    code : int; (** the taxa code containing this genome list *)

    chrom_pam : Data.dyna_pam_t (** user-defined parameters used for creating genome medians *)   
}

(** [init_med genome chrom_pam tcode num_taxa] 
* returns a genome list with only one element
* created from genome [genome] *) 
val init_med :
  Sequence.s Data.dyna_data -> Data.dyna_pam_t -> int -> int -> meds_t

(** [keep chrom_pam med_ls] selects randomly a subset 
* of medians kept to process further *)
val keep : Data.dyna_pam_t -> 'a list -> 'a list

(** [update_cost_mat meds1 meds2] creates the median
* between genome [meds1] and [meds2]
* and update their approximate cost, recost
* arrays if they are not yet computed *)
val update_cost_mat : meds_t -> meds_t -> unit

(** [find_meds2 keep_all_meds meds1 meds2] finds median list
 * between two lists of medians [meds1=(x1,...,xk)] and [meds2=(y1,...,yt)]
 * where xi and yj are medians. For each pair (xi, yj) we have 
 * a list of medians z_ij with the same cost c_ij. 
 * Find z*_ij = minargv(z_ij )(c_ij) *)
val find_meds2 : ?keep_all_meds:bool -> meds_t -> meds_t -> meds_t

(** [find_meds3 medsp meds1 meds2] creates the median list
 * of three lists of medians [medsp=(x1,...,xk)], [meds1=(y1,...,yt)]
 * and [meds2=(z1,...,zq)] where xi, yj, and zp are medians. 
 * For each triplet (xi, yj, zp) we have 
 * a list of medians w_ijp with the same cost c_ijp. 
 * Find w*ijp = minargv_(w_ijp) (c_ijp) *)
val find_meds3 : meds_t -> meds_t -> meds_t -> meds_t

(** [cmp_min_pair_cost] computes the min median cost
 * between two lists of medians [meds1=(x1,...,xk)] and [meds2=(y1,...,yt)]
 * where xi and yj are medians. For each pair (xi, yj) we have 
 * a list of medians z_ij with the same cost c_ij. 
 * returns c*_ij = min (c_ij) *)
val cmp_min_pair_cost : meds_t -> meds_t -> int * int

(** [cmp_max_pair_cost] computes the max median cost
 * between two lists of medians [meds1=(x1,...,xk)] and [meds2=(y1,...,yt)]
 * where xi and yj are medians. For each pair (xi, yj) we have 
 * a list of medians z_ij with the same cost c_ij. 
 * returns c*_ij = max (c_ij) *)
val cmp_max_pair_cost : meds_t -> meds_t -> int * int

(** [Compare meds1 meds2] returns 0 if these 
* two lists [meds1] and [meds2] are the same, 
* otherwise (-1) or 1 *)
val compare : meds_t -> meds_t -> int

val clean_median : 'a -> 'b -> 'a * 'b

(** [get_active_ref_code meds] returns active 
* reference codes of genome medians [meds] *)
val get_active_ref_code : meds_t -> int * int * int

(** [copy_chrom_map s_ch d_ch] copies genome map
* from median [s_ch] into median [d_ch] *)
val copy_chrom_map : meds_t -> meds_t -> meds_t

(** [readjust_3d ch1 ch2 mine c2 c3 parent] readjusts
* the current median [mine] of three medians [ch1],
* [ch2], and [parent] using three dimentional alignments*)
val readjust_3d :
  meds_t ->
  meds_t ->
  meds_t ->
  Cost_matrix.Two_D.m ->
  Cost_matrix.Three_D.m -> meds_t -> int * meds_t * bool
