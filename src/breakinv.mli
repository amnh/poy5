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


(** Median module contains functions to create medians
*    between two lists of breakinv chracters *)


val fprintf : out_channel -> ('a, out_channel, unit) format -> 'a
type breakinv_t = BreakinvAli.breakinv_t

(** [meds_t] is a data structure for a list medians 
* between two breakinv character lists *)
type meds_t = {
    med_ls : breakinv_t list; (** breakinv list *)
    num_med : int; (** number of breakinv characters *)
    total_cost : int;   (** the cost to create this breakinv list *)
    total_recost : int; (** the recost to create this breakinv list *)
    breakinv_pam : Data.dyna_pam_t; (** breakinv paramenters used to create breakinv median *)
    gen_cost_mat : Cost_matrix.Two_D.m;
    pure_gen_cost_mat : int array array;
    alpha : Alphabet.a 
}


(** init_med seq gen_cost_mat alpha breakinv_pam] returns
* a breakinv character list with only one element 
* created from a sequence of general character [seq]*)
val init_med :
  Sequence.s ->
  Cost_matrix.Two_D.m -> Alphabet.a -> Data.dyna_pam_t -> meds_t

(** [keep chrom_pam med_ls] returns a sublist of median list
* [med_ls] to be kept to process further based on the customs's defined paramaters *)
val keep : Data.dyna_pam_t -> 'a list -> 'a list

(** [find_meds2 meds1 meds2] returns a list of 
* breakinv character medians created for two lists of medians 
* [meds1=(x1,...,xk)] and [meds2=(y1,...,yt)]
* where xi and yj are medians. For each pair (xi, yj) we have 
* a list of medians z_ij with the same cost c_ij. 
* Find z*_ij = minargv(z_ij )(c_ij) *)
val find_meds2 : meds_t -> meds_t -> meds_t


(** [cmp_min_pair_cost] returns the minimum cost
* between two lists of medians [meds1=(x1,...,xk)] and [meds2=(y1,...,yt)]
 * where xi and yj are medians. For each pair (xi, yj) we have 
 * a list of medians z_ij with the same cost c_ij. 
 * Find c*_ij = min (c_ij) *)
val cmp_min_pair_cost : meds_t -> meds_t -> int * int

(** [cmp_max_pair_cost] returns the maximum cost
* between two lists of medians [meds1=(x1,...,xk)] and [meds2=(y1,...,yt)]
 * where xi and yj are medians. For each pair (xi, yj) we have 
 * a list of medians z_ij with the same cost c_ij. 
 * Find c*_ij = min (c_ij) *)
val cmp_max_pair_cost : meds_t -> meds_t -> int * int


(** [find_meds3 medsp meds1 meds2] returns 
* the median of three lists of medians [medsp=(x1,...,xk)], [meds1=(y1,...,yt)]
 * and [meds2=(z1,...,zq)] where xi, yj, and zp are medians. 
 * For each triplet (xi, yj, zp) we have 
 * a list of medians w_ijp with the same cost c_ijp. 
 * Find w*ijp = minargv_(w_ijp) (c_ijp) *)
val find_meds3 : meds_t -> meds_t -> meds_t -> meds_t


(** [readjust_3d ch1 ch2 mine c2 c3 parent] readjusts
* the breakinv median [mine] of three breakinv medians 
* [ch1], [ch2] and [parent] *) 
val readjust_3d : meds_t -> meds_t -> meds_t -> 'a -> 'b -> meds_t -> int * meds_t * bool

(** [compare meds1 meds2] returns 0 if breakinv list [meds1]
* is the same as breakinv list [meds2], otherwise (-1) or (1) *)
val compare : meds_t -> meds_t -> int

(** [get_active_ref_code meds] return active reference codes
* of breakinv medians [meds] *)
val get_active_ref_code : meds_t -> int * int * int


