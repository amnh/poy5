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

type order_t = [ `BothSeq | `First | `Second ]
type direction_t = [ `BothDir | `Negative | `Positive ]
type re_meth_t = Data.re_meth_t
type chromPairAliPam_t = {
  k : int;
  sig_k : int;
  sig_block_len : int;
  max_gap : int;
  rearranged_len : int;
  min_pos1 : int;
  max_pos1 : int;
  min_pos2 : int;
  max_pos2 : int;
  gap_opening_cost : int;
  gap_ext_cost : int;
  mat_cost : int;
  mismat_cost : int;
  re_meth : re_meth_t;
  chrom_breakpoint : int;
  keep_median : int;
  swap_med : int;
  approx : order_t;
  symmetric : bool;
  negative : bool;
  locus_indel_cost : int * int;
  chrom_indel_cost : int * int;
  chrom_hom : int;
  circular : int;
  max_3d_len : int;
}
val chromPairAliPam_default : chromPairAliPam_t
val get_chrom_pam : Data.dyna_pam_t -> chromPairAliPam_t
val cloneChromPairPam : chromPairAliPam_t -> chromPairAliPam_t
val print_pair_ali_pam : chromPairAliPam_t -> unit
val locus_indel_cost_default : int * int
