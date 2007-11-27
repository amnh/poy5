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


val fprintf : out_channel -> ('a, out_channel, unit) format -> 'a
type seq_t = {
  seq : Sequence.s;
  seq_ref_code : int;
  alied_med : Sequence.s;
  seq_ord1 : int;
  alied_seq1 : Sequence.s;
  seq_ord2 : int;
  alied_seq2 : Sequence.s;
}
type annchrom_t = {
  seq_arr : seq_t array;
  ref_code : int;
  ref_code1 : int;
  ref_code2 : int;
  cost1 : int;
  recost1 : int;
  cost2 : int;
  recost2 : int;
}
val clone_seq : seq_t -> seq_t
val clone_med : annchrom_t -> annchrom_t
type annchromPam_t = {
  re_meth : Data.re_meth_t;
  keep_median : int;
  circular : int;
  swap_med : int;
  approx : ChromPam.order_t;
  symmetric : bool;
  locus_indel_cost : int * int;
}
val annchromPam_default : annchromPam_t
val init_seq_t : Sequence.s * int -> seq_t
val init : (Sequence.s * int) array -> annchrom_t
val printMap : seq_t array -> unit
val get_seq_arr : annchrom_t -> Sequence.s array
val convert_map : annchrom_t -> (int * int * int * int * int * int) list
val get_annchrom_pam : Data.dyna_pam_t -> annchromPam_t
val print : annchrom_t -> Alphabet.a -> unit
val split : annchrom_t -> Sequence.s array * int array
val create_pure_gen_cost_mat :
  Sequence.s array ->
  Sequence.s array ->
  Cost_matrix.Two_D.m ->
  annchromPam_t -> int array array * int array * int array * int
val cmp_cost :
  annchrom_t ->
  annchrom_t -> Cost_matrix.Two_D.m -> 'a -> Data.dyna_pam_t -> int * int
val find_med2_ls :
  annchrom_t ->
  annchrom_t ->
  Cost_matrix.Two_D.m -> 'a -> Data.dyna_pam_t -> int * int * annchrom_t list

val find_med3 :
  annchrom_t ->
  annchrom_t ->
  annchrom_t ->
  annchrom_t ->
  Cost_matrix.Two_D.m ->
  Cost_matrix.Three_D.m -> 'a -> Data.dyna_pam_t -> int * annchrom_t

val compare : annchrom_t -> annchrom_t -> int
val find_approx_med2 : annchrom_t -> annchrom_t -> annchrom_t -> annchrom_t
val assign_seq_ref : annchrom_t -> int -> annchrom_t * int
val create_map : annchrom_t -> int -> int * int * Tags.output
val create_single_map : annchrom_t -> Tags.output


val to_single : annchrom_t -> int -> Cost_matrix.Two_D.m -> Sequence.s array
val to_single_root :
  annchrom_t -> int -> Cost_matrix.Two_D.m -> Sequence.s array

val change_to_single : annchrom_t -> Sequence.s array -> annchrom_t
val to_formater : annchrom_t -> Alphabet.a -> string
val copy_chrom_map : annchrom_t -> annchrom_t -> annchrom_t
