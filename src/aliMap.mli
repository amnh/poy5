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

type chromPairAliPam_t = ChromPam.chromPairAliPam_t
type block_pam_t = Block.blockPam_t
type seed_t = Seed.seed_t
type block_t = Block.block_t
type order_t = ChromPam.order_t
type subseq_t = Subseq.subseq_t
val fprintf : out_channel -> ('a, out_channel, unit) format -> 'a
val create_general_ali :
  [> `Breakinv ] ->
  Block.block_t list ->
  Sequence.s ->
  Sequence.s ->
  Cost_matrix.Two_D.m ->
  Block.pairChromPam_t ->
  Block.subseq_t list * Block.subseq_t list * Block.block_t list *
  (Sequence.s * Sequence.s) array array * int array * int array * int *
  (int * int)
