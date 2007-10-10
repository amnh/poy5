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
type blockPam_t = { max_connecting_dis : int; max_connecting_shift : int; }
type block_t = {
  mutable id : int;
  mutable is_dum : bool;
  mutable sta1 : int;
  mutable sta2 : int;
  mutable en1 : int;
  mutable en2 : int;
  mutable direction : direction_t;
  mutable cost : int;
  mutable alied_seq1 : Sequence.s option;
  mutable alied_seq2 : Sequence.s option;
  mutable seed_ls : seed_t list;
  mutable subseq1_id : int;
  mutable subseq2_id : int;
}
val blockPam_default : blockPam_t
val cloneBlockPam : blockPam_t -> blockPam_t
val create_from_seed : int -> seed_t -> block_t
val create_simple_block : int -> int -> int -> int -> int -> block_t
val get_dum_first_block : ChromPam.chromPairAliPam_t -> block_t
val get_dum_last_block : ChromPam.chromPairAliPam_t -> block_t
val max_len : block_t -> int
val invert : block_t -> int -> int -> unit
val get_pos : block_t -> order_t -> int * int
val cmp_dia_dis : block_t -> block_t -> int
val print : block_t -> unit
val add_seed : block_t -> seed_t -> unit
val cmp_cost_based_seed : block_t -> pairChromPam_t -> int
val create_from_seed_ls : int -> seed_t list -> pairChromPam_t -> block_t
val cmp_ali_cost :
  Sequence.s -> Sequence.s -> direction_t -> pairChromPam_t -> int
val find_local_block : seed_t list -> pairChromPam_t -> block_t list
val is_free : block_t -> int array -> int array -> bool
val assign : block_t -> int array -> int array -> unit
val create_sig_arr :
  block_t list -> ChromPam.chromPairAliPam_t -> int array * int array
val select_separated_block :
  block_t list -> ChromPam.chromPairAliPam_t -> block_t list
val create_pos_alied_block :
  block_t ->
  Sequence.s -> Sequence.s -> Cost_matrix.Two_D.m -> pairChromPam_t -> unit
val is_inside : seed_t -> block_t -> bool
val determine_separated_subseq :
  block_t list ->
  order_t -> int -> [> `Alied | `Both | `Deleted ] -> subseq_t list
val create_alied_block_ls :
  block_t list ->
  pairChromPam_t -> Sequence.s -> Sequence.s -> Cost_matrix.Two_D.m -> unit
val check_sep_block : block_t list -> unit
val prepen :
  block_t ->
  block_t ->
  blockPam_t ->
  Sequence.s -> Sequence.s -> Cost_matrix.Two_D.m -> pairChromPam_t -> unit
val connect_pos_consecutive_block :
  block_t list ->
  blockPam_t ->
  Sequence.s ->
  Sequence.s -> Cost_matrix.Two_D.m -> pairChromPam_t -> block_t list
val connect_consecutive_block :
  block_t list ->
  blockPam_t ->
  Sequence.s ->
  Sequence.s -> Cost_matrix.Two_D.m -> pairChromPam_t -> block_t list
val create_subseq_id :
  [> `Alied | `Both | `Deleted ] ->
  block_t list ->
  pairChromPam_t -> block_t list * subseq_t list * subseq_t list
val create_median : ?approx:ChromPam.order_t -> block_t -> Cost_matrix.Two_D.m -> Sequence.s * int
val find_block : block_t list -> int -> int -> block_t option
val find_subseq1 : block_t list -> int -> block_t option
