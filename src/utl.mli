(* POY 5.1.1. A phylogenetic analysis program using Dynamic Homologies.       *)
(* Copyright (C) 2014 Andrés Varón, Lin Hong, Nicholas Lucaroni, Ward Wheeler,*)
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

val large_int : int

val break_code : int -> int array

val max_seq_len : int

val gen_chrom_ref_code : int ref

val gen_seq_ref_code : int ref

val gen_genome_ref_code : int ref

val get_new_chrom_ref_code : unit -> int

val get_new_genome_ref_code : unit -> int

val get_new_seq_ref_code : unit -> int

val deref : 'a option -> 'a

val compare_non_dec_list : int list -> int list -> bool

val get_sum_arr : int array -> int -> int -> int

val invert_subarr : 'a array -> int -> int -> unit

val binary_index_search : int array -> int -> int

val find_index : 'a array -> 'b -> ('b -> 'a -> int) -> int

val get_common :
  'a array -> 'a array -> ('a -> 'a -> int) -> 'a array * 'a array

val get_common3 :
  'a array -> 'a array -> 'a array -> ('a -> 'a -> int) -> 'a array * 'a array * 'a array

val insert : 'a array -> int -> 'a -> 'a array

val swap_item : int -> int -> 'a array -> 'a array

val bigger_int : int -> int -> int

val float_to_int_mat : float array array -> int array array

val int_to_int32_arr : int array -> int32 array

val int_to_int32_mat : int array array -> int32 array array

val printIntArr : int array -> unit

val printIntArrWithIdx : int array -> unit

val printIntMat : int array array -> unit

val printIntMatWithIdx : int array array -> unit

(** printIntList print a list of int, with a newline.*)
val printIntList : int list -> unit

(** printIntList print a list of int, without a newline.*)
val printIntList2 : int list -> unit

val printIntListList : int list list -> unit

val printIntListList2 : int list list -> unit

val printIntListListList : int list list list -> unit

val printIntListToFile : out_channel -> int list -> unit

val create_ls : int -> 'a -> 'a list

val get_neg_rev_intlst : int list -> int list

val get_abs_intlst : int list -> int list

val remove_nth : ?acc:'a list -> 'a list -> int -> 'a * 'a list

val insert_arr : int array -> int array -> int -> int array

val pairwisep : ('a -> 'a -> bool) -> 'a list -> bool

val get_k_random_elem : 'a list -> int -> 'a list

val isEqualArr : 'a array -> 'b array -> ('a -> 'b -> int) -> bool

val break_array : 'a array -> (int * int) list -> 'a array list

val max_arr : 'a array -> 'a

val min_arr : 'a array -> 'a

val get_dir : [> `Negative | `Positive ] -> string

val factorial :int -> int

val p_m_n : int -> int -> int

val get_avg_of_intlst : int list -> float

val get_min_of_lst : 'a list -> 'a

val get_avg_of_floatlst : float list -> float

