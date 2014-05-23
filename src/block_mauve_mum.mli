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


val get_proper_seedlen : int -> int

type m_i = {  (* the ith element M_i of a local mum *)
    sequence_NO : int ; (*which sequence is this m_i in *)
    left_end : int;  (*left end coordination of M_i*)
    right_end : int;
    orientation : int; (*1 or -1, 1 is positive, -1 is reverse complement*)
    }

type mum = {
    seedNO : int; (*the # of seed that contruct this mum. start from
    only one seed, than extend*)
    mumseq : int array; (*subsequence of seq0. NOTE: seq0*) 
    positions : m_i list;  (* the positions in sequence this mum shows up*)
    mumkey : int;
    size : int; (* size of this mum, also the size of position list *)
    neighborhood_lst: ( int * int * int * int * int * int ) list; 
    (* list of neighborhood (seqNO,j_seedNO,i_ori,j_ori,distance left,distance right) 
     * we can get seqNO from positions.(i_idx).seqNO*)
    subsuming_pointer : int ; 
    extendable : int ; 
    (*extendable= 0,1,2,3: when we say extendable, we mean "extendable from"
    0. this seed is extendable (from) 
    1 ~ 3 : not extendable (from),but might be extend to.
    1. during the scan_seqlst, it means this mum shows up more then
    once in a sequence, also shows up in every sequence,
    we should not extend the seed from this mum(but
    we can extend from other mum to this one),
    2 .during the resolve_overlap, it means this mum does not
    show up in one or more sequence. we don't extend to this one. 
    3. during add to/remove from position2seed table, only one seed is kept as
    extendable, others are marked as unextendable. but if we remove the
    extendable one for some reason, we can upgrade one of this kind to be
    extendable.
    *)
    priority_lst : int list; 
    (*priority_lst keep a list of seedNO that make this mum un-extendable (which
    * means this mum has extendable=3 because of those seeds) we can upgrade
    * this mum to ext=0 if this priority_lst is empty.
    * Note: each seedNO can show up more than once in priority_lst, since we can
    * have two MUM share more than one start/end positions*)
    mumscore : int ;
    mumalgn : Sequence.s list;
}

type mum_btree = ((int array),mum) BinaryTree.b_tree


val print_pos_list : m_i list -> unit

val print_pos_list_to_file : out_channel -> m_i list -> unit

val get_mum_from_mumtbl : int -> (int, (int array, mum) BinaryTree.b_tree) Hashtbl.t -> 
    (int, (int * int * int) list * int array * int) Hashtbl.t -> mum

val get_position_by_seqNO : m_i list -> int -> m_i

val print_mum : bool -> bool -> mum -> unit

val get_seed_weight : int -> int

val build_seed_and_position_tbl :
           int array array ->
           int ->
           (int , int list ) Hashtbl.t ->
           (int * int , (int * int * int) list) Hashtbl.t ->
           (int * int, (int * int * int) list) Hashtbl.t ->
           (int, (int * int * int) list * int array * int) Hashtbl.t ->
           (int, (int array, mum) BinaryTree.b_tree) Hashtbl.t ->
           int array ref -> bool -> int

val build_local_mums2  :
    (int, (int array, mum) BinaryTree.b_tree) Hashtbl.t ->
    (int, (int * int * int) list * int array * int) Hashtbl.t ->
    (int * int, (int * int * int) list) Hashtbl.t ->
    (int * int, (int * int * int) list) Hashtbl.t ->
    bool -> bool -> unit

val transpose_mum_back :
    (int*int*int list*int) list list ->
    (int, (int array, mum) BinaryTree.b_tree) Hashtbl.t ->
    (int , int list ) Hashtbl.t ->
    (int, (int array, mum) BinaryTree.b_tree) Hashtbl.t ->
    (int * int, (int * int * int) list) Hashtbl.t ->
    (int * int, (int * int * int) list) Hashtbl.t->
    (int, (int * int * int) list * int array * int) Hashtbl.t ->
    (int array) ref -> unit

val extend_seeds :
(int, (int array, mum) BinaryTree.b_tree) Hashtbl.t ->
(int, (int * int * int) list * int array * int) Hashtbl.t ->
(int * int, (int * int * int) list) Hashtbl.t ->
(int * int, (int * int * int) list) Hashtbl.t -> unit

val resolve_overlap_mum :
(int, (int array, mum) BinaryTree.b_tree) Hashtbl.t ->
int ->
(int * int, (int * int * int) list) Hashtbl.t ->
(int * int, (int * int * int) list) Hashtbl.t ->
(int, (int * int * int) list * int array * int) Hashtbl.t ->
(int array) ref -> bool -> bool -> int

val update_score_for_each_mum :
(int, (int array, mum) BinaryTree.b_tree) Hashtbl.t ->
int array array -> unit


val get_mum_lst_for_each_seq :
(int, (int array, mum) BinaryTree.b_tree) Hashtbl.t ->
(int, (int * int * int) list * int array * int) Hashtbl.t ->
(int * int, (int * int * int) list) Hashtbl.t ->
int  -> int list -> int list list list


val print_mumtbl : 
(int, (int array, mum) BinaryTree.b_tree) Hashtbl.t ->
bool -> unit



(*
* seed2pos_tbl
* (int, (int * int * int) list * int array * int) Hashtbl.t
* mumtbl 
* (int, (int array, mum) BinaryTree.b_tree) Hashtbl.t
* pos2seed_tbl
* (int * int, (int * int * int) list) Hashtbl.t
* seedNO2seq_tbl
* (int , int list ) Hashtbl.t
*)
