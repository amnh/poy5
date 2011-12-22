(* POY 4.0 Beta. A phylogenetic analysis program using Dynamic Homologies.    *\
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
\* USA                                                                        *)


type m_i = {  (* the ith element M_i of a local mum *)
    sequence_NO : int ; (*which sequence is this m_i in *)
    left_end : int;  (*left end coordination of M_i*)
    right_end : int;
    orientation : int; (*1 or -1, 1 is positive, -1 is reverse complement*)
    }

type mum = { 
    seedNO : int; (*the # of seed that contruct this mum. start from
    only one seed, than extend*)
    mumseq : int array;
    positions : m_i list;  (* the positions in sequence this mum shows up*)
    mumkey : int;
    size : int; (* size of this mum, also the size of position list *)
    (* we don't extend seed with subset/superset for multi-sequence.
    left_superset : int list  ; (* seedNOlst of another mum*)
    right_superset : int list;
    left_subset : subset list  ; (* seedNOlst of another mum*)
    right_subset : subset list;
    *)
    neighborhood_lst: ( int * int * int * int * int * int ) list; 
    (* list of neighborhood (seqNO,j_seedNO,i_ori,j_ori,distance left,distance right) 
     * we can get seqNO from positions.(i_idx).seqNO*)
    subsuming_pointer : int ; 
    extendable : int ; 
    (*extendable= 0,1,2,3:
    0. this seed is extendable (from). 
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
    * this mum to ext=0 if this priority_lst is empty*)
    mumscore : int ; 
    (*mumscore is not only used during build mum table, 
    * after we have the lcbs, if this mum is part of a huge lcb block, mumscore
    * is set to alignment score in function search_inside_lcb.*)
    mumalgn : Sequence.s list;
    (*alignment record the alignment result in function search_inside_lcb*)
}

type lcb = {
    seedNOlst : int list; (*each lcb is a list of mums,also the key to lcb_tbl*)
    range_lst : m_i list; (*lcb_range we need is exactly the same truct as m_i,
                           no need to create a new type here*)
    ratio : float; (* [score/length] of seq in this lcb*)
    ref_code : int; (*just a code for bk/rearr purpose*)
    avg_range_len : int; (*average length from range_lst*)
    (*score between subsequence contained by range_lst during lcb building.
    * after we have lcb_tbl, huge lcb blocks are aligned seperately, score in
    * function search_inside_each_lcb is set to alignment cost of seq in this lcb*)
    score : int; 
    alignment : Sequence.s array;
}


(** [get_matcharr_and_costmatrix] is the main function of this module. it take two
* sequences, set of parameters (min lcb ratio, min lcb length, etc) and 
* a costmatrix between charactors as input,
* output two code array of matching blocks, a cost_mat between lcb blocks and non-lcb
* blocks, ali_mat is the cost and alignment sequence between lcb blocks.
* NOTE: edit_cost is the total editing cost between lcb blocks (which is not included
* is cost_mat), full_cost_listlst is a little bit redundant here, for it
* contains the two code array, but it also has the begin and end point of each
* block. *)
val get_matcharr_and_costmatrix : Sequence.s -> Sequence.s -> float -> float ->
    float -> float ->
int*int -> Cost_matrix.Two_D.m ->  bool -> 
    int array * int array * int array array *
(int*Sequence.s*Sequence.s) array array * int * int * 
(int * (int * int)) list list * int

val output2mauvefile : string -> int -> (int option) -> int array -> int array
-> (int * (int * int)) list list -> (int*Sequence.s*Sequence.s) array
array ->int -> int -> int -> int -> unit 

    
(* we don't call this directly outside of mauve any more
val create_lcb_tbl : int array array -> float -> float -> int ->
    (int list, lcb) Hashtbl.t * int list list list * int list list * 
    (int * int) list list
*)
val print_mum : bool -> bool -> mum  -> unit

val print_lcb : lcb -> unit

val print_int_list : int list -> unit

val print_int_lstlst : int list list -> unit

val print_lcblst : int list list list -> unit

val get_abs_lst : int list -> int list

val update_lcb_ref_code : (int list, lcb) Hashtbl.t -> int list -> int -> unit

val get_position_by_seqNO : m_i list -> int -> m_i

