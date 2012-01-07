(* POY 5.0 Alpha. A phylogenetic analysis program using Dynamic Homologies.   *\
(* Copyright (C) 2011  Andrés Varón, Lin Hong, Nicholas Lucaroni, Ward Wheeler*)
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


type lcb = {
    seedNOlst : int list; (*each lcb is a list of mums,also the key to lcb_tbl*)
    range_lst : Block_mauve_mum.m_i list; (*lcb_range we need is exactly the same truct as m_i,
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
(int * (int * int) * int list) list list * int

val output2mauvefile : string -> int -> (int option) -> int array -> int array
-> (int * (int * int) * int list) list list -> (int*Sequence.s*Sequence.s) array
array ->int -> int -> int -> int -> unit 

    
(* we don't call this directly outside of mauve any more
val create_lcb_tbl : int array array -> float -> float -> int ->
    (int list, lcb) Hashtbl.t * int list list list * int list list * 
    (int * int) list list
*)

val print_lcb : lcb -> unit

val print_lcblst : int list list list -> unit

val get_abs_lst : int list -> int list

val update_lcb_ref_code : (int list, lcb) Hashtbl.t -> int list -> int -> unit


