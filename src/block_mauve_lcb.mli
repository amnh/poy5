type lcb = {
    seedNOlst : int list; (*each lcb is a list of mums*)
    range_lst : Block_mauve_mum.m_i list; 
    (* range_lst = [range of this lcb in seq0; range of this lcb in seq1].
    lcb_range we need is exactly the same truct as m_i,
    no need to create a new type here
    we sort it in update_lcbs_range_and_score, in decrease range size order*)
    ratio : float; (* [score/length] of seq in this lcb*)
    ref_code : int; (*just a code for bk/rearr purpose*)
    avg_range_len : int; (*average of range length from range_lst *)
    (*score between subsequence contained by range_lst during lcb building.
    * after we have lcb_tbl, huge lcb blocks are aligned seperately, score in
    * function search_inside_each_lcb is set to alignment cost of seq in this lcb*)
    score : int; 
    alignment : Sequence.s array;   
}

(** [print_lcb lcb] print one lcb*)
val print_lcb : lcb -> unit
(** [print_lcblst lcbkeylst] print out a list of lcb key, which is a intlstlst*)
val print_lcblst : int list list list -> unit

val search_inside_a_lcb : 
lcb -> Sequence.s -> Sequence.s -> 
int -> int -> 
(int, (int array, Block_mauve_mum.mum) BinaryTree.b_tree) Hashtbl.t ->
(int, (int * int * int) list * int array * int) Hashtbl.t ->
Cost_matrix.Two_D.m -> bool -> Sequence.s array * int


val create_lcb_tbl : int array array -> float -> int -> float -> int ->
Cost_matrix.Two_D.m -> bool ->
 (int list,lcb) Hashtbl.t * 
All_sets.IntegerListMap.key list list *
int list list * 
(int * int * int list * int) list list *
(int, (int array, Block_mauve_mum.mum) BinaryTree.b_tree) Hashtbl.t *
(int, (int * int * int) list * int array * int) Hashtbl.t

val update_lcb_ref_code : (int list, lcb) Hashtbl.t -> int list -> int -> unit
