


type m_i = {  (* the ith element M_i of a local mum *)
    sequence_NO : int ; (*which sequence is this m_i in *)
    left_end : int;  (*left end coordination of M_i*)
    right_end : int;
    orientation : int; (*1 or -1, 1 is positive, -1 is reverse complement*)
    }

type mum = { 
    seedNO : int; (*the # of seed that contruct this mum. start from
    only one seed, than extend*)
    mumseq : int list;
    positions : m_i list;  (* the positions in sequence this mum shows up*)
    mumkey : int;
    size : int; (* size of this mum, also the size of position list *)
    (* we don't extend seed with subset/superset for multi-sequence.
    left_superset : int list  ; (* seedNOlst of another mum*)
    right_superset : int list;
    left_subset : subset list  ; (* seedNOlst of another mum*)
    right_subset : subset list;
    *)
    neighborhood_lst: ( int * int * int * int * int ) list; 
    (* list of neighborhood (seqNO,j_seedNO,i_ori,j_ori,distance) 
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
    mumscore : int ;
}

type lcb = {
    seedNOlst : int list; (*each lcb is a list of mums,also the key to lcb_tbl*)
    range_lst : m_i list; (*lcb_range we need is exactly the same truct as m_i,
                           no need to create a new type here*)
    score : int; (*sum of individual MUM scores, not subsequence contained by
                  the range*)
    ratio : float; (* [score/length] of seq in this lcb*)
    ref_code : int; (*just a code for bk/rearr purpose*)
    avg_range_len : int; (*average length from range_lst*)
}


val create_lcb_tbl : int list list -> float -> int -> int ->
    (int list, lcb) Hashtbl.t * int list list list * int list list * 
    (int * int) list list

val print_mum : mum -> bool -> bool -> unit

val print_lcb : lcb -> unit

val print_int_list : int list -> unit

val print_int_lstlst : int list list -> unit

val print_lcblst : int list list list -> unit

val get_abs_lst : int list -> int list

val update_lcb_ref_code : (int list, lcb) Hashtbl.t -> int list -> int -> unit

val get_position_by_seqNO : m_i list -> int -> m_i

val get_lcb_key_by_range : int -> (int * int) -> (int list, lcb) Hashtbl.t 
                          -> int list * int
