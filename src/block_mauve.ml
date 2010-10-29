open Printf

(* W = weight, R = ratio *)
(*let sampleseq =
    [
        [1;2;4;8;1;2;4;8;1;2;4;8;4;4;4;4;]; 
        [4;4;4;1;2;4;8;1;2;4;8;1;2;4;8;4;]; 
        [8;8;8;1;2;4;8;1;2;4;8;1;2;4;8;4;]
    ]
let samplepattern = [|1;2;3|]
let seedlen = Array.length samplepattern 
let seedweight = samplepattern.(seedlen-1) 
*)

(*this is a cost matrix between ATGC and ATGC. default matrix used by Mauve *)
let hodx_matrix = [| 
(*           (*A*) (*T*) (*G*) (*C*)
  (*A*)  [|  91; -123;  -31; -114|];
  (*T*)  [|-123;   91; -114;  -31|];
  (*G*)  [| -31; -114;  100; -125|];
  (*C*)  [|-114;  -31; -125;  100|];
*)
           (*A*) (*C*) (*G*) (*T*)
  (*A*)  [|  91; -114;  -31; -123|];
  (*C*)  [|-114;  100; -125;  -31|];
  (*G*)  [| -31; -125;  100; -114|];
  (*T*)  [|-123;  -31; -114;   91|];

|]

let palindromic_spaced_seed_tbl =
    let res = Hashtbl.create 15 in
    let seed5 = [ [1;2;4;6;7];
        [1;4;5;6;9];
        [1;2;5;8;9]; ] in
    Hashtbl.add res 5 seed5;
    Hashtbl.add res 7
    [
        [1;2;5;7;9;12;13];
        [1;3;4;8;12;13;15];
        [1;2;4;8;12;14;15];
    ];
    Hashtbl.add res 9
    [
        [1;2;3;5;8;11;13;14;15];
        [1;2;3;6;9;12;15;16;17];
        [1;2;3;6;8;10;13;14;15];
    ];
    Hashtbl.add res 11
    [
        [1;2;3;5;7;10;13;15;17;18;19];
        [1;2;3;6;9;11;13;16;19;20;21];
        [1;2;3;4;7;9;11;14;15;16;17];
    ];
    (*Hashtbl.add res 13
    [
        1111**1**1*1*1**1**1111
        [1;2;3;4]
    ];
    add 13 15 19 21 later*)
    res

type m_i = {  (* the ith element M_i of a local mum *)
    sequence_NO : int ; (*which sequence is this m_i in *)
    left_end : int;  (*left end coordination of M_i*)
    right_end : int;
    orientation : int; (*1 or -1, 1 is positive, -1 is reverse complement*)
    }

type mum = { 
    seedNO : int; (*the # of seed that contruct this mum. start from
    only one seed, than extend*)
    positions : m_i list;  (* the positions in sequence this mum shows up*)
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
    seedNOlst : int list; (*each lcb is a list of mums*)
    range_lst : m_i list; (*lcb_range we need is exactly the same truct as m_i,
                           no need to create a new type here*)
    score : int; (*score between subsequence contained by range_lst*)
    ratio : float; (* [score/length] of seq in this lcb*)
    ref_code : int; (*just a code for bk/rearr purpose*)
}

let minimum_seed_num = ref 10
(*if the number of match we found is smaller than
    this, we should reduce the size of seedlen*)
let minimum_lcb_ratio = ref 30.0
let min_break_point_penalty = ref 4000 

(*lcb with lower weight should be considered removing*)
let minimum_lcb_weight = 2 
let init_tb_size = 50
let init_seed_size = 50
let max_gap_num = 0  (*the w *)
let seedNO_available_arr = ref (Array.make init_seed_size 1)

let bigger_int a b = 
    if (a>b) then a else b

let print_int_list inlist = 
    Printf.printf "[%!";
    List.iter (fun x -> Printf.printf "%d,%!" x) inlist;
    Printf.printf "] %!"

let print_int_list_to_file oc inlist =
    fprintf oc "[%!";
    List.iter (fun x -> fprintf oc "%d,%!" x) inlist;
    fprintf oc "] %!"

let print_int_lstlst inlstlst =
    List.iter (fun lst ->
        print_int_list lst;
        Printf.printf "\n%!";
    ) inlstlst

let print_lcblst lcbs =
    List.iter (fun lcblst ->
        List.iter (fun record -> print_int_list record) lcblst;
        Printf.printf "\n%!";
    ) lcbs

let print_int_list2 inlist  = 
    List.iter (fun x -> Printf.printf "%3i%!" x) inlist;
    Printf.printf " \n%!"

let print_int_list3 inlist = 
    let idx = ref 0 in
    List.iter (fun x -> Printf.printf "[%d]:%d; " !idx x; idx:= !idx +1; ) inlist;
    Printf.printf " \n%!"

let print_int_lstlst3 inlstlst =
    List.iter (fun lst -> print_int_list3 lst) inlstlst

let print_int_arr arr = 
    Array.iteri (fun idx item -> Printf.printf "[%d]:%d,%!"  idx item) arr;
    Printf.printf "\n%!"

let print_int_matrix m =
    Array.iter (fun arr ->
            print_int_arr arr ;
    ) m

let print_pos_list in_lst = 
    List.iter (fun pos ->
        Printf.printf 
        "[seqNO: %d, left end: %d, right end: %d, orientation: %d ] \n%!"
        pos.sequence_NO pos.left_end pos.right_end pos.orientation
    ) in_lst

let print_pos_list_to_file oc in_lst = 
    List.iter (fun pos ->
        fprintf oc 
        "[seqNO: %d, left end: %d, right end: %d, orientation: %d ] \n%!"
        pos.sequence_NO pos.left_end pos.right_end pos.orientation
    ) in_lst


let print_lcbs_range in_lstlst =
    Printf.printf "lcb range lst lst :\n%!";
    List.iter (fun lcb_range_lst ->
        Printf.printf "{ %!";
        List.iter (fun (leftend,rightend) ->
            Printf.printf "(%d,%d),%!" leftend rightend;
        )lcb_range_lst;
        Printf.printf " } \n%!";
    ) in_lstlst


let print_neighborhood_list in_lst = 
    List.iter (fun (seqNO, j,i_ori,j_ori,d) ->
        Printf.printf "{seqNO.%d, seed.#%d, ori:%d/%d,d:%d }\n%!" seqNO j i_ori
        j_ori d 
    ) in_lst

let print_mum in_SeedNO mum_tbl print_neighborhood print_unextendable =
    let in_mum = try (Hashtbl.find mum_tbl in_SeedNO) with 
    | Not_found -> failwith ("Print MUMs, not found err ") in
    if (in_mum.extendable = 0) ||
    ((in_mum.extendable !=0 )&&(print_unextendable=true) ) then begin
        Printf.printf "MUM# %d\n%!" in_SeedNO;
        Printf.printf "seedNO = %d,%!" in_mum.seedNO;
        Printf.printf " size:%d ,\n" in_mum.size;
        Printf.printf " subsuming_pointer = %d ,%!" in_mum.subsuming_pointer; 
        Printf.printf " extendable = %d \n%!" in_mum.extendable;
        Printf.printf "positions: \n%!";
        print_pos_list in_mum.positions;
        if print_neighborhood then begin
        Printf.printf " neighborhood list : \n%!";
        print_neighborhood_list in_mum.neighborhood_lst;
        end;
        Printf.printf "\n%!";
    end

(*print_pos_list2 print out one record of pos2seed_tbl, also the mum info*)
let print_pos_list2 seqNO pos in_lst mum_tbl printmum =
    Printf.printf "check record on (%d,%d) of pos2seed tbl(size=%d):\n%!"
    seqNO pos (List.length in_lst);
    List.iter (fun (seedNO,weight,orientation) ->
        Printf.printf "[seed: %d, weight: %d, ori: %d ] \n%!"
        seedNO weight orientation;
        if printmum then
        print_mum seedNO mum_tbl false true;
    ) in_lst;
    Printf.printf "end of checking \n%!"

let print_lcb x = 
    Printf.printf "LCB.%!";
    print_int_list x.seedNOlst; Printf.printf "\n%!";
    print_pos_list x.range_lst;
    Printf.printf "score = %d;%!" x.score;
    Printf.printf "ratio = %f;%!" x.ratio;
    Printf.printf "ref code = %d\n%!" x.ref_code

let print_lcb_to_file oc x = 
    fprintf oc "LCB.%!";
    print_int_list_to_file oc x.seedNOlst; fprintf oc "\n%!";
    print_pos_list_to_file oc x.range_lst;
    fprintf oc "score = %d\n%!" x.score;
    fprintf oc "ratio = %f\n%!" x.ratio 


let update_lcb_ref_code lcb_tbl key refcode = 
    let debug = false in
    let record = try (Hashtbl.find lcb_tbl key)
    with |Not_found -> failwith "not found updateing lcb refcode" in
    if debug then begin
        Printf.printf "update refcode %d to %!" refcode;
        print_int_list key;
        Printf.printf "\n%!";
    end;
    Hashtbl.remove lcb_tbl key;
    Hashtbl.add lcb_tbl key {record with ref_code = refcode}

let get_shortest_len seqlst =
    let initlen = (List.length (List.hd seqlst)) in
    List.fold_left (fun acc seq ->
        let newlen = List.length seq in
        if (newlen<acc) then newlen else acc
    ) initlen seqlst

(*return average number(float) of a int list *)
let get_avg_of_intlst in_lst =
    let acc = ref 0 in
    List.iter (fun size ->
        acc := size + !acc
    ) in_lst;
    (float !acc) /.(float (List.length in_lst))

let get_avg_of_floatlst in_lst =
    let acc = ref 0. in
    List.iter (fun size ->
        acc := size +. !acc
    ) in_lst;
    (!acc) /.(float (List.length in_lst))

let get_cons_id in_seq_size_lst num_of_mums = 
    let acc = ref 0.0 in
    List.iter (fun size ->
        acc := (1.0 /. (float size)) +. !acc 
    ) in_seq_size_lst;
    !acc *. (float num_of_mums) /. (float (List.length in_seq_size_lst))

let expand_arr add_len = 
    let add_arr = Array.make add_len 1 in
    let add_lst = Array.to_list add_arr in
    let ori_lst = Array.to_list !seedNO_available_arr in
    let new_lst = ori_lst@add_lst in
    seedNO_available_arr := Array.of_list new_lst

let return_a_seedNO idx =
    let seedNO_arr = !seedNO_available_arr in
    seedNO_arr.(idx)<-1;
    seedNO_available_arr := seedNO_arr

let get_a_seedNO () = 
    let found = ref (-1) in
    let idx = ref 1 in
    let seedNO_arr = !seedNO_available_arr in
    let arrlen = Array.length seedNO_arr in
    while (!found=(-1))&&(!idx<arrlen) do
        if (seedNO_arr.(!idx)=1) then begin 
            found := !idx;
            seedNO_arr.(!idx)<-0;
            seedNO_available_arr := seedNO_arr; 
        end;
        idx := !idx +1;
    done;
    if (!found = (arrlen-1)) then expand_arr init_seed_size;
    !found
     
let print_seedNO2seq_tbl in_tbl = 
    Hashtbl.iter (fun key record ->
        Printf.printf "seed#%d:%!" key;
        print_int_list record;
        Printf.printf "\n%!"
    )

let print_seed2postbl in_tbl =
    Hashtbl.iter (fun key record ->
    Printf.printf "seed# %d, { %!" key;
        List.iter (fun (seqNO,idx)->
        Printf.printf " (%d,%d);%!" seqNO idx
        )record;
        Printf.printf " }\n%!"
    ) in_tbl


let print_position2seedtbl position2seed_tbl = 
    Hashtbl.iter (fun (seqNO,pos) record_lst ->
        Printf.printf "pos.(%d,%d) ==> " seqNO pos;
        List.iter (fun (seed,seed_weight,dir) ->
        Printf.printf "seed.#%d,weight=%d,dir=%d; %!" seed seed_weight dir;
        ) record_lst;
        Printf.printf "\n%!";
    ) position2seed_tbl;
    Printf.printf "\n%!"

let print_seedNO2seqtbl in_tbl =
    Hashtbl.iter (fun seedNO seq ->
        Printf.printf "{ seed#%d:%!" seedNO;
        print_int_list seq;
        Printf.printf "}\n%!";
    ) in_tbl

(* if lst2 is a sub list of lst1, return 1, else 0 *)
let is_sub_list lst1 lst2 =
    let res = ref 1 in
    if (List.length lst2)>(List.length lst1) then 0
    else begin 
        List.iter ( fun y -> 
            if (List.mem y lst1) then () else res := 0
        ) lst2;
        !res
    end

(** Given an array [arr], a [looking_item] and a compare
* function [cmp_fun], this function return an integer as
* an  index of the [looking_item]. If the [looking_item] is not
* in the array [arr], return (-1) *)
let find_index arr looking_item cmp_fun = 
    let len = Array.length arr in 
    let rec find pos = 
        if pos = len then -1
        else
            if cmp_fun looking_item arr.(pos) = 0 then pos
            else find (pos + 1)
    in
    find 0

let get_abs_lst in_lst =
    List.sort (fun x y -> compare x y) (List.map (fun x -> abs x) in_lst)

(*NOTE: rev_comp just give us the complement, not reverse, use List.rev to do that*)
let rev_comp seqlst =
    List.map (fun x -> 
        if (x=1) then 2
        else if (x=2) then 1
        else if (x=4) then 8
        else if (x=8) then 4
        else 
            begin
            Printf.printf "reverse complement with x = %d,%!" x;
            failwith "we are dealing DNA sequence%!"
            end
    ) seqlst

let mark_unextendable_seed2pos_tbl seedNO ext_type seed2pos_tbl =
    let debug = false in
    if debug then 
        Printf.printf "mark seedNO#%d as UNextendab=%d in seed2pos tbl\n%!" 
        seedNO ext_type;
    let old_record = Hashtbl.find seed2pos_tbl seedNO in
    Hashtbl.remove seed2pos_tbl seedNO;
    Hashtbl.add seed2pos_tbl seedNO ( ((-1),ext_type)::old_record )

(*extendable= 0,1,2,3:
0. this seed is extendable (from). 
1 ~ 3 : not extendable (from),but might be extend to.
1. during the scan_seqlst, it means this mum shows up more then
once in a sequence, also shows up in every sequence,we should not 
extend the seed from this mum(but we can extend from other mum to this one),
2 .during the resolve_overlap, it means this mum does not
show up in one or more sequence. we don't extend to this one. 
3. during add to/remove from position2seed table, only one seed is kept as
extendable, others are marked as unextendable. but if we remove the
extendable one for some reason, we can upgrade one of this kind to be
extendable. 
*)
let mark_unextendable_mum_tbl seedNO ext_type mum_tbl = 
    let debug = false in
    if debug then Printf.printf "mark seedNO#%d as UNextendable=%d in mum_tbl\n%!" 
    seedNO ext_type;
    let oldmum = Hashtbl.find mum_tbl seedNO in
    Hashtbl.remove mum_tbl seedNO;
    Hashtbl.add mum_tbl seedNO {oldmum with extendable = ext_type}

let mark_extendable_mum_tbl seedNO mum_tbl = 
    let debug = false in
    if debug then Printf.printf "mark seedNO#%d as extendable in mum_tbl\n%!" seedNO;
    let oldmum = Hashtbl.find mum_tbl seedNO in
    Hashtbl.remove mum_tbl seedNO;
    Hashtbl.add mum_tbl seedNO {oldmum with extendable = 0}


let remove_seed seedNO pos2seedtbl_left poslst pos2seedtbl_right seedweight seedseqtbl
seedNOtbl =
    List.iter (fun (seqNO,pos) ->
        let pos = abs pos in
        if (Hashtbl.mem pos2seedtbl_left (seqNO,pos)) then
        Hashtbl.remove pos2seedtbl_left (seqNO,pos);
        if (Hashtbl.mem pos2seedtbl_right (seqNO,pos+seedweight-1)) then
        Hashtbl.remove pos2seedtbl_right (seqNO,pos+seedweight-1);
        if( Hashtbl.mem seedNOtbl seedNO) then begin 
        let seedseq =Hashtbl.find seedNOtbl seedNO
        in
        Hashtbl.remove seedseqtbl seedseq;
        Hashtbl.remove seedNOtbl seedNO;
        end 
        else (); (*we already remove seedNO and seedseq before *)
    ) poslst



(*when some seed shows up at most once in any sequence, doesn't mean it shows up in every
* sequence of input_seq*)
let match_at_most_once_in_each_seq poslst = 
    let sign = ref true in
    let hittbl = Hashtbl.create init_tb_size in 
    List.iter (fun record ->
        let seqNO = record.sequence_NO in
        if (Hashtbl.mem hittbl seqNO) then  sign := false
        else 
            Hashtbl.add hittbl seqNO 1 
    ) poslst;
    !sign

(*different poslst input type, this one is (seqNO,pos) list *)
let match_at_most_once_in_each_seq2 poslst = 
    let sign = ref true in
    let hittbl = Hashtbl.create init_tb_size in 
    List.iter (fun (seqNO,pos) ->
        if (Hashtbl.mem hittbl seqNO) then  sign := false
        else 
            Hashtbl.add hittbl seqNO 1 
    ) poslst;
    !sign

let match_at_least_once_in_every_seq poslst input_seqlst_size debug = 
    let hittbl = Hashtbl.create init_tb_size in 
    List.iter (fun record ->
        if debug then Printf.printf "%d,%d,%d\n%!" record.sequence_NO record.left_end
        record.right_end;
        let seqNO = record.sequence_NO in
        if (Hashtbl.mem hittbl seqNO) then ()
        else 
            Hashtbl.add hittbl seqNO 1 
    ) poslst;
    if (Hashtbl.length hittbl = input_seqlst_size) then true
    else false

let match_at_least_once_in_every_seq2 poslst input_seqlst_size  = 
    let hittbl = Hashtbl.create init_tb_size in 
    List.iter (fun (seqNO,pos) ->
        if (Hashtbl.mem hittbl seqNO) then ()
        else 
            Hashtbl.add hittbl seqNO 1 
    ) poslst;
    if (Hashtbl.length hittbl = input_seqlst_size) then true
    else false

let get_extendable_type poslst input_seqlst_size debug =
    let shows_up_in_every_seq = 
        match_at_least_once_in_every_seq poslst input_seqlst_size debug
    in
    let shows_up_at_most_once = match_at_most_once_in_each_seq poslst in
    if debug then Printf.printf "get extendable type, %b,%b\n%!" 
    shows_up_in_every_seq shows_up_at_most_once;
    let res:int =
        if (shows_up_in_every_seq=false) then 2
        else if shows_up_at_most_once then 0
        else 1
    in
    res
        

let remove_useless_seed seed_weight seq2seedNO_tbl seedNO2seq_tbl pos2seed_tbl_left pos2seed_tbl_right seed2pos_tbl =
    let debug = false in
    let input_seqlst_size = 2 in (*now we only have 2 input sequences*)
    if debug then print_seed2postbl seed2pos_tbl;
    Hashtbl.iter (fun seedNO poslst ->
        if debug then Printf.printf "clear up seed#.%d\n%!" seedNO;
        if (List.length poslst)<input_seqlst_size then begin 
            remove_seed seedNO pos2seed_tbl_left poslst pos2seed_tbl_right 
            seed_weight seq2seedNO_tbl seedNO2seq_tbl;
            Hashtbl.remove seed2pos_tbl seedNO;
            return_a_seedNO seedNO;
        end
        (*for match that shows up more than once in a sequence, we still need
        * them in seed extension, we can extend qualified match to them*)
        else begin
            let sign1 = match_at_most_once_in_each_seq2 poslst in
            let sign2 = match_at_least_once_in_every_seq2 poslst 
            input_seqlst_size  in
            let ext_type = 
                if (sign1)&&(sign2) then 0
                else if (sign1=false)&&sign2 then 1
                else if (sign2=false)&&sign1 then 2
                else 2
            in
            if (ext_type=2) then begin
                remove_seed seedNO pos2seed_tbl_left poslst pos2seed_tbl_right 
                seed_weight seq2seedNO_tbl seedNO2seq_tbl;
                Hashtbl.remove seed2pos_tbl seedNO;
                return_a_seedNO seedNO;
            end;
            if (ext_type=1) then
                mark_unextendable_seed2pos_tbl seedNO ext_type seed2pos_tbl;
        end; 
    ) seed2pos_tbl;
    Hashtbl.length seed2pos_tbl

(*get extendable record out of in_lst, in_lst is a recordlist of a position
* (seqNO,pos), from a pos2seed_tbl *)    
let get_extendable_recordlst in_lst mum_tbl = 
    List.filter (fun (seedNO,weight,ori) ->
        let mum = try (Hashtbl.find mum_tbl seedNO)
        with | Not_found -> failwith "get extendable subset not found 1" in
        (mum.extendable=0)
    ) in_lst 

(* this only allows at most one extendable seed start/end from on position*)
(* when there are two extendable seed, break the tie 
*  Note: remember to reload mum after this function, for its "extendable" sign
*  could be changed. *)
let add_to_pos2seed_tbl pos2seed_tbl seqNO pos seedNO seedweight orientation
mum_tbl = 
    let debug = false in
    if debug then Printf.printf "add seed#%d to pos (%d,%d)\n%!" seedNO seqNO
    pos;
    let add_to_table key newrecord tbl =
        let old_record = Hashtbl.find tbl key in
        let new_record = newrecord::old_record in
        Hashtbl.remove tbl key;
        Hashtbl.add tbl key new_record
    in
    let replace_to_table key newrecord tbl = 
        Hashtbl.remove tbl key;  Hashtbl.add tbl key newrecord
    in
    if (Hashtbl.mem pos2seed_tbl (seqNO,pos) ) then begin
        let new_mum = try (Hashtbl.find mum_tbl seedNO) 
        with | Not_found -> failwith "add_to_pos2seed_tbl not found 1" in
        if (new_mum.extendable!=0) then begin
            add_to_table (seqNO,pos) (seedNO,seedweight,orientation) pos2seed_tbl;
        end
        else begin
            let old_record = Hashtbl.find pos2seed_tbl (seqNO,pos) in
            if debug then print_pos_list2 seqNO pos old_record mum_tbl true;
            let old_record = get_extendable_recordlst old_record mum_tbl in
            assert((List.length old_record) <= 1);
            if debug then 
                Printf.printf "there is %d ext-record on this pos\n%!" 
                (List.length old_record) ;
            if ((List.length old_record) = 0 ) then 
                add_to_table (seqNO,pos) (seedNO,seedweight,orientation) pos2seed_tbl
            else begin
                let (old_seedNO,old_weight,old_ori) = List.hd old_record in
                if debug then 
                    Printf.printf "there was other mum on position (%d,%d)->(%d,%d,%d)\n%!"
                    seqNO pos old_seedNO old_weight old_ori;
                let old_mum = try (Hashtbl.find mum_tbl old_seedNO)
                with | Not_found -> failwith "add_to_pos2seed_tbl not found 1" in
                let old_multi = List.length old_mum.positions in
                let new_multi = List.length new_mum.positions in
                if (new_multi>old_multi)||((new_multi=old_multi)&&(seedweight>old_weight)) 
                then begin (*newseed become the only extendable seed on this pos*)
                    if debug then Printf.printf "newseed become the only extendable seed \n%!";
                    replace_to_table (seqNO,pos) [(seedNO,seedweight,orientation)] pos2seed_tbl;
                    mark_unextendable_mum_tbl old_seedNO 3 mum_tbl;
                    (*replace_to_table old_seedNO {old_mum with
                    * extendable=false} mum_tbl;*)
                end
                else begin (*add newseed to postbl, mark it as untendable*)
                    if debug then 
                    Printf.printf "seed#%d is add to and mark as untendatbl\n%!" seedNO;
                    add_to_table (seqNO,pos) (seedNO,seedweight,orientation) pos2seed_tbl;
                    replace_to_table seedNO {new_mum with extendable=3} mum_tbl;
                end;
            end;
        end (*end of if new mum is extendable*)
    end
    else begin
        if debug then Printf.printf "first record on this position\n%!";
        Hashtbl.add pos2seed_tbl (seqNO,pos) [(seedNO,seedweight,orientation)];
    end

(*if we remove the extendable seed from a position, in case there are
* other seeds on current pos, we can upgrade one of them to be extendable. *)
let remove_from_pos2seed_tbl pos2seed_tbl seqNO old_pos seedNO mum_tbl = 
    let debug = false in
    if debug then Printf.printf "remove seed#%d from (%d,%d) \n%!" seedNO seqNO old_pos;
    let old_record = try (Hashtbl.find pos2seed_tbl (seqNO,old_pos)) with
    | Not_found -> failwith ("remove pos2seedtbl, not found 1") in
    let new_record = List.filter (fun (x,_,_) -> x<>seedNO) old_record in
    let mum_i = try (Hashtbl.find mum_tbl seedNO) with
    | Not_found -> failwith ("remove pos2seedtbl, not found 2") in
    if (mum_i.extendable=0) then begin
        if debug then 
            Printf.printf "seed#%d is the only extendable one on this pos\n%!" seedNO;
        let sorted_new_record = List.sort (fun (_,j_weight,_) (_,k_weight,_) ->
            compare j_weight k_weight
        ) new_record in
        let seedNO_to_upgrade,_,_ = 
        try (List.find (fun (j_seedNO,_,_) ->
            let j_mum = try (Hashtbl.find mum_tbl seedNO) with
            | Not_found -> failwith ("remove pos2seedtbl, not found 3") in
            (j_mum.extendable = 3)
        ) sorted_new_record)
        with | Not_found -> (-1,-1,-1)
        in
        if (seedNO_to_upgrade != (-1)) then begin
            if debug then 
                Printf.printf "seed# %d become the extendable one\n%!" seedNO_to_upgrade;
            mark_extendable_mum_tbl seedNO_to_upgrade mum_tbl;
        end;
    end;
    Hashtbl.remove pos2seed_tbl (seqNO,old_pos);
    Hashtbl.add pos2seed_tbl (seqNO,old_pos) new_record

let modify_record_in_pos2seed_tbl pos2seed_tbl seqNO pos seedNO weight ori =
    let old_record = try (Hashtbl.find pos2seed_tbl (seqNO,pos)) with
    | Not_found -> failwith ("remove pos2seedtbl, not found 1") in
    let unchanged = List.filter (fun (x,_,_) -> x<>seedNO) old_record in
    let new_record = (seedNO,weight,ori)::unchanged in
    Hashtbl.add pos2seed_tbl (seqNO,pos) new_record


let add_to_seed2pos_tbl seed2pos_tbl seqNO pos seedNO seedweight orientation =
    if( Hashtbl.mem seed2pos_tbl seedNO ) then begin
        let oldposlst = Hashtbl.find seed2pos_tbl seedNO in
        Hashtbl.remove seed2pos_tbl seedNO;
        Hashtbl.add seed2pos_tbl seedNO ((seqNO,pos)::oldposlst);
    end
    else
        Hashtbl.add seed2pos_tbl seedNO [(seqNO,pos)]
                    
let remove_from_seed2pos_tbl seed2pos_tbl seqNO pos seedNO = 
    assert(Hashtbl.mem seed2pos_tbl seedNO);
    let oldposlst = Hashtbl.find seed2pos_tbl seedNO in
    let newposlst = List.filter (fun item ->
        let (old_seqNO,old_pos) = item in
        if (old_seqNO=seqNO)&&(old_pos=pos) then false else true
    ) oldposlst in
    Hashtbl.remove seed2pos_tbl seedNO;
    Hashtbl.add seed2pos_tbl seedNO newposlst


let scan_seqlst (inseqlst:int list list) patternarr 
seq2seedNO_tbl seedNO2seq_tbl 
position2seed_tbl_left position2seed_tbl_right 
seed2position_tbl debug =
    let inseqarr = Array.of_list (List.map (fun lst -> Array.of_list lst
    ) inseqlst) in
    let seedweight = patternarr.((Array.length patternarr)-1) in
    if debug then Printf.printf "seed weight = %d\n" seedweight;
    Array.iteri (
        fun sequenceNO seq ->
            let seqlen = Array.length seq in
            if debug then Printf.printf "seed len = %d\n%!" seqlen;
            Array.iteri (fun idx item ->
            if ( (idx+seedweight-1) < seqlen) then begin
                let subseqlst = Array.fold_right 
                     (fun p acc -> seq.(idx+p-1)::acc) patternarr [] in
                if (Hashtbl.mem seq2seedNO_tbl subseqlst) then begin
                    let seedNO = Hashtbl.find seq2seedNO_tbl subseqlst in
                    Hashtbl.add position2seed_tbl_left (sequenceNO,idx)
                    [(seedNO,seedweight,1)]; (*1 is orientation*)
                    let oldposlst = 
                        try (Hashtbl.find seed2position_tbl seedNO) with
                        | Not_found -> failwith ("scan seq Not found 1")
                    in
                    Hashtbl.remove seed2position_tbl seedNO;
                    Hashtbl.add seed2position_tbl seedNO ((sequenceNO,idx)::oldposlst);
                end
                else
                    begin
                    let seedNO = get_a_seedNO () in
                    Hashtbl.add seq2seedNO_tbl subseqlst seedNO;
                    Hashtbl.add position2seed_tbl_left (sequenceNO,idx)
                    [(seedNO,seedweight,1)];
                    Hashtbl.add seed2position_tbl seedNO [(sequenceNO,idx)];
                    end
            end
            else ()
        ) seq
    ) inseqarr;
    if debug then Printf.printf "check the reverse complement\n%!";
    let rev_inseqlst = (List.map (fun lst -> List.rev lst) inseqlst) in
    let revseqarr = Array.of_list (List.map (fun lst -> Array.of_list lst
    ) rev_inseqlst) in
    Array.iteri (
        fun sequenceNO revseq ->
        let seqlen = Array.length revseq in
        Array.iteri (
            fun idx item ->
                let realidx = seqlen-idx-1 in
                if ( (idx+seedweight-1) < seqlen) then begin
                    let subseqlst = Array.fold_right 
                         (fun p acc -> revseq.(idx+p-1)::acc) patternarr [] in
                    let rev_subseqlst = rev_comp subseqlst in
                    if (Hashtbl.mem seq2seedNO_tbl rev_subseqlst) then
                        let new_seedNO = Hashtbl.find seq2seedNO_tbl rev_subseqlst in
                        let old_record_lst = 
                            Hashtbl.find position2seed_tbl_left
                            (sequenceNO,(realidx-seedweight+1)) in
                        assert((List.length old_record_lst)=1);
                        let old_seedNO,seed_weight,dir = List.hd old_record_lst in 
                        (*if our seed is not DNA-palindromic, old_seedNO = new_seedNO will never happen*)
                        if (old_seedNO<>(-1))&&(old_seedNO<>new_seedNO) then begin 
                            Hashtbl.remove seq2seedNO_tbl (List.rev subseqlst);
                            let poslst = try (Hashtbl.find seed2position_tbl
                            old_seedNO) with |Not_found -> failwith ("scan seq, Not found 5") in
                            List.iter (fun (seqNO,pos) ->
                                Hashtbl.remove position2seed_tbl_left (seqNO,pos);
                                Hashtbl.add position2seed_tbl_left (seqNO,pos)
                                [(new_seedNO,seed_weight,(-1))];
                                if (Hashtbl.mem seed2position_tbl new_seedNO) then
                                begin
                                let oldposlst = Hashtbl.find seed2position_tbl (new_seedNO) in
                                Hashtbl.remove seed2position_tbl new_seedNO;
                                Hashtbl.add seed2position_tbl new_seedNO ((seqNO,-pos)::oldposlst);
                                end
                                else Hashtbl.add seed2position_tbl new_seedNO [(seqNO,-pos)]
                            ) poslst;
                            Hashtbl.remove seed2position_tbl old_seedNO;
                            return_a_seedNO old_seedNO;
                        end
                    else ();
                end
                else ()
        ) revseq;
    ) revseqarr;
    Hashtbl.iter (fun (seqNO,pos) recordlist ->
        Hashtbl.add position2seed_tbl_right (seqNO,(pos+seedweight-1)) recordlist;
    )position2seed_tbl_left;
    Hashtbl.iter (fun seedseq seedNO ->
        Hashtbl.add seedNO2seq_tbl seedNO seedseq;
    ) seq2seedNO_tbl;
    if debug then print_seedNO2seqtbl seedNO2seq_tbl;
    if debug then 
        Printf.printf "end of scan sequence, return seedweight = %d \n%!" seedweight;
    seedweight



let build_seed_and_position_tbl inseqlst init_seedlen min_seednumber 
seq2seedNO_tbl seedNO2seq_tbl pos2seed_tbl_left pos2seed_tbl_right
seed2pos_tbl debug
=
    if debug then Printf.printf "build_seed_and_position_tbl, init_seedlen = %d\n%! "
    init_seedlen;
    let seed_len = ref init_seedlen in
    let seedlst = 
        try (Hashtbl.find palindromic_spaced_seed_tbl init_seedlen)
        with | Not_found -> failwith ("cannot find seed pattern1")
    in
    let seedpattern = List.hd seedlst in
    if debug then begin 
        Printf.printf "seed pattern list = %!";
        print_int_list (seedpattern); Printf.printf "\n%!";
    end;
    let patternarr = Array.of_list seedpattern in
    let seedweight = scan_seqlst inseqlst patternarr seq2seedNO_tbl seedNO2seq_tbl 
        pos2seed_tbl_left pos2seed_tbl_right seed2pos_tbl false in
    if debug then
        Printf.printf "seed number before remove useless seed = %d\n%!"
        (Hashtbl.length seedNO2seq_tbl);
    let seed_number = ref
    (remove_useless_seed seedweight seq2seedNO_tbl seedNO2seq_tbl 
    pos2seed_tbl_left pos2seed_tbl_right seed2pos_tbl) 
    in
    if debug then 
        Printf.printf "seed number after remove useless seed :%d\n%!" !seed_number;
    while (!seed_number<min_seednumber)&&(!seed_len>5 ) do 
        seed_len := !seed_len - 2;
        let seedlst = 
        try (Hashtbl.find palindromic_spaced_seed_tbl !seed_len)
        with | Not_found -> failwith ("cannot find seed pattern 2")
        in
        let seedpattern = List.hd seedlst in
        Printf.printf "reduce seed weight, new seed pattern list = %!";
        print_int_list (seedpattern); Printf.printf "\n%!";
        let patternarr = Array.of_list seedpattern in
        let seedweight = scan_seqlst inseqlst patternarr 
        seq2seedNO_tbl seedNO2seq_tbl
        pos2seed_tbl_left pos2seed_tbl_right seed2pos_tbl false in
        seed_number := (remove_useless_seed seedweight seq2seedNO_tbl seedNO2seq_tbl 
        pos2seed_tbl_left pos2seed_tbl_right seed2pos_tbl);
    done;
    if debug then begin
    (*Printf.printf "check position2seed tbl left and right\n%!";
    print_position2seedtbl pos2seed_tbl_left;
    print_position2seedtbl pos2seed_tbl_right;
    Printf.printf "check seed2position tbl :\n%!";
    print_seed2postbl seed2pos_tbl;
    Printf.printf "check seedNO2seq_tbl :\n%!";
    print_seedNO2seqtbl seedNO2seq_tbl; *)
    Printf.printf "end of build_seed_and_position_tbl,seed_len=%d\n%!" !seed_len;
    end;
    seedweight

let sort_positions_by_seqNO in_lst =
    List.sort (fun x y ->
       compare x.sequence_NO y.sequence_NO 
    ) in_lst

(*this works if there is at most one match in that sequence*)
let get_position_by_seqNO in_lst seqNO = 
    let index = find_index (Array.of_list in_lst) seqNO 
    (fun looking_item item ->
        compare item.sequence_NO looking_item ) in
    let len = List.length in_lst in
    if (index<0)||(index>=len) then begin
        Printf.printf "try to get seqNO = %d from \n%!" seqNO;
        print_pos_list in_lst;
    end;
    assert(index>=0);
    assert(index<len);
    List.nth in_lst index

let get_largest_range_of_mi_lst in_lst = 
    let largest_range = ref 0 in
    let acc_idx = ref 0 
    and idx = ref 0 in
    List.iter (fun record ->
        let range = record.right_end - record.left_end in
        if (!largest_range<range) then begin
            acc_idx := !idx;
            largest_range := range;
        end;
        idx := !idx + 1;
    ) in_lst;
    let res = List.nth in_lst !acc_idx in
    res.right_end - res.left_end

(*
let match_only_once_in_every_seq2 mi_lst = 
    let sign = ref true in
    let hittbl = Hashtbl.create init_tb_size in 
    List.iter (fun mi ->
        let seqNO = mi.sequence_NO in
        if (Hashtbl.mem hittbl seqNO) then  sign := false
        else Hashtbl.add hittbl seqNO 1 
    ) mi_lst;
    !sign
*)

(*this works for adding new position to a pos_list, or update content of a
* position in the pos_list*)
let update_positions_by_seqNO seedNO seqNO new_pos mum_tbl = 
    (*this works only for input two sequences*)
    let debug =  false in
    let input_seqlst_size = 2 in
    let mi = try (Hashtbl.find mum_tbl seedNO) with
    | Not_found -> failwith "update positions by seqNO not found 1" in 
    let old_poslen = List.length mi.positions in
    let positions_unchanged = List.filter (fun record ->
        if(record.sequence_NO <> seqNO) then true
        else false
    ) mi.positions in 
    let new_positions = new_pos :: positions_unchanged in
    let new_poslen = List.length new_positions in
    let ext_sign = 
        if (new_poslen = old_poslen) then mi.extendable
        else
            get_extendable_type new_positions input_seqlst_size debug 
    in
    let new_mi = {mi with positions = new_positions; size = List.length
    new_positions; extendable = ext_sign } in
    Hashtbl.remove mum_tbl mi.seedNO;
    Hashtbl.add mum_tbl mi.seedNO new_mi;
    if debug then print_mum mi.seedNO mum_tbl false true

let remove_positions_by_seqNO seedNO seqNO mum_tbl = 
    let mi = try (Hashtbl.find mum_tbl seedNO) with
    | Not_found -> failwith "update positions by seqNO not found 1" in 
    let positions_unchanged = List.filter (fun record ->
        if(record.sequence_NO <> seqNO) then true
        else false
    ) mi.positions in 
    let new_positions = positions_unchanged in
    let new_mi = {mi with positions = new_positions; size = List.length
    new_positions } in
    Hashtbl.remove mum_tbl mi.seedNO;
    Hashtbl.add mum_tbl mi.seedNO new_mi


(*remove_seed2 remove j_seedNO frompos2seedtbl_left/right, also mark j_seedNO as
* unextendable -- maybe we should remove j_seedNO from mum_tbl completely. also
* remove k_seedNO, we create k_seedNO during the trimp of overlap-mums*)
let remove_seed2 j_seedNO seqNO pos2seed_tbl_left pos2seed_tbl_right jleft jright
mum_tbl k_seedNO = 
    let debug = false in
    if debug then Printf.printf "remove match on seq#%d of mum#%d\n%!" seqNO j_seedNO;
    remove_from_pos2seed_tbl pos2seed_tbl_left seqNO
    jleft j_seedNO mum_tbl;
    remove_from_pos2seed_tbl pos2seed_tbl_right seqNO
    jright j_seedNO mum_tbl ;
    remove_positions_by_seqNO j_seedNO seqNO mum_tbl;
    mark_unextendable_mum_tbl j_seedNO 2 mum_tbl;
    return_a_seedNO k_seedNO


(*our idx starts from 0, while mauve starts from 1*)
let format_output mum_tbl = 
    let res = ref [[]] in
    Hashtbl.iter (fun seedNO mum ->
        if (mum.subsuming_pointer=(-1))&&(mum.extendable=0) then begin
            let poslst = sort_positions_by_seqNO mum.positions in
            let poshd = List.hd poslst in
            let size = poshd.right_end - poshd.left_end + 1 in
            let tmp = ref [] in
            List.iter (fun mi ->
                let add = if (mi.orientation = 1) then mi.left_end+1
                else -(mi.left_end+1) in
                tmp := !tmp@[add]
            ) poslst;
            tmp := size::!tmp;
            res := !tmp::!res;
        end 
        else ();
    ) mum_tbl;
    !res

let is_subsumed_by j_seed i_seed mum_tbl =
    let j_mum = try (Hashtbl.find mum_tbl j_seed) with 
    | Not_found -> failwith ("is subsumed , not found err ") in
    if (j_mum.subsuming_pointer=i_seed) then true
    else false


let is_subsuming i_seedNO j_seedNO mum_tbl =
    let i_mum = try (Hashtbl.find mum_tbl i_seedNO) with 
            | Not_found -> failwith ("is subsuming, not found err 1") in
    let j_mum = try (Hashtbl.find mum_tbl j_seedNO) with 
            | Not_found -> failwith ("is subsuming, not found err 2") in
    let i_pos_lst = i_mum.positions in
    let j_pos_lst = j_mum.positions in
    if (List.length i_pos_lst)>=(List.length j_pos_lst) then begin
        let sign = ref true in
        List.iter (fun j_pos ->
            let i_idx = find_index (Array.of_list i_pos_lst) j_pos 
            (fun j_pos i_pos ->
                let j_le = j_pos.left_end and j_re = j_pos.right_end in
                let i_le = i_pos.left_end and i_re = i_pos.right_end in
                if (j_le>=i_le)&&(j_re<=i_re) then 0
                else (-1) )
            in
            if (i_idx<0) then sign := false
        ) j_pos_lst;
        !sign
    end
    else false 



let sort_tbl_by_size in_tbl =
    let p_queue = ref [] in
    Hashtbl.iter (fun seedNO mum ->
        let size = mum.size in
        p_queue := (seedNO, size)::!p_queue;
    ) in_tbl;
    let p_queue = List.sort (fun (_,sizex) (_,sizey) ->
        compare sizex sizey ) !p_queue 
    in
    List.rev p_queue


let sort_neighborhood_lst_by_dis in_lst = 
    List.sort(fun x y ->
        let (_,_,_,_,dx) = x and (_,_,_,_,dy) = y in
        compare (abs dx) (abs dy)
    ) in_lst

(*in_lst is neighborhood list of one j_seedNO to i_seedNO *)
let group_neighborhood_lst_by_dis in_lst debug =
    if debug then Printf.printf "group neighborhood lst\n%!";
    let in_lst = sort_neighborhood_lst_by_dis in_lst in
    let first_record = List.hd in_lst in
    let (_,_,first_i_ori,first_j_ori,first_distance) = first_record in
    let res = ref [] in
    if debug then Printf.printf "start with (%d,%d,%d)\n%!" first_i_ori first_j_ori
    first_distance;
    let _,_,last_dis, last_acc = List.fold_right (fun record (pre_iori,pre_jori,previous_dis, acc) ->
        let (_,_,i_ori,j_ori,distance) = record in
        if debug then Printf.printf "work on (%d,%d,%d) with pre=(%d,%d,%d)\n%!" i_ori j_ori distance pre_iori pre_jori previous_dis;
        if (pre_iori*pre_jori=i_ori*j_ori)&&(pre_iori*previous_dis=i_ori*distance) then 
            (i_ori,j_ori,distance,record::acc)
        else begin
            if debug then Printf.printf "different distance,wrap up\n%!";
            res := (previous_dis,acc) :: !res;
            (i_ori,j_ori,distance, [record])
        end
    ) (List.tl in_lst) (first_i_ori,first_j_ori,first_distance,[first_record]) in
    res := (last_dis,last_acc) :: !res;
    (*make sure the first record of res is the one group with smallest distance*)
    List.rev !res

(* this allows more than one seed starts/ends from one positon*)
(* we call this to expand the position to left or/and right*)
let update_position2seedtbl seqNO old_pos new_pos seed_weight ori seedNO positiontbl =
    let old_record = try (Hashtbl.find positiontbl (seqNO,old_pos)) with
    | Not_found -> failwith ("update pos2seedtbl, not found 1") in
    let new_record = List.filter (fun (x,_,_) -> x<>seedNO) old_record in
    Hashtbl.remove positiontbl (seqNO,old_pos);
    Hashtbl.add positiontbl (seqNO,old_pos) new_record;
    if (Hashtbl.mem positiontbl (seqNO,new_pos)) then begin
        (*set the subsuming pointer here?*)
        let old_record = Hashtbl.find positiontbl (seqNO,new_pos) in
        let new_record = (seedNO,seed_weight,ori)::old_record in
        Hashtbl.remove positiontbl (seqNO,new_pos);
        Hashtbl.add positiontbl (seqNO,new_pos) new_record;
    end
    else Hashtbl.add positiontbl (seqNO,new_pos) [(seedNO,seed_weight,ori)]


let deal_with_tandem i_seedNO repeat_units_lst mum_tbl =
    Printf.printf "figure out how to do this...\n%!"

let categorize_neighboorhood_lst in_list i_seedNO i_size mum_tbl =
    let debug = false in
    if debug then begin 
        Printf.printf " categorize_neighboorhood_lst, check in_list : \n%!";
        print_neighborhood_list in_list;
    end;
    (*let chainable_lst = ref [] in  we will need this later*)
    let overlap_lst = ref [] in
    (*we don't need to worry about j_mum.seqNO = i_mum.seqNO
    * mum shows up only onece in a sequence,*)
    let tbj = Hashtbl.create init_tb_size in
    List.iter (fun item ->
        let (seqNO,j,i_ori,j_ori,d) = item in
        if (is_subsumed_by j i_seedNO mum_tbl) then ()
        else begin 
            if (Hashtbl.mem tbj j) then begin
                let old_record = Hashtbl.find tbj j in
                Hashtbl.remove tbj j;
                Hashtbl.add tbj j (item::old_record);
            end
            else Hashtbl.add tbj j [item];
        end
    ) in_list;
    Hashtbl.iter (fun jseedNO record_list  ->
        (*we only extend to overlap seed now, and only to those have same
        * distance to current mum in each sequence. *)
        (*group list is neighborhood of iseedNO with jseedNO, group them by
        * distance, the first item in group_list is the best seed-extension
        * candidates regarding jseedNO*)
        let group_list = group_neighborhood_lst_by_dis record_list false in
        if debug then Printf.printf "check grouped list:\n%!";
        if debug then List.iter (fun (dis,sublst) ->
            Printf.printf "with dis=%d,we have \n%!" dis;
            print_neighborhood_list sublst;
        ) group_list;
        (*we only care about the group with min distance*)
        let find_one = ref false in
        List.iter (fun (dis,sublst)->
            if ((!find_one=false)&&((List.length sublst)=i_size)) then 
                begin
                overlap_lst := (dis,sublst)::(!overlap_lst);
                find_one := true;
                end
        ) group_list ;
        (* this part is for chainable list, keep it here .
        let j_mum = try (Hashtbl.find mum_tbl jseedNO) with 
        | Not_found -> failwith ("categorize neighboorhood, not found err 1") in
        let j_size = j_mum.size in
        let record_lst1,record_lst2 = 
            List.fold_right (fun record (acc1,acc2) ->
                let (_,_,i_ori,j_ori,_) = record in
                if (i_ori,j_ori=1) then record::acc1,acc2
                else acc1, record::acc2
            ) record_list ([],[])
        in
        let r_size1 = List.length record_lst1 
        and r_size2 = List.length record_lst2 in
        if debug then 
        Printf.printf "Categorize on i & j = %d & %d, r_size1/2=%d/%d, i_size=%d,j_size=%d\n%!" i_seedNO j_mum.seedNO r_size1 r_size2 i_size j_size;
        if (j_size=r_size1)&&(i_size=r_size1) then
            chainable_lst := record_lst1::(!chainable_lst) 
        else if (j_size=r_size2)&&(i_size=r_size2) then
            chainable_lst := record_lst2::(!chainable_lst)
        else ();
        *)
    )tbj;
     let res = List.sort (fun (disa,lsta) (disb,lstb) ->
         compare (abs disa) (abs disb)
     ) !overlap_lst in
     res

(*update neighborhood of MUM#.seedNO, if MUM#seedNO is extendable (=0 or 3),
* other matches subsumed by seedNO will not be included in neighborhood. *)     
let update_neighborhood seedNO mum_tbl pos_tbl_left pos_tbl_right debug = 
    if debug then 
    Printf.printf "update neighborhood on seedNO# %d\n%!" seedNO;
    (*print_position2seedtbl pos_tbl_left; print_position2seedtbl pos_tbl_right; *)
    let i_mum = try (Hashtbl.find mum_tbl seedNO) with 
    | Not_found -> failwith ("build_neighborhood, not found err 1") in
    if debug then print_mum seedNO mum_tbl true true;
    (*what if mi or mj is subsuming by other mum ?*)
    (*we don't care neighborhood of unextendable (from) seed, as long as we
    * don't do the chainable list, this should be fine *)
    if (i_mum.extendable=0||i_mum.extendable=3) then begin
        let neighborhood_lst = ref [] in (*start from empty list*)
        let m_i_lst = i_mum.positions in
        let i_idx = ref 0 and d = max_gap_num in
        List.iter (fun mi ->
            let seqNO, i_left_end,i_right_end,i_ori = 
                mi.sequence_NO, mi.left_end, mi.right_end, mi.orientation in
            if debug then Printf.printf "work on seqNO=%d,left=%d,right=%d\n%!"
            seqNO i_left_end i_right_end ;
            let stop = ref 0 in
            (*"next to" also count as "overlap", we add it into neighborhood, too*)
            let pos = ref (i_left_end-d-1) in
            while (!stop<>1)&&( (!pos) <= (i_right_end+d+1)) do
                if ((!pos<>i_left_end)&&(Hashtbl.mem pos_tbl_right (seqNO,!pos)))then begin
                    let record_list = Hashtbl.find pos_tbl_right (seqNO,!pos) in
                    List.iter (fun (j_seedNO,seed_weight,j_ori) ->
                        if (Hashtbl.mem mum_tbl j_seedNO)=false then
                            Printf.printf "cannot find j_seedNO:%d%!" j_seedNO;
                        let j_mum = try (Hashtbl.find mum_tbl j_seedNO) with 
                        | Not_found -> failwith ("build_neighborhood, not found 2") in
                        if debug then begin 
                            Printf.printf "right tbl,pos=(%d,%d),j_seedNO=\n%!" seqNO !pos;
                            print_mum j_seedNO mum_tbl false true;
                        end;
                        (*don't care about j_mum subsumed by other mum*)
                        if (j_mum.subsuming_pointer<>(-1)) then () 
                        (*stop adding neighborhoodlst if we reach another match
                        * of mum_i in the same seq.
                        * Note. current i_left/right_end will show up in pos_tbl,
                        * since mauve also consider "next to" as "overlap", we need
                        * to check the position left to i_left_end, and right to
                        * i_right_end, don't stop when we reach i_left/right_end*)
                        else if (j_seedNO<>seedNO) then begin
                            let j_idx = 
                                find_index (Array.of_list j_mum.positions)
                                !pos (fun p item ->
                                    if (item.right_end=p) then 0 else (-1) )
                            in
                            if debug then print_mum j_seedNO mum_tbl false true;
                            if (j_idx<0) then
                                Printf.printf "j_idx<0,j=%d,pos=%d\n%!" j_seedNO !pos;
                            assert(j_idx>=0);
                            let j_pos = List.nth j_mum.positions j_idx in
                            let j_left_end = j_pos.left_end in
                            (*assert( (i_left_end-j_left_end)<>0 );it happens*)
                            let newneighbor =
                                (seqNO,j_seedNO,i_ori,j_ori,i_left_end-j_left_end)
                            in
                            if (j_idx>=0) then
                                if ((List.mem newneighbor !neighborhood_lst)=false) then 
                                neighborhood_lst := newneighbor::!neighborhood_lst
                                else ()
                            else failwith ("build neighborhood,idx not found 4");
                        end
                        else if ((j_seedNO=seedNO))&&(!pos=i_right_end) then ()
                        else stop := 1;
                    )record_list;
                end;
                if ((!pos<>i_left_end)&&(Hashtbl.mem pos_tbl_left (seqNO,!pos))) 
                then begin
                    let record_list = Hashtbl.find pos_tbl_left (seqNO,!pos) in
                    List.iter (fun (j_seedNO,seed_weight,j_ori) ->
                        if debug then begin 
                            Printf.printf "left tbl,pos=(%d,%d),j_seedNO=\n%!" seqNO !pos;
                            print_mum j_seedNO mum_tbl false true;
                        end;
                        let j_mum = try (Hashtbl.find mum_tbl j_seedNO) with 
                        | Not_found -> failwith ("build_neighborhood, not found 3") in
                        if (j_mum.subsuming_pointer<>(-1)) then ()
                        else if (j_seedNO<>seedNO) then begin 
                            if debug then print_mum j_seedNO mum_tbl false true;
                            let j_idx = 
                                find_index (Array.of_list j_mum.positions)
                                !pos (fun p item ->
                                    if (item.left_end=p) then 0 else (-1) )
                            in
                            assert(j_idx>=0);
                            assert( (i_left_end-(!pos))<>0);
                            let newneighbor =
                                (seqNO,j_seedNO,i_ori,j_ori,i_left_end-(!pos)) in
                            if (j_idx>=0) then 
                                if ((List.mem newneighbor !neighborhood_lst)=false) then
                                neighborhood_lst := newneighbor ::!neighborhood_lst 
                                else ()
                            else 
                                failwith ("build neighborhood,idx not found 5");
                        end
                        else if ((j_seedNO=seedNO))&&(!pos=i_left_end) then ()
                        else stop := 1;
                    )record_list;
                end;
                pos := !pos + 1;
            done; (* for pos = (i_left_end-(!d)) to (i_right_end+(!d)) do*)
            i_idx := !i_idx + 1;
        ) m_i_lst;
        let res = sort_neighborhood_lst_by_dis !neighborhood_lst in
        let new_i_mum = {i_mum with neighborhood_lst = res } in
        Hashtbl.remove mum_tbl seedNO;
        Hashtbl.add mum_tbl seedNO new_i_mum;
        if debug then begin
            Printf.printf "end of update_neighborhood check mum again\n%!";
            print_mum seedNO mum_tbl true true;
        end
    end
    else () (*if mi.extendable = 0 or 3*)
  (*  end *)
        

(*dir: direction , left = 1, right = (-1)*)
let build_neighborhood p_queue mum_tbl pos_tbl_left pos_tbl_right debug =
    if debug then Printf.printf "build neighborhood start \n%!";
    let idx_p_queue = ref 0 and size_p_queue = List.length p_queue in 
    while ( (!idx_p_queue<size_p_queue) ) do
        let (seedNO, size) = List.nth p_queue !idx_p_queue in
        update_neighborhood seedNO mum_tbl pos_tbl_left pos_tbl_right false;
        let i_mum = Hashtbl.find mum_tbl seedNO in
        let nei_list = i_mum.neighborhood_lst in
        Hashtbl.remove mum_tbl seedNO;
        Hashtbl.add mum_tbl seedNO {i_mum with neighborhood_lst = nei_list};
        idx_p_queue := !idx_p_queue + 1;
    done;
    if debug then Printf.printf "build neighborhood done\n%!"

(*make a chain chain i_seed and j_seed together. chainable_lst is a sublist of
* i_seed's neighborhood, each of them is a record of j_seed. 
* Note: number of j_seed in chainable_lst might be larger then number of
* matches of i_seed, for j_seed could show up more than once in a sequence -- in
* that case, j_seed's extendable is 1.*)
let make_a_chain i_seedNO chainable_lst mum_tbl pos2seed_tbl_left pos2seed_tbl_right =
    let debug = false in
    let (distance,nei_list) = List.hd chainable_lst in
    let (_,j_seedNO,_,_,_) = List.hd nei_list in
    let j_mum = try (Hashtbl.find mum_tbl j_seedNO) with 
    | Not_found -> failwith ("in make a chain, not found err 2") in
    let j_pos_lst = j_mum.positions in
    if debug then Printf.printf "make a chain start \n%!";
    if debug then print_neighborhood_list nei_list;
    let i_mum = try (Hashtbl.find mum_tbl i_seedNO) with 
    | Not_found -> failwith ("in make a chain, not found err 1") in
    let i_pos_lst = i_mum.positions in
    if (j_mum.subsuming_pointer = i_seedNO) then ()
    else begin
        List.iter (fun (seqNO,j_seedNO,i_ori,j_ori,d) ->
            if debug then begin
            Printf.printf "chain seed#%d and #%d together,d=%d,i_ori,j_ori=%d/%d \n%!"
            i_seedNO j_seedNO d i_ori j_ori;
            Printf.printf "before chainning:\n%!";
            print_mum i_seedNO mum_tbl false true;
            print_mum j_seedNO mum_tbl false true;
            end;
            (*we don't need to reload mi and mj here, for we only modify
            * positions of another match in different sequence of mi *)
            let mi = get_position_by_seqNO i_pos_lst seqNO in
            let i_seqNO, i_left_end,i_right_end = 
                mi.sequence_NO, mi.left_end,mi.right_end in
            let j_idx = find_index (Array.of_list j_pos_lst) i_seqNO 
            (fun looking_item item ->
                if(looking_item = item.sequence_NO)&&
                ( (i_left_end-item.left_end)= d )
                then 0 else (-1) ) in
            if (j_idx<0) then 
                Printf.printf "could not find j_pos for i:(%d,%d,%d)\n%!"
                i_seqNO i_left_end d;
            assert(j_idx>=0);
            let mj = List.nth j_pos_lst j_idx in
            let expand_to_left =  
                if (mj.left_end<i_left_end) then true else false
            in
            let expand_to_right =
                if (mj.right_end>i_right_end) then true else false
            in
            let new_left_end = 
                if expand_to_left  then mj.left_end else i_left_end 
            in
            let new_right_end = 
                if expand_to_right then mj.right_end else i_right_end 
            in
            if debug then 
                Printf.printf "new_l_end/r_end=%d/%d\n%!" new_left_end new_right_end;
            let seed_weight = new_right_end - new_left_end+1 in
            if expand_to_left then begin
                (*do we need to modify seedNO2seq_tbl, seq2seedNO_tbl?*)
                remove_from_pos2seed_tbl pos2seed_tbl_left i_seqNO
                i_left_end i_seedNO mum_tbl;
                add_to_pos2seed_tbl pos2seed_tbl_left i_seqNO mj.left_end
                i_seedNO seed_weight mi.orientation mum_tbl;
            end ;
            if expand_to_right then begin
                remove_from_pos2seed_tbl pos2seed_tbl_right i_seqNO
                i_right_end i_seedNO mum_tbl;
                add_to_pos2seed_tbl pos2seed_tbl_right i_seqNO mj.right_end
                i_seedNO seed_weight mi.orientation mum_tbl;
            end;
            let new_pos = {mi with left_end = new_left_end; right_end =
            new_right_end }in
            update_positions_by_seqNO i_seedNO i_seqNO new_pos mum_tbl ;  
        ) nei_list;
    end;
    (*set j's subsuming_pointer to i*)
    let j_mum = try (Hashtbl.find mum_tbl j_seedNO) with 
    | Not_found -> failwith ("in make a chain, not found err 3") in
    Hashtbl.remove mum_tbl j_seedNO;
    Hashtbl.add mum_tbl j_seedNO {j_mum with subsuming_pointer = i_seedNO };
    if debug then begin
        Printf.printf "after chain: \n%!";
        print_mum i_seedNO mum_tbl false true;
        print_mum j_seedNO mum_tbl false true;
    end
(* chainable_lst *)



let extend_seeds procrastination_queue mum_tbl pos2seed_tbl_left pos2seed_tbl_right =
    let debug = false in
    if debug then Printf.printf "\n ===========   extend seeds ============= \n%!";
    List.iter (fun (i_seedNO, i_size) ->
        if debug then Printf.printf "+++++++ extend seeds on seed#%d +++++++\n%!" 
        i_seedNO;
        let i_mum = try (Hashtbl.find mum_tbl i_seedNO) with 
        | Not_found -> failwith ("in extend seeds, not found err 1") in
        (*this should be a pop from neighborhood_list*)
        if ( i_mum.subsuming_pointer <> (-1)) then begin 
            if debug then 
            Printf.printf " is subsumed by seed# %d\n%!" i_mum.subsuming_pointer; 
        end
        else if ( i_mum.extendable != 0 ) then begin
            if debug then Printf.printf "this one is not extendable from\n%!";
        end
        else begin
            let find_a_chain = ref 1 in
            while ( !find_a_chain = 1 ) do
                if debug then begin
                Printf.printf "=== begin of while ===> check pos2seed tbl \n%!";
                (*print_position2seedtbl pos2seed_tbl_left;
                print_position2seedtbl pos2seed_tbl_right; *)
                end; 
                update_neighborhood i_seedNO mum_tbl pos2seed_tbl_left
                pos2seed_tbl_right false;
                (* build_neighborhood procrastination_queue mum_tbl
                pos2seed_tbl_left pos2seed_tbl_right; *)
                (*update_subset_list mum_tbl;*)
                let i_mum = try (Hashtbl.find mum_tbl i_seedNO) with 
                | Not_found -> failwith ("in extend seeds, not found err 2") in
                if debug then begin
                    Printf.printf "extend seed# %d\n%!" i_seedNO;
                    print_mum i_seedNO mum_tbl true true;
                end;
                    if debug then Printf.printf "categorize neighborhood :\n%!";
                    let nei_list = i_mum.neighborhood_lst in
                    let overlap_list = (*nei list for each j_seedNO*)
                        categorize_neighboorhood_lst nei_list i_seedNO i_size mum_tbl in
                    find_a_chain := 0;
                    if (List.length overlap_list)>0 then begin
                        find_a_chain := 1;
                        (*we only care about the closest overlap-lst to i_mum*)
                        if debug then 
                            List.iter ( fun (dis,nei_lst) -> 
                            Printf.printf "closest nei list with distance = %d\n%!" dis;
                            print_neighborhood_list nei_lst;
                            ) overlap_list;
                        make_a_chain i_seedNO overlap_list mum_tbl pos2seed_tbl_left
                        pos2seed_tbl_right ;
                    end
                    else ();
                    (*
                    let c_l,c_r = List.fold_right (fun item (accl,accr)->
                        let (_,_,_,distance) = item in
                        if (distance>0) then item::accl,accr
                        else accl,item::accr
                    ) chainable_list ([],[]) in
                    if debug then print_mum i_seedNO mum_tbl true true;
                    if debug then Printf.printf "chain left : \n%!";
                    let chainable_lst = sort_neighborhood_lst_by_dis c_l in
                    find_a_chain := 0;
                    if (List.length chainable_lst)>0 then begin
                        find_a_chain := 1; 
                        make_a_chain i_seedNO chainable_lst mum_tbl pos2seed_tbl_left pos2seed_tbl_right 1;
                    end;
                   if debug then  Printf.printf "chain right : \n%!";
                    let chainable_lst =  sort_neighborhood_lst_by_dis c_r in
                    if (List.length chainable_lst)>0 then begin
                        find_a_chain := 1; 
                        make_a_chain i_seedNO chainable_lst mum_tbl pos2seed_tbl_left pos2seed_tbl_right (-1);
                    end; *) 
               (*if debug then begin
                   Printf.printf "check mum_tbl before end of while\n%!";
                   print_mum i_seedNO mum_tbl true true;
                   Printf.printf " === end of while ===\n%!";
               end;  *)
            done; (* end of while ( !find_a_chain = 1 ) *)
        end;
    ) procrastination_queue

    
let print_result_position seq poslst seedNO =
    List.iter (fun pos ->
        let l = pos.left_end 
        and r = pos.right_end 
        and ori=pos.orientation in
        let sign = if (ori=1) then '+' else '-' in
        for i=0 to (l-1) do
            Printf.printf "   %!"
        done;
        for i=0 to (r-l) do
            Printf.printf "%c%2i%!" sign seedNO
        done;
        Printf.printf "\n%!";
    ) poslst


let is_a_number chr =  
        if (((chr>='0')&&(chr<='9'))||(chr='-')) then 
            true 
        else 
            false
    
let charlst_to_int str start_idx =
    let res = ref 0 in
    let code0 = Char.code '0' in
    let nextchr = ref 'a' in
    let idx = ref start_idx in
    nextchr := str.[!idx] ;
    let ori = 
        if (!nextchr='-') then (-1) else 1 in
    if (!nextchr='-') then begin
               idx := !idx + 1;
              nextchr := str.[!idx] 
    end;
    while (is_a_number !nextchr) do
           let num = (Char.code !nextchr) - code0 in
           res := !res * 10 + num ;
           idx := !idx + 1;
           nextchr := str.[!idx];
    done;
    !res,!idx,ori

let get_mauve_result filename = 
    let lines = ref [] in
    let chan = open_in filename in
    let _ = 
        try 
        while true do
        lines := input_line chan :: !lines
        done;
        with End_of_file -> close_in chan
    in
    let lines = List.rev !lines in
    (* List.iter (fun str -> Printf.printf "%s\n%!" str) lines; *)
    let mauve_result = List.fold_right (fun str acc ->
        if (is_a_number str.[0]) then begin
            let reslst = ref [] in 
            let blank = ref 0 in
            let idx = ref 0 in 
            while (!idx <(String.length str)&&(!blank<3)) do
                if (!blank=0)&&(is_a_number str.[!idx]) then begin
                    let size,newidx,_ = charlst_to_int str !idx in
                    idx := newidx ;
                    blank :=1;
                    reslst := !reslst@[size];
                end
                else if (!blank=1)&&(is_a_number str.[!idx]) then begin
                    let leftpos,newidx,ori = charlst_to_int str !idx in
                    idx := newidx ;
                    reslst := !reslst@[ori*(leftpos)];
                    blank := 2;
                end
                else if (!blank=2)&&(is_a_number str.[!idx]) then begin
                    let rightpos,newidx,ori = charlst_to_int str !idx in
                    idx := newidx ;
                    reslst := !reslst@[ori*(rightpos)];
                    blank := 3;
                end
                else idx := !idx +1
            done;
            !reslst::acc
        end
        else
            acc
    ) lines [[]] in
    mauve_result

let print_sub_seq seqlstlst seqNO leftend size =
    let leftend = (abs leftend)-1 in
    assert((List.length seqlstlst)>seqNO);
    let seq = List.nth seqlstlst seqNO in
    assert( (List.length seq)>=(leftend+size) );
    let subseq = Array.to_list (Array.sub (Array.of_list seq) leftend size) in
    print_int_list subseq;
    Printf.printf "\n%!"

let get_sub_seq seq leftend size = 
    assert( (List.length seq)>=(leftend+size) );
     Array.to_list (Array.sub (Array.of_list seq) leftend size)

let get_diff in_lstlst1 in_lstlst2 seqlst =
    let in_lst1 = List.map (fun item ->
        assert((List.length item)=3);
        let size = List.nth item 0
        and lefta = List.nth item 1
        and leftb = List.nth item 2 in
        (size,lefta,leftb)
    ) in_lstlst1 in
    let in_lst2 = List.map (fun item ->
        assert((List.length item)=3);
        let size = List.nth item 0
        and lefta = List.nth item 1
        and leftb = List.nth item 2 in
        (size,lefta,leftb)
    ) in_lstlst2 in
    List.iter (fun item ->
        let (size1,lefta,leftb) = item in
            if ( List.mem item in_lst2)
            || ( List.mem (size1,-lefta,-leftb) in_lst2) then ()
            else begin
                Printf.printf "[%d,%d,%d] is not in our record\n%!" size1 lefta leftb;
                print_sub_seq seqlst 0 lefta size1;
                print_sub_seq seqlst 1 leftb size1;
            end
    ) in_lst1;
    List.iter (fun item ->
        let (size2,lefta,leftb) = item in
            if ( List.mem item in_lst1)
            || ( List.mem (size2,-lefta,-leftb) in_lst1) then ()
            else begin
                Printf.printf "[%d,%d,%d] is not in mauve record\n%!" size2 lefta leftb;
                print_sub_seq seqlst 0 lefta size2;
                print_sub_seq seqlst 1 leftb size2;
            end
    ) in_lst2;
    Printf.printf "end of comparing our mums with mauve's result\n%!"


let build_local_mums mum_tbl seedseq_tbl seedNO_tbl pos2seed_tbl_left
pos2seed_tbl_right seed2pos_tbl debug debug_neighborhood =
    if debug then Printf.printf "build local mums start ... \n%! ";
    Hashtbl.iter (fun (seqNO,index) record_list  ->
        List.iter (fun (seedNO,seed_weight,dir) ->
(* we did mark unextendable seeds in "remove_useless_seed" with pair of (-1) *)
            let ext_sign = 
                let (sign1,sign2) = List.hd ( Hashtbl.find seed2pos_tbl seedNO ) in
                if (sign1!=(-1)) then 0
                else sign2
            in
            if (ext_sign = 0)||(ext_sign=1)||(ext_sign=3) then begin
            let left_end,right_end = index, index + seed_weight - 1 
            and orientation = dir in
            let new_pos = { sequence_NO = seqNO;left_end = left_end; right_end = right_end;
            orientation = orientation; } in
            if (Hashtbl.mem mum_tbl seedNO) then
                begin
                    let old_mum = Hashtbl.find mum_tbl seedNO in
                    let old_pos_lst = old_mum.positions in
                    let new_pos_lst = sort_positions_by_seqNO old_pos_lst@[new_pos] in
                    let old_size = old_mum.size in
                    let new_mum = { old_mum with 
                    positions = new_pos_lst; size = old_size +1 } in
                    Hashtbl.remove mum_tbl seedNO;
                    Hashtbl.add mum_tbl seedNO new_mum;
                end
            else 
                Hashtbl.add mum_tbl seedNO 
                { seedNO = seedNO ; positions = [new_pos]; size = 1; 
                  neighborhood_lst = [];
                  subsuming_pointer = (-1);
                  extendable = ext_sign;
                  mumscore = 0;
                }
            end
            else () (*ext_sign=2, no need to add this one to mum_tbl*)
        )record_list
    ) pos2seed_tbl_left;
    (*extension start*)
    let procrastination_queue = sort_tbl_by_size mum_tbl in
    (*if debug then begin 
    Printf.printf "check procrastination queue: \n%!";
    List.iter (fun (seedNO,size) ->
        Printf.printf "seed.# %d (size = %d)\n%!" seedNO size;
    ) procrastination_queue;
    end; *)
    build_neighborhood procrastination_queue mum_tbl pos2seed_tbl_left
    pos2seed_tbl_right debug_neighborhood;
    if debug then begin
    Printf.printf "check mum tbl after build neighborhood, size=%d (un-extendable included)\n%!"
    (Hashtbl.length mum_tbl);
    (* Hashtbl.iter (fun seedNO mum -> print_mum seedNO mum_tbl false false) mum_tbl;
    *)
    end; 
    (*extend seeds*)
    extend_seeds procrastination_queue mum_tbl pos2seed_tbl_left pos2seed_tbl_right;
    if debug then begin
        Printf.printf "check mum tbl after extend seeds, size=%d (un-extendable included)\n%!"
        (Hashtbl.length mum_tbl);
       (* Hashtbl.iter (fun seedNO mum -> 
            if (mum.subsuming_pointer=(-1)) then 
                print_mum seedNO mum_tbl false false
        ) mum_tbl; *)
    end   


let to_ori_position position lcb_range_lstlst =
    let debug = false in
    let leftend = position.left_end and rightend = position.right_end in
    let len = rightend - leftend +1 in
    let seqNO = position.sequence_NO in
    if debug then
        Printf.printf "to_ori_position, seqNO=%d,l/r = %d/%d, %!" seqNO leftend rightend;
    let lcb_range_lst = List.nth lcb_range_lstlst seqNO in
    let (last_acc_size,last_right,ori_left,ori_right) = 
        List.fold_left (fun (acc_size,pre_r,res_left,res_right) (l,r) ->
            if debug then 
                Printf.printf "acc_size=%d,pre_r=%d,res_left=%d,res_right=%d,\
               range:[%d,%d]\n%!" acc_size pre_r res_left res_right l r; 
            if (res_left <> (-1)) then 
                (0,0,res_left,res_right)
            else begin
                let current_acc_size = acc_size+(l-1-(pre_r+1)+1) in
                if debug then Printf.printf "current_acc_size=%d\n%!" current_acc_size;
                (*what if (leftend,rightend) is split by a LCB?
                * -- we discard this new mum,*)
                if (leftend>=acc_size)&&(rightend<current_acc_size) then
                    let res_l = leftend-acc_size+1+pre_r in
                    let res_r = res_l + len -1 in
                    (0,0,res_l,res_r)
                else if (rightend>=current_acc_size) then (*not there yet*)
                    (current_acc_size, r, (-1),(-1))
                else (*split by existing lcb*)
                    (current_acc_size, r, (-2),(-2))
            end
        ) (0,(-1),(-1),(-1)) lcb_range_lst 
    in
    if debug then 
        Printf.printf "ori left/right = %d/%d\n%!" ori_left ori_right;
    let ori_left =
        if ori_left=(-1) then 
            last_right+leftend-last_acc_size+1
        else 
            ori_left in
    let ori_right =
        if ori_right=(-1) then
            ori_left + rightend-leftend
        else 
            ori_right in
    if debug then Printf.printf "ori left/right = %d/%d \n%!" ori_left ori_right;
    (ori_left,ori_right)

(*transpose *)
let transpose_back lcb_range_lstlst new_mum_tbl new_seedNO2seq_tbl
mum_tbl pos2seed_tbl_left pos2seed_tbl_right seed2pos_tbl
= 
    let debug = false in
    if debug then begin
        Printf.printf "transpose back, range list is :\n%!";
        List.iter (fun rangelst ->
            List.iter (fun (left,right) ->
                Printf.printf "(%d,%d),%!" left right
            ) rangelst;
            Printf.printf "\n%!";
        ) lcb_range_lstlst;
    end;
    Hashtbl.iter (fun new_seedNO new_mum ->
            if debug then Printf.printf "mum before transpose: %!";
            if debug then print_mum new_seedNO new_mum_tbl false true;
            let sign = ref true in (*true means we can add the new mum to mum_tbl*)
            let ori_positions = List.map (fun mi ->
                let ori_l,ori_r = to_ori_position mi lcb_range_lstlst in
                if (ori_l=(-2)) then sign:=false;
                {mi with left_end = ori_l; right_end = ori_r}
            ) new_mum.positions in
            if (!sign = true) then begin
                
                let seedNO = get_a_seedNO () in
                assert( (Hashtbl.mem mum_tbl seedNO)=false);
                Hashtbl.add mum_tbl seedNO 
                { new_mum with 
                        seedNO = seedNO; 
                        positions = ori_positions; 
                        neighborhood_lst = [] };
                let debug2 = false in
                if debug2
                then begin
                    Printf.printf "add new seed to mum tbl:%!";
                    print_mum seedNO mum_tbl false true
                end;
                List.iter ( fun mi ->
                    let seqNO = mi.sequence_NO and ori = mi.orientation in
                    let left = mi.left_end and right = mi.right_end in
                    let weight = right-left+1 in
                    add_to_pos2seed_tbl pos2seed_tbl_left seqNO left seedNO
                    weight ori mum_tbl;
                    add_to_pos2seed_tbl pos2seed_tbl_right seqNO right seedNO
                    weight ori mum_tbl;
                    add_to_seed2pos_tbl seed2pos_tbl seqNO left seedNO weight ori;
                )ori_positions;
                if debug then Printf.printf "mum after transpose :%!";
                if debug then print_mum seedNO mum_tbl false true;
        end;
        return_a_seedNO new_seedNO;
    ) new_mum_tbl;
    Hashtbl.iter (fun seedNO mum ->
        let debug_nei = false in
        update_neighborhood seedNO mum_tbl pos2seed_tbl_left pos2seed_tbl_right
        debug_nei;
    ) mum_tbl;
    if debug then Printf.printf "end of transposition back\n%!"
    


(* we only allow at most one mum starts each position, 
* when a position has no qualified mum, we return (-1,-1,-1) 
* (mum that subsumed by other mum, or shows up more than once
* in a sequence - which is unextendable, is ignored)*)
let get_extendable_record_pos2seed_tbl pos2seed_tbl pos mum_tbl =
    let recordlst = Hashtbl.find pos2seed_tbl pos in
    let recordlst = List.filter (fun record ->
        let seedNO,weight,ori = record in
        let mum = try (Hashtbl.find mum_tbl seedNO)
        with | Not_found -> failwith "not found in get lcb 1" in
        if (mum.subsuming_pointer=(-1))&&(mum.extendable=0) then true
        else false
    ) recordlst in
    assert( (List.length recordlst)<=1 );
    if (List.length recordlst)=0 then ((-1),(-1),(-1))
    else   List.hd recordlst
    (* no longer need to break the tie between extendable mums, since we only keep
    * one of them extendable now.
    else begin
        let _, highest_multi_lst = 
            List.fold_right (fun recordseed (max_multi, acc) ->
                let seedNO,weight,ori = recordseed in
                let multiplicity = List.length mum.positions in
                if (multiplicity>max_multi) then
                    (multiplicity,[recordseed])
                else if (multiplicity=max_multi) then
                    (max_multi, recordseed::acc)
                else 
                    (max_multi,acc)
            ) recordlst (0,[])
        in
        if (List.length highest_multi_lst)=0 then ((-1),(-1),(-1))
        else if (List.length highest_multi_lst)=1 then
            List.hd highest_multi_lst
        else
            List.hd (
            List.sort (fun x y ->
                let (_,weightx,_) = x and (_,weighty,_) = y in
                compare (-weightx) (-weighty)
            ) highest_multi_lst )
    end *)


    

let get_mum_lst_for_each_seq mum_tbl pos2seed_tbl seqlst_size seq_size_lst 
= 
    let debug = false in
    if debug then Printf.printf "get lcb lst for each seq\n%!";
    let seedNOlstlst = ref [[]] in
    for i = 0 to (seqlst_size-1) do
        if debug then Printf.printf "work on sequence %d\n%!" i;
        let reslst = ref [] in
        for j=0 to (List.nth seq_size_lst i)-1 do
            if (Hashtbl.mem pos2seed_tbl (i,j)) then begin
                let (seedNO,seed_weight,ori) = 
                    get_extendable_record_pos2seed_tbl pos2seed_tbl (i,j) mum_tbl in
                if (seedNO<>(-1)) then
                    begin
                    if debug then
                        Printf.printf "seed#%d at pos (%d,%d)\n%!" seedNO i j;
                    reslst := !reslst@[ (ori*seedNO) ]
                    end
                else ();
            end
            else ();
        done;
        seedNOlstlst := !reslst::!seedNOlstlst
    done;
    let seedNOlstlst = List.filter (fun item -> (List.length item)>0) 
    (List.rev !seedNOlstlst) in
    if debug then begin
    Printf.printf "seedNOlst for each seq is :\n%!"; print_int_lstlst seedNOlstlst;
    end;
    seedNOlstlst


(*this works for two sequences input only*)
let get_break_point_matrix inarrarr debug =
    let sizex = Array.length inarrarr 
    and sizey = Array.length (inarrarr.(0)) in
    assert(sizex>0); 
    assert(sizey>0);
    if debug then 
        Printf.printf "get_break_point_matrix, sizex=%d,sizey=%d\n%!" sizex sizey;
    let bkmatrix = Array.make_matrix sizex sizey (-1) in 
    for i=0 to (sizey-1) do
        bkmatrix.(0).(i) <- (i+1);
    done;
    for i=1 to (sizex-1) do
        for j=0 to (sizey-1) do
            let item = inarrarr.(0).(j) in
            let idx,ori = 
                let tmp = find_index inarrarr.(i) item compare in
                if (tmp<0) then
                    find_index inarrarr.(i) (-item) compare ,-1
                else tmp,1
            in
            if (idx<0) then begin
                Printf.printf "i=%d,j=%d,item=%d\n%!" i j item;
            end;
            assert(idx>=0);
            bkmatrix.(i).(idx)<-(ori*(j+1));
        done;
    done;
    if debug then  print_int_matrix bkmatrix;
    let bk_number = ref 0 in
    for i=1 to (sizex-1) do
       for j=(sizey-1) downto 1 do
            let b = bkmatrix.(i).(j) and a = bkmatrix.(i).(j-1) in
            if ((b*a>0)) && ( (b-a)=1 ) then begin
                bkmatrix.(i).(j) <- 0;
                let idx = bigger_int (abs a) (abs b) in
                bkmatrix.(0).(idx-1) <- 0;
            end
            else begin
                bk_number := !bk_number + 1;
                bkmatrix.(i).(j) <- 1;
            end
        done;
        bkmatrix.(i).(0) <- 1;
        bk_number := !bk_number + 1;
    done;
    for j=0 to (sizey-1) do
        if bkmatrix.(0).(j)<>0 then bkmatrix.(0).(j)<-1
        else ();
    done;
    bkmatrix.(0).(0) <- 1;
    bkmatrix, !bk_number

(*get range of current lcb in each sequence, seedNOlst is the list of seedNO of that lcb*)
let get_range_of_a_lcb seedNOlst seqNO mum_tbl = 
    let min_abs_code = List.hd (List.sort (fun x y -> compare (abs x) (abs y) )
    seedNOlst) in
    let ori =  if (min_abs_code>0) then 1 else (-1) in 
    let weight = List.length seedNOlst in
    let seedNO_first = abs (List.hd seedNOlst) in
    let seedNO_last = abs (List.nth seedNOlst (weight-1)) in
    let mumfirst = Hashtbl.find mum_tbl seedNO_first 
    and mumlast = Hashtbl.find mum_tbl seedNO_last in
    let posfirst = get_position_by_seqNO mumfirst.positions seqNO in
    let poslast = get_position_by_seqNO mumlast.positions seqNO in
    let leftend = 
        if (posfirst.left_end < poslast.left_end) then posfirst.left_end
        else poslast.left_end 
    in
    let rightend = 
        if (posfirst.right_end < poslast.right_end) then poslast.right_end
        else posfirst.right_end 
    in
    let debug = false in
    if debug then begin
    Printf.printf "get range of lcb: \n%!";
    List.iter (fun x -> print_mum x mum_tbl false true) seedNOlst;
    Printf.printf "=> %d,%d\n %!" leftend rightend;
    end;
    (leftend,rightend),ori

(*convert (a,t,g,c) to (0,1,2,3), then we can lookup hodx table for their cost.*)
let atgc_2_idx x = 
    if x=1 then 0 
    else if x=2 then 1
    else if x=4 then 2
    else if x=8 then 3
    else failwith "unkown gene code"

(*score bwteen two sequences based on hodx matrix *)
let get_score_from_2seq sequence1 sequence2 ori = 
    let debug = false in
    let seq1,seq2 = 
        if ori=1 then sequence1,sequence2
        else 
            List.rev (rev_comp sequence1), sequence2
    in      
    if debug then begin
    print_int_list seq1; print_int_list seq2;
    end;
    List.fold_left2 (fun acc x y -> 
        acc + hodx_matrix.(atgc_2_idx x).(atgc_2_idx y) 
        ) 0 seq1 seq2
        
(* get_score returns score of mum, and yes this works for more than 2
* input sequence, but we only care about lcbs between 2 sequences *)
let get_mum_score seedNO mum_tbl old_seqlst =
    let mum = try (Hashtbl.find mum_tbl seedNO) with
    | Not_found -> failwith "not found get_mum_score\n%!" in 
    let poslst = mum.positions in
    let first_range = List.hd poslst in
    let first_seq = 
        get_sub_seq (List.nth old_seqlst first_range.sequence_NO) 
        first_range.left_end (first_range.right_end-first_range.left_end+1) in
    let first_ori = first_range.orientation in
    let score,_,_ = 
        List.fold_left (fun (sum_score,previous_seq,previous_ori) lcbrange ->
        let seqNO,lend,rend,current_ori =
            lcbrange.sequence_NO,lcbrange.left_end,
            lcbrange.right_end,lcbrange.orientation
        in
        let current_seq = 
            get_sub_seq (List.nth old_seqlst seqNO) lend (rend-lend+1) in
        let ori = if (current_ori=previous_ori) then 1 else (-1) in
        let score = get_score_from_2seq previous_seq current_seq ori in
        sum_score + score, current_seq, current_ori
    ) (0,first_seq,first_ori) (List.tl poslst)
    in
    score

let get_sorted_lcb_by_ratio lcb_tbl =
    let lst = ref [] in
    Hashtbl.iter (fun key record ->
        lst := record :: (!lst);
        lst := List.sort (fun x y -> compare x.ratio y.ratio) (!lst)
    ) lcb_tbl;
    !lst

let get_num_of_low_ratio_lcb lcb_tbl min_ratio check_weight =
    let res = ref 0 in
    Hashtbl.iter (fun key record ->
        let check1 =
            if (check_weight) then
                if ((List.length key)>=minimum_lcb_weight) then true
                else false
            else true
        in
        if (check1)&&(record.ratio < min_ratio) 
        then begin
            res:= !res +1;
        end;
    ) lcb_tbl;
    !res

(*this is a test we do, maybe we can replace the old way we score lcb*)
let get_new_score lcb_tbl seqsize_lst =
    let seqlst_len = 2 in (*we only deal with two sequences now*)
    let sum_score = Array.make seqlst_len 0 and sum_len = Array.make seqlst_len 0 in
    Hashtbl.iter (fun key record ->
        let range_lst = record.range_lst in
        List.iter (fun range->
            let idx = range.sequence_NO in
            let length = range.right_end - range.left_end in
            sum_len.(idx) <- sum_len.(idx)+length;
            sum_score.(idx) <- sum_score.(idx)+record.score;
        )range_lst;
    ) lcb_tbl;
    let sum_len = Array.to_list sum_len 
    and sum_score = Array.to_list sum_score in
    let res_R = List.map2 (fun x y -> (float x)/.(float y)) sum_score sum_len in
    let avg_R = get_avg_of_floatlst res_R in
    let coverage = List.map2 (fun len size ->
        (float len) /. (float size) ) sum_len seqsize_lst in
    let avg_cov = get_avg_of_floatlst coverage in
    let bk_number = Hashtbl.length lcb_tbl in
    Printf.printf "avg_R = %f,%!" avg_R;
    Printf.printf "avg_cov = %f, R*cov=%f, bk_number = %d\n%!" 
    avg_cov (avg_cov*.avg_R) bk_number



(*this only works for two sequences ,for now*)
(*  List.nth lcbs n is the lcb list, each lcb is a list of mums. 
* so lcb list looks like this; [muma1,muma2,...],[mumb1,mumb2,...],,..
* Instead of getting rid of all low-weight lcbs, we only get rid one of them.
* "all of the LCBs get scored, then the sum-of-pairs
breakpoint score gets calculated.  Then each one of the LCBs is
evaluated for removal by calculating how much larger the sum-of-pairs
LCB score would be if that block were removed.  The LCB that would
create the largest increase in the sum-of-pairs breakpoint score gets
removed."  
" The total score would be calculated as the sum of the individual MUM scores, minus
(number of breakpoints)*(breakpoint penalty)". *)
(*lcb_tbl is empty before this function*)
let update_lcbs_range_and_score lcbs mum_tbl lcb_tbl input_seqlst = 
    let debug = false in
   (*fill in range of each lcb with mum_tbl*)
    let idx = ref (-1) in
    List.iter (fun lcblist ->
        idx := !idx + 1;
         List.iter (fun record ->
            let (lend,rend),ori = get_range_of_a_lcb record !idx mum_tbl in
            let abs_record = get_abs_lst record in
            if (Hashtbl.mem lcb_tbl abs_record) then begin
                let old_lcb = Hashtbl.find lcb_tbl abs_record in
                let newlcbrange = 
                { sequence_NO = (!idx); left_end=lend; right_end=rend; orientation=ori;} 
                in
                let old_range_lst = old_lcb.range_lst in
                Hashtbl.remove lcb_tbl abs_record;
                Hashtbl.add lcb_tbl abs_record 
                {old_lcb with range_lst = newlcbrange::old_range_lst};
            end
            else begin
                assert( !idx = 0 );
                let firstitem = 
                    [{sequence_NO= (!idx);left_end=lend;right_end=rend;
                    orientation=ori}]
                in
                Hashtbl.add lcb_tbl abs_record 
                {seedNOlst = abs_record; range_lst = firstitem; 
                score=0; ratio = 0.0; ref_code = (-1) }
            end;
        ) lcblist
    )lcbs;
    (* Printf.printf "check lcb range lstlst:\n%!";  print_lcbs_range lcb_range_lst_lst; *)
   (*fill in score for each lcb*)
    Hashtbl.iter (fun seedNOlst oldlcb ->
        let total_score = 
            List.fold_left (fun acc seedNO ->
                 let mumi = try (Hashtbl.find mum_tbl seedNO) with
                 | Not_found -> failwith "not found 1 update_lcbs_range_and_score\n%!" 
                 in
                 acc + mumi.mumscore
            ) 0 seedNOlst
        in
        let largest_range = get_largest_range_of_mi_lst oldlcb.range_lst in 
        let ratio = (float total_score)/.(float largest_range) in
        if debug then begin
            Printf.printf "update lcb_tbl on %!" ; print_lcb oldlcb;
            Printf.printf "score <-- %d,ratio <-- %f\n%!" total_score ratio;
        end; 
        Hashtbl.remove lcb_tbl seedNOlst;
        Hashtbl.add lcb_tbl seedNOlst 
        {oldlcb with score = total_score; ratio = ratio };
    ) lcb_tbl


let get_lcbs_range lcb_tbl old_seqlst special_case = 
    let lcb_range_arr = Array.make (List.length old_seqlst) [] in
    let min_lcb_w = 
        if special_case then 0 else minimum_lcb_weight 
    in
    Hashtbl.iter (fun seedNOlst record ->
        if (List.length seedNOlst)>= min_lcb_w then
            List.iter (fun mi ->
                let oldlst = lcb_range_arr.(mi.sequence_NO) in
                let newlst = (mi.left_end,mi.right_end)::oldlst in
                lcb_range_arr.(mi.sequence_NO) <- newlst ;
            ) record.range_lst;
    ) lcb_tbl;
    let lcb_range_lst_lst = Array.to_list lcb_range_arr in
    let lcb_range_lst_lst = (*sort the range by leftend*)
        List.map (fun x -> 
            List.sort (fun (aleft,_) (bleft,_) -> compare aleft bleft) x
        ) lcb_range_lst_lst
    in
    lcb_range_lst_lst

let get_lcb_key_by_range seqNO range lcb_tbl =
    let (leftend,rightend) = range in
    let reskey= ref [] and rescode = ref 0 in (*-1 is not a good idea for the
    rescode could be (code=1)*(ori=-1)=-1 *)
    Hashtbl.iter (fun key record ->
        List.iter (fun mi ->
            if (leftend=mi.left_end)&&(rightend=mi.right_end)&&
            (mi.sequence_NO=seqNO) then begin
                reskey := key;
                rescode := record.ref_code*mi.orientation;
            end
        ) record.range_lst;
    ) lcb_tbl;
    !reskey,!rescode

(*given a lcbs range list list, get_full_range_lstlst returns a range list list
* of both lcbs and area outside lcbs 
* what happens if lcbs range start from 0 of the sequence ?
* (0,-1) will become the first range of full_range_lstlst, does not make 
* much sense, get rid of this range. also take care of the case when one lcb
* ends at the end of a sequence. *)
let get_full_range_lstlst lcb_range_lstlst in_seqsize_lst =
    List.map2 (fun lcb_range_list total_size ->
        let last_rightend,acc=
            List.fold_left (fun (pre_rightend,acc) (leftend,rightend) ->
                assert(pre_rightend<=leftend);
                let newacc = 
                    if (leftend>0) then acc@[((pre_rightend+1),(leftend-1))] 
                    else acc
                in
                let newacc = newacc@[(leftend,rightend)] in
                rightend,newacc 
                ) (-1,[]) lcb_range_list
        in
        let res = 
            if (last_rightend+1=total_size) then acc
            else acc@[(last_rightend+1,total_size-1)]
        in
        res
    ) lcb_range_lstlst in_seqsize_lst


    
let get_seq_outside_lcbs old_seqlst lcb_tbl =
    let debug = false in
    if debug then Printf.printf "get seq outside lcbs :\n%!";
    let lcb_range_lst_lst = get_lcbs_range lcb_tbl old_seqlst false in
    if debug then print_lcbs_range lcb_range_lst_lst;
    let seq_outside_lcb_lst = List.map2 (fun seq lcb_range ->
        let seqarr = Array.of_list seq in
        let seqlen = List.length seq in
        let last_rightend,tmpseq = 
            List.fold_left (
                fun (pre_rightend,subseq) (leftend,rightend) ->
                assert(pre_rightend<=leftend);
                let size = leftend-pre_rightend-1 in
                let subarr =
                    if size>0 then
                    Array.sub seqarr (pre_rightend+1) (leftend-pre_rightend-1) 
                    else [||]
                in
                rightend,subseq@(Array.to_list subarr)
            ) (-1,[]) lcb_range 
        in
        let last_size = seqlen-1-last_rightend in
        let subend = 
            if last_size>0 then
            Array.sub seqarr (last_rightend+1) (seqlen-1-last_rightend) 
            else [||]
        in
        tmpseq@(Array.to_list subend)
    ) old_seqlst lcb_range_lst_lst in
    lcb_range_lst_lst, seq_outside_lcb_lst
    


(* local collinear blocks, works only for two sequences now *)
let build_LCBs seedNOlstlst  =
    let debug = false in
    if debug then Printf.printf "build LCBS \n %!";
    let resmatrix = Array.of_list (List.map (fun lst -> Array.of_list lst)
    seedNOlstlst) in 
    let debug_bk_matrix = false in
    let bkmatrix,_ = get_break_point_matrix resmatrix debug_bk_matrix in
    if debug then print_int_matrix bkmatrix;
    let lcbs = Array.to_list (Array.mapi (fun i bkarr ->
        let res = ref [[]] in
        let j = ref (-1) in
        let lastacc = List.fold_left(fun acc record ->
            j := !j+1;
            if (record=0) then 
                acc@[resmatrix.(i).(!j)]
            else if (record=1) then begin
                res := (!res)@[acc];
                [resmatrix.(i).(!j)]
            end
            else 
                failwith "what's this? we only record 1 or 0 in bk_matrix";
        ) [] (Array.to_list bkarr) in
        List.filter (fun item -> (List.length item)>0) (!res@[lastacc]);
    ) bkmatrix )
    in
    if debug then begin
        Printf.printf "check lcbs :\n%!";
        List.iter (fun seqlcbs ->
        List.iter (fun lcb -> print_int_list lcb ) seqlcbs;  Printf.printf "\n%!";
        ) lcbs;
    end;
    if debug then Printf.printf "qualified lcbs now: %!";
    let sum_weight =
        List.fold_left (fun sum record ->
        let len = List.length record in
        if (len>=minimum_lcb_weight) then begin 
            if debug then print_int_list record;
            sum + (List.length record)
        end
        else sum
        ) 0 (List.nth lcbs 1)
    in
    if debug then Printf.printf "\n sum_weight=%d, end of building LCBs\n%!" sum_weight;
    lcbs,sum_weight


let get_worst_lcb_from_light_lcbtbl light_lcb_tbl lcbs mum_tbl in_seqlst
bk_penalty total_mum_score init_higher_score init_low_ratio_lcb_num debug =
    let lcb_to_remove = ref [] in
    let lcbs_after_remove = ref [[[]]] in
    let higher_score = ref init_higher_score in
    let lcb_tbl_after_remove = ref (Hashtbl.create init_tb_size) in
    let tmp_tbl = Hashtbl.create init_tb_size in
    Hashtbl.iter (fun seedNOlst lcb ->
        let removed_lcb_score = lcb.score in
        if debug then begin
            Printf.printf "if we remove lcb:%!"; print_lcb lcb;
        end;
        let new_lcbs = List.map (fun lcblst ->
            List.filter (fun record ->( (get_abs_lst record) <> seedNOlst )
             ) lcblst
        ) lcbs in
        let seedNOlstlst = List.map (fun lcblst -> List.flatten lcblst) new_lcbs in
        let resmatrix = Array.of_list (List.map (fun lst -> Array.of_list lst)
        seedNOlstlst) in 
        let debug_bk_matrix = false in
        let _,bk_number = get_break_point_matrix resmatrix debug_bk_matrix in
        let new_score = 
            total_mum_score - removed_lcb_score - bk_number*bk_penalty in
        if debug then 
            Printf.printf " bk_number=%d,score will be %d-%d-%d=%d\n%!" 
            bk_number total_mum_score removed_lcb_score (bk_number*bk_penalty)
            new_score;
    (*just a test -- the old way of scoring sucks, try some new one*)
        let new_seedNOlstlst = 
            List.map (fun lcblst -> List.flatten lcblst) new_lcbs in
        let new_lcbs,new_weight = build_LCBs new_seedNOlstlst in
        Hashtbl.clear tmp_tbl;
        update_lcbs_range_and_score new_lcbs mum_tbl tmp_tbl in_seqlst;
       (* get_new_score tmp_tbl in_seq_size_lst; this is the test*)
        let num_low_ratio_lcb = get_num_of_low_ratio_lcb tmp_tbl
        !minimum_lcb_ratio true in
        if debug then
            Printf.printf "num of low ratio lcb = %d\n%!" num_low_ratio_lcb;
        if (new_score> (!higher_score))&&
        (init_low_ratio_lcb_num >= num_low_ratio_lcb)
        then begin
            higher_score := new_score;
            lcb_to_remove := seedNOlst;
            lcbs_after_remove := new_lcbs;
            lcb_tbl_after_remove := Hashtbl.copy tmp_tbl;
        end;
    ) light_lcb_tbl;
    !higher_score,!lcb_to_remove,!lcbs_after_remove,!lcb_tbl_after_remove
    
     

(*get_worst_lcb returns the lcb key so we can remove it to get a better score.
* also returns the score if we remove the lcb, number of low W lcb left if we
* remove the lcb*)
let get_worst_lcb lcbs lcb_tbl bk_penalty in_seq_size_lst in_seqlst
mum_tbl(*mum_tbl is for the test*) 
=
    let debug = false in
    let total_num = Hashtbl.length lcb_tbl in
    let light_num = get_num_of_low_ratio_lcb lcb_tbl !minimum_lcb_ratio true in
    let total_mum_score = ref 0 in
    Hashtbl.iter (fun seedNOlst record ->
            total_mum_score := !total_mum_score + record.score
    ) lcb_tbl;
    let total_mum_score = !total_mum_score in
    let init_bk_num = List.length (List.hd lcbs) in
    let higher_score = ref (total_mum_score - init_bk_num*bk_penalty)  in
    if debug then begin
            Printf.printf "\n<------get worst lcb, init_bk_num = %d,init_higher_score=%d;%!" 
            init_bk_num !higher_score;
            Printf.printf "total_mum_score = %d,%!" total_mum_score;
            Printf.printf "init low ratio lcb num = %d/%d\n%!" light_num total_num;
    end;
    let lcb_to_remove = ref [] in
    Hashtbl.iter (fun key record -> (*remove high W lcb with low ratio first*)
            if (!lcb_to_remove=[])&&((List.length key)>=minimum_lcb_weight)
            && (record.ratio< !minimum_lcb_ratio) then
                lcb_to_remove := key
    ) lcb_tbl;
    (*we only need light_lcb_tbl when we remove lcb with low W high R, but the size of
    * light_lcb_tbl is useful when we remove lcb with high W and low R*)
    let light_lcb_tbl = Hashtbl.create init_tb_size in
    Hashtbl.iter (fun key record ->
            if ((List.length key)<minimum_lcb_weight) then 
                Hashtbl.add light_lcb_tbl key record
    ) lcb_tbl;
    if (!lcb_to_remove<>[]) then begin (*always remove high w but low R lcb first*)
        let lcb = try(Hashtbl.find lcb_tbl !lcb_to_remove) 
        with |Not_found -> failwith "not found in get worst lcb" in
        let removed_lcb_score = lcb.score in
        if debug then begin
            Printf.printf "remove high weight but low ratio lcb:%!"; 
            print_lcb lcb;
        end;
        let seedNOlst = !lcb_to_remove in
        let new_lcbs = List.map (fun lcblst ->
            List.filter (fun record ->( (get_abs_lst record) <> seedNOlst )
             ) lcblst
        ) lcbs in
        let seedNOlstlst = List.map (fun lcblst -> List.flatten lcblst) new_lcbs in
        let resmatrix = 
            Array.of_list (List.map (fun lst -> Array.of_list lst) seedNOlstlst) in 
        let debug_bk_matrix = false in
        let _,bk_number = get_break_point_matrix resmatrix debug_bk_matrix in
        let new_score = 
            total_mum_score - removed_lcb_score - bk_number*bk_penalty in
        let lcb_tbl_after_remove = Hashtbl.copy lcb_tbl in
        Hashtbl.remove lcb_tbl_after_remove seedNOlst;
        seedNOlst,new_score, (Hashtbl.length light_lcb_tbl),
        new_lcbs,lcb_tbl_after_remove
    end
    else if (Hashtbl.length light_lcb_tbl)>1 then begin
        let init_low_ratio_lcb_num = light_num in
        let higher_score, lcb_to_remove,lcbs_after_remove,lcb_tbl_after_remove = 
            get_worst_lcb_from_light_lcbtbl light_lcb_tbl lcbs mum_tbl in_seqlst
            bk_penalty total_mum_score !higher_score init_low_ratio_lcb_num debug
        in
        if debug then begin
            Printf.printf "++++== RESULT: remove lcb %!"; 
            print_int_list lcb_to_remove;
            Printf.printf " score after remove lcb = %d\n%!" higher_score;
        end;
        let light_lcb_num = 
            if (lcb_to_remove = []) then (Hashtbl.length light_lcb_tbl) 
            (* maybe return 0 to end removing light lcb?*)
            else (Hashtbl.length light_lcb_tbl)-1
        in
        lcb_to_remove,higher_score, light_lcb_num,
        lcbs_after_remove,lcb_tbl_after_remove
    end
    else
        [],!higher_score,1 (*1 to make sure lcb_tbl is update in outer
        function*), lcbs, lcb_tbl



(*break point penaly formula is from Mauve.*)
let get_break_point_penalty input_seq_size_lst num_of_mums =
    let debug = false in
    if debug then 
        Printf.printf "get break point penalty, num_of_mums = %d, %!" num_of_mums;
    let avg_seq_len = get_avg_of_intlst input_seq_size_lst in
    let	avg_seq_len = log( avg_seq_len ) /. log( 2.0 ) in
    let avg_seq_len = 7000. *. avg_seq_len in
    let cons_id = 1. -. (get_cons_id input_seq_size_lst num_of_mums) in
    let bp_estimate = int_of_float (3. *.avg_seq_len*.cons_id*.cons_id*.cons_id*.cons_id) in
    if debug then Printf.printf "bp_estimate = %d \n%!" bp_estimate;
    let res = 
        if (bp_estimate> !min_break_point_penalty) then bp_estimate
        else !min_break_point_penalty
    in
    res

let remove_lcb_from_lcblstlst key_to_remove lcbs =
    List.map (fun lcblst ->
        List.filter (fun record ->
           ( (get_abs_lst record) <> (key_to_remove) )
        ) lcblst
    ) lcbs


(*remove light weight lcbs will remove  
* -- high W low R lcbs
* -- low W lcbs
* (W is weight, R is ratio)
* the result lcb_tbl could still have low W lcbs, if remove them will create
* high W low R lcb, we won't do it*)
let remove_light_weight_lcbs lcbs lcb_tbl mum_tbl old_seqlst in_seq_size_lst
num_of_mums =
    let debug = false in
    let bk_penalty = get_break_point_penalty in_seq_size_lst num_of_mums in 
    if debug then begin
        Printf.printf "start to remove light weight lcbs, bk_penalty = %d\n
        input lcbs is%!" bk_penalty;
        print_lcblst lcbs;
    end;
    let init_lcb_weight = ref 0 in
    (*light_lcb_tbl is for low W and low R lcbs*)
    let light_lcb_tbl = Hashtbl.create init_tb_size in
    Hashtbl.iter (fun key record ->
        if ((List.length key)< minimum_lcb_weight)||(record.ratio< !minimum_lcb_ratio)
        then 
            Hashtbl.add light_lcb_tbl key record
        else begin
            if debug then print_lcb record;
            init_lcb_weight := !init_lcb_weight + (List.length key)
        end
    ) lcb_tbl;
    let total_mum_score = ref 0 in
    Hashtbl.iter (fun seedNOlst record ->
            total_mum_score := !total_mum_score + record.score
    ) lcb_tbl;
    let total_mum_score = !total_mum_score in
    let init_bk_num = List.length (List.hd lcbs) in
    let init_score = total_mum_score - init_bk_num*bk_penalty in
    (*pick the first lcb to remove*)
    (*let seedNOlst_to_remove,_,_ = 
        get_worst_lcb lcbs lcb_tbl bk_penalty in_seq_size_lst old_seqlst mum_tbl
    in*)
    if debug then begin
        Printf.printf "init light weight lcb_tbl size = %d/(total size = %d) \n%!"
        (Hashtbl.length light_lcb_tbl) (Hashtbl.length lcb_tbl); 
       (* print_int_list seedNOlst_to_remove;
        Printf.printf "\n%!";*)
    end;
    (*if (seedNOlst_to_remove<>[]) then begin*) 
        let res_score = ref init_score in
        let res_lcb_tbl = ref lcb_tbl in
        let res_lcbs = ref lcbs in
        (*and key_to_remove = ref seedNOlst_to_remove in*)
        let sign = ref true in
        while (!sign) do
            if debug then 
            Printf.printf "current score = %d\n %!"  !res_score;
            let seedNOlst_to_remove,score,light_lcb_left,lcbs_after_remove,lcb_tbl_after_remove = 
                get_worst_lcb !res_lcbs !res_lcb_tbl bk_penalty in_seq_size_lst
                old_seqlst mum_tbl in
            if debug then begin
                Printf.printf "try to remove lcb : %!";
                print_int_list seedNOlst_to_remove; 
           Printf.printf "new score will be %d, light lcb left=%d\n%!" score light_lcb_left;
            end;
            (*mauve stops when there is only one light weight lcb left*)
            (*when we only have 1 low W lcb left, light_lcb_left will be 1,
            * seedNOlst_to_remove will be [], end of loop*)
            if (!res_score <= score) && (light_lcb_left>1)
            && (seedNOlst_to_remove <> [])
            then begin
                sign := true; 
                res_lcb_tbl := Hashtbl.copy lcb_tbl_after_remove;
                res_lcbs := lcbs_after_remove;
                if debug then begin
                    Printf.printf "current lcbs is :%!";
                    print_lcblst !res_lcbs;
                end;
                res_score := score;
                if debug then 
                    Printf.printf "new score >= old score, continue with loop\n%!";
            end
            else begin
                sign := false;
                if debug then 
                   (* if (light_lcb_left=1) then begin
                        Printf.printf "reach the last one of light weight lcbs,\
                        get out of loop\n%!";
                        res_lcb_tbl := new_lcb_tbl;
                        current_lcbs := new_lcbs;
                        old_score := score;
                    end
                    else *)
                    Printf.printf "no improvement in score,get out of loop\n%!";
            end;
        done; (* while(!sign) *)
        let qualified_lcb_weight = ref 0 in
        Hashtbl.iter (fun key record ->
            let w = List.length key in
            if w>=minimum_lcb_weight then qualified_lcb_weight := !qualified_lcb_weight + w
        ) !res_lcb_tbl;
        if debug then 
            Printf.printf "end of remove_light_weight_lcbs,q_lcb_weight = %d\n%!"
            !qualified_lcb_weight;
        !res_lcbs,!qualified_lcb_weight, !res_lcb_tbl
    (*end
    else begin
        let qualified_lcb_weight = ref 0 in
        Hashtbl.iter (fun key record ->
            if (List.length key)>=minimum_lcb_weight then
                qualified_lcb_weight := !qualified_lcb_weight + 1
        ) lcb_tbl;
        lcbs,!init_lcb_weight,lcb_tbl
    end*)
    


(*for jseedNO, we remove the old record of oldpos of pos2seed_tbl1, add newpos the the
* new position. for pos2seed_tbl2, we only update the weight of the old record
* on pos2 .
* for kseedNO, we create it if it's not already there, and update
* pos2seed_tbl1&2, seed2pos_tbl, also it's own mum record in mum_tbl.*)
let update_tables mum_tbl 
jseedNO seqNO kseedNO pos2seed_tbl1 oldpos newpos pos2seed_tbl2
pos2 newweight seed2pos_tbl new_mi oldori weight_reduce trim_from_left add_to_new_mum = 
    let debug = false in
    if debug then 
    Printf.printf "update tables,seqNO=%d,jseedNO=%d,kseedNO=%d,oldpos=%d,newpos=%d,newweight=%d\
    weight_reduce=%d,pos2=%d trim_from_left:%b,add_to_new_mum:%b\n%!"
    seqNO jseedNO kseedNO oldpos newpos newweight weight_reduce pos2 trim_from_left add_to_new_mum;
    (*trim mumj, update positions of mumj *)
    remove_from_pos2seed_tbl pos2seed_tbl1 seqNO oldpos jseedNO mum_tbl;
    add_to_pos2seed_tbl pos2seed_tbl1 seqNO newpos jseedNO newweight oldori mum_tbl;
    (*the weigh in pos2seed_tbl2 need to be updated*)
    modify_record_in_pos2seed_tbl pos2seed_tbl2 seqNO pos2 jseedNO newweight oldori;
    (*update position list of j_mum *)
    update_positions_by_seqNO jseedNO seqNO new_mi mum_tbl;
    (*no need to modify seed2pos_tbl if tabl1 is the right table, for
    * seed2pos_tbl just record left end of each seed. *)
    if (trim_from_left) then begin
    remove_from_seed2pos_tbl seed2pos_tbl seqNO newpos jseedNO;
    add_to_seed2pos_tbl seed2pos_tbl seqNO newpos jseedNO newweight oldori;
    end;
    (*add positions to the new mum_k*)
    if (add_to_new_mum=true) then begin
        let k_left,k_right = 
            if (trim_from_left) then oldpos,oldpos+weight_reduce-1
            else oldpos-weight_reduce+1,oldpos 
        in
        let k_pos = { sequence_NO = seqNO; left_end = k_left; right_end =
            k_right; orientation = oldori } in
        if ((Hashtbl.mem mum_tbl kseedNO)=false) then
            (*this is the first position of k, create the k_mum.*)
            Hashtbl.add mum_tbl kseedNO 
            {seedNO=kseedNO; positions=[k_pos]; size=1;  neighborhood_lst = []; 
            subsuming_pointer= -1; extendable = 2; mumscore = 0 ;} 
        else  (*update position list of k_mum *)
            update_positions_by_seqNO kseedNO seqNO k_pos mum_tbl
        ;
        if debug then print_mum kseedNO mum_tbl false true;
        (*update seed2pos table in any case*)
        add_to_seed2pos_tbl seed2pos_tbl seqNO k_left kseedNO weight_reduce oldori;
        let k_pos_tbl1,k_pos_tbl2 =
            if (trim_from_left) then k_left,k_right
            else k_right,k_left 
        in
        add_to_pos2seed_tbl pos2seed_tbl1 seqNO k_pos_tbl1 kseedNO weight_reduce oldori mum_tbl;
        add_to_pos2seed_tbl pos2seed_tbl2 seqNO k_pos_tbl2 kseedNO weight_reduce oldori mum_tbl;
    end;
    if debug then Printf.printf "end of update tables \n%!"
    
    
(*we only consider mums with extendable=true here*)  
(* if two mums are overlap with each other in one sequence, we trim the one with
* lower multiplicity, if there is a tie, we trim the one with shorter length, if
* there still a tie, pick any one.
* say mumA and mumB has match in sequence1 at (starta,enda) and (startb,endb) 
* General case is , starta < startb < enda < endb. If we decide to trim mumB in
* this case, new match of mumA in sequence1 remains the same, match of mumB in
* sequence1 become (enda+1,endb). 
* If mumB has another match in sequence2, like (startb2,endb2) with the same
* orientation. we also break this matche of mumB in two parts, one is 
* (startb2,startb2 + (enda-startb+1) ), make this match a new MUM, say mumC,
* the second part (startb2 + (enda-startb+1) +1, endb2) will be the new match of
* mumB. Do this for each sequence, add new matches to mumC. 
* The result : we have a shorter mumB, also have a new MUM, mumC.
* *)
(*mums smaller than min_weight will be discarded*)
let resolve_overlap_mum mum_tbl min_weight pos2seed_tbl_left pos2seed_tbl_right 
seed2pos_tbl debug debug_neighborhood =
    if debug then Printf.printf "resolve overlap mums with min_weight = %d\n%!" min_weight;
    Hashtbl.iter (fun key mum ->
        update_neighborhood key mum_tbl pos2seed_tbl_left
        pos2seed_tbl_right false
    ) mum_tbl;
    let keep_i multi_i multi_j weight_i weight_j =
        if (multi_i>multi_j)||((multi_i=multi_j)&&(weight_i>=weight_j)) then true
        else false
    in
    let get_three_pos trim_from_left record weight_reduce =
        if (trim_from_left) then 
            record.left_end, record.left_end+weight_reduce, record.right_end
        else 
            record.right_end,record.right_end - weight_reduce, record.left_end
    in
    let get_new_mi trim_from_left record newpos = 
        if (trim_from_left) then { record with left_end = newpos }
        else {record with right_end = newpos} 
    in
    (*when we trim mumj from left, we remove the old record for mumj pos2seed_tbl_left
    * and add the new one. for pos2seed_tbl_right, just update the old record
    * with new j_weight. So we mark tbl_left as tbl1 and tbl_right as tbl2.
    * When we trim mumj from right, it's another way around. I don't want to
    * write two "update table" functions for each direction, that why we have this function here. *)
    let update_table_trim_right_or_left trim_from_left add_to_new_mum 
        j_seedNO current_mi new_mi k_seedNO oldpos newpos pos2 new_jweight
        weight_reduce =
        let record = current_mi in
        if (trim_from_left) then
            update_tables mum_tbl 
            j_seedNO record.sequence_NO k_seedNO pos2seed_tbl_left 
            oldpos newpos pos2seed_tbl_right pos2 new_jweight 
            seed2pos_tbl new_mi record.orientation weight_reduce 
            true add_to_new_mum 
        else 
            update_tables mum_tbl 
            j_seedNO record.sequence_NO k_seedNO pos2seed_tbl_right
            oldpos newpos pos2seed_tbl_left pos2 new_jweight 
            seed2pos_tbl new_mi record.orientation weight_reduce 
            false add_to_new_mum
    in
    (* print_position2seedtbl pos2seed_tbl_left; *)
    let seedlst = ref [] in 
    Hashtbl.iter (fun key mum -> seedlst := key::!seedlst ) mum_tbl;
    let seedlst = List.rev !seedlst in
    List.iter (fun i_seedNO -> 
        let debug2 = false in
        if debug2 then Printf.printf "work on i_seed: %d\n%!" i_seedNO;
        let mi = try (Hashtbl.find mum_tbl i_seedNO)
            with | Not_found -> failwith "not found in resolve overlap 0" 
        in
        if debug2 then print_mum i_seedNO mum_tbl true true;
        if (mi.extendable=0)&&(mi.subsuming_pointer=(-1)) then begin 
        if (mi.neighborhood_lst=[]) then
            (*in case this mum is a new mum added by function "transpose_back" 
            * when we search area outside lcbs*)
            update_neighborhood i_seedNO mum_tbl pos2seed_tbl_left
            pos2seed_tbl_right false;
        let idx = ref 0 and sign = ref true in
        while (!sign)&&((List.length mi.neighborhood_lst)>0) do
            (*mi's positions could be different from last iteration with another
            * mj, that's why we need to load mi again.*)
            let mi = try (Hashtbl.find mum_tbl i_seedNO)
            with | Not_found -> failwith "not found in resolve overlap 0" in
            if debug2 then begin
                Printf.printf "i_mum is: %!"; print_mum i_seedNO mum_tbl true true;
            end;
            let (seqNO,j_seedNO,i_ori,j_ori,d) = List.nth mi.neighborhood_lst !idx in
            if debug2 then Printf.printf "idx=%d,i=%d,j=%d,seqNO=%d\n%!" 
            !idx i_seedNO j_seedNO seqNO;
            let mj = try (Hashtbl.find mum_tbl j_seedNO) 
            with | Not_found -> failwith "not found in resolve overlap 2" in
            if debug2 then begin 
                Printf.printf "j_mum is :%!"; print_mum j_seedNO mum_tbl false true;
            end;
            if (mj.extendable=0)&&(mi.extendable=0)&&(mj.subsuming_pointer=(-1)) 
            then begin
                let poslst_i = mi.positions in
                let position_i = get_position_by_seqNO poslst_i seqNO in
                let ileft = position_i.left_end and iright = position_i.right_end
                and iori = position_i.orientation  in
                let poslst_j = mj.positions in
                let position_j = get_position_by_seqNO poslst_j seqNO in
                let jleft = position_j.left_end and jright = position_j.right_end 
                and jori = position_j.orientation in
                let sizei = List.length poslst_i and sizej = List.length poslst_j in
                let leni= iright-ileft+1 and lenj= jright-jleft+1 in
                if debug2 then Printf.printf "i=%d(%d,%d) j=%d(%d,%d),d=%d\n%!" 
                        i_seedNO ileft iright j_seedNO jleft jright d;
                if (jleft > iright)|| (jright < ileft) then begin
                    idx := !idx + 1;
                    if debug2 then Printf.printf "idx ++ = %d\n%!" !idx;
                end
                else begin
                    idx := 0 ; (*neighborhood of i_mum will change, start over again*)
                    assert ((ileft-jleft)=d);
                    let i_seedNO,j_seedNO,poslst_j,position_j,mj,
                    ileft,iright,jleft,jright,jori,d = 
                        if (keep_i sizei sizej leni lenj) then begin
                            i_seedNO,j_seedNO,poslst_j,position_j,mj,
                            ileft,iright,jleft,jright,jori,d
                        end
                        else begin (*we keep j_mum,trim i_mum instead,*) 
                             j_seedNO,i_seedNO,poslst_i,position_i,mi,
                             jleft,jright,ileft,iright,iori,(0-d)
                        end
                    in
                    let k_seedNO = get_a_seedNO () in 
                    if debug2 then begin
                        Printf.printf "add new mum:%d, before trim mums \n%!" k_seedNO;
                        print_mum i_seedNO mum_tbl false false;
                        print_mum j_seedNO mum_tbl false false; end;
                    if (d<0) then begin
                        let new_jleft = iright + 1 in
                        let new_jweight = jright-new_jleft+1 in
                        let weight_reduce = iright - jleft + 1 in
                        if debug2 then Printf.printf "d<0,%!";
                        if (new_jweight>min_weight) then 
                            List.iter (fun record ->
                                let trim_from_left = 
                                    if (record.orientation=jori) then true else false 
                                in
                                let add_to_new_mum =
                                    if (seqNO<>record.sequence_NO) then true else false in
                                let oldpos,newpos,pos2 =
                                    get_three_pos trim_from_left record weight_reduce in
                                let new_mi = get_new_mi trim_from_left record newpos
                                in
                                update_table_trim_right_or_left trim_from_left add_to_new_mum 
                                j_seedNO record new_mi k_seedNO oldpos newpos pos2 new_jweight
                                weight_reduce;
                            ) poslst_j
                        else (*no need to carry j_seed or k_seed with us*) 
                            remove_seed2 j_seedNO seqNO pos2seed_tbl_left
                            pos2seed_tbl_right jleft jright mum_tbl k_seedNO;
                    end
                    else if (d>0) then begin
                        let new_jright = ileft-1 in
                        let new_jweight = new_jright-jleft+1 in
                        let weight_reduce = jright - ileft + 1 in
                        if debug2 then Printf.printf "d>0,%!";
                        if (new_jweight<=0) then
                            Printf.printf "ileft=%d,newj_right=%d,jleft=%d\n%!"
                            ileft new_jright jleft;
                        assert(new_jweight>0);
                        if (new_jweight>min_weight) then 
                        List.iter (fun record ->
                            let trim_from_left = 
                                if (record.orientation=jori) then false else true 
                            in
                            let add_to_new_mum =
                                if (seqNO<>record.sequence_NO) then true else false in
                            let oldpos,newpos,pos2 =
                                get_three_pos trim_from_left record weight_reduce in
                            let new_mi = 
                                get_new_mi trim_from_left record newpos in
                            update_table_trim_right_or_left trim_from_left add_to_new_mum 
                            j_seedNO record new_mi k_seedNO oldpos newpos pos2 new_jweight
                            weight_reduce;
                        ) poslst_j       
                        else 
                            remove_seed2 j_seedNO seqNO pos2seed_tbl_left
                            pos2seed_tbl_right jleft jright mum_tbl k_seedNO;
                    end; (* if d<0 else if d>0..*)
                    if debug2 then begin 
                        Printf.printf "after trim mums \n%!";
                        print_mum i_seedNO mum_tbl false true;
                        print_mum j_seedNO mum_tbl false true;
                        if (Hashtbl.mem mum_tbl k_seedNO) then 
                        print_mum k_seedNO mum_tbl false true;
                    end;
                    (*we only upate neighborhood when trim happens*)
                    Hashtbl.iter (fun key mum ->
                      update_neighborhood key mum_tbl pos2seed_tbl_left
                      pos2seed_tbl_right debug2
                    ) mum_tbl;
                end; (* if trimable -- mi and mj are overlaped*)
            end (*if (mj and mi.extendable=true)*)
            else
                idx := !idx +1;
            if (!idx >= ((List.length mi.neighborhood_lst)-1)) then 
                sign := false
            else sign := true;
        done; (*while sign*)
        end
        else () (*mi.extendable=false*)
    ) seedlst;
    let res = ref 0 in
    Hashtbl.iter (fun key record ->
        if (record.extendable = 0) then res := !res + 1
    ) mum_tbl;
    if debug then 
        Printf.printf "end of resolve overlap mums, extendable mum = %d\n %!" !res;
    !res
    
let update_score_for_each_mum mum_tbl in_seqlst = 
    Hashtbl.iter (fun key record ->
        if (record.extendable = 0) then begin
        let score = get_mum_score key mum_tbl in_seqlst in
        Hashtbl.remove mum_tbl key;
        Hashtbl.add mum_tbl key {record with mumscore = score };
        end
    ) mum_tbl

let search_inside_each_lcb lcbs input_seqlst mum_tbl min_weight = 
    Printf.printf "Search inside each lcb ,min_weight = %d\n%!" min_weight;
    let lcb_tbl = Hashtbl.create init_tb_size in
    let len = List.length input_seqlst in
    let seqlcbs = List.nth lcbs 1 in
    let init_size = init_tb_size in 
    for i = 0 to (len-1) do
        List.iter (fun seedNOlst ->
                let (left,right),_ = get_range_of_a_lcb seedNOlst i mum_tbl in
                Hashtbl.add lcb_tbl (seedNOlst,i) (left,right)
        )seqlcbs
    done;
    List.iter (fun seedNOlst ->
        Printf.printf "work on seedNOlst=%!";
        print_int_list seedNOlst;
        let (l,r) = Hashtbl.find lcb_tbl (seedNOlst,0) in 
        let res_seqlst = ref [[]] in
        if ((r-l+1)>min_weight) then begin
        for i = 0 to (len-1) do
            let (lefti,righti) = 
                try (Hashtbl.find lcb_tbl (seedNOlst,i)) with
                |Not_found -> failwith "not found1 inside_each_lcb"
            in
            let seqi = get_sub_seq (List.nth input_seqlst i) lefti (righti-lefti+1) in
            res_seqlst := seqi::!res_seqlst;
        done;
        let new_seq2seedNO_tbl = Hashtbl.create init_size in
        let new_seedNO2seq_tbl = Hashtbl.create init_size in
        let new_pos2seed_tbl_left = Hashtbl.create init_size in
        let new_pos2seed_tbl_right = Hashtbl.create init_size in
        let new_seed2pos_tbl = Hashtbl.create init_size in
        let lcb_seqlst = List.rev !res_seqlst in
        let shorted_seqlen = get_shortest_len lcb_seqlst in
        let seedlen = int_of_float ( ceil (log (float shorted_seqlen))) in
        let seedlen = if (seedlen mod 2)=0 then seedlen+1 else seedlen in 
        let seedlen = if (seedlen<5) then 5 else seedlen in
        let seedweight = build_seed_and_position_tbl lcb_seqlst seedlen 2 
        new_seq2seedNO_tbl new_seedNO2seq_tbl 
        new_pos2seed_tbl_left new_pos2seed_tbl_right 
        new_seed2pos_tbl false in
        Printf.printf "; seedweight = %d\n%!" seedweight;
        let new_mum_tbl = Hashtbl.create init_tb_size in
        build_local_mums new_mum_tbl new_seq2seedNO_tbl new_seedNO2seq_tbl 
        new_pos2seed_tbl_left new_pos2seed_tbl_right new_seed2pos_tbl false false;
        let _ = resolve_overlap_mum new_mum_tbl seedweight 
        new_pos2seed_tbl_left new_pos2seed_tbl_right new_seed2pos_tbl false false in
        Printf.printf 
        "search inside each lcb, check mum table after resolve_overlap_mum:\n%!";
        Hashtbl.iter 
        (fun seedNO mum -> print_mum seedNO new_mum_tbl false false) new_mum_tbl; 
        end
    ) seqlcbs

let search_outside_lcbs inner_lcbs mum_tbl
pos2seed_tbl_left pos2seed_tbl_right seed2pos_tbl
in_seqlst in_seq_size_lst  =
let debug = false in
    let lcb_tbl = Hashtbl.create init_tb_size in
    (*fill in empty lcb_tbl with current_lcbs and mum_tbl*)
    update_lcbs_range_and_score inner_lcbs mum_tbl lcb_tbl in_seqlst; 
    (*fprintf oc "inner loop, lcb are :\n%!";
    Hashtbl.iter (fun key record ->
     print_lcb_to_file oc record
    ) lcb_tbl; *)
    if debug then
        Printf.printf "\n Start working on sequence outside existing lcbs\n%!";
    (*this works even when the seq outside is shorter than seed's length*)
    let lcb_range_lstlst, seq_outside_lcbs_lst = 
        get_seq_outside_lcbs in_seqlst lcb_tbl
    in
    let current_seq_size_lst = List.fold_right (fun seq acc ->
        (List.length seq)::acc ) seq_outside_lcbs_lst [] 
    in
    if debug then begin
        Printf.printf "current_seq_size_lst :%!";
        print_int_list current_seq_size_lst;
    end;
    let shorted_seqlen = get_shortest_len seq_outside_lcbs_lst in
    let new_seedlen = int_of_float ( ceil (log (float shorted_seqlen))) in
    let new_seedlen = if (new_seedlen mod 2)=0 then new_seedlen+1 else new_seedlen in 
    let new_seedlen = if (new_seedlen<5) then 5 else new_seedlen in
    if debug then 
        Printf.printf "seedlen=%d,%!" new_seedlen;
    (*these hashtable are based on the sub-seq list*)
    let new_seq2seedNO_tbl = Hashtbl.create init_tb_size in
    let new_seedNO2seq_tbl = Hashtbl.create init_tb_size in
    let new_pos2seed_tbl_left = Hashtbl.create init_tb_size in
    let new_pos2seed_tbl_right = Hashtbl.create init_tb_size in
    let new_seed2pos_tbl = Hashtbl.create init_tb_size in
    let new_seedweight = 
        build_seed_and_position_tbl seq_outside_lcbs_lst new_seedlen 
        !minimum_seed_num new_seq2seedNO_tbl new_seedNO2seq_tbl 
        new_pos2seed_tbl_left new_pos2seed_tbl_right 
        new_seed2pos_tbl false in
    if debug then 
        Printf.printf "seedweight=%d\n%!" new_seedweight;
    let new_mum_tbl = Hashtbl.create init_tb_size in
    build_local_mums new_mum_tbl new_seq2seedNO_tbl new_seedNO2seq_tbl 
    new_pos2seed_tbl_left new_pos2seed_tbl_right new_seed2pos_tbl false false;
    (*transposed new mums and add them to the old mum_tbl we have*)
    if debug then
        Printf.printf "\n ============= Transpose back: ================\n%!";
    let res_mum_tbl = Hashtbl.copy mum_tbl in
    let res_pos2seed_tbl_left = Hashtbl.copy pos2seed_tbl_left in
    let res_pos2seed_tbl_right = Hashtbl.copy pos2seed_tbl_right in
    let res_seed2pos_tbl = Hashtbl.copy seed2pos_tbl in
    transpose_back lcb_range_lstlst new_mum_tbl new_seedNO2seq_tbl
    res_mum_tbl res_pos2seed_tbl_left res_pos2seed_tbl_right res_seed2pos_tbl;
    if debug then 
        Printf.printf "\n ============= Resolve overlap mum: =============\n%!";
    let num_of_mums = resolve_overlap_mum res_mum_tbl (new_seedweight-1)
    res_pos2seed_tbl_left res_pos2seed_tbl_right res_seed2pos_tbl false false in
    update_score_for_each_mum res_mum_tbl in_seqlst;
    (*get the new weight of lcbs*)
    let new_seedNOlstlst = get_mum_lst_for_each_seq res_mum_tbl res_pos2seed_tbl_left
    (List.length in_seqlst) in_seq_size_lst in
    if debug then begin
        Printf.printf "new seedNOlstlst is :\n%!";
        print_int_lstlst3  new_seedNOlstlst;
    end;
    let new_lcbs,new_weight = build_LCBs new_seedNOlstlst in
    let previous_lcb_tbl = lcb_tbl in 
    (*this is the lcb_tbl before adding any new mums*)
    new_lcbs,new_weight,previous_lcb_tbl,res_mum_tbl, 
    res_pos2seed_tbl_left, res_pos2seed_tbl_right, res_seed2pos_tbl, num_of_mums

(*main function here*)
let create_lcb_tbl in_seqlst min_lcb_ratio min_seed_num min_bk_penalty =
    minimum_seed_num := min_seed_num ; 
    minimum_lcb_ratio := min_lcb_ratio ;
    min_break_point_penalty := min_bk_penalty ;
    let debug = false in
    seedNO_available_arr := Array.make init_seed_size 1;
    (*output result to file ...
    let outfile = "outfile.txt" in
    let oc = open_out outfile in*)
    let in_seqlst_size = (List.length in_seqlst) in
    let in_seq_size_lst = List.fold_right (fun seq acc ->
         (List.length seq)::acc ) in_seqlst [] 
    in
    let shorted_seqlen = get_shortest_len in_seqlst in
    let seedlen = int_of_float ( ceil (log (float shorted_seqlen))) in
    let seedlen = if (seedlen mod 2)=0 then seedlen+1 else seedlen in 
    let seedlen = if (seedlen<5) then 5 else seedlen in
    if debug then Printf.printf "seedlen=%d\n%!" seedlen;
    let init_size = init_tb_size in 
    let seq2seedNO_tbl = Hashtbl.create init_size in
    let seedNO2seq_tbl = Hashtbl.create init_size in
    (* Hashtbl { pos(seqNO,idx) -> seedNO,seedweight,orientation }, 
    * seqNO is the NO of sequence, idx is the position of that
    * sequence.  pos2seed_tbl_left records the left_end of each seed,
    * pos2seed_tbl_right records the right_end of each seed *)
    let pos2seed_tbl_left = Hashtbl.create init_size in
    let pos2seed_tbl_right = Hashtbl.create init_size in
    let seed2pos_tbl = Hashtbl.create init_size in
    let seedweight = build_seed_and_position_tbl in_seqlst seedlen 2 
    seq2seedNO_tbl seedNO2seq_tbl pos2seed_tbl_left pos2seed_tbl_right
    seed2pos_tbl false in
    let mum_tbl = Hashtbl.create init_tb_size in
    build_local_mums mum_tbl seq2seedNO_tbl seedNO2seq_tbl 
    pos2seed_tbl_left pos2seed_tbl_right seed2pos_tbl false false;
    if debug then
        Printf.printf "++++++++ init seedweight=%d, end of mum table, mumtbl size=%d\
    (un-extendable included)\n%!" seedweight (Hashtbl.length mum_tbl);
    (* Hashtbl.iter (fun key record -> print_mum key !mum_tbl false true) !mum_tbl; *)
    (*print_position2seedtbl pos2seed_tbl_left;
    print_position2seedtbl pos2seed_tbl_right;*)
    let init_num_mums = resolve_overlap_mum mum_tbl (seedweight-1) 
    pos2seed_tbl_left pos2seed_tbl_right seed2pos_tbl false false in
    update_score_for_each_mum mum_tbl in_seqlst;
    if debug then
        Printf.printf "++++++++++ check mum table after resolve_overlap_mum +++++++ \n%!";
    (*Hashtbl.iter (fun seedNO mum -> print_mum seedNO mum_tbl false false) mum_tbl;*)
    let seedNOlstlst = get_mum_lst_for_each_seq mum_tbl pos2seed_tbl_left
    in_seqlst_size in_seq_size_lst in 
    if debug then begin
        Printf.printf "init seedNOlstlst is :\n%!";
        print_int_lstlst3  seedNOlstlst;
    end;
    let init_lcbs,init_weight = build_LCBs seedNOlstlst in
    (*we might have lcb with high W low R here, these will be removed, so the
    * init_weight could be higher than it should be *)
    if debug then Printf.printf "init lcb weight = %d\n%! " init_weight;
    (*init inner tbl*)
    let inner_old_weight = ref init_weight in
    let inner_lcbs = ref init_lcbs in
    let inner_lcb_tbl = ref (Hashtbl.create init_tb_size) in
    let inner_mum_tbl = ref mum_tbl in
    let inner_pos2seed_tbl_left = ref pos2seed_tbl_left in
    let inner_pos2seed_tbl_right = ref pos2seed_tbl_right in
    let inner_seed2pos_tbl = ref seed2pos_tbl in
    (*init outer tbl*)
    let outer_sign = ref true in
    let outer_old_weight = ref 0 in
    let outer_lcbs = ref init_lcbs in
    let outer_lcb_tbl = ref (Hashtbl.create init_tb_size) in
    let outer_mum_tbl = ref mum_tbl in
    let current_num_of_mums = ref init_num_mums in
    while (*(init_weight>0)&&*)(!outer_sign) do
        if debug then
        Printf.printf "\n begin of outer while, old_weight = %d\n%!" !outer_old_weight;
    (*****work on regions outside of lcbs******)
    (* we do need to update inner lcbs with the outer one, so we can work on new
    * range of subsequence, but inner_old_weight should remain the same, 
    * since inner_weight is the weight before any low-weight lcb removing *)
        inner_lcbs := !outer_lcbs;
        inner_mum_tbl := !outer_mum_tbl; (*this one is not necessary*)
        let inner_sign = ref true in
        while (!inner_sign) do
            if debug then
                Printf.printf "\n ----------- Inner loop -------------- \n%!";
            let inner_new_lcbs,inner_new_weight,previous_lcb_tbl,res_mum_tbl,
            res_pos2seed_tbl_left, res_pos2seed_tbl_right, res_seed2pos_tbl, 
            num_of_mums = 
            search_outside_lcbs !inner_lcbs !inner_mum_tbl !inner_pos2seed_tbl_left 
            !inner_pos2seed_tbl_right !inner_seed2pos_tbl in_seqlst in_seq_size_lst 
            in
            if (num_of_mums>0)&&(!inner_old_weight < inner_new_weight) then begin
                current_num_of_mums := num_of_mums;
                if debug then
                    Printf.printf "we found more lcbs (%d>%d), update num_of_mums to %d,\
                continue with inner loop\n%!" inner_new_weight !inner_old_weight 
                !current_num_of_mums;
                inner_sign := true;
                inner_old_weight := inner_new_weight;
                inner_lcbs := inner_new_lcbs;
                inner_mum_tbl := res_mum_tbl;
                inner_pos2seed_tbl_left := res_pos2seed_tbl_left;
                inner_pos2seed_tbl_right := res_pos2seed_tbl_right;
                inner_seed2pos_tbl := res_seed2pos_tbl;
            end
            else begin
                inner_sign := false;
                inner_lcb_tbl := previous_lcb_tbl;
                if debug then Printf.printf "no improve on lcb weight (%d<=%d) \
                get out of inner loop\n%!" inner_new_weight !inner_old_weight ;
            end;
            if debug then
                Printf.printf "end of inner while, current_weight=%d\n%!" inner_new_weight;
        done; 
        (*end of inner while loop*)
        if debug then
            Printf.printf "\n --------  Outer loop, remove light lcbs ----- \n%!";
        let new_outer_lcbs, new_outer_weight, new_outer_lcb_tbl = 
            remove_light_weight_lcbs !inner_lcbs !inner_lcb_tbl !inner_mum_tbl 
            in_seqlst in_seq_size_lst !current_num_of_mums in
        (*we could still have low W lcb here -- if remove them will create high W
        * low R lcbs, "remove_light_weight_lcbs" won't do it, get rid of those
        * lcbs here
        * Note. mauve keep the last low W lcb, we will do the same*)
        Hashtbl.iter (fun key record ->
            if ((Hashtbl.length new_outer_lcb_tbl)>1)&&(List.length key)<minimum_lcb_weight 
            then  Hashtbl.remove new_outer_lcb_tbl key;
        ) new_outer_lcb_tbl;
        let new_outer_weight = Hashtbl.length new_outer_lcb_tbl in
        (* remove corresponding item in new_outer_lcbs*)
        let new_outer_lcbs = 
                List.map (fun lcblst ->
                List.filter (fun record ->
                 Hashtbl.mem new_outer_lcb_tbl record
                ) lcblst
                ) new_outer_lcbs 
        in
        if debug&&(new_outer_weight=0) then 
            Printf.printf "remove lcb did not give us any high W lcbs\n%!";
        if debug then begin
            Printf.printf "\n after remove light lcbs, we have :\n%!";
            Hashtbl.iter (fun key record ->
           (* if (List.length key) >= minimum_lcb_weight then*) 
                print_lcb record 
            ) new_outer_lcb_tbl;
        end;
        (*
        fprintf oc "outer loop, lcb :\n%!";
        Hashtbl.iter (fun key record ->
             print_lcb_to_file oc record
        ) new_outer_lcb_tbl;*)
        if (!outer_old_weight < new_outer_weight) then begin
            if debug then 
                Printf.printf "\n we have a new lcb weight = %d>old one=%d, \
            continue with outer loop again\n%!" new_outer_weight !outer_old_weight;
            outer_sign := true;
            outer_old_weight := new_outer_weight;
            outer_lcbs := new_outer_lcbs;
            outer_lcb_tbl := new_outer_lcb_tbl;
            outer_mum_tbl := !inner_mum_tbl;
        end
        else begin
            outer_sign := false;
            if debug then 
            Printf.printf "\n no improve on lcb weight (%d>=%d),\
            get out of outer loop\n%!" !outer_old_weight new_outer_weight;
        end
    done; (*end of outer while loop*)
    if debug then
        Printf.printf "init_weight=%d,outer old weight = %d, outer_lcb_tbl len=%d\n%!"
        init_weight !outer_old_weight (Hashtbl.length !outer_lcb_tbl);
    (*if no qualified lcb are found, just make the whole sequence as one lcb*)
    if (Hashtbl.length !outer_lcb_tbl)=0 then begin 
        if debug then Printf.printf "we didn't find any qualified lcb\n%!";
        (*let lcb_range_lstlst = get_lcbs_range !outer_lcb_tbl in_seqlst true in
        print_lcbs_range lcb_range_lstlst;*)
        let code_list = [[1];[1]] in
        let range0 =  (0,(List.hd in_seq_size_lst)-1)
        and range1 =  (0,(List.nth in_seq_size_lst 1)-1)
        in
        let m0 = {sequence_NO = 0; left_end = 0; right_end = (List.hd
        in_seq_size_lst)-1; orientation = 1}
        and m1 = {sequence_NO = 1; left_end =0; right_end = (List.nth
        in_seq_size_lst 1)-1; orientation = 1}
        in
        let full_range_lstlst = [ [range0];[range1] ] in
        let single_lcb = 
            {
                seedNOlst = [1];
                range_lst = [m0;m1];
                score = 0;
                ratio = 0.0;
                ref_code = 1;
            }
        in
        let outer_lcb_tbl = !outer_lcb_tbl in
        Hashtbl.clear outer_lcb_tbl;(*not necessary,outer_lcb_tbl should be empty here*)
        Hashtbl.add outer_lcb_tbl [1] single_lcb;
        let outer_lcbs = [[[1]];[[1]]] in
        outer_lcb_tbl,outer_lcbs,code_list,full_range_lstlst
    end
    (* we do find some high W high R lcbs, return them to aliMap  *)
    else begin
        (*fprintf oc "Final result :\n%!";
        Hashtbl.iter (fun key record ->     print_lcb_to_file oc record
        ) !outer_lcb_tbl;  close_out oc;*)
        let outer_lcbs = !outer_lcbs in
        let outer_lcb_tbl = !outer_lcb_tbl in
        if debug then
            Hashtbl.iter (fun key record ->
             print_lcb record
        ) outer_lcb_tbl; 
        let lcb_range_lstlst = get_lcbs_range outer_lcb_tbl in_seqlst true in
        let full_range_lstlst = 
            get_full_range_lstlst lcb_range_lstlst in_seq_size_lst in
        if debug then begin
            Printf.printf "here is LCB range lstlst %!";
            print_lcbs_range lcb_range_lstlst;
            Printf.printf "here is the full range lstlst (includes non-LCB blocks %!";
            print_lcbs_range full_range_lstlst;
        end;
        let refcode = ref 1 in
        Hashtbl.iter (fun key record -> 
            refcode := !refcode;
            update_lcb_ref_code outer_lcb_tbl key (!refcode);
            refcode := !refcode + 1;
        ) outer_lcb_tbl;
        let seqNO = ref (-1) in
        let code_list = List.map (fun lcblst -> 
            seqNO := !seqNO + 1;
            let tmplst = List.map (fun key ->
                let record = try (Hashtbl.find outer_lcb_tbl (get_abs_lst key)) 
                with |Not_found -> failwith "not found in creating code list" in
                let mi_lst = record.range_lst in
                let mi = get_position_by_seqNO mi_lst !seqNO in
                record.ref_code * mi.orientation ) lcblst 
            in
            if debug then print_int_list tmplst;
            tmplst
        ) outer_lcbs in
        (*clear up*)
        seedNO_available_arr := Array.make init_seed_size 1;
        (*return lcb table, lcbs, code list and range list*)
        outer_lcb_tbl,outer_lcbs,code_list,full_range_lstlst
    end

    

