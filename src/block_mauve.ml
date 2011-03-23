open Printf
let error_user_message format = Printf.ksprintf (Status.user_message Status.Error) format
let info_user_message format = Printf.ksprintf (Status.user_message Status.Information) format

(* W = weight, R = ratio *)

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
    Hashtbl.add res 13
    [
        [1;2;3;5;7;8;11;14;15;17;19;20;21];
        [1;2;3;5;8;9;11;13;14;17;19;20;21];
        (* 1111**1**1*1*1**1**1111 this one sucks *)
    ];
    Hashtbl.add res 15
    [
        [1;2;3;4;6;8;9;12;15;16;18;20;21;22;23]
    ];
    Hashtbl.add res 19
    [
        [1;2;3;4;6;7;9;10;11;13;15;16;17;19;20;22;23;24;25]
    ];
    Hashtbl.add res 21
    [
        [1;2;3;4;5;7;8;9;11;12;14;16;17;19;20;21;23;24;25;26;27]
    ];
    (*add 15 19 21 later*)
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
}

type b_tree = 
    | Leaf of int array * mum
    | Node of int array * b_tree * b_tree

type lcb = {
    seedNOlst : int list; (*each lcb is a list of mums*)
    range_lst : m_i list; (*lcb_range we need is exactly the same truct as m_i,
                           no need to create a new type here
    we sort it in update_lcbs_range_and_score, in decrease range size order*)
    score : int; (*score between subsequence contained by range_lst*)
    ratio : float; (* [score/length] of seq in this lcb*)
    ref_code : int; (*just a code for bk/rearr purpose*)
    avg_range_len : int; (*average of range length from range_lst *)
}

let empty_lcb = {seedNOlst=[];range_lst=[];score=0;ratio=0.;ref_code=0;avg_range_len=0;}

(** for a single lcb, if the coverR = (length of subseq in lcb/total length of seq)
 * is lower than this, it's a light-weight lcb*)
let minimum_cover_ratio = ref 0.0020
(** for all lcbs. if they don't cover more than this of input sequence, adjust
* parameters like minimum_cover_ratio and minimum_lcb_ratio to find more matches  
*)
let minimum_full_cover_ratio = ref 0.4
 
(*we stop looking for more lcbs when cover ratio is bigger than this number*) 
let maximum_cover_ratio = ref 0.9
(*we should stop looking when the improvement is too small*)
let minimum_improvement = ref 0.01
(*lcb with ratio smaller than this is considered low ratio lcb.*)
let minimum_lcb_ratio = ref 30.0
(** we lower the minimum_lcb_ratio when initial step could not find enough
* match(minimum_full_cover_ratio). we stop doing that until we have enough match, 
* or minimum_lcb_ratio is lower than this number*)
let stop_lcb_ratio = ref 5.

(*minimum value of break point penalty*)
let min_break_point_penalty = ref 4000 

(*lcb with lower weight should be considered removing*)
let minimum_lcb_weight = 2 
let init_tb_size = 50
let init_seed_size = 50
let max_gap_num = 0  (*the w *)
let seedNO_available_arr = ref (Array.make init_seed_size 1)

(** some operation on int/int list/int array/etc.. we gonna need *)
let to_ori_code code =
    assert(code<>0);
    if( code mod 2 == 0) then -(code/2)  else (code+1)/2 

let from_ori_code code =
    assert(code<>0);
    if (code<0) then (abs code)*2
    else (code*2-1)
    
let bigger_int a b =  if (a>b) then a else b

let get_shorter_len seqarr =
    let len0 = Array.length seqarr.(0) 
    and len1 = Array.length seqarr.(1) in
    if len0<len1 then len0 else len1

(*return average of a int list, return float*)
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

(*if lcb cover ratio is too small, we may need to align super big chunk of
* sequences after mauve, user might want to adjust pameter like 
* minimum_cover_ratio and minimum_lcb_ratio to find more matches*)
let get_adjusted_parameter full_covR in_seq_size_lst minlcbR mincovR = 
    let e = 2.718 in 
    let divided = 
        if full_covR > 0. then 
            !minimum_full_cover_ratio /. full_covR 
        else e *. e (*if we didn't find any good matches before, cut the
        requirement by half*)
    in
    Printf.printf "adjust parameter,divided = %f\n%!" divided;
    if divided <= 1. then 
        minlcbR,mincovR
    else if divided <= e then 
        minlcbR *. (0.9), mincovR *. (0.9)
    else begin
        minlcbR /. (log divided), mincovR *. 0.9
    end


let tiny_improvement newresult oldresult =
    if (oldresult = 0.) then begin
        Printf.printf "we didn't find any good matches before, any improvment \
        this round is good to us\n%!";
        false
    end
    else begin
        if (abs_float (newresult -. oldresult) ) /. oldresult < (!minimum_improvement) 
        then begin
            Printf.printf " tiny improvement, no need to continue\n%!";
            true
        end
        else false;
    end

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

let print_int_arr2 arr = 
    Array.iteri (fun idx item -> Printf.printf "[%d]:%d,%!"  idx item) arr;
    Printf.printf "\n%!"

let print_int_arr arr =
    Array.iter (fun item -> 
        if item=Utl.large_int then Printf.printf "  L%!"
        else 
            Printf.printf "%3i%!" item) arr;
    Printf.printf "\n%!"

let print_int_matrix m =
    Array.iter (fun arr ->
            print_int_arr2 arr ;
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
    List.iter (fun (seqNO, j,i_ori,j_ori,dl,dr) ->
        Printf.printf "{seqNO.%d, seed#%d, ori:%d/%d,d:%d,%d }\n%!" seqNO j i_ori
        j_ori dl dr
    ) in_lst

let print_mum in_mum  print_neighborhood print_unextendable =
    (*let in_mum = try (Hashtbl.find mum_tbl in_SeedNO) with 
    | Not_found -> failwith ("Print MUMs, not found err ") in*)
    if (in_mum.extendable = 0) ||
    ((in_mum.extendable !=0 )&&(print_unextendable=true) ) then begin
        Printf.printf "MUM#%d\n%!" in_mum.seedNO;
        Printf.printf " mumkey : %d\n%!" in_mum.mumkey;
        Printf.printf " mumseq = %!";
        print_int_arr in_mum.mumseq; 
        Printf.printf " size: %d ,\n" in_mum.size;
        Printf.printf " subsuming_pointer(is subsumed by) = %d ,%!" in_mum.subsuming_pointer; 
        Printf.printf " extendable = %d \n%!" in_mum.extendable;
        Printf.printf " priority_lst = %!";
        print_int_list in_mum.priority_lst;
        Printf.printf " positions: \n%!";
        print_pos_list in_mum.positions;
        if print_neighborhood then begin
        Printf.printf " neighborhood list : \n%!";
        print_neighborhood_list in_mum.neighborhood_lst;
        end;
    end

let print_lcb x = 
    Printf.printf "LCB.%!";
    print_int_list x.seedNOlst; Printf.printf "\n%!";
    print_pos_list x.range_lst;
    Printf.printf "score = %d;%!" x.score;
    Printf.printf "ratio = %f;%!" x.ratio;
    Printf.printf "avg range = %d;%!" x.avg_range_len;
    Printf.printf "ref code = %d\n%!" x.ref_code

let print_lcb_to_file oc x = 
    fprintf oc "LCB.%!";
    print_int_list_to_file oc x.seedNOlst; fprintf oc "\n%!";
    print_pos_list_to_file oc x.range_lst;
    fprintf oc "score = %d\n%!" x.score;
    fprintf oc "ratio = %f\n%!" x.ratio 

(***************** function deal with sub_seq start *********************)
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

let get_sub_seq2 seq leftend size = 
    Array.sub seq leftend size 

(****************** function deal with sub_seq end ***********************)

(***************** binary tree function begins ****************************)

let create_btree seedseq newmum = Leaf (seedseq,newmum)

(*add new mum_node to b_tree, if same mum already exist, replace it with new one
* if other mum with same match sequence already exist, do nothing, return false*)
let add_to_btree seedseq newmum old_bt =
    let debug = false in
    if debug then  Printf.printf "add mum#%d to btree\n%!" newmum.seedNO;
    let sign = ref true in
    let rec insert_node sub_bt key mum =
        match sub_bt with
        | Node (middle_key, left_bt, right_bt) ->
                let insert_dir = compare middle_key key in
                if insert_dir=1 then (*middle_key>key*)
                    let new_left_bt = insert_node left_bt key mum in
                    Node (middle_key, new_left_bt, right_bt)
                else 
                    let new_right_bt = insert_node right_bt key mum in
                    Node (middle_key,
                    left_bt,
                    new_right_bt)
        | Leaf (leafkey, leafmum) ->
                assert ((Array.length leafkey)>0);
                let insert_dir = compare key leafkey in
                if debug then begin
                    Utl.printIntArr key;
                    Utl.printIntArr leafkey;
                    Printf.printf " insert_dir=%d,leafmum.seedNO=%d\n%!" 
                    insert_dir leafmum.seedNO;
                end;
                if insert_dir=1 then (*key>leafkey*)
                    Node (key, 
                          Leaf (leafkey, leafmum),
                          Leaf (key,mum))
                else if insert_dir=(-1) then
                    Node (leafkey, 
                          Leaf (key,mum),
                          Leaf (leafkey, leafmum))
                else 
                    if newmum.seedNO=leafmum.seedNO then  
                        Leaf (leafkey,newmum)
                    else begin
                        sign := false;
                        Leaf (leafkey,leafmum)
                    end
    in
    let res_bt = insert_node old_bt seedseq newmum in
    res_bt,!sign

let search_in_btree key bt = (*mumseq is the key inside b_tree*)
    let rec search_node sub_bt key =
        match sub_bt with
        | Node (middle_key, left_bt, right_bt) ->
                let search_dir = compare middle_key key in
                if search_dir=1 then(*key<middle_key*)
                    search_node left_bt key
                else 
                    search_node right_bt key
        | Leaf (leafkey, leafmum) ->
                if (compare leafkey key)=0 then leafmum
                else begin
                        Printf.printf "we are searching for :%!";
                        print_int_arr key;
                        Printf.printf "search failed, we reach\n:%!";
                        print_mum leafmum false true;
                    failwith "search reach a dead end in btree";
                end
    in
    search_node bt key

let remove_from_btree seedseq bt = (*key in b_tree is seedseq*)
    let rec remove_node sub_bt key = 
        match sub_bt with 
        | Node (middle_key, left_bt, right_bt) ->
                let remove_dir = compare middle_key key in
                if remove_dir=1 then (*key<middle_key*)
                    let new_left_child = remove_node left_bt key in
                    (match new_left_child with
                    | Leaf ([||],_) ->  right_bt
                    | x ->  
                            Node (middle_key,x,right_bt))
                else 
                    let new_right_child = remove_node right_bt key in
                    (match new_right_child with
                    | Leaf ([||],_) -> left_bt
                    | x -> 
                            Node (middle_key,left_bt, x))
        | Leaf (leafkey, leafmum) ->
                if (compare leafkey key)=0 then  Leaf ([||],leafmum)
                else begin
                    Printf.printf "we are searching for :%!";
                        print_int_arr key;
                        Printf.printf "search failed, we reach\n:%!";
                        print_mum leafmum false false;
                    failwith "remove node reach a dead end in btree"
                end
    in
    remove_node bt seedseq 
    
let iter_b_tree bt debug f =
    let rec iter_sub_tree sub_bt f =
        match sub_bt with
        | Node (middle_key, left_bt, right_bt) ->
                if debug then begin
                    Printf.printf "{ mkey=%!";
                    print_int_arr middle_key;
                    Printf.printf "\n ( %!";
                end;
                iter_sub_tree left_bt f;
                iter_sub_tree right_bt f;
                if debug then Printf.printf " ) }\n%!";
        | Leaf (leafkey, leafmum) -> f leafmum 
    in
    iter_sub_tree bt f


let print_b_tree bt debug =
    let ourf leafmum = 
        if debug then print_mum leafmum false true
        else 
            Printf.printf " MUM#.%d %!" leafmum.seedNO
    in
    iter_b_tree bt true ourf

(******************** binary tree function ends *********************)

(**************** b_tree mum_tbl function starts *********************)

let print_mumtbl mum_tbl debug = 
    Hashtbl.iter (fun key record ->
        Printf.printf "mumkey = %d \n%!" key;
        print_b_tree record debug
    ) mum_tbl


(*get_mumkey_from_milst returns the new key for mum_tbl with positions list of a mum*)
let get_mumkey_from_milst milst =
    let headmi = List.find (fun record -> record.sequence_NO=0) milst in
    let s1 = headmi.left_end in
     List.fold_left (fun acc mi ->
        acc+abs(mi.left_end-s1)
    ) 0 milst 

(*get_mumkey_from_poslst and get_mumkey_and_seedseq_from_seed2postbl returns the
* key of mum_tbl with poslst from seed2pos_tbl*)
let get_mumkey_from_poslst poslst = 
    let _,s1,_ = List.hd poslst in
    List.fold_left (fun acc (seqNO,pos,_) ->
        acc+abs(pos-s1)
    ) 0 poslst

let get_mumkey_and_seedseq_from_seed2postbl seedNO seed2pos_tbl =
    let poslst,seedseq,_ = try (Hashtbl.find seed2pos_tbl seedNO) with
    | Not_found -> 
        begin 
            Printf.printf "seedNO=%d\n%!" seedNO;
            failwith "not found, get key and seedseq";
        end
    in
    get_mumkey_from_poslst poslst,seedseq

let get_mum_from_mumtbl seedNO mum_tbl seed2pos_tbl =
    let mumkey,seedseq = 
        get_mumkey_and_seedseq_from_seed2postbl seedNO seed2pos_tbl
    in
    let bt = try (Hashtbl.find mum_tbl mumkey) with
    | Not_found ->
            begin
                Printf.printf "seedNO=%d,mumkey=%d\n%!" seedNO mumkey;
            failwith "not found, get mum from mumtbl" 
            end 
    in
    let res = search_in_btree seedseq bt in
    if res.seedNO<>seedNO then begin
        let _,seedseq1,w1 = Hashtbl.find seed2pos_tbl seedNO 
        and _,seedseq2,w2 = Hashtbl.find seed2pos_tbl res.seedNO in
        Printf.printf "seed1=%d,w=%d,%!" seedNO w1;
        print_int_arr seedseq1;
        Printf.printf "; seed2=%d,w=%d,%!" res.seedNO w2;
        print_int_arr seedseq2;
        Printf.printf "\n%!";
        assert(false);
    end;
    res
        

(*add newmum to mum_tbl, create a b_tree if there is no entry of newmum.mumkey,
* if some b_tree already exist, add the new node to it. if old record of this
* mum exist in the that b_tree, replace it with the new one*)    
let add_mum_to_mumtbl newmum mum_tbl =
    let debug = false in
    if debug then 
        Printf.printf "add seed#%d with mumkey=%d to mumtbl\n%!"
        newmum.seedNO newmum.mumkey;
    let mumkey = newmum.mumkey in
    if (Hashtbl.mem mum_tbl mumkey) then begin
        let old_bt = Hashtbl.find mum_tbl mumkey in
        let new_bt,sign_newbt = add_to_btree newmum.mumseq newmum old_bt in
        if sign_newbt then Hashtbl.replace mum_tbl mumkey new_bt;
        if debug then Printf.printf "sign_newbt=%b\n%!" sign_newbt;
        sign_newbt
    end
    else begin
        let bt = create_btree newmum.mumseq newmum in
        Hashtbl.add mum_tbl mumkey bt;
        true
    end

let remove_mum_from_mumtbl mum2remove mum_tbl =
    let debug = false in
    if debug then Printf.printf "remove mum#%d from mumtbl\n%!"
    mum2remove.seedNO;
    let mumkey = mum2remove.mumkey in
    let bt = try (Hashtbl.find mum_tbl mumkey) 
    with |Not_found -> failwith "not found, remove mum from mumtbl" in
    match bt with 
    | Leaf (leafkey,leafmum) ->
            Hashtbl.remove mum_tbl mumkey
    | _ -> 
        let treekey = mum2remove.mumseq in
        let newbt = remove_from_btree treekey bt in
        match newbt with 
        | Leaf ([||],leafmum) ->
                Hashtbl.remove mum_tbl mumkey
        | _ -> Hashtbl.replace mum_tbl mumkey newbt
    


(*update_mum_to_mumtbl, update mum_tbl with newmum. 
* In case the key for mum_tbl changes -- like when we change mum.positions, 
* oldkey should be passed to remove the node from old tree. 
* if we have a new position for some mum, but "add_mum_to_mumtbl" return 
* a false sign, this means we already have the same match. *)    
let update_mum_to_mumtbl (oldmum:mum option) newmum mum_tbl debug =
    match oldmum with
    | Some oldmum -> 
            if debug then Printf.printf 
            "update_mum_to_mumtbl: remove oldmum#.%d,mumkey=%d;%!" 
            oldmum.seedNO oldmum.mumkey;
            remove_mum_from_mumtbl oldmum mum_tbl;
            if debug then Printf.printf "add mum,mumkey=%d\n%!" newmum.mumkey;
            let sign_newmum = add_mum_to_mumtbl newmum mum_tbl in
            sign_newmum
    | None -> 
            let sign_newmum = add_mum_to_mumtbl newmum mum_tbl in
            if debug then 
                Printf.printf "update_mum_to_mumtbl,sign_newmum=%b\n%!" sign_newmum;
            sign_newmum

(***************** b_tree mum_tbl function ends *********************)


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
    Hashtbl.iter (fun key (record,seedseq,seedweight) ->
        Printf.printf "seed#%d, { %!" key;
        List.iter (fun (seqNO,idx,dir)->
            Printf.printf " (%d,%d,%d);%!" seqNO idx dir
        )record;
        Printf.printf " }\n%!"
    ) in_tbl


(*print_pos_list2 print out one record of pos2seed_tbl, also the mum info*)
let print_pos_list2 seqNO pos in_lst mum_tbl seed2pos_tbl printmum =
    Printf.printf "check record on (%d,%d) of pos2seed tbl(size=%d):\n%!"
    seqNO pos (List.length in_lst);
    List.iter (fun (seedNO,weight,orientation) ->
        Printf.printf "[seed: %d, weight: %d, ori: %d ] \n%!"
        seedNO weight orientation;
        if printmum then begin
        let mum2print = get_mum_from_mumtbl seedNO mum_tbl seed2pos_tbl in
        print_mum mum2print false true;
        end;
    ) in_lst;
    Printf.printf "end of record in pos2seed tbl \n%!"


let print_position2seedtbl position2seed_tbl somepos = 
    let print_record record_lst seqNO pos= 
        Printf.printf "pos.(%d,%d) ==> " seqNO pos;
        List.iter (fun (seed,seed_weight,dir) ->
        Printf.printf "seed#%d,weight=%d,dir=%d; %!" seed seed_weight dir;
        ) record_lst;
        Printf.printf "\n%!"
    in
    (match somepos with
    | Some (seqNO,pos) ->
            if (Hashtbl.mem position2seed_tbl (seqNO,pos) ) then
                let record_lst = 
                    Hashtbl.find position2seed_tbl (seqNO,pos)
                in
                print_record record_lst seqNO pos
            else
                Printf.printf "there is no (%d,%d) in pos2seed tbl\n%!"
                seqNO pos 
    | None ->
        Hashtbl.iter (fun (seqNO,pos) record_lst ->
            print_record record_lst seqNO pos
        ) position2seed_tbl;
        Printf.printf "\n%!" )

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

let break_code in_code =
    let resarr = Array.make 4 0 in
    let tmp_code = ref in_code in
    for i = 3 downto 0 do
        let base = int_of_float ( 2. ** (float_of_int i) )  in
        if (!tmp_code >= base) then begin
            resarr.(i)<-1;
            tmp_code := !tmp_code - base;
        end;
    done;
    assert( !tmp_code=0);
    resarr
    

let rev_comp_chr x pos = 
    match (Alphabet.complement x Alphabet.nucleotides) with 
    | Some res -> res
    | None -> failwith "cannot find complement, block_mauve.ml"

(*NOTE: rev_comp_lst just give us the complement seq, not reverse, use List.rev to do that*)
let rev_comp_lst seqlst =
    List.map (fun x ->  rev_comp_chr x (-1)
    ) seqlst

(*rev_comp_arr give us the reverse complement of seqarr*)
let rev_comp_arr seqarr = 
    let size = Array.length seqarr in
    Array.mapi (fun idx x -> rev_comp_chr seqarr.(size-idx-1) (size-idx-1)) seqarr


(* subsumed related function start *)
let is_subsumed_by j_seed i_seed mum_tbl seed2pos_tbl =
    let j_mum = get_mum_from_mumtbl j_seed mum_tbl seed2pos_tbl in
    if (j_mum.subsuming_pointer=i_seed) then true else false

(*i_mum is subsuming j_mum means, for each mach (jseqNO,jleft,jright) in j_mum,
* there is a match in i_mum (iseqNO,ileft,iright), so that iseqNO=jseqNO &&
* ileft<=jleft && iright>=jright. 
* this function works for any ext_type mums, also for more than 2 sequences*)
let is_subsuming i_mum j_mum =
    let i_pos_lst = i_mum.positions in
    let j_pos_lst = j_mum.positions in
    if (List.length i_pos_lst)>=(List.length j_pos_lst) then begin
        let sign = ref true in
        List.iter (fun j_pos ->
            let i_idx = find_index (Array.of_list i_pos_lst) j_pos 
            (fun j_pos i_pos ->
                let j_le = j_pos.left_end and j_re = j_pos.right_end in
                let i_le = i_pos.left_end and i_re = i_pos.right_end in
                if (j_pos.sequence_NO=i_pos.sequence_NO)&&
                (j_le>=i_le)&&(j_re<=i_re) then 0
                else (-1) )
            in
            if (i_idx<0) then sign := false
        ) j_pos_lst;
        !sign
    end
    else false 
(* subsumed related function end *)

(***************** extendable related functions start ************************)
(*extendable= 0,1,2,3:
0. this seed is extendable (from). 
1 ~ 3 : not extendable (from),but might be extend to.
1. during the scan_seqlst, it means this mum shows up more then
once in a sequence, also shows up in every sequence,we should not 
extend the seed from this mum(but we can extend from other mum to this one),
2 .during tke resolve_overlap, it means this mum does not
show up in one or more sequence. we don't extend to this one. 
3. during add to/remove from position2seed table, only one seed is kept as
extendable, others are marked as unextendable. but if we remove the
extendable one for some reason, we can upgrade one of this kind to be
extendable. 
*)
let mark_unextendable_mum_tbl seedNO (priority_seedNO:int option) 
(subsuming_p : int option) ext_type seed2pos_tbl mum_tbl = 
    let debug = false in
    if debug then 
        Printf.printf "mark seedNO#%d as UNextendable=%d in mum_tbl\n%!"  seedNO ext_type;
    let oldmum = get_mum_from_mumtbl seedNO mum_tbl seed2pos_tbl in
    let new_prilst =
        match priority_seedNO with
        | Some pri_seedNO -> pri_seedNO::oldmum.priority_lst
        | _ -> oldmum.priority_lst
    in
    let subsumed_by = 
        match subsuming_p with
        | Some x -> x
        | None -> oldmum.subsuming_pointer
    in
    let sign_newmum = update_mum_to_mumtbl None 
    {oldmum with subsuming_pointer = subsumed_by; 
    extendable = ext_type; priority_lst=new_prilst} mum_tbl false in
    assert(sign_newmum)

let mark_extendable_mum_tbl seedNO seed2pos_tbl mum_tbl = 
    let debug = false in
    if debug then Printf.printf "mark seedNO#%d as extendable in mum_tbl\n%!" seedNO;
    let oldmum = get_mum_from_mumtbl seedNO mum_tbl seed2pos_tbl in
    let sign_newmum = update_mum_to_mumtbl None 
    {oldmum with extendable = 0; subsuming_pointer=(-1); priority_lst=[]} mum_tbl false in
    if debug then print_mum oldmum false true;
    assert(sign_newmum)


(** [update_priority_lst] add seedNO_to_remove to priority_lst, and remove
* seedNO_to_add to priority_lst. return if priority_lst is empty.
* same seedNO in priority_lst can show up more than one time.*)
let update_priority_lst in_mum (seedNO_to_remove:int option) (seedNO_to_add:int
option) mum_tbl seed2pos_tbl=
    let debug = false in
    if debug then begin
        Printf.printf "update prilst on seed#%d," in_mum.seedNO;
        print_mum in_mum false true;
    end;
    let old_lst = in_mum.priority_lst in
    let lst_after_remove = 
    match seedNO_to_remove with
    | Some seedNO->
        let debug = 
            if (List.mem seedNO old_lst)=false then true
            else debug
        in
        if debug then print_mum in_mum false true;
        if debug then Printf.printf "remove seed#%d," seedNO;
        assert (List.mem seedNO old_lst); (*or these must be something wrong*)
        let find1 = ref false in
        List.filter (fun x -> 
            let res = (x<>seedNO)||(x=seedNO && !find1=true) in
            if x=seedNO then find1 := true;
            res
            ) old_lst 
    | None -> old_lst
    in
    let reslst = 
    match seedNO_to_add with
    | Some seedNO-> if debug then Printf.printf "add seed#%d %!" seedNO; 
                 seedNO::lst_after_remove
    | None -> lst_after_remove
    in
    if debug then Printf.printf "prilst len = %d\n%!" (List.length reslst);
    let sign_newmum = 
        update_mum_to_mumtbl None {in_mum with priority_lst=reslst} mum_tbl false
    in
    assert(sign_newmum);
    List.length reslst

let upgrade_extendable_pos2seed_tbl seedNO2upgrade_lst pos2seed_tbl_left
pos2seed_tbl_right seed2pos_tbl mum_tbl =
    let debug = false in
    let add_upgrade_seedNO_2_prilst seedNO_to_upgrade seqNO pos pos2seed_tbl =
        let old_record_k = try (Hashtbl.find pos2seed_tbl (seqNO,pos)) with
        | Not_found -> failwith ("upgrade_extendable_pos2seed_tbl, not found") in
        List.iter (fun (k_seedNO,_,_) ->
            if k_seedNO<>seedNO_to_upgrade then begin
                if debug then Printf.printf "add mum#%d to mum#%d's prilst,%!" 
                seedNO_to_upgrade k_seedNO;
                let k_mum = get_mum_from_mumtbl k_seedNO mum_tbl seed2pos_tbl in
                let len_prilst = 
                update_priority_lst k_mum None (Some seedNO_to_upgrade) mum_tbl 
                seed2pos_tbl in
                if debug then Printf.printf "len prilst<-%d\n%!" len_prilst;
            end;
        )old_record_k;
    in
    List.iter (fun seedNO ->
        let mum2upgrade = get_mum_from_mumtbl seedNO mum_tbl seed2pos_tbl in
        let prilst = mum2upgrade.priority_lst 
        and ext = mum2upgrade.extendable in
        if (prilst=[])&&(ext=3) then begin
            let poslst = mum2upgrade.positions in
            List.iter (fun pos ->
            let left,right = pos.left_end,pos.right_end in
            let seqNO = pos.sequence_NO in
            add_upgrade_seedNO_2_prilst seedNO seqNO left pos2seed_tbl_left;
            add_upgrade_seedNO_2_prilst seedNO seqNO right pos2seed_tbl_right;
            ) poslst;
        end;
    ) seedNO2upgrade_lst


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
    List.iter (fun (seqNO,pos,dir) ->
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
    List.iter (fun (seqNO,pos,dir) ->
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

(* no longer need this
* get extendable record out of in_lst, in_lst is a recordlist of a position
* (seqNO,pos), from a pos2seed_tbl 
let get_extendable_recordlst in_lst mum_tbl seed2pos_tbl = 
    List.filter (fun (seedNO,weight,ori) ->
        let mum = get_mum_from_mumtbl seedNO mum_tbl seed2pos_tbl in
        (mum.extendable=0)
    ) in_lst
*)

(* [get_extendable_record] return the extendable mum record on the input
* position (pos=(seqNO,position)), when a position has no qualified mum, 
* we return (-1,-1,-1) *)
let get_extendable_record pos2seed_tbl pos seed2pos_tbl mum_tbl =
    let debug = false in
    let recordlst = Hashtbl.find pos2seed_tbl pos in
    let recordlst = List.filter (fun record ->
        let seedNO,weight,ori = record in
        let mumi = get_mum_from_mumtbl seedNO mum_tbl seed2pos_tbl in
        if (mumi.subsuming_pointer=(-1))&&(mumi.extendable=0) then true
        else false
    ) recordlst in
    if debug&&(List.length recordlst)>1 then
        List.iter (fun (seedNO,weight,ori) ->
            let seqNO,position = pos in
            Printf.printf "seqNO=%d,pos=%d,seedNO=%d,weight=%d,ori=%d\n%!" 
            seqNO position seedNO weight ori;
    ) recordlst;
    assert( (List.length recordlst)<=1 );
    if (List.length recordlst)=0 then ((-1),(-1),(-1))
    else List.hd recordlst


(***************** extendable related functions end ************************)

(******************* seed removing function start  ************************)
(*[remove_from_pos2seed_tbl] remove mum record from position (seqNO,old_pos),
* mum_i is the mum we want to remove, mum_tbl&seed2pos_tbl is for other mums 
* on the same position,so mum_i could be removed from mum_tbl before this 
* function. return a list of mum to upgrade.
* if we remove the ext-mum from this position, and there are
* other mums on current pos, we can upgrade one of them to be extendable
* after this function by [upgrade_extendable_pos2seed_tbl]. *)
let remove_from_pos2seed_tbl pos2seed_tbl seqNO old_pos mum_i mum_tbl seed2pos_tbl= 
    let res_lst = ref [] in
    let seedNO = mum_i.seedNO in
    let debug = false in
    if debug then Printf.printf "remove seed#%d from pos2seedtbl (%d,%d) \n%!"
    seedNO seqNO old_pos;
    let old_record = try (Hashtbl.find pos2seed_tbl (seqNO,old_pos)) with
    | Not_found -> failwith ("remove pos2seedtbl, not found 1") in
    if debug then print_position2seedtbl pos2seed_tbl (Some (seqNO,old_pos));
    let new_record = List.filter (fun (x,_,_) -> x<>seedNO) old_record in
    if (mum_i.extendable=0)&&(List.length new_record)>0 then begin
        if debug then 
            Printf.printf "seed#%d is the only ext-mum on this pos\n%!" seedNO;
        let seedNO_to_upgrade,_ = List.fold_left (fun (best_seedNO,best_w)
        (j_seedNO,j_weight,_) ->
            let j_mum = get_mum_from_mumtbl j_seedNO mum_tbl seed2pos_tbl in
            let len_prilst = 
                update_priority_lst j_mum (Some seedNO) None mum_tbl seed2pos_tbl
            in
            if debug then Printf.printf "seed#%d,len_prilst=%d;" j_seedNO len_prilst;
            if (j_weight>best_w)&&(len_prilst=0)&&(j_mum.extendable=3) then
                j_seedNO,j_weight
            else 
                best_seedNO,best_w
        ) ((-1),0) new_record
        in
        if (seedNO_to_upgrade != (-1)) then
            res_lst := seedNO_to_upgrade :: !res_lst; 
    end
    else begin 
        if debug then Printf.printf "seed#%d is not extendable, or is the only \
        seed on this position,just remove it\n%!" seedNO;
    end;
    if (List.length new_record)>0 then
    Hashtbl.replace pos2seed_tbl (seqNO,old_pos) new_record
    else
        Hashtbl.remove pos2seed_tbl (seqNO,old_pos);
    if debug then Printf.printf "END of remove seed from pos2seed tbl\n%!";
    !res_lst
    

(* remove_seed2 remove j_seedNO from mum_tbl/seed2pos/pos2seed completely.*)
let remove_seed2 mumj pos2seed_tbl_left pos2seed_tbl_right mum_tbl seed2pos_tbl 
remove_from_mumtbl = 
    let seedNO2upgrade_lst = ref [] in
    let j_seedNO = mumj.seedNO in
    let debug = false in
    if debug then Printf.printf "remove mum#%d \n%!" j_seedNO;
    let j_positions = mumj.positions in
    List.iter (fun record ->
        let seqNO = record.sequence_NO 
        and left = record.left_end 
        and right = record.right_end in
        let reslst1 = remove_from_pos2seed_tbl pos2seed_tbl_left seqNO left 
        mumj mum_tbl seed2pos_tbl in
        let reslst2 = remove_from_pos2seed_tbl pos2seed_tbl_right seqNO right 
        mumj mum_tbl seed2pos_tbl in
        if reslst1<>[] then seedNO2upgrade_lst := !seedNO2upgrade_lst@reslst1;
        if reslst2<>[] then seedNO2upgrade_lst := !seedNO2upgrade_lst@reslst2;
    ) j_positions ;
    if !seedNO2upgrade_lst<>[] then
    upgrade_extendable_pos2seed_tbl !seedNO2upgrade_lst pos2seed_tbl_left
pos2seed_tbl_right seed2pos_tbl mum_tbl;
    if remove_from_mumtbl then remove_mum_from_mumtbl mumj mum_tbl;
    Hashtbl.remove seed2pos_tbl j_seedNO
(******************* seed removing function end  ************************)

(*************** functions for pos2seed table start ************************)
(* [add_to_pos2seed_tbl] add new positions of mum, might modify mum's 
* extendable sign and update it to mum_tbl.
*  Note1: pos2seed_tbl only allows at most one extendable seed start/end 
 * from on position, when there are two extendable seed, break the tie 
*  Note2: remember to reload mum after this function, for its "extendable" sign
*  could be changed. 
*  Note3: **REQUIRE** seed2pos_tbl&mum_tbl must be updated before pos2seed_tbl *)
let add_to_pos2seed_tbl main_left pos2seed_tbl_left pos2seed_tbl_right  
seqNO pos seedNO seedweight orientation ext_sign new_multi seed2pos_tbl 
mum_tbl = 
    let debug = false in
    if debug then Printf.printf "add seed#%d(weight=%d) to pos \
    (%d,%d),main_left=%b\n%!" seedNO seedweight seqNO pos main_left;
    let main_tbl = 
        if main_left then pos2seed_tbl_left  else pos2seed_tbl_right in
    (*call add_to_table to add extendable seed, call add_to_table2 to add all
    * other seed. make sure the extendable record is the head of recordlst*)
    let add_to_table key newrecord tbl =
        let old_record = Hashtbl.find tbl key in
        let new_record = newrecord::old_record in
        Hashtbl.replace tbl key new_record;
    in
    let add_to_table2 key newrecord tbl =
        let old_record = Hashtbl.find tbl key in
        let new_record = old_record@[newrecord] in
        Hashtbl.replace tbl key new_record;
    in
    if (Hashtbl.mem main_tbl (seqNO,pos) ) then begin
        let new_mum = get_mum_from_mumtbl seedNO mum_tbl seed2pos_tbl in
        let old_record = Hashtbl.find main_tbl (seqNO,pos) in
        if debug then print_position2seedtbl main_tbl (Some (seqNO,pos));
        (*the extendable record is the head of recordlst*)
        let h_seedNO,h_weight,h_ori = List.hd old_record in
        let h_mum = get_mum_from_mumtbl h_seedNO mum_tbl seed2pos_tbl in
        if (ext_sign!=0) then begin (*non-ext mum, just add it*)
            if debug then 
                Printf.printf "non-ext(ext=%d) mum, just add it.\n%!" ext_sign;
            let seedNO_to_add = 
                if h_mum.extendable<>0 then None
                else (Some h_seedNO)
            in
            let _ = update_priority_lst new_mum None seedNO_to_add mum_tbl
            seed2pos_tbl in
            add_to_table2 (seqNO,pos) (seedNO,seedweight,orientation) main_tbl;
        end
        else begin (*extendable seed,compare with old ext-record if need*)
            if h_mum.extendable<>0(*ext_record_len=0*) then begin
                if debug then Printf.printf "no ext-record before,become first\n%!";
                List.iter (fun (oldseedNO,_,_) ->
                    let oldmum = get_mum_from_mumtbl oldseedNO mum_tbl seed2pos_tbl in
                    let pri_len = 
                        update_priority_lst oldmum None (Some seedNO) mum_tbl seed2pos_tbl
                    in
                    if debug then Printf.printf "len of prilst:%d\n%!" pri_len;
                ) old_record; 
                add_to_table (seqNO,pos) (seedNO,seedweight,orientation) main_tbl;
            end
            else begin
                (*let (old_ext_seedNO,old_ext_weight,old_ext_ori) = List.hd
                * ext_old_record in*)
                let old_ext_seedNO,old_ext_weight,old_ext_ori = 
                    h_seedNO,h_weight,h_ori in
                if debug then Printf.printf "(%d,%d,%d) is the ext record on this pos\n%!"
                    old_ext_seedNO old_ext_weight old_ext_ori;
                let old_ext_mum = h_mum in
                    (*get_mum_from_mumtbl old_ext_seedNO mum_tbl seed2pos_tbl in*)
                let old_multi = List.length old_ext_mum.positions in
                if (new_multi>old_multi)||
                ((new_multi=old_multi)&&(seedweight>old_ext_weight))
                then begin (*newseed become the ext one on this pos*)
                    if debug then Printf.printf "newseed become the ext record \n%!";
                    let sign_subsume = is_subsuming new_mum old_ext_mum in
                    let subsumed_by = 
                        if sign_subsume then (Some seedNO)
                        else None in
                    mark_unextendable_mum_tbl old_ext_seedNO (Some seedNO) 
                    subsumed_by 3 seed2pos_tbl mum_tbl;
                    List.iter (fun (oldseedNO,_,_) ->
                        if oldseedNO<>old_ext_seedNO then begin
                        let oldmum = get_mum_from_mumtbl oldseedNO mum_tbl seed2pos_tbl in
                        let pri_len = update_priority_lst oldmum (Some old_ext_seedNO)
                        (Some seedNO) mum_tbl seed2pos_tbl in
                        if debug then Printf.printf "len of prilst:%d\n%!" pri_len;
                        end;
                    ) old_record;
                    add_to_table (seqNO,pos) (seedNO,seedweight,orientation) main_tbl;
                end
                else begin (*add newseed to postbl, mark it as untendable*)
                    if debug then 
                    Printf.printf "seed#%d is add to as untendatbl,%!" seedNO;
                    let sign_subsume = is_subsuming old_ext_mum new_mum in
                    let subsumed_by =
                        if sign_subsume then (Some old_ext_mum.seedNO) 
                        else None
                    in
                    mark_unextendable_mum_tbl seedNO (Some old_ext_seedNO)
                    subsumed_by 3 seed2pos_tbl mum_tbl;
                    add_to_table2 (seqNO,pos) (seedNO,seedweight,orientation) main_tbl;
                end;
            end;(*end of weather the new mum is the first ext record on this pos*)
        end (*end of if new mum is extendable*)
    end
    else begin
        if debug then Printf.printf "first record on this position\n%!";
        Hashtbl.add main_tbl (seqNO,pos) [(seedNO,seedweight,orientation)];
    end;
    if debug then begin
        print_position2seedtbl main_tbl (Some (seqNO,pos));
        Printf.printf "END of add 2 pos2seed tbl\n%!";
    end

(*[modify_record_in_pos2seed_tbl] modify weight and/or orientation of seedNO in
* pos2seed_tbl, on position (seqNO,pos), don't care about ext-sign *)    
let modify_record_in_pos2seed_tbl pos2seed_tbl seqNO pos seedNO weight ori =
    let debug = false in
    let old_record = try (Hashtbl.find pos2seed_tbl (seqNO,pos)) with
    | Not_found -> failwith ("remove pos2seedtbl, not found 1") in
    let h_seedNO,h_w,h_ori = List.hd old_record in
    let unchanged = List.filter (fun (x,_,_) -> x<>seedNO) old_record in
    if debug then Printf.printf "modify pos2seedtbl for seed#%d at pos(%d,%d) \
    w=%d,ori=%d\n%!"  seedNO seqNO pos weight ori;
    let new_record = (*make sure ext-mum is the head of recordlst*)
        if h_seedNO=seedNO then (seedNO,weight,ori)::unchanged 
        else unchanged @ [(seedNO,weight,ori)]
    in
    Hashtbl.replace pos2seed_tbl (seqNO,pos) new_record
(****************** functions for pos2seed table end ************************)

(***************** function for seed2position table start *******************)
let add_to_seed2pos_tbl seed2pos_tbl seqNO pos seedNO seedweight orientation =
    if( Hashtbl.mem seed2pos_tbl seedNO ) then begin
        let oldposlst,seedseq,_ = Hashtbl.find seed2pos_tbl seedNO in
        Hashtbl.replace seed2pos_tbl seedNO        
        (((seqNO,pos,orientation)::oldposlst),seedseq,seedweight);
    end
    else
        Hashtbl.add seed2pos_tbl seedNO ([(seqNO,pos,orientation)],[],seedweight)
                    
let remove_from_seed2pos_tbl seed2pos_tbl seqNO pos seedNO = 
    let oldposlst,seedseq,seedweight = Hashtbl.find seed2pos_tbl seedNO in
    let newposlst = List.filter (fun item ->
        let (old_seqNO,old_pos,old_dir) = item in
        if (old_seqNO=seqNO)&&(old_pos=pos) then false else true
    ) oldposlst in
    Hashtbl.replace seed2pos_tbl seedNO (newposlst,seedseq,seedweight)

(*[update_seed2pos_tbl update mumseq, new_pos_set(optional) for mum#.seedNO*)
let update_seed2pos_tbl seedNO new_pos_set new_mumseq seed2pos_tbl =
    let oldposlst,seedseq,seedweight = Hashtbl.find seed2pos_tbl seedNO in
    let newposlst = 
    match new_pos_set with
        | Some (seqNO,pos,newpos) ->
                List.map (fun (oldseqNO,oldpos,olddir) ->
                if (oldseqNO=seqNO)&&(oldpos=pos) 
                then (seqNO,newpos,olddir)
                else (oldseqNO,oldpos,olddir)
                ) oldposlst 
        | None -> oldposlst
    in
    Hashtbl.replace seed2pos_tbl seedNO 
    (newposlst,new_mumseq,(Array.length new_mumseq))
(*************** function for seed2position table end *******************)

let intlst2int inlst =
    List.fold_left (fun acc x -> acc*10+x ) 0 inlst

let radix_sort inarr =
    let max_size = 15 in
    let debug = false in
    let get_idx x = x-1  in
    let size_digits = 
        let _,_,_,first_seq,_ = inarr.(0) in
        Array.length first_seq in
    if debug then Printf.printf "digits size is %d\n%!" size_digits;
    let sort_by_digit inarr digit =
       (*List.sort (fun (_,_,_,lst1) (_,_,_,lst2) ->
            compare (List.nth lst1 digit) (List.nth lst2 digit)
        )inlst*)
       if debug then Printf.printf "sort by digit %d\n%!" digit;
        let count_arr = Array.make max_size 0 in
        Array.iteri (fun i (seqNO,pos,dir,subseq,_) ->
            let x = subseq.(digit) in
            let idx = get_idx x in
            assert((idx>=0)&&(idx<15)); (*chr not in [1,2,4,8]*)
            count_arr.(idx) <- count_arr.(idx)+1; 
            inarr.(i)<-(seqNO,pos,dir,subseq,idx);
        ) inarr;
        let index = ref 0 in
        let full_count_arr = Array.make max_size (0,0,0) in
        let _ = Array.fold_left ( fun pre_count count ->
            full_count_arr.(!index) <- (count,0,pre_count); 
            index := !index +1;
            count+pre_count
        ) 0 count_arr in
        let sorted_inarr = Array.make (Array.length inarr) (0,0,0,[||],0) in
        Array.iteri (fun i (seqNO,pos,dir,arr,idx) ->
            let total_num,hit_num,bigger_than = full_count_arr.(idx) in
            assert(hit_num<total_num);
            let index2 = bigger_than + hit_num in
            sorted_inarr.(index2)<-(seqNO,pos,dir,arr,idx);
            full_count_arr.(idx) <- (total_num,hit_num+1,bigger_than);
        )inarr;
        sorted_inarr
    in
    let digit_arr = Array.init size_digits (fun x->size_digits-x-1) in
    Array.fold_left (fun current_arrarr digit ->
        sort_by_digit current_arrarr digit
    ) inarr digit_arr 
    (*
    Array.fast_sort (fun x y ->
                let (_,_,_,lstx) = x and (_,_,_,lsty)=y in
                compare lstx lsty ) inarr;
    inarr*)
    
    

let extend_seq ext_dir init_matcharr seedweight inseqarr =
    let debug = false in
    if debug then Printf.printf "extend seq in dir=%d,seedW=%d\n%!" ext_dir
    seedweight;
    let new_seedweight = ref seedweight in
    let posarr = init_matcharr in
    let stillmatch = ref true in
    while !stillmatch do
        (*we modify posarr when we have a possible extension between current
        * match and previous match, but this extension might not work for next
        * match, in that case we need to undo the change, 
        * "posarr_copy" is here for the way back*)
        let posarr_copy = Array.init (Array.length posarr)  (fun i -> posarr.(i) ) in
        let idx = ref (-1) in
        let last_chr,last_dir = 
            Array.fold_left (fun (pre_chr,pre_dir) (seqNO,pos,dir) ->
            idx := !idx +1;
            if debug then Printf.printf "check (%d,%d,%d) with (%d,%d)\n%!" 
            seqNO pos dir pre_chr pre_dir;
            if pre_chr<0 then (pre_chr,pre_dir) (*mismatch already happened*)
            else begin
                let inseq = inseqarr.(seqNO) in
                let inseqlen = Array.length inseq in
                let ext2right = if dir*ext_dir=1 then true else false in
                let nexpos = 
                    if ext2right then !new_seedweight+pos
                    else pos-1
                in
                if debug then Printf.printf "nextpos=%d,%!" nexpos;
                if nexpos>=0&&nexpos<inseqlen then begin
                    let next_chr = 
                        if pre_chr=0||dir=pre_dir then inseq.(nexpos)
                        else rev_comp_chr inseq.(nexpos) nexpos
                    in
                    if debug then Printf.printf "chr=%d/%d\n%!" next_chr pre_chr;
                    if pre_chr=0(*the first item*)||next_chr=pre_chr then 
                        begin
                            let _,oldpos,_ = posarr.(!idx) in
                            let newpos =
                                if ext2right then oldpos
                                else oldpos-1
                            in
                            posarr.(!idx) <- (seqNO,newpos,dir);
                            (next_chr,dir)
                        end
                    else (-1,-1)
                end
                else (-1,-1)
            end
            ) (0,0) posarr in
        if last_chr<0 then begin (*undo extension, get out of while loop*)
            Array.iteri (fun i record_copy ->
                posarr.(i) <- record_copy ) posarr_copy;
            stillmatch := false
        end
        else new_seedweight := !new_seedweight + 1;
    done;
    if debug then
        Printf.printf "end of extend_seq, new seedweight = %d\n%!" !new_seedweight;
    posarr,
    !new_seedweight


let extend_seq_in_both_dir matcharr seedweight inseqarr = 
    let extend_matcharr,new_seedweight = 
        extend_seq 1 matcharr seedweight inseqarr in
    let res_matcharr,res_seedweight = 
        extend_seq (-1) extend_matcharr new_seedweight inseqarr in
    Array.sort (fun (seqNOx,_,_) (seqNOy,_,_) -> compare seqNOx seqNOy) res_matcharr;
    res_matcharr,res_seedweight

(*[add_seed_to_tbl] add match we found during [find_SML] to
* mum_tbl/seed2pos/pos2seed tbl. expend the match to both direction before
* adding it. if the result of extension is already in those tbls, skip this one*)
let add_seed_to_tbl init_accarr init_seedweight inseqarr input_seqlst_size
position2seed_tbl_left position2seed_tbl_right seed2position_tbl mum_tbl =
    let debug = false in
    if debug then begin
        Printf.printf "add seed to tbl,size of match arr=%d\n%!" 
        (Array.length init_accarr);
        (*Array.iter (fun (seqNO,pos,dir) ->
            Printf.printf "(%d,%d,%d) %!" seqNO pos dir
        ) init_accarr;
        Printf.printf "\n%!"*)
    end;
    let poslst =  Array.to_list init_accarr in
    let sign1 = match_at_most_once_in_each_seq2 poslst in
    let sign2 = match_at_least_once_in_every_seq2 poslst 
    input_seqlst_size  in
    let ext_type = 
        if (sign1)&&(sign2) then 0
        else if (sign1=false)&&sign2 then 1
        else if (sign2=false)&&sign1 then 2
        else 2
    in
    (*when ext=0,we have two match ,one in seq0, another in seq1*)
    if ext_type=0(*||ext_type=1 for now we only care about extendable ones*) then begin
        let check_old_record position2seed_tbl (seqNO,pos) ext_2dir new_seedweight
        pre_sign acclst newseedNO =
            (*a shortcut: we only allow one extendable seed start/end 
            * from one position. 
            * when adding a new seed, and there is an old seed start/or/end with
            * the same position. we could
            * 1.get rid of old seed, add the new one
            * 2.don't add new one
            * 3.get rid of old seed, and don't add the new one.
            * case 3 happens when old seed is same as new seed,just 
            * with different direction.*)
            if (Hashtbl.mem position2seed_tbl (seqNO,pos)) then begin
                    (*already have this seed match. this could happen when 
                    1.we have both seedseq and (rev seedseq) in one sequence. 
                    since we append this sequence with its own
                    revseq, the same match will show up 2n(n>1) times.
                    2.the extended subseq of other seed is the same as this
                    one, we keep the longer one.*)
                    let recordlst = Hashtbl.find position2seed_tbl (seqNO,pos) in
                    assert((List.length recordlst)=1);
                    if debug then begin
                        Printf.printf "new seedweight=%d(2dir=%d),old record:%!"
                        new_seedweight ext_2dir;
                        print_position2seedtbl position2seed_tbl (Some (seqNO,pos)); 
                    end;
                    let (oldseedNO,oldweight,olddir) = List.hd recordlst in
                    if oldweight>new_seedweight then
                        2,[]
                    else if oldweight=new_seedweight then begin
                        (*we find the same match again, if it's just with
                        * different dir, then it has ext=1*)
                        let old_poslst,_,_ = try (Hashtbl.find seed2position_tbl
                        oldseedNO) with |Not_found -> 
                            failwith "not found seed2postbl, check old record" in
                        let _,_,old_dir0 = List.hd old_poslst in
                        let _,_,old_dir1 = List.nth old_poslst 1 in
                        let old_2dir = old_dir0*old_dir1 in
                        if debug then Printf.printf "old_2dir=%d,%!" old_2dir;
                        if old_2dir<>ext_2dir then
                            let newacclst = 
                                if (List.mem oldseedNO acclst) then acclst
                                else oldseedNO::acclst in
                            3,newacclst
                        else
                            2,[]
                    end
                    else begin
                        if debug then 
                        Printf.printf "newseed=%d,remove oldseed#%d(%d<%d)\n%!" 
                        newseedNO oldseedNO oldweight new_seedweight;
                        let newacclst = 
                            if (List.mem oldseedNO acclst) then acclst
                            else oldseedNO::acclst in
                        pre_sign,newacclst;
                    end;
            end
            else pre_sign,acclst;
        in
        let ext_posarr,ext_seedweight = 
            extend_seq_in_both_dir init_accarr init_seedweight inseqarr in
        let seqNO,pos,dir0 = ext_posarr.(0) in 
        let _,_,dir1 = ext_posarr.(1) in         
        let ext_2dir = dir0*dir1 in
        assert (seqNO=0);
        let seedNO = get_a_seedNO () in
        let newsign,seed2remove_lst = 
            Array.fold_left (fun (pre_sign,acclst) (ext_seqNO,ext_pos,_) -> 
            if pre_sign<>2 then begin
                (*check left table first*)
                let sign_left,acclst_left = check_old_record position2seed_tbl_left 
                (ext_seqNO,ext_pos) ext_2dir ext_seedweight pre_sign acclst seedNO in
                if sign_left<>2 then (*then check right table*)
                 check_old_record position2seed_tbl_right 
                 (ext_seqNO,ext_pos+ext_seedweight-1) ext_2dir 
                 ext_seedweight sign_left acclst_left seedNO
                else
                    2,[]
            end
            else
                pre_sign,[]
            ;
        ) (1,[]) ext_posarr in
        if debug then Printf.printf "newsign=%d,\n%!" newsign;
        if seed2remove_lst<>[] then 
            List.iter (fun seed2remove ->
                if debug then Printf.printf "remove seed#%d\n%!" seed2remove;
                let poslst,_,seedweight = try (Hashtbl.find seed2position_tbl
                seed2remove) with |Not_found -> 
                    let _ = Printf.printf  "cannot find seed %d in seed2pos tbl"
                    seed2remove in
                    failwith "add to seed, not found in seed2pos tbl"
                in
                List.iter (fun (seqNO,pos,dir) ->
                    Hashtbl.remove position2seed_tbl_left (seqNO,pos);
                    Hashtbl.remove position2seed_tbl_right (seqNO,pos+seedweight-1);
                ) poslst;
                let mum2remove = get_mum_from_mumtbl seed2remove mum_tbl
                seed2position_tbl in
                remove_mum_from_mumtbl mum2remove mum_tbl;
                Hashtbl.remove seed2position_tbl seed2remove;
                return_a_seedNO seed2remove;
            ) seed2remove_lst;
        let ext_poslst = Array.to_list ext_posarr in
        if newsign=1 then begin
            let mumseq = 
            get_sub_seq2 inseqarr.(seqNO) pos ext_seedweight in
            Hashtbl.add seed2position_tbl seedNO (ext_poslst,mumseq,ext_seedweight);
            let key = get_mumkey_from_poslst ext_poslst in
            let newpositions = List.map (fun (seqNO,pos,dir) ->
                Hashtbl.add position2seed_tbl_left (seqNO,pos) 
                [(seedNO,ext_seedweight,dir)];
                Hashtbl.add position2seed_tbl_right (seqNO,(pos+ext_seedweight-1))
                [(seedNO,ext_seedweight,dir)];
                let left_end = pos and right_end = pos + ext_seedweight - 1 in
                {
                    sequence_NO = seqNO; 
                    left_end = left_end; 
                    right_end = right_end;
                    orientation = dir;
                } 
            ) ext_poslst in
            let mumsize = List.length ext_poslst in
            let newmum =
                {   seedNO = seedNO;
                    mumseq = mumseq;
                    positions = newpositions;
                    mumkey = key;
                    size = mumsize;
                    neighborhood_lst = [];
                    subsuming_pointer = (-1);
                    extendable = if mumsize=2 then 0 else 1;
                    priority_lst=[];
                    mumscore = 0;
                } in
            if debug then begin
                Printf.printf "add seed %d with mumkey %d,treekey = %!" seedNO key;
                print_mum newmum false true;
                Printf.printf " to mum_tbl\n%!";
            end;
            let sign_newmum = add_mum_to_mumtbl newmum mum_tbl in
            assert(sign_newmum); 
        end
        else
            return_a_seedNO seedNO;
    end


let find_SML patternarr init_seedweight  inseqarr input_seqlst_size
position2seed_tbl_left position2seed_tbl_right seed2position_tbl mum_tbl = 
    let debug = false and debug2 = false in
    if debug then 
        Printf.printf "get sequences with index and possible seed subseq\n%!";
    let inseqarr_w_idxarr = Array.mapi (fun seqNO inseq->
        let seqlen = Array.length inseq in
        if debug then Printf.printf "seq len = %d\n%!" seqlen;
        let resarr =  
            Array.mapi (fun idx x ->
            let sequenceNO,pos,dir,subseqlst =
                if ( (idx+init_seedweight-1) < seqlen) then 
                    let subseq = 
                     Array.fold_right 
                     (fun p acc ->
                         let chr=inseq.(idx+p-1) in
                         chr::acc ) patternarr [] 
                    in
                    seqNO,idx,1,subseq
                else 0,0,0,[]
            in
            if debug2&&(idx mod 100000)=0 then
                Printf.printf "reach idx=%d\n%!" idx;
            (*0,false are for function radix_sort later.*)
            sequenceNO,pos,dir,Array.of_list subseqlst,0
        ) inseq in
        Array.sub resarr 0 (seqlen-init_seedweight+1)
    ) inseqarr
    in
    if debug then 
        Printf.printf "append it with reverse complement of first sequence\n%!";
    let rev_seq_w_idx = Array.map (fun (sequenceNO,pos,dir,subseq,_) ->
        let rev_subseq = rev_comp_arr subseq in
        let rev_pos = pos in
        (sequenceNO,rev_pos,-1,rev_subseq,0)
    ) inseqarr_w_idxarr.(0) in
    if debug then Printf.printf "array append\n%!";
    let inseqarr_w_idxarr_w_rev = 
        Array.append (Array.append inseqarr_w_idxarr.(0) inseqarr_w_idxarr.(1)) rev_seq_w_idx
    in
    if debug then Printf.printf "call radix sort\n%!";
    let sorted_arr = radix_sort inseqarr_w_idxarr_w_rev in
    if debug then Printf.printf "end of radix sort\n%!";
    let pre_subseq = ref [||] and pre_sign = ref 0 in
    let sorted_arr = Array.mapi (fun idx (seqNO,pos,dir,(subseq:int array),_) ->
        if (compare !pre_subseq subseq)=0 then begin
            pre_sign := !pre_sign+1;
        end
        else begin
            pre_subseq := subseq;
            pre_sign := 0;
        end;
        (!pre_sign,seqNO,pos,dir,subseq)
    ) sorted_arr in
    if debug2 then begin 
        Printf.printf "sorted arr is:\n%!";
        Array.iter (fun (count,sequenceNO,pos,dir,subseqarr) ->
        Printf.printf "%2i,%2i,%3i,%2i,[%!" count sequenceNO pos dir;
        Utl.printIntArr subseqarr;
        Printf.printf "]\n%!";
        ) sorted_arr;
    end;
    if debug then Printf.printf "start adding seeds to tables\n%!";
    let last_acclst = Array.fold_left (fun acc item ->
        let (sign,seqNO,pos,dir,subseq) = item in
        let newacc =
            if sign>0 then  acc@[(seqNO,pos,dir)] 
            else [(seqNO,pos,dir)] in
        if sign=0&&(List.length acc)>1 then 
            add_seed_to_tbl (Array.of_list acc) init_seedweight inseqarr input_seqlst_size
            position2seed_tbl_left position2seed_tbl_right seed2position_tbl mum_tbl;
        newacc
    ) [] sorted_arr in
    if (List.length last_acclst)>1 then 
        add_seed_to_tbl (Array.of_list last_acclst) init_seedweight inseqarr input_seqlst_size
position2seed_tbl_left position2seed_tbl_right seed2position_tbl mum_tbl
        

(*try to deal with huge dataset, remember to use 64bits*)
let scan_seqlst2 (inseqarr:int array array) patternarr mum_tbl
position2seed_tbl_left position2seed_tbl_right seed2position_tbl debug=
    let input_seqarr_size = 2 in (*now we only have 2 input sequences*)
    let init_seedweight = patternarr.((Array.length patternarr)-1) in
    if debug then Printf.printf "scan_seqlst2, seedweight=%d, find_SML with inseqlst_w_idx_w_rev%!" 
    init_seedweight;
    find_SML patternarr init_seedweight inseqarr input_seqarr_size 
    position2seed_tbl_left position2seed_tbl_right seed2position_tbl mum_tbl ; 
    if debug then Printf.printf "find_SML for reverse complement done\n%!";
    let debug2 = false in
    if debug2 then begin
        Printf.printf "check seed2pos tbl :\n%!";
        print_seed2postbl seed2position_tbl;
        Printf.printf "check pos2seed tbl :\n%!";
        print_position2seedtbl position2seed_tbl_left None;
    end;
    if debug then begin
        Printf.printf "end of scan_seq2 \n%!";
    end;
    Hashtbl.length seed2position_tbl



let build_seed_and_position_tbl inseqarr init_seedlen min_seednumber 
seq2seedNO_tbl seedNO2seq_tbl pos2seed_tbl_left pos2seed_tbl_right
seed2pos_tbl mum_tbl debug
=
    if debug then Printf.printf "build_seed_and_position_tbl, init_seedlen = %d\n%! "
    init_seedlen;
    let seed_len = ref init_seedlen in
    let seedlst = 
        try (Hashtbl.find palindromic_spaced_seed_tbl init_seedlen)
        with | Not_found -> 
            let _ = Printf.printf "looking for %d\n%!" init_seedlen in
            failwith ("cannot find seed pattern1")
    in
    let seedpattern = List.hd seedlst in
    if debug then begin 
        Printf.printf "seed pattern list = %!";
        print_int_list (seedpattern); Printf.printf "\n%!";
    end;
    let patternarr = Array.of_list seedpattern in
    let _ = scan_seqlst2 inseqarr patternarr  
        mum_tbl pos2seed_tbl_left pos2seed_tbl_right seed2pos_tbl false in
    (*
    let seed_number = ref init_seed_number in
    while (!seed_number<min_seednumber)&&(!seed_len>5 ) do 
        seed_len := !seed_len - 2;
        let seed_len = if !seed_len = 17 then 15 else !seed_len in
        let seedlst = 
        try (Hashtbl.find palindromic_spaced_seed_tbl seed_len)
        with | Not_found -> failwith ("cannot find seed pattern 2")
        in
        let seedpattern = List.hd seedlst in
        Printf.printf "reduce seed weight (%d<%d), new seed pattern list = %!"
        !seed_number min_seednumber;
        print_int_list (seedpattern); Printf.printf "\n%!";
        let patternarr = Array.of_list seedpattern in
        let new_seed_number = scan_seqlst2 inseqarr patternarr 
        mum_tbl pos2seed_tbl_left pos2seed_tbl_right seed2pos_tbl true in
        seed_number := new_seed_number;
    done;*)
    if debug then begin
    (*Printf.printf "check position2seed tbl left and right\n%!";
    print_position2seedtbl pos2seed_tbl_left None;
    print_position2seedtbl pos2seed_tbl_right None;
    Printf.printf "check seed2position tbl :\n%!";
    print_seed2postbl seed2pos_tbl;
    Printf.printf "check seedNO2seq_tbl :\n%!";
    print_seedNO2seqtbl seedNO2seq_tbl; *)
    Printf.printf "end of build_seed_and_position_tbl,seed_len=%d\n%!" !seed_len;
    end;
    patternarr.((Array.length patternarr)-1)

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
    
let get_max_range in_lst = 
    let tmp = List.hd in_lst in
    tmp.right_end - tmp.left_end

(* [update_positions_by_seqNO] update positions of mum with newpos=(seqNO,pos)
  this works for adding new position to a pos_list, or update content of a
* position in the pos_list, in mum_tbl*)
(*this works only for input two sequences*)
let update_positions_by_seqNO seedNO seqNO oldmum new_pos new_mumseq mum_tbl = 
    let debug = false in
    let input_seqlst_size = 2 in
    if debug then Printf.printf "update_positions_by_seqNO,seedNO:%d\n%!" seedNO;
    let old_poslen = List.length oldmum.positions in
    let positions_unchanged = List.filter (fun record ->
        if(record.sequence_NO <> seqNO) then true
        else false
    ) oldmum.positions in 
    let new_positions = new_pos :: positions_unchanged in
    let new_mumkey = get_mumkey_from_milst new_positions in
    let new_poslen = List.length new_positions in
    let ext_sign = 
        if (new_poslen = old_poslen) then oldmum.extendable
        else
            get_extendable_type new_positions input_seqlst_size debug 
    in
    let new_mum = 
        {oldmum with positions = new_positions; 
        mumkey = new_mumkey;
        mumseq = new_mumseq;
        size = List.length new_positions; extendable = ext_sign } 
    in
    let sign_newmum = update_mum_to_mumtbl (Some oldmum) new_mum mum_tbl false in
    if debug then begin
        Printf.printf "update_positions_by_seqNO END, (sign_newmum=%b):\n%!" sign_newmum;
        if sign_newmum then print_mum new_mum false true;
    end;
    sign_newmum

let remove_positions_by_seqNO seedNO seqNO mum_tbl seed2pos_tbl = 
    let debug = false in
    let mi = get_mum_from_mumtbl seedNO mum_tbl seed2pos_tbl in 
    let old_mumkey = get_mumkey_from_milst mi.positions in
    let positions_unchanged = List.filter (fun record ->
        if(record.sequence_NO <> seqNO) then true
        else false
    ) mi.positions in 
    let new_positions = positions_unchanged in
    let new_mumkey = get_mumkey_from_milst new_positions in
    let new_mi = {mi with 
    positions = new_positions; 
    mumkey = new_mumkey;
    size = List.length new_positions } in
    if debug then Printf.printf "update mum#%d,new_mumkey = %d, old_mumkey = %d.%!"
    seedNO new_mumkey old_mumkey;
    let sign_newmum = update_mum_to_mumtbl (Some mi) new_mi mum_tbl false in
                        assert(sign_newmum);
    if debug then Printf.printf "sign=%b\n%!" sign_newmum


(*our idx starts from 0, while mauve starts from 1*)
let format_output mum_tbl = 
    let res = ref [[]] in
    let format_mum mum =
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
        else ()
    in
    Hashtbl.iter (fun mumkey bt ->
        iter_b_tree bt false format_mum
    ) mum_tbl;
    !res


let sort_tbl_by_size in_tbl =
    let p_queue = ref [] in
    Hashtbl.iter (fun seedNO mum ->
        let size = mum.size in
        p_queue := (seedNO, size)::!p_queue;
    ) in_tbl;
    List.sort (fun (_,sizex) (_,sizey) ->
        compare sizey sizex ) !p_queue 


let sort_neighborhood_lst_by_dis in_lst = 
    List.sort(fun x y ->
        let (_,_,_,_,dxl,dxr) = x and (_,_,_,_,dyl,dyr) = y in
        compare (abs dxl) (abs dyl)
    ) in_lst

(*in_lst is neighborhood list of one j_seedNO to i_seedNO *)
let group_neighborhood_lst_by_dis in_lst debug =
    if debug then Printf.printf "group neighborhood lst\n%!";
    let in_lst = sort_neighborhood_lst_by_dis in_lst in
    let first_record = List.hd in_lst in
    let (_,_,first_i_ori,first_j_ori,first_distanceL,first_distanceR) = first_record in
    let res = ref [] in
    if debug then Printf.printf "start with (%d,%d,%d,%d)\n%!" first_i_ori first_j_ori
    first_distanceL first_distanceR;
    let _,_,last_disL,last_distR,last_acc = List.fold_right (fun 
        record (pre_iori,pre_jori,previous_disL,previous_disR, acc) ->
        let (_,_,i_ori,j_ori,distanceL,distanceR) = record in
        if debug then Printf.printf "work on (%d,%d,%d,%d) with pre=(%d,%d,%d,%d)\n%!" 
        i_ori j_ori distanceL distanceR pre_iori pre_jori previous_disL previous_disR;
        if (pre_iori*pre_jori=i_ori*j_ori)&&
        (pre_iori*previous_disL=i_ori*distanceL)&&
        (pre_iori*previous_disR=i_ori*distanceR) then 
            (i_ori,j_ori,distanceL,distanceR,record::acc)
        else begin
            if debug then Printf.printf "different distance,wrap up\n%!";
            res := (previous_disL,acc) :: !res;
            (i_ori,j_ori,distanceL,distanceR,[record])
        end
    ) (List.tl in_lst) 
    (first_i_ori,first_j_ori,first_distanceL,first_distanceR,[first_record]) in
    res := (last_disL,last_acc) :: !res;
    (*make sure the first record of res is the one group with smallest distance*)
    List.rev !res




(* this allows more than one seed starts/ends from one positon*)
(* we call this to expand the position to left or/and right*)
let update_position2seedtbl seqNO old_pos new_pos seed_weight ori seedNO positiontbl =
    let debug = false in
    if debug then
        Printf.printf "update pos2seedtbl on seed#%d,pos=(%d,old=%d,new=%d),\
        weight=%d\n%!" seedNO seqNO old_pos new_pos seed_weight;
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

let categorize_neighboorhood_lst in_list i_seedNO i_size mum_tbl seed2pos_tbl =
    let debug = false in
    if debug then begin 
        Printf.printf " categorize_neighboorhood_lst, check in_list : \n%!";
        print_neighborhood_list in_list;
    end;
    let overlap_lst = ref [] in
    (*we don't need to worry about j_mum.seqNO = i_mum.seqNO
    * mum shows up only onece in a sequence,*)
    let tbj = Hashtbl.create init_tb_size in
    List.iter (fun item ->
        let (seqNO,j,i_ori,j_ori,d,dr) = item in
        if (is_subsumed_by j i_seedNO mum_tbl seed2pos_tbl) then ()
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
        (* this part is for chainable list, we don't need them now.
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
let update_neighborhood seedNO mum_tbl seed2pos_tbl pos_tbl_left pos_tbl_right debug =
    let debug2 = false in
    if debug then Printf.printf "update neighborhood on seedNO# %d\n%!" seedNO;
    let i_mum = get_mum_from_mumtbl seedNO mum_tbl seed2pos_tbl in
    if debug then print_mum i_mum true true;
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
            if debug then Printf.printf "work on mi{seqNO=%d,left=%d,right=%d}\n%!"
            seqNO i_left_end i_right_end ;
            let stop = ref 0 in
            (*"next to" also count as "overlap", we add it into neighborhood, too*)
            let pos = ref (i_left_end-d-1) in
            while (!stop<>1)&&( (!pos) <= (i_right_end+d+1)) do
                if ((!pos<>i_left_end)&&(Hashtbl.mem pos_tbl_right (seqNO,!pos)))then begin
                    let record_list = Hashtbl.find pos_tbl_right (seqNO,!pos) in
                    List.iter (fun (j_seedNO,seed_weight,j_ori) ->
                        if debug2 then  
                            Printf.printf "right tbl,pos=(%d,%d),i_seedNO=%d,\
                            j_seedNO#%d\n%!" seqNO !pos seedNO j_seedNO;
                        let j_mum = 
                            get_mum_from_mumtbl j_seedNO mum_tbl seed2pos_tbl in
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
                            if debug2 then print_mum j_mum false true;
                            if (j_idx<0) then
                                Printf.printf "r-tbl,j_idx<0,j=%d,pos=%d\n%!" j_seedNO !pos;
                            assert(j_idx>=0);
                            let j_pos = List.nth j_mum.positions j_idx in
                            let j_left_end = j_pos.left_end 
                            and j_right_end = j_pos.right_end in
                            let newneighbor =
                                (seqNO,j_seedNO,i_ori,j_ori,i_left_end-j_left_end,
                                i_right_end-j_right_end)
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
                        if debug2 then  
                            Printf.printf "left tbl,pos=(%d,%d),j_seedNO#.%d\n%!" 
                            seqNO !pos j_seedNO;
                        let j_mum = 
                            get_mum_from_mumtbl j_seedNO mum_tbl seed2pos_tbl in
                        if debug2 then print_mum j_mum false true;
                        if (j_mum.subsuming_pointer<>(-1)) then ()
                        else if (j_seedNO<>seedNO) then begin 
                            let j_idx = 
                                find_index (Array.of_list j_mum.positions)
                                !pos (fun p item ->
                                    if (item.left_end=p) then 0 else (-1) )
                            in
                            if (j_idx<0) then
                                Printf.printf "j_idx<0,j=%d,pos=%d\n%!" j_seedNO !pos;
                            assert(j_idx>=0);
                            assert( (i_left_end-(!pos))<>0);
                            let j_pos = List.nth j_mum.positions j_idx in
                            let j_right_end = j_pos.right_end in
                            let newneighbor =
                                (seqNO,j_seedNO,i_ori,j_ori,i_left_end-(!pos),
                                i_right_end-j_right_end) in
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
        let sign_newmum = update_mum_to_mumtbl None new_i_mum mum_tbl false in
        assert(sign_newmum);
        if debug then begin
            Printf.printf "end of update_neighborhood check mum again\n%!";
            let imum =  get_mum_from_mumtbl seedNO mum_tbl seed2pos_tbl in
            print_mum imum true true;
        end
    end
    else () (*if mi.extendable = 0 or 3*)
  (*  end *)
 

let build_neighborhood mum_tbl pos_tbl_left pos_tbl_right seed2pos_tbl debug =
    if debug then Printf.printf "build neighborhood start \n%!";
    Hashtbl.iter (fun seedNO record ->
        update_neighborhood seedNO mum_tbl seed2pos_tbl
        pos_tbl_left pos_tbl_right false;
    ) seed2pos_tbl;
    if debug then Printf.printf "build neighborhood done\n%!"

(*make a chain chain i_seed and j_seed together. chainable_lst is a sublist of
* i_seed's neighborhood, each of them is a record of j_seed. 
* Note: number of j_seed in chainable_lst might be larger then number of
* matches of i_seed, for j_seed could show up more than once in a sequence -- in
* that case, j_seed's extendable is 1.*)
let make_a_chain i_seedNO chainable_lst mum_tbl seed2pos_tbl 
pos2seed_tbl_left pos2seed_tbl_right find_a_chain =
    let debug = false in
    (*we only chain the closest j_mum*)
    let (distance,nei_list) = List.hd chainable_lst in
    let (_,j_seedNO,_,_,_,_) = List.hd nei_list in
    let j_mum = get_mum_from_mumtbl j_seedNO mum_tbl seed2pos_tbl in
    let j_pos_lst = j_mum.positions in
    if debug then Printf.printf "make a chain start \n%!";
    if debug then print_neighborhood_list nei_list;
    let i_mum = get_mum_from_mumtbl i_seedNO mum_tbl seed2pos_tbl in
    let i_pos_lst = i_mum.positions in
    if (j_mum.subsuming_pointer = i_seedNO) then ()
    else begin
        let debug = false in
        let seedNO2upgrade_lst = ref [] in
        List.iter (fun (seqNO,j_seedNO,i_ori,j_ori,d,dr) ->
            let i_mum = get_mum_from_mumtbl i_seedNO mum_tbl seed2pos_tbl in
            let j_mum = get_mum_from_mumtbl j_seedNO mum_tbl seed2pos_tbl in
            if debug then begin
            Printf.printf "chain seed#%d and #%d together,d=%d,i_ori,j_ori=%d/%d \n%!"
            i_seedNO j_seedNO d i_ori j_ori;
            Printf.printf "before chainning:\n%!";
            print_mum i_mum false true;
            print_mum j_mum false true;
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
            let jmumseq = j_mum.mumseq and imumseq = i_mum.mumseq in 
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
            let new_imumseq = (*we only take subseq in seq0 as mumseq*)
                if (seqNO=0)&&expand_to_left then
                    Array.append (get_sub_seq2 jmumseq 0 (mi.left_end-mj.left_end)) imumseq
                else if (seqNO=0)&&expand_to_right then
                    Array.append imumseq
                    (get_sub_seq2 jmumseq (mi.right_end-mj.left_end+1)
                    (mj.right_end-mi.right_end))
                else imumseq
            in
            if debug then 
                Printf.printf "new_l_end/r_end=%d/%d\n%!" new_left_end new_right_end;
            let seed_weight = new_right_end - new_left_end+1 in
            let new_pos = 
                {mi with left_end = new_left_end; right_end = new_right_end }
            in
            if expand_to_left then
               update_seed2pos_tbl i_seedNO (Some(i_seqNO,i_left_end,new_left_end)) 
               new_imumseq seed2pos_tbl
            else  
                update_seed2pos_tbl i_seedNO None new_imumseq seed2pos_tbl;
            let sign_newmum = update_positions_by_seqNO i_seedNO i_seqNO i_mum 
            new_pos new_imumseq mum_tbl in
            if sign_newmum=false then begin
                if debug then Printf.printf "remove seed#%d from tbls\n%!" i_seedNO;
                find_a_chain := 0;
                (*no need to remove old i_mum, it's done in update_mum_to_mumtbl
                * called by update_positions_by_seqNO *)
                remove_seed2 i_mum pos2seed_tbl_left pos2seed_tbl_right mum_tbl
                seed2pos_tbl false;
            end
            else begin
            if expand_to_left then begin
                (*do we need to modify seedNO2seq_tbl, seq2seedNO_tbl?*)
                let res_lst = remove_from_pos2seed_tbl pos2seed_tbl_left i_seqNO
                i_left_end i_mum mum_tbl seed2pos_tbl in
                if res_lst <> [] then
                    seedNO2upgrade_lst := !seedNO2upgrade_lst@res_lst;
                add_to_pos2seed_tbl true pos2seed_tbl_left pos2seed_tbl_right i_seqNO 
                mj.left_end i_seedNO seed_weight mi.orientation i_mum.extendable
                i_mum.size seed2pos_tbl mum_tbl;
                modify_record_in_pos2seed_tbl pos2seed_tbl_right i_seqNO
                new_right_end i_seedNO seed_weight mi.orientation;
            end ;
            if expand_to_right then begin
                let res_lst = remove_from_pos2seed_tbl pos2seed_tbl_right i_seqNO
                i_right_end i_mum mum_tbl seed2pos_tbl in
                if res_lst <> [] then
                    seedNO2upgrade_lst := !seedNO2upgrade_lst@res_lst;
                add_to_pos2seed_tbl false pos2seed_tbl_left pos2seed_tbl_right i_seqNO 
                mj.right_end i_seedNO seed_weight mi.orientation i_mum.extendable 
                i_mum.size seed2pos_tbl mum_tbl;
                modify_record_in_pos2seed_tbl pos2seed_tbl_left i_seqNO
                new_left_end i_seedNO seed_weight mi.orientation;
            end;      
            end
        )nei_list;
        if !seedNO2upgrade_lst <> [] then 
            upgrade_extendable_pos2seed_tbl !seedNO2upgrade_lst pos2seed_tbl_left
pos2seed_tbl_right seed2pos_tbl mum_tbl;
    end;
    (*set j's subsuming_pointer to i*)
    let j_mum = get_mum_from_mumtbl j_seedNO mum_tbl seed2pos_tbl in
    let sign_newmum = update_mum_to_mumtbl None {j_mum with subsuming_pointer = i_seedNO } mum_tbl false in
    assert(sign_newmum);
    if debug then begin
        Printf.printf "after chain: \n%!";
        let i_mum = get_mum_from_mumtbl i_seedNO mum_tbl seed2pos_tbl in
        let j_mum = get_mum_from_mumtbl j_seedNO mum_tbl seed2pos_tbl in
        print_mum i_mum false true;
        print_mum j_mum false true;
    end
(* chainable_lst *)


let extend_seeds mum_tbl seed2pos_tbl pos2seed_tbl_left pos2seed_tbl_right =
    let debug = false and debug2 = false in
    if debug then Printf.printf "\n ===========   extend seeds ============= \n%!";
    Hashtbl.iter (fun i_seedNO _ ->
        let debug = false in
        if debug then Printf.printf "+++++++ extend seeds on seed#%d \n%!" i_seedNO;
        let i_mum = get_mum_from_mumtbl i_seedNO mum_tbl seed2pos_tbl in
        if ( i_mum.subsuming_pointer <> (-1)) then begin 
            if debug then 
            Printf.printf " is subsumed by seed#%d\n%!" i_mum.subsuming_pointer; 
        end
        else if ( i_mum.extendable != 0 ) then begin
            if debug then Printf.printf "this one is not extendable from\n%!";
        end
        else begin
            let find_a_chain = ref 1 in
            while ( !find_a_chain = 1 ) do
                if debug then begin
                Printf.printf "=== begin of while ===> check pos2seed tbl \n%!";
                if debug2 then print_position2seedtbl pos2seed_tbl_left None;
                if debug2 then print_position2seedtbl pos2seed_tbl_right None;
                end; 
                update_neighborhood i_seedNO mum_tbl seed2pos_tbl
                pos2seed_tbl_left pos2seed_tbl_right false;
                let i_mum = get_mum_from_mumtbl i_seedNO mum_tbl seed2pos_tbl in
                if debug then begin
                    Printf.printf "extend seed#%d\n%!" i_seedNO;
                    print_mum i_mum false true;
                end;
                    if debug then Printf.printf "categorize neighborhood :\n%!";
                    let nei_list = i_mum.neighborhood_lst in
                    let overlap_list = (*nei list for each j_seedNO*)
                        categorize_neighboorhood_lst nei_list i_seedNO
                        i_mum.size mum_tbl seed2pos_tbl in
                    find_a_chain := 0;
                    if (List.length overlap_list)>0 then begin
                        find_a_chain := 1;
                        (*we only care about the closest overlap-lst to i_mum*)
                        if debug then 
                            List.iter ( fun (dis,nei_lst) -> 
                            Printf.printf "closest nei list with distance = %d\n%!" dis;
                            print_neighborhood_list nei_lst;
                            ) overlap_list;
                        make_a_chain i_seedNO overlap_list mum_tbl seed2pos_tbl
                        pos2seed_tbl_left pos2seed_tbl_right find_a_chain;
                        if debug then 
                            if (Hashtbl.mem pos2seed_tbl_right (0,13497)) then
                                print_position2seedtbl pos2seed_tbl_right
                                (Some (0,13497));
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
    ) seed2pos_tbl

    
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

    
(*this is for comparison of our result and mauve's*)
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




let build_local_mums2 mum_tbl seed2pos_tbl pos2seed_tbl_left pos2seed_tbl_right
 debug debug_neighborhood =
    if debug then Printf.printf "build local mums2 start,seed number=%d\n%!"
     (Hashtbl.length seed2pos_tbl);
    build_neighborhood mum_tbl pos2seed_tbl_left
    pos2seed_tbl_right seed2pos_tbl debug_neighborhood;
    if debug then begin
    Printf.printf "check mum tbl after build neighborhood, size=%d \
    (un-extendable included)\n%!" (Hashtbl.length mum_tbl);
    end; 
    extend_seeds mum_tbl seed2pos_tbl  pos2seed_tbl_left pos2seed_tbl_right;
    if debug then Printf.printf "build_local_mums2 finished\n%!"

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
    let ourfunction new_mum = 
        if debug then Printf.printf "mum before transpose: %!";
        if debug then print_mum new_mum false true;
        if new_mum.extendable<>2 then begin
            let sign = ref true in (*true means we can add the new mum to mum_tbl*)
            let ori_positions = List.map (fun mi ->
                let ori_l,ori_r = to_ori_position mi lcb_range_lstlst in
                if (ori_l=(-2)) then sign:=false;
                {mi with left_end = ori_l; right_end = ori_r}
            ) new_mum.positions in
            if (!sign = true) then begin
                    let seedNO = get_a_seedNO () in
                    let mumkey = get_mumkey_from_milst ori_positions in
                    let mum2add = { new_mum with 
                                    seedNO = seedNO;
                                    positions = ori_positions; 
                                    mumkey = mumkey;
                                    subsuming_pointer = (-1);
                                    extendable = 0;(*we will change extendable
                                    sign when we add this mum to pos2seed_tbl*)
                                    priority_lst = [];
                                    neighborhood_lst = [] } in
                    if debug then Printf.printf "mum after transpose :%!";
                    if debug then print_mum mum2add false true;
                    let sign_newmum = add_mum_to_mumtbl mum2add mum_tbl in
                    if debug then 
                        Printf.printf "sign_newmum=%b\n%!" sign_newmum;
                    if sign_newmum then begin
                        let poslst = List.map (fun mi ->
                            (mi.sequence_NO,mi.left_end,mi.orientation) ) ori_positions 
                        in
                        Hashtbl.add seed2pos_tbl seedNO (poslst,mum2add.mumseq,
                        (Array.length mum2add.mumseq));
                        List.iter ( fun mi ->
                        let seqNO = mi.sequence_NO and ori = mi.orientation in
                        let left = mi.left_end and right = mi.right_end in
                        let weight = right-left+1 and ext_sign = 0 in
                        let mumsize = new_mum.size in
                        add_to_pos2seed_tbl true pos2seed_tbl_left
                        pos2seed_tbl_right seqNO left seedNO
                        weight ori ext_sign mumsize seed2pos_tbl mum_tbl;
                        add_to_pos2seed_tbl false pos2seed_tbl_left 
                        pos2seed_tbl_right seqNO right seedNO
                        weight ori ext_sign mumsize seed2pos_tbl mum_tbl;
                        )ori_positions;
                    end;
            end;
        end;
        (*we return the seedNO of original newmum anyway*)
        return_a_seedNO new_mum.seedNO;
    in
    Hashtbl.iter (fun mumkey bt ->
        iter_b_tree bt false ourfunction
    ) new_mum_tbl;
    Hashtbl.iter (fun seedNO record ->
        let debug_nei = false in
        update_neighborhood seedNO mum_tbl seed2pos_tbl 
        pos2seed_tbl_left pos2seed_tbl_right debug_nei;
    ) seed2pos_tbl;
    if debug then Printf.printf "end of transpose back\n%!"
    


    
(*[get_mum_lst_for_each_seq] returns the lcb_lst_lst we need to build lcb blocks 
* to do: rewrite this with seed2pos tbl *)
let get_mum_lst_for_each_seq mum_tbl seed2pos_tbl pos2seed_tbl seqlst_size seq_size_lst 
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
                    get_extendable_record pos2seed_tbl (i,j)
                    seed2pos_tbl mum_tbl in
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
    let seedNOlst1 = List.hd seedNOlstlst 
    and seedNOlst2 = List.nth seedNOlstlst 1 in
    List.iter (fun seedNO ->
        let negseedNO = (-1)*seedNO in
        if (List.mem seedNO seedNOlst1 )||(List.mem negseedNO seedNOlst1 ) 
        then ()
        else begin 
            Printf.printf "cannot find seedNO.#%d \n%!" seedNO;
            let badmum = get_mum_from_mumtbl seedNO mum_tbl seed2pos_tbl in
            print_mum badmum false true;
            let poslst = badmum.positions in
            List.iter (fun mi ->
                let seqNO = mi.sequence_NO 
                and pos = mi.left_end in
                let record_lst = Hashtbl.find pos2seed_tbl (seqNO,pos) in
                Printf.printf "(seqNO=%d,pos=%d):%!" seqNO pos;
                List.iter (fun (seed,seed_weight,dir) ->
                Printf.printf "seed#%d,weight=%d,dir=%d; \n%!" seed seed_weight dir;
                ) record_lst;
            )poslst;
        end;
    ) seedNOlst2;
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
let get_range_of_a_lcb seedNOlst seqNO mum_tbl seed2pos_tbl= 
    let debug = false in
    let min_abs_code = List.hd (List.sort (fun x y -> compare (abs x) (abs y) )
    seedNOlst) in
    let ori =  if (min_abs_code>0) then 1 else (-1) in 
    let weight = List.length seedNOlst in
    let seedNO_first = abs (List.hd seedNOlst) in
    let seedNO_last = abs (List.nth seedNOlst (weight-1)) in
    let mumfirst = get_mum_from_mumtbl seedNO_first mum_tbl seed2pos_tbl 
    and mumlast = get_mum_from_mumtbl seedNO_last mum_tbl seed2pos_tbl in
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
    if debug then begin
    Printf.printf "get range of lcb: \n%!";
    Printf.printf "=> %d,%d\n %!" leftend rightend;
    end;
    (leftend,rightend),ori

(*score bwteen two sequences based on hodx matrix *)
let get_score_from_2seq sequence1 sequence2 ori = 
    let debug = false in
    let seq1,seq2 = 
        if ori=1 then sequence1,sequence2
        else 
            rev_comp_arr sequence1, sequence2
    in      
    assert ( (Array.length seq1)=(Array.length seq2) );
    if debug then begin
    Printf.printf "get score from 2seq\n%!";
    print_int_arr seq1; print_int_arr seq2;
    end;
    Array_ops.fold_right_2 (fun acc x y -> 
        let codearr1 = break_code x 
        and codearr2 = break_code y in
        if debug then begin 
            Printf.printf "x=%d,y=%d,code1arr1/codearr2=\n%!" x y;
        Utl.printIntArr codearr1; Utl.printIntArr codearr2;
        end;
        let maxv = ref (-125) in
        Array.iteri (fun idx1 code1 ->
            Array.iteri (fun idx2 code2 ->
                if code1>0 && code2>0 then
                    if !maxv< hodx_matrix.(idx1).(idx2) then
                        maxv := hodx_matrix.(idx1).(idx2)
            )codearr2
        )codearr1;
        if debug then Printf.printf "acc(%d)+%d=%d\n%!" acc !maxv (!maxv+acc);
        acc + !maxv (*hodx_matrix.(atgc_2_idx x).(atgc_2_idx y) *)
    ) 0 seq1 seq2
        
(* get_score returns score of mum, and yes this works for more than 2
* input sequence, but we only care about lcbs between 2 sequences *)
let get_mum_score old_seqarr mum =
    let poslst = mum.positions in
    let first_range = List.hd poslst in
    let first_seq = 
        get_sub_seq2 old_seqarr.(first_range.sequence_NO) 
        first_range.left_end (first_range.right_end-first_range.left_end+1) in
    let first_ori = first_range.orientation in
    let score,_,_ = 
        List.fold_left (fun (sum_score,previous_seq,previous_ori) lcbrange ->
        let seqNO,lend,rend,current_ori =
            lcbrange.sequence_NO,lcbrange.left_end,
            lcbrange.right_end,lcbrange.orientation
        in
        let current_seq = 
            get_sub_seq2 old_seqarr.(seqNO) lend (rend-lend+1) in
        let ori = if (current_ori=previous_ori) then 1 else (-1) in
        if (Array.length previous_seq)<>(Array.length current_seq) then
            print_mum mum false true;
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

(*return true if lcb_record is a light weight lcb.*)    
let is_light_weight_lcb lcb_record in_seq_size_lst =
    (*if (List.length lcb_record.seedNOlst)<minimum_lcb_weight then
        true   else false*)
    let covR = !minimum_cover_ratio in
    let avg_seq_len = get_avg_of_intlst in_seq_size_lst in
    if (lcb_record.avg_range_len < (int_of_float (avg_seq_len *. covR))) 
    then true
    else false

let is_low_score_lcb item =
    if item.ratio < !minimum_lcb_ratio then true
    else false


(*[get_num_of_low_ratio_lcb] returns number of low lcb_ratio. is has_high_covR
* is set to true, only cares about lcbs with high covR and low lcb_ratio*)
let get_num_of_low_ratio_lcb lcb_tbl min_ratio in_seq_size_lst has_high_covR =
    let res = ref 0 in
    Hashtbl.iter (fun key record ->
        let check_c = 
            is_light_weight_lcb record in_seq_size_lst in
        if (check_c=false)&&(record.ratio < min_ratio) 
        then begin
            res:= !res +1;
        end;
    ) lcb_tbl;
    !res

(*this is a test to replace the old way we score lcb*)
let get_new_score lcb_tbl seqsize_lst =
    let debug = false in
    if debug then Printf.printf "Test .... new score : %!";
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
    if debug then 
        Printf.printf "avg_R = %f, avg_cov = %f, R*cov=%f, bk_number = %d\n%!" 
        avg_R avg_cov (avg_cov*.avg_R) bk_number



(*this only works for two sequences ,for now*)
(*  List.nth lcbs n is the lcb list, each lcb is a list of mums. 
* so lcb list looks like this; [mum a1,mum a2,...],[mum b1,mum b2,...],,..
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
let update_lcbs_range_and_score lcbs mum_tbl seed2pos_tbl lcb_tbl  = 
    let debug = false in
   (*fill in range of each lcb with mum_tbl*)
    let idx = ref (-1) in
    List.iter (fun lcblist ->
        idx := !idx + 1;
         List.iter (fun record ->
            let (lend,rend),ori = get_range_of_a_lcb record !idx mum_tbl
            seed2pos_tbl in
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
                score=0; ratio = 0.0; ref_code = (-1); avg_range_len = 0 }
            end;
        ) lcblist
    )lcbs;
    (* Printf.printf "check lcb range lstlst:\n%!";  print_lcbs_range lcb_range_lst_lst; *)
   (*fill in score for each lcb*)
    Hashtbl.iter (fun seedNOlst oldlcb ->
        let total_score = 
            List.fold_left (fun acc seedNO ->
                 let mumi = get_mum_from_mumtbl seedNO mum_tbl seed2pos_tbl in
                 acc + mumi.mumscore
            ) 0 seedNOlst
        in
        let sorted_range_lst = List.rev ( List.sort (fun x y ->
            compare (x.right_end-x.left_end) (y.right_end-y.left_end)
        ) oldlcb.range_lst ) in
        let avg_range = 
            get_avg_of_intlst 
            (List.map (fun x -> x.right_end-x.left_end+1) sorted_range_lst) in
        let ratio = (float total_score)/.(avg_range) in
        if debug then begin
            Printf.printf "update lcb_tbl on %!" ; print_lcb oldlcb;
            Printf.printf "score <-- %d,ratio <-- %f\n%!" total_score ratio;
        end; 
        Hashtbl.remove lcb_tbl seedNOlst;
        Hashtbl.add lcb_tbl seedNOlst 
        {oldlcb with 
                     range_lst = sorted_range_lst;
                     score = total_score; 
                     ratio = ratio;
                     avg_range_len = int_of_float avg_range};

    ) lcb_tbl

(*[get_lcbs_range returns the range lst of current lcb, 
* if special_case=true, we don't care the lcb weight.]*)
let get_lcbs_range lcb_tbl in_seq_size_lst = 
    let lcb_range_lst_arr = Array.make 2 [] in
    Hashtbl.iter (fun seedNOlst record ->
        if (is_light_weight_lcb record in_seq_size_lst)=false then
            List.iter (fun mi ->
                let oldlst = lcb_range_lst_arr.(mi.sequence_NO) in
                let newlst = (mi.left_end,mi.right_end)::oldlst in
                lcb_range_lst_arr.(mi.sequence_NO) <- newlst ;
            ) record.range_lst;
    ) lcb_tbl;
    let lcb_range_lst_arr = (*sort the range by leftend*)
        Array.map (fun x -> 
            List.sort (fun (aleft,_) (bleft,_) -> compare aleft bleft) x
        ) lcb_range_lst_arr
    in
    Array.to_list lcb_range_lst_arr

(* [get_lcb_key_by_range] get the lcb_key and lcb_code by its range. if no lcb is found,
* return [],0. 
* remember what we did in [fill_in_indel]. the range in full_range_lstlst could
* be different from the original lcb range. *)
let get_lcb_key_by_range seqNO range lcb_tbl =
    let (leftend,rightend) = range in
    let reskey= ref [] and rescode = ref 0 in (*-1 is not a good idea for the
    rescode could be (code=1)*(ori=-1)=-1 *)
    Hashtbl.iter (fun key record ->
        List.iter (fun mi ->
            let left = mi.left_end and right = mi.right_end in
            if ((leftend-1=left)||(leftend=left))&&
            ((rightend+1=right)||(rightend=right))&&
            (mi.sequence_NO=seqNO) then begin
                reskey := key;
                rescode := record.ref_code*mi.orientation;
            end
        ) record.range_lst;
    ) lcb_tbl;
    !reskey,!rescode


(* no worry about this, alignment function will take care of this.
* what if just in one sequence, all lcb are adjacent to each another? 
* say we have N lcb blocks in one sequence, in this case, one sequence 
* will have just N blocks, while the other one has 2*N+1 blocks 
let fill_in_indel full_range_lstlst = 
    let rangelst1 = List.hd full_range_lstlst 
    and rangelst2 = List.nth full_range_lstlst 1 in
    let size1 = List.length rangelst1 in
    let size2 = List.length rangelst2 in
    let fill_in_shorter_rangelst lsty sizex sizey =
        if (sizex-sizey)=2 then
            let (lastl,lastr) = List.hd (List.rev lsty) in
            let (firstl,firstr) = List.hd lsty in
            let lst2keep = List.tl (List.rev (List.tl lsty)) in
            ([(0,0);(1,firstr)] @ lst2keep @ [(lastl,lastr-1);(lastr,lastr)])
        else if (sizex-sizey)=1 then
            let (firstl,firstr) = List.hd lsty in
            let lst2keep = List.tl lsty in
           ( [(0,0);(1,firstr)]@lst2keep )
        else begin
            Printf.printf "fill_in_gap_seq, sizex=%d,sizey=%d\n%!;" sizex sizey;
            assert(false);
        end
    in
    if size1=size2 then full_range_lstlst
    else if size1>size2 then
        let newlst2 = fill_in_shorter_rangelst rangelst2 size1 size2 in
        [rangelst1;newlst2]
    else 
        let newlst1 = fill_in_shorter_rangelst rangelst1 size2 size1 in
        [newlst1;rangelst2]
*)
(*given a lcbs range list list, get_full_range_lstlst returns a range list list
* of both lcbs and area outside lcbs. *)
let get_full_range_lstlst lcb_range_lstlst in_seqsize_lst =
    let res = List.map2 (fun lcb_range_list total_size ->
        let last_rightend,acc=
            List.fold_left (fun (pre_rightend,acc) (leftend,rightend) ->
                if pre_rightend>leftend then
                    Printf.printf "%d>%d\n%!" pre_rightend leftend;
                assert(pre_rightend<=leftend);
                let newacc = 
                    (*if leftend=0 then acc@[(0,0)]
                    else*) if leftend>0 then   
                    acc@[((pre_rightend+1),(leftend-1))] 
                    else acc
                in
                (*let newleftend = if leftend=0 then 1 else leftend in*)
                let newacc = newacc@[(leftend,rightend)] in
                rightend,newacc 
                ) (-1,[]) lcb_range_list
        in
        let res = 
            if (last_rightend+1=total_size) then acc
                (*let last_l,_ = List.hd (List.rev acc) in
                let acc2keep = List.rev(List.tl (List.rev acc)) in
                acc2keep@[(last_l,last_rightend-1);(last_rightend,last_rightend)]
                *) 
            else acc@[(last_rightend+1,total_size-1)]
        in
        res
    ) lcb_range_lstlst in_seqsize_lst in
    (*fill_in_indel*) res

    
let get_seq_outside_lcbs old_seqarr old_seq_size_lst lcb_tbl =
    let debug = false in
    if debug then Printf.printf "get seq outside lcbs :\n%!";
    let lcb_range_lst_lst = get_lcbs_range lcb_tbl old_seq_size_lst in
    if debug then print_lcbs_range lcb_range_lst_lst;
    let seq_outside_lcb_arr = Array_ops.map_2 (fun seq lcb_range ->
        let seqlen = Array.length seq in
        let last_rightend,tmpseq = 
            List.fold_left (
                fun (pre_rightend,accseq) (leftend,rightend) ->
                assert(pre_rightend<=leftend);
                let size = leftend-pre_rightend-1 in
                let newsubseq =
                    if size>0 then
                    Array.sub seq (pre_rightend+1) (leftend-pre_rightend-1) 
                    else [||]
                in
                rightend,Array.append accseq newsubseq
            ) (-1,[||]) lcb_range 
        in
        let last_size = seqlen-1-last_rightend in
        let subend = 
            if last_size>0 then
            Array.sub seq (last_rightend+1) (seqlen-1-last_rightend) 
            else [||]
        in
        Array.append tmpseq subend
    ) old_seqarr (Array.of_list lcb_range_lst_lst) in
    lcb_range_lst_lst, seq_outside_lcb_arr
    

(*[get_lcb_covR_num_badlcb] return the cover rate of qualified lcb, 
* number of un-qualified lcb.  
* but if we only have one lcb in lcb_tbl, qualified or not, 
* we get the cover rate of that lcb, num_badlcb doen't matter in this case*)    
let get_lcb_covR_num_badlcb lcb_tbl in_seq_size_lst =
    let q_cov_len = ref 0 in
    let num_badlcb = ref 0 in
    let lcb_tbl_size = Hashtbl.length lcb_tbl in
    Hashtbl.iter (fun key lcb_record ->
        if (lcb_tbl_size=1) (*we keep the last bad lcb*) ||
        (((is_light_weight_lcb lcb_record in_seq_size_lst)=false)
        && (is_low_score_lcb lcb_record)=false)
        then q_cov_len := !q_cov_len + lcb_record.avg_range_len
        else 
            num_badlcb := !num_badlcb + 1
    ) lcb_tbl;
    let avg_seq_len = get_avg_of_intlst in_seq_size_lst in
    (float !q_cov_len)/.avg_seq_len , !num_badlcb
    

(* LCB = local collinear blocks, works only for two sequences now *)
let build_LCBs seedNOlstlst mum_tbl seed2pos_tbl  in_seq_size_lst
special_case =
    let debug = false and debug_bk_matrix = false in
    if debug then Printf.printf "build LCBS \n %!";
    let resmatrix = Array.of_list (List.map (fun lst -> Array.of_list lst)
    seedNOlstlst) in 
    let bkmatrix,_ = get_break_point_matrix resmatrix debug_bk_matrix in
    if debug_bk_matrix then print_int_matrix bkmatrix;
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
                failwith "we only record 1 or 0 in bk_matrix";
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
    let lcb_tbl = Hashtbl.create init_tb_size in
    update_lcbs_range_and_score lcbs mum_tbl seed2pos_tbl lcb_tbl ;
    let q_cov_rate,num_badlcb = 
        get_lcb_covR_num_badlcb lcb_tbl in_seq_size_lst in
    if debug then 
        Printf.printf "q_cov_rate = %f,end of building LCBs\n%!" 
        q_cov_rate;
    lcbs,q_cov_rate,num_badlcb,lcb_tbl

(* remove one lcb at the time, the rest of lcbs might collapse into each other
* return the score, cover_rate, lcbs and lcb_tbl*)
let filter_out_one_lcb seedNOlst removed_lcb_score lcbs total_mum_score bk_penalty
in_seq_size_lst seed2pos_tbl mum_tbl =
    assert (seedNOlst<>[]);
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
    let new_lcbs, new_lcb_cov_rate,num_badlcb,new_lcb_tbl = 
         build_LCBs seedNOlstlst mum_tbl seed2pos_tbl in_seq_size_lst false in
    new_score,new_lcb_cov_rate,num_badlcb,new_lcbs,new_lcb_tbl
    


let get_worst_lcb_from_light_lcbtbl param_shortcut light_lcb_tbl lcbs mum_tbl 
seed2pos_tbl in_seq_size_lst param_init_value =
    let debug = false in
    let shortcut,lcb_w_worstR,worst_score = param_shortcut in
    let bk_penalty,total_mum_score,init_higher_score,init_better_covR,init_low_ratio_lcb_num = param_init_value in
    let lcb_to_remove = ref [] in
    let lcbs_after_remove = ref [[[]]] in
    let higher_score = ref init_higher_score in
    let better_covR = ref init_better_covR in
    let num_badlcb = ref 0 in
    let lcb_tbl_after_remove = ref (Hashtbl.create init_tb_size) in
        if shortcut then begin (*Case 1: shortcut, just remove this one*)
        if debug then begin
            Printf.printf "take a shortcut, get rid of lcb:%!";
            print_int_list lcb_w_worstR; 
            Printf.printf "\n%!";
        end;
        let new_score,new_lcb_cov_rate,new_num_badlcb,new_lcbs,new_lcb_tbl =
            filter_out_one_lcb lcb_w_worstR worst_score lcbs total_mum_score bk_penalty
            in_seq_size_lst seed2pos_tbl mum_tbl
        in
        higher_score := new_score;
        lcb_to_remove := lcb_w_worstR;
        lcbs_after_remove := new_lcbs;
        better_covR := new_lcb_cov_rate;
        num_badlcb := new_num_badlcb;
        lcb_tbl_after_remove := Hashtbl.copy new_lcb_tbl;
    end
    else  (*Case 2: iter all low W lcb to find the worst one.*)
        Hashtbl.iter (fun seedNOlst lcb ->
        if debug then begin
            Printf.printf "if we remove lcb:%!"; print_lcb lcb;
        end;
        let new_score,new_lcb_cov_rate,new_num_badlcb,new_lcbs,new_lcb_tbl =
            filter_out_one_lcb lcb.seedNOlst lcb.score lcbs total_mum_score bk_penalty
            in_seq_size_lst seed2pos_tbl mum_tbl
        in
        let num_low_ratio_lcb = get_num_of_low_ratio_lcb new_lcb_tbl
        !minimum_lcb_ratio in_seq_size_lst true in
        if debug then
            Printf.printf "new score=%d,new covR=%f,num_lowR_lcb(%d) should not increase\n%!" 
            new_score new_lcb_cov_rate num_low_ratio_lcb;
        (* rules of iteration:
        *  1. don't create more lowR lcbs 
        *  2. if cover rate increase, a decrease in score is ok
        *  3. if cover rate does not increase, we take the biggest increase in score
        *  4. don't let cover rate drop*)
        if (init_low_ratio_lcb_num >= num_low_ratio_lcb) &&
        ( (new_lcb_cov_rate > !better_covR)
        ||(new_lcb_cov_rate = !better_covR)&&(new_score > (!higher_score)))
        then begin
            higher_score := new_score;
            better_covR := new_lcb_cov_rate;
            lcb_to_remove := seedNOlst;
            lcbs_after_remove := new_lcbs;
            num_badlcb := new_num_badlcb;
            lcb_tbl_after_remove := Hashtbl.copy new_lcb_tbl;
        end;
        ) light_lcb_tbl
    ;
    if debug then 
        Printf.printf "end of get_worst_lcb_from_light_lcbtbl,tbl size=%d\n%!"
        (Hashtbl.length !lcb_tbl_after_remove);
    !higher_score,!better_covR,!num_badlcb,!lcb_to_remove,
    !lcbs_after_remove,!lcb_tbl_after_remove
    
 

(*get_worst_lcb returns the worst lcb so we can remove it to get a better score.
* also returns the score, number of low W lcb left if we
* remove the lcb*)
let get_worst_lcb lcbs lcb_tbl bk_penalty in_seq_size_lst mum_tbl seed2pos_tbl
shortcut init_cov_rate =
    let debug = false in
    let total_num = Hashtbl.length lcb_tbl in
    (*num_high_W_low_R : make sure we won't create more low ratio lcb.*)
    let num_high_W_low_R = ref 0 and total_mum_score = ref 0 in
    Hashtbl.iter (fun seedNOlst record ->
            total_mum_score := !total_mum_score + record.score;
            if (is_light_weight_lcb record in_seq_size_lst)=false &&
               (is_low_score_lcb record) then
                   num_high_W_low_R := !num_high_W_low_R + 1;
    ) lcb_tbl;
    let total_mum_score = !total_mum_score in
    let init_bk_num = List.length (List.hd lcbs) in
    let init_score =  total_mum_score - init_bk_num*bk_penalty in
    if debug then begin
        Printf.printf "\n<-------get worst lcb, bk_num=%d,init score=%d,init covR =%f;%!" 
        init_bk_num init_score init_cov_rate;
        Printf.printf "total_mum_score = %d,lowR&highW lcb num = %d;total num = %d\n%!" 
        total_mum_score !num_high_W_low_R total_num;
    end;
    let lcb_to_remove = ref [] and tmp_score = ref 0 in
    Hashtbl.iter (fun key record -> (*remove highW lcb with lowR first*)
            if (!lcb_to_remove=[])&&(is_low_score_lcb record)&&
            ((is_light_weight_lcb record in_seq_size_lst)=false)
            then begin
                lcb_to_remove := key; tmp_score := record.score;
            end
    ) lcb_tbl;
    (*Case 1: remove high w but low R lcb first*)
    if (!lcb_to_remove<>[]) then begin 
        if debug then Printf.printf "remove a lowR&highW lcb first\n%!";
        let new_score,new_cov_rate,num_badlcb,new_lcbs,lcb_tbl_after_remove =
            filter_out_one_lcb !lcb_to_remove !tmp_score lcbs total_mum_score 
            bk_penalty in_seq_size_lst seed2pos_tbl mum_tbl
        in
        !lcb_to_remove,new_score, new_cov_rate, 
        num_badlcb,
        new_lcbs,lcb_tbl_after_remove
    end
    (*Case 2: remove the worst lcb from the rest low w lcbs *)
    else begin
        let worst_ratio = ref !minimum_lcb_ratio 
        and lcb_w_worstR = ref [] 
        and worst_score = ref 0 in 
        (*we don't really need a tbl, just keep a list of lowW lcb key and score
        * should be enough*)
        let light_lcb_tbl = Hashtbl.create init_tb_size in
        Hashtbl.iter (fun key record ->
            if shortcut && (record.ratio < !worst_ratio) then begin
                    lcb_w_worstR := record.seedNOlst;
                    worst_ratio := record.ratio;
                    worst_score := record.score;
            end;
            if (is_light_weight_lcb record in_seq_size_lst) then 
                Hashtbl.add light_lcb_tbl key record
        ) lcb_tbl;
        if shortcut&&debug then
            Printf.printf "worst_ratio=%f,worst_score=%d\n%!" !worst_ratio !worst_score;
        if (Hashtbl.length light_lcb_tbl)>1 then begin
            if debug then Printf.printf "choose one from light_lcb_tbl to remove\n%!";
            let shortcut =
                if shortcut&&(!lcb_w_worstR=[]) then begin
                    if debug then Printf.printf 
                    "no shortcut to take, try all bad lcb 1 by 1.\n%!";
                    false
                end
                else shortcut
            in
            let param_shorcut = (shortcut,!lcb_w_worstR,!worst_score) in
            let param_init_value = (bk_penalty,total_mum_score,init_score,
            init_cov_rate,!num_high_W_low_R) in
            let higher_score,higher_cov_rate,num_badlcb,lcb_to_remove,lcbs_after_remove,
            lcb_tbl_after_remove = 
                get_worst_lcb_from_light_lcbtbl param_shorcut 
                light_lcb_tbl lcbs mum_tbl seed2pos_tbl in_seq_size_lst 
                param_init_value
            in
            if debug then begin
                Printf.printf "RESULT: remove lcb %!"; 
                print_int_list lcb_to_remove;
                Printf.printf "score after remove lcb is %d, cov rate = %f\n%!" 
                higher_score higher_cov_rate;
            end;
            (*let light_lcb_num = Hashtbl.length light_lcb_tbl in*)
            lcb_to_remove,higher_score,higher_cov_rate,num_badlcb,(*light_lcb_num,*)
            lcbs_after_remove,lcb_tbl_after_remove
        end
        else 
            (*in this case cov_rate does not matter, light_lcb_left does not matter*)
            [],init_score, 0., 
            1 (*so we update lcb_tbl in remove_light_weight_lcbs*), 
            lcbs, lcb_tbl
        ; (*end of if (Hashtbl.length light_lcb_tbl)>1 *)
    end



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

(*remove bad lcbs (W is weight, R is ratio)  
* -- high W low R lcbs
* -- low W lcbs
* remove highW lowR lcb first 
* the result lcb_tbl could still have low W lcbs, if remove them will create
* more high W low R lcb, we don't do it -- this is different from Mauve *)
let remove_bad_lcbs lcbs lcb_tbl mum_tbl seed2pos_tbl in_seq_size_lst
num_of_mums =
    let debug = true and debug2 = false in
    let bk_penalty = get_break_point_penalty in_seq_size_lst num_of_mums in 
    if debug then begin
        Printf.printf "start to remove lowW or lowR lcbs (min_cov=%f,min_R=%f) \
        inputlcbs is\n%!" !minimum_cover_ratio !minimum_lcb_ratio;
        if debug2 then print_lcblst lcbs;
        if debug2 then Printf.printf "qualified lcb are :\n%!";
    end;
    (*light_lcb_tbl is for low W and low R lcbs
    let light_lcb_tbl = Hashtbl.create init_tb_size in *)
    let init_num_badlcb = ref 0 in
    let total_mum_score = ref 0 in (*we need to pass an init score*)
    Hashtbl.iter (fun key record ->
        total_mum_score := !total_mum_score + record.score;
        if (is_light_weight_lcb record in_seq_size_lst)||
        (is_low_score_lcb record)
        then  
            init_num_badlcb := !init_num_badlcb + 1
        else  
            if debug2 then print_lcb record
    ) lcb_tbl;
    (*pass the init score/bk_num/etc to while loop*)
    let total_mum_score = !total_mum_score in
    let init_bk_num = List.length (List.hd lcbs) in
    let init_score = total_mum_score - init_bk_num*bk_penalty in
    let res_score = ref init_score in
    let res_cov_rate,res_num_badlcb = 
        get_lcb_covR_num_badlcb lcb_tbl in_seq_size_lst in
    let res_cov_rate = ref res_cov_rate in
    let res_lcb_tbl = ref lcb_tbl in
    let res_lcbs = ref lcbs in
    let res_num_badlcb = ref res_num_badlcb in
    let sign = ref ( !res_num_badlcb > 1 ) in
    while (!sign) do
        if debug then 
        Printf.printf "\n in while loop,current score=%d,covR=%f,num_badlcb=%d\n %!" 
        !res_score !res_cov_rate !res_num_badlcb;
    (*if we have too many bad matches, we mneed a shortcut to remove them quickly*)
    let num_alllcb = (Hashtbl.length !res_lcb_tbl) in
    let shortcut = 
            if !res_num_badlcb=0 then false
            else if (num_alllcb / !res_num_badlcb)<2 then true
            else false
    in
    if debug then begin
        Printf.printf "we have %d lightW or lowR lcb (among %d),shortcut<-%b \n%!"
            !res_num_badlcb num_alllcb shortcut; 
    end;
        let seedNOlst_to_remove,score,cov_rate,num_badlcb,lcbs_after_remove,
        lcb_tbl_after_remove = 
            get_worst_lcb !res_lcbs !res_lcb_tbl bk_penalty in_seq_size_lst
            mum_tbl seed2pos_tbl shortcut !res_cov_rate in
        if debug then begin
            Printf.printf "if we remove lcb : %!";
            print_int_list seedNOlst_to_remove; 
            Printf.printf "new score will be %d,cov_rate=%f, bad lcb num=%d\n%!" 
            score cov_rate num_badlcb;
        end;
        (*mauve stops when there is only one light weight lcb left*)
        if (seedNOlst_to_remove<>[])&&(num_badlcb>1)&&(!res_score <= score) 
        then begin
            sign := true; 
            res_lcb_tbl := Hashtbl.copy lcb_tbl_after_remove;
            res_lcbs := lcbs_after_remove;
            if debug2 then begin
                Printf.printf "current lcbs is :%!";
                print_lcblst !res_lcbs;
            end; 
            res_score := score;
            res_cov_rate := cov_rate;
            res_num_badlcb := num_badlcb;
            if debug then 
                Printf.printf "new score >= old score (cov_R set to %f), \
                continue with loop\n%!" !res_cov_rate;
        end
        else begin
            sign := false;
            if (seedNOlst_to_remove<>[])&&(num_badlcb=1)&&(!res_score <= score)
            then begin
                if debug then Printf.printf "only one lowW lcb left,%!";
                res_lcb_tbl := Hashtbl.copy lcb_tbl_after_remove;
                res_lcbs := lcbs_after_remove;
                if debug2 then begin
                    Printf.printf "current lcbs is :%!";
                    print_lcblst !res_lcbs;
                end; 
                res_score := score;
                res_cov_rate := cov_rate;
                res_num_badlcb := num_badlcb;
            end;
            if debug then 
                Printf.printf ",get out of loop\n%!";
        end;
    done; (* while(!sign) *)
    if debug then 
        Printf.printf "end of remove_light_weight_lcbs,res_covR = %f,\n%!"
        !res_cov_rate;
    !res_lcbs,!res_cov_rate, !res_lcb_tbl
    

let update_k_seed_lst pre_jseedNO previous_k_lst pre_kposlst k_pos j_mumseq trim_from_left
weight_reduce pos2seed_tbl1 pos2seed_tbl2 seed2pos_tbl mum_tbl =
let debug = false in
    let k_seedNO = get_a_seedNO () in
    let k_extsign = 0 in (*since mi and mj are both ext=0*)
    let k_milst = k_pos::pre_kposlst in
    let k_size = (List.length k_milst) in
    if debug then 
        Printf.printf "we have a new seedk : %d,W=%d\n%!" k_seedNO weight_reduce;
    let k_mumkey = get_mumkey_from_milst k_milst in
    if debug then begin
        Printf.printf "kmumkey = %d\n%!" k_mumkey;
        print_int_arr j_mumseq;
    end;
    let start_pos = 
        if trim_from_left then 0 
        else (Array.length j_mumseq)-weight_reduce
    in
    let k_mumseq = get_sub_seq2 j_mumseq start_pos weight_reduce in  
    if debug then print_int_arr k_mumseq;
    let new_kmum =  
        {seedNO=k_seedNO; 
        positions=k_milst; 
        mumkey = k_mumkey;
        mumseq = k_mumseq; 
        size= k_size;  
        neighborhood_lst = []; 
        subsuming_pointer= -1; 
        extendable = k_extsign; 
        priority_lst=[];
        mumscore = 0 } in  
    if debug then print_mum new_kmum false true;
    let sign_newmum = add_mum_to_mumtbl new_kmum mum_tbl in
    if sign_newmum then begin
    let poslst = List.map (fun record ->
        record.sequence_NO,record.left_end,record.orientation
    )  k_milst in
    Hashtbl.add seed2pos_tbl k_seedNO (poslst,k_mumseq,weight_reduce);
    List.iter (fun record ->
        let left=record.left_end and right = record.right_end in
        let klen = right - left + 1 in
        let seqNO = record.sequence_NO 
        and ori = record.orientation in
        let pos2seed_tbl_l,pos2seed_tbl_r =
            if (trim_from_left) then pos2seed_tbl1,pos2seed_tbl2
            else pos2seed_tbl2,pos2seed_tbl1
        in
        add_to_pos2seed_tbl trim_from_left pos2seed_tbl_l pos2seed_tbl_r 
        seqNO left k_seedNO klen 
        ori k_extsign k_size seed2pos_tbl mum_tbl;
        add_to_pos2seed_tbl trim_from_left pos2seed_tbl_l pos2seed_tbl_r 
        seqNO right k_seedNO klen 
        ori k_extsign k_size seed2pos_tbl mum_tbl;
    ) k_milst;
    end;
    List.filter (fun (pre_j,_,_)-> pre_j<>pre_jseedNO) previous_k_lst





(* [update_tables] trim mumj, and create mumk out of the trimed parts of
* mumj if necessary. return the list of other trimed parts, also if jmum become
* redundant after trim (which means it's the same as some mum we already have),
* return true to remove it, else return false.
* for mumj, we remove the old record of oldpos of pos2seed_tbl1, add newpos
* as  new position. in pos2seed_tbl2, we only update the weight of the old record
* on pos2 .
* for mumk, we wait until we have both subseqs, then create it.*)
let update_tables mum_tbl mumj seqNO pos2seed_tbl1 oldpos newpos pos2seed_tbl2
pos2 newweight seed2pos_tbl new_mi oldori weight_reduce min_weight trim_from_left 
previous_k_lst = 
    let seedNO2upgrade_lst = ref [] in
    let jseedNO,j_mumseq,j_extsign,j_size =
        mumj.seedNO,mumj.mumseq,mumj.extendable,mumj.size in
    let debug = false in
    if debug then 
    Printf.printf "update tables on seed#%d,seqNO=%d,oldpos=%d,newpos=%d,newweight=%d weight_reduce=%d,pos2=%d trim_from_left:%b\n%!"
    jseedNO seqNO oldpos newpos newweight weight_reduce pos2 trim_from_left ;
    let new_mumseq =
        if (seqNO=0)&&trim_from_left then
            get_sub_seq2 j_mumseq weight_reduce (pos2-newpos+1)
        else if (seqNO=0) then
            get_sub_seq2 j_mumseq 0 (newpos-pos2+1)
        else j_mumseq
    in 
    let sign_newmum = 
        update_positions_by_seqNO jseedNO seqNO mumj new_mi new_mumseq mum_tbl
    in
    let left_tbl,right_tbl = 
        if trim_from_left then pos2seed_tbl1,pos2seed_tbl2
        else pos2seed_tbl2,pos2seed_tbl1
    in
    let upgrade_mums seedNO2upgrade_lst =
        upgrade_extendable_pos2seed_tbl seedNO2upgrade_lst left_tbl
        right_tbl seed2pos_tbl mum_tbl;
    in
    if sign_newmum then begin (*trim mumj, update positions of mumj *)
        if trim_from_left then  (*update mumseq and position*)
            update_seed2pos_tbl jseedNO (Some(seqNO,oldpos,newpos)) new_mumseq
            seed2pos_tbl
        else (*only update mumseq*)
            update_seed2pos_tbl jseedNO None new_mumseq seed2pos_tbl;
        let reslst = remove_from_pos2seed_tbl pos2seed_tbl1 seqNO oldpos mumj mum_tbl 
        seed2pos_tbl in
        if reslst<>[] then 
            seedNO2upgrade_lst := !seedNO2upgrade_lst@reslst;
        if !seedNO2upgrade_lst <> [] then  
            upgrade_mums !seedNO2upgrade_lst;
        add_to_pos2seed_tbl trim_from_left left_tbl right_tbl seqNO newpos jseedNO newweight oldori j_extsign j_size seed2pos_tbl mum_tbl;
        modify_record_in_pos2seed_tbl pos2seed_tbl2 seqNO pos2 jseedNO newweight oldori;
        if debug then Printf.printf "start working on seedK lst:\n%!";
        if weight_reduce>min_weight then begin
            let k_left,k_right = 
                if (trim_from_left) then oldpos,oldpos+weight_reduce-1
                else oldpos-weight_reduce+1,oldpos 
            in
            let k_pos = { sequence_NO = seqNO; left_end = k_left; right_end =
                k_right; orientation = oldori } in
            (* mi and mj are both extendable=0, there will be at most two subseq 
            * for mk, one in seqNO=0, one in seqNO=1 . 
            * also there won't be more than 1 prev_k record fit this*)
            let pre_jseedNO,pre_kposlst,pre_jmumseq = 
                try (List.find (fun (pre_j,_,_) -> pre_j = jseedNO)) previous_k_lst 
                with | Not_found -> (0,[],[||])
            in
            if debug then
                    List.iter (fun (pre_j,kposlst,_) ->
                            Printf.printf "jseed=%d,%!" pre_j;
                            List.iter (fun kpos ->
                                Printf.printf "(%d,%d,%d,%d)%!"
                                kpos.sequence_NO kpos.left_end kpos.right_end 
                                kpos.orientation
                            ) kposlst;
                            Printf.printf "\n%!";
                    ) previous_k_lst;
            let ori_jmumseq =
                if (Array.length pre_jmumseq)>(Array.length j_mumseq) then
                    pre_jmumseq else j_mumseq
            in
            let current_k_lst =
                if (pre_jseedNO<>0) then 
                update_k_seed_lst pre_jseedNO previous_k_lst 
                pre_kposlst k_pos ori_jmumseq trim_from_left weight_reduce 
                pos2seed_tbl1 pos2seed_tbl2 seed2pos_tbl mum_tbl
                else  
                    previous_k_lst@[(jseedNO,[k_pos],ori_jmumseq)] 
            in
            if debug then Printf.printf "end of update tables,add new kseed 2 acclst \n%!";
            current_k_lst,false
        end
        else begin
            if debug then Printf.printf "end of update tables,no new kseed \n%!";
            previous_k_lst,false
        end; (*end of if weight_reduce>min_weight*)
    end (*if sign_newmum *)
    else begin (*jmum is redundant after trimming. remove it completely*)
        if debug then 
            Printf.printf "mum#%d is redundant after trimming.\n%!" mumj.seedNO;
        List.iter (fun record ->
        let seqNO = record.sequence_NO and left = record.left_end 
        and right = record.right_end in
        let reslst1 = remove_from_pos2seed_tbl left_tbl seqNO left mumj mum_tbl
        seed2pos_tbl in
        let reslst2 = remove_from_pos2seed_tbl right_tbl seqNO right mumj mum_tbl 
        seed2pos_tbl in
        if reslst1<>[] then seedNO2upgrade_lst := !seedNO2upgrade_lst@reslst1;
        if reslst2<>[] then seedNO2upgrade_lst := !seedNO2upgrade_lst@reslst2;
        ) mumj.positions;
        if !seedNO2upgrade_lst <> [] then 
            upgrade_mums !seedNO2upgrade_lst;
        Hashtbl.remove seed2pos_tbl mumj.seedNO;
        previous_k_lst,true
    end
    
    
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
let resolve_overlap_mum mum_tbl min_weight pos2seed_tbl_left 
pos2seed_tbl_right seed2pos_tbl debug debug_neighborhood =
    if debug then 
        Printf.printf "resolve overlap mums with min_weight = %d\n%!" min_weight;
    Hashtbl.iter (fun seedNO (_,_,_) ->
        update_neighborhood seedNO mum_tbl seed2pos_tbl pos2seed_tbl_left
        pos2seed_tbl_right debug_neighborhood
    ) seed2pos_tbl;   
    if debug then Printf.printf "update neighborhood first, done\n%!";
    let keep_i multi_i multi_j weight_i weight_j =
        if (multi_i>multi_j)||((multi_i=multi_j)&&(weight_i>=weight_j)) then true
        else false
    in
    (*[get_three_pos] return a.the old pos for tbl1, b.the new pos for tbl1. 
    * c. the pos for tbl2.*)
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
    (*I know function [update_table...] has a long input list, here is why.
    * when we trim mumj from left, we remove the old record(oldpos) for mumj 
    * in pos2seed_tbl_left and add the new one(newpos). 
    * for pos2seed_tbl_right, just update the old record (pos2) with new j_weight. 
    * So here we mark tbl_left as tbl1 and tbl_right as tbl2.
    * When we trim mumj from right, it's another way around. I don't want to
    * write two "update_table" functions for each dir.*)
    let update_table_trim_right_or_left trim_from_left mumj current_mi new_mi
    oldpos newpos pos2 new_jweight weight_reduce min_weight previous_k_lst =
        let record = current_mi in
        if (trim_from_left) then
            update_tables mum_tbl mumj record.sequence_NO pos2seed_tbl_left 
            oldpos newpos pos2seed_tbl_right pos2 new_jweight seed2pos_tbl 
            new_mi record.orientation weight_reduce min_weight true previous_k_lst 
        else 
            update_tables mum_tbl mumj record.sequence_NO pos2seed_tbl_right 
            oldpos newpos pos2seed_tbl_left pos2 new_jweight seed2pos_tbl 
            new_mi record.orientation weight_reduce min_weight false previous_k_lst
    in
    let seedlst = ref [] and removed_seedlst = ref [] in 
    Hashtbl.iter (fun key record -> seedlst := key::!seedlst ) seed2pos_tbl;
    let seedlst = List.rev !seedlst in
    List.iter (fun i_seedNO -> 
        if (List.mem i_seedNO !removed_seedlst)=false then begin
        let debug2 = false in
        if debug2 then Printf.printf "resolve overlap,work on seed#%d\n%!" i_seedNO;
        let mi = get_mum_from_mumtbl i_seedNO mum_tbl seed2pos_tbl in
        if debug2 then print_mum mi true true;
        assert(mi.seedNO=i_seedNO);
        if (mi.extendable=0)&&(mi.subsuming_pointer=(-1)) then begin 
        if (mi.neighborhood_lst=[]) then
            (*in case this mum is a new mum added by function "transpose_back" 
            * when we search area outside lcbs*)
            update_neighborhood i_seedNO mum_tbl seed2pos_tbl pos2seed_tbl_left
            pos2seed_tbl_right false;
        let idx = ref 0 and sign = ref true in
        while (!sign)&&(mi.neighborhood_lst<>[]) do
            (*mi's positions/neighbors could be different from last iteration 
            * trimed by other mj, that's why we need to load mi here.*)
            let mi = get_mum_from_mumtbl i_seedNO mum_tbl seed2pos_tbl in
            if debug2 then begin
                Printf.printf "idx = %d, i_mum is: %!" !idx; print_mum mi true false;
            end;
            let (seqNO,j_seedNO,i_ori,j_ori,d,dr) = List.nth mi.neighborhood_lst !idx in
            if debug2 then Printf.printf "i=%d,j=%d,seqNO=%d\n%!" i_seedNO j_seedNO seqNO;
            let mj = get_mum_from_mumtbl j_seedNO mum_tbl seed2pos_tbl in
            if debug2 then begin 
                Printf.printf "j_mum is :%!"; print_mum mj false true;
            end;
            assert(mj.seedNO=j_seedNO);
            if (mj.extendable=0)&&(mi.extendable=0)&&(mj.subsuming_pointer=(-1)) 
            then begin
                let poslst_i = mi.positions in
                let position_i = get_position_by_seqNO poslst_i seqNO in
                let ileft = position_i.left_end and iright = position_i.right_end
                and iori = position_i.orientation and i_size = mi.size in
                let poslst_j = mj.positions in
                let position_j = get_position_by_seqNO poslst_j seqNO in
                let jleft = position_j.left_end and jright = position_j.right_end 
                and jori = position_j.orientation and j_size = mj.size in
                let leni= iright-ileft+1 and lenj= jright-jleft+1 in
                if debug2 then Printf.printf "i=%d(%d,%d) j=%d(%d,%d),d=%d\n%!" 
                        i_seedNO ileft iright j_seedNO jleft jright d;
                if (jleft > iright)|| (jright < ileft) then begin
                    idx := !idx + 1;
                    if debug2 then Printf.printf "idx ++ = %d\n%!" !idx;
                end
                else begin
                    if d=0 then print_mum mi false true;
                    if d=0 then print_mum mj false true;
                    assert(d<>0);
                    let acc_k_lst = ref [] in
                    idx := 0 ; (*neighborhood of i_mum will change, start over again*)
                    assert ((ileft-jleft)=d);
                    let i_seedNO,j_seedNO,poslst_j,position_j,mj,
                    ileft,iright,jleft,jright,jori,d = 
                        if (keep_i i_size j_size leni lenj) then begin
                            if debug2 then Printf.printf "keep i\n%!";
                            i_seedNO,j_seedNO,poslst_j,position_j,mj,
                            ileft,iright,jleft,jright,jori,d
                        end
                        else begin (*we keep j_mum,trim i_mum instead,*) 
                            if debug2 then Printf.printf "keep j\n%!";
                             j_seedNO,i_seedNO,poslst_i,position_i,mi,
                             jleft,jright,ileft,iright,iori,(0-d)
                        end
                    in
                    if (d<0) then begin
                        let new_jleft = iright + 1 in
                        let new_jweight = jright-new_jleft+1 in
                        let weight_reduce = iright - jleft + 1 in
                        if debug2 then Printf.printf "d<0,new_j_w=%d,%!"
                        new_jweight;
                        if (new_jweight>min_weight) then 
                            let continue_sign = ref true in
                            List.iter (fun record ->
                                if !continue_sign then begin
                                let trim_from_left = 
                                    if (record.orientation=jori) then true else false 
                                in
                                let oldpos,newpos,pos2 =
                                    get_three_pos trim_from_left record weight_reduce in
                                let new_mi = 
                                    get_new_mi trim_from_left record newpos
                                in
                                (*mj might be diff because the previous inter*)
                                let new_mj = get_mum_from_mumtbl j_seedNO
                                mum_tbl seed2pos_tbl in
                                let current_klst,removej = 
                                update_table_trim_right_or_left trim_from_left new_mj 
                                record new_mi oldpos newpos pos2 new_jweight 
                                weight_reduce min_weight !acc_k_lst in
                                acc_k_lst := current_klst;
                                if removej then begin
                                    continue_sign := false;
                                removed_seedlst := mj.seedNO::(!removed_seedlst);
                                end;
                                if debug2 then
                                Printf.printf "add seed#%d to removelst1\n%!" mj.seedNO;
                                end;
                            ) poslst_j
                        else 
                            begin (*no need to carry j_seed with us*) 
                            remove_seed2 mj pos2seed_tbl_left 
                            pos2seed_tbl_right mum_tbl seed2pos_tbl true;
                            removed_seedlst := mj.seedNO::(!removed_seedlst);
                            if debug2 then
                            Printf.printf "add seed#%d to removelst2\n%!" mj.seedNO;
                            end;
                    end
                    else if (d>0) then begin
                        let new_jright = ileft-1 in
                        let new_jweight = new_jright-jleft+1 in
                        let weight_reduce = jright - ileft + 1 in
                        if debug2 then Printf.printf "d>0\n%!";
                        if (new_jweight<=0) then
                            Printf.printf "ileft=%d,newj_right=%d,jleft=%d\n%!"
                            ileft new_jright jleft;
                        assert(new_jweight>0);
                        if (new_jweight>min_weight) then begin
                            let continue_sign = ref true in
                            List.iter (fun record ->
                            if !continue_sign then begin
                            let trim_from_left = 
                                if (record.orientation=jori) then false else true 
                            in
                            let oldpos,newpos,pos2 =
                                get_three_pos trim_from_left record weight_reduce in
                            let new_mi = 
                                get_new_mi trim_from_left record newpos in
                            let new_mj = get_mum_from_mumtbl j_seedNO
                                mum_tbl seed2pos_tbl in
                            let current_klst,removej =
                            update_table_trim_right_or_left trim_from_left  
                            new_mj record new_mi oldpos newpos pos2 
                            new_jweight weight_reduce min_weight !acc_k_lst
                            in
                            acc_k_lst := current_klst;
                            if removej then begin
                                continue_sign := false;
                                removed_seedlst := mj.seedNO::(!removed_seedlst);
                            end;
                            if debug2 then
                            Printf.printf "add seed#%d to removelst3\n%!" mj.seedNO;
                            end
                            ) poslst_j;
                        end
                        else begin
                            remove_seed2 mj pos2seed_tbl_left
                            pos2seed_tbl_right mum_tbl seed2pos_tbl true;
                            removed_seedlst := mj.seedNO::(!removed_seedlst);
                            if debug2 then
                            Printf.printf "add seed#%d to removelst4\n%!" mj.seedNO;
                        end;
                    end; (* if d<0 else if d>0..*)
                    (*we only upate neighborhood when trim happens*)
                    if debug2 then Printf.printf "trim happends,update neighbor:%!";
                    Hashtbl.iter (fun key record ->
                      update_neighborhood key mum_tbl seed2pos_tbl 
                      pos2seed_tbl_left pos2seed_tbl_right false
                    ) seed2pos_tbl;
                    if debug2 then Printf.printf "update neighbor done\n%!"
                end; (* if trimable -- mi and mj are overlaped*)
            end (*if (mj and mi.extendable=0)*)
            else
                idx := !idx +1;
            (*mi's neighborhood list could be diferent due to trimming, we need
            * to load mi again*)
            let mi =
                if (List.mem i_seedNO !removed_seedlst) then 
                    {mi with neighborhood_lst=[]}
                else
                    get_mum_from_mumtbl i_seedNO mum_tbl seed2pos_tbl 
            in
            if (!idx > ((List.length mi.neighborhood_lst)-1)) then 
                sign := false
            else sign := true;
        done; (*while sign*)
        end
        else () (*mi.extendable=false*)
        end (*if i_seedNO is not in removed_seedlst *)
    ) seedlst;
    let res = ref 0 in
    let count_f mum = 
        if (mum.extendable = 0) then res := !res + 1 in
    Hashtbl.iter (fun key bt ->
        iter_b_tree bt false count_f
    ) mum_tbl;
    if debug then 
        Printf.printf "end of resolve overlap mums, extendable mum = %d\n %!" !res;
    !res
    
let update_score_for_each_mum mum_tbl in_seqarr = 
    let update_mum_score record = 
        if (record.extendable = 0) then begin
        let score = get_mum_score in_seqarr record in
        let sign_newmum = update_mum_to_mumtbl None 
        {record with mumscore = score } mum_tbl false in
        assert(sign_newmum);
        end
    in
    Hashtbl.iter (fun key bt ->
        iter_b_tree bt false update_mum_score 
    ) mum_tbl


let search_inside_each_lcb lcb_tbl in_seqarr in_seq_size_lst min_len max_len = 
    Printf.printf "Search inside each lcb ,min and max len = %d/%d\n%!" min_len max_len;
    let new_seq2seedNO_tbl = Hashtbl.create init_tb_size in
    let new_seedNO2seq_tbl = Hashtbl.create init_tb_size in
    let new_pos2seed_tbl_left = Hashtbl.create init_tb_size in
    let new_pos2seed_tbl_right = Hashtbl.create init_tb_size in
    let new_seed2pos_tbl = Hashtbl.create init_tb_size in
    let new_mum_tbl = Hashtbl.create init_tb_size in
    Hashtbl.iter (fun key lcbrecord ->
        Printf.printf "work on lcb:%!";
        print_lcb lcbrecord;
        let lcblen = lcbrecord.avg_range_len in        
        let rangelst = List.sort (fun x y -> compare x.sequence_NO
        y.sequence_NO) lcbrecord.range_lst in
        let lcb_seqarr = 
            Array_ops.map_2 (fun seq mi->
                  get_sub_seq2 seq mi.left_end (mi.right_end-mi.left_end+1) 
            ) in_seqarr (Array.of_list rangelst)
        in
        if (lcblen>min_len)&&(lcblen<max_len) then begin
            Hashtbl.clear new_seq2seedNO_tbl;
            Hashtbl.clear new_seedNO2seq_tbl;
            Hashtbl.clear new_pos2seed_tbl_left;
            Hashtbl.clear new_pos2seed_tbl_right;
            Hashtbl.clear new_seed2pos_tbl;
            Hashtbl.clear new_mum_tbl;
            let avglen = get_avg_of_intlst in_seq_size_lst in
            let seedlen = int_of_float ( ceil (log avglen)) in
            let seedlen = if (seedlen mod 2)=0 then seedlen+1 else seedlen in 
            let seedlen = if (seedlen<5) then 5 else seedlen in
            let seedlen = if (seedlen=17) then 19 else seedlen in
            let seedweight = build_seed_and_position_tbl lcb_seqarr seedlen 
            !minimum_cover_ratio  new_seq2seedNO_tbl new_seedNO2seq_tbl 
            new_pos2seed_tbl_left new_pos2seed_tbl_right 
            new_seed2pos_tbl new_mum_tbl false in
            Printf.printf "seedweight = %d\n%!" seedweight;
            build_local_mums2 new_mum_tbl new_seed2pos_tbl 
            new_pos2seed_tbl_left new_pos2seed_tbl_right  false false;
            let num_of_mums = resolve_overlap_mum new_mum_tbl seedweight 
            new_pos2seed_tbl_left new_pos2seed_tbl_right new_seed2pos_tbl false false in
            Printf.printf "num of mums = %d\n%!" num_of_mums;
        end
        else
            Printf.printf "this lcb is too large or too small to work with\n%!";
        
    ) lcb_tbl;
    Printf.printf "Search inside each lcb, done\n%!"


(* search area outside existing lcb block. *)    
let search_outside_lcbs inner_lcbs lcb_tbl mum_tbl
pos2seed_tbl_left pos2seed_tbl_right seed2pos_tbl
in_seqarr in_seq_size_lst  =
    let debug = false and debug2 = false in
    if debug then begin 
        Printf.printf "\n search outside lcbs, lcb tbl size = %d\n%!" 
        (Hashtbl.length lcb_tbl);
        if debug2 then Hashtbl.iter (fun key record -> print_lcb record ) lcb_tbl;
    end;
    (*this works even when the seq outside is shorter than seed's length*)
    let lcb_range_lstlst, seq_outside_lcbs_arr = 
        get_seq_outside_lcbs in_seqarr in_seq_size_lst lcb_tbl in
    let current_seq_size_lst = Array.fold_right (fun seq acc ->
        (Array.length seq)::acc ) seq_outside_lcbs_arr [] 
    in
    if debug then begin
        Printf.printf "current_seq_size_lst :%!";
        print_int_list current_seq_size_lst;
    end;
    let shorted_seqlen = get_shorter_len seq_outside_lcbs_arr in
    let new_seedlen = int_of_float ( ceil (log (float shorted_seqlen))) in
    let new_seedlen = if (new_seedlen mod 2)=0 then new_seedlen+1 else new_seedlen in 
    let new_seedlen = 
        if new_seedlen=17 then 15 else new_seedlen in
    let new_seedlen = if (new_seedlen<5) then 5 else new_seedlen in
    if debug then 
        Printf.printf "seedlen=%d,%!" new_seedlen;
    (*these hashtable are based on the sub-seq list*)
    let new_seq2seedNO_tbl = Hashtbl.create init_tb_size in
    let new_seedNO2seq_tbl = Hashtbl.create init_tb_size in
    let new_pos2seed_tbl_left = Hashtbl.create init_tb_size in
    let new_pos2seed_tbl_right = Hashtbl.create init_tb_size in
    let new_seed2pos_tbl = Hashtbl.create init_tb_size in
    let new_mum_tbl = Hashtbl.create init_tb_size in
    let new_seedweight = 
        build_seed_and_position_tbl seq_outside_lcbs_arr new_seedlen 
        !minimum_cover_ratio new_seq2seedNO_tbl new_seedNO2seq_tbl 
        new_pos2seed_tbl_left new_pos2seed_tbl_right 
        new_seed2pos_tbl new_mum_tbl false in
    if debug then 
        Printf.printf "base seedweight=%d\n%!" new_seedweight;
    build_local_mums2 new_mum_tbl new_seed2pos_tbl 
    new_pos2seed_tbl_left new_pos2seed_tbl_right false false;
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
        Printf.printf "\n ================ Extend seeds: ================\n%!";
    extend_seeds res_mum_tbl res_seed2pos_tbl res_pos2seed_tbl_left res_pos2seed_tbl_right;
    if debug then 
        Printf.printf "\n ============= Resolve overlap mum: =============\n%!";
    let num_of_mums = resolve_overlap_mum res_mum_tbl (new_seedweight-1)
    res_pos2seed_tbl_left res_pos2seed_tbl_right res_seed2pos_tbl false false in
    update_score_for_each_mum res_mum_tbl in_seqarr;
    (*get the new weight of lcbs*)
    let new_seedNOlstlst = 
        get_mum_lst_for_each_seq res_mum_tbl res_seed2pos_tbl res_pos2seed_tbl_left
        (Array.length in_seqarr) in_seq_size_lst in
    if debug then begin
        Printf.printf "new seedNOlstlst is :\n%!";
        print_int_lstlst3  new_seedNOlstlst;
    end;
    let new_lcbs,new_covR,num_badlcb,new_lcb_tbl = 
        build_LCBs new_seedNOlstlst res_mum_tbl res_seed2pos_tbl 
         in_seq_size_lst false in
    if debug then 
        Printf.printf "new covR = %f, end of searching outside lcbs,check tbl \
        size (%d,%d,%d,%d)\n%!" new_covR (Hashtbl.length res_mum_tbl)
        (Hashtbl.length res_seed2pos_tbl) (Hashtbl.length res_pos2seed_tbl_left)
        (Hashtbl.length res_pos2seed_tbl_right);
    new_lcbs,new_covR,new_lcb_tbl,res_mum_tbl, 
    res_pos2seed_tbl_left, res_pos2seed_tbl_right, res_seed2pos_tbl, num_of_mums

let rec get_init_lcbs seedNOlstlst 
seed2pos_tbl mum_tbl 
in_seq_size_lst init_num_mums 
min_lcb_ratio min_cover_ratio previous_fullcovR=
    let debug = true in
    if debug then Printf.printf "get initial lcbs with min_lcbR=%f,min_covR=%f\n%!"
    min_lcb_ratio min_cover_ratio;
    assert(min_lcb_ratio > 0. );
    minimum_cover_ratio := min_cover_ratio ; 
    minimum_lcb_ratio := min_lcb_ratio ;
    let init_lcbs,init_covR,num_badlcb,init_lcb_tbl = 
        build_LCBs seedNOlstlst mum_tbl seed2pos_tbl in_seq_size_lst
        false in
    let init_lcbs, init_covR, init_lcb_tbl = 
        remove_bad_lcbs init_lcbs init_lcb_tbl mum_tbl 
        seed2pos_tbl in_seq_size_lst init_num_mums in
    if debug then begin
        Printf.printf 
        "after remove light lcbs, we have (don't show lightW or lowR lcbs) :\n%!";
        Hashtbl.iter (fun key record ->
            if (is_light_weight_lcb record in_seq_size_lst) then ()
            else
                print_lcb record 
        ) init_lcb_tbl;
        Printf.printf "init lcb covR = %f\n%! " init_covR;
    end;
    if init_covR < !minimum_full_cover_ratio then begin
        let new_lcbR,new_covR = get_adjusted_parameter init_covR in_seq_size_lst
        min_lcb_ratio min_cover_ratio in
        if (tiny_improvement previous_fullcovR init_covR)||(new_lcbR <
        !stop_lcb_ratio) then begin
            if debug then 
                Printf.printf "Warning: not much improve (%f<=>%f), or lcb ratio
                \ drops too low(%f<=>%f), just return initial lcbs.\n%!" 
                previous_fullcovR init_covR new_lcbR !stop_lcb_ratio;
            init_lcbs,init_covR,init_lcb_tbl
        end
        else
            get_init_lcbs seedNOlstlst 
            seed2pos_tbl mum_tbl 
            in_seq_size_lst init_num_mums
            new_lcbR new_covR init_covR;
    end
    else 
        init_lcbs,init_covR,init_lcb_tbl
    

    
let create_lcb_tbl in_seqarr min_lcb_ratio min_cover_ratio min_bk_penalty =
    min_break_point_penalty := min_bk_penalty ;
    let debug = false and debug2 = false in
    seedNO_available_arr := Array.make init_seed_size 1;
    (*output result to file ...    
    * let outfile = "outfile.txt" in let oc = open_out outfile in*)
    let in_seqarr_size = 2 in
    let in_seq_size_lst = Array.fold_right 
    (fun seq acc -> (Array.length seq)::acc ) in_seqarr [] 
    in
    let shorted_seqlen = get_shorter_len in_seqarr in
    let seedlen = int_of_float ( ceil (log (float shorted_seqlen))) in
    let seedlen = if (seedlen mod 2)=0 then seedlen+1 else seedlen in 
    let seedlen = if (seedlen<5) then 5 else seedlen in
    let seedlen = if (seedlen=17) then 19 else seedlen in
    if debug then 
        Printf.printf "\n **********Block_mauve.create_lcb_tbl;seqlen=%d,%d;seedlen=%d \
        ****************\n%!" (List.hd in_seq_size_lst) (List.nth
        in_seq_size_lst 1 ) seedlen;
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
    let mum_tbl = Hashtbl.create init_tb_size in
    (*find initial mums*)
    let seedweight = build_seed_and_position_tbl in_seqarr seedlen 
    !minimum_cover_ratio seq2seedNO_tbl seedNO2seq_tbl pos2seed_tbl_left
    pos2seed_tbl_right seed2pos_tbl mum_tbl false in
    if debug then 
        Printf.printf "base seedweight=%d, call build_local_mums2 \n%!" seedweight;
    build_local_mums2 mum_tbl seed2pos_tbl 
    pos2seed_tbl_left pos2seed_tbl_right false false;
    if debug then
        Printf.printf "++++++++ init seedweight=%d, end of building mum\
        table\n%!" seedweight;
    let init_num_mums = resolve_overlap_mum mum_tbl (seedweight-1) 
    pos2seed_tbl_left pos2seed_tbl_right seed2pos_tbl false false in
    update_score_for_each_mum mum_tbl in_seqarr;
    if debug then
        Printf.printf "++++++++++ resolve_overlap_mum done (tbl size\
        check:%d,%d,%d,%d) ++++ \n%!" (Hashtbl.length mum_tbl) (Hashtbl.length
        seed2pos_tbl) (Hashtbl.length pos2seed_tbl_left) (Hashtbl.length
        pos2seed_tbl_right);
    if debug2 then print_mumtbl mum_tbl true;
    let seedNOlstlst = get_mum_lst_for_each_seq mum_tbl seed2pos_tbl 
    pos2seed_tbl_left in_seqarr_size in_seq_size_lst in 
    if debug2 then begin
        Printf.printf "init seedNOlstlst is :\n%!";
        print_int_lstlst3  seedNOlstlst;
    end;
    (*get initial lcbs*)
    let init_lcbs,init_covR,init_lcb_tbl = 
        get_init_lcbs seedNOlstlst seed2pos_tbl mum_tbl in_seq_size_lst
        init_num_mums min_lcb_ratio min_cover_ratio 0. in
    (* sign of any improvement in outer&inner while loops*)
    let any_improvement = ref false in
    (*init inner tbl*)
    let inner_old_covR = ref !minimum_cover_ratio in (*this is not the init_covR
    we get from get_init_lcbs, that number is the covR after remove light weight
    lcbs, this one is not*)
    let inner_lcbs = ref init_lcbs in
    let inner_lcb_tbl = ref (Hashtbl.copy init_lcb_tbl) (* not necessary *)in
    let inner_mum_tbl = ref mum_tbl in
    let inner_pos2seed_tbl_left = ref pos2seed_tbl_left in
    let inner_pos2seed_tbl_right = ref pos2seed_tbl_right in
    let inner_seed2pos_tbl = ref seed2pos_tbl in
    (*init outer tbl*)
    let outer_old_covR = ref init_covR in
    let outer_lcbs = ref init_lcbs in
    let outer_lcb_tbl = ref (Hashtbl.copy init_lcb_tbl) (* necessary *)in
    let outer_mum_tbl = ref mum_tbl in
    let current_num_of_mums = ref init_num_mums in
    let outer_sign = 
        ref ((init_covR< !maximum_cover_ratio)&&(init_covR>= !minimum_cover_ratio)) in
    while (!outer_sign) do
        if debug then
        Printf.printf "\n begin of outer while, old_covR = %f\n%!" !outer_old_covR;
    (*****work on regions outside of lcbs******)
    (* we do need to update inner lcbs with the outer one, so we can work on new
    * range of subsequence, but inner_old_weight should remain the same, 
    * since inner_weight is the weight before any low-weight lcb removing *)
        inner_lcbs := !outer_lcbs;
        inner_mum_tbl := !outer_mum_tbl; (*this one is not necessary*)
        inner_lcb_tbl := Hashtbl.copy !outer_lcb_tbl;
        let inner_sign = ref true in
        while (!inner_sign) do
            if debug then
                Printf.printf "\n ----------- Inner loop -------------- \n%!";
            let inner_new_lcbs,inner_new_covR,res_lcb_tbl,res_mum_tbl,
            res_pos2seed_tbl_left, res_pos2seed_tbl_right, res_seed2pos_tbl, 
            num_of_mums = 
            search_outside_lcbs !inner_lcbs !inner_lcb_tbl !inner_mum_tbl 
            !inner_pos2seed_tbl_left !inner_pos2seed_tbl_right 
            !inner_seed2pos_tbl in_seqarr in_seq_size_lst 
            in
            if (num_of_mums>0)&&(!inner_old_covR < inner_new_covR) then begin
                if debug then
                    Printf.printf "we found better lcbs (%f>%f), update \
                    num_of_mums from %d to %d,continue with inner loop\n%!" 
                    inner_new_covR !inner_old_covR !current_num_of_mums num_of_mums;
                if (inner_new_covR > !maximum_cover_ratio)||(tiny_improvement
                inner_new_covR !inner_old_covR) 
                    then inner_sign :=false
                else inner_sign := true;
                current_num_of_mums := num_of_mums;
                inner_old_covR := inner_new_covR;
                inner_lcbs := inner_new_lcbs;
                inner_lcb_tbl := Hashtbl.copy res_lcb_tbl;
                inner_mum_tbl := res_mum_tbl;
                inner_pos2seed_tbl_left := res_pos2seed_tbl_left;
                inner_pos2seed_tbl_right := res_pos2seed_tbl_right;
                inner_seed2pos_tbl := res_seed2pos_tbl;
            end
            else begin
                inner_sign := false;
                if debug then Printf.printf "no improve on lcb (%f<=%f) \
                get out of inner loop\n%!" inner_new_covR !inner_old_covR ;
            end;
        done; 
        (*end of inner while loop*)
        if debug then
            Printf.printf "\n --------  Outer loop, remove bad lcbs ----- \n%!";
        let new_outer_lcbs, new_outer_covR, new_outer_lcb_tbl = 
            remove_bad_lcbs !inner_lcbs !inner_lcb_tbl !inner_mum_tbl 
            !inner_seed2pos_tbl in_seq_size_lst !current_num_of_mums in
        (*we could still have low W lcb here -- if remove them will create high W
        * low R lcbs, "remove_light_weight_lcbs" won't do it, get rid of those
        * lcbs here
        * Note. mauve keep the last low W lcb, we do the same, so
        * new_outer_weight will at least be 1*)
        Hashtbl.iter (fun key record ->
            if ((Hashtbl.length new_outer_lcb_tbl)>1)&&
            (is_light_weight_lcb record in_seq_size_lst) 
            then  
                Hashtbl.remove new_outer_lcb_tbl key;
        ) new_outer_lcb_tbl;
        (*let new_outer_covR = Hashtbl.length new_outer_lcb_tbl in*)
        let new_outer_covR,_ = get_lcb_covR_num_badlcb new_outer_lcb_tbl in_seq_size_lst in
        (* remove corresponding item in new_outer_lcbs*)
        let new_outer_lcbs = 
                List.map (fun lcblst ->
                List.filter (fun record ->
                 Hashtbl.mem new_outer_lcb_tbl (get_abs_lst record)
                ) lcblst
                ) new_outer_lcbs 
        in
        if debug&&(new_outer_covR=0.) then 
            Printf.printf "remove lcb did not give us any high W lcbs\n%!";
        if debug then begin
            Printf.printf "\n after remove light lcbs, we have :\n%!";
            Hashtbl.iter (fun key record ->
                print_lcb record 
            ) new_outer_lcb_tbl;
        end;
        if (!outer_old_covR < new_outer_covR) then begin
            if debug then 
                Printf.printf "\n we have a new lcb cov = %f>old one=%f, \
            continue with outer loop again\n%!" new_outer_covR !outer_old_covR;
            any_improvement := true;
            if (new_outer_covR> !maximum_cover_ratio)||(tiny_improvement
                new_outer_covR !outer_old_covR) then outer_sign := false
            else outer_sign := true;
            outer_old_covR := new_outer_covR;
            outer_lcbs := new_outer_lcbs;
            outer_lcb_tbl := Hashtbl.copy new_outer_lcb_tbl;
            outer_mum_tbl := !inner_mum_tbl;
        end
        else begin
            outer_sign := false;
            if debug then 
            Printf.printf "\n no improve on lcb weight (%f>=%f),\
            get out of outer loop\n%!" !outer_old_covR new_outer_covR;
        end
    done; (*end of outer while loop*)
    (*when outer&inner while did not find any qualified lcb, outer_lcb_tbl still
    * have the lcb from initial function, remove lightW ones*)
    let outer_lcb_tbl = !outer_lcb_tbl in
    if !any_improvement=false then 
        (*Hashtbl.clear outer_lcb_tbl; *)
        Hashtbl.iter (fun key record ->
            if (is_light_weight_lcb record in_seq_size_lst) then begin
                Hashtbl.remove outer_lcb_tbl key;
                outer_lcbs := List.map (fun lcblst ->
                List.filter (fun x -> (get_abs_lst x)<>key ) lcblst
            ) !outer_lcbs;
            end;
        ) outer_lcb_tbl;
    if debug then
        Printf.printf "init_covR=%f,outer old covR = %f, outer_lcb_tbl len=%d\n%!"
        init_covR !outer_old_covR (Hashtbl.length outer_lcb_tbl);
    if (Hashtbl.length outer_lcb_tbl)=1 then begin
        let lightW = ref false in
        Hashtbl.iter (fun key record ->
            if (is_light_weight_lcb record in_seq_size_lst)
                || (is_low_score_lcb record) then  lightW := true
        ) outer_lcb_tbl;
        if !lightW then  Hashtbl.clear outer_lcb_tbl;
    end;
    (*if no qualified lcb are found, just make the whole sequence as one lcb*)
    if (Hashtbl.length outer_lcb_tbl)=0 then begin 
        if debug then Printf.printf "we didn't find any qualified lcb\n%!";
        (*let lcb_range_lstlst = get_lcbs_range outer_lcb_tbl in_seqlst true in
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
                avg_range_len = 0;
            }
        in
        Hashtbl.clear outer_lcb_tbl;(*not necessary,outer_lcb_tbl should be empty here*)
        Hashtbl.add outer_lcb_tbl [1] single_lcb;
        let outer_lcbs = [[[1]];[[1]]] in
        outer_lcb_tbl,outer_lcbs,code_list,full_range_lstlst
    end
    else begin
        let outer_lcbs = !outer_lcbs in
        if debug then begin
            Printf.printf 
            "we do find some high W high R lcbs(or the last low RorW lcb), \
            return them to aliMap\n%!";
            Hashtbl.iter (fun key record -> print_lcb record) outer_lcb_tbl;
            Printf.printf "lcb list list is:%!";
            print_lcblst outer_lcbs;
        end;
        let lcb_range_lstlst = get_lcbs_range outer_lcb_tbl in_seq_size_lst in
        let full_range_lstlst = 
            get_full_range_lstlst lcb_range_lstlst in_seq_size_lst in
        if debug then begin
            Printf.printf "LCB range lstlst %!";
            print_lcbs_range lcb_range_lstlst;
            Printf.printf "the full range lstlst (includes non-LCB blocks)%!";
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
        if debug then Printf.printf "end of create lcb tbl\n%!";
        (*return lcb table, lcbs, code list and range list*)
        outer_lcb_tbl,outer_lcbs,code_list,full_range_lstlst
    end


(*main function here*)
(** [get_matcharr_and_costmatrix] take seq1 and seq2 as input sequence, devide 
genAli.create_gen_ali_new need codearr out of
[from_ori_code],not the original ones. while Utlgrappa.inv_med
requires the original codearr, call [to_ori_code] to transform 
* *)
let get_matcharr_and_costmatrix seq1 seq2 min_lcb_ratio min_cover_ratio min_bk_penalty 
locus_indel_cost cost_mat  =
    let debug = false and debug2 = false in
    let seq1arr = Sequence.to_array seq1 
    and seq2arr = Sequence.to_array seq2 in
    let in_seqarr = [|seq1arr;seq2arr|] in
    let lcb_tbl,lcbs,code_list,full_range_lstlst = 
        create_lcb_tbl in_seqarr min_lcb_ratio min_cover_ratio min_bk_penalty in
    if debug then begin
        Hashtbl.iter (fun key record ->
        print_lcb record ) lcb_tbl;
        (* Printf.printf "lcbs is \n%!";  Block_mauve.print_lcblst lcbs;*)
        Printf.printf "original code list is \n%!";
        print_int_lstlst code_list;
    end;
    let base = List.length (List.hd code_list) in (*start number of non-lcb block*)
    let len_lst1 = List.length (List.hd full_range_lstlst) in
    let len_lst2 = List.length (List.nth full_range_lstlst 1) in
    let gen_gap_code = (len_lst1 + len_lst2) * 2 + 1 in
    let matlen = gen_gap_code + 1 in
    if debug then
        Printf.printf "make empty matrix with size = %d,base=%d\n%!" matlen base;
    let gen_cost_mat = Array.make_matrix matlen matlen Utl.large_int in
    let set_cost code1 code2 cost = gen_cost_mat.(code1).(code2) <- cost in
    let empty_seq = Sequence.get_empty_seq () in
    let ali_mat = Array.make_matrix matlen matlen (0,empty_seq, empty_seq) in
    let block_gap_cost = Utl.large_int/len_lst1 in
    let seqNO = ref (-1) in
    let full_code_lstlst =
        List.map (fun full_range_lst ->
            seqNO := !seqNO +1 ;
            let coderef = ref base and start_num = !seqNO*len_lst1*2 in
            List.map(fun (left,right) ->
                let lcbkey,lcb_refcode = 
                    get_lcb_key_by_range !seqNO (left,right) lcb_tbl
                in
                if lcb_refcode<>0 then (*a lcb block, use lcb refcode*)
                    ((from_ori_code lcb_refcode)+start_num,(left,right))
                else begin (*not lcb block, give it a new code*)
                    coderef := !coderef +1;
                    ((from_ori_code !coderef)+start_num,(left,right))
                end
            )full_range_lst
        )full_range_lstlst 
    in
    if debug then begin
    Printf.printf "full (code,range) list is :\n%!";
    List.iter(fun full_code_lst ->
        List.iter (fun (code,(l,r)) ->
            Printf.printf "[%d,(%d,%d)],%!" code l r
        ) full_code_lst;
        print_newline();
    )full_code_lstlst;
    end;
    let edit_cost = ref 0 in
    let set_cost_code_arr code1_arr code2_arr gen_gap_code block_gap_cost =
        List.iter (fun (code1,(left1,right1)) ->
        List.iter (fun (code2,(left2,right2)) ->
            let subseq1,subseq2 = 
                    Sequence.sub seq1 left1 (right1-left1+1),
                    Sequence.sub seq2 left2 (right2-left2+1)
            in
            let ori_code1 = to_ori_code code1 
            and ori_code2 = to_ori_code (code2 - len_lst1*2) in
            if debug2 then 
                Printf.printf "work on %d,%d (ori=%d,%d) ->\n %!" 
                code1 code2 ori_code1 ori_code2;
            if ((abs ori_code2)<=base) && ((abs ori_code1)<=base) 
            && ((abs ori_code1)=(abs ori_code2)) then begin
                set_cost code1 code2 0;
                set_cost code2 code1 0; 
                (*Note: we only set cost(code1,code2) to 0 to make sure
                code1 is alied to code2 later in "GenAli.create_gen_ali_new",
                the real alignment cost is kept in ali_mat, then pass to
                function create_median_mauve in chromAli.ml.*)
                set_cost code1 gen_gap_code block_gap_cost; 
                set_cost gen_gap_code code1 block_gap_cost;
                set_cost code2 gen_gap_code block_gap_cost; 
                set_cost gen_gap_code code2 block_gap_cost;
                let alied_seq1, alied_seq2, cost, _ =  
                    Sequence.align2 subseq1 subseq2 cost_mat 
                in
                edit_cost := !edit_cost + cost; (*acc the real edit cost*)
                ali_mat.(code1).(code2) <- (cost,alied_seq1, alied_seq2);
                ali_mat.(code2).(code1) <- (cost,alied_seq2, alied_seq1);
                if debug2 then 
                    Printf.printf "set %d,%d with cost=%d,%d,(real ali cost=%d)\n%!"
                     code1 code2 0 block_gap_cost cost;
            end
            else if ((abs ori_code2)>base) && ((abs ori_code1)>base) then
                begin
                let alied_seq1, alied_seq2, cost, _ =  
                    Sequence.align2 subseq1 subseq2 cost_mat 
                in
                ali_mat.(code1).(code2) <- (cost,alied_seq1, alied_seq2);
                ali_mat.(code2).(code1) <- (cost,alied_seq2, alied_seq1);
                set_cost code1 code2 cost;
                set_cost code2 code1 cost;
                let del_cost1 = 
                    Sequence.cmp_gap_cost locus_indel_cost subseq1 
                and del_cost2 = 
                    Sequence.cmp_gap_cost locus_indel_cost subseq2
                in
                let len1 = Sequence.length subseq1 
                and len2 = Sequence.length subseq2 in
                let indel1 = Sequence.create_gap_seq len1  
                and indel2 = Sequence.create_gap_seq len2 in
                set_cost code1 gen_gap_code del_cost1;
                set_cost gen_gap_code code1 del_cost1;
                set_cost code2 gen_gap_code del_cost2;
                set_cost gen_gap_code code2 del_cost2;
                ali_mat.(code1).(gen_gap_code) <- (del_cost1,subseq1,indel1);
                (*cost in ali_mat.(gap).(code1) is for later looking up in this ali_mat*)
                ali_mat.(gen_gap_code).(code1) <- (del_cost1,indel1,subseq1);
                ali_mat.(code2).(gen_gap_code) <- (del_cost2,subseq2,indel2);
                ali_mat.(gen_gap_code).(code2) <- (del_cost2,indel2,subseq2);
                if debug2 then
                    Printf.printf "set %d,%d with cost=%d,del_cost=(%d|%d)\n%!" 
                    code1 code2 cost del_cost1 del_cost2;
            end
            else if ((abs ori_code2)>base) && ((abs ori_code1)<=base) then begin
                let del_cost2 = 
                    Sequence.cmp_gap_cost locus_indel_cost subseq2 in
                let len2 = Sequence.length subseq2 in
                let indel2 = Sequence.create_gap_seq len2 in
                ali_mat.(code2).(gen_gap_code) <- (del_cost2,subseq2,indel2);
                ali_mat.(gen_gap_code).(code2) <- (del_cost2,indel2,subseq2);
                set_cost code2 gen_gap_code del_cost2;
                if debug2 then 
                    Printf.printf "set %d,%d(gap) with cost=%d\n%!" code2
                    gen_gap_code del_cost2;
            end
            else if ((abs ori_code1)>base) && ((abs ori_code2)<=base) then begin
                let del_cost1 = 
                    Sequence.cmp_gap_cost locus_indel_cost subseq1 in
                let len1 = Sequence.length subseq1 in
                let indel1 = Sequence.create_gap_seq len1 in
                ali_mat.(code1).(gen_gap_code) <- (del_cost1,subseq1,indel1);
                ali_mat.(gen_gap_code).(code1) <- (del_cost1,indel1,subseq1);
                set_cost code1 gen_gap_code del_cost1;
                if debug2 then 
                    Printf.printf "set %d,%d(gap) with cost=%d\n%!" code1
                    gen_gap_code del_cost1;
            end
            else () (*we don't need to ali mauve-recognized block with other blocks*)
        ) (List.nth full_code_lstlst 1)
    ) (List.hd full_code_lstlst); 
    in
    let get_code_arr_from_fullcode_lst fullcode_lst = 
        Array.of_list ( List.map (fun (code,(_,_)) -> code) fullcode_lst )
    in
    let code1_arr = 
        get_code_arr_from_fullcode_lst (List.hd full_code_lstlst)
    and code2_arr =
        get_code_arr_from_fullcode_lst (List.nth full_code_lstlst 1)
    in
    set_cost_code_arr code1_arr code2_arr gen_gap_code block_gap_cost;
    code1_arr,code2_arr,gen_cost_mat,ali_mat,gen_gap_code,!edit_cost,full_code_lstlst,len_lst1

(** [get_range_with_code] return the range of match block code1 and block code2. if
    * a block is all_gap_seq, return (totalsize,totalsize+len_of_block-1) so that
    * graphic_mauve can work on it. totalsize starts as the size of the sequence,
    * increase each time when meet an all_gap block. *)
let get_range_with_code code1 code2 full_code_lstlst gen_gap_code totalsize1 totalsize2 =
    let debug = false in
    if debug then Printf.printf "get range with code %d,%d,gapcode=%d,total size1,2=%d,%d%!" 
    code1 code2 gen_gap_code totalsize1 totalsize2;
    let left1,right1 = 
        List.fold_left (fun acc (codex,(left,right)) ->
        if code1=codex then (left,right)
        else acc ) (-1,-1) (List.hd full_code_lstlst)
    in
    let left2,right2 = 
        List.fold_left (fun acc (codey,(left,right)) ->
        if code2=codey then (left,right)
        else acc ) (-1,-1) (List.nth full_code_lstlst 1)
    in
    if debug then Printf.printf "left1,right1=%d,%d,left2,right2=%d,%d\n%!" 
    left1 right1 left2 right2 ;
    if (code1=gen_gap_code)&&(code2<>gen_gap_code) then
        totalsize1,totalsize1,left2,right2
    else if (code2=gen_gap_code)&&(code1<>gen_gap_code) then
        left1,right1,totalsize2,totalsize2
    else if (code1<>gen_gap_code)&&(code2<>gen_gap_code) then
        left1,right1,left2,right2
    else
        failwith "Block_mauve.get_range_with_code, \
        we have indel match up with indel, something is wrong"

let get_seqlst_for_mauve in_seq =
    let totallen = Sequence.length in_seq in
    let iteration = totallen/80 in
    let reslst = ref [] in
    for i=0 to iteration do
        let startpos = i*80 in 
        let seqlen = 
            if (startpos+80)<totallen then 80
            else (totallen-startpos) in
        reslst := !reslst@[(Sequence.subseq in_seq startpos seqlen)]
    done;
    !reslst

    
let output2mauvefile filename cost old_cost alied_gen_seq1 alied_gen_seq2 full_code_lstlst
ali_mat gen_gap_code len_lst1 seqsize1 seqsize2 = 
let debug = false in
        (*let oc = open_out_gen [Open_creat(*;Open_append*)] 0o666 filename in
        let oc = open_out filename in*)
        let rewrite = match old_cost with
        | Some oldcost -> 
                if cost<oldcost then true else false
        | None -> true 
        in
        if rewrite then begin
            let oc = match filename with 
            | "" -> stdout
            | fname -> open_out fname in
            fprintf oc "#FormatVersion Mauve1\n";
            fprintf oc "#Sequence1File	somefile.in\n";
            fprintf oc "#Sequence1Entry	1\n";
            fprintf oc "#Sequence1Format	FastA\n";
            fprintf oc "#Sequence2File	somefile.in\n";
            fprintf oc "#Sequence2Entry	2\n";
            fprintf oc "#Sequence2Format	FastA\n";
            let totalsize1,totalsize2 = seqsize1, seqsize2 in
            List.iter2 ( fun alied_code1 alied_code2 ->
            let ori_code1 = to_ori_code alied_code1 
            and ori_code2 = to_ori_code (alied_code2 - len_lst1*2) in
            let cost,alied_seq1,alied_seq2 = 
                ali_mat.(alied_code1).(alied_code2)
            in
            let dir1 = if ori_code1>0 then "+" else "-"
            and dir2 = if ori_code2>0 then "+" else "-"
            in
            let left1,right1,left2,right2 = 
                get_range_with_code alied_code1 alied_code2 full_code_lstlst
                gen_gap_code totalsize1 totalsize2
            in
            if debug then Printf.printf 
            "output2mauve: left/right1=%d,%d, left/right2=%d,%d\n%!"
            left1 right1 left2 right2;
            (*let oc = open_out_gen [Open_creat(*;Open_append*)] 0o666 filename in*)
            let seqlst1 = get_seqlst_for_mauve alied_seq1 in
            let seqlst2 = get_seqlst_for_mauve alied_seq2 in
            fprintf oc "> 1:%d-%d %s c=%d\n" (left1+1) (right1+1) dir1 cost;
            if (right1-left1+1)<500 then 
                info_user_message "this sequence is too short for mauve graphic \
                output (length<500)";
            if left1=right1 then fprintf oc "X";
            List.iter (fun seq1 ->
                Sequence.print oc seq1 Alphabet.nucleotides;
                fprintf oc "\n";
            ) seqlst1;
            fprintf oc "> 2:%d-%d %s c=%d\n" (left2+1) (right2+1) dir2 cost;
            if left2=right2 then fprintf oc "X";
            List.iter (fun seq2 ->
                Sequence.print oc seq2 Alphabet.nucleotides;
                fprintf oc "\n";
            ) seqlst2;
            fprintf oc "=\n";
            ) (Array.to_list alied_gen_seq1) (Array.to_list alied_gen_seq2);
            if filename<>"" then close_out oc;
        end

