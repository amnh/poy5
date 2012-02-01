(* POY 5.0 Alpha. A phylogenetic analysis program using Dynamic Homologies.   *)
(* Copyright (C) 2011 Andrés Varón, Lin Hong, Nicholas Lucaroni, Ward Wheeler,*)
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

let () = SadmanOutput.register "Block_mauve" "$Revision: 2219 $"
(* A.D. = Aaron E. Darling*)
(* W = weight, R = ratio *)

open Printf
open Block_mauve_seed
let printIntArr = Utl.printIntArr
let printIntList = Utl.printIntList
let printIntList2 = Utl.printIntList
let printIntListListList = Utl.printIntListListList
let get_neg_rev_intlst = Utl.get_neg_rev_intlst
let find_index = Utl.find_index
let error_user_message format = Printf.ksprintf (Status.user_message Status.Error) format
let info_user_message format = Printf.ksprintf (Status.user_message Status.Information) format

let get_sub_seq2 = Array.sub

let init_tb_size = 50

let max_gap_num = 0  (*the w *)

(*[get_proper_seedlen inlen] return the proper seedlen.
* seedlen is the key to palindromic_spaced_seed_tbl, 
* seedlen cannot be bigger than 21 , or smaller than 5. 
* we only have entry for odd number, because even length of palidromic bring us
* problems, check out A.D.'s paper "Procrastination leads to efficient
* filtration for local multiple alignment".
* also, there is no entry for 17.*)
let get_proper_seedlen avg_seqlen =
    let avg_seqlen = float_of_int avg_seqlen in
    let seedlen = int_of_float ( ceil (log (avg_seqlen))) in
    let seedlen = if (seedlen mod 2)=0 then seedlen+1 else seedlen in 
    let seedlen = if (seedlen<5) then 5 else seedlen in
    let seedlen = if (seedlen=17) then 19 else seedlen in
    let seedlen = if (seedlen>21) then 21 else seedlen in
    seedlen


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

let compare_mums mum1 mum2 =
    if mum1.seedNO = mum2.seedNO then
        true
    else false
    
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


(*function to print mum*)
let print_neighborhood_list in_lst = 
    List.iter (fun (seqNO, j,i_ori,j_ori,dl,dr) ->
        Printf.printf "{seqNO.%d, seed#%d, ori:%d/%d,d:%d,%d }\n%!" seqNO j i_ori
        j_ori dl dr
    ) in_lst

let print_mum print_neighborhood print_unextendable in_mum =
    if (in_mum.extendable = 0) ||
    ((in_mum.extendable !=0 )&&(print_unextendable=true) ) then begin
        Printf.printf "MUM#%d\n%!" in_mum.seedNO;
        Printf.printf " mumkey : %d\n%!" in_mum.mumkey;
        Printf.printf " mumseq = %!";
        printIntArr in_mum.mumseq; 
        Printf.printf " size: %d ,\n" in_mum.size;
        Printf.printf " subsuming_pointer(is subsumed by) = %d ,%!" in_mum.subsuming_pointer;
        Printf.printf " extendable = %d \n%!" in_mum.extendable;
        Printf.printf " priority_lst = %!";
        printIntList2 in_mum.priority_lst;
        Printf.printf " positions: \n%!";
        print_pos_list in_mum.positions;
        if print_neighborhood then begin
        Printf.printf " neighborhood list : \n%!";
        print_neighborhood_list in_mum.neighborhood_lst;
        end;
    end


(**************** b_tree mum_tbl function starts *********************)

let print_mumtbl mum_tbl debug = 
    Hashtbl.iter (fun key record ->
        Printf.printf "mumkey = %d \n%!" key;
        BinaryTree.print_b_tree record (print_mum false true) printIntArr
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
    let res,newbt = BinaryTree.search_in_btree seedseq bt printIntArr 
    (print_mum false true) in
    Hashtbl.replace mum_tbl mumkey newbt;
    if res.seedNO<>seedNO then begin
        let _,seedseq1,w1 = Hashtbl.find seed2pos_tbl seedNO 
        and _,seedseq2,w2 = Hashtbl.find seed2pos_tbl res.seedNO in
        Printf.printf "seed1=%d,w=%d,%!" seedNO w1;
        printIntArr seedseq1;
        Printf.printf "; seed2=%d,w=%d,%!" res.seedNO w2;
        printIntArr seedseq2;
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
        assert(Array.length newmum.mumseq > 0);
        let new_bt,sign_newbt = BinaryTree.add_to_btree newmum.mumseq newmum 
        old_bt printIntArr (print_mum false true) compare_mums in
        if sign_newbt then Hashtbl.replace mum_tbl mumkey new_bt;
        if debug then Printf.printf "sign_newbt=%b\n%!" sign_newbt;
        sign_newbt
    end
    else begin
        let bt = BinaryTree.create_btree newmum.mumseq newmum in
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
    if (BinaryTree.just_a_leaf bt) then
            Hashtbl.remove mum_tbl mumkey
    else 
        let treekey = mum2remove.mumseq in
        let newbt = BinaryTree.remove_from_btree treekey bt printIntArr
        (print_mum false true) in
        if (BinaryTree.just_an_empty_leaf newbt) then
            Hashtbl.remove mum_tbl mumkey
        else Hashtbl.replace mum_tbl mumkey newbt
    


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

(***************** print function for seed2pos/pos2seed tbl *************)
let print_seedNO2seq_tbl in_tbl = 
    Hashtbl.iter (fun key record ->
        Printf.printf "seed#%d:%!" key;
        printIntList record;
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
        print_mum false true mum2print;
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
        printIntList2 seq;
        Printf.printf "}\n%!";
    ) in_tbl


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
    if debug then print_mum false true oldmum;
    assert(sign_newmum)


(** [update_priority_lst] add seedNO_to_remove to priority_lst, and remove
* seedNO_to_add to priority_lst. return if priority_lst is empty.
* same seedNO in priority_lst can show up more than one time.*)
let update_priority_lst in_mum (seedNO_to_remove:int option) (seedNO_to_add:int
option) mum_tbl seed2pos_tbl=
    let debug = false in
    if debug then begin
        Printf.printf "update prilst on seed#%d," in_mum.seedNO;
        print_mum false true in_mum;
    end;
    let old_lst = in_mum.priority_lst in
    let lst_after_remove = 
    match seedNO_to_remove with
    | Some seedNO->
        let debug = 
            if (List.mem seedNO old_lst)=false then true
            else debug
        in
        if debug then print_mum false true in_mum;
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

(*************** add seed to mum/seed2pos/pos2seed tbl function start  ***************)
(*[add_seed_to_tbl] add match we found during [find_SML] to
* mum_tbl/seed2pos/pos2seed tbl. expend the match to both direction before
* adding it. if the result of extension is already in those tbls, skip this one*)
let add_seed_to_tbl init_accarr init_seedweight inseqarr input_seqlst_size
position2seed_tbl_left position2seed_tbl_right seed2position_tbl mum_tbl seedNO_available_arr =
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
        let seedNO = get_a_seedNO seedNO_available_arr in
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
                return_a_seedNO seed2remove seedNO_available_arr;
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
                    mumalgn = [];
                } in
            if debug then begin
                Printf.printf "add seed %d with mumkey %d,treekey = %!" seedNO key;
                print_mum false true newmum;
                Printf.printf " to mum_tbl\n%!";
            end;
            let sign_newmum = add_mum_to_mumtbl newmum mum_tbl in
            assert(sign_newmum); 
        end
        else
            return_a_seedNO seedNO seedNO_available_arr;
    end


let find_SML patternarr init_seedweight  inseqarr input_seqlst_size
position2seed_tbl_left position2seed_tbl_right seed2position_tbl mum_tbl seedNO_available_arr = 
    let debug = false and debug2 = false in
    if debug then 
        Printf.printf "get sequences with index and possible seed
        subseq,seedweight=%d\n%!" init_seedweight;
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
        if debug then
            Printf.printf "seqlen=%d,init_seedweight=%d\n%!" seqlen init_seedweight;
        Array.sub resarr 0 (seqlen-init_seedweight+1)
    ) inseqarr
    in
    if debug then 
        Printf.printf "append it with reverse complement of first sequence\n%!";
    let rev_seq_w_idx = Array.map (fun (sequenceNO,pos,dir,subseq,_) ->
        let rev_subseq = Alphabet.rev_comp_arr subseq Alphabet.nucleotides in
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
            position2seed_tbl_left position2seed_tbl_right seed2position_tbl
            mum_tbl seedNO_available_arr;
        newacc
    ) [] sorted_arr in
    if (List.length last_acclst)>1 then 
        add_seed_to_tbl (Array.of_list last_acclst) init_seedweight inseqarr input_seqlst_size
position2seed_tbl_left position2seed_tbl_right seed2position_tbl mum_tbl seedNO_available_arr
        

(*try to deal with huge dataset, remember to use 64bits*)
let scan_seqlst2 (inseqarr:int array array) patternarr mum_tbl
position2seed_tbl_left position2seed_tbl_right seed2position_tbl seedNO_available_arr debug=
    let input_seqarr_size = 2 in (*now we only have 2 input sequences*)
    let init_seedweight = patternarr.((Array.length patternarr)-1) in
    if debug then Printf.printf "scan_seqlst2, seedweight=%d, find_SML with inseqlst_w_idx_w_rev%!" 
    init_seedweight;
    find_SML patternarr init_seedweight inseqarr input_seqarr_size 
    position2seed_tbl_left position2seed_tbl_right seed2position_tbl mum_tbl seedNO_available_arr; 
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



let build_seed_and_position_tbl inseqarr init_seedlen 
(*seq2seedNO_tbl*) seedNO2seq_tbl pos2seed_tbl_left pos2seed_tbl_right
seed2pos_tbl mum_tbl seedNO_available_arr debug
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
        printIntList (seedpattern); 
    end;
    let patternarr = Array.of_list seedpattern in
    let _ = scan_seqlst2 inseqarr patternarr  
        mum_tbl pos2seed_tbl_left pos2seed_tbl_right seed2pos_tbl seedNO_available_arr false in
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

(*************** add seed to mum/seed2pos/pos2seed  tbl function start  ***************)

(******************************* positions and seqNO **************************)

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
        if sign_newmum then print_mum false true new_mum;
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



(******************************* positions and seqNO function end **************************)

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
        BinaryTree.iter_b_tree bt format_mum printIntArr false
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


(************************* neighborhood function start ********************)

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
    if debug then print_mum true true i_mum;
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
                            if debug2 then print_mum false true j_mum;
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
                        if debug2 then print_mum false true j_mum;
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
            print_mum true true imum;
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

(************************* neighborhood function end ********************)

(************************* seed extension funnction start ********************)
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
            print_mum false true i_mum;
            print_mum false true j_mum;
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
        print_mum false true i_mum;
        print_mum false true j_mum;
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
                    print_mum false true i_mum;
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
                    if debug then print_mum true true i_seedNO mum_tbl;
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
                   print_mum true true i_seedNO mum_tbl;
                   Printf.printf " === end of while ===\n%!";
               end;  *)
            done; (* end of while ( !find_a_chain = 1 ) *)
        end;
    ) seed2pos_tbl


(************************* seed extension funnction start ********************)


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


(********************** add mum with new pos system to old mum tbl function start **********)
let to_ori_position position lcb_range_lstlst =
    let debug = false in
    let leftend = position.left_end and rightend = position.right_end in
    let len = rightend - leftend +1 in
    let seqNO = position.sequence_NO in
    if debug then
        Printf.printf "to_ori_position, seqNO=%d,l/r = %d/%d, %!" seqNO leftend rightend;
    let lcb_range_lst = List.nth lcb_range_lstlst seqNO in
    let (last_acc_size,last_right,ori_left,ori_right) = 
        List.fold_left (fun (acc_size,pre_r,res_left,res_right) (l,r,_,_) ->
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

(*transpose positions of new_mum, join it with existing mums.*)
let transpose_positions lcb_range_lstlst pos2seed_tbl_left pos2seed_tbl_right
seed2pos_tbl mum_tbl seedNO_available_arr new_mum  = 
let debug = false in
    if debug then Printf.printf "mum before transpose: %!";
    if debug then print_mum false true new_mum;
    if new_mum.extendable<>2 then begin
        let sign = ref true in (*true means we can add the new mum to mum_tbl*)
        let ori_positions = List.map (fun mi ->
            let ori_l,ori_r = to_ori_position mi lcb_range_lstlst in
            if (ori_l=(-2)) then sign:=false;
            {mi with left_end = ori_l; right_end = ori_r}
        ) new_mum.positions in
        if (!sign = true) then begin
                let seedNO = get_a_seedNO seedNO_available_arr in
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
                if debug then print_mum false true mum2add;
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
    return_a_seedNO new_mum.seedNO seedNO_available_arr


(*transpose mums found outside current ones, join them to current ones.*)
let transpose_mum_back lcb_range_lstlst new_mum_tbl new_seedNO2seq_tbl
mum_tbl pos2seed_tbl_left pos2seed_tbl_right seed2pos_tbl seedNO_available_arr
= 
    let debug = false in
    if debug then begin
        Printf.printf "transpose back, range list is :\n%!";
        List.iter (fun rangelst ->
            List.iter (fun (left,right,_,_) ->
                Printf.printf "(%d,%d),%!" left right
            ) rangelst;
            Printf.printf "\n%!";
        ) lcb_range_lstlst;
    end;
    let ourf = transpose_positions lcb_range_lstlst pos2seed_tbl_left 
    pos2seed_tbl_right seed2pos_tbl mum_tbl seedNO_available_arr in
    Hashtbl.iter (fun mumkey bt ->
        BinaryTree.iter_b_tree bt ourf printIntArr false
    ) new_mum_tbl;
    Hashtbl.iter (fun seedNO record ->
        let debug_nei = false in
        update_neighborhood seedNO mum_tbl seed2pos_tbl 
        pos2seed_tbl_left pos2seed_tbl_right debug_nei;
    ) seed2pos_tbl;
    if debug then Printf.printf "end of transpose back\n%!"

(********************** add mum with new pos system to old mum tbl function end **********)

    
(*[get_mum_lst_for_each_seq] returns the lcb_lst_lst we need to build lcb blocks 
* to do: rewrite this with seed2pos tbl *)
let get_mum_lst_for_each_seq mum_tbl seed2pos_tbl pos2seed_tbl seqlst_size seq_size_lst 
= 
    let debug = false in
    if debug then Printf.printf "get lcb lst for each seq\n%!";
    let seedNOlstlst = ref [[]] in
    for i = 0 to (seqlst_size-1) do
        if debug then Printf.printf "work on sequence %d\n%!" i;
        let reslst = ref [[]] in
        for j=0 to (List.nth seq_size_lst i)-1 do
            if (Hashtbl.mem pos2seed_tbl (i,j)) then begin
                let (seedNO,seed_weight,ori) = 
                    get_extendable_record pos2seed_tbl (i,j)
                    seed2pos_tbl mum_tbl in
                if (seedNO<>(-1)) then
                    begin
                    if debug then
                        Printf.printf "seed#%d at pos (%d,%d)\n%!" seedNO i j;
                        reslst := [ori*seedNO]::!reslst
                    end
                else ();
            end
            else ();
        done;
        seedNOlstlst := !reslst::!seedNOlstlst
    done;
    let seedNOlstlst = List.filter (fun item -> (List.length item)>0) 
    (List.rev !seedNOlstlst) in
    if debug then begin (*sanity check*)
    Printf.printf "seedNOlst for each seq is :\n%!"; 
    printIntListListList seedNOlstlst;
    let seedNOlst1 = List.hd seedNOlstlst 
    and seedNOlst2 = List.nth seedNOlstlst 1 in
    List.iter (fun seedNO ->
        let negseedNO = get_neg_rev_intlst seedNO in
        if (List.mem seedNO seedNOlst1 )||(List.mem negseedNO seedNOlst1 ) 
        then ()
        else begin 
             assert((List.length seedNO)=1);
            let seedNO = List.hd seedNO in
            Printf.printf "cannot find seedNO.#%d \n%!" seedNO;
            let badmum = get_mum_from_mumtbl seedNO mum_tbl seed2pos_tbl in
            print_mum false true badmum;
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
            print_mum false true mum;
        let score = get_score_from_2seq previous_seq current_seq ori in
        sum_score + score, current_seq, current_ori
    ) (0,first_seq,first_ori) (List.tl poslst)
    in
    score

(********************* resolve overlap mum function start **************)

 
let update_k_seed_lst pre_jseedNO previous_k_lst pre_kposlst k_pos j_mumseq trim_from_left
weight_reduce pos2seed_tbl1 pos2seed_tbl2 seed2pos_tbl mum_tbl seedNO_available_arr =
let debug = false in
    let k_seedNO = get_a_seedNO seedNO_available_arr in
    let k_extsign = 0 in (*since mi and mj are both ext=0*)
    let k_milst = k_pos::pre_kposlst in
    let k_size = (List.length k_milst) in
    if debug then 
        Printf.printf "we have a new seedk : %d,W=%d\n%!" k_seedNO weight_reduce;
    let k_mumkey = get_mumkey_from_milst k_milst in
    if debug then begin
        Printf.printf "kmumkey = %d\n%!" k_mumkey;
        printIntArr j_mumseq;
    end;
    let start_pos = 
        if trim_from_left then 0 
        else (Array.length j_mumseq)-weight_reduce
    in
    let k_mumseq = get_sub_seq2 j_mumseq start_pos weight_reduce in  
    if debug then printIntArr k_mumseq;
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
        mumscore = 0;
        mumalgn = []} in  
    if debug then print_mum false true new_kmum;
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
        add_to_pos2seed_tbl true pos2seed_tbl_l pos2seed_tbl_r 
        seqNO left k_seedNO klen 
        ori k_extsign k_size seed2pos_tbl mum_tbl;
        add_to_pos2seed_tbl false pos2seed_tbl_l pos2seed_tbl_r 
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
previous_k_lst seedNO_available_arr = 
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
                pos2seed_tbl1 pos2seed_tbl2 seed2pos_tbl mum_tbl seedNO_available_arr
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
pos2seed_tbl_right seed2pos_tbl seedNO_available_arr debug debug_neighborhood =
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
    oldpos newpos pos2 new_jweight weight_reduce min_weight previous_k_lst seedNO_available_arr =
        let record = current_mi in
        if (trim_from_left) then
            update_tables mum_tbl mumj record.sequence_NO pos2seed_tbl_left 
            oldpos newpos pos2seed_tbl_right pos2 new_jweight seed2pos_tbl 
            new_mi record.orientation weight_reduce min_weight true previous_k_lst 
            seedNO_available_arr
        else 
            update_tables mum_tbl mumj record.sequence_NO pos2seed_tbl_right 
            oldpos newpos pos2seed_tbl_left pos2 new_jweight seed2pos_tbl 
            new_mi record.orientation weight_reduce min_weight false previous_k_lst
            seedNO_available_arr
    in
    let seedlst = ref [] and removed_seedlst = ref [] in 
    Hashtbl.iter (fun key record -> seedlst := key::!seedlst ) seed2pos_tbl;
    let seedlst = List.rev !seedlst in
    List.iter (fun i_seedNO -> 
        if (List.mem i_seedNO !removed_seedlst)=false then begin
        let debug2 = false in
        if debug2 then Printf.printf "resolve overlap,work on seed#%d\n%!" i_seedNO;
        let mi = get_mum_from_mumtbl i_seedNO mum_tbl seed2pos_tbl in
        if debug2 then print_mum true true mi;
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
                Printf.printf "idx = %d, i_mum is: %!" !idx; print_mum true false mi;
            end;
            let (seqNO,j_seedNO,i_ori,j_ori,d,dr) = List.nth mi.neighborhood_lst !idx in
            if debug2 then Printf.printf "i=%d,j=%d,seqNO=%d\n%!" i_seedNO j_seedNO seqNO;
            let mj = get_mum_from_mumtbl j_seedNO mum_tbl seed2pos_tbl in
            if debug2 then begin 
                Printf.printf "j_mum is :%!"; print_mum false true mj;
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
                    if d=0 then print_mum false true mi;
                    if d=0 then print_mum false true mj;
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
                                weight_reduce min_weight !acc_k_lst seedNO_available_arr in
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
                            seedNO_available_arr 
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
        BinaryTree.iter_b_tree bt count_f printIntArr false
    ) mum_tbl;
    if debug then 
        Printf.printf "end of resolve overlap mums, extendable mum = %d\n %!" !res;
    !res
   

    
(********************* resolve overlap mum function end **************)

    
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
        BinaryTree.iter_b_tree bt update_mum_score printIntArr false
    ) mum_tbl


