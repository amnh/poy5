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

let () = SadmanOutput.register "Block_mauve" "$Revision: 2823 $"

(* A.D. = Aaron E. Darling*)
(* W = weight, R = ratio *)

let debug_main = false 

let debug_major_slice = false

open Printf
open Block_mauve_seed
open Block_mauve_mum
open Block_mauve_lcb
let printIntArr = Utl.printIntArr
let printIntMatWithIdx = Utl.printIntMatWithIdx
let printIntMat = Utl.printIntMatWithIdx
let printIntList = Utl.printIntList
let printIntList2 = Utl.printIntList2
let printIntListList = Utl.printIntListList
let printIntListListList = Utl.printIntListListList
let printIntListToFile = Utl.printIntListToFile
let get_neg_rev_intlst = Utl.get_neg_rev_intlst
let get_avg_of_intlst = Utl.get_avg_of_intlst
let get_abs_lst in_lst = List.sort compare (Utl.get_abs_intlst in_lst)

let error_user_message format = Printf.ksprintf (Status.user_message Status.Error) format
let info_user_message format = Printf.ksprintf (Status.user_message Status.Information) format


let skip_huge_nonlcb_block = true




(*lcb with lower weight should be considered removing*)
let minimum_lcb_weight = 2 
let init_tb_size = 50
let init_seed_size = 50

let checkpoint_major_slice msg =
    prerr_endline ("block_mauve.ml checkpoint_major_slice at position "^ msg);
    let slice_number = Gc.major_slice 0 in
    prerr_endline ("number of major slice = "^(string_of_int slice_number))


(** some operation on int/int list/int array/etc.. we gonna need *)
let to_ori_code code =
    assert(code<>0);
    if( code mod 2 == 0) then -(code/2)  else (code+1)/2 

let from_ori_code code =
    assert(code<>0);
    if (code<0) then (abs code)*2
    else (code*2-1)
    





let print_int_list2 inlist  = 
    List.iter (fun x -> Printf.printf "%3i%!" x) inlist;
    Printf.printf " \n%!"

let print_int_list3 inlist = 
    let idx = ref 0 in
    List.iter (fun x -> Printf.printf "[%d]:%d; " !idx x; idx:= !idx +1; ) inlist;
    Printf.printf " \n%!"

let print_int_lstlst3 inlstlst =
    List.iter (fun lst -> print_int_list3 lst) inlstlst

let get_sub_seq2 = Array.sub 

(*function for parsing mauve result file, just for comparing result of our mauve
* and original mauve*)

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


(*this is a test to replace the old way we score lcb
no one calls this now
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

*)

    


(* no longer useful [get_lcb_key_by_range] get the lcb_key and lcb_code by its range. if no lcb is found,
* return [],0. 
* remember what we did in [fill_in_indel]. the range in full_range_lstlst could
* be different from the original lcb range. 
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
    *)


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



(*[fill_in_cost_ali_mat] fill in cost between blocks .
* for cost matrix ali_mat. [Sequence.align2] is being called to
get the editing cost between two lcb blocks. [Sequence.cmp_gap_cost] is being
called to get the deletion cost of a lcb block*)
let fill_in_cost_ali_mat maximum_lcb_len cost_mat seq1 seq2 in_seqarr code_range_lst1 code_range_lst2 gen_gap_code
locus_indel_cost ali_mat gen_cost_mat base len_lst1 use_ukk 
mum_tbl seed2pos_tbl lcb_tbl =
    let debug = false and debug2 = false in
    if debug then Printf.printf "fill_in_cost_ali_mat with basecode=%d\n%!"
    base;
    let set_cost code1 code2 cost = gen_cost_mat.(code1).(code2) <- cost in
    let add_indel_cost code seq left right acc_cost =
        (*we are not align non-lcb block any more, they become indel blocks*)
        if debug then 
            Printf.printf "add indel cost to non-lcb block,code=%d,(left=%d,right=%d)," code left right;
        if (right>left) then begin
        let subseq = Sequence.sub seq left (right-left+1) in
        let del_cost = 
            Sequence.cmp_gap_cost locus_indel_cost subseq in
        if debug then 
            Printf.printf "gen_gap_code=%d,add indel_cost=%d \n%!" 
            gen_gap_code del_cost;
        let len = Sequence.length subseq in
        let indel = Sequence.create_gap_seq len in
        set_cost code gen_gap_code del_cost;
        set_cost gen_gap_code code del_cost;
        ali_mat.(code).(gen_gap_code) <- (del_cost,subseq,indel);
        (*cost in ali_mat.(gap).(code1) is for later looking up in this ali_mat*)
        ali_mat.(gen_gap_code).(code) <- (del_cost,indel,subseq);
        acc_cost + del_cost
        end
        else begin
            if debug then Printf.printf "empty non-lcb block,set cost to 0\n%!";
            set_cost code gen_gap_code 0;
            set_cost gen_gap_code code 0;
            acc_cost
        end
    in
    (*acc cost for edit and indel*)
    let edit_cost = ref 0 and indel_cost = ref 0 in
    (*collect edit cost of lcb blocks between seq1 and seq2, also indel cost from non-lcb block of seq1*)
    List.iter (fun (code1,(left1,right1),lcbkey1) ->
        if code1<=base*2 then begin
            let ori_code1 = to_ori_code code1 in
            let code2,(left2,right2),lcbkey2 = 
                try List.find (fun (c2,(_,_),_) ->
                    abs (to_ori_code (c2-len_lst1*2)) = abs ori_code1
                ) code_range_lst2 
                with | Not_found -> failwith "could not find matching lcb in second seq"
            in
            let avglen = get_avg_of_intlst [right1-left1+1;right2-left2+1] in
            let avglen = int_of_float avglen in
            set_cost code1 code2 0;
            set_cost code2 code1 0; 
            (*Note: we only set cost(code1,code2) to 0 to make sure
            code1 is alied to code2 later in "GenAli.create_gen_ali_new",
            the real alignment cost is kept in ali_mat, then pass to
            function create_median_mauve in chromAli.ml.*)
            let subseq1,subseq2 = 
                Sequence.sub seq1 left1 (right1-left1+1),
                Sequence.sub seq2 left2 (right2-left2+1)
            in
            let alied_seq1, alied_seq2, cost, _ =
            if (avglen > maximum_lcb_len) then begin
                if debug then Printf.printf "search inside a huge lcb block\n%!";
                let lcbrecord = Hashtbl.find lcb_tbl (get_abs_lst lcbkey1) in    
                let alied_seqarr,cost =
                    search_inside_a_lcb lcbrecord seq1 seq2 
                    maximum_lcb_len max_int mum_tbl seed2pos_tbl cost_mat use_ukk
                in
                if debug2 then Printf.printf "cost from search inside a lcb = %d\n%!" cost;
                if cost>0 then
                    alied_seqarr.(0),alied_seqarr.(1),cost,-1
                else
                    Sequence.align2 subseq1 subseq2 cost_mat use_ukk
            end
            else
                Sequence.align2 subseq1 subseq2 cost_mat use_ukk
            in
            edit_cost := !edit_cost + cost; (*acc the real edit cost*)
            ali_mat.(code1).(code2) <- (cost,alied_seq1, alied_seq2);
            ali_mat.(code2).(code1) <- (cost,alied_seq2, alied_seq1);
            if debug then 
                Printf.printf "set %d,%d with cost=%d,(real ali cost=%d)\n%!"
                 code1 code2 0 cost;
        end
        else 
            (*we are not align non-lcb block any more, they are indel blocks*)
            indel_cost := add_indel_cost code1 seq1 left1 right1 !indel_cost;
        ) code_range_lst1;
    (*work on non-lcb block of seq2*)
    List.iter (fun (code2,(left2,right2),lcbkey2) ->
        if (code2-len_lst1*2) <= (base*2) then ()
        else 
            indel_cost := add_indel_cost code2 seq2 left2 right2 !indel_cost;
    ) code_range_lst2;
    !edit_cost,!indel_cost

    
(*main function here*)
(** [get_matcharr_and_costmatrix] is the main function of this module. it take two
* sequences, set of parameters (min lcb ratio, min lcb length, etc) and 
* a costmatrix between charactors as input,
* output two code array of matching blocks, a cost_mat between lcb blocks and non-lcb
* blocks, ali_mat is the cost and alignment sequence between lcb blocks.
* NOTE: edit_cost is the total editing cost between lcb blocks (which is not included
* is cost_mat), full_cost_lstlst is a little bit redundant here, for it
* contains the two code array, but it also has the begin and end point of each
* block. *)
let get_matcharr_and_costmatrix seq1 seq2 min_lcb_ratio min_cover_ratio min_lcb_len max_lcb_len 
locus_indel_cost cost_mat use_ukk =
    let debug2 = false in
    (*transform poy costmatrix into hodx_matrix in mauve*)
    if debug_major_slice then checkpoint_major_slice "begin of get_matcharr_and_costmatrix";
    fill_in_hmatrix cost_mat;
    let len1 = Sequence.length seq1 and len2 = Sequence.length seq2 in
    let minlen = if (len1>len2) then len2 else len1 in
    let min_lcb_len,max_lcb_len = 
        (float_of_int minlen) *. min_lcb_len,
        (float_of_int minlen)*. max_lcb_len in
    let min_lcb_len,max_lcb_len = 
        int_of_float min_lcb_len, int_of_float max_lcb_len in
    let seq1arr = Sequence.to_array seq1 
    and seq2arr = Sequence.to_array seq2 in
    let in_seqarr = [|seq1arr;seq2arr|] in
    if debug_main then
        Printf.printf "min_lcb_len,max_lcb_len=%d,%d,call create_lcb_tbl\n%!" min_lcb_len max_lcb_len;
    let lcb_tbl,lcbs,code_list,full_range_lstlst, mum_tbl, seed2pos_tbl = 
        create_lcb_tbl in_seqarr min_lcb_ratio min_lcb_len min_cover_ratio
        max_lcb_len cost_mat use_ukk in
    if debug_main then begin
        if debug2 then 
            Hashtbl.iter (fun key record ->
            print_lcb record ) lcb_tbl;
        (* Printf.printf "lcbs is \n%!";  Block_mauve.print_lcblst lcbs;*)
        Printf.printf "original code list is \n%!";
        printIntListList code_list;
    end;
    if debug_major_slice then checkpoint_major_slice "got lcbs, start to fill in cost matrix";
    let base = List.length (List.hd code_list) in (*start number of non-lcb block*)
    let len_lst1 = List.length (List.hd full_range_lstlst) in
    let len_lst2 = List.length (List.nth full_range_lstlst 1) in
    let gen_gap_code = (len_lst1 + len_lst2) * 2 + 1 in
    (*we use ali_mat start from  (1).(1)*)
    let matlen = gen_gap_code + 1 in
    let gen_cost_mat = Array.make_matrix matlen matlen Utl.large_int in
    let empty_seq = Sequence.get_empty_seq () in
    let ali_mat = Array.make_matrix matlen matlen (0,empty_seq, empty_seq) in
    let seqNO = ref (-1) in
    let full_code_lstlst =
        List.map (fun full_range_lst ->
            seqNO := !seqNO +1 ;
            let start_num = !seqNO*len_lst1*2 in
            let coderef = ref base in
            List.map(fun (left,right,lcbkey,lcb_refcode) ->
                if lcb_refcode<>0 then (*a lcb block, use lcb refcode*)
                    ((from_ori_code lcb_refcode)+start_num,(left,right),lcbkey)
                else begin (*not lcb block, give it gapcode,lcbkey=[] anyway*)
                    coderef := !coderef +1;    
                    (from_ori_code !coderef)+start_num,(left,right),[]
                end
            )full_range_lst
        )full_range_lstlst 
    in
    if debug_main then begin
    Printf.printf "full (code,range) list is :\n%!";
    List.iter(fun full_code_lst ->
        List.iter (fun (code,(l,r),lcbkey) ->
            Printf.printf "[%d,(%d,%d),%!" code l r;
            printIntList2 lcbkey;
            Printf.printf "]; %!";
        ) full_code_lst;
        print_newline();
    )full_code_lstlst;
    end;
    let get_code_arr_from_fullcode_lst fullcode_lst = 
        Array.of_list ( (List.map (fun (code,(_,_),_) -> code) fullcode_lst) )
    in
    let code1_arr = 
        get_code_arr_from_fullcode_lst (List.hd full_code_lstlst)
    and code2_arr =
        get_code_arr_from_fullcode_lst (List.nth full_code_lstlst 1)
    in
    let edit_cost,indel_cost = fill_in_cost_ali_mat max_lcb_len cost_mat seq1 seq2 in_seqarr 
    (List.hd full_code_lstlst) (List.nth full_code_lstlst 1) 
    gen_gap_code locus_indel_cost ali_mat
    gen_cost_mat base len_lst1 use_ukk mum_tbl seed2pos_tbl lcb_tbl in
    (*let full_code_lstlst = List.map (fun full_code_lst ->
        List.map (fun (code,(l,r),_ ) -> code,(l,r)
        ) full_code_lst;
    ) full_code_lstlst
    in*)
    if debug_main then Printf.printf "end of main function in block_mauve. return with cost=(edit:%d,indel=%d)\n%!"
    edit_cost indel_cost;
    if debug_major_slice then
        checkpoint_major_slice "end of get_matcharr_and_costmatrix";
    code1_arr,code2_arr,gen_cost_mat,ali_mat,gen_gap_code,edit_cost,indel_cost,full_code_lstlst

(** [get_range_with_code] return the range of match block code1 and block code2. if
    * a block is all_gap_seq, return (totalsize,totalsize+len_of_block-1) so that
    * graphic_mauve can work on it. totalsize starts as the size of the sequence,
    * increase each time when meet an all_gap block. *)
let get_range_with_code code1 code2 full_code_lstlst gen_gap_code totalsize1 totalsize2 =
    let debug = false in
    if debug then Printf.printf "get range with code %d,%d,gapcode=%d,total size1,2=%d,%d%!" 
    code1 code2 gen_gap_code totalsize1 totalsize2;
    let left1,right1 = 
        List.fold_left (fun acc (codex,(left,right),_) ->
        if code1=codex then (left,right)
        else acc ) (-1,-1) (List.hd full_code_lstlst)
    in
    let left2,right2 = 
        List.fold_left (fun acc (codey,(left,right),_) ->
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
ali_mat gen_gap_code seqsize1 seqsize2 = 
        let debug = false in
        let len_lst1 = Array.length alied_gen_seq1 in
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
                "output2mauve: code1=%d,code2=%d,left/right1=%d,%d, left/right2=%d,%d\n%!"
                alied_code1 alied_code2 left1 right1 left2 right2;
                (*let oc = open_out_gen [Open_creat(*;Open_append*)] 0o666 filename in*)
                let seqlst1 = get_seqlst_for_mauve alied_seq1 in
                let seqlst2 = get_seqlst_for_mauve alied_seq2 in
                fprintf oc "> 1:%d-%d %s c=%d\n" (left1+1) (right1+1) dir1 cost;
                if debug && (right1-left1+1)<500 then 
                info_user_message "this sequence alone might too short for mauve graphic output (length<500)";
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

