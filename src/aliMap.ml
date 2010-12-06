(* POY 4.0 Beta. A phylogenetic analysis program using Dynamic Homologies.    *)
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
(* USA                                                                        *)

let () = SadmanOutput.register "AliMap" "$Revision: 2845 $"


(** AliMap module implements methods to align two general 
* character sequence allowing rearrangements *)



type chromPairAliPam_t = ChromPam.chromPairAliPam_t
type block_pam_t = Block.blockPam_t
type seed_t = Seed.seed_t
type block_t = Block.block_t
type order_t = ChromPam.order_t
type subseq_t = Subseq.subseq_t

let fprintf = Printf.fprintf

let print_intarr inarr = 
    Array.iter (fun x -> Printf.printf "%d %!" x ) inarr; 
    print_newline()

(** [create_gen_cost_mat subseq1_ls subseq2_ls global_map gen_gap_code 
*        seq1 seq2 cost_mat ali_pam] creates a general cost matrix 
* between [subseq1_ls] and [subseq2_ls. Each subseq is considered
* as a character state *)
let create_gen_cost_mat subseq1_ls subseq2_ls global_map gen_gap_code 
        seq1 seq2 cost_mat ali_pam =
    let debug = false in
    let len1 = List.length subseq1_ls in 
    let len2 = List.length subseq2_ls in 
    let len = 2 * (len1 + len2 + 1) in 
    let gen_cost_mat = Array.make_matrix len len Utl.large_int in 
    let set_cost code1 code2 cost = 
        gen_cost_mat.(code1).(code2) <- cost in 
    set_cost gen_gap_code gen_gap_code 0;
    let num_block = List.length global_map in 
    let block_gap_cost = Utl.large_int / (num_block + 1) in 
    if debug then
        Printf.printf "create_gen_cost_mat,len1/len2=%d/%d, block gap cost = %d(%d)\n%!"
    len1 len2 block_gap_cost num_block;
    List.iter (fun b -> 
                   let code1_id = b.Block.subseq1_id in 
                   let code2_id = b.Block.subseq2_id in 
                   let cost = b.Block.cost in 
                   if debug then begin
                       Block.print b;
                       Printf.printf " set cost for code1=%d,code2=%d,gap=%d,cost=%d\n%!" 
                       code1_id code2_id gen_gap_code cost;
                   end;
                   set_cost code1_id code2_id cost;
                   set_cost code2_id code1_id cost;
                   set_cost code1_id gen_gap_code block_gap_cost; 
                   set_cost gen_gap_code code1_id block_gap_cost;
                   set_cost code2_id gen_gap_code block_gap_cost; 
                   set_cost gen_gap_code code2_id block_gap_cost;
    ) global_map;
    let del_subseq1_ls = List.filter 
        (fun subseq -> List.length subseq.Subseq.block_id_ls = 0) subseq1_ls in 
    let del_subseq2_ls = List.filter 
        (fun subseq -> List.length subseq.Subseq.block_id_ls = 0) subseq2_ls in 
    let empty_seq = Sequence.get_empty_seq () in 
    let ali_mat = Array.make_matrix len len (empty_seq, empty_seq) in 
    let pair_gap subseq seq = 
       let id = subseq.Subseq.id in 
       (*del_cost is the cost for replacing un-recognized block with all_gaps*)
       let del_cost = Sequence.cmp_gap_cost ali_pam.ChromPam.locus_indel_cost
           (Subseq.get_subseq seq subseq)
       in
       if debug then
           Printf.printf "pair gap for %d,%d with %d\n%!" id gen_gap_code del_cost;
       set_cost id gen_gap_code  del_cost;
       set_cost gen_gap_code id  del_cost;
    in
    if debug then 
        Printf.printf "pair gap begin... \n%!";
    List.iter (fun subseq -> pair_gap subseq seq1) del_subseq1_ls;
    List.iter (fun subseq -> pair_gap subseq seq2) del_subseq2_ls;
    List.iter 
        (fun subseq1 ->
             let id1 = subseq1.Subseq.id in 
             let s1 = subseq1.Subseq.sta in 
             let e1 = subseq1.Subseq.en in 
	         let subseq1 = Sequence.sub seq1 s1 (e1 - s1 + 1) in  
             List.iter 
                 (fun subseq2 ->
                      let id2 = subseq2.Subseq.id in 
                      let s2 = subseq2.Subseq.sta in 
                      let e2 = subseq2.Subseq.en in 	              
	                  let subseq2 = Sequence.sub seq2 s2 (e2 - s2 + 1) in 
                      let alied_seq1, alied_seq2, cost, _ =  
                          Sequence.align2 subseq1 subseq2 cost_mat 
                      in
                      if debug then
                          Printf.printf "ali %d and %d, cost =%d\n%!" id1 id2 cost;
                      set_cost id1 id2  cost;          
                      set_cost id2 id1  cost;    
                      ali_mat.(id1).(id2) <- (alied_seq1, alied_seq2);
                      ali_mat.(id2).(id1) <- (alied_seq2, alied_seq1)
                 ) del_subseq2_ls
        ) del_subseq1_ls;
    gen_cost_mat, ali_mat

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

let get_range_with_code code1 code2 full_code_lstlst gen_gap_code =
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
    if (code1=gen_gap_code)&&(code2<>gen_gap_code) then
        left2,right2,left2,right2
    else if (code2=gen_gap_code)&&(code1<>gen_gap_code) then
        left1,right1,left1,right1
    else if (code1<>gen_gap_code)&&(code2<>gen_gap_code) then
        left1,right1,left2,right2
    else
        failwith "aliMap.get_range_with_code, \
        we have indel match up with indel, something is wrong"

let create_general_ali_mauve seq1 seq2 cost_mat ali_pam outputtofile old_cost =
    let debug = true and debug2 = false in
    let min_seed_num = ChromPam.get_min_seed_num ali_pam
    and min_lcb_ratio = ChromPam.get_min_lcb_ratio ali_pam
    and min_bk_penalty = ChromPam.get_min_bk_penalty ali_pam
    in
    if debug then 
        Printf.printf "====  create general ali with mauve, len1=%d,len2=%d\
         min lcb ratio = %f, min seed num = %d, min bk penalty = %d\n%!"
    (Sequence.length seq1) (Sequence.length seq2) min_lcb_ratio min_seed_num min_bk_penalty;
    let seq1arr = Sequence.to_array seq1 
    and seq2arr = Sequence.to_array seq2 in
    let in_seqarr = [|seq1arr;seq2arr|] in
    let lcb_tbl,lcbs,code_list,full_range_lstlst = 
    Block_mauve.create_lcb_tbl in_seqarr min_lcb_ratio min_seed_num min_bk_penalty in
    if debug then begin
        Hashtbl.iter (fun key record ->
        Block_mauve.print_lcb record
        ) lcb_tbl;
        (* Printf.printf "lcbs is \n%!";  Block_mauve.print_lcblst lcbs;*)
        Printf.printf "original code list is \n%!";
        Block_mauve.print_int_lstlst code_list;
    end;
    assert((List.length code_list)=2);
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
                    Block_mauve.get_lcb_key_by_range !seqNO (left,right) lcb_tbl
                in
                if lcb_refcode<>0 then (*a lcb block, use lcb refcode*)
                    ((GenAli.from_ori_code lcb_refcode)+start_num,(left,right))
                else begin (*not lcb block, give it a new code*)
                    coderef := !coderef +1;
                    ((GenAli.from_ori_code !coderef)+start_num,(left,right))
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
            let ori_code1 = GenAli.to_ori_code code1 
            and ori_code2 = GenAli.to_ori_code (code2 - len_lst1*2) in
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
                    Sequence.cmp_gap_cost ali_pam.ChromPam.locus_indel_cost subseq1 
                and del_cost2 = 
                    Sequence.cmp_gap_cost ali_pam.ChromPam.locus_indel_cost subseq2
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
    let cost, rc, alied_gen_seq1, alied_gen_seq2 = 
        GenAli.create_gen_ali_new code1_arr code2_arr gen_cost_mat gen_gap_code 
        ali_pam.ChromPam.re_meth ali_pam.ChromPam.circular false
    in
    if debug then begin
        Printf.printf "alied code1/code2 arr =\n%!";
        Utl.printIntArr alied_gen_seq1; 
        Utl.printIntArr alied_gen_seq2; 
    end;
    let _ = match outputtofile with 
    | Some filename ->
        (*let oc = open_out_gen [Open_creat(*;Open_append*)] 0o666 filename in
        let oc = open_out filename in*)
        let rewrite = match old_cost with
        | Some oldcost -> 
                if cost<oldcost then true else false
        | None -> true 
        in
        if rewrite then begin
            let oc = open_out filename in
            fprintf oc "#FormatVersion Mauve1\n";
            fprintf oc "#Sequence1File	blabla.in\n";
            fprintf oc "#Sequence1Entry	1\n";
            fprintf oc "#Sequence1Format	FastA\n";
            fprintf oc "#Sequence2File	blabla.in\n";
            fprintf oc "#Sequence2Entry	2\n";
            fprintf oc "#Sequence2Format	FastA\n";
            List.iter2 ( fun alied_code1 alied_code2 ->
            let ori_code1 = GenAli.to_ori_code alied_code1 
            and ori_code2 = GenAli.to_ori_code (alied_code2 - len_lst1*2) in
            let cost,alied_seq1,alied_seq2 = 
                ali_mat.(alied_code1).(alied_code2)
            in
            let dir1 = if ori_code1>0 then "+" else "-"
            and dir2 = if ori_code2>0 then "+" else "-"
            in
            let left1,right1,left2,right2 = 
                get_range_with_code alied_code1 alied_code2 full_code_lstlst gen_gap_code
            in
            assert(left1>=0); assert(left2>=0);
            (*let oc = open_out_gen [Open_creat(*;Open_append*)] 0o666 filename in*)
            fprintf oc "> 1:%d-%d %s bla.in,c=%d\n" (left1+1) (right1+1) dir1 cost;
            let seqlst1 = get_seqlst_for_mauve alied_seq1 in
            let seqlst2 = get_seqlst_for_mauve alied_seq2 in
            List.iter (fun seq1 -> 
                Sequence.print oc seq1 Alphabet.nucleotides;
                fprintf oc "\n";
            ) seqlst1;
            fprintf oc "> 2:%d-%d %s bla.in,c=%d\n" (left2+1) (right2+1) dir2 cost;
            List.iter (fun seq2 ->
                Sequence.print oc seq2 Alphabet.nucleotides;
                fprintf oc "\n";
            ) seqlst2;
            fprintf oc "=\n";
            ) (Array.to_list alied_gen_seq1) (Array.to_list alied_gen_seq2);
            close_out oc;
        end;
    | None -> ()
    in
    if debug then begin
        Printf.printf "mauve ali res : cost=%d(+%d),rc=%d,alied_gen_seq1/2 = %!" 
        cost !edit_cost rc;
        Block_mauve.print_int_list (Array.to_list alied_gen_seq1);
        Block_mauve.print_int_list (Array.to_list alied_gen_seq2);
        print_newline();
    end;
    let cost = cost + !edit_cost in
    full_code_lstlst, gen_gap_code, ali_mat,
    (Array.to_list alied_gen_seq1),(Array.to_list alied_gen_seq2), cost,(0,rc)



(** [create_general_ali state global_map seq1 seq2 cost_mat ali_pam] 
* returns a general alignement between [seq1] and [seq2] allowing rearrangements *)    
let create_general_ali state global_map seq1 seq2 cost_mat ali_pam =    
    let global_map, subseq1_ls, subseq2_ls = 
        Block.create_subseq_id `Both global_map ali_pam in
    let len1 = List.length subseq1_ls in 
    let subseq1_ls = List.map (fun sub -> 
                                    {sub with Subseq.id = sub.Subseq.id * 2 - 1}
                              ) subseq1_ls   
    in 
    let subseq2_ls = List.map (fun sub ->                                    
                                   {sub with Subseq.id = len1 * 2 + sub.Subseq.id * 2 - 1}
                              ) subseq2_ls   
    in 
    List.iter (fun b -> 
                   b.Block.subseq1_id <- b.Block.subseq1_id * 2 - 1;
                   b.Block.subseq2_id <- 
                       match b.Block.direction with 
                       | `Positive -> len1 * 2 + b.Block.subseq2_id * 2 - 1  
                       | _ -> len1 * 2 + b.Block.subseq2_id * 2

              ) global_map;
    let gen_gap_code = ((List.length subseq1_ls) + (List.length subseq2_ls)) * 2 + 1 in  
    let gen_cost_mat, ali_mat = create_gen_cost_mat subseq1_ls subseq2_ls 
        global_map gen_gap_code seq1 seq2 cost_mat ali_pam 
    in 
    let gen_seq1 = Array.map (fun sub -> sub.Subseq.id) (Array.of_list subseq1_ls) in 
    let gen_seq2 = Array.map 
        (fun sub ->                                   
             let id = 
                 try 
                     let b = List.find (fun b ->
                                            (b.Block.subseq2_id + 1)/2 = (sub.Subseq.id + 1)/2
                                       ) global_map
                     in 
                     b.Block.subseq2_id
                 with not_found -> sub.Subseq.id
             in 
             id 
        ) (Array.of_list subseq2_ls) 
    in 
    let cost, rc, alied_gen_seq1, alied_gen_seq2 = 
        GenAli.create_gen_ali_new gen_seq1 gen_seq2 gen_cost_mat gen_gap_code 
        ali_pam.ChromPam.re_meth ali_pam.ChromPam.circular false
    in   
    let recost = (0,rc) in
    subseq1_ls, subseq2_ls, gen_gap_code, global_map, ali_mat, 
    alied_gen_seq1, alied_gen_seq2, cost, recost
