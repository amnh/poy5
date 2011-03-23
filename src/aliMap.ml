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

let create_general_ali_mauve seq1 seq2 cost_mat ali_pam outputtofile old_cost =
    let debug = true in
    let min_cover_ratio:float = ChromPam.get_min_cover_ratio ali_pam
    and min_lcb_ratio = ChromPam.get_min_lcb_ratio ali_pam
    and min_bk_penalty = ChromPam.get_min_bk_penalty ali_pam
    in
    if debug then 
        Printf.printf "====  create general ali with mauve, len1=%d,len2=%d\
         min lcb ratio = %f, min cover R = %f, min bk penalty = %d\n%!"
    (Sequence.length seq1) (Sequence.length seq2) min_lcb_ratio min_cover_ratio min_bk_penalty;
    let code1_arr,code2_arr,gen_cost_mat,ali_mat,gen_gap_code,edit_cost,
    full_code_lstlst,len_lst1 =
        Block_mauve.get_matcharr_and_costmatrix seq1 seq2 min_lcb_ratio min_cover_ratio
        min_bk_penalty ali_pam.ChromPam.locus_indel_cost cost_mat  in
    let cost, rc, alied_gen_seq1, alied_gen_seq2 = 
        GenAli.create_gen_ali_new code1_arr code2_arr gen_cost_mat gen_gap_code 
        ali_pam.ChromPam.re_meth ali_pam.ChromPam.circular false in
    if debug then begin
        Printf.printf "alied code1/code2 arr =\n%!";
        Utl.printIntArr alied_gen_seq1; 
        Utl.printIntArr alied_gen_seq2; 
    end;
    let _ = match outputtofile with 
    | Some filename ->
        Block_mauve.output2mauvefile filename cost old_cost alied_gen_seq1 alied_gen_seq2 
        full_code_lstlst ali_mat gen_gap_code len_lst1 (Sequence.length seq1)
        (Sequence.length seq2)
    | None -> ()
    in
    if debug then begin
        Printf.printf "mauve ali res : cost=%d(+edit cost=%d),rc=%d,alied_gen_seq1/2 = %!" 
        cost edit_cost rc;
        Block_mauve.print_int_list (Array.to_list alied_gen_seq1);
        Block_mauve.print_int_list (Array.to_list alied_gen_seq2);
        print_newline();
    end;
    let cost = cost + edit_cost in
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
