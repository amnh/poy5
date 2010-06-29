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

let () = SadmanOutput.register "GenAli" "$Revision: 1616 $"

(** This module implements methods to align two general
* characters allowing rearrangements *)

let fprintf = Printf.fprintf;;
type dyna_state_t = Data.dyna_state_t


let sanity_check arr size =
    let mark_arr = Array.make size 0 in
    Array.iter (fun x ->
        if (x<size) then
        mark_arr.(x) <- mark_arr.(x)+1
    ) arr;
    Array.fold_left (fun sign x -> if (x>1) then (sign+1) else sign ) 0 mark_arr

(** [oritentation code] returns the direction of the 
* character [code]. Characters whose codes are even numbers are
* considered as negative characters, otherwise postive *)
let get_orientated_code code = 
    if code mod 2 = 0 then -(code / 2)
    else (code + 1) /2

let to_ori_arr arr =
    Array.init (Array.length arr) (fun index->
        assert((arr.(index)<>0));
        if( arr.(index) mod 2 == 0) then -(arr.(index)/2)
        else (arr.(index)+1)/2 )

let from_ori_arr arr = 
    Array.init (Array.length arr) (fun index->
        assert((arr.(index)<>0));
        if (arr.(index)<0) then (abs arr.(index))*2
        else (arr.(index))*2-1
        )

let equal_orientation code1 code2 = compare (abs code1) (abs code2) 

let equal_content arr1 arr2 =
    if (Array.length arr1)=(Array.length arr2) then begin
        let oriarr1 = to_ori_arr (arr1)
        and oriarr2 = to_ori_arr (arr2) in
        let common_ori_arr1,common_ori_arr2 = Utl.get_common oriarr1 oriarr2
        equal_orientation in
        if (Array.length common_ori_arr1)=(Array.length arr1) then true
        else false
    end
    else false

(* [cmp_recost_simple seq1 seq2] is similar to [cmp_recost], except it only take
* two sequence seq1 and seq2 as input, no 'reseq2' *)
let cmp_recost_simple seq1 seq2 re_meth circular orientation =
    let arr1, arr2 = Utl.get_common seq1 seq2 equal_orientation in 
    let res = 
        match re_meth with 
        | `Locus_Inversion cost -> 
            (UtlGrappa.cmp_inversion_dis arr1 arr2 circular) * cost  
        | `Locus_Breakpoint cost ->       
            (UtlGrappa.cmp_oriented_breakpoint_dis arr1 arr2 circular) * cost
    in
    (*debug msg 
    Printf.printf "cmp_recost_simple on seq/reseq =\n%!";
    Utl.printIntArr seq1;
    Utl.printIntArr seq2;
     Printf.printf "res = %d\n%!" res;
    debug msg*)
    res
 


(** [cmp_recost state seq1 seq2 reseq2 re_meth circular] returns
* the rearrangement distance between two sequence [seq1] and [seq2] *)
let cmp_recost state seq1 seq2 reseq2 re_meth circular orientation = 
    let seq1, seq2, reseq2 = match orientation with
    | true -> to_ori_arr seq1, to_ori_arr seq2, to_ori_arr reseq2 
    | false -> seq1, seq2, reseq2
    in
    if Array.length seq2 = 0 then 0, 0
    else begin
        let recost1 = 
            match state with
            | `Breakinv ->                                              
                  let com_seq1_arr, com_reseq2_arr = Utl.get_common seq1 reseq2 equal_orientation in 
                  (match re_meth with 
                  | `Locus_Inversion cost -> 
                        (UtlGrappa.cmp_inversion_dis com_seq1_arr com_reseq2_arr circular) * cost  
                  | `Locus_Breakpoint cost ->       
                        (UtlGrappa.cmp_oriented_breakpoint_dis com_seq1_arr com_reseq2_arr circular) * cost)                     
            | _ -> 0
        in 
        let recost2  = 
            let seq2, reseq2 =  Utl.get_common seq2 reseq2 equal_orientation in
            (match re_meth with 
            | `Locus_Inversion cost -> 
                  (UtlGrappa.cmp_inversion_dis seq2 reseq2 circular) * cost  
            | `Locus_Breakpoint cost -> 
                  (UtlGrappa.cmp_oriented_breakpoint_dis seq2 reseq2 circular) * cost   
            )
        in  
        (* debug msg 
            Printf.printf "cmp_recost: [%!";
            Array.iter (Printf.printf "%d,") seq1; Printf.printf "];[%!";
            Array.iter (Printf.printf "%d,") seq2; Printf.printf "];[%!";
            Array.iter (Printf.printf "%d,") reseq2; Printf.printf "]: %!";
            Printf.printf "recost1=%d,recost2=%d \n%!" recost1 recost2;
         debug msg *)
        recost1, recost2
    end 


(** [cmp_cost state code1_arr code2_arr recode2_arr 
*              cost_mat gap re_meth circular] returns
* the total cost between [seq1] and [reseq2]. Precisely,
* total cost = editing cost ([seq1], [reseq2]) + rearrangement cost ([seq2], [reseq2]) *)
let cmp_cost state code1_arr code2_arr recode2_arr 
        (cost_mat : Cost_matrix.Two_D.m) gap re_meth circular orientation = 
    let seq1 = Sequence.of_array code1_arr
    and reseq2 = Sequence.of_array recode2_arr in
    let alied_seq1, alied_reseq2, editing_cost =  
        Sequence.Align.align_2 ~first_gap:false seq1 reseq2 cost_mat
            Matrix.default  
    in  
    let alied_code1_arr = Sequence.to_array alied_seq1 in 
    let alied_recode2_arr = Sequence.to_array alied_reseq2 in 
    let recost1, recost2 = 
        cmp_recost state code1_arr code2_arr recode2_arr re_meth circular 
        orientation 
    in 
   (editing_cost + recost1 + recost2), (recost1, recost2), alied_code1_arr, alied_recode2_arr

(** [find_wagner_ali state seq1 seq2 gen_cost_mat gap re_meth circular]
 * returns rearranged sequence [reseq2] of sequence [seq2] using stepwise addition method 
 * such that the total cost is minimum where 
 * total cost = editing cost ([seq1], [reseq2]) + rearrangement cost ([seq2], [reseq2]) *)
let find_wagner_ali (kept_wag : int) state seq1 seq2 gen_cost_mat gap re_meth circular orientation =
    (*debug msg
    Printf.printf "find_wagner_ali seq1/seq2=[%!"; 
    Array.iter (Printf.printf "%d,") seq1;
    Printf.printf "][%!";
    Array.iter (Printf.printf "%d,") seq2;
    Printf.printf "]\n%!";
    debug msg*)
    let rec add (best_wagner_seq2_arr : int array array) added_seq2_ls rem_seq2_ls =
        match rem_seq2_ls with
        | [] -> best_wagner_seq2_arr
        | code2 :: tl ->
              let added_seq2_ls = added_seq2_ls @ [code2] in 
              let wagner_seq2_ls = ref [] in 
              let update partial_seq2 =
                  let cost, (_, _), _, _  = 
                      cmp_cost state seq1 (Array.of_list added_seq2_ls) partial_seq2 
                          gen_cost_mat gap re_meth circular orientation
                  in
                  wagner_seq2_ls := (partial_seq2, cost)::!wagner_seq2_ls;
              in 
              
              let num_w = Array.length best_wagner_seq2_arr in  
              let len2 = Array.length best_wagner_seq2_arr.(0) in
              for w = 0 to num_w - 1 do                     
                  for pos = 0 to len2 do
                      let partial_seq2 = 
                          Utl.insert best_wagner_seq2_arr.(w) pos code2 in 
                      update partial_seq2;
                       if (state = `Annotated) && (code2 mod 2 = 1) then begin
                            let partial_seq2 = 
                                Utl.insert best_wagner_seq2_arr.(w) pos (code2 + 1) 
                            in 
                            update partial_seq2;
                      end 
                  done;     
              done; 
              let subseq2 = Array.sub seq2 0 (len2 + 1) in
              update subseq2; 
              let wagner_seq2_ls = Sort.list (fun (_, c1) (_, c2) -> c1 < c2) !wagner_seq2_ls in
              let best_w_arr : int array array = Array.init (min kept_wag (List.length wagner_seq2_ls)) 
                        (fun p -> 
                            let w2, c = List.nth wagner_seq2_ls p in 
                            w2
                        )  
              in
    (*debug msg
              Printf.printf "best_w_arr = \n%!";
              Array.iter (fun x -> 
                  Printf.printf "[%!"; 
                  Array.iter (Printf.printf "%d,%!") x; 
              Printf.printf "] %!")best_w_arr;
    debug msg*)
              add best_w_arr added_seq2_ls tl
    in
    let wagner_seq2_arr = add [|[||]|] [] (Array.to_list seq2) in 
    (*debug msg
    Printf.printf "result is [:%!";
    Array.iter (Printf.printf "%d ") wagner_seq2_arr.(0);
    Printf.printf "]\n%!";
    debug msg*)
    wagner_seq2_arr.(0)




(** [multi_swap_locus state seq1 seq2 best_seq2 best_cost 
*                     gen_cost_mat gap re_meth max_swap_med circular num_done_swap] 
* swaps [reseq2] in order to minimize the total cost between [seq1] and [reseq2] where 
* total cost = editing cost ([seq1], [reseq2]) + rearrangement cost ([seq2], [reseq2]) *)
let rec multi_swap_locus state seq1 seq2 best_seq2 best_cost 
        gen_cost_mat gap re_meth max_swap_med
        circular orientation num_done_swap =             
    let len2 = Array.length best_seq2 in  
    let swap_ls = ref [] in 
    for donor_pos = 0 to len2 - 2 do
        for rev_pos = donor_pos + 1 to len2 - 1 do 
            let new_seq2 = Utl.swap_item donor_pos rev_pos best_seq2 in 

            let new_cost, (_, _), _, _ = cmp_cost state seq1 seq2 new_seq2 
                gen_cost_mat gap re_meth circular orientation
            in                 
            if new_cost < best_cost then 
                swap_ls := (donor_pos, rev_pos, new_cost)::!swap_ls
       done 
    done;
    let swap_ls = 
        List.sort (fun (_, _, s1) (_, _, s2) -> compare s1 s2) !swap_ls
    in 
    let is_depend (s1, e1, _) (s2, e2, _) = 
        if (s1 > e2) || (e1 < s2) 
            || ( (s1 > s2) && (e1 < e2) ) 
              ||( (s2 > s1) && (e2 < s1) )  then false
        else true
    in 
    let ind_swap_ls = List.fold_left 
        (fun ind_swap_ls swap ->
             if List.exists (fun swap2 -> is_depend swap swap2) ind_swap_ls 
             then ind_swap_ls 
             else swap::ind_swap_ls) [] swap_ls 
    in 
    let ind_swap_arr = Array.of_list (List. rev ind_swap_ls) in 
    let rec swap num_swap = 
        let new_seq2 = Array.copy best_seq2 in
        for pos = 0 to num_swap - 1 do 
            let sta, en, _ = ind_swap_arr.(pos) in 
            let tmp = new_seq2.(sta) in
            new_seq2.(sta) <- new_seq2.(en);
            new_seq2.(en) <- tmp;
        done;
            
        let new_cost, (_, _), _, _ = cmp_cost state seq1 seq2 new_seq2 
            gen_cost_mat gap re_meth circular orientation
        in                  
        if new_cost < best_cost then begin
            new_cost, new_seq2
        end 
        else swap (num_swap  / 2)
    in 
    let num_swap = Array.length ind_swap_arr in 
    if (num_swap = 0) then  best_cost, best_seq2
    else begin
        let new_cost, new_seq2 = swap num_swap in 
        if num_done_swap + 1 >=  max_swap_med then new_cost, new_seq2 
        else 
            multi_swap_locus state seq1 seq2 new_seq2 new_cost gen_cost_mat gap 
                re_meth max_swap_med circular orientation (num_done_swap + 1)
    end 

(* I add these function and structure "matched_loci_array" for annotated chromosome, 
* They should work for other data type too, 
* but we need to verify that -- verifying....maybe we
* should change the name "loci" to something else....*)
type matched_loci_array = {
        cost : int;
        size : int; (* size of matched array *)
        (* since we are going to take gap/indel into match, size above won't be
        * able to tell us whether the match is done, we need another two number
        * here, arr1_left/arr2_left is the number of un-matched item in array1/array2
        * *)
        arr1_left : int; 
        arr2_left : int;
        marked_matrix : int array array  ; (* 1 means not-available, 0 is ok *)
        cost_array : (int * int * int) array ; (* ( index i, index j, cost) *)
        matched_loci_list : (int * int * int) list ;(* (index i, index j, cost)*)
}

let mark_matrix (in_matrix:int array array) indexi indexj mark_i_or_j_only =
    if (mark_i_or_j_only=0) then
        Array.mapi (fun idi arr ->
            Array.mapi ( fun idj item ->
                if (idi = indexi) then 1
                else if (idj = indexj) then 1
                else item
            )arr
        )in_matrix
    else if (mark_i_or_j_only=1) then
        Array.mapi (fun idi arr ->
            Array.mapi ( fun idj item ->
                if (idi = indexi) then 1
                else item
            ) arr
        )in_matrix
    else (* if (mark_i_or_j_only=2) then *)
        Array.mapi (fun idi arr ->
            Array.mapi ( fun idj item ->
                if (idj = indexj) then 1
                else item
            ) arr
        )in_matrix



let make_new_mark_matrix in_matrix indexi indexj gapcode =
    if (indexi=gapcode)&&(indexj<>gapcode) then
        let leftidxj = 
            if (indexj mod 2)=0 then indexj-1 else indexj
        in
        let mm = mark_matrix in_matrix indexi leftidxj 2 in
        let mm = mark_matrix mm indexi (leftidxj+1) 2 in
        mm
    else if (indexj=gapcode)&&(indexi<>gapcode) then
        let upidxi =
            if (indexi mod 2)=0 then indexi-1 else indexi 
        in
        let mm = mark_matrix in_matrix upidxi indexj 1 in
        let mm = mark_matrix mm (upidxi+1) indexj 1 in
        mm
    else if (indexi<>gapcode)&&(indexj<>gapcode) then
        let lefttopi, lefttopj = 
            match (indexi mod 2),(indexj mod 2) with
            | 1,1 -> indexi, indexj
            | 1,0 -> indexi, indexj-1
            | 0,1 -> indexi-1, indexj
            | 0,0 -> indexi-1, indexj-1
            | _,_ -> failwith ("what? besides 0 and 1, you find different result\
            for (x mod 2)?\n")
        in
        let mm = mark_matrix in_matrix lefttopi lefttopj 0 in
        let mm = mark_matrix mm lefttopi (lefttopj+1) 0 in
        let mm = mark_matrix mm (lefttopi+1) lefttopj 0 in
        let mm = mark_matrix mm (lefttopi+1) (lefttopj+1) 0 in
        mm
    else
        failwith ("why we are matching gap with gap?")

let array_filter_nth in_array index =
    let len = Array.length in_array in
    assert(index>=0); assert(index<len);
    let arra = 
        if (index=0) then [||]
        else Array.sub in_array 0 index 
    in
    let arrb = 
    if (index+1 < len) then Array.sub in_array (index+1) (len-index-1)
    else [||] 
    in
    Array.append arra arrb 
    
let get_best_n (in_list:matched_loci_array list) size =
    assert(size <= List.length in_list);
    let (in_array:matched_loci_array array) = Array.of_list in_list in
    Array.sort (fun item1 item2 ->
        let cost1 = item1.cost and cost2 = item2.cost in
        compare cost1 cost2
    ) in_array ;
    let out_array = Array.sub in_array 0 size in
    Array.to_list out_array

(* this function takes the result cost matrix from 
*  function [create_pure_gen_cost_mat] as in_cost_mat. 
*  get the part of matrix we need as out_cost_mat, also get a sorted cost array out of it. *)
let make_cost_matrix_and_array in_cost_mat code1_arr code2_arr
gapcode re_meth
=    
    (* we will consider gap(indel) in the match *)
    let lst1 = (List.sort compare (Array.to_list code1_arr)) 
    and lst2 = (List.sort compare (Array.to_list code2_arr)) in
    let end1 = List.nth lst1 ((List.length lst1)-1)
    and end2 = List.nth lst2 ((List.length lst2)-1)  in
    let len1 = if (end1 mod 2)=0 then end1+2 else end1+3
    and len2 = if (end2 mod 2)=0 then end2+2 else end2+3 in
    (*debug msg
    Printf.printf "make cost matrix and array, code1_arr/code2_arr = \n %!";
    Utl.printIntArr code1_arr; Utl.printIntArr code2_arr; 
    debug msg*)
    let cost_list = ref [] in
    let out_cost_mat = 
        Array.make_matrix len1 len2 Utl.large_int in
    (* the first column and line of above matrix will not be used *)
    Array.iter(fun code1 ->
        Array.iter (fun code2 ->
            let i = if (code1 mod 2)<>0 then code1
            else (code1-1) in
            let j = if (code2 mod 2)<>0 then code2
            else (code2-1) in
            let cost = in_cost_mat.(i).(j) in
            out_cost_mat.(i).(j) <- cost;
            cost_list := (i,j,cost) :: !cost_list;
            let cost = in_cost_mat.(i+1).(j) in
            out_cost_mat.(i+1).(j) <- cost;
            cost_list := (i+1,j,cost) :: !cost_list;
            let cost = in_cost_mat.(i).(j+1) in
            out_cost_mat.(i).(j+1) <- cost;
            cost_list := (i,j+1,cost) :: !cost_list;
            let cost = in_cost_mat.(i+1).(j+1) in
            out_cost_mat.(i+1).(j+1) <- cost;
            cost_list := (i+1,j+1,cost) :: !cost_list;
        ) code2_arr;
    )code1_arr;
    (*match with gap*)
    Array.iter (fun code1 ->
            let i = if (code1 mod 2)<>0 then code1
            else (code1-1) in
            let cost = in_cost_mat.(i).(gapcode) in
            cost_list := (i,gapcode,cost) :: !cost_list;
            let cost = in_cost_mat.(i+1).(gapcode) in
            cost_list := (i+1,gapcode,cost) :: !cost_list;
    ) code1_arr;
    Array.iter (fun code1 ->
            let i = if (code1 mod 2)<>0 then code1
            else (code1-1) in
            let cost = in_cost_mat.(i).(gapcode) in
            cost_list := (gapcode,i,cost) :: !cost_list;
            let cost = in_cost_mat.(i+1).(gapcode) in
            cost_list := (gapcode,i+1,cost) :: !cost_list;
    ) code2_arr;
    let cost_array = Array.of_list ( 
        List.sort (fun (_,_,cost1) (_,_,cost2) -> compare cost1 cost2
        ) !cost_list )
    in
(*debug msg 
  (*  Printf.printf "check in cost matrix:\n%!"; Utl.printIntMat in_cost_mat;*)
    Printf.printf "check cost array:\n%!";
    Array.iter (fun (id1,id2,cost) ->  Printf.printf "(%d,%d,%d) " id1 id2 cost;
    ) cost_array;   Printf.printf "\n\n%!";
 debug msg*)
    out_cost_mat, cost_array

(* heuristic function to match pair of input loci arrays.
*  h_size is the heuristic size -- the number of best choice we keep in each
*  round. *)
let match_pair_heuristic full_cost_array marked_matrix size1 size2 h_size
gapcode =
    let candidate_list = ref [] in
    let cost_array = full_cost_array in
    for x = 0 to (h_size-1) do
        let (indexi, indexj, costij) = cost_array.(x) in
        let costarr_left = array_filter_nth cost_array x in
        let (init_match: matched_loci_array) = {
            cost = costij; 
            size = 1;
            arr1_left = if (indexi<gapcode) then size1-1 else size1 ;
            arr2_left = if (indexj<gapcode) then size2-1 else size2 ;
            matched_loci_list = [(indexi, indexj, costij)];
            marked_matrix = make_new_mark_matrix marked_matrix indexi indexj gapcode;
            cost_array = costarr_left;
        } in
        candidate_list := (!candidate_list) @ [init_match] ;
    done;
    let somethingleft = ref 1 in
    while ( !somethingleft = 1 ) do
        let current_list = ref [] in
        let (tmp_list:matched_loci_array list) = (!candidate_list) in
        List.iter (fun (item:matched_loci_array) ->
            let base_size = item.size and base_cost = item.cost in
            let base_left1 = item.arr1_left and base_left2 = item.arr2_left in
            let matched_list = item.matched_loci_list in
            let marked_matrix = item.marked_matrix in
            let cost_array = item.cost_array in
(*debug msg
            Printf.printf "left1=%d,left2=%d, gapcode = %d, check matched_list:\n%!"
            base_left1 base_left2 gapcode;
            List.iter (fun (id1,id2,cost) ->  Printf.printf "(%d,%d,%d) " id1 id2 cost;
            ) matched_list; Printf.printf "\n%!";
(*            Printf.printf "check cost array:\n%!";
            Array.iter (fun (id1,id2,cost) ->  
                Printf.printf "(%d,%d,%d) " id1 id2 cost;
            ) cost_array; Printf.printf "\n%!";
*)
 (*       Printf.printf "check marked matrix:\n%!";  Utl.printIntMat
 *       marked_matrix;*)
debug msg*)
            if( (base_left1=0)&&(base_left2=0) ) then begin
                current_list := !current_list @ [item];
                somethingleft := 0
            end 
            else begin
                let picked = ref 0 and left_len = ref (Array.length cost_array)  
                and z = ref 0 in
                while ( (!picked < h_size)&&(!z < !left_len)&&(!left_len >0) ) do
                    let (indexi, indexj, costij) =  cost_array.(!z) in
                    let costarr_left = array_filter_nth cost_array !z in
                    if (indexi >= Array.length(marked_matrix) ) then
                    Printf.printf "idi=%d > Array.sizei=%d\n%!" 
                    indexi (Array.length marked_matrix);
                    if (indexj >= Array.length(marked_matrix.(indexi))) then
                        Printf.printf "idj=%d > Array.sizej=%d\n%!" indexj 
                        ( Array.length marked_matrix.(indexi));
                    if ( marked_matrix.(indexi).(indexj) = 0 ) then begin
                        let new_matched_list = (indexi, indexj, costij) :: matched_list in
                        let new_marked_matrix = 
                            make_new_mark_matrix marked_matrix indexi indexj
                            gapcode in
                        let new_match_item = {
                            cost=costij+base_cost;
                            size = base_size +1;
                            arr1_left = 
                                if (indexi<gapcode) then base_left1-1
                                else base_left1;
                            arr2_left = 
                                if (indexj<gapcode) then base_left2-1
                                else base_left2;
                            matched_loci_list = new_matched_list;
                            marked_matrix = new_marked_matrix;
                            cost_array = costarr_left;
                        } in
                        current_list := (!current_list) @ [new_match_item];
                        picked := !picked + 1;
                    end else ();
                    left_len := Array.length costarr_left; 
                    z := !z +1;
                done;
            end 
            ;
        ) tmp_list ;
        candidate_list := 
            (*if (!somethingleft = 1) then *)
            if ((List.length !current_list)>0) then 
                get_best_n !current_list h_size 
            else
                tmp_list ;
    done;
    !candidate_list

let get_neg_code code = 
    if ( (code mod 2)=1 ) then code + 1 else code - 1

(* match arr1=[i1,i2,i3...] with arr2=[j1,j2,j3....], try to get a lower sum of 
*  all cost(i,j), i=i1,i2.., j=j1,j2... 
*  we don't consider adding indel(all gaps sequence) here, call algn function in
*  sequence.ml to take care of that.
*  *)
let match_pair arr1 arr2 cost_array sizex sizey gapcode =
    (* heuristic size: how many different match we keep at every step*)
    let heuristic_size = 2 in 
    let size1 = Array.length arr1 
    and size2 = Array.length arr2 in
    let marked_matrix = Array.make_matrix sizex sizey 0 in
    let (res_list:matched_loci_array list) = 
        match_pair_heuristic cost_array marked_matrix size1 size2 heuristic_size gapcode in
    let best = List.hd res_list in
    let matched_cost = best.cost in
    let matched_list = best.matched_loci_list in
    let matched_list = List.map (fun (code1,codem,_) ->  code1, codem
    ) matched_list in
    (*debug msg
    Printf.printf "match pair with arr1/arr2=\n%!";
    Utl.printIntArr arr1; Utl.printIntArr arr2;
    Printf.printf "editng cost =%d, matched list is [%!" matched_cost;
    List.iter (fun (idx1,idx2) -> 
        Printf.printf "(%d,%d) " idx1 idx2 
    ) matched_list;
    Printf.printf "]\n  %!"; 
    debug msg*)
    matched_list, 
    matched_cost


(* [create_gen_ali_new] take code1_arr codem_arr as input code array, c2 is the
* 2d matrix and cost_matrix is the "pure_cost_matrix" -- which is the "int array
* array" version of c2. It then calls "make_cost_matrix_and_array" to get the
* cost_array for the function "match_pair", which will give us a match, or you
* can call it alignment between codes in the two input array, also the editing
* cost. To get the rerrangement cost, we call "cmp_recost_simple" *)
let create_gen_ali_new state code1_arr codem_arr c2 (*we don't need c2 in this
function , remove it later*) 
cost_matrix gapcode re_meth circular orientation  = 
    let sizex = Array.length cost_matrix in
    let sizey = Array.length cost_matrix.(0) in
    let cost_mat, cost_array = 
        make_cost_matrix_and_array cost_matrix code1_arr codem_arr gapcode re_meth
    in
    let matched_list, matched_cost = 
        match_pair code1_arr codem_arr cost_array sizex sizey gapcode in
    let matched_array = Array.of_list matched_list in
    let nongap_matched_lst = ref [] in
    let alied_code1_lst = ref [] and alied_codem_lst = ref [] in
    Array.iter (fun item ->
        let trypos = 
                Utl.find_index matched_array item 
                (fun looking_item (indexi,_) -> compare looking_item indexi ) 
        in
        if (trypos<0) then begin 
                let neg_item = get_neg_code item in
                alied_code1_lst := !alied_code1_lst @ [item];
                let pos = 
                    Utl.find_index matched_array neg_item
                    (fun looking_item (indexi,_) -> compare looking_item indexi)
                in
                assert(pos>=0);
                let _,matched_item =  matched_array.(pos) in
                let neg_matched_item = 
                    if (matched_item = gapcode) then gapcode 
                    else  get_neg_code matched_item in
                alied_codem_lst := !alied_codem_lst @ [neg_matched_item];
                if (neg_matched_item<gapcode) then
                 nongap_matched_lst :=
                  (item,neg_matched_item)::(!nongap_matched_lst)
        end
        else begin
                alied_code1_lst := !alied_code1_lst @ [item];
                let _,matched_item =  matched_array.(trypos) in
                alied_codem_lst := !alied_codem_lst @ [matched_item];
                if (matched_item<gapcode) then
                 nongap_matched_lst :=
                  (item,matched_item)::(!nongap_matched_lst)
        end;
    ) code1_arr;
    Array.iteri (fun idx codem ->
        if (List.mem codem !alied_codem_lst)||
        (List.mem (get_neg_code codem) !alied_codem_lst)  then ()
        else begin 
            alied_code1_lst := !alied_code1_lst @ [gapcode];
            alied_codem_lst := !alied_codem_lst @ [codem];
        end
    )codem_arr;
    let alied_code1_lst = !alied_code1_lst 
    and alied_codem_lst = !alied_codem_lst in
    let nongap_matched_lst = List.rev !nongap_matched_lst in
    let codem_matched_with_nongap = Array.map ( fun (_ ,codem) -> codem ) 
    (Array.of_list nongap_matched_lst) in
    let arr1, arr2 = to_ori_arr codem_matched_with_nongap,
    to_ori_arr codem_arr in
    let recost =  
        cmp_recost_simple arr1 arr2 re_meth circular orientation in
    let editingcost = matched_cost in
    let alied_code1,alied_codem = 
        Array.of_list alied_code1_lst, Array.of_list alied_codem_lst in
    editingcost+recost, recost, alied_code1, alied_codem


  


let create_gen_ali kept_wag state (seq1 : Sequence.s) (seq2 : Sequence.s) 
        (gen_cost_mat : Cost_matrix.Two_D.m) pure_gen_cost_mat alpha re_meth 
        max_swap_med circular orientation =
    let gapcode = Alphabet.get_gap alpha in
  (* debug msg 
    Printf.printf "\n GenAli.ml create_gen_ali: gapcode=%d,seq1/seq2=\n%!"
    gapcode ;
    Sequence.printseqcode seq1; Sequence.printseqcode seq2;
  (*  Printf.printf "check cost matrix:\n%!";
    Array.iter(fun arr -> Utl.printIntArr arr) pure_gen_cost_mat; *)
   debug msg *)
    let arr1 = Sequence.to_array seq1 in 
    let arr2 = Sequence.to_array seq2 in
    let tc, rc, alied_code1, alied_code2 = 
        create_gen_ali_new state arr1 arr2 gen_cost_mat
        pure_gen_cost_mat gapcode re_meth circular orientation
    in
    (*debug msg
    Printf.printf "tc=%d,rc=%d, alied1/alied2 = \n%!" tc rc;
    Utl.printIntArr alied_code1;  Utl.printIntArr alied_code2; 
    debug msg*)
    let alied_seq1 = Sequence.of_array alied_code1 in
    let alied_seq2 = Sequence.of_array alied_code2 in
    tc, (0,rc), alied_seq1, alied_seq2


(* Note: seqXX here is result sequence of alignment before. 
*  when seq21 and seq22 are different, we need to rearrange seq11 and seq12, so
* that each item in seq11 match that in seq12. 
* For example, seq11 = [1,2,3,-,-], seq12 = [1,-,-,4,5], 
* seq21 = [1,-,6,7], seq22=[1,4,5,-] // I know this is a bad alignment, but still an
* alignment//. Now, we need to match [1,-,-,4,5] with [1,4,5,-], along with
* seq11 and seq21.
* *)
let re_align seq11 seq12 seq21 seq22 gapcode = 
    if (Sequence.compare seq12 seq22) = 0 then
        seq11,seq12,seq21
    else
        let len1 = Sequence.length seq12 and len2 = Sequence.length seq22 in
        let arr11,arr12,arr21,arr22 =
            if len1>=len2 then 
                Sequence.to_array seq11, Sequence.to_array seq12,
                Sequence.to_array seq21, Sequence.to_array seq22
            else Sequence.to_array seq21,Sequence.to_array seq22,
                 Sequence.to_array seq11,Sequence.to_array seq12
        in
        let finditem_in_array arr item num =
            let arrlen = Array.length arr in
            let res = ref arrlen and count = ref num in
            for i=0 to (arrlen-1) do
                if (arr.(i) = item) then 
                    if (!count)=1 then res := i
                    else count := (!count)-1
            done;
            (!res)
        in
        let arrlen = Array.length arr12 and arrlen_short = Array.length arr22 in
        let new_arr = Array.init arrlen 
        (fun x ->
            let gapcount = ref 0 in
            let item1 = Array.get arr12 x in
            let index2 = 
                if item1=gapcode then begin 
                    gapcount := (!gapcount)+1;
                    finditem_in_array arr22 item1 (!gapcount);
                end
                else finditem_in_array arr22 item1 1;
            in
            if (index2>=0)&&(index2<arrlen_short) then
                arr21.(index2)
            else
                gapcode
                ;
        ) in
        let new_seq = Sequence.of_array new_arr in
        if len1>=len2 then seq11,new_seq,seq12
        else new_seq,seq21,seq22

let create_gen_ali3_by_medsov medsov kept_wag (seq1 : Sequence.s) (seq2 :
    Sequence.s) (seq3 : Sequence.s) (delimiter_lstlst : int list list) (gen_cost_mat: Cost_matrix.Two_D.m) alpha re_meth  max_swap_med circular orientation sym =
    (* debug msg 
    Printf.printf "create_gen_ali3_by_medsov, seq1,seq2,seq3=\n%!";
    Sequence.printseqcode seq1; Sequence.printseqcode seq2;
    Sequence.printseqcode seq3;
    Printf.printf "delimiters lst =\n%!";
    List.iter (fun lst -> 
        Printf.printf "[%!"; 
        List.iter (Printf.printf "%d,") lst; 
        Printf.printf "];%!") delimiter_lstlst;
    Printf.printf "\n%!";
     debug msg *)
    let gapcode = Alphabet.get_gap alpha in
    let arr1 = Sequence.to_array seq1
    and arr2 = Sequence.to_array seq2 
    and arr3 = Sequence.to_array seq3 in
    let oriarr1 = to_ori_arr (arr1)
    and oriarr2 = to_ori_arr (arr2) 
    and oriarr3 = to_ori_arr (arr3) in
    let comoriarr1,comoriarr2,comoriarr3 = Utl.get_common3 oriarr1 oriarr2 oriarr3 equal_orientation in
    let ori_arr_med3,delimiter_arr =
        match medsov with
        |`Vinh ->
                failwith "Vinh median solver is not in grappa"
        |`MGR
        |`SimpleLK
        |`ChainedLK 
        |`COALESTSP
        |`BBTSP
        |`Albert
        |`Siepel ->
            ( UtlGrappa.inv_med medsov comoriarr1 comoriarr2 comoriarr3
            delimiter_lstlst circular)
    in
    let arr_med3 = from_ori_arr ori_arr_med3 in
    let comarr1 = from_ori_arr comoriarr1 in
    let comarr2 = from_ori_arr comoriarr2 in
    let comarr3 = from_ori_arr comoriarr3 in
    let totalcost1, (recost11,recost12), alied_arr11,alied_arr12  =
        cmp_cost `Breakinv arr1 arr_med3 comarr1 gen_cost_mat gapcode re_meth circular orientation
    in
    (* debug msg 
    Printf.printf "GenAli.ml: delimiters :\n%!";
    Array.iter (Printf.printf "%d,") delimiter_arr; Printf.printf "\n%!";
    Printf.printf "cost1 between [%!";  Array.iter (Printf.printf "%d,") arr1; 
    Printf.printf " ] and [ %!";  Array.iter (Printf.printf "%d," ) arr_med3;
    Printf.printf " ] and [ %!";  Array.iter (Printf.printf "%d," ) comarr1;
    Printf.printf " ] is %d\n%!" totalcost1;
    debug msg *)
    let totalcost2, (recost21,recost22), alied_arr21, alied_arr22 =
        cmp_cost `Breakinv arr2 arr_med3 comarr2 gen_cost_mat gapcode re_meth circular orientation
    in
    let totalcost3, (recost31,recost32), alied_arr31, alied_arr32 =
        cmp_cost `Breakinv arr3 arr_med3 comarr3 gen_cost_mat gapcode re_meth circular orientation
    in
    (Sequence.of_array arr_med3), (Array.to_list delimiter_arr),
    totalcost1,totalcost2,totalcost3,  
    recost11+recost12,recost21+recost22,recost31+recost32,
    alied_arr11,alied_arr21,alied_arr31,
    alied_arr12,alied_arr22,alied_arr32


(** [create_gen_ali_code state seq1 seq2 gen_cost_mat gen_gap_code 
*        re_meth max_swap_med circular] creates the general 
* alignment between [seq1] and [seq2] with minimum total cost
* where total cost = editing cost + rearrangement cost *)
let create_gen_ali_code kept_wag state (seq1 : int array) (seq2 : int array) 
        (gen_cost_mat : int array array) gen_gap_code re_meth max_swap_med circular orientation =   
(*debug msg*)
    Printf.printf "create_gen_ali_code, seq1 and seq2 = \n%!"; 
    Utl.printIntArr seq1; Utl.printIntArr seq2;
    Printf.printf "check matrix \n%!";
    Array.iter (fun code1 ->
        Array.iter (fun code2 ->
               Printf.printf "%10i" gen_cost_mat.(code1).(code2); )seq2;
        Printf.printf "%10i" gen_cost_mat.(code1).(gen_gap_code);
        Printf.printf "\n%!"; )seq1;
    Array.iter (fun code2 -> Printf.printf "%10i" gen_cost_mat.(code2).(gen_gap_code);
    )seq2; Printf.printf "\n%!";
(*debug msg*)
    let size = Array.length gen_cost_mat in 
    let gen_cost_mat = Array.init (size - 1) 
        (fun i -> Array.init (size - 1) (fun j -> gen_cost_mat.(i + 1).(j + 1))) 
    in 
    let gen_cost_ls = List.map (fun arr -> Array.to_list arr) (Array.to_list gen_cost_mat) in       
    let gen_cost_mat = Cost_matrix.Two_D.of_list ~use_comb:false gen_cost_ls (-1) in
    Cost_matrix.Two_D.set_gap gen_cost_mat gen_gap_code; 
    let wag_seq2 = find_wagner_ali kept_wag state seq1 seq2 gen_cost_mat 
        gen_gap_code re_meth circular orientation
    in 
    let init_cost, recost, alied_seq1, alied_seq2 =  
        cmp_cost state seq1 seq2 wag_seq2 gen_cost_mat gen_gap_code re_meth circular orientation
    in 
    let _, best_seq2 = 
        match max_swap_med with 
        | 0 -> init_cost, wag_seq2
        | _ ->
              multi_swap_locus state seq1 seq2 wag_seq2 init_cost  
                  gen_cost_mat gen_gap_code re_meth max_swap_med circular orientation 0  
    in   
    let final_cost, recost, alied_seq1, alied_seq2 =   
        cmp_cost state seq1 seq2 best_seq2 gen_cost_mat 
            gen_gap_code re_meth circular orientation
    in   
    final_cost, recost, alied_seq1, alied_seq2  


(**[cmp_cost3 seq1 seq2 seq3 med cost_mat gap re_meth cir sym] returns
* the total cost between [med] and three sequences [seq1], [seq2], [seq3] *)

let cmp_cost3 kept_wag seq1 seq2 seq3 med cost_mat gap re_meth cir orientation sym = 
    let max_swap_med = 1 in 
    let cmp_cost2 code1_arr code2_arr = 
        match sym with 
        | true ->
              let cost1, _, _, _ = create_gen_ali_code kept_wag `Breakinv code1_arr code2_arr
                  cost_mat gap re_meth max_swap_med cir orientation
              in 
              let cost2, _, _, _ = create_gen_ali_code kept_wag `Breakinv code2_arr code1_arr
                  cost_mat gap re_meth max_swap_med cir orientation
              in 
              min cost1 cost2
        | false ->
              let seq1 = Sequence.init (fun id -> code1_arr.(id)) (Array.length code1_arr) in 
              let seq2 = Sequence.init (fun id -> code2_arr.(id)) (Array.length code2_arr) in 
              let cost, _, _, _ = 
                  if Sequence.compare seq1 seq2 < 0 then 
                      create_gen_ali_code kept_wag `Breakinv code1_arr code2_arr 
                          cost_mat gap re_meth max_swap_med cir orientation
                  else 
                      create_gen_ali_code kept_wag  `Breakinv code2_arr code1_arr 
                          cost_mat gap re_meth max_swap_med cir orientation
              in 
              cost
    in 
    let cost1:int = cmp_cost2 seq1 med in 
    let cost2:int = cmp_cost2 seq2 med in 
    let cost3:int = cmp_cost2 seq3 med in 
    (cost1,cost2,cost3) 


(** [find_wagner_ali3 seq1 seq2 seq3 gen_cost_mat gap re_meth circular sym]
* finds the best median sequence of [seq1], [seq2] and [seq3] according
* to wagner-based algorithm *)
let find_wagner_ali3  kept_wag seq1 seq2 seq3 gen_cost_mat gap re_meth circular orientation sym = 
    let max_code = max (max (Utl.max_arr seq1) (Utl.max_arr seq2)) 
        (Utl.max_arr seq3) 
    in         
    let rec add (best_wagner : int array) cost code =
        if code > max_code then best_wagner, cost
        else begin
            let min_cost = ref cost in
            let new_best_wagner = ref best_wagner in 
            for p = 0 to Array.length best_wagner do
                let check acode =
                    let wagner = Utl.insert best_wagner p acode in 
                    let cost1,cost2,cost3 = cmp_cost3  kept_wag  seq1 seq2 seq3 wagner gen_cost_mat
                        gap re_meth circular orientation sym
                    in
                    let cost = cost1+cost2+cost3 in
                    if cost < !min_cost then begin
                        min_cost := cost;
                        new_best_wagner := wagner;
                    end 
                in 
                check code;
                check (code + 1);
            done;
            add !new_best_wagner !min_cost (code + 2)
        end 
    in  
    let init_cost1, init_cost2, init_cost3 = cmp_cost3  kept_wag seq1 seq2 seq3 [||] gen_cost_mat gap re_meth
        circular  orientation sym in 
    let init_cost = init_cost1+init_cost2+init_cost3 in
    let best_wagner, cost = add [||] init_cost 1 in 
    
    best_wagner, cost

(** [swap3 seq1 seq2 seq3 med gen_cost_mat gap re_meth circular sym]
    swaps the best found [med] to improve its quality *)
let rec swap3  kept_wag seq1 seq2 seq3 med gen_cost_mat gap re_meth circular orientation sym = 
    let best_cost1,best_cost2,best_cost3  =  (cmp_cost3 kept_wag seq1 seq2 seq3 med gen_cost_mat gap re_meth circular orientation sym) in 
    let best_cost1 = ref best_cost1 in
    let best_cost2 = ref best_cost2 in
    let best_cost3 = ref best_cost3 in
    let best_med = ref med in 
    let continue = ref false in 
    let len = Array.length med in 
    for p1 = 0 to len - 2 do
        for p2 = p1 + 1 to len - 1 do
            let new_med = Utl.swap_item p1 p2 med in
            let new_cost1,new_cost2,new_cost3 = cmp_cost3  kept_wag  seq1 seq2 seq3 new_med 
                gen_cost_mat gap re_meth circular orientation sym in 
            let new_cost = new_cost1+new_cost2+new_cost3 in
            let best_cost = !best_cost1 + !best_cost2 + !best_cost3 in
            if new_cost < best_cost then begin
                best_cost1 := new_cost1;
                best_cost2 := new_cost2;
                best_cost3 := new_cost3;
                best_med := new_med;
                continue := true;
            end 
        done
    done;
    if !continue then swap3 kept_wag seq1 seq2 seq3 !best_med gen_cost_mat gap re_meth
        circular orientation sym
    else  !best_med, !best_cost1, !best_cost2, !best_cost3


(** [create_gen_ali3 seq1 seq2 seq3 med gen_cost_mat 
*     alpha re_meth  max_swap_med circular sym] creates
* the general alignment among [seq1], [seq2], and [seq3] 
* such that total cost = editing cost + rearrangement cost is minimized *)
let create_gen_ali3  kept_wag  (seq1 : Sequence.s) (seq2 : Sequence.s) (seq3 : Sequence.s)
        (med : Sequence.s) gen_cost_mat alpha re_meth  max_swap_med circular orientation sym =

    let gap = Alphabet.get_gap alpha in 
    let seq1 : int array = Sequence.to_array seq1 in 
    let seq2 : int array = Sequence.to_array seq2 in 
    let seq3 : int array = Sequence.to_array seq3 in 
    let med : int array = Sequence.to_array med in 
    let med, cost1,cost2,cost3 = swap3 kept_wag seq1 seq2 seq3 med gen_cost_mat gap re_meth
        circular orientation sym in 
    let med_len = Array.length med in 
    let med_seq = Sequence.init (fun idx -> med.(idx)) med_len in
    med_seq, cost1,cost2,cost3 



(** [create_gen_ali_code3 state seq1 seq2 seq3 med gen_cost_mat 
        gen_gap_code re_meth max_swap_med circular sym] creates
* the general alignment among [seq1], [seq2], and [seq3] 
* such that total cost = editing cost + rearrangement cost is minimized *)
let create_gen_ali_code3 kept_wag state (seq1 : int array) (seq2 : int array) 
        (seq3 : int array) (med : int array) 
        (gen_cost_mat : int array array) gen_gap_code 
        re_meth max_swap_med circular orientation sym =
(*
    let med, cost = find_wagner_ali3 seq1 seq2 seq3 gen_cost_mat 
        gen_gap_code re_meth circular sym
    in 
*)
    let med, cost1,cost2,cost3 = swap3 kept_wag seq1 seq2 seq3 med gen_cost_mat gen_gap_code re_meth
        circular orientation sym in 
    med, cost1+cost2+cost3
