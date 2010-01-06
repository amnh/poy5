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

(** [cmp_recost state seq1 seq2 reseq2 re_meth circular] returns
* the rearrangement distance between two sequence [seq1] and [seq2] *)
let cmp_recost state seq1 seq2 reseq2 re_meth circular orientation = 
    let seq1, seq2, reseq2 = match orientation with
    | true -> to_ori_arr seq1, to_ori_arr seq2, to_ori_arr reseq2 
            (*(Array.map get_orientated_code seq1),
                    (Array.map get_orientated_code seq2),
                    (Array.map get_orientated_code reseq2)*)
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
            match re_meth with 
            | `Locus_Inversion cost -> 
                  (UtlGrappa.cmp_inversion_dis seq2 reseq2 circular) * cost  
            | `Locus_Breakpoint cost -> begin               
                  (UtlGrappa.cmp_oriented_breakpoint_dis seq2 reseq2 circular) * cost   
            end;
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
                      let partial_seq2 = Utl.insert best_wagner_seq2_arr.(w) pos code2 in 
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
              add best_w_arr added_seq2_ls tl
    in
    let wagner_seq2_arr = add [|[||]|] [] (Array.to_list seq2) in  
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




(** [create_gen_ali state seq1 seq1 gen_cost_mat alpha re_meth max_swap_med circular]
* creates the general alignment between [seq1] and [seq2] with minimum total cost 
* where total cost = editing cost + rearrangement cost *)
let create_gen_ali kept_wag state (seq1 : Sequence.s) (seq2 : Sequence.s) 
        (gen_cost_mat : Cost_matrix.Two_D.m) alpha re_meth 
        max_swap_med circular orientation =
(*debug msg  
    Printf.printf " create_gen_ali: seq1/seq2=\n%!";
    Sequence.printseqcode seq1; Sequence.printseqcode seq2;
  debug msg*)
    let gap = Alphabet.get_gap alpha in 
    let seq1 = Sequence.to_array seq1 in 
    let seq2 = Sequence.to_array seq2 in
    let equal_content =
        let oriarr1 = to_ori_arr (seq1)
        and oriarr2 = to_ori_arr (seq2) in
        let common_ori_arr1,common_ori_arr2 = Utl.get_common oriarr1 oriarr2
        equal_orientation in
        if (Array.length common_ori_arr1)=(Array.length seq1) then true
        else false
    in
    let wag_seq2 = 
        match equal_content with
        | false ->
        find_wagner_ali kept_wag state seq1 seq2 gen_cost_mat gap re_meth circular orientation
        | true -> seq2
    in 
 (*  debug msg 
       Printf.printf "wag_seq2= [ %!"; 
      Array.iter (Printf.printf "%d,") wag_seq2; Printf.printf " ] \n%!"; 
   debug msg*)
    let init_cost, recost, alied_seq1, alied_seq2 = 
        cmp_cost state seq1 seq2 wag_seq2 gen_cost_mat gap re_meth circular orientation
    in 
    let _, best_seq2 = 
        match max_swap_med with 
        | 0 -> init_cost, wag_seq2
        | _ ->
             (
                match equal_content with
                |true -> init_cost, seq2
                |false ->
                        multi_swap_locus state seq1 seq2 wag_seq2 init_cost  
                        gen_cost_mat gap re_meth max_swap_med circular orientation 0
             )
    in   
    let final_cost, (recost1,recost2), alied_seq1, alied_seq2 =  
        cmp_cost state seq1 seq2 best_seq2 gen_cost_mat gap re_meth circular orientation
    in   
    let alied_seq1 = Sequence.of_array alied_seq1 in
    let alied_seq2 = Sequence.of_array alied_seq2 in
(*debug msg 
   let (recost1,recost2) = recost in
   Printf.printf "cost=%d, recost1/recost2=%d/%d, alied_seq1/alied_seq2=\n%!"
       final_cost recost1 recost2;
    Sequence.printseqcode alied_seq1; Sequence.printseqcode alied_seq2;
 debug msg*)

    match equal_content with
    | true ->
          recost1, (recost1,recost2), alied_seq1,alied_seq2
    | false ->
           final_cost, (recost1,recost2), alied_seq1, alied_seq2  

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
    Sequence.s) (seq3 : Sequence.s) (gen_cost_mat: Cost_matrix.Two_D.m) alpha re_meth  max_swap_med circular orientation sym =
    (* debug msg 
    Printf.printf "create_gen_ali3_by_medsov, seq1,seq2,seq3=\n%!";
    Sequence.printseqcode seq1; Sequence.printseqcode seq2;
    Sequence.printseqcode seq3;
    debug msg *)
    let gapcode = Alphabet.get_gap alpha in
 (*   let size = Array.length gen_cost_mat in
    let gen_cost_mat = Array.init (size - 1) 
        (fun i -> Array.init (size - 1) (fun j -> gen_cost_mat.(i + 1).(j + 1))) 
    in 
  
    let gen_cost_ls = List.map (fun arr -> Array.to_list arr) (Array.to_list gen_cost_mat) in  
    let gen_cost_mat = Cost_matrix.Two_D.of_list ~use_comb:false gen_cost_ls (-1) in
    Cost_matrix.Two_D.set_gap gen_cost_mat gapcode; 
 *)
    let arr1 = Sequence.to_array seq1
    and arr2 = Sequence.to_array seq2 
    and arr3 = Sequence.to_array seq3 in
    let oriarr1 = to_ori_arr (arr1)
    and oriarr2 = to_ori_arr (arr2) 
    and oriarr3 = to_ori_arr (arr3) in
    let comoriarr1,comoriarr2,comoriarr3 = Utl.get_common3 oriarr1 oriarr2 oriarr3 equal_orientation in
    let ori_arr_med3 =
        match medsov with
        |`Vinh ->
                failwith "Vinh median solver is not in grappa"
        |`SimpleLK
        |`ChainedLK 
        |`COALESTSP
        |`BBTSP
        |`Albert
        |`Siepel ->
            ( UtlGrappa.inv_med medsov comoriarr1 comoriarr2 comoriarr3 circular)

    in
    let arr_med3 = from_ori_arr ori_arr_med3 in
    let comarr1 = from_ori_arr comoriarr1 in
    let comarr2 = from_ori_arr comoriarr2 in
    let comarr3 = from_ori_arr comoriarr3 in
    let totalcost1, (recost11,recost12), alied_arr11,alied_arr12  =
        cmp_cost `Breakinv arr1 arr_med3 comarr1 gen_cost_mat gapcode re_meth circular orientation
    in
    (* debug msg
    Printf.printf "~~ cost1 between [%!"; 
    Array.iter (Printf.printf "%d,") arr1; 
    Printf.printf " ] and [ %!";
    Array.iter (Printf.printf "%d," ) arr_med3;
    Printf.printf " ] and [ %!";
    Array.iter (Printf.printf "%d," ) comarr1;
    Printf.printf " ] is %d\n%!" totalcost1;
    debug msg *)
    let totalcost2, (recost21,recost22), alied_arr21, alied_arr22 =
        cmp_cost `Breakinv arr2 arr_med3 comarr2 gen_cost_mat gapcode re_meth circular orientation
    in
    (* debug msg
    Printf.printf "~~ cost2 between [%!"; 
    Array.iter (Printf.printf "%d,") arr2; 
    Printf.printf " ] and [ %!";
    Array.iter (Printf.printf "%d," ) arr_med3;
    Printf.printf " ] and [ %!";
    Array.iter (Printf.printf "%d," ) comarr2;
    Printf.printf " ] is %d\n%!" totalcost2;
     debug msg *)
    let totalcost3, (recost31,recost32), alied_arr31, alied_arr32 =
        cmp_cost `Breakinv arr3 arr_med3 comarr3 gen_cost_mat gapcode re_meth circular orientation
    in
    (* debug msg
    Printf.printf "~~ cost3 between [%!"; 
    Array.iter (Printf.printf "%d,") arr3; 
    Printf.printf " ] and [ %!";
    Array.iter (Printf.printf "%d," ) arr_med3;
    Printf.printf " ] and [ %!";
    Array.iter (Printf.printf "%d," ) comarr3;
    Printf.printf " ] is %d\n%!" totalcost3;
    debug msg*)
    (Sequence.of_array arr_med3), 
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
