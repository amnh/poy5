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

let fprintf = Printf.fprintf;;

type dyna_state_t = Data.dyna_state_t


let orientation code = 
    if code mod 2 = 0 then -(code / 2)
    else (code + 1) /2

let equal_orientation code1 code2 = compare (abs code1) (abs code2) 

(** Given two arrays of orders [seq2] and [re_seq2], 
 * compute their rearrangement distance*)
let cmp_recost state seq1 seq2 reseq2 re_meth circular =     
    let seq1 = Array.map orientation seq1 in 
    let seq2 = Array.map orientation seq2 in 
    let reseq2 = Array.map orientation reseq2 in 

    if Array.length seq2 = 0 then 0, 0
    else begin
        let recost1 = 
            match state with
            | `Breakinv ->                            
                  
                  let com_seq1_arr, com_reseq2_arr = Utl.get_common seq1 reseq2 equal_orientation in 

                  (match re_meth with 
                  | `Inversion cost -> 
                        (UtlGrappa.cmp_inversion_dis com_seq1_arr com_reseq2_arr circular) * cost  
                  | `Breakpoint cost ->                                       
                        (UtlGrappa.cmp_oriented_breakpoint_dis com_seq1_arr com_reseq2_arr circular) * cost)                     
            | _ -> 0
        in 

        let recost2  = 
            match re_meth with 
            | `Inversion cost -> 
                  (UtlGrappa.cmp_inversion_dis seq2 reseq2 circular) * cost  
            | `Breakpoint cost ->               
                  (UtlGrappa.cmp_oriented_breakpoint_dis seq2 reseq2 circular) * cost   
        in  

        recost1, recost2
    end 


(** Given two sequences [seq1], [seq2] and rearranged sequence 
 * [reseq2] of [seq2], compute the total cost between [seq1] and [reseq2] where 
 * total cost = editing cost (seq1, reseq2) plus rearrangement cost (seq2, reseq2) *)
let cmp_cost state code1_arr code2_arr recode2_arr 
        (cost_mat : Cost_matrix.Two_D.m) gap re_meth circular = 

    let seq1 = Sequence.init (fun idx -> code1_arr.(idx)) (Array.length code1_arr) in
    let reseq2 = Sequence.init (fun idx -> recode2_arr.(idx)) (Array.length recode2_arr) in


    let alied_seq1, alied_reseq2, editing_cost =  
        Sequence.Align.align_2 ~first_gap:false seq1 reseq2 cost_mat
            Matrix.default  
    in   

    let alied_code1_arr = Sequence.to_array alied_seq1 in 
    let alied_recode2_arr = Sequence.to_array alied_reseq2 in 

    let recost1, recost2 = cmp_recost state code1_arr code2_arr recode2_arr re_meth circular in 
   (editing_cost + recost1 + recost2), (recost1, recost2), alied_code1_arr, alied_recode2_arr



(** Given two sequences [seq1], [seq2], find rearranged sequence [reseq2]
 * of sequence [seq2] using stepwise addition method 
 * such that the total cost is minimum where 
 * total cost = editing cost (seq1, reseq2) plus rearrangement cost (seq2, reseq2) *)
let find_wagner_ali state seq1 seq2 gen_cost_mat gap re_meth circular = 
    let rec add (best_wagner_seq2 : int array) added_seq2_ls rem_seq2_ls = 
        match rem_seq2_ls with
        | [] -> best_wagner_seq2 
        | code2 :: tl ->
              let added_seq2_ls = added_seq2_ls @ [code2] in 

              let wagner_cost = ref Utl.large_int in 
              let wagner_seq2 = ref [||] in 
              let update partial_seq2 =
                  let cost, (_, _), _, _  = 
                      cmp_cost state seq1 (Array.of_list added_seq2_ls) partial_seq2 
                          gen_cost_mat gap re_meth circular
                  in

                  if cost < !wagner_cost then begin
                      wagner_cost := cost; 
                      wagner_seq2 := partial_seq2 
                  end;

              in 

              for pos = 0 to Array.length best_wagner_seq2 do 
                  let partial_seq2 = Utl.insert best_wagner_seq2 pos code2 in 
                  update partial_seq2;

                  if code2 mod 2 = 1 then begin
                      let partial_seq2 = Utl.insert best_wagner_seq2 pos (code2 + 1) in 
                      update partial_seq2;
                  end 
              done;  
              add !wagner_seq2 added_seq2_ls tl
    in
              
    let wagner_seq2 = add [||] [] (Array.to_list seq2) in  
    wagner_seq2




(** Given two sequences [seq1], [seq2] and rearranged sequence 
 * [reseq2] of [seq2], swap [reseq2] in order to minimize the total
 * cost between [seq1] and [reseq2] where 
 * total cost = editing cost (seq1, reseq2) plus rearrangement cost (seq2, reseq2) *)
let rec multi_swap_locus state seq1 seq2 best_seq2 best_cost 
        gen_cost_mat gap re_meth max_swap_med
        circular num_done_swap =             

    let len2 = Array.length best_seq2 in  
    let swap_ls = ref [] in 

    for donor_pos = 0 to len2 - 2 do
        for rev_pos = donor_pos + 1 to len2 - 1 do 
            let new_seq2 = Utl.swap_item donor_pos rev_pos best_seq2 in 

            let new_cost, (_, _), _, _ = cmp_cost state seq1 seq2 new_seq2 
                gen_cost_mat gap re_meth circular
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
            gen_cost_mat gap re_meth circular
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
                re_meth max_swap_med circular (num_done_swap + 1)
    end 

        
(** Given two sequences [seq1] and [seq2] of general character, create 
 * the general alignment between [seq1] and [seq2] with minimum total cost 
 * where total cost = editing cost + rearrangement cost *)
let create_gen_ali state (seq1 : Sequence.s) (seq2 : Sequence.s) 
        (gen_cost_mat : Cost_matrix.Two_D.m) alpha re_meth 
        max_swap_med circular =

    let gap = Alphabet.get_gap alpha in 
    let seq1 : int array = Sequence.to_array seq1 in 
    let seq2 : int array = Sequence.to_array seq2 in 

    let wag_seq2 = find_wagner_ali state seq1 seq2 gen_cost_mat 
        gap re_meth circular
    in 

    let init_cost, recost, alied_seq1, alied_seq2 = 
        cmp_cost state seq1 seq2 wag_seq2 gen_cost_mat gap re_meth circular
    in 

    let _, best_seq2 = 
        match max_swap_med with 
        | 0 -> init_cost, wag_seq2
        | _ -> 
              multi_swap_locus state seq1 seq2 wag_seq2 init_cost  
                  gen_cost_mat gap re_meth max_swap_med circular 0
    in   

    let final_cost, recost, alied_seq1, alied_seq2 =   
        cmp_cost state seq1 seq2 best_seq2 gen_cost_mat gap re_meth circular 
    in   

    let ali_len = Array.length alied_seq1 in 
    let alied_seq1 = Sequence.init (fun idx -> alied_seq1.(idx)) ali_len in
    let alied_seq2 = Sequence.init (fun idx -> alied_seq2.(idx)) ali_len in
    final_cost, recost, alied_seq1, alied_seq2  





(** Given two sequence [seq1] and [seq2] of character codes, create 
 * the general alignment between [seq1] and [seq2] with minimum total cost 
 * where total cost = editing cost + rearrangement cost *)
let create_gen_ali_code state (seq1 : int array) (seq2 : int array) 
        (gen_cost_mat : int array array) gen_gap_code 
        re_meth max_swap_med circular =

    let size = Array.length gen_cost_mat in 
    let gen_cost_mat = Array.init (size - 1) 
        (fun i -> Array.init (size - 1) (fun j -> gen_cost_mat.(i + 1).(j + 1))) 
    in 

    let gen_cost_ls = List.map (fun arr -> Array.to_list arr) (Array.to_list gen_cost_mat) in       
    let gen_cost_mat = Cost_matrix.Two_D.of_list ~use_comb:false gen_cost_ls in
    Cost_matrix.Two_D.set_gap gen_cost_mat gen_gap_code; 

    let wag_seq2 = find_wagner_ali state seq1 seq2 gen_cost_mat 
        gen_gap_code re_meth circular
    in 

    let init_cost, recost, alied_seq1, alied_seq2 =  
        cmp_cost state seq1 seq2 wag_seq2 gen_cost_mat gen_gap_code re_meth circular
    in 


    let _, best_seq2 = 
        match max_swap_med with 
        | 0 -> init_cost, wag_seq2
        | _ ->
              multi_swap_locus state seq1 seq2 wag_seq2 init_cost  
                  gen_cost_mat gen_gap_code re_meth max_swap_med circular 0  
    in   
    


    let final_cost, recost, alied_seq1, alied_seq2 =   
        cmp_cost state seq1 seq2 best_seq2 gen_cost_mat 
            gen_gap_code re_meth circular 
    in   

    final_cost, recost, alied_seq1, alied_seq2  





let cmp_cost3 seq1 seq2 seq3 med cost_mat gap re_meth cir sym = 
    let max_swap_med = 1 in 
    let cmp_cost2 code1_arr code2_arr = 
        match sym with 
        | true ->
              let cost1, _, _, _ = create_gen_ali_code `Breakinv code1_arr code2_arr
                  cost_mat gap re_meth max_swap_med cir
              in 


              let cost2, _, _, _ = create_gen_ali_code `Breakinv code2_arr code1_arr
                  cost_mat gap re_meth max_swap_med cir
              in 
              min cost1 cost2
        | false ->
              let seq1 = Sequence.init (fun id -> code1_arr.(id)) (Array.length code1_arr) in 
              let seq2 = Sequence.init (fun id -> code2_arr.(id)) (Array.length code2_arr) in 
              
              let cost, _, _, _ = 
                  if Sequence.compare seq1 seq2 < 0 then 
                      create_gen_ali_code `Breakinv code1_arr code2_arr 
                          cost_mat gap re_meth max_swap_med cir 
                  else 
                      create_gen_ali_code `Breakinv code2_arr code1_arr 
                          cost_mat gap re_meth max_swap_med cir 
              in 
              cost
    in 

    let cost1 = cmp_cost2 seq1 med in 
    let cost2 = cmp_cost2 seq2 med in 
    let cost3 = cmp_cost2 seq3 med in 
(*    fprintf stdout "Cost3: %i %i %i\n" cost1 cost2 cost3;*)
    cost1 + cost2 + cost3
    


(** Given three sequences [seq1], [seq2], [seq3], find the best median of these
 * three sequences *)
let find_wagner_ali3 seq1 seq2 seq3 gen_cost_mat gap re_meth circular sym = 
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
                    let cost = cmp_cost3 seq1 seq2 seq3 wagner gen_cost_mat
                        gap re_meth circular sym
                    in
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
    let init_cost = cmp_cost3 seq1 seq2 seq3 [||] gen_cost_mat gap re_meth
        circular sym in 
    let best_wagner, cost = add [||] init_cost 1 in 
    
    best_wagner, cost


(** Given three sequences [seq1], [seq2], [seq3], find the best median of these
 * three sequences *)
let rec swap3 seq1 seq2 seq3 med gen_cost_mat gap re_meth circular sym = 
    let best_cost = ref (cmp_cost3 seq1 seq2 seq3 med gen_cost_mat gap re_meth
                             circular sym) in 
    let best_med = ref med in 
    let continue = ref false in 

    let len = Array.length med in 
    for p1 = 0 to len - 2 do
        for p2 = p1 + 1 to len - 1 do
            let new_med = Utl.swap_item p1 p2 med in
            let new_cost = cmp_cost3 seq1 seq2 seq3 new_med 
                gen_cost_mat gap re_meth circular sym in 

            if new_cost < !best_cost then begin
                best_cost := new_cost;
                best_med := new_med;
                continue := true;
            end 
        done
    done;

    if !continue then swap3 seq1 seq2 seq3 !best_med gen_cost_mat gap re_meth
        circular sym
    else  !best_med, !best_cost


(** Given two sequences [seq1] and [seq2] of general character, create 
 * the general alignment between [seq1] and [seq2] with minimum total cost 
 * where total cost = editing cost + rearrangement cost *)
let create_gen_ali3  (seq1 : Sequence.s) (seq2 : Sequence.s) (seq3 : Sequence.s)
        (med : Sequence.s) 
        gen_cost_mat alpha re_meth  max_swap_med circular sym =

    let gap = Alphabet.get_gap alpha in 
    let seq1 : int array = Sequence.to_array seq1 in 
    let seq2 : int array = Sequence.to_array seq2 in 
    let seq3 : int array = Sequence.to_array seq3 in 
    let med : int array = Sequence.to_array med in 

(*
    let med, _  = find_wagner_ali3 seq1 seq2 seq3 gen_cost_mat 
        gap re_meth circular sym
    in 
*)
    let med, cost = swap3 seq1 seq2 seq3 med gen_cost_mat gap re_meth
        circular sym in 



    let med_len = Array.length med in 
    let med_seq = Sequence.init (fun idx -> med.(idx)) med_len in

    med_seq, cost 



(** Given two sequence [seq1] and [seq2] [seq3] of character codes, create 
 * the general alignment between [seq1] and [seq2] [seq3] with minimum total cost 
 * where total cost = editing cost + rearrangement cost *)
let create_gen_ali_code3 state (seq1 : int array) (seq2 : int array) 
        (seq3 : int array) (med : int array) 
        (gen_cost_mat : int array array) gen_gap_code 
        re_meth max_swap_med circular sym =
(*
    let med, cost = find_wagner_ali3 seq1 seq2 seq3 gen_cost_mat 
        gen_gap_code re_meth circular sym
    in 
*)
    let med, cost = swap3 seq1 seq2 seq3 med gen_cost_mat gen_gap_code re_meth
        circular sym in 
    med, cost
