(* POY 5.0 Alpha. A phylogenetic analysis program using Dynamic Homologies.   *)
(* Copyright (C) 2011  Andrés Varón, Le Sy Vinh, Illya Bomash, Ward Wheeler,  *)
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

let () = SadmanOutput.register "Seed" "$Revision: 2782 $"
type pairChromPam_t = ChromPam.chromPairAliPam_t

(** The seed module contains data and funtions to 
* determine all seeds between two chromosomes *)

(** seed_t is a data structure presenting a segment of nucleotides
* which is in both chromosomes *)
type seed_t = {
    mutable id : int; (* Each seed is assigned an [id]*)
    mutable sta1 : int; (* Start position of the seed on the chromosome 1 *)
    mutable sta2 : int; (* Start position of the seed on the chromosome 2 *)
    mutable len : int; (* Length of the seed *)
    mutable score : int; (* Score of the seed. A better seed is longer and has higer score *)
    mutable is_chained : bool (* is_changed = true if this seed can be chained 
                                with another seed, otherwise false *)
}

type sufTree_t = SufTree.sufTree_t
type incList_t = seed_t IncList.incList_t

(** The number of memory cells we allocate each time to contain found seeds *)
let num_located_cell = 100000

let num_used_cell = ref num_located_cell

(** An array of memory cells which were allocated to contain seeds *)
let located_cell_arr : incList_t array ref = ref [||]

let fprintf = Printf.fprintf
let deref = Utl.deref


let copy (seed : seed_t) = {
    id = seed.id;
    sta1 = seed.sta1;
    sta2 = seed.sta2;
    len = seed.len;
    score = seed.score;
    is_chained = seed.is_chained;
}

let swap (seed : seed_t) = 
    let sta1 = seed.sta1 in 
    seed.sta1 <- seed.sta2;
    seed.sta2 <- sta1


let create_empty_seed () = {
    id = -1; 
    sta1 = -1; sta2 = -1; len = 0; 
    score = 0; is_chained = false}

(** This function returns a dummy seed of length 1
* which is used as a starting point for dynamic programming to 
* chain seeds together *)
let get_dum_first_seed (min_pos1 : int) (min_pos2 : int) = {
    id = -1; 
    sta1 = min_pos1 - 1; sta2 = min_pos2 - 1; 
    len = 1; score = 0; is_chained = false}

(** This function returns a dummy seed of length 1
* which is used as an end point for dynamic programming to 
* chain seeds together *)
let get_dum_last_seed (max_pos1 : int) (max_pos2 : int) = {
    id = -1; 
    sta1 = max_pos1 + 1; sta2 = max_pos2 + 1; len = 1; 
    score = 0; is_chained = false}

let set_seed (seed : seed_t) (new_sta1 : int) (new_sta2 : int) 
        (new_len : int) (new_score : int) (new_is_chained : bool) = 
    seed.sta1 <- new_sta1;
    seed.sta2 <- new_sta2;
    seed.len <- new_len;
    seed.score <- new_score;
    seed.is_chained <- new_is_chained


let create_new_seed (new_sta1 : int) (new_sta2 : int) 
        (new_len : int) (new_score : int) = {
    id = -1; sta1 = new_sta1; sta2 = new_sta2; 
    len = new_len; score = new_score; is_chained = false}

(** Given two seeds [s1] and [s2], this function 
* returns an integer number as the distance between 
* [s1] and [s2] on the first chromosome*)
let cmp_dis1 (s1 : seed_t) (s2 : seed_t) = 
    s2.sta1 - (s1.sta1 + s1.len - 1) - 1

(** Given two seeds [s1] and [s2], this function 
* returns an integer number as the distance between 
* [s1] and [s2] on the second chromosome*)
let cmp_dis2 (s1 : seed_t) (s2 : seed_t) = 
    s2.sta2 - (s1.sta2 + s1.len - 1) - 1

(** Given two seeds [s1] and [s2], this function 
* returns an integer number as the distance between 
* [s1] and [s2] which is the minimum distance
* on the both chromosomes *)
let cmp_dis (s1 : seed_t) (s2 : seed_t) = 
    min (cmp_dis1 s1 s2) (cmp_dis2 s1 s2)

(** Given two seeds [s1] and [s2], this function 
* returns an integer number as the diagonal distance between 
* [s1] and [s2] which is the distance difference on the
* first and second chromosomes *)
let cmp_dia_dis (s1 : seed_t) (s2 : seed_t) = 
    abs ( (cmp_dis1 s1 s2) - (cmp_dis2 s1 s2) )

(** Given a [seed], [min_pos2] and [max_pos2]
* which are start and end of the second chromosome, 
* this function inverts the [sta2] of the seed 
* due to the inversion of the second chromosome *)
let invert (min_pos2 : int) (max_pos2 : int) (seed : seed_t) = 
    let end_pos2 = seed.sta2 + seed.len - 1 in    
    seed.sta2 <- min_pos2 + (max_pos2 - end_pos2) 

let print (seed : seed_t) = 
    fprintf stdout "id: %i --> (%i, %i), len: %i, score: %i \n" 
        seed.id seed.sta1 seed.sta2 seed.len seed.score;
    flush stdout

let compare_end1 (seed1 : seed_t) (seed2 : seed_t) = 
    compare (seed1.sta1 + seed1.len) (seed2.sta1 + seed2.len)

let compare_start1 (seed1 : seed_t) (seed2 : seed_t) = 
    compare seed1.sta1 seed2.sta1 
    
(** Given two seeds [seed1], [seed2] and pairwise chromosome
* parameters [ali_pam], this function returns an integer number
* as the cost to connect two seeds *)
let cmp_connecting_cost (seed1 : seed_t) (seed2 : seed_t) 
        (ali_pam :  pairChromPam_t) = 
    let dis1  = seed2.sta1 - (seed1.sta1 + seed1.len - 1) - 1 in 
    let dis2  = seed2.sta2 - (seed1.sta2 + seed1.len - 1) - 1 in 
    let digDis = abs (dis1 - dis2) in
    match digDis with 
    | 0 -> 0
    | _ ->  ali_pam.ChromPam.gap_opening_cost + 
          (digDis - 1) * ali_pam.ChromPam.gap_ext_cost 
    
(** Given a [seed_ls] and pairwise chromosome parameters [ali_pam], 
* this function does dynamic programming to create local connections
* between seeds in order to chain seeds into blocks *)
let create_local_connection (seed_ls : seed_t list) (ali_pam : pairChromPam_t) =    
    let num_seed = List.length seed_ls in
    let sorted_end_seed_arr = Array.of_list seed_ls in                                                          
    Array.sort compare_end1 sorted_end_seed_arr;
    Array.iteri (fun index seed -> seed.id <- index) sorted_end_seed_arr;
        
    let seed_cost_arr = Array.init num_seed 
        (fun index -> sorted_end_seed_arr.(index).score * ali_pam.ChromPam.mat_cost) in

    let acc_cost_arr = Array.copy seed_cost_arr in 
    let back_arr = Array.make num_seed (-1) in 
    
    
    let rearranged_len = ChromPam.get_min_rearrangement_len ali_pam in 
    for seed_no = 0 to num_seed - 2 do
            
        let donor_seed = sorted_end_seed_arr.(seed_no) in
        let donor_cost = acc_cost_arr.(seed_no) in 
        let donor_end_pos1 = donor_seed.sta1 + donor_seed.len - 1 in 
        let donor_end_pos2 = donor_seed.sta2 + donor_seed.len - 1 in 

        let rec move (rec_seed_no : int) =
            let rec_seed = sorted_end_seed_arr.(rec_seed_no) in                 
            if rec_seed.sta1 - donor_end_pos1 - 1 <= rearranged_len then 
            begin               
                if (donor_end_pos1 < rec_seed.sta1) 
                    && (donor_end_pos2 < rec_seed.sta2) 
                    && (rec_seed.sta2 - donor_end_pos2 - 1 <= rearranged_len) then
                begin
                    let connecting_cost = 
                        cmp_connecting_cost donor_seed rec_seed ali_pam in                
                    if acc_cost_arr.(rec_seed_no) >= donor_cost + 
                        connecting_cost + seed_cost_arr.(rec_seed_no) then 
                    begin
                        acc_cost_arr.(rec_seed_no) <- donor_cost + 
                            connecting_cost + seed_cost_arr.(rec_seed_no);
                        back_arr.(rec_seed_no) <- seed_no
                    end
                end;
                
                if rec_seed_no + 1 < num_seed then 
                    move (rec_seed_no + 1)
            end
        in
        move (seed_no + 1)
    done;   
    acc_cost_arr, back_arr, sorted_end_seed_arr


(** This function allocate a bunch of [num_located_cell] memory cells
* to contain found seeds. The reason is that it reduces
* a signigicant amount of time compared to allocate
* one memory cell each time *)
let locate_mem () = 
    num_used_cell := 0;
    located_cell_arr := Array.init num_located_cell 
        (fun _ -> {IncList.data = create_empty_seed (); IncList.next = None})
    

(** Given [sta1], [sta2], [len], [is_chained],this function returns 
* a memory cell which contains a seed made from these information. 
* Note that the type of this memory cell is incList*)
let create_cell (sta1 : int) (sta2 : int) 
        (len : int) (is_chained : bool) = 
    (if (!num_used_cell = num_located_cell) then 
        locate_mem ()); 
    let new_cell = !located_cell_arr.(!num_used_cell) in    
    num_used_cell := !num_used_cell + 1;
    
    set_seed new_cell.IncList.data sta1 sta2 len len is_chained;      
    new_cell
    

(** Given
* - [sta2]:  A position in the second chromosome
* - [sta1_ls]: A list of positions in the first chromosomes
*              where the nucleotide at [sta2] in the second chromosome 
*              appears in the first chromosome), 
* - [evolving_seed_ls]: A nascent seed list
* - [ali_pam]: Pairwise chromosome parameters
* - [code1_arr], [code2_arr]: chromosomes on in array presentations
* This function update the [evolving_seed_ls] by prolonging 
* these evolving seeds if possible*)
let connect (sta2 : int) (sta1_ls : int list) 
    (evolving_seed_ls : incList_t) (ali_pam : pairChromPam_t) 
    (code1_arr : int array) (code2_arr : int array) = 
    
    let k = ChromPam.get_min_seed_len ali_pam in 
    let max_gap = ali_pam.ChromPam.max_gap in 

    let rec chain_insert (sta1 : int) 
            (new_cell : incList_t) (cell : incList_t) = 

        match cell.IncList.next with 
        | None -> (* last cell -> insert *)              
              cell.IncList.next <- Some (ref new_cell); 
        | Some next_cell ->  
            let next_seed = !next_cell.IncList.data in 
            match compare (sta1 + k - 1) (next_seed.sta1 + next_seed.len) with
            | 0 ->  
                if (sta2 + k - 1 = next_seed.sta2 + next_seed.len) then 
                begin
                    next_seed.len <- next_seed.len + 1;
                    cell.IncList.next <- !next_cell.IncList.next; (* delete and reinsert to keep the same increasing order *)
                    IncList.insert cell compare_end1 !next_cell;
                end 
                else chain_insert sta1 new_cell !next_cell
                    
            | cmp when cmp > 0 -> (* delete *)
                let end_pos1_seed = next_seed.sta1 + next_seed.len - 1 in 
                let end_pos2_seed = next_seed.sta2 + next_seed.len - 1 in 
                let dis1 = sta1 - end_pos1_seed - 1 in  
                let dis2 = sta2 - end_pos2_seed - 1 in 
                if (dis1 >= 0) && (dis2 >= 0) && (dis1 <= max_gap) then
                begin
                    next_seed.is_chained <- true;
                    new_cell.IncList.data.is_chained <- true;
                end;
                chain_insert sta1 new_cell !next_cell                     

            | _ -> (* cmp < 0 *)
                cell.IncList.next <- Some (ref new_cell);
                new_cell.IncList.next <- Some next_cell;      
    in

    let new_seed_ls = ref [] in     
    let is_reached_last_cell = ref false in
    let sig_k = ali_pam.ChromPam.sig_k in 

    let rec delete (sta1 : int) (cell : incList_t) = 
        match cell.IncList.next with
        | None ->  
              is_reached_last_cell := true;    
              cell 
        | Some next_cell -> 
              let next_seed = !next_cell.IncList.data in  
              if (sta1 + 1 - (next_seed.sta1 + next_seed.len - 1) - 1) > max_gap then begin
                  if (sta2  + 1 - (next_seed.sta2 +  next_seed.len - 1) - 1) > max_gap then  
                  begin 
                      cell.IncList.next <- !next_cell.IncList.next; 
                      next_seed.score <- next_seed.len; 
                        
                      (if next_seed.is_chained || next_seed.score >= sig_k then
                          new_seed_ls := next_seed::!new_seed_ls);                
                      delete sta1 cell 
                  end 
                  else delete sta1 !next_cell 
              end
              else cell 
    in              
    

    let last_cell = ref evolving_seed_ls in     
    let rec move_forward (sta1 : int) =       
        let new_cell = create_cell sta1 sta2  k false in            
        chain_insert sta1 new_cell !last_cell;
        last_cell := delete sta1 !last_cell;
    in
    List.iter move_forward sta1_ls;
    
    (if !is_reached_last_cell = false then 
        last_cell := delete max_int !last_cell); 

    !new_seed_ls




(** Given a suffix tree [suf_tree] which encodes 
* the first chromosome, the second chromsome [chrom2],
* pairwise chromosome parameters [ali_pam], this function
* returns a list of positive seeds *)
let determine_pos_seed (suf_tree : sufTree_t) (chrom2 : Sequence.s) ali_pam =        
(* First, create the first k_mer from position 0..(k_mer - 1) *)
    let code1_arr = suf_tree.SufTree.code_arr in 

    let code2_arr = Sequence.to_array chrom2 in 
    let code2_arr = Array.map (fun code -> code land 15) code2_arr in

    let k = ChromPam.get_min_seed_len ali_pam in         
    let seed_ls = ref [] in 
    let root = suf_tree.SufTree.root in 
    
    let start_seed = 
        {id = -1; sta1 = -1; sta2 = -1; 
         len = 0; score = 0; is_chained = false} 
    in
    let evolving_seed_ls = {IncList.data = start_seed; IncList.next = None} in

    let full_major_size = 2000000 in 
    (if (ali_pam.ChromPam.max_pos2 - ali_pam.ChromPam.min_pos2 > full_major_size) then Gc.full_major () );


    for sta2 = ali_pam.ChromPam.min_pos2 to (ali_pam.ChromPam.max_pos2 - k + 1) do
                                
        let k_mer = Array.to_list (Array.sub code2_arr sta2 k) in   
        let stop_node, res_k_mer = SufNode.find_stop_node ~exact:true root code1_arr k_mer in
        let sta1_ls = 
            if (res_k_mer = []) then stop_node.SufNode.leaf_ls
            else []
        in


        let new_seed_ls = connect sta2 sta1_ls evolving_seed_ls
            ali_pam code1_arr code2_arr in 
        
        
        seed_ls := new_seed_ls @ !seed_ls;

    done; (* End of position sta2_*)  

    (* Get the remaining seeds *)
    let new_seed_ls = connect max_int [] evolving_seed_ls 
        ali_pam code1_arr code2_arr in

    seed_ls := new_seed_ls @ !seed_ls;

    !seed_ls


(** Given a suffix tree [suf_tree] which encodes 
* the first chromosome, the second chromsome [chrom2],
* pairwise chromosome parameters [ali_pam], this function
* returns a list of negative seeds *)
let determine_neg_seed (suf_tree : sufTree_t) (chrom2 : Sequence.s) 
        (ali_pam : pairChromPam_t) = 
    let com_chrom2 = Sequence.complement_chrom Alphabet.nucleotides chrom2 in 
    let neg_seed_ls = determine_pos_seed suf_tree com_chrom2 ali_pam in 
    neg_seed_ls  

    
let print_sta (seed_ls : seed_t list) = 
    let max_len = ref 0 in 
    let find_max_len (seed : seed_t) = max_len := max !max_len seed.len in 
    List.iter find_max_len seed_ls;
    
    let frq_arr = Array.make (!max_len + 1) 0 in 
    List.iter (fun seed -> frq_arr.(seed.len) <- 
                   frq_arr.( seed.len) + 1) seed_ls;

    for len = 0 to !max_len do
        if frq_arr.(len) > 0 then    
            fprintf stdout "%i -> %i \n" len frq_arr.(len)
    done;   

    let num_chained_seed = ref 0 in 
    List.iter (fun seed -> if seed.is_chained then 
                   num_chained_seed := !num_chained_seed + 1) seed_ls;

    fprintf stdout "Number seeds: %i, number chained seds: %i\n" 
        (List.length seed_ls) !num_chained_seed

    

(** Given two chromosomes [chrom1] and [chrom2],
* pairwise chromosome parameters [ali_pam], seed direction [direction],
* this function returns two lists: positive seed list and negative seed_ls.
* Note that:
* - if [direction] is positive, negative seed list is empty,
* - if [direction] is negative, postive seed list is empty *)
let determine_seed (chrom1 : Sequence.s) (chrom2 : Sequence.s) 
        (ali_pam : pairChromPam_t) (direction : ChromPam.direction_t) = 


    let code1_arr = Sequence.to_array chrom1 in
    let code1_arr = Array.map (fun code -> code land 15) code1_arr in 

    let k_mer_tree = SufTree.create_k_mer_tree code1_arr 
        (ChromPam.get_min_seed_len ali_pam) in     
    
    let pos_seed_ls = 
        match direction with        
        | `Positive | `BothDir -> 
              determine_pos_seed k_mer_tree chrom2 ali_pam 
        | `Negative -> []
    in

    let neg_seed_ls = 
        match direction, ali_pam.ChromPam.negative with        
        | `Negative, true | `BothDir, true ->  
              determine_neg_seed k_mer_tree chrom2 ali_pam 
        | _, _ -> [] 
    in  
    pos_seed_ls, neg_seed_ls



(** Given a [seed], two chromosomes [chrom1] and [chrom2], 
* transformation cost matrix [cost_mat], this function 
* returns a pairwise sequence alignment coresponding to [seed] *)
let create_alied_seed (seed : seed_t) (chrom1 : Sequence.s) 
        (chrom2 : Sequence.s) (cost_mat : Cost_matrix.Two_D.m) =
    match seed.len > 1 with
    | true ->
        let subseq1 = Sequence.sub chrom1 seed.sta1 seed.len in 
        let subseq2 = Sequence.sub chrom2 seed.sta2 seed.len in         
        let cost = Sequence.cmp_ali_cost subseq1 subseq2 `Positive cost_mat in 
        subseq1, subseq2, cost
    | false -> (* dummy seed*)
        Sequence.get_empty_seq (), Sequence.get_empty_seq (), 0
    

(** Given two seeds [seed1] and [seed2], chromosomes [chrom1], [chrom2],
* and transformation cost matrix, this function returns 
* a pairwise sequence alignments which lies between two seeds [seed1] and [seed2] *)
let create_alied_subseq (seed1 : seed_t) (seed2 : seed_t) 
    (chrom1 : Sequence.s) (chrom2 : Sequence.s) (cost_mat : Cost_matrix.Two_D.m) = 
    let sta1 = seed1.sta1 + seed1.len in (* the first position is served fo gap *)
    let sta2 = seed1.sta2 + seed1.len in
        
    let len1 = seed2.sta1 - sta1 in 
    let len2 = seed2.sta2 - sta2 in                 
    let subseq1 = Sequence.sub chrom1 sta1 len1 in 
    let subseq2 = Sequence.sub chrom2 sta2 len2 in
    let use_ukk = false in (*module Seed and Block are with old annotation
    method, we don't use them anymore, it doens't matter if we use use_ukk here
    or not.*)
    let alied_subseq1, alied_subseq2, cost, ali_len = 
        Sequence.align2 subseq1 subseq2 cost_mat use_ukk in         
    alied_subseq1, alied_subseq2, cost
    
