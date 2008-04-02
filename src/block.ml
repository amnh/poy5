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

let () = SadmanOutput.register "Block" "$Revision: 2658 $"
(** The module contains default parameters and 
    funtions to create blocks between two chromosomes. *)


type pairChromPam_t = ChromPam.chromPairAliPam_t
type seed_t = Seed.seed_t
type direction_t = ChromPam.direction_t
type order_t = ChromPam.order_t
type subseq_t = Subseq.subseq_t

let deref = Utl.deref
let fprintf = Printf.fprintf
module IntSet = All_sets.Integers


(** Parameters are used to create blocks between two chromosomes*)
type blockPam_t = {
    (** Two consecutive blocks are connected if their distance and shift are
        smaller than thresholds *)
    max_connecting_dis : int;
    max_connecting_shift: int;
}   



type block_t = {
    mutable id : int; (** The block id *)
    mutable is_dum : bool; (** Dum blocks are used as boundary *)

    (** Seld-documented*)
    mutable sta1 : int;
    mutable sta2 : int;
    mutable en1 : int;
    mutable en2 : int;
    mutable direction : direction_t;
    mutable cost : int;
    mutable alied_seq1 : Sequence.s option;
    mutable alied_seq2 : Sequence.s option;

    (** A block contains a list of seeds *)
    mutable seed_ls : seed_t list;

    (** A chromosome is divided into consecutive sub-sequences *)
    mutable subseq1_id : int;
    mutable subseq2_id : int;
}
        

let blockPam_default = {
    max_connecting_dis = 30000; 
    max_connecting_shift = 3000;
}


let cloneBlockPam (donor : blockPam_t) = {    
    max_connecting_dis = donor.max_connecting_dis;
    max_connecting_shift = donor.max_connecting_shift;
}



let create_from_seed (block_id : int) (seed : seed_t) = {
    id = block_id;
    is_dum = false;
    sta1 = seed.Seed.sta1;
    sta2 = seed.Seed.sta2;
    en1 = seed.Seed.sta1 + seed.Seed.len - 1;
    en2 = seed.Seed.sta2 + seed.Seed.len - 1;
    direction = `Positive;
    cost = 0;
    seed_ls = [seed];
    alied_seq1 = None;
    alied_seq2 = None;
    
    subseq1_id = -1;
    subseq2_id = -1;
}       


let create_simple_block (block_id : int) (new_sta1 : int) 
        (new_en1 : int) (new_sta2 : int) (new_en2 : int) = 

{   id = block_id;
    is_dum = false;
    sta1 = new_sta1;
    sta2 = new_sta2;
    en1 = new_en1;
    en2 = new_en2;
    direction = `Positive;
    cost = 0;
    seed_ls = [];
    alied_seq1 = None;
    alied_seq2 = None;
    
    subseq1_id = -1;
    subseq2_id = -1;
}       

let get_dum_first_block ali_pam = {
    id = -1; 
    is_dum = true;
    sta1 = ali_pam.ChromPam.min_pos1 - 1; en1 = ali_pam.ChromPam.min_pos1 - 1; 
    sta2 = ali_pam.ChromPam.min_pos2 - 1; en2 = ali_pam.ChromPam.min_pos2 - 1;

    cost = ali_pam.ChromPam.sig_k * ali_pam.ChromPam.mat_cost;
    direction = `Positive;
    seed_ls = [Seed.get_dum_first_seed ali_pam.ChromPam.min_pos1 ali_pam.ChromPam.min_pos2];
    alied_seq1 = None;
    alied_seq2 = None;


    subseq1_id = -1;
    subseq2_id = -1;    
}

let get_dum_last_block ali_pam = {
    id = -1; 
    is_dum = true;
    sta1 = ali_pam.ChromPam.max_pos1 + 1; en1 = ali_pam.ChromPam.max_pos1 + 1; 
    sta2 = ali_pam.ChromPam.max_pos2 + 1; en2 = ali_pam.ChromPam.max_pos2 + 1;

    cost = ali_pam.ChromPam.sig_k * ali_pam.ChromPam.mat_cost;
    direction = `Positive;
    seed_ls = [Seed.get_dum_last_seed ali_pam.ChromPam.max_pos1 ali_pam.ChromPam.max_pos2];
    alied_seq1 = None;
    alied_seq2 = None;

    subseq1_id = -1;
    subseq2_id = -1;
}

let max_len (b : block_t) = 
    max (b.en1 - b.sta1 + 1) (b.en2 - b.sta2 + 1)


let invert (block : block_t) (min_pos2 : int) (max_pos2 : int) = 
    let new_sta2 = min_pos2 + (max_pos2 - block.en2) in 
    let new_en2 = min_pos2 + (max_pos2 - block.sta2) in 
    block.sta2 <- new_sta2;
    block.en2 <- new_en2;
    
    (match block.direction with
        | `Positive -> block.direction <- `Negative
        | _ -> block.direction <- `Positive);
    
    (match block.alied_seq2 with
        | None -> ()
        | Some seq -> Sequence.reverse_ip seq);
        
    List.iter (fun seed -> Seed.invert min_pos2 max_pos2 seed) block.seed_ls   
    

let get_pos (block : block_t) (order : order_t) = 
    match order with
        | `First -> block.sta1, block.en1
        | _ -> block.sta2, block.en2

let cmp_dia_dis (b1 : block_t) (b2 : block_t) =         
    abs ( (b1.en1 - b1.en2) - (b2.sta1 - b2.sta2) )



let print (block :  block_t) = 
    let dir = 
        match block.direction with
            | `Positive -> 1
            | _ -> -1
    in
    
    fprintf stdout "id: %i, (%i, %i) <--> (%i, %i) --> len: (%i, %i), "
        block.id block.sta1 block.en1  
        block.sta2 block.en2 
        (block.en1 - block.sta1 + 1) 
        (block.en2 - block.sta2 + 1);

    fprintf stdout "num_seed: %i, subseq: (%i, %i), cost: %i, direction: %i \n" 
        (List.length block.seed_ls) block.subseq1_id block.subseq2_id block.cost
        dir;
    List.iter Seed.print block.seed_ls;
    print_endline "-------------------------------------"
    

    

let add_seed (block : block_t) (seed : seed_t) = 
    block.sta1 <- seed.Seed.sta1;
    block.sta2 <- seed.Seed.sta2;
    block.seed_ls <- seed::block.seed_ls
                                

(** Compute the cost of a block which is a list of seeds *)
let cmp_cost_based_seed (block : block_t) (ali_pam : pairChromPam_t) =  
    let rec cmp (cur_seed : seed_t) (res_ls : seed_t list) (cost : int) = 
        match res_ls with 
            | [] -> cost
            | next_seed::tail ->
                let connecting_cost = 
                    Seed.cmp_connecting_cost cur_seed next_seed ali_pam in 
                cmp next_seed tail (cost + connecting_cost + 
                    next_seed.Seed.len *ali_pam.ChromPam.mat_cost) 
    in
    let seed_ls = block.seed_ls in 
    match seed_ls with 
        | [] -> 0 
        | first_seed::tail ->
            cmp first_seed tail 
                (first_seed.Seed.len * ali_pam.ChromPam.mat_cost)

    
(** Create a block from a list of seeds *)
let create_from_seed_ls (block_id : int) (seed_ls : seed_t list) 
        (ali_pam : pairChromPam_t) = 

    let head = List.hd seed_ls in  
    let tail = List.hd (List.rev seed_ls) in  

    let new_block = {      
        id = block_id;   
        is_dum = false;    
        sta1 = head.Seed.sta1;    
        sta2 = head.Seed.sta2;    
        en1 = tail.Seed.sta1 + tail.Seed.len - 1;    
        en2 = tail.Seed.sta2 + tail.Seed.len - 1;    
        direction = `Positive;    
        cost = 0;    
        seed_ls = seed_ls;    
        alied_seq1 = None;    
        alied_seq2 = None;    
        subseq1_id = -1;    
        subseq2_id = -1}     
    in      
    new_block.cost <- cmp_cost_based_seed new_block ali_pam;     
    new_block     

(*=======================================================*)


(** Compute the cost between two aligned sequences using cost parameters in
    ali_pam (ChromPam) *)
let cmp_ali_cost (seq1 : Sequence.s) (seq2 : Sequence.s) 
        (direction : direction_t) (ali_pam : pairChromPam_t) = 

    let len = Sequence.length seq1 in 
    let code1_arr = Sequence.to_array seq1 in          
    let code2_arr = Sequence.to_array seq2 in
    (if direction = `Negative then Utl.invert_subarr code2_arr 0 (Array.length code2_arr));

    let mat_cost = ali_pam.ChromPam.mat_cost in 
    let mismat_cost = ali_pam.ChromPam.mismat_cost in 
    let gap_opening_cost = ali_pam.ChromPam.gap_opening_cost in 
    let gap_ext_cost = ali_pam.ChromPam.gap_ext_cost in 
    
    let gap_code = Alphabet.get_gap Alphabet.nucleotides in 
    let rec count pos cost = 
        match pos >= len with
        | true -> cost
        | false ->
              match code1_arr.(pos) land code2_arr.(pos) > 0 with
              | true -> count (pos + 1) (cost + mat_cost)
              | false ->
                   match (code1_arr.(pos) = gap_code) || (code2_arr.(pos) = gap_code) with                       
                   | false ->  count (pos + 1) (cost + mismat_cost)
                   | true ->
                         if ( (code1_arr.(pos) = gap_code) 
                              && (code1_arr.(pos - 1) = gap_code) ) ||
                             ( (code2_arr.(pos) = gap_code) 
                               && (code2_arr.(pos - 1) = gap_code) ) then
                                 count (pos + 1) (cost + gap_ext_cost)

                         else count (pos + 1) (cost + gap_opening_cost)
    in

    if len = 0 then 0 
    else 
        match code1_arr.(0) = code2_arr.(0) with
        | true -> count 1 mat_cost
        | false ->
              if (code1_arr.(0) = gap_code) || (code2_arr.(0) = gap_code) then 
                  count 1 gap_opening_cost
              else count 1 mismat_cost
        
    
    


(** Seeds which are near each other are connected into
    local blocks *)
let find_local_block (seed_ls : seed_t list) (ali_pam : pairChromPam_t) = 

    let acc_cost_arr, back_arr, sorted_end_seed_arr = 
        Seed.create_local_connection seed_ls ali_pam in       
    
    let num_seed = List.length seed_ls in 
    let sorted_best_seed_arr = Array.copy sorted_end_seed_arr in 
    
    let cmp_seed (seed1 : seed_t) (seed2: seed_t) = 
        let id1 = seed1.Seed.id in
        let id2 = seed2.Seed.id in 
        compare acc_cost_arr.(id1) acc_cost_arr.(id2)
    in
        
    Array.sort cmp_seed sorted_best_seed_arr; 
    let avail_seed_arr = Array.make num_seed true in 
    
    let num_block = ref 0  in 
    let rev_block_ls = ref [] in
    
    for index = 0 to num_seed - 1 do
        let seed = sorted_best_seed_arr.(index) in 
        let seed_id = seed.Seed.id in 
        if avail_seed_arr.(seed_id) then 
        begin
            let seed_ls = ref [] in 

            let rec chase cur_seed_id = 
                if avail_seed_arr.(cur_seed_id) then 
                begin
                    let cur_seed = sorted_end_seed_arr.(cur_seed_id) in 
                    seed_ls := cur_seed::!seed_ls;
                    avail_seed_arr.(cur_seed_id) <- false;
                                
                    if back_arr.(cur_seed_id) != -1 then 
                        chase back_arr.(cur_seed_id)
                end 
            in
            chase seed_id;
            
            avail_seed_arr.(seed_id) <- false;

            let new_block = create_from_seed_ls !num_block !seed_ls ali_pam in          
            num_block := !num_block + 1; 
            rev_block_ls := new_block::!rev_block_ls                        

        end 
    done;   
    
    List.rev !rev_block_ls


let is_free (b : block_t) (sig1_arr : int array) (sig2_arr : int array) = 
    let rec travel sig_arr pos en =
        if pos > en then true
        else
            match sig_arr.(pos) = -1 with
                | true -> travel sig_arr (pos + 1) en
                | false -> false
    in
    if travel sig1_arr b.sta1 b.en1 = false then false
    else travel sig2_arr b.sta2 b.en2



let assign block sig1_arr sig2_arr = 
    let assign_subseq block_id sig_arr (sta : int) (en : int) = 
        for pos = sta to en do
            sig_arr.(pos) <- block_id
        done
    in 
    assign_subseq block.id sig1_arr block.sta1 block.en1;
    assign_subseq block.id sig2_arr block.sta2 block.en2


let create_sig_arr (block_ls : block_t list) ali_pam = 
    let _ = List.fold_left 
        (fun index block -> block.id <- index; index + 1) 0 block_ls in  

    let sig1_arr = Array.make (ali_pam.ChromPam.max_pos1 + 1) (-1) in
    let sig2_arr = Array.make (ali_pam.ChromPam.max_pos2 + 1) (-1) in

    List.iter (fun b -> assign b sig1_arr sig2_arr) block_ls;
    sig1_arr, sig2_arr




let select_separated_block b_ls ali_pam = 
    let b_arr = Array.of_list b_ls in 
    Array.sort (fun b1 b2 -> 
                    let len1 : float = float (max_len b1) in
                    let len2 : float = float (max_len b2) in
                    let cost1 = float b1.cost in 
                    let cost2 = float b2.cost in 
                    int_of_float ( len1 *. log(-.cost1) -. len2 *. log(-.cost2) )
               ) b_arr;
    
    let sig1_arr = Array.make (ali_pam.ChromPam.max_pos1 + 1) (-1) in
    let sig2_arr = Array.make (ali_pam.ChromPam.max_pos2 + 1) (-1) in

    List.fold_right 
        (fun b sep_bl_ls ->
             match is_free b sig1_arr sig2_arr with
             | true ->
                 assign b sig1_arr sig2_arr;
                 b::sep_bl_ls
             | false -> sep_bl_ls
        ) (Array.to_list b_arr) []

    
(** Given a list of seeds, align them to create the block*)        
let create_pos_alied_block (block : block_t) (seq1 : Sequence.s) 
        (seq2 : Sequence.s) (cost_mat : Cost_matrix.Two_D.m) ali_pam =

    let rev_alied_subseq1_ls = ref [] in 
    let rev_alied_subseq2_ls = ref [] in 
    let total_ali_cost = ref 0  in

    let rec create (pre_seed : seed_t) (cur_map_ls : seed_t list) = 
        let alied_pre_subseq1, alied_pre_subseq2, pre_cost = 
            Seed.create_alied_seed pre_seed seq1 seq2 cost_mat in 
        rev_alied_subseq1_ls := alied_pre_subseq1:: !rev_alied_subseq1_ls;
        rev_alied_subseq2_ls := alied_pre_subseq2:: !rev_alied_subseq2_ls;
        total_ali_cost := !total_ali_cost + pre_cost;
                
        match cur_map_ls with
            | [] -> ()
            | head::tail ->     
                let alied_subseq1, alied_subseq2, cost = 
                    Seed.create_alied_subseq pre_seed head seq1 seq2 cost_mat in 

                rev_alied_subseq1_ls := alied_subseq1::!rev_alied_subseq1_ls;
                rev_alied_subseq2_ls := alied_subseq2::!rev_alied_subseq2_ls;
                total_ali_cost := !total_ali_cost + cost;                               
                create head tail
    in
    create (List.hd block.seed_ls) (List.tl block.seed_ls);

    let alied_seq1 = UtlPoy.concat (List.rev !rev_alied_subseq1_ls) in  
    let alied_seq2 = UtlPoy.concat (List.rev !rev_alied_subseq2_ls) in
    
    block.alied_seq1 <- Some alied_seq1;
    block.alied_seq2 <- Some alied_seq2;
    block.cost <- cmp_ali_cost alied_seq1 alied_seq2 block.direction ali_pam
    




    
let is_inside (seed : seed_t) (block : block_t) = 
    (seed.Seed.sta1 >= block.sta1) && 
    (seed.Seed.sta1 + seed.Seed.len - 1 <= block.en1) && 
    (seed.Seed.sta2 >= block.sta2) && 
    (seed.Seed.sta2 + seed.Seed.len - 1 <= block.en2)

        

(** create an array of separated sub-sequences *)
let determine_separated_subseq (block_ls : block_t list) (order : order_t)
        (max_pos : int) subseq_type : subseq_t list = 
    (* create the label_arr.(pos) -> list of blocks containing this position*)
    let label_arr = Array.create (max_pos + 2) [] in        
    let rec create_label (block : block_t) = 
        if block.is_dum = false then 
        begin
            let sta, en = get_pos block order in 
            let block_id = block.id in
            (*fprintf stdout "%i %i %i" block_id sta en; print_newline ();*)
            for pos = sta to en do
                label_arr.(pos) <- block_id::label_arr.(pos) 
            done    
        end
    in              
    List.iter create_label block_ls;
    
    
    let rev_sep_subseq_ls = ref [] in
    let sta = ref 0 in
    let num_sep_subseq = ref 0 in 
    for pos = 1 to max_pos + 1 do               
        let create_new_subseq (sta : int) (en : int) = 
            let new_subseq = {Subseq.id = !num_sep_subseq + 1; 
                              Subseq.sta = sta; 
                              Subseq.en = en; 
                              Subseq.block_id_ls = label_arr.(sta)} in               

            rev_sep_subseq_ls := new_subseq::!rev_sep_subseq_ls;               
            num_sep_subseq := !num_sep_subseq + 1;              
        in
        
        if (pos = max_pos + 1) || 
        ( (Utl.compare_non_dec_list label_arr.(pos) label_arr.(pos - 1)) = false) then 
        begin
            if (  subseq_type = `Both) || 
               ( (subseq_type = `Alied) &&  (label_arr.(!sta) != []) ) || 
               ( (subseq_type = `Deleted) &&  (label_arr.(!sta) = []) )
            then create_new_subseq !sta (pos - 1);

            sta := pos
        end
    done;
    
    List.rev !rev_sep_subseq_ls
(************ End of create separeated subseq list *****************)
(*******************************************************************)            



        

let create_alied_block_ls (block_ls : block_t list) (ali_pam : pairChromPam_t) 
        (seq1 : Sequence.s) (seq2 : Sequence.s) cost_mat =
    
    List.iter 
        (fun block -> 
             if block.direction = `Positive then    
                 create_pos_alied_block block seq1 seq2 cost_mat ali_pam
        ) block_ls; 
    
    
    let min_pos2 = ali_pam.ChromPam.min_pos2 in 
    let max_pos2 = ali_pam.ChromPam.max_pos2 in 
    let com_seq2 = Sequence.complement_chrom Alphabet.nucleotides seq2 in
    List.iter (fun block -> 
                   if block.direction = `Negative then begin
                       invert block min_pos2 max_pos2;
                       create_pos_alied_block block seq1 com_seq2 cost_mat ali_pam;
                       invert block min_pos2 max_pos2;
                    end) block_ls


(*==========================================================*)
 
let check_sep_block (sep_block_ls : block_t list) = 
    let sep_block_arr = Array.of_list sep_block_ls in 
    let num_block = Array.length sep_block_arr in 
    for no1 = 0 to num_block - 2 do
        for no2 = no1 + 1 to num_block - 1 do
            if ((sep_block_arr.(no1).en1 < 
                     sep_block_arr.(no2).sta1) ||
               (sep_block_arr.(no2).en1 < 
                    sep_block_arr.(no1).sta1) ) &&
               ( (sep_block_arr.(no1).en2 < 
                      sep_block_arr.(no2).sta2) ||
               (sep_block_arr.(no2).en2 < 
                    sep_block_arr.(no1).sta2) )
            then ()
            else begin
                print sep_block_arr.(no1);
                print sep_block_arr.(no2);
                failwith "Fucking diving block, they are still overlapped"
            end
        done
    done;
    print_endline "All block are separated!!!"


let prepen (b1 : block_t) (b2 : block_t) (block_pam : blockPam_t) 
        (seq1 : Sequence.s) (seq2 : Sequence.s) (cost_mat : Cost_matrix.Two_D.m) 
        (ali_pam : pairChromPam_t) =

    
    let alied_between_seq1, alied_between_seq2, cost = 
        UtlPoy.create_subalign2 seq1 seq2 cost_mat (b1.en1 + 1)
        (b2.sta1 - 1) (b1.en2 + 1) (b2.sta2 - 1) 
    in 
    

    let between_cost = cmp_ali_cost alied_between_seq1 alied_between_seq2
        `Positive ali_pam 
    in 

    
    
    b2.sta1 <- b1.sta1;
    b2.sta2 <- b1.sta2;
    b2.cost <- b1.cost + between_cost + b2.cost;
    b2.seed_ls <- b1.seed_ls @ b2.seed_ls;
    b2.alied_seq1 <- Some (UtlPoy.concat [ (deref b1.alied_seq1);
                                           alied_between_seq1; 
                                           (deref b2.alied_seq1)] );

    b2.alied_seq2 <- Some (UtlPoy.concat [ (deref b1.alied_seq2); 
                                           alied_between_seq2; 
                                           (deref b2.alied_seq2)] )

    
let connect_pos_consecutive_block (block_ls : block_t list) 
        (block_pam : blockPam_t) (seq1 : Sequence.s) (seq2 : Sequence.s) 
        (cost_mat : Cost_matrix.Two_D.m) (ali_pam : pairChromPam_t)  = 

    let sorted_block_arr = Array.of_list block_ls in
    Array.sort (fun b1 b2 -> compare b1.en1 b2.en1) sorted_block_arr;
    let _ = Array.fold_left 
        (fun id b -> b.id <- id; id + 1) 0 sorted_block_arr in 

    
    let marker1_set = ref IntSet.empty in
    let marker2_set = ref IntSet.empty in
    List.iter (fun b -> marker1_set := IntSet.add b.sta1 !marker1_set;
                        marker1_set := IntSet.add b.en1 !marker1_set;
                        marker2_set := IntSet.add b.sta2 !marker2_set;
                        marker2_set := IntSet.add b.en2 !marker2_set
              ) block_ls;

    let marker1_arr = Array.of_list (IntSet.elements !marker1_set) in 
    let marker2_arr = Array.of_list (IntSet.elements !marker2_set) in 

    let max_pos1 = ali_pam.ChromPam.max_pos1 in 
    let max_pos2 = ali_pam.ChromPam.max_pos2 in 
    let sig1_arr = Array.make (max_pos1 + 1) [] in 
    let sig2_arr = Array.make (max_pos2 + 1) [] in 
    List.iter (fun b -> if b.direction = `Positive then 
                begin
                    sig1_arr.(b.sta1) <- b.id::sig1_arr.(b.sta1);
                    sig2_arr.(b.sta2) <- b.id::sig2_arr.(b.sta2);
                end) block_ls;



    let num_marker1 = Array.length marker1_arr in 
    let num_marker2 = Array.length marker2_arr in 

    let do_connection (b : block_t) =  
        let m1 = Utl.binary_index_search marker1_arr b.en1 in 
        let m2 = Utl.binary_index_search marker2_arr b.en2 in 
        if (m1 >= 0) && (m1 + 1 < num_marker1) && 
            (m2 >= 0) && (m2 + 1 < num_marker2) then 
            begin
                let e1 = marker1_arr.(m1 + 1) in 
                let e2 = marker2_arr.(m2 + 1) in 
                let dis1 = e1 - b.en1 in
                let dis2 = e2 - b.en2 in 

 
                if (max dis1 dis2 <= block_pam.max_connecting_dis) &&
                    (abs (dis1 - dis2) <= block_pam.max_connecting_shift) then
                begin
                    let com_block_id_ls = List.filter 
                        (fun id -> List.mem id sig2_arr.(e2) ) sig1_arr.(e1) in                    
                    if List.length com_block_id_ls > 0 then
                    begin
                        List.iter (fun id -> prepen b sorted_block_arr.(id) 
                                       block_pam seq1 seq2 cost_mat ali_pam
                                  ) com_block_id_ls;
                        b.is_dum <- true
                    end
                end
            end
    in

    Array.iter (fun b -> if (b.is_dum = false) && (b.direction = `Positive) then 
                    do_connection b) sorted_block_arr;
    let pos_con_block_ls = Array.fold_left 
        (fun l b -> if b.is_dum then l else b::l) [] sorted_block_arr in 
    pos_con_block_ls




let connect_consecutive_block (block_ls : block_t list) (block_pam : blockPam_t)
        (seq1 : Sequence.s) (seq2 : Sequence.s) (cost_mat : Cost_matrix.Two_D.m)
        (ali_pam : pairChromPam_t) = 

    let pos_con_block_ls = connect_pos_consecutive_block block_ls block_pam seq1
        seq2 cost_mat ali_pam in 
   
    pos_con_block_ls



let create_subseq_id subseq_type (sep_block_ls : block_t list) 
        (ali_pam : pairChromPam_t) =  
                                                                            
    let _ = List.fold_left 
        (fun index block -> block.id <- index; index + 1) 0 sep_block_ls in 
    let sep_block_arr = Array.of_list sep_block_ls in 
    

    let subseq1_ls = determine_separated_subseq sep_block_ls 
        `First ali_pam.ChromPam.max_pos1 subseq_type in 
    let subseq2_ls = determine_separated_subseq sep_block_ls 
        `Second ali_pam.ChromPam.max_pos2 subseq_type in 


    List.iter (fun subseq -> List.iter 
       (fun block_id -> sep_block_arr.(block_id).subseq1_id <- subseq.Subseq.id
       ) subseq.Subseq.block_id_ls) subseq1_ls;

    List.iter (fun subseq -> List.iter 
       (fun block_id -> sep_block_arr.(block_id).subseq2_id <- subseq.Subseq.id
       ) subseq.Subseq.block_id_ls) subseq2_ls;
    
    sep_block_ls, subseq1_ls, subseq2_ls



let create_median ?(approx=`BothSeq) (block : block_t) cost_mat = 
    let alied_seq1 = Utl.deref block.alied_seq1 in 
    let alied_seq2 = Utl.deref block.alied_seq2 in 


    match block.direction = `Positive with
    | true -> UtlPoy.create_median_seq ~approx:approx alied_seq1 alied_seq2 cost_mat 
    | false ->  
          Sequence.reverse_ip alied_seq2;  
          let med, cost = UtlPoy.create_median_seq ~approx:approx alied_seq1 alied_seq2 cost_mat in  
          Sequence.reverse_ip alied_seq2; 
          med, cost 


let find_block  block_ls subseq1_id subseq2_id = 
    let rec check cur_block_ls =
        match cur_block_ls with
        | [] -> None
        | hd::tl ->
              if (hd.subseq1_id = subseq1_id) 
                  && (hd.subseq2_id = subseq2_id) then Some hd
              else check tl
    in
    check block_ls


let find_subseq1  block_ls subseq1_id = 
    let rec check cur_block_ls =
        match cur_block_ls with
        | [] -> None
        | hd::tl ->
              if (hd.subseq1_id = subseq1_id) then Some hd
              else check tl
    in
    check block_ls


