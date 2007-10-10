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

let () = SadmanOutput.register "AnnchromAli" "$Revision: 911 $"

(** The implementation of funtions to calculate the cost, alignments and medians
    between annotated chromosomes where both point mutations and rearrangement operations
    are considered *)
let fprintf = Printf.fprintf

type seq_t = {
    seq : Sequence.s;
    seq_ref_code : int;
    alied_med : Sequence.s;
    seq_ord1 : int;
    alied_seq1 : Sequence.s;

    seq_ord2 : int;
    alied_seq2 : Sequence.s;

}

type annchrom_t = {
    seq_arr : seq_t array; 
    ref_code : int; 
    ref_code1 : int;    (** Child's code *)    
    ref_code2 : int;  (** Child's code *)
    cost1 : int;
    recost1 : int;
    cost2 : int;
    recost2 : int;
}


let clone_seq s = {
    seq = Sequence.clone s.seq;
    seq_ref_code = s.seq_ref_code;
    alied_med = Sequence.clone s.alied_med;
    seq_ord1 = s.seq_ord1;
    alied_seq1 = Sequence.clone s.alied_seq1;
    seq_ord2 = s.seq_ord2;
    alied_seq2 = Sequence.clone s.alied_seq2;
}


let clone_med m = {
    seq_arr = Array.map clone_seq m.seq_arr;
    ref_code = m.ref_code;
    ref_code1 = m.ref_code1;
    ref_code2 = m.ref_code2;
    cost1 = m.cost1;
    cost2 = m.cost2;
    recost1 = m.recost1;
    recost2 = m.recost2;
}


(** Parameters used to align two general character sequences *)
type annchromPam_t = {
    re_meth : Data.re_meth_t;
    keep_median : int;
    circular : int;
    swap_med : int;
    approx : ChromPam.order_t;
    symmetric : bool;
    locus_indel_cost : (int * int);
}

let annchromPam_default = {
    re_meth = `Breakpoint 10;
    keep_median = 1;
    circular = 0;
    swap_med = 1;
    approx = `BothSeq;
    symmetric = false;
    locus_indel_cost = (10, 100);
}

let init_seq_t (seq, code) = {
    seq = seq;
    seq_ref_code = code;
    alied_med = seq;
    seq_ord1 = -1;
    alied_seq1 = seq;
    seq_ord2 = -1;
    alied_seq2 = seq;
    
}

let init seq_arr = {
    seq_arr = Array.map init_seq_t seq_arr;
    ref_code = Utl.get_new_chrom_ref_code ();
    ref_code1 = -1;
    ref_code2 = -1;
    cost1 = 0;
    cost2 = 0;
    recost1 = 0;
    recost2 = 0;
}

let printMap seq_arr =
    print_endline "Start printing alied_annchrom ";
    Array.iter (fun m -> 
                    fprintf stdout "Order1: %i, Order2: %i\n"  m.seq_ord1 m.seq_ord2;
                    UtlPoy.printDNA m.alied_seq1;
                    UtlPoy.printDNA m.alied_seq2;
                    UtlPoy.printDNA m.alied_med;
               ) seq_arr;
    print_endline "End printing alied_annchrom ";
    print_newline ()


let swap_seq s = 
    {
        s with seq_ord1 = s.seq_ord2;
            alied_seq1 = s.alied_seq2;
            seq_ord2 = s.seq_ord1;
            alied_seq2 = s.alied_seq1                
    } 



let swap_med m = 
     {m with ref_code1 = m.ref_code2;
         ref_code2 = m.ref_code1;
         cost1 = m.cost2;
         recost1 = m.recost2;
         cost2 = m.cost1;
         recost2 = m.recost1;
         seq_arr = Array.map swap_seq m.seq_arr
     }



let get_seq_arr t =
    Array.map (fun seq -> seq.seq) t.seq_arr

(** for implied alignments *)
let convert_map med = 
    let gap = Alphabet.gap in 
    let num_frag = Array.length med.seq_arr in

    let alied_med_arr = Array.map (fun seg -> seg.alied_med) med.seq_arr in
    let alied_seq1_arr = Array.init num_frag 
        (fun idx1 ->
             let seg = List.find (fun seg -> seg.seq_ord1 = idx1)
                 (Array.to_list med.seq_arr) 
             in
             seg.alied_seq1)
    in 

    let alied_seq2_arr = Array.init num_frag 
        (fun idx2 ->
             let seg = List.find (fun seg -> seg.seq_ord2 = idx2)
                 (Array.to_list med.seq_arr) 
             in
             seg.alied_seq2)
    in 

    let create_pos seq_arr = 
        let pos = ref (-1) in
        Array.map
            (fun s ->         
                 Array.init (Sequence.length s) 
                     (fun idx -> 
                          let code = Sequence.get s idx in
                          match code = gap with
                          | true -> -1
                          | false ->                                
                              incr pos; !pos)

            ) seq_arr 
    in 
    let pos_mat = create_pos alied_med_arr in
    let pos1_mat = create_pos alied_seq1_arr in 
    let pos2_mat = create_pos alied_seq2_arr in 


    let rev_map = ref [] in 
    for seg_id = 0 to num_frag - 1 do
        let seg = med.seq_arr.(seg_id) in 
        let len = Array.length pos_mat.(seg_id) in
        for idx = 0 to len - 1 do
            let p = pos_mat.(seg_id).(idx) in 
            let code = Sequence.get seg.alied_med idx in 
            let p1 = pos1_mat.(seg.seq_ord1).(idx) in 
            let code1 = Sequence.get seg.alied_seq1 idx in 
            let p2 = pos2_mat.(seg.seq_ord2).(idx) in 
            let code2 = Sequence.get seg.alied_seq2 idx in 
            rev_map := (p, code, p1, code1, p2, code2)::!rev_map
        done
    done; 

    List.rev !rev_map  


let get_annchrom_pam user_annchrom_pam = 
    let chrom_pam = annchromPam_default in  
    let chrom_pam = 
        match user_annchrom_pam.Data.re_meth with
        | None -> chrom_pam
        | Some re_meth -> {chrom_pam with re_meth = re_meth}
    in 

    let chrom_pam = 
        match user_annchrom_pam.Data.circular with
        | None -> chrom_pam
        | Some circular -> {chrom_pam with circular = circular}
    in 

    let chrom_pam = 
        match user_annchrom_pam.Data.keep_median with
        | None -> chrom_pam
        | Some keep_median -> {chrom_pam with keep_median = keep_median}
    in 
    
    let chrom_pam = 
        match user_annchrom_pam.Data.swap_med with
        | None -> chrom_pam
        | Some swap_med -> {chrom_pam with swap_med = swap_med}
    in 

    let chrom_pam = 
        match user_annchrom_pam.Data.approx with
        | None -> chrom_pam
        | Some approx -> 
              if approx then {chrom_pam with approx = `First}
              else {chrom_pam with approx = `BothSeq}
    in 


    let chrom_pam = 
        match user_annchrom_pam.Data.symmetric with
        | None -> chrom_pam
        | Some sym -> {chrom_pam with symmetric = sym}
    in 

    let chrom_pam = 
        match user_annchrom_pam.Data.locus_indel_cost with
        | None -> chrom_pam
        | Some locus_indel_cost -> {chrom_pam with locus_indel_cost = locus_indel_cost}
    in 
    chrom_pam

let print annchrom alpha = 
    Array.iter (fun s -> 
                    let seq = Sequence.to_string s.seq alpha in 
                    fprintf stdout "(%i, %s) | " s.seq_ref_code seq) annchrom.seq_arr;
    print_newline ()

let split chrom =  
    let seq_arr = Array.map (fun s -> s.seq) chrom.seq_arr in  
    let code_arr = Array.map (fun s -> s.seq_ref_code) chrom.seq_arr in  
    seq_arr, code_arr 




(** Given two arrays of sequences [seq1_arr] and [seq2_arr], 
 *  create the general cost matrix and corresponding code arrays  *)
let create_pure_gen_cost_mat seq1_arr seq2_arr cost_mat ali_pam =        
    let seq1_arr = Array.mapi (fun ith seq -> seq, (ith * 2 + 1) ) seq1_arr in 

    let len1 = Array.length seq1_arr in 

    let seq2_arr = Array.mapi 
        (fun ith seq -> seq, (ith * 2 + 1 + len1 * 2) ) seq2_arr 
    in 

    let len2 = Array.length seq2_arr in 
    let gen_gap_code = (len1 + len2) * 2 + 1 in 


    let pure_gen_cost_mat = Array.make_matrix (gen_gap_code + 1) (gen_gap_code + 1)
        Utl.infinity 
    in 

    pure_gen_cost_mat.(gen_gap_code).(gen_gap_code) <- 0;
    
    let update (seq1, code1) (seq2, code2) =
        let com_seq1 = Sequence.complement_chrom Alphabet.nucleotides seq1 in 
        let com_seq2 = Sequence.complement_chrom Alphabet.nucleotides seq2 in 
(*
        UtlPoy.printDNA seq1;
        UtlPoy.printDNA com_seq1;
        print_newline ();
        UtlPoy.printDNA seq2;
        UtlPoy.printDNA com_seq2;
        print_endline "----------------------------------";
*)        
        let _, _, cost, _ = UtlPoy.align2 seq1 seq2 cost_mat in 
  (*      fprintf stdout "%i %i -> %i\n" code1 code2 cost;*)
        pure_gen_cost_mat.(code1).(code2) <- cost;
        pure_gen_cost_mat.(code2).(code1) <- cost;

        let _, _, cost, _ = UtlPoy.align2 seq1 com_seq2 cost_mat in 
(*        fprintf stdout "%i %i -> %i\n" code1 (-code2) cost; *)
        pure_gen_cost_mat.(code1).(code2 + 1) <- cost;
        pure_gen_cost_mat.(code2 + 1).(code1) <- cost;

        let _, _, cost, _ = UtlPoy.align2 com_seq1 seq2 cost_mat in 
(*        fprintf stdout "%i %i -> %i\n" (-code1) code2 cost; *)
        pure_gen_cost_mat.(code1 + 1).(code2) <- cost;
        pure_gen_cost_mat.(code2).(code1 + 1) <- cost;

        let _, _, cost, _ = UtlPoy.align2 com_seq1 com_seq2 cost_mat in 
(*        fprintf stdout "%i %i -> %i\n" (-code1) (-code2) cost;*)
        pure_gen_cost_mat.(code1 + 1).(code2 + 1) <- cost;
        pure_gen_cost_mat.(code2 + 1).(code1 + 1) <- cost;

    in 

        
    Array.iter (fun (seq1, code1) ->
                    Array.iter (fun (seq2, code2) -> 
                                    update (seq1, code1) (seq2, code2)
                               ) seq2_arr) seq1_arr;


    let update_gap (seq, code) = 
        pure_gen_cost_mat.(gen_gap_code).(code) <- UtlPoy.cmp_gap_cost ali_pam.locus_indel_cost seq;
        pure_gen_cost_mat.(code).(gen_gap_code) <- pure_gen_cost_mat.(gen_gap_code).(code);
    in 
    Array.iter update_gap seq1_arr;
    Array.iter update_gap seq2_arr;


    let code1_arr = Array.map (fun (seq, code) -> code) seq1_arr in 
    let code2_arr = Array.map (fun (seq, code) -> code) seq2_arr in 

    pure_gen_cost_mat, code1_arr, code2_arr, gen_gap_code



(** Given two annotated chromosomes [chrom1] and [chrom2], compute 
 * the total cost between them which is comprised of editing cost and 
 * rearrangement cost *)
let cmp_simple_cost (chrom1: annchrom_t) (chrom2 : annchrom_t) 
        cost_mat alpha ali_pam =


    let chrom_len1 = Array.fold_left (fun len s -> len + Sequence.length s.seq) 0 chrom1.seq_arr in     
    let chrom_len2 = Array.fold_left (fun len s -> len + Sequence.length s.seq) 0 chrom2.seq_arr in 
    
    if (chrom_len1 < 2) || (chrom_len2 < 2) then 0, 0
    else begin
        let seq1_arr, _ = split chrom1 in  
        let seq2_arr, _ = split chrom2 in  
    
        let pure_gen_cost_mat, code1_arr, code2_arr, gen_gap_code = 
            create_pure_gen_cost_mat seq1_arr seq2_arr cost_mat ali_pam  
        in 
    
        let total_cost, (recost1, recost2), _, _ = 
            GenAli.create_gen_ali_code  `Annchrom code1_arr code2_arr 
                pure_gen_cost_mat gen_gap_code  
                ali_pam.re_meth ali_pam.swap_med 
                ali_pam.circular  
        in 
        total_cost, (recost1 + recost2)
    end 



let cmp_cost (chrom1: annchrom_t) (chrom2 : annchrom_t) 
        cost_mat alpha annchrom_pam = 
    let ali_pam = get_annchrom_pam annchrom_pam in     
    match ali_pam.symmetric with
    | true ->
          let cost12 = cmp_simple_cost chrom1 chrom2 cost_mat alpha ali_pam in
          let cost21 = cmp_simple_cost chrom2 chrom1 cost_mat alpha ali_pam in
          min cost12 cost21
    | false ->
          if compare chrom1 chrom2 < 0 then
              cmp_simple_cost chrom1 chrom2 cost_mat alpha ali_pam
          else cmp_simple_cost chrom2 chrom1 cost_mat alpha ali_pam

              
(** Given two annotated chromosomes [chrom1] and [chrom2], 
 * find all median chromoromes between [chrom1] and [chrom2]. 
 * Rearrangements are allowed *) 
let find_simple_med2_ls (chrom1: annchrom_t) (chrom2 : annchrom_t) 
        (cost_mat : Cost_matrix.Two_D.m) alpha ali_pam = 
    
    let chrom_len1 = Array.fold_left (fun len s -> len + Sequence.length s.seq) 0 chrom1.seq_arr in     
    let chrom_len2 = Array.fold_left (fun len s -> len + Sequence.length s.seq) 0 chrom2.seq_arr in 
    if (chrom_len1 < 2) then 0,0, [chrom2]
    else if chrom_len2 < 2 then 0,0, [chrom1]
    else begin    
        let approx = ali_pam.approx in 

        let seq1_arr, _ = split chrom1 in  
        let seq2_arr, _ = split chrom2 in    
        
        let pure_gen_cost_mat, code1_arr, code2_arr, gen_gap_code = 
            create_pure_gen_cost_mat seq1_arr seq2_arr cost_mat ali_pam  
        in 
    

    
        let total_cost, (recost1, recost2), alied_code1_arr, alied_code2_arr = 
            GenAli.create_gen_ali_code  `Annchrom code1_arr code2_arr 
                pure_gen_cost_mat gen_gap_code  
                ali_pam.re_meth ali_pam.swap_med 
                ali_pam.circular  
        in 


        let ali_len = Array.length alied_code1_arr in 
    
        let ali_chrom = Array.init ali_len
            (fun idx ->   
                 let idx1 = Utl.find_index code1_arr alied_code1_arr.(idx) compare in   
                 let seq1 =  
                     match idx1 with  
                     | -1 -> None
                     | _ -> Some chrom1.seq_arr.(idx1).seq
                 in                             
    
                 let pos_code2 = 
                     if alied_code2_arr.(idx) mod 2 = 0 then alied_code2_arr.(idx) - 1
                     else alied_code2_arr.(idx)
                 in 
                 let idx2 = Utl.find_index code2_arr pos_code2 compare in  
                 let seq2 =
                     match idx2 with 
                     | -1 -> None
                     | _ -> 
                           if pos_code2 = alied_code2_arr.(idx) then Some chrom2.seq_arr.(idx2).seq
                           else Some (Sequence.complement_chrom Alphabet.nucleotides chrom2.seq_arr.(idx2).seq)
                 in                            
                     
                 let alied_med_seq, alied_seq1, alied_seq2 =
                     match seq1, seq2 with   
                     | Some seq1, Some seq2 -> 
                           let med, alied_seq1, alied_seq2, _  =
                               UtlPoy.create_median  ~approx:approx seq1 seq2 cost_mat  
                           in
                           med, alied_seq1, alied_seq2
                     | Some seq1, None -> 
                           let len1 = Sequence.length seq1 in 
                           let seq2 = UtlPoy.create_gap_seq len1 in 
                           let med, _ = UtlPoy.create_median_seq ~approx:approx seq1 seq2 cost_mat in 
                           med, seq1, seq2
                     | None, Some seq2 -> 
                           let len2 = Sequence.length seq2 in 
                           let seq1 = UtlPoy.create_gap_seq len2 in 
                           let med, _ = UtlPoy.create_median_seq ~approx:approx seq1 seq2 cost_mat in
                           med, seq1, seq2
                     | _, _ -> UtlPoy.get_empty_seq (), UtlPoy.get_empty_seq (), UtlPoy.get_empty_seq ()
                 in 
                 (alied_med_seq, idx1, alied_seq1, idx2, alied_seq2)
            )
        in 
    
    
    
        let create_dynamic_med (re_seq2, recost1, recost2) = 
            let rec follow index_ls index =
                if (index = ali_len) ||  
                    (alied_code2_arr.(index) != gen_gap_code) then index_ls  
                else follow (index::index_ls) (index + 1)
            in  
    
            let get_index re_seq2 = 
                let start_index_ls = follow [] 0 in
                let rev_index_ls = 
                    Array.fold_left  
                    (fun index_ls code2 ->
                         let index = Utl.find_index alied_code2_arr (abs code2) compare in
                         
                         match code2 > 0 with   
                         | true -> follow (index::index_ls) (index + 1)
                         | false -> follow (-index::index_ls) (index + 1)                     
                    ) start_index_ls re_seq2
                in 
                List.rev rev_index_ls
            in 
                                     
            let index_ls = get_index re_seq2 in         
    
            let med = 
                List.fold_right  
                    (fun index med ->                     
                         let alied_med, seq_ord1, alied_seq1, seq_ord2, alied_seq2 =  ali_chrom.(index) in
                         {seq= UtlPoy.delete_gap alied_med;  
                          seq_ref_code = Utl.get_new_seq_ref_code();  
                          alied_med = alied_med;
                          seq_ord1 = seq_ord1; 
                          alied_seq1 = alied_seq1;
                          seq_ord2 = seq_ord2;
                          alied_seq2 = alied_seq2;
                         }::med 
                    ) index_ls []
            in   


            {seq_arr = Array.of_list med; 
             ref_code = Utl.get_new_chrom_ref_code ();
             ref_code1 = chrom1.ref_code;
             ref_code2 = chrom2.ref_code;
             cost1 = total_cost - recost2;
             cost2 = total_cost - recost1;
             recost1 = recost1;
             recost2 = recost2;
            }
        in  
    
    
        let re_code2_arr = Utl.filterArr alied_code2_arr  
            (fun code2 -> code2 != gen_gap_code)  
        in   
        let code2_arr = UtlGrappa.get_ordered_permutation re_code2_arr in 
             
        let all_order_ls =   
            if (Utl.equalArr code2_arr re_code2_arr compare) ||  
                (ali_pam.keep_median = 1) ||
                (ali_pam.approx = `First) then [re_code2_arr, recost1, recost2]   
            else [(re_code2_arr, recost1, recost2); (code2_arr, recost2, recost1)]   
        in   
    

        let med_ls =  List.fold_left   
            (fun med_ls (re_seq2, recost1, recost2)  ->   
                 let med = create_dynamic_med (re_seq2, recost1, recost2) in 
                 med::med_ls  
            ) [] all_order_ls     
        in 
        total_cost, (recost1 + recost2), med_ls           
    end




(** Given two annotated chromosomes [chrom1] and [chrom2], 
 * find all median chromoromes between [chrom1] and [chrom2]. 
 * Rearrangements are allowed *) 
let find_med2_ls (chrom1: annchrom_t) (chrom2 : annchrom_t) 
        (cost_mat : Cost_matrix.Two_D.m) alpha annchrom_pam = 
    
    let ali_pam = get_annchrom_pam annchrom_pam in          
    match ali_pam.symmetric with
    | true ->
          let cost12, recost12, med12_ls = find_simple_med2_ls chrom1 chrom2
              cost_mat alpha ali_pam in 
          
          let ali_pam = 
              if ali_pam.approx = `First then {ali_pam with approx = `Second}
              else ali_pam
          in 
          let cost21, recost21, med21_ls = find_simple_med2_ls chrom2 chrom1
              cost_mat alpha ali_pam in 
          if cost12 <= cost21 then cost12, recost12, med12_ls
          else begin 
              let med12_ls = List.map swap_med med21_ls in 
              cost21, recost21, med12_ls
          end 

    | false ->
          let med1, med2, ali_pam, swaped = 
              match compare chrom1 chrom2 < 0 with 
              | true ->  chrom1, chrom2, ali_pam, false
              | false -> 
                    let ali_pam = 
                        if ali_pam.approx = `First then {ali_pam with approx = `Second}
                        else ali_pam
                    in 
                    chrom2, chrom1, ali_pam, true
          in 

          let cost, recost, med_ls = find_simple_med2_ls med1 med2 cost_mat
              alpha ali_pam in
          let med_ls = 
              match swaped with
              | false -> med_ls
              | true -> List.map swap_med med_ls 
          in 
          cost, recost, med_ls


 



(** Given two annotated chromosomes [chrom1] and [chrom2], 
 * compare chrom1 and chrom2 *)
let compare annchrom1 annchrom2 = 
    
    let seq1_arr, _ = split annchrom1 in 
    let seq2_arr, _ = split annchrom2 in 
    let len1 = Array.length seq1_arr in
    let len2 = Array.length seq2_arr in 
    match len1 = len2 with
    | false -> len1 - len2
    | true ->
          let rec compare_seq index = 
              if index = len1 then 0 
              else begin 
                  let cmp = Sequence.compare seq1_arr.(index) seq2_arr.(index) in
                  match cmp with 
                  | 0 -> compare_seq (index + 1)
                  | _ -> cmp
              end 
          in 
          compare_seq 0


let find_approx_med2 med1 med2 med12 = 
    let new_med12 = clone_med med12 in 
    let new_seq_arr = Array.map (fun seq -> 
                                     {seq with seq_ref_code = Utl.get_new_seq_ref_code ()}
                                ) new_med12.seq_arr
    in 
    let ref_code1 = med1.ref_code in 
    let ref_code2 = med2.ref_code in 
    let ref_code = Utl.get_new_chrom_ref_code () in 
    {new_med12 with 
         seq_arr = new_seq_arr;
         ref_code = ref_code; ref_code1 = ref_code1; ref_code2 = ref_code2}



(** Assign ids for all sequenes of this chromosome *)
let assign_seq_ref annchrom seq_ref_code = 
    let seq_ls, seq_id = Array.fold_right 
        (fun seq (annchrom, seq_id) ->
             match seq.seq_ref_code with 
             | -1 ->
                   let seq = {seq with seq_ref_code = seq_ref_code} in 
                   (seq::annchrom), (seq_ref_code + 1) 
             | _ -> (seq::annchrom), seq_ref_code
        ) annchrom.seq_arr ([], seq_ref_code)
    in 
    let seq_arr = Array.of_list seq_ls in 
    {annchrom with seq_arr = seq_arr}, seq_ref_code



let create_map med child_ref = 
    let str = string_of_int in  


    let seq_arr = Array.mapi 
        (fun med_id m ->
             let p_ref_code, p_seq_ord = (str med.ref_code), (str med_id) in

             let c_ref_code, c_seq_ord = 
                 match child_ref = med.ref_code1 with
                 | true -> (str med.ref_code1), (str m.seq_ord1)
                 | false -> (str med.ref_code2), (str m.seq_ord2)
             in 
             let attributes = [(Tags.GenomeMap.a_ref_code, p_ref_code);
                               (Tags.GenomeMap.a_seq_order, p_seq_ord);
                               (Tags.GenomeMap.a_dir_seg, "+");
                               (Tags.GenomeMap.d_ref_code, c_ref_code); 
                               (Tags.GenomeMap.d_seq_order, c_seq_ord);
                               (Tags.GenomeMap.d_dir_seg, "+")
                              ] 
             in 
             let m : Tags.output = (Tags.GenomeMap.seg, attributes, `String "") in 
             `Single m
        ) med.seq_arr
    in 


    let chrom_map : Tags.output = 
        (Tags.GenomeMap.chrom, [], `Structured (`Set  (Array.to_list seq_arr))) 
    in 
    match child_ref = med.ref_code1 with
    | true -> med.cost1, med.recost1, chrom_map
    | false -> med.cost2, med.recost2, chrom_map


let create_single_map med =
    let str = string_of_int in  
    let seq_arr = Array.map 
        (fun m ->
             let p_ref_code, p_seq_ord = (str med.ref_code2), (str m.seq_ord2) in

             let c_ref_code, c_seq_ord = 
                 (str med.ref_code1), (str m.seq_ord1)
             in 

             let attributes = [(Tags.GenomeMap.a_ref_code, p_ref_code);
                               (Tags.GenomeMap.a_seq_order, p_seq_ord);
                               (Tags.GenomeMap.a_dir_seg, "+");
                               (Tags.GenomeMap.d_ref_code, c_ref_code); 
                               (Tags.GenomeMap.d_seq_order, c_seq_ord);
                               (Tags.GenomeMap.d_dir_seg, "+")
                              ] 
             in 
             let m : Tags.output = (Tags.GenomeMap.seg, attributes, `String "") in 
             `Single m
        ) med.seq_arr
    in 


    let chrom_map : Tags.output = 
        (Tags.GenomeMap.chrom, [], `Structured (`Set  (Array.to_list seq_arr))) 
    in 
    chrom_map


let to_single single_parent child_ref c2 =
    let gap = Cost_matrix.Two_D.gap c2 in 
    let map = 
        Array.map (fun m ->
                       let alied_single_seq, alied_child_seq  =
                          match child_ref = single_parent.ref_code1 with
                          | true ->
                                let single, _ = UtlPoy.closest_alied_seq
                                    m.alied_med m.alied_seq1 c2
                                in 
                                single, m.alied_seq1
                          | false ->
                                let single, _ = UtlPoy.closest_alied_seq
                                    m.alied_med m.alied_seq2 c2
                                in 
                                single, m.alied_seq2
                      in 
                      let order = 
                          match child_ref = single_parent.ref_code1 with
                          | true -> m.seq_ord1
                          | false -> m.seq_ord2
                      in 
 
                      let ungap_alied_med = Sequence.fold_righti 
                          (fun ungap_alied_med p code ->
                               match code = gap with
                               | false -> code::ungap_alied_med
                               | true ->
                                     if Sequence.get alied_child_seq p = gap then ungap_alied_med
                                     else code::ungap_alied_med
                          ) [] alied_single_seq
                      in                       
  
                      let ungap_alied_med = UtlPoy.of_array (Array.of_list  ungap_alied_med) in                      
                      ungap_alied_med, order
                 ) single_parent.seq_arr
    in 

    Array.sort (fun seg1 seg2 ->
                    let _,  ord1 = seg1 in 
                    let _,  ord2 = seg2 in
                    ord1 - ord2
               ) map;

    
    let seq_ls = Array.fold_right 
        (fun seg seq_ls -> 
             let seq,  ord = seg in 
             if ord < 0 then seq_ls
             else seq::seq_ls
        ) map []
    in 

    Array.of_list seq_ls





let to_single_root root other_code c2 =
    let gap = Cost_matrix.Two_D.gap c2 in 
    let map = 
        Array.map (fun m ->
                       let child_alied_seq, other_alied_seq = 
                           match root.ref_code2 = other_code with 
                           | true -> m.alied_seq1, m.alied_seq2
                           | false -> m.alied_seq1, m.alied_seq2 
                       in 
                       let single_seq, _ = UtlPoy.closest_alied_seq
                           other_alied_seq child_alied_seq c2 
                       in 
                      let ungap_alied_med = Sequence.fold_righti 
                          (fun ungap_alied_med p code ->
                               match code = gap with
                               | false -> code::ungap_alied_med
                               | true ->
                                     if Sequence.get m.alied_med p = gap then ungap_alied_med
                                     else code::ungap_alied_med
                          ) [] single_seq
                      in                       
  
                      let ungap_alied_med = UtlPoy.of_array (Array.of_list  ungap_alied_med) in                      
                      ungap_alied_med
                 ) root.seq_arr
    in 
    map
    



let change_to_single med single_seq_arr = 
    let gap = Alphabet.gap in 

    let new_seq_arr = Array.mapi 
        (fun idx m ->
             let single_seq = single_seq_arr.(idx) in
             let num_dna = ref 0 in 
             let single_alied_med = UtlPoy.map
                 (fun code ->
                      if code = gap then gap
                      else begin
                          let single_code = Sequence.get single_seq !num_dna in 
                          incr num_dna;
                          single_code
                      end 
                 ) m.alied_med
             in 
             {m with alied_med = single_alied_med;
                  seq = (UtlPoy.delete_gap single_alied_med) }
        ) med.seq_arr
    in 
    {med with seq_arr = new_seq_arr}




let to_formater med alph = 
    let seq_str_arr = 
        Array.map (fun seg -> Sequence.to_formater seg.seq alph) med.seq_arr 
    in 
    String.concat "|" (Array.to_list seq_str_arr)

let copy_chrom_map s d = {d with ref_code = s.ref_code; 
                              ref_code1 = s.ref_code1;
                              ref_code2 = s.ref_code2;
                              seq_arr = s.seq_arr}
