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

(** The implementation of funtions to calculate the cost, 
* alignments and medians between annotated chromosomes 
* where both point mutations and rearrangement operations 
* are considered *)
let fprintf = Printf.fprintf
type direction_t = ChromPam.direction_t

(** [seq_t] is data structure to contain a segment of an annotated chromosome *)
type seq_t = {
    seq : Sequence.s; (** the segment sequence *)
    seq_ref_code : int; (** the reference code of this segment *)
    alied_med : Sequence.s; (** the aligned sequence this segement *)

    seq_ord1 : int; (** the segment order of its first child *)
    alied_seq1 : Sequence.s; (** the aligned sequence of its first child *)
    seq_ord2 : int; (** the segment order of its second child *)
    alied_seq2 : Sequence.s; (** the aligned sequence of its second child *)
    dir2 : direction_t; (** the orientation of this segment on the second chromosome *) 
}

(** [annchrom_t] is a data structure to contain an annoated chromosome *)
type annchrom_t = {
    seq_arr : seq_t array;  (** annotated chromosome *)
    ref_code : int;  (** reference code of this chromosome *)
    ref_code1 : int;    (** first child's code *)    
    ref_code2 : int;  (** second child's code *)
    cost1 : int; (** cost from this chromosome to its first child *)
    recost1 : int; (** recost from this chromosome to its first child *)
    cost2 : int; (** cost from this chromosome to its second child *)
    recost2 : int; (** recost from this chromosome to its second child *) 
}

(** [clone_seq] returns a fresh clone of segment [s] *)
let clone_seq s = {
    seq = Sequence.clone s.seq;
    seq_ref_code = s.seq_ref_code;
    alied_med = Sequence.clone s.alied_med;
    seq_ord1 = s.seq_ord1;
    alied_seq1 = Sequence.clone s.alied_seq1;
    seq_ord2 = s.seq_ord2;
    alied_seq2 = Sequence.clone s.alied_seq2;
    dir2 = s.dir2;
}

(** [clone_med] returns a fresh clone of annotated chromosome [m] *)
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


(** [annchromPam_t] is data structure to 
* contains parameters used to align two annotated chromosomes *)
type annchromPam_t = {
    re_meth : Data.re_meth_t;
    keep_median : int;
    circular : int;
    swap_med : int;
    approx : ChromPam.order_t;
    symmetric : bool;
    locus_indel_cost : (int * int);
    kept_wag : int;
}

let annchromPam_default = {
    re_meth = `Locus_Breakpoint 10;
    keep_median = 1;
    circular = 0;
    swap_med = 1;
    approx = `BothSeq;
    symmetric = false;
    locus_indel_cost = (10, 100);
    kept_wag = 1;
}

(** [init_seq_t (seq, code)] returns a segment from [seq] and [code]*)
let init_seq_t (seq, code) = {
    seq = seq;
    seq_ref_code = code;
    alied_med = seq;
    seq_ord1 = -1;
    alied_seq1 = seq;
    seq_ord2 = -1;
    alied_seq2 = seq;
    dir2 = `Positive;
    
}

(** [init seq_arr] returns an annotated chromosome 
* from sequence array [seq_arr] *)
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
                    Sequence.printDNA m.alied_seq1;
                    Sequence.printDNA m.alied_seq2;
                    Sequence.printDNA m.alied_med;
               ) seq_arr;
    print_endline "End printing alied_annchrom ";
    print_newline ()

let printSeq x = Printf.printf "%s\n" (Sequence.to_string x Alphabet.nucleotides)


    
(** [swap_seq s] swaps the first child and second child of 
* this segment [s] *)
let swap_seq s =  {
    s with seq_ord1 = s.seq_ord2;
        alied_seq1 = s.alied_seq2;
        seq_ord2 = s.seq_ord1;
        alied_seq2 = s.alied_seq1                
} 



(** [swap_med m] swaps the first child and second child 
of this annotated chromosome [m]*)
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

(** [get_annchrom_pam] returns user defined parameters
* to align two annotated chromosomes *)
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

    let chrom_pam = 
        match user_annchrom_pam.Data.max_kept_wag with
        | None -> chrom_pam
        | Some max_kept_wag -> {chrom_pam with kept_wag = max_kept_wag}
    in 

    chrom_pam

let print annchrom alpha = 
    Array.iter (fun s -> 
                    let seq = Sequence.to_string s.seq alpha in 
                    fprintf stdout "(%i, %s) | " s.seq_ref_code seq) annchrom.seq_arr;
    print_newline ()

(** [split chrom] returns an array code sequences and an array
* of codes coresponding to the sequence array*)
let split chrom =  
    let seq_arr = Array.map (fun s -> s.seq) chrom.seq_arr in  
    let code_arr = Array.map (fun s -> s.seq_ref_code) chrom.seq_arr in  
    seq_arr, code_arr 



(** [create_pure_gen_cost_mat seq1_arr seq2_arr cost_mat ali_pam]
* returns a cost matrix between segments of [seq1_arr] and segments of [seq2_arr] *)
let create_pure_gen_cost_mat seq1_arr seq2_arr cost_mat ali_pam =        
    let seq1_arr = Array.mapi (fun ith seq -> seq, (ith * 2 + 1) ) seq1_arr in 

    let len1 = Array.length seq1_arr in 

    let seq2_arr = Array.mapi 
        (fun ith seq -> seq, (ith * 2 + 1 + len1 * 2) ) seq2_arr 
    in 

    let len2 = Array.length seq2_arr in 
    let gen_gap_code = (len1 + len2) * 2 + 1 in 


    let pure_gen_cost_mat = Array.make_matrix (gen_gap_code + 1) (gen_gap_code + 1)
        Utl.large_int 
    in 

    pure_gen_cost_mat.(gen_gap_code).(gen_gap_code) <- 0;
    
    let update (seq1, code1) (seq2, code2) =
        let com_seq1 = Sequence.complement_chrom Alphabet.nucleotides seq1 in 
        let com_seq2 = Sequence.complement_chrom Alphabet.nucleotides seq2 in 
            let _, _, cost, _ = Sequence.align2 seq1 seq2 cost_mat in 

        pure_gen_cost_mat.(code1).(code2) <- cost;
        pure_gen_cost_mat.(code2).(code1) <- cost;


        let _, _, cost, _ = Sequence.align2 seq1 com_seq2 cost_mat in 

        pure_gen_cost_mat.(code1).(code2 + 1) <- cost;
        pure_gen_cost_mat.(code2 + 1).(code1) <- cost;

        let _, _, cost, _ = Sequence.align2 com_seq1 seq2 cost_mat in 

        pure_gen_cost_mat.(code1 + 1).(code2) <- cost;
        pure_gen_cost_mat.(code2).(code1 + 1) <- cost;

        let _, _, cost, _ = Sequence.align2 com_seq1 com_seq2 cost_mat in 

        pure_gen_cost_mat.(code1 + 1).(code2 + 1) <- cost;
        pure_gen_cost_mat.(code2 + 1).(code1 + 1) <- cost;

    in 

        
    Array.iter (fun (seq1, code1) ->
                    Array.iter (fun (seq2, code2) -> 
                                    update (seq1, code1) (seq2, code2)
                               ) seq2_arr) seq1_arr;


    let update_gap (seq, code) = 

(*       pure_gen_cost_mat.(gen_gap_code).(code) <- Sequence.cmp_locus_indel_cost seq cost_mat ali_pam.locus_indel_cost; *)
        pure_gen_cost_mat.(gen_gap_code).(code)<- Sequence.cmp_gap_cost ali_pam.locus_indel_cost seq;


        pure_gen_cost_mat.(code).(gen_gap_code) <- pure_gen_cost_mat.(gen_gap_code).(code);
        pure_gen_cost_mat.(gen_gap_code).(code + 1) <- pure_gen_cost_mat.(gen_gap_code).(code);
        pure_gen_cost_mat.(code + 1).(gen_gap_code) <- pure_gen_cost_mat.(gen_gap_code).(code);
    in 

    Array.iter update_gap seq1_arr; 
    Array.iter update_gap seq2_arr;
    let code1_arr = Array.map (fun (seq, code) -> code) seq1_arr in 
    let code2_arr = Array.map (fun (seq, code) -> code) seq2_arr in 
        pure_gen_cost_mat, code1_arr, code2_arr, gen_gap_code





(** [create_pure_gen_cost_mat_3 seq1_arr seq2_arr seq3_arr seqm_arr c2 ali_pam]
* creates three general cost matrices corresponding to
* the distances from [seq1_arr], [seq2_arr] and [seq3_arr] to [seqm_arr] *)
let create_pure_gen_cost_mat_3 seq1_arr seq2_arr seq3_arr seqm_arr c2 ali_pam =        
    let code = ref 0 in 
    let codem_arr = Array.map (fun seq -> code := !code + 2; seq, (!code - 1)) seqm_arr in
    let code1_arr = Array.map (fun seq -> code := !code + 2; seq, (!code - 1)) seq1_arr in
    let code2_arr = Array.map (fun seq -> code := !code + 2; seq, (!code - 1)) seq2_arr in
    let code3_arr = Array.map (fun seq -> code := !code + 2; seq, (!code - 1)) seq3_arr in
    let gen_gap_code = !code + 1 in 


    let cost1_mat = Array.make_matrix (gen_gap_code + 1) (gen_gap_code + 1) Utl.large_int  in 
    let cost2_mat = Array.make_matrix (gen_gap_code + 1) (gen_gap_code + 1) Utl.large_int  in 
    let cost3_mat = Array.make_matrix (gen_gap_code + 1) (gen_gap_code + 1) Utl.large_int  in 

    cost1_mat.(gen_gap_code).(gen_gap_code) <- 0;
    cost2_mat.(gen_gap_code).(gen_gap_code) <- 0;
    cost3_mat.(gen_gap_code).(gen_gap_code) <- 0;
    
    let update cost_mat (seq1, code1) (seq2, code2) =
        let com_seq1 = Sequence.complement_chrom Alphabet.nucleotides seq1 in 
        let com_seq2 = Sequence.complement_chrom Alphabet.nucleotides seq2 in 
        let _, _, cost, _ = Sequence.align2 seq1 seq2 c2 in 
        cost_mat.(code1).(code2) <- cost;
        cost_mat.(code2).(code1) <- cost;

        let _, _, cost, _ = Sequence.align2 seq1 com_seq2 c2 in 
        cost_mat.(code1).(code2 + 1) <- cost;
        cost_mat.(code2 + 1).(code1) <- cost;

        let _, _, cost, _ = Sequence.align2 com_seq1 seq2 c2 in 
        cost_mat.(code1 + 1).(code2) <- cost;
        cost_mat.(code2).(code1 + 1) <- cost;

        let _, _, cost, _ = Sequence.align2 com_seq1 com_seq2 c2 in 
        cost_mat.(code1 + 1).(code2 + 1) <- cost;
        cost_mat.(code2 + 1).(code1 + 1) <- cost;
    in 

        
    Array.iter 
        (fun (seqm, codem) ->
             Array.iter (fun (seq1, code1) -> 
                             update cost1_mat (seqm, codem) (seq1, code1)
                        ) code1_arr;

             Array.iter (fun (seq2, code2) -> 
                             update cost2_mat (seqm, codem) (seq2, code2)
                        ) code2_arr;

             Array.iter (fun (seq3, code3) -> 
                             update cost3_mat (seqm, codem) (seq3, code3)
                        ) code3_arr;
        ) codem_arr;


    let update_gap cost_mat (seq, code) = 
        cost_mat.(gen_gap_code).(code) <- Sequence.cmp_gap_cost ali_pam.locus_indel_cost seq;
        cost_mat.(code).(gen_gap_code) <- cost_mat.(gen_gap_code).(code);
        cost_mat.(gen_gap_code).(code + 1) <- cost_mat.(gen_gap_code).(code);
        cost_mat.(code + 1).(gen_gap_code) <- cost_mat.(gen_gap_code).(code);
    in 

    Array.iter (update_gap cost1_mat) codem_arr;
    Array.iter (update_gap cost1_mat) code1_arr;

    Array.iter (update_gap cost2_mat) codem_arr;    
    Array.iter (update_gap cost2_mat) code2_arr;

    Array.iter (update_gap cost3_mat) codem_arr;    
    Array.iter (update_gap cost3_mat) code3_arr;




    let codem_arr = Array.map (fun (seq, code) -> code) codem_arr in 
    let code1_arr = Array.map (fun (seq, code) -> code) code1_arr in 
    let code2_arr = Array.map (fun (seq, code) -> code) code2_arr in 
    let code3_arr = Array.map (fun (seq, code) -> code) code3_arr in 

    cost1_mat, cost2_mat, cost3_mat, code1_arr, code2_arr, code3_arr, codem_arr, gen_gap_code


(** [cmp_simple_cost chrom1 chrom2 cost_mat alpha ali_pam]
* computes the cost between annotated chromosome [chrom1]
* and annotated chromosome [chrom2] allowing rearrangements *)
let cmp_simple_cost (chrom1: annchrom_t) (chrom2 : annchrom_t) 
        cost_mat alpha ali_pam =


    let chrom_len1 = Array.fold_left (fun len s -> len + Sequence.length s.seq) 0 chrom1.seq_arr in     
    let chrom_len2 = Array.fold_left (fun len s -> len + Sequence.length s.seq) 0 chrom2.seq_arr in 
(*Printf.printf "chrom_len1 = %d, chrom_len2 = %d \n" chrom_len1 chrom_len2;*)
    if (chrom_len1 < 2) || (chrom_len2 < 2) then 0, 0
    else begin
        let seq1_arr, _ = split chrom1 in  
        let seq2_arr, _ = split chrom2 in  
(*let printseq x = Printf.printf "%s \n" (Sequence.to_string x Alphabet.nucleotides) in
Printf.printf "seq1_arr =\n"; Array.iter printseq seq1_arr;
Printf.printf "seq2_arr =\n"; Array.iter printseq seq2_arr;*)
        let pure_gen_cost_mat, code1_arr, code2_arr, gen_gap_code = 
            create_pure_gen_cost_mat seq1_arr seq2_arr cost_mat ali_pam  
        in
        let total_cost, (recost1, recost2), _, _ = 
            GenAli.create_gen_ali_code ali_pam.kept_wag `Annotated code1_arr code2_arr 
                pure_gen_cost_mat gen_gap_code  
                ali_pam.re_meth ali_pam.swap_med 
                ali_pam.circular true
        in 
        total_cost, (recost1 + recost2)    
    end 


(** [cmp_cost chrom1 chrom2 cost_mat alpha annchrom_pam] 
* compute the minimum cost from annotated chromosome [chrom1]
* to annotated chromosome [chrom2], and from [chrom2] to [chrom1] *)
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

              
(** [find_simple_med2_ls chrom1 chrom2 cost_mat alpha ali_pam]
* finds medians between [chrom1] and [chrom2] allowing rearrangements *)
let find_simple_med2_ls (chrom1: annchrom_t) (chrom2 : annchrom_t) 
        (cost_mat : Cost_matrix.Two_D.m) alpha ali_pam =

(*Printf.printf "annchromAli.ml, find_simple_med2_ls between chrom1=%d and chrom2=%d.\n" chrom1.ref_code chrom2.ref_code;
*)
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
        (*Printf.printf "seq1_arr = \n";
        Array.iter printSeq seq1_arr;
        Printf.printf "seq2_arr = \n"; Array.iter printSeq seq2_arr;i*)
        let total_cost, (recost1, recost2), alied_code1_arr, alied_code2_arr = 
            GenAli.create_gen_ali_code ali_pam.kept_wag `Annotated code1_arr code2_arr 
                pure_gen_cost_mat gen_gap_code  
                ali_pam.re_meth ali_pam.swap_med 
                ali_pam.circular  true
        in 
  (*       Printf.printf "code1_arr and code2_arr =\n";
         Array.iter (Printf.printf "%d ") code1_arr;
         print_newline();
         Array.iter (Printf.printf "%d ") code2_arr;
         print_newline();
        Printf.printf "total_cost = %d, alied_code1 and alied_code2 = \n" total_cost;
        Array.iter (Printf.printf "%d ") alied_code1_arr;
         print_newline();
        Array.iter (Printf.printf "%d ") alied_code2_arr;
        print_newline();
  *)
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
                 (*Printf.printf "pos_code2 = %d, " pos_code2;*)
                 let idx2 = Utl.find_index code2_arr pos_code2 compare in  
                 let seq2, dir2 =
                     match idx2 with 
                     | -1 -> None, `Positive
                     | _ -> 
                           if pos_code2 = alied_code2_arr.(idx) then 
                               begin
                                Some chrom2.seq_arr.(idx2).seq, `Positive
                               end
                           else
                               begin
                               Some (Sequence.complement_chrom
                           Alphabet.nucleotides chrom2.seq_arr.(idx2).seq),`Negative
                               end
                 in                            
                 let alied_med_seq, alied_seq1, alied_seq2 =
                     match seq1, seq2 with   
                     | Some seq1, Some seq2 -> 
                           let med, alied_seq1, alied_seq2, _  =
                               Sequence.create_median  ~approx:approx seq1 seq2 cost_mat  
                           in
                         med, alied_seq1, alied_seq2
                     | Some seq1, None -> 
                           let len1 = Sequence.length seq1 in 
                           let seq2 = Sequence.create_gap_seq len1 in 
                           let med, _ = Sequence.create_median_seq ~approx:approx seq1 seq2 cost_mat in 
                           med, seq1, seq2
                     | None, Some seq2 -> 
                           let len2 = Sequence.length seq2 in 
                           let seq1 = Sequence.create_gap_seq len2 in 
                           let med, _ = Sequence.create_median_seq ~approx:approx seq1 seq2 cost_mat in
                           med, seq1, seq2
                     | _, _ -> Sequence.get_empty_seq (), Sequence.get_empty_seq (), Sequence.get_empty_seq ()
                 in
                                    
                  (*Printf.printf "alied_seq1/alied_seq2=\n"; printSeq alied_seq1;
                 printSeq alied_seq2; *)
                 (alied_med_seq, idx1, alied_seq1, idx2, alied_seq2, dir2)
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
                         let alied_med, seq_ord1, alied_seq1, seq_ord2, alied_seq2, dir2 =  ali_chrom.(index) in
                         {seq= Sequence.delete_gap alied_med;  
                          seq_ref_code = Utl.get_new_seq_ref_code();  
                          alied_med = alied_med;
                          seq_ord1 = seq_ord1; 
                          alied_seq1 = alied_seq1;
                          seq_ord2 = seq_ord2;
                          alied_seq2 = alied_seq2;
                          dir2 = dir2;
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
            if (Utl.isEqualArr code2_arr re_code2_arr compare) ||  
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
        total_cost, (recost1 + recost2),med_ls           
    end




(** [find_med2_ls chrom1 chrom2 cost_mat alpha annchrom_pam]
* finds medians between [chrom1] and [chrom2] and vice versa
* allowing rearrangements *)
let find_med2_ls (chrom1: annchrom_t) (chrom2 : annchrom_t) 
        (cost_mat : Cost_matrix.Two_D.m) alpha annchrom_pam = 
    let ali_pam = get_annchrom_pam annchrom_pam in          
    match ali_pam.symmetric with
    | true ->
    (*          Printf.printf "find_med2_ls, get cost12:\n" ;*)
          let cost12, recost12, med12_ls = find_simple_med2_ls chrom1 chrom2
              cost_mat alpha ali_pam in 
          let ali_pam = 
              if ali_pam.approx = `First then {ali_pam with approx = `Second}
              else ali_pam
          in 
        (* Printf.printf "cost12 = %d, get cost21 => \n" cost12;*)
           let cost21, recost21, med21_ls = find_simple_med2_ls chrom2 chrom1
              cost_mat alpha ali_pam in 
         (* Printf.printf "cost21 = %d \n" cost21;*)
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


(** [find_med3 ch1 ch2 ch3 mine c2 c3 alpha annchrom_pam]
* finds the median of annotated chromosome [ch1], [ch2],
* and [ch3] based on the current median [mine] *) 
let find_med3 ch1 ch2 ch3 mine c2 c3 alpha annchrom_pam = 
    let ali_pam = get_annchrom_pam annchrom_pam in 
    let seq1_arr, _ = split ch1 in 
    let seq2_arr, _ = split ch2 in 
    let seq3_arr, _ = split ch3 in 
    let seqm_arr, _ = split mine in 

    let cost1_mat, cost2_mat, cost3_mat, code1_arr, code2_arr, code3_arr, codem_arr, gen_gap_code = 
        create_pure_gen_cost_mat_3 seq1_arr seq2_arr seq3_arr seqm_arr c2 ali_pam
    in 


    let _, _, alied_code1_arr, alied_code1m_arr = 
        GenAli.create_gen_ali_code ali_pam.kept_wag `Annotated code1_arr codem_arr cost1_mat
            gen_gap_code ali_pam.re_meth ali_pam.swap_med ali_pam.circular true
    in 


    let _, _, alied_code2_arr, alied_code2m_arr = 
        GenAli.create_gen_ali_code ali_pam.kept_wag `Annotated code2_arr codem_arr cost2_mat
            gen_gap_code ali_pam.re_meth ali_pam.swap_med ali_pam.circular true
    in 


    let _, _, alied_code3_arr, alied_code3m_arr = 
        GenAli.create_gen_ali_code ali_pam.kept_wag `Annotated code3_arr codem_arr cost3_mat
            gen_gap_code ali_pam.re_meth ali_pam.swap_med ali_pam.circular true
    in 

    let seq_arr = Array.mapi 
        (fun ith seqt -> 
             let seqm = seqt.seq in 
             let codem = codem_arr.(ith) in 

             let comp v1 v2 = (v1 + 1) / 2 -  (v2 + 1) / 2 in


             let get_seq seq_arr code_arr alied_code_arr alied_codem_arr  = 
                 let indexm = Utl.find_index alied_codem_arr codem comp in 
                 let code = alied_code_arr.(indexm) in 
                 let index = Utl.find_index code_arr code comp in 
                 match index with 
                 | -1 -> Sequence.get_empty_seq ()
                 | _ -> seq_arr.(index).seq 
             in 

             let seq1 = get_seq ch1.seq_arr code1_arr alied_code1_arr
                 alied_code1m_arr 
             in 

             let seq2 = get_seq ch2.seq_arr code2_arr alied_code2_arr
                 alied_code2m_arr 
             in 

             let seq3 = get_seq ch3.seq_arr code3_arr alied_code3_arr
                 alied_code3m_arr 
             in 
             
             if (Sequence.length seq1 = 0) ||  (Sequence.length seq2 = 0) ||
                 (Sequence.length seq3 = 0) then {seqt with seq = seqm}
             else begin
                 let _, median_seq, _ = Sequence.Align.readjust_3d
                     ~first_gap:false seq1 seq2 seqm c2 c3 seq3 in
                 {seqt with seq = median_seq}
             end               
        ) mine.seq_arr
    in 

    let cost1, recost1 = cmp_cost ch1 mine c2 alpha annchrom_pam in 
    let cost2, recost2 = cmp_cost ch2 mine c2 alpha annchrom_pam in 
    let cost3, recost3 = cmp_cost ch3 mine c2 alpha annchrom_pam in 
    let total_cost = cost1 + cost2 + cost3 in

    let mine3 = {mine with seq_arr = seq_arr} in 
    let cost1, recost1 = cmp_cost ch1 mine3 c2 alpha annchrom_pam in 
    let cost2, recost2 = cmp_cost ch2 mine3 c2 alpha annchrom_pam in 
    let cost3, recost3 = cmp_cost ch3 mine3 c2 alpha annchrom_pam in 
    let total_cost3 = cost1 + cost2 + cost3 in
(*    fprintf stdout "%i -> %i\n" total_cost total_cost3; flush stdout;*)
    if total_cost3 < total_cost then total_cost3, mine3
    else total_cost, mine






(** [compare annchrom1 annchrom2] compares
* annotated chromosome [annchrom1] and annotated chromosome [annchrom2] *)
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


(** [find_approx_med2 med1 med2 med12] returns an approximate
* median between annotated chromosomes [med1] and [med2] based 
* on computed median [med12] *)
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



(** [assign_seq_ref annchrom seq_ref_code] assigns
* ID for all segments of annotated chromosome [annchrom] *)
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


(** [create_map med child_ref] returns the map
* between chromosome [med] and its child [child_ref] *)
let create_map med child_ref = 
    let str x = `Int x in
    let seq_arr = Array.mapi 
        (fun med_id m ->
             let p_ref_code, p_seq_ord = (str med.ref_code), (str med_id) in
             let c_ref_code, c_seq_ord = 
                 match child_ref = med.ref_code1 with
                 | true -> (str med.ref_code1), (str m.seq_ord1)
                 | false -> (str med.ref_code2), (str m.seq_ord2)
             in 
             let attributes = [(Xml.GenomeMap.a_ref_code, p_ref_code);
                               (Xml.GenomeMap.a_seq_order, p_seq_ord);
                               (Xml.GenomeMap.a_dir_seg, `String "+");
                               (Xml.GenomeMap.d_ref_code, c_ref_code); 
                               (Xml.GenomeMap.d_seq_order, c_seq_ord);
                               (Xml.GenomeMap.d_dir_seg, `String "+")
                              ] 
             in 
             let m : Xml.xml = (Xml.GenomeMap.seg, attributes, `String "") in 
             `Single m
        ) med.seq_arr
    in 


    let chrom_map : Xml.xml = 
        (Xml.GenomeMap.chrom, [], (`Set  (Array.to_list seq_arr))) 
    in 
    match child_ref = med.ref_code1 with
    | true -> med.cost1, med.recost1, chrom_map
    | false -> med.cost2, med.recost2, chrom_map


let create_single_map med =
    let str x = `Int x in  
    let seq_arr = Array.map 
        (fun m ->
             let p_ref_code, p_seq_ord = (str med.ref_code2), (str m.seq_ord2) in
             let c_ref_code, c_seq_ord = 
                 (str med.ref_code1), (str m.seq_ord1)
             in 

             let attributes = [(Xml.GenomeMap.a_ref_code, p_ref_code);
                               (Xml.GenomeMap.a_seq_order, p_seq_ord);
                               (Xml.GenomeMap.a_dir_seg, `String "+");
                               (Xml.GenomeMap.d_ref_code, c_ref_code); 
                               (Xml.GenomeMap.d_seq_order, c_seq_ord);
                               (Xml.GenomeMap.d_dir_seg, `String "+")
                              ] 
             in 
             let m : Xml.xml = (Xml.GenomeMap.seg, attributes, `String "") in 
             `Single m
        ) med.seq_arr
    in 


    let chrom_map : Xml.xml = 
        (Xml.GenomeMap.chrom, [], (`Set  (Array.to_list seq_arr))) 
    in 
    chrom_map

let get_alied_seq_code parent_chrom child_chrom c2 annchrom_pam =
    let ali_pam = get_annchrom_pam annchrom_pam in 
    let seq1_code, _ = split parent_chrom in  
    let seq2_code, _ = split child_chrom in  
    let pure_gen_cost_mat, code1_arr, code2_arr, gen_gap_code = 
            create_pure_gen_cost_mat seq1_code seq2_code c2 ali_pam  
    in
    let _, (_, _), ali_seq1_code, ali_seq2_code = 
         GenAli.create_gen_ali_code ali_pam.kept_wag `Annotated code1_arr code2_arr 
                pure_gen_cost_mat gen_gap_code  
                ali_pam.re_meth ali_pam.swap_med 
                ali_pam.circular true
    in 
    (code1_arr, code2_arr, ali_seq1_code,ali_seq2_code)


let compare_int looking_item array_item =
    if(looking_item = array_item) then 0
    else 1


let seqcode_to_seq chrom1 chrom2 code1_arr code2_arr alied_code1_arr alied_code2_arr =
    let full_len = Array.length alied_code1_arr in 
    let full_seq_t_pair_arr = Array.init full_len
        (fun idx ->   
         let idx1 = Utl.find_index code1_arr alied_code1_arr.(idx) compare_int in   
         let seq1 =  
             match idx1 with  
             | -1 -> None
             | _ -> Some chrom1.seq_arr.(idx1).seq
         in                             
         let pos_code2 = 
                 if alied_code2_arr.(idx) mod 2 = 0 then alied_code2_arr.(idx) - 1
                 else alied_code2_arr.(idx)
         in 
         let idx2 = Utl.find_index code2_arr pos_code2 compare_int in  
         let seq2 =
                 match idx2 with 
                 | -1 -> None
                 | _ -> 
                       if pos_code2 = alied_code2_arr.(idx) then 
                            Some chrom2.seq_arr.(idx2).seq
                       else Some (Sequence.complement_chrom Alphabet.nucleotides chrom2.seq_arr.(idx2).seq)
         in 
         let full_seq_t1, full_seq_t2 =
                 match seq1, seq2 with   
                 | Some seq1, Some seq2 -> 
                         chrom1.seq_arr.(idx1),chrom2.seq_arr.(idx2) 
                 | Some seq1, None -> 
                       let len1 = Sequence.length seq1 in 
                       let seq2 = Sequence.create_gap_seq len1 in
                       let tmp_seq_t = { chrom2.seq_arr.(0) with seq = seq2;
                       alied_med = seq2; alied_seq1 = seq2; alied_seq2 = seq2 }  
                       in
                       chrom1.seq_arr.(idx1), tmp_seq_t
                 | None, Some seq2 -> 
                       let len2 = Sequence.length seq2 in 
                       let seq1 = Sequence.create_gap_seq len2 in 
                       let tmp_seq_t = { chrom1.seq_arr.(0) with seq = seq1;
                       alied_med = seq1; alied_seq1 = seq1; alied_seq2 = seq1 }  
                       in
                       tmp_seq_t, chrom2.seq_arr.(idx2)
                 | _, _ -> 
                         let empty_seq = Sequence.get_empty_seq () in
                         let tmp_seq_t1 = { chrom1.seq_arr.(0) with seq = empty_seq;
                         alied_med = empty_seq; alied_seq1 = empty_seq; 
                         alied_seq2 = empty_seq}  
                         in
                         let tmp_seq_t2 = { chrom2.seq_arr.(0) with seq = empty_seq;
                         alied_med = empty_seq; alied_seq1 = empty_seq;
                         alied_seq2 = empty_seq }  
                         in
                         tmp_seq_t1, tmp_seq_t2
         in
         (full_seq_t1, full_seq_t2)
        )
    in
    let full_seq_t1_arr = Array.init full_len
    (fun idx ->
           match full_seq_t_pair_arr.(idx) with
           | (full_seq_t1, full_seq_t2) ->  full_seq_t1
    ) in
    let full_seq_t2_arr = Array.init full_len
    (fun idx ->
            match full_seq_t_pair_arr.(idx) with
            | (full_seq_t1, full_seq_t2) -> full_seq_t2
    ) in       
    {chrom1 with seq_arr = full_seq_t1_arr},{ chrom2 with seq_arr = full_seq_t2_arr}    
    

    (** [to_single_root root other_code c2] returns the single
* state sequence for the [root] *)
let to_single_root root other_code c2 =
    match (root.ref_code1 = -1) && (root.ref_code2 = -1) with
    | true -> Array.map (fun s -> s.seq) root.seq_arr
    | false ->   
          let gap = Cost_matrix.Two_D.gap c2 in 
          let map = 
              Array.map 
                  (fun m ->
                    let child_alied_seq, other_alied_seq =  m.alied_seq1, m.alied_seq2 in 
                       (*Printf.printf "At root, alied_seq1/alied_seq2 =\n";
                       printSeq child_alied_seq; printSeq other_alied_seq;*)
                       let single_seq, _ = Sequence.closest_alied_seq
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
                      let ungap_alied_med = Sequence.of_array (Array.of_list  ungap_alied_med) 
                      in    
                      ungap_alied_med
                 ) root.seq_arr
          in 
          map 
    


let assign_single_nonroot parent child child_ref c2 annchrom_pam = 
    let parent_seq_code, child_seq_code, parent_ali_seq_code, child_ali_seq_code = 
        get_alied_seq_code parent child c2 annchrom_pam
    in
    (* though we don't need alied sequence here, we do need to fill in the gaps
    * for the child sequence, so that two sequences have the same size *)
    let parent_with_full_seq, child_with_full_seq = seqcode_to_seq 
    parent child parent_seq_code child_seq_code parent_ali_seq_code
    child_ali_seq_code in  
    let parent = parent_with_full_seq in 
    let child = child_with_full_seq in
    match (parent.ref_code1 = -1) && (parent.ref_code2 = -1) with
    | true -> Array.map (fun s -> s.seq) parent.seq_arr
    | false -> 
          let gap = Cost_matrix.Two_D.gap c2 in 
          let map = 
              Array_ops.map_2 
              (fun parent_seq_t child_seq_t ->
                let alied_single_seq,_ = 
                  Sequence.Align.closest parent_seq_t.alied_med child_seq_t.seq c2 Matrix.default
                in
                 (*Printf.printf "child's seq = "; printSeq child_seq_t.seq;
                 Printf.printf "parent's alied_seq "; printSeq parent_seq_t.alied_med;  *)
                let alied_child_seq,order  =
                   match child_ref = parent.ref_code1 with
                  | true -> 
                          parent_seq_t.alied_seq1,parent_seq_t.seq_ord1
                  | false ->
                          parent_seq_t.alied_seq2,parent_seq_t.seq_ord2
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
                let ungap_alied_med = Sequence.of_array (Array.of_list  ungap_alied_med) in
                ungap_alied_med, order
            ) parent.seq_arr child.seq_arr
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
              
(*    
(** [to_single single_parent child_ref c2] returns
* thesingle state sequence for annotated chromosome [child_ref] *)
let to_single single_parent child_ref c2 =
    Printf.printf "annchromAli.ml to_single, parent's ref = %d, two children's
    ref= %d/%d, child's ref = %d\n %!" single_parent.ref_code single_parent.ref_code1 single_parent.ref_code2 child_ref;
    match (single_parent.ref_code1 = -1) && (single_parent.ref_code2 = -1) with
    | true -> Array.map (fun s -> s.seq) single_parent.seq_arr
    | false -> 
          let gap = Cost_matrix.Two_D.gap c2 in 
          let map = 
              Array.map 
                  (fun m ->
                       let alied_single_seq, alied_child_seq  =
                           match child_ref = single_parent.ref_code1 with
                          | true ->
                                 let single, _ = Sequence.closest_alied_seq
                                    m.alied_med m.alied_seq1 c2
                                in 
                                single, m.alied_seq1
                          | false ->
                                 let single, _ = Sequence.closest_alied_seq
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
                      let ungap_alied_med = Sequence.of_array (Array.of_list  ungap_alied_med) in
                      Printf.printf "+++alied_single_seq/ungap_alied_med=\n"; 
                      printSeq alied_single_seq; printSeq ungap_alied_med;
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
*)              
              
(** [change_to_single med single_seq_arr c2] changes the
* single state sequence of [med] by [single_seq] *)
let change_to_single med single_seq_arr c2 = 
    let gap = Alphabet.gap in 
    let new_seq_arr = Array.mapi 
        (fun idx m ->
             let single_seq = single_seq_arr.(idx) in
             let len = Sequence.length single_seq in               
             let num_dna = ref 0 in 
             let num_not_gap = Sequence.cmp_num_not_gap m.alied_med in
             let single_alied_med = 
                if len = num_not_gap then begin
                    Sequence.map
                        (fun code ->    
                            if code = gap then gap
                            else begin
                                let single_code = Sequence.get single_seq !num_dna in 
                                incr num_dna;
                                single_code
                            end 
                         ) m.alied_med
                end else Sequence.get_single_seq m.alied_med c2
             in 

             {m with alied_med = single_alied_med;
                  seq = (Sequence.delete_gap single_alied_med) }
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

(** [convert_map med] converts the map of the 
* median [med] into the format suitable for
* creating implied alignments *)
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


