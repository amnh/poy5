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

let () = SadmanOutput.register "ChromAli" "$Revision: 2403 $"

(** The implementation of funtions to calculate the cost, alignments and medians
    between chromosomes where both point mutations and rearrangement operations
    are considered *)
type seed_t = Seed.seed_t
type pairChromPam_t = ChromPam.chromPairAliPam_t
type block_t = Block.block_t
type subseq_t = Subseq.subseq_t
let  fprintf = Printf.fprintf
let  deref = Utl.deref
let  debug = false

type direction_t = ChromPam.direction_t

    
type seg_t = {
    sta : int;
    en : int;
    cost : int; 
    alied_med : Sequence.s;

    sta1 : int;
    en1 : int;
    alied_seq1 : Sequence.s;
    dir1 : direction_t;

    sta2 : int;
    en2 : int;
    alied_seq2 : Sequence.s; (** alied_seq1 <-> reversed alied_seq2 
                                if dir2 is Negative *)
    dir2 : direction_t;
}

type med_t = {
    seq : Sequence.s;
    ref_code : int;
    ref_code1 : int;    (** Child's code *)    
    ref_code2 : int;  (** Child's code *)
    cost1 : int;
    recost1 : int;
    cost2 : int;
    recost2 : int;
    chrom_map : seg_t list;    
}


let create_med seq = {
    seq = seq; 
    ref_code = -1; 
    ref_code1 = -1;
    ref_code2 = -1;
    chrom_map = [];
    cost1 = 0;
    recost1 = 0;
    cost2 = 0;
    recost2 =0;
}

let init_med seq = 
    let med = {seq = seq; 
               ref_code = Utl.get_new_chrom_ref_code ();  
               ref_code1 = -1;
               ref_code2 = -1;
               chrom_map =[];
               cost1 = 0;
               recost1 = 0;
               cost2 = 0;
               recost2 =0;
              } 
    in  
    med


let clone_seg s = {
    sta = s.sta;
    en = s.en;
    cost = s.cost;
    alied_med = Sequence.clone s.alied_med;
    
    sta1 = s.sta1;
    en1 = s.en1;
    alied_seq1 = Sequence.clone s.alied_seq1;
    dir1 = s.dir1;

    sta2 = s.sta2;
    en2 = s.en2;
    alied_seq2 = Sequence.clone s.alied_seq2;
    dir2 = s.dir2;    
}


let clone_med m = {
    seq = Sequence.clone m.seq;
    ref_code = m.ref_code;
    ref_code1 = m.ref_code1;
    ref_code2 = m.ref_code2;
    cost1 = m.cost1;
    cost2 =m.cost2;
    recost1 = m.recost1;
    recost2 = m.recost2;
    chrom_map = List.map clone_seg m.chrom_map
}


let swap_seg s = 
    {s with sta1 = s.sta2; en1 = s.en2; alied_seq1 = s.alied_seq2; dir1 = s.dir2;
         sta2 = s.sta1; en2 = s.en1; alied_seq2 = s.alied_seq1; dir2 = s.dir1}
            
let swap_chrom_map m = 
    List.map swap_seg m

let swap_med m = 
    {m with ref_code1 = m.ref_code2;
         ref_code2 = m.ref_code1;
         cost1 = m.cost2;
         recost1 = m.recost2;
         cost2 = m.cost1;
         recost2 = m.recost1;
         chrom_map = swap_chrom_map m.chrom_map}


let print med =
    fprintf stdout "%i -> %i %i\n " med.ref_code
        med.ref_code1 med.ref_code2;
    flush stdout  


let print_map map = 
    print_endline "Start the chromosome map";
    let f = stdout in
    List.iter 
        (fun seg ->
                 fprintf f "length: %i, cost: %i, " (Sequence.length seg.alied_seq1) seg.cost;
                 fprintf f "sta: %i, end: %i, " seg.sta seg.en;
                 (match seg.dir1 with
                 | `Positive -> fprintf f "dir1: positive, "
                 | _ -> fprintf f "dir1: negative, ");

                 fprintf f "sta1: %i, end1: %i, " seg.sta1 seg.en1;
                 (match seg.dir1 with
                 | `Positive -> fprintf f "dir1: positive, "
                 | _ -> fprintf f "dir1: negative, ");

                 fprintf f "sta2: %i, end2: %i, " seg.sta2 seg.en2;
                 (match seg.dir2 with
                 | `Positive -> fprintf f "dir2: positive\n"
                 | _ -> fprintf f "dir2: negative\n");

                 fprintf f "%s\n" (Sequence.to_string seg.alied_seq1 Alphabet.nucleotides);
                 fprintf f "%s\n" (Sequence.to_string seg.alied_seq2 Alphabet.nucleotides);
                 fprintf f "%s\n" (Sequence.to_string seg.alied_med Alphabet.nucleotides);
                 
        ) map;
    print_endline "End of chromosome map";
    print_newline ()

  

let print_median med_ls outfile = 
    let f = open_out outfile in 

    let print_med med = 
        List.iter 
            (fun seg ->
                 fprintf f "length: %i, cost: %i, " (Sequence.length seg.alied_seq1) seg.cost;
                 fprintf f "sta: %i, end: %i, " seg.sta seg.en;
                 (match seg.dir1 with
                 | `Positive -> fprintf f "dir1: positive, "
                 | _ -> fprintf f "dir1: negative, ");

                 fprintf f "sta1: %i, end1: %i, " seg.sta1 seg.en1;
                 (match seg.dir1 with
                 | `Positive -> fprintf f "dir1: positive, "
                 | _ -> fprintf f "dir1: negative, ");

                 fprintf f "sta2: %i, end2: %i, " seg.sta2 seg.en2;
                 (match seg.dir2 with
                 | `Positive -> fprintf f "dir2: positive\n"
                 | _ -> fprintf f "dir2: negative\n");

                 fprintf f "%s\n" (Sequence.to_string seg.alied_seq1 Alphabet.nucleotides);
                 fprintf f "%s\n" (Sequence.to_string seg.alied_seq2 Alphabet.nucleotides);
                 
            ) (List.rev med.chrom_map);
    in 


    print_med (List.hd med_ls);
    close_out f


let create_map anc_med des_ref : (int * int * Tags.output) = 
    let str = string_of_int in  
    let seg_ls = List.map 
        (fun m -> 
             let a_ref_code, a_sta, a_en, a_dir  = 
                 (str anc_med.ref_code), (str m.sta), (str m.en), "+" 
             in 
             let d_ref_code, d_sta, d_en, d_dir = 
                 match des_ref = anc_med.ref_code1 with
                 | true ->
                       (str anc_med.ref_code1), (str m.sta1), (str m.en1), (Utl.get_dir m.dir1) 
                 | false ->
                       (str anc_med.ref_code2), (str m.sta2), (str m.en2), (Utl.get_dir m.dir2) 
             in 
             let attributes = [(Tags.GenomeMap.a_ref_code, a_ref_code);
                               (Tags.GenomeMap.a_start_seg, a_sta);
                               (Tags.GenomeMap.a_end_seg, a_en );
                               (Tags.GenomeMap.a_dir_seg, a_dir );
                               (Tags.GenomeMap.d_ref_code, d_ref_code);
                               (Tags.GenomeMap.d_start_seg, d_sta);
                               (Tags.GenomeMap.d_end_seg, d_en );
                               (Tags.GenomeMap.d_dir_seg, d_dir )
                              ] 
             in 
             let m : Tags.output = (Tags.GenomeMap.seg, attributes, `String "") in 
             `Single m
        ) anc_med.chrom_map 
    in 


    let chrom_map : Tags.output = 
        (Tags.GenomeMap.chrom, [], `Structured (`Set  seg_ls)) 
    in 

    match des_ref = anc_med.ref_code1 with
    | true -> anc_med.cost1, anc_med.recost1, chrom_map
    | false -> anc_med.cost2, anc_med.recost2, chrom_map




let create_single_map med : Tags.output = 
    let str = string_of_int in  
    let seg_ls = List.map 
        (fun m -> 
             let a_ref_code, a_sta, a_en, a_dir  = 
                 (str med.ref_code2), (str m.sta2), (str m.en2), "+" 
             in 
             let d_ref_code, d_sta, d_en, d_dir = 
                 (str med.ref_code1), (str m.sta1), (str m.en1), (Utl.get_dir m.dir1) 
             in 
             let attributes = [(Tags.GenomeMap.a_ref_code, a_ref_code);
                               (Tags.GenomeMap.a_start_seg, a_sta);
                               (Tags.GenomeMap.a_end_seg, a_en );
                               (Tags.GenomeMap.a_dir_seg, a_dir );
                               (Tags.GenomeMap.d_ref_code, d_ref_code);
                               (Tags.GenomeMap.d_start_seg, d_sta);
                               (Tags.GenomeMap.d_end_seg, d_en );
                               (Tags.GenomeMap.d_dir_seg, d_dir )
                              ] 
             in 
             let m : Tags.output = (Tags.GenomeMap.seg, attributes, `String "") in 
             `Single m
        ) med.chrom_map 
    in 


    let chrom_map : Tags.output = 
        (Tags.GenomeMap.chrom, [], `Structured (`Set  seg_ls)) 
    in 
    chrom_map








(** Create a global map between two chromsomes. Rearrangements are taken into account *)
let rec create_global_map (seq1 : Sequence.s) (seq2 : Sequence.s) cost_mat ali_pam =
    let pos_seed_ls, neg_seed_ls =
        Seed.determine_seed seq1 seq2 ali_pam `BothDir
    in 

    let pos_block_ls = Block.find_local_block pos_seed_ls ali_pam in
    let neg_block_ls = Block.find_local_block neg_seed_ls ali_pam in
(*    print_endline "Positive blocks";
    List.iter Block.print pos_block_ls;
    print_endline "Negative blocks";
    List.iter Block.print neg_block_ls;
*)    

    List.iter (fun block -> Block.invert block ali_pam.ChromPam.min_pos2 ali_pam.ChromPam.max_pos2) neg_block_ls;     
    let all_b_ls = pos_block_ls @ neg_block_ls in 
    let all_b_ls = List.filter 
        (fun b -> Block.max_len b >= ali_pam.ChromPam.sig_block_len) all_b_ls
    in 
    
    
    let sep_b_ls = Block.select_separated_block all_b_ls ali_pam in  
(*
    List.iter Block.print sep_b_ls;
    print_endline "End of sep blocks list"; 
*)

    Block.create_alied_block_ls sep_b_ls ali_pam seq1 seq2 cost_mat;

(*    let _ = List.fold_left (fun id b -> b.Block.id <- id; id + 1) 0 sep_b_ls in *)
        
    let block_pam = Block.blockPam_default in
    let con_sep_b_ls = Block.connect_consecutive_block 
        sep_b_ls block_pam seq1 seq2 cost_mat ali_pam 
    in
    
    let sig_map_ls, subseq1_ls, subseq2_ls = 
        Block.create_subseq_id `Alied con_sep_b_ls ali_pam in   


    List.iter (fun b ->
                   let alied_seq1 = deref b.Block.alied_seq1 in
                   let alied_seq2 = deref b.Block.alied_seq2 in
                   let cost = UtlPoy.cmp_ali_cost alied_seq1 alied_seq2 
                       b.Block.direction cost_mat in                    
                   b.Block.cost <- cost;     
              ) sig_map_ls;
                   
(*
    List.iter Block.print sig_map_ls;
    print_endline "End of con_sep blocks list";  
*)

    (List.rev sig_map_ls),  subseq1_ls, subseq2_ls
      
 
let check_chrom_map seq1 seq2 chrom_map =
    let len1 = Sequence.length seq1 in
    let len2 = Sequence.length seq2 in 
    let mark1_arr = Array.make len1 0 in 
    let mark2_arr = Array.make len2 0 in 
    List.iter (fun seg ->
                   (if seg.sta1 >= 0 && seg.en1 >= 0 then 
                       for p = seg.sta1 to seg.en1 do
                           mark1_arr.(p) <- mark1_arr.(p) + 1
                       done);

                   (if seg.sta2 >= 0 && seg.en2 >= 0 then 
                       for p = seg.sta2 to seg.en2 do
                           mark2_arr.(p) <- mark2_arr.(p) + 1
                       done)
              ) chrom_map;
 
    for p = 0 to len1 - 1 do
        if mark1_arr.(p) != 1 then begin
             fprintf stdout "The position %i in sequence1 is covered %i times" p mark1_arr.(p);
             failwith "Create median in ChromAli"
         end 
    done; 

    for p = 0 to len2 - 1 do
        if mark2_arr.(p) != 1 then begin
            fprintf stdout "The position %i in sequence2 is covered %i times" p mark2_arr.(p);
             failwith "Create median in ChromAli"
        end
    done

        
   

(** Create the median between two chromosomes which are divided into subseqs lists *)
let create_median subseq1_ls subseq2_ls (seq1, chrom1_id) (seq2, chrom2_id) global_map 
        ali_mat alied_gen_seq1 alied_gen_seq2 
        (order2_arr, total_cost, recost1, recost2) cost_mat ali_pam = 

    let approx = ali_pam.ChromPam.approx in


    let locus_indel_cost = ali_pam.ChromPam.locus_indel_cost in 
    
    let adder (submed_ls, nascent_len, chrom_map) ali_pos  = 
        let gen_code1 = alied_gen_seq1.(ali_pos) in 
        let gen_code2 = alied_gen_seq2.(ali_pos) in 

        match gen_code1, gen_code2 with
        | 0, 0 -> submed_ls, nascent_len, chrom_map
        | 0, _ -> 
              let subseq2 = List.find 
                  (fun sq -> sq.Subseq.id = gen_code2) subseq2_ls in
              let sta2 = subseq2.Subseq.sta in
              let en2 = subseq2.Subseq.en in 
 
              let len2 = en2 - sta2 + 1 in
              let subseq1 = UtlPoy.create_gap_seq len2 in 
              let subseq2 = Sequence.sub seq2 sta2 len2 in
              let submed, _ = UtlPoy.create_median_seq ~approx:approx subseq1 subseq2 cost_mat in

              let med_len = UtlPoy.cmp_num_not_gap submed in
              let sta, en, nascent_len = match med_len with
              | 0 -> -1, -1, nascent_len
              | _ -> nascent_len + 1, nascent_len + med_len, nascent_len + med_len
              in  

              let map = {sta=sta; en = en; alied_med = submed;
                         cost = UtlPoy.cmp_gap_cost locus_indel_cost subseq2;
                         sta1 = -1; en1 = -1; alied_seq1 = subseq1; dir1 = `Positive;                         
                         sta2 = sta2; en2 = en2; alied_seq2 = subseq2; dir2 = `Positive
                        }
              in 
              List.append submed_ls [submed], nascent_len, List.append chrom_map [map]

        | _, 0 ->
              let subseq1 = List.find 
                  (fun sq -> sq.Subseq.id = gen_code1) subseq1_ls in

              let sta1 = subseq1.Subseq.sta in
              let en1 = subseq1.Subseq.en in 
              let len1 = en1 - sta1 + 1 in

              let subseq1 = Sequence.sub seq1 sta1 len1 in
              let subseq2 = UtlPoy.create_gap_seq len1 in 
              let submed, _ = UtlPoy.create_median_seq ~approx:approx subseq1 subseq2 cost_mat in

              let med_len = UtlPoy.cmp_num_not_gap submed in 
              let sta, en, new_nascent_len = match med_len with
              | 0 -> -1, -1, nascent_len
              | _ -> nascent_len + 1, nascent_len + med_len, nascent_len + med_len
              in  
              let map = {sta=sta; en = en; alied_med = submed;
                         cost = UtlPoy.cmp_gap_cost locus_indel_cost subseq1;
                         sta1 = sta1; en1 = en1; alied_seq1 = subseq1; dir1 = `Positive;                         
                         sta2 = -1; en2 = -1; alied_seq2 = subseq2; dir2 = `Positive;
                        }
              in 
              List.append submed_ls [submed], new_nascent_len, List.append chrom_map [map]
                  


        | _, _ -> 
              let block_opt = Block.find_block global_map gen_code1 gen_code2 in 
              match block_opt with
              | Some b -> 
                    let submed, cost = Block.create_median ~approx:approx b cost_mat in
                    let med_len = UtlPoy.cmp_num_not_gap submed in 

                    let sta, en, nascent_len = match med_len with
                    | 0 -> -1, -1, nascent_len
                    | _ -> nascent_len + 1, nascent_len + med_len, nascent_len + med_len
                    in  
                    
                    let map = {sta=sta; en = en; alied_med = submed;
                               cost = cost;
                               sta1 = b.Block.sta1; en1 = b.Block.en1; 
                               alied_seq1 = deref b.Block.alied_seq1; dir1 = `Positive;                         

                               sta2 = b.Block.sta2; en2 = b.Block.en2; 
                               alied_seq2 = deref b.Block.alied_seq2; dir2 = b.Block.direction;
                              }
                    in 
                    List.append submed_ls [submed], nascent_len, List.append chrom_map [map]


              | None -> 
                    let alied_seq1, alied_seq2 = ali_mat.(gen_code1).(gen_code2) in 
                    
                    let subseq1 = List.find  
                        (fun sq -> sq.Subseq.id = gen_code1) subseq1_ls 
                    in
                    let subseq2 = List.find  
                        (fun sq -> sq.Subseq.id = gen_code2) subseq2_ls 
                    in

                    let submed, cost = UtlPoy.create_median_seq 
                        ~approx:approx alied_seq1 alied_seq2 cost_mat 
                    in 


                    let med_len = UtlPoy.cmp_num_not_gap submed in 
                    let sta, en, nascent_len = match med_len with
                    | 0 -> -1, -1, nascent_len
                    | _ -> nascent_len + 1, nascent_len + med_len, nascent_len + med_len 
                    in  

                    let map = {sta=sta; en = en; alied_med = submed;
                               cost = cost;
                               sta1 = subseq1.Subseq.sta; en1 = subseq1.Subseq.en; 
                               alied_seq1 = alied_seq1; dir1 = `Positive;                         

                               sta2 = subseq2.Subseq.sta; en2 = subseq2.Subseq.en; 
                               alied_seq2 = alied_seq2; dir2 = `Positive;                         
                              }
                    in 
                    List.append submed_ls [submed], nascent_len, List.append chrom_map [map]
    in
    


    let create_ali_order order_arr alied_gen_seq =
        let gen_gap_code = 0 in
        let ali_len = Array.length alied_gen_seq  in 
        let rec collect pos_ls pos = 
            let new_pos_ls = pos::pos_ls in
            if (pos + 1 < ali_len) && (alied_gen_seq.(pos + 1) = gen_gap_code) then 
                collect new_pos_ls (pos + 1)
            else new_pos_ls
    
        in
    
        let rec travel pos pos_ls = 
            match pos = Array.length order_arr with
            | true -> List.rev pos_ls
            | false ->
                  let order = order_arr.(pos) in
                  let index = Utl.find_index alied_gen_seq (abs order) compare in 
                  let new_pos_ls = collect pos_ls index in
                  travel (pos + 1) new_pos_ls
        in
    
        let sta_ls = 
            match alied_gen_seq.(0)= gen_gap_code with
            | true -> collect [] 0
            | false -> []
        in
        let ali_order_ls = travel 0 sta_ls in     
        ali_order_ls
    in

    let ali_order_ls : int list= create_ali_order order2_arr alied_gen_seq2 in 
    
    let submed_ls, med_len, chrom_map = List.fold_left adder ([], -1, [])  ali_order_ls  in
    if (List.length chrom_map = 0) then begin
        failwith "Chrom_map length is Zero";
    end;
    
    let seq = UtlPoy.delete_gap (UtlPoy.concat submed_ls) in 
    let ref_code = Utl.get_new_chrom_ref_code() in 

    (if Sequence.length seq = 0 then  failwith "Sequence length is Zero");
(*    check_chrom_map seq1 seq2 chrom_map; *)
    {seq = seq; ref_code = ref_code;
     ref_code1 = chrom1_id;
     ref_code2 = chrom2_id;
     chrom_map = chrom_map;
     cost1 = total_cost - recost2;
     cost2 = total_cost - recost1;
     recost1 = recost1;
     recost2 = recost2;
    }




(** Compute the cost between two sequences with rearrangement operations *)
let cmp_simple_cost med1 med2 cost_mat ali_pam = 
    if debug = true then begin
        let len1 = Sequence.length med1.seq in 
        let len2 = Sequence.length med2.seq in 
        fprintf stdout "Cmp_cost with lens %i %i: " len1 len2;
        flush stdout;
        let seqfile = open_out "seq12_cost" in 
        fprintf seqfile ">seq1\n";  
        Sequence.print seqfile med1.seq Alphabet.nucleotides;  
        fprintf seqfile "\n";  

        fprintf seqfile ">seq2\n";  
        Sequence.print seqfile med2.seq Alphabet.nucleotides;  
        close_out seqfile;  
        print_endline "End of printing seq12";  
    end;   


    let seq1 = med1.seq and seq2 = med2.seq in    
    let len1 = Sequence.length seq1 in
    let len2 = Sequence.length seq2 in 

    if (len1 < 2) or (len2 < 2) then 0, 0
    else begin
        let ali_pam = {ali_pam with  
                           ChromPam.min_pos1 = 0; 
                           ChromPam.max_pos1 = len1 - 1; 
                           ChromPam.min_pos2 = 0; 
                           ChromPam.max_pos2 = len2 - 1; 
                      } 
        in  
        
        let global_map, _, _ = create_global_map seq1 seq2 cost_mat ali_pam in  
        let _, _, _, _, _, _, total_cost, (recost1, recost2) =
            AliMap.create_general_ali `Chromosome global_map seq1 seq2 cost_mat ali_pam
        in     



        total_cost, (recost1 + recost2)
    end  



let cmp_cost med1 med2 cost_mat chrom_pams state = 
    if debug = true then begin
        let len1 = Sequence.length med1.seq in 
        let len2 = Sequence.length med2.seq in 
        fprintf stdout "Cmp_cost with lens %i %i: " len1 len2;
        flush stdout;
        let seqfile = open_out "seq12_cost" in 
        fprintf seqfile ">seq1\n";  
        Sequence.print seqfile med1.seq Alphabet.nucleotides;  
        fprintf seqfile "\n";  

        fprintf seqfile ">seq2\n";  
        Sequence.print seqfile med2.seq Alphabet.nucleotides;  
        close_out seqfile;  
        print_endline "End of printing seq12";  
    end;   
    let ali_pam = ChromPam.get_chrom_pam chrom_pams in 
    let ali_pam = match state with 
    | `Genome -> {ali_pam with ChromPam.negative = false}
    | `Chromosome -> ali_pam
    in 

    let cost, recost =  
        match ali_pam.ChromPam.symmetric with
        | true ->
              let cost12 = cmp_simple_cost med1 med2 cost_mat ali_pam in 
              let cost21 = cmp_simple_cost med2 med1 cost_mat ali_pam in 
              min cost12 cost21
        | false -> 
              if Sequence.compare med1.seq med2.seq < 0 then 
                  cmp_simple_cost med1 med2 cost_mat ali_pam  
              else cmp_simple_cost med2 med1 cost_mat ali_pam  

    in 
    cost, recost

        

let find_simple_med2_ls (med1 : med_t) (med2 : med_t) cost_mat ali_pam =
    if debug = true then begin
        let len1 = Sequence.length med1.seq in 
        let len2 = Sequence.length med2.seq in 
        fprintf stdout "Find median list with lens: %i %i" len1 len2; print_newline ();
        let seqfile = open_out "seq12" in
        fprintf seqfile ">seq1\n";
        Sequence.print seqfile med1.seq Alphabet.nucleotides;
        fprintf seqfile "\n";

        fprintf seqfile ">seq2\n";

        Sequence.print seqfile med2.seq Alphabet.nucleotides;
        close_out seqfile;
    end;

    let seq1 = med1.seq and seq2 = med2.seq in    
    let len1 = Sequence.length seq1 in
    let len2 = Sequence.length seq2 in 

    let ali_pam = {ali_pam with 
                       ChromPam.min_pos1 = 0;
                       ChromPam.max_pos1 = len1 - 1;
                       ChromPam.min_pos2 = 0;
                       ChromPam.max_pos2 = len2 - 1;
                  }
    in 
    
    if (len1 < 2) then begin
        let med = {seq = Sequence.clone seq2; 
                   ref_code = Utl.get_new_chrom_ref_code (); 
                   ref_code1 = -1;
                   ref_code2 = med2.ref_code; 
                   chrom_map = [];
                   cost1 = 0; recost1 = 0;
                   cost2 = 0; recost2 = 0;
                  } 
        in 
        0, 0, [med]
    end  else if len2 < 2 then begin 
        let med = {seq = Sequence.clone seq1; 
                   ref_code = Utl.get_new_chrom_ref_code (); 
                   ref_code1 = med1.ref_code;
                   ref_code2 = -1;
                   chrom_map = [];
                   cost1 = 0; recost1 = 0;
                   cost2 = 0; recost2 = 0;                  
                  } in 
        0, 0, [med]
    end else begin
        let global_map, _, _ = create_global_map seq1 seq2 cost_mat ali_pam in 
(*       List.iter Block.print global_map; *)

        let subseq1_ls, subseq2_ls, global_map, ali_mat, alied_gen_seq1,
            alied_gen_seq2, total_cost, (_, recost)  = 
            AliMap.create_general_ali `Chromosome global_map seq1 seq2 cost_mat ali_pam 
        in


        let re_gen_seq2 = Utl.filterArray (fun code2 -> code2 > 0) alied_gen_seq2 in 
        let gen_seq2 = UtlGrappa.get_ordered_permutation re_gen_seq2 in 

        let all_order_ls = 
            if ali_pam.ChromPam.approx = `Second then [(gen_seq2, recost, 0)]
            else 
            if (Utl.equalArr gen_seq2 re_gen_seq2 compare) ||
                (ali_pam.ChromPam.keep_median = 1) || 
                (ali_pam.ChromPam.approx = `First) then [(re_gen_seq2, 0, recost)]
            else [(re_gen_seq2, 0, recost); (gen_seq2, recost, 0)]
        in 


        let med_ls = List.fold_right
            (fun (order_arr, recost1, recost2) med_ls ->
                 let med = 
                     create_median subseq1_ls subseq2_ls 
                         (seq1, med1.ref_code) (seq2, med2.ref_code) global_map
                         ali_mat alied_gen_seq1 alied_gen_seq2 
                         (order_arr, total_cost, recost1, recost2) cost_mat ali_pam
                 in
                 med::med_ls
            ) all_order_ls []
        in
        total_cost, recost, med_ls
    end 




let find_med2_ls (med1 : med_t) (med2 : med_t) cost_mat user_chrom_pam = 
    let ali_pam = ChromPam.get_chrom_pam user_chrom_pam in 
    match ali_pam.ChromPam.symmetric with
    | true ->
          let cost12, recost12, med12_ls = find_simple_med2_ls med1 med2 cost_mat ali_pam in 
          
          let ali_pam = 
              if ali_pam.ChromPam.approx = `First then {ali_pam with ChromPam.approx = `Second}
              else ali_pam
          in 
          let cost21, recost21, med21_ls = find_simple_med2_ls med2 med1 cost_mat ali_pam in 
          if cost12 <= cost21 then cost12, recost12, med12_ls
          else begin 
              let med12_ls = List.map swap_med med21_ls in 
              cost21, recost21, med12_ls
          end 

    | false ->
          let med1, med2, ali_pam, swaped = 
              match Sequence.compare med1.seq med2.seq < 0 with 
              | true ->  med1, med2, ali_pam, false
              | false -> 
                    let ali_pam = 
                        if ali_pam.ChromPam.approx = `First then {ali_pam with ChromPam.approx = `Second}
                        else ali_pam
                    in 
                    med2, med1, ali_pam, true
          in 

          let cost, recost, med_ls = find_simple_med2_ls med1 med2 cost_mat ali_pam in
          let med_ls = 
              match swaped with
              | false -> med_ls
              | true -> List.map swap_med med_ls 
          in 
          cost, recost, med_ls





let find_approx_med2 (med1 : med_t) (med2 : med_t) (med12 : med_t) =
    let new_med12 = clone_med med12 in 
    let ref_code1 = med1.ref_code in 
    let ref_code2 = med2.ref_code in 
    let ref_code = Utl.get_new_chrom_ref_code () in 
    {new_med12 with ref_code = ref_code; ref_code1 = ref_code1; ref_code2 = ref_code2}



let test () =
    let file_name = Sys.argv.(1) in
    let seq_file = open_in file_name in
    let rev_ls = Parser.Fasta.of_channel Parser.Nucleic_Acids seq_file in
    let rev_seq_ls : Sequence.s list = 
        List.map (fun (a, b) -> List.hd (List.hd (List.hd a))) rev_ls in


    let seq_ls = List.rev rev_seq_ls in    
    let len = List.length seq_ls in 
    for id0 = 0 to len - 2 do
        for id1 = id0 + 1 to len - 1 do
            let seq1 = List.nth seq_ls id0 in 
            let len1 = Sequence.length seq1 in  
            let seq1 = Sequence.sub seq1 1 (len1 - 1) in  

            let seq2 = List.nth seq_ls id1 in 
            let len2 = Sequence.length seq2 in  
            let seq2 = Sequence.sub seq2 1 (len2 - 1) in  

            let cost_mat = Cost_matrix.Two_D.default in  
            let chrom_pam = Data.dyna_pam_default in   
            let chrom_pam = {chrom_pam with Data.approx = Some true} in 

            let med1 = {seq = seq1; ref_code = -1; 
                        ref_code1 = -1; ref_code2 = -1; chrom_map = [];
                        cost1 = 0; recost1 = 0;
                        cost2 = 0; recost2 = 0;
                       } in  
            let med2 = {seq = seq2; ref_code = -1; 
                        ref_code1 = -1; ref_code2 = -1; chrom_map = [];
                        cost1 = 0; recost1 = 0; 
                        cost2 = 0; recost2 = 0;
                       } in  

            
            let total_cost, _, med_ls = find_med2_ls med1 med2 cost_mat chrom_pam in    
            print_median med_ls (file_name ^ ".ali");
            UtlPoy.printDNA (List.hd med_ls).seq;
            fprintf stdout "Total cost: %i \n End of testing in ChromAli!!!" total_cost;
            print_newline ();
        done
    done;
    failwith "Finish testing chromAli"
        


let copy_chrom_map s d = 
    {d with ref_code = s.ref_code; 
         ref_code1 = s.ref_code1;
         ref_code2 = s.ref_code2;
         chrom_map = s.chrom_map}
    

let to_single single_parent child_ref c2 pam = 
    if (single_parent.ref_code1 = -1) && (single_parent.ref_code2 = -1) then single_parent.seq
    else begin
        let gap = Cost_matrix.Two_D.gap c2 in    
        let map = List.map 
            (fun seg -> 
                 let alied_single_seq, alied_child_seq = 
                     match child_ref = single_parent.ref_code1 with
                     | true ->
                           let single, _ = UtlPoy.closest_alied_seq
                               seg.alied_med seg.alied_seq1 c2
                           in 
                           single, seg.alied_seq1
                     | false ->
                           let single, _  = 
                               if seg.dir2 = `Positive then 
                                   UtlPoy.closest_alied_seq seg.alied_med seg.alied_seq2 c2 
                               else  begin
                                   let rev_alied_seq2 = Sequence.reverse seg.alied_seq2 in 
                                   let rev_single_alied_seq2, _ = UtlPoy.closest_alied_seq
                                       seg.alied_med rev_alied_seq2 c2
                                   in 
                                   (Sequence.reverse rev_single_alied_seq2), 0
                               end 
                           in 
                           single, seg.alied_seq2 
                 in  
                 let sta = 
                     match child_ref = single_parent.ref_code1 with 
                     | true -> seg.sta1
                     | false -> seg.sta2
                 in 
                 let ungap_alied_med = Sequence.fold_righti  
                     (fun ungap_alied_med p code ->
                          match code = gap with
                          | false -> code::ungap_alied_med
                          | true ->                            
                                if Sequence.get alied_child_seq p = gap then  ungap_alied_med
                                else code::ungap_alied_med
                     ) [] alied_single_seq
                 in                       
                 
                 let ungap_alied_med = UtlPoy.of_array (Array.of_list  ungap_alied_med) in
                 ungap_alied_med, sta
            ) single_parent.chrom_map 
        in  
        
        let sorted_map = List.sort 
            (fun seg1 seg2 ->
                 let _,  sta1 = seg1 in 
                 let _,  sta2 = seg2 in
                 sta1 - sta2
            ) map 
        in     
        
        let seq_ls = List.fold_right 
            (fun seg seq_ls -> 
                 let seq, sta = seg in 
                 if sta < 0 then seq_ls
                 else seq::seq_ls
            ) sorted_map []
        in  
        let single_seq = UtlPoy.concat seq_ls in   
        single_seq
    end 



let to_single_root root other_code c2 =
    if (root.ref_code1 = -1) && (root.ref_code2 = -1) then root.seq
    else begin
        let gap = Cost_matrix.Two_D.gap c2 in 
        let map = List.map 
            (fun seg ->  
                 let child_alied_seq = seg.alied_seq1 in 
                 let other_alied_seq = match seg.dir2 = `Positive with 
                 | true -> seg.alied_seq2
                 | false -> Sequence.reverse seg.alied_seq2 
                 in 
                 let single_seg, _ = UtlPoy.closest_alied_seq 
                     other_alied_seq child_alied_seq c2 in 
                 let ungap_alied_med = Sequence.fold_righti  
                     (fun ungap_alied_med p code ->
                          match code = gap with
                          | false -> code::ungap_alied_med
                          | true ->
                                if Sequence.get seg.alied_med p = gap then ungap_alied_med
                                else code::ungap_alied_med
                     ) [] single_seg
                 in                       
                 let ungap_alied_med = UtlPoy.of_array (Array.of_list ungap_alied_med) in         
                 ungap_alied_med
            ) root.chrom_map 
        in  
        let single_seq = UtlPoy.concat map in 
        single_seq
    end 

let change_to_single med single_seq c2 = 
    (if Sequence.length single_seq != Sequence.length med.seq then begin
        fprintf stdout "single_len: %i, med_len:%i\n" (Sequence.length single_seq)
            (Sequence.length med.seq); flush stdout;
        failwith "XXXXXXXXXXXXXXXXXXXXXX at change_to_single @ chromAli.ml"
    end);


    let gap = Cost_matrix.Two_D.gap c2 in
    let single_pos = ref 0 in 
    let new_map = List.map 
        (fun seg ->
             let single_alied_med = UtlPoy.map 
                 (fun code ->
                      let single_code = 
                          if code = gap then gap
                          else begin
                              let single_code = Sequence.get single_seq !single_pos in
(*                              (if (single_code land code = 0) then begin
                                   fprintf stdout "Code: %i, single_code: %i" code single_code;                               
                                   failwith "The code does not include the single_code";
                               end);*)
                              incr single_pos;
                              single_code
                          end 
                      in 
                      single_code
                 ) seg.alied_med
             in 
             let new_seg = {seg with alied_med = single_alied_med} in 
             new_seg
        ) med.chrom_map
    in 
    let gapless_single_seq = UtlPoy.delete_gap single_seq in 
    {med with seq = gapless_single_seq;
         chrom_map = new_map}



    
