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

let () = SadmanOutput.register "ChromAli" "$Revision: 2869 $"

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
let  no_single_error = true

type direction_t = ChromPam.direction_t

(** A pair of segements on both chromosomes and the median sequence *)
type seg_t = {
    sta : int; (** start position of this segment on the median sequence *)
    en : int; (** end position of this segment on the median sequence *)
    cost : int; (** the cost of aligning the segment *)
    alied_med : Sequence.s; (** the median sequence of the segment *)

    sta1 : int; (** start position of the segment on the first chromosome *)
    en1 : int; (** end position of the segment on the first chromosome *)
    alied_seq1 : Sequence.s; (** aligned sequence of this segment on the first chromosome *)
    dir1 : direction_t; (** the orientation of this segment on the first chromosome *) 

    sta2 : int; (** start position of the segment on the second chromosome *)
    en2 : int;  (** end position of the segment on the second chromosome *)
    alied_seq2 : Sequence.s; (** alied_seq1 <-> reversed alied_seq2 
                                if dir2 is Negative *)
    dir2 : direction_t; (** the orientation of this segment on the first chromosome *) 
}

(** the median data structure of two chromosomes *)
type med_t = {
    seq : Sequence.s; (** the median sequence *)
    ref_code : int; (** reference code of this median *)
    ref_code1 : int;    (** Child's code *)    
    ref_code2 : int;  (** Child's code *)
    cost1 : int;  
    recost1 : int;
    cost2 : int;
    recost2 : int;
    chrom_map : seg_t list;    
}


(** [create_med seq] return a new median from a seq *)
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


let get_recost user_pams = 
    match user_pams.Data.re_meth with
    | None -> failwith "The rearrangement cost is not specified"
    | Some re_meth ->
        match re_meth with
            | `Locus_Breakpoint c -> c
            | `Locus_Inversion c -> c



(** [clone_seg s] return a fresh clone of segment [s] *)
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

(** [clone_med m] returns a fresh clone of median [m] *)
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

(** [swap_seg s] swaps the first segment and 
* the second segment of this segment pair *)
let swap_seg s = 
    {s with sta1 = s.sta2; en1 = s.en2; alied_seq1 = s.alied_seq2; dir1 = s.dir2;
         sta2 = s.sta1; en2 = s.en1; alied_seq2 = s.alied_seq1; dir2 = s.dir1}
            
let swap_chrom_map m = 
    List.map swap_seg m

(** [swap_med m] swaps the first chromosome 
* and the second chromosome of this median *)
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


(** [create_map anc_med des_ref] creates a map from the 
* ancestor sequence  [anc_med] to descendant sequence whose ref_code is [des_ref]*)
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


(** [create_global_map seq1 seq2 cost_mat ali_pam] creates the 
*  global map between two chromosomes [seq1] and [seq2] 
* Rearrangements are allowed in the global map *)
let rec create_global_map (seq1 : Sequence.s) (seq2 : Sequence.s) cost_mat ali_pam =
    let pos_seed_ls, neg_seed_ls =
        Seed.determine_seed seq1 seq2 ali_pam `BothDir
    in 

    let pos_block_ls = Block.find_local_block pos_seed_ls ali_pam in
    let neg_block_ls = Block.find_local_block neg_seed_ls ali_pam in


    List.iter (fun block -> Block.invert block ali_pam.ChromPam.min_pos2 ali_pam.ChromPam.max_pos2) neg_block_ls;     
    let all_b_ls = pos_block_ls @ neg_block_ls in 
    let all_b_ls = List.filter 
        (fun b -> Block.max_len b >= ali_pam.ChromPam.sig_block_len) all_b_ls
    in 
    
    
    let sep_b_ls = Block.select_separated_block all_b_ls ali_pam in  

    Block.create_alied_block_ls sep_b_ls ali_pam seq1 seq2 cost_mat;
        
    let block_pam = Block.blockPam_default in
    let con_sep_b_ls = Block.connect_consecutive_block 
        sep_b_ls block_pam seq1 seq2 cost_mat ali_pam 
    in
    
    let sig_map_ls, subseq1_ls, subseq2_ls = 
        Block.create_subseq_id `Alied con_sep_b_ls ali_pam in   


    List.iter (fun b ->
                   let alied_seq1 = deref b.Block.alied_seq1 in
                   let alied_seq2 = deref b.Block.alied_seq2 in
                   let cost = Sequence.cmp_ali_cost alied_seq1 alied_seq2 
                       b.Block.direction cost_mat in                    
                   b.Block.cost <- cost;     
              ) sig_map_ls;
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

        
   

(** [create_median subseq1_ls subseq2_ls gen_gap_code (seq1, chrom1_id) (seq2, chrom2_id) global_map 
*                  ali_mat alied_gen_seq1 alied_gen_seq2 (order2_arr, total_cost, recost1, recost2) 
*                   cost_mat ali_pam] creates the median between 
* two chromosome [seq1] and [seq2] given a [global_map] *)
let create_median subseq1_ls subseq2_ls gen_gap_code (seq1, chrom1_id) (seq2, chrom2_id) global_map 
        ali_mat alied_gen_seq1 alied_gen_seq2 
        (order2_arr, total_cost, recost1, recost2) cost_mat ali_pam = 

    let approx = ali_pam.ChromPam.approx in

    let locus_indel_cost = ali_pam.ChromPam.locus_indel_cost in 
    
    let adder (submed_ls, nascent_len, chrom_map) ali_pos  = 
        let gen_code1 = alied_gen_seq1.(ali_pos) in 
        let gen_code2 = alied_gen_seq2.(ali_pos) in 

        match gen_code1 = gen_gap_code, gen_code2 = gen_gap_code with
        | true, true -> submed_ls, nascent_len, chrom_map
        | true, _ -> 
              let subseq2 = List.find 
                  (fun sq -> sq.Subseq.id = gen_code2) subseq2_ls in
              let sta2 = subseq2.Subseq.sta in
              let en2 = subseq2.Subseq.en in 
 
              let len2 = en2 - sta2 + 1 in
              let subseq1 = Sequence.create_gap_seq len2 in 
              let subseq2 = Sequence.sub seq2 sta2 len2 in
              let submed, _ = Sequence.create_median_seq ~approx:approx subseq1 subseq2 cost_mat in

              let med_len = Sequence.cmp_num_not_gap submed in
              let sta, en, nascent_len = match med_len with
              | 0 -> -1, -1, nascent_len
              | _ -> nascent_len + 1, nascent_len + med_len, nascent_len + med_len
              in  

              let map = {sta=sta; en = en; alied_med = submed;
                         cost = Sequence.cmp_gap_cost locus_indel_cost subseq2;
                         sta1 = -1; en1 = -1; alied_seq1 = subseq1; dir1 = `Positive;                         
                         sta2 = sta2; en2 = en2; alied_seq2 = subseq2; dir2 = `Positive
                        }
              in 
              List.append submed_ls [submed], nascent_len, List.append chrom_map [map]

        | _, true ->
              let subseq1 = List.find 
                  (fun sq -> sq.Subseq.id = gen_code1) subseq1_ls in

              let sta1 = subseq1.Subseq.sta in
              let en1 = subseq1.Subseq.en in 
              let len1 = en1 - sta1 + 1 in

              let subseq1 = Sequence.sub seq1 sta1 len1 in
              let subseq2 = Sequence.create_gap_seq len1 in 
              let submed, _ = Sequence.create_median_seq ~approx:approx subseq1 subseq2 cost_mat in

              let med_len = Sequence.cmp_num_not_gap submed in 
              let sta, en, new_nascent_len = match med_len with
              | 0 -> -1, -1, nascent_len
              | _ -> nascent_len + 1, nascent_len + med_len, nascent_len + med_len
              in  
              let map = {sta=sta; en = en; alied_med = submed;
                         cost = Sequence.cmp_gap_cost locus_indel_cost subseq1;
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
                    let med_len = Sequence.cmp_num_not_gap submed in 

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

                    let submed, cost = Sequence.create_median_seq 
                        ~approx:approx alied_seq1 alied_seq2 cost_mat 
                    in 


                    let med_len = Sequence.cmp_num_not_gap submed in 
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
    
    let seq = Sequence.delete_gap (Sequence.concat submed_ls) in 
    let ref_code = Utl.get_new_chrom_ref_code() in 

    (if Sequence.length seq = 0 then  failwith "Sequence length is Zero");

    {seq = seq; ref_code = ref_code;
     ref_code1 = chrom1_id;
     ref_code2 = chrom2_id;
     chrom_map = chrom_map;
     cost1 = total_cost - recost2;
     cost2 = total_cost - recost1;
     recost1 = recost1;
     recost2 = recost2;
    }




(** [cmp_simple_cost med1 med2 cost_mat ali_pam] computes
* the cost between sequence [med1] and sequence [med2]. 
* Rearrangement operations are allowed *)
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

        let _, _, _, _, _, _, _, total_cost, (recost1, recost2) =
            AliMap.create_general_ali `Chromosome global_map seq1 seq2 cost_mat ali_pam
        in     


        total_cost, (recost1 + recost2)
    end  



(** [cmp_cost med1 med2 cost_mat chrom_pams state] computes
* the cost which is the min of cost between [med1] and [med2],  
* and cost between [med2] and [med1]. Rearrangement operations are allowed *)
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

        
(** [find_simple_med2_ls med1 med2d cost_mat ali_pam] finds 
* the median list between chromosome [med1] and chromosome [med2]. 
* Rearrangements are allowed *)
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


        let subseq1_ls, subseq2_ls, gen_gap_code, global_map, ali_mat, alied_gen_seq1,
            alied_gen_seq2, total_cost, (_, recost)  = 
            AliMap.create_general_ali `Chromosome global_map seq1 seq2 cost_mat ali_pam 
        in

        let re_gen_seq2 = Utl.filterArray (fun code2 -> code2 != gen_gap_code) alied_gen_seq2 in 
        let gen_seq2 = UtlGrappa.get_ordered_permutation re_gen_seq2 in 

        let all_order_ls = 
            if ali_pam.ChromPam.approx = `Second then [(gen_seq2, recost, 0)]
            else 
            if (Utl.isEqualArr gen_seq2 re_gen_seq2 compare) ||
                (ali_pam.ChromPam.keep_median = 1) || 
                (ali_pam.ChromPam.approx = `First) then [(re_gen_seq2, 0, recost)]
            else [(re_gen_seq2, 0, recost); (gen_seq2, recost, 0)]
        in 


        let med_ls = List.fold_right
            (fun (order_arr, recost1, recost2) med_ls ->
                 let med = 
                     create_median subseq1_ls subseq2_ls gen_gap_code
                         (seq1, med1.ref_code) (seq2, med2.ref_code) global_map
                         ali_mat alied_gen_seq1 alied_gen_seq2 
                         (order_arr, total_cost, recost1, recost2) cost_mat ali_pam
                 in
                 med::med_ls
            ) all_order_ls []
        in
(*
       List.iter Block.print global_map; 
        fprintf stdout "Total_cost, recost: %i %i\n" total_cost recost; flush stdout;
*)  
      total_cost, recost, med_ls
    end 



(** [find_med2_ls med1 med2 cost_mat user_chrom_pam] find the
* median list whose cost is minimum cost between [med1] and [med2]
* and between [med2] and [med1] *) 
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




(** [find_approx_med2 med1 med2 med12] returns the median 
* between [med1] and [med2] based on all information 
* from [med12] which is a computed median of [med1] and [med2] *)
let find_approx_med2 (med1 : med_t) (med2 : med_t) (med12 : med_t) =
    let new_med12 = clone_med med12 in 
    let ref_code1 = med1.ref_code in 
    let ref_code2 = med2.ref_code in 
    let ref_code = Utl.get_new_chrom_ref_code () in 
    {new_med12 with ref_code = ref_code; ref_code1 = ref_code1; ref_code2 = ref_code2}


(** [find_med3 ch1 ch2 ch3 mine c2 c3 pam] returns
* the median of [ch1], [ch2] and [ch3] where [mine]
* is the current median of [ch1], [ch2] and [ch3] *) 
let find_med3 ch1 ch2 ch3 mine c2 c3 pam = 
    let ali_pam = ChromPam.get_chrom_pam pam in 
    let _, _, med1m_ls = find_med2_ls ch1 mine c2 pam in
    let _, _, med2m_ls = find_med2_ls ch2 mine c2 pam in
    let _, _, med3m_ls = find_med2_ls ch3 mine c2 pam in

    let med1m = List.hd med1m_ls in 
    let med2m = List.hd med2m_ls in 
    let med3m = List.hd med3m_ls in 
    
    let gap = Cost_matrix.Two_D.gap c2 in 
    let create_pos chrom_map =
        let sorted_chrom_map = List.sort 
            (fun seg1 seg2 -> seg1.sta2 - seg2.sta2) chrom_map
        in 
        let pos_ls = List.fold_left 
            (fun pos_ls seg -> 
                 let alied_arr1 = Sequence.to_array seg.alied_seq1 in 
                 let alied_arr2 = Sequence.to_array seg.alied_seq2 in 
                 let p1 = ref 0 in 
                 let num_base1 = ref 0 in 
                 let pos_ls = Array.fold_left 
                     (fun (pos_ls : int list) base2 -> 
                          let base1 = alied_arr1.(!p1) in 
                          p1:= !p1 + 1;
                          (if base1 != gap then 
                              num_base1 := !num_base1 + 1);
                          if base2 = gap then pos_ls
                          else 
                              if (seg.dir2 = `Negative) || (seg.sta1 = -1) then (-2)::pos_ls
                              else 
                                  if base1 = gap then (-1)::pos_ls
                                  else (seg.sta1 + !num_base1 - 1)::pos_ls
                     ) pos_ls alied_arr2
                 in 
                 pos_ls
            ) [] sorted_chrom_map 
        in 
        Array.of_list (List.rev pos_ls)
    in  

    let pos1_arr = create_pos med1m.chrom_map in 
    let pos2_arr = create_pos med2m.chrom_map in 
    let pos3_arr = create_pos med3m.chrom_map in 

    let max_3d_len = min ali_pam.ChromPam.detected_3d_len ali_pam.ChromPam.max_3d_len in 
    let mine_len = Sequence.length mine.seq in

    let rec detect_change new_med f_p =
        if f_p >= mine_len then new_med
        else begin
            let rec find_first f_p = 
                if f_p = mine_len then -1, -1, -1, -1
                else begin 
                    let f_p1 = pos1_arr.(f_p) 
                    and f_p2 = pos2_arr.(f_p)
                    and f_p3 = pos3_arr.(f_p)                         
                    in  
                    if f_p1 >=0 && f_p2 >=0 && f_p3 >=0 then f_p, f_p1, f_p2, f_p3
                    else find_first (f_p + 1)
                end
            in 

            let f_p, f_p1, f_p2, f_p3 = find_first f_p in 
            let rec find_last l_p = 
                if l_p = f_p then -1, -1, -1, -1
                else begin
                    let l_p1 = pos1_arr.(l_p) 
                    and l_p2 = pos2_arr.(l_p)
                    and l_p3 = pos3_arr.(l_p)                         
                    in  
                    if (l_p1 >= 0 ) && (l_p1 > f_p1) && (l_p1 - f_p1 <= max_3d_len)
                        && (l_p2 >= 0 ) && (l_p2 > f_p2) && (l_p2 - f_p2 <= max_3d_len)
                        && (l_p3 >= 0 ) && (l_p3 > f_p3) && (l_p3 - f_p3 <= max_3d_len) 
                    then l_p, l_p1, l_p2, l_p3
                    else find_last (l_p - 1)
               end 
            in
            
            if f_p = -1 then new_med 
            else begin
                let max_l_p = min (f_p + max_3d_len) (mine_len - 1)  in
                let l_p, l_p1, l_p2, l_p3 = find_last max_l_p in 
                if l_p = -1 then 
                    detect_change new_med (max_l_p + 1)
                else begin
(*                    fprintf stdout "(%i, %i), (%i, %i), (%i, %i), (%i, %i)\n"
                        f_p l_p f_p1 l_p1 f_p2 l_p2 f_p3 l_p3; flush stdout;
*)
                    let sub_seq1 = Sequence.sub ch1.seq f_p1 (l_p1 - f_p1 + 1) in
                    let sub_seq2 = Sequence.sub ch2.seq f_p2 (l_p2 - f_p2 + 1) in
                    let sub_seq3 = Sequence.sub ch3.seq f_p3(l_p3 - f_p3 + 1) in
                    let sub_seqm = Sequence.sub mine.seq f_p (l_p - f_p + 1) in


                    if (Sequence.length sub_seq1 = 0) ||  (Sequence.length sub_seq2 = 0) ||
                        (Sequence.length sub_seq3 = 0) then detect_change new_med (l_p + 1)
                    else begin                        

                        let _, median_seq, _ = Sequence.Align.readjust_3d
                            ~first_gap:false sub_seq1 sub_seq2 sub_seqm c2 c3 sub_seq3 
                        in

                        detect_change ((f_p, l_p, median_seq)::new_med) (l_p + 1)
                    end  
                end                                     
            end 
        end 
    in 

    
    let new_med_rev = detect_change [] 0 in
    let l_p, seq_ls = List.fold_left 
        (fun (l_p, seq_ls) (s, e, med) ->
             if s > l_p then begin
                 let sub_seq = Sequence.sub mine.seq l_p (s - l_p) in 
                 let seq_ls = med::(sub_seq::seq_ls) in 
                 (e + 1), seq_ls
             end else (e + 1), med::seq_ls 
        ) (0, []) (List.rev new_med_rev) 
    in
    let seq_ls = 
        if l_p < mine_len then begin
            let sub_seq = Sequence.sub mine.seq l_p (mine_len - l_p) in  
            sub_seq::seq_ls
        end else seq_ls
    in 
    
    let new_mine_seq = Sequence.concat (List.rev seq_ls) in
    let new_mine = {mine with seq = new_mine_seq} in

    let cost1, _ = cmp_cost ch1 new_mine c2 pam `Chromosome in 
    let cost2, _ = cmp_cost ch2 new_mine c2 pam `Chromosome in 
    let cost3, _ = cmp_cost ch3 new_mine c2 pam `Chromosome in 

    cost1 + cost2 + cost3, new_mine





(** [copy_chrom_map s d] copy the chromosome map
* of chromosome [s] into chromosome map of chromosome [d] *)
let copy_chrom_map s d = 
    {d with ref_code = s.ref_code; 
         ref_code1 = s.ref_code1;
         ref_code2 = s.ref_code2;
         chrom_map = s.chrom_map}
    
(** [to_single single_parent child_ref c2 pam] returns the
* single state sequence for the node [child_ref] *)
let to_single single_parent child_ref c2 pam = 
    if (single_parent.ref_code1 = -1) && (single_parent.ref_code2 = -1) then single_parent.seq
    else begin
        let gap = Cost_matrix.Two_D.gap c2 in    
        let map = List.map 
            (fun seg -> 
                 let alied_single_seq, alied_child_seq = 
                     match child_ref = single_parent.ref_code1 with
                     | true ->
                           let single, _ = Sequence.closest_alied_seq
                               seg.alied_med seg.alied_seq1 c2
                           in 
                           single, seg.alied_seq1
                     | false ->
                           let single, _  = 
                               if seg.dir2 = `Positive then 
                                   Sequence.closest_alied_seq seg.alied_med seg.alied_seq2 c2 
                               else  begin
                                   let rev_alied_seq2 = Sequence.reverse seg.alied_seq2 in 
                                   let rev_single_alied_seq2, _ = Sequence.closest_alied_seq
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
                 
                 let ungap_alied_med = Sequence.of_array (Array.of_list  ungap_alied_med) in
                 

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
        let single_seq = Sequence.concat seq_ls in   
        single_seq
    end 


(** [to_single_root root other_code c2] returns a 
* the map of the root for single state *)
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
                 let single_seg, _ = Sequence.closest_alied_seq 
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
                 let ungap_alied_med = Sequence.of_array (Array.of_list ungap_alied_med) in         
                 ungap_alied_med
            ) root.chrom_map 
        in  
        let single_seq = Sequence.concat map in 
        single_seq
    end 

(** [change_to_single med single_seq c2] 
* changes the single state of the [med] by [single_seq] *)
let change_to_single med single_seq c2 = 
    let single_seq = 
        if Sequence.length single_seq != Sequence.length med.seq then 
            if no_single_error then 
                Sequence.get_single_seq med.seq c2 
            else begin
                fprintf stdout "single_len: %i, med_len:%i\n" (Sequence.length single_seq)
                (Sequence.length med.seq); flush stdout;
                 print_endline "Single_seq";
                 Sequence.printDNA single_seq;
                 print_endline "med_seq";
                 Sequence.printDNA med.seq;
                failwith "Fail at change_to_single @ chromAli.ml" 
            end 
        else single_seq
    in  
    


    let gap = Cost_matrix.Two_D.gap c2 in
    let single_pos = ref 0 in 
    let new_map = List.map 
        (fun seg ->
             let single_alied_med = Sequence.map 
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
    let gapless_single_seq = Sequence.delete_gap single_seq in 
    {med with seq = gapless_single_seq;
         chrom_map = new_map}



    
