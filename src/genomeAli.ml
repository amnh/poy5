(* POY 5.1.1. A phylogenetic analysis program using Dynamic Homologies.       *)
(* Copyright (C) 2014 Andrés Varón, Lin Hong, Nicholas Lucaroni, Ward Wheeler,*)
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
let () = SadmanOutput.register "GenomeAli" "$Revision: 3663 $"

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
let  debug_median = false
let  debug_distance = false 
let  gap = Alphabet.gap

type direction_t = ChromPam.direction_t

(** [seg_t] is a data structure to contain a pair of segments 
* on a chromosome of the genome *)
type seg_t = {
    sta : int; (** start position on the median sequence *)
    en : int; (** end position on the median sequence *)
    cost : int; (** the cost segment *)
    med_chrom_id : int; (** median chromosome id *)
    alied_med : Sequence.s; (** aligned median sequence *)


    ref_code1 : int;    (** first child's code *)
    sta1 : int; (** start position on the first child *)
    en1 : int; (** end position on the first child *)
    chi1_chrom_id : int; (** first child chromosome id *)
    alied_seq1 : Sequence.s; (** aligned sequence of the first child *)
    dir1 : direction_t; (** sequence orientation of the first child *)

    ref_code2 : int;    (** second child's code *)
    sta2 : int; (** start position on the second child *)
    en2 : int; (** end position on the second child *)
    chi2_chrom_id : int; (** second child chromosome id *)
    alied_seq2 : Sequence.s; (** aligned sequence of the second child *)
    dir2 : direction_t; (** sequence orientation of the second child *)
}


let create_new_seg_t sta en cost med_chrom_id alied_med ref_code1
    ref_code2 sta1 sta2 en1 en2 chi1_chrom_id chi2_chrom_id alied_seq1
    alied_seq2 dir1 dir2 =
    {    
     sta=sta; en = en; 
     cost = cost; med_chrom_id = med_chrom_id; alied_med = alied_med;
     ref_code1 = ref_code1; 
     sta1 = sta1; en1 = en1; chi1_chrom_id = chi1_chrom_id;
     alied_seq1 = alied_seq1; dir1 = dir1;
     ref_code2 = ref_code2;
     sta2 = sta2; en2 = en2; chi2_chrom_id = chi2_chrom_id;
     alied_seq2 = alied_seq2; dir2 = dir2
    }


(** [chrom_t] is a data structure to contain a chromosome *)
type chrom_t = {
    chrom_id : int ref; (** homologous chromosome id *)
    main_chrom1_id : int; (** main first child's chromosome's id *)
    main_chrom2_id : int; (** main second child's chromosome's id *)

    chrom_ref_code : int; (** chromosome  reference code *)
    map : seg_t list;    
    seq : Sequence.s;   (** the chromosome sequence *)
}

let create_new_chromt sta en cost med_chrom_id alied_med ref_code1
    ref_code2 sta1 sta2 en1 en2 chi1_chrom_id chi2_chrom_id alied_seq1
    alied_seq2 dir1 dir2 =
    let newmap =
        create_new_seg_t sta en cost med_chrom_id alied_med ref_code1
        ref_code2 sta1 sta2 en1 en2 chi1_chrom_id chi2_chrom_id alied_seq1
        alied_seq2 dir1 dir2
    in
    let new_chrom_ref_code = Utl.get_new_chrom_ref_code () in
    { 
    chrom_id = ref med_chrom_id;
    main_chrom1_id = chi1_chrom_id;
    main_chrom2_id = chi2_chrom_id;
    chrom_ref_code = new_chrom_ref_code;
    map = [newmap]; 
    seq = alied_med;} 
    

(** [med_t] is data structure to contain a genome*)
type med_t = {
    chrom_arr : chrom_t array; (** chromosome array *)
    genome_ref_code : int; (** genome reference code *)
    genome_ref_code1 : int; (** the first genome child's code *)
    genome_ref_code2 : int; (** the second genome child's code *)
    cost1 :int; (** the cost from this genome to its first child *)
    recost1 : int; (** the recost from this genome to its first child *)
    cost2 : int; (** the cost from this genome to its second child *)
    recost2 : int; (** the recost from this genome to its second child *)
    (*store cost3 and sumcost in each med_t*)
    cost3 : int; (** cost(mine,ch1)+cost(mine,ch2)+cost(mine,parent) *)
    sum_cost : int; (** cost of subtree root at this med_t of this node*) 
}

(** [genome_block_t] is data structure for a block
* which is a pair of high similar sequences on both
* genomes*)
type genome_block_t = {
    block_id : int; (** the genome block ID *)
    chrom1_id : int; (** chromosome id on the first genome *)
    chrom2_id : int; (** chromosome id on the second genome *)
    block : Block.block_t; (** block information *)
}


(** [swap_seg s] swaps the first child's segment and 
* the second child's segment of segment [s] *)
let swap_seg s = {
    s with ref_code1 = s.ref_code2;
        sta1 = s.sta2;
        en1 = s.en2;
        chi1_chrom_id = s.chi2_chrom_id;
        alied_seq1 = s.alied_seq2;
        dir1 = s.dir2;

        ref_code2 = s.ref_code1;
        sta2 = s.sta1;
        en2 = s.en1;
        chi2_chrom_id = s.chi1_chrom_id;
        alied_seq2 = s.alied_seq1;
        dir2 = s.dir1;
}

(** [swap_chrom s] swaps the first child's chromosome and 
* the second child's chromosome of chromosome [ch] *)
let swap_chrom ch = {
    ch with main_chrom1_id = ch.main_chrom2_id;
        main_chrom2_id = ch.main_chrom1_id;
        map = List.map swap_seg ch.map
}

(** [swap_med s] swaps the first child's genome and 
* the second child's genome of genome [m] *)
let swap_med m = {
    m with genome_ref_code1 = m.genome_ref_code2;
        genome_ref_code2 = m.genome_ref_code1;
        cost1 = m.cost2;
        cost2 = m.cost1;
        recost1 = m.recost2;
        recost2 = m.recost1;
        chrom_arr = Array.map swap_chrom m.chrom_arr
}

let print_block b = 
    fprintf stdout "BlockID: %i, chrom1_id: %i, chrom2_id: %i\n" b.block_id
        b.chrom1_id b.chrom2_id;
    Block.print b.block;
    print_newline (); 
    flush stdout

(** [to_string med] converts all information in 
* median [med] into a string *) 
let to_string med = 
    let chromStr_arr = Array.map (fun chrom -> Sequence.to_string chrom.seq Alphabet.nucleotides)
        med.chrom_arr 
    in
    String.concat "@@" (Array.to_list chromStr_arr)


let print_genome genome = 
    fprintf stdout "Genome code: %i\n" genome.genome_ref_code;
    Array.iter (fun chrom -> 
                    fprintf stdout "chrom_ref_code=%i, chrom_id=%i:" chrom.chrom_ref_code !(chrom.chrom_id);
                    Sequence.print stdout chrom.seq Alphabet.nucleotides; 
                    fprintf stdout " @ ") genome.chrom_arr;
    print_newline ()


let get_chroms med = Array.map (fun chrom -> chrom.seq) med.chrom_arr

let print_seg seg = 
    fprintf stdout "sta: %i, en: %i, med_chrom_id: %i\n" seg.sta seg.en seg.med_chrom_id;
    Sequence.printDNA seg.alied_med;
    fprintf stdout "sta1: %i, en1: %i, chi1_chrom_id: %i\n" seg.sta1 seg.en1 seg.chi1_chrom_id;
    Sequence.printDNA seg.alied_seq1;
    fprintf stdout "sta2: %i, en2: %i, chi2_chrom_id: %i\n" seg.sta2 seg.en2 seg.chi2_chrom_id;
    Sequence.printDNA seg.alied_seq2;
    print_newline (); 
    flush stdout        

let print_chrom chrom =
    fprintf stdout "ChromId: %i, chrom_ref_code: %i, main_chrom1_id: %i, main_chrom2_id: %i \n" 
        !(chrom.chrom_id) chrom.chrom_ref_code chrom.main_chrom1_id chrom.main_chrom2_id;
    Sequence.printDNA chrom.seq;
    List.iter (fun seg -> 
                   print_seg seg;
                   print_newline ();
              ) chrom.map;
    print_endline "------------------------------------------------------------";
    flush stdout



let print_genome genome =
    Array.iter print_chrom genome.chrom_arr;
    print_endline "=============================================================";
    flush stdout

(** [chrom_map_to_string chrom_map] converts the information
* in the chromosome map [chrom_ma] into a string *)
let chrom_map_to_string chrom_map = 
    let get_dir dir = 
        match dir with
        | `Positive -> "Positive"
        | `Negative -> "Negative"
        | _ -> ""
    in 

    let str_ls = List.map 
        (fun m -> 
             let node_str =  "Sequence (chrom_id:" ^ (string_of_int m.med_chrom_id)
                 ^ ", start:" ^ (string_of_int m.sta) 
                 ^ ", end: " ^ (string_of_int m.en) ^ ")"
             in 

             let child1_str = 
                 match m.ref_code1 = -1 with
                 | true -> ""
                 | false -> 
                       ", Child1 (genome_ref:" ^ (string_of_int m.ref_code1) 
                       ^ ", chrom_id:" ^ (string_of_int m.chi1_chrom_id) 
                       ^ ", start:" ^ (string_of_int m.sta1)  
                       ^ ", end: " ^ (string_of_int m.en1)  
                       ^ ", direction:" ^ (get_dir m.dir1) ^ ")" 
             in 

             let child2_str = 
                 match m.ref_code2 = -1 with
                 | true -> ""
                 | false -> 
                       ", Child2 (genome_ref:" ^ (string_of_int m.ref_code2) 
                       ^ ", chrom_id:" ^ (string_of_int m.chi2_chrom_id) 
                       ^ ", start:" ^ (string_of_int m.sta2)  
                       ^ ", end: " ^ (string_of_int m.en2)  
                       ^ ", direction:" ^ (get_dir m.dir2) ^ ")" 
             in
             "(" ^ node_str ^ child1_str ^ child2_str ^ ")"
        ) chrom_map 
    in 
    String.concat " | " str_ls

(** [genome_map_to_string genome] converts all information
* in the genome [genome] into a string *)
let genome_map_to_string genome = 
    let chrom_map_str = Array.map 
        (fun chrom -> chrom_map_to_string chrom.map
        ) genome.chrom_arr 
    in 
    String.concat " @@ " (Array.to_list chrom_map_str) 


let ref_genome = ref {chrom_arr = [||]; 
                      genome_ref_code = -1; 
                      genome_ref_code1 = -1;
                      genome_ref_code2 = -1;
                      cost1 = 0;
                      recost1 = 0;
                      cost2 = 0;
                      recost2 = 0;
                      cost3 = 0;
                      sum_cost = 0;
}

(** [assign_hom_chrom med cost_mat user_chrom_pams]
* assigns homologous Id chrom all chromosomes in
* genome [med] *)
let assign_hom_chrom med cost_mat user_chrom_pams = 
    if !(med.chrom_arr.(0).chrom_id) = -1 then begin
        let num_chrom1 = Array.length !ref_genome.chrom_arr in
        let num_chrom2 = Array.length med.chrom_arr in
        Printf.printf "assign hom chrom, num_chrom1=%d,num_chrom2=%d\n%!"
        num_chrom1 num_chrom2;
        let chrom_cost_mat = Array.make_matrix num_chrom1 num_chrom2 0.0 in 
        for c1 = 0 to num_chrom1 - 1 do
            for c2 = 0 to num_chrom2 - 1 do
                let chrom1 = !ref_genome.chrom_arr.(c1).seq in 
                let chrom2 = med.chrom_arr.(c2).seq in 
                let chrom_med1 = ChromAli.create_med chrom1 in 
                let chrom_med2 = ChromAli.create_med chrom2 in 
                let cost, recost = (ChromAli.cmp_cost chrom_med1 chrom_med2
                                        cost_mat user_chrom_pams `Genome) 
                in
                let len = max (Sequence.length chrom1) (Sequence.length chrom2) in 
                chrom_cost_mat.(c1).(c2) <- (float cost) /. (float len);
                Printf.printf "chrom_cost_mat.(%d).(%d) <- %d/%d\n%!" c1 c2 cost len;
            done;
        done;
        (*chrom_cost_mat.(best_c1).(best_c2) has the min cost. *)
        let map1_arr = Array.make num_chrom1 false in 
        let map2_arr = Array.make num_chrom2 false in 
        for times = 0 to num_chrom2 - 1 do
            let best_c1 = ref (-1) in 
            let best_c2 = ref (-1) in
            let min_cost = ref (float Utl.large_int) in
            for c1 = 0 to num_chrom1 - 1 do
                for c2 = 0 to num_chrom2 - 1 do
                    if (map1_arr.(c1) = false) && (map2_arr.(c2) = false) 
                        && (chrom_cost_mat.(c1).(c2) < !min_cost) then begin 
                        min_cost := chrom_cost_mat.(c1).(c2); 
                        best_c1 := c1; 
                        best_c2 := c2; 
                    end;                               
                done;    
            done;  
(* chrom hom:(FLOAT) Specifies the lower limit of distance between two
chromosomes beyond which the chromosomes are not considered
to be homologous. The default value of chrom hom is 0.75*)
            let chrom_pams = ChromPam.get_chrom_pam user_chrom_pams in        
            let chrom_hom = float chrom_pams.ChromPam.chrom_hom in 
            if  (!min_cost < chrom_hom /. 100.) then begin                
                med.chrom_arr.(!best_c2).chrom_id := !best_c1;
                map1_arr.(!best_c1) <- true;
                map2_arr.(!best_c2) <- true;
                Printf.printf "min_cost=%f<%f,map1_arr.%d <- true, map2_arr.%d <- true\n%!"
                !min_cost chrom_hom !best_c1 !best_c2;
            end 
        done;
        Array.iteri 
            (fun pos mapped -> 
                 if mapped = false then begin 
                     let chrom_id = Array.length !ref_genome.chrom_arr in  
                     let chrom2 = med.chrom_arr.(pos) in  
                     chrom2.chrom_id := chrom_id; 
                     let chrom_ls = Array.to_list !ref_genome.chrom_arr in  
                     let new_chrom_arr = Array.of_list (chrom_ls @ [chrom2]) in  
                     ref_genome := {!ref_genome with chrom_arr = new_chrom_arr};
                 end 
            ) map2_arr;
        (*    
        Printf.printf "check med\n%!";
        print_genome med; Printf.printf "check ref_genome\n%!";
        print_genome !ref_genome;*)
        print_endline "End of assigning homologous chromosome in GenomeAli"; 
    end 



(** [init genome] inits a new genome
* from an array of chromosomes [genome] passed 
* from Data module *)
let init (genome  : Sequence.s Data.dyna_data) = 
    let chrom_arr = Array.map 
        (fun chrom -> {chrom_ref_code = chrom.Data.code; 
                       main_chrom1_id = -1;
                       main_chrom2_id = -1;
                       chrom_id = ref (-1);
                       map = [];
                       seq = chrom.Data.seq}
        ) genome.Data.seq_arr
    in 
    let clean_chrom_ls = List.filter 
        (fun chrom -> Sequence.length chrom.seq > 0
        ) (Array.to_list chrom_arr)
    in 
    let med = {chrom_arr = Array.of_list clean_chrom_ls; 
               genome_ref_code = Utl.get_new_genome_ref_code ();
               genome_ref_code1 = -1;
               genome_ref_code2 = -1;
               cost1 = 0;
               recost1 = 0;
               cost2 = 0;
               recost2 = 0;
                cost3 = 0;
                sum_cost = 0;
              }
    in 
    med

(** [is_leaf med] returns true if genome [med] is at 
* a leaf on the tree, othewise false *)
let is_leaf med = (med.genome_ref_code1 = -1) && (med.genome_ref_code2 = -1)

(** [create_med_from_seq chrom_arr] returns a genome
* created from an array chromosomes [chrom_arr] *)
let create_med_from_seq chrom_arr = 
    let chrom_arr = Array.map (fun chrom -> 
                                   {chrom_ref_code = -1;
                                    main_chrom1_id = -1;
                                    main_chrom2_id = -1;
                                    chrom_id = ref (-1);
                                    map = [];
                                    seq = chrom;
                                   }
                              ) chrom_arr
    in 
    let med = {chrom_arr = chrom_arr;
               genome_ref_code = Utl.get_new_genome_ref_code ();
               genome_ref_code1 = -1;
               genome_ref_code2 = -1;
               cost1 = 0;
               recost1 = 0;
               cost2 = 0;
               recost2 = 0;
               cost3 = 0;
                sum_cost = 0;
              }
    in 
    med    


(** [find_conserved_areas seq1 seq2 cost_mat chrom_pams]
* returns a global map containing conserved areas between
* two sequences [seq1] and [seq2] *)
let find_conserved_areas seq1 seq2 cost_mat chrom_pams =
    let ali_pam = ChromPam.get_chrom_pam chrom_pams in      
    let len1 = Sequence.length seq1 in
    let len2 = Sequence.length seq2 in 
    let ali_pam = {ali_pam with   
                       ChromPam.min_pos1 = 0;  
                       ChromPam.max_pos1 = len1 - 1;  
                       ChromPam.min_pos2 = 0;  
                       ChromPam.max_pos2 = len2 - 1;  
                       negative = false;
                  }  
    in   
    let global_map, _, _ = ChromAli.create_global_map seq1 seq2 cost_mat ali_pam in   
    global_map

(** [find_chrom chrom_id genome] returns the chromosome [chrom_id]
* in the genome [genome] *)
let find_chrom chrom_id genome =
    try
        Some (List.find 
                  (fun chrom ->  !(chrom.chrom_id) = chrom_id) (Array.to_list genome.chrom_arr) 
             ) 
    with | Not_found -> None

let get_chrom_from_arr idx chrom_arr =
    let len = Array.length chrom_arr in
    if idx<len then
        Some (chrom_arr.(idx))
    else
        None

(** [create_loci seq subseq_ls] returns a list of 
* loci where ends of subsequence list [subseq_ls] 
* is used as milestones to divide the sequence [seq] *)
let create_loci seq subseq_ls =
    let subseq_arr = Array.of_list subseq_ls in 
    Array.sort (fun s1 s2 -> compare s1.Subseq.sta s2.Subseq.sta) subseq_arr;

    let loci_ls = ref [] in 
    let loci_id = ref 1 in 
    let sta = ref 0 in 
    Array.iter (fun s -> 
                    if !sta < s.Subseq.sta  then begin
                        let ns = {Subseq.id = !loci_id;
                                  Subseq.sta = !sta;
                                  Subseq.en = s.Subseq.sta - 1;
                                  Subseq.block_id_ls = []
                                 } 
                        in 
                        incr loci_id;
                        loci_ls := ns::!loci_ls;
                    end;
                    loci_ls := {s with Subseq.id = !loci_id}::!loci_ls;
                    incr loci_id;                     
                    sta := s.Subseq.en + 1;
               ) subseq_arr;  

    let seq_len = Sequence.length seq in 
    if !sta < seq_len then begin
        let ns = {Subseq.id = !loci_id; 
                  Subseq.sta = !sta; 
                  Subseq.en = seq_len - 1;
                  Subseq.block_id_ls = [] 
                 }   
        in 
        incr loci_id;
        loci_ls := ns::!loci_ls
    end;
    List.rev (!loci_ls)
    

(** create_fast_general_ali chrom_id genome1_ref_code chrom1_seq loci1_ls
*        genome2_ref_code chrom2_seq loci2_ls gb_ls cost_mat ali_pam] 
* returns the median sequence between chromosomes [chrom_id] of genome
* [genome1_ref_code] and [genome2_ref_code] *)
let create_fast_general_ali chrom_id genome1_ref_code chrom1_seq loci1_ls
        genome2_ref_code chrom2_seq loci2_ls gb_ls cost_mat ali_pam =
    let loci1_arr = Array.of_list loci1_ls in 
    let loci2_arr = Array.of_list loci2_ls in 
    let len1 = Array.length loci1_arr in 
    let len2 = Array.length loci2_arr in 
    Printf.printf "create_fast_general_ali loci arr size = %d,%d\n%!" len1 len2;
    let use_ukk = 
                match !Methods.algn_mode with
                        | `Algn_Newkk -> true
                                | _ -> false
                                    in
    Array.iter (fun sq1 -> sq1.Subseq.id <- sq1.Subseq.id  * 2 - 1) loci1_arr; 
    Array.iter (fun sq2 -> sq2.Subseq.id <- len1 * 2 + sq2.Subseq.id  * 2 - 1) loci2_arr;
    let gen_gap_code = (len1 + len2 + 1) * 2 - 1 in     
    let gen_c2 = Array.make_matrix  (gen_gap_code + 1) (gen_gap_code + 1) Utl.large_int in 
    for r1 = 0 to len1 - 1 do
        let sq1 = loci1_arr.(r1) in 
        let sq1_id = sq1.Subseq.id in 
        let s1 = sq1.Subseq.sta and e1 = sq1.Subseq.en in 
        let sq1_seq = Sequence.sub chrom1_seq s1 (e1 - s1 + 1) in 
        let free1 = Subseq.is_free sq1 in 
        for r2 = 0 to len2 - 1 do
            let sq2 = loci2_arr.(r2) in 
            let sq2_id = sq2.Subseq.id in 
            let s2 = sq2.Subseq.sta and e2 = sq2.Subseq.en in 
            let sq2_seq = Sequence.sub chrom2_seq s2 (e2 - s2 + 1) in 
            let free2 = Subseq.is_free sq2 in
            if free1 && free2 then begin
                let alied_sq1, alied_sq2, cost, _ = 
                    Sequence.align2 sq1_seq sq2_seq cost_mat use_ukk
                in 
                gen_c2.(sq1_id).(sq2_id) <- cost;
                gen_c2.(sq2_id).(sq1_id) <- cost;
            end; 
        done
    done;
    let pair_gap subseq chrom_seq = 
       let id = subseq.Subseq.id in 
       let ss = Subseq.get_subseq chrom_seq subseq in
       let del_cost = Sequence.cmp_gap_cost
           ali_pam.ChromPam.chrom_indel_cost ss
       in
       gen_c2.(id).(gen_gap_code) <-  del_cost;
       gen_c2.(gen_gap_code).(id) <-  del_cost;
    in
    gen_c2.(gen_gap_code).(gen_gap_code) <- 0;
    Array.iter (fun sq1 -> if Subseq.is_free sq1 then pair_gap sq1 chrom1_seq) loci1_arr;
    Array.iter (fun sq2 -> if Subseq.is_free sq2 then pair_gap sq2 chrom2_seq) loci2_arr;
    let free_id1_ls = List.fold_right 
        (fun locus free_id1_ls -> 
             if Subseq.is_free locus then locus.Subseq.id::free_id1_ls
             else free_id1_ls
        ) loci1_ls []
    in 
    let free_id2_ls = List.fold_right 
        (fun locus free_id2_ls -> 
             if Subseq.is_free locus then locus.Subseq.id::free_id2_ls
             else free_id2_ls
        ) loci2_ls []
    in 
    let free_id1_arr = Array.of_list free_id1_ls in
    let free_id2_arr = Array.of_list free_id2_ls in
    let swap_med = ali_pam.ChromPam.swap_med in 
    let total_cost, (recost1, recost2), alied_free_id1, alied_free_id2 = 
        GenAli.create_gen_ali_code         
        ali_pam.ChromPam.kept_wag
        `Genome free_id1_arr free_id2_arr gen_c2 gen_gap_code 
        ali_pam.ChromPam.re_meth swap_med ali_pam.ChromPam.circular false
        use_ukk
    in   
    let max_sq2_id = gen_gap_code - 2 in 
    let mark2_arr = Array.make (max_sq2_id + 1) true in 
    Array.iter (fun locus2 -> mark2_arr.(locus2.Subseq.id) <- false) loci2_arr;
    Array.iteri (fun idx _ -> 
                     (if (alied_free_id1.(idx) != gen_gap_code) && (alied_free_id2.(idx) != gen_gap_code) then begin
                          mark2_arr.(alied_free_id2.(idx)) <- true
                      end );
                ) alied_free_id1;    
    Array.iter (fun sq2 -> 
                    let is_free = Subseq.is_free sq2 in
                    (if is_free = false then begin
                         let gb_id = List.hd sq2.Subseq.block_id_ls in 
                         let gb1 = List.find (fun gb1 -> gb1.block_id = gb_id) gb_ls in 
                         (if gb1.chrom2_id = chrom_id then 
                              mark2_arr.(sq2.Subseq.id) <- true);
                     end);  
               ) loci2_arr;
    let rec add_missing_sq2 sq2_id (seg_ls : seg_t list) =
        if (sq2_id > len1 * 2 ) && (mark2_arr.(sq2_id) = false) then begin
            let sq2 = loci2_arr.( (sq2_id - len1 * 2) / 2 ) in 
            let sta2 = sq2.Subseq.sta and en2 = sq2.Subseq.en in 
            let sq2_seq = Sequence.sub chrom2_seq sta2 (en2 - sta2 + 1) in  
            let seg = 
                {sta = -1; en = -1;
                 cost = gen_c2.(gen_gap_code).(sq2_id);
                 med_chrom_id = chrom_id;  
                 alied_med = Sequence.create_gap_seq (Sequence.length sq2_seq);
                 ref_code1 = genome1_ref_code; sta1 = -1;  en1 = -1;  
                 chi1_chrom_id = chrom_id;  
                 alied_seq1 = Sequence.create_gap_seq (Sequence.length sq2_seq); dir1 = `Positive;  
                 ref_code2 = genome2_ref_code; sta2 = sta2;  en2 = en2;  
                 chi2_chrom_id = chrom_id;  alied_seq2 = sq2_seq; dir2 = `Positive
                }                         
            in 
            assert (Sequence.length seg.alied_seq1 > 0);
            assert (Sequence.length seg.alied_seq2 > 0); 
            add_missing_sq2 (sq2_id - 2) (seg::seg_ls)
        end else seg_ls
    in 
    let alied_free_len = Array.length alied_free_id1 in     
    let search_sq2 sq1 = 
        let sq1_id = sq1.Subseq.id in
        let rec move p =             
            if p = alied_free_len then (-1)
            else if alied_free_id1.(p) = sq1_id then alied_free_id2.(p)
            else move (p+1)
        in 
        let sq2_id =  move 0 in
        if sq2_id = gen_gap_code then chrom_id, gen_gap_code, `Positive
        else if sq2_id > 0 then chrom_id, sq2_id, `Positive
        else begin (* sq2_id = -1 *)
            let gb_id = List.hd sq1.Subseq.block_id_ls in  
            let gb1 = List.find (fun gb1 -> gb1.block_id = gb_id) gb_ls in              
            let b = gb1.block in 
            if gb1.chrom2_id = chrom_id then begin
                let sq2 = List.find 
                    (fun locus2 -> locus2.Subseq.sta = b.Block.sta2) loci2_ls 
                in
                chrom_id, sq2.Subseq.id, b.Block.direction
            end else gb1.chrom2_id, -1, `Positive
        end 
    in   
    let chrom_med_len = ref 0 in 
    let med_ls, seg_ls, chrom_cost = 
        Array.fold_left (fun (med_ls, seg_ls, chrom_cost) sq1 ->
             let sq1_id = sq1.Subseq.id in
             let s1 = sq1.Subseq.sta and e1 = sq1.Subseq.en in 
             let sq1_seq = Sequence.sub chrom1_seq s1 (e1 - s1 + 1) in 
             match Subseq.is_free sq1 with
             | true ->
                   let _, sq2_id, _ = search_sq2 sq1 in
                   if sq2_id = gen_gap_code then begin 
                       let sq2_seq = Sequence.create_gap_seq (Sequence.length sq1_seq) in 
                       let med = Sequence.create_median_gap sq1_seq cost_mat in
                       let med_len = Sequence.cmp_num_not_gap med in  
                       let sta, en = 
                           match med_len with 
                           | 0 -> -1, -1
                           | _ -> !chrom_med_len, !chrom_med_len + med_len - 1
                       in 
                       let seg = 
                           {sta = sta; en = en;
                            cost = gen_c2.(sq1_id).(gen_gap_code);  
                            med_chrom_id = chrom_id;  
                            alied_med = med;
                            ref_code1 = genome1_ref_code; sta1 = s1;  en1 = e1;  
                            chi1_chrom_id = chrom_id;  alied_seq1 = sq1_seq;  dir1 = `Positive;  
                            ref_code2 = -1; sta2 = -1;  en2 = -1;  
                            chi2_chrom_id = chrom_id;  alied_seq2 = sq2_seq; dir2 = `Positive
                           }                        
                       in  
                       assert (Sequence.length seg.alied_seq1 > 0); 
                       assert (Sequence.length seg.alied_seq2 > 0); 
                       chrom_med_len := !chrom_med_len + med_len;
                       med::med_ls, seg::seg_ls, chrom_cost
                    end else begin
                        let sq2 = loci2_arr.( (sq2_id - len1 * 2) / 2 ) in 
                        let s2 = sq2.Subseq.sta and e2 = sq2.Subseq.en in 
                        let sq2_seq = Sequence.sub chrom2_seq s2 (e2 - s2 + 1) in  
                        let alied_seq1, alied_seq2, _, _ = Sequence.align2 sq1_seq
                            sq2_seq cost_mat use_ukk
                        in 
                        let med, cost = Sequence.create_median_seq alied_seq1 alied_seq2 cost_mat in
                        let med_len = Sequence.cmp_num_not_gap med in  
                        let sta, en = 
                            match med_len with 
                            | 0 -> -1, -1
                            | _ -> !chrom_med_len, !chrom_med_len + med_len - 1
                        in 
                        let seg = 
                            {sta = sta; en = en;
                             cost = cost; med_chrom_id = chrom_id; 
                             alied_med = med;
                             ref_code1 = genome1_ref_code; sta1 = s1;  en1 = e1; 
                             chi1_chrom_id = chrom_id; alied_seq1 = alied_seq1; dir1 = `Positive;  
                             ref_code2 = genome2_ref_code; sta2 = s2; en2 = e2; 
                             chi2_chrom_id = chrom_id;  alied_seq2 = alied_seq2; dir2 = `Positive;                          
                            }
                        in 
                        assert (Sequence.length seg.alied_seq1 > 0);
                        assert (Sequence.length seg.alied_seq2 > 0); 
                        chrom_med_len := !chrom_med_len + med_len; 
                        let mis_seg2_ls = add_missing_sq2 (sq2_id - 1) [] in 
                        med::med_ls, (List.append mis_seg2_ls (seg::seg_ls)), chrom_cost
                    end 
             | false ->
                   let gb_id = List.hd sq1.Subseq.block_id_ls in 
                   let gbcp4 = List.find (fun gbcp4 -> gbcp4.block_id = gb_id) gb_ls in 
                   let b = gbcp4.block in 
                   let alied_seq1 =  Utl.deref b.Block.alied_seq1 in
                   let alied_seq2=  Utl.deref b.Block.alied_seq2 in
                   let med, cost = Sequence.create_median_seq alied_seq1 
                       alied_seq2 cost_mat  
                   in  
                   let med_len = Sequence.cmp_num_not_gap med in  
                   let sta, en = 
                       match med_len with 
                       | 0 -> -1, -1
                       | _ -> !chrom_med_len, !chrom_med_len + med_len - 1
                   in 
                   let seg =  
                      {sta = sta; en = en;
                       cost = cost; med_chrom_id = chrom_id; 
                       alied_med = med;
                       ref_code1 = genome1_ref_code;
                       sta1 = b.Block.sta1; en1 = b.Block.en1; chi1_chrom_id = gbcp4.chrom1_id;
                       alied_seq1 = alied_seq1; dir1 = `Positive;
                       ref_code2 = genome2_ref_code;
                       sta2 = b.Block.sta2; en2 = b.Block.en2; chi2_chrom_id = gbcp4.chrom2_id;
                       alied_seq2 = alied_seq2; dir2 = b.Block.direction;
                      }
                   in 
                   assert (Sequence.length seg.alied_seq1 > 0);
                   assert (Sequence.length seg.alied_seq2 > 0); 
                   chrom_med_len := !chrom_med_len + med_len; 
                   med::med_ls, seg::seg_ls, chrom_cost + cost 
        ) ([], (add_missing_sq2 max_sq2_id []), total_cost) loci1_arr 
    in 
    (*calculate chromosome translocation cost -- we used to call it
    * chrom_breakpoint cost.*)
    let chrom_bp = ali_pam.ChromPam.translocation in 
    let bp =
        match ali_pam.ChromPam.re_meth with
        | `Locus_Breakpoint bp -> bp
        | _ -> failwith "Can not calculate inversion distance for multiple chromosome!"
    in 
    let cmp_breakpoint_cost = 
        let total_bp = ref 0 in
        for p = 0 to len1 - 2 do
            let chrom21_id, sq21_id, dir21 = search_sq2 loci1_arr.(p) in 
            let chrom22_id, sq22_id, dir22 = search_sq2 loci1_arr.(p+1) in
            let one_bp = 
                if (chrom21_id != chrom_id) || (chrom22_id != chrom_id) then chrom_bp
                else if (sq21_id = gen_gap_code) || (sq22_id = gen_gap_code) then 0
                else if ( (sq21_id + 1 = sq22_id) && (dir21 = `Positive) && (dir22 = `Positive) ) ||
                    ( (sq21_id - 1 = sq22_id) && (dir21 = `Negative)  && (dir22 = `Negative) ) then 0
                else  bp
            in 
            total_bp := !total_bp + one_bp
        done;
        !total_bp
    in 
    let total_bp_cost = cmp_breakpoint_cost in
    let med_ls = List.rev med_ls in 
    let seg_ls = List.rev seg_ls in 
    let med = Sequence.concat med_ls in 
    let edit_cost = chrom_cost - (recost1 + recost2) in
    Printf.printf "return edit_cost = %d, total_bp_cost = %d\n%!" edit_cost
    total_bp_cost;
    med, seg_ls, (edit_cost + total_bp_cost), 0, total_bp_cost



(** [create_chrom_med (genome1_ref_code, chrom1) 
*    (genome2_ref_code, chrom2) gb_ls cost_mat chrom_pams]
* return a chromosme median between chromosome [chrom1]
* of genome [genome1_ref_code] and chromosome [chrom2] of
* genome [genome2_ref_code] *)
let create_chrom_med (genome1_ref_code, chrom1) (genome2_ref_code, chrom2) gb_ls cost_mat chrom_pams = 
    Printf.printf "create_chrom_med\n%!";
    let subseq1_ls = 
        List.fold_left  
        (fun subseq1_ls gbcp4 -> 
             if gbcp4.chrom1_id = !(chrom1.chrom_id) then begin
                 let subseq1 = 
                     {Subseq.id = -1;
                      Subseq.sta = gbcp4.block.Block.sta1; 
                      Subseq.en = gbcp4.block.Block.en1; 
                      Subseq.block_id_ls = [gbcp4.block_id]} 
                 in 
                 subseq1::subseq1_ls
             end else subseq1_ls
        ) [] gb_ls
    in 
    let subseq2_ls = 
        List.fold_left  
        (fun subseq2_ls gbcp4 -> 
             if gbcp4.chrom2_id = !(chrom2.chrom_id) then begin
                 let subseq2 = 
                     {Subseq.id = -1;  
                      Subseq.sta = gbcp4.block.Block.sta2; 
                      Subseq.en = gbcp4.block.Block.en2; 
                      Subseq.block_id_ls = [gbcp4.block_id]} 
                 in 
                 subseq2::subseq2_ls
             end else subseq2_ls
        ) [] gb_ls        
    in 
    let loci1_ls = create_loci chrom1.seq subseq1_ls in 
    let loci2_ls = create_loci chrom2.seq subseq2_ls in 
    let chrom_id = !(chrom1.chrom_id) in
    let med, seg_ls, chrom_cost, recost1, recost2 = create_fast_general_ali chrom_id
        genome1_ref_code chrom1.seq loci1_ls genome2_ref_code chrom2.seq loci2_ls gb_ls cost_mat chrom_pams 
    in 
    let chrom_ref_code = Utl.get_new_chrom_ref_code () in 
    let chrom_med = {
        chrom_id = ref chrom_id; 
        main_chrom1_id = !(chrom1.chrom_id);
        main_chrom2_id = !(chrom2.chrom_id);
        chrom_ref_code = chrom_ref_code; 
        seq = med; 
        map = seg_ls}
    in 
    chrom_med, chrom_cost, recost1, recost2


    



(** [create_genome_blocks med1 med2 cost_mat chrom_pams]
* returns a list of detected genome blocks between 
* genome [med1] and genome [med2] *)
let create_genome_blocks med1 med2 cost_mat chrom_pams =
    let num_chrom1 = Array.length med1.chrom_arr in  
    let num_chrom2 = Array.length med2.chrom_arr in 
    let max_chrom = Array.length !ref_genome.chrom_arr in 
    let gb_ls = ref [] in 
    let gb_id = ref 0 in 
    for chrom_id = 0 to max_chrom - 1 do
        let chrom1_opt = find_chrom chrom_id med1 in 
        let chrom2_opt = find_chrom chrom_id med2 in
        Printf.printf "find chrom id#%d,%!" chrom_id;
        match chrom1_opt, chrom2_opt with
        | Some chrom1, Some chrom2 ->
              Printf.printf "create global map,add new block to gb_ls with block_id = %d\n%!" !gb_id;
              let blocks = find_conserved_areas chrom1.seq chrom2.seq cost_mat chrom_pams in
              gb_ls := List.fold_left 
                  (fun gb_ls block -> 
                       let gbcp4 = {block_id = !gb_id; chrom1_id  = chrom_id;
                                 chrom2_id = chrom_id; block = block}
                       in 
                       incr gb_id;
                       gbcp4::gb_ls
                  ) !gb_ls blocks;                               
        | _, _ -> ()    
    done;
    let is_overlaped gbcp4 gb_ls = 
        let loci_overlaped s1 e1 s2 e2 = not (e1 < s2) || (e2 < s1) in 
        List.fold_left
            (fun overlaped gb2 -> 
                 let b1 = gbcp4.block and b2 = gb2.block in 
                 if (gbcp4.chrom1_id = gb2.chrom1_id) &&
                     (loci_overlaped b1.Block.sta1 b1.Block.en1 b2.Block.sta1 b2.Block.en1) then true
                 else if (gbcp4.chrom2_id = gb2.chrom2_id) && 
                     (loci_overlaped b1.Block.sta2 b1.Block.en2 b2.Block.sta2 b2.Block.en2) then true
                 else overlaped
            ) false gb_ls  
    in 
    let mapped_arr = Array.make max_chrom 0 in 
    Array.iter (fun chrom -> 
                    let id = !(chrom.chrom_id) in 
                    mapped_arr.(id) <- mapped_arr.(id) + 1
               ) med1.chrom_arr;
    Array.iter (fun chrom -> 
                    let id = !(chrom.chrom_id) in 
                    mapped_arr.(id) <- mapped_arr.(id) + 1
               ) med2.chrom_arr;
    for c1 = 0 to num_chrom1 -  1 do
        let chrom1 = med1.chrom_arr.(c1) in 
        let chrom1_id = !(chrom1.chrom_id) in 
        for c2 = 0 to num_chrom2 - 1 do
            let chrom2 = med2.chrom_arr.(c2) in 
            let chrom2_id = !(chrom2.chrom_id) in 
            if (chrom1_id != chrom2_id) && (mapped_arr.(chrom1_id) = 2) &&
                (mapped_arr.(chrom2_id) = 2) then begin
                let blocks = find_conserved_areas chrom1.seq chrom2.seq 
                    cost_mat chrom_pams
                in  
                gb_ls := List.fold_left  
                  (fun gb_ls block -> 
                       let gbcp4 = {block_id = !gb_id; chrom1_id = !(chrom1.chrom_id);
                                 chrom2_id = !(chrom2.chrom_id); block = block}
                       in                        
                       incr gb_id;
                       if is_overlaped gbcp4 gb_ls then gb_ls
                       else gbcp4::gb_ls
                  ) !gb_ls blocks                
            end 
        done
    done;
    !gb_ls
    
        
let create_med_vinh_annotator med1 med2 cost_mat user_chrom_pams =
    assign_hom_chrom med1 cost_mat user_chrom_pams;
    assign_hom_chrom med2 cost_mat user_chrom_pams;
    let gb_ls = create_genome_blocks med1 med2 cost_mat user_chrom_pams in 
    let max_chrom = Array.length !ref_genome.chrom_arr in 
    let rev_chrom_med_ls = ref [] in 
    let g_cost = ref 0 in 
    let g_recost1 = ref 0 in 
    let g_recost2 = ref 0 in 
    let ali_pam = ChromPam.get_chrom_pam user_chrom_pams in 
    let child1_genome_ref_code = med1.genome_ref_code in 
    let child2_genome_ref_code = med2.genome_ref_code in 
    for chrom_id = 0 to max_chrom - 1 do 
        let chrom1_opt = find_chrom chrom_id med1 in 
        let chrom2_opt = find_chrom chrom_id med2 in 
        Printf.printf "work on chromid # %d," chrom_id;
        match chrom1_opt, chrom2_opt with
        | Some chrom1, Some chrom2 ->
                Printf.printf "create chrom median for chrom1 and chrom2\n%!";
              let chrom_med,  chrom_cost, recost1, recost2 = create_chrom_med
                  (child1_genome_ref_code, chrom1)
                  (child2_genome_ref_code, chrom2) gb_ls cost_mat ali_pam
              in 
              rev_chrom_med_ls := chrom_med::!rev_chrom_med_ls; 
              g_cost := !g_cost + chrom_cost;
              g_recost1 := !g_recost1 + recost1;
              g_recost2 := !g_recost2 + recost2;
        | Some chrom1, None -> 
                Printf.printf "add indel chrom for chrom2\n%!";
              let med_seq = chrom1.seq in 
              let med_len = Sequence.cmp_num_not_gap med_seq in 
              let chrom1_len = Sequence.length chrom1.seq in  
              let indel_cost = Sequence.cmp_gap_cost
                  ali_pam.ChromPam.chrom_indel_cost chrom1.seq
              in
              let seg =   
                  {sta = 0; en = med_len - 1; 
                   cost = indel_cost; med_chrom_id = chrom_id; 
                   alied_med = med_seq;
                   ref_code1 = child1_genome_ref_code; sta1 = 0;  en1 = chrom1_len - 1;
                   chi1_chrom_id = !(chrom1.chrom_id); alied_seq1 = chrom1.seq; dir1 = `Positive;  
                   ref_code2 = -1; sta2 = -1; en2 = -1;  
                   chi2_chrom_id = -1;  alied_seq2 = Sequence.create_gap_seq (Sequence.length chrom1.seq); 
                    dir2 = `Positive;                          
                  } 
              in  
              assert (Sequence.length seg.alied_seq1 > 0);
              assert (Sequence.length seg.alied_seq2 > 0); 
              let chrom_med = {chrom_ref_code = Utl.get_new_chrom_ref_code ();
                               main_chrom1_id = !(chrom1.chrom_id);
                               main_chrom2_id = -1;
                               chrom_id = ref chrom_id; 
                               map = [seg];
                               seq = med_seq}
              in 
              g_cost := !g_cost + indel_cost;
              rev_chrom_med_ls := chrom_med::!rev_chrom_med_ls;                   
        | None, Some chrom2 -> 
                Printf.printf "add indel chrom for chrom1\n%!";
              let med_seq = chrom2.seq in
              let med_len = Sequence.cmp_num_not_gap med_seq in 
              let chrom2_len = Sequence.length chrom2.seq in 
              let indel_cost = Sequence.cmp_gap_cost
                  ali_pam.ChromPam.chrom_indel_cost chrom2.seq
              in
              let seg =  
                  {sta = 0; en = med_len - 1;
                   cost = indel_cost; med_chrom_id = chrom_id; 
                   alied_med = med_seq;
                   ref_code1 = -1; sta1 = -1;  en1 = -1;
                   chi1_chrom_id = -1; alied_seq1 = Sequence.create_gap_seq (Sequence.length chrom2.seq); dir1 = `Positive;  
                   ref_code2 = child2_genome_ref_code; sta2 = 0;  en2 = chrom2_len - 1;
                   chi2_chrom_id = chrom_id;  alied_seq2 = chrom2.seq; dir2 = `Positive;                          
                  } 
              in  
              assert (Sequence.length seg.alied_seq1 > 0);
              assert (Sequence.length seg.alied_seq2 > 0); 
              let chrom_med = {chrom_ref_code = Utl.get_new_chrom_ref_code ();
                               main_chrom1_id = -1;
                               main_chrom2_id = !(chrom2.chrom_id);
                               chrom_id = ref chrom_id; 
                               map = [seg];
                               seq = med_seq}
              in 
              g_cost := !g_cost + indel_cost;
              rev_chrom_med_ls := chrom_med::!rev_chrom_med_ls;               
        | None, None -> ()
    done;
    let chrom_med_ls = 
        List.map (fun chrom -> 
                      {chrom with seq = (Sequence.delete_gap chrom.seq)}
                 ) (List.rev !rev_chrom_med_ls)
    in 
    let chrom_med_ls = List.filter 
        (fun chrom -> Sequence.length chrom.seq > 0) chrom_med_ls
    in 
    let genome_med = {chrom_arr = Array.of_list chrom_med_ls; 
                      genome_ref_code = Utl.get_new_genome_ref_code ();
                      genome_ref_code1 = med1.genome_ref_code;
                      genome_ref_code2 = med2.genome_ref_code;
                      cost1 = !g_cost - !g_recost2;
                      recost1 = !g_recost1;
                      cost2 = !g_cost - !g_recost1;
                      recost2 = !g_recost2;
                      cost3 = 0; (* cost3 is 0 for median2, until we update it with median3 function*)
                      sum_cost = !g_cost + med1.sum_cost + med2.sum_cost;
                      }
    in 
    if (List.length chrom_med_ls = 0) then failwith "Created an empty genome";
    genome_med, !g_cost, (!g_recost1, !g_recost2)
 


let concat_chrom_seq chrom_arr =
    let revseqlst,revrangelst =
    Array.fold_left (fun (accseq,accrange) chromt ->
        let seqlen = Sequence.length chromt.seq in
        let resrange = match accrange with
        | (startpos,endpos)::t ->
                (endpos+1,endpos+1+seqlen-1)::accrange
        | [] ->
                [(0,0+seqlen-1)]
        in
        chromt.seq::accseq, resrange
    ) ([Sequence.get_empty_seq ()],[]) chrom_arr 
    in
    Sequence.concat (List.rev revseqlst), List.rev revrangelst


let create_med_mauve_annotator med1 med2 cost_mat ali_pam =
    let num_chrom1 = Array.length med1.chrom_arr in  
    let num_chrom2 = Array.length med2.chrom_arr in
    let max_num_chrom = 
        if num_chrom1>num_chrom2 then num_chrom1 else num_chrom2
    in
    if debug_median then begin
        Printf.printf "create_med_mauve_annotator,num_chrom = %d,%d,med1/med2=\n%!" 
        num_chrom1 num_chrom2;
        (*print_genome med1; print_genome med2;*)
        Printf.printf "call concat_chrom_seq on med1 and on med2\n%!";
    end;
    let seq1,rangelst1 = concat_chrom_seq med1.chrom_arr
    and seq2,rangelst2 = concat_chrom_seq med2.chrom_arr 
    in
    if debug_median then begin 
        let print_rangelst lst = 
        List.iter (fun (l,r) ->
        Printf.printf "(l=%d,r=%d);%!" l r ) lst;
        print_newline()
        in
        Printf.printf "check rangelst after concat chrom seq\n%!";
        print_rangelst rangelst1;
        print_rangelst rangelst2;
    end;
    let len1,len2 = Sequence.length seq1,Sequence.length seq2 in
    let ali_pam = {ali_pam with  
                       ChromPam.min_pos1 = 0; 
                       ChromPam.max_pos1 = len1 - 1; 
                       ChromPam.min_pos2 = 0; 
                       ChromPam.max_pos2 = len2 - 1; 
                  } 
    in
    let code_range_lstlst,gen_gap_code, alignment_matrix, 
    alied_code1_lst, alied_code2_lst, total_cost, (recost1, recost2) =
    AliMap.create_general_ali_mauve seq1 seq2 cost_mat ali_pam None None
    in
    let edit_and_indel_cost = total_cost - recost1 - recost2 in
    if debug_median then begin
        Printf.printf "create_general_ali_mauve return with cost=%d(rc1=%d,rc2=%d)\n%!"
        total_cost recost1 recost2;
        List.iter (fun lst ->
            List.iter (fun ( code,(l,r),lcbkey ) ->
                Printf.printf " code=%d,(l=%d,r=%d)," code l r;
                Utl.printIntArr (Array.of_list lcbkey);
            ) lst;
            print_newline();
        ) code_range_lstlst;
    end;
    (*function for id translocation of loci between chromosomes.*)
    let in_which_chromosome leftend rightend rangelst = 
        let idx = ref (-1) in
        List.fold_left (fun accset (l,r) ->
            idx := !idx + 1;
            if (leftend>=l && rightend<=r) || 
            (leftend<l && rightend>=l && rightend<=r) ||
            (leftend>=l && leftend<=r && rightend>r)  
            then All_sets.Integers.add !idx accset 
            else accset
        ) All_sets.Integers.empty rangelst 
    in        
    let translocation,lcb2remove = List.fold_left (fun (trloc,lcb2rmv) ( code,(left,right),lcbkey ) ->
        if (List.length lcbkey)>0 then begin
            let show_up_in_seq1 =  
                in_which_chromosome left right rangelst1 in
            let _,(l,r),_ = 
            try 
                List.find (fun (_,(_,_),lk) -> (lk=lcbkey) ) (List.nth code_range_lstlst 1)
            with | Not_found -> failwith "could not find lcb in second code lst"
            in
            let show_up_in_seq2 = in_which_chromosome l r rangelst2 in
            if debug_median then begin
                Printf.printf "lcb with code=%d,lcbkey=%!" code;
                Utl.printIntArr (Array.of_list lcbkey); 
                Printf.printf "in seq1 shows up in chrom:";
                All_sets.Integers.iter (Printf.printf " %d,") show_up_in_seq1;
                Printf.printf "; in seq2 shows up in chrom:";
                All_sets.Integers.iter (Printf.printf " %d,") show_up_in_seq2;
                print_newline();
            end;
            let uniset = All_sets.Integers.union show_up_in_seq2 show_up_in_seq1
            and interset = All_sets.Integers.inter show_up_in_seq2 show_up_in_seq1
            in
            let uniset = All_sets.Integers.fold (fun chromNO accset ->
                All_sets.Integers.remove chromNO accset
            ) interset uniset in
            All_sets.Integers.fold (fun crhomNO (acc,acclcb) ->
                (*if (List.mem chromNO show_up_in_seq2) then (acc,acclcb)*)
                if (All_sets.Integers.cardinal uniset)=0 then acc,acclcb
                else acc+1,lcbkey::acclcb
            ) show_up_in_seq1 (trloc,lcb2rmv) 
        end
        else trloc,lcb2rmv;
    ) (0,[]) (List.hd code_range_lstlst) in
    (* recost might be result of translocation,but we already counting
    * translocation into result.
    * for example, we have four locus 1,2,3,4. loci2 is being translocated into
    * second chromosome:
    * [1;2;@3;4] and [1;3;@2;4] 
    * recost2 will be bigger than 0, but if we remove loci 2, there is no
    * breakpoint/breakinv cost at all.
    * also the code in code_range_lstlst is not a permutation of each other, 
    * like [1;3;5;7] and [9;11;13;15] , we need to standardize it. there is a
    * "standardize" funcion in utlGrapp, but since the code array we pass to
    * utlGrapp does not include translocation lcb, that function won't work for us. 
    * we have to do it here.*)
    let codearr1_len = List.length (List.hd code_range_lstlst) in
    let base = codearr1_len * 2 in
    let codelst_after_rmv = List.map (fun lst ->
        let reslst = List.filter (fun (_,_,lcbkey) -> 
           (List.length lcbkey)>0 && (List.mem lcbkey lcb2remove)=false 
        ) lst in
        List.map (fun (x,_,_) -> 
            if x>base then x-base else x ) reslst;
    ) code_range_lstlst in
    let arr1,arr2 = 
        Array.of_list (List.hd codelst_after_rmv),
        Array.of_list (List.nth codelst_after_rmv 1) in
    let real_recost = GenAli.cmp_recost_simple arr1 arr2 ali_pam.ChromPam.re_meth 
    ali_pam.ChromPam.circular in
    let real_recost2 = real_recost in 
    if debug_median then begin
        Printf.printf "translocation number = %d,lcb code 2 remove = %!" translocation;
        List.iter (fun lcbkey -> Utl.printIntArr (Array.of_list lcbkey)) lcb2remove;  print_newline();
        Printf.printf "after removing, recost = %d, code arr1/code arr2=\n%!" real_recost2;
        Utl.printIntArr arr1; Utl.printIntArr arr2;
    end;
    (*about the median, what are we going to do here? the alignment seq and cost
    * we got from mauve are for the whole genome, not for each chromosome. 
    * cut it? how? does it even make any sense? 
    * anyway, I just put the longer seq of two child as median and alignment for now.*)
    let new_med_chrom_arr = 
        let res = ref [] in
        for chrom_id = 0 to max_num_chrom - 1 do
            if debug_median then Printf.printf "create new med arr, work on chrom_id = %d\n%!" chrom_id;
            let chrom1_opt = get_chrom_from_arr chrom_id med1.chrom_arr in 
            let chrom2_opt = get_chrom_from_arr chrom_id med2.chrom_arr in
            match chrom1_opt, chrom2_opt with
            | Some chrom1, Some chrom2 ->
                if debug_median then Printf.printf "create med with chrom1 and chrom2\n%!";
                (*Here we are filling in seq from one of chromosome as median
                * seq. yes we can get the alied seq between genome1 and genome2, but
                * that's for concatenated chromsome, not for each pair of them.
                * cutting it does not make much sense.*)
                let chrom1,chrom2 =
                    if Sequence.length chrom1.seq >= Sequence.length chrom2.seq
                    then chrom1,chrom2
                    else chrom2,chrom1
                in
                let medseq = chrom1.seq in
                let sta = 0 and en = Sequence.length medseq in 
                let new_med_chromt = create_new_chromt sta en 0 chrom_id medseq
                med1.genome_ref_code med2.genome_ref_code 0 0 en (Sequence.length
                chrom2.seq) !(chrom1.chrom_id) !(chrom2.chrom_id) medseq medseq
                `Positive `Positive in
                res := new_med_chromt :: !res;
            | Some chrom1, None ->
                if debug_median then Printf.printf "create med with chrom1 and None\n%!";
                let medseq = chrom1.seq in
                let sta = 0 and en = Sequence.length medseq in
                let indel_cost = Sequence.cmp_gap_cost
                  ali_pam.ChromPam.chrom_indel_cost chrom1.seq
                in
                let new_med_chromt = create_new_chromt sta en indel_cost
                chrom_id medseq med1.genome_ref_code med2.genome_ref_code 0 (-1)
                en (-1) !(chrom1.chrom_id) (-1) medseq 
                (Sequence.create_gap_seq (Sequence.length medseq)) `Positive `Positive 
                in 
                res := new_med_chromt :: !res;
            | None, Some chrom2 ->
                if debug_median then Printf.printf "create med with None and chrom1\n%!";
                let medseq = chrom2.seq in
                let sta = 0 and en = Sequence.length medseq in
                let indel_cost = Sequence.cmp_gap_cost
                  ali_pam.ChromPam.chrom_indel_cost chrom2.seq
                in
                let new_med_chromt = create_new_chromt sta en indel_cost
                chrom_id medseq med1.genome_ref_code med2.genome_ref_code (-1) 0
                (-1) en (-1) !(chrom2.chrom_id) (Sequence.create_gap_seq
                (Sequence.length medseq)) medseq `Positive `Positive in
                res := new_med_chromt :: !res;
            | None, None -> failwith "genomeAli.create_med_mauve_annotator,both
            chrom1 and chrom2 are empty, impossible"
        done;
        List.rev !res
    in
    let total_cost = real_recost2+edit_and_indel_cost+ali_pam.ChromPam.translocation * translocation in
    if debug_median then Printf.printf "update total cost with recost(%d) + \
    edit&indel(%d),translocation number(%d)*trlcost(%d) = %d \n%!" real_recost2
    edit_and_indel_cost translocation ali_pam.ChromPam.translocation total_cost;
    let genome_med = { 
            chrom_arr = Array.of_list new_med_chrom_arr;
            genome_ref_code = Utl.get_new_genome_ref_code ();
            genome_ref_code1 = med1.genome_ref_code;
            genome_ref_code2 = med2.genome_ref_code;
            cost1 = total_cost - recost2;
            recost1 = recost1;
            cost2 = total_cost - recost1;
            recost2 = recost2;
            cost3 = 0; (* cost3 is 0 for median2, until we update it with median3 function*)
            sum_cost = total_cost + med1.sum_cost + med2.sum_cost;
        }
    in
    genome_med, total_cost,(recost1,recost2) 


let get_extra_cost_for_root x cost_mat =
    Array.fold_left (fun acc1 chromt ->
        let c = List.fold_left (fun acc2 segt ->
            let seq1 = segt.alied_seq1 
            and seq2 = segt.alied_seq2 in
            let newcost = Sequence.Align.cost_2 seq1 seq2 cost_mat Matrix.default in
            if debug_distance then Printf.printf "get_extra_cost_for_root, acc(%d)\
            += %d - %d\n%!" acc2 segt.cost newcost;
            segt.cost - newcost + acc2
        ) 0 chromt.map
        in
        acc1 + c
    ) 0 x.chrom_arr

(** [create_med med1 med2 cost_mat user_chrom_pams] 
* creates a median between two genomes [med1] and [med2] *)
let create_med med1 med2 cost_mat user_chrom_pams = 
    if ( Data.use_mauve_annotator user_chrom_pams ) then
        let ali_pam =  ChromPam.get_chrom_pam user_chrom_pams in
        create_med_mauve_annotator med1 med2 cost_mat ali_pam
    else 
        create_med_vinh_annotator med1 med2 cost_mat user_chrom_pams
          

(** [cmp_cost med1 med2 cost_mat user_chrom_pams] 
* returns the cost between genome [med1] and genome [med2] *)
let cmp_cost med1 med2 cost_mat user_chrom_pams = 
    if debug = true then begin
        Printf.printf "genomeAli.cmp_cost\n%!";
        let genomeFile = open_out "genome12" in 
        fprintf genomeFile ">genome1\n";  
        Array.iter (fun chrom -> 
                    Sequence.print genomeFile chrom.seq Alphabet.nucleotides;
                    fprintf genomeFile " @ ") med1.chrom_arr;
        fprintf genomeFile "\n";  
        fprintf genomeFile ">genome2\n";  
        Array.iter (fun chrom -> 
                    Sequence.print genomeFile chrom.seq Alphabet.nucleotides;
                    fprintf genomeFile " @ ") med2.chrom_arr;
        fprintf genomeFile "\n";  
        close_out genomeFile;
    end;   
    let ali_pam = ChromPam.get_chrom_pam user_chrom_pams in 
    match ali_pam.ChromPam.symmetric with
    | true ->
          let _, cost12, recost12 = create_med med1 med2 cost_mat user_chrom_pams in 
          let _, cost21, recost21 = create_med med2 med1 cost_mat user_chrom_pams in 
          if cost12 <= cost21 then cost12, recost12
          else cost21, recost21
    | false ->
          let med1, med2 = 
              match compare med1 med2 < 0 with 
              | true ->  med1, med2
              | false ->  med2, med1
          in 
          let med, cost, recost = create_med med1 med2 cost_mat user_chrom_pams in 
          cost, recost


(** [find_med2_ls med1 med2 cost_mat user_chrom_pams]
* return the median list between genome [med1] and genome [med2] *)
let find_med2_ls med1 med2 cost_mat user_chrom_pams = 
    if debug = true then begin
        let genomeFile = open_out "genome12" in 
        fprintf genomeFile ">genome1\n";  
        Array.iter (fun chrom -> 
                    Sequence.print genomeFile chrom.seq Alphabet.nucleotides;
                    fprintf genomeFile " @ ") med1.chrom_arr;
        fprintf genomeFile "\n";  

        fprintf genomeFile ">genome2\n";  
        Array.iter (fun chrom -> 
                    Sequence.print genomeFile chrom.seq Alphabet.nucleotides;
                    fprintf genomeFile " @ ") med2.chrom_arr;
        fprintf genomeFile "\n";  
        close_out genomeFile;
    end;   
    let ali_pam = ChromPam.get_chrom_pam user_chrom_pams in 
    let med, cost, recost = 
        match ali_pam.ChromPam.symmetric with
        | true ->
              let med12, cost12, recost12 = create_med med1 med2 cost_mat user_chrom_pams in 
              let med21, cost21, recost21 = create_med med2 med1 cost_mat user_chrom_pams in 
              if cost12 <= cost21 then med12, cost12, recost12
              else begin 
                  let med12 = swap_med med21 in 
                  med12, cost21, recost21
              end 

        | false ->
              let med1, med2, ali_pam, swaped = 
                  match compare med1 med2 < 0 with 
                  | true ->  med1, med2, ali_pam, false
                  | false ->  med2, med1, ali_pam, true
              in 
              
              let med, cost, recost = create_med med1 med2 cost_mat user_chrom_pams in 
              let med = 
                  match swaped with
                  | false -> med
                  | true -> swap_med med 
              in 
              med, cost, recost
    
    in                   
    cost, recost, [med]
    
(** [compare med1 med2] returns 0 if genome [med1] is
* identical to genome [med2], otherwise (-1) or 1 *)
let compare med1 med2 =
    let max_chrom = Array.length (!ref_genome).chrom_arr in 
    let chrom_ids = Array.make max_chrom 0 in 
    Array.iter (fun chrom -> 
                    let id = !(chrom.chrom_id) in 
                    chrom_ids.(id) <- chrom_ids.(id) + 1
               ) med1.chrom_arr;

    Array.iter (fun chrom -> 
                    let id = !(chrom.chrom_id) in 
                    chrom_ids.(id) <- chrom_ids.(id) + 1
               ) med2.chrom_arr;

    let rec is_compared chrom_id =
        if chrom_id = max_chrom then true  
        else if (chrom_ids.(chrom_id) = 0) || (chrom_ids.(chrom_id) = 2) then
            is_compared (chrom_id + 1)
        else false
    in 

    match is_compared 0 with
    | false -> compare (Array.length med1.chrom_arr) (Array.length med2.chrom_arr)
    | true -> 
          let rec compare_chrom chrom_id = 
              if chrom_id = max_chrom then 0
              else begin
                  let chrom1_opt = find_chrom chrom_id med1 in   
                  let chrom2_opt = find_chrom chrom_id med2 in   
                  match chrom1_opt, chrom2_opt with  
                  | Some chrom1, Some chrom2 ->  
                        let cmp = Sequence.compare chrom1.seq chrom2.seq in
                        if cmp = 0 then compare_chrom (chrom_id + 1) 
                        else cmp
                  | _, _ -> compare_chrom (chrom_id + 1) 
              end 
          in 
          compare_chrom 0
    




(** [create_map anc_med des_ref] returns the map
* from ancestor genome [anc_med] to descendant genome [des_ref] *)
let create_map anc_med des_ref : (int * int * Xml.xml) = 
    let str x = `Int x in  
    let create_chrom_map m =      
        let a_chrom_id, a_sta, a_en, a_dir  = 
            (str m.med_chrom_id), (str m.sta), (str m.en), `String "+" 
        in 
        let d_chrom_id, d_sta, d_en, d_dir = 
            match des_ref = anc_med.genome_ref_code1 with
            | true ->
                  (str m.chi1_chrom_id), (str m.sta1), (str m.en1), 
                  `String (Utl.get_dir m.dir1) 
            | false ->
                  (str m.chi2_chrom_id), (str m.sta2), (str m.en2), 
                    `String (Utl.get_dir m.dir2) 
        in 
        let attributes = [(Xml.GenomeMap.a_chrom_id, a_chrom_id);
                          (Xml.GenomeMap.a_start_seg, a_sta);
                          (Xml.GenomeMap.a_end_seg, a_en);
                          (Xml.GenomeMap.a_dir_seg, a_dir);
                          (Xml.GenomeMap.d_chrom_id, d_chrom_id);
                          (Xml.GenomeMap.d_start_seg, d_sta);
                          (Xml.GenomeMap.d_end_seg, d_en);
                          (Xml.GenomeMap.d_dir_seg, d_dir)
                         ] 
        in 
        let m : Xml.xml = (Xml.GenomeMap.seg, attributes, `String "") in 
        `Single m
    in   
    
    let chrom_map_ls = List.map 
        (fun chrom -> 
             let seg_ls = List.map create_chrom_map chrom.map in
             let chrom_map : Xml.xml = 
                 (Xml.GenomeMap.chrom, [], (`Set seg_ls))
             in 
             `Single chrom_map
        ) (Array.to_list anc_med.chrom_arr) 
    in 

    let attributes = [(Xml.GenomeMap.a_ref_code, str anc_med.genome_ref_code); 
                      (Xml.GenomeMap.d_ref_code, str des_ref)] 
    in 

    let genome_map : Xml.xml = 
        (Xml.GenomeMap.genome, attributes, (`Set chrom_map_ls)) 
    in 

    match des_ref = anc_med.genome_ref_code1 with
    | true -> anc_med.cost1, anc_med.recost1, genome_map
    | false -> anc_med.cost2, anc_med.recost2, genome_map


(** [create_single_map med] returns the map 
* of single states of genome [med] in Tag.Output format *)
let create_single_map med : Xml.xml = 
    let str x = `Int x in  
 
    let create_chrom_map m =      
        let a_chrom_id, a_sta, a_en, a_dir  = 
            (str m.chi2_chrom_id), (str m.sta2), (str m.en2), 
            `String (Utl.get_dir m.dir2)
        in 
        let d_chrom_id, d_sta, d_en, d_dir = 
            (str m.chi1_chrom_id), (str m.sta1), (str m.en1), 
            `String (Utl.get_dir m.dir1) 
        in 
        let attributes = [(Xml.GenomeMap.a_chrom_id, a_chrom_id);
                          (Xml.GenomeMap.a_start_seg, a_sta);
                          (Xml.GenomeMap.a_end_seg, a_en);
                          (Xml.GenomeMap.a_dir_seg, a_dir);
                          (Xml.GenomeMap.d_chrom_id, d_chrom_id);
                          (Xml.GenomeMap.d_start_seg, d_sta);
                          (Xml.GenomeMap.d_end_seg, d_en);
                          (Xml.GenomeMap.d_dir_seg, d_dir)
                         ] 
        in 
        let m : Xml.xml = (Xml.GenomeMap.seg, attributes, `String "") in 
        `Single m
    in   
    
    let chrom_map_ls = List.map 
        (fun chrom -> 
             let seg_ls = List.map create_chrom_map chrom.map in
             let chrom_map : Xml.xml = 
                 (Xml.GenomeMap.chrom, [], (`Set seg_ls))
             in 
             `Single chrom_map
        ) (Array.to_list med.chrom_arr) 
    in 

    let attributes = [(Xml.GenomeMap.a_ref_code, str med.genome_ref_code2); 
                      (Xml.GenomeMap.d_ref_code, str med.genome_ref_code1)] 
    in 

    let genome_map : Xml.xml = 
        (Xml.GenomeMap.genome, attributes, (`Set chrom_map_ls)) 
    in 

    genome_map 





(** [to_single single_parent med c2 pam] creates
* the single states from single state genome parent 
* [single_parent] to genome [med] *)
let to_single single_parent med c2 pam =
    if debug then begin
    Printf.printf "genomeAli.to_single, single parent = \n%!";
    print_genome single_parent;
    Printf.printf "med = \n%!";
    print_genome med;
    end;
    let is_first_child = single_parent.genome_ref_code1 = med.genome_ref_code in    
    let gap = Cost_matrix.Two_D.gap c2 in
    let single_genome = Array.map 
        (fun chromt -> 
             let seg_ls = Array.fold_left 
                 (fun seg_ls anc_chrom ->                  
                      List.fold_left 
                          (fun seg_ls seg ->                              
                               match is_first_child with
                               | true ->
                                     if seg.chi1_chrom_id != !(chromt.chrom_id) then seg_ls
                                     else begin                                         
                                         let single, _ = Sequence.closest_alied_seq
                                             seg.alied_med seg.alied_seq1 c2
                                         in  
                                         (single, seg.alied_seq1, seg.sta1, seg.en1)::seg_ls
                                     end   
                               | false  -> 
                                     if seg.chi2_chrom_id != !(chromt.chrom_id) then seg_ls
                                     else begin
                                         let single, _ = Sequence.closest_alied_seq
                                             seg.alied_med seg.alied_seq2 c2
                                         in  
                                         (single, seg.alied_seq2, seg.sta2, seg.en2)::seg_ls
                                     end 
                          ) seg_ls anc_chrom.map                           
                 ) [] single_parent.chrom_arr                  
             in   
             let sorted_seg_ls = List.sort 
                 (fun seg1 seg2 ->
                      let _,  _, sta1, _ = seg1 in 
                      let _,  _, sta2, _ = seg2 in
                      sta1 - sta2 
                 ) (List.filter (fun (_,  _, sta, _) -> sta != -1) seg_ls)
             in      
             let process_indel_locus sta en = 
                 match sta > en with
                 | true -> Sequence.get_empty_seq () 
                 | false ->                       
                       let ss = Sequence.sub chromt.seq sta (en - sta +1 ) in 
                       let single = Sequence.get_single_seq ss c2 in 
                       single
             in        
             let seq_len = Sequence.length chromt.seq in 
             let seg_ls, last_p =
                 List.fold_right 
                     (fun seg (seg_ls, last_p) ->
                          let alied_single, alied_seq, sta, en = seg in 
                          let gapless_alied_single = Sequence.fold_righti 
                              (fun gapless_alied_single p code ->
                                   match code = gap with
                                   | false -> code::gapless_alied_single
                                   | true ->
                                         if Sequence.get alied_seq p = gap then gapless_alied_single
                                         else code::gapless_alied_single
                              ) [] alied_single
                          in                       
                          let gapless_alied_single = Sequence.of_array (Array.of_list  gapless_alied_single) in
                          let indel_single = process_indel_locus (en + 1) last_p in 
                          (indel_single::gapless_alied_single::seg_ls), sta - 1
                     ) sorted_seg_ls ([], seq_len - 1) 
             in 
             let seg_ls = match last_p >= 0 with
             | false -> seg_ls
             | true ->
                   let indel_single = process_indel_locus 0 last_p in 
                   indel_single::seg_ls
             in 
            let single_chrom = Sequence.concat seg_ls in 
             single_chrom
        ) med.chrom_arr 
    in
    single_genome






(** [to_single_root root c2] create the single states
* for the genome at [root] *)
let to_single_root root c2 = 
    if debug then begin
        Printf.printf "genomeAli.to_single_root, root is\n%!";
        print_genome root;
    end;
    if (root.genome_ref_code1 = -1) && (root.genome_ref_code2 = -1) then 
        Array.map (fun ch -> ch.seq) root.chrom_arr
    else begin
        let gap = Cost_matrix.Two_D.gap c2 in
        let single_genome = Array.map 
            (fun chromt -> 
                 let map = 
                     List.map
                     (fun seg ->                              
                          let single_seg, _ = Sequence.closest_alied_seq
                              seg.alied_seq2 seg.alied_seq1 c2
                          in  
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
                     ) chromt.map  
                 in 
                 Sequence.concat map
            ) root.chrom_arr 
        in     
        single_genome
    end 

(** [change_to_single med single_genome] assigns
* the single states of genome [med] by [single_genome] *)
let change_to_single med single_genome = 
    if debug then begin
        Printf.printf "genomeAli.change_to_single,single_genome is\n%!";
    Array.iter(fun s -> Sequence.printseqcode s) single_genome;
    Printf.printf "med is \n%!";
    print_genome med;
    end;
    let gap = Alphabet.gap in 
    let new_chrom_arr = Array.mapi 
        (fun idx chromt -> 
             let single_seq = single_genome.(idx) in 
             let num_dna = ref 0 in  
             let new_map = match is_leaf med with 
             | true -> chromt.map
             | false -> 
                   List.map 
                 (fun seg ->
                     if debug then 
                         Printf.printf "genomeAli.change to single, idx=%d,single seq len=%d, alied_med seq len=%d\n%!"
                     idx (Sequence.length single_seq) (Sequence.length seg.alied_med);
                      let single_alied_med = Sequence.map 
                          (fun code ->
                               if code = gap then gap
                               else begin
                                   let single_code = Sequence.get single_seq !num_dna in 
                                   incr num_dna;
                                   single_code
                               end 
                          ) seg.alied_med
                      in         
                      {seg with alied_med = single_alied_med} 
                 ) chromt.map
             in 
             (if Sequence.length (Sequence.delete_gap single_seq) = 0 then 
                 failwith "The created single sequence is EMPTY");
             {chromt with map = new_map; seq = (Sequence.delete_gap single_seq)}
        ) med.chrom_arr 
    in 
    {med with chrom_arr = new_chrom_arr}


(** [copy_chrom_map s d] copys the chromosome map
* from genome source [s] to genome destination [d] *) 
let copy_chrom_map s d = 
    {d with genome_ref_code = s.genome_ref_code; 
         genome_ref_code1 = s.genome_ref_code1;
         genome_ref_code2 = s.genome_ref_code2;
         chrom_arr = s.chrom_arr}
        




(** [find_med3 ch1 ch2 ch3 mine cost_mat cost_cube pam]
* create the median sequence of [ch1], [ch2], [ch3] based
* on the current median [min] *)
let find_med3 ch1 ch2 ch3 mine cost_mat ori_cost_mat cost_cube pam = 
    let _, _, med1m_ls = find_med2_ls  ch1 mine cost_mat pam in
    let _, _, med2m_ls = find_med2_ls  ch2 mine cost_mat pam in
    let _, _, med3m_ls = find_med2_ls  ch3 mine cost_mat pam in

    let med1m = List.hd med1m_ls in
    let med2m = List.hd med2m_ls in
    let med3m = List.hd med3m_ls in

    let gap = Cost_matrix.Two_D.gap cost_mat in

    let num_mine_chrom = Array.length mine.chrom_arr in

    let get_chrom_id_map genome = 
        let max_chrom_id = Array.fold_left 
            (fun max_id chrom -> max max_id  !(chrom.chrom_id)) 0 genome.chrom_arr
        in 
        let map_id = Array.make (max_chrom_id + 1) 0 in 
        Array.iteri (fun idx chrom -> map_id.(!(chrom.chrom_id)) <- idx) genome.chrom_arr;
        map_id;
    in 

    let map1_id = get_chrom_id_map ch1 in 
    let map2_id = get_chrom_id_map ch2 in 
    let map3_id = get_chrom_id_map ch3 in 
    let mapm_id = get_chrom_id_map mine in 


    let create_pos map_id genome_med =
        let pos_mat = Array.init num_mine_chrom 
            (fun chi -> 
                 let li = Sequence.length mine.chrom_arr.(chi).seq in 
                 Array.init li (fun _ -> (-1, -2))
            )
        in 
        Array.iter 
            (fun chrom -> 
                 List.iter 
                     (fun seg ->
                          let c1 = map_id.(seg.chi1_chrom_id) in 
                          let c2 = mapm_id.(seg.chi2_chrom_id) in 
                          if (c1 >= 0) && (c2 >= 0)  && (seg.sta1 >= 0) && (seg.sta2 >= 0) then begin
                              let num1 = ref 0 in 
                              let num2 = ref 0 in 
                              let ali_len = Sequence.length seg.alied_med in 
                              let ali1 = seg.alied_seq1 in
                              let ali2 = seg.alied_seq2 in  

                              for p = 0 to ali_len -1 do
                                  let b1 = Sequence.get ali1 p in
                                  let b2 = Sequence.get ali2 p in 
                                  (if b1 != gap then incr num1);
                                  (if b2 != gap then begin
                                       incr num2;
                                       (pos_mat.(c2).(seg.sta2 + !num2 - 1) <- 
                                            if b1 = gap then (c1, -1)
                                            else (c1, seg.sta1 + !num1 - 1));
                                   end);
                                  
                              done;
                          end  
                     ) (List.filter (fun seg -> seg.dir2 = `Positive) chrom.map);
            ) genome_med.chrom_arr;
        pos_mat
    in 

    let pos1_mat = create_pos map1_id med1m in 
    let pos2_mat = create_pos map2_id med2m in 
    let pos3_mat = create_pos map3_id med3m in 

    let aliPam = ChromPam.get_chrom_pam pam in 
    let max_3d_len = min aliPam.ChromPam.detected_3d_len aliPam.ChromPam.max_3d_len in
    
    
    let rec detect_change chrom new_med f_p =
        let mine_len = Sequence.length chrom.seq  in 
        let c_id = !(chrom.chrom_id) in
        if f_p >= mine_len then new_med 
        else begin
            let rec find_first f_p = 
                if f_p = mine_len then (-1, -1), (-1, -1),(-1, -1), (-1, -1)
                else begin 
                    let c1, f_p1 = pos1_mat.(c_id).(f_p) 
                    and c2, f_p2 = pos2_mat.(c_id).(f_p)
                    and c3, f_p3 = pos3_mat.(c_id).(f_p)                         
                    in  
                    if (c1 >= 0 && f_p1 >=0) && (c2 >= 0 && f_p2 >=0) && 
                        (c3 >= 0 && f_p3 >=0) then (c_id, f_p), (c1, f_p1), (c2, f_p2), (c3, f_p3)
                    else find_first (f_p + 1)
                end
            in 

            let (_, f_p), (c1, f_p1), (c2, f_p2), (c3, f_p3) = find_first f_p in 
            let rec find_last l_p = 
                if l_p = f_p then (-1, -1), (-1, -1), (-1, -1), (-1, -1)
                else begin
                    let lc1, l_p1 = pos1_mat.(c_id).(l_p) 
                    and lc2, l_p2 = pos2_mat.(c_id).(l_p)
                    and lc3, l_p3 = pos3_mat.(c_id).(l_p)                         
                    in  

                    if (c1 = lc1) && (l_p1 >= 0 ) && (l_p1 > f_p1) && (l_p1 - f_p1 <= max_3d_len)
                        && (c2 = lc2) && (l_p2 >= 0 ) && (l_p2 > f_p2) && (l_p2 - f_p2 <= max_3d_len)
                        && (c3 = lc3) && (l_p3 >= 0 ) && (l_p3 > f_p3) && (l_p3 - f_p3 <= max_3d_len) 
                    then (c_id, l_p), (c1, l_p1), (c2, l_p2), (c3, l_p3)
                    else find_last (l_p - 1)
               end 
            in
            
            if f_p = -1 then new_med 
            else begin
                let max_l_p = min (f_p + max_3d_len) (mine_len - 1)  in
                let (_, l_p), (_, l_p1), (_, l_p2), (_, l_p3) = find_last max_l_p in 
                if l_p = -1 then 
                    detect_change chrom new_med (max_l_p + 1)
                else begin
                    let sub_seq1 = Sequence.sub ch1.chrom_arr.(c1).seq f_p1 (l_p1 - f_p1 + 1) in
                    let sub_seq2 = Sequence.sub ch2.chrom_arr.(c2).seq f_p2 (l_p2 - f_p2 + 1) in
                    let sub_seq3 = Sequence.sub ch3.chrom_arr.(c3).seq f_p3 (l_p3 - f_p3 + 1) in
                    let sub_seqm = Sequence.sub mine.chrom_arr.(c_id).seq f_p (l_p - f_p + 1) in



                    if (Sequence.length sub_seq1 = 0) ||  (Sequence.length sub_seq2 = 0) ||
                        (Sequence.length sub_seq3 = 0) then detect_change chrom new_med (l_p + 1)
                    else begin                        

                        let _, _, median_seq, _, _ , _  = Sequence.Align.readjust_3d
                            ~first_gap:false sub_seq1 sub_seq2 sub_seqm cost_mat
                            ori_cost_mat cost_cube sub_seq3 
                            mine.cost2 mine.cost3
                        in

                        detect_change chrom ((c_id, f_p, l_p, median_seq)::new_med) (l_p + 1)
                    end  
                end                                     
            end 
        end 
    in 
    let change_chrom chrom = 
        let new_med_rev = detect_change chrom [] 0 in         
        let l_p, seq_ls = List.fold_left 
            (fun (l_p, seq_ls) (c_id, s, e, med) ->
                 if s > l_p then begin
                     let sub_seq = Sequence.sub chrom.seq l_p (s - l_p) in 
                     let seq_ls = med::(sub_seq::seq_ls) in 
                     (e + 1), seq_ls
                 end else (e + 1), med::seq_ls 
            ) (0, []) (List.rev new_med_rev) 
        in

        let mine_len = Sequence.length chrom.seq  in 
        let seq_ls = 
            if l_p < mine_len then begin
                let sub_seq = Sequence.sub chrom.seq l_p (mine_len - l_p) in  
                sub_seq::seq_ls 
            end else seq_ls 
        in  
    
        let new_seq = Sequence.concat (List.rev seq_ls) in
        let new_chrom = {chrom with seq = new_seq} in
        new_chrom
    in
(* no need to calculate cost3 again for we carry it with us now
    let cost1, _ = cmp_cost ch1 mine cost_mat pam in
    let cost2, _ = cmp_cost ch2 mine cost_mat pam in
    let cost3, _ = cmp_cost ch3 mine cost_mat pam in
    let total_cost = cost1 + cost2 + cost3 in
*)
    (*get old cost12 and cost3*)
    let old_cost3 = mine.cost3 in
    let old_cost12 = mine.cost1 + mine.cost2 in
    let old_sumcost = mine.sum_cost in
    (*get new cost12 and cost3*)
    let new_chrom_arr = Array.map change_chrom mine.chrom_arr in 
    let new_mine = {mine with chrom_arr = new_chrom_arr} in 
    let new_cost_mine_ch1, (rc1,rc2) = cmp_cost ch1 new_mine cost_mat pam in
    let new_recost_mine_ch1 = rc1 + rc2 in
    let new_cost_mine_ch2, (rc1,rc2) = cmp_cost ch2 new_mine cost_mat pam in
    let new_recost_mine_ch2 = rc1 + rc2 in
    let new_cost12 = new_cost_mine_ch1 + new_cost_mine_ch2 in
    let new_cost_mine_parent, _ = cmp_cost ch3 new_mine cost_mat pam in
    let new_cost3 = new_cost_mine_ch1 + new_cost_mine_ch2 + new_cost_mine_parent in
    let new_sumcost = ch1.sum_cost + ch2.sum_cost + new_cost12 in
    (*update new_mine with cost12 and cost3*)
    let new_mine = {new_mine with
            cost1 = new_cost_mine_ch1;
            recost1 = new_recost_mine_ch1;
            cost2 = new_cost_mine_ch2;
            recost2 = new_recost_mine_ch2;
            cost3 = new_cost3; 
            sum_cost = new_sumcost;
    } in
    (*return new median with cost12, cost3*)
    if (new_cost12 < old_cost12)||(new_cost3 < old_cost3)||(new_sumcost <
    old_sumcost)(*we don't compare the median sequence, for they are array of
    sequences, and they might be huge, ||(0 <> compare new_chrom_arr chrom_arr)*) 
    then 
            new_cost12, new_cost3, new_mine, true
    else 
        old_cost12, old_cost3, mine, false
