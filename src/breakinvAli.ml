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

let () = SadmanOutput.register "BreakinvAli" "$Revision: 911 $"

(** The implementation of funtions to calculate the cost, 
* alignments and medians between general sequences 
* where both point mutations and rearrangement operations
* are considered *)

let debug = false
let fprintf = Printf.fprintf

(** Data structure to contain a breakinv sequence *)
type breakinv_t = {
    seq : Sequence.s; (** the breakinv sequence *)

    alied_med : Sequence.s; (** the aligned median sequence  *)
    alied_seq1 : Sequence.s; (** the aligned sequence of its first child *)
    alied_seq2 : Sequence.s; (** the aligned sequence of its second child *)

    ref_code : int; (** the reference code of this breakinv *)
    ref_code1 : int; (** the reference code of its first child *)
    ref_code2 : int; (** the reference code of its second child *)
    cost1 : int; (** the cost between this breakinv and its first child *)
    cost2 : int; (** the cost between this breakinv and its second child *)
    recost1 : int; (** the recost between this breakinv and its first child *)
    recost2 : int; (** the recost between this breakinv and its second child *)

    delimiter_lst : int list (* delimiter list for multichromosome *)
}


(** Data structure to contain parameters 
* used to align two breakinv sequences *)
type breakinvPam_t = {
    median_solver: Data.median_solver_t; 
    re_meth : Data.re_meth_t;
    keep_median : int;
    circular : int;
    swap_med : int;
    symmetric : bool;
    kept_wag : int;
}

(** [swap_med m] swaps the first child and second child
* of breakinv sequence [m] *)
let swap_med m = {
    m with alied_seq1 = m.alied_seq2;
        alied_seq2 = m.alied_seq1;
        ref_code1 = m.ref_code2;
        ref_code2 = m.ref_code1;
        cost1 = m.cost2;
        cost2 = m.cost1;
        recost1 = m.recost2;
        recost2 = m.recost1        
}


let breakinvPam_default = {
    median_solver = `Albert;
    re_meth = `Locus_Breakpoint 10;
    keep_median = 1;
    circular = 0;
    swap_med = 1;
    symmetric = true;   
    kept_wag = 1;
}


(** [init seq] creates a new breakinv sequence from [seq] *)
let init seq =  {        
    seq = seq;  
    alied_med = Sequence.get_empty_seq ();
    alied_seq1 = Sequence.get_empty_seq ();
    alied_seq2 = Sequence.get_empty_seq ();

    ref_code = Utl.get_new_chrom_ref_code ();
    ref_code1 = -1;
    ref_code2 = -1;
    cost1 = 0;
    cost2 = 0;
    recost1 = 0;
    recost2 = 0;    

    delimiter_lst = [];
}

let equal_orientation code1 code2 = compare (abs code1) (abs code2) 

let get_common_seq (med1:breakinv_t) (med2:breakinv_t) (med3:breakinv_t) = 
    let arr1 = Sequence.to_array med1.seq
    and arr2 = Sequence.to_array med2.seq
    and arr3 = Sequence.to_array med3.seq in
    let res = Utl.get_common3 arr1 arr2 arr3 equal_orientation in
    res

(** [get_breakinv_pam user_breakinv_pam] returns 
* user defined parameters used to align two breakinv sequences *)
let get_breakinv_pam user_breakinv_pam = 
    let chrom_pam = breakinvPam_default in  
    let chrom_pam =
         match user_breakinv_pam.Data.median_solver with
        | None -> chrom_pam
        | Some median_solver -> {chrom_pam with median_solver = median_solver}
    in 
    let chrom_pam = 
        match user_breakinv_pam.Data.re_meth with
        | None -> chrom_pam
        | Some re_meth -> {chrom_pam with re_meth = re_meth}
    in 
    let chrom_pam = 
        match user_breakinv_pam.Data.circular with
        | None -> chrom_pam
        | Some circular -> {chrom_pam with circular = circular}
    in 
    let chrom_pam = 
        match user_breakinv_pam.Data.keep_median with
        | None -> chrom_pam
        | Some keep_median -> {chrom_pam with keep_median = keep_median}
    in 
    let chrom_pam = 
        match user_breakinv_pam.Data.swap_med with
        | None -> chrom_pam
        | Some swap_med -> {chrom_pam with swap_med = swap_med}
    in 
    let chrom_pam = 
        match user_breakinv_pam.Data.symmetric with
        | None -> chrom_pam
        | Some sym -> {chrom_pam with symmetric = sym}
    in 
    let chrom_pam = 
        match user_breakinv_pam.Data.max_kept_wag with
        | None -> chrom_pam
        | Some max_kept_wag -> {chrom_pam with kept_wag = max_kept_wag}
    in
    chrom_pam

        
let get_recost user_pams = 
    match user_pams.Data.re_meth with
    | None -> failwith "The rearrangement cost is not specified"
    | Some re_meth ->
        match re_meth with
            | `Locus_Breakpoint c -> c
            | `Locus_Inversion c -> c



(** [cmp_cost med1 med2 gen_cost_mat pure_gen_cost_mat alpha breakinv_pam]
* computes total cost between two breakinv sequences [med1] and [med2].
* the total cost = editing cost + rearrangement cost *)
let cmp_cost med1 med2 gen_cost_mat pure_gen_cost_mat alpha breakinv_pam =
(*debug msg 
    Printf.printf "cmp_cost: %!"; 
    Sequence.printseqcode med1.seq; Sequence.printseqcode med2.seq;
 debug msg*)
    let ali_pam = get_breakinv_pam breakinv_pam in     
    let len1 = Sequence.length med1.seq in 
    let len2 = Sequence.length med2.seq in 
    let orientation = Alphabet.get_orientation alpha in
    if (len1 < 1) || (len2 < 1) then 0, (0, 0)
    else begin
        match ali_pam.symmetric with
        | true ->
                let cost12, (recost12a,recost12b), _, _ =
                  GenAli.create_gen_ali ali_pam.kept_wag `Breakinv med1.seq med2.seq gen_cost_mat 
                      alpha ali_pam.re_meth ali_pam.swap_med ali_pam.circular orientation
              in 
              let cost21, (recost21a,recost21b), _, _ = 
                  GenAli.create_gen_ali ali_pam.kept_wag `Breakinv med2.seq med1.seq gen_cost_mat 
                      alpha ali_pam.re_meth ali_pam.swap_med ali_pam.circular orientation
              in 
              if cost12 <= cost21 then cost12, (recost12a,recost12b)
              else cost21, (recost21a,recost21b)
        | false ->
              let cost, recost, _, _ = 
                  if Sequence.compare med1.seq med2.seq < 0 then                       
                      GenAli.create_gen_ali  ali_pam.kept_wag `Breakinv med1.seq med2.seq gen_cost_mat 
                          alpha ali_pam.re_meth ali_pam.swap_med ali_pam.circular orientation 
                  else 
                      GenAli.create_gen_ali  ali_pam.kept_wag `Breakinv med2.seq med1.seq gen_cost_mat 
                          alpha ali_pam.re_meth ali_pam.swap_med ali_pam.circular orientation
              in      
              cost , recost
    end 

let better_capping med1_seq med2_seq delimiters =
    let arr1 = Sequence.to_array med1_seq 
    and arr2 = Sequence.to_array med2_seq in
    assert( (Array.length arr1) = (Array.length arr2) );
    let deli1 = Array.of_list delimiters in 
    let deli2 = deli1 in 
    (*this is true since we choose delimiters of med1 as
    the delimiter for med2. btw, med1 is one of the two 
    parents of med2 in "pick_delimiters" *)
    let new_medseq2, new_deli2 = 
        UtlGrappa.find_better_capping arr1 arr2 deli1 deli2
    in
    (Sequence.of_array new_medseq2),(Array.to_list new_deli2)


let pick_delimiters med1 med2 med_seq = 
    let arr1 = Sequence.to_array med1.seq 
    and arr2 = Sequence.to_array med2.seq
    and arr3 = Sequence.to_array med_seq
    in
    assert( (Array.length arr1) = (Array.length arr2) );
    assert( (Array.length arr3) = (Array.length arr2) );
    let num_genes = Array.length arr3 in
    let deli1 = Array.of_list med1.delimiter_lst 
    and deli2 = Array.of_list med2.delimiter_lst
    in
    (* debug msg 
    let print_intarr arr = 
        Printf.printf "[%!";
        Array.iter (Printf.printf "%d,%!") arr;
        Printf.printf "],%!";
    in
    Printf.printf "pick delimiters :\n { %!";
    Sequence.printseqcode med1.seq;
    Sequence.printseqcode med2.seq;
    Sequence.printseqcode med_seq;
    print_intarr deli1; print_intarr deli2; 
    Printf.printf " } num_genes = %d \n%!" num_genes;
     debug msg *)
    let dis11 = UtlGrappa.cmp_inversion_dis_multichrom 
                arr1 arr3 deli1 deli1 num_genes
    in
    let dis12 = UtlGrappa.cmp_inversion_dis_multichrom 
                arr2 arr3 deli2 deli1 num_genes
    in
    let dist21 = UtlGrappa.cmp_inversion_dis_multichrom 
               arr1 arr3 deli1 deli2 num_genes  
    in
    let dist22 = UtlGrappa.cmp_inversion_dis_multichrom 
               arr2 arr3 deli2 deli2 num_genes
    in
   (*debug msg
    if (dis11+dis12)>(dist21+dist22) then Printf.printf "pick delimiter 1 \n%!"
    else Printf.printf "pick delimiter 2\n%!";
    debug msg*)
    if (dis11+dis12)>(dist21+dist22) then med2,med2.delimiter_lst
    else med1,med1.delimiter_lst


(** find_simple_med2_ls med1 med2 gen_cost_mat pure_gen_cost_mat alpha ali_pam]
* finds all medians between breakinv sequence [med1] and [med2]
* allowing rearrangements *) 
let find_simple_med2_ls med1 med2 gen_cost_mat pure_gen_cost_mat alpha ali_pam = 
    (* debug msg 
    let print_intlist lst = 
        Printf.printf "[%!";
        List.iter (Printf.printf "%d,") lst;
        Printf.printf "]\n%!"
    in
     debug msg *)
    (* debug msg 
    Printf.printf "find_simple_med2_ls in breakinvAli.ml :\n";
    Sequence.printseqcode med1.seq; Sequence.printseqcode med2.seq;
     debug msg *)
    let len1 = Sequence.length med1.seq in 
    let len2 = Sequence.length med2.seq in
    let orientation = Alphabet.get_orientation alpha in
    if len1 < 1 then 0, (0, 0), [med2]
    else if len2 < 1 then 0, (0, 0), [med1] 
    else begin        
        let total_cost, (recost1, recost2), alied_gen_seq1, alied_gen_seq2 = 
            GenAli.create_gen_ali  ali_pam.kept_wag `Breakinv med1.seq med2.seq gen_cost_mat 
                alpha ali_pam.re_meth ali_pam.swap_med ali_pam.circular orientation 
        in 
        let re_seq2 =
            Utl.filterArr (Sequence.to_array alied_gen_seq2) (fun code2 -> code2 != Alphabet.get_gap alpha)
        in     
        let all_order_ls =  [(re_seq2, recost1, recost2)]  in 
        let med_ls = List.fold_left 
            (fun med_ls (re_seq2, recost1, recost2) ->
                 let med_seq = Sequence.init 
                     (fun index ->
                          if re_seq2.(index) > 0 then 
                              begin
                                  Printf.printf "";
                                  re_seq2.(index)
                              end
                          else 
                              if re_seq2.(index) mod 2  = 0 then -re_seq2.(index) - 1
                              else -re_seq2.(index) + 1
                     ) (Array.length re_seq2)
                 in
                 (* create delimiter for the new median *)
                 (* just pick the better one from its two parents, for now *)
                 (* then we do "find better capping" after this*)
                 let chosen_parent,newdelimiters = pick_delimiters med1 med2 med_seq in
                 let med_seq, newdelimiters = 
                     better_capping chosen_parent.seq med_seq newdelimiters in
                 (* debug msg 
                 Printf.printf "seqcode with better capping:\n%!";
                 Sequence.printseqcode med_seq; print_intlist newdelimiters; 
                 debug msg *)
                 let newrefcode =  Utl.get_new_chrom_ref_code () in
                 let med = 
                     {seq = med_seq; 
                      alied_med = alied_gen_seq2;
                      alied_seq1 = alied_gen_seq1;
                      alied_seq2 = alied_gen_seq2;
                      ref_code = newrefcode;
                      ref_code1 = med1.ref_code;
                      ref_code2 = med2.ref_code;
                      cost1 = total_cost - recost2;
                      cost2 = total_cost - recost1;
                      recost1 = recost1;
                      recost2 = recost2;
                      delimiter_lst = newdelimiters 
                     }
                 in    
               med::med_ls
            ) [] all_order_ls 
        in
        total_cost, (recost1, recost2), med_ls
    end



(** [find_med2_ls med1 med2 gen_cost_mat pure_gen_cost_mat alpha breakinv_pam]
* finds all median sequences between [med1] and [med2] and vice versa.
* Rearrangements are allowed *) 
let find_med2_ls med1 med2 gen_cost_mat pure_gen_cost_mat alpha breakinv_pam = 
    if med1 == med2 then
        begin
        0, (0,0), [
        { med1 with 
        alied_seq1 = med1.seq;
        alied_seq2 = med2.seq; 
        ref_code = Utl.get_new_chrom_ref_code ();
        ref_code1 = med1.ref_code;
        ref_code2 = med2.ref_code;
        cost1 = 0;
        cost2 = 0;
        recost1 = 0;
        recost2 = 0;
        delimiter_lst = med1.delimiter_lst
        }]
        end
    else
    let ali_pam = get_breakinv_pam breakinv_pam in          
    match ali_pam.symmetric with 
    | true ->
          let cost12, recost12, med12_ls = find_simple_med2_ls med1 med2
              gen_cost_mat pure_gen_cost_mat alpha ali_pam in
          let cost21, recost21, med21_ls = find_simple_med2_ls med2 med1
              gen_cost_mat pure_gen_cost_mat alpha ali_pam in
          if cost12 <= cost21 then 
              begin
                  cost12, recost12, med12_ls
              end
          else begin 
              let med12_ls = List.map swap_med med21_ls in 
              cost21, recost21, med12_ls
          end 
    | false ->
          let med1, med2, swaped = 
              match Sequence.compare med1.seq med2.seq < 0 with 
              | true ->  med1, med2, false
              | false -> med2, med1, true
          in 
          let cost, recost, med_ls = find_simple_med2_ls med1 med2 
              gen_cost_mat pure_gen_cost_mat alpha ali_pam in
          let med_ls = 
              match swaped with
              | false -> med_ls
              | true -> List.map swap_med med_ls
          in
    (* debug msg 
    Sequence.printseqcode med1.seq; Sequence.printseqcode med2.seq;
    Sequence.printseqcode (List.hd med_ls).seq;
    debug msg*)
    cost, recost, med_ls

(** [get_costs med child_ref] returns the cost
* from this breakinv median to its child [child_ref] *)
let get_costs med child_ref = 
    match child_ref = med.ref_code1 with
    | true -> med.cost1, med.recost1
    | false -> med.cost2, med.recost2


let create_breakinv_test =
    if debug then
        let char_arr = [|"A"; "B"; "C"; "D"; "E"; "F"; "G"; "H"; "I"; "J"; "K"; "N"; "M";
                         "L"; "P"; "Q"; "Z"; "X"; "Y"; "R"; "S"; "V"; "O"; "W"; "T";
                         "U"|] in 
        let num_char = 27 in
        let gap_cost = 10 in 
        let cost_mat = Array.create_matrix num_char num_char gap_cost in 
        cost_mat.(num_char - 1).(num_char - 1) <- 0;
        for i = 0 to num_char - 2 do
            for j = 0 to num_char - 2 do
                cost_mat.(i).(j) <-
                    if i = j then 0
                    else (Random.int (gap_cost - 1)) + 1 
            done;
        done;
        Random.self_init ();
        let num_tax = 22 in 
        let seqfile = open_out "taxa.bk" in 
        for t = 0 to num_tax - 1 do
            fprintf seqfile ">T%i\n" t;
            let mask_arr = Array.make num_char true in 
            for c = 0 to num_char - 2 do 
                let ch = Random.int (num_char - 1) in 
                if mask_arr.(ch) = true then begin
                    (if Random.int 2 = 0 then fprintf seqfile "%s " char_arr.(ch)
                    else fprintf seqfile "%s " char_arr.(ch));
                    mask_arr.(ch) <- false;
                end 
            done;
            fprintf seqfile "\n";
        done;
        close_out seqfile;
        
        let matfile = open_out "matrix.bk" in 
        Array.iter (fun s -> fprintf matfile "%s  " s) char_arr;
        fprintf matfile "\n"; 
        for i = 0 to num_char - 1 do
            for j = 0 to num_char - 1 do
                fprintf matfile "%i  " cost_mat.(i).(j);
            done;
            fprintf matfile "\n";
        done;
        close_out matfile;
        print_endline "End of create_breakinv_test"
    else ()

let update_bkinv_t old_bkinvt newseq delimiters =
    { old_bkinvt with seq = newseq; delimiter_lst = delimiters }

let single_to_multi bkinvt =
    let seq = bkinvt.seq and delimiters = bkinvt.delimiter_lst in
    assert( (List.length delimiters)>0 );
    let count = ref 0 in
    Printf.printf "single to multi in breakinvAli.ml, check delimiter list:\n [%!";
    List.iter (Printf.printf "%d,") delimiters;
    Printf.printf "]\n%!";
    let seq_lst = List.map ( fun delimiter ->
        let empty: Sequence.s = Sequence.get_empty_seq () in
        for i = !count to (!count + delimiter)-1 do
            Sequence.prepend empty (Sequence.get seq i)
        done;
        count := !count + delimiter;
        (Sequence.reverse empty)
    ) delimiters in
    Printf.printf "break single seq into seq list:\n%!";
    List.map (fun s ->
        Sequence.printseqcode s;
        { bkinvt with seq = s; delimiter_lst = [] }
    ) seq_lst

