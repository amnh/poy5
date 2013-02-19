(* POY 5.0 Beta. A phylogenetic analysis program using Dynamic Homologies.    *)
(* Copyright (C) 2013 Andrés Varón, Lin Hong, Nicholas Lucaroni, Ward Wheeler,*)
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

let () = SadmanOutput.register "BreakinvAli" "$Revision: 3160 $"

(** The implementation of funtions to calculate the cost, 
* alignments and medians between general sequences 
* where both point mutations and rearrangement operations
* are considered *)

let debug = false
let debug_distance = false
let fprintf = Printf.fprintf
let error_user_message format = Printf.ksprintf (Status.user_message Status.Error) format
let info_user_message format = Printf.ksprintf (Status.user_message Status.Information) format

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

    (*store cost3 and sumcost in each breakinv_t*)
    cost3 : int; (** cost(mine,ch1)+cost(mine,ch2)+cost(mine,parent) *)
    sum_cost : int; (** cost of subtree root at this breakinv_t of this node*)

    delimiter_lst : int list (* delimiter list for multichromosome *)
}


let print bkt =
    Printf.printf "[ seq=%!";
    Sequence.printseqcode bkt.seq;
    Printf.printf "alied_med=%!";
    Sequence.printseqcode bkt.alied_med;
    Printf.printf "alied_seq1=%!";
    Sequence.printseqcode bkt.alied_seq1;
    Printf.printf "alied_seq2=%!";
    Sequence.printseqcode bkt.alied_seq2;
    Printf.printf "cost1=%d,cost2=%d,recost1=%d,recost2=%d\n%!"
    bkt.cost1 bkt.cost2 bkt.recost1 bkt.recost2;
    Printf.printf "]\n%!"




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
let init seq delimiters =  
    {        
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

    cost3 = 0;
    sum_cost = 0;

    delimiter_lst = delimiters;
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

let get_extra_cost_for_root x cost_mat pure_gen_cost_mat alpha breakinv_pam = 
    let ali_pam = get_breakinv_pam breakinv_pam in 
    let gapcode = Alphabet.get_gap alpha in
    let seq1 = x.alied_seq1 and seq2 = x.alied_seq2 in
    let seq1 = Sequence.remove_gaps seq1 ~gapcode:gapcode false in
    let seq2 = Sequence.remove_gaps seq2 ~gapcode:gapcode false in
    let len1 = Sequence.length seq1 
    and len2 = Sequence.length seq2 in
    let orientation = Alphabet.get_orientation alpha in
    let newcost = 
        if (len1 < 1) || (len2 < 1) then 0
        else begin
            match ali_pam.symmetric with
            | true ->
                   let cost12, (recost12a,recost12b), _, _ =
                      GenAli.create_gen_ali ali_pam.kept_wag `Breakinv seq1
                      seq2 pure_gen_cost_mat alpha ali_pam.re_meth 
                      ali_pam.swap_med ali_pam.circular orientation
                  in 
                  let cost21, (recost21a,recost21b), _, _ = 
                      GenAli.create_gen_ali ali_pam.kept_wag `Breakinv seq2
                      seq1 pure_gen_cost_mat alpha ali_pam.re_meth 
                      ali_pam.swap_med ali_pam.circular orientation
                  in 
                  if cost12 <= cost21 then cost12
                  else cost21
            | false ->
                  let cost, _ , _, _ = 
                      if Sequence.compare x.alied_seq1 x.alied_seq2 < 0 then                       
                          GenAli.create_gen_ali  ali_pam.kept_wag `Breakinv seq1
                          seq2 pure_gen_cost_mat alpha ali_pam.re_meth 
                          ali_pam.swap_med ali_pam.circular orientation 
                      else 
                          GenAli.create_gen_ali  ali_pam.kept_wag `Breakinv seq2
                          seq1 pure_gen_cost_mat alpha ali_pam.re_meth 
                          ali_pam.swap_med ali_pam.circular orientation
                  in      
                  cost 
        end 
    in
    (*let newcost = Sequence.Align.cost_2 seq1 seq2 cost_mat Matrix.default in
* *)
    let oldcost = x.cost1 + x.recost1 in 
    if debug_distance then Printf.printf "get_extra_cost_for_root = %d - %d\n%!"  oldcost newcost;
    oldcost - newcost 



(** [cmp_cost med1 med2 gen_cost_mat pure_gen_cost_mat alpha breakinv_pam]
* computes total cost between two breakinv sequences [med1] and [med2].
* the total cost = editing cost + rearrangement cost *)
let cmp_cost med1 med2 gen_cost_mat pure_gen_cost_mat alpha breakinv_pam =
    let debug = false in
    let len1wg = Sequence.length med1.seq 
    and len2wg = Sequence.length med2.seq in
    (*remove gaps if there is any from medseq*)
    let gapcode = Alphabet.get_gap alpha in
    let seq1 = Sequence.remove_gaps med1.seq ~gapcode:gapcode false in
    let seq2 = Sequence.remove_gaps med2.seq ~gapcode:gapcode false in
    if debug then begin
    Printf.printf "breakinvAli.cmp_cost,gapcod=%d: \n%!" gapcode; 
    Sequence.printseqcode med1.seq; Sequence.printseqcode med2.seq;
    end;
    let ali_pam = get_breakinv_pam breakinv_pam in     
    let len1 = Sequence.length seq1 in 
    let len2 = Sequence.length seq2 in 
    let orientation = Alphabet.get_orientation alpha in
    if len1<>len2 then
        failwith ("cmp_cost funciton in breakinvAli.ml.ERROR:input sequence without gap have different length");
    let gaplen1 = len1wg - len1 
    and gaplen2 = len2wg - len2 in
    if debug && (gaplen1>0||gaplen2>0) then
        info_user_message "cmp_cost funciton in breakinvAli.ml.\
        Warning:we have gaps in input, ignore them";
    if (len1 < 1) || (len2 < 1) then 0, (0, 0)
    else begin
        match ali_pam.symmetric with
        | true ->
               let cost12, (recost12a,recost12b), _, _ =
                  GenAli.create_gen_ali ali_pam.kept_wag `Breakinv seq1
                  seq2 pure_gen_cost_mat alpha ali_pam.re_meth 
                  ali_pam.swap_med ali_pam.circular orientation
              in 
              let cost21, (recost21a,recost21b), _, _ = 
                  GenAli.create_gen_ali ali_pam.kept_wag `Breakinv seq2
                  seq1 pure_gen_cost_mat alpha ali_pam.re_meth 
                  ali_pam.swap_med ali_pam.circular orientation
              in 
              if cost12 <= cost21 then cost12, (recost12a,recost12b)
              else cost21, (recost21a,recost21b)
        | false ->
              let cost, recost, _, _ = 
                  if Sequence.compare med1.seq med2.seq < 0 then                       
                      GenAli.create_gen_ali  ali_pam.kept_wag `Breakinv seq1
                      seq2 pure_gen_cost_mat alpha ali_pam.re_meth 
                      ali_pam.swap_med ali_pam.circular orientation 
                  else 
                      GenAli.create_gen_ali  ali_pam.kept_wag `Breakinv seq2
                      seq1 pure_gen_cost_mat alpha ali_pam.re_meth 
                      ali_pam.swap_med ali_pam.circular orientation
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
    Printf.printf "find_simple_med2_ls in breakinvAli.ml,med1/med2.seq = :\n";
    Sequence.printseqcode med1.seq; Sequence.printseqcode med2.seq;
    Printf.printf " delimiter list = \n%!";
    print_intlist med1.delimiter_lst; print_intlist med2.delimiter_lst; 
     debug msg *)
    let len1wg = Sequence.length med1.seq 
    and len2wg = Sequence.length med2.seq in
    if (len1wg<>len2wg) then
        failwith ("median funciton in breakinvAli.ml.ERROR:input sequence have different length");
    (*I add this for multipul data file input.
    * file1 has t1 and t2
        * >t1
        * abc
        * >t2
        * bca
    * file2 has only t3 and t4
        * >t3
        * cab
        * >t4
        * bac
    * then we have for t1 and t2, med1.seq=[abc][-] and med2.seq=[bca][-] 
    How about this?
    * file1 
        * >t1 
        * abc
        * >t2
        * a
    * file2:
        * >t1
        * there is no t1 in file2
        * >t2
        * bc
    * here we end up having med1.seq=[abc][-] and med2.seq=[a][bc], is this
    * allowed?*)
    let gapcode = Alphabet.get_gap alpha in
    let seq1 = Sequence.remove_gaps med1.seq ~gapcode:gapcode false in
    let seq2 = Sequence.remove_gaps med2.seq ~gapcode:gapcode false in
    let len1 = Sequence.length seq1 in 
    let len2 = Sequence.length seq2 in
    if len1<>len2 then
        failwith ("median funciton in breakinvAli.ml.ERROR:input sequence without gap have different length");
    let gaplen = len1wg - len1 in
    if gaplen>0 then
        info_user_message "median funciton in breakinvAli.ml.\
        Warning:we have gaps in input, median function will ignore the gaps "; 
    let orientation = Alphabet.get_orientation alpha in
    if len1 < 1 then 0, (0, 0), [med2]
    else if len2 < 1 then 0, (0, 0), [med1] 
    else begin        
        let total_cost, (recost1, recost2), alied_gen_seq1, alied_gen_seq2 = 
            GenAli.create_gen_ali  ali_pam.kept_wag `Breakinv seq1 seq2
            pure_gen_cost_mat alpha ali_pam.re_meth ali_pam.swap_med 
            ali_pam.circular orientation 
        in
        let re_seq2 =
            Array_ops.filter (fun code2 -> code2 != Alphabet.get_gap alpha) (Sequence.to_array alied_gen_seq2)
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
                 let med_seq = 
                     if gaplen>0 then
                         let gapseq = Sequence.create_gap_seq ~gap:gapcode gaplen in
                         Sequence.concat [med_seq;gapseq]
                     else med_seq
                 in
                 (* create delimiter for the new median *)
                 (* just pick the better one from its two parents, for now *)
                 (* then we do "find better capping" after this*)
                 let len1 = List.length med1.delimiter_lst 
                 and len2 =  List.length med2.delimiter_lst in
                 let med_seq, newdelimiters =
                     if gaplen>0 then 
                         med_seq,med1.delimiter_lst
                     else 
                     match med1.delimiter_lst,med2.delimiter_lst with
                     | h1::t1, h2::t2 -> 
                         if (len1=1)&&(len2=1) then
                             med_seq, [Sequence.length med_seq]
                         else begin
                            let chosen_parent,newdelimiters = 
                             pick_delimiters med1 med2 med_seq in
                            better_capping chosen_parent.seq med_seq newdelimiters 
                         end
                     | [] , [] -> med_seq, []
                     | _ -> 
                             failwith "uncompatible delimiters from two parents' medians"
                 in
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
                      cost1 = total_cost - recost1;
                      cost2 = total_cost - recost2;
                      recost1 = recost1;
                      recost2 = recost2;
                      cost3 = 0; (*for median2, set cost3 to 0. median3 will update this*)
                      sum_cost = total_cost + med1.sum_cost + med2.sum_cost;
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
        cost3 = 0;
        sum_cost = med1.sum_cost + med2.sum_cost;
        delimiter_lst = med1.delimiter_lst
        }]
        end
    else
    let ali_pam = get_breakinv_pam breakinv_pam in          
    match ali_pam.symmetric with 
    | true ->
          let cost12, recost12, med12_ls = find_simple_med2_ls med1 med2
              gen_cost_mat pure_gen_cost_mat alpha ali_pam in
    (*
          let cost21, recost21, med21_ls = find_simple_med2_ls med2 med1
              gen_cost_mat pure_gen_cost_mat alpha ali_pam in
          if cost12 <= cost21 then 
    *)
               cost12, recost12, med12_ls
    (*
           else begin 
              let med12_ls = List.map swap_med med21_ls in 
              cost21, recost21, med12_ls
          end 
    *)
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
    (*
    Printf.printf "single to multi in breakinvAli.ml, check delimiter list:\n [%!";
    List.iter (Printf.printf "%d,") delimiters;
    Printf.printf "]\n%!";
    *)
    let seq_lst = List.map ( fun delimiter ->
        let empty: Sequence.s = Sequence.get_empty_seq () in
        for i = !count to (!count + delimiter)-1 do
            Sequence.prepend empty (Sequence.get seq i)
        done;
        count := !count + delimiter;
        (Sequence.reverse empty)
    ) delimiters in
    List.map (fun s ->
        Sequence.printseqcode s;
        { bkinvt with seq = s; delimiter_lst = [] }
    ) seq_lst

