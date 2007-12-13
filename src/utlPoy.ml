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

let () = SadmanOutput.register "UtlPoy" "$Revision: 2495 $"

let fprintf = Printf.fprintf

let is_existed_code code seq = 
    Sequence.fold (fun existed c -> if c = code then true
                   else existed) false seq

let is_existed_char ch seq = 
    let code = Alphabet.match_base ch Alphabet.nucleotides in 
    is_existed_code code seq

let printDNA seq = 
    Sequence.print stdout seq Alphabet.nucleotides; 
    print_newline ()


let create_gap_seq ?(gap=Alphabet.gap) len = 
    Sequence.init (fun _ -> gap) len



let cmp_num_all_DNA seq = 
    let len = Sequence.length seq in 
    let gap = Alphabet.gap in 
    let num_nu = ref 0 in
    for p = 0 to len - 1 do 
        if (Sequence.get seq p) land gap != gap then num_nu := !num_nu + 1;
    done;
    !num_nu


let cmp_num_not_gap seq = 
    let len = Sequence.length seq in 
    let gap = Alphabet.gap in 
    let num_nu = ref 0 in
    for p = 0 to len - 1 do 
        if (Sequence.get seq p) != gap then num_nu := !num_nu + 1;
    done;
    !num_nu



let cmp_gap_cost indel seq = 
    let num_char = cmp_num_all_DNA seq in 
    match num_char with 
    | 0 -> 0
    | _ ->
          let o, e = indel in 
          o + num_char * e / 100


let cmp_ali_cost (alied_seq1 : Sequence.s) (alied_seq2 : Sequence.s) 
        direction (cost_mat : Cost_matrix.Two_D.m) =

    let opening_cost = 
        match Cost_matrix.Two_D.affine cost_mat with
        | Cost_matrix.Affine o -> o
        | _ -> 0
    in 

    
	let len = Sequence.length alied_seq1 in  
    let gap = Cost_matrix.Two_D.gap cost_mat in 

    let rec sum_up (pre_b1, pre_b2, total_cost) p  = 
        if p >= len then total_cost 
        else begin
		    let b1 = Sequence.get alied_seq1 p in 
		    let b2 =   
                match direction = `Positive with  
                | true -> Sequence.get alied_seq2 p   
                | false -> Sequence.get alied_seq2 (len - 1 - p)  
            in  
		    let cost = Cost_matrix.Two_D.cost b1 b2 cost_mat in   
            let total_cost = 
                if ( (pre_b1 != gap) && (b1 = gap) )  
                    || ( (pre_b2 != gap && b2 = gap)) then total_cost + cost + opening_cost
                else total_cost + cost
            in 
            sum_up (b1, b2, total_cost) (p + 1)

        end
    in 
    sum_up (1, 1, 0) 0



	
let get_empty_seq () = Sequence.create 0
						
let subseq seq start len = 
	match len < 1 with
		| true -> get_empty_seq ()
		| false -> Sequence.sub seq start len




(*==================================*)			
let align2 (seq1 : Sequence.s) (seq2 : Sequence.s) 
        (cost_mat : Cost_matrix.Two_D.m) =

	let len1 = Sequence.length seq1 in
	let len2 = Sequence.length seq2 in
	

	let gap_code = Cost_matrix.Two_D.gap cost_mat in 
	let ext_seq1 = Sequence.init (fun pos -> if pos = 0 then gap_code 
						else Sequence.get seq1 (pos - 1)) (len1 + 1) in 
	let ext_seq2 = Sequence.init (fun pos -> if pos = 0 then gap_code 
						else Sequence.get seq2 (pos - 1)) (len2 + 1) in 

(*	print_endline "Start aligning2";
 	Sequence.print stdout ext_seq1 Alphabet.nucleotides; print_newline ();
	Sequence.print stdout ext_seq2 Alphabet.nucleotides; print_newline ();  *)

	(*Cost_matrix.Two_D.output stdout cost_mat; print_newline (); *)

(*	fprintf stdout "The seq 1 length is %d and the sequence 2 length is %d"
    (Sequence.length ext_seq1) (Sequence.length ext_seq2); print_newline ();  *)


	let ext_alied_seq1, ext_alied_seq2, cost = Sequence.Align.align_2 
        ext_seq1 ext_seq2 cost_mat Matrix.default in 		
(*	Sequence.print stdout ext_alied_seq1 Alphabet.nucleotides; print_newline ();
	Sequence.print stdout ext_alied_seq2 Alphabet.nucleotides; print_newline (); 
	print_endline "End of aligning2";  *)

	let ali_len = Sequence.length ext_alied_seq1 - 1 in 
	let alied_seq1 = subseq ext_alied_seq1 1 ali_len in 
	let alied_seq2 = subseq ext_alied_seq2 1 ali_len in 

	alied_seq1, alied_seq2, cost, ali_len


(*====================================================================*)			
let align3 (seq1 : Sequence.s) (seq2 : Sequence.s) (seq3 : Sequence.s) 
        (cost_cube : Cost_matrix.Three_D.m)= 

	let len1 = Sequence.length seq1 in
	let len2 = Sequence.length seq2 in
	let len3 = Sequence.length seq3 in
(*	fprintf stdout "Aligning sequences: %i %i %i\n" len1 len2 len3; 
    flush stdout; *)
	
	let gap_code = Alphabet.gap in 
	let ext_seq1 = Sequence.init (fun pos -> if pos = 0 then gap_code 
						else Sequence.get seq1 (pos - 1)) (len1 + 1) in 
	let ext_seq2 = Sequence.init (fun pos -> if pos = 0 then gap_code 
						else Sequence.get seq2 (pos - 1)) (len2 + 1) in 
	let ext_seq3 = Sequence.init (fun pos -> if pos = 0 then gap_code 
						else Sequence.get seq3 (pos - 1)) (len3 + 1) in 

    let nuc = Alphabet.nucleotides in 
	Sequence.print stdout ext_seq1 nuc; print_newline ();
	Sequence.print stdout ext_seq2 nuc; print_newline ();
	Sequence.print stdout ext_seq3 nuc; print_newline ();
	
	let ext_alied_seq1, ext_alied_seq2, ext_alied_seq3, cost = 
        Sequence.Align.align_3 ext_seq1 ext_seq2 ext_seq3 cost_cube 
            Matrix.default in 		
	
	Sequence.print stdout ext_alied_seq1 nuc; print_newline ();
	Sequence.print stdout ext_alied_seq2 nuc; print_newline ();
	Sequence.print stdout ext_alied_seq3 nuc; print_newline ();

	print_endline "End of POY align_3";
	let ali_len = Sequence.length ext_alied_seq1 - 1 in 
	let alied_seq1 = subseq ext_alied_seq1 1 ali_len in 
	let alied_seq2 = subseq ext_alied_seq2 1 ali_len in 
	let alied_seq3 = subseq ext_alied_seq3 1 ali_len in 
	alied_seq1, alied_seq2, alied_seq3, cost, ali_len
	
	

let closest_alied_seq alied_parent alied_child c2 = 
    let len = Sequence.length alied_parent in 
    let single_seq = Sequence.init 
        (fun p -> 
             let p_code = Sequence.get alied_parent p in 
             let c_code = Sequence.get alied_child p in 
             Cost_matrix.Two_D.get_closest c2 p_code c_code
        ) len 
    in 
    let cost = Sequence.foldi  
        (fun cost p single_code ->
             let p_code = Sequence.get alied_parent p in 
             cost + (Cost_matrix.Two_D.cost single_code p_code c2)
        ) 0 single_seq 
    in 
    single_seq, cost

	
	
(*==================================*)			
let concat (seq_ls : Sequence.s list) = 
    let total_len = List.fold_left (fun acc_len seq -> 
                        acc_len + Sequence.length seq) 0 seq_ls in

    let concat_seq = Sequence.init (fun index -> -1) total_len in            
	let concat_pos = ref 0 in    
    let copier seq = 
		let len = Sequence.length seq in
       	for pos = 0 to len - 1 do
       		Sequence.set concat_seq !concat_pos (Sequence.get seq pos);
			concat_pos := !concat_pos + 1
        done;
    in
    List.iter copier seq_ls;
    concat_seq
(*==================================*)			



	
let create_subalign2 (seq1 : Sequence.s) (seq2 : Sequence.s) 
        (cost_mat : Cost_matrix.Two_D.m) (start_pos1 : int) (end_pos1 : int) 
        (start_pos2 : int) (end_pos2 : int) = 	

	let len1 = end_pos1 - start_pos1 + 1 in
	let len2 = end_pos2 - start_pos2 + 1 in 
	let subseq1 = Sequence.sub seq1 start_pos1 len1 in 
	let subseq2 = Sequence.sub seq2 start_pos2 len2 in	
	let alied_subseq1, alied_subseq2, cost, ali_len = 
        align2 subseq1 subseq2 cost_mat in 	
	alied_subseq1, alied_subseq2, cost



let dna_gap = Alphabet.get_gap Alphabet.nucleotides 

let get_num_base (seq : Sequence.s) = 
	Sequence.fold (fun num_code code -> 
                       if code != dna_gap then num_code + 1 
                       else num_code) 0 seq
	
	

let delete_gap ?(gap_code = dna_gap) seq = 
	let new_len = Sequence.fold 
        (fun len code -> 
             if code = gap_code then len 
             else (len + 1) ) 0 seq in 

	let new_seq = Sequence.init (fun _ -> -1) new_len in 
	let _ = Sequence.fold (fun new_pos code -> 
                               match code = gap_code with
							   | true -> new_pos
							   | false -> Sequence.set new_seq new_pos code;
								     new_pos + 1) 0 seq in
	new_seq
	


let create_median_gap seq ?(start_pos=(-1)) ?(end_pos=(-1)) cost_mat =
    let start_pos, end_pos = 
        match start_pos with
        | -1 -> 
              let len = Sequence.length seq in
              0, (len - 1)
        | _ -> start_pos, end_pos
    in 
    let gap_code = Cost_matrix.Two_D.gap cost_mat in 

    let med = Sequence.init 
        (fun pos -> 
             let code1 = Sequence.get seq (start_pos + pos) in
             let med_code = Cost_matrix.Two_D.median code1 gap_code cost_mat in 
             med_code) (end_pos - start_pos + 1)
    in
    med


let create_median_seq ?(approx=`BothSeq) alied_seq1 alied_seq2 cost_mat =
    let len = Sequence.length alied_seq1 in 
    let get_median_code pos = 
        let code1 = Sequence.get alied_seq1 pos in 
        let code2 = Sequence.get alied_seq2 pos in          
        match approx with 
        | `First -> code1
        | `Second -> code2
        | `BothSeq ->              
              Cost_matrix.Two_D.median code1 code2 cost_mat
    in

    let median = Sequence.init (fun pos -> get_median_code pos) len in


    let cost = ref 0 in 
    for p = 0 to len - 1 do 
        let code1 = Sequence.get alied_seq1 p in 
        let code2 = Sequence.get alied_seq2 p in         
        cost := !cost + (Cost_matrix.Two_D.cost code1 code2 cost_mat)
    done;
    median, !cost



let create_median_deled_seq seq cost_mat =
    let len = Sequence.length seq in 
    let gap = Alphabet.gap in
    let get_median_code pos = 
        let code = Sequence.get seq pos in 
        if code land gap = 0 then code
        else gap
    in

    let median = Sequence.init (fun pos -> get_median_code pos) len in

    median


let create_median ?(approx=`BothSeq) seq1 seq2 
        ?(s1=(-1)) ?(e1=(-1)) ?(s2=(-1)) ?(e2=(-1)) cost_mat = 

    let s1, e1, s2, e2 =
        match s1 with  
        | -1 ->  
              let len1 = Sequence.length seq1 in 
              let len2 = Sequence.length seq2 in 
              0, (len1 - 1), 0, (len2 -1)        
        | _ -> s1, e1, s2, e2
    in 

    let alied_seq1, alied_seq2, _ = 
        create_subalign2 seq1 seq2 cost_mat s1 e1 s2 e2 
    in
    let alied_med, cost = create_median_seq ~approx:approx alied_seq1 alied_seq2 cost_mat in 
    alied_med, alied_seq1, alied_seq2, cost



let check_repeated_char seq alpha =  
    let len = Sequence.length seq in  
    let rec check_char p1 p2 =  
        if p1 = len then () 
        else if p2 = len then check_char (p1 + 1) (p1 + 2) 
        else if Sequence.get seq p1 = Sequence.get seq p2 then begin   
            let ch = Alphabet.match_code (Sequence.get seq p1) alpha in  
            print_endline ("Character " ^ ch ^ " appears twice in sequence"); 
            Sequence.print stdout seq alpha;  
            print_newline ();  
            failwith "In Breakinv, characters MUST BE NOT duplicated"; 
        end else check_char p1 (p2 + 1)                       
    in  
    check_char 0 1 

    

let create_general_ali code1_arr code2_arr gap_code cost_mat =
(*    print_endline "Create general alignment"; *)

    let len1 = Array.length code1_arr in
    let ext_seq1 = Sequence.init 
        (fun index -> 
             match index with
             | 0 -> gap_code
             | _ -> code1_arr.(index - 1) ) (len1 + 1) 
    in 

    let len2 = Array.length code2_arr in 
    let ext_seq2 = Sequence.init 
        (fun index -> 
             match index with
             | 0 -> gap_code 
             | _ -> code2_arr.(index - 1) ) (len2 + 1) 
    in 
    
	let ext_alied_seq1, ext_alied_seq2, cost = Sequence.Align.align_2 
        ext_seq1 ext_seq2 cost_mat Matrix.default in 		


    let ali_len = Sequence.length ext_alied_seq1 in 
    let alied_seq1 = Sequence.to_array ext_alied_seq1 in 
    let alied_seq1 = Array.sub alied_seq1 1 (ali_len - 1) in 


    let alied_seq2 = Sequence.to_array ext_alied_seq2 in 
    let alied_seq2 = Array.sub alied_seq2 1 (ali_len - 1) in 
    alied_seq1, alied_seq2, cost


(** do not use map of sequences, because it i is running up-down*)
let map f s =
    let len = Sequence.length s in 
    let m_s = Sequence.init (fun _ -> 0) len in
    for p = 0 to len - 1 do
        Sequence.set m_s p (f (Sequence.get s p));
    done; 
    m_s


let of_array code_arr = 
    let len = Array.length code_arr in 
    Sequence.init (fun idx -> code_arr.(idx)) len    

let test_general_ali  () = 
    print_endline "Testing general alignmnet";
    let size1 = 5 in 
    let code1_arr = [|1; 2; 3; 4; 5|] in 

    let size2 = 3 in 

    let code2_arr = [|6; 7; 8|] in 

    let size = size1 + size2 + 1 in 
    let cost_mat = Array.make_matrix size size (999) in 

    let gap_code = 9 in 
    let gap_cost = 100 in
        
    for p1 = 0 to size1 - 1 do
        for p2 = 0 to size2 - 1 do
            cost_mat.( code1_arr.(p1) - 1). ( code2_arr.(p2) - 1) <- 
                code1_arr.(p1) * 10 + code2_arr.(p2)
        done 
    done;
    
    for p1 = 0 to size1 - 1 do
        cost_mat.(code1_arr.(p1) - 1).(gap_code - 1) <- gap_cost;
        cost_mat.(gap_code - 1).(code1_arr.(p1) - 1)<- gap_cost;
    done;

    for p2= 0 to size2 - 1 do
        cost_mat.(code2_arr.(p2) - 1).(gap_code - 1) <- gap_cost;
        cost_mat.(gap_code - 1).(code2_arr.(p2) - 1) <- gap_cost;  
    done;



    let cost_ls = Array.to_list cost_mat in 
    let cost_ls = List.map (fun arr -> Array.to_list arr) cost_ls in       
    let cost_mat = Cost_matrix.Two_D.of_list ~use_comb:false cost_ls in 
    Cost_matrix.Two_D.set_gap cost_mat gap_code;


    let alied_seq1, alied_seq2, cost = 
        create_general_ali code1_arr code2_arr gap_code cost_mat
    in 
    fprintf stdout "General ali cost: %i" cost; print_newline (); 
    Utl.printIntArr alied_seq1;
    Utl.printIntArr alied_seq2



        
let get_single_seq seq c2 = Sequence.select_one seq c2

