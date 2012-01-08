(* POY 5.0 Alpha. A phylogenetic analysis program using Dynamic Homologies.   *)
(* Copyright (C) 2011 Andrés Varón, Lin Hong, Nicholas Lucaroni, Ward Wheeler,*)
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

let () = SadmanOutput.register "Block_mauve" "$Revision: 2209 $"

(* A.D. = Aaron E. Darling*)
(* W = weight, R = ratio *)

open Printf
let printIntArr = Utl.printIntArr
let break_code = Utl.break_code
let error_user_message format = Printf.ksprintf (Status.user_message Status.Error) format
let info_user_message format = Printf.ksprintf (Status.user_message Status.Information) format


(*this is a cost matrix between ATGC and ATGC. default matrix used by Mauve *)
let hodx_matrix = [| 
           (*A*) (*C*) (*G*) (*T*) 
  (*A*)  [|  91; -114;  -31; -123|];
  (*C*)  [|-114;  100; -125;  -31|];
  (*G*)  [| -31; -125;  100; -114|];
  (*T*)  [|-123;  -31; -114;   91|];

|]

(*transform poy cost matrix into hodx_matrix in mauve*)
let fill_in_hmatrix cost_mat = 
    let debug = false in
    let ac = Cost_matrix.Two_D.cost 1 2 cost_mat in
    let ag = Cost_matrix.Two_D.cost 1 4 cost_mat in
    let at = Cost_matrix.Two_D.cost 1 8 cost_mat in
    let cg = Cost_matrix.Two_D.cost 2 4 cost_mat in
    let ct = Cost_matrix.Two_D.cost 2 8 cost_mat in
    let gt = Cost_matrix.Two_D.cost 4 8 cost_mat in
    let aa = hodx_matrix.(0).(0) 
    and cc = hodx_matrix.(1).(1)
    and gg = hodx_matrix.(2).(2)
    and tt = hodx_matrix.(3).(3) in
    let r = (float_of_int ((aa+cc+gg+tt)*6)) /. 
    (float_of_int (4*(ac+ag+at+cg+ct+gt))) in
    let ac = int_of_float ((float_of_int ac) *. r) in
    let ag = int_of_float ((float_of_int ag) *. r) in
    let at = int_of_float ((float_of_int at) *. r) in
    let cg = int_of_float ((float_of_int cg) *. r) in
    let ct = int_of_float ((float_of_int ct) *. r) in
    let gt = int_of_float ((float_of_int gt) *. r) in
    hodx_matrix.(0).(1) <- -ac;
    hodx_matrix.(1).(0) <- -ac;
    hodx_matrix.(0).(2) <- -ag;
    hodx_matrix.(2).(0) <- -ag;
    hodx_matrix.(0).(3) <- -at;
    hodx_matrix.(3).(0) <- -at;
    hodx_matrix.(1).(2) <- -cg;
    hodx_matrix.(2).(1) <- -cg;
    hodx_matrix.(1).(3) <- -ct;
    hodx_matrix.(3).(1) <- -ct;
    hodx_matrix.(3).(2) <- -gt;
    hodx_matrix.(2).(3) <- -gt;
    if debug then
        Utl.printIntMat hodx_matrix

(*this is the seed table we are going to use.*)    
let palindromic_spaced_seed_tbl =
    let res = Hashtbl.create 15 in
    let seed5 = [ [1;2;4;6;7];
        [1;4;5;6;9];
        [1;2;5;8;9]; ] in
    Hashtbl.add res 5 seed5;
    Hashtbl.add res 7
    [
        [1;2;5;7;9;12;13];
        [1;3;4;8;12;13;15];
        [1;2;4;8;12;14;15];
    ];
    Hashtbl.add res 9
    [
        [1;2;3;5;8;11;13;14;15];
        [1;2;3;6;9;12;15;16;17];
        [1;2;3;6;8;10;13;14;15];
    ];
    Hashtbl.add res 11
    [
        [1;2;3;5;7;10;13;15;17;18;19];
        [1;2;3;6;9;11;13;16;19;20;21];
        [1;2;3;4;7;9;11;14;15;16;17];
    ];
    Hashtbl.add res 13
    [
        [1;2;3;5;7;8;11;14;15;17;19;20;21];
        [1;2;3;5;8;9;11;13;14;17;19;20;21];
        (* 1111**1**1*1*1**1**1111 this one sucks *)
    ];
    Hashtbl.add res 15
    [
        [1;2;3;4;6;8;9;12;15;16;18;20;21;22;23]
    ];
    Hashtbl.add res 19
    [
        [1;2;3;4;6;7;9;10;11;13;15;16;17;19;20;22;23;24;25]
    ];
    Hashtbl.add res 21
    [
        [1;2;3;4;5;7;8;9;11;12;14;16;17;19;20;21;23;24;25;26;27]
    ];
    res



(******************** seedNO_available_arr and its functions ********************)
(*we recreate this arr every time with [create_lcb_tbl]*)

(*when all seedNO in seedNO_available_arr is in use, we expand it*)
let expand_arr add_len seedNO_available_arr = 
    let add_arr = Array.make add_len 1 in
    let add_lst = Array.to_list add_arr in
    let ori_lst = Array.to_list !seedNO_available_arr in
    let new_lst = ori_lst@add_lst in
    seedNO_available_arr := Array.of_list new_lst

let return_a_seedNO idx seedNO_available_arr =
    let seedNO_arr = !seedNO_available_arr in
    seedNO_arr.(idx)<-1;
    seedNO_available_arr := seedNO_arr

let get_a_seedNO seedNO_available_arr = 
    let found = ref (-1) in
    let idx = ref 1 in
    let seedNO_arr = !seedNO_available_arr in
    let arrlen = Array.length seedNO_arr in
    while (!found=(-1))&&(!idx<arrlen) do
        if (seedNO_arr.(!idx)=1) then begin 
            found := !idx;
            seedNO_arr.(!idx)<-0;
            seedNO_available_arr := seedNO_arr; 
        end;
        idx := !idx +1;
    done;
    if (!found = (arrlen-1)) then 
expand_arr (Array.length !seedNO_available_arr) seedNO_available_arr;
    !found


let radix_sort inarr =
    let max_size = 15 in
    let debug = false in
    let get_idx x = x-1  in
    let size_digits = 
        let _,_,_,first_seq,_ = inarr.(0) in
        Array.length first_seq in
    if debug then Printf.printf "digits size is %d\n%!" size_digits;
    let sort_by_digit inarr digit =
       (*List.sort (fun (_,_,_,lst1) (_,_,_,lst2) ->
            compare (List.nth lst1 digit) (List.nth lst2 digit)
        )inlst*)
       if debug then Printf.printf "sort by digit %d\n%!" digit;
        let count_arr = Array.make max_size 0 in
        Array.iteri (fun i (seqNO,pos,dir,subseq,_) ->
            let x = subseq.(digit) in
            let idx = get_idx x in
            assert((idx>=0)&&(idx<15)); (*chr not in [1,2,4,8]*)
            count_arr.(idx) <- count_arr.(idx)+1; 
            inarr.(i)<-(seqNO,pos,dir,subseq,idx);
        ) inarr;
        let index = ref 0 in
        let full_count_arr = Array.make max_size (0,0,0) in
        let _ = Array.fold_left ( fun pre_count count ->
            full_count_arr.(!index) <- (count,0,pre_count); 
            index := !index +1;
            count+pre_count
        ) 0 count_arr in
        let sorted_inarr = Array.make (Array.length inarr) (0,0,0,[||],0) in
        Array.iteri (fun i (seqNO,pos,dir,arr,idx) ->
            let total_num,hit_num,bigger_than = full_count_arr.(idx) in
            assert(hit_num<total_num);
            let index2 = bigger_than + hit_num in
            sorted_inarr.(index2)<-(seqNO,pos,dir,arr,idx);
            full_count_arr.(idx) <- (total_num,hit_num+1,bigger_than);
        )inarr;
        sorted_inarr
    in
    let digit_arr = Array.init size_digits (fun x->size_digits-x-1) in
    Array.fold_left (fun current_arrarr digit ->
        sort_by_digit current_arrarr digit
    ) inarr digit_arr 
    (*
    Array.fast_sort (fun x y ->
                let (_,_,_,lstx) = x and (_,_,_,lsty)=y in
                compare lstx lsty ) inarr;
    inarr*)
    
    

let extend_seq ext_dir init_matcharr seedweight inseqarr =
    let debug = false in
    if debug then Printf.printf "extend seq in dir=%d,seedW=%d\n%!" ext_dir
    seedweight;
    let new_seedweight = ref seedweight in
    let posarr = init_matcharr in
    let stillmatch = ref true in
    while !stillmatch do
        (*we modify posarr when we have a possible extension between current
        * match and previous match, but this extension might not work for next
        * match, in that case we need to undo the change, 
        * "posarr_copy" is here for the way back*)
        let posarr_copy = Array.init (Array.length posarr)  (fun i -> posarr.(i) ) in
        let idx = ref (-1) in
        let last_chr,last_dir = 
            Array.fold_left (fun (pre_chr,pre_dir) (seqNO,pos,dir) ->
            idx := !idx +1;
            if debug then Printf.printf "check (%d,%d,%d) with (%d,%d)\n%!" 
            seqNO pos dir pre_chr pre_dir;
            if pre_chr<0 then (pre_chr,pre_dir) (*mismatch already happened*)
            else begin
                let inseq = inseqarr.(seqNO) in
                let inseqlen = Array.length inseq in
                let ext2right = if dir*ext_dir=1 then true else false in
                let nexpos = 
                    if ext2right then !new_seedweight+pos
                    else pos-1
                in
                if debug then Printf.printf "nextpos=%d,%!" nexpos;
                if nexpos>=0&&nexpos<inseqlen then begin
                    let next_chr = 
                        if pre_chr=0||dir=pre_dir then inseq.(nexpos)
                        else Alphabet.complement2 inseq.(nexpos) Alphabet.nucleotides
                    in
                    if debug then Printf.printf "chr=%d/%d\n%!" next_chr pre_chr;
                    if pre_chr=0(*the first item*)||next_chr=pre_chr then 
                        begin
                            let _,oldpos,_ = posarr.(!idx) in
                            let newpos =
                                if ext2right then oldpos
                                else oldpos-1
                            in
                            posarr.(!idx) <- (seqNO,newpos,dir);
                            (next_chr,dir)
                        end
                    else (-1,-1)
                end
                else (-1,-1)
            end
            ) (0,0) posarr in
        if last_chr<0 then begin (*undo extension, get out of while loop*)
            Array.iteri (fun i record_copy ->
                posarr.(i) <- record_copy ) posarr_copy;
            stillmatch := false
        end
        else new_seedweight := !new_seedweight + 1;
    done;
    if debug then
        Printf.printf "end of extend_seq, new seedweight = %d\n%!" !new_seedweight;
    posarr,
    !new_seedweight


let extend_seq_in_both_dir matcharr seedweight inseqarr = 
    let extend_matcharr,new_seedweight = 
        extend_seq 1 matcharr seedweight inseqarr in
    let res_matcharr,res_seedweight = 
        extend_seq (-1) extend_matcharr new_seedweight inseqarr in
    Array.sort (fun (seqNOx,_,_) (seqNOy,_,_) -> compare seqNOx seqNOy) res_matcharr;
    res_matcharr,res_seedweight

(*score bwteen two sequences based on hodx matrix *)
let get_score_from_2seq sequence1 sequence2 ori = 
    let debug = false in
    let seq1,seq2 = 
        if ori=1 then sequence1,sequence2
        else 
            Alphabet.rev_comp_arr sequence1 Alphabet.nucleotides, sequence2
    in      
    assert ( (Array.length seq1)=(Array.length seq2) );
    if debug then begin
    Printf.printf "get score from 2seq\n%!";
    printIntArr seq1; printIntArr seq2;
    end;
    Array_ops.fold_right_2 (fun acc x y -> 
        let codearr1 = break_code x 
        and codearr2 = break_code y in
        if debug then begin 
            Printf.printf "x=%d,y=%d,code1arr1/codearr2=\n%!" x y;
        Utl.printIntArr codearr1; Utl.printIntArr codearr2;
        end;
        let maxv = ref (-125) in
        Array.iteri (fun idx1 code1 ->
            Array.iteri (fun idx2 code2 ->
                if code1>0 && code2>0 then
                    if !maxv< hodx_matrix.(idx1).(idx2) then
                        maxv := hodx_matrix.(idx1).(idx2)
            )codearr2
        )codearr1;
        if debug then Printf.printf "acc(%d)+%d=%d\n%!" acc !maxv (!maxv+acc);
        acc + !maxv (*hodx_matrix.(atgc_2_idx x).(atgc_2_idx y) *)
    ) 0 seq1 seq2

