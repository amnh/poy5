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
(** This module implements basic functions *)

let () = SadmanOutput.register "Utl" "$Revision: 3160 $"
module IntSet = All_sets.Integers

let large_int = 100000000

let max_seq_len = 50000000

(** Each chromosome is assigned an unique chromosome ID, 
* whenever new chromosome is loaded from input or created during
* the run,gen_chrom_ref_code is increated *)
let gen_chrom_ref_code = ref 0

(** Each sequence or segment is assigned an unique sequence ID, 
* whenever new sequence is loaded from input or created during
* the run,gen_seq_ref_code is increated *)
let gen_seq_ref_code = ref 1000000 

(** Each genome is assigned an unique genome ID, whenever
* new genome is loaded from input or created during
* the run,gen_genome_ref_code is increated *)
let gen_genome_ref_code = ref 0

(** When a new chromosome is loaded from input or created
* during the run, this function returns an integer number 
* as an unique chromosome ID for this chromosome *)
let get_new_chrom_ref_code () = 
    let new_chrom_ref_code = !gen_chrom_ref_code in 
    incr gen_chrom_ref_code;
    new_chrom_ref_code


(** When a new genome is loaded from input or created
* during the run, this function returns an integer number 
* as an unique genome ID for this genome *)
let get_new_genome_ref_code () = 
    let new_genome_ref_code = !gen_genome_ref_code in 
    incr gen_genome_ref_code;
    new_genome_ref_code


(** When a new sequence or segment is loaded from input 
* or created during the run, this function returns 
* an integer number as an unique seq ID for this sequence *)
let get_new_seq_ref_code () = 
    let new_seq_ref_code = !gen_seq_ref_code in 
    incr gen_seq_ref_code;
    new_seq_ref_code

let deref ptr = match ptr with
    | Some content -> content 
    | None -> failwith "It is a null pointer" 
            

(*int to bit array, int = 0~15*)
let break_code in_code =
    if in_code>15 || in_code<0 then
        failwith ("break code only take int from 0 to 15");
    let resarr = Array.make 4 0 in
    let tmp_code = ref in_code in
    for i = 3 downto 0 do
        let base = int_of_float ( 2. ** (float_of_int i) )  in
        if (!tmp_code >= base) then begin
            resarr.(i)<-1;
            tmp_code := !tmp_code - base;
        end;
    done;
    assert( !tmp_code=0);
    resarr

        
(** Given two integer lists [l1] and [l2]
* which are sorted non-decreasing. This function returns 
* true if [l1] is identical [l2], otherwise false *)
let compare_non_dec_list (l1 : int list) (l2 : int list) : bool = 
    if (List.length l1) != (List.length l2) then 
        false
    else begin
        let rec compare rem_l1 rem_l2 = 
            match rem_l1, rem_l2 with
                | [], [] -> true
                | h1::t1, h2::t2 ->
                    if (h1 != h2) then 
                        false
                    else
                        compare t1 t2
                | _ -> false
        in        
        compare l1 l2
    end
        
(** Given an integer array [a], two array indices [pos1] 
* and [pos2], the function computes an integer number 
* as  the sum of array elements from [pos1] to [pos2]. 
* Note that [pos1] is not necessary smaller than [pos2] *)
let get_sum_arr (a : int array) (pos1 : int) (pos2 : int) = 
    let rec add pos en sum = match pos > en with
        | true  -> sum
        | false -> add (pos + 1) en (sum + a.(pos))
    in
    match pos1 < pos2 with
    | true  -> add pos1 pos2 0 
    | false -> add pos2 pos1 0 
    
    

(** Given an array [arr], two indices [start_pos] and [end_pos],
* this function returns an inverted array [arr'] where 
* orders of elements from [start_pos] to [end_pos] are inverted.
* Note that element directions (+/-) are kept unchanged *)
let invert_subarr (arr :  'a array) (start_pos : int) (end_pos : int) = 
    for offset = 0 to (end_pos - start_pos + 1) / 2 - 1 do
        let tmp = arr.(start_pos + offset) in
        arr.(start_pos + offset) <- arr.(end_pos - offset);
        arr.(end_pos - offset) <- tmp;
    done    


(** Given a non-decreasing array [arr] and a [looking_val], this
* function returns an integer as an index of [looking_val]. 
* if [looking_val] is not in the [arr], return (-1). 
* To this end, binary search is employed  *)
let binary_index_search (arr : int array) (looking_val : int) = 
    let rec search l u =
        if l > u then -1
        else begin
            let m = (l + u) / 2 in
            if arr.(m) = looking_val then m
            else
                if arr.(m) < looking_val then search (m + 1) u
                else search l (m - 1)
        end
    in
    search 0 ( (Array.length arr) - 1)


(** Given an array [arr], a [looking_item] and a compare
* function [cmp_fun], this function return an integer as
* an  index of the [looking_item]. If the [looking_item] is not
* in the array [arr], return (-1) *)
let find_index arr looking_item cmp_fun = 
    let len = Array.length arr in 
    let rec find pos = 
        if pos = len then -1
        else
            if cmp_fun looking_item arr.(pos) = 0 then pos
            else find (pos + 1)
    in
    find 0


(** Given two arrays [arr1] and [arr2] and compare function [equal], 
* this function returns [arr1'] and [arr2'] where
* - [arr1'] includes elements of [arr1] which are also in [arr2],
* - [arr2'] includes elements of [arr2] which are also in [arr1] *)
let get_common arr1 arr2 equal = 
    let ls1 = Array.to_list arr1 in
    let ls2 = Array.to_list arr2 in 
    let ls1 = List.filter  (fun code1 -> (find_index arr2 code1 equal) != -1) ls1 in
    let ls2 = List.filter  (fun code2 -> (find_index arr1 code2 equal) != -1) ls2 in
    Array.of_list ls1, Array.of_list ls2

let get_common3 arr1 arr2 arr3 equal =
    let ls1 =
        List.filter
            (fun code1 ->
                ((find_index arr2 code1 equal) != -1) && ((find_index arr3 code1 equal) != -1))
            (Array.to_list arr1)
    and ls2 =
        List.filter
            (fun code2 ->
                ((find_index arr1 code2 equal) != -1) && ((find_index arr3 code2 equal) != -1))
            (Array.to_list arr2)
    and ls3 =
        List.filter
            (fun code3 ->
                ((find_index arr1 code3 equal) != -1) && ((find_index arr2 code3 equal) != -1))
            (Array.to_list arr3)
    in
    Array.of_list ls1, Array.of_list ls2, Array.of_list ls3
    

(** Given an array [arr], an array index [pos] and a new element 
* [new_item], this function returns [arr'] where [new_item] 
* is inserted into array [arr] at position [pos] *)
let insert arr pos new_item = 
    let len = Array.length arr in
    Array.init (len + 1) (fun index -> 
      if index < pos then arr.(index)
      else if index = pos then new_item 
      else arr.(index - 1)) 
    let swap_item pos1 pos2 arr =     
    let new_arr = Array.copy arr in  
    new_arr.(pos1) <- arr.(pos2);
    new_arr.(pos2) <- arr.(pos1);
    new_arr

let int_to_int32_arr (arr :int array) = 
    Array.map (fun x -> Int32.of_int x) arr

let int_to_int32_mat (mat : int array array) = 
    Array.map (fun arr -> int_to_int32_arr arr) mat

let float_to_int_arr (arr: float array) =
    Array.map (fun x -> int_of_float x) arr

let float_to_int_mat (mat : float array array) =
    Array.map (fun arr-> float_to_int_arr arr) mat

let printIntArr (arr : int array) = 
    Array.iter
        (fun x ->
            if x =large_int then Printf.fprintf stdout "  L"
                            else Printf.fprintf stdout "%d," x)
        arr;
    print_newline ()

let printIntArrWithIdx (arr :int array) =
    Array.iteri (fun idx item -> Printf.printf "[%d]:%d,%!"  idx item) arr;
    print_newline()

let printIntMat (arr : int array array) = 
    Array.iter printIntArr arr;
    print_newline ()

let printIntMatWithIdx (arr : int array array) =
    Array.iter printIntArrWithIdx arr;
    print_newline ()

let printIntList2 inlist = 
    Printf.printf "[%!";
    List.iter (fun x -> Printf.printf "%d,%!" x) inlist;
    Printf.printf "] %!"

let printIntList inlist = 
    Printf.printf "[%!";
    List.iter (fun x -> Printf.printf "%d,%!" x) inlist;
    Printf.printf "]\n%!"

let printIntListToFile oc inlist =
    Printf.fprintf oc "[%!";
    List.iter (fun x -> Printf.fprintf oc "%d,%!" x) inlist;
    Printf.fprintf oc "]\n%!"

let printIntListList inlstlst =
    List.iter (fun lst ->
        printIntList lst;
    ) inlstlst

let printIntListList2 inlstlst =
    List.iter (fun lst ->
        printIntList2 lst;
    ) inlstlst;
    print_newline()

let printIntListListList inlstlstlst = 
    List.iter (fun lstlst ->
       printIntListList2 lstlst
    ) inlstlstlst

let create_ls len value = 
    Array.to_list (Array.init len (fun _ -> value))

let get_neg_rev_intlst in_lst =
    List.rev (List.map (fun x -> -x ) in_lst)

let get_abs_intlst in_lst =
    List.map (fun x -> abs x) in_lst

(*return average of a int list, return float*)
let get_avg_of_intlst in_lst =
    let acc = ref 0 in
    List.iter (fun size ->
        acc := size + !acc
    ) in_lst;
    (float !acc) /.(float (List.length in_lst))

let get_min_of_lst in_lst =
    match in_lst with
    | h::t ->
        List.fold_left (fun accmin x ->
            if x<accmin then x
            else accmin
        ) h t
    | [] -> failwith "ERROR: Utl.get_min_of_intlst, empty list "


let get_avg_of_floatlst in_lst =
    let acc = ref 0. in
    List.iter (fun size ->
        acc := size +. !acc
    ) in_lst;
    (!acc) /.(float (List.length in_lst))

(** [remove_nth list n] returns the [n]th element of [list] and [list] with the
    [n]th element removed. *)
let rec remove_nth ?(acc=[]) list n =
    if n = 0
    then (List.hd list, List.rev_append acc (List.tl list))
    else if n > 0
    then remove_nth ~acc:(List.hd list :: acc) (List.tl list) (pred n)
    else raise (Invalid_argument "remove_nth")

    
(** Given two arrays [src_arr], [des_arr], and an array index [pos],
    this function returns [des_arr'] where [src_arr] is inserted
    into [des_arr] at position [pos] *)
let insert_arr src_arr des_arr pos =
    let src_len = Array.length src_arr in 
    let des_len = Array.length des_arr in 
    let arr = Array.make (src_len + des_len) 0 in 
    let cur_pos = ref 0 in 
    let add act_arr p1 p2 =
        for p = p1 to p2 do
            arr.(!cur_pos) <- act_arr.(p);
            cur_pos := !cur_pos + 1;
        done 
    in
    add des_arr 0 (pos - 1);
    add src_arr 0 (src_len - 1);
    add des_arr pos (des_len - 1);
    arr

(** [pairwisep p l] Apply a predicate [p] to a list with the head of a list and
    the lists tail.  Thus, a list of, [a;b;c;d], would apply (p a [b;c;d]), (p b
    [c;d]) ... until one of these is false, which we then return. *)
let rec pairwisep p list = match list with
    | l :: ls ->
          let bools = List.map (p l) ls in
          if List.fold_left (&&) true bools
          then pairwisep p ls
          else false
    | [] -> true

(** Given a list [elem_ls] and an integer number [k],
    this function return a list of [k] elements chosen
    randomly from the list [elem_ls] *)
let get_k_random_elem elem_ls k = 
    let elem_arr = Array.of_list elem_ls in 
    let len = Array.length elem_arr in 
    if (k = 1) || (len = 1)  then 
        [elem_arr.(0)]
    else if (k = 2) then 
        [elem_arr.(0); elem_arr.(len - 1)]
    else begin
        let chosen_ls = ref [elem_arr.(0); elem_arr.(len - 1)] in 
        let rest = ref (len - 2) in 
        for pos = 1 to len - 2 do
            if !rest > 0 then begin
                let chosen_elem =                 
                    match len - 2 - pos + 1 = !rest with
                    | true -> Some elem_arr.(pos)
                    | false ->
                          let ran_num = Random.int (len - 2 - pos + 1) in 
                          if ran_num < !rest then Some elem_arr.(pos)
                          else None
                in 
                match chosen_elem with
                | Some elem -> 
                      chosen_ls := elem::!chosen_ls;
                      rest := !rest - 1
                | None -> ()
            end 
        done;
        !chosen_ls
    end 


(** Given two arrays [arr1], [arr2] and compare function [cmp_fun],
    this function returns true if [arr1] and [arr2] are identical, otherwise
    false *)
let isEqualArr arr1 arr2 cmp_fun = 
    if Array.length arr1 != Array.length arr2 then false
    else begin
        let len = Array.length arr1 in 
        let rec check pos = 
            if pos = len then true
            else 
                if cmp_fun arr1.(pos) arr2.(pos) != 0 then false
                else check (pos + 1)
        in 
        check 0
    end 
        
(** Given an array [arr] and a list [break_ls] : (int * int) list,
    this function returns a list of segments which are broken arcording to the
    [break_ls] *)
let break_array arr break_ls = 
    let rev_seg_ls = 
        List.fold_left 
            (fun rev_seg_ls (sta, en) -> 
                 let len = en - sta + 1 in
                 let seg = Array.init len (fun idx -> arr.(sta + idx)) in
                 seg::rev_seg_ls
            ) [] break_ls
    in 
    List.rev rev_seg_ls

let bigger_int a b =  if (a>b) then a else b

let max_arr arr = 
    Array.fold_left (fun v max_v -> if v > max_v then v else max_v) arr.(0) arr

let min_arr arr = 
    Array.fold_left (fun v min_v -> if v < min_v then v else min_v) arr.(0) arr
    

let printIntSet s = 
    IntSet.iter (fun v -> Printf.fprintf stdout "%i " v) s;
    print_newline (); flush stdout


let get_dir dir =  match dir with 
    | `Positive -> "+" 
    | `Negative -> "-" 
    | _ -> "" 

let factorial n =
        let rec func t acc = 
           if (t>1) then   func (t-1) t*acc
           else acc
        in
        func n 1

let p_m_n m n =
        let rec func t acc =
            if ( t < n ) then func (t+1) (m-t)*acc
            else acc
        in
        func 0 1 

