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

let () = SadmanOutput.register "Block_mauve" "$Revision: 2836 $"

(* A.D. = Aaron E. Darling*)
(* W = weight, R = ratio *)
let debug_main = false
let debug_remove_bad_match =  false 
let debug_remove_bad_match2 = false
let debug_build_lcbs = false
let debug_search_outside = false
let debug_join_3_cell = false 
(*now we have 3 different function for removing bad lcbs: greedy,
* dyn1,dyn2. greedy is from A.D.'s paper, dyn1 and dyn2 are my
* idea. 
* greedy is O(n^3), dyn1 and dyn2 are both O(n^2).
* if faster_remove=false, the greedy one will be used. if faster_remove=true,
* we call the dynamic programming one. there is another dynamic
* one "remove_bad_lcbs_dyn2" ,which is slower than "remove_bad_lcbs_dyn"
* my guess is, it's doing bit arr operations, which should be
* faster in c_side. anyway, we are not using dyn1 now.*)
let faster_remove = true


open Printf
open Block_mauve_mum

let error_user_message format = Printf.ksprintf (Status.user_message Status.Error) format
let info_user_message format = Printf.ksprintf (Status.user_message Status.Information) format

let printIntArr = Utl.printIntArr
let printIntMatWithIdx = Utl.printIntMatWithIdx
let printIntMat = Utl.printIntMatWithIdx
let printIntList = Utl.printIntList
let printIntList2 = Utl.printIntList2
let printIntListList = Utl.printIntListList
let printIntListListList = Utl.printIntListListList
let printIntListToFile = Utl.printIntListToFile
let get_abs_lst in_lst = List.sort compare (Utl.get_abs_intlst in_lst)
let get_neg_rev_intlst = Utl.get_neg_rev_intlst
let bigger_int = Utl.bigger_int
let get_avg_of_floatlst = Utl.get_avg_of_floatlst
let get_avg_of_intlst = Utl.get_avg_of_intlst
let get_min_of_lst = Utl.get_min_of_lst
(** for all lcbs. if they all together don't cover more than this of input sequence, adjust
* parameters like minimum_cover_ratio and minimum_lcb_ratio to find more matches  
*)

let checkpoint_lcb p =
  (*
Gc.compact ();
  *)
  prerr_endline ("block_mauve_lcb.ml checkpoint at position " ^ p);
  let slice_number = Gc.major_slice 0 in
  prerr_endline ("major slice number = "^(string_of_int slice_number))


let minimum_cover_ratio = ref 0.4
(*if after outer loop, a lcb is bigger than this, cut it in "search_inside_each_lcb"*)
let maximum_lcb_len = ref 500
(** lcb shorter than this is a light-weight lcb*)
let minimum_lcb_len = ref 50
(*lcb with ratio smaller than this is considered low ratio lcb.*)
let minimum_lcb_ratio = ref 30.0
(** we lower the minimum_lcb_ratio when initial step could not find enough
* match(minimum_full_cover_ratio). we stop doing that until we have enough match, 
* or minimum_lcb_ratio is lower than this number*)

(*we should stop looking when the improvement is too small*)
let minimum_improvement = ref 0.05
(*we stop looking for more lcbs when cover ratio is bigger than this number*) 
let maximum_cover_ratio = ref 0.9
let stop_lcb_ratio = ref 5.

let init_tb_size = 50

let init_seed_size = 50

(*minimum value of break point penalty*)
let min_break_point_penalty = ref 4000


let sanity_check_is_dna_seq seqarr = 
    let res = ref true in
    Array.iter (fun x ->
        if x>=0 && x<=16 then ()
        else res := false;
    ) seqarr;
    if (!res) then ()
    else begin
        Printf.printf "this is not a dna sequence:\n%!";
        printIntArr seqarr;
    end;
    !res

let tiny_improvement newresult oldresult =
    if (oldresult = 0.) then begin
        Printf.printf "we didn't find any good matches before, any improvment \
        this round is good to us\n%!";
        false
    end
    else begin
        if (abs_float (newresult -. oldresult) ) /. oldresult < (!minimum_improvement) 
        then begin
            if debug_main then Printf.printf " tiny improvement, no need to continue\n%!";
            true
        end
        else false;
    end



let get_cons_id in_seq_size_lst num_of_mums = 
    let acc = ref 0.0 in
    List.iter (fun size ->
        acc := (1.0 /. (float size)) +. !acc 
    ) in_seq_size_lst;
    !acc *. (float num_of_mums) /. (float (List.length in_seq_size_lst))


(*if lcb cover ratio is too small, we may need to align super big chunk of
* sequences after mauve, user might want to adjust pameter like 
* minimum_cover_ratio and minimum_lcb_ratio to find more matches*)
let get_adjusted_parameter full_covR in_seq_size_lst minlcbR minlcbL = 
    let e = 2.718 in 
    let divided = 
        if full_covR > 0. then 
            !minimum_cover_ratio /. full_covR 
        else e *. e (*if we didn't find any good matches before, cut the
        requirement by half*)
    in
    if debug_main then Printf.printf "adjust parameter,divided = %f\n%!" divided;
    if divided <= 1. then 
        minlcbR,minlcbL
    else if divided <= e then 
        minlcbR *. (0.9), int_of_float ((float_of_int minlcbL)*.(0.9)) 
    else begin
        minlcbR /. (log divided), int_of_float ((float_of_int minlcbL)*.0.9)
    end


type lcb = {
    seedNOlst : int list; (*each lcb is a list of mums*)
    range_lst : Block_mauve_mum.m_i list; 
    (* range_lst = [range of this lcb in seq0; range of this lcb in seq1].
    lcb_range we need is exactly the same truct as m_i,
    no need to create a new type here
    we sort it in update_lcbs_range_and_score, in decrease range size order*)
    ratio : float; (* [score/length] of seq in this lcb*)
    ref_code : int; (*just a code for bk/rearr purpose*)
    avg_range_len : int; (*average of range length from range_lst *)
    (*score between subsequence contained by range_lst during lcb building.
    * after we have lcb_tbl, huge lcb blocks are aligned seperately, score in
    * function search_inside_each_lcb is set to alignment cost of seq in this lcb*)
    score : int; 
    alignment : Sequence.s array;   
}

let print_lcblst lcbs = (*lcbs should be int list list list.*)
    List.iter (fun lcblst ->
        List.iter (fun record -> printIntList2 record) lcblst;
        Printf.printf "\n%!";
    ) lcbs

let print_lcbs_range in_lstlst =
    Printf.printf "lcb range lst lst :\n%!";
    List.iter (fun lcb_range_lst ->
        Printf.printf "{ %!";
        List.iter (fun (leftend,rightend,lcbkey,lcb_refcode) ->
            Printf.printf "(%d,%d,ref=%d,%!" leftend rightend lcb_refcode;
            printIntList2 lcbkey; 
            Printf.printf ")%!";
        )lcb_range_lst;
        Printf.printf " } \n%!";
    ) in_lstlst

(*function to print lcb*)
let print_lcb x = 
    Printf.printf "LCB.%!";
    printIntList2 x.seedNOlst; Printf.printf "\n%!";
    print_pos_list x.range_lst;
    Printf.printf "score = %d;%!" x.score;
    Printf.printf "ratio = %f;%!" x.ratio;
    Printf.printf "avg range = %d;%!" x.avg_range_len;
    Printf.printf "ref code = %d\n%!" x.ref_code

let print_lcb_to_file oc x = 
    fprintf oc "LCB.%!";
    printIntListToFile oc x.seedNOlst; 
    print_pos_list_to_file oc x.range_lst;
    fprintf oc "score = %d\n%!" x.score;
    fprintf oc "ratio = %f\n%!" x.ratio 

let update_lcb_ref_code lcb_tbl key refcode = 
    let debug = false in
    let record = try (Hashtbl.find lcb_tbl key)
    with |Not_found -> failwith "not found updateing lcb refcode" in
    if debug then begin
        Printf.printf "update refcode %d to %!" refcode;
        printIntList2 key;
        Printf.printf "\n%!";
    end;
    Hashtbl.remove lcb_tbl key;
    Hashtbl.add lcb_tbl key {record with ref_code = refcode}

(*return true if lcb_record is a light weight lcb.*)    
let is_light_weight_lcb lcb_record =
    if (lcb_record.avg_range_len < !minimum_lcb_len) 
    then true
    else false

let is_low_score_lcb item =
    if item.ratio < !minimum_lcb_ratio then true
    else false

(*0: good ratio long len, 1: good ratio short len, 2: poor ratio long len. 3: poor ratio short len.*)
let lcb_quality_test lcb_record =
    let resl = is_light_weight_lcb lcb_record  
    and resr = is_low_score_lcb lcb_record in
    if resl&&resr then 3
    else if resl&&(resr=false) then 1
    else if (resl=false)&&resr then 2
    else 0

(*return the weight of low weight or low ratio lcbs in lcbarr, also return the number of qualified
* lcb in lcbarr *)
let get_bad_block_weight lcbarr lcb_tbl debug = 
        Array.fold_left (fun (acc_w,acc_count) lcbkey ->
            let lcbrecord = try Hashtbl.find lcb_tbl (get_abs_lst lcbkey)
            with 
            |Not_found -> 
                error_user_message " get_light_block_weight,could not find lcb";
                printIntList2 lcbkey;
                assert(false) 
            in
            if debug then print_lcb lcbrecord;
            if (is_light_weight_lcb lcbrecord)||(is_low_score_lcb lcbrecord) then 
                lcbrecord.avg_range_len+acc_w,acc_count
            else
                acc_w,acc_count+1
        ) (0,0) lcbarr 

let both_are_high_score_lcbs lcbkey1 lcbkey2 lcb_tbl debug =
    let lcbrecord1 = try Hashtbl.find lcb_tbl (get_abs_lst lcbkey1)
    with  |Not_found -> 
                error_user_message "could not find lcb in both_are_high_score_lcbs";
                printIntList2 lcbkey1;
                assert(false) 
    in
    let q1 = lcb_quality_test lcbrecord1 in
    if debug then print_lcb lcbrecord1;
    let lcbrecord2 = try Hashtbl.find lcb_tbl (get_abs_lst lcbkey2)
        with  |Not_found -> 
                error_user_message "could not find lcb in both_are_high_score_lcbs";
                printIntList2 lcbkey2;
                assert(false) 
    in
    let q2 = lcb_quality_test lcbrecord2 in
    if q1=0||q1=1 then begin
               if debug then print_lcb lcbrecord2;
                if debug then Printf.printf "q1=%d,q2=%d\n%!" q1 q2;
        (q2=0 || q2=1),q1,q2
    end
    else
        false,q1,q2


let get_from_code2lcbmap (codelst:int list) m =
    let ori = 
        if (List.hd codelst)<0 then (-1) else 1 in
    let codelst = if ori>0 then codelst
    else get_neg_rev_intlst codelst in
    let abs_lcbkeylst = 
        List.map (fun code ->
        try (All_sets.IntegerMap.find code m) 
        with | Not_found -> 
            Printf.printf "code %d is not in code2lcb map\n%!" code;
            assert(false);
    ) codelst in
    abs_lcbkeylst,ori
    (*if ori>0 then abs_lcbkeylst 
    else (get_neg_rev_intlst abs_lcbkeylst)*) 

let print_codemap codemap =
    All_sets.IntegerListMap.iter (fun lcbkey codelist ->
        Printf.printf "lcb:%!";
        printIntList2 lcbkey;
        Printf.printf " --> %!";
        printIntList2 codelist;
        Printf.printf "\n%!";
    ) codemap

let get_code_from_codemap lcbkey codemap = 
    let ori = ref 1 in
    let code =  
    try All_sets.IntegerListMap.find lcbkey codemap 
     with | Not_found -> 
        ( 
            ori := (-1);
            try All_sets.IntegerListMap.find (get_neg_rev_intlst lcbkey) codemap
            with | Not_found -> 
            Printf.printf "could not find lcb in codemap,%!";
            printIntList2 lcbkey;
            assert(false);
        )
    in
    code,!ori

(*[remove_from_codemap lcbkey codemap lcbarr0 debug]. return the codemap after
* removing lcbkey. we have two lcb arrays. after get breakpoint analy*)
let remove_from_codemap lcbkey codemap lcbarr0 debug = 
    let map_after_remove = 
    All_sets.IntegerListMap.remove lcbkey codemap in
    let map_after_remove = 
    All_sets.IntegerListMap.remove (get_neg_rev_intlst lcbkey)
    map_after_remove in
    if debug then begin
    Printf.printf "codemap after remove:\n%!";
    print_codemap map_after_remove;
    end;
    (*say if we remove 6th,7th and 8th lcb so that 1th and 2th lcb could collapse
    * into each other, then 5th and 9th lcb are now adjancent, we should
    * add code6,7,8 to 5th lcb in codemap, so "is_adjacent" on 5th and
    * 9th lcb return true.
    * NOTE: code in codemap is idx+1, so 1th lcb has code 1, but its idx
    * in lcblst is 0, that the reason of -1-1*)
    let codelst2remove,_ = get_code_from_codemap lcbkey codemap in
    let mincode2remove = List.hd (List.sort compare codelst2remove) in
    let idx_ahead = ref(mincode2remove-1-1) 
    and lcbkey_ahead = ref [] in
    while (!idx_ahead>=0)&&(!lcbkey_ahead <> []) do
        let tmplcb = lcbarr0.(!idx_ahead) in
        if All_sets.IntegerListMap.mem tmplcb map_after_remove then
            lcbkey_ahead := tmplcb
        else
            idx_ahead := !idx_ahead -1;
    done;
    let lcbkey_ahead = !lcbkey_ahead in
    if lcbkey_ahead <> [] then begin
        let codelst_ahead,_ = get_code_from_codemap lcbkey_ahead map_after_remove in
        let newlst_ahead = codelst_ahead @ codelst2remove in
        let codemap_ahead = All_sets.IntegerListMap.add lcbkey_ahead
        newlst_ahead map_after_remove in
        if debug then begin
            Printf.printf "mincode2remove=%d,idx_ahead=%d,lcb_ahead:%!" 
            mincode2remove !idx_ahead;
            printIntList2 lcbkey_ahead;
            Printf.printf "codemap ahead:\n%!";
            print_codemap codemap_ahead;
        end;
        codemap_ahead
    end
    (*if there is no available lcb ahead, no need to worry*)
    else map_after_remove

let get_lcblst_from_codemap codemap =
    let debug = false in
    let lst = All_sets.IntegerListMap.fold (fun lcbkey codelst acc ->
        (codelst,lcbkey)::acc;
    ) codemap [] in
    if debug then Printf.printf "get lcblst from codemap, size of map=%d\n%!"
    (List.length lst);
    if debug then
        List.iter (fun (codelst,lcbkey) ->
        printIntList2 codelst;
        Printf.printf " -- > lcblst=%!";
        printIntList2 lcbkey;
        print_newline();
    ) lst;
    let sorted_lst = List.sort(fun (codelst1,lcbkey1) (codelst2,lcbkey2) ->
        let mincode1 = List.hd codelst1 
        and mincode2 = List.hd codelst2 in
        compare mincode1 mincode2
    ) lst in
    List.fold_left (fun acc (codelst,lcbkey) -> 
        if acc=[[]] then [lcbkey] else acc @ [lcbkey]
    ) [[]] sorted_lst 


let get_bkmatrix_codemap inarrarr get_neg_item =
    let debug = false and debug2 = false in
    let sizex = Array.length inarrarr 
    and sizey = Array.length (inarrarr.(0)) in
    if debug then 
        Printf.printf "get_bkmatrix_codemap, sizex=%d,sizey=%d\n%!" sizex sizey;
    assert(sizex=2);  assert(sizey>0);
    
    let bkmatrix = Array.make_matrix sizex sizey (-1) in
    (*codemap : lcb key list -> code list. here code list has only one int, but
    * later in removing bad lcbs, we need to merge two lcblist into one, we
    * merge their code into one code list as well.*)
    let emptycodemap = All_sets.IntegerListMap.empty in 
    let emptycodemap2 = All_sets.IntegerMap.empty in 
    let i = ref 0 in
    let codemap,codemap2 = Array.fold_left (fun (accmap,accmap2) item ->
        let newmap = All_sets.IntegerListMap.add item [(!i+1)] accmap in
        let newmap2 = All_sets.IntegerMap.add (!i+1) item accmap2 in
        bkmatrix.(0).(!i) <- (!i+1);
        i := !i + 1;
        newmap,newmap2 ) (emptycodemap,emptycodemap2) inarrarr.(0) 
    in
    let i = 1 in
    for j=0 to (sizey-1) do
        let item = inarrarr.(1).(j) in
        let ori = ref 1 in 
        let singlecodelst = try All_sets.IntegerListMap.find item codemap 
        with | Not_found ->
            (
            ori := (-1);
            try All_sets.IntegerListMap.find (get_neg_item item) codemap
            with | Not_found -> 
                Printf.printf "i=%d,j=%d,could not find item\n%!" i j;
                printIntList2 item;
                assert(false);
            )
        in
        assert((List.length singlecodelst)=1);
        let code = List.hd singlecodelst in
        let ori = !ori in
        bkmatrix.(i).(j)<-(ori*code);
    done;
    if debug then printIntMatWithIdx bkmatrix;
    if debug2 then print_codemap codemap;
    bkmatrix,codemap,codemap2


let get_break_point_matrix_faster inarrarr get_neg_item lcb_param (*optional*) debug =
    let debug = false in
    let sizex = Array.length inarrarr 
    and sizey = Array.length (inarrarr.(0)) in
    assert(sizex=2);  assert(sizey>0);
    if debug then 
        Printf.printf "get_break_point_matrix_faster, sizex=%d,sizey=%d\n%!" sizex sizey;
    let bkmatrix = Array.make_matrix sizex sizey (-1) in
    let emptycodemap = All_sets.IntegerListMap.empty in 
    let i = ref 0 in
    let codemap = Array.fold_left (fun accmap item ->
        let newmap = All_sets.IntegerListMap.add item (!i+1) accmap in
        bkmatrix.(0).(!i) <- (!i+1);
        i := !i + 1;
        newmap ) emptycodemap inarrarr.(0) 
    in
    let light_block_lst = ref [] (*record bad lcb blocks between good ratio pairs*)
    and worst_weight = ref (!minimum_lcb_len) (*we don't need this now*)
    and pre_pos = ref (-1) 
    and q3_lst = ref [] (*record quality=3 lcbs,we remove this kind of match first*)
    in
    let i = 1 in
    for j=0 to (sizey-1) do
        let item = inarrarr.(1).(j) in
        let ori = ref 1 in 
        let code = try All_sets.IntegerListMap.find item codemap 
        with | Not_found ->
            (
            ori := (-1);
            try All_sets.IntegerListMap.find (get_neg_item item) codemap
            with | Not_found -> 
                Printf.printf "i=%d,j=%d,could not find item\n%!" i j;
                printIntList2 item;
                assert(false);
            )
        in
        let ori = !ori in
        bkmatrix.(i).(j)<-(ori*code);
        let idx = code - 1 in
        (match lcb_param with
            | Some lcb_tbl ->
                if (!pre_pos>0) then begin
                    let a,b = 
                        if !pre_pos > idx then idx,!pre_pos else !pre_pos,idx 
                    in
                    let debug = false in
                    let both_are_good,q_a,q_b = both_are_high_score_lcbs inarrarr.(i).(a) inarrarr.(i).(b)
                    lcb_tbl debug in
                    if q_a=3||q_a=2 then q3_lst := inarrarr.(i).(a) :: !q3_lst ;
                    if q_b=3||q_b=2 then q3_lst := inarrarr.(i).(b) :: !q3_lst ;
                     if (b-a>1)&&both_are_good then 
                        let middle_block_full = 
                            Array.sub inarrarr.(i) (a) (b-a+1) in
                        let middle_block = Array.sub middle_block_full 1
                        ((Array.length middle_block_full)-2) in
                        let w,c = get_bad_block_weight middle_block lcb_tbl
                        debug in
                        if debug then Printf.printf "c=%d,w=%d\n%!" c w;
                        if c=0 && w>0 && !worst_weight<= w then begin
                        light_block_lst := middle_block_full :: !light_block_lst;
                        end
                        else if c=0 && w>0 && !worst_weight>w  then begin 
                        light_block_lst := !light_block_lst @ [middle_block_full];
                        worst_weight := w;
                        end
                end;
                pre_pos := idx;
            | None -> () )
    done;
    if debug then  printIntMatWithIdx bkmatrix;
    let bk_number = ref 0 in
    let i = 1 in
       for j=(sizey-1) downto 1 do
            let b = bkmatrix.(i).(j) and a = bkmatrix.(i).(j-1) in
            if ((b*a>0)) && ( (b-a)=1 ) then begin
                bkmatrix.(i).(j) <- 0;
                let idx = bigger_int (abs a) (abs b) in
                bkmatrix.(0).(idx-1) <- 0;
            end
            else begin
                bk_number := !bk_number + 1;
                bkmatrix.(i).(j) <- 1;
                
            end
        done;
        bkmatrix.(i).(0) <- 1;
        bk_number := !bk_number + 1;
    for j=0 to (sizey-1) do
        if bkmatrix.(0).(j)<>0 then bkmatrix.(0).(j)<-1
        else ();
    done;
    bkmatrix.(0).(0) <- 1;
    bkmatrix, !bk_number, (!light_block_lst,!worst_weight,!q3_lst)
    
(*get range of current lcb , seedNOlst is the list of seedNO of
* that lcb, this returns the leftmost position and right most position of that
* lcb in each seq*)
let get_range_of_a_lcb seedNOlst seqNO mum_tbl seed2pos_tbl= 
    let debug = false in
    if debug then Printf.printf "get range of a lcb :\n%!";
    let min_abs_code = List.hd (List.sort (fun x y -> compare (abs x) (abs y) )
    seedNOlst) in
    let ori =  if (min_abs_code>0) then 1 else (-1) in 
    let leftend,rightend,orirangelst =  
        List.fold_left (fun (le,re,orl) seedNO -> 
            if debug then Printf.printf "add mum#%d,%!" seedNO;
        let thismum = get_mum_from_mumtbl (abs seedNO) mum_tbl seed2pos_tbl in
        let thispos = get_position_by_seqNO thismum.positions seqNO in
        let newle = 
            if thispos.left_end<le then thispos.left_end 
            else le
        and newre =
            if thispos.right_end>re then thispos.right_end
            else re
        in
        newle,newre,
        thismum.positions @ orl
    ) (max_int/2,- max_int/2,[]) seedNOlst in
    if debug then print_pos_list orirangelst;
    (leftend,rightend),ori,orirangelst


let get_sorted_lcb_by_ratio lcb_tbl =
    let lst = ref [] in
    Hashtbl.iter (fun key record ->
        lst := record :: (!lst);
        lst := List.sort (fun x y -> compare x.ratio y.ratio) (!lst)
    ) lcb_tbl;
    !lst


(*[get_num_of_low_ratio_lcb] returns number of low lcb_ratio. is has_high_covR
* is set to true, only cares about lcbs with high covR and low lcb_ratio*)
let get_num_of_low_ratio_lcb lcb_tbl min_ratio in_seq_size_lst has_high_covR =
    let res = ref 0 in
    Hashtbl.iter (fun key record ->
        let check_c = 
            is_light_weight_lcb record in
        if (check_c=false)&&(record.ratio < min_ratio) 
        then begin
            res:= !res +1;
        end;
    ) lcb_tbl;
    !res

(*this only works for two sequences ,for now*)
(*  List.nth lcbs n is the lcb list, each lcb is a list of mums. 
* so lcb list looks like this; [mum a1,mum a2,...],[mum b1,mum b2,...],,..
* Instead of getting rid of all low-weight lcbs, we only get rid one of them.
* "all of the LCBs get scored, then the sum-of-pairs
breakpoint score gets calculated.  Then each one of the LCBs is
evaluated for removal by calculating how much larger the sum-of-pairs
LCB score would be if that block were removed.  The LCB that would
create the largest increase in the sum-of-pairs breakpoint score gets
removed."  
" The total score would be calculated as the sum of the individual MUM scores, minus
(number of breakpoints)*(breakpoint penalty)". *)
(*NOTE: lcb_tbl is EMPTY before this function*)
let update_lcbs_range_and_score lcbs mum_tbl seed2pos_tbl lcb_tbl  = 
    let debug = false in
   (*fill in range of each lcb with mum_tbl*)
    let idx = ref (-1) in
    List.iter (fun lcblist ->
        idx := !idx + 1;
         List.iter (fun record ->
            let (lend,rend),ori,ori_rlst = 
                get_range_of_a_lcb record !idx mum_tbl seed2pos_tbl in
            let abs_record = get_abs_lst record in
            if (Hashtbl.mem lcb_tbl abs_record) then begin
                (*we add range for one seq already, add info for another seq*)
                let old_lcb = Hashtbl.find lcb_tbl abs_record in
                let newlcbrange = 
                { sequence_NO = (!idx); left_end=lend; right_end=rend; orientation=ori;} 
                in
                let old_range_lst = old_lcb.range_lst in
                Hashtbl.remove lcb_tbl abs_record;
                Hashtbl.add lcb_tbl abs_record 
                {old_lcb with 
                range_lst = newlcbrange::old_range_lst;
                };
            end
            else begin (*this is the first range info*)
                assert( !idx = 0 );
                let firstitem = 
                    [{sequence_NO= (!idx);left_end=lend;right_end=rend;
                    orientation=ori}]
                in
                Hashtbl.add lcb_tbl abs_record 
                {seedNOlst = abs_record; 
                range_lst = firstitem; 
                score=0; ratio = 0.0; ref_code = (-1); avg_range_len = 0;
                alignment = [||]}
            end;
        ) lcblist
    )lcbs;
    (* Printf.printf "check lcb range lstlst:\n%!";  print_lcbs_range lcb_range_lst_lst; *)
   (*fill in score for each lcb*)
    Hashtbl.iter (fun seedNOlst oldlcb ->
        let total_score = 
            List.fold_left (fun acc seedNO ->
                 let mumi = get_mum_from_mumtbl seedNO mum_tbl seed2pos_tbl in
                 acc + mumi.mumscore
            ) 0 seedNOlst
        in
        let sorted_range_lst = List.rev ( List.sort (fun x y ->
            compare (x.right_end-x.left_end) (y.right_end-y.left_end)
        ) oldlcb.range_lst ) in
        let avg_range = 
            get_avg_of_intlst 
            (List.map (fun x -> x.right_end-x.left_end+1) sorted_range_lst) in
        let ratio = (float total_score)/.(avg_range) in
        if debug then begin
            Printf.printf "update lcb_tbl on %!" ; print_lcb oldlcb;
            Printf.printf "score <-- %d,ratio <-- %f,avg_range<-%f\n%!" total_score ratio avg_range;
        end; 
        Hashtbl.remove lcb_tbl seedNOlst;
        Hashtbl.add lcb_tbl seedNOlst 
        {oldlcb with 
                     range_lst = sorted_range_lst;
                     score = total_score; 
                     ratio = ratio;
                     avg_range_len = int_of_float avg_range};

    ) lcb_tbl


(*[get_lcbs_range] returns the range lst of current lcb, the list is sorted by leftend. *)
let get_lcbs_range lcb_tbl in_seq_size_lst = 
    let lcb_range_lst_arr = Array.make 2 [] in
    Hashtbl.iter (fun seedNOlst record ->
        if (is_light_weight_lcb record)=false then
            List.iter (fun mi ->
                let oldlst = lcb_range_lst_arr.(mi.sequence_NO) in
                let newlst =
                    (mi.left_end,mi.right_end,seedNOlst,
                    record.ref_code * mi.orientation)
                    ::oldlst in
                lcb_range_lst_arr.(mi.sequence_NO) <- newlst ;
            ) record.range_lst;
    ) lcb_tbl;
    let lcb_range_lst_arr = (*sort the range by leftend*)
        Array.map (fun x -> 
            List.sort (fun (aleft,_,_,_) (bleft,_,_,_) -> compare aleft bleft) x
        ) lcb_range_lst_arr
    in
    Array.to_list lcb_range_lst_arr

(*[get_lcblst_len] return sum of avg_reange_len and score of all qualified lcb 
* in lcbkeylst*)
let get_lcblst_len_and_score lcbkeylst lcb_tbl =
    List.fold_left (fun (acclen,accscore) lcbkey ->
        let tmp = try (Hashtbl.find lcb_tbl (get_abs_lst lcbkey)) 
        with |Not_found -> failwith "not found get lcblst range" in
        let lcb_q = (lcb_quality_test tmp) in
        if lcb_q=1 || lcb_q=0 then
        tmp.avg_range_len + acclen,tmp.score + accscore
        else acclen,accscore
    ) (0,0) lcbkeylst


(*given a lcbs range list list, get_full_range_lstlst returns a range list list
* of both lcbs and area outside lcbs. *)
let get_full_range_lstlst lcb_range_lstlst in_seqsize_lst =
    let res = List.map2 (fun lcb_range_list total_size ->
        let last_rightend,acc=
            List.fold_left (fun (pre_rightend,acc)
            (leftend,rightend,lcbkey,lcb_refcode) ->
                if pre_rightend>leftend then
                    Printf.printf "%d>%d\n%!" pre_rightend leftend;
                assert(pre_rightend<=leftend);
                let newacc = 
                    if leftend>0 then   
                    acc@[((pre_rightend+1),(leftend-1),[],0)] 
                    else acc
                in
                let newacc = newacc@[(leftend,rightend,lcbkey,lcb_refcode)] in
                rightend,newacc 
                ) (-1,[]) lcb_range_list
        in
        let res = 
            if (last_rightend+1=total_size) then acc
            else acc@[(last_rightend+1,total_size-1,[],0)]
        in
        res
    ) lcb_range_lstlst in_seqsize_lst in
    res

(*[get_seq_ouside_lcbs] get sequence outside of list of lcbs, which range is in
* lcb_range_lst_lst. old_seq_range_lst is the range of total sequence. If we are
* building lcb_tbl, old_seq_range_lst is range of input sequences, but if we
* are working inside a huge lcb, old_seq_range_lst is the range of that huge lcb*)    
let get_seq_outside_lcbs old_seqarr lcb_range_lst_lst old_seq_range_lst =
    let debug = false and debug2 = false and debug3 = false in
    if debug then Printf.printf "get seq outside lcbs\n%!";
    if debug2 then checkpoint_lcb "beginning of get seq outside lcbs";
    (*let lcb_range_lst_lst = get_lcbs_range lcb_tbl old_seq_size_lst in*)
    if debug then print_lcbs_range lcb_range_lst_lst;
    let seqNO = ref 0 in
    let seq_outside_lcb_arr = Array_ops.map_2 (fun seq lcb_range ->
        (*let seqlen = Array.length seq in*)
        let leftmost,rightmost = List.nth old_seq_range_lst (!seqNO) in
        if debug then Printf.printf "leftmost=%d,rightmost=%d\n%!" leftmost rightmost;
        seqNO := !seqNO + 1;
        let last_rightend,tmpseq = 
            List.fold_left (
                fun (pre_rightend,accseq) (leftend,rightend,_,_) ->
                assert(pre_rightend<=leftend);
                let size = leftend-pre_rightend-1 in
                let newsubseq =
                    if size>0 then
                    Array.sub seq (pre_rightend+1) (leftend-pre_rightend-1) 
                    else [||]
                in
                if debug3 then assert(sanity_check_is_dna_seq newsubseq);
                if debug2 then checkpoint_lcb ("before append ["^(string_of_int
                leftend)^","^(string_of_int rightend)^"]");
                rightend,Array_ops.array_append accseq newsubseq
            ) (leftmost-1,[||]) lcb_range 
        in
        let last_size = rightmost-last_rightend in
        let subend = 
            if last_size>0 then
            Array.sub seq (last_rightend+1) last_size 
            else [||]
        in
        if debug2 then begin
                Printf.printf "tmpsq size = %d, last_size =%d-%d=%d,[%!" 
		(Array.length tmpseq) rightmost last_rightend last_size;
		Array.iter (fun x -> Printf.printf "%d,%!" x) subend;
		Printf.printf "]\n%!";
		checkpoint_lcb "before append possible last piece";
	end;
        Array_ops.array_append tmpseq subend
    ) old_seqarr (Array.of_list lcb_range_lst_lst) in
    let size_lst = Array.fold_right (fun seq acc ->
        (Array.length seq)::acc ) seq_outside_lcb_arr [] 
    in
    if debug then begin
        Printf.printf "return seq arr with size = %!";
        printIntList2 size_lst;
        if debug3 then
            Array.iter (fun s -> 
                assert(sanity_check_is_dna_seq s);
            ) seq_outside_lcb_arr;
        Printf.printf "return seq outside lcbs\n%!";
    end;
    seq_outside_lcb_arr,size_lst
    

(*[get_lcb_covR_num_badlcb] return the cover rate of qualified lcb, 
* number of un-qualified lcb.  
* max_q : if return of lcb_quality_test is bigger than max_q, it's not a
* qualified lcb -- just put max_q=0 if you need highR&highL lcb; max_q=1 if you
* can accept lowL lcb with highR, etc 
* but if we only have one lcb in lcb_tbl, qualified or not, 
* we get the cover rate of that lcb, num_badlcb doen't matter in this case*)    
let get_lcb_covR_num_badlcb lcb_tbl in_seq_size_lst max_q =
    let q_cov_len = ref 0 in
    let num_badlcb = ref 0 in
    let lcb_tbl_size = Hashtbl.length lcb_tbl in
    Hashtbl.iter (fun key lcb_record ->
        if (lcb_tbl_size=1) (*we keep the last bad lcb*) ||
        (lcb_quality_test lcb_record)<=max_q 
            then q_cov_len := !q_cov_len + lcb_record.avg_range_len
        else 
            num_badlcb := !num_badlcb + 1
    ) lcb_tbl;
    let avg_seq_len = get_avg_of_intlst in_seq_size_lst in
    (float !q_cov_len)/.avg_seq_len , !num_badlcb


(* build_LCBs build lcb table from mum_tb, seed2pos_tbl. 
*  previous_result is the breakpoint number from previous lcbs. it's optional.
*  this is useful when we remove only one lcb from lcbs. we only create a new 
*  lcb table when breakpoint number decrease by at least 2. *)
let build_LCBs seedNOlstlst mum_tbl seed2pos_tbl in_seq_size_lst max_q
(previous_result:int option) =
    let debug = false in
    if debug then Printf.printf "build LCBS \n %!";
    let resmatrix = Array.of_list (List.map (fun lst -> Array.of_list lst)
    seedNOlstlst) in 
    let bkmatrix,_,_ = get_break_point_matrix_faster resmatrix get_neg_rev_intlst None
    false in
    let lcbs = Array.to_list (Array.mapi (fun i bkarr ->
        let res = ref [[]] in
        let j = ref (-1) in
        let lastacc = List.fold_left(fun acc record ->
            j := !j+1;
            if (record=0) then 
                acc @ resmatrix.(i).(!j)
            else if (record=1) then begin
                res := (!res)@[acc];
                resmatrix.(i).(!j)
            end
            else 
                failwith "we only record 1 or 0 in bk_matrix";
        ) [] (Array.to_list bkarr) in
        List.filter (fun item -> (List.length item)>0) (!res@[lastacc]);
    ) bkmatrix )
    in
    if debug then begin
        Printf.printf "build_LCBs,lcblstlst :\n%!";
        assert( ( List.length (List.hd lcbs) ) = (List.length (List.nth lcbs 1))
        );
        printIntListListList lcbs;
    end;
    let new_bk_num = List.length (List.hd lcbs) in
    let old_bk_num = 
        match previous_result with
        | Some x -> x  (*we could use the old lcb_tbl*)
        | None -> new_bk_num+2 (*we need to rebuild lcb_tbl*)
    in
    let lcb_tbl = Hashtbl.create init_tb_size in
    (*only update lcb range and score when two lcb collapse into each other*)
    if (old_bk_num-new_bk_num)>=2 then begin
        update_lcbs_range_and_score lcbs mum_tbl seed2pos_tbl lcb_tbl ;
        let q_cov_rate,num_badlcb = 
            get_lcb_covR_num_badlcb lcb_tbl in_seq_size_lst max_q in
        if debug then 
            Printf.printf "q_cov_rate = %f,end of building LCBs\n%!" 
            q_cov_rate;
        lcbs,q_cov_rate,num_badlcb,lcb_tbl
    end
    else (*no need to build lcb_tbl*)
        lcbs,0.0,-1,lcb_tbl (*this lcb_tbl is empty, we don't need it*)

(* remove one lcb at the time, the rest of lcbs might collapse into each other
* return the score, cover_rate, lcbs and lcb_tbl*)
let filter_out_one_lcb seedNOlst removed_lcb_score lcbs total_mum_score bk_penalty
in_seq_size_lst seed2pos_tbl mum_tbl old_lcb_tbl old_lcb_cov_rate old_badlcb_num =
    assert (seedNOlst<>[]);
    let old_bk_number = List.length (List.hd lcbs) in
    let new_lcbs = List.map (fun lcblst ->
        List.filter (fun record ->( (get_abs_lst record) <> seedNOlst )
         ) lcblst
    ) lcbs in
    let resmatrix = 
        Array.of_list (List.map (fun lst -> Array.of_list lst) new_lcbs) in 
    let debug_bk_matrix = false in
    let _,bk_number,_ = get_break_point_matrix_faster resmatrix get_neg_rev_intlst None debug_bk_matrix in
    let new_score = 
        total_mum_score - removed_lcb_score - bk_number*bk_penalty in
    let new_lcbs, new_lcb_cov_rate,num_badlcb,new_lcb_tbl = 
         build_LCBs new_lcbs mum_tbl seed2pos_tbl in_seq_size_lst 0 (Some old_bk_number ) in
    if new_lcb_cov_rate > 0. then 
        (*after remove the bad lcb, two or more lcbs clapse into each other*)
        new_score,new_lcb_cov_rate,num_badlcb,new_lcbs,new_lcb_tbl
    else begin
        (*nothing new happened, just remove the lcb, and get out*)
        Hashtbl.remove old_lcb_tbl seedNOlst;
        new_score,old_lcb_cov_rate,old_badlcb_num-1,new_lcbs,old_lcb_tbl
    end
    
let get_worst_lcb_from_light_lcbtbl param_shortcut light_lcb_tbl old_lcb_tbl lcbs mum_tbl 
seed2pos_tbl in_seq_size_lst param_init_value =
    let shortcut,lcb_w_worstR,worst_score = param_shortcut in
    let bk_penalty,total_mum_score,init_higher_score,init_better_covR,
    init_num_badlcb,init_low_ratio_lcb_num = param_init_value in
    (*we gonna store best result in these, and return them*)
    let lcb_to_remove = ref [] in
    let lcbs_after_remove = ref [[[]]] in
    let higher_score = ref init_higher_score in
    let better_covR = ref init_better_covR in
    let num_badlcb = ref init_num_badlcb in
    let lcb_tbl_after_remove = ref (Hashtbl.create init_tb_size) in
    if shortcut then begin (*Case 1: shortcut, just remove this one*)
        if debug_remove_bad_match2 then begin
            info_user_message "take a shortcut, get rid of lcb:%!";
            printIntList2 lcb_w_worstR; 
        end;
        let new_score,new_lcb_cov_rate,new_num_badlcb,new_lcbs,new_lcb_tbl =
            filter_out_one_lcb lcb_w_worstR worst_score lcbs total_mum_score
            bk_penalty in_seq_size_lst seed2pos_tbl mum_tbl old_lcb_tbl
            !better_covR !num_badlcb
        in
        higher_score := new_score;
        lcb_to_remove := lcb_w_worstR;
        lcbs_after_remove := new_lcbs;
        better_covR := new_lcb_cov_rate;
        num_badlcb := new_num_badlcb;
        lcb_tbl_after_remove := new_lcb_tbl; (*copy?*)
    end
    else  (*Case 2: iter all low W lcb to find the worst one.*)
        Hashtbl.iter (fun seedNOlst lcb ->
        let old_lcb_tbl_copy = Hashtbl.copy old_lcb_tbl in
        if debug_remove_bad_match2 then begin
            info_user_message "if we remove lcb:"; print_lcb lcb;
        end;
        let new_score,new_lcb_cov_rate,new_num_badlcb,new_lcbs,new_lcb_tbl =
            filter_out_one_lcb lcb.seedNOlst lcb.score lcbs total_mum_score bk_penalty
            in_seq_size_lst seed2pos_tbl mum_tbl old_lcb_tbl_copy !better_covR !num_badlcb
        in
        let num_low_ratio_lcb = get_num_of_low_ratio_lcb new_lcb_tbl
        !minimum_lcb_ratio in_seq_size_lst true in
        if debug_remove_bad_match2 then
            info_user_message "new score=%d,new covR=%f,num_lowR_lcb(%d) should not increase" 
            new_score new_lcb_cov_rate num_low_ratio_lcb;
        (* rules of iteration:
        *  1. don't create more lowR lcbs 
        *  2. if cover rate increase, a decrease in score is ok
        *  3. if cover rate does not increase, we take the biggest increase in score
        *  4. don't let cover rate drop*)
        if (init_low_ratio_lcb_num >= num_low_ratio_lcb) &&
        ( (new_lcb_cov_rate > !better_covR)
        ||(new_lcb_cov_rate = !better_covR)&&(new_score > (!higher_score)))
        then begin
            higher_score := new_score;
            better_covR := new_lcb_cov_rate;
            lcb_to_remove := seedNOlst;
            lcbs_after_remove := new_lcbs;
            num_badlcb := new_num_badlcb;
            lcb_tbl_after_remove := new_lcb_tbl;(*copy?*)
        end;
        ) light_lcb_tbl
    ;
    if debug_remove_bad_match2 then 
        info_user_message "end of get_worst_lcb_from_light_lcbtbl";
    !higher_score,!better_covR,!num_badlcb,!lcb_to_remove,
    !lcbs_after_remove,!lcb_tbl_after_remove
    
 
(*get_worst_lcb returns the worst lcb so we can remove it to get a better score.
* also returns the score, number of low W lcb left if we
* remove the lcb*)
let get_worst_lcb lcbs lcb_tbl bk_penalty in_seq_size_lst mum_tbl seed2pos_tbl
shortcut init_cov_rate init_num_badlcb =
    let total_num = Hashtbl.length lcb_tbl in
    (*num_high_W_low_R : make sure we won't create more low ratio lcb.*)
    let num_high_W_low_R = ref 0 and total_mum_score = ref 0 in
    Hashtbl.iter (fun seedNOlst record ->
            total_mum_score := !total_mum_score + record.score;
            if (is_light_weight_lcb record)=false &&
               (is_low_score_lcb record) then
                   num_high_W_low_R := !num_high_W_low_R + 1;
    ) lcb_tbl;
    let total_mum_score = !total_mum_score in
    let init_bk_num = List.length (List.hd lcbs) in
    let init_score =  total_mum_score - init_bk_num*bk_penalty in
    if debug_remove_bad_match then begin
        info_user_message "get worst lcb, bk_num=%d,init score=%d,init covR =%f;" 
        init_bk_num init_score init_cov_rate;
        info_user_message "total_mum_score = %d,lowR&highW lcb num = %d;total num = %d" 
        total_mum_score !num_high_W_low_R total_num;
    end;
    let lcb_to_remove = ref [] and tmp_score = ref 0 in
    Hashtbl.iter (fun key record -> (*remove highW lcb with lowR first*)
            if (!lcb_to_remove=[])&&(is_low_score_lcb record)&&
            ((is_light_weight_lcb record)=false)
            then begin
                lcb_to_remove := key; tmp_score := record.score;
            end
    ) lcb_tbl;
    (*Case 1: remove high w but low R lcb first*)
    if (!lcb_to_remove<>[]) then begin 
        if debug_remove_bad_match then 
            info_user_message "remove a lowR&highW lcb first";
        let new_score,new_cov_rate,num_badlcb,new_lcbs,lcb_tbl_after_remove =
            filter_out_one_lcb !lcb_to_remove !tmp_score lcbs total_mum_score 
            bk_penalty in_seq_size_lst seed2pos_tbl mum_tbl lcb_tbl 
            init_cov_rate init_num_badlcb
        in
        !lcb_to_remove,new_score, new_cov_rate, 
        num_badlcb,
        new_lcbs,lcb_tbl_after_remove
    end
    (*Case 2: remove the worst lcb from the rest low w lcbs *)
    else begin
        let worst_ratio = ref !minimum_lcb_ratio 
        and lcb_w_worstR = ref [] 
        and worst_score = ref 0 in 
        (*we don't really need a tbl, just keep a list of lowW lcb key and score
        * should be enough*)
        let light_lcb_tbl = Hashtbl.create init_tb_size in
        Hashtbl.iter (fun key record ->
            if shortcut && (record.ratio < !worst_ratio) then begin
                    lcb_w_worstR := record.seedNOlst;
                    worst_ratio := record.ratio;
                    worst_score := record.score;
            end;
            if (is_light_weight_lcb record) then 
                Hashtbl.add light_lcb_tbl key record
        ) lcb_tbl;
        if shortcut&&debug_remove_bad_match then
            info_user_message "worst_ratio=%f,worst_score=%d" !worst_ratio !worst_score;
        if (Hashtbl.length light_lcb_tbl)>1 then begin
            if debug_remove_bad_match then 
                info_user_message "choose one from light_lcb_tbl to remove%!";
            let shortcut =
                if shortcut&&(!lcb_w_worstR=[]) then begin
                    if debug_remove_bad_match then 
                        info_user_message "no shortcut to take, try all 1 by 1";
                    false
                end
                else shortcut
            in
            let param_shorcut = (shortcut,!lcb_w_worstR,!worst_score) in
            let param_init_value = (bk_penalty,total_mum_score,init_score,
            init_cov_rate,init_num_badlcb,!num_high_W_low_R) in
            let higher_score,higher_cov_rate,num_badlcb,lcb_to_remove,lcbs_after_remove,
            lcb_tbl_after_remove = 
                get_worst_lcb_from_light_lcbtbl param_shorcut 
                light_lcb_tbl lcb_tbl lcbs mum_tbl seed2pos_tbl in_seq_size_lst 
                param_init_value
            in
            if debug_remove_bad_match then begin
                info_user_message "remove lcb"; 
                printIntList2 lcb_to_remove;
                info_user_message "score after remove lcb is %d, cov rate=%f,num_badlcb=%d/%d" 
                higher_score higher_cov_rate num_badlcb (Hashtbl.length lcb_tbl_after_remove);
            end;
            (*let light_lcb_num = Hashtbl.length light_lcb_tbl in*)
            lcb_to_remove,higher_score,higher_cov_rate,num_badlcb,(*light_lcb_num,*)
            lcbs_after_remove,lcb_tbl_after_remove
        end
        else begin
            if debug_remove_bad_match then 
                info_user_message "no light weight lcb left, nothing to remove";
            (*in this case cov_rate does not matter, light_lcb_left does not matter*)
            [],init_score, 0., 
            1 (*so we update lcb_tbl in remove_light_weight_lcbs*), 
            lcbs, lcb_tbl
        end
        ; (*end of if (Hashtbl.length light_lcb_tbl)>1 *)
    end




(*break point penaly formula is from Mauve.*)
let get_break_point_penalty input_seq_size_lst num_of_mums =
    let debug = false in
    if debug then 
        Printf.printf "get break point penalty, num_of_mums = %d, %!" num_of_mums;
    let avg_seq_len = get_avg_of_intlst input_seq_size_lst in
    let	avg_seq_len = log( avg_seq_len ) /. log( 2.0 ) in
    let avg_seq_len = 7000. *. avg_seq_len in
    let cons_id = 1. -. (get_cons_id input_seq_size_lst num_of_mums) in
    let bp_estimate = int_of_float (3. *.avg_seq_len*.cons_id*.cons_id*.cons_id*.cons_id) in
    if debug then Printf.printf "bp_estimate = %d \n%!" bp_estimate;
    let res = 
        if (bp_estimate> !min_break_point_penalty) then bp_estimate
        else !min_break_point_penalty
    in
    res

let remove_lcb_from_lcblstlst key_to_remove lcbs =
    List.map (fun lcblst ->
        List.filter (fun record ->
           ( (get_abs_lst record) <> (key_to_remove) )
        ) lcblst
    ) lcbs

(*remove bad lcbs (W is weight, R is ratio)  
* -- high W low R lcbs
* -- low W lcbs
* remove highW lowR lcb first 
* the result lcb_tbl could still have low W lcbs, if remove them will create
* more high W low R lcb, we don't do it -- this is different from Mauve *)
let remove_bad_lcbs lcbs lcb_tbl mum_tbl seed2pos_tbl in_seq_size_lst
num_of_mums =
    let bk_penalty = get_break_point_penalty in_seq_size_lst num_of_mums in 
    if debug_remove_bad_match then begin
        info_user_message "start to remove lowW or lowR lcbs (min_L=%d,min_R=%f)" 
        !minimum_lcb_len !minimum_lcb_ratio;
        if debug_remove_bad_match2 then print_lcblst lcbs;
        if debug_remove_bad_match2 then info_user_message "qualified lcb are:";
    end;
    let init_num_badlcb = ref 0 in
    let total_mum_score = ref 0 in 
    Hashtbl.iter (fun key record ->
        total_mum_score := !total_mum_score + record.score;
        if (is_light_weight_lcb record )||
        (is_low_score_lcb record)
        then  
            init_num_badlcb := !init_num_badlcb + 1
        else  
            if debug_remove_bad_match2 then print_lcb record
    ) lcb_tbl;
    (*pass the init score/bk_num/etc to while loop*)
    let total_mum_score = !total_mum_score in
    let init_bk_num = List.length (List.hd lcbs) in
    let init_score = total_mum_score - init_bk_num*bk_penalty in
    let res_score = ref init_score in
    let res_cov_rate,res_num_badlcb = 
        get_lcb_covR_num_badlcb lcb_tbl in_seq_size_lst 0 in
    let res_cov_rate = ref res_cov_rate in
    let res_lcb_tbl = ref lcb_tbl in
    let res_lcbs = ref lcbs in
    let res_num_badlcb = ref res_num_badlcb in
    let sign = ref ( !res_num_badlcb > 1 ) in
    while (!sign) do
        if debug_remove_bad_match then 
        info_user_message "current score=%d,covR=%f,num_badlcb=%d" 
        !res_score !res_cov_rate !res_num_badlcb;
    (*if we have too many bad matches, we mneed a shortcut to remove them quickly*)
        let num_alllcb = (Hashtbl.length !res_lcb_tbl) in
        assert (num_alllcb>= !res_num_badlcb);
        let shortcut = 
                if !res_num_badlcb=0 then false
                else if (num_alllcb / !res_num_badlcb)<2 then true
                else false
        in
        if debug_remove_bad_match then begin
            info_user_message "we have %d lightW or lowR lcb (among %d),shortcut<-%b "
            !res_num_badlcb num_alllcb shortcut; 
        end;
        let seedNOlst_to_remove,score,cov_rate,num_badlcb,lcbs_after_remove,
        lcb_tbl_after_remove = 
            get_worst_lcb !res_lcbs !res_lcb_tbl bk_penalty in_seq_size_lst
            mum_tbl seed2pos_tbl shortcut !res_cov_rate !res_num_badlcb in
        if (*debug_remove_bad_match&&*)(cov_rate <> !res_cov_rate) then begin
            Printf.printf "remove worst lcb: [%d,..]" (List.hd seedNOlst_to_remove);
            Printf.printf "new score <- %d,cov_rate=%f(>%f), bad lcb num=%d\n" 
            score cov_rate !res_cov_rate num_badlcb;
            print_lcblst lcbs_after_remove;
            (*Printf.printf "old lcbs is:%!";
            print_lcblst !res_lcbs;
            Printf.printf "new lcbs is:%!";
            print_lcblst lcbs_after_remove;*)
        end;
        (*mauve stops when there is only one light weight lcb left*)
        if (seedNOlst_to_remove<>[])&&(num_badlcb>1)&&(!res_score <= score) 
        then begin
            sign := true; 
            res_lcb_tbl := Hashtbl.copy lcb_tbl_after_remove; (*copy?*)
            res_lcbs := lcbs_after_remove;
            if debug_remove_bad_match2 then begin
                info_user_message "current lcbs is:";
                print_lcblst !res_lcbs;
            end; 
            res_score := score;
            res_cov_rate := cov_rate;
            res_num_badlcb := num_badlcb;
            if debug_remove_bad_match then info_user_message "continue with while loop";
        end
        else begin
            sign := false;
            if (seedNOlst_to_remove<>[])&&(num_badlcb=1)&&(!res_score <= score)
            then begin
                if debug_remove_bad_match then 
                    info_user_message "only one lowW lcb left";
                res_lcb_tbl := Hashtbl.copy lcb_tbl_after_remove;
                res_lcbs := lcbs_after_remove;
                if debug_remove_bad_match2 then begin
                    Printf.printf "current lcbs is :%!";
                    print_lcblst !res_lcbs;
                end; 
                res_score := score;
                res_cov_rate := cov_rate;
                res_num_badlcb := num_badlcb;
            end;
            if debug_remove_bad_match then 
                Printf.printf "end of while loop";
        end;
    done; (* while(!sign) *)
     info_user_message "end of remove_bad_lcbs,covR=%f, (%d,%d,%d)"
        !res_cov_rate !res_num_badlcb (List.length (List.hd !res_lcbs))
        (Hashtbl.length !res_lcb_tbl);
    if debug_remove_bad_match then 
        info_user_message "end of remove_bad_lcbs,covR=%f, (%d,%d,%d)"
        !res_cov_rate !res_num_badlcb (List.length (List.hd !res_lcbs))
        (Hashtbl.length !res_lcb_tbl);
    !res_lcbs,!res_cov_rate, !res_lcb_tbl


 

let remove_bad_lcbs_faster lcbs lcb_tbl mum_tbl seed2pos_tbl in_seq_size_lst
num_of_mums old_cov_rate = 
    if debug_remove_bad_match then
    info_user_message "remove bad lcbs faster,lcb_tbl size=%d" (Hashtbl.length lcb_tbl);
    let improving = ref true in
    let res_lcbs = ref lcbs 
    and res_cov_rate = ref old_cov_rate 
    and res_lcb_tbl = ref lcb_tbl in
    while (!improving) do
    if debug_remove_bad_match then
    info_user_message "improving: previous cov rate = %f" !res_cov_rate;
    let resmatrix = 
        Array.of_list (List.map (fun lst -> Array.of_list lst) !res_lcbs) in
    let _,_,(light_block_lst,worst_w,q3_lcb_lst) = 
        get_break_point_matrix_faster resmatrix get_neg_rev_intlst (Some !res_lcb_tbl) false in
    if debug_remove_bad_match then info_user_message 
    "candidate list size:q3=%d,light=%d,worst weight=%d,go through them 1 by 1" 
    (List.length q3_lcb_lst) (List.length light_block_lst) worst_w;
    if q3_lcb_lst<>[] then begin
        if debug_remove_bad_match then 
            let _ = Printf.printf "blocks in q3 lst :%!" in
            List.iter (fun x -> printIntList2 x) q3_lcb_lst;
        let new_lcbs = List.map (fun lcblst ->
            List.filter (fun record ->
                (*get rid of lcbrecord shows up in q3_list*)
                List.fold_left (fun sign lcbkey ->
                    if sign && ((get_abs_lst record)<>(get_abs_lst lcbkey)) then
                    true else false  
                    ) true q3_lcb_lst
            ) lcblst 
        ) !res_lcbs in
        res_lcbs := new_lcbs;
    end
    ;
    let block2remove = ref [] in (*this is just for debug output*)
    if (light_block_lst)<>[] then begin
        if debug_remove_bad_match then
            info_user_message "done with q3 lst, try other light blocks 1 by 1";
        let better_result = 
            (*greedy here: remove the first on in light_block_lst that can increase covR.*)
        List.fold_left (fun better_result light_block ->
            let _,better_cov_rate,_ = better_result in
            if better_cov_rate > !res_cov_rate then
                better_result
            else
                let size = Array.length light_block in
                let good_lcb1,good_lcb2,light_block =
                light_block.(0),light_block.(size-1),
                (Array.sub light_block 1 (size-2)) in
                let light_block = Array.to_list light_block in
                block2remove := light_block;
                let new_lcbs = List.map (fun lcblst ->
                List.filter (fun record ->
                    (*see if lcbrecord shows up in light_block*)
                    List.fold_left (fun sign lcbkey ->
                        if sign && ((get_abs_lst record)<>(get_abs_lst lcbkey)) then
                        true else false  
                        ) true light_block
                ) lcblst 
                ) !res_lcbs in
                Printf.printf "\nif we remove light lcb block:%!" ;
            List.iter (fun lcbkey -> printIntList2 lcbkey) light_block;
                let new_lcbs, new_lcb_cov_rate,num_badlcb,new_lcb_tbl = 
                build_LCBs new_lcbs mum_tbl seed2pos_tbl in_seq_size_lst 0 None
                in 
                (new_lcbs,new_lcb_cov_rate,new_lcb_tbl)
            (*
            let _,better_cov_rate,_ = better_result in
            let size = Array.length light_block in
            let good_lcb1,good_lcb2,light_block =
                light_block.(0),light_block.(size-1),
                (Array.sub light_block 1 (size-2)) in
            let light_block = Array.to_list light_block in
            let new_lcbs = List.map (fun lcblst ->
                List.filter (fun record ->
                    (*see if lcbrecord shows up in light_block*)
                    List.fold_left (fun sign lcbkey ->
                        if sign && ((get_abs_lst record)<>(get_abs_lst lcbkey)) then
                        true else false  
                        ) true light_block
                ) lcblst 
            ) !res_lcbs in
            if debug_remove_bad_match then info_user_message "rebuild lcb table again";
            let new_lcbs, new_lcb_cov_rate,num_badlcb,new_lcb_tbl = 
                build_LCBs new_lcbs mum_tbl seed2pos_tbl in_seq_size_lst None
            in 
            if new_lcb_cov_rate > better_cov_rate then begin
                if debug_remove_bad_match2 then begin
                Printf.printf "remove bad lcb block:%!" ;
                List.iter (fun lcbkey -> printIntList2 lcbkey) light_block;
                Printf.printf "new lcb_tbl,size=%d,covR=%f(pre covR=%f),#badlcb=%d\n%!"
                (Hashtbl.length new_lcb_tbl) new_lcb_cov_rate better_cov_rate num_badlcb;
                end;
                block2remove := light_block;
                (new_lcbs,new_lcb_cov_rate,new_lcb_tbl)
            end
            else better_result
            *)
        ) (!res_lcbs,!res_cov_rate,!res_lcb_tbl) light_block_lst in
        let better_lcbs,better_cov_rate,better_lcb_tbl = better_result in
        if better_cov_rate> !res_cov_rate then begin
            if debug_remove_bad_match then begin
            Printf.printf "remove worst lcb block:%!" ;
            List.iter (fun lcbkey -> printIntList2 lcbkey) !block2remove;
            Printf.printf "new lcb_tbl,size=%d,covR=%f(pre covR=%f),\n%!"
            (Hashtbl.length better_lcb_tbl) better_cov_rate !res_cov_rate;
            Printf.printf "lcb arrarr after removing:%!";
            print_lcblst better_lcbs;
            end;
            res_lcbs := better_lcbs;
            res_cov_rate := better_cov_rate;
            res_lcb_tbl := better_lcb_tbl;
        end
        else begin
            if debug_remove_bad_match then info_user_message 
            "covR drops from %f, improving set to false." !res_cov_rate;
            improving := false;
           
        end; 
    end 
    else begin
        if debug_remove_bad_match then
        info_user_message "could not improve covR by removing lcbs";
        improving := false;
        (*lcbs,old_cov_rate,lcb_tbl*)
    end;
    if !improving=false && q3_lcb_lst<>[] then begin
            let new_lcbs, new_lcb_cov_rate,num_badlcb,new_lcb_tbl = 
            build_LCBs !res_lcbs mum_tbl seed2pos_tbl in_seq_size_lst 0 None
            in
            res_lcbs := new_lcbs;
            res_cov_rate := new_lcb_cov_rate;
            res_lcb_tbl := new_lcb_tbl;
    end;
    done;
    !res_lcbs,!res_cov_rate,!res_lcb_tbl


(*remove bad lcbs with dynamic prog*)
let remove_bad_lcbs_dyn (lcbs:int list list list) lcb_tbl mum_tbl seed2pos_tbl in_seq_size_lst
num_of_mums old_cov_rate = 
let debug = false and debug2 = false in
    if debug then Printf.printf "remove bad lcbs dyn\n%!";
    if debug2 then Hashtbl.iter (fun key record ->
                    print_lcb record 
                ) lcb_tbl;
    let size = List.length (List.hd lcbs) in 
    (*ok, we don't need the whole matrix, just the diagonal cells*)
        
   (* don't alloc the whole matrix 
   * let cov_mat = Array.make_matrix size size
    (([[]],0,0,All_sets.IntegerListMap.empty),(-1)) in *)
    let resmatrix =  Array.of_list 
        (List.map (fun lst -> Array.of_list lst) lcbs) in
    let bkmatrix,codemap,_ = get_bkmatrix_codemap resmatrix get_neg_rev_intlst in
    (*functions start*)
    let is_adjacent codelst1 codelst2 ori1 ori2 =
        if ori1 * ori2 < 0 then false
        else begin
            let code1 = List.hd (List.rev codelst1) 
            and code2 = List.hd codelst2 in
            (code2-code1) = 1
        end
    in
    let print_cell i j ((lcbkeylst,score,len,codemap),dir) printcodemap =
        Printf.printf "cell.%d.%d:%!" i j;
        print_lcblst [lcbkeylst];
        Printf.printf "score = %d,len=%d,dir=%d\n%!" score len dir;
        if printcodemap then
            print_codemap codemap;
    in
    let join_2_cell ((lcbkeylst: int list list),score,len,codemap) (lcbkey:int list) leftjoin =
        let debug_join_2_cell = false in 
        let addlen,addscore = get_lcblst_len_and_score [lcbkey] lcb_tbl in
        if debug_join_2_cell then begin
            Printf.printf "join 2 cell,addscore=%d,addlen=%d;%!" addscore addlen;
            printIntList2 lcbkey;
            print_newline();
        end;
        let reslcblst =
            if addscore > 0 then 
                match leftjoin with 
                | true -> lcbkeylst @ [lcbkey]
                | false -> [lcbkey] @ lcbkeylst 
            else
                lcbkeylst
        in
        let sumscore = addscore + score and sumlen = addlen + len in
        let rescodemap = 
            if addscore > 0 then codemap
            else 
                remove_from_codemap lcbkey codemap resmatrix.(0) false
        in
        reslcblst,sumscore,sumlen,rescodemap
    in
    let merge_two_lcb lcbkey1 lcbkey2 midlcblst code1lst code2lst midcodemap debug_merge_two_lcb =
        let key12 = List.append lcbkey1 lcbkey2 in
        (*get range including midlcblst*)
        let (lend0,rend0),_,_ = get_range_of_a_lcb key12 0 mum_tbl seed2pos_tbl
        and (lend1,rend1),_,_ = get_range_of_a_lcb key12 1 mum_tbl seed2pos_tbl
        in
        let _,score12 = get_lcblst_len_and_score [lcbkey1;lcbkey2] lcb_tbl in
        let range12 = get_avg_of_intlst [(rend0-lend0);(rend1-lend1)]  in
        let ratio = (float score12)/.(range12) in
        if debug_merge_two_lcb then begin
            Printf.printf "merge_two_lcb:%!";
            printIntList2 lcbkey1; printIntList2 lcbkey2;
            Printf.printf "l0,r0=%d,%d;l1,r1=%d,%d;score=%d,range=%f;
            midlcblst to remove:%!" 
            lend0 rend0 lend1 rend1 score12 range12;
            print_lcblst [midlcblst];
        end;
        if ratio >= !minimum_lcb_ratio then begin
            let newcodemap = List.fold_left (fun accmap lcbkey ->
                remove_from_codemap lcbkey accmap resmatrix.(0) false
            ) (midcodemap) midlcblst in
            (*update code12*)
            let code1lst,_ = 
                get_code_from_codemap lcbkey1 newcodemap
            and code2lst,_ = 
                get_code_from_codemap lcbkey2 newcodemap
            in
            let code12lst = get_abs_lst (code1lst @ code2lst) in
            let codemap_remove1 = All_sets.IntegerListMap.remove lcbkey1 newcodemap in
            let codemap_remove2 = All_sets.IntegerListMap.remove
            (get_neg_rev_intlst lcbkey1) codemap_remove1 in
            let codemap_remove3 = All_sets.IntegerListMap.remove lcbkey2
            codemap_remove2 in
            let codemap_remove4 = All_sets.IntegerListMap.remove
            (get_neg_rev_intlst lcbkey2) codemap_remove3 in
            let codemap12 = All_sets.IntegerListMap.add (get_abs_lst key12)
            code12lst codemap_remove4 in
            if debug_merge_two_lcb then begin
                Printf.printf "codemap12 :\n%!";
                print_codemap codemap12;
            end;
            true,(get_abs_lst key12),score12,int_of_float range12,codemap12
        end
        else
            false,[],0,0,All_sets.IntegerListMap.empty
    in
    let join_3_cell midcell lcbkey1 lcbkey2 debug_join_3_cell=
        if debug_join_3_cell then Printf.printf "join 3 cell \n%!";
        let midlcblst,_,_,midcodemap = midcell in
        let code1lst,ori1 = get_code_from_codemap lcbkey1 midcodemap in
        let code2lst,ori2 = get_code_from_codemap lcbkey2 midcodemap in
        if is_adjacent code1lst code2lst ori1 ori2 then begin
            if debug_join_3_cell then Printf.printf "try to merge left and right lcb\n%!";
            match merge_two_lcb lcbkey1 lcbkey2 midlcblst code1lst code2lst
            midcodemap debug_join_3_cell with
            | true, key12,score12,range12,codemap12 ->
                    [key12],score12,range12,codemap12
            | false, _,_,_,_ ->
            if debug_join_3_cell then Printf.printf "did not merge, join2 instead\n%!";
                    let (keym1,scorem1,lenm1,codemapm1) as tmpjoin = 
                        join_2_cell midcell lcbkey1 false in
                    join_2_cell tmpjoin lcbkey2 true 
        end
        else
            let tmpjoin = join_2_cell midcell lcbkey1 false in
            join_2_cell tmpjoin lcbkey2 true
    in
(*build the matrix from bottom to top*)
    (*init the first two diagonal*)
    let diag0 = Array.make (size)
        (([[]],0,0,All_sets.IntegerListMap.empty),(-1)) in
    let diag1 = Array.make (size-1)
        (([[]],0,0,All_sets.IntegerListMap.empty),(-1)) in
    for i = 0 to size-1 do
        let lcbkey : int list = resmatrix.(1).(i) in 
        let lcblen,lcbscore = get_lcblst_len_and_score [lcbkey] lcb_tbl in
        diag0.(i) <- ([lcbkey],lcbscore,lcblen,codemap),(-1);
        (*cov_mat.(i).(i) <- ([lcbkey],lcbscore,lcblen,codemap),(-1);*)
    done;
    for i = 0 to size-2 do
        let cellii,_ = diag0.(i) in
        let couplejoin = join_2_cell cellii resmatrix.(1).(i+1) true
        in
        (*let couplelcb = [resmatrix.(1).(i);resmatrix.(1).(i+1)] in
        let couplelen,couplescore = get_lcblst_len_and_score couplelcb lcb_tbl
        in*)
        diag1.(i) <- couplejoin,(-1);
        (* cov_mat.(i).(i+1) <-
            couplejoin(*couplelcb,couplescore,couplelen,codemap*),(-1);*)
    done;
    let rec walk_up_matrix k diagD diagLR =
        if debug2 then Printf.printf "walk up matrix, k=%d,sizeD=%d,sizeLR=%d\n%!"
        k (Array.length diagD) (Array.length diagLR);
        assert((Array.length diagLR)>=1);
        assert((Array.length diagD)>=2);
        if (k=size) then  (*we are done with righttop cell of the matrix*)
        begin
            assert((Array.length diagLR)=1);
            diagLR.(0)
        end
        else begin
            let diagC = Array.make (size-k)
            (([[]],0,0,All_sets.IntegerListMap.empty),(-1)) in
            for j = 0 to size-k-1 do
if debug2 then Printf.printf "work on cell.%d.%d(size=%d)\n%!" j (j+k) size;
                (*if we need mat.2.5 : solution for 2th,3th,4th,5th lcb in lcblst of
                * sequence1*)
                let cell_left,_ = diagLR.(j) 
                    (* cov_mat.(j).(j+i-1) join mat.2.4 and 5th lcb*) 
                and cell_right,_ = diagLR.(j+1)
                    (*cov_mat.(j+1).(j+i) join 2th lcb and mat.3.5*)
                and cell_down,_ = diagD.(j+1)
                    (* cov_mat.(j+1).(j+i-1) see if we remove lcb
                3th and 4th, lcb 2th and 5th can collapse into one lcb. if not
                just join 2th lcb, mat.3.4 and 5th lcb*)
                in
                let (resL,scoreL,lenL,codemapL) as from_cell_left = 
                    join_2_cell cell_left resmatrix.(1).(j+k) true
                in
                let (resR,scoreR,lenR,codemapR) as from_cell_right = 
                    join_2_cell cell_right resmatrix.(1).(j) false
                in
                let (resD,scoreD,lenD,newcodemap) as from_cell_down = 
                    join_3_cell cell_down resmatrix.(1).(j) resmatrix.(1).(j+k) debug_join_3_cell in
                let best_of_three = (*we should also compare score here?*) 
                    (*get_best_out_of_three lenL lenR lenD scoreL scoreR scoreD*)
                    if lenL>=lenR && lenL>=lenD then 
                        from_cell_left,0
                    else if lenR>=lenL && lenR>=lenD then 
                        from_cell_right,1
                    else
                        from_cell_down,2
                in
                if debug2 then Printf.printf "lenL=%d,lenR=%d,lenD = %d\n%!" lenL
                lenR lenD;
                diagC.(j) <- best_of_three;
                (*cov_mat.(j).(j+i) <- best_of_three;*)
                if debug2 then print_cell j (j+k) best_of_three false;
            done;
            walk_up_matrix (k+1) diagLR diagC
        end
    in
    let (lcblst1, resscore, reslen, resmap),_ = walk_up_matrix 2 diag0 diag1 in
    let avg_in_seq_size = get_avg_of_intlst in_seq_size_lst in
    if debug then begin
        Printf.printf "end of filling in cov_mat,res score=%d,len=%d,covR=%f \
        (lcblst size=%d)\n %!" 
        resscore reslen ((float reslen)/.avg_in_seq_size) (List.length lcblst1);
        if debug2 then print_codemap resmap;
    end;
    let lcblst0 = get_lcblst_from_codemap resmap in
    (*All_sets.IntegerListMap.iter(fun key item ->
        if (List.mem key lcblst1)||(List.mem (get_neg_rev_intlst key) lcblst1)
        then ()
        else begin
            Printf.printf "this is not in lcblst1:%!";
            printIntList2 key; 
            Printf.printf " --> ";
            printIntList2 item;
        end;
    ) resmap;*)       
    assert((List.length lcblst0)=(List.length lcblst1));
    let res_lcbs = [lcblst0;lcblst1] in
    if debug then print_lcblst res_lcbs;
    let res_lcbtbl = Hashtbl.create init_tb_size in
    update_lcbs_range_and_score res_lcbs mum_tbl seed2pos_tbl res_lcbtbl;
    res_lcbs,(float reslen)/.avg_in_seq_size,res_lcbtbl

(*
let remove_bad_lcbs_dyn2 (lcbs:int list list list) lcb_tbl mum_tbl seed2pos_tbl in_seq_size_lst num_of_mums old_cov_rate = 
let debug = false and debug2 = false in
    if debug then Printf.printf "remove bad lcbs dyn2\n%!";
    if debug then print_lcblst lcbs;
    if debug2 then Hashtbl.iter (fun key record ->
                    print_lcb record 
                ) lcb_tbl;
    let size = List.length (List.hd lcbs) in 
    (*ok, we don't need the whole matrix, just the diagonal cells*)
    let cov_mat = Array.make_matrix size size
    (([[]],0,0,Array.create (size+1) 1),(-1)) in
    let resmatrix =  Array.of_list 
        (List.map (fun lst -> Array.of_list lst) lcbs) in
    let bkmatrix,lcb2code,code2lcb = 
        get_bkmatrix_codemap resmatrix get_neg_rev_intlst in
    let print_cell i j ((codelst,score,len,maskarr),dir) printmaskarr =
        Printf.printf "cell.%d.%d:%!" i j;
        print_lcblst [codelst];
        Printf.printf "score = %d,len=%d,dir=%d\n%!" score len dir;
        if printmaskarr then Utl.printIntArr maskarr;
    in 
    let join_2_cell cell1 cell2 trivial leftjoin debug_join_2_cell =
        let codelst1,score1,len1,maskarr1 = cell1 
        and codelst2,score2,len2,maskarr2 = cell2 in
        if debug_join_2_cell then begin
            Printf.printf "join 2 cell ,trivial=%b,leftjoin=%b,cell1/cell2= %!"
            trivial leftjoin;
            print_int_lstlst2 codelst1;
            Printf.printf "score1=%d,len1=%d,maskarr1=\n%!" score1 len1;
            Utl.printIntArr maskarr1;
            print_int_lstlst2 codelst2;
            Printf.printf "score2=%d,len2=%d,maskarr2=\n%!" score2 len2;
            Utl.printIntArr maskarr2;
        end;
        if trivial then begin
            assert((List.length codelst1)=1);
            let singlecode1 = List.hd codelst1 in
            let singlelcb1,ori = get_from_code2lcbmap singlecode1 code2lcb in
            assert((List.length singlecode1)=1);
            let code1 = List.hd singlecode1 in
            let addlen,addscore = get_lcblst_len_and_score singlelcb1 lcb_tbl in
            let res = if addlen>0 then
                match leftjoin with
                | true ->
                        let codelst12 = codelst1@codelst2 in
                        codelst12,score2+addscore,len2+addlen,maskarr2
                | false -> 
                        let codelst21 = codelst2@codelst1 in
                        codelst21,score2+addscore,len2+addlen,maskarr2
            else begin
                if debug_join_2_cell then Printf.printf "----------- set maskarr.%d to 0\n%!" code1;
                let resmaskarr = Array.copy maskarr2 in
                resmaskarr.(abs code1)<-0;
                codelst2,score2,len2,resmaskarr
            end
            in
            if debug_join_2_cell then begin
                let (reslst,resscore,reslen,resarr) = res in
                Printf.printf "res of join2,score12=%d,len12=%d,codelst/maskarr12=%!" resscore reslen;
                print_int_lstlst2 reslst;
                Utl.printIntArr resarr;
            end;
            res
        end
        else begin
            let maskarr12 = Array_ops.map_2 (fun a b -> a land b) maskarr1 maskarr2
            in
            let tmplst12 = codelst1 @ codelst2 in
            let _,_,_,codelst12,score12,len12 = 
                List.fold_right (fun this (pre,pre_score,pre_len,acclst,accscore,acclen) ->
                let lcblstt,_ = get_from_code2lcbmap this code2lcb in
                let len_t,score_t = get_lcblst_len_and_score lcblstt lcb_tbl in
                if pre = [] then 
                    this,score_t,len_t,[this],score_t,len_t
                else begin
                    (*join [codeth,codet2,...codett] and [codeph,codep2....codepp]
                     * codetx is list of int code *)
                    let code_ph = List.hd pre 
                    and code_tt = List.hd (List.rev this) in
                    if debug_join_2_cell then begin
                        Printf.printf "(acc=%d,%d;pre=%d,%d;this=%d,%d),pre.head=%d,this.tail=%d -> %!" 
                        accscore acclen pre_score pre_len len_t score_t code_ph code_tt;
                        Printf.printf "this="; printIntList2 this;
                        Printf.printf "pre="; printIntList2 pre;
                        Printf.printf "acc="; print_int_lstlst2 acclst;
                    end;
                    let count = 
                        if code_ph * code_tt < 0 then 3 
                        else begin
                        let absp = abs code_ph and abst = abs code_tt in
                        let startcode,endcode = if abst>absp then 
                            absp,abst else abst,absp in
                        let adjtp = Array.make (abs(absp-abst)+1) 1 in
                        let masktp = Array.sub maskarr12 startcode
                        (endcode-startcode+1) in
                        assert ((Array.length adjtp)=(Array.length masktp));
                        let contp = Array_ops.map_2 (fun a b->a land b) 
                        adjtp masktp in
                        let res = Array.fold_left (fun acc x -> acc+x) 0 contp
                        in
                        if res<2 then begin
                        Printf.printf "startcode=%d,endcode=%d\n%!" startcode
                        endcode;
                        Utl.printIntArr masktp;
                        end;
                        res
                        end;
                    in
                                        assert(count>=2);
                    if count=2 then begin
                        if debug_join_2_cell then
                            Printf.printf " are adjacent to each other,%!";
                        let lcblstp,_ =
                            get_from_code2lcbmap pre code2lcb
                        in
                        let lcbtp_spr = lcblstt @ lcblstp in
                        let lcbtp = List.flatten lcbtp_spr in
                        let (lend0,rend0),_,_ = get_range_of_a_lcb lcbtp 0
                        mum_tbl seed2pos_tbl 
                        and (lend1,rend1),_,_ = get_range_of_a_lcb lcbtp 1
                        mum_tbl seed2pos_tbl in
                        let _,score_tp = get_lcblst_len_and_score lcbtp_spr
                        lcb_tbl in
                        let range_tp =  
                            get_avg_of_intlst [(rend0-lend0)+1;(rend1-lend1)+1] in
                        let range_tp = range_tp in
                        let ratio_tp = (float score_tp)/.(range_tp) in
                        if debug_join_2_cell then 
                            Printf.printf "tp=%f(%f,%d)\n%!" ratio_tp range_tp score_tp;
                        if ratio_tp >= !minimum_lcb_ratio then
                            let rescodelst = (this@pre) :: (List.tl acclst) in
                            this@pre,score_t,len_t,
                            rescodelst,
                            accscore-pre_score+score_tp,
                            acclen-pre_len+(int_of_float range_tp)
                        else
                            this,score_t,len_t,
                            this :: acclst, 
                            accscore+score_t,
                            acclen+len_t
                    end
                    else
                        this,score_t,len_t,
                        this :: acclst,
                        accscore+score_t,
                        acclen+len_t 
                end
            ) tmplst12 ([],0,0,[[]],0,0) in
            if debug_join_2_cell then begin
                Printf.printf "res of join2,score12=%d,len12=%d,codelst/maskarr12=%!" score12 len12;
                print_int_lstlst2 codelst12;
                Utl.printIntArr maskarr12;
            end;
             codelst12,score12,len12,maskarr12 
        end;
    in
    let join_3_cell midcell leftcell rightcell debug_join_3_cell =
        if debug_join_3_cell then Printf.printf "join 3 cell start\n%!";
        let (codelstM,scoreM,lenM,maskarrM) = midcell in
        let (codelstL,scoreL,lenL,maskarrL) = leftcell 
        and (codelstR,scoreR,lenR,maskarrR) = rightcell in
        (*normal join is redundant
        let join_ml = join_2_cell mc lc true false debug_join_3_cell in
        let res_join = join_2_cell join_ml rc false false debug_join_3_cell in*)
        let arrMcopy = Array.copy maskarrM in
        List.iter (fun codeMid ->
            let c = abs codeMid in
            arrMcopy.(c) <- 0; 
            Printf.printf "set code%d to 0;%!" c;
        ) (List.flatten codelstM);
        let arrLcopy = Array_ops.map_2 (fun a b -> a land b) maskarrL arrMcopy in
        let arrRcopy = Array_ops.map_2 (fun a b -> a land b) maskarrR arrMcopy in
        let leftcell_removeMid = codelstL,scoreL,lenL,arrLcopy 
        and rightcell_removeMid = codelstR,scoreR,lenR,arrRcopy in
        let res_removeMid = join_2_cell leftcell_removeMid rightcell_removeMid
        false false true in
        let rescodlst,resscore,reslen,_ = res_removeMid in
        if debug_join_3_cell then Printf.printf "resscore=%d,reslen=%d\n%!" resscore reslen;
        if debug_join_3_cell then print_int_lstlst2 rescodlst;
        res_removeMid
    in
    for i = 0 to size-1 do
        let lcbkey : int list = resmatrix.(1).(i) in 
        let lcblen,lcbscore = get_lcblst_len_and_score [lcbkey] lcb_tbl in
        if debug2 then Printf.printf "fill in cell.%d.%d,score=%d,len=%d\n%!" i i lcbscore lcblen;
        cov_mat.(i).(i) <- ([[bkmatrix.(1).(i)]],lcbscore,lcblen,
        Array.create (size+1) 1),(-1);
    done;
    for i = 0 to size-2 do
        if debug2 then Printf.printf "fill in cell.%d.%d\n%!" i (i+1);
        let cell1,_ = cov_mat.(i).(i) in
        let cell2,_ = cov_mat.(i+1).(i+1) in
        let couplejoin = join_2_cell cell1 cell2 true true false in
        cov_mat.(i).(i+1) <- (couplejoin,(-1));
    done;
    for i = 2 to size-1 do
        for j = 0 to size-1-i do (*join c.j.k-1, c.k.k and c.k+1.j+i-1*)
            if debug then Printf.printf "+++++++++++ work on cell.%d.%d ++++++++++\n%!" j (j+i);
            let bestres = ref ([],0,0,[||]) and bestk1 = ref (-1) and bestk2 = ref (-1) in
            (*for k = j to j+i-1 do
                let debug_join = true in
                let l1,l2 = j,k in
                let r1,r2 = k+1,j+i
                in
                let cellleft,_ = cov_mat.(l1).(l2)
                and cellright,_ = cov_mat.(r1).(r2) in
                if debug then Printf.printf 
                "========== k=%d, join cell.%d.%d and cell.%d.%d %! ======= \n%!" k l1 l2 r1 r2;
                let (codelstk,scorek,lenk,maskarrk) as resk =
                    if k=j then
                    join_2_cell cellleft cellright true true debug_join
                    else if k=j+i-1 then
                    join_2_cell cellright cellleft true false debug_join
                    else
                    join_2_cell cellleft cellright false false debug_join
                in
                *)
            let debug_join = false in
            let lenperbk = ref 0 and basebknum = ref 0 in 
            let cellleft,_ = cov_mat.(j).(j) and cellright,_ = cov_mat.(j+1).(j+i) in
            let (codelstk,scorek,lenk,maskarrk) = 
                join_2_cell cellleft cellright true true debug_join 
            in
            let bknumk = List.length codelstk in
            let lpk = lenk/bknumk in
            let _,_,bestlen,_ = !bestres in
            if bestlen=0 || lenk>bestlen then begin (*bestlen is 0 here*)
                    bestres := (codelstk,scorek,lenk,maskarrk); 
                    bestk1 := j; bestk2 :=j;
                    lenperbk := lpk;
                    basebknum := bknumk;
            end;
            let cellleft,_ = cov_mat.(j).(j+i-1) and cellright,_ = cov_mat.(j+i).(j+i) in
            let (codelstk,scorek,lenk,maskarrk) =
            join_2_cell cellright cellleft true false debug_join in
            let bknumk = List.length codelstk in
            let _,_,bestlen,_ = !bestres in
            if bestlen=0 || lenk>bestlen then begin
                    bestres := (codelstk,scorek,lenk,maskarrk); 
                    bestk1 := j+i; bestk2 :=j+i;
                    lenperbk := lenk/bknumk;
                    basebknum := bknumk;
            end;
            let bkbonus = !lenperbk * 3/2 in
            if debug then Printf.printf "bkbonus update to %d\n%!" bkbonus;
            for k1 = j+1 to (*j+1*)j+i-2 do
                for k2 = k1(*j+i-1*) to j+i-1 do 
                let l1,l2 = 
                    j,k1-1 
                and r1,r2 =
                    k2+1,j+i in
                let m1,m2 = k1,k2 in
                let cellleft,_ = cov_mat.(l1).(l2)
                and cellright,_ = cov_mat.(r1).(r2) in
                if debug then Printf.printf 
                "========== k=%d/%d, join cell.%d.%d and cell.%d.%d %! =======\n"
                k1 k2 l1 l2 r1 r2;
                let (codelstk,scorek,lenk,maskarrk),bk_num_mid =
                    if k1=j then
                    join_2_cell cellleft cellright true true debug_join,0
                    else if k1=j+i then
                    join_2_cell cellright cellleft true false debug_join,0
                    else
                    let cellmid,_ = cov_mat.(m1).(m2) in
                    let codelstmid,_,_,_ = cellmid in
                    (*let _ =
                        print_cell l1 l2 cov_mat.(l1).(l2) true;
                        print_cell m1 m2 cov_mat.(m1).(m2) true;
                        print_cell r1 r2 cov_mat.(r1).(r2) true;
                    in*)
                    join_3_cell cellmid cellleft cellright debug_join,
                    (List.length codelstmid)
                in
                (*******)
                let _,_,bestlen,_ = !bestres in
                let bkreduce = !basebknum - (List.length codelstk) -
                bk_num_mid in
                let lenk = lenk + bkreduce * bkbonus in
                if debug then Printf.printf "basebknum=%d,bkreduce=%d\n%!" !basebknum bkreduce;
                if debug then Printf.printf "lenk = %d (bestlen=%d)\n%!" lenk
                bestlen;
                if bestlen=0 || lenk>bestlen then begin
                    bestres := (codelstk,scorek,lenk,maskarrk); 
                    bestk1 := k1; bestk2 :=k2;
                end;
                done;
            done;
            cov_mat.(j).(j+i) <- (!bestres,!bestk1); 
            if debug then begin
                Printf.printf "++++++++++++++++++++++++++ bestk=%d/%d,fill in cell:%!" !bestk1 !bestk2;
                print_cell j (j+i) cov_mat.(j).(j+i) false;
            end;
        done;
    done;
    (*right top cell of matrix has lcblst1 now, just get lcblst0 from codemap
    * and right top cell, build lcb_tbl, return it*)
    let (rescodelst, resscore, reslen, resmaskarr),_ = cov_mat.(0).(size-1) in
    let avg_in_seq_size = get_avg_of_intlst in_seq_size_lst in
    if debug then begin
        Printf.printf "end of filling in cov_mat,res score=%d,len=%d,fakecovR=%f \
        rescodelst (size=%d)\n %!" 
        resscore reslen ((float reslen)/.avg_in_seq_size) (List.length rescodelst);
        if debug then print_lcblst [rescodelst];
        if debug2 then Utl.printIntArr resmaskarr;
    end;
    let lcblst0 = 
        let idx = ref (size+1) in
        List.fold_right (fun lcbkey acc->
            idx := !idx - 1;
                if resmaskarr.(!idx)>0 then 
                    if acc = [[]] then [lcbkey] else lcbkey :: acc
                else acc
        ) (List.hd lcbs) [[]]
    in
    let lcblst1 = List.fold_right (fun lcbkey acc ->
        let ori = ref 1 in
        let code = try All_sets.IntegerListMap.find lcbkey lcb2code 
        with | Not_found ->
            (
            ori := (-1);
            try All_sets.IntegerListMap.find (get_neg_rev_intlst lcbkey) lcb2code
            with | Not_found ->
            printIntList2 lcbkey;
            failwith "no found lcb2code map" 
            )
        in
            assert((List.length code)=1);
            let code = List.hd code in
            if resmaskarr.(code)=1 then 
                if acc = [[]] then [lcbkey] else lcbkey :: acc
            else acc
    ) (List.nth lcbs 1) [[]]
    in
    assert((List.length lcblst0)=(List.length lcblst1));
    let newlcbs = [lcblst0;lcblst1] in
    let reslcbs,rescovR,_,reslcbtbl = 
    build_LCBs newlcbs mum_tbl seed2pos_tbl in_seq_size_lst 1 None in
    if debug then begin
        Printf.printf "return covR=%f, lcbs:" rescovR;
        print_lcblst reslcbs;
    end;
    reslcbs,rescovR,reslcbtbl
*) 


(* search area outside existing lcb block. *)    
let search_outside_lcbs inner_lcbs lcb_tbl mum_tbl
pos2seed_tbl_left pos2seed_tbl_right seed2pos_tbl
in_seqarr in_seq_size_lst seedNO_available_arr=
    let debug2 = false in
    if debug_search_outside then begin 
        info_user_message "search outside lcbs, lcb_tbl size = %d" 
        (Hashtbl.length lcb_tbl);
        if debug2 then Hashtbl.iter (fun key record -> print_lcb record ) lcb_tbl;
    end;
    (*this works even when the seq outside is shorter than seed's length*)
    let lcb_range_lstlst = get_lcbs_range lcb_tbl in_seq_size_lst in
    let in_seq_range_lst = List.map (fun size -> (0,size-1) ) in_seq_size_lst in
    let seq_outside_lcbs_arr,current_seq_size_lst = 
        get_seq_outside_lcbs in_seqarr lcb_range_lstlst in_seq_range_lst in
    (*let current_seq_size_lst = Array.fold_right (fun seq acc ->
        (Array.length seq)::acc ) seq_outside_lcbs_arr [] 
    in*)
    if debug_search_outside then begin
        Printf.printf "seq outside matches, size :";
        printIntList2 current_seq_size_lst;
    end;
    let get_shorter_len seqarr =
	    let len0 = Array.length seqarr.(0) 
	    and len1 = Array.length seqarr.(1) in
	    if len0<len1 then len0 else len1
    in
    let shorted_seqlen = get_shorter_len seq_outside_lcbs_arr in
    if shorted_seqlen < !minimum_lcb_len then
        inner_lcbs, (- 0.0),lcb_tbl, mum_tbl, pos2seed_tbl_left, pos2seed_tbl_right,
        seed2pos_tbl, (-1)
    else begin
    let new_seedlen = 
        get_proper_seedlen (get_min_of_lst current_seq_size_lst) in
    if debug_search_outside then 
        info_user_message "use seedlen:%d" new_seedlen;
    (*these hashtable are based on the sub-seq list*)
    let new_seedNO2seq_tbl = Hashtbl.create init_tb_size in
    let new_pos2seed_tbl_left = Hashtbl.create init_tb_size in
    let new_pos2seed_tbl_right = Hashtbl.create init_tb_size in
    let new_seed2pos_tbl = Hashtbl.create init_tb_size in
    let new_mum_tbl = Hashtbl.create init_tb_size in
    let new_seedweight = 
        build_seed_and_position_tbl seq_outside_lcbs_arr new_seedlen 
        new_seedNO2seq_tbl new_pos2seed_tbl_left new_pos2seed_tbl_right 
        new_seed2pos_tbl new_mum_tbl seedNO_available_arr false in
    if debug_search_outside then 
        info_user_message "build_seed_and_position_tbl done(newseed weight=%d)" new_seedweight;
    build_local_mums2 new_mum_tbl new_seed2pos_tbl 
    new_pos2seed_tbl_left new_pos2seed_tbl_right false false;
    (*transposed new mums and add them to the old mum_tbl we have*)
    if debug_search_outside then
        info_user_message "transpose new matches, join them with old matches";
    let res_mum_tbl = Hashtbl.copy mum_tbl in
    let res_pos2seed_tbl_left = Hashtbl.copy pos2seed_tbl_left in
    let res_pos2seed_tbl_right = Hashtbl.copy pos2seed_tbl_right in
    let res_seed2pos_tbl = Hashtbl.copy seed2pos_tbl in
    transpose_mum_back lcb_range_lstlst new_mum_tbl new_seedNO2seq_tbl
    res_mum_tbl res_pos2seed_tbl_left res_pos2seed_tbl_right res_seed2pos_tbl seedNO_available_arr;
    if debug_search_outside then 
        info_user_message "extend seeds in both direction,(total seed number = %d)" (Hashtbl.length res_seed2pos_tbl);
    extend_seeds res_mum_tbl res_seed2pos_tbl res_pos2seed_tbl_left res_pos2seed_tbl_right;
    if debug_search_outside then 
        info_user_message "resolve overlap mum";
    let num_of_mums = resolve_overlap_mum res_mum_tbl  (new_seedweight-1)
    res_pos2seed_tbl_left res_pos2seed_tbl_right res_seed2pos_tbl seedNO_available_arr false false in
    update_score_for_each_mum res_mum_tbl in_seqarr;
    (*get the new weight of lcbs*)
    let new_seedNOlstlst = 
        get_mum_lst_for_each_seq res_mum_tbl res_seed2pos_tbl res_pos2seed_tbl_left
        (Array.length in_seqarr) in_seq_size_lst in
    let maxq = 1 in (*this is inner loop, before any remove_badlcb functions,
    we are ok with lowL&highR lcbs.*)
    let new_lcbs,new_covR,num_badlcb,new_lcb_tbl = 
        build_LCBs new_seedNOlstlst res_mum_tbl res_seed2pos_tbl 
        in_seq_size_lst maxq None in
    if debug_search_outside then 
        info_user_message "new covR raise to %f, end of searching outside lcbs,check tbl \
        size (%d,%d,%d,%d)" new_covR (Hashtbl.length res_mum_tbl)
        (Hashtbl.length res_seed2pos_tbl) (Hashtbl.length res_pos2seed_tbl_left)
        (Hashtbl.length res_pos2seed_tbl_right);
    new_lcbs,new_covR,new_lcb_tbl,res_mum_tbl, 
    res_pos2seed_tbl_left, res_pos2seed_tbl_right, res_seed2pos_tbl, num_of_mums
    end

let rec get_init_lcbs seedNOlstlst 
seed2pos_tbl mum_tbl 
in_seq_size_lst init_num_mums 
min_lcb_ratio min_lcb_len previous_fullcovR=
    let debug = false in
    if debug then Printf.printf "get initial lcbs with min_lcbR=%f,min_lcblen=%d\n%!"
    min_lcb_ratio min_lcb_len;
    assert(min_lcb_ratio > 0. );
    minimum_lcb_len := min_lcb_len ; 
    minimum_lcb_ratio := min_lcb_ratio ;
    let init_lcbs,init_covR,num_badlcb,init_lcb_tbl = 
        build_LCBs seedNOlstlst mum_tbl seed2pos_tbl in_seq_size_lst 1 None in
    let lcbnum = List.length (List.hd init_lcbs) in
    if debug then Printf.printf "before remove bad lcbs, init_covR = %f\n%!" init_covR;
    let better_lcbs, better_covR, better_lcb_tbl = 
    if (lcbnum>2) then begin
        if faster_remove then
            remove_bad_lcbs_dyn init_lcbs init_lcb_tbl mum_tbl 
            seed2pos_tbl in_seq_size_lst init_num_mums init_covR 
        else
            remove_bad_lcbs init_lcbs init_lcb_tbl mum_tbl 
            seed2pos_tbl in_seq_size_lst init_num_mums ;
    end
    else
        init_lcbs,init_covR,init_lcb_tbl
    in
    if debug then begin
        Printf.printf 
        "after remove bad lcbs, we have %!";
        (*Hashtbl.iter (fun key record ->
            if (is_light_weight_lcb record ) then ()
            else
                print_lcb record 
        ) better_lcb_tbl;*)
        Printf.printf "better lcb covR = %f\n%! " better_covR;
    end;
    if better_covR < !minimum_cover_ratio then begin
        let new_lcbR,new_lcbL = get_adjusted_parameter better_covR in_seq_size_lst
        min_lcb_ratio min_lcb_len in
        if (tiny_improvement previous_fullcovR better_covR)||(new_lcbR <
        !stop_lcb_ratio) then begin
            if debug then 
                Printf.printf "Warning: not much improve (%f<=>%f), or lcb ratio\
                drops too low(%f<=>%f), just return initial lcbs.\n%!" 
                previous_fullcovR better_covR new_lcbR !stop_lcb_ratio;
            better_lcbs,init_covR,better_covR,better_lcb_tbl
        end
        else
            get_init_lcbs seedNOlstlst 
            seed2pos_tbl mum_tbl 
            in_seq_size_lst init_num_mums
            new_lcbR new_lcbL better_covR;
    end
    else 
        begin
        let lcb2remove = ref [] in
        Hashtbl.iter (fun key record ->
            if ((Hashtbl.length better_lcb_tbl)>1)&&
            (is_light_weight_lcb record ) 
            then begin
                lcb2remove := key :: !lcb2remove;
                Hashtbl.remove better_lcb_tbl key;
            end
        ) better_lcb_tbl;
        let lcb2remove = !lcb2remove in
        let better_lcbs = List.map (fun lcblst ->
            List.filter (fun record -> 
                (List.mem (get_abs_lst record) lcb2remove)=false
            ) lcblst
        ) better_lcbs in
        better_lcbs,init_covR(*we need covR before remove bad lcbs*),
        better_covR,better_lcb_tbl
        end



(*[merge_ali_seq] is for function [search_inside_a_lcb] working on a big lcb block. 
* after we have aligned seq for smaller lcbs, we need to combine it with seq
* outside those lcbs*)
let merge_ali_seq in_lst0 in_lst1 seq0 seq1 cost_mat use_ukk
total_range_lst  =
    let debug = false and debug2 = false in
    let leftmost0,rightmost0 = List.hd total_range_lst 
    and leftmost1,rightmost1 = List.nth total_range_lst 1 in 
    if debug then
        Printf.printf "merge_ali_seq,[lmost,rmost] = [%d,%d];[%d,%d]\n%!" 
        leftmost0 rightmost0 leftmost1 rightmost1;
    let out_cost,last_rightend0,last_rightend1,resseq0,resseq1 =
    List.fold_left2 (fun (acc_cost,pre_rightend0,pre_rightend1,accseq0,accseq1) 
    ((leftend0,rightend0,_,_),aliseq0) ((leftend1,rightend1,_,_),aliseq1) ->
        let outlen0 = leftend0-pre_rightend0-1
        and outlen1 = leftend1-pre_rightend1-1 in
        if debug then Printf.printf
        "[pre_rightend=%d,le=%d,re=%d],[pre_rightend=%d,le=%d,re=%d]\n%!" 
         pre_rightend0 leftend0 rightend0 pre_rightend1 leftend1 rightend1;
        let out_s0,out_s1,out_c,_ = 
            let os0 = Sequence.subseq seq0 (pre_rightend0+1) outlen0 in
            let os1 =Sequence.subseq seq1 (pre_rightend1+1) outlen1 in
            Sequence.align2 os0 os1 cost_mat use_ukk
        in
        acc_cost+out_c,rightend0,rightend1,
        Sequence.concat [accseq0;out_s0;aliseq0],
        Sequence.concat [accseq1;out_s1;aliseq1]
    ) (0,leftmost0-1,leftmost1-1,Sequence.get_empty_seq (),Sequence.get_empty_seq ()) 
    in_lst0 in_lst1
    in
    if debug then Printf.printf 
    "last_rightend0=%d,last_rightend1=%d\n%!" last_rightend0 last_rightend1 ;
    assert(rightmost0=last_rightend0);
    assert(rightmost1=last_rightend1);
    let resarr = [|resseq0;resseq1|] in
    let len0,len1 = (Sequence.length resarr.(0)),(Sequence.length resarr.(1)) in
    if debug then begin
        Printf.printf "merge_ali_seq ends, check res seq (size=%d,%d):\n%!" len0 len1;
        if debug2 then Sequence.printseqcode resarr.(0);
        if debug2 then Sequence.printseqcode resarr.(1);
    end;
    assert(len0=len1);
    resarr,out_cost


let search_inside_a_lcb lcbrecord seq0 seq1 min_len max_len mum_tbl seed2pos_tbl cost_mat use_ukk =
    let debug = false in
    if debug then 
        Printf.printf "Search inside a lcb ,min and max len = %d/%d\n%!" min_len max_len;
        if debug then Printf.printf "work on lcb:%!";
        if debug then print_lcb lcbrecord;
        let key = lcbrecord.seedNOlst in
        let lcblen = lcbrecord.avg_range_len in        
        let accinlen0 = ref 0 and accinlen1 = ref 0 in
        if (lcblen>min_len)&&(lcblen<max_len) then begin
            let sorted_total_rangelst = List.sort (fun x y -> compare x.sequence_NO
            y.sequence_NO) lcbrecord.range_lst in
            let rlst0,rlst1,in_cost = 
            List.fold_left (fun (rlist0,rlist1,acc_cost) seedNO ->
                let thismum = get_mum_from_mumtbl (abs seedNO) mum_tbl seed2pos_tbl in
                let mi0 = get_position_by_seqNO thismum.positions 0 in
                let mi1 = get_position_by_seqNO thismum.positions 1 in
                let le0,re0 = mi0.left_end,mi0.right_end 
                and le1,re1 = mi1.left_end,mi1.right_end in
                let subseq0,subseq1 = 
                    Sequence.sub seq0 le0 (re0-le0+1),
                    Sequence.sub seq1 le1 (re1-le1+1)
                in
                accinlen1 := !accinlen1 + re1-le1+1;
                accinlen0 := !accinlen0 + re0-le0+1;
                let alied_seq0, alied_seq1, cost, _  =  
                Sequence.align2 subseq0 subseq1 cost_mat use_ukk
                in
                if debug then begin
                let len1,len2 =Sequence.length alied_seq0,(Sequence.length
                alied_seq1) in
                assert(len1=len2);
                end;
                ((le0,re0,[seedNO],mi0.orientation),alied_seq0)::rlist0,
                ((le1,re1,[seedNO],mi1.orientation),alied_seq1)::rlist1,
                acc_cost + cost
            ) ([],[],0) key in
            if debug then Printf.printf "cost from smaller lcbs = %d\n%!" in_cost;
            let sort_by_leftend lst = 
                (*this sorting might use bunch of memory, do we really need to
                * sort? mums in lcbkey are suppose to be sorted by leftend.*)
                List.sort (fun ((aleft,_,_,_),_) ((bleft,_,_,_),_) -> 
                    compare aleft bleft) lst
            in
            let rlst0,rlst1 = sort_by_leftend rlst0,sort_by_leftend rlst1 in
            let range_lst0,_ = List.split rlst0 
            and range_lst1,_ = List.split rlst1 in 
            if debug then print_lcbs_range [range_lst0;range_lst1];
            let total_range_lst =
                List.map (fun mi -> (mi.left_end,mi.right_end) ) sorted_total_rangelst
            in
            let aliseqarr,out_cost = 
                merge_ali_seq rlst0 rlst1 seq0 seq1 cost_mat use_ukk total_range_lst 
            in
            if debug then 
            Printf.printf "cost from seq outside smaller lcbs = %d\n%!" out_cost;
            aliseqarr,
            in_cost+out_cost
        end
        else begin
            if debug then 
                Printf.printf "this lcb is too large(>%d) or too small(<%d) to work with\n%!"
                max_len min_len;
            [||],0
        end



let create_lcb_tbl in_seqarr min_lcb_ratio min_lcb_len min_cover_ratio
max_lcb_len cost_mat use_ukk =
    maximum_lcb_len := max_lcb_len;
    minimum_cover_ratio := min_cover_ratio;
    let debug2 = false in
    (*we keep array of available seed#, mark each not-inused # "1"*)
    let seedNO_available_arr = ref (Array.make init_seed_size 1) in
    (*seedNO_available_arr := Array.make init_seed_size 1;*)
    (*output result to file ...    
    * let outfile = "outfile.txt" in let oc = open_out outfile in*)
    let in_seqarr_size = 2 in
    let in_seq_size_lst = Array.fold_right 
    (fun seq acc -> (Array.length seq)::acc ) in_seqarr [] 
    in
    let min_seqlen = get_min_of_lst in_seq_size_lst in
    let seedlen = get_proper_seedlen min_seqlen in 
    if debug_main then 
        Printf.printf "\n **********Block_mauve.create_lcb_tbl;seqlen=%d,%d;seedlen=%d \
        ****************\n%!" (List.hd in_seq_size_lst) (List.nth
        in_seq_size_lst 1 ) seedlen;
    let init_size = init_tb_size in 
    (*let seq2seedNO_tbl = Hashtbl.create init_size in*)
    let seedNO2seq_tbl = Hashtbl.create init_size in
    (* Hashtbl { pos(seqNO,idx) -> seedNO,seedweight,orientation }, 
    * seqNO is the NO of sequence, idx is the position of that
    * sequence.  pos2seed_tbl_left records the left_end of each seed,
    * pos2seed_tbl_right records the right_end of each seed *)
    let pos2seed_tbl_left = Hashtbl.create init_size in
    let pos2seed_tbl_right = Hashtbl.create init_size in
    let seed2pos_tbl = Hashtbl.create init_size in
    let mum_tbl = Hashtbl.create init_tb_size in
(*i move these two out because we need them after we have lcb tbl, in search
* inside a huge lcb*)
    let outer_mum_tbl = ref (Hashtbl.create init_tb_size ) in 
    let inner_seed2pos_tbl = ref (Hashtbl.create init_tb_size ) in
    (*if seedweight is longer than shortest seq, no need to find seeds*)
    let init_num_mums =
        let seedweight = get_seed_weight seedlen in
        if seedweight>min_seqlen then begin
            if debug_main then 
                Printf.printf "seedweight=%d>min seqlen=%d,get out of mauve" seedweight min_seqlen;
            0
        end
        else begin
            (*find initial mums*)
            let seedweight = build_seed_and_position_tbl in_seqarr seedlen 
            (*seq2seedNO_tbl*) seedNO2seq_tbl pos2seed_tbl_left
            pos2seed_tbl_right seed2pos_tbl mum_tbl seedNO_available_arr false in
            if debug_main then 
                Printf.printf "base seedweight=%d, call build_local_mums2 \n%!" seedweight;
            build_local_mums2 mum_tbl seed2pos_tbl 
            pos2seed_tbl_left pos2seed_tbl_right false false;
            if debug_main then
                Printf.printf "++++++++ init seedweight=%d, end of building mum\
                table\n%!" seedweight;
            let init_num_mums = resolve_overlap_mum mum_tbl (seedweight-1) 
            pos2seed_tbl_left pos2seed_tbl_right seed2pos_tbl seedNO_available_arr false false in
            update_score_for_each_mum mum_tbl in_seqarr;
            init_num_mums
        end
    in
    (*to do : too many if else here, do something*)
    let outer_lcb_tbl = ref (Hashtbl.create 0) in
    let outer_lcbs = ref [[[]]] in
    let outer_old_covR = ref 0. in
    (*trivial case, no mum was found at all*)
    if (init_num_mums=0) then 
        info_user_message "cannot find any mums, algn the sequence at once.\n%!"
    else begin
    if debug_main then
        Printf.printf "++++++++++ resolve_overlap_mum done (tbl size\
        check:%d,%d,%d,%d) ++++ \n%!" (Hashtbl.length mum_tbl) (Hashtbl.length
        seed2pos_tbl) (Hashtbl.length pos2seed_tbl_left) (Hashtbl.length
        pos2seed_tbl_right);
    if debug2 then print_mumtbl mum_tbl true;
    let seedNOlstlst = get_mum_lst_for_each_seq mum_tbl seed2pos_tbl 
    pos2seed_tbl_left in_seqarr_size in_seq_size_lst in 
    if debug2 then begin
        Printf.printf "init seedNOlstlst is :\n%!";
        printIntListListList  seedNOlstlst;
    end;
    (*get initial lcbs, and init cover ratio, before and after remove bad lcbs*)
    let init_lcbs,init_covR_before,init_covR_after,init_lcb_tbl = 
        get_init_lcbs seedNOlstlst seed2pos_tbl mum_tbl in_seq_size_lst
        init_num_mums min_lcb_ratio min_lcb_len 0.0 in
    (* sign of any improvement in outer&inner while loops*)
    let any_improvement_outer = ref false in
    (*init inner tbl*)
    let inner_old_covR = ref init_covR_before in
    let inner_lcbs = ref init_lcbs in
    let inner_lcb_tbl = ref (Hashtbl.copy init_lcb_tbl) (* not necessary *)in
    let inner_mum_tbl = ref mum_tbl in
    let inner_pos2seed_tbl_left = ref pos2seed_tbl_left in
    let inner_pos2seed_tbl_right = ref pos2seed_tbl_right in
    inner_seed2pos_tbl := seed2pos_tbl;
    (*init outer tbl*)
    outer_old_covR := init_covR_after ;
    outer_lcbs := init_lcbs ;
    outer_lcb_tbl := (Hashtbl.copy init_lcb_tbl) (* necessary *) ;
    outer_mum_tbl := mum_tbl ;
    let current_num_of_mums = ref init_num_mums in
    let outer_sign = 
        ref ((init_covR_after< !maximum_cover_ratio)&&(init_covR_after>= !minimum_cover_ratio)) in
    if debug_main && (init_covR_after< !minimum_cover_ratio) then
        Printf.printf "the init round could not cover more than %f of all \
        sequence, skip the part of removing bad lcbs\n%!" !minimum_cover_ratio;
    (*to do : create an inner function and an outer one, move this part out*)
    while (!outer_sign) do
        if debug_main then
        Printf.printf "\n begin of outer while, old_covR = %f\n%!" !outer_old_covR;
    (*****work on regions outside of lcbs******)
    (* we do need to update inner lcbs with the outer one, so we can work on new
    * range of subsequence, but inner_old_weight should remain the same, 
    * since inner_weight is the weight before any low-weight lcb removing *)
        inner_lcbs := !outer_lcbs;
        inner_mum_tbl := !outer_mum_tbl; (*this one is not necessary*)
        inner_lcb_tbl := Hashtbl.copy !outer_lcb_tbl;
        let inner_sign = ref true in
        let any_improvement_inner = ref false in
        while (!inner_sign) do
            if debug_main then
                Printf.printf "\n ----------- Inner loop -------------- \n%!";
            let inner_new_lcbs,inner_new_covR,res_lcb_tbl,res_mum_tbl,
            res_pos2seed_tbl_left, res_pos2seed_tbl_right, res_seed2pos_tbl, 
            num_of_mums = 
            search_outside_lcbs !inner_lcbs !inner_lcb_tbl !inner_mum_tbl 
            !inner_pos2seed_tbl_left !inner_pos2seed_tbl_right 
            !inner_seed2pos_tbl in_seqarr in_seq_size_lst seedNO_available_arr 
            in
            if (num_of_mums>0)&&(!inner_old_covR < inner_new_covR) then begin
                if debug_main then
                    Printf.printf "we found better lcbs (%f>%f), update \
                    num_of_mums from %d to %d,continue with inner loop\n%!" 
                    inner_new_covR !inner_old_covR !current_num_of_mums num_of_mums;
                if (inner_new_covR > !maximum_cover_ratio)||(tiny_improvement
                inner_new_covR !inner_old_covR) 
                    then inner_sign :=false
                else inner_sign := true;
                current_num_of_mums := num_of_mums;
                inner_old_covR := inner_new_covR;
                inner_lcbs := inner_new_lcbs;
                inner_lcb_tbl := Hashtbl.copy res_lcb_tbl;
                inner_mum_tbl := res_mum_tbl;
                inner_pos2seed_tbl_left := res_pos2seed_tbl_left;
                inner_pos2seed_tbl_right := res_pos2seed_tbl_right;
                inner_seed2pos_tbl := res_seed2pos_tbl;
                any_improvement_inner := true;
            end
            else begin
                inner_sign := false;
                if debug_main then Printf.printf "no improve on lcb (%f<=%f) \
                get out of inner loop\n%!" inner_new_covR !inner_old_covR ;
            end;
        done; 
        (*end of inner while loop*)
        if debug_main then Printf.printf "\n --------  Outer loop, remove bad lcbs ----- \n%!";
        if !any_improvement_inner && (List.length (List.hd !inner_lcbs))>2 then begin
            let new_outer_lcbs, new_outer_covR, new_outer_lcb_tbl = 
                if faster_remove then
                    remove_bad_lcbs_dyn !inner_lcbs !inner_lcb_tbl !inner_mum_tbl 
                !inner_seed2pos_tbl in_seq_size_lst !current_num_of_mums
                !inner_old_covR 
                else 
                    remove_bad_lcbs !inner_lcbs !inner_lcb_tbl !inner_mum_tbl 
                !inner_seed2pos_tbl in_seq_size_lst !current_num_of_mums
            in
            (*we could still have low W lcb here -- if remove them will create high W
            * low R lcbs, "remove_light_weight_lcbs" won't do it, get rid of those
            * lcbs here
            * Note. mauve keep the last low W lcb, we do the same, so
            * new_outer_weight will at least be 1*)
            Hashtbl.iter (fun key record ->
                if ((Hashtbl.length new_outer_lcb_tbl)>1)&&
                (is_light_weight_lcb record ) 
                then begin
                    (*if debug_main then begin
                        Printf.printf "remove lowW lcb:%!";
                        printIntList2 key;
                    end;*)
                    Hashtbl.remove new_outer_lcb_tbl key;
                end
            ) new_outer_lcb_tbl;
            let maxq = 0 in (*we only count highR&highL lcbs here*)
            let new_outer_covR2,_ = get_lcb_covR_num_badlcb new_outer_lcb_tbl
            in_seq_size_lst maxq in
            if debug_main then Printf.printf "new_outer_covR=%f,new_outer_covR2=%f\n%!"
            new_outer_covR new_outer_covR2;
            (* remove corresponding item in new_outer_lcbs*)
            let new_outer_lcbs = 
                    List.map (fun lcblst ->
                    List.filter (fun record ->
                     Hashtbl.mem new_outer_lcb_tbl (get_abs_lst record)
                    ) lcblst
                    ) new_outer_lcbs 
            in
            if debug_main&&(new_outer_covR=0.) then 
                Printf.printf "remove lcb did not give us any high W lcbs\n%!";
            (*if debug_main then begin
                Printf.printf "\n after remove light lcbs, we have :\n%!";
                Hashtbl.iter (fun key record ->
                    print_lcb record 
                ) new_outer_lcb_tbl;
            end;*)
            if (!outer_old_covR < new_outer_covR) then begin
                if debug_main then 
                    Printf.printf "\n we have a new lcb cov = %f>old one=%f, \
                continue with outer loop again\n%!" new_outer_covR !outer_old_covR;
                any_improvement_outer := true;
                if (new_outer_covR> !maximum_cover_ratio)||(tiny_improvement
                    new_outer_covR !outer_old_covR) then outer_sign := false
                else outer_sign := true;
                outer_old_covR := new_outer_covR;
                outer_lcbs := new_outer_lcbs;
                outer_lcb_tbl := Hashtbl.copy new_outer_lcb_tbl;
                outer_mum_tbl := !inner_mum_tbl;
            end
            else begin
                outer_sign := false;
                if debug_main then 
                Printf.printf "\n no improve on covR (%f>=%f), get out of outer loop\n%!" 
                !outer_old_covR new_outer_covR;
            end;
        end
        else begin
            outer_sign := false;
            if debug_main then
                Printf.printf "no improve in innerloop, no need to redo outerloop again\n%!";
        end;
    done; (*end of outer while loop*)
    (*search_inside_each_lcb !outer_lcb_tbl in_seqarr !maximum_lcb_len max_int
    !outer_mum_tbl !inner_seed2pos_tbl cost_mat use_ukk ;*)
    if !any_improvement_outer = false then 
    (*when outer&inner while did not find any qualified lcb, outer_lcb_tbl still
    * have the lcb from initial function, remove lightW ones*)
        Hashtbl.iter (fun key record ->
            if (is_light_weight_lcb record ) then begin
                Hashtbl.remove !outer_lcb_tbl key;
                outer_lcbs := List.map (fun lcblst ->
                List.filter (fun x -> (get_abs_lst x)<>key ) lcblst
            ) !outer_lcbs;
            end;
        ) !outer_lcb_tbl;
    if debug_main then
        Printf.printf "outer covR = %f, outer_lcb_tbl len=%d\n%!"
        !outer_old_covR (Hashtbl.length !outer_lcb_tbl);
    if (Hashtbl.length !outer_lcb_tbl)=1 then begin
        let lightW = ref false in
        Hashtbl.iter (fun key record ->
            if (is_light_weight_lcb record )
                || (is_low_score_lcb record) then  lightW := true
        ) !outer_lcb_tbl;
        (*Say at the end of loops, we have only one lcb.
        * if it's light weight(W), we remove it.
        * if total coverage is lower than minimum_cover_ratio, we remove it.
        *)
        if !lightW then Hashtbl.clear !outer_lcb_tbl;
        if !outer_old_covR< !minimum_cover_ratio then 
            Hashtbl.clear !outer_lcb_tbl;
    end;
    end;(*end of non trivial case -- when init mum number > 0*)
    let outer_lcb_tbl = !outer_lcb_tbl in
    (*if no qualified lcb are found, just make the whole sequence as one lcb*)
    if (Hashtbl.length outer_lcb_tbl)=0 then begin 
        if debug_main then Printf.printf "we didn't find any qualified lcb\n%!";
        (*let lcb_range_lstlst = get_lcbs_range outer_lcb_tbl in_seqlst true in
        print_lcbs_range lcb_range_lstlst;*)
        let code_list = [[1];[1]] in
        let range0 =  (0,(List.hd in_seq_size_lst)-1,[1],1)
        and range1 =  (0,(List.nth in_seq_size_lst 1)-1,[1],1)
        in
        let m0 = {sequence_NO = 0; left_end = 0; right_end = (List.hd
        in_seq_size_lst)-1; orientation = 1}
        and m1 = {sequence_NO = 1; left_end =0; right_end = (List.nth
        in_seq_size_lst 1)-1; orientation = 1}
        in
        let full_range_lstlst = [ [range0];[range1] ] in
        let single_lcb = 
            {
                seedNOlst = [1];
                range_lst = [m0;m1];
                score = 0;
                ratio = 0.0;
                ref_code = 1;
                avg_range_len = 0;
                alignment = [||]
            }
        in
        Hashtbl.clear outer_lcb_tbl;(*not necessary,outer_lcb_tbl should be empty here*)
        Hashtbl.add outer_lcb_tbl [1] single_lcb;
        let outer_lcbs = [[[1]];[[1]]] in
        outer_lcb_tbl,outer_lcbs,code_list,full_range_lstlst,!outer_mum_tbl,!inner_seed2pos_tbl
    end
    else begin
        let outer_lcbs = !outer_lcbs in
        if debug_main then begin
            Printf.printf 
            "we do find some high W high R lcbs(or the last low RorW lcb), \
            return them to aliMap (final covR=%f)\n%!" !outer_old_covR;
            if debug2 then begin
                Hashtbl.iter (fun key record -> print_lcb record) outer_lcb_tbl;
                Printf.printf "lcb list list is:%!";
                print_lcblst outer_lcbs;
            end;
        end;
        let refcode = ref 1 in
        Hashtbl.iter (fun key record -> 
            refcode := !refcode;
            update_lcb_ref_code outer_lcb_tbl key (!refcode);
            refcode := !refcode + 1;
        ) outer_lcb_tbl;
        let lcb_range_lstlst = get_lcbs_range outer_lcb_tbl in_seq_size_lst in
        let full_range_lstlst = 
            get_full_range_lstlst lcb_range_lstlst in_seq_size_lst in
        if debug2 then begin
            Printf.printf "LCB range lstlst %!";
            print_lcbs_range lcb_range_lstlst;
            Printf.printf "the full range lstlst (includes non-LCB blocks)%!";
            print_lcbs_range full_range_lstlst; 
        end; 
        let seqNO = ref (-1) in
        let code_list = List.map (fun lcblst -> 
            seqNO := !seqNO + 1;
            let tmplst = List.map (fun key ->
                let record = try (Hashtbl.find outer_lcb_tbl (get_abs_lst key)) 
                with |Not_found -> 
                    Printf.printf "cannot find lcb:%!"; 
                    printIntList2 key;
                    failwith "not found in creating code list" in
                let mi_lst = record.range_lst in
                let mi = get_position_by_seqNO mi_lst !seqNO in
                record.ref_code * mi.orientation ) lcblst 
            in
            if debug2 then printIntList2 tmplst;
            tmplst
        ) outer_lcbs in
        (*clear up,may not be since we create seedNO_available_arr at the beginning of this function*)
        seedNO_available_arr := Array.make init_seed_size 1;
        if debug_main then Printf.printf "end of create lcb tbl\n%!";
        (*return lcb table, lcbs, code list and range list*)
        outer_lcb_tbl,outer_lcbs,code_list,full_range_lstlst,!outer_mum_tbl,!inner_seed2pos_tbl
    end


