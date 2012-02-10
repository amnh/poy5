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

(* This implementation of sequence can only handle DNA sequences now, it
* will be easy to extend it if necessary, but the payoff in speed can be 
* considerable, so the best option is, if needed, to implement a second 
* library for extended alphabet sizes. *)
exception Invalid_Argument of string;;
exception Invalid_Sequence of (string * string * int);; 

let () = SadmanOutput.register "Sequence" "$Revision: 2871 $"

external register : unit -> unit = "seq_CAML_register"
let () = register ()

let ndebug = true

type s

external c_reverse : (s -> s -> unit) = "seq_CAML_reverse"

external reverse_ip : (s -> unit) = "seq_CAML_reverse_ip"

external capacity : (s -> int) = "seq_CAML_get_cap"

external create : int -> s = "seq_CAML_create"

external copy : (s -> s -> unit) = "seq_CAML_copy"

external length : (s -> int) = "seq_CAML_length"

external get : (s -> int -> int) = "seq_CAML_get"

external c_median_2 :
    s -> s -> Cost_matrix.Two_D.m -> s -> unit = "seq_CAML_median_2_no_gaps"

external c_median_2_with_gaps :
    s -> s -> Cost_matrix.Two_D.m -> s -> unit = "seq_CAML_median_2_with_gaps"

let get s x =
    assert (x >= 0);
    if (x>= (length s)) then
        Printf.printf "error in sequence.get: x=%d while seqlen=%d \n%!" x (length s);
    assert (x < length s);
    get s x

external count : (int -> int -> s -> int ) = "seq_CAML_count";;

external set : (s -> int -> int -> unit) = "seq_CAML_set";;

external prepend : s -> int -> unit = "seq_CAML_prepend";;


let check_level cm =
    let level = Cost_matrix.Two_D.get_level cm in
    let size = Cost_matrix.Two_D.get_ori_a_sz cm in
    (level>1) && (level<=size)


let make_empty a =
    let s = create 1 in
    prepend s (Alphabet.get_gap a);
    s


let map_to_array f s = 
    Array.init (length s) (fun x -> f (get s x)) 


let mapi f s =
    let len = length s in
    let new_s = create len in
    for i = len - 1 downto 0 do
        prepend new_s (f (get s i) i);
    done;
    new_s


let map f s =
    let len = length s in
    let new_s = create len in
    for i = len - 1 downto 0 do
        prepend new_s (f (get s i));
    done;
    new_s


let fold f acc s =
    let len = length s in
    let rec folder acc pos =
        if pos < len then 
            folder (f acc (get s pos)) (pos + 1)
        else acc
    in
    folder acc 0


let foldi f acc s =
    let len = length s in
    let rec folder acc pos = 
        if pos < len then
            folder (f acc pos (get s pos)) (pos + 1)
        else acc
    in
    folder acc 0


let fold_right f acc s =
    let len = length s in
    let rec folder acc pos =
        if pos > (-1) then 
            folder (f acc (get s pos)) (pos - 1)
        else acc
    in
    folder acc (len - 1)


let fold_righti f acc s =
    let len = length s in
    let rec folder acc pos =
        if pos > (-1) then
            folder (f acc pos (get s pos)) (pos - 1)
        else acc
    in
    folder acc (len - 1)


let iter f s = 
    let len = length s in
    for i = len - 1 downto 0 do 
        f (get s i);
    done;
    ()


let init f len = 
    let seq = create len in
    for i = len - 1 downto 0 do
        prepend seq (f i);
    done;
    seq


let resize n ns = 
    let v = create ns in
    copy !n v;
    n := v;;


let clone n =
    let sz = length n in
    let res = create sz in
    copy n res;
    res;;


let reverse s1 = 
    let sp = create (length s1) in
    c_reverse s1 sp;
    sp;;


let safe_reverse x = 
    let gap = get x 0 in
    let y = create (length x) in
    for i = 1 to (length x) - 1 do
        prepend y (get x i);
    done;
    prepend y gap;
    y


let rec aux_of_list seq l it = match l with 
    | h :: t -> 
        set !seq it h;
        aux_of_list seq t (it + 1);
    | [] -> seq;;


let of_list l =
    let length = List.length l in
    let seq = create length in
    aux_of_list (ref seq) l 0;;


let remove_gaps2 s2' cm =
    let remove_gaps gap seq base = 
        if base <> gap then 
            let () = prepend seq base in
            seq
        else seq
    in
    let res = 
        fold_right (remove_gaps (Cost_matrix.Two_D.gap cm)) 
                   (create (length s2'))
                   s2'
    in
    prepend res (Cost_matrix.Two_D.gap cm);
    res


let remove_gaps s2' ?(gapcode=Alphabet.gap) prependgap =
    let remove_gap gap seq base = 
        if base <> gap then 
            let () = prepend seq base in
            seq
        else seq
    in
    let res = 
        fold_right (remove_gap gapcode) 
                   (create (length s2'))
                   s2'
    in
    if prependgap then prepend res gapcode;
    res


let is_empty seq gap =
    let length = length seq in
    let rec check p =
        if p = length then true
        else 
            if gap <> get seq p then false
            else check (p + 1) 
    in
    check 0


let prepend_gap s m =
    let () = prepend s (Cost_matrix.Two_D.gap m) in
    s


(* Why isn't ocaml making this function a loop rather than a recursive call? *)
let rec aux_parse str seq alph = function 
    | -1 -> seq
    | it -> 
        let rec get_base cur_str cur_it =
            let ch = Char.escaped str.[cur_it] in
            let new_str = ch ^ cur_str in 
            begin try
                let code = Alphabet.match_base new_str alph in
                prepend seq code;
                cur_it
            with | Alphabet.Illegal_Character c -> 
                match cur_it with
                    | 0 -> raise (Invalid_Sequence (str, c, it))
                    | _ -> get_base new_str (cur_it - 1)
            end
        in
        let new_it = get_base "" it in 
        aux_parse str seq alph (new_it - 1) 
              

let aux_parse_imp str seq alph =
    let len = String.length str in
    for i = len - 1 downto 0 do
        try
            let v = Char.escaped str.[i] in
            prepend seq (Alphabet.match_base v alph);
        with
        | Alphabet.Illegal_Character c ->
            raise (Invalid_Sequence (str, c, i));
    done;
    seq


let of_string str alph =
    let len = String.length str in
    let seq = create len in
    aux_parse str seq alph (len - 1)


let of_string_ls str_ls alph = 
    let len = List.length str_ls in 
    let seq = create len in
    let rec aux_parse_ls str_ls alph = match str_ls with
        | [] -> seq
        | hd::tl -> 
            begin try
                let code = Alphabet.match_base hd alph in
                prepend seq code;              
            with | Alphabet.Illegal_Character c -> 
                failwith ("The character state " ^ c ^ " does not exists!!!")
            end;
            aux_parse_ls tl alph
    in
    aux_parse_ls str_ls alph 


let to_string seq alph =
    let len = length seq in
    let seq,len = 
        if len=0 then
            (make_empty alph),1    
        else
            seq,len
    in
    let b = Buffer.create len in
    for i = 0 to len - 1 do
        Buffer.add_string b (Alphabet.match_code (get seq i) alph);
    done;
    Buffer.contents b


let to_formater seq alph =
    let len = length seq in
    let b = Buffer.create len in
    for i = 0 to len - 1 do
        let tmp = get seq i in
        if(tmp=0) then
            Printf.printf "Warning: try to get code 0 in alphabet map, check the rest of sequence code: \n%!";
        let str = 
            let v = Alphabet.match_code (get seq i) alph in
            if v = "@" then "@@"
            else if v = "%" then "%%"
            else v
        in
        Buffer.add_string b str
    done;
    Buffer.contents b

let print chn seq alph =
    let str = to_string seq alph in
    Pervasives.output_string chn str;;

(** [printDNA seq] prints the DNA sequence [seq] into stdout *)
let printDNA seq = 
    print stdout seq Alphabet.nucleotides; 
    print_newline ()

let printseqcode seq =
    let len = length seq in
    Printf.printf "[%!";
    for i = 0 to (len-1) do  
        Printf.printf "%d,%!" (get seq i) 
    done;
    Printf.printf "]\n%!"


let copy_from_in x y z u =
    let to_copy = get x z in
    set y u to_copy


let concat x = 
    let len = List.fold_left (fun x y -> x + length y) 0 x in
    let ns = init (fun _x -> 0) len in 
            
	let pos = ref 0 in    
    let copier x = 
        let len = length x in
       for i = 0 to len - 1 do
			copy_from_in x ns i !pos;
            Pervasives.incr pos;            
        done;
    in
    List.iter (copier) x;
    ns


let to_array s =
    Array.init (length s) (fun x -> get s x)


(** [of_array code_arr] converts the code array [code_arr]
* into a sequence *)
let of_array code_arr = 
    let len = Array.length code_arr in 
    init (fun idx -> code_arr.(idx)) len  


let sub s st len =
    init (fun x -> get s (x + st)) len


let sub_ignore_gap ?(gap=Alphabet.gap) s st len = 
    let debug = false in
    let idx = ref st in
    let count = ref 0 in
    let slen = length s in
    let acclst = ref [] in
    if debug then Printf.printf "sub ignore gap(=%d), slen=%d,start idx=%d,len without gap = %d -> %!" 
    gap slen st len;
    while (!count < len) do
        assert( !idx < slen);
        let add = get s !idx in
        acclst := !acclst @ [add];
        idx := !idx +1;
        if add<>gap then 
            count := !count + 1
        else 
            if debug then Printf.printf "hit a gap;%!"
    done;
    if debug then Printf.printf "end of sub ignore gap,next idx set to %d\n%!" !idx;
    of_array (Array.of_list !acclst) , !idx

let length_without_gap ?(gap=Alphabet.gap) s =
    let count = ref 0 in
    let slen = length s in
    for i = 0 to slen-1 do
        let tmp = get s i in
        if tmp<>gap then count := !count + 1
        else ()
        ;
    done;
    !count

let prepend_char seq element =  
    let len = length seq in
    let ext_seq = init
        (fun index ->
             match  index with
             | 0  -> element
             | _ -> get seq (index - 1)
        )  (len + 1) 
    in
    ext_seq


let del_first_char seq =
    let len = length seq in  
    sub seq 1 (len - 1)


let compare a b = 
    let la = length a
    and lb = length b in
    let lc = min la lb in
    let rec comparator cnt =
        if cnt < lc then begin
            let ca = get a cnt 
            and cb = get b cnt in
            match ca - cb with
            | 0 -> comparator (cnt + 1)
            | v -> v
        end else 0
    in
    match comparator 0 with
    | 0 -> 
            if la = lb then 0
            else la - lb
    | v -> v

let median_2_with_gaps s1 s2 c =
    let sz1 = length s1 
    and sz2 = length s2 in
    if (sz1 = sz2) then begin
        let sp = create sz1 in
        c_median_2_with_gaps s1 s2 c sp;
        sp
    end else 
        raise (Invalid_Argument "The size of the sequences is not the same.")

let median_2 s1 s2 c =
    let sz1 = length s1 
    and sz2 = length s2 in
    if (sz1 = sz2) then begin
        let sp = create (sz1 + 1) in
        c_median_2 s1 s2 c sp;
        sp
    end else 
        raise 
        (Invalid_Argument "The size of the sequences is not the same.")

module SerSequence = struct
    type v = int
    type ser = s
    let equal (a:v) (b:v) = a = b
    let get_pos = get
    let to_int v = v
    let of_int v = v
    let length = length
    let sub = sub
end

(* A finger print module for sequences *)
module SeqFp = FingerPrint.Make (SerSequence)
module SeqMFp = FingerPrint.MakeMult (SerSequence)

module Anchors = struct
    (* Builds a map containing the last substring of length [len] with a unique 
    * hash value. It's a quick and dirty index generator. The hashes are
    * returned in the original order of the sequence 0, 1, 2, 3, ... *)
    let gen_subhash s p len alph ffp ifp = 
        let a = (Alphabet.size alph) in
        let inital = ffp s p 0 len a in
        let rec builder it max acc =
            if it = max then acc 
            else begin match acc with
            | (prev, _) :: _ -> 
                    let next = ifp s p (it - 1) len prev a in
                    builder (it + 1) max ((next, it) :: acc)
            | _ -> failwith "Unexpected"
            end
        in
        let res = builder 1 ((length s) - len) [(inital, 0)] in
        List.rev res

    (* Generates the list of shared values between two sequences [a] and [b] 
     * with corresponding hash lists [a_hs] and [b_hs] as generated by
     * [gen_subhash]. The [a_mp] parameter is a mapping of the association list
     * [a_hs] and mapmem and mapfind are functions to perform the mem and find
     * operations in the map [a_mp]. The function returns an association list of 
     * [(x, y)] where [x] is the position in the sequence [b] that has the
     * same hash value of the position [y] in [a]. *)
    let gen_shared_subseqs p a_hs b_hs a_mp mapmem mapfind = 
        let builder acc (it, cnt) =
            if mapmem it a_mp then begin
                let pos = mapfind it a_mp in 
                (cnt, pos) :: acc
            end else acc
        in
        List.fold_left builder [] b_hs

    let gen_make_map_of_subhash a_hs addmap emptymap =
        let a_mp, _ = List.fold_left 
            (fun (acc, cnt) x -> addmap x cnt acc, cnt + 1)
            (emptymap, 0) a_hs
        in
        a_mp

    let subhash s p len alpha =
        gen_subhash s p len alpha SeqFp.fp SeqFp.fp_incremental

    let shared_subseqs p a_hs b_hs a_mp =
        gen_shared_subseqs p a_hs b_hs a_mp 
        All_sets.IntegerMap.mem All_sets.IntegerMap.find 

    let subhash_mult s p len alpha = 
        gen_subhash s p len alpha SeqMFp.fp SeqMFp.fp_incremental

    let shared_subseqs_mult p a_hs b_hs a_mp l alph =
        gen_shared_subseqs p a_hs b_hs a_mp 
        All_sets.IntegerListMap.mem All_sets.IntegerListMap.find

    let consistent_groups lst = 
        let rec add_to_proper_group acc ((x, y) as res) = 
            match acc with
            | (((_, hb) :: _) as h) :: t when y < hb -> 
                    (res :: h) :: t
            | h :: t -> 
                    h :: (add_to_proper_group t res)
            | [] -> [[res]]
        in
        List.fold_left add_to_proper_group [] lst

    let consecutive_matches lst = 
        let prepender (acc, (a, b), (x, y)) (c, d) =
            ((((x, a), (y, b)) :: acc), (c, d), (c, d))
        in
        let classifier ((acc, (a, b), (x, y)) as res) (c, d) =
            if (a = (c - 1)) && (b = (d - 1)) then
                (acc, (c, d), (x, y))
            else prepender res (c, d)
        in
        match lst with
        | hd :: tl -> 
                let res = List.fold_left classifier ([], hd, hd) tl in
                let res, _, _ = prepender res (max_int, max_int) in
                res
        | [] -> []

end



module Align = struct

    let code2list = Cost_matrix.Two_D.combcode_to_comblist
    let list2code = Cost_matrix.Two_D.comblist_to_combcode
    let print_intlist = Cost_matrix.Two_D.print_intlist
    let gap_filter_for_combcode = Cost_matrix.Two_D.gap_filter_for_combcode
    let get_level = Cost_matrix.Two_D.get_level
    let get_ori_a_sz = Cost_matrix.Two_D.get_ori_a_sz
    let get_full_a_sz = Cost_matrix.Two_D.get_map_sz
    let clear_duplication_in_list = Cost_matrix.Two_D.clear_duplication_in_list
    let calc_comb_with_gap = Cost_matrix.Two_D.calc_num_of_comb_with_gap

    external c_max_cost_2 : s -> s -> Cost_matrix.Two_D.m -> int = "algn_CAML_worst_2"

    (*[c_cost_2] and [c_cost_2_limit] return the alignment cost, the result of
    * alignment is stored in Matrix.m, we can retrieve it with
    * [extract_edited_2] *)
    external c_cost_2 :
        s -> s -> Cost_matrix.Two_D.m -> Matrix.m -> int -> int =
            "algn_CAML_simple_2"

    external c_cost_2_limit :
        s -> s -> Cost_matrix.Two_D.m -> Matrix.m -> int -> int -> int ->
            int -> int -> int -> int = "algn_CAML_limit_2_bc" "algn_CAML_limit_2"

    (*[cost_2_affine] and [c_align_affine_3] for affine alignment,
    * besides cost, [c_align_affine_3] also give
    * us alignment, which is called traceback in algn.c *)
    external c_align_affine_3 : s -> s -> Cost_matrix.Two_D.m -> Matrix.m -> 
        s -> s -> s -> s -> int = "algn_CAML_align_affine_3_bc"
        "algn_CAML_align_affine_3"

    external cost_2_affine : s -> s -> Cost_matrix.Two_D.m -> Matrix.m -> int =
        "algn_CAML_cost_affine_3" 

    let align_affine_3 si sj cm =
        let len = length si + length sj + 2 in
        let resi = create len
        and resj = create len 
        and median = create len 
        and medianwg = create len in
        let cost = 
            c_align_affine_3 si sj cm Matrix.default resi resj median
            medianwg in
        (median, resi, resj, cost, medianwg)
    
    let max_cost_2 a b c =
        let gap = Cost_matrix.Two_D.gap c in
        if is_empty a gap || is_empty b gap then 0
        else begin
            assert ((length a) = (length b));
            c_max_cost_2 a b c
        end

    type gap_side = NoGap | AGap | BGap

    let verify_cost_2 expected_cost a b c =
        let gap = Cost_matrix.Two_D.gap c in
        let do_combine = 1 = Cost_matrix.Two_D.combine c in
        if is_empty a gap || is_empty b gap then 0
        else begin
            let go = Cost_matrix.Two_D.gap_opening c in
            let lena = length a
            and alph = Cost_matrix.Two_D.alphabet_size c in
            assert (lena = (length b));
            let rec verify_cost_recursively a b pos cur_best gap_block_side acc =
                if cur_best < acc then cur_best
                else if pos = lena then
                    if cur_best <= expected_cost then raise Exit
                    else if cur_best <= acc then cur_best
                    else acc
                else
                    let get_list x =
                        let base = get x pos in
                        if do_combine then
                            List.sort (fun a b -> b - a)
                                      (BitSet.Int.list_of_packed_max base alph)
                        else
                            [base]
                    in
                    let put_shared_first lst1 lst2 =
                        let shared =
                            List.filter (fun x -> List.exists (fun y -> y = x) lst1) lst2
                        in
                        match shared with
                        | [] -> lst1, lst2
                        | _ ->
                            let filter lst = 
                                shared @ (List.filter (fun x -> not
                                (List.exists (fun y -> y = x) shared)) lst)
                            in
                            (filter lst1), (filter lst2)
                    in
                    let basea_list = get_list a
                    and baseb_list = get_list b in
                    let basea_list, base_blist = put_shared_first basea_list
                    baseb_list in
                    let process_pairwise basea baseb (acc, gap_block_side) =
                        let cost = Cost_matrix.Two_D.cost basea baseb c in
                        let acc = cost + acc in
                        match gap_block_side with
                        | NoGap -> 
                                if basea = gap then
                                    if baseb = gap then 
                                        acc, NoGap
                                    else 
                                        (acc + go), AGap
                                else if baseb = gap then
                                    acc + go, BGap
                                else 
                                    acc, NoGap
                        | AGap ->
                                if basea = gap then
                                    acc, AGap
                                else if baseb = gap then
                                    acc + go, BGap
                                else acc, NoGap
                        | BGap ->
                                if baseb = gap then
                                    acc, BGap
                                else if basea = gap then
                                     acc + go, AGap
                                else acc, NoGap
                    in
                    let cur_best, _ =
                        List.fold_left (fun acc basea ->
                            List.fold_left (fun (cur_best, acc) baseb ->
                                let nacc, gap_block_side = 
                                    process_pairwise basea baseb 
                                    (acc, gap_block_side) 
                                in
                                (verify_cost_recursively a b (pos + 1) cur_best
                                gap_block_side nacc), acc) 
                            acc baseb_list) 
                        (cur_best, acc) basea_list
                    in
                    cur_best
            in
            try verify_cost_recursively a b 0 max_int NoGap 0 with
            | Exit -> expected_cost
        end

    let default_length_calculation w h c d =
        let default_length = 
            min d (2 * (abs (d - c))) in
        let w = 
            match w with
            | Some v -> v
            | None -> default_length
        and h = 
            match h with
            | Some v -> v
            | None -> default_length
        in
        w, h

    let counter = ref 0

    let cost_2_limit ?w ?h s1 s2 m1 m2 a b c d =
        let w, h = default_length_calculation w h c d in
        let w = if w > c then c else w
        and h = if h > d then d else h in
        if not ndebug then 
            Printf.printf "I will request positions %d and %d in sequences \
            of length %d and %d respectively.\n%!" a b (length s1) (length s2);
        let bs1 = get s1 a
        and bs2 = get s2 b 
        and gap = Cost_matrix.Two_D.gap m1 in
        set s1 a gap;
        set s2 b gap;
        let res = 
            if c >= d then 
                c_cost_2_limit s1 s2 m1 m2 w h a b c d
            else 
                c_cost_2_limit s2 s1 m1 m2 h w b a d c
        in
        set s1 a bs1;
        set s2 b bs2;
        res

    let cost_2_stitched p w alph s1 s2 m1 m2 =
        if not ndebug then 
            Printf.printf "Going to align two sequences using stitched \
            algorithm.\n%!";
        let s1_lst = Anchors.subhash_mult s1 p w alph
        and s2_lst = Anchors.subhash_mult s2 p w alph in
        if not ndebug then begin
            Printf.printf "\nSequence s1 hashed values : ";
            List.iter (fun (x, y) -> 
                match x with
                | [x1; x2] -> Printf.printf "(%d:%d, %d) " x1 x2 y
                | _ -> failwith "Unexpected") s1_lst;
            Printf.printf "\nSequence s2 hashed values : ";
            List.iter  (fun (x, y) -> 
                match x with
                | [x1; x2] -> Printf.printf "(%d:%d, %d) " x1 x2 y
                | _ -> failwith "Unexpected") s2_lst;
        end;
        let generate_map lst = List.fold_left (fun acc (x, y) ->
            All_sets.IntegerListMap.add x y acc)
            All_sets.IntegerListMap.empty lst
        in
        let s1_map = generate_map s1_lst in
        let res = Anchors.shared_subseqs_mult p s1_lst s2_lst s1_map w alph in
        let res = Anchors.consistent_groups res in
        let res, _ = List.fold_left (fun ((x, len) as acc) y -> 
            let ylen = List.length y in
            if len >= ylen then acc
            else (y, ylen)) ([], 0) res
        in
        let lst = Anchors.consecutive_matches res in
        if not ndebug then begin
            Printf.printf "\nSequences consecutive matches: ";
            List.iter (fun ((x, y), (u, v)) ->
                Printf.printf "((%d, %d), (%d, %d)) %!" x y u v) lst;
            print_newline ();
            Printf.printf "The total number of shared subsequences is %d\n%!" 
            (List.length lst);
        end;
        let distance_of_two_subseqs 
            (((a, b), (c, d)), cost) (((e, f), (g, h)) as res) =
                if not ndebug then begin
                    Printf.printf "The current cost is %d\n%!" cost;
                    Printf.printf "First shared position is ((%d, %d), \
                    (%d, %d)) for s2, s1 and the second is \
                    ((%d, %d), (%d, %d)) for s2, s1\n%!" 
                    a b c d e f g h
                end;
                let lens2 = a - f
                and lens1 = c - h in
                if not ndebug then begin
                    Printf.printf "The length of the sequences is %d for s1 \
                    and %d for s2\n%!" lens1 lens2;
                end;
                res, cost + (cost_2_limit s1 s2 m1 m2 h f lens1 lens2) 
        in
        let lens1 = length s1 
        and lens2 = length s2 in
        let all_distances_excepting_first_segment = 
            let last_segment = ((lens2, lens2), (lens1, lens1)) in
            List.fold_left distance_of_two_subseqs (last_segment, 0) lst 
        in
        let _, total = 
            let first_segment = ((0, 0), (0, 0)) in
            distance_of_two_subseqs all_distances_excepting_first_segment 
            first_segment in
        total

    let count_gaps s m = 
        let gap = Cost_matrix.Two_D.gap m in
        let uselevel = check_level m in
        if uselevel then
            let level = Cost_matrix.Two_D.get_level m in
            let orisize = Cost_matrix.Two_D.get_ori_a_sz m in
            let totalnum = Cost_matrix.Two_D.get_map_sz m in
            let numwithgap = 
                Cost_matrix.Two_D.calc_num_of_comb_with_gap orisize level in
            count gap (totalnum-numwithgap+1) s 
        else
           count gap 0 s 

    let cost_2 ?deltaw s1 s2 m1 m2 =
        let deltaw_calc s1len s2len = 
            let dif = s1len - s2len in
            let lower_limit = int_of_float ((float_of_int s1len) *.  0.10) in
            match deltaw with 
            | None -> 
                    if dif < lower_limit then (lower_limit / 2)
                    else 2
            | Some v -> 
                    if dif < lower_limit then lower_limit
                    else v
        in
        let ls1 = length s1 and ls2 = length s2 in
        let s1,ls1 = 
          if ls1=0 then (prepend_gap s1 m1),1
          else s1,ls1
        in
        let s2,ls2 =
            if ls2=0 then (prepend_gap s2 m1),1 
            else s2,ls2
        in
        (*assert (ls1 <> 0);  assert (ls2 <> 0);*)
        let g1 = count_gaps s1 m1
        and g2 = count_gaps s2 m1 in
        let gaps = max g1 g2 in
        (*note that we pass longer seq as seq1 to the cside, this is different
        * in newkkonen/affine alignement/likelihood*)
        if ls1 >= ls2 then
            let deltaw = gaps + deltaw_calc ls1 ls2 in
            c_cost_2 s1 s2 m1 m2 deltaw
        else 
            let deltaw = gaps + deltaw_calc ls2 ls1 in
            c_cost_2 s2 s1 m1 m2 deltaw


    let cost_2 ?deltaw s1 s2 m1 m2 =
        match Cost_matrix.Two_D.affine m1 with
        | Cost_matrix.Affine _ ->
                cost_2_affine s1 s2 m1 m2
        | _ -> 
                match deltaw with 
                | None -> cost_2 s1 s2 m1 m2
                | Some deltaw -> cost_2 ~deltaw s1 s2 m1 m2

    external myers : s -> s -> int = "algn_CAML_myers"

    external cost_3 :
        s -> s -> s -> Cost_matrix.Three_D.m -> Matrix.m -> int -> int =
            "algn_CAML_simple_3_bc" "algn_CAML_simple_3"

    external extract_edited_2 :
        s -> s -> s -> s -> Matrix.m -> Cost_matrix.Two_D.m -> bool -> unit =
            "algn_CAML_backtrace_2d_bc" "algn_CAML_backtrace_2d"

    external extract_edited_2_limit :
        s -> s -> s -> s -> Matrix.m -> Cost_matrix.Two_D.m -> 
            int -> int -> int -> int -> bool -> unit =
            "algn_CAML_backtrace_2d_limit_bc" "algn_CAML_backtrace_2d_limit"

    external extract_edited_3 :
        s -> s -> s -> s -> s -> s -> Matrix.m -> Cost_matrix.Three_D.m -> 
            unit = "algn_CAML_backtrack_3d_bc" "algn_CAML_backtrack_3d"

    external c_align_2 :
        s -> s -> Cost_matrix.Two_D.m -> Matrix.m -> s -> s -> int -> int ->
            bool =
            "algn_CAML_align_2d_bc" "algn_CAML_align_2d"

    external c_align_3 :
        s -> s -> s -> Cost_matrix.Three_D.m -> Matrix.m -> s -> s -> s
        -> int = "algn_CAML_align_3d_bc" "algn_CAML_align_3d"

    external c_median_3 :
        s -> s -> s -> Cost_matrix.Three_D.m -> s -> unit =
            "algn_CAML_median_3"

    external c_print_bcktrack :
        s -> s -> Matrix.m -> unit = "algn_CAML_print_bcktrck"

    let print_backtrack a b m =
        if length a >= length b then c_print_bcktrack a b m
        else c_print_bcktrack b a m

    external c_create_backtrack : 
        int -> int -> Matrix.m -> int array array -> unit =
            "algn_CAML_create_backtrack"

    let make_backtrack a b m =
        let la = length a 
        and lb = length b in
        let comparison = la >= lb in
        let matrix = 
            if comparison then Array.make_matrix la lb 0 
            else Array.make_matrix lb la 0 
        in
        if comparison then c_create_backtrack la lb m matrix
        else c_create_backtrack lb la m matrix;
        matrix

    let rec check_one x y v mtx =
        if v = 0 then 0
        else aux_count_paths x y mtx

    and aux_count_paths x y mtx =
        if x = 0 || y = 0 then 1
        else
            let v = mtx.(x).(y) in
            let a = v land 1
            and b = v land 2
            and c = v land 4 in
            let first = (check_one (x - 1) (y - 1) a mtx) in
            if first > 10000 then first
            else
                let first = first + (check_one x (y - 1) b mtx) in
                if first > 10000 then first 
                else first + (check_one (x - 1) y c mtx)

    let count_paths mtx = 
        let x = Array.length mtx in
        if x = 0 then 1
        else
            let y = Array.length mtx.(0) in
            if y = 0 || y = 1 then 1
            else aux_count_paths (x - 1) (y - 1) mtx

    let create_edited_2 s1 s2 tm c =
        let debug = false in
        let sz1 = length s1
        and sz2 = length s2 in
        let s1p = create (sz1 + sz2)
        and s2p = create (sz1 + sz2) in
        let size_compared = sz1 >= sz2 in
        if debug then 
            Printf.printf "create_edited_2,sz1=%d,sz2=%d,size_compared=%b\n%!"
            sz1 sz2 size_compared;
        if size_compared then 
            extract_edited_2 s1 s2 s1p s2p tm c size_compared
        else 
            extract_edited_2 s2 s1 s2p s1p tm c size_compared;
        s1p, s2p

    let create_edited_2_limit s1 s2 tm c st1 st2 len1 len2 =
        let s1p = create (len1 + len2)
        and s2p = create (len1 + len2) in
        let size_compared = len1 >= len2 in
        Printf.printf "create_edited_2_limit,sz1=%d,sz2=%d,size_compared=%b\n%!"
        len1 len2 size_compared;
        if size_compared then 
            extract_edited_2_limit s1 s2 s1p s2p tm c st1 st2 len1 len2
            size_compared
        else 
            extract_edited_2_limit s2 s1 s2p s1p tm c st2 st1 len2 len1
            size_compared;
        s1p, s2p

    let create_edited_3 s1 s2 s3 tm cm =
        let sz1 = length s1
        and sz2 = length s2
        and sz3 = length s3 in
        let total_len = sz1 + sz2 + sz3 in
        let s1p = create total_len
        and s2p = create total_len 
        and s3p = create total_len in
        extract_edited_3 s1 s2 s3 s1p s2p s3p tm cm;
        s1p, s2p, s3p


    let align_2 ?(first_gap=true) s1 s2 c m =
        let cmp s1 s2 =
            match Cost_matrix.Two_D.affine c with
            | Cost_matrix.Affine _ ->
                    let _, s1p, s2p, tc, _ = align_affine_3 s1 s2 c in
                    s1p, s2p, tc
            | _ ->
                    let tc = cost_2 s1 s2 c m in   
                    let s1p, s2p = create_edited_2 s1 s2 m c in   
                    s1p, s2p, tc   
        in 
        match first_gap with
        | true -> cmp s1 s2 
        | false ->              
              let gap = Cost_matrix.Two_D.gap c in 
              let s1 = prepend_char s1 gap in 
              let s2 = prepend_char s2 gap in 
              let s1p, s2p, tc = cmp s1 s2 in
              let s1p = del_first_char s1p in 
              let s2p = del_first_char s2p in 
              s1p, s2p, tc

    let align_3 ?(first_gap = true) s1 s2 s3 c m =
        let align s1 s2 s3 =
            let sz1 = length s1
            and sz2 = length s2 
            and sz3 = length s3 in
            let sum = sz1 + sz2 + sz3 in
            let s1p = create sum
            and s2p = create sum  
            and s3p = create sum in 
            let c = c_align_3 s1 s2 s3 c m s1p s2p s3p in
            s1p, s2p, s3p, c
        in 
        match first_gap with
        | true -> align s1 s2 s3
        | false ->
              let gap = Cost_matrix.Three_D.gap c in 
              let s1 = prepend_char s1 gap in 
              let s2 = prepend_char s2 gap in 
              let s3 = prepend_char s3 gap in 
              let s1p, s2p, s3p, tc = align s1 s2 s3 in 
              let s1p = del_first_char s1p in 
              let s2p = del_first_char s2p in 
              let s3p = del_first_char s3p in 
              s1p, s2p, s3p, tc


    

    external c_ancestor_2 : 
        s -> s -> Cost_matrix.Two_D.m -> s -> unit = "algn_CAML_ancestor_2" 

    let ancestor_2 s1 s2 c =
        let sz1 = length s1 
        and sz2 = length s2 in
        if (sz1 = sz2) then begin
            let sp = create (sz2 + 1) in
            c_ancestor_2 s1 s2 c sp;
            sp
        end else 
            raise 
            (Invalid_Argument "The size of the sequences is not the same.")

    let median_3 s1 s2 s3 c =
        let sz1 = length s1
        and sz2 = length s2 
        and sz3 = length s3 in 
        if (sz1 = sz2) && (sz2 = sz3) then begin
            (* Leave space for the initial gap *)
            let sp = create (sz1 + 1) in 
            c_median_3 s1 s2 s3 c sp;
            sp;
        end else 
            raise 
            (Invalid_Argument "The size of the sequences is not the same.")

    let full_median_2 a b cm m = 
        match Cost_matrix.Two_D.affine cm with
        | Cost_matrix.Affine _ ->
                let m, _, _, _, _ = align_affine_3 a b cm in
                m
        | _ ->
                let a, b, _ = align_2 a b cm m in
                median_2 a b cm

    let full_median_3 a b c cm m =
        let a, b, c, _ = align_3 a b c cm m in
        median_3 a b c cm

    external c_union : s -> s -> s -> unit = "algn_CAML_union"
 
    let union a b cm =
        let lena = length a and lenb = length b in
        assert(lena == lenb);
        let len = lena in
        let res = create (len + 1) in
        let level = get_level cm in
        if (check_level cm) then begin
            let break_union newlist listitem =
            (code2list listitem cm) @ newlist
            in
            for i = len-1 downto 0 do
                let ai = get a i and bi = get b i in
                let alst = code2list ai cm and blst = code2list bi cm in
                let alstlen = (List.length alst) and
                blstlen = (List.length blst) in
                assert (alstlen<=level); assert(blstlen<=level);
                let union_ai_bi =
                    let tmplst = List.fold_left break_union [] ([ai]@[bi]) in
                    let newlist = List.rev (clear_duplication_in_list tmplst) in
                    if (List.length newlist)<= level then
                        list2code ([ai]@[bi]) cm
                    else
                        ai
                    (*if(alstlen+blstlen)<=level then
                        list2code ([ai]@[bi]) cm 
                    else
                        ai*)
                in
                prepend res union_ai_bi
            done;
        end
        else 
            c_union a b res
        ;
        res
        

    let closest s1 s2 cm m =
        let uselevel = check_level cm in
        if is_empty s2 (Cost_matrix.Two_D.gap cm) then
            s2, 0
        else
        let (s, _) as res = 
        if 0 = Cost_matrix.Two_D.combine cm then 
            (* We always have just one sequence per type s *)
            let s1', s2', cst = align_2 s1 s2 cm m in
            let get_closest v i =
                let v' = get s1' i in
                let all = Cost_matrix.Two_D.get_all_elements cm in
                if v = all then 
                    if v' = all then 1 (* We pick one, any will be fine *)
                    else v'
                else v
            in
            (*Printf.printf "; s1' and s2' = %!"; printseqcode s1'; printseqcode s2';
            Printf.printf "end of sequence.ml \n%!";*)
            remove_gaps2 (mapi get_closest s2') cm, cst
        else
            let s1', s2', comp =
                if 0 = compare s1 s2 then
                    (* We remove all the gaps if we are using combinations *)
                    if 0 = Cost_matrix.Two_D.combine cm then 
                        s1, s2, true
                    else
                        if uselevel then 
                            let gap_filter x =
                                gap_filter_for_combcode x cm
                            in
                            let news1 = map gap_filter s1 and
                            news2 = map gap_filter s2 in
                            news1, news2, true
                        else
                            let mask = lnot (Cost_matrix.Two_D.gap cm) in
                            mapi (fun x pos -> 
                                if pos > 0 then x land mask else x) s1, 
                            mapi (fun x pos -> 
                                if pos > 0 then x land mask else x) s2, 
                            true
                else
                    let s1', s2', _ = align_2 s1 s2 cm m in
                    s1', s2', false
            in
            let s2' =
                let s2' = 
                    let get_closest v i =
                        let v' = get s1' i in
                        Cost_matrix.Two_D.get_closest cm v' v 
                    in
                    mapi get_closest s2' 
                in
                (*let _ =
                Printf.printf ", s2' = %!"; printseqcode s2'
                in*)
                remove_gaps2 s2' cm
            in
            (*Printf.printf ", s2' = %!"; printseqcode s2';
            Printf.printf "end of sequence.ml \n%!";*)
            (* We must recalculate the distance between sequences because the
            * set ditance calculation is an upper bound in the affine gap cost
            * model *)
            if comp then 
                s2', 0
            else
                s2', cost_2 s1 s2' cm m
        in
        res


    let recost a b cm =
        assert (length a = length b);
        let len = length a 
        and gap = Cost_matrix.Two_D.gap cm 
        and gap_opening =
            match Cost_matrix.Two_D.affine cm with
            | Cost_matrix.Affine x -> x
            | _ -> 0
        in
        if 1 = Cost_matrix.Two_D.combine cm then
            let rec process_cost it cost is_gap_block =
                if it = len then cost
                else 
                    let ba = get a it
                    and bb = get b it in
                    let c = Cost_matrix.Two_D.cost ba bb cm in
                    if not is_gap_block && ((0 <> ba land gap) || (0 <> bb land
                    gap)) then
                        process_cost (it + 1) (cost + c + gap_opening) true
                    else if (0 = ba land gap) && (0 = bb land gap) then
                        process_cost (it + 1) (cost + c) false
                    else process_cost (it + 1) (cost + c) true
            in
            process_cost 0 0 false
        else
            let rec process_cost it cost is_gap_block =
                if it = len then cost
                else
                    let ba = get a it
                    and bb = get b it in
                    let c = Cost_matrix.Two_D.cost ba bb cm in
                    if not is_gap_block && ((ba = gap) || (bb = gap)) then
                        process_cost (it + 1) (cost + c + gap_opening) true
                    else if (ba <> gap) && (bb <> gap) then
                        process_cost (it + 1) (cost + c) false
                    else process_cost (it + 1) (cost + c) true
            in
            process_cost 0 0 false

    external powell_3d : s -> s -> s -> s -> s -> s -> int -> int -> int -> int =
        "powell_3D_align_bc" "powell_3D_align"

    let align_3_powell s1 s2 s3 mm go ge =
        let sz1 = length s1
        and sz2 = length s2 
        and sz3 = length s3 in
        let sum = sz1 + sz2 + sz3 in
        let s1p = create sum
        and s2p = create sum  
        and s3p = create sum in 
        let c = powell_3d s1 s2 s3 s1p s2p s3p mm go ge in
        s1p, s2p, s3p, c

    let align_3_powell_inter s1 s2 s3 cm cm3 =
        let gap = Cost_matrix.Two_D.gap cm in
        let mm, go, ge =
            match Cost_matrix.Two_D.affine cm with
            | Cost_matrix.Linnear
            | Cost_matrix.No_Alignment -> 
                    Cost_matrix.Two_D.cost 1 2 cm,
                    0,
                    Cost_matrix.Two_D.cost 1 16 cm
            | Cost_matrix.Affine x ->
                    Cost_matrix.Two_D.cost 1 2 cm,
                    x,
                    Cost_matrix.Two_D.cost 1 16 cm
        in
        let a, b, c, cost = align_3_powell s1 s2 s3 mm go ge in
        let len = length a in
        let median = create (len + 1) in
        for i = len - 1 downto 0 do
            let a = get a i in
            let b = get b i in
            let c = get c i in
            let to_prepend = Cost_matrix.Three_D.median a b c cm3 in
            if to_prepend <> gap then prepend median to_prepend
            else ();
        done;
        prepend median gap;
        median, cost

    let readjust_3d ?(first_gap = true) s1 s2 m cm cm3 p =
        if length s1 = length s2 && length m = length p && length s1 = length m then
            0, m, false
        else 
            let gap = Cost_matrix.Two_D.gap cm in
            if is_empty s1 gap && not (is_empty s2 gap) then
                0, s2, 0 <> compare m s2
            else if is_empty s2 gap && not (is_empty s1 gap) then 
                0, s1, 0 <> compare m s1
            else if is_empty p gap then 
                0, p, 0 <> compare m p
            else
                match first_gap with 
                | true -> 
                      let res, x = align_3_powell_inter s1 s2 p cm cm3 in
                      x, res, 0 <> compare m res
                | false ->
                      let s1 = prepend_char s1 gap in 
                      let s2 = prepend_char s2 gap in 
                      let p = prepend_char p gap in 
                      let res, x = align_3_powell_inter s1 s2 p cm cm3 in
                      let res = del_first_char res in 
                      x, res, 0 <> compare m res

end

module NewkkAlign = struct 
    
    type ukkm
    external create_ukkm : unit -> ukkm = "newkkonen_CAML_create_general"
    (*
    external print_2d : m -> int -> int -> unit = "mat_CAML_print_algn_2d"
    external print_3d : m -> int -> int -> int -> unit = "mat_CAML_print_algn_3d"
    external get_dir : m -> int -> int -> int -> int -> int = "mat_CAML_get_value"
    *)
    (*c side function*)
    external init : unit -> unit = "newkkonen_CAML_initialize"
    external newkk_cost2 : s -> s -> Cost_matrix.Two_D.m -> ukkm -> int =
         "newkkonen_CAML_algn"
    external newkk_cost2_affine : s -> s -> Cost_matrix.Two_D.m -> ukkm -> int =
         "newkkonen_CAML_algn_affine"
    external newkk_backtrace : s -> s -> s -> s -> Cost_matrix.Two_D.m -> ukkm -> unit = 
        "newkkonen_CAML_backtrace_bc" "newkkonen_CAML_backtrace"
    external newkk_backtrace_affine : s -> s -> s -> s -> Cost_matrix.Two_D.m -> ukkm -> unit = 
        "newkkonen_CAML_backtrace_affine_bc" "newkkonen_CAML_backtrace_affine"
    external get_k : ukkm -> int = "newkkonen_CAML_get_k"
    (*functions from cost_matrix.2d*)
    let gap_filter_for_combcode = Cost_matrix.Two_D.gap_filter_for_combcode

    let _ = init ()
    let default_ukkm = create_ukkm ()
(*median_2 and median_2_with_gaps are totally the same as in module Align
    let median_2 s1 s2 c =
        let sz1 = length s1 
        and sz2 = length s2 in
        if (sz1 = sz2) then begin
            let sp = create (sz1 + 1) in
            c_median_2 s1 s2 c sp;
            sp
        end else 
            raise 
            (Invalid_Argument "Newkkonen.median_2: size of two seq must be the same.")

    let median_2_with_gaps s1 s2 c =
        let sz1 = length s1 
        and sz2 = length s2 in
        if (sz1 = sz2) then begin
            let sp = create sz1 in
            c_median_2_with_gaps s1 s2 c sp;
            sp
        end else 
            raise 
            (Invalid_Argument "The size of the sequences is not the same.")
*)


    let call_newkk_backtrace s1 s2 s1p s2p c m affine =
        if affine then
            newkk_backtrace_affine s1 s2 s1p s2p c m
        else
            newkk_backtrace s1 s2 s1p s2p c m 


    let get_alignment s1 s2 c m affine =
        let debug = false in
        let sz1 = length s1
        and sz2 = length s2 in
        let s1p = create (sz1 + sz2)
        and s2p = create (sz1 + sz2) in
        if debug then Printf.printf "Sequence.Newkkonen.get_alignment, len1=%d,len2=%d\n%!" sz1 sz2;
        (*call traceback function*)
        let size_compared = (sz1 <= sz2) in
        if size_compared then 
            call_newkk_backtrace s1 s2 s1p s2p c m affine
        else 
            call_newkk_backtrace s2 s1 s2p s1p c m affine;
        if debug then begin 
            Printf.printf " seq1:%!"; printseqcode s1;
            Printf.printf " seq2:%!"; printseqcode s2;
            Printf.printf " alied seq1:%!"; printseqcode s1p;
            Printf.printf " alied seq2:%!"; printseqcode s2p;
        end;
        s1p, s2p

    let align_2 ?(first_gap=true) s1 s2 c m =
        let debug = false in
        let ls1 = length s1 and ls2 = length s2 in
        if debug then Printf.printf "Sequence.Newkkonen.align_2,len1=%d,len2=%d\n%!" ls1 ls2;
        let s1,s2,ls1,ls2,exchange = 
            if ls1<=ls2 then s1,s2,ls1,ls2,false
            else s2,s1,ls2,ls1,true
        in
        let cmp s1 s2 =
            match Cost_matrix.Two_D.affine c with
            | Cost_matrix.Affine _ ->
                    (*printseqcode s1; printseqcode s2;*)
                    let tc = newkk_cost2_affine s1 s2 c m in   
                    (*Printf.printf " editing cost = %d,call traceback\n%!"
                    * tc;*)
                    let s1p, s2p = get_alignment s1 s2 c m true in
                    if exchange then s2p,s1p,tc
                    else
                    s1p, s2p, tc
            | _ ->
                    (*printseqcode s1; printseqcode s2;*)
                    let tc = newkk_cost2 s1 s2 c m in   
                    (*Printf.printf " editing cost = %d,call traceback\n%!"
                    * tc;*)
                    let s1p, s2p = get_alignment s1 s2 c m false in
                    if exchange then s2p,s1p,tc
                    else
                    s1p, s2p, tc   
        in 
        match first_gap with
        | true -> cmp s1 s2 
        | false ->           
            let gap = Cost_matrix.Two_D.gap c in 
              let s1 = prepend_char s1 gap in 
              let s2 = prepend_char s2 gap in 
              let s1p, s2p, tc = cmp s1 s2 in
              let s1p = del_first_char s1p in 
              let s2p = del_first_char s2p in 
              if debug then Printf.printf "end of Newkkonen.alied len1=%d, len2=%d\n%!" (length s1p) (length
              s2p);
              s1p, s2p, tc

    let cost_2 ?deltaw s1 s2 m1 m2 =
        let debug = false in
        if debug then Printf.printf "Newkkonen.cost_2\n%!";
        let ls1 = length s1 and ls2 = length s2 in
        let s1,s2,ls1,ls2 = 
            if ls1<=ls2 then s1,s2,ls1,ls2
            else s2,s1,ls2,ls1
        in
        let s1,ls1 = 
          if ls1=0 then (prepend_gap s1 m1),1
          else s1,ls1
        in
        let s2,ls2 =
            if ls2=0 then (prepend_gap s2 m1),1 
            else s2,ls2
        in
        (*assert (ls1 <> 0);  assert (ls2 <> 0);*)
        if ls1 <= ls2 then
            newkk_cost2 s1 s2 m1 m2 
        else 
            newkk_cost2 s2 s1 m1 m2 

    let full_median_2 a b cm m = 
        match Cost_matrix.Two_D.affine cm with
        | Cost_matrix.Affine _ (*->
                let a, b, _ = align_2 a b cm m in
                median_2 a b cm*)
        | _ ->
                let a, b, _ = align_2 a b cm m in
                median_2 a b cm


(*closest s1 s2 cm ukkm*)
    let closest s1 s2 cm m =
        let uselevel = check_level cm in
        if is_empty s2 (Cost_matrix.Two_D.gap cm) then
            s2, 0
        else
        let (s, _) as res = 
        if 0 = Cost_matrix.Two_D.combine cm then 
            (* We always have just one sequence per type s *)
            let s1', s2', cst = align_2 s1 s2 cm m in
            let get_closest v i =
                let v' = get s1' i in
                let all = Cost_matrix.Two_D.get_all_elements cm in
                if v = all then 
                    if v' = all then 1 (* We pick one, any will be fine *)
                    else v'
                else v
            in
            remove_gaps2 (mapi get_closest s2') cm, cst
        else
            let s1', s2', comp =
                if 0 = compare s1 s2 then
                    (* We remove all the gaps if we are using combinations *)
                    if 0 = Cost_matrix.Two_D.combine cm then 
                        s1, s2, true
                    else
                        if uselevel then 
                            let gap_filter x =
                                gap_filter_for_combcode x cm
                            in
                            let news1 = map gap_filter s1 and
                            news2 = map gap_filter s2 in
                            news1, news2, true
                        else
                            let mask = lnot (Cost_matrix.Two_D.gap cm) in
                            mapi (fun x pos -> 
                                if pos > 0 then x land mask else x) s1, 
                            mapi (fun x pos -> 
                                if pos > 0 then x land mask else x) s2, 
                            true
                else
                    let s1', s2', _ = align_2 s1 s2 cm m in
                    s1', s2', false
            in
            let s2' =
                let s2' = 
                    let get_closest v i =
                        let v' = get s1' i in
                        Cost_matrix.Two_D.get_closest cm v' v 
                    in
                    mapi get_closest s2' 
                in
                remove_gaps2 s2' cm
            in
            (* We must recalculate the distance between sequences because the
            * set ditance calculation is an upper bound in the affine gap cost
            * model *)
            if comp then 
                s2', 0
            else
                s2', cost_2 s1 s2' cm m
        in
        res

end


let select_one_generic get_one_item s cm =
    if 0 = Cost_matrix.Two_D.combine cm then s
    else 
        let asz = Cost_matrix.Two_D.alphabet_size cm in
        map (get_one_item asz) s

let select_one s cm =
    let sort_list = List.sort (fun a b -> a - b) in
    let get_one_item asz b =
        match sort_list (BitSet.Int.list_of_packed_max b asz) with
        | h :: _ -> h
        | [] -> failwith "Nothing?"
    in
    select_one_generic get_one_item s cm

let select_one_randomized s cm =
    let get_one_item asz b =
        let lst = BitSet.Int.list_of_packed_max b asz in
        let len = List.length lst in
        List.nth lst (Random.int len)
    in
    select_one_generic get_one_item s cm

let readjust a b m cm parent use_ukk =
    let matr = Matrix.default in
    let algn s1 s2 =
        match Cost_matrix.Two_D.affine cm with
        | Cost_matrix.Affine _ ->  
                if use_ukk then 
                    let s1', s2', c =
                    NewkkAlign.align_2 s1 s2 cm NewkkAlign.default_ukkm
                    in
                    let median = median_2 s1' s2' cm in
                    c,median
                else
                    let m, _, _, cost, _ = Align.align_affine_3 s1 s2 cm in
                    cost, m
        | _ ->
                let s1', s2', c = 
                    if use_ukk then 
                    NewkkAlign.align_2 s1 s2 cm NewkkAlign.default_ukkm
                    else
                    Align.align_2 ~first_gap:true s1 s2 cm Matrix.default
                in
                let median = median_2 s1' s2' cm in
                c, median
    in
    let make_center a b c =
        let cab, ab = algn a b
        and cbc, bc = algn b c
        and cac, ac = algn a c in
        let cabc, abc = algn ab c
        and cbca, bca = algn bc a
        and cacb, acb = algn ac b in
        let cabc = cabc + cab
        and cbca = cbca + cbc
        and cacb = cacb + cac in
        if cabc <= cbca then
            if cabc <= cacb then false, cabc, Align.closest c ab cm matr, cabc
            else true, cacb, Align.closest b ac cm matr, cabc
        else if cbca < cacb then
            true, cbca, Align.closest a bc cm matr, cabc 
        else true, cacb, Align.closest b ac cm matr, cabc 
    in
    let has_to_print, c, (s, _), previous = make_center a b parent in
    c, s, has_to_print



module CamlAlign = struct

    type e = Align | Delete | Insert
    type ep = Alignp of int | Deletep of int | Insertp of int

    type s = int array

(** Give two code arrays and the cost matrix of translations between
 *  codes. Align two code arrays with minimum cost using deletion, insertion and
 *  substitution *)
    let create_pair_align code1_arr code2_arr cost_mat gap_code = 
        let infinity = 100000000 in    
        let len1 = Array.length code1_arr + 2 in 
        let len2 = Array.length code2_arr + 2 in 

        let row_code_arr = Array.init len1 
            (fun index -> 
                 if (index = 0) || (index = len1 - 1) then gap_code
                 else code1_arr.(index - 1) )
        in

        let col_code_arr = Array.init len2 
            (fun index -> 
                 if (index = 0) || (index = len2 - 1) then gap_code
                 else code2_arr.(index - 1) )
        in
        
        let tc_mat = Array.make_matrix len1 len2 infinity in 
        tc_mat.(0).(0) <- 0;

        let cost r c = 
        if (r < 0) || (c < 0) then infinity
        else tc_mat.(r).(c)
        in

        for r = 0 to len1 - 1 do
            for c = 0 to len2 - 1 do
                let row_code = row_code_arr.(r) in
                let col_code = col_code_arr.(c) in            
            
                let cost1 = cost (r - 1) (c -  1) + 
                    cost_mat.(row_code).(col_code) in 
                
                let cost2 = cost (r - 1) c + cost_mat.(row_code).(gap_code) in 
                let cost3 = cost r (c - 1) + cost_mat.(gap_code).(col_code) in 
                if tc_mat.(r).(c) > min (min cost1 cost2) cost3 then
                    tc_mat.(r).(c) <- min (min cost1 cost2) cost3
            done
        done;
        
        let rec trace ali r c  = 
            if (r = 0) && (c = 0) then ali
            else begin            
                let row_code = row_code_arr.(r) in
                let col_code = col_code_arr.(c) in            
                
                let cost1 = cost (r - 1) (c -  1) + 
                    cost_mat.(row_code).(col_code) in 
                
                let cost2 = cost (r - 1) c + cost_mat.(row_code).(gap_code) in 
                
                
                if cost1 = tc_mat.(r).(c) then 
                    trace ((row_code, col_code)::ali) (r-1) (c-1)
                else
                    if cost2 = tc_mat.(r).(c) then 
                        trace ((row_code, gap_code)::ali) (r-1) c
                    else trace ((gap_code, col_code)::ali) r (c - 1)
            end
        in
        
        match tc_mat.(len1 - 1).(len2 - 1) = infinity with
        | true -> [||], [||], infinity
        | false -> 
              let ali = trace [] (len1 - 1) (len2 - 1) in 
              let ali = List.filter 
                  (fun (c1, c2) -> (c1 != gap_code) || (c2 != gap_code)) ali in 
              
              let alied_seq1, alied_seq2 = List.split ali in 
              let ali_cost = tc_mat.(len1 - 1).(len2 - 1) in 
              Array.of_list alied_seq1, Array.of_list alied_seq2, ali_cost
                  
                  
            
            
 (** Give two code arrays and the cost matrix of translations between
  *  codes. Align two code arrays with minimum cost using deletion, insertion and
  *   substitution *)
    let create_triple_align code1_arr code2_arr code3_arr cost_mat gap = 
(* 
    printIntArr code1_arr;
    printIntArr code2_arr;
    printIntArr code3_arr;
    print_newline ();
    printIntMat cost_mat;
    print_newline ();
*)
        let infinity = 100000000 in    
        let len1 = Array.length code1_arr + 2 in 
        let len2 = Array.length code2_arr + 2 in 
        let len3 = Array.length code3_arr + 2 in 

        let code1_arr = Array.init len1 
            (fun index -> 
                 if (index = 0) || (index = len1 - 1) then gap
                 else code1_arr.(index - 1) )
        in
        
        let code2_arr = Array.init len2 
            (fun index -> 
                 if (index = 0) || (index = len2 - 1) then gap
                 else code2_arr.(index - 1) )
        in
        
        let code3_arr = Array.init len3 
            (fun index -> 
                 if (index = 0) || (index = len3 - 1) then gap
                 else code3_arr.(index - 1) )
        in
        
        let tc_cub = Array.init len1
            (fun _ -> Array.init len2 (fun _ -> Array.make len3 infinity))
        in 
        let cmp_cost3 code1 code2 code3 = cost_mat.(code1).(code2) +
            cost_mat.(code1).(code3) + cost_mat.(code2).(code3) 
        in 
        
        let find_min ls = 
            List.fold_left (fun v min_v -> 
                                if v < min_v then v
                                else min_v) 
                (List.hd ls) ls
        in 

        tc_cub.(0).(0).(0) <- 0;
        for i1 = 1 to len1 - 1 do
            for i2 = 1 to len2 - 1 do
                for i3 = 1 to len3 - 1 do
                    let c1 = code1_arr.(i1) in 
                    let c2= code2_arr.(i2)in 
                    let c3= code3_arr.(i3) in
                    
                    
                    let cost1 = tc_cub.(i1 - 1).(i2 - 1).(i3 - 1) + cmp_cost3 c1 c2 c3 in
                    let cost2 = tc_cub.(i1 - 1).(i2 - 1).(i3) + cmp_cost3 c1 c2 gap in                
                    let cost3 = tc_cub.(i1 - 1).(i2).(i3 - 1) + cmp_cost3 c1 gap c3 in                
                    let cost4 = tc_cub.(i1 - 1).(i2).(i3) + cmp_cost3 c1 gap gap in                
                    let cost5 = tc_cub.(i1).(i2 - 1).(i3 - 1) + cmp_cost3 gap c2 c3 in
                    let cost6 = tc_cub.(i1).(i2 - 1).(i3) + cmp_cost3 gap c2 gap in
                    let cost7 = tc_cub.(i1).(i2).(i3 - 1) + cmp_cost3 gap gap c3 in
                    let cost_ls = [cost1; cost2; cost3; cost4; cost5; cost6; cost7] in 
                    tc_cub.(i1).(i2).(i3) <- find_min cost_ls; 
                    
                done;
            done;
        done;
        

        let rec trace ali i1 i2 i3 = 
            if (i1 + i2 + i3 = 0) then ali
            else begin
                let c1 = code1_arr.(i1) in 
                let c2= code2_arr.(i2) in 
                let c3= code3_arr.(i3) in
                let min_cost = tc_cub.(i1).(i2).(i3) in 
                
                let cost1 = tc_cub.(i1 - 1).(i2 - 1).(i3 - 1) + cmp_cost3 c1 c2 c3 in
                let cost2 = tc_cub.(i1 - 1).(i2 - 1).(i3) + cmp_cost3 c1 c2 gap in                
                let cost3 = tc_cub.(i1 - 1).(i2).(i3 - 1) + cmp_cost3 c1 gap c3 in                
                let cost4 = tc_cub.(i1 - 1).(i2).(i3) + cmp_cost3 c1 gap gap in                
                let cost5 = tc_cub.(i1).(i2 - 1).(i3 - 1) + cmp_cost3 gap c2 c3 in
                let cost6 = tc_cub.(i1).(i2 - 1).(i3) + cmp_cost3 gap c2 gap in

                
                if min_cost = cost1 then trace ((c1, c2, c3)::ali) (i1 - 1) (i2 - 1) (i3 - 1)
                else if min_cost = cost2 then trace ((c1, c2, gap)::ali) (i1 - 1) (i2 - 1) i3
                else if min_cost = cost3 then trace ((c1, gap, c3)::ali) (i1 - 1) i2 (i3 - 1)        
                else if min_cost = cost4 then trace ((c1, gap, gap)::ali) (i1 - 1) i2 i3        
                else if min_cost = cost5 then trace ((gap, c2, c3)::ali) i1 (i2 - 1) (i3 - 1)        
                else if min_cost = cost6 then trace ((gap, c2, gap)::ali) i1 (i2 - 1) i3        
                else trace ((gap, gap, c3)::ali) i1 i2 (i3 - 1)        
            end             
        in 
        
        if tc_cub.(len1 - 1).(len2 - 1).(len3 - 1) >= infinity then 
            [||],[||], [||], infinity
        else begin
            let ali = trace [] (len1 - 1) (len2 - 1) (len3 - 1) in
            let ali = List.filter (fun (c1, c2, c3) -> c1 + c2 + c3 != 3 * gap) ali in 
            let alied_code1_ls, alied_code2_ls, alied_code3_ls = 
                List.fold_right (fun (c1, c2, c3) (l1, l2, l3) -> c1::l1, c2::l2,
                                     c3::l3) ali ([], [], [])
            in 
            
            let cost = List.fold_left (fun ali_cost (c1, c2, c3) -> ali_cost + cmp_cost3 c1 c2 c3) 0 ali
            in 
(*            
            (if cost != tc_cub.(len1 - 1).(len2 - 1).(len3 -1) || cost >= infinity then 
                 failwith "You are working with the triple alignment cost"
             else begin 
                 fprintf stdout "Triple_alignment_cost: %i\n" cost;
                 printIntArr (Array.of_list alied_code1_ls);
                 printIntArr (Array.of_list alied_code2_ls);
                 printIntArr (Array.of_list alied_code3_ls);
                 print_newline (); 
             end  
            );  
            
*)            
            Array.of_list alied_code1_ls, 
            Array.of_list alied_code2_ls, 
            Array.of_list alied_code2_ls, 
            cost        
        end 
            

 
(*
    let align a b =
        let lena = Array.length a
        and lenb = Array.length b in
        let cost_matrix = Array.make_matrix (lenb + 1) (lena + 1) 0 
        and even_matrix = Array.make_matrix (lenb + 1) (lena + 1) [] in
        (* Initialize the first column and row *)
        let rec initialization prev it assign max =
            if it <= max then begin
                let res = prev + 1 in
                assign it res;
                initialization res (it + 1) assign max;
            end else ()
        in
        let for_row it v = 
            cost_matrix.(0).(it) <- v;
            even_matrix.(0).(it) <- [Insert];
        in
        let for_column it v = 
            cost_matrix.(it).(0) <- v;
            even_matrix.(it).(0) <- [Delete];
        in
        (* To fill the cost cost_matrix *)
        let match_cost x y =
            if a.(x - 1) = b.(y - 1) then 0
            else 1
        in
        let fill_item x y =
            let align_cost = cost_matrix.(y - 1).(x - 1) + (match_cost x y) in
            let insert_cost = cost_matrix.(y).(x - 1) + 1
            and delete_cost = cost_matrix.(y - 1).(x) + 1 in
            let cost, event = 
                if align_cost < insert_cost then begin
                    if align_cost < delete_cost then 
                        align_cost, [Align]
                    else if delete_cost < align_cost then 
                        delete_cost, [Delete]
                    else delete_cost, [Delete; Align]
                end else if insert_cost < delete_cost then 
                    insert_cost, [Insert]
                else if delete_cost < insert_cost then 
                    delete_cost, [Delete]
                else if 
                    (insert_cost = delete_cost) 
                    && (align_cost = insert_cost)
                then
                    delete_cost, [Align; Insert; Delete]
                else if align_cost = insert_cost then
                    align_cost, [Insert; Align]
                else insert_cost, [Insert; Delete]
            in
            cost_matrix.(y).(x) <- cost;
            even_matrix.(y).(x) <- event;
        in
        let fill_matrix max_x max_y =
            for i = 1 to max_y do
                for j = 1 to max_x do
                    fill_item j i;
                done;
            done;
            ()
        in
        initialization 0 1 for_row lena;
        initialization 0 1 for_column lenb;
        fill_matrix lena lenb;
        cost_matrix, even_matrix

    let backtrace (mtx : e list array array) (a : s) (b : s) =
        let lena = Array.length a 
        and lenb = Array.length b in
        let rec builder acc x y = 
            match x, y with
            | 0, 0 -> [acc]
            | x, y ->
                    let local_builder myacc it = 
                        myacc @ (match it with
                        | Align  -> builder (a.(x - 1) :: acc) (x - 1) (y - 1)
                        | Delete -> builder (b.(y - 1) :: acc) x (y - 1)
                        | Insert -> builder (a.(x - 1) :: acc) (x - 1) y)
                    in
                    List.fold_left local_builder [] mtx.(y).(x)
        in
        List.map Array.of_list (builder [] lena lenb)
*)







    exception Illegal_argument

    (* The parametric complexity of a particular indel, sub and match model with
    * parameters [n], [m] and [alpha], where [n] is the length of one sequence,
    * [m] the length of a second sequence and [alpha] is the maximum number of
    * matches between the sequences. *)
    let pc n m alpha = 
        let fm = float_of_int m in
        (* A staged factorial calculation *)
        let rec do_division a b res =  (* a numerator, b denominator *)
            match a, b with
            | ha :: ta, hb :: tb -> 
                    do_division ta tb ((ha /. hb) *. res)
            | ha :: ta, [] ->
                    do_division ta [] (ha *. res)
            | [], hb :: tb ->
                    do_division [] tb (res /. hb)
            | [], [] -> res
        in
        let rec list_of_int acc it =
            if it > 0.0 then list_of_int (it :: acc) (it -. 1.0)
            else acc
        in
        let element_calculation m h_s =
            match h_s with
            | h1 :: h2 :: h' ->
                    let sum = List.fold_left (+.) 0.0 h' in
                    let total_sum = sum +. h1 +. h2 in
                    (* First calculate the numerator *)
                    let num = list_of_int [] m in
                    let num = list_of_int num sum in
                    let den = list_of_int [] (h1 +. h2) in
                    let den = list_of_int den (m -. (h1 +. h2)) in
                    let all_res = 
                        List.map (fun x -> (x /. total_sum) ** x) h_s 
                    in
                    let den = List.fold_left list_of_int den h' in
                    let num = num @ all_res in
                    do_division num den 1.0
            | _ -> raise Illegal_argument
        in
        let summation () =
            let result = ref 0.0 in
            for i = 0 to n do
                for j = 0 to m do
                    for k = 0 to n do
                        for l = 0 to alpha do
                            if (2 * (k + l)) + i + j = n + m && (l + k) <= n 
                            then begin
                                if not ndebug then
                                    print_endline (string_of_int l ^ "\t" ^
                                    string_of_int k ^ "\t" ^ string_of_int i ^
                                    "\t" ^ string_of_int j);
                                let i = float_of_int i
                                and j = float_of_int j 
                                and k = float_of_int k
                                and l = float_of_int l in
                                let res = 
                                    element_calculation fm [l; k; i; j] 
                                in
                                result := !result +. res;
                            end else ()
                        done;
                    done;
                done;
            done;
            !result
        in
        (log (summation ()))

        let costm = ref (Array.make_matrix 1 1 0.0)
        let dirm = ref (Array.make_matrix 1 1 (Alignp 0))

        let create_of_tcm tcm =
            let gap = Cost_matrix.Two_D.gap tcm in
            (fun cm x y s1 s2 ->
                cm.(x - 1).(y - 1) +. 
                (float_of_int (Cost_matrix.Two_D.cost (get s1 x) (get s2 y)
                tcm)),
                Alignp 1),
            (fun cm x y s ->
                cm.(x - 1).(y) +.
                (float_of_int (Cost_matrix.Two_D.cost gap (get s x) tcm)),
                Insertp 1),
            (fun cm x y s ->
                cm.(x).(y - 1) +.
                (float_of_int (Cost_matrix.Two_D.cost (get s y) gap tcm)),
                Deletep 1),
            (fun cm x s ->
                cm.(x - 1).(0) +.
                (float_of_int (Cost_matrix.Two_D.cost gap (get s x) tcm)),
                Insertp 1),
            (fun cm y s ->
                cm.(0).(y - 1) +.
                (float_of_int (Cost_matrix.Two_D.cost (get s y) gap tcm)),
                Deletep 1)

        let create_of_non_linnear_indel f tcm =
            let a = 
                (fun cm x y s1 s2 ->
                    cm.(x - 1).(y - 1) +. 
                    (float_of_int (Cost_matrix.Two_D.cost (get s1 x) (get s2 y)
                    tcm)),
                    Alignp 1)
            and b =
                (fun cm x y s ->
                    let min = ref max_float
                    and len = ref 0 in
                    for i = 0 to x - 1 do
                        let cost = cm.(x - i).(y) +. (f i x s) in
                        if cost < !min then begin 
                            min := cost;
                            len := x - i;
                        end;
                    done;
                    !min, Insertp !len)
            and c =
                (fun cm x y s ->
                    let min = ref max_float
                    and len = ref 0 in
                    for i = 0 to y - 1 do
                        let cost = cm.(x).(y - i) +. (f i y s) in
                        if cost < !min then begin 
                            min := cost;
                            len := y - i;
                        end;
                    done;
                    !min, Insertp !len) 
            in
            let d = fun cm x s -> b cm x 0 s
            and e = fun cm y s -> c cm 0 y s in
            a, b, c, d, e

        let cost_2 cf insf delf rowf colf s1 s2 = 
            let s1l = length s1
            and s2l = length s2 in
            let max = max s1l s2l in
            if max > Array.length !costm.(0) then begin
                costm := Array.make_matrix max max 0.0;
                dirm := Array.make_matrix max max (Alignp 0);
            end;
            let costm = !costm
            and dirm = !dirm in
            for i = 1 to s1l - 1 do
                let a, b = rowf costm i s1 in
                costm.(i).(0) <- a;
                dirm.(i).(0) <- b;
            done;
            for i = 1 to s2l - 1 do
                let a, b = colf costm i s2 in
                costm.(0).(i) <- a;
                dirm.(0).(i) <- b;
            done;
            for i = 1 to s1l - 1 do 
                for j = 1 to s2l - 1 do
                    let insa, dirins = insf costm i j s1 in
                    let dela, dirdel = delf costm i j s2 in
                    let suba, dirsub = cf costm i j s1 s2 in
                    if suba <= insa then
                        if suba <= dela then begin
                            costm.(i).(j) <- suba;
                            dirm.(i).(j) <- dirsub;
                        end else begin
                            costm.(i).(j) <- dela;
                            dirm.(i).(j) <- dirdel;
                        end
                    else if insa <= dela then begin
                        costm.(i).(j) <- insa;
                        dirm.(i).(j) <- dirins;
                    end else begin
                        costm.(i).(j) <- dela;
                        dirm.(i).(j) <- dirdel;
                    end;
                done;
            done;
            costm.(s1l - 1).(s2l - 1) 

        let cost_2_of_tcm tcm =
            let a, b, c, d, e = create_of_tcm tcm in
            cost_2 a b c d e

        let non_linnear f tcm =
            let a, b, c, d, e = create_of_non_linnear_indel f tcm in
            cost_2 a b c d e

        let backtrace_ext s1 s2 gap = 
            let s1l = length s1 
            and s2l = length s2 in
            let i = ref (s1l - 1)
            and j = ref (s2l - 1) 
            and seq1 = create (s1l + s2l + 1) 
            and seq2 = create (s1l + s2l + 1)
            and dirm = !dirm in
            while (!i <> 0) && (!j <> 0) do
                match dirm.(!i).(!j) with
                | Alignp len -> 
                        for k = 1 to len do
                            prepend seq1 (get s1 !i);
                            prepend seq2 (get s2 !j);
                            decr i;
                            decr j;
                        done;
                | Deletep len ->
                        for k = 1 to len do
                            prepend seq2 (get s2 !j);
                            prepend seq1 gap;
                            decr j;
                        done;
                | Insertp len ->
                        for k = 1 to len do
                            prepend seq1 (get s1 !i);
                            prepend seq2 gap;
                            decr i;
                        done;
            done;
            if gap <> get seq1 0 then prepend seq1 gap;
            if gap <> get seq2 0 then prepend seq2 gap;
            seq1, seq2

end

(*external fix_randomly : s -> unit = "seq_CAML_fix_randomly"*)

let proportion a b = (float_of_int a) /. (float_of_int b)

let gap_saturation seq alph = 
    let gap = Alphabet.get_gap alph in
    let len = length seq
    and gaps = 
        fold
        (fun acc base -> if 0 <> base land gap then acc + 1 else acc)
        0
        seq
    in
    proportion gaps len

(* Calculate the fraction of positions in the sequence that show
* polymorphism. *)
let poly_saturation sequence n =
    let len = length sequence 
    and poly =
        fold 
        (fun acc base -> if n = (BitSet.Int.count_bits base) then acc + 1 else acc)
        0
        sequence
    in
    proportion poly len

let split positions s alph = 
    let gap = Alphabet.get_gap alph in
    let len = (length s) - 1 in
    let positions = (0, 0) :: positions in
    let do_one_pair a b acc =
        let first = a
        and last = b in
        let total = 1 + (last - first) in
        let seq = create (total + 1)  in
        for i = (last - 1) downto first do
            prepend seq (get s i);
        done;
        if first <> 0 then prepend seq gap;
        seq :: acc
    in
    let rec splitter acc = function
        | (a, _) :: (((c, b) :: _) as t) -> 
                assert (
                    if a <= b then true
                    else 
                        let _ = print_endline 
                        ("Trying to split " ^ string_of_int a ^ 
                        " and " ^ string_of_int b) in
                        false
                );
                splitter (do_one_pair a c acc) t
        | (a, _) :: [] ->
                (* We add one at the end because we remove one in do_one_pair *)
                List.rev (do_one_pair a (len + 1) acc)
        | [] -> []
    in
    let res = splitter [] positions in
    assert (
        let res = List.map (fun x -> to_string x alph) res in
        let res = String.concat "" res in
        let initial = to_string s alph in
        let res = Str.global_replace (Str.regexp "-") "" res in
        let assertion = ("-" ^ res) = initial in
        if not assertion then Printf.printf "Initial: %s\nFinal: %s\n%!" initial
        res;
        assertion);
    res

module Unions = struct
    (* 
    * ATTENTION, WARNING: The following structure is used AS IS in the C side,
    * don't change it unless you change also the C side (union.c and union.h). 
    * *)
IFDEF USE_LONG_SEQUENCES THEN
        type off_type = 
            (int32, Bigarray.int32_elt, Bigarray.c_layout) Bigarray.Array1.t
        let zero = Int32.zero
        let to_int = Int32.to_int
        let of_int = Int32.of_int
        let constructor = Bigarray.int32
ELSE
        type off_type = 
            (int, Bigarray.int16_signed_elt, Bigarray.c_layout) Bigarray.Array1.t
        let zero = 0
        let to_int x = x
        let of_int x = x
        let constructor = Bigarray.int16_signed
END

        type u = {
            seq : s;
            offset : off_type;
            union_c1 : off_type;
            union_c2 : off_type; } 
        let create_seq = create 

        let create_union do_init s = 
            let create cap =
                    Bigarray.Array1.create 
                    constructor
                    Bigarray.c_layout 
                    cap
            in
            let init_0 res =
                let _ = 
                    let cap = Bigarray.Array1.dim res in
                    for i = cap - 1 downto 0 do
                        res.{i} <- zero;
                    done
                in
                res
            in
            let init_x len res = 
                let _ = 
                    let cap = Bigarray.Array1.dim res in
                    let offset = cap - len in
                    for i = 0 to len - 1 do
                        res.{i + offset} <- of_int i;
                    done;
                in
                res
            in
            let cap = capacity s 
            and len = length s in
            { 
                seq = s; 
                offset = 
                    if do_init then init_0 (create cap)
                    else create cap;
                union_c1 = 
                    if do_init then init_x len (create cap)
                    else create cap;
                union_c2 = 
                    if do_init then init_x len (create cap)
                    else create cap;
            }

        let leaf s =
            let s = clone s in
            assert (length s == capacity s);
            create_union true s 

        external make_union : 
            s -> s -> s -> u -> u -> u -> Cost_matrix.Two_D.m -> unit = 
                "union_CAML_make_b" "union_CAML_make"

        let union a b m ua ub cm =
            let new_seq, len = 
                let c = (length ua.seq) + (length ub.seq) + 2 in
                create c, c
            in
            let u = create_union false new_seq in
            make_union a b m ua ub u cm;
            u

        let get_seq u = u.seq

        let get_positions ua positions = 
            let len = (Bigarray.Array1.dim ua.offset) in
            let seq_len = length ua.seq in
            let rec pos_finder modifier arr x =
                if x = 0 then arr.{0}
                else if x = len then pos_finder ( - ) arr (x - 1)
                else if (x < len) && (zero <> arr.{x}) then arr.{x}
                else pos_finder modifier arr (modifier x 1)
            in
            let get_position arr poslen pos (x, y) = 
                try
                    let xmod = 
                        if pos = poslen then ( + )
                        else ( - )
                    and ymod = 
                        if pos = 0 then ( - )
                        else ( + ) 
                    in
                    let posx = len - seq_len + x 
                    and posy = len - seq_len + y in
                    to_int (pos_finder xmod arr posx), 
                    to_int (pos_finder ymod arr posy)
                with
                | err ->
                        print_endline ("I have an error with len " ^
                        string_of_int len ^ " and sequence length " ^ 
                        string_of_int seq_len ^ " with " ^ 
                        string_of_int x ^ " and " ^ string_of_int y ^ 
                        " for the sequence " ^ to_string ua.seq
                        Alphabet.nucleotides);
                        raise err;
            in
            let positions = Array.of_list positions in
            let len = (Array.length positions) - 1 in
            let a = Array.mapi (get_position ua.union_c1 len) positions 
            and b = Array.mapi (get_position ua.union_c2 len) positions in
            Array.to_list a, Array.to_list b

        let compare a b = 
            match compare a.seq b.seq with
            | 0 -> Pervasives.compare a.offset b.offset
            | x -> x

end

let of_code_arr code_arr gap = 
    let num_nus = Array.fold_left   
        (fun num_nus dna_code ->  
             if (dna_code = gap) or (dna_code = 0) then num_nus 
             else num_nus + 1  
        ) 0 code_arr 
    in  

    let seq = init (fun _ -> gap) num_nus in 

    let _ =  
        Array.fold_left   
            (fun num_nus code ->                                              
                 match (code=gap) or (code = 0) with 
                 | true -> num_nus  
                 | false ->  
                       set seq num_nus code;  
                       num_nus + 1 
            ) 0 code_arr 
    in   
    seq


external encoding : 
    (float, Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array1.t -> 
        s -> float = "seq_CAML_encoding"

let aux_complement start a s = 
    let res =
        let acc = (create (length s)) in
        for i = start to (length s) - 1 do
            match Alphabet.complement (get s i) a with
            | Some x -> prepend acc x
            | None -> failwith "I can't complement this alphabet"
        done;
        acc
    in
    res

let complement a s = 
    let gap = Alphabet.get_gap a in
    let res = aux_complement 1 a s in
    prepend res gap;
    res

let complement_chrom = aux_complement 0 

(** [is_existed_code code seq] returns true
* if the [code] is in the sequence [seq], otherwise false *)
let is_existed_code code seq = 
    fold (fun existed c -> if c = code then true
                   else existed) false seq

(** [is_existed_char code seq] returns true
* if the code of character [ch] is in the sequence [seq], otherwise false *)
let is_existed_char ch seq = 
    let code = Alphabet.match_base ch Alphabet.nucleotides in 
    is_existed_code code seq


(** [create_gap_seq gap len] create a 
* sequence of  [len] gaps *)
let create_gap_seq ?(gap=Alphabet.gap) len = 
    init (fun _ -> gap) len

let is_gap_seq ?(gap=Alphabet.gap) seq =
    let num_gap = count gap 0 seq in
    let len = length seq in
    if ( num_gap <> len ) then 0
    else 1

(** [cmp_num_all_DNA] returns number of
* codes in [seq] which do not include gap *)
let cmp_num_all_DNA seq = 
    let len = length seq in 
    let gap = Alphabet.gap in 
    let num_nu = ref 0 in
    for p = 0 to len - 1 do 
        if (get seq p) land gap != gap then num_nu := !num_nu + 1;
    done;
    !num_nu

(** [cmp_num_not_all_DNA] returns number of
* codes in [seq] which are not gaps *)
let cmp_num_not_gap seq = 
    let len = length seq in 
    let gap = Alphabet.gap in 
    let num_nu = ref 0 in
    for p = 0 to len - 1 do 
        if (get seq p) != gap then num_nu := !num_nu + 1;
    done;
    !num_nu

(** [cmp_gap_cost indel seq] returns the indel cost
* of sequence [seq] *)
let cmp_gap_cost indel seq = 
    let num_char = cmp_num_all_DNA seq in 
    match num_char with 
    | 0 -> 0
    | _ ->
          let o, e = indel in 
          o + num_char * e / 100


(** [cmp_ali_cost alied_seq1 alied_seq2 direction cost_mat]
* returns the alignment cost for alignment [alied_seq1] and [alied_seq2] *)
let cmp_ali_cost (alied_seq1 : s) (alied_seq2 : s) 
        direction (cost_mat : Cost_matrix.Two_D.m) =

    let opening_cost = 
        match Cost_matrix.Two_D.affine cost_mat with
        | Cost_matrix.Affine o -> o
        | _ -> 0
    in 

    
	let len = length alied_seq1 in  
    let gap = Cost_matrix.Two_D.gap cost_mat in 

    let rec sum_up (pre_b1, pre_b2, total_cost) p  = 
        if p >= len then total_cost 
        else begin
		    let b1 = get alied_seq1 p in 
		    let b2 =   
                match direction = `Positive with  
                | true -> get alied_seq2 p   
                | false -> get alied_seq2 (len - 1 - p)  
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
(** [get_empty_seq] returns an empty sequence *)
let get_empty_seq () = create 0


(** [subseq seq start len] returns a subsequence
* of sequence [seq] start from position [start] 
* for a length of [len]. If length [len] < 1, return 
* an empty sequence *)
let subseq seq start len = 
	match len < 1 with
	| true -> get_empty_seq ()
    | false -> sub seq start len

(** [subseq_ignore_gap ?gap seq start len] get subseq of seq, start at idx
* 'start', subseq contains number of 'len' non-gap charactor. 
* since gap does not count in 'len'. ending idx is returned as well.
 * *)
let subseq_ignore_gap ?(gap=Alphabet.gap) seq start len =
    match len < 1 with
	| true -> get_empty_seq (),start
    | false -> sub_ignore_gap ~gap:gap seq start len  


(** [align2 seq1 seq2 cost_mat] aligns 
* two sequences [seq1], [seq2] without conditioning 
* that two first characters are gaps *)
let align2 (seq1 : s) (seq2 : s) 
        (cost_mat : Cost_matrix.Two_D.m) use_ukk =

	let len1 = length seq1 in
	let len2 = length seq2 in
	

	let gap_code = Cost_matrix.Two_D.gap cost_mat in 
	let ext_seq1 = init (fun pos -> if pos = 0 then gap_code 
						else get seq1 (pos - 1)) (len1 + 1) in 
	let ext_seq2 = init (fun pos -> if pos = 0 then gap_code 
						else get seq2 (pos - 1)) (len2 + 1) in
(*this calls the ukkonen space saver version. in ml. I'm moving it to C
    let alied_arr1,alied_arr2, cost = 
        Ukk_space_save.ukkonen_align (to_array ext_seq1) (to_array ext_seq2) cost_mat in
    let ext_alied_seq1, ext_alied_seq2 = 
        of_array alied_arr1, of_array alied_arr2 in
*)

    let ext_alied_seq1, ext_alied_seq2, cost =
(*this calls the c version of ukk_space_save*)
    if use_ukk then
        NewkkAlign.align_2 ext_seq1 ext_seq2 cost_mat NewkkAlign.default_ukkm 
    else
(* this is calls ukkonen full memory version, in C, which allocation the whole matrix*)
        Align.align_2 ext_seq1 ext_seq2 cost_mat Matrix.default 
    in 		


(*just for test, we call them both and compare the cost
let ext_alied_seq3, ext_alied_seq4, cost2 =
Align.align_2 ext_seq1 ext_seq2 cost_mat Matrix.default 
in
Printf.printf "cost from align = %d\n%!" cost2;
let ext_alied_seq1, ext_alied_seq2, cost1 =
        NewkkAlign.align_2 ext_seq1 ext_seq2 cost_mat NewkkAlign.default_ukkm 
in
Printf.printf "cost from newkk = %d\n%!" cost1;

assert(cost1=cost2);
let cost = cost1 in
just for test, we call them both and compare the cost*)

	let ali_len = length ext_alied_seq1 - 1 in 
	let alied_seq1 = subseq ext_alied_seq1 1 ali_len in 
	let alied_seq2 = subseq ext_alied_seq2 1 ali_len in 

	alied_seq1, alied_seq2, cost, ali_len


(** [align3] returns the alignment of three sequences
* [seq1], [seq2], and [seq3] without requiring 
* that three first letters of three sequences to be gaps *)
let align3 (seq1 : s) (seq2 : s) (seq3 : s) 
           (cost_cube : Cost_matrix.Three_D.m)= 

	let len1 = length seq1 in
	let len2 = length seq2 in
	let len3 = length seq3 in
	
	let gap_code = Alphabet.gap in 
	let ext_seq1 = init (fun pos -> if pos = 0 then gap_code 
						else get seq1 (pos - 1)) (len1 + 1) in 
	let ext_seq2 = init (fun pos -> if pos = 0 then gap_code 
						else get seq2 (pos - 1)) (len2 + 1) in 
	let ext_seq3 = init (fun pos -> if pos = 0 then gap_code 
						else get seq3 (pos - 1)) (len3 + 1) in 

    let nuc = Alphabet.nucleotides in 
	print stdout ext_seq1 nuc; print_newline ();
	print stdout ext_seq2 nuc; print_newline ();
	print stdout ext_seq3 nuc; print_newline ();
	
	let ext_alied_seq1, ext_alied_seq2, ext_alied_seq3, cost = 
        Align.align_3 ext_seq1 ext_seq2 ext_seq3 cost_cube 
            Matrix.default in 		
	
	print stdout ext_alied_seq1 nuc; print_newline ();
	print stdout ext_alied_seq2 nuc; print_newline ();
	print stdout ext_alied_seq3 nuc; print_newline ();

	print_endline "End of POY align_3";
	let ali_len = length ext_alied_seq1 - 1 in 
	let alied_seq1 = subseq ext_alied_seq1 1 ali_len in 
	let alied_seq2 = subseq ext_alied_seq2 1 ali_len in 
	let alied_seq3 = subseq ext_alied_seq3 1 ali_len in 
	alied_seq1, alied_seq2, alied_seq3, cost, ali_len


(** [closest_alied_seq alied_parent alied_child c2]
* returns the single sequence of sequence [alied_child]
* which is closest to the aligned parent sequence [alied_parent] *)
let closest_alied_seq alied_parent alied_child c2 = 
    let len = length alied_parent in
    let single_seq = init 
        (fun p -> 
             let p_code = get alied_parent p in 
             let c_code = get alied_child p in 
             Cost_matrix.Two_D.get_closest c2 p_code c_code
        ) len 
    in 
    let cost = foldi  
        (fun cost p single_code ->
             let p_code = get alied_parent p in 
             cost + (Cost_matrix.Two_D.cost single_code p_code c2)
        ) 0 single_seq 
    in 
    single_seq, cost
    
(** [to_single_seq inseq cost_mat] resolve polymorphism in inseq with cost_mat.
* now if there is a combination code of a/t/g/c/gap, we just pick one closest to
* gap*)
let to_single_seq inseq c2 =
    let len = length inseq in
    let single_seq = init 
    (fun p -> 
             let c_code = get inseq p in 
             Cost_matrix.Two_D.get_closest c2 (Cost_matrix.Two_D.gap c2) c_code
    ) len in
    single_seq
	
	
(** [concat seq_ls] returns a concatination 
* of sequences in the sequence list [seq_ls] *)
let concat (seq_ls : s list) = 
    let total_len = List.fold_left (fun acc_len seq -> 
                        acc_len + length seq) 0 seq_ls in

    let concat_seq = init (fun index -> -1) total_len in            
	let concat_pos = ref 0 in    
    let copier seq = 
		let len = length seq in
       	for pos = 0 to len - 1 do
       		set concat_seq !concat_pos (get seq pos);
			concat_pos := !concat_pos + 1
        done;
    in
    List.iter copier seq_ls;
    concat_seq

(** [create_subalign2 seq1 seq2 cost_mat 
* start_pos1 end_pos1 start_pos2 end_pos2] returns
* the alignment of subsequence from [start_pos1] to [end_pos1] 
* in the first sequence [seq1] with subsequence from [start_pos2] to [end_pos2] 
* in the second sequence [seq2]  *)
let create_subalign2 (seq1 : s) (seq2 : s) 
        (cost_mat : Cost_matrix.Two_D.m) (start_pos1 : int) (end_pos1 : int) 
        (start_pos2 : int) (end_pos2 : int) (use_ukk : bool) = 	

	let len1 = end_pos1 - start_pos1 + 1 in
	let len2 = end_pos2 - start_pos2 + 1 in 
	let subseq1 = sub seq1 start_pos1 len1 in 
	let subseq2 = sub seq2 start_pos2 len2 in	
	let alied_subseq1, alied_subseq2, cost, ali_len = 
        align2 subseq1 subseq2 cost_mat use_ukk in 	
	alied_subseq1, alied_subseq2, cost


(** [dna_gap] returns the gap code in tha nucleotide alphabet *)
let dna_gap = Alphabet.get_gap Alphabet.nucleotides 

(** [get_num_base seq] returns number bases, not gaps, of sequence [seq] *)
let get_num_base (seq : s) = 
	fold (fun num_code code -> 
                       if code != dna_gap then num_code + 1 
                       else num_code) 0 seq
	
	
(** [delete_gap] deletes all gaps in the sequence [seq] *)
let delete_gap ?(gap_code = dna_gap) seq = 
    let new_len = fold 
        (fun len code -> 
             if code = gap_code then len 
             else (len + 1) ) 0 seq in 
	let new_seq = init (fun _ -> -1) new_len in 
	let _ = fold (fun new_pos code -> 
                               match code = gap_code with
							   | true -> new_pos
							   | false -> set new_seq new_pos code;
								     new_pos + 1) 0 seq in
	new_seq
	

(** [create_median_gap seq start_pos end_pos cost_mat] 
* returns the median sequence between sequence [seq] and gaps *)
let create_median_gap seq ?(start_pos=(-1)) ?(end_pos=(-1)) cost_mat =
    let start_pos, end_pos = 
        match start_pos with
        | -1 -> 
              let len = length seq in
              0, (len - 1)
        | _ -> start_pos, end_pos
    in 
    let gap_code = Cost_matrix.Two_D.gap cost_mat in 

    let med = init 
        (fun pos -> 
             let code1 = get seq (start_pos + pos) in
             let med_code = Cost_matrix.Two_D.median code1 gap_code cost_mat in 
             med_code) (end_pos - start_pos + 1)
    in
    med


(** [create_median_seq approx alied_seq1 alied_seq2 cost_mat]
* returns the median sequence between aligned sequence [alied_seq1] 
* and aligned sequence [alied_seq2] *)
(* Note : this only works for two non-indel alied sequence.*)
let create_median_seq ?(approx=`BothSeq) alied_seq1 alied_seq2 cost_mat =
    let len = length alied_seq1 in 
    let get_median_code pos = 
        let code1 = get alied_seq1 pos in 
        let code2 = get alied_seq2 pos in          
        match approx with 
        | `First -> code1
        | `Second -> code2
        | `BothSeq ->              
              Cost_matrix.Two_D.median code1 code2 cost_mat
    in
    let median = init (fun pos -> get_median_code pos) len in
    let cost = ref 0 in 
    for p = 0 to len - 1 do 
        let code1 = get alied_seq1 p in 
        let code2 = get alied_seq2 p in         
        (*this cost is wrong if one of input seq is a gap seq*)
        cost := !cost + (Cost_matrix.Two_D.cost code1 code2 cost_mat)
    done;
    median, !cost



let create_median_deled_seq seq cost_mat =
    let len = length seq in 
    let gap = Alphabet.gap in
    let get_median_code pos = 
        let code = get seq pos in 
        if code land gap = 0 then code
        else gap
    in

    let median = init (fun pos -> get_median_code pos) len in

    median

(** [create_median approx seq1 seq2 s1 e1 s2 e2 cost_mat]
* returns the median sequence between subsequence
* from [s1] to [e1] in the first sequence [seq1] and subsequence 
* from [s2] tp [e2] in the second sequence [seq2] *)
let create_median ?(approx=`BothSeq) seq1 seq2 
        ?(s1=(-1)) ?(e1=(-1)) ?(s2=(-1)) ?(e2=(-1)) cost_mat use_ukk = 

    let s1, e1, s2, e2 =
        match s1 with  
        | -1 ->  
              let len1 = length seq1 in 
              let len2 = length seq2 in 
              0, (len1 - 1), 0, (len2 -1)        
        | _ -> s1, e1, s2, e2
    in 

    let alied_seq1, alied_seq2, _ = 
        create_subalign2 seq1 seq2 cost_mat s1 e1 s2 e2 use_ukk
    in
    let alied_med, cost = create_median_seq ~approx:approx alied_seq1 alied_seq2 cost_mat in 
    alied_med, alied_seq1, alied_seq2, cost



let check_repeated_char seq alpha =  
    let len = length seq in  
    let rec check_char p1 p2 =  
        if p1 = len then () 
        else if p2 = len then check_char (p1 + 1) (p1 + 2) 
        else if get seq p1 = get seq p2 then begin   
            let ch = Alphabet.match_code (get seq p1) alpha in  
            print_endline ("Character " ^ ch ^ " appears twice in sequence"); 
            print stdout seq alpha;  
            print_newline ();  
            failwith "In Breakinv, characters MUST BE NOT duplicated"; 
        end else check_char p1 (p2 + 1)                       
    in  
    check_char 0 1 

    
(** [create_general_ali code1_arr code2_arr gap_code cost_mat]
* returns the general alignment between sequence [code1_arr] and 
* sequence [code2_arr] *)
let create_general_ali code1_arr code2_arr gap_code cost_mat use_ukk =
(*    print_endline "Create general alignment"; *)

    let len1 = Array.length code1_arr in
    let ext_seq1 = init 
        (fun index -> 
             match index with
             | 0 -> gap_code
             | _ -> code1_arr.(index - 1) ) (len1 + 1) 
    in 

    let len2 = Array.length code2_arr in 
    let ext_seq2 = init 
        (fun index -> 
             match index with
             | 0 -> gap_code 
             | _ -> code2_arr.(index - 1) ) (len2 + 1) 
    in 
	let ext_alied_seq1, ext_alied_seq2, cost = 
        if use_ukk then
        NewkkAlign.align_2 ext_seq1 ext_seq2 cost_mat
         NewkkAlign.default_ukkm
        else
        Align.align_2 
        ext_seq1 ext_seq2 cost_mat Matrix.default in 		

    let ali_len = length ext_alied_seq1 in 
    let alied_seq1 = to_array ext_alied_seq1 in 
    let alied_seq1 = Array.sub alied_seq1 1 (ali_len - 1) in 
    let alied_seq2 = to_array ext_alied_seq2 in 
    let alied_seq2 = Array.sub alied_seq2 1 (ali_len - 1) in 
    alied_seq1, alied_seq2, cost


(** do not use map of sequences, because it i is running up-down*)
let map f s =
    let len = length s in 
    let m_s = init (fun _ -> 0) len in
    for p = 0 to len - 1 do
        set m_s p (f (get s p));
    done; 
    m_s

(** [get_single_seq seq c2] returns a signle sequence
* of sequence [seq] *)        
let get_single_seq seq c2 = select_one seq c2

(** [cmp_locus_indel_cost s c2 locus_indel]
* returns the locus indel cost of locus [s] *)
let cmp_locus_indel_cost s c2 locus_indel =
    let locus_open, locus_ext = locus_indel in
    let gap_open = 
        match Cost_matrix.Two_D.affine c2 with
        | Cost_matrix.Affine o -> o
        | _ -> 0
    in 
    let gap = Cost_matrix.Two_D.gap c2 in 
    let seq_len = length s in


    let cmp_indel_cost p = 
        let dna = get s p in
        if (dna land gap) = gap then 0
        else Cost_matrix.Two_D.cost dna gap c2 
    in 

    let f1 = locus_open + locus_ext in 
    let f2 = locus_open + gap_open + 
             (cmp_indel_cost 0)
    in  

    let rec cmp p f1 f2 =
        if p = seq_len then min f1 f2
        else begin
            let new_f1 = (min f1 f2) + locus_ext in 
            let indel_cost = cmp_indel_cost p in 
            let new_f2 = min (f1 + gap_open + indel_cost)
                             (f2 + indel_cost)
            in 
            cmp (p + 1) new_f1 new_f2
        end 
    in
    cmp 1 f1 f2

module Clip = struct

    let count_tip_gaps gap seqa seqb =
        let rec aux pos acc = 
            let abase = get seqa pos
            and bbase = get seqb pos in
            if gap = abase || gap = bbase then
                if gap = abase && gap = bbase then 
                    aux (pos + 1) acc
                else aux (pos + 1) (acc + 1)
            else acc
        in
        aux 0 0

    let correct_distance gap c2 seqa seqb cost =
        let gap_opening = Cost_matrix.Two_D.gap_opening c2 in
        let is_gap seq pos = gap = get seq pos in
        let rec add_initial pos acc =
            if is_gap seqa pos || is_gap seqb pos then
                add_initial (pos + 1)
                (acc + (Cost_matrix.Two_D.cost (get seqa pos) 
                (get seqb pos) c2))
            else acc
        in
        let subst = add_initial 0 0 in
        if subst = 0 then cost
        else cost - (subst + gap_opening)

    let corrected_distance c2 missing_distance a b use_ukk =
        let gap = Cost_matrix.Two_D.gap c2 in
        if is_empty a gap || is_empty b gap then 
            missing_distance
        else
            let seqa, seqb, cost = 
                if use_ukk then 
                    NewkkAlign.align_2 a b c2 NewkkAlign.default_ukkm
                else 
                    Align.align_2 a b c2 Matrix.default
            in
            correct_distance gap c2 seqa seqb cost

    type old_s = s
    type s = [ `DO of old_s | `First of old_s | `Last of old_s ]

    let extract_s s = 
        match s with
        | `DO x | `First x | `Last x -> x

    let fraction = 1.00

    let clip_n_fix s1 s2 = 
        (* We return the sequences with their corrections, and functions 
        * to correct the median and the pairwise alignments *)
        let identity x = x in
        let aux s1 s2 s1len s2len =
            assert (s1len >= s2len);
            if s1len > (s2len *. fraction) then
                (* We eliminate enough to make the length difference 
                * only 10% of the overall length *)
                let maxlen = truncate (s2len *. fraction) in
                let dif = (truncate s1len) - maxlen in
                let cnt = ref (truncate s1len) in
                let s1' = init (fun x ->
                    decr cnt;
                    if x = 0 then get s1 0
                    else get s1 !cnt) maxlen 
                in
                let add_bottom x =
                    let len = length x in
                    let newx = create (len + dif) in
                    for i = (len - 1) downto 1 do
                        prepend newx (get x i);
                    done;
                    newx
                in
                let fixs1 x =
                    let newx = add_bottom x in
                    for i = dif downto 0 do
                        prepend newx (get s1 i);
                    done;
                    newx
                in
                let fixs2 x =
                    let newx = add_bottom x in
                    let gap = get s1 0 in
                    for i = dif downto 0 do
                        prepend newx gap;
                    done;
                    newx
                in
                s1', s2, fixs1, fixs2, dif
            else 
                s1, s2, identity, identity, 0
        in
        let s1len = float_of_int (length s1)
        and s2len = float_of_int (length s2) in
        if s1len >= s2len then aux s1 s2 s1len s2len
        else 
            let a, b, c, d, dif = aux s2 s1 s2len s1len in
            b, a, d, c, dif

    module Align = struct
        let align_2 ?(first_gap=true) (s1 : s) (s2 : s) c m =
            let gap = Cost_matrix.Two_D.gap c in
            let algn_and_fix s1 s2 = 
                let s1, s2, f1, f2, dif = clip_n_fix s1 s2 in
                let s1, s2, cost = Align.align_2 ~first_gap s1 s2 c m in
                let cost = correct_distance gap c s1 s2 cost in
                let tip = count_tip_gaps gap s1 s2 in
                (f1 s1), (f2 s2), cost, dif + tip, s1, s1
            in
            match s1, s2 with
            | `DO s1, `DO s2 ->
                    let s1, s2, cost = Align.align_2 ~first_gap s1 s2 c m in
                    `DO s1, `DO s2, cost, 0, `DO s1, `DO s2
            | `Last s1, `DO s2
            | `DO s1, `Last s2
            | `Last s1, `Last s2 ->
                    let s1 = safe_reverse s1
                    and s2 = safe_reverse s2 in
                    let s1, s2, cost, dif, s1o, s2o = 
                        algn_and_fix s1 s2 in
                    `Last (safe_reverse s1), `Last (safe_reverse s2), 
                    cost, dif, `Last (safe_reverse s1o),
                    `Last (safe_reverse s2o)
            | `First s1, `DO s2
            | `DO s1, `First s2
            | `First s1, `First s2 ->
                    let s1, s2, cost, dif, s1o, s2o = 
                        algn_and_fix s1 s2 in
                    `First s1, `First s2, cost, dif, 
                    `First s1o, `First s2o
            | `First s1, `Last s2 -> assert false
            | `Last s1, `First s2 -> assert false
    end

    let cmp_ali_cost a b c d =
        cmp_ali_cost (extract_s a) (extract_s b) c d
    let cmp_gap_cost a s = cmp_gap_cost a (extract_s s)
    let complement_chrom alph s =
        match s with
        | `DO s -> `DO (complement_chrom alph s)
        | `First s -> `First (complement_chrom alph s)
        | `Last s -> `Last (complement_chrom alph s)

    let concat s = 
        `DO (concat (List.map extract_s s))

    let create kind len = 
        let s = create len in
        match kind with
        | `DO -> `DO s
        | `Last -> `Last s
        | `First -> `First s

    let delete_gap ?(gap_code = dna_gap) s = 
        match s with
        | `DO s -> `DO (delete_gap ~gap_code s)
        | `Last s -> `Last (delete_gap ~gap_code s)
        | `First s -> `First (delete_gap ~gap_code s)
    let foldi f acc s = foldi f acc (extract_s s)
    let get s i = get (extract_s s) i
    let get_empty_seq kind = 
        let empty = get_empty_seq () in
        match kind with
        | `DO -> `DO empty
        | `First -> `First empty
        | `Last -> `Last empty
    let init kind f l = 
        let s = init f l in
        match kind with
        | `DO -> `DO s
        | `First -> `First s
        | `Last -> `Last s
    let is_empty s gap = is_empty (extract_s s) gap
    let length s = length (extract_s s)
    let prepend s i = prepend (extract_s s) i
    let sub s a b = 
        match s with
        | `DO s -> `DO (sub s a b)
        | `First s -> `First (sub s a b)
        | `Last s -> `Last (sub s a b)

    let to_array s = to_array (extract_s s)
    let to_formater s alph = to_formater (extract_s s) alph
    let to_string s alph = to_string (extract_s s) alph
    let safe_reverse s =
        match s with
        | `DO s -> `DO (safe_reverse s)
        | `First s -> `First (safe_reverse s)
        | `Last s -> `Last (safe_reverse s)
end
