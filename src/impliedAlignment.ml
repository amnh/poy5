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

let () = SadmanOutput.register "ImpliedAlignment" "$Revision: 2871 $"

exception NotASequence of int

let debug = false

module Codes = All_sets.IntegerMap
module Handles = All_sets.Integers
module IntSet = All_sets.Integers


type dyna_state_t = Data.dyna_state_t


(* An implied alignment of raw sequences *)
type ias = {
    seq : Sequence.s;
    codes : (int, int) Hashtbl.t; (* (key=pos -> code) Hashtble *)
    homologous: (int, int Sexpr.t) Hashtbl.t; (* (code, hom_code list) Hashtbl *)
    indels: (int * string * int * [ `Insertion | `Deletion ] * int Sexpr.t) Sexpr.t; 
         (** The location and contents of an insertion block 
          * (p * s * l * t * c) Sexpr.t where
          *  p: start indel position 
          *  s: content of this indel block
          *  l: indel block length 
          *  t: block type either insertion or deletion
          *  c: children of this subtree *)
          
    dum_chars : (int * int Sexpr.t) Sexpr.t;
              (** (c * children) where 
               * c: the dummy character cost
               * children: children of this subtree  *)

    order : int list; (* codes list in reverse order *)
}


(* t is the presentation of a dynamic set (DynamicCS) in order to create implied
   alignments.
   Note: sequences : ias_arr Codes.t where 
   isa_arr is an array of ias in the case of annotated chromosomes (ach ias
   presents a locus)
*)
type t = {
    sequences : ias array Codes.t;
    c2 : Cost_matrix.Two_D.m;
    chrom_pam : Data.dyna_pam_t;
    state : dyna_state_t;
    code : int;
    children : int Sexpr.t; 
    cannonic_code : int;
    alpha : Alphabet.a
}

type pairs = int * t list (* taxon_id * character list *)

module OrderedTuple = struct
    type t = pairs
    let compare (a, _) (b, _) = a - b
end

module AssList = Set.Make (OrderedTuple)

type cg = (unit -> int)


let fprintf = Printf.fprintf

(** [code_generator] returns an unique code for a position in a
* sequence. Thus, codes created by [code_genrator] function
* are all different *)
let code_generator () =
    let counter = ref (-1) in
    fun () ->
        incr counter;
        !counter

(** [create_ias state s code cg] returns an initiated 
* implied alignment for raw sequence [s] where [cg] is 
* code generator function, [state] indicates if sequence [s] 
* is a sequence or a chromosome character.
* The parameter [code] is redundant for now *)
let create_ias (state : dyna_state_t) s code cg =
    let add_codes ((code_acc, hom_acc, order_lst) as res) pos _ =
        if (pos = 0) && (state = `Seq) then res 
        else 
            let code = cg () in
            Hashtbl.add code_acc pos code;
            let single = `Single code in 
            Hashtbl.add hom_acc code single;
            code_acc, hom_acc, code :: order_lst
    in
    let c = Hashtbl.create 1667
    and h = Hashtbl.create 1667 in
    let c, h, o = Sequence.foldi add_codes (c, h, []) s in
    { seq = s; codes = c; homologous = h; indels = `Empty; dum_chars = `Empty; order = o }

let rec prepend_until_shared tgt src it = 
    match src with
    | h :: t when h = it -> tgt, t
    | h :: t -> prepend_until_shared (h :: tgt) t it
    | [] -> failwith "prepend_until_shared"

let rec prepend_all tgt = function
    | h :: t -> prepend_all (h :: tgt) t
    | [] -> tgt

let print_debug a' b' a b m =
    let printem = Status.user_message (Status.Output (Some "ia_dia", false, [])) in
    printem "For implied alignment:\n";
    printem (Sequence.to_formater a' Alphabet.nucleotides );
    printem "\n";
    printem (Sequence.to_formater b' Alphabet.nucleotides);
    printem "\n";
    printem "This is the alignment matrix\n";
    Sequence.Align.print_backtrack a.seq b.seq m;
    printem "\n"

let print_algn_debug = false
let print_anc_debug = false

let calculate_indels a b alph b_children = 
    (* Create a list with all the starting and ending positions of indels, and
    * the indel string, with their homologous beginning and ending base
    * codes *)
    let in_indel_row = ref `None
    and result_list = ref [] in
    let gap = Alphabet.get_gap alph in
    assert ( (Sequence.length a > 1) || (Sequence.get a 0 != gap) ||
                 (Sequence.get b 0 != gap) );
    assert (Sequence.length b = Sequence.length a);


    let len = Sequence.length a in
    let assign_in_indel_row i a_gap b_gap =
        (* Is any of them a gap? *)
        if a_gap && b_gap then 
            in_indel_row := `None
        else if b_gap then 
            in_indel_row := `A i
        else 
            in_indel_row := `B i
    in
    for i = 1 to len - 1 do
        let a_base = Sequence.get a i
        and b_base = Sequence.get b i in
        let a_gap = a_base <> gap (* a_gap = true if a_base is not a gap, otherwise false *)
        and b_gap = b_base <> gap in
        match !in_indel_row with
        | `None -> assign_in_indel_row i a_gap b_gap;
        | `A pos ->
                if a_gap then begin
                    let len = i - pos in
                    let seq = 
                        let seq = Sequence.sub b pos len in
                        Sequence.to_string seq alph 
                    in
                    result_list :=
                        (`Single (pos, seq, len, `Insertion, b_children)) :: 
                            !result_list;
                    assign_in_indel_row i a_gap b_gap;
                end else ()
        | `B pos ->
                if b_gap then begin
                    let len = i - pos in
                    let seq = 
                        let seq = Sequence.sub a pos len in
                        Sequence.to_string seq alph 
                    in
                    result_list :=
                        (`Single (pos, seq, len, `Deletion, b_children)) :: 
                            !result_list;
                    assign_in_indel_row i a_gap b_gap;
                end else ()
    done;
    let _ = 
        match !in_indel_row with
        | `None -> ()
        | `A pos ->
                let len = len - pos in
                let seq = 
                    let seq = Sequence.sub b pos len in
                    Sequence.to_string seq alph 
                    in
                    result_list :=
                        (`Single (pos, seq, len, `Insertion, b_children)) :: 
                            !result_list;
        | `B pos ->
                let len = len - pos in
                let seq = 
                    let seq = Sequence.sub a pos len in
                    Sequence.to_string seq alph 
                in
                result_list :=
                    (`Single (pos, seq, len, `Deletion, b_children)) :: 
                        !result_list;
    in
    match !result_list with
    | [] -> `Empty
    | x -> `Set x

(** [ancestor calc_m a b cm m] creates a common ancestor for sequences [a] and [b]
 * using the cost matrix [cm] and the alignment matrix [m] 
 * The resulting common ancestor holds the homology
 * relationships of the codes assigned in [a] and [b]. 
 * state indicates the type of sequence, i.e., Sequence, Chromosome, Annotated,
 * Genome, Breakinv....
 * [calc_m] if true, will ask the function to automatically calculate the median
 * between the two input sequences as the ancestor, otherwise, it will assume
 * that [a] is the ancestor of [b].
*)
(* TODO CHILDREN *)
let ancestor calculate_median state prealigned all_minus_gap a b 
            codea codeb cm alph achld bchld = 

    if print_anc_debug then
        Status.user_message Status.Information
        ("The ancestors of " ^ string_of_int codea ^ " and " ^ string_of_int codeb);
    let a, b, mincode = 
        if calculate_median then
            if codea < codeb then a, b, codea
            else b, a, codeb
        else a, b, codea
    in
    let lena = Sequence.length a.seq
    and lenb = Sequence.length b.seq 
    and gap = Cost_matrix.Two_D.gap cm in
    let create_gaps len = Sequence.init (fun _ -> gap) len 
    and aempty = (Sequence.is_empty a.seq gap) && (state = `Seq)
    and bempty = (Sequence.is_empty b.seq gap) && (state = `Seq) in
    let a', b', nogap, indels =
        let anb_indels = `Set [a.indels; b.indels] in
        let a', b', _, nogap, indels =
            if aempty && bempty then
                if lena > lenb then a.seq, a.seq, 0, `A, anb_indels
                else b.seq, b.seq, 0, `B, anb_indels
            else if aempty then
                    (create_gaps lenb), b.seq, 0, `A, anb_indels
            else if bempty then
                a.seq, (create_gaps lena), 0, `B, anb_indels
            else begin
                if prealigned then  begin
                    let inds = match lena with
                    | 0 -> `Empty
                    | _ -> calculate_indels a.seq b.seq alph bchld 
                    in
                    a.seq, b.seq, 0, `Both, `Set [inds; anb_indels]
                
                end else
                    let a, b, c = 
                        Sequence.Align.align_2 a.seq b.seq cm
                        Matrix.default 
                    in
                    let inds = 
                        calculate_indels a b alph bchld 
                    in
                    a, b, c, `Both, `Set [inds; anb_indels]
            end 
        in
        let nogap = 
            (* if we are not calculating the median, we better always
            *  pick whatever is assigned to the true ancestor *)
            if calculate_median then nogap
            else `B 
        in
        a', b', nogap, indels
    in
    if print_algn_debug then print_debug a' b' a b Matrix.default;
    let lena' = Sequence.length a' in
    let anc = Sequence.create (lena' + 1) 
    and a_ord = a.order
    and b_ord = b.order in
    let correct_gaps_in_sequences a =
        (* For some of the sequence classes we really don't use the leading gap.
        * So we better clean that up when required *)
        match state with 
        | `Chromosome | `Annotated | `Breakinv | `Genome -> 
                let gapless_seqa = Sequence.delete_gap ~gap_code:gap a.seq in 
                {a with seq = gapless_seqa}, Sequence.length gapless_seqa
        | _ -> a, Sequence.length a.seq
    in 
    let a, lena = correct_gaps_in_sequences a in
    let b, lenb = correct_gaps_in_sequences b in
    let rec builder = 
        fun position a_pos b_pos anc_pos codes hom a_hom b_hom a_or b_or res_or ->
        if position > (-1) then begin
            let it_a = Sequence.get a' position
            and it_b = Sequence.get b' position in
            let med = 
                match nogap with
                | `A -> it_b
                | `B -> it_a
                | `Both -> 
                      match state with 
                      | `Breakinv -> it_b
                      | _ -> Cost_matrix.Two_D.median it_a it_b cm 
            in
            let is_gap_median =
                if it_a <> gap && it_b <> gap then
                    (((Cost_matrix.Two_D.cost it_a it_b cm) <
                    (Cost_matrix.Two_D.cost (all_minus_gap it_a)
                    (all_minus_gap it_b) cm)) || (med = gap))
                else (med = gap)
            in
            let code, hom, n_a_pos, n_b_pos, na_hom, nb_hom, a_or, b_or, res_or =
                match is_gap_median, it_a = gap, it_b = gap with
                | false, false, false ->
                        let codea = 
                            try Hashtbl.find a.codes a_pos with
                            | Not_found ->
                                  failwith 
                                      (Printf.sprintf "Could not find %d with gap
                                     %d and it_a %d and it_b %d\n" a_pos gap it_a it_b)
                        in
                        let codeb = Hashtbl.find b.codes b_pos in
                        let hom_a = Hashtbl.find a_hom codea 
                        and hom_b = Hashtbl.find b_hom codeb in
                        Hashtbl.remove a_hom codea;
                        Hashtbl.remove b_hom codeb;
                        let new_ab = `Set [hom_a; hom_b] in
                        let new_res_or, new_a_or = 
                            prepend_until_shared res_or a_or codea in
                        let new_res_or, new_b_or = 
                            prepend_until_shared new_res_or b_or codeb in
                        Hashtbl.replace hom codea new_ab;
                        codea, hom, a_pos - 1, b_pos - 1, a_hom, 
                        b_hom, new_a_or, new_b_or, (codea :: new_res_or)
                | true, false, false ->
                        let codeb = Hashtbl.find b.codes b_pos 
                        and codea = Hashtbl.find a.codes a_pos in
                        let hom_b = Hashtbl.find b_hom codeb 
                        and hom_a = Hashtbl.find a_hom codea in
                        Hashtbl.remove b_hom codeb;
                        Hashtbl.remove a_hom codea;
                        let new_res_or, new_b_or =
                            prepend_until_shared res_or b_or codeb 
                        in
                        let new_res_or, new_a_or = 
                            prepend_until_shared new_res_or a_or codea 
                        in
                        Hashtbl.add hom codeb hom_b;
                        Hashtbl.add hom codea hom_a;
                        let new_res_or = (codea :: codeb :: new_res_or) in
                        codea, hom, a_pos - 1, b_pos - 1, a_hom, b_hom,
                        new_a_or, new_b_or, new_res_or
                | _, true, false ->
                        let codeb = Hashtbl.find b.codes b_pos in
                        let hom_b = Hashtbl.find b_hom codeb in
                        Hashtbl.remove b_hom codeb;
                        Hashtbl.replace hom codeb hom_b;
                        let new_res_or, new_b_or =
                            prepend_until_shared res_or b_or codeb in
                        codeb, hom, a_pos, b_pos - 1, a_hom, b_hom,
                        a_or, new_b_or, (codeb :: new_res_or)
                | _, false, true ->
                        let codea = Hashtbl.find a.codes a_pos in
                        let hom_a = Hashtbl.find a_hom codea in
                        Hashtbl.remove a_hom codea;
                        Hashtbl.replace hom codea hom_a;
                        let new_res_or, new_a_or = 
                            prepend_until_shared res_or a_or codea in
                        codea, hom, a_pos - 1, b_pos, a_hom, b_hom,
                        new_a_or, b_or, (codea :: new_res_or)
                | _, true, true ->
(*                        assert (a_pos = 0);*)
(*                        assert (b_pos = 0);*)
                        (-1), hom, a_pos, b_pos, a_hom, b_hom, a_or, b_or,
                        res_or
            in
            let n_anc_pos = 
                if not is_gap_median then begin
                    Sequence.prepend anc med;
                    Hashtbl.replace codes anc_pos code;
                    anc_pos + 1
                end else anc_pos
            in
            
(*            List.iter (fun ord -> fprintf stdout "%i " ord) res_or;
            print_newline (); *)
            builder (position - 1) n_a_pos n_b_pos n_anc_pos codes hom na_hom
            nb_hom a_or b_or res_or 
        end else begin
            (* We have to prepend a gap to the
            *  ancestor, though we don't include it in the set of homologies. *)
            Sequence.prepend anc gap; 
            (* The codes of the ancestor positions are shifted and are wrong,
            * depending on weather or not the median in a particular position is
            * a gap. And so, we fix the codes *)
            let codes =
                let res = Hashtbl.create 1667 in
                Hashtbl.iter 
                (fun k c -> Hashtbl.add res (anc_pos - k) c) codes;
                res
            in
            (* Now we add the homologies that are in the tree above us, after
            * the current sequence. *)
            Hashtbl.iter (Hashtbl.replace a_hom) hom;
            Hashtbl.iter (Hashtbl.replace a_hom) b_hom;
            let new_res_or = prepend_all res_or a_or in
            let new_res_or = prepend_all new_res_or b_or in
            codes, a_hom, new_res_or
        end
    in
    let initial_codes = Hashtbl.create 1667
    and initial_hom = Hashtbl.create 1667 in
    let codes, hom, order = 
        builder (lena' - 1) (lena - 1) (lenb - 1) 0 initial_codes
        initial_hom a.homologous b.homologous a_ord b_ord []
    in
    let order = List.rev order in
    { seq = anc; codes = codes; homologous = hom; indels = indels; dum_chars = `Empty; order = order }


let assign_act_order sta en codes ord_arr act_ord_arr =
    if sta != -1 then begin
        let sta_code = Hashtbl.find codes sta in 
        let en_code = Hashtbl.find codes en in 
        let sta = Utl.find_index ord_arr sta_code compare in 
        let en = Utl.find_index ord_arr en_code compare in 
        Array.fill act_ord_arr sta (en - sta + 1) true;
    end else ()
        

let get_suborder is_main_chrom sta en codes ord_arr act_ord_arr =
    if sta = -1 then [||] 
    else  begin
        let sta_code = Hashtbl.find codes sta in 
        let en_code = Hashtbl.find codes en in 
        let sta = Utl.find_index ord_arr sta_code compare in 
        let en = Utl.find_index ord_arr en_code compare in 

        let sta = ref sta in
        while is_main_chrom &&(!sta > 0) && (act_ord_arr.(!sta - 1) = false) do
            act_ord_arr.(!sta - 1) <- true;
            sta := !sta - 1;
        done;
        
        let num_ord = Array.length ord_arr in 
        let en = ref en in
        while is_main_chrom && (!en + 1 < num_ord) && (act_ord_arr.(!en + 1) = false) do
            act_ord_arr.(!en + 1) <- true;
            en := !en + 1;
        done;
        
        let sub_ord_arr = Array.sub ord_arr !sta (!en - !sta + 1) in 
        sub_ord_arr
    end 

let ancestor_sequence prealigned calculate_median all_minus_gap 
        acode bcode achld bchld a_ls b_ls cm alpha chrom_pam =
    let get_one x =
        match x with
        | [|x|] -> x
        | _ -> assert false
    in
    let res =
        ancestor calculate_median `Seq prealigned all_minus_gap
        (get_one a_ls) (get_one b_ls) acode bcode cm alpha achld bchld
    in
    [|res|]


(* [ancestor_chrom prealigned calculate_median all_minus_gap acode bcode 
*  achld bchld a b cm alpha chrom_pam] merges the implied alignments of two clades and 
* their respective roots into one common ancestor. This is similar to [ancestor]
* function, but for chromosome characters  *)
let ancestor_chrom prealigned calculate_median all_minus_gap acode bcode 
        achld bchld a b cm alpha chrom_pam = 
    let a, b, min_can_code = 
        if calculate_median then 
            if acode < bcode then a, b, acode
            else b, a, bcode
        else a, b, acode
    in 
    let a = a.(0) and b = b.(0) in 
    let lena = Sequence.length a.seq
    and lenb = Sequence.length b.seq 
    and gap = Cost_matrix.Two_D.gap cm in
    let create_gaps len = Sequence.init (fun _ -> gap) len 
    and aempty = Sequence.is_empty a.seq gap
    and bempty = Sequence.is_empty b.seq gap in
    let seq_a, seq_b =
        if aempty && bempty then
            if lena > lenb then a.seq, a.seq
            else b.seq, b.seq
        else if aempty then
            (create_gaps lenb), b.seq
        else if bempty then
            a.seq, (create_gaps lena)
        else a.seq, b.seq
    in  
    let med_a = ChromAli.create_med seq_a in 
    let med_b = ChromAli.create_med seq_b in 
    let chrom_pam = {chrom_pam with Data.approx = Some true} in 
    let _, _, med_ls = ChromAli.find_med2_ls med_a med_b cm chrom_pam in 
    let med = List.hd med_ls in 

    let ordera_arr = Array.of_list (List.rev a.order) in 
    let orderb_arr = Array.of_list (List.rev b.order) in 

    let act_orda_arr = Array.map (fun _ -> false) ordera_arr in
    let act_ordb_arr = Array.map (fun _ -> false) orderb_arr in        

    List.iter 
        (fun seg ->
             let staa = seg.ChromAli.sta1 in 
             let ena = seg.ChromAli.en1 in                  
             let stab = seg.ChromAli.sta2 in 
             let enb = seg.ChromAli.en2 in 
             (if seg.ChromAli.dir2 = `Negative then 
                  failwith 
                      "The impliedAlignment function has 
                       not handled the negative segment");
             assign_act_order staa ena a.codes ordera_arr act_orda_arr;
             assign_act_order stab enb b.codes orderb_arr act_ordb_arr
        ) med.ChromAli.chrom_map;

    let init_ias = {seq = med.ChromAli.seq; 
                    codes = Hashtbl.create 1667;
                    homologous = Hashtbl.create 1667;
                    order = [];
                    indels = `Empty;
                    dum_chars = `Empty}
    in 
    
    let added_locus_indel_cost = ref 0 in 
    let num_locus_indel = ref 0 in
    let new_ias, rev_anc_seq_ls = List.fold_left 
        (fun (nascent_ias, rev_anc_seq_ls) seg ->
             let sta = seg.ChromAli.sta in 
             let staa = seg.ChromAli.sta1 in 
             let ena = seg.ChromAli.en1 in 
             let sub_orda_arr = get_suborder true staa ena a.codes ordera_arr act_orda_arr in 
             let stab = seg.ChromAli.sta2 in 
             let enb = seg.ChromAli.en2 in 
             let sub_ordb_arr = get_suborder true stab enb b.codes orderb_arr act_ordb_arr in 

             let sub_codes_a = Hashtbl.create 1667 in
             (if staa != -1 then 
                  for p = staa to ena do
                      let code = Hashtbl.find a.codes p in 
                      Hashtbl.add sub_codes_a (p - staa) code;
                  done);
             let sub_codes_b = Hashtbl.create 1667 in
             (if (stab != -1) then 
                  for p = stab to enb do
                      let code = Hashtbl.find b.codes p in 
                      Hashtbl.add sub_codes_b (p - stab) code;
                  done);
                  
             let sub_a = {a with 
                              seq = seg.ChromAli.alied_seq1;
                              codes = sub_codes_a;
                              order = List.rev (Array.to_list sub_orda_arr);                                  
                         } 
             in 
             let sub_b = {b with 
                              seq =  seg.ChromAli.alied_seq2;
                              codes = sub_codes_b;
                              order = List.rev (Array.to_list sub_ordb_arr);
                         } 
             in  
             let ans = 
                 ancestor calculate_median `Chromosome prealigned all_minus_gap sub_a
                     sub_b acode bcode cm alpha achld bchld
             in 
             (if (sta != -1) then 
                  Hashtbl.iter (fun p code ->                                     
                                    Hashtbl.add nascent_ias.codes (p + sta - 1) code) ans.codes);


             Hashtbl.iter (fun code hom -> 
                               Hashtbl.replace nascent_ias.homologous code hom) ans.homologous; 
 
             (if ( (staa = -1) || (stab = -1) ) then begin

                let gap_cost = Sequence.cmp_ali_cost 
                    seg.ChromAli.alied_seq1 seg.ChromAli.alied_seq2 seg.ChromAli.dir2 cm
                in 
                incr num_locus_indel;
                added_locus_indel_cost := !added_locus_indel_cost + (seg.ChromAli.cost-gap_cost);
             end);

             let nascent_ias = {nascent_ias with order = List.append nascent_ias.order (List.rev ans.order)} in
             (nascent_ias, ans.seq::rev_anc_seq_ls)
        ) (init_ias, []) med.ChromAli.chrom_map
    in 

    
    let recost1 = med.ChromAli.recost1 in
    let recost2 = med.ChromAli.recost2 in

    let rea = `Single (recost1, achld) in 
    let reb = `Single (recost2, bchld) in
    let locus_indel = 
        if (List.hd (Sexpr.to_list achld)) < (List.hd (Sexpr.to_list bchld)) then 
            `Single (!added_locus_indel_cost,  achld) 
        else
            `Single (!added_locus_indel_cost,  bchld) 
    in
    
    let dum_chars = ref [a.dum_chars; b.dum_chars] in 
    (if recost1 > 0 then dum_chars:= rea::!dum_chars);
    (if recost2 > 0 then dum_chars:= reb::!dum_chars);
    (if !num_locus_indel > 0 then dum_chars:= locus_indel::!dum_chars);

    let anc_chrom = Sequence.concat (List.rev rev_anc_seq_ls) in
    let anc_chrom = Sequence.delete_gap anc_chrom in 
    let new_ias = {new_ias with order = List.rev new_ias.order;
                                seq = anc_chrom;
                                dum_chars = `Set !dum_chars} in     


    [|new_ias|]



(* [ancestor_annchrom prealigned calculate_median all_minus_gap acode bcode 
*  achld bchld a b cm alpha annchrom_pam] merges the implied alignments of two clades and 
* their respective roots into one common ancestor. This is similar to [ancestor]
* function, but for annotated chromosome characters  *)
let ancestor_annchrom prealigned calculate_median all_minus_gap acode bcode
        achld bchld a b cm alpha annchrom_pam  = 

    let a, b, min_can_code = 
        if calculate_median then 
            if acode < bcode then a, b, acode
            else b, a, bcode
        else a, b, acode
    in 
    let seqa_arr = Array.map (fun ias -> (ias.seq, -1)) a in
    let seqb_arr = Array.map (fun ias -> (ias.seq, -1)) b in
    let ana = AnnchromAli.init seqa_arr in 
    let anb = AnnchromAli.init seqb_arr in 
    let annchrom_pam = {annchrom_pam with Data.approx = Some true} in 
    let _, _, med_ls = 
        AnnchromAli.find_med2_ls ana anb cm Alphabet.nucleotides annchrom_pam 
    in
    let med = List.hd med_ls in 

    let added_locus_indel_cost = ref 0 in 
    let num_locus_indel = ref 0 in
    let merge_ias seqt = 
        let orda = seqt.AnnchromAli.seq_ord1 in 
        let ordb = seqt.AnnchromAli.seq_ord2 in 
        let isa = 
            match orda = -1 with
            | false -> {a.(orda) with 
                            seq = seqt.AnnchromAli.alied_seq1 }
            | true -> 
                    { seq = seqt.AnnchromAli.alied_seq1;
                    codes = Hashtbl.create 1667;
                    homologous = Hashtbl.create 1667;
                    order = [];
                    indels = `Empty;
                    dum_chars = `Empty}
        in 
        let isb = 
            match ordb = -1 with
            | false -> {b.(ordb) with 
                            seq = seqt.AnnchromAli.alied_seq2 }
            | true -> {seq = seqt.AnnchromAli.alied_seq2;
                       codes = Hashtbl.create 1667;
                       homologous = Hashtbl.create 1667;
                       order = []; indels = `Empty; dum_chars = `Empty}
        in 
        let ans = 
            ancestor calculate_median `Annotated prealigned
            all_minus_gap isa isb acode bcode cm alpha achld bchld
        in
        (if (orda = -1) || (ordb = -1) then begin
            let gap_cost = Sequence.cmp_ali_cost 
                seqt.AnnchromAli.alied_seq1 seqt.AnnchromAli.alied_seq2 `Positive cm
            in 

            let indel_cost = Data.get_locus_indel_cost annchrom_pam in 
            let seq_cost = 
                if orda = -1 then 
                    Sequence.cmp_gap_cost indel_cost seqt.AnnchromAli.alied_seq2
                else 
                    Sequence.cmp_gap_cost indel_cost seqt.AnnchromAli.alied_seq1
            in
                    
            incr num_locus_indel;
            added_locus_indel_cost := !added_locus_indel_cost + (seq_cost - gap_cost);
        end);


        let new_codes = Hashtbl.create 1667 in 
        Hashtbl.iter (fun p code -> Hashtbl.add new_codes (p-1) code) ans.codes;
        {ans with seq = seqt.AnnchromAli.seq; codes = new_codes}
    in 

    let new_ias_arr = Array.map merge_ias med.AnnchromAli.seq_arr in 

    let recost1 = med.AnnchromAli.recost1 in
    let recost2 = med.AnnchromAli.recost2 in

    let rea = `Single (recost1, achld) in 
    let reb = `Single (recost2, bchld) in
    let locus_indel = 
        if (List.hd (Sexpr.to_list achld)) < (List.hd (Sexpr.to_list bchld)) then 
            `Single (!added_locus_indel_cost,  achld) 
        else
            `Single (!added_locus_indel_cost,  bchld) 
    in
    

    let dum_chars = ref [a.(0).dum_chars; b.(0).dum_chars] in 
    (if recost1 > 0 then dum_chars:= rea::!dum_chars);
    (if recost2 > 0 then dum_chars:= reb::!dum_chars);
    (if !num_locus_indel > 0 then dum_chars:= locus_indel::!dum_chars);
    new_ias_arr.(0) <- {new_ias_arr.(0) with dum_chars = `Set !dum_chars}; 

    new_ias_arr



(* [ancestor_breakinv prealigned calculate_median all_minus_gap acode bcode 
*  achld bchld a b cm alpha breakinv_pam] merges the implied alignments of two clades and 
* their respective roots into one common ancestor. This is similar to [ancestor]
* function, but for breakinv characters  *)
let ancestor_breakinv prealigned calculate_median all_minus_gap acode bcode
        achld bchld a b cm alpha breakinv_pam = 

    let a, b, min_can_code = 
        if calculate_median then 
            if acode < bcode then a, b, acode
            else b, a, bcode
        else a, b, acode
    in 
    let gap_code = Cost_matrix.Two_D.gap cm in 

    let a, b = a.(0), b.(0) in 
    let meda = BreakinvAli.init a.seq in 
    let medb = BreakinvAli.init b.seq in 
    let pure_cm = Cost_matrix.Two_D.get_pure_cost_mat cm in 
    let _, _, med_ls = BreakinvAli.find_med2_ls meda medb cm pure_cm alpha breakinv_pam in
    let med = List.hd med_ls in


    let seqa_arr = Sequence.to_array a.seq in 
    let re_seqa = Sequence.delete_gap ~gap_code:gap_code med.BreakinvAli.alied_seq1 in
    let new_codes_a = Hashtbl.create 1667 in 
    Array.iteri 
    (fun re_pos re_base -> 
        let ori_pos = Utl.find_index seqa_arr re_base compare in 
        let code = Hashtbl.find a.codes ori_pos in         
        Hashtbl.add new_codes_a re_pos code) (Sequence.to_array re_seqa);

    let order_a_arr = Array.of_list (List.rev a.order) in 
    let num_codes = Array.length order_a_arr in 
    let exist_code_a = Hashtbl.fold 
        (fun p code code_set -> IntSet.add code code_set) new_codes_a IntSet.empty
    in 
    let new_orders_a = ref [] in 
    let rec add_deled_code pos = 
        if pos < num_codes then begin
            let code = order_a_arr.(pos) in 
            if IntSet.mem code exist_code_a = false then begin
                new_orders_a := List.append !new_orders_a [code];
                add_deled_code (pos + 1)
            end 
        end            
    in 
    Array.iteri 
    (fun pos _ -> 
        let code = Hashtbl.find new_codes_a pos in 
        let ori_pos = Utl.find_index order_a_arr code compare in 
        new_orders_a := List.append !new_orders_a [code];                     
        add_deled_code (ori_pos + 1)) 
    seqa_arr;
    add_deled_code 0;



    let seqb_arr = Sequence.to_array b.seq in 
    let re_seqb = Sequence.delete_gap ~gap_code:gap_code med.BreakinvAli.alied_seq2 in
    let new_codes_b = Hashtbl.create 1667 in 
    Array.iteri 
    (fun re_pos re_base -> 
        let ori_pos = Utl.find_index seqb_arr re_base compare in 
        let code = Hashtbl.find b.codes ori_pos in         
        Hashtbl.add new_codes_b re_pos code) (Sequence.to_array re_seqb);

    let order_b_arr = Array.of_list (List.rev b.order) in 
    let num_codes = Array.length order_b_arr in 
    let exist_code_b = Hashtbl.fold 
        (fun p code code_set -> IntSet.add code code_set) new_codes_b IntSet.empty
    in 
    let new_orders_b = ref [] in 
    let rec add_deled_code pos = 
        if pos < num_codes then begin
            let code = order_b_arr.(pos) in 
            if IntSet.mem code exist_code_b = false then begin
                new_orders_b := List.append !new_orders_b [code];
                add_deled_code (pos + 1)
            end 
        end            
    in 
    Array.iteri 
    (fun pos _ -> 
        let code = Hashtbl.find new_codes_b pos in 
        let ori_pos = Utl.find_index order_b_arr code compare in 
        new_orders_b := List.append !new_orders_b [code];                     
        add_deled_code (ori_pos + 1)) 
    seqb_arr;
    add_deled_code 0;


    let isa = {a with seq = med.BreakinvAli.alied_seq1;
                   codes = new_codes_a;
                   order = List.rev  !new_orders_a} in 

    let isb = {b with seq = med.BreakinvAli.alied_seq2;
                   codes = new_codes_b;
                   order = List.rev  !new_orders_b} in 
    let ans = 
        ancestor calculate_median `Breakinv prealigned all_minus_gap isa isb
        acode bcode cm alpha achld bchld
    in


    let recost1 = med.BreakinvAli.recost1 in
    let recost2 = med.BreakinvAli.recost2 in

    let rea = `Single (recost1, achld) in 
    let reb = `Single (recost2, bchld) in
    
    let dum_chars = ref [a.dum_chars; b.dum_chars] in 
    (if recost1 > 0 then dum_chars:= rea::!dum_chars);
    (if recost2 > 0 then dum_chars:= reb::!dum_chars);

    let ans = {ans with dum_chars = `Set !dum_chars} in

    let new_codes = Hashtbl.create 1667 in 
    Hashtbl.iter (fun p code -> Hashtbl.add new_codes (p-1) code) ans.codes;
    [|{ans with seq = med.BreakinvAli.seq; codes = new_codes}|]


(* [ancestor_genome prealigned calculate_median all_minus_gap acode bcode 
*  achld bchld a b cm alpha chrom_pam] merges the implied alignments of two clades and 
* their respective roots into one common ancestor. This is similar to [ancestor]
* function, but for genome characters  *)
let ancestor_genome prealigned calculate_median all_minus_gap acode bcode achld
        bchld a b cm alpha chrom_pam = 

    let ias1_arr, ias2_arr, min_can_code, max_can_code =
        if acode < bcode then a, b, acode, bcode
        else b, a, bcode, acode
    in
    let chroma_arr = Array.map (fun ias -> ias.seq) ias1_arr in
    let chromb_arr = Array.map (fun ias -> ias.seq) ias2_arr in
    let med1 = GenomeAli.create_med_from_seq chroma_arr in 
    let med2 = GenomeAli.create_med_from_seq chromb_arr in
    let _, _, med_ls = GenomeAli.find_med2_ls med1 med2 cm chrom_pam in 
    let med = List.hd med_ls in 

    let order1_mat = Array.map (fun ias1 -> Array.of_list (List.rev ias1.order) ) ias1_arr in  
    let order2_mat = Array.map (fun ias2 -> Array.of_list (List.rev ias2.order) ) ias2_arr in  
    let act_ord1 = Array.map (fun ord_arr -> Array.map (fun _ -> false) ord_arr) order1_mat in
    let act_ord2 = Array.map (fun ord_arr -> Array.map (fun _ -> false) ord_arr) order2_mat in
    let id_to_index chrom_id med = 
        let chrom_id_arr = Array.map (fun chrom -> !(chrom.GenomeAli.chrom_id)) 
            med.GenomeAli.chrom_arr
        in 
        Utl.find_index chrom_id_arr chrom_id compare
    in 
    Array.iter 
        (fun chrom -> 
             List.iter 
                 (fun seg -> 
                      let sta1 = seg.GenomeAli.sta1 in 
                      let en1 = seg.GenomeAli.en1 in 
                      let idx1 = id_to_index seg.GenomeAli.chi1_chrom_id med1 in 
                      (if idx1 != -1 then 
                          assign_act_order sta1 en1 ias1_arr.(idx1).codes order1_mat.(idx1) act_ord1.(idx1));
                      let sta2 = seg.GenomeAli.sta2 in 
                      let en2 = seg.GenomeAli.en2 in 
                      let idx2 = id_to_index seg.GenomeAli.chi2_chrom_id med2 in 
                      (if idx2 != -1 then 
                          assign_act_order sta2 en2 ias2_arr.(idx2).codes order2_mat.(idx2) act_ord2.(idx2));
                 ) chrom.GenomeAli.map
        ) med.GenomeAli.chrom_arr; 

    let added_locus_indel_cost = ref 0 in 
    let num_locus_indel = ref 0 in

    let new_ias_arr = Array.map 
        (fun chrom ->
             let init_ias = {seq = Sequence.get_empty_seq (); 
                             codes = Hashtbl.create 1667;
                             homologous = Hashtbl.create 1667;
                             indels = `Empty; dum_chars = `Empty;
                             order = []}
             in 
             let main_idx1 = id_to_index chrom.GenomeAli.main_chrom1_id med1 in 
             let main_idx2 = id_to_index chrom.GenomeAli.main_chrom2_id med2 in 
             let med_len = ref 0 in 
             let new_ias, rev_anc_seq_ls = List.fold_left
                 ( fun (nascent_ias, rev_anc_seq_ls) seg ->
                       let sta = seg.GenomeAli.sta in 
                       let sta1 = seg.GenomeAli.sta1 in 
                       let en1 = seg.GenomeAli.en1 in 
                       let idx1 = id_to_index seg.GenomeAli.chi1_chrom_id med1 in 
                       let sub_ord1_arr = match idx1 = -1 with
                       | true -> [||]
                       | false ->
                             get_suborder (idx1 = main_idx1) sta1 en1 
                                 ias1_arr.(idx1).codes order1_mat.(idx1) act_ord1.(idx1)
                       in 
                       
                       let sta2 = seg.GenomeAli.sta2 in 
                       let en2 = seg.GenomeAli.en2 in 
                       let idx2 = id_to_index seg.GenomeAli.chi2_chrom_id med2 in 
                       let sub_ord2_arr = match idx2 = -1 with
                       | true -> [||]
                       | false ->
                             get_suborder (idx2 = main_idx2) sta2 en2 
                                 ias2_arr.(idx2).codes order2_mat.(idx2) act_ord2.(idx2)
                       in 
                       
                       let sub_codes1 = Hashtbl.create 1667 in
                       (if sta1 != -1 then 
                            for p = sta1 to en1 do
                                let code = Hashtbl.find ias1_arr.(idx1).codes p in 
                                Hashtbl.add sub_codes1 (p - sta1) code;
                            done);
                       let sub_codes2 = Hashtbl.create 1667 in
                       (if sta2 != -1 then 
                            for p = sta2 to en2 do
                                let code = Hashtbl.find ias2_arr.(idx2).codes p in 
                                Hashtbl.add sub_codes2 (p - sta2) code;
                            done);

                       (if (sta1 = -1) || (sta2 = -1) then begin
                                let gap_cost = Sequence.cmp_ali_cost 
                                seg.GenomeAli.alied_seq1 seg.GenomeAli.alied_seq2 seg.GenomeAli.dir2 cm
                            in 
                            incr num_locus_indel;
                            added_locus_indel_cost := !added_locus_indel_cost + (seg.GenomeAli.cost-gap_cost);
                        end);

                       let hom1 = match idx1 = -1 with
                       | true -> Hashtbl.create 1668 
                       | false -> ias1_arr.(idx1).homologous
                       in 
                       let hom2 = match idx2 = -1 with
                       | true -> Hashtbl.create 1668 
                       | false -> ias2_arr.(idx2).homologous
                       in 
                       
                       let sub1 = {
                           seq = seg.GenomeAli.alied_seq1;
                           codes = sub_codes1;
                           order = List.rev (Array.to_list sub_ord1_arr);                                  
                           homologous = hom1;
                            indels = `Empty;
                            dum_chars = `Empty
                        }  
                       in 
                       let sub2 = {
                           seq = seg.GenomeAli.alied_seq2;
                           codes = sub_codes2;
                           order = List.rev (Array.to_list sub_ord2_arr);                                  
                           homologous = hom2;
                            indels = `Empty;
                            dum_chars = `Empty
                       }  
                       in  
                       let ans_ias = 
                           ancestor calculate_median 
                           `Genome prealigned all_minus_gap
                           sub1 sub2 acode bcode cm alpha achld bchld
                       in 
                       (if (sta != -1) then 
                            Hashtbl.iter 
                                (fun p code -> 
                                     Hashtbl.add nascent_ias.codes (p + !med_len - 1) code
                                ) ans_ias.codes); 
                       Hashtbl.iter 
                           (fun code hom -> 
                                Hashtbl.replace nascent_ias.homologous code hom
                           ) ans_ias.homologous; 
                       let nascent_ias = 
                            {nascent_ias with order = 
                               List.append nascent_ias.order (List.rev ans_ias.order)}                 
                       in
                       med_len := !med_len + (Sequence.length ans_ias.seq) - 1;
                       (nascent_ias, ans_ias.seq::rev_anc_seq_ls) 
                 ) (init_ias, []) chrom.GenomeAli.map
             in 
             let chrom_seq = Sequence.concat (List.rev rev_anc_seq_ls) in
             let chrom_seq = Sequence.delete_gap chrom_seq in
             {new_ias with order = List.rev new_ias.order;
                seq = chrom_seq}

        ) med.GenomeAli.chrom_arr
    in 

    let recost1 = med.GenomeAli.recost1 in
    let recost2 = med.GenomeAli.recost2 in

    let rea = `Single (recost1, achld) in 
    let reb = `Single (recost2, bchld) in
    let locus_indel = 
        if (List.hd (Sexpr.to_list achld)) < (List.hd (Sexpr.to_list bchld)) then 
            `Single (!added_locus_indel_cost,  achld) 
        else
            `Single (!added_locus_indel_cost,  bchld) 
    in
    
    let dum_chars = ref [a.(0).dum_chars; b.(0).dum_chars] in 
    (if recost1 > 0 then dum_chars:= rea::!dum_chars);
    (if recost2 > 0 then dum_chars:= reb::!dum_chars);
    (if !num_locus_indel > 0 then dum_chars:= locus_indel::!dum_chars);
    new_ias_arr.(0) <- {new_ias_arr.(0) with dum_chars = `Set !dum_chars}; 

    new_ias_arr


exception IsSankoff

type matrix_class = 
    | AllOne of int
    | AllOneGapSame of (int * int)
    | AffinePartition of (int * int * int)
    | AllSankoff of (string -> int) option 
    (* If using affine gap cost or non metric tcm (where gaps and substitutions
    * need to be split), we pass a function to compute
    * the cost of an indel block to deal with affine. *)

let present_absent_alph = 
    Alphabet.list_to_a 
    [("present", 1, None); ("absent", 2, None)] 
    "absent" None Alphabet.Sequential

(* A function that analyzes a cost matrix and an alphabet and
* generates a pair of functions f and g, such that f converts 
* a state into a list of character states, and g converts a state into it's
* appropriate Parser.Hennig.Encoding.s *)
let analyze_tcm tcm alph =
    let gap = Alphabet.get_gap alph 
    and all = Alphabet.get_all alph in
    let for_sankoff =
        let go = 
            match Cost_matrix.Two_D.affine tcm with
            | Cost_matrix.No_Alignment 
            | Cost_matrix.Linnear -> 0
            | Cost_matrix.Affine go -> go
        in
        let to_string string pos = String.make 1 string.[pos] in
        let rec processor pos max cost string =
            if max = pos then cost
            else 
                let base = Alphabet.find_base (to_string string pos) alph in
                processor (pos + 1) max 
                ((Cost_matrix.Two_D.cost gap base tcm) + cost) string
        in
        if Cost_matrix.Two_D.is_metric tcm then
            fun _ -> go
        else
            fun string -> processor 0 (String.length string) go string
    in
    let alph = Alphabet.simplify alph in
    let single_compare (_, a) res (_, b) =
        match res with
        | None -> 
                let cost = (Cost_matrix.Two_D.cost a b tcm) in
                Some cost
        | Some y ->
                let x = Cost_matrix.Two_D.cost a b tcm in
                if x = y then res
                else raise IsSankoff
    in
    let rec compare_costs l1 l2 res =
        match l1, l2 with
        | _, []
        | [], _ -> res
        | h1 :: t1, _ :: t2 ->
                compare_costs t1 t2 (List.fold_left (single_compare h1) res l2)
    in
    let check_x x =
        match all with
        | Some y -> y <> x
        | None -> true
    in
    let get_cost_of_all_subs () =
        match 
            List.filter (fun (_, x) -> (x <> gap) && (check_x x))
            (Alphabet.to_list alph) 
        with
        | [] -> failwith "An empty alphabet?"
        | (_ :: t) as res ->
                match compare_costs res t None with
                | Some v -> v
                | None -> failwith "No costs?"
    in
    let get_cost_of_gap () =
        match 
            List.filter (fun (_, x) -> (x <> gap) && (check_x x))
            (Alphabet.to_list alph) 
        with
        | [] -> failwith "An empty alphabet?"
        | res ->
                match compare_costs [("", gap)] res None with
                | Some v -> v
                | None -> failwith "No costs?"
    in
    let get_gap_opening tcm =
        match Cost_matrix.Two_D.affine tcm with
        | Cost_matrix.No_Alignment 
        | Cost_matrix.Linnear -> failwith "not affine"
        | Cost_matrix.Affine go -> go
    in
    let is_affine tcm =
        match Cost_matrix.Two_D.affine tcm with
        | Cost_matrix.Affine _ -> true
        | _ -> false
    in
    let all_same_affine () =
        try let _ = get_gap_opening tcm in true with
        | _ -> false
    in
    let get_case =
        try
            let all_excepting_gap = get_cost_of_all_subs ()
            and all_and_gap = get_cost_of_gap () in
            match Alphabet.kind alph with
            | Alphabet.Simple_Bit_Flags ->
                    if 32 > Alphabet.distinct_size alph then
                        if all_same_affine () then
                            AffinePartition 
                            (all_excepting_gap, all_and_gap, 
                            get_gap_opening tcm)
                        else if all_excepting_gap = all_and_gap then 
                            AllOne all_excepting_gap
                        else if not (is_affine tcm) then
                            AllOneGapSame 
                            (all_excepting_gap, all_and_gap)
                        else 
                            AllSankoff (Some for_sankoff)
                    else if is_affine tcm then
                        AllSankoff (Some for_sankoff)
                    else AllSankoff None
            | _ -> 
                    if is_affine tcm then 
                        AllSankoff (Some for_sankoff)
                    else AllSankoff None
        with
        | IsSankoff -> 
                    AllSankoff (Some for_sankoff)
    in
    let extract_all all =
        match all with
        | Some all -> all
        | None -> assert false
    in
    let table = Hashtbl.create 67 in
    let find_item it =
        if Hashtbl.mem table it then
            Hashtbl.find table it
        else begin
            let r = Parser.Unordered_Character (it, false) in
            Hashtbl.add table it r;
            r
        end
    in
    match get_case with
    | AllOne weight ->
            (* We assume that we have dna sequences *)
            let all = extract_all all in
            let encoding = 
                alph,
                (Parser.OldHennig.Encoding.set_weight
                Parser.OldHennig.Encoding.dna_encoding weight)
            in
            let to_parser is_missing states acc = 
                match is_missing, states with
                | `Missing, _ -> (find_item all) :: acc
                | `Exists, 0 -> (find_item gap) :: acc
                | `Exists, x -> (find_item x) :: acc
            and to_encoding _ acc = encoding :: acc in
            get_case, to_parser, to_encoding
    | AllOneGapSame (subsc, gapcost) ->
            let present_absent = 
                present_absent_alph,
                (Parser.OldHennig.Encoding.gap_encoding gapcost)
            and subs = 
                alph,
                (Parser.OldHennig.Encoding.set_weight
                Parser.OldHennig.Encoding.dna_encoding subsc)
            in
            let notgap = lnot gap in
            (* We assume we have dna sequences *)
            let all = notgap land (extract_all all) in
            let to_parser is_missing states acc =
                match is_missing, states with
                | `Missing, _ ->
                        (find_item all) :: (find_item (1 lor 2)) :: acc
                | `Exists, 0 ->
                        (* All characters, and the gap itself, in other words,
                        * we treat the gap as a separate character, and the
                        * state as missing data *)
                        (find_item all) :: (find_item 1) :: acc
                | `Exists, x ->
                        let r = if x = all then 3 else 2 in
                        (find_item (x land notgap)) :: (find_item r) :: acc
            and to_encoding _ acc = 
                subs :: present_absent :: acc
            in
            get_case, to_parser, to_encoding
    | AffinePartition (subsc, gapcost, gapopening) ->
            (* We have to partition the column in three columns, each
            * corresponding to gap opening, gap extension, and substitution.
            * We will have to filter out columns that are not gap opening
            * but only extension.
            * *)
            let subs = 
                alph,
                (Parser.OldHennig.Encoding.set_weight
                Parser.OldHennig.Encoding.dna_encoding
                subsc)
            in
            let notgap = lnot gap in
            let all = notgap land (extract_all all) in
            let to_parser is_missing states acc =
                match is_missing, states with
                | `Missing, _ 
                | `Exists, 0 -> 
                        (* We have a gap, so we assign both gap opening and
                        * gap extension, we will later cleaunup when gap
                        * opening is not needed *)
                        (find_item all) :: acc
                | `Exists, x -> (find_item (x land notgap)) :: acc
            in 
            let to_encoding _ acc = subs :: acc in
            get_case, to_parser, to_encoding
    | AllSankoff _ ->
            let size = 
                (* We remove one from the all elements representation *)
                match Alphabet.get_all alph with
                | Some _ -> (Alphabet.distinct_size alph) - 1 
                | None -> 
                      Alphabet.distinct_size alph
            in
            let is_metric = Cost_matrix.Two_D.is_metric tcm in
            let make_tcm () =
                let size = 
                    (* We will assume that the source of the non metricity is
                    * the gap *)
                    if is_metric then size else size - 1
                in
                match Alphabet.kind alph with
                | Alphabet.Simple_Bit_Flags ->
                        Array.init size (fun x -> Array.init size 
                        (fun y -> 
                            Cost_matrix.Two_D.cost (1 lsl x) (1 lsl y) tcm)) 
                | Alphabet.Sequential ->
                      let tcm_size = Cost_matrix.Two_D.alphabet_size tcm in
                      let all = 
                          match Alphabet.get_all alph with
                          | None -> (-1)
                          | Some x -> x
                      in
                      Array.init size 
                          (fun x -> 
                               Array.init size 
                                   (fun y -> 
                                        let x = min (x + 1) tcm_size in 
                                        let y = min (y + 1) tcm_size in
                                        let x =
                                            if x = all then x + 1
                                            else x
                                        and y =
                                            if y = all then y + 1
                                            else y
                                        in
                                        Cost_matrix.Two_D.cost x y tcm))
                | Alphabet.Extended_Bit_Flags -> 
                        failwith "Impliedalignment.make_tcm"
            in
            let enc = 
                let alph = Alphabet.to_sequential alph in
                let res = Parser.OldHennig.Encoding.default () in
                let res = Parser.OldHennig.Encoding.set_min res 0 in
                let res = Parser.OldHennig.Encoding.set_max res (size - 1) in
                let set = 
                    let rec add_consecutive_integers cur max acc = 
                        if cur = max then acc
                        else 
                            add_consecutive_integers (cur + 1) max 
                            (All_sets.Integers.add cur acc)
                    in
                    add_consecutive_integers 0 size All_sets.Integers.empty
                in
                let res = Parser.OldHennig.Encoding.set_set res set in
                alph, Parser.OldHennig.Encoding.set_sankoff res (make_tcm ())
            in
            let rec generate_all acc size = 
                if size < 0 then acc
                else generate_all (size :: acc) (size - 1)
            in
            let convert_to_list x =                
                match Alphabet.kind alph with
                | Alphabet.Simple_Bit_Flags ->
                        let rec match_bit v pos mask acc = 
                            if pos = 6 then acc
                            else if 0 <> (v land mask) then
                                let acc =
                                    if is_metric || mask <> gap then
                                        ((pos - 1) :: acc)
                                    else acc
                                in
                                match_bit v (pos + 1) (mask lsl 1) acc
                            else match_bit v (pos + 1) (mask lsl 1) acc
                        in
                        match_bit x 1 1 []
                | Alphabet.Sequential -> 
                        [x - 1]
                | Alphabet.Extended_Bit_Flags -> 
                        failwith "Impliedalignment.convert_to_list"
            in
            let all = 
                (* The size minus 1 for the codes (starting in 0), and another
                * one for the gap which we code separated *)
                generate_all [] (size - 2) in
            let gap_holder = 
                (* We always use all in the gap because we code it separated for
                * Sankoff characters 
                if is_metric then [gap_code] else all 
                *) all
            in
            let table = Hashtbl.create 67 in
            let find_item it =
                if Hashtbl.mem table it then
                    Hashtbl.find table it
                else begin
                    let r = Parser.Sankoff_Character (it, false) in
                    Hashtbl.add table it r;
                    r
                end
            in
            let to_parser is_missing states acc = 
                match is_missing, states with
                | `Missing, _ -> 
                        (find_item all) :: acc
                | `Exists, 0 -> 
                        (find_item gap_holder) :: acc
                | `Exists, x -> 
                        (find_item (convert_to_list x)) :: acc
            and to_encoding _ acc = 
                enc :: acc 
            in
            get_case, to_parser, to_encoding


module type S = sig
    type a 
    type b
    type tree = (a, b) Ptree.p_tree

            
    (** [of_tree t] generates the implied alignment of all the sequences in the tree
    * [t]. *)
    val of_tree : ((int -> int) * tree) -> Methods.implied_alignment


    val create : (tree -> int list -> tree) ->
        int list -> Data.d ->
        tree -> Methods.implied_alignment list

    val to_static_homologies : string -> bool ->
        (tree -> int list -> tree) ->
            bool  -> Methods.characters -> Data.d -> tree -> Data.d

end
module Make (Node : NodeSig.S) (Edge : Edge.EdgeSig with type n = Node.n) = struct
    type a = Node.n
    type b =  Edge.e
    type tree = (a, b) Ptree.p_tree

    (** return (taxon_id, character_ls) list (of taxa) * (final ias for each
        character set) list (of characters) *)
    let of_tree_handle all_minus_gap cg handle ptree =
        let calculate_median = false
        (* We turn this off because `Normal and `Exact would produce the wrong
        * implied alignment when using affine gap costs 
            match !Methods.cost with
            | `Normal | `Exact -> true
            | `Iterative -> false *)
        in
        let get_dynamic_data = 
            if calculate_median then Node.get_dynamic_preliminary 
            else Node.get_dynamic_adjusted 
        in
        let vertices = 
            try
                let root = Ptree.get_component_root handle ptree in
                match root.Ptree.root_median with
                | Some (_, v) -> Some (((Node.num_otus None v) * 2) - 1)
                | None -> None
            with
            | Not_found -> None
        in
        let st = 
            Status.create "Implied Alignments" vertices "vertices calculated"
        in
        let convert_data taxon_id data =
            let data = 
                List.map 
                (fun dyn ->
                    let sequences = DynamicCS.leaf_sequences dyn in
                    let state = DynamicCS.state dyn in 
                    let new_sequences = 
                        Codes.fold 
                        (fun code seq_arr acc ->
                            let ias_arr = 
                                Array.map 
                                (fun seq -> create_ias state seq taxon_id cg) 
                                seq_arr
                            in 
                            Codes.add code ias_arr acc
                        ) sequences Codes.empty  
                    in                                   
                    {sequences = new_sequences;   
                     c2 = DynamicCS.c2 dyn;   
                     chrom_pam = DynamicCS.chrom_pam dyn;  
                     state = DynamicCS.state dyn; 
                     alpha = DynamicCS.alpha dyn;
                     code = DynamicCS.code dyn; 
                     cannonic_code = taxon_id;
                     children = `Single taxon_id}
                ) 
                    data  
            in
            AssList.singleton (taxon_id, data), data
        in
        let convert_node parent ptree _ id _ =
            let data = Ptree.get_node_data id ptree in
            let taxon_id = Node.taxon_code data in
            let data = 
                let par =                     
                    match parent with
                    | Some _ -> parent
                    | None -> Some (Ptree.get_parent taxon_id ptree)
                in
                if Tree.is_leaf id ptree.Ptree.tree then
                    (* In a leaf we have to do something more complex, if we are
                    * dealing with simplified alphabets, not bitsets, we must
                    * pick the dynamic adjusted *)
                    let pre = Node.get_dynamic_preliminary par data
                    and adj = Node.get_dynamic_adjusted par data in
                    List.map2 (fun pre adj ->
                        if 0 = Cost_matrix.Two_D.combine (DynamicCS.c2 pre)
                        then
                            adj
                        else pre) pre adj
                else get_dynamic_data par data
            in
            let data = convert_data taxon_id data in
            let did = Status.get_achieved st in
            Status.full_report ~adv:(did + 1) st;
            data
        in
        let join_2_nodes _ _ (ac, a) (bc, b) =
            let t_ancestor x y =
                let ancestor_f =
                    assert (x.state = y.state);
                    match x.state with
                    | `Seq -> ancestor_sequence false
                    | `Chromosome -> ancestor_chrom true
                    | `Annotated -> ancestor_annchrom true
                    | `Breakinv -> ancestor_breakinv true
                    | `Genome -> ancestor_genome true
                in
                Codes.fold 
                (fun u v acc ->
                    let homs = Codes.find u y.sequences in

                    let ancestor = 
                        ancestor_f calculate_median all_minus_gap
                        x.cannonic_code y.cannonic_code x.children y.children 
                        v homs x.c2 x.alpha x.chrom_pam;
                    in
                    Codes.add u ancestor acc) x.sequences Codes.empty
            in
            let rec ancestor_builder hx hy =
                { hx with 
                cannonic_code = min hx.cannonic_code hy.cannonic_code;
                children = `Set [hx.children; hy.children];
                sequences = t_ancestor hx hy } 
            in
            AssList.union ac bc, List.map2 ancestor_builder a b
        in 
        match Tree.get_node handle ptree.Ptree.tree with
        | Tree.Single self -> 
                let a, b = convert_node None ptree () self ([], []) in
                Status.finished st;
                AssList.elements a, b
        | _ ->
              let self, other, root  =      
                  let root = Ptree.get_component_root handle ptree in  
                  match root.Ptree.root_median with
                  | Some ((`Edge (a, b)), the_root) -> a, b, the_root
                  | _ -> failwith "no root?"
              in
              let x, y =
                  if calculate_median then begin
                        let a, b = 
                            Ptree.post_order_node_with_edge_visit 
                            (convert_node None ptree) join_2_nodes 
                            (Tree.Edge (self, other)) ptree
                            (AssList.empty, [])
                        in 
                        let  a' = 
                            let new_ptree = 
                                let self_data = Ptree.get_node_data self ptree 
                                and other_data = Ptree.get_node_data other ptree
                                in
                                let single = 
                                    Node.to_single (Some root) (Some self) 
                                    other_data (Some other) self_data 
                                in
                                Ptree.add_node_data self single ptree
                            in
                            convert_node (Some other) new_ptree () self ([], []) 
                        in
                        let a = join_2_nodes () () a a' in
                        let tmp = join_2_nodes () () a b in
                        tmp
                  end 
                  else 
                        (* We will do a different tree traversal, converting each
                        interior vertex as we move down the tree starting with the
                        actual root *)
                        let rec tree_traverser parent_ia parent_code my_code =
                            let my_ia =
                                let my_ia = 
                                    convert_node parent_code ptree () my_code
                                    ([], [])
                                in
                                match Ptree.get_node my_code ptree with
                                | (Tree.Interior _) as n -> 
                                        let parent_code = 
                                            match parent_code with
                                            | Some parent_code -> parent_code
                                            | None -> assert false
                                        in
                                        let a, b = Tree.other_two_nbrs parent_code n in
                                        let my_ia = 
                                            tree_traverser my_ia (Some my_code) a 
                                        in
                                        tree_traverser my_ia (Some my_code) b 
                                | _ -> my_ia
                            in
                            join_2_nodes () () parent_ia my_ia
                        in
                        match Ptree.get_node self ptree with
                        | Tree.Interior (_, a, b, c) ->
                                let myd = Ptree.get_node_data  self ptree in
                                let my_ia = 
                                    convert_data self
                                    (get_dynamic_data (Some other) myd)
                                in
                                let self = Some self in
                                let my_ia = tree_traverser my_ia self a in
                                let my_ia = tree_traverser my_ia self b in
                                tree_traverser my_ia self c
                        | Tree.Leaf (_, b) ->
                                let myd = Ptree.get_node_data self ptree in
                                let my_ia = 
                                    convert_data self (get_dynamic_data (Some b) myd)
                                in
                                tree_traverser my_ia (Some self) b
                        | Tree.Single _ ->
                                let data = Ptree.get_node_data self ptree in
                                convert_data (Node.taxon_code data) 
                                (get_dynamic_data None data)
              in
              let _ = Status.finished st in
              let cleanedup = 
                  List.fold_left (fun acc ((code, y) as assoc) ->
                      if Tree.is_leaf code ptree.Ptree.tree then
                          assoc :: acc
                      else acc) []
                      (AssList.elements x)
              in
              cleanedup, y

    (** t is a final ias for a character set   
        =>
        (number columns, (com_hom_code -> column) map, (code -> com_hom_code) map) array map (a character set) *)
    let invert_codes t = 
        Codes.fold 
            (fun char_code iat_arr acc ->
                 let res_arr = 
                     Array.map 
                     (fun iat ->
                     (* number columns, (code -> column) map, (hom_code -> code) map *)
                     List.fold_left   
                     (fun (c, remap, acc) code ->  
                         let hom = Hashtbl.find iat.homologous code in 
                         let acc =   
                             Sexpr.fold_left 
                             (fun acc hom_code -> 
                                 Codes.add hom_code code acc) 
                             acc hom 
                         in  
                         c + 1, Codes.add code c remap, acc) 
                     (0, Codes.empty, Codes.empty) iat.order) 
                     iat_arr
                 in 
                 Codes.add char_code res_arr acc  
            ) t.sequences Codes.empty  
 (**  (number columns, (com_hom_code -> column) map, (code -> com_hom_code) map
  *  return arrays of aligned bases and aligned positions
  *) 
    let convert_a_taxon fi_ias_arr tax_ias_arr =
        let num_ias = Array.length fi_ias_arr in 
        match num_ias with 
        | 1 ->
              let len, remap, recode = fi_ias_arr.(0) in 
              let ias = tax_ias_arr.(0) in 
              let column code = 
                  let fin_code =  
                      try   
                          Codes.find code recode  
                      with | Not_found ->  
                          failwith ("Not_found hom_code ->  find_code " ^ (string_of_int code))  
                  in 
                  try 
                      Codes.find fin_code remap 
                  with Not_found ->
                      failwith ("Not_found in fin_code -> column" ^ (string_of_int fin_code))  
              in
              let results = Array.make (len + 1) 0 in
              let pos_results = Array.make (len + 1) (-1) in
              let add_result ias =
                  Hashtbl.iter (fun pos code -> 
                      let base = Sequence.get ias.seq pos in
                      let col = 
                          let col = column code in
                          len - col 
                      in
                      results.(col) <- base;
                      pos_results.(col) <- pos) 
                  ias.codes
              in
              add_result ias;
              [|results|], [|pos_results|]
        | _ -> (** that is for annotated or multiple chromosomes *)
              let pos_mat = Array.map 
                  (fun (len, _, _) -> Array.make (len + 1) (-1)) fi_ias_arr 
              in 
              let base_mat = Array.map 
                  (fun (len, _, _) -> Array.make (len + 1) 0) fi_ias_arr 
              in 
              let adder fi_ias pos_arr base_arr = 
                  let len, remap, recode = fi_ias in 
                  Array.iteri 
                      (fun chrom_idx tax_ias -> 
                           let seq = Sequence.to_array tax_ias.seq in                                
                           Array.iteri 
                               (fun (pos : int) (base : int) ->
                                    let code = Hashtbl.find tax_ias.codes pos in
                                    let com_hom_code = 
                                        try 
                                            Codes.find code recode 
                                        with  Not_found -> (-1)
                                    in 
                                    if (com_hom_code != -1) then begin
                                        let col = Codes.find com_hom_code remap in
                                        base_arr.(len - col) <- base;
                                        pos_arr.(len-col) <- Utl.max_seq_len * chrom_idx + pos;
                                    end 
                               ) seq;                               
                      ) tax_ias_arr;
              in 
              Array.iteri (fun idx fi_ias -> adder fi_ias pos_mat.(idx) base_mat.(idx)) fi_ias_arr;
              base_mat, pos_mat
   (** (taxon_id * list of alignments for each character 
       in the character set) list (of characters) ) list (of taxa) *)
    let of_tree (all_but_gap, tree) = 
        let cg = code_generator () in
    (** return ( (taxon_id, character_ls) list (of taxa) * (final ias for each character)
        list (of characters) ) list (of handles) *)
        let res = 
            Handles.fold 
            (fun (x : int) acc -> 
                (of_tree_handle all_but_gap cg x tree) :: acc) 
            (Tree.get_handles tree.Ptree.tree) []
        in
        let indel_blocks = 
            List.map (fun (_, x) -> 
                (List.map (fun x -> 
                    Codes.fold (fun _ arr acc ->
                        Array.fold_left 
                            (fun acc y ->  y.indels :: acc) acc arr) x.sequences []) x)) res
        in

        let dum_chars = 
            List.map (fun (_, x) -> 
                (List.map (fun x -> 
                    Codes.fold (fun _ arr acc ->
                        Array.fold_left 
                            (fun acc y ->  y.dum_chars :: acc) acc arr) x.sequences []) x)) res
        in

   (** ((taxon_id * (aligned_code arrays for each character
      set) list (of characters) ) list (of taxa) ) of list (of handles)*)
        let ali = List.map 
            (fun ((a:pairs list), (b:t list)) -> 
                (** (char_code -> (int * int Codes.t * 
                int Codes.t)) Codes.t list 
                (number columns, (com_hom_code -> column) map, (code -> com_hom_code) map) 
                map (a character set) list (of character sets)  *) 
                let inv_codes_b = List.map invert_codes b in
                      (** (taxon_id * (aligned_code arrays for each character
                          set) list (of characters) ) list (of taxa) *)
                 let char_states = ref (Codes.empty) in 
                 let new_a = List.map  
                      (fun ((tc : int), (tit : t list)) ->
                           (** a= (number_colum * (com_hom_code -> colum) map *
                               (code -> com_hom_code) map ) 
                                array map (a chracter set) list (of character sets) 
                               b = t list *)
                           let rec builder a b = 
                               match a, b with 
                               | ha :: ta, hb :: tb -> 
                                     let result =  Codes.fold 
                                         (fun char_code fi_hom (acc_bases, acc_pos) ->  
                                             (** fi_hom = 
                                                 * (number_colum * (com_hom_code -> colum) map 
                                                 * * (code -> com_hom_code) map ) array *)
                                              let state = hb.state in
                                              char_states := 
                                                  (Codes.add char_code state !char_states);
                                              let seq_arr = Codes.find char_code hb.sequences in 
                                              let bases, pos = convert_a_taxon fi_hom seq_arr in 
                                              (Codes.add char_code bases acc_bases),
                                              (Codes.add char_code pos acc_pos)) 
                                         ha (Codes.empty, Codes.empty)
                                     in 
                                     (result) :: (builder ta tb) 
                               | [], [] -> [] 
                               | _, _ -> failwith "Unexpected" 
                           in 
                           tc, builder inv_codes_b tit 
                      ) a 
                 in 
                 let get_ali_pos looking_char_code = 
                     List.fold_left 
                         (fun alied_seq_ls (tx_code, ali) ->
                              List.fold_left 
                                  (fun alied_seq_ls (_, alied_pos) ->
                                       let alied_seq = Codes.find
                                           looking_char_code alied_pos 
                                       in        
                                       alied_seq::alied_seq_ls
                                  ) alied_seq_ls ali
                         ) [] new_a
                 in 
                 let find_break_map () = 
                     Codes.fold 
                         (fun char_code state break_pos_map ->
                              let ali_pos_ls = get_ali_pos char_code in
                              
                              let ali_pos_ls = List.map 
                                  (fun ali_pos_mat -> 
                                       let pos_arr = Array.concat (Array.to_list ali_pos_mat) in
                                       pos_arr
                                  ) ali_pos_ls
                              in 
                              let alied_pos_mat = Array.of_list ali_pos_ls in 
                              let num_taxa = Array.length alied_pos_mat in 
                              let num_col = Array.length (alied_pos_mat.(0)) in
                              let break_ls = ref [] in                               
                              let sta = ref 0 in 
                              let cur_pos_arr = Array.init num_taxa (fun ti -> alied_pos_mat.(ti).(0)) in 
                              for col = 1 to num_col - 1 do
                                  let continue = ref true in 
                                  for t = 0 to num_taxa - 1 do 
                                      let pre_pos = cur_pos_arr.(t) in 
                                      let pos = alied_pos_mat.(t).(col) in
                                      (if (pre_pos != -1) && (pos != -1) 
                                           && (abs (pos - pre_pos) > 1) then continue := false);
                                  done;  
                                  if !continue then begin 
                                      for t = 0 to num_taxa - 1 do
                                          if alied_pos_mat.(t).(col) != -1 then 
                                              cur_pos_arr.(t) <- alied_pos_mat.(t).(col);
                                      done;
                                  end  else begin
                                      break_ls:= List.append !break_ls [(!sta, (col - 1))] ;
                                      sta := col;
                                      for t = 0 to num_taxa - 1 do
                                          cur_pos_arr.(t) <- alied_pos_mat.(t).(col)
                                      done;
                                  end;
                              done;    
                              break_ls:= List.append  !break_ls [(!sta, (num_col - 1))]; 
                              Codes.add char_code !break_ls break_pos_map
                         ) !char_states (Codes.empty)  
                 in
                 List.map 
                 (fun (taxa_code, new_a_taxa) ->
                     let new_a_taxa = List.map 
                     (fun (base_map, _) -> Codes.mapi 
                     (fun char_code alied_seq ->
                         let char_state = Codes.find char_code !char_states in 
                         match char_state with 
                         | `Seq | `Annotated | `Breakinv -> alied_seq
                         | `Chromosome | `Genome ->
                                 let break_map = find_break_map () in
                                 let alied_seq = Array.concat (Array.to_list alied_seq) in
                                 let break_ls = Codes.find char_code break_map in 
                                 let seg_ls = Utl.break_array alied_seq break_ls in 
                                 let seg_arr = Array.of_list seg_ls in 
                                 for idx = 1 to (Array.length seg_arr) - 1 do
                                     seg_arr.(idx) <- 
                                         Array.of_list( 0 :: 
                                             (Array.to_list seg_arr.(idx)))
                                 done;
                                 seg_arr) 
                     base_map) new_a_taxa in 
                     (taxa_code, new_a_taxa)) 
                 new_a) 
            res  
        in 
        List.combine (List.combine ali indel_blocks) dum_chars
(** End of of_tree function *)

    let post_process_affine_gap_cost f data (enc, taxa)
    all_blocks=
        let all_blocks = `Set all_blocks in
        let process_indel (enc, taxa) (loc, string, length, clas, taxa_list) =
            let present_absent_alph = 
                Alphabet.list_to_a 
                [(string, 1, None); ("-", 2, None)] 
                "-" None Alphabet.Sequential
            in
            let in_taxa, not_in_taxa =
                match clas with 
                | `Insertion -> 
                        Parser.Unordered_Character (1, false), 
                        Parser.Unordered_Character (2, false)
                | `Deletion -> 
                        Parser.Unordered_Character (2, false),
                        Parser.Unordered_Character (1, false)
            in
            let taxa_list : All_sets.Integers.t = 
                Sexpr.fold_left 
                (fun acc x -> All_sets.Integers.add x acc) 
                All_sets.Integers.empty taxa_list
            in
            let newenc = 
                Parser.OldHennig.Encoding.gap_encoding (f string) 
            in
            ((present_absent_alph, newenc) :: enc),
            List.map (fun (characters, taxon) ->
                let code = Data.taxon_code taxon data in
                (if All_sets.Integers.mem code taxa_list then
                    in_taxa :: characters
                else not_in_taxa :: characters), taxon) taxa
        in
        let acc = Array.to_list enc, 
        List.map (fun (y, x) -> Array.to_list y, x) taxa 
        in
        let enc, taxa = Sexpr.fold_left process_indel acc all_blocks in
        Array.of_list enc, List.map (fun (y, x) -> Array.of_list y, x) taxa



    let post_process_dum_cost data (enc, taxa) dum_chars =
        let dum_chars = `Set dum_chars in
        let process_rearr (enc, taxa) (dum_cost, taxa_list) =
            let present_absent_alph = 
                Alphabet.list_to_a 
                [("R", 1, None); ("r", 2, None)] 
                "r" None Alphabet.Sequential
            in
            let in_taxa, not_in_taxa = 
                        Parser.Unordered_Character (1, false), 
                        Parser.Unordered_Character (2, false)
            in 
            let taxa_list : All_sets.Integers.t = 
                Sexpr.fold_left 
                (fun acc x -> All_sets.Integers.add x acc) 
                All_sets.Integers.empty taxa_list
            in
            let newenc = 
                Parser.OldHennig.Encoding.rearr_encoding (dum_cost)
            in
            ((present_absent_alph, newenc) :: enc),
            List.map (fun (characters, taxon) ->
                let code = Data.taxon_code taxon data in
                (if All_sets.Integers.mem code taxa_list then
                    in_taxa :: characters
                else not_in_taxa :: characters), taxon) taxa
        in
        let acc = Array.to_list enc, 
        List.map (fun (y, x) -> Array.to_list y, x) taxa 
        in
        let enc, taxa = Sexpr.fold_left process_rearr acc dum_chars in
        Array.of_list enc, List.map (fun (y, x) -> Array.of_list y, x) taxa



    let ia_to_parser_compatible data imtx =
        match imtx with
        | [(all_taxa, all_blocks), dum_chars] ->
              let process_each = fun (acc, enc, clas) (taxcode, sequence) ->
                    let preprocess_sequence alph x =
                        let len = Array.length x in
                        let rec check it =
                            if it = len then true
                            else if x.(it) = 0 then check (it + 1)
                            else false
                        in
                        if check 0 then begin
                            for i = len - 1 downto 0 do
                                x.(i) <- 0;
                            done;
                            `Missing
                        end else `Exists
                    in
                    (* Fold over every sequence, and return a list containing all of
                    * them in a tuple with their code and the sequence itself in a
                    * preprocessed way (if the sequence is missing data, then all
                    * the states are on), and a map of the sequence code to the
                    * function that will convert an observed state into it's
                    * appropriate Parser.t *)
                    let sequence, transform_functions = 
                        List.fold_left (fun acc s ->
                            All_sets.IntegerMap.fold
                                (fun c (s_arr : int array array) (acc, funs) -> 
                                    let alph = Data.get_alphabet data c in
                                    let funs =
                                        if All_sets.IntegerMap.mem c funs then funs
                                        else begin
                                            let tcm = Data.get_sequence_tcm c data in
                                            let mc, f, g = analyze_tcm tcm alph in
                                            All_sets.IntegerMap.add c (mc, f, g) funs
                                        end 
                                    in
                                    let mis_sta = preprocess_sequence alph 
                                                 (Array.concat (Array.to_list s_arr)) 
                                    in
                                    let s_arr' = 
                                        Array.map (fun s ->                                            
                                            c, (mis_sta, s)) 
                                        s_arr 
                                    in 
                                    List.append (Array.to_list s_arr') acc, funs)  
                                s acc)
                        ([], All_sets.IntegerMap.empty) sequence 
                    in
                    let clas, res, encf = 
                        List.fold_left 
                        (fun (_, acc, acc2) (code, (is_missing, s)) -> 
                            let clas, to_parser, to_encoding = 
                                All_sets.IntegerMap.find code transform_functions 
                            in
                            let is_missing = 
                                match is_missing with 
                                | `Exists -> `Exists
                                | `Missing -> begin
                                      let state = Data.get_character_state data code in 
                                      match state with 
                                      | `Chromosome | `Annotated ->
                                            let seq_len = Array.length s in 
                                            let seq_op, seq_ex = 
                                                match clas with
                                                | AllOne gapc -> 0, gapc
                                                | AllOneGapSame (_, gapc) ->  0, gapc
                                                | AffinePartition (_, ex, op) -> op, ex
                                                | AllSankoff _ -> 
                                                        (* This is an error, but
                                                        * I am confused with the
                                                        * handling of seq_ex
                                                        * below. Somehow it
                                                        * assumes that
                                                        * ex = seq_ex? And if
                                                        * sankoff is used this
                                                        * should raise an error
                                                        * ... *)
                                                        0, 0
                                            in 
                                            let pam = Data.get_pam data code in 
                                            let op, ex = 
                                                match pam.Data.locus_indel_cost with
                                                | Some (op, ex) -> op, ex 
                                                | None -> ChromPam.locus_indel_cost_default 
                                            in
                                            let locus_indel_cost = op + ex * seq_len / 100 in 
                                            let num_gaps = (locus_indel_cost - seq_op)/seq_ex in
                                            let alph = Data.get_alphabet data code in 
                                            let all = Utl.deref (Alphabet.get_all alph) in 
                                            for p = 0 to seq_len - num_gaps - 1 do
                                                s.(p) <- all;
                                            done;
                                            `Exists
                                      | _ -> `Missing
                                  end 
                            in  
                            clas,
                            (Array.fold_right (to_parser is_missing) s acc), 
                            (Array.fold_right to_encoding s acc2))
                        (AllSankoff None, [], []) sequence 
                    in 
                    let name = 
                        try Data.code_taxon taxcode data with
                        | Not_found -> (string_of_int taxcode) 
                    in
                    match enc with
                    | Some _ -> (res, name) :: acc, enc, clas
                    | None ->
                            (res, name) :: acc,
                                (let rec apply_map l1 l2 =
                                    match l1, l2 with
                                    | f1 :: t1, it2 :: t2 ->
                                            (f1 it2) :: (apply_map t1 t2)
                                    | [], [] -> []
                                    | _, _ -> failwith "Not matching numbers?"
                                in
                                Some encf), clas
                in
                (match List.fold_left process_each ([], None, AllSankoff None) all_taxa with
                | r, Some enc, clas -> 
                        let arr = 
                            Array.of_list enc, (List.map (fun (x, y) ->
                                Array.of_list x, y) r)
                        in
                        let (a, b) = 
                            match clas with
                            | AffinePartition (subs, gapcost, gapopening) ->
                                    (* We have to postprocess and check by
                                    * groups of three whether or not we have a
                                    * gap opening indeed *)
                                    let f string =
                                        (gapopening + (gapcost * (String.length string))) 
                                    in
                                    List.fold_left
                                    (post_process_affine_gap_cost f data) arr all_blocks
                            | AllSankoff (Some f) -> 
                                    List.fold_left 
                                    (post_process_affine_gap_cost f data) arr
                                    all_blocks
                            | _ -> arr
                        in

                        let (a,b) = List.fold_left
                                    (post_process_dum_cost data) (a,b) dum_chars
                        in

                        a, b, []
                | [], _, _ -> [||], [], []
                | _, None, _ -> failwith "How is this possible?")
        | _ -> failwith "ImpliedAlignment.ia_to_parser_compatible"

    let update_ia_encodings (encs, species, trees) =
        let add_states int acc =
            match int with
            | Parser.Unordered_Character (int, _) ->
                    let a = 1 land int
                    and b = 2 land int 
                    and c = 4 land int 
                    and d = 8 land int 
                    and e = 16 land int in
                    let addit acc x =
                        if x <> 0 then All_sets.Integers.add x acc
                        else acc
                    in
                    addit (addit (addit (addit (addit acc a) b) c) d) e
            | _ -> failwith "Unordered here?"
        in
        let arr = Array.of_list species in
        let arr = Array.map (fun (a, _) -> a) arr in
        let updater pos (alph, enc) = 
            if Parser.OldHennig.Encoding.is_sankoff enc ||
                Parser.OldHennig.Encoding.is_ordered enc then alph, enc
            else 
                let ns = Array.fold_left (fun acc taxon ->
                    add_states taxon.(pos) acc) All_sets.Integers.empty 
                    arr 
                in
                alph, Parser.OldHennig.Encoding.set_set enc ns
        in
        let arr = Array.mapi updater encs in
        arr, species, trees

    let to_static_character remove_non_informative character iamtx data = 
        let st = 
            Status.create "Static Approximation" None 
            "Converting implied alignments to static characters"
        in
        let res = (ia_to_parser_compatible data iamtx) in
        let (a, b, c) = 
            if remove_non_informative then update_ia_encodings res 
            else res
        in
        let alphabets = Array.map fst a
        and encodings = Array.map snd a in
        let res = 
            Parser.SC.of_old_parser character (Some alphabets) (encodings, b, c)
        in
        Status.finished st;
        character, res

   (** (sequence code list), ( (taxon_id * (aligned_code arrays for each character
       set) list (of characters) ) list (of taxa) ) of list (of trees) *)
    let aux_create_implied_alignment filter_fn codes data tree = 
        let operate_on_tree tree =
            let filtered_trees = 
                List.map 
                    (fun x -> 
                         let alph = Data.get_alphabet data x in
                         let gap = Alphabet.get_gap alph in
                         let tcm = Data.get_sequence_tcm x data in
                         let all = 
                             if 1 = Cost_matrix.Two_D.combine tcm then
                                 match Alphabet.get_all alph with
                                 | Some all -> all
                                 | None -> assert false
                             else (-1) (* we won't use it anyway *)
                         in 
                         (if 1 = Cost_matrix.Two_D.combine tcm then
                            fun x -> x land ((lnot gap) land all)
                         else fun x -> x), filter_fn tree [x]) codes
            in
            List.map of_tree filtered_trees
        in
        codes, operate_on_tree tree



   (** (sequence code list), ( (taxon_id * (aligned_code arrays for each character
       set) list (of characters) ) list (of taxa) ) of list (of trees) *)
    let create filter_fn codes data tree = 
        let codes = (* Check if the codes are sequence codes or not *) 
            List.filter (fun x -> 
                if (List.exists (fun y -> x = y) data.Data.dynamics) then true
                else begin
                    Status.user_message Status.Error
                    ("The character with code " ^ string_of_int x ^
                    " is not a sequence character. You have requested an "
                    ^ "implied alignment of such thing, I will ignore that "
                    ^ "character for the implied alignment.");
                    false
                end) codes
        in
        let _, ia = aux_create_implied_alignment filter_fn codes data tree in
        ia

    let get_char_codes (chars : Methods.characters)  data =
        let codes = 
            let codes = 
                match chars with
                | `Some (dont_complement, codes) ->
                        let codes = Data.get_chars_codes data (`Some codes) in
                        if dont_complement then `Some codes
                        else Data.complement_characters data (`Some codes)
                | `Names (dont_complement, names) ->
                        let codes = Data.get_chars_codes data (`Names names) in
                        if dont_complement then `Some codes
                        else Data.complement_characters data (`Some codes)
                | `Random _ | `Missing _ | `All | `AllStatic | `AllDynamic as x ->
                        `Some (Data.get_chars_codes data x)
            in
            let codes = 
                match codes with
                | `Some x -> x
                | _ -> failwith "Impossible?"
            in
            (* Ensure we are not removing anything that is not affected by this
            * transformation *)
            let dyn = 
                (Data.get_code_from_characters_restricted `Dynamic data 
                (`Some codes))
            in
            List.filter (fun x -> List.exists (fun y -> y = x) dyn) codes 
        in
        codes
    
    let to_static_homologies prefix ignore filter_fn remove_noninformative 
            (chars : Methods.characters)  data tree = 
        let codes = get_char_codes chars data in
        let names = List.map (fun x -> Data.code_character x data) codes in
        let all_to_add = 
            List.fold_left (fun acc code -> 
                let _, ia = 
                    aux_create_implied_alignment filter_fn [code] data tree 
                in
                assert (1 = List.length ia);
                let ia = List.hd ia in

                let name =
                    let code = 
                        incr data.Data.character_code_gen;
                        string_of_int !(data.Data.character_code_gen)
                    in
                    prefix ^ code
                in
                (to_static_character remove_noninformative name ia data) :: acc) 
            []
            codes
        in
        let d = Data.add_multiple_static_parsed_file data all_to_add in
        if ignore then Data.process_ignore_characters false d (`Names names)
        else d

end 
