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

let () = SadmanOutput.register "Diagnosis" "$Revision: 2871 $"

let debug = true

let sort_using_tree tree all_taxa =
    let set = 
        List.fold_left 
        (fun acc ((x, y) as z) -> All_sets.IntegerMap.add x z acc) 
        All_sets.IntegerMap.empty all_taxa
    in
    let res = 
        All_sets.Integers.fold (fun handle acc ->
            Tree.post_order_node_visit 
            (fun _ x acc -> 
                Tree.Continue,
                if All_sets.IntegerMap.mem x set then 
                    (All_sets.IntegerMap.find x set) :: acc
                else acc)
            handle tree acc)
        (Tree.get_handles tree) []
    in
    List.rev res

let all_missing_codes indels =
    let add_sexpr set sexpr = 
        Sexpr.fold_left (fun x y -> All_sets.Integers.add y x) set sexpr
    in
    List.fold_left
        (fun acc x -> List.fold_left
            (fun acc x -> Sexpr.fold_left
                (fun acc (_,_,_,x,cs) -> match x with 
                    | `Missing -> add_sexpr acc cs 
                    | `Insertion 
                    | `Deletion -> acc)
                acc x)
            acc x)
        (All_sets.Integers.empty)
        (indels)

let replace_with_missing missing a =
    let s = ref "" in 
    for i = (String.length a)-1 downto 0 do 
        s := missing ^ !s 
    done; 
    !s

(** Functions to output immediately all the taxa without memory consumption *)
module O = struct
    let header_f filename seqname () = 
        Status.user_message (Status.Output (filename, false, [])) 
        ("@[<v>@,@[New Tree@ for@ sequence@ " ^ StatusCommon.escape seqname ^ 
        "@]@,@[")

    let rec sequence_f ?(break=true) alphabet gap (taxon, result, pos, cnt) base =
        if cnt = 0 then (taxon, result, pos, 1) 
        else if break && 80 = pos then 
            sequence_f alphabet gap (taxon, result ^ "@,", 0, cnt) base
        else if base = 0 then 
            taxon, result ^ gap, pos + 1, cnt
        else 
            let item = 
                Alphabet.match_code base alphabet 
        in
        (taxon, result ^ item, pos + 1, cnt) 

    let taxonf name = (name, "", 0, 0)

    let taxon_closef filename (name, result, _, _) =
        let fo = Status.Output (filename, false, []) in
        Status.user_message fo ("@,@[<v>@,>" ^ StatusCommon.escape name ^
        "@,");
        Status.user_message fo result;
        Status.user_message fo "@]%!"

    let closef filename () = 
        Status.user_message (Status.Output (filename, false, [])) "@]@]"

end

(**************)

(** Functions to concatenate all of the sequences belonging to the same taxon 
* in one long line *)

module C = struct
    let break = false
    let header_f _ _ () = ()
    let sequence_f missings code a g data b = 
        if All_sets.Integers.mem code missings then
            let t,r,p,c = O.sequence_f ~break a g data b in
            (t, replace_with_missing (Alphabet.get_missing a) r, p, c)
        else
            O.sequence_f ~break a g data b

    let taxonf = O.taxonf
    let taxon_closef hash filename (name, result, _, _) =
        Hashtbl.add hash name result
    let closef _ () = ()
    let output_all_taxa include_header filename seqname hash =
        if include_header then
            O.header_f filename seqname ()
        else 
            Status.user_message (Status.Output (filename, false, [])) 
            ("@[<v>@,@[");
        let names = ref All_sets.Strings.empty in
        Hashtbl.iter (fun x b -> 
            if not (All_sets.Strings.mem x !names) then
                names := All_sets.Strings.add x !names
            else ()) hash;
        All_sets.Strings.iter (fun taxon ->
            let terms = List.rev (Hashtbl.find_all hash taxon) in
            let fo = Status.user_message (Status.Output (filename, false, [])) in
            fo ("@,@[<v>@,>" ^ StatusCommon.escape taxon ^ "@,");
            List.iter (fun x -> fo x; fo " ") terms;
            fo "@]%!") !names 
end

(**************)
let output_implied_alignment (tree, seqname) headerf taxonf taxon_closef seqf 
closef data to_process = 
    headerf ();
    (* This function expects only one element *)
    match to_process with
    | [(all_taxa, indels), _] ->
        let all_taxa = sort_using_tree tree all_taxa in
        let all_missing = all_missing_codes indels in
        let process_each acc (taxcode, sequence) =
            let name = 
                try Data.code_taxon taxcode data 
                with | Not_found -> (string_of_int taxcode)
            in
            match sequence with
            | hd_sequence :: tl ->
                let res = 
                    All_sets.IntegerMap.fold (fun c s acc -> Some (c, s)) hd_sequence None
                in
                begin match res with 
                    | Some (seqcode, sequence_arr) ->
                        let sequence = sequence_arr.(0) in 
                        let alphabet = Data.get_sequence_alphabet seqcode data in
                        let gapcode = Alphabet.get_gap alphabet in
                        let gap = Alphabet.match_code gapcode alphabet in
                        (* Check if the sequence is missing data *)
                        let preprocess_sequence x = match x with
                            | `DO x -> x
                            | `Last x
                            | `First x as m ->
                                let all = Utl.deref (Alphabet.get_all alphabet) in 
                                let len = Array.length x in
                                try match m with
                                    | `Last x ->
                                        for i = len - 1 downto 0 do
                                            if x.(i) <> 0 then
                                                raise Exit
                                            else x.(i) <- all;
                                        done;
                                        x
                                    | `First x ->
                                        for i = 0 to len - 1 do
                                            if x.(i) <> 0 then
                                                raise Exit
                                            else x.(i) <- all;
                                        done;
                                        x
                                with | Exit -> x
                        in
                        let sequence = preprocess_sequence sequence in
                        let accumulator = 
                            Array.fold_left (seqf all_missing taxcode alphabet gap) 
                                            (taxonf name) 
                                            (sequence)
                        in
                        let () = taxon_closef accumulator in
                        if Array.length sequence_arr = 1 then 
                            (taxcode, tl) :: acc
                        else begin
                            let hd_sequence = 
                                All_sets.IntegerMap.map
                                    (fun sequence_arr -> 
                                        Array.of_list (List.tl (Array.to_list sequence_arr)))
                                    hd_sequence
                            in 
                            (taxcode, hd_sequence::tl)::acc
                        end 
                    | None -> (taxcode, tl) :: acc
                end
            | [] -> acc
        in
        let rec process_all_sequences = function
            | (_, []) :: _ -> ()
            | x -> let acc = List.fold_left process_each [] x in
                   let acc = List.rev acc in
                   process_all_sequences acc
        in
        process_all_sequences all_taxa;
        closef ();
    | _ -> failwith "Diagnosis.output_implied_alignment 1"

module type S = sig
    type a 
    type b

    val diagnosis :
      Data.d -> (a, b) Ptree.p_tree -> Methods.diagnosis -> unit
end

module Make 
    (Node : NodeSig.S with type other_n = Node.Standard.n) (Edge : Edge.EdgeSig with type n = Node.n) 
    (TreeOps : 
        Ptree.Tree_Operations 
                        with type a = Node.n with type b = Edge.e) = struct

    type a = Node.n
    type b = Edge.e

    module IA = ImpliedAlignment.Make (Node) (Edge)
    module CT = CharTransform.Make (Node) (Edge) (TreeOps)
    module TO = TreeOps 

    let report_all_roots fo tree =
        let report_root ((Tree.Edge (a, b)), cost) =
            fo ("@[" ^ string_of_int a ^ "-" ^ string_of_int b ^ ": " ^
            string_of_float cost ^ "@]@\n");
        in
        fo "@[<v 2>@{<u>Tree@}@,@[<v>";
        List.iter report_root (TO.root_costs tree);
        fo "@]";
        fo "@]@,"

    let rec diagnosis (data : Data.d) (tree : (a, b) Ptree.p_tree) = 
        function
        | `AllRootsCost filename ->
                let fo = 
                    Status.user_message (Status.Output (filename, false, [])) 
                in
                fo "@[<v>@,@,@{<b>All Roots Cost@}@,";
                report_all_roots fo tree;
                fo "@]@\n%!"
        | `Implied_Alignment (filename, chars, include_header) ->
                let char_codes = 
                    Data.get_code_from_characters_restricted_comp
                    `AllDynamic data chars 
                in
                let char_codes = List.sort compare char_codes in
                let res = 
                    List.map (fun code -> 
                        let seqname = Data.code_character code data in
                        let ia = IA.create CT.filter_characters [code] data tree in
                        (tree.Ptree.tree, seqname) , ia) char_codes  
                in
                let extract_name name = 
                    match Str.split (Str.regexp ":") name with
                    | h :: _ -> h
                    | [] -> name
                in
                let hash = ref (Hashtbl.create 97) 
                and name = ref "" in
                List.iter (fun (x, y) ->
                    if include_header then
                        Status.user_message 
                        (Status.Output (filename, false, [])) 
                        "@[<v>@,@,@{<b>Implied Alignments@}@,"
                    else ();
                    let (_, seqname) = x in
                    if (extract_name seqname) = !name then ()
                    else begin
                        C.output_all_taxa include_header filename !name !hash;
                        hash := Hashtbl.create 97;
                        name := extract_name seqname;
                    end;
                    let header_f = 
                        if include_header then 
                            C.header_f filename seqname
                        else (fun () -> ())
                    in
                    List.iter 
                    (output_implied_alignment x 
                    header_f C.taxonf 
                    (C.taxon_closef !hash filename) (C.sequence_f) 
                    (C.closef filename) data) y;
                    Status.user_message 
                    (Status.Output (filename, false, [])) "@]%!")
                res;
                C.output_all_taxa include_header filename !name !hash


end
    
