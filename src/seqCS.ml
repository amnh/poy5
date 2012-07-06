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

(** A Sequence Character Set implementation *)
exception Illegal_Arguments
let () = SadmanOutput.register "SeqCS" "$Revision: 2871 $"

let debug = false

module Codes = All_sets.IntegerMap

type union_element = 
    | Single of Sequence.Unions.u
    | Array of Sequence.Unions.u array


type cost_tuple = {
    min : float;
    max : float;
}

type packed_algn = 
    | Packed of (int * BitSet.t * packed_algn)
    | Raw of Sequence.s


type heuristic = {
    seqs: int;
    c2 : Cost_matrix.Two_D.m;
    c3 : Cost_matrix.Three_D.m;
}

let make_default_heuristic ?c3 c2 =
    let c3 =
        match c3 with
        | None -> Cost_matrix.Three_D.of_two_dim c2 
        | Some x -> x
    in
    { seqs = 1; c2 = c2; c3 = c3 }

module ProtAff = struct

    type s = Sequence.s
    module Unions = Sequence.Unions
    let encoding = Sequence.encoding
    let compare = Sequence.compare
    let to_formater = Sequence.to_formater
    let create = Sequence.create
    let length = Sequence.length
    let get = Sequence.get
    let prepend = Sequence.prepend
    let is_empty = Sequence.is_empty
    let close_block_diagonal = ref [|[|0|]|]
    let extend_block_diagonal = ref [|[|0|]|]
    let extend_vertical = ref [|[|0|]|]
    let extend_horizontal = ref [|[|0|]|]
    let final_cost_matrix = ref [|[|0|]|]
    let direction_matrix = ref [|[|0|]|]

    let create_arrays si sj =
        let is = Array.length !close_block_diagonal in
        let js = Array.length !close_block_diagonal.(0) in
        let new_is = max is (Array.length si)
        and new_js = max js (Array.length sj) in
        if new_is > is || new_js > js then begin
                let make w = 
                    w := Array.make_matrix new_is new_js 0 
                in
                make close_block_diagonal;
                make extend_block_diagonal;
                make extend_vertical;
                make extend_horizontal;
                make final_cost_matrix;
                make direction_matrix;
        end

    let ( --> ) a b = b a

    let max_int = 10000

    let has_gap_extension ge a b =
        assert (ge > 0);
        if 0 <> (16 land b) then 0
        else ge

    let align = 1
    let align_to_vertical = 2
    let align_to_horizontal = 4
    let align_to_diagonal = 8
    let begin_block = 16
    let end_block = 32
    let begin_vertical = 64
    let end_vertical = 128
    let begin_horizontal = 256
    let end_horizontal = 512
    let do_align = 1024
    let do_vertical = 2048
    let do_horizontal = 4096
    let do_diagonal = 8192

    let a = 1
    let c = 2
    let g = 4
    let t = 8
    let p = 16


    let align_sequences gap_open gap_extension substitution si sj =
        let min (a : int) (b : int) = if a <= b then a else b in
        let close_block_diagonal = !close_block_diagonal
        and extend_block_diagonal = !extend_block_diagonal
        and extend_vertical = !extend_vertical
        and extend_horizontal = !extend_horizontal
        and final_cost_matrix = !final_cost_matrix
        and direction_matrix = !direction_matrix in
        let lor_with_direction_matrix i j v =
            direction_matrix.(i).(j) <- direction_matrix.(i).(j) lor v
        in
        let assign_minimum i j =
            let mask = ref do_horizontal in
            final_cost_matrix.(i).(j) <- extend_horizontal.(i).(j);
            if final_cost_matrix.(i).(j) >= extend_vertical.(i).(j) then
                if final_cost_matrix.(i).(j) > extend_vertical.(i).(j) then
                    begin
                    final_cost_matrix.(i).(j) <- extend_vertical.(i).(j);
                    mask := do_vertical;
                    end
                else begin
                    mask := !mask lor do_vertical;
                end;
            if final_cost_matrix.(i).(j) >= extend_block_diagonal.(i).(j) then
                if final_cost_matrix.(i).(j) > extend_block_diagonal.(i).(j) then
                    begin
                    final_cost_matrix.(i).(j) <- extend_block_diagonal.(i).(j);
                    mask := do_diagonal;
                    end
                else begin
                    mask := !mask lor do_diagonal;
                end;
            if final_cost_matrix.(i).(j) >= close_block_diagonal.(i).(j) then
                if final_cost_matrix.(i).(j) > close_block_diagonal.(i).(j) then
                    begin
                    final_cost_matrix.(i).(j) <- close_block_diagonal.(i).(j);
                    mask := do_align;
                    end
                else begin
                    mask := !mask lor do_align;
                end;
            lor_with_direction_matrix i j !mask
        in
        let has_gap_opening a b =
            assert (gap_open > 0);
            if 0 = (16 land a) && 0 <> (16 land b) then 0
            else gap_open
        in
        let fill_extend_vertical i j =
            let res = 
                let ext_cost = 
                    (extend_vertical.(i - 1).(j) + 
                    (has_gap_extension gap_extension (si.(i - 1)) 
                    (si.(i))))
                and open_cost =
                    (close_block_diagonal.(i - 1).(j) + 
                    (has_gap_opening (si.(i - 1)) (si.(i))) +
                    (has_gap_extension gap_extension (si.(i - 1)) 
                    (si.(i))))
                in
                let () =
                    assert (ext_cost >= extend_vertical.(i - 1).(j));
                    assert (open_cost >= close_block_diagonal.(i - 1).(j));
                in
                if ext_cost < open_cost then
                    let () = lor_with_direction_matrix i j begin_vertical in
                    ext_cost
                else 
                    let () = lor_with_direction_matrix i j end_vertical in
                    open_cost
            in
            extend_vertical.(i).(j) <- res;
        in
        let fill_extend_horizontal i j =
            let res =
                let ext_cost =
                    (extend_horizontal.(i).(j - 1) +
                    (has_gap_extension gap_extension 
                    (sj.(j - 1)) (sj.(j))))
                and open_cost = 
                    (close_block_diagonal.(i).(j - 1) + 
                    (has_gap_opening (sj.(j - 1)) 
                    (sj.(j))) +
                    (has_gap_extension gap_extension 
                    (sj.(j - 1)) (sj.(j))))
                in
                let () =
                    assert (ext_cost >= extend_horizontal.(i).(j - 1));
                    assert (open_cost >= close_block_diagonal.(i).(j - 1));
                in
                if ext_cost < open_cost then
                    let () = lor_with_direction_matrix i j begin_horizontal in
                    ext_cost
                else 
                    let () = lor_with_direction_matrix i j end_horizontal in
                    open_cost
            in
            extend_horizontal.(i).(j) <- res;
        in
        let has_extend_block_diagonal i j =
            let i = si.(i)
            and j = sj.(j) in
            if (0 <> (16 land i)) && (0 <> (16 land j)) then
                0
            else max_int
        in
        let no_gap_cost i j =
            let i = si.(i)
            and j = sj.(j) in
            let ai = (lnot 16) land i
            and aj = (lnot 16) land j in
            if 0 = (ai land aj) then substitution
            else 0
        in
        let fill_extend_block_diagonal i j =
            let diag = has_extend_block_diagonal i j in
            let res =
                let ext_cost = (extend_block_diagonal.(i - 1).(j - 1) + diag)
                and open_cost = (close_block_diagonal.(i - 1).(j - 1) + diag) in
                if ext_cost < open_cost then
                    let () = lor_with_direction_matrix i j begin_block in
                    ext_cost
                else 
                    let () = lor_with_direction_matrix i j end_block in
                    open_cost
            in
            extend_block_diagonal.(i).(j) <- res;
        in
        let fill_close_block_diagonal i j =
            let diag = no_gap_cost i j in
            let extra_gap_opening =
                let i = si.(i)
                and j = sj.(j)
                and ip = si.(i - 1)
                and jp = sj.(j - 1) in
                (has_gap_opening ip i) +
                (has_gap_opening jp j)
            in
                let algn = close_block_diagonal.(i - 1).(j - 1) + diag
                and from_vertical = 
                    extend_vertical.(i - 1).(j - 1) + diag + extra_gap_opening
                and from_horizontal =
                    extend_horizontal.(i - 1).(j - 1) + diag + extra_gap_opening
                and from_diagonal = 
                    extend_block_diagonal.(i - 1).(j - 1) + diag +
                    extra_gap_opening
                in
                let mask = ref align in
                close_block_diagonal.(i).(j) <- algn;
                if close_block_diagonal.(i).(j) >= from_vertical then
                    if close_block_diagonal.(i).(j) > from_vertical then
                        begin
                            close_block_diagonal.(i).(j) <- from_vertical;
                            mask := align_to_vertical;
                        end
                    else begin
                        mask := !mask lor align_to_vertical;
                    end;
                if close_block_diagonal.(i).(j) >= from_horizontal then
                    if close_block_diagonal.(i).(j) > from_horizontal then
                        begin
                            close_block_diagonal.(i).(j) <- from_horizontal;
                            mask := align_to_horizontal;
                        end
                    else begin
                        mask := !mask lor from_horizontal;
                    end;
                if close_block_diagonal.(i).(j) >= from_diagonal then
                    if close_block_diagonal.(i).(j) > from_diagonal then
                        begin
                            close_block_diagonal.(i).(j) <- from_diagonal;
                            mask := align_to_diagonal;
                        end
                    else begin
                        mask := !mask lor from_diagonal;
                    end;
                lor_with_direction_matrix i j !mask
                (*
                let () = if smallest = 0 then
                    Printf.printf "The position %d %d has 0 cost with algn \
                    %d from_vertical %d from_horizontal %d and from_diagonal \
                    %d\n%!" i j algn from_vertical from_horizontal from_diagonal
                in
                *)
        in
        let leni = (Array.length si) - 1
        and lenj = (Array.length sj) - 1  in
        assert (lenj >= leni);
        let start_pos = ref 1
        and end_pos = ref (min lenj (max 40 ((lenj - leni) + 8))) 
        and start_v = 40 in
        for i = 1 to leni do
            if i > start_v then incr start_pos;
            direction_matrix.(i).(!start_pos - 1) <- do_horizontal lor end_horizontal;
            extend_horizontal.(i).(!start_pos - 1) <- max_int;
                close_block_diagonal.(i).(!start_pos - 1) <- max_int;
            for j = !start_pos to !end_pos do
                direction_matrix.(i).(j) <- 0;
                fill_extend_horizontal i j;
                fill_extend_vertical i j;
                fill_extend_block_diagonal i j;
                fill_close_block_diagonal i j;
                assign_minimum i j;
            done;
            if !end_pos < lenj then begin
                direction_matrix.(i).(!end_pos + 1) <- do_vertical lor
                end_vertical;
                extend_vertical.(i).(!end_pos + 1) <- max_int;
                close_block_diagonal.(i).(!end_pos + 1) <- max_int;
                incr end_pos;
            end;
        done;
        ()

    let initialize_matrices go ge s si sj =
        let close_block_diagonal = !close_block_diagonal
        and extend_block_diagonal = !extend_block_diagonal
        and extend_vertical = !extend_vertical
        and extend_horizontal = !extend_horizontal
        and final_cost_matrix = !final_cost_matrix
        and direction_matrix = !direction_matrix in
        final_cost_matrix.(0).(0) <- 0;
        close_block_diagonal.(0).(0) <- 0;
        extend_block_diagonal.(0).(0) <- 0;
        extend_horizontal.(0).(0) <- go;
        extend_vertical.(0).(0) <- go;
        direction_matrix.(0).(0) <- (2 * do_diagonal) - 1;
        let leni = (Array.length si) - 1
        and lenj = (Array.length sj) - 1  in
        for j = 1 to lenj do
            let jc = sj.(j)
            and jp = sj.(j - 1) in
            let r =
                extend_horizontal.(0).(j - 1) + 
                (has_gap_extension ge jp jc)
            in
            extend_horizontal.(0).(j) <- r;
            close_block_diagonal.(0).(j) <- r;
            final_cost_matrix.(0).(j) <- r;
            extend_block_diagonal.(0).(j) <- max_int;
            extend_vertical.(0).(j) <- max_int;
            direction_matrix.(0).(j) <- do_horizontal lor end_horizontal;
        done;
        direction_matrix.(0).(1) <- do_horizontal lor end_horizontal;
        for i = 1 to leni do
            let ic = si.(i)
            and ip = si.(i - 1) in
            let r =
                extend_vertical.(i - 1).(0) +
                (has_gap_extension ge ip ic)
            in
            extend_horizontal.(i).(0) <- max_int;
            close_block_diagonal.(i).(0) <- r;
            final_cost_matrix.(i).(0) <- r;
            extend_block_diagonal.(i).(0) <- max_int;
            extend_vertical.(i).(0) <- r;
            direction_matrix.(i).(0) <- do_vertical lor end_vertical;
        done;
        ()

    type mode = Todo | Vertical | Horizontal | Diagonal | Align

    exception Wrong_Position of (int * int * string)

    let to_s = function
        | Todo -> "Todo"
        | Vertical -> "Vertical"
        | Horizontal -> "Horizontal"
        | Diagonal -> "Diagonal"
        | Align -> "Align"

    let backtrace mtx median resi resj i j si sj mode =
        (*
        Printf.printf "Back %d %d\n%!" i j;
        *)
        let i = ref i
        and j = ref j in
        let ic = ref (Sequence.get si !i)
        and jc = ref (Sequence.get sj !j) in
        let mode = ref mode in
        let has_flag flag = 0 <> (mtx.(!i).(!j) land flag) in
        while (!i <> 0) && (!j <> 0) do
            match !mode with
            | Todo -> 
                    (* do_* branches and and don't move *)
                    mode :=
                        if has_flag do_align then Align
                        else if has_flag do_vertical then Vertical
                        else if has_flag do_horizontal then Horizontal
                        else if has_flag do_diagonal then Diagonal
                        else raise (Wrong_Position (!i, !j, to_s !mode))
            | Vertical -> 
                    (* begin_vertical and end_vertical *)
                    if has_flag end_vertical then mode := Todo;
                    if 0 = !ic land p then Sequence.prepend median (!ic lor p);
                    Sequence.prepend resi !ic;
                    Sequence.prepend resj p;
                    decr i;
                    ic  := Sequence.get si !i;
            | Horizontal -> 
                    (* begin_horizontal and end_horizontal *)
                    if has_flag end_horizontal then mode := Todo ;
                    if 0 = !jc land p then Sequence.prepend median (!jc lor p);
                    Sequence.prepend resi p;
                    Sequence.prepend resj !jc;
                    decr j;
                    jc  := Sequence.get sj !j;
            | Diagonal ->
                    (* begin_block and end_block *)
                    if has_flag end_block then mode := Todo;
                    Sequence.prepend resi !ic;
                    Sequence.prepend resj !jc;
                    decr i;
                    decr j;
                    jc  := Sequence.get sj !j;
                    ic  := Sequence.get si !i;
            | Align ->
                    (* align_* *)
                    if has_flag align_to_vertical then mode := Vertical
                    else if has_flag align_to_horizontal then mode := Horizontal
                    else if has_flag align_to_diagonal then mode :=
                        Diagonal;
                    let prep = 
                        if ((!ic land !jc) land (lnot p)) <> 0 then (!ic land !jc)
                        else (!ic lor !jc)
                    in
                    Sequence.prepend median prep;
                    Sequence.prepend resi !ic;
                    Sequence.prepend resj !jc;
                    decr i;
                    decr j;
                    jc  := Sequence.get sj !j;
                    ic  := Sequence.get si !i;
        done;
        while !i <> 0 do
            (* begin_vertical and end_vertical *)
            if 0 = !ic land p then Sequence.prepend median (!ic lor p);
            Sequence.prepend resi !ic;
            Sequence.prepend resj p;
            decr i;
            ic := Sequence.get si !i;
        done;
        while (!j <> 0) do
            (* begin_horizontal and end_horizontal *)
            if has_flag end_horizontal then mode := Todo ;
            if 0 = !jc land p then Sequence.prepend median (!jc lor p);
            Sequence.prepend resi p;
            Sequence.prepend resj !jc;
            decr j;
            jc := Sequence.get sj !j;
        done;
        Sequence.prepend resi 16;
        Sequence.prepend resj 16;
        if 16 <> Sequence.get median 0 then Sequence.prepend median p

    let print_matrix w h ttl mtx =
        print_endline ttl;
        for i = 0 to h - 1 do
            for j = 0 to w - 1 do
                Printf.printf "%d " mtx.(i).(j)
            done;
            print_newline ();
        done

    let backtrace si sj =
        let w = (Sequence.length sj) - 1
        and h = (Sequence.length si) - 1 in
        let create () = Sequence.create (w + h + 2) in
        let a = create ()
        and b = create () 
        and c = create () in
        assert (w >= h);
        try 
            backtrace !direction_matrix a b c h w si sj Todo;
            a, b, c
        with
        | err ->
                print_matrix w h "Final Matrix" !final_cost_matrix;
                print_matrix w h "Close block" !close_block_diagonal;
                print_matrix w h "Extend horizontal" !extend_horizontal;
                print_matrix w h "Extend vertical" !extend_vertical;
                print_matrix w h "Extend diagonal" !extend_block_diagonal;
                print_endline (Sequence.to_string si Alphabet.nucleotides);
                print_endline (Sequence.to_string sj Alphabet.nucleotides);
                raise err

    let seq_to_array s =
        Array.init (Sequence.length s) (fun x -> Sequence.get s x)

    let do_alignment go ge s si sj =
        let si' = seq_to_array si
        and sj' = seq_to_array sj in
        create_arrays si' sj';
        initialize_matrices go ge s si' sj';
        align_sequences go ge s si' sj';
        backtrace si sj

    let print_list ttl lst = 
        print_endline ttl;
        List.iter (Printf.printf "%d ") lst;
        print_newline ()

        (*
    let () = 
        let si = 
            [|p; a; t; t; g lor a; p lor g lor a; p lor t; a; t; p lor t; c|]
        and sj = [|p; a; g; p lor t; p lor t; c|] in
        let close_block_diagonal, extend_block_diagonal, extend_vertical,
        extend_horizontal, final_cost_matrix, direction_matrix = 
            do_alignment 6 1 1 si sj
        in
        print_matrix "Final Matrix" final_cost_matrix;
        print_matrix "Close block" close_block_diagonal;
        print_matrix "Extend horizontal" extend_horizontal;
        print_matrix "Extend vertical" extend_vertical;
        print_matrix "Extend diagonal" extend_block_diagonal;
        let median, bt1, bt2 = backtrace direction_matrix si sj in
        print_list "Median" median;
        print_list "Sequence 1" bt1;
        print_list "Sequence 2" bt2
        *)
    module Align = struct
        let of_list x = Sequence.of_array (Array.of_list x)
        let readjust_3d _ _ _ _ _ _ = failwith "Not implemented"
        let align_2 a b cm =
            let lena = Sequence.length a
            and lenb = Sequence.length b in
            let swap = lena > lenb in
            let go = Cost_matrix.Two_D.gap_opening cm in
            (*
            Printf.printf "Aligning\n%s\n%s\n%!" (Sequence.to_string
            a Alphabet.nucleotides ) (Sequence.to_string b Alphabet.nucleotides );
            *)
            let lst, a, b = 
                let lst, a, b =
                    if swap then do_alignment go 1 go b a 
                    else do_alignment go 1 go a b 
                in
                if swap then lst, b, a else lst, a, b
            in
            (*
            Printf.printf "Results are\n%s\n%s\n%s\nwith cost%d\n%!" 
            (Sequence.to_string lst Alphabet.nucleotides ) (Sequence.to_string
            a Alphabet.nucleotides ) (Sequence.to_string b Alphabet.nucleotides )
            (!final_cost_matrix.(lena - 1).(lenb - 1));
            print_matrix "Final Matrix" !final_cost_matrix;
            print_matrix "Close block" !close_block_diagonal;
            print_matrix "Extend horizontal" !extend_horizontal;
            print_matrix "Extend vertical" !extend_vertical;
            print_matrix "Extend diagonal" !extend_block_diagonal;
            *)
            lst, a, b,
            (!final_cost_matrix.((min lena lenb) - 1).((max lena lenb) - 1))
         let union _ _ = failwith "Not implemented"
     let cost_2 si sj cm = 
            let lena = Sequence.length si
            and lenb = Sequence.length sj in
            let si = seq_to_array si
            and sj = seq_to_array sj in
            let si, sj = 
                if lena > lenb then sj, si
                else si, sj
            in
            let go = Cost_matrix.Two_D.gap_opening cm in
            create_arrays si sj;
            initialize_matrices go 1 go si sj;
            align_sequences go 1 go si sj;
            (!final_cost_matrix.(min (lena - 1) (lenb - 1)).(max (lenb - 1) (lena - 1)))

        let closest parent mine c2 _ =
            let _, s1', s2', cst = align_2 parent mine c2 in
            let remove_gaps s2' =
                (* We first define a function to eliminate gaps from the 
                * final selection *)
                let remove_gaps gap seq base = 
                    if base <> gap then 
                        let _ = prepend seq base in
                        seq
                    else seq
                in
                let res = 
                    Sequence.fold_right (remove_gaps (Cost_matrix.Two_D.gap c2)) 
                    (create (length s2')) s2'
                in
                prepend res (Cost_matrix.Two_D.gap c2);
                res
            in
            let get_closest v i =
                let v' = get s1' i in
                Cost_matrix.Two_D.get_closest c2 v' v 
            in
            remove_gaps (Sequence.mapi get_closest s2'), cst

     let full_median_2 a b cm _ =
         let m,_,_,_ = align_2 a b cm in
         m
    end
end

module DOS = struct


    let rec bitset_to_seq gap set =
        match set with
        | Raw x -> x
        | Packed (len, bitset, seq) ->
                let res = Sequence.create len in
                let seq = bitset_to_seq gap seq in
                let cnt = ref (Sequence.length seq) in
                for i = len - 1 downto 0 do
                    let to_prepend =
                        if BitSet.is_set bitset i then 
                            let () = decr cnt in
                            Sequence.get seq !cnt
                        else gap
                    in
                    Sequence.prepend res to_prepend
                done;
                (*
                assert (
                    if not (res = x) then begin
                        Printf.printf "I am getting \n%s\n%s\n%!"
                        (Sequence.to_string res Alphabet.nucleotides)
                        (Sequence.to_string x Alphabet.nucleotides);
                        false end else true);
                *)
                res

    let seq_to_bitset gap seq seqo =
        let len = Sequence.length seq in
        let set = BitSet.create len in
        for i = 0 to len - 1 do
            if gap <> Sequence.get seq i then BitSet.set set i;
        done;
        Packed (len, set, seqo)

    type do_single_sequence = {
            sequence : Sequence.s;
            aligned_children : packed_algn * packed_algn * packed_algn;
            costs : cost_tuple;
            position : int;
            delimiters: int list;
        }

    let update_do_single_sequence oldone newseq newdeli =
        { oldone with sequence = newseq; delimiters = newdeli }
    
    let safe_reverse x = 
        { x with sequence = Sequence.safe_reverse x.sequence }

    let create seq = {
        sequence = seq;
        aligned_children = Raw seq, Raw seq, Raw seq;
        costs = { min = 0.0; max = 0.0 };
        position = 0;
        delimiters = []
    }

    let print do_single_seq = 
        Sequence.print stdout do_single_seq.sequence Alphabet.nucleotides

    let to_union a = Sequence.Unions.leaf a.sequence

    let make_cost tmpcost = 
        let tmpcost = float_of_int tmpcost in
        {min = tmpcost; max= tmpcost}

    let readjust mode h ch1 ch2 parent mine use_ukk =
        let c2 = h.c2 in
        let gap = Cost_matrix.Two_D.gap c2 in
        let res, cost = 
            let algn s1 s2 =
                let cm = c2 in
                match Cost_matrix.Two_D.affine cm with
                | Cost_matrix.Affine _ ->
                        let median,cost =
                            if use_ukk then 
                                let s1',s2',c = 
                                    Sequence.NewkkAlign.align_2 s1 s2 cm
                                    Sequence.NewkkAlign.default_ukkm
                                in
                                Sequence.median_2 s1' s2' cm, c
                            else
                                let m, _, _, cost, _ = 
                                    Sequence.Align.align_affine_3 s1 s2 cm in
                                m,cost
                        in
                        let m = Sequence.select_one median cm in
                        create m, cost
                | _ ->
                    let s1', s2', c =
                            if use_ukk then
                            Sequence.NewkkAlign.align_2 s1 s2 cm
                            Sequence.NewkkAlign.default_ukkm
                            else
                        Sequence.Align.align_2 ~first_gap:true s1 s2 cm Matrix.default
                    in
                    let median = Sequence.median_2 s1' s2' cm in
                    let a = Sequence.select_one median cm in 
                    create a, c 
            in
            let empty1 = Sequence.is_empty ch1.sequence gap 
            and empty2 = Sequence.is_empty ch2.sequence gap 
            and emptypar = Sequence.is_empty parent.sequence gap in
            match empty1, empty2, emptypar with
            | false, false, false ->
                    let tmpcost, seqm, changed =
                        match mode with
                        | `ThreeD _ -> 
                                Printf.printf "iter = 3d, call readjust_3d\n%!";
                                Sequence.Align.readjust_3d ch1.sequence ch2.sequence
                                mine.sequence h.c2 h.c3 parent.sequence
                        | `ApproxD _ ->
                                Printf.printf "iter=approx, call readjust\n%!";
                                Sequence.readjust ch1.sequence ch2.sequence
                                mine.sequence h.c2 parent.sequence use_ukk
                    in
                    let rescost = make_cost tmpcost in
                    { mine with sequence = seqm; costs = rescost },
                    tmpcost
            | true, true, _ -> ch1, 0
            | true, _, true -> ch1, 0
            | _, true, true -> ch2, 0
            | false, false, true ->
                    algn ch1.sequence ch2.sequence
            | false, true, false ->
                    algn ch1.sequence parent.sequence
            | _, false, false ->
                    algn ch2.sequence parent.sequence
        in
        0 <> compare res.sequence mine.sequence, res, cost

    let to_single h parent mine =
        let debug = false in
        if debug then Printf.printf "seqCS.DOS.to_single,%!";
        let gap = Cost_matrix.Two_D.gap h.c2 in
        let use_ukk = match !Methods.algn_mode with
            | `Algn_Newkk  -> true
            | `Algn_Normal -> false
        in
        if Sequence.is_empty mine.sequence gap then
            create mine.sequence, 0
        else
            let parent =
                if Sequence.is_empty parent.sequence gap then 
                    mine.sequence
                else parent.sequence 
            in
            let seqm, tmpcost = 
                if use_ukk then
                    Sequence.NewkkAlign.closest parent mine.sequence h.c2 Sequence.NewkkAlign.default_ukkm
                else
                    Sequence.Align.closest parent mine.sequence h.c2 Matrix.default 
            in
            let rescost = make_cost tmpcost in
            if debug then begin
                Printf.printf "parent seq = %!";
                Sequence.printseqcode parent;
                Printf.printf " seq with ambiguity = %!";
                Sequence.printseqcode mine.sequence;
                Printf.printf " ==> seq without ambiguity = %!";
                Sequence.printseqcode seqm;
            end;
            { mine with sequence = seqm; costs = rescost }, tmpcost

    let median alph code h a b use_ukk =
        let debug = false in
        if debug then Printf.printf "seqCS.DOS.median,use_ukk=%b\n%!" use_ukk; 
        let gap = Cost_matrix.Two_D.gap h.c2 in
        (* above are debug functions *)
        if Sequence.is_empty a.sequence gap then
            create b.sequence, 0 
        else if Sequence.is_empty b.sequence gap then
            create a.sequence, 0
        else 
            let seqm, tmpa, tmpb, tmpcost, seqmwg =
                match Cost_matrix.Two_D.affine h.c2 with
                | Cost_matrix.Affine _ ->
                        (**)
                        if use_ukk then
                            let tmpa,tmpb,tmpcost = 
                                Sequence.NewkkAlign.align_2 a.sequence b.sequence 
                                h.c2 Sequence.NewkkAlign.default_ukkm in
                            let seqm = Sequence.Align.ancestor_2 tmpa tmpb h.c2 in
                            let seqmwg = 
                                Sequence.median_2_with_gaps tmpa tmpb h.c2 
                            in
                            seqm, tmpa, tmpb, tmpcost, seqmwg
                        else
                            Sequence.Align.align_affine_3 a.sequence b.sequence h.c2
                        (**)
                        (* compare cost and alignment
                        let notthesame s1 s2 =
                            if (Sequence.compare s1 s2)<>0 then true
                            else false
                        in
                        let tmpa,tmpb,tmpcost = 
                                Sequence.NewkkAlign.align_2 a.sequence b.sequence 
                                h.c2 Sequence.NewkkAlign.default_ukkm in
                        let seqm = Sequence.Align.ancestor_2 tmpa tmpb h.c2 in
                        let seqmwg = Sequence.median_2_with_gaps tmpa tmpb h.c2 
                        in
                        let seqm2, tmpa2, tmpb2, tmpcost2, seqmwg2 =
                        Sequence.Align.align_affine_3 a.sequence b.sequence h.c2
                        in
                        (*just compare the cost*)
                        let oc1 =  open_out "algn.cost" in
                        Printf.fprintf oc1 "%d%!" tmpcost;
                        close_out oc1;
                        let oc2 =  open_out "newkkonen.cost" in
                        Printf.fprintf oc2 "%d%!" tmpcost2;
                        close_out oc2;
                        (* compare the traceback*)
                        if (notthesame seqm seqm2)||(notthesame seqmwg
                        seqmwg2)||(notthesame tmpa tmpa2)||(notthesame tmpb
                        tmpb2) then begin
                            Printf.printf "diff between newkkonen and normal affine alignment:\n%!";
                            Printf.printf "input seq(len1=%d,len2=%d) :\n%!"
                            (Sequence.length a.sequence) (Sequence.length b.sequence);
                            Sequence.print stdout a.sequence Alphabet.nucleotides;print_newline();
                            Sequence.print stdout b.sequence Alphabet.nucleotides;print_newline();
                            Printf.printf "tmpa from algn.c, tmpa2 from newkkonen.c:\n%!";
                            Sequence.print stdout tmpa Alphabet.nucleotides;print_newline();
                            Sequence.print stdout tmpa2 Alphabet.nucleotides;print_newline();
                            Printf.printf "tmpb,tmpb2:\n%!";
                            Sequence.print stdout tmpb Alphabet.nucleotides;print_newline();
                            Sequence.print stdout tmpb2 Alphabet.nucleotides;print_newline();
                            Printf.printf "seqm with gap :\n%!";
                            Sequence.print stdout seqmwg Alphabet.nucleotides;print_newline();
                            Sequence.print stdout seqmwg2 Alphabet.nucleotides;print_newline();
                            Printf.printf "seqm : \n%!";
                            Sequence.print stdout seqm Alphabet.nucleotides;print_newline();
                            Sequence.print stdout seqm2 Alphabet.nucleotides;print_newline();
                            assert(false);
                        end;
                        seqm, tmpa, tmpb, tmpcost, seqmwg
                        compare cost and alignment *)
                | _ ->
                        let tmpa, tmpb, tmpcost =
                           if use_ukk then 
                                Sequence.NewkkAlign.align_2 a.sequence b.sequence
                            h.c2 Sequence.NewkkAlign.default_ukkm
                        (*in
                        let tmpa2, tmpb2, tmpcost2 =*)
                            else 
                                Sequence.Align.align_2 a.sequence b.sequence h.c2
                            Matrix.default
                        in
                        let seqm = Sequence.Align.ancestor_2 tmpa tmpb h.c2 in
                        let seqmwg = 
                            Sequence.median_2_with_gaps tmpa tmpb h.c2 
                        in
                        (*
                        let seqm2 = Sequence.Align.ancestor_2 tmpa tmpb h.c2 in
                        let seqmwg2 = 
                            Sequence.median_2_with_gaps tmpa2 tmpb2 h.c2 
                        in
                        let notthesame s1 s2 =
                            if (Sequence.compare s1 s2)<>0 then true
                            else false
                        in
                        if (notthesame seqm seqm2)||(notthesame seqmwg
                        seqmwg2)||(notthesame tmpa tmpa2)||(notthesame tmpb
                        tmpb2) then begin
                            Printf.printf "diff between newkkonen and normal alignment:\n%!";
                            Printf.printf "input seq :\n%!";
                            Sequence.print stdout a.sequence Alphabet.nucleotides;print_newline();
                            Sequence.print stdout b.sequence Alphabet.nucleotides;print_newline();
                            Printf.printf "tmpa from algn.c, tmpa2 from newkkonen.c:\n%!";
                            Sequence.print stdout tmpa Alphabet.nucleotides;print_newline();
                            Sequence.print stdout tmpa2 Alphabet.nucleotides;print_newline();
                            Printf.printf "tmpb,tmpb2:\n%!";
                            Sequence.print stdout tmpb Alphabet.nucleotides;print_newline();
                            Sequence.print stdout tmpb2 Alphabet.nucleotides;print_newline();
                            Printf.printf "seqm with gap :\n%!";
                            Sequence.print stdout seqmwg Alphabet.nucleotides;print_newline();
                            Sequence.print stdout seqmwg2 Alphabet.nucleotides;print_newline();
                            Printf.printf "seqm : \n%!";
                            Sequence.print stdout seqm Alphabet.nucleotides;print_newline();
                            Sequence.print stdout seqm2 Alphabet.nucleotides;print_newline();
                            assert(false);
                        end;*)
                        seqm, tmpa, tmpb, tmpcost, seqmwg
            in
            let rescost = make_cost tmpcost in
            if debug then begin 
                let print_seqlist seq = 
                    (*Sequence.print stdout seq Alphabet.nucleotides;*)
                    Sequence.printseqcode seq;
                    print_newline();
                in
                Printf.printf "costs = (%f,%f), a/b= %!" rescost.min rescost.max; 
                print_seqlist a.sequence; print_seqlist b.sequence;
                print_seqlist tmpa; print_seqlist tmpb;
                Printf.printf "seqm: %!"; print_seqlist seqm;
            end;
            let ba = seq_to_bitset gap tmpa (Raw a.sequence)
            and bb = seq_to_bitset gap tmpb (Raw b.sequence) 
            and bm = seq_to_bitset gap seqmwg (Raw seqm) in
            { sequence = seqm; aligned_children = (ba, bb, bm); costs = rescost;
            position = 0; delimiters = a.delimiters}, tmpcost

    let median_3_no_union h p n c1 c2 use_ukk =
        let with_parent c =
            let s1, s2, costs =
                if use_ukk then 
                    Sequence.NewkkAlign.align_2 p.sequence c.sequence h.c2 Sequence.NewkkAlign.default_ukkm
                else 
                    Sequence.Align.align_2 p.sequence c.sequence h.c2 Matrix.default
            in
            { n with 
                sequence = Sequence.median_2 s1 s2 h.c2;
                costs = 
                    { min = float_of_int costs; 
                    max = float_of_int (Sequence.Align.max_cost_2 s1 s2 h.c2) }
            }, costs
        in
        let (res1, cost1) = with_parent c1
        and (res2, cost2) = with_parent c2 in
        let res = if cost1 < cost2 then res1 else res2 in
        let gap = Cost_matrix.Two_D.gap h.c2 in
        if gap <> (Sequence.get res.sequence 0) then 
            Sequence.prepend res.sequence gap;
        res

    let median_3_union h p n c1 c2 use_ukk =
        let gap = Cost_matrix.Two_D.gap h.c2 in
        let a, b, _ = n.aligned_children in
        let a = bitset_to_seq gap a
        and b = bitset_to_seq gap b in
        assert (Sequence.length a = Sequence.length b);
        let res = Sequence.Align.union a b h.c2 in
        let a, b, cost = 
            if use_ukk then 
                Sequence.NewkkAlign.align_2 p.sequence res h.c2 Sequence.NewkkAlign.default_ukkm
            else
                Sequence.Align.align_2 p.sequence res h.c2 Matrix.default 
        in
        let res = 
            let res = Sequence.median_2 a b h.c2 in
            if gap <> Sequence.get res 0 then
                Sequence.prepend res gap;
            res
        in
        let rescost = 
            { min = float_of_int cost; 
            max = float_of_int (Sequence.Align.max_cost_2 a b h.c2) } in
        { n with sequence = res; costs = rescost }

    let distance alph h missing_distance a b use_ukk =
        let gap = Cost_matrix.Two_D.gap h.c2 in 
        if Sequence.is_empty a.sequence gap || 
            Sequence.is_empty b.sequence gap then missing_distance
        else
IFDEF USE_VERIFY_COSTS THEN
            let seqa, seqb, cost = 
                if use_ukk then
                Sequence.NewkkAlign.align_2 a.sequence b.sequence h.c2 Sequence.NewkkAlign.default_ukkm
                else
                    Sequence.Align.align_2 a.sequence b.sequence h.c2 
                Matrix.default 
            in
            let () = 
                assert (
                    let real_cost = 
                        Sequence.Align.verify_cost_2 cost seqa seqb h.c2 
                    in
                    if cost < real_cost then
                        let () = 
                            Printf.printf 
                            "Failed alignment between \
                            \n%s\nand\n%s\nwith claimed cost %d \
                            and real cost %d\n%!" 
                            (Sequence.to_string seqa alph)
                            (Sequence.to_string seqb alph)
                            cost
                            real_cost
                        in
                        false
                    else if real_cost < cost then begin
                        Printf.printf "Check this case\n%s\n%s\n has real cost %d and expected cost is %d\n%!"
                        (Sequence.to_string seqa alph) (Sequence.to_string seqb
                        alph) real_cost cost;
                        true
                    end else true
            ) 
            in
            cost
ELSE 
            let deltaw = 
                let tmp = 
                    (max (Sequence.length a.sequence) (Sequence.length
                    b.sequence)) 
                    - (min (Sequence.length a.sequence) (Sequence.length
                    b.sequence)) in
                if tmp > 8 then tmp 
                else 8
            in
            if use_ukk then
            Sequence.NewkkAlign.cost_2  ~deltaw a.sequence b.sequence h.c2
            Sequence.NewkkAlign.default_ukkm
            else
                Sequence.Align.cost_2 ~deltaw a.sequence b.sequence h.c2 
            Matrix.default
END

        let dist_2 h n a b use_ukk =
            let gap = Cost_matrix.Two_D.gap h.c2 in
            let tmp =
                if Sequence.is_empty a.sequence gap then
                    n.sequence
                else 
                    Sequence.Align.full_median_2 a.sequence
                    b.sequence h.c2 Matrix.default
            in
            let cost =
                if use_ukk then
                Sequence.NewkkAlign.cost_2 n.sequence tmp h.c2 Sequence.NewkkAlign.default_ukkm
                else
                    Sequence.Align.cost_2 n.sequence tmp h.c2
                Matrix.default
            in
            cost

        let compare a b =
            Sequence.compare a.sequence b.sequence
end

module PartitionedDOS = struct
    type fragment = 
        | DO of DOS.do_single_sequence
        | First of DOS.do_single_sequence
        | Last of DOS.do_single_sequence

    type partitioned_sequence = fragment array

    let create mode clip alph =
        let seqs =
            match mode with
            | `Partitioned seqs -> seqs
            | `AutoPartitioned (splits, seq) ->
                    let seqs = Sequence.split splits seq alph in
                    Array.of_list seqs 
        in
        let last = (Array.length seqs) - 1 in
        Array.mapi (fun pos x -> 
            if pos = 0 && clip = Data.Clip then
                First (DOS.create x)
            else if pos = last && clip = Data.Clip then
                Last (DOS.create x)
            else DO (DOS.create x)) seqs

    let empty clip size alph = 
        let s = Sequence.create 1 in
        Sequence.prepend s (Alphabet.get_gap alph);
        let cr = DOS.create s in
        Array.init size (fun x -> 
            if x = 0 && clip = Data.Clip then
                First cr
            else if x = size - 1 && Data.Clip = clip then
                Last cr
            else DO cr)

    let to_union a = Array.map (function 
        | Last x
        | First x
        | DO x -> DOS.to_union x) a

    let clip_n_fix s1 s2 = 
        (* We return the sequences with their corrections, and functions 
        * to correct the median and the pairwise alignments *)
        let identity x = x in
        let aux s1 s2 s1len s2len =
            assert (s1len >= s2len);
            if s1len > (s2len *. Sequence.Clip.fraction) then
                (* We eliminate enough to make the length difference 
                * only 10% of the overall length *)
                let maxlen = truncate (s2len *. Sequence.Clip.fraction) in
                let dif = (truncate s1len) - maxlen in
                let cnt = ref (truncate s1len) in
                let s1' = Sequence.init (fun x ->
                    decr cnt;
                    if x = 0 then Sequence.get s1.DOS.sequence 0
                    else Sequence.get s1.DOS.sequence !cnt) maxlen 
                in
                let s1' = { s1 with DOS.sequence = s1' } in
                let add_bottom x =
                    let len = Sequence.length x.DOS.sequence in
                    let newx = Sequence.create (len + dif) in
                    for i = (len - 1) downto 1 do
                        Sequence.prepend newx 
                        (Sequence.get x.DOS.sequence i);
                    done;
                    { x with DOS.sequence = newx }
                in
                let fixs1 x =
                    let newx = add_bottom x in
                    for i = dif downto 0 do
                        Sequence.prepend newx.DOS.sequence 
                        (Sequence.get s1.DOS.sequence i);
                    done;
                    newx
                in
                let fixs2 x =
                    let newx = add_bottom x in
                    let gap = Sequence.get s1.DOS.sequence 0 in
                    for i = dif downto 0 do
                        Sequence.prepend newx.DOS.sequence gap;
                    done;
                    newx
                in
                let fixmedian x = 
                    let newx = add_bottom x in
                    for i = dif downto 0 do
                        Sequence.prepend newx.DOS.sequence
                        (Sequence.get s1.DOS.sequence i)
                    done;
                    newx
                in
                assert (s1.DOS.sequence = (fixs1 s1').DOS.sequence);
                s1', s2, fixs1, fixs2, fixmedian, fixs1, identity
            else 
                s1, s2, identity, identity, identity, identity, identity
        in
        let s1len = float_of_int (Sequence.length s1.DOS.sequence)
        and s2len = float_of_int (Sequence.length s2.DOS.sequence) in
        if s1len >= s2len then aux s1 s2 s1len s2len
        else 
            let a, b, c, d, e, f, g = aux s2 s1 s2len s1len in
            b, a, d, c, e, g, f

    let readjust mode h ch1 ch2 parent mine use_ukk =
        let cost = ref 0 in
        let is_0 = ref true in
        let res =
            Array.init (Array.length ch1) (fun i ->
                match ch1.(i), ch2.(i), parent.(i), mine.(i) with
                | _, _, _, First _ 
                | _, _, First _, _ 
                | _, First _, _, _ 
                | First _, _, _, _ 
                | _, _, _, Last _ 
                | _, _, Last _, _ 
                | _, Last _, _, _ 
                | Last _, _, _, _  ->
                        mine.(i)
                | DO ch1, DO ch2, DO parent, DO mine ->
                        let x, y, z = 
                            DOS.readjust mode h ch1 ch2 parent mine use_ukk
                        in
                        is_0 := !is_0 && x;
                        cost := !cost + z;
                        DO y)
        in
        !is_0, res, !cost

    let identity x = x

    let to_single h parent mine =
        let tmp_cost = ref 0 in
        let gap = Cost_matrix.Two_D.gap h.c2 in
        let to_single_tip reverse parent mine =
            if Sequence.is_empty mine.DOS.sequence gap then
                mine
            else 
                let parent = 
                    if Sequence.is_empty parent.DOS.sequence gap then
                        mine
                    else parent
                in
                let parent = reverse parent
                and mine = reverse mine in
                let parent', mine', _, _, _, _, fixmine = 
                    clip_n_fix parent mine 
                in
                let x, y = DOS.to_single h parent' mine' in
                tmp_cost := y + !tmp_cost;
                reverse (fixmine x)
        in
        let res = Array.init (Array.length mine) (fun i ->
            match parent.(i), mine.(i) with
            | DO parent, First mine
            | First parent, First mine -> 
                    First (to_single_tip identity parent mine)
            | First parent, DO mine ->
                    DO (to_single_tip identity parent mine)
            | DO parent, Last mine
            | Last parent, Last mine -> 
                    Last (to_single_tip DOS.safe_reverse parent mine)
            | Last parent, DO mine ->
                    DO (to_single_tip DOS.safe_reverse parent mine)
            | DO parent, DO mine ->
                    let x, y = DOS.to_single h parent mine in
                    tmp_cost := y + !tmp_cost;
                    DO x
            | Last _, First _
            | First _, Last _ -> assert false)
        in
        res, !tmp_cost



    let median alph code h a b use_ukk =
        let tmpcost = ref 0 in
        let median_of_tip reverse a b =
            let gap = Cost_matrix.Two_D.gap h.c2 in 
            if Sequence.is_empty a.DOS.sequence gap then
                DOS.create b.DOS.sequence
            else if Sequence.is_empty b.DOS.sequence gap then
                DOS.create a.DOS.sequence
            else 
                let a = reverse a
                and b = reverse b in
                let a, b, fixma, fixmb, fixmedian, fixa, fixb = 
                    clip_n_fix a b in
                let seqm, tmpa, tmpb, cost, seqmwg =
                    match Cost_matrix.Two_D.affine h.c2 with
                    | Cost_matrix.Affine _ ->
                            if use_ukk then
                                let tmpa,tmpb,tmpcost =
                                    Sequence.NewkkAlign.align_2 a.DOS.sequence 
                                b.DOS.sequence h.c2 Sequence.NewkkAlign.default_ukkm in
                                let seqm = Sequence.Align.ancestor_2 tmpa tmpb h.c2 in
                                let seqmwg = 
                                    Sequence.median_2_with_gaps tmpa tmpb h.c2 
                                in
                                seqm, tmpa, tmpb, tmpcost, seqmwg
                            else
                                Sequence.Align.align_affine_3 a.DOS.sequence b.DOS.sequence h.c2
                    | _ ->
                            let tmpa, tmpb, cost = 
                                if use_ukk then
                                Sequence.NewkkAlign.align_2 a.DOS.sequence 
                                b.DOS.sequence h.c2 Sequence.NewkkAlign.default_ukkm
                                else
                                    Sequence.Align.align_2 a.DOS.sequence 
                                b.DOS.sequence h.c2 Matrix.default
                            in
                            let seqm = Sequence.Align.ancestor_2 tmpa tmpb h.c2 in
                            let seqmwg = 
                                Sequence.median_2_with_gaps tmpa tmpb h.c2 
                            in
                            seqm, tmpa, tmpb, cost, seqmwg
                in
                let rescost = 
                    DOS.make_cost 
                    (Sequence.Clip.correct_distance gap h.c2 tmpa tmpb cost) 
                in
                let a = fixa a in
                let b = fixb b in
                let a = reverse a
                and b = reverse b in
                let tmpa = reverse (fixma (DOS.create tmpa))
                and tmpb = reverse (fixmb (DOS.create tmpb))
                and seqm = reverse (fixmedian (DOS.create seqm))
                and seqmwg = reverse (fixmedian (DOS.create seqmwg)) in
                let ba = 
                    DOS.seq_to_bitset gap tmpa.DOS.sequence (Raw a.DOS.sequence)
                and bb = 
                    DOS.seq_to_bitset gap tmpb.DOS.sequence (Raw b.DOS.sequence) 
                and bm = 
                    DOS.seq_to_bitset gap seqmwg.DOS.sequence 
                    (Raw seqm.DOS.sequence) 
                in
                assert (
                    if not (tmpa.DOS.sequence = DOS.bitset_to_seq gap ba) then
                        begin
                            Sequence.printDNA tmpa.DOS.sequence;
                            Sequence.printDNA (DOS.bitset_to_seq gap ba);
                            false
                        end else true);
                assert (
                    if not (tmpb.DOS.sequence = DOS.bitset_to_seq gap bb) then
                        begin
                            Sequence.printDNA tmpb.DOS.sequence;
                            Sequence.printDNA (DOS.bitset_to_seq gap bb);
                            false
                        end else true);
                assert (seqmwg.DOS.sequence = DOS.bitset_to_seq gap bm);
                tmpcost := !tmpcost + cost;
                { DOS.sequence = seqm.DOS.sequence; 
                aligned_children = (ba, bb, bm); 
                costs = rescost; position = 0; 
                (*tmp solution, fix this later*)
                delimiters = a.DOS.delimiters }
        in
        let res = 
            Array.init (Array.length a) (fun i ->
                match a.(i), b.(i) with
                | DO a, First b
                | First a, DO b -> DO (median_of_tip identity a b)
                | First a, First b -> First (median_of_tip identity a b)
                | Last a, DO b
                | DO a, Last b -> DO (median_of_tip DOS.safe_reverse a b)
                | Last a, Last b -> Last (median_of_tip DOS.safe_reverse a b)
                | DO a, DO b ->
                        let x, y = DOS.median alph code h a b use_ukk in
                        tmpcost := y + !tmpcost;
                        DO x
                | First _, _
                | Last _, _ -> assert false
                )
        in
        res, !tmpcost

    let median_3_no_union h p n c1 c2 use_ukk =
        Array.init (Array.length n) (fun i ->
            match p.(i), n.(i), c1.(i), c2.(i) with
            | Last p, Last n, Last c1, Last c2 
            | First p, First n, First c1, First c2 ->
                    assert (false) (* TODO *)
            | DO p, DO n, DO c1, DO c2 ->
                    DO (DOS.median_3_no_union h p n c1 c2 use_ukk)
            | Last _, _, _, _
            | First _, _, _, _
            | DO _, _, _, _ -> assert false)

    let median_3_union h p n c1 c2 use_ukk=
        Array.init (Array.length n) (fun i ->
            match p.(i), n.(i), c1.(i), c2.(i) with
            | First p, First n, First c1, First c2 -> 
                    First n
            | Last p, Last n, Last c1, Last c2 ->
                    Last n
            | DO p, DO n, DO c1, DO c2 ->
                    DO (DOS.median_3_union h p n c1 c2 use_ukk)
            | First _, _, _, _
            | Last _, _, _, _
            | DO _, _, _, _ -> assert false)


    let distance alph h missing_distance a b use_ukk =
        let dist = ref 0 in
        let distance_tip reverse a b =
            let a = reverse a
            and b = reverse b in
            let a, b, _, _, _, _, _ = clip_n_fix a b in
            dist := !dist + 
                Sequence.Clip.corrected_distance h.c2 missing_distance
                a.DOS.sequence b.DOS.sequence use_ukk
        in
        for i = (Array.length a) - 1 downto 0 do
            match a.(i), b.(i) with
            | First a, DO b
            | DO a, First b
            | First a, First b -> distance_tip identity a b
            | DO a, Last b 
            | Last a, DO b
            | Last a, Last b -> distance_tip DOS.safe_reverse a b
            | DO a, DO b -> 
                    dist := !dist + DOS.distance alph h missing_distance a b use_ukk
            | Last _, First _
            | First _, Last _ -> assert false
        done;
        !dist

    let merge arr = 
        let total_len, min, max = 
            Array.fold_left (fun (total_len, mincost, maxcost) x -> 
                match x with
                | Last x 
                | First x
                | DO x ->
                        total_len + Sequence.length x.DOS.sequence,
                        mincost +. x.DOS.costs.min,
                        maxcost +. x.DOS.costs.max)
            (0, 0., 0.) arr
        in
        let res = Sequence.create total_len in
        for i = (Array.length arr) - 1 downto 0 do
            match arr.(i) with
            | Last seq
            | First seq
            | DO seq ->
                    let seq = seq.DOS.sequence in
                    let len = Sequence.length seq in
                    for j = (len - 1) downto 1 do
                        Sequence.prepend  res (Sequence.get seq j);
                    done;
                    if i = 0 then Sequence.prepend res (Sequence.get seq 0) ;
        done;
        { DOS.sequence = res;
        aligned_children = Raw res, Raw res, Raw res;
        costs = {min = min; max = max};
        position = 0;
        (* fix this later *)
        delimiters = []
        }


    let dist_2 h n a b use_ukk =
        let total = ref 0 in
        let do_one n a b =
            let a, b, _, _, fixmedian, _, _ = clip_n_fix a b in
            let gap = Cost_matrix.Two_D.gap h.c2 in
            let tmp =
                if Sequence.is_empty a.DOS.sequence gap then
                    n.DOS.sequence
                else 
                    Sequence.Align.full_median_2 a.DOS.sequence
                    b.DOS.sequence h.c2 Matrix.default
            in
            let tmp = fixmedian { a with DOS.sequence = tmp } in
            let cost = 
                let n, tmp, _, _, _, _, _ = clip_n_fix n tmp in
                if use_ukk then
                Sequence.NewkkAlign.cost_2 n.DOS.sequence 
                tmp.DOS.sequence h.c2 Sequence.NewkkAlign.default_ukkm
                else
                    Sequence.Align.cost_2 n.DOS.sequence 
                tmp.DOS.sequence h.c2 Matrix.default
            in
            total := !total + cost
        in
        for i = (Array.length n) - 1 downto 0 do
            match n.(i), a.(i), b.(i) with
            | DO n, First a, First b
            | First n, DO a, First b
            | First n, First a, DO b
            | DO n, DO a, First b
            | DO n, First a, DO b
            | First n, DO a, DO b
            | First n, First a, First b -> do_one n a b
            | DO n, Last a, Last b
            | Last n, DO a, Last b
            | Last n, Last a, DO b
            | DO n, DO a, Last b
            | DO n, Last a, DO b
            | Last n, DO a, DO b
            | Last n, Last a, Last b ->
                    let n = DOS.safe_reverse n
                    and a = DOS.safe_reverse a
                    and b = DOS.safe_reverse b in
                    do_one n a b
            | DO n, DO a, DO b ->
                    total := !total + (DOS.dist_2 h n a b use_ukk)
            | First _, _, _
            | Last _, _, _
            | DO _, _, _ -> assert false
        done;
        !total

    let compare a b =
        let compare = ref 0 in
        try
            for i = 0 to (Array.length a) - 1 do
                match a.(i), b.(i) with
                | First a, First b 
                | Last a, Last b
                | DO a, DO b -> 
                        compare := DOS.compare a b;
                        if !compare <> 0 then raise Exit
                | _, DO _ 
                | First _, _ -> 
                        compare := -1;
                        raise Exit
                | DO _, _
                | Last _, First _ -> 
                        compare := 1;
                        raise Exit
            done;
            0
        with 
        | Exit -> !compare

    let tabu_distance a =
        Array.fold_left (fun sum y -> 
            match y with
            | First y
            | Last y
            | DO y -> y.DOS.costs.max +. sum) 0.0 a

    let encoding enc arr =
        Array.fold_left (fun sum x ->
            match x with
            | First x
            | Last x
            | DO x -> sum +. (Sequence.encoding enc x.DOS.sequence)) 0. arr


end


let max_float = float_of_int (max_int / 20) 


type sequence_characters =
    | General_Prealigned of GenNonAdd.gnonadd_sequence
    | Heuristic_Selection of DOS.do_single_sequence
    | Partitioned of PartitionedDOS.partitioned_sequence
    (*| Relaxed_Lifted of RL.t_w_seqtbl*)

(** A sequence character type. *)
type t = { 
    (** The sequences that form the set *)
    characters : sequence_characters array;
    codes : int array;
    total_cost : float;             (** The total cost of the character set *)
    alph : Alphabet.a;              (** The alphabet of the sequence set *)
    code : int;                     (** The set code *)
    heuristic : heuristic;          (** The heuristic to be used *)
    priority : int list;            (** The information ordering *)
}
(*
let is_fixedstates x =
    match x.characters.(0) with
    | Relaxed_Lifted _ -> true
    | _ -> false
*)

let print in_data =
    let seq_chr_arr = in_data.characters in
    Array.iter (fun item ->
     match item with
        | General_Prealigned x -> 
                Printf.printf "General_Prealigned,no print function yet%!"
        | Heuristic_Selection x -> 
                Printf.printf "Heuristic_Selection,%!";
                DOS.print x;
        (*| Relaxed_Lifted x -> Printf.printf "Relaxed_Lifted,no print functionyet%!"*)
        | Partitioned x -> Printf.printf "Partitioned,no print function yet%!"
    ) seq_chr_arr;
    Printf.printf "\n%!"

(* return 1 if we are dealing with this kind of SeqCS data for multi-chromosome.*)
let is_available in_data =
    match (in_data.characters).(0) with
    | General_Prealigned _ -> 0
    | Heuristic_Selection _ -> 0
    | Partitioned _ -> 0
    (*| Relaxed_Lifted _ -> 0*)
    
let flatten t_lst = 
    List.map ( fun x -> 
    let gapcode = Alphabet.get_gap x.alph in
    let seqchar_arr =  x.characters in
    List.map (fun x -> match x with
        | Heuristic_Selection dos_single_seq ->
                let ori_seq = dos_single_seq.DOS.sequence in
                (* for multi-chromosome,
                * delete the "beginning position" gap of sequence 
                * we are going to concat these seq into one later in 
                * "transform_multi_chromosome" of node.ml,
                * *)
                let new_seq = Sequence.fold_righti 
                (fun acc pos item ->
                    if (pos=0)&&(item=gapcode) then acc
                    else
                        let _ = Sequence.prepend acc item in
                        acc
                ) (Sequence.create ((Sequence.length ori_seq)-1)) ori_seq 
                in
                new_seq   
        | _ -> failwith "seqCS.flatten,only work for Heuristic_Selection type"
        ) (Array.to_list seqchar_arr)
    ) t_lst 

let update_t oldt (newseqlst:Sequence.s list) (delimiterslst: int list list)  =
    let i = ref 0 in
    (* do we need to patch a gap to the begging of the newseq? *)
    (* if we do so, the delimiter list should be modified as well*)
    (* why do we have a gap at the beginning of each sequence for SeqCS, 
    * but not for BreakinvCS? *)
    (* let gapcode = Alphabet.get_gap x.alph in *)
    let old_characters_list = Array.to_list oldt.characters in
    let new_characters_list = List.map (fun seq_char ->
        match seq_char with
        | Heuristic_Selection oldone ->
                let newone = DOS.update_do_single_sequence oldone  
                (List.nth newseqlst !i) (List.nth delimiterslst !i) in
                i := !i +1;
                Heuristic_Selection newone
        | _ -> failwith "not there yet"
    ) old_characters_list
    in
    let new_characters = Array.of_list new_characters_list in
    {oldt with characters = new_characters }

module Union = struct
    (* The union of sequences *)
    type ustr = {
        unions : union_element option array;
        u_c2 : Cost_matrix.Two_D.m;
        u_alph : Alphabet.a;
        u_codes : int array;
    }


    type u = ustr option

    let compare_union a b =
        match a, b with
        | Some a, Some b ->
                Array_ops.fold_right_2 
                (fun acc x y ->
                    match acc with
                    | 0 -> 
                            (match x, y with
                            | Some (Single x), Some (Single y) -> Sequence.Unions.compare x y
                            | Some (Array x), Some (Array y) -> 
                                    (let res = ref 0 in
                                    try
                                        for i = 0 to (Array.length x) - 1 do
                                            res := Sequence.Unions.compare x.(i)
                                            y.(i);
                                            if !res <> 0 then raise Exit
                                        done;
                                        0
                                    with
                                    Exit -> !res)
                            | Some (Single _), Some (Array _) -> -1
                            | Some (Array _), Some (Single _) -> 1
                            | None, None -> 0
                            | Some _, None -> 1
                            | None, Some _ -> -1)
                    | x -> x)
                0
                a.unions
                b.unions
        | None, None -> 0
        | Some _, None -> -1
        | None, Some _ -> 1

    let cardinal_union ua = 
        match ua with
        | None -> 0
        | Some ua -> Array.length ua.unions

    let poly_saturation x v =
        match x with
        | None -> 0.0
        | Some x ->
                let single_saturation (acc, len) x =
                    let nlen = Sequence.length x.Sequence.Unions.seq 
                    and sat = 
                        Sequence.poly_saturation x.Sequence.Unions.seq v 
                    in
                    acc +. (sat *. (float_of_int nlen)), len + nlen
                in
                let sat, len =
                    Array.fold_left (fun ((acc, len) as acc1) x ->
                        match x with
                        | None -> acc1 
                        | Some x -> 
                                match x with
                                | Single x -> single_saturation acc1 x
                                | Array x -> 
                                        Array.fold_left single_saturation 
                                        acc1 x)
                    (0.0, 0) x.unions 
                in
                sat /. (float_of_int len)

    let union self ua ub =
        match ua, ub with
        | Some ua, Some ub ->
                let c2 = ua.u_c2 in
                let gap = Cost_matrix.Two_D.gap c2 in
                let do_one self uniona unionb =
                    if Sequence.is_empty uniona.Sequence.Unions.seq gap then
                        unionb 
                    else if Sequence.is_empty unionb.Sequence.Unions.seq gap 
                    then
                        uniona 
                    else 
                        let tmpa, tmpb, tmpc = self.DOS.aligned_children in
                        let tmpa = DOS.bitset_to_seq gap tmpa 
                        and tmpb = DOS.bitset_to_seq gap tmpb 
                        and tmpc = DOS.bitset_to_seq gap tmpc in
                        (Sequence.Unions.union tmpa tmpb tmpc uniona unionb c2)
                in
                let union uniona unionb self =
                    match self with
                    | General_Prealigned _ -> failwith "write union later"
                    | Heuristic_Selection self ->
                        (match uniona, unionb with
                        | Some (Single uniona), Some (Single unionb) ->
                                Some (Single (do_one self uniona unionb))
                        | _ -> assert false)
                    (*| Relaxed_Lifted _ -> None*)
                    | Partitioned self -> 
                            (match uniona, unionb with
                            | Some (Array uniona), Some (Array unionb) ->
                                    let arr = 
                                        Array_ops.map_3 (fun a b c ->
                                            let a = 
                                                match a with
                                                | PartitionedDOS.Last a 
                                                | PartitionedDOS.First a
                                                | PartitionedDOS.DO a -> a
                                            in
                                            do_one a b c)
                                            self uniona unionb
                                    in
                                    Some (Array arr)
                            | _ -> assert false)
                in
                let union = 
                    Array_ops.map_3 union ua.unions ub.unions
                    self.characters
                in
                Some { ua with unions = union }
        | None, None -> None
        | _, _ -> failwith "SeqCS.union"

    let distance_union a b = 
        match a, b with
        | Some a, Some b ->
                let use_ukk = match !Methods.algn_mode with
                    | `Algn_Newkk  -> true
                    | `Algn_Normal -> false
                in
                let sub_factor = 
                    match Cost_matrix.Two_D.affine a.u_c2 with
                    | Cost_matrix.Affine _ -> 0.8
                    | _ -> 1.0 
                in
                let gap = Cost_matrix.Two_D.gap a.u_c2 in
                let distance =
                    let one_distance acc seqa seqb = 
                        let seqa = seqa.Sequence.Unions.seq
                        and seqb = seqb.Sequence.Unions.seq in
                        if Sequence.is_empty seqa gap || 
                        Sequence.is_empty seqb gap then
                            acc
                        else
                            let deltaw = 
                                let tmp = 
                                    (max (Sequence.length seqa) 
                                    (Sequence.length seqb)) -
                                    (min (Sequence.length seqa) 
                                    (Sequence.length seqb)) 
                                in
                                if tmp > 8 then tmp 
                                else 8
                            in
                            acc +.
                            (sub_factor *. 
                            (let d = 
                                if use_ukk then 
                                    Sequence.NewkkAlign.cost_2 ~deltaw:deltaw
                                seqa seqb a.u_c2
                                Sequence.NewkkAlign.default_ukkm
                                else
                                    Sequence.Align.cost_2 ~deltaw:deltaw
                                seqa seqb a.u_c2 Matrix.default in
                            float_of_int d))
                    in
                    Array_ops.fold_right_2 (fun acc seqa seqb ->
                        match seqa, seqb with
                        | Some (Single seqa), Some (Single seqb) -> 
                                one_distance acc seqa seqb
                        | Some (Array seqa), Some (Array seqb) ->
                                Array_ops.fold_right_2 one_distance acc seqa seqb
                        | None, None -> acc
                        | _ -> assert false)
                    0. a.unions b.unions
                in
                distance
        | None, None -> 0.0
        | Some _, _ 
        | None, _ -> failwith "Impossible?"

    let get_sequence_union code x =
        match x with
        | None -> None
        | Some x ->
                let res = ref (-1) in
                try 
                    for i = (Array.length x.u_codes) - 1 downto 0 do
                        if x.u_codes.(i) = code then begin
                            res := i;
                            raise Exit
                        end;
                    done;
                    None
                with
                | Exit -> (x.unions.(!res))

end

let cardinal x = Array.length x.codes

let empty code c2 alph =
    let set = 
        {
            characters = [||];
            codes = [||];
            total_cost = 0.0;
            alph = alph;
            code = code;
            priority = [];
            heuristic = make_default_heuristic c2;
        }
    in
    set

let to_union a = 
    if a.alph = Alphabet.nucleotides then
        let new_unions = Array.map (function
            | General_Prealigned _ -> failwith "write to_union later"
            | Heuristic_Selection x -> Some (Single (DOS.to_union x))
            (*| Relaxed_Lifted _ -> None*)
            | Partitioned s -> 
                    Some (Array (Array.map (function
                        | PartitionedDOS.Last x
                        | PartitionedDOS.First x
                        | PartitionedDOS.DO x -> DOS.to_union x) s))) 
        a.characters in
        Some { Union.unions = new_unions;
            u_c2 = a.heuristic.c2;
            u_alph = a.alph;
            u_codes = a.codes; 
        }
    else None

let to_string a =
    let builder acc code seq =
        let code = string_of_int code 
        and seq = 
            match seq with
            | General_Prealigned seq ->
                    Sequence.to_formater seq.GenNonAdd.seq a.alph
            | Partitioned seq ->
                    let seq = PartitionedDOS.merge seq in
                    Sequence.to_formater seq.DOS.sequence a.alph
            | Heuristic_Selection seq ->
                    Sequence.to_formater seq.DOS.sequence a.alph
            (*| Relaxed_Lifted x ->
                    let b_states = RL.get_states x 1 in
                    let b_leftchild = RL.get_states x 2 in
                    let b_rightchild = RL.get_states x 3 in
                    let res = 
                        Array_ops.map_3 (fun a b c ->
                            Printf.sprintf "Cost: %d - left: %d - right: %d" a b
                            c) b_states b_leftchild b_rightchild 
                    in
                    String.concat "--" (Array.to_list res) *)
        in
        acc ^ code ^ ": " ^ seq ^ "; "
    in
    Array_ops.fold_right_2 builder "" a.codes a.characters 

let of_array spec sc code taxon =
    let c3 = spec.Data.tcm3d in
    let heur = make_default_heuristic ~c3 spec.Data.tcm2d in
    let create_item (x, _) =
        match spec.Data.initial_assignment, x with
        | `Partitioned clip, x ->
                Partitioned
                (PartitionedDOS.create (`Partitioned x) clip spec.Data.alph)
        | `AutoPartitioned (clip, size, table), [|x|] ->
                Partitioned 
                    (try
                        let partition = Hashtbl.find table taxon in
                        (PartitionedDOS.create (`AutoPartitioned
                        (partition, x)) clip 
                        spec.Data.alph)
                    with
                    | Not_found ->  
                            PartitionedDOS.empty clip size spec.Data.alph)
        | `GeneralNonAdd, [|x|] -> General_Prealigned (GenNonAdd.init_gnonadd_t x)
        | `DO, [|x|] -> Heuristic_Selection (DOS.create x)
        (*| `FS fs, [|x|] ->
                let res = RL.off_array fs taxon in
                Relaxed_Lifted res*)
        | `AutoPartitioned _, _ 
        | `GeneralNonAdd, _ 
        | `DO, _ -> assert false
        (*| `FS _, _ -> assert false*)

    in
    let codes = Array.map snd sc in
    let characters = Array.map create_item sc in
    let res = 
        { characters = characters; codes = codes; total_cost = 0.0; 
        alph = spec.Data.alph; code = code; heuristic = heur;
        priority = Array.to_list codes;} 
    in
    (*
    Status.user_message Status.Information (to_string res);
    *)
    res

let of_list spec lst code = of_array spec (Array.of_list lst) code

let same_codes a b =
    Array_ops.fold_right_2 (fun acc a b -> acc && (a = b)) true a.codes b.codes

(** [readjust ch1 ch2 par mine] returns a tuple [(a, b)], where [b] is the 
* set of sequences generated from (heuristically) readjusting [mine] to 
* somewhere in between [ch1], [ch2], and [par] (the two children and 
* parent of [mine] respectively, and [a] is the new cost of [b] as 
* parent of [ch1] and [ch2]. *)
let readjust mode to_adjust modified ch1 ch2 parent mine =
    assert (parent.alph = Alphabet.nucleotides);
    let use_ukk = match !Methods.algn_mode with
        | `Algn_Newkk  -> true
        | `Algn_Normal -> false
    in
    let new_modified = ref [] 
    and total_cost = ref 0 in
    let adjusted = 
        Array_ops.map_5 (fun code a b c d ->
            let skip_it =
                match to_adjust with
                | None -> false
                | Some to_adjust ->
                        not (All_sets.Integers.mem code to_adjust)
            in
            if skip_it then d
            else
                match a, b, c, d with
                | Heuristic_Selection a, Heuristic_Selection b, 
                    Heuristic_Selection c, Heuristic_Selection d ->
                        let changed, res, cost = 
                            DOS.readjust mode mine.heuristic a b c d use_ukk in
                        if changed then begin
                            new_modified := code :: !new_modified;
                            total_cost := cost + !total_cost;
                        end;
                        Heuristic_Selection res
                | _ -> assert false)
        mine.codes ch1.characters ch2.characters parent.characters
        mine.characters
    in
    let modified = 
        List.fold_left (fun acc x -> All_sets.Integers.add x acc)
        modified !new_modified
    in
    let total_cost = float_of_int !total_cost in
    modified, total_cost,
    { mine with characters = adjusted; total_cost = total_cost }

let to_single parent mine =
    let total_cost = ref 0 in
    let characters =
        Array_ops.map_2 (fun a b ->
            match a, b with
            | General_Prealigned a, General_Prealigned b ->
                    let res, c = GenNonAdd.to_single mine.alph mine.heuristic.c2 a b in
                    total_cost := c + !total_cost;
                    General_Prealigned res
            | Partitioned a, Partitioned b ->
                    let res, c = PartitionedDOS.to_single mine.heuristic a b in
                    total_cost := c + !total_cost;
                    Partitioned res
            | Heuristic_Selection a, Heuristic_Selection b ->
                    let res, c = DOS.to_single mine.heuristic a b in
                    total_cost := c + !total_cost;
                    Heuristic_Selection res
            | Partitioned _, _
            | _, Partitioned _
            | General_Prealigned _, _ 
            | Heuristic_Selection _, General_Prealigned _ -> assert false) parent.characters
                    mine.characters
    in
    let total_cost = float_of_int !total_cost in
    mine.total_cost, total_cost, 
    { mine with characters = characters; total_cost = total_cost }


let median code a b =
    let debug = false in
    if debug then Printf.printf "seqCS.median -> \n%!";
    let total_cost = ref 0 in
    let h = a.heuristic in
    let alph = a.alph in
    let use_ukk = match !Methods.algn_mode with
        | `Algn_Newkk  -> true
        | `Algn_Normal -> false
    in
    let characters =
        Array_ops.map_2 (fun a b ->
            match a, b with
            | General_Prealigned a, General_Prealigned b ->
                    let res, c = GenNonAdd.median h.c2 a b in
                    total_cost := c + !total_cost;
                    General_Prealigned res
            | Partitioned a, Partitioned b ->
                    let res, c = PartitionedDOS.median alph code h a b use_ukk in
                    total_cost := c + !total_cost;
                    Partitioned res
            | Heuristic_Selection a, Heuristic_Selection b ->
                    let res, c = DOS.median alph code h a b use_ukk in
                    total_cost := c + !total_cost;
                    Heuristic_Selection res
            (*| Relaxed_Lifted a, Relaxed_Lifted b ->
                    let res,tc = RL.median code a b in
                    total_cost := (int_of_float tc) + !total_cost;
                    Relaxed_Lifted res *)
            | Partitioned _, _
            | _, Partitioned _
            (*| Relaxed_Lifted _, _
            | _, Relaxed_Lifted _ *)
            | Heuristic_Selection _, _
            | General_Prealigned _ , _ -> assert false) a.characters b.characters
    in
    let res = { a with characters = characters; total_cost = float_of_int
    !total_cost } in
    if debug then Printf.printf "<- seqCS.median, tocal_cost<-%d\n%!" !total_cost;
    (*
    Status.user_message Status.Information (to_string res);
    *)
    res


let median_3 p n c1 c2 =
    let h = n.heuristic in
    let generic_map_4 f1 f3 a b c d use_ukk =
        Array_ops.map_4 (fun a b c d ->
            match a, b, c, d with
            | Heuristic_Selection a, Heuristic_Selection b, 
            Heuristic_Selection c, Heuristic_Selection d ->
                Heuristic_Selection (f1 h a b c d use_ukk)
            | General_Prealigned a, General_Prealigned b,
            General_Prealigned c, General_Prealigned d ->
                General_Prealigned (f3 h.c2 a b c d) 
            (*| Relaxed_Lifted a, Relaxed_Lifted b, 
                Relaxed_Lifted c, Relaxed_Lifted d ->
                    Relaxed_Lifted (f2 h a b c d) *)
            | Partitioned _, Partitioned _, Partitioned _, Partitioned _ ->
                    b
            | _ -> assert false)
        a b c d
    in
    (* A function to calculate the uppass values if the alphabet cannot 
    * handle the union of the items inside *)
    let median_no_union () = 
        let use_ukk = match !Methods.algn_mode with
            | `Algn_Newkk  -> true
            | `Algn_Normal -> false
        in
        generic_map_4 DOS.median_3_no_union GenNonAdd.median_3
                p.characters n.characters c1.characters c2.characters use_ukk
    in
    (* A function to calculate the uppass values if the alphabet does handle
    * properly the union of the items inside. *)
    let median_union () =
        let use_ukk = match !Methods.algn_mode with
            | `Algn_Newkk  -> true
            | `Algn_Normal -> false
        in
        generic_map_4 DOS.median_3_union GenNonAdd.median_3
                p.characters n.characters c1.characters c2.characters use_ukk
    in
    let characters = 
        let has_combinations = 1 = Cost_matrix.Two_D.combine n.heuristic.c2 in
        if has_combinations then median_union () else median_no_union ()
    in
    { n with characters = characters }

let distance missing_distance a b = 
    let missing_distance = int_of_float missing_distance in
    let h = a.heuristic in
    let alph = a.alph in
    let use_ukk = match !Methods.algn_mode with
        | `Algn_Newkk  -> true
        | `Algn_Normal -> false
    in
    float_of_int
        (Array_ops.fold_right_2
            (fun acc a b -> match a, b with
                | Partitioned a, Partitioned b ->
                    acc + (PartitionedDOS.distance alph h missing_distance a b use_ukk)
                | General_Prealigned a, General_Prealigned b ->
                    acc + (GenNonAdd.distance a b h.c2)
                | Heuristic_Selection a, Heuristic_Selection b ->
                    acc + (DOS.distance alph h missing_distance a b use_ukk) 
                (* | Relaxed_Lifted a, Relaxed_Lifted b ->
                    0  acc + (int_of_float (RL.distance a b))*)
                | Partitioned _, _
                (*| Relaxed_Lifted _, _*)
                | Heuristic_Selection _, _ 
                | General_Prealigned _, _ -> assert false)
            0 a.characters b.characters)

let dist_2 delta n a b =
    let h = n.heuristic in
    let delta = int_of_float delta in
    let use_ukk = match !Methods.algn_mode with
        | `Algn_Newkk  -> true
        | `Algn_Normal -> false
    in
    let x, deltaleft =
        Array_ops.fold_right_3 (fun (acc, deltaleft) n a b ->
            if deltaleft < 0 then (max_int / 10, deltaleft)
            else
                match n, a, b with
                | Partitioned n, Partitioned a, Partitioned b ->
                        let cost = PartitionedDOS.dist_2 h n a b use_ukk in
                        (acc + cost), (deltaleft - cost)
                | General_Prealigned n, General_Prealigned a, 
                General_Prealigned b ->
                    failwith "we don't have dist_2 for general prealigned data-type yet"
                | Heuristic_Selection n, 
                    Heuristic_Selection a, Heuristic_Selection b ->
                        let cost = DOS.dist_2 h n a b use_ukk in
                        (acc + cost), (deltaleft - cost)
                (*| Relaxed_Lifted n, Relaxed_Lifted a, Relaxed_Lifted b ->
                        let cost = RL.dist_2 n a b in
                        (acc + cost), (deltaleft - cost) *)
                | Partitioned _, _, _
                | _, Partitioned _, _
                | _, _, Partitioned _
                (*| Relaxed_Lifted _, _ , _
                | _, Relaxed_Lifted _ , _
                | _, _, Relaxed_Lifted _  *)
                | General_Prealigned _ , _, _
                | Heuristic_Selection _ , _ , _ -> assert false) (0, delta) 
                        n.characters a.characters b.characters
    in
    float_of_int x 


let f_codes memf s c = 
    let check x = memf x c in
    let codes = ref []
    and characters = ref [] in
    for i = (cardinal s) - 1 downto 0 do
        if check s.codes.(i) then begin
            codes := s.codes.(i) :: !codes;
            characters := s.characters.(i) :: !characters;
        end;
    done;
    { s with codes = Array.of_list !codes; characters = Array.of_list
    !characters }

let f_codes_comp s c = 
    f_codes (fun a b -> not (All_sets.Integers.mem a b)) s c

let f_codes s c =
    f_codes (fun a b -> (All_sets.Integers.mem a b)) s c

let compare_data a b =
    Array_ops.fold_right_2 (fun acc a b ->
        if acc <> 0 then acc
        else
            match a, b with
            | General_Prealigned a, General_Prealigned b ->
                    GenNonAdd.compare a b
            | Heuristic_Selection a, Heuristic_Selection b ->
                    DOS.compare a b
            | Partitioned a, Partitioned b ->
                    PartitionedDOS.compare a b
            | Partitioned _, _
            | _, Partitioned _
            | General_Prealigned _ , _ 
            | Heuristic_Selection _ , _ -> assert false)
    0 a.characters b.characters

let ( --> ) a b = b a 

let to_formatter report_type attr t do_to_single d : Xml.xml Sexpr.t list = 
    let use_ukk = match !Methods.algn_mode with
        | `Algn_Newkk  -> true
        | `Algn_Normal -> false
    in
    let h = t.heuristic in
    let res_state = ref (-1) in
    let is_fs = ref false in
    let rec output_sequence acc code seq do_to_single =
        let one_sequence (cmin, cmax, ccost, seqs) par seq =
            let cost = seq.DOS.costs in
            match par with
            | None -> 
                    (cmin +. cost.min), (cmax +. cost.max), 
                    ccost +. cost.max, seq.DOS.sequence :: seqs
            | Some par -> 
                    let par = par.DOS.sequence in
                    let s1, s2, min =
                        if use_ukk then 
                        Sequence.NewkkAlign.align_2 seq.DOS.sequence 
                        par h.c2 Sequence.NewkkAlign.default_ukkm
                        else
                            Sequence.Align.align_2 seq.DOS.sequence 
                        par h.c2 Matrix.default
                    in
                    let max = Sequence.Align.max_cost_2 s1 s2 h.c2 in
                    (cmin +. (float_of_int min)), 
                    (cmax +. (float_of_int max)),
                    ccost +. cost.max, seq.DOS.sequence :: seqs
        in
        let cost, costb, max, seq =
            match seq with
            | Partitioned seqs -> 
                    let min, max, total, seqs =
                        match do_to_single with
                        | None -> Array.fold_right (fun a acc ->
                                match a with
                                | PartitionedDOS.Last a 
                                | PartitionedDOS.First a
                                | PartitionedDOS.DO a -> 
                                        one_sequence acc None a) 
                                seqs (0.,0.,0.,[]) 
                        | Some (Partitioned par) ->
                                Array_ops.fold_right_2 (fun acc a b ->
                                    match a, b with
                                    | PartitionedDOS.Last a, 
                                        PartitionedDOS.Last b 
                                    | PartitionedDOS.First a, 
                                        PartitionedDOS.First b
                                    | PartitionedDOS.DO a, 
                                        PartitionedDOS.DO b -> 
                                            one_sequence acc (Some b) a
                                    | _ -> assert false)
                                (0., 0., 0., []) seqs par
                        | Some _ -> assert false
                    in
                    { min = min; max = max }, (`FloatFloatTuple (min, max)),
                    max, seqs
            | General_Prealigned seq ->
                    if debug then Printf.printf "seqCS.GP.to_formatter\n%!";
                    let cost = seq.GenNonAdd.costs in
                    let mincost,maxcost = 
                        cost.GenNonAdd.min,cost.GenNonAdd.max
                    in
                    let cost = { min = mincost;  max = maxcost } in
                    let costb,max = `FloatFloatTuple (mincost,maxcost), maxcost in 
                    cost, costb, max, [seq.GenNonAdd.seq]
            | Heuristic_Selection seq ->
                    if debug then  Printf.printf "seqCS.HS.to_formatter\n%!";
                    let cost = seq.DOS.costs in
                    let costb, max = 
                        match do_to_single with
                        | None -> 
                            `FloatFloatTuple (cost.min, cost.max), cost.max
                        | Some (Heuristic_Selection par) ->
                                let par = par.DOS.sequence in
                                let s1, s2, min =
                                    if use_ukk then
                                        Sequence.NewkkAlign.align_2 seq.DOS.sequence 
                                    par h.c2 Sequence.NewkkAlign.default_ukkm
                                    else 
                                        Sequence.Align.align_2 seq.DOS.sequence 
                                    par h.c2 Matrix.default
                                in
                                let max = 
                                    Sequence.Align.max_cost_2 s1 s2 h.c2 
                                in
                                `IntTuple (min, max), float_of_int max
                        | Some (General_Prealigned _)
                        | Some (Partitioned _) -> assert false 
                    in
                    cost, costb, max, [seq.DOS.sequence]
            (*| Relaxed_Lifted x  ->
                    is_fs := true;
                    if debug then  Printf.printf "seqCS.RL.to_formatter\n%!";
                    let bestc,best_state = RL.get_min_states x in
                    res_state := best_state;
                    DOS.make_cost bestc, 
                    `IntTuple (bestc, bestc),float_of_int bestc,
                    [x.RL.sequence_table.(best_state)] *)
        in
        let statename = (string_of_int !res_state) in
        let seq () = 
            match seq with
            | [] -> assert false
            | seq -> 
                    match report_type,!is_fs with
                    | `Normal,_ 
                    | `StateOnly,false ->
                    "state:"^statename^",seq:"^
                    (String.concat "#" 
                    (List.map (fun x -> 
                        let x = Sequence.del_first_char x in
                        Sequence.to_formater x t.alph) seq))
                    | `StateOnly,true -> statename
        in
        let definite_str = 
            if max > 0. then  `String "true"
            else `String "false"
        in 
        let module T = Xml.Characters in
        (PXML 
            -[T.sequence] 
                (* Attributes *)
                ([T.name] = [`String ((Data.code_character code
                d))])
                ([T.cost] = [costb])
                ([T.definite] = [definite_str])
                ([attr])
                
                { `Fun seq }
            --) :: acc
    in
    let parent = 
        match do_to_single with
        | None -> Array.map (fun _ -> None) t.characters
        | Some x -> Array.map (fun x -> Some x) x.characters
    in
    Array_ops.fold_right_3 output_sequence [] t.codes t.characters parent

let get_sequences (data:t) : Sequence.s array array =
    Array.map 
        (function
            | General_Prealigned x -> [| x.GenNonAdd.seq |]
            | Heuristic_Selection x -> [| x.DOS.sequence |]
            (*| Relaxed_Lifted rl -> rl.RL.sequence_table*)
            | Partitioned y ->
                Array.map 
                    (function | PartitionedDOS.DO x | PartitionedDOS.First x 
                              | PartitionedDOS.Last x -> x.DOS.sequence)
                    y)
        data.characters

let align_2 (one:t) (two:t) =
    let ones = get_sequences one and twos = get_sequences two in
    assert( (Array.length ones) = (Array.length twos));
    if 0 = Array.length ones then [| |]
    else begin
        Array.mapi (fun i _ -> 
            Array.mapi (fun j _ ->
                let a,b,m = 
(*                    if one.use_ukk then *)
(*                    Sequence.NewkkAlign.align_2 ones.(i).(j) twos.(i).(j) *)
(*                    one.heuristic.c2 Sequence.NewkkAlign.default_ukkm*)
(*                    else*)
                    Sequence.Align.align_2 ones.(i).(j) twos.(i).(j) one.heuristic.c2 Matrix.default
                in
                (a,b,m))
            ones.(i))
        ones
    end

let tabu_distance a = 
    Array.fold_left
        (fun sum -> function
            | Partitioned x         -> sum +. (PartitionedDOS.tabu_distance x)
            (*| Relaxed_Lifted _      -> sum*)
            | General_Prealigned y  -> (GenNonAdd.get_max_cost y.GenNonAdd.costs) +. sum
            | Heuristic_Selection y -> y.DOS.costs.max +. sum)
        0.0 a.characters


let explode cs = 
    let h = cs.heuristic in
    Array_ops.fold_right_2
        (fun acc code -> function
            | Partitioned _ -> failwith "TODO 12345"
            (*| Relaxed_Lifted _ -> failwith "TODO 12345"*)
            | General_Prealigned seq ->
                (code, seq.GenNonAdd.seq, h.c2, h.c3, cs.alph) :: acc
            | Heuristic_Selection seq ->
                (code, seq.DOS.sequence, h.c2, h.c3, cs.alph) :: acc)
        []
        cs.codes
        cs.characters


let classify_transformations leafa seqa leafb seqb acc =
    let cm = seqa.heuristic.c2 in
    let align_2_function = match !Methods.algn_mode with
        | `Algn_Newkk ->
            (fun x y -> Sequence.NewkkAlign.align_2 x y cm Sequence.NewkkAlign.default_ukkm)
        | `Algn_Normal ->
            (fun x y -> Sequence.Align.align_2 x y cm Matrix.default)
    and add_comb t v pair =
        let v =
            if All_sets.FullTupleMap.mem pair t
                then v +. (All_sets.FullTupleMap.find pair t)
                else v
        in
        All_sets.FullTupleMap.add pair v t
    and add_freq l e v f =
        if not l then f
        else begin
            let v =
                if All_sets.IntegerMap.mem e f
                    then (All_sets.IntegerMap.find e f) +. v
                    else v
            in
            All_sets.IntegerMap.add e v f
        end
    in
    let classify_sequence_pair useq1 useq2 acc =
        let aseq1, aseq2, _ = align_2_function useq1 useq2 in
        Sequence.foldi_2
            (fun acc i c1 c2 -> 
                let es1 = BitSet.Int.list_of_packed c1
                and es2 = BitSet.Int.list_of_packed c2 in
                let v1 = float_of_int (List.length es1)
                and v2 = float_of_int (List.length es2) in
                let v  = (v1 *. v2) in
                List.fold_left2
                    (fun (t,f) e1 e2 -> 
                        let t = add_comb t v (e1,e2) in
                        let f = f --> add_freq leafa e1 v1
                                  --> add_freq leafb e2 v2 in
                        t,f)
                    acc es1 es2)
            acc aseq1 aseq2
    in
    let a_rr = get_sequences seqb
    and b_rr = get_sequences seqa in
    Array_ops.fold_right_2
        (fun acc a_r b_r ->
            Array_ops.fold_right_2
                (fun acc a_seq b_seq -> classify_sequence_pair a_seq b_seq acc)
                acc a_r b_r)
        acc a_rr b_rr


let encoding enc x =
    Array.fold_left (fun acc x -> 
        match x with
        | Partitioned x -> acc +. (PartitionedDOS.encoding enc x)
        (*| Relaxed_Lifted _ -> acc*)
        | General_Prealigned x -> 
                acc +. (Sequence.encoding enc x.GenNonAdd.seq)
        | Heuristic_Selection x -> 
                acc +. (Sequence.encoding enc x.DOS.sequence)) 
    0.0 x.characters

module Kolmogorov = struct

    let correct_cost t m = 
        match m.Data.ks.Data.kolmo_spec.Data.mo with
        | Data.InDels _ 
        | Data.InDelSub _ 
        | Data.Subs _
        | Data.AffInDelSub _ ->
                { t with total_cost = t.total_cost /. Data.kolmo_round_factor }
        | Data.AffInDelAffSub _ -> failwith "Not implemented"
end
