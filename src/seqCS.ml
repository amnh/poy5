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

(** A Sequence Character Set implementation *)
exception Illegal_Arguments
let () = SadmanOutput.register "SeqCS" "$Revision: 2174 $"


module Codes = All_sets.IntegerMap

type cost_tuple = {
    min : float;
    max : float;
}

(** A sequence character type. *)
type t = { 
    (** The sequences that form the set *)
    sequences : Sequence.s Codes.t;
    aligned_children : (Sequence.s * Sequence.s) Codes.t;
    costs : cost_tuple Codes.t; (** The cost of each median *)
    total_cost : float;             (** The total cost of the character set *)
    c2 : Cost_matrix.Two_D.m;       (** The two dimensional cost matrix to 
                                    be used in the character set *)
    c3 : Cost_matrix.Three_D.m;     (** The three dimensional cost matrix to be 
                                    used in the character set *)
    alph : Alphabet.a;              (** The alphabet of the sequence set *)
    code : int;                     (** The set code *)
    pool : Sequence.Pool.p;         (** The pool of sequences these are produced
                                    from *)
    priority : int list;            (** The information ordering *)
}

(* The union of sequences *)
type ustr = {
    unions : Sequence.Unions.u Codes.t;
    u_c2 : Cost_matrix.Two_D.m;
    u_alph : Alphabet.a;
}

type u = ustr option

let cardinal x = Codes.fold (fun _ _ acc -> acc + 1) x.sequences 0

let reprioritize a b =
    { b with priority = a.priority }

let compare_union a b =
    match a, b with
    | Some a, Some b ->
            let single_compare x y =
                let z = Codes.find x b.unions in
                Sequence.Unions.compare y z
            in
            Codes.fold (fun x y acc ->
                match acc with
                | 0 -> single_compare x y 
                | x -> x)
            a.unions
            0
    | None, None -> 0
    | Some _, None -> -1
    | None, Some _ -> 1

let prioritize a = 
    let lst = Codes.fold (fun x y acc -> (x, y) :: acc) a.costs [] in
    let lst = List.sort (fun (_, a) (_, b) -> compare b a) lst in
    let x, _ = List.split lst in
    { a with priority = x }

let empty code c2 alph pool = 
    let c3 = Cost_matrix.Three_D.of_two_dim c2 in
    let emp = Codes.empty in
    {
        sequences = emp;
        aligned_children = emp;
        costs = emp;
        total_cost = 0.0;
        c2 = c2;
        c3 = c3;
        alph = alph;
        code = code;
        pool = pool;
        priority = [];
    }

let to_union a = 
    if a.alph = Alphabet.nucleotides then
        let folder code seq acc = 
            Codes.add code (Sequence.Unions.leaf seq) acc 
        in
        let new_unions = Codes.fold folder a.sequences Codes.empty in
        Some { unions = new_unions; u_c2 = a.c2; u_alph = a.alph }
    else None

let add t (seq, code) =
    let ns = Codes.add code seq t.sequences in
    let na = Codes.add code (seq, seq) t.aligned_children in
    { t with sequences = ns; aligned_children = na;
    priority = code :: t.priority }

let of_array spec sc code = 
    let no_cost = { min = 0.0; max = 0.0 } in
    let adder (x, y, z, acc) (a, b) = 
        let a = Sequence.clone_pool spec.Data.pool a in
        (Codes.add b a x), (Codes.add b no_cost y), (Codes.add b (a, a) z), 
        (b :: acc)
    in
    let empty = Codes.empty in
    let seqs, costs, algn_chld, priority = 
        Array.fold_left adder (empty, empty, empty, []) sc 
    in
    {
        sequences = seqs;
        aligned_children = algn_chld;
        costs = costs;
        total_cost = 0.0;
        c2 = spec.Data.tcm2d;
        c3 = spec.Data.tcm3d;
        alph = spec.Data.alph;
        code = code;
        pool = spec.Data.pool;
        priority = priority;
    }

let of_list spec lst code = 
    let arr = Array.of_list lst in
    of_array spec arr code 

let to_list t =
    Codes.fold (fun code seq acc -> (seq, code) :: acc) t.sequences []

let same_codes a b =
    let checker x _ res = res && (Codes.mem x b) in
    Codes.fold checker a true

let sequence_median c2 m gap code seqa seqb (res_medians, res_algn_chld, res_costs, total) =
    let no_cost = { min = 0.0; max = 0.0 } in
    if Sequence.is_empty seqa gap then
        let new_costs = Codes.add code no_cost res_costs 
        and new_median = Codes.add code seqb res_medians 
        and new_algn_chld = Codes.add code (seqb, seqb) res_algn_chld in
        new_median, new_algn_chld, new_costs, total
    else if Sequence.is_empty seqb gap then
        let new_costs = Codes.add code no_cost res_costs 
        and new_median = Codes.add code seqa res_medians 
        and new_algn_chld = Codes.add code (seqa, seqa) res_algn_chld in
        new_median, new_algn_chld, new_costs, total
    else 
        let tmpa, tmpb, tmpcost = Sequence.Align.align_2 seqa seqb c2 m in
        let seqm = Sequence.Align.ancestor_2 tmpa tmpb c2 
        and rescost = 
            let tmp_cost = float_of_int tmpcost in
            { min = tmp_cost; max = tmp_cost } 
        in
        let new_median = Codes.add code seqm res_medians
        and new_costs = Codes.add code rescost res_costs 
        and new_algn_chld = Codes.add code (tmpa, tmpb) res_algn_chld
        and new_total = total + tmpcost in
        new_median, new_algn_chld, new_costs, new_total

let cardinal_union ua =
    match ua with
    | None -> 0
    | Some ua ->
            Codes.fold (fun _ x acc -> 
                Sequence.length x.Sequence.Unions.seq + acc) ua.unions 0

let poly_saturation x v =
    match x with
    | None -> 0.0
    | Some x ->
            let sat, len =
                Codes.fold (fun _ x (acc, len) ->
                    let nlen = Sequence.length x.Sequence.Unions.seq 
                    and sat = Sequence.poly_saturation x.Sequence.Unions.seq v in
                    acc +. (sat *. (float_of_int nlen)), len + nlen) 
                x.unions (0.0, 0)
            in
            sat /. (float_of_int len)

let union self ua ub =
    match ua, ub with
    | Some ua, Some ub ->
            let empty = Codes.empty 
            and c2 = ua.u_c2 in
            let gap = Cost_matrix.Two_D.gap c2 in
            let union code (tmpa, tmpb) res_union =
                let uniona = Codes.find code ua.unions
                and unionb = Codes.find code ub.unions in
                if Sequence.is_empty uniona.Sequence.Unions.seq gap then
                    Codes.add code unionb res_union
                else if Sequence.is_empty unionb.Sequence.Unions.seq gap then
                    Codes.add code uniona res_union
                else 
                    let union = Sequence.Unions.union tmpa tmpb uniona unionb c2 in
                    Codes.add code union res_union
            in
            let union = Codes.fold union self.aligned_children empty in
            Some { ua with unions = union }
    | None, None -> None
    | _, _ -> failwith "SeqCS.union"

let to_single_root parent =
    let empty = Codes.empty 
    and no_cost = { min = 0.0; max = 0.0 }
    and c2 = parent.c2 in
    let median code parent_seq acc =
        let (res_medians, res_costs) = acc in
        let seqm = Sequence.select_one parent_seq c2 in
        let rescost = no_cost in
        let new_single = Codes.add code seqm res_medians
        and new_costs = Codes.add code rescost res_costs in
        new_single, new_costs
    in
    let sequences, costs = 
        Codes.fold median parent.sequences (empty, empty)
    in
    { parent with sequences = sequences; costs = costs; total_cost = 0.0 }

(** [readjust ch1 ch2 par mine] returns a tuple [(a, b)], where [b] is the set
* of sequences generated from (heuristically) readjusting [mine] to somewhere in
* between [ch1], [ch2], and [par] (the two children and parent of [mine]
* respectively, and [a] is the new cost of [b] as parent of [ch1] and [ch2]. *)
let readjust to_adjust modified ch1 ch2 parent mine =
    assert (parent.alph = Alphabet.nucleotides);
    let empty = Codes.empty 
    and no_cost = { min = 0.0; max = 0.0 }
    and c2 = parent.c2 
    and c3 = parent.c3 in
    let gap = Cost_matrix.Two_D.gap c2 in
    let adjusted code parent_seq acc =
        let to_adjust =
            match to_adjust with
            | None -> All_sets.Integers.singleton code
            | Some x -> x
        in
        let (modified, res_medians, res_costs, total) = acc in
        let my_sequence = Codes.find code mine.sequences 
        and ch1_sequence = Codes.find code ch1.sequences
        and ch2_sequence = Codes.find code ch2.sequences in
        if (not (All_sets.Integers.mem code to_adjust)) then 
            let new_costs = Codes.add code no_cost res_costs 
            and new_single = Codes.add code my_sequence res_medians in
            modified, new_single, new_costs, total
        else if Sequence.is_empty ch1_sequence gap then 
            let new_costs = Codes.add code no_cost res_costs 
            and new_single = Codes.add code ch2_sequence res_medians in
            let modified = 
                if 0 <> compare my_sequence ch2_sequence then
                    All_sets.Integers.add code modified
                else modified
            in
            modified, new_single, new_costs, total
        else if Sequence.is_empty ch2_sequence gap then 
            let new_costs = Codes.add code no_cost res_costs 
            and new_single = Codes.add code ch1_sequence res_medians in
            let modified = 
                if 0 <> compare my_sequence ch2_sequence then
                    All_sets.Integers.add code modified
                else modified
            in
            modified, new_single, new_costs, total
        else 
            let tmpcost, seqm, changed = 
                Sequence.Align.readjust_3d ch1_sequence ch2_sequence my_sequence
                c2 c3 parent_seq
            in
            let rescost = 
                let tmpcost = float_of_int tmpcost in
                { min = tmpcost; max = tmpcost } 
            in
            let new_single = Codes.add code seqm res_medians
            and new_costs = Codes.add code rescost res_costs 
            and new_total = total + tmpcost in
            let modified = 
                if changed then All_sets.Integers.add code modified
                else modified
            in
            modified, new_single, new_costs, new_total
    in
    let modified, sequences, costs, total_cost = 
        Codes.fold adjusted parent.sequences (modified, empty, empty, 0)
    in
    let tc = float_of_int total_cost in
    modified,
    tc,
    { mine with sequences = sequences; costs = costs; total_cost = tc }

let to_single parent mine =
    let previous_cost = mine.total_cost in
    let empty = Codes.empty 
    and no_cost = { min = 0.0; max = 0.0 }
    and c2 = parent.c2
    and m = Matrix.default in
    let gap = Cost_matrix.Two_D.gap c2 in
    let median code parent_seq acc =
        let (res_medians, res_costs, total) = acc in
        let my_sequence = Codes.find code mine.sequences in
        if Sequence.is_empty parent_seq gap then
            let new_costs = Codes.add code no_cost res_costs 
            and new_single = Codes.add code my_sequence res_medians in
            new_single, new_costs, total
        else if Sequence.is_empty my_sequence gap then
            let new_costs = Codes.add code no_cost res_costs 
            and new_single = Codes.add code parent_seq res_medians in
            new_single, new_costs, total
        else begin 
            let seqm, tmpcost = 
                Sequence.Align.closest parent_seq my_sequence c2 m 
            in
            let rescost = 
                let tmpcost = float_of_int tmpcost in
                { min = tmpcost; max = tmpcost } 
            in
            let new_single = Codes.add code seqm res_medians
            and new_costs = Codes.add code rescost res_costs 
            and new_total = total + tmpcost in
            new_single, new_costs, new_total
        end 
    in
    let sequences, costs, total_cost = 
        Codes.fold median parent.sequences (empty, empty, 0)
    in
    let tc = float_of_int total_cost in
    previous_cost, tc, { mine with sequences = sequences; costs = costs; total_cost = tc }

let median a b =
    let empty = Codes.empty 
    and no_cost = { min = 0.0; max = 0.0 }
    and c2 = a.c2
    and m = Matrix.default in
    let gap = Cost_matrix.Two_D.gap c2 in
    let median code seqa acc =
        let (res_medians, res_algn_chld, res_costs, total) = acc in
        let seqb = Codes.find code b.sequences  in
        if Sequence.is_empty seqa gap then
            let new_costs = Codes.add code no_cost res_costs 
            and new_median = Codes.add code seqb res_medians 
            and new_algn_chld = Codes.add code (seqb, seqb) res_algn_chld in
            new_median, new_algn_chld, new_costs, total
        else if Sequence.is_empty seqb gap then
            let new_costs = Codes.add code no_cost res_costs 
            and new_median = Codes.add code seqa res_medians 
            and new_algn_chld = Codes.add code (seqa, seqa) res_algn_chld in
            new_median, new_algn_chld, new_costs, total
        else 
            let tmpa, tmpb, tmpcost = Sequence.Align.align_2 seqa seqb c2 m in
            let seqm = Sequence.Align.ancestor_2 tmpa tmpb c2 
            and rescost = 
                let tmp_cost = float_of_int tmpcost in
                { min = tmp_cost; max = tmp_cost } 
            in
            let new_median = Codes.add code seqm res_medians
            and new_costs = Codes.add code rescost res_costs 
            and new_algn_chld = Codes.add code (tmpa, tmpb) res_algn_chld
            and new_total = total + tmpcost in
            new_median, new_algn_chld, new_costs, new_total
    in
    let sequences, algn_chld, costs, total_cost = 
        Codes.fold median a.sequences (empty, empty, empty, 0)
    in
    { a with sequences = sequences; costs = costs; total_cost = float_of_int
    total_cost; aligned_children = algn_chld; }

let median_3 p n c1 c2 =
    let gap = Cost_matrix.Two_D.gap p.c2 in
    (* A function to calculate the uppass values if the alphabet cannot handle
    * the union of the items inside *)
    let median_no_union code sp (res_medians, costs) = 
        let get_all s =
            let s1, s2, costs = 
                Sequence.Align.align_2 sp s p.c2 Matrix.default
            in
            Sequence.Align.median_2 s1 s2 p.c2, costs,
            Sequence.Align.max_cost_2 s1 s2 p.c2
        in
        let res, cost = 
            let a = Codes.find code c1.sequences
            and b = Codes.find code c2.sequences in
            let res, minc, maxc = 
                let (_, costa, _) as resa = get_all a
                and (_, costb, _) as resb = get_all b in
                if costa < costb then resa
                else resb
            in
            (if gap = Sequence.get res 0 then res
            else let _ = Sequence.prepend res gap in res), 
            { min = float_of_int minc; 
            max = float_of_int maxc }
        in
        (Codes.add code res res_medians), 
        (Codes.add code cost costs)
    in
    (* A function to calculate the uppass values if the alphabet does handle
    * properly the union of the items inside. *)
    let median_union code sp (res_medians, costs) = 
        let res,cost = 
            let a, b = Codes.find code n.aligned_children in
            assert (Sequence.length a = Sequence.length b);
            let res =
                Sequence.Align.union a b in
            let a, b, cost = 
                Sequence.Align.align_2 sp res p.c2 Matrix.default 
            in
            let res = Sequence.Align.median_2 a b p.c2 in
            (if gap = Sequence.get res 0 then res
            else let _ = Sequence.prepend res gap in res), 
            { min = float_of_int cost; 
            max = float_of_int (Sequence.Align.max_cost_2 a b p.c2) }
        in
        (Codes.add code res res_medians), 
        (Codes.add code cost costs)
    in
    let acc = Codes.empty in
    let seqs, costs = 
        let has_combinations = 1 = Cost_matrix.Two_D.combine p.c2 in
        if has_combinations then
            Codes.fold median_union p.sequences (acc, acc) 
        else Codes.fold median_no_union p.sequences (acc, acc)
    in
    { n with sequences = seqs; costs = costs }

let distance a b = 
    let gap = Cost_matrix.Two_D.gap a.c2 in
    let single_distance code seqa acc =
        let seqb = Codes.find code b.sequences in
        if Sequence.is_empty seqa gap || Sequence.is_empty seqb gap then
            acc
        else
            let cost =
IFDEF USE_VERIFY_COSTS THEN
                let seqa, seqb, cost = 
                    Sequence.Align.align_2 seqa seqb a.c2 Matrix.default 
                in
                let () = 
                    assert (
                        let real_cost = Sequence.Align.verify_cost_2 cost seqa seqb
                        a.c2 in
                        if cost < real_cost then
                            let () = 
                                Printf.printf "Failed alignment between
                                \n%s\nand\n%s\nwith claimed cost %d and real cost %d\n%!" 
                                (Sequence.to_string seqa a.alph)
                                (Sequence.to_string seqb b.alph)
                                cost
                                real_cost
                            in
                            false
                        else true
                ) 
                in
                cost
ELSE 
                let deltaw = 
                    let tmp = (max (Sequence.length seqa) (Sequence.length seqb)) - (min
                    (Sequence.length seqa) (Sequence.length seqb)) in
                    if tmp > 8 then tmp 
                    else 8
                in
                (Sequence.Align.cost_2 ~deltaw seqa seqb a.c2 Matrix.default)
END
            in
            acc + cost
    in
    float_of_int (Codes.fold single_distance a.sequences 0)

let distance_union a b = 
    match a, b with
    | Some a, Some b ->
            let gap = Cost_matrix.Two_D.gap a.u_c2 in
            let single_distance code seqa acc =
                let seqa = seqa.Sequence.Unions.seq in
                let seqb =
                    let seqb = Codes.find code b.unions in
                    seqb.Sequence.Unions.seq 
                in
                if Sequence.is_empty seqa gap || Sequence.is_empty seqb gap then
                    acc
                else
                    let deltaw = 
                        let tmp = 
                            (max (Sequence.length seqa) (Sequence.length seqb)) -
                            (min (Sequence.length seqa) (Sequence.length seqb)) 
                        in
                        if tmp > 8 then tmp 
                        else 8
                    in
                    acc + 
                    (Sequence.Align.cost_2 ~deltaw:deltaw seqa seqb a.u_c2 Matrix.default)
            in
            float_of_int (Codes.fold single_distance a.unions 0)
    | None, None -> 0.0
    | Some _, _ 
    | None, _ -> failwith "Impossible?"

let to_string a =
    let builder code seq acc =
        let code = string_of_int code 
        and seq = Sequence.to_formater seq a.alph in
        acc ^ code ^ ": " ^ seq ^ "; "
    in
    Codes.fold builder a.sequences ""

let dist_2 delta n a b =
    let delta = int_of_float delta 
    and gap = Cost_matrix.Two_D.gap n.c2 in
    let cost_calculator (acc, deltaleft) code =
        let seqb = Codes.find code a.sequences in
        if deltaleft < (-1) then begin
            Status.user_message Status.Information "Breaking early";
            (max_int, deltaleft)
        end else begin
            let seqn = Codes.find code n.sequences
            and seqa = Codes.find code a.sequences in
            if Sequence.is_empty seqn gap then acc, deltaleft
            else 
                let tmp = 
                    if Sequence.is_empty seqa gap then seqn
                    else 
                        Sequence.Align.full_median_2 seqa seqb b.c2 
                        Matrix.default
                in
                let cost= Sequence.Align.cost_2 seqn tmp b.c2 Matrix.default in
                acc + cost, deltaleft - cost
        end
    in
    let x, deltaleft = List.fold_left cost_calculator (0, delta) b.priority in
    float_of_int x 

let f_codes s c = 
    let check x = All_sets.Integers.mem x c in
    let adder x y acc = 
        if check x then Codes.add x y acc 
        else acc
    in
    let n_seq_semap = Codes.fold adder s.sequences Codes.empty
    and n_seq_costs = Codes.fold adder s.costs Codes.empty in
    { s with sequences = n_seq_semap; costs = n_seq_costs }

let f_codes_comp s c = 
    let check x = not (All_sets.Integers.mem x c) in
    let adder x y acc = 
        if check x then Codes.add x y acc 
        else acc
    in
    let n_seq_semap = Codes.fold adder s.sequences Codes.empty
    and n_seq_costs = Codes.fold adder s.costs Codes.empty in
    { s with sequences = n_seq_semap; costs = n_seq_costs }

let compare_data a b =
    let comparator_seqs code seqb acc =
        if acc = 0 then begin
            let seqa = Codes.find code a.sequences in
            Sequence.compare seqb seqa 
        end else acc
    in
    Codes.fold comparator_seqs b.sequences 0 

let ( --> ) a b = b a 

let to_formatter attr t to_single d : Tags.output list = 
    let output_sequence code seq acc =
        let cost = Codes.find code t.costs in
        let costb, max = 
            match to_single with
            | None -> 
                    (string_of_float cost.min) ^ " - " ^ (string_of_float
                    cost.max), cost.max
            | Some par ->
                    let par = Codes.find code par.sequences in
                    let s1, s2, min = 
                        Sequence.Align.align_2 seq par t.c2 Matrix.default
                    in
                    let max = Sequence.Align.max_cost_2 s1 s2 t.c2 in
                    (string_of_int min) ^ " - " ^ (string_of_int max), 
                    float_of_int max
        in
        let definite_str = 
            if max > 0. then  "true"
            else "false"
        in 
        let seq = Sequence.to_formater seq t.alph in
        let attributes = 
            (Tags.Characters.name, (Data.code_character code d)) ::
                (Tags.Characters.cost, costb) ::
                (Tags.Characters.definite, definite_str) :: attr
        in
        let contents = `String seq in
        (Tags.Characters.sequence, attributes, contents) :: acc
    in
    Codes.fold output_sequence t.sequences []

let tabu_distance a = 
    All_sets.IntegerMap.fold (fun _ y sum -> y.max +. sum) a.costs 0.0

let explode cs = 
    All_sets.IntegerMap.fold 
    (fun code seq acc ->
        (code, seq, cs.c2, cs.c3, cs.alph) :: acc)
    cs.sequences
    []

let get_sequence_union code x =
    match x with
    | None -> None
    | Some x ->
            try Some (Codes.find code x.unions) with
            | _ -> None

let encoding enc x =
    Codes.fold (fun _ x acc -> acc +. (Sequence.encoding enc x)) x.sequences 0.0

module Kolmogorov = struct

    let create_function distr =
        let init = distr.Data.selfp in
        match distr.Data.distr with
        | Data.Items (IntSpec.Constant arr) ->
                fun len -> init +. ((float_of_int len) *. arr.(0))
        | Data.Items (IntSpec.Variable arr) ->
                fun len -> init +. arr.(len)
        | Data.MaxLength (len, IntSpec.Constant arr) ->
                fun clen ->
                    let c = float_of_int (clen / len)
                    and r = clen mod len in
                    (init *. (c +. 1.)) +. (arr.(0) *. c) +. ((float_of_int r) *.
                    arr.(0))
        | Data.MaxLength (len, IntSpec.Variable arr) ->
                fun clen ->
                    let c = float_of_int (clen / len)
                    and r = clen mod len in
                    (init *. (c +. 1.)) +. (arr.(len) *. c) +. arr.(r)

    let correct_cost t m = 
        let tcm = m.Data.dhs.Data.tcm2d in
        match m.Data.ks.Data.kolmo_spec.Data.mo with
        | Data.InDels _ 
        | Data.InDelSub _ 
        | Data.Subs _ -> 
                { t with total_cost = t.total_cost /. Data.kolmo_round_factor }
        | Data.AffInDelSub (ins, del, sub) ->
                let ins = create_function ins
                and del = create_function del 
                and gap = Cost_matrix.Two_D.gap tcm in
                let pair_cost _ (s1, s2) acc =
                    let len = Sequence.length s1 
                    and i = ref 0 in
                    while !i < len do
                        let b1 = Sequence.get s1 !i
                        and b2 = Sequence.get s2 !i in
                        if b1 <> gap && b2 <> gap then begin
                            acc := 
                                !acc +. ((float_of_int (Cost_matrix.Two_D.cost
                                b1 b2 tcm)) /. Data.kolmo_round_factor);
                            incr i;
                        end else if b1 = gap && b2 <> gap then
                            let total = 
                                let len = 
                                    incr i;
                                    ref 1 
                                in
                                while (!i < !len) && (gap <> Sequence.get s1 !i) do
                                    incr i;
                                    incr len;
                                done;
                                !len
                            in
                            acc := !acc +. (ins total)
                        else if b2 = gap && b1 <> gap then
                            let total = 
                                let len = 
                                    incr i;
                                    ref 0 
                                in
                                while (!i < !len) && (gap <> Sequence.get s2 !i) do
                                    let b2 = Sequence.get s2 !i in
                                    incr i;
                                    incr len;
                                    acc :=
                                        !acc +. ((float_of_int
                                        (Cost_matrix.Two_D.cost gap b2 tcm)) /.
                                        Data.kolmo_round_factor);
                                done;
                                !len
                            in
                            acc := !acc +. (del total)
                        else incr i
                    done;
                    acc
                in
                let cost = Codes.fold pair_cost t.aligned_children (ref 0.0) in
                { t with total_cost = !cost }
        | _ -> failwith "Still programming it"
end
