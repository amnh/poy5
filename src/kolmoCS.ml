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

let () = SadmanOutput.register "KolmoCS" "$Revision: 1616 $"

(** The Kolmogorov complexity class of character *)

type t = {
    model : Data.kolmo_spec;
    matrix : Kolmo.Align.matrix;
    characters : SeqCS.t;
    branch_cost : float;
    leaf_cost : float;
}

type u = SeqCS.Union.u

let get_model set = 
    match set.model.Data.ks.Data.kolmo_spec.Data.mo with
    | Data.InDels _ -> Xml.KolSpecs.indelsonly
    | Data.InDelSub _ -> Xml.KolSpecs.indelsub
    | Data.Subs _ -> Xml.KolSpecs.subsonly
    | Data.AffInDelAffSub _ -> Xml.KolSpecs.affineindelsaffsub
    | Data.AffInDelSub _ -> Xml.KolSpecs.affineindelssub

let to_string set = 
    let model = get_model set in
    let rest = SeqCS.to_string set.characters in
    model ^ " - " ^ rest

let seqCS_median matrix a b =
    let total_cost = ref 0. in
    let gap = Cost_matrix.Two_D.gap matrix.Kolmo.Align.matrix in
    let characters =
        Array_ops.map_2 (fun a b ->
            match a, b with
            | SeqCS.Heuristic_Selection a, SeqCS.Heuristic_Selection b ->
                    let c, a_algn, b_algn, res = 
                        Kolmo.Align.align a.SeqCS.DOS.sequence 
                        b.SeqCS.DOS.sequence matrix 
                    in
                    let to_bitset a b = 
                        SeqCS.DOS.seq_to_bitset gap a (SeqCS.Raw b.SeqCS.DOS.sequence)
                    in
                    total_cost := c +. !total_cost;
                    let res = 
                        { SeqCS.DOS.sequence = res;
                        aligned_children = 
                            (to_bitset a_algn a), (to_bitset b_algn b), 
                            (SeqCS.Raw res);
                        costs = SeqCS.DOS.make_cost (int_of_float (ceil c));
                        position = 0; }
                    in
                    SeqCS.Heuristic_Selection res
            | SeqCS.Relaxed_Lifted _, SeqCS.Relaxed_Lifted _ 
            | SeqCS.Partitioned _, SeqCS.Partitioned _
            | SeqCS.Partitioned _, _
            | _, SeqCS.Partitioned _
            | SeqCS.Relaxed_Lifted _, _
            | _, SeqCS.Relaxed_Lifted _ -> assert false) 
        a.SeqCS.characters b.SeqCS.characters
    in
    let total_cost = !total_cost /. Data.kolmo_round_factor in
    let res = { a with SeqCS.characters = characters; total_cost = total_cost} in
    res, total_cost

let median code a b = 
    { a with characters = fst (seqCS_median a.matrix a.characters b.characters) }

let total_cost a = a.characters.SeqCS.total_cost 

let median_3 a b c d =
    (* We don't do anything here *)
    b

let distance c d =
    snd (seqCS_median c.matrix c.characters d.characters)
    
let dist_2 a b c d = distance c d

let to_formatter ref_codes attr t d =
    let attr = (Xml.KolSpecs.model, `String (get_model t)) :: attr in
    SeqCS.to_formatter attr t.characters None d

let f_codes a b = 
    { a with characters = SeqCS.f_codes a.characters b } 

let f_codes_comp a b =
    { a with characters = SeqCS.f_codes_comp a.characters b }

let cardinal a = SeqCS.cardinal a.characters

let root_cost x = 
    x.model.Data.ks.Data.kolmo_spec.Data.root_cost +.
    (SeqCS.encoding x.model.Data.ks.Data.kolmo_spec.Data.be x.characters)
    +. (total_cost x)

let of_array spec code taxon num_taxa =
    let c = SeqCS.of_array spec.Data.dhs code taxon num_taxa in
    let prob = spec.Data.ks.Data.kolmo_spec.Data.event_prob in
    let log2 x = (log x) /. (log 2.) in
    let first = ~-. (log2 prob)
    and second = ~-. (log2 (prob *. (1. -. prob))) in
    let diff = second -. first in
    let initial = first -. diff in
    Printf.printf "Diff is %f and initial is %f\n%!" diff initial;
    let matrix = { Kolmo.Align.event_cost =  diff *. Data.kolmo_round_factor;
    first_event_cost = initial *. Data.kolmo_round_factor;
        matrix = c.SeqCS.heuristic.SeqCS.c2; } in
    let other_spec = spec.Data.ks.Data.kolmo_spec in
    let bc = other_spec.Data.branch_cost 
    and lc = other_spec.Data.leaf_cost in
    { model = spec; characters = c; branch_cost = bc; leaf_cost = lc; matrix =
        matrix }

let to_union x = SeqCS.to_union x.characters

let get_sequence_union code x = SeqCS.Union.get_sequence_union code x

let union a b c = SeqCS.Union.union a.characters b c

let compare_data a b = SeqCS.compare_data a.characters b.characters

let tabu_distance a b = 0.

let distance_union a b =
    (SeqCS.Union.distance_union a b) /. Data.kolmo_round_factor

let get_dynamic_preliminary d = DynamicCS.SeqCS d.characters

let to_single a b c d = 
    let a, _, r = SeqCS.to_single c.characters d.characters in
    let r = { d with characters = r } in
    let cost = distance c r in
    a, cost, r
