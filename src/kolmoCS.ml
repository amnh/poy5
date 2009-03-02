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
    characters : DynamicCS.t;
    branch_cost : float;
    leaf_cost : float;
}

type u = DynamicCS.u

let get_model set = 
    match set.model.Data.ks.Data.kolmo_spec.Data.mo with
    | Data.InDels _ -> Xml.KolSpecs.indelsonly
    | Data.InDelSub _ -> Xml.KolSpecs.indelsub
    | Data.Subs _ -> Xml.KolSpecs.subsonly
    | Data.AffInDelAffSub _ -> Xml.KolSpecs.affineindelsaffsub
    | Data.AffInDelSub _ -> Xml.KolSpecs.affineindelssub

let to_string set = 
    let model = get_model set in
    let rest = DynamicCS.to_string set.characters in
    model ^ " - " ^ rest

let median code a b = 
    let res = DynamicCS.median code a.characters b.characters in
    let res = DynamicCS.Kolmogorov.correct_cost res b.model in
    { a with characters = res }

let total_cost a = DynamicCS.total_cost a.characters

let median_3 a b c d = 
    let m =
        DynamicCS.median_3 a.characters b.characters c.characters d.characters
    in
    let m = DynamicCS.Kolmogorov.correct_cost m c.model in
    { b with characters = m }

let distance c d =
    let res = DynamicCS.median (-1) c.characters d.characters in
    let res = DynamicCS.Kolmogorov.correct_cost res c.model in
    (DynamicCS.total_cost res) 
    
let dist_2 a b c d = distance c d

let to_formatter ref_codes attr t d =
    let attr = (Xml.KolSpecs.model, `String (get_model t)) :: attr in
    DynamicCS.to_formatter ref_codes attr  t.characters None d

let f_codes a b = 
    { a with characters = DynamicCS.f_codes a.characters b } 

let f_codes_comp a b =
    { a with characters = DynamicCS.f_codes_comp a.characters b }

let cardinal a = DynamicCS.cardinal a.characters

let root_cost x = 
    x.model.Data.ks.Data.kolmo_spec.Data.root_cost +.
    (DynamicCS.encoding x.model.Data.ks.Data.kolmo_spec.Data.be x.characters)
    +. (total_cost x)

let of_array spec c code taxon num_taxa =
    let c = DynamicCS.of_array spec.Data.dhs c code taxon num_taxa in
    let other_spec = spec.Data.ks.Data.kolmo_spec in
    let bc = other_spec.Data.branch_cost 
    and lc = other_spec.Data.leaf_cost in
    { model = spec; characters = c; branch_cost = bc; leaf_cost = lc; }

let to_union x = DynamicCS.to_union x.characters

let get_sequence_union code x = DynamicCS.get_sequence_union code x

let union a b c = DynamicCS.union a.characters b c

let compare_data a b = DynamicCS.compare_data a.characters b.characters

let tabu_distance a b = 
    let d = DynamicCS.tabu_distance a.characters b.characters in
    d /. Data.kolmo_round_factor

let distance_union a b =
    (DynamicCS.distance_union a b) /. Data.kolmo_round_factor

let get_dynamic_preliminary d = d.characters

let to_single a b c d = 
    let b = 
        match b with
        | None -> None
        | Some x -> Some x.characters
    in
    let a, _, r = DynamicCS.to_single a b c.characters d.characters in
    let r = { d with characters = r } in
    let cost = distance c r in
    a, cost, r
