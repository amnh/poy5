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

let () = SadmanOutput.register "SearchInformation" "$Revision: 2047 $"

let handle_tree_information trees acc = function
    | `Number -> 
            let counter = ref 0 in
            Sexpr.leaf_iter (fun _ -> incr counter) trees;
            acc ^ ("@[<hov>Number:@ " ^ string_of_int !counter ^ "@]@,")
    | `Minimum ->
            let min = 
                let minimizer = fun acc x ->
                    let x = Ptree.get_cost `Adjusted x in
                    if x < acc then x else acc
                in
                Sexpr.fold_left minimizer infinity trees 
            in
            acc ^ ("@[<hov>Minimum@ lenght:@ " ^ string_of_float min ^ "@]@,")
    | `Maximum ->
            let max = 
                let maximizer = fun acc x ->
                    let x = Ptree.get_cost `Adjusted x in
                    if x > acc then x else acc
                in
                Sexpr.fold_left maximizer 0.0 trees 
            in
            acc ^ ("@[<hov>Maximum@ length:@ " ^ string_of_float max ^ "@]@,")
    | `Summary ->
            let counter = ref 0 in
            Sexpr.leaf_iter (fun _ -> incr counter) trees;
            let min, cnt = 
                let minimizer = fun (acc, cnt) x ->
                    let x = Ptree.get_cost `Adjusted x in
                    if x < acc then (x, 1) 
                    else if x = acc then (x, cnt + 1)
                    else (acc, cnt)
                in
                Sexpr.fold_left minimizer (infinity, 0) trees 
            in
            let max = 
                let maximizer = fun acc x ->
                    let x = Ptree.get_cost `Adjusted x in
                    if x > acc then x else acc
                in
                Sexpr.fold_left maximizer (-. 100000000.) trees 
            in
            let storing =
                acc ^ "@[<hov>Storing@ " ^ string_of_int !counter
                ^ "@ tree" in
            (match !counter with
            | 0 -> storing ^ "s@]"
            | 1 -> storing ^ "@ with@ cost@ " ^ string_of_float min ^ "@]"
            | _ -> storing ^ "s@ with@ costs@ "
                  ^ string_of_float min ^ " to " ^ string_of_float max
                    ^ "@]" ^ 
                    if cnt > 0 then "@,@[<hov>Best@ cost@ was@ found@ " ^
                        if cnt > 1 then 
                            string_of_int cnt ^ "@ times@]" 
                        else "once@]" 
                    else "")

let append acc h =  acc ^ ",@ " ^ h

(** [with_commas false "" list] comma-separates a list *)
let rec with_commas start acc = function
    | h :: [] -> 
            if start then h
            else append acc h
    | h :: t ->
            with_commas false (append acc h) t
    | [] -> acc

let handle_taxa_information data acc = 
    let taxa = Data.get_taxa data in
    acc ^ "@[Taxa@ being@ processed:@ @[" ^ (with_commas true "" taxa) ^ "@]@," 

let handle_character_information data acc = function
    | `Type -> 
            let printer acc str lst =
                match lst with
                | [] -> acc
                | lst ->
                        let find = fun x ->
                            match Hashtbl.find data.Data.character_specs x with
                            | Data.Static enc -> enc.Parser.SC.st_name
                            | Data.Kolmogorov d -> d.Data.dhs.Data.filename
                            | Data.Dynamic d -> d.Data.filename
                            | Data.Set -> ""
                        in
                        let names = List.map find lst in
                        acc ^ "@,@[" ^ str ^ ": @[" ^ with_commas true "" names
                        ^ "@]@]"
            in
            let acc = printer acc "Non Additive between 2 and 8 states"
            data.Data.non_additive_8 in
            let acc = printer acc "Non@ Additive@ between@ 9@ and@ 16@ states"
            data.Data.non_additive_16 in
            let acc = printer acc "Non@ Additive@ between@ 17@ and@ 32@ states"
            data.Data.non_additive_32 in
            let acc = printer acc "Non@ Additive@ with@ more@ then@ 32@ states"
            data.Data.non_additive_33 in
            let acc = printer acc "Additive" data.Data.additive in
            let acc = printer acc "Unaligned Sequences" data.Data.dynamics in
            acc

let handle_file_information data acc = function
    | `Filename ->
            let files = List.map (fun (x, _) -> x) data.Data.files in
            acc ^ "@[Loaded: @[" ^ with_commas true "" files ^ "@]@]"

let show_information trees data timer acc = function
    | `TreeInformation items ->
            (match trees with
            | Some trees ->
                let acc = acc ^ "@,@[<v 2>Trees:@ " in
                let acc = List.fold_left (handle_tree_information trees) acc items in
                acc ^ "@] "
            | None -> acc)
    | `TaxonInformation -> 
            (match data with
            | Some data -> 
                    (handle_taxa_information data acc) ^ " "
            | None -> acc)
    | `CharacterInformation lst ->
            (match data with
            | Some data ->
                let to_insert = List.fold_left (handle_character_information data) ""
                lst in
                acc ^ "@,@[<v 2>Characters: " ^ to_insert ^ "@]"
            | None -> acc)
    | `FileInformation lst ->
            (match data with
            | Some data ->
                let to_insert = List.fold_left (handle_file_information data) "" lst
                in
                acc ^ "@,@[<v 2>Files: " ^ to_insert ^ "@]"
            | None -> acc)
    | `SearchInformation _ -> acc
    | `Timer -> 
            match timer with
            | Some timer ->
                    acc ^ "@,@[<v 2>Timer:@, @[" ^ timer ^ " seconds@]@]"
            | None -> acc

let show_information trees data timer show =
    let show =
        match show with
        | Some v -> v
        | None -> 
                match data with
                | Some data -> data.Data.search_information
                | None -> []
    in
    let acc = 
        List.fold_left (show_information trees data timer) "" show
    in
    "@[<v 2>" ^ acc ^ "@]"
