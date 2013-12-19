(* POY 5.1.1. A phylogenetic analysis program using Dynamic Homologies.       *)
(* Copyright (C) 2011  Andrés Varón, Lin Hong, Nicholas Lucaroni, Ward Wheeler*)
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

open Nexus.File     (* only used for Nexus.File.nexus *)
    
let failwithf format = Printf.ksprintf failwith format

(** Creates the default_hennig character spec from the nexus.file spec *)
let default_hennig gap_handling alph equates file pos = 
    let gaps =  match gap_handling with
        | None        -> equates @ [(Alphabet.gap_repr, [])]
        | Some `Nogap -> equates @ [(Alphabet.gap_repr, [])]
        | Some `Gap   -> equates 
    in
    { st_filesource = file;
      st_name = file ^ ":" ^ string_of_int pos;
      st_alph = alph;
      st_observed = [];
      st_normal = Some 0;
      st_labels = [];
      st_weight = 1.0;
      st_type = STOrdered;
      st_equivalents = gaps;
      st_missing = "?"; 
      st_matchstate = None;
      st_gap = Alphabet.gap_repr; (* Someting that can't come from the input *)
      st_eliminate = false;
      st_case = true;
      st_used_observed = None;
      st_observed_used = None; }

(** Create a a list of length [n] of unique characters *)
let hennig_upto n = 
    let all_list =
        ["0"; "1"; "2"; "3"; "4"; "5"; "6"; "7"; "8"; "9";
         "A"; "B"; "C"; "D"; "E"; "F"; "G"; "H"; "I"; "J";
         "K"; "L"; "M"; "N"; "O"; "P"; "Q"; "R"; "S"; "T";
         "U"; "V"]
    in
    if n < 0 then failwith "Illegal number of states requested"
    else
        let rec get_up_to lst x = match x, lst with
            | 0, _      -> []
            | 1, h :: _ -> [h]
            | n, h :: t -> h :: get_up_to t (n - 1)
            | _         -> failwith "Illegal number"
        in
        get_up_to all_list n

(** Create an alphabet from a symbol table and length of alphabet *)
let hennig_for_upto is_gap_state file n pos =
        let lst = hennig_upto n in
        let gap = Alphabet.gap_repr in
        let alph, equates= 
            Nexus.File.make_symbol_alphabet gap lst []
                            [Nexus.P.Datatype Nexus.P.DStandard]
        in
        default_hennig is_gap_state alph equates file pos

(** Generate the default character spec from type of input *)
let rec generate_default of_type file pos = match of_type with
    | Some (`Continuous) ->
            let gap = Alphabet.gap_repr in
            let alph, equates = 
                Nexus.File.make_symbol_alphabet gap [] []
                                [Nexus.P.Datatype Nexus.P.Continuous]
            in
            default_hennig None alph equates file pos
    | Some (`Dna x) -> 
            let gap = Alphabet.gap_repr in
            let equ = [("0",["A"]); ("1",["C"]); ("2",["G"]);
                       ("3",["T"]); ("4",[gap])] in
            let alph, equates = 
                Nexus.File.make_symbol_alphabet gap [] equ
                                [Nexus.P.Datatype Nexus.P.Dna]
            in
            default_hennig x alph equates file pos
    | Some (`Rna x) ->
            let gap = Alphabet.gap_repr in
            let equ = [("0",["A"]); ("1",["C"]); ("2",["G"]);
                       ("3",["U"]); ("4",[gap])] in
            let alph, equates = 
                Nexus.File.make_symbol_alphabet gap [] equ
                                [Nexus.P.Datatype Nexus.P.Rna]
            in
            default_hennig x alph equates file pos
    | Some (`Protein x) ->
            let gap = Alphabet.gap_repr in
            let alph, equates =
                Nexus.File.make_symbol_alphabet gap [] []
                                [Nexus.P.Datatype Nexus.P.Protein]
            in
            default_hennig x alph equates file pos
    | Some (`Number x) when x <  9 -> hennig_for_upto None file  8 pos
    | Some (`Number x) when x < 17 -> hennig_for_upto None file 16 pos
    | Some (`Number x) when x < 33 -> hennig_for_upto None file 32 pos
    | Some (`Number _) -> generate_default (Some `Continuous) file pos
    | None             -> generate_default (Some (`Number 32)) file pos


(** Seperate the name of a taxa and the characters specification *)
let assign_names characters name = match name with
    | char :: char_name :: states_names ->
            let pos = int_of_string char in
            characters.(pos) <- 
                { characters.(pos) with st_name   = char_name;
                                        st_labels = states_names }
    | _ -> failwith "illegal character name specification"

(** Create Set of characters; range, single or polymorphic of all *)
let get_chars max chars = 
    List.flatten
        (List.map 
            (fun char ->
                let rec sequence a b acc =
                    if a > b || a >= max then acc
                    else sequence (a + 1) b (a :: acc)
                in
                match char with
                | P.All -> sequence 0 (max - 1) []
                | P.Range (a, b) -> sequence a b []
                | P.Single x -> [x])
            chars)


(** make a default 1,1 matrix; used for sankoff characters *)
let make_sankoff_matrix spec =
    let len = List.length spec.st_observed in
    Array.init len (fun x -> Array.init len (fun y -> if x = y then 0 else 1))


(** normalize columns in a continuous character set so we can fit in a
    constrained character size for vectorization *)
let normalize_continuous_characters 
        (taxa  : Nexus.File.taxon option array)
        (matrix: Nexus.File.static_state array array)
        (chars : Nexus.File.static_spec array) : bool =
    let state_max = 255 in
    (* find minimum character value for character c; ignore missing *)
    let find_by_f_character f i c =
        let ret = ref i in
        for t = 0 to ((Array.length matrix)-1) do
            begin match matrix.(t).(c) with
                | Some `List [x;y] -> ret := f !ret (f x y)
                | Some `List [-1]  -> ()
                | Some `List [x]   -> ret := f !ret x
                | None             -> ()
                | _                -> assert false
            end
        done;
        !ret
    (* normalize character states by minimum of column; also insert missing *)
    and normalize_by_minimum state_min c : unit =
        for t = 0 to ((Array.length matrix)-1) do
            begin match matrix.(t).(c) with
                | Some `List [-1]  ->
                    matrix.(t).(c) <- Some (`List [0;state_max])
                | Some `List [x]   ->
                    assert( x >= state_min );
                    assert( (x-state_min) <= state_max );
                    matrix.(t).(c) <- Some (`List [x-state_min])
                | Some `List [x;y] ->
                    assert( (x >= state_min) && ( y>= state_min) );
                    assert( ((max x y) - (min x y)) <= state_max );
                    matrix.(t).(c) <- Some (`List [x-state_min;y-state_min])
                | None -> assert false
                | _    -> assert false
            end
        done;
        chars.(c) <- { chars.(c) with st_normal = Some state_min; }
    and denormalize_by_minimum c : unit =
        chars.(c) <- { chars.(c) with st_normal = None; }
    in
    let can_normalize = ref true in
    for c = 0 to ((Array.length chars)-1) do
        let smin = find_by_f_character min max_int c in
        let smax = find_by_f_character max 0 c in
        if smax - smin > state_max then can_normalize := false
    done;
    let () =
        if !can_normalize then begin
            for c = 0 to ((Array.length chars)-1) do
                let smin = find_by_f_character min max_int c in
                normalize_by_minimum smin c
            done
        end else begin
            for c = 0 to ((Array.length chars)-1) do
                denormalize_by_minimum c
            done
        end
    in
    !can_normalize

(** process a matrix with spaces seperating the characters. Polymorphisms are
    seperated by spaces in brackets *)
let process_matrix matrix taxa characters get_row_number assign_item to_parse =
    let data = Parser.OldHennig.process_matrix to_parse
                            true (Array.length taxa) (Array.length characters)
    in
    List.iter
        (fun (name,data) ->
            let x = get_row_number name in
            Array.iteri
                (fun i states -> assign_item x i states)
                (Array.of_list data))
        (data);
    ()

(** Process a command from the lexer; Order matters, so mode can be set to
    process the character matrix properly **)
let process_command file (mode, (acc:nexus)) = function
    |P.Nstates x  -> (x, acc)
    |P.Xread data ->
        let lex = Lexing.from_string data in
        let (nch, ntaxa, to_parse) = Grammar.xread Lexer.xread lex in
        let taxa, characters, matrix =
            (* In hennig files we only allow one xread per file, that's
               common for this kind of files *)
            match acc.taxa, acc.characters, acc.matrix with
            | [||], [||], [||] ->
                Array.make ntaxa None,
                Array.init nch (generate_default mode file),
                Array.init ntaxa (fun _ -> Array.make nch None)
            | _ -> 
                failwith "We only allow one xread command per hennig file"
        (* Update the spec; only two elements, min and max are necessary *)
        and add_observed _ new_states spec =
            let new_states =
                if List.mem (Alphabet.get_gap spec.st_alph) new_states then begin
                    []
                end else begin
                    new_states @ spec.st_observed
                end
            in
            let min_observed,max_observed =
                List.fold_left min max_int new_states, List.fold_left max 0 new_states
            in
            if min_observed = max_observed then
                { spec with st_observed = [min_observed]; }
            else
                { spec with st_observed = [min_observed;max_observed]; }
        (* set the character state from a list; we don't assume two elements,
           although we only save the high and low elements in the list *)
        and add_taxastate x y new_states spec =
            if List.mem (Alphabet.get_gap spec.st_alph) new_states then begin
                Some (`List [~-1])
            end else begin
                let states = List.sort Pervasives.compare new_states in
                Some (`List states)
            end
        in
        begin match mode with
            | Some (`Continuous) ->
                process_matrix matrix taxa characters
                    (fun name -> Nexus.File.find_taxon taxa name)
                    (fun x y v ->
                        characters.(y) <- add_observed y v characters.(y);
                        matrix.(x).(y) <- add_taxastate x y v characters.(y))
                    to_parse;
                (* We determine the normalization from st_normal; can ignore
                   results returned from this function. *)
                ignore (normalize_continuous_characters taxa matrix characters);
                mode,{acc with
                        taxa = taxa;
                        characters = characters;
                        matrix = matrix; }
            | Some (`Number x) when x >= 33 ->
                process_matrix matrix taxa characters
                    (fun name -> Nexus.File.find_taxon taxa name)
                    (fun x y v ->
                        characters.(y) <- add_observed y v characters.(y);
                        matrix.(x).(y) <- add_taxastate x y v characters.(y))
                    to_parse;
                mode,{acc with
                        taxa = taxa;
                        characters = characters;
                        matrix = matrix; }
            | None
            | Some (`Dna _)
            | Some (`Number _)
            | Some (`Protein _) ->
                Nexus.File.process_matrix true `Hennig matrix taxa characters
                    (fun name -> Nexus.File.find_taxon taxa name)
                    (fun x y v -> matrix.(x).(y) <- v) to_parse;
                mode,{acc with
                        taxa = taxa;
                        characters = characters;
                        matrix = matrix; }
        end
    |P.Charname name_list ->
        let name_list = List.map (Str.split (Str.regexp "[ \t\n]+")) name_list in
        List.iter (assign_names acc.characters) name_list;
        mode, acc
    |P.Ccode char_changes ->
        List.iter
            (fun char_change ->
                let modifier, chars = match char_change with
                    |P.Additive chars ->
                        (fun x -> { x with st_type = STOrdered}), chars
                    |P.NonAdditive chars ->
                        (fun x -> { x with st_type = STUnordered}), chars
                    |P.Active chars ->
                        (fun x -> { x with st_eliminate = false }), chars
                    |P.Inactive chars ->
                        (fun x -> { x with st_eliminate = true }), chars
                    |P.Sankoff chars ->
                        (fun x ->
                            let mat = make_sankoff_matrix x in
                            { x with st_type = STSankoff mat }), chars
                    |P.Weight (v, chars) ->
                        (fun x -> { x with st_weight = float_of_int v } ), chars
                in
                List.iter
                    (fun x -> 
                        acc.characters.(x) <- modifier acc.characters.(x))
                    (get_chars (Array.length acc.characters) chars))
            char_changes;
        mode, acc
    |P.Tread new_trees ->
        let new_trees = 
            let trees = Tree.Parse.of_string new_trees in
            try let m t =
                Tree.Parse.map_tree
                    (fun str ->
                        try match acc.taxa.(int_of_string str) with
                            | Some x -> x
                            | None -> str 
                        with _ -> str)
                    t
                in
                List.rev_map (fun x -> None,List.rev_map m x) trees
            with _ -> []
        in
        mode, {acc with trees= acc.trees @ new_trees}
    |P.Cost _ -> mode, acc
    |P.EOF    -> mode, acc
    |P.Ignore -> mode, acc

(** Return true if the last command is not End of File *)
let last_command_is_not_eof = function
    | P.EOF::tl -> false
    | _         -> true

(** Load a Hennig file from a channel based on a file *)
let of_channel ch (file : string) =
    let parsed =
        let res = ref [] in
        let lex = Lexing.from_channel ch in
        while last_command_is_not_eof !res do
                let command = Grammar.command Lexer.token lex in
                res := command :: !res;
        done;
        List.rev !res
    in
    let res =
        let _, data =
            List.fold_left (process_command file) (None,(empty_parsed ())) parsed
        in
        (* Now it is time to correct the order of the terminals to 
           guarantee the default rooting of the tree. *)
        let tlen = Array.length data.taxa
        and mlen = Array.length data.matrix in
        assert (tlen >= mlen);
        let taxa = Array.init tlen (fun x -> data.taxa.(tlen - x - 1)) 
        and matrix = Array.init mlen (fun x -> data.matrix.(mlen - x - 1)) in
        {data with taxa=taxa; matrix=matrix}
    in
    res
