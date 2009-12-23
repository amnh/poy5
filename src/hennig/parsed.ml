open Nexus.Parsed 

let default_hennig gap_handling alph equates file pos = 
    let gaps =  match gap_handling with
        | None -> equates @ [(Alphabet.gap_repr, [])]
        | Some `Nogap -> equates @ [(Alphabet.gap_repr, [])]
        | Some `Gap -> equates 
    in

    { st_filesource = file;
      st_name = file ^ ":" ^ string_of_int pos;
      st_alph = alph;
      st_observed = [];
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
      st_observed_used = None}

let hennig_upto n = 
    let all_list = ["0"; "1"; "2"; "3"; "4"; "5"; "6"; "7"; "8"; "9";
    "A"; "B"; "C"; "D"; "E"; "F"; "G"; "H"; "I"; "J"; "K"; "L"; "M";
    "N"; "O"; "P"; "Q"; "R"; "S"; "T"; "U"; "V"] in
    if n < 0 then failwith "Illegal number of states requested"
    else
        let rec get_up_to lst x =
            match x, lst with
            | 0, _ -> []
            | 1, h :: _ -> [h]
            | n, h :: t -> h :: get_up_to t (n - 1)
            | _ -> failwith "Illegal number"
        in
        get_up_to all_list n

let hennig_for_upto is_gap_state file n pos =
        let lst = hennig_upto n in
        let alph, equates= 
            Nexus.Parsed.make_symbol_alphabet "7" lst [] [Nexus.File.Datatype
            Nexus.File.DStandard]
        in
        default_hennig is_gap_state alph equates file pos

let rec generate_default of_type file pos =
    
    match of_type with
    | Some (`Dna x) -> 
            let gap = "-" in
            let equ =  [("0",["A"]); ("1",["C"]);
                        ("2",["G"]); ("3",["T"]); ("4",[gap])] in
            let alph, equates = 
                Nexus.Parsed.make_symbol_alphabet gap [] equ [Nexus.File.Datatype
                Nexus.File.Dna]
            in
            default_hennig x alph equates file pos

    | Some (`Rna x) ->
            let gap = "-" in
            let equ =  [("0",["A"]); ("1",["C"]);
                        ("2",["G"]); ("3",["U"]); ("4",[gap])]  in
            let alph, equates = 
                Nexus.Parsed.make_symbol_alphabet "-" [] equ [Nexus.File.Datatype
                Nexus.File.Rna]
            in
            default_hennig x alph equates file pos
    | Some (`Protein x) ->
            let alph, equates = 
                Nexus.Parsed.make_symbol_alphabet "-" [] [] [Nexus.File.Datatype
                Nexus.File.Protein]
            in
            default_hennig x alph equates file pos
    | Some (`Number x) ->
            if x < 9 then
                hennig_for_upto None file 8 pos
            else if x < 17 then
                hennig_for_upto None file 16 pos
            else if x < 33 then
                hennig_for_upto None file 32 pos
            else failwith "We can't hold it, can you?"
    | None -> generate_default (Some (`Number 32)) file pos

let assign_names characters name =
    match name with
    | char :: char_name :: states_names ->
            let pos = int_of_string char in
            characters.(pos) <- 
                { characters.(pos) with st_name = char_name; st_labels =
                    states_names }
    | _ -> failwith "illegal character name specification"

let get_chars max chars = 
    List.flatten (List.map (fun char ->
        let rec sequence a b acc =
            if a > b || a >= max then acc
            else sequence (a + 1) b (a :: acc)
        in
        match char with
        | File.All -> sequence 0 (max - 1) []
        | File.Range (a, b) -> sequence a b []
        | File.Single x -> [x]) chars)

let make_sankoff_matrix spec =
    let len = List.length spec.st_observed in
    Array.init len (fun x ->
        Array.init len (fun y ->
            if x = y then 0
            else 1))

let process_command file (mode, (acc:nexus)) = function
    | File.Nstates x -> (x, acc)
    | File.Xread data ->
            let lex = Lexing.from_string data in
            let (nch, ntaxa, to_parse) = 
                Grammar.xread Lexer.xread lex
            in
            let taxa, characters, matrix =
                (* In hennig files we only allow one xread per file, that's
                * common for this kind of files *)
                match acc.taxa, acc.characters, acc.matrix with
                | [||], [||], [||] ->
                        (* We are OK to continue *)
                        Array.make ntaxa None, Array.init nch
                        (generate_default mode file), 
                        Array.init ntaxa (fun _ -> Array.make nch None)
                | _ -> 
                        failwith 
                        "We only allow one xread command per hennig file"
            in
            (* Now we can parse the contents using the default parser *)
            Nexus.Parsed.process_matrix true `Hennig matrix taxa characters 
            (fun name -> Nexus.Parsed.find_taxon taxa name)
            (fun x y v -> matrix.(x).(y) <- v) to_parse;
            mode, {acc with
                        taxa = taxa;
                        characters = characters;
                        matrix = matrix}
    | File.Charname name_list ->
            (* We need to parse each of the charname entries *)
            let name_list = 
                List.map (Str.split (Str.regexp "[ \t\n]+"))
                name_list
            in
            List.iter (assign_names acc.characters) name_list;
            mode, acc
    | File.Ignore -> mode, acc
    | File.Ccode char_changes ->
            List.iter (fun char_change ->
                let modifier, chars = 
                    match char_change with
                    | File.Additive chars -> 
                            (fun x -> { x with st_type = STOrdered}), chars
                    | File.NonAdditive chars -> 
                            (fun x -> 
                                { x with st_type = STUnordered}), chars
                    | File.Active chars ->
                            (fun x -> { x with st_eliminate = false }), chars
                    | File.Inactive chars ->
                            (fun x -> { x with st_eliminate = true }), chars
                    | File.Sankoff chars ->
                            (fun x -> { x with st_type = STSankoff
                            (make_sankoff_matrix x) }), chars
                    | File.Weight (v, chars) ->
                            (fun x -> { x with st_weight = float_of_int v } ), 
                            chars
                in
                List.iter (fun x -> 
                    acc.characters.(x) <- modifier acc.characters.(x))
                (get_chars (Array.length acc.characters) chars)) char_changes;
            mode, acc
    | File.Tread new_trees ->
            let new_trees = 
                let trees = Tree.Parse.of_string new_trees in
                (* convert them to taxa *)
                try
                    let m t = Tree.Parse.map_tree
                        (fun str ->
                            try 
                                match acc.taxa.(int_of_string str) with
                                | Some x -> x
                                | None -> str 
                            with _ -> str)
                        t
                    in
                    List.rev_map (fun x -> None,List.rev_map m x) trees
                with _ -> []
            in
            mode, {acc with trees= acc.trees @ new_trees}
    | _ -> mode, acc

let of_channel ch (file : string) =
    let parsed = 
        let res = ref [] in
        let lex = Lexing.from_channel ch in
        try 
            while true do
                let command = Grammar.command Lexer.token lex in
                res := command :: !res;
            done;
            []
        with
        | Lexer.Eof -> List.rev !res
    in
    let res =
        let _, data = List.fold_left (process_command file)
                                     (None,(empty_parsed ()))
                                     parsed
        in
        (* Now it is time to correct the order of the terminals to 
        * guarantee the default rooting of the tree. *)
        let tlen = Array.length data.taxa
        and mlen = Array.length data.matrix in
        assert (tlen >= mlen);
        let taxa = Array.init tlen (fun x -> data.taxa.(tlen - x - 1)) 
        and matrix = Array.init mlen (fun x -> data.matrix.(mlen - x - 1)) in
        {data with taxa=taxa; matrix=matrix}
    in
    res
