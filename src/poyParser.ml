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

let () = SadmanOutput.register "PoyParser" "$Revision: 2669 $"

open StdLabels

(* The necessary types to produce the tree of the parsed input. *)
exception Illegal_Grammar
exception Name_Collision of (CharacSpec.t * string)
exception Illegal_argument
type verbosity = Min | Med | Max
type range = int * int
type distr = Poisson of float
type ulang = 
    | CInt of int * string
    | CString of string * string
    | CBase of string * string
    | Variable of string
    | Predecessor of ulang
    | Successor of ulang
    | Prepend of ulang * ulang
    | Tail of ulang
    | Head of ulang
    | Bool of string * ulang
    | IfThen of ulang * ulang * ulang
    | Application of string * ulang list
type ufuncs = string * (string * string) list * ulang
type filename = string
type spec =
    | Alphabet of string * string list * spec option
    | WordSet of string * string * range * spec option
    | IntSet of string * range * spec option
    | FProbability of distr
    | EProbability of string * (string * float) list
    | Character of string * ufuncs list * spec option
    | Protein_Sequences of (string list * (Data.dynhom_opts list option)) list
    | Synonyms of string list
    | SynonymsList of (string * string) list 
    | Static of string list
    | Ignore_files of string list
    | Ignore of string list
    | Trees of string list
    | TreesList of string list

open Camlp4.PreCast
let verbosity = ref Max
let expr = Gram.Entry.mk "expr"

let print_verbosity a b = 
    match a, !verbosity with
    | Max, Max | Med, Med | Med, Max | Min, _ -> print_endline b
    | _ -> ()

EXTEND Gram
    GLOBAL: expr;
    expr: [ [a = LIST1 [x = recexpr -> x]; `EOI -> a] ];
    recexpr: 
        [ 
            [ "define"; "character"; name = UIDENT; "as"; 
            specs = LIST1 [x = character_funcs -> x]; opt = OPT options;
            "end_def"
            -> Character (name, specs, opt) ] |
            [ "define"; "set_def"; name = UIDENT; "as"; "alphabet"; "with"; "{"; 
            alph = LIST1 [x = UIDENT -> x] SEP ";"; "}"; opt = OPT options ;
            "end_def" -> 
                Alphabet (name, alph, opt) ] |
            [ "define"; "set_def"; name = UIDENT; "as"; "wordset"; "with"; 
            alph = UIDENT; "with"; "["; min = INT; max = INT; "]";
            opt = OPT options; "end_def" -> 
                let range = int_of_string min, int_of_string max in
                WordSet (name, alph, range, opt) ] |
            [ "define"; "set_def"; name = UIDENT; "as"; "integers"; "with"; "["; 
            min = INT; max = INT; "]"; opt = OPT options; "end_def" -> 
                let range = int_of_string min, int_of_string max in
                IntSet (name, range, opt) ] 
        ];
    heading_char_funcs:
        [ 
            [ "let"; name = LIDENT; params = LIST0 [x = params -> x] -> name, params]
        ];
    character_funcs:
        [
            [ n = heading_char_funcs; "="; res = u_lang_exprs  -> 
                let (a, b) = n in a, b, res ]
        ];
    params:
        [ 
            [ name = LIDENT; ":"; set = UIDENT -> (name, set) ]
        ];
    u_function_application:
        [
            [ "prep"; b = u_lang_apps; s = u_lang_apps -> Prepend (b, s) ] 
            | [ "head"; s = u_lang_apps -> Head s ] 
            | [ "tail"; s = u_lang_apps -> Tail s ] 
            | [ "pre"; s = u_lang_apps -> Predecessor s ] 
            | [ "suc"; s = u_lang_apps -> Successor s ] 
            | [ fname = LIDENT; res = LIST0 [x = u_lang_apps -> x] -> 
                    Application (fname, res) ] (* Careful with this one! *)
        ];
    u_lang_apps:
        [
            [ a = LIDENT -> Variable a ] |
            [ a = constant -> a ] |
            [ "("; a = u_lang_exprs_1; ")" -> a ]
        ];
    u_lang_exprs:
        [ 
            [ "if"; a = u_lang_bool; "then"; b = u_lang_exprs; "else"; 
            c = u_lang_exprs -> IfThen (a, b, c) ] |
            [ r = u_lang_exprs_1 -> r ] |
            [ a = constant -> a ] 
        ];
    u_lang_exprs_1:
        [
            [ "("; r = u_lang_exprs_1; ")" -> r ] |
            [ r = u_function_application -> r ] |
            [ r = constant -> r ]
        ];
    options: 
        [ 
            [ "with"; "probability"; "{"; probs = LIST1 [x = assigned_prob -> x] SEP ";";
            "}" -> EProbability ("", probs) ] |
            [ "with"; "probability"; f = distributions -> FProbability f ] 
        ];
    distributions:
        [
            [ "poisson"; x = FLOAT -> Poisson (float_of_string x) ]
        ];
    assigned_prob: 
        [ 
            [ el = fname_for_prob; "="; prb = FLOAT -> (el, float_of_string prb)]
        ];
    fname_for_prob:
        [
            [ el = "prep" -> Gram.Token.to_string el ] |
            [ el = "head" -> Gram.Token.to_string el ] |
            [ el = "tail" -> Gram.Token.to_string el ] |
            [ el = "pre" -> Gram.Token.to_string el ] |
            [ el = "suc" -> Gram.Token.to_string el ] |
            [ el = UIDENT -> el] |
            [ el = LIDENT -> el ]
        ];
    constant:
        [
            [ r = STRING; ":"; t = UIDENT -> CString (r, t) ] | 
            [ r = INT; ":"; t = UIDENT -> CInt ((int_of_string r), t) ] |
            [ r = CHAR; ":"; t = UIDENT -> CBase (r, t) ]
        ];
    u_lang_bool:
        [
            [ v = LIDENT; "=="; r = constant -> Bool (v, r)] 
        ];
END

let of_channel ch = 
    let st = Stream.of_channel ch in
    Gram.parse expr (Loc.mk "<stream>") st

(* Process the output of the [of_channel] function above. This is a rather long
* thing, so lets see how this documentation goes.  *)
let process_tree prev lst = 
    let do_bool charspec var cnst =
        try
            let _, _, _, varspec = CharacSpec.find_variable charspec var in
            match varspec, cnst with
            | SpecIndex.Alph x, CBase (cnst, _) -> 
                    let len = AlphSpec.length ~alph:x ~elem:cnst in
                    CharacSpec.add_decoder charspec len
            | SpecIndex.Int x, CInt (cnst, _) -> 
                    let len = IntSpec.length x cnst in
                    CharacSpec.add_decoder charspec len
            | SpecIndex.Word x, CString (cnst, _) ->
                    let len = 
                        WordSpec.length ~word:x ~length:(String.length cnst) 
                    in
                    CharacSpec.add_decoder charspec len
            | _, _ -> raise Illegal_Grammar
        with
        | e ->
                prerr_string "PoyParser.do_bool\n";
                raise e
    in
    let rec do_count index charspec = function
        | CInt (i, t) ->
                begin try match SpecIndex.find index t with
                | SpecIndex.Int s -> 
                        let encoded = IntSpec.length s i in
                        CharacSpec.add_decoder (CharacSpec.count charspec t 1)
                        encoded
                | _ -> raise Illegal_Grammar
                with
                | Not_found -> assert (false);
                end;
        | CString (s, t) ->  
                begin try match SpecIndex.find index t with
                | SpecIndex.Word w ->
                        let encoded = WordSpec.wordlength w s in
                        CharacSpec.add_decoder (CharacSpec.count charspec t 1) 
                        encoded
                | SpecIndex.Alph a ->
                        let encoded = AlphSpec.length a s in
                        CharacSpec.add_decoder (CharacSpec.count charspec t 1) 
                        encoded;
                | _ -> raise Illegal_Grammar
                with
                | e -> assert (false);
                end;
        | CBase (b, t) -> 
                begin try match SpecIndex.find index t with
                | SpecIndex.Word w ->
                        let encoded = WordSpec.wordlength w b in
                        CharacSpec.add_decoder (CharacSpec.count charspec t 1) 
                        encoded;
                | SpecIndex.Alph a ->
                        let encoded = AlphSpec.length a b in
                        CharacSpec.add_decoder (CharacSpec.count charspec t 1) 
                        encoded;
                | _ -> raise Illegal_Grammar
                with 
                | e -> assert (false)
                end;
        | Variable v -> 
                let _, t, _, _ = CharacSpec.find_variable charspec v in
                CharacSpec.count charspec t 1
        | Predecessor x -> 
                let f = CharacSpec.get_string CharacSpec.Pre in
                do_count index (CharacSpec.count charspec f 1) x
        | Successor x -> 
                let f = CharacSpec.get_string CharacSpec.Suc in
                do_count index (CharacSpec.count charspec f 1) x
        | Prepend (x, y) -> 
                let f = CharacSpec.get_string CharacSpec.Prep in
                let tmp = do_count index (CharacSpec.count charspec f 1) x in
                do_count index tmp y
        | Tail x -> 
                let f = CharacSpec.get_string CharacSpec.Tl in
                let charspec = CharacSpec.count charspec f 1 in
                do_count index charspec x;
        | Head x -> 
                let f = CharacSpec.get_string CharacSpec.Hd in
                let charspec = CharacSpec.count charspec f 1 in
                do_count index charspec x
        | Bool (var, cnst) -> 
                do_bool charspec var cnst
        | IfThen (Bool (var, cnst), x, y) -> 
                let charspec = do_bool charspec var cnst in
                let charspec = do_count index charspec x in
                do_count index charspec y
        | Application (name, items) -> 
                let charspec = CharacSpec.count charspec name 1 in
                List.fold_left ~f:(do_count index) ~init:charspec items
        | IfThen _ -> raise Illegal_Grammar
    in
    let process_functions index charspec (fname, fparams, fulang) =
        (* Check if the names are being overriden and raise an exception if
        * it is the case *)
        if not (CharacSpec.exists charspec fname) then begin 
            (* Store the type of each parameter in the function *)
            let charspec = CharacSpec.add charspec fname fparams index in
            (* Count the function tree while checking the types *)
            do_count index charspec fulang
        end else raise (Name_Collision (charspec, fname))
    in
    let process_specs data item = 
        let index = data.Data.specification_index
        and chars = data.Data.character_index  in
        (* U LANGUAGE DEFINITIONS *)
        match item with
        | Alphabet (name, items, probs) -> 
                let alph = AlphSpec.create items in
                (try begin match probs with
                | None -> 
                        { data with Data.specification_index = 
                        SpecIndex.add index name (SpecIndex.Alph alph) }
                | Some (EProbability ("", enumeration)) -> 
                        let alph = 
                            AlphSpec.probs alph (AlphSpec.Enum enumeration)
                        in
                        { data with Data.specification_index = 
                        SpecIndex.add index name (SpecIndex.Alph alph) }
                | _ -> raise Illegal_Grammar
                end;
                with
                | e -> 
                        prerr_string "PoyParser.process_specs";
                        raise e)
        | WordSet (name, alphabet, (min, max), probs) -> 
                (try begin match SpecIndex.find ~index:index ~name:alphabet with
                | SpecIndex.Alph alph -> 
                        let words = WordSpec.create min max alph in
                        begin match probs with
                        | None -> 
                                let to_add = SpecIndex.Word words in
                                { data with Data.specification_index = 
                                (SpecIndex.add index name to_add) }
                        | _ -> raise Illegal_Grammar
                        end;
                | _ -> raise Illegal_Grammar
                end;
                with
                | e ->
                        raise e)
        | IntSet (name, (min, max), probs) -> 
                let ints = IntSpec.create min max in
                (try begin match probs with
                | None -> 
                        { data with Data.specification_index = 
                        SpecIndex.add index name (SpecIndex.Int ints) }
                | Some (FProbability (Poisson lambda)) -> 
                        (* We have some function that handles it.*)
                        let poisson x lambda =
                            let rec fact x res = 
                                if x = 0 then res
                                else fact (x - 1) (x * res)
                            in
                            let x' = float_of_int x in 
                            (exp (-. lambda)) *. (lambda ** x') /. (float_of_int
                            (fact x 1))
                        in
                        let deff x = poisson x lambda in
                        let ints = IntSpec.codes ints (IntSpec.Function deff) in
                        { data with Data.specification_index = 
                        SpecIndex.add index name (SpecIndex.Int ints) }
                | _ -> raise Illegal_Grammar
                end;
                with 
                | e -> 
                        prerr_string "PoyParser.process_specs";
                        raise e)
        | Character (name, flist, probs) -> 
                let res = 
                    let res = 
                        List.fold_left ~f:(process_functions index) 
                        ~init:CharacSpec.empty flist 
                    in
                    match probs with 
                    | None -> CharacSpec.prob res CharacSpec.Counter
                    | Some (EProbability ("", lst)) -> 
                            CharacSpec.prob res (CharacSpec.Assigned lst)
                    | _  -> raise Illegal_Grammar
                in
                { data with Data.specification_index = index; character_index = 
                    (name, res) :: chars }
        (* INPUT DATA PROCESSING FOR ANALYSIS *)
        | Synonyms syns ->
                (* Load a list of files containing synonyms *)
                let syns = List.map (fun x -> `Remote x) syns in
                List.fold_left ~f:Data.add_synonyms_file ~init:data syns 
        | SynonymsList lst -> 
                (* Load a list of synonyms directly *)
                List.fold_left ~f:Data.add_synonym ~init:data lst
        | Static files ->
                let files = List.map (fun x -> `Remote x) files in
                (* Load a hennig or dpread file containing static characters *)
                List.fold_left ~f:(Data.add_static_file `Hennig) ~init:data files 
        | Ignore_files files ->
                (* Load a file containing a list of taxa to be ignored in the
                * analysis *)
                let files = List.map (fun x -> `Remote x) files in
                List.fold_left ~f:Data.process_ignore_file ~init:data files 
        | Ignore names ->
                (* Load a list of taxa to be ignored *)
                List.fold_left ~f:Data.process_ignore_taxon ~init:data names 
        | Trees treefiles ->
                let process_tree_file data file =
                    let file = `Remote file in
                    let tr = Parser.Tree.of_file file in
                    { data with Data.trees = data.Data.trees @ tr }
                in
                List.fold_left ~f:process_tree_file ~init:data
                treefiles 
        | TreesList trees ->
                List.fold_left ~f:(fun acc x -> { acc with Data.trees =
                    (Parser.Tree.of_string x) @ acc.Data.trees })
                ~init:data trees 
        | Protein_Sequences _ | EProbability _ | FProbability _ ->
                failwith "Grammar not yet supported. "
    in
    let c = 
        List.fold_left ~f:process_specs ~init:prev lst
    in
    Data.categorize (Data.remove_taxa_to_ignore c)

let of_file data file =
    try
        let ch = open_in file in
        let r = of_channel ch in
        process_tree data r
    with
    | (Sys_error err) as e ->
            let msg = "System error while attempting to open the file " ^ 
            StatusCommon.escape file ^ ". The error message is \"" ^
            StatusCommon.escape err ^ "\"" in
            Status.user_message Status.Error msg;
            raise e

let set_character_weight data (c, w) = Data.transform_weight (`ReWeight (`Some
(true, [c]), w)) data

let guess_class_and_add_file annotated is_prealigned data filename =
    if Data.file_exists data filename then
        let _ =
            let filename = FileStream.filename filename in
            let msg = 
                "@[A@ file@ with@ name@ " ^ StatusCommon.escape filename ^ 
                "@ has@ previously@ " 
                ^ "been@ loaded.@ Sorry,@ I@ will@ cowardly@ refuse@ to@ "
                ^ "load@ its@ contents@ again.@ However,@ I@ will@ continue@ "
                ^ "loading@ any@ files@ remaining.@]"
            in
            Status.user_message Status.Error msg
        in
        data
    else
        let file_type_message str = 
            let msg =
                let filename = FileStream.filename filename in
                "@[Reading@ file@ " ^ StatusCommon.escape filename ^ 
                "@ of@ type@ " ^ str ^ "@]@." 
            in
            Status.user_message Status.Information msg
        in
        let add_file contents = Data.add_file data contents filename in
        let res = 
            match Parser.test_file filename with
            | Parser.Is_Poy -> 
                    let data = add_file [] in
                    file_type_message "POY";
                    (match filename with
                    | `Local filename
                    | `Remote filename -> of_file data filename)
            | Parser.Is_Clustal
            | Parser.Is_TinySeq
            | Parser.Is_Fasta | Parser.Is_Genome | Parser.Is_ASN1
            | Parser.Is_Genbank | Parser.Is_INSDSeq | Parser.Is_GBSeq
            | Parser.Is_XML | Parser.Is_NewSeq ->
                    let data = add_file [Data.Characters] in
                    file_type_message "input@ sequences";
                    Data.process_molecular_file "tcm:(1,2)" 
                    Cost_matrix.Two_D.default Cost_matrix.Three_D.default 
                    annotated Alphabet.nucleotides is_prealigned `Seq data filename
            | Parser.Is_Phylip | Parser.Is_Hennig -> 
                    let data = add_file [Data.Characters; Data.Trees] in
                    file_type_message "hennig86/Nona";
                    Data.add_static_file `Hennig data filename
            | Parser.Is_Dpread ->
                    let data = add_file [Data.Characters; Data.Trees] in
                    file_type_message "dpread file";
                    let parsed = Parser.OldHennig.of_file filename in
                    let fn = FileStream.filename filename in
                    let converted = Parser.SC.of_old_parser fn None parsed in
                    Data.add_static_parsed_file data fn converted
            | Parser.Is_Fixed_States_Dictionary ->
                    let data = add_file [] in
                    file_type_message "Fixed@ States@ Dictionary";
                    Data.process_fixed_states data (Some filename)
            | Parser.Is_Dictionary ->
                    let data = add_file [Data.Characters] in
                    file_type_message "Synonyms@ Dictionary";
                    Data.add_synonyms_file data filename
            | Parser.Is_Trees ->
                    let data = add_file [Data.Trees] in
                    file_type_message "Tree@ List";
                    Data.process_trees data filename
            | Parser.Is_Nexus -> 
                    file_type_message "Nexus@ File";
                    Data.add_static_file `Nexus data filename 
            | Parser.Is_Unknown ->
                    let data = 
                        add_file [Data.Characters; Data.Trees;
                        Data.CostMatrix] 
                    in
                    file_type_message "input@ sequences@ (default)";
                    Data.process_molecular_file 
                    "tcm:(1,2)"
                    Cost_matrix.Two_D.default Cost_matrix.Three_D.default
                    annotated Alphabet.nucleotides false `Seq data filename
            | Parser.Is_ComplexTerminals ->
                    let data = add_file [Data.Characters] in
                    file_type_message "Complex@ terminals@ definition@ file";
                    Data.process_complex_terminals data filename
        in
        Data.categorize (Data.remove_taxa_to_ignore res)

let explode_filenames files = 
    let explode_filename file = 
        let file = FileStream.filename file in
        let ch = 
            let file = 
                if Sys.os_type = "Win32" then Filename.quote file
                else if Sys.os_type = "Unix" then
                    Str.global_replace (Str.regexp "\\\\ ") "\\ " file
                else file
            in
            let line = 
                match Sys.os_type with
                | "Win32" -> ("dir /B " ^ file)
                | _ -> "ls -1 " ^ file ^ " 2> /dev/null"
            in
            Unix.open_process_in line 
            in
        let res = Parser.IgnoreList.of_channel ch in
        close_in ch;
        match res with
        | [] -> 
                let msg = "@[No@ file@ matching@ @{<b>" ^ StatusCommon.escape file ^ 
                "@}@ found.@]" in
                Status.user_message Status.Error msg;
                failwith "File not found"
        | _ -> res
    in
    List.flatten (List.map explode_filename files)

