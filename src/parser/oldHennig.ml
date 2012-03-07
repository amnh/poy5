(* POY 5.0 Alpha. A phylogenetic analysis program using Dynamic Homologies.   *)
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

(* A hennig86 file parser *)

(* The default range for an ordered type *)
let default_ordered = ( 9999, -10 )

let strict = true   (* Do not ignore errors *)
let ndebug = false  (* no debug *)

type options = 
    | Weight of (int * int) 
    | Deactivate of int 
    | Activate of int
    | Ordered of int 
    | Unordered of int
    | Sankoff of (int array array * int)
    | Trees of Tree.Parse.tree_types list
    | Unkown_option of string

exception IncompleteMatrix

(* Set of regexes to match all the supported commands in a Hennig86 file *)
let tree_re = Str.regexp "tread +\\((.*)\\)"

(* The general idea of ccost is:
   ccode ((costspec)+(charspec)+)+; *)
let multi_cc = Str.regexp "cc[a-z]* *\\([^a-z].*\\)"
let costs = Str.regexp "[Cc][Oo][Ss][a-zA-Z]* *\\[ *\\([0-9.]+\\) *\\$ *\\([0-9]+\\) *\\([0-9 ]+\\)"
let proc = Str.regexp "[Pp][Rr][Oo][Cc]/ *"

(** [build_list_between_intervals lower higher []] builds a list of
    integers from [lower] to [higher], inclusive *)
let rec build_list_between_intervals lower higher thelist =
    if higher < lower then thelist
    else build_list_between_intervals lower (higher-1) (higher::thelist)

let build_list_from_option str characters =
    if str = "." then (* All the values should be changed *)
        build_list_between_intervals 0 (characters - 1) []
    else if String.contains str '.' then begin (*range of values *) 
        let temp = Str.split (Str.regexp "\\.") str in
        let lower = int_of_string (List.nth temp 0) 
        and higher = int_of_string (List.nth temp 1) in
        build_list_between_intervals lower higher []
    end else begin
        let res = Str.split (Str.regexp " +") str in
        List.map (int_of_string) res 
    end

let rec process_an_option str characters res =
    let res = build_list_from_option str characters in
    if List.for_all (fun x -> x < characters && (-1) < x) res then res
    else 
        raise (E.Illegal_hennig86_format ("Wrong character index in " ^ str))
    
let match_num = 
    Str.regexp " *\\([1-90\\-]+\\) *\\(.*\\)"

(* used for dpread files to read in the cost matrix - read the string
*  and store as a list of integers *)   
let rec load_all_integers str l =
    match Str.string_match match_num str 0 with
    | true ->
            let v = Str.matched_group 1 str
            and r = Str.matched_group 2 str in
            let v = Pervasives.int_of_string v in
            load_all_integers r (v :: l);
    | false ->
            List.rev l 

(* used for dpread files this function takes the list of integers from
*  load_all_integers function and converts to a matrix (cost matrix) *)
let convert_list_to_matrix size lst =
    try let matrix = Array.make_matrix size size 0 in
        for i = 0 to size - 1 do 
            for j = 0 to size - 1 do
                matrix.(i).(j) <- List.nth lst (size + i*size + j);
            done;
        done;
        matrix
    with _ -> raise IncompleteMatrix

let process_ccode string n_chars =
    let r = new FileStream.string_reader string in
    let read_num def =
        try
            let str = r#read_while FileStream.is_num in
            if str = ""
            then def
            else int_of_string str
        with End_of_file -> def in
    let rec read_verb ~act ~add ~w ~acc =
        match r#getch_safe with
        | None -> List.rev acc
        | Some ch -> match ch with
          | '[' -> read_verb ~act:`Active ~add ~w ~acc
          | ']' -> read_verb ~act:`Inactive ~add ~w ~acc
          | '+' -> read_verb ~act ~add:`Additive ~w ~acc
          | '-' -> read_verb ~act ~add:`Nonadditive ~w ~acc
          | '*' -> read_verb ~act:`None ~add:`None ~w:`None ~acc
          | '/' ->
                (* see whether we can read an integer weight argument *)
                (* if not, default to 1... *)
                let weight =
                    try r#skip_ws_nl; r#read_int
                    with _ -> 1 in
                read_verb ~act ~add ~w:(`Set weight) ~acc
          | ' ' | '\t' | '\010' | '\013' ->
                read_verb ~act ~add ~w ~acc
          | c -> begin
                r#putback c;
                read_object ~act ~add ~w ~acc
            end
    and read_object ~act ~add ~w ~acc =
        (* this can be [int].[int], where the first int defaults to 0, and
           the last defaults to n_chars - 1.
           to complicate things further, we should ignore opening and
           closing parentheses... *)
        let acc = try
            r#skip_ws_nl;
            while r#match_prefix "(" do r#skip_ws_nl done;
            let nfrom = read_num 0 in
            let nto =
                if r#match_prefix "."
                then read_num (n_chars - 1)
                else nfrom in
            let chars = build_list_between_intervals nfrom nto [] in
            (* add all the properties *)
            let acc = match act with
                | `None -> acc
                | `Active -> 
                    List.fold_left (fun acc i -> (Activate i) :: acc) acc chars
                | `Inactive -> 
                    List.fold_left (fun acc i -> (Deactivate i) :: acc) acc chars
            in
            let acc = match add with
                | `None -> acc
                | `Additive -> 
                    List.fold_left (fun acc i -> (Ordered i) :: acc) acc chars
                | `Nonadditive -> 
                    List.fold_left (fun acc i -> (Unordered i) :: acc) acc chars
            in
            let acc = match w with
                | `None -> acc
                | `Set w ->
                    List.fold_left (fun acc i -> (Weight (w, i)) :: acc) acc chars
            in
            (* skip closing paren *)
            let () = 
                try ignore (r#read_while (FileStream.is_or [FileStream.is_ws_nl; FileStream.is_char ')']))
                with End_of_file -> ()
            in
            acc
        with End_of_file -> acc in
        read_verb ~act ~add ~w ~acc
    in read_verb ~act:`None ~add:`None ~w:`None ~acc:[]
            
let is_unordered_matrix matrix =
    let height = Array.length matrix in
    if height = 0 then true
    else 
        let width = Array.length matrix.(0) in
        if width = 0 then true
        else 
            try
                for x = 0 to height - 1 do
                    for y = 0 to width - 1 do
                        if x = y && matrix.(x).(y) <> 0 then raise Exit
                        else if x <> y && matrix.(x).(y) <> 1 then raise Exit
                        else ()
                    done;
                done;
                true
            with
            | Exit -> false


let process_single_command taxa_data x characters =
    try
        if not ndebug then Printf.printf "Processing: %s\n%!" x;
        if Str.string_match multi_cc x 0 then begin
            let res = Str.matched_group 1 x in
            let res = process_ccode res characters in
            res
        end
        else if Str.string_match costs x 0 then begin
            let res = Str.matched_group 1 x in
            let size = int_of_string (Str.matched_group 2 x) in
            let matrix_string = Str.matched_group 3 x in
            if not ndebug then print_endline matrix_string;
            let matrix_list = load_all_integers matrix_string [] in
            let matrix = convert_list_to_matrix size matrix_list in 
            let res = process_an_option res characters [] in 
            let processor =
                if is_unordered_matrix matrix then (fun x -> Unordered x)
                else (fun x -> Sankoff (matrix, x))
            in
            List.map processor res
        end
        else if Str.string_match tree_re x 0 then begin
            let tree = Str.matched_group 1 x in
            let trees = Tree.Parse.of_string tree in
            (* convert them to taxa *)
            try
                let m (t:Tree.Parse.tree_types):Tree.Parse.tree_types =
                    Tree.Parse.map_tree
                        (fun str ->
                            let name =
                                try fst (List.nth taxa_data (int_of_string str))
                                with _ -> str in
                            name)
                        t
                in
                let trees = List.rev_map (List.rev_map m) trees in
                List.rev_map (fun t -> Trees t) trees
            with _ -> []
        end 
        else if Str.string_match proc x 0 then begin
            []
        end else 
            [Unkown_option x]
    with 
    | IncompleteMatrix -> 
        let msg = "Incomplete Matrix in Hennig86 file. " ^ x in
        raise (E.Illegal_hennig86_format msg)
    | _ -> 
        let msg = "Illegal command in Hennig86 file. " ^ x in
        raise (E.Illegal_hennig86_format msg)


let process_options taxa_data opts y =
    let single_option_processor x =
        match process_single_command taxa_data x y with
        | [Unkown_option str] when strict -> 
            Status.user_message Status.Error ("@[Parser: Unknown Hennig86 \
                command:@ @[" ^ StatusCommon.escape str ^ "@]@]");
            assert false
        | [Unkown_option str] as res -> 
            Status.user_message Status.Error ("@[Parser: Unknown Hennig86 \
                command:@ @[" ^ StatusCommon.escape str ^ "@]@ Ignoring@]");
            res
        | res -> 
            res
    in
    let my_filter = function | [Unkown_option _] -> false | _ -> true in
    let res = List.map (single_option_processor) opts in
    let res = List.filter (my_filter) res in
    List.flatten res


(* returns int list list *)
let read_data_xread r n_chars =
    let rec line acc n_reading =
        if n_reading <= n_chars then
            match r#getch with
            | ' ' | '\n' | '\t' | '.' -> line acc n_reading
            | '[' -> multi acc [] n_reading
            | '?' | '-' -> line ([-1] :: acc) (succ n_reading)
            | c ->
                line ([(Char.code c) - (Char.code '0')] :: acc) (succ n_reading)
        else
            List.rev acc
    and multi std_acc multi_acc n_reading =
        match r#getch with
        | ' ' | '\t' | '\n' | '.' -> multi std_acc multi_acc n_reading
        | ']' -> line (List.rev multi_acc :: std_acc) (succ n_reading)
        | c ->
            multi std_acc (((Char.code c) - (Char.code '0')) :: multi_acc) n_reading
    in
    line [] 1


let read_data_dpread r n_chars =
    let clear fn n_reading acc str = match str with
        | ""  -> fn n_reading acc ""
        | str ->
            let int =
                try (int_of_string str)
                with _ -> raise (E.Illegal_hennig86_format
                                        ("Not an integer in dpread: " ^ str))
            in
            fn (succ n_reading) ([int] :: acc) ""
    and clear' fn n_reading acc str = match str with
        | ""  -> fn n_reading acc ""
        | str ->
            let int =
                try (int_of_string str)
                with _ -> raise (E.Illegal_hennig86_format
                                        ("Not an integer in dpread: " ^ str))
            in
            fn (succ n_reading) (int :: acc) ""
    in
    let rec line n_reading acc str =
        if n_reading <= n_chars then begin
            match r#getch with
                | ' ' | '\n' | '\t' | '.' ->
                    clear line n_reading acc str
                | '[' ->
                    clear
                        (fun n_reading acc str -> multi n_reading acc 0 [] str)
                        n_reading acc str
                | '-' | '?' ->
                    clear
                        (fun n_reading acc str ->
                            line (succ n_reading) ([-1] :: acc) "")
                        n_reading acc str
                | c ->
                    line n_reading acc (str ^ Char.escaped c)
        end else begin
            List.rev acc
        end
    and multi r_n r_acc m_n m_acc str =
        match r#getch with
        | ' ' | '\n' | '\t' | '.' ->
            clear' (multi r_n r_acc) m_n m_acc str
        | ']' ->
            clear'
                (fun m_n m_acc str ->
                    line (succ r_n) ((List.rev m_acc) :: r_acc) str)
                m_n m_acc str
        | c   ->
            multi r_n r_acc m_n m_acc (str ^ Char.escaped c)
    in
    line 1 [] ""


let rec read_taxa ?(acc=[]) r is_dpread taxa chars =
    r#skip_ws_nl;
    (* if we're done, check for correctness *)
    if r#match_prefix ";" then begin
        if taxa = List.length acc then begin
            List.rev acc
        end else begin
            raise (E.Illegal_hennig86_format "Number of taxa in data matrix\
                    doesn't match reported number in heading of Hennig86 file")
        end
    end else begin
        (* read the name *)
        let name = 
            r#read_while (FileStream.is_or 
                            [FileStream.is_alpha;    FileStream.is_num;
                             FileStream.is_char '_'; FileStream.is_char '-';
                             FileStream.is_char '.'; FileStream.is_char '"'; ])
        in
        if not ndebug then print_endline ("Name: " ^ name);
        r#skip_ws_nl;
        let data : int list list =
            if is_dpread
                then read_data_dpread r chars
                else read_data_xread r chars
        in
        let acc = (name, data) :: acc in
        match List.length acc with
        | x when x < taxa -> read_taxa ~acc r is_dpread taxa chars
        | x when x = taxa -> acc
        | _ ->
            raise (E.Illegal_hennig86_format "Number of taxa in data matrix\
                    doesn't match reported number in heading of Hennig86 file")
    end


let process_matrix s = read_taxa (new FileStream.string_reader s)


let rec extract_options ?(acc=[]) r =
    let maybe_cons s acc = if s = "" then acc else s :: acc in
    let rec rs str =
        match try Some r#getch with End_of_file -> None with
        | Some ';' -> str, true
        | Some '\010'
        | Some '\013' -> rs (str ^ " ")
        | Some c -> rs (str ^ Char.escaped c)
        | None -> str, false
    in
    let option, more =
        try r#skip_ws_nl; rs ""
        with End_of_file -> "", false
    in
    if not ndebug then print_endline ("Read option " ^ option);
    if more
        then extract_options ~acc:(maybe_cons option acc) r
        else List.rev (maybe_cons option acc)

(* Parses the dataset from a Hennig86 file. This is the basic parsing
* procedure of the data contents of the file *)
let parse_file characters taxa line is_dpread = 
    let taxa_data = read_taxa line is_dpread taxa characters in
    let names, data = List.split taxa_data in
    let opts = extract_options line in
    let options = process_options taxa_data opts characters in
    let trees, options =
        List.partition (function Trees _ -> true | _ -> false) options
    in
    let trees =
        List.map (function Trees a -> a | _ -> assert false) trees
    in
    names, data, options, trees

type ordtype = Is_ordered | Is_unordered | Is_sankoff

let print_ordtype x =
    match x with
    | Is_ordered -> print_string "Ordered"
    | Is_unordered -> print_string "Unordered"
    | Is_sankoff -> print_string "Sankoff"
    
type mapping_type = 
    | Do_Nothing 
    | Do_Ordered of (int, int) Hashtbl.t 
    | Do_Unordered of (int, int) Hashtbl.t
    | Do_Sankoff 
    (* Each Hashtbl is a mapping of observed numbers in the file and 
    * integer codes assigned by POY for internal representation *)

type encoding_spec = 
    { (* Encoding specifications for a given character *)
        max : int;
        min : int;
        set : All_sets.Integers.t;
        weight : int;
        active : bool;
        ordered : ordtype;
        cost_matrix : int array array;
        likelihood_model : MlModel.model option;
        observed_used : mapping_type;
        used_observed : mapping_type;
    }

let print_encoding_spec s = 
    let str_of_mapping = function
        | Do_Nothing -> "nothing"
        | Do_Ordered _ -> "ordered"
        | Do_Unordered _ -> "unordered"
        | Do_Sankoff -> "sankoff"
    in
    Printf.printf ("Max: %d\t Min: %d\n"^^
                   "Weight: %d\t Active: %b\n"^^
                   "Type: %s\tModel?: %b\n"^^
                   "ObservedUsed?: %s\n"^^
                   "UsedObserved?: %s\nSets: ")
        s.max s.min
        s.weight s.active
        (match s.ordered with | Is_ordered -> "ordered"
            | Is_unordered -> "unordered" | Is_sankoff -> "sankoff")
        (match s.likelihood_model with | Some _ -> true | None -> false)
        (str_of_mapping s.observed_used)
        (str_of_mapping s.used_observed);
        All_sets.Integers.iter (fun k -> Printf.printf "%d, " k) s.set;
        print_newline ()

let dna_encoding = 
    let codes = 
        [(1, 1); (2, 2); (4, 3); (8, 4) ; (16, 5)]
    in
    let dna_set = 
        List.fold_left (fun x (y, _) -> 
            All_sets.Integers.add y x) 
        All_sets.Integers.empty codes
    and used_observed = 
        let htb = Hashtbl.create 5 in
        List.iter (fun (a, b) -> Hashtbl.add htb a b) codes;
        htb
    in
    {
        max = 16;
        min = 1;
        set = dna_set;
        weight = 1;
        active = true;
        ordered = Is_unordered;
        cost_matrix = [||];
        likelihood_model = None;
        observed_used = Do_Nothing;
        used_observed = Do_Unordered used_observed;
    }

(* A taxon name code counter and an association list for taxa names and
* codes *)
let counter_names = ref 0 
let names = ref []


(* A character counter and an association list for characters and codes *)
let counter_characters = ref 0 
let characters = ref [] 

(* A homologous character counter and association list as in the previous *)
let counter_characters_hom = ref 0
let characters_hom = ref []

let clear_characters () =
    counter_characters_hom := 0;
    counter_characters := 0;
    characters := [];
    characters_hom := [];
    ()

let clear_taxa () =
    counter_names := 0;
    names := [];
    ()
    
let print_character_specs ch =
    output_string ch "Character\tMax\tMin\tSet\tWeight\tActive\tOrdered\n";
    let printer (b, a) =
        let print = output_string ch in
        let tab () = print "\t" in
        let print_int x = print (string_of_int x); tab () in
        let print_bool x = print (if x then "true" else "false"); tab () in
        print_int a;
        print_int b.max;
        print_int b.min;
        All_sets.Integers.iter print_int b.set;
        print_int b.weight;
        print_bool b.active;
        print_ordtype b.ordered;
        print "\n"
    in
    List.iter printer !characters

(* A generic function to associate a name with a code, used for either
* characters or taxa names *)
let some_code name lst count f =
        try
            f name !lst
        with 
        | Not_found -> 
                let n_code = !count in
                incr count;
                lst := (name, n_code) :: !lst;
                n_code

(* Symmetric to the previous. *)
let code_some code lst  =
    let rec finder = function
        | (name, mcode) :: tl when code = mcode -> name
        | _ :: tl -> finder tl
        | [] -> raise Not_found
    in
    finder !lst

let taxon_code name = 
    some_code name names counter_names List.assoc

let code_taxon code =
    code_some code names

let character_code spec =
    let a = some_code spec characters_hom counter_characters_hom List.assq
    and b = some_code spec characters counter_characters List.assoc in
    a, b

let code_character code =
    code_some code characters

let code_character_hom code =
    code_some code characters_hom

let default_encoding_specs _ =
    { 
        max = -1;
        min = Pervasives.max_int;
        set = All_sets.Integers.empty;
        weight = 1;
        active = true;
        ordered = Is_ordered;
        cost_matrix = [||];
        likelihood_model = None;
        observed_used = Do_Nothing;
        used_observed = Do_Nothing;
    }

let update_encoding_specs curr_specs data =
    let update_data x y = 
        let max = max x.max y
        and min = if (y < 0 ) then x.min else min x.min y
        and set = if (y < 0 ) then x.set else All_sets.Integers.add y x.set 
        in
        { x with max = max; min = min; set = set; }
    in
    let do_update_encoding_specs x y = List.fold_left (update_data) x y
    in
    List.map2 (do_update_encoding_specs) curr_specs data

let handle_array_position_error_and_message x = 
    let position, command = 
        match x with
        | Weight (_, a) -> a, "reweight"
        | Deactivate a -> a, "deactivate"
        | Activate a -> a, "activate"
        | Ordered a -> a, "make ordered"
        | Unordered a -> a, "make unordered"
        | Sankoff (_, a) -> a, "make sankoff"
        | _ -> failwith "Impossible state?"
    in
    Status.user_message Status.Error
    ("@[<v 2>@[Illegal@ character@ position:@]@,@[You@ asked@ me@ to@ " ^
    command ^ "@ the@ character@ " ^ string_of_int position ^ ",@ but@ " ^
    "@ that@ position@ does@ not@ exists@ in@ the@ input@ matrix.@ " ^
    "I@ am@ cancelling@ the@ processing@ of@ this@ file.@]@]");
    failwith "Illegal character number in input file"

let update_options_in_specs arr x = 
    try
        match x with
        | Weight (a, b) -> arr.(b) <- { arr.(b) with weight = a }
        | Deactivate a -> arr.(a) <- { arr.(a) with active = false }
        | Activate a -> arr.(a) <- { arr.(a) with active = true }
        | Ordered a -> arr.(a) <- { arr.(a) with ordered = Is_ordered }
        | Unordered a -> arr.(a) <- { arr.(a) with ordered = Is_unordered }
        | Sankoff (b, a) -> arr.(a) <- 
            { arr.(a) with ordered = Is_sankoff; cost_matrix = b }
        | Trees _ -> ()
        | Unkown_option _ -> 
                let msg = "Yet another impossible error" in
                raise (E.Illegal_hennig86_format msg)
    with
    | Invalid_argument _ -> handle_array_position_error_and_message x

(* Given a data set that was read from the Hennig86 matrix, the encoding
 * options of the characters and the character ranges are considered to build a set
 * of encoding parameters for each character *)
let calc_encoding_specs data_matrix options =
    match data_matrix with
    | hd :: _ -> 
          let empty = List.map (default_encoding_specs) hd in
          let param = 
              List.fold_left (update_encoding_specs) empty data_matrix 
          in
          let param = Array.of_list param in
          List.iter (update_options_in_specs param) options;
          param
    | [] -> [||]

let list_list_to_matrix x = Array.of_list (List.map (Array.of_list) x)

let objectify_channel ch = FileStream.stream_reader ch

(* Checks whether the given channel contains a possible hennig86 file and
 * extracts the header information (number of characters, taxa and the rest of the
 * text in a single long string). The regular expression reg checks the general
 * format and breaks down the different parts of the file. Note that the options at
 * the end of the file are not processed yet. *)
let is_hennig ch =
    let r = objectify_channel ch in
    r#skip_ws_nl;
    let is_dpread =
        if r#match_prefix "xread"
        then false
        else if r#match_prefix "dpread"
        then true
        else raise (E.Illegal_hennig86_format "Illegal heading") in
    r#skip_ws_nl;
    let comment =
        if r#match_prefix "'"
        then begin
        (* skip the comment *)
            let c = (r#read_excl ['\'']) in
            ignore (r#match_prefix "'");
            Some c
        end
        else None
    in
    (* read the taxon and character counts *)
    if not ndebug then begin
        print_endline ("Reading " ^ (if is_dpread then "dpread" else "xread")
                       ^ " file with comment "
                       ^ (match comment with | None -> "" | Some c -> c))
    end;
    r#skip_ws_nl; let a = r#read_int in
    r#skip_ws_nl; let b = r#read_int in
    if not ndebug then begin
        print_endline ("Reading " ^ (if is_dpread then "dpread" else "xread")
                       ^ " file with " ^ string_of_int a ^ " characters and "
                       ^ string_of_int b ^ " taxa, with comment "
                       ^ (match comment with | None -> "" | Some c -> c))
    end;
    a, b, r, is_dpread


(* Like a regular array mapping but taking two arrays a *)
let map2 a b c = 
    let b = Array.to_list b
    and c = Array.to_list c in
    let res = List.map2 a b c in
    Array.of_list res


let _architecture = 32 

let make_mapping x = 
    match x with
    | { active = true; ordered = Is_unordered; weight = w; set = s } -> 
            let size = All_sets.Integers.cardinal s in
            if size <= _architecture then begin
                let lst = All_sets.Integers.elements s in
                let lst = List.sort (compare) lst in
                let observed_used = Hashtbl.create size 
                and used_observed = Hashtbl.create size in
                let process x y =
                    Hashtbl.add observed_used y x;
                    Hashtbl.add used_observed x y;
                    x * 2
                in
                let _ = List.fold_left (process) 1 lst in
                (Do_Unordered observed_used), (Do_Unordered used_observed)
            end else begin
                let msg = "The number of states is bigger than POY's \
                           supported size." in
                raise (E.Illegal_hennig86_format msg)
            end
    | { active = true; ordered = Is_ordered; weight = w; set = s } -> 
            let size = All_sets.Integers.cardinal s in
            let lst = All_sets.Integers.elements s in
            let lst = List.sort (compare) lst in 
            let res = Hashtbl.create size in
            let process y = Hashtbl.add res y y in
            List.iter (process) lst;
            (Do_Ordered res), (Do_Ordered res)
    | { active = true; ordered = Is_sankoff; weight = w; set = s } ->
            Do_Sankoff, Do_Sankoff  
    | { active = false } -> Do_Nothing, Do_Nothing

(** [lor_list_withhash l hash] returns the logical or of the hash values of all
    of the values in [l] *)
let lor_list_withhash l hash =
    let process f =
        let proc x y = x lor (f y) in
        List.fold_left proc 0 l
    in
    match hash with
    | Some hash -> process (fun x -> Hashtbl.find hash x)
    | None -> process (fun x -> x)

let single_encoding_appl spec (data, unknown) =
    try 
        if not spec.active then FileContents.Inactive_Character
        else 
            match spec.ordered with
            | Is_unordered ->
                    let a = 
                        match spec.observed_used with
                        | Do_Unordered a -> a
                        | _ -> failwith "Unexpected"
                    in
                    FileContents.Unordered_Character (lor_list_withhash data (Some a), unknown)
            | Is_ordered ->
                let processor ((min_x, max_x) as x) y = 
                    let y_is_less_than_min = y < min_x
                    and y_is_greater_than_max = y > max_x in
                    match y_is_less_than_min, y_is_greater_than_max with
                    | false, false -> x
                    | true, false -> y, max_x
                    | false, true -> min_x, y
                    | true, true -> y, y (* This is impossible! *)
                in
                let st = default_ordered in
                let min, max = List.fold_left (processor) st data in
                FileContents.Ordered_Character (min, max, unknown)
            | Is_sankoff ->
                    FileContents.Sankoff_Character (data, unknown)
    with
    | _ -> 
            let msg = "This truly is an unexpected error." in
            raise (E.Illegal_hennig86_format msg)

let encode specs data = 
    Array.map (map2 (single_encoding_appl) specs) data

let make_list names data = 
    let res = Array.to_list data in
    List.combine res names

let merge_arrays a b =
    assert (Array.length a == Array.length b);
    let len = Array.length a in
    let first (x, _) = x
    and second (_, x) = x in
    let a = Array.init len (fun x -> first a.(x))
    and c = Array.init len (fun x -> second a.(x)) in
    Array.init len 
        (fun i -> {b.(i) with observed_used = a.(i); 
        used_observed = c.(i)})

let correct_wildcard data specs =
    let nd = Array.init (Array.length data)
        (fun i ->
             let a = data.(i) in
             Array.init (Array.length a)
                 (fun j ->
                      match a.(j) with
                      | [-1] -> 
                            let x = specs.(j) in
                            let lst = match x.ordered with
                            | Is_ordered -> x.min :: x.max :: []
                            | Is_sankoff
                            | Is_unordered ->
                                  All_sets.Integers.elements x.set in
                            lst, true
                      | a -> a, false)) in
    nd
    
let of_channel ch =
    let parser_status = Status.create "Parser" (Some 7) 
    "Parsing input files" in
    Status.report parser_status;
    let characters,taxa,text,is_dpread = is_hennig ch in
    Status.full_report ~adv:1 parser_status;
    let names,data,options,trees = parse_file characters taxa text is_dpread in
    let trees = List.map (fun x -> None,x) trees in
    Status.full_report ~adv:2 parser_status;
    let encoding_specs = calc_encoding_specs data options
    and data_matrix = list_list_to_matrix data in
    Status.full_report ~adv:3 parser_status;
    let specs = Array.map (make_mapping) encoding_specs in
    Status.full_report ~adv:4 parser_status;
    let encoding_specs = merge_arrays specs encoding_specs in
    Status.full_report ~adv:5 parser_status;
    let new_data_matrix = correct_wildcard data_matrix encoding_specs in
    Status.full_report ~adv:6 parser_status;
    let res = encode encoding_specs new_data_matrix in
    Status.full_report ~adv:7 parser_status;
    if not ndebug then begin
        print_string ("I loaded a hennig file with " ^ string_of_int
        (Array.length encoding_specs) ^ " characters.")
    end;
    Status.finished parser_status;
    encoding_specs, make_list names res,trees

let of_file f =
    let ch = FileStream.open_in f in
    let res = of_channel ch in
    close_in ch;
    res

    


let split_ordered encoding_specs taxa = 
    (* ordered and unordered characters index location in the 
    * encoding_specs array and each taxon in the list of taxa *)
    let indices_ordered, indices_unordered = 
        let ordered = ref []
        and unordered = ref []
        and len = Array.length encoding_specs in
        for i = len - 1 downto 0 do
            if encoding_specs.(i).ordered = Is_ordered then 
                ordered := i :: !ordered
            else unordered := i :: !unordered
        done;
        (Array.of_list !ordered), (Array.of_list !unordered)
    in
    let number_ordered = Array.length indices_ordered
    and number_unordered = Array.length indices_unordered in
    (* Given a single taxon, divide its characters in ordered and unordered
    * *)
    let splitter (characters, name) = 
        let ordered_chars = 
            try Array.make number_ordered characters.(0) with _ -> [||]
        in
        let unordered_chars = 
            try Array.make number_unordered characters.(0) with _ -> [||]
        in
        for i = number_ordered - 1 downto 0 do
            ordered_chars.(i) <- characters.(indices_ordered.(i));
        done;
        for i = number_unordered - 1 downto 0 do
            unordered_chars.(i) <- characters.(indices_unordered.(i));
        done;
        (ordered_chars, name), (unordered_chars, name)
    in
    (* Extract an array elements from specs stored in its index 
    * positions as stored in the indices array *)
    let extract indices specs = 
        try
            let len = Array.length indices in
            let res = Array.make len specs.(0) in
            for i = len - 1 downto 0 do
                res.(i) <- specs.(i);
            done;
            res
        with
        | _ -> [||]
    in
    let divided = List.map (splitter) taxa 
    and ordered_specs = extract indices_ordered encoding_specs
    and unordered_specs = extract indices_unordered encoding_specs in
    let all_ordered, all_unordered = List.split divided in
    (ordered_specs, all_ordered), (unordered_specs, all_unordered)

    
module Encoding = struct
    type s = encoding_spec
    let dna_encoding = dna_encoding
    (* Encoding specification *)
    let default = default_encoding_specs
    let print = print_encoding_spec
    let get_min x = x.min
    let get_max x = x.max
    let set_min x y = { x with min = y }
    let set_max x y = { x with max = y }
    let get_set x = x.set
    let set_set x s = { x with set = s }
    let set_tcm x m = { x with cost_matrix = m }
    let get_weight x = x.weight
    let set_weight x v = { x with weight = v }
    let is_active x = x.active
    let is_ordered x = x.ordered = Is_ordered
    let set_unordered x = { x with ordered = Is_unordered }
    let set_sankoff x mat = { 
        x with 
        cost_matrix = mat;
        ordered = Is_sankoff 
    }
    let is_sankoff x = x.ordered = Is_sankoff
    let get_observed_used x = 
        match x.observed_used with 
        | Do_Ordered a | Do_Unordered a -> Some a 
        | Do_Nothing -> None
        | Do_Sankoff -> None
    let get_used_observed x = 
        match x.used_observed with 
        | Do_Ordered a | Do_Unordered a -> Some a 
        | Do_Nothing -> None
        | Do_Sankoff -> None
    let set_likelihood_model model x =
        {x with likelihood_model = model;}
    let set_used_observed item x = 
        {x with used_observed = item }
    let to_string enc = 
        let active = if enc.active then "Active, " else "Inactive, "
        and ordered = if enc.ordered = Is_ordered then "Additive, " 
        else "Non Additive, "
        and weight = "Weight: " ^ string_of_int enc.weight in
        active ^ ordered ^ weight
    let has_states min max x = 
        let cardinal = All_sets.Integers.cardinal x.set in
        cardinal >= min && cardinal <= max
    let get_tcm enc = enc.cost_matrix
    let get_set enc = enc.set
    let gap_encoding gapcost =
        let res = default () in
        let res = set_min res 1 in
        let res = set_max res 2 in
        let res = set_weight res gapcost in
        let set = 
            List.fold_left (fun acc x -> All_sets.Integers.add x acc) 
            All_sets.Integers.empty [1;2] 
        in
        let res = set_set res set in
        let res = set_unordered res in
        let used_observed = 
            let htb = Hashtbl.create 2 in
            Hashtbl.add htb 1 1;
            Hashtbl.add htb 2 2;
            htb
        in
        { res with used_observed = Do_Unordered used_observed }
    let rearr_encoding recost = gap_encoding recost
end

let convert_to_Phylip_format_file hennig_filename phylip_filename = 
   let fprintf = Printf.fprintf in 
   let ch = FileStream.open_in hennig_filename in 
   let encode_arr, taxon_ls, _ = of_channel ch in 
   close_in ch;

   let phy_file = open_out phylip_filename in
   let num_taxon = List.length taxon_ls in 
   let num_char = Array.length encode_arr in 
   fprintf phy_file "%i %i\n" num_taxon num_char;


    let elt_iter elt =
        let rec elt_iter acc elt c =
            if elt = 0 then acc 
            else if (elt land 1) <> 0 then 
                elt_iter (c :: acc) (elt lsr 1) (c * 2)
            else elt_iter acc (elt lsr 1)  (c * 2)
        in
        elt_iter [] elt 1
    in

    let name_f = ref 0 in 
    List.iter (fun  (cont_arr, name) -> 
        fprintf phy_file "T%i     " (!name_f + 10000);                   
       name_f := !name_f + 1;
       Array.iteri (fun index state ->
            match state with
            | FileContents.Unordered_Character (key, _) ->       
                  let key_ls = elt_iter key in 
                  let ht_opt = Encoding.get_used_observed encode_arr.(index) in 
                  let state_ls = 
                      match ht_opt with
                      | None -> failwith "Hashtable is NONE"
                      | Some ht -> 
                            List.map (fun state_key -> 
                                          Hashtbl.find ht state_key) key_ls
                  in 

                  let dna_code = List.fold_left 
                      (fun acc state ->
                           match state with 
                           | 0 -> acc + 1 
                           | 1 -> acc + 2 
                           | 2 -> acc + 4 
                           | 3 -> acc + 8 
                           | _ -> acc + 16) 0 state_ls 
                  in
                  let new_dna_code = 
                      match dna_code > 16 with
                      | true -> dna_code - 16
                      | false ->  dna_code
                  in 

                  if new_dna_code = 16 then fprintf phy_file "-"
                  else fprintf phy_file "%s"  
                      (Alphabet.match_code new_dna_code Alphabet.nucleotides)
            | _ -> () ) cont_arr;  
           
       fprintf phy_file "\n") taxon_ls; 
    
    close_out phy_file;
    print_endline "End of converting from Hennig format to Phtlip format"


(* If taxa lists are added to the functionality (and they will), this is the
* place. Simply add a parameter with the set of taxa that will be used in
* the analysis and match it in the execution (below all these functions). I
* don't need it now, so Ron can attack that feature later. *)
let merger (all_hennig_files : (Encoding.s array * (FileContents.t array * string) list)
list) =
    (* Some internal functions necessary *)
    (* [codify_taxa files] checks the parsed hennig files contained in the
    * [files] list, see if all the files are actually shared or not, and
    * request codes for each one from the general parser taxon code
    * generator. The taxon names are then replaced by integers. *)
    let all_exist name =
        try 
            let _ = taxon_code name in
            ()
        with
        | Not_found -> raise (E.Unknown_taxon name)
    in
    let codify_taxa all_hennig_files =
        let extract_taxa (_, it) = 
            let _, b = List.split it in 
            b
        in
        let replace_names_with_codes (a, it) =
            let b, c = List.split it in
            a, List.combine b (List.map taxon_code c)
        in
        match List.map extract_taxa all_hennig_files with
        | hd :: tl ->
                (* First verify that all the files share the same taxa
                * *)
                let _ = List.map taxon_code hd in
                List.iter (List.iter all_exist) tl;
                (* Assign to each taxon a particular code from the
                * module taxon code assigned. *)
                List.map replace_names_with_codes all_hennig_files
        | [] -> []
    in
    (* Merges all the hennig files in one long list array with all the
    * characters. The function assumes that all the taxa already have codes
    * and all the files have the same taxa. *)
    let merge_all_files files = 
        (* We need to guarantee that all the files have the taxa sorted *)
        let sort_taxa (a, b) = 
            let b = Array.of_list b in
            Array.sort (fun (_, x) (_, y) -> x - y) b;
            a, Array.to_list b
        in
        (* now we need something able to merge a pair of files once they are
        * sorted. *)
        let merge_pair (fa, lsta) (fb, lstb) =
            let rec merger_of_taxa a b =
                match a, b with
                | (ch_a, nm_a) :: tla, (ch_b, nm_b) :: tlb 
                    when nm_a = nm_b ->
                        let ntl = merger_of_taxa tla tlb in
                        ((Array.append ch_a ch_b), nm_a) :: ntl
                | [], [] -> []
                | _, _ -> 
                        raise (E.Unexpected "merge_all_files with unequal \
                        length lists.");
            in
            Array.append fa fb, merger_of_taxa lsta lstb
        in
        (* So sort and merge them all now *)
        match List.map sort_taxa files with
        | hd :: tl -> List.fold_left merge_pair hd tl
        | [] -> [||], []
    in
    let assign_character_code (chars, b) =
        (Array.map character_code chars), b
    in
    (* We don't have an index of equivalent taxa, so the names
    must match in all the files. *)
    let all_files = codify_taxa all_hennig_files in
    let all_files = merge_all_files all_files in
    assign_character_code all_files

let general_character_filter f = 
    let res = List.filter (fun (a, b) -> f a) !characters in
    let _, b = List.split res in
    b

(* Return the character codes with states between min and max *)
let character_states_minmax min max = 
    let filter = fun x -> x.min >= min && x.max <= max in
    general_character_filter filter

(* Filter characters with the number of states between min and max *)
let character_states_number min max = 
    let filter = 
        fun x -> 
            let cardinal = All_sets.Integers.cardinal x.set in
            cardinal >= min && cardinal <= max
    in
    general_character_filter filter

let character_additive () = 
    let filter = fun x -> x.ordered = Is_ordered in
    general_character_filter filter

let character_nonadditive () =
    let filter = fun x ->  x.ordered = Is_unordered in
    general_character_filter filter

let character_active () =
    let filter = fun x -> x.active in
    general_character_filter filter


let character_inactive () =
    let filter = fun x -> not x.active in
    general_character_filter filter

let character_complement lst = 
    let filter = 
        fun (_, c) -> not (List.exists (fun x -> c = x) lst)
    in
    let res = List.filter filter !characters in
    let _, a = List.split res in
    a

let character_additive_minmax min max=
    let filter = fun x -> x.ordered = Is_ordered
    &&  x.min >= min && x.max <= max in
    general_character_filter filter

let character_nonadditive_minmax min max=
    let filter = fun x -> x.ordered = Is_unordered
    &&  x.min >= min && x.max <= max in
    general_character_filter filter

let character_active_minmax min max=
    let filter = fun x -> x.active &&  x.min >= min && x.max <= max in
    general_character_filter filter

let character_inactive_minmax min max=
    let filter = fun x -> not x.active &&  x.min >= min && x.max <= max in
    general_character_filter filter

let filter_matrix its m = 
    let (codes, characters) = m in
    let _, selected = 
        Array.fold_left begin
            fun (pos, accu) (hom_code, spec_code) ->
                if List.exists (fun x -> x = spec_code) its then (pos + 1, pos ::
                    accu)
                else (pos + 1, accu)
        end (0, []) codes 
    in
    let selected = List.rev selected in
    let len = List.length selected in
    (* Now we need a couple of functions, one to select the appriate columns
    * in a particular taxon row. *)
    let filter_taxon = fun (chars, tcode) ->
        if 0 < Array.length chars then begin
            let new_chars = Array.make len chars.(0) in
            let _ = List.fold_left begin
                    fun x y -> new_chars.(x) <- chars.(y); x + 1 
                end 0 selected  
            in
            new_chars, tcode
        end else ([||], tcode)
    in
    (* In the same way we now need to filter out the character codes *)
    let filter_codes = fun arr ->
        if 0 < Array.length arr then begin
            let new_codes = Array.make len arr.(0) in
            let _ = List.fold_left begin
                fun x y -> new_codes.(x) <- arr.(y); x + 1
            end 0 selected 
            in
            new_codes
        end else [||]
    in
    filter_codes codes, List.map filter_taxon characters

let flatten_matrix (codes, mtx) =
    let codes = Array.map (fun (x, _) -> x) codes in
    List.map begin
        fun (arr, taxon) ->
            Array.mapi begin
                fun pos c ->
                    c, codes.(pos)
            end arr, taxon
    end mtx

let categorize_chars (chararray, taxon_index) =
    let characters_hom = !characters_hom in
    let has_states min max x = 
        let cardinal = All_sets.Integers.cardinal x.set in
        cardinal >= min && cardinal <= max
    in
    let has_more_states min x =
        let cardinal = All_sets.Integers.cardinal x.set in
        cardinal >= min
    in
    let is_additive x =
        (x.ordered = Is_ordered) in
    let is_nonadditive x =
        (x.ordered = Is_unordered) in
    let filter_characters filter list =
        List.filter (fun (a, b) -> filter a) list in
    let get_codes list =
        let _, b = List.split list in
        b in
    
    let chars_additive = filter_characters is_additive characters_hom in
    let chars_nonadd = filter_characters is_nonadditive characters_hom in
    let chars_nonadd_8 = filter_characters (has_states 0 8) chars_nonadd in
    let chars_nonadd_16 = filter_characters (has_states 9 16) chars_nonadd
    in
    let chars_nonadd_32 = filter_characters (has_states 17 31) chars_nonadd
    in
    let chars_nonadd_33plus = filter_characters (has_more_states 32)
        chars_nonadd in

    let codes_add, codes_nonadd8, codes_nonadd16, codes_nonadd32, codes_nonadd33 =
        (get_codes chars_additive,
         get_codes chars_nonadd_8,
         get_codes chars_nonadd_16,
         get_codes chars_nonadd_32,
         get_codes chars_nonadd_33plus) in

    let array_filter array filter =
        let list = Array.to_list array in
        let list = List.filter filter list in
        let array = Array.of_list list in
        array in

    let filter_array_codes codes =
        array_filter chararray (fun (char, code) -> List.mem code codes) in

    let final_filter codes =
        (filter_array_codes codes, taxon_index) in

    (final_filter codes_add,
     final_filter codes_nonadd8,
     final_filter codes_nonadd16,
     final_filter codes_nonadd32,
     final_filter codes_nonadd33)
     
let generate_alphabet alph spec =
    let alph, gap = 
        match alph with
        | None ->
                List.map string_of_int
                (All_sets.Integers.elements spec.set),
                Alphabet.gap_repr
        | Some alph ->
                let alph = Alphabet.to_sequential alph in
                fst (List.split (Alphabet.to_list alph)), 
                Alphabet.match_code (Alphabet.get_gap alph) alph
    in
    Nexus.File.generate_alphabet alph gap

let to_new_spec ?(separator=":") filename alph spec pos =
    let newspec = 
        Nexus.File.spec_of_alph alph filename (filename ^ separator ^ string_of_int pos)
    in
    let newspec =
        if not spec.active then 
            { newspec with Nexus.File.st_eliminate = true }
        else newspec
    in
    let newspec = 
        { newspec with Nexus.File.st_weight = float_of_int
        spec.weight } 
    in
    match spec.ordered,spec.likelihood_model with
    | _, Some model ->
        { newspec with
            Nexus.File.st_type = Nexus.File.STLikelihood model }        
    | Is_unordered, _ -> newspec
    | Is_ordered, _ -> 
            { newspec with 
                Nexus.File.st_type = Nexus.File.STOrdered }
    | Is_sankoff, _ -> 
            { newspec with 
                Nexus.File.st_type = Nexus.File.STSankoff spec.cost_matrix }

let to_new_atom table_of_atoms 
(newspec : Nexus.File.static_spec) (oldspec : encoding_spec)  
data : Nexus.File.static_state =
    if Hashtbl.mem table_of_atoms data then 
        Hashtbl.find table_of_atoms data
    else
        let res = 
            match data with
            | FileContents.Ordered_Character (_, _, true)
            | FileContents.Unordered_Character (_, true)
            | FileContents.Sankoff_Character (_, true) -> None
            | FileContents.Ordered_Character (min, max, false) -> Some (`List [min; max])
            | FileContents.Unordered_Character (x, false) ->
                    let set = BitSet.create 31 in
                    let () = 
                        let rec process_bits pos v =
                            if v = 0 || pos > 32 then ()
                            else 
                                let mask = 1 lsl pos in
                                let () = 
                                    if 0 <> v land mask then
                                        BitSet.set set pos 
                                    else ()
                                in 
                                let next = (v land (lnot mask)) in
                                process_bits (pos + 1) next
                        in
                        process_bits 0 x
                    in
                    Some (`Bits set)
            | FileContents.Sankoff_Character (s, false) -> Some (`List s)
            | _ -> assert false
        in
        let () = Hashtbl.add table_of_atoms data res in
        res

let to_new_parser ?(separator=":") filename alphabet (specs, data, trees) :
    Nexus.File.nexus =
    let taxa, data = 
        let data, taxa = List.split data in
        Array.map (fun x -> Some x) (Array.of_list taxa),
        Array.of_list data
    in
    let nchars = Array.length specs 
    and ntaxa = Array.length taxa in
    let new_specs = 
        match alphabet with
        | None ->
                let table = Hashtbl.create 1667 in
                Array.iter (fun x ->
                    if not (Hashtbl.mem table x) then 
                        let a = generate_alphabet None x in
                        Hashtbl.replace table x a) specs;
                let alph = Array.map (Hashtbl.find table) specs in
                Array.init nchars (fun x ->
                    to_new_spec ~separator filename alph.(x) specs.(x) x)
        | Some alphs ->
                let table = Hashtbl.create 1667 in
                Array.iteri (fun pos x ->
                    if not (Hashtbl.mem table x) then 
                        Hashtbl.replace table x (generate_alphabet (Some x)
                        specs.(pos))) alphs;
                let alphs = Array.map (Hashtbl.find table) alphs in
                Array.init nchars (fun x ->
                    to_new_spec ~separator filename 
                    alphs.(x) 
                    specs.(x) x)
    in
    let table_of_atoms = Hashtbl.create 1667 in
    let new_data : Nexus.File.static_state array array =
        Array.init ntaxa (fun x ->
            Array.init nchars (fun y ->
                to_new_atom table_of_atoms 
                new_specs.(y) specs.(y) data.(x).(y)))
    in
    Nexus.File.fill_observed new_specs new_data;
    { (Nexus.File.empty_parsed ()) with
        Nexus.File.taxa = taxa;
        characters = new_specs;
        matrix = new_data;
        trees = trees;
    }

