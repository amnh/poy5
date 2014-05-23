(* POY 5.1.1. A phylogenetic analysis program using Dynamic Homologies.       *)
(* Copyright (C) 2014 Andrés Varón, Lin Hong, Nicholas Lucaroni, Ward Wheeler,*)
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

let (-->) b a = a b
let failwithf format = Printf.ksprintf (failwith) format
let debug = false

let debug_fo fo = 
    (fun s -> output_string stdout s; fo s)

let make_tcm_name a b = 
    "TCM" ^ string_of_int a ^ "_" ^ string_of_int b

(** LIN REMOVE THIS; DANGER ZONE!! *)
external convert_dyn :
    Data.dyna_pam_t -> Dyn_pam.dyna_pam_t = "%identity"

let make_name name =
    if name.[0] = '\'' then name else "'" ^ name ^ "'" 


let create_name_and_matrix_for_nexus tcm alph =
    let create_array a b = 
        let alph = Alphabet.to_sequential alph in
        let size = Alphabet.size alph in
        Array.init size (fun y ->
            Array.init size (fun x ->
                if x = y then 0
                else if x = (size - 1) || y = (size - 1) then b
                else a))
    in
    let rec entry tcm = match tcm with
        | Data.Substitution_Indel (a, b) -> 
            make_tcm_name a b, create_array a b
        | Data.Substitution_Indel_GapOpening (a, b, _) -> 
            make_tcm_name a b, create_array a b
        | Data.Input_file (name, mtx) ->
            make_name name, Array.of_list (List.map Array.of_list mtx)
        | Data.Input_file_GapOpening (name, mtx, _) ->
            make_name name, Array.of_list (List.map Array.of_list mtx)
        | Data.Level (d,_) -> entry d
    in
    entry tcm


let output_range format x =
    let fixit int = if format = `Hennig then int else int + 1 in
    match x with
    | `Single min -> string_of_int (fixit min)
    | `Pair (min, max) -> 
        (string_of_int (fixit min)) ^
        (if format = `Hennig then "." else "-") ^
        (string_of_int (fixit max))


let state_to_string ambig_open ambig_close sep a resolve code t =
    let f =
        let a = Alphabet.to_sequential a in
        if resolve
            then (fun x -> try Alphabet.match_code x a 
                           with | e -> Alphabet.print a; raise e)
            else string_of_int
    in
    match t with
    | None -> "?@?"
    | Some lst ->
        match Nexus.File.static_state_to_list lst with
        | []     -> "-@?"
        | [item] -> (f item) ^ "@?" 
        | lst    -> let lst = List.sort ( - ) lst in
                    ambig_open ^ String.concat sep (List.map f lst) ^ ambig_close ^ "@?"


let static_character_to_string sep ambig_open ambig_close fo resolve alph charset code =
    let () =
        try match Hashtbl.find charset code with
            | (_, `Unknown) -> fo ("?@?")
            | (spec, _) ->
                match spec with
                | Data.Stat (_,t) ->
                            fo (state_to_string ambig_open ambig_close sep alph resolve code t)
                | Data.FS code -> failwith ("compOut.static_character_to_string do we need to print fixedstates here?")
                | Data.Dyna _ -> assert false
        with
        | Not_found -> fo ("?@?")
    in
    fo sep

let create_set_data_pairs (fo : string -> unit) data code_char_pairs =
    (* associate data to each category; right now associated with likelihood *)
    let add_data charset = match charset with
        | []    -> failwith "CompOut.create_set_data_pairs; Empty Category?"
        | hd::_ ->
            let setname = Data.get_set_of_character data hd in
            begin match Hashtbl.find data.Data.character_specs hd with
                | Data.Static x ->
                    begin match x with 
                        | Data.NexusFile spec ->
                            begin match spec.Nexus.File.st_type with
                                | Nexus.File.STLikelihood m -> (setname,charset,Some m)
                                | Nexus.File.STNCM _
                                | Nexus.File.STOrdered
                                | Nexus.File.STUnordered
                                | Nexus.File.STSankoff _ -> (setname,charset,None)
                            end
                        | _ -> failwith "CompOut.create_set_data_pairs is not for fixedstates"
                    end
                | Data.Dynamic x -> (setname,charset,x.Data.lk_model)
                | Data.Kolmogorov _
                | Data.Set -> (setname,charset,None)
            end
    and print_set fo set = match set with
        | [] -> assert false
        | xs ->
            List.iter
                (fun x ->
                    let assoc = List.assoc x code_char_pairs in
                    fo (" "^(string_of_int assoc)) )
                xs
    in
    let set_pairs =
        List.map
            (fun xs ->
                try Data.categorize_likelihood_chars_by_model data (`Some xs)
                with _ -> [xs])
            (Data.categorize_characters_comp data `All)
        --> List.flatten
        --> List.map add_data
    in
    fo "@[BEGIN SETS;@\n";
    let count = ref ~-1 in
    let set_pairs =
        List.map
            (fun (n,s,m) ->
                let set = Buffer.create 1000 in
                let foset = Buffer.add_string set in
                let n = match n with
                    | None   -> incr count; "poy_"^(string_of_int !count) 
                    | Some n -> n
                in
                try foset "@[charset "; foset n; foset "= "; print_set foset s; foset ";@]@,"; 
                    fo (Buffer.contents set); (n,s,m)
                with | Not_found -> (n,s,m))
            set_pairs
    in
    if !count > 0 then fo ";@]@,END;@]@\n" else fo "END;@]@\n";
    set_pairs


let output_tree_nexus_block fo trees ic : unit =
    let leaf_only = not (List.exists (function `Cost -> true | _ -> false) ic) in
    let output cnt (name,tree) =
        let name = match name with
            | Some x -> x
            (* we gotta HOPE that the labeling it's going to be an issue *)
            | None   -> "POY_" ^ string_of_int cnt
        in
        fo ("TREE " ^ name ^ " = ");
        fo "@[";
        fo (AsciiTree.for_formatter true true leaf_only tree);
        fo ";@]@,";
        cnt + 1
    in
    match trees with
    | [] ->  ()
    | xs -> 
        fo "@[<v>@[BEGIN TREES;@]@\n";
        ignore (List.fold_left output 0 xs);
        fo "@[END;@]@\n@]"

(* The labeling generated is geared to building up the tree and it's labeling
 * and not printing the data; here we transform the structure so the code for
 * printing the results is simple *)
let transform_labeling labeling set_names =
    let find_set_from_chars ch =
        if Array.length ch > 0 then
            let n,_,_ = List.find (fun (_,s,_) -> List.mem ch.(0) s) set_names in
            n
        else 
            ""
    in
    let table = Hashtbl.create 1373 in
    Hashtbl.iter
        (fun (t,n) data -> match data with
         | [] -> ()
         | xs ->
            List.iter
                (fun (char,branch) -> match branch with
                 | None   -> ()
                 | Some x ->
                    let c = find_set_from_chars char in
                    let set = 
                        if Hashtbl.mem table (t,c)
                            then Hashtbl.find table (t,c)
                            else All_sets.StringMap.empty
                    in
                    let set = All_sets.StringMap.add n x set in
                    Hashtbl.replace table (t,c) set)
                data)
            labeling;
    table


let output_characterbranch fo in_poy labeling set_names =
    let fof format = Printf.ksprintf fo format
    and process_name n = Str.replace_first (Str.regexp "^&") "" n in
    let print_single (t,c) nbset =
        fo "@.@[@[CharacterBranch@]@\n";
        if t = "" then () else fof "@[TREES = %s;@] @\n" t;
        if c = "" then () else fof "@[CHARSET = %s;@] @\n" c;
        let count =
            All_sets.StringMap.fold
                (fun n b k ->
                    if k = 0 then fof "@[MAP @[%s %f" (process_name n) b
                             else fof ",@]@, @[%s %f" (process_name n) b;
                    k+1)
                nbset
                0
        in
        fo (if count > 0 then "@];@\n@];@\n" else "@];@\n")
    in
    if not in_poy then fo "@[@[BEGIN POY;@]@\n";
    Hashtbl.iter print_single (transform_labeling labeling set_names);
    if not in_poy then fo ";END;@]@\n"

let output_poy_nexus_block (fo : string -> unit) data labeling code_char_pairs : unit =
    let output_nexus_model (name,_,model) = match model with
        | Some m -> MlModel.output_model fo None `Nexus m (Some [name])
        | _      -> ()
    and is_empty lst = match lst with
        | [] -> true
        | _  -> false
    in
    let assoc_sets = create_set_data_pairs fo data code_char_pairs in
    let add_tcm_assumption fo code name =
        let tcm,_  = Data.get_tcm2d data code in
        let list_mat = Cost_matrix.Two_D.ori_cm_to_list tcm
        and s = Cost_matrix.Two_D.get_ori_a_sz tcm
        and a = Alphabet.to_sequential (Data.get_alphabet data code) in
        fo ("@[UserType "^name^" (StepMatrix) = "^(string_of_int s)^"@]@\n");
        fo "@[";
        for i = 0 to s-1 do
            try fo ((Alphabet.match_code i a)^" ")
            with _ -> ()
(*                begin*)
(*                    Alphabet.print a;*)
(*                    failwith ("Could not find "^string_of_int i^" in alphabet")*)
(*                end*)
        done;
        fo "@]@\n@[";
        List.fold_left
            (fun i x ->
                if i = s then (fo ("@]@\n@["^(string_of_int x)^" ");1)
                         else (fo ((string_of_int x)^" "); i+1))
            0
            list_mat --> ignore;
        fo "@];@\n";
    in
    if not (is_empty data.Data.dynamics) then begin
        fo "@[BEGIN POY;@]@\n";
        let dynamics = Array.of_list data.Data.dynamics in
        let go = Buffer.create 1000 and use_go = ref false
        and weights = Buffer.create 1000 and use_weights = ref false
        and level = Buffer.create 1000 and use_level = ref false
        and dynpam = Buffer.create 1000 and use_pam = ref false
        and tcm = Buffer.create 1000 and use_tcm = ref false
        and assump = Buffer.create 1000 and use_assump = ref false in
        Buffer.add_string go " @[GAPOPENING * POYGENERATED = ";
        Buffer.add_string tcm " @[TCM * POYGENERATED = ";
        Buffer.add_string weights " @[WTSET * POYWEIGH = ";
        Buffer.add_string level " @[LEVEL * POYLEVEL = ";
        let len = (Array.length dynamics) - 1 in
        let add_ab code pos posstr a b =
            use_tcm := true;
            use_assump := true;
            let name = make_tcm_name a b in
            add_tcm_assumption (Buffer.add_string assump) code name;
            Buffer.add_string tcm name;
            Buffer.add_string tcm ":";
            Buffer.add_string tcm posstr;
            if pos < len then Buffer.add_string tcm ","
            else Buffer.add_string tcm ";@]@,"
        and add_name code pos posstr name =
            use_tcm := true;
            use_assump := true;
            let name = make_name name in
            add_tcm_assumption (Buffer.add_string assump) code name;
            Buffer.add_string tcm name;
            Buffer.add_string tcm ":";
            Buffer.add_string tcm posstr;
            if pos < len then Buffer.add_string tcm ","
            else Buffer.add_string tcm ";@]@,"
        and add_go code pos posstr x =
            use_go := true;
            Buffer.add_string go (string_of_int x);
            Buffer.add_string go ":";
            Buffer.add_string go posstr;
            if pos < len then Buffer.add_string go ","
            else Buffer.add_string go ";@]@,";
        and add_weight code pos posstr x =
            use_weights := true;
            Buffer.add_string weights (string_of_float x);
            Buffer.add_string weights ":";
            Buffer.add_string weights posstr;
            if pos < len then Buffer.add_string weights ","
            else Buffer.add_string weights ";@]@,";
        and add_level code pos posstr x =
            use_level := true;
            Buffer.add_string level (string_of_int x);
            Buffer.add_string level ":";
            Buffer.add_string level posstr;
            if pos < len then Buffer.add_string level ","
            else Buffer.add_string level ";@]@,";
        and add_dynpam fo code data posstr : unit =
            match Data.get_dyn_state data code with
            | `Chromosome | `Genome -> 
                use_pam := true;
                use_assump := true;
                add_tcm_assumption (Buffer.add_string assump) code "";
                Dyn_pam.to_nexus fo (convert_dyn (Data.get_pam data code)) posstr
            | `SeqPrealigned | `CustomAlphabet | `Seq | `Ml | `Annotated | `Breakinv -> ()
        in
        let add_data pos code = 
            let posstr = string_of_int (pos + 1) in
            let weight = Data.get_weight code data in
            add_dynpam (Buffer.add_string dynpam) code data [posstr];
            add_weight code pos posstr weight;
            let rec add_ = function
                | Data.Substitution_Indel (a, b) -> 
                    add_ab code pos posstr a b;
                | Data.Input_file (name, _) -> 
                    add_name code pos posstr name;
                | Data.Substitution_Indel_GapOpening (a, b, x) ->
                    add_ab code pos posstr a b;
                    add_go code pos posstr x
                | Data.Input_file_GapOpening (name, _, x) ->
                    add_name code pos posstr name;
                    add_go code pos posstr x
                | Data.Level (d,x) ->
                    add_level code pos posstr x;
                    add_ d
            in
            add_ (Data.get_tcmfile data code)
        in
        Array.iteri add_data dynamics;
        fo (if !use_tcm then Buffer.contents tcm else "");
        fo (if !use_go then Buffer.contents go else "");
        fo (if !use_weights then Buffer.contents weights else "");
        fo (if !use_level then Buffer.contents level else "");
        fo (if !use_pam then Buffer.contents dynpam else "");
        List.iter output_nexus_model assoc_sets;
        output_characterbranch fo true labeling assoc_sets;
        fo "END;@]@\n";

        fo "@[BEGIN ASSUMPTIONS;@]@\n";
        fo (if !use_assump then Buffer.contents assump else "");
        fo "END;@]@\n";
        ()
    end else begin
        fo "@[BEGIN POY;@]@\n";
        List.iter output_nexus_model assoc_sets;
        output_characterbranch fo true labeling assoc_sets;
        fo "END;@\n@]"
    end


let output_character_types fo output_format resolve_a data all_of_static =
    (* We first output the non additive character types *)
    if output_format = `Nexus then fo "@[BEGIN ASSUMPTIONS;@]@\n" else ();
    let output_element name position tcm =
        let output_matrix m = 
            let buffer = Buffer.create 100 in
            Buffer.add_string buffer "@[<v 0>";
            Array.iter
                (fun x ->
                    Buffer.add_string buffer "@[<h>";
                    Array.iter
                        (fun y ->
                            let to_add =
                                if y > max_int / 8 && output_format = `Nexus then "i"
                                else (string_of_int y)
                            in
                            Buffer.add_string buffer to_add;
                            Buffer.add_string buffer " ")
                        x;
                    Buffer.add_string buffer "@]@,")
                m;
            Buffer.add_string buffer ";@,@]";
            Buffer.contents buffer
        in
        let output_codes m =
            let buffer = Buffer.create 100 in
            Buffer.add_string buffer "@[<h>";
            Array.iteri 
                (fun pos _ -> 
                    Buffer.add_string buffer (string_of_int pos);
                    Buffer.add_string buffer " ")
                m.(0);
            Buffer.add_string buffer "@]@,";
            Buffer.contents buffer;
        in
        let codes = output_codes tcm in
        let matrix = output_matrix tcm in
        if output_format = `Hennig then
            "@[<v 0>costs [ " ^ string_of_int position ^ " $" ^
            string_of_int (Array.length tcm) ^ "@," ^ codes ^ matrix
        else begin
            fo ("@[USERTYPE " ^ name ^ " STEPMATRIX =" ^ 
                string_of_int (Array.length tcm) ^ "@," ^ codes ^ matrix);
            ""
        end
    in
    let fixit int = if output_format = `Hennig then int else int + 1 in
    let table_of_matrices = Hashtbl.create 97 in
    let get_name_of_matrix suggested_name matrix =
        if Hashtbl.mem table_of_matrices matrix then
            Some (Hashtbl.find table_of_matrices matrix)
        else begin
            Hashtbl.add table_of_matrices matrix suggested_name;
            None
        end
    in
    let print_type is_last cnt
        (x : (([`Pair of (int * int) | `Single of int]) * Nexus.File.st_type) option) = 
        match x with
        | None -> ""
        | Some (range, Nexus.File.STNCM _)
        | Some (range, Nexus.File.STLikelihood _)
        | Some (range, Nexus.File.STUnordered) ->
                if output_format = `Hennig then begin
                    "@[<v 0>cc - " ^ output_range output_format range ^ ";@]@,"
                end else begin
                    "UNORD: " ^ output_range output_format range ^
                    (if not is_last then ", " else "")
                end
        | Some (range, Nexus.File.STOrdered) ->
                if output_format = `Hennig then begin
                    "@[<v 0>cc + " ^ output_range output_format range ^ ";@]@,"
                end else begin
                    "ORD: " ^ output_range output_format range ^
                    (if not is_last then ", " else "")
                end
        | Some ((`Single min) as range, Nexus.File.STSankoff matrix) ->
                let name, element = 
                    let name = "MATRIX" ^ string_of_int min in
                    match get_name_of_matrix name matrix, output_format with
                    | None, _            ->  name, output_element name min matrix
                    | Some nname,`Hennig -> nname, output_element name min matrix
                    | Some nname,`Nexus  -> nname, ""
                in
                if output_format = `Hennig then element
                else "@[" ^ name ^ ":" ^ output_range output_format range ^
                     (if not is_last then ", " else "") ^ "@]@,"
        | Some (((`Pair (min, max)) as range), Nexus.File.STSankoff matrix) ->
                if output_format = `Hennig then
                    let rec output acc i = 
                        if i > max then acc
                        else 
                            let next = 
                                output_element ("MATRIX" ^ string_of_int i) i matrix
                            in
                            output (acc ^ next) (i + 1)
                    in
                    output "" min
                else
                    let name = "MATRIX" ^ string_of_int min in
                    let res = output_element name min matrix in
                    res ^ "@[" ^ name ^ ":" ^ output_range output_format range ^
                    (if not is_last then ", " else "") ^ "@]@,"
    in
    fo "@[<v 0>";
    let acc, last, cnt =
        List.fold_left 
            (fun (acc, previous, cnt) code ->
                let spec = 
                    match Hashtbl.find data.Data.character_specs code with
                    | Data.Static x -> 
                            (match x with 
                            | Data.NexusFile enc -> enc.Nexus.File.st_type
                            | _ -> failwith "compOut.output_character_types is not for
                            fixedstates" )
                    | _ -> assert false
                in
                match previous with
                | None -> (acc, Some ((`Single cnt), spec), cnt + 1)
                | Some ((`Single min), spec') ->
                    if spec' = spec then begin
                        (acc, Some ((`Pair (min, cnt)), spec), cnt + 1)
                    end else begin
                        let acc = acc ^ (print_type false cnt previous) in
                        (acc, Some ((`Single cnt), spec), cnt + 1)
                    end
                | Some ((`Pair (min, max)), spec') ->
                    if spec' = spec then begin
                        (acc, Some ((`Pair (min, cnt)), spec), cnt + 1)
                    end else begin
                        let acc = acc ^ print_type false cnt previous in
                        (acc, Some ((`Single cnt), spec), cnt + 1)
                    end)
            ("",None, 0)
            all_of_static
    in
    let acc = acc ^ print_type true cnt last in
    if output_format = `Nexus then fo ("@[<h>TYPESET * POY = ");
    fo acc;
    if output_format = `Nexus then fo ";@]@,";
    fo "@]";
    let reweight_command, weight_separator = 
        if output_format = `Hennig then "ccode /", " "
        else "", ": "
    in
    let pos = ref ~-1 in
    let output_weights code = 
        incr pos;
        match Hashtbl.find data.Data.character_specs code with
        | Data.Static x ->
                (match x with 
                | Data.NexusFile enc ->
                    let weight = enc.Nexus.File.st_weight in 
                    if weight = 1. then ""
                    else 
                        ("@[<v 0>" ^ reweight_command ^
                        string_of_int (truncate weight) ^
                        weight_separator ^ 
                        string_of_int (fixit !pos) ^ "@]")
                | _ -> failwith "compOut.output_character_types is not for fixedstates" )
        | _ -> failwith "Sequence characters are not supported in fastwinclad"
    in
    let weights = List.map output_weights all_of_static in
    let weights = List.filter ((<>) "") weights in
    let weights = 
        String.concat (if output_format = `Nexus then ", " else ";@,")
                      weights
    in
    let weights = if output_format = `Hennig then weights ^ ";@," else weights in
    if output_format = `Nexus then fo "@[<h>WTSET * WEIGHT = " else ();
    fo weights;
    if output_format = `Nexus then fo ";@]@\n@[END;@]@,@\n" else ()


let output_character_names fo output_format resolve_a data all_of_static =
    let character_begining, state_names_beginning, character_separator, name_enclosing, to_add =
        match output_format with
        | `Hennig -> "{", " ", ";", "", 0
        | `Nexus -> " ", " / ", ",", "'", 1
    in
    let output_name (position,prev) code =
        let name = Hashtbl.find data.Data.character_codes code in
        if (output_format = `Nexus) && position > 0 then fo character_separator;
        fo ("@[<h>" ^ character_begining ^ string_of_int (position + to_add)^ " " ^
        name_enclosing ^ name ^ name_enclosing ^ state_names_beginning);
        let labels = 
            let string_labels, number_labels = 
                match Hashtbl.find data.Data.character_specs code with
                | Data.Kolmogorov _ | Data.Dynamic _ | Data.Set -> assert false
                | Data.Static x ->
                    (match x with 
                    | Data.NexusFile spec ->
                        (match spec.Nexus.File.st_labels with
                        | [] -> 
                            let f = 
                                if resolve_a then fst
                                else (fun x -> string_of_int (snd x))
                            in
                            (Data.get_alphabet data code)
                                --> Alphabet.to_list
                                --> List.sort (fun (_,a) (_,b) -> a - b)
                                --> List.map f
                                --> List.partition 
                                        (fun x -> try ignore (int_of_string x); false
                                                  with _ -> true)
                        | lst -> lst, [] )
                        | _ -> failwith "compOut.output_character_names is not for fixedstates"
                    )
            in
            string_labels @ number_labels
        in
        List.iter (fun x -> fo "@[<h>"; fo name_enclosing; fo x; 
                            fo name_enclosing; fo "@]"; fo " ")
                  labels;
        if output_format = `Hennig then fo character_separator;
        fo "@]@,";
        (position + 1,(code,position + to_add)::prev)
    in
    (* The mysterious commands for old programs *)
    match output_format with
    | `Hennig ->
            fo "@[<v 0>proc /;@,#@,";
            fo "$@,";
            fo ";@,";
            fo "@[<v 0>cn ";
            let _,codes = List.fold_left output_name (0,[]) all_of_static in
            fo ";@]@]@,";
            List.rev codes
    | `Nexus ->
            fo "@[CHARSTATELABELS@,";
            let _,codes = List.fold_left output_name (0,[]) all_of_static in
            fo ";@]@,";
            List.rev codes


let all_of_static data = 
    List.sort 
        ( - )
        (Hashtbl.fold 
            (fun c s acc -> match s with
                | Data.Static _ -> c :: acc(*do we need to filter fixed states out?*)
                | _ -> acc)
            data.Data.character_specs 
            [])


let to_nexus (data:Data.d) (trees:(string option * Tree.Parse.tree_types) list)
             (labeling:(string*string,(int array * float option) list) Hashtbl.t)
             (options:Methods.information_contained list) (filename:string option) =
    let resolve_a = true in
    let fo = Status.user_message (Status.Output (filename, false, [])) in
    let fo = if debug then debug_fo fo else fo in
    let all_of_static =  all_of_static data in

    let terminals_not_ignored =
        All_sets.IntegerMap.fold
            (fun code name acc ->
                if All_sets.Strings.mem name data.Data.ignore_taxa_set
                    then acc
                    else All_sets.IntegerMap.add code name acc)
            data.Data.taxon_codes
            All_sets.IntegerMap.empty
    in
    let terminals_sorted =
        List.sort
            (fun a b -> compare (fst a) (fst b))
            (All_sets.IntegerMap.fold (fun a b acc -> (a, b) :: acc) terminals_not_ignored [])
    in
    let output_nexus_header () = fo "#NEXUS@\n"
    and output_taxa_block () =
        fo "@[BEGIN TAXA;@]@\n";
        fo "@[DIMENSIONS NTAX=";
        fo (string_of_int (List.length terminals_sorted));
        fo ";@]@,@[TAXLABELS ";
        List.iter (fun (i, name) -> fo name; fo " ";) terminals_sorted;
        fo ";@]@,END;@,"
    in
    let output_characters_blocks () : unit =
        let output_taxa_sequences alph character_code =
            Hashtbl.iter
              (fun code characters ->
                if not (All_sets.IntegerMap.mem code terminals_not_ignored) then
                    ()
                else
                    let name = All_sets.IntegerMap.find code data.Data.taxon_codes in
                    if Hashtbl.mem characters character_code then
                        match Hashtbl.find characters character_code with
                        | Data.Stat _, _        -> assert false
                        | Data.FS _, _ -> assert false
                        | Data.Dyna _, `Unknown -> ()
                        | Data.Dyna (code, data), `Specified ->
                            assert (code = character_code);
                            fo "@["; fo name; fo " ";
                            Array.iter
                                (fun x -> fo (Sequence.to_string x.Data.seq alph))
                                data.Data.seq_arr;
                            fo "@]@\n"
                    else ())
                data.Data.taxon_characters
        in
        let output_dynamic_homology character_code = 
            match Hashtbl.find data.Data.character_specs character_code with
            | Data.Static _ | Data.Kolmogorov _ | Data.Set -> assert false
            | Data.Dynamic spec ->
                fo "@[BEGIN UNALIGNED;@]@\n";
                fo "[CHARACTER NAME: ";
                fo (Data.code_character character_code data);
                fo "]@,";
                let alphabet, symbols = 
                    if spec.Data.alph = Alphabet.nucleotides then 
                        "NUCLEOTIDE", []
                    else if spec.Data.alph = Alphabet.dna then 
                        "DNA",[]
                    else if (spec.Data.alph = Alphabet.aminoacids) 
                         || (spec.Data.alph = Alphabet.aminoacids_use_3d) then 
                        "PROTEIN", []
                    else 
                        "STANDARD", 
                        List.map fst (Alphabet.to_list spec.Data.alph)
                in
                fo " @[FORMAT DATATYPE=";
                fo alphabet;
                fo "@]@,";
                let () = match symbols with
                    | [] -> ()
                    | lst -> fo "@[SYMBOLS=\"";
                             List.iter (fun x -> fo x; fo " ") lst;
                             fo "\"@]@,"
                in
                fo ";@\n";
                fo "@[MATRIX@]@\n";
                output_taxa_sequences spec.Data.alph character_code;
                fo ";@]@[END;@]@\n"
        in
        let output_likelihood_symbols fo data codes =
            let alph,inc_gap =
                let rep = List.hd (List.tl codes) in
                try match Hashtbl.find data.Data.character_specs rep with
                    | Data.Static y ->
                        begin match y with
                            | Data.NexusFile x ->
                                let inc_gap = match x.Nexus.File.st_type with
                                    | Nexus.File.STLikelihood m ->
                                        m.MlModel.spec.MlModel.use_gap
                                    | _ -> failwith "Impossible"
                                in
                                x.Nexus.File.st_alph,inc_gap
                            | _ -> failwith "compOut.output_likelihood_symbols is not for fixedstates"
                        end
                    | _ -> failwith "Impossible"
                with | _ ->
                    Printf.printf "Codes: ";
                    Hashtbl.iter (fun c _ -> Printf.printf "%d, " c) data.Data.character_specs;
                    print_newline ();
                    failwithf "Cannot find code %d" rep
            in
            let name, symbols =
                let filter_gap lst =
                    let g = Alphabet.get_gap alph in
                    List.filter (fun (a,b) -> not (b = g)) lst
                in
                match alph with
                | x when x = Alphabet.nucleotides -> "NUCLEOTIDE", []
                | x when x = Alphabet.dna         -> "DNA",        []
                | x when x = Alphabet.aminoacids_use_3d
                                                  -> "PROTEIN",    []
                | x when x = Alphabet.aminoacids  -> "PROTEIN",    []
                | x ->
                    let f =
                        if resolve_a then (fun x-> fst x)
                                     else (fun x-> string_of_int (snd x))
                    in
                    "STANDARD", List.map f (filter_gap (Alphabet.to_list x))
            in
            fo "@[FORMAT DATATYPE="; fo name;
            fo " GAP=";
            let () =
                if resolve_a then
                    fo (Alphabet.match_code (Alphabet.get_gap alph) alph)
                else
                    fo (string_of_int (Alphabet.get_gap alph))
            in
            let () = match symbols with
                | []  -> ()
                | lst ->
                    fo " @[SYMBOLS=\"";
                    List.iter (fun x -> fo x; fo " ") lst;
                    fo "\"@]@,"
            in
            fo "@]@,"
        in
        (* Now the static homology characters, in one big matrix *)
        let number_of_static_characters = List.length all_of_static in
        let all_of_static_pairs =
            if 0 < number_of_static_characters then begin
                fo "@[BEGIN CHARACTERS;@]@\n";
                fo "@[DIMENSIONS NCHAR=";
                fo (string_of_int number_of_static_characters);
                fo ";@]@,";
                let () =
                    if Data.has_static_likelihood data then
                        output_likelihood_symbols fo data all_of_static
                    else begin
                        fo "@[FORMAT ";
                        fo "@[SYMBOLS=\"0 1 2 3 4 5 6 7 8 9 A B C D E F G H I J K L";
                        fo " M N O P Q R S T U V\"@]@,"
                    end
                in
                fo ";@]@\n";
                let all_of_static_pairs =
                    output_character_names fo `Nexus resolve_a data all_of_static
                in
                fo "@[MATRIX @\n";
                List.iter
                    (fun (code, name) ->
                        let specs = Hashtbl.find data.Data.taxon_characters code in
                        fo "@["; fo name; fo " ";
                        List.iter
                            (fun character_code ->
                                let a = Data.get_alphabet data character_code in
                                let sep = if resolve_a then "" else " " in
                                static_character_to_string
                                    sep "(" ")" fo resolve_a a specs character_code)
                            all_of_static;
                        fo "@]@\n")
                    terminals_sorted;
                fo ";@]@\n@[END;@]@\n";
                output_character_types fo `Nexus resolve_a data all_of_static;
                all_of_static_pairs
            end else begin
                []
            end
        in
        (* We print the dynamic homology characters first *)
        List.iter output_dynamic_homology data.Data.dynamics;
        let () = output_tree_nexus_block fo trees options in
        let () = output_poy_nexus_block fo data labeling all_of_static_pairs in
        ()
    in
    output_nexus_header ();
    fo "@\n";
    output_taxa_block ();
    fo "@\n";
    output_characters_blocks ();
    fo "@\n@."

let to_faswincladfile data filename =
    let has_sankoff = match data.Data.sankoff with
        | [] -> false
        | _ -> true
    in
    let fo = Status.user_message (Status.Output (filename, false, [])) in
    let all_of_static = all_of_static data in
    let number_of_characters = List.length all_of_static in
    let number_of_taxa =
        Hashtbl.fold (fun _ _ x -> x + 1) data.Data.taxon_characters 0
    in
    let sep = if has_sankoff then " " else "" in
    let output_taxon tid name =
        if All_sets.Strings.mem name data.Data.ignore_taxa_set then
            ()
        else begin
            fo name;
            fo " ";
            let () =
                let charset = Hashtbl.find data.Data.taxon_characters tid in
                List.iter 
                    (fun cc ->
                        let a = Data.get_alphabet data cc in
                        static_character_to_string sep "[" "]" fo false a charset cc)
                    all_of_static
            in
            fo "@\n@?";
        end
    in
    let output_all_taxa () = 
        All_sets.IntegerMap.iter output_taxon data.Data.taxon_codes;
        fo ";@,";
    in
    let output_header () = 
        fo (if has_sankoff then "@[<v 0>dpread " else "@[<v 0>xread ");
        fo "'Generated by POY 4.0' ";
        fo (string_of_int number_of_characters);
        fo " ";
        fo (string_of_int number_of_taxa);
        fo "@]@,@?";
    in
    fo "@[<v 0>";
    output_header ();
    output_all_taxa ();
    output_character_types fo `Hennig false data all_of_static;
    ignore (output_character_names fo `Hennig false data all_of_static);
    fo "@,";
    fo "@]@."

