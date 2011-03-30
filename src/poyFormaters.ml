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

let () = SadmanOutput.register "PoyFormaters" "$Revision: 2400 $"

exception Illegal_formater of string

let user_messagef st format = Printf.ksprintf (Status.user_message st) format

let is_xml_filename filename =
    match filename with 
    | None -> false
    | Some filename -> Filename.check_suffix filename ".xml" 

let sort_matrix mtx = 
    let splitter str = 
        match Str.split (Str.regexp ":") str with
        | [_] | [] -> str, ""
        | lst ->
                let rec last = function
                    | [h] -> [], h
                    | h :: t -> 
                            let a, b = last t in
                            h :: a, b
                    | [] -> 
                            failwith "this case was already filtered"
                in
                let a, b = last lst in 
                String.concat ":" a, b
    in
    let comparator a b =
        let _, r = 
            Array.fold_left (fun (p, c) v -> 
                if c <> 0 then (p + 1, c)
                else 
                    let vpre, vpost = splitter v 
                    and bpre, bpost = splitter b.(p) in
                    let comparison = 
                        let init = String.compare vpre bpre in
                        if  0 <> init then init
                        else 
                            try 
                                let v = int_of_string vpost
                                and b = int_of_string bpost in
                                v - b
                            with
                            | _ -> 0
                    in
                    (p + 1, comparison)) (0, 0) a
        in
        r
    in
    Array.iteri (fun x str -> mtx.(0).(x) <- " " ^ str) mtx.(0);
    Array.sort comparator mtx

let build_contents_row ((_, attributes, contents) : Xml.xml)  =
    let mapper (_, value) = StatusCommon.escape (Xml.value_to_string value) in
    match contents with
    | #Xml.structured_xml -> (* We will ignore the structured contents for a table *)
            let res = Array.of_list attributes in
            Array.map mapper res
    | #Xml.unstructured as v ->
            let res = 
                Array.of_list (List.rev ( ("", v) :: (List.rev attributes)))
            in
            Array.map mapper res

let build_names_row ?(contents_name = " ") (_, attributes, _) =
    let mapper (name, _) = "@{<b>" ^ name ^ "@}" in
    let res = 
        Array.of_list (List.rev ((contents_name, `String contents_name) :: 
            (List.rev attributes))) 
    in
    Array.map mapper res

type t = [ Xml.unstructured | Xml.xml Xml.structured_xml ]

let build_set (contents : t) = 
    let res = Buffer.create 16 in
    let rec add_items (_, _, x) =
        match x with
        | #Xml.unstructured as v -> 
                Buffer.add_string res (Xml.value_to_string v);
                Buffer.add_string res ", ";
        | #Xml.structured as x ->
                let x = Xml.eagerly_compute x in
                Sexpr.leaf_iter add_items x
        | `CDATA _ -> failwith "Can't handle CDATA contents in poyFormaters"
    in
    let _ =
        match contents with
        | #Xml.structured as todo ->
                let todo = Xml.eagerly_compute todo in
                Buffer.add_string res "{";
                Sexpr.leaf_iter add_items todo;
                Buffer.add_string res "}";
        | #Xml.unstructured as v -> 
                Buffer.add_string res (Xml.value_to_string v)
        | `CDATA _ -> failwith "Can't handle CDATA contents in poyFormaters"
    in
    (Buffer.contents res)

let build_table_with_contents_as_set ?(cn = " ") (lst :
    Xml.xml list) = 
    let mapper ((a, b, contents) : Xml.xml) : Xml.xml = 
        (a, b, `String (build_set contents))
    in
    match lst with
    | h :: _ ->
            let res = List.map (function x -> build_contents_row (mapper x)) lst 
            and names = build_names_row ~contents_name:cn h in
            Array.of_list (names :: res)
    | [] -> [|[||]|]

let character_type_as_attribute (tag, attributes, contents) =
    let t = "Type" in
    let attr = 
        if tag = Xml.Characters.additive then 
            "Additive"
        else if tag = Xml.Characters.nonadditive then 
            "Non Additive"
        else if tag = Xml.Characters.sankoff then 
            "Sankoff"
        else if tag = Xml.Characters.molecular then 
            "Molecular"
        else if tag = Xml.Characters.kolmogorov then 
            "Kolmogorov"
        else failwith "Unexpected"
    in
    (tag, ((t, attr) :: attributes), contents)

let filter_tag tag (item : Xml.xml) : Xml.xml list =
    let rec build acc ((mytag, _, contents) as item) =
        let nacc = 
            if tag  = mytag then 
                item :: acc
            else acc
        in
        match contents with
        | #Xml.structured as x ->
                let x = Xml.eagerly_compute x in
                Sexpr.fold_left build nacc x
        | #Xml.unstructured -> nacc
        | `CDATA _ -> failwith "Can't handle CDATA contents in poyFormaters"
    in
    List.rev (build [] item)

let rec build_values_as_list st (_, _, contents) =
    match contents with
    | #Xml.structured as x -> 
            let x = Xml.eagerly_compute x in
            Sexpr.leaf_iter (build_values_as_list st) x
    | #Xml.unstructured as v -> 
            Status.user_message st (StatusCommon.escape (Xml.value_to_string v));
            Status.user_message st ";@ "
    | `CDATA _ -> failwith "Can't handle CDATA contents in poyFormaters"
    
let output_rows st matrix = 
    Status.user_message st (string_of_int ((Array.length matrix) - 1))

let output_characters st (characters : Xml.xml)=
    let output_table_of_list title cn (lst : Xml.xml list) =
        let matrix = build_table_with_contents_as_set ~cn:cn lst in
        sort_matrix matrix;
        Status.user_message st ("@[<v 4>@{<b>" ^ title ^ 
        "@}@\n@[<v 0>@\n@{<u>Total@} ");
        output_rows st matrix;
        Status.user_message st "@\n";
        Status.output_table st matrix;
        Status.user_message st "@]@]@\n"
    in
    (* For each tag we will output the appropriate character *)
    let additive = filter_tag Xml.Characters.additive characters
    and nonadditive = filter_tag Xml.Characters.nonadditive characters
    and sankoff = filter_tag Xml.Characters.sankoff characters
    and molecular = filter_tag Xml.Characters.molecular characters
    and kolmogorov = filter_tag Xml.Characters.kolmogorov characters in
    Status.user_message st "@[<v 4>@{<b>Characters@}@\n@[<v 0>";
    output_table_of_list "Non Additive" "Set" nonadditive;
    output_table_of_list "Additive" "Range" additive;
    output_table_of_list "Sankoff" "Set" sankoff;
    output_table_of_list "Molecular" " " molecular;
    output_table_of_list "Kolmogorov" " " kolmogorov;
    Status.user_message st "@]@]@\n"

let output_taxa st (_, _, taxa) =
    match taxa with
    | #Xml.structured as x ->
            let x = Xml.eagerly_compute x in
            let lst = Sexpr.to_list x in
            let mtx = build_table_with_contents_as_set lst in
            sort_matrix mtx;
            Status.user_message st 
            "@[<v 4>@{<b>Taxa@}@\n@[<v 0>@\n@{<u>Total@} ";
            output_rows st mtx;
            Status.user_message st "@\n";
            Status.output_table st mtx;
            Status.user_message st "@]@]@\n"
    | #Xml.unstructured -> ()
    | `CDATA _ -> failwith "Can't handle CDATA contents in poyFormaters"

(* [output_files st c] outputs the list of files from the output c as a list of
* elements separated by semicolon. All the contents, all the values are
* converted into a list. *)
let output_list title st c =
    Status.user_message st ("@[<v 4>@{<b>" ^ title ^ "@}@\n@[");
    build_values_as_list st c;
    Status.user_message st "@]@]@\n"

(* This is a general purpose formatter for attributes.*)
let format_attributes st attributes = 
    let format_attribute (a, b) =
        Status.user_message st "@[";
        Status.user_message st (a);
        Status.user_message st " : ";
        Status.user_message st (Xml.value_to_string b);
        Status.user_message st "@]@\n";
    in
    Status.user_message st "@[<v 0>@\n";
    List.iter format_attribute attributes;
    Status.user_message st "@]"

let rec aux_data_to_status st ((tag, attributes, contents) as c : Xml.xml) =
    if tag = Xml.Data.characters then output_characters st c
    else if tag = Xml.Data.taxa then output_taxa st c
    (*
    else if tag = Xml.Data.files then output_list "Files" st c
    *)
    else if tag = Xml.Data.ignored_taxa then output_list "Ignored taxa" st c
    else if tag = Xml.Data.ignored_characters then output_list "Ignored characters" st c
    (* TODO else if tag = Xml.Data.synonyms then output_synonyms st c *)
    else begin
        let str = 
            match contents with
            | #Xml.structured ->  "@[<v 4>@ "
            | #Xml.unstructured -> "@[@ "
            | `CDATA _ -> failwith "Can't handle CDATA contents in poyFormaters"
        in
        Status.user_message st "@{<b>";
        Status.user_message st tag;
        Status.user_message st "@}";
        Status.user_message st str;
        Status.user_message st "@\n";
        format_attributes st attributes;
        Status.user_message st "@ ";
        begin match contents with
        | #Xml.structured as sexpr ->
                let sexpr = Xml.eagerly_compute sexpr in
                Sexpr.leaf_iter (aux_data_to_status st) sexpr
        | #Xml.unstructured as v ->
                let value = StatusCommon.escape (Xml.value_to_string v) in
                Status.user_message st value;
        | `CDATA _ -> failwith "Can't handle CDATA contents in poyFormaters"
        end;
        Status.user_message st "@]@\n@]@\n"
    end

(* [data_to_status filename tag] outputs the contents of Data.to_formatter in a table
* format presented in the ncurses and non-ncurses interface. *)
let data_to_status filename tag =
    StatusCommon.Files.set_margin filename 0;
    let st = Status.Output (filename, false, []) in
    if is_xml_filename filename then
        let filename = 
            match filename with | Some x -> x | None -> assert false in
        let ch = StatusCommon.Files.channel filename in
        match ch with
        | `NotCompressed ch -> Xml.to_file ch tag
        | `Zlib ch -> 
                failwith "Output fo XML in compressed format is not supported"
    else begin
        aux_data_to_status st tag;
        Status.user_message st "%!"
    end

let get_name_class_and_cost attr = 
    let get_them ((a, cclass, b) as acc) (tag, value) =
        if tag = Xml.Characters.name then
            (Some value), cclass, b
        else if tag = Xml.Characters.cclass then 
            a, (Some value), b
        else if tag = Xml.Characters.cost then
            a, cclass, (Some value)
        else acc
    in
    match List.fold_left get_them (None, None, None) attr with
    | Some name, Some cclass, Some cost -> name, cclass, cost 
    | _, _, _ -> raise (Illegal_formater "get_name_class_and_cost")

let find_item tag lst = (List.assoc tag lst)

let get_recost attr = find_item Xml.Characters.recost attr

let get_ref_code attr = find_item Xml.Characters.ref_code attr

let get_map attr = `String "Available only in XML format"


let min_and_max ((x, y) as acc) (a, _, c) = 
    match c with
    | #Xml.unstructured as c ->
            let c = StatusCommon.escape (Xml.value_to_string c) in
            if a = Xml.Characters.min then
                (Some (c), y)
            else if a = Xml.Characters.max then
                (x, Some (c))
            else acc
    | #Xml.structured_xml ->
            raise (Illegal_formater "min_and_max")

let addcs_to_formater ((tag, attr, cont) : Xml.xml) = 
    if tag = Xml.Characters.additive then begin
        let name, cclass, cost = get_name_class_and_cost attr in
        let recost = `Float 0. in 
        let chrom_ref = `String "-" in 
        let map = chrom_ref in 
        let minmax = 
            match cont with
            | #Xml.structured as x -> 
                    let x = Xml.eagerly_compute x in
                    Sexpr.fold_left min_and_max (None, None) x
            | _ -> raise (Illegal_formater "addcs_to_formater")
        in
        match minmax with
        | Some min, Some max -> 
                let c = Buffer.create 10 in
                Buffer.add_string c "[";
                Buffer.add_string c min;
                Buffer.add_string c ", ";
                Buffer.add_string c max;
                Buffer.add_string c "]";
                [|name; cclass; cost; recost;
                chrom_ref; map; 
                `String (Buffer.contents c)|]
        | _ -> raise (Illegal_formater "addcs_to_formater 2")
    end else raise (Illegal_formater "addcs_to_formater 3")

let nonaddcs_to_formater ((tag, attr, cont) : Xml.xml) =
    if tag = Xml.Characters.nonadditive then begin
        let (_, name) = 
            List.find (fun (a, _) -> a = Xml.Characters.name) attr
        in
        let (_, cclass) = 
            List.find (fun (a, _) -> a = Xml.Characters.cclass) attr
        in
        let (_, cost) =
            try List.find (fun (a, _) -> a = Xml.Characters.cost) attr 
            with Not_found -> 
                let str = " " in
                str, `String str
        in
        let recost = `String "0." in 
        let chrom_ref = `String "0." in 
        let map = `String "-" in 
        let states = 
            match cont with
            | #Xml.structured as x -> (* This is what we expect *)
                    let x = Xml.eagerly_compute x in
                    let c = Buffer.create 10 in
                    Buffer.add_string c "{";
                    Sexpr.leaf_iter (fun (_, _, x) ->
                        match x with
                        | #Xml.unstructured as x -> 
                                Buffer.add_string c (Xml.value_to_string x);
                                Buffer.add_string c ",";
                        | _ ->
                                raise (Illegal_formater "nonaddcs_to_formater"))
                        x;
                    Buffer.add_string c "}";
                    StatusCommon.escape (Buffer.contents c)
            | _ -> raise (Illegal_formater "nonaddcs_to_formater 2")
        in
        [| name; cclass; cost; recost; chrom_ref; map; `String states |]
    end else raise (Illegal_formater "nonaddcs_to_formater")

let sankcs_to_formater ((tag, attr, cont) : Xml.xml) =
    if tag = Xml.Characters.sankoff then
        nonaddcs_to_formater (Xml.Characters.nonadditive, attr,cont)
    else raise (Illegal_formater "sankcs_to_formater")

let seq_to_formater ((tag, attr, cont) : Xml.xml) : Xml.unstructured array  =
    if tag = Xml.Characters.sequence then begin
        let name, cclass, cost = get_name_class_and_cost attr in
        let chrom_ref = `String "-" in 
        let recost = `Float 0. in 
        let map = chrom_ref in
        match cont with
            | #Xml.unstructured as v ->
                let v = `Fun (fun () -> StatusCommon.escape (Xml.value_to_string v)) in
                [|name; cclass; cost; `String "-,-"; recost; chrom_ref; map; v|]
            | #Xml.structured_xml -> raise (Illegal_formater "seq_to_formater")
    end else raise (Illegal_formater ("seq_to_formater"))

let breakinv_to_formater ((tag, attr, cont) : Xml.xml) =
    if tag = Xml.Characters.breakinv then begin
        let name, cclass, cost = get_name_class_and_cost attr in
        let breakinv_ref = `String "-" in 
        let recost = get_recost attr in
        let map = breakinv_ref in
        let cont = match cont with 
            | #Xml.structured as cont ->
                let cont = Xml.eagerly_compute cont in
                (match List.rev (Sexpr.to_list cont) with
                    | (_, _, cont) :: _ -> cont
                    | [] -> `Empty)
            | _  -> cont
        in 
        match cont with
            | #Xml.unstructured as v ->
                let v = `Fun (fun () -> StatusCommon.escape (Xml.value_to_string v)) in
                [|name; cclass; cost; `String "-,-"; recost; breakinv_ref; map; v|]
            | #Xml.structured_xml ->
                raise (Illegal_formater "breakinv_to_formater")
    end else raise (Illegal_formater ("breakinv_to_formater"))

let chrom_to_formater ((tag, attr, cont) : Xml.xml) =
    if tag = Xml.Characters.chromosome then begin
        let name, cclass, cost = get_name_class_and_cost attr in
        let recost = get_recost attr in 
        let chrom_ref = get_ref_code attr in 
        let map = get_map attr in 
        let cont = match cont with 
            | #Xml.structured as cont -> 
                let cont = Xml.eagerly_compute cont in
                (match List.rev (Sexpr.to_list cont) with
                    | (_, _, h) :: _ -> h
                    | [] -> `Empty)
            | _ -> cont
        in 
        match cont with
            | #Xml.unstructured as v -> 
                let v = `Fun (fun () -> StatusCommon.escape (Xml.value_to_string v)) in
                [|name; cclass; cost; `String "-,-"; recost; chrom_ref; map; v|]
            | #Xml.structured_xml ->
                raise (Illegal_formater "chrom_to_formater")
    end else raise (Illegal_formater ("chrom_to_formater"))

let genome_to_formater ((tag, attr, cont) : Xml.xml) =
    if tag = Xml.Characters.genome then begin
        let name, cclass, cost = get_name_class_and_cost attr in
        let recost = get_recost attr in 
        let genome_ref = get_ref_code attr in 
        let map = get_map attr in 
        let cont = match cont with 
            | #Xml.structured as cont -> 
                let cont = Xml.eagerly_compute cont in
                (match (List.rev (Sexpr.to_list cont))  with
                    | (_, _, h) :: _ -> h
                    | [] -> `Empty)
            | _ -> cont
        in 
        match cont with
            | #Xml.unstructured as v ->
                let v = `Fun (fun () -> StatusCommon.escape (Xml.value_to_string v)) in
                [|name; cclass; cost; `String "-,-"; recost; genome_ref; map; v|]
            | #Xml.structured_xml ->
                raise (Illegal_formater "genome_to_formater")
    end else raise (Illegal_formater ("genome_to_formater"))

let annchrom_to_formater ((tag, attr, cont) : Xml.xml) =
    if tag = Xml.Characters.annchrom then begin
        let name, cclass, cost = get_name_class_and_cost attr in
        let recost = get_recost attr in 
        let chrom_ref = get_ref_code attr in
        let map = get_map attr in 
        let cont = match cont with 
            | #Xml.structured as cont -> 
                let cont = Xml.eagerly_compute cont in
                (match (List.rev (Sexpr.to_list cont))  with
                    | (_, _, h) :: _ -> h
                    | [] -> `Empty)
            | _ -> cont
        in 
        match cont with
            | #Xml.unstructured as v ->
                let v = `Fun (fun () -> StatusCommon.escape (Xml.value_to_string v)) in
                [|name; cclass; cost; `String "-,-"; recost; chrom_ref; map; v|]
            | #Xml.structured_xml ->
                raise (Illegal_formater "annchrom_to_formater")
    end else raise (Illegal_formater ("annchrom_to_formater"))

let likelihood_to_formater ((tag, attr, cont) as xml: Xml.xml) =
    if tag = Xml.Characters.likelihood then begin
        (* Basic data *)
        let name    = Xml.find_attr xml Xml.Characters.name
        and cost    = Xml.find_attr xml Xml.Characters.cost
        (* empty data *)
        and cclass  = `String Xml.Nodes.preliminary
        and recost  = `Float 0.
        and chrom   = `String "-"
        and map     = `String "-"
        and branch1 = Xml.value_to_string (Xml.find_attr xml Xml.Nodes.min_time)
        and branch2 = Xml.value_to_string (Xml.find_attr xml Xml.Nodes.oth_time)
        (* state data --since we have a complex vector; leave blank *)
        and v       = `String "-" in
        let branches = `String ((branch1)^","^(branch2)) in
        [| name;cclass; cost; branches; recost; chrom; map; v |]
    end else
        raise (Illegal_formater "likelihood_to_formater")

let dyn_likelihood_to_formater ((tag, attr, cont) as xml: Xml.xml) = 
    if tag = Xml.Characters.dlikelihood then begin
        let ((_,_,c) as c_xml) = Xml.find_tag xml Xml.Characters.characters in
        (* Basic data *)
        let name    = Xml.find_attr c_xml Xml.Characters.name
        and cclass  = Xml.find_attr xml Xml.Characters.cclass
        and cost    = Xml.find_attr xml Xml.Characters.cost
        (* empty data *)
        and recost  = `Float 0.
        and chrom   = `String "-"
        and map     = `String "-"
        and branch1 = Xml.value_to_string (Xml.find_attr xml Xml.Nodes.min_time)
        and branch2 = Xml.value_to_string (Xml.find_attr xml Xml.Nodes.oth_time)
        (* state data --pull from sequence data *)
        and v = match c with
            | #Xml.unstructured as c ->
                `Fun (fun () -> StatusCommon.escape (Xml.value_to_string c))
            | _ -> failwith "Unexpected format in Dynamic Likelihood output"
        in
        let branches = `String ((branch1)^","^(branch2)) in
        [| name;cclass; cost; branches; recost; chrom; map; v |]
    end else
        raise (Illegal_formater "dyn_likelihood_to_formater")

let node_character_to_formater ((tag, _, _) as v) =
         if tag = Xml.Characters.sequence     then seq_to_formater v
    else if tag = Xml.Characters.chromosome   then chrom_to_formater v
    else if tag = Xml.Characters.genome       then genome_to_formater v
    else if tag = Xml.Characters.breakinv     then breakinv_to_formater v
    else if tag = Xml.Characters.annchrom     then annchrom_to_formater v
    else if tag = Xml.Characters.sankoff      then sankcs_to_formater v
    else if tag = Xml.Characters.nonadditive  then nonaddcs_to_formater v
    else if tag = Xml.Characters.additive     then addcs_to_formater v
    else if tag = Xml.Characters.likelihood   then likelihood_to_formater v
    else if tag = Xml.Characters.dlikelihood  then dyn_likelihood_to_formater v
    else raise (Illegal_formater ("node_character_to_formater: " ^ tag) )

let node_to_formater st (tag, attr, cont) =
    let assoc x y = Xml.value_to_string (List.assoc x y) in
    if tag = Xml.Nodes.node then begin
        let name = assoc Xml.Characters.name attr
        and child1_name = assoc Xml.Nodes.child1_name attr
        and child2_name = assoc Xml.Nodes.child2_name attr in
     (*  we gonna get cost and recost out of the lst *)
        let cost =
            try assoc Xml.Characters.cost attr 
            with | Not_found -> " "
        in
        let recost =
            try assoc Xml.Characters.recost attr 
            with | Not_found -> " "
        in
        let lst = match cont with
            | #Xml.structured as x -> Sexpr.to_list (Xml.eagerly_compute x)
            | _ -> raise (Illegal_formater "node_to_formater 2")
        in
        let lst = List.map node_character_to_formater lst in
        let lst =
            let arr = [|"@{<u>Characters@}"; "@{<u>Class@}"; "@{<u>Cost@}";
                        "@{<u>Child Branch Lengths@}"; "@{<u>Rearrangement Cost@}";
                        "@{<u>Chrom Ref@}"; "@{<u>Median Map@}"; "@{<u>States@}"|]
            in
            arr :: (List.map (Array.map Xml.value_to_string) lst)
        in
        user_messagef st "@\n@\n@[<v 0>@{<b>%s@}@\n" name;
        user_messagef st "@[<v 0>@{<u>Cost %s@}@\n" cost;
        user_messagef st "@[<v 0>@{<u>Rearrangement cost %s@}@\n" recost;
        user_messagef st "@[<v 0>@{<u>Children: %s %s@}@\n" child1_name child2_name;
        Status.output_table st (Array.of_list lst);
        user_messagef st "@\n";
        user_messagef st "@]@]";
    end else raise (Illegal_formater "node_to_formater")

let forest_to_formater st ((tag, attr, cont) as v) =
    let find x y = 
        try Xml.value_to_string (List.assoc x y) with Not_found -> " " 
    in
    let cost = find Xml.Characters.cost attr in
    let recost = find Xml.Characters.recost attr in 
    match cont with
    | #Xml.structured ->
           (* let cost2 = ref 0.0 and recost2 = ref 0.0 in *)
            let rec do_nodes ((tag, _, cont) as x) =
                if tag = Xml.Nodes.node then 
                    node_to_formater st x (* cost2 recost2 *)
                else (*what's else part here doing?*) 
                    match cont with
                    | #Xml.structured as v -> 
                            let v = Xml.eagerly_compute v in
                            Sexpr.leaf_iter do_nodes v
                    | _ -> ()
            in
      (*
            let cost,recost = 
                (string_of_float !cost), (string_of_float !recost) in
      *)
            Status.user_message st ("@[<v 4>@{<b>Tree@}@\n@{<u>Tree cost: @}");
            Status.user_message st (cost ^ "@\n");
            Status.user_message st ("@{<u>Tree rearrangement cost: @}");
            Status.user_message st (recost ^ "@\n@\n");
            do_nodes v;
            Status.user_message st ("@]");
    |  _ -> 
            raise (Illegal_formater "forest_to_formater")

let trees_to_formater st ((tag, _, _) as r) =
    if tag = Xml.Trees.forest || tag = Xml.Trees.tree then
        forest_to_formater st r
    else raise (Illegal_formater "trees_to_formater")


let trees_to_formater filename fo_ls tree = 
    StatusCommon.Files.set_margin filename 0;
    let st = Status.Output (filename, false, fo_ls) in
    if is_xml_filename filename then
        let filename = 
            match filename with | Some x -> x | None -> assert false in
        let ch = StatusCommon.Files.channel filename in
        match ch with
        | `NotCompressed ch -> Xml.to_file ch tree
        | `Zlib ch -> 
                failwith "Output fo XML in compressed format is not supported"
    else trees_to_formater st tree
