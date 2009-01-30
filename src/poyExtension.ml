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

exception Exit 

(* The necessary types to produce the tree of the parsed input. *)
open Camlp4.Sig

module Id = struct
    let name = "POYLanguage"
    let version = "0.1"
end

module POYLanguage (Syntax : Camlp4Syntax) = struct
    open Syntax

    let to_local x = List.map (fun x -> 
                let _loc = Loc.ghost in
            <:expr<`Local $str:x$>>) x

    let expr_poy = Gram.Entry.mk "expr_poy" 
    let arg_poy = Gram.Entry.mk "arg_poy"
    let single_poy = Gram.Entry.mk "single_poy"
    let xml_poy = Gram.Entry.mk "xml_poy"
    let attr_poy = Gram.Entry.mk "attr_poy"
    let root_xml = Gram.Entry.mk "root_xml"

    let rec exSem_of_list = function
        | [] -> 
                let _loc = Loc.ghost in
                <:expr<[]>>
        | [x] -> 
                let _loc = Ast.loc_of_expr x in
                <:expr< [$x$]>>
        | x :: xs ->
                let _loc = Ast.loc_of_expr x in
                <:expr< [$x$ :: $exSem_of_list xs$] >> 

    let rec exSemCom_of_list = function
        | [] -> 
                let _loc = Loc.ghost in
                <:expr<[]>>
        | [`S x] -> 
                let _loc = Ast.loc_of_expr x in
                <:expr< [$x$]>>
        | [`L x] -> 
                let _loc = Ast.loc_of_expr x in
                <:expr< $x$>>
        | (`S x) :: xs ->
                let _loc = Ast.loc_of_expr x in
                <:expr< [$x$ :: $exSemCom_of_list xs$] >> 
        | (`L x) :: xs ->
                let _loc = Ast.loc_of_expr x in
                <:expr< $x$ @ $exSemCom_of_list xs$ >> 

    let handle_optional x = 
        let _loc = Loc.ghost in
        match x with
        | Some x -> <:expr<Some $x$>>
        | None -> <:expr<None>>

    let handle_optional_lst x = 
        let _loc = Loc.ghost in
        match x with
        | Some x -> <:expr<$exSemCom_of_list x$>>
        | None -> <:expr<[]>>

    EXTEND Gram
        GLOBAL: expr_poy arg_poy single_poy xml_poy attr_poy root_xml;
        xml_poy : [
            [ "---" -> <:expr<`Empty>> ] |
            [ "-"; "--" -> <:expr<`Empty>> ] |
            [ "{"; LIDENT "single"; x = expr; "}" -> <:expr<`Single $x$>> ] |
            [ "{"; LIDENT "set"; x = expr; "}" -> <:expr<`Set $x$>> ] |
            [ "{"; LIDENT "string"; x = expr; "}" -> <:expr<`String $x$>> ] |
            [ "{"; LIDENT "int"; x = expr; "}" -> <:expr<`Int $x$>> ] |
            [ "{"; LIDENT "float"; x = expr; "}" -> <:expr<`Float $x$>> ] |
            [ "{"; LIDENT "delayed"; x = expr; "}" -> <:expr<`Delayed $x$>> ] |
            [ "-"; LIDENT "set"; x = LIST1 [ x = xml_poy -> x]; "--" -> 
                <:expr<`Set $exSem_of_list x$>> ] |
            [ "{"; x = expr; "}" -> x ] |
            [ x = root_xml -> <:expr<`Single $x$>> ]
        ];
        root_xml : [
            [ "-"; t = tag; a = OPT [ x = attr_poy -> x ]; 
                c = OPT contents; "--" ->
                    match c with
                    | None -> 
                            <:expr<($t$, $handle_optional_lst a$, `Empty)>>
                    | Some c -> 
                            <:expr<($t$, $handle_optional_lst a$, $c$)>>]
        ];
        tag : [
            [ x = LIDENT -> <:expr<$str:x$>> ] | 
            [ x = UIDENT -> <:expr<$str:x$>> ] |
            [ x = STRING -> <:expr<$str:x$>> ] |
            [ "["; x = expr; "]" -> x ]
        ];
        attr_poy : [
            [ x = LIST0 [ 
                "("; a = att_label; v = OPT ["="; v = xml_value -> v]; 
                ")" -> 
                    match v with
                    | None -> `L <:expr<$a$>>
                    | Some v -> `S <:expr< ($a$, $v$)>> ] -> x ]
        ];
        att_label : [
            [ x = LIDENT -> <:expr<$str:x$>> ] |
            [ x = UIDENT -> <:expr<$str:x$>> ] |
            [ x = STRING -> <:expr<$str:x$>> ] |
            [ "["; x = expr; "]" -> <:expr<$x$>> ]
        ];
        xml_value : [
            [ x = LIDENT -> <:expr<`String $str:x$>> ] |
            [ x = UIDENT -> <:expr<`String $str:x$>> ] |
            [ x = STRING -> <:expr<`String $str:x$>> ] |
            [ x = INT -> <:expr<`Int $int:x$>> ] |
            [ x = FLOAT -> <:expr<`Float $flo:x$>> ] |
            [ "["; LIDENT "string"; x = expr; "]" -> <:expr<`String $x$>> ] |
            [ "["; LIDENT "int"; x = expr; "]" -> <:expr<`Int $x$>> ] |
            [ "["; LIDENT "float"; x = expr; "]" -> <:expr<`Float $x$>> ] |
            [ "["; x = expr; "]" -> <:expr<$x$>> ]
        ];
        contents : [
            [ x = xml_value -> <:expr<$x$>> ] |
            [ x = LIST1 [ x = xml_poy -> x ] -> 
                match x with
                | [] -> assert false
                | [x] -> x
                | lst -> <:expr<`Set $exSem_of_list x$>> ]
        ];
        expr_poy: [ [ 
            x = LIST1 [a = transform -> `S a 
                    | a = fuse -> `S a 
                    | a = calculate_support -> `S a 
                    | a = report -> `S a 
                    | a = select -> `S a 
                    | a = perturb -> `S a 
                    | a = swap -> `S a 
                    | a = build -> `S a 
                    | a = application_command -> `S a 
                    | x = read -> `S x 
                    | x = rename -> `S x
                    | x = search -> `S x
                    | x = simp_expr -> `S x
                    | x = cur_expr -> `L x 
                    | x = plugin -> `S x
            ]
        -> x] ];
        single_poy : [
            [ a = transform -> a ] 
            | [ a = fuse -> a ]
            | [ a = calculate_support -> a ]
            | [ a = report -> a ]
            | [ a = select -> a ]
            | [ a = perturb -> a ]
            | [ a = swap -> a ]
            | [ a = build -> a ]
            | [ a = application_command -> a ]
            | [ x = read -> x ]
            | [ x = rename -> x ]
            | [ x = search -> x ]
            | [ x = plugin -> x ]
            | [ x = simp_expr -> x ]
            | [ x = cur_expr -> x ]
        ];
        arg_poy: [ 
            [ x = swap_argument -> x ] |
            [ x = fuse_argument -> x ] |
            [ x = build_argument -> x ] |
            [ x = select_argument -> x ] |
            [ x = transform_argument -> x ] |
            [ x = read_argument -> x ] | 
            [ x = perturb_argument -> x ] |
            [ x = support_argument -> x ] |
            [ x = std_search_argument -> x ] 
        ];
        (* Application commands *)
        plugin: 
            [
                [ x = LIDENT; left_parenthesis; 
                    y = LIST0 [ x = plugin_args -> x] SEP ",";
                        right_parenthesis -> 
                    match y with
                    | [] -> <:expr<`Plugin ($str:x$, `Empty)>>
                    | [y] -> <:expr<`Plugin ($str:x$, $y$)>>
                    | y -> <:expr<`Plugin ($str:x$, `List $exSem_of_list y$)>> ]
            ];
        plugin_args:
            [   
                [ x = FLOAT -> <:expr<`Float ($flo:x$)>> ] |
                [ x = INT -> <:expr<`Int ($int:x$)>> ] |
                [ x = STRING -> <:expr<`String $str:x$>> ] |
                [ x = LIDENT; ":"; y = plugin_args -> 
                    <:expr<`Labled ($str:x$, $y$)>>] |
                [ x = LIDENT -> <:expr<`Lident $str:x$>> ] | 
                [ left_parenthesis; x = 
                    LIST0 [ x = plugin_args -> x] SEP ",";
                right_parenthesis -> <:expr<`List $exSem_of_list x$>> ]
            ];
        cur_expr: [[ "["; x = expr; "]" -> x ]];
        simp_expr: [[ "{"; x = expr; "}" -> x ]];
        (* Transforming taxa or characters *)
        setting:
            [
                [ LIDENT "timer"; ":"; x = flex_integer -> 
                    <:expr<`TimerInterval $x$>> ] |
                [ LIDENT "history"; ":"; x = flex_integer -> 
                    <:expr<`HistorySize $x$>> ] |
                [ LIDENT "log"; ":"; x = flex_string -> 
                    <:expr<`Logfile (Some $x$)>> ] |
                [ LIDENT "nolog" -> <:expr<`Logfile None>> ] |
                [ LIDENT "seed"; ":"; x = flex_integer -> 
                    <:expr<`SetSeed $x$>> ] |
                [ LIDENT "root"; ":"; x = flex_string -> 
                    <:expr<`RootName $x$>> ] |
                [ LIDENT "exhaustive_do" -> <:expr<`Exhaustive_Weak>> ] |
                [ LIDENT "iterative"; ":"; x = iterative_mode -> 
                    <:expr<`Iterative $x$>> ] |
                [ LIDENT "normal_do_plus" -> <:expr<`Normal_plus_Vitamines>> ] |
                [ LIDENT "normal_do" -> <:expr<`Normal>> ]
            ];
        iterative_mode:
            [
                [ LIDENT "exact" -> <:expr<`ThreeD None>> ] |
                [ LIDENT "approximate" -> <:expr<`ApproxD None>> ]
            ];
        application_command:
            [
                [ LIDENT "version"; left_parenthesis; right_parenthesis ->
                    <:expr<`Version>> ] |
                [ LIDENT "exit"; left_parenthesis; right_parenthesis ->
                    <:expr<`Exit>> ] |
                [ LIDENT "recover"; left_parenthesis; right_parenthesis -> 
                    <:expr<`Recover>> ] |
                [ LIDENT "clear_recovered"; left_parenthesis; 
                    right_parenthesis -> <:expr<`ClearRecovered>> ] |
                [ LIDENT "quit" ; left_parenthesis; right_parenthesis ->
                    <:expr<`Exit>> ] |
                [ LIDENT "echo"; left_parenthesis; a = flex_string; OPT ",";
                    x = LIST0 [x = output_class -> x ] SEP ","; 
                    right_parenthesis -> 
                        <:expr<`Echo ($a$, $exSem_of_list x$)>> ] |
                [ LIDENT "help"; left_parenthesis; a = OPT optional_string; 
                    right_parenthesis -> <:expr<`Help $handle_optional a$>> ] |
                [ LIDENT "set"; left_parenthesis; 
                    b = LIST0 [x = setting -> x] SEP ","; 
                    right_parenthesis -> <:expr<`Set ($exSem_of_list b$)>> ] |
                [ LIDENT "redraw"; left_parenthesis; right_parenthesis ->
                    <:expr<`Redraw>> ] |
                [ LIDENT "wipe"; left_parenthesis; right_parenthesis ->
                    <:expr<`Wipe>> ] |
                [ LIDENT "clear_memory"; left_parenthesis; x = LIST0 [x =
                    clear_options -> x]
                    SEP ","; right_parenthesis -> <:expr<`ClearMemory
                    $exSem_of_list x$>> ] |
                [ LIDENT "load"; left_parenthesis; a = flex_string; 
                    right_parenthesis -> <:expr<`Load $a$>> ] |
                [ LIDENT "save"; left_parenthesis; a = flex_string; 
                    b = OPT file_comment; right_parenthesis -> <:expr<`Save
                    ($a$, $handle_optional b$)>> ] |
                [ LIDENT "inspect"; left_parenthesis; a = flex_string; 
                    right_parenthesis -> <:expr<`InspectFile $a$>> ] |
                [ LIDENT "rediagnose"; left_parenthesis; 
                    right_parenthesis -> <:expr<`ReDiagnose>> ] |
                [ LIDENT "run"; left_parenthesis; 
                    a = LIST0 [x = flex_string -> x] SEP ","; 
                    right_parenthesis -> 
                        <:expr<`ReadScript $exSem_of_list a$>> ] |
                [ LIDENT "cd"; left_parenthesis; a = flex_string; 
                    right_parenthesis ->
                        <:expr<`ChangeWDir $a$>> ] |
                [ LIDENT "pwd"; left_parenthesis; right_parenthesis -> 
                    <:expr<`PrintWDir>> ]
            ];
        clear_options:
            [
                [ LIDENT "m" -> <:expr<`Matrices>> ] |
                [ LIDENT "s" -> <:expr<`SequencePool>> ]
            ];
        file_comment:
            [
                [ ","; x = flex_string -> x ]
            ];
        output_class:
            [
                [ LIDENT "info" -> <:expr<`Information>> ] |
                [ LIDENT "error" -> <:expr<`Error>> ] |
                [ LIDENT "output"; x = OPT optional_string -> <:expr<`Output
                $handle_optional x$>> ]
            ];
        read:
            [
                [ LIDENT "read"; "{"; x = expr; "}" -> <:expr< `Read $x$>> ] | 
                [ LIDENT "read"; left_parenthesis; a = LIST0 [x = read_argument
                -> x] SEP ","; 
                    right_parenthesis -> <:expr<`Read $exSem_of_list a$>> ] 
            ];
        read_argument:
            [ 
                [ "{"; x = expr; "}" -> x ] |
                [ LIDENT "annotated"; ":"; left_parenthesis; a = LIST1 [x =
                    otherfiles -> x] SEP ","; 
                    right_parenthesis -> <:expr<`AnnotatedFiles $exSem_of_list
                    a$>> ] |
                [ LIDENT "prealigned"; ":"; left_parenthesis; a = otherfiles;
                ","; b = prealigned_costs; ","; LIDENT "gap_opening"; c =
                    optional_integer; 
                right_parenthesis ->
                    <:expr<`Prealigned ($a$, $b$, $c$)>> ] |
                [ x = otherfiles -> x ]
            ];
        prealigned_costs:
            [
                [ LIDENT "tcm"; ":";  x = STRING ->
                    <:expr<(`Assign_Transformation_Cost_Matrix (`Local
                    $str:x$))>> ] |
                [ LIDENT "tcm"; ":"; left_parenthesis; x = INT; ","; y = INT; 
                    right_parenthesis -> 
                        <:expr<`Create_Transformation_Cost_Matrix ($int:x$,
                        $int:y$)>> ] |
                [ LIDENT "tcm"; ":"; left_parenthesis; x = cur_expr; ","; y = cur_expr; 
                    right_parenthesis -> 
                        <:expr<`Create_Transformation_Cost_Matrix ($x$,
                        $y$)>> ] | 
                [ LIDENT "tcm"; ":";  x = cur_expr ->
                    <:expr<(`Assign_Transformation_Cost_Matrix (`Local
                    $x$))>> ]
            ];
        otherfiles:
            [
                [ f = STRING -> <:expr<`AutoDetect $exSem_of_list (to_local [f])$>> ] |
                [ f = cur_expr -> <:expr<`AutoDetect [`Local $f$]>> ] |
                [ LIDENT "partitioned"; ":"; left_parenthesis; a = LIST1 [x =
                    STRING -> x] SEP ","; 
                    right_parenthesis -> <:expr<`PartitionedFile $exSem_of_list (to_local
                    a)$>> ] |
                [ LIDENT "nucleotides"; ":"; left_parenthesis; a = LIST1 [x =
                    STRING -> x] SEP ","; 
                    right_parenthesis -> <:expr<`Nucleotides $exSem_of_list (to_local
                    a)$>> ] |
                [ LIDENT "nucleotides"; ":"; left_parenthesis; a = LIST1 [x =
                    expr -> x] SEP ","; 
                    right_parenthesis -> <:expr<`Nucleotides $exSem_of_list
                    (List.map (fun x -> <:expr<`Local $x$>>) a)$>> ] |
                [ LIDENT "chromosome"; ":"; left_parenthesis; a = LIST1 [x =
                    STRING -> x] SEP ","; 
                    right_parenthesis -> <:expr<`Chromosome $exSem_of_list (to_local a)$>> ] |
                [ LIDENT "genome"; ":"; left_parenthesis; a = LIST1 [x = STRING
                -> x]SEP ","; 
                    right_parenthesis -> <:expr<`Genome $exSem_of_list (to_local a)$>> ] |
                [ LIDENT "aminoacids"; ":"; left_parenthesis; a = LIST1 [x =
                    STRING -> x] SEP ","; 
                    right_parenthesis -> <:expr<`Aminoacids $exSem_of_list
                    (to_local a)$>> ] |
                [ LIDENT "breakinv"; ":"; left_parenthesis; seq = STRING; ","; cost_mat = STRING; OPT ",";
                  read_options = LIST0 [x = read_optiona -> x] SEP ","; right_parenthesis 
                      -> <:expr<`Breakinv (`Local $str:seq$, `Local
                      $str:cost_mat$, $exSem_of_list read_options$)>>  ] 
            ];
        read_optiona:
            [
                [LIDENT "init3D"; ":"; init3D = STRING -> <:expr<`Init3D
                $str:init3D$>>] |
                [LIDENT "orientation"; ":"; ori = STRING -> 
                        <:expr<`Orientation $str:ori$>>] |
                [LIDENT "init3D"; ":"; init3D = cur_expr -> <:expr<`Init3D
                $init3D$>>] |
                [LIDENT "orientation"; ":"; ori = cur_expr -> 
                        <:expr<`Orientation $ori$>>] 
            ];
        rename:
            [
                [ rename_cmd; left_parenthesis; x = LIST0 [x = rename_argument
                -> x] SEP ","; 
                    right_parenthesis -> <:expr<`Rename $exSem_of_list x$>> ]
            ];
        rename_cmd:
            [
                [ LIDENT "synonymize" ] |
                [ LIDENT "rename" ]
            ];
        charortax:
            [
                [ LIDENT "characters" -> <:expr<`Characters>> ] |
                [ LIDENT "terminals" -> <:expr<`Taxa>> ]
            ];
        rename_argument:
            [
                [ "{"; x = expr; "}" -> x ] |
                [ x = STRING -> <:expr<`File $str:x$>> ] | 
                [ x = charortax -> <:expr<$x$>>  ] |
                [ left_parenthesis; a = STRING; ","; b = STRING; right_parenthesis -> 
                    <:expr<`Syn ($str:a$, $str:b$)>> ] |
                [ x = cur_expr -> <:expr<`File $x$>> ] | 
                [ left_parenthesis; a = cur_expr; ","; b = cur_expr; right_parenthesis -> 
                    <:expr<`Syn ($a$, $b$)>> ]
            ];
        build:
            [
                [ LIDENT "build"; "{"; a = expr; "}" -> <:expr<`Build $a$>> ] |
                [ LIDENT "build"; left_parenthesis; a = LIST0 [x =
                    build_argument -> x] SEP ","; 
                    right_parenthesis -> <:expr<`Build $exSem_of_list a$ >>]
            ];
        join_method:
            [
                [ LIDENT "sectorial"; x = OPT optional_integer -> 
                    <:expr<`UnionBased $handle_optional x$>> ]|
                [ LIDENT "all"; x = OPT optional_integer -> <:expr<`AllBased
                $handle_optional x$>> ] |
                [ LIDENT "constraint_p"; ":"; "("; 
                    x = LIST1 [x = constraint_options -> x] SEP ","; right_parenthesis
                    -> <:expr<`Partition $exSem_of_list x$>> ] |
                [ LIDENT "constraint_p"; ":"; x = flex_integer ->
                    <:expr<`Partition [`MaxDepth $x$]>> ] |
                [ LIDENT "sets"; ":"; x = cur_expr -> <:expr<`Partition [`Sets $x$]>> ] |
                [ LIDENT "constraint_p" -> <:expr<`Partition []>> ]
            ];
        build_argument:
            [
                [ "{"; x = cur_expr; "}" -> x ] |
                [ x = threshold_and_trees -> x ] |
                [ x = build_method -> x ] |
                [ x = join_method -> x ] |
                [ x = keep_method -> x ] |
                [ x = cost_calculation -> x ] |
                [ LIDENT "lookahead"; ":"; x = flex_integer -> 
                    <:expr<`Lookahead $x$>> ]
            ];
        flex_float:
            [
                [ x = INT -> <:expr<$flo:x$>> ] |
                [ x = FLOAT -> <:expr<$flo:x$>> ] |
                [ x = cur_expr -> x ] 
            ];
        flex_integer:
            [ [x = INT -> <:expr<$int:x$>> ] |
            [ x = cur_expr -> x ] ];
        flex_string:
            [ [x = STRING -> <:expr<$str:x$>> ] |
            [ x = cur_expr -> x ] ];
        optional_integer_or_float:
            [ [ ":"; x = flex_float -> x ] ];
        optional_string:
            [ [ ":"; x = flex_string -> x ] ];
        optional_integer:
            [ [ ":"; x = flex_integer -> x ] ];
        threshold_and_trees:
            [
                [ LIDENT "threshold"; ":"; x = flex_float -> 
                    <:expr<`Threshold $x$>> ] |
                [ LIDENT "trees"; ":"; x = flex_integer -> <:expr<`Trees $x$>> ] |
                [ x = flex_integer -> <:expr<`Trees $x$>> ]
            ];
        build_method:
            [
                [ x = flex_string -> <:expr<`Prebuilt (`Local $x$)>> ] |
                [ LIDENT "of_file"; ":"; x = flex_string -> <:expr<`Prebuilt
                (`Local $x$)>> ] |
                [ LIDENT "randomized" -> <:expr<`Random>> ] |
                [ LIDENT "random" -> <:expr<`RandomTree>> ] |
                [ LIDENT "as_is" -> <:expr<`Ordered>> ] |
                [ LIDENT "branch_and_bound"; x = OPT optional_integer_or_float -> 
                    let thresh = handle_optional x in
                    <:expr<`Branch_and_Bound $thresh$>> ] |
                [ LIDENT "constraint_p"; x = OPT optional_string -> 
                    <:expr<`Constraint $handle_optional x$>> ] |
                [ LIDENT "nj" -> <:expr<`Nj>> ] | 
                [ LIDENT "_mst" -> <:expr<`Mst>> ] |
                [ LIDENT "_distances" -> <:expr<`DistancesRnd>> ]
            ];
        keep_method:
            [ 
                [ LIDENT "last" -> <:expr<`Last>> ] |
                [ LIDENT "first" -> <:expr<`First>> ] | 
                [ LIDENT "at_random" -> <:expr<`Keep_Random>> ] 
            ];
        perturb:
            [
                [ LIDENT "perturb"; "{"; x = expr; "}" -> <:expr<`Perturb $x$>>]
                | [ LIDENT "perturb"; left_parenthesis; x = LIST0 [x =
                    perturb_argument -> x] SEP ","; 
                    right_parenthesis -> <:expr<`Perturb $exSem_of_list x$>> ]
            ];
        perturb_argument:
            [
                [ "{"; x = expr; "}" -> x ] |
                [ x = ratchet -> x ]  |
                [ x = resample -> x ] |
                [ x = swap -> x ] |
                [ x = transform -> x ] |
                [ LIDENT "iterations"; ":"; x = flex_integer -> 
                    <:expr<`Iterations $x$>> ] |
                [ LIDENT "timeout"; ":"; x = flex_float -> 
                    <:expr<`TimeOut $x$>> ]
            ];
        ratchet:
            [
                [ LIDENT "ratchet"; ":"; left_parenthesis; x = flex_float; ",";
                y = flex_integer; right_parenthesis -> <:expr<`Ratchet ($x$, $y$)>> ] |
                [ LIDENT "ratchet" -> <:expr<perturb_default_perturb>> ]
            ];
        resample:
            [
                [ LIDENT "resample"; ":"; left_parenthesis; x = flex_integer; ","; y = charortax; 
                    right_parenthesis -> <:expr<`Resample ($x$, $y$)>> ]
            ];
        report:
            [
                [ LIDENT "report"; "{"; x = expr; "}" -> <:expr<`Report $x$>> ] |
                [ LIDENT "report"; left_parenthesis; a = LIST0 [x =
                    report_argument -> x] SEP ","; 
                    right_parenthesis -> <:expr<`Report $exSem_of_list a$>> ]
            ];
        boolean: 
            [
                [ x = cur_expr -> <:expr<$x$>> ] |
                [ LIDENT "true" -> <:expr<$str:"True"$>> ] |   
                [ LIDENT "false" -> <:expr<$str:"False"$>> ]    
            ];
        optional_boolean:
            [
                [ ":"; x = boolean -> x ]
            ];
        collapse:
            [
                [ LIDENT "collapse"; y = OPT optional_boolean ->
                    <:expr<`Collapse $handle_optional y$>> ]
            ];
        optional_collapse:
            [ 
                [ ":"; x = collapse -> x ]
            ];
        tree_information_list:
            [   
                [ ":"; "("; x = LIST0 [x = tree_information -> x] SEP ","; ")"
                -> exSem_of_list x ]
            ];
        tree_information:
            [
                [ LIDENT "_cost" -> <:expr<`Cost>> ] |
                [ LIDENT "hennig" -> <:expr<`HennigStyle>> ] |
                [ LIDENT "total" -> <:expr<`Total>> ] |
                [ LIDENT "newick" -> <:expr<`Newick>> ] |
                [ LIDENT "margin"; ":"; m = flex_integer -> <:expr<`Margin $m$>> ] |
                [ LIDENT "nomargin" -> <:expr<`Margin (1000000010 - 1)>> ] |
                [ x = collapse -> x  ]
            ];
        optional_old_identifiers:
            [ [":"; x = old_identifiers -> x] ];
        opt_support_names:
            [
                [ ":"; y = support_names -> y ]
            ];
        summary_class:
            [ 
                [ ":"; LIDENT "individual" -> <:expr<`Individual>> ] | 
                [ ":"; LIDENT "consensus" -> <:expr<`Consensus>> ] |
                [ ":"; x = flex_string -> <:expr<`InputFile $x$ >> ]
            ];
        support_names:
            [
                [ LIDENT "bremer"; ":"; x = flex_string -> <:expr<`Bremer (Some
            [(`Local $x$)])>>] |
                [ LIDENT "bremer" -> <:expr<`Bremer None>> ] |
                [ LIDENT "jackknife"; y = OPT [ x = summary_class -> x ] -> 
                    match y with
                    | None -> <:expr<`Jackknife `Individual>>
                    | Some y -> <:expr<`Jackknife $y$>>] |
                [ LIDENT "bootstrap"; y = OPT [ x = summary_class -> x ] -> 
                    match y with
                    | None -> <:expr<`Bootstrap `Individual>>
                    | Some y -> <:expr<`Bootstrap $y$>>]
            ];
        report_argument:
            [
                [ "{"; x = expr; "}" -> x ] |
                [ x = flex_string -> <:expr<`File $x$>> ] |
                [ LIDENT "kml"; ":"; left_parenthesis; 
                    plugin = LIDENT; ","; csv = flex_string;
                right_parenthesis -> 
                    <:expr<`KML ($str:plugin$, $csv$)>> ] |
                [ LIDENT "asciitrees" ; y = OPT optional_collapse -> 
                    match y with
                    | None -> <:expr<`Ascii False>> 
                    | Some x -> <:expr<`Ascii $x$>> ] |
                [ LIDENT "memory" -> <:expr<`Memory>> ] | 
                [ LIDENT "graphtrees" -> <:expr<`Graph True>> ] |
                [ LIDENT "trees"; x = OPT tree_information_list -> 
                    match x with
                    | Some x -> <:expr<`Trees $x$>> | None -> <:expr<`Trees []>> ] |
                [ LIDENT "treestats" -> <:expr<`TreesStats>> ] |
                [ LIDENT "treecosts" -> <:expr<`TreeCosts>> ] |
                [ LIDENT "timer"; ":"; x = flex_string -> <:expr<`TimeDelta $x$>>] |
                [ LIDENT "_mst" -> <:expr<`MstR>> ] | 
                [ LIDENT "consensus"; x = OPT optional_integer_or_float -> 
                    <:expr<`Consensus $handle_optional x$>> ] | 
                [ LIDENT "graphconsensus"; x = OPT optional_integer_or_float -> 
                    <:expr<`GraphicConsensus $handle_optional x$>> ] | 
                [ LIDENT "clades" -> <:expr<`Clades>> ] |
                [ LIDENT "phastwinclad" -> <:expr<`FasWinClad>> ] | 
                [ LIDENT "seq_stats"; ":"; ch = old_identifiers ->
                    <:expr<`SequenceStats $ch$>> ] |
                [ LIDENT "ci"; x = OPT optional_old_identifiers -> <:expr<`Ci
                $handle_optional x$>> ] |
                [ LIDENT "ri"; x = OPT optional_old_identifiers -> <:expr<`Ri
                $handle_optional x$>> ] |
                [ LIDENT "compare"; ":"; left_parenthesis; complement = boolean;
                ","; ch1 = old_identifiers; ","; ch2 = old_identifiers; right_parenthesis ->
                    <:expr<`CompareSequences ($complement$, $ch1$, $ch2$)>> ] |
                [ LIDENT "script_analysis"; ":"; x = flex_string ->
                    <:expr<`ExplainScript $x$>> ] |
                [ LIDENT "supports"; y = OPT opt_support_names ->
                    <:expr<`Supports $handle_optional y$>> ] |
                [ LIDENT "graphsupports"; y = OPT opt_support_names -> 
                    <:expr<`GraphicSupports $handle_optional y$>> ] |
                [ LIDENT "diagnosis" -> <:expr<`Diagnosis>> ] |
                [ LIDENT "data" -> <:expr<`Data >>] |
                [ LIDENT "xslt"; ":"; "("; a = flex_string; ","; b = flex_string; ")" ->
                    <:expr<`Xslt ($a$, $b$)>> ] |
                [ LIDENT "implied_alignments"; ":"; x = identifiers ->
                    <:expr<`Implied_Alignments ($x$, True)>> ] |
                [ LIDENT "fasta"; ":"; x = identifiers ->
                    <:expr<`Implied_Alignments ($x$, False)>> ] |
                [ LIDENT "all_roots" -> <:expr<`AllRootsCost >>] |
                [ LIDENT "implied_alignments" -> 
                    <:expr<`Implied_Alignments (`All, True)>>] |
                [ LIDENT "ia"; ":"; x = identifiers -> 
                    <:expr<`Implied_Alignments ($x$, True)>> ] | 
                [ LIDENT "ia" -> <:expr<`Implied_Alignments (`All, True)>> ] |
                [ LIDENT "nodes" -> <:expr<`Nodes >>] |
                [ LIDENT "cross_references"; ":"; x = old_identifiers -> 
                    <:expr<`CrossReferences (Some $x$)>> ] |
                [ LIDENT "terminals" -> <:expr<`TerminalsFiles >>] | 
                [ LIDENT "cross_references" -> <:expr<`CrossReferences None >>]
            ];
        transform:
            [
                [ LIDENT "transform"; "{"; x = expr; "}" -> <:expr<`Transform
                $x$>> ] |
                [ LIDENT "transform"; left_parenthesis; 
                    x = LIST0 [ x = transform_argument -> x] SEP ","; right_parenthesis ->
                        <:expr<(`Transform $exSem_of_list x$ )>> ] 
            ];
        transform_argument:
            [
                [ "{"; x = expr; "}" -> x ] |
                [ left_parenthesis; x = identifiers; ","; t = transform_method; 
                    right_parenthesis -> <:expr<(x, t)>> ] |
                [ t = transform_method -> <:expr<(`All, $t$)>> ]
            |   [ LIDENT "origin_cost"; ":"; x = flex_float
                        -> <:expr<(`All, `OriginCost x)>> ] |
                [ LIDENT "prioritize" -> <:expr<(`All, `Prioritize)>> ] 
            ];
        transform_method:
            [
                [ LIDENT "prealigned" -> <:expr<`Prealigned_Transform >>] |
                [ LIDENT "randomize_terminals" -> <:expr<`RandomizedTerminals >>] |
                [ LIDENT "alphabetic_terminals" -> <:expr<`AlphabeticTerminals
                >>] |
                [ LIDENT "tcm"; ":"; "("; x = flex_integer; ","; y
                = flex_integer; ")" -> <:expr<`Gap ($x$, $y$)>> ] |
                [ LIDENT "tcm"; ":";  x = flex_string -> <:expr<`Tcm $x$>> ] |
                [ LIDENT "partitioned" -> <:expr<`Partitioned>> ] |
                [ LIDENT "fixed_states" -> <:expr<`Fixed_States>> ] |
                [ LIDENT "direct_optimization" -> <:expr<`Direct_Optimization>> ] |
                [ LIDENT "gap_opening"; ":"; x = flex_integer -> <:expr<`AffGap $x$>> ] |
                [ LIDENT "static_approx"; x = OPT informative_characters -> 
                    match x with 
                    | None -> <:expr<`StaticApproximation True>>
                    | Some v -> <:expr<`StaticApproximation $v$>>] |
                [ LIDENT "multi_static_approx"; x = OPT informative_characters -> 
                    match x with 
                    | None -> <:expr<`MultiStaticApproximation True>> 
                    | Some v -> <:expr<`MultiStaticApproximation $v$ >>] |
                [ LIDENT "auto_static_approx"; x = OPT optional_boolean -> 
                    match x with
                    | None -> <:expr<`Automatic_Static_Aprox False>> 
                    | Some x -> <:expr<`Automatic_Static_Aprox $x$ >>] |
                [ LIDENT "auto_sequence_partition"; x = OPT optional_boolean -> 
                    match x with
                    | None -> <:expr<`Automatic_Sequence_Partition (False,
                    None)>>
                    | Some x -> <:expr<`Automatic_Sequence_Partition ($x$, None)>> ] |
                [ LIDENT "sequence_partition"; ":"; x = flex_integer -> 
                    <:expr<`Automatic_Sequence_Partition (False, Some $x$)>> ] |
                [ LIDENT "weight"; ":"; x = neg_integer_or_float -> 
                    <:expr<`ReWeight $x$>> ] |
                [ LIDENT "weightfactor"; ":"; x = neg_integer_or_float ->
                    <:expr<`WeightFactor $x$>> ] |
                [ LIDENT "search_based" -> <:expr<`SearchBased>> ] |
                [ LIDENT "seq_to_chrom"; ":"; left_parenthesis; x = LIST0
                        [ x = chromosome_argument -> x] SEP ",";
                        right_parenthesis -> <:expr<`SeqToChrom $exSem_of_list x$>> ] | 
                [ LIDENT "custom_to_breakinv"; ":"; left_parenthesis; x = LIST0
                        [ x = chromosome_argument -> x] SEP ","; right_parenthesis -> 
                            <:expr<`CustomToBreakinv $exSem_of_list x$>> ] | 
                [ LIDENT "annchrom_to_breakinv"; ":"; left_parenthesis; x = LIST0
                        [x = chromosome_argument -> x] SEP ","; right_parenthesis -> 
                            <:expr<`AnnchromToBreakinv $exSem_of_list x$>> ] | 
                [ LIDENT "dynamic_pam"; ":"; left_parenthesis; x = LIST0 
                        [ x = chromosome_argument -> x] SEP ","; right_parenthesis -> 
                            <:expr<`ChangeDynPam $exSem_of_list x$>> ] | 
                [ LIDENT "chrom_to_seq" -> <:expr<`ChromToSeq []>> ] |
                [ LIDENT "breakinv_to_custom" -> <:expr<`BreakinvToSeq []>> ] 
            ];
        time:
            [
                [ days = flex_float; ":"; hours = flex_float; ":";
                minutes = flex_float ->
                    <:expr<(((($days$) *. 60. *. 60. *. 24.) +.
                    (($hours$) *. 60. *. 60.) +.
                    (($minutes$) *. 60. )))>> ] 
            ];
        memory:
            [
                [ LIDENT "gb"; ":"; x = flex_float -> 
                    <:expr<int_of_float ($x$ *. 1000. *. 1000. *. 1000. /. (float_of_int
                    (Sys.word_size / 8)))>> ] |
                [ LIDENT "mb"; ":"; x = flex_float ->
                    <:expr<int_of_float ($x$ *. 1000. *. 1000. /. (float_of_int
                    (Sys.word_size / 8)))>> ]
            ];
        std_search_argument:
            [   
                [ "{"; x = expr; "}" -> x ] |
                [ LIDENT "memory"; ":"; x = memory -> <:expr<`MaxRam $x$>> ] |
                [ LIDENT "hits"; ":"; x = flex_integer -> <:expr<`MinHits $x$>> ] |
                [ LIDENT "target_cost"; ":"; x = flex_float -> <:expr<`Target $x$>> ] |
                [ LIDENT "max_time"; ":"; x = time -> <:expr<`MaxTime $x$>> ] |
                [ LIDENT "min_time"; ":"; x = time -> <:expr<`MinTime $x$>> ] |
                [ LIDENT "visited"; ":"; x = flex_string -> <:expr<`Visited
                (Some $x$)>> ] |
                [ LIDENT "visited" -> <:expr<`Visited None>> ] |
                [ LIDENT "constraint_s"; ":"; x = flex_string -> 
                    <:expr<`ConstraintFile $x$>> ]
            ];
        search:
            [
                [ LIDENT "search"; "{"; x = expr; "}" -> <:expr<`StandardSearch $x$>> ] |
                [ LIDENT "search"; left_parenthesis; a = LIST0 [ x =
                    std_search_argument -> x ] SEP ",";  right_parenthesis ->
                    <:expr<`StandardSearch $exSem_of_list a$>>] 
            ];
        neg_integer_or_float:
            [
                [ x = flex_float -> x ]
            ];
        informative_characters:
            [
                [ ":"; LIDENT "keep" -> <:expr<$str:"False"$>> ] |
                [ ":"; LIDENT "remove" -> <:expr<$str:"True"$>> ]
            ];
        chromosome_argument:
            [
                [ LIDENT "locus_inversion"; ":"; c = flex_integer -> 
                      <:expr<`Locus_Inversion $c$>> ]  |
                [ LIDENT "locus_breakpoint"; ":"; c = flex_integer -> 
                      <:expr<`Locus_Breakpoint $c$>> ]  |
                [ LIDENT "chrom_breakpoint"; ":"; c = flex_integer -> 
                      <:expr<`Chrom_Breakpoint $c$>> ]  |
                [ LIDENT "circular"; ":"; e = boolean -> <:expr<`Circular $e$>>] |
                [ LIDENT "locus_indel"; ":"; left_parenthesis; o = flex_integer; 
                    ","; e = flex_float; right_parenthesis ->
                      <:expr<`Locus_Indel_Cost ( $o$, 
                                      int_of_float ($e$ *. 100.0) )>>] | 
                [ LIDENT "chrom_indel"; ":"; left_parenthesis; o = flex_integer; 
                ","; e = flex_float; right_parenthesis ->
                      <:expr<`Chrom_Indel_Cost ( $o$, int_of_float ($e$ *.
                      100.0) )>> ] | 
                [ LIDENT "chrom_hom"; ":"; c = flex_float -> 
                      <:expr<`Chrom_Hom (int_of_float ($c$ *. 100.))>> ] | 
                [ LIDENT "min_block_len"; ":"; c = flex_integer -> 
                      <:expr<`Sig_Block_Len $c$>> ] | 
                [ LIDENT "min_rearrangement_len"; ":"; c = flex_integer -> 
                      <:expr<`Rearranged_Len $c$ >>] | 
                [ LIDENT "min_seed_length"; ":"; c = flex_integer -> 
                      <:expr<`Seed_Len $c$>> ] | 
                [ LIDENT "median"; ":"; c = flex_integer ->
                      <:expr<`Keep_Median $c$>> ] |
                [ LIDENT "swap_med"; ":"; iters = flex_integer ->
                    <:expr<`SwapMed $iters$>>] | 
                [ LIDENT "med_approx"; ":"; ans = boolean -> <:expr<`Approx $ans$>>] |
                [ LIDENT "symmetric"; ":"; ans = boolean -> <:expr<`Symmetric
                $ans$>>] |
                [ LIDENT "max_3d_len"; ":"; l = flex_integer ->
                    <:expr<`Max_3D_Len $l$>> ]  |
                [ LIDENT "max_kept_wag"; ":"; l = flex_integer ->
                    <:expr<`Max_kept_wag $l$>> ]  
            ];
        calculate_support:
            [
                [ LIDENT "calculate_support"; "{"; x = expr; "}" ->
                    <:expr<`Supports $x$>> ] |
                [ LIDENT "calculate_support"; left_parenthesis; a = LIST0
                [x = support_argument -> x]  SEP ","; 
                right_parenthesis -> <:expr<`Support $exSem_of_list a$>> ]
            ];
        fuse:
            [
                [ LIDENT "fuse"; "{"; a = expr; "}" -> <:expr<`Fuse $a$>> ] |
                [ LIDENT "fuse"; left_parenthesis; a = LIST0 [x = fuse_argument
                -> x] SEP ","; 
                    right_parenthesis -> <:expr<(`Fuse $exSem_of_list a$)>> ]
            ];
        fuse_argument:
            [
                [ "{"; x = expr; "}" -> x ] |
                [ LIDENT "keep"; ":"; i = flex_integer -> <:expr<`Keep $i$>> ]
            |   [ LIDENT "iterations"; ":"; i = flex_integer ->
                    <:expr<`Iterations $i$>> ]
            |   [ LIDENT "replace"; ":"; r = fuseareplace -> <:expr<`Replace $r$>>]
            |   [ x = swap -> x ]
            |   [ LIDENT "weighting"; ":"; w = fuseaweighting ->
                    <:expr<`Weighting $w$>>]
            |   [ LIDENT "clades"; ":"; cfrom = flex_integer; cto = OPT fusea_cto ->
                      <:expr<`Clades ($cfrom$, $handle_optional cto$)>>]
            ];
        fuseareplace:
            [ [ LIDENT "better" -> <:expr<`Better >>]
            | [ LIDENT "best" -> <:expr<`Best >>] ];
        fuseaweighting:
            [ [ LIDENT "uniform" -> <:expr<`Uniform >>] ];
        fusea_cto:
            [ [ "-"; cto = flex_integer -> cto ] ];
        support_argument:
            [
                [ "{"; x = expr; "}" -> x ] |
                [ x = build -> x ] |
                [ x = swap -> x ] |
                [ x = support_method -> x ]
            ];
        list_of_jackknifea:
            [
                [ ":"; left_parenthesis; x = LIST0 [x = jackknifea -> x] SEP ","; 
                    right_parenthesis -> exSem_of_list x ]
            ];
        jackknifea:
            [
                [ LIDENT "remove"; ":"; x = flex_float -> <:expr<`Select $x$>> ] |
                [ LIDENT "resample"; ":"; x = flex_integer -> <:expr<`Resample
                $x$>> ]
            ];
        support_method:
            [
                [ LIDENT "bremer" -> <:expr<`Bremer>> ] |
                [ LIDENT "jackknife"; x = OPT list_of_jackknifea ->
                    match x with
                    | None -> <:expr<`Jackknife []>>
                    | Some v -> <:expr<`Jackknife $v$>> ] |
                [ LIDENT "bootstrap"; x = OPT optional_integer -> <:expr<`Bootstrap
                $handle_optional x$>>]
            ];
        (* Selecting characters or taxa *)
        select:
            [
                [ LIDENT "select"; "{"; x = expr; "}" -> <:expr<`Select $x$>>]
                | [ LIDENT "select"; left_parenthesis; x = LIST0 [x =
                    select_argument -> x] SEP ","; 
                    right_parenthesis -> <:expr<`Select $exSem_of_list x$>> ]
            ];
        select_argument:
            [
                [ "{"; x = expr; "}" -> x ] |
                [ x = identifiers -> x ] |
                [ x = charortax -> x ] |
                [ x = seltrees -> x ] |
                [ x = flex_string -> <:expr<`Files (True, $exSem_of_list [x]$)>> ] 
            ];
        identifiers:
            [
                [ LIDENT "not" ; LIDENT "files"; ":"; left_parenthesis; x =
                    LIST0 [x = flex_string -> x] SEP ","; 
                    right_parenthesis -> <:expr<`Files (False, $exSem_of_list x$)>> ] |
                [ x = old_identifiers -> x ] |
                [ LIDENT "files"; ":"; left_parenthesis; x = LIST0 [x =
                    flex_string -> x] SEP ","; 
                    right_parenthesis -> <:expr<`Files (True, $exSem_of_list x$)>> ] 
            ];
        old_identifiers:
            [
                [ LIDENT "all" -> <:expr<`All>> ] |
                [ LIDENT "not"; LIDENT "names"; ":"; left_parenthesis; x = LIST0
                [x = flex_string -> x] SEP ","; 
                    right_parenthesis -> <:expr<`Names (False, $exSem_of_list x$)>> ] |
                [ LIDENT "not"; LIDENT "codes"; ":"; left_parenthesis; x = LIST0
                [x = flex_integer -> x] SEP ","; 
                    right_parenthesis -> <:expr<`Some (False, $exSem_of_list x$)>> ] |
                [ LIDENT "names"; ":"; left_parenthesis; x = LIST0 [x =
                    flex_string -> x] SEP ","; 
                    right_parenthesis -> <:expr<`Names (True, $exSem_of_list
                    x$)>> ] |
                [ LIDENT "codes"; ":"; left_parenthesis; x = LIST0 [x =
                    flex_integer ->
                    x] SEP ","; right_parenthesis -> <:expr<`Some (True,
                    $exSem_of_list x$)>> ] |
                [ LIDENT "static" -> <:expr<`AllStatic>> ] | 
                [ LIDENT "dynamic" -> <:expr<`AllDynamic >>] |
                [ LIDENT "missing"; ":"; x = flex_integer -> 
                    <:expr<`Missing (True, 100 - $x$)>> ] |
                [ LIDENT "not"; LIDENT "missing"; ":"; x = flex_integer -> 
                    <:expr<`Missing (False, 100 - $x$)>> ] |
                [ LIDENT "_random"; ":"; x = flex_float -> 
                    <:expr<`Random (100. -. $x$)>> ]
            ];
        seltrees:
            [
                [ LIDENT "optimal" -> <:expr<`BestN None>> ] |
                [ LIDENT "unique" -> <:expr<`Unique>> ] |
                [ LIDENT "best"; ":"; x = flex_integer -> 
                    <:expr<`BestN (Some $x$)>> ] |
                [ LIDENT "within"; ":"; x = flex_float -> 
                    <:expr<`BestWithin $x$>> ] |
                [ LIDENT "random"; ":"; x = flex_integer -> 
                    <:expr<`RandomTrees $x$>> ]
            ];
        swap:
            [
                [ LIDENT "swap"; "{"; x = expr; "}" -> <:expr<`Swap $x$>> ] |
                [ LIDENT "swap"; left_parenthesis; 
                    a = LIST0 [x = swap_argument -> x] SEP ","; 
                    right_parenthesis -> <:expr<`Swap $exSem_of_list a$>> ] 
            ];
        swap_argument:
            [
                [ x = sample_method -> x ] |
                [ x = swap_method -> x ] |
                [ x = keep_method -> x ] |
                [ x = threshold_and_trees -> x ] |
                [ x = cost_calculation -> x ] |
                [ LIDENT "forest"; a = OPT optional_integer_or_float -> 
                    match a with
                    | None -> <:expr<`Forest $flo:"0."$>>
                    | Some a -> <:expr<`Forest $a$>> ] |
                [ a = trajectory_method -> a ] |
                [ a = break_method -> a ] |
                [ a = reroot_method -> a ] |
                [ a = join_method -> a ] |
                [ "{"; x = expr; "}" -> x ]
            ];
        trajectory_method:
            [
                [ LIDENT "annealing"; ":"; left_parenthesis; 
                  coeff = flex_float;
                  ","; exp = flex_float; right_parenthesis ->
                      <:expr<`Annealing ($coeff$, $exp$)>>] |
                [ LIDENT "drifting"; ":"; left_parenthesis;
                    equalprob = flex_float; ","; 
                    worstfactor = flex_float; right_parenthesis ->
                        <:expr<`PoyDrifting ($equalprob$, $worstfactor$)>> ] |
                [ LIDENT "current_neighborhood"; f = OPT optional_string -> 
                    <:expr<`AllAround $handle_optional f$>> ] |
            [ LIDENT "around" -> <:expr<`AllThenChoose>> ]
            ];
        sample_method:
            [
                [ LIDENT "_maxtrees"; ":"; x = flex_integer ->
                    <:expr<`MaxTreesEvaluated $x$>> ] |
                [ LIDENT "timeout"; ":"; x = flex_float -> 
                    <:expr<`TimeOut $x$>> ] |
                [ LIDENT "timedprint"; ":"; left_parenthesis; x = flex_float; ","; 
                    y = flex_string; right_parenthesis -> 
                        <:expr<`TimedPrint ($x$, Some $y$)>> ] |
                [ LIDENT "visited"; x = OPT optional_string ->
                    <:expr<`AllVisited $handle_optional x$>> ] |
                [ LIDENT "trajectory"; x = OPT optional_string -> 
                        <:expr<`PrintTrajectory $handle_optional x$>> ] |
                [ LIDENT "recover" -> <:expr<`KeepBestTrees>> ] |
                [ LIDENT "_unionstats"; ":"; left_parenthesis; 
                    x = flex_string; ","; y = flex_integer; right_parenthesis ->
                    <:expr<`UnionStats (Some $x$, $y$)>> ] |
                [ LIDENT "_rootuniondistr"; ":"; x = flex_string -> 
                    <:expr<`RootUnionDistr (Some $x$)>> ] |
                [ LIDENT "_attemptsdistr"; ":"; x = flex_string ->
                    <:expr<`AttemptsDistr (Some $x$)>> ] |
                [ LIDENT "_breakvsjoin"; x = OPT optional_string ->
                    <:expr<`BreakVsJoin $handle_optional x$>> ]
            ];
        swap_method:
            [
                [ LIDENT "spr"; y = OPT single_option -> 
                    match y with
                    | None -> <:expr<`ChainNeighborhoods `Spr>>
                    | Some _ -> <:expr<`SingleNeighborhood `Spr>>
                ] |
                [ LIDENT "tbr"; y = OPT single_option -> 
                    match y with
                    | None -> <:expr<`ChainNeighborhoods `Tbr>> 
                    | Some _ -> <:expr<`SingleNeighborhood `Tbr>>] |
                [ LIDENT "alternate" -> <:expr<`Alternate (`Spr, `Tbr)>> ]
            ];
        single_option:
            [
                [ ":"; LIDENT "once" -> 1 ]
            ];
        break_method:
            [
                [ LIDENT "randomized" -> <:expr<`Randomized>>] |
                [ LIDENT "distance"; x = OPT optional_boolean -> 
                    match x with
                    | None -> <:expr<`DistanceSorted False>>
                    | Some x -> <:expr<`DistanceSorted $x$>>] |
                [ LIDENT "once" -> <:expr<`OnlyOnce>> ]
            ];
        constraint_options:
            [
                [ x = flex_integer -> <:expr<`MaxDepth $x$>> ] |
                [ LIDENT "depth"; ":"; x = flex_integer -> 
                    <:expr<`MaxDepth $x$>> ] |
                [ x = flex_string -> <:expr<`ConstraintFile (`Local $x$)>> ] |
                [ LIDENT "file"; ":"; x = flex_string -> <:expr<`ConstraintFile
                (`Local $x$)>>]
            ];
        reroot_method:
            [
                [ LIDENT "bfs"; x = OPT optional_integer -> <:expr<`Bfs
                $handle_optional x$>> ]
            ];
        cost_calculation:
            [
                [ x = transform -> x ]
            ];
        left_parenthesis:
            [ ["("] ];
        right_parenthesis:
            [ [")"] ];
    END;;

    EXTEND Gram
    Syntax.expr : LEVEL "top" [
        [ "PXML"; s = xml_poy -> s ] | 
        [ "PCDATA"; s = xml_poy -> <:expr<`CDATA $s$>> ] |
        [ "AXML"; s = attr_poy -> exSemCom_of_list s ] | 
        [ "RXML"; s = root_xml -> s ] | 
        [ "CPOY"; s = expr_poy -> exSemCom_of_list s ] |
        [ "SPOY"; s = single_poy -> s ] |
        [ "APOY"; s = arg_poy -> s ] |
        [ "POY"; s = expr_poy -> 
            <:expr<Phylo.parsed_run (PoyCommand.of_parsed True 
            $exSemCom_of_list s$)>> ] |
        [ "NPOY"; s = expr_poy -> 
            <:expr<Phylo.parsed_run (PoyCommand.of_parsed False
            $exSemCom_of_list s$)>> ] |
        [ "GPOY" -> <:expr<Phylo.get_console_run ()>> ] 
    ];
    END;;

    include Syntax
end

let () =
    let module M = Camlp4.Register.OCamlSyntaxExtension (Id) (POYLanguage) in
    ()
