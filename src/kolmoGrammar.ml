(* A File that defines OCaml extensions and a DSL for SK operations and
* definitions *)
open Camlp4.Sig

type primitives = 
    [ `S 
    | `K 
    | `Label of string 
    | `Node of primitives list 
    | `Debugger of string 
    | `Lazy of primitives Lazy.t ]

type 'a kolmo_function =
    [ `Module of 'a kolmo_function list
    | `LetVal of (name * arguments * 'a definition)
    | `RecVal of (name * arguments * 'a definition) 
    | `Let of (name * arguments * 'a definition)
    | `Rec of (name * arguments * 'a definition) ]
and name = string
and arguments = string list 
and 'a definition =
    [ `Value of 'a values
    | `IfElse of ('a condition * 'a definition * 'a definition)
    | `Letin of ('a kolmo_function * 'a definition) ]
and 'a condition = [ `Condition of 'a definition ]
and 'a values =
    [ `Apply of (string list * string * ('a definition list))
    | `Integer of string 
    | `Expr of 'a ]


module Id = struct
    let name = "SKLanguage"
    let version = "0.1"
end

module SKLanguage (Gram : Grammar.Static) (Syntax : Camlp4Syntax) = struct
    open Syntax

    (* An SK expression *)
    type 'loc expr = 
        | S of 'loc 
        | K of 'loc 
        | Label of string * 'loc 
        | Node of 'loc expr list * 'loc
        | Expr of Ast.expr * 'loc
        | Debg of Ast.expr * 'loc

    (* A valid statement in our extension to OCaml *)
    type 'loc stmt = 
        | Def of string * string list * 'loc expr * 'loc
        | Evaluate of 'loc expr * 'loc
        | BitEncode of 'loc expr * 'loc
        | BitDecode of string * string * 'loc

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

    let rec bare_expression_converter = function
        | S _loc -> <:expr< (`S) >>
        | K _loc -> <:expr< (`K) >>
        | Label (v, _loc) -> 
                <:expr< `Label $str:v$ >> 
        | Node (items, _loc) -> 
                let es = List.map bare_expression_converter items in
                <:expr< (`Node ( $exSem_of_list es$)) >>
        | Expr (item, _loc) ->
                <:expr< $item$ >>
        | Debg (item, _loc) ->
                <:expr< (`Debugger ( $item$ )) >> 

    let rec sk_to_ocaml x : primitives = 
        match x with
        | S _loc -> `S
        | K _loc -> `K
        | Label (v, _loc) -> `Label v
        | Expr _
        | Debg _ -> assert false 
        | Node (items, _loc) -> 
                let es = List.map sk_to_ocaml items in
                match es with
                | [es] -> es
                | _ -> `Node es


    let expression_converter _loc x =
        let res = bare_expression_converter x in
        <:expr< `Processed $res$>>

    let lst_to_expr _loc lst = 
        exSem_of_list (List.map (fun x -> <:expr<$str:x$>>) lst)

    let sklst_to_expr _loc lst =
        let rec aux_to_expr lst = 
            match lst with
            | [] -> <:expr<`Label "EmptyList">>
            | h :: t -> 
                    let at = aux_to_expr t in
                    <:expr<`Node [`Label "Pair"; $h$; $at$]>>
        in
        aux_to_expr (List.map bare_expression_converter lst)

    let create_grammar () =
        let expr_sk = Gram.Entry.mk "expr_sk" in
        EXTEND Gram 
        expr_sk: [
            [ UIDENT "S" -> S _loc ] | 
            [ UIDENT "K" -> K _loc ] | 
            [ v = UIDENT -> Label (v, _loc) ] | 
            [ v = LIDENT -> Label (v, _loc) ] |
            [ "("; x = LIST1 [ x = expr_sk -> x] ; ")" -> Node (x, _loc) ] ];
        END;
        expr_sk

(*            Gram.Entry.clear Syntax.str_item *)
    let init () =
        let expr_sk = create_grammar () in
        EXTEND Gram 
        expr_sk: [
            [ "["; x = Syntax.expr; "]" -> Expr (x, _loc) ] |
            [ "{"; x = Syntax.expr; "}" -> Debg (x, _loc) ] 
        ];
        END;
        EXTEND Gram
        Syntax.expr : LEVEL "top" [ 
            [ UIDENT "SK"; s = expr_sk -> bare_expression_converter s 
            | UIDENT "LS"; s = LIST0 [x = UIDENT -> x | x = LIDENT -> x] -> 
                    lst_to_expr _loc s 
            | UIDENT "SKLS"; s = LIST0 [ x = expr_sk -> x ] -> 
                    sklst_to_expr _loc s ]];
        END;;

    include Syntax
end

module IdSKOcaml = struct
    let name = "SKOcaml"
    let version = "0.1"
end

module SKOcamlLanguage (Gram : Grammar.Static) (Syntax : Camlp4Syntax) = struct
    open Syntax

    module MySKLanguage = SKLanguage (Gram) (Syntax)

    let () = MySKLanguage.init ()

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

    type 'a kolmo_function =
        | Module of (name * 'a kolmo_function list * Syntax.Ast.loc)
        | LetVal of (name * arguments * 'a definition * Syntax.Ast.loc)
        | RecVal of (name * arguments * 'a definition * Syntax.Ast.loc) 
        | Let of (name * arguments * 'a definition * Syntax.Ast.loc)
        | Rec of (name * arguments * 'a definition * Syntax.Ast.loc) 
    and name = (string * Syntax.Ast.loc)
    and arguments = name list 
    and 'a definition =
        | Value of ('a values * Syntax.Ast.loc)
        | IfElse of ('a condition * 'a definition * 'a definition *
        Syntax.Ast.loc)
        | Letin of ('a kolmo_function * 'a definition * Syntax.Ast.loc)
    and 'a condition = Condition of ('a definition * Syntax.Ast.loc)
    and 'a values =
        | Apply of (name list * name * ('a definition list) * Syntax.Ast.loc)
        | Integer of name 
        | Expr of 'a

    let rec to_ast = function
        | Module ((n, _), x, _loc) ->
                <:expr<`Module ($str:n$, 
                    $exSem_of_list (List.map to_ast x)$)>>
        | LetVal ((name, _), arguments, definition, _loc) ->
                <:expr<`LetVal($str:name$, 
                    $exSem_of_list (args_to_ast arguments)$, 
                    $def_to_ast definition$)>>
        | RecVal ((name, _), arguments, definition, _loc) ->
                <:expr<`RecVal($str:name$, 
                    $exSem_of_list (args_to_ast arguments)$, 
                    $def_to_ast definition$)>>
        | Let ((name, _), arguments, definition, _loc) ->
                <:expr<`Let($str:name$, 
                    $exSem_of_list (args_to_ast arguments)$, 
                    $def_to_ast definition$)>>
        | Rec ((name, _), arguments, definition, _loc) ->
                <:expr<`Rec($str:name$, 
                    $exSem_of_list (args_to_ast arguments)$, 
                    $def_to_ast definition$)>>
    and def_to_ast = function 
        | Value (values, _loc) -> <:expr<`Value $values_to_ast values$>>
        | IfElse (condition, def1, def2, _loc) ->
                <:expr<`IfElse ($condition_to_ast condition$,  
                    $def_to_ast def1$, $def_to_ast def2$)>>
        | Letin (f, d, _loc) ->
                <:expr<`Letin($to_ast f$, $def_to_ast d$)>>
    and string_to_ast (x, _loc) = <:expr<$str:x$>> 
    and args_to_ast args = List.map string_to_ast args 
    and condition_to_ast = function
        | Condition (def, _loc) -> <:expr<`Condition ($def_to_ast def$)>>
    and values_to_ast = function
        | Apply (a, b, defs, _loc) -> 
                <:expr<`Apply($exSem_of_list (args_to_ast a)$, 
                    $string_to_ast b$, 
                    $exSem_of_list (List.map def_to_ast defs)$)>>
        | Integer ((_, _loc) as x) -> <:expr<`Integer $string_to_ast x$>>
        | Expr (x, _loc) -> <:expr<`Expr $x$>>

    let rec to_ocaml = function
        | Module (name, x, _loc) ->
                `Module (fst name, List.map to_ocaml x) 
        | LetVal (name, arguments, definition, _loc) ->
                `LetVal (fst name, List.map fst arguments, def_to_ocaml definition)
        | RecVal (name, arguments, definition, _loc) ->
                `RecVal (fst name, List.map fst arguments, 
                    def_to_ocaml definition)
        | Let (name, arguments, definition, _loc) ->
                `Let (fst name, List.map fst arguments, def_to_ocaml definition)
        | Rec (name, arguments, definition, _loc) ->
                `Rec (fst name, List.map fst arguments, def_to_ocaml definition)
    and def_to_ocaml = function 
        | Value (values, _loc) -> 
                `Value (values_to_ocaml values)
        | IfElse (condition, def1, def2, _loc) ->
                `IfElse (condition_to_ocaml condition,  
                    def_to_ocaml def1, def_to_ocaml def2)
        | Letin (f, d, _loc) ->
                `Letin(to_ocaml f, def_to_ocaml d)
    and string_to_ocaml (x, _loc) = x
    and args_to_ocaml args = List.map string_to_ocaml args 
    and condition_to_ocaml = function
        | Condition (def, _loc) -> `Condition (def_to_ocaml def)
    and values_to_ocaml = function
        | Apply (a, b, defs, _loc) -> 
                `Apply((args_to_ocaml a), string_to_ocaml b, 
                    (List.map def_to_ocaml defs))
        | Integer ((_, _loc) as x) -> `Integer (string_to_ocaml x)
        | Expr (x, _) -> `Expr (MySKLanguage.sk_to_ocaml x)



    type loc =Syntax.Ast.loc

    type range = (string * string)

    type pair = (string * string * loc)
    type distribution = (string * (pair list) * loc)
    type options = 
        | EProbability of (pair list * loc)
        | FProbability of (distribution * loc)

    type 'a kolmo_specification =
        | Alphabet of (string * string list * options option * loc)
        | Character of (string * 'a * options option * loc)
        | WordSet of (string * string * range * options option * loc)
        | IntSet of (string * range * options option * loc)

    let pair_to_ocaml (a, b, _) = (a, float_of_string b)

    let options_to_ocaml = function
        | Some (EProbability (probs, _)) -> 
                Some (`EProbability (List.map pair_to_ocaml probs))
        | Some (FProbability ((name, parameters, _), _)) ->
                Some (`FProbability (name, List.map pair_to_ocaml parameters))
        | None -> None

    let range_to_ocaml (a, b) = (float_of_string a, float_of_string b)

    let rec spec_to_ocaml = function
        | Alphabet (name, elements, options, _) ->
                `Alphabet (name, elements, options_to_ocaml options)
        | Character (name, definition, options, _) ->
                `Character 
                    (name, List.map to_ocaml definition, 
                    options_to_ocaml options)
        | WordSet (name, alphabet, range, options, _) ->
                let range = range_to_ocaml range in
                `WordSet (name, alphabet, range, options_to_ocaml options)
        | IntSet (name, range, options, _) ->
                `IntSet (name, range_to_ocaml range, options_to_ocaml options)


    let create_specification machine_parser =
        let spec = Gram.Entry.mk "kolmo_spec" in 
        EXTEND Gram
        GLOBAL: spec;
        spec: [
            [ LIDENT "character"; name = UIDENT; ":"; specs = 
                LIST1 [x = machine_parser -> x]; opt = OPT options; 
                KEYWORD "end" -> Character (name, specs, opt, _loc) ] |
            [ LIDENT "alphabet"; name = UIDENT; ":"; alph = LIST1 [x = UIDENT ->
                x] ; opt = OPT options; KEYWORD "end" -> 
                    Alphabet (name, alph, opt, _loc) ] |
            [ LIDENT "wordset"; name = UIDENT; LIDENT "of"; alph = UIDENT; 
                "["; min = INT; max = INT; "]"; opt = OPT
                options; KEYWORD "end" -> 
                    WordSet (name, alph, (min, max), opt, _loc) ] |
            [ LIDENT "integers"; name = UIDENT; "["; min = INT; max = INT; 
            opt = OPT options; KEYWORD "end" -> IntSet (name, (min, max), opt,
            _loc) ] 
        ];
        options :[
            [ KEYWORD "with"; LIDENT "probability"; "{"; x = LIST1 [ x =
                assigned_prob -> x ]; "}" -> EProbability (x, _loc) ] |
            [ KEYWORD "with"; d = distribution -> FProbability (d, _loc) ] 
        ];
        assigned_prob :[
            [ x = LIDENT; "="; y = FLOAT -> (x, y, _loc) ]
        ];
        distribution :[
            [ x = LIDENT; parameters =  LIST1 [ x = assigned_prob -> x ] -> 
                (x, parameters, _loc) ]
        ];
        END;
        spec

    let create_grammar my_expression =
        let ocaml_sk = Gram.Entry.mk "ocaml_sk" in
        let definition = Gram.Entry.mk "definition" in
        EXTEND Gram 
        GLOBAL: ocaml_sk definition;
        ocaml_sk: [
            [ KEYWORD "module"; n = UIDENT; "="; KEYWORD "struct";
                x = LIST1 [ x = ocaml_sk -> x]; KEYWORD "end" ->
                    Module ((n, _loc), x, _loc) ] |
            [ KEYWORD "val";  KEYWORD "rec";
                (name, arguments, definition) = sk_definition ->
                    RecVal (name, arguments, definition, _loc) ] |
            [ KEYWORD "val"; 
                (name, arguments, definition) = sk_definition ->
                    LetVal (name, arguments, definition, _loc) ] |
            [ KEYWORD "let";  KEYWORD "rec";
                (name, arguments, definition) = sk_definition ->
                    Rec (name, arguments, definition, _loc) ] |
            [ KEYWORD "let"; 
                (name, arguments, definition) = sk_definition ->
                    Let (name, arguments, definition, _loc) ] 
        ];
        sk_definition: [
            [ arguments = LIST1 [ x = LIDENT -> (x, _loc)];
                "="; 
                def = definition -> 
                    match arguments with
                    | name :: arguments ->
                            (name, arguments, def) 
                    | [] -> assert false]
        ];
        definition: [
            [ "("; x = definition; ")" -> x ] |
            [ KEYWORD "if"; x = condition; KEYWORD "then"; y = definition; 
                KEYWORD "else"; z = definition -> 
                    IfElse (x, y, z, _loc) ] |
            [ f = ocaml_sk; KEYWORD "in"; d = definition -> 
                    Letin (f, d, _loc) ] |
            [ v = values -> Value (v, _loc) ] 
        ];
        condition: [
            [ a = definition -> Condition (a, _loc) ]
        ];
        values: [
            [ a = LIDENT; b = LIST0 
                [ x = sub_definition -> x] -> Apply ([], (a, _loc), b, _loc) ] |
            [ mods = LIST1 [x = UIDENT; "." -> (x, _loc)];
                a = LIDENT; b = LIST0 
                [ x = sub_definition -> x] -> Apply (mods, (a, _loc), b, _loc) ] |
            [ b = INT -> Integer (b, _loc) ] 
                | [ "["; e = my_expression; "]" -> Expr (e, _loc) ] 
        ];
        sub_definition: [
            [ "("; x = definition; ")" -> x ] |
            [ x = LIDENT -> Value (Apply ([], (x, _loc), [], _loc), _loc) ] |
            [ b = INT -> Value (Integer (b, _loc), _loc) ] |
            [ mods = LIST1 [x = UIDENT; "." -> (x, _loc) ];
                a = LIDENT; b = LIST0 
                [ x = sub_definition -> x] -> 
                Value (Apply (mods, (a, _loc), b, _loc), _loc) ] 
            | [ "["; e = my_expression; "]" -> Value (Expr (e, _loc), _loc) ] 
        ];
        END;
        ocaml_sk, definition

    let second_extend () =
        let ocaml_sk, definition = create_grammar expr in
        EXTEND Gram
        Syntax.expr : LEVEL "top" [ 
            [ UIDENT "OCAMLSK"; s = LIST1 [ x = ocaml_sk -> x] -> 
                <:expr<$exSem_of_list (List.map to_ast s)$>> ] |
            [ UIDENT "OCAMLSKV"; s = definition -> <:expr<$def_to_ast s$>> ] 
        ];
        END;;

    include Syntax
end

module SKOcamlLanguageExt (Syntax : Camlp4Syntax) = struct
    module S = SKOcamlLanguage (Camlp4.PreCast.Gram) (Syntax)
    open S
    let () = 
        let () = second_extend () in
        ()
    include S
end

module Gram = Camlp4.PreCast.Gram 
module MySKOcaml = SKOcamlLanguage (Gram) (Camlp4.PreCast.Syntax)

let of_stream str =
    try
        let sk_expr = MySKOcaml.MySKLanguage.create_grammar () in
        let ocaml_sk, _ = MySKOcaml.create_grammar sk_expr in
        let kolmo_spec = MySKOcaml.create_specification ocaml_sk in
        MySKOcaml.spec_to_ocaml (Gram.parse kolmo_spec (Camlp4.PreCast.Loc.mk "<stream>")
        str)
    with
    | (Camlp4.PreCast.Loc.Exc_located (err, _)) as error ->
            print_endline (Camlp4.PreCast.Loc.to_string err);
            raise error

let of_channel ch =
    of_stream (Stream.of_channel ch)

let of_string str =
    of_stream (Stream.of_string str)
