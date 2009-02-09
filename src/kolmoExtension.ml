(* A File that defines OCaml extensions and a DSL for SK operations and
* definitions *)


open Camlp4.Sig

type primitives = [ `S | `K | `Label of string | `Node of primitives list | `Debugger of string | `Lazy of primitives Lazy.t]

type 'a kolmo_function =
    [ `Module of 'a kolmo_function list
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

module SKLanguage (Syntax : Camlp4Syntax) = struct
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

    let expression_converter _loc x =
        let res = bare_expression_converter x in
        <:expr< `Processed $res$>>

    let lst_to_expr _loc lst = 
        exSem_of_list (List.map (fun x -> <:expr<$str:x$>>) lst)

    let expr_sk = Gram.Entry.mk "expr_sk"

    let sklst_to_expr _loc lst =
        let rec aux_to_expr lst = 
            match lst with
            | [] -> <:expr<`Label "EmptyList">>
            | h :: t -> 
                    let at = aux_to_expr t in
                    <:expr<`Node [`Label "Pair"; $h$; $at$]>>
        in
        aux_to_expr (List.map bare_expression_converter lst)
    EXTEND Gram 
    expr_sk: [
        [ UIDENT "S" -> S _loc ] | 
        [ UIDENT "K" -> K _loc ] | 
        [ v = UIDENT -> Label (v, _loc) ] | 
        [ v = LIDENT -> Label (v, _loc) ] |
        [ "["; x = expr; "]" -> Expr (x, _loc) ] |
        [ "{"; x = expr; "}" -> Debg (x, _loc) ] |
        [ "("; x = LIST1 [ x = expr_sk -> x] ; ")" -> Node (x, _loc) ] ];
    END;;

(*            Gram.Entry.clear Syntax.str_item *)

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

let () =
    let module M = Camlp4.Register.OCamlSyntaxExtension (Id) (SKLanguage) in 
    ()

module IdSKOcaml = struct
    let name = "SKOcaml"
    let version = "0.1"
end

module SKOcamlLanguage (Syntax : Camlp4Syntax) = struct
    open Syntax

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

    let ocaml_sk = Gram.Entry.mk "ocaml_sk"
    let definition = Gram.Entry.mk "definition"

    EXTEND Gram 
    GLOBAL: ocaml_sk definition;
    ocaml_sk: [
        [ KEYWORD "module"; n = UIDENT; "="; KEYWORD "struct";
            x = LIST1 [ x = ocaml_sk -> x]; KEYWORD "end" ->
                <:expr<`Module ($str:n$, $exSem_of_list x$)>> ] |
        [ KEYWORD "let";  KEYWORD "rec";
            (name, arguments, definition) = sk_definition ->
                <:expr<`Rec ($name$, $arguments$, $definition$)>> ] |
        [ KEYWORD "let"; 
            (name, arguments, definition) = sk_definition ->
                <:expr<`Let ($name$, $arguments$, $definition$)>> ] 
    ];
    sk_definition: [
        [ arguments = LIST1 [ x = LIDENT -> <:expr<$str:x$>>];
            "="; 
            def = definition -> 
                match arguments with
                | name :: arguments ->
                        (<:expr<$name$>>, exSem_of_list arguments,
                        <:expr<$def$>>) 
                | [] -> assert false]
    ];
    definition: [
        [ "("; x = definition; ")" -> x ] |
        [ KEYWORD "if"; x = condition; KEYWORD "then"; y = definition; 
            KEYWORD "else"; z = definition -> 
                <:expr<`IfElse ($x$, $y$, $z$)>> ] |
        [ f = ocaml_sk; KEYWORD "in"; d = definition -> 
                <:expr<`Letin ($f$, $d$)>> ] |
        [ v = values -> <:expr<`Value $v$>> ] 
    ];
    condition: [
        [ a = definition -> <:expr<`Condition $a$>> ]
    ];
    values: [
        [ a = LIDENT; b = LIST0 
            [ x = sub_definition -> x] -> 
            <:expr<`Apply ([], $str:a$, $exSem_of_list b$)>> ] |
        [ mods = LIST1 [x = UIDENT; "." -> <:expr<$str:x$>>];
            a = LIDENT; b = LIST0 
            [ x = sub_definition -> x] -> 
            <:expr<`Apply ($exSem_of_list mods$, $str:a$, $exSem_of_list b$)>> ] |
        [ b = INT -> <:expr<`Integer $str:b$>> ] |
        [ "["; e = expr; "]" -> <:expr<`Expr $e$>> ] 
    ];
    sub_definition: [
        [ "("; x = definition; ")" -> x ] |
        [ x = LIDENT -> <:expr<`Value (`Apply ([], $str:x$, []))>> ] |
        [ b = INT -> <:expr<`Value (`Integer $str:b$)>> ] |
        [ "["; e = expr; "]" -> <:expr<`Value (`Expr $e$)>> ] 
    ];
    END;;

    EXTEND Gram
    Syntax.expr : LEVEL "top" [ 
        [ UIDENT "OCAMLSK"; s = LIST1 [ x = ocaml_sk -> x] -> 
            <:expr<List.iter Kolmo.Compiler.compile $exSem_of_list s$>> ] |
        [ UIDENT "OCAMLSKV"; s = definition -> 
            <:expr<Kolmo.S_K.eval (Kolmo.Compiler.evaluate $s$)>> ] 
    ];
    END;;

    include Syntax
end

let () =
    let module M = Camlp4.Register.OCamlSyntaxExtension (IdSKOcaml) (SKOcamlLanguage) in 
    ()
