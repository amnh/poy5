(* A File that defines OCaml extensions and a DSL for SK operations and
* definitions *)


open Camlp4.Sig

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
