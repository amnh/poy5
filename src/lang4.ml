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

let () = SadmanOutput.register "Lang4" "$Revision: 1644 $"

exception Inconsistent_types

type var = [ `Var of string ]       (* A Variable *)
type constant =
    [ `CInt of int 
    | `CBase of string 
    | `CSeq of string ]

(* The only comparison allowed in the u language is a variable vs. a constant *)
type boolean = 
    [`Boolean of var * constant ]   

(* A function application, could be either a simple variable or a complete
* function with it's parameters. The variant type will be fixed in the
* expression definition. This applies to all the following variant types. *)
type 'a func = 
    [ var 
    | `Function of string * 'a list ]
(* The only three valid types in the u language, integer, base and sequence. *)
type 'a integer = 
    [ `Successor of 'a integer 
    | `Predeccessor of 'a integer 
    | `CInt of int 
    | 'a func ]
type 'a base = 
    [ `Head of 'a sequence 
    | `CBase of string 
    | 'a func ]
and 'a sequence = 
    [ `Prepend of 'a base * 'a sequence 
    |`Tail of 'a sequence 
    | `CSeq of string 
    | 'a func ]

(* Finally, a valid expression in the u language for function specification *)
type expression = 
    [ expression func 
    | expression integer 
    | expression base
    | expression sequence 
    | `If_then of boolean * expression * expression ]

(* A function in the u language, specified with the triple (a, x, y), where [a]
 * is the name of the function, [x] is the list of names of the parameters and
 * [y] is the expression that specifies the function. *)
type a_ufunc =  [ `Aufunc of var * var list * expression ]
(* A character definition requires a tuple (a, b) where [a] is the name of the
* character and [b] is a list of functions that specify its valid transformation
* events *)
type a_poy_macro = [`Acharacter of string * a_ufunc list ] list

let addpar seq = print_string (" " ^ seq ^ " ") 
let inpar f el txt = 
    print_string ("( " ^ txt ^ " ");
    List.iter f el;
    print_string " )"

let rec print_expression (x : [< expression ]) : unit = 
    let print_function : [< expression func ] -> unit = function
            | `Var name -> print_string (name ^ " ");
            | `Function (name, param) -> 
                    print_string ("( " ^ name ^ " ");
                    List.iter (fun x -> print_expression x; print_string " ";) param;
    in
    let rec print_integer : [< expression integer ] -> unit = function
            | `Successor i -> inpar (print_integer) [i] "succ";
            | `Predeccessor i -> inpar (print_integer) [i] "pred";
            | `CInt i -> print_int i;
            | #func as f -> print_function f
    in
    let rec print_base x =
        let rec print_sequence = function
            | `Tail seq -> inpar (print_sequence) [seq] "tail";
            | `Prepend (b, seq) -> 
                    print_string "prepend ( ";
                    print_base b;
                    print_string ", ";
                    print_sequence seq;
                    print_string ")";
            | `CSeq seq -> addpar seq;
            | #func as f -> print_function f;
        in
        match x with
            | `Head seq -> inpar (print_sequence) [seq] "head";
            | `CBase base -> addpar base; 
            | #func as f -> print_function f;
    in
    let rec print_sequence = function
            | `Tail seq -> inpar (print_sequence) [seq] "tail";
            | `Prepend (b, seq) -> 
                    print_string "prepend ( ";
                    print_base b;
                    print_string ", ";
                    print_sequence seq;
                    print_string ")";
            | `CSeq seq -> addpar seq;
            | #func as f -> print_function f;
    in
    let print_constant = function
        | `CInt i -> print_int i;
        | `CSeq s | `CBase s -> print_string s
    in
    match x with
    | #func as f -> print_function f
    | #integer as i -> print_integer i
    | #base as b -> print_base b
    | #sequence as s -> print_sequence s
    | `If_then (`Boolean (`Var a, b), c, d) ->
            print_string ("( if " ^ a ^ " == ");
            print_constant b;
            print_string " then ";
            print_expression c;
            print_string " else ";
            print_expression d;
            print_string ")"

let rec print = function 
    | #a_ufunc as x ->
            begin match x with
            | `Aufunc (`Var name, parameters, definition) ->
                    print_string ("let " ^ name ^ " ");
                    List.iter 
                        (fun (`Var x) -> print_string (x ^ " ")) parameters;
                    print_expression definition;
                    print_string ("\n end of function \n");
            | `Aprev_ufunc name -> print_string ("Character " ^ name ^ " is \
            already defined.\n")
            end

let print_macro = function
    | `Acharacter (name, its) -> List.iter print its

type ht = 
    [ `Var of string 
    | `Function of string 
    | `Successor 
    | `Predeccessor
    | `CInt
    | `CSeq
    | `CBase
    | `Tail
    | `Head
    | `Prepend 
    | `If_then ]

let _counter = ref ((Hashtbl.create 10) : (ht, int) Hashtbl.t)

let clear () = _counter := Hashtbl.create 10

let do_count x =
    try
        let res = Hashtbl.find!_counter x in
        Hashtbl.replace !_counter x (res + 1);
    with
    | Not_found -> Hashtbl.add !_counter x 1

let output_counter ch counter =
    let do_print a b =
        let printer_match = function
            | `Var name -> name 
            | `Function name -> name 
            | `Successor -> "successor"
            | `Predeccessor -> "predeccessor"
            | `CInt -> "integer"
            | `Tail -> "tail"
            | `Prepend -> "prepend"
            | `CSeq -> "sequence"
            | `Head -> "head"
            | `CBase -> "base"
            | `If_then -> "If_then"
        in
        let res = printer_match a 
        and co = string_of_int b in
        output_string ch (res ^ ": " ^ co ^ "\n");
    in
    Hashtbl.iter do_print counter 

let count_if_then_else x = ()

let rec count_expression (x : expression ) : unit = 
    let it, (lst : [< expression] list) = 
    match x with
    | `Var name -> `Var name, []
    | `Function (name, param) -> (`Function name), param
    | `Successor i -> `Successor, [(i :> expression) ]
    | `Predeccessor i -> `Predeccessor, [(i :> expression)]
    | `CInt i -> `CInt, ([] : expression list)
    | `Tail seq -> `Tail, [(seq :> expression)]
    | `Prepend (b, seq) ->  `Prepend, [(seq :> expression); (b :> expression)]
    | `CSeq seq -> `CSeq, []
    | `Head seq -> `Head, [(seq :> expression)]
    | `CBase base -> `CBase, []
    | `If_then (`Boolean (`Var a, b), c, d) -> 
            count_if_then_else b;
            `If_then, [c; d]
    in
    do_count it;
    List.iter count_expression lst

let count_character_func (x : a_ufunc) = 
    match x with
    | `Aufunc (name, parameters, exp) -> count_expression exp

let count_character (x : [`Acharacter of string * a_ufunc list ]) = 
    match x with
    | `Acharacter (name, defs) ->
            List.iter count_character_func defs;
            let res = name, !_counter in
            clear ();
            res

let count_poy_macro (x : a_poy_macro) = 
    List.map count_character x

type to_store = [`Var of string | `Function of string]
type types = [`Unknown | `CBase | `CSeq | `CInt ]
let _types : (to_store, types) Hashtbl.t  ref = ref (Hashtbl.create 10)

let clear_types () = _types := Hashtbl.create 10

let set_type t v =
    try
        match Hashtbl.find !_types v, t with 
        | a, `Unknown -> ()
        | `Unknown, _ -> Hashtbl.replace !_types v t
        | a, t when a = t -> ()
        | _, _ -> raise Inconsistent_types
    with
    | Not_found -> Hashtbl.add !_types v t

let get_type v =
    try Hashtbl.find !_types v with Not_found -> `Unknown

let _par_types = ref (Hashtbl.create 10)

let clear_par_types () = _par_types := Hashtbl.create 10

let get_par_types v par = 
    try Hashtbl.find !_par_types v 
    with Not_found -> List.map (fun x -> `Unknown) par

let check_parameter_types nl pl var =
    let rec checker a b =
        match a, b with
        | ha :: ta, hb :: tb when ha = hb -> ha :: checker ta tb
        | ha :: ta, `Unknown :: tb -> ha :: checker ta tb
        | [], [] -> []
        | _, _ -> raise Inconsistent_types
    in
    let res = checker nl pl in
    Hashtbl.replace !_par_types var res

let rec update_types t (x : expression ) = 
    match x with
    | `Var name as x-> 
            set_type t (`Var name);
            get_type x 
    | `Function (name, param) -> 
            set_type t (`Function name);
            let parameter_types = get_par_types name param in
            let to_apply = List.map update_types parameter_types in
            let rec apply_pairs a b =
                match a, b with
                | ha :: ta, hb :: tb -> 
                        (ha hb) :: apply_pairs ta tb
                | [], [] -> []
                | _, _ -> raise Inconsistent_types
            in
            let res = apply_pairs to_apply param in
            check_parameter_types res parameter_types name;
            t
    | `Successor i -> 
            begin match t with
            | `Unknown | `CInt -> update_types (`CInt) (i :> expression)
            | _ -> raise Inconsistent_types
            end 
    | `Predeccessor i ->
            begin match t with
            | `Unknown | `CInt -> update_types (`CInt) (i :> expression)
            | _ -> raise Inconsistent_types
            end
    | `CInt i -> `CInt
    | `Tail seq ->
            begin match t with
            | `Unknown | `CSeq -> update_types (`CSeq) (seq :> expression)
            | _ -> raise Inconsistent_types
            end
    | `Prepend (b, seq) ->
            begin match t with
            | `Unknown | `CSeq -> 
                    begin match update_types (`CBase) (b :> expression) with
                    | `CBase -> update_types (`CSeq) (seq :> expression)
                    | _ -> raise Inconsistent_types
                    end;
            | _ -> raise Inconsistent_types
            end
    | `CSeq seq -> `CSeq
    | `Head seq -> 
            begin match t with
            | `Unknown | `CBase -> 
                    let _ = update_types (`CSeq) (seq :> expression) in
                    `CBase
            | _ -> raise Inconsistent_types
            end 
    | `CBase base -> `CBase
    | `If_then (`Boolean (`Var a, b), c, d) -> 
            let t = update_types (`Unknown) (c :> expression) in
            update_types t (d :> expression)

let check_char_func (x : a_ufunc) =
    match x with
    | `Aufunc (a, b, c) ->
            clear_par_types ();
            clear_types ();
            let res1 = update_types `Unknown c in
            update_types
            res1 c 

let check_character (x : [`Acharacter of string * a_ufunc list]) = 
    match x with
    | `Acharacter (_, lst) -> List.map check_char_func lst

let check_poy_macro (x : a_poy_macro) = List.map check_character x
