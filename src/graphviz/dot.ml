(**************************************************************************)
(*                                                                        *)
(*  Ocamlgraph: a generic graph library for OCaml                         *)
(*  Copyright (C) 2004-2010                                               *)
(*  Sylvain Conchon, Jean-Christophe Filliatre and Julien Signoles        *)
(*                                                                        *)
(*  This software is free software; you can redistribute it and/or        *)
(*  modify it under the terms of the GNU Library General Public           *)
(*  License version 2.1, with the special exception on linking            *)
(*  described in file LICENSE.                                            *)
(*                                                                        *)
(*  This software is distributed in the hope that it will be useful,      *)
(*  but WITHOUT ANY WARRANTY; without even the implied warranty of        *)
(*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.                  *)
(*                                                                        *)
(**************************************************************************)

(** Parser for DOT file format *)

open Dot_ast

let of_channel c =
    let lb = Lexing.from_channel c in
    let dot =
        try Dot_parser.file Dot_lexer.token lb
        with Parsing.Parse_error ->
            let n = Lexing.lexeme_start lb in
            failwith (Printf.sprintf "Dot.parse: parse error character %d" n)
    in
    dot

let of_file f =
    let chan = open_in f in
    let res = of_channel chan in
    close_in chan;
    res
