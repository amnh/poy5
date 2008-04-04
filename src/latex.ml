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

type parsed = 
    | Word of string
    | WordNoSpace of string
    | Blank
    | Text of parsed list
    | Command of (string * parsed list)

let channel = ref stdout

let o str = output_string !channel str

let rec produce_element = function
    | Command ("begin", (Word h :: tl)) ->
            (match h with
            | "command" -> 
                    (match tl with
                    | h :: _ -> 
                            o "@]\n\n";
                            produce_element h;
                            o "\n.\n@[<v 2>"
                    | [] -> failwith "command with no args?")
            | "description" ->
                    o "@,@[<v 2>@,@["
            | "poydescription" ->
                    o "@,@[<v 2>@{<c:cyan>Description@}@,@,@["
            | "arguments" ->
                    o "@,@,@[<v 2>@{<c:cyan>Arguments@}@,@,@[<v>"
            | "argumentgroup" ->
                    (match tl with
                    | [title ; description] ->
                            o "@,@[<v 2>@{<c:cyan>@[";
                            produce_element title;
                            o "@]@}@,@[";
                            produce_element description;
                    | [title] ->
                            o "@,@[<v 2>@{<c:cyan>@[";
                            produce_element title;
                            o "@]@}@,@[";
                    | _ -> failwith "argumentgroup without the necessary args?")
            | "statement" -> o "@,@[@,@["
            | "poyexamples" -> o "@,@[<v 2>@{<c:cyan>Examples@}@,@,@[<v>"
            | "poyalso" -> o "@,@[<v 2>@[@{<c:cyan>See Also@}@]@,@,@["
            | "flushleft" 
            | "center" -> o "@[<v 2>@,@[" 
            | "atsymbol" -> o "@@"
            | _ -> ())
    | Command ("%", []) -> o "%"
    | Command (h, []) when h = "" -> ()
    | Command (h, []) when h.[0] = '_' -> o h
    | Command ("end", [Word h]) -> o "@]@]@,"
    | Command ("argumentdefinition", [com ; args ; definition ; cross_reference]) ->
            o "@[<v 2>@,@[@{<c:yellow>";
            produce_element com;
            produce_element args;
            o "@}@]@,@[";
            produce_element definition;
            o "@]@]@,"
    | Command ("poydefaults", [args; descr]) ->
            o "@,@,@[<v 2>@{<c:cyan>Defaults@}@,@,@[";
            produce_element args;
            o "@]@,@[";
            produce_element descr;
            o "@]@]@,"
    | Command ("obligatory", [arg]) -> 
            o ": ";
            produce_element arg
    | Command ("optional", [arg]) -> 
            o "[: ";
            produce_element arg;
            o "]";
    | Command ("poyexample", [example; explanation]) ->
            o "@[<v 2>@[@{<c:green>";
            produce_element example;
            o "@}@]@,@[";
            produce_element explanation;
            o "@]@]@,@,";
    | Command ("ncross", [arg; _])
    | Command ("cross", [arg]) ->
            o "@[";
            produce_element arg;
            o "@]@,";
    | Command ("poycommand", arg) ->
            o "@[<h>";
            List.iter produce_element arg;
            o "@]";
    | Command ("ccross", [Word arg])
    | Command ("nccross", [Word arg; _]) ->
            o (arg ^ "(see help (" ^ arg ^ ")) ");
    | Command _ -> ()
    | Text lst -> 
            List.iter produce_element lst;
    | Word x -> 
            if x <> "~" then o x;
            (*
            if x <> "(" && x <> "[" then o " ";
            *)
    | WordNoSpace x -> 
            if x <> "~" then o x
    | Blank -> o " "

let rec collapse = function
    | Blank
    | (WordNoSpace _)
    | (Word _) as x -> x
    | Text [(Text _) as y] -> collapse y
    | Text [(Word _) as y] -> y
    | Text x ->
            Text (List.map collapse x)
    | Command ("poy", _) -> Word "POY"
    | Command ("poybool", _) -> Word "BOOL"
    | Command ("poyfloat", _) -> Word "FLOAT"
    | Command ("poyint", _) -> Word "INTEGER"
    | Command ("poystring", _) -> Word "STRING"
    | Command ("poylident", _) -> Word "LIDENT"
    | Command ("poycommand", [arg]) ->
            Command ("poycommand", [collapse arg])
    | Command ("poyargument", [arg])
    | Command ("texttt", [arg])
    | Command ("emph", [arg]) ->
            collapse arg
    | Command (a, b) ->
            Command (a, List.map collapse b)

let rec flatten x =
    match x with
    | ((WordNoSpace _) as h) :: t 
    | ((Word _) as h) :: t -> h :: (flatten t)
    | (Text x) :: t -> (flatten x) @ (flatten t)
    | h :: t -> h :: (flatten t)
    | [] -> []

let rec collapse2 = function
    | (Text a) :: t ->
            collapse2 ((flatten a) @ t)
    | (Word a) :: (Command (h, [])) :: t when h = "" ->
            Word h :: (collapse2 t)
    | (Word a) :: (Command (h, [])) :: t when h.[0] = '_' ->
            Word (a ^ h) :: (collapse2 t)
    | ((Word a) as h) :: (((Word b) :: _) as t) ->
            if b = "." || b = "," || b = ")" || b = "]" || b = "," 
                || b = ";" then
                (WordNoSpace a) :: (collapse2 t)
            else h :: (collapse2 t)
    | h :: t -> h :: (collapse2 t)
    | [] -> []

let rec the_parser fstream = 
    let is_newline fstream =
        fstream#skip_ws_nl;
        fstream#match_prefix "\\\\"
    in
    let is_command fstream = 
        fstream#skip_ws_nl;
        fstream#match_prefix "\\"
    in
    let is_word fstream = 
        fstream#skip_ws_nl;
        not (fstream#match_prefix "\\")
    in
    let is_enclosed fstream = 
        fstream#skip_ws_nl;
        fstream#match_prefix "{"
    in
    let is_close_command fstream = 
        fstream#skip_ws_nl;
        fstream#match_prefix "}"
    in
    let is_end_of_file fstream = 
        match fstream#getch_safe with
        | Some x -> fstream#putback x; false
        | None -> true
    in
    let is_blank fstream =
        let res = fstream#read_incl [' '; '\t'; '\010'; '\013'] in
        if res = "" then false
        else true
    in
    let rec get_param acc fstream =
        fstream#skip_ws_nl;
        if fstream#match_prefix "{" then begin
            let res = the_parser fstream in
            get_param (acc @ [(Text res)]) fstream
        end else acc
    in
    let skip1 = [' '; '\t'; '\\'; '\010'; '\013'; '{'; '}'; ','; '.'; ')'; '(';
            '['; ']'; '~'] 
    and noskip = ['('; ')'; ','; '.'; '['; ']'; '~'] in
    let get_comm fstream =
        fstream#skip_ws_nl;
        let cmd = fstream#read_excl skip1 in
        let param = get_param [] fstream in
        Command (cmd, param)
    in
    let get_word fstream = 
        fstream#skip_ws_nl;
        let cmd = fstream#read_excl skip1 in
        let cmd = 
            if cmd = "" then fstream#read_incl noskip
            else cmd
        in
        Text [Word cmd]
    in
    let rec split_on_commands acc fstream = 
        if is_end_of_file fstream then
            List.rev acc
        else if is_blank fstream then
            split_on_commands (Blank :: acc) fstream
        else if is_close_command fstream then 
            List.rev acc
        else if is_newline fstream then
            split_on_commands ((Text [Word "@\n"]) :: acc) fstream
        else if is_command fstream then
            split_on_commands ((get_comm fstream) :: acc) fstream
        else if is_enclosed fstream then
            split_on_commands ((Text (the_parser fstream)) :: acc) fstream
        else if is_word fstream then
            split_on_commands ((get_word fstream) :: acc) fstream 
        else failwith "Huh?";
    in
    let res = split_on_commands [] fstream in
    collapse2 (List.map collapse res)

let process fstree = 
    let res = the_parser fstree in
    List.iter produce_element res

let process_file filename = 
    let ch = FileStream.Pervasives.open_in (`Local filename) in 
    process ch;
    o "\n\n\n"

let _ = process_file "../doc/allcommands.tex";;
