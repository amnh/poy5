(* POY 5.0 Alpha. A phylogenetic analysis program using Dynamic Homologies.   *)
(* Copyright (C) 2011 Andrés Varón, Lin Hong, Nicholas Lucaroni, Ward Wheeler,*)
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

let debug = false

type parsed = 
    | Word of string
    | WordNoSpace of string
    | Blank
    | Text of parsed list
    | Command of (string * parsed list)

let channel = ref stdout

let print_parsed l =
    let rec print_parsed d x = 
        let () = print_depth d in match x with
        | Word x        -> Printf.printf "Word:%s\n" x
        | WordNoSpace x -> Printf.printf "Word:%s\n" x
        | Blank         -> Printf.printf "' '\n";
        | Text l        -> List.iter (print_parsed (d+1)) l;
        | Command (x,l) -> Printf.printf "CMD:%s\n" x;
                           List.iter (print_parsed (d+1)) l
    and print_depth x =
        if x = 0 
            then ()
            else begin
                Printf.printf "\t";
                print_depth (x-1)
            end
    in
    print_parsed 0 l

let rec produce_latex channel data =
    let o str = output_string channel str in
    let rec produce_latex = function
        | (Command ("begin", (Word h :: tl))) as command ->
            begin match h with
                | "command" ->
                    if debug then print_parsed command;
                    begin match tl with
                        | h::_ -> o "@]\n\n.\n"; produce_latex h; o "\n@[<v 2>"
                        | []   -> failwith "command with no args?"
                    end
                | "description" -> o "@,@[<v 2>@,@["
                | "poydescription" ->
                    o "@,@[<v 2>@{<c:cyan>Description@}@,@,@["
                | "arguments" ->
                    o "@,@,@[<v 2>@{<c:cyan>Arguments@}@,@,@[<v>"
                | "argumentgroup" ->
                    begin match tl with
                        | [title ; description] ->
                            o "@,@[<v 2>@{<c:cyan>@[";
                            produce_latex title;
                            o "@]@}@,@[";
                            produce_latex description
                        | [title] ->
                            o "@,@[<v 2>@{<c:cyan>@[";
                            produce_latex title;
                            o "@]@}@,@["
                        | _ -> 
                            failwith "argumentgroup without the necessary args?"
                    end
                | "statement" -> o "@,@[@,@["
                | "poyexamples" -> o "@,@[<v 2>@{<c:cyan>Examples@}@,@,@[<v>"
                | "poyalso" -> o "@,@[<v 2>@[@{<c:cyan>See Also@}@]@,@,@[<v>"
                | "flushleft"
                | "center" -> o "@[<v 2>@,@[" 
                | "atsymbol" -> o "@@"
                | _ -> ()
            end
        | Command ("%", []) -> o "%"
        | Command (h, []) when h = "" -> ()
        | Command (h, []) when h.[0] = '_' -> o h
        | Command ("end", [Word h]) -> o "@]@]@,"
        | Command ("argumentdefinition", [com ; args ; definition ; cross_reference]) ->
            o "@[<v 2>@,@[@{<c:yellow>";
            produce_latex com;
            produce_latex args;
            o "@}@]@,@[";
            produce_latex definition;
            o "@]@]@,"
        | Command ("poydefaults", [args; descr]) ->
            o "@,@,@[<v 2>@{<c:cyan>Defaults@}@,@,@[";
            produce_latex args;
            o "@]@,@[";
            produce_latex descr;
            o "@]@]@,"
        | Command ("obligatory", [arg]) -> 
            o ": ";
            produce_latex arg
        | Command ("optional", [arg]) -> 
            o "[: ";
            produce_latex arg;
            o "]";
        | Command ("poyexample", [example; explanation]) ->
            o "@[<v 2>@[@{<c:green>";
            produce_latex example;
            o "@}@]@,@[";
            produce_latex explanation;
            o " @]@]@,@,@?"
        | Command ("ncross", [arg; _])
        | Command ("cross", [arg]) ->
            o "@["; produce_latex arg; o "@]@,"
        | Command ("poycommand", arg) ->
            o "@[<h>"; List.iter produce_latex arg; o " @]"
        | Command ("ccross", [Word arg])
        | Command ("nccross", [Word arg; _]) ->
            o (arg ^ " (see help (" ^ arg ^ ")) ");
        | Command _ -> ()
        | Text lst -> 
            List.iter produce_latex lst;
        | Word x -> 
            if x <> "~" then o x;
        | WordNoSpace x -> 
            if x <> "~" then o x
        | Blank ->
            o " "
    in
    produce_latex data

let rec produce_troff channel data =
    let o str = output_string channel str in
    let rec produce_troff = function
        | Command ("begin", (Word h :: tl)) ->
            begin match h with
                | "command" -> 
                    begin match tl with
                        | h::_ -> o "\n.SH "; produce_troff h; o "\n.P\n"
                        | []   -> failwith "command with no args?"
                    end
                | "description" ->
                        o "\n.P\n"
                | "poydescription" ->
                        o "\n.SS Description\n"
                | "arguments" ->
                        o "\n.SS Arguments\n"
                | "argumentgroup" ->
                    begin match tl with
                        | [title ; description] ->
                            o "\n.SS \"";
                            produce_troff title;
                            o "\"\n.P\n";
                            produce_troff description
                        | [title] ->
                            o "\n.SS \"";
                            produce_troff title;
                            o "\"\n.P\n"
                        | _ ->
                            failwith "argumentgroup without the necessary args?"
                    end
                | "statement" -> o "\n.P\n"
                | "poyexamples" -> 
                        o "\n.SS Examples\n"
                | "poyalso" -> 
                        o "\n.SS \"See Also\"\n"
                | "flushleft" 
                | "center" -> o "\n.P\n"
                | "atsymbol" -> o "@"
                | _ -> ()
            end
        | Command ("%", []) -> o "%"
        | Command (h, []) when h = "" -> ()
        | Command (h, []) when h.[0] = '_' -> o h
        | Command ("end", [Word h]) -> o "\n.P\n"
        | Command ("argumentdefinition", [com ; args ; definition ; cross_reference]) ->
            o "\n.IP \"";
            produce_troff com;
            produce_troff args;
            o "\"\n";
            produce_troff definition;
            o "\n"
        | Command ("poydefaults", [args; descr]) ->
            o "\n.SS Defaults\n.P\n";
            produce_troff args;
            o "\n.P\n";
            produce_troff descr;
            o "\n";
        | Command ("obligatory", [arg]) -> 
            o ": ";
            produce_troff arg
        | Command ("optional", [arg]) -> 
            o "[: ";
            produce_troff arg;
            o "]";
        | Command ("poyexample", [example; explanation]) ->
            o "\n.P\n";
            produce_troff example;
            o "\n.P\n";
            produce_troff explanation;
            o "\n";
        | Command ("ncross", [arg; _])
        | Command ("cross", [arg]) ->
            produce_troff arg;
            o "\n.br\n";
        | Command ("poycommand", arg) ->
            o " ";
            List.iter produce_troff arg;
            o " ";
        | Command ("ccross", [Word arg])
        | Command ("nccross", [Word arg; _]) ->
            o (arg ^ "(see help (" ^ arg ^ ")) ");
        | Command _ -> ()
        | Text lst -> 
            List.iter produce_troff lst;
        | Word x -> 
            if x <> "~" then o x;
        | WordNoSpace x -> 
            if x <> "~" then o x
        | Blank -> o " "
    in
    produce_troff data

let rec collapse = function
    | Blank
    | (WordNoSpace _)
    | (Word _) as x -> x
    | Text [(Text _) as y] -> collapse y
    | Text [(Word _) as y] -> y
    | Text x ->
            Text (List.map collapse x)
    | Command ("poy", _) -> Word "POY "
    | Command ("poybool", _) -> Word "BOOL "
    | Command ("poyfloat", _) -> Word "FLOAT "
    | Command ("poyint", _) -> Word "INTEGER "
    | Command ("poystring", _) -> Word "STRING "
    | Command ("poylident", _) -> Word "LIDENT "
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
            if b = "." || b = "," || b = ")" || b = "]" || b = "," || b = ";" then
                (WordNoSpace a) :: (collapse2 t)
            else h :: (collapse2 t)
    | h :: t -> h :: (collapse2 t)
    | [] -> []

let rec the_parser mode fstream = 
    let brstr = match mode with
        | `OnlineHelp -> "\n\n"
        | `Troff -> "\n.br\n"
    in
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
            let res = the_parser mode fstream in
            get_param (acc @ [(Text res)]) fstream
        end else acc
    in
    let skip1 = [' '; '\t'; '\\'; '\010'; '\013'; '{'; '}'; ','; '.'; ')'; '('; '['; ']'; '~'] 
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
            let data = List.rev acc in
            let () = if debug then List.iter print_parsed data in
            data
        else if is_blank fstream then
            split_on_commands (Blank :: acc) fstream
        else if is_close_command fstream then 
            List.rev acc
        else if is_newline fstream then
            split_on_commands ((Text [Word brstr]) :: acc) fstream
        else if is_command fstream then
            split_on_commands ((get_comm fstream) :: acc) fstream
        else if is_enclosed fstream then
            split_on_commands ((Text (the_parser mode fstream)) :: acc) fstream
        else if is_word fstream then
            split_on_commands ((get_word fstream) :: acc) fstream 
        else failwith "Huh?";
    in
    let res = split_on_commands [] fstream in
    collapse2 (List.map collapse res)

let process parsed mode channel =
    let generator = match mode with
        | `OnlineHelp -> produce_latex channel
        | `Troff      -> produce_troff channel
    in
    List.iter generator parsed

let process_file parsed mode output_file = 
    let channel = open_out output_file in
    let o str = output_string channel str in
    match mode with
    | `Troff ->
        o (".TH POY 1 LOCAL\n\n");
        o (".SH NAME\npoy \\- A phylogenetic analysis program using dynamic \
            homologies\n.SH SYNOPSIS\n.B poy [options] filename.\n\n\
            .SH DESCRIPTION\npoy is a phylogenetic analysis program for morphological \
            and molecular characters with support for dynamic homology characters: \
            that is, supports the analysis of unaligned sequences.\n\n.SH OPTIONS\n\n\
            .TP 5\n\
            -w\n\
            Run poy in the specified working directory. \n\
            .TP\n\
            -e\n\
            Exit upon error. \n\n\
            .TP\n\
            -d\n\
            Dump filename in case of error. \n\
            .TP\n\
            -q \n\
            Don't wait for input other than the program argument script. \n\
            .TP\n\
            -no-output-xml \n\
            Do not generate the output.xml file. \n\
            .TP\n\
            -help\n\
            Display this list of options. \n\
            .TP\n\
            --help\n\
            Display this list of options. \n\
            \n.SH VERSION\n 4.0." ^ (Str.global_replace (Str.regexp " ") "" BuildNumber.build) ^ "\n\
            .SH COMMANDS\n.P\n\
            For complete documentation go to \
            http://research.amnh.org/scicomp/projects/poy.php.\n\
            The following are the valid commands for \n\
            .B poy.");
        process parsed mode channel;
        o "\n.SH AUTHOR\nPOY was written by Andres Varon, Lin Hong, \
            NIcholas Lucaroni, and Ward Wheeler.\n.PP\nThis manual page was written \
            by Andres Varon, Lin Hong, Nicholas Lucaroni, Ward Wheeler\
            Ilya Temkin, Megan Cevasco, Kurt M. Pickett, Julian Faivovich, \
            Taran Grant, and William Leo Smith.\n.RS\n\n"
    | `OnlineHelp ->
        o (".\nauthors\n@[<v 2> POY was written by Andres Varon, Lin Hong, \
            Nicholas LUcaroni, and Ward Wheeler.This manual page was written \
            by Andres Varon, Lin Hong, Nicholas Lucaroni, Ward Wheeler\
            Ilya Temkin, Megan Cevasco, Kurt M. Pickett, Julian Faivovich, \
            Taran Grant, and William Leo Smith.@]\n\n");
        o (".\ncopyright\n@[<v 2>"^ Version.copyright_authors ^ Version.warrenty_information  ^"@]");
        process parsed mode channel;
        o ("\n\n.\nversion\n@[<v 2> "^ Version.copyright_authors ^ Version.compile_information ^ "@]");
        ()

let () = 
    let ch = FileStream.Pervasives.open_in (`Local "../doc/allcommands.tex") in 
    let parsed = the_parser `Troff ch in
    let () = process_file parsed `OnlineHelp "help.txt" in
    let () = process_file parsed `Troff      "poy.1" in
    ()
