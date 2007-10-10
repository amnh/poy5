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

(* $Id: helpRead.ml 1644 2007-02-14 19:05:47Z andres $ *)
(* Created Tue Apr 25 13:55:43 2006 (Illya Bomash) *)

(** [HelpRead] contains functions to parse our help file format and generate
    the structures for our online help.  It also contains functions to output
    to HTML. *)

let debug = true

let do_html = ref false
let () =
    Arg.parse
        [ ("-html", Arg.Set do_html, "Output help.html from help.txt") ]
        (fun _ -> ())
        "Process POY's help.txt file into POY source code"

module OrderedString = struct
    type t = string 
    let compare a b = Pervasives.compare a b
end
module Strings = Set.Make (OrderedString)

let build_cross_ref prev item =
    prev ^ "@ " ^ "@{<u>@{<help-xref:" ^ item ^ ">" ^ item ^ "@}@}"

let replace str = Str.global_replace (Str.regexp " ") "@ " str
let another_replace str = Str.global_replace (Str.regexp "<v@ 2>") "<v 2>" str
let yet_another_replace str = Str.global_replace (Str.regexp "%@") "%" str

let replace str = yet_another_replace (another_replace (replace str))


let make_help_line title usage description cross_ref =
    "@[<v 2>@{<b>" ^ title
    ^ (if usage = "" then "" else " " ^ usage)
    ^ ":@}@ @[<hov 2>" ^ replace description ^ 
    "@]" ^ 
    (match cross_ref with [] -> "" | _ -> "@,@[<v 2>See also:" ^ List.fold_left
    build_cross_ref "" cross_ref ^ "@]") ^ "@]@,"

let rec read_help_file_structure ?acc file =
    let input_line ?(skip_blank=false) file =
        let rec il file =
            let l = input_line file in
            if 0 = String.compare l ""
            then (if skip_blank then il file else l)
            else if l.[0] = '#'
            then il file
            else l in
        il file in
    let acc = match acc with
    | Some a -> a
    | None -> [] in
    try
        let name = input_line ~skip_blank:true file in
        let desc = input_line file in
        let desc = if 0 = String.compare desc "." then "" else desc in
        let xref = ref [] in
        let rec get file =
            let l = input_line file in
            if l = "" then []
                (* TODO: xref support *)
            else if l.[0] = '['
            then begin
                let split = Str.regexp ";[ \t]*" in
                let l = Str.global_replace (Str.regexp "^\\[[ \t]*") "" l in
                let l = Str.global_replace (Str.regexp "[ \t]*\\]$") "" l in
                xref := Str.split split l;
                (get file)
            end
            else l :: (get file) in
        let lines = get file in
        let text = String.concat "\n" lines in
        (* TODO: replace with "@ " *)
        read_help_file_structure ~acc:((name, desc, text, !xref) :: acc)
            file
    with _ -> acc

let read_help_file file =
    let map (name, desc, text, xref) =
        (name, make_help_line name desc text xref) in
    (* ignore first line *)
    let _ = input_line file in
    let structure = read_help_file_structure file in

    if debug then begin
        let all_nodes =
            List.fold_left
                (fun set (_, _, _, xrefs) ->
                     List.fold_left
                         (fun set xref -> Strings.add xref set)
                         set xrefs)
                Strings.empty structure in
        (* check whether all the nodes exist *)
        let pr_exists n =
            if List.exists (fun (name, _, _, _) -> n = name) structure
            then ()
            else print_endline ("No help entry for " ^ n)
        in
        Strings.iter pr_exists all_nodes
    end;

    List.map map structure

let helpfile = open_in "help.txt"
let index : (string * string) list = List.rev (read_help_file helpfile)
let () = close_in helpfile

(* Output to help.ml *)
let () = match !do_html with
| false ->
      let ch = open_out "help.ml" in
      let print_endline s = output_string ch s; output_string ch "\n"; flush ch in
      print_endline "(** [help.ml] is an automatically-generated source file.";
      print_endline "    It is generated by helpRead.ml from help.txt. *)";
      print_endline "";
      print_endline "let index = [";
      List.iter
          (fun (f, t) ->
               print_endline ("    (\"" ^ String.escaped f ^ "\", \""
                              ^ String.escaped t ^ "\");"))
          index;
      print_endline "]"
| true ->
      let ch = open_out "help.html" in
      let print_endline s = output_string ch s; output_string ch "\n"; flush ch
      in
      let command_re = Str.regexp "^\\([^:]*\\): " in
      let subs =
          [ ("@\\[<[^>]*>", "");
            ("@]", "");
            ("@{<help-xref:\\([^>]*\\)>\\([^}]*\\)@}", "<a href=\"#\\1\">\\2</a>");
            ("@{<[^>]*>", "");
            ("@}", "");
            ("@,", "\n<p />");
            ("\n", "\n<p />");
            ("POY", "P<font size=\"-1\">OY</font>");
          ] in
      let subs = List.map
          (fun (a, b) -> (Str.regexp a, b)) subs in
      let unfilter (t : string) : string =
          List.fold_left
              (fun t (filter, rep) ->
                   Str.global_replace filter rep t)
              t subs in
      let break (f, t) =
          let oldt = t in
          let t = Str.global_replace (Str.regexp "@ ") " " (unfilter t) in
          let title =
              if Str.string_match command_re t 0
              then Str.matched_group 1 t
              else f in
          let t = Str.replace_first command_re "" t in
          f, title, t ^ "\n<!-- " ^ oldt ^ " -->" in
      print_endline "<html>";
      print_endline "<head>";
      print_endline "<title>POY cheat sheet</title>";
      print_endline "<style type=\"text/css\">";
      print_endline "  body {";
      print_endline "  width: 600px;";
      print_endline "  margin-left: 3em;";
      print_endline "  text-align: justify;";
      print_endline "  font-family: sans;";
      print_endline "  }";
      print_endline "  h1 {";
      print_endline "  background: #aaf;";
      print_endline "  }";
      print_endline "  h2 {";
      print_endline "  background: #eef;";
      print_endline "  text-align: left;";
      print_endline "  font-size: 14pt;";
      print_endline "  font-weight: bold;";
      print_endline "  margin-top: 2em;";
      print_endline "  }";
      print_endline "  p {";
      print_endline "  }";
      print_endline "</style>";
      print_endline "</head>";
      print_endline "<body>";
      print_endline "<h1>POY Cheat Sheet</h1>";
      List.iter
          (fun entry ->
               let short, title, descr = break entry in
               print_endline ("<h2><a name=\""
                              ^ short
                              ^ "\">"
                              ^ title ^ "</a></h2>");
               print_endline ("<p />" ^ descr))
          index;
      print_endline "</body>";
      print_endline "</html>"
      
