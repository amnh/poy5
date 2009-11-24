(* POY 4.0 Beta. A phylogenetic analysis program using Dynamic Homologies.    *)
(* Copyright (C) 2007  Andr�s Var�n, Le Sy Vinh, Illya Bomash, Ward Wheeler,  *)
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

let () = SadmanOutput.register "Arguments" "$Revision: 1959 $"

let just_exit = ref false
let only_run_argument_script = ref false
let dump_file = ref "ft_output.poy"
let input : [`Inlined of string | `Filename of string] list ref = ref []

let change_working_directory str =
    try Sys.chdir str with
    | err ->
            prerr_endline 
            ("Attemting to change to directory " ^ str ^ 
            " failed misserably. I am cancelling this run, check your " ^
            "path.\n");
            exit 1

let process_poy_plugin plugin =
IFDEF USE_NATIVEDYNLINK THEN
    if Sys.file_exists plugin then 
        let extension = if Dynlink.is_native then "cmxs" else "cmo" in
        if Filename.check_suffix plugin extension then Dynlink.loadfile plugin
        else 
            failwith ("A plugin for this executable must have extension " ^
            extension)
    else failwith ("Could not find KML plugin " ^ plugin)
ELSE
    failwith "This version of POY was compiled without plugin support"
END

let anon_fun kind str =
    match kind with
    | `Inlined -> input := (`Inlined str) :: !input
    | `Filename -> input := (`Filename str) :: !input


let parse_list = [
    ("-w", Arg.String change_working_directory,
    "Run poy in the specified working directory"); 
    ("-e", Arg.Unit (fun () -> just_exit := true), "Exit upon error");
    ("-d", Arg.String (fun str -> dump_file := str), "Dump filename in case of \
    error");
    ("-q", Arg.Unit (fun () -> only_run_argument_script := true), 
    "Don't wait for input other than the program argument script.");
    ("-no-output-xml", Arg.Unit (fun () -> SadmanOutput.do_sadman := false),
    "Do not generate the output.xml file.");
    ("-plugin", Arg.String process_poy_plugin, 
    "Load the selected plugins.");
    ("-script", Arg.String (anon_fun `Inlined), 
    "Inlined script to be included in the analysis.")
]

let usage = 
    "poy [OPTIONS] filename ..."
