(* POY 5.0 Beta. A phylogenetic analysis program using Dynamic Homologies.    *)
(* Copyright (C) 2011  Andrés Varón, Lin Hong, Nicholas Lucaroni, Ward Wheeler*)
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

let anywhere_match regex line =
    try
        let _ = Str.search_forward regex line 0 in
        true
    with
    | _ -> false

let explode_filenames files = 
    let explode_filename file = 
        let file = FileStream.filename file in
        (* Warning: For some strange reason, if we run file_exists in
        * windows, then we can not catch the interrupt signal in the 
        * application. That suck for end users (ctr-C doesn't work). *)
        if Sys.os_type <> "Win32" && Sys.file_exists file then [file]
        else
            let ch = 
                let file = 
                    if Sys.os_type = "Win32" then file
                    else if Sys.os_type = "Unix" then
                        Str.global_replace (Str.regexp "\\\\ ") "\\ " file
                    else file
                in
                let line = 
                    match Sys.os_type with
                    | "Win32" -> ("dir /B \"" ^ file ^ "\" 2> NUL")
                    | _ -> "ls -1 " ^ file ^ " 2> /dev/null"
                in
                Unix.open_process_in line 
                in
            let res = IgnoreList.of_channel ch in
            close_in ch;
            match res with
            | [] -> 
                    let msg = "@[No@ file@ matching@ @{<b>" ^ StatusCommon.escape file ^ 
                    "@}@ found.@]" in
                    Status.user_message Status.Error msg;
                    failwith "File not found"
            | _ -> res
        in
        List.flatten (List.map explode_filename files)

