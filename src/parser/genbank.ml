(* POY 5.1.1. A phylogenetic analysis program using Dynamic Homologies.       *)
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

let convert_to_fasta ?filename file =
    let ch, file = FileStream.channel_n_filename file 
    and outName, chout = 
    match filename with
    | None -> Filename.open_temp_file "fasta" ".tmp"
    | Some f -> f, open_out f
in
    try
        while true do
            let line = input_line ch in
            (*
            if Str.string_match (Str.regexp 
            " *\\(VERSION+ *\\) \\([a-zA-Z0-9.]+\\) *GI:\\([0-9]+\\)") 
            line 0 then 
                output_string chout 
                (">gi|" ^ (Str.matched_group 3 line) ^ "|" ^
                (Str.matched_group 2 line) ^ "|" ^ !definition_str)
            *)
            if Str.string_match (Str.regexp 
            " *\\(ORGANISM+ *\\) \\([a-zA-Z0-9 ,]+\\)") line 0 then begin
                let str = Str.matched_group 2 line in
                let str = Str.global_replace (Str.regexp " ") "_" str in
                output_string chout ">";
                output_string chout str;
                output_string chout "\n";
            end else
                if Str.string_match 
                (Str.regexp " *\\([0-9]+\\) \\([a-zA-Z ]+\\)+") line 0 then
                    let temp = Str.matched_group 2 line in
                    let temp2 = Str.global_replace (Str.regexp "[ ]+")
                    "" temp in
                    output_string chout (String.uppercase temp2) ;
                    output_string chout "\n";
        done;
        open_in outName;
    with
    | End_of_file -> (open_in outName) 

