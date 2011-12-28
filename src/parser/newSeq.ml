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

let to_fasta file =
    let ch, file = FileStream.channel_n_filename file and
    start_string = ref ">" and
    outName, chout = Filename.open_temp_file "fasta" ".tmp" in
    try
        while true do
            let line = input_line ch in
            if Str.string_match (Str.regexp 
            "\\(^[^ ]+\\)") line 0 then 
                output_string chout 
                (!start_string ^ (Str.matched_group 1 line) ^ "\n")
            else
                if Str.string_match 
                (Str.regexp " *\\([0-9]+\\) *\\([a-zA-Z ]+\\)+") line 0 then
                    let temp = Str.matched_group 2 line in
                    let temp2 = Str.global_replace (Str.regexp "[ ]+")
                    "" temp in
                    output_string chout ((String.uppercase temp2) ^ "\n");
                    start_string := "\n>"
        done;
        close_out chout;
        outName
    with
    | End_of_file -> 
            close_out chout;
            outName

let convert_to_fasta file =
    let file = to_fasta file in
    open_in file

