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

let convert_to_fasta file =
    let ch, file = FileStream.channel_n_filename file and
    definition_str = ref "" and
    accession_str = ref "" and
    outName, chout = Filename.open_temp_file "fasta" ".tmp" in
    try
        while true do
            let line = input_line ch in
            if Str.string_match (Str.regexp 
            " *\\(<GBSeq_accession-version>+\\)\\([a-zA-Z0-9. ,]+\\)") 
            line 0 then
                accession_str := 
                ((Str.matched_group 2 line) ^ "|"  ^ !definition_str)
            else if Str.string_match (Str.regexp 
            " *\\(<GBSeqid>gi|+\\)\\([a-zA-Z0-9 ,]+\\)") line 0 then
                output_string chout (">gi|" ^ (Str.matched_group 2 line) ^ 
                "|" ^ !accession_str)
            else if Str.string_match (Str.regexp 
            " *\\(<GBSeq_definition>+\\)\\([a-zA-Z0-9 ,]+\\)") line 0 then
                definition_str := (Str.matched_group 2 line) ^ "\n"
            else if Str.string_match (Str.regexp 
            " *\\(<GBSeq_sequence>+\\)\\([a-zA-Z]+\\)") line 0 then
                output_string chout ( (Str.matched_group 2 line) ^ "\n\n");
        done;
        open_in outName;
    with
    | End_of_file -> (open_in outName) 
