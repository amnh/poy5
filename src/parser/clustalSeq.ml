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

let process_line filename hash taxa line_number line =
    if 0 = String.length line || ' ' = line.[0] then 
        taxa
    else 
        match Str.split (Str.regexp " +") line with
        | [taxon; sequence] -> 
                Hashtbl.add hash taxon sequence;
                if All_sets.Strings.mem taxon taxa then taxa
                else All_sets.Strings.add taxon taxa 
        | _ -> 
                raise (E.Unexpected
                ("Unexpected line " ^ string_of_int line_number ^ 
                " in Clustal file " ^ filename ^ ": " ^ line))

let merge_lines hash ch taxon = 
    let seqs  = Hashtbl.find_all hash taxon in
    output_string ch ">";
    output_string ch taxon;
    output_string ch "\n";
    List.fold_right (fun x () ->
        output_string ch x; ()) 
    seqs ();
    output_string ch "\n\n"

let convert_to_fasta file =
    let ch, file = FileStream.channel_n_filename file in
    let first_line = input_line ch 
    and outName, chout = Filename.open_temp_file "fasta" ".tmp" in
    if "CLUSTAL" = String.sub first_line 0 7 then
        let hash = Hashtbl.create 97 
        and line_number = ref 0 
        and taxa = ref All_sets.Strings.empty in
        try
            while true do
                let nt = 
                    process_line file hash !taxa !line_number
                    (input_line ch)
                in
                taxa := nt;
                incr line_number;
            done;
            open_in outName
        with
        | End_of_file ->
                All_sets.Strings.iter (merge_lines hash chout) !taxa;
                close_out chout;
                open_in outName
    else 
        raise (E.Unexpected 
        ("The first line of a fasta file has to begin with the word \
        CLUSTAL"))

