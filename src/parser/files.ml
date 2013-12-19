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

(** [is_unknown t] tells observers whether this present character was listed as
    unknown by the user *)
let is_unknown t = match t with
    | FileContents.Nucleic_Acids
    | FileContents.Proteins
    | FileContents.Prealigned_Alphabet _
    | FileContents.Genes _
    | FileContents.AlphSeq _
    | FileContents.Inactive_Character -> false
    | FileContents.Unordered_Character (_, bool)
    | FileContents.Sankoff_Character (_, bool)
    | FileContents.Ordered_Character (_, _, bool) -> bool

type ft = 
    | Is_Hennig 
    | Is_Dpread
    | Is_Clustal
    | Is_Fasta 
    | Is_Poy 
    | Is_Genome 
    | Is_ASN1 
    | Is_Genbank 
    | Is_INSDSeq  
    | Is_GBSeq 
    | Is_TinySeq 
    | Is_XML 
    | Is_Nexus 
    | Is_NewSeq 
    | Is_Dictionary
    | Is_Fixed_States_Dictionary
    | Is_Phylip
    | Is_Unknown 
    | Is_Trees
    | Is_ComplexTerminals


let poy_file_regex = Str.regexp "\\(define\\|load\\|trees\\|names\\|ignore\\)"


let is_hennig file =
    let ch = FileStream.Pervasives.open_in file in
    try 
        while true do
            let line = FileStream.Pervasives.input_line ch in
            if String.length line > 4 then
                let lst = [ line.[0]; line.[1]; line.[2]; line.[3]; line.[4] ] in
                let lst = List.map Char.uppercase lst in
                (* Does the line start with xread or tread? *)
                if lst = ['X';'R';'E';'A';'D'] || lst = ['T';'R';'E';'A';'D'] then 
                    raise Exit;
        done;
        assert false
    with
    | Exit -> true
    | End_of_file -> false


(* Some general utilities for the parser *)
let test_file file = 
    let ch = FileStream.Pervasives.open_in file in
    let line = FileStream.Pervasives.input_line ch in
    let line2 = 
        try FileStream.Pervasives.input_line ch with
        | _ -> ""
    in
    FileStream.Pervasives.close_in ch;
    if Wildcard.anywhere_match (Str.regexp "^CLUSTAL") line then Is_Clustal
    else if is_hennig file then Is_Hennig
    (* treat dpread as a hennig file *)
    else if Wildcard.anywhere_match (Str.regexp "COMPLEX") line then
        Is_ComplexTerminals
    else if Wildcard.anywhere_match (Str.regexp "dpread") line then Is_Dpread
    else if Wildcard.anywhere_match poy_file_regex line then
        Is_Poy
    else if Wildcard.anywhere_match (Str.regexp "^>") line then 
        begin
            if Wildcard.anywhere_match (Str.regexp "[a-zA-Z]") line2 then
                Is_Fasta 
            else
                Is_Genome
        end
    else if Wildcard.anywhere_match (Str.regexp "Seq") line then Is_ASN1
    else if Wildcard.anywhere_match (Str.regexp "LOCUS") line then Is_Genbank
    else if Wildcard.anywhere_match (Str.regexp "#NEXUS") (String.uppercase line) then Is_Nexus
    else if Wildcard.anywhere_match (Str.regexp "<\\?xml") line then
        begin
            if Wildcard.anywhere_match (Str.regexp "<!DOCTYPE INSDSeq") line2 then
                Is_INSDSeq
            else if Wildcard.anywhere_match (Str.regexp "<!DOCTYPE GB") 
            line2 then
                Is_GBSeq
            else if Wildcard.anywhere_match (Str.regexp "<!DOCTYPE TSeq") 
            line2 then
                Is_TinySeq
            else 
                Is_XML
        end
    else if Wildcard.anywhere_match (Str.regexp "^[0-9]+[ \t]*[0-9]+") line then
        Is_Phylip
    else if Wildcard.anywhere_match (Str.regexp "^[a-zA-Z0-9_]+") line 
            && Wildcard.anywhere_match (Str.regexp " *[0-9]+") line2 then
                Is_NewSeq
    else if Wildcard.anywhere_match (Str.regexp "^[a-zA-Z._-]+ +[a-zA-Z._-]+\\s*") line then
        Is_Dictionary
    else if Wildcard.anywhere_match (Str.regexp "^ *(") line then
        Is_Trees
    else Is_Unknown

let molecular_to_fasta file =
    match test_file file with
    | Is_Clustal -> ClustalSeq.convert_to_fasta file
    | Is_GBSeq -> GBSeq.convert_to_fasta file
    | Is_TinySeq -> TinySeq.convert_to_fasta file
    | Is_INSDSeq -> INSDSeq.convert_to_fasta file
    | Is_XML -> XMLGB.convert_to_fasta file
    | Is_ASN1 -> Asn1.convert_to_fasta file
    | Is_Genbank -> Genbank.convert_to_fasta file
    | Is_NewSeq -> NewSeq.convert_to_fasta file
    | _ -> FileStream.open_in file

