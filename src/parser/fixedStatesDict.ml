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

let of_channel t ch = 
    let rec reader lexer alph acc counter =
        try
            let line = input_line ch in
            let acc = (Fasta.process_sequence true lexer alph [line]) :: acc in
            reader lexer alph acc (counter + 1)
        with
        | End_of_file -> acc
    in
    match t with
    | FileContents.Nucleic_Acids -> 
            let lexer = Alphabet.Lexer.make_lexer true false Alphabet.nucleotides in
            reader lexer Alphabet.nucleotides [] 1
    | FileContents.Proteins -> 
            let msg = "The current list of protein codes is empty. \
            Request an update or contact the POY maintainer." in
            raise (E.Unsupported_file_format msg)
    | FileContents.AlphSeq x ->
            let lexer = Alphabet.Lexer.make_lexer true true x in
            reader lexer x [] 1 
    | _ -> 
            let msg = 
                "Unknown file format for a Fixed states dictionary. "
            in
            raise (E.Unsupported_file_format msg)

let of_file t f =
    let ch = FileStream.open_in f in
    let res = of_channel t ch in
    close_in ch;
    res

module OderedSequence = struct
    type t = Sequence.s
    let compare = Sequence.compare
end

module SeqMap = Map.Make (OderedSequence)

let create_dp_read_file ?filename file taxa lst tcm = 
    (* We report to the user what we intend to do *)
    let outName, ch = 
        match filename with
        | Some v -> v, open_out v
        | None -> Filename.open_temp_file file ".tmp" 
    in
    let len = List.length taxa in
    (* Print the header of the file *)
    Printf.fprintf ch "dpread\n";
    Printf.fprintf ch "\'Poy file generated for fixed states analysis \
    from\nthe original file %s.\'\n" file;
    Printf.fprintf ch "1 %d\n" len;
    (* We first load the fixed states *)
    let res = List.fold_left (fun (code, acc) seq ->
        match SeqMap.mem seq acc with
        | false -> code + 1, SeqMap.add seq code acc
        | true -> code, acc) (0, SeqMap.empty) lst 
    in
    (* Now we load the observations and print them out in the file *)
    let len, map = List.fold_left (fun (code, acc) (seq, taxon) ->
        match SeqMap.mem seq acc with
        | false -> 
                Printf.fprintf ch "%s %d\n" taxon code;
                code + 1, SeqMap.add seq code acc 
        | true -> 
                let tc = SeqMap.find seq acc in
                Printf.fprintf ch "%s %d\n" taxon tc;
                code, acc) res taxa 
    in
    (* Print the transformation cost matrix *)
    Printf.fprintf ch ";\ncost [ 0 $ %d\n" len;
    let asc_lst = SeqMap.fold (fun a b acc -> (b, a) :: acc) map [] in
    let asc_lst = List.sort (fun (a, _) (b, _) -> a - b) asc_lst in
    (* A function to print out the necessary header of the transformation
    * cost matrix in the dpread file *)
    let rec print_table_header x =
        if x < len then begin
            Printf.fprintf ch "%d " (x mod 10);
            print_table_header (x + 1)
        end else 
            Printf.fprintf ch "\n";
    in
    print_table_header 0;
    let len = List.length asc_lst in
    let st = Status.create "Parser" (Some len) ("Converting " ^ file ^ 
    " to dpread file \"" ^ outName ^ "\"") in
    List.iter (fun (_, x) ->
        List.iter (fun (_, y) ->
            let c = Sequence.Align.cost_2 x y tcm Matrix.default in
            Printf.fprintf ch "%d " c;) asc_lst;
            (* Report to the user *)
            Status.full_report ~adv:(1 + Status.get_achieved st) st;
            Printf.fprintf ch "\n";) asc_lst;
    Status.finished st;
    close_out ch;
    outName
