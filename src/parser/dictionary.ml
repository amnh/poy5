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

let gen_of_channel adder ch = 
    let input_handler = FileStream.stream_reader ch in
    let rec reader counter =
        try let line = input_handler#read_line in
            match Str.split (Str.regexp "[\t ]+") line with
            | [] | [_] -> 
                let msg = ("Line " ^ string_of_int counter ^ ": " ^ line) in
                raise (E.Illegal_dictionary msg)
            | hd :: tl ->
                List.iter (fun x -> adder x hd) tl;
                reader (counter + 1)
        with | End_of_file -> ()
    in
    reader 1

let of_channel ch =
    let table = Hashtbl.create 1667 in
    gen_of_channel (Hashtbl.add table) ch;
    table

let of_channel_assoc ch =
    let table = ref [] in
    gen_of_channel (fun a b -> table := (a, b) :: !table) ch;
    List.rev !table

