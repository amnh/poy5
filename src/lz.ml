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

type compress = Single of int | Pair of (int * int)

type table = { 
    encoding_table : (compress, int) Hashtbl.t; (* the table *)
    decoding_table : (int, compress) Hashtbl.t;
    cnt : int ref; (* the next available code *)
    state : (int option) ref;
}

let initial_table () = 
    let encoding_table = Hashtbl.create 1667 
    and decoding_table = Hashtbl.create 1667 in
    for i = 0 to 255 do
        Hashtbl.add encoding_table (Single i) i;
        Hashtbl.add decoding_table i (Single i);
    done;
    { 
        encoding_table = encoding_table; 
        decoding_table = decoding_table; 
        cnt = ref 256; 
        state = ref None 
    }

let add_to_tables table code =
    (* Because of the way we print files, using only shorts, we can support a
    * dictionary with size at most 0xFFFF, therefore, we check if we have
    * reached the limit, in which case we don't do anything *)
    if !(table.cnt) > 0xFFFF then ()
    else begin
        Hashtbl.replace table.encoding_table code !(table.cnt); 
        Hashtbl.replace table.decoding_table !(table.cnt) code; 
        incr table.cnt;
        table.state := None;
    end

let rec compress ?initial ?length table string acc =
    let initial =
        match initial with
        | None -> 0
        | Some x -> x
    and max =
        match length with
        | None -> String.length string 
        | Some x -> x
    in
    let rec aux_compress pos acc =
        if pos < max then
            let code = Char.code string.[pos] in
            match !(table.state) with
            | None -> 
                    table.state := Some code;
                    aux_compress (pos + 1) acc
            | Some previous ->
                    let potential_state = Pair (previous, code) in
                    if Hashtbl.mem table.encoding_table potential_state then
                        begin
                        table.state := 
                            Some (Hashtbl.find table.encoding_table
                            potential_state);
                        aux_compress (pos + 1) acc
                    end else begin
                        add_to_tables table potential_state;
                        let acc = code :: previous :: acc in
                        aux_compress (pos + 1) acc
                    end 
        else acc
    in
    aux_compress initial acc

let rec decode table string h = 
    let handle_encoding i =
        let () = 
            match !(table.state) with
            | None ->  table.state := Some i;
            | Some previous ->
                    let potential = Pair (previous, i) in
                    if Hashtbl.mem table.encoding_table potential then
                        table.state := 
                            Some (Hashtbl.find table.encoding_table potential)
                    else add_to_tables table potential;
        in
        i
    in
    let i =
        match Hashtbl.find table.decoding_table h with
        | Single i -> handle_encoding i
        | Pair (prev, i) ->
                decode table string prev;
                handle_encoding i
    in
    Buffer.add_char string (Char.chr i)

let rec decompress table compressed string = 
    List.iter (decode table string) compressed

let full_compress string = 
    let table = initial_table () in
    let acc = compress table string [] in
    match !(table.state) with
    | None -> List.rev acc
    | Some x -> List.rev (x :: acc)

let full_decompress lst = 
    let table = initial_table () in
    let res = Buffer.create 200 in
    decompress table lst res;
    Buffer.contents res
