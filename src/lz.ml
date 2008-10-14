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

type protocol = NoCompression | Compressed of int
(** 
* NoCompression has no compression turned on.
* Compressed uses compression in different versions:
*   1 Has a bug in the compression and so it doesn't really compress very well,
*   moreover, we have to fix the basic procedure to decompress to workaround a
*   bug. Compression version 1 starts with a 0 byte.
*   2 The previous bug was corrected. This kind of files start with a 1 byte.
**)


type table = { 
    encoding_table : (compress, int) Hashtbl.t; (* the table *)
    decoding_table : (int, compress) Hashtbl.t;
    cnt : int ref; (* the next available code *)
    state : (int option) ref;
    sticky_state : (int option) ref; (* For bug workaround in Compressed 1 *)
    reached_limit : bool ref;
    did_exit : int ref;
}

let initial_table () = 
    let encoding_table = Hashtbl.create 99991 
    and decoding_table = Hashtbl.create 99991 in
    for i = 0 to 255 do
        Hashtbl.add encoding_table (Single i) i;
        Hashtbl.add decoding_table i (Single i);
    done;
    { 
        encoding_table = encoding_table; 
        decoding_table = decoding_table; 
        cnt = ref 256; 
        state = ref None;
        sticky_state = ref None;
        reached_limit = ref false;
        did_exit = ref 0;
    }

(* Always compress with the latest protocol version *)
let latest = Compressed 2

let add_to_tables protocol table code =
    (* Because of the way we print files, using only shorts, we can support a
    * dictionary with size at most 0xFFFF, therefore, we check if we have
    * reached the limit, in which case we don't do anything *)
    if !(table.cnt) > 0xFFFF then 
        let () =
            if not !(table.reached_limit) then begin
                table.reached_limit := true;
                table.sticky_state := !(table.state);
            end else ()
        in
        match protocol with
        | Compressed 1 | NoCompression -> ()
        | Compressed 2 -> table.state := None
        | Compressed x -> 
                failwith ("Unknown protocol format " ^ string_of_int x)
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
                        add_to_tables latest table potential_state;
                        let acc = code :: previous :: acc in
                        aux_compress (pos + 1) acc
                    end 
        else acc
    in
    aux_compress initial acc


let decode protocol table string h =
    let rec aux_decode stack table string h = 
            match Hashtbl.find table.decoding_table h with
            | Single i -> i :: stack
            | Pair (prev, i) ->
                    aux_decode (i :: stack) table string prev;
    in
    let handle_encoding i =
        let () = 
            match !(table.state) with
            | None ->  table.state := Some i;
            | Some previous ->
                    let potential = Pair (previous, i) in
                    if Hashtbl.mem table.encoding_table potential then begin
                        if !(table.cnt) > 0xFFFF && protocol = Compressed 1 then
                            begin
                                incr table.did_exit;
                            if !(table.did_exit) = 20 then raise Exit;
                        end;
                        table.state := 
                            Some (Hashtbl.find table.encoding_table potential);
                        table.sticky_state := None;
                        table.reached_limit := false;
                    end else add_to_tables protocol table potential;
        in
        Buffer.add_char string (Char.chr i)
    in
    if (protocol = Compressed 1) && (Some h) = !(table.sticky_state) then ()
    else List.iter handle_encoding (aux_decode [] table string h)

let rec decompress protocol table compressed string = 
    List.iter (decode protocol table string) compressed

let full_compress string = 
    let table = initial_table () in
    let acc = compress table string [] in
    match !(table.state) with
    | None -> List.rev acc
    | Some x -> List.rev (x :: acc)

let full_decompress protocol lst = 
    let table = initial_table () in
    let res = Buffer.create 200 in
    decompress protocol table lst res;
    Buffer.contents res

let magic_number1 = 0x6FD (* The french revolution begins! *)
let magic_number2 = 0x707 (* And the french revolution ends ... *)
let poy_major_version = 0x0004
let poy_minor_version = 0x0000
let protocol_major_version = 0x0000
let protocol_minor_version = 0x0000

let valid_headers = 
    [(magic_number1, magic_number2, poy_major_version, poy_minor_version, 
        protocol_major_version, protocol_minor_version), Compressed 2 ]

let valid_protocols = List.map (fun (a, b) -> b, a) valid_headers

let output_header protocol ch =
    let a, b, c, d, e, f = List.assoc protocol valid_protocols in
    output_byte ch 0x1;
    output_binary_int ch a;
    output_binary_int ch b;
    output_binary_int ch c;
    output_binary_int ch d;
    output_binary_int ch e;
    output_binary_int ch f


let detect_type ch = 
    match Pervasives.input_byte ch with
    | 0 -> Compressed 1
    | 1 -> 
            (let the_header = PoyFile.get_header ch in
            try List.assoc the_header valid_headers with
            | Not_found -> raise PoyFile.InvalidFile)
    | x -> NoCompression

let skip_header ch protocol =
    match protocol with
    | NoCompression | Compressed 1 -> ()
    | Compressed 2 ->
            (* We get rid of the protocol flag, and then the header itself *)
            let _ = input_byte ch in
            let _ = PoyFile.get_header ch in
            ()
    | _ -> ()
