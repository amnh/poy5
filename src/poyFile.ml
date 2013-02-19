(* POY 5.0 Beta. A phylogenetic analysis program using Dynamic Homologies.    *)
(* Copyright (C) 2013 Andrés Varón, Lin Hong, Nicholas Lucaroni, Ward Wheeler,*)
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

let () = SadmanOutput.register "PoyFile" "$Revision: 3160 $"

(* A poy file format magic number *)

exception InvalidFile

let magic_number1 = 0x6FD (* The french revolution begins! *)
let magic_number2 = 0x707 (* And the french revolution ends ... *)
let poy_major_version = 0x0004
let poy_minor_version = 0x0000
let data_structure_major_version = 0x0002
let data_structure_minor_version = 0x0014

(* [get_header ch] gets the poy header numbers. If the numbers can't be properly
* read, raise an InvalidFile exception. *)
let get_header ch =
    try
        let m1 = input_binary_int ch in
        let m2 = input_binary_int ch in
        let maj_v = input_binary_int ch in
        let min_v = input_binary_int ch in
        let dsmaj_v = input_binary_int ch in
        let dsmin_v = input_binary_int ch in
        m1, m2, maj_v, min_v, dsmaj_v, dsmin_v
    with
    | End_of_file ->
            raise InvalidFile

(* [output_header ch] outputs the POY magic and version numbers to the out
* channel [ch]. *)
let output_header ch =
    output_binary_int ch magic_number1;
    output_binary_int ch magic_number2;
    output_binary_int ch poy_major_version;
    output_binary_int ch poy_minor_version;
    output_binary_int ch data_structure_major_version;
    output_binary_int ch data_structure_minor_version

(* [is_valid_header (m1, m2, maj_v, min_v, dsmaj_v, dsmin_v)] checks weather or
* not the 6-tuple holds a set of magic numbers which is the same as those stored
* in this module (magic_number1, magic_number2, poy_major_version,
* poy_minor_version, data_structure_major_version,
* data_structure_minor_version). *)
let is_valid_header (m1, m2, maj_v, min_v, dsmaj_v, dsmin_v) =
    m1 = magic_number1 &&
    m2 = magic_number2 && 
    maj_v = poy_major_version &&
    min_v = poy_minor_version &&
    dsmaj_v = data_structure_major_version &&
    dsmin_v = data_structure_minor_version

(* [check_for_magic_number ch] true if the channel holds a valid POY file header
* and it's readable by this version of the module, false otherwise. Raises an
* InvalidFile error if the file does not contain a valid POY header. *)
let check_for_magic_number ch = 
    is_valid_header (get_header ch)

(* [read_file filename] outputs the contents of the POY file [filename]. If the
* file is not a valid POY format or it's not readable by this version of the
* module, raises InvalidFile. *)
let read_file filename = 
    let ch = open_in_bin filename in
    if check_for_magic_number ch then begin
        let res = Marshal.from_channel ch in
        close_in ch;
        res
    end else raise InvalidFile

(* [store_file data filename] outputs the [data] to the [filename] with the header
* values of the current module. *)
let store_file data filename =
    let ch = open_out_bin filename in
    output_header ch;
    Marshal.to_channel ch data [Marshal.Closures];
    close_out ch
