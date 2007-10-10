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

let () = SadmanOutput.register "HashedFuncs" "$Revision: 1644 $"

let same_10000 _ = 
    log 10000.0

let same_100 _ =
    log 100.0

let defined_functions = [
    ("same_10000", same_10000);
    ("same_100", same_100)
]

let index = ref (Hashtbl.create 10)

let add name f = Hashtbl.add !index name f 

let remove name = Hashtbl.remove !index name

let find name = 
    try 
        Hashtbl.find !index name 
    with 
    | Not_found as e ->
            prerr_string "HashedFuncs.find";
            raise e

let _ =
    List.iter (fun (x, y) -> add x y ) defined_functions
