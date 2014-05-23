(* POY 5.1.1. A phylogenetic analysis program using Dynamic Homologies.       *)
(* Copyright (C) 2014 Andrés Varón, Lin Hong, Nicholas Lucaroni, Ward Wheeler,*)
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

(** Type of data contained in a file (either Fasta, POY or Hennig86
 * formats *)
type t =
    | Nucleic_Acids
    | Proteins
    | Prealigned_Alphabet of Alphabet.a
    | AlphSeq of Alphabet.a
    | Inactive_Character
    | Ordered_Character of (int * int * bool)
    | Unordered_Character of (int * bool)
    | Sankoff_Character of (int list * bool)
    | Genes of int array

let to_string x =
    let pp_ilst () x =
        List.fold_left (fun acc x -> acc ^ "," ^ (string_of_int x)) "" x
    and pp_iray () x =
        Array.fold_left (fun acc x -> acc ^ "," ^ (string_of_int x)) "" x in
    match x with
    | Nucleic_Acids -> "dna"
    | Proteins      -> "aa"
    | Prealigned_Alphabet a -> "Pre"
    | AlphSeq a             -> "Seq"
    | Inactive_Character    -> "inactive"
    | Ordered_Character (a,b,c) -> Printf.sprintf "ord:(%d,%d,%B)" a b c
    | Unordered_Character (a,b) -> Printf.sprintf "unord:(%d,%B)" a b
    | Sankoff_Character (is,b)  -> Printf.sprintf "sank:(%a,%B)" pp_ilst is b
    | Genes ray                 -> Printf.sprintf "gene:(%a)" pp_iray ray

