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

exception Illegal_Grammar
exception Name_Collision of (CharacSpec.t * string)
type range = int * int
type distr = Poisson of float
(* The types returned by the parser to express the functions in the U 
* language *)
type ulang =
  | CInt of int * string
  | CString of string * string
  | CBase of string * string
  | Variable of string
  | Predecessor of ulang
  | Successor of ulang
  | Prepend of ulang * ulang
  | Tail of ulang
  | Head of ulang
  | Bool of string * ulang
  | IfThen of ulang * ulang * ulang
  | Application of string * ulang list

(* A function in the U language is composed of he following:
    * Name
    * Parameters with their corresponding type 
    * The function definition 
*)
type ufuncs = string * (string * string) list * ulang

type filename = string

(* A specification in the input file. *)
type spec =
    | Alphabet of string * string list * spec option 
        (** An alphabet specification with [(name, items, probability
        * distribution)] *)
    | WordSet of string * string * range * spec option
        (** A word set specification with [(name, alphabet, length range, 
        * probablity distribution)] *)
    | IntSet of string * range * spec option
        (** An integer set specification with [(name, valid range, probability
        * distribution)]. *)
    | FProbability of distr
        (** A function probability with [(name of prob, name of function)]. *)
    | EProbability of string * (string * float) list
        (** An enumeration probablity with [(name of probability, pairs of items
        * and probability)] where each item is a member of the enumeration. *)
    | Character of string * ufuncs list * spec option
        (** A character specification, with [(name, functions, function
            * probabilit)]. *)
    | Protein_Sequences of (filename list * (Data.dynhom_opts list option)) list
        (** A protein sequence dataset defined in the same way as Dna_Sequences.
        * *)
    | Synonyms of string list
        (** A list of files containing taxon synonyms *)
    | SynonymsList of (string * string) list
        (** A list of pairs of synonyms, where [(a, b)] means that [a] should be
        * treated as taxon [b]. *)
    | Static of filename list
        (** A list of files containing static homology character types. *)
    | Ignore_files of string list
        (** A list of files with lists of taxa to be ignored. *)
    | Ignore of string list
        (** A list of taxa to be ignored. *)
    | Trees of string list
    | TreesList of string list

(* [process_tree lst] process the specifications defined in [lst] for an 
* analysis. The output is a triple [(a, b, c)] where [a] is the index of
* specifications that define [b], the list of characters defined in the file,
* where each item [b'] in [b] is a tuple with the name and the specification of
* a character. Finally, [c] is the actuall data loaded in the file. *)
val process_tree : Data.d -> spec list -> Data.d

val of_file : Data.d -> string -> Data.d

(* [of_channel c] parses and checks the grammar of a Poy script loaded from [c].
* The output is valid input for [process_tree]. *)
val of_channel : Pervasives.in_channel -> spec list

val guess_class_and_add_file : bool -> bool -> Data.d -> Parser.filename -> Data.d

val explode_filenames : [`Local of string | `Remote of string ] list -> string list
