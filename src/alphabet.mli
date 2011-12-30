(* POY 5.0 Alpha. A phylogenetic analysis program using Dynamic Homologies.   *)
(* Copyright (C) 2011 Andrés Varón, Lin Hong, Nicholas Lucaroni, Ward Wheeler,*)
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

(* $Id: alphabet.mli 2845 2008-05-14 15:12:35Z vinh $ *)
(* Alphabet.
*
* Description and handling of different kinds of alphabets for the analysis of
* molecular data. *)

(** Alphabet specification for sequences.
*
* Sequences are usually composed of nucleotides or aminoacides. Other kind of
* sequences are valid though. This library provides the interface for such
* specifications, as well as the default ones (nucleotides and aminoacids). *)

(** {2 Exceptions} *)

exception Illegal_Character of string
exception Illegal_Code of int

(** {2 Types} *)

(** The type of the alphabet, either a sequence of integers, up to 32 elements,
    each represented by a bit, or each represented by a bit including all the
    combinations. *)
type kind = 
    | Sequential
    | Simple_Bit_Flags
    | Extended_Bit_Flags

(** Alphabet type *)
type a

(** {2 Alphabets} *)

(** [present_absent] contains binary data; used to calculate costs of indels *)
val present_absent : a

(** [present_absent_io] contains binary data; for IO; the string is a character *)
val present_absent_io : string -> a

(** [dna] contains codes for A, G, C, and T *)
val dna : a

(** [nucleotides] contains codes for A, C, G, and T, and all their possible
    combinations, including gaps, represented as "_". The alphabet used follows
    the IUPAC specification. *)
val nucleotides : a

(** [aminoacids] contains aminoacids alphabet according to the IUPAC.
    Combinations are not considered in this alphabet. *)
val aminoacids : a 
(** [aminoacids_use_3d] is like aminoacids, but we calculate 3d costmatrix for
* it, which might take a long time.*)
val aminoacids_use_3d : a

(** [is_aminoacids] return true is alph is aminoacids, it might be with
* any level, it doesn't matter*)
val is_aminoacids : a -> bool

(** [of_string l] creates an encoding for the string list [l] to produce an 
    alphabet. *)
val of_string : ?respect_case:bool -> ?orientation:bool -> ?init3D:bool -> string list -> string -> string option -> a

(** {2 Finding} *)

(** [match_base b a] finds the code assigned in an alphabet [a] to an element
    [b]. If the element is not found raises an [Illegal_Character] exception. *)
val match_base : string -> a -> int

val gap : int

(** The obligatory representation of a gap for any alphabet *)
val gap_repr : string

(** The obligatory representation of the complement of an alphabet element. *)
val elt_complement : string

(** Same as [match_base] *)
val find_base : string -> a -> int

(** [match_base c a] finds the string representation of code [c] in the alphabet
    specification [a]. If not found raises an [Illegal_Code] exception. *)
val match_code : int -> a -> string 

(** Same as match_code *)
val find_code : int -> a -> string

(** find the code and return the expanded list of associated states **)
val find_codelist : int -> a -> int list

(** find the character associated with the polymorphism in the list **)
val find_comb : int list -> a -> int

(** Find  the specified complement of the element in the alphabet. If no
    complement is specified, return None. *)
val complement : int -> a -> int option

(** [rnd a] creates a function that generates random elements in the alphabet [a] *)
val rnd : a -> (unit -> int)

(** {2 Alphabets Properties} *)

(** [size a] gets the size of the alphabet [a]. *)
val size : a -> int

(** [get_ori_size a] returns the value of original alphabet size of alphabet
* a here *)
val get_ori_size : a -> int

(** [get_orientation a] gets the orientation of the alphabet [a]. *)
val get_orientation : a -> bool

val use_3d : a -> bool

(** [get_all a] returns the assigned code to represent all the elements in the
* alphabet [a]. *)
val get_all : a -> int option

(** [get_gap a] returns the assigned code to represent the gap element in the
* alphabet [a]. *)
val get_gap : a -> int

(** [get_missing a] return the string representing missing data **)
val get_missing : a -> string

(** [get_level a] returns the level value of alphabet a *)
val get_level: a -> int

(** [kind a] returns the kind of the alphabet [a]. *)
val kind : a -> kind

(** [distinct_size a] returns the number of distinguishable elements in [a]. Two
    elements are distinguishable if their codes are different. *)
val distinct_size : a -> int

(** [print a] debug function to print the contents of an alphabet **)
val print : a -> unit

(** [check_level a] returns true if level>1 && level<=a_sz*)
val check_level : a -> bool

(** [set_level newlevel] returns the alphabet with new level value*)
val set_level : a -> int -> a

(** {2 Extracting and Generating Alphabets} *)

(** [list_to_a l g a k] generate an alphabet using the association
    list of strings, codes, and optional complements [l], with gap
    code [g] and all code [a], to create * an alphabet of kind [k] *)
val list_to_a : ?respect_case:bool -> ?orientation:bool -> ?init3D:bool ->
  (string * int * int option) list -> string -> string option -> kind -> a

(** [simplify a] return an alphabet with the following conditions:
    If the kind of [a] is [Sequential] or [Simple_Bit_Flags], then the same 
    alphabet is returned, otherwise, only the bit flags, and the all elements are
    returned in a fresh alphabet, of type [Simple_Bit_Flags]. *)
val simplify : a -> a 

(** [to_sequential a] returns an alphabet of any kind, with its elements
    represented in the simpli:e alfied Sequential kind representation. *)
val to_sequential : a -> a

(** [create_alph_by_level alph level] creates a new alphabet based on the new
    level value *)
val create_alph_by_level : a -> int -> int -> a


(** [explote a level ori_sz] takes an alphabet of any [kind] and generates an
    [Extended_Bit_Flags] alphabet, where every combination is represented within
    square brackets. *)
val explote : a -> int -> int -> a

(** [to_list a] convert the alphabet to a list of string*code pairs *)
val to_list : a -> (string * int) list

(** [to_formatter a] produce a text-only representation of the alphabet [a].
    Useful for uniform end-user reporting. *)
val to_formatter : a -> Xml.xml

module Lexer : sig

    (** The returned list is in the inverse order of the stream. It is done like
        this because the sequences are loaded by prepending, so there is no need
        to hit the performance by keeping the order of the sequences being read.
        The boolean argument is used to issue error warnings (true) or not
        (false). *)
    val make_lexer : 
        bool -> bool -> a -> (char Stream.t -> int list -> int -> int list * int)

    (** Simmilar to [make_lexer] with two main differences: 0) first argument
        specifies the polymorphic style according to the file format i) second argument is
        a boolean for whether or not it should respect the case of the alphabet,
        and ii) it doesn't consume the entire stream but returns one element at a
        time. The list of elements is returned because a small set of alphabet
        elements can be enclosed inside [{{}}] or [{()}]. *)
    val make_simplified_lexer : [ `Nexus | `Hennig | `None ] -> 
        bool -> bool -> a -> (char Stream.t -> int list)
end

(** [of_file stream o 3d] parse an alphabet using orientation [o],
    and optionally initialize it to 3d dimensions [3d]. *)
val of_file : FileStream.f -> bool -> bool -> int -> bool ->
                a * (Cost_matrix.Two_D.m * int list list) * Cost_matrix.Three_D.m

