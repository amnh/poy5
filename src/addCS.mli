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

(** Additive Character Sets. 
*
* The additive character type has some properties that make it nicer to
* implement in group in C. For this reason a specific library for this kind of
* character was implemented in POY. The current interface follows the complete
* Set.Make specification so that it can be used as a functor parameter for set
* in scalable character sets. *)

(** {2 Exceptions} *)

exception Exists
exception Illegal_Arguments
exception Duplicated
exception Illegal_State
exception Not_Found

(** {2 Types} *)

(** A set of additive characters. Each character is defined with a minimum a
* maximum and a code. In a set, not two elements share a code, and most of the
* operations are only allowed if the sets involved share the set of codes of
* their respective characters. *)
type t

(**[min * max * code] defines the set of valid possible states in a character
* assignation, where [min] is the lower limit, [max] is the upper limit of the
* segment of states assigned and [min] <= [max], and code is the assigned code
* to the character. *)
type c = int * int * int

(** {2 Creation and Conversion } *)

(**[of_array a k] creates a new set of characters containing the characters in
 * [x] with code [k]. The function [raise Duplicated] if there is a duplicated code in the
 * array, and [raise Illegal_State] if [max] < [min] (See {!Char_addc.c}.).*)
val of_array : c array -> int -> t

(**[of_list a k] creates a new set of characters containing the characters in
* [a] with code [k].
* the function [raise Duplicated] if there is a duplicated code in the array,
* and [raise Illegal_State] if [max] < [min] (See {!Char_addc.c}.). *)
val of_list : c list -> int -> t

(**[to_list a] is symmetric to {!of_list}. *)
val to_list : t -> c list

val to_list_with_cost : t -> (int * int * int * float) list

(** [copy a b] copies the contents of the additive character set [a] to the additive 
* character set [b]. If the two character sets have a different length then
* [raise Illegal_Arguments]. *)
val copy : t -> t -> unit 

(** [clone a] creates a fresh copy of the additive character [a]. *)
val clone : t -> t 

(** {2 Distance and Medians} *)

(** [median a b] creates a set of medians between the elements in the two sets.
* The two sets must contain the same number of elements and the same codes,
* otherwise [raise Non_Equal]. *)
val median : t option -> t -> t -> t

val reroot_median : t -> t -> t
val median_3 : t -> t -> t -> t -> t 

(** [distance a b] calculates the distance between the additive character 
* set [a] and the additive character set [b]. The function raises the same
* errors as {!Char_addc.median}. *)
val distance : t -> t -> float 
val distance_2 : t -> t -> t -> float 
val dist_2 : t -> t -> t -> float 

(** [distance_median a b] is equivalent to [distance a b, median a b] but faster.
* Look in {!Char_addc.distance} and {!Char_addc.median}. *)
val distance_median : t -> t -> float * t 

(** [median_cost a] calculates the sum of distance between the two parents that gave
* rise to the additive character set [a]. If [a] belongs to an OTU then the
* function returns [0]. *)
val median_cost : t -> float 

(** [compare_codes a b] returns [0] if a and b share the same character codes in their
* sets, [-1] if for the first character with code in [a] different from [b] the
* code in [a] is less than in [b], otherwise returns [1]. This is a horrible
* description, so in other words, they are compared in lexicographic order
* according to the codes of the characters they contain. *)
val compare_codes : t -> t -> int 

(** [compare_data a b] performs a lexicographic comparison of sequences [a] and
* [b], returning [0] if [a = b], [>0] if [a > b] and otherwise [<0]. *)
val compare_data : t -> t -> int

(** {2 Particular Character Information } *)

(** [set_state a x y c] creates a fresh set of states that set the value of the
* character with code [c] with minimum [x] and maximum [y] (See {!Char_addc.c}),
* in the additive character set [a]. If the character is not found [raise
* Failed "Not_Found"]. *)
val set_state : t -> int -> int -> int -> t

(** [get_max a c] gets the maximum of the state of the character with code [c]
 * in the additive character set [a]. If the character doesn't exist in the set
 * [raise Failed "Not_Found"]. *)
val get_max : t -> int -> int 

(** [get_min a c] is symmetric to {!Char_addc.get_max}. *)
val get_min : t -> int -> int 

(** [get_cost a c] is symmetric to {!Char_addc.get_max}. *)
val get_cost : t -> int -> float

(** [find_pos a c] returns the position in the set of additive characters [a] where
* the character with code [c] is located. *)
val find_pos : t -> int -> int

(** [codes c] returns the list of codes stored in [c] *)
val codes : t -> int list

(** {2 Set.Make Interface } *)

(** [cardinal a] calculates the cardinality of [a]. *)
val cardinal : t -> int
val deep_cardinal : t -> int

(** [get_state a c] returns the current state of the character with code [c] in
 * the character set [a]. *)
val get_state : t -> int -> c

(** [elt_code c] returns the code of the additive character [c]. *)
val elt_code : c -> int

(** [set_code a] returns the code of the additive character set [a]. Note that
 * this function doesn't {i set} the value of the set code, it only returns it. *)
val set_code : t -> int
val code : t -> int

(** [elements a]. *)
val elements : t -> c list

(** [of_list l] creates a set containing the additive characters in [l]. *)
(* val of_list : c list -> t *)

(** [map f a] maps an additive character set to another additive character set.
* by applyting the [f] to each element in [a]. The resulting set has the same
* set code of [a]. *)
val map : (c -> c) -> t -> t

(** [empty k] is the empty additive character set with set code [k]. *)
val empty : int -> t

(** [is_empty a] tests wheather or not [a] is empty. *)
val is_empty : t -> bool

(** [mem i a] is true if and only if [i] is an element of [a]. It compares the
* code, minimum and maximum of each element of [a] to decide weather it is the
* same or not. *)
val mem : c -> t -> bool

(** [code_exists c a] checks weather or not the code [c] exists in the additive
* character [a]. *)
val code_exists : int -> t -> bool

(** [add c a] creates a fresh additive character set composed of the union of
* [a] and \{c\}. If the code assigned to [c] already exists in the set [raise
* Exists]. *)
val add : c -> t -> t

(** [singleton c k] creates a fresh additive character with the only element [c]
 * and set code [k].*)
val singleton : c -> int -> t

(** [remove c a] creates a fresh additive character set with all the elements in
 * [a] excepting [c]. The function compares the code, minimum and maximum of the
 * state of the characters to perform the remotion. *)
val remove : c -> t -> t

(** [union a b] creates a fresh additive character set containing the elements
* of [a] and [b] and set code of b. If a particular code is contained in both [a] 
* and [b] the state of the corresponding code in the output set is the one found 
* in [a]. *)
val union : t -> t -> t  

(** [inter a b k] creates a fresh additive character set containing the
* intersectoin of [a] and [b] and set code [k]. Uses {!Char_addc.mem} to perform the
* comparisons. *)
val inter : t -> t -> int -> t

(** [diff a b] creates a fresh additive character set containing the difference
* between [a] and [b] (elements in [a] not present in [b]) and set code [k]. Uses
* {!Char_addc.mem} to perform the comparisons. *)
val diff : t -> t -> int -> t

(** [subset a b] is true if every element in [a] is in [b]. Uses the
* {!Char_addc.mem} to perform the comparisons. *)
val subset : t -> t -> bool

(** [equal a b] is true if every element in [a] is in [b] and every element in
 * [b] is in [a]. Uses the {!Char_addc.mem} to perform the comparisons. *)
val equal : t -> t -> bool

(** [iter f a] applies the function [f] to every additive character in [a]. *)
val iter : (c -> unit) -> t -> unit

(** [fold f a x] is analogous to List.fold_left, where the list of additive
* characters is sorted by code. *)
val fold : ('a -> c -> 'a) -> t -> 'a -> 'a

(** [for_all f a] is true iff [f x] is true for every element [x] in [a]. *)
val for_all : (c -> bool) -> t -> bool

(** [exists f a] is true iff [f x] is true for some element [x] in [a]. *)
val exists : (c -> bool) -> t -> bool

(** [filter f a] creates a fresh set contining every element [x] in [a] that
* evaluates true for [f x]. *)
val filter : (c -> bool) -> t -> t

val f_codes : t -> All_sets.Integers.t -> t
val f_codes_comp : t -> All_sets.Integers.t -> t

(** [partition f a] is equivalent (but faster) to 
* [filter f a (set_code a), filter (not f) a (set_code a)].
* *)
val partition : (c -> bool) -> t -> t * t

(** [min_elt a] returns the character with the minimum code in [a]. *)
val min_elt : t -> c

(** [max_elt a] returns the character with the maximum code in [a]. *)
val max_elt : t -> c

(** [choose a] selects one element from [a]. It guarantees that for two sets with
* the same contents, both choose the same element. *)
val choose : t -> c

(** [full_union a b c] stores in the character [c] the full union of the
* character sets [a] and [b]. The full union of two additive characters [(x_min,
* x_max)] and [(y_min, y_max)] is [((min x_min y_min), (max x_max y_max))], and
* this is performed for every pair of corresponding characters in [a] and [b].
* *)
val full_union : t -> t -> t -> unit 

(** [split c a] creates a triple [(x, y, z)] where [x] is the set of additive
* characters with code less than [c], [y] is true iff element [c] exists in [a]
* and [z] containes those elements with code greater than [c]. [x] and [z] 
* have the same set code of [a].*)
val split : c -> t -> t * bool * t

(** [of_parser p] generates a character set from the output of the PoyParser.
* For more information look at its usage in the {!Node} module. *)
val of_parser : Data.d -> ((int list option * int) array * int) -> int -> t * int
val is_potentially_informative : int list option list -> bool
val min_possible_cost : int list option list -> float

(** [to_string a] generates a string representation of the additive character
* set [a] *)
val to_string : t -> string 

(** [state_to_xml ch x d] generates a XML representation of the character set
* [x] using the data of [d], with the output dumped in the channel [ch]. *)
val state_to_xml : 
    Pervasives.out_channel -> t -> Data.d -> unit

(** [to_formatter attrs c parent d : Tags.output list] returns the formatter for
    node c where parent is optional parent of c if available *)
val to_formatter : Tags.attributes -> t -> t option -> Data.d -> Tags.output list

(** {2 Convenient Modules} *)

module Imperative : sig
    (** An imperative implementation of the Additive Character Set 
    *
    * This is a mirror of the functional definition of an additive character set
    * in an imperative form for efficiency purposes, if the application shows
    * the necessity. See {!AddCS} for the documentation of the respective
    * funcitons. *)
    type it 
    type ic
    val create : t -> it
    val of_array : ic array -> int -> it
    val of_list : ic list -> int -> it
    val to_list : it -> ic list
    val copy : it -> it -> unit
    val clone : it -> it
    val median : it -> it -> it -> unit
    val distance : it -> it -> float 
    val median_cost : it -> float
    val compare : it -> it -> int
    val set_state : it -> ic -> unit
    val get_max : it -> int -> int
    val get_min : it -> int -> int
    val get_cost : it -> int -> float 
    val cardinal : it -> int
    val get_state : it -> int -> ic
    val elt_code : ic -> int
    val get_set_code : it -> int
end

module Test : sig
    (** An Additive Character Set test module. 
    *
    * A module that can be turned on with a debug flag in the source code to
    * self check the consistency of the Additive Character Set library.
    * It also provides functionality for self tests at runtime.
    *)

    (** [random a b c d] creates a fresh set of additive characters with [a]
     * characters, with states between [b] and [c] and set code [d]. *)
    val random : int -> int -> int -> int -> t

    (** [run t] runs a test [t] which is a tuple of generating function and the
    * output check. All the following [test_*] functions generate valid output
    * for [run] *)
    val run : ((unit -> 'a) * ('a -> 'b)) -> 'b

    val test_full_union : int -> int -> int ->
        (unit -> (t * t * t)) * ((t * t * t) -> (t * t * t))

    val test_of_array : int -> int -> int -> 
        (unit -> c array * t) * ((c array * t) -> (c array * t))

    val test_get_functions : int -> int -> int ->
        (unit -> c array * t) * ((c array * t) -> (c array * t))
end
