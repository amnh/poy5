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

(** This type interacts with the tree traversal code to continue or stop traversal *)
type update_res
val rContinue : update_res
val rStop : update_res

(** Utility function to create a new, unique code for a character module *)
val new_code : unit -> int

(** Colors for characters *)
type c = White | Black | Blue | Yellow | Red | Grey | Green | Brown
(** Heuristic values *)
type h = Fast | Medium | Exact

(** The [CHARACTER] signature is the general signature for different types of
    characters.  This signature is currently missing several functions:
    * [median_3]
    * [update] for performing a 3median and returning an [update_res] value
    Note that characters are meant to be immutable values. *)
module type CHARACTER =
  sig
    type t

    (** Every character has a code; characters with the same code are
        homologous. *)
    val code : t -> int

    (** Type for generating random instances of this character.  Used for
        testing. *)
(*     type gen *)

    (** Make a random set of parameters for generating mutually consistent and
        comparable characters *)
(*     val rand_gen : unit -> gen *)

    (** Given a set of parameters, make a random character *)
(*     val make_rand : gen -> t *)

    (** [median prev_median a b] returns the median of [a] and [b].  Optionally,
        you may specify [Some prev] as [prev_median], the median of the previous
        versions of [a] and [b].  This is useful if [a] and [b] themselves keep
        a list of changes from their previous incarnations.  See the module
        [Charset] for an example of how this is used. *)
    val median : t option -> t -> t -> t

    val median_3 : t -> t -> t -> t -> t
    val reroot_median : t -> t -> t
    val dist_2 : t -> t -> t -> float

    (** Distance between two character instances. *)
    val distance : t -> t -> float

    (** [compare_codes a b] imposes a total ordering on the set of characters
        based on their homologies.  For homologous characters (those with the
        same codes), this function should return 0. *)
    val compare_codes : t -> t -> int

    (** [compare_data a b] imposes a total ordering on the set of characters
        based on their codes and data.  Only in cases where both the codes and
        the data are the same will this function return 0. *)    
    val compare_data : t -> t -> int

    val cardinal : t -> int
    val deep_cardinal : t -> int

    (** Returns a string representation of a character instance. *)
    val to_string : t -> string

  end

(** [CharacterSet]s are [CHARACTER]s in their own right, but also have
    additional functionality relevant to sets of characters.  The [Charset]
    module defines a functor that implements sets of a generic character. *)
module type CharacterSet = sig
    type e                              (** contents of the set *)
    type t                              (** set type *)

(*     (\** Type for generating random instances of this character.  Used for *)
(*         testing. *\) *)
(*     type gen *)
(*     (\** Make a random set of parameters for generating mutually consistent and *)
(*         comparable characters *\) *)
(*     val rand_gen : unit -> gen *)
(*     (\** Given a set of parameters, make a random character *\) *)
(*     val make_rand : gen -> t *)

    (** An empty set of characters *)
    val empty : t
    (** Find the cardinality (number of elements) of a set *)
    val cardinal : t -> int
    val deep_cardinal : t -> int

    val code : t -> int
    val set_code : t -> int -> t
    val elt_code : e -> int

    (** [add code elt cost set] adds an element [elt] with a given code and cost
        to the character set [set] *)
    val add : e -> float -> t -> t
    (** [del code set] returns a new set with [code] deleted *)
    val del : int -> t -> t

    (** return the list of codes *)
    val codes : t -> int list
    (** return the list of code/cost pairs *)
    val costs : t -> (int * float) list
    (** [get_elt_withcode code set] returns the element with code [code], if it exists *)
    val get_elt_withcode : int -> t -> e option
        
    (** [substitute a b] takes any element in [a] with a corresponding element
        in [b] and replaces the element in [b] with the one in [a].  The set of
        codes in [b] is not changed, but elements in [b] are overwritten by ones
        from [a].  Of course, the sets are immutable, so the resulting set is
        returned. *)
    val substitute : t -> t -> t
    (** [merge a b] returns the merge of [a] and [b] *)    
    val merge : t -> t -> t
    (** set subtraction *)
    val minus : t -> t -> t
    (** [random p set] returns a random subset of [set], given a random boolean
        generator [p]. *)
    val random : (unit -> bool) -> t -> t

    (** Turns a set into a [(code, elt, cost)] list *)    
    val to_list : t -> (int * e * float) list
    (** Takes a list of [(code, elt, cost)] elements and returns a set *)
    val of_list : (int * e * float) list -> t

    (** gets the heuristic currently used on a given set *)    
    val set_heu : h -> t -> t
    (** sets the heuristic to be used *)    
    val get_heu : t -> h

    (** [fold fun start set] calls [fun elt accum] to fold the set *)    
    val fold : (e -> 'b -> 'b) -> 'b -> t -> 'b
    (** [f_codes set list] returns a subset with only codes that appear in list
    *)    
    val f_codes : t -> int list -> t

    val iter : (e -> int -> unit) -> t -> unit
    (** [iter fun set] calls [fun elt code color] for each element *)    

    (** [map fun set] returns a new set, made by calling [fun elt code color],
        which should return a new element *)    
    val map : (e -> int -> e) -> t -> t

    val is_empty : t -> bool
        
    (** [median prev left right] should calculate a median of [left] and
        [right].  [prev] is provided as a set of medians of the "previous"
        versions of [left] and [right].  However, under certain tree operations,
        [prev] will be the median of other nodes;  the program must check for
        this case.  See charset.ml for an example. *)
    val median : t option -> t -> t -> t

    val to_string : t -> string
    val distance_list : t -> t -> (int * float) list
    val distance : t -> t -> float
    val compare_codes : t -> t -> int
    val compare_data : t -> t -> int


(*     val update : t option -> t -> t -> t -> (t * update_res) *)
(*     val median_3 : t -> t -> t -> t *)
end

