(* POY 5.0. A phylogenetic analysis program using Dynamic Homologies.         *)
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

(** This module implements sets of equally-weighted non-additive characters
    in C. These sets are immutable but can share data through reference
    counting on the C side. C is used mainly to ensure fast median-taking;
    medians are made with vectorizable bitwise operators.

    This implementation can handle sets of size up to the word size in bits,
    i.e., 32 characters on 32-bit machines, 64 characters on 64-bit machines.
    Using n-1 characters (31, 63) ensures complete interoperability with OCaml;
    using 32 or 64 allows you to use all the standard functions, but some of the
    raw operations included will fail.  However, these should only be used for
    debugging in the first place.

    They are implemented as custom blocks in C, allowing for serialization,
    deserialization, and comparison.  The functions listed here mainly come from
    the interfaces mentioned above; any additional functions are explained here. *)


(** {2 Types} *)

(** Abstract type that defines a set of NonAdditive characters *)
type ct

(** The Abstract type for union based non additive characters *)
type cu

(** The type of a set of nonadditive characters *)
type t = {
    codes : int array;
    data : ct;
}

 (** The type of individual set elements *)
type e = int * int


(** {2 Accessing Codes and Querying Codes} *)

val mem : int list option -> t -> bool

val code : ct -> int



(** {2 Union operations} *)

val union : ct -> cu -> cu -> cu

val to_union : ct -> cu

val elt_code : e -> int


(** {2 Creation} *)

val make_new : int -> int -> ct
(** Makes a new nonadditive character set of a given size and code.  All
    elements are initialized to zero, which is an inconsistent state.  (Please
    see this module's TeX file for details.) *)

val set_elt : ct -> int -> e -> unit
(** [set_elt set loc val] sets element number [loc] to value [val].  [val] is
    treated as bit-encoding its possible states. *)

val set_elt_bit : ct -> int -> int -> unit
(** [set_elt_bit set loc bit] sets bit [bit] of element [loc], signifying that
    element [loc] is / can be in state [bit]. *)


(** {2 Medians and distances} *)

val median : 'a -> t -> t -> t
(** [median _ child1 child2] determine the median of two non-additive characters *)

val median_3 : t -> t -> t -> t -> t
(** [median_3 parent node child1 child2] returns node_final *)

val distance : ct -> ct -> float

val distance_list : ct -> ct -> (int * float) list

val dist_2 : t -> t -> t -> int
(** [dist_2 a b c] returns the smallest distance from [a] to either [b] or [c],
    computed element-wise *)

val reroot_median : ct -> ct -> ct
(** [reroot_median a b] performs a special median for finding the root value
    between two nodes for which the final states are known.  This is done by
    taking the union of the states, as described in Goloboff 1993. *)


(** {2 State, listing, parsing} *)

val cardinal : ct -> int

val cardinal_union : cu -> int
(** number of elements in set *)

val poly_items : cu -> int -> int

val to_int : ct -> int -> int

(** [elt_to_list set eltnum] returns a list of states that element [eltnum]
    might be in. ~Deprecated *)
val elt_to_list : ct -> int -> int list

(* this is the prefered method to inspect set bits. *)
val e_to_list : e -> int list

val to_list : ct -> (int * e * float) list

val of_list : (int * e * float) list -> t

val of_parser : Data.d -> (Parser.SC.static_state * int) array * 'a -> int -> t * 'a

val is_potentially_informative : int list option list -> bool

val max_possible_cost : int list option list -> float

val min_possible_cost : int list option list -> float

val to_simple_list : t -> (int * e) list

val to_string : t -> string
(** [to_formatter attrs c parent d : Xml.xml list] returns the formatter for
    node c where parent is optional parent of c if available *)
val to_formatter : Xml.xml Sexpr.t list -> 
    Xml.attributes -> t -> t option -> Data.d -> 
        Xml.xml Sexpr.t list


(** {2 Other standard functions} *)

val compare_codes : t -> t -> int
val compare_data : t -> t -> int
val empty : t
val add : e -> float -> t -> t
val del : int -> t -> t
val codes : t -> int list
val costs : t -> (int * float) list
val get_elt_withcode : int -> t -> e option
val pair_map : ('a -> 'b) -> 'a * 'a -> 'b * 'b
val substitute : t -> t -> t
val merge : t -> t -> t
val minus : t -> t -> t
val random : (unit -> bool) -> t -> t
val fold : (e -> 'a -> 'a) -> 'a -> t -> 'a
val filter : (int * e * float -> bool) -> t -> t
val f_codes : t -> All_sets.Integers.t -> t
val f_codes_comp : t -> All_sets.Integers.t -> t
val iter : (e -> int -> unit) -> t -> unit
val map : (e -> int -> e) -> t -> t
val is_empty : t -> bool
