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

(* $Id: charset.mli 1644 2007-02-14 19:05:47Z andres $ *)

(** Note: [Charset] implements homogeneous sets of characters in a general
    (functorized) way.  It is currently unused; sets are implemented directly in
    [Node].  Original comments follow.

    Charset.ml implements a generic set-of-characters type, as specified in
    character.mli.  Set members are identified by integer codes, i.e. chromosome
    number.  Characters are compared pairwise by element, and the median of two
    character sets will be the set containing the pairwise medians of the sets'
    elements.

    The datatype defined here is immutable.  This is easily accomplished, since
    the standard OCaml set is immutable.

    Note that this is a type of {i homogenous} sets.  All contained elements
    must be of the same class of characters.  There will be a separate module
    for heterogenous sets.

    This type implements an optimization for recalculating medians.  Whenever
    a set is made as a median, it keeps track of the sets from which it was
    made.  When this median needs to be updated due to changes in these sets,
    the old median can be passed as the first argument to [median], and the
    character set code will automatically figure out the fewest number of
    changes we need to consider.

    The problem with this optimization is that characters implicitly chain
    themselves together from most recent to oldest, so that more memory is
    technically "accessible" than will ever be accessed.  This prevents the
    garbage collector from recouping the memory.

    We will (TODO:) soon provide a functional way of removing this previous
    history to one or two ply manually, or will implement a forced "horizon" of
    one link in history.
*)

(** Character instances should always have the same codes under this
    implementation.  If this isn't true, we raise an Unmatched_code
    exception. *)
exception Unmatched_code of int

(** These describe ways of performing operations in parallel *)
type parscheme =
    | ParCodewise                       (** Send different codes to different
                                           processes *)
    | Nopar                             (** No parallelization *)

(** Make a character set out of any character.  Many of these functions come
    from the [Character.CHARACTER] signature, so we do not describe them again,
    unless there are notes about the implementation. *)
module Of :
    functor (Elt : Character.CHARACTER) ->
sig
    type e = Elt.t                      (** Type of set element *)
    type t                              (** Type of set *)

    module CharSet : Set.S
    module Elt : Character.CHARACTER

(*     type gen                            (\** Type to make random sets *\) *)
(*     val rand_gen : unit -> gen          (\** Make a random generator *\) *)
(*     val rand_gen_parallel : unit -> gen (\** Make a random generator with *)
(*                                             parallel support *\) *)
(*     val make_rand : gen -> t            (\** Make a random set *\) *)
        
    val code : t -> int
    val set_code : t -> int -> t
    val elt_code : e -> int
    val compare_codes : t -> t -> int
    val compare_data : t -> t -> int
        
    val empty : t
    val is_empty : t -> bool

    (** [derived_copy set from1 from2 from3] returns a new set with an empty
        changeset, and which claims to be derived from [set] by aligning [from1],
        [from2], and [from3].  Using [update] below, we can make changes to this set,
        one code at a time. *)
    val derived_copy : t -> t option -> t option -> t option -> t
        
    val add : e -> float -> t -> t
    val del : int -> t -> t

      (** [update code elt set1] replaces [set1]'s instance of [code] with the
          new element [elt].  If the new [elt] differs from the old one, we add
          [code] to the internal changeset; otherwise, we don't.  We use the
          cost from the previous instance of [code].  Use in conjunction with
          [fresh_copy] above. *)
    val update : int -> e -> t -> t
    val cardinal : t -> int
    val deep_cardinal : t -> int
        
    val codes : t -> int list
    val costs : t -> (int * float) list
        
    (** [merge a b] merges the characters in the given sets (by elt code) *)
    val merge : t -> t -> t

    (** [minus a b] returns the set difference [a - b] *)  
    val minus : t -> t -> t

    (** [random rboolfn set] selects elements from [set] based on the random
        boolean generator [rboolfn] *)    
    val random : (unit -> bool) -> t -> t
    val get_elt_withcode : int -> t -> e option

    (** [substitute a b]: any elements in [b] with the same codes as in [a] get
        replaced with elements from [a] *)    
    val substitute : t -> t -> t
        
    val to_list : t -> (int * e * float) list
    val of_list : (int * e * float) list -> t
        
    val set_heu : Character.h -> t -> t
    val get_heu : t -> Character.h
        
    val fold : (e -> 'a -> 'a) -> 'a -> t -> 'a

    (** [f_codes set codes] selects just the elements with codes in the list
        [codes] *)    
    val f_codes : t -> int list -> t

    (** [iter fn set] iterates over set, calling [fn elt code] for each
        element *)    
    val iter : (e -> int -> unit) -> t -> unit

    (** [map fn set] maps over set, calling [fn elt code] which should
        return a new element *)    
    val map : (e -> int -> e) -> t -> t
        
    val to_string : t -> string
        
    (** distance between two sets *)
    val distance : t -> t -> float

    (** list of (code, distance) pairs: distances between pairs with the same
        codes *)    
    val distance_list : t -> t -> (int * float) list
    
    (** median between two sets, including an optional previous median *)
    val median : t option -> t -> t -> t

    (** 3-median between four sets *)
    val median_3 : t -> t -> t -> t -> t

    (** [median_cost t] returns the cost of this median *)
    val median_cost : t -> float

    val reroot_median : t -> t -> t
    val dist_2 : t -> t -> t -> float
    val to_formatter : Tags.attributes -> t -> Data.d -> Tags.xml list
end
