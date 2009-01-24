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

(** This module implements sets of Sankoff characters *)

(** A cost matrix is a matrix of integers.  Costs cannot be infinite. *)
type cm = int array array

type t                                  (** The type of the set *)
type elt                                (** The type of set elements *)

(** {2 Standard functions} *)
val color : Character.c
val code : t -> int                     (** The code of a set *)
val to_list : t -> (int * int array) list
val ecode : elt -> int                  (** The code of an element *)
val codes : t -> int list               (** Set of ecodes stored here *)
val tcm : t -> cm                       (** Get the cost matrix of a set *)
val set_tcm : cm -> t -> t              (** Set the cost matrix *)
val median : t option -> t -> t -> t    (** See {!Character.CharacterSet} *)
val distance : t -> t -> float          (** See {!Character.CharacterSet} *)
val reroot : t -> t -> t -> t
    (** Sankoff has its own clade-rerooting procedure *)
val compare_codes : t -> t -> int       (** See {!Character.CharacterSet} *)
val compare_data : t -> t -> int        (** See {!Character.CharacterSet} *)
val parse : 'a -> (string * t) list     (** See {!Character.CharacterSet} *)
val to_string : t -> string             (** See {!Character.CharacterSet} *)
val median_3 : t -> t -> t -> t -> t    (** See {!Character.CharacterSet} *)
val dist_2 : t -> t -> t -> float       (** See {!Character.CharacterSet} *)
val cardinal : t -> int

(** {2 Creation} *)
val make_random_tcm : ?max:int -> int -> cm
    (** makes a random cost matrix with given dimension *)
val make_sank : int -> cm -> (int * int list) list -> t
    (** [make_sank code tcm eltlist] makes a new Sankoff set; the list of
        elements should be made up of [(code, states)] pairs, giving the code of
        the element and the list of observed states. *)

(** {2 Output and other functions} *)
val of_parser : cm -> ((Parser.SC.static_state * int) array * int) -> int -> t * int

val print_tcm : cm -> unit              (** prints to standard output *)

val get_minstates : elt -> int list
    (** get the list of optimal states at the node *)
val get_elt_array : t -> elt array
    (** return the set as an array of elements *)
val f_codes : t -> All_sets.Integers.t -> t
    (** filter the Sankoff set given a list of codes *)
val f_codes_comp : t -> All_sets.Integers.t -> t
    (** filter the Sankoff set given a list of codes to exclude *)

(** [to_formatter attrs c parent d : Xml.xml list] returns the formatter for
    node c where parent is optional parent of c if available *)
val to_formatter : Xml.attributes -> t -> t option -> Data.d -> Xml.xml
Sexpr.t list

val get_all_possible_assignments : int list option list -> int list list 

val min_possible_cost : int array array -> Parser.SC.static_state list -> float

val max_possible_cost : int array array -> Parser.SC.static_state list -> float
