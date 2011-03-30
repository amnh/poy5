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

(** Operations to Fingerprint (generate Hash values) using the Karp-Rabi
* technique *)

(** {2 Types} *)

exception Illegal_argument
exception Not_match

(** {2 Functors} *)

module type SERIALIZABLE = sig
    (** A serializable type for fingerprinting using the Karp-Rabi techinque *)

    (** The type of the serialized type *)
    type ser

    (** The type of each element in the serialized type *)
    type v

    (** Compares two elements of the serialized type *)
    val equal : v -> v -> bool

    (** [get_pos ser p] Gets the contents in position p of the serialized ser. If the
    * position is less than 0 or greated than the length of ser a Failure is
    * raised. *)
    val get_pos : ser -> int -> v

    (** [to_int v] gets an integer representation of the value v. *)
    val to_int : v -> int

    (** [of_int v] converts an integer to the internal representation of the
    * value v. *)
    val of_int : int -> v

    (** [length ser] gets the length of the serialized value ser *)
    val length : ser -> int

    (** [sub ser t l] gets a slice of the serialized ser starting in position t with
    * length l. If ser is not long enough to get the slice, an Illegal_argument
    * error is raised. *)
    val sub : ser -> int -> int -> ser

end

module Make : functor (Elt : SERIALIZABLE) -> sig
    (** Fingerprinting functor using the Karp-Rabi algorithm for serializable 
    * types *)

    (** [fp t q i n a] gets the fingerprint of a subsection of the serialized
    * t using module q (a prime number) for the fingerprinting, starting in
    * position i of t, for a slice of length n using an alphabet size a. *)
    val fp : Elt.ser -> int -> int -> int -> int -> int

    (** [complete_fp t q a] calculates the fingerprint of the complete
    * serialized t using module q and an alphabet size a *)
    val complete_fp : Elt.ser -> int -> int -> int

    (** [fp_incremental t q i n cur a] performs an incremental recalculation of
    * a fingerprint of the slice starting in position i + 1 of the serialized
    * t with length n, using the cur fingerprint of the slice of length n
    * starting in position i of t. Uses q as the module for the calculations
    * (should be the same module used to calculate q), and uses the alphabet
    * size a. *)
    val fp_incremental : Elt.ser -> int -> int -> int -> int -> int -> int

    (* [fp_match t p a q] finds the first match of the pattern p in the text t
    * using alphabet size a and calculations module q. *)
    val fp_match : Elt.ser -> Elt.ser -> int -> int -> int

    val do_match : Elt.ser -> Elt.ser -> int -> int

end

module MakeMult : functor (Elt : SERIALIZABLE) -> sig
    (** Fingerprinting functor using the Karp-Rabi algorithm for serializable 
    * types, using lists of values as hash values as opposed to single values as
    * in {!Make}. 
    *
    * The functions specifications are symmetric to those in {!Make}. *)

    val fp : Elt.ser -> int list -> int -> int -> int -> int list

    val complete_fp : Elt.ser -> int list -> int -> int -> int list

    val fp_incremental : 
        Elt.ser -> int list -> int -> int -> int list -> int -> int list

end

(** {2 Concrete Modules} *)

(** Some useful instances of the {!Make} and {!MakeMult} *)

module SerGenArray : sig
    type v = int
    type ser = v array
    val equal : v -> v -> bool
    val get_pos : ser -> int -> v
    val to_int : v -> int
    val of_int : int -> v
    val length : ser -> int
    val sub : ser -> int -> int -> ser
end

module IntArrayFp : sig
    val fp : SerGenArray.ser -> int -> int -> int -> int -> int
    val complete_fp : SerGenArray.ser -> int -> int -> int
    val fp_incremental :
        SerGenArray.ser -> int -> int -> int -> int -> int -> int
    val fp_match : SerGenArray.ser -> SerGenArray.ser -> int -> int -> int
    val do_match : SerGenArray.ser -> SerGenArray.ser -> int -> int
end

module SerString : sig
    type v = char
    type ser = string
    val equal : v -> v -> bool
    val get_pos : ser -> int -> v
    val to_int : v -> int
    val of_int : int -> v
    val length : ser -> int
    val sub : ser -> int -> int -> ser
end

module StringFp : sig
    val fp : SerString.ser -> int -> int -> int -> int -> int
    val complete_fp : SerString.ser -> int -> int -> int
    val fp_incremental :
        SerString.ser -> int -> int -> int -> int -> int -> int
    val fp_match : SerString.ser -> SerString.ser -> int -> int -> int
    val do_match : SerString.ser -> SerString.ser -> int -> int
end
