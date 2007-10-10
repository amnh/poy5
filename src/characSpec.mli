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

(** Character specifications in the U language *)

(** {2 Exceptions} *)
exception Illegal_Assigned_Probabilities of float
exception Undefined_Variable of string
exception Undefined_Function of string
exception Illegal_Combination
exception Precission_Failure

(** {2 Types} *)

(** The type of a function and builtin functions of U *)
type c = 
    | Base of string    (** A base *)
    | Seq of string     (** A sequence *)
    | Int of string     (** An integer *)
    | Prep              (** Prepending to a string *)
    | Hd                (** The head of a string *)
    | Tl                (** The tail of a string *)
    | Pre               (** The predeccessor of an integer *)
    | Suc               (** The successor of an integer *)

(** A character specification index *)
type t = { 
    functions: string list;     (** The list of character transformation
                                            functions *)
    variables: (string * string * c * SpecIndex.s) list; (** The list of
    variables, in each function of the character. The list holds cuadruples [(a,
    b, c, d)] where [a] is the name of the variable, [b] is the type of the
    variable, [c] is the contents of the variable in U language types, and [d]
    is the type as stored in {!SpecIndex}. *)
    counter : int All_sets.StringMap.t;  (** Counts the number of times each
    function appears in a character specification. *)
    probabilities: (string * float) list; (** Is an association list
    with the probability of each function in the character. *)
    decoders_length : float; (** Holds the calculated length of the
    functions decoder in nats. *)
}

(** The type of probability assignment *)
type p = 
    | Counter  (** The probability of a function is defined by it's count number
    in the specification. This is typical of only internal functions *)
    | Assigned of (string * float) list (** The probability of each function is
    assigned using the association list. *)

(** {2 Creation, Conversion, and Builtins} *)

(** [builtin lst] returns the list of builtin functions as specified in the U
* language. *)
val builtin : string list

(** [empty] is the empty specification *)
val empty : t

(** [to_list t] returns an association list with the functions and their count
* numbers according to the character specification [t]. *)
val to_list : t -> (string * int) list

(** [join a b c d] creates a fresh character specification that merges the
* contents of the character specifications [a] and [b], and prepends the strings
* [c] and [d] to the functions and variables in [a] and [b] to avoid name
* clashes. *)
val join : t -> t -> string -> string -> SpecIndex.t -> t

(** {2 Adding Information} *)

(** [add spec name params index] adds the function [name] with parameters and types
* [params] and using the specification index [index] to the specification
* [spec]. *)
val add : t -> string -> (string * string) list -> SpecIndex.t -> t

(** [count spec name ammount] sets the total number of times the function [name]
 * was counted in the character specification [spec] to [ammount]. *)
val count : t -> string -> int -> t

(** [estimate_prob t] updates the probabilities of the functions in [t] using
* the internal counters. *)
val estimate_prob : t -> t

(** [prob spec p] sets the probability of the functions in [spec] to the
* specification of [p]. The functions must sum less than or equal to [1.0]. If
* some functions have no assigned probability, their probability will be
* proportional to their [count] and will add up to [1]. If one of the functions
* doesn't occurr in [spec], and [Undefined_Function] exception is raised. If the
* functions probabilities don't add up to [1.0], an
* [Illegal_Assigned_Probabilities] exception is raised. *)
val prob : t -> p -> t

(** [add_decoder t f] sets the length of the decoder of the functions in [t] to
* [f]. *)
val add_decoder : t -> float -> t

(** {2 Gathering Information} *)

(** [exists a b] checks weather or not a function [b] is specified in the
* character [a] *)
val exists : t -> string -> bool

(** [find_variable spec name] collects the specification information for a
* variable [name] in [spec]. The output has the same type of
* {!CharacSpec.t.variables}. *)
val find_variable : t -> string -> string * string * c * SpecIndex.s

(** [get_string c] gets the internal string representation of a function [c]. *)
val get_string : c -> string

(** {2 Complexity} *)

(** [length spec name] calculates the encoding length in nats of a function
* [name] in the character specification [spec]. *)
val length : t -> string -> float

(** [get_decoder t] gets the length of the decoder of the functions in [t]. *)
val get_decoder : t -> float

(** [k t] calcualtes the complexity of the character specification [t] in nats.
* *)
val k : t -> float

val to_formatter : (string * t) list -> Tags.output

val non_standard : t -> string list

