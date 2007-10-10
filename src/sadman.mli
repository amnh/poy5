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

(* $Id: sadman.mli 1800 2007-05-08 01:02:29Z andres $ *)

(** The [Sadman] module is an OCaml interface to the Sadman status printing
    system (which is currently under development).  It facilitates the logging
    of status information whenever an operation is performed.  Such "annotated"
    functions return key-value pairs called "features" that describe their input
    (both the identities of their arguments and whatever higher-level features
    of the input might be interesting) and the algorithm used,  and key-value
    pairs called "results" that describe properties of their results.

    The goal is to store these values in a database searchable by features and
    results.  If the descriptions are sufficiently good and sufficiently
    automated, then these results can be stored in the database automatically,
    which provides a convenient logging feature for the program's
    development. *)

(** A property has a name and a value *)
type prop = (string * string)
(*

(** [print_features ch prop_list] prints an XML version of the property list on
    the given channel.  It prints the properties as "features" of that
    particular run of the program. *)
val print_features : out_channel -> prop list -> unit

(** [print_results ch prop_list] prints an XML version of the property list on
    the given channel.  It prints the properties as results of that
    particular run of the program. *)
val print_results : out_channel -> prop list -> unit

*)
(** [prefix str prop_list] affixes a prefix to all the names of the properties
    in the list.  They are separated by periods. *)
val prefix : string -> prop list -> prop list

val start : string -> prop list -> unit
val finish : prop list -> unit
