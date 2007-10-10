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

(** Convenient extensions of the Array module in the Standard library. *)

(* $Id: array_ops.mli 1644 2007-02-14 19:05:47Z andres $ *)

(** [flatten_array a] converts an array array ref referenced by a into a
 single array of the same type of a, pointed by tgt.*)
val flatten_array : 'a array array -> 'a array 

(** [fold_min x y z] applies to each element i of z the functions xy and x
(that is the composition of x and y, and x itself). The x i of the minima
of xy i are returned in a fresh array. *)
val fold_min : ('a -> 'b) -> ('b -> int) -> 'a array -> 'b array 

(** Same as fold_min but returns the set that provides the maximum result *)
val fold_max : ('a -> 'b) -> ('b -> int) -> 'a array -> 'b array 

(** [map_ip x y] performs an inplace mapping, that is, applies the function x to
* each element in y and replace it *)
val map_ip : ('a -> 'a) -> 'a array -> unit

(** [mapi_ip x y] performs an inplace mapping, that is, applies the function x to
* each element in y and it's position in y *)
val mapi_ip : (int -> 'a -> 'a) -> 'a array -> unit

(** [randomize x] randomize the elements in the array x in place *)
val randomize : 'a array -> unit

(** [filter f x] creates a fresh array with the elements that are true for [f]
 * in [x]. *)
val filter : ('a -> bool) -> 'a array -> 'a array

(* vim:sw=4 et tw=80 : *)
