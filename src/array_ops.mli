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

(** Convenient extensions of the Array module in the Standard library. *)

(* $Id: array_ops.mli 2871 2008-05-23 17:48:34Z andres $ *)

(** [array_append_caml3 arr1 arr2] return Array.append arr1 arr2.
* this is a walk around of seg fault problem causing by Array.append of Ocaml4.0.X. 
* here we create our version of Array.append, does the job just like append function from Ocaml3.12.X. *)
val array_append : 'a array -> 'a array -> 'a array


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

val mem : 'a array -> 'a -> bool

val map_2 : ('a -> 'b -> 'c) -> 'a array -> 'b array -> 'c array

val map_3 : ('a -> 'b -> 'c -> 'd) -> 'a array -> 'b array -> 
    'c array -> 'd array

val map_4 : ('a -> 'b -> 'c -> 'd -> 'e) -> 'a array -> 'b array -> 
    'c array -> 'd array -> 'e array

val map_5 : ('a -> 'b -> 'c -> 'd -> 'e -> 'f) -> 'a array -> 'b array -> 
    'c array -> 'd array -> 'e array -> 'f array

val fold_right_2 : ('a -> 'b -> 'c -> 'a) -> 'a -> 'b array -> 'c array -> 'a
val fold_right_3 : ('a -> 'b -> 'c -> 'd -> 'a) -> 'a -> 'b array -> 'c array ->
    'd array -> 'a
val fold_right_4 : ('a -> 'b -> 'c -> 'd -> 'e -> 'a) -> 'a -> 'b array -> 'c array ->
    'd array -> 'e array -> 'a

val fold_righti : (int -> 'a -> 'b -> 'b) -> 'b -> 'a array -> 'b

val split : int -> 'a array -> 'a array array

(** [is_identical2 arr1 arr2] returns 0 if arr1 and arr2 are not the same, 1 if
    * they are the same.*)
val is_identical2 : int array -> int array -> int

val is_identical3 : int array -> int array -> int array -> int

val fill_symmetric_square_matrix : ?status : string option ->
    (int -> int -> 'a -> 'a -> float) -> 'a array -> float array array -> unit

