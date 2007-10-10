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

(** Asynchronous non fault tolerant operations
*
* Provides parallel operations in asynchronous environment using group
* communications (the main difference with asynch.ml). It assumes that every
* process is running the same code (no slave loop as in asynch.ml. 
* The interface is the same as in {!Paroper} and {!Asynch}. *)

(* $Id: paroperas.mli 1644 2007-02-14 19:05:47Z andres $ *)

val all_minimize : ('a -> 'b) -> ('b -> float) -> 'a array ref -> float
val all_maximize : ('a -> 'b) -> ('b -> float) -> 'a array ref -> float
val gather_minimum :
  ('a -> 'b) -> ('b -> float) -> 'a array ref -> float -> 'b array ref 
val gather_maximum :
  ('a -> 'b) -> ('b -> float) -> 'a array ref -> float -> 'b array ref
val map : ('a -> 'b) -> 'a array ref -> 'b array ref
val scan : ('a -> bool) -> 'a array ref -> 'a array ref

(* vim:sw=4 et tw=80
 *)
