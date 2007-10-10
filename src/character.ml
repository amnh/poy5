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

(* $Id: character.ml 1644 2007-02-14 19:05:47Z andres $ *)
let () = SadmanOutput.register "Character" "$Revision: 1644 $"


(* To later interface with Vamsi's code *)
type update_res = Continue | Stop
let rContinue = Continue
let rStop = Stop

let autocode = ref 5
let new_code () =
    incr autocode;
    !autocode

type t = Indels

type c = White | Black | Blue | Yellow | Red | Grey | Green | Brown

type h = Fast | Medium | Exact

(* Why functors and not object signatures:
 *
 * We have many more static operations than object-embedded operations.
 * Furthermore, our predominant data types will contain homogenous objects which
 * should interact with each other but not with objects of other types.  Thus, I
 * think polymorphism should happen on the level of code bases (modules,
 * functors) rather than as a class hierarchy.
 *)

module type CHARACTER = sig
    type t
        
(*     val color : c *)
    val code : t -> int

(*     type gen *)
(*     val rand_gen : unit -> gen *)
(*     val make_rand : gen -> t *)

    (* basic median-making functions *)
    val median : t option -> t -> t -> t
    val median_3 : t -> t -> t -> t -> t
    val reroot_median : t -> t -> t
    val dist_2 : t -> t -> t -> float

    val distance : t -> t -> float

    val compare_codes : t -> t -> int
    val compare_data : t -> t -> int

    val cardinal : t -> int
    val deep_cardinal : t -> int

    (* input and output *)
    val to_string : t -> string
end

module type CharacterSet = sig
    type e                              (* contents of the set *)
    type t                              (* set type *)

(*     type gen *)
(*     val rand_gen : unit -> gen *)
(*     val make_rand : gen -> t *)

(*     val stype : t *)
    val empty : t

    val cardinal : t -> int
    val deep_cardinal : t -> int

(*     val color : c *)
    val code : t -> int
    val set_code : t -> int -> t
    val elt_code : e -> int

    val add : e -> float -> t -> t
    val del : int -> t -> t

(*     val colors : t -> (int * c) list *)
    val codes : t -> int list
    val costs : t -> (int * float) list
    val get_elt_withcode : int -> t -> e option
        
    val substitute : t -> t -> t
    val merge : t -> t -> t
    val minus : t -> t -> t
    val random : (unit -> bool) -> t -> t

    val to_list : t -> (int * e * float) list
    val of_list : (int * e * float) list -> t

    val set_heu : h -> t -> t
    val get_heu : t -> h

    val fold : (e -> 'b -> 'b) -> 'b -> t -> 'b
    val f_codes : t -> int list -> t
(*     val f_colors : t -> c -> t *)
    val iter : (e -> int -> unit) -> t -> unit
    val map : (e -> int -> e) -> t -> t

    val is_empty : t -> bool
    val median : t option -> t -> t -> t
(*     val median_3 : t -> t -> t -> t *)
    val to_string : t -> string
    val distance_list : t -> t -> (int * float) list
    val distance : t -> t -> float
    val compare_codes : t -> t -> int
    val compare_data : t -> t -> int

(*     val update : t option -> t -> t -> t -> (t * update_res) *)
end

