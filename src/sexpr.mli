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

type 'a t = [ `Empty | `Set of 'a t list | `Single of 'a ]

val fold_left :
  ('a -> 'b -> 'a) -> 'a -> 'b t -> 'a

val fold_right :
  ('a -> 'b -> 'a) -> 'a -> 'b t -> 'a

val rev: 'a t -> 'a t

val map :
  ('a -> 'b) -> 'a t -> 'b t

val map_fold_left :
  ('a -> 'b -> 'a * 'c) ->
  'a -> ([< `Empty | `Set of 'd list | `Single of 'b ] as 'd) -> 'c t * 'a

val cardinal : 'a t -> int

val map_feedback :
  (int -> unit) -> ('a -> 'b) -> 'a t -> 'b t

val map_status :
    string -> ?eta:bool -> ('a -> 'b) -> 'a t -> 'b t

val fold_status :
    string -> ?eta:bool -> ('a -> 'b -> 'a) -> 'a -> 'b t -> 'a

val compose_status :
    string -> ?eta:bool -> (int -> 'a -> 'a) -> int -> 'a -> 'a

val singleton : 'a -> 'a t

val of_list : 'a list -> 'a t

val to_list : 'a t -> 'a list

val first : 'a t -> 'a

val all_to_all :
  ('a -> 'b -> 'c) -> 'a t -> 'b t -> 'c t

val flatten : 'a t t -> 'a t

val iter : ('a t -> unit) -> 'a t -> unit

val shallow_all_to_all : ('a t -> 'a t -> 'b) -> 'a t -> 'a t -> 'b t

val leaf_iter : ('a -> unit) -> 'a t -> unit

val to_array : 'a t -> 'a array t

val length : 'a t -> int

val nth : int -> 'a t -> 'a

val filter : ('b -> bool) -> 'b t -> 'b t

val split : ('b -> bool) -> 'b t -> ('b t * 'b t)

val combine : 'a t * 'b t -> ('a * 'b t) t

val full_combine : 'a t * 'b t -> ('a * 'b) t

val map_insexpr : ('a -> ('b t) list) -> 'a t -> 'b t

val choose_random : 'a t -> 'a option

val union : 'a t -> 'a t -> 'a t
