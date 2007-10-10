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

(* $Id: dequeue.ml 1644 2007-02-14 19:05:47Z andres $ *)
(* Created Fri Mar 31 17:57:47 2006 (Illya Bomash) *)

(** [Dequeue] implements double-ended queues (as pairs of linked lists) *)

type 'a dequeue = 'a list * 'a list
let empty = (([], []) : 'a dequeue)
let is_empty ((a, b) : 'a dequeue) = a = [] && b = []

(** [cons q] returns a pair of the first element of [q] and [q] with that
    element removed. *)
let cons ((a, b) : 'a dequeue) =
    match a with
    | [] -> let rev = List.rev b in
      (List.hd rev, ((List.tl rev, []) : 'a dequeue))
    | elt :: elts -> (elt, (elts, b))

(** [pop_front q] is a synonym for [cons q] *)
let pop_front = cons

(** [pop_back q] returns [(e, q')], where [e] is the last element in [q] and
    [q'] is [q] with that element removed *)
let pop_back ((a, b) : 'a dequeue) =
    match b with
    | [] -> let rev = List.rev a in
      (List.hd rev, (([], List.tl rev) : 'a dequeue))
    | elt :: elts -> (elt, (a, elts))

(** [head dq] returns the head element of [dq] *)
let head dq = let head, _ = cons dq in head

(** [tail dq] returns the tail element of [dq] *)
let tail dq = let _, tail = cons dq in tail

let push_back ((a, b) : 'a dequeue) elt =
    ((a, elt :: b) : 'a dequeue)
let push_front ((a, b) : 'a dequeue) elt =
    ((elt :: a, b) : 'a dequeue)
let length ((a, b) : 'a dequeue) = List.length a + List.length b
let oflist l = ((l, []) : 'a dequeue)
let remove ((a, b) : 'a dequeue) elt =
    let filter = List.filter (fun x -> elt <> x) in
    ((filter a, filter b) : 'a dequeue)
