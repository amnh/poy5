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

let () = SadmanOutput.register "Profiles" "$Revision: 1644 $"

(* Book keeping of the current state of the processes *)
type condition =
    | Busy
    | Idle;;

let _notags = -1
and _nochekpt = ""
and _waiting_parameter = 1;;

let avail_proc = ref Queue.create ()
and true_ranks = ref [||]
and all_proc = ref [||];;

let rec init_proc c it =
    match c with 
    | [] -> ();
    | h :: t ->
            Queue.push h !avail_proc;
            begin match it with
            | 0 ->
                    all_proc := Array.create (List.length c) h;
                    true_ranks := Array.create (List.length c) 0;
            | _ ->
                    let rank, comm = h in
                    !all_proc.(it) <- (rank, comm, Idle, Alive, _notags,
                    _nochekpt);
                    !true_ranks.(it) <- c;
            end;
            init_proc t (it + 1);;

let get_position r = 
    !true_ranks.(r);;

let rec get_proc tag chekpt =
    match Queue.is_empty !avail_proc with
    | false ->
            let true_rank = get_position p in
            let rank, comm, _, _, _, _ = !all_proc.(true_rank) in
            !all_proc.(true_rank) <- rank, comm, Idle, Alive, _notags,
            _nochekpt;
            Queue.push (p, comm);;

let clean_chekpt p =
    let true_rank = get_position p in
    let r, c, s, l, t, _ = !all_proc.(true_rank) in
    !all_proc.(true_rank) <- r, c, s, l, t, _nochekpt;;

let clean_tag p =
    let true_rank = get_position p in
    let r, c, s, l, _, cp = !all_proc.(true_rank) in
    !all_proc.(true_rank) <- r, c, s, l, _notags, cp;;

let clean_all p =
    clean_tag p;
    clean_chekpt p;;
