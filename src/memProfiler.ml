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

(** A module to profile the program's real memory usage. As we use C structures,
* we can't really trust the garbage collector, so we are going to use the vsz
* and rss reports from the standard unix ps program. *)

let do_profile = false

let current_snapshot = 
    if do_profile then
        let wd = Unix.getcwd () in
        let pid = Unix.getpid () in
        fun funct ->
            let comd = 
                Printf.sprintf "echo %s >> %s/poy_memory_profile.txt; ps -A -o vsz,rss,pid | grep %d >> %s/poy_memory_profile.txt"
                funct wd pid wd
            in
            match Sys.command comd with
            | 0 -> ()
            | _ -> 
                    Printf.printf 
                    "Error while attempting to generate memory profile\n%!"
    else fun _ -> ()

let _ = 
    if do_profile then
        Sys.command "rm -f ./poy_memory_profile.txt"
    else 0
