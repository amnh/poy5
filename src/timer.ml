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

IFDEF USEWIN32 THEN
type rt = {
    utime : float;
    stime : float;
}

let c_start _ =
    let time = Unix.gettimeofday () in
    { utime = time; stime = time }

ELSE
type rt = {
    utime : float;
    stime : float;
    maxrss : int;
    ixrss : int;
    idrss : int;
    isrss : int;
    minflt : int;
    majflt : int;
    nswap : int;
    inblock : int;
    oublock : int;
    msgsnd : int;
    msgrcv : int;
    nsignals : int;
    nvcsw : int;
    nivcsw : int;
}

external c_start : int -> rt = "CAML_getrusage"

END
    

type t = {
    rt : rt;
    st_time : float;
}

let start () = 
    let a = c_start 0 in
    let b = Unix.gettimeofday () in
    { rt = a; st_time = b }


let get_user st = 
    let cur = start () in
    cur.rt.utime -. st.rt.utime

let get_system st = 
    let cur = start () in
    cur.rt.stime -. st.rt.stime

let wall st = 
    let cur = start () in
    cur.st_time -. st.st_time

external nanosleep : int -> float -> unit = "sleep_CAML_nanosleep"

let to_string time =
    let present remain msg = 
        (Printf.sprintf "%.0f" remain) ^ " " ^ msg
    in
    let mine_remain v factor = 
        let mine = floor (v /. factor) in
        let remain = v -. (mine *. factor) in
        mine, remain
    in
    let rec full_calculation remain =
        if (remain) < 60. then present remain "s"
        else if (remain) < 3600. then (* Show minutes *)
            let mine, remain = mine_remain remain 60. in
            (present mine "m, ") ^ full_calculation remain
        else if (remain) < 3600. *. 24. then (* Show days *)
            let mine, remain = mine_remain remain 3600. in
            (present mine "h, ") ^ full_calculation remain
        else 
            let mine, remain = mine_remain remain (3600. *. 24.) in
            (present mine "d, ") ^ full_calculation remain
    in
    full_calculation time

let status_msg elapsed_wall adv n =
    let tpi = elapsed_wall /. (float_of_int adv) in
    let remain = tpi *. (float_of_int (n - adv)) in
    let remain = to_string remain in
    "Estimated finish in " ^ remain

