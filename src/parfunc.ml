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

(* $Id: parfunc.ml 1644 2007-02-14 19:05:47Z andres $ *)
let () = SadmanOutput.register "Parfunc" "$Revision: 1644 $"


let all_minimize f c a =
    match Comm.get_exec_mode () with
 (*   | Comm.Non_Fault_Tolerant, Comm.Synch -> Paroper.all_minimize f c a;
    | Comm.Non_Fault_Tolerant, Comm.Asynch -> Paroperas.all_minimize f c a;*)
    | Comm.Fault_Tolerant, Comm.Asynch -> Asynch.all_minimize f c a;
    | _, _ -> raise  (Comm.Error "The execution mode is not supported. \
        And this is not supposed to be called!");;

let all_maximize f c a =
    match Comm.get_exec_mode () with
 (*   | Comm.Non_Fault_Tolerant, Comm.Synch -> Paroper.all_maximize f c a;
    | Comm.Non_Fault_Tolerant, Comm.Asynch -> Paroperas.all_maximize f c a; *)
    | Comm.Fault_Tolerant, Comm.Asynch -> Asynch.all_maximize f c a;
    | _, _ -> raise  (Comm.Error "The execution mode is not supported. \
        And this is not supposed to be called!");;

let gather_minimum f c a th =
    match Comm.get_exec_mode () with
(*  | Comm.Non_Fault_Tolerant, Comm.Synch -> Paroper.gather_minimum f c a th;
    | Comm.Non_Fault_Tolerant, Comm.Asynch -> Paroperas.gather_minimum f c a th;*)
    | Comm.Fault_Tolerant, Comm.Asynch -> Asynch.gather_minimum f c a th;
    | _, _ -> raise  (Comm.Error "The execution mode is not supported. \
        And this is not supposed to be called!");;

let gather_maximum f c a th =
    match Comm.get_exec_mode () with
(*  | Comm.Non_Fault_Tolerant, Comm.Synch -> Paroper.gather_maximum f c a th;
    | Comm.Non_Fault_Tolerant, Comm.Asynch -> Paroperas.gather_maximum f c a th;*)
    | Comm.Fault_Tolerant, Comm.Asynch -> Asynch.gather_maximum f c a th;
    | _, _ -> raise  (Comm.Error "The execution mode is not supported. \
        And this is not supposed to be called!");;

let map f a =
    match Comm.get_exec_mode () with
(*  | Comm.Non_Fault_Tolerant, Comm.Synch -> Paroper.map f a;
    | Comm.Non_Fault_Tolerant, Comm.Asynch -> Paroperas.map f a;*)
    | Comm.Fault_Tolerant, Comm.Asynch -> Asynch.map f a;
    | _, _ -> raise  (Comm.Error "The execution mode is not supported. \
        And this is not supposed to be called!");;

let scan f a =
    match Comm.get_exec_mode () with
(*  | Comm.Non_Fault_Tolerant, Comm.Synch -> Paroper.scan f a;
    | Comm.Non_Fault_Tolerant, Comm.Asynch -> Paroperas.scan f a; *)
    | Comm.Fault_Tolerant, Comm.Asynch -> Asynch.scan f a;
    | _, _ -> raise  (Comm.Error "The execution mode is not supported. \
    And this is not supposed to be called!");;



(* vim:sw=4 et tw=80
 *)
