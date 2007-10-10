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

(* $Id: parparam.mli 1644 2007-02-14 19:05:47Z andres $ *)

(** Parallelization parameters calculation. 
*
* Operations to calculate (and in the future to learn) parameters for jobs
* parallelization. *)

exception Impossible of string

(** [cal_range s] calculates the portion of a job of size s that is reponsability
* of the calling process. The function returns a Impossible s condition if s is
* less than the total number of processors. Note that this is a recoverable
* error by applying a partition in the first s processes, using the cal_spread
* function. The function returns the index of the first element to be processed
* and the number of elements to be processed. If the size of the array is less
* than the number of processes, an Impossible exception is raised.*)
val cal_range: int -> (int * int)

(** [cal_spread l] returns a list of l / processes elements to be executed by the
* calling process. The difference with cal_range is that there are not error
* messages if the total number of processes is more than the number of tasks l
* to be performed.*)
val cal_spread: int -> Mpi.rank list 

(* vim:sw=4 et tw=80
 *)
