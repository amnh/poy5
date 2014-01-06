(* POY 5.1.1. A phylogenetic analysis program using Dynamic Homologies.       *)
(* Copyright (C) 2014 Andrés Varón, Lin Hong, Nicholas Lucaroni, Ward Wheeler,*)
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

(*this is the seed table we are going to use.*)    
val palindromic_spaced_seed_tbl : (int,int list list) Hashtbl.t
   
val fill_in_hmatrix : Cost_matrix.Two_D.m -> unit


val return_a_seedNO : int -> (int array) ref -> unit

val get_a_seedNO : (int array) ref -> int

val radix_sort : (int*int*int*int array*int) array -> (int*int*int*int array*int) array

val extend_seq_in_both_dir : (int*int*int) array -> int -> int array array -> (int*int*int) array * int

val get_score_from_2seq : int array -> int array -> int -> int
