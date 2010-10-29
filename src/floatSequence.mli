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

val to_char_list : string -> char list
val sequence_of_string : string -> Alphabet.a -> Sequence.s

module FloatAlign : sig

    type floatmem
    type s

    val get_mem : s -> s -> floatmem
    val s_of_seq: Sequence.s -> s 
    val seq_of_s: s -> Sequence.s

    val cost_2          : ?debug:bool -> ?deltaw:int -> s -> s -> MlModel.model -> float -> floatmem -> float
    val verify_cost_2   : float -> s -> s -> MlModel.model -> float -> floatmem -> float
    val c_cost_2        : s -> s -> MlModel.model -> float -> floatmem -> int -> float
    val create_edited_2 : s -> s -> MlModel.model -> float -> floatmem -> s * s
    val align_2         : ?first_gap:bool -> s -> s -> MlModel.model -> float -> floatmem -> s * s * float
    val median_2        : s -> s -> MlModel.model -> float -> floatmem -> s
    val median_2_cost   : s -> s -> MlModel.model -> float -> floatmem -> float * s
    val full_median_2   : s -> s -> MlModel.model -> float -> floatmem -> s
    val closest         : s -> s -> MlModel.model -> float -> floatmem -> s * float
    val readjust        : s -> s -> s -> MlModel.model -> float -> float -> float -> floatmem -> float * s * bool

    (* delete these later *)
    val print_mem : floatmem -> unit
    val print_s   : s -> unit
    val clear_mem : floatmem -> unit
    val print_cm  : MlModel.model -> float -> unit

end
