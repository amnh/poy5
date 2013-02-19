(* POY 5.0 Beta. A phylogenetic analysis program using Dynamic Homologies.    *)
(* Copyright (C) 2013 Andrés Varón, Lin Hong, Nicholas Lucaroni, Ward Wheeler,*)
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

type elt = SankCS.elt (*cside*)
type t = SankCS.t (*cside*)

type t_w_seqtbl = {
    states_and_distbl : t;
    sequence_table: Sequence.s array; (*we carry original sequence before
    transform(fixed_state) with us*)
    alph : Alphabet.a;
}

val cardinal : t_w_seqtbl -> int

val f_codes : t_w_seqtbl -> All_sets.Integers.t -> t_w_seqtbl

val f_codes_comp : t_w_seqtbl -> All_sets.Integers.t -> t_w_seqtbl

val get_states : t_w_seqtbl -> int -> int array

val get_earray : t_w_seqtbl -> int array

val get_min_states : t_w_seqtbl -> int * int

val of_array : Data.fixed_state_spec -> int -> int -> t_w_seqtbl 

val distance : t_w_seqtbl -> t_w_seqtbl -> float

val median : int -> t_w_seqtbl -> t_w_seqtbl -> t_w_seqtbl * float

val median_3 : t_w_seqtbl ->  t_w_seqtbl ->  t_w_seqtbl ->  t_w_seqtbl -> t_w_seqtbl 

val dist_2 : t_w_seqtbl -> t_w_seqtbl -> t_w_seqtbl -> float 

val to_single : t_w_seqtbl -> t_w_seqtbl -> t_w_seqtbl * int

val to_string : t_w_seqtbl -> string

val compare_data : t_w_seqtbl -> t_w_seqtbl -> int 

val to_formatter :
    Methods.diagnosis_report_type -> Xml.attributes -> t_w_seqtbl -> Data.d 
    -> Xml.xml Sexpr.t list 


