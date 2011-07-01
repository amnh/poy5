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

(** Module to align and create medians for dynamic likelihood characters under a
    numbner of different criteras *)

(** {6 Types} Representation of different dynamic likelihood characters *)

type 'a align = { ss : 'a array; }

(** Define the union of all the types of dynamic likelihood characters *)
type r = | FPAlign      of FloatSequence.FloatAlign.s align
(** Characters that transform across a combined time. The assignment is
    restricted to gaps and the states of the children. This model is equivlent
    to the MPL model when dealing with two sequences. The cost is defined as, 
        $P(X,Y|t_1+t_2)$, and the assignment is $X \vee Y$ *)

         | MPLAlign     of FloatSequence.MPLAlign.s align
(** Define characters for maximum parsimonious likelihood or ancestral
    likelihood. These characters find the best assignment of a transformation
    from the two children. If X and Y are the states of the children, $P(a,b|t)$
    is the probability of state b over branch length t, then MPL is defined as,
        $\max_{\alpha = \{ACTG-\}} P(X,\alpha|t_1) * P(Y,\alpha|t_2)$ *)

type t = { model  : FloatSequence.dyn_model;
            data  : r;
            codes : int array;
            code  : int;
            cost  : float; }
(** The type for dynamic likelihood characters; encompasses integerized
    likelihood, floating point alignment, maximum parsimonious likelihood and
    eventually, fixed state likelihood and maximum average likelihood. *)

(** {6 Basic Functions} Implemented in similar modules. *)

val alph : t -> Alphabet.a
(** return the alphabet for the dynamic likelihood characters. *)

val total_cost : t -> float
(** Return the loglikelihood of the sequence; priors are obtained from model in [t] *)

val prior : t -> float
(** Calculate the prior cost of the sequence [t] *)

val get_cm : t -> Cost_matrix.Two_D.m
(** Return the cost matrix used to align [t]; this will fail, as the dynamic
    likelihood charaters use a floating point alignment matrix not compatible
    with the Cost_matrix module *)

val model : t -> FloatSequence.dyn_model
(** Return the likelihood model for the characters *)

val code : t -> int
(** Return the code for the character sequence. *)

val get_codes : t -> int array
(** Return the array of codes for the character partitions. *)

val compare : t -> t -> bool
(** Compare the model and sequence data of two dynamic likelihood characters *)

val to_string : t -> string
(** Simple to_string function for the data-type t *)

val name_string : t -> string
(** Name of the cost function and that we are using dynamic likelihood *)

val cardinal : t -> int
(** The number of partitions in the character set *)

val make : Alphabet.a -> SeqCS.t -> MlModel.model -> t
(** Create the type t used in the module. MlModel contains information regarding
    the type of dynamic likelihood, the alphabet could be different (for
    example, mlModel stores the singular representation, where the code is the
    index, this is different then the bitset version which this alphabet would
    represent. *)

val encoding : (float, Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array1.t -> t -> float
(** Determine the encoding of the sequence characters *)

val leaf_sequences : t -> 
    [> `DO of Sequence.s | `First of Sequence.s | `Last of Sequence.s ] array All_sets.IntegerMap.t
(** Extract the sequence data from the type *)

val to_formatter : Xml.attributes -> t -> t option -> 
            float option * float option -> Data.d -> Xml.xml Sexpr.t list
(** Create XML output of the datatype [t] *)

(** {6 Median Functions} *)

val median : int -> t -> t -> float option -> float option -> t
(** [median c a b ta tb] find the median between [a] and [b] seperated by [ta]
    and [tb] respectively; applying code [c]. *)

val median_i : int -> t -> t -> float -> float -> t * float * float
(** [median c a b ta tb] find the median between [a] and [b] seperated by [ta]
    and [tb] respectively, and optimize each pair by brents method. *)

val median_3 : t -> t -> t -> t -> t
(** [median c a b ta tb] find the median between [a], [b] and [c] seperated by
    [ta], [tb], and [tc] respectively. Not implemented. *)

val readjust : t -> t -> t -> float -> float -> bool * float * float * (float * float) * t
(** [readjust p a b t1 t2] optimize the branch lengths [t1] and [t2] from nodes
    [a] and [b] respectively, updating parent [p]. *)

val readjust3 : t -> t -> t -> t -> float -> float -> float
    -> bool * float * float * (float * float * float) * t

val to_single : t -> t -> float -> float * float * t
(** [to_single p m t] Find the sequence for m as a sequence of assigned states,
    based on [p] and the branch length [t]. *)

(** {6 Distance Functions} *)

val distance : float -> t -> t -> float option -> float
(** [distance f a b] Distance metric for two sequence data; calls the equivlent
    functions in sequences data. *)

val dist_2 : float -> t -> t -> t -> float
(** Heuristic distance of two sequences and a median *)

val tabu_distance : t -> float

val estimate_time : t -> t -> float * float 
(** Estimate the branch lengths between the two sequences using the closed form
    Jukes-Cantor distance solution. *)

(** {6 Filter Functions} *)

val f_codes : t -> All_sets.Integers.t -> t
(** [f_codes t a] Filter the codes in [a] from [t] *)

val f_codes_comp : t -> All_sets.Integers.t -> t
(** [f_codes t a] Filter the codes not in [a] from [t] *)
