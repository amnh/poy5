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

(** Floating Point Alignment Modules *)

val sequence_of_string : ?filter_gap:bool -> string -> Alphabet.a -> Sequence.s
(** Convert a string to a sequence based on the provided alphabet *)

val remove_gaps : int -> Sequence.s -> Sequence.s
(** filter out gaps and prepend sequence with an opening gap *)

type dyn_model = { static : MlModel.model; alph : Alphabet.a; }
(** type that builds on the static likelihood model. The alphabet contained in
    the static likelihood is sequential, representing the indexes of the
    characters in the likelihood matrix. This type adds the unsequential
    alphabet to the model for use as bitsets with the vectors. *)

val make_model : Alphabet.a -> MlModel.model -> dyn_model
(** simple function to wrap up composing a dynamic likelihood model *)

val cost_fn : dyn_model -> Methods.ml_costfn
(** Return the cost function for the dynamic likelihood model. This helps in
    retrieving from other modules since the name would be quite long. *)

val test_all : bool -> out_channel -> Sequence.s -> Sequence.s -> float -> float -> dyn_model -> unit
(** Run an alignment with each of the methods below with the specified branch
    lengths and model. Also, optimize the alignment branch lengths and print out
    all the scores. The optimized scores for MPL and FLK should be equal; this
    is based on the pully principle. *)

module type A = sig

    (** {6 Types *)

    type floatmem
    (** [floatmem] The representation of the memory scratch space used in
        alignment of the sequence data *)

    type s
    (** [s] The type that represents the sequence data of the module; this can
        be converted to the parsimony characters by [seq_of_s], and converted
        from it conversely. *)

(** {6 Auxiliary/helper functions *)

    val get_mem     : s -> s -> floatmem
    (** Return scratch space to be used. The returned space is checked to fix [s] and [s]. *)

    val create_mem  : int -> int -> floatmem
    (** [create_mem w h] Initially create a scratch space of [w] [h] *)

    val clear_mem   : floatmem -> unit
    (** [clear_mem f] Clear the memory by setting all the values in the matrix
        to -1.0. This can be used with debug settings to ensure cells being used
        in the calculation of an alignment are all positive. *)

    (* {6 Conversion of datatype *)

    val s_of_seq    : Sequence.s -> s 
    (** [s_of_seq seq] Convert the Sequence data [seq] to the type s defined in the module *)

    val seq_of_s    : s -> Sequence.s
    (** [seq_of_s s] Convert the FloatSequence data [s] to the Sequence.s type *)

    val compare     : s -> s -> bool
    (** [compare a b] Compare the two sequences data *)

    (** {6 Functions for testing externally *)

    val print_mem   : floatmem -> unit
    (** [print_mem mem] Print the memory --alignment matrix-- to stdout *)

    val print_s     : s -> Alphabet.a -> unit
    (** [print_s s a] Print the sequence resolving the characters with the alphabet defined by [a]. *)

    val print_raw   : s -> unit
    (** [print_raw s] Print the sequence to stdout. *)

    val print_cm    : dyn_model -> float -> unit
    (** [print_cm m t] Compose and print the cost matrix defined by the model
        [m] and the branch length [t]. *)

(* {6 Cost Function/Matrix functions *)

    val get_cf : dyn_model -> float -> float -> (int -> int -> float * int)
    (** [get_cf m at bt] -> [cost i j] Create a function that defines the cost
        of converting from character [i] to character [j], based on the composed
        matrix from the model [m] and branch lengths [at] and [bt]. *)

    val get_cm : dyn_model -> float -> float -> (float * int) array array
    (** [get_cm m at bt] -> [cost i j] Create a matrix that defines the cost
        of converting from character [i] to character [j], based on the composed
        matrix from the model [m] and branch lengths [at] and [bt]. *)

    val cost   : dyn_model -> float -> float -> (int -> int -> float * int)
    (** [cost m t1 t2 -> x y -> cost] create a function that returns the cost
        between bases as integers. *)

    val aln_cost_2      : s -> s -> dyn_model -> float -> float
    (** [alncost_2 a b m t] Determines the cost of the branch from [a] to [b]
        and branch length [t], of model [m]. [a] and [b] must be aligned. *)

(* {6 2D Alignment Operations *)

    val cost_2          : ?deltaw:int -> s -> s -> dyn_model -> float -> float -> floatmem -> float
    (** [cost_2 ?delta a b m at bt mem] Return the cost of the median, from the
        alignment of [a] and [b] with branch lengths [at] and [bt] respectively,
        and model [m]. The [mem] is updated, and can be used to return a
        backtrace or the edited sequences, although it is recommended that
        another function that returns that data be called instead. *)

    val full_cost_2     : float -> s -> s -> dyn_model -> float -> float -> floatmem -> float
    (** [verify_cost_2 cost a b m at bt mem] Return the cost of the alignment of
        [a] and [b] with model [m], branch lengths [at] and [bt], respectivly.
        This alignment does a full alignment, and tests the final cost against
        [cost]; no action is done to if they are different if the debug
        parameter isn't set. *)

    val create_edited_2 : s -> s -> dyn_model -> float -> float -> floatmem -> s * s
    (** [create_edited_2 a b m at bt mem] Create the edited sequences of [a] and
        [b] from their alignment with branch lengths [at] and [bt], model [m]. *)

    val align_2         : ?first_gap:bool -> s -> s -> dyn_model -> float -> float -> floatmem -> s * s * float
    (** [align_2 a b m at bt mem] Create the edited sequences of [a] and [b]
       from their alignment with branch lengths [at] and [bt], model [m]. Also
       return the cost of the alignment. *)

    val clip_align_2    : ?first_gap:bool -> Sequence.Clip.s -> Sequence.Clip.s
                            -> dyn_model -> float -> float 
                                -> Sequence.Clip.s * Sequence.Clip.s * float * int * Sequence.Clip.s * Sequence.Clip.s
    (** [clip_align_2 a b m at bt] Convert the clip data to S data, and
        align, returning cost, length, and clipped and unclipped alignments.
        Although, this module does not deal with clips, but all input should be
        `DO, and output should be the same, thus a - anoclip. Used in the
        implied alignment module to help with processing *)

    val median_2        : s -> s -> dyn_model -> float -> float -> floatmem -> s
    (** [median_2 a b m at bt mem] Create the backtrace of the mem, from an
        alignment of [a] and [b]. *)

    val median_2_cost   : s -> s -> dyn_model -> float -> float -> floatmem -> float * s
    (** [median_2_cost a b m t1 t2 mem] Create the backtrace of the mem, from an
        alignment of [a] and [b]. Return the cost of the median as well. *)
    
    val full_median_2   : s -> s -> dyn_model -> float -> float -> floatmem -> s
    (** [full_median_2 a b m t1 t2 mem] Create the backtrace through a full
        alignment of the data; avoids the Ukkonen Approximate String Matching
        Algorithm. Similar to [full_cost_2]. *)

    val gen_all_2       : s -> s -> dyn_model -> float -> float -> floatmem -> s * s * float * s
    (** [gen_all_2 a b m t1 t2 mem] Create the edited sequences of [a] and [b]
        from the alignment along the branches of length [at] and [bt]
        respectively, and the median associated with the edited distances and
        costs. The length of all the returned sequences are equal. *)

    val backtrace       : ?filter_gap:bool -> dyn_model -> floatmem -> s -> s -> s
    (** [backtrace ?fg m mem a b] Produce the backtrace from an alignment. This
        can be used mostly in testing and probably should not be called
        directly. More appropriate functions are median_2/full_median_2/gen_all_2 *)

(** {6 Pseudo 3d operations *)

    val closest  : p:s -> m:s -> dyn_model -> float -> floatmem -> s * float
    (** [closest p m model t mem] Create a sequence from [m] using [p] as its
        parent that is the closest assignment of [m] from [p] over length [t].
        The returned sequence will have no polymorphisms. *)

    val get_closest : dyn_model -> float -> i:int -> p:int -> m:int -> int * float
    (** [get_closest model t i p m] Find the state of [m] that would create
        the minimum cost to [p]; return the optimal state, and cost associated
        with that transformation based on the model and branch length, [t]. *)

    val readjust : s -> s -> s -> dyn_model -> float -> float -> float -> floatmem -> float * s * bool
    (** [readjust a b c m at bt ct mem] Perform a pseudo 3D alignment by using
        the best score of any pair, and then performing [closest] on that median
        with the third sequence. *)

    val optimize : s -> s -> dyn_model -> float -> floatmem -> float * float
    (** [optimize a b m t mem] find the minimum cost by adjusting the branch
        length between [a] and [b]; with initial value at [t]. Return the pair,
        cost * branch length. **)

end 
(** The sequence alignment module for floating point cost matrices/regimes. *)

module FloatAlign : A
(** This module calculates the alignment (and the cost of the alignment), by
    combining the branch lengths, and creating the transition matrix from that
    over-all distance. The assignment does not deviate from the domain of the
    children --polymorphisms happen from a substitution only-- and is assigned
    the most likely state on that assumption. *)

module MPLAlign   : A
(** This module determines two costs matrices, and finds the median by the
    likelihood of transform either of the children over their cost matrix to an
    assignment median. Each assignment is considered. 

    PROD (x,y=0 to n,m) of MAX (i in A) of ( P(t1)_xi * P(t2)_yi )
    Where A is the set of characters in the alphabet, and P(t) = e^Qt. *)

module MALAlign   : A
(** This is a partially implemented module to define the distance between two
    sequences by the criteria of total likelihood ->dynamic MAL. Only the cost
    function is defined along with the memory and cost matrix functions.
    Alignments are not defined (as the total cost is the sum; finding an
    alignment would be equivlent to the FloatAlign module. *)

val pair_distance : dyn_model -> Sequence.s -> Sequence.s -> float * float
