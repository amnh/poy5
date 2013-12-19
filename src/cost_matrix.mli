(* POY 5.1.1. A phylogenetic analysis program using Dynamic Homologies.       *)
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

(** Cost_matrix

    This Module handles all the operations on two-dimensional and three
    dimensional cost matrices. Given a set of n possible states A = a1 .... an,
    a cost matrix keeps precaulcated costs for tuples (two dimensional cost
    matrix) and triples (three dimensional cost matrix) of these elements.

    The type m in {!Two_D} and {!Three_D} are imperative.  *)


(** {2 Exceptions} *)

exception Illegal_Cm_Format


(** {2 Types} *)

type gap_opening = int

type cost_model = 
    | No_Alignment
    | Linnear
    | Affine of gap_opening 


(** {2 Two Dimensional and Three Dimensional Matrices} *)

module Two_D : sig

    (** Two Dimensional sequence alignment cost matrices. *)


    (** {2 Types} *)

    (** A two dimensional cost matrix type. *)
    type m


    (** {2 Creating} *)

    (** [create a com aff go ge] creates a cost matrix for an alphabet of size a,
        claculating all the possible combinations (effective size 2^a) iff com is
        true, using the affine cost model aff and gap opening and gap extension *)
    val create : int -> bool -> int -> int -> int -> int -> int -> int -> int -> m
    
    (** [clone x] creates a fresh copy of the cost matrix x *)
    val clone : m -> m

    (** [perturbe x y z] perturbates the cost matrix x with severity y and
        probability z *)
    val perturbe : m -> int -> int -> m

    (** [of_list a] creates a fresh transformation cost matrix with the values
        contained in the squared matrix [a]. *)
    val of_list : ?use_comb:bool -> ?level:int -> ?suppress:bool -> 
                    int list list -> int -> m * m

    (** [of_transformations_and_gaps uc as t g] creates a fresh two dimensional
        transformation cost matrix for an alphabet of size [as], with
        transformation cost [t] and indel cost [g]. If [uc] is true, then all the
        combinations for the sequences are calculated, otherwise they are not. *)
    val of_transformations_and_gaps : bool -> int -> int -> int -> int -> m * m

    (** [default] is the default cost matrix in POY for pairwise, nucleotide
        sequence alignments. *)
    val default : m

    val default_nucleotides : m

    val default_aminoacids : m


    (** {2 IO} *)

    (** [of_channel file] parse the file containing a cost matrix and returns
        the processed data. Raise an Illegal_Cm_Format if the format can't be
        parsed. *)
    val of_channel: ?tie_breaker:Methods.keep_method -> ?orientation:bool -> ?use_comb:bool -> ?level:int -> int ->
        FileStream.greader -> m * m * int list list

    (** [of_channel_nocomb file] parse the file containing a cost matrix and
        returns the processed data, but without calculating combinations. Raise an
        Illegal_Cm_Format if the format can't be parsed. *)
    val of_channel_nocomb: 
        ?orientation:bool -> int -> FileStream.greader -> m * m * int list list

    (** [print ma] prints the matrix ma with alphabet size a in stdout. *)
    val output : out_channel -> m -> unit

    val output_constrained : out_channel -> int -> m -> unit


    (** {2 Setting Values} *)

    (** [set_gap cm v] sets the gap representation value to v in the cost matrix cm. *)
    val set_gap : m -> int -> unit

    (* [set_level m v] sets the level value of cost matrix *)
    val set_level : m -> int -> unit

    (** [set_alphabet_size cm v] sets the alphabet size of the cost matrix cm to v.  *)
    val set_alphabet_size : int -> m -> unit

    (** [set_lcm cm v] sets the log2 of the alphabet size in the cost matrix cm to v. *)
    val set_lcm : m -> int -> unit

    (** [set_cost_model cm v] sets the cost model to the appropriate value. *)
    val set_cost_model : m -> cost_model -> unit

    (** [create_cm_by_level m level oldlevel all_elements] creates a new cost
        matrix based on the original matrix and new level. *)
    val create_cm_by_level : m -> int -> int -> int -> Methods.keep_method -> m

    (** [set_cost x y cm v] sets the cost of transforming element x into y in cost
        matrix cm to v*)
    val set_cost : int -> int -> m -> int -> unit

    (** [set_median x y z u] sets the value of the median of states x and y to u
        in the cost matrix z *)
    val set_median : int -> int -> m -> int -> unit


    (** {2 Getting Values} *)

    (** [is_metric mtx] checks if the cost matrix specified by the list [lst]
        (which should be a valid input for [of_list]), is a metric matrix or not. *)
    val is_metric : m -> bool
    
    (*[is_idntity] return true is there is no cost between same states in matrix m*)
    val is_identity : m -> bool

    val get_all_elements : m -> int

    (** [alphabet_size cm] gets the total alphabet size in the cost matrix cm. *)
    val alphabet_size : m -> int

    (** [gap cm] retrieves the gap representation value in cm *)
    val gap : m -> int

    (* [get_ori_a_sz cm] returns the original alphabet size*)
    val get_ori_a_sz : m -> int

    (* [get_level cm] returns the level value of current alphabet *)
    val get_level : m -> int

     (* [get_map_sz cm] returns the map size, thus the number of combiantions of
        current alphabet *)
    val get_map_sz : m -> int

    (** [lcm cm] retrieves the celing of the log2 of the alphabet size of the
        cost matrix cm. *)
    val lcm : m -> int

    val get_tie_breaker : m -> int
    
    val load_file_as_list : FileStream.greader -> int list

    val fill_tail : int array -> m -> unit

    val fill_prepend : int array -> m -> unit

    (** [combine cm] gets the combinations flag in the cost matrix cm. The
        combinations flag stablishes if the median calculation will consider all
        or some of the possible combinations of elements in an original alphabet
        (1) or not (0). *)
    val combine : m -> int

    (** [get_combination] return true if combine=1, false if combine=0 *)
    val get_combination : m -> bool 

    (** [affine cm] retrieves the cost model set in the cost matrix cm. *)
    val get_cost_model : m -> cost_model

    (** [gap_opening cm] retrieves the gap opening cost as set in the cost matrix * cm. *)
    val gap_opening : m -> int

    (** [cost x y cm] retrieves the cost of transforming alphabet element x into y
        according to the cost matrix cm*)
    val cost : int -> int -> m -> int

    (** [worst_cost x y cm] retrieves the WORST cost of transforming alphabet
        element x into y according to the cost matrix cm. *)
    val worst_cost : int -> int -> m -> int

    (** [median x y cm] retrieves the median for transforming character x into y
        according to the cost matrix cm. *)
    val median : int -> int -> m -> int

    val get_gap_startNO : m -> int

    (** [states_of_code code cm] return the list of states in an alphabet from a
        state code from a median/leaf *)
    val states_of_code : int -> m -> int list

    (** [code_of_states states cm] return the code from a list of states. This
        may not work if the level is not high enough. *)
    val code_of_states : int list -> m -> int

    (** [get_pure_cost_mat cost_mat] return the cost matrix in an 
        int alphabet_size * alphabet_size matrix format *)  
    val get_pure_cost_mat : m -> int array array

    (** return a reduced matrix; it's essence, the individual costs *)
    val ori_cm_to_list : m -> int list

    (** [get_closest m x y], finds the closest element from the bitset [y] to
        the element [x] according to the two dimensional cost matrix [m]. [x]
        can be a bitset, and [y] need not to be a bitset. *)
    val get_closest : m -> int -> int -> int

    (** [calc_number_of_combinations_by_level a_sz level] returns the number of
        combinations based on alphabet size:a_sz and level value:level *)
    val calc_number_of_combinations_by_level: int -> int -> int 

    (** [calc_num_of_comb_with_gap ori_a_sz level] returns the number of
        combination code that contains gap *)
    val calc_num_of_comb_with_gap : int -> int -> int
    
    (** [gap_filter_for_combcode combcode level ori_a_sz], get rid of "gap" in a
        combination code if the code contains gap code. for example, if gap is
        "-", then input [a-] will get [a] *)
    val gap_filter_for_combcode : int -> m -> int

    (** [comblist_to_combcode lst m], returns the combination code given a
        combiantion codelist. This function should be more efficent if we use
        Bigarray -- add Bigarray later. *)
    val comblist_to_combcode: int list  -> m -> int

    (** [combcode_to_comblist code m], returns combination codelist given a
        combination code *)
    val combcode_to_comblist: int -> m -> int list

    (** [clear_duplication_in_list lst] clears up the duplicate element in the
        given list *)
    val clear_duplication_in_list: int list -> int list

    (** [print_intlist list] prints out the int list, for debug....*)
    val print_intlist: int list -> unit

    val of_file : ?tie_breaker:Methods.keep_method -> ?orientation:bool ->
        ?use_comb:bool -> ?level:int -> FileStream.f -> int -> bool ->
            m * m * int list list

    (* [matrix_of_file fn file] Read a file into an array array, and map a
       function over the values; we ensure that the matrix is rectangular. *)
    val matrix_of_file: (string -> 'a) -> FileStream.f -> 'a list list

    (** [check_level return true if cost matrix is using level]*)
    val check_level : m -> bool
end


module Three_D : sig

    (** Three Dimensional sequence alignment cost matrices. *)

    (** {2 Types} *)

    (** A two dimensional cost matrix type. *)
    type m


    (** {2 Creating} *)

    (** [create a_sz comb aff go ge dim] creates a new cost matrix with alphabet
        size a_sz, considering all the possible combinations of members of the
        alphabet (total efective size 2^a_sz - 1), if comb is true with affine
        gap cost model if aff is not 0, in which case will use gap opening cost
        go and extension gap ge. The gap is represented as the next available
        integer depending on the a_sz. IF dim is true the matrix will be three
        dimensional, otherwise it will be two dimensional. *)
    val create : int -> bool -> int -> int -> int -> int -> int -> int -> int -> m

    (** [clone x] creates a fresh copy of the cost matrix x *)
    val clone : m -> m

    (** [of_two_dim cm] creates a fresh three dimensional cost matrix using the
        values stored in the two dimensional matrix cm. This is the only way to
        create a three dimensional cost matrix. *)
    val of_two_dim: Two_D.m -> m

    (** The default three dimensional cost matrix in Poy. It's equivalent to
        [of_two_dim Two_D.default]. *)
    val default : m

    val default_nucleotides : m

    val default_aminoacids : m Lazy.t

    (** [perturbe x y z] perturbates the cost matrix x with severity y and
        probability z *)
    val perturbe : m -> int -> int -> m


    (** {2 IO } *)

    (** [print ma] prints the matrix ma with alphabet size a in stdout. *)
    val output : out_channel -> m -> unit


    (** {2 Setting Values} *)

    (** [set_gap cm v] sets the gap representation value to v in the cost matrix cm. *)
    val set_gap : m -> int -> unit

    (** [set_alphabet_size cm v] sets the alphabet size of the cost matrix cm to v. *)
    val set_alphabet_size : int -> m -> unit

    (** [set_lcm cm v] sets the log2 of the alphabet size in the cost matrix cm to v. *)
    val set_lcm : m -> int -> unit

    (** [set_cost_model cm v] sets the cost model to the appropriate value. *)
    val set_cost_model : m -> cost_model -> unit

    (** [set_cost x y z cm v] sets the cost of transforming element x, y and z 
        in cost matrix cm to v *)
    val set_cost : int -> int -> int -> m -> int -> unit

    (** [set_median x y z u v] sets the value of the median of states x, y and z 
        to v in the cost matrix u *)
    val set_median : int -> int -> int -> m -> int -> unit


    (** {2 Getting Values} *)

    (** [alphabet_size cm] gets the total alphabet size in the cost matrix cm. *)
    val alphabet_size : m -> int

    (** [gap cm] retrieves the gap representation value in cm *)
    val gap : m -> int

    (** [lcm cm] retrieves the celing of the log2 of the alphabet size of the cost
        matrix cm. *)
    val lcm : m -> int

    (** [combine cm] gets the combinations flag in the cost matrix cm. The
        combinations flag stablishes if the median calculation will consider all
        the possible combinations of elements in an original alphabet (1) or not *)
    val combine : m -> int

    (** [affine cm] retrieves the cost model set in the cost matrix cm. *)
    val affine : m -> int

    (** [gap_opening cm] retrieves the gap opening cost as set in the cost matrix cm. *)
    val gap_opening : m -> int

    (** [cost x y z cm] retrieves the cost of transforming alphabet element x, y
        and z simultaneusly according to the cost matrix cm. *)
    val cost : int -> int -> int -> m -> int

    (** [median x y z cm] retrieves the median for transforming character x, y and
        z simultaneusly according to the cost matrix cm. *)
    val median : int -> int -> int -> m -> int
end
