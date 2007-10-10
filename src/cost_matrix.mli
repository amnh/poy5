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

(** Cost_matrix 
*
* This Module handles all the operations on two-dimensional and three
* dimensional cost matrices. Given a set of n possible states A = a1 .... an, 
* a cost matrix keeps precaulcated costs for tuples (two dimensional cost 
* matrix) and triples (three dimensional cost matrix) of these elements. 
*
* The type m in {!Two_D} and {!Three_D} are imperative.
* *)

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
    * claculating all the possible combinations (effective size 2^a) iff com is
    * true, using the affine cost model aff and gap opening go and gap 
    * extension ge.  *)
    external create : int -> bool -> int -> int -> m = "cm_CAML_create"

    (** [clone x] creates a fresh copy of the cost matrix x *)
    external clone : m -> m = "cm_CAML_clone"

    (** [perturbe x y z] perturbates the cost matrix x with severity y and
    * probability z *)
    val perturbe : m -> int -> int -> m

    (** [of_list a] creates a fresh transformation cost matrix with the values
    * contained in the squared matrix [a]. *)
    val of_list : ?use_comb:bool -> int list list -> m

    (** [of_transformations_and_gaps uc as t g] creates a fresh two dimensional
    * transformation cost matrix for an alphabet of size [as], with
    * transformation cost [t] and indel cost [g]. If [uc] is true, then all the
    * combinations for the sequences are calculated, otherwise they are not. *)
    val of_transformations_and_gaps : bool -> int -> int -> int -> m

    (** [default] is the default cost matrix in POY for pairwise, nucleotide
    * sequence alignments. *)
    val default : m
    val default_nucleotides : m

    val default_aminoacids : m

    (** {3 IO} *)

    (** [of_channel file] parse the file containing a cost matrix and returns
    the processed data. Raise an Illegal_Cm_Format if the format can't be
    parsed. *)
    val of_channel: ?orientation:bool -> ?use_comb:bool ->
        FileStream.greader -> m

    (** [of_channel_nocomb file] parse the file containing a cost matrix and
        returns the processed data, but without calculating combinations. Raise an
        Illegal_Cm_Format if the format can't be parsed. *)
    val of_channel_nocomb: ?orientation:bool -> FileStream.greader -> m

    (** [print ma] prints the matrix ma with alphabet size a in stdout. *)
    val output : out_channel -> m -> unit

    (** {2 Setting Values} *)

    (** [set_gap cm v] sets the gap representation value to v in the cost matrix
    * cm. *)
    external set_gap : m -> int -> unit = "cm_CAML_set_gap"

    (** [set_alphabet_size cm v] sets the alphabet size of the cost matrix cm to
    * v.  *)
    external set_alphabet_size : int -> m -> unit = "cm_CAML_set_a_sz"

    (** [set_lcm cm v] sets the log2 of the alphabet size in the cost matrix cm to
    * v. *)
    external set_lcm : m -> int -> unit = "cm_CAML_set_lcm"

    (** [set_affine cm v] sets the cost model to the appropriate value. *)
    val set_affine : m -> cost_model -> unit

    (** [set_cost x y cm v] sets the cost of transforming element x into y in cost
    * matrix cm to v*)
    external set_cost : int -> int -> m -> int -> unit = "cm_CAML_set_cost"

    (** [set_median x y z u] sets the value of the median of states x and y to u
    * in the cost matrix z *)
    external set_median : 
        int -> int -> m -> int -> unit = "cm_CAML_set_median"

    (** {2 Getting Values} *)
    (** [is_metric mtx] checks if the cost matrix specified by the list
    * [lst] (which should be a valid input for [of_list]), is a metric matrix or
    * not. *)
    val is_metric : m -> bool

    (** [alphabet_size cm] gets the total alphabet size in the cost matrix cm. *)
    external alphabet_size : m -> int = "cm_CAML_get_a_sz"

    (** [gap cm] retrieves the gap representation value in cm *)
    external gap : m -> int = "cm_CAML_get_gap"

    (** [lcm cm] retrieves the celing of the log2 of the alphabet size of the cost
    * matrix cm. *)
    external lcm : m -> int = "cm_CAML_get_lcm"

    val load_file_as_list : FileStream.greader -> int list
    val fill_tail : int array -> m -> unit 
    val fill_prepend : int array -> m -> unit 
    (** [combine cm] gets the combinations flag in the cost matrix cm. The
    * combinations flag stablishes if the median calculation will consider all
    * the possible combinations of elements in an original alphabet (1) or not
    * (0). *)
    external combine : m -> int = "cm_CAML_get_combinations"

    (** [affine cm] retrieves the cost model set in the cost matrix cm. *)
    val affine : m -> cost_model

    (** [gap_opening cm] retrieves the gap opening cost as set in the cost matrix
    * cm. *)
    external gap_opening : m -> int = "cm_CAML_get_gap_opening"

    (** [cost x y cm] retrieves the cost of transforming alphabet element x into y
    * according to the cost matrix cm*)
    external cost : int -> int -> m -> int = "cm_CAML_get_cost"

    (** [median x y cm] retrieves the median for transforming character x into y
    * according to the cost matrix cm. *)
    external median : int -> int -> m -> int = "cm_CAML_get_median"


    (** [get_pure_cost_mat cost_mat] return the cost matrix in an 
     * int alphabet_size * alphabet_size matrix format *)  
    val get_pure_cost_mat : m -> int array array
    val list_of_bits : int -> int -> int list

    (* [get_closest m x y], finds the closest element from the bitset [y] to the
    * element [x] according to the two dimensional cost matrix [m]. [x] can be a
    * bitset, and [y] need not to be a bitset. *)
    val get_closest : m -> int -> int -> int
end

module Three_D : sig
    (** Three Dimensional sequence alignment cost matrices. *)

    (** {2 Types} *)

    (** A two dimensional cost matrix type. *)
    type m

    (** {2 Creating} *)

    (** [create a_sz comb aff go ge dim] creates a new cost matrix with
    * alphabet size a_sz, considering all the possible combinations
    * of members of the alphabet (total efective size 2^a_sz - 1),
    * if comb is true with affine gap cost model if aff is not 0, in
    * which case will use gap opening cost go and extension gap ge. The
    * gap is represented as the next available integer depending on
    * the a_sz. IF dim is true the matrix will be three dimensional,
    * otherwise it will be two dimensional. *)
    external create : 
        int -> bool -> int -> int -> int -> m = "cm_CAML_create_3d"

    (** [clone x] creates a fresh copy of the cost matrix x *)
    external clone : m -> m = "cm_CAML_clone_3d"

    (** [of_two_dim cm] creates a fresh three dimensional cost matrix using the
    * values stored in the two dimensional matrix cm. This is the only way to
    * create a three dimensional cost matrix. *)
    val of_two_dim: Two_D.m -> m

    (** The default three dimensional cost matrix in Poy. It's equivalent to
    * [of_two_dim Two_D.default]. *)
    val default : m
    val default_nucleotides : m
    val default_aminoacids : m Lazy.t
    (** [perturbe x y z] perturbates the cost matrix x with severity y and
    * probability z *)
    val perturbe : m -> int -> int -> m

    (** {3 IO } *)

    (** [print ma] prints the matrix ma with alphabet size a in stdout. 
    * *)
    val output : out_channel -> m -> unit

    (** {2 Setting Values} *)

    (** [set_gap cm v] sets the gap representation value to v in the cost matrix
    * cm. *)
    external set_gap : m -> int -> unit = "cm_CAML_set_gap_3d"

    (** [set_alphabet_size cm v] sets the alphabet size of the cost matrix cm to
    * v. *)
    external set_alphabet_size : int -> m -> unit = "cm_CAML_set_a_sz_3d"

    (** [set_lcm cm v] sets the log2 of the alphabet size in the cost matrix cm to
    * v. *)
    external set_lcm : m -> int -> unit = "cm_CAML_set_lcm_3d"

    (** [set_affine cm v] sets the cost model to the appropriate value. *)
    val set_affine : m -> cost_model -> unit

    (** [set_cost x y z cm v] sets the cost of transforming element x, y and z 
    * in cost matrix cm to v *)
    external set_cost :
        int -> int -> int -> m -> int -> unit = "cm_CAML_set_cost_3d"

    (** [set_median x y z u v] sets the value of the median of states x, y and z 
    * to v in the cost matrix u *)
    external set_median : 
        int -> int -> int -> m -> int -> unit = "cm_CAML_set_median_3d"

    (** {2 Getting Values{ *)

    (** [alphabet_size cm] gets the total alphabet size in the cost matrix cm. *)
    external alphabet_size : m -> int = "cm_CAML_get_a_sz_3d"

    (** [gap cm] retrieves the gap representation value in cm *)
    external gap : m -> int = "cm_CAML_get_gap_3d"

    (** [lcm cm] retrieves the celing of the log2 of the alphabet size of the cost
    * matrix cm. *)
    external lcm : m -> int = "cm_CAML_get_lcm_3d"

    (** [combine cm] gets the combinations flag in the cost matrix cm. The
    * combinations flag stablishes if the median calculation will consider all
    * the possible combinations of elements in an original alphabet (1) or not
    * (0). *)
    external combine : m -> int = "cm_CAML_get_combinations_3d"

    (** [affine cm] retrieves the cost model set in the cost matrix cm. *)
    val affine : m -> int

    (** [gap_opening cm] retrieves the gap opening cost as set in the cost matrix
    * cm. *)
    external gap_opening : m -> int = "cm_CAML_get_gap_opening_3d"

    (** [cost x y z cm] retrieves the cost of transforming alphabet element x, y
    * and z simultaneusly according to the cost matrix cm. *)
    external cost : int -> int -> int -> m -> int = "cm_CAML_get_cost_3d"

    (** [median x y z cm] retrieves the median for transforming character x, y and
    * z simultaneusly according to the cost matrix cm. *)
    external median : int -> int -> int -> m -> int = "cm_CAML_get_median_3d"


end
