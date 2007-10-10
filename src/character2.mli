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

exception Illegal_argument of string

val assign_codes : 
    int -> 'a array list -> 'b list list -> int array list * int list * int

(** Operations on sets of non additive characters. *)
module NonAdditive : sig

    (** A character state *)
    type c = int

    (** The set of non additive characters. *)
    type s 

    (** [create x y z] creates a fresh non additive character set using the 
    * encoding specifications [x] and the codes assigned to each correspoding 
    * element in [y], for the taxon character set [z]. *)
    val create : 
        Parser.Hennig.Encoding.s array -> int array -> Parser.t array -> s

    (** [set_preliminary x y z] sets the preliminary state of the character in
     * position [x] to [y] in the set [z]. *)
    val set_preliminary : int -> c -> s -> unit

    (** [get_preliminary x y] get's the preliminary state of the character in
     * position [x] in the set [y]. *)
    val get_preliminary : int -> s -> c

    (** [set_final x y z] sets the final state of the character in
     * position [x] to [y] in the set [z]. *)
    val set_final : int -> c -> s -> unit

    (** [get_final x y] get's the final state of the character in
     * position [x] in the set [y]. *)
    val get_final : int -> s -> c

    (** [set_cost x y z] sets the cost of the transformation of the state in
     * position [x] to [y] in the two children of the node (HTU) that holds the
     * set [z]. *)
    val set_cost : int -> int -> s -> unit

    (** [get_cost x y] gets the cost of the transformation of the state in
     * position [x] in the two children of the node (HTU) that holds the
     * set [y]. *)
    val get_cost : int -> s -> int

    (** [set_weight x y z] sets the weight of the characater in position [x] to
    * [y] in the set [z]. *)
    val set_weight : int -> int -> s -> unit

    (** [get_weight x y z] gets the weight of the characater in position [x] 
    * in the set [z]. *)
    val get_weight : int -> s -> int

    (** [set_total x y z] sets the total cost of the transformation of the 
    * state in position [x] to [y] in the two children trees of the node (HTU) 
    * that holds the set [z]. *)
    val set_total : int -> int -> s -> unit

    (** [get_total x y] gets the total cost of the transformation of the state in
     * position [x] in the two children trees of the node (HTU) that holds the
     * set [y]. *)
    val get_total : int -> s -> int

    (** [set_char_code x y z] sets to [y] the character code of the character in
     * position [x] of the set [z]. *)
    val set_char_code : int -> int -> s -> unit

    (** [get_char_code x y] gets the character code of the character in
     * position [x] of the set [z]. *)
    val get_char_code : int -> s -> int

    (** [set_changed x y] sets to true the changed flag of the character in
     * position [x] of the set [y]. *)
    val set_changed : int -> s -> unit

    (** [set_changed x y] sets to false the changed flag of the character in
     * position [x] of the set [y]. *)
    val set_no_changed : int -> s -> unit

    (** [get_changed x y] gets the changed flag of the character in position [x]
     * of the character set [y]. *)
    val get_changed : int -> s -> bool

    (** [get_index_of_code x y] gets the index of the character with code [x] in
     * the character set [y] *)
    val get_index_of_code : int -> s -> int

    (** [downpass x y z] performs a downpass operation on the character set [x]
     * and [y] and stores the result in the character set [z]. The function
     * returns the total cost of the downpass step *)
    val downpass : s -> s -> s -> int

    (** [uppass x y z u] performs an uppass that sets the final states of the
    * character set [u] whith [x] being it's parent and [y] and [z] it's
    * children. The function assumes that the final states of [x] are already
    * set.  *)
    val uppass : s -> s -> s -> s -> unit

    (* [preliminary x] sets the preliminary states of [x] as the final. This
    * function should be called on the root of a tree before a full uppass.*)
    val preliminary_is_final : s -> unit

    (* [join x y] takes two disjoint sets of characters [x] and [y] and outputs
    * a fresh character set. *)
    val join : s -> s -> s

    (* [split x] takes a single character ste [x] and outputs two character sets
    * with half of it (or the floor and ceiling of the half of it) *)
    val split : s -> s * s

end

(** Operations on sets of additive characters. Same functionality of NonAdditive
* *)
module Additive : sig

    (** A character state *)
    type c = int * int

    type s 

    val create : 
        Parser.Hennig.Encoding.s array -> int array -> Parser.t array -> s 

    val set_preliminary : int -> int * int -> s -> unit

    val get_preliminary : int -> s -> int * int

    val set_final : int -> int * int -> s -> unit

    val get_final : int -> s -> int * int

    val set_cost : int -> int -> s -> unit

    val get_cost : int -> s -> int

    val set_weight : int -> int -> s -> unit

    val get_weight : int -> s -> int

    val set_total : int -> int -> s -> unit

    val get_total : int -> s -> int

    val set_code : int -> int -> s -> unit

    val get_code : int -> s -> int

    val set_changed : int -> s -> unit

    val set_no_changed : int -> s -> unit

    val get_changed : int -> s -> bool

    val get_index_of_code : int -> s -> int

    val downpass : s -> s -> s -> unit

    val uppass : s -> s -> s -> s -> unit

    val preliminary_is_final : s -> unit

    val set_all_weights : int -> s -> unit

    val remove_weights : s -> unit

    val join : s -> s -> s

    val split : s -> s * s

end
(*
module Sankoff_Precomputed : sig

    (** A character state *)
    type c = int

    (** The set of sankof characters. *)
    type s 

    (** [create x y z] creates a fresh sankoff character set using the 
    * encoding specifications [x] and the codes assigned to each correspoding 
    * element in [y], for the taxon character set [z]. *)
    val create : 
        Parser.Hennig.Encoding.s array -> int array -> Parser.t array -> s

    (** [set_preliminary x y z] sets the preliminary state of the character in
     * position [x] to [y] in the set [z]. *)
    val set_preliminary : int -> c -> s -> unit

    (** [get_preliminary x y] get's the preliminary state of the character in
     * position [x] in the set [y]. *)
    val get_preliminary : int -> s -> c

    (** [set_final x y z] sets the final state of the character in
     * position [x] to [y] in the set [z]. *)
    val set_final : int -> c -> s -> unit

    (** [get_final x y] get's the final state of the character in
     * position [x] in the set [y]. *)
    val get_final : int -> s -> c

    (** [set_cost x y z] sets the cost of the transformation of the state in
     * position [x] to [y] in the two children of the node (HTU) that holds the
     * set [z]. *)
    val set_cost : int -> int -> s -> unit

    (** [get_cost x y] gets the cost of the transformation of the state in
     * position [x] in the two children of the node (HTU) that holds the
     * set [y]. *)
    val get_cost : int -> s -> int

    (** [set_weight x y z] sets the weight of the characater in position [x] to
    * [y] in the set [z]. *)
    val set_weight : int -> int -> s -> unit

    (** [get_weight x y z] gets the weight of the characater in position [x] 
    * in the set [z]. *)
    val get_weight : int -> s -> int

    (** [set_total x y z] sets the total cost of the transformation of the 
    * state in position [x] to [y] in the two children trees of the node (HTU) 
    * that holds the set [z]. *)
    val set_total : int -> int -> s -> unit

    (** [get_total x y] gets the total cost of the transformation of the state in
     * position [x] in the two children trees of the node (HTU) that holds the
     * set [y]. *)
    val get_total : int -> s -> int

    (** [set_char_code x y z] sets to [y] the character code of the character in
     * position [x] of the set [z]. *)
    val set_char_code : int -> int -> s -> unit

    (** [get_char_code x y] gets the character code of the character in
     * position [x] of the set [z]. *)
    val get_char_code : int -> s -> int

    (** [set_changed x y] sets to true the changed flag of the character in
     * position [x] of the set [y]. *)
    val set_changed : int -> s -> unit

    (** [set_changed x y] sets to false the changed flag of the character in
     * position [x] of the set [y]. *)
    val set_no_changed : int -> s -> unit

    (** [get_changed x y] gets the changed flag of the character in position [x]
     * of the character set [y]. *)
    val get_changed : int -> s -> bool

    (** [get_index_of_code x y] gets the index of the character with code [x] in
     * the character set [y] *)
    val get_index_of_code : int -> s -> int

    (** [downpass x y z] performs a downpass operation on the character set [x]
     * and [y] and stores the result in the character set [z]. The function
     * returns the total cost of the downpass step *)
    val downpass : s -> s -> s -> Cost_matrix.Two_D.m -> int

    (** [uppass x y z u] performs an uppass that sets the final states of the
    * character set [u] whith [x] being it's parent and [y] and [z] it's
    * children. The function assumes that the final states of [x] are already
    * set.  *)
    val uppass : s -> s -> s -> s -> unit

    (* [preliminary x] sets the preliminary states of [x] as the final. This
    * function should be called on the root of a tree before a full uppass.*)
    val preliminary_is_final : s -> unit

    (* [join x y] takes two disjoint sets of characters [x] and [y] and outputs
    * a fresh character set. *)
    val join : s -> s -> s

    (* [split x] takes a single character ste [x] and outputs two character sets
    * with half of it (or the floor and ceiling of the half of it) *)
    val split : s -> s * s

end
*)

(************************************************************************)
module Chromosome :
  sig
    type c
    type s = c array
    val cost : int
    val weight : int
    val total : int
    val code : int
    val changed : int
    external c_read_input : in_channel -> int -> int -> s
      = "read_input2_CAML"
    external c_calc_num_genomes_genes : in_channel -> int * int
      = "calc_num_genomes_genes2_CAML"
    external c_print_gene_matrix : s -> unit = "print_gene_matrix2_CAML"
    external c_print_gene_array : c -> unit = "print_gene_array2_CAML"
    external c_set_preliminary : int -> c -> s -> unit
      = "chromosome_CAML_set_preliminary"
    external c_get_preliminary : int -> s -> c
      = "chromosome_CAML_get_preliminary"
    external c_set_final : int -> c -> s -> unit
      = "chromosome_CAML_set_final" 
    external c_get_final : int -> s -> c = "chromosome_CAML_get_final"
    val set_preliminary : int -> c -> s -> unit
    val get_preliminary : int -> s -> c
    val set_final : int -> c -> s -> unit 
    val get_final : int -> s -> c
    external c_set_pre_to_fin : s ->  unit
      = "chromosome_CAML_set_pre_to_fin"
    val preliminary_is_final : s -> unit 
    external c_set_val : int -> int -> int -> s -> unit
      = "chromosome_CAML_set_val"
    external c_get_val : int -> int -> s -> int = "chromosome_CAML_get_val"
    val set_cost : int -> int -> s -> unit
    val set_all_cost : int -> s -> unit
    val get_cost : int -> s -> int
    val set_weight : int -> int -> s -> unit
    val get_weight : int -> s -> int
    val set_total : int -> int -> s -> unit
    val get_total : int -> s -> int
    val set_char_code : int -> int -> s -> unit
    val get_char_code : int -> s -> int
    val set_changed : int -> s -> unit
    val set_no_changed : int -> s -> unit
    val get_changed : int -> s -> bool
    external c_get_index_of_code : int -> s -> int
      = "chromosome_CAML_get_index_of_code"
    val get_index_of_code : int -> s -> int
    external c_downpass : s -> s -> s -> int -> int
      = "chromosome_CAML_downpass"
    val downpass : s -> s -> s -> int -> int
    external c_uppass : s -> s -> s -> s -> int -> unit
      = "chromosome_CAML_uppass"
    external c_set_all_weights : int -> s -> unit
      = "chromosome_CAML_set_all_weights"
    val set_all_weights : int -> s -> unit
    val remove_weights : s -> unit
    external c_join : s -> s -> s = "chromosome_CAML_join"
    val join : s -> s -> s
    external c_split : s -> s * s = "chromosome_CAML_split"
    val split : s -> s * s
    external c_create_empty_chrom : int -> int -> s = 
        "chromosome_CAML_create_empty"
    val create : 
        Parser.Hennig.Encoding.s array -> int array -> Parser.t array -> s
  end
  
(************************************************************************)
(* Replacing missing character stuff                                    *)

type types = | NonAdditive | Additive | Chromosome

type characters = 
    | NonAdditive_Char of NonAdditive.c | Additive_Char of Additive.c |
    Chromosome_Char of Chromosome.c

(* [same_type a b] is true iff a and b are of the same type. *)
val same_type : characters -> characters -> bool


(* [get_type_index t] produce the unique index value assigned to characters of
* type t. *)
val get_type_index : types -> int

(* [get_types c] produce a list of the character types to which c belongs *)
val get_types : characters -> types list

val create_types_array : 'a -> 'a array ref

(* [initial ()] creates a new set of characters array structure to be filled
* with the data for a given taxon. *)
val initial : unit -> characters array ref 
