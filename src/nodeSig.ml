(* POY 4.0 Beta. A phylogenetic analysis program using Dynamic Homologies.    *)
(* Copyright (C) 2007 Andrés Varón, Le Sy Vinh, Illya Bomash, Ward Wheeler,   *)
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

let () = SadmanOutput.register "NodeSig" "$Revision: 1616 $"

type direction = Without of int | Median of (int * int)

(** Specification of the signature for every node implementation. *)
module type S = sig

    (** {2 Types} *)

    (** The type of the nodes in the tree. *)
    type n 
    type other_n

    (** The calculation of the nodes in a tree start in a particular vertex. For
    * this reason we have to do two passes: the downpass and the uppass. The
    * downpass defined states are named the "preliminary" states in most
    * parsimony related literature. At the beginning, there is one vertex that
    * has to use the same preliminary and final states (usually the root), so
    * the fix_preliminary function assigns the preliminary states as final too. *)
    val fix_preliminary : n -> n

    (** [distance a b] calculates the distance between the nodes [a] and [b]. *)
    val distance : ?para:int option -> ?parb:int option -> float -> n -> n -> float

    (** [median par code prev a b] calculates a median between the nodes [a] and [b].
    * [prev] holds the previously calculated value of that median (could be used
    * for heuristic or speedup purposes). *)
    val median : ?branches:(int,(int,float) Hashtbl.t) Hashtbl.t ->
                        int option -> n option -> n -> n -> n

    (** [final_states granpa par cur a b] creates new node with [cur], parent
    * [par], children [a] and [b]. and parent of parent [grandpa]. This only
    * updates the final_assignment for the specified direction. *)
    val final_states : int option -> n -> n -> n -> n -> n

    (** [uppass_heuristic pardata partime cur a b] updates a node in three
     * directions in AllDirF. In One Direction and Standard we call the
     * final_states functions --which used to be `median_3`. This function will
     * fill in the data from the downpass in the other directions --passing the
     * data from [cur] to [a] and [b] as it's children, and recieving the data
     * from [par] as it's parent. This data, in likelihood, are the `times` or
     * branch lengths. The [pardata] and [partime] is split when we have no data
     * in the case of two disjoint nodes over a root. Time is taken from [partime]
     * and data from [pardata] *)
    val uppass_heuristic : n -> n option -> n -> n -> n -> n

    (** [apply_time given curr par] applies the time in par to cur --used for leafs **)
    val apply_time : bool -> n -> n -> n

    (** [extract_states a r b] extract parsimony states at node [b] toward [a],
     * exclude all characters not present in [r] if it exists, else all of them **)
    val extract_states : Alphabet.a -> Data.d -> int option -> int array option -> n ->
                            (float * int * MlModel.chars) list

    (** [get_times_between left right] time/branch lengths between two nodes **)
    val get_times_between : n -> n option -> (int array * (float option)) list

    (** [to_string n] produces a string representation of the node. This is used
    * for debugging purposes. There is no particular format requirement. *)
    val to_string : n -> string

    (** [total_cost par n] calculates the cost of the tree rooted by the node [n]. 
    * when the parent has code [par].
    * If [n] was created using [load_data] (see below), [tree_cost n = 0],
    * otherwise, if it was created using [median a b] (or subsequently the
    * [final_states] of such a median), then [total_cost n = (total_cost a) +.
    * (total_cost b) +. (distance a b)]. *)
    val total_cost : int option -> int list option -> n -> float

    val tree_cost : int option -> n -> float

    (** return the minimum prior in one of directions (if more then one) held
     * by the node on the median sequence *)
    val min_prior : n -> float

    (** [node_cost n] calculates the cost of generating [n]. If [n] was created
    * using [load_data], then [node_cost n = 0], otherwise if it was created
    * using [median a b] (or subsequently the [final_states] of such a median), then
    * [node_cost n = distance a b]. *)
    val node_cost : int option -> n -> float

    (** [update_leaf n] sets as the final state of [n] its preliminary state of
    * [n]. *)
    val update_leaf : n -> n

    (** [taxon_code n] gets the code assigned to the node [n]. Whenever a node
    * is created from [load_data], the [taxon_code n] should correspond to the
    * [Data.d] code assigned to the corresponding taxon. *)
    val taxon_code : n -> int

    (** [union_distance m n] calculates the distance between the unions
    * contained in [m] and [n]. 
    *
    * TODO: This function should be removed from this module *)
    val union_distance : n -> n -> float

    (** [is_collapsable m n] returns true if the edge defined by [m] and [n] can
    * be collapsed in a tree safely (for I/O purposes). *)
    val is_collapsable : [ `Dynamic | `Static | `Any ] -> n -> n -> bool

    (** [to_xml d chan n] outputs in an (obscure) XML format the contents of the
    * ndoe. This functionality is not currently used, and should be removed from
    * the source code. 
    *
    * TODO: This function should be removed from this module *)
    val to_xml : Data.d -> Pervasives.out_channel -> n -> unit

    (** [num_height n] calculates the maximum depth of the vertex from which
    * information is being requested. *)
    val num_height : int option -> n -> int

    (** [num_otus n] calculates the number of leaves contained in the tree
    * rooted by [n]. *)
    val num_otus : int option -> n -> int

    (** return all the sequence information contained in the node (including
    * the unions right now). This function is used in the auto sequence partition 
    * functions. Watch out! *)
    val get_sequences : 
	int option ->
        n -> 
            (int * Sequence.s * 
            Cost_matrix.Two_D.m * Cost_matrix.Three_D.m * Alphabet.a) list

    (** [get_dynamic_preliminary par n] returns a list with all the dynamic homology
    * characters contained in the node [n] in the preliminary field, in the
    * direction specified by the parental [par] . *)
    val get_dynamic_preliminary : int option -> n -> DynamicCS.t list

    (** [get_dynamic_adjusted n] is the same as [get_dynamic_preliminary] but
    * returning all the dynamic homology characters contained in the node in the
    * single assignment (if available). If not available, a failure is raised.
    * *)
    val get_dynamic_adjusted : int option -> n -> DynamicCS.t list

    (** [edge_distance m n] calculates the maximum distance between the vertices
    * [m] and [n]. *)
    val edge_distance : n -> n -> float

    (** [support_chars starting code n] returns the cost of the support 
    * characters. See the usage of the function in the Support module *)
    val support_chars : int -> int option -> n -> (float * int) list

    (** [load d] loads the data contained in [d], and creates a tuple 
    * [(a, b)], where [b] is the list of nodes contained in [d], and [d] is an
    * updated [d] that holds the information of the nodes created. The function
    * can optionally take two labled arguments, the list of codes of the [taxa]
    * to be loaded (then the list only contains those taxa), or the list of the
    * [codes] of the characters that are to be loaded.
    * MAKE SURE THE DATA IS CATEGORIZED by, (Data.categorize data) *)
    val load_data : 
        ?is_fixedstates:bool -> ?silent:bool -> (* ?taxa:(int list) -> ?codes:(int list) -> *)
            ?classify:bool -> Data.d -> Data.d * (n list)

    (** [n_chars ?acc n] gets the number of characters stored in [n], plus the
    * optional argument [acc]. If [acc] is not provided, the default value is
    * [0]. *)
    val n_chars : ?acc:int -> n -> int

    (** [prioritize n] assigns priorities to the characters in [n], so that
    * those with larger information contents, and those that are faster, are
    * processed first. *)
    val prioritize : n -> n

    (** [reprioritize a b], creates a fresh node that has the priorities
    * assigned to [a] but is otherwise undistinguishable from [b]. *)
    val reprioritize : n -> n -> n

    (* [f_codes c n] creates a fresh node [n] of which the codes [c] have been
    * filtered out. *)
    val f_codes : int list -> n -> n 

    (** [min_child_code par n] gets the smallest code assigned to any of it's
    * children when the parental of [n] is [par]. If [par = None] then it has to
    * be the case the [n] holds information in only one direction.*)
    val min_child_code : int option -> n -> int

    (** {3 Support Calculation} *)
    (** The type of the exclude information for the support value calculation. *)
    type e 

    (** The type of the non additive characters, for the support value
    * calculation. *)
    type nad8

    val new_characters :
        int ->
            nad8 All_sets.IntegerMap.t ->
                (All_sets.IntegerMap.key * int array All_sets.IntegerMap.t list) list ->
                    nad8 All_sets.IntegerMap.t
    val build_node :
        nad8 All_sets.IntegerMap.t ->
            All_sets.Integers.elt list -> n -> n

    val set_exclude_info : e -> n -> n
    val excludes_median : int option -> n -> n -> e 
    val has_excluded : e -> bool

    module T : sig
        val add_exclude : All_sets.Integers.t -> n -> n
        val remove_exclude : n -> n
    end

    (* This module handles the union operations of different characters *)

    module Union : sig
        (* The type of a union of nodes *)
        type u 
        (* [union x y z] creates a union of [x], [y], and [z]. Remember that 
        * this function is used on vertices, building the union of the 
        * union of the two children of a vertex ([y] and [z]), and the 
        * contents of the vertex itself [x]. *)
        val union : int option -> n -> u -> u -> u

        (* [union_preliminary x y z] adds the contents of the preliminary states
        * assigned to [z] to the union [y] with optional parent [x]. This 
        * function is only safe for static homology characters, all other 
        * characters will not modify the union in any form. Use the standard 
        * union for them. See also [union_final]. *)
        val union_preliminary : int option -> u -> n -> u

        (* [union_final x y z] is similar to [union_preliminary] excepting that
        * the final instead of the preliminary states of [y] are added to [z].
        * The same restrictions of [union_preliminary] apply. *)
        val union_final : int option -> u -> n -> u

        (* [leaf code par x] creates a union that only has the information in 
        * [x]. This is used in the leaves and the parental of the root of each 
        * union cluster. [code] is the optional assigned code for the leaf, 
        * otherwise the code of [x] is assigned. [par] is the optional 
        * code of the parental of [x]. If it is not assigned, it is assumed that
        * [x] has only one direction (this should be enforced by an
        * implementation). *)
        val leaf : int option -> int option -> n -> u

        (* [distance x y] calculates the distance between [x] and [y] *)
        val distance : u -> u -> float

        (* [saturation x] calculates the fraction of characters that have 
        * nonpolymorphic positions  *)
        val saturation : u -> float

        (* [distance_node x y] calculates the distance between the node 
        * [x] and the union [y]. *)
        val distance_node : int option -> n -> u -> float

        (** [compare a b] returns [0] iff [a] and [b] are the same, 
        * returns a positive integer iff [a > b], otherwise returns a negative
        * integer. This function complies with the [Pervasives.compare]
        * function. *)
        val compare : u -> u -> int

        (* [get_sequence code union] returns the character [code] of the
        * [union], provided that [code] refers to a sequence character. *)
        val get_sequence : int option -> int -> u -> SeqCS.union_element

    end

    (** [compare a b] compares the preliminary states of the nodes [a] and [b]
     * only. This function can not be used for the general comparison of the
     * complete contents of a node. *)
    val compare : n -> n -> int

    val for_support : 
        int -> (int * n) list -> int list -> int list -> n list

    val root_cost : n -> float
    val to_single : n option -> int option -> n -> int option -> n -> 
                        All_sets.Integers.t * All_sets.Integers.t -> n
                    
    val character_costs : int option -> n -> ([`NonAdd | `Add | `Sank] * int * float) list
    (* The following group of functions return the preliminary and final
    * assignments for a given vertex for each type, in the sets as grouped for
    * the node for faster computation *)
    val get_nonadd_8 : int option -> n -> (NonaddCS8.t * NonaddCS8.t) list
    val get_nonadd_16 : int option -> n -> (NonaddCS16.t * NonaddCS16.t) list
    val get_nonadd_32 : int option -> n -> (NonaddCS32.t * NonaddCS32.t) list
    val get_add : int option -> n -> (AddCS.t * AddCS.t) list
    val get_sank : int option -> n -> (SankCS.t * SankCS.t) list
    val get_dynamic : int option -> n -> (DynamicCS.t * DynamicCS.t) list
    val get_mlstatic : int option -> n -> (MlStaticCS.t * MlStaticCS.t) list

    (* Map all the internal codes of a node using the function *)
    val recode : (int -> int) -> n -> n

    val to_other : n -> other_n

    (** [edge_iterator mine child1 child2] -> mine,child1,child2
     * iterates the edges of child1 and child2 to find a minimum. main node is
     * modified, and in AllDirF the two children are also updated.
     *)
    val edge_iterator : n option -> n -> n -> n -> n

    (** [readjust dir mode to_adjust ch1 ch2 par mine] -> mine * to_adjust
     *
     * Adjusts [mine] with children [ch1],[ch2] and parent [par], using
     * direction [dir] --as grand parent, if needed. Returns an updated [mine]
     * with a set heuristically determined in the same spirit of [to_adjust].
     *
     * Dynamic uses [mode] to determine if approximations or three dimension
     * medians are calculated.
    **)
    val readjust : [`ThreeD of int option | `ApproxD of int option ] ->
        All_sets.Integers.t option -> n -> n-> n-> n -> n * All_sets.Integers.t

    val force : n -> n

end
