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



(** A character set with all the information it needs for fast calculations *)
type 'a r = { 
    preliminary : 'a;    (** Characters of the downpass *)
    final       : 'a;    (** Characters of the uppass *)
    cost        : float; (** Total cost of the set of characters *)
    sum_cost    : float; (** Total cost up to this node in tree (inc'l) *)
    weight      : float; (** The weight of the set of characters *)
    time        : float;  (** The time length of the edge connecting this node
                            and its parent *)
}

(** Schemes for origin and loss costs *)
type origin_loss_cost = [
| `Flat of int * int ]                 (** Flat cost for each origin and loss *)

(** Methods for complex terminal alignments *)
type complex_term_method = [
| `Strictly_Same  (** Don't allow recombination between elements of a set:  sets
                      are only used for grouping *)

| `Any_Of of ((int * int * (int * int) list) * float)
(** Allow choosing of any of the elements to use for the median *)

(* | `Origin_Loss of origin_loss_cost *)
(** Allow simple origin and loss costs *)
]   



(** Sets of characters (or other items) *)
type 'a css = {
    sid   : int;                        (** ID (code) of this subset *)
    set   : 'a list;                    (** Contents *)
    smethod : complex_term_method; (** How to treat complex terminals *)
}

type cs =
    | Nonadd8 of NonaddCS8.t r  (** A set of non additive characters 
    with at most 8 states *)
    | Nonadd16 of NonaddCS16.t r  (** A set of non additive characters 
    with at most 16 states *)
    | Nonadd32 of NonaddCS32.t r  (** A set of non additive characters 
    with at most 32 states *)
    | Add of AddCS.t r                  (** A set of additive characters *)
    | Sank of SankCS.t r                (** A set of sankoff characters *)
    | Dynamic of DynamicCS.t r                  (** A set of dynamics *)
    | Kolmo of KolmoCS.t r              (** A set of kolmogorov characters *)
    | Set of cs css r                   (** A set of other characters *)
    | StaticMl of MlStaticCS.t r        (** A set of static ML characters *)

type exclude = ([`Excluded | `NotExcluded | `Either] * int * int * int) list

type node_data = { 
    characters : cs list;             (** The set of characters of the [ho]tu *)
    total_cost : float;     (** The total cost of the tree rooted by the node *)
    node_cost : float;
    taxon_code : int;       (** The code of the taxon associated with the 
                                node *)
    min_child_code : int;
    num_child_edges : int;
    num_height : int;
    num_otus : int;                (** How many OTUs are a child of this node *)
    exclude_sets : All_sets.Integers.t list;
    exclude_info : exclude;
    (** This allows us to count how many taxa from a set are children of the
        given node *)
}

(** Compares the final states information between two nodes. Follows the
* Pervasives.compare definition. *)
val compare_uppass : node_data -> node_data -> int

(** [compare_downpass a b] Compares the preliminary states information 
* between the nodes [a] and [b]. Follows the Pervasives.compare 
* definition. *)
val compare_downpass : node_data -> node_data -> int

val has_to_single : [ `Add | `Annchrom | `Breakinv | `Chrom | `Genome | `Kolmo
            | `Nonadd | `Sank | `Seq | `StaticMl ] list

(* Like [distance] but calculates it only for the characters that match the
* listed types. *)
val distance_of_type :   
    ?para:int option ->
      ?parb:int option ->
          [ `Add | `Annchrom | `Breakinv | `Chrom | `Genome | `Kolmo
            | `Nonadd | `Sank | `Seq | `StaticMl ] list -> 
                    node_data -> node_data -> float

(** [dist_2 delta a b c] calculates the cost of joining a vertex [n] between
* the vertices [a] and [b]. If the function finds that the cost is greated
* than [delta], it can return any value that is already greater than [delta]. 
* In that case, the actual cost of the join is unknown. *)
val dist_2 : float -> node_data -> node_data -> node_data -> float

(** [median_3 par x a b] calculates a median between [par], [a], and [b].
* This median does not need to be optimal, but should be as good as
* possible. *)
val median_3 : node_data -> node_data -> node_data -> node_data -> node_data

(** [set_node_cost c n] creates a fresh node [e] such that 
* [node_cost e = c], but otherwise undistinguishable from [n]. *)
val set_node_cost : float -> node_data -> node_data


(** [to_formatter_single] is a horrible function, horrible, horrible; it outputs in
* a horrible format. Check the Tags module for further information. *)
(* Compute the total rearrangement cost of the subtree rooted by node_data *)
val cmp_subtree_recost : node_data -> float

(**
*  The function take accumulated formatter (acc), data (d), node_data,
 *  node_id, parent_node_data as an option
*)
val to_formatter_single :
    ChromCS.IntSet.t * ChromCS.IntSet.t ->    
        Tags.attributes ->
            Data.d -> (node_data * node_data) -> int -> (node_data * node_data) option -> Tags.output


(**
 * [to_formatter_subtree (final, prel) b c d e f g h] creates Tags.output of the contents 
* of the node [d] with code [e], using as attributed [b], eliminating the
* chromosomal characters [final] and [prel], with children of code and node
* content [f] and [g], and with [(parent node * parent single assignment node)]
* [h]. The single assignment is a unique assignment per character, as opposed to
* the multiple valid assignments that are possible in a tree *)
val to_formatter_subtree : ChromCS.IntSet.t * ChromCS.IntSet.t -> Tags.attributes ->
  Data.d -> (node_data * node_data) -> int ->
    int *  node_data -> int *  node_data -> 
        (node_data * node_data) option -> Tags.output


(** [to_single (pre_ref_code, fi_ref_code) root p n] returns a node that contains per character a single state
 * which is closest to [p] among those available in [n]. Useful for tree length
 * verification.  Pre_ref_code and fi_ref_code contain active codes for chromosome characters. 
 * Inactive codes are eliminated from diagnosis. 
 * If p is the handle, alied_map is the root containing the aligned map between p
 * and n needed for chromosome stuff, else alied_map is assigned by p *)
val to_single : ChromCS.IntSet.t * ChromCS.IntSet.t -> node_data option ->
    node_data -> node_data -> node_data

(** [to_single_root n] is equivalent to [to_single n n]. *)
val to_single_root : ChromCS.IntSet.t * ChromCS.IntSet.t -> node_data -> node_data


(** [readjust ch1 ch2 par mine] returns a heuristically selected set of nodes which
* is located somewhere in between [ch1], [ch2], and [par], which are the two
* children and parent of [mine]. If no better node than [mine] can be found,
* then [mine] itself is returned.
*
* When only parsimony characters are being adjusted, the returned first element
* in the tupe holds [`Mine x], where [x] is the readjusted [mine]. When
* likelihood characters are inside the node, the first element in the tuple
* containes [`MineNChildren (x, y, z)] where [x] containes the readjusted
* [mine], [y] the readjusted [ch1] and [z] the readjusted [ch2]. For static
* homology characters, the only parameter readjusted in [y] and [z] is their
* branch length specification. *)
val readjust : All_sets.Integers.t option -> node_data -> node_data -> node_data ->
    node_data -> [`Mine of node_data | `MineNChildren of (node_data * node_data
    * node_data) ] * All_sets.Integers.t

val get_active_ref_code : node_data -> All_sets.Integers.t * All_sets.Integers.t *
    All_sets.Integers.t * All_sets.Integers.t 

module Standard : NodeSig.S with type e = exclude and type n = node_data

val merge : node_data -> node_data -> node_data

val empty : unit -> node_data

(** [total_cost_of_type t n] extracts the sum of the total cost of the node [n]
 * for all the characters of the type [t], as listed below *)
val total_cost_of_type :
  [> `Add | `Annchrom | `Breakinv | `Chrom | `Genome | `Nonadd | `Sank | `Seq ] ->
  node_data -> float

val to_string : node_data -> string
val print : node_data -> unit
val copy_chrom_map : node_data -> node_data -> node_data


