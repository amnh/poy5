(* POY 5.0 Alpha. A phylogenetic analysis program using Dynamic Homologies.   *)
(* Copyright (C) 2011 Andrés Varón, Lin Hong, Nicholas Lucaroni, Ward Wheeler,*)
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

type to_single = 
    [ `Add | `Annchrom | `Breakinv | `Chrom | `Genome | `Kolmo | `Nonadd 
    | `FixedStates | `Sank | `Seq | `StaticMl | `Ml ]
(** A character set with all the information it needs for fast calculations *)
type 'a r = { 
    preliminary : 'a;    (** Characters of the downpass *)
    final       : 'a;    (** Characters of the uppass *)
    cost        : float; (** Total cost of the set of characters *)
    sum_cost    : float; (** Total cost up to this node in tree (inc'l) *)
    weight      : float; (** The weight of the set of characters *)
    time        : float option * float option * float option;
                         (** time that this node connecst to it's children,
                          * first is the min_taxon_code, second is the other *)
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

IFDEF USE_LIKELIHOOD THEN
type ml_rep = MlStaticCS.t r
ELSE
type ml_rep = unit
END

type cs =
    | Nonadd8 of NonaddCS8.t r      (** non additive chars w/ <=8 states  *)
    | Nonadd16 of NonaddCS16.t r    (** non additive chars w/ <=16 states *)
    | Nonadd32 of NonaddCS32.t r    (** non additive chars w/ <=32 states *)
    | AddVec of AddCS.Vector.t r    (** additive characters   *)
    | AddGen of AddCS.General.t r   (** additive characters   *)
    | Sank of SankCS.t r            (** sankoff characters    *)
    | FixedStates of Fixed_states.t_w_seqtbl r (** fixed states characters *)
    | Dynamic of DynamicCS.t r      (** dynamics              *)
    | Kolmo of KolmoCS.t r          (** kolmogorov characters *)
    | Set of cs css r               (** other characters      *)
    | StaticMl of ml_rep            (** static ML characters  *)

type exclude = ([`Excluded | `NotExcluded | `Either] * int * int * int) list

type node_data = { 
    characters : cs list;   (** The set of characters of the [ho]tu *)
    total_cost : float;     (** The total cost of the tree rooted by the node *)
    node_cost : float;
    taxon_code : int;       (** The code of the taxon associated with the node *)
    min_child_code : int;
    num_child_edges : int;
    num_height : int;
    num_otus : int;         (** How many OTUs are a child of this node *)
    exclude_sets : All_sets.Integers.t list;
    exclude_info : exclude;
    cost_mode : [ `Likelihood | `NotLikelihood | `SumLikelihood  ]; 
}

(** Compares the final states information between two nodes. Follows the
* Pervasives.compare definition. *)
val compare_uppass : node_data -> node_data -> int

(** [compare_downpass a b] Compares the preliminary states information 
* between the nodes [a] and [b]. Follows the Pervasives.compare 
* definition. *)
val compare_downpass : node_data -> node_data -> int

val has_to_single : to_single list

val not_to_single : to_single list

val print_times : node_data -> unit

val total_cost : int option -> node_data -> float

(* Like [distance] but calculates it only for the characters that match the
* listed types. *)
val distance_of_type :   
  ?branches:(int,(int,float) Hashtbl.t) Hashtbl.t ->
    ?para:int option ->
      ?parb:int option ->
        to_single list -> 
          float -> node_data -> node_data -> float

(** [dist_2 delta a b c] calculates the cost of joining a vertex [n] between
* the vertices [a] and [b]. If the function finds that the cost is greated
* than [delta], it can return any value that is already greater than [delta]. 
* In that case, the actual cost of the join is unknown. *)
val dist_2 : float -> node_data -> node_data -> node_data -> float

(** [update_cost_only mine child1 child2] update mine with new sum_cost,
* calculated by new sum_cost of child1 or/and child2*)
val update_cost_only : node_data -> node_data -> node_data -> node_data

(** [final_states par x a b] calculates a median between [par], [a], and [b].
* This median does not need to be optimal, but should be as good as
* possible.
val final_states : node_data -> node_data -> node_data -> node_data -> node_data
*)

(** [set_node_cost c n] creates a fresh node [e] such that 
* [node_cost e = c], but otherwise undistinguishable from [n]. *)
val set_node_cost : float -> node_data -> node_data

(* Compute the total rearrangement cost of the subtree rooted by node_data *)
val cmp_subtree_recost : node_data -> float

(** [to_formatter_single] Check the Xml module for further information on the
    format.The function take accumulated formatter (acc), data (d), node_data,
    node_id, parent_node_data as an option *)
val to_formatter_single :
    Methods.diagnosis_report_type -> ChromCS.IntSet.t * ChromCS.IntSet.t ->    
        Xml.attributes ->
            Data.d -> (node_data * node_data) -> int -> (node_data * node_data)
            option -> Xml.xml 

(**
 * [to_formatter_subtree (final, prel) b c d e f g h] creates Xml.xml of the contents 
* of the node [d] with code [e], using as attributed [b], eliminating the
* chromosomal characters [final] and [prel], with children of code and node
* content [f] and [g], and with [(parent node * parent single assignment node)]
* [h]. The single assignment is a unique assignment per character, as opposed to
* the multiple valid assignments that are possible in a tree *)
val to_formatter_subtree : Methods.diagnosis_report_type -> ChromCS.IntSet.t * ChromCS.IntSet.t -> Xml.attributes ->
  Data.d -> (node_data * node_data) -> int ->
    int *  node_data -> int *  node_data -> 
        (node_data * node_data) option -> Xml.xml 

(** [to_single (pre_ref_code, fi_ref_code) root p n] returns a node that contains per character a single state
 * which is closest to [p] among those available in [n]. Useful for tree length
 * verification.  Pre_ref_code and fi_ref_code contain active codes for chromosome characters. 
 * Inactive codes are eliminated from diagnosis. 
 * If p is the handle, alied_map is the root containing the aligned map between p
 * and n needed for chromosome stuff, else alied_map is assigned by p *)
val to_single :
    ChromCS.IntSet.t * ChromCS.IntSet.t -> bool -> node_data option
        -> node_data -> node_data -> node_data

(** [to_single_root n] is equivalent to [to_single n n]. *)
val to_single_root : ChromCS.IntSet.t * ChromCS.IntSet.t -> node_data -> node_data

(** [apply_time r a b] applies the time [a] to the minimum child node and [b] to the
* other. If [r] is true, b is a root, and we add the edges together. *)
val apply_time : bool -> node_data -> node_data -> node_data

(** [edge_iterator gp p a b] is a function to iterate the branch lengths of
 * likelihood characters of parent [p], with children [a] and [b], and grand
 * parent [gp], which is optional. This function will return new node_data for 
 * [p], the branch lengths should be transfered to other directions via
 * uppass_heuristic
val edge_iterator : node_data option -> node_data -> node_data -> node_data ->
                    node_data
*)

(** [readjust mode to_adjust ch1 ch2 par mine] returns a heuristically selected
* set of nodes which is located somewhere in between [ch1], [ch2], and [par], 
* which are the two children and parent of [mine]. If no better node than [mine] 
* can be found, then [mine] itself is returned.
*
* Characters to adjust are help in [to_adjust], and characters that have changed
* are added into the second return value. [mode] is used in Dynamic Homology.
val readjust : [`ThreeD | `ApproxD ] -> All_sets.Integers.t option -> node_data ->
    node_data -> node_data -> node_data -> node_data * All_sets.Integers.t
*)

(* [median_w_times code node_1 node_2 times_1 times_2] 
 * uses the time data of nodes in different directions, [times_1] and [times_2]
 * to fill in the other directions with children [node_1] and [node_2].
*)
val median_w_times: 
    int option -> node_data option -> node_data -> node_data -> 
        float option list -> float option list -> float option list option ->
            node_data

(* [replace_parent_time node t] replace the third time in [node] with the time *)
val replace_parent_time : node_data -> float option list -> node_data

(* [median_of_child_branch child parent] creates a median between a child and
    a parent (from the persepective of an unrooted tree this makes sense), and
    uses the branch length from the child. This is intended to be used to update
    the parent with the branch length from a readjustment in three directions *)
val median_of_child_branch : int option -> node_data -> node_data -> node_data

(**[get_times_between(_tbl tbl) nd child_code]
 * helper functions to unify the distribution of times in three directions
 * --either from input tree or from the post_order downpass *)
val get_times_between : node_data -> int option -> float option list

(** [using_likelihood kind x] check if we are using likelihood characters; used
    for assertions in the allDirNode module to reduce errors and yet still keep
    code consistent with multi-character data **)
val using_likelihood : [`Static | `Dynamic | `Either] -> node_data -> bool

val get_active_ref_code : node_data -> All_sets.Integers.t * All_sets.Integers.t *
    All_sets.Integers.t * All_sets.Integers.t 

module Standard : NodeSig.S with type e = exclude
                            and type n = node_data
                            and type other_n = node_data

val merge : node_data -> node_data -> node_data

val empty : [ `Likelihood | `SumLikelihood | `NotLikelihood ] -> node_data

(** [total_cost_of_type t n] extracts the sum of the total cost of the node [n]
 * for all the characters of the type [t], as listed below *)
val total_cost_of_type : to_single -> node_data -> float

val extra_cost_from_root : node_data -> float

val get_cost_mode : node_data -> [ `Likelihood | `SumLikelihood | `NotLikelihood ]

(*val classify_data : bool -> node_data -> bool -> node_data -> int list option ->*)
(*        (float All_sets.FullTupleMap.t) * (float All_sets.IntegerMap.t) ->*)
(*            (float All_sets.FullTupleMap.t) * (float All_sets.IntegerMap.t)*)

val to_string : node_data -> string

val print : node_data -> unit

val copy_chrom_map : node_data -> node_data -> node_data

val median_counter : int ref

