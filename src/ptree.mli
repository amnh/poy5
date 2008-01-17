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

(** ptree.mli - Interface to the phylogenetic tree/forest module. Phylogentic
* trees are unrooted trees with data associated to both interior and leaf nodes.
* Data could potentially be attached to the edges of these trees. Note that the
* direction of the edges is determined by the location of handles in the
* underlying unrooted tree. So when associating data with the edges, care must
* be taken to ensure that the data gets updated when handles move around in the
* underlying rooted tree. *)

type id = Tree.id
(** polymorphic variant to distinguish handles from ids. *)

type incremental = [
    | `Children of int
    | `No_Children of int
    | `HandleC of (int * int)
    | `HandleNC of (int * int)
]

type clade_cost = NoCost | Cost of float
(** The cost of pruning a clade. *)

type node = Tree.node
type edge = Tree.edge
(** The tree node and tree edge which are the node and edge types of the
* underlying unrooted tree. *)

type t_status = Tree.t_status
(** The traversal status. See Tree.t_status *)

type 'a root_node = ([ `Edge of (int * int) | `Single of int ] * 'a) option

type 'a root = {
    root_median : 'a root_node;
    component_cost : float;
    adjusted_component_cost : float;
}

type ('a, 'b) p_tree = {
  node_data : 'a All_sets.IntegerMap.t;
  edge_data : 'b Tree.EdgeMap.t;
  tree : Tree.u_tree;
  component_root : 'a root All_sets.IntegerMap.t;
  origin_cost : float;
}

type cost_type = [ `Adjusted | `Unadjusted ]

val get_cost : cost_type -> ('a, 'b) p_tree -> float

val set_origin_cost : float -> ('a, 'b) p_tree -> ('a, 'b) p_tree

type phylogeny = (Node.node_data, unit) p_tree
(** The phylogentic tree (p_tree). 'a and 'b are node and edge data types. *)

val empty : ('a, 'b) p_tree
(** The empty phylogenetic tree. *)

type ('a, 'b) break_fn = Tree.break_jxn -> ('a, 'b) p_tree ->
    (('a, 'b) p_tree * Tree.break_delta * float * int * 'a * incremental list )
(** type of function that breaks the tree at a given break jxn and returns the
* tree with the edge broken, change in the tree topology, any additional data
* that might be useful to the caller, the cost of breaking the tree at this jxn
* and finally the node that was the root of the clade. *)
    
type ('a, 'b) join_fn =   
    incremental list ->
    Tree.join_jxn ->
    Tree.join_jxn ->
    ('a, 'b) p_tree ->
    (('a, 'b) p_tree * Tree.join_delta)
(** type of function that performs a join of the clade and the tree. It takes as
* input the join jxns in the tree and the clade, the change in tree topology due
* to the pruning of the clade from the tree and returns the tree with the clade
* joined at the specified junction, the cost of the new tree and the topology
* change due to the join operation. Note that the tree_delta that is passed in
* is just the tree_delta that was returned as a result of the edge break
* operation. *)

type ('a, 'b) cost_fn =
    Tree.join_jxn -> Tree.join_jxn ->
    float ->
    'a ->
    ('a, 'b) p_tree ->
    clade_cost
(** type of the function that determines the cost of a join. It takes both the
* tree and clade junctions as arguments, the cost of pruning the clade and the
* node that was the root of the clade and returns the cost of joining the clade
* at this location. *)

type ('a, 'b) reroot_fn =
    bool ->
    Tree.edge ->
    ('a, 'b) p_tree ->
    ('a, 'b) p_tree * incremental list
(** type of the function to reroot the clade before joining it in a TBR
* operation. *)

type ('a, 'b) print_fn =
    string ->
    ('a, 'b) p_tree -> 
    unit
(** type of the function used to print debugging information etc during
* searches. *)

(** The tree operations needed to parameterize the Search operations. *)
module type Tree_Operations =
  sig
    type a 
    type b
    val break_fn : (a, b) break_fn
    val join_fn : (a, b) join_fn
    val cost_fn : (a, b) cost_fn
    val reroot_fn : (a, b) reroot_fn
    val string_of_node : a -> string
    val features : Methods.local_optimum -> (string * string) list -> (string * string) list
    val clear_internals : (a, b) p_tree -> (a, b) p_tree
    val downpass : (a, b) p_tree -> (a, b) p_tree
    val uppass : (a, b) p_tree -> (a, b) p_tree
    val incremental_uppass : (a, b) p_tree -> incremental list -> (a, b) p_tree
    val to_formatter :  
        Tags.attributes -> Data.d -> (a, b) p_tree -> Tags.output 
    (** [root_costs t] returns all possible roots in a tree (eg. every edge)
    * and the respective tree cost associated with it. *)
    val root_costs : (a, b) p_tree -> (Tree.edge * float) list
    val unadjust : (a, b) p_tree -> (a, b) p_tree
  end

class type ['a, 'b] wagner_edges_mgr = object
    method break_distance : float -> unit
    method next_edge : Tree.edge option
    method next_clade : 'a -> unit
    method new_delta : float -> unit
    method update_join : ('a, 'b) p_tree -> Tree.join_delta -> unit
    method clone : ('a, 'b) wagner_edges_mgr
    method exclude : Tree.edge list -> unit
end

  class type ['a, 'b] wagner_mgr =
    object
        method any_trees : bool
        method clone : ('a, 'b) wagner_mgr
        method init : 
            (('a, 'b) p_tree * float * clade_cost * 
                ('a, 'b) wagner_edges_mgr) list -> unit
        method next_tree : ('a, 'b) p_tree * float * ('a, 'b) wagner_edges_mgr
        method process :
            ('a, 'b) cost_fn ->
                float ->
                    'a ->
                        ('a, 'b) join_fn ->
                            Tree.join_jxn ->
                                Tree.join_jxn -> 
                                    ('a, 'b) p_tree -> 
                                        ('a, 'b) wagner_edges_mgr ->
                                            t_status
    method evaluate : unit
    method results : (('a, 'b) p_tree * float) list
    end


class type ['a, 'b] tabu_mgr = object

    method clone : ('a, 'b) tabu_mgr
    (** Function to create a deep-copy of the tabu. *)
        
    method break_edge : Tree.edge option
    method break_distance : float -> unit
    method features : (string * string) list -> (string * string) list
    method join_edge : [ `Left | `Right ] -> Tree.edge option
    method reroot_edge : [`Left | `Right] -> Tree.edge option
    method update_break :
        ('a, 'b) p_tree -> Tree.break_delta -> int -> int -> int -> unit
    method update_join : 
        ('a, 'b) p_tree -> Tree.join_delta -> unit

    method break_edges : Tree.edge list

end 

    (** The search manager object that parameterizes the search. This
    * allows us to incorporate various heuristics into the search easily. The
    * parameterized types are the types from the p_tree. *)
      class type ['a, 'b] search_mgr = object
          method features : (string * string) list -> (string * string) list

          method init : (('a, 'b) p_tree * float * clade_cost * ('a, 'b) tabu_mgr) list 
              -> unit
        (** Function to initialize the list of trees to be searched
         * and their individual costs and break_deltas associated with
         * them. *)
              
          method clone : ('a, 'b) search_mgr
        (** Function to create a fresh instance of a search mgr object
            with all data initialized to default values. *)
              
        (** [process cost_fn join_fn
         *       join_1_jxn join_2_jxn tree_delta 
         *       broken_tree -> Travesal Status *)
          method process : ('a, 'b) cost_fn -> float -> 'a ->
              ('a, 'b) join_fn -> incremental list ->
              Tree.join_jxn -> Tree.join_jxn -> 
              ('a, 'b) tabu_mgr -> ('a, 'b) p_tree -> t_status 
        (** This function decides whether to perform a join operation 
         * and add the tree to the queue of trees to be searched. *)
              
          method any_trees : bool 
        (** Function to return whether there are anymore trees to be
         * searched. *)
              
          method next_tree : (('a, 'b) p_tree * float * ('a, 'b) tabu_mgr)
        (** Function to return the next tree to be searched and its cost *)

          method results : (('a, 'b) p_tree * float * ('a, 'b) tabu_mgr) list
        (** Function to return the results of the search. *)

          method breakin : Tree.edge -> unit
          method should_repeat : bool
      end 

module type SEARCH = sig

    type a 
    type b
        
    val features : Methods.local_optimum -> (string
    * string) list -> (string * string) list

  val make_wagner_tree :
      ?sequence:(int list) ->
    (a, b) p_tree ->
    (a, b) wagner_mgr ->
    ((a, b) p_tree -> int -> (a, b) wagner_edges_mgr) ->
    (a, b) wagner_mgr

  val trees_considered : int ref

      type searcher =
          (a, b) search_mgr ->
          (a, b) search_mgr

      type search_step = 
          (a, b) p_tree ->
          (a, b) tabu_mgr ->
          searcher

      (** [spr_step ptree tabu search] performs one round of SPR searching on
          tree [ptree] using a given tabu manager and search manager. *)
      val spr_step : search_step

      (** [spr_simple search] takes each tree in search manager [search] and
          performs rounds of SPR until there is no further improvement *)
      val spr_simple : searcher
      
      (** [tbr_step ptree tabu search] performs one round of TBR searching on
          tree [ptree] using a given tabu manager and search manager. *)
      val tbr_step : search_step

      (** [tbr_single search] performs one step of TBR on each tree in the
          search manager *)
      val tbr_single : searcher
      val spr_single : searcher

      (** [tbr_simple search] takes each tree in search manager [search] and
          performs rounds of TBR until there is no further improvement *)
      val tbr_simple : searcher

      val tbr_join :
          (a, b) search_mgr ->
          ?updt:incremental list ->
          (a, b) tabu_mgr ->
          ?rerooted:bool ->
          (a, b) p_tree ->
          Tree.join_jxn -> a -> float -> Tree.t_status

      (** [alternate spr_searcher tbr_searcher search] takes each tree in search manager [search]
          and performs rounds of alternating SPR and TBR until there is no further
          improvement, using the [spr_searcher] and [tbr_searcher] passed as
          arugment. *)
      val alternate : searcher -> searcher -> searcher

      (** [repeat_until_no_more f s sm] iterates a search using whatever
       * searcher [s] is selected until no better tree is found, using the search
       * manager [sm], and creating a new tabu manager for each repetition using
       * the function [f]. This function is used in the [swap (around)] kind of
       * commands. *)
      val repeat_until_no_more : 
          ((a, b) p_tree -> (a, b) tabu_mgr) -> 
          searcher -> (a, b) search_mgr -> (a, b) search_mgr

      val get_trees_considered : unit -> int
      val reset_trees_considered : unit -> unit
      val uppass : (a, b) p_tree -> (a, b) p_tree
      val downpass : (a, b) p_tree -> (a, b) p_tree
      val diagnosis : (a, b) p_tree -> (a, b) p_tree

      (** [fuse_generations trees max_trees tree_weight tree_keep iterations
          process] runs a genetic algorithm-style search using tree fusing.  The function
          takes a list of trees to start with, the max number of trees and a method for
          keeping trees, a method for weighting trees, a number of iterations to perform,
          and a function to process new trees *)
      val fuse_generations :
          (a, b) p_tree list ->
          int ->
          ((a, b) p_tree -> float) ->
          Methods.fusing_keep_method ->
          int ->
          ((a, b) p_tree -> (a, b) p_tree
               list) ->
          (int * int) ->
          (a, b) p_tree list

      val search_local_next_best : (search_step * string) -> searcher
      val search : (search_step * string) -> searcher

        val convert_to :
          string Parser.Tree.t list ->
          Data.d * a list -> (a, b) p_tree

        val build_trees: Tree.u_tree -> 
            (int -> string) -> 
                (int -> int -> bool) -> 
                    string -> string Parser.Tree.t list
        val build_tree : Tree.u_tree -> 
            (int -> string) -> 
                (int -> int -> bool) -> string -> 
                    string Parser.Tree.t

        val never_collapse :  (a, b) p_tree -> int -> int -> bool

        val collapse_as_needed : (a, b) p_tree -> int -> int -> bool

        val get_unique : (a, b) p_tree list -> (a, b) p_tree list 

        val build_tree_with_names :
        bool -> (a, b) p_tree -> Data.d -> string Parser.Tree.t

        val build_tree_with_names_n_costs :
        bool -> (a, b) p_tree -> Data.d -> string -> string Parser.Tree.t
        val build_forest :
            bool -> (a, b) p_tree ->
          Data.d -> string -> string Parser.Tree.t list
        val build_forest_as_tree :
            bool -> (a, b) p_tree -> Data.d -> string -> string Parser.Tree.t

        val build_forest_with_names :
            bool -> (a, b) p_tree -> Data.d -> string Parser.Tree.t list
        val build_forest_with_names_n_costs :
            bool -> (a, b) p_tree -> Data.d -> string -> string Parser.Tree.t list


        val to_xml : 
            Pervasives.out_channel -> (a, b) p_tree -> Data.d -> unit
        val disp_trees : 
            string -> 
                (a, b) p_tree -> 
                    ((a, b) p_tree -> int -> string) -> 
                            string -> unit

    end
    
val set_avail_start : ('a, 'b) p_tree -> ('a, 'b) p_tree
val int_of_id : Tree.id -> int
val get_id : Tree.node -> int
val is_handle : int -> ('a, 'b) p_tree -> bool
val is_edge : Tree.edge -> ('a, 'b) p_tree -> bool
val get_node_id : int -> ('a, 'b) p_tree -> Tree.id
val get_handle_id : int -> ('a, 'b) p_tree -> Tree.id
val get_handles : ('a, 'b) p_tree -> All_sets.Integers.t

(** [components forest] returns the number of components in the forest *)
val components :
    ('a, 'b) p_tree -> int

val get_node : int -> ('a, 'b) p_tree -> Tree.node
val get_edge : int * int -> ('a, 'b) p_tree -> Tree.edge
val get_parent : int -> ('a, 'b) p_tree -> int
val other_two_nbrs : int -> Tree.node -> int * int
val get_nodes : ('a, 'b) p_tree -> Tree.node list
val get_pre_order_edges : int -> ('a, 'b) p_tree -> edge list
val get_edges_tree : ('a, 'b) p_tree -> edge list
val get_node_ids : ('a, 'b) p_tree -> int list
val add_node_data :
  int -> 'a -> ('a, 'b) p_tree -> ('a, 'b) p_tree
val add_edge_data :
  Tree.EdgeMap.key -> 'a -> ('b, 'a) p_tree -> ('b, 'a) p_tree
val remove_node_data :
  int -> ('a, 'b) p_tree -> ('a, 'b) p_tree
val remove_edge_data : Tree.EdgeMap.key -> ('a, 'b) p_tree -> ('a, 'b) p_tree
val get_node_data : int -> ('a, 'b) p_tree -> 'a
val get_edge_data : Tree.EdgeMap.key -> ('a, 'b) p_tree -> 'b
val handle_of : id -> ('a, 'b) p_tree -> id

(** [get_leaves h tree] returns a list of the leaves in [tree] with handle [h] *)
val get_leaves : ?init:('a list) -> int -> ('a, 'b) p_tree -> 'a list

(** [get_all_leaves tree] returns a list of all the leaves in [tree] *)
val get_all_leaves : ('a, 'b) p_tree -> 'a list

(** [get_leaves h tree] returns a list of the leaves in [tree] with handle [h] *)
val get_leaves_ids : ?init:(int list) -> int -> ('a, 'b) p_tree -> int list

(** [get_all_leaves tree] returns a list of all the leaves in [tree] *)
val get_all_leaves_ids : ('a, 'b) p_tree -> int list

val move_handle : int -> ('a, 'b) p_tree -> ('a, 'b) p_tree * int list
val pre_order_node_visit :
  (int option -> int -> 'a -> Tree.t_status * 'a) ->
  int -> ('b, 'c) p_tree -> 'a -> 'a
val post_order_node_visit :
  (int option -> int -> 'a -> Tree.t_status * 'a) ->
  int -> ('b, 'c) p_tree -> 'a -> 'a

(* [post_order_downpass_style l i c t] applies the functions [l] and [i] to the
* leafs (or singles), and [i] to the interior nodes, following a downpass style
* of pass. The function requires the root being assigned to the requested
* component [c] from the tree [t]. The arguments for [l] and [i] are the code of
* the parent (optional), and the code of the vertex being applied (optional in
* the [i] function, as the root may have code [None]). *)
val post_order_downpass_style : (int option -> int -> 'a) ->
    (int option -> int option -> 'a -> 'a -> 'a) ->
        int -> ('b, 'c) p_tree -> 'a

val post_order_node_with_edge_visit :
    (int -> int -> 'a -> 'a) -> (int -> int -> 'a -> 'a -> 'a) ->
        edge -> ('b, 'c) p_tree -> 'a -> 'a * 'a

val pre_order_edge_visit :
(*   ?max_d:int ->  *)(Tree.edge -> 'a -> Tree.t_status * 'a) ->
  int -> ('b, 'c) p_tree -> 'a -> 'a
val print_tree : int -> ('a, 'b) p_tree -> unit
val print_forest : ('a, 'b) p_tree -> unit
val make_disjoint_tree : 'a All_sets.IntegerMap.t -> ('a, 'b) p_tree
val build_consensus :
  int ->
  int ->
  int ->
  ('a, 'b) p_tree list ->
  int * Gen_rtree.g_rtree * (int, int) Hashtbl.t *
  (int, int) Hashtbl.t


module Search (Node : NodeSig.S) 
    (Edge : Edge.EdgeSig with type n = Node.n) (Tree_Ops : Tree_Operations with type a = Node.n with
    type b = Edge.e) : SEARCH with type a = Tree_Ops.a with type b = Tree_Ops.b

(** Module to fingerprint trees and compare them *)
module Fingerprint : sig
    type t
    val fingerprint : ('a, 'b) p_tree -> t
    val compare : t -> t -> int
    val to_string : t -> string
    val empty : t
    end

(** [assign_root_to_connected_component h n ec ac t] assigns a root to a
* connected component of the tree with handle [h], assigning the node specified
* in the edge or single by [n], with estimated component cost [ec] and adjusted
* cost [ac] to the tree [t]. *)
val assign_root_to_connected_component : 
    int -> 'a root_node -> float -> float option ->
    ('a, 'b) p_tree -> ('a, 'b) p_tree

(** [get_cost c t] gets the adjusted or unadjusted cost (as specified by [c]),
* of the tree [t]. *)
val get_cost : cost_type -> ('a, 'b) p_tree -> float

(** [get_handle_cost c t h] gets the adjusted or unadjusted cost (as specified by [c]),
* of the connected component with handle [h] of the tree [t]. *)
val get_handle_cost : cost_type -> ('a, 'b) p_tree -> int -> float

(** the same as assign_root_to_connected_component with different argument
* ordering *)
val set_component_cost : 
    float -> float option -> 'a root_node -> 
        int -> ('a, 'b) p_tree -> ('a, 'b) p_tree

val wipe_costs : ('a, 'b) p_tree -> ('a, 'b) p_tree 

val move_cost_n_root : int -> int -> ('a, 'b) p_tree -> ('a, 'b) p_tree

val get_component_root : int -> ('a, 'b) p_tree -> 'a root 

val get_roots : ('a, 'b) p_tree -> 'a root list

val change_component_root :
  All_sets.IntegerMap.key -> 'a root -> ('a, 'b) p_tree -> ('a, 'b) p_tree


val remove_root_of_component : int -> ('a, 'b) p_tree -> ('a, 'b) p_tree

val fix_handle_neighbor : int -> int -> ('a, 'b) p_tree -> ('a, 'b) p_tree

val break_median_edge :
    (int option -> int option -> 'a -> 'a -> 'b) -> int -> ('a, 'c) p_tree -> Tree.edge

(** [jxn_of_handle ptree handle] returns the junction corresponding to a
    handle's position in the tree (forest) *)
val jxn_of_handle :
    ('a, 'b) p_tree -> int -> Tree.join_jxn

val get_parent_or_root_data :
    id -> ('a, 'b) p_tree -> 'a

val choose_leaf : ('a, 'b) p_tree -> int

val cannonize :  string Parser.Tree.t -> string * string Parser.Tree.t

val consensus :
    (('a, 'b) p_tree -> int -> int -> bool) -> (int -> string) -> int -> 
    ('a, 'b) p_tree list -> int -> string Parser.Tree.t

val add_tree_to_counters :
    (int -> int -> bool) -> int Tree.CladeFPMap.t ->
    Tree.u_tree -> int Tree.CladeFPMap.t

val add_consensus_to_counters :
    int Tree.CladeFPMap.t ->
    ((int -> int -> bool) * Tree.u_tree) list -> int Tree.CladeFPMap.t

val supports :
    (int -> string) -> int -> float -> Tree.u_tree -> int Tree.CladeFPMap.t ->
    string Parser.Tree.t

val extract_bremer :  (All_sets.Integers.elt -> string) ->
  int Tree.CladeFPMap.t -> string Parser.Tree.t

(** [bremer to_string cost t conversion file] calculates a bremer support tree (for
* printing purposes) of the tree [t] (which must be properly rooted) which has
* cost [cost] (the cost must be an integer, as bremer is only used in parsimony
* and this improves function reusage inside the library), using the
* function [to_string] to convert leaves in the tree to strings (taxon names),
* and the sexpr of clades and costs contained in the file [file], that is, a
* file containing trees and cost of those trees [(a, * b)], where [b] is the cost of the tree [a]. The function [conversion] takes
* care of converting each touple in [sets] to a cost and the corresponding set
* of clades for the actual bremer calculations. This is done to reduce memory
* consumption.
*
* The resulting tree assigns to each branch the minimum cost found for a tree
* not containing the child clade of the branch within [sets]. *)
val bremer :
    (int -> string) -> int -> Tree.u_tree -> 
              ((string Parser.Tree.t * string) -> (int * Tree.CladeFP.CladeSet.t)) ->
          Parser.filename -> 
              string Parser.Tree.t
val preprocessed_consensus :
  (All_sets.Integers.elt -> string) ->
  int -> int -> int Tree.CladeFPMap.t -> string Parser.Tree.t

