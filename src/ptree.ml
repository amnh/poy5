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

let () = SadmanOutput.register "Ptree" "$Revision: 2329 $"

let ndebug = false
let ndebug_break_delta = false
let ndebug_traject = false
let ndebug_traject_spr = false || ndebug_traject
let ndebug_traject_tbr = false || ndebug_traject
let ndebug_traject_summary = false      (* summary of the trajectory *)
let ndebug_jxn_of_handle = false
let debug_wagner_traject = false
let odebug = Status.user_message Status.Information

type id = Tree.id
type incremental = [
    | `Children of int
    | `No_Children of int
    | `HandleC of (int * int)
    | `HandleNC of (int * int)
]

type clade_cost = NoCost | Cost of float
type node = Tree.node
type edge = Tree.edge
type t_status = Tree.t_status

type 'a root_node = ([ `Edge of (int * int) | `Single of int ] * 'a) option

type 'a root = {
    root_median : 'a root_node;
    component_cost : float;
    adjusted_component_cost : float;
}

type ('a, 'b) p_tree = { 
    node_data : 'a All_sets.IntegerMap.t ;
    edge_data : 'b Tree.EdgeMap.t ;
    tree : Tree.u_tree;
    component_root : 'a root All_sets.IntegerMap.t;
    origin_cost : float;
}


type phylogeny = (Node.node_data, unit) p_tree

    let rec cannonize a = 
        match a with
        | Parser.Tree.Leaf x -> x, a
        | Parser.Tree.Node (lst, _) ->
                let res = List.rev_map cannonize lst in
                match List.sort (fun ((x : string), _) ((y : string), _) ->
                    Pervasives.compare x y ) res  with
                    | ((h, _) :: _) as lst ->
                            let _, res = List.split lst in
                            h, (Parser.Tree.Node (res, h))
                    | [] -> failwith "Ptree.cannonize"

type cost_type = [ `Adjusted | `Unadjusted ]

let get_cost clas ptree =
    let get_cost = 
        match clas with
        | `Adjusted -> fun x -> x.adjusted_component_cost
        | `Unadjusted -> fun x -> x.component_cost
    in
    if ptree.origin_cost = infinity
    then
          let adder = fun _ v acc -> (get_cost v) +. acc in
          All_sets.IntegerMap.fold adder ptree.component_root
              0.
    else
          let adder = fun _ v acc -> (get_cost v) +. acc +. ptree.origin_cost in
          All_sets.IntegerMap.fold adder ptree.component_root
              (-. ptree.origin_cost)

let set_origin_cost cost ptree =
    { ptree with origin_cost = cost }

let remove_root_of_component node ptree = 
    { ptree with component_root = 
        All_sets.IntegerMap.remove node ptree.component_root }

let empty = { 
    node_data = All_sets.IntegerMap.empty ;
    edge_data = Tree.EdgeMap.empty ;
    tree = Tree.empty ();
    component_root = All_sets.IntegerMap.empty;
    origin_cost = 0.;
}

(** [set_avail_start ptree id] tells [ptree] to start creating new nodes with
    ID [id] *)
let set_avail_start ptree =
    { ptree with tree = Tree.set_avail_start ptree.tree }

type ('a, 'b) break_fn = Tree.break_jxn -> ('a, 'b) p_tree ->
    (('a, 'b) p_tree * Tree.break_delta * float * int * 'a * incremental list)

type ('a, 'b) join_fn =   
    incremental list ->
    Tree.join_jxn ->
    Tree.join_jxn ->
    ('a, 'b) p_tree ->
    ('a, 'b) p_tree * Tree.join_delta

type ('a, 'b) cost_fn =
    Tree.join_jxn -> Tree.join_jxn ->
    float ->
    'a ->
    ('a, 'b) p_tree ->
    clade_cost
    
type ('a, 'b) reroot_fn =
    bool ->
    Tree.edge ->
    ('a, 'b) p_tree ->
    ('a, 'b) p_tree * incremental list

type ('a, 'b) print_fn =
    string ->
    ('a, 'b) p_tree -> 
    unit

module type Tree_Operations = 
    sig
        type a
        type b
        val break_fn : (a, b) break_fn
        val join_fn : (a, b) join_fn 
        val cost_fn : (a, b) cost_fn
        val reroot_fn : (a, b) reroot_fn
        val string_of_node : a -> string
        val features : 
            Methods.local_optimum -> (string * string) list -> 
                (string * string) list
        val clear_internals : (a, b) p_tree -> (a, b) p_tree
        val downpass : (a, b) p_tree -> (a, b) p_tree
        val uppass : (a, b) p_tree -> (a, b) p_tree
        val incremental_uppass : 
            (a, b) p_tree -> incremental list -> (a, b) p_tree


        val to_formatter :  
            Tags.attributes -> Data.d -> (a, b) p_tree ->
            Tags.output

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

(** The tabu manager object that parameterizes the both the edges and the
* order in which the edges are visited by the searches. The invariant to be
* maintained by any tabu_mgr is that the edges in the tabu are always present in
* the tree. The directionality is flexible though. i.e. If the tree has edge(e1,
* e2) the tabu could have edge(e2, e1). *)
class type ['a, 'b] tabu_mgr = object

    method break_edge : Tree.edge option

    method break_distance : float -> unit

    method join_edge : [`Left | `Right] -> Tree.edge option

    method reroot_edge : [`Left | `Right] -> Tree.edge option

    method clone : ('a, 'b) tabu_mgr
    (** Function to create a deep-copy of the tabu. *)
    
    method update_break :
        ('a, 'b) p_tree -> Tree.break_delta -> int -> int -> int
        -> unit
    (** Function to update the tabu after a break operation. This function
     * should ensure the invariant that edges in the tabu and the edges in the tree are
     * in sync. Note that directionality could be reversed. *)

    method update_join : 
        ('a, 'b) p_tree -> Tree.join_delta -> unit
        
    method features : (string * string) list -> (string * string) list
    (** What does this method do? *)

    method break_edges : Tree.edge list
    
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

(** The search manager object that parameterizes the search. This
* allows us to incorporate various heuristics into the search easily. The
* parameterized types are the types from the p_tree. *)
  class type ['a, 'b] search_mgr = object
      method features : (string * string) list -> (string * string) list

      method init : 
          (('a, 'b) p_tree * float * clade_cost * ('a, 'b) tabu_mgr) list -> unit
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

      (** [alternate_spr_tbr search] takes each tree in search manager [search]
          and performs rounds of alternating SPR and TBR until there is no further
          improvement *)
      val alternate : searcher -> searcher -> searcher

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
    
(** The internal debug flag. *)
let debug = false
let debug_verify_costfn = false
let debug_verify_costfn_printnodes = false
let debug_verify_costfn_except_wagner = false
let debug_cost_print = true 
    
(** [int_of_id id]
    @param id is the id with associated type info `Node or `Handle.
    @return the int form of the id. See {!Tree.int_of_id} *)
let int_of_id = Tree.int_of_id

let handle_of id {tree=tree} = Tree.handle_of id tree

(** [get_id node]
    @return the id associated with the node. *)
let get_id = Tree.get_id 

(** [is_handle hnd ptree]
    @param hnd the handle id being verified.
    @param ptree the tree in which the id is either a handle or not. 
    @return true if the id corresponds to a handle. false, otherwise. *)
let is_handle hnd ptree =
    (Tree.is_handle hnd ptree.tree)

(** [is_edge (x,y) ptree]
    @param (x,y) the tuple that corresponds to an edge or not.
    @param ptree the tree in which the tuple, corresponds to an edge or not.
    @return true if the tuple corresponds to an edge, false if otherwise. *)
let is_edge e ptree = 
    (Tree.is_edge e ptree.tree)

(** [get_node_id id ptree]
    @param id the int whose id form is desired. id form tacks on a `Node
              to the integer.
    @param ptree the tree in which the int represents a node. 
    @raise Tree.Invalid_Node_Id
    @return The id form or throws an exception if the node is invalid. *)
let get_node_id id ptree = 
    (Tree.get_node_id id ptree.tree)

(** [get_handle_id id ptree]
    @param id the int whose id form is desired. id form tacks on a `Handle
              to the integer.
    @param ptree the tree in which the int represents a handle. 
    @raise Tree.Invalid_Handle_Id
    @return the handle form of the id. *)
let get_handle_id id ptree = 
    (Tree.get_handle_id id ptree.tree)
   
(** [get_handles ptree]
    @return the handle set of the tree. *)
let get_handles ptree = 
    (Tree.get_handles ptree.tree)

(** [components forest] returns the number of components in the forest *)
let components forest =
    All_sets.Integers.cardinal (get_handles forest)


(** [get_node id ptree]
    @param id the int-id of the the node.
    @param ptree the tree from which the node of the given id is being
                retrieved.
    @return the node associated with the given id.
    @raise Tree.Invalid_Node_Id when the id is invalid. *)
let get_node id ptree =
    (Tree.get_node id ptree.tree)

(** [get_edge (x, y) ptree]
    @param (x, y) the pair that could correspond to an edge in the tree.
    @param ptree the tree from which the edge is being retrieved.
    @return the edge corresponding to the pair.
    @raise Tree.Invalid_Edge if the pair doesnt correspond to an edge. *)
let get_edge (x, y) ptree = 
    (Tree.get_edge (x, y) ptree.tree)

(** [get_parent id ptree]
    @param id the id of the node whose parent is desired. Parent is the node
    that would be visited immediately before the present node in a
    pre-order-traversal from the handle of the component to which the node
    belongs. An exception is thrown when the id corresponds to a handle.
    @param ptree the tree in which the parent of node with id=id is determined.
    @return the parent of the node if it exists, otherwise an exception is
    thrown. *)
let get_parent id ptree = 
    (Tree.get_parent id ptree.tree)

(** [other_two_nbrs nbr node]
    @param nbr a nbr of the node
    @param node the interior node whose other two nbrs are desired.
    @return the other two nbrs of the node. 
    @raise Invalid_argument if the node is not an interior node. *)
let other_two_nbrs = Tree.other_two_nbrs
    
(** [get_nodes ptree]
    @return the list of all the nodes of the tree. *)
let get_nodes ptree = 
    (Tree.get_nodes ptree.tree)
   
(** [get_pre_order_edges hs ptree]
    @return the list of edges when visiting the tree in a pre-order
            traversal starting at the handle hs. *)
let get_pre_order_edges hs ptree = 
    (Tree.get_pre_order_edges hs ptree.tree)
    
let get_edges_tree ptree = 
    (Tree.get_edges_tree ptree.tree)
    
(** [get_node_ids ptree]
    @return the list of all node ids of the tree. *)
let get_node_ids ptree = 
    (Tree.get_node_ids ptree.tree)

(** [add_node_data id data ptree]
    @param id to which node data is added.
    @param data being added to the node
    @param ptree to which node data is being added. 
    @return new ptree with the data added to node with id=id.
            Any old data is silently overwritten. *)
let add_node_data id data ptree = 
    let new_node_data = (All_sets.IntegerMap.add id data ptree.node_data) in
        { ptree with node_data = new_node_data }

(** [add_edge_data edge data ptree]
    @param edge to which data is being added.
    @param data being added to the edge.
    @param ptree to which the edge data is being added.
    @return new ptree with the edge data added. Any old data
            will be silently overwritten. *)
let add_edge_data edge data ptree = 
    let edge = Tree.normalize_edge edge ptree.tree in 
    let new_edge_data = (Tree.EdgeMap.add edge data ptree.edge_data) in
        { ptree with edge_data = new_edge_data }

(** [remove_node_data id ptree]
    @return a new ptree with data associated with node_id=id erased. *)
let remove_node_data id ptree =
    let new_node_data = (All_sets.IntegerMap.remove id ptree.node_data) in
        { ptree with node_data = new_node_data }

(** [remove_edge_data edge ptree]
    @return a new ptree with data associated with edge removed. *)
let remove_edge_data edge ptree = 
    let edge = Tree.normalize_edge edge ptree.tree in
    let new_edge_data = (Tree.EdgeMap.remove edge ptree.edge_data) in
        { ptree with edge_data = new_edge_data }

(** [get_node_data id ptree]
    @return data associated with the node. *)
let get_node_data id ptree = 
    (All_sets.IntegerMap.find id ptree.node_data)

(** [get_parent_or_root_data id ptree] returns the data of the node's "literal"
    parent if the node is in the tree, or returns the root data if the root is
    the node's real edge *)
let get_parent_or_root_data id ptree =
    let get_root parent =
        match
            (All_sets.IntegerMap.find parent
                 ptree.component_root).root_median
        with
        | None -> failwith "get_parent_or_root_data"
        | Some (_, data) -> data in
    if is_handle id ptree
    then get_root id
    else
        let parent = get_parent id ptree in
        if is_handle parent ptree
        then
            match get_node parent ptree with
            | Tree.Leaf (_, par) when par = id ->
                  get_root parent
            | Tree.Interior (_, par, _, _) when par = id ->
                  get_root parent
            | _ -> get_node_data parent ptree
        else get_node_data parent ptree
    
(** [get_edge_data edge ptree]
    @return data associated with the edge. *)
let get_edge_data edge ptree = 
    try (Tree.EdgeMap.find edge ptree.edge_data) with
    | Not_found -> 
            let Tree.Edge (a, b) = edge in
            Tree.EdgeMap.find (Tree.Edge (b, a)) ptree.edge_data
    
let move_cost_n_root hid id ptree =
    let comp_cost = ptree.component_root in
    try 
        let x = All_sets.IntegerMap.find hid comp_cost in
        let res = All_sets.IntegerMap.remove hid comp_cost in
        let res = 
            All_sets.IntegerMap.add id { x with root_median = None } 
            res
        in
        { ptree with component_root = res }
    with
    | Not_found -> ptree

(** [move_handle id ptree]
    @return new ptree with the current id as a handle. THE DATA
    ASSOCIATED WITH edges Edge(e2, e1) IS NOT UPDATED. *)
let move_handle id ptree = 
    let hid, _ = Tree.get_path_to_handle id ptree.tree in
    let bt, path = (Tree.move_handle id ptree.tree) in
    move_cost_n_root hid id { ptree with tree = bt }, path

let fix_handle_neighbor h n ptree =
    { ptree with tree = Tree.fix_handle_neighbor h n ptree.tree }

(* 
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    How do we specify that edges have no data and no need to update
    edge data ? 
    let rev_edge pt edge = 
        let Tree.Edge(e1, e2) = edge in
        let rev = Tree.Edge(e2, e1) in
        let ed = (get_edge_data edge pt) in
        let pt = (remove_edge_data edge pt) in
        let pt = (add_edge_data rev ed pt) in
            pt
    in
        (List.fold_left rev_edge ptree edges)
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
*)

(** [pre_order_node_visit f id ptree acc]
    Function to perform a pre order node visit on the tree.
    See {!Tree.pre_order_node_visit} *)
let pre_order_node_visit f id ptree acc = 
    Tree.pre_order_node_visit f id ptree.tree acc 

(** [post_order_node_visit f id ptree acc]
    Function to perform a post order node visit on the tree.
    See {!Tree.post_order_node_visit} *)
let post_order_node_visit f id ptree acc = 
    Tree.post_order_node_visit f id ptree.tree acc 

(** [post_order_node_with_edge_visit]
    @return What does this function do? *)
let post_order_node_with_edge_visit f g e ptree acc = 
    Tree.post_order_node_with_edge_visit f g e ptree.tree acc
    
(** [pre_order_edge_visit f id ptree acc]
    Function to perform a pre order edge visit on the tree.
    See {!Tree.pre_order_edge_visit} *)
let pre_order_edge_visit f id ptree acc = 
    Tree.pre_order_edge_visit f id ptree.tree acc
    
(** [print_tree hnd ptree]
    Function to print a tree rooted at the given handle *)
let print_tree hnd ptree = 
    Tree.print_tree hnd ptree.tree
    
(** [print_forest ptree]
    Function to print the entire forest. *)
let print_forest ptree = 
    Tree.print_forest ptree.tree



    
(** [make_disjoint_tree n leaf_data_map]
    @param n the number of leaves.
    @param leaf_data_map a map from {1, 2, ... n} to the leaf data.
    @return the disjointed tree with n nodes. 0 edges and data
            associated with the nodes. *)
let make_disjoint_tree leaf_data_map = 
    let f k d acc = k :: acc in
    let nodes = (All_sets.IntegerMap.fold f leaf_data_map []) in
    let bt = Tree.make_disjoint_tree nodes in
        { empty with tree = bt ; node_data = leaf_data_map }

(** [build_consensus out_grp rt_id ptrees]
    @param maj the majority desired. 
    @param out_grp the leaf node at which the ptrees are rooted to get 
        the rooted topologies used for building consensus trees. 
    @param rt_id the id of the added root node. 
    @param ptrees the list of ptrees whose consensus is desired. 
    @return the consensus tree and an assoc-list of the node-ids and the node
        data. *)
let build_consensus maj out_grp rt_id ptrees = 
    (* extract the topologies of the ptrees. *)
    let utrees = List.map (fun ptree -> ptree.tree) ptrees in
    (* root the topologies using out-group to get rooted versions. *)
    let rtrees = List.map 
                    (fun utree -> (Rtree.root_at out_grp rt_id utree))
                    utrees in
    assert (List.for_all (fun x -> Rtree.is_root rt_id x) rtrees);
    (* build the consensus tree out of those rooted trees. *)
    let num_nodes = Rtree.get_num_nodes rt_id (List.hd rtrees) in
    let rtrees = List.map (fun rtree -> (rt_id, rtree)) rtrees in
    let croot_id, ctree, nd_indx_tbl, indx_nd_tbl = 
        Gen_rtree.build_ctree maj num_nodes rtrees in
        croot_id, ctree, nd_indx_tbl, indx_nd_tbl

let get_component_root handle ptree = 
    All_sets.IntegerMap.find handle ptree.component_root

let change_component_root handle root ptree = 
    let com_root = All_sets.IntegerMap.add handle root ptree.component_root in
    {ptree with component_root = com_root}


let choose_leaf tree = 
    match 
    All_sets.IntegerMap.fold (fun x n acc -> match n with Tree.Leaf _ -> Some x | _
    -> acc) tree.tree.Tree.u_topo None
    with
    | Some v -> v
    | None -> failwith "A tree with no leafs?"

(** Functor to allow for SPR/TBR searches over phylogenetics trees of 
    different types of characters. *)
module Search (Node : NodeSig.S) (Edge : Edge.EdgeSig with type n = Node.n) 
    (Tree_Ops : Tree_Operations with type a = Node.n with type b = Edge.e) = struct

        type a = Tree_Ops.a
        type b = Tree_Ops.b
    let features meth lst = Tree_Ops.features meth lst

    let downpass = Tree_Ops.downpass
    let uppass = Tree_Ops.uppass
    let diagnosis x = uppass (downpass x)

    (** Function used for debugging purposes...*)
    let verify_cost old_cost old_delta cdnd new_delta j1 j2 t_delta pt =
        let t, _ = Tree_Ops.join_fn [] j1 j2 pt in
        let new_cost = get_cost `Adjusted t in
        match new_delta with
        | Cost new_delta ->
              let supposed_new_cost = old_cost -. old_delta +. new_delta in
              if supposed_new_cost <> new_cost
              then begin
                  if new_cost > supposed_new_cost
                  then print_endline ("cost_fn optimistic: "
                                      ^ string_of_float
                                      (new_cost -. supposed_new_cost)
                                      ^ " (delta reported was " ^
                                      string_of_float new_delta ^ ")")
                  else print_endline ("cost_fn pessimistic: "
                                      ^ string_of_float
                                      (supposed_new_cost -. new_cost));
                  if debug_verify_costfn_printnodes
                  then begin
                      let clade_node = get_node_data cdnd pt in
                      print_endline ("clade node: "
                                     ^ Tree_Ops.string_of_node clade_node);
                      begin
                          match j1 with
                          | Tree.Single_Jxn h ->
                                print_endline ("trying to join at "
                                               ^ Tree_Ops.string_of_node
                                               (get_node_data
                                                    (Tree.int_of_id h) pt))
                          | Tree.Edge_Jxn (n1, n2) ->
                                print_endline ("trying to join at "
                                               ^ Tree_Ops.string_of_node
                                               (get_node_data
                                                    (Tree.int_of_id n1) pt));
                                print_endline ("              and "
                                               ^ Tree_Ops.string_of_node
                                               (get_node_data
                                                    (Tree.int_of_id n2) pt))
                      end
                  end
              end
        | NoCost ->
              (* Just check that the actual delta > old delta *)
              if new_cost < old_cost
              then print_endline ("..cost_fn skips incorrectly!")
        

    (** [make_wagner_tree ptree join_fn cost_fn]
        @param ptree The initial ptree that is just a bunch of single nodes
                     with node data associated.
        @param join_fn function used to join the nodes to build the wagner 
                       tree.
        @param cost_fn function to determine the cost of the tree.
        @return the wagner tree. i.e. best spr tree over the given data. *)
    let make_wagner_tree ?(sequence) ptree 
            (srch_mgr : (Tree_Ops.a, Tree_Ops.b) wagner_mgr) 
            (create_tabu_mgr : (('a, b) p_tree -> int -> (Tree_Ops.a, Tree_Ops.b)
            wagner_edges_mgr))
            = 
        let nodes = 
            match sequence with
            | None -> Tree.handle_list ptree.tree
            | Some r -> List.map (fun x -> handle_of x ptree) r 
        in
        (* make sure you have atleast two nodes to build a tree 
         *  (one edge only) of type a -- b *)
        match nodes with
        | n1 :: n2 :: rest ->
                let status = 
                    Status.create "Wagner" (Some (2 + List.length rest)) ""
                in
                (* build one edge tree with h1 h2 *)
                let h1 = Tree.get_handle_id n1 ptree.tree
                and h2 = Tree.get_handle_id n2 ptree.tree in
                let j1, j2 =
                    let get_corrected_jnx nd =
                        match get_node nd ptree with
                        | Tree.Single _ -> Tree.Single_Jxn nd
                        | Tree.Leaf (x, y)
                        | Tree.Interior (x, y, _, _) -> Tree.Edge_Jxn (x, y)
                    in
                    (get_corrected_jnx h1), (get_corrected_jnx h2)
                in
                let ptree, tree_delta = (Tree_Ops.join_fn [] j1 j2 ptree) in
                (* Now we ensure that the root is located in between the two 
                * handles that we just joined. This is needed for constrained
                * building. *)
                let ptree  = 
                    let l, r, _ = tree_delta in
                    let new_vertex x = 
                        match x with
                        | `Single (x, _)
                        | `Edge (x, _, _, _) -> x
                    in
                    let l = new_vertex l 
                    and r = new_vertex r in
                    let tree, inc = 
                        Tree_Ops.reroot_fn true (Tree.Edge (l, r)) ptree 
                    in
                    Tree_Ops.incremental_uppass tree inc
                in
                let cst = get_cost `Adjusted ptree in
                (* function adds the given nd to each of the edges of pt and
                * picks the tree/s according to some optimality criterion. *)
                let add_node_everywhere (pt, cst, tabu_mgr) nd srch_mgr =
                    let j2, nd_data =
                        match (get_component_root nd ptree).root_median with
                        | None -> assert false
                        | Some ((`Edge (x, y)), z) -> (Tree.Edge_Jxn (x, y)), z
                        | Some ((`Single x), z) -> (Tree.Single_Jxn x), z
                    in
                    tabu_mgr#next_clade nd_data;
                    (* function to add a node to an edge and determine the 
                    * optimality of the resulting tree. *)
                    let add_node_to_edge e srch_mgr tabu_mgr = 
                        let Tree.Edge(e1, e2) = e in
                        let h1 = (Tree.get_node_id e1 pt.tree) 
                        and h2 = (Tree.get_node_id e2 pt.tree) in
                        let j1 = Tree.Edge_Jxn(h1, h2) in
                        let status:t_status = 
                            (srch_mgr#process Tree_Ops.cost_fn infinity 
                                 nd_data Tree_Ops.join_fn j1 j2 pt tabu_mgr) 
                        in
                        status, srch_mgr
                    in
                    (* Sequentially add rest of the nodes keeping
                    * the best tree/s *)
                    let srch_mgr = 
                        let rec do_all_edges srch_mgr tabu_mgr =
                            match tabu_mgr#next_edge with
                            | None -> srch_mgr
                            | Some e ->
                                    let _, mgr = 
                                        add_node_to_edge e srch_mgr tabu_mgr
                                    in
                                    do_all_edges mgr tabu_mgr
                        in
                        do_all_edges srch_mgr tabu_mgr
                    in
                    (* There has to be at least one new tree *)
                    assert(srch_mgr#any_trees) ;
                    () 
                in
              (* sequentially add rest of the nodes to the
              * tree *)
              let rec seq_add nodes srch_mgr added = 
                  match nodes with
                  | [] -> srch_mgr
                  | nd :: rest ->
                          let n_srch_mgr = (srch_mgr#clone) in
                          assert(n_srch_mgr#any_trees = false) ;
                          while srch_mgr#any_trees do 
                              let (_, cst, _) as it = srch_mgr#next_tree in
                              Status.full_report ~msg:("Wagner tree with cost "
                              ^ string_of_float (cst))
                              ~adv:(added) status;
                              add_node_everywhere it nd n_srch_mgr;
                              n_srch_mgr#evaluate;
                          done;
                          (seq_add rest n_srch_mgr (added + 1))
              in
              let tabu_mgr = create_tabu_mgr ptree h1 in
              srch_mgr#init [(ptree, cst, Cost(infinity), tabu_mgr)] ;
              let result = (seq_add rest srch_mgr 2) in
              Status.finished status;
              result
              (* need at least two nodes *)
      | _ -> 
              let msg = "POY requires at least two taxa in order to start a \
              search. Please load some data." in
              Status.user_message Status.Error msg;
              failwith "Wagner building  requires a tree with at least two \
              vertices."

    
let trees_considered = ref 0

let fix_edge ptree ((Tree.Edge (e1, e2)) as edge) =
    (* if the edge is in tree, return the edge. *)
    if (is_edge edge ptree) then
        edge
    (* otherwise, check whether the reversed-edge is in the tree. *)
    else if (is_edge (Tree.Edge (e2, e1)) ptree) then
        Tree.Edge(e2, e1) 
    (* neither the edge nor its reverse was found in tree i.e. tabu and
     * tree are out of sync. *)
    else failwith "fix_edge"


(** [joins ptree loc tabu joiner] is a utility function that takes every [loc]
    (i.e. [`Left] or [`Right]) join edge in [tabu] and passes the join junction to
    [joiner], a function that will perform the join and return a [Tree] status
    code *)
let rec joins ?(counter=0) ptree loc tabu joiner =
    match tabu#join_edge loc with
    | None ->
          if ndebug
          then odebug ("no more places to join");
          if ndebug then
              odebug ("Total number of joins performed is : " ^ 
              string_of_int counter);
          (Tree.Continue, counter)
    | Some edge ->
          let Tree.Edge (je1, je2) = fix_edge ptree edge in
          if ndebug
          then odebug ("trying to join at "
                       ^ string_of_int je1
                       ^ ","
                       ^ string_of_int je2);
          let j1 = 
              try
                  Tree.Edge_Jxn
                  (Tree.get_node_id je1 ptree.tree,
                  Tree.get_node_id je2 ptree.tree) 
              with
              | err ->
                      print_endline "Ptree.joins";
                      raise err
          in
          let status = joiner j1 in
          match status with
          | Tree.Break -> (Tree.Break, succ counter)
          | Tree.Continue -> joins ~counter:(succ counter) ptree loc tabu joiner
          | Tree.Skip -> (Tree.Skip, succ counter)

(** [spr_join search break_delta clade_node clade_id incremental j2 tabu
    ptree] performs a set of SPR joins.  The parameters are:
    @param search The search manager to which to feed the results
    @param break_delta Cost difference from breaking the previous edge
    @param clade_node Data for clade root node
    @param clade_id Code of clade root node
    @param incremental Incremental update information obtained from the break,
    which should be updated on the join
    @param j2 Junction in clade on which to join
    @param tabu Tabu manager to provide join edges
    @param ptree Tree upon which to operate (already broken)
*)
let spr_join search
        (break_delta : float)
        (clade_node : Tree_Ops.a) clade_id incremental j2 tabu ptree =
    let joiner j1 =
        if ndebug_traject_spr
        then odebug ("traj: spr join: " ^ Tree.string_of_jxn j1 ^ " and " ^
                         Tree.string_of_jxn j2);
        search#process
            Tree_Ops.cost_fn
            break_delta                 (* break delta *)
            clade_node                  (* clade node *)
            (* clade_id *)
            Tree_Ops.join_fn
            incremental
            j1 j2                       (* junctions *)
            tabu                        (* current tabu..? *)
            ptree in                    (* post-break tree *)
    joins ptree `Left tabu joiner


let spr_step
        ptree
        (tabu : (Tree_Ops.a, Tree_Ops.b) tabu_mgr)
        (search : (Tree_Ops.a, Tree_Ops.b) search_mgr) =

    let nbreaks = ref 0 in
    let njoins = ref 0 in
    let njoins_last = ref 0 in
    if ndebug_traject_spr
    then odebug "traj: spr: begin";
    let swap_tree_and_clade (ptree, tree_delta, break_delta, clade_id, 
        clade_node, incremental) : ((Tree_Ops.a, Tree_Ops.b) p_tree * 
        Tree.break_delta * float * int * 'a * incremental list)=
        let old_tree_delta, old_clade_delta = tree_delta in
        match old_tree_delta with
        | `Edge (v, l1, l2, lst) ->
                let edge = Tree.Edge (l1, l2) in
                let ptree, newincr = Tree_Ops.reroot_fn false edge ptree in
                let new_tree_delta = (old_clade_delta, old_tree_delta) in
                let ptree = Tree_Ops.incremental_uppass ptree incremental in
                let clade_root = get_component_root l1 ptree in
                (match clade_root.root_median with
                | Some (_, clade_node) ->
                        (ptree, new_tree_delta, break_delta, v, clade_node, 
                        newincr)
                | None -> failwith "Unexpected swap_tree_and_clade")
        | `Single (v, _) -> 
                let new_tree_delta = old_clade_delta, old_tree_delta in
                let ptree = Tree_Ops.incremental_uppass ptree incremental in
                let clade_node = get_node_data v ptree in
                (ptree, new_tree_delta, break_delta, v, clade_node, [])
    in
    (* break edge (e1, e2) and treat the e2 part as the clade *)
    (* then try to reattach the clade in each tabu#join_edge location *)
    let break_side (ptree, tree_delta, break_delta, clade_id, 
        clade_node, incremental) =
        if ndebug_break_delta
        then odebug ("break_delta is " ^ string_of_float break_delta);
        let t1_h, t2_h = Tree.get_break_handles tree_delta (ptree.tree) in
        let t1_h, t2_h = handle_of t1_h ptree, handle_of t2_h ptree in
        let j2 = Tree.side_to_jxn (let (l, r) = tree_delta in r) in

        let tabu = tabu#clone in
        tabu#update_break
            ptree tree_delta t1_h t2_h clade_id;

        let res, joins =
            spr_join search break_delta clade_node clade_id incremental j2 tabu
                ptree in
        njoins := !njoins + joins;
        njoins_last := joins;
        if ndebug_traject_summary
        then odebug (string_of_int joins ^ " joins");
        res
    in

    (* Try to break in each possible location *)
    let rec break () =
        match tabu#break_edge with
        | None -> if ndebug then odebug "SPR: no more edges to break"
        | Some e ->
              let Tree.Edge (e1, e2) = 
                  let break = fix_edge ptree e in
                  search#breakin break;
                  break
              in
              if ndebug then odebug ("SPR: breaking at edge "
                                     ^ string_of_int e1
                                     ^ ","
                                     ^ string_of_int e2);
(*              let test = Tree.depths ptree.tree in*)
(*              odebug ("SPR: depth " ^ string_of_int (All_sets.TupleMap.find (e2,*)
(*              -1)*)
(*              test));*)
              incr nbreaks;
              if ndebug || ndebug_traject_summary
              then odebug "breaking";
              let (_, _, bdelta, _, _, _) as breakage = 
                  Tree_Ops.break_fn (e1, e2) ptree 
              in
              if bdelta > 0. then
                  let res =
                      if not (Tree.is_leaf e1 ptree.tree) then begin
                              if ndebug || ndebug_traject_summary
                              then odebug "checking left...";
                              break_side breakage
                          end else Tree.Continue
                  in
                  match res with
                  | Tree.Break -> ()
                  | _ ->
                          if not (Tree.is_leaf e2 ptree.tree) then begin
                              if ndebug || ndebug_traject_summary
                              then odebug "checking right...";
                              let breakage = swap_tree_and_clade breakage in
                              match break_side breakage with
                              | Tree.Break -> ()
                              | _ -> break ()
                          end else break ()
            else break ()
    in

    break ();
    if ndebug_traject_summary
    then odebug (string_of_int !nbreaks ^ " breaks, "
                 ^ string_of_int !njoins ^ " joins");
    search(* , (!nbreaks, !njoins, !njoins_last) *)

(** [tbr_join search break_delta clade_node clade_id incremental j2 tabu
    ptree] performs a set of TBR joins.  The parameters are:
    @param search The search manager to which to feed the results
    @param tabu Tabu manager to provide join edges
    @param ptree Tree upon which to operate (already broken)
    @param j2 Junction in clade on which to join first
    @param clade_node Data for clade root node
    @param break_delta Cost difference from breaking the previous edge
    @param clade_id Code of clade root node
*)
let rec tbr_join search ?(updt=[]) tabu ?(rerooted=false) ptree j2 clade_node
        break_delta =
    let reroots = ref 0 in
    let original_ptree = ptree in
    let rec reroot updt tabu rerooted ptree j2 the_node =
        let joiner j1 =
            if ndebug then odebug "TBR: Joining";
            if ndebug_traject_tbr
            then odebug ("traj: tbr join: " ^ Tree.string_of_jxn j1 ^ " and " ^
                             Tree.string_of_jxn j2);
            let tabu = tabu#clone in
            search#process
                Tree_Ops.cost_fn
                break_delta             (* break delta *)
                the_node                (* clade node *)
                (* clade_id *)
                Tree_Ops.join_fn
                updt              (* Incremental uppass has been done already *)
                j1 j2                   (* junctions *)
                tabu                    (* current tabu..? *)
                ptree                   (* post-break tree *)
        in
        match joins ptree `Left tabu#clone joiner with
        | Tree.Break, _ -> Tree.Break
        | Tree.Continue, _ 
        | Tree.Skip, _ ->
              match tabu#reroot_edge `Right with
              | None -> Tree.Continue
              | Some reroot_edge -> begin
                    if ndebug then odebug ("TBR: Rerooting");
                    let ptree, updt = 
                        Tree_Ops.reroot_fn false reroot_edge original_ptree 
                    in
                    let j2, the_node =
                        let Tree.Edge (e1, e2) = reroot_edge in
                        let (e1, e2) = (get_handle_id e1 ptree,
                                        get_node_id e2 ptree) in
                        let root = get_component_root e1 ptree in
                        match root.root_median with
                        | Some (`Edge (l, m), v) ->
                              Tree.Edge_Jxn (e1, e2), v
                        | _ -> failwith "Ptree.tbr_step.reroot"
                    in
                    incr reroots;
                    reroot updt tabu true ptree j2 the_node
                end
    in
    let res = reroot updt tabu rerooted ptree j2 clade_node in
    if ndebug then 
        odebug ("Total number of reroots performed " ^ string_of_int
                    !reroots);
    res


let tbr_step
        ptree
        (tabu : (Tree_Ops.a, Tree_Ops.b) tabu_mgr)
        (search : (Tree_Ops.a, Tree_Ops.b) search_mgr) =

    if ndebug_traject_tbr
    then odebug "traj: tbr: begin";
    (* break edge (e1, e2) and treat the e2 part as the clade *)
    (* then try to reattach the clade in each tabu#join_edge location *)
    let break_side_reroot tabu e1 e2 =
        search#breakin (Tree.Edge (e1, e2));
        let (ptree, tree_delta, break_delta, clade_id, clade_node, incremental) = 
            Tree_Ops.break_fn (e1, e2) ptree in
        (* NOTE: we don't update based on [incremental], so the tree is in a
           slightly inconsistent state.  This is OK in our particular case,
           since we then proceed to reroot the tree. *)

        let t1_h, t2_h = Tree.get_break_handles tree_delta (ptree.tree) in
        let j2 = Tree.side_to_jxn (let (l, r) = tree_delta in r) in

        if ndebug then odebug ("Breaking at "
                               ^ string_of_int t1_h
                               ^ ","
                               ^ string_of_int t2_h);
        if ndebug then begin
            let edges = get_pre_order_edges t2_h ptree in
            let len = string_of_int (List.length edges) in
            odebug ("Subtree with " ^ len ^ " edges");
        end;
        tabu#update_break
            ptree tree_delta t1_h t2_h clade_id;

        tbr_join search tabu ptree j2 clade_node break_delta
    in

    (* Try to break in each possible location *)
    let breaks = ref 0 in
    let rec break () =
        incr breaks;
        match tabu#break_edge with
        | None -> ()
        | Some e ->
              let Tree.Edge (e1, e2) = e in
              if ndebug_traject_tbr
              then odebug ("traj: tbr break: " ^ string_of_int e1 ^ " -> " ^
                               string_of_int e2);
              let res = break_side_reroot tabu#clone e1 e2 in
              match res with
              | Tree.Break -> ()
              | _ -> break ()
    in

    break ();
    if ndebug then begin
        odebug ("TBR: done. number of breaks: " ^ string_of_int !breaks);
    end;
    search(* , !breaks *)

let tbr_single search =
    while search#any_trees do
        let (tree, cost, tabu) = search#next_tree in
        ignore(tbr_step tree tabu search)
    done;
    search

type searcher =
        (Tree_Ops.a, Tree_Ops.b) search_mgr ->
            (Tree_Ops.a, Tree_Ops.b) search_mgr
type search_step = 
        (Tree_Ops.a, Tree_Ops.b) p_tree ->
            (Tree_Ops.a, Tree_Ops.b) tabu_mgr ->
                searcher

let search (searcher, name) search =
    let status = Status.create name None ("Searching") in
    while search#any_trees do
        let (ptree, cost, tabu) = search#next_tree in
        Status.full_report ~adv:(int_of_float cost) status;
        searcher ptree tabu#clone search;
    done;
    Status.finished status;
    search

(* This function will not find the local optimum, it will return as soon as a
* better tree is found. *)
let search_local_next_best (searcher, name) (search : (Tree_Ops.a, Tree_Ops.b)
                                                 search_mgr)
        : (Tree_Ops.a, Tree_Ops.b) search_mgr =
    let status = Status.create ("Single " ^ name) None ("Searching") in
    let ptree, cost, tabu = search#next_tree in
    searcher ptree tabu#clone search;
    Status.finished status;
    search

let spr_simple = search (spr_step, "SPR")
let tbr_simple = search (tbr_step, "TBR")

let spr_single = search_local_next_best (spr_step, "SPR")
let tbr_single = search_local_next_best (tbr_step, "TBR")

let alternate spr tbr search =
    let find_best_cost lst = 
              List.fold_left (fun best (_, cost, _) -> if best < cost then best
              else cost) max_float lst
    in
    let status = Status.create "Alternate" None ("") in
    let () = Status.full_report ~msg:("Beginning search") status in
    let rec try_spr prev_best search = match search#any_trees with
    | false -> search
    | true ->
          Status.full_report ~msg:("SPR search") status;
          let search = spr search in
          (* SPR is done---run TBR steps on the results *)
          Status.full_report
              ~msg:"Performing TBR swapping" status;
          let results = search#results in
          let best_cost = find_best_cost results in
          let search = search#clone in
          let () = search#init
              (List.map (fun (tree, cost, tabu) -> tree, cost, NoCost, tabu)
                   results) in
          let search = tbr search in
(*           let search = Sexpr.fold_status *)
(*               "TBR swapping" *)
(*               (fun search (tree, cost, tabu) -> *)
(*                    tbr_step tree tabu search) search (Sexpr.of_list results) in *)
            let new_cost = find_best_cost search#results in
            if (new_cost < best_cost && new_cost < prev_best) ||
            search#should_repeat then 
                let search = search#clone in
                let () = search#init
                    (List.map (fun (tree, cost, tabu) -> tree, cost, NoCost, tabu)
                    results) 
                in
                try_spr new_cost search
            else search
              
    in
    let search = try_spr max_float search in
    Status.finished status;
    search

let repeat_until_no_more tabu_creator neighborhood queue =
    Status.user_message Status.Information "Starting on tree";
    let rec go queue =
        let queue = neighborhood queue in
        if queue#should_repeat then 
            let results = queue#results in
            let queue = queue#clone in
            let _ =
                queue#init
                (List.map (fun (tree, cost, tbu) -> tree, cost,
                NoCost, tabu_creator tree)
                results);
                Status.user_message Status.Information "Iterating";
            in
            go queue
        else queue
    in
    let res = go queue in
    Status.user_message Status.Information "Finished on tree";
    res

    (** @return the number of trees considered during the search. *)
    let get_trees_considered () = !trees_considered

    (** @return Sets the number of trees considered to zero. *)
    let reset_trees_considered () = trees_considered := 0

(** {2 Tree Fusing} *)
type ('a, 'b) fuse_locations =
        (('a, 'b) p_tree * Tree.u_tree * Tree.edge) list Sexpr.t

let fuse_all_locations ?min ?max trees =
    let min = Some 3 in
    let filter = match min, max with
    | None, None -> (fun _ -> true)
    | Some min, None -> (fun (_, s) -> s >= min)
    | None, Some max -> (fun (_, s) -> s <= max)
    | Some min, Some max -> Tree.fuse_cladesize ~min ~max in
    let trees = List.map (fun t -> (t, t.tree)) trees in
    Tree.fuse_all_locations ~filter trees

let fuse source_arg target_arg =
    (* reroot if necessary *)
    let maybe_reroot ((tree, utree, (Tree.Edge(efrom, eto) as edge)) as arg) =
        if is_edge edge tree
        then arg
        else
            let tree, updt = Tree_Ops.reroot_fn false edge tree in
            let tree = Tree_Ops.incremental_uppass tree updt in
            tree, tree.tree, edge in
    let source, source_u, sedge = maybe_reroot source_arg in
    let target, target_u, tedge = maybe_reroot target_arg in

    let u_tree =
        Tree.fuse ~source:(source_u, sedge) ~target:(target_u, tedge) in
    let res = { target with tree = u_tree } in
    diagnosis res

(** [fuse_generations trees max_trees tree_weight tree_keep iterations process]
    runs a genetic algorithm-style search using tree fusing.  The function takes a
    list of trees to start with, the max number of trees and a method for keeping
    trees, a method for weighting trees, a number of iterations to perform, and a
    function to process new trees *)
let fuse_generations trees max_trees tree_weight tree_keep iterations
        process (cmin, cmax) =
    if 2 > List.length trees
    then failwith "Tree fusing: must have at least two trees";
    let status = Status.create "Tree Fusing" (Some iterations)
        "" in
    let () = Status.full_report status in
    
    let remove_worst_tree trees =
        let maxcost =
            List.fold_left (fun m tree -> max m (get_cost `Adjusted tree)) 0. trees
        in
        let rec rem lst = match lst with
        | tree :: trees -> if get_cost `Adjusted tree = maxcost then trees
          else tree :: rem trees
        | [] -> assert false in
        rem trees
    in
    let limit_num trees =
        let len = List.length trees in
        let trees = ref trees in
        for i = (len - max_trees) downto 1 do
            trees := remove_worst_tree !trees
        done; !trees in

    let keeper new_trees source target trees' = match tree_keep with
    | `Best ->
          let old_trees = source :: target :: trees' in
          limit_num (List.rev_append new_trees old_trees)
    | `Better ->
          let target_cost = get_cost `Adjusted target in
          let new_trees =
              List.filter (fun t -> target_cost >= get_cost `Adjusted t) new_trees in
          match new_trees with
          | [] -> source :: target :: trees'
          | new_trees ->
                (* omit target if we're running out of space *)
                let old_trees =
                    if List.length new_trees + List.length trees' + 1
                        >= max_trees
                    then source :: trees'
                    else source :: target :: trees'
                in
                limit_num (List.rev_append new_trees old_trees)
    in

    let rec choose_remove ?(i=1) f weights items = match weights with
    | w :: ws ->
          if w >= f
          then List.hd items, i, ws, List.tl items
          else let c, i, ws, is =
              choose_remove ~i:(succ i) (f -. w) ws (List.tl items) in
          c, i, w :: ws, (List.hd items) :: is
    | [] -> failwith "choose_remove" in
    let choose_remove f weights items =
        (* DEBUG: these assertions should not fail; if they do, we should have
           some sort of failsafe behavior. *)
        (* TODO: make sure there are at least two items (see below)... *)
        assert (List.length weights = List.length items);
        assert (f <= List.fold_left ( +. ) 0. weights);
        try
            choose_remove f weights items
        with Failure "choose_remove" ->
            begin
                let msg =
                    "Warning: choose_remove: "
                    ^ "total weight "
                    ^ string_of_float (List.fold_left ( +. ) 0. weights)
                    ^ "; f = "
                    ^ string_of_float f
                    ^ "; len = "
                    ^ string_of_int (List.length weights)
                in
                odebug msg;
                List.hd items, 1, List.tl weights, List.tl items
            end
    in

    let rec gen trees iter = if iter > iterations then trees else begin
        let trees = limit_num trees in
        Status.full_report ~adv:iter status;
        let weights = List.map tree_weight trees in
        let wsum = List.fold_left (+.) 0. weights in
        let source, snum, weights, trees' =
            choose_remove (Random.float wsum) weights trees in

        let wsum = List.fold_left (+.) 0. weights in
        let target, tnum, weights, trees' =
            choose_remove (Random.float wsum) weights trees' in

        let msg = ("tree #" ^ string_of_int snum
                   ^ " [" ^ string_of_float (get_cost `Adjusted source) ^ "]"
                   ^ " -> " ^ "tree #" ^ string_of_int tnum
                   ^ " [" ^ string_of_float (get_cost `Adjusted target) ^ "]") in

        let locations = fuse_all_locations ~min:cmin ~max:cmax [source; target] in
        let location = Sexpr.choose_random locations in
        match location with
        | None ->                       (* no good location *)
              Status.message status (msg ^ ": no join locations");
              gen trees (succ iter)
        | Some [t1; t2] ->
              let t1, t2 =
                  match t1 with
                  | (tree, _, _) when tree == source -> t1, t2
                  | (tree, _, _) when tree == target -> t1, t2
                  | _ -> assert false in
              let new_tree = fuse t1 t2 in
              Status.full_report ~msg:(msg ^ ": tree with cost "
                                       ^ string_of_float (get_cost `Adjusted new_tree))
                  status;
              let new_trees = process new_tree in

              let trees = keeper new_trees source target trees' in
              gen trees (succ iter)
        | Some _ ->                     (* why?! *)
              assert false

    end in
    let res = gen trees 1 in
    Status.finished status;
    limit_num res

(** [convert_to tree d]
    @param tree the Parser.Tree.t that is being used to build trees.
    @param d Data.d the data associated with the parser tree.
    @return p_tree that corresponds to the Parser.Tree.t *)
let convert_to tree (d, nd_data_lst) = 
    (* convert the Parser.Tree.t to Tree.u_tree *)
    let ut = Tree.convert_to tree d in
    let pt = { empty with tree = ut } in
    (* function to add leaf-node data to the ptree. *)
    let data_adder ptree nd = 
        (add_node_data (Node.taxon_code nd) nd ptree) in
        (List.fold_left data_adder pt nd_data_lst)

(** [build_trees tree]
    @param tree the ptree which is being converted into a Parser.Tree.t
    @param str_gen is a function that generates a string for each vertex in the
    tree
    @param collapse is a function that check weather or not a branch can be
    collapsed.
    @return the ptree in the form of a Parser.Tree.t *)
let build_trees (tree : Tree.u_tree) str_gen collapse root =
    let sortthem a b ao bo data ad bd =
        match String.compare ao bo with
        | 0 | 1 -> 
                Parser.Tree.Node (b @ a, data), bd + 1, bo
        | _ -> 
                Parser.Tree.Node (a @ b, data), bd + 1, ao
    in
    let get_children = function
        | Parser.Tree.Leaf _ as x -> [x]
        | Parser.Tree.Node (lst, _) -> lst
    in
    let rec rec_down node prev_node =
        match node with
        | Tree.Leaf (self, parent) -> 
              let data = str_gen self in
              Parser.Tree.Leaf data, 0, data
        | Tree.Interior (our_id, _, _, _) ->
              let (ch1, ch2) = 
                  assert (prev_node <> Tree.get_id node);
                  Tree.other_two_nbrs prev_node node in
              let a, ad, ao = rec_down (Tree.get_node ch1 tree) our_id in
              let b, bd, bo = rec_down (Tree.get_node ch2 tree) our_id in
              let a =
                  if collapse our_id ch1 then 
                      get_children a
                  else [a]
              and b =
                  if collapse our_id ch2 then 
                      get_children b
                  else [b]
              in
              let data = 
                  try str_gen our_id 
                  with Not_found -> ""
              in
              if bd > ad then
                  Parser.Tree.Node (a @ b, data), bd + 1, bo
              else if bd = ad then sortthem a b ao bo data ad bd
              else Parser.Tree.Node (b @ a, data), ad + 1, ao
        | Tree.Single _ -> failwith "Unexpected single"
    in
    let map : string Parser.Tree.t list =
        List.map
        (fun handle ->
             match Tree.get_node handle tree with
             | Tree.Leaf (self, parent) ->
                   let acc, _, _ = rec_down (Tree.get_node parent
                                           tree) handle in
                   let str = str_gen self in 
                   let nd_data = 
                       [(Parser.Tree.Leaf str); acc]
                   in
                   Parser.Tree.Node (nd_data, root)
             | Tree.Interior (self, par, ch1, ch2) ->
                   let par = Tree.get_node par tree in
                   let ch1 = Tree.get_node ch1 tree in
                   let ch2 = Tree.get_node ch2 tree in
                   let acc, accd, acco = rec_down par handle in
                   let acc1, acc1d, acc1o = rec_down ch1 handle in
                   let acc2, acc2d, acc2o = rec_down ch2 handle in
                   let data = 
                       try str_gen self with
                       | Not_found -> ""
                   in
                   if acc2d > acc1d then
                       if acc1d > accd then
                           Parser.Tree.Node ([Parser.Tree.Node ([acc; acc1], 
                           data); acc2], root)
                       else 
                           Parser.Tree.Node ([Parser.Tree.Node ([acc1; acc], 
                           data); acc2], root)
                   else if acc1d > accd then
                       if acc2d > accd then 
                           Parser.Tree.Node ([Parser.Tree.Node ([acc; acc2], 
                           (try str_gen self with Not_found -> "")); acc1], root)
                       else 
                           Parser.Tree.Node ([Parser.Tree.Node ([acc2; acc], 
                           (try str_gen self with Not_found -> "")); acc1], root)
                   else 
                       if acc2d > acc1d then
                           Parser.Tree.Node ([Parser.Tree.Node ([acc1; acc2], 
                           (try str_gen self with Not_found -> "")); acc], root)
                       else 
                           Parser.Tree.Node ([Parser.Tree.Node ([acc2; acc1], 
                           (try str_gen self with Not_found -> "")); acc], root)
             | Tree.Single self -> 
                     Parser.Tree.Node
                   ([(Parser.Tree.Leaf (str_gen self))], root))
            (All_sets.Integers.elements (Tree.get_handles tree))
    in
    map


    (** [build_tree tree]
        @param tree the ptree being converted into a Parser.Tree.t
        @return the tree with the smallest handle_id in ptree *)
    let build_tree tree strgen collapse root =
        let map = build_trees tree strgen collapse root in
        List.hd map

    let never_collapse a b c = false

    let collapse_as_needed tree code chld = 
        let data = get_node_data code tree in
        let datac = get_node_data chld tree in
        Node.is_collapsable `Any data datac

let extract_names pd ptree code = 
    let data = get_node_data code ptree in
    try Data.code_taxon (Node.taxon_code data) pd
    with _ -> ""

let extract_codes pd ptree code = string_of_int code

let build_tree_with_codes tree pd = 
    build_tree tree.tree (extract_codes pd tree) (collapse_as_needed tree) ""


let rec fold_2 f acc a b =
    match a, b with
    | (ha :: ta), (hb :: tb) -> fold_2 f (f acc ha hb) ta tb
    | [], [] -> acc
    | _, _ -> false

let rec compare acc a b = 
    acc &&
    (match a, b with
    | Parser.Tree.Leaf x, Parser.Tree.Leaf y ->  x = y
    | Parser.Tree.Node (ca, _), Parser.Tree.Node (cb, _) ->
            fold_2 compare true ca cb
    | _, _ -> false)


    let get_unique trees =
        match trees with 
        | tree :: _ ->
                let a, _ = Tree.choose_leaf tree.tree in
                let trees =
                    List.rev_map (fun x -> 
                        x, { x with tree = Tree.cannonize_on_leaf a x.tree }) trees 
                in
                let trees = 
                    List.rev_map 
                    (fun (x, y) -> x, (let _, z = cannonize (build_tree_with_codes
                    y "") in z)) 
                    trees
                in
                let rec remove_duplicated acc = function
                    |  (x, y) :: t -> 
                            let are_different (_, z) = not (compare true y z) in
                            remove_duplicated (x :: acc) (List.filter are_different
                            t) 
                    | [] -> acc
                in
                remove_duplicated [] trees
        | x -> x

let handle_collapse bool = 
    if bool then collapse_as_needed
    else (fun _ _ _ -> false)

(** [build_tree_with_names tree pd]
    @return What does this function do? *)
let build_tree_with_names collapse tree pd =
    let collapse_f = handle_collapse collapse tree in
    build_tree tree.tree (extract_names pd tree) collapse_f ""

let build_forest_with_names collapse tree pd =
    let collapse_f = handle_collapse collapse tree in
    build_trees tree.tree (extract_names pd tree) collapse_f ""

let build_tree_with_names_n_costs collapse tree pd cost = 
    let collapse_f = handle_collapse collapse tree in
    let extract_names code =
        let data = get_node_data code tree in
        match get_node code tree with
        | Tree.Interior (_, par, _, _) ->
                string_of_float (Node.total_cost (Some par) data)
        | _ ->
                try Data.code_taxon (Node.taxon_code data) pd
                with _ -> ""
    in
    build_tree tree.tree extract_names collapse_f cost

let build_forest collapse tree pd cost =
    let collapse_f = handle_collapse collapse tree in
    let extract_names code =
        let data = get_node_data code tree in
        match get_node code tree with
        | Tree.Interior (_, par, _, _) ->
(*                string_of_float (Node.total_cost (Some par) data)*)
                    ""
        | _ ->
                try Data.code_taxon (Node.taxon_code data) pd
                with _ -> ""
    in
    let trees = build_trees tree.tree extract_names collapse_f cost in
    trees

let build_forest_as_tree collapse tree pd cost =
    match build_forest collapse tree pd cost with
    | [tree] -> tree
    | trees ->
          Parser.Tree.Node (trees, "forest")

let build_forest_with_names_n_costs collapse tree pd cost = 
    let collapse_f = handle_collapse collapse tree in
    let extract_names code =
        let data = get_node_data code tree in
        match get_node code tree with
        | Tree.Interior (_, par, _, _) ->
                string_of_float (Node.total_cost (Some par) data)
        | _ ->
                try Data.code_taxon (Node.taxon_code data) pd
                with _ -> ""
    in
    build_trees tree.tree extract_names collapse_f cost

(** [disp_trees tree]
    @param str the title of the image.
    @param tree the ptree being displayed. 
    @return the tree is drawn either in graphical format or ascii format. *)
let disp_trees str tree strgen root =
    let trees = ref (build_trees tree.tree (strgen tree) (collapse_as_needed tree) root) in
    let ntrees = List.length !trees in
    for i = 1 to ntrees do
        print_endline (str ^ ": " ^ string_of_int i
        ^ "/" ^ string_of_int ntrees);
        AsciiTree.draw false stdout (List.hd !trees);
        trees := List.tl !trees
    done
(** [to_xml tree data f]
@param tree the ptree which is being converted into a Parser.Tree.t
@return the ptree in the form of a Parser.Tree.t *)
let to_xml ch tree d =
    let rec rec_down node prev_node =
        match node with
        | Tree.Leaf (self, parent) -> 
                output_string ch "<otu>\n";
                flush ch;
                let data = get_node_data self tree in
                Node.to_xml d ch data;
                output_string ch "</otu>\n";
                flush ch;
        | Tree.Interior (our_id, _, _, _) ->
                let data = get_node_data our_id tree in
                let (ch1, ch2) = 
                    assert (prev_node <> Tree.get_id node);
                    Tree.other_two_nbrs prev_node node in
                output_string ch "<htu>\n";
                flush ch;
                Node.to_xml d ch data;
                flush ch;
                rec_down (get_node ch1 tree) our_id;
                rec_down (get_node ch2 tree) our_id;
                output_string ch "</htu>\n";
                flush ch;
        | Tree.Single _ -> failwith "Unexpected single"
    in
    let print_item handle =
        output_string ch "<tree>\n";
        begin match get_node handle tree with
        | Tree.Leaf (self, parent) ->
                rec_down (get_node parent tree) handle
        | Tree.Interior (self, par, ch1, ch2) ->
                (* rec_down (get_node par tree) self; *)
                rec_down (get_node ch1 tree) self;
                rec_down (get_node ch2 tree) self;
                Node.to_xml d ch (get_node_data self tree);
        | Tree.Single self -> ()
        end;
        output_string ch "</tree>\n"
    in
    All_sets.Integers.iter print_item (get_handles tree)

end                                     (* module Search *)



module Fingerprint = struct
    type t = Tree.Fingerprint.t
    let fingerprint {tree=t} = Tree.Fingerprint.fingerprint t
    let compare = Tree.Fingerprint.compare
    let to_string = Tree.Fingerprint.to_string
    let empty = Tree.Fingerprint.empty
end

let get_leaves ?(init=[]) root t =
    pre_order_node_visit
        (fun parent node_id acc ->
            try
             match get_node node_id t with
             | Tree.Leaf _ -> (Tree.Continue,
                          (get_node_data node_id t) :: acc)
             | _ -> (Tree.Continue, acc)
            with 
            | Not_found as err -> 
                    Status.user_message Status.Information 
                    ("Failed with code " ^ string_of_int node_id);
                    raise err)
        root t init

let get_all_leaves t =
    All_sets.Integers.fold
        (fun h init -> get_leaves ~init h t)
        t.tree.Tree.handles
        []

let select_default adjusted cost =
    match adjusted with
    | Some v -> v
    | None -> cost

let assign_root_to_connected_component handle item cost adjusted ptree =
    let adjusted = select_default adjusted cost in
    let root = { 
        root_median = item;
        component_cost = cost;
        adjusted_component_cost = adjusted 
    } in
    let new_component = 
        All_sets.IntegerMap.add handle root ptree.component_root 
    in
    { ptree with component_root = new_component; }

let get_handle_cost adjusted ptree handle_id =
    let res = All_sets.IntegerMap.find handle_id ptree.component_root in
    match adjusted with
    | `Unadjusted -> res.component_cost
    | `Adjusted -> res.adjusted_component_cost

let set_component_cost cost adjusted root handle_id ptree = 
    let adjusted = select_default adjusted cost in
    try
        let root = { root_median = root; component_cost = cost;
        adjusted_component_cost = adjusted } in 
        let new_comp = 
            All_sets.IntegerMap.add handle_id root ptree.component_root 
        in
        { ptree with component_root = new_comp }
    with
    | Not_found ->
            let root = { root_median = root; component_cost = cost;
            adjusted_component_cost = adjusted } in
            let new_comp =
                All_sets.IntegerMap.add handle_id root ptree.component_root 
            in
            { ptree with component_root = new_comp; }

let get_leaves_ids ?(init=[]) root t =
    pre_order_node_visit
        (fun parent node_id acc ->
             match get_node node_id t with
             | Tree.Leaf _ -> (Tree.Continue,
                               node_id :: acc)
             | _ -> (Tree.Continue, acc))
        root t init

let get_all_leaves_ids t =
    All_sets.Integers.fold
        (fun h init -> get_leaves_ids ~init h t)
        t.tree.Tree.handles
        []

let wipe_costs ptree = 
    { 
        ptree with
        component_root = All_sets.IntegerMap.empty;
    }

(** [break_median ptree id] creates two trees from one by removing one of the
    edges incident to [id].  ([id] must be an internal node.)  If the neighbors
    of [id] are a, b, and c, it does this by comparing the medians (a, b), (b,
    c), and (a, c); it keeps the shortest edge. *)
let break_median_edge medianfn id ptree =
    match get_node id ptree with
    | Tree.Interior (self, a, b, c) ->
          let da, db, dc =
              get_parent_or_root_data self ptree,
              get_node_data b ptree,
              get_node_data c ptree in
          (* calculate the medians *)
            let self = Some self in
          let mab = medianfn self self da db in
          let mac = medianfn self self da dc in
          let mbc = medianfn self self db dc in
          (* break the longest edge *)
          if mab <= mac && mab <= mbc
          then Tree.Edge (id, c)
          else if mac <= mab && mac <= mbc
          then Tree.Edge (id, b)
          else if mbc <= mab && mbc <= mac
          then Tree.Edge (id, a)
          else failwith "break_median_edge"
    | _ -> failwith "break_median_edge1"

(** [jxn_of_handle ptree handle] returns the junction corresponding to a
    handle's position in the tree (forest) *)
let jxn_of_handle ptree handle =
    Tree.test_tree ptree.tree;
    match get_node handle ptree with
    | Tree.Single n -> Tree.Single_Jxn n
    | Tree.Interior (id, par, _, _)
    | Tree.Leaf (id, par) ->
          if ndebug_jxn_of_handle
          then (odebug ("jxn of handle with id " ^ string_of_int id
                        ^ " and parent " ^ string_of_int par);
                ignore(get_node par ptree);
                odebug "Parent node exists");
          Tree.Edge_Jxn (id, par)

(* A function that takes a map of counts of sets of terminals, and a tree, and
* adds the sets defined by tree to the counters *)
let add_tree_to_counters is_collapsable counters (tree : Tree.u_tree) = 
    let add_or_not addme set counters = 
        if (not addme) && Tree.CladeFPMap.mem set counters then
            let res = Tree.CladeFPMap.find set counters in
            set, Tree.CladeFPMap.add set (res + 1) counters
        else if not addme then
            set, Tree.CladeFPMap.add set 1 counters
        else set, counters
    in
    let rec add_all_clades node prev counters =
        let addme = is_collapsable prev node in
        let add_singleton () =
            let single = All_sets.Integers.singleton node in
            add_or_not false single counters
        in
        match Tree.get_node node tree with
        | Tree.Leaf _ 
        | Tree.Single _ -> add_singleton ()
        | Tree.Interior (_, a, b, c) ->
                let a, b = 
                    if a = prev then b, c
                    else if b = prev then a, c
                    else if c = prev then a, b
                    else failwith "Tree.consensus"
                in
                let sidea, counters = add_all_clades a node counters in
                let sideb, counters = add_all_clades b node counters in
                let mine = All_sets.Integers.union sidea sideb in
                add_or_not addme mine counters
    in
    let add_handle tree handle counters = 
        match Tree.get_node handle tree with
        | Tree.Leaf (a, b) 
        | Tree.Interior (a, b, _, _) ->
                let sidea, counters = add_all_clades a b counters in
                let sideb, counters = add_all_clades b a counters in
                let mine = All_sets.Integers.union sidea sideb in
                let _, res = add_or_not false mine counters in
                res
        | Tree.Single _ -> 
                let _, res = add_all_clades handle handle counters in
                res
    in
    All_sets.Integers.fold (add_handle tree) (Tree.get_handles tree) counters

let build_a_tree to_string denominator print_frequency coder trees (set, cnt) = 
    let module CladeFPSet = Set.Make (Tree.CladeFP.Ordered) in
    let all_trees, _, trees = 
        All_sets.Integers.fold (fun x (my_trees, sets, builttrees) ->
            try
                let trees, clades = All_sets.IntegerMap.find x builttrees in
                if CladeFPSet.mem clades sets then my_trees, sets, builttrees
                else
                    let sets = CladeFPSet.add clades sets in
                    trees :: my_trees, sets, builttrees
            with
            | Not_found ->
                    assert (1 = All_sets.Integers.cardinal set);
                    let code = All_sets.Integers.choose set in
                    coder := code;
                    let name = to_string code in
                    let newtree = Parser.Tree.Leaf name in
                    let trees = 
                        All_sets.IntegerMap.add code (newtree, set) builttrees 
                    in
                    newtree :: my_trees, sets, trees) 
        set 
        ([], CladeFPSet.empty, trees)
    in
    let tree = 
        match all_trees with
        | [] -> failwith "Tree.consensus2"
        | [all_trees] -> all_trees
        | all_trees ->
                let cnt = float_of_int cnt in
                let msg = 
                    if print_frequency then
                        Printf.sprintf "%.2f" (cnt /. denominator)
                    else ""
                in
                Parser.Tree.Node (all_trees, msg)
    in
    All_sets.Integers.fold (fun x trees ->
        All_sets.IntegerMap.add x (tree, set) trees)
    set
    trees

let make_tree majority_cutoff coder builder map = 
    let sets = 
        Tree.CladeFPMap.fold (fun set cnt lst ->
        if cnt >= majority_cutoff then (set, cnt) :: lst 
        else lst) 
        map
        []
    in
    let sets = 
        List.sort (fun (x, _) (y, _) -> 
        (All_sets.Integers.cardinal x) - 
        (All_sets.Integers.cardinal y)) sets
    in
    let trees = List.fold_left builder All_sets.IntegerMap.empty sets in
    let tree, _ = All_sets.IntegerMap.find !coder trees in
    tree

let (-->) a b = b a

let consensus is_collapsable to_string maj trees root =
    let number_of_trees = List.length trees in
    let print_frequency = maj <> number_of_trees 
    and number_of_trees = float_of_int number_of_trees 
    and coder = ref 0 in
    let tree_builder = 
        build_a_tree to_string number_of_trees print_frequency coder
    in
    trees 
    -->
        List.map (fun x -> 
            let y = get_parent root x in
            Tree.reroot (root, y) x.tree, is_collapsable x)
    --> List.fold_left (fun acc (x, y) ->
        add_tree_to_counters y acc x) Tree.CladeFPMap.empty 
    --> make_tree maj coder tree_builder

let preprocessed_consensus to_string maj num_trees map =
    let coder = ref 0 in
    let tree_builder = 
        build_a_tree to_string (float_of_int num_trees) true coder 
    in
    make_tree maj coder tree_builder map

let extract_counters sets set =
    Tree.CladeFPMap.fold (fun a b acc ->
        let freq = 
            try Tree.CladeFPMap.find a sets with
            | Not_found -> 0
        in
        Tree.CladeFPMap.add a freq acc ) 
    set
    Tree.CladeFPMap.empty

let supports to_string maj number_of_samples tree sets =
    let coder = ref 0 in
    let tree_builder = 
        build_a_tree to_string number_of_samples true coder
    in
    tree 
    --> add_tree_to_counters (fun _ _ -> false) Tree.CladeFPMap.empty 
    --> extract_counters sets
    --> make_tree maj coder tree_builder

(* A function that returns the bremer support tree based on the set of (costs, 
* tree) of sets, for the input tree *)
let bremer to_string cost tree generator file =
    let tree_generator = Parser.Tree.stream_of_file file in
    (* We first create a function that takes a map of clades and best cost found
    * for a tree _not_ containing the set, and a set of clades belonging to a
    * tree, with it's associated cost, and update the map according to the cost
    * for the set of clades, only if better. *)
    let replace_when_smaller map =
        let map = ref map in
        let cntr = ref 1 in
        let status = 
            Status.create "Bremer Estimation" None 
            "Comparing tree with trees in file" 
        in
        try
            while true do
                try 
                    Status.full_report ~adv:!cntr status;
                    let input_tree = tree_generator () in
                    let new_cost, sets = generator input_tree in
                    map :=
                        Tree.CladeFPMap.fold (fun my_clade best_cost acc ->
                        if (not (Tree.CladeFP.CladeSet.mem my_clade sets)) &&
                            ((new_cost - cost) < best_cost) then
                            Tree.CladeFPMap.add my_clade (new_cost - cost) acc
                        else acc) !map !map;
                    incr cntr;
                with
                | End_of_file as err -> raise err 
                | _ -> ()
            done;
            !map
        with
        | End_of_file -> 
                Status.finished status;
                !map
    in
    (** We create a map with all the sets of clades in the input tree *)
    let map : int Tree.CladeFPMap.t = 
        let map = 
            add_tree_to_counters (fun _ _ -> false) Tree.CladeFPMap.empty
            tree
        in
        Tree.CladeFPMap.map
        (fun _ -> max_int) map 
    in
    (* We now update that map with the best length found for each tree not
    * having the clade *)
    let coder = ref 0 in
    let tree_builder =
        build_a_tree to_string 1. true coder
    in
    map
    --> replace_when_smaller
    --> make_tree (-1) coder tree_builder
