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

(** The Tree module that provides support for local search of phylogenetic
    trees using SPR and TBR on unrooted trees. *)



(** {2 Types that define an unrooted tree topology **)

(** Type for id's in the node information *)
type id = int

(** The nodes have the following form: NodeType(id, _...) where id is the 
    node id. In these unrooted trees, the locations of the tuples have no special
    significance. The nodes are interpreted as 
        - Interior (id, nbr1, nbr2, nbr3)
        - Leaf (id, nbr1)
        - Single (id) *)
type node =
  | Single of id
  | Leaf of id * id
  | Interior of (id * id * id * id)

(** In an rooted tree, edges are implicitly directed according to parent ->
    child relationships. The edges in unrooted trees are directed based on a
    DFS from a handle. These edges are used only by unrooted trees to represent
    the directionality imposed on an undirected edge when traversing from a
    handle. *)
type edge = Edge of (id * id)

(** The edge comparator, allows comparison of edges so we
can create edge -> {% $\alpha$ %} maps. *)
module EdgeComparator :
  sig 
    type t = edge 
    val compare : edge -> edge -> int 
  end

(** The set of edges *)
module EdgeSet : Set.S with type elt = edge

(** The edge -> {% $\alpha$ %} map. *)
module EdgeMap : Map.S with type key = edge

(** The unrooted tree. Just as the rooted tree, the un-rooted 
    tree's topology is represented as a node_id -> node map. 
    [d_edges] are the directed edges imposed by handles in each
    of the connected components (when the tree is not fully connected). 
    There is one handle per connected component in the tree. *)
type u_tree = {
    tree_name : string option;
    u_topo : node All_sets.IntegerMap.t; (** Stores the tree topology as node
                                             records *)
    d_edges : EdgeSet.t;                (** Set of edges *)
    handles : All_sets.Integers.t;      (** Set of handle ids *)
    avail_ids : id list;                (** Unused ids available for reuse *)
    new_ids : id;                       (** Where to start numbering new ids *)
}



(** {2 Exceptions for handling errors *)

exception Invalid_Node_Id of id

exception Invalid_Handle_Id

exception Invalid_Edge

exception Invalid_End_Nodes

exception Node_Is_Handle



(** {2 Break and Join Type Information **)

(** A simple record to keep track of added and removed edges *)
type edge_delta = {
    added : edge list;
    removed : edge list;
}

(** The junctions where a tree can be broken to create two subtrees. *)
type break_jxn = id * id

(** Summarizes what a tree break or join operation did to a side of the tree
    being created or broken *)
type side_delta =
        [ `Single of id * bool
            (** the code of the single, whether or not a
                handle was created/removed here *)
        | `Edge of id * id * id * id option ]
            (** l (the node created/destroyed), l1, l2,
                id of the handle created/removed, if any *)

(** Joining can be done at an edge or with a single isolated node *)
type join_jxn = 
    | Single_Jxn of id
    | Edge_Jxn of id * id

(** Moving the root means changing the nodes in a path *)
type reroot_delta = id list

(** A break operates on two sides of an edge *)
type break_delta = side_delta * side_delta

(** A join joins two components, and removes a root on the right-hand
    component, resulting in a "rerooting" operation on that side *)
type join_delta = side_delta * side_delta * reroot_delta



(** {2 Define Basic Tree creation functions *)

(** The empty u_tree. *)
val empty : unit -> u_tree

(** Define a tree from a list of id's *)
val make_disjoint_tree : id list -> u_tree

(** Create a uniformly random tree from a list of ID's *)
val random : id list -> u_tree



(** {2 General Accessor Functions for the tree *)

(** [get_id n] get node id from node; the first int in the tuple *)
val get_id : node -> int

(** [int_of_id i] transform an id to an integer; identity function here *)
val int_of_id : id -> int

(** [is_handle i] checks if i is a handle in the tree *)
val is_handle : int -> u_tree -> bool

(** [is_edge e t] checks if e is an edge in t *)
val is_edge : edge -> u_tree -> bool

(** [normalized_edge e t] flips elements of the edge to ensure that it is
    contained in the tree,t *)
val normalize_edge : edge -> u_tree -> edge

(** [get_node_id i t] get node id from it's integer; ensures that it is valid  **)
val get_node_id : int -> u_tree -> id

(** [get_handle_id i t] get handle id from it's integer; ensures that it is valid  **)
val get_handle_id : int -> u_tree -> id

(** [get_handle_id i t] gets the handle attached to the subtree containing i *)
val handle_of : id -> u_tree -> id

(** [get_handles t] get the set of handles for the tree *)
val get_handles : u_tree -> All_sets.Integers.t

(** [handle_list t] get the list of handles in t *)
val handle_list : u_tree -> int list

(** [get_nodes t] get list of nodes *)
val get_nodes : u_tree -> node list

(** [get_node_ids t] get the integer keys of the nodes *)
val get_node_ids : u_tree -> id list

(** [get_node i t] get the node for i; neighbors information *)
val get_node : id -> u_tree -> node

(** [get_all_leaves t] return all the leaf nodes *)
val get_all_leaves : u_tree -> id list

(** [is_leaf i t] return if i is a Leaf(i,_) *)
val is_leaf : id -> u_tree -> bool

(** [is_single i t] return if i is a Single(i) *)
val is_single : id -> u_tree -> bool

(** [get_edge (a,b) t] Return an edge from a pair of ID's *)
val get_edge : id * id -> u_tree -> edge

(** [get_parent] get the parent of a node. In an unrooted context this may not
    be what is expected. *)
val get_parent : id -> u_tree -> id

(** Get a pre-order list of edges from a node *)
val get_pre_order_edges : id -> u_tree -> edge list

(** Return a list of edges in a tree *)
val get_edges_tree : u_tree -> edge list

(** [other_two_nbrs] get the other two neighbors of an interior node excluding
    the id passed. *)
val other_two_nbrs : id -> node -> id * id

(** Choose a leaf and return it and it's parent. The choice is done by the Set
    module. Although not specified in the implementation, it is not random *)
val choose_leaf : u_tree -> (id * id)

(** Verify that an edge exists in the tree *)
val verify_edge : edge -> u_tree -> bool


(** Get a Map of Handle Edges to depths of a tree *)
val depths : u_tree -> int All_sets.TupleMap.t



(** {2 Handle Manipulation functions *)

(** Return the path from a node to a handle. Return the handle id and an edge
    list that defines the path *)
val get_path_to_handle : id -> u_tree -> id * edge list

(** Move a handle to a node; return the tree with the new handle, and a path
    from the old handle to the new one *)
val move_handle : id -> u_tree -> u_tree * id list

(** Return a list of nodes from an id to a handle *)
val get_vertices_to_handle :  id -> u_tree -> id list

val fix_handle_neighbor : id -> id -> u_tree -> u_tree



(** {2 Break and Join Functions *)

val break_to_edge_delta : break_delta -> edge_delta
(** [break_to_edge_delta delta] returns the set of edges created and removed by
    a break operation *)

val join_to_edge_delta : join_delta -> edge_delta
(** [join_to_edge_delta delta] returns the set of edges created and removed by
    a join operation *)

val get_break_handles : break_delta -> u_tree -> id * id
(** [get_break_handles delta tree] returns the ids of the handles of the left
    and right components of a tree after a break *)

val break : break_jxn -> u_tree -> u_tree * break_delta
(** Break an edge defined by a jxn point *)

val join : join_jxn -> join_jxn -> u_tree -> u_tree * join_delta
(** Joint two join jxns as defined above. *)

val reroot : (id * id) -> u_tree -> u_tree
(** reroot a tree on an edge defined by the pair of ID's *)

val print_break_jxn : break_jxn -> unit
(** Debugging function for a break jxn *)

val print_join_1_jxn : join_jxn -> unit
(** Print join function labeled with a "1" *)

val print_join_2_jxn : join_jxn -> unit
(** Print join function labeled with a "2" *)



(** {2 Traversal Functions **)

(** Status of a traversal. *)
type t_status = 
    | Continue  (** Continue the current traversal. *)
    | Skip      (** Skip the subtree rooted at this node, for pre-order only *)
    | Break     (** Stop traversal and return current val. *)

(** [pre_order_edge_visit f id ptree accum]
    @param f function to be applied to each edge of the subtree
           of ptree rooted at id.
    @param id all the edges that are descendents of this node are visited
           in pre-order manner.
    @param ptree the ptree whose edges are being visited.
    @param accum the accumulator, a list.
    @return the return values of the function being applied on the edges. *)
val pre_order_edge_visit :
    (edge -> 'a -> t_status * 'a) -> id -> u_tree -> 'a -> 'a

(** [pre_order_node_visit f id bt ad acc]
    @param f function to applied to all the nodes in pre-order.
    @param id the node_id from where the traversal is started.
    @param bt the tree whose nodes are visited.
    @param acc the result of the function application on the edges.
    @return the function applied to the subtree of id as imposed by the
        handle in the component, the results are returned as a list. *)
val pre_order_node_visit :
    (id option -> id -> 'a -> t_status * 'a) -> id -> u_tree -> 'a -> 'a

(** [post_order_node_with_edge_visit] A function to traverse a tree from an
    edge, applying a function on this edge in each direction, applying [f] if
    the node is a leaf, or g if the node is interior, with two accumulators, one
    from the left and one from the right. Single nodes are ignored. *)
val post_order_node_with_edge_visit :
    (id -> id -> 'a -> 'a) -> (id -> id -> 'a -> 'a -> 'a) -> edge -> u_tree -> 'a -> 'a * 'a

(** [post_order_node_with_edge_visit_simple f e t a] is a simplified visitor
    function [f], which is applied on every (non single) vertex, starting in 
    the (hypothetical) root located between the two vertices of edge [e], over
    tree [t] with accumulator [a]. *)
val post_order_node_with_edge_visit_simple :
    (id -> id -> 'a -> 'a) -> edge -> u_tree -> 'a -> 'a

(** [pre_order_node_with_edge_visit_simple f e t a] is a simplified visitor
    function [f], which is applied on every (non single) vertex, starting in 
    the (hypothetical) root located between the two vertices of edge [e], over
    tree [t] with accumulator [a]. *)
val pre_order_node_with_edge_visit_simple :
    (id -> id -> 'a -> 'a) -> edge -> u_tree -> 'a -> 'a

(** [pre_order_node_with_edge_visit_simple_root rf f e t a] is a simplified
    visitor function [f], which is applied on every (non single) vertex, with a
    special function [rf] applied to the (hypothetical) root location between
    two vertices of edge [e], over tree [t] with accumulator [a]. *)
val pre_order_node_with_edge_visit_simple_root :
    (int -> int -> 'a -> 'a) -> edge -> u_tree -> 'a -> 'a

(** [post_order_node_visit f id bt ad acc]
    @param f function to applied to all the nodes in post-order.
    @param id the node_id from where the traversal is started.
    @param bt the tree whose nodes are visited.
    @param acc the result of the function application on the edges.
    @return the function applied to the subtree of id as imposed by the
            handle in the component, the results are returned as a list. *)
val post_order_node_visit :
    (id option -> id -> 'a -> t_status * 'a) -> id -> u_tree -> 'a -> 'a



(** {2 Distance Functions for trees *)

(* Create two sets defining the partition accross an edge *)
val create_partition : u_tree -> edge -> All_sets.Integers.t * All_sets.Integers.t

(* Calculate the Robinson Foulds distance between two trees *)
val robinson_foulds : u_tree -> u_tree -> int



(** {2 Functions for Re-Mapping Trees *)

val cannonize_on_edge : (int * int) -> u_tree -> u_tree

val cannonize_on_leaf : int -> u_tree -> u_tree

val exchange_codes : int -> int -> u_tree -> u_tree

(** A function to replace the codes on a tree from a given function. The
    function is called multiple times with the same code, thus, it should not,
    have side-effects or should take that into account. For example, although an
    identity function will operate properly, a hidden incremental reference will
    not, unless the input is memoized to be onto. *)
val replace_codes : (int -> int) -> u_tree -> u_tree

val destroy_component : int -> u_tree -> u_tree

val copy_component : int -> u_tree -> u_tree -> u_tree



(** {2 Finger Print Modules *)

(** Module to fingerprint trees and compare them *)
module Fingerprint : sig
    type t
    val fingerprint : u_tree -> t
    val compare : t -> t -> int
end

(** [CladeFP] implements fingerprinting of clades.  All clades with the same
    set of terminal nodes will have the same fingerprint. *)
module CladeFP :
sig
    type t

    type fp = All_sets.Integers.t

    module Ordered : Map.OrderedType with type t = fp

    (* A module of sets of clades *)
    module CladeSet : Set.S with type elt = fp

    (** [fpcompare] compares two fingerprints *)
    val fpcompare : fp -> fp -> int

    (** [calc tree] calculates the set of fingerprints for [tree] *)
    val calc : u_tree -> t

    (* [sets t] Returns the set of all the clades contained in [t]. *)
    val sets : u_tree -> CladeSet.t

    (** [query edge fps] retrieves the fingerprint for [edge] *)
    val query : edge -> t -> fp

    (** [return the number of leaves in a fingerprint *)
    val num_leaves : fp -> int

    (** fold over a tree *)
    val fold : (edge -> fp -> 'a -> 'a) -> t -> 'a -> 'a
end

module CladeFPMap : (Map.S with type key = CladeFP.fp)



(** {2 Tree Fusing Functions *)

type 'a fuse_locations = ('a * u_tree * edge) list Sexpr.t

val source_to_target : u_tree * int -> u_tree * int -> u_tree

val fuse_cladesize : min:int -> max:int -> 'a * int -> bool

val fuse_locations :
    ?filter:(CladeFP.fp * int -> bool) -> ('a * u_tree) list -> ('a * u_tree) -> 'a fuse_locations

val fuse_all_locations :
    ?filter:(CladeFP.fp * int -> bool) -> ('a * u_tree) list -> 'a fuse_locations

val fuse :
    source:u_tree * edge -> target:u_tree * edge -> u_tree



(** {2 Parse Module / IO *)

(** Define a Parsing module for trees; we also include some intermediary
    types of trees used when we have annotations on nodes, branch lengths on
    nodes, et cetera *)
module Parse : sig

    exception Illegal_argument
    exception Illegal_tree_format of string

    type 'a t = Leafp of 'a | Nodep of 'a t list * 'a

    type tree_types =
        | Flat of string t
        | Annotated of (tree_types * string)
        | Branches of ((string * float option) t)
        | Characters of ((string * string option) t)

    (** Strip the tree of an annotation; thus we have a flat, branch, or
        character result in the tree. *)
    val remove_annotations : tree_types -> tree_types

    (** [print_tree t] prints a tree **)
    val print_tree : string option -> tree_types -> unit

    (** [of_string x] given a string x of a tree in the form (a (b c)), returns
        its representation as an internal t type. If an error occurs, an
        Illegal_tree_format error is raised. *)
    val of_string : string -> tree_types list list

    (** [of_channel x] creates a list of all the trees contained in the input
        channel x in ascii format. Each tree should be in a single line. In case
        of error the function raises Illegal_molecular_format. *)
    val of_channel : in_channel -> tree_types list list

    (** [of_file] is a shortcut of [of_channel], when the channel is an opened
        file. *)
    val of_file : FileStream.f -> tree_types list list

    (** [stream_of_file f] produces a function that returns on each call one of
        the trees in the input file [f]. If no more trees are found, an
        End_of_file exception is raised. The function _requires_ that the trees
        be separated with associated information, semicolons, or stars.*)
    val stream_of_file : bool -> FileStream.f -> ((unit -> tree_types) * (unit -> unit))

    (** Order the nodes alphabetically *)
    val cannonic_order : tree_types -> tree_types

    (** [cleanup ~newroot f t] removes from the tree [t] the leaves that contain
        infromation [x] such that [f x = true]. If there is the need for a new
        root, then [newroot] must be provided. If the [newroot] is required and
        not provided, the function raises an [Illegal_argument] exception. *)
    val cleanup : ?newroot:string  -> (string -> bool) -> tree_types -> tree_types option

    (** [post_process t] takes a basic tree and converts it to Flat,Branches or
        Annotated based on it's properties. *)
    val post_process : (string * (float option * string option)) t * string -> tree_types

    (** Map on the node data of a tree *)
    val map : ('a -> 'b) -> 'a t -> 'b t

    (** Map a function on the taxa of the tree *)
    val map_tree : (string -> string) -> tree_types -> tree_types

    (** Strip the tree of all node annotations and extra information to a basic
        tree; when we only care about the topology *)
    val strip_tree : tree_types -> string t

    val maximize_tree : tree_types -> (string * string option) t

    (** Convert a tree to the u-tree interface. Many details are stipped from
        the conversion, and thus need to be dealt with in other ways. *)
    val convert_to : string option * tree_types list -> (string -> int) -> u_tree
end



(** {2 Functions to aid in debugging *)

(** return string representation of a node *)
val string_of_node : node -> string

(** return a string of a join-junction *)
val string_of_jxn : join_jxn -> string

(** Run a tree through a battery of consistency tests *)
val test_tree : u_tree -> unit

(** Print an edge; to aid debuging *)
val print_edge : edge -> unit

(** Print tree from id as the starting point. backward-in-order traversal *)
val print_tree : id -> u_tree -> unit

(** Generalization of above, applied to each handle id *)
val print_forest : u_tree -> unit




