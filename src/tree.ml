(* POY 5.0 Beta. A phylogenetic analysis program using Dynamic Homologies.    *)
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

let () = SadmanOutput.register "Tree" "$Revision: 3235 $"

exception Invalid_Node_Id of int
exception Invalid_Handle_Id
exception Invalid_Edge
exception Invalid_End_Nodes
exception Node_Is_Handle

let debug_print_edge_removal = false
let debug_print_edge_addition = false
let debug_print_node_removal = false
let debug_print_node_addition = false
let debug_tests = false
let debug_fusing = false
let debug_handle_of = false
let odebug = Status.user_message Status.Information

type id = int

type node =
    | Single of id
    | Leaf of id * id
    | Interior of (id * id * id * id)

type edge =
    | Edge of (id * id)

type t_status = Continue | Skip | Break

module EdgeComparator = struct
    type t = edge
    let compare x y =
        match x, y with
        | Edge(e11, e12), Edge(e21, e22) ->
            match (e11 - e21) with
            | 0 -> (e12 - e22)
            | x -> x
end


module EdgeSet = Set.Make(EdgeComparator)

module EdgeMap =  Map.Make(EdgeComparator)

type edge_delta = {
    added : edge list;
    removed : edge list;
}

let empty_ed = { added = []; removed = [] }

let merge_edged a b =
    { added = b.added @ a.added; removed = b.removed @ a.removed }

type side_delta =
    [   | `Single of int * bool
                (** l, l1, l2, handle added/removed *)
        | `Edge of int * int * int * int option
    ]

type join_jxn =
    | Single_Jxn of id
    | Edge_Jxn of id * id

let string_of_jxn = function
    | Single_Jxn id -> string_of_int id
    | Edge_Jxn (a, b) -> string_of_int a ^ " -- " ^ string_of_int b

let jxn_choose_node = function
    | Single_Jxn id -> id
    | Edge_Jxn (a, b) -> a

type reroot_delta = id list

type break_delta = side_delta * side_delta

type join_delta = side_delta * side_delta * reroot_delta

let join_to_edge_delta (sd1, sd2, rrd) =
    let a, r, n1 = match sd1 with
        | `Single (id, _) -> [], [], id
        | `Edge (l, l1, l2, _) ->
            [Edge(l, l1); Edge(l, l2);], [Edge(l1, l2);], l
    in
    let a, r, n2 = match sd2 with
        | `Single (id, _) -> a, r, id
        | `Edge (l, l1, l2, _) ->
            Edge(l, l1) :: Edge(l, l2) :: a, Edge(l1, l2) :: r, l
    in
    let a = Edge(n1, n2) :: a in
    { added = a; removed = r; }

let break_to_edge_delta (sd1, sd2) =
    let { added = a; removed = r; } = join_to_edge_delta (sd1, sd2, []) in
    { added = r; removed = a; }

(** [int_of_id id]
    @return the corresponding int of the id. *)
let int_of_id id = id

(** [print_break_jxn jxn]
    @return prints the break jxn. *)
let print_break_jxn (h1, h2) =
    let h1 = (int_of_id h1)
    and h2 = (int_of_id h2) in
    let h1_s = (string_of_int h1)
    and h2_s = (string_of_int h2) in
    print_string ("(" ^ h1_s ^ "," ^ h2_s ^ ")")

(** [string_of_node node]
    @return the node as a string. *)
let string_of_node = function
    | Single x           -> Printf.sprintf "S:%d" x
    | Leaf (x,y)         -> Printf.sprintf "L:%d->%d" x y
    | Interior (a,x,y,z) -> Printf.sprintf "I:%d->(%d,%d,%d)" a x y z

(** [print_join_1_jxn jxn]
    @return prints join_1_jxn. *)
let print_join_jxn jxn = match jxn with
    | Single_Jxn(s) ->
          print_string ("Single_1_Jxn " ^ (string_of_int s))
    | Edge_Jxn(e1, e2) ->
          let e1 = (int_of_id e1)
          and e2 = (int_of_id e2) in
          let e1_s = (string_of_int e1)
          and e2_s = (string_of_int e2) in
          print_string ("Edge_1_Jxn (" ^ e1_s ^ "," ^ e2_s ^ ")")


type u_tree = {
    tree_name : string option;
    u_topo  : node All_sets.IntegerMap.t;
    d_edges : EdgeSet.t;
    handles : All_sets.Integers.t;
    avail_ids : id list;
    new_ids : id;
}

let replace_codes f tree =
    let topo = 
        let fix_node = function
            | Single x -> Single (f x)
            | Leaf (x, y) -> Leaf ((f x), (f y))
            | Interior (a, b, c, d) -> Interior ((f a), (f b), (f c), (f d))
        in
        All_sets.IntegerMap.fold 
            (fun code node acc -> 
                All_sets.IntegerMap.add (f code) (fix_node node) acc)
            tree.u_topo
            All_sets.IntegerMap.empty
    and edges =
        let fix_edge (Edge (a, b)) = Edge ((f a), (f b)) in
        EdgeSet.fold 
            (fun edge acc -> EdgeSet.add (fix_edge edge) acc)
            tree.d_edges
            EdgeSet.empty
    and handles = 
        All_sets.Integers.fold
            (fun code acc -> All_sets.Integers.add (f code) acc) 
            tree.handles
            All_sets.Integers.empty
    and avail_ids = List.map f tree.avail_ids in
    { tree with
        u_topo = topo; d_edges = edges;
        handles = handles; avail_ids = avail_ids }

(** Return the next available code to apply to a tree. Check that the code is
    able to be a new node and continue. 

    TODO: Although, originally the codes should never be duplicated in avail_ids,
    an issue occured where duplication, and these checks were failing. I'm unsure
    why this is happening, it needs to be looked at. Although, it is not a major
    error and is probably due to join/break/build routines or re-coding routines
    being called in non-uniform ways. *)
let rec get_available tree = match tree.avail_ids with
    | id :: ids when List.mem id ids ->
        get_available {tree with avail_ids = ids; }
    | id :: ids when All_sets.IntegerMap.mem id tree.u_topo ->
        get_available {tree with avail_ids = ids; }
    | id :: ids ->
        id, { tree with avail_ids = ids }
    | [] ->
        assert (not (All_sets.IntegerMap.mem (succ tree.new_ids) tree.u_topo));
        tree.new_ids, { tree with new_ids = succ tree.new_ids }


type break_jxn = id * id

(** [get_id node]
    @param node the node whose id is desired
    @return id of the node. *)
let get_id node = match node with
    | Single(id) 
    | Leaf(id, _)
    | Interior(id, _, _, _) -> id


(** The empty unrooted tree. *)
let empty () = 
    {
        tree_name = None;
        u_topo = All_sets.IntegerMap.empty;
        d_edges = EdgeSet.empty;
        handles = All_sets.Integers.empty;
        avail_ids = [];
        new_ids = 0;
    }

(** [is_handle id tree]
    @param id node id.
    @param tree tree in which the node exists.
    @return true, when the node corresponding to id is a handle.
            false, otherwise. *)
let is_handle id tree =
    (All_sets.Integers.mem id tree.handles)

(** [is_edge e tree]
    @param e=Edge (x, y) tuple indicating end vertices x, y.
    @param tree tree in which Edge(x,y) exists or not.
    @return true if Edge(x, y) exists, else false. *)
let is_edge e tree =
    (EdgeSet.mem e tree.d_edges)

(**************** GET FUNCTIONS **********************************)

(** [get_node_ids tree]
    @return all the node ids of the tree in a list. *)
let get_node_ids tree =
    let first x y z = x::z in
    (All_sets.IntegerMap.fold first tree.u_topo [])

(** [get_nodes tree]
    @return all the nodes of the tree in a list. *)
let get_nodes tree =
    let second x y z = y::z in
    (All_sets.IntegerMap.fold second tree.u_topo [])

(** [get_handles tree]
    @return the handles of the tree. *)
let get_handles tree =
    tree.handles

(** [handle_list tree] returns a list of the handles *)
let handle_list tree =
    All_sets.Integers.elements tree.handles

(** [get_node_id id btree]
    @param id the id being verified
    @param tree the tree in which id is a valid node.
    @return returns the id with a [`Node] tag
    after verifying it.
    @raise Invalid_Node_Id when the id does not have a node associated
    with it. *)
let get_node_id id tree = match (All_sets.IntegerMap.mem id tree.u_topo) with
    | true -> id
    | false -> raise (Invalid_Node_Id id)

(** [get_handle_id id btree]
    @param id the id being verified
    @param tree the tree in which id is a valid handle.
    @return returns the id with a [`Handle] tag
    after verifying it.
    @raise Invalid_Handle_Id when the id does not have a node associated
    with it. *)
let get_handle_id id tree = match (All_sets.Integers.mem id tree.handles) with
    | true -> id
    | false -> raise Invalid_Handle_Id

(** [get_node id btree]
    @param id id of the node desired.
    @param tree tree to which the node belongs
    @return returns the node associated with the id.
    @raise Invalid_Node_Id, if the id is not a valid id. *)
let get_node id tree =
    try (All_sets.IntegerMap.find id tree.u_topo)
    with | Not_found -> raise (Invalid_Node_Id id)

(** [get_leaves ?acc tree handle] returns a list of the leaves in [tree]
    below node [handle] *)
let rec get_leaves ?(acc=[]) tree handle =
    if is_handle handle tree
    then begin match get_node handle tree with
    | Leaf (id, par) -> get_leaves ~acc:(id :: acc) tree par
    | Single id -> id :: acc
    | Interior (id, p, a, b) ->
          let acc = get_leaves ~acc tree p in
          let acc = get_leaves ~acc tree a in
          let acc = get_leaves ~acc tree b in
          acc
    end else begin match get_node handle tree with
    | Leaf (id, _) -> id :: acc
    | Single id -> failwith "Tree.get_leaves"
    | Interior (id, _, a, b) ->
          let acc = get_leaves ~acc tree a in
          let acc = get_leaves ~acc tree b in
          acc
    end

(** [get_all_leaves tree] returns a list of all the leaves in the tree *)
let get_all_leaves tree =
    let fn handle acc = get_leaves ~acc tree handle in
    All_sets.Integers.fold fn tree.handles []

(** [is_leaf id tree]
    @param id the id of the node.
    @param tree the tree in which the id corresponds to  the leaf or not.
    @return whether the node corresponding to the id is a leaf or not. *)
let is_leaf id tree = match get_node id tree with
    | Leaf(_, _) -> true
    | _ -> false

let is_single id tree = match get_node id tree with
    | Single _ -> true
    | _ -> false

(** [normalize_edge edge tree] returns either [edge] or the reverse-direction
    [edge'] if one of those is in the tree.  Otherwise, it fails with exception
    [Invalid_Edge]. *)
let normalize_edge e tree =
    let Edge (x,y) = e in
    if EdgeSet.mem e tree.d_edges
        then e
        else begin
            let e' = Edge (y,x) in
            if EdgeSet.mem e' tree.d_edges
                then e'
                else 
                    let () = prerr_endline (string_of_int x ^ "-" ^ string_of_int y) in
                    raise Invalid_Edge
        end

(** [get_edge edge btree]
    @param edge the edge as a tuple (x,y).
    @param btree the tree.
    @return the edge Edge(x,y) after verifying the tuple.
    @raise Invalid edge when the edge does not belong to the tree. *)
let get_edge (x, y) btree =
    normalize_edge (Edge (x, y)) btree

(** [nbr_of node1 node2]
    @param node1 a node
    @param node2 another node
    @return true if node2 is a neighbor of node1. *)
let nbr_of node1 node2 = match node1, node2 with
    | Leaf(id1, nbr1), Leaf(id2, nbr2) ->
        if (id1 = nbr2) && (id2 = nbr1) then
            true
        else begin
            Printf.printf "The neighbors of %d and %d are %d and %d"
            id1 id2 nbr1 nbr2;
            false
        end;
    | Leaf(id1, nbr1), Interior(id2, nbr21, nbr22, nbr23)
    | Interior(id2, nbr21, nbr22, nbr23), Leaf(id1, nbr1) ->
        if (id2 = nbr1) && (List.mem id1 [nbr21; nbr22; nbr23]) then
            true
        else begin
            Printf.printf "The neighbors of %d and %d are %d %d %d and \
            %d\n%!" id2 id1 nbr21 nbr22 nbr23 nbr1;
            false
        end
    | Interior(id1, nbr11, nbr12, nbr13),
      Interior(id2, nbr21, nbr22, nbr23) ->
        if ((List.mem id1 [nbr21; nbr22; nbr23]) &&
            (List.mem id2 [nbr11; nbr12; nbr13])) then
            true
        else begin
            Printf.printf "The neighbors of %d and %d are %d %d %d and %d %d %d\n%!"
                          id1 id2 nbr11 nbr12 nbr13 nbr21 nbr22 nbr23;
            false
        end
    | _ -> raise Invalid_End_Nodes

(** [verify_edge edge btree]
    @param edge the edge which is being verified i.e. whether it
                it indeed exists in the tree or not.
    @param btree the tree in which the edge exists or not.
    @return true if the edge exists in the tree, false otherwise. *)
let verify_edge edge btree =
    let Edge(e1, e2) = edge in
    let node1 = (get_node e1 btree) in
    let node2 = (get_node e2 btree) in
    let edge = normalize_edge edge btree in
    (EdgeSet.mem edge btree.d_edges) && (nbr_of node1 node2)

(**************** BEGIN - UNSAFE OPS ***************************************)

(** The following functions are unsafe i.e. they allow us to modify the
structure of the tree in a manner that does not maintain the invariants
of the tree data-structure. These functions should never be exposed to the
outside user. *)

(************************ ADD FUNCTIONS *************************)

(** [add_node node btree]
    @param node the node that needs to be added to the tree.
    @param tree the tree to which the node is being added.
    @return the tree with the node added.
    This function does not check whether the node is already
    present in the tree or not. If the node is already present then
    the tree is returned unmodified. *)
let add_node node btree =
    let id = (get_id node) in
    if debug_print_node_addition then
        Printf.fprintf stderr "Adding node %d\n%!" id;
    assert (if All_sets.IntegerMap.mem id btree.u_topo then begin print_int id;
    print_newline (); false end else true);
    let avail_ids = List.filter (fun i -> i <> id) btree.avail_ids in
    let new_topo = (All_sets.IntegerMap.add id node btree.u_topo) in
    { btree with u_topo = new_topo; avail_ids = avail_ids; }

(** [add_nodes nodes btree]
    @param nodes list of nodes that need to be added to the tree.
    @param btree tree to which the nodes are being added.
    @return tree with the list of nodes added. Older nodes with the
    same ids will be overwritten. *)
let add_nodes nodes btree =
    let f btree node = (add_node node btree) in
        (List.fold_left f btree nodes)

(** [add_edge edge btree]
    @param edge the edge that needs to be added to the tree.
    @param btree the tree to which the edge is being added.
    @return tree with the given edge added. If the edge already
    exists, tree is returned as is. *)
let add_edge edge btree =
    if debug_print_edge_addition then begin
        let Edge (e1, e2) = edge in
        Printf.fprintf stderr "Adding edge %d %d.\n%!" e1 e2;
    end;
    let new_edges = (EdgeSet.add edge btree.d_edges) in
        { btree with d_edges = new_edges }

(** [add_edges edges btree]
    @param edges the edges that need to be added to the tree.
    @param btree the tree to which edges need to be added.
    @return a new tree with the edges added. *)
let add_edges edges btree =
    let f btree edge = (add_edge edge btree) in
    List.fold_left f btree edges

(** [add_handle id btree]
    @param id node id to be added as a handle.
    @param btree tree to which the node is being added as a handle.
    @return tree with the node added as a handle. Any previous node
    or handle with the id of node is overwritten. *)
let add_handle id btree =
    let new_handles = (All_sets.Integers.add id btree.handles) in
    { btree with handles = new_handles }

(** [add_handles nodes btree]
    @param nodes node ids to be added as handles.
    @param btree tree to which the nodes are added as handles.
    @return btree with nodes added as handles. *)
let add_handles handles btree =
    let f btree handle = (add_handle handle btree) in
    List.fold_left f btree handles

(************************ REMOVE FUNCTIONS ****************************)

(** [get_node id btree]
    @param id the id of the node that belongs to the tree.
    @param tree tree from whch the node with the given id is
                removed.
    @return tree with the node removed.
    This operation only removes nodes that are not handles.
    To remove handles, use the remove_handles function.
    This operation does not check whether the node is actually
    present in the tree or not. *)
let remove_node id btree =
    assert (All_sets.IntegerMap.mem id btree.u_topo);
    if debug_print_node_removal then
        Printf.fprintf stderr "Removing node %d\n%!" id;
    let new_topo = (All_sets.IntegerMap.remove id btree.u_topo) in
    assert (not (All_sets.IntegerMap.mem id new_topo));
    { btree with u_topo = new_topo; avail_ids = id :: btree.avail_ids }

(* [remove_reserved_node id btree] 
 * Removes a node from a tree, without making it an available node for other
 * tree operations in [avail_ids].
 * @param id is the node to be removed but reserved.
 * @param btree is the tree that holds the node
 * @return the tree with the node removed but the available ids unmodified. *)
let remove_reserved_node id btree =
    let new_topo = All_sets.IntegerMap.remove id btree.u_topo in
    { btree with u_topo = new_topo }

(** [remove_nodes ids btree]
    @param ids the list of ids of nodes that belong to the tree.
    @param btree tree from which nodes with ids given the ids list
    will be removed from.
    @return tree with nodes with id in ids list removed. *)
let remove_nodes ids btree =
    let f btree id =  (remove_node id btree) in
        (List.fold_left f btree ids)

(** [remove_handle id btree]
    @param id id of the handle.
    @param btree tree in which the node corresponding to id is a handle.
    @return tree with handle removed. *)
let remove_handle id btree =
    let new_handles = (All_sets.Integers.remove id btree.handles) in
    { btree with handles = new_handles }

(** [remove_handles ids btree]
    @param ids of the handles to be removed.
    @param btree the tree whose handles are to be removed.
    @return tree with the handles removed. *)
let remove_handles ids btree =
    let f btree handle = (remove_handle handle btree) in
    List.fold_left f btree ids

(** [remove_edge edge btree]
    @param edge edge to be removed
    @param btree tree from which the edge is being removed.
    @return tree with the edge removed. *)
let remove_edge edge btree =
    if debug_print_edge_removal then begin
        let Edge (e1, e2) = edge in
        Printf.printf "Removing edge %d %d.\n%!" e1 e2;
    end;
    let new_d_edges = (EdgeSet.remove edge btree.d_edges) in
        { btree with d_edges = new_d_edges }

(** [remove_edges edges btree]
    @param edges the edges to be removed.
    @param btree tree from which the edges are being removed.
    @return tree with the edges removed. *)
let remove_edges edges btree =
    let f btree edge = (remove_edge edge btree) in
        (List.fold_left f btree edges)

(** [remove_edge_either edge tree] removes either the edge given or its reverse,
    depending on which exists in the tree.
    @raise Invalid_edge if neither edge exists *)
let remove_edge_either e t = remove_edge (normalize_edge e t) t

(** [remove_incident_edges id tree] returns a tree with all the edges coming
    into or out of node [id] removed.
let remove_incident_edges id tree =
    let rem a b t =
        if verify_edge (Edge (a, b)) t
        then remove_edge_either (Edge (a, b)) t
        else t in
    match get_node id tree with
    | Leaf (id, nbr) -> rem id nbr tree
    | Interior (id, n1, n2, n3) ->
          rem id n1
              (rem id n2
                   (rem id n3 tree))
    | Single _ -> tree *)

(********************** OTHER FNS *************************************)

(** [replace_nbr nbr new_nbr node]
    @return replace the nbr with new_nbr and returns the new node. *)
let replace_nbr nbr new_nbr node = match node with
    | Leaf(id1, nbr1) when (nbr1 = nbr) -> Leaf(id1, new_nbr)
    | Interior(id1, nbr1, nbr2, nbr3) when (nbr1 = nbr) ->
        Interior(id1, new_nbr, nbr2, nbr3)
    | Interior(id1, nbr1, nbr2, nbr3) when (nbr2 = nbr) ->
        Interior(id1, nbr1, new_nbr, nbr3)
    | Interior(id1, nbr1, nbr2, nbr3) when (nbr3 = nbr) ->
        Interior(id1, nbr1, nbr2, new_nbr)
    | _ -> raise
            (Invalid_argument "Single node or Incorrect nbr passed in")

(** [replace_nbr_in_tree nbr new_nbr id tree]
    @return returns the tree with node(id) changed. nbr of node(id) is
            replaced with new_nbr. *)
let replace_nbr_in_tree nbr new_nbr id btree =
    let node = (get_node id btree) in
    let new_node = (replace_nbr nbr new_nbr node) in
    let btree = remove_node id btree in
    let btree = add_node new_node btree in
    btree

(** [get_parent id tree]
    @return the nbr, such that Edge(nbr, id) exists in tree. *)
let get_parent id tree = match get_node id tree with
    | Leaf(a, nb) -> nb
    | Interior(a, p, l, r) -> p
    | _ -> raise (Invalid_argument "Single Handle passed in")

(** [test_component tree handle] performs consistency tests on a component of a
    tree. *)
let test_component tree handle =
    let test_parent e1 e2 =
        let parent = get_parent e2 tree in
        assert (if not (parent = e1) then 
            Printf.fprintf stderr "Failed for %d, %d, %d\n%!" parent e1 e2;
        parent = e1); () in
    let test_edge e1 e2 =
        assert (let res = is_edge (Edge (e1, e2)) tree in
            if not res then 
            Printf.fprintf stderr "Failed on %d %d\n%!" e1 e2; 
            res);
        test_parent e1 e2 in
    let test_node e1 = match get_node e1 tree with
    | Single id -> assert (id = e1); ()
    | Leaf (id, nbr) -> assert (nbr <> id); assert (id = e1); ()
    | Interior (id, n1, n2, n3) ->
          assert (id = e1);
          assert (id <> n1);
          assert (id <> n2);
          assert (id <> n3);
          assert (n1 <> n2);
          assert (n1 <> n3);
          assert (n2 <> n3); () in
    let rec r prev e1 = match get_node e1 tree with
    | Single _ -> assert false
    | Leaf (id, n) ->
          assert (prev = n);
          assert (All_sets.IntegerMap.mem n tree.u_topo);
          test_node e1;
          test_edge n id
    | Interior (id, n1, n2, n3) ->
          assert (All_sets.IntegerMap.mem n1 tree.u_topo);
          assert (All_sets.IntegerMap.mem n2 tree.u_topo);
          assert (All_sets.IntegerMap.mem n3 tree.u_topo);
          assert (if not (prev = n1) then 
            Printf.fprintf stderr "Failed test with prev:%d and n1=%d in %d"
            prev n1 id;
            prev = n1);
          test_node e1;
          test_edge n1 id;
          r id n2;
          r id n3 in
    match get_node handle tree with
    | Single id -> test_node handle
    | Leaf (id, nbr) -> test_node handle; r handle nbr
    | Interior (id, n1, n2, n3) ->
          assert (handle = id);
          test_node handle;
          r handle n1;
          r handle n2;
          r handle n3

(** [test_tree tree] performs consistency tests on a tree *)
let test_tree tree =
    All_sets.Integers.iter (test_component tree) tree.handles

(******************* END - UNSAFE OPS ****************************************)

(** [a --> b] performs [b a], allowing for chaining of tree-modifying functions *)
let (-->) a b = b a

(** [make_disjoint_tree n]
    @return a disjointed tree with the given nodes and 0 edges *)
let make_disjoint_tree nodes =
    let nodes = List.sort (fun a b -> b - a) nodes in
    let tree = 
        let t = empty () in
        let next_avail = 
            match nodes with
            | [] -> 0
            | h :: _ -> h + 1
        in
        { t with new_ids = next_avail }
    in
    let f tr m =
        let tr = (add_node (Single(m)) tr) in
        let tr = (add_handle m tr) in
            tr
    in
    (List.fold_left f tree nodes)

(** [handle_of n tree] returns the handle of the component in [tree] containing
    node [n] *)
let rec handle_of n tree =
    if debug_handle_of then odebug ("handle of: " ^ string_of_int n);
    if is_handle n tree then n else handle_of (get_parent n tree) tree

(** [get_path_to_handle id tree]
    @return this function returns the path from node to the handle.
            The path is a list of directed edges. *)
let get_path_to_handle id tree =
    let rec aux id elist =
        if (is_handle id tree) then
            id, elist
        else
            let nbr = (get_parent id tree) in
            let edge = (Edge(nbr, id)) in
                (aux nbr (edge :: elist))
    in
    aux id []

(** [get_vertices_to_handle id tree]
    @return A list of vertices from an [id] to the handle that tree contains *)
let get_vertices_to_handle id tree =
    let rec aux id acc =
        if is_handle id tree then 
            List.rev (id :: acc)
        else 
            let nbr = get_parent id tree in
            aux nbr (id :: acc)
    in
    aux id []

(** [path_up ?acc nfrom nto tree] returns the path from [nfrom] towards the
    root to [nto]; the path is returned in reverse order, with [nto] first *)
let rec path_up ?(acc=[]) nfrom nto tree =
    try if nfrom = nto then nto :: acc
        else path_up ~acc:(nfrom :: acc) (get_parent nfrom tree) nto tree
    with | err ->
        Status.user_message Status.Error
            ("nfrom: " ^ string_of_int nfrom ^ " nto: " ^ string_of_int nto);
        raise err

(** [reorient_node parent child tree] changes the !Interior structore of
    [child] so that [parent] is its parent (i.e., is in the first position) *)
let reorient_node parent child tree =
    match get_node child tree with
    | Interior (node, c1, c2, c3) ->
        if parent = c1 then
            tree
        else if parent = c2 then begin
            let tree = remove_node child tree in
            let tree = add_node (Interior (child, parent, c1, c3)) tree in
            tree
        end else if parent = c3 then begin
            let tree = remove_node child tree in
            let tree = add_node (Interior (child, parent, c1, c2)) tree in
            tree
        end else 
            failwith "reorient_node"
    | Leaf (node, c) when parent = c ->
        tree
    | Leaf _
    | Single _ -> failwith "reorient_node"

(** [fix_path_to_handle tree list] changes nodes and reorients edges to provide
    consistency when moving a handle *)
let rec fix_path_to_handle tree list = match list with
    | [] -> tree
    | [new_handle] -> tree
    | p :: c :: list ->
        let tree = 
            tree --> remove_edge (Edge (p, c))
                 --> add_edge (Edge (c, p))
                 --> reorient_node c p
        in
        fix_path_to_handle tree (c :: list)

(** Function to move the node handle from one node to another.
    This will change the directions of the edges to account
    for this change. *)
let move_handle node_id tree =
    let handle = handle_of node_id tree in
    let path = path_up node_id handle tree in
    let tree = fix_path_to_handle tree path in
    let tree = tree --> remove_handle handle --> add_handle node_id in
    if debug_tests then test_tree tree;
    tree, path

(** [other_two_nbrs nbr int_node]
    @return returns the other two nbrs of the interior node. *)
let other_two_nbrs nbr node =
    match node with
    | Interior (id, nbr1, nbr2, nbr3) ->
            if (nbr = nbr1) then (nbr2, nbr3)
            else if (nbr = nbr2) then (nbr1, nbr3)
            else begin
                assert (
                    if nbr = nbr3 then true
                    else begin
                        let mst = 
                            Printf.sprintf "This is the failure point: %d but the \
                            neighbors of %d are %d %d and %d" nbr id nbr1 nbr2
                            nbr3 
                        in
                        Status.user_message Status.Error mst;
                        false
                    end);
                (nbr1, nbr2)
            end
    | _ -> raise (Invalid_argument "Not interior node")

(** [break edge btree] breaks the tree at the jxn.
    @param edge the junction {!Tree.break_jxn} where the tree should
    be broken.
    @param btree the binary tree whose edge is being broken.
    @return an {!Tree.aux_tree}, the new binary tree with the edge
            removed and the changes to the difference between the old
            and new trees as a {!Tree.tree_delta} structure. *)
let break jxn tree =
    let (e1, e2) = jxn in
    assert(verify_edge (Edge (e1, e2)) tree);
    (* Check whether root is on left side or right side *)
    let handle_left = is_edge (Edge (e1, e2)) tree in
    (* Handle the left side *)
    let node1 = get_node e1 tree in
    let node2 = get_node e2 tree in
    let handle_side handle_our_side node other_vertex tree = match node with
        | Leaf (id, _) ->
            (* if it's a leaf: make it a single *)
            let tree = tree --> remove_node id --> add_node (Single id) in
            (* add a handle if it's not already *)
            let tree =
                if not handle_our_side
                then add_handles [id] tree
                else tree
            in
            tree, id, Single_Jxn id, `Single (id, not handle_our_side)
        | Interior (id, n1, n2, n3) ->
            let l = id in
            let l1, l2 = other_two_nbrs other_vertex node in
            let new_root = not handle_our_side || (is_handle l tree) in
            let tree, jxn, root =
                if new_root then begin
                    let tree = if is_handle l tree 
                        then remove_handle l tree
                        else tree
                    in
                    let tree = 
                        tree --> add_handle l1
                             --> add_edge (Edge (l1, l2))
                             --> replace_nbr_in_tree l l2 l1
                             --> replace_nbr_in_tree l l1 l2
                             --> reorient_node l1 l2
                             --> reorient_node l2 l1
                    in
                    tree, Edge_Jxn (l1, l2),Some l1
                end else begin
                    let li, lj = if is_edge (Edge (l1, l)) tree 
                        then l1, l2
                        else l2, l1
                    in
                    let tree = 
                        tree --> add_edge (Edge (li, lj))
                             --> replace_nbr_in_tree l li lj
                             --> replace_nbr_in_tree l lj li
                             --> reorient_node li lj
                    in
                    tree, Edge_Jxn (li, lj), None
                end
            in
            let tree = 
                tree --> remove_edge_either (Edge(l, l1))
                     --> remove_edge_either (Edge(l, l2))
                     --> remove_node l
            in
            tree, id, jxn, `Edge (l, l1, l2, root)
        | Single _ -> assert false
    in
    let tree,t1_handle,t1_jxn,ld = handle_side handle_left node1 e2 tree in
    let tree,t2_handle,t2_jxn,rd = handle_side (not handle_left) node2 e1 tree in
    let tree = remove_edge_either (Edge(e1, e2)) tree in
    if debug_tests then test_tree tree;
    tree, (ld, rd)

(** [join jxn1 jxn2 atree]
    @param jxn1 jxn in the first tree to which the second tree will
                be grafted onto.
    @param jxn2 jxn in the second tree which will be grafted to a jxn
                in the first tree.
    @param atree aux_tree.
    @return an aux tree with the two smaller trees grafted onto a larger
            tree. *)
let join jxn1 jxn2 tree =
    let debug = false in
    if debug then Printf.printf "tree.join\n%!";
    test_tree tree;
    (* The junctions should come from different components *)
    assert (handle_of (jxn_choose_node jxn1) tree
            <> handle_of (jxn_choose_node jxn2) tree);
    (* [split node edge tree] *)
    let split ex_nd edge bt =
        let Edge (e1, e2) = edge in
        bt --> replace_nbr_in_tree e1 ex_nd e2
           --> replace_nbr_in_tree e2 ex_nd e1
           --> remove_edge edge
           --> add_edge (Edge (e1, ex_nd))
           --> add_edge (Edge (ex_nd, e2))
    in
    (* changes and return tree, left new node id, and a function to make
       the actual node given its neighbor's id *)
    let tree, lid, lnode, ldel = match jxn1 with
        | Single_Jxn id ->
                tree --> remove_reserved_node id, id, (fun nbr -> Leaf (id, nbr)), 
                `Single (id, false)
        | Edge_Jxn (e1, e2) ->
                let Edge (e1, e2) = normalize_edge (Edge (e1, e2)) tree in
                let n, tree = get_available tree in
                let tree = split n (Edge (e1, e2)) tree in
                tree, n, (fun nbr -> Interior (n, e1, e2, nbr)), `Edge (n, e1, e2, None) 
    in
    (* make changes and return tree, right new node id, a function to make the
       actual node given its neighbor's id, and a list of nodes that have to be
       changed in direction *)
    let tree, rid, rnode, rdel, rpath = match jxn2 with
        | Single_Jxn id ->
                let tree = tree --> remove_handle id --> remove_reserved_node id in
                tree, id, (fun nbr -> Leaf (id, nbr)), `Single (id, true), [id]
        | Edge_Jxn (e1, e2) ->
                (* This case is the most difficult, as it puts the tree on the RHS
                   into a bad state.  We proceed as follows.   *)
                (* make edge e1->e2 so that e1 is closer to the root *)
                let Edge (e1, e2) = normalize_edge (Edge (e1, e2)) tree in
                let handle = handle_of e1 tree in
                let n, tree = get_available tree in
                let tree = remove_handle handle tree in
                (* fix nodes along the path *)
                (* (this returns path with handle as the first element) *)
                let path = 
                    try path_up ~acc:[n] e1 handle tree 
                    with err ->
                        Status.user_message Status.Error
                            (Printf.sprintf "The error is with new node %d\n%!" n);
                        raise err
                in
                let Edge (e1, e2) = normalize_edge (Edge (e1, e2)) tree in
                let tree = split n (Edge (e1, e2)) tree in
                let tree = tree --> reorient_node n e1 --> reorient_node n e2 in
                let tree = fix_path_to_handle tree path in
                tree, n, (fun nbr -> Interior (n, nbr, e1, e2)),
                    `Edge (n, e1, e2, Some handle), path 
    in
    let n1 = lnode rid in
    let n2 = rnode lid in
    let tree = 
        tree --> add_edge (Edge (lid, rid)) --> add_node n1 --> add_node n2 
    in
    let path = rpath in
    if debug_tests then test_tree tree;
    tree, (ldel, rdel, path)

(** [random t]
    @param t - nodes to create new random tree from
    @return  - a uniformally random tree *)
let random (nodes : int list)  : u_tree =
    let random_edge t =
        let ray = Array.of_list (EdgeSet.elements t.d_edges) in
        let i = Random.int (Array.length ray) in
        ray.(i)
    in
    let add_node t n =
        let (Edge (a,b)) = random_edge t in
        let t,_ = join (Single_Jxn n) (Edge_Jxn (a,b)) t in
        t
    in
    let tree = make_disjoint_tree nodes in
    match nodes with
    | a :: b :: t ->
        let acc,_ = join (Single_Jxn a) (Single_Jxn b) tree in
        List.fold_left (add_node) acc t
    | _ :: [] | [] -> tree


(** [pre_order_edge_visit f id ptree accum]
@param f function to be applied to each edge of the subtree
         of ptree rooted at id.
@param id all the edges that are descendents of this node are visited
          in pre-order manner.
@param ptree the ptree whose edges are being visited.
@param accum the accumulator, a list.
@return the return values of the function being applied on the edges. *)
let pre_order_edge_visit f id bt accum =
    let break_or_continue pred nd ret c_with s_flag =
        let status, accum = ret in
        match status with
        | Continue -> (c_with pred nd accum)
        | Break    -> (Break, accum)
        | Skip     ->
              begin match s_flag with
                  | true -> (c_with pred nd accum)
                  | false -> (Continue, accum)
              end
    in
    let visit_node pred nd acc = (f (Edge (pred,nd)) acc) in
    let rec process_children pred nd accum =
        let node = (get_node nd bt) in
        let (x, y) = other_two_nbrs pred node in
        let x, y = match Random.int 2 with
            | 0 -> x, y
            | _ -> y, x
        in
        let ret = (traverse nd x accum) in
        (break_or_continue nd y ret traverse true)
    and traverse pred id accum = match (get_node id bt) with
        | Leaf(nd, nbr) ->
              assert(pred = nbr);
              (visit_node pred id accum)
        | Interior(nd, nbr1, nbr2, nbr3) ->
              let ret = (visit_node pred nd accum) in
              (break_or_continue pred nd ret process_children false)
        | _ -> raise (Invalid_argument "poev 2: Check node type")
    in
    match (get_node id bt) with
    | Leaf(nd, nbr) ->
          assert( is_edge (Edge (nd, nbr)) bt);
          let _, accum = (traverse nd nbr accum) in
          accum
    | Interior(nd, nbr1, nbr2, nbr3) ->
          begin
              match (is_handle nd bt) with
              | true ->
                    let ret = (traverse nd nbr1 accum) in
                    let _, accum =
                        (break_or_continue nbr1 nd ret process_children true) in
                    accum
              | false ->
                    let pred = (get_parent nd bt) in
                    let _, accum = (process_children pred nd accum) in
                    accum
          end
    | Single(_) -> accum

(** [pre_order_node_visit f id bt ad acc]
    @param f function to applied to all the nodes in pre-order.
    @param id the node_id from where the traversal is started.
    @param bt the tree whose nodes are visited.
    @param acc the result of the function application on the edges.
    @return the function applied to the subtree of id as imposed by the
        handle in the component, the results are returned as a list. *)
let pre_order_node_visit f id bt accum =
    let visit_node pred nd acc =
        (f pred nd acc)
    in
    let break_or_continue pred nd ret c_with s_flag =
        let status, accum = ret in
        match status with
        | Continue ->
              (c_with pred nd accum)
        | Break -> (Break, accum)
        | Skip ->
              begin
                  match s_flag with
                  | true -> (c_with pred nd accum)
                  | false -> (Continue, accum)
              end
    in
    let rec process_children pred nd acc =
        match pred with
        | Some p ->
              let node = (get_node nd bt) in
              let (x, y) = (other_two_nbrs p node) in
              let ret = (traverse (Some nd) x acc) in
              (break_or_continue (Some nd) y ret traverse true)
        | _ -> raise (Invalid_argument "Needs valid predecessor.")
    and traverse pred id accum =
        match (get_node id bt) with
        | Leaf(nd, nbr) ->
              begin
                  match pred with
                  | None -> assert(is_handle id bt)
                  | Some(p) -> assert(p = nbr)
              end;
              (visit_node pred nd accum)
        | Interior(nd, nbr1, nbr2, nbr3) ->
              let ret = (visit_node pred nd accum) in
              (break_or_continue pred nd ret process_children false)
        | _ -> raise (Invalid_argument "ponv 1: Check node type")
    in
    match (get_node id bt) with
    | Leaf(nd, nbr) ->
        (* if the node is a leaf and the edge is pointing into the node,
         * then we just visit the one node and stop. *)
          let ret = (visit_node None nd accum) in
          if (is_edge (Edge (nd, nbr)) bt) then
              let _, acc =
                  (break_or_continue (Some nd) nbr ret traverse false) in
              acc
          else
              let _, acc = ret in
              acc
    | Interior(nd, nbr1, nbr2, nbr3) ->
          begin
              match (is_handle nd bt) with
              | true ->
                    let ret = (visit_node None nd accum) in
                    let ret =
                        (break_or_continue (Some nd) nbr1 ret traverse true) in
                    let _, acc =
                        (break_or_continue (Some nbr1) nd
                             ret process_children true) in
                    acc
              | false ->
                    let pred = (get_parent nd bt) in
                    let ret = (visit_node (Some pred) nd accum) in
                    let _, acc =
                        (break_or_continue (Some pred) nd
                             ret process_children false) in
                    acc
          end
    | Single(_) -> accum

(** [post_order_node_with_edge_visit] A function to traverse a tree from an
    edge, applying a function on this edge in each direction, applying [f] if
    the node is a leaf, or g if the node is interior, with two accumulators, one
    from the left and one from the right. Single nodes are ignored. *)
let post_order_node_with_edge_visit f g (Edge (a, b)) bt accum =
    let rec processor prev curr accum =
        match get_node curr bt with
        | Leaf (nd, nbr) ->
                f prev curr accum
        | Interior (nd, nbr1, nbr2, nbr3) as node ->
                let a, b = other_two_nbrs prev node in
                let aacc = processor nd a accum
                and bacc = processor nd b accum in
                g prev curr aacc bacc
        | Single _ -> accum
    in
    let a = processor b a accum
    and b = processor a b accum in
    a, b

(** [post_order_node_with_edge_visit_simple f e t a] is a simplified visitor
    function [f], which is applied on every (non single) vertex, starting in 
    the (hypothetical) root located between the two vertices of edge [e], over
    tree [t] with accumulator [a]. *)
let post_order_node_with_edge_visit_simple f (Edge (a, b)) bt acc =
    let rec processor prev curr acc =
        match get_node curr bt with
        | Leaf (nd, nbr) ->
                f prev curr acc
        | Interior (nd, nbr1, nbr2, nbr3) as node ->
                let a, b = other_two_nbrs prev node in
                acc --> processor curr a 
                    --> processor curr b 
                    --> f prev curr
        | Single _ -> acc
    in
    acc --> processor b a --> processor a b

(** [pre_order_node_with_edge_visit_simple f e t a] is a simplified visitor
    function [f], which is applied on every (non single) vertex, starting in 
    the (hypothetical) root located between the two vertices of edge [e], over
    tree [t] with accumulator [a]. *)
let pre_order_node_with_edge_visit_simple f (Edge (a, b)) bt acc =
    let rec processor prev curr acc =
        match get_node curr bt with
        | Leaf (nd, nbr) ->
                f prev curr acc
        | Interior (nd, nbr1, nbr2, nbr3) as node ->
                let a, b = other_two_nbrs prev node in
                acc 
                    --> f prev curr
                    --> processor curr a 
                    --> processor curr b 
        | Single _ -> acc
    in
    acc --> processor b a --> processor a b

(** [pre_order_node_with_edge_visit_simple_root rf f e t a] is a simplified
    visitor function [f], which is applied on every (non single) vertex, with a
    special function [rf] applied to the (hypothetical) root location between
    two vertices of edge [e], over tree [t] with accumulator [a]. *)
let pre_order_node_with_edge_visit_simple_root f (Edge (a, b)) bt acc =
    let rec processor prev curr acc = match get_node curr bt with
        | Leaf (nd, nbr) ->
                f prev curr acc
        | Interior (nd, nbr1, nbr2, nbr3) as node ->
                let a, b = other_two_nbrs prev node in
                acc --> f prev curr
                    --> processor curr a
                    --> processor curr b
        | Single _ -> acc
    and processor_skip prev curr acc = match get_node curr bt with
        | Leaf (nd, nbr) -> acc
        | Interior (nd, nbr1, nbr2, nbr3) as node ->
                let a, b = other_two_nbrs prev node in
                acc --> processor curr a
                    --> processor curr b
        | Single _ -> acc
    in
    acc --> processor_skip b a --> processor_skip a b


(** [create_partition t l r]; returns the partitions that that seperate the edge
 * between the [l] and [r] nodes using the tree topology [t]. *)
let create_partition tree edge = 
    post_order_node_with_edge_visit
        (fun p x a -> All_sets.Integers.add x a)
        (fun p x acc_l acc_r -> All_sets.Integers.union acc_l acc_r)
        edge
        tree
        All_sets.Integers.empty


(** determine the Robinson-Foulds distance between two trees. **)
let robinson_foulds tree1 tree2 : int =
    (* compare two sets; or complement for partition equality. *)
    let compare_set all a b =
        let other_b = All_sets.Integers.diff all b in
        let equality_one = All_sets.Integers.compare a b
        and equality_two = All_sets.Integers.compare a other_b in
        (equality_one = 0) || (equality_two = 0) in
    (* ensure partitions are unique; removes extra partition on root *)
    let unique all lst =
        let rec unique acc = function
            | [] -> acc
            | hd::tl ->
                if List.exists (compare_set all hd) acc 
                    then unique acc tl
                    else unique (hd::acc) tl
        in
        unique [] lst in
    (* generate list of all partitions in a tree *)
    let partitions_of_tree all tree =
        let lst = 
            EdgeSet.fold
                (fun e acc ->
                    (fst (create_partition tree e)) :: acc)
                tree.d_edges
                []
        in
        unique all lst in
    (* return set of all from a tree; we also assert one handle *)
    let get_all tree = 
        List.fold_left
            (fun acc x -> All_sets.Integers.add x acc)
            (All_sets.Integers.empty)
            (get_all_leaves tree) in
    (* remove partitions of one list from the other *)
    let compare_partitions all p1 p2 = 
        let remove_partition p ps = List.filter (fun x -> not (compare_set all p x)) ps in
        let p_diff21 = List.fold_left (fun acc x -> remove_partition x acc) p2 p1 in
        let p_diff12 = List.fold_left (fun acc x -> remove_partition x acc) p1 p2 in
        (List.length p_diff21) + (List.length p_diff12) in
    (* assertions; ensure both trees have the same nodes *)
    let all = get_all tree1 in
    assert( All_sets.Integers.is_empty (All_sets.Integers.diff all (get_all tree2)) );
    let p1 = partitions_of_tree all tree1 
    and p2 = partitions_of_tree all tree2 in
    compare_partitions all p1 p2


(** [post_order_node_visit f id bt ad acc]
@param f function to applied to all the nodes in post-order.
@param id the node_id from where the traversal is started.
@param bt the tree whose nodes are visited.
@param acc the result of the function application on the edges.
@return the function applied to the subtree of id as imposed by the
        handle in the component, the results are returned as a list. *)
let post_order_node_visit f id bt accum =
    let visit_node pred nd acc = (f pred nd acc) in
    let break_or_continue pred nd ret c_with =
        let status, acc = ret in
        match status with
        | Continue -> (c_with pred nd acc)
        | Break    -> (Break, acc)
        | _ -> (failwith "Invalid status returned")
    in
    let rec process_children pred nd acc = match pred with
        | Some p ->
              let node = (get_node nd bt) in
              let (x, y) = (other_two_nbrs p node) in
              let ret = (traverse (Some nd) x acc) in
              (break_or_continue (Some nd) y ret traverse)
        | _ -> raise (Invalid_argument "Needs a predecessor")
    and traverse pred id accum =
        match (get_node id bt) with
        | Leaf(nd, nbr) ->
              begin match pred with
                  | None -> assert(is_handle id bt)
                  | Some p -> assert(p = nbr);
              end;
              (visit_node pred nd accum)
        | Interior(nd, nbr1, nbr2, nbr3) ->
              let ret = (Continue, accum) in
              let ret = (break_or_continue pred nd ret process_children) in
              (break_or_continue pred nd ret visit_node)
        | _ -> raise (Invalid_argument "ponv 1: Check node type")
    in
    match (get_node id bt) with
    | Leaf (nd, nbr) ->
        let ret = (Continue, accum) in
        let ret = (break_or_continue (Some nd) nbr ret traverse) in
        let _, acc = (break_or_continue (Some nbr) nd ret visit_node) in
        acc
    | Interior (nd, nbr1, nbr2, nbr3) ->
        begin match (is_handle nd bt) with
            | true ->
                let ret = (Continue, accum) in
                let ret = break_or_continue (Some nd) nbr1 ret traverse in
                let ret = break_or_continue (Some nbr1) nd ret process_children in
                let _, acc = break_or_continue (Some nbr1) nd ret visit_node in
                acc
            | false ->
                let pred = (get_parent nd bt) in
                let ret = (Continue, accum) in
                let ret = break_or_continue (Some pred) nd ret process_children in
                let _, acc = break_or_continue (Some pred) nd ret visit_node in
                acc
          end
    | Single (_) -> accum

(** [get_pre_order_edges hs tree]
    @param hs the start handle.
    @param tree the tree whose edges in pre-order are desired.
    @return the edges of the tree as visited by a pre_order traversal
            in list form. *)
let get_pre_order_edges hs tree =
    let get_edge e acc =  Continue, e :: acc in
    let rev_edges = (pre_order_edge_visit get_edge hs tree []) in
       List.rev rev_edges

(** [get_edges tree hs] returns the edges in a tree component in an arbitrary
    order *)
let get_edges ?(acc=[]) tree hs =
    let get_edge e acc =  Continue, e :: acc in
    pre_order_edge_visit get_edge hs tree acc

(** [get_edges_tree tree] returns all the edges in the tree in an arbitrary
    order *)
let get_edges_tree tree =
    List.fold_left
        (fun acc h -> get_edges ~acc tree h)
        []
        (handle_list tree)

(*****************************************************************************)

(** [print_edge edge]
    @return prints the edge as a tuple. *)
let print_edge (Edge(e1, e2)) =
    let e1_s = (string_of_int e1) and
        e2_s = (string_of_int e2) in
    print_string ("(" ^ e1_s ^ "," ^ e2_s ^ ")")


(** [print_tree tree]
    @param id id of the node that will be used as a starting point.
    @param tree tree being printed.
    @return () - this function prints the binary tree to the screen
            using a backward-in-order traversal. *)
let print_tree id tree =
    let pipe = "|"
    and empty_string = ""
    and suffix = " " in
    let update_prefix prefix =
        (prefix ^ suffix ^ pipe) in
    let already_visited node_id visited =
        (All_sets.Integers.mem node_id visited) in
    let mark_visited node_id visited =
        (All_sets.Integers.add node_id visited) in
    let visit_node node_id prefix =
        (prerr_endline (prefix ^ (string_of_int node_id))) in
    let rec visit_nbr nbr_id prefix visited =
        if (already_visited nbr_id visited)
            then visited
            else (print_tree_aux nbr_id prefix visited)
    and print_tree_aux node_id prefix visited =
        match (get_node node_id tree) with
        | Single(id) ->
            (visit_node id prefix);
            (mark_visited id visited)
        | Leaf(id, pid) ->
            let visited_n = (mark_visited node_id visited) in
            let new_prefix = (update_prefix prefix) in
            let visited_n_p = (visit_nbr pid new_prefix visited_n) in
                (visit_node id prefix);
                visited_n_p
        | Interior(id, pid, lid, rid) ->
            let visited_n = (mark_visited node_id visited) in
            let new_prefix = (update_prefix prefix) in
            let visited_p = (visit_nbr pid new_prefix visited_n) in
            let visited_n_l = (visit_nbr lid new_prefix visited_p) in
                (visit_node id prefix);
                (visit_nbr rid new_prefix visited_n_l)
    in
    ignore (print_tree_aux id empty_string All_sets.Integers.empty)

(** [print_forest forest]
    @return () - prints all the trees in the forest. *)
let print_forest forest =
    let handles = (All_sets.Integers.elements (get_handles forest)) in
    (ignore (List.map (fun x -> prerr_newline ();
                                print_tree x forest) handles))
module Parse = struct
    (* A simple representation of a tree *)
    type 'a t = Leafp of 'a | Nodep of 'a t list * 'a

    (* types of trees that can be returned *)
    type tree_types =
        | Flat of string t
        | Annotated of (tree_types * string)
        | Branches of (string * float option) t
        | Characters of (string * string option) t

    let remove_annotations = function
        | Annotated (t,_) -> t
        | t               -> t

    let rec print_tree (output:string option) (t:tree_types) : unit = 
        let c = Status.user_message (Status.Output (output, false, [])) in
        let fprint_out format = Printf.ksprintf c format in
        let rec interior f t : unit = match t with
            | []     -> ()
            | hd::[] -> f hd
            | hd::tl -> f hd; fprint_out ","; interior f tl
        in
        let rec print_ n_to_string l_to_string t : unit = match t with
            | Leafp d     -> l_to_string d
            | Nodep (n,d) ->
                fprint_out "(";
                interior (print_ n_to_string l_to_string) n;
                fprint_out ")";
                n_to_string d;
                ()
        in
        match t with
            | Annotated (t,_) -> print_tree output t
            | Flat t ->
                print_ (fun _ -> ()) (fun x -> fprint_out "%s" x) t
            | Branches t -> 
                print_ 
                    (fun (_,d) -> 
                        match d with | Some x -> fprint_out ":%g" x 
                                     | None   -> ())
                    (fun (n,d) -> 
                        match d with | Some x -> fprint_out "%s:%g" n x 
                                     | None   -> ())
                    t
            | Characters t ->
                print_ 
                    (fun (_,d) -> 
                        match d with | Some x -> fprint_out ":%s" x 
                                     | None   -> ())
                    (fun (n,d) -> 
                        match d with | Some x -> fprint_out "%s:%s" n x 
                                     | None   -> fprint_out "%s" n)
                    t

    exception Trailing_characters

    exception Illegal_tree_format of string

    (* map on tree w/out annotations *)
    let rec map fn = function
        | (Leafp d) -> Leafp (fn d)
        | Nodep (list,d) ->
              Nodep (List.map (map fn) list, fn d)
    (* map on post_processed tree *)
    let rec map_tree fn = function
        | Annotated (t,a) -> Annotated (map_tree fn t, a)
        | Branches t -> Branches (map (fun (x,d) -> (fn x,d)) t)
        | Flat t -> Flat (map fn t)
        | Characters t -> Characters (map (fun (x,d) -> (fn x,d)) t)

    (* fold_left on tree data *)
    let fold_left_data fn a t =
        let rec fold_left_data2 fn a t = match t with
            | Leafp d -> fn a d
            | Nodep (list,d) ->
                fn (List.fold_left (fold_left_data2 fn) a list) d
        in
        match t with
        | (Leafp d),str -> (fn a d)
        | ((Nodep (list, d)),str) ->
            fn (List.fold_left (fold_left_data2 fn) a list) d

    (* fold_left on tree annotations *)
    let fold_left_annot fn a t = match t with
        | (Leafp d),str -> fn a str
        | (Nodep (list, d)),str -> fn a str

    let gen_aux_of_stream_gen do_stream stream =
        let taxon_name x = not (FileStream.is_taxon_delimiter x) in
        let close_squared_parenthesis = [ ']' ] in
        (*
        let ignore_cost_bracket () =
            ignore (stream#read_excl close_squared_parenthesis);
            ignore (stream#getch);
            ()
        in
        *)
        let get_cost_bracket () =
            let res = stream#read_excl close_squared_parenthesis in
            ignore (stream#getch);
            res
        in
        let read_taxon_name () = 
            stream#read_while taxon_name
        in

        let rec read_branch_length acc : float option=
            match stream#getch with
            | v when not (taxon_name v) -> 
                stream#putback v;
                let str = 
                    acc --> List.rev
                        --> List.map (String.make 1)
                        --> String.concat ""
                in
                Some (float_of_string str)
            | ('0' .. '9') as v -> read_branch_length (v :: acc)
            | 'e' -> read_branch_length ('e' :: acc)
            | 'E' -> read_branch_length ('E' :: acc)
            | '.' -> read_branch_length ('.' :: acc)
            | '-' -> read_branch_length ('-' :: acc)
            |  v   -> raise (Illegal_tree_format ("Unexpected Char "^(Char.escaped v)))
        in

        let rec read_branch acc : (string * (float option * string option)) t =
            stream#skip_ws_nl;
            match stream#getch with
            | '(' -> 
                    let res = read_branch [] in
                    read_branch (res :: acc)
            | ')' -> 
                    Nodep (acc, ("",(None,None)))
            | '[' -> 
                    let bracket = Some (get_cost_bracket ()) in
                    let acc = (match acc with
                     | Leafp (t,(branch,_)) ::tl      -> Leafp ((t,(branch,bracket))) :: tl
                     | Nodep (a, (t,(branch,_))) ::tl -> Nodep ((a, (t,(branch,bracket)))) :: tl
                     | _ -> raise (Illegal_tree_format "Unexpected Character")
                    ) in
                    read_branch acc
            | ':' ->
                    let branch = read_branch_length [] in
                    let acc = (match acc with
                     | Leafp (t,(_,bracket)) ::tl      -> Leafp ((t,(branch,bracket))) :: tl
                     | Nodep (a,(t,(_,bracket))) ::tl -> Nodep ((a, (t,(branch,bracket)))) :: tl
                     | _ -> raise (Illegal_tree_format "Unexpected Character")
                    ) in
                    read_branch acc
            | v when taxon_name v ->
                    stream#putback v;
                    let taxon = read_taxon_name () in
                    read_branch ((Leafp (taxon,(None,None))) :: acc)
            | v ->
                    let character = stream#get_position in
                    let ch = Char.escaped v in
                    if List.exists (fun x -> x = ch) [","; ";"] then
                        read_branch acc
                    else
                        let () = 
                            let msg = 
                                "I@ will@ use@ the@ character@ " ^ ch ^ 
                                " in@ position@ " ^ string_of_int character ^ 
                                "@ as@ a@ taxon@ name@ separator" 
                            in
                            Status.user_message Status.Information msg;
                        in
                        read_branch acc
        in
        let rec read_tree acc1 acc2 =
            let prepend acc2 acc1 =
                (* We have to do this to avoid bogus warnings when we read a
                * separator that closes the list of trees (for example a
                * trailing semicolon) *)
                match acc2 with
                | [] -> acc1
                | acc2 -> acc2 :: acc1 
            in
            let () = try stream#skip_ws_nl; with | End_of_file -> () in
            match stream#getch_safe with
            | None -> prepend acc2 acc1
            | Some v -> (match v with
                | '(' -> 
                        let res = 
                            try read_branch [] with
                            | End_of_file -> 
                                    let msg = "Unexpected end of file" in
                                    raise (Illegal_tree_format msg)
                        in
                        read_tree acc1 ((res, "") :: acc2)
                | '*'
                | ';' -> 
                        let acc1 = prepend acc2 acc1 in
                        read_tree acc1 []
                | ':' ->
                        let dist = read_branch_length [] in
                        let acc2 = (match acc2 with
                         | ((Nodep (lst,(name,(None,bracket)))),q) :: tl ->
                                 ((Nodep (lst,(name,(dist,bracket)))),q) :: tl
                         | ((Leafp (name,(None,bracket))),q) :: tl ->
                                 ((Leafp (name,(dist,bracket))),q) :: tl
                         | _ -> raise 
                            (Illegal_tree_format "Unexpected Character")
                        ) in
                        read_tree acc1 acc2
                | '[' -> 
                        let contents = 
                            try get_cost_bracket () with
                            | End_of_file ->
                                    let msg = "Unexpected end of file" in
                                    raise (Illegal_tree_format msg)
                        in
                        let acc2 = 
                            match acc2 with
                            | (h, _) :: t -> (h, contents) :: t
                            | [] -> 
                                    let msg = "Unexpected cost spec" in
                                    raise (Illegal_tree_format msg)
                        in
                        read_tree acc1 acc2
                | v -> 
                        let character = stream#get_position in
                        let ch = Char.escaped v in
                        let message = "Unexpected character " ^ ch ^ 
                        " in position " ^ string_of_int character in
                        failwith message
                )
        in
        let read_tree_str =
            let acc2 = ref None in
            let rec tree_generator () =
                stream#skip_ws_nl;
                match stream#getch with
                | '(' -> 
                        let res = 
                            try read_branch [] with
                            | End_of_file -> 
                                    let msg = "Unexpected end of file" in
                                    raise (Illegal_tree_format msg)
                        in
                        (match !acc2 with
                        | None -> 
                                acc2 := Some (res, "");
                                tree_generator ()
                        | Some _ -> raise (Illegal_tree_format "Trees ignored!"))
                | '*'
                | ';' -> 
                        (match !acc2 with
                        | None -> raise (Illegal_tree_format "No trees to read?")
                        | Some tree ->
                                acc2 := None;
                                tree)
                | ':' -> let dist = read_branch_length [] in
                         (match !acc2 with
                          | Some ((Nodep (lst,(name,(_,bracket)))),q) ->
                                  acc2 := None;
                                 ((Nodep (lst,(name,(dist,bracket)))),q)
                          | Some ((Leafp (name,(_,bracket))),q) ->
                                  acc2 := None;
                                 ((Leafp (name,(dist,bracket))),q)
                          | _ -> raise 
                            (Illegal_tree_format "Unexpected Character")
                         )
                | '[' -> 
                        let contents = 
                            try get_cost_bracket () with
                            | End_of_file ->
                                    let msg = "Unexpected end of file" in
                                    raise (Illegal_tree_format msg)
                        in
                        (match !acc2 with
                        | Some (h, _) -> 
                                acc2 := None;
                                (h, contents) 
                        | None -> 
                                let msg = "Unexpected cost spec" in
                                raise (Illegal_tree_format msg))
                | v -> 
                        let character = stream#get_position in
                        let ch = Char.escaped v in
                        let message = "Unexpected character " ^ ch ^ 
                        " in position " ^ string_of_int character in
                        failwith message
            in
            tree_generator
        in
        if not do_stream then `Treess (read_tree [] [])
        else `Stream (read_tree_str)

    (* post_processing function for trees. Will figure out what type of tree we
    * are parsing and apply the correct varient type *)
    let post_process (t:(string * (float option * string option)) t * string):tree_types = 
        let branches_exist t = 
                fold_left_data
                   (fun a (s,(n,c)) -> match n with | None -> a | _ -> true)
                    false t
        and characters_exist t = 
                fold_left_data
                    (fun a (s,(n,c)) -> match c with | None -> a | _ -> true)
                    false t
        and annotations_exist (t,str) = match str with
            | "" -> false
            | _ -> true
        in
        let remove_annot_1branches (t,s) = map (fun (d,(n,c)) -> d,c) t
        and remove_annot_2branches (t,s) = map (fun (d,(n,c)) -> d,n) t
        and remove_all (t,s)  = map (fst) t in
        match annotations_exist t,branches_exist t,characters_exist t with
        | true , true, _     -> Annotated (Branches (remove_annot_2branches t),(snd t))
        | true ,  _  , true  -> Annotated (Characters (remove_annot_1branches t),(snd t))
        | true ,false, false -> Annotated (Flat (remove_all t),(snd t))
        | false, true, false -> Branches (remove_annot_2branches t)
        |   _  , true,  _    -> Branches (remove_annot_2branches t)
        |   _  ,  _  , true  -> Characters (remove_annot_1branches t)
        | false,false, false -> Flat (remove_all t)

    (** general function for trees *)
    let gen_aux_of_stream str = match gen_aux_of_stream_gen false str with
        | `Treess t -> List.rev_map (List.rev_map (post_process)) t
        | `Stream _ -> assert false

    (** [(gen_)of_string ...] **)
    let gen_of_string f str = 
        try
            let stream = new FileStream.string_reader str in
            f stream
        with
        | Trailing_characters -> 
                raise (Illegal_tree_format "Trailing characters in tree.")
    let of_string str = gen_of_string gen_aux_of_stream str

    (** [(gen_)of_channel ...] **)
    let gen_of_channel f ch =
        try
            let stream = FileStream.stream_reader ch in
            f stream
        with
        | End_of_file -> failwith "Unexpected end of file"
    let of_channel ch = gen_of_channel gen_aux_of_stream ch

    (** [(gen_)of_file ...] **)
    let gen_of_file f file =
        try
            let ch = FileStream.open_in file in
            let x = f ch in
            close_in ch;
            x
        with
        | Failure msg ->
                let file = FileStream.filename file in
                let msg = file ^ ": " ^ msg in
                raise (Illegal_tree_format msg)
        | (Sys_error err) as e ->
                let file = FileStream.filename file in
                let msg = "Couldn't@ open@ the@ trees@ file@ " ^ 
                StatusCommon.escape file ^ 
                ".@ The@ error@ message@ is@ @[" ^ err ^ "@]"in
                Status.user_message Status.Error msg;
                raise e
    let of_file file = gen_of_file of_channel file

    (** [stream of_file .. ] *)
    let stream_of_file is_compressed file =
        let ch, to_close = 
            if is_compressed then
                let protocol = 
                    let ch = FileStream.open_in_bin file in
                    let res = Lz.detect_type ch in
                    close_in ch;
                    res
                in
                match protocol with
                | `Zlib -> 
                        let real_ch = FileStream.open_in_gz file in
                        FileStream.gz_reader real_ch,
                        `Zlib real_ch
                | #Lz.protocol as protocol ->
                        match protocol with
                        | `Compressed x -> 
                                let real_ch = FileStream.open_in_bin file in
                                new FileStream.compressed_reader protocol real_ch, 
                                `Regular real_ch
                        | `NoCompression -> 
                                let real_ch = FileStream.open_in_bin file in
                                FileStream.stream_reader real_ch,
                                `Regular real_ch
            else
                let real_ch = FileStream.open_in_bin file in
                FileStream.stream_reader real_ch, `Regular real_ch
        in
        match gen_aux_of_stream_gen true ch with
        | `Treess _ -> assert false
        | `Stream s -> 
                (fun () -> post_process (s ())), 
                match to_close with
                | `Regular real_ch -> (fun () -> close_in real_ch)
                | `Zlib real_ch -> (fun () -> Gz.close_in real_ch)

    let rec cannonic_order tree =
        let rec build_cannonic_order fn = function
            | Leafp d -> (Leafp d, d)
            | Nodep (chld , cnt) ->
                    let nch = List.sort fn (List.map (build_cannonic_order fn) chld) in
                    let _, b = List.hd nch in
                    let nch = List.map (fun (x, _) -> x) nch in
                    Nodep (nch, cnt), b
        in
        match tree with
            | Annotated (t,str) ->
                Annotated (cannonic_order t, str)
            | Flat t ->
                let tree,_ = 
                    build_cannonic_order (fun (_,a) (_,b) -> compare a b) t
                in Flat tree
            | Branches t ->
                let tree,_ = build_cannonic_order 
                                (fun (_,(d1,_)) (_,(d2,_)) -> compare d1 d2)
                                t
                in Branches tree
            | Characters t ->
                let tree,_ = build_cannonic_order 
                                (fun (_,(d1,_)) (_,(d2,_)) -> compare d1 d2)
                                t
                in Characters tree

    exception Illegal_argument

    let rec strip_tree : tree_types -> string t = function
        | Annotated (t,_) -> strip_tree t
        | Flat t -> t
        | Characters t -> map (fst) t
        | Branches t -> map (fst) t

    let fold_left_tree_data (f: 'a -> string -> 'a) (acc:'a) (t:tree_types) : 'a =
        fold_left_data f acc ((strip_tree t),None)

    let rec maximize_tree = function
        | Characters t -> t
        | Branches t -> map (fun (x,_) -> x,None) t
        | Annotated (t,_) -> maximize_tree t
        | Flat t -> map (fun x -> x,None) t

    let rec cleanup ?newroot f tree =
        let rec aux_cleanup f tree = 
            match tree with
            | Leafp x -> 
                    if f x then []
                    else [tree]
            | Nodep (chld, x) ->
                    let res = List.flatten (List.map (aux_cleanup f) chld) in
                    match res with
                    | [] | [_] -> res
                    | y -> [Nodep (y, x)]
        in
        match tree with
        | Annotated (t,str) ->
            begin match cleanup ?newroot f t with
                | None -> None
                | Some t -> Some (Annotated (t,str))
            end
        | Branches t ->
            begin match aux_cleanup (fun (s,d) -> f s) t with
                | []  -> None
                | [x] -> Some (Branches x)
                | x   ->
                    begin match newroot with
                        | None -> raise Illegal_argument
                        | Some y -> Some (Branches (Nodep (x,(y,None))))
                    end
            end
        | Characters t -> 
            begin match aux_cleanup (fun (s,d) -> f s) t with
                | []  -> None
                | [x] -> Some (Characters x)
                | x   ->
                    begin match newroot with
                        | None -> raise Illegal_argument
                        | Some y -> Some (Characters (Nodep (x,(y,None))))
                    end
            end
        | Flat t -> 
            begin match aux_cleanup f t with
                | []  -> None
                | [x] -> Some (Flat x)
                | x   ->
                    begin match newroot with
                        | None -> raise Illegal_argument
                        | Some y -> Some (Flat (Nodep (x, y)))
                    end
            end

    let add_tree_to taxon_code add_to tree =
        (* below, fill in tree names for consistent tree *)
        let tree = maximize_tree tree in
        let rec assign_codes tree parent = function
            | Leafp (name,nname) ->
                let tc =
                    try taxon_code name
                    with | Not_found -> 
                        try int_of_string name
                        with | Not_found as err ->
                            Status.user_message Status.Error
                                ("Could@ not@ find@ data@ loaded@ for@ taxon@ "^
                                 StatusCommon.escape name ^ "@ in@ a@ loaded@ tree.");
                            raise err
                in
                tree, Leafp (Leaf (tc, parent)), tc
            | Nodep (child_nodes, (txt,nname)) ->
                let rec resolve_more_children (chil:(string * string option) t list) = match chil with
                    | [a; b] as x -> x
                    | [taxon] ->
                        Status.user_message Status.Error
                            ("Your@ tree@ file@ has@ a@ subtree@ or@ taxon@ " ^
                            "occurring@ as@ the@ unique@ member@ of@ an@ " ^
                            "internal@ node@ in@ the@ tree@ " ^
                            "(maybe@ there@ is@ a@ space@ missing@ between@ " ^
                            "taxa?).@ I@ am@ cancelling@ this@ read@ command");
                            failwith "Illegal Tree file format"
                    | [] -> 
                        Status.user_message Status.Error
                            ("Your@ tree@ file@ has@ a@ (),@ that@ is,@ "
                            ^"an@ opening@ parentheses@ followed@ immediately@ "
                            ^ "by@ a@ closing@ parentheses.@ Empty@ trees@ "
                            ^ "are@ not@ allowed@ in@ a@ tree@ file,@ so@ I@ "
                            ^"am@ cancelling@ this@ read@ command");
                        failwith "Illegal Tree file format"
                    | a :: b :: t ->
                        resolve_more_children ((Nodep ([a; b], (txt,nname))) :: t)
                in
                let child_nodes = resolve_more_children child_nodes in
                let sc,tree = get_available tree in
                match child_nodes with
                | [a; b] ->
                    let tree, ta, ca = assign_codes tree sc a in
                    let tree, tb, cb = assign_codes tree sc b in
                    tree, Nodep ([ta; tb], Interior (sc, parent, ca, cb)), sc
                | _ -> failwith "Tree.assign_codes"
        in
        let add_edge a b = EdgeSet.add (Edge (a, b)) in
        let remove_edge a b = EdgeSet.remove (Edge (a, b)) in
        let replace_parent vertices v par = 
            let vertex = All_sets.IntegerMap.find v vertices in
            let vertex = match vertex with
                | Interior (a, _, b, c) -> Interior (a, par, b, c)
                | Leaf (a, _) -> Leaf (a, par)
                | Single _ -> failwith "Unexpected Tree.replace_parent"
            in
            All_sets.IntegerMap.add v vertex vertices
        in
        let rec add_edges_n_vertices tree (edges, vertices) = match tree with
            | Nodep ([x; y], ((Interior (a, b, c, d)) as v)) ->
                let vertices = All_sets.IntegerMap.add a v vertices
                and edges =
                    edges --> add_edge b a
                          --> add_edge a c
                          --> add_edge a d
                in
                add_edges_n_vertices x (add_edges_n_vertices y (edges,vertices))
            | Leafp ((Leaf (a, b)) as v) ->
                let vertices = All_sets.IntegerMap.add a v vertices
                and edges = add_edge b a edges in
                (edges, vertices)
            | _ ->
                failwith "Unexpected Tree.add_edges_n_vertices"
        in
        let tree =  match tree with
            (* We clean up a tree that has only one leaf, this is valid when
               processing forest, but invalid inside a tree. *)
            | (Leafp _) as x
            | Nodep ([((Leafp _) as x)], _) -> x
            | y -> y
        in
        match assign_codes add_to (-1) tree with
        | _, Nodep ([a; b], (Interior (sc, _, ca, cb))), _ ->
                let edges, vertices =
                    add_edges_n_vertices b
                        (add_edges_n_vertices a (add_to.d_edges, add_to.u_topo))
                in
                let vertices = replace_parent vertices ca cb in
                let vertices = replace_parent vertices cb ca in
                let edges =
                    edges --> remove_edge sc ca
                          --> remove_edge sc cb
                          --> add_edge ca cb
                in
                let handles = All_sets.Integers.add ca add_to.handles in
                { add_to with
                    u_topo = vertices;
                    d_edges = edges;
                    handles = handles;
                }
        | _, Leafp (Leaf (tc, _)), _ ->
                let vertices =
                    All_sets.IntegerMap.add tc (Single tc) add_to.u_topo in
                let handles = All_sets.Integers.add tc add_to.handles in
                { add_to with
                        u_topo = vertices;
                        handles = handles;
                }
        | _ -> failwith "We need trees with more than two taxa"

    let convert_to ((name,trees): string option * tree_types list) taxon_code =
        let get_id_information taxon_code_fn tree : int list * int =
            let taxon_ids =
                List.fold_left
                    (fun acc t ->
                        fold_left_tree_data
                            (fun acc d ->
                                if d = "" then acc else (taxon_code_fn d) :: acc)
                            acc
                            t)
                    []
                    tree
            in
            let avail_id = ref [] in
            let maxim_id = List.fold_left (max) 0 taxon_ids in
            for i = 1 to maxim_id do
                if not (List.mem i taxon_ids) then
                    avail_id := i :: !avail_id
            done;
            !avail_id, (maxim_id+1)
        in
        let avail,next = get_id_information taxon_code trees in
        let tree =
            { (empty ()) with tree_name = name;
                              avail_ids = avail;
                                new_ids = next; }
        in
        List.fold_left (add_tree_to taxon_code) tree trees

end

module Fingerprint = struct
    type t = int Parse.t

    let find_smallest_leaf {u_topo = topo} =
        All_sets.IntegerMap.fold
            (fun k _ m -> min k m) topo max_int
    let move_handle_to_smallest t =
        let r, _ = move_handle (find_smallest_leaf t) t in
        r
    let fingerprint t =
        let smallest = find_smallest_leaf t in
        let t, _ = move_handle smallest t in
        let parent_id = 
            match get_node smallest t with
            | Leaf (_, p) -> p
            | _ -> failwith "Bad leaf node" 
        in
        let rec fp p = function
            | Single _ -> failwith "Strange place to leave a single"
            | Leaf (myid, _) -> myid, Parse.Leafp myid
            | (Interior (myid, c1, c2, c3)) as node ->
                  let c1, c2 = other_two_nbrs p node in
                  let c1id, c1l = fp myid (get_node c1 t) in
                  let c2id, c2l = fp myid (get_node c2 t) in
                  if c1id < c2id then c1id, Parse.Nodep ([c1l; c2l], c1id)
                  else c2id, Parse.Nodep ([c2l; c1l], c2id)
        in 
        let _, left = fp smallest (get_node parent_id t) in
        Parse.Nodep ([Parse.Leafp smallest; left], smallest)

    let compare = compare

end


module CladeFP = struct
    type fp = All_sets.Integers.t
    type t = fp EdgeMap.t
    let fpcompare = All_sets.Integers.compare

    (** [Ordered] allows for maps and sets of fingerprints *)
    module Ordered = struct
        type t = fp
        let compare = fpcompare
    end

    let calc tree =
        let module Set = All_sets.Integers in
        let module Map = EdgeMap in
        let calc_component h m =
            let comp_leaves = get_leaves tree h in
            let comp_leaves =
                List.fold_right Set.add comp_leaves Set.empty in
            (* complement *)
            let comp s = Set.diff comp_leaves s in
            let add p n set m =
                let m = Map.add (Edge (p, n)) set m in
                let m = Map.add (Edge (n, p)) (comp set) m in
                m in
            (* [r p n m] is parent, node, map *)
            let rec r p n m = 
                match get_node n tree with
                | Leaf (id, _) ->
                        let set = Set.singleton id in
                        set, add p n set m
                | Single _ -> assert (false)
                | Interior (_, _, a, b) ->
                        let aset, m = r n a m in
                        let bset, m = r n b m in
                        let set = Set.union aset bset in
                        set, add p n set m 
            in
            match get_node h tree with
            | Leaf (id, par) ->
                    let set = Set.singleton id in
                    let m = add par id set m in
                    let _, m = r id par m in
                    m
            | Single _ -> m
            | Interior (id, a, b, c) ->
                    let _, m = r id a m in
                    let _, m = r id b m in
                    let _, m = r id c m in
                    m 
        in
        (Set.fold calc_component tree.handles Map.empty : t)

    module CladeSet = Set.Make (Ordered)

    let sets t =
        let module Set = All_sets.Integers in
        let module Map = EdgeMap in
        Map.fold (fun _ set acc -> CladeSet.add set acc) (calc t) 
        CladeSet.empty

    let query = (EdgeMap.find : edge -> t -> fp)

    let num_leaves (fp : fp) = All_sets.Integers.cardinal fp

    let fold (fn : edge -> fp -> 'a -> 'a) t a =
        EdgeMap.fold fn t a

end

(** [get_break_handles delta tree] returns the ids of the handles of the left
    and right components of a tree after a break *)
let get_break_handles (a, b) tree =
    let h1 = match a with
    | `Single (i, _) -> i
    | `Edge (_, l1, l2, m) -> begin match m with
      | None -> handle_of l1 tree
      | Some m -> m
      end in
    let h2 = match b with
    | `Single (i, _) -> i
    | `Edge (_, l1, l2, m) -> begin match m with
      | None -> handle_of l1 tree
      | Some m -> m
      end in
    h1, h2

let fix_handle_neighbor h n tree = match get_node h tree with
    | Interior (_, par, a, b) ->
        if par = n then tree
        else if n = a then 
            let tree = remove_node h tree in
            add_node (Interior (h, a, par, b)) tree
        else if n = b then
            let tree = remove_node h tree in
            add_node (Interior (h, b, a, par)) tree
        else failwith "Tree.fix_handle_neighbor"
    | Leaf _ -> tree
    | Single _ -> failwith "Tree.fix_handle_neighbor 2"



(** {2 Tree Fusing} *)

(** [edge_summary tree] takes the edges in [tree] and returns a map from all
    nodes to their parents and a map from all nodes to a list of their children *)
let edge_summary tree =
    let parents  = ref All_sets.IntegerMap.empty in
    let children = ref All_sets.IntegerMap.empty in
    EdgeSet.iter
        (fun (Edge(f, t)) ->
            assert (false = All_sets.IntegerMap.mem t !parents);
            parents  := All_sets.IntegerMap.add t f !parents;
            children := match All_sets.IntegerMap.mem f !children with
                | false -> All_sets.IntegerMap.add f [t] !children
                | true  ->
                    All_sets.IntegerMap.add f
                        (t :: All_sets.IntegerMap.find f !children) !children)
        tree.d_edges;
    (!parents, !children)

(** [update_clade_with_edges node tree] returns [tree] with all nodes under
    [node] updated to respect the edges in [tree]
let update_clade_with_edges node tree =
    let (parents, children) = edge_summary tree in
    let rec r node tree' = match get_node node tree with
        | Leaf _ -> tree'                   (* leaves stay leaves *)
        | Interior (_, par, a, b) ->
            let par' = All_sets.IntegerMap.find node parents in
            let a', b' = match All_sets.IntegerMap.find node children with
                | [ a; b ] -> a, b
                | _ -> failwith "Tree.update_clade_with_edges too many/few children"
            in
            let tree' = 
                tree' --> remove_node node
                      --> add_node (Interior (node, par', a', b')) in
                      --> r a'
                      --> r b'
            tree'
        | Single _ -> tree'
    in
    r node tree *)

(** [destroy_under_node node tree] removes all the edges and nodes underneath
    [node] in [tree].  WARNING: this leaves the tree in an inconsistent state.
    Leaves are left as-is, and the node itself is untouched, both of which are
    inconsistent.  For some measure of usefulness, we set the node's parent to
    itself... *)
let rec destroy_under_node node tree =
    match get_node node tree with
    | Leaf _ -> tree                    (* leave leaves: we don't want their
                                           IDs to be recycled... *)
          (* TODO: set node's parent to itself *)
    | Single _ -> tree                  (* same with singles *)
    | Interior (_, par, a, b) ->
          (* remove edges *)
          let tree = tree --> remove_edge (Edge (node, a))
              --> remove_edge (Edge (node, b)) in
          (* recursively remove *)
          let tree = tree --> destroy_under_node a
              --> destroy_under_node b in
          (* remove self *)
          let tree = tree --> remove_node node in
          tree

(** [reconstruct_node (source, snode) (target, tparent)] reconstructs the
    topology in [source] tree into [target] tree.  The clade under node [snode]
    will be reconstructed in the target tree as the child of node [ *)
let rec reconstruct_node (source, snode) (target, tparent) =
    match get_node snode source with
    | Leaf (sid, _) ->
          let target = target --> remove_node sid
              --> add_edge (Edge (tparent, sid))
              --> add_node (Leaf (sid, tparent)) in
          if debug_fusing then print_endline ("leaf: " ^ string_of_int sid
                                              ^ " with parent "
                                              ^ string_of_int tparent);
          target, sid
    | Single _ -> failwith "Tree.reconstruct_node"
    | Interior (sid, spar, sa, sb) ->
          let tid, target = get_available target in
          if debug_fusing then print_endline ("interior: " ^ string_of_int tid
                                              ^ " with parent "
                                              ^ string_of_int tparent);
          let target, t1 =
              reconstruct_node (source, sa) (target, tid) in
          let target, t2 =
              reconstruct_node (source, sb) (target, tid) in
          let target = target --> add_node (Interior (tid, tparent, t1, t2))
              --> add_edge (Edge (tparent, tid)) in
          target, tid

(** [source_to_target (source, snode) (target, tnode)] fuses the clade in
    [source] rooted at [snode] into [target], beginning at node [tnode] *)
let source_to_target (source, snode) (target, tnode) =
    assert (not (is_handle tnode target));
    assert (not (is_leaf tnode target));

    (* make sure the set of leaves is the same *)
    if debug_fusing then begin
        let sleaves = get_leaves source snode in
        let tleaves = get_leaves target tnode in
        let sleaves = List.sort compare sleaves in
        let tleaves = List.sort compare tleaves in
        print_string "sleaves: ";
        print_endline (String.concat "; " (List.map string_of_int sleaves));
        print_string "tleaves: ";
        print_endline (String.concat "; " (List.map string_of_int tleaves));
        if sleaves <> tleaves
        then begin
            print_endline "sets of leaves disagree!\n";
            assert false
        end
    end;
    let tparent = get_parent tnode target in
    let target = target --> destroy_under_node tnode
        --> remove_edge (Edge (tparent, tnode)) in
    let (target, tnewid) =
        reconstruct_node (source, snode) (target, tparent) in
    assert (tnewid = tnode);            (* TODO: if the node has changed,
                                           update the parent to refer to the
                                           new node (...) *)
    if debug_fusing then test_tree target;
    target


module CladeSet = Set.Make (CladeFP.Ordered)
module CladeFPMap = Map.Make (CladeFP.Ordered)

let tree_fps_map ?(map=CladeFPMap.empty) tree =
    let fps = CladeFP.calc tree in
    let store edge fp map =
        try
            let list = CladeFPMap.find fp map in
            let map = CladeFPMap.remove fp map in
            CladeFPMap.add fp ((tree, edge) :: list) map
        with Not_found ->
            CladeFPMap.add fp [(tree, edge)] map in
    CladeFP.fold store fps map

let tree_fps_map_withaux ?(map=CladeFPMap.empty) (aux, tree) =
    let fps = CladeFP.calc tree in
    let store edge fp map =
        try
            let list = CladeFPMap.find fp map in
            let map = CladeFPMap.remove fp map in
            CladeFPMap.add fp ((aux, tree, edge) :: list) map
        with Not_found ->
            CladeFPMap.add fp [(aux, tree, edge)] map in
    CladeFP.fold store fps map

type 'a fuse_locations = ('a * u_tree * edge) list Sexpr.t

let fuse_locations ?(filter=fun _ -> true) sources (taux, target) =
    let target_fps = tree_fps_map target in
    let source_fps =
        let folder map tree = tree_fps_map_withaux ~map tree in
        List.fold_left folder CladeFPMap.empty sources in
    let locations =
        CladeFPMap.fold
            (fun key value sexpr ->
                try let _, tedge = match value with
                        | [e] -> e
                        | _ -> failwith "Tree.fuse_locations.locations"
                    in
                    let value = CladeFPMap.find key source_fps in
                    if filter (key, CladeFP.num_leaves key)
                        then (`Single ((taux, target, tedge) :: value)) :: sexpr
                        else sexpr
                with Not_found -> sexpr)
            target_fps []
    in
    match locations with
    | [] -> `Empty
    | l -> `Set l

let exchange_codes a b tree =
    let add_edge ?(force=false) ((Edge (a, b)) as e) tree =
        if is_handle b tree && not force then add_edge (Edge (b, a)) tree
        else add_edge e tree
    in
    assert (All_sets.IntegerMap.fold (fun key node acc ->
        match node with
        | Single _ -> acc
        | Leaf (a, b) -> (a <> b) && acc
        | Interior (a, b, c, d) ->
                (a <> b) && (a <> c) && (a <> d) && acc) tree.u_topo true);
    let replace_edge org_tree a prev next tree =
        let ((Edge (a, b)) as e) = normalize_edge (Edge (a, prev)) org_tree in
        let tree = remove_edge e tree in
        if a = prev then 
            add_edge (Edge (next, b)) tree
        else 
            add_edge (Edge (a, next)) tree
    in
    let map_quadruple (u, v, w, x) a b =
        let exchange x = if x = b then a else if x = a then b else x in
        let u = exchange u and v = exchange v 
        and w = exchange w and x = exchange x in
        (u, v, w, x)
    in
    let replace_neighbor org_tree node prev next tree =
        assert (prev <> node);
        match get_node node org_tree with
        | Leaf (x, y) ->
            assert (x = node);
            assert (y = prev);
            tree
                --> remove_node node
                --> replace_edge org_tree node prev next
                --> add_node (Leaf (x, next))
        | Single _ -> assert false
        | Interior (u, v, w, x) ->
            tree
                --> remove_node node
                --> replace_edge org_tree node prev next
                --> add_node (Interior (map_quadruple (u, v, w, x) prev next))
    in
    let replace_handles a b tree =
        if is_handle a tree || is_handle b tree then
            if not (is_handle a tree && is_handle b tree) then
                let proc_handle x y = tree --> remove_handle x --> add_handle y in
                if is_handle a tree
                    then proc_handle a b
                    else proc_handle b a
            else tree
        else tree
    in
    if a = b then tree
    else
        let pick_to_do f ow c t = if c then f t else ow t in
        let verify_shared_neighs a b c x =
            Printf.fprintf stderr "Verifying %d %d %d %d\n%!" a b c x;
            x <> a && x <> b && x <> c
        in
        let an = get_node a tree in
        if All_sets.IntegerMap.mem b tree.u_topo then
            let bn = get_node b tree in
            (* Both are valid codes and we must exchange their respective codes *)
            match an, bn with
            | Single _, _
            | _, Single _
            | Leaf _, _
            | _, Leaf _ -> failwith "Illegal argument"
            | Interior ((_, n, o, p) as h), Interior ((_, r, s, t) as i) ->
                (* We have to first check if the two vertices are neighbors, if yes
                * we better do something different or things will be just incorrect *)
                assert (a <> b);
                assert (a <> n);
                assert (a <> o);
                assert (a <> p);
                assert (b <> r);
                assert (b <> s);
                assert (b <> t);
                assert (verify_edge (Edge (b, r)) tree);  
                assert (verify_edge (Edge (b, s)) tree);  
                assert (verify_edge (Edge (b, t)) tree);  
                assert (verify_edge (Edge (a, n)) tree);  
                assert (verify_edge (Edge (a, o)) tree);  
                assert (verify_edge (Edge (a, p)) tree);  
                try
                    let (Edge (x, y)) as e = normalize_edge (Edge (a, b)) tree in
                    let a_n1, a_n2 = other_two_nbrs b an in
                    let b_n1, b_n2 = other_two_nbrs a bn in
                    tree
                        --> remove_edge e
                        --> remove_node a
                        --> remove_node b
                        --> replace_neighbor tree a_n1 a b
                        --> replace_neighbor tree a_n2 a b
                        --> replace_neighbor tree b_n1 b a
                        --> replace_neighbor tree b_n2 b a
                        --> add_node (Interior (map_quadruple h a b))
                        --> add_node (Interior (map_quadruple i a b))
                        --> add_edge ~force:true (Edge (y, x))
                        --> replace_handles a b
                with Invalid_Edge ->
                    (* OK they are not neighbors, let's keep going with the standard
                    * replacement *)
                    tree
                        --> remove_node b
                        --> remove_node a
                        --> remove_edge (Edge (n, a))
                        --> remove_edge (Edge (r, b))
                        --> add_node (Interior (map_quadruple h a b))
                        --> add_node (Interior (map_quadruple i a b))
                        --> replace_neighbor tree n a b
                        --> replace_neighbor tree o a b
                        --> replace_neighbor tree p a b
                        --> pick_to_do (replace_neighbor tree r b a)
                                       (fun x -> x)
                                       (verify_shared_neighs n o p r)
                        --> pick_to_do (replace_neighbor tree s b a)
                                       (add_edge (Edge (a, s)))
                                       (verify_shared_neighs n o p s)
                        --> pick_to_do (replace_neighbor tree t b a)
                                       (add_edge (Edge (a, t)))
                                       (verify_shared_neighs n o p t)
                        --> add_edge (Edge (r, a))
                        --> add_edge (Edge (n, b))
                        --> replace_handles a b
        else
            match an with
            | Single _
            | Leaf _ -> failwith "Illegal argument"
            | Interior ((_, x, y, z) as h) ->
                tree
                    --> remove_node a
                    --> remove_edge (Edge (x, a))
                    --> remove_edge (Edge (a, x))
                    --> add_node (Interior (map_quadruple h a b))
                    --> replace_neighbor tree x a b
                    --> replace_neighbor tree y a b
                    --> replace_neighbor tree z a b
                    --> replace_handles a b

let exchange_codes a b tree = 
    test_tree tree;
    let res = exchange_codes a b tree in
    test_tree res;
    res

let fuse_all_locations ?(filter=fun _ -> true) trees =
    let fps =
        let folder map tree = tree_fps_map_withaux ~map tree in
        List.fold_left folder CladeFPMap.empty trees 
    in
    (* Now we convert them in a list and we sort them by size, so that we can
    * decide if the clade that we are observing occurs in both and has a
    * different resolution efficiently *)
    let fps = 
        let sort_function (a, _) (b, _) = 
            (CladeFP.num_leaves a) - (CladeFP.num_leaves b)
        in
        [] --> CladeFPMap.fold (fun key value acc -> (key, value) :: acc) fps
        --> List.sort sort_function 
    in
    let atleast2 = function
        | a :: b :: _ -> true
        | _ -> false 
    in
    let _, locations = List.fold_left
        (fun (potential, sexpr) (key, value) ->
            let my_filter = filter (key, CladeFP.num_leaves key) in
            if atleast2 value && my_filter then 
                let has_potential, potential =
                    CladeSet.fold 
                    (fun pot (acc, new_pot) -> 
                        assert (CladeFP.num_leaves pot <= CladeFP.num_leaves
                        key);
                        let has_it = All_sets.Integers.subset pot key in
                        if has_it then true, CladeSet.remove pot new_pot
                        else acc, new_pot) potential
                        (false, potential) 
                in
                let sexpr = 
                    if has_potential then 
                        (`Single value) :: sexpr
                    else sexpr 
                in
                potential, sexpr
             else if not (atleast2 value) && my_filter then 
                 (CladeSet.add key potential), sexpr
             else potential, sexpr) (CladeSet.empty, []) fps  in
    match locations with
    | [] -> `Empty
    | l -> `Set l

let fuse_cladesize ~min ~max (_, size) = min <= size && size <= max

let fuse ~source ~target =
    let (stree, sedge) = source in
    let (ttree, tedge) = target in
    let Edge(sf, st) = sedge in
    let Edge(tf, tt) = tedge in
    let stree, spath =
        if is_edge sedge stree
        then stree, []
        else let Edge(f, t) = sedge in move_handle f stree in
    let ttree, tpath =
        if is_edge tedge ttree
        then ttree, []
        else let Edge(f, t) = tedge in move_handle f ttree in
    source_to_target (stree, st) (ttree, tt)


let rec update_vertex tree depth vertex acc =
    let cur_depth = All_sets.TupleMap.find (vertex, -1) acc in
    if cur_depth > depth + 1 then
        match get_node vertex tree with
        | Interior (_, par, c1, c2) ->
                let my_depth = depth + 1 in
                let acc = All_sets.TupleMap.add (vertex, par) my_depth acc in
                let acc = update_vertex tree my_depth c1 acc in
                update_vertex tree my_depth c2 acc
        | _ -> acc
    else acc

(** Return the best depth from a node to the vertex given; a general case to the
    ones based on a handle **)
let rec depth_of_vertex tree best_known vertex acc =
    match get_node vertex tree with
    | Interior (_, par, c1, c2) ->
            let dc1, acc = depth_of_vertex tree (best_known + 1) c1 acc in
            let dc2, acc = 
                let new_best_known = min (best_known + 1) (dc1 + 1) in
                depth_of_vertex tree new_best_known c2 acc 
            in
            let my_depth, acc = 
                if dc2 < dc1 && dc2 < best_known then
                    let my_depth = dc2 + 1 in
                    my_depth, update_vertex tree my_depth c1 acc
                else if dc1 < best_known then
                    dc1 + 1, acc
                else best_known + 1, acc
            in
            my_depth, All_sets.TupleMap.add (vertex, par) my_depth acc
    | Leaf (_, par) ->
            1, All_sets.TupleMap.add (vertex, par) 1 acc
    | Single _ -> 0, acc

(** Find the maximum depth of a tree with the handle mentioned *)
let depth_of_handle tree handle acc =
    let big_number = 10000000 in
    match get_node handle tree with
    | Interior (_, par, c1, c2) ->
            let dp, acc = depth_of_vertex tree big_number par acc in
            let dc1, acc = depth_of_vertex tree dp c1 acc in
            let dc2, acc = depth_of_vertex tree (min dp dc1) c2 acc in
            let mydepth, acc = 
                if dc2 < dp && dc2 < dc1 then 
                    let acc = update_vertex tree dc2 par acc in
                    dc2 + 1, update_vertex tree dc2 c1 acc
                else if dc1 < dp then
                    dc1 + 1, update_vertex tree dc1 par acc
                else dp, acc
            in
            All_sets.TupleMap.add (handle,par) mydepth acc
    | Leaf (_, par) ->
            let dp, acc = depth_of_vertex tree 1 par acc in
            All_sets.TupleMap.add (handle, par) dp acc
    | Single _ -> acc

(** Return a map from edges of handles to the depth of the tree.*)
let depths tree =
    All_sets.Integers.fold
        (depth_of_handle tree) (get_handles tree) All_sets.TupleMap.empty

(* Functions to compare two trees weather or not they are the same *)
let reroot (a, b) tree =
    let bt, _ = move_handle a tree in
    fix_handle_neighbor a b bt

(* Find a leaf *)
let choose_leaf tree =
    let (Edge (a, _)) = EdgeSet.choose tree.d_edges in
    let rec find_leaf x = match get_node x tree with
        | Leaf (a, b)           ->  a, b
        | Interior (_, _, c, _) -> find_leaf c
        | Single _              -> failwith "Tree.choose_leaf"
    in
    find_leaf a

let cannonize_on_edge ((a, b) as edge) tree = 
    assert (is_leaf a tree);
    let tree = reroot edge tree in
    let rec my_cannonizer vertex res = 
        match get_node vertex tree with
        | (Interior (a, b, c, d)) as vx ->
                (* Recursively cannonize the children and then myself *)
                let res, minc = my_cannonizer c res in
                let res, mind = my_cannonizer d res in
                if minc < mind then
                    All_sets.IntegerMap.add a vx res, minc
                else 
                    let v = (Interior (a, b, d, c)) in
                    All_sets.IntegerMap.add a v res, mind
        | vx ->
                (* There is nothing to cannonize in leafs and singles *)
                All_sets.IntegerMap.add vertex vx res, vertex
    in
    let res, _ = my_cannonizer b All_sets.IntegerMap.empty in
    let res = All_sets.IntegerMap.add a (Leaf (a, b)) res in
    { tree with u_topo = res }

let cannonize_on_leaf a tree =
    try match get_node a tree with
        | Leaf (a, b) -> cannonize_on_edge (a, b) tree
        | _ -> failwith "Tree.cannonize_on_leaf: the vertex is not a leaf"
    with | Not_found as err ->
        Status.user_message Status.Error
            ("Could not find " ^ string_of_int a);
        raise err

let destroy_component handle tree =
    assert (is_handle handle tree);
    let remove_node a tree =
        if All_sets.IntegerMap.mem a tree.u_topo then
            if (is_leaf a tree) then 
                tree --> remove_node a --> add_node (Single a) --> add_handle a
            else 
                if is_handle a tree then
                    tree --> remove_handle a --> remove_node a 
                else remove_node a tree
        else tree
    in
    let remove_edge e tree =
        if EdgeSet.mem e tree.d_edges then remove_edge e tree
        else tree
    in
    let edges = get_edges tree handle in
    let tree =
        List.fold_left (fun tree ((Edge (a, b)) as e) ->
            tree --> remove_node a --> remove_node b --> remove_edge e)
            tree edges
    in
    test_tree tree;
    tree

let copy_component handle source target =
    let edges = get_edges source handle in
    let add_node a tree =
        let tree =
            if All_sets.IntegerMap.mem a tree.u_topo
                then tree --> remove_node a
                else tree
        in
        add_node (get_node a source) tree
    in
    let target =
        List.fold_left
            (fun tree ((Edge (a, b)) as e) ->
                tree --> add_node a --> add_node b --> add_edge e)
            target edges
    in
    let avail, handles =
        List.fold_left
            (fun (lst, hdls) (Edge (a, b)) ->
                List.filter (fun x -> x <> a && x <> b) lst,
                (hdls --> All_sets.Integers.remove a --> All_sets.Integers.remove b))
            (target.avail_ids, target.handles)
            edges
    in
    let handles = All_sets.Integers.add handle handles in
    let res = { target with avail_ids = avail; handles = handles } in
    test_tree res;
    res

