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

(** rtree.ml - implements a rooted forest. *)
let () = SadmanOutput.register "Rtree" "$Revision: 1644 $"


(** Exn raised when node id is invalid. *)
exception Invalid_Node_Id

(** The node type of the rooted tree. The conventions are described in the
* comments below. *)
type node = Root of int * int * int 
            (** Root(id, left_id, right_id) *)
            | Interior of int * int * int * int 
            (** Interior(id, parent_id, left_id, right_id) *)
            | Leaf of int * int 
            (** Leaf(id, parent_id) *)

type r_tree = { r_topo : node All_sets.IntegerMap.t ;
                (** The topology of the rooted forest. *)
                roots : All_sets.Integers.t }
                (** Set of roots of the rooted forest. *)

(** The empty rooted tree. *)
let empty = { r_topo = All_sets.IntegerMap.empty ;
              roots = All_sets.Integers.empty }

(*************************** GET FUNCTIONS ***************************)

(** [get_id node] 
    @param node the node whose id is desired
    @return id of the node. *)
let get_id node =
    match node with
    | Root(id, _, _) -> id
    | Leaf(id, _) -> id 
    | Interior(id, _, _, _) -> id 

(** [get_node id rtree]
    @param id id of the node desired.
    @param tree tree to which the node belongs
    @return returns the node associated with the id.
    @raise Invalid_Node_Id, if the id is not a valid id. *)
let get_node id tree = 
    try 
        (All_sets.IntegerMap.find id tree.r_topo)
    with 
        | Not_found -> raise Invalid_Node_Id

(** [get_roots tree]
    @return the handles of the tree. *)
let get_roots tree = 
    tree.roots     

        
(** [is_root id rtree]
    @param id the id being checked whether it belongs to the root or not. 
    @param rtree the rooted tree in which id is either a root or not. *)
let is_root id rtree = 
    (* get the node corresponding to the id. *)
    match (get_node id rtree) with
    (* if it is a root node, make sure it exists in the root set and 
    * return true. *)
    | Root(_) -> 
        assert(All_sets.Integers.mem id rtree.roots) ;
        true 
    (* otherwise, return false. *)
    | _ -> false 

(** [is_leaf id rtree] 
    @param id the id of the node being checked to determine whether it
            is a leaf or not.
    @param rtree the rooted tree in which the id belongs to a leaf node or
            not. *)
let is_leaf id rtree = 
    (* get the node. *)
    match (get_node id rtree) with
    (* if the node is a leaf, then return true. *)
    | Leaf(_) -> true 
    (* otherwise false. *)
    | _ -> false 

(********************* BEGIN - UNSAFE OPS ****************************)

(** [add_node node rtree]
    @param node the node that needs to be added to the tree.
    @param tree the tree to which the node is being added.
    @return the tree with the node added.
    This function does not check whether the node is already
    present in the tree or not. If the node is already present then
    the tree is returned unmodified. *) 
let add_node node rtree = 
    let id = (get_id node) in
    let new_topo = (All_sets.IntegerMap.add id node rtree.r_topo) in
    let rtree = { rtree with r_topo = new_topo } in
    (* check whether the added node was a root. *)
    match node with
    | Root(_) -> 
        (* if it is, add it to the root set. *)
        { rtree with roots = All_sets.Integers.add id rtree.roots }
        (* else, do nothing. *)
    | _ -> rtree 

(** [add_nodes nodes rtree]
    @param nodes list of nodes that need to be added to the tree.
    @param btree tree to which the nodes are being added.
    @return tree with the list of nodes added. Older nodes with the
    same ids will be overwritten. *)
let add_nodes nodes rtree = 
    let f rtree node = (add_node node rtree) in
        (List.fold_left f rtree nodes)

(********************** END - UNSAFE OPS *****************************)

(** [pre_order_node_visitor visitor_op node_id rtree accum]
    @param visit_op function to applied to each node as it is visited.
    @param node_id the node whose children are visited in pre-order. When
        node_id is root, then all the nodes of the tree are visited; 
        otherwise only the children of the node are visited.
    @param rtree the tree whose nodes are being visited.
    @param accum the accumulator that accumulates the results of all the 
                 visits. 
    @return the accumulator with the results of the node visit. *)
let pre_order_node_visitor visit_op node_id rtree accum = 
    (* function to visit the node. *)
    let visit node_id accum = 
        (visit_op node_id accum) 
    in
    (* function to visit the node and process the children. *)
    let rec visit_node_n_children node_id left_id right_id accum =
        (* visit the node. *)
        let accum = (visit node_id accum) in
        (* visit the left sub-tree. *)
        let accum = (traverse left_id accum) in
        (* visit the right sub-tree. *)
        let accum = (traverse right_id accum) in
            (* return the accumulator. *)
            accum 
    (* function to traverse the tree in pre-order. *)
    and traverse node_id accum = 
        (* get the node corresponing to the start node id. *)
        match (get_node node_id rtree) with
        (* if it is the root node. *)
        | Root(root_id, left_id, right_id) -> 
            (* make sure the root is included in the roots set. *)
            assert(is_root root_id rtree) ;
            (* visit the node and its children. *)
            (visit_node_n_children root_id left_id right_id accum)
        (* if it is an interior node. *)
        | Interior(int_id, _, left_id, right_id) ->
            (* visit the node and its children. *)
            (visit_node_n_children int_id left_id right_id accum)
        (* if it is a leaf node *)
        | Leaf(leaf_id, _) -> 
            (* just, visit the node. *)
            (visit leaf_id accum) 
    in
        (* traverse the tree starting at the node_id. *)
        (traverse node_id accum) 

(** [post_order_node_visitor visit_op node_id rtree accum]
    @param visit_op function to applied to each node as it is visited.
    @param node_id the node whose children are visited in pre-order. When
        node_id is root, then all the nodes of the tree are visited; 
        otherwise only the children of the node are visited.
    @param rtree the tree whose nodes are being visited.
    @param accum the accumulator that accumulates the results of all the 
                 visits. 
    @return the accumulator with the results of the node visit. *)
let post_order_node_visitor visit_op node_id rtree accum = 
    (* function to visit the node. *)
    let visit node_id accum = 
        (visit_op node_id accum) 
    in
    (* function to visit the node and process the children. *)
    let rec visit_node_n_children node_id left_id right_id accum =
        (* visit the left sub-tree. *)
        let accum = (traverse left_id accum) in
        (* visit the right sub-tree. *)
        let accum = (traverse right_id accum) in
        (* visit the node. *)
        let accum = (visit node_id accum) in
            (* return the accumulator. *)
            accum 
    (* function to traverse the tree in pre-order. *)
    and traverse node_id accum = 
        (* get the node corresponing to the start node id. *)
        match (get_node node_id rtree) with
        (* if it is the root node. *)
        | Root(root_id, left_id, right_id) -> 
            (* make sure the root is included in the roots set. *)
            assert(is_root root_id rtree) ;
            (* visit the node and its children. *)
            (visit_node_n_children root_id left_id right_id accum)
        (* if it is an interior node. *)
        | Interior(int_id, _, left_id, right_id) ->
            (* visit the node and its children. *)
            (visit_node_n_children int_id left_id right_id accum)
        (* if it is a leaf node *)
        | Leaf(leaf_id, _) -> 
            (* just, visit the node. *)
            (visit leaf_id accum) 
    in
        (* traverse the tree starting at the node_id. *)
        (traverse node_id accum) 

(** [get_num_nodes rtree]
    @return the number of nodes in the subtree rooted at root_id. *)
let get_num_nodes root_id rtree = 
    let visit_op node_id num_nodes = 
        (num_nodes + 1) 
    in
        (pre_order_node_visitor visit_op root_id rtree 0) 
        
(** [index_nodes start_index root_id rtree]
    @param start_index the index at which the indexing starts.
    @param root_id identifies the tree 
    @param rtree the rooted forest.
    @return a pair of hash-tables, 1. indx to node 2. node to index.
*) 
let index_nodes (li, ii) num_nodes root_id rtree = 
    (* Function to add the current node_id to the index table. *)
    let visit_op node_id ((cli, cii), indx_nd_tbl, nd_indx_tbl) =
        (* if the node is a leaf *)
        if (is_leaf node_id rtree) then begin
            (* add it to the table. *)
            (Hashtbl.add nd_indx_tbl node_id cli) ;
            (Hashtbl.add indx_nd_tbl cli node_id) ;
            (* increment leaf index. *)
            ((cli + 1, cii), indx_nd_tbl, nd_indx_tbl)
        (* else, if it is an interior node or root node. *)
        end else begin
            (* add it to the table. *)
            (Hashtbl.add nd_indx_tbl node_id cii) ;
            (Hashtbl.add indx_nd_tbl cii node_id) ;
            (* increment the interior index. *)
            ((cli, cii + 1), indx_nd_tbl, nd_indx_tbl)
        end
    in
    let nd_indx_tbl = (Hashtbl.create num_nodes) in
    let indx_nd_tbl = (Hashtbl.create num_nodes) in
    let (_, indx_nd_tbl, nd_indx_tbl) = 
        (pre_order_node_visitor visit_op 
                                root_id 
                                rtree 
                                ((li, ii), indx_nd_tbl, nd_indx_tbl)) in
        (indx_nd_tbl, nd_indx_tbl)

(** [of_tree rt_id rt_jxn utree]
    @param rt_id the id of the new root.
    @param rt_jxn the jxn where the new root would be added.
    @param utree the unrooted tree that would be rooted at the root_jxn.
    @return a rooted version of the unrooted tree. *)
let of_tree rt_id rt_jxn utree = 
    (* Get the end-points of the root_jxn *)
    assert (not (All_sets.IntegerMap.mem rt_id utree.Tree.u_topo));
    let (hd_id, nd_id) = rt_jxn in
    let hd = Tree.int_of_id hd_id in
    let nd = Tree.int_of_id nd_id in
    (* Create an rtree with the root node. *)
    let rtree = (add_node (Root(rt_id, hd, nd)) empty) in
    (* function to add a utree node into the rtree. *)
    let add_to_rtree par_id nd_id (utree, rtree) = 
        (* get the parent of the node. *)
        let parent = 
            match par_id with
            | Some p -> p
            | None -> rt_id in
        (* Special case when the parent happens to be the root. *)
        let par_rt = 
            if ((parent = hd) || (parent = nd)) then
                rt_id
            else
                parent
        in
        (* get the node. *)
        let node = (Tree.get_node nd_id utree) in
            match node with
            | Tree.Leaf(_, p_id) ->
                (* if the node is a leaf, make sure the parents are
                * consistent. *)
                (Tree.Continue, 
                    (utree, (add_node (Leaf(nd_id, par_rt)) rtree)))
            | Tree.Interior(_, nb1, nb2, nb3) ->
                (* if the node is an interior node. get the other two nbrs*)
                let (x, y) = Tree.other_two_nbrs parent node in
                    (Tree.Continue, 
                     (utree, (add_node (Interior(nd_id, par_rt, x, y)) rtree)))
            | Tree.Single(_) ->
                raise (Failure "Single node found.") 
    in
    let _, rtree = 
        Tree.pre_order_node_visit add_to_rtree nd utree (utree, rtree) in
    let utree, _ = Tree.move_handle nd utree in
    let _, rtree = 
        Tree.pre_order_node_visit add_to_rtree hd utree (utree, rtree) in
        rtree

(** [root_at out_grp rt_id utree]
    @param out_grp the out group.
    @param rt_id the id of the root node.
    @param utree the unrooted tree being rooted at the out-group. 
    @return a rooted version of the unrooted tree, rooted at the given out
    group. *)
let root_at out_grp rt_id utree = 
    (* get the nbr of the out group. *)
    let nbr = match Tree.get_node out_grp utree with
                | Tree.Leaf(_, nbr) -> nbr 
                | _ -> raise (Failure "Outgroup must be a leaf") in
    let utree, _ = (Tree.move_handle out_grp utree) in
    let hnd_id = Tree.get_handle_id out_grp utree in
    let nbr_id = Tree.get_node_id nbr utree in
    let brk_jxn = (hnd_id, nbr_id) in
        (of_tree rt_id brk_jxn utree)
        
(** [distance num_nodes (rt_id1, rt1) (rt_id2, rt2)] 
    @param num_nodes number of nodes in each of the trees. 
    @param (rt_id, rt) root, rooted_tree pairs.
    @return the robinson-foulds distance between the trees on the same set 
    of leaves. *)
let distance num_nodes (rt_id1, rt1) (rt_id2, rt2) = ()
    
        
(** [print_tree tree]
    @param id id of the node that will be used as a starting point.
    @param tree tree being printed.
    @return () - this function prints the binary tree to the screen
            using a backward-in-order traversal. *)
let print_tree id tree = 
    let pipe = "|" and 
        empty_string = "" and
        suffix = " " in
    let update_prefix prefix = 
        (prefix ^ suffix ^ pipe) in
    let already_visited node_id visited = 
        (All_sets.Integers.mem node_id visited) in
    let mark_visited node_id visited =
        (All_sets.Integers.add node_id visited) in
    let visit_node node_id prefix = 
        (print_endline (prefix ^ (string_of_int node_id))) in
    let rec visit_nbr nbr_id prefix visited = 
        if (already_visited nbr_id visited) then
            visited
        else
            (print_tree_aux nbr_id prefix visited)
    and print_tree_aux node_id prefix visited = 
        match (get_node node_id tree) with
        | Leaf(id, pid) -> 
            let visited_n = (mark_visited node_id visited) in
            let new_prefix = (update_prefix prefix) in
            let visited_n_p = (visit_nbr pid new_prefix visited_n) in
                (visit_node id prefix) ;
                visited_n_p
        | Root(id, lid, rid) -> 
            let visited_n = (mark_visited node_id visited) in
            let new_prefix = (update_prefix prefix) in
            let visited_n_l = (visit_nbr lid new_prefix visited_n) in
                (visit_node id prefix) ;
                (visit_nbr rid new_prefix visited_n_l) 
        | Interior(id, pid, lid, rid) ->
            let visited_n = (mark_visited node_id visited) in
            let new_prefix = (update_prefix prefix) in
            let visited_p = (visit_nbr pid new_prefix visited_n) in
            let visited_n_l = (visit_nbr lid new_prefix visited_p) in
                (visit_node id prefix) ;
                (visit_nbr rid new_prefix visited_n_l) 
    in
        let _ = 
            (print_tree_aux id empty_string All_sets.Integers.empty) in
            ()

(** [print_forest forest]
    @return () - prints all the trees in the forest. *)
let print_forest forest = 
    let handles = (All_sets.Integers.elements (get_roots forest)) in
        (ignore (List.map (fun x -> print_newline () ;
                                    print_tree x forest) handles))
