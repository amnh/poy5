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

let () = SadmanOutput.register "Gen_rtree" "$Revision: 1644 $"

(* Module to represent the clusters of a generalized rooted tree. 
   Clusters are represented using BitSets. *) 

module Cluster = struct
    
    (** [universal_cluster size] 
        @param size size of the cluster. 
        @return a universal cluster of the given size. *)
    let universal_cluster size = 
        let uc = BitSet.create size in
            for i = 0 to size - 1 do
                BitSet.set uc i
            done ;
            uc 
    
    (** [disjoint c1 c2] 
        @param first cluster.
        @param second cluster.
        @return whether the two clusters are disjoint or not. *)
    let disjoint c1 c2 = 
        let intsxn = (BitSet.inter c1 c2) in
            ((BitSet.count intsxn) = 0)

    (** [complement c uc]
        @param c the cluster whose complement in the the passed in
            universal cluster is desired.
        @param size the size of the cluster.
        @return the complement of the cluster. *)
    let complement c size = 
        let copy = BitSet.copy c in
            for i = 0 to (size - 1) do
                BitSet.toggle copy i
            done ;
            copy
        
    (** [compatible c1 c2 size] 
        @param c1 the first cluster 
        @param c2 the second cluster
        @param size the size of the bit set. 
        @return whether the clusters are compatible or not. *)
    let compatible c1 c2 size =  
        if ( disjoint c1 c2 ) then
            true
        else 
            let nc1 = (complement c1 size) in
                if ( disjoint c2 nc1 ) then
                    true
                else 
                    let nc2 = (complement c2 size) in
                        if ( disjoint c1 nc2 ) then
                            true
                        else
                            (disjoint nc1 nc2)
    
    (** [includes c_sub c_sup]
        @param c_sub the cluster in question.
        @param c_sup does this cluster include c_sub?
        @param true if c_sup includes c_sub, false otherwise. *)
    let includes c_sub c_sup = 
        let diff = BitSet.diff c_sub c_sup in
            ((BitSet.count diff) = 0) 

    (** [merge clusters]
        @param clusters the list of clusters that need to be merged into one
                cluster. 
        @return function to merge the clusters into one. *)
    let merge clusters = 
        let fold_fn (acc_index, acc_size, acc_set) (index, size, set) =
            (-1, acc_size + size, BitSet.union acc_set set)
        in
        let (_, size, _) as cluster = 
            List.fold_left fold_fn (-1, 0, (BitSet.empty ())) clusters 
        in
            assert(size > 1) ;
            cluster
end

(* Indicates an invalid node id. *)
exception Invalid_Node_Id

(* A flexible node type to allow for trees with polytomies. *)
type node = Root of int * int list 
            (* (id, child id list) *)
          | Interior of int * int * int list 
            (* (id, parent id, child id list) *)
          | Leaf of int * int 
            (* (id, parent) *)

(* The forest of consensus trees or better still a consensus forest
* could be represented using this type. *)
type g_rtree = { g_rtopo : node All_sets.IntegerMap.t ;
                (* the topology of the consensus forest. *)
                roots : All_sets.Integers.t 
                (* the set of roots of the consensus forest. *) }

(* the empty consensus forest. *)
let empty = { g_rtopo = All_sets.IntegerMap.empty ;
              roots = All_sets.Integers.empty }

(*************************** GET FUNCTIONS ***************************)

(** [get_id node] 
    @param node the node whose id is desired
    @return id of the node. *)
let get_id node =
    match node with
    | Root(id, _) -> id
    | Leaf(id, _) -> id 
    | Interior(id, _, _) -> id 

(** [get_node id ctree]
    @param id id of the node desired.
    @param tree tree to which the node belongs
    @return returns the node associated with the id.
    @raise Invalid_Node_Id, if the id is not a valid id. *)
let get_node id tree = 
    try 
        (All_sets.IntegerMap.find id tree.g_rtopo)
    with 
        | Not_found -> raise Invalid_Node_Id

(** [get_parent nd_id ctree]
    @return the id of the parent node. *)
let get_parent nd_id ctree = 
    match get_node nd_id ctree with
    | Leaf(_, pid) -> pid 
    | Root(_, _) -> raise (Failure "Cannot get parent of root node.")
    | Interior(_, pid, _) -> pid
        
(** [is_node id ctree]
    @param id id of the node desired.
    @param tree tree to which the node belongs
    @return returns true when a node with the given id exists in the tree,
        returns false otherwise. *)
let is_node id tree = 
    (All_sets.IntegerMap.mem id tree.g_rtopo)
    
        
(** [is_root id ctree]
    @param id the id being checked whether it belongs to the root or not. 
    @param ctree the rooted tree in which id is either a root or not. *)
let is_root id ctree = 
    (* get the node corresponding to the id. *)
    match (get_node id ctree) with
    (* if it is a root node, make sure it exists in the root set and 
    * return true. *)
    | Root(_) -> 
        assert(All_sets.Integers.mem id ctree.roots) ;
        true 
    (* otherwise, return false. *)
    | _ -> false 

(** [is_leaf id ctree] 
    @param id the id of the node being checked to determine whether it
            is a leaf or not.
    @param ctree the rooted tree in which the id belongs to a leaf node or
            not. *)
let is_leaf id ctree = 
    (* get the node. *)
    match (get_node id ctree) with
    (* if the node is a leaf, then return true. *)
    | Leaf(_) -> true 
    (* otherwise false. *)
    | _ -> false 

(** [get_num_children nd_id g_rtree] 
    @param nd_id the id of the node whose number of children is desired.
    @param g_rtree the tree to which the node belongs.
    @return the number of children of this node. *)
let get_num_children nd_id g_rtree = 
    match get_node nd_id g_rtree with
    | Interior(_, _, chld_lst) 
    | Root(_, chld_lst) ->
        List.length chld_lst
    | Leaf(_, _) -> 0 


(** [other_children nd_id chld_id g_rtree]
    @param nd_id the node whose children other than chld_id are desired.
    @param chld_id the child that is to be ignored when returning the child
        list.
    @param g_rtree the tree in which all the children of nd_id, other than
        chld_id are desired.
    @return all the children of nd_id excepting chld_id. *)
let other_children nd_id chld_id g_rtree = 
    match get_node nd_id g_rtree with
    | Interior(_, _, chld_lst)
    | Root(_, chld_lst) ->
        List.filter (fun x -> (not (x = chld_id))) chld_lst
    | _ -> raise (Failure "Need interior or root node.")
    
(********************* BEGIN - UNSAFE OPS ****************************)

(** [add_node node ctree]
    @param node the node that needs to be added to the tree.
    @param tree the tree to which the node is being added.
    @return the tree with the node added.
    This function does not check whether the node is already
    present in the tree or not. If the node is already present then
    the tree is returned unmodified. *) 
let add_node node ctree = 
    let id = (get_id node) in
    let new_topo = (All_sets.IntegerMap.add id node ctree.g_rtopo) in
    let ctree = { ctree with g_rtopo = new_topo } in
    (* check whether the added node was a root. *)
    match node with
    | Root(_) -> 
        (* if it is, add it to the root set. *)
        { ctree with roots = All_sets.Integers.add id ctree.roots }
        (* else, do nothing. *)
    | _ -> ctree 

(** [remove_node nd_id ctree]
    @param nd_id id of the node that is being removed from the tree. 
    @param ctree tree from which the node is being removed. 
    @return the tree with the node removed. This function does not update
        any connectivity of the node's parent, children etc. All such
        connectivity information should be updated by the caller prior to
        calling this fn. *)
let remove_node nd_id ctree = 
    (* remove the node from the topology. *)
    let is_root = is_root nd_id ctree in
    let new_topo = (All_sets.IntegerMap.remove nd_id ctree.g_rtopo) in
    let ctree = { ctree with g_rtopo = new_topo } in
        (* if the removed node was a root, then, remove it from the 
        * set of roots. *)
        if ( is_root ) then
            { ctree with roots = All_sets.Integers.remove nd_id ctree.roots }
        else
            ctree

(** [add_nodes nodes ctree]
    @param nodes list of nodes that need to be added to the tree.
    @param btree tree to which the nodes are being added.
    @return tree with the list of nodes added. Older nodes with the
    same ids will be overwritten. *)
let add_nodes nodes ctree = 
    let f ctree node = (add_node node ctree) in
        (List.fold_left f ctree nodes)

(** [add_child nd_id child_id ctree]
    @return a ctree in which nd_id had a child child_id. *)
let add_child nd_id child_id ctree = 
    let node = (get_node nd_id ctree) in
    match node with
    | Root(_, chld_lst) ->
        (add_node (Root(nd_id, child_id :: chld_lst)) ctree)
    | Interior(_, pid, chld_lst) ->
        (add_node (Interior(nd_id, pid, child_id :: chld_lst)) ctree)
    | Leaf(_, pid) -> 
        (add_node (Interior(nd_id, pid, [child_id])) ctree)

(** [set_parent nd_id pid ctree]
    @return a ctree in which nd_id has parent set to pid. *)
let set_parent nd_id pid ctree = 
    let node = get_node nd_id ctree in
    match node with
    | Root(_, chld_lst) ->
        raise (Failure "Cannot set parent of root node.")
    | Interior(_, _, chld_lst) ->
        (add_node (Interior(nd_id, pid, chld_lst)) ctree)
    | Leaf(_, pid) -> 
        (add_node (Leaf(nd_id, pid)) ctree)

(** [remove_child nd_id chld_id ctree]
    @return a new ctree in which nd_id does not have chld_id as its child. *)
let remove_child nd_id chld_id ctree = 
    let node = get_node nd_id ctree in
    let filter_fn nd_id = (nd_id != chld_id) in
    match node with
    | Root(_, chld_lst) ->
        let n_chld_lst = List.filter filter_fn chld_lst in
            (add_node (Root(nd_id, n_chld_lst)) ctree) 
    | Interior(_, pid, chld_lst) ->
        let n_chld_lst = List.filter filter_fn chld_lst in
            (add_node (Interior(nd_id, pid, n_chld_lst)) ctree) 
    | Leaf(_, _) ->
        raise (Failure "Cannot replace child of leaf node.")
        
(********************** END - UNSAFE OPS *****************************)

(** [of_rtree rtree]
    @param rt_id the id of the root node.
    @param rtree the rooted tree. 
    @return a generalized rooted tree corresponding to the rooted tree. *)
let of_tree rtree = 
    (* function to convert an rtree node to a gen_rtree node. *)
    let convert_node rtree_nd = 
        match rtree_nd with
        | Rtree.Leaf (nd, p) -> Leaf(nd, p)
        | Rtree.Interior(nd, p, l, r) -> Interior(nd, p, [l; r])
        | Rtree.Root(nd, l, r) -> Root(nd, [l; r])
    in
    let { Rtree.roots = roots ; Rtree.r_topo = r_topo } = rtree in
    let g_rtree = { empty with roots = roots } in
    let g_rtree = { g_rtree with g_rtopo = 
                        All_sets.IntegerMap.map convert_node r_topo } in
        g_rtree 

(** [pre_order_node_visitor visitor_op node_id ctree accum]
    @param visit_op function to applied to each node as it is visited.
    @param node_id the node whose children are visited in pre-order. When
        node_id is root, then all the nodes of the tree are visited; 
        otherwise only the children of the node are visited.
    @param ctree the tree whose nodes are being visited.
    @param accum the accumulator that accumulates the results of all the 
                 visits. 
    @return the accumulator with the results of the node visit. *)
let pre_order_node_visitor visit_op node_id ctree accum = 
    (* function to visit the node. *)
    let visit node_id accum = 
        (visit_op node_id accum) 
    in
    (* function to visit the node and process the children. *)
    let rec visit_node_n_children node_id chld_lst accum =
        (* visit the node. *)
        let accum = (visit node_id accum) in
        (* visit the children. *)
        let accum = (List.fold_left traverse accum chld_lst) in
            (* return the accumulator. *)
            accum 
    (* function to traverse the tree in pre-order. *)
    and traverse accum node_id = 
        (* get the node corresponing to the start node id. *)
        match (get_node node_id ctree) with
        (* if it is the root node. *)
        | Root(root_id, chld_lst) -> 
            (* make sure the root is included in the roots set. *)
            assert(is_root root_id ctree) ;
            (* visit the node and its children. *)
            (visit_node_n_children root_id chld_lst accum)
        (* if it is an interior node. *)
        | Interior(int_id, _, chld_lst) ->
            (* visit the node and its children. *)
            (visit_node_n_children int_id chld_lst accum)
        (* if it is a leaf node *)
        | Leaf(leaf_id, _) -> 
            (* just, visit the node. *)
            (visit leaf_id accum) 
    in
        (* traverse the tree starting at the node_id. *)
        (traverse accum node_id)

(** [post_order_node_visitor visit_op node_id ctree accum]
    @param visit_op function to applied to each node as it is visited.
    @param node_id the node whose children are visited in pre-order. When
        node_id is root, then all the nodes of the tree are visited; 
        otherwise only the children of the node are visited.
    @param ctree the tree whose nodes are being visited.
    @param accum the accumulator that accumulates the results of all the 
                 visits. 
    @return the accumulator with the results of the node visit. *)
let post_order_node_visitor visit_op node_id ctree accum = 
    (* function to visit the node. *)
    let visit node_id accum = 
        (visit_op node_id accum) 
    in
    (* function to visit the node and process the children. *)
    let rec visit_node_n_children node_id chld_lst accum =
        (* visit the children. *)
        let accum = (List.fold_left traverse accum chld_lst) in
        (* visit the node. *)
        let accum = (visit node_id accum) in
            (* return the accumulator. *)
            accum 
    (* function to traverse the tree in pre-order. *)
    and traverse accum node_id = 
        (* get the node corresponing to the start node id. *)
        match (get_node node_id ctree) with
        (* if it is the root node. *)
        | Root(root_id, chld_lst) -> 
            (* make sure the root is included in the roots set. *)
            assert(is_root root_id ctree) ;
            (* visit the node and its children. *)
            (visit_node_n_children root_id chld_lst accum)
        (* if it is an interior node. *)
        | Interior(int_id, _, chld_lst) ->
            (* visit the node and its children. *)
            (visit_node_n_children int_id chld_lst accum)
        (* if it is a leaf node *)
        | Leaf(leaf_id, _) -> 
            (* just, visit the node. *)
            (visit leaf_id accum) 
    in
        (* traverse the tree starting at the node_id. *)
        (traverse accum node_id) 

(** [get_num_nodes ctree]
    @return the number of nodes in the subtree rooted at root_id. *)
let get_num_nodes root_id ctree = 
    let visit_op node_id (num_leaves, num_nodes) = 
        if (is_leaf node_id ctree) then
            (num_leaves + 1, num_nodes + 1)
        else
            (num_leaves, num_nodes + 1)
    in
        (pre_order_node_visitor visit_op root_id ctree (0, 0)) 

(** [print_node node]
    @return prints the node. *)
let print_node node =
    let print_lst lst = 
        let print_fn item = 
            print_string ((string_of_int item) ^ " ") in
            List.iter print_fn lst
    in
    match node with
    | Leaf(id, pid) ->
         print_endline ("Leaf(" ^ (string_of_int id) ^ ", " ^
                               (string_of_int pid) ^ ")") 
    | Root(id, chld_lst) ->
         print_string ("Root(" ^ (string_of_int id) ^ ", ") ;
         print_lst chld_lst ;
         print_endline ")" 
    | Interior(id, pid, chld_lst) ->           
         print_string ("Interior(" ^ (string_of_int id) ^ ", ") ;
         print_string ((string_of_int pid) ^ " ") ;
         print_lst chld_lst ;
         print_endline ")" 

(** [move_root nd_id g_rtree]
    @param nd_id the node to which the root is being moved. 
    @param g_rtree the generalized rooted tree/forest in which the
                   root is being moved to nd_id. 
    @return the g_rtree with the root moved to nd_id. If the root had
            only two children then it would be contracted out of the tree. *)
let move_root nd_id g_rtree = 
    (* function to fix the parent child relationships after a re-root
    * operation. *)
    let rec fix_tree nd_id prev_pid g_rtree = 
        (* get the previous parent. *)
        match get_node prev_pid g_rtree with
        (* if it is an interior node. *)
        | Interior(_, prev_ppid, chld_lst) ->
            (* remove nd_id from previous parent's child list. *)
            let g_rtree = (remove_child prev_pid nd_id g_rtree) in
            (* add previous parent's parent to the child list. *)
            let g_rtree = (add_child prev_pid prev_ppid g_rtree) in
            (* set the nd_id as the parent of prev parent. *)
            let g_rtree = (set_parent prev_pid nd_id g_rtree) in
                (* continue fixing the tree untill you run into the 
                * previous root. *)
                (fix_tree prev_pid prev_ppid g_rtree) 
        (* if it is a root. *)
        | Root(_, chld_lst) ->
            (* get the number of children. *)
            let num_children = (get_num_children prev_pid g_rtree) in
            (* if the root has 2 children, contract the root out. *)
            if ( num_children = 2 ) then
                (* get the other child of the root node. *)
                let o_childs = (other_children prev_pid nd_id g_rtree) in
                assert((List.length o_childs) = 1) ;
                let o_c = List.hd o_childs in
                (* remove the root from nd_id's child list. *)
                let g_rtree = (remove_child nd_id prev_pid g_rtree) in
                (* set the other child's parent to nd_id. *)
                let g_rtree = (set_parent o_c nd_id g_rtree) in
                (* add other child as a child to nd_id. *)
                let g_rtree = (add_child nd_id o_c g_rtree) in
                (* remove the old root. *)
                let g_rtree = (remove_node prev_pid g_rtree) in
                    g_rtree 
            (* if the root has more than 2 children, convert it into an
            * interior node. *)
            else if ( num_children > 2 ) then
                let g_rtree = 
                    (add_node (Interior(prev_pid, nd_id, chld_lst)) g_rtree) in
                    (remove_child prev_pid nd_id g_rtree) 
            else
                (* if the root has less than two children, then it is
                * an error. *)
                raise (Failure "Root has less than two children.")
        | _ -> raise (Failure "Expecting an interior or root node.")
    in
    (* get the node *)
    match get_node nd_id g_rtree with
    (* if it is an interior node, replace it with a root node. *)
    | Interior(_, pid, chld_lst) ->
        (* note that the previous parent now becomes a child. *)
        let g_rtree = (add_node (Root(nd_id, pid :: chld_lst)) g_rtree) in
            (* fix the parent-child relationships after a reroot. *)
            (fix_tree nd_id pid g_rtree)
    (* if the node is already a root, do nothing. *)
    | Root(_) -> g_rtree 
    | _ -> raise (Failure "expecting an interior node.")

(** [set_out_group out_grp g_rtree]
    @param out_grp the out group. 
    @param g_rtree the tree in which out_grp will be made the out group. 
    @return the tree with root moved to the only nbr of the out group. *)
let set_out_group out_grp g_rtree = 
    match get_node out_grp g_rtree with
    | Leaf(_, nbr) -> 
         (nbr, (move_root nbr g_rtree))
    | _ -> raise (Failure "Expecting leaf node.") 
    
(** [attempt_ctree rtree_nd_hash_lst index_tbl sp_tbl]
    @param rtree_nd_hash_lst the rooted tree and corresponding node_hash 
        list.
    @param index_tbl the indexes of the leaf nodes.
    @param sp_tbl the hash_table with the split counts.
    @return the consensus tree of the trees in the rtree_nd_hash_lst. *)
let attempt_ctree maj num_nodes rtree_nd_hash_lst sp_tbl =
    (* function to generate the ids of the nodes in the ctree. *)
    let id_fn = 
        let x = ref ((num_nodes + 1) / 2) in
            (fun () -> (incr x; !x))
    in
    (* The table of node hashes to the index in the ctree. Used to 
    * determine whether a node already exists in the ctree or not. *)
    let nd_hashes_index_tbl = Hashtbl.create num_nodes in
    (* Hash table to store the size of paritions at consensus tree
    * nodes. *)
    let nd_size_tbl = Hashtbl.create num_nodes in
    (* function to determine whether a node exists in the ctree 
    * corresponding to the node_hashes. *)
    let has_node (nd_hash, sp_hash) = 
        (Hashtbl.mem nd_hashes_index_tbl (nd_hash, sp_hash)) 
    in
    (* function to get the index in the ctree, given node and 
    * split hashes. *)
    let get_node_index (nd_hash, sp_hash) = 
        (Hashtbl.find nd_hashes_index_tbl (nd_hash, sp_hash))
    in
    (* function to save a node's index in the ctree. The node is 
    * identified using its node and split hashes. *)
    let save_node_index (nd_hash, sp_hash) index = 
        (Hashtbl.add nd_hashes_index_tbl (nd_hash, sp_hash) index)
    in
    (* function to incorporate the data of one input tree into the 
    * consensus tree. *)
    let update_consensus_with_tree (croot_id, ctree)
                                   (root_id, rtree, nd_hashes) = 
        (* function to process the node. *)
        let visit_op nd_id (pid, ctree) = 
            let (nd_hash, sp_hash, _) = Hashtbl.find nd_hashes nd_id in
            let (count, size, index) = 
                Hash_tbl.get_data (nd_hash, sp_hash) sp_tbl in
            (* if the node is a majority node. *)
            if ( count >= maj ) then begin
                (* Is the node in the consensus tree? *)
                if (has_node (nd_hash, sp_hash)) then begin
                    let index = get_node_index (nd_hash, sp_hash) in
                    (* update the parent. *)
                    if (not (is_root index ctree)) then begin
                        let cpid = (get_parent index ctree) in
                        let cpsize = Hashtbl.find nd_size_tbl cpid in
                        let psize = Hashtbl.find nd_size_tbl pid in
                            match psize < cpsize with
                            | true ->
                                (* better parent found. *)
                                assert(not (pid = cpid)) ;
                                let ctree = (set_parent index pid ctree) in
                                let ctree = (remove_child cpid index ctree) in
                                let ctree = (add_child pid index ctree) in
                                    (index, ctree)
                            | false ->
                                (* do nothing. *)
                                (index, ctree)
                    end else
                        (index, ctree)
                end else begin
                    (* create a new node for this majority node. *)
                    assert(size > 0) ;
                    match size with
                    | 1 -> 
                        let _ = (save_node_index (nd_hash, sp_hash) index) in
                        let ctree = (add_node (Leaf(index, pid)) ctree) in
                        let ctree = (add_child pid index ctree) in
                        let _ = Hashtbl.add nd_size_tbl index 1 in
                            (index, ctree)
                    | _ ->
                        let index = id_fn () in
                        let _ = (save_node_index (nd_hash, sp_hash) index) in
                        let ctree = 
                            (add_node (Interior(index, pid, [])) ctree) in
                        let ctree = (add_child pid index ctree) in
                        let _ = Hashtbl.add nd_size_tbl index size in
                            (index, ctree)
                end
            end else
                (* do nothing. *)
                (pid, ctree)
        in
        (* function to perform a pre-order-traversal of the rtree. *)
        let rec traverse nd_id rtree (pid, ctree) =
            (* visit the node. *)
            let (pid, ctree) = visit_op nd_id (pid, ctree) in
            (* visit its children. *)
            match Rtree.get_node nd_id rtree with
            | Rtree.Root(_, lid, rid) 
            | Rtree.Interior(_, _, lid, rid) ->
                let _, ctree = traverse lid rtree (pid, ctree) in
                let pid, ctree = traverse rid rtree (pid, ctree) in
                    pid, ctree
            | Rtree.Leaf(_, _) ->
                (pid, ctree) 
        in
        (* make sure root_id is a root id. *)
        assert(Rtree.is_root root_id rtree) ;
        (* visit the root. *)
        let (_, ctree) = traverse root_id rtree (croot_id, ctree) in
            (croot_id, ctree)
    in
        (* Build the root of the consensus tree. *)
        match rtree_nd_hash_lst with
        | [] -> (-1, empty, nd_size_tbl)
        | (root_id, rtree, nd_hashes) :: _ -> 
            (* get the hash for the root_id. *)
            let (nd_hash, sp_hash, sp_size) = Hashtbl.find nd_hashes root_id in
            let croot_id = (id_fn ()) in
            let ctree = (add_node (Root(croot_id, [])) empty) in
            let _ = Hashtbl.add nd_size_tbl croot_id sp_size in
            let _ = (save_node_index (nd_hash, sp_hash) croot_id) in
            (* Build the rest of the consensus tree. *)
            let croot_id, ctree =
                List.fold_left update_consensus_with_tree 
                               (croot_id, ctree) 
                               rtree_nd_hash_lst in
                (croot_id, ctree, nd_size_tbl) 

(** [verify_ctree (croot_id, ctree, nd_size_tbl)] 
    @param croot_id the root of the ctree.
    @param ctree the consensus tree.
    @nd_size_tbl the table containing the sizes of the partitions
        of each node in the original-hash-table built by consensus
        tree algorithm.
    @return a boolean to indicate whether the consensus tree was built 
        correctly or has a double-collision. *) 
let verify_ctree (croot_id, ctree, nd_size_tbl) = 
    (* function to make a post-order-traversal and verify the 
    * consensus tree. *)    
    let rec traverse nd_id = 
        (* get the node. *)
        let node = get_node nd_id ctree in
        match node with 
        (* if the node is a root or an interior node, compute the size of the
        * split and check against the nd_size_tbl. *)
        | Root(_, chld_lst)
        | Interior(_, _, chld_lst) ->
            (* get the results. *)
            let res = List.map traverse chld_lst in
            let statuses, sizes = List.split res in
            let status = List.fold_left (&&) true statuses in
                if (not status) then
                    (raise Hash_tbl.Double_Collision)
                else
                    let actual_size = List.fold_left (+) 0 sizes in
                    let expected_size = (Hashtbl.find nd_size_tbl nd_id) in
                        (actual_size = expected_size, actual_size)
        (* if the node is a leaf, make sure it has size 1. *)
        | Leaf(_, _) -> 
            assert((Hashtbl.find nd_size_tbl nd_id) = 1) ;
            (true, 1) 
    in
        try
            (* Verify the tree. *)
            let (status, _) = traverse croot_id in
                assert(status) ;
                status
        with
        (* Got double collison, return false. *)
        | Hash_tbl.Double_Collision -> false
        (* Unexpected exn, raise it as is. *)
        | exn -> raise exn


(** [build_ctree maj num_nodes rtree_lst]
    @param maj the number of trees in which a partition must appear, in order 
    to appear in the final tree, 
    @param num_nodes number of nodes in each of the rooted trees. Note that
    all these trees need to have the same number of nodes.
    @param rtree_lst = (root_id, rtree) the list of the rooted trees on the 
    same set of leaves. 
    @return A consensus tree of all the trees in the rtree_lst is returned. *)
let build_ctree  maj num_nodes rtree_lst = 
    (* maximum number of tries before a good consensus tree is found. *)
    let max_tries = 150 in
    (* Initialize the random number generator. *)
    let _ = Random.self_init () in
    (* This factor determines the probability of a double collision. 
    * prob = 1 / safety_factor. Large safety factor is good, but it could lead
    * to overflow issues. *)
    let safety_factor = 500 in
    (* get the number of trees. *)
    let num_trees = List.length rtree_lst in
    (* make sure maj is correct. *)
    if (maj <= num_trees / 2) then begin
        Status.user_message Status.Error ("You@ have@ requested@ a@ consensus@ "
        ^ "with@ majority@ below@ 0.5.@ Are@ you@ sure@ you'd@ like@ to@ do@ "
        ^ "that?@ That@ is@ no@ majority@ at@ all!.@ I@ can@ not@ calculate@ "
        ^ "such@ consensus.");
        failwith "Consensus: Illegal argument";
    end;
    if (num_trees <= 1) then begin
        Status.user_message Status.Information ("There@ " ^
        (if num_trees = 1 then 
            "is@ " ^ string_of_int num_trees ^ " tree@ "
        else "are@ no@ trees@ ")
        ^ "in@ memory.");
        failwith "Consensus: Illegal argument";
    end;
    (* make sure there are atleast two trees. *)
    assert(num_trees > 1) ;
    (* get the prime numbers m1 and m2 of the hash functions. *)
    let total_nodes = (num_nodes * num_trees) in
    let m1 = Primes.Probable.next_prime_gt total_nodes in
    (* potential overflow. *)
    let m2 = Primes.Probable.next_prime_gt (safety_factor * total_nodes) in
    (* get the number of leaves. *)
    let num_leaves = (num_nodes + 1) / 2 in
    (* create the index table. *)
    let (rt_id, rtree) = List.hd rtree_lst in
    let indx_nd_tbl, nd_indx_tbl = 
        Rtree.index_nodes (0, num_leaves) num_nodes rt_id rtree in
    (* Attempt to create a consensus tree. *)
    let rec attempt trial_num = 
        (* create the parameters for the hash_table. *)
        let h1_p = Array.make num_leaves 0 in
        let h2_p = Array.make num_leaves 0 in
        let _ = for i = 0 to num_leaves - 1 do
                    h1_p.(i) <- (Random.int m1) ;
                    h2_p.(i) <- (Random.int m2)
                done 
        in
        (* create the actual split-hash table. *)
        let sp_hash_tbl = Hash_tbl.make h1_p h2_p m1 m2 in
        (* populate the split-hash table. *)
        let populate (rt_id, rtree) = 
            (* create the node_hashes tbl. *)
            let node_hashes = Hashtbl.create num_nodes in
            let node_hashes, _ = Hash_tbl.populate rt_id 
                                  (rtree, nd_indx_tbl, node_hashes)
                                  sp_hash_tbl in
                (rt_id, rtree, node_hashes) 
        in
        (* Attempt to build the ctree. *)
        try
            let rt_nd_hash_lst = List.map populate rtree_lst in
            let croot_id, ctree, nd_sz_tbl = 
                attempt_ctree maj num_nodes rt_nd_hash_lst sp_hash_tbl in
            let verified = (verify_ctree (croot_id, ctree, nd_sz_tbl)) in
                begin
                    match verified with
                    | true -> (croot_id, ctree, indx_nd_tbl, nd_indx_tbl) 
                    | false -> 
                        if ( trial_num < max_tries ) then
                            (attempt (trial_num + 1)) 
                        else
                            raise (Failure "max_tries exceeded.")
                end
        with
        | Hash_tbl.Double_Collision -> 
                Status.user_message Status.Error ("Hrm,@ a@ hashtable@ "
                ^ "collition@ while@ calculating@ the@ consensus.@ "
                ^ "This@ is@ normal,@ I@ will@ try@ again@ "
                ^ "(up@ to@ 150@ times,@ so@ don't@ get@ desperate).");
                (attempt (trial_num + 1))
        | xn -> (raise xn)
    in
        (attempt 1)

(** [index_nodes start_index root_id rtree]
    @param start_index the index at which the indexing starts.
    @param root_id identifies the tree 
    @param rtree the generalized rooted forest.
    @return a hash-table where the indexes are hashed with node_ids as keys. 
*) 
let index_nodes (li, ii) num_nodes root_id rtree = 
    (* Function to add the current node_id to the index table. *)
    let visit_op node_id ((cli, cii), index_tbl) =
        (* if the node is a leaf *)
        if (is_leaf node_id rtree) then begin
            (* add it to the table. *)
            (Hashtbl.add index_tbl node_id cli) ;
            (* increment leaf index. *)
            ((cli + 1, cii), index_tbl)
        (* else, if it is an interior node or root node. *)
        end else begin
            (* add it to the table. *)
            (Hashtbl.add index_tbl node_id cii) ;
            (* increment the interior index. *)
            ((cli, cii + 1), index_tbl)
        end
    in
    let index_tbl = (Hashtbl.create num_nodes) in
    let (_, index_tbl) = (pre_order_node_visitor visit_op root_id rtree
                        ((li, ii), index_tbl)) in
        index_tbl 

(** [get_cluster_lst index_tbl (rt_id, g_rtree)]
    @param index_tbl the table that indexes the nodes such that leaf_nodes 
            start at 0. 
    @param num_leaves the number of leaves in the tree. 
    @param rt_id the id of the root node. 
    @param g_rtree the generalized rooted forest. 
    @return the list of clusters of the tree. *)
let get_cluster_lst index_tbl num_leaves (rt_id, g_rtree) = 
    (* the list of the clusters. *)
    let cluster_lst = ref [] in
    (* function to populate the cluster table. *)
    let rec traverse nd_id = 
        match (get_node nd_id g_rtree) with
        | Leaf(_, _) ->
            let index = (Hashtbl.find index_tbl nd_id) in
            let leaf_set = BitSet.create num_leaves in
            let _ = BitSet.set leaf_set index in
            let cluster = (index, 1, leaf_set) in
                cluster_lst := cluster :: (!cluster_lst) ;
                cluster
        | Interior(_, _, chld_lst) 
        | Root(_, chld_lst) ->
            let chld_clusters = List.map traverse chld_lst in
            let cluster = Cluster.merge chld_clusters in
                cluster_lst := cluster :: (!cluster_lst) ;
                cluster
    in
    let _ = (traverse rt_id) in
        !cluster_lst

(** [get_cluster_lsts index_tbl leaf_hashes <(rt_id, g_rtree)>]
    @param index_tbl the hash-tbl that contains the indexes of nodes. 
    @param leaf_hashes the hashes of the leaf-nodes. 
    @param <rt_id, g_rtree> list of the generalized rooted trees and their
            nodes. *)
let get_cluster_lsts index_tbl leaf_hashes g_rtree_lst = 
    let get_clusters (rt_id, g_rtree) = 
        get_cluster_lst index_tbl leaf_hashes (rt_id, g_rtree)
    in
        List.map get_clusters g_rtree_lst

(** [get_pairwise_compatible_clusters cluster_lsts]
    @param num_leaves the number of leaves in the trees to which these
            clusters correspond. 
    @param cluster_lsts the list of cluster lists' of each tree.
    @return a list of pairwise consistent clusters. *)
let get_pairwise_compatible_clusters num_leaves cluster_lsts = 
    (* function to filter out the leaves and root nodes. *)
    let filter_leaves_n_root (_, size, _) =
        not ((size = 1) || (size = num_leaves))
    in
    (* function to keep root and leaves *)
    let keep_leaves_n_root (_, size, _) = 
        ((size = 1) || (size = num_leaves))
    in
    (* function to keep only the root.
    let keep_root (_, size, _) = 
        (size = num_leaves)
    in
    *)
    (* function to get the list of all interior clusters of the 
    * trees in a single list. *)
    let folder_fn accum cluster_lst = 
        let filtered_lst = List.filter filter_leaves_n_root cluster_lst in
            List.rev_append filtered_lst accum 
    in
    (* the leaf and root clusters. *)
    let root_n_leaf_clusters = 
        List.filter keep_leaves_n_root (List.hd cluster_lsts)
    in
    (* long list of interior clusters, only these could potentially be
    * incompatible with each other. *)
    let interior_clusters = List.fold_left folder_fn [] cluster_lsts in
    (* function to check whether the cluster is compatible with
    * every item in the given list. *)
    let rec compatible_with_all_items ((_, _, hash) as cluster) cluster_lst = 
        match cluster_lst with
        | [] ->
            true 
        | (c_index, c_size, c_hash) :: rest ->
            begin
                match Cluster.compatible hash c_hash num_leaves with
                | true ->
                    (compatible_with_all_items cluster rest)
                | false ->
                    false
            end
    in
    (* function to compute the pair-wise compatible interior clusters. *)
    let rec pairwise_compatible cluster_lst results = 
        match cluster_lst with
        | [] -> results
        | hd :: rest ->
            begin
                match compatible_with_all_items hd interior_clusters with
                | true -> 
                    (pairwise_compatible rest (hd :: results))
                | false -> 
                    (pairwise_compatible rest results)
            end
    in
    let pairwise_compatible_interior_clusters = 
        (pairwise_compatible interior_clusters []) in
        List.rev_append root_n_leaf_clusters
                        pairwise_compatible_interior_clusters 

(** [build_semi_strict_consensus leaf_hashes clusters]
    @param leaf_hashes the hashes of the leaf_nodes in the original trees.
    @param clusters the semi-strict consensus clusters. 
    @return the semi-strict consensus tree. *)
let build_semi_strict_consensus num_leaves clusters = 
    (* get the number of nodes. *)
    let num_nodes = List.length clusters in
    (* the node cluster sets needed to build the ssct. *)
    let node_cluster_sets = Array.init num_nodes 
        (fun x -> (BitSet.empty ())) 
    in
    (* function to generate the ids of the nodes in the ctree. *)
    let id_fn = 
        let x = ref num_leaves in
            (fun () -> (let ret = !x in incr x; ret))
    in
    (* sort the clusters in descending order of size. *)
    let compare (_, size1, _) (_, size2, _) = 
        (size2 - size1)
    in
    (* function to sort the clusters in descreasing order. *)
    let sorted_clusters = List.sort compare clusters in
    match sorted_clusters with
    | (index, size, hash) :: rest -> 
        (* make sure the first cluster is the root. *)
        assert(size = num_leaves) ;
        let rt_id = id_fn () in
        let ssct = add_node (Root(rt_id, [])) empty in
        let root_cluster_set = Cluster.universal_cluster num_leaves in
        let _ = node_cluster_sets.(rt_id) <- root_cluster_set in
        (* function to add cluster to semi-strict consensus tree. *)
        let rec add_cluster (index, size, cluster_set) (nd_id, ssct) = 
            (* get the node's cluster set. *)
            let node_cluster_set = node_cluster_sets.(nd_id) in
                (* if the cluster is included in the node. *)
                if Cluster.includes cluster_set node_cluster_set then
                begin
                    match (get_node nd_id ssct) with
                    | Root(_, chld_lst)
                    | Interior(_, _, chld_lst) ->
                        (* attempt to add cluster to any of the child nodes. *)
                        let status, ssct =
                            try_cluster_with_child_nodes 
                                        (index, size, cluster_set)
                                                         chld_lst
                                                         ssct in
                            (* if the cluster could be added. *)
                            if status then
                                (* return the tree. *)
                                true, ssct
                            (* otherwise, add the cluster as a new child. *)
                            else begin
                                (* if the cluster has size 1, then it is a leaf.
                                * *)
                                match size with
                                | 1 -> 
                                    let ssct = 
                                        (add_node (Leaf(index, nd_id)) ssct) in
                                    let ssct = (add_child nd_id index ssct) in
                                        true, ssct
                                (* otherwise it is an interior cluster. *)
                                | n -> 
                                    let index = id_fn () in
                                    let ssct = 
                                        (add_node 
                                            (Interior(index, nd_id, [])) ssct)
                                    in
                                    let ssct = (add_child nd_id index ssct) in
                                    let _ = 
                                        node_cluster_sets.(index) <- 
                                            cluster_set in
                                        true, ssct
                            end
                    | Leaf(_, _) -> 
                        raise (Failure "cluster couldn't be added.")
                (* if the cluster is not included, return false status and 
                * the tree built so far. *)
                end else
                    false, ssct
        (* function to attempt to add cluster to one of the child nodes. *)
        and try_cluster_with_child_nodes (index, hash, size) chld_lst ssct =
        begin
            match chld_lst with
            (* if the child list is empty then the the cluster could not 
            * be added to any of the children. *)
            | [] -> 
                (false, ssct)
            | nd_id :: rest ->
                let status, ssct = 
                    add_cluster (index, hash, size) (nd_id, ssct) in
                    if status then
                        status, ssct
                    else
                        try_cluster_with_child_nodes (index, hash, size) 
                                                     rest 
                                                     ssct
        end in
        (* function to add clusters to the semi-strict consensus tree. *)
        let rec add_clusters clusters (rt_id, ssct) = 
        begin
            match clusters with 
            (* if the clusters list is empty, return the semi-strict consensus
            * tree built so far. *)
            | [] -> ssct 
            (* otherwise, add a cluster to the semi-strict consensus tree built
            * so far. then, recursively add the rest. *)
            | cluster :: rest -> 
                (* add the cluster to the ssct. *)
                let status, ssct = add_cluster cluster (rt_id, ssct) in
                    assert(status) ;
                    (add_clusters rest (rt_id, ssct))
        end in
            (add_clusters rest (rt_id, ssct))
    | _ -> raise (Failure "found no clusters.")

(** [build_semi_strict_ctree out_grp g_rtrees]
    @param out_grp the out group where all the trees will be matched at.
    @param g_rtrees the generalized rooted trees whose semi-strict consensus is
        needed. 
    @return the semi-strict consensus tree of the input set of trees. *)
let build_semi_strict_ctree out_grp g_rtrees = 
    (* get the number of leaves and nodes. *)
    let (rt_id, g_rtree) = (List.hd g_rtrees) in
    let (num_leaves, num_nodes) = get_num_nodes rt_id g_rtree in
    (* root all the trees at the out_grp. *)
    let g_rtrees = 
        List.map (fun (rt_id, g_rtree) -> (set_out_group out_grp g_rtree))
                 g_rtrees
    in
    (* index the leaf nodes. *)
    let index_tbl = index_nodes (0, num_leaves) num_nodes rt_id g_rtree in
    (* get the cluster lists. *)
    let cluster_lsts = get_cluster_lsts index_tbl num_leaves g_rtrees in
    (* get pairwise compatible clusters. *)
    let pwcc_lst = get_pairwise_compatible_clusters num_leaves cluster_lsts in
    (* build the semi-strict consensus tree. *)
    let ssct = build_semi_strict_consensus num_leaves pwcc_lst in
        ssct 

(* Construct a parser friendly consensus tree *)
let to_parser c_topo to_string =
    let rec builder item = 
        match get_node item c_topo with
        | Root (_, children) ->
                let res = List.map builder children in
                Parser.Tree.Node (res, to_string item)
        | Interior (_, _, children) ->
                let res = List.map builder children in
                Parser.Tree.Node (res, to_string item)
        | Leaf (_, _) ->
                Parser.Tree.Leaf (to_string item)
    in
    let roots = All_sets.Integers.elements c_topo.roots in
    let trees = List.map builder roots in
    Parser.Tree.Node (trees, "")
