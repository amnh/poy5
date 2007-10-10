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

let () = SadmanOutput.register "UnionTree" "$Revision: 1616 $"

(* This module handles operations on trees with Node.Union.n data in their
* vertices *)

let debug = false 

let union_size = 10

module type S = sig

    (** The type of a node in a tree *)
    type n 

    (* The type of an edge in the tree *)
    type e 

    (** The type of a node union *)
    type u

    (** The data holding the union of a set of Node.node_data vertices *)
    type node

    val leaf : int list -> int -> u -> node

    (** A comparison function for total ordering of trees. The function 
    * complies with the specifications required by Pervasives.compare.*)
    val compare_trees_data : 
      (node, 'b) Ptree.p_tree ->
          (node, 'b) Ptree.p_tree ->
              int list

    (* [create_initial_tree h x y] 
     * Creates an initial tree out of a regular phylogenetic tree [y],
     * optionally only in the connected components with handle [h], 
     * restarting a new union cluster whenever the saturation of 
     * polymorphic positions is greaier than [x]. Outputs a tuple [(a, b)], 
     * whera [a] is the new union holding tree with the same topology as [x], 
     * and the map [b] of vertices that hold unions of clusters of vertices. 
     * *)
    val create_initial_tree : int option ->
      float ->
      (n, e) Ptree.p_tree ->
      (node, e) Ptree.p_tree * node All_sets.IntegerMap.t 

    (* [update v w x y z] 
    * Updates the vertices of the union tree down to the root, defined in the
    * set [z]. It is a requirement that the vertices in [z] form a connected
    * component that includes the handle of the tree [x].
    * The function uses the original phylogenetic tree [w],
    * making clusters with saturation of polymorphic positions below [v], and
    * updating the map of unions [y]. The function returns a quadruple 
    * [(a, b, c, d)], where [a] is the new union tree, [b] is the new map of 
    * union vertices to the union data, [c] is the vertices that have been 
    * added to the union set, and [d] the vertices that where removed from 
    * the union set. *)
    val update :
        int option 
        -> float -> (n, e) Ptree.p_tree -> (node, e) Ptree.p_tree 
        -> node All_sets.IntegerMap.t -> All_sets.Integers.t 
        -> (node, e) Ptree.p_tree * node All_sets.IntegerMap.t *
            All_sets.Integers.t * All_sets.Integers.t

    (* Calculates the distance between two nodes *)
    val distance : node -> node -> float

    (* [get_union_vertex x y z]
    * Find the code of the vertex in [y] that holds in it's union the vertex 
    * [x]. The set of vertices is contained in the map [y]. *)
    val get_union_vertex : 
        int -> (node, e) Ptree.p_tree -> node All_sets.IntegerMap.t -> int

    (* Some useful functions to deal with clusters of vertices *)
    module Clusters : sig
        (* The type of the maps of clusters, from vertices to clusters *)
        type c = node All_sets.IntegerMap.t

        (* [distance a b c] is the distance between the union vertices [b] 
         * and [c] using the unions defined in the map [a] *)
        val distance : c -> int -> int -> float

        (* [update_matrix w x y z] updates the contents of the distance matrix
        * [w] using the information of the union map [x], when the set [y] 
        * are new vertices in [x] and [z] are removed vertices in [x]. The 
        * function returns an updates matrix, where the members of [z] have 
        * been added and the members of [y] have been removed. *)
        val update_matrix : SparceMatrix.m -> c -> All_sets.Integers.t ->
            All_sets.Integers.t -> SparceMatrix.m

        (* [initial_matrix x] creates the initial all to all distance matrix 
        * for the union map [x] *)
        val initial_matrix : c -> SparceMatrix.m

        module Test : sig
            val unity : SparceMatrix.m -> c -> All_sets.Integers.t ->
                All_sets.Integers.t -> bool
        end
    end

end

module Codes = All_sets.Integers
module IntMap = All_sets.IntegerMap

let (-->) a b = b a

type 'a node1 = {
        union : 'a;                   
            (* The union *)
        children_unions : All_sets.Integers.t; 
            (* The vertices holding the next unions *)
        code : int;
            (* The code of the vertex that holds this node *)
        size : int;
    }

module Make (Node : NodeSig.S) (Edge : Edge.EdgeSig) : 
    S with type n = Node.n with type u = Node.Union.u 
    with type node = Node.Union.u node1 with type e = Edge.e = struct 

        type n = Node.n
        type u = Node.Union.u
        type e = Edge.e

        type node = u node1

        let create a b c z =
            { 
                union = a;
                children_unions = b;
                code = c;
                size = z;
            }

        let leaf items code a = {
            union = a;
            children_unions = 
                List.fold_left (fun acc x -> Codes.add x acc)
                Codes.empty
                items;
            code = code;
            size = 1;
        }

        let compare_trees_data a b =
            let comparator x y acc =
                let bitem = IntMap.find x b.Ptree.node_data in
                match Node.Union.compare y.union bitem.union with
                | 0 -> acc
                | _ -> x :: acc
            in
            IntMap.fold comparator a.Ptree.node_data []

        (* A function to add a union assigned to a single vertex to a cluster;
        * the vertex is the root of the cluster  *)
        let add_to_cluster = IntMap.add 

        (* Create a new pair of clusters, given that the saturation of the 
        * current cluster has reached more than the acceptable value *)
        let handle_new_cluster a b ca cb clusters =
            add_to_cluster b cb (add_to_cluster a ca clusters)

        (* We define a function that can handle a single vertex at a time, and
        * recursively process the complete subtree rooted by the [vertex]
        * argument and create a union of it's two children and itself. 
        * It returns a tuple, [((c, d), b)], where [c] is a map of vertices 
        * and their associated union values, [d] is the set of vertices that 
        * conform a the root of a cluster, and [d] is the union created for 
        * the processed [vertex]. 
        * If the function finds that the data of a particular vertex has already
        * been calculated, then it will not attempt to continue that subtree and
        * assumes that such information is correct *)
        let rec process_subtree saturation_parameter ptree 
            ((acc, clusters) as allacc) vertex =
            if IntMap.mem vertex acc then allacc, (IntMap.find vertex acc)
            else
                let my_node = Ptree.get_node_data vertex ptree in
                match Ptree.get_node vertex ptree with
                | Tree.Interior (_, par, a, b) ->
                        let (allacc, child1) = 
                            process_subtree saturation_parameter ptree allacc a
                        in
                        let ((acc, clusters), child2) = 
                            process_subtree saturation_parameter ptree allacc b 
                        in
                        if union_size < child1.size + child2.size then
                            let res = 
                                let code = 
                                    Node.min_child_code (Some par) my_node 
                                in
                                my_node 
                                --> Node.Union.leaf (Some code) (Some par) 
                                --> (fun x -> 
                                    Node.Union.union_final (Some par) 
                                    x my_node)
                                --> leaf [a; b] vertex
                            and clusters = 
                                handle_new_cluster a b child1 child2 clusters 
                            in
                            ((IntMap.add vertex res acc), clusters), res
                        else 
                            let res = 
                                let res = 
                                    Node.Union.union (Some par) my_node 
                                    child1.union child2.union 
                                    -->
                                        (fun x ->
                                            Node.Union.union_final (Some par)
                                            x my_node)
                                and children = 
                                    Codes.union 
                                    child1.children_unions 
                                    child2.children_unions
                                in
                                create res children vertex 
                                (child1.size + child2.size + 1) 
                            in
                            (IntMap.add vertex res acc, clusters), res
                | _ ->
                        (* Create the leaf union and add it to the set of 
                        * vertices *)
                        let res = 
                            my_node 
                            --> Node.Union.leaf None None 
                            --> leaf [vertex] vertex 
                        in
                        (IntMap.add vertex res acc, clusters), res

        (* A function that can process a single handle at a time, and generate
        * the unions tree and the clusters based on that tree *)
        let process_handle saturation_parameter ptree (handle : int) 
            (acc, clusters) =
            match Ptree.get_node handle ptree with
            | Tree.Interior (a, b, _, _) 
            | Tree.Leaf (a, b) -> 
                    let (allacc, child1) = 
                        process_subtree saturation_parameter ptree 
                        (acc, clusters) a 
                    in
                    let (acc, clusters), child2 = 
                        process_subtree saturation_parameter ptree allacc b 
                    in
                    acc, handle_new_cluster a b child1 child2 clusters
            | Tree.Single _ -> 
                    let (acc, cluster), child1 = 
                        process_subtree saturation_parameter ptree 
                        (acc, clusters) handle 
                    in
                    acc, add_to_cluster handle child1 cluster

        let internal_create_initial_tree handle acc saturation_parameter ptree = 
            (* We fold over every handle in the tree and then put the new tree
            * together with the unions of the vertices and the set of 
            * clusters *)
            let new_tree_data, clusters = 
                match handle with
                | None ->
                        Codes.fold 
                        (process_handle saturation_parameter ptree) 
                        (Ptree.get_handles ptree)
                        acc
                | Some v ->
                        process_handle saturation_parameter ptree v acc
            in
            { 
                Ptree.empty with 
                tree = ptree.Ptree.tree;
                Ptree.node_data = new_tree_data 
            },
            clusters

        (* Calculate a tree with Node.Union.n data in the vertices, and a map 
        * of vertices and clusters assigned to groups of vertices in the tree 
        * *)
        let create_initial_tree handle saturation_parameter ptree =
            internal_create_initial_tree handle (IntMap.empty, IntMap.empty) 
            saturation_parameter ptree


        let distance a b = Node.Union.distance a.union b.union

        let rec get_union_vertex a b c =
            if IntMap.mem a c then a
            else if Tree.is_handle a b.Ptree.tree then raise Not_found
            else get_union_vertex (Tree.get_parent a b.Ptree.tree) b c

        module Clusters = struct

            type c = node IntMap.t

            let distance map a b =
                let find x = 
                    try IntMap.find x map with 
                    | Not_found -> 
                            failwith ("UnionTree.Make.Clusters.distance:"
                                ^ "Not_found " ^ string_of_int x)
                in
                distance (find a) (find b)

            let update_matrix matrix clusters added removed =
                let matrix =
                    Codes.fold
                    (fun x acc -> SparceMatrix.remove x acc)
                    removed
                    matrix
                in
                let distance = distance clusters in
                Codes.fold
                (fun x acc -> SparceMatrix.add x distance acc)
                added 
                matrix

            let initial_matrix map = 
                let distance = distance map in
                let add_item_to_matrix x _ acc =
                    SparceMatrix.add x distance acc
                in
                let res = 
                    IntMap.fold add_item_to_matrix map SparceMatrix.empty
                in
                res

            module Test = struct
                let unity matrix cls added removed =
                    let items = SparceMatrix.rows matrix in
                    let is_ok x acc =
                        acc && 
                        ((All_sets.Integers.mem x removed) || 
                        (All_sets.IntegerMap.mem x cls))
                    in
                    All_sets.Integers.fold is_ok 
                    (All_sets.Integers.union items added) true
            end
        end

        (* Update the unions of all the vertices down to the root, starting 
        * in the current vertex. *)
        let update handle saturation_parameter ptree union_tree clusters 
            vertices =
            let removed, union_tree, clusters = 
                (* The clusters that have been removed *)
                Codes.fold (fun x (removed, union_tree, clusters) ->
                    let is_cluster = IntMap.mem x clusters in
                    (if is_cluster then Codes.add x removed
                    else removed),
                    (Ptree.remove_node_data x union_tree),
                    (if is_cluster then IntMap.remove x clusters
                    else clusters))
                vertices
                (Codes.empty, union_tree, clusters)
            in
            let union_tree, new_clusters = 
                let data = union_tree.Ptree.node_data in
                internal_create_initial_tree handle (data, clusters)
                saturation_parameter ptree 
            in
            let added = 
                let compare_clusters old newc = 
                    IntMap.fold (fun x _ acc ->
                        if IntMap.mem x old then acc
                        else Codes.add x acc)
                    newc 
                    Codes.empty
                in
                compare_clusters clusters new_clusters 
            in
            union_tree, new_clusters, added, removed

end

module Make2 (Node : NodeSig.S) (Edge : Edge.EdgeSig with type n = Node.n) : 
    S with type n = Node.n with type u = Node.Union.u 
    with type node = Node.Union.u with type e = Edge.e = struct 

        type n = Node.n
        type u = Node.Union.u
        type e = Edge.e

        type node = u

        let leaf _ _ a = a

        let compare_trees_data a b =
            let comparator x y acc =
                let bitem = IntMap.find x b.Ptree.node_data in
                match Node.Union.compare y bitem with
                | 0 -> acc
                | _ -> x :: acc
            in
            IntMap.fold comparator a.Ptree.node_data []

        (* A function to add a union assigned to a single vertex to a cluster;
        * the vertex is the root of the cluster  *)
        let add_to_cluster = IntMap.add 

        (* Create a new pair of clusters, given that the saturation of the 
        * current cluster has reached more than the acceptable value *)
        let handle_new_cluster = IntMap.add 

        let rec process_subtree ptree vertex parent clusters =
            let my_node = Ptree.get_node_data vertex ptree in
            let is_limit =  Codes.mem vertex clusters in
            match Ptree.get_node vertex ptree with
            | (Tree.Interior _) as node when not is_limit ->
                    let child1, child2 = 
                        let a, b = Tree.other_two_nbrs parent node in
                        process_subtree ptree a vertex clusters,
                        process_subtree ptree b vertex clusters 
                    in
                    Node.Union.union (Some parent) my_node child1 child2
            | _ ->
                    (* Create the leaf union and add it to the set of 
                    * vertices *)
                    let code = Node.min_child_code (Some parent) my_node in
                    Node.Union.leaf (Some code) (Some parent) my_node

        let process_edge code (Tree.Edge (a, b)) ptree clusters =
            let e = 
                (Ptree.get_edge_data (Tree.Edge (a, b)) ptree)
                --> Edge.to_node code (a, b)
            and ua = process_subtree ptree a b clusters
            and ub = process_subtree ptree b a clusters in
            Node.Union.union None e ua ub

        let build_all_clusters clusters ptree =
            (* A function that will fold over all the clusters and will add the
            * calculated cluster into the clusters map. *)
            let process_cluster clusters ptree vertex unions = 
                let union = 
                    match Ptree.get_node vertex ptree with
                    | Tree.Interior (_, _, a, _)
                    | Tree.Leaf (_, a) ->
                            let edge = 
                                Tree.normalize_edge (Tree.Edge (vertex, a))
                                ptree.Ptree.tree
                            in
                            process_edge vertex edge ptree clusters
                    | Tree.Single _ ->
                            let node = Ptree.get_node_data vertex ptree in
                            Node.Union.leaf None None node
                in
                IntMap.add vertex union unions
            in
            Codes.fold (process_cluster clusters ptree) clusters IntMap.empty 

        let rec define_clusters ptree vertex (assignments, clusters) =
            match Ptree.get_node vertex ptree with
            | Tree.Interior (_, par, a, b) ->
                    let u, v, w, x = 
                        define_clusters ptree a (assignments, clusters)
                    in
                    let u, v, y, z = define_clusters ptree b (u, v) in
                    if w + y < union_size then 
                        let assigned = Codes.union x z in
                        assignments, clusters, w + y + 1, 
                        Codes.add vertex assigned
                    else
                        (* We define each child as a root of a cluster *)
                        let assignments = IntMap.add b z (IntMap.add a x u) 
                        and clusters = Codes.add a (Codes.add b clusters) in
                        assignments, clusters, 1, Codes.singleton vertex
            | Tree.Leaf (a, b) ->
                    assignments, clusters, 1, Codes.singleton vertex
            | Tree.Single _ -> 
                    assignments, Codes.add vertex clusters, 1,
                    Codes.singleton vertex

        let define_clusters_handle ptree handle ((assgn, clst) as acc) = 
            match Ptree.get_node handle ptree with
            | Tree.Interior (a, b, _, _) 
            | Tree.Leaf (a, b) ->
                    let u, v, w, x = define_clusters ptree a acc in
                    let u, v, y, z = define_clusters ptree b (u, v) in
                    let assignments = IntMap.add b z (IntMap.add a x u) 
                    and clusters = Codes.add a (Codes.add b v) in
                    assignments, clusters
            | Tree.Single _ ->
                    IntMap.add handle (Codes.singleton handle) assgn,
                    Codes.add handle clst

        (* Calculate a tree with Node.Union.n data in the vertices, and a map 
        * of vertices and clusters assigned to groups of vertices in the tree 
        * *)
        let create_initial_tree _ _ ptree =
            let _, clst = 
                Codes.fold (define_clusters_handle ptree)  
                ptree.Ptree.tree.Tree.handles (IntMap.empty, Codes.empty)
            in
            let unions = build_all_clusters clst ptree in
            { 
                Ptree.empty with 
                Ptree.node_data = unions;
                tree = ptree.Ptree.tree;
            }, unions


        let distance a b = Node.Union.distance a b

        let rec get_union_vertex a b c =
            if IntMap.mem a c then a
            else if Tree.is_handle a b.Ptree.tree then raise Not_found
            else get_union_vertex (Tree.get_parent a b.Ptree.tree) b c

        module Clusters = struct

            type c = node IntMap.t

            let distance map a b =
                let find x = 
                    try IntMap.find x map with 
                    | Not_found -> 
                            failwith ("UnionTree.Make.Clusters.distance:"
                                ^ "Not_found " ^ string_of_int x)
                in
                distance (find a) (find b)

            let update_matrix matrix clusters added removed =
                let matrix =
                    Codes.fold
                    (fun x acc -> SparceMatrix.remove x acc)
                    removed
                    matrix
                in
                let distance = distance clusters in
                Codes.fold
                (fun x acc -> SparceMatrix.add x distance acc)
                added 
                matrix

            let initial_matrix map = 
                let distance = distance map in
                let add_item_to_matrix x _ acc =
                    SparceMatrix.add x distance acc
                in
                let res = 
                    IntMap.fold add_item_to_matrix map SparceMatrix.empty
                in
                res

            module Test = struct
                let unity matrix cls added removed =
                    let items = SparceMatrix.rows matrix in
                    let is_ok x acc =
                        acc && 
                        ((All_sets.Integers.mem x removed) || 
                        (All_sets.IntegerMap.mem x cls))
                    in
                    All_sets.Integers.fold is_ok 
                    (All_sets.Integers.union items added) true
            end
        end

        (* Update the unions of all the vertices down to the root, starting 
        * in the current vertex. *)
        let update _ _ ptree _ prev _ =
            let union_tree, clusters = create_initial_tree None 0.0 ptree in
            let get_elements map = 
                IntMap.fold (fun x _ acc -> Codes.add x acc) map
                Codes.empty 
            in
            union_tree, clusters, get_elements clusters, get_elements prev
end

