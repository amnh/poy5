(* POY 5.1.1. A phylogenetic analysis program using Dynamic Homologies.       *)
(* Copyright (C) 2014 Andrés Varón, Lin Hong, Nicholas Lucaroni, Ward Wheeler,*)
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

let () = SadmanOutput.register "Build" "$Revision: 3649 $"

let debug_profile_memory = false

let (-->) a b = b a

let current_snapshot x = 
    if debug_profile_memory then 
        let () = Printf.printf "%s\n%!" x in
        MemProfiler.current_snapshot x
    else ()

let rec build_features (meth:Methods.build) =
    match meth with
    | `Prebuilt fn ->
            let fn = 
                match fn with
                | `Local x | `Remote x -> x
            in
            ("type", "prebuilt") :: 
            ("filename", fn) :: []
    | `Nj -> ("type", "neighbor joining") :: []
    | `Build (n, meth, _, _) -> (* TODO: report type of iteration patterns; here, BB and RT *)
            let other_features str max_n = 
                [ ("type", str); 
                ("number of trees to keep", string_of_int max_n)]
            in
            ("number of trees", string_of_int n) ::
                (match meth with 
                    | `Wagner_Rnd (max_n, _, _, _, _) -> 
                            other_features "randomized wagner" max_n
                    | `Wagner_Ordered (max_n, _, _, _, _) ->
                            other_features "ordered wagner" max_n
                    | `Wagner_Mst (max_n, _, _, _, _) ->
                            other_features "minimum spanning tree wagner" max_n
                    | `Wagner_Distances (max_n, _,  _, _, _) ->
                            other_features "distances based tree wagner" max_n
                    | `Nj
                    | (`Branch_and_Bound _) 
                    | (`Prebuilt _) as x ->
                            build_features x
                    | _ -> [] )
    | `Branch_and_Bound ((bound, threshold, keep_method, max_trees, _),_) ->
            [("type", "Branch and bound"); ("initial_bound", match bound with | None ->
                "none" | Some x -> string_of_float x); 
                ("threshold", match threshold with
                | None -> "0.0" | Some x -> string_of_float x);
                ("number of trees to keep", string_of_int max_trees);
                ("keep method", match keep_method with
                | `Last -> "last"
                | `First -> "first"
                | `Keep_Random -> "random")]
    | `Build_Random _ ->  
            [("type", "Uniformly at random")]


let remove_exact (meth : Methods.cost_calculation) (acc : Methods.transform
list) : Methods.transform list =
    match meth with
    | #Methods.transform as meth -> meth :: acc

let rec get_transformations (meth : Methods.build) : Methods.transform list =
    match meth with
    | `Nj
    | `Prebuilt _ -> []
    | `Build (_, build_meth, trans,_) ->
            List.fold_right remove_exact (trans @ 
            (match build_meth with
            | `Branch_and_Bound ((_, _, _, _, trans),_)
            | `Wagner_Distances (_, _, _, trans, _)
            | `Wagner_Mst (_, _, _, trans, _) 
            | `Wagner_Rnd (_, _, _, trans, _) 
            | `Wagner_Ordered (_, _, _, trans, _)
            | `Constraint (_, _, _, trans)
            | `Build_Random ((_, _, _, trans, _),_) -> trans
            | `Nj
            | `Prebuilt _ -> [])) []
    | `Branch_and_Bound ((_, _, _, _, trans),_)
    | `Build_Random ((_, _, _, trans, _),_) -> 
            List.fold_right remove_exact trans []

module type S = sig
      type a 
      type b

      val report_mst : Data.d -> a list -> string option -> unit

        (** [prebuilt a b] generates a list of trees using the list of trees [a] and the
    * node data [b]. It is required that the names of the taxa in [a] have data
    * associated in [b]. [b] is generated by the {!Node.load_data} function. *)
      val prebuilt :
        (string option * Tree.Parse.tree_types list) list ->
        Data.d * a list ->
        [> `Set of [> `Single of (a, b) Ptree.p_tree ] list ]

    val build_initial_trees :
        (a, b) Ptree.p_tree Sexpr.t -> Data.d ->
      a list -> Methods.build -> (a, b) Ptree.p_tree Sexpr.t

    end

module MakeNormal (Node : NodeSig.S) (Edge : Edge.EdgeSig with type n = Node.n) 
    (TreeOps : Ptree.Tree_Operations with type a = Node.n with type b = Edge.e)
    = struct

    module PtreeSearch = Ptree.Search (Node) (Edge) (TreeOps)
    module BuildTabus = Tabus.Make (Node) (Edge)
    module Queues = Queues.Make (Node) (Edge)

    type a = PtreeSearch.a
    type b = PtreeSearch.b
    type phylogeny = (a, b) Ptree.p_tree

    let map_of_list f x =
        List.fold_left (fun acc x -> All_sets.IntegerMap.add (f x) x acc)
        All_sets.IntegerMap.empty x


    let randomize lst =
        let arr = Array.of_list lst in
        Array_ops.randomize arr;
        Array.to_list arr


    let single_wagner_search_manager = new Queues.wagner_srch_mgr true 1 0.0


    let set_of_leafs data = 
        map_of_list (fun x -> Node.taxon_code x) data 


    let disjoin_tree data node =
        let leafs = set_of_leafs node in
        let tree = Ptree.make_disjoint_tree data leafs in
        tree --> PtreeSearch.downpass --> PtreeSearch.uppass


    let edges_of_tree tree =
        Tree.EdgeSet.fold (fun x acc -> x :: acc) tree.Ptree.tree.Tree.d_edges []


    let random_tree (data : Data.d) (nodes : a list) adj_mgr =
        let tree : phylogeny = 
            { (Ptree.empty data) with
                Ptree.tree = Tree.random (List.map Node.taxon_code nodes);
                Ptree.node_data = map_of_list Node.taxon_code nodes; }
        in
        tree --> PtreeSearch.downpass --> PtreeSearch.uppass


    let branch_and_bound keep_method max_trees threshold data nodes bound adj_mgr =
        let select_appropriate bound_plus_threshold lst =
            let lst = List.filter (fun x -> bound_plus_threshold >=
                Ptree.get_cost `Adjusted x) lst in 
            let len = List.length lst in
            if len <= max_trees then lst
            else
                let () = assert (len = max_trees + 1) in
                match keep_method with
                | `Last ->
                        (match List.rev lst with
                        | h :: lst -> List.rev lst
                        | [] -> assert false)
                | `First ->
                        (match lst with
                        | h :: lst -> lst
                        | [] -> assert false)
                | `Keep_Random ->
                        (let arr = Array.of_list lst in
                        Array_ops.randomize arr;
                        match Array.to_list arr with
                        | h :: lst -> lst
                        | [] -> assert false)
        in
        (* calculate the initial bound if nto provided; /2 for threshold *)
        let bound = match bound with
            | None -> max_float /. 2.
            | Some x -> x
        (* create status for ncurses *)
        and st = Status.create "Branch and Bound Build" (Some 100) "percent complete" in
        let () = Status.full_report ~adv:(0) st in
        (* We need to present some output; find a decent depth and that percentage done *)
        let report_depth,report_percent =
            (* the number of possibilities at level n *)
            let n t =
                let rec n acc t = match t with
                    | 0 | 1 | 2 | 3 -> acc
                    | t             -> n (acc*(2*t-5)) (t-1)
                in
                n 1 t
            (* depth=6 will prune ~1% of the tree (having 105 possibilities) *)
            and depth = max 1 (min ((List.length nodes)-1) 6) and p = ref 0.0 in
            depth,
            (fun depth ->
                let p_incr = (1.0 /. float (n depth)) *. 100.0 in
                p := p_incr +. !p;
                Status.full_report ~adv:(int_of_float !p) st)
        in
        let rec aux_branch_and_bound depth ((bound, best_trees) as acc) tree
                    edges cur_handle other_handles =
            match edges with
            | (Tree.Edge (x, y)) :: t ->
                if report_depth = depth then report_percent depth;
                let new_tree, _ =
                    TreeOps.join_fn adj_mgr [] (Tree.Edge_Jxn (x,y)) cur_handle tree
                in
                let new_cost = Ptree.get_cost `Adjusted new_tree in
                if new_cost >  bound +. threshold then begin
                    if depth < report_depth then report_percent depth;
                    aux_branch_and_bound depth acc tree t cur_handle other_handles
                end else begin
                    let acc = match other_handles with
                        | nh :: oh ->
                            aux_branch_and_bound (depth+1) acc new_tree
                                (edges_of_tree new_tree) (Tree.Single_Jxn nh) oh
                        | [] -> 
                            let realbound = bound +. threshold in
                            if new_cost < bound then
                                new_cost, select_appropriate (new_cost +. threshold) (new_tree :: best_trees)
                            else if new_cost <= realbound then
                                (bound, select_appropriate realbound (new_tree::best_trees))
                            else
                                acc
                    in
                    aux_branch_and_bound depth acc tree t cur_handle other_handles
                end
            | [] -> acc
        in
        let initial_tree = disjoin_tree data nodes in
        let _, trees =
            if max_trees < 1 then 0., []
            else begin match List.map Node.taxon_code nodes with
                | f :: s :: tail ->
                    let new_tree, _ =
                        TreeOps.join_fn adj_mgr [] (Tree.Single_Jxn f) (Tree.Single_Jxn s) initial_tree
                    in
                    begin match tail with
                        | t :: tail ->
                            let edges = edges_of_tree new_tree in
                            aux_branch_and_bound 3 (bound,[]) new_tree edges (Tree.Single_Jxn t) tail
                        | [] ->
                            Ptree.get_cost `Adjusted new_tree, [new_tree]
                    end
                | _ -> 0., [initial_tree]
            end
        in
        let () = Status.finished st in
        Sexpr.of_list (List.map PtreeSearch.uppass trees)

    let sort_list_of_trees ptrees = 
        let cost = Ptree.get_cost `Adjusted in
        List.sort (fun a b -> compare (cost a) (cost b)) ptrees

    let constrained_build cg data n constraint_tree nodes adj_mgr = 
        let rec randomize_tree tree = 
            match tree with
            | Tree.Parse.Leafp _ -> tree
            | Tree.Parse.Nodep (lst, res) ->
                    let arr = Array.map randomize_tree (Array.of_list lst) in
                    Array_ops.randomize arr;
                    Tree.Parse.Nodep ((Array.to_list arr), res)
        in
        let ptree = disjoin_tree data nodes in
        let rec aux_constructor fn ptree tree = 
            match tree with
            | Tree.Parse.Leafp x -> 
                    ptree, (Data.taxon_code (fn x) data)
            | Tree.Parse.Nodep ([x], _) ->
                    aux_constructor fn ptree x
            | Tree.Parse.Nodep (lst, _) ->
                    let ptree, handles = 
                        List.fold_left (fun (ptree, acc) item ->
                            let ptree, nh = aux_constructor fn ptree item in
                            ptree, (nh :: acc)) (ptree, []) lst 
                    in
                    let handles = 
                        List.map 
                        (fun (x : int) -> Ptree.handle_of x ptree) 
                        handles
                    in
                    let constraints = List.fold_left (fun acc x ->
                        let acc = All_sets.Integers.add x acc in
                        try 
                            let parent = Ptree.get_parent x ptree in
                            All_sets.Integers.add parent acc
                        with
                        | _ -> acc) All_sets.Integers.empty handles 
                    in
                    let ptrees = 
                        PtreeSearch.make_wagner_tree ~sequence:handles
                        ptree adj_mgr single_wagner_search_manager 
                        (BuildTabus.wagner_constraint
                        constraints)
                    in
                    let ptrees = sort_list_of_trees ptrees in
                    match ptrees, handles with
                    | (ptree :: _), (h :: _) -> 
                             ptree, h
                    | _ -> assert false
        in
        let st = Status.create "Constrained Wagner Replicate" (Some n)
        "Constrained Wagner tree replicate building" in
        let rec total_builder res blt =
            Status.full_report ~adv:(n - blt) st;
            if blt = n then 
                let () = Status.finished st in
                `Set res
            else 
                let rec deal_with_tree = function
                    | Tree.Parse.Annotated (t,_) ->
                        deal_with_tree t
                    | Tree.Parse.Flat t -> 
                        let tree,_= aux_constructor (fun x-> x) ptree (randomize_tree t) in
                        tree
                    | Tree.Parse.Branches t ->
                        let tree,_= aux_constructor (fst) ptree (randomize_tree t) in
                        tree
                    | Tree.Parse.Characters t ->
                        let tree,_ = aux_constructor (fst) ptree (randomize_tree t) in
                        tree
                in
                let tree = deal_with_tree constraint_tree in
                let tree = TreeOps.uppass tree in
                total_builder ((`Single tree) :: res) (blt + 1)
        in
        if n < 0 then `Empty 
        else
            total_builder [] 0

    let single_wagner data tabu_mgr wmgr cg nodes adj_mgr = 
        let disjoin_tree = disjoin_tree data nodes
        and nodes = List.map (fun x -> Node.taxon_code x) nodes in
        let ptrees = 
            PtreeSearch.make_wagner_tree ~sequence:nodes 
                                         disjoin_tree adj_mgr wmgr tabu_mgr
        in
        let ptrees = sort_list_of_trees ptrees in
        match ptrees with
        | hd :: _ -> `Single (PtreeSearch.uppass hd)
        | _ -> failwith "No wagner trees built!"

    let wagner data tabu_mgr cg nodes adj_mgr =
        single_wagner data tabu_mgr single_wagner_search_manager cg nodes adj_mgr

    let randomized_single_wagner data tabu_mgr randomize wmgr cg adj_mgr = 
        let nodes = randomize () in
        single_wagner data tabu_mgr wmgr cg nodes adj_mgr

    (** [rand_wagner a b] creates a fresh wagner tree with its 
    * corresponding cost using
    * the initial set of leaves as generated by [b ()] and the code 
    * generator [a]. The addition
    * sequence depends on the list as produced by [b] 
    * before creating the tree. *)
    let rand_wagner data tabu_mgr data_generator cg =
        randomized_single_wagner data tabu_mgr data_generator 
        single_wagner_search_manager cg 

    let n_independent_wagner data_generator tabu_mgr cg data n =    
        (* This is equivalent to having multiple replicates *)
        let st = Status.create "Wagner Replicate" (Some n) "Wagner tree replicate \
        building" in
        let mgr = single_wagner_search_manager in
        let rec builder cnt acc =
            if cnt > 0 then begin
                Status.full_report ~adv:(n - cnt) st;
                let next = 
                    (randomized_single_wagner data tabu_mgr data_generator mgr cg) :: acc 
                in
                builder (cnt - 1) next
            end else begin
                Status.finished st;
                `Set acc
            end
        in
        builder n []

    let create_adjust_manager (m,b) = 
        let thrsh = match m with 
            | `Threshold f 
            | `Both (f,_) -> Some f
            | `Always     -> Some 0.0
            | `Null 
            | `MaxCount _ -> None
        and count =  match m with
            | `MaxCount m 
            | `Both (_,m) -> Some m
            | `Always     -> Some 0
            | `Null 
            | `Threshold _ -> None
        in
        let mgr = match b with
            | `Null           -> BuildTabus.simple_nm_none count thrsh
            | `AllBranches    -> BuildTabus.simple_nm_all count thrsh
            | `JoinDelta      -> BuildTabus.complex_nm_delta count thrsh
            | `Neighborhood x -> BuildTabus.complex_nm_neighborhood x count thrsh
        in
        Some mgr

    (* compose the iteration manager tabu *)
    let pick_tabu_manager = function
        | `UnionBased _ -> (* Maximum distance is useless in this case *)
                BuildTabus.wagner_tabu
        | `AllBased _
        | `Partition _ ->
                BuildTabus.distance_dfs_wagner

        (** [max_n_wagner a b n] creates a list of at most [n] wagner trees using those
    * other trees that have the same cost of the current best tree. It uses the
    * original addition sequence as generated by [b ()].*)
    let max_n_wagner data threshold tabu_mgr data_generator cg n adj_mgr =
        let wmgr = new Queues.wagner_srch_mgr true n threshold in
        let nodes = data_generator () in
        let nodesl = List.map (fun x -> Node.taxon_code x) nodes in
        let res = 
            PtreeSearch.make_wagner_tree ~sequence:nodesl 
                (disjoin_tree data nodes) adj_mgr wmgr tabu_mgr 
        in
        `Set (List.map (fun x -> `Single (TreeOps.uppass x)) res)

    let make_distances_table ?(both=true) nodes = 
        let tmp = Hashtbl.create 99991 in
        List.iter (fun a ->
            List.iter (fun b ->
                let ca = Node.taxon_code a
                and cb = Node.taxon_code b in
                if ca = cb then ()
                else
                    if Hashtbl.mem tmp (ca, cb) then ()
                    else begin
                        let d = Node.distance 100000. a b in
                        Hashtbl.add tmp (ca, cb) d;
                        if both then Hashtbl.add tmp (cb, ca) d;
            end) nodes) nodes;
            tmp

    module OrderedPairs = struct
        type t = (float * (int * int))
        let compare a b = compare ((fst a) : float) (fst b)
    end

    module H = Heap.Make (OrderedPairs)

    let table_of_trees nodes = 
        let tbl = Hashtbl.create 1667 in
        List.iter (fun n -> Hashtbl.add tbl (Node.taxon_code n)
        (Tree.Parse.Leafp (Node.taxon_code n))) nodes;
        tbl

    let set_of_trees nodes =
        List.fold_left (fun acc x ->
            All_sets.Integers.add (Node.taxon_code x) acc) 
        All_sets.Integers.empty nodes

    let (-->) a b = b a 

    let nj_qtable table terminals =
        let d x y = Hashtbl.find table ((min x y), (max x y)) in
        let r = float_of_int (All_sets.Integers.cardinal terminals) in
        let q_table = Hashtbl.create 99991 in
        Hashtbl.iter (fun ((i, j) as p) dist ->
            let sum =
                All_sets.Integers.fold (fun c sum ->
                    sum -.
                    (if c = i then 0. else d c i) -.
                    (if c = j then 0.  else d c j)) 
                terminals 0.
            in
            Hashtbl.add q_table p (((r -. 2.) *. dist) +. sum)) table;
        q_table

    let nj_distance_to_ancestor table terminals f g =
        let d x y = Hashtbl.find table ((min x y), (max x y)) in
        let sum_of_distance_to x =
            All_sets.Integers.fold (fun y sum -> sum +. (d x y)) terminals 0.
        in
        let r = float_of_int (All_sets.Integers.cardinal terminals) in
        ((0.5 *. (d f g)) 
        +.  ((1. /. (2. *. (r -. 2.))) *.
            ((sum_of_distance_to f) -. (sum_of_distance_to g))))

    let nj_new_distance table distance_fu distance_gu f g k =
        let d x y = Hashtbl.find table ((min x y), (max x y)) in
        ((0.5 *. ((d f k) -.  distance_fu)) +.
        (0.5 *. ((d g k) -.  distance_gu)))

    let join_trees code distance_table a b tree_table trees heap =
        let trees = 
            trees
            --> All_sets.Integers.remove a 
            --> All_sets.Integers.remove b 
        in
        let ta = Hashtbl.find tree_table a
        and tb = Hashtbl.find tree_table b in
        let tab = Tree.Parse.Nodep ([ta; tb], code) in
        Hashtbl.add tree_table code tab;
        let heap = 
            let distance_acode = 
                nj_distance_to_ancestor distance_table trees a b
            and distance_bcode = 
                nj_distance_to_ancestor distance_table trees b a
            in
            All_sets.Integers.fold (fun c heap ->
                let dc = 
                    nj_new_distance distance_table distance_acode
                    distance_bcode a b c
                in
                let pair = (code, c) in
                Hashtbl.add distance_table pair dc;
                H.insert (dc, pair) heap) trees heap
        in
        let trees = All_sets.Integers.add code trees in
        trees, code - 1, heap


    let rec merge_nj_trees code distance_table tree_table trees heap =
        let (_, (a, b)) = H.findMin heap in
        let heap = H.deleteMin heap in
        if All_sets.Integers.mem a trees && 
            All_sets.Integers.mem b trees then
                join_trees code distance_table a b tree_table trees heap
        else 
            merge_nj_trees code distance_table tree_table trees heap

    let nj data nodes = 
        let distance_table = 
            let both = false in
            make_distances_table ~both nodes 
        in
        let heap =
            Hashtbl.fold 
            (fun x y acc ->
                H.insert (y, x) acc) 
            distance_table H.empty
        in
        let tree_table = table_of_trees nodes in
        let trees = set_of_trees nodes in
        let rec complete_merge code trees heap =
            if 1 = All_sets.Integers.cardinal trees then
                Hashtbl.find tree_table (All_sets.Integers.choose trees)
            else 
                let trees, code, heap = 
                    merge_nj_trees code distance_table tree_table trees heap
                in
                complete_merge code trees heap
        in
        let tree = complete_merge (-1) trees heap in
        Tree.Parse.map 
                (fun x -> 
                    if x >= 0 then Data.code_taxon x data
                        else "") tree


    let distances_ordered nodes = 
        let distances_table = make_distances_table nodes in
        let distances_list = 
            let tmp =
                Hashtbl.fold (fun a b acc ->
                    (a, b) :: acc) distances_table [] 
            in
            List.sort (fun (_, a) (_, b) -> compare a b) tmp
        in
        let _, addition_list = 
            let add_one visited acc x =
                if All_sets.Integers.mem x visited then visited, acc
                else All_sets.Integers.add x visited, x :: acc
            in
            List.fold_left (fun (visited, acc) ((x, y), _) ->
                let v, a = add_one visited acc x in
                add_one v a y) (All_sets.Integers.empty, []) distances_list
        in
        let addition_list = List.rev addition_list in
        let rec addition_function lst =
            match lst with
            | h1 :: ((h2 :: t) as rest) ->
                    if 0 = Random.int 2 then 
                        h1 :: (addition_function rest)
                    else h2 :: (addition_function (h1 :: t))
            | _ -> lst
        in
        let create_list () =
            let lst = addition_function addition_list in
            List.rev (List.map (fun x ->
                List.find (fun y -> 
                    x = Node.taxon_code y) nodes) lst)
        in
        create_list

    let mst data nodes =
        let distances_table = make_distances_table nodes in
        let distance_fn a b = Hashtbl.find distances_table (a, b) 
        and codes = List.map Node.taxon_code nodes in
        let mst = Mst.kruskal Mst.Closest distance_fn codes in
        let do_mst () =
            let data = Mst.bfs_traversal Mst.Closest2 mst in
            List.map (fun x ->
                List.find (fun y -> 
                    x = Node.taxon_code y) nodes) data 
        in
        do_mst

    let report_mst data nodes filename =
        let distances_table = make_distances_table nodes in
        let distance_fn a b = Hashtbl.find distances_table (a, b) 
        and codes = List.map Node.taxon_code nodes in
        let mst = Mst.kruskal Mst.Closest distance_fn codes in
        Mst.print_mst_tree (fun x -> Data.code_taxon x data) mst filename

    let max_n_dg_p_wagner data threshold tabu_mgr data_generator cg n p adj_mgr: phylogeny Sexpr.t =
        match p, n with
        | 0, _ -> `Empty
        | 1, 1 -> 
                let st = Status.create "Building Wagner Tree" None "" in
                let res = rand_wagner data tabu_mgr data_generator cg adj_mgr in
                Status.finished st;
                res
        | 1, n -> 
                let st = Status.create "Building Wagner Tree" None "" in
                let res = max_n_wagner data threshold tabu_mgr data_generator cg n adj_mgr in
                Status.finished st;
                res
        | p, _ ->
              let builder cnt acc =
                  let next = max_n_wagner data threshold tabu_mgr data_generator cg n adj_mgr in
                  next :: acc
              in
              `Set (Sexpr.compose_status "Wagner build" builder p [])

    (** [max_n_randomized_p_wagner a b n p] generates [p] independent
    * wagner trees, on each keeping the best [n] trees found on each step,
    * following the addition sequence as specified by the node generating
    * function [b ()]. *)
    let max_n_randomized_p_wagner data threshold tabu_mgr cg nodes n p adj =
        let data_generator () = randomize nodes in
        max_n_dg_p_wagner data threshold tabu_mgr data_generator cg n p adj

    let max_n_mst_p_wagner data threshold tabu_mgr node_data cg nodes n p adj =
        let mst = mst node_data nodes in
        max_n_dg_p_wagner data threshold tabu_mgr mst cg n p adj

    let max_n_distances_p_wagner data threshold tabu_mgr cg nodes n p adj =
        let dord = distances_ordered nodes in
        max_n_dg_p_wagner data threshold tabu_mgr dord cg n p adj

    let split_in_forests trees = 
        let make_tree_set_of_taxa acc x = 
            let rec make_tree_set_of_taxa acc x = 
                match x with
                | Tree.Parse.Leafp name -> All_sets.Strings.add name acc
                | Tree.Parse.Nodep (chld, _) ->
                    List.fold_left make_tree_set_of_taxa acc chld
            in
            make_tree_set_of_taxa acc (Tree.Parse.strip_tree x)
        in
        let are_different a acc b =
            acc && (All_sets.Strings.is_empty (All_sets.Strings.inter a b))
        in
        let are_same a acc b =
            acc && (0 = (All_sets.Strings.compare a b))
        in
        let rec are_something pairwise_comparison acc lst = 
            match lst with
            | h :: t ->
                    let acc = List.fold_left (pairwise_comparison h) acc t in
                    are_something pairwise_comparison acc t
            | [] -> acc
        in
        let are_all_the_same_set lst =
            are_something are_same true lst
        and are_all_different lst =
            are_something are_different true lst
        in
        List.fold_left (fun acc (name,x) ->
            if 1 = List.length x then (name,x) :: acc
            else 
                let taxa = 
                    List.map (make_tree_set_of_taxa All_sets.Strings.empty) x 
                in
                if are_all_the_same_set taxa then 
                    let x = List.map (fun x -> name,[x]) x in
                    x @ acc
                else if are_all_different taxa then (name,x) :: acc
                else 
                    let _ =
                        Status.user_message Status.Error 
                        ("While@ trying@ to@ read@ the@ trees@ from@ the@ "
                        ^ "input@ files@ I@ have@ found@ some@ internal@ "
                        ^ "inconsistencies:@ POY@ can@ read@ either@ forests@ "
                        ^ "or@ trees,@ and@ recognize@ each@ by@ comparing@ "
                        ^ "the@ trees@ in@ memory,@ either@ all@ the@ trees@ "
                        ^ "between@ separators@ (, or ;)@ share@ same@ taxa@ "
                        ^ "in@ which@ case@ I@ treat@ them@ as@ just@ trees@ "
                        ^ "or@ they@ are@ disjoint@, and@ I@ treat@ them@ as@ "
                        ^ "a@ forest.@ Your@ input@ don't@ have@ one@ of@ "
                        ^ "those@ properties.@ I@ think@ you@ intend@ to@ "
                        ^ "read@ just@ trees@, but@ there@ is@ some@ tree@ "
                        ^ "with@ taxa@ that@ doesn't@ appear@ in@ some@ other@ "
                        ^ "tree.@ Sorry,@ I@ can't@ recover@ from@ this,@ and@ "
                        ^ "won't@ load@ the@ trees@ you@ gave@ me.")
                    in
                    failwith "Illegal tree input")
                [] trees

    let prebuilt (trees: (string option * Tree.Parse.tree_types list) list) ((data,_) as sumdata) =
        let trees = split_in_forests trees in
        let st = Status.create "Loading Trees" (Some (List.length trees)) "" in
        let constructor (cnt, lst) x =
            Status.full_report ~adv:cnt st;
            let t = 
                current_snapshot "Build.prebuilt.constructor begin";
                let tree = PtreeSearch.convert_to x sumdata in
                current_snapshot "Build.prebuilt.constructor converted";
                let tree = PtreeSearch.downpass tree in
                current_snapshot "Build.prebuilt.constructor downpass";
                let tree = PtreeSearch.uppass tree in
                current_snapshot "Build.prebuilt.constructor uppass";
                tree
            in
            cnt + 1, (`Single t) :: lst
        in
        let res = 
            let _, res = List.fold_left constructor (1, []) trees in
            `Set res 
        in
        Status.finished st;
        res

    let prebuilt trees sumdata = match trees with
        | [] -> `Set []
        | xs -> prebuilt xs sumdata

    let rec build_initial_trees trees data nodes (meth : Methods.build) =
        let d = (data, nodes) in
        let cg = 
            let code = ref data.Data.number_of_taxa in
            fun () -> incr code; !code
        in
        let built_tree_report acc trees =
            let builder (acc, cnt) t =
                let cost = Ptree.get_cost `Adjusted t in
                let hd = ("tree_" ^ string_of_int cnt ^ "_cost", string_of_float
                    cost) in
                hd :: acc, cnt + 1
            in
            builder acc trees
        in
        let do_constraint file = match file with
            | None -> 
                let hd, tree_list = match Sexpr.to_list trees with
                    | (h :: _) as t -> h, t
                    | [] -> failwith "No trees for constraint"
                in
                let maj = float_of_int (List.length tree_list) in
                Ptree.consensus
                    (PtreeSearch.get_collapse_function None)
                    (fun code -> Data.code_taxon code data)
                    (maj)
                    (Sexpr.to_list trees)
                    (match data.Data.root_at with
                        | Some v -> v
                        | None ->
                            let f = Sexpr.first trees in
                            Ptree.choose_leaf f)
            | Some file ->
                begin match (Data.process_trees data file).Data.trees with
                    | [((_,[t]), _, _) as one] when Data.verify_trees data one -> t
                    | _ -> failwith "Illegal input constraint file"
                end
        in
        let perform_build () = match meth with
            | `Branch_and_Bound ((bound, threshold, keep_method, max_trees, _),adj_meth) ->
                let threshold = match threshold with
                    | None -> 0.
                    | Some x -> x
                and adj_mgr = create_adjust_manager adj_meth in
                branch_and_bound keep_method max_trees threshold data nodes bound adj_mgr
            | `Prebuilt file -> 
                let data = Data.process_trees data file in
                let trees = List.filter (Data.verify_trees data) data.Data.trees in
                let trees = List.map (fun (a, _, id) -> a) trees in
                prebuilt trees d
            | `Nj ->
                let tree = None, [Tree.Parse.Flat (nj data nodes)] in
                prebuilt [tree] d
            | `Build (n, build_meth, lst, adj_meth) ->
                let new_nodes = nodes
                and adj_mgr = create_adjust_manager adj_meth in
                (** TODO: Add different cost calculation heuristic methods *)
                (** TODO: Add different keep methods *)
                if n < 1 then trees
                else
                    begin match build_meth with
                        | `Constraint (_, threshold, file, _) ->
                            let constraint_tree = do_constraint file in
                            constrained_build cg data n constraint_tree nodes adj_mgr
                        | `Branch_and_Bound ((bound, threshold, keep_method, max_trees, _),_) ->
                            let threshold = 
                                match threshold with
                                | None -> 0.
                                | Some x -> x
                            in
                            branch_and_bound keep_method (n * max_trees) threshold data nodes bound adj_mgr
                        | `Wagner_Rnd (max_n, threshold, _, lst, tabu_mgr) ->
                            let tabu_mgr = pick_tabu_manager tabu_mgr in
                            let res = 
                                max_n_randomized_p_wagner data threshold 
                                        tabu_mgr cg new_nodes max_n n adj_mgr
                            in
                            Sexpr.of_list (Sexpr.to_list res)
                        | `Wagner_Ordered (max_n, threshold, keep_method, lst, tabu_mgr) ->
                            let tabu_mgr = pick_tabu_manager tabu_mgr in
                            let status = Status.create "Wagner ordered build" None "" in
                            let () =
                                Status.full_report ~msg:"Building ordered tree" status 
                            in
                            let hd = 
                                max_n_wagner data threshold tabu_mgr 
                                    (fun () -> new_nodes) cg max_n adj_mgr
                            in
                            let () =
                                Status.full_report ~msg:"Building random trees"
                                status in
                            let tl = 
                                max_n_randomized_p_wagner data threshold tabu_mgr
                                            cg new_nodes max_n (n - 1) adj_mgr
                            in
                            let () = Status.finished status in
                            Sexpr.of_list (Sexpr.to_list (`Set [hd; tl]))
                        | `Wagner_Distances (max_n, threshold, keep_method, lst, tabu_mgr) ->
                            let tabu_mgr = pick_tabu_manager tabu_mgr in
                            let status = 
                                Status.create "Wagner Distances-ordered build" 
                                None ""
                            in
                            let () =
                                Status.full_report ~msg:"Building trees"
                                status in
                            let lst = 
                                max_n_distances_p_wagner data threshold tabu_mgr
                                                    cg new_nodes max_n n adj_mgr
                            in
                            let () = Status.finished status in
                            Sexpr.of_list (Sexpr.to_list lst)
                        | `Wagner_Mst (max_n, threshold, keep_method, lst, tabu_mgr) ->
                            let tabu_mgr = pick_tabu_manager tabu_mgr in
                            let status = 
                                Status.create "Wagner MST-ordered build" 
                                None ""
                            in
                            let () =
                                Status.full_report ~msg:"Building trees"
                                status in
                            let lst = 
                                max_n_mst_p_wagner data threshold tabu_mgr data 
                                                    cg new_nodes max_n n adj_mgr
                            in
                            let () = Status.finished status in
                            Sexpr.of_list (Sexpr.to_list lst)
                        | `Nj
                        | (`Prebuilt _) as x -> build_initial_trees trees data nodes x
                        | `Build_Random _ ->
                            let st = Status.create "Random Trees build" (Some n) "" in
                            let arr = 
                                Array.init n 
                                    (fun x ->
                                        Status.full_report ~adv:x st;
                                        random_tree data nodes adj_mgr)
                            in
                            Status.finished st;
                            Sexpr.of_list (Array.to_list arr)
                end;
            | `Build_Random ((n, _, _, _, _),adj_meth) -> 
                let st = Status.create "Random Trees build" (Some n) ""
                and adj_mgr = create_adjust_manager adj_meth in
                let arr = 
                    Array.init n 
                        (fun x ->
                            Status.full_report ~adv:x st;
                            random_tree data nodes adj_mgr)
                in
                Status.finished st;
                Sexpr.of_list (Array.to_list arr)
        in
        Sadman.start "build" (build_features meth);
        let timer = Timer.start () in
        let res = perform_build () in
        let time = Timer.get_user timer in
        let report, n = 
            Sexpr.fold_left (fun acc x -> built_tree_report acc x) ([], 0) res
        in
        Sadman.finish ((TreeSearch.search_time_and_trees_considered time n) @ report);
        res

end

module Make (NodeH : NodeSig.S with type other_n = Node.Standard.n)
            (EdgeH : Edge.EdgeSig with type n = NodeH.n) 
            (TreeOps : Ptree.Tree_Operations with type a = NodeH.n with type b = EdgeH.e) =
    struct

        type a = NodeH.n
        type b = EdgeH.e

        module TOH = TreeOps
        module TOS = Chartree.TreeOps 
        module NodeS = Node.Standard
        module DH = MakeNormal (NodeH) (EdgeH) (TreeOps)
        module SH = MakeNormal (NodeS) (Edge.SelfEdge) (TOS)

        let report_mst = DH.report_mst

        let prebuilt = DH.prebuilt

        let replace_contents downpass uppass get_code nodes data ptree =
            let nt = { (Ptree.empty data) with Ptree.tree = ptree.Ptree.tree } in
            nodes 
                --> List.fold_left
                        (fun nt node ->
                            Ptree.add_node_data (get_code node) node nt) nt
                --> downpass
                --> uppass

        let from_s_to_h = replace_contents TOH.downpass TOH.uppass NodeH.taxon_code 
        let from_h_to_s = replace_contents TOS.downpass TOS.uppass NodeS.taxon_code

        let build_initial_trees trees data n b =
            let has_dyn = Data.has_dynamic data in
            let has_lik = Data.has_static_likelihood data in
            if has_dyn || has_lik then
                DH.build_initial_trees trees data n b
            else
                let s_nodes = List.map NodeH.to_other n in
                let trees = Sexpr.map (from_h_to_s s_nodes data) trees in
                let trees = SH.build_initial_trees trees data s_nodes b in
                Sexpr.map (from_s_to_h n data) trees

    end

(* Fast heuristics used *)


