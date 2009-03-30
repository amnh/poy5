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

let () = SadmanOutput.register "AllDirChar" "$Revision: 1616 $"

module IntSet = All_sets.Integers

let debug_profile_memory = false
let debug_node_fn = false
let debug_adjust_fn = false
let debug_clear_subtree = false
let debug_join_fn = false
let debug_branch_fn = false
let debug_cost_fn = false
let debug_uppass_fn = false

let current_snapshot x = 
    if debug_profile_memory then MemProfiler.current_snapshot x
    else ()

module F : Ptree.Tree_Operations with 
type a = AllDirNode.AllDirF.n
with type b = AllDirNode.OneDirF.n = struct

    type a = AllDirNode.AllDirF.n
    type b = AllDirNode.OneDirF.n
    type phylogeny = (a, b) Ptree.p_tree

    let (-->) a b = b a
    
    let force_node x = AllDirNode.force_val x.AllDirNode.lazy_node

    let error_user_message format = 
        Printf.ksprintf (Status.user_message Status.Error) format
    let info_user_message format = 
        Printf.ksprintf (Status.user_message Status.Information) format
    let failwithf format = 
        Printf.ksprintf (failwith) format

    (* process tree data to find branch lengthes *)
    let hashdoublefind data tree node_ids = 
        let transform_keys ntable keylation table =
            Hashtbl.iter 
                (fun name length ->
                    try Hashtbl.add ntable (Hashtbl.find keylation name) length
                    with | Not_found ->
                        if debug_branch_fn then
                            error_user_message "Couldn't find character name %s" name
                        else ()
                ) table;
            ntable
        in
        match tree.Tree.tree_name with
        | Some tree_n ->
            let ntable = Hashtbl.create 27 in
            let tree_n = String.uppercase tree_n in
           (try let tree_t = Hashtbl.find data.Data.branches (tree_n) in
                let res = List.fold_left
                    (fun acc node_id ->
                        try let node_n = Hashtbl.find (tree.Tree.names) node_id in
                            let node_n = String.uppercase node_n in
                            try let char_n = Hashtbl.find tree_t node_n in
                                let tbl = transform_keys (Hashtbl.create 27)
                                          data.Data.character_names char_n in
                                let () = Hashtbl.replace ntable node_id tbl in
                                true
                            with | Not_found ->
                                if debug_branch_fn then
                                    error_user_message "Couldn't find node name %s" node_n;
                                    false || acc
                        with | Not_found ->
                            if debug_branch_fn then
                                error_user_message "Couldn't find node id %d" node_id;
                                false || acc)
                    false
                    node_ids
                in
                if res then Some ntable else None

            with | Not_found ->
                if debug_branch_fn then
                    error_user_message "Couldn't find tree name %s" tree_n; None)
        | None -> None

(*  Creates a lazy edge which is the median between the data of the vertices
    with codes [a] and [b] in the tree [ptree].
    
    times in opposite directions SHOULD be set except when [root] is true. In
    this case, median_w_times is used and the tree is updated with the new
    directions calculated.
    
    root_node is used when the edge data is precomputed (during a break)
*)
    let create_lazy_edge ?branches root root_node adjusted ptree a b = 
        if debug_node_fn then
            info_user_message "Creating lazy edge between %d and %d" a b
        else ();

        if root then 
            let aa,ab = match Ptree.get_node a ptree with
                | Tree.Leaf (_,a_p) -> assert(a_p = b);
                    None,None
                | (Tree.Interior (_,a_p,ac1,ac2)) as an -> 
                    let ch1, ch2 = Ptree.other_two_nbrs b an in
                    Some (Ptree.get_node_data ch1 ptree),
                    Some (Ptree.get_node_data ch2 ptree)
                | Tree.Single _ -> failwith "create_lazy_edge of Singleton"
            and ba,bb = match Ptree.get_node b ptree with
                | Tree.Leaf (_,b_p) -> assert (b_p = a);
                    None,None
                | (Tree.Interior (_,b_p,bc1,bc2)) as bn ->
                    let ch1, ch2 = Ptree.other_two_nbrs a bn in
                    Some (Ptree.get_node_data ch1 ptree),
                    Some (Ptree.get_node_data ch2 ptree)
                | Tree.Single _ -> failwith "create_lazy_edge of Singleton"
            in
            let a_node,b_node,lazy_node =
                    AllDirNode.create_root ?branches
                            (Ptree.get_node_data a ptree) aa ab
                            (Ptree.get_node_data b ptree) ba bb
                            root_node
            in
            lazy_node,
                ptree --> Ptree.add_node_data a a_node
                      --> Ptree.add_node_data b b_node
        else begin
            (* creates root of an edge with time data *)
            assert( match root_node with | None -> true | _ -> false );
            let dataa = Ptree.get_node_data a ptree
            and datab = Ptree.get_node_data b ptree in
            AllDirNode.create_root_w_times adjusted dataa datab,ptree
        end

    (* Creates a valid vertex that only has the downpass information *)
    let create_lazy_interior_down ?branches ptree code a b =
        if debug_node_fn then
            match code with
            | Some x ->
                info_user_message
                    "Creating lazy interior down (%d) between %d and %d" x a b
            | None ->
                info_user_message
                    "Creating lazy interior down (?) between %d and %d" a b
        else ();
        let a_nd = Ptree.get_node_data a ptree 
        and b_nd = Ptree.get_node_data b ptree in
        AllDirNode.AllDirF.median ?branches code None a_nd b_nd

    (* Creates a valid vertex that has the downpass and uppass information.
    * [a] and [b] have to be the component of the currently calculated subtree
    * that is still valid. *)
    let create_lazy_interior_up ptree code a b c =
        if debug_node_fn then
            info_user_message 
                "Creating lazy interior up of %d with c:%d, c:%d and p:%d" code a b c
        else ();
        (* interior nodes should have all this information *)
        let cur_data = try Ptree.get_node_data code ptree with
            | Not_found -> create_lazy_interior_down ptree (Some code) a b
        and a_data = try Ptree.get_node_data a ptree with
            | Not_found -> failwith "AllDirChar.create_lazy_interior_up: need child data1"
        and b_data = try Ptree.get_node_data b ptree with
            | Not_found -> failwith "AllDirChar.create_lazy_interior_up: need child data2"
        and c_data = try Ptree.get_node_data c ptree with
            | Not_found -> failwith "AllDirChar.create_lazy_interior_up: need parent"
        in
        AllDirNode.AllDirF.uppass_heuristic c cur_data a_data b_data c_data

    let create_root a b (ptree : phylogeny) =
        let make_internal () = 
            let norm = 
                Tree.normalize_edge (Tree.Edge (a, b)) ptree.Ptree.tree 
            in
            let node = Ptree.get_edge_data norm ptree in
            let nnode = AllDirNode.force_val node in
            let node = {AllDirNode.lazy_node = node; code = 0; dir = Some (a,b)} in
            let node = { AllDirNode.unadjusted = [node]; adjusted = [node] } in
            { 
                Ptree.root_median = Some ((`Edge (a, b)), node);
(**)            component_cost = Node.Standard.tree_cost None nnode;
                adjusted_component_cost = Node.Standard.tree_cost None nnode;
            }
        in
        match Ptree.get_node a ptree with
        | Tree.Leaf (_, x) ->
                assert (x = b);
                make_internal ();
        | Tree.Interior (_, x, y, z) ->
                assert ((x = b) || (y = b) || (z = b));
                make_internal ();
        | Tree.Single _ ->
                let root = Ptree.get_node_data a ptree in
                { 
                    Ptree.root_median = Some ((`Single a), root);
                    component_cost = AllDirNode.AllDirF.tree_cost None root;
                    adjusted_component_cost = AllDirNode.AllDirF.tree_cost None root;
                }

    let print_roots ptree = 
        let roots = Ptree.get_roots ptree in
        let _ = Printf.printf "\nCurrent Root ids: " in
        let _ = List.iter (fun x -> match x.Ptree.root_median with
                | None -> Printf.printf "none"
                | Some e -> ( match e with
                        | ((`Edge (a,b)),_) -> Printf.printf "[%d -- %d]\t" a b
                        | ((`Single (a)),_) -> Printf.printf "[%d]\t" a )) roots
        in
        let _ = Printf.printf "\n%!" in ()

    let check_three_directions ptree = 
        let verify_one node = match Ptree.get_node node ptree with
            | Tree.Single _
            | Tree.Leaf _ ->
                let dat1 = (Ptree.get_node_data node ptree).AllDirNode.adjusted in
                let dat2 = (Ptree.get_node_data node ptree).AllDirNode.unadjusted in
                ((List.length dat1) = 1) && ((List.length dat2) = 1) &&
                    (match (List.hd dat1).AllDirNode.dir,
                           (List.hd dat2).AllDirNode.dir with 
                        | None,None -> true | _,_ -> false)
            | Tree.Interior _ ->
                let dat1 = (Ptree.get_node_data node ptree).AllDirNode.adjusted in
                let dat2 = (Ptree.get_node_data node ptree).AllDirNode.unadjusted in
                ((List.length dat1) = 3) && ((List.length dat2) = 3)
        in
        All_sets.Integers.fold
            (fun h ac1 ->
                match Ptree.get_node h ptree with
                | Tree.Leaf (a,b)
                | Tree.Interior (a,b,_,_) ->
                    Tree.pre_order_node_with_edge_visit_simple
                        (fun prev curr ac2 -> (verify_one curr) && ac2)
                        (Tree.Edge (a,b)) ptree.Ptree.tree ac1
                | Tree.Single a -> verify_one a
            )
            (Ptree.get_handles ptree)
            (true)


    (* We first define a function to check the total cost of a tree for
    * sequences and likelihood only! *)
    let check_cost new_tree handle new_root =
        (* Some characters are computed by the downpass, we extract that
        * portion of the cost, which is contained in the root.
        * not_single_character_cost holds the total cost of those characters *)
        let not_single_character_cost, root, root_edge = 
            (* set initial value to subtract from -- cost of likelihood *)
            let three_root, root_edge = 
                match new_root with
                | Some a -> 
                        a
                | None ->  
                        let x = Ptree.get_component_root handle new_tree in
                        match x.Ptree.root_median with
                        | Some (root_edge, root) -> root, root_edge
                        | None -> failwith "No root 1?"
            in
            match three_root.AllDirNode.unadjusted with
            | [root] ->
                    let root = AllDirNode.force_val root.AllDirNode.lazy_node in
                    (List.fold_left 
                    (fun acc y ->
                        acc +. 
                        (Node.total_cost_of_type y root)) 0. Node.not_to_single),
                    three_root, root_edge
            | _ -> failwith "What?"
        in 
        (* Other characters have their cost computed by adding up the length of
        * all of the branches. single_characters_cost is exactly that. *)
        let distance a b acc =
            let nda =
                (List.hd ((Ptree.get_node_data a
                new_tree).AllDirNode.adjusted)).AllDirNode.lazy_node
            and ndb = (List.hd ((Ptree.get_node_data b
            new_tree).AllDirNode.adjusted)).AllDirNode.lazy_node
            in
            let dist = 
                (Node.distance_of_type Node.has_to_single 100000.
                (AllDirNode.force_val nda) (AllDirNode.force_val ndb)) 
            in
            dist +. acc
        in
        let single_characters_cost =  
            match root_edge with
            | `Single _ -> 0.0
            | `Edge (a, b) ->
                    Tree.post_order_node_with_edge_visit_simple 
                    distance
                    (Tree.Edge (a, b))
                    new_tree.Ptree.tree 
                    (~-. (distance b a 0.0))
        in

        if debug_cost_fn then begin
            info_user_message "Single Character Cost: %f" single_characters_cost;
            info_user_message "Other Character Cost: %f" not_single_character_cost;
            info_user_message "Root Cost: %f" (AllDirNode.AllDirF.root_cost root)
        end;

        let res = 
            single_characters_cost +. not_single_character_cost +. 
            (AllDirNode.AllDirF.root_cost root)
        in
        res

    let check_cost_all_handles ptree = 
        let new_cost = 
            All_sets.Integers.fold (fun handle cost ->
            (check_cost ptree handle None) +. cost) 
            (Ptree.get_handles ptree) 0.0
        in
        new_cost

    let check_assertion_two_nbrs a b c =
        if a <> Tree.get_id b then true
        else 
            let _ = Status.user_message Status.Error c in
            false

    let get_pre_active_ref_code ptree = 
        let rec get_subtree parent current acc_pre_codes = 
            let pre_codes = 
                try                      
                    let a, b = 
                        let currentn = Ptree.get_node current ptree in 
                        assert (check_assertion_two_nbrs parent currentn "1");
                        Tree.other_two_nbrs parent currentn
                    in
                    let current_d = 
                        let current_3d = Ptree.get_node_data current ptree in
                        AllDirNode.not_with parent current_3d.AllDirNode.unadjusted
                    in

                    let _, pre_codes, _, _ = Node.get_active_ref_code 
                        (AllDirNode.force_val current_d.AllDirNode.lazy_node)
                    in 
                    let pre_child1 = get_subtree current a IntSet.empty in 
                    let pre_child2 = get_subtree current b IntSet.empty in
                    IntSet.union pre_codes (IntSet.union pre_child1 pre_child2)
                with
                | Invalid_argument _ -> IntSet.empty
            in 
            IntSet.union pre_codes acc_pre_codes
        in

        (* Now we define a function that can assign single sequences to the
        * connected component of a handle *)
        let get_handle handle pre_codes =
            let get_root_direction root = 
                match root.AllDirNode.unadjusted with
                | [x] -> AllDirNode.force_val (x.AllDirNode.lazy_node), x
                | _ -> failwith "get_handle at allDirChar"
            in
            let comp = Ptree.get_component_root handle ptree in
            match comp.Ptree.root_median with
            | Some ((`Edge (a, b)), rootg) ->
                  let root, rooth = get_root_direction rootg in
                  let r_pre, r_pre_child, _, __ = Node.get_active_ref_code root in
                  let prea_codes = get_subtree a b IntSet.empty in 
                  let preb_codes = get_subtree b a IntSet.empty in 
                  let new_pref_codes = IntSet.union (IntSet.union prea_codes preb_codes)
                      (IntSet.union r_pre r_pre_child)
                  in
                  IntSet.union pre_codes new_pref_codes
            | Some ((`Single a), rootg) ->
                  let root, rooth = get_root_direction rootg in 
                  let new_pref_codes, _, _, _ = Node.get_active_ref_code root in 
                  IntSet.union pre_codes new_pref_codes                      
            | _ -> failwith "Get_active_ref_code in allDirChar.ml"
        in 
        let pre_codes = 
            All_sets.Integers.fold get_handle
                ptree.Ptree.tree.Tree.handles IntSet.empty
        in
        pre_codes

    (* A function to assign a unique sequence on each vertex of the ptree in the
    * [AllDirNode.adjusted] field of the node. *)
    let assign_single keep_three ptree = 
        (* We first define a function that can traverse the tree and assign
        * a single sequence to each vertex on it. *)
        let pre_ref_codes = get_pre_active_ref_code ptree in  
        let fi_ref_codes = pre_ref_codes in 
        let rec assign_single_subtree parentd parent current ptree = 
            let current_d, initial_d =
                let tmp = Ptree.get_node_data current ptree in
                AllDirNode.not_with parent  tmp.AllDirNode.unadjusted, tmp
            in
            let nd, original = 
                current_d.AllDirNode.lazy_node
                --> AllDirNode.force_val 
(**)                --> fun x -> 
                        Node.to_single 
                        (pre_ref_codes, fi_ref_codes)
                        None parentd x, x
            in
            let nnd = 
                { current_d with AllDirNode.lazy_node = AllDirNode.lazy_from_val nd }
            in
            let oths =  if keep_three then
                            nnd::(List.filter (fun x -> 
                                      not (AllDirNode.has_code parent x))
                                      initial_d.AllDirNode.unadjusted)
                        else [nnd]
            in
            let final_d = { initial_d with AllDirNode.adjusted = oths } in

            assert ( ( (List.length initial_d.AllDirNode.unadjusted) ==
                       (List.length final_d.AllDirNode.adjusted) )
                     || (not keep_three) );

            let ptree = Ptree.add_node_data current final_d ptree in
            try 
                let a, b =
                    let currentn = Ptree.get_node current ptree in 
                    assert (check_assertion_two_nbrs parent currentn "2");
                    Tree.other_two_nbrs parent currentn 
                in
                ptree
                --> assign_single_subtree nd current a
                --> assign_single_subtree nd current b 
            with
            | Invalid_argument _ -> ptree
        in
        (* Now we define a function that can assign single sequences to the
        * connected component of a handle *)
        let assign_single_handle handle ptree =
            let get_root_direction root = 
                match root.AllDirNode.unadjusted with
                | [x] -> AllDirNode.force_val (x.AllDirNode.lazy_node), x
                | _ -> 
                        failwith "more than one root? \
                        AllDirChar.assign_single_handle 2"
            in
            let generate_root_and_assign_it rootg edge ptree =
                let a, b =
                    match edge with
                    | `Edge x -> x
                    | `Single a -> a, a
                in
                let root, rooth = get_root_direction rootg in
                let handle_node = 
                    (Ptree.get_node_data a ptree).AllDirNode.unadjusted 
                        --> AllDirNode.not_with b
                        --> (fun x -> AllDirNode.force_val x.AllDirNode.lazy_node)
                and other_node = 
                    (AllDirNode.not_with a (Ptree.get_node_data b
                    ptree).AllDirNode.unadjusted)
                    --> (fun x -> AllDirNode.force_val x.AllDirNode.lazy_node)
                in
                let root = 
                    Node.to_single (pre_ref_codes, fi_ref_codes) 
                                   (Some root) other_node handle_node
                in
                let rooti = 
                    [{ rooth with
                        AllDirNode.lazy_node = AllDirNode.lazy_from_val (root)
                    }]
                in
                let readjusted = { rootg with AllDirNode.adjusted = rooti} in

                let ptree = Ptree.assign_root_to_connected_component 
                        handle (Some (edge, readjusted)) 
                        (Node.Standard.tree_cost None root) None ptree
                in
                ptree,root,readjusted
            in
            let comp = Ptree.get_component_root handle ptree in
            match comp.Ptree.root_median with
            | Some ((`Edge (a, b)) as edge, rootg) ->
                    let ptree, root, readjusted = 
                        generate_root_and_assign_it rootg edge ptree 
                    in
                    ptree
                    --> assign_single_subtree root b a 
                    --> assign_single_subtree root a b 
                    --> (fun ptree ->
                        Ptree.assign_root_to_connected_component 
                        handle 
                        (Some ((`Edge (a, b)), readjusted))
                        comp.Ptree.component_cost
                        (Some (check_cost ptree handle None))
                        (* handle is already added in generate_root_and_assign *)
                        ptree)
            | Some ((`Single a) as edge, rootg) ->
                    let ptree, root, readjusted = 
                        generate_root_and_assign_it rootg edge ptree 
                    in
                    Ptree.add_node_data a readjusted ptree
            | None -> failwith "no root? AllDirChar.assign_single_handle"
        in
        (* Finally, we are ready to proceed on all the handles available *)
        let res =
            All_sets.Integers.fold assign_single_handle
            ptree.Ptree.tree.Tree.handles ptree 
        in
        (*
        assert ((Ptree.get_cost `Adjusted res) = check_cost res handle);
        *)
        res

    let unadjust ptree = ptree

    let refresh_all_edges ?data adjusted root_opt do_roots start_edge_opt ptree =

        (* refresh a single edge with root value *)
        let refresh_edge rhandle root_opt ((Tree.Edge (a,b)) as e) (acc,ptree) =
            if debug_uppass_fn then
                Printf.printf "Refreshing %d--%d as %s\n%!" a b 
                                (if rhandle then "a root edge" else "an edge")
            else ();
            let data,ptree = 
                try match data with
                    | Some data ->
                        begin 
                            match hashdoublefind data ptree.Ptree.tree [a;b] with
                            | Some x -> create_lazy_edge ~branches:x
                                            rhandle root_opt adjusted ptree a b 
                            | None -> create_lazy_edge rhandle root_opt adjusted ptree a b 
                        end
                    | None -> create_lazy_edge rhandle root_opt adjusted ptree a b
                with | err ->
                    let print_node = function
                        | Tree.Interior (u,v,w,x) ->
                            error_user_message "with neighbors %d, %d, %d," v w x;
                        | Tree.Leaf (x,y) -> 
                            error_user_message "with neighbor %d" x;
                        | Tree.Single _ -> ()
                    in
                    error_user_message "Failure attempting lazy_edge of %d and %d" a b;
                    print_node (Ptree.get_node a ptree);
                    print_node (Ptree.get_node b ptree);
                    raise err
            in
            (Tree.EdgeMap.add e data acc,ptree)
        (* perform uppass heuristic on a node *)
        and add_vertex_pre_order prev code (ptree:phylogeny) =
            match Ptree.get_node code ptree with
            | Tree.Single _ -> ptree
            | Tree.Leaf (_, p) ->
                    let this = Ptree.get_node_data code ptree
                    and parn = Ptree.get_node_data p ptree in
                    let leaf = AllDirNode.AllDirF.apply_time false this parn in
                    Ptree.add_node_data code leaf ptree
            | (Tree.Interior (_, par,a ,b)) as v -> 
                    let a,b = Tree.other_two_nbrs prev v in
                    let interior = create_lazy_interior_up ptree code a b prev in
                    Ptree.add_node_data code interior ptree
        in

        (* perform refresh on root node for uppass, to hold invariant that all
        * nodes have a parent with relevant data --both nodes will have all dirs *)
        if debug_uppass_fn then
            info_user_message "Performing Calculation on Root\n%!"
        else ();
        let ptree = 
            All_sets.Integers.fold 
                (fun h ptree ->
                    match start_edge_opt with
                    | Some (a,b) ->
                        let _,t = refresh_edge true root_opt
                                    (Tree.Edge (a,b)) 
                                    (Tree.EdgeMap.empty,ptree) in
                        t
                    | None ->
                        try 
                            begin match (Ptree.get_component_root h ptree).Ptree.root_median with
                                | Some ((`Edge (a,b)),c) -> 
                                    let _,t = refresh_edge true root_opt
                                                (Tree.Edge (a,b)) 
                                                (Tree.EdgeMap.empty,ptree) in
                                    t
                                | None
                                | Some _ -> ptree
                            end
                        with | Not_found -> 
                            begin match Ptree.get_node h ptree with
                                | Tree.Leaf (a,b)
                                | Tree.Interior (a,b,_,_) -> 
                                    let _,t = refresh_edge true root_opt
                                                (Tree.Edge (a,b)) 
                                                (Tree.EdgeMap.empty,ptree) in
                                    t
                                | Tree.Single _ -> ptree
                            end
                )
                ptree.Ptree.tree.Tree.handles
                ptree
        in
        (* perform uppass heuristic --fill all directions *)
        current_snapshot "AllDirChar.refresh_all_edges uppass heuristic";
        if debug_uppass_fn then
            info_user_message "Performing Uppass Heurisitic\n%!"
        else ();
        let ptree =
            All_sets.Integers.fold
                (fun h ptree ->
                    match start_edge_opt with
                    | Some (a,b) ->
                        Tree.pre_order_node_with_edge_visit_simple
                                add_vertex_pre_order
                                (Tree.Edge (a,b))
                                ptree.Ptree.tree ptree
                    | None ->
                        try 
                            begin match (Ptree.get_component_root h ptree).Ptree.root_median with
                                | Some ((`Edge (a,b)),c) -> 
                                    Tree.pre_order_node_with_edge_visit_simple
                                            add_vertex_pre_order
                                            (Tree.Edge (a,b))
                                            ptree.Ptree.tree ptree
                                | None
                                | Some _ -> ptree
                            end
                        with | Not_found -> 
                            begin match Ptree.get_node h ptree with
                                | Tree.Leaf (a,b)
                                | Tree.Interior (a,b,_,_) -> 
                                    Tree.pre_order_node_with_edge_visit_simple
                                            add_vertex_pre_order
                                            (Tree.Edge (a,b))
                                            ptree.Ptree.tree ptree
                                | Tree.Single _ -> ptree
                            end
                )
                ptree.Ptree.tree.Tree.handles
                ptree in
        (* fill in roots for all edges *)
        current_snapshot "AllDirChar.refresh_all_edges internal fold";
        if do_roots then begin
            if debug_uppass_fn then
                info_user_message "Performing Refresh on all edges\n"
            else ();
            let new_edges,ptree =
                Tree.EdgeSet.fold
                        (refresh_edge false None)
                        ptree.Ptree.tree.Tree.d_edges
                        (Tree.EdgeMap.empty,ptree) in
            { ptree with Ptree.edge_data = new_edges }
        end else ptree

    (** functions to return the adjusted and unadjusted values **)
    let get_single, get_unadjusted =
        let general_get f parent node =
            let nd = AllDirNode.not_with parent (f node) in
            AllDirNode.force_val nd.AllDirNode.lazy_node
        in
        (general_get (fun x -> x.AllDirNode.adjusted)),
        (general_get (fun x -> x.AllDirNode.unadjusted))

    (** refresh root of all trees *)
    let refresh_roots ptree =
        let new_roots =
            All_sets.Integers.fold (fun x acc ->
                let root = 
                    match Ptree.get_node x ptree with
                    | Tree.Leaf (a, b)
                    | Tree.Interior (a, b, _, _) ->
                            create_root a b ptree
                    | Tree.Single _ ->
                            create_root x x ptree
                in
                All_sets.IntegerMap.add x root acc
            ) ptree.Ptree.tree.Tree.handles All_sets.IntegerMap.empty
        in
        {
            ptree with 
            Ptree.component_root = new_roots;
        }

    let get_active_ref_code tree =
        let get_active_ref_code parent node = 
            Node.get_active_ref_code (get_unadjusted parent node)
        in
        let get_active_ref_code_handle handle (pre, fi) =
            let leaf_handler parent node _ =
                let node_data = Ptree.get_node_data node tree in
                let _, _, fi, fi_child = 
                    get_active_ref_code parent node_data 
                in
                Some (IntSet.empty, (IntSet.union fi fi_child))
            and node_handler par node a b =
                let extract = function
                    | Some x -> x
                    | None -> assert false
                in
                let (apre, afi) = extract a 
                and (bpre, bfi) = extract b in
                let node_data = Ptree.get_node_data node tree in
                let _, pre, _, fi = get_active_ref_code par node_data in
                Some (IntSet.union (IntSet.union apre bpre) pre,
                IntSet.union (IntSet.union afi bfi) fi)
            in
            let pre, fi, root =
                match (Ptree.get_component_root handle tree).Ptree.root_median 
                with
                | None -> assert false
                | Some ((`Single _), root) ->
                        IntSet.empty, IntSet.empty, root
                | Some ((`Edge (a, b)), root) ->
                        match Ptree.post_order_node_with_edge_visit
                        leaf_handler node_handler (Tree.Edge (a, b)) tree None
                        with
                        | Some (apre, afi), Some (bpre, bfi) -> 
                                IntSet.union apre bpre, IntSet.union afi bfi,
                                root
                        | _ -> assert false
            in

            let fi = IntSet.filter (fun x -> not (IntSet.mem x pre)) fi in
            let rpre, rprech, rfi, _ = get_active_ref_code (-1) root in

            IntSet.union (IntSet.union pre rprech) rpre,
            IntSet.union fi rfi
        in
        All_sets.Integers.fold
        get_active_ref_code_handle
        (Ptree.get_handles tree)
        (All_sets.Integers.empty, All_sets.Integers.empty)

    let dump_tree f x ptree =
        let printf format = Printf.ksprintf (f) format in
        let traversal prev code _ = match Ptree.get_node code ptree with
            | Tree.Single x -> Tree.Continue, ()
            | Tree.Leaf (x,p) ->
                printf "Leaf %d -- %d:" x p;
                AllDirNode.AllDirF.dump_node (f) (Ptree.get_node_data x ptree)
                                                 (Ptree.get_node_data p ptree);
                printf "\n";
                Tree.Continue, ()
            | Tree.Interior (x,p,_,_) ->
                printf "Node %d -- %d:" x p;
                AllDirNode.AllDirF.dump_node (f) (Ptree.get_node_data x ptree)
                                                 (Ptree.get_node_data p ptree);
                printf "\n";
                Tree.Continue, ()
        in
        Ptree.pre_order_node_visit traversal x ptree ()



    (* ------------------------------------------------------------------------ *)
    (* Now we define a function that can adjust all the vertices in the tree
    * to improve the overall cost of the tree, using only the
    * [AllDirNode.adjusted] field of each. *)
    let rec adjust_tree max_count nodes ptree =

        let mode = match !Methods.cost with
            | `Iterative x -> x
            | _ -> assert false
        and max_count = 
            match max_count with
            | Some x -> x
            | None -> max_int
        in

        (* We start by defining a function to adjust one node *)
        (* ptree has been changed to a ref, so this function is NOT THREAD SAFE *)
        let adjust_node chars_to_check ch1_k ch2_k parent_k mine_k p_ref =

            if debug_adjust_fn then
                info_user_message "Adjusting %d with %d,%d then %d" 
                                        mine_k ch1_k ch2_k parent_k 
            else ();
            
            (* we just have the keys, so we need to get the node data *)
            let gnd x = Ptree.get_node_data x !p_ref in
            let mine,modified = AllDirNode.AllDirF.readjust
                                            None (* adjusted nodes *)
                                            mode
                                            chars_to_check
                                            (gnd ch1_k)
                                            (gnd ch2_k)
                                            (gnd parent_k)
                                            (gnd mine_k)
            in
            if All_sets.Integers.is_empty modified then 
                modified,[]
            else
                let () = p_ref := !p_ref --> Ptree.add_node_data mine_k mine in
                modified, [ch1_k;ch2_k;parent_k]
        in
        
        (* add modified vertices in node_list to the set *)
        let add_vertices_affected node_list codes affected = 
            let add_one code affected =
                if All_sets.IntegerMap.mem code affected then
                    match All_sets.IntegerMap.find code affected with
                    | Some conts ->
                            let res = All_sets.Integers.union conts codes in
                            All_sets.IntegerMap.add code (Some res) affected
                    | None -> assert false
                else All_sets.IntegerMap.add code (Some codes) affected
            in
            List.fold_right (fun x acc -> add_one x acc) node_list affected
        in

        (* adjust the root and return the best cost *)
        let adjust_root handle root a b ptree = ptree (*
    
            if debug_adjust_fn then
                info_user_message "Adjusting root with %d,%d then None" a b
            else ();
            let sets = get_active_ref_code ptree (*pre_ref_codes,fi_ref_codes*)
            and a_nd = Ptree.get_node_data a ptree
            and b_nd = Ptree.get_node_data b ptree in
            
            let new_root = 
                (AllDirNode.AllDirF.median None (Some root) a_nd b_nd)
                --> (fun x -> AllDirNode.AllDirF.edge_iterator None x a_nd b_nd)
                --> (fun x -> AllDirNode.AllDirF.to_single (Some x) None a_nd None b_nd sets)
            in

            Ptree.assign_root_to_connected_component
                    handle 
                    (Some ((`Edge (a, b)),new_root ))
                    (check_cost ptree handle (Some (new_root , `Edge (a, b))))
                    None
                    ptree
        *)
        in

        (* Now we need to be able to adjust the root of the tree, and it's
        * cost, once we have finished adjusting every vertex in the tree *)
        let adjust_until_nothing_changes max_count start_ptree = 
            (* empty set *)
            let first_affected = 
                match nodes with
                | None -> 
                        All_sets.IntegerMap.map (fun x -> None)
                        ptree.Ptree.node_data 
                | Some items ->
                        items
            in

            (* adjusts vertices *)
            let adjust_vertices p_ref c2c affect prev curr a_ch b_ch =

                if not (All_sets.IntegerMap.mem curr c2c) then
                    a_ch || b_ch
                else (* it was modified, so get it and modify it *)
                    let c2c = All_sets.IntegerMap.find curr c2c in
                    let a,b = 
                        let currc = Ptree.get_node curr !p_ref in
                        match Ptree.get_node curr !p_ref with
                        | Tree.Interior (_,a,b,c) -> 
                            Tree.other_two_nbrs prev currc
                        | _ -> failwith "this is impossible"
                    in
                    
                    let codes, new_affected = adjust_node c2c a b prev curr p_ref in
                    let () = 
                        affect := add_vertices_affected new_affected codes !affect
                    in
                    (0 != (List.length new_affected)) || a_ch || b_ch
            in
            
            (* Post order traversal of internal nodes *)
            let adjust_loop affected p_ref a_ref handle changed = try
                let a,b,comp_root = 
                    match (Ptree.get_component_root handle !p_ref).Ptree.root_median with
                    | Some ((`Edge (a,b)),c) -> a,b,c
                    | None
                    | Some _ -> raise Not_found
                in
                let (a_ch),(b_ch) =
                    Ptree.post_order_node_with_edge_visit
                        (* F - leaf function: PREV CURR ACC *)
                            (fun prev curr acc -> acc)
                        (* G - interior function: PREV CURR aACC bACC *)
                            (adjust_vertices p_ref affected a_ref)
                            (Tree.Edge(a,b))
                             !p_ref
                             false
                in
                (* we need to iterate the Edge(a,b) and then update cost/ptree *)
                let () = p_ref := adjust_root handle comp_root a b !p_ref in
                changed || a_ch || b_ch

                with | Not_found -> false
            in

            (* recursive loop of for changes *)
            let rec iterator count prev_cost affected ptree =
                let changed,new_affected,new_ptree = 
                    let pref = ref ptree in    (* can't merge two ptrees *)
                    let aref = ref All_sets.IntegerMap.empty in (* map - no union operator *)
                    (* perform on each tree *)
                    let changed = All_sets.Integers.fold
                        (adjust_loop affected pref aref)
                        ptree.Ptree.tree.Tree.handles
                        (false)
                    in
                    changed,!aref,!pref
                in
                (* now ptree can be used normally *)
                let new_cost = check_cost_all_handles new_ptree in
                if (not changed) || (count = 1) then
                    ptree
                else if new_cost < prev_cost then
                    iterator (count - 1) new_cost new_affected new_ptree
                else
                    ptree
            in
            iterator max_count (check_cost_all_handles ptree) first_affected ptree
        in

        let set_handle_n_root_n_cost handle ptree =
            match (Ptree.get_component_root handle ptree).Ptree.root_median with
            | None -> failwith "Huh? AllDirChar.adjust_handle"
            | Some ((`Edge (a, b)), root) ->
                let sets = get_active_ref_code ptree
                and ad = Ptree.get_node_data a ptree
                and bd = Ptree.get_node_data b ptree in

                let root = 
                    let x = AllDirNode.AllDirF.median None None ad bd in
                    AllDirNode.AllDirF.to_single (Some x) None ad None bd sets
                in
                Ptree.assign_root_to_connected_component handle
                    (Some ((`Edge (a,b), root)))
                    (check_cost ptree handle None)
                    None
                    ptree
            | Some _ -> ptree
        in
        All_sets.Integers.fold (set_handle_n_root_n_cost)
                               (ptree.Ptree.tree.Tree.handles)
                               (adjust_until_nothing_changes max_count ptree)

    (* ------------------------------------------------------------------------ *)

    let verify_downpass x ptree : bool =
        let traversal prev code acc =
           match Ptree.get_node code ptree with
            | Tree.Single x 
            | Tree.Leaf (x, _) -> (Tree.Continue, acc)
            | (Tree.Interior (_, par, a, b)) as v ->
                let a, b, c, m =
                    match prev with
                    | Some prev ->
                        assert (check_assertion_two_nbrs prev v "4");
                        let a, b = Tree.other_two_nbrs prev v in
                        (Ptree.get_node_data a ptree,
                         Ptree.get_node_data b ptree,
                         Ptree.get_node_data prev ptree,
                         Ptree.get_node_data code ptree)
                    | None -> 
                        (Ptree.get_node_data a ptree,
                         Ptree.get_node_data b ptree,
                         Ptree.get_node_data par ptree,
                         Ptree.get_node_data code ptree)
                in
                (Tree.Continue, acc &&
                        (AllDirNode.verify_branch_lengths a b c m))
        in
        info_user_message "Verifying Branch Lengths";
        Ptree.pre_order_node_visit traversal x ptree true

    (** [internal_downpass] Traverse every vertex in the tree and assign the
     * downpass and uppass information using the lazy all direction nodes *)
    let internal_downpass ?data do_roots (ptree : phylogeny) : phylogeny =
        info_user_message "Internal Downpass"; 

        (* function to process tree->node->charactername to int->float hashtbl
        * for all the node ids passed --(multiple node_id capability for uppass) *)

        let ptree = 
            (* A function to add  the vertices using a post order traversal 
            * from the Ptree library. *)
            let add_vertex_post_order prev code ptree =
                current_snapshot
                "AllDirChar.internal_downpass.add_vertex_post_order";
                match Ptree.get_node code ptree with
                | Tree.Single _
                | Tree.Leaf (_, _) -> 
                        assert (All_sets.IntegerMap.mem code ptree.Ptree.node_data);
                        ptree
                | (Tree.Interior (_, par, a, b)) as v ->
                        let a,b = Tree.other_two_nbrs prev v in
                        if debug_branch_fn then
                            info_user_message 
                                "Adding Vertex %d post Order: (%d,%d) and %d%!" 
                                                code a b prev
                        else ();

                        let interior = match data with
                            | Some data ->
                                (match hashdoublefind data ptree.Ptree.tree [a;b] with
                                 | Some x -> create_lazy_interior_down 
                                                ~branches:x ptree (Some code) a b
                                 | None ->
                                    create_lazy_interior_down ptree (Some code) a b)
                            | None ->
                                create_lazy_interior_down ptree (Some code) a b
                        in
                        Ptree.add_node_data code interior ptree
            in
            (* A function to add the vertices using a post order traversal from
            * the Ptree library *)
            All_sets.Integers.fold 
                (fun x (ptree:phylogeny) ->
                    try begin
                        match (Ptree.get_component_root x ptree).Ptree.root_median with
                        | Some ((`Edge (a,b)),c) ->
                            Tree.post_order_node_with_edge_visit_simple
                                add_vertex_post_order
                                (Tree.Edge (a,b)) 
                                ptree.Ptree.tree ptree
                        | None | Some _ -> ptree
                    end with | Not_found -> 
                        begin match Ptree.get_node x ptree with
                        | Tree.Leaf (a,b)
                        | Tree.Interior (a,b,_,_) ->
                            Tree.post_order_node_with_edge_visit_simple
                                add_vertex_post_order
                                (Tree.Edge (a,b)) 
                                ptree.Ptree.tree ptree
                        | Tree.Single _ -> ptree
                    end
                )
                ptree.Ptree.tree.Tree.handles
                ptree
        in
        if do_roots then
            ptree --> refresh_all_edges ?data false None true None --> refresh_roots
        else 
            refresh_all_edges ?data false None true None ptree

    let clear_internals t = 
        let tree = internal_downpass false t in
        assert ( check_three_directions tree );
        tree

    let blindly_trust_downpass ptree 
        (edges, handle) (cost, cbt) ((Tree.Edge (a, b)) as e) =
        let data = Ptree.get_edge_data e ptree in
        let c = AllDirNode.OneDirF.tree_cost None data in
        if abs_float cost > abs_float c then 
            let data = 
                [{ AllDirNode.lazy_node = data; dir = None; code = -1 }] 
 (**)       in
            let data = { AllDirNode.unadjusted = data; adjusted = data } in
            let comp = Some ((`Edge (a, b)), data) in
            c, 
            Lazy.lazy_from_fun (fun () ->
                Ptree.set_component_cost c None comp handle ptree)
        else (cost, cbt)

(*
    let blindly_trust_single ptree 
        (edges, handle) (cost, cbt) ((Tree.Edge (a, b)) as e) =
        let data = Ptree.get_edge_data e ptree in
        let tree =
            let c = AllDirNode.OneDirF.tree_cost None data in
            let data = 
(**)            [{ AllDirNode.lazy_node = data; dir = None; code = -1 }] 
            in
            let data = { AllDirNode.unadjusted = data; adjusted = data } in
            let comp = Some ((`Edge (a, b)), data) in
            Ptree.set_component_cost c None comp handle ptree
        in
        let tree = assign_single false tree in
        let c = Ptree.get_cost `Adjusted tree in
        if abs_float cost > abs_float c then c, lazy tree
        else (cost, cbt)

    let blindly_trust_adjusted ptree 
        (edges, handle) (cost, cbt) ((Tree.Edge (a, b)) as e) =
        let data = Ptree.get_edge_data e ptree in
        let tree =
            let c = AllDirNode.OneDirF.tree_cost None data in
            let data = 
                [{ AllDirNode.lazy_node = data; dir = None; code = -1 }] 
(**)        in
            let data = { AllDirNode.unadjusted = data; adjusted = data } in
            let comp = Some ((`Edge (a, b)), data) in
            Ptree.set_component_cost c None comp handle ptree
        in
        let tree = tree --> assign_single false --> 
                    adjust_tree None None --> refresh_all_edges false in
        let c = Ptree.get_cost `Adjusted tree in
        if abs_float cost > abs_float c then c, lazy tree
        else (cost, cbt)
*)

    let general_pick_best_root selection_method ptree = 
        let edgesnhandles = 
            All_sets.Integers.fold 
            (fun handle acc ->
                ((Ptree.get_pre_order_edges handle ptree), handle) :: acc)
            ptree.Ptree.tree.Tree.handles 
            []
        in
        let process ptree (edges, handle) =
            let current_root_of_tree =
                let r = Ptree.get_component_root handle ptree in
                match r.Ptree.root_median with
                | Some (`Single _, _) 
                | None -> 0., lazy ptree
                | Some ((`Edge e), n) ->
                        r.Ptree.adjusted_component_cost, lazy ptree
            in
            let _, ptree =
                List.fold_left (selection_method ptree (edges, handle)) current_root_of_tree
                    (List.sort (fun (Tree.Edge (a, b)) (Tree.Edge (c, d)) ->
                        match c - a with
                        | 0 -> d - b
                        | x -> x)
                        edges)
            in
            Lazy.force_val ptree
        in 
        List.fold_left process ptree edgesnhandles 

    let pick_best_root ptree = general_pick_best_root blindly_trust_downpass ptree

    let downpass ?data ptree = 
        current_snapshot "AllDirChar.downpass a";
        let res = 
            match !Methods.cost with
            | `Exhaustive_Strong
            | `Exhaustive_Weak
            | `Normal_plus_Vitamines
            | `Normal -> internal_downpass ?data true ptree
            | `Iterative (`ApproxD iterations) ->
                    ptree --> internal_downpass ?data true -->
                        pick_best_root -->
                        assign_single false --> 
                        adjust_tree iterations None -->
                        refresh_all_edges true None true None
            | `Iterative (`ThreeD iterations) ->
                    ptree
                    --> internal_downpass ?data true
                    --> pick_best_root
                    --> assign_single false
                    --> adjust_tree iterations None
                    --> refresh_all_edges true None true None
        in
        assert( check_three_directions res );
        current_snapshot "AllDirChar.downpass b";
        res

    let uppass ptree = 
        let tree = match !Methods.cost with
            | `Exhaustive_Strong
            | `Exhaustive_Weak
            | `Normal_plus_Vitamines
            | `Normal -> 
                    ptree --> pick_best_root 
                          --> assign_single true
            | `Iterative (`ApproxD _)
            | `Iterative (`ThreeD _) -> ptree
        in
        assert( check_three_directions tree );
        tree

    let rec clear_subtree v p ptree = 
        if debug_clear_subtree then
            info_user_message "Clearing vertex %d with parent %d." v p
        else ();
        match Ptree.get_node v ptree with
        | Tree.Leaf _ | Tree.Single _ -> ptree
        | (Tree.Interior (_, a, b, c)) as vn ->
                let uno,dos = Tree.other_two_nbrs p vn in
                p --> create_lazy_interior_up ptree v uno dos
                  --> (fun x -> Ptree.add_node_data v x ptree)
                  --> clear_subtree uno v
                  --> clear_subtree dos v

    (* reset the data in the removed direction, by running essentially an uppass
     * heuristic to fill in the other data starting at (a,b). edge data is a
     * previous node data that can be used for edge data, and direction
     * information, IT MUST HAVE (a,b) IN IT'S DIRECTIONS. *)
    let clear_up_over_edge (a, b) edge_data_opt ptree =
        assert( match edge_data_opt with
                | None -> true
                | Some x ->
                    try let _ = AllDirNode.with_both a b x.AllDirNode.adjusted in true
                    with | Not_found -> false );

        (* lets move the root to this edge, that way a simple refresh all edges
        * will take care of the missing node and update all the uppass data *)
    match edge_data_opt with
    | None -> ptree --> clear_subtree b a --> clear_subtree a b
    | Some edge ->
        let edge = (* use this as edge data for likelihood *)
            let single = [AllDirNode.with_both a b edge.AllDirNode.adjusted] in
            { AllDirNode.unadjusted = single; adjusted = single }
        in
        refresh_all_edges true None false (Some (a,b)) ptree
        (* ONLY LIFT DATA FOR LIKELIHOOD 
            refresh_all_edges true (Some edge) false (Some (a,b)) ptree
        *)

    let set_clade_root (ptree : phylogeny) ex_sister root handle = 
        match Ptree.get_node handle ptree with
        | Tree.Interior (handle, sister, _, _)
        | Tree.Leaf (handle, sister) ->
                let node, snode = 
                    let a = AllDirNode.not_with ex_sister root.AllDirNode.unadjusted in
                    { AllDirNode.unadjusted = [a]; adjusted = [a]}, a
                in
                let ptree = 
                    Ptree.add_edge_data (Tree.Edge (handle, sister))
                    snode.AllDirNode.lazy_node ptree
                in
(**)            let snode = force_node snode in
                let cost = Node.Standard.tree_cost None snode in
                Ptree.set_component_cost cost None
                (Some (`Edge (handle, sister), node)) handle ptree
        | Tree.Single _ ->
                let cost = AllDirNode.AllDirF.tree_cost None root in
                Ptree.set_component_cost cost None (Some (`Single handle, root)) 
                handle ptree

    let clean_ex_neighbor a b ptree = 
        let data = Ptree.get_node_data a ptree in
        let notwith_un = AllDirNode.not_with b data.AllDirNode.unadjusted in
        let node = { AllDirNode.unadjusted = [notwith_un]; adjusted = [notwith_un] } 
        in
        Ptree.add_node_data a node ptree

    let get_edge_n_force a b ptree =
        let data = Ptree.get_edge_data (Tree.Edge (a, b)) ptree in
        AllDirNode.force_val data

    let replace_topology tree ptree = { ptree with Ptree.tree = tree } 

    let add_component_root ptree handle root = 
        { ptree with 
        Ptree.component_root = All_sets.IntegerMap.add handle root
        ptree.Ptree.component_root }

    (* break_fn has type handle * int (node) -> tree -> tree * delta * aux_data *)
    let break_fn (tree_node, clade_node_id) (ptree : phylogeny) =
        let (Tree.Edge (tree_node, clade_node_id)) as edge = 
            Tree.normalize_edge (Tree.Edge (tree_node, clade_node_id)) 
            ptree.Ptree.tree 
        in
        (* We want to know the cost of the tree, so we force the calculation of
        * the downpass all the way down to the place of the breakage *)
        let prev_cost = 
            let edge_res = get_edge_n_force tree_node clade_node_id ptree in
            Node.Standard.tree_cost None edge_res
        in

        (* Figure out the cost of the broken tree *)
        let new_cost = 
            let clade_node = Ptree.get_node_data clade_node_id ptree
            and tree_node_dir = 
                ptree --> Ptree.get_node_data tree_node 
                    --> fun x -> x.AllDirNode.unadjusted
                    --> AllDirNode.not_with clade_node_id 
                    --> force_node
            in
            let clade_node_dir = 
                clade_node --> fun x -> x.AllDirNode.unadjusted
                    --> AllDirNode.not_with tree_node 
                    --> force_node
            in

            if debug_join_fn then begin
                info_user_message "Previous Cost: %f" prev_cost;
                info_user_message "New Costs around with %f and %f = %f"
                    clade_node_dir.Node.total_cost tree_node_dir.Node.total_cost
                    (clade_node_dir.Node.total_cost +. tree_node_dir.Node.total_cost)
            end else ();
            
            clade_node_dir.Node.total_cost +. tree_node_dir.Node.total_cost
        in

        (* Break the topology and update the data *)
        let ptree, tree_delta, clade_handle, tree_handle =
            (* A function that takes one side of a tree delta and updates the
            * tree's data using that information *)
            let update_break_delta delta ptree = 
                match delta with
                | `Edge (rem, l1, l2, _) ->
                    assert ((tree_node = rem) || (clade_node_id = rem));
                    let old_data = Ptree.get_node_data rem ptree in
                    ptree --> clean_ex_neighbor l1 rem
                          --> clean_ex_neighbor l2 rem
                          --> Ptree.remove_node_data rem
                          --> clear_up_over_edge (l1, l2) (Some old_data)
                | _ -> ptree
            in
            (* Perform the topology break *)
            let nbt, ((left_delta, right_delta) as tree_delta) = 
                Tree.break (tree_node, clade_node_id) ptree.Ptree.tree 
            in
            let tree_handle, clade_handle = 
                Tree.get_break_handles tree_delta nbt 
            in

            (* Update the actual contents of the tree *)
            let ptree =
                ptree --> Ptree.remove_root_of_component tree_handle
                    --> Ptree.remove_root_of_component clade_handle
                    --> Ptree.remove_edge_data edge
                    --> replace_topology nbt
                    --> (fun x -> assert( check_three_directions x); x)
                    --> update_break_delta left_delta
                    --> (fun x -> assert( check_three_directions x); x)
                    --> update_break_delta right_delta
                    --> (fun x -> assert( check_three_directions x); x)
                    --> refresh_all_edges true None true None
                    --> refresh_roots
            in
            ptree, tree_delta, clade_handle, tree_handle
        in
        (* Compare costs, and calculate the break delta *)
        let b_delta =
            if prev_cost = infinity && new_cost = infinity then 0.
            else 
                let rc, tc = 
                    let clade_root = 
                        let c = Ptree.get_component_root clade_handle ptree in
                        match c.Ptree.root_median with
                        | Some (_, b) -> b
                        | None -> failwith "AllDirChar.break_fn Huh?"
                    in
                    AllDirNode.AllDirF.root_cost clade_root, 
                    AllDirNode.AllDirF.total_cost None clade_root
                in
                (prev_cost -. (new_cost -. (rc +. ptree.Ptree.origin_cost))) -. tc
        in

        if debug_join_fn then
            info_user_message "Break Delta: %f" b_delta;

        let left, right =
            let extract_side x side =
                let component_root x =
                    let cr = Ptree.get_component_root x ptree in
                    match cr.Ptree.root_median with
                    | Some (_, x) -> x
                    | None -> assert false
                in
                { Ptree.clade_id = x; 
                clade_node = component_root x;
                topology_delta = side;}
            in
            let (left, right) = tree_delta in
            extract_side tree_handle left, extract_side clade_handle right
        in
        assert (left.Ptree.topology_delta = fst tree_delta);
        assert (right.Ptree.topology_delta = snd tree_delta);
        assert (
            let get_handle side = 
                match side.Ptree.topology_delta with
                | `Edge (_, a, _, _) -> 
                        Ptree.handle_of a ptree
                | `Single (a, _) ->
                        let res = Ptree.handle_of a ptree in
                        assert (a = res);
                        res
            in
            get_handle left <> get_handle right);
        {
            Ptree.ptree = ptree;
            tree_delta = tree_delta;
            break_delta = b_delta;
            left = left;
            right = right;
            incremental = [];
        }

    let get_other_neighbors (a, b) tree acc = 
        let add_one a b acc =
            match Ptree.get_node a tree with
            | Tree.Interior (_, x, y, z) ->
                    if x <> b then All_sets.IntegerMap.add x None acc
                    else All_sets.IntegerMap.add y None acc
            | _ -> acc
        in
        let beg = 
            match acc with
            | None -> All_sets.IntegerMap.empty 
            | Some x -> x
        in
        Some (beg --> add_one a b --> add_one b a)

    let reroot_fn force edge ptree =
        let Tree.Edge (h, n) = edge in
        let my_handle = Ptree.handle_of h ptree in
        let root = Ptree.get_component_root my_handle ptree in
        let ptree, _ = 
            ptree 
            --> Ptree.remove_root_of_component my_handle 
            --> Ptree.move_handle h 
        in
        let ptree = Ptree.fix_handle_neighbor h n ptree in
        let tree,inc = match !Methods.cost with
            | `Exhaustive_Strong
            | `Exhaustive_Weak
            | `Normal_plus_Vitamines
            | `Iterative `ApproxD _
            | `Normal -> 
                    let root = 
                        let new_roots = create_root h n ptree in
                        if force || 
                            (abs_float new_roots.Ptree.component_cost) < 
                            (abs_float root.Ptree.component_cost) then
                            new_roots
                        else root
                    in
                    add_component_root ptree h root, []
            | `Iterative `ThreeD _ -> 
                    add_component_root ptree h root, []
        in
        assert (check_three_directions tree);
        tree,inc

    let break_fn ((s1, s2) as a) b =
        assert( check_three_directions b );
        let res = match !Methods.cost with
        | `Iterative (`ApproxD _)
        | `Iterative (`ThreeD _)
        | `Exhaustive_Weak
        | `Normal_plus_Vitamines
        | `Normal ->
                let breakage = break_fn a b in
                let nt = refresh_all_edges true None true None
                                           (breakage.Ptree.ptree) in
                { breakage with 
                    Ptree.ptree = nt; }
        | `Exhaustive_Strong ->
                let breakage = break_fn a b in
                let nt = refresh_all_edges true None true None
                                           (breakage.Ptree.ptree) in
                { breakage with 
                    Ptree.ptree = nt;
                    incremental = []; 
                    break_delta = (Ptree.get_cost `Adjusted b) -. 
                                  (Ptree.get_cost `Adjusted nt);
                }
        in
        assert( check_three_directions res.Ptree.ptree );
        res

    (* join_fn must have type join_1_jxn -> join_2_jxn -> delta -> tree -> tree
    * *)
    let join_fn _ jxn1 jxn2 ptree =
        if debug_join_fn then
            info_user_message "Time to join! (%s) and (%s)"
                (match jxn1 with
                    | Tree.Single_Jxn x -> string_of_int x
                    | Tree.Edge_Jxn (x,y) -> (string_of_int x) ^","^ (string_of_int y)
                )
                (match jxn2 with
                    | Tree.Single_Jxn x -> string_of_int x
                    | Tree.Edge_Jxn (x,y) -> (string_of_int x) ^","^ (string_of_int y)
                )
        else ();

        let lift_data edge_l edge_r i_code ptree = 
            
            let node = AllDirNode.AllDirF.median (Some i_code) None
                            (Ptree.get_node_data edge_l ptree)
                            (Ptree.get_node_data edge_r ptree) in
            Ptree.add_node_data i_code node ptree
            (* ONLY FOR LIKELIHOOD 
            try
                let lr = Ptree.get_edge_data (Tree.Edge (edge_l,edge_r)) ptree in
                let node = 
                    [{
                        AllDirNode.dir = Some (edge_l,edge_r);
                        code = i_code;
                        lazy_node = lr;
                    }] in
                let node = 
                    {AllDirNode.adjusted = node; AllDirNode.unadjusted = node}
                in
                Ptree.add_node_data i_code node ptree
            with | Not_found -> failwithf "Cannot lift %d -- %d to %d"
                                                    edge_l edge_r i_code
            *)
        in

        let ret, ((td1,td2,rrd) as tree_delta) = Tree.join jxn1 jxn2 ptree.Ptree.tree in

        let v, h, ptree = match tree_delta with
            | (`Edge (v, a, b, _)), (`Single (h, true)), _ ->
                    let ptree = 
                        ptree --> Ptree.remove_node_data v 
                            --> clean_ex_neighbor a b
                            --> clean_ex_neighbor b a
                            --> Ptree.remove_root_of_component h
                            --> lift_data a b v
                    in
                    v, h, ptree
            | (`Single (v, _)), (`Single (h, true)), _ ->
                    v, h, Ptree.remove_root_of_component h ptree
            | (`Edge (v, c, d, _)), (`Edge (r, a, b, Some h)), _ ->
                    let ptree = 
                        ptree --> Ptree.remove_root_of_component h 
                            --> Ptree.remove_node_data r 
                            --> Ptree.remove_node_data v
                            --> clean_ex_neighbor c d
                            --> clean_ex_neighbor d c
                            --> clean_ex_neighbor a b
                            --> clean_ex_neighbor b a
                            --> lift_data a b r
                            --> lift_data c d v
                    in
                    r, v, ptree
            | (`Single (v, _)), (`Edge (r, a, b, Some h)), _ ->
                    let ptree = 
                        ptree --> Ptree.remove_root_of_component h 
                            --> Ptree.remove_node_data r 
                            --> clean_ex_neighbor a b
                            --> clean_ex_neighbor b a
                            --> lift_data a b r
                    in
                    r, h, ptree
            | _ -> failwith "Unexpected Chartree.join_fn"
        in
        let ptree = { ptree with Ptree.tree = ret } in
        let handle, parent = 
            let handle = Ptree.handle_of v ptree in 
            let parent = Ptree.get_parent handle ptree in
            handle, parent
        in
        let ptree = 
            ptree --> Ptree.remove_root_of_component handle 
                (* taken care of by uppass --> clear_up_over_edge (v, h) None *)
                --> refresh_all_edges false None true (Some (v,h))
        in
        let ptree = 
            add_component_root ptree handle (create_root handle parent ptree)
        in
        ptree, tree_delta

    let get_one side = 
        match side with
        | `Single (x, _) | `Edge (x, _, _, _) -> x

    let join_fn a b c d =
        let ptree, tdel = match !Methods.cost with
            | `Normal -> join_fn a b c d 
            | `Iterative (`ThreeD iterations)
            | `Iterative (`ApproxD iterations) ->
                let tree, delta = join_fn a b c d in
                let tree = adjust_tree iterations None
                           (assign_single false (pick_best_root tree)) in
                refresh_all_edges true None true None tree, delta
            | `Normal_plus_Vitamines
            | `Exhaustive_Weak
            | `Exhaustive_Strong ->
                let tree, delta = join_fn a b c d in
                uppass tree, delta 
        in
        assert( check_three_directions ptree );
        ptree,tdel

    type tmp = Edge of (int * int) | Clade of a 

    let cost_fn jxn1 jxn2 delta clade_data (tree : phylogeny) =
        let rec forcer edge =
            match edge with
            | Edge (a, b) ->
                    AllDirNode.force_val (Ptree.get_edge_data (Tree.Edge (a, b)) tree)
            | Clade x -> 
                    (match x.AllDirNode.unadjusted with
                    | [x] -> force_node x
                    | _ -> failwith "AllDirChar.cost_fn")
        in
        let clade_data = 
            match !Methods.cost with
            | `Iterative (`ThreeD _) ->
                    (match jxn2 with
                    | Tree.Single_Jxn _ -> forcer (Clade clade_data)
                    | Tree.Edge_Jxn (h, n) ->
                            let (Tree.Edge (h, n)) = 
                                Tree.normalize_edge 
                                (Tree.Edge (h, n)) tree.Ptree.tree
                            in
                            forcer (Edge (h, n)))
            | _ -> forcer (Clade clade_data)
        in
        match jxn1 with
        | Tree.Single_Jxn h ->
                let d = 
                    Node.Standard.distance 0.
                    (forcer (Clade (Ptree.get_node_data (Tree.int_of_id h)
                    tree)))
                    clade_data
                in
                Ptree.Cost d
        | Tree.Edge_Jxn (h, n) ->
                let (Tree.Edge (h, n)) = 
                    Tree.normalize_edge (Tree.Edge (h, n)) tree.Ptree.tree
                in
                let ndata = forcer (Edge (h, n)) in
                let d = Node.Standard.distance 0. clade_data ndata in
                Ptree.Cost d

    let cost_fn a b c d e =
        match !Methods.cost with
            | `Iterative (`ApproxD _) -> 
                    (match cost_fn a b c d e with 
                    | Ptree.Cost x -> Ptree.Cost (0.85 *. x)
                    | x -> x)
            | `Iterative `ThreeD _
            | `Exhaustive_Weak
            | `Normal_plus_Vitamines
            | `Normal -> cost_fn a b c d e 
            | `Exhaustive_Strong ->
                    let pc = Ptree.get_cost `Adjusted e in
                    let (nt, _) = join_fn [] a b e in
                    assert( check_three_directions nt );
                    Ptree.Cost ((Ptree.get_cost `Adjusted nt) -. pc)

 
    let root_costs tree = 
        let collect_edge_data edge node acc =
            let cost = AllDirNode.OneDirF.root_cost node in
            (edge, cost) :: acc
        in
        Tree.EdgeMap.fold collect_edge_data tree.Ptree.edge_data []

    let string_of_node _ = ""

    let features meth lst = 
        Chartree.features meth (("all directions", "true") :: lst)

    let incremental_uppass tree _ = tree

    let assign_final_states ptree =
        let assign_final_states_handle handle ptree =
            try
                let root_data, a, b = 
                    let rt = Ptree.get_component_root handle ptree in
                    match rt.Ptree.root_median with
                    | Some ((`Edge (a, b)), root) -> root, a, b
                    | Some _ -> failwith "Single vertex" (* Used down below *)
                    | None -> failwith "No root?"
                in
                let root_data c = 
                    (* We prepare a function to replace the taxon code for a
                    * meaningful one to start the uppass with on each side *)
                    match root_data.AllDirNode.unadjusted with
                    | [x] -> 
                            { root_data with 
                            AllDirNode.unadjusted = 
                                [{ x with AllDirNode.code = c }] }
                    | _ -> assert false
                in
                (* We move recursively up on and b calculating their final 
                * states *)
                let rec uppass grandparent_code parent_code parent_final vertex
                acc =
                    let my_data = Ptree.get_node_data vertex ptree in
                    match Ptree.get_node vertex acc with
                    | (Tree.Interior _) as nd ->
                            let a, b = Tree.other_two_nbrs parent_code nd in
                            let nda = Ptree.get_node_data a ptree
                            and ndb = Ptree.get_node_data b ptree in
                            let my_data =
                                AllDirNode.AllDirF.final_states grandparent_code 
                                parent_final my_data nda ndb 
                            in
                            acc
                            --> Ptree.add_node_data vertex my_data 
                            --> uppass (Some parent_code) vertex my_data a 
                            --> uppass (Some parent_code) vertex my_data b
                    | Tree.Leaf _ ->
                            let my_data = 
                                AllDirNode.AllDirF.final_states grandparent_code 
                                parent_final my_data my_data my_data 
                            in
                            Ptree.add_node_data vertex my_data acc
                    | Tree.Single _ -> acc
                in
                ptree --> uppass None a (root_data a) b 
                --> uppass None b (root_data b) a
            with
            | Failure "Single vertex" -> ptree
        in
        All_sets.Integers.fold assign_final_states_handle (Ptree.get_handles ptree) 
        ptree
    
    (** [create_branch_table table ptree] 
     *
     * Creates a hashtable with all the branch data like in Data.branches
     * Branch names are not used at this point, since they are probably
     * unreliable, and the nodes code is used in preference. The 
    *)
    let branch_table ptree =
        let tree_table = Hashtbl.create 13 in
        let create_branch_table handle _ = 
            (* let tree_code = (*  = handle *)
                match ptree.Ptree.tree.Tree.tree_name with
                | Some name -> name
                | None ->   incr tree_num;
                            "tree_"^(string_of_int !tree_num)
            in *)
            let single_node prev curr =
                let ntable =
                    if Hashtbl.mem tree_table 0 then
                        Hashtbl.find tree_table 0
                    else
                        let tbl = Hashtbl.create 13 in
                        Hashtbl.add tree_table 0 tbl;
                        tbl
                in 

                (* get the branch length of parent/prev for current *)
                let dat = AllDirNode.AllDirF.get_times_between 
                            (Ptree.get_node_data curr ptree)
                            (Ptree.get_node_data prev ptree) in
                let name_it x = match dat with | [_] -> `Single x | _ -> `Name in
                let () = List.iter
                    (fun (code,length) -> match length with
                        | Some length -> 
                           (* let name = (* name of char set *)
                                try Hashtbl.find (data.Data.character_nsets)
                                           (Hashtbl.find (data.Data.character_codes) code)
                                with | Not_found -> "node_"^(string_of_int code)
                              in
                           *)
                            Hashtbl.add ntable curr (name_it length)
                        | None -> () (* no data, forgetaboutit *)
                    ) dat in
                ()
            in
            let a,b =
                match (Ptree.get_component_root handle ptree).Ptree.root_median with
                | Some ((`Edge (a,b)),_) -> a,b
                | None | Some _ -> raise Not_found
            in
            let (),() = Ptree.post_order_node_with_edge_visit
                                (fun prev curr _ -> single_node prev curr)
                                (fun prev curr _ _ -> single_node prev curr)
                                (Tree.Edge (a,b))
                                ptree
                                ()
            in
            (* added everything already -- root? *) 
            ()
        in
        let () = All_sets.Integers.fold
                    (create_branch_table)
                    (Ptree.get_handles ptree)
                    ()
        in

        tree_table

    let to_formatter (atr : Xml.attributes) (data : Data.d) 
            (tree : (a, b) Ptree.p_tree) : Xml.xml =

        let tree = assign_final_states tree in
        let pre_ref_codes, fi_ref_codes = get_active_ref_code tree in 
(*
        Utl.printIntSet pre_ref_codes;
        Utl.printIntSet fi_ref_codes;
*)
        let get_simplified parent x = 
            let nd = Ptree.get_node_data x tree in
            nd, get_unadjusted parent nd, get_single parent nd
        in
        let merger a b root = 
            (`Set [`Single root; `Single a; `Single b]) 
        and splitter parent a = get_unadjusted parent a, get_single parent a in
        (* Now we are ready to process the contents of the tree *)
        let rec subtree_to_formatter (pre, fi) cur par 
                ((node_parent, single_parent) as tmp2) : Xml.xml =
            match Ptree.get_node cur tree with
            | (Tree.Interior _) as nd ->
                    let cur_data = Ptree.get_node_data cur tree in
                    let ch1, ch2 = Ptree.other_two_nbrs par nd in
                    let ch1d, ch1u, ch1s = get_simplified cur ch1 
                    and ch2d, ch2u, ch2s = get_simplified cur ch2 in
                    let ((cur_data, cur_single) as tmp) = 
                        splitter par cur_data 
                    in 
                    let mine = 
                        Node.to_formatter_subtree (pre, fi) [] data tmp cur 
                        (ch1, ch1u) (ch2, ch2u) (Some tmp2)
                    in
                    let ch1 = subtree_to_formatter (pre, fi) ch1 cur tmp in
                    let ch2 = subtree_to_formatter (pre, fi) ch2 cur tmp in
                    ((RXML 
                        -[Xml.Trees.tree] 
                            {single mine} { single ch1 } 
                            { single ch2 } --) : Xml.xml)
            | (Tree.Leaf (_, par)) ->
                    let node_data = Ptree.get_node_data cur tree in
                    
                    let nodest = 
                        Node.to_formatter_single
                        (pre, fi) [] data 
                        (splitter par node_data) cur (Some tmp2)
                    in
                    (RXML -[Xml.Trees.tree] { single nodest }--)
            | (Tree.Single _) ->
                    let node_data = Ptree.get_node_data cur tree in
                    let nodest = 
                        Node.to_formatter_single
                        (pre, fi) [] data (splitter (-1) node_data) cur None
                    in
                    (RXML -[Xml.Trees.tree] { single nodest } --)
        in
        let handle_to_formatter (pre, fi) handle (recost, trees) =
            let r = Ptree.get_component_root handle tree in
            let recost, contents, attr =
                match r.Ptree.root_median with
                | Some ((`Edge (a, b)), root) -> 
                        let recost = 
                            let root = get_unadjusted (-1) root in
                            (Node.cmp_subtree_recost root) +. recost 
                        in
                        (* We override the root now to continue using the single
                        * assignment of the handle *)
                        let sroot, sa = 
                            let a = Ptree.get_node_data a tree in
                            let s = get_single (-1) a in
                            let root = get_unadjusted (-1) root in 
                            let s_root = Node.copy_chrom_map root s in 
                            (root, s_root), s
                        in
                        let a : Xml.xml = 
                            subtree_to_formatter (pre, fi) a b sroot
                        and b : Xml.xml = 
                            subtree_to_formatter (pre, fi) b a sroot
                        and froot : Xml.xml =
                            let handle = Ptree.get_node_data a tree 
                            and parent = Ptree.get_node_data b tree in
                            Node.to_formatter_subtree 
                            (pre, fi) [] data 
                            ((get_unadjusted (-1) root), sa) a (a, get_unadjusted b handle)
                            (b, get_unadjusted a parent) None
                        in
                        recost, (merger a b froot), 
                        [Xml.Trees.cost, `Float r.Ptree.component_cost]
                | Some ((`Single a), root) ->
                        let c1 : Xml.xml = 
                            let nd = splitter (-1) root in
                            subtree_to_formatter (pre, fi) a a nd
                        in
                        recost, (`Single c1),
                        [Xml.Trees.cost, `Float r.Ptree.component_cost]
                | None -> assert false
            in
            recost, 
            (((PXML -[Xml.Trees.tree] ([attr]) { contents }--)) ::
                trees)
        in
        let recost, trees =
            All_sets.Integers.fold 
            (handle_to_formatter (pre_ref_codes, fi_ref_codes)) 
            (Ptree.get_handles tree)
            (0., [])
        in
        let cost = Ptree.get_cost `Adjusted tree in
        (RXML -[Xml.Trees.forest] 
            ([Xml.Trees.recost] = [`Float recost])
            ([Xml.Trees.cost] = [`Float cost])
            ([atr])
            { set trees } --)

end

module CharScripting : CharacterScripting.S with type cs =
    CharacterScripting.Standard.cs with type n = AllDirNode.AllDirF.n = struct
        type cs = CharacterScripting.Standard.cs
        type n = AllDirNode.AllDirF.n

        type character_input_output = [ 
            | `Characters of cs Sexpr.t
            | `Floats of float Sexpr.t
        ]

        let distance = CharacterScripting.Standard.distance
        let median = CharacterScripting.Standard.median
        let convert =
            List.map 
            (fun x -> AllDirNode.force_val (List.hd x.AllDirNode.unadjusted).AllDirNode.lazy_node)

        let scriptchar_operations nodes = 
            CharacterScripting.Standard.scriptchar_operations (convert nodes)

        let filter_char_operations nodes =
            CharacterScripting.Standard.filter_char_operations (convert nodes)

        let extract_character x y =
            CharacterScripting.Standard.extract_character x 
            (List.hd (convert [y]))

        let character_operations = 
            CharacterScripting.Standard.character_operations
end
