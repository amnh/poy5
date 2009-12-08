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
let debug_branch_fn = false (* downpass and given branch lengths *)
let debug_cost_fn = false
let debug_uppass_fn = false
let debug_downpass_fn = false
let debug_single_assignment = false

let current_snapshot x = 
    if debug_profile_memory then MemProfiler.current_snapshot x
    else ()

module F : Ptree.Tree_Operations 
    with type a = AllDirNode.AllDirF.n
        with type b = AllDirNode.OneDirF.n = struct

    type a = AllDirNode.AllDirF.n
    type b = AllDirNode.OneDirF.n
    type phylogeny = (a, b) Ptree.p_tree

    let (-->) a b = b a

    let (=.) a b = 
        match classify_float ( a -. b ) with
        | FP_infinite | FP_nan | FP_normal -> false
        | FP_subnormal | FP_zero -> true
    
    let force_node x = AllDirNode.force_val x.AllDirNode.lazy_node

    let error_user_message format = 
        Printf.ksprintf (Status.user_message Status.Error) format
    let info_user_message format = 
        Printf.ksprintf (Status.user_message Status.Information) format
    let failwithf format = 
        Printf.ksprintf (failwith) format

    let create_partition ptree left_node right_node =
        let left,right =
            Ptree.post_order_node_with_edge_visit
                (fun p x a -> All_sets.Integers.add x a)
                (fun p x acc_l acc_r -> All_sets.Integers.union acc_l acc_r)
                (Tree.Edge (left_node,right_node) )
                ptree
                All_sets.Integers.empty
        in
        left (* either side can return, no reason for choosing left *)

    (* build table of branch lengths like the table in Data.branches for tree *)
    let update_branches ptree =
        let get_codestable data node1 node2 = 
            let codestimes = 
                AllDirNode.AllDirF.get_times_between
                    (Ptree.get_node_data node1 ptree)
                    (Ptree.get_node_data node2 ptree)
            and table = Hashtbl.create 1227
            and insert_set codes table time = 
                Array.iter
                    (fun code -> 
                        let name = Hashtbl.find data.Data.character_codes code in
                        Hashtbl.add table name time)
                    codes
            in
            List.iter
                (fun (codes,times) -> match times with
                    | Some x -> insert_set codes table x
                    | None -> () )
                codestimes;
            table
        in
        let name = match ptree.Ptree.tree.Tree.tree_name with
            | Some x -> String.uppercase x | None -> ""
        and treebranches = Hashtbl.create 1 (* doing this for only one tree *)
        and insert_function setmap edge =
            let Tree.Edge (left,right) = edge in
            let p1 = create_partition ptree left right
            and p2 = create_partition ptree right left in
            let codestable = get_codestable ptree.Ptree.data left right in
            setmap
                --> All_sets.IntSetMap.add p1 codestable
                --> All_sets.IntSetMap.add p2 codestable 
        in
        let () =
            Ptree.get_edges_tree ptree
                --> List.fold_left insert_function All_sets.IntSetMap.empty
                --> Hashtbl.add treebranches name
        in
        {ptree.Ptree.data with Data.branches = treebranches }
            --> (fun x -> {ptree with Ptree.data = x})

    (* process tree data to find branch lengths *)
    let hashdoublefind tree partitions : (int,(int,float) Hashtbl.t) Hashtbl.t option =
        (* converts table with character names to character ids *)
        let transform_keys to_table from_table char_name_tbl =
            Hashtbl.iter 
                (fun name length ->
                    try Hashtbl.add to_table (Hashtbl.find char_name_tbl name) length
                    with | Not_found ->
                        if debug_branch_fn then
                            error_user_message "Couldn't find character name %s" name;)
                from_table;
            to_table
        (* retrieve the name of the tree, or no name *)
        and t_name = match (tree.Ptree.tree).Tree.tree_name with
            | Some tree_n -> String.uppercase tree_n | None -> ""
        in
        (* test if that table exists, and convert each partition to the node id
            it belongs from, and the charcter_names to character ids *) 
        try
            let ret_table = Hashtbl.create 27 in
            let tree_tbl = Hashtbl.find (tree.Ptree.data).Data.branches t_name in
            let res = List.fold_left
                (fun acc (partition,node_id) ->
                    try 
                        let node_n = All_sets.IntSetMap.find partition tree_tbl in
                        let tbl = transform_keys (Hashtbl.create 27)
                                                 node_n
                                                 (tree.Ptree.data).Data.character_names
                        in
                        let () = Hashtbl.replace ret_table node_id tbl in
                        true
                    with | Not_found -> (false or acc) )
                false
                partitions
            in
            if res then Some ret_table else None
        (* return nothing if the node wasn't found *)
        with | Not_found ->
            if debug_branch_fn then
                error_user_message "Couldn't find tree name %s" t_name;
            None
    
    let using_likelihood ptree =
        All_sets.Integers.fold
            (fun handle acc ->
                acc && (AllDirNode.AllDirF.using_likelihood 
                                (Ptree.get_node_data handle ptree)))
            ptree.Ptree.tree.Tree.handles
            true

(*  Creates a lazy edge which is the median between the data of the vertices
    with codes [a] and [b] in the tree [ptree].
    
    times in opposite directions SHOULD be set except when [root] is true. In
    this case, median_w_times is used and the tree is updated with the new
    directions calculated.
    
    root_node is used when the edge data is precomputed (during a break)
*)
    let create_lazy_edge ?branches root root_node adjusted ptree a b = 
        if debug_node_fn then
            info_user_message 
                "Creating lazy edge between %d and %d with root %s"
                a b (if root then "true" else "false")
        else ();
        if root then 
            let aa,ab = match Ptree.get_node a ptree with
                | Tree.Leaf (a,a_p) -> assert(a_p = b);
                    None,None
                | (Tree.Interior (a,a_p,ac1,ac2)) as an -> 
                    let ch1, ch2 = Ptree.other_two_nbrs b an in
                    Some (Ptree.get_node_data ch1 ptree),
                    Some (Ptree.get_node_data ch2 ptree)
                | Tree.Single _ -> failwith "create_lazy_edge of Singleton"
            and ba,bb = match Ptree.get_node b ptree with
                | Tree.Leaf (b,b_p) -> assert (b_p = a);
                    None,None
                | (Tree.Interior (b,b_p,bc1,bc2)) as bn ->
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
            assert (root_node = None);
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
        let cur_data = Ptree.get_node_data code ptree
        and a_data = try Ptree.get_node_data a ptree with
            | Not_found -> failwith "AllDirChar.create_lazy_interior_up: need child data1"
        and b_data = try Ptree.get_node_data b ptree with
            | Not_found -> failwith "AllDirChar.create_lazy_interior_up: need child data2"
        and c_data = try Ptree.get_node_data c ptree with
            | Not_found -> failwith "AllDirChar.create_lazy_interior_up: need parent"
        in
        AllDirNode.AllDirF.uppass_heuristic c_data None cur_data a_data b_data

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

    let tree_size ptree = 
        let traversal a b =
            let l,r = Ptree.post_order_node_with_edge_visit
                        (fun p x a -> a) (* leafs *)
                        (fun par node accleft accright ->
                            let nd = Ptree.get_node_data node ptree in
                            let curr = AllDirNode.AllDirF.tree_size (Some par) nd in
                            accleft +. accright +. curr )
                        (Tree.Edge (a,b))
                        ptree
                        0.0
            and final_edge = Ptree.get_edge_data (Tree.Edge (a,b)) ptree in
            l +. r +. (AllDirNode.OneDirF.tree_size None final_edge)
        in
        let single_tree handle =
            try match (Ptree.get_component_root handle ptree).Ptree.root_median with
                | Some ((`Edge (a,b)),_) -> traversal a b
                | None | Some _ -> raise Not_found
            with | Not_found ->
                begin match Ptree.get_node handle ptree with
                    | Tree.Leaf (a,b)
                    | Tree.Interior (a,b,_,_) -> traversal a b
                    | Tree.Single _ -> 0.0
                end
        in
        let sum =
            try All_sets.Integers.fold
                    (fun h x -> x +. (single_tree h))
                    (Ptree.get_handles ptree)
                    0.0
            with | Not_found -> -1.0
        in
        sum

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
                | Some a -> a
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
            (* debug msg 
            Printf.printf " distance of (%d,%d) is %!" a b;
             debug msg*)
            let nda =
                (List.hd ((Ptree.get_node_data a
                new_tree).AllDirNode.adjusted)).AllDirNode.lazy_node
            and ndb = (List.hd ((Ptree.get_node_data b
            new_tree).AllDirNode.adjusted)).AllDirNode.lazy_node
            in
            let dist = 
                Node.distance_of_type Node.has_to_single 0.
                (AllDirNode.force_val nda) (AllDirNode.force_val ndb)
            in
            (* debug msg
            * Printf.printf " [%f/%f]\n %!" dist acc;
            * debug msg*)
            let tmp = dist +. acc in
            tmp
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
            info_user_message "Root Cost: %f" (AllDirNode.AllDirF.root_cost root);
            info_user_message "Size of Tree: %f" (tree_size new_tree)
        end;
        let res = 
            single_characters_cost +. not_single_character_cost +. 
            (AllDirNode.AllDirF.root_cost root)
        in
        res

    let check_cost_all_handles ptree =
        All_sets.Integers.fold 
            (fun handle cost -> (check_cost ptree handle None) +. cost)
            (Ptree.get_handles ptree)
            0.0

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
            if debug_single_assignment then
            info_user_message "assign signle subtree on node %d,parent=%d" 
            current parent;
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
                | _   -> failwith "more than one root? AllDirChar.assign_single_handle 2"
            in
            let generate_root_and_assign_it rootg edge ptree =
                let a, b =
                    match edge with
                    | `Edge x ->  x
                    | `Single a -> a, a
                in
                let root, rooth = get_root_direction rootg in
                let handle_node = 
                    (Ptree.get_node_data a ptree).AllDirNode.unadjusted 
                        --> AllDirNode.not_with b
                        --> (fun x -> AllDirNode.force_val x.AllDirNode.lazy_node)
                and other_node = 
                    (AllDirNode.not_with a (Ptree.get_node_data b ptree).AllDirNode.unadjusted)
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
                    if debug_uppass_fn then 
                          Printf.printf "root_median is (%d,%d)\n%!" a b; 
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
        res

    let unadjust ptree = ptree

    let refresh_all_edges adjusted root_opt do_roots start_edge_opt ptree =
        let refresh_edge rhandle root_opt ((Tree.Edge (a,b)) as e) (acc,ptree) =
            if debug_uppass_fn then
                info_user_message "Refreshing %d--%d as %s" a b 
                                (if rhandle then "a root edge" else "an edge")
            else ();
            let data,ptree = 
                if rhandle then
                    let p1 = (create_partition ptree b a,b)
                    and p2 = (create_partition ptree a b,a) in
                    match hashdoublefind ptree [p1;p2] with
                    | Some x -> create_lazy_edge ~branches:x rhandle root_opt adjusted ptree a b 
                    | None -> create_lazy_edge rhandle root_opt adjusted ptree a b 
                else
                    create_lazy_edge rhandle root_opt adjusted ptree a b 
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
            info_user_message "Performing Calculation on Root" else ();
        let ptree = 
            match start_edge_opt with
            | Some (a,b) ->
          (*      info_user_message "start_edge_opt=%d,%d,call refresh_edge" a
          *      b; *)
                let _,t = refresh_edge true root_opt
                            (Tree.Edge (a,b)) 
                            (Tree.EdgeMap.empty,ptree) in
               (* info_user_message "end of refresh_edge"; *)
                t
            | None ->
               (* info_user_message "no start_edge_opt"; *)
                All_sets.Integers.fold 
                    (fun h ptree ->
                        try begin
                        match (Ptree.get_component_root h ptree).Ptree.root_median with
                            | Some ((`Edge (a,b)),c) -> 
                                let _,t = refresh_edge true root_opt
                                            (Tree.Edge (a,b)) 
                                            (Tree.EdgeMap.empty,ptree) in
                                t
                            | None
                            | Some _ -> ptree
                        end with | Not_found -> 
                            begin match Ptree.get_node h ptree with
                                | Tree.Leaf (a,b)
                                | Tree.Interior (a,b,_,_) -> 
                                    let _,t = refresh_edge true root_opt
                                                (Tree.Edge (a,b)) 
                                                (Tree.EdgeMap.empty,ptree) in
                                    t
                                | Tree.Single _ -> ptree
                            end)
                    ptree.Ptree.tree.Tree.handles
                    ptree
        in
        (* perform uppass heuristic --fill all directions *)
        current_snapshot "AllDirChar refresh_all_edges uppass heuristic";
        if debug_uppass_fn then
            info_user_message "Performing Uppass Heurisitic"
        else ();
        let ptree =
            match start_edge_opt with
            | Some (a,b) ->
            (*    info_user_message "start_edge_opt=%d,%d,call pre_order_node.."
            *    a b; *)
                let res =   Tree.pre_order_node_with_edge_visit_simple_root
                            add_vertex_pre_order
                            (Tree.Edge (a,b))
                            ptree.Ptree.tree ptree
                in
           (*     info_user_message "End of pre_order_node..."; *)
                res
            | None ->
               (*     info_user_message "no start_edge_opt "; *)
                    All_sets.Integers.fold
                        (fun h ptree ->
                                try 
                                    begin match (Ptree.get_component_root h ptree).Ptree.root_median with
                                        | Some ((`Edge (a,b)),c) -> 
                                            Tree.pre_order_node_with_edge_visit_simple_root
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
                                            Tree.pre_order_node_with_edge_visit_simple_root
                                                    add_vertex_pre_order
                                                    (Tree.Edge (a,b))
                                                    ptree.Ptree.tree ptree
                                        | Tree.Single _ -> ptree
                                    end
                        )
                        ptree.Ptree.tree.Tree.handles
                        ptree
        in
        (* fill in roots for all edges *)
        current_snapshot "AllDirChar refresh_all_edges internal fold";
        if do_roots then begin
            if debug_uppass_fn then
                info_user_message "Performing Refresh on all edges"
            else ();
            let new_edges,ptree =
                Tree.EdgeSet.fold
                        (refresh_edge false None)
                        ptree.Ptree.tree.Tree.d_edges
                        (Tree.EdgeMap.empty,ptree)
            in 
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
                    (* although this section is large, it only chooses the
                     * starting edge for the starting edge to traverse *)
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
                printf "\n%!";
                Tree.Continue, ()
            | Tree.Interior (x,p,_,_) ->
                printf "Node %d -- %d:" x p;
                AllDirNode.AllDirF.dump_node (f) (Ptree.get_node_data x ptree)
                                                 (Ptree.get_node_data p ptree);
                printf "\n%!";
                Tree.Continue, ()
        in
        Ptree.post_order_node_visit traversal x ptree ()

    let add_component_root ptree handle root = 
        { ptree with 
        Ptree.component_root = All_sets.IntegerMap.add handle root
        ptree.Ptree.component_root }

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
        (tree,inc)

    (* ------------------------------------------------------------------------ *)
    (* Now we define a function that can adjust all the vertices in the tree
    * to improve the overall cost of the tree, using only the
    * [AllDirNode.adjusted] field of each. *)
    type adjust_acc = bool * All_sets.Integers.t option All_sets.IntegerMap.t * phylogeny
    let adjust_tree max_count nodes ptree =
        let mode = match !Methods.cost with
            | `Iterative x -> x
            | _ -> assert false
        and max_count = match max_count with
            | Some x -> x
            | None -> max_int
        in
        let all_edges =
            Tree.EdgeMap.fold
                (* change the above to tuples and this as well *)
                (fun (Tree.Edge (a,b)) _ acc -> (a,b)::acc)
                ptree.Ptree.edge_data
                []
        in
    (* We start by defining a function to adjust one node *)
        let adjust_node chars_to_check ch1_k ch2_k parent_k mine_k ptree =
            current_snapshot 
                (Printf.sprintf "AllDirChar.adjust_node %d" mine_k);
            if debug_adjust_fn then
                info_user_message "AllDirCahr.adjust_node, on mine=%d with c1=%d,c2=%d p=%d" 
                                        mine_k ch1_k ch2_k parent_k;
            let gnd x = Ptree.get_node_data x ptree in
            let mine,modified =
                AllDirNode.AllDirF.readjust mode chars_to_check (gnd ch1_k)
                                            (gnd ch2_k) (gnd parent_k) (gnd mine_k)
            in
            if All_sets.Integers.is_empty modified
                then modified,[],ptree
                else begin
                    let ptree = Ptree.add_node_data mine_k mine ptree in
   (*                 Printf.printf "\n 4 Try to test refresh all edges..\n%!";
                    refresh_all_edges true None true None ptree;
                    Printf.printf "\n End of refresh all edges... \n\n%!";
   *)
                    modified, [ch1_k;ch2_k;mine_k;parent_k], ptree
                end
    (* adjust root --for likelihood *)
        and adjust_root (changed,affected,ptree) c2c handle a b : adjust_acc = 
            if debug_adjust_fn then
                    info_user_message "Adjusting root with %d,%d then None" a b;
            (* find edge root, create otherwise, and iterate *)
            let new_root = 
                let a_nd = Ptree.get_node_data a ptree 
                and b_nd = Ptree.get_node_data b ptree in
                let o_nd =
                    try
                        let one = Ptree.get_edge_data (Tree.Edge (a,b)) ptree in
                        let tmp = [{ AllDirNode.lazy_node = one;
                                                      dir = Some (a,b);
                                                     code = ~-1; }]
                        in
                        {AllDirNode.adjusted=tmp;AllDirNode.unadjusted=tmp}
                    with | Not_found ->
                        AllDirNode.AllDirF.median None None a_nd b_nd
                in
                AllDirNode.AllDirF.edge_iterator None o_nd a_nd b_nd
            in 
            (* below we apply the new branch length data to the left and right
             * and add the result to the ptree. n_root is striped of direction. *)
            let n_root,ptree = 
                create_lazy_edge true (Some new_root) false ptree a b
            in
            let n_root =
                let tmp = [{ AllDirNode.lazy_node = n_root;
                                              dir = Some (a,b);
                                             code = ~-1; }]
                in
                {AllDirNode.adjusted=tmp;AllDirNode.unadjusted=tmp}
            in
            (* assign the root and cost *)
            let ptree = 
                let check_cost = check_cost ptree handle (Some (n_root,`Edge (a,b))) in
                Ptree.assign_root_to_connected_component 
                    handle
                    (Some (`Edge (a,b),n_root))
                    check_cost
                    None
                    ptree
            in
            let ptree =
                refresh_all_edges true (Some n_root) true (Some (a,b)) ptree
            in
            (changed,affected,ptree)
    (* add modified vertices in node_list to the set *)
        and add_vertices_affected node_list codes affected = 
            let add_one code affected =
                if All_sets.IntegerMap.mem code affected then
                    match All_sets.IntegerMap.find code affected with
                    | Some conts ->
                            let res = All_sets.Integers.union conts codes in
                            All_sets.IntegerMap.add code (Some res) affected
                    | None -> (*assert false*) affected
                else All_sets.IntegerMap.add code (Some codes) affected
            in
            List.fold_right (fun x acc -> add_one x acc) node_list affected
        in
    (* compose the above functions to adjust and modify the affected nodes *)
        let adjust_vertices_affected ((modified,affected_nodes,ptree) as acc) c2c prev curr =
            if not (All_sets.IntegerMap.mem curr c2c) then acc
            else match Ptree.get_node curr ptree with 
                | (Tree.Interior (c,p,c1,c2)) as nd ->
                                    (* it was modified before, so iterate it *)
                    let c2c = All_sets.IntegerMap.find curr c2c in
                    let a,b = Tree.other_two_nbrs prev nd in
                    (* modified codes, affected nodes, tree *)
                    let ccodes,affected,n_ptree = adjust_node c2c a b prev curr ptree in
                    let new_affected = add_vertices_affected affected ccodes affected_nodes in
                    let modified = (0 != (List.length affected)) || modified in
                    (modified, new_affected, n_ptree)
                | Tree.Leaf _ 
                | Tree.Single _ ->  acc
        in
    (* loop to adjust a tree and *)
        let adjust_until_nothing_changes max_count start_ptree = 
            (* empty set *)
            let first_affected = 
                match nodes with
                | None -> 
                    All_sets.IntegerMap.map (fun _ -> None) ptree.Ptree.node_data
                | Some items -> items
            in
            (* Post order traversal of internal nodes *)
            let adjust_loop prev_affected handle adjust_acc =
                match (Ptree.get_component_root handle ptree).Ptree.root_median with
                | Some ((`Edge(a,b)),c) -> 
                        let start_edge = Tree.Edge (a,b) in
                        Tree.post_order_node_with_edge_visit_simple (* f e ptree acc *)
                        (fun prev curr acc ->
                           adjust_vertices_affected acc prev_affected prev curr)
                        start_edge
                        ptree.Ptree.tree
                        adjust_acc
                | Some ((`Single a), rootg) ->  false,prev_affected,ptree
                | None -> false,prev_affected,ptree
            (* loop for rerooting and applying iterative on the resultant path *)
            and adjust_reroot_loop affected (modified,aff_n,ptree) (a,b) =
                (* a simple reroot, since the reroot_fn requires incremental as
                 * a return type, but nodes are integers, this can be changed by
                 * maping over set to convert, but this is here to quickly get
                 * things going --it is the same functionality as DO with force *)
                let simple_reroot edge ptree = 
                    let Tree.Edge (h, n) = edge in
                    let ol_handle = Ptree.handle_of h ptree in
                    let ptree,inc =
                        ptree --> Ptree.remove_root_of_component ol_handle 
                              --> Ptree.move_handle h
                    in
                    let ptree =
                        ptree --> create_root h n
                              --> add_component_root ptree h
                    in
                    (ptree,inc,h)
                in
                (* main portion of reroot -- iterate steps *)
                let ptree,incs,handle = simple_reroot (Tree.Edge (a,b)) ptree in
                adjust_root (modified,aff_n,ptree) affected handle a b
            in
            (* recursive loop of for changes *)
            let rec iterator count prev_cost affected ptree =
                let (changed,new_affected,new_ptree : adjust_acc) = 
                    let none_affected = All_sets.IntegerMap.empty in
                    (* perform on each tree *)
                    if using_likelihood ptree
                        then begin
                            List.fold_left
                                (adjust_reroot_loop affected)
                                (true,none_affected,ptree)
                                (all_edges)
                        end else
                            All_sets.Integers.fold
                                (adjust_loop affected)
                                (ptree.Ptree.tree.Tree.handles)
                                (true,none_affected,ptree)
                in
                (* now ptree can be used normaliy *)
                let new_cost = check_cost_all_handles new_ptree in
                if debug_adjust_fn then
                    info_user_message "Iteration %d completed: %f --> %f (%b)" 
                                      (max_count - count) prev_cost new_cost changed;
                if (not changed) || (count = 1) || (prev_cost =. new_cost) || (new_cost > prev_cost)
                    then ptree
                    else iterator (count - 1) new_cost new_affected new_ptree
            in
            let initial_cost = check_cost_all_handles ptree in
            iterator max_count initial_cost 
                               first_affected
                               ptree
        in
        let set_handle_n_root_n_cost handle ptree =
            if using_likelihood ptree then ptree 
            else begin
                let comp_root = Ptree.get_component_root handle ptree in
                match comp_root.Ptree.root_median with
                | None -> failwith "Huh? AllDirChar.set_handle_n_root"
                | Some ((`Edge (a, b)) as edge, root) ->
                    let sets = get_active_ref_code ptree
                    and ad = Ptree.get_node_data a ptree
                    and bd = Ptree.get_node_data b ptree in
                    let root = 
                        let n = AllDirNode.AllDirF.median None None ad bd in
                        AllDirNode.AllDirF.to_single (Some n) None ad None bd sets
                    in
                    ptree
                        --> Ptree.assign_root_to_connected_component handle
                                (Some (edge, root))
                                (check_cost ptree handle None)
                                None
                        --> refresh_all_edges true None true None
                | Some _ -> ptree
            end
        in
        let newtree = adjust_until_nothing_changes max_count ptree in   
        let ptree = All_sets.Integers.fold (set_handle_n_root_n_cost)
                               (ptree.Ptree.tree.Tree.handles)
                               (newtree)
        in
        ptree
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
    let internal_downpass do_roots (ptree : phylogeny) : phylogeny =
        (* function to process tree->node->charactername to int->float hashtbl
         * for all the node ids passed --(multiple node_id capability for uppass) *)
        let ptree = 
            (* A function to add  the vertices using a post order traversal 
             * from the Ptree library. *)
            let add_vertex_post_order prev code ptree =
                current_snapshot "AllDirChar.internal_downpass.add_vertex_post_order";
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
                                                code a b prev;
                        let interior = 
                            let p1 = create_partition ptree b code,b
                            and p2 = create_partition ptree a code,a in
                            match hashdoublefind ptree [p1;p2] with
                            | Some x -> create_lazy_interior_down ~branches:x ptree (Some code) a b
                            | None   -> create_lazy_interior_down ptree (Some code) a b
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
                        | None 
                        | Some _ -> ptree
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
        let ptree = refresh_all_edges false None true None ptree in
        if do_roots then refresh_roots ptree else ptree

    let clear_internals t = 
        {t with Ptree.data = Data.remove_bl t.Ptree.data; }

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


    let general_pick_best_root selection_method ptree =
        (* debug msg 
        Printf.printf "\n #2. general_pick_best_root\n%!";
        let newcost = check_cost_all_handles ptree in
         debug msg*)
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

    let pick_best_root ptree =
        if using_likelihood ptree then ptree
        else general_pick_best_root blindly_trust_downpass ptree

    (* ----------------- *)
    (* function to adjust the likelihood model of a tree using BFGS --quasi
     * newtons method. Function requires three directions. *)
    let model_fn tree = 
        (* replace nodes in a tree, copying relevent data structures *)
        let substitute_nodes nodes tree =
            let adder acc x = All_sets.IntegerMap.add (AllDirNode.AllDirF.taxon_code x) x acc in
            let node_data = List.fold_left adder All_sets.IntegerMap.empty nodes in
            internal_downpass true {tree with Ptree.node_data = node_data}
        (* get all characters to iterate *)
        and chars = 
            let chars = `Some (Data.get_chars_codes_comp tree.Ptree.data `All) in
            Data.get_code_from_characters_restricted `AllStatic tree.Ptree.data chars
        (* get current model *)
        and current_model data chars = 
            let get_model x = match Hashtbl.find data.Data.character_specs x with
                | Data.Static dat ->
                    begin match dat.Parser.SC.st_type with
                        | Parser.SC.STLikelihood model -> model
                        | _ -> failwith "unsupported static character"
                    end
                | _ -> failwith "unsupported character type"
            in
            match List.map get_model chars with
            | h :: t ->
                if List.fold_left (fun acc x -> acc && (x = h)) true t then
                    h
                else failwith "Inconsistent Model over characters"
            | _ -> failwith "No Characters found"
        in
        (* function for processing a model and applying to a tree --inner loop *)
        let f_likelihood f tree chars current_model new_values =
            let ntree =
                Parser.SC.STLikelihood (f current_model new_values)
                    --> Data.apply_on_static_chars tree.Ptree.data chars
                    --> AllDirNode.AllDirF.load_data
                    --> (fun (x,y) -> substitute_nodes y {tree with Ptree.data = x} )
            in
            let ncost = Ptree.get_cost `Adjusted ntree in
            ntree,ncost
        in
        (* iterate the model only if using likelihood *)
        if using_likelihood tree then begin
            let current_model = current_model tree.Ptree.data chars
            and current_cost = Ptree.get_cost `Adjusted tree in
            let current_array =
                match MlModel.get_current_parameters_for_model current_model with
                | Some x -> x | None -> failwith "No parameters to modify" in 
            let params, (best_tree, best_cost) = 
                match MlModel.get_update_function_for_model current_model with
                | Some func -> 
                    let tree = update_branches tree in
                    MlModel.bfgs_method (f_likelihood func tree chars current_model)
                                        current_array 
                                        (tree,current_cost)
                | None -> failwith "No function to optimize with"
            in
            if best_cost < current_cost then best_tree else tree
        end else begin
            tree
        end

    let adjust_fn ?(epsilon=1.0e-4) ?(max_iter=20) iterations tree = 
        if using_likelihood tree then begin
            let rec loop_ iter curr tree =
                if iter = max_iter then tree
                else begin
                    let mtree = model_fn tree in
                    let mcost = Ptree.get_cost `Adjusted mtree in
                    if (abs_float (curr -. mcost)) <= epsilon then mtree
                    else begin
                        let btree = adjust_tree iterations None mtree in
                        let bcost = Ptree.get_cost `Adjusted btree in
                        if (abs_float (mcost -. bcost)) <= epsilon then btree
                        else loop_ (iter+1) bcost btree
                    end
                end
            in
            (* ensures that we modify the branch lengths *)
            let first = adjust_tree iterations None tree in
            loop_ 0 (Ptree.get_cost `Adjusted first) first
        end else begin
            adjust_tree iterations None tree
        end

    (* ---------- *)
    let downpass ptree =
        if debug_downpass_fn then Printf.printf "DownPass Begins\n%!";
        current_snapshot "AllDirChar.downpass a";
        let res = 
            match !Methods.cost with
            | `Exhaustive_Strong
            | `Exhaustive_Weak
            | `Normal_plus_Vitamines
            | `Normal -> internal_downpass true ptree
            | `Iterative (`ApproxD iterations)
            | `Iterative (`ThreeD  iterations) ->
                ptree
                    --> clear_internals (* remove BLs *)
                    --> internal_downpass true
                    --> pick_best_root
                    --> assign_single true
                    --> adjust_fn iterations
        in
        current_snapshot "AllDirChar.downpass b";
        if debug_downpass_fn then Printf.printf "Downpass Ends\n%!";
        res

    let uppass ptree = 
        if debug_uppass_fn then Printf.printf "UPPASS begin: \n%!";
        let tree = match !Methods.cost with
            | `Exhaustive_Strong
            | `Exhaustive_Weak
            | `Normal_plus_Vitamines
            | `Normal -> 
                ptree --> pick_best_root --> assign_single true
            | `Iterative (`ApproxD _)
            | `Iterative (`ThreeD _) -> ptree
        in
        if debug_uppass_fn then Printf.printf "UPPASS ends. \n%!";
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
    | None ->
        refresh_all_edges true None false (Some (a,b)) ptree
    | Some edge ->
        if not (using_likelihood ptree) then
            refresh_all_edges true None false (Some (a,b)) ptree
        else begin
            let edge = (* use this as edge data for likelihood *)
                let single = [AllDirNode.with_both a b edge.AllDirNode.adjusted] in
                { AllDirNode.unadjusted = single; adjusted = single }
            in
            refresh_all_edges true (Some edge) false (Some (a,b)) ptree
        end

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

    (* break_fn has type handle * int (node) -> tree -> tree * delta * aux_data *)
    let break_fn (tree_node, clade_node_id) (ptree : phylogeny) =
        (* remove branch length table *)
        let () = 
            if Hashtbl.length (ptree.Ptree.data.Data.branches) > 0 
                then begin
                    Hashtbl.clear (ptree.Ptree.data.Data.branches);
                    info_user_message "Clearing supplied branch lengths."
                end else ()
        in
        (* -------------------------- *)
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
                    --> update_break_delta left_delta
                    --> update_break_delta right_delta
                    --> refresh_all_edges false None true None
                    --> refresh_roots
                    --> uppass
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
                abs_float 
                ( (prev_cost -. (new_cost -. (rc +. ptree.Ptree.origin_cost))) -.  tc)
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

    let break_fn ((s1, s2) as a) b =
        let res = match !Methods.cost with
        | `Normal ->
            b   --> clear_internals
                --> break_fn a
        | `Iterative (`ApproxD _)
        | `Iterative (`ThreeD _)
        | `Exhaustive_Weak
        | `Normal_plus_Vitamines ->
                let breakage = break_fn a (clear_internals b) in
                let nt =
                    refresh_all_edges true None true None
                                           (breakage.Ptree.ptree)
                in
                { breakage with 
                    Ptree.ptree = nt; }
        | `Exhaustive_Strong ->
                let breakage = break_fn a (clear_internals b) in
                let nt = 
                    refresh_all_edges true None true None
                                           (breakage.Ptree.ptree) in
                { breakage with 
                    Ptree.ptree = nt;
                    incremental = []; 
                    break_delta = (Ptree.get_cost `Adjusted b) -. 
                                  (Ptree.get_cost `Adjusted nt);
                }
        in
        res

    let equal_float = (* Up to three significant positions *)
        let positions = 3. in
        let factor = 10. ** positions in
        fun a b ->
            let truncate x = truncate (x *. factor) in
            (truncate a) = (truncate b)

    (* ----------------- *)
    (* join_fn must have type join_1_jxn -> join_2_jxn -> delta -> tree -> tree *)
    let join_fn _ jxn1 jxn2 ptree =
        if debug_join_fn then
            info_user_message "Time to join! (%s) and (%s)"
                (match jxn1 with
                    | Tree.Single_Jxn x -> string_of_int x
                    | Tree.Edge_Jxn (x,y) -> (string_of_int x) ^","^ (string_of_int y))
                (match jxn2 with
                    | Tree.Single_Jxn x -> string_of_int x
                    | Tree.Edge_Jxn (x,y) -> (string_of_int x) ^","^ (string_of_int y))
        else ();
        let lift_data edge_l edge_r i_code ptree = 
            if not (using_likelihood ptree) then begin
                let node = AllDirNode.AllDirF.median (Some i_code) None
                            (Ptree.get_node_data edge_l ptree)
                            (Ptree.get_node_data edge_r ptree) in
                Ptree.add_node_data i_code node ptree
            end else begin
                try let lr = Ptree.get_edge_data (Tree.Edge (edge_l,edge_r)) ptree in
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
                with | Not_found ->
                    failwithf "Cannot lift %d -- %d to %d" edge_l edge_r i_code
            end
        in
        let ret, ((td1,td2,rrd) as tree_delta) = 
            Tree.join jxn1 jxn2 ptree.Ptree.tree 
        in
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
            add_component_root ptree handle (create_root v h ptree)
        in
(*
        assert (
            let ptree, _ = reroot_fn true (Tree.Edge (v, h)) ptree in
            let cost = Ptree.get_cost `Unadjusted ptree in
            let size = tree_size ptree in
            Printf.printf "REDIAGNOSING THE TREE:\n\n%!";
            let ptree, _ = ptree --> downpass --> uppass -->
                           reroot_fn true (Tree.Edge (v, h)) in
            let res = equal_float cost (Ptree.get_cost `Unadjusted ptree) in
            if not res then 
                Printf.printf ("The old cost: %f\t new cost: %f\n%!"^^ 
                               "The old size: %f\t new size: %f\n%!")
                    cost (Ptree.get_cost `Unadjusted ptree)
                    size (tree_size ptree);
            res);
*)
        ptree, tree_delta

    let get_one side = 
        match side with
        | `Single (x, _) | `Edge (x, _, _, _) -> x

    let join_fn a b c d =
        let ptree, tdel = match !Methods.cost with
            | `Normal -> 
                d --> clear_internals --> join_fn a b c 
            | `Iterative (`ThreeD iterations)
            | `Iterative (`ApproxD iterations) ->
                let tree, delta = join_fn a b c (clear_internals d) in
                let tree = 
                   tree --> pick_best_root
                        --> assign_single true 
                        --> adjust_fn iterations
                in
                tree, delta
            | `Normal_plus_Vitamines
            | `Exhaustive_Weak
            | `Exhaustive_Strong ->
                let tree, delta = join_fn a b c (clear_internals d) in
                uppass tree, delta 
        in
        ptree, tdel

    type tmp = Edge of (int * int) | Clade of a 

    let cost_fn jxn1 jxn2 _ (*delta*) clade_data (tree : phylogeny) =
        let rec forcer edge =
            match edge with
            | Edge (a, b) ->
                    AllDirNode.force_val (Ptree.get_edge_data (Tree.Edge (a, b)) tree)
            | Clade x -> 
                    (match x.AllDirNode.unadjusted with (* leaf *)
                    | [x] -> force_node x
                    | _ -> failwith "AllDirChar.cost_fn")
        in
        let clade_data = 
            match !Methods.cost with
            | `Iterative (`ThreeD _) ->
                    (match jxn2 with
                    | Tree.Single_Jxn h -> (* forcer (Clade clade_data) *)
                        forcer (Clade (Ptree.get_node_data (Tree.int_of_id h) tree))
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
                    | Ptree.Cost x -> Ptree.Cost (abs_float (0.85 *. x))
                    | x -> x)
            | `Iterative `ThreeD _
            | `Exhaustive_Weak
            | `Normal_plus_Vitamines
            | `Normal ->
                    (match cost_fn a b c d e with 
                          | Ptree.Cost x -> Ptree.Cost (abs_float x)
                          | x -> x)
            | `Exhaustive_Strong ->
                    let pc = Ptree.get_cost `Adjusted e in
                    let (nt, _) = join_fn [] a b e in
                    Ptree.Cost (abs_float (((Ptree.get_cost `Adjusted nt) -.  pc)))

 
    let root_costs tree = 
        let collect_edge_data edge node acc =
            let cost = AllDirNode.OneDirF.tree_cost None node in
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
     * Creates a hashtable with all the branch data. The key is the pair of
     * branch lengths and the value is either a single length or a list of
     * branch lengths in the case of multiple character sets. *)
    let branch_table ptree =
        let trees_table = Hashtbl.create 13 in
        let create_branch_table handle _ = 
            let single_node prev curr =
                let pair = (min curr prev, max curr prev) in
                let dat = AllDirNode.AllDirF.get_times_between 
                            (Ptree.get_node_data curr ptree)
                            (Ptree.get_node_data prev ptree) in
                let name_it x = match dat with | [_] -> `Single x 
                                               | []  -> failwith "No character Sets"
                                               | _   -> `Name
                in
                let () = List.iter
                    (fun (code,length) -> match length with
                        | Some length -> Hashtbl.add trees_table pair (name_it length)
                        | None -> ()
                    ) dat in
                ()
            in
            let a,b =
                match (Ptree.get_component_root handle ptree).Ptree.root_median with
                | Some ((`Edge (a,b)),_) -> a,b
                | None | Some _ -> raise Not_found
            in
            let (),() =
                Ptree.post_order_node_with_edge_visit
                    (fun prev curr _ -> single_node prev curr)
                    (fun prev curr _ _ -> single_node prev curr)
                    (Tree.Edge (a,b))
                    ptree
                    ()
            in
            ()
        in
        let () = All_sets.Integers.fold
                    (create_branch_table)
                    (Ptree.get_handles ptree)
                    ()
        in
        trees_table


    let to_formatter (atr : Xml.attributes)  
            (tree : (a, b) Ptree.p_tree) : Xml.xml =
        let data = tree.Ptree.data in
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
        let merger a b root = (`Set [`Single root; `Single a; `Single b]) 
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
                     (*   Printf.printf "root median at (%d,%d) : %!" a b;*)
                        let recost = 
                            let root = get_unadjusted (-1) root in
                            (Node.cmp_subtree_recost root) +. recost 
                        in
                       (* Printf.printf "recost = %f \n%!" recost;*)
                        (* We override the root now to continue using the single
                        * assignment of the handle *)
                        let sroot, sa = 
                            let a = Ptree.get_node_data a tree in
                            let s = get_single b a in
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
        (*Printf.printf "ALLDIRCHAR.TO_FORMATTER: recost = %f, cost = %f \n%!"
        recost cost; *)
        (RXML -[Xml.Trees.forest] 
            ([Xml.Trees.recost] = [`Float recost])
            ([Xml.Trees.cost] = [`Float cost])
            ([atr])
            { set trees } --)

end
