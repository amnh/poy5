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

let () = SadmanOutput.register "AllDirChar" "$Revision: 3349 $"

module IntSet = All_sets.Integers
module IntMap = All_sets.IntegerMap
module IntSetMap = All_sets.IntSetMap

let debug_profile_memory    = false
let debug_node_fn           = false
let debug_model_fn          = false
let debug_adjust_fn         = false
let debug_join_fn           = false
let debug_break_fn          = false
let debug_branch_fn         = false
let debug_cost_fn           = false
let debug_uppass_fn         = false
let debug_downpass_fn       = false
let debug_single_assignment = false
let debug_diagnosis         = false
let debug_reroot            = false
let debug_create_root       = false 
let debug_assign_single_handle = false

let current_snapshot x = 
    if debug_profile_memory then MemProfiler.current_snapshot x

let update_node_manager ptree f d : unit = match d with
    | Some node_mgr -> node_mgr#update_iterate ptree f
    | None          -> ()

let (-->) a b = b a

open Numerical.FPInfix      (* open up fuzzy floating point comparison *)

let error_user_message format =
    Printf.ksprintf (Status.user_message Status.Error) format
let info_user_message format =
    Printf.ksprintf (Status.user_message Status.Information) format
let warning_user_message format =
    Printf.ksprintf (Status.user_message Status.Warning) format
let failwithf format =
    Printf.ksprintf (failwith) format

module F : Ptree.Tree_Operations
    with type a = AllDirNode.AllDirF.n
        with type b = AllDirNode.OneDirF.n = struct

    type a = AllDirNode.AllDirF.n
    type b = AllDirNode.OneDirF.n
    type phylogeny = (a, b) Ptree.p_tree

    (** Type used to define an accumulator to adjust the nodes/edges of a tree.
        This stores a set of modified nodes, the tree, and if things were
        updated in the adjustment function. *)
    type adjust_acc = bool * IntSet.t option IntMap.t * phylogeny

    let force_node x = AllDirNode.force_val x.AllDirNode.lazy_node

    (* process tree data to find branch lengths *)
    let hashdoublefind tree partitions : (int,(int,float) Hashtbl.t) Hashtbl.t option =
        (* converts table with character names to character ids *)
        let transform_keys to_table from_table char_name_tbl =
            Hashtbl.iter 
                (fun name length ->
                    try Hashtbl.add to_table (Hashtbl.find char_name_tbl name) length
                    with | Not_found ->
                        if debug_branch_fn then
                            error_user_message "Couldn't find character name %s" name)
                from_table;
            to_table
        (* retrieve the name of the tree, or no name *)
        and t_name = match (tree.Ptree.tree).Tree.tree_name with
            | Some tree_n -> String.uppercase tree_n | None -> ""
        in
        (* test if that table exists, and convert each partition to the node id
            it belongs from, and the charcter_names to character ids *) 
        try match tree.Ptree.data.Data.branches with
            | Some btable ->
                let ret_table = Hashtbl.create 27 in
                let tree_tbl = Hashtbl.find btable t_name in
                let res = List.fold_left
                    (fun acc (partition,node_id) ->
                        try let node_n = All_sets.IntSetMap.find partition tree_tbl in
                            let tbl = 
                                transform_keys (Hashtbl.create 27) node_n
                                               (tree.Ptree.data).Data.character_names
                            in
                            let () = Hashtbl.replace ret_table node_id tbl in
                            true
                        with | Not_found -> (false or acc))
                    false
                    partitions
                in
                if res then Some ret_table else None
            | None -> None
        (* return nothing if the node wasn't found *)
        with | Not_found ->
            if debug_branch_fn then
                error_user_message "Couldn't find tree name %s" t_name;
            None


    (* check if the ptree has likelihood characters by it's data; obviously this
     * means that data must be up to date --an invariant that we hold. The old
     * version of this function, that checks the ptree is left below for
     * posterity. *)
    let rec using_likelihood types ptree =
        let data_test = match types with
            | `Static   -> Data.has_static_likelihood ptree.Ptree.data
            | `Dynamic  ->
                begin match Data.type_of_dynamic_likelihood ptree.Ptree.data with
                    | Some _ -> true
                    | None   -> false
                end
            | `OnlyStatic ->
                (using_likelihood `Static ptree ) && not (using_likelihood `Dynamic ptree)
            | `Either   ->
                (using_likelihood `Static ptree ) || (using_likelihood `Dynamic ptree)
        in
        data_test


    (** [create_branch_table table ptree] 
     * Creates a hashtable with all the branch data. The key is the pair of
     * nodes lengths and the value is either a single length or a list of
     * branch lengths in the case of multiple character sets. *)
    let branch_table (opts : [`Single | `Final | `Max ] option) ptree =
        let trees_table = Hashtbl.create 13 in
        let inc_parsimony = match opts with
            | Some _ -> (true, opts)
            | None   -> (false,None)
        in
        let adjusted = match opts with
            | Some (`Final | `Max) | None -> false
            | _ when using_likelihood `Either ptree -> false
            | Some `Single -> true
        in
        let create_branch_table handle () =
            let rec single_node prev curr =
                let pair = (min curr prev, max curr prev) in
                let dat =
                    AllDirNode.AllDirF.get_times_between ~adjusted ~inc_parsimony
                            (Ptree.get_node_data curr ptree)
                            (Some (Ptree.get_node_data prev ptree))
                in
                let name_it x = match dat with
                    | [_] -> `Single x
                    | []  -> failwith "No character Sets"
                    | xs  -> `Name xs
                in
                List.iter
                    (fun (code,length) -> match length with
                        | Some length -> Hashtbl.add trees_table pair (name_it length)
                        | None -> ())
                    dat
            and traversal a b =
                Ptree.post_order_node_with_edge_visit
                    (fun prev curr _ -> single_node prev curr)
                    (fun prev curr _ _ -> single_node prev curr)
                    (Tree.Edge (a,b))
                    ptree
                    ()
            in
            let (),() =
                try match (Ptree.get_component_root handle ptree).Ptree.root_median with
                    | Some ((`Edge (a,b)),_) -> traversal a b
                    | None | Some _ -> raise Not_found
                with | Not_found ->
                    begin match Ptree.get_node handle ptree with
                        | Tree.Leaf (a,b)
                        | Tree.Interior (a,b,_,_) -> traversal a b
                        | Tree.Single _ -> (),()
                    end
            in
            ()
        in
        IntSet.fold (create_branch_table) (Ptree.get_handles ptree) ();
        trees_table


    (* Update Data.d in ptree with branch data. Used to transfer data between
     * dynamic and static likelihood; used under Dynamic Likelihood only. *)
    let update_branches ptree =
        let get_codestable data node1 node2 =
            let codestimes =
                AllDirNode.AllDirF.get_times_between
                    (Ptree.get_node_data node1 ptree)
                    (Some (Ptree.get_node_data node2 ptree))
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
                    | None   -> ())
                codestimes;
            table
        in
        if using_likelihood `Either ptree then begin
            let name = match ptree.Ptree.tree.Tree.tree_name with
                | Some x -> String.uppercase x | None -> ""
            and treebranches = Hashtbl.create 1
            and insert_function setmap edge =
                let Tree.Edge (left,right) = edge in
                let p1,p2 = Ptree.create_partition ptree edge in
                let codestable = get_codestable ptree.Ptree.data left right in
                setmap --> IntSetMap.add p1 codestable
                       --> IntSetMap.add p2 codestable
            in
            let () =
                Ptree.get_edges_tree ptree
                    --> List.fold_left insert_function All_sets.IntSetMap.empty
                    --> Hashtbl.add treebranches name
            in
            {ptree with
                Ptree.data = {ptree.Ptree.data with Data.branches = Some treebranches;}}
        end else begin
            ptree
        end

(*  Creates a lazy edge which is the median between the data of the vertices
    with codes [a] and [b] in the tree [ptree]. Times in opposite directions
    SHOULD be set except when [root] is true. In this case, median_w_times is
    used and the tree is updated with the new directions calculated.
    [root_node] is used when the edge data is precomputed (during a break) *)
    let create_lazy_edge ?branches root root_node ptree a b = 
        if debug_node_fn then
            info_user_message 
                "Creating lazy edge between %d and %d with root %s"
                a b (if root then "true" else "false");
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
            AllDirNode.create_root_w_times dataa datab, ptree
        end


    (* Creates a valid vertex that only has the downpass information *)
    let create_lazy_interior_down ?branches ptree code a b =
        let () =
            if debug_node_fn then
                begin match code with
                | Some x ->
                    info_user_message
                        "Creating lazy interior down (%d) between %d and %d" x a b
                | None ->
                    info_user_message
                        "Creating lazy interior down (?) between %d and %d" a b;
                end
        in
        let a_nd = Ptree.get_node_data a ptree 
        and b_nd = Ptree.get_node_data b ptree in
        AllDirNode.AllDirF.median ?branches code None a_nd b_nd


    (* Creates a valid vertex that has the downpass and uppass information.
    * [a] and [b] have to be the component of the currently calculated subtree
    * that is still valid. *)
    let create_lazy_interior_up ptree code a b c =
        if debug_node_fn then
            info_user_message 
                "Creating lazy interior up of %d with c:%d, c:%d and p:%d" code a b c;
        (* interior nodes should have all this information *)
        let m_data = Ptree.get_node_data code ptree
        and a_data = Ptree.get_node_data a ptree
        and b_data = Ptree.get_node_data b ptree
        and c_data = Ptree.get_node_data c ptree in
        AllDirNode.AllDirF.uppass_heuristic c_data None m_data a_data b_data


    (* Create the edge data on (a,b) and use to update root information in tree *)
    let create_root a b (ptree : phylogeny) =
        let make_internal a b ptree =
            if debug_create_root then 
                Printf.printf "create root with edge (%d,%d)\n%!" a b;
            let norm = Tree.normalize_edge (Tree.Edge (a, b)) ptree.Ptree.tree in
            let lazy_edge_data = Ptree.get_edge_data norm ptree in
            let edge_data = AllDirNode.force_val lazy_edge_data in
            let adj_nodedata = {AllDirNode.lazy_node = lazy_edge_data; code = 0; dir = Some (a,b)} in
            let nodedata = { AllDirNode.unadjusted = [adj_nodedata]; adjusted = Some adj_nodedata} in
            let cost = Node.Standard.tree_cost None edge_data in
            if debug_create_root then begin
                Printf.printf "treecost=%f, root node data:\n%!" cost;
                AllDirNode.print_node_data nodedata true;
            end;
            {
                Ptree.root_median = Some ((`Edge (a, b)), nodedata);
                component_cost = cost;
                adjusted_component_cost = cost;
            }
        in
        match Ptree.get_node a ptree with
        | Tree.Leaf (_, x) ->
                assert(x = b);
                make_internal a b ptree
        | Tree.Interior (_, x, y, z) ->
                assert((x = b) || (y = b) || (z = b));
                make_internal a b ptree
        | Tree.Single _ ->
                let root = Ptree.get_node_data a ptree in
                let cost = AllDirNode.AllDirF.tree_cost None root in
                {
                    Ptree.root_median = Some ((`Single a), root);
                    component_cost = cost;
                    adjusted_component_cost = cost;
                }

    let prior_cost tree chars =
        let dyn_cost =
            if using_likelihood `Dynamic tree then
                All_sets.IntegerMap.fold
                    (fun k v min_prior ->
                        let new_prior = AllDirNode.AllDirF.min_prior chars v in
                        min new_prior min_prior)
                    (tree.Ptree.node_data)
                    (infinity)
            else
                0.0
        in
        dyn_cost +. (MlStaticCS.ncm_priors tree.Ptree.data chars)


    let total_cost ptree adj chars = Ptree.get_cost `Adjusted ptree


    (* calculate the size of a tree under likelihood; this is the sum of all the
     * branch lengths. The sum of each character set is calculated; this might
     * be better if each one is seperate(?). 0 is returned under other critera. *)
    let tree_size ptree chars =
        let has_code codes chars = match chars with
            | None    -> true
            | Some [] -> false
            | Some xs -> 
                let rec found i = 
                    if i = Array.length codes 
                        then false
                        else begin
                            if List.mem codes.(i) xs
                                then true
                                else found (i+1)
                        end
                in
                found 0
        in
        let edge_cost (Tree.Edge (a,b)) =
            let lst =
                try AllDirNode.AllDirF.get_times_between
                        (Ptree.get_node_data a ptree) (Some (Ptree.get_node_data b ptree))
                with | _ ->
                    AllDirNode.AllDirF.get_times_between
                        (Ptree.get_node_data b ptree) (Some (Ptree.get_node_data a ptree))
            in
            List.fold_left
                (fun acc x -> match x with 
                    | c,Some x when has_code c chars -> x +. acc 
                    | _,_ -> acc)
                0.0 lst
        in
        Tree.EdgeMap.fold
            (fun k _ acc -> acc +. (edge_cost k))
            ptree.Ptree.edge_data
            0.0

    (* Determine the cost of a tree from the handle. A optional root can be
     * passed if the tree requires it for applying the root. *)
    let check_cost new_tree handle new_root =
        let debug = false in
        if debug then Printf.printf "alldirchar.check_cost ->\n%!";
        (* Some characters are computed by the downpass, we extract that
         * portion of the cost, which is contained in the root.
         * not_single_character_cost holds the total cost of those characters *)
        let not_single_character_cost, root, root_edge =
            (* set initial value to subtract from -- cost of likelihood *)
            let tree_root, root_edge =
                match new_root with
                | Some a -> a
                | None ->
                    let x = Ptree.get_component_root handle new_tree in
                    match x.Ptree.root_median with
                    | Some (root_edge, root) -> root, root_edge
                    | None -> failwith "No root 1?"
            in
            match tree_root.AllDirNode.unadjusted with
            | [root] ->
                    let root = AllDirNode.force_val root.AllDirNode.lazy_node in
                    let cost =
                        List.fold_left
                            (fun acc y -> 
                                let add = Node.total_cost_of_type y root in
                                acc +. add 
                            )
                            (0.0)
                            (Node.not_to_single)
                    in
                    (*when there is cost between same states, we might have
                    * extra cost on the root because we 'create' the root node
                    * left child is state.1    right child is state.1
                    *                 |            |
                    *               root node take state1
                    * if the cost between same states is c, above tree give us
                    * cost = 2c. but in fact there was no root node, we create
                    * that to make things clear, the tree is just like this:
                    * left node is state.1 ---- right node is state.1  
                    * this give us cost = c only.*)
                    if debug then 
                        Printf.printf "not_single_character_cost, cost  = %f\n%!" cost;
                    cost, tree_root, root_edge
            | _ -> failwith "What?"
        in
        (* Other characters have their cost computed by adding up the costs of
           all of the branches. single_characters_cost is exactly that. *)
        let distance a b acc =
            let nodeA = Ptree.get_node_data a new_tree in
            let nodeB = Ptree.get_node_data b new_tree in
            let nda =
                let node = AllDirNode.get_adjusted_nodedata nodeA
                                            "AllDirChar.distance,no adjdata"
                in
                node.AllDirNode.lazy_node
            and ndb =
                let node = AllDirNode.get_adjusted_nodedata nodeB
                                            "AllDirChar.distance,no adjdata"
                in
                node.AllDirNode.lazy_node
            in
            (*AllDirNode.print_node_data nodeA true; 
            AllDirNode.print_node_data nodeB true; *)
            let dist =
                let p1 = fst (Ptree.create_partition new_tree (Tree.Edge (a,b))), a in
                match hashdoublefind new_tree [p1] with
                | Some _ as branches -> 
                    Node.distance_of_type ?branches (Node.has_to_single) (0.0)
                            (AllDirNode.force_val nda) (AllDirNode.force_val ndb)
                | None   ->
                    Node.distance_of_type (Node.has_to_single) (0.0)
                            (AllDirNode.force_val nda) (AllDirNode.force_val ndb)
            in
            if debug then
                Printf.printf "distance between node.%d and node.%d = %f(%f)\n%!"
                a b dist acc;
            dist +. acc
        in
        let single_characters_cost = match root_edge with
            | `Single _    -> 0.0
            | `Edge (a, b) ->
                if using_likelihood `Either new_tree
                    then 0.0
                else begin
                    Tree.post_order_node_with_edge_visit_simple
                            distance (Tree.Edge (a, b)) new_tree.Ptree.tree
                            (~-. (distance b a 0.0))
                end
        in
        let pi_cost = prior_cost new_tree None in (* all chars = None *)
        let root_cost = AllDirNode.AllDirF.root_cost root in
        if debug then begin
            let () = match root_edge with
                | `Single x   -> info_user_message "Cost of: %d" x
                | `Edge (a,b) -> info_user_message "Cost from: (%d,%d)" a b
            in
            info_user_message "Single Character Cost: %f" single_characters_cost;
            info_user_message "Other Character Cost: %f" not_single_character_cost;
            info_user_message "Root Cost: %f" (root_cost);
            info_user_message "Prior Cost: %f" (pi_cost);
            info_user_message "Size of Tree: %f" (tree_size new_tree None);
        end;
        let res =
            single_characters_cost +. not_single_character_cost +. root_cost +. pi_cost
        in
        res

    (* above function over all handles *)
    let check_cost_all_handles ptree =
        IntSet.fold 
            (fun handle cost -> (check_cost ptree handle None) +. cost)
            (Ptree.get_handles ptree)
            0.0

    let root_costs tree =
        let prior = prior_cost tree None in
        let collect_edge_data edge node acc =
            let cost = AllDirNode.OneDirF.tree_cost None node in
            (edge, cost +. prior) :: acc
        in
        Tree.EdgeMap.fold collect_edge_data tree.Ptree.edge_data []

    let report_all_roots tree =
        List.iter
            (fun (Tree.Edge (a,b),c) ->
                Printf.printf "root at (%d,%d) -- %f\n" a b c)
            (root_costs tree)


    let check_assertion_two_nbrs a b c =
        if a <> Tree.get_id b then true
        else 
            let () = Status.user_message Status.Error c in
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
            IntSet.fold get_handle
                ptree.Ptree.tree.Tree.handles IntSet.empty
        in
        pre_codes


    let clear_internals force t = t

    let print_all_nodes ptree =
        All_sets.IntegerMap.iter
            (fun k v ->
                Printf.printf "Component Root %d - (%s) - " k
                    (match v.Ptree.root_median with
                        | Some (`Edge (x,y), n) ->
                            string_of_int x ^ "," ^ string_of_int y
                        | Some (`Single x,n) ->
                            string_of_int x
                        | None -> "none"))
            ptree.Ptree.component_root;
        All_sets.IntegerMap.iter
            (fun k v -> AllDirNode.q_print v)
            (ptree.Ptree.node_data)

    (* A function to assign a unique sequence on each vertex of the ptree in the
    * [AllDirNode.adjusted] field of the node. *)
    let assign_single ptree = 
        (* We first define a function that can traverse the tree and assign
        * a single sequence to each vertex on it. *)
        let pre_ref_codes = get_pre_active_ref_code ptree in  
        let fi_ref_codes = pre_ref_codes in 
        let rec assign_single_subtree root_edge parentd parent current ptree =
            if debug_single_assignment then
                info_user_message "assign single subtree on node %d,parent=%d" current parent;
            let current_d, initial_d =
                let tmp = Ptree.get_node_data current ptree in
                AllDirNode.not_with parent  tmp.AllDirNode.unadjusted, tmp
            in
            let nd, original =
                let x = AllDirNode.force_val current_d.AllDirNode.lazy_node in
                let n = Node.to_single (pre_ref_codes, fi_ref_codes) root_edge None parentd x in
                n, x
            in
            let nnd = 
                { current_d with AllDirNode.lazy_node = AllDirNode.lazy_from_val nd }
            in
            let final_d = { initial_d with AllDirNode.adjusted = Some nnd } in
            let ptree = Ptree.add_node_data current final_d ptree in
            try 
                let a, b =
                    let currentn = Ptree.get_node current ptree in 
                    assert (check_assertion_two_nbrs parent currentn "2");
                    Tree.other_two_nbrs parent currentn 
                in
                ptree
                    --> assign_single_subtree false nd current a
                    --> assign_single_subtree false nd current b 
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
                let a, b = match edge with
                    | `Edge x ->  x
                    | `Single a -> a, a
                in
                if debug_single_assignment then
                    info_user_message "assign single subtree on root node %d,%d" a b;
                (*get root unadjusted node data as root*)
                let root, rooth = get_root_direction rootg in
                (*get handle node and other node*)
                let handle_node =
                    let hnode = Ptree.get_node_data a ptree in
                    hnode.AllDirNode.unadjusted
                        --> AllDirNode.not_with b
                        --> (fun x -> AllDirNode.force_val x.AllDirNode.lazy_node)
                and other_node =
                    (Ptree.get_node_data b ptree).AllDirNode.unadjusted 
                        --> AllDirNode.not_with a
                        --> (fun x -> AllDirNode.force_val x.AllDirNode.lazy_node)
                in
		        if debug_assign_single_handle then begin
                    Printf.printf "handle node before to single:\n%!";
                    Node.print handle_node;
			        Printf.printf "root node before to single\n%!";
		            Node.print root; 
                end;
                (*update single assignment to root*)
                let root_tosingle =
                    Node.to_single (pre_ref_codes, fi_ref_codes) true
                                   (Some root) other_node handle_node
                in
                (*update rooth with new lazy_node*)
                let rooti = 
                    { rooth with
                        AllDirNode.lazy_node = AllDirNode.lazy_from_val (root_tosingle) }
                in
                (*update rootg with new adjusted*)
                let readjusted = { rootg with AllDirNode.adjusted = Some rooti} in
                if debug_assign_single_handle then begin
                    Printf.printf "handle node after to single\n%!";
                    AllDirNode.print_node_data readjusted false; 
                end;
                let treecost = Node.Standard.tree_cost None root_tosingle in
                if debug_assign_single_handle then  
                    info_user_message "assign_single, assign cost to root:%f" treecost;
                let ptree = 
                    Ptree.assign_root_to_connected_component 
                        handle (Some (edge, readjusted)) 
                        treecost None ptree
                in
                ptree, root, readjusted
            in
            let comp = Ptree.get_component_root handle ptree in
            match comp.Ptree.root_median with
            | Some ((`Edge (a, b)) as edge, rootg) ->
                    let ptree, root, readjusted =
                        generate_root_and_assign_it rootg edge ptree 
                    in
                    ptree
                        --> assign_single_subtree true root b a
                        --> assign_single_subtree true root a b
                        --> (fun ptree ->
                if debug_assign_single_handle then  info_user_message "assign cost to root(%d,%d) with %f and check_cost function" a b comp.Ptree.component_cost;
                                Ptree.assign_root_to_connected_component 
                                    handle
                                    (Some ((`Edge (a, b)), readjusted))
                                    comp.Ptree.component_cost
                                    (Some (check_cost ptree handle None))
                                    ptree)
            | Some ((`Single a) as edge, rootg) ->
                    let ptree, root, readjusted =
                        generate_root_and_assign_it rootg edge ptree
                    in
                    Ptree.add_node_data a readjusted ptree
            | None -> 
                    failwith "no root? AllDirChar.assign_single_handle"
        in
        (* Finally, we are ready to proceed on all the handles available *)
        let res = IntSet.fold assign_single_handle (Ptree.get_handles ptree) ptree in
        res

    let assign_single ptree =
        if using_likelihood `OnlyStatic ptree
            then ptree
            else assign_single ptree

    let unadjust ptree = ptree
 
    let refresh_edge_data ptree = 
        let refresh_edge ((Tree.Edge (a,b)) as e) (acc,ptree) =
            if debug_uppass_fn then
                info_user_message "Refreshing %d--%d as an edge" a b;
            let data,ptree = create_lazy_edge false None ptree a b in
            (Tree.EdgeMap.add e data acc, ptree)
        in
        if debug_uppass_fn then
            info_user_message "Performing Refresh on all edges";
        let new_edges,ptree =
            Tree.EdgeSet.fold
                    (refresh_edge)
                    (ptree.Ptree.tree.Tree.d_edges)
                    (Tree.EdgeMap.empty, ptree)
        in
        { ptree with Ptree.edge_data = new_edges }


    (** Performs much of the node functions after a downpass is completed.
     * Essentially the uppass; this function fills in other directions of all the
     * nodes; their roots, and transfer branch lengths or other data. This
     * function requires only downpass data, and can be used to update a tree if
     * any nodes change --based on a change in the traversal *)
    let refresh_all_edges root_opt do_roots start_edge_opt ptree =
        if debug_uppass_fn then info_user_message "Refresh all edges, do_roots = %b" do_roots;
        (* A function to refresh the data on a edge *)
        let refresh_edge rhandle root_opt ((Tree.Edge (a,b)) as e) (acc,ptree) =
            if debug_uppass_fn then
                info_user_message "Refreshing %d--%d as %s" a b 
                                (if rhandle then "a root edge" else "an edge");
            let data,ptree =
                if rhandle then
                    let p1,p2 = Ptree.create_partition ptree e in
                    match hashdoublefind ptree [(p1,a);(p2,b)] with
                    | Some x -> create_lazy_edge ~branches:x rhandle root_opt ptree a b 
                    | None   -> create_lazy_edge rhandle root_opt ptree a b 
                else
                    create_lazy_edge rhandle root_opt ptree a b 
            in
            (Tree.EdgeMap.add e data acc,ptree)
        (* perform uppass heuristic on a node; fill in three directions *)
        and add_vertex_pre_order prev code (ptree:phylogeny) =
            match Ptree.get_node code ptree with
            | Tree.Single _ -> ptree
            | Tree.Leaf (_, p) ->
                    let this = Ptree.get_node_data code ptree
                    and parn = Ptree.get_node_data p ptree in
                    let leaf = AllDirNode.AllDirF.apply_time false this parn in
                    Ptree.add_node_data code leaf ptree
            | (Tree.Interior (_, par,a ,b)) as v -> 
                    if debug_uppass_fn then info_user_message "add vertex with node#.%d(par=%d,a=%d,b=%d,prev=%d)" code par a b prev;
                    let a,b = Tree.other_two_nbrs prev v in
                    let interior = create_lazy_interior_up ptree code a b prev in
                    if debug_uppass_fn then info_user_message "add internal node#.%d(a=%d,b=%d,prev=%d as par) to ptree"  code a b prev;
                    Ptree.add_node_data code interior ptree
        in
        (* Because the tree is currently disjoint; the final edge has no data;
         * and in likelihood terms, no length, we must fill this edge in first.
         * As an improvement, we can create an edge set to avoid duplicating this
         * median later *)
        if debug_uppass_fn then
            info_user_message "Performing Calculation on Root" else ();
        let _,ptree = match start_edge_opt with
            | Some (a,b) ->
		        if debug_uppass_fn then 
		            info_user_message "call refresh_edge with start_edge(%d,%d)" a b;
                refresh_edge true root_opt (Tree.Edge (a,b)) (Tree.EdgeMap.empty,ptree) 
            | None ->
		        if debug_uppass_fn then 
		            info_user_message "no start edge, iter over all handles to refresh_edges";
                IntSet.fold 
                    (fun h (edge_map,ptree) ->
                        try begin
                        match (Ptree.get_component_root h ptree).Ptree.root_median with
                            | Some ((`Edge (a,b)),c) -> 
                                refresh_edge true root_opt (Tree.Edge (a,b))
                                             (edge_map,ptree) 
                            | None
                            | Some _ -> (edge_map,ptree)
                        end with | Not_found -> 
                            begin match Ptree.get_node h ptree with
                                | Tree.Leaf (a,b)
                                | Tree.Interior (a,b,_,_) -> 
                                    refresh_edge true root_opt (Tree.Edge (a,b)) (edge_map,ptree)
                                | Tree.Single _ -> 
                                    (edge_map,ptree) 
                            end)
                    ptree.Ptree.tree.Tree.handles
                    (Tree.EdgeMap.empty,ptree)
        in
        (* perform uppass heuristic --fill all directions *)
        current_snapshot "AllDirChar refresh_all_edges uppass heuristic";
        if debug_uppass_fn then
            info_user_message "Performing Uppass Heuristic";
        let ptree = match start_edge_opt with
            | Some (a,b) ->
                Tree.pre_order_node_with_edge_visit_simple_root
                    add_vertex_pre_order (Tree.Edge (a,b)) ptree.Ptree.tree ptree
            | None ->
                IntSet.fold
                    (fun h ptree ->
                        try begin
                        match (Ptree.get_component_root h ptree).Ptree.root_median with
                            | Some ((`Edge (a,b)),c) ->
                                Tree.pre_order_node_with_edge_visit_simple_root
                                        add_vertex_pre_order (Tree.Edge (a,b))
                                        ptree.Ptree.tree ptree
                            | None
                            | Some _ -> ptree
                        end with | Not_found ->
                            begin match Ptree.get_node h ptree with
                                | Tree.Leaf (a,b)
                                | Tree.Interior (a,b,_,_) ->
                                    Tree.pre_order_node_with_edge_visit_simple_root
                                            add_vertex_pre_order (Tree.Edge (a,b))
                                            ptree.Ptree.tree ptree
                                | Tree.Single _ -> ptree
                            end)
                    ptree.Ptree.tree.Tree.handles
                    ptree
        in
        (* fill in roots for all edges *)
        current_snapshot "AllDirChar refresh_all_edges internal fold";
        if do_roots then
            refresh_edge_data ptree
        else 
            ptree

    (* return the unadjusted value in the direction to parent in one dir *)
    let get_unadjusted par node =
        let nd = AllDirNode.not_with par (node.AllDirNode.unadjusted) in
        AllDirNode.force_val nd.AllDirNode.lazy_node

    (** refresh root of all trees *)
    let refresh_roots ptree =
        let new_roots = 
            IntSet.fold
                (fun x acc ->
                    let root = match Ptree.get_node x ptree with
                        | Tree.Leaf (a, b)
                        | Tree.Interior (a, b, _, _) -> create_root a b ptree
                        | Tree.Single _ -> create_root x x ptree
                    in
                    IntMap.add x root acc)
                (Ptree.get_handles ptree)
                (IntMap.empty)
        in
        { ptree with Ptree.component_root = new_roots; }

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
        IntSet.fold
            get_active_ref_code_handle
            (Ptree.get_handles tree)
            (IntSet.empty, IntSet.empty)

    let add_component_root ptree handle root = 
        { ptree with
            Ptree.component_root = IntMap.add handle root ptree.Ptree.component_root }

    let reroot_fn n_mgr force edge ptree =
        let Tree.Edge (h, n) = edge in
        let my_handle = Ptree.handle_of h ptree in
        let root = Ptree.get_component_root my_handle ptree in
        let ptree, _ = 
            ptree --> Ptree.remove_root_of_component my_handle 
                  --> Ptree.move_handle h 
        in
        let ptree = Ptree.fix_handle_neighbor h n ptree in
        let tree,inc = match !Methods.cost with
            | `Exhaustive_Strong
            | `Exhaustive_Weak
            | `Normal_plus_Vitamines
            | `Iterative `ApproxD _
            | `Normal ->
		        if debug_reroot then Printf.printf "reroot_fn, Normal,Iterative2D,etc..call create_root\n%!";
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
		        if debug_reroot then Printf.printf "reroot_fn, Iterative3D, just add component root\n%!";
                add_component_root ptree h root, []
        in
        update_node_manager tree (`Reroot inc) n_mgr;
        (tree,inc)

    let print_tree_times ptree = 
        let str_dir x = match x.AllDirNode.dir with
            | Some (x,y) -> (string_of_int x)^","^(string_of_int y)
            | None       -> "none"
        in
        let get_node x = AllDirNode.force_val x.AllDirNode.lazy_node in
        All_sets.IntegerMap.iter
            (fun k v -> 
                List.iter
                    (fun ndir ->
                        Printf.printf "u- (%s)\t--" (str_dir ndir);
                        Node.print_times (get_node ndir) )
                    v.AllDirNode.unadjusted;
                match v.AllDirNode.adjusted with
                | Some v -> Printf.printf "a- (%s)\t--" (str_dir v);
                            Node.print_times (get_node v)
                | None   -> Printf.printf "a- %d no adjusted\n" k)
            ptree.Ptree.node_data

    let update_edge_from_child_node single cn pn (tree : phylogeny) : phylogeny = 
        (* the child contains branch information to parent, but parent does not.
           We need to update this information, then call, refresh_all_edges with
           this new node data to update directions. *)
        let child  = Ptree.get_node_data cn tree
        and parent = Ptree.get_node_data pn tree in
        let edgenode = AllDirNode.create_root_from_child_branch child parent in
        let tree = refresh_all_edges (Some edgenode) false (Some (cn,pn)) tree in
        let tree = if single then assign_single tree else tree in
        tree

    (** We define a function that can optimize a subset of branches in the tree
        to improve the overall cost of the tree. The arguments for this function
        are as follows,
            - max_count : number of passes to do.
            - branches  : the branches in the tree to iterate (likelihood only).
            - ptree     : the tree to adjust. *)
    let branch_fn max_count branches ptree =
        let max_count = match max_count with
            | Some x -> x
            | None   -> max_int
        and first_affected = IntMap.map (fun _ -> None) ptree.Ptree.node_data in
        (* adjust root --for likelihood; after completion of iteration, we also
           need to iterate the final edge, the root, of the tree. *)
        let adjust_branches (changed,affected,ptree) c2c handle a b : adjust_acc =
            if debug_adjust_fn then
                    info_user_message "Adjusting root with %d,%d then None" a b;
            (* find edge data and package for AllDir optimization *)
            let new_root =
                let a_nd = Ptree.get_node_data a ptree
                and b_nd = Ptree.get_node_data b ptree in
                let o_nd =
                    let one = Ptree.get_edge_data (Tree.Edge (a,b)) ptree in
                    let tmp = { AllDirNode.lazy_node = one;
                                                 dir = Some (a,b);
                                                code = ~-1; }
                    in
                    { AllDirNode.adjusted= None; AllDirNode.unadjusted=[tmp] }
                in
                AllDirNode.AllDirF.edge_iterator None o_nd a_nd b_nd
            in
            (* below we apply the new branch length data to the left and right
               and add the result to the ptree. n_root is striped of direction.*)
            let e_root,ptree = create_lazy_edge true (Some new_root) ptree a b in
            let n_root = { AllDirNode.lazy_node = e_root;
                                          dir = Some (a,b);
                                         code = ~-1; } in
            let n_root = { AllDirNode.adjusted= None;
                           AllDirNode.unadjusted=[n_root] }
            in
            (* assign the root and cost *)
            let ptree = refresh_all_edges (Some n_root) true (Some (a,b)) ptree in
            let treecost = AllDirNode.OneDirF.tree_cost None e_root in
		let debug = false in
	    if debug then Printf.printf "adjust_branches,asign root with (%f,None)\n%!" treecost;
            let ptree =
                let root_edge = (Some (`Edge (a,b),n_root)) in
                Ptree.assign_root_to_connected_component handle root_edge
                                            treecost None ptree
            in
            let ptree = assign_single ptree in
            (changed,affected,ptree)
        in
        (* loop to adjust a tree and *)
        let adjust_until_nothing_changes max_count start_ptree =
            (* loop for rerooting and applying iterative on the resultant path *)
            let adjust_reroot_loop affected (modified,aff_n,ptree) (a,b) =
                (* a simple reroot, since the reroot_fn requires incremental as
                 * a return type, and a nodes_manager; which this function does *)
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
                adjust_branches (modified,aff_n,ptree) affected handle a b
            in
            let all_edges = match branches with
                | Some x -> (* strip the constructor; Tree.Edge *)
                    List.map (fun (Tree.Edge (a,b)) -> (a,b)) x
                    --> Array.of_list
                | None   -> (* all edges *)
                    Tree.EdgeMap.fold
                        (fun (Tree.Edge (a,b)) _ acc -> (a,b)::acc)
                        ptree.Ptree.edge_data
                        []
                    --> Array.of_list
            in
            (* recursive loop of for changes *)
            let rec iterator count prev_cost affected ptree =
                let (changed,new_affected,new_ptree : adjust_acc) = 
                    let none_affected = IntMap.empty in
                    let () = Array_ops.randomize all_edges in
                    Array.fold_left
                        (adjust_reroot_loop affected)
                        (true,none_affected,ptree)
                        (all_edges)
                in
                (* now ptree can be used normaliy *)
                let new_cost = check_cost_all_handles new_ptree in
                if debug_adjust_fn then
                    info_user_message "Iteration %d completed: %f --> %f (%b)" 
                                      (max_count - count) prev_cost new_cost changed;
                if (not changed) || (count = 1) || (prev_cost =. new_cost) || (new_cost > prev_cost)
                    then begin
                        if debug_adjust_fn then Printf.printf "return OLD ptree\n%!";
                        ptree
                    end
                    else begin
                        if debug_adjust_fn then Printf.printf "call iterator again with new ptree and new_cost=%f\n%!" new_cost;
                        iterator (count - 1) new_cost new_affected new_ptree
                    end
            in
            let initial_cost = check_cost_all_handles ptree in
            if debug_adjust_fn then Printf.printf "call iterator with init cost:%f\n%!" initial_cost;
            iterator max_count initial_cost first_affected ptree
        in
        adjust_until_nothing_changes max_count ptree

    (** We define a function that can adjust all the vertices in the tree
        to improve the overall cost of the tree, using only the
        [AllDirNode.adjusted] field of each. This function is aimed at
        optimization of the assignments at the internal nodes of the tree. The
        arguments for this function are as follows,
            - max_count : option for maximum number of optimization passes
            - nodes     : the tested nodes/characters to avoid
            - ptree     : the tree to adjust *)
    let adjust_assignment max_count nodes ptree =
        let mode = match !Methods.cost with
            | `Iterative x -> x
            | _            -> assert false
        and max_count = match max_count with
            | Some x -> x
            | None   -> max_int
        and first_affected = match nodes with
            | None   -> IntMap.map (fun _ -> None) ptree.Ptree.node_data
            | Some i -> i
        in
        (* We start by defining a function to adjust one node *)
        let adjust_node chars_to_check ch1_k ch2_k parent_k mine_k ptree =
            current_snapshot (Printf.sprintf "AllDirChar.adjust_node %d" mine_k);
            if debug_adjust_fn then
                info_user_message "AllDirChar.adjust_node, on %d, (%d,%d), %d"
                                        mine_k ch1_k ch2_k parent_k;
            let gnd x = Ptree.get_node_data x ptree in
            let mine,modified =
                AllDirNode.AllDirF.readjust mode chars_to_check (gnd ch1_k)
                                        (gnd ch2_k) (gnd parent_k) (gnd mine_k)
            in
            if IntSet.is_empty modified
                then modified,[],ptree
                else begin
                    let ptree = Ptree.add_node_data mine_k mine ptree in
                    modified, [ch1_k;ch2_k;mine_k;parent_k], ptree
                end
        in
        (* add modified vertices in node_list to the set *)
        (*node_list : the id list of node we need to update, 
        * codes : an IntSet of 
        * int list(character id we need to update) 
        * character codes 
        * affected i:an IntMap of previous affected nodes.
        * int(node id) -> (Some IntSet of character id we need to update) *)
        let add_vertices_affected node_list codes affected =
            let add_one code affected =
                if IntMap.mem code affected then
                    begin match IntMap.find code affected with
                        | Some conts ->
                            let res = IntSet.union conts codes in
                            IntMap.add code (Some res) affected
                        | None ->
                                IntMap.add code (Some codes) affected
                                (*since we init adjust function with
                                * first_affected: a IntMap of
                                * nodeid -> None, we should not assert false here
                                * assert false*)
                    end
                else
                    IntMap.add code (Some codes) affected
            in
            List.fold_right (fun x acc -> add_one x acc) node_list affected
        in
        (* compose the above functions to adjust and modify the affected nodes *)
        let adjust_vertices_affected ((modified,affected_nodes,ptree) as acc) c2c prev curr =
            if not (IntMap.mem curr c2c) then acc
            else match Ptree.get_node curr ptree with 
                | (Tree.Interior (c,p,c1,c2)) as nd ->
                        if debug_adjust_fn then begin
                            Printf.printf "adjust_vertices_affected, affected_nodes:\n%!";
                            IntMap.iter(fun key v ->
                                Printf.printf "key:%d --> chacter set:" key;
                                match v with
                                | Some cset -> 
                                        Printf.printf "[%!";
                                        IntSet.iter (fun c -> Printf.printf "%d," c)
                                        cset;
                                        Printf.printf "]\n%!"
                                | None ->
                                        Printf.printf "None !\n%!"
                            ) affected_nodes;
                        end;
                    let c2c = IntMap.find curr c2c in
                    let a,b = Tree.other_two_nbrs prev nd in
                    let ccodes,affected,n_ptree = adjust_node c2c a b prev curr ptree in
                    if debug_adjust_fn then begin
                        Printf.printf "end of adjust_node(m:%d,p:%d,c1:%d,c2:%d), \
                                character codes we need to udpate:[%!" c p c1 c2;
                        IntSet.iter (fun x -> Printf.printf "%d," x) ccodes;
                        Printf.printf "]\n%!";
                    end;
                    let n_ptree = match using_likelihood `Dynamic n_ptree with
                        (* we need to update the branch length from optimization *)
                        | true  -> update_edge_from_child_node true curr prev n_ptree
                        | false -> n_ptree
                    in
                    let new_affected = add_vertices_affected affected ccodes affected_nodes
                    and modified = (0 != (List.length affected)) || modified in
                    (modified, new_affected, n_ptree)
                | Tree.Leaf _ 
                | Tree.Single _ ->  acc
        in
        (* loop to adjust a tree and *)
        let adjust_until_nothing_changes max_count start_ptree =
            (* Post order traversal of internal nodes *)
            let debug = false in
            if debug then Printf.printf "adjust until nothing changes begin\n%!";
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
            in
            (* recursive loop of for changes *)
            let rec iterator count prev_cost affected ptree =
                let (changed,new_affected,new_ptree : adjust_acc) = 
                    let none_affected = IntMap.empty in
                    IntSet.fold
                            (adjust_loop affected)
                            (ptree.Ptree.tree.Tree.handles)
                            (true,none_affected,ptree)
                in
                (* now ptree can be used normaliy *)
                if debug_adjust_fn then Printf.printf "call check_cost_all_handles :\n%!";
                let new_cost = check_cost_all_handles new_ptree in
                if debug_adjust_fn then
                    Printf.printf "Iteration %d completed: prev cost:%f --> new cost:%f anything changed:%b\n%!"
                                  (max_count - count) prev_cost new_cost changed;
                if (not changed) || (count = 1) || (prev_cost =. new_cost) || (new_cost > prev_cost)
                    then begin
                        if debug_adjust_fn then
                            Printf.printf "end of adjust_until_nothing_changes; return old\n%!";
                        ptree
                    end else begin
                        iterator (count - 1) new_cost new_affected new_ptree
                    end
            in
            let initial_cost = check_cost_all_handles ptree in
            if debug_adjust_fn then
                Printf.printf "call iterator with initial cost = %f\n%!" initial_cost;
            iterator max_count initial_cost first_affected ptree
        in
        let set_handle_n_root_n_cost handle ptree =
            (*move AllDirF.to_single here,we are going to move to_single
            * related function in allDirChar to allDirNode *)
            let quick_2single root b d set =
                let b',d' =
                    let one = match b.AllDirNode.adjusted with
                        | None   -> assert false
                        | Some x -> x
                    and two = match d.AllDirNode.adjusted with
                        | None   -> assert false
                        | Some y -> y
                    in
                    one,two
                in
                let root = match root.AllDirNode.unadjusted with
                    | [x] -> Some x.AllDirNode.lazy_node
                    |  _  -> assert false
                in
                let lazy_node =
                    AllDirNode.OneDirF.to_single root None 
                        b'.AllDirNode.lazy_node None d'.AllDirNode.lazy_node set
                in
                let node = { d' with AllDirNode.lazy_node = lazy_node } in
                { AllDirNode.unadjusted = [node]; AllDirNode.adjusted = Some node }
            in
            (* quick_2single was AllDirF.to_single in allDirNode.ml; should be
               moved back when things are stable. *)
            let comp_root = Ptree.get_component_root handle ptree in
            match comp_root.Ptree.root_median with
            | None -> assert false
            | Some ((`Edge (a, b)) as edge, root) ->
                let sets = get_active_ref_code ptree
                and ad = Ptree.get_node_data a ptree
                and bd = Ptree.get_node_data b ptree in
                let root = 
                    let n = AllDirNode.AllDirF.median None None ad bd in
                    quick_2single n ad bd sets 
                    (*AllDirNode.AllDirF.to_single (Some n) None ad None bd sets*)
                in
                Ptree.assign_root_to_connected_component 
                        handle (Some (edge, root)) 
                        (check_cost ptree handle None) None ptree
            | Some _ -> ptree
        in
        let newtree = adjust_until_nothing_changes max_count ptree in
        let ptree =
            IntSet.fold (set_handle_n_root_n_cost)
                        (ptree.Ptree.tree.Tree.handles)
                        (newtree)
        in
        ptree

    let adjust_assignment max_count nodes ptree =
        if using_likelihood `Either ptree
            then ptree
            else adjust_assignment max_count nodes ptree


    (* ------------------------------------------------------------------------ *)
    (** [internal_downpass] Traverse every vertex in the tree and assign the
     * downpass and uppass information using the lazy all direction nodes *)
    let internal_downpass do_roots (ptree : phylogeny) : phylogeny =
        let add_vertex_post_order prev code ptree =
            current_snapshot "AllDirChar.internal_downpass.add_vertex_post_order";
            match Ptree.get_node code ptree with
            | Tree.Single _
            | Tree.Leaf (_, _) -> 
                assert (IntMap.mem code ptree.Ptree.node_data);
                if debug_downpass_fn then
                    info_user_message "Skipping Leaf/Single %d%!" code; 
                ptree
            | (Tree.Interior (_, par, a, b)) as v ->
                let a,b = Tree.other_two_nbrs prev v in
                if debug_downpass_fn then
                    info_user_message "Adding Vertex %d post Order: (%d,%d) and %d%!"
                                        code a b prev;
                let interior =
                    let p1 = fst (Ptree.create_partition ptree (Tree.Edge (b,code))), b
                    and p2 = fst (Ptree.create_partition ptree (Tree.Edge (a,code))), a in
                    match hashdoublefind ptree [p1;p2] with
                    | Some x -> create_lazy_interior_down ~branches:x ptree (Some code) a b
                    | None   -> create_lazy_interior_down ptree (Some code) a b
                in
                Ptree.add_node_data code interior ptree
        in
        let ptree =
            IntSet.fold
                (fun x (ptree:phylogeny) ->
                    try begin
                        match (Ptree.get_component_root x ptree).Ptree.root_median with
                        | Some ((`Edge (a,b)),c) ->
                            if debug_downpass_fn then
                                info_user_message "Downpass from Given (%d,%d)" a b;
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
                            if debug_downpass_fn then
                                info_user_message "Downpass from Handle (%d,%d)" a b;
                            Tree.post_order_node_with_edge_visit_simple
                                add_vertex_post_order
                                (Tree.Edge (a,b))
                                ptree.Ptree.tree ptree
                        | Tree.Single _ -> ptree
                    end)
                ptree.Ptree.tree.Tree.handles
                ptree
        in
        let ptree = refresh_all_edges None true None ptree in
        if do_roots then refresh_roots ptree else ptree


    let blindly_trust_downpass ptree (edges, handle) (cost, cbt) ((Tree.Edge (a, b)) as e) =
        let data = Ptree.get_edge_data e ptree in
        let c = AllDirNode.OneDirF.tree_cost None data in
        if abs_float cost > abs_float c then
            let data =  { AllDirNode.lazy_node = data; dir = None; code = -1 } in
            let data = { AllDirNode.unadjusted = [data]; adjusted = Some data } in
            let comp = Some ((`Edge (a, b)), data) in
            c, 
            Lazy.lazy_from_fun
                (fun () ->
                    Ptree.assign_root_to_connected_component handle comp c None ptree)
        else 
            (cost, cbt)


    let general_pick_best_root selection_method ptree =
        let edgesnhandles = 
            IntSet.fold 
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
                    Ptree.get_cost `Adjusted ptree, lazy ptree
            in
            let _, ptree =
                List.fold_left 
                    (selection_method ptree (edges, handle))
                     current_root_of_tree
                     (List.sort 
                            (fun (Tree.Edge (a, b)) (Tree.Edge (c, d)) ->
                                match c - a with
                                | 0 -> d - b
                                | x -> x)
                            edges)
            in
            Lazy.force_val ptree
        in 
        List.fold_left process ptree edgesnhandles 


    let pick_best_root ptree =
        if using_likelihood `OnlyStatic ptree then ptree
        else general_pick_best_root blindly_trust_downpass ptree


    (* ----------------- *)
    (* function to adjust the likelihood model ;for static characters of a tree
     * using BFGS --quasi newtons method. Function requires three directions. *)
    let static_model_chars_fn chars tree = 
        assert( using_likelihood `Static tree );
        (* replace nodes in a tree, copying relevent data structures *)
        let substitute_nodes nodes tree =
            let adder acc x = IntMap.add (AllDirNode.AllDirF.taxon_code x) x acc in
            let node_data = List.fold_left adder IntMap.empty nodes in
            internal_downpass true {tree with Ptree.node_data = node_data}
        in
        (* function for processing a model and applying to a tree --inner loop *)
        let f_likelihood f tree chars current_model new_values =
            let ntree =
                new_values
                    --> f current_model
                    --> Data.apply_likelihood_model_on_chars tree.Ptree.data chars
                    --> AllDirNode.AllDirF.load_data
                    --> (fun (x,y) -> substitute_nodes y {tree with Ptree.data = x})
            in
            let ncost = Ptree.get_cost `Adjusted ntree in
            ntree,ncost
        and get_some = function | Some x -> x | None -> raise Not_found in
        (* compose above functions to initiate adjustments *)
        let current_model = Data.get_likelihood_model tree.Ptree.data chars
        and current_cost = Ptree.get_cost `Adjusted tree in
        (* optimize Transition Model parameters *)
        let best_tree, best_cost =
            let params = MlModel.get_current_parameters_for_model current_model in
            match MlModel.get_update_function_for_model current_model with
            | Some func ->
                let opt = Numerical.default_numerical_optimization_strategy
                                    !Methods.opt_mode (Array.length params) in
                let _,results =
                    Numerical.run_method opt (f_likelihood func tree chars current_model)
                                             (params,(tree,current_cost))
                in
                if debug_model_fn then
                    info_user_message "\tOptimized Rates to %f --> %f"
                                      current_cost (snd results);
                results
            | None -> (tree,current_cost)
        in
        (* Optimize Alpha *)
        let current_model = Data.get_likelihood_model best_tree.Ptree.data chars in
        let best_tree, best_cost = 
            match MlModel.get_update_function_for_alpha current_model with
            | None      -> best_tree,best_cost
            | Some func ->
                let current_a = MlModel.get_current_parameters_for_alpha current_model in
                let _,results = 
                    Numerical.brents_method (f_likelihood func best_tree chars current_model)
                                            ((get_some current_a),(best_tree,best_cost))
                                          
                in
                if debug_model_fn then
                    info_user_message "\tOptimized Alpha to %f --> %f"
                                      best_cost (snd results);
                results
        in
        if best_cost < current_cost then best_tree else tree

    (* Group all the characters and optimize each with function above *)
    let static_model_fn tree = 
        List.fold_left
            (fun tree -> function
                | []  -> tree
                | xs  -> static_model_chars_fn xs tree)
            (tree)
            (Data.categorize_likelihood_chars_by_model tree.Ptree.data `AllStatic)


    module IA = ImpliedAlignment.Make (AllDirNode.AllDirF) (Edge.LazyEdge)
    (* Optimize model of dynamic likelihood characters by converting to an
     * implied alignment. We do ONE pass; ensuring the improvement of costs
     * exists after the call. *)
    let dynamic_model_fn old_tree =
        (* define a function to convert and optimize static model *)
        let optimize_static_tree ptree =
            let old_verbosity = Status.get_verbosity () in
            let dlk_categories = Data.categorize_likelihood_chars_by_model
                                                ptree.Ptree.data `AllDynamic in
            List.fold_left
                (fun ptree chars ->
                    Status.set_verbosity `None;
                    let data,chars =
                        IA.to_static_homologies true IA.filter_characters true
                            false (`Some (true,chars)) ptree.Ptree.data ptree
                    in
                    let data,nodes = AllDirNode.AllDirF.load_data 
                                            ~silent:true ~classify:false data
                    in
                    let nodes =
                        List.fold_left
                            (fun acc x ->
                                IntMap.add (AllDirNode.AllDirF.taxon_code x)
                                           x acc)
                            IntMap.empty
                            nodes
                    in
                    Status.set_verbosity old_verbosity;
                    let tree = 
                        { ptree with Ptree.data = data;
                                     Ptree.node_data = nodes; }
                                --> internal_downpass true
                                --> static_model_chars_fn chars
                    in
                    tree)
                ptree
                dlk_categories
        (* define a function to apply model from static to dynamic *)
        and static_model_to_dyn_chars static_tree dyn_tree =
            let data, nodes =
                Data.sync_static_to_dynamic_model
                    ~src:static_tree.Ptree.data ~dest:dyn_tree.Ptree.data
                --> AllDirNode.AllDirF.load_data ~silent:true ~classify:false
            in
            let node_data =
                List.fold_left
                    (fun acc x -> IntMap.add (AllDirNode.AllDirF.taxon_code x) x acc)
                    IntMap.empty
                    nodes
            in
            { dyn_tree with Ptree.data      = data;
                            Ptree.node_data = node_data; }
        in
        let old_cost = Ptree.get_cost `Adjusted old_tree in
        let new_tree =
            let static_tree = optimize_static_tree old_tree in
            old_tree
                --> static_model_to_dyn_chars static_tree
                (** Below is a downpass + uppass *)
                --> internal_downpass true
                --> pick_best_root
                --> assign_single
        in
        let new_cost = Ptree.get_cost `Adjusted new_tree in
        if new_cost < old_cost then begin
            if debug_model_fn then
                info_user_message "Updated Dynamic Likelihood Score: %f --> %f"
                                  old_cost new_cost;
            new_tree
        end else begin
            old_tree
        end

    (* adjust the assignment of internal nodes, and models for a tree *)
    let model_fn ?max_iter node_man tree =
        let rates_fn tree =
            let tree =
                if using_likelihood `Static tree then static_model_fn tree else tree in
            let tree =
                if using_likelihood `Dynamic tree then dynamic_model_fn tree else tree in
            tree
        and max_iter = match max_iter with
            | Some i -> i
            | None   -> Numerical.default_number_of_passes !Methods.opt_mode
        in
        (* adjust model and branches -- for likelihood *)
        let adjust_ do_model do_branches branches iterations first_tree = 
            let rec loop_ iter icost itree =
                let mcost,mtree,iter =
                    if do_model then begin
                        let mtree = rates_fn itree in
                        let mcost = Ptree.get_cost `Adjusted mtree in
                        if debug_model_fn then
                            info_user_message "Step %d; Optimized Model %f --> %f" iter icost mcost;
                        mcost,mtree,iter+1
                    end else begin
                        icost,itree,iter
                    end
                in
                let bcost,btree,iter =
                    if do_branches then begin
                        let btree = branch_fn iterations branches mtree in
                        let btree = update_branches btree in
                        let bcost = Ptree.get_cost `Adjusted btree in
                        if debug_model_fn then
                            info_user_message "Step %d; Optimized Branches %f --> %f" iter mcost bcost;
                        bcost,btree,iter+1
                    end else begin
                        mcost,mtree,iter
                    end
                in
                if icost =. bcost || iter >= max_iter
                    then btree
                    else loop_ iter bcost btree
            in
            let first_cost = Ptree.get_cost `Adjusted first_tree in
            let first_tree = update_branches first_tree in
            loop_ 0 first_cost first_tree
        in
        let tree =
            let do_branches,branches,do_model = match node_man with
                | Some node_man ->
                    let do_branches = match node_man#branches with
                        | Some [] -> false
                        | _       -> true
                    in
                    do_branches, node_man#branches, node_man#model
                | None -> true,None,true
            in
            let do_model,do_branches =
                if `None = !Methods.opt_mode then false,false else do_model,do_branches
            in
            if (not do_branches) && (not do_model) then
                tree
            else
                let n_tree = adjust_ do_model do_branches branches None tree in
                if debug_model_fn then
                    info_user_message "Optimized Likelihood Params: %f to %f"
                        (Ptree.get_cost `Adjusted tree) (Ptree.get_cost `Adjusted n_tree);
                n_tree
        in
        tree

    let model_fn ?max_iter node_man ptree =
        if using_likelihood `Either ptree
            then model_fn ?max_iter node_man (assign_single ptree)
            else ptree


    (* ---------- *)
    let downpass ptree =
        if debug_downpass_fn then info_user_message "Downpass Begins,%!";
        current_snapshot "AllDirChar.downpass start";
        let ptree = match !Methods.cost with
            | `Exhaustive_Strong
            | `Exhaustive_Weak
            | `Normal_plus_Vitamines
            | `Normal -> internal_downpass true ptree
            | `Iterative (`ThreeD  iterations)
            | `Iterative (`ApproxD iterations) ->
                ptree --> clear_internals false
                      --> internal_downpass true
                      --> pick_best_root
                      (*we need assign_single before adjust function because
                      *we call [check_cost] to get an initial value*)
                      --> assign_single 
                      --> adjust_assignment iterations None
                      (*we need assign_single after adjust function or report(diagnosis) won't work *)
                      --> assign_single
        in
        current_snapshot "AllDirChar.downpass end";
        if debug_downpass_fn then info_user_message "Downpass Ends\n%!";
        ptree

    (* ----------------- *)
    let uppass ptree =
        if debug_uppass_fn then info_user_message "UPPASS begin:%!";
        current_snapshot "AllDirChar.uppass start";
        let ptree = match !Methods.cost with
            | `Exhaustive_Strong
            | `Exhaustive_Weak
            | `Normal_plus_Vitamines
            | `Normal ->
                ptree --> pick_best_root --> assign_single --> model_fn None
            | `Iterative (`ApproxD _)
            | `Iterative (`ThreeD _) -> ptree --> model_fn None
        in
        current_snapshot "AllDirChar.uppass end";
        if debug_uppass_fn then info_user_message "UPPASS ends.%!";
        ptree

    (* reset the data in the removed direction, by running essentially an uppass
     * heuristic to fill in the other data starting at (a,b). edge data is a
     * previous node data that can be used for edge data, and direction
     * information, IT MUST HAVE (a,b) IN ITS DIRECTIONS. *)
    let clear_up_over_edge (a, b) edge_data_opt ptree =
        match edge_data_opt with
        (* Under static likelihood, because of the pully principle, we can keep
           the keep use the node data as the edge data/cost of subtree. *)
        | Some edge when using_likelihood `OnlyStatic ptree ->
            let edge =
                let single = AllDirNode.with_both a b edge.AllDirNode.unadjusted in
                { AllDirNode.unadjusted = [single]; adjusted = None }
            in
            refresh_all_edges (Some edge) false (Some (a,b)) ptree
        (* lets move the root to this edge, that way a simple refresh all edges
         * will take care of the missing node and update all the uppass data *)
        | Some _
        | None -> refresh_all_edges None false (Some (a,b)) ptree


    let create_edge edge_l edge_r i_code ptree = 
        let node =
            AllDirNode.AllDirF.median (Some i_code) None
                    (Ptree.get_node_data edge_l ptree)
                    (Ptree.get_node_data edge_r ptree) in
        Ptree.add_node_data i_code node ptree


    let clean_ex_neighbor a b ptree = 
        let data = Ptree.get_node_data a ptree in
        let notwith_un = AllDirNode.not_with b data.AllDirNode.unadjusted in
        let node = { AllDirNode.unadjusted = [notwith_un]; adjusted = None; } in
        if debug_node_fn then
            info_user_message "Cleaning vertex %d with parent %d." a b;
        Ptree.add_node_data a node ptree

    let get_edge_n_force a b ptree =
        let data = Ptree.get_edge_data (Tree.Edge (a, b)) ptree in
        AllDirNode.force_val data

    let replace_topology tree ptree = { ptree with Ptree.tree = tree }

    (* break_fn has type handle * int (node) -> tree -> tree * delta * aux_data *)
    let break_fn (tree_node_id, clade_node_id) (ptree : phylogeny) =
        if debug_break_fn then
            info_user_message "break_fn, Breaking %d -- %d" tree_node_id clade_node_id;
        let ptree = clear_internals true ptree in
        let (Tree.Edge (tree_node_id, clade_node_id)) as edge =
            Tree.normalize_edge (Tree.Edge (tree_node_id, clade_node_id)) ptree.Ptree.tree
        in
        (* We want to know the cost of the tree, so we force the calculation of
           the downpass all the way down to the place of the breakage *)
        let prev_cost =
            let edge_res = get_edge_n_force tree_node_id clade_node_id ptree in
            Node.Standard.tree_cost None edge_res
        in
        (* Figure out the cost of the broken tree from subtrees *)
        let new_cost =
            let tree_cost =
                Ptree.get_node_data tree_node_id ptree
                    --> AllDirNode.AllDirF.total_cost (Some clade_node_id)
            and clade_cost =
                Ptree.get_node_data clade_node_id ptree
                    --> AllDirNode.AllDirF.total_cost (Some tree_node_id)
            in
            if debug_join_fn then begin
                info_user_message "Previous Cost: %f" prev_cost;
                info_user_message "New Cost = clade_cost(%f) + tree_cost(%f) = %f"
                    clade_cost tree_cost (clade_cost+.tree_cost)
            end;
            clade_cost +. tree_cost
        in
        (* Break the topology and update the data *)
        let ptree, tree_delta, clade_handle, tree_handle =
            (* A function that takes one side of a tree delta and updates the
            * tree's data using that information *)
            let update_break_delta delta ptree =
                match delta with
                | `Edge (rem, l1, l2, _) ->
                    assert ((tree_node_id = rem) || (clade_node_id = rem));
                    let old_data = Ptree.get_node_data rem ptree in
                    ptree --> clean_ex_neighbor l1 rem
                          --> clean_ex_neighbor l2 rem
                          --> Ptree.remove_node_data rem
                          --> clear_up_over_edge (l1, l2) (Some old_data)
                | _ -> ptree
            in
            (* Perform the topology break *)
            let nbt, ((left_delta, right_delta) as tree_delta) =
                Tree.break (tree_node_id, clade_node_id) ptree.Ptree.tree
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
                    --> refresh_edge_data
                    --> refresh_roots
            in
            (* The uppass can change the tree/clade_handles; update based on previous handles *)
            ptree, tree_delta, clade_handle, tree_handle
        in
        (* Compare costs, and calculate the break delta *)
        let b_delta =
            if prev_cost = infinity && new_cost = infinity then 0.
            else begin
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
                let bd = prev_cost -. (new_cost -. (rc +. ptree.Ptree.origin_cost)) in
                if debug_break_fn then begin
                    info_user_message "Root (%d) Cost: %f" clade_handle rc;
                    info_user_message "Total Cost: %f" tc;
                    info_user_message "Origin Cost: %f" ptree.Ptree.origin_cost;
                    info_user_message "New Cost: %f" new_cost;
                    info_user_message "Prev Cost: %f" new_cost;
                    info_user_message "Break Delta: %f" bd
                end;
                abs_float bd (* likelihood break-delta is negative *)
            end
        in
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
                  topology_delta = side; }
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


    (* ----------------- *)
    let break_fn n_mgr ((s1, s2) as a) b =
        let res = match !Methods.cost with
            | `Iterative (`ApproxD _)
            | `Iterative (`ThreeD _)
            | `Exhaustive_Weak
            | `Normal_plus_Vitamines
            | `Normal -> break_fn a b
            | `Exhaustive_Strong ->
                let breakage = break_fn a b in
                { breakage with 
                    Ptree.break_delta = (Ptree.get_cost `Adjusted b) -.
                                        (Ptree.get_cost `Adjusted breakage.Ptree.ptree);
                    Ptree.incremental = []; }
        in
        update_node_manager (res.Ptree.ptree) (`Break res) n_mgr;
        {res with
            Ptree.ptree = model_fn n_mgr res.Ptree.ptree; }


    (* ----------------- *)
    let join_fn _ jxn1 jxn2 ptree =
        if debug_join_fn then
            info_user_message "Joining: (%s) and (%s)"
                (match jxn1 with
                    | Tree.Single_Jxn x -> string_of_int x
                    | Tree.Edge_Jxn (x,y) -> (string_of_int x) ^","^ (string_of_int y))
                (match jxn2 with
                    | Tree.Single_Jxn x -> string_of_int x
                    | Tree.Edge_Jxn (x,y) -> (string_of_int x) ^","^ (string_of_int y));
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
                        --> create_edge a b v
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
                        --> create_edge a b r
                        --> create_edge c d v
                in
                r, v, ptree
            | (`Single (v, _)), (`Edge (r, a, b, Some h)), _ ->
                let ptree =
                    ptree --> Ptree.remove_root_of_component h
                        --> Ptree.remove_node_data r
                        --> clean_ex_neighbor a b
                        --> clean_ex_neighbor b a
                        --> create_edge a b r
                in
                v, r, ptree
            | _ -> failwith "Unexpected AllDirChar.join_fn"
        in
        let ptree = { ptree with Ptree.tree = ret } in
        let handle, parent =
            let handle = Ptree.handle_of v ptree in
            let parent = Ptree.get_parent handle ptree in
            handle, parent
        in
        let ptree =
            ptree --> Ptree.remove_root_of_component handle
                  --> refresh_all_edges None true (Some (v,h))
                  --> refresh_roots
        in
        if debug_join_fn then
            info_user_message "end of Joining, return ptree with cost adj:%f(unadj:%f)"
                (Ptree.get_cost `Adjusted ptree) (Ptree.get_cost `Unadjusted ptree);
        ptree, tree_delta

    let get_one side = match side with
        | `Single (x, _) | `Edge (x, _, _, _) -> x

    let join_fn n_mgr a b c d =
        let d = clear_internals true d in
        let (ptree, tdel) = match !Methods.cost with
            | `Normal ->
                let tree,delta =join_fn a b c d in
                update_node_manager tree (`Join delta) n_mgr;
                tree, delta
            | `Iterative (`ThreeD iterations)
            | `Iterative (`ApproxD iterations) ->
                let tree, delta = join_fn a b c d in
                update_node_manager tree (`Join delta) n_mgr;
                let tree =
                   tree --> pick_best_root
                        --> assign_single
                        --> adjust_assignment iterations None
                in
                tree, delta
            | `Normal_plus_Vitamines
            | `Exhaustive_Weak
            | `Exhaustive_Strong ->
                let tree, delta = join_fn a b c d in
                update_node_manager tree (`Join delta) n_mgr;
                uppass tree, delta
        in
        let ptree = model_fn n_mgr ptree in
        if debug_join_fn then
            info_user_message "Joined with cost: %f (%f)" 
                    (Ptree.get_cost `Adjusted ptree) (Ptree.get_cost `Unadjusted ptree);
        (ptree,tdel)


    type tmp = Edge of (int * int) | Clade of a
    let cost_fn jxn1 jxn2 delta clade (tree : phylogeny) =
        if debug_cost_fn then Printf.printf "alldirchar.cost_fn 2 -> %!";
        let rec forcer edge = match edge with
            | Edge (a, b) ->
                AllDirNode.force_val (Ptree.get_edge_data (Tree.Edge (a, b)) tree)
            | Clade x ->
                begin match x.AllDirNode.unadjusted with (* Leaf *)
                    | [x] -> force_node x
                    | _   -> assert false
                end
        in
        let clade_data = match !Methods.cost with
            | `Iterative (`ThreeD _) ->
                begin match jxn2 with
                    | Tree.Single_Jxn h    -> forcer (Clade clade)
                    | Tree.Edge_Jxn (h, n) ->
                        let Tree.Edge (h, n) =
                            Tree.normalize_edge (Tree.Edge (h, n)) tree.Ptree.tree
                        in
                        forcer (Edge (h, n))
                end
            | _ -> forcer (Clade clade)
        in
        let res = match jxn1 with
            | Tree.Single_Jxn h ->
                let d =
                    Node.Standard.distance 0.
                        (forcer (Clade (Ptree.get_node_data (Tree.int_of_id h) tree)))
                        clade_data
                in
                if debug_cost_fn then Printf.printf "single jxn,cost=%f\n%!" d;
                Ptree.Cost d
            | Tree.Edge_Jxn (h, n) ->
                let (Tree.Edge (h, n)) = 
                    Tree.normalize_edge (Tree.Edge (h, n)) tree.Ptree.tree
                in
                let ndata = forcer (Edge (h, n)) in
                let c = Node.Standard.distance 0. clade_data ndata in
                if debug_cost_fn then Printf.printf "edge jxn, cost=%f\n%!" c;
                Ptree.Cost c
        in
        res


    let cost_fn n_mgr a b c d e =
        if debug_cost_fn then Printf.printf "alldirchar.cost_fn 1 ->%!";
        let cost = match !Methods.cost with
            | `Iterative (`ApproxD _) ->
                begin match cost_fn a b c d e with 
                    | Ptree.Cost x -> Ptree.Cost (abs_float (0.85 *. x))
                    | x -> x
                end
            | `Iterative `ThreeD _
            | `Exhaustive_Weak
            | `Normal_plus_Vitamines
            | `Normal -> cost_fn a b c d e
            | `Exhaustive_Strong ->
                let pc = Ptree.get_cost `Adjusted e in
                let (nt, _) = join_fn n_mgr [] a b e in
                let res = Ptree.get_cost `Adjusted nt in
                Ptree.Cost (res -. pc)
        in
        if debug_cost_fn then begin 
            match cost with 
            | Ptree.Cost x ->  Printf.printf "cost = %f\n%!" x
            | Ptree.NoCost -> Printf.printf "NoCost\n%!"
        end;
        update_node_manager e (`Cost) n_mgr;
        cost

    let string_of_node _ = ""

    let features meth lst = 
        Chartree.features meth (("all directions", "true") :: lst)

    let incremental_uppass tree _ = tree

    let assign_final_states ptree =
        let assign_final_states_handle handle ptree =
            try let root_data, a, b =
                    let rt = Ptree.get_component_root handle ptree in
                    match rt.Ptree.root_median with
                    | Some ((`Edge (a, b)), root) -> root, a, b
                    | Some _ -> failwith "Single vertex" (* Used down below *)
                    | None -> failwith "No root?"
                in
                let root_data c = match root_data.AllDirNode.unadjusted with
                    (* We prepare a function to replace the taxon code for a
                       meaningful one to start the uppass with on each side *)
                    | [x] ->
                        { root_data with
                            AllDirNode.unadjusted = [{ x with AllDirNode.code = c }] }
                    | _ -> assert false
                in
                (* We move recursively up on a/b calculating their final states *)
                let rec uppass grandparent_code parent_code parent_final vertex acc =
                    let my_data = Ptree.get_node_data vertex ptree in
                    match Ptree.get_node vertex acc with
                    | Tree.Leaf _ -> acc
                    | Tree.Single _ -> acc
                    | (Tree.Interior _) as nd ->
                            let a, b = Tree.other_two_nbrs parent_code nd in
                            let nda = Ptree.get_node_data a ptree
                            and ndb = Ptree.get_node_data b ptree in
                            let my_data =
                                AllDirNode.AllDirF.final_states grandparent_code
                                                    parent_final my_data nda ndb
                            in
                            acc --> Ptree.add_node_data vertex my_data
                                --> uppass (Some parent_code) vertex my_data a
                                --> uppass (Some parent_code) vertex my_data b
                in
                ptree --> uppass None a (root_data a) b
                      --> uppass None b (root_data b) a
            with | Failure "Single vertex" -> ptree
        in
        IntSet.fold assign_final_states_handle (Ptree.get_handles ptree) ptree


    let to_formatter diag_report_type (atr : Xml.attributes) tree : Xml.xml =
        if debug_diagnosis then
            Printf.printf "AllDirChar.to_formatter \n%!";
        let get_single par node =
            if using_likelihood `OnlyStatic tree then
                get_unadjusted par node
            else match node.AllDirNode.adjusted with
                | Some x -> AllDirNode.force_val x.AllDirNode.lazy_node
                | None   -> failwithf "AllDirChar.to_formatter; no single data: %d -- %d"
                                        par (AllDirNode.AllDirF.taxon_code node)
        in
(*        let tree = assign_final_states tree in*)
        let pre_ref_codes, fi_ref_codes = get_active_ref_code tree in
        let get_simplified parent x =
            let nd = Ptree.get_node_data x tree in
            nd, get_unadjusted parent nd, get_single parent nd
        in
        let merger a b root = (`Set [`Single root; `Single a; `Single b])
        and splitter parent a = get_unadjusted parent a, get_single parent a in
        (* Now we are ready to process the contents of the tree *)
        let rec subtree_to_formatter diag_report_type (pre, fi) cur par ((node_parent, single_parent) as tmp2) : Xml.xml =
            if debug_diagnosis then
                Printf.printf "alldirchar.to_formatter cur:%d, par:%d\n%!" cur par;
            match Ptree.get_node cur tree with
            | (Tree.Interior (me,nbr1,nbr2,nbr3)) as nd ->
                if debug_diagnosis then
                    Printf.printf "Is an interior (%d,%d,%d,%d) => \n%!" me nbr1 nbr2 nbr3;
                let cur_data = Ptree.get_node_data cur tree in
                let ch1, ch2 = Ptree.other_two_nbrs par nd in
                let ch1d, ch1u, ch1s = get_simplified cur ch1
                and ch2d, ch2u, ch2s = get_simplified cur ch2 in
                let ((cur_data, cur_single) as tmp) = splitter par cur_data in
                let mine =
                    Node.to_formatter_subtree diag_report_type (pre, fi) [] tree.Ptree.data tmp cur
                                              (ch1, ch1u) (ch2, ch2u) (Some tmp2)
                in
                let ch1 = subtree_to_formatter diag_report_type (pre, fi) ch1 cur tmp in
                let ch2 = subtree_to_formatter diag_report_type (pre, fi) ch2 cur tmp in
                ((RXML
                    -[Xml.Trees.tree]
                        {single mine} { single ch1 }
                        { single ch2 } --) : Xml.xml)
            | (Tree.Leaf (me, par)) ->
                if debug_diagnosis then Printf.printf "Is a leaf (%d,%d) => %!" me par;
                let node_data = splitter par (Ptree.get_node_data cur tree) in
                let nodest =
                    Node.to_formatter_single diag_report_type (pre, fi) [] tree.Ptree.data
                                             node_data cur (Some tmp2)
                in
                (RXML -[Xml.Trees.tree] { single nodest }--)
            | (Tree.Single me) ->
                if debug_diagnosis then Printf.printf "Is a single: %d => %!" me;
                let node_data = splitter (-1) (Ptree.get_node_data cur tree) in
                let nodest =
                    Node.to_formatter_single diag_report_type (pre, fi) [] tree.Ptree.data
                                             (node_data) cur None
                in
                (RXML -[Xml.Trees.tree] { single nodest } --)
        in
        let handle_to_formatter diag_report_type (pre, fi) handle (recost, trees) =
            let r = Ptree.get_component_root handle tree in
            let recost, contents, attr = match r.Ptree.root_median with
                | Some ((`Edge (a, b)), root) ->
                    if debug_diagnosis then
                        Printf.printf "root median at (%d,%d) : %!" a b;
                    let recost =
                        let root = get_unadjusted (-1) root in
                        (Node.cmp_subtree_recost root) +. recost
                    in
                    (* We override the root now to continue using the single
                       assignment of the handle *)
                    let sroot, sa =
                        let a = Ptree.get_node_data a tree in
                        let s = get_single b a in
                        let root = get_unadjusted (-1) root in
                        let s_root = Node.copy_chrom_map root s in
                        (root, s_root), s
                    in
                    let a : Xml.xml = subtree_to_formatter diag_report_type (pre, fi) a b sroot
                    and b : Xml.xml = subtree_to_formatter diag_report_type (pre, fi) b a sroot
                    and froot : Xml.xml =
                        let handle = Ptree.get_node_data a tree
                        and parent = Ptree.get_node_data b tree in
                        Node.to_formatter_subtree diag_report_type
                            (pre, fi) [] tree.Ptree.data (get_unadjusted (-1) root, sa) a
                            (a, get_unadjusted b handle) (b, get_unadjusted a parent) None
                    in
                    recost, (merger a b froot),
                        [Xml.Trees.cost, `Float r.Ptree.component_cost]
                | Some ((`Single a), root) ->
                    let c1 : Xml.xml =
                        let nd = splitter (-1) root in
                        subtree_to_formatter diag_report_type (pre, fi) a a nd
                    in
                    recost, (`Single c1),
                        [Xml.Trees.cost, `Float r.Ptree.component_cost]
                | None -> assert false
            in
            recost,
                (((PXML -[Xml.Trees.tree] ([attr]) { contents }--)) :: trees)
        in
        let recost, trees =
            IntSet.fold
                (handle_to_formatter diag_report_type (pre_ref_codes, fi_ref_codes))
                (Ptree.get_handles tree)
                (0., [])
        in
        let cost = Ptree.get_cost `Adjusted tree in
        if debug_diagnosis then
            Printf.printf "end of alldirchar.to_formatter: recost = %f, cost = %f \n%!"
        recost cost;
        (RXML -[Xml.Trees.forest]
            ([Xml.Trees.recost] = [`Float recost])
            ([Xml.Trees.cost] = [`Float cost])
            ([atr])
            { set trees } --)

end
