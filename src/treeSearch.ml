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

(** [TreeSearch] contains high-level functions to perform tree searches *) 
let () = SadmanOutput.register "TreeSearch" "$Revision: 2592 $"

let has_something something (`LocalOptimum (cost_calculation)) =
    let cost_calculation = cost_calculation.Methods.cc in
    List.exists (fun x -> x = something) cost_calculation

    (* This needs to be moved out, so that the module is constructed accoding to
    * this decition 
    let has_exact = has_exact cost_calculation in
    *)

module type S = sig
      type a
      type b

      val report_trees :
          Methods.information_contained list -> string option ->  
        Data.d ->
        (a, b) Ptree.p_tree Sexpr.t -> unit

      val forest_break_search_tree :
        float ->
        (a, b) Ptree.p_tree ->
        (a, b) Ptree.p_tree

      val diagnose :
        (a, b) Ptree.p_tree ->
        (a, b) Ptree.p_tree

      val find_local_optimum :
          ?base_sampler:(a, b) Sampler.search_manager_sampler ->  ?queue : (float array * (int * int) list * int * Status.status *
          int ref * float) ->
          Data.d ->
              Sampler.ft_queue ->
                  (a, b) Ptree.p_tree Sexpr.t ->
                      All_sets.IntSet.t Lazy.t ->
                          Methods.local_optimum -> (a, b) Ptree.p_tree Sexpr.t

      val forest_search :
          Data.d ->
              Sampler.ft_queue ->
                  float ->
                      Methods.local_optimum ->
                          (a, b) Ptree.p_tree
                          Sexpr.t ->
                              (a, b) Ptree.p_tree Sexpr.t

      val fusing :
          Data.d ->
              Sampler.ft_queue ->
                  (a, b) Ptree.p_tree Sexpr.t ->
                      int option * int option * Methods.tree_weights * 
                      'a * Methods.local_optimum *
                      (int * int) -> (a, b) Ptree.p_tree Sexpr.t

        val output_consensus : Data.d -> (a, b) Ptree.p_tree Sexpr.t -> string
        option -> float option -> bool -> unit
    end

let get_transformations (`LocalOptimum (l_opt)) = 
    let clist = l_opt.Methods.cc in
    let remover meth acc =
        match meth with
        | #Methods.transform as x -> x :: acc
    in 
    List.fold_right remover clist []

let sets_of_consensus trees  = 
    Lazy.lazy_from_fun 
    (fun () ->
        let len = Sexpr.length trees 
        and trees = 
            Sexpr.fold_left 
            (fun acc x -> x.Ptree.tree :: acc)
            []
            trees 
        in
        let counters = 
            List.fold_left
            (Ptree.add_tree_to_counters (fun _ _ -> false))
            Tree.CladeFPMap.empty
            trees
        in
        Tree.CladeFPMap.fold 
        (fun set cnt acc ->
            if cnt = len then All_sets.IntSet.add set acc
            else acc)
        counters
        All_sets.IntSet.empty)

let sets_of_parser data tree =
    let get_code x = Data.taxon_code x data in
    let rec process tree acc =
        match tree with
        | Parser.Tree.Leaf x ->
                let res = All_sets.Integers.singleton (get_code x) in
                All_sets.IntSet.add res acc, res
        | Parser.Tree.Node (chld, _) ->
                let acc, lst = List.fold_left (fun (acc, res) x ->
                    let a, b = process x acc in
                    a, (b :: res)) (acc, []) chld
                in
                let union =
                    List.fold_left All_sets.Integers.union
                    All_sets.Integers.empty lst
                in
                All_sets.IntSet.add union acc, union
    in
    let sets, _ = process tree All_sets.IntSet.empty in
    sets

let sets meth data trees = 
    let `LocalOptimum (meth) = meth in
    match meth.Methods.tabu_join with
    | `Partition options ->
        (match 
            List.fold_left (fun acc x -> 
                match x with
                | `ConstraintFile file -> Some file
                | _ -> acc) None options
        with
        | None -> sets_of_consensus trees
        | Some filename ->
            try match Parser.Tree.of_file filename with
            | [[tree]] -> lazy (sets_of_parser data tree)
            | _ -> 
                    Status.user_message Status.Error
                    ("To@ use@ constraint@ files@ you@ must@ provide@ a@ " ^
                    "single@ tree,@ not@ more,@ no@ forests@ are@ allowed.");
                    failwith "Illegal input file"
            with
            | err ->
                    Status.user_message Status.Error
                    "Error@ reading@ constraint@ file";
                    raise err)
    | _ -> lazy (All_sets.IntSet.empty)

let search_time_and_trees_considered a b = 
    [ ("search-time", string_of_float a) ; 
    ("trees-considered", string_of_int b)]

module MakeNormal
    (Node : NodeSig.S) 
    (Edge : Edge.EdgeSig with type n = Node.n) 
    (TreeOps : 
            Ptree.Tree_Operations with 
            type a = Node.n with type b = Edge.e) = struct

    module PtreeSearch = Ptree.Search (Node) (Edge) (TreeOps)
    module SamplerApp = Sampler.MakeApp (Node) (Edge)
    module SamplerRes = Sampler.MakeRes (Node) (Edge)
    module PhyloQueues = Queues.Make (Node) (Edge)
    module PhyloTabus = Tabus.Make (Node) (Edge)

    type a = Node.n
    type b = Edge.e


    let debug_forest = false
    let trace_forest = true
    (** Whether to output forest search actions to the output log (top window) *)

    let odebug = Status.user_message Status.Information

    let simplified_report_trees filename data (tree, cost, _) =
        let fo = Status.Output (filename, false, []) in
        let output tree = 
            let cost = string_of_float cost in
            let tree = 
                PtreeSearch.build_trees 
                tree
                (fun x -> Data.code_taxon x data) 
                (fun _ _ -> false)
                ""
            in
            let output tree =
                Status.user_message fo "@[";
                Status.user_message fo (AsciiTree.for_formatter false true tree);
                Status.user_message fo ("[" ^ cost ^ "]");
                Status.user_message fo "@]@," in
            List.iter output tree
        in
        Status.user_message fo "@[<v>";
        output tree;
        Status.user_message fo "@]";
        Status.user_message fo "%!"

    let report_trees ic filename data trees =
        let leafsonly = 
            not (List.exists (function `Cost -> true | _ -> false) ic)
        in
        let use_hennig_style =
            List.exists (function `HennigStyle -> true | _ -> false) ic
        in
        let report_tree_len = 
            List.exists (function `Total -> true | _ -> false) ic
        in
        let newick = 
            (List.exists (function `Newick -> true | _ -> false) ic)
        in
        let collapse =
            List.exists (function `Collapse x -> x | _ -> false) ic
        in
        let ori_margin = StatusCommon.Files.get_margin filename in 
        let fo_ls = ref [] in 
        if (List.exists (function `Margin _ -> true | _ -> false) ic) then begin
            let margin = List.find (function `Margin _ -> true | _ -> false) ic in 
            match margin with
            | `Margin m -> 
                  fo_ls := (StatusCommon.Margin m)::!fo_ls;
                  StatusCommon.Files.set_margin filename m 
            | _ -> failwith "This never happens"
        end; 
        let fo = Status.Output (filename, false, !fo_ls) in
        let is_first = ref true in
        let output tree = 
            let cost = Ptree.get_cost `Adjusted tree in
            let cost = string_of_float cost in
            let tree = 
                PtreeSearch.build_forest_with_names_n_costs 
                collapse tree data cost 
            in
            let output tree =
                if use_hennig_style && not !is_first then 
                    Status.user_message fo "@,*@,"
                else is_first := false;
                Status.user_message fo "@[";
                Status.user_message fo (AsciiTree.for_formatter newick leafsonly tree);
                if leafsonly && report_tree_len then
                    Status.user_message fo ("[" ^ cost ^ "]");
                if not use_hennig_style then
                    Status.user_message fo ";";
                Status.user_message fo "@]@\n" in
            List.iter output tree
        in
        Status.user_message fo "@[<v>";
        if use_hennig_style then Status.user_message fo "tread ";
        Sexpr.leaf_iter (output) trees;
        if use_hennig_style then Status.user_message fo ";";
        Status.user_message fo "@]";
        Status.user_message (Status.Output (filename, false, !fo_ls)) "%!";
        StatusCommon.Files.set_margin filename ori_margin 
              

    let get_search_function tabu_creator trajectory meth =
        let stepfn = function
            | `Spr ->  PtreeSearch.spr_step, "SPR"
            | `Tbr -> PtreeSearch.tbr_step, "TBR"
        in
            match trajectory with
            | `AllThenChoose ->
                    (match meth with
                    | `Alternate _ -> 
                            PtreeSearch.alternate
                            (PtreeSearch.repeat_until_no_more tabu_creator 
                            (PtreeSearch.search (stepfn `Spr)))
                            (PtreeSearch.search_local_next_best (stepfn `Tbr))
                    | `SingleNeighborhood x -> 
                            PtreeSearch.repeat_until_no_more tabu_creator
                            (PtreeSearch.search_local_next_best (stepfn x))
                    | `ChainNeighborhoods x -> 
                            PtreeSearch.repeat_until_no_more tabu_creator
                            (PtreeSearch.search (stepfn x))
                    | `None -> (fun a -> a))
            | _ ->
                    (match meth with
                    | `Alternate _ -> 
                            PtreeSearch.alternate
                            (PtreeSearch.search (stepfn `Spr))
                            (PtreeSearch.search_local_next_best (stepfn `Tbr))
                    | `SingleNeighborhood x -> 
                            PtreeSearch.search_local_next_best (stepfn x)
                    | `ChainNeighborhoods x -> 
                            PtreeSearch.search (stepfn x)
                    | `None -> (fun a -> a))

    (** [forest_break_search_tree origin_cost tree] attempts to break all edges
        in [tree] with length greater than [origin_cost].  It attempts to be
        liberal about breaking edges---it may break more edges than is strictly
        correct.  However, this is fine for the forest searching algorithm.

        Note that in our trees, weights are on nodes rather than edges; thus,
        what we actually do is check the median cost at each node, and then, for
        candidate nodes, we pick the best of its edges to break.
    *)
    let forest_break_search_tree origin_cost tree =
        (** [median_cost_fn a b] returns the cost of taking the median of [a] and
            [b]---needed for breaking at a median *)
        let median_cost_fn x y a b = Node.distance ~para:x ~parb:y a b in

        (** [break edge tree] breaks an edge in the tree and updates the tree data
            accordingly *)
        let break edge tree =
            let Tree.Edge(bfrom, bto) = edge in
            let tree, break_delta, float, int, data, incr =
                TreeOps.break_fn (bfrom, bto) tree in
            let tree = TreeOps.uppass tree in
            tree in

        let tree = Ptree.set_origin_cost origin_cost tree in

        (* Iterate over all the nodes, possibly breaking... *)
        let tree = All_sets.IntegerMap.fold
            (fun id _ tree ->
                 (* This is where we decide whether or not to break at a given
                    node *)
(*                 try*)
(*                    let parent = Ptree.get_parent id tree*)
(*                    and median = Ptree.get_node_data id tree in*)
(*                    let cost = Node.node_cost (Some parent) median in*)
(*                 if cost >= origin_cost && not (Tree.is_leaf id tree.Ptree.tree)*)
(*                 then begin*)
(*                         if trace_forest*)
(*                         then odebug ("Forest search: Breaking at node"*)
(*                                      ^ " with cost "*)
(*                                      ^ string_of_float cost)*)
(*                         else ();*)
                         (* To add more breaks: this is the place *)
                         try if (not (Tree.is_leaf id tree.Ptree.tree)) then
                             let new_tree = 
                                 break
                                 (Ptree.break_median_edge
                                 median_cost_fn
                                 id tree) tree
                             in
                             let new_cost = Ptree.get_cost `Adjusted new_tree
                             and old_cost = Ptree.get_cost `Adjusted tree in
                             Status.user_message Status.Information
                             ("The break new old is " ^ string_of_float new_cost
                             ^ " - " ^ string_of_float old_cost);
                             if new_cost < old_cost then new_tree
                             else tree
                        else tree 
                         with Tree.Invalid_Node_Id _ -> tree )
(*                end else*)
(*                     (if debug_forest*)
(*                      then odebug ("Not breaking at node "*)
(*                                   ^ string_of_int id*)
(*                                   ^ " with cost "*)
(*                                   ^ string_of_float cost);*)
(*                            tree)*)
(*                     with Tree.Invalid_Node_Id id ->*)
(*                         (if debug_forest*)
(*                          then odebug ("Node doesn't exist in tree: "*)
(*                                       ^ string_of_int id);*)
(*                          tree))*)
            tree.Ptree.node_data tree in
        (* (TODO: also check the roots for breaking!!) *)

        tree

    (** [forest_joins forest] attempts to join pairs of forest components using
        TBR join.  Those whose join cost is less than the origin/loss cost will
        be kept. *)
    let rec forest_joins forest =
        let components = Ptree.components forest in
        let join_tabu = PhyloTabus.join_to_tree_in_forest forest in

        let status = Status.create "Attempting to join forest components"
            (Some components) "" in

        (** [tbr_joins component] tries to join [component] to all other
            components in the tree *)
        let tbr_joins component =
            Status.full_report ~adv:component status;
            let tabu, right = join_tabu component in
            let mgr = new PhyloQueues.first_best_srch_mgr (new Sampler.do_nothing) in
            let () =
                mgr#init [(forest, Ptree.get_cost `Adjusted forest,
                           Ptree.NoCost, tabu)] in

            (* get the current `Right junction *)
            let j2 = Ptree.jxn_of_handle forest right in
            
            (* get the `Right (clade) root node *)
            let clade_node =
                let root =
                    All_sets.IntegerMap.find right forest.Ptree.component_root in
                let root = root.Ptree.root_median in
                match root with
                    (* This shouldn't happen;  our passed-in forest should have
                       been evaluated, etc., and all the roots should be 
                       defiend. *)
                | None -> assert false
                | Some (_, clade_node) -> clade_node in

            let status =
                PtreeSearch.tbr_join mgr tabu forest j2 clade_node
                    forest.Ptree.origin_cost in

            match status with
            | Tree.Break ->
                  let results = mgr#results in
                  (* Only interested in the first one *)
                  let (forest, _, _) = List.hd results in
                  Some forest
            | Tree.Continue
            | Tree.Skip -> None 
        in
        let rec try_comp component =
            if component = components
            then None
            else match tbr_joins component with
            | None -> try_comp (succ component)
            | Some forest -> Some forest in

        let res = try_comp 0 in
        Status.finished status;

        match res with
        | None -> forest
        | Some forest -> forest_joins forest

    let diagnose tree =
        PtreeSearch.uppass (PtreeSearch.downpass tree)

    let queue_manager max th keep sampler =
        fun () ->
            match max, th with
            | 1, _ -> new PhyloQueues.first_best_srch_mgr (sampler ())
            | n, th -> new PhyloQueues.hold_n_threshold_srch_mgr n keep th (sampler ())

    let create_sampler data queue previous item = 
        let ob = 
            match item with
            | `PrintTrajectory filename -> 
                  (new SamplerApp.print_next_tree 
                  (simplified_report_trees filename data))
            | `KeepBestTrees ->
                    (new SamplerApp.local_optimum_holder queue)
            | `TimeOut time ->
                    (new SamplerApp.timed_cancellation time) 
            | `TimedPrint (time, filename) ->
                    (new SamplerApp.timed_printout queue time 
                     (simplified_report_trees filename data))
            | `UnionStats (filename, depth) ->
                    new SamplerRes.union_table depth
                    (Status.user_message (Status.Output (filename, false, [])))
            | `RootUnionDistr filename ->
                    Status.user_message Status.Error 
                    ("Sorry@ the@ root@ union@ distribution@ sampler@ is@ "
                    ^ "currently@ unsupported");
                    new Sampler.do_nothing 
                    (* TODO: This sampler is off for now *)
(*                    new SamplerRes.union_root_distribution *)
(*                    (Status.user_message (Status.Output (filename, false)))*)
            | `AttemptsDistr filename ->
                    new SamplerRes.tests_before_next 
                    (Status.user_message (Status.Output (filename, false, [])))
            | `BreakVsJoin filename ->
                    new SamplerRes.break_n_join_distances 
                    TreeOps.join_fn 
                    (Status.user_message (Status.Output (filename, false, [])))
            | `AllVisited filename ->
                    let join_fn a b c = 
                        let a, _ = TreeOps.join_fn [] a b c in
                        a
                    in
                    (new SamplerApp.visited join_fn 
                    (simplified_report_trees filename data))
        in
        new Sampler.composer previous ob

    let sampler sampler data queue lst () =
        let sampler = 
            match sampler with
            | Some x -> x
            | None -> new Sampler.do_nothing
        in
        List.fold_left 
        (fun prev item -> create_sampler data queue prev item) 
        sampler
        lst

let rec find_local_optimum ?base_sampler ?queue data emergency_queue
        (trees : (a, b) Ptree.p_tree Sexpr.t)
        (sets :   All_sets.IntSet.t Lazy.t)
        (meth : Methods.local_optimum) :
    (a, b) Ptree.p_tree Sexpr.t =
    let trees = Sexpr.map TreeOps.unadjust trees in
    let local_search_results_reporting lst =
        let builder (acc, cnt) (_, cost) =
            let hd = 
                ("tree_" ^ string_of_int cnt ^ 
                "_cost", string_of_float cost) 
            in
            hd :: acc, cnt + 1
        in
        let acc, _ = 
            Sexpr.fold_left 
            (fun acc x -> List.fold_left builder acc x) 
            ([], 0) 
            lst 
        in
        acc
    in

    (* let `LocalOptimum
            (search_space, th, max, keep, cost_calculation, origin, 
            trajectory, break_tabu, join_tabu, reroot_tabu, nodes_tabu, samples)
            = meth in *)
    let `LocalOptimum (l_opt) = meth in
    let samplerf = sampler base_sampler data emergency_queue l_opt.Methods.samples in
    let queue_manager =
        match queue with
        | Some
            (best_vals, node_indices, starting, status, nbhood_count,
            orig_cost) -> 
                (fun () -> new PhyloQueues.supports_manager best_vals node_indices 
                starting status nbhood_count orig_cost (new Sampler.do_nothing))
        | None ->
                match l_opt.Methods.tm with (* trajectory *)
                | `AllAround f -> 
                        (fun () -> new PhyloQueues.all_possible_joins f
                        (samplerf ()))
                | `BestFirst -> 
                            queue_manager l_opt.Methods.num_keep
                            l_opt.Methods.threshold l_opt.Methods.keep samplerf
                | `AllThenChoose -> 
                        fun () -> 
                            new PhyloQueues.all_neighbors_srch_mgr
                            TreeOps.join_fn TreeOps.break_fn TreeOps.reroot_fn 
                            TreeOps.incremental_uppass
                            (samplerf ())
                | `Annealing (a, b) -> 
                        fun () ->
                            new PhyloQueues.annealing_srch_mgr 1 `Last 0.0 a b
                            (samplerf ())
                | `PoyDrifting (a, b) -> 
                        fun () ->
                            new PhyloQueues.classic_poy_drifting_srch_mgr 1 `Last
                            a b (samplerf ())
    in
    let tabu_manager ptree = 
        let get_depth = function
            | None -> max_int 
            | Some v -> v
        in
        let breakfn =
            match l_opt.Methods.tabu_break with
            | `Randomized -> PhyloTabus.random_break
            | `DistanceSorted -> PhyloTabus.sorted_break
            | `OnlyOnce -> PhyloTabus.only_once_break
        in
        let joinfn =
            match l_opt.Methods.tabu_join with
            | `UnionBased depth ->
                    PhyloTabus.union_join (get_depth depth)
            | `AllBased depth -> 
                    PhyloTabus.distance_join (get_depth depth)
            | `Partition options ->
                    let depth = 
                        List.fold_left (fun acc x ->
                            match x with
                            | `MaxDepth x -> Some x
                            | _ -> acc)
                        None options
                    in
                    PhyloTabus.partitioned_join (Lazy.force sets) (get_depth depth)
        in
        let iterfn =
            match l_opt.Methods.tabu_nodes with
            | `Null -> PhyloTabus.simple_nm_none 
            | `Leaves -> PhyloTabus.simple_nm_leaves
            | `Randomized -> PhyloTabus.random_nm
            | `Longest -> PhyloTabus.sorted_longest_nm
            | `Averaged ->PhyloTabus.sorted_average_nm
        in
        let rerootfn =
            match l_opt.Methods.tabu_reroot with
            | `Bfs depth ->
                    PhyloTabus.reroot (get_depth depth)
        in
        new PhyloTabus.standard_tabu ptree joinfn rerootfn breakfn iterfn
    in
    let search_fn = get_search_function tabu_manager l_opt.Methods.tm l_opt.Methods.ss in
    let search_features =
            PtreeSearch.features meth ((queue_manager ())#features
            (* TODO: tabu features *)
            (([])))
    in
    Sadman.start "search" search_features;
    let timer = Timer.start () in
    let result = 
            PhyloQueues.reset_trees_considered ();
            let process_tree tree = 
                let cost = Ptree.get_cost `Adjusted tree in
                let tabu_manager = tabu_manager tree in
                let queue_manager = queue_manager () in
                queue_manager#init [(tree, cost, Ptree.NoCost, tabu_manager)];
                try
                    let res = (search_fn queue_manager)#results in
                    List.map (fun (a, _, _) ->
                        let a = PtreeSearch.uppass a in
                        (a, Ptree.get_cost `Adjusted a)) res
                with
                | Sampler.TimedOut -> [(tree, Ptree.get_cost `Adjusted tree)]
            in
            Sexpr.map_status "Tree search" process_tree trees 
    in
    let time = Timer.get_user timer in
    let trees_considered = PhyloQueues.get_trees_considered () in
    let sadman_results = local_search_results_reporting result in
    Sadman.finish 
        (search_time_and_trees_considered time trees_considered @
            sadman_results);
    Sexpr.map_insexpr (List.map (fun (a, _) -> `Single a)) result


let forest_search data queue origin_cost search trees =
    (* Forest search happens in three steps.  First, we find candidate high-cost
       medians for breaking.  Once those are broken, we search on the entire
       forest of trees.  We then attempt to TBR join to improve our overall
       score... *)

    let trees =
        Sexpr.map_status "Breaking trees"
            (forest_break_search_tree origin_cost)
            trees in

    let trees = find_local_optimum data queue trees (sets search data trees) search in

    (* TBR joins *)
    let trees = Sexpr.map_status "TBR joining trees" forest_joins trees in

    Sexpr.leaf_iter
        (fun trees ->
    let cost = Ptree.get_cost `Adjusted trees in
    Status.user_message Status.Information ("Forest search found forest with "
                                            ^ string_of_int (Ptree.components
                                                                 trees)
                                            ^ " components and cost "
                                            ^ string_of_float cost)
        ) trees;

    let trees = Sexpr.map
        (Ptree.set_origin_cost 0.) trees in

    trees

    (** [fusing trees params] performs tree fusing with the given parameters *)
    let fusing data queue trees 
            (iterations, max_trees, weighting, keep_method, local_optimum, (min,
                                                                            max)) =
        (* default max_trees to number of trees *)
        let max_trees = match max_trees with
        | Some m -> m
        | None -> Sexpr.cardinal trees in
        (* default iterations to....... 4 * max_trees ?? *)
        let iterations = match iterations with
        | Some i -> i
        | None -> max_trees * 4 in

        Sadman.start "tree-fusing"
            [("iterations", string_of_int iterations);
             ("max_trees", string_of_int max_trees);
            ];

        let weighting = match weighting with
        | `Uniform -> (fun _ -> 1.)
        in
        let keep_method = `Best in
        (* TODO: take only best of these *)
        let process t = Sexpr.to_list (find_local_optimum data queue
        (Sexpr.singleton t) (sets_of_consensus trees) local_optimum) in
        let trees = 
            PtreeSearch.fuse_generations
            (Sexpr.to_list trees)
            max_trees
            weighting
            keep_method
            iterations
            process
            (min, max)
        in
        Sadman.finish [];
        Sexpr.of_list trees

    let output_consensus data trees filename v graphic = 
        let _ = incr Data.median_code_count in
        (* We will use a negative number for the root code to avoid
        * any clash *)
        let ntrees = Sexpr.length trees in
        let majority, majority_text = 
            match v with
            | None -> ntrees, "Strict"
            | Some v -> 
                    if v > 100.0 then begin
                        Status.user_message Status.Error 
                        ("You@ have@ requested@ a@ consensus@ " ^
                        "with@ majority@ rule@ of@ more@ than@ " ^
                        "100@ percent,@ I@ don't@ see@ how@ to@ " ^
                        "do@ that.");
                        failwith "Illegal Consensus Parameter"
                    end else if v <= 50.0 then begin
                        Status.user_message Status.Error 
                        ("You@ have@ requested@ a@ consensus@ " ^
                        "with@ majority@ rule@ less@ than@ or@ equal@ to" ^
                        "50@ percent;@ I@ can@ not@ do@ that@ because@ " ^
                        "that@ percentage@ is@ not@ a@ majority!.@ " ^
                        "Either@ you@ made@ a@ typo@ " ^
                        "or@ you@ should@ reconsider@ your@ " ^
                        "parameter@ selection.");
                        failwith "Illegal Consensus Parameter";
                    end else 
                        int_of_float 
                        (ceil ((v *. (float_of_int ntrees) /.
                        100.0))), ((string_of_float v) ^ " percent")
        in
        let lst_trees = Sexpr.to_list trees in
        let rooting_leaf, tree = 
            match data.Data.root_at with
            | Some v -> 
                    (match lst_trees with
                    | hd :: _ -> v, hd
                    | _ -> failwith "No trees in memory")
            | None -> 
                    match lst_trees with
                    | hd :: _ -> Ptree.choose_leaf hd, hd
                    | _ -> failwith "No trees in memory"
        in
        let res = Ptree.consensus PtreeSearch.collapse_as_needed (fun code ->
            Data.code_taxon code data) majority 
            (Sexpr.to_list trees) rooting_leaf
        in
        let fo = Status.Output (filename, false, []) in
        if not graphic then begin
            Status.user_message fo "@[<v>";
            Status.user_message fo 
            ("@[" ^ majority_text ^ "@ Majority@ Consensus Tree@]@,@[");
            Status.user_message fo (AsciiTree.for_formatter false false res);
            Status.user_message fo "@]@]\n%!";
        end else
            match filename with
            | Some filename -> 
                    let title = majority_text ^ " Majority Consensus Tree" in
                    GraphicsPs.display title filename [|(0.0, res)|]
            | None -> 
                    let r = AsciiTree.to_string ~sep:2 ~bd:2 true res in
                    Status.user_message fo 
                    ("@[@[<v>@[" ^ majority_text ^ "@ Majority@ Consensus Tree@]@,@[");
                    Status.user_message (Status.Output (filename,false, [])) r;
                    Status.user_message (Status.Output (filename, false, [])) "@]@]@]%!";

end

module Make
    (NodeH : NodeSig.S) 
    (EdgeH : Edge.EdgeSig with type n = NodeH.n) 
    (TreeOpsH : 
            Ptree.Tree_Operations with 
            type a = NodeH.n with type b = EdgeH.e)
     = struct

    module NodeS = Node.Standard
    module EdgeS = Edge.SelfEdge
    module TreeOpsS = Chartree.TreeOps

    type a = NodeH.n
    type b = EdgeH.e

    module SH = MakeNormal (NodeS) (EdgeS) (TreeOpsS)
    module DH = MakeNormal (NodeH) (EdgeH) (TreeOpsH)

    module TOS = TreeOpsS 
    module TOH = TreeOpsH 

    let report_trees = DH.report_trees
    let forest_break_search_tree = DH.forest_break_search_tree
    let diagnose = DH.diagnose

    let replace_contents downpass uppass get_code nodes ptree =
        let nt = { Ptree.empty with Ptree.tree = ptree.Ptree.tree } in
        uppass (downpass 
        (List.fold_left (fun nt node ->
            Ptree.add_node_data (get_code node) node nt) 
        nt nodes))

    let from_s_to_h = replace_contents TOH.downpass TOH.uppass NodeH.taxon_code 
    let from_h_to_s = replace_contents TOS.downpass TOS.uppass NodeS.taxon_code

    let find_local_optimum ?base_sampler ?queue data b trees d e =
        let sampler = 
            match base_sampler with
            | None -> new Sampler.do_nothing
            | Some x -> x 
        in
        match  Data.has_likelihood data || Data.has_dynamic data, queue with
        | true, None -> DH.find_local_optimum ~base_sampler:sampler data b trees d e
        | true, Some queue -> 
                DH.find_local_optimum ~base_sampler:sampler ~queue data b trees d e
        | false, queue ->
                let data, nodes = NodeS.load_data data in
                let trees = Sexpr.map (from_h_to_s nodes) trees in
                let trees = 
                    match queue with
                    | None -> SH.find_local_optimum data b trees d e
                    | Some queue -> 
                            SH.find_local_optimum ~queue data b trees d e
                in
                let data, nodes = NodeH.load_data data in
                Sexpr.map (from_s_to_h nodes) trees

    let forest_search data b c d trees =
        if Data.has_dynamic data then
            DH.forest_search data b c d trees 
        else
            let data, nodes = NodeS.load_data data in
            let trees = Sexpr.map (from_h_to_s nodes) trees in
            let trees = SH.forest_search data b c d trees in
            let data, nodes = NodeH.load_data data  in
            Sexpr.map (from_s_to_h nodes) trees

    let fusing data a trees b =
        if Data.has_dynamic data then
            DH.fusing data a trees b
        else
            let data, nodes = NodeS.load_data data in
            let trees = Sexpr.map (from_h_to_s nodes) trees in
            let trees = SH.fusing data a trees b in
            let data, nodes = NodeH.load_data data in
            Sexpr.map (from_s_to_h nodes) trees


    let output_consensus = DH.output_consensus
end
