(* POY 4.0 Beta. A phylogenetic analysis program using Dynamic Homologies.    *)
(* Copyright (C) 2007  Andr�s Var�n, Le Sy Vinh, Illya Bomash, Ward Wheeler,  *)
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
let () = SadmanOutput.register "TreeSearch" "$Revision: 2871 $"

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
          ?base_sampler:(a, b) Sampler.search_manager_sampler -> 
              ?queue : (float array * (int * int) list * int * Status.status * int ref * float) ->
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
        | Tree.Parse.Leafp x ->
                let res = All_sets.Integers.singleton (get_code x) in
                All_sets.IntSet.add res acc, res
        | Tree.Parse.Nodep (chld, _) ->
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
    fst (process (Tree.Parse.strip_tree tree) All_sets.IntSet.empty)

let sets meth data trees = 
    match meth with
    | `Partition options ->
        (match 
            List.fold_left (fun acc x -> 
                match x with
                | `ConstraintFile file -> Some file
                | _ -> acc) None options
        with
        | None -> sets_of_consensus trees
        | Some filename ->
            try match Tree.Parse.of_file filename with
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
    module SamplerRes = Sampler.MakeRes (Node) (Edge) (TreeOps)
    module PhyloQueues = Queues.Make (Node) (Edge)
    module PhyloTabus = Tabus.Make (Node) (Edge)

    type a = Node.n
    type b = Edge.e


    let debug_forest = false
    let trace_forest = true
    (** Whether to output forest search actions to the output log (top window) *)

    let odebug = Status.user_message Status.Information

    let simplified_report_trees compress filename data (tree, cost, _) =
        let fo = 
            let lst = if compress then [StatusCommon.Compress] else [] in
            Status.Output (filename, false, lst) 
        in
        let output tree = 
            let cost = string_of_float cost in
            let tree = 
                PtreeSearch.build_trees 
                tree
                (fun x -> Data.code_taxon x data) 
                (fun _ _ -> false)
                None
                ""
            in
            let output tree =
                Status.user_message fo "@[";
                Status.user_message fo 
                (AsciiTree.for_formatter false true true tree);
                Status.user_message fo ("[" ^ cost ^ "]");
                Status.user_message fo "@]@," in
            List.iter output tree
        in
        Status.user_message fo "@[<v>";
        output tree;
        Status.user_message fo "@]";
        Status.user_message fo "%!"

    let report_trees_and_branches compress filename data branches ptree : unit =
        let fo = 
            let lst = if compress then [StatusCommon.Compress] else [] in
            Status.Output (filename, false, lst)

        and trees = 
            PtreeSearch.build_trees
                ptree.Ptree.tree
                (fun x -> Data.code_taxon x data)
                (fun _ _ -> false) (* don't collapse *)
                (Some branches)
                ""
        in
        let adj_cost = string_of_float (Ptree.get_cost `Adjusted ptree)
        and unadj_cost = string_of_float (Ptree.get_cost `Unadjusted ptree) in
        Status.user_message fo "@[<v>";
        List.iter 
            (fun x -> 
                Status.user_message fo "@[";
                Status.user_message fo (AsciiTree.for_formatter false true true x);
                Status.user_message fo ("[" ^ adj_cost ^ "|" ^ unadj_cost ^ "]");
                Status.user_message fo "@]@,"
            ) trees;
        Status.user_message fo "@]%!"

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
        let newline = if use_hennig_style then "" else "@\n" in
        let ic = 
            if use_hennig_style then (`Margin (1000000010 - 1)) :: ic
            else ic
        in
        let branches = List.exists (function `Branches -> true | _ -> false) ic in
        (*
        let newick = 
            (List.exists (function `Newick -> true | _ -> false) ic)
        in
        *)
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
                collapse tree cost branches
            in
            let output tree =
                if use_hennig_style && not !is_first then 
                    Status.user_message fo " * "
                else is_first := false;
                Status.user_message fo "@[";
                Status.user_message fo 
                (AsciiTree.for_formatter (not use_hennig_style ) 
                (not use_hennig_style) leafsonly tree);
                if leafsonly && report_tree_len then
                    Status.user_message fo ("[" ^ cost ^ "]");
                if not use_hennig_style then Status.user_message fo ";"
                else Status.user_message fo "@?";
                Status.user_message fo "@]";
                Status.user_message fo newline;
            in
            List.iter output tree
        in
        Status.user_message fo (if use_hennig_style then "@[<h>" else "@[<v>");
        if use_hennig_style then Status.user_message fo "tread ";
        Sexpr.leaf_iter (output) trees;
        if use_hennig_style then Status.user_message fo ";";
        Status.user_message fo (if use_hennig_style then  "@]" else "@]@\n" );
        Status.user_message (Status.Output (filename, false, !fo_ls)) "@\n%!";
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
                            (PtreeSearch.search true (stepfn `Spr)))
                            (PtreeSearch.search_local_next_best (stepfn `Tbr))
                    | `SingleNeighborhood x -> 
                            PtreeSearch.repeat_until_no_more tabu_creator
                            (PtreeSearch.search_local_next_best (stepfn x))
                    | `ChainNeighborhoods x -> 
                            PtreeSearch.repeat_until_no_more tabu_creator
                            (PtreeSearch.search false (stepfn x))
                    | `None -> (fun a -> a))
            | _ ->
                    (match meth with
                    | `Alternate _ -> 
                            PtreeSearch.alternate
                            (PtreeSearch.search true (stepfn `Spr))
                            (PtreeSearch.search_local_next_best (stepfn `Tbr))
                    | `SingleNeighborhood x -> 
                            PtreeSearch.search_local_next_best (stepfn x)
                    | `ChainNeighborhoods x -> 
                            PtreeSearch.search false (stepfn x)
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
        (** [break edge tree] breaks an edge in the tree and updates the tree data
            accordingly *)
        let break edge tree =
            let Tree.Edge(bfrom, bto) = edge in
            let breakage = TreeOps.break_fn (bfrom, bto) tree in
            TreeOps.uppass breakage.Ptree.ptree 
        in
        let tree = Ptree.set_origin_cost origin_cost tree in
        (* Iterate over all the edges, possibly breaking... *)
        let tree = 
            Tree.EdgeMap.fold
            (fun edge _ tree ->
                try
                    let cost = Ptree.get_cost `Adjusted tree in
                    let new_tree = break edge tree in
                    let new_cost = Ptree.get_cost `Adjusted new_tree in
                    if new_cost < cost then new_tree else tree
                with _ -> tree)
            tree.Ptree.edge_data tree in
        (* (TODO: also check the roots for breaking!!) *)

        tree

    (** [forest_joins forest] attempts to join pairs of forest components using
        TBR join.  Those whose join cost is less than the origin/loss cost will
        be kept. *)
    let rec forest_joins forest =
        let components = Ptree.components forest in
        (*
        let join_tabu = PhyloTabus.join_to_tree_in_forest forest in
        *)
        let status = Status.create "Attempting to join forest components"
            (Some components) "" in

        (** [tbr_joins component] tries to join [component] to all other
            components in the tree *)
        let tbr_joins component = 
            failwith "Forest searches are off in this release"
            (* TODO:
                * To get the forest search working again, we need to modify the
                * tabu managers so that instead of using left and right uses a
                * code assigned to each individual component. That's the only
                * way to get the necessary way to connect multiple elements in a
                * forest using the tabu managers. Right now there is no nice way
                * to do it, and this is a low priority issue, therefore, I am
                * leaving a note and doing it later. *)
            (*
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
            *)
        in
        let rec try_comp component =
            if component = components
            then None
            else match tbr_joins component with
            | None -> try_comp (succ component)
            | Some forest -> Some forest in

        let res = None in
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
                  (simplified_report_trees false filename data))
            | `KeepBestTrees ->
                    (new SamplerApp.local_optimum_holder queue)
            | `MaxTreesEvaluated trees ->
                    (new SamplerApp.counted_cancellation trees) 
            | `TimeOut time ->
                    (new SamplerApp.timed_cancellation time) 
            | `TimedPrint (time, filename) ->
                    (new SamplerApp.timed_printout queue time 
                     (simplified_report_trees false filename data))
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
                    (Status.user_message (Status.Output (filename, false, [])))
            | `Likelihood filename ->
                    let do_compress = None <> filename in
                    new SamplerRes.likelihood_verification
                        (Status.user_message (Status.Output (filename, false, [])))
                        (report_trees_and_branches do_compress filename data)
            | `AllVisited filename ->
                    let join_fn incr a b c = 
                        let a, _ = TreeOps.join_fn incr a b c in
                        a
                    in
                    let do_compress = None <> filename in
                    (new SamplerApp.visited join_fn 
                    (simplified_report_trees do_compress filename data))
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
    let sets =
        try
            match l_opt.Methods.tabu_join with
            | `Partition opts -> 
                    (match 
                        List.find (function `Sets _ -> true | _ -> false)
                        opts
                    with
                    | `Sets x -> x
                    | _ -> assert false)
            | _ -> sets
        with
        | Not_found -> sets
    in
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
    let partition_for_other_tabus =
        match l_opt.Methods.tabu_join with
        | `Partition _ -> 
                Some (`Sets (Lazy.force sets))
                (* TMP
                Some (`Height 2)
                *)
        | _ -> None
    in
    let tabu_manager ptree = 
        let get_depth = function
            | None -> max_int 
            | Some v -> v
        in
        let breakfn =
            match l_opt.Methods.tabu_break with
            | `Randomized -> PhyloTabus.random_break
            | `DistanceSorted early -> 
                    PhyloTabus.sorted_break early partition_for_other_tabus
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
                    PhyloTabus.partitioned_join (`Sets (Lazy.force sets)) 
                    (get_depth depth)
                    (* TMP 
                    PhyloTabus.partitioned_join (`Height 4) (get_depth depth)
                    *)
        in
        let iterfn =
            match l_opt.Methods.tabu_nodes with
            | `Null -> PhyloTabus.simple_nm_none 
            | `Leaves -> PhyloTabus.simple_nm_leaves
        in
        let rerootfn =
            match l_opt.Methods.tabu_reroot with
            | `Bfs depth ->
                    PhyloTabus.reroot partition_for_other_tabus (get_depth depth)
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
                        if a.Ptree.tree <> tree.Ptree.tree then
                            let a = PtreeSearch.uppass a in 
                            (a, Ptree.get_cost `Adjusted a)
                        else (a, Ptree.get_cost `Adjusted a)) res
                with
                | Methods.TimedOut -> [(tree, Ptree.get_cost `Adjusted tree)]
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

    let trees = 
        let `LocalOptimum s = search in
        find_local_optimum data queue trees 
        (sets s.Methods.tabu_join data trees) search
    in 
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
        let max_code, cntr = 
            All_sets.IntegerMap.fold (fun code _ (acc, cnt) -> 
                Pervasives.max code acc, cnt + 1)
            data.Data.taxon_codes (0, 0)
        in
        (* TODO: take only best of these *)
        let process t = Sexpr.to_list (find_local_optimum data queue
        (Sexpr.singleton t) (sets_of_consensus trees) local_optimum) in
        let trees = 
            PtreeSearch.fuse_generations
            (Sexpr.to_list trees)
            max_code
            max_trees
            weighting
            keep_method
            iterations
            process
            (min, Pervasives.min (cntr - 3) max)
        in
        Sadman.finish [];
        Sexpr.of_list trees

    let output_consensus data trees filename v graphic = 
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
            Status.user_message fo (AsciiTree.for_formatter false true false 
            res);
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
    (NodeH : NodeSig.S with type other_n = Node.Standard.n) 
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

    let collect_nodes data trees =
        let tree = Sexpr.first trees in
        let aux_collect key node taxaacc = 
            match node with
            | Tree.Single _ | Tree.Leaf _ -> 
                    (Ptree.get_node_data key tree) :: taxaacc
            | _ -> taxaacc
        in
        let nodes = All_sets.IntegerMap.fold aux_collect 
            tree.Ptree.tree.Tree.u_topo [] in
        nodes, List.map NodeH.to_other nodes

    let replace_contents downpass uppass get_code nodes ptree =
        let nt = { (Ptree.empty ptree.Ptree.data) with Ptree.tree = ptree.Ptree.tree } in
        uppass (downpass 
        (List.fold_left (fun nt node ->
            Ptree.add_node_data (get_code node) node nt) 
        nt nodes))

    let from_s_to_h = replace_contents TOH.downpass TOH.uppass NodeH.taxon_code 
    let from_h_to_s = replace_contents TOS.downpass TOS.uppass NodeS.taxon_code

    let find_local_optimum ?base_sampler ?queue data b trees d e =
        if 0 = Sexpr.length trees then trees
        else
            let sampler =
                match base_sampler with
                | None -> new Sampler.do_nothing
                | Some x -> x 
            in
        match  Data.has_likelihood data || Data.has_dynamic data, queue with
            | true, None -> 
                    DH.find_local_optimum ~base_sampler:sampler data b trees d e
            | true, Some queue -> 
                    DH.find_local_optimum ~base_sampler:sampler ~queue data b 
                    trees d e
            | false, queue ->
                    let nodeh, nodes = collect_nodes data trees in
                    let trees = Sexpr.map (from_h_to_s nodes) trees in
                    let trees = 
                        match queue with
                        | None -> SH.find_local_optimum data b trees d e
                        | Some queue -> 
                                SH.find_local_optimum ~queue data b trees d e
                    in
                    Sexpr.map (from_s_to_h nodeh) trees

    let forest_search data b c d trees =
        if 0 = Sexpr.length trees then trees
        else
            if Data.has_dynamic data then
                DH.forest_search data b c d trees 
            else
                let nodeh, nodes = collect_nodes data trees in
                let trees = Sexpr.map (from_h_to_s nodes) trees in
                let trees = SH.forest_search data b c d trees in
                Sexpr.map (from_s_to_h nodeh) trees

    let fusing data a trees b =
        if Data.has_dynamic data then
            DH.fusing data a trees b
        else
            let nodeh, nodes = collect_nodes data trees in
            let trees = Sexpr.map (from_h_to_s nodes) trees in
            let trees = SH.fusing data a trees b in
            Sexpr.map (from_s_to_h nodeh) trees


    let output_consensus = DH.output_consensus
end
