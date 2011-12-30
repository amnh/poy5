(* POY 5.0 Alpha. A phylogenetic analysis program using Dynamic Homologies.   *)
(* Copyright (C) 2011 Andrés Varón, Lin Hong, Nicholas Lucaroni, Ward Wheeler,*)
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

(* $Id: supports.ml 2871 2008-05-23 17:48:34Z andres $ *)
(* Created Tue Jan 31 16:39:25 2006 (Illya Bomash) *)
let () = SadmanOutput.register "Support" "$Revision: 2871 $"

let infinity = float_of_int (max_int / 4)

module type S = sig
        type a 
        type b

    val support : 
        (a, b) Ptree.p_tree Sexpr.t ->
        a list -> 
        Methods.support_method -> Data.d -> 
         Sampler.ft_queue ->
             int Tree.CladeFPMap.t

    val bremer_support : 
        (a, b) Ptree.p_tree Sexpr.t ->
        int -> int ->
        a list -> (a, b) Ptree.p_tree Sexpr.t -> Methods.local_optimum
        -> Methods.build -> Data.d -> Sampler.ft_queue ->
            Methods.support_tree Sexpr.t

    val support_to_string_tree : 
        Data.d -> Methods.support_tree -> Tree.Parse.tree_types

    (** [join_support_trees trees] takes a list of [(iterations, support_tree)]
        pairs and combines them into a single support tree *)
    val join_support_trees : 
        (int * Methods.support_tree) list -> Methods.support_tree

val bremer_of_input_file :
    (Tree.u_tree -> string -> int) -> int ->
        (int -> string) -> Data.d -> Methods.filename list -> 
            (a, b) Ptree.p_tree Sexpr.t -> Tree.Parse.tree_types Sexpr.t

(** Like [bremer_of_input_file] but trust whatever input cost is provided with
* each tree .*)
val bremer_of_input_file_but_trust_input_cost : int ->
    (int -> string) -> Data.d -> Methods.filename list -> 
        (a, b) Ptree.p_tree Sexpr.t -> Tree.Parse.tree_types Sexpr.t

end
(** support.ml *)

module MakeNormal (Node : NodeSig.S with type other_n = Node.Standard.n) (Edge : Edge.EdgeSig with type n = Node.n) 
    (TreeOps : 
        Ptree.Tree_Operations with type a = Node.n with type b = Edge.e) = struct
        type a = Node.n
        type b = Edge.e

    let ndebug_bremer_best_costs = false
    let ndebug_edgevisit = false
    let ndebug_bremertbr = false
    let odebug = Status.user_message Status.Information

    type support_tree = Methods.support_tree

    module CompHashes = struct
        type t = Hash_tbl.Interface.fp
        let compare = compare
    end

    module HashCount = Map.Make(CompHashes)
    module PTS = TreeSearch.MakeNormal (Node) (Edge) (TreeOps)
    module PQueue = Queues.Make (Node) (Edge)
    module CT = CharTransform.Make (Node) (Edge) (TreeOps)
    module B = Build.MakeNormal (Node) (Edge) (TreeOps)
    module TO = Ptree.Search (Node) (Edge) 
    (TreeOps) 

    (** [combine_counts acc] combines a list of counts (maps) by adding them
        together. *)
    let rec combine_counts ?(acc=HashCount.empty) count_list = match count_list with
    | [] -> acc
    | count :: counts ->
          combine_counts 
              ~acc:(HashCount.fold
                        (fun k v map ->
                             let prev =
                                 try HashCount.find k map
                                 with Not_found -> 0 in
                             HashCount.add k (prev + v) map)
                        count
                        acc)
              counts

    let resample_support trees data queue iterations otus perturb search build
        outgroup counters =
        let n, data = 
            let cnt = ref 0 in
            let chars = ref [] in
            Hashtbl.iter (fun _ code -> 
                if Data.apply_boolean
                NonaddCS8.is_potentially_informative
                AddCS.is_potentially_informative data code then incr cnt 
                else chars := code :: !chars) 
            data.Data.character_names;
            let data = Data.process_ignore_characters false data (`Some !chars) in
            !cnt, data
        in
        let new_data = 
            match perturb with
            | `Bootstrap -> Data.Sample.generate data `Bootstrap
            | `Jackknife m -> 
                    if m > 100.0 || m < 0.0 then begin
                        Status.user_message Status.Error 
                        "Illegal@ Jackknife@ sample@ size:@ The@ sample@ \
                        is@ the@ percentage@ of@ characters@ that@ should@ \
                        be@ sampled,@ therefore@ it's@ valid@ values@ are@ \
                        between@ 0.0@ and@ 100.0";
                        failwith "Illegal Jackknife sample";
                    end;
                    Data.Sample.generate data (`Jackknife (m /. 100.))
            | _ -> 
                    (* TODO *)
                    failwith 
                    "TODO: The only resamples supported are Bootstrap and
                    Jackknife"
        in
        let data, new_otus = Node.load_data new_data in
        let trees = 
            let node_data = 
                List.fold_left (fun acc x ->
                All_sets.IntegerMap.add (Node.taxon_code x) x acc)
                All_sets.IntegerMap.empty
                new_otus
            in
            Sexpr.map (fun x ->
                    TO.uppass (TO.downpass 
                    { x with Ptree.node_data = node_data })) trees
        in
        let perturbed = B.build_initial_trees trees data new_otus build in
        let set = 
            let `LocalOptimum tabu = search in
            TreeSearch.sets tabu.Methods.tabu_join data perturbed 
        in
        let res_trees = 
            PTS.find_local_optimum data queue perturbed set search 
        in
        let tree, _ = 
            match Sexpr.to_list res_trees with
            | tree :: res_trees ->
                    List.fold_left 
                    (fun ((trees, cost) as acc) x ->
                        let cx = Ptree.get_cost `Adjusted x in
                        if cx < cost then ([x], cx)
                        else if cx = cost then ((x :: trees), cx)
                        else acc) 
                    ([tree], Ptree.get_cost `Adjusted tree) 
                    res_trees
            | [] -> assert false
        in
        let ntree = 
            List.map (fun tree -> 
                (TO.collapse_as_needed tree),
                let tree = tree.Ptree.tree in
                Tree.reroot (outgroup, Tree.get_parent outgroup tree) tree)
            tree
        in
        Ptree.add_consensus_to_counters counters ntree

    let bremer_to_sexpr ?edges tree =
        (* function returns None if edge should not be added, Some (edge, otu_set)
           for good edges *)
        let visit (Tree.Edge (h, n) as edge) =
            let otus = Ptree.get_leaves n tree in
            let n_otus = List.length otus in
            if n_otus = 1
            then None
            else
                let otu_set = List.fold_left
                    (fun set elt -> All_sets.Integers.add
                         (Node.taxon_code elt)
                         set) All_sets.Integers.empty otus in
                Some (`Single (edge, otu_set))
        in

        (* Fold over passed-in edges or iterate over all edges *)
        match edges with
        | Some edges -> begin
              let visit acc edge =
                  match visit edge with
                  | None -> acc
                  | Some s -> s :: acc in
              `Set (List.fold_left visit [] edges)
          end
        | None -> begin
              let visit edge list =
                  match visit edge with
                  | None -> (Tree.Continue, list)
                  | Some s -> (Tree.Continue, s :: list) in
              let list_of_handle h acc =
                  Ptree.pre_order_edge_visit visit h tree acc in
              let handles = tree.Ptree.tree.Tree.handles in
              `Set (All_sets.Integers.fold list_of_handle handles [])
          end

    (** [bremer_item ?status tree (edge, otu_set)] makes the tree based on an entry
        in the sexpr as produced by !bremer_to_sexpr *)
    let bremer_item ?status original_tree tree build nodes (edge, otu_set) =
        let () = match status with
        | None -> ()
        | Some status ->
              Status.full_report ~msg:("Making Bremer tree for "
                                       ^ "clade with "
                                       ^ string_of_int
                                       (All_sets.Integers.cardinal otu_set)
                                       ^ " OTUs")
                  status in
        let treecost = Ptree.get_cost `Adjusted original_tree in
        if ndebug_bremer_best_costs
        then odebug ("Tree cost is " ^ string_of_float treecost);
        tree

    let minimum item =
        let comparison acc ((n, tree) as res) =
            match acc with
            | None -> Some res
            | Some (n, othertree) ->
                    if (Ptree.get_cost othertree) < (Ptree.get_cost tree) then
                        acc
                    else Some res
        in
        match Sexpr.fold_left comparison None item with
        | None -> failwith "Empty sexpr?"
        | Some v -> v

    let bremer_run_sexpr trees data queue nodes search build tree sexpr =
        let tree_cost = Ptree.get_cost `Adjusted tree 
        and number_of_nodes = List.length nodes in
        let calculate_one_tree ((Tree.Edge (_, n)), otu_set) =
            if number_of_nodes = All_sets.Integers.cardinal otu_set then
                n, infinity
            else
                let initial_trees = B.build_initial_trees trees data nodes build in
                let nodes = List.map (Node.T.add_exclude otu_set) nodes in
                let res = 
                    match 
                    Sexpr.fold_left (fun acc x ->
                        let x = 
                            let x = CT.substitute_nodes nodes x in
                            let set = 
                                let `LocalOptimum tabu = search in
                                TreeSearch.sets tabu.Methods.tabu_join 
                                data (`Single x) 
                            in
                            PTS.find_local_optimum data queue (`Single x) set search 
                        in
                        Sexpr.fold_left (fun acc x ->
                            let xc = Ptree.get_cost `Adjusted x in
                            if xc < tree_cost then begin
                                Status.user_message Status.Error 
                                ("The following tree was found@ during a Bremer@ " ^
                                "support@ search@ and has better cost:");
                                PTS.report_trees [`Cost] None (`Single x);
                            end else ();
                            match acc with
                            | None -> Some xc
                            | Some c ->
                                    if c < xc then acc
                                    else Some xc) 
                        acc
                        x)
                    None
                    initial_trees
                    with
                    | Some x -> x
                    | None -> failwith "no trees?"
                in
                n, res
        in
        Sexpr.map_status "Calculating Bremer Supports" ~eta:true 
        calculate_one_tree sexpr

    let bremer_sexpr_to_support tree handle tree_cost sexpr =
        let alist = Sexpr.to_list sexpr in

        let rec descend parent id =
            match All_sets.IntegerMap.find id (tree.Ptree.tree.Tree.u_topo)
            with
            | Tree.Leaf (id, _) ->
                    let node = Ptree.get_node_data id tree in
                    let code = Node.taxon_code node in
                    Methods.Leaf code
            | (Tree.Interior (id, a, b, c) as node) ->
                  let left, right = Ptree.other_two_nbrs parent node in
                  let cost =
                      try 
                          let res = (List.assoc id alist) in
                          if res >= (infinity /. 2.) then Pervasives.infinity
                          else res -. tree_cost
                      with Not_found ->
                          Pervasives.infinity in
                  Methods.Node (cost,
                                descend id left, descend id right)
            | Tree.Single _ -> assert false
        in

        let tree = match Ptree.get_node handle tree with
        | Tree.Single i ->
                let node = Ptree.get_node_data i tree in
                let code = Node.taxon_code node in
                Methods.Leaf code
        | Tree.Leaf (id, parent) ->
                let node = Ptree.get_node_data id tree in
                let code = Node.taxon_code node in
                let leaf = Methods.Leaf code in
              Methods.Node (Pervasives.infinity, leaf, descend id parent)
        | Tree.Interior (id, parent, left, right) ->
              Methods.Node (Pervasives.infinity,
                            descend id parent,
                            descend parent id)
        in

        tree

    (** [bremer_singletbr_fps tree handle cost param] uses a fingerprint technique
        to do a single-pass TBR calculation of Bremer supports.  This is
        inefficient due to an inefficient implementation of Tree.CladeFP *)
    let bremer_singletbr_nonadds tree handle cost (param : Methods.local_optimum)
    data emergency_queue =
        let search = param in

        (** [enumerate_list list] returns an assoc. list from list elements to
            indices 0..n-1, and also returns the length of the original list *)
        let enumerate_list ?(starting=0) list =
            let counter = ref starting in
            let ind =
                List.map
                    (fun id ->
                         let c = !counter in
                         incr counter;
                         (id, c))
                    list
            in
            (ind, !counter) in

        let orig_cost = Ptree.get_cost `Adjusted tree in

        (* assoc list: node id -> array index *)
        let node_indices, num_nodes =
            enumerate_list ~starting:(succ (!(data.Data.character_code_gen)))
                (Ptree.get_node_ids tree) in
        let num_leaves = List.length (Ptree.get_all_leaves_ids tree) in

        (* we make a list of leaf data *)
        let starting = succ !(data.Data.character_code_gen) in
        let node_data = 
            let node_data, leaves, nodes =
                let leaves = Ptree.get_all_leaves_ids tree in
                let node_data = 
                    List.map (fun x -> x, Ptree.get_node_data x tree) leaves
                and nodes = Ptree.get_node_ids tree in
                node_data, leaves, nodes
            in
            Node.for_support starting node_data leaves nodes 
        in
        (* make the new tree *)
        let tree' = CT.substitute_nodes node_data tree in
        let orig_cost' = Ptree.get_cost `Adjusted tree' in
        assert (orig_cost = orig_cost');

        (* make the queue and perform the search *)
        let best_vals = Array.make num_nodes infinity in
        (* assume one component. TODO: multiple components. *)
        let status =
            let n = num_leaves in
            (* (2n-3)! / ( 2^(n-2) * (n-2) ) *)
        let tmp1 =  Utl.factorial (2*n-3) in
        let tmp1 = float_of_int tmp1 in
        let nbhd =  tmp1 /. ( (2.**(float_of_int (n-2))) *. (float_of_int (n-2)) )
        in
            Status.create "Bremer TBR neighborhood" (Some (int_of_float nbhd)) "processing" in
        let () = Status.full_report status in
        let nbhood_count = ref 0 in
        let queue = 
            (best_vals, node_indices, starting,
            status, nbhood_count, orig_cost')
        in
        Status.full_report ~msg:"Performing neighborhood search" status;
        let set = 
            let `LocalOptimum tabu = search in
            TreeSearch.sets tabu.Methods.tabu_join data (`Single tree') in
        let _ = 
            PTS.find_local_optimum ~queue data emergency_queue (`Single
            tree') set search in
        Status.full_report ~msg:"Gathering results" status;

        let find_cost parent id =
            best_vals.(List.assoc id node_indices) in

        let rec descend parent id =
            match All_sets.IntegerMap.find id (tree.Ptree.tree.Tree.u_topo)
            with
            | Tree.Leaf (id, _) ->
                let node = Ptree.get_node_data id tree in
                let code = Node.taxon_code node in
                Methods.Leaf code
            | (Tree.Interior (id, a, b, c) as node) ->
                  let left, right = Ptree.other_two_nbrs parent node in
                  let cost =
                      try (find_cost parent id) -. orig_cost
                      with Not_found ->
                          (try 
                              let res = (find_cost id parent) in
                              if res >= infinity then Pervasives.infinity
                              else res -. orig_cost
                           with Not_found -> Pervasives.infinity) in
                  Methods.Node (cost,
                                descend id left, descend id right)
            | Tree.Single _ -> assert false
        in

        let tree = match Ptree.get_node handle tree with
        | Tree.Single id ->
                let node = Ptree.get_node_data id tree in
                let code = Node.taxon_code node in
                Methods.Leaf code
        | Tree.Leaf (id, parent) ->
                let node = Ptree.get_node_data id tree in
                let code = Node.taxon_code node in
                let leaf = Methods.Leaf code in
              Methods.Node (Pervasives.infinity, leaf, descend id parent)
        | Tree.Interior (id, parent, left, right) ->
              Methods.Node (Pervasives.infinity,
                            descend id parent,
                            descend parent id)
        in
        Status.finished status;
        tree
        

    (** [bremer_support tree handle cost search] calculates a Bremer support for the
        given tree.  If the search strategy is [single_tbr], then a special method is
        used that runs a single search to calculate all supports.  If the search
        strategy is anything else, a search will be run for each edge (excluding the
        clade below that edge in the original tree); the results of that search will be
        used to assess the support for that edge. *)
    let bremer_support trees my_rank modul tree handle cost 
    (search : Methods.local_optimum) (build :
        Methods.build) data nodes queue =
        let `LocalOptimum (search_type) = search in
        match search_type.Methods.ss with
        | `SingleNeighborhood `Tbr ->
              bremer_singletbr_nonadds tree handle cost search data queue
        | _ ->
              let sexpr = 
                  if my_rank = 0 && modul = 1 then
                      bremer_to_sexpr tree 
                  else 
                      let all_edges = Ptree.get_edges_tree tree in
                      let all_edges = 
                          List.sort Tree.EdgeComparator.compare all_edges 
                      in
                      let _, edges = 
                          List.fold_left (fun (cnt, acc) edge ->
                              if my_rank = cnt mod modul then 
                                  (cnt + 1), (edge :: acc)
                              else (cnt + 1), acc) (0, []) all_edges
                      in
                      let edges = List.rev edges in
                      bremer_to_sexpr ~edges tree
              in
              let sexpr = bremer_run_sexpr trees data queue nodes search build tree sexpr in
              bremer_sexpr_to_support tree handle cost sexpr

    let support trees nodes meth data queue =
        try
            let get_outgroup = function
                | None -> Node.taxon_code (List.hd nodes)
                | Some v -> v
            in
            let sampler_function, iterations = 
                match meth with
                | `Bootstrap (iterations, search_meth, build_meth, out) -> 
                      let outgroup = get_outgroup out in
                      resample_support trees data queue iterations nodes `Bootstrap
                      search_meth build_meth outgroup,
                      iterations
                | `Jackknife (jack_par, iterations, search_meth, build_meth, out) ->
                      let outgroup = get_outgroup out in
                      resample_support trees data queue iterations nodes 
                      (`Jackknife jack_par) search_meth build_meth outgroup, 
                      iterations
            in
            let st = 
                Status.create "Calculating Support" (Some iterations)
                "Iterations"
            in
            let rec iterator cnt max f acc =
                Status.full_report ~adv:cnt st;
                if cnt > max then acc
                else iterator (cnt + 1) max f (f acc)
            in
            let res = 
                iterator 1 iterations sampler_function Tree.CladeFPMap.empty
            in
            Status.finished st;
            res
        with
        | Failure "hd" -> failwith "Illegal support: of an emty node set?"

    let bremer_support build_trees my_rank modul nodes trees search_method build_method data queue =
        let f tree =
            let cost = Ptree.get_cost `Adjusted tree in
            bremer_support
            build_trees
            my_rank
            modul
            tree
            (All_sets.Integers.choose (tree.Ptree.tree.Tree.handles))
            cost
            search_method build_method
            data nodes queue 
        in
        Sexpr.map_status "Support for trees" f trees

    let support_to_xml out diag =
        let print = output_string out in
        let rec d = function
            | Methods.Leaf i -> print ("<leaf id=\"" ^ string_of_int i ^ "\" />\n")
            | Methods.Node (f, left, right) ->
                  print ("<node support=\"" ^ string_of_float f ^ "\">\n");
                  d left;
                  d right;
                  print "</node>\n"
        in d diag

    let support_to_parser_tree to_a diag =
        let rec d = function
            | Methods.Leaf i -> Tree.Parse.Leafp (to_a (Methods.Leaf i))
            | Methods.Node (f, left, right) ->
                  Tree.Parse.Nodep ([d left; d right;],
                                    to_a (Methods.Node (f, left, right)))
        in Tree.Parse.Flat (d diag)

    let support_to_string_tree data = support_to_parser_tree
        (function
             | Methods.Leaf i -> Data.code_taxon i data
             | Methods.Node (f, _, _) -> string_of_float f)

    (** [join_support_trees trees] takes a list of [(iterations, support_tree)]
        pairs and combines them into a single support tree *)
    let join_support_trees trees =
        let join2 (i1i, t1) (i2i, t2) =
            let i1 = float_of_int i1i in
            let i2 = float_of_int i2i in
            let rec join2 t1 t2 =
                match t1, t2 with
                | Methods.Leaf i, Methods.Leaf j when i = j -> Methods.Leaf i
                | Methods.Node (w1, t1l, t1r), Methods.Node (w2, t2l, t2r) ->
                      Methods.Node (((w1 *. i1) +. (w2 *. i2)) /. (i1 +. i2),
                            join2 t1l t2l,
                            join2 t1r t2r)
                | _ -> failwith "join_support_trees/join2: incompatible trees"
            in
            (i1i + i2i, join2 t1 t2) in
        match trees with
        | [] -> failwith "join_support_trees: no trees!"
        | tree :: trees ->
              let (iter, tree) = List.fold_left join2 tree trees in
              tree

    let generate_sets process_cost data root tree =
        let pair tree = 
            match Tree.get_node root tree with
            | Tree.Leaf (a, b) 
            | Tree.Interior (a, b, _, _) -> (a, b)
            | Tree.Single _ -> failwith "Not a tree?"
        in
        let expected_cost = match tree with
            | Tree.Parse.Annotated (_,str) -> str
            | _ -> failwith "No Expected Cost" in
        let tree = 
            Tree.Parse.convert_to (None,[tree]) 
            (fun taxon -> Data.taxon_code taxon data)
            (Data.number_of_taxa data)
        in
        let tree = Tree.reroot (pair tree) tree in
        process_cost tree expected_cost, Tree.CladeFP.sets tree

    let generate_sets_for_bremer process_cost root data file =
        let trees = Tree.Parse.of_file file in
        let for_brem = List.flatten (List.map (List.map (generate_sets
        process_cost data root)) trees) in
        Sexpr.of_list for_brem

    let bremer_of_input_file process_cost root to_string data file 
    (trees : (a, b) Ptree.p_tree Sexpr.t) =
        let process_one_tree tree : Tree.Parse.tree_types =
            Ptree.bremer to_string
            (int_of_float (Ptree.get_cost `Adjusted tree))
            tree.Ptree.tree (generate_sets process_cost data root) file
        in
        Sexpr.map_status "Bremer of input tree" ~eta:true process_one_tree trees


    let bremer_of_input_file_but_trust_input_cost =
        bremer_of_input_file (fun _ b -> try int_of_float (float_of_string b)
        with _ -> max_int)

end

module Make (NodeH : NodeSig.S with type other_n = Node.Standard.n) (EdgeH : Edge.EdgeSig with type n = NodeH.n) 
    (TreeOpsH : 
        Ptree.Tree_Operations with type a = NodeH.n with type b = EdgeH.e) = struct

    type a = NodeH.n
    type b = EdgeH.e

    module SH = MakeNormal (Node.Standard) (Edge.SelfEdge) (Chartree.TreeOps)
    module DH = MakeNormal (NodeH) (EdgeH) (TreeOpsH)

    module TOS = Chartree.TreeOps 
    module TOH = TreeOpsH 

    let replace_contents downpass uppass get_code nodes ptree =
        let nt = { (Ptree.empty ptree.Ptree.data) with 
            Ptree.tree = ptree.Ptree.tree } in
        uppass (downpass 
        (List.fold_left (fun nt node ->
            Ptree.add_node_data (get_code node) node nt) 
        nt nodes))

    let from_s_to_h = replace_contents TOH.downpass TOH.uppass NodeH.taxon_code 
    let from_h_to_s = replace_contents TOS.downpass TOS.uppass Node.Standard.taxon_code

    let support trees a b data c =
        if Data.has_dynamic data then
            DH.support trees a b data c
        else
            let data, nodes = Node.Standard.load_data data in
            let trees = Sexpr.map (from_h_to_s nodes) trees in
            SH.support trees nodes b data c

    let bremer_support build_trees a b c itrees d e data f =
        if Data.has_dynamic data then
            DH.bremer_support build_trees a b c itrees d e data f
        else
            let data, nodes = Node.Standard.load_data data in
            let build_trees = Sexpr.map (from_h_to_s nodes) build_trees 
            and itrees = Sexpr.map (from_h_to_s nodes) itrees in
            SH.bremer_support build_trees a b nodes itrees d e data f

    let support_to_string_tree = DH.support_to_string_tree

    let join_support_trees = DH.join_support_trees

    let bremer_of_input_file = DH.bremer_of_input_file

    let bremer_of_input_file_but_trust_input_cost = 
        DH.bremer_of_input_file_but_trust_input_cost

end
