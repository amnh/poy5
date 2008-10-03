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

(** queues.ml - contains the serach managers that manage the queue of trees to
 * be searched. *)

(* $Id: queues.ml 2272 2007-10-05 15:03:07Z andres $ *)
let () = SadmanOutput.register "Queues" "$Revision: 2272 $"

(** {1 Types} *)

(** The incemental optimization type. When set of vertices are modified and
* might be subject to a second pass of modifications (in general an uppass),
* that is missing to complete their update, each constructor holds information
* for weather or not a vertex should be modified, and weather or not the list
* holding the items has already it's children or no children at all.  *)
type incremental = [ 
    | `Children of int
    | `HandleC of int * int
    | `HandleNC of int * int
    | `No_Children of int ]

(** The type of the module that is built by the functor {!Queues.Make}. *)
module type S = sig

    (** {2 Types} *)
    type a (* The type of the data in the vertices of the tree *)
    type b (* The type of the data in the edges of the tree *)

    (** Get the number of trees that have been evaluated after the last
    * [reset_trees_considered ()] *)
    val get_trees_considered : unit -> int

    (** Set the counter of trees evaluated to 0 *)
    val reset_trees_considered : unit -> unit

    (** {2 Wagner Build Managers} *)

    class type wagner_mgr = [a, b] Ptree.wagner_mgr

    (** Constructor for the wagner search manager used in all the building
    * procedures. The constructor takes two arguments, the maximum number of
    * trees that are accepted in the wagner building step, and the maximum
    * threshold above the current optimum that is acceptable for a tree to be
    * elgible for the next build round (next vertex addition). *)
    class wagner_srch_mgr : bool -> int -> float -> wagner_mgr

    (** {2 Local Search Managers *)
    (** All the local search managers take a Sampler.search_manager as an
    * argument. That sampler allows the collection of data about the searches
    * for research puproses. *)

    class type search_mgr = [a, b] Ptree.search_mgr

    (** A search manager that chooses the first tree that shows a better cost
    * than the current best tree *)
    class first_best_srch_mgr : 
        (a, b) Sampler.search_manager_sampler -> 
            search_mgr

    (** A search manager that keeps up to n trees that have the same cost as the
    * current best tree. The preference order is given by the keep_method
    * argument. *)
    class hold_n_fb_srch_mgr :
        int ->
            Methods.keep_method ->
                (a, b) Sampler.search_manager_sampler -> 
                    search_mgr

    (** A search manager that keeps up to n trees that have cost within the
    * threshold of the current tree. The preference order is given by the
    * keep_method argumen *)
    class hold_n_threshold_srch_mgr :
        int ->
            Methods.keep_method ->
                float ->
                    (a, b) Sampler.search_manager_sampler ->
                        search_mgr

    class annealing_srch_mgr :
        int ->
            Methods.keep_method ->
                float ->
                    float ->
                        float ->
                            (a, b) Sampler.search_manager_sampler ->
                                search_mgr

    class classic_poy_drifting_srch_mgr :
        int -> Methods.keep_method -> float -> float ->
                            (a, b) Sampler.search_manager_sampler ->
                                search_mgr

    (** A search manager that instead of choosing the first tree that shows a
    * better cost, evaluates all the trees that are in the neighborhood of the
    * current best to continue the search. *)
    class all_possible_joins :
        string option ->
            (a, b) Sampler.search_manager_sampler ->
                search_mgr

    (* A support manager class that handles the searches for the support values
    * calucations. The constructor takes arguments [a b c d e f g], where [a]
    * is the array of indices of clades and their corresponding best cost found
    * so far, [b] is the association list that connects character codes with their
    * corresponding clade (we mark clades using a specially added character),
    * [c] is the initial code of the special characters, [d] is the status
    * holder for the user interface updates, [e] is the counter of the number of
    * rearrangementes of the support calculation (for the user interface
    * update), [f] is the origin_cost for forest searches, and [g] is the
    * sampler search manager. For more information check it's usage in the
    * {!Support} module. *)
    class supports_manager : 
        float array -> (int * int) list -> int ->
            Status.status -> int ref -> float ->
                (a, b) Sampler.search_manager_sampler ->
                    search_mgr

    class all_neighbors_srch_mgr : 
        (a, b) Ptree.join_fn ->
        (a, b) Ptree.break_fn ->
        (a, b) Ptree.reroot_fn ->
        ((a, b) Ptree.p_tree -> Ptree.incremental list
        -> (a, b) Ptree.p_tree) ->
        (a, b) Sampler.search_manager_sampler -> search_mgr
end

(** {1 Functor} *)
(** The functor to construct a module of queue managers for a tree with nodes
* and edges*)
module Make (Node : NodeSig.S) (Edge : Edge.EdgeSig with type n = Node.n) 
    : S with type a = Node.n with type b = Edge.e = struct

    type a = Node.n
    type b = Edge.e
    class type wagner_mgr = [a, b] Ptree.wagner_mgr
    class type search_mgr = [a, b] Ptree.search_mgr
    let data_mining = false
    let debug_costfn = false
    let debug_costfn_callback = None
    let debug_costfn_opportunistic = true
    let debug_fb_infcost_result = true
    let odebug = Status.user_message Status.Information

    let total_trees_considered = ref 0
    let get_trees_considered () = !total_trees_considered
    let reset_trees_considered () = total_trees_considered := 0
    let incr_trees_considered () = incr total_trees_considered
    let add_trees_considered i = 
        total_trees_considered := !total_trees_considered + i


    (* This is a debugging function for when a cost is found to be incorrect.  It
       prints the nodes involved in the attachment. *)
    let debug_costfn_fn expected actual j1 j2 cladenode prevtree newtree =
        print_endline "Joining at:";
        match j1 with
        | Tree.Single_Jxn h ->
              let h = Tree.int_of_id h in
              print_endline (Node.to_string
                                 (Ptree.get_node_data h prevtree))
        | Tree.Edge_Jxn (h, n) ->
              let h = Tree.int_of_id h in
              let n = Tree.int_of_id n in
              (print_endline
                  (Node.to_string
                       (Ptree.get_node_data h prevtree) ^ " ->");
               print_endline
                   (Node.to_string
                        (Ptree.get_node_data n prevtree)));
        let j2data = Ptree.get_node_data cladenode prevtree in
        print_endline "Join node:";
        print_endline (Node.to_string j2data)
    (* let debug_costfn_callback = Some debug_costfn_fn *)
        

    class counting_srch_mgr =
    object (self)
        val mutable tree_count = 0

        method count =
            tree_count <- tree_count + 1;
            incr_trees_considered ()

        method tried_trees = tree_count
    end

    (** a search manager of the wagner tree construction. This keeps track of the
    * best tree at each round of node addition. *)
    class wagner_srch_mgr verify_cost max thresh : wagner_mgr = 
        object (self)
        
        val mutable c_delta = Ptree.NoCost 
            (** the current delta *)
            
        val mutable m_max = max

        val true_max = max

        val m_thr = thresh 

        val mutable srch_trees = [] 
            (** the current best tree and cost *)

            (** Function to initialize the object explicitly. *)
        method init tree_cost_delta_lst =
            match tree_cost_delta_lst with
            | (_, _, d, _) :: _ -> 
                    let rec get_upto cnt cur acc = 
                        match cnt, cur with
                        | 0, _ -> 0, List.rev acc
                        | _, [] -> cnt, List.rev acc
                        | _, ((_, cst, c, tabu_mgr) as h) :: t -> 
                                let acc = (((Lazy.lazy_from_val h), cst, c) :: acc) in
                                get_upto (cnt - 1) t acc
                    in
                    let cnt, tst = get_upto max tree_cost_delta_lst [] in
                    srch_trees <- tst;
                    c_delta <- d;
                    m_max <- true_max - cnt;
            | [] ->
                    raise (Invalid_argument "The initialization must contain at \
                    least 1 element.")

            (** Function to clone the current instance and fill it with default
             * values. *)
        method clone = 
            new wagner_srch_mgr verify_cost true_max m_thr
                
            (** Function to return whether there are any more trees to search. *)
        method any_trees =
            match srch_trees with
            | _ :: _ -> true 
            | [] -> false 

            (** Function to get the next tree. *)
        method next_tree = 
            match srch_trees with
            | (x, _, c) :: tl ->
                    srch_trees <- tl;
                    c_delta <- c;
                    m_max <- m_max + 1;
                    let a, b, _, c = Lazy.force x in
                    a, b, c
            | _ -> 
                    failwith "There are no more trees to search in the wagner \
                    manager."
                
        method private filterout m_max trees real_cost new_cost old_best =
            let filter_maximum (worst_cost, tree, res) ((_, cur_cost, _) as v) =
                if cur_cost > worst_cost then begin
                    match tree with
                    | None -> (cur_cost, Some v, res)
                    | Some x -> (cur_cost, Some v, (x :: res))
                end else (worst_cost, tree, (v :: res))
            in
            (*
            let do_filter (counter, results) ((_, cur_cost, _) as v) =
                if cur_cost > (min new_cost old_best) +. m_thr then (counter + 1, results)
                else (counter, v :: results)
            in
            match List.fold_left do_filter (0, []) trees with
            | (0, _) as x -> 
                    *)
                    if m_max = 0 then 
                        let acc = (real_cost, None, []) in
                        let worst_cost, tree, res = 
                            List.fold_left filter_maximum acc trees
                        in
                        match tree with
                        | None -> (false, 0, res)
                        | Some _ -> 
                                (true, 1, res)
                    else if m_max < 0 then assert false
                    else true, 0, trees
                    (*
            | x -> x
                    *)
        method evaluate =
            let checker ((m_max, cur_best, results) as r) (v, estimate, _) =
                if estimate > cur_best +. m_thr then  r
                else begin
                    let a, b, c, tabu_mgr, cost_for_comparisons = 
                        let a, b, c, tabu_mgr = Lazy.force v in
                        if verify_cost then a, b, c, tabu_mgr, b
                        else a, b, c, tabu_mgr, estimate
                    in
                    if ( cost_for_comparisons <= cur_best ) then begin
                        (* We are ok, we must add this one *)
                        let is_there_room, counter, results = 
                            self#filterout m_max results cost_for_comparisons 
                            cost_for_comparisons cur_best 
                        in
                        if is_there_room then
                            m_max - 1 + counter, cost_for_comparisons,
                            ((Lazy.lazy_from_val (a, b, c, tabu_mgr)), b, c) :: 
                                results
                        else m_max, cur_best, results
                    end else if cost_for_comparisons < cur_best +. m_thr then begin
                        (* We are more or less ok, we make the decision based on
                        * m_max *)
                        if m_max > 0 then begin
                            m_max - 1, cur_best,
                            ((Lazy.lazy_from_val (a, b, c, tabu_mgr)), b, c) ::
                                results
                        end else begin
                            let is_there_room, counter, results = 
                                self#filterout m_max results
                                cost_for_comparisons cost_for_comparisons 
                                cur_best 
                            in              
                            if is_there_room then
                                m_max - 1 + counter, cur_best,
                                ((Lazy.lazy_from_val (a, b, c, tabu_mgr)), b, c) ::
                                    results
                            else m_max, cur_best, results
                        end
                    end else 
                        m_max, cur_best, results
                end
            in
            let new_max, _, res = 
                List.fold_left checker (true_max, infinity, []) srch_trees
            in
            assert (true_max - new_max = List.length res);
            m_max <- new_max;
            srch_trees <- res

             (** [process cost_fn b_delta cd_nd join_fn j1 j2 t_delta]
             @param cost_fn function to return the cost of the attempted join. 
             @param b_delta the cost of pruning the clade. 
             @param cd_nd the node that was on the clade side of the break_jxn. 
             @param join_fn function to perform the join operation. 
             @param j1, j2 the join junctions. 
             @param t_delta the change in the tree topology due to clade prune
                        operation. 
             @return Function to decide whether to perform a join
             and possibly store the resulting tree for further searching.  *)
        method process (cost_fn : (a, b) Ptree.cost_fn)
            (b_delta : float)
            (cd_nd : a)
            (join_fn : (a, b) Ptree.join_fn)
            (j1 : Tree.join_jxn)
            (j2 : Tree.join_jxn) 
            pt 
            (tabu_mgr : (a, b) Ptree.wagner_edges_mgr) = 
                (* Ignores the b_delta as it does not make sense for 
                 * wagner trees algo *)
            let cld_cst = 
                match c_delta with
                | Ptree.NoCost -> infinity
                | Ptree.Cost(x) -> x
            in
            match cost_fn j1 j2 cld_cst cd_nd pt with
            | Ptree.Cost cc ->
                    let (_, real_cost, _) as v = 
                        let tabu_mgr = tabu_mgr#clone in
                        (Lazy.lazy_from_fun (fun () ->
                        let nt, dlt = join_fn [] j1 j2 pt in
                        let cst = Ptree.get_cost `Adjusted nt in
                        tabu_mgr#update_join nt dlt;
                        (nt, cst, c_delta, tabu_mgr)), 
                        (Ptree.get_cost `Adjusted pt +. cc), c_delta)
                    in
                    tabu_mgr#break_distance cc;
                    if ( cc <= cld_cst +. m_thr ) then begin
                        (* We are ok, we must add this one *)
                        let is_there_room, counter, new_search_trees = 
                            self#filterout m_max srch_trees real_cost cc cld_cst 
                        in
                        if is_there_room then begin
                            tabu_mgr#new_delta cc;
                            c_delta <- Ptree.Cost cc;
                            srch_trees <- v :: new_search_trees;
                            m_max <- (m_max - 1 + counter) + counter;
                        end;
                        Tree.Continue 
                    end else Tree.Continue
            | Ptree.NoCost -> Tree.Skip

            (** Function to return the results. *)
        method results = 
            let mapper = fun (x, _, _) ->
                let (a, b, _, _) = Lazy.force x in
                a,b
            in
            List.map mapper srch_trees
            
    end 

    (** A simple search manager that picks the first tree that is better than the
     * current one and continues search with that tree. *)
    class first_best_srch_mgr (sampler : (Node.n, Edge.e) Sampler.search_manager_sampler) = object 

        val mutable c_delta = Ptree.NoCost 
        (** the current delta *)
            
        val mutable srch_trees = [] 
        (** the current best tree and cost *)

        val mutable results = [] 
        (** the results are stored here. *)

        val mutable cur_best_cost = Pervasives.max_float

        method features lst =
            ("type", "first best") :: lst

        val mutable do_repeat = false

        method should_repeat = do_repeat

        (** [init list] initializes the contents of the queue.  [list] must be a
            list of triplets [(tree, cost, delta)]. *)
        method init (tree_cost_delta_lst : 
                         (('a, 'b) Ptree.p_tree * float * Ptree.clade_cost *
                              ('a, 'b) Ptree.tabu_mgr) list) = 
            sampler#init (List.map (fun (a, b, c, _) -> (a, b, c))
            tree_cost_delta_lst);
            match tree_cost_delta_lst with
            | [(tr, cst, delta, tabu)] -> 
                  begin
                      srch_trees <- [(tr, cst, tabu)] ;
                      results <- [(tr, cst, tabu)] ;
                      c_delta <- delta ;
                      cur_best_cost <- cst
                  end
            | _ -> raise (Invalid_argument "List must be of length 1")

        (** Function to clone the current instance and fill it with default
         * values. *)
        method clone = 
            let nsampler = sampler#clone in
            new first_best_srch_mgr nsampler
                
        (** Function to return whether there are any more trees to search. *)
        method any_trees =
            let res = 
                match srch_trees with
                | [] -> false 
                | [(_, _, _)] -> true 
                | _ -> raise (Failure "More than one tree stored in list") 
            in
            sampler#any_trees res;
            res

        (** Function to get the next tree. *)
        method next_tree = 
            let ret = List.hd srch_trees in
            srch_trees <- [] ;
            sampler#next_tree (let (a, b, _) = ret in (a, b));
            ret
                
        (** Function to decide whether to decide whether to perform a join and
         possibly store the resulting tree for further searching. This function
         has to update the tabu manager and keep the list of edges in the tree
         in sync with list of edges in the tabu. *)
        method process
            (cost_fn : (Node.n, Edge.e) Ptree.cost_fn)
            (b_delta : float)
            (cd_nd : 'a)
            (join_fn : (Node.n, Edge.e) Ptree.join_fn)
            incremental
            (j1 : Tree.join_jxn)
            (j2 : Tree.join_jxn) 
            (tabu_mgr : (Node.n, Edge.e) Ptree.tabu_mgr)
            pt = 
            incr_trees_considered ();
            let clade_cost = (cost_fn j1 j2 b_delta cd_nd pt) in
            match clade_cost with
            | Ptree.NoCost -> Tree.Skip
            | Ptree.Cost(cc) ->
                    tabu_mgr#break_distance cc;
                  if debug_costfn then begin
                      if cur_best_cost <> max_float then begin
                          let expected_new_cost =
                              cur_best_cost -. b_delta +. cc in
                          let newtree, tree_delta =
                              join_fn incremental j1 j2 pt in
                          let actual_new_cost = 
                              Ptree.get_cost `Adjusted newtree 
                          in
                          if expected_new_cost <> actual_new_cost
                          then begin
                              Status.user_message
                                  Status.Error
                                  ("Cost function not exact. Claimed cost "
                                   ^ string_of_float expected_new_cost
                                   ^ " (delta = " ^ string_of_float cc ^ ")"
                                   ^ ", actual cost "
                                   ^ string_of_float actual_new_cost
                                   ^ ", off by "
                                   ^ string_of_float (expected_new_cost -.
                                                          actual_new_cost)
                                  );
                              match debug_costfn_callback with
                              | None -> ()
                              | Some f ->
                                    f
                                        expected_new_cost
                                        actual_new_cost
                                        j1 j2 cd_nd pt
                                        newtree
                          end
    (*                       else print_endline "Good: Cost function is exact!" *)
                      end
                  end;
                  if ( cc < b_delta ) then begin
                      c_delta <- (Ptree.Cost(cc)) ;
                      let nt, j_delta = (join_fn incremental j1 j2 pt) in
                      let cst = Ptree.get_cost `Adjusted nt in
                      if cst < cur_best_cost then begin
                          let new_tabu = tabu_mgr#clone in
                          sampler#process incremental j1 j2 cd_nd pt (Some nt) 
                          b_delta cc (Some cst); 
                          cur_best_cost <- cst;
                          new_tabu#update_join nt j_delta;
                          srch_trees <- [(nt, cst, new_tabu)] ;
                          results <- [(nt, cst, new_tabu)] ;
                          Tree.Break 
                      end else begin
                          sampler#process incremental j1 j2 cd_nd pt None b_delta cc None; 
                          Tree.Continue
                      end
                  end else begin
                      sampler#process incremental j1 j2 cd_nd pt None b_delta cc None; 
                      Tree.Continue
                  end

        (** Function to return the results. *)
        method results =
            sampler#results (List.map (fun (a, b, _) -> (a, b)) results);
            results

        val mutable current_break_position = None

        method breakin (edge : Tree.edge) = current_break_position <- Some edge
            
    end




    type fp = Ptree.Fingerprint.t Lazy.t
    let get_fp tree = lazy (Ptree.Fingerprint.fingerprint tree)
    let compare_fp a b =
    (*     let fp_to_string = Ptree.Fingerprint.to_string in *)
        let a = Lazy.force a in
        let b = Lazy.force b in
        let cmp = Ptree.Fingerprint.compare a b in
        cmp

    module OrderTreeCost = struct
        type t = float * fp
                * (Node.n, Edge.e) Ptree.p_tree
                * (Node.n, Edge.e) Ptree.tabu_mgr
        let compare (a, fpa, _, _) (b, fpb, _, _) =
            let cmp =
                Pervasives.compare (a : float) (b : float) in
            if cmp <> 0
            then cmp
            else compare_fp fpa fpb
    end
    module TreeCostSet = Set.Make(OrderTreeCost)

    class ['a] dont_repeat =
    object (self)
        val mutable fingerprints = All_sets.IntegerMap.empty
        method private get_fps_of_cost cost =
            try All_sets.IntegerMap.find (int_of_float cost) fingerprints
            with Not_found -> []
        method private set_fps_of_cost cost list =
            fingerprints <-
                All_sets.IntegerMap.add (int_of_float cost) list fingerprints

        method remember_tree (tree, cost, fp) =
            let _ = (tree : 'a) in
            let list = self#get_fps_of_cost cost in
            if List.exists (fun f -> 0 = compare_fp fp f) list
            then ()
            else self#set_fps_of_cost cost (fp :: list)

        method forget_tree (tree, cost, fp) =
            let _ = (tree : 'a) in
            let list = self#get_fps_of_cost cost in
            let list = List.filter (fun f -> 0 <> compare_fp fp f) list in
            self#set_fps_of_cost cost list

        method new_tree (tree, cost, fp) =
            let _ = (tree : 'a) in
            let list =
                try All_sets.IntegerMap.find (int_of_float cost) fingerprints
                with Not_found -> [] in
            if List.exists (fun f -> 0 = compare_fp fp f) list
            then false
            else true
    end

    class ['a] repeat =
    object
        method remember_tree (_ : 'a) = ()
        method new_tree (_ : 'a) = true
        method forget_tree (_ : 'a) = ()
    end


    (** A simple search manager that picks the first tree that is better than the
     * current one and continues search with that tree. *)
    class hold_n_fb_srch_mgr n (strategy : Methods.keep_method)
    (sampler : (Node.n, Edge.e) Sampler.search_manager_sampler)  =
    object (self)

        val sampler = sampler

        val mutable srch_trees = []
        (** the current best tree and cost *)

        val mutable results = Array.of_list []
        (** the results are stored here. *)

        val mutable ins_index = 0
        (** where in the array to insert new items; meaning is strategy-specific:
            - for `First: have we read up to n?  if not, this is where to insert the
            next tree.  if yes, this will just be n and we won't insert any more.
            - for `Last: this is the array index mod n where we should insert the
            next element.  if we've already inserted at least n elements, we also
            set all_stored to true.
            - for `Keep_Random: total number of trees we've seen.  we use this to
            calculate the correct position for the new element.
        *)

        (** [all_stored]: whether we've filled up our buffer already (nec. for the
            Last strategy) *)
        val mutable all_stored = false

        (** best cost found so far *)
        val mutable cur_best_cost = Pervasives.max_float

        val mutable do_repeat = false

        method should_repeat = do_repeat

        (** [filter t] is a private virtual method that allows you to specify a
            filter for the trees considered for storing *)
        method private filter t cost = cost <= cur_best_cost

        (** [insert tree] is a private method that inserts a new tree in the proper
            place. *)    
        method private insert tree = match strategy with
        | `First ->
              if ins_index < n
              then begin
                      results.(ins_index) <- tree;
                      ins_index <- ins_index + 1
                  end
        | `Last ->
              results.(ins_index) <- tree;
              ins_index <- (ins_index + 1) mod n;
              if ins_index = 0 then all_stored <- true
        | `Keep_Random ->
              if ins_index < n
              then begin
                      results.(ins_index) <- tree;
                      ins_index <- ins_index + 1
                  end
              else begin
                      ins_index <- ins_index + 1;
                      let outside = ins_index - n in
                      let rand = Random.int ins_index in
                      if rand < outside
                      then ()
                      else begin
                              let rand = Random.int n in
                              results.(rand) <- tree
                          end
                  end

        method private clear =
            ins_index <- 0;
            all_stored <- false

        (** Use the current filter to refilter the stored trees *)
        method private refilter =
            let stored = self#results in
            let stored = List.map
                (fun (tr, c, ta) -> (lazy tr, c, lazy ta)) stored in
            let filtered = List.map (fun (t, cst, _) -> self#filter t cst) stored in
            if List.exists not filtered
            then begin
                    self#clear;
                    List.iter2
                        (fun i b -> if b then self#insert i)
                        stored filtered
                end

        method features lst =
            ("type", "hold n first best")
            :: ("holding", string_of_int n)
            :: ("strategy", match strategy with
                | `First -> "keep first"
                | `Last -> "keep last"
                | `Keep_Random -> "keep random set")
            :: lst

        (** [init list] initializes the contents of the queue.  [list] must be a
            list of triplets [(tree, cost, delta)]. *)
        method init (tree_cost_delta_lst :
                         (('a, 'b) Ptree.p_tree * float * Ptree.clade_cost *
                              ('a, 'b) Ptree.tabu_mgr) list) = 
            sampler#init (List.map (fun (a, b, c, _) -> (a, b, c))
            tree_cost_delta_lst);
            match tree_cost_delta_lst with
            | (tr, cst, delta, tabu) :: lst ->
                  begin
                      srch_trees <- List.map
                          (fun (a, b, _, c) -> (a, b, c)) tree_cost_delta_lst;
                      let (tr, cst, delta, tabu) = List.fold_left
                          (fun ((_, cst1, _, _) as t1) ((_, cst2, _, _) as t2) ->
                               if cst1 <= cst2
                               then t1 else t2)
                          (tr, cst, delta, tabu) lst in
                      (* find the best *)
                      cur_best_cost <- cst;
                      let tr = lazy tr in
                      let tabu = lazy tabu in
                      results <- Array.make n (tr, cst, tabu);
                      self#insert (tr, cst, tabu)
                  end
            | _ -> raise (Invalid_argument "List must be of length 1")

        (** Function to clone the current instance and fill it with default
         * values. *)
        method clone = ({< sampler = sampler#clone >} :> (Node.n, Edge.e) Ptree.search_mgr)
                
        (** Function to return whether there are any more trees to search. *)
        method any_trees =
            let res = 
                match srch_trees with
                | [] -> false 
                | _ -> true 
            in
            sampler#any_trees res;
            res

        (** Function to get the next tree. *)
        method next_tree = 
            let (tree, b, tab) as ret = List.hd srch_trees in
            self#insert (lazy tree, b, lazy tab);
            srch_trees <- List.tl srch_trees;
            sampler#next_tree (tree,b);
            ret
                
        (** Function to decide whether to decide whether to perform a join and
         * possibly store the resulting tree for further searching. *)
        method process
            (cost_fn : ('a, 'b) Ptree.cost_fn)
            (b_delta : float)
            (cd_nd : 'a)
            (join_fn : ('a, 'b) Ptree.join_fn)
            incremental
            (j1 : Tree.join_jxn)
            (j2 : Tree.join_jxn) 
            (tabu_mgr : ('a, 'b) Ptree.tabu_mgr)
            pt = 
            incr_trees_considered ();
            let clade_cost = (cost_fn j1 j2 b_delta cd_nd pt) in
            match clade_cost with
            | Ptree.NoCost -> Tree.Skip
            | Ptree.Cost(cc) ->
                    tabu_mgr#break_distance cc;
                    let new_tabu = tabu_mgr#clone in
                    let cost = Ptree.get_cost `Adjusted pt +. cc in
                    let ljoin = lazy (join_fn incremental j1 j2 pt) in
                    let ltree = lazy (let (t, _) = Lazy.force ljoin in t) in
                    let ltabu =
                        lazy (let (t, j) = Lazy.force ljoin in
                        new_tabu#update_join t j;
                        new_tabu) in
                    if self#filter ltree cost then 
                        let cost = 
                            let nt = Lazy.force ltree in
                            Ptree.get_cost `Adjusted nt 
                        in
                        if self#filter ltree cost then
                            self#insert (ltree, cost, ltabu)
                        else ()
                    else ();
                    (* insert self-call here *)
                    self#trajectory incremental cc cd_nd b_delta j1 j2 pt ljoin ltabu cost

        (** [trajectory cc b_delta ljoin ltabu cost] is used to control the
            trajectory of the search.  This version implements first-best search;
            see later for an example of annealing search implemented this way. *)
        method private trajectory incremental cc cd_nd b_delta j1 j2 pt ljoin ltabu cost =
            if ( cc < b_delta ) then begin
                let nt, j_delta = Lazy.force ljoin in
                let cst = Ptree.get_cost `Adjusted nt in
                sampler#process incremental j1 j2 cd_nd pt (Some nt) b_delta 
                cc (Some cst);
                (*
                if debug_costfn_opportunistic && (cst <> cost)
                then (odebug "Wrong cost!"; assert false);
                *)
                if cst < cur_best_cost then begin
                        let new_tabu = Lazy.force ltabu in
                        cur_best_cost <- cst;
                        srch_trees <- [(nt, cst, new_tabu)];
                        self#refilter;
                        Tree.Break 
                    end else begin
                        sampler#process incremental j1 j2 cd_nd pt (Some nt) 
                        b_delta cc None;
                        Tree.Continue
                    end
                end else begin
                    sampler#process incremental j1 j2 cd_nd pt None b_delta cc None;
                    Tree.Continue
                end

        (** Function to return the results.  We don't currently do this in the
            correct order. *)
        method results =
            let res = match strategy with
            | `Last ->
                  if all_stored
                  then Array.to_list results
                  else Array.to_list (Array.sub results 0 (ins_index + 1))
            | _ ->
                  if ins_index >= n
                  then Array.to_list results
                  else Array.to_list (Array.sub results 0 (ins_index + 1))
            in
            let res = 
                List.map
                (fun (t, c, tabu) ->
                     let t = Lazy.force t in
                     (t, Ptree.get_cost `Adjusted t, Lazy.force tabu))
                res
            in
            sampler#results (List.map (fun (a, b, _) -> (a, b)) res);
            assert (0 < List.length res);
            res

        val mutable current_break_position = None

        method breakin (e : Tree.edge) = current_break_position <- Some e
    end

    class all_neighbors_srch_mgr  
        (join_fn : (a, b) Ptree.join_fn) (break_fn : (a, b) Ptree.break_fn)
        (reroot_fn : (a, b) Ptree.reroot_fn) 
        (incremental_uppasss : ((a, b) Ptree.p_tree -> Ptree.incremental list
        -> (a, b) Ptree.p_tree))
        sampler = object (self)
        inherit first_best_srch_mgr sampler as super

        val mutable all_places_to_join = []

        val mutable current_best_gain = None

        method init lst = 
            all_places_to_join <- [];
            current_best_gain <- None;
            super#init lst

        method breakin edge = 
            let _ =
                match current_best_gain, current_break_position with
                | None, _ -> ()
                | Some (gain, r1, r2), Some cbp -> 
                        all_places_to_join <-
                            (gain, cbp, r1, r2) :: all_places_to_join;
                        current_best_gain <- None
                | _ -> failwith "impossible!"
            in
            super#breakin edge;

        method clone = 
            let nsampler = sampler#clone in
            new all_neighbors_srch_mgr join_fn break_fn reroot_fn
            incremental_uppasss nsampler

        val mutable current_tree_in_search = None 

        method private check_results = 
            let attempt_to_make_them_true tree tabu (acc_tree, (acc_tabu : ('a, 'b)
            Ptree.tabu_mgr), acc_todo) (_, break, r1, r2) =
                let icost = Ptree.get_cost `Adjusted acc_tree in
                let ntabu = tabu#clone in
                let Tree.Edge break = 
                    Tree.normalize_edge break tree.Ptree.tree 
                in
                let breakage = break_fn break tree in
                let () = ntabu#update_break breakage in
                let tree, treed = 
                    join_fn breakage.Ptree.incremental 
                    r1 r2 breakage.Ptree.ptree 
                in
                let () = ntabu#update_join tree treed in
                let cost = Ptree.get_cost `Adjusted tree in
                if cost < icost then 
                    tree, ntabu, true
                else acc_tree, acc_tabu, acc_todo
            in
            let what_to_do = 
                List.sort (fun (a, _, _, _) (b, _, _, _) -> compare b a)
                all_places_to_join 
            in
            match current_tree_in_search with
            | Some (tr, tabu) ->
                    let tree, tabu, do_again = 
                        List.fold_left (attempt_to_make_them_true tr tabu) 
                        (tr, tabu, false) what_to_do
                    in
                    do_repeat <- do_again;
                    all_places_to_join <- [];
                    let cost = Ptree.get_cost `Adjusted tree in
                    results <- [(tree, cost, tabu)]
            | None -> ()

        method any_trees =
            self#check_results;
            super#any_trees 

        method next_tree =
            let ((t, _, tabu) as tree)= super#next_tree in
            current_tree_in_search <- Some (t, tabu#clone);
            tree

        method process
            (cost_fn : ('a, 'b) Ptree.cost_fn)
            (b_delta : float)
            (cd_nd : 'a)
            (join_fn : ('a, 'b) Ptree.join_fn)
            incremental
            (j1 : Tree.join_jxn)
            (j2 : Tree.join_jxn) 
            (tabu_mgr : ('a, 'b) Ptree.tabu_mgr)
            pt = 
            incr_trees_considered ();
            match cost_fn j1 j2 b_delta cd_nd pt with
            | Ptree.NoCost -> Tree.Skip
            | Ptree.Cost cc ->
                    tabu_mgr#break_distance cc;
                    sampler#process incremental j1 j2 cd_nd pt None b_delta cc None; 
                    if cc < b_delta then 
                      let nt, j_delta = (join_fn incremental j1 j2 pt) in
                      let cst = Ptree.get_cost `Adjusted nt in
                      if cst < cur_best_cost then 
                        let _ =
                            (* If the gain is better, store it *)
                            let gain = cur_best_cost -. cst in
                            match current_best_gain with
                            | None -> 
                                    current_best_gain <- 
                                        Some (gain, j1, j2)
                            | Some (tmp, _, _) ->
                                    if tmp > gain then ()
                                    else 
                                        current_best_gain <- 
                                            Some (gain, j1, j2)
                        in
                        Tree.Continue
                      else Tree.Continue
                  else Tree.Continue
    end

    (** ['a, 'b] hold_n_threshold_srch_mgr n strategy threshold] is the same as
    above, but keeps a threshold of trees above the current best *)
    class hold_n_threshold_srch_mgr n strategy threshold sampler =
    object (self)
        inherit (hold_n_fb_srch_mgr n strategy sampler) as super

        method private filter t cost =
            cost <= (cur_best_cost +. threshold)

        method features l =
            super#features (("threshold", string_of_float threshold) :: l)
    end

    (** ['a, 'b] annealing_srch_mgr n strategy threshold coeff e] uses simulated
        annealing to control the trajectory.  It uses the standard temperature
        function with coefficient [coeff] and exponent [e]. *)
    class virtual probabilistic_search_mgr n strategy threshold sampler =
    object (self)
        inherit (hold_n_threshold_srch_mgr n strategy threshold sampler) as super

        val mutable tree_count = 0

        method virtual private is_better : float -> float -> bool

        method private trajectory incremental cc cd_nd b_delta j1 j2 pt ljoin ltabu cost =
            tree_count <- succ tree_count;
            if self#is_better cc b_delta then begin
                    let nt, j_delta = Lazy.force ljoin in
                    let cst = Ptree.get_cost `Adjusted nt in
                    if self#is_better cst cur_best_cost then begin
                        cur_best_cost <- cst;
                        sampler#process incremental j1 j2 cd_nd pt (Some nt) b_delta cc 
                        (Some cst);
                        let new_tabu = Lazy.force ltabu in
                        srch_trees <- (nt, cst, new_tabu) :: srch_trees;
                        Tree.Break
                    end else begin
                        sampler#process incremental j1 j2 cd_nd pt None b_delta cc None;
                        Tree.Continue
                    end
                end else begin
                    sampler#process incremental j1 j2 cd_nd pt None b_delta cc None;
                    Tree.Continue
                end
    end

    class annealing_srch_mgr n strategy threshold coeff e sampler =
    object (self)
        inherit (probabilistic_search_mgr n strategy threshold sampler)

        method clone = ({< sampler = sampler#clone >} :> (Node.n, Edge.e) Ptree.search_mgr)

        method features lst = 
            ( "type", "simmulated annealing") :: ("temp.type", "standard") ::
                ("temp.coefficient", string_of_float coeff) :: 
                ("temp.exponent", string_of_float e) :: lst

        method private temp_fn iter = coeff *. (exp (-. (float_of_int iter) /. e))

        method private is_better a b =
            if a < b then true
            else if a = b then false
            else
                let temp = self#temp_fn tree_count in
                if temp = 0. then false
                else
                    let prob = (exp (-. (a -. b) /. temp)) in
                    (Random.float 1.) < prob
    end

    class classic_poy_drifting_srch_mgr n strategy equalprob 
    worstfactor sampler =
        object (self)
            inherit (probabilistic_search_mgr n strategy 0.0 sampler)

        method clone = ({< sampler = sampler#clone >} :> (Node.n, Edge.e) Ptree.search_mgr)

        method features lst = 
            ( "type", "tree drifting") :: ("temp.type", "POY 3.0") ::
                ("temp.equal_cost_prob", string_of_float equalprob) :: 
                ("temp.worst_prob_factor", string_of_float worstfactor) :: lst

        method private is_better a b =
            if a < b then true
            else if a = b then (Random.float 1.) < equalprob
            else (Random.float 1.) < (1. /. (worstfactor +. a -. b))
    end

    class all_possible_joins file sampler = 
        let print_jxn x = 
            "(" ^ 
            (match x with 
            | Tree.Single_Jxn x -> 
                    string_of_int x
            | Tree.Edge_Jxn (x, y) -> 
                    string_of_int x ^ "," ^ string_of_int y) ^ ")"
        in
        let length x pt =
            "(" ^
            match x with
            | Tree.Single_Jxn x -> "Single"
            | Tree.Edge_Jxn (x, y) ->
                    let a = Ptree.get_node_data x pt
                    and b = Ptree.get_node_data y pt in
                    string_of_float (Node.edge_distance a b)
        in
        object
        inherit first_best_srch_mgr sampler

        val mutable current_jxnts = None

        method process 
            (cost_fn : (Node.n, Edge.e) Ptree.cost_fn)
            (b_delta : float)
            (cd_nd : Node.n)
            (join_fn : (Node.n, Edge.e) Ptree.join_fn)
            incremental
            (j1 : Tree.join_jxn)
            (j2 : Tree.join_jxn) 
            (tabu_mgr : (Node.n, Edge.e) Ptree.tabu_mgr)
            pt =
        incr_trees_considered ();
            let clade_cost = (cost_fn j1 j2 b_delta cd_nd pt) in
            let union_cost = 
                match file with
                | Some _ ->
                    let base = 
                        match j1 with
                        | Tree.Edge_Jxn (a, _) 
                        | Tree.Single_Jxn a -> Ptree.get_node_data a pt 
                    in
                    Node.union_distance cd_nd base
                | None -> 0.0
            in
            match clade_cost with
            | Ptree.Cost(cc) ->
                    tabu_mgr#break_distance cc;
                    let _ =
                        match file with
                        | Some file ->
                                let msg = (print_jxn j1) ^ "-" ^ (print_jxn j2) ^ 
                                "=" ^ string_of_float b_delta ^ " - " ^ 
                                string_of_float cc ^ " - " ^ string_of_float union_cost ^ 
                                length j1 pt ^ "\n" in
                                Status.user_message (Status.Output ((Some file), false, [])) msg;
                        | None -> ()
                    in
                    if cc < b_delta then
                        let _ = 
                            sampler#process incremental j1 j2 cd_nd pt 
                            None b_delta cc None 
                        in
                        let _ = 
                            match file with
                            | Some file ->
                                    let msg = 
                                        "Joint with cost_delta " ^ string_of_float
                                        (b_delta -. cc)
                                    in
                                    Status.user_message (Status.Output ((Some file), false, []))
                                    msg
                            | None -> ()
                        in
                        ()
                    else sampler#process incremental j1 j2 cd_nd pt None b_delta cc None;
                    Tree.Continue
            | Ptree.NoCost -> Tree.Skip

        method clone = new all_possible_joins file sampler#clone
    end

    class supports_manager best_vals node_indices starting status  nbhood_count orig_cost' sampler = 
        object
        inherit first_best_srch_mgr sampler as super

        val nbhood_report = 500 (* The module of the report update *)
        method process cost_fn b_delta cd_nd join_fn incremental j1 j2
            tabu_mgr tree =
            let cost = match cost_fn j1 j2 infinity cd_nd tree with
            | Ptree.NoCost -> infinity
            | Ptree.Cost c -> orig_cost' +. c -. b_delta in
            (*** update values of all clades not in this tree *)
            (* update the status *)
            incr nbhood_count;
            if (!nbhood_count mod nbhood_report) = 0
            then
                Status.full_report
                    ~msg:("Examining tree with cost " ^ string_of_float cost)
                    ~adv:(!nbhood_count)
                    status;
            (* for now, we restrict our attention to a single component.  this
               is insufficiently general.  TODO: fix.(?) *)
            let root_median =
                match j1 with
                | Tree.Single_Jxn h ->
                      Node.median None None (Ptree.get_node_data h tree) cd_nd
                | Tree.Edge_Jxn (h, n) ->
                      let j1median = Node.median None None 
                                        (Ptree.get_node_data h tree)
                                        (Ptree.get_node_data n tree) in
                      Node.median None None j1median cd_nd
            in
            (* run through each character... *)
            (* (extract those we've created) *)
            let char_list = Node.support_chars starting None root_median in
            List.iter
                (fun (cost, charcode) ->
                     (* clade not present iff cost > 1 *)
                     if cost > 1.  then begin
                         try
                             let index = List.assoc charcode node_indices in
                             if cost < best_vals.(index)
                             then best_vals.(index) <- cost
                         with 
                         | Not_found -> () (* this is for non-clade characters *)
                     end
                )
                char_list;

            (* continue *)
            Tree.Continue
    end
end 
