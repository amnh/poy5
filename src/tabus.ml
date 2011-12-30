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

(** tabus.ml : This is the file where different kinds of tabu_mgr objects are
 * implemented. The tabu manager specifies the order in which edges are broken by
 * the SPR and TBR search procedures. The list of edges in the tabu should always
 * match the edges in the tree. *)
let () = SadmanOutput.register "Tabus" "$Revision: 2871 $"

(* A module that provides the managers for a local search (rerooting, edge
* breaking and joining. A tabu manager controls what edges are next ina series
* of requests for rerooting, breaking and joining *)

(** The type produced by the {!Tabus.Make} functor *)
module type S = sig

    (** {2 Types} *)

    (** The type of the data of a vertex in a tree *)
    type a 

    (** The type of the data of an edge in a tree *)
    type b

    (** A phylogenetic tree, holding data of type {!a} and {!b} *)
    type phylogeny = (a, b) Ptree.p_tree

    (** A wagner edge manager with concrete types a and b in the tree *)
    class type wem = [a, b] Ptree.wagner_edges_mgr

    (** A function that constructs a wagner edge manager for the phylogeny. *)
    val wagner_tabu : phylogeny -> int -> wem

    val wagner_constraint : All_sets.Integers.t -> phylogeny -> int -> wem

    (* A function that constructs a wagner edge manager for the phylogeny, that
     * uses the unions of the vertices to reduce the number of attempts. *)
    val wagner_union : phylogeny -> int -> wem

    val distance_dfs_wagner : phylogeny -> int -> wem

    (** {2 Edge and 2 Node Manager Constructors} *)

    (** An edge/node manager that handles a particular strategy to process the edges of a tree *)
    type emc 
    class type nm = [a, b] Ptree.nodes_manager

    (** An edge/node manager for a side of a tree *)
    type semc

    (** An edge manager that uses the union of the vertices contained in a
    * subtree to decide weather or not it is necessary to continue down that
    * path. The edge manager will traverse up to distance [int] from the last
    * breaking position. It uses Depth First Search to select the next edge to
    * evaluate. *)
    val union_join : int -> semc

    (** An edge manager that returns all the edges that are at distance up to
    * [int] from the last breaking position *)
    val distance_join : int -> semc

    (** An edge manager that will traverse in a using depth first search
    * starting in the root of the tree, but will not pass the vertices included
    * in the argument set. The procedure will also return edges that occur at
    * distance [int] from the last break position. *)
    val partitioned_join : 
        [ `Sets of All_sets.IntSet.t | `Height of int ] -> int -> semc

    (** A simple depth first search manager to choose the next rerooting
    * position. The manager will return edges up to distance [int] from the last
    * break position *)
    val reroot : 
        [ `Sets of All_sets.IntSet.t | `Height of int ] option -> int -> semc

    (** Choose at random any non previously selected edge *)
    val random_break : emc

    (** Choose the edges to break, starting with those that show greates 
    * length.  *)
    val sorted_break : 
        bool -> [ `Sets of All_sets.IntSet.t | `Height of int ] option ->  emc

    (** Break an edge once and only once, never again *)
    val only_once_break : emc

    (** iteration managers -- from doing nothing to everything and something
     * between. The first two arguments define the model iteration pattern *)
    val simple_nm_none : int option -> float option -> nm
    val simple_nm_all : int option -> float option -> nm
    val complex_nm_delta : int option -> float option -> nm
    val complex_nm_neighborhood : int option -> float option -> nm

    (** {2 Tabu Managers} *)

    class type tabu_mgr = [a, b] Ptree.tabu_mgr

    (** Simple tabu manager to direct the search on a forest. *)
    val join_to_tree_in_forest :
        phylogeny -> int -> tabu_mgr * All_sets.Integers.elt

    (* A standard tabu manager that composes three edge managers on a starting
    * phylogeny: the join, reroot, iterate and break managers. *)
    class standard_tabu : phylogeny -> semc -> semc -> emc -> nm -> tabu_mgr

end


let timer_interval = ref 1

module Make  (Node : NodeSig.S) (Edge : Edge.EdgeSig with type n = Node.n) : S with type a = Node.n with type b = Edge.e = struct

    type a = Node.n
    type b = Edge.e
    type phylogeny = (Node.n, Edge.e) Ptree.p_tree


    let ndebug = false
    let debug_num_edges = false             (* Print how many edges are in the
                                               to-break and to-join queues 
                                               whenever an update is called *)
    let debug_tabu_join_once = false
    let debug_join_to_tree_in_forest = false
    let odebug = Status.user_message Status.Information

    let data_mining = false

    let (-->) a b = b a

    let print_time timer current_time = 
        if !timer_interval = 0 then 0
        else
            let wall = Timer.wall timer in
            let nt = int_of_float wall in
            if Status.is_interactive () then
                if !timer_interval <= nt - current_time then begin
                    Status.user_message Status.SearchReport 
                    (SearchInformation.show_information 
                    None 
                    None 
                    (Some (string_of_int nt))
                    (Some [`Timer]));
                    nt
                end else current_time
            else 0

    class type ['a, 'b] nodes_manager = object
        method update_iterate :
            ('a, 'b) Ptree.p_tree -> 
            ([ `Break of ('a, 'b) Ptree.breakage 
             | `Join of Tree.join_delta 
             | `Reroot of Ptree.incremental list
             | `Cost ]) -> unit
        method clone : ('a, 'b) nodes_manager
        method branches : Tree.edge list option
        method model : bool
        method to_string : string 
        method fuse : ('a,'b) nodes_manager -> ('a,'b) nodes_manager 
    end

    class virtual ['a, 'b] tabu_base = object (self)
        val my_features : (string * string) list = []

        method private inherited_features (l : (string * string) list) = l

        method features l = self#inherited_features l

        method virtual update_break : ('a, 'b) Ptree.breakage -> unit

        method update_reroot (_ :  ('a, 'b) Ptree.breakage) = ()

        method virtual update_join :
            ('a, 'b) Ptree.p_tree -> Tree.join_delta -> unit

        method virtual break_edge : Tree.edge option

        method virtual join_edge : [`Left | `Right] -> Tree.edge option

        method virtual break_edges : Tree.edge list

        method reroot_edge a = self#join_edge a

        method break_distance (_ : float) = ()

        method virtual get_node_manager : ('a,'b) nodes_manager option
    end

    let rand arr = 
        let len = Array.length arr in
        let random pos _ = 
            if pos < len - 1 then begin
                let other = pos + (Random.int (len - 1 - pos)) in
                let tmp = arr.(pos) in
                arr.(pos) <- arr.(other);
                arr.(other) <- tmp;
            end else ()
        in
        Array.iteri random arr

    let randomize lst =
        let arr = Array.of_list lst in
        rand arr;
        Array.to_list arr

    class type ['a, 'b] edges_manager = object
        method break_distance : float -> unit
        method next_edge : Tree.edge option
        method update_break   : ('a, 'b) Ptree.breakage -> unit
        method update_reroot  : ('a, 'b) Ptree.breakage -> unit
        method update_join    : ('a, 'b) Ptree.p_tree -> Tree.join_delta -> unit
        method clone : ('a, 'b) edges_manager
        method exclude : Tree.edge list -> unit
    end

    type tip_path = (int list * int)
    type middle_path = (tip_path * tip_path * int * int)

    let rec longest_path tree vertex = 
        match Tree.get_node vertex tree with
        | Tree.Interior (_, _, l, r) ->
                let (ltp, llen), (ltp1, ltp2, lcode, llenm) = longest_path tree l
                and (rtp, rlen), (rtp1, rtp2, rcode, rlenm) = 
                    longest_path tree r
                in
                let merged =
                    if (llen + rlen + 1) > (max llenm rlenm) then
                        (ltp, rtp, vertex, llenm + rlen + 1)
                    else if llenm > rlenm then
                        (ltp1, ltp2, lcode, llenm)
                    else 
                        (rtp1, rtp2, rcode, rlenm)
                in
                let path =
                    if llen > rlen then
                        (vertex :: ltp, llen + 1)
                    else 
                        (vertex :: rtp, rlen + 1)
                in
                path, merged
        | _ -> ([vertex], 0), ([], [], vertex, 0)

    let longest_path tree = 
        let vertex = All_sets.Integers.choose (Tree.get_handles tree) in
        match Tree.get_node vertex tree with
        | Tree.Interior (_, b, c, d) ->
                let (btp, blen), (btp1, btp2, bcode, blenm) = longest_path tree b
                and (ctp, clen), (ctp1, ctp2, ccode, clenm) = longest_path tree c
                and (dtp, dlen), (dtp1, dtp2, dcode, dlenm) = longest_path tree d
                in
                let (_, _, _, len) as res =
                    if blenm > max clenm dlenm then
                        (btp1, btp2, bcode, blenm)
                    else if clenm > dlenm then
                        (ctp1, ctp2, ccode, clenm)
                    else 
                        (dtp1, dtp2, dcode, dlenm)
                in
                let (lp1, llen1), (lp2, llen2) = 
                    if blen > max clen dlen then
                        if clen > dlen then (btp, blen), (dtp, dlen)
                        else (btp, blen), (dtp, dlen)
                    else if clen > max blen dlen then
                        if blen > dlen then (ctp, clen), (btp, blen)
                        else (ctp, clen), (dtp, dlen)
                    else 
                        if clen > blen then (dtp, dlen), (ctp, clen)
                        else (dtp, dlen), (btp, blen)
                in
                if len > (llen2 + llen1 + 1) then res
                else (lp1, lp2, vertex, llen2 + llen1 + 1)
        | Tree.Leaf (_, b) ->
                let (btp, blen), (btp1, btp2, bcode, blenm) = longest_path tree b
                in
                if blenm > blen + 1 then
                    (btp1, btp2, bcode, blenm)
                else ([], btp, vertex, blen + 1)
        | Tree.Single _ -> failwith "Single once again?"

    let middle_of_path (a, b, vertex, len) =
        let rec find_in_list lst pos =
            match pos, lst with
            | 1, (a :: b :: _) 
            | 0, (a :: b :: _) -> a
            | n, (_ :: t) -> find_in_list t (n - 1)
            | _, _ -> failwith "Empty list?"
        in
        let la = List.length a in
        let lh = len / 2 in
        if la < lh then 
            find_in_list (vertex :: b) (lh - la)
        else if la = lh then vertex
        else find_in_list (vertex :: a) (lh - List.length b)

    let undirected_edge (a, b) (Tree.Edge (c, d)) =
        (a = c && b = d) || (a = d && b = c)

    class exclude_edges = object
        val mutable exclude_edges = []

        method exclude (edges : Tree.edge list) =
            exclude_edges <- edges

        method private add_to_stack ((Tree.Edge (a, b), (_ : int)) as item) 
        stack =
            if List.exists (undirected_edge (a, b)) exclude_edges then ()
            else Stack.push item stack

    end

    class type em = [Node.n, Edge.e] edges_manager
    class type nm = [Node.n, Edge.e] nodes_manager
    class type tabu_mgr = [Node.n, Edge.e] Ptree.tabu_mgr
    class type wem = [Node.n, Edge.e] Ptree.wagner_edges_mgr

    class tabu_join_once left right : [Node.n, Edge.e] Ptree.tabu_mgr =
        let get_list l =
            match l with
            | [] -> None, []
            | l :: ls -> Some l, ls in
        let print e = (match e with
        | Some Tree.Edge(e1, e2) ->
            if debug_tabu_join_once
            then odebug ("Returning edge " ^ string_of_int e1
                         ^ " -> " ^ string_of_int e2)
        | None -> ());
              e in
    object (self)
        inherit [Node.n, Edge.e] tabu_base as super

        val mutable left = left
        val mutable right = right

        method clone = ({<>} :> (Node.n, Edge.e) Ptree.tabu_mgr)

        method private clear =
            left <- []; right <- []

        method private get_left =
            let r, list = get_list left in
            left <- list;
            print r

        method private get_right =
            let r, list = get_list right in
            right <- list;
            print r

        method features l = ("tabu", "join-once") :: l

        method update_break _ =
            self#clear;
            ()

        method update_join _ _ =
            self#clear;
            ()

        method get_node_manager = None

        method break_edge = None

        method join_edge = function
            | `Left -> self#get_left
            | `Right -> self#get_right

        method break_edges = []

    end

    let join_to_tree_in_forest forest =
        let handles =
            All_sets.Integers.elements (Ptree.get_handles forest) in
        fun num ->
            (* We need to make `Left and `Right lists.  The specified tree becomes
               the `Right, and we choose a junction from it for j2.  We then use
               all the other trees for `Left. *)
            let right, left = Utl.remove_nth handles num in
            let left_edges =
                List.map (fun i -> Ptree.get_pre_order_edges i forest) left in
            let left_edges = List.fold_left List.rev_append [] left_edges in
                (* TODO: remove right junction from the list of `Right edges... *)
            let right_edges = Ptree.get_pre_order_edges right forest in 
            (new tabu_join_once left_edges right_edges, right)

    class virtual wagner_dfs_distance_based constraints max_distance (ptree : phylogeny) 
        (handle : int) = 
            let to_calculate_for_compare =
                let to_calculate_for_compare = Stack.create () in
                match Ptree.get_node handle ptree with
                | Tree.Interior (s1, s2, _, _)
                | Tree.Leaf (s1, s2) ->
                        let add_edges parent child =
                            let node = Tree.get_node child ptree.Ptree.tree in
                            try
                                let (a, b) = Tree.other_two_nbrs parent node in
                                Stack.push (Tree.Edge (child, a), 1) to_calculate_for_compare;
                                Stack.push (Tree.Edge (child, b), 1) to_calculate_for_compare
                            with
                            | _ -> ()
                        in
                        if All_sets.Integers.mem s2 constraints || 
                            All_sets.Integers.mem s1 constraints then ()
                        else add_edges s2 s1;
                        Stack.push (Tree.Edge (s1, s2), 0) to_calculate_for_compare;
                        to_calculate_for_compare
                | Tree.Single _ -> failwith "Joining with a single?"
            in
            object (self)

        val timer = Timer.start ()
        val mutable current_time = 0
        val to_compare = Stack.create ()
        val to_calculate_for_compare = to_calculate_for_compare
        val to_do_later = Stack.create ()
        val mutable currently_doing = None
        val mutable current_handle = handle
        val mutable current_depth = 0
        val mutable current_joined_ptree = ptree
        val mutable current_broken_ptree = ptree
        val mutable current_delta = max_float
        val mutable current_clade : Node.n option = None

        method next_clade x = 
            current_clade <- Some x

        method private initialize_children_edges child parent depth =
            let node = Tree.get_node child current_broken_ptree.Ptree.tree in
            try
                let new_depth = depth + 1 in
                let (a, b) = Tree.other_two_nbrs parent node in
                Stack.push (Tree.Edge (child, a), new_depth) to_calculate_for_compare;
                Stack.push (Tree.Edge (child, b), new_depth) to_calculate_for_compare;
            with
            | _ -> ()

        method update_join (tree : phylogeny) (_ : Tree.join_delta) = 
            current_joined_ptree <- tree;
            let handle, _ = 
                Tree.get_path_to_handle current_handle tree.Ptree.tree
            in
            current_handle <- handle;
            current_delta <- 0.0;
            current_broken_ptree <- tree;
            current_delta <- max_float;
            self#clear_all_joins;
            match Ptree.get_node current_handle tree with
            | Tree.Interior (s1, s2, _, _)
            | Tree.Leaf (s1, s2) ->
                    self#initialize_children_edges s1 s2 0;
                    Stack.push (Tree.Edge (s1, s2), 0) to_calculate_for_compare
            | Tree.Single _ -> failwith "Joining with a single?"
            
        method private clear_all_joins =
            Stack.clear to_compare;
            Stack.clear to_calculate_for_compare;
            Stack.clear to_do_later

        method private select_existant_and_update =
            if Stack.is_empty to_do_later then None
            else 
                let edge, depth = Stack.pop to_do_later in
                let Tree.Edge (a, b) = edge in
                let _ = self#initialize_children_edges b a depth in
                self#next_edge

        method private sort_and_choose_join_edge =
            let rec get_all_elements acc =
                if Stack.is_empty to_compare then acc
                else 
                    get_all_elements ((Stack.pop to_compare) :: acc)
            in
            let lst = get_all_elements [] in
            let sorter (_, _, (a : float)) (_, _, (b : float)) = compare b a in
            let lst = List.sort sorter lst in
            List.iter (fun (x, y, _) -> Stack.push (x, y) to_do_later) lst;
            self#select_existant_and_update 

        val mutable has_potential = false

        method virtual private should_continue_this_path : Tree.edge option -> bool -> bool

        method next_edge = 
            current_time <- print_time timer current_time;
            has_potential <- false;
            let edge, do_check, depth = 
                if Stack.is_empty to_calculate_for_compare then 
                    if Stack.is_empty to_compare then
                        if Stack.is_empty to_do_later then None, false, 0
                        else self#select_existant_and_update, false, 0
                    else self#sort_and_choose_join_edge, false, 0
                else begin 
                    has_potential <- true;
                    let edge, depth = Stack.pop to_calculate_for_compare in
                    Some edge, true, depth
                end
            in
            if depth > max_distance then self#next_edge
            else begin
                if do_check then begin
                    if self#should_continue_this_path edge has_potential then begin
                        currently_doing <- edge;
                        current_depth <- depth;
                    end else begin
                        currently_doing <- None;
                        current_depth <- 0;
                    end
                end else ();
                edge
            end

        method break_distance (x : float) =
            if current_delta > x then 
                current_delta <- x
            else ();
            begin match currently_doing with
            | Some edge -> 
                    Stack.push (edge, current_depth, x) to_compare;
            | None -> ()
            end;
            currently_doing <- None

        method exclude (_ : Tree.edge list) = ()

    end

    class union_dfs_wagner max_distance ptree handle : wem = object

        inherit wagner_dfs_distance_based All_sets.Integers.empty max_distance ptree handle

        method new_delta _ = ()

        method clone =  ({<
            to_compare = Stack.copy to_compare;
            to_calculate_for_compare = Stack.copy to_calculate_for_compare;
            to_do_later = Stack.copy to_do_later;
             >} :> wem)

        method private should_continue_this_path edge x =
            x &&
                (match edge with
                | None -> false
                | Some (Tree.Edge (a, b)) ->
                        let tree = current_broken_ptree.Ptree.tree in
                        if a = Tree.get_parent b tree then
                            if Tree.is_leaf b tree then false
                            else 
                                let clade =
                                    match current_clade with
                                    | Some clade -> clade
                                    | None -> failwith "No root?"
                                in
                                let ca = 
                                    Ptree.get_node_data a current_broken_ptree 
                                in
                                let union_delta = Node.union_distance clade ca in
                                if current_delta >= union_delta then true
                                else false
                        else true)
    end

    class constrained_dfs_wagner constrains max_distance ptree handle : wem =
        object
            inherit wagner_dfs_distance_based constrains max_distance ptree
            handle as super

            method new_delta _ = ()
        method clone =  ({<
            to_compare = Stack.copy to_compare;
            to_calculate_for_compare = Stack.copy to_calculate_for_compare;
            to_do_later = Stack.copy to_do_later;
             >} :> wem)

            method private should_continue_this_path edge x =
                x &&
                    (match edge with
                    | None -> false
                    | Some (Tree.Edge (a, b)) ->
                            if All_sets.Integers.mem b constrains then 
                                false
                            else if All_sets.Integers.mem a constrains then
                                false
                            else true)
        end

    class probabilitisic_union_dfs_wagner max_distance ptree handle 
        : wem = object 
        inherit wagner_dfs_distance_based All_sets.Integers.empty max_distance ptree handle

        method new_delta _ = ()

        method clone = ({<
            to_compare = Stack.copy to_compare;
            to_calculate_for_compare = Stack.copy to_calculate_for_compare;
            to_do_later = Stack.copy to_do_later;
        >} :> wem)

        method private should_continue_this_path edge x =
            x &&
                (match edge with
                | None -> false
                | Some (Tree.Edge (a, b)) ->
                        let tree = current_broken_ptree.Ptree.tree in
                        if a = Tree.get_parent b tree then
                            if Tree.is_leaf b tree then false
                            else 
                                if Random.bool () then
                                    let clade =
                                        match current_clade with
                                        | Some clade -> clade
                                        | None -> failwith "No root?"
                                    in
                                    let ca = Ptree.get_node_data a current_broken_ptree in
                                    let union_delta = Node.union_distance clade ca in
                                    if current_delta >= union_delta then true
                                    else false
                                else true
                        else true)
    end 

    class distance_dfs_wagner max_distance ptree handle : wem = object (self)

        inherit wagner_dfs_distance_based All_sets.Integers.empty max_distance ptree handle

        method new_delta _ = ()

        method private should_continue_this_path (y : Tree.edge option) x = 
            if not x then x
            else 
                match y with
                | None -> false
                | _ -> true

        method clone = ({<
            to_compare = Stack.copy to_compare;
            to_calculate_for_compare = Stack.copy to_calculate_for_compare;
            to_do_later = Stack.copy to_do_later;
        >} :> wem)

    end

    (* root class; only updates the model using threshold and count. subclasses
     * must provide branch functionality *)
    class virtual nm_simple_base mcount mthreshold = object (self)

            val threshold = mthreshold
            val maxcount  = mcount
            val mutable count = 0
            val mutable score = None
            val mutable iterate = false

            method private update_model (tree: ('a,'b) Ptree.p_tree)
                    (fn : ([ `Break of ('a, 'b) Ptree.breakage 
                           | `Join of Tree.join_delta 
                           | `Reroot of Ptree.incremental list
                           | `Cost ]) ) = match fn with
                    | `Break _ | `Reroot _ | `Cost -> ()
                    | `Join _ -> begin
                        let new_score = Ptree.get_cost `Adjusted tree
                        and _count () = match maxcount with
                            | Some c when count > c -> count <- 0; true
                            | Some _ | None   -> false
                        in
                        count <- count + 1;
                        match score, threshold with
                            | Some x, Some t when x *. (1.0 +. t) < new_score ->
                                iterate <- true;
                                score <- Some new_score;
                            | None, _ ->
                                iterate <- _count ();
                                score <- Some new_score;
                            | Some x, _ ->
                                iterate <- _count ()
                    end

            method model    =
                if not iterate then false
                else begin iterate <- false; true end

            method clone = ({< >} :> (Node.n, Edge.e) nodes_manager)

            method fuse _ = self#clone

            (* these much be implemented *)
            method virtual branches       : Tree.edge list option
            method virtual update_iterate : ('a, 'b) Ptree.p_tree ->
                                            ([ `Break of ('a, 'b) Ptree.breakage 
                                             | `Join of Tree.join_delta 
                                             | `Reroot of Ptree.incremental list
                                             | `Cost ]) -> unit
            method virtual to_string : string
        end

    class nm_simple_none mcount mthreshold : [Node.n,Edge.e] nodes_manager = object (self)
            inherit nm_simple_base mcount mthreshold as super

            method branches = Some []
            method update_iterate t m = self#update_model t m
            method clone = ({< >} :> (Node.n, Edge.e) nodes_manager)
            method to_string =
                let opt_ to_str = function | Some x -> to_str x | None -> "none" in
                Printf.sprintf 
                    "Nodes Manager:\n\tmodel- counts: (%s,%d) threshold:(%s,%s)\n\tBranches:%s\n%!"
                    (opt_ string_of_int mcount) count (opt_ string_of_float mthreshold)
                    (opt_ string_of_float score) "None"
(*
            method fuse x = 
                let r = x#clone in
                let r#count = x#count + self#count in
                let r#score = min self#count x#score in
                let r#iterate = x#iterate || self#iterate in
                r
*)
    end

    class nm_simple_all mcount mthreshold : [Node.n,Edge.e] nodes_manager = object (self)
            inherit nm_simple_base mcount mthreshold as super

            method branches = None
            method update_iterate t m = self#update_model t m
            method clone = ({< >} :> (Node.n, Edge.e) nodes_manager)
            method to_string =
                let opt_ to_str = function | Some x -> to_str x | None -> "none" in
                Printf.sprintf 
                    "Nodes Manager:\n\tmodel- counts: (%s,%d) threshold:(%s,%s)\n\tBranches:%s\n%!"
                    (opt_ string_of_int mcount) count (opt_ string_of_float mthreshold)
                    (opt_ string_of_float score) "Always"
    end

    class nm_complex_delta mcount mthreshold : [Node.n,Edge.e] nodes_manager = object (self)
            inherit nm_simple_base mcount mthreshold as super

            val mutable b_list = None

            method private edges_of_path tree delta = 
                let make_edge a b = 
                    let new_edge = Tree.Edge (a,b) in
                    if Ptree.is_edge new_edge tree then
                        new_edge
                    else begin
                        let oth_edge = Tree.Edge (b,a) in
                        assert( Ptree.is_edge oth_edge tree );
                        oth_edge
                    end
                in
                let rec loop_ acc p_id oth_ids = match oth_ids with
                    | []    -> acc
                    | x::xs -> loop_ ((make_edge p_id x) :: acc) x xs
                in match delta with
                | [_] | [] -> []
                | x::xs -> loop_ [] x xs

            method update_iterate tree fn = match fn with
                | `Break _ | `Cost | `Reroot _ -> ()
                | `Join (_,_,delta) ->
                    let () = self#update_model tree fn in
                    b_list <- Some (self#edges_of_path tree delta)

            method branches = match b_list with
                | Some [] -> b_list (* not iterating anything *)
                | Some x -> 
(*                    Printf.printf "Iterating Edges: ";*)
(*                    List.iter (fun (Tree.Edge (a,b)) ->*)
(*                                    Printf.printf "%d -- %d, " a b) x;*)
(*                    print_newline ();*)
                    let copy = x in
                    b_list <- Some [];
                    Some copy
                | None -> 
(*                    Printf.printf "Iterating Edges: All\n%!";*)
                    None

            method to_string =
                let opt_ to_str = function | Some x -> to_str x | None -> "none" in
                Printf.sprintf 
                    "Nodes Manager:\n\tmodel- counts: (%s,%d) threshold:(%s,%s)\n\tBranches:%s\n%!"
                    (opt_ string_of_int mcount) count (opt_ string_of_float mthreshold)
                    (opt_ string_of_float score) "Join Delta"
    end

    class nm_complex_neighborhood mcount mthreshold : [Node.n,Edge.e] nodes_manager = object (self)
            inherit nm_simple_base mcount mthreshold as super
    
            val mutable b_list = []

            method private edges_of_delta tree deltal deltar = 
                let rec make_edge a b = 
                    let edge =
                        let new_edge = Tree.Edge (a,b) in
                        if Ptree.is_edge new_edge tree then
                            new_edge
                        else begin
                            let oth_edge = Tree.Edge (b,a) in 
                            assert( Ptree.is_edge oth_edge tree );
                            oth_edge
                        end
                    in
                    edge
                and get_edges_on_side = function
                    | `Single (i,_)      -> [], i
                    | `Edge (i,l,r,_)  -> [(make_edge i l); (make_edge i r)], i
                in
                let acc1,il = get_edges_on_side deltal
                and acc2,ir = get_edges_on_side deltar in
                (make_edge il ir) :: (acc1 @ acc2)

            method update_iterate tree fn = match fn with
                | `Break _ | `Cost | `Reroot _ -> ()
                | `Join (d1,d2,_) ->
                    let () = self#update_model tree fn in
                    b_list <- self#edges_of_delta tree d1 d2

            method branches = Some b_list

            method to_string =
                let opt_ to_str = function | Some x -> to_str x | None -> "none" in
                Printf.sprintf 
                    "Nodes Manager:\n\tmodel- counts: (%s,%d) threshold:(%s,%s)\n\tBranches:%s\n%!"
                    (opt_ string_of_int mcount) count (opt_ string_of_float mthreshold)
                    (opt_ string_of_float score) "Join Neighborhood"
    end

    class simple_dfs (ptree : (Node.n, Edge.e) Ptree.p_tree) : 
        [Node.n, Edge.e] edges_manager = object (self)

            val mutable starting = true
            val to_do_later = Stack.create ()
            val mutable tree = ptree.Ptree.tree
            val mutable exclude_edges = []

            method update_join ptree _ =
                tree <- ptree.Ptree.tree;
                starting <- true

            method update_break _ = ()
            method update_reroot _ = ()

            method clone = ({< >} :> (Node.n, Edge.e) edges_manager)

            method break_distance _ = ()

            method private add_both_children ((parent, child) as e) =
                if List.exists (undirected_edge e) exclude_edges then ()
                else
                    try
                        let node = Tree.get_node child tree in
                        let (a, b) = Tree.other_two_nbrs parent node in
                        Stack.push (child, a) to_do_later;
                        Stack.push (child, b) to_do_later;
                    with
                    | _ -> ()

            method next_edge =
                if starting then
                    let vertex = 
                        let edges = tree.Tree.d_edges in
                        let choose_rnd item ((pos, _) as res) =
                            if pos = 0 then (-1, Some item)
                            else if pos > 0 then (pos - 1, None)
                            else res
                        in
                        let item = Random.int (Tree.EdgeSet.cardinal edges) in
                        match Tree.EdgeSet.fold choose_rnd edges (item, None) with
                        | _, Some (Tree.Edge (x, _)) -> x
                        | _, None -> failwith "No edges?"
                    and _ = starting <- false in
                    let _ =
                        match Tree.get_node vertex tree with
                        | Tree.Interior (v, a, b, c) ->
                                Stack.push (v, a) to_do_later;
                                Stack.push (v, b) to_do_later;
                                Stack.push (v, c) to_do_later;
                        | Tree.Leaf (v, a) ->
                                Stack.push (v, a) to_do_later;
                        | Tree.Single _ -> failwith "Single again?"
                    in
                    self#next_edge
                else
                    let (a, b) as res = Stack.pop to_do_later in
                    let _ = self#add_both_children res in
                    Some (Tree.Edge (a, b))

            method exclude edges = 
                exclude_edges <- edges

    end

    let get_side side breakage = 
        match side with
        | `Left -> breakage.Ptree.left.Ptree.topology_delta
        | `Right -> breakage.Ptree.right.Ptree.topology_delta
        | `Both -> failwith "Both?"

    let get_both side delta = 
        match side, delta with
        | `Both, (l, r) -> [l; r]
        | _ -> failwith "No both?"


    module type S = sig
        type 'a t
        val pop : 'a t -> 'a
        val push : 'a -> 'a t -> unit
        val create : unit -> 'a t
        val is_empty : 'a t -> bool
        val copy : 'a t -> 'a t 
        val clear : 'a t -> unit
        val iter : ('a -> unit) -> 'a t -> unit
    end

    module SQ (M : S) = struct
        class generic_distance side max_distance ptree : [Node.n, Edge.e]
        edges_manager = object (self) 
            val items = M.create ()
            val mutable ptree = ptree
            val mutable exclude = []

            method break_distance _ = ()
            method next_edge =
                let debug = false in
                if debug then begin
                    Printf.printf "SQ contains\n%!";
                    M.iter (fun ((Tree.Edge (a, b)), d) ->
                            Printf.printf "%d,%d - %d\n%!" a b d) items;
                end;
                if M.is_empty items then None
                else 
                    let ((Tree.Edge (a, b)) as e, depth) = M.pop items in
                    if depth < max_distance then begin
                        let node = Ptree.get_node a ptree in
                        match node with
                        | Tree.Interior _ -> 
                                let newdepth = depth + 1 in
                                let c, d = Ptree.other_two_nbrs b node in
                                M.push ((Tree.Edge (c, a)), newdepth) items;
                                M.push ((Tree.Edge (d, a)), newdepth) items;
                        | Tree.Leaf _ | Tree.Single _ -> ()
                    end;
                    Some e

            method update_break breakage = 
                ptree <- breakage.Ptree.ptree;
                let clade =
                    match side with
                    | `Left -> breakage.Ptree.left
                    | `Right -> breakage.Ptree.right
                in
                M.clear items;
                match clade.Ptree.topology_delta with
                | `Single _ -> ()
                | `Edge (_, a, b, _) ->
                        let handle = Ptree.handle_of a breakage.Ptree.ptree in
                        match Ptree.get_node handle breakage.Ptree.ptree with
                        | Tree.Interior (a, b, c, d) ->
                                M.push ((Tree.Edge (b, a)), 0) items;
                                M.push ((Tree.Edge (c, a)), 1) items;
                                M.push ((Tree.Edge (d, a)), 1) items;
                        | Tree.Leaf (a, b) -> M.push ((Tree.Edge (b, a)), 0) items;
                        | Tree.Single _ -> ()

            method update_reroot _  = ()
            method update_join _ _ = ()

            method clone = 
                ({< items = M.copy items;
                ptree = ptree; >} :> (Node.n, Edge.e) edges_manager)
            method exclude edges = exclude <- edges
        end
    end

    module DFS = SQ (Stack)
    module BFS = SQ (Queue)

    class virtual dfs_distance_based side max_distance 
        (initial_ptree : (Node.n, Edge.e) Ptree.p_tree) = object (self)
            inherit exclude_edges as super

        val to_compare = Stack.create ()
        val to_calculate_for_compare = Stack.create ()
        val to_do_later = Stack.create ()
        val mutable currently_doing = None
        val mutable current_depth = 0
        val mutable current_joined_ptree = initial_ptree
        val mutable current_broken_ptree = initial_ptree
        val mutable current_delta = 0.0
        val mutable current_clade = None

        method update_reroot (_  : (Node.n, Edge.e) Ptree.breakage) = ()

        method update_join (tree : (Node.n, Edge.e) Ptree.p_tree) 
            (_ : Tree.join_delta) =
            current_joined_ptree <- tree;
            current_delta <- 0.0

        method private initialize_children_edges child parent depth =
            let node = Tree.get_node child current_broken_ptree.Ptree.tree in
            try
                let new_depth = depth + 1 in
                let (a, b) = Tree.other_two_nbrs parent node in
                super#add_to_stack (Tree.Edge (child, a), new_depth)
                to_calculate_for_compare; 
                super#add_to_stack (Tree.Edge (child, b), new_depth) 
                to_calculate_for_compare;
            with
            | _ -> ()

        method private clear_all_joins =
            Stack.clear to_compare;
            Stack.clear to_calculate_for_compare;
            Stack.clear to_do_later

        method update_break breakage =
            let ptree = breakage.Ptree.ptree in
            let delta = get_side side breakage in
            current_broken_ptree <- ptree;
            current_clade <- 
                Some (match side with
                (* We have to select the oposite side *)
                | `Left -> breakage.Ptree.right.Ptree.clade_node
                | `Right -> breakage.Ptree.left.Ptree.clade_node);
            current_delta <- 
                (Ptree.get_cost `Adjusted current_joined_ptree) -. 
                (Ptree.get_cost `Adjusted current_broken_ptree);
            self#clear_all_joins;
            match delta with
            | `Single (_, _) -> ()
            | `Edge (r, s1, s2, _) ->
                    let add x = 
                        let e = Tree.Edge (s1, x) in
                        if self#should_continue_this_path (Some e) true then
                            let d = if x = s2 then 0 else 1 in
                            super#add_to_stack (e, d) 
                            to_calculate_for_compare 
                        else ()
                    in
                    match Ptree.get_node s1 breakage.Ptree.ptree with
                    | Tree.Interior (s1, a, b, c) ->
                            add a;
                            add b;
                            add c;
                    | Tree.Leaf (_, a) -> add a;
                    | Tree.Single _ -> ()

        method private select_existant_and_update =
            if Stack.is_empty to_do_later then None
            else 
                let edge, depth = Stack.pop to_do_later in
                let Tree.Edge (a, b) = edge in
                let _ = self#initialize_children_edges b a depth in
                self#next_edge

        method private sort_and_choose_join_edge =
            let rec get_all_elements acc =
                if Stack.is_empty to_compare then acc
                else 
                    get_all_elements ((Stack.pop to_compare) :: acc)
            in
            let lst = get_all_elements [] in
            let sorter (_, _, (a : float)) (_, _, (b : float)) = compare b a in
            let lst = List.sort sorter lst in
            List.iter (fun (x, y, _) -> super#add_to_stack (x, y) to_do_later) lst;
            self#select_existant_and_update 

        val mutable has_potential = false

        method virtual private should_continue_this_path : Tree.edge option -> bool -> bool

        method next_edge = 
            has_potential <- false;
            let edge, do_check, depth = 
                if Stack.is_empty to_calculate_for_compare then 
                    if Stack.is_empty to_compare then
                        if Stack.is_empty to_do_later then None, false, 0
                        else self#select_existant_and_update, false, 0
                    else self#sort_and_choose_join_edge, false, 0
                else begin 
                    has_potential <- true;
                    let edge, depth = Stack.pop to_calculate_for_compare in
                    Some edge, true, depth
                end
            in
            if depth > max_distance then begin
                self#next_edge
            end else begin
                if do_check then begin
                    if self#should_continue_this_path edge has_potential then 
                        begin
                        currently_doing <- edge;
                        current_depth <- depth;
                    end else begin
                        currently_doing <- None;
                        current_depth <- 0;
                    end
                end else ();
                edge
            end

        method break_distance (x : float) =
            begin match currently_doing with
            | Some edge -> Stack.push (edge, current_depth, x) to_compare;
            | None -> ()
            end;
            currently_doing <- None

    end

    class distance_dfs side max_distance ptree : [Node.n, Edge.e] edges_manager = 
        object (self)
        inherit dfs_distance_based side max_distance ptree

        method private should_continue_this_path (y : Tree.edge option) x = 
            if not x then x
            else 
                match y with
                | None -> false
                | _ -> true

        method clone = 
            ({<
            currently_doing = currently_doing;
            to_compare = Stack.copy to_compare;
            to_calculate_for_compare = Stack.copy to_calculate_for_compare;
            to_do_later = Stack.copy to_do_later;
        >} :> (Node.n, Edge.e) edges_manager)

    end

    let order_edge (a, b) = min a b, max a b

    class union_dfs side max_distance edges ptree : [Node.n, Edge.e] edges_manager =
        object (self)
        inherit dfs_distance_based side max_distance ptree as super

        val mutable remove_first = 0

        val mutable do_not_cross = 
            List.fold_left 
            (fun acc (Tree.Edge ((a, b) as edge)) -> 
                acc 
                --> All_sets.FullTuples.add edge 
                --> All_sets.FullTuples.add (b, a))
            All_sets.FullTuples.empty
            edges

        val mutable extra_constraints = []

        method clone = ({<
            currently_doing = currently_doing;
            to_compare = Stack.copy to_compare;
            to_calculate_for_compare = Stack.copy to_calculate_for_compare;
            to_do_later = Stack.copy to_do_later;
        >} :> (Node.n, Edge.e) edges_manager)

        method update_break breakage =
            super#update_break breakage;
            self#clearup_remove_first;
            match get_side side breakage with
            | `Edge (o, s1, s2, _) ->
                    ()
                    (*
                    let hass1 = All_sets.FullTuples.mem o do_not_cross in
                    if hass1 then extra_constraints <- [s1; s2];
                    *)
            | `Single (_, _) -> ()

        method next_edge = 
            match super#next_edge with
            | None -> None
            | (Some (Tree.Edge ((u, v) as edge))) as res ->
                    if (All_sets.FullTuples.mem edge do_not_cross) then 
                        self#next_edge
                    else res

        method private clearup_remove_first =
            extra_constraints <- []

        method private should_continue_this_path edge x =
            x &&
                (match edge with
                | None -> false
                | Some (Tree.Edge ((a, _) as edge)) ->
                        let edge = order_edge edge in
                    if not (All_sets.FullTuples.mem edge do_not_cross) then
                        let clade =
                            match current_clade with
                            | Some clade -> clade
                            | None -> failwith "No root?"
                        in
                        let ca = Ptree.get_node_data a current_broken_ptree in
                        current_delta >= (Node.union_distance clade ca)
                    else false)
    end 

    module UnionTree = UnionTree.Make (Node) (Edge)

    let dsp = 0.7 (* default saturation parameter *)

    (* A class to choose the next join edge based on the cost of joining with the
    * union of the subtree rooted by each node *)
    class partitioned_dfs side max_distance edges ptree = object 
        inherit (dfs_distance_based side max_distance ptree)  as super

        val mutable extra_constraints = [] 

        val do_not_cross = 
            List.fold_left 
            (fun acc (Tree.Edge (a, _)) -> All_sets.Integers.add a acc) 
            All_sets.Integers.empty
            edges

        method update_break breakage =
            extra_constraints <- [];
            let () = 
                match get_side side breakage with
                | `Edge (o, s1, s2, _) ->
                        let hass1 = All_sets.Integers.mem o do_not_cross in
                        if hass1 then extra_constraints <- [s1; s2];
                | `Single (_, _) -> ()
            in
            super#update_break breakage

        method clone = ({<
            to_compare = Stack.copy to_compare;
            to_calculate_for_compare = Stack.copy to_calculate_for_compare;
            to_do_later = Stack.copy to_do_later;
        >} :> (Node.n, Edge.e) edges_manager)

        method private should_continue_this_path edge x =
            x && 
            (match edge with
            | None -> false
            | Some (Tree.Edge (a, b)) ->
                    let tree = current_broken_ptree.Ptree.tree in
                    let base = 
                        if a = Tree.get_parent b tree then a
                        else b
                    in
                    not ((All_sets.Integers.mem base do_not_cross) &&
                    (List.exists (fun x -> x = a || x = b) extra_constraints)))
    end

    class union_dfs_print_limits side max_distance ptree : 
        [Node.n, Edge.e] edges_manager = object
        inherit dfs_distance_based side max_distance ptree

        method clone = ({<
            to_compare = Stack.copy to_compare;
            to_calculate_for_compare = Stack.copy to_calculate_for_compare;
            to_do_later = Stack.copy to_do_later;
        >} :> (Node.n, Edge.e) edges_manager)

        method private should_continue_this_path edge x =
            x &&
                (match edge with
                | None -> false
                | Some (Tree.Edge (a, b)) ->
                        let tree = current_broken_ptree.Ptree.tree in
                        if a = Tree.get_parent b tree then
                            if Tree.is_leaf b tree then false
                            else 
                                let clade =
                                    match current_clade with
                                    | Some clade -> clade
                                    | None -> failwith "No root?"
                                in
                                let ca = 
                                    Ptree.get_node_data a current_broken_ptree 
                                in
                                let union_delta = 
                                    Node.union_distance clade ca 
                                in
                                current_delta >= union_delta 
                        else true)
    end

    class lazy_union_calculations (ptree : (Node.n, Edge.e) Ptree.p_tree) = 
            let empty = 
                { (Ptree.empty ptree.Ptree.data) with 
                    Ptree.tree = ptree.Ptree.tree } in
        object (self)

        val mutable tree : 
            ([`Known of int | `Beyond | `Unknown | `FirstTime], unit) Ptree.p_tree = 
                empty

        val mutable barriers_table = Hashtbl.create 157

        val mutable broken_edge = None

        method is_not_neighbor (Tree.Edge (a, b)) (Tree.Edge (c, d)) =
            (a <> c) && (a <> d) && (b <> c) && (b <> d)

        method store_current_edge =
            match broken_edge with
            | None -> ()
            | Some (Tree.Edge (a, b)) ->
                    Hashtbl.replace barriers_table a tree;
                    Hashtbl.replace barriers_table b tree

        method get_a_neighbor_calculation (Tree.Edge (a, b)) =
            let res = 
                try Hashtbl.find barriers_table a with
                | Not_found ->
                        try Hashtbl.find barriers_table b with
                        | Not_found -> 
                                empty
            in
            tree <- res

        method update_break (edge : Tree.edge) =
            let _ =
                match broken_edge with
                | Some broken_edge ->
                        self#store_current_edge;
                        if self#is_not_neighbor edge broken_edge then
                            self#get_a_neighbor_calculation edge;
                | None -> ()
            in
            broken_edge <- Some edge
                    

        method get_vertex_distance_to_barrier b =
            try Ptree.get_node_data b tree with
            | _ -> `FirstTime

        method get_parent node =
            if Tree.is_handle node tree.Ptree.tree then raise Not_found 
            else
                match Tree.get_node node tree.Ptree.tree with
                | Tree.Interior (_, p, _, _) 
                | Tree.Leaf (_, p) -> p
                | Tree.Single _ -> failwith "Single?"

        method set_current_vertex v state =
            tree <- Ptree.add_node_data v state tree

        method set_children state b b_parent =
            try
                let node = Ptree.get_node b tree in
                let (a, b) = Tree.other_two_nbrs b_parent node in
                self#set_current_vertex a state;
                self#set_current_vertex b state;
            with
            | _ -> ()

        method check_if_this_is_barrier 
            (broken_tree : (Node.n, Edge.e) Ptree.p_tree) b b_parent
            clade delta cb =
                let union_delta = Node.union_distance clade cb in
                if delta < union_delta then 
                    let parent = 
                            match 
                                try Some (Ptree.get_node_data b_parent broken_tree)
                                with
                                | Not_found -> None
                            with
                            | Some cb ->
                                    let parent = self#get_parent b_parent in
                                    self#check_if_this_is_barrier broken_tree b_parent 
                                    parent clade delta cb
                            | None -> true
                    in
                    let _ = 
                        if parent then 
                            let _ = self#set_current_vertex b (`Known 0) in
                            self#set_children `Beyond b b_parent
                        else 
                            let _ = self#set_current_vertex b (`Beyond) in
                            self#set_children `Beyond b b_parent
                    in
                    false
                else 
                    let _ = self#set_children `Unknown b b_parent in
                    true

        method check_depth_to_make_barrier_check tree b a clade delta cb =
            if (Node.num_otus (Some (Ptree.get_parent (Node.taxon_code cb) tree)) cb) < 3 then 
                self#check_if_this_is_barrier tree b a clade delta cb
            else 
                true

        method should_continue_this_path broken_tree edge clade delta =
            try
                match edge with
                | None -> false
                | Some (Tree.Edge (a, b)) ->
                        if a = Tree.get_parent b tree.Ptree.tree then 
                            if Tree.is_leaf b tree.Ptree.tree then false
                            else 
                                match self#get_vertex_distance_to_barrier b with
                                | `FirstTime ->
                                        let cb = Ptree.get_node_data b broken_tree in
                                        self#check_if_this_is_barrier broken_tree b a
                                        clade delta cb
                                | `Unknown
                                | `Known 0 ->
                                        let cb = Ptree.get_node_data b broken_tree in
                                        self#check_if_this_is_barrier broken_tree b a
                                        clade delta cb
                                | `Beyond -> false
                                | `Known _ -> true
                        else true
            with
            | Not_found -> true

    end

    class lazy_union_dfs side max_distance ptree : [Node.n, Edge.e] edges_manager = 
        object 
        inherit dfs_distance_based side max_distance ptree as super

        val mutable union = new lazy_union_calculations ptree

        method clone = ({<
            to_compare = Stack.copy to_compare;
            to_calculate_for_compare = Stack.copy to_calculate_for_compare;
            to_do_later = Stack.copy to_do_later;
        >} :> (Node.n, Edge.e) edges_manager)

        method private should_continue_this_path edge x =
            let clade = 
                match current_clade with
                | Some c -> c
                | None -> failwith "No root?"
            in
            x && union#should_continue_this_path current_broken_ptree edge clade current_delta

        method next_edge =
            let edge = super#next_edge in
            match edge with
            | Some e ->
                    union#update_break e;
                    edge
            | None -> None

        method update_join tree r =
            union <- new lazy_union_calculations tree;
            super#update_join tree r

    end 

    (* A class to choose where to root next using bfs, and with maximum depth
    * max_depth *)
    class bfs_based side limits (max_depth : int) (ptree : (Node.n, Edge.e) Ptree.p_tree) : 
        [Node.n, Edge.e] edges_manager =
        object (self)

        val to_do = Queue.create ()
        val mutable current_tree = ptree
        val mutable exclude_edges = limits

        val mutable remove_first = 0

        method update_reroot (_  : (Node.n, Edge.e) Ptree.breakage) = ()
        method exclude edges = 
            exclude_edges <- edges

        method clone = ({<
            to_do = Queue.copy to_do;
        >} :> (Node.n, Edge.e) edges_manager)

        method break_distance (_ : float) = ()

        method private add_children_to_queue child parent depth =
            if depth >= max_depth then ()
            else 
                try
                    match Ptree.get_node child current_tree with
                    | Tree.Single _ | Tree.Leaf _ -> ()
                    | node ->
                            let exists x = 
                                List.exists (undirected_edge (x,
                                child)) exclude_edges 
                            in
                            let (a, b) = Tree.other_two_nbrs parent node in
                            if exists a then ()
                            else Queue.push (a, child, depth + 1) to_do;
                            if exists b then ()
                            else Queue.push (b, child, depth + 1) to_do;
                with
                | err -> 
                        Printf.fprintf stderr "There is an error?\n%!";
                        raise err

        method private clear_queues = Queue.clear to_do

        method private get_children_edges node ptree =
            match Ptree.get_node node ptree with
            | Tree.Interior (_, _, a, b) -> 
                    remove_first <- remove_first + 2;
                    [Tree.Edge (node, a); Tree.Edge (node, b)]
            | _ -> []

        method private clearup_remove_first =
            if remove_first > 0 then begin
                exclude_edges <-
                    (match exclude_edges with
                    | _ :: _ :: t -> t
                    | _ -> []);
                remove_first <- remove_first - 2;
                self#clearup_remove_first 
            end else ()

        method update_break breakage =
                let ptree = breakage.Ptree.ptree in
                self#clearup_remove_first;
                current_tree <- ptree;
                self#clear_queues;
                match get_side side breakage with
                | `Edge (o, s1, s2, _) ->
                        let exists x =
                            List.exists (undirected_edge (o, x)) exclude_edges
                        in
                        let hass1 = exists s1
                        and hass2 = exists s2 in
                        if hass1 then 
                            exclude_edges <- 
                                (self#get_children_edges s1 ptree) 
                                @ exclude_edges;
                        if hass2 then 
                            exclude_edges <- 
                                (self#get_children_edges s2 ptree) 
                                @ exclude_edges;
                        Queue.push (s1, s2, 0) to_do;
                        self#add_children_to_queue s2 s1 0;
                | `Single (_, _) -> ()

        method update_join (tree : (Node.n, Edge.e) Ptree.p_tree) (_ : Tree.join_delta) = 
            current_tree <- tree

        method next_edge =
            if Queue.is_empty to_do then None
            else begin
                let (child, par, depth) = Queue.pop to_do in
                self#add_children_to_queue child par depth;
                Some (Tree.Edge (child, par))
        end

    end

    class virtual all_edges side initial_edges (ptree : (Node.n, Edge.e) Ptree.p_tree) =
        object (self)

        val mutable edges = initial_edges
        val mutable next_edge = -1
        val mutable end_of_check = (Array.length initial_edges) - 1
        val mutable break_deltas = None
        val mutable this_is_the_end = false

        method exclude (_ : Tree.edge list) = ()

        method clone = 
            ({<
            edges = Array.copy edges;
            break_deltas = break_deltas;
            next_edge = next_edge;
            end_of_check = end_of_check;
            this_is_the_end = this_is_the_end;
        >} :> (Node.n, Edge.e) edges_manager)

        method private check_consistency tree edges =
            let test = ref false in
            Array.iter (fun ((Tree.Edge (a, b)) as x) -> 
                try 
                    let _ = Tree.normalize_edge x tree.Ptree.tree in
                    ()
                with
                | _ -> 
                        test := true;
                        Status.user_message Status.Error ("Found an inconsistent
                edge "  ^ string_of_int a ^ ", " ^ string_of_int b)) edges;
            if !test then failwith "Failed"
            else ()

        method virtual private update_with_deltas : (Node.n, Edge.e)
        Ptree.p_tree ->
            Tree.join_delta -> unit

        method break_distance (_ : float) = ()

        method next_edge =
            if 0 = Array.length edges then begin
                this_is_the_end <- true;
                None;
            end else if this_is_the_end then None
            else 
                let dges = Array.length edges in
                let () = next_edge <- (next_edge + 1) mod dges in
                if next_edge = end_of_check then 
                    let () = this_is_the_end <- true in
                    Some (edges.(next_edge))
                else Some (edges.(next_edge))

        method update_break breakage = 
               break_deltas <- Some breakage.Ptree.tree_delta

        method update_join (tree : (Node.n, Edge.e) Ptree.p_tree) (join_delta : Tree.join_delta) = 
            self#update_with_deltas tree join_delta;
            break_deltas <- None

    end

    let calculate_edges ptree = 
        let res = Tree.get_edges_tree ptree.Ptree.tree in
        let res = Array.of_list res in
        rand res;
        res

    class randomized_edges (ptree : (Node.n, Edge.e) Ptree.p_tree) :
        [Node.n, Edge.e] edges_manager =
        object (self)
            inherit all_edges `Both (calculate_edges ptree) ptree

        method update_reroot (_  : (Node.n, Edge.e) Ptree.breakage) = ()
        method private update_with_deltas tree tdelta =
            let removed = 
                match break_deltas with
                | Some (a, b) -> (Tree.break_to_edge_delta (a, b))
                | None -> failwith "No deltas?"
            and gen_created = Tree.join_to_edge_delta tdelta in
            let created = ref (gen_created.Tree.added @ removed.Tree.added) in
            let is_the_same (Tree.Edge (a, b)) (Tree.Edge (c, d)) = 
                ((a = c) && (b = d)) || ((a = d) && (b = c))
            in
            let removed = removed.Tree.removed @ gen_created.Tree.removed in
            let rec filter x = 
                if List.exists (is_the_same x) removed then 
                    match !created with
                    | h :: t ->
                            created := t;
                            if 
                                List.exists (is_the_same h) gen_created.Tree.removed
                            then filter x
                            else h
                    | [] -> failwith "Different number?"
                else x
            in
            edges <- Array.map filter edges;
            (*assert (0 = List.length !created);*)
            end_of_check <- 
                if next_edge < 0 then (Array.length edges) - 1 else next_edge;
            this_is_the_end <- false

    end

    module EdgeComparator = struct
        type t = Tree.edge
        let compare x y =
            match x, y with
            | Tree.Edge (x1, x2), Tree.Edge (y1, y2) ->
                    let x1 = min x1 x2 
                    and x2 = max x1 x2
                    and y1 = min y1 y2
                    and y2 = max y1 y2 in
                    match x1 - y1 with
                    | 0 -> x2 - y2
                    | x -> x
    end

    module EdgeSet = Set.Make (EdgeComparator)

    class only_once_break (ptree : (Node.n, Edge.e) Ptree.p_tree) :
        [Node.n, Edge.e] edges_manager = object (self)

            val mutable edges = 
                Tree.EdgeSet.fold (fun x acc ->
                    EdgeSet.add x acc) 
                ptree.Ptree.tree.Tree.d_edges
                EdgeSet.empty;

        method update_reroot (_  : (Node.n, Edge.e) Ptree.breakage) = ()
            val mutable break_delta = None
            val mutable current_tree = ptree.Ptree.tree

            method break_distance _ = ()
            method next_edge =
                try
                    let next = EdgeSet.choose edges in 
                    edges <- EdgeSet.remove next edges;
                    let next = Tree.normalize_edge next current_tree in
                    Some next
                with
                | _ -> None

            method update_break breakage = 
                break_delta <- Some breakage.Ptree.tree_delta

            method update_join nt (a, b, _) =
                current_tree <- nt.Ptree.tree;
                let add_update_side side =
                    match side with
                    | `Single a -> ()
                    | `Edge (a, b, c, _) ->
                            edges <- 
                                (edges 
                                --> EdgeSet.add (Tree.Edge (a, b))
                                --> EdgeSet.add (Tree.Edge (a, c)))
                in
                let remove_update_side side =
                    match side with
                    | `Single _ -> ()
                    | `Edge (a, b, c, _) ->
                            edges <-
                                (edges 
                                --> EdgeSet.add (Tree.Edge (b, c))
                                --> EdgeSet.remove (Tree.Edge (a, b))
                                --> EdgeSet.remove (Tree.Edge (a, c)))
                in
                let _ = 
                    match break_delta with
                    | Some (a, b) -> 
                            remove_update_side a;
                            remove_update_side b;
                    | None -> ()
                in
                add_update_side a;
                add_update_side b

            method clone = ({<>} :> (Node.n, Edge.e) edges_manager)
            method exclude _ = ()

    end

    let make_edge_list dont_cut ptree =
        let cnt = 
            let cnt = Hashtbl.create 1667 in
            List.iter 
            (fun ((Tree.Edge ((a, b) as e))) ->
                Hashtbl.add cnt e 1;
                Hashtbl.add cnt (b, a) 1)
            dont_cut;
            cnt
            (*
            let res = Hashtbl.create 1667 in
            let cnt = Hashtbl.create 1667 in
            let add x = 
                if Hashtbl.mem cnt x then
                    let c = Hashtbl.find cnt x in
                    if c = 2 then Hashtbl.add res x 3
                    else Hashtbl.replace cnt x (1 + (Hashtbl.find cnt x))
                else Hashtbl.add cnt x 1
            in
            List.iter (fun (Tree.Edge (a, b)) -> add (a, b); add (b, a)) dont_cut;
            cnt
            *)
        in
        let tabu_distance2 ptree (Tree.Edge (x, y)) =
            let xn = Ptree.get_node_data x ptree
            and yn = Ptree.get_node_data y ptree in
            Node.edge_distance xn yn
        in
        let filter (a, (Tree.Edge x)) = not (Hashtbl.mem cnt x) in
        let handle = All_sets.Integers.choose (Ptree.get_handles ptree) in
        let edge_lst = Ptree.get_pre_order_edges handle ptree in
        let id_cost = List.map (fun x -> tabu_distance2 ptree x, x) edge_lst in
        let id_cost = List.filter filter id_cost in
        let id_cost = List.sort (fun (a, _) (b, _) -> compare b a) id_cost in
        id_cost

    let create_initial_edges dont_cut ptree = 
        let edges = make_edge_list dont_cut ptree in
        let initial_edges = Array.of_list edges in
        Array.map (fun (_, x) -> x) initial_edges

    class sort_edges_by_cost early_stop dont_cut ptree =
        object (self)
            inherit all_edges `Both (create_initial_edges dont_cut ptree) ptree
            as super

        val mutable es = early_stop
        method update_reroot (_  : (Node.n, Edge.e) Ptree.breakage) = ()

        method next_edge =
            if es && (next_edge > (Array.length edges) / 2) then None
            else super#next_edge

        method private update_with_deltas tree tdelta =
            es <- false;
            let removed = 
                match break_deltas with
                | Some (a, b) -> (Tree.break_to_edge_delta (a, b))
                | None -> failwith "No deltas?"
            and gen_created = Tree.join_to_edge_delta tdelta in
            let created = ref (gen_created.Tree.added @ removed.Tree.added) in
            let is_the_same (Tree.Edge (a, b)) (Tree.Edge (c, d)) = 
                ((a = c) && (b = d)) || ((a = d) && (b = c))
            in
            let removed = removed.Tree.removed @ gen_created.Tree.removed in
            let rec filter x = 
                if List.exists (is_the_same x) removed then 
                    match !created with
                    | h :: t ->
                            created := t;
                            if
                                List.exists (is_the_same h) 
                                gen_created.Tree.removed
                            then 
                                filter x
                            else h
                    | [] -> failwith "Different number?"
                else x
            in
            edges <- Array.map filter edges;
            (*assert (0 = List.length !created);*)
            end_of_check <- 
                if next_edge < 0 then (Array.length edges) - 1 else next_edge;
            this_is_the_end <- false

    end

    type emc  =((Node.n, Edge.e) Ptree.p_tree -> (Node.n, Edge.e) edges_manager)
    type semc = [`Left | `Right] -> emc

    class side_manager side ptree (join : semc) (reroot : semc) = object
        val join = join side ptree
        val reroot = reroot side ptree 

        method break_distance x =
            join#break_distance x

        method join_edge = join#next_edge

        method reroot_edge = reroot#next_edge

        method update_break breakage =
            join#update_break breakage;
            reroot#update_break breakage

        method update_join tree jd =
            join#update_join tree jd;
            reroot#update_join tree jd

        method update_reroot breakage = 
            join#update_reroot breakage;
            reroot#update_reroot breakage

        method clone = {< 
            join = join#clone;
            reroot = reroot#clone;
        >}

    end

    (* A submodule to centralize all the union related operations to avoid further
    * confussion of this tabus module *)
    module Unions = struct
        (* Let tree be a tree, with vertex vertex, if the closest ancestor
        * present in the union of true_vertices and false_vertices is in 
        * true_vertices, then the edge is added to queue, otherwise it is not.
        * The cur value is used to carry on the value of the latest decision. *)
        let rec get_tf_edges tree vertex queue true_vertices false_vertices cur = 
            match Tree.get_node vertex tree with
            | Tree.Interior (_, _, c, d) ->
                    let cur =
                        (cur || (All_sets.Integers.mem vertex true_vertices))
                        && (not (All_sets.Integers.mem vertex false_vertices)) 
                    in
                    if cur then begin
                        Queue.push (Tree.Edge (vertex, c)) queue;
                        Queue.push (Tree.Edge (vertex, d)) queue;
                    end;
                    get_tf_edges tree c queue true_vertices false_vertices cur;
                    get_tf_edges tree d queue true_vertices false_vertices cur;
            | _ -> ()

        (* This function adds all the edges according to the definition of the
        * get_tf_edges function, but performs it on the edge of the handle *)
        let get_edges add_handle queue tree handle true_vertices false_vertices =
            match Tree.get_node handle tree with
            | Tree.Interior (a, b, _, _)
            | Tree.Leaf (a, b) ->
                    get_tf_edges tree a queue true_vertices false_vertices false;
                    get_tf_edges tree b queue true_vertices false_vertices false;
                    if add_handle then Queue.push (Tree.Edge (a, b)) queue;
            | _ -> ()

        (* Find possible unions between the union with code [union_code] and the other
        * unions in [all_unions], when their mutual distance is less than [delta], 
        * using the [matrix] of distances. The function returns a tuple of sets 
        * [(a, b)], where [a] are the unions that have distance to [union_code] less
        * than [delta], and [b] are the unions that doesn't show this property. *)
        let find_possible_unions matrix union_code all_unions delta = 
            let folder y _ (t, f) =
                let dist = 
                    try SparceMatrix.find union_code y matrix with
                    | err ->
                            Status.user_message Status.Error
                            ("Failed to find " ^ string_of_int union_code);
                            raise err
                in
                if dist <= delta then 
                    All_sets.Integers.add y t, f
                else t, All_sets.Integers.add y f
            in
            let empty = All_sets.Integers.empty in
            All_sets.IntegerMap.fold folder all_unions (empty, empty)
            
        (* A class to join vertices using information from the unions of the
        * vertices *)
        class wagner_joiner sat_par ptree handle : wem = 
            let initial_union_tree, initial_clusters = 
                UnionTree.create_initial_tree (Some handle) sat_par ptree 
            in
            let initial_matrix = 
                UnionTree.Clusters.initial_matrix initial_clusters 
            in
            object (self)
            val queue = Queue.create ()
            val mutable current_tree_cost = Ptree.get_cost `Adjusted ptree
            val mutable matrix = initial_matrix
            val mutable union_tree = initial_union_tree
            val mutable clusters = initial_clusters
            val mutable current_delta = max_float
            val mutable current_clade = None
            val unions = Queue.create ()
            val mutable cluster_set = All_sets.IntegerMap.fold 
                (fun x _ acc -> All_sets.Integers.add x acc) 
                initial_clusters All_sets.Integers.empty

            method break_distance _ = ()

            method next_edge =
                try Some (Queue.pop queue) with
                | Queue.Empty -> 
                        try 
                            let un = Queue.pop unions in
                            let clust = 
                                let r = All_sets.IntegerMap.find un clusters
                                in
                                All_sets.IntegerMap.add un r 
                                (All_sets.IntegerMap.empty)
                            in
                            let t, f =
                                let vertex = 
                                    match current_clade with
                                    | None -> failwith "No node?"
                                    | Some v -> v
                                in
                                find_possible_unions matrix vertex clust
                                current_delta
                            in
                            if All_sets.Integers.is_empty t then 
                                self#next_edge
                            else begin
                                get_edges true queue union_tree.Ptree.tree handle 
                                t (All_sets.Integers.remove un cluster_set);
                                self#next_edge
                            end
                        with
                        | Queue.Empty -> None



            method private create_update_break_vertices 
                (delta1, delta2) tree acc =
                    let update_side delta acc = 
                        match delta with
                        | (`Single _) -> acc
                        | (`Edge (h, x, v, _)) ->
                                let to_handle = 
                                    Tree.get_vertices_to_handle v tree 
                                in
                                let set = 
                                    List.fold_right (fun x acc -> 
                                        All_sets.Integers.add x acc)
                                    to_handle
                                    (All_sets.Integers.add h acc)
                                in
                                All_sets.Integers.add x set
                    in
                    acc
                    --> update_side delta1 
                    --> update_side delta2

            method new_delta f = current_delta <- f

            method next_clade clade = 
                (* We have to add a new cluster, that only contains this clade.
                * It is equivalent to breaking in a leaf. *)
(*                let t,f =*)
                    (* We add the vertex to the clusters, then add it to the
                    * matrix, and we only give the element in that list that has
                    * lowest cost as a suitable cluster *)
                    let vertex = Node.taxon_code clade in
                    current_clade <- Some vertex;
                    let uclade = 
                        clade 
                        --> Node.Union.leaf None None 
                        --> UnionTree.leaf [] vertex 
                    in
                    All_sets.IntegerMap.iter 
                        (fun x _ -> Queue.add x unions) clusters;
                    clusters <- All_sets.IntegerMap.add vertex uclade clusters;
                    matrix <- 
                        UnionTree.Clusters.update_matrix matrix clusters 
                        (All_sets.Integers.singleton vertex)
                        All_sets.Integers.empty;
                        (*
                    let delta_cost = SparceMatrix.minimum vertex matrix in
                    find_possible_unions matrix vertex clusters (delta_cost +.
                    1.)
                        *)
(*                in*)
                Queue.clear queue
(*                get_edges true queue union_tree.Ptree.tree handle t f;*)

            method update_join tree (delt1, delt2, ids) =
                let new_union_tree = 
                    { union_tree with Ptree.tree = tree.Ptree.tree } 
                and my_update_down = 
                    All_sets.Integers.empty
                    --> self#create_update_break_vertices 
                        (delt1, delt2) union_tree.Ptree.tree
                    -->
                        (fun acc ->
                            List.fold_left (fun acc x ->
                                All_sets.Integers.add x acc)
                            acc ids)
                in
                let nt, cls, added, removed = 
                    UnionTree.update (Some handle) sat_par tree new_union_tree 
                    clusters my_update_down
                in
                let mtx = 
                    UnionTree.Clusters.update_matrix matrix cls added removed
                in
                current_delta <- max_float;
                matrix <- mtx;
                union_tree <- nt;
                clusters <- cls;
                cluster_set <- 
                    All_sets.Integers.union (All_sets.Integers.diff cluster_set
                    removed) added;
                current_tree_cost <- Ptree.get_cost `Adjusted tree;
                Queue.clear queue

            method clone = ({< >} :> wem)

            method exclude _ = ()

        end

        (* A class to join vertices using information from the unions of the
        * vertices *)
        class joiner side sat_par threshold ptree : [a, b] edges_manager = 
            let initial_union_tree, initial_clusters = 
                UnionTree.create_initial_tree None sat_par ptree 
            in
            let initial_matrix = 
                UnionTree.Clusters.initial_matrix initial_clusters 
            in
            let get_cost_without_encoding ptree =
                let adder _ b (its_first, acc) =
                    if its_first then 
                        match b.Ptree.root_median with
                        | Some (_, b) -> false, acc +. Node.total_cost None None b
                        | _ -> failwith "No root?"
                    else false, b.Ptree.component_cost +. acc
                in
                let _, res =
                    All_sets.IntegerMap.fold adder ptree.Ptree.component_root
                    (true, 0.)
                in
                res
            in
            object (self)
            val queue = Queue.create ()
            val mutable current_tree_cost = Ptree.get_cost `Adjusted ptree
            val mutable matrix = initial_matrix
            val mutable union_tree = initial_union_tree
            val mutable clusters = initial_clusters
            val mutable to_update_down = None 

            method update_reroot breakage = self#update_break breakage

            method break_distance _ = ()

            method next_edge =
                try Some (Queue.pop queue) with
                | Queue.Empty -> None

            method private create_update_break_vertices 
                (delta1, delta2) tree acc =
                    let update_side delta acc = 
                        match delta with
                        | (`Single _) -> acc
                        | (`Edge (h, x, v, _)) ->
                                let to_handle = 
                                    Tree.get_vertices_to_handle v tree 
                                in
                                let set = 
                                    List.fold_right (fun x acc -> 
                                        All_sets.Integers.add x acc)
                                    to_handle
                                    (All_sets.Integers.add h acc)
                                in
                                All_sets.Integers.add x set
                    in
                    acc
                    --> update_side delta1 
                    --> update_side delta2

            method update_break breakage =
                let ptree = breakage.Ptree.ptree 
                and my_delta = breakage.Ptree.tree_delta
                and handle, node_data = 
                    let clade_info =
                        match side with
                        | `Left -> breakage.Ptree.left
                        | `Right -> breakage.Ptree.right
                    in
                    clade_info.Ptree.clade_id, clade_info.Ptree.clade_node
                in
                let t,f =
                    let other_side = 
                        match side with
                        | `Left -> `Right
                        | `Right -> `Left 
                    in
                    let vertex =
                        match get_side other_side breakage with
                        | `Single (a, _) 
                        | `Edge (a, _, _, _) -> a
                    in
                    let delta_cost = 
                        current_tree_cost -. (get_cost_without_encoding ptree)
                    and vertex = 
                        UnionTree.get_union_vertex vertex union_tree clusters 
                    in
                    find_possible_unions matrix vertex clusters delta_cost
                in
                to_update_down <-
                    Some 
                    (self#create_update_break_vertices 
                    my_delta ptree.Ptree.tree All_sets.Integers.empty);
                Queue.clear queue;
                get_edges false queue ptree.Ptree.tree handle t f

            method update_join tree (delt1, delt2, ids) =
                let new_union_tree = 
                    { union_tree with Ptree.tree = tree.Ptree.tree } 
                and my_update_down = 
                    match to_update_down with
                    | Some v -> 
                            v 
                            --> self#create_update_break_vertices 
                                (delt1, delt2) union_tree.Ptree.tree
                            -->
                                (fun acc ->
                                    List.fold_left (fun acc x ->
                                        All_sets.Integers.add x acc)
                                    acc ids)
                    | None -> failwith "Nothing broken?"
                in
                let nt, cls, added, removed = 
                    UnionTree.update None sat_par tree new_union_tree 
                    clusters my_update_down
                in
                let mtx = 
                    UnionTree.Clusters.update_matrix matrix cls added removed
                in
                matrix <- mtx;
                union_tree <- nt;
                clusters <- cls;
                current_tree_cost <- Ptree.get_cost `Adjusted tree;
                to_update_down <- None;
                Queue.clear queue

            method clone = ({<>} :> (a, b) edges_manager)

            method exclude _ = ()

        end
    end

(* The functions to feed the standard tabu with *)

(* This tabu is used by andres for testing purposes *)
    let wagner_tabu ptree handle = new Unions.wagner_joiner dsp ptree handle

    let wagner_constraint constrain ptree handle = 
        new constrained_dfs_wagner constrain max_int ptree handle

    let distance_dfs_wagner ptree handle = 
        new distance_dfs_wagner max_int ptree handle

(* The official wagner tabu *)
    let wagner_union ptree handle = new union_dfs_wagner max_int ptree handle

(* These two union_join are currently being used by andres in his tests. They
* will be removed when the current union_dfs is officially deprecated *)
    let union_join max_distance side ptree = 
        if max_distance = max_int then new Unions.joiner side dsp 0.7 ptree
        else new union_dfs side max_distance [] ptree

    let distance_join max_distance side ptree = 
        new DFS.generic_distance side max_distance ptree

    let interior_vertex_sets_partition sets leaves ptree =
        (* Produce a list holding all the edges such that:
            *
            * If internal, the four partitions that the incident vertices
            * generate exist in sets
            *
            * If one is a leaf, then the three partitions implied by the other
            * incident vertex 
            *
            * The resultsing edges should not be violated during the search.
            * *)
        let is_mem x = 
            (All_sets.IntSet.mem x sets) || (All_sets.IntSet.mem
            (All_sets.Integers.diff leaves x) sets) 
        in
        let add_to_list do_add item lst = 
            if do_add then (Tree.Edge item) :: lst else lst
        in
        let rec tree_traverser parent vertex acc =
            match Ptree.get_node vertex ptree with
            | Tree.Leaf _ -> 
                    true, acc, All_sets.Integers.singleton vertex
            | Tree.Interior (_, _, a, b) ->
                    let a_chld_belong_to_set, acc, a_group = 
                        tree_traverser vertex a acc in
                    let b_chld_belong_to_set, acc, b_group =
                        tree_traverser vertex b acc in
                    let ab_group = All_sets.Integers.union a_group b_group in
                    let par_group = All_sets.Integers.diff leaves ab_group in
                    let res =
                        let par_group = is_mem par_group in
                        acc 
                        --> add_to_list 
                            (a_chld_belong_to_set && is_mem b_group && par_group) 
                            (vertex, a)
                        --> add_to_list 
                            (b_chld_belong_to_set && is_mem a_group && par_group)
                            (vertex, b)
                    in
                    is_mem b_group && is_mem a_group, res, ab_group
            | Tree.Single _ -> assert false
        in
        let of_handle ptree handle acc =
            match Ptree.get_node handle ptree with
            | Tree.Leaf (mine, other) -> 
                    let _, res, _ = tree_traverser mine other acc in
                    res
            | Tree.Single _ -> acc
            | Tree.Interior (mine, a, b, c) -> 
                    let a_add, acc, a_group = tree_traverser mine a acc in
                    let b_add, acc, b_group = tree_traverser mine b acc in
                    let c_add, acc, c_group = tree_traverser mine c acc in
                    let b_mem = is_mem b_group 
                    and a_mem = is_mem a_group
                    and c_mem = is_mem c_group in
                    acc 
                    --> add_to_list (a_add && b_mem && c_mem) (mine, a)
                    --> add_to_list (b_add && a_mem && c_mem) (mine, b)
                    --> add_to_list (c_add && b_mem && a_mem) (mine, c)
        in
        All_sets.Integers.fold (of_handle ptree) (Ptree.get_handles ptree) []

    let generate_partition_edges sets ptree = 
        let tree_leaves = 
            let leaves = Tree.get_all_leaves ptree.Ptree.tree in
            List.fold_left (fun acc x -> All_sets.Integers.add x acc)
            All_sets.Integers.empty leaves 
        in
        let res = interior_vertex_sets_partition sets tree_leaves ptree in
        res, res

    let max_height ptree =
        let best = ref max_int in
        let e = 
            Tree.EdgeMap.fold (fun ((Tree.Edge (a, b)) as edge) data acc ->
                let data = Edge.to_node (-1) (a, b) data in
                if !best > Node.num_height None data then 
                    let () = best := Node.num_height None data in
                    Some edge
                else acc) ptree.Ptree.edge_data None
        in
        match e with
        | None -> assert false
        | Some x -> [x]

    let partitioned_join sets max_distance side ptree = 
        let edges =
            match sets with
            | `Sets sets ->
                    fst (generate_partition_edges sets ptree)
            | _ -> max_height ptree
        in
        new union_dfs side max_distance edges ptree 

    let reroot sets max_int side ptree = 
        let dont_cross =
            match sets with
            | None -> []
            | Some (`Sets x) -> fst (generate_partition_edges x ptree) 
            | Some (`Height x) -> max_height ptree
        in
        new bfs_based side dont_cross max_int ptree

    let random_break ptree = new randomized_edges ptree

    let sorted_break early_stop sets ptree = 
        let dont_cut =
            match sets with
            | None -> []
            | Some (`Sets x) ->
                    snd (generate_partition_edges x ptree)
            | Some (`Height x) -> max_height ptree
        in
        new sort_edges_by_cost early_stop dont_cut ptree

    let only_once_break ptree = new only_once_break ptree
    let simple_dfs_from_middle ptree = new simple_dfs ptree
    
    let simple_nm_none c t = new nm_simple_none c t
    let simple_nm_all c t = new nm_simple_all c t
    let complex_nm_delta c t = new nm_complex_delta c t
    let complex_nm_neighborhood c t = new nm_complex_neighborhood c t

    let add_to_banned_list neigh1 neigh2 b1 b2 map_of_banned tree =
        let add_to_map_banned a neigh1 b1 b2 map_of_banned =
            Tree.EdgeMap.add (Tree.Edge (a, neigh1)) 
            (Tree.EdgeSet.add (Tree.Edge (b1, b2)) 
            (try Tree.EdgeMap.find (Tree.Edge (a, neigh1)) map_of_banned with
             Not_found -> Tree.EdgeSet.empty))
            map_of_banned
        in
        let add_one_direction neigh1 neight2 map_of_banned =
            match Ptree.get_node neigh1 tree with
            | Tree.Single _ | Tree.Leaf _ -> map_of_banned
            | node ->
                    let a, b = Tree.other_two_nbrs neight2 node in
                    map_of_banned 
                    --> add_to_map_banned a neigh1 b1 b2 
                    --> add_to_map_banned b neigh1 b1 b2
        in
        map_of_banned
        --> add_one_direction neigh1 neigh2
        --> add_one_direction neigh2 neigh1

    let total_banned_size map =
        Tree.EdgeMap.fold (fun _ set acc ->
            acc + (Tree.EdgeSet.cardinal set)) map 0

    let add_to_banned_list map_of_banned broken left right tree =
        match broken with
        | Some (Tree.Edge (b1, b2)) -> 
                (*
                Status.user_message Status.Information ("The initial size is " ^
                string_of_int (total_banned_size map_of_banned));
                *)
                let process_side side map_of_banned =
                    match side with
                    | None -> map_of_banned
                    | Some (Tree.Edge (neigh1, neigh2)) -> 
                            add_to_banned_list neigh1 neigh2 b1 b2 map_of_banned
                            tree
                in
                let res = map_of_banned --> process_side left --> process_side
                right in
                (*
                Status.user_message Status.Information ("The new size is " ^
                string_of_int (total_banned_size res));
                *)
                res
        | None -> map_of_banned

    let find_banned_for_broke edge map =
        match edge with
        | Some ((Tree.Edge (a, b)) as e) -> 
                let baned = 
                    Tree.EdgeSet.union 
                    (try Tree.EdgeMap.find e map with Not_found ->
                        Tree.EdgeSet.empty)
                    (try Tree.EdgeMap.find (Tree.Edge (b, a)) map with Not_found ->
                        Tree.EdgeSet.empty)
                in
                (*
                Status.user_message Status.Information 
                ("The current baned has size " ^ string_of_int
                (Tree.EdgeSet.cardinal baned));
                *)
                baned
        | None -> Tree.EdgeSet.empty

    let is_banned did_reroot e set =
        if not did_reroot then
            match e with
            | Some ((Tree.Edge (a, b)) as e) ->
                    Tree.EdgeSet.mem e set || Tree.EdgeSet.mem (Tree.Edge (b, a)) set
            | None -> false
        else false

        (* see, ('a, 'b) Ptree.tabu_mgr *)
    class standard_tabu (ptree : (Node.n, Edge.e) Ptree.p_tree) (join : semc) 
    (reroot : semc) (break : emc) (iterate : nm) = object (self)
        val timer = Timer.start ()
        val mutable current_time = 0
        val left = new side_manager `Left ptree join reroot
        val right = new side_manager `Right ptree join reroot
        val break = break ptree
        val adjust_mgr = Some iterate
        val mutable last_side : [`Left | `Right] = `Left
        val banned_map = ref Tree.EdgeMap.empty
        val mutable current_banned = Tree.EdgeSet.empty
        val mutable last_edge_to_break = None

        val mutable left_edge = None
        val mutable right_edge = None

        val mutable did_reroot = false

        method clone = ({< 
            last_edge_to_break = last_edge_to_break;
            left = left#clone;
            right = right#clone;
            break = break#clone;
        >} :> (Node.n, Edge.e) Ptree.tabu_mgr)

        method break_edge = 
            current_time <- print_time timer current_time;
            let r = break#next_edge in
            last_edge_to_break <- r;
            r

        method break_distance x = 
            let mgr = 
                match last_side with
                | `Left -> left
                | `Right -> right
            in
            mgr#break_distance x

        method join_edge side =
            last_side <- side;
            let to_join =
                match side with
                | `Left -> left#join_edge
                | `Right -> right#join_edge
            in
            if is_banned did_reroot to_join current_banned then 
                self#join_edge side
            else to_join


        method reroot_edge (side : [`Left | `Right]) =
            did_reroot <- true;
            match side with
            | `Left ->
                    (match left#reroot_edge with
                    | None -> None
                    | new_root ->
                            if new_root = left_edge then 
                                self#reroot_edge `Left
                            else new_root)
            | `Right ->
                    match right#reroot_edge with
                    | None -> None
                    | new_root ->
                            if new_root = right_edge then 
                                self#reroot_edge `Right
                            else new_root

        method update_break breakage =
            let tree = breakage.Ptree.ptree in
            did_reroot <- false;
            let () =
                let my_side = function
                    | (`Edge (_, a, b, _)) -> Some (Tree.Edge (a, b))
                    | _ -> None
                in
                left_edge <- my_side (get_side `Left breakage);
                right_edge <- my_side (get_side `Right breakage)
            in
            banned_map :=
                add_to_banned_list !banned_map last_edge_to_break
                left_edge right_edge tree;
            current_banned <- 
                find_banned_for_broke last_edge_to_break
                !banned_map;
            left#update_break breakage;
            right#update_break breakage;
            break#update_break breakage

        method update_reroot breakage =
            left#update_reroot breakage;
            right#update_reroot breakage;
            break#update_reroot breakage

        method update_join tree jd =
            banned_map := !banned_map;
            left#update_join tree jd;
            right#update_join tree jd;
            break#update_join tree jd
        
        method get_node_manager = adjust_mgr

        method break_edges : Tree.edge list = []

        method features lst : (string * string) list = lst

    end

end
