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

let () = SadmanOutput.register "Sampler" "$Revision: 1616 $"

(* The sampler module is divided in two submodules, one containing te
* application related samplers (App), and one with research intended samplers
* (Res). The former are visible for external users, the later only for people
* interested on exploring new heuristics *)

type ft_queue = {
        queue : (Tree.u_tree * float * Ptree.clade_cost) list Queue.t;
        stack : (Tree.u_tree * float * Ptree.clade_cost) list Stack.t;
    }

let create () = {
    queue = Queue.create ();
    stack = Stack.create ();
}

class type ['a, 'b] search_manager_sampler = object
    method init : 
        (('a, 'b) Ptree.p_tree * float * Ptree.clade_cost) list -> unit
    method clone : ('a, 'b) search_manager_sampler
    (* process a b c d e f joins the tree c in positions a and b, with delta
    * cost d, which heuristic join cost e while the actual join cost (if appears
    * better), f. *)
    method process : 
            Ptree.incremental list ->
                Tree.join_jxn -> Tree.join_jxn -> 'a ->
                ('a, 'b) Ptree.p_tree -> ('a, 'b) Ptree.p_tree option -> 
                    float -> float -> float option -> unit
    method any_trees : bool -> unit
    method next_tree : (('a, 'b) Ptree.p_tree * float) -> unit
    method evaluate : (('a, 'b) Ptree.p_tree * float) list -> unit
    method results : (('a, 'b) Ptree.p_tree * float) list -> unit
end

class ['a, 'b] do_nothing : ['a, 'b] search_manager_sampler = object
    method init _ = ()
    method clone = ({<>} :> ('a, 'b) search_manager_sampler)
    method process _ _ _ _ _ _ _ _ _ = 
        StatusCommon.process_parallel_messages Status.user_message
    method any_trees _ = ()
    method next_tree _  = ()
    method evaluate _ = ()
    method results _ = ()
end

class ['a, 'b] composer (f : ('a, 'b) search_manager_sampler) 
    (s : ('a, 'b) search_manager_sampler) : 
        ['a, 'b] search_manager_sampler = object

    val f = f
    val s = s

    method init a = 
        f#init a;
        s#init a

    method clone = ({<f = f#clone; s = s#clone>} :> ('a, 'b) search_manager_sampler)

    method process z a b n c d y e x = 
        f#process z a b n c d y e x;
        s#process z a b n c d y e x;

    method any_trees a =
        f#any_trees a;
        s#any_trees a

    method next_tree a =
        f#next_tree a;
        s#next_tree a

    method evaluate a =
        f#evaluate a;
        s#evaluate a

    method results a =
        f#results a;
        s#results a

end

module MakeApp (Node : NodeSig.S) 
    (Edge : Edge.EdgeSig with type n = Node.n) 
    = struct

    class print_next_tree printer : 
        [Node.n, Edge.e] search_manager_sampler = object
        inherit [Node.n, Edge.e] do_nothing
        method next_tree (tree, cost) = printer (tree.Ptree.tree, cost, ())
    end

    class ['a, 'b] general_local_optimum_holder 
        (preprocess : (Tree.u_tree * float * Ptree.clade_cost) list ->
            unit) (postprocess : unit -> unit) 
            :
        ['a, 'b] search_manager_sampler = object
        val mutable current_best = None
        val mutable current_cost = None

        method init trees =
            let _ = 
                let t =
                    List.map (fun (a, b, c) -> a.Ptree.tree, b, c)
                        (match current_best with
                        | Some t -> t 
                        | None ->  trees)
                in
                preprocess t
            in
            current_best <- Some trees;
            current_cost <- None

        method clone = ({<>} :> ('a, 'b) search_manager_sampler)

        method process _ _ _ _ _ tree _ _ cost = 
            match tree with
            | None -> ()
            | Some tree ->
                    match cost with
                    | None -> ()
                    | Some c ->
                           match current_cost with
                           | None -> 
                                   current_cost <- cost;
                                   current_best <- Some [(tree, c, Ptree.NoCost)]
                           | Some cc ->
                                   if c < cc then begin
                                       current_cost <- cost;
                                       current_best <- 
                                           Some [(tree, c, Ptree.NoCost)];
                                       postprocess ();
                                       preprocess [(tree.Ptree.tree, c,
                                       Ptree.NoCost)];
                                   end else ()

        method any_trees _ = ()
        method next_tree _ = ()
        method evaluate _ = ()
        method results _ = ()

    end

    class ['a, 'b] local_optimum_holder results = object (self)
        inherit ['a, 'b] general_local_optimum_holder 
        (fun x -> Stack.push x results.stack) 
        (fun () -> ignore (Stack.pop results.stack))

        method clone = ({<>} :> (Node.n, Edge.e) search_manager_sampler)
    end

    (** A class that executes the function [printer] on every visited tree,
    * which is generated by joining the necessary positions using the function
    * [join_tree]. This sampler is expensive to use. *)
    class ['a, 'b] visited join_tree printer : ['a, 'b] search_manager_sampler = 
        object (self)
        inherit ['a, 'b] general_local_optimum_holder 
        (fun _ -> ())
        (fun _ -> ()) 

        method clone = ({<>} :> ('a, 'b) search_manager_sampler)

        method process incremental j1 j2 _ tree _ _ _ _ = 
            let tree = join_tree incremental j1 j2 tree in
            let cost = Ptree.get_cost `Adjusted tree in
            printer (tree.Ptree.tree, cost, ())
    end

    class virtual timed_trajectory time = object (self)
            inherit [Node.n, Edge.e] do_nothing 

            val mutable timer = Timer.start ()
            method virtual private timed_operation : 
                (Node.n, Edge.e) Ptree.p_tree * float -> unit

            method next_tree input =
                if time < Timer.wall timer then begin
                    self#timed_operation input;
                    timer <- Timer.start ();
                end else ()
    end

    class timed_printout queue time print :
        [Node.n, Edge.e] search_manager_sampler = object (self)
            inherit timed_trajectory time

            method private print_stack =
                Stack.iter (List.iter print) queue.stack;
                timer <- Timer.start ();

            method clone = ({<>} :> ('a, 'b) search_manager_sampler)

            method private timed_operation _ =
                self#print_stack

            method process _ _ _ _ _ _ _ _ _ =
                if time < Timer.wall timer then 
                    self#print_stack
                else ()
    end

    class ['a] timed_cancellation time :
        [Node.n, Edge.e] search_manager_sampler = object 
            inherit timed_trajectory time

            method process _ _ _ _ _ _ _ _ _ =
                if time < Timer.wall timer then 
                    raise Methods.TimedOut
                else ()
            
            method private timed_operation _ =
                raise Methods.TimedOut
    end

end

module MakeRes (Node : NodeSig.S) (Edge : Edge.EdgeSig with type n = Node.n) = struct

    let proportion a b = (float_of_int a) /. (float_of_int b)

    (* The following function assumes that the sequence is a nucleotide sequence
    * with all combinations turned on *)
    let gap_saturation sequence =
        (* Assume that a gap is 16 *)
        let len = Sequence.length sequence
        and gaps = 
            Sequence.fold 
            (fun acc base -> if 0 <> base land 16 then acc + 1 else acc)
            0
            sequence
        in
        proportion gaps len

    (* Count the number of bits in the integer n *)
    let count_bits n =
        assert (n < 32);
        let rec counter cnt n =
            if n = 0 then cnt
            else counter (cnt + (n land 1)) (n lsr 1)
        in
        counter 0 n

    (* Calculate the fraction of positions in the sequence that show
    * polymorphism. *)
    let poly_saturation sequence n =
        let len = Sequence.length sequence 
        and poly =
            Sequence.fold 
            (fun acc base -> if n = count_bits base then acc + 1 else acc)
            0
            sequence
        in
        proportion poly len

    (* Calculate the fraction of positions in the sequence that show more than
    * one base *)
    let total_saturation sequence = 
        let len = Sequence.length sequence 
        and tot =
            Sequence.fold
            (fun acc base -> if 1 < count_bits base then acc + 1 else acc)
            0
            sequence
        in
        proportion tot len

    (* The type of the unions *)
    type union_node = {
        rooted : Node.Union.u; (* The node of the subtree I root *)
        self : Node.Union.u;   (* The node of myself for lower unions *)
        depth : int;
    }

    let union_depth ptree depth = 
        let create new_tree nd =
            match Ptree.get_node nd ptree with
            | Tree.Single _ -> new_tree
            | Tree.Leaf _ -> 
                    let node = Ptree.get_node_data nd ptree in
                    let node = Node.Union.leaf None None node in
                    let new_node = { rooted = node; self = node; depth = 0 } in
                    Ptree.add_node_data nd new_node new_tree
            | Tree.Interior (_, par, c1, c2) ->
                    let nc1 = Ptree.get_node_data c1 new_tree
                    and nc2 = Ptree.get_node_data c2 new_tree
                    and ndd = Ptree.get_node_data nd ptree in
                    let median = 
                        Node.Union.union (Some par) ndd nc1.self nc2.self 
                    in
                    let mod_depth = (Node.num_height (Some par) ndd) mod depth in
                    let self_median =
                        if 0 <> mod_depth then median
                        else Node.Union.leaf None (Some par) ndd
                    in
                    let new_node = {
                        rooted = median;
                        self = self_median;
                        depth = mod_depth;
                    } in
                    Ptree.add_node_data nd new_node new_tree
        in
        let rec post_order_fold ptree vertex f acc =
            match Ptree.get_node vertex ptree with
            | Tree.Leaf _ ->
                    f acc vertex
            | Tree.Interior (_, _, c1, c2) ->
                    let acc = post_order_fold ptree c1 f acc in
                    let acc = post_order_fold ptree c2 f acc in
                    f acc vertex
            | Tree.Single _ -> f acc vertex
        in
        let post_order_fold_handle ptree f handle acc = 
            match Ptree.get_node handle ptree with
            | Tree.Leaf (_, b) ->
                    let acc = f acc handle in
                    post_order_fold ptree b f acc
            | Tree.Interior (a, b, _, _) ->
                    let acc = post_order_fold ptree a f acc in
                    post_order_fold ptree b f acc
            | Tree.Single _ -> f acc handle
        in
        let post_order_fold_tree ptree f acc = 
            let handles = Ptree.get_handles ptree in
            All_sets.Integers.fold (post_order_fold_handle ptree f) handles acc
        in
        let initial_tree = { 
            Ptree.node_data = All_sets.IntegerMap.empty;
            edge_data = Tree.EdgeMap.empty;
            tree = ptree.Ptree.tree;
            component_root = All_sets.IntegerMap.empty;
            origin_cost = 0.0;
        }
        in
        post_order_fold_tree ptree create initial_tree

    let union_distance parent node tree code =
        let cnode = Ptree.get_node_data code tree in
        Node.Union.distance_node parent node cnode.rooted

    let is_handle_or_other_handle node ptree =
        Ptree.is_handle node ptree ||
        (match Ptree.get_node node ptree with
        | Tree.Leaf (_, b)
        | Tree.Interior (_, b, _, _) ->
                Ptree.is_handle b ptree
        | Tree.Single _ -> failwith "IMpossible")

    let rec cluster_union node ptree depth =
        if is_handle_or_other_handle node ptree then node
        else 
            let nd = Ptree.get_node_data node ptree in
            if nd.depth = depth then node
            else 
                match Ptree.get_node node ptree with
                | Tree.Leaf (_, b) 
                | Tree.Interior (_, b, _, _) -> cluster_union b ptree depth
                | Tree.Single _ -> failwith "IMpossible"

    let union_distance_cluster parent node depth tree code =
        union_distance parent node tree (cluster_union code tree depth)

    let get_all_sequences_data node = 
        let sequences = Node.get_sequences None node in
        List.map (fun (_, a, _, _, _) -> a) sequences

    let get_sequence_data node = 
        match get_all_sequences_data node with
        | hd :: _ -> hd
        | _ -> failwith "No sequence data loaded"

    let preliminary_gap_saturation node = 
        let seq = get_sequence_data node in
        string_of_float (gap_saturation seq)

    let preliminary_poly_saturation node =
        let seq = get_sequence_data node in
        string_of_float (total_saturation seq)

    class ['a] union_table depth print = object (self)
        inherit [Node.n, 'a] do_nothing 

        val mutable current_tree = None
        val results = Array.create 20 ""
        val mutable union_tree = None

        method private update_tree new_tree = 
            match union_tree with
            | None -> 
                    let res = union_depth new_tree depth in
                    union_tree <- Some res;
                    res
            | Some v ->
                        if v.Ptree.tree == new_tree.Ptree.tree then v
                        else 
                            let res = union_depth new_tree depth in
                            let _ = union_tree <- Some res in
                            res

        method process _ edge _ vertex tree _ delta cost real_cost = 
            let str_delta = string_of_float delta in
            let union_tree = self#update_tree tree in
            let _ =
                match edge with
                | Tree.Single_Jxn vtx -> (* There is nothing to print here *)
                        for i = 0 to 19 do
                            results.(i) <- "";
                        done;
                | Tree.Edge_Jxn (a, b) ->
                        let e = Tree.Edge (a, b) in
                        let (Tree.Edge (a, b)) = 
                            Tree.normalize_edge e tree.Ptree.tree 
                        in
                        let adata = Ptree.get_node_data a tree 
                        and bdata = Ptree.get_node_data b tree
                        and a_union = Ptree.get_node_data a union_tree
                        and b_union = Ptree.get_node_data b union_tree in
                        let vertex_to_a_self =
                            Node.Union.distance_node (Some b) vertex 
                            a_union.self
                        and vertex_to_b_self =
                            Node.Union.distance_node (Some a) vertex 
                            b_union.self
                        and vertex_to_a_root =
                            Node.Union.distance_node (Some b) vertex 
                            a_union.rooted 
                        and vertex_to_b_root = 
                            Node.Union.distance_node (Some a) vertex 
                            b_union.rooted
                        in
                        results.(0) <-
                            (match real_cost with
                            | None -> ""
                            | Some v -> string_of_float v);
                        results.(1) <- string_of_float cost;
                        results.(2) <- string_of_int 
                            (Node.num_height (Some (Ptree.get_parent b tree)) bdata);
                        results.(3) <- string_of_int b_union.depth;
                        results.(4) <- preliminary_gap_saturation vertex;
                        results.(7) <- preliminary_gap_saturation adata;
                        results.(8) <- preliminary_gap_saturation bdata;
                        results.(9) <- preliminary_poly_saturation vertex;
                        results.(12) <- preliminary_poly_saturation adata;
                        results.(13) <- preliminary_poly_saturation bdata;
                        results.(14) <- str_delta;
                        results.(15) <- string_of_float vertex_to_a_self;
                        results.(16) <- string_of_float vertex_to_a_root;
                        results.(17) <- string_of_float vertex_to_b_self;
                        results.(18) <- string_of_float vertex_to_b_root;
                        results.(19) <- string_of_int depth;
            in
            Array.iter (fun x -> print x; print "\t") results;
            print "\n";

    end

    (*
    class ['a] union_root_distribution print = object (self)
        inherit [Node.n, 'a] do_nothing 
        method next_tree (tree, cost) =
            let roots = Ptree.get_handles tree in
            match All_sets.Integers.elements roots with
            | [hd] ->
                    let root = Ptree.get_component_root hd tree in
                    (match root.Ptree.root_median with
                    | Some (_, b) ->
                            let seqb = get_all_sequences_data b in
                            let seqb = 
                                List.map (fun x -> x.Sequence.Unions.seq) seqb 
                            in
                            List.iter (fun seqb ->
                                Sequence.iter (fun x -> print (string_of_int
                                (count_bits x))) seqb;
                                print "\t";
                                Sequence.iter (fun x -> print (string_of_int 
                                (if 16 = (x land 16) then 1 else 0))) seqb;
                                print "\n";)
                            seqb
                    | None -> failwith "no root?")
            | _ -> failwith "Please update the union_root_distribution to
            handle forests. class Sampler.Res.union_root_distribution"
    end
    *)

    class ['a] tests_before_next print = 
        object (self) 
        inherit [ Node.n, 'a] do_nothing

        val mutable counter = 0

        method any_trees x =
            if x = false then begin
                print (string_of_int counter);
                print "\n%!"
            end else ()

        method next_tree _ =
                    print (string_of_int counter ^ "\t");
                    counter <- 0

        method process _ edge _ _ _ _ _ _ x = 
            match x with
            | Some _ ->
                    print (string_of_int counter ^ "\t");
                    counter <- 0
            | None ->
                    counter <- counter + 1
    end

    class ['a] break_n_join_distances (join_fn : (Node.n, 'a) Ptree.join_fn) 
    print =
        object (self)
            inherit [Node.n, 'a] do_nothing

            val mutable cost = "unknown"

            method clone = ({<>} :> (Node.n, 'b) search_manager_sampler)

            method init lst =
                match lst with
                | [(h, _, _)] -> 
                        cost <- string_of_float (Ptree.get_cost `Adjusted h)
                | _ -> ()

            method process incremental j1 j2 _ pt _ bd _ _ =
                let total_cost = 
                    let joined, _ = join_fn incremental j1 j2 pt in
                    let total_cost = Ptree.get_cost `Adjusted joined in
                    string_of_float total_cost
                in
                print ((string_of_float bd) ^ "\t" ^ cost ^ "\t" ^ total_cost ^ "\n")
        end

end

