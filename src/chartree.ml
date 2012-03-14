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

(* $Id: chartree.ml 2871 2008-05-23 17:48:34Z andres $ *)
let () = SadmanOutput.register "Chartree" "$Revision: 2871 $"

let info_user_message format =
    Printf.ksprintf (Status.user_message Status.Information) format

(** chartree.ml *)
type 'a p_tree = (Node.node_data, 'a) Ptree.p_tree

type incremental = [
    | `Children of int
    | `No_Children of int
    | `HandleC of (int * int)
    | `HandleNC of (int * int)
]

module Handles = All_sets.Integers

let debug_print_on_break = false
let debug_exclude = false
let debug_breakfn_deltas = false
let debug_joinfn = false
let debug_costfn = true
let debug_uppass_which_handle = false
let debug_no_incremental = false
let debug_diagnosis = false
let odebug = Status.user_message Status.Information

let check_assertion_two_nbrs a b c =
    if a <> Tree.get_id b then true
    else 
        let _ = Status.user_message Status.Error c in
        false

module IntSet = All_sets.Integers
let fprintf = Printf.fprintf 
let (-->) a b = b a

let tree_iter fn {Ptree.tree=tree} =
    Handles.fold
        (fun handle ptree ->
             Tree.pre_order_node_visit
                 (fun parent_id self_id curr_tree ->
                      fn parent_id self_id;
                      (Tree.Continue, ()))
                 handle
                 tree
                 ptree)
        (Tree.get_handles tree)
        ()

let tree_iter_depth fn {Ptree.tree=tree} =
    let rec get_depth self_id =
        if Tree.is_handle self_id tree
        then 0
        else 1 + get_depth (Tree.get_parent self_id tree)
    in
    Handles.fold
        (fun handle () ->
             Tree.pre_order_node_visit
                 (fun parent_id self_id () ->
                      let curr = get_depth self_id in
                      fn parent_id self_id curr;
                      (Tree.Continue, ()))
                 handle
                 tree
                 ();
             ())
        (Tree.get_handles tree)
        ()


let print_nodes tree =
    tree_iter (fun par_id self_id ->
                   let node = Ptree.get_node_data self_id tree in
                   let node_type = match Ptree.get_node self_id tree with
                       | Tree.Leaf (x,a)         -> Printf.sprintf "leaf(%d)-%d" x a
                       | Tree.Interior (x,a,b,c) -> Printf.sprintf "intr(%d)-(%d,%d,%d)" x a b c
                       | Tree.Single x           -> Printf.sprintf "sing(%d)" x
                   in
                   let par_string = match par_id with
                       | None -> "none"
                       | Some a -> string_of_int a 
                   in
                   print_endline ("Node " ^ (string_of_int self_id)
                                  ^ " is " ^ node_type ^ " with"
                                  ^ " parent " ^ par_string
                                  ^ " cost " ^ string_of_float
                                  (Node.Standard.total_cost par_id node)))
                tree

let print_node_data_indent ?(indent = 3) trees =
    tree_iter_depth
        (fun par_id self_id depth ->
             print_string (String.make (depth * indent) ' ');
             print_endline (Node.Standard.to_string
                                (Ptree.get_node_data
                                     self_id
                                     trees)))
        trees


let hashdoublefind tree node_ids = None 
(* we don't use likelihood in chartree, so no extra data is neccessary *)

let downpass_handle handle ({Ptree.tree=tree} as ptree) =
    Tree.post_order_node_visit
        (fun parent_id self_id curr_tree ->
            (* assumption: tree topology has not changed *)
            let node = Tree.get_node self_id tree in
            match node with
            | Tree.Single self
            | Tree.Leaf (self, _) ->
                let selfdata = Ptree.get_node_data self curr_tree in
                let selfdata = Node.Standard.update_leaf selfdata in
                (Tree.Continue, Ptree.add_node_data self selfdata curr_tree)
            | Tree.Interior (nid, ch1id, ch2id, ch3id) ->
                let ch1id, ch2id = match parent_id with
                    | None -> ch2id, ch3id
                    | Some parent_id ->
                        assert (check_assertion_two_nbrs parent_id node "9");
                        Tree.other_two_nbrs parent_id node
                in
                Printf.printf "Calculating median between %d and %d\n%!" ch1id ch2id;
                let ch1data = Ptree.get_node_data ch1id curr_tree in
                let ch2data = Ptree.get_node_data ch2id curr_tree in
                let selfdata =
                    try Some (Ptree.get_node_data nid curr_tree)
                    with Not_found -> None
                in
                let median = match hashdoublefind tree [ch1id;ch2id] with
                    | Some x -> Node.Standard.median ~branches:x (Some nid) selfdata ch1data ch2data
                    | None -> Node.Standard.median (Some nid) selfdata ch1data ch2data
                in
                (Tree.Continue, Ptree.add_node_data nid median curr_tree))
        handle
        tree
        ptree

let downpass ({Ptree.tree=tree} as ptree) =
    if debug_diagnosis then Printf.printf "Downpass Begins\n%!";
    Handles.fold
        (downpass_handle)
        (Tree.get_handles tree)
        ptree

let uppass_handle handle ({Ptree.tree=tree} as ptree) =
    if debug_diagnosis then Printf.printf "Uppass Begins\n%!";
    let ptree = 
        let d =
            match Ptree.get_node handle ptree with
            | Tree.Single _ -> 
                    let d = Ptree.get_node_data handle ptree in
                    Some ((`Single handle), d)
            | _ -> None
        in
        Ptree.assign_root_to_connected_component handle d 0. None ptree 
    in
    let virt_root : Node.node_data Ptree.root_node ref = ref None in
    let get_parent myid id tree =
        match !virt_root with
        | Some (`Edge (rid, rchid), rnode) ->
              if rid = id && rchid = myid then rnode
              else Ptree.get_node_data id tree
        | None -> failwith "Out-of-sequence uppass"
        | Some (`Single node, rnode) ->
                failwith "Unexpected Chartree.uppass_handle.get_parent"
    in
    let get_node = Ptree.get_node_data in
    Tree.pre_order_node_visit
        (fun parent_id self_id ptree ->
             (* assumption: tree topology has not changed *)
             let node = Tree.get_node self_id tree in
             match node with
             | Tree.Single _ -> 
                     let data = Ptree.get_node_data self_id ptree in
                     let ptree = 
                         Ptree.assign_root_to_connected_component
                            self_id (Some ((`Single self_id), data)) 0. None ptree
                     in
                     (Tree.Continue, ptree)
             | Tree.Leaf (selfid, otherid) -> begin
                   match parent_id with
                   | Some p ->
                         assert (p = otherid);
                         (* For OTUs, we perform a 3-median to handle the case
                            of polymorphic states... *)
                         let mydata = Ptree.get_node_data selfid ptree in
                         let otherdata = Ptree.get_node_data otherid ptree in
                         if debug_diagnosis then
                             info_user_message "Final States Leaf: %d, %d" selfid otherid;
                         let median3 = Node.Standard.final_states None otherdata mydata mydata mydata in
                         let ptree = Ptree.add_node_data selfid median3 ptree in
                         (Tree.Continue, ptree) (* OTU *)
                   | None -> begin
                         (* We're dealing with a handle node that's also a leaf.
                            We want to calculate the last median and store it. *)
                         let mydata = Ptree.get_node_data selfid ptree in
                         let otherdata = Ptree.get_node_data otherid ptree in
                         let root_prelim = Node.Standard.median None None mydata otherdata in
                         virt_root := Some (`Edge (selfid, otherid), root_prelim);
                         let tree_cost = Node.Standard.tree_cost None root_prelim in
                         if debug_uppass_which_handle then odebug "uppass: handle is leaf";
                         if debug_diagnosis then
                             info_user_message "Final States Leaf: %d, %d" selfid otherid;
                         let mydata =
                             Node.Standard.final_states None root_prelim mydata mydata mydata 
                         in
                         let ptree = Ptree.add_node_data selfid mydata ptree in
                         let ptree = 
                             Ptree.assign_root_to_connected_component 
                                        self_id !virt_root tree_cost None ptree
                         in
                         (Tree.Continue, ptree)
                     end
               end
             | Tree.Interior (nid, ch1id, ch2id, ch3id) ->
                   match parent_id with
                   | None -> begin
                         (* We're dealing with a handle node.  We calculate a
                            median for the handle and its "parent"; then, we
                            calculate the handle's final state.  WARNING: my
                            parent does not get assigned its final state. *)
                         let mydata = get_node nid ptree in
                         let pardata = get_node ch1id ptree in
                         let ch1data = get_node ch2id ptree in
                         let ch2data = get_node ch3id ptree in
                         if debug_uppass_which_handle then odebug "uppass: handle is internal";
                         let root_prelim = Node.Standard.median None None mydata pardata in
                         let tree_cost = Node.Standard.tree_cost None root_prelim in
                         virt_root := Some (`Edge (nid, ch1id), root_prelim);
                         if debug_diagnosis then
                             info_user_message "Final States Node: %d -(%d,%d,%d)" nid ch1id ch2id ch3id;
                         let median3 = 
                             Node.Standard.final_states None root_prelim mydata ch1data ch2data
                         in
                         let ptree = Ptree.add_node_data nid median3 ptree in
                         let ptree = Ptree.assign_root_to_connected_component 
                                            self_id !virt_root tree_cost None ptree
                         in
                         (Tree.Continue, ptree)
                     end
                   | Some parent_id -> begin
                         (* We should treat the case of the node which is the
                            handle's parent separately; we now do this. *)
                         let ch1id, ch2id = 
                             assert (check_assertion_two_nbrs parent_id node "0");
                             Tree.other_two_nbrs parent_id node 
                         in
                         let mydata = get_node nid ptree in
                         let pardata = get_parent nid parent_id ptree in
                         let ch1data = get_node ch1id ptree in
                         let ch2data = get_node ch2id ptree in
                         if debug_diagnosis then
                             info_user_message "Final States Node: %d -(%d,%d,%d)" nid ch1id ch2id ch3id;
                         let median3 =
                             Node.Standard.final_states None pardata mydata ch1data ch2data
                         in
                         (Tree.Continue, Ptree.add_node_data nid median3 ptree)
                    end)
        handle
        tree
        ptree

let uppass ({Ptree.tree=tree} as ptree) =
    Handles.fold
        uppass_handle
        (Tree.get_handles tree)
        ptree

(** [downpass_cost_step ptree node_id c1 c2] will perform a downpass like step
 * in the node with id [node_id], which has children [c1] and [c2], where the only
 * data element from its Node.Standard.node_data to be updated is its cost, not the actual
 * preliminary character states.  The function assumes that [node_id] has
 * information already associated, and that the median it has assigned is a valid
 * one, already assigned. *)
let downpass_cost_step ptree node_id c1 c2 =
    let c1d = 
        try Ptree.get_node_data c1 ptree
        with Not_found -> failwith "Chartree.downpass_step c1d not found"
    and c2d = 
        try Ptree.get_node_data c2 ptree 
        with Not_found -> failwith "Chartree.downpass_step c2d not found"
    in
    let selfdata = 
        try 
            let data = (Ptree.get_node_data node_id ptree) in
            let new_cost = 
                data.Node.node_cost +.
                c1d.Node.total_cost +.
                c2d.Node.total_cost 
            and new_exclude_info = 
                (* We can use none because Node.Standard do not use this
                * information *)
                Node.Standard.excludes_median None c1d c2d 
            in
            let new_excluded = Node.Standard.has_excluded new_exclude_info in
            let new_cost = 
                if new_excluded then infinity
                else new_cost 
            in
            Ptree.add_node_data node_id { data with 
                Node.total_cost = new_cost;
                Node.exclude_info = new_exclude_info } ptree
        with 
        | err -> 
              let msg = 
                  ("Failed to complete on " ^ string_of_int node_id)
              in
              Status.user_message Status.Error msg;
              raise err
    in
    selfdata

(** [downpass_step ptree node_id c1 c2] recalculates the preliminary states of
 * the [node_id] data using the data from it's two children [c1] and [c2].  The
 * function requires that both children have already been assigned data (a downpass
 * has been performed over them). The function returns a tuple, the newly
 * calculated data and the old data assigned to node_id if present, otherwise the
 * second element of the tuple is None *)
let downpass_step ptree node_id c1 c2 =
    let c1d = 
        try Ptree.get_node_data c1 ptree
        with Not_found -> failwith "Chartree.downpass_step c1d not found"
    and c2d = 
        try Ptree.get_node_data c2 ptree 
        with Not_found -> failwith "Chartree.downpass_step c2d not found"
    and selfdata = 
        match node_id with
        | None -> None
        | Some x -> 
              try 
                  Some (Ptree.get_node_data x ptree)
              with 
              | Not_found -> None
              | err -> 
                    let msg = 
                        ("Failed to complete on " ^ string_of_int x)
                    in
                    Status.user_message Status.Error msg;
                    raise err
    in
    (Node.Standard.median node_id selfdata c1d c2d), selfdata

(** [iterate e (a,b)r t] iterates the branch lengths, for likelihood characters.
 * [e] is returned as the edges for the nodes that need to be updated because of
 * changes in the tree. [r] is the root of the tree, with [a] and [b] being the
 * nodes surrounding it. [t] is the tree. and this function returns a new ptree
 * with the updates applied. *)
let iterate edges ((Tree.Edge (a, b)) as root) ptree =
    let visitor c1 c2 (ptree, edgeset) =
        let par = Ptree.get_parent c1 ptree in
        assert (par = Ptree.get_parent c2 ptree);
        let gpar = Ptree.get_parent par ptree in
        assert( gpar <> c1 && gpar <> c2 );
        if All_sets.Integers.mem par edgeset then
            let gd id = Ptree.get_node_data id ptree in
            let pard = 
                Node.Standard.edge_iterator (Some (gd gpar)) (gd par) (gd c1) (gd c2)
            in
            let ptree = Ptree.add_node_data par pard ptree in
            let edgeset =
                All_sets.Integers.add (Ptree.get_parent par ptree) edgeset
            in
            (ptree, edgeset)
        else 
            (ptree, edgeset) in
    let tree = ptree.Ptree.tree in
    let ptree, _ = 
        Tree.post_order_node_with_edge_visit_simple visitor root tree (ptree,edges)
    in
    ptree 

(** The incremental downpass type. If a function returns [Same], then the
 * downpass over a particular node didn't cause any change in the tree. If the
 * returned value is [Cost tree], then a new tree where only the cost of that node
 * was modified is returned. If the data of the node was affected, then an [All
 * tree] is returned. *)
type 'a id = 
    | Same 
    | Cost of (Node.node_data, 'a) Ptree.p_tree 
    | All of (Node.node_data, 'a) Ptree.p_tree 

(** A set of flags used to pass messages between functions, so that they know
 * weather or not they should do nothing ([MSame]), update cost ([MCost]), or
 * update all data ([MAll]). *)
type tdp = MSame | MCost | MAll

(** [simple_downpass continuation ptree node_id c1 c2] is the main incremental
 * downpass function. Provides functionality to signal calling functions weather or
 * not it is necessary to continue with a downpass. It will attempt to perform a
 * downpass on [node_id] with children [c1] and [c2] for which downpass data has
 * already been calculated. According to the [continuation] flag, the appropriate
 * action is performed.  *)
let simple_downpass ?(debug_joinfn=debug_joinfn) continuation ptree node_id c1 c2 =
    match continuation with
    | MSame -> Same
    | MCost -> Cost (downpass_cost_step ptree node_id c1 c2) 
    | _ ->
            (* The function to compare data *)
          let cd = Node.compare_downpass in
          let median, selfd = 
              downpass_step ptree (Some node_id) c1 c2 
          in
          let costa = median.Node.total_cost in
          match selfd with
          | Some s when ((0 = cd median s) 
                         && (costa = s.Node.total_cost)
                         && (0 = compare median.Node.exclude_info s.Node.exclude_info)) -> 
                Same
          | Some s when (0 <> cd median s) ->
                if debug_joinfn then begin
                    Printf.printf "The new some cost of vertex %d with\
                        chidren %d and %d is %f\n%!" node_id c1 c2 costa;
                    Printf.printf "The old some cost of vertex %d is %f\n%!" node_id
                        s.Node.total_cost;
                end;
                All (Ptree.add_node_data node_id median ptree)
          | Some s ->
                if debug_joinfn then begin
                    Printf.printf "The new cost cost of vertex %d with\
                        chidren %d and %d is %f\n%!" node_id c1 c2 costa;
                    Printf.printf "The old cost cost of vertex %d is %f\n%!" node_id
                        s.Node.total_cost;
                    Printf.printf "Only the cost has changed";
                end;
                Cost (Ptree.add_node_data node_id median ptree)
          | None ->
                if debug_joinfn then
                    Printf.printf "The new none cost of vertex %d is %f\n%!" node_id
                        costa;
                All (Ptree.add_node_data node_id median ptree)

(** [calculate_root ptree node_id neighbor] calculates a new virtual root for
 * the tree, when [node_id] or [neighbor] is the handle. The function requires that
 * both node_id and neighbor have already calculated their downpass data. *)
let calculate_root ptree node_id neighbor =
    match downpass_step ptree None node_id neighbor with
    | median, None ->
          if Ptree.is_handle node_id ptree then
              (let root_median = Some (`Edge (node_id, neighbor), median) in
               if debug_joinfn then 
                   Printf.printf "The cost of this root is %f\n%!"
                       median.Node.total_cost;
               let cost = Node.Standard.tree_cost None median in
               Ptree.assign_root_to_connected_component node_id root_median
                   cost None ptree)
          else if Ptree.is_handle neighbor ptree then
              (if debug_joinfn then 
                   Printf.printf "The cost of this root is %f\n%!"
                       median.Node.total_cost;
               let root_median = Some (`Edge (neighbor, node_id), median) in
               let cost = Node.Standard.tree_cost None median in
               Ptree.assign_root_to_connected_component neighbor root_median
                    cost None ptree)
          else failwith "Chartree.calculate_root 1"
    | median, Some _ -> 
            (* As downpass_step received a None as node_id data, we don't expect
             * any data to come from any node_id in the second element in the
             * tuple *)
            failwith "Chartree.calculate_root"

(** [incremental_downpass_step continuation ptree node_id] will attempt a
 * downpass step over the vertex with id [node_id] in the tree [ptree]. Depending
 * on the kind of vertex [node_id] is, a different operation is performed. *)
let incremental_downpass_step continuation ptree node_id = 
    match Tree.get_node node_id ptree.Ptree.tree with
    | Tree.Interior (_, par, c1, c2) -> 
            (* We should try to continue, as an interior node could change if it
             * has new data in it's children *)
            if debug_joinfn then 
                Printf.printf "The children of %d are %d and %d and the parent \
                is %d\n%!" node_id c1 c2 par;
            simple_downpass continuation ptree node_id c1 c2
      | _ -> Same

(** [add_to_accumlator acc c1 c2 nd previous] will add the appropriate vertices
 * to the incremental uppass list of vertices, when an update in the node [nd] with
 * children [c1] and [c2] is updated. If [previous] is None, then none of the
 * children of [nd] are already present in [acc], in which case they are both
 * added, otherwise, [previous] is compared with [c1] and [c2] to add the missing
 * one to [acc] with the [`No_Children] flag. The function assumes that if
 * [previous] has been already added, this is due to the fact that it was also
 * evaluated during an incremental downpass, therefore, it's children are already
 * part of [acc]. *)
let add_to_accumlator cont acc c1 c2 nd previous =
    match previous with
    | Some v ->
          (match cont with
           | MSame 
           | MCost -> (`Children nd) ::  acc
           | _ ->
                 (`Children nd) ::
                     (if v = c1 then (`No_Children c2) :: acc
                      else if v = c2 then (`No_Children c1) :: acc
                      else (`No_Children c1) :: (`No_Children c2) :: acc))
    | None ->
          (`Children nd) ::
              (`No_Children c1) :: 
              (`No_Children c2) :: 
              acc

(** Finally we are ready for the set of functions that perform downpasses over
 * trees, we will start with those that force a full downpass over it. These are
 * used for consistency checks more than anything else. *)
let force_downpass ptree node_id =
    match Tree.get_node node_id ptree.Ptree.tree with
    | Tree.Single _ | Tree.Leaf _ -> ptree
    | Tree.Interior (_, parent, c1, c2) ->
          match simple_downpass MAll ptree node_id c1 c2 with
          | Same -> ptree
          | Cost ptree 
          | All ptree -> ptree

let rec force_downpass_node ptree node_id =
    match Tree.get_node node_id ptree.Ptree.tree with
    | Tree.Interior (_, par, c1, c2) ->
          let ptree = force_downpass_node ptree c1 in
          let ptree = force_downpass_node ptree c2 in
          let ptree = force_downpass ptree node_id in
          ptree
    | Tree.Leaf _ | Tree.Single _ -> ptree

let force_downpass_handle handle ptree =
    assert (Tree.is_handle handle ptree.Ptree.tree);
    match Tree.get_node handle ptree.Ptree.tree with
    | Tree.Interior (_, neigh, c1, c2) ->
          let ptree = force_downpass_node ptree handle in
          let ptree = force_downpass_node ptree neigh in
          calculate_root ptree handle neigh 
    | Tree.Leaf (_, neigh) ->
          let ptree = force_downpass_node ptree neigh in
          calculate_root ptree handle neigh
    | Tree.Single _ -> ptree

let force_downpass ptree = 
    Handles.fold force_downpass_handle (Tree.get_handles ptree.Ptree.tree) ptree
        

let is_root_with_handle node ptree =
    match Tree.get_node node ptree.Ptree.tree with
    | Tree.Interior (_, par, _, _)
    | Tree.Leaf (_, par) ->
          if Tree.is_handle par ptree.Ptree.tree then begin
              match Tree.get_node par ptree.Ptree.tree with
              | Tree.Interior (_, par, _, _)
              | Tree.Leaf (_, par) -> node = par
              | Tree.Single _ -> assert (false);
          end else false
    | Tree.Single _ -> assert (false)


(* True incremental downpass functions *)

(** [aux_incremental_downpass continuation acc prev ptree node_id] is the main
* function for incremental downpass. It is called aux_incremental_downpass only
* because it is not visible from outside and should never be (if it is, is just
* for debugging purposes!). 
*
* This function will attempt to perform an incremental downpass in the node
* [node_id] of the tree [ptree] using the kind of modification defined in
* [continuation]. The uppass list for the incremental uppass is held in [acc],
* and the previously calculated vertex (if this is not the first call to
* aux_incremental_downpass but a recursive one), is [prev]. 
* The function returns a tuple, holding the resulting tree from the incremental
* downpass, and the list of vertices that must be visited by an incremental
* uppass. *)
let rec aux_incremental_downpass continuation acc prev ptree node_id =
    let ptree, continue =
        (* We just pass the continuation parameter to incremental_downpass_step
        * to make the appropriate thing about it. *)
        match incremental_downpass_step continuation ptree node_id with
        | Same -> ptree, MSame
        | Cost ptree -> ptree, MCost
        | All ptree -> ptree, MAll
    in
    if Tree.is_handle node_id ptree.Ptree.tree then begin 
        match Tree.get_node node_id ptree.Ptree.tree with
        | Tree.Interior (_, neigh, c1, c2) ->
                let acc = add_to_accumlator continue acc c1 c2 node_id prev in
                (calculate_root ptree node_id neigh,
                (match acc, prev with
                | (_ :: tl), None ->
                        (`HandleNC (node_id, node_id) ::
                            `HandleNC (neigh, node_id) ::
                                tl)
                | (_ :: tl), Some _ ->
                        (`HandleC (node_id, node_id) ::
                            `HandleNC (neigh, node_id) ::
                                tl)
                | [], _ -> 
                        failwith "Unexpected aux_incremental_downpass"))
        | Tree.Leaf (_, neigh) ->
                assert (prev = None);
                calculate_root ptree node_id neigh,
                (`HandleNC (node_id, node_id)) ::
                    (`HandleNC (neigh, node_id)) ::
                        acc
        | Tree.Single _ -> ptree, acc
    end else if is_root_with_handle node_id ptree then begin
        match Tree.get_node node_id ptree.Ptree.tree with
        | Tree.Interior (_, handle, c1, c2) ->
                let acc = add_to_accumlator continue acc c1 c2 node_id prev in
                (calculate_root ptree node_id handle,
                (match acc, prev with
                | (_ :: tl), None ->
                        (`HandleNC (node_id, handle)) ::
                            (`HandleNC (handle, handle)) ::
                                tl
                | (_ :: tl), Some _ ->
                        (`HandleC (node_id, handle)) ::
                            (`HandleNC (handle, handle)) ::
                                tl
                | [], _ -> 
                        failwith "Unexpected aux_incremental_downpass"))
        | Tree.Leaf (_, handle) ->
                assert (prev = None);
                aux_incremental_downpass continuation acc None ptree handle
        | Tree.Single _ -> assert (false);
    end else begin
        match Tree.get_node node_id ptree.Ptree.tree with
        | Tree.Interior (_, next, c1, c2) -> 
                let acc = add_to_accumlator continue acc c1 c2 node_id prev in
                (match continue with
                | MSame -> ptree, acc
                | _ ->
                        aux_incremental_downpass continue acc (Some node_id) 
                        ptree next)
        | Tree.Single _
        | Tree.Leaf _ -> ptree, acc
    end

let rec incremental_downpass_path_to_handle acc prev ptree = function
    | [] -> ptree, acc
    | [h] ->
            assert (Tree.is_handle h ptree.Ptree.tree);
            (match Tree.get_node h ptree.Ptree.tree with
            | Tree.Interior (_, neigh, c1, c2) ->
                    let new_h, _ = downpass_step ptree None c1 c2 in
                    let ptree = Ptree.add_node_data h new_h ptree 
                    and nacc = add_to_accumlator MAll acc c1 c2 h prev in
                    calculate_root ptree h neigh, 
                    (match nacc, prev with
                    | (_ :: tl), Some _ ->
                            (`HandleC (h, h) :: 
                                `HandleNC (neigh, h) ::
                                    tl)
                    | (_ :: tl), None ->
                            (`HandleNC (h, h) :: 
                                `HandleNC (neigh, h) ::
                                    tl)
                    | [], _ -> 
                            failwith "Unexpected \
                            incremental_downpass_path_to_handle")
            | Tree.Leaf (_, neigh) ->
                    calculate_root ptree h neigh, 
                    (`HandleC (h, h) :: 
                        `HandleNC (neigh, h) ::
                            acc)
            | Tree.Single _ -> ptree, acc)
    | [r; h] when (is_root_with_handle r ptree)-> (* h has to be the handle *)
            assert (Tree.is_handle h ptree.Ptree.tree);
            let r_node = Tree.get_node r ptree.Ptree.tree
            and h_node = Tree.get_node h ptree.Ptree.tree in
            let ptree = 
                match r_node, h_node with
                | Tree.Interior (_, _, c1r, c2r), 
                    Tree.Interior (_, _, c1h, c2h) ->
                        let new_r, _ = downpass_step ptree None c1r c2r in
                        let new_h, _ = downpass_step ptree None c1h c2h in
                        ptree --> Ptree.add_node_data r new_r -->
                            Ptree.add_node_data h new_h
                | _, Tree.Interior (node_id, _, c1, c2)
                | Tree.Interior (node_id, _, c1, c2), _ ->
                        let new_node, _ = downpass_step ptree None c1 c2 in
                        Ptree.add_node_data node_id new_node ptree 
                | _, _ -> ptree
            in
            calculate_root ptree h r,
            (match prev with
            | Some _ ->
                    (`HandleC (r, h) :: 
                        (`HandleNC (h, h)) ::
                            acc)
            | None ->
                    (`HandleNC (r, h) :: 
                        (`HandleNC (h, h)) ::
                            acc))
    | h :: t ->
            assert (not (Tree.is_handle h ptree.Ptree.tree));
            let new_prev = Some h in
            match Tree.get_node h ptree.Ptree.tree with
            | Tree.Single _ -> assert (false);
            | Tree.Leaf _ -> 
                    let acc = (`Children h) :: acc in
                    incremental_downpass_path_to_handle acc new_prev ptree t
            | Tree.Interior (_, _, c1, c2) ->
                    let new_h, _ = downpass_step ptree None c1 c2 in
                    let ptree = Ptree.add_node_data h new_h ptree in
                    let acc = add_to_accumlator MAll acc c1 c2 h prev in
                    incremental_downpass_path_to_handle acc new_prev ptree t


let simple_uppass ptree parentd selfd c1d c2d node_id =
    let median3 = Node.Standard.final_states None parentd selfd c1d c2d in
    if 0 = Node.compare_uppass selfd median3 then 
        None 
    else 
        Some (Ptree.add_node_data node_id median3 ptree)

let incremental_downpass node_id ptree =
    aux_incremental_downpass MAll [] None ptree node_id
    
let simple_uppass_in_leaf ptree parentd selfd node_id =
    let median3 = Node.Standard.final_states None parentd selfd selfd selfd in
    Ptree.add_node_data node_id median3 ptree

let incremental_uppass_step ptree node_id =
    let selfd = Ptree.get_node_data node_id ptree in
    match Tree.get_node node_id ptree.Ptree.tree with
    | Tree.Interior (nid, parent, c1, c2) ->
            let parentd = Ptree.get_node_data parent ptree 
            and c1d = Ptree.get_node_data c1 ptree 
            and c2d = Ptree.get_node_data c2 ptree in
            simple_uppass ptree parentd selfd c1d c2d node_id
    | Tree.Leaf (selfid, parent) ->
            let parentd = Ptree.get_node_data parent ptree in
            Some (simple_uppass_in_leaf ptree parentd selfd selfid)
    | Tree.Single _ -> None

let add_children node_id ptree lst =
    match Tree.get_node node_id ptree.Ptree.tree with
    | Tree.Interior (_, _, c1, c2) ->
            (`No_Children c1) :: (`No_Children c2) :: lst
    | Tree.Leaf _ | Tree.Single _ -> lst

let incr_uppass_step_in_root_child ptree node_id handle_id =
    try
        let root = 
            assert (Tree.is_handle handle_id ptree.Ptree.tree);
            try Ptree.get_component_root handle_id ptree
            with Not_found -> 
                failwith  "Unexpected \
                Chartree.incr_uppass_step_in_root_child a"
        and selfd =
            try Ptree.get_node_data node_id ptree
            with Not_found ->
                 failwith "Unexpected \
                 Chartree.incr_uppass_step_in_root_child b"
        in
        match root.Ptree.root_median with
        | Some (`Edge (a, b), root) ->
                assert ((a = node_id) || (b = node_id));
                (match Ptree.get_node node_id ptree with
                | Tree.Interior (_, _, c1, c2) ->
                        let c1d = Ptree.get_node_data c1 ptree
                        and c2d = Ptree.get_node_data c2 ptree in
                        simple_uppass ptree root selfd c1d c2d node_id
                | Tree.Leaf (_, _) ->
                        Some (simple_uppass_in_leaf ptree root selfd node_id)
                | Tree.Single _ -> Some ptree)
        | Some ((`Single id), root) ->
                Some (simple_uppass_in_leaf ptree selfd selfd node_id)
        | None -> failwith "Unexpected Chartree.incr_uppass_step_in_root_child"
    with
    | Not_found -> 
            failwith "Unexpected Chartree.incr_uppass_step_in_root_child 2"

let counter_iu = ref 0 
let show_stat = false

let print_id id =
    if show_stat then 
        Status.user_message Status.Information ("Uppass on " ^ string_of_int
        id)
    else ()

let rec incremental_uppass ptree = function
    | (`No_Children node_id) :: t ->
            incr counter_iu;
            print_id node_id;
            (match incremental_uppass_step ptree node_id with
            | None -> incremental_uppass ptree t
            | Some ptree -> 
                    let nt = add_children node_id ptree t in
                    incremental_uppass ptree nt)
    | (`Children node_id) :: t ->
            incr counter_iu;
            print_id node_id;
            (match incremental_uppass_step ptree node_id with
            | None -> incremental_uppass ptree t
            | Some ptree -> incremental_uppass ptree t)
    | (`HandleC (node_id, handle_id)) :: t ->
            incr counter_iu;
            print_id node_id;
            (match incr_uppass_step_in_root_child ptree node_id handle_id with
            | None -> incremental_uppass ptree t
            | Some ptree -> incremental_uppass ptree t)
    | (`HandleNC (node_id, handle_id)) :: t ->
            incr counter_iu;
            print_id node_id;
            (match incr_uppass_step_in_root_child ptree node_id handle_id with
            | None -> incremental_uppass ptree t
            | Some ptree ->
                    let nt = add_children node_id ptree t in
                    incremental_uppass ptree nt)
    | [] -> 
            if show_stat then begin
                let edges = Tree.EdgeSet.cardinal ptree.Ptree.tree.Tree.d_edges in
                let vertices = (edges + 3) / 2 in
                let fraction = (0.5 *. (float_of_int !counter_iu)) /. (float_of_int vertices) in
                Status.user_message Status.Information ("total fraction affected by
                uppass was " ^ string_of_float fraction ^ " with vertices " ^
                string_of_int vertices ^ " and uppass " ^ string_of_int
                !counter_iu);
                counter_iu := 0;
            end;
            ptree

let incremental_uppass a b =
    counter_iu := 0;
    incremental_uppass a b

(* cost_fn must have type key -> tree -> cost *)
let set_clade_root ptree clade_node new_clade_handle =
    let handle = Tree.int_of_id new_clade_handle 
    and cost = Node.Standard.tree_cost None clade_node in
    match Tree.get_node handle ptree.Ptree.tree with
    | Tree.Leaf (_, neigh)
    | Tree.Interior (_, neigh, _, _) ->
            Ptree.assign_root_to_connected_component handle 
            (Some (`Edge (handle, neigh), clade_node)) cost None ptree
    | Tree.Single _->
            Ptree.assign_root_to_connected_component handle None
            cost None ptree

let update_tree_data_break doup delta ptree =
    match delta with
    | `Edge (nid, l1, l2, handle) ->
            let ptree = Ptree.remove_node_data nid ptree in
            (match handle with
            | Some v ->
                    (match Tree.get_node v ptree.Ptree.tree with
                    | Tree.Interior (_, neigh, _, _)
                    | Tree.Leaf (_, neigh) ->
                            if debug_joinfn then begin
                                Printf.printf "I am going to do incremental downpass in \
                                break of vertex %d.\n%!" l1;
                                Printf.printf "My current handle is %d\n%!" (Ptree.handle_of
                                l1 ptree);
                            end;
                            let nt, todo = incremental_downpass l1 ptree in
                            if not doup then nt, todo 
                            else incremental_uppass nt todo, []
                    | Tree.Single _ -> ptree, [])
            | None ->
                    if debug_joinfn then begin
                        Printf.printf "I am going to do incremental downpass in \
                        break of vertex %d.\n%!" l1;
                        Printf.printf "My current handle is %d\n%!" (Ptree.handle_of
                        l1 ptree);
                    end;
                    let nt, updlt = incremental_downpass l1 ptree in
                    if not doup then nt, updlt
                    else incremental_uppass nt updlt, [])
    | `Single (nid, newhandle) -> 
            assert (Tree.is_handle nid ptree.Ptree.tree);
            let data = Ptree.get_node_data nid ptree in
            let cost = Node.Standard.tree_cost None data in
            let data = Node.Standard.fix_preliminary data in
            let item = Some ((`Single nid), data) in
            Ptree.assign_root_to_connected_component nid item cost None ptree, []

let update_tree_data_break (tree_delta, clade_delta) ptree =
    let ptree, _ = update_tree_data_break true tree_delta ptree in
    update_tree_data_break false clade_delta ptree

let reroot_fn (_ : bool) edge ptree =
    let Tree.Edge (h, n) = edge in
    let my_handle = Ptree.handle_of h ptree in
    let ptree = Ptree.remove_root_of_component my_handle ptree in
    let ptree, path = Ptree.move_handle h ptree in
    let ptree = Ptree.fix_handle_neighbor h n ptree in
    incremental_downpass_path_to_handle [] None ptree path

(* break_fn must have type handle * node -> tree -> tree * delta * aux_data *)
let break_fn (tree_node, clade_node_id) ptree =
    let Tree.Edge (tree_node, clade_node_id) = 
        Tree.normalize_edge (Tree.Edge (tree_node, clade_node_id)) 
        ptree.Ptree.tree 
    in
    let jxn = (tree_node, clade_node_id) in
    let prev_cost = Ptree.get_cost `Adjusted ptree in
    (* Figure out the clade root node *)
    let clade_node = 
        (Ptree.get_node_data (Tree.int_of_id clade_node_id) ptree)
    in
    (* Get rid of old uppass data *)
    let clade_node = Node.Standard.fix_preliminary clade_node in
    (* Print the tree if asked *)
    if debug_print_on_break then print_node_data_indent ptree;
    (* Perform the break *)
    let nbt, tree_delta = Tree.break jxn ptree.Ptree.tree in
    let tree_handle, clade_handle =
        Tree.get_break_handles tree_delta nbt 
    in
    let ptree = 
        ptree 
            --> Ptree.remove_root_of_component tree_node
            --> Ptree.remove_root_of_component clade_node_id
    in
    let ptree = { ptree with Ptree.tree = nbt } in
    let ptree = set_clade_root ptree clade_node clade_handle in
    let ptree, udlt = update_tree_data_break tree_delta ptree in
    let new_cost = Ptree.get_cost `Adjusted ptree in
    (* Compare costs *)
    let b_delta =
        if prev_cost = infinity && new_cost = infinity
        then 0.
        else prev_cost -. (new_cost -. ptree.Ptree.origin_cost)
    in
    if debug_breakfn_deltas then 
        odebug ("old cost " ^ string_of_float prev_cost
                 ^ " and new cost " ^ string_of_float new_cost);
    if debug_joinfn then begin
        let tc = Ptree.get_cost `Adjusted ptree in
        let nt = force_downpass ptree in
        let nc = Ptree.get_cost `Adjusted nt in
        if nc <> tc then begin
            Printf.printf "The calculated cost is %f but the real cost is \
            %f.\n%!" tc nc;
            failwith "Break update is failing";
        end;
    end;
    let left, right =
        let component_root x =
            let cr = Ptree.get_component_root x ptree in
            match cr.Ptree.root_median with
            | Some (_, b) -> b
            | None -> assert false
        in
        let extract_side x side =
            { Ptree.clade_id = x; 
            clade_node = component_root x;
            topology_delta = side;}
        in
        let (left, right) = tree_delta in
        extract_side tree_handle left, extract_side clade_handle right
    in
    {
        Ptree.ptree = ptree; 
        tree_delta = tree_delta;
        break_delta = b_delta;
        left = left;
        right = right;
        incremental = udlt;
    }

let prepare_tree_for_downpass ptree tree_delta =
    match tree_delta with
    | (`Edge (v, _, _, _)), (`Single (h, true)), _ ->
            if debug_joinfn then Printf.printf "Removing the vertex %d\n%!" v;
            let ptree = Ptree.remove_node_data v ptree in
            v, Ptree.remove_root_of_component h ptree
    | (`Single (v, _)), (`Single (h, true)), _ ->
            v, Ptree.remove_root_of_component h ptree
    | (`Edge (v, _, _, _)), (`Edge (r, a, b, Some h)), _ ->
            if debug_joinfn then 
                Printf.printf "Removing the vertices %d and %d\n%!" v r;
            (let ptree = Ptree.remove_root_of_component h ptree in
            let ptree = Ptree.remove_node_data r ptree in
            let ptree = Ptree.remove_node_data v ptree in
            r, ptree)
    | (`Single (v, _)), (`Edge (r, a, b, Some h)), _ ->
            if debug_joinfn then Printf.printf "Removing the vertex %d\n%!" r;
            let ptree = Ptree.remove_root_of_component h ptree in
            let ptree = Ptree.remove_node_data r ptree in
            r, ptree
    | _ -> failwith "Unexpected Chartree.join_fn"

let join_topologies_and_data jxn1 jxn2 ptree = 
    let ptree = 
        match jxn2 with
        | Tree.Single_Jxn _ -> ptree
        | Tree.Edge_Jxn (a, b) -> 
                let do_reroot () = 
                    let ptree, incr = 
                        reroot_fn true (Tree.Edge (a, b)) ptree in
                    incremental_uppass ptree incr
                in
                let handle = Ptree.handle_of a ptree in
                let root = Ptree.get_component_root handle ptree in
                match root.Ptree.root_median with
                | Some ((`Single _), _)
                | None -> do_reroot ()
                | Some ((`Edge (x, y)), _) when 
                    (x <> a || y <> b) && (x <> b || y <> a) -> do_reroot ()
                | _ -> ptree
    in
    let ret, tree_delta = Tree.join jxn1 jxn2 ptree.Ptree.tree in
    let ptree = { ptree with Ptree.tree = ret } in
    let v, ptree = prepare_tree_for_downpass ptree tree_delta in
    let ptree, updt = 
        if debug_joinfn then 
            Printf.printf "I will do the downpass\n%!";
        incremental_downpass v ptree
    in
    ptree, v, tree_delta, updt

(* join_fn must have type join_1_jxn -> join_2_jxn -> delta -> tree -> tree *)
let join_fn incremental jxn1 jxn2 ptree =
    let ptree = incremental_uppass ptree incremental in
    let ptree, v, tree_delta, updt = join_topologies_and_data jxn1 jxn2 ptree in
    let ptree = incremental_uppass (force_downpass ptree) updt in
    if debug_joinfn then begin
        Printf.printf "Processing join\n%!";
        let current_cost = Ptree.get_cost `Adjusted ptree in
        Printf.printf "Forcing a downpass.\n%!";
        let tmp_tree = force_downpass ptree in
        let new_cost = Ptree.get_cost `Adjusted tmp_tree in
        if new_cost <> current_cost then 
            (Tree.print_join_1_jxn jxn1;
            Tree.print_join_2_jxn jxn2;
            Printf.printf "The vertex for the downpass is %d and has as \
            parental the vertex %d\n. The handle of v is %d\n. The calculated cost \
            is %f but the real cost is %f\n" 
            v (Ptree.get_parent v ptree) (Ptree.handle_of v ptree) current_cost
            new_cost;
            Printf.printf "The handles I have recorded are: \n%!";
            All_sets.IntegerMap.iter (fun code _ -> 
                Printf.printf "%d - " code) ptree.Ptree.component_root;
            print_newline ();
            );
        Printf.printf "End of join\n\n\n%!";
        end;
    ptree, tree_delta

let cost_fn jxn1 jxn2 delta clade_data tree =
    if debug_costfn
    then begin
        (* check clade_data against jxn2 *)
        match jxn2 with
        | Tree.Single_Jxn h ->
              assert (0 = Node.compare_downpass clade_data
                      (Ptree.get_node_data h tree))
        | Tree.Edge_Jxn (h, n) ->
              assert (0 = Node.compare_downpass clade_data
                      (Node.Standard.median None None
                           (Ptree.get_node_data h tree)
                           (Ptree.get_node_data n tree)))
    end;
    match jxn1 with
    | Tree.Single_Jxn h ->
          Ptree.Cost
              (Node.Standard.distance 0.
                   (Ptree.get_node_data (Tree.int_of_id h) tree)
                   clade_data)
    | Tree.Edge_Jxn (h, n) ->
          let h, n =
              if Ptree.is_edge (Tree.Edge (h, n)) tree
              then h, n
              else n, h in
          let ndata = Ptree.get_node_data n tree in

         (* We must check whether our parent is the handle, and, if we are
            the handle's "parent," then our real parent is the virtual
            root, which we calculate by taking a median. *)
         (* Note that [h] is guaranteed to be closer to the handle than
            [n]. *)
          let hdata =
              Ptree.get_parent_or_root_data n tree in

          Ptree.Cost (Node.dist_2 delta clade_data hdata ndata)

let within_threshold b a = b < (a +. (a *. 0.05))

let exact_cost_fn jxn1 jxn2 delta clade_data ptree =
    match cost_fn jxn1 jxn2 delta clade_data ptree with
    | (Ptree.Cost v) as res ->
            if within_threshold v delta then
                let previous_cost = Ptree.get_cost `Adjusted ptree in
                let ptree, _, _, _ = join_topologies_and_data jxn1 jxn2 ptree in
                Ptree.Cost ((Ptree.get_cost `Adjusted ptree) -. previous_cost)
            else res
    | x -> x

let incremental_downpass_path ptree path =
    incremental_downpass_path_to_handle [] None ptree path

let rec features meth lst =
    let `LocalOptimum (meth) = meth in
    
    let hood =
        let t = function
            | `Spr -> "spr"
            | `Tbr -> "tbr" in
        match meth.Methods.ss with
        | `SingleNeighborhood x -> "single-" ^ t x
        | `ChainNeighborhoods x -> t x
        | `Alternate (x, y) -> "alternate-" ^ t x ^ "-" ^ t y
        | `None -> "none" in
    let lst = ("type", hood) :: lst in
    let lst = ("origin-cost", match meth.Methods.oo with None -> "none"
               | Some cost -> string_of_float cost) :: lst in
    let lst =
        match meth.Methods.tm with
        | `BestFirst -> ("trajectory", "best-first") :: lst
        | `AllAround _ -> ("trajectory", "all-around") :: lst
        | `AllThenChoose -> ("trajectory", "All around then choose") :: lst
        | `PoyDrifting (equalprob, worstfactor) ->
              ("trajectory", "poy 3.0 drifting")
              :: ("equal-probability", string_of_float equalprob)
              :: ("worst-factor", string_of_float worstfactor)
              :: lst 
        | `Annealing (coeff, exp) ->
              ("trajectory", "simulated annealing")
              :: ("annealing-coefficient", string_of_float coeff)
              :: ("annealing-exponent", string_of_float exp)
              :: lst in
    let lst = 
        let lst = 
            match meth.Methods.tabu_break with
            | `Randomized -> ("tabu.break", "randomized") :: lst
            | `DistanceSorted true -> 
                    ("tabu.break", "sorted edges by distance") :: lst
            | `DistanceSorted false -> 
                    ("tabu.break", "sorted edges by distance with early stop") :: lst
            | `OnlyOnce ->
                    ("tabu.break", "only break once each edge, never again") :: lst
        in
        let add_depth depth lst = 
            (match depth with
            | Some v -> ("tabu.join.depth", string_of_int v) :: lst
            | None -> lst)
        in
        let lst = 
            match meth.Methods.tabu_join with
            | `UnionBased depth ->
                    ("tabu.join", "Union of children states prunes the tree") ::
                        add_depth depth lst
            | `AllBased depth ->
                    ("tabu.join", "All edges are joined") ::
                        add_depth depth lst
            | `Partition options ->
                    List.fold_left (fun acc x ->
                        match x with
                        | `Sets sets ->
                                ("tbu.join", 
                                "Programmatically provided constraint") :: acc
                        | `MaxDepth depth ->
                                ("tabu.join", "Consensus based partition") ::
                                    add_depth (Some depth) acc
                        | `ConstraintFile file ->
                                ("tabu.join", "File based partition") :: 
                                ("tabu.join.constraint", (FileStream.filename
                                file)) :: 
                                    acc) 
                    lst options
        in
        match meth.Methods.tabu_reroot with
        | `Bfs depth ->
                ("tabu.reroot", "Breath first search") ::
                    (match depth with
                    | Some v -> ("tabu.reroot.depth", string_of_int v) :: lst
                    | None -> lst)
    in
    lst

let rec get_active_ref_code_node tree node_id = 
    let node_data = Ptree.get_node_data node_id tree in

    match Ptree.get_node node_id tree with 
    | Tree.Interior (_, parent, child1_id, child2_id) ->

          let acc_pre_child1, acc_fi_child1 = get_active_ref_code_node tree child1_id in 
          let acc_pre_child2, acc_fi_child2 = get_active_ref_code_node tree child2_id in 
          let _, pre_child, _, fi_child = Node.get_active_ref_code node_data in 

          IntSet.union (IntSet.union acc_pre_child1 acc_pre_child2) pre_child, 
          IntSet.union (IntSet.union acc_fi_child1 acc_fi_child2) fi_child 

    | Tree.Leaf (_, parent) ->
          let _, _, fi, fi_child = Node.get_active_ref_code node_data in 
          IntSet.empty, (IntSet.union fi fi_child)

    | Tree.Single _ ->
          let pre, _, fi, _ = Node.get_active_ref_code node_data in 
          pre, fi

let get_active_ref_code tree = 
    let handles = All_sets.Integers.elements tree.Ptree.tree.Tree.handles in
    let pre_act, fi_act = List.fold_left
        (fun (acc_pre_act, acc_fi_act) handle_id ->
             let pre, fi =
                 match Ptree.get_node handle_id tree with
                 | Tree.Interior (_, parent, _, _)
                 | Tree.Leaf (_, parent) ->
                       let pre, fi = get_active_ref_code_node tree handle_id in
                       let pre2, fi2 = get_active_ref_code_node tree parent in
                       IntSet.union pre pre2, IntSet.union fi fi2
                 | Tree.Single _ -> IntSet.empty, IntSet.empty
             in
             let fi = IntSet.filter (fun fi_code -> not (IntSet.mem fi_code pre)) fi
             in
             let root = Ptree.get_component_root handle_id tree in
             match root.Ptree.root_median with  
             | None -> pre, fi 
             | Some (_, root_data) ->  
                   let r_pre, r_pre_child, r_fi, _ = 
                       Node.get_active_ref_code root_data 
                   in
                   IntSet.union (IntSet.union pre r_pre_child) r_pre,
                   IntSet.union fi r_fi) (IntSet.empty, IntSet.empty) handles
    in 
    pre_act, fi_act


(** Note that parent_data = root_data in case the 
    anaylzing node is the handle or handle's ancestor *)
let rec subtree_to_formatter diag_report_type (pre_ref_codes, fi_ref_codes) 
        attr tree node_id 
        (parent_data : (Node.node_data * Node.node_data) option) : Xml.xml =
    match Ptree.get_node node_id tree with
    | Tree.Interior (_, parent_id, c1, c2)  ->
            let node_data = Ptree.get_node_data node_id tree
            and child1_node_data = Ptree.get_node_data c1 tree
            and child2_node_data = Ptree.get_node_data c2 tree
            and parent_node_data, single_parent  = Utl.deref parent_data 
            in 
            let my_single_assignment =
                Node.to_single (pre_ref_codes, fi_ref_codes) false
                    None single_parent node_data 
            in
            (* We recursively call for the current vertex, then the two children
            * *)
            let nodest =
                try
                    let node_formatter = 
                        Node.to_formatter_subtree diag_report_type (pre_ref_codes, fi_ref_codes) [] 
                        tree.Ptree.data (node_data, node_data (* TODO: possible bug, but
                        this feature is not used  in this branch of the code *) )
                        node_id (c1, child1_node_data) (c2, child2_node_data)
                        (Some (parent_node_data, single_parent))
                    in
                    (`Single node_formatter)
                with Not_found -> `Empty
            in 
            let child1_formatter = 
                subtree_to_formatter diag_report_type
                (pre_ref_codes, fi_ref_codes) 
                [] tree c1 (Some (node_data, my_single_assignment))
            in
            let child2_formatter = 
                subtree_to_formatter diag_report_type
                (pre_ref_codes, fi_ref_codes)
                [] tree c2 (Some (node_data, my_single_assignment))
            in
            let c1st = `Single child1_formatter and c2st = `Single child2_formatter in
            (Xml.Trees.tree, [], (`Set [nodest; c1st; c2st]))
    | Tree.Leaf (_, parent_id) -> 
            let node_data = Ptree.get_node_data node_id tree in
            assert (match parent_data with 
                | Some _ -> true
                | None -> false);
            let nodest =
                try
                    let node_formatter = 
                        Node.to_formatter_single
                        diag_report_type (pre_ref_codes, fi_ref_codes) []  
                        tree.Ptree.data 
                        (node_data, node_data (* TODO Bug *)) node_id parent_data
                    in
                    `Single node_formatter
                with Not_found -> `Empty
            in 
          (Xml.Trees.tree, [], nodest)
    | Tree.Single _ ->
          let nodest = 
              let node_data = Ptree.get_node_data node_id tree in 
              let node_formatter = 
                  Node.to_formatter_single diag_report_type (pre_ref_codes, fi_ref_codes) 
                  [] tree.Ptree.data (node_data, node_data (* TODO: Bug *)) node_id None 
              in
              `Single node_formatter
          in 
          (Xml.Trees.tree, [], nodest)

let handle_to_formatter diag_report_type (pre_ref_codes, fi_ref_codes)
        attr tree handle_id : Xml.xml =
    let root = Ptree.get_component_root handle_id tree in
    let data = 
        match Ptree.get_node handle_id tree with
        | Tree.Interior (_, parent, _, _)
        | Tree.Leaf (_, parent) ->
              let root_formatter, tree_root, handle_node_data, parent_node_data =
                  match root.Ptree.root_median with 
                  | Some ((`Edge (handle_id, parent)), root) -> 
                          let handle_node_data = Ptree.get_node_data handle_id tree
                          and parent_node_data = Ptree.get_node_data parent tree in
                          let root_single = 
                              Node.to_single  (pre_ref_codes, fi_ref_codes) false
                                  (Some root) parent_node_data handle_node_data
                          in
                          let root_f = 
                              Node.to_formatter_subtree diag_report_type (pre_ref_codes, fi_ref_codes) 
                                  [] tree.Ptree.data  (root, root_single (* TODO: Bug *))  handle_id (handle_id,  handle_node_data) 
                                  (parent, parent_node_data) None
                          in 
                          root_f, Some (root_single, root_single),
                          handle_node_data, parent_node_data
                  | _ -> failwith "How is it possible we have no root?"
              in               
              let handle_f = 
                  subtree_to_formatter diag_report_type (pre_ref_codes, fi_ref_codes) 
                  [] tree handle_id tree_root 
              in 
              let parent_f = 
                  subtree_to_formatter diag_report_type (pre_ref_codes, fi_ref_codes) 
                  [] tree parent tree_root
              in 
              let c1 = `Single handle_f and c2 = `Single parent_f in 
              let root = `Single root_formatter in 
              (`Set [root; c1; c2])
        | Tree.Single _ ->
                let node_single = 
                    let nd = Ptree.get_node_data handle_id tree in
                    Some (nd, Node.to_single (pre_ref_codes, fi_ref_codes) false None nd nd)
                in
                let c1 = 
                    subtree_to_formatter diag_report_type (pre_ref_codes, fi_ref_codes) 
                    [] tree handle_id node_single 
                in
                (`Single c1)
    in
    let attr = 
        (Xml.Trees.cost, `Float root.Ptree.component_cost) :: attr
    in
    (Xml.Trees.tree, attr, data)

let to_formatter diag_report_type atr tree : Xml.xml =
        (* We don't include the cost of the tree because it comes from the three
        * directional tree attributes. *)
    let (pre_ref_codes, fi_ref_codes) = get_active_ref_code tree in
    let tag = Xml.Trees.forest in
    let handles = All_sets.Integers.elements tree.Ptree.tree.Tree.handles in
    let root_recost = ref 0. in 
    let data = 
        List.map 
        (fun handle_id -> 
             let handle_formatter = handle_to_formatter diag_report_type (pre_ref_codes, fi_ref_codes)
                 [] tree handle_id in 
             let root = Ptree.get_component_root handle_id tree in 
             (match root.Ptree.root_median with  
             | None -> ()
             | Some (_, root_data) ->  
                   let recost = Node.cmp_subtree_recost root_data in
                   root_recost := !root_recost +. recost);
             `Single handle_formatter) handles 
    in
    let atr = (Xml.Trees.recost, `Float !root_recost) :: atr in
    (tag, atr, (`Set data))


let root_costs tree = 
    All_sets.Integers.fold (fun handle acc ->
        let get_cost edge acc = 
            let ntree, _ = reroot_fn false edge tree in
            (Tree.Continue, (((edge, Ptree.get_cost `Adjusted ntree) :: acc)))
        in
        Ptree.pre_order_edge_visit get_cost handle tree acc)
    (Tree.get_handles tree.Ptree.tree) []

let total_cost ptree adj chars =
    let total_handle_cost acc h =
        match (Ptree.get_component_root h ptree).Ptree.root_median with
            | Some ((`Edge (a,b)),c) -> 
                acc +. (Node.total_cost None c)
            | None
            | Some _ -> acc
    in
    All_sets.Integers.fold
        (fun handle acc_cost -> total_handle_cost acc_cost handle)
        (Ptree.get_handles ptree)
        0.0

let refresh_all_edges a b c d = d

module TreeOps = struct
    type a = Node.node_data
    type b = Node.node_data
    let clear_internals _ ptree = ptree
    let break_fn _ = break_fn
    let join_fn _ = join_fn
    let model_fn x = x
    let adjust_fn ?(max_iter) _ t = t
    let cost_fn _ = 
        match !Methods.cost with
        | `Normal -> cost_fn
        | _ -> exact_cost_fn
    let reroot_fn _ = reroot_fn
    let string_of_node = Node.Standard.to_string
    let features meth lst = features meth (("skipping", "false") :: lst)
    let downpass = downpass
    let uppass = uppass
    let incremental_uppass = incremental_uppass   
    let to_formatter = to_formatter
    let get_active_ref_code = get_active_ref_code
    let branch_table _ = Hashtbl.create 1
    let root_costs = root_costs
    let unadjust ptree = ptree
    
    let refresh_all_edges = refresh_all_edges
    (* not implemented in this module *)
    let tree_size _ _ = 0.0
    let prior_cost _ _ = 0.0
    let total_cost = total_cost
end
