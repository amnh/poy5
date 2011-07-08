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

(** Threading/sequence function for a point-free programming experience *)
let (-->) a b = b a

(** Use continuations to produce a formatted failwith message *)
let failwithf format = Printf.ksprintf (failwith) format

(** Id is used for keys to label nodes *)
type id = int

(** Edge are a pair of IDs; no direction is specified *)
type edge = id * id

(** A phantom type to label a graph checked for cycles and other invariants *)
type checked

(** The graph needs to be checked to be used by diagnosis functions *)
type unchecked

(** A description of all the nodes; the first id in each is the node id for
*     itself, the rest depend on context, but parents come before childen. *)
type node =
    | Single     of id                (** single; included in root set *)
    | Leaf       of id * id           (** standard leaf; one parent *)
    | Root       of id * id * id      (** provides two children *)
    | Interior   of id * id * id * id (** one parent and two children  *)
    | Reticulate of id * id * id * id (** two parents and one child    *)

(** Record to define the entire graph; roots are kept in a seperate structure
    for quick lookup (use assert_roots function to verify these are correct). This
    digraph may have cycles, and phantom types are used to ensure we are do not
    have infinite loops; [checked/unchecked]. **)
type 'a di_graph =
    { name  : string option;
      roots : All_sets.Integers.t;
      graph : node All_sets.IntegerMap.t;
      avail_ids : id list;
      new_ids : id;
    }

(* an empty di_graph to start computation off *)
let empty : checked di_graph = 
    {   name = None;
        roots = All_sets.Integers.empty;
        graph = All_sets.IntegerMap.empty;
        avail_ids = [];
        new_ids = 0;
    }

(** get the next id for the graph and return the updated information for later
    retreval of new graph node ids. *)
let pop_id (g: 'a di_graph) : id * 'a di_graph = 
    match g.avail_ids with
    | []     -> g.new_ids + 1, { g with new_ids = g.new_ids + 1; }
    | hd::tl -> hd, { g with avail_ids = tl; }

(** push ids to a graph; used after breaking edges. *)
let push_id (n: id) (g: 'a di_graph) : 'a di_graph =
    if n = (g.new_ids - 1) then
        { g with new_ids = n; }
    else 
        {g with avail_ids = n :: g.avail_ids; }

(** functions to control phantom types; mark checked *)
external check_graph : unchecked di_graph -> checked di_graph = "%identity"

(** add the phantom type 'unchecked' to di_graph *)
external uncheck_graph : checked di_graph -> unchecked di_graph = "%identity"

(** verify roots checks that all roots are in the roots structure, and no roots
    are missing in the integermap. *)
let verify_roots g : bool =
    All_sets.IntegerMap.fold
        (fun k v a -> match v with
            | Leaf _
            | Interior _
            | Reticulate _ -> (not (All_sets.Integers.mem k g.roots)) && a
            | Single _
            | Root _ -> (All_sets.Integers.mem k g.roots) && a)
        (g.graph)
        (true)

(** Post-Order node visit; fold ovear directed graph with funciton f in a
    post-order fashion; types ensure we do not have ridiculous cycles. *)
let post_order_node_visit ~f (r: id) (acc: 'a) (g: checked di_graph): 'a =
    assert (All_sets.Integers.mem r g.roots);
    let rec eval cur (acc:'a): 'a =
        match All_sets.IntegerMap.find cur g.graph with
        | Root (n,c1,c2)         -> acc --> eval c1 --> eval c2 --> f n
        | Interior (n,p,c1,c2)   -> acc --> eval c1 --> eval c2 --> f n
        | Reticulate (n,p1,p2,c) -> acc --> eval c  --> f n
        | Leaf (n,p)             -> acc --> f n
        | Single (n)             -> acc --> f n
    in
    eval r acc


let post_order_node_visit_once ~f r acc (g: checked di_graph) =
    let augmented_f cur ((oacc,set) as acc) =
        if All_sets.Integers.mem cur set then acc
        else (f cur oacc, All_sets.Integers.add cur set)
    in
    post_order_node_visit ~f:augmented_f (acc,All_sets.Integers.empty) g


let post_order_node_visit_once_short_circuit ~f (r: id) (acc: 'a) (g: checked di_graph): 'a =
    assert (All_sets.Integers.mem r g.roots);
    let f n (facc,set) = f n facc, All_sets.Integers.add n set in
    let rec eval cur ((_,set) as acc) =
        if All_sets.Integers.mem n set then
            acc
        else match All_sets.IntegerMap.find cur g.graph with
            | Root (n,c1,c2)         -> acc --> eval c1 --> eval c2 --> f n
            | Interior (n,p,c1,c2)   -> acc --> eval c1 --> eval c2 --> f n
            | Reticulate (n,p1,p2,c) -> acc --> eval c --> f n
            | Leaf (n,p)             -> acc --> f n
            | Single (n)             -> acc --> f n
    in
    eval r acc


(** Post_order edges; fold over the edges of the tree calling f on each edge. The
    parent is always the first argument in the function call. If no edges in the
    graph exist; we return the accumulator --no exceptions, no failures. *)
let post_order_edge_visit ~f ~rf (r: id) (acc: 'a) (g: checked di_graph): 'a =
    assert (All_sets.Integers.mem r g.roots);
    let rec eval (par: id) (cur: id) (acc:'a): 'a =
        match All_sets.IntegerMap.find cur g.graph with
        | Leaf (n,p)             -> f p n acc
        | Interior (n,p,c1,c2)   -> acc --> eval n c1 --> eval n c2 --> f p n
        | Reticulate (n,p1,p2,c) -> acc --> eval n c --> f p1 n --> f p2 n
        | (Root _ | Single _ )   -> assert false
    and eval_r (cur: id) (acc:'a): 'a =
        match All_sets.IntegerMap.find cur g.graph with
        | Root (n,c1,c2)         -> acc --> eval n c1 --> eval n c2 --> rf n c1 c2
        | Single (n)             -> acc
        | Leaf _ | Interior _ | Reticulate _ -> assert false
    in
    eval_r r acc

open Graphviz
(** Transform Dot_ast to di_graph format; unchecked **)
let process_file (d : Data.d) (g : Dot_ast.file) : unchecked di_graph =
    (* process any node information *)
(*    let names : Dot_ast.id Dot_ast.IdMap.t =*)
    (* process edge information and normalize the names *)
(*    let edges,names : (Dot_ast.id list) Dot_ast.IdMap.t * (Dot_ast.id Dot_ast.IdMap.t) =*)
    (* build associations with names against Data.d *)
(*    let normal_names : id Dot_ast.IdMap.t =*)
    (* find the root node(s); a node with two children; and ensure that they
       aren't apart of the same leaf-set. *)
(*    let roots =*)
    (* traverse from the root out determining if the new node is reticulate, a
       leaf or interior *)
(*    let graph =*)
    (* add missing taxon nodes as Single nodes *)
(*    let roots,graph =*)
    let graph_name = match g.Dot_ast.id with
        | Some (Dot_ast.Ident s) -> Some s
        | Some (Dot_ast.Number s) -> Some s
        | Some (Dot_ast.String s) -> Some s
        | Some (Dot_ast.Html s) -> None
        | None -> None
    in
    { empty with
        name = graph_name;
    }

(** transform a tree to a graph *)
let process_tree (tree: Tree.u_tree) : checked di_graph =
    failwith "not_done"
