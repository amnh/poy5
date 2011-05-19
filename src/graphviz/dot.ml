(**************************************************************************)
(*                                                                        *)
(*  Ocamlgraph: a generic graph library for OCaml                         *)
(*  Copyright (C) 2004-2010                                               *)
(*  Sylvain Conchon, Jean-Christophe Filliatre and Julien Signoles        *)
(*                                                                        *)
(*  This software is free software; you can redistribute it and/or        *)
(*  modify it under the terms of the GNU Library General Public           *)
(*  License version 2.1, with the special exception on linking            *)
(*  described in file LICENSE.                                            *)
(*                                                                        *)
(*  This software is distributed in the hope that it will be useful,      *)
(*  but WITHOUT ANY WARRANTY; without even the implied warranty of        *)
(*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.                  *)
(*                                                                        *)
(**************************************************************************)

(** Parser for DOT file formats and a stripped down representation for our use *)

(** {6 Types *)

(** A striped down subset of details from Dot_ast for easier usage. We keep the
    equality list and not normlize the values in the map since it is unknown which
    is necessary to use as a name. *)
type basic =  
    { name  : string option;
      nodes : Dot_ast.IdSet.t;
      children : Dot_ast.IdSet.t Dot_ast.IdMap.t;
      parents  : Dot_ast.IdSet.t Dot_ast.IdMap.t;
      equal : (Dot_ast.id * Dot_ast.id) list;
    }

(** empty basic information *)
let empty = 
    {   name  = None;
        nodes = Dot_ast.IdSet.empty;
        children = Dot_ast.IdMap.empty;
        parents  = Dot_ast.IdMap.empty;
        equal = []; 
    }

(** parse a channel of the Dot Language *)
let of_channel c =
    let lb = Lexing.from_channel c in
    let dot =
        try Dot_parser.file Dot_lexer.token lb
        with Parsing.Parse_error ->
            let n = Lexing.lexeme_start lb in
            failwith (Printf.sprintf "Dot.parse: parse error character %d" n)
    in
    dot

(** parse a file of the Dot Language *)
let of_file f =
    let chan = open_in f in
    let res = of_channel chan in
    close_in chan;
    res

(** find all the name equivalences to a name *)
let equivalences name basic =
    let rec follow_left acc name =
        let neigh = 
            List.find_all (fun (a,_) -> 0 = Dot_ast.compare a name) basic.equal
        in
        List.fold_left
            (fun acc (a,b) ->
                if Dot_ast.IdSet.mem b acc then acc
                else follow_left (Dot_ast.IdSet.add b acc) b)
            acc
            neigh
    and follow_right acc name =
        let neigh = 
            List.find_all (fun (_,b) -> 0 = Dot_ast.compare b name) basic.equal
        in
        List.fold_left
            (fun acc (a,b) ->
                if Dot_ast.IdSet.mem a acc then acc
                else follow_right (Dot_ast.IdSet.add a acc) a)
            acc
            neigh
    in
    follow_right (follow_left (Dot_ast.IdSet.singleton name) name) name

(** Convert Dot_ast to a basic stripped down type mentioned above *)
let to_basic (g : Dot_ast.file) : basic =
    assert( g.Dot_ast.digraph ); (* only parse di-graphs *)
    (* find all the equality pairs *)
    let rec get_equals acc = function
        | [] -> 
            acc
        | (Dot_ast.Equal (data1,data2)) :: oth -> 
            get_equals ((data1,data2) :: acc) oth
        | (Dot_ast.Node_stmt _) :: oth
        | (Dot_ast.Subgraph (Dot_ast.SubgraphId _)) :: oth
        | (Dot_ast.Attr_edge _) :: oth
        | (Dot_ast.Attr_node _) :: oth
        | (Dot_ast.Attr_graph _):: oth -> get_equals acc oth
        | (Dot_ast.Subgraph (Dot_ast.SubgraphDef (_,data))) :: oth ->
            get_equals (get_equals acc data) oth
        | (Dot_ast.Edge_stmt (name,names,_)) :: oth ->
            let acc =
                begin match name with
                    | Dot_ast.NodeId (name,_) -> acc
                    | Dot_ast.NodeSub sub -> get_equals acc [Dot_ast.Subgraph sub]
                end
            in
            get_equals acc oth
    in
    (* since we are going to need the equality for functions below... *)
    let almost =
        { empty with
            equal = get_equals [] g.Dot_ast.stmts;
            name  = match g.Dot_ast.id with
                    | Some x -> Some (Dot_ast.string_of_id x);
                    | None   -> None;
        }
    in
    (* add an edge to the map *)
    let rec add_edge a b (nodes,(cedges,pedges)) =
        let a = normalize_name a in
        let b = normalize_name b in
        let cset =
            if Dot_ast.IdMap.mem a cedges
                then Dot_ast.IdSet.add b (Dot_ast.IdMap.find a cedges)
                else Dot_ast.IdSet.singleton b
        and pset =
            if Dot_ast.IdMap.mem b pedges
                then Dot_ast.IdSet.add a (Dot_ast.IdMap.find b pedges)
                else Dot_ast.IdSet.singleton a
        in
        let nodes = add_node a (add_node b nodes) in
        nodes, (Dot_ast.IdMap.add a cset cedges,
                Dot_ast.IdMap.add b pset pedges)
    (* add a node name to a map *)
    and add_node name nodes =
        if Dot_ast.IdSet.mem name nodes then nodes
        else Dot_ast.IdSet.add name nodes
    (* normalize a name to the FIRST argument (this is so we can easily use the
       assoc functions in list later *)
    and normalize_name name = (* to the minimum of equaivlences *)
        Dot_ast.IdSet.min_elt (equivalences name almost)
    in 
    (* Process a single element in the graph; skip over attributes, recurse
       through subgraph, add node names when they come across, and process edges
       in the peculiar manor that they need to be processed. *)
    let rec process_nodes ((nodes,edges) as acc) = function
        | [] -> acc
        | (Dot_ast.Subgraph (Dot_ast.SubgraphId name)) :: oth ->
            let resolved = Dot_ast.Subgraph (Dot_ast.find_subgraph name g) in
            process_nodes acc (resolved::oth)
        | (Dot_ast.Subgraph (Dot_ast.SubgraphDef (_,data))) :: oth ->
            process_nodes (process_nodes acc data) oth
        | (Dot_ast.Node_stmt ((name,_),_)) :: oth ->
            process_nodes (add_node name nodes,edges) oth
        | (Dot_ast.Edge_stmt (name,names,_)) :: oth ->
            let parents,edges =
                begin match name with
                    | Dot_ast.NodeId (name,_) ->
                        Dot_ast.IdSet.singleton name,edges
                    | Dot_ast.NodeSub sub ->
                        process_nodes (Dot_ast.IdSet.empty,edges) [Dot_ast.Subgraph sub]
                end
            in
            let nodes = Dot_ast.IdSet.union parents nodes in
            process_nodes (process_edges parents names (nodes,edges)) oth
        | (Dot_ast.Equal _) :: oth
        | (Dot_ast.Attr_edge _) :: oth
        | (Dot_ast.Attr_node _) :: oth
        | (Dot_ast.Attr_graph _):: oth -> process_nodes acc oth
    (* process an edge starting from [ps] to the path [ns]; each element in ns
       could be a single node or a subgraph. This adds the interesting property
       that the parent could be a set of nodes as in, 4 -> {2 -> 3, 1} -> 5. In
       this case, all nodes contained in the subgraph are directed to the next *)
    and process_edges ps ns ((nodes,edges) as acc) = match ns with
        | [] -> acc
        | (Dot_ast.NodeId (name,_)) :: xs ->
            let acc = Dot_ast.IdSet.fold (fun par iacc -> add_edge par name iacc) ps acc in
            process_edges (Dot_ast.IdSet.singleton name) xs acc
        | (Dot_ast.NodeSub (Dot_ast.SubgraphId name)) :: xs ->
            let resolved = Dot_ast.NodeSub (Dot_ast.find_subgraph name g) in
            process_edges ps (resolved::xs) acc
        | (Dot_ast.NodeSub (Dot_ast.SubgraphDef (_,stmts))) :: xs ->
            let subnodes,edges = process_nodes (Dot_ast.IdSet.empty,edges) stmts in
            let acc = (* attach subgraph nodes to parents *)
                Dot_ast.IdSet.fold (fun par acc ->
                    Dot_ast.IdSet.fold (fun chi iacc -> add_edge par chi iacc)
                                        subnodes
                                        acc)
                    ps
                    (nodes,edges)
            in
            process_edges subnodes xs acc
    in
    let acc = Dot_ast.IdSet.empty,(Dot_ast.IdMap.empty,Dot_ast.IdMap.empty) in
    let nodes,(cedges,pedges) = process_nodes acc g.Dot_ast.stmts in
    { empty with
        nodes = nodes;
        children = cedges;
        parents = pedges;
    }

(** helper function to return the leaf set **)
let leaves base : Dot_ast.IdSet.t =
    Dot_ast.IdSet.fold
        (fun k acc ->
            if Dot_ast.IdMap.mem k base.children then
                match Dot_ast.IdSet.cardinal (Dot_ast.IdMap.find k base.children) with
                | 0 -> Dot_ast.IdSet.add k acc
                | _ -> acc
            else
                Dot_ast.IdSet.add k acc)
        base.nodes
        Dot_ast.IdSet.empty

(** write representation to a Dot_ast readable format in a string**)
let basic_to_string ?(node_attr="") ?(edge_attr="") base =
    let buffer = Buffer.create 1000 in
    let bufferf format =
        Printf.ksprintf (Buffer.add_string buffer) format in
    let () = match base.name with
        | Some x -> bufferf "digraph %s {\n" x
        | None   -> bufferf "digraph {"
    in
    let () =
        Dot_ast.IdSet.iter
            (fun x ->
                bufferf "\t\"%s\" [%s];\n" (Dot_ast.string_of_id x) node_attr)
            base.nodes
    in
    let () =
        Dot_ast.IdMap.iter
            (fun k vs ->
                Dot_ast.IdSet.iter
                    (fun v ->
                        bufferf "\t\"%s\" -> \"%s\" [%s];\n" 
                            (Dot_ast.string_of_id k)
                            (Dot_ast.string_of_id v) edge_attr)
                    vs)
            base.children
    in
    bufferf "}";
    Buffer.contents buffer

