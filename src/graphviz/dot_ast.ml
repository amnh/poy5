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

(* $Id:$ *)

(** AST for DOT file format. *)

type id = 
  | Ident of string
  | Number of string
  | String of string
  | Html of string

type attr = (id * id option) list

type compass_pt = N | Ne | E | Se | S | Sw | W | Nw

type port = 
  | PortId of id * compass_pt option
  | PortC of compass_pt

type node_id = id * port option

type subgraph = 
  | SubgraphId of id
  | SubgraphDef of id option * stmt list

and node =
  | NodeId of node_id
  | NodeSub of subgraph

and stmt = 
  | Node_stmt of node_id * attr list
  | Edge_stmt of node * node list * attr list
  | Attr_graph of attr list
  | Attr_node of attr list
  | Attr_edge of attr list
  | Equal of id * id
  | Subgraph of subgraph

type file =
  { strict : bool;
    digraph : bool;
    id : id option;
    stmts : stmt list; }

let string_of_id = function
    | Ident s
    | Number s
    | String s
    | Html s -> s

let compare a b = compare (string_of_id a) (string_of_id b)

let find_subgraph name g =
    let rec find_ x = match x with
        | (Equal _ ) ::xs
        | (Attr_edge _ ) :: xs
        | (Attr_node _) :: xs
        | (Attr_graph _):: xs
        | (Subgraph (SubgraphId _)) :: xs
        | (Node_stmt _) :: xs
        | (Subgraph (SubgraphDef (None,_))) :: xs
        | (Edge_stmt _) :: xs -> find_ xs
        | (Subgraph ((SubgraphDef (Some sub,graph)) as def)) :: xs ->
            if 0 = compare name sub then def else find_ graph
        | [] -> raise Not_found
    in
    find_ g.stmts

module OrderedId = struct
    type t = id
    let compare a b = compare a b
end

module IdSet = Set.Make (OrderedId)
module IdMap = Map.Make (OrderedId)
