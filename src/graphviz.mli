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

(** For parsing DOT files --courtesy of ocamlgraph *)
module Dot_ast : sig

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
        stmts : stmt list }

    val string_of_id : id -> string

    val find_subgraph : id -> file -> subgraph

    module IdSet : Set.S with type elt = id
    module IdMap : Map.S with type key = id

end

module Dot : sig

    (** a basic representation for ease of use *)
    type basic =  
        { name  : string option;
          nodes : Dot_ast.IdSet.t;
          children : Dot_ast.IdSet.t Dot_ast.IdMap.t;
          parents  : Dot_ast.IdSet.t Dot_ast.IdMap.t;
          equal : (Dot_ast.id * Dot_ast.id) list; }

    (** parse a Dot_ast from a channel *)
    val of_channel : in_channel -> Dot_ast.file

    (** parse a Dot_ast from a file *)
    val of_file : string -> Dot_ast.file

    (** Convert a Dot_ast to a Basic map of sets which define directed edges. *)
    val to_basic : Dot_ast.file -> basic

    (** helper function to return the leaf set *)
    val leaves : basic -> Dot_ast.IdSet.t

    (** return a set of names that are equivlent *)
    val equivalences : Dot_ast.id -> basic -> Dot_ast.IdSet.t
    
    (** convert the basic to a Dot_ast parsable string *)
    val basic_to_string : ?node_attr:string -> ?edge_attr:string -> basic -> string
end
