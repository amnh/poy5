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

exception IllegalEdgeConversion

module type EdgeSig = sig
    (** The type of the edge data *)
    type e
    type n
    (** Weather or not an edge holds data (false when e = unit), or the edge
    * type is not to be used *)
    val has_information : bool
    (*
    (** calculates the edge defined by the two nodes *)
    val calculate : n -> n -> e
    *)
    (* Convert the contents of an edge into a node (if it is possible). If it is
    * not possible due to the nature of the edge information, raise
    * IllegalEdgeConversion *)
    val to_node : int -> (int * int) -> e -> n
    val of_node : int option -> n -> e
end

module Edge : EdgeSig with type e = unit with type n = Node.node_data
module SelfEdge : EdgeSig with type e = Node.node_data with type n = Node.node_data

module LazyEdge : EdgeSig with type e = AllDirNode.OneDirF.n with type n =
    AllDirNode.AllDirF.n 

(*
module SuperRoot (Node : NodeSig.S) : 
    EdgeSig with type e = Node.n with type n = Node.n

module LazyRoot (Node : NodeSig.S) : 
    EdgeSig with type e = Node.n Lazy.t with type n = Node.n
*)
