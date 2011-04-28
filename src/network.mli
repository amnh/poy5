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

(** Module that defines directed graphs for phylogenetic networks *)

(** {6 Types used in this module to define Networks *)

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
    for quick lookup This digraph may have cycles, and phantom types are used to
    ensure we do not have any issues with our invariants. **)
type 'a di_graph =
    { name  : string option;
      roots : All_sets.Integers.t;
      graph : node All_sets.IntegerMap.t;
      avail_ids : id list;
      new_ids : id;
    }


(** {6 Basic Graph Functions / Values **)

(* an empty di_graph to start computation off *)
val empty : checked di_graph


(** {6 Traversal Functions **)

(** post order on nodes; the only argument for the function is the current node *)
val post_order_node_visit : f:(id -> 'a -> 'a) -> id -> 'a -> checked di_graph -> 'a

(** post order on all the edges; call rf on the root (pass its id and its two
    children), and f on the subsequent edges (parent and child). No edges are
    repeated; thus f will never be called with the first argument a leaf *)
val post_order_edge_visit : f:(id -> id -> 'a -> 'a) -> rf:(id -> id -> id -> 'a -> 'a) -> id -> 'a -> checked di_graph -> 'a


(** {6 IO / transformation Functions *)

(** Transform Dot_ast to di_graph format; unchecked **)
val process_file : Data.d -> Graphviz.Dot_ast.file -> unchecked di_graph

(** convert a Tree.u_tree to a di_graph; preserves leaf ids *)
val process_tree : Tree.u_tree -> checked di_graph

