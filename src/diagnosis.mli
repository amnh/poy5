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

module type S = sig
    type a 
    type b
(*    val assign_ref : Data.d -> int -> int -> (a, b) Ptree.p_tree Sexpr.t -> 
        (a, b) Ptree.p_tree Sexpr.t * int * int
*)
    val diagnosis :
      Data.d -> (a, b) Ptree.p_tree Sexpr.t -> Methods.diagnosis -> unit
end

module Make 
    (Node : NodeSig.S) (Edge : Edge.EdgeSig with type n = Node.n)
    (TreeOps : 
        Ptree.Tree_Operations with type a = Node.n with type b = Edge.e) 
    : S with type a = Node.n with type b = Edge.e
