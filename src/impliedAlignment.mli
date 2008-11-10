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

(** A module to build implied alignments of sequences over trees. *)

(** The implied alignments are multiple sequence alignments that follow a
* homology hypothesis implied by the best hypothesis available to explain their
* origin. *)

exception NotASequence of int

(** A sequence representation used in the generation of an implied alignment *)
type ias = {
    seq : Sequence.s;    (** The sequence contained *)
    codes : (int, int) Hashtbl.t;
    (** The code of each position in the sequence *)
    homologous: (int, int Sexpr.t) Hashtbl.t;
    (** The other assigned codes that are homologous to a particular code 
    * assigned to [seq] or to other sequences. *)
    indels: (int * string * int * [`Insertion | `Deletion] * int Sexpr.t)
    Sexpr.t; 
        (* The location and contents of an insertion block *)
    
     dum_chars : (int * int Sexpr.t) Sexpr.t;

    order : int list; (* The order of the codes stored in homologous *)

    dir : int; (* dir indicates the orientation of this ias. 
                * dir = 1 if positive, dir = -1 if negative *)
}

type dyna_state_t = Data.dyna_state_t

(* t is the presentation of a dynamic set (DynamicCS) in order to create implied
   alignments.
   Note: sequences : ias_arr Codes.t where 
   isa_arr is an array of ias in the case of annotated chromosomes (ach ias
   presents a locus)
*)

type t = {

    sequences : ias array All_sets.IntegerMap.t;
    c2 : Cost_matrix.Two_D.m;
    chrom_pam : Data.dyna_pam_t;
    state : dyna_state_t;
    code : int;
    children : int Sexpr.t;
    cannonic_code : int;
    alpha : Alphabet.a;
}

type pairs = int * t list

(** A unique code generation function *)
type cg = (unit -> int)

(** Creates a fresh code generation functions *)
val code_generator : unit -> cg

(** [create_ias status s cg] creates a fresh sequence representation for implied
 * alignments of sequence [s] using the code generation function [cg].  
 * state indicate the type of sequence, i.e., Sequence, Chromosome, Annotated,
 * Genome, Breakinv....
 *)
val create_ias : dyna_state_t -> Sequence.s -> int -> cg -> ias

exception IsSankoff

type matrix_class = 
    | AllOne of int
    | AllOneGapSame of (int * int)
    | AffinePartition of (int * int * int)
    | AllSankoff of (string -> int) option

val analyze_tcm :
    Cost_matrix.Two_D.m -> Alphabet.a ->
        matrix_class *
        ([`Exists | `Missing ] -> int -> Parser.t list -> Parser.t list) *
        (int -> (Alphabet.a * Parser.OldHennig.Encoding.s) list -> 
            (Alphabet.a * Parser.OldHennig.Encoding.s) list)

module type S = sig
    type a 
    type b
    type tree = (a, b) Ptree.p_tree


    (** [of_tree t] generates the implied alignment of all the sequences in the tree
    * [t], where the first element of the tuple [t] is a function that removes a
    * gap from the alphabet set of the sequences contained in tree (it implies
    * that [t] only has one sequence). *)
    val of_tree : ((int -> int) * tree) -> Methods.implied_alignment

    val create : (tree -> int list -> tree) ->
        int list -> Data.d ->
        tree -> Methods.implied_alignment list
    

    val to_static_homologies : bool ->
        (tree -> int list -> tree) ->
            bool  -> Methods.characters -> Data.d -> tree -> Data.d

end

module Make (Node : NodeSig.S) (Edge : Edge.EdgeSig with type n = Node.n) : 
    S with type a = Node.n with type b = Edge.e
