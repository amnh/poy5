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

(** {1 Scripting POY} *)
type support_class = (int * int Tree.CladeFPMap.t) option

type 'a str_htbl = (string, 'a) Hashtbl.t


type search_results = {
    tree_costs_found : int All_sets.FloatMap.t;
    total_builds : int;
    total_fuse : int;
    total_ratchet : int;
}

type ('a, 'b, 'c) run = {
    description : string option;
    trees : ('a, 'b) Ptree.p_tree Sexpr.t;
    data : Data.d;
    nodes : 'a list;
    characters : ('c Sexpr.t, float Sexpr.t) 
        Methods.character_input_output;
    bremer_support : Methods.support_tree Sexpr.t;
    jackknife_support : support_class;
    bootstrap_support : support_class;
    runtime_store : (('a, 'b, 'c) run) str_htbl;
    data_store : (Data.d * 'a list) str_htbl;
    bremer_store : Methods.support_tree Sexpr.t str_htbl;
    bootstrap_store : support_class str_htbl;
    jackknife_store : support_class str_htbl;
    tree_store : ('a, 'b) Ptree.p_tree Sexpr.t str_htbl;
    queue : Sampler.ft_queue;
    stored_trees : ('a, 'b) Ptree.p_tree Sexpr.t;
    original_trees : ('a, 'b) Ptree.p_tree Sexpr.t;
    search_results : search_results;
}

val build_has : Methods.cost_calculation -> Methods.build -> bool

module type S = sig
        type a 
        type b
        type c

type tree = (a, b) Ptree.p_tree 

type r = (a, b, c) run

type minimum_spanning_tree = tree 
type build = minimum_spanning_tree list
type minimum_spanning_family = minimum_spanning_tree list
type build_optimum = tree list

type script = Methods.script

val empty : unit -> r

val args : string array

val run : 
    ?folder:(r -> script -> r) ->
    ?output_file:string -> ?start:r -> script list -> r

val update_mergingscript : (r -> script -> r) -> script list -> r -> r -> r

val process_input : r -> 
    Methods.input -> r

val get_dump : ?file:string -> unit -> r * script list

val restart : ?file:string -> unit -> r

val process_random_seed_set : r -> int -> r

val console_run : string -> unit

val parsed_run : script list -> unit

val channel_run : in_channel -> unit

val get_console_run : unit -> r

val update_trees_to_data : ?classify:bool -> bool -> bool -> r -> r

val set_console_run : r -> unit

    module PhyloTree : sig
        type phylogeny = (a, b) Ptree.p_tree
        val get_cost :  phylogeny -> float
        (* Operating on edges and vertices *)
        val fold_edges : ('a -> Tree.edge -> 'a) -> 'a -> (a, b) Ptree.p_tree -> 'a
        val fold_nodes : ('a -> Tree.node -> 'a) -> 'a -> (a, b) Ptree.p_tree -> 'a
        val fold_vertices : ('a -> int -> 'a) -> 'a -> (a, b) Ptree.p_tree -> 'a
        val add_node_data : int -> a -> phylogeny -> phylogeny
        val get_node_data : int -> phylogeny -> a
        val add_edge_data : Tree.edge -> b -> phylogeny -> phylogeny
        val get_edge_data  : Tree.edge -> phylogeny -> b
        val get_parent : int -> phylogeny -> int
        val get_neighs : int -> phylogeny -> int list

        (* Modifying a tree *)
        val join : 
            Tree.join_jxn -> Tree.join_jxn -> phylogeny -> 
                phylogeny * Tree.join_delta
        val break : Tree.break_jxn -> phylogeny -> phylogeny * Tree.break_delta
        val reroot : Tree.edge -> phylogeny -> phylogeny

        (* Recomputing the contents of a tree *)
        val downpass : phylogeny -> phylogeny
        val uppass : phylogeny -> phylogeny

        (* Tree conversion for IO *)
        val of_string : string -> Data.d -> a list -> phylogeny list
        val to_string : bool -> phylogeny -> Data.d -> string list
        val of_file : string -> Data.d -> a list -> phylogeny list
        val of_nodes : Data.d -> a list -> phylogeny

        (* Swapping a tree *)
        val build : Data.d -> a list -> phylogeny list
        val spr : ((phylogeny * float) list -> unit) -> Data.d -> phylogeny ->
            phylogeny list 
        val tbr : ((phylogeny * float) list -> unit) -> Data.d -> phylogeny ->
            phylogeny list 
    end
    module Runtime : sig
        type phylogeny = (a, b) Ptree.p_tree
        val min_cost : unit -> float option
        val max_cost : unit -> float option
        val all_costs : unit -> float list
        val trees : unit -> phylogeny list
        val set_trees : phylogeny list -> unit
        val data : unit -> Data.d
        val to_string : bool -> string list list 
        val of_string : string -> unit
    end


end

module Make 
    (Node : NodeSig.S with type other_n = Node.Standard.n) 
    (Edge : Edge.EdgeSig with type n = Node.n)
    (TreeOps : 
        Ptree.Tree_Operations with type a = Node.n with type b = Edge.e)
    (CScrp : CharacterScripting.S with type n = Node.n) : 
        S with type a = Node.n with type b = Edge.e with type c = CScrp.cs

(** {2 Scripting in Ocaml} 
 *
 * The following utility functions let us handle efficiently some common
 * operations required when writting test or experimentation scripts. 
 * *)

(** Handling multiple files efficiently *)
module FILES : sig
    (** A simplified interface for efficient file handling when doing scripts *)

    (** [explode pattern] generates a list of all the files matchin the (bash)
     * [pattern] *)
    val explode : string -> string list

    (** [open_all_in p] opens all files matching the pattern [p] for reading,
    * returning their list, containing tuples: the filename and the input
    * channel *)
    val open_all_in : string -> (string * in_channel) list

    (** [open_all_out p] is simmilar to [open_all_in p] but the resulting
    * channel is an output channel *)
    val open_all_out : string -> (string * out_channel) list

    (** [run_n_close n f] runs the function f in the input channel opened for
    * the file [n], closes the channel and returns the result of [f]. *)
    val run_n_close : string -> (in_channel -> 'a) -> 'a
end

module DNA : sig
    (** Simplified interface to deal with DNA or RNA sequences *)

    module CM : sig
        (** Handling transformation cost matrices for DNA or RNA sequences *)

        (** The type of a two dimensional transformation cost matrix *)
        type cm2 = Cost_matrix.Two_D.m

        (** [of_list lst] converts the 5x5 matrix [lst] to the corresponding
        * transformation cost matrix. The function assumes that rows and columns
        * are sorted as follows: A, C, G, T/U, Gap. *)
        val of_list : int list list -> cm2

        (** [of_array arr] is the same as [of_list] but for arrays. *)
        val of_array : int array array -> cm2

        (** [of_sub_indel a b] assigns the same cost [a] to all substitutions,
        * and the same cost [b] to all indels. *)
        val of_sub_indel : int -> int -> cm2

        (** [of_sub_indel a b c] assigns the same cost [a] to all substitutions,
        * and the same cost [b] to all indels, and [c] cost to opening a block
        * of indels. *)
        val of_sub_indel_affine : int -> int -> int -> cm2

        (** [of_file str] reads the input file [str], containing a sequence of
        * 25 elements, each corresponding to the left-right, top-down 25
        * elements that would be required for [of_list]. *)
        val of_file : string -> cm2

        (** [all_ones ()] is equivalent to [of_sub_indel 1 1] *)
        val all_ones : unit -> cm2

        (** [median a b cm] claculates the median between the bases [a] and [b]
         * according to the cost matrix [cm]. *)
        val median : char -> char -> cm2 -> char

        (** [cost a b cm] claculates the cost between the bases [a] and [b]
         * according to the cost matrix [cm]. *)
        val cost : char -> char -> cm2 -> int
    end

    module Seq : sig
        (** The type of a DNA or RNA sequence *)
        type s = Sequence.s

        (** The type of a DNA or RNA base *)
        type base = char

        (** The internal representation of a gap *)
        val gap : base

        (** [of_string str] ~onverts a valid IUPAC string that represents a
        * sequence [str] into the internal sequence representation. *)
        val of_string : string -> s

        (** [to_string s] converts the sequence [s] in the valid IUPAC
        * representation string. *)
        val to_string : s -> string

        (** [length s] returns the length of the sequence s. *)
        val length : s -> int

        (** [prepend s b] prepends the base [b] to the sequence [s], returning a
        * fresh sequence. *)
        val prepend : s -> base -> s

        (** [merge a b] merges the sequences [a] and [b] in a fresh sequence. *)
        val merge : s -> s -> s

        (** [get s p] returns the base in position [p] of sequence [s]. *)
        val get : s -> int -> base

        (** [set s b p] returns a fresh sequence which is a copy of [s]
         * excepting that that the base in position [p] is set to the base [b].
         * *)
        val set : s -> base -> int -> s 

        (** [delete s p] returns a fresh sequence that has deleted the base in
         * position [p] of [s]. *)
        val delete : s -> int -> s

        (* [slice s b e] returns a fresh sequence with the subsequence contained
        * between positions [b] and [e] (inclusive). *)
        val slice : s -> int -> int -> s
    end

    module Fasta : sig
        (** Simplified parsing and generation of FASTA files. *)

        (** The type of lists of sequences from a FASTA file, as tuples of names
        * and sequences. *)
        type seqs = (string * Seq.s) list
        (** The type of sequences from a file when it has been broken using
        * pipes and pound signs (POY format). Each column of those files
        * correspond to a list.*)
        type multi_seqs = seqs list

        (** [of_channel ch] returns the sequences contained in the FASTA channel
        * [ch]. *)
        val of_channel : bool -> in_channel -> seqs

        (** [multi_of_channel ch] is the same as [of_channel] but returns
        * [multi_seqs]. *)
        val multi_of_channel : in_channel -> multi_seqs

        (** [to_channel ch seqs] outputs all the sequences in [seqs] to the
        * channel [ch] in FASTA format. *)
        val to_channel : out_channel -> seqs -> unit

        (** [of_file f] reads the file [f] containing fasta sequences. *)
        val of_file : bool -> string -> seqs

        (** [multi_of_file ch] is the same as [of_file] but returns
        * [multi_seqs]. *)
        val multi_of_file : string -> multi_seqs

        (** [to_file f s] creates a fresh file [f] containing the sequences [s]
         * in FASTA format. *)
        val to_file : string -> seqs -> unit


        (** [random_sample seqs name number size] generates [number] random samples 
        * of the sequences [seqs], each simple of size [size] without
        * repetitions, and outputs them in the filename produced by [name ()],
        * sequentially. *)
        val random_sample : seqs -> (unit -> string) -> int -> int -> unit

        (** [multi_sample seqs column name number size] is the same as
        * [random_sample (List.nth seqs column) name number size]. *)
        val multi_sample : multi_seqs -> int -> (unit -> string) -> int -> int -> unit
    end

    module Generic : sig
        (** Generic parsing of files (auto guess the type). *)

        (** [molecular f] attempts to guess the type of [f] and if successful
        * returns the parsed sequences. *)
        val molecular : string -> Fasta.seqs
    end

    module Align : sig
        (** Simplified sequence alignment *)

        (** [algn a b cm] returns a quadruple [(a', b', c, med)] where [a'] and
        * [b'] are the edited sequences of [a] and [b] respectively, using the
        * cost matrix [cm], that achieve the optimal edition cost [c]. [med] is
        * a valid median located between [a] and [b]. *)
        val algn : 
            Seq.s -> Seq.s -> CM.cm2 -> Seq.s * Seq.s * int * Seq.s

        (** [algn_and_print a b cm] is the same as [algn] excepting that the
        * sequences are not returned in the abstract representation but as OCaml
        * strings. *)
        val algn_and_print :
            Seq.s -> Seq.s -> CM.cm2 -> string * string * int * string

        (** [algn_all seqs cm] takes the list of sequences [seqs] and align them
        * all to all, producing a matrix of alignments as generated by [algn].*)
        val algn_all : Fasta.seqs -> CM.cm2 -> 
            ((string * Seq.s) * (string * Seq.s) * int * Seq.s) list list

        (** [algn_all_and_print seqs cm] is the same as [algn_all] excepting
        * that the sequences are not returned in the abstract representation but
        * as OCaml strings. *)
        val algn_all_and_print : Fasta.seqs -> CM.cm2 -> ((string * string) *
        (string * string) * int * string) list list

    end

end
