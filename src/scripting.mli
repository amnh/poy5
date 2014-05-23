(* POY 5.1.1. A phylogenetic analysis program using Dynamic Homologies.       *)
(* Copyright (C) 2014 Andrés Varón, Lin Hong, Nicholas Lucaroni, Ward Wheeler,*)
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

type model_data =
  (int list * MlModel.spec) list * 
    (string, ((string, float) Hashtbl.t) All_sets.IntSetMap.t) Hashtbl.t option

type ('a, 'b) run = {
    description : string option;
    trees : ('a, 'b) Ptree.p_tree Sexpr.t;
    data : Data.d;
    nodes : 'a list;
    bremer_support : Methods.support_tree Sexpr.t;
    jackknife_support : support_class;
    bootstrap_support : support_class;
    runtime_store : (('a, 'b) run) str_htbl;
    data_store : (Data.d * 'a list) str_htbl;
    bremer_store : Methods.support_tree Sexpr.t str_htbl;
    bootstrap_store : support_class str_htbl;
    jackknife_store : support_class str_htbl;
    tree_store : ('a, 'b) Ptree.p_tree Sexpr.t str_htbl;
    queue : Sampler.ft_queue;
    stored_trees : (float * model_data * Tree.u_tree) Sexpr.t;
    original_trees : ('a, 'b) Ptree.p_tree Sexpr.t;
    search_results : search_results;
}

val build_has : Methods.cost_calculation -> Methods.build -> bool

module type S = sig
        type a 
        type b

type tree = (a, b) Ptree.p_tree 

    module Kml : sig
        type phylogeny = tree

        module GIS : sig 

            type point = { 
                latitude : float;
                longitude : float;
                altitude : float;
            }

            type triangle = (point * point * point)

            (* Estimate the horizontal distance between points, in meters *)
            val horizontal_distance : point -> point -> float
            
            (* Calculate a point located in the center between two points *)
            val center_points : point -> point -> point

            (* Find the centroid of a triangle *)
            val center_triangle : triangle -> point

        end

        module TemporalGIS : sig
            type date = (int * int * int) option (* Year month day *)

            type sample = {
                coordinates : GIS.point;
                date : date;
            }

            (* Compare a pair of dates and return the minimum *)
            val min_date : date -> date -> date

            (* Produce a hash table of terminals and their corresponding samples, as
            * read from a csv file *)
            val csv : string -> (Xml.unstructured, sample) Hashtbl.t
        end

        module KTree : sig
            (* A module to easily process trees for KML generation *)

            (* POY uses internally a relatively complex data structure to hold the
            * trees, this is a sipmlified version that has all the information in XML
            * like format. Only binary trees are alowed, and each is a tuple, consisting
            * of the contents of the vertex in the phylogenetic tree, and the temporal
            * and geographic information associated with it. The leaves are exactly the
            * input data, while the interior vertices are computed by POY or user
            * provided plugins. *)
            type simplified_topology = (Xml.xml * TemporalGIS.sample) Tree.Parse.t

            (* The represenatation of the name of a node. We don't use plain strings
            * because they would make the generation of the XML a little bit too verbose
            * *)
            type node_name = Xml.unstructured 

            (* The topology of a tree, with the simplified version of the topology, and
            * quick access to the nodes of the tree, and the ancestors. This suplies
            * some functions that the simplified_topology can nos perform (like finding
            * the ancestor of a tree or quickly reaching a particular vertex of the
            * tree). *)
            type topology =
                { ancestors : (Xml.unstructured, Xml.xml Xml.contents option) Hashtbl.t;
                nodes : (Xml.unstructured, Xml.xml) Hashtbl.t;
                topo : simplified_topology }


            (* The default tree adjustment function *)
            val adjust_tree : simplified_topology -> simplified_topology

            (* [process data csv phylogeny] producess a topology consisting of the
            * contents computed in the [phylogeny] tree, with temporal and geographic
            * information contained in the CSV file [csv], and all the data
            * representation [data]. *)
            val process : Data.d -> string -> phylogeny -> topology

            (* [ancestor topology vertex] gets the ancestor of the [vertex] in the
            * [topology]. The output is optional as the root of the tree has no
            * ancestor. *)
            val ancestor : topology -> node_name -> Xml.xml Xml.contents option

            (* [children topology vertex] gets the pair of the vertex [vertex] in the
            * [topology]. The output is optional as the leaves of the tree have no
            * children. *)
            val children : topology -> node_name -> (node_name * node_name) option

            (* [sister topology vertex] gets the sister group of the [vertex] in the
            * [topology] (that is, the other child of the ancestor of [vertex]). 
            * The output is optional as the root has no sister. *)
            val sister : topology -> node_name ->  node_name option

            (** [node topology vertex] extracts all the data about [vertex] contained 
            * in the original phylogeny as stored in [topology]. *)
            val node : topology -> node_name -> Xml.xml


            (** [is_root topology vertex] is true iff [vertex] is the root of the
            * [topology] *)
            val is_root : topology -> node_name -> bool

            (** [extract_gis simplified_topology] extracts the temporal and gis
            * information stored in the root vertex of the [simplified_topology]. *)
            val extract_gis : simplified_topology -> TemporalGIS.sample

        end

        module KFile : sig
            (** The following types are needed to produce a Plugin for POY. *)

            (** [node_information data topology vertex] produces an HTML-equivalent
            * structure with the information that should be printed about the [vertex] in
            * the [topology]. [data] is provided in case the specification of some of
            * the characters in [vertex] or the [vertex] itself is needed. *)
            type node_information = 
                Data.d -> KTree.topology -> Xml.unstructured -> 
                    [ Xml.unstructured | Xml.xml Xml.structured ]

            (** [create_node node_information data topology vertex parent_sample
            *   child1_sample child2_sample vertex_sample] produces the XML structure
            *   with the representation of [vertex] in the KML file. The information
            *   about the vertex should be generated and enclosed in a CDATA using the
            *   [node_information] function provided in the argument. The
            *   [parent_sample], [child1_sample], and [child2_sample] are provided for
            *   convenience, as the main goal of this function is to print the node and
            *   edges connected with it. *)
            type create_node =
                    node_information -> Data.d -> KTree.topology ->
                    Xml.xml -> TemporalGIS.sample option -> 
                        TemporalGIS.sample option ->
                        TemporalGIS.sample option ->TemporalGIS.sample -> 
                            Xml.xml Sexpr.t


            (** [adjust_tree simple_topology] beautifies the location of the vertices in
             * the tree *)
            type adjust_tree = KTree.simplified_topology -> KTree.simplified_topology

            (** [styles ()] produces all the styles used in the KML. *)
            type styles = unit -> Xml.xml Sexpr.t

            type folder = {
                name : string;
                node_information : node_information;
                create_node : create_node option;
            }

            type plugin = {
                folders : folder list;
                adjust_tree : adjust_tree;
                styles : styles;
            }

            (** [default ] is the default plugin *)
            val default : plugin

            (** [register_plugin name plugin] registers the [plugin] under the [name]
             * provided. This plugin will be usable in the user interface using the
             * command report (kml:name). *)
            val register_plugin : string -> plugin -> unit

            (** [has_plugin name] is [true] iff [register_plugin name plugin] has been
            * called before *)
            val has_plugin : string -> bool

            (** [kml ?plugin name output data csv tree] dumps in the file [output] a
            * KML file using the [plugin] selected for the tree [tree] and using the
            * [csv] file with the geographic and temporal information. The KML will be
            * registerd with the [name] provided. If not [plugin] is given, then
            * [default] is selected .*)
            val kml : ?plugin:string -> string -> string -> Data.d -> string ->
                phylogeny Sexpr.t -> unit


            val create_line : TemporalGIS.sample option -> TemporalGIS.sample -> 
                Xml.xml Xml.contents
        end
    end

type r = (a, b) run

    val register_function : 
        string -> (Methods.script Methods.plugin_arguments -> r -> r) -> unit

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

val rediagnose_trees : ?classify:bool -> bool -> r -> r

val set_console_run : r -> unit

    module PhyloTree : sig
        type phylogeny = (a, b) Ptree.p_tree
        val get_cost :  phylogeny -> float
        val get_root_costs : phylogeny -> (Tree.edge * float) list
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
        val cost : Tree.join_jxn -> Tree.join_jxn -> a -> phylogeny -> Ptree.clade_cost

        (* Recomputing the contents of a tree *)
        val downpass : phylogeny -> phylogeny
        val uppass : phylogeny -> phylogeny

        (* Tree conversion for IO *)
        val of_string : string -> Data.d -> a list -> phylogeny list
        val to_string : 
            Methods.report_branch option -> phylogeny -> Data.d -> string list
        val of_file : string -> Data.d -> a list -> phylogeny list
        val of_nodes : Data.d -> a list -> phylogeny
        val branch_table : Methods.report_branch option -> phylogeny ->
            ((int * int),[ `Name of (int array * float option) list | `Single of float ]) Hashtbl.t

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
        val nodes : unit -> a list
        val to_string : Methods.report_branch option -> string list list
        val of_string : string -> unit
    end

    module Node : NodeSig.S with type n = a


end

module Make 
    (Node : NodeSig.S with type other_n = Node.Standard.n) 
    (Edge : Edge.EdgeSig with type n = Node.n)
    (TreeOps : 
        Ptree.Tree_Operations with type a = Node.n with type b = Edge.e) :
        S with type a = Node.n with type b = Edge.e 

