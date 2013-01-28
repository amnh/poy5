(* POY 5.0 Alpha. A phylogenetic analysis program using Dynamic Homologies.   *)
(* Copyright (C) 2011 Andrés Varón, Lin Hong, Nicholas Lucaroni, Ward Wheeler,*)
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

let () = SadmanOutput.register "Scripting" "$Revision: 3048 $"

module IntSet = All_sets.Integers

type support_class = (int * int Tree.CladeFPMap.t) option

type 'a str_htbl = (string, 'a) Hashtbl.t

type search_results = {
    tree_costs_found : int All_sets.FloatMap.t;
    total_builds : int;
    total_fuse : int;
    total_ratchet : int;
}

(* Before we begin, set the correct error printing functions in the parsers *)
let () = Nexus.P.print_error := Status.user_message Status.Error

(* We define a few functions to deal with the search results type *)
let empty_search_results = {
    tree_costs_found = All_sets.FloatMap.empty;
    total_builds = 0;
    total_fuse = 0;
    total_ratchet = 0;
}

let incr_search_results t cnt sr = match t with
    | `Builds   -> { sr with total_builds   = sr.total_builds + cnt }
    | `Fuse     -> { sr with total_fuse     = sr.total_fuse + cnt   }
    | `Ratchet  -> { sr with total_ratchet  = sr.total_ratchet + cnt}

let add_results cost sr =
    let cnt = 
        if All_sets.FloatMap.mem cost sr.tree_costs_found then
            1 + (All_sets.FloatMap.find cost sr.tree_costs_found)
        else 1
    in
    let set = All_sets.FloatMap.add cost cnt sr.tree_costs_found in
    { sr with tree_costs_found = set }

type ('a, 'b) run = {
    description : string option;
    trees : ('a, 'b) Ptree.p_tree Sexpr.t;
    data : Data.d;
    nodes : 'a list;
    bremer_support : Methods.support_tree Sexpr.t;
    jackknife_support : support_class;
    bootstrap_support : support_class;
    runtime_store : (('a, 'b) run) str_htbl;
    data_store : (Data.d * ('a list)) str_htbl;
    bremer_store : Methods.support_tree Sexpr.t str_htbl;
    bootstrap_store : support_class str_htbl;
    jackknife_store : support_class str_htbl;
    tree_store : ('a, 'b) Ptree.p_tree Sexpr.t str_htbl;
    queue : Sampler.ft_queue;
    stored_trees : ('a, 'b) Ptree.p_tree Sexpr.t;
    original_trees : ('a, 'b) Ptree.p_tree Sexpr.t;
    search_results : search_results;
}

let is_forest = function
    | `LocalOptimum (_, _, _, _, _, x, _, _, _, _, _) -> x

let is_something y x = x = y

let has_something x = List.exists (is_something x)

let build_has item = function
    | `Mst _ 
    | `Nj
    | `Prebuilt _ -> false
    | `Branch_and_Bound ((_, _, _, _, l),_) 
    | `Build (_, _, l,_)
    | `Build_Random ((_, _, _, l, _),_) -> has_something item l


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

    val register_function : string -> (Methods.script Methods.plugin_arguments -> r -> r) -> unit

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
    
    val rediagnose_trees : ?classify:bool -> bool -> r -> r

    val set_console_run : r -> unit

    module PhyloTree : sig
        type phylogeny = (a, b) Ptree.p_tree
        val get_cost : phylogeny -> float
        val get_root_costs : phylogeny -> (Tree.edge * float) list
        val fold_edges : ('a -> Tree.edge -> 'a) -> 'a -> (a, b) Ptree.p_tree -> 'a
        val fold_nodes : ('a -> Tree.node -> 'a) -> 'a -> (a, b) Ptree.p_tree -> 'a
        val fold_vertices : ('a -> int -> 'a) -> 'a -> (a, b) Ptree.p_tree -> 'a
        val add_node_data : int -> a -> phylogeny -> phylogeny
        val get_node_data : int -> phylogeny -> a
        val add_edge_data : Tree.edge -> b -> phylogeny -> phylogeny
        val get_edge_data  : Tree.edge -> phylogeny -> b
        val get_parent : int -> phylogeny -> int
        val get_neighs : int -> phylogeny -> int list
        val join : Tree.join_jxn -> Tree.join_jxn -> phylogeny -> phylogeny * Tree.join_delta
        val break : Tree.break_jxn -> phylogeny -> phylogeny * Tree.break_delta
        val reroot : Tree.edge -> phylogeny -> phylogeny
        val downpass : phylogeny -> phylogeny
        val branch_table : phylogeny -> ((int * int),[ `Name of (int array * float option) list | `Single of float ]) Hashtbl.t
        val uppass : phylogeny -> phylogeny
        val of_string : string -> Data.d -> a list -> phylogeny list
        val to_string : bool -> phylogeny -> Data.d -> string list
        val of_file : string -> Data.d -> a list -> phylogeny list
        val of_nodes : Data.d -> a list -> phylogeny
        val build : Data.d -> a list -> phylogeny list
        val spr : ((phylogeny * float) list -> unit) -> Data.d -> phylogeny -> phylogeny list 
        val tbr : ((phylogeny * float) list -> unit) -> Data.d -> phylogeny -> phylogeny list 
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
        val to_string : bool -> string list list 
        val of_string : string -> unit
    end

    module Node : NodeSig.S with type n = a

end


module Make (Node : NodeSig.S with type other_n = Node.Standard.n)
            (Edge : Edge.EdgeSig with type n = Node.n)
            (TreeOps : Ptree.Tree_Operations with type a =
                                      Node.n with type b = Edge.e)
    = struct
    type a = Node.n
    type b = Edge.e
    type tree = (a, b) Ptree.p_tree 

    module Kml = struct
        exception IllegalCVS

        type phylogeny = tree

        module GIS = struct

            type point = { 
                latitude : float; 
                longitude : float;
                altitude : float;
            }

            let horizontal_distance a b =
                let to_rad x = x /. 57.2958 (* 180 / pi *) in
                let convert a = 
                    { a with latitude = to_rad a.latitude; 
                    longitude = to_rad a.longitude }
                in
                let a = convert a
                and b = convert b in
                let r = 63787000. (* Kilometers *) in
                let dist = 
                    r *. acos (((sin a.latitude) *. (sin b.latitude)) +.
                    ((cos a.latitude) *. (cos b.latitude) *. (cos (b.longitude -.
                    a.longitude))))
                in
                dist

            type triangle = (point * point * point)

            let latitude_middle a b = (a.latitude +. b.latitude) /. 2. 

            let longitude_middle a b = 
                let dif = a.longitude -. b.longitude in
                let mid = (a.longitude +. b.longitude) /. 2. in
                if dif > 180. then mid -. 180.
                else if dif < (-. 180.) then mid +. 180.
                else mid

            let altitude_middle a b = 
                (a.altitude +. b.altitude) /. 2.

            let altitude_delta a b =
                let d = horizontal_distance a b in
                if d < 50000. then d
                else 50000.

            let center_points a b =
                { latitude = latitude_middle a b; 
                longitude = longitude_middle a b;
                altitude = altitude_middle a b}

            let point_mid a b = 
                let max_altitude = 
                    (altitude_delta a b) +. (max a.altitude b.altitude)
                in
                { latitude = latitude_middle a b; 
                longitude = longitude_middle a b;
                altitude = max_altitude; }

            let rec point_mid2 i a b =
                if i = 0 then b
                else 
                    let b = point_mid a b in
                    point_mid2 (i - 1) a b


            let center_triangle (a, b, c) =
                let rec aux iterations (a, b, c) =
                    if iterations = 0 then a
                    else
                        let iterations = iterations - 1 in
                        if a = b && b = c then a
                        else aux iterations 
                            (point_mid a b, point_mid b c, point_mid a c)
                in
                aux 5 (a, b, c)

            let rec triangle_mid iterations (a, b, c) = 
                if iterations = 0 then a
                else
                    let iterations = iterations - 1 in
                    if a = b && b = c then a
                    else triangle_mid iterations 
                        (point_mid a b, point_mid b c, point_mid a c)

            let triangle_mid (a, b, c) =
                let ab = horizontal_distance a b
                and bc = horizontal_distance b c
                and ac = horizontal_distance a c in
                if 10000. > min (min ab bc) ac then
                    if ab < bc then
                        if ab < ac then (point_mid a b) 
                        else (point_mid a c) 
                    else if bc < ac then (point_mid b c) 
                    else (point_mid a c) 
                else triangle_mid 5 (a, b, c)

        end

        module TemporalGIS = struct

            open GIS 

            type date = (int * int * int) option 
            type sample = { coordinates : GIS.point; 
            date : date}

            let min_date a b =
                match a, b with
                | x, None
                | None, x -> x
                | Some (yeara, ma, da), Some (yearb, mb, db) ->
                        if yeara = yearb then
                            if ma = mb then
                                if da = db then a
                                else if da < db then a
                                else b
                            else if ma < mb then a
                            else b
                        else if yeara < yearb then a
                        else b

            let initial_ancestor x y =
                { coordinates = point_mid x.coordinates y.coordinates;
                date = min_date x.date y.date }

            let adjust_ancestor mine x y z =
                let altitude = mine.coordinates.altitude in
                let coord = 
                    triangle_mid (x.coordinates, y.coordinates, z.coordinates)
                in
                { mine with coordinates = { coord with altitude = altitude }}

            let csv file =
                (* We need to read the file and output the hash table
                with all of * the coordinates and time if available *)
                let ch = open_in file in
                let res = Hashtbl.create 1667 in
                try
                    (* We read the first line to check what is the input 
                    * contents *)
                    let line = input_line ch in
                    let comma = Str.regexp ","
                    and space = Str.regexp " " in
                    let with_date name lat long year month day =
                        Hashtbl.add res (`String name)
                        { coordinates = 
                            { latitude = lat; longitude = long; altitude = 0.}; 
                            date = Some (year, month, day) }
                    and without_date name lat long =
                        Hashtbl.add res (`String name)
                        {coordinates = 
                            {latitude = lat; longitude = long; altitude = 0.};
                        date = None}
                    in
                    let has_date = 
                        match Str.split comma (Str.global_replace space "" line) with
                        | [a; b; c; d] ->
                                let aref = "strain_name" in
                                if aref = aref then 
                                    if b = "latitude" then 
                                        if c = "longitude" then
                                            if d = "date" then true
                                            else raise IllegalCVS
                                        else raise IllegalCVS
                                    else raise IllegalCVS
                                else 
                                    raise IllegalCVS
                        | [a; b; c] ->
                                if 
                                    a = "strain_name" && b = "latitude" && 
                                    c = "longitude" then false
                                else raise IllegalCVS
                        | "strain_name" :: "latitude" :: "longitude" :: tl  ->
                                (match tl with
                                | [] -> false
                                | ["date"] -> true
                                | _ -> raise IllegalCVS)
                        | lst -> raise IllegalCVS
                    in
                    while true do
                        let line = input_line ch in
                        match Str.split comma line with
                        | [a; b; c] when not has_date ->
                                without_date a (float_of_string b) (float_of_string c)
                        | [a; b; c; d] when has_date ->
                                (match Str.split (Str.regexp "-") d with
                                | [x; y; z] ->
                                        with_date a (float_of_string b) 
                                        (float_of_string c) (int_of_string x) 
                                        (int_of_string y) (int_of_string z)
                                | _ -> 
                                        prerr_string line;
                                        raise IllegalCVS)
                        | _ -> 
                                prerr_string line;
                                raise IllegalCVS
                    done;
                    assert false
                with
                | End_of_file -> 
                        close_in ch;
                        res
        end


        module KTree = struct

            type node_name = Xml.unstructured 
            let index_formatter tree =
                let tree = TreeOps.to_formatter `Normal [] tree in
                let nodes =
                    let rec get_nodes acc ((x, _, sub) as n) =
                        match x with
                        | "Node" -> (n :: acc)
                        | _ ->
                                match sub with
                                | #Sexpr.t as sub ->
                                        Sexpr.fold_left get_nodes  acc sub
                                | _ -> acc
                    in
                    get_nodes [] tree
                in
                let tbl = Hashtbl.create 1667 in
                List.iter (fun node -> 
                    let name = Xml.attribute Xml.Nodes.name node in
                    Hashtbl.add tbl name node) nodes;
                tbl

            let extract_gis x = 
                snd 
                (match x with  Tree.Parse.Leafp x | Tree.Parse.Nodep (_, x) -> x)

            let empty = 
                { TemporalGIS.coordinates = 
                    { GIS.latitude = 0.; longitude = 0.; altitude = 0.};
                date = None }

            let create_simplified_topology gis_table nodes data tree =
                let get_name code =
                    try `String (Data.code_taxon code data) with
                    | Not_found -> `String (string_of_int code)
                in
                let get_node code =
                    try Hashtbl.find nodes (get_name code)
                    with er -> raise er
                in
                let get_gis code =
                    let name = get_name code in
                    try Hashtbl.find gis_table name
                    with er -> raise er
                in
                match Ptree.get_roots tree with
                | [x] ->
                    begin match x.Ptree.root_median with
                        | Some ((`Single x), root) ->
                            let node = get_node x in
                            let gis = get_gis x in
                            Tree.Parse.Leafp (node, gis)
                        | Some ((`Edge edge), root) ->
                            let leaf _ cur _ =
                                let node = get_node cur in
                                let gis = get_gis cur in
                                Tree.Parse.Leafp (node, gis)
                            and internal _ cur x y =
                                let xgis = extract_gis x
                                and ygis = extract_gis y in
                                let my_gis = TemporalGIS.initial_ancestor xgis ygis in
                                let node = get_node cur in
                                Tree.Parse.Nodep ([x; y], (node, my_gis))
                            in
                            let node =
                                try Hashtbl.find nodes (`String "root")
                                with err -> raise err
                            in
                            let left, right = 
                                Tree.post_order_node_with_edge_visit leaf 
                                internal (Tree.Edge edge) tree.Ptree.tree 
                                (Tree.Parse.Leafp (node, empty))
                            in
                            let gis = 
                                TemporalGIS.initial_ancestor (extract_gis left) 
                                (extract_gis right)
                            in
                            Tree.Parse.Nodep ([left; right], (node, gis))
                        | None -> assert false
                    end
                | _ -> failwith "Kml generation requires forest of only one tree"

            let adjust_tree tree =
                let rec downpass iterations tree = 
                    if iterations = 0 then tree
                    else
                        let rec aux parent_gis tree =
                            match tree with
                            | Tree.Parse.Leafp _ -> (tree, false)
                            | Tree.Parse.Nodep ([x; y], (node, gis)) ->
                                    let x, x_changed = aux gis x
                                    and y, y_changed = aux gis y in
                                    let my_gis = 
                                        TemporalGIS.adjust_ancestor gis parent_gis 
                                        (extract_gis x) (extract_gis y) 
                                    in
                                    Tree.Parse.Nodep ([x; y], (node, my_gis)), 
                                    x_changed || y_changed || my_gis <> gis
                            | _ -> assert false
                        in
                        match tree with
                        | Tree.Parse.Leafp _ -> tree
                        | Tree.Parse.Nodep ([x; y], (node, root_gis)) ->
                                let x_gis = extract_gis x in
                                let y, y_changed = aux x_gis y in
                                let y_gis = extract_gis y in
                                let x, x_changed = aux y_gis x in
                                let root_gis = 
                                    TemporalGIS.initial_ancestor (extract_gis x) 
                                    (extract_gis y)
                                in
                                let res = 
                                    Tree.Parse.Nodep ([y; x], (node, root_gis)) in
                                if x_changed || y_changed then 
                                    downpass (iterations - 1) res
                                else res
                        | _ -> assert false
                in
                downpass 30 tree

            let parent_table tree =
                let res = Hashtbl.create 1667 in
                let rec aux parent tree =
                    let name = 
                        match tree with
                        | Tree.Parse.Leafp (x, _) ->
                                Xml.attribute Xml.Nodes.name x
                        | Tree.Parse.Nodep (chld, (x, _)) ->
                                let name = Xml.attribute Xml.Nodes.name x in
                                List.iter (aux (Some (name :> Xml.xml Xml.contents))) chld;
                                name
                    in
                    Hashtbl.add res name parent;
                in 
                aux None tree;
                res

            type simplified_topology = (Xml.xml * TemporalGIS.sample) Tree.Parse.t
            type topology = 
                { ancestors : (Xml.unstructured, Xml.xml Xml.contents option) Hashtbl.t;
                  nodes : (Xml.unstructured, Xml.xml) Hashtbl.t; 
                  topo : simplified_topology; } 


            let process data csv tree =
                let csv = TemporalGIS.csv csv in
                let nodes = index_formatter tree in
                let topo = create_simplified_topology csv nodes data tree in
                let ancestors = parent_table topo in
                { ancestors = ancestors; nodes = nodes; topo = topo }

            let ancestor topo node =
                try Hashtbl.find topo.ancestors node with
                | Not_found -> None

            let children topo node =
                let (_, attrs, _) = Hashtbl.find topo.nodes node in
                match 
                List.filter (fun (x, _) ->
                    x = Xml.Nodes.child1_name || x = Xml.Nodes.child2_name) attrs
                with
                | [] -> None
                | [(_, `String ""); (_, `String "")] -> None
                | [(_, x); (_, y)] -> Some (x, y)
                | _ -> assert false

            let sister topo node =
                let ance = ancestor topo node in
                match ance with
                | None -> None
                | Some (#Xml.unstructured as ance) ->
                    let children = children topo ance in
                    (match children with
                    | Some (a, b) ->
                            if a = node then Some b
                            else begin
                                assert (b = node);
                                Some a
                            end
                    | None -> assert false)
                | Some _ -> assert false

            let node topo node =
                Hashtbl.find topo.nodes node 

            let is_root topo node = None = ancestor topo node
        end

        module KFile = struct

            (* Apply the function on the value if not None, otherwise produce the empty
            * contents *)
            let optional contents v =
                match v with
                | None -> (PXML ---)
                | Some v -> contents v

            (* A function to generate the summary of the immediate neighbors of a node
            * *)
            let family_summary topology node () : Xml.xml Sexpr.t =
                let do_row title contents =
                    let c = 
                        match contents with
                        | `String contents -> contents
                        | _ -> ""
                    in
                    (PXML 
                        -tr -td (align=right) -font (size="+1") [`String title] -- --
                            -td -a 
                                (href=[`String ("#pm_" ^ c)])
                            [(contents :> Xml.xml Xml.contents)] -- --
                        --)
                in
                let make_pair (l, r) = `Set [do_row "Descendant:" l; do_row "Descendant:" r] in
                let ancestor = KTree.ancestor topology node
                and sister = KTree.sister topology node
                and children = KTree.children topology node in
                (PXML
                    -table 
                        -tr -td (colspan=2) (align=center)
                            -b -font (size="+2") (color="#000000")
                                "Family Summary"
                            -- --
                        -- --
                        { optional (do_row "Ancestor:") ancestor }
                        { optional (do_row "Sister:") sister }
                        { optional make_pair children } --)

            let ( --> ) a b = b a

            let is_transition str1 str2 =
                match str1, str2 with
                | "A", "G"
                | "G", "A"
                | "C", "T"
                | "T", "C" -> true
                | _ -> false

            let is_transversion str1 str2 = not (is_transition str1 str2)

            let white = `String "#FFFFFF"
            let gray = `String "#E1E1E1"

            let transforms_text data tree node () =
                (* How will the transformations look like *)
                let ancestor = KTree.ancestor tree node in
                let coherce x  = (x :> Xml.xml Xml.contents) in
                match ancestor with
                | None -> 
                        (* We are at the root, we don't show any transforms *)
                        (PXML ---)
                | Some ancestor ->
                        let ancestor = 
                            match ancestor with
                            | #Xml.unstructured as ancestor -> KTree.node tree ancestor
                            | _ -> assert false
                        and node = KTree.node tree node in
                        (* Good, so we have to traverse the list of characters and do
                        * accordingly *)
                        let process_pair node ancestor =
                            let tag = Xml.tag node in
                            assert (tag = (Xml.tag ancestor));
                            (* We define a function that processes a pair of characters,
                            * we filter what we are interested on *)
                            let use_final = 
                                tag = Xml.Characters.nonadditive || 
                                tag = Xml.Characters.additive || 
                                tag = Xml.Characters.sankoff
                            in
                            let cost = Xml.attribute Xml.Characters.cost node
                            and name = Xml.attribute Xml.Characters.name node
                            and clas = Xml.attribute Xml.Characters.cclass node in
                            assert (clas = Xml.attribute Xml.Characters.cclass ancestor);
                            let costn = match cost with
                                | `FloatFloatTuple (costn, _)
                                | `Float costn -> costn
                                | `IntTuple (costn, _) -> float_of_int costn
                                | _ -> 0.
                            in
                            if costn > 0. && use_final && 
                            ((`String Xml.Nodes.final) = clas) then
                                (* A function to ge the list of states contained 
                                * in a node *)
                                let get_chars cha =
                                    cha --> Xml.children Xml.Characters.value 
                                        --> Sexpr.map Xml.value
                                        --> Sexpr.map Xml.value_to_string
                                        --> Sexpr.to_list
                                        --> String.concat ";" 
                                in
                                let node_chars = get_chars node
                                and ance_chars = get_chars ancestor in
                                let t = 
                                    let name = Xml.value_to_string name in
                                    let code = Data.character_code name data in
                                    let alph = Data.get_alphabet data code in
                                    let alph = Alphabet.to_sequential alph in
                                    let alph = Alphabet.to_list alph in
                                    let alph = List.map fst alph in
                                    let res =
                                        if alph = ["A"; "C"; "G"; "T"] ||
                                            alph = ["A"; "C"; "G"; "T"; "-"] then
                                            if ance_chars = "-" then "Ins"
                                            else if node_chars = "-" then "Del"
                                            else if is_transition node_chars
                                            ance_chars then "Ti"
                                            else if is_transversion node_chars
                                            ance_chars then "Tv"
                                            else "ABC"
                                        else "ABC"
                                    in
                                    `String res
                                    in
                                (PXML
                                -tr -td (align=left) (bgcolor=[white]) { coherce name } --
                                    -td (align=left) (bgcolor=[gray]) -font (size="+1") 
                                        { `String ance_chars } -- --
                                    -td (align=left) (bgcolor=[white]) -font (size="+1") 
                                        { `String node_chars } -- --
                                    -td (align=left) (bgcolor=[gray]) { coherce t } --
                                    -td (align=left) (bgcolor=[white]) { coherce cost } -- --)
                            else if costn > 0. && 
                                (`String Xml.Nodes.single) = clas then
                                let get_chars x = 
                                    x --> Xml.value --> Xml.value_to_string
                                in
                                let ance_chars = get_chars ancestor
                                and node_chars = get_chars node in
                                (PXML
                                -tr -td (align=left) (bgcolor=[white]) { coherce name } --
                                    -td (align=left) (bgcolor=[gray])
                                        -font (size="+1") { `String ance_chars } -- --
                                    -td (align=left) (bgcolor=[white]) -font (size="+1") 
                                        { `String node_chars } -- --
                                    -td (align=left) (bgcolor=[gray]) 
                                        { `String "Sequence" } --
                                    -td (align=left) (bgcolor=[white]) 
                                        { coherce cost } -- --)
                            else 
                                (PXML ---)
                        in
                        let get_list_of_all_characters (_, _, n) =
                            match n with
                            | `Delayed f -> Sexpr.to_list (f ())
                            | #Sexpr.t as n -> Sexpr.to_list n
                            | `CDATA _
                            | #Xml.unstructured -> []
                        in
                        let node_characters = get_list_of_all_characters node 
                        and ance_characters = get_list_of_all_characters ancestor in
                        let title color content = 
                            (PXML -td (align=left) (bgcolor=[color]) 
                                -i { `String content } -- --) 
                        in
                        (PXML 
                            -table (border=0) (align=left) -font (color="#A4BFDB")
                                -tr -td (colspan=3) (align=left) -b -font 
                                    (size="+2") (color="#000000") Transformations
                                -- -- --
                                -tr
                                    { title white "Position" }
                                    { title gray "Ancestor" }
                                    { title white "Descendant" }
                                    { title gray "Type" }
                                    { title white "Cost" }
                                --
                                { set List.map2 process_pair node_characters 
                                ance_characters }
                            -- -- --)

            let create_coords point =
                (string_of_float point.GIS.longitude ^ "," ^ 
                string_of_float point.GIS.latitude ^ "," ^
                string_of_float point.GIS.altitude)

            let line_to_curve iterations a b =
                let accumulator = Buffer.create 1667 in
                let append x =
                    Buffer.add_string accumulator " ";
                    Buffer.add_string accumulator x
                in
                let rec aux iterations a b =
                    let center = GIS.center_points a b in
                    let coords = create_coords center in
                    if iterations = 0 then begin
                        append coords;
                    end else begin
                        let iterations = iterations - 1 in
                        aux iterations a center;
                        append coords;
                        aux iterations center b
                    end
                in
                Buffer.add_string accumulator (create_coords a);
                aux iterations a b;
                append (create_coords b);
                Buffer.contents accumulator

            let create_line parent_gis gis =
                (PXML -LineString
                    -altitudeMode relativeToGround --
                    -coordinates 
                        [ match parent_gis with
                        | None -> `Empty
                        | Some x -> 
                                `String 
                                    (line_to_curve 5 gis.TemporalGIS.coordinates x.TemporalGIS.coordinates)]
                    --
                --)


            let associated_info data tree name =
                let is_root = KTree.is_root tree name in
                (PXML -text -body (bgcolor="#ffffff") 
                    -table (width=400) 
                        -tr 
                            -td { `Delayed (family_summary tree name) } -- 
                            {if is_root then (PXML ---) 
                            else
                               (PXML -td 
                               { `Delayed (transforms_text data tree name) } --)}
                        -- -- -- --)

            (* A function that generates the contents that should be printed in the node
            * ballon box of the KML files.  *)
            type node_information = 
                Data.d -> KTree.topology -> Xml.unstructured -> 
                    [ Xml.unstructured | Xml.xml Xml.structured ]

            let place_node (associated_info : node_information) data tree node parent_gis left_gis right_gis gis 
            : Xml.xml Sexpr.t =
                let name = Xml.attribute Xml.Nodes.name node in
                let is_root = KTree.is_root tree name in
                let id = Xml.value_to_string name in
                let summary : Xml.xml Xml.contents =
                    (PXML { PCDATA {associated_info data tree name} }) 
                in
                (PXML 
                -Placemark 
                        (* Only one attribute *)
                        (id=[`String ("pm_" ^ id)])
                        -visibility --
                        -name { `String id } --
                        -Snippet (maxLines=0) empty --
                        -LookAt
                            - longitude [`Float gis.TemporalGIS.coordinates.GIS.longitude] --
                            - latitude [`Float gis.TemporalGIS.coordinates.GIS.latitude] --
                            - altitude [`Float gis.TemporalGIS.coordinates.GIS.altitude] --
                            - range 500 --
                            - tilt 10 --
                            - heading 0 --
                        --
                        { match gis.TemporalGIS.date with
                        | None -> (PXML ---)
                        | Some (y, m, d) -> 
                                let str = 
                                    String.concat "-" 
                                    (List.map string_of_int [y; m; d])
                                in
                                (PXML -TimeSpan -"begin" {`String str} -- --)}
                        - description { summary } --
                        - styleUrl "#00" --
                        - Style
                            { if is_root then 
                                (PXML 
                                    - BallonStyle 
                                        -textColor "FF000000" --
                                        -bgColor "FFffffff" --
                                    --)
                            else (PXML ---) }
                            - LineStyle 
                                (id=[`String ("ls_" ^ id)])
                                -color "FFffffff" --
                            --
                        --
                        -MultiGeometry
                            -Point (id=[`String id])
                                -altitudeMode relativeToGround --
                                -coordinates 
                                    [ `String (create_coords gis.TemporalGIS.coordinates)  ] --
                            --
                            { create_line parent_gis gis }
                            { create_line left_gis gis }
                            { create_line right_gis gis }
                        --
                    --)

            let default_styles () = 
                let supramap_url suffix =
                    `String ("http://supramap.osu.edu/supramap/" ^ suffix) in
                let a = (PXML 
                    -Style (id="00a") 
                        -IconStyle -scale 0.35 -- 
                            -Icon -href { supramap_url "images/circle.png" } --
                        -- --
                        -LabelStyle -scale 0 -- --
                        -LineStyle -width 1.5 -- --
                    --)
                and b = (PXML
                    -Style (id="00b")
                        -IconStyle -scale 0.75 -- 
                            -Icon -href { supramap_url "images/circle.png" } --
                        -- --
                        -LabelStyle --
                        -LineStyle -width 1.5 -- --
                    --) 
                and c = (PXML 
                    -Style (id="00c")
                        -IconStyle -scale 0.95 -- 
                            -Icon -href { supramap_url "images/circle.png" } --
                        -- --
                        -LabelStyle --
                        -LineStyle -width 4 -- --
                    --)
                and d = (PXML 
                    -StyleMap (id="00") 
                        -Pair -key normal -- -styleUrl "#00a" -- --
                        -Pair -key highlight -- -styleUrl "#00c" -- --
                    --)
                and e = (PXML
                    -StyleMap (id="01") 
                        -Pair -key normal -- -styleUrl "#00b" -- --
                        -Pair -key highlight -- -styleUrl "#00c" -- --
                    --)
                in
                (PXML { set [a; b; c; d; e] })

            let base_ml style name contents : Xml.xml Sexpr.t =
                (PXML 
                    -kml (xmlns="http://earth.google.com/kml/2.0")
                        -Document
                            -Name {`String name} --
                            { set (style () :: contents) }
                        --
                    --)


            (* A function to produce a node in the sphere. *)
            type create_node =
                    node_information -> Data.d -> KTree.topology ->
                    Xml.xml -> TemporalGIS.sample option -> 
                        TemporalGIS.sample option ->
                        TemporalGIS.sample option ->TemporalGIS.sample -> 
                            Xml.xml Sexpr.t

            type adjust_tree = KTree.simplified_topology -> KTree.simplified_topology

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

            let plugins : (string, plugin) Hashtbl.t = Hashtbl.create 1667 

            let register_plugin name plugin =
                if Hashtbl.mem plugins name then 
                    failwith "Plugin already exists"
                else Hashtbl.add plugins name plugin


            let get_plugin name =
                try Hashtbl.find plugins name with
                | Not_found -> failwith "Plugin not found"

            let contents adjust_tree data csv tree folder =
                let place_node =
                    match folder.create_node with
                    | None -> place_node
                    | Some x -> x
                in
                let topology = KTree.process data csv tree in
                let topology = 
                    { topology with KTree.topo = adjust_tree topology.KTree.topo }
                in
                let res = ref [] in
                let rec aux parent_gis tree =
                    match tree with
                    | Tree.Parse.Leafp (node, gis) -> 
                            res := 
                                (place_node folder.node_information 
                                data topology node parent_gis None None gis) :: !res
                    | Tree.Parse.Nodep ([l;r] as chld, (node, gis)) ->
                            let l = Some (KTree.extract_gis l)
                            and r = Some (KTree.extract_gis r) in
                            res := 
                                (place_node folder.node_information
                                data topology node parent_gis l r gis) 
                                    :: !res;
                            List.iter (aux (Some gis)) chld
                    | _ -> assert false
                in
                aux None topology.KTree.topo;
                (PXML -Folder -name {`String folder.name} -- { set !res } --)

            let default =
                let folder = 
                    { name = "Tree";
                    node_information = associated_info; 
                    create_node = None }
                in
                { folders = [folder]; 
                 adjust_tree = KTree.adjust_tree; 
                 styles = default_styles }

            let () = register_plugin "default" default

            let generate_kml name plugin data csv trees : Xml.xml Sexpr.t = 
                let cnt = ref 0 in
                let trees =
                    Sexpr.map (fun tree ->
                        incr cnt;
                        let res =
                            List.map (contents plugin.adjust_tree data csv tree) 
                            plugin.folders 
                        in
                        (PXML -Folder -name {`Int !cnt} -- { set res } --)) trees
                in
                base_ml plugin.styles name (Sexpr.to_list trees)

            let has_plugin name = Hashtbl.mem plugins name

            let kml ?(plugin="default") name file data csv trees =
                let ch = open_out file in
                let kml = generate_kml name (get_plugin plugin) data csv trees in
                output_string ch "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
                let () =
                    match kml with
                    | `Single kml -> Xml.to_file ch kml;
                    | _ -> assert false
                in
                close_out ch


        end

    end


    module MainBuild = Build
    module Build = Build.Make (Node) (Edge) (TreeOps)
    module CT = CharTransform.Make (Node) (Edge) (TreeOps)
    module TS = Ptree.Search (Node) (Edge) (TreeOps)
    module PTS = TreeSearch.Make (Node) (Edge) (TreeOps)
    module D = Diagnosis.Make (Node) (Edge) (TreeOps)
    module S = Supports.Make (Node) (Edge) (TreeOps)

    type r = (a, b) run

    let args = 
        IFDEF USEPARALLEL THEN
            Mpi.init Sys.argv
        ELSE
            Sys.argv
        END

    type plugin_function = Methods.script Methods.plugin_arguments -> r -> r

    let plugins : (string, plugin_function) Hashtbl.t = Hashtbl.create 97

    let register_function name f =
        if Hashtbl.mem plugins name then begin
            Status.user_message Status.Error 
                ("There is already a plugin registering the function " ^ name);
            failwith "Function already exists"
        end else
            Hashtbl.add plugins name f

type minimum_spanning_tree = tree 
type build = minimum_spanning_tree list
type minimum_spanning_family = minimum_spanning_tree list
type build_optimum = tree list

let ndebug = true
let ndebug_no_catch = true

(** [reroot_at_outgroup data ptree] reroots [ptree] at the root specified in
    [data].  If the root is not present, the tree will not be rerooted. *)
let reroot_at_outgroup run =
    let reroot_at_outgroup data ptree =
        match data.Data.root_at with
        | None          -> ptree
        | Some outgroup ->
           try let nbr = Ptree.get_parent outgroup ptree in
               let ptree, update =
                   TreeOps.reroot_fn None false (Tree.Edge (outgroup, nbr)) ptree in
               let ptree = TreeOps.incremental_uppass ptree update in
               ptree
           with _ -> ptree
    in
    { run with
        trees = Sexpr.map (reroot_at_outgroup run.data) run.trees }

let report_memory () = 
    (* A function that reports in a formatter-aware style the statistics of the
    * Garbage Collector *)
    let append name value to_string acc = 
        (* Add an element to the table of garbage collection stats *)
        acc ^ "@,@[@{<u>" ^ name ^ "@}: " ^ to_string value ^ "@]" 
    and ( --> ) a b = b a 
    and stat = Gc.stat () in
    "@[<v 2>@{<b>Memory usage:@}@,@[<v>" 
        --> append "Minor Words" stat.Gc.minor_words string_of_float
        --> append "Promoted Words" stat.Gc.promoted_words string_of_float
        --> append "Major Words" stat.Gc.major_words string_of_float
        --> append "Major Collections" stat.Gc.major_collections string_of_int
        --> append "Heap Words" stat.Gc.heap_words string_of_int
        --> append "Heap Chunks" stat.Gc.heap_chunks string_of_int
        --> append "Live Words" stat.Gc.live_words string_of_int
        --> append "Free Words" stat.Gc.free_words string_of_int
        --> append "Free Blocks" stat.Gc.free_blocks string_of_int
        --> append "Largest Free" stat.Gc.largest_free string_of_int
        --> append "Fragments" stat.Gc.fragments string_of_int
        --> append "Compactions" stat.Gc.compactions string_of_int
        --> append "Top Heap Words" stat.Gc.top_heap_words string_of_int 
        --> fun x -> x ^ "@]@]@,%!"


let explode_filenames files =
   List.map (fun x -> `Local x)  (Parser.Wildcard.explode_filenames files)


(* this function differs from update_trees_to_data in that we do not modify the
 * data of the trees, or update the data to the tree (if load_data is false).
 * This is done when cost modes and other settings have been changed, and does
 * not affect the data or how the data is created. *)
let rediagnose_trees ?(classify=true) load_data run =
    let replacer nodes nd = 
        let code = Node.taxon_code nd in
        try All_sets.IntegerMap.find code nodes
        with | Not_found as err ->
            Status.user_message Status.Error ("Couldn't find code " ^ string_of_int code);
            raise err
    and create_node_set ns =
        List.fold_left
            (fun acc x -> All_sets.IntegerMap.add (Node.taxon_code x) x acc)
            All_sets.IntegerMap.empty ns
    in
    let doit st tree =
        let data, nodes =
            if load_data then
                let d,n = Node.load_data ~classify tree.Ptree.data in
                (d, create_node_set n)
            else
                (tree.Ptree.data, tree.Ptree.node_data)
        in
        let tree = { tree with Ptree.node_data = nodes; data = data; } in
        let res = CT.transform_tree (replacer nodes) tree in
        let ach = Status.get_achieved st in
        let () = Status.full_report ~adv:(ach + 1) st in
        res
    in
    (* start it up... *)
    let len = (Sexpr.length run.trees) + (Sexpr.length run.stored_trees) in
    if len > 0 then
        let st = Status.create "Diagnosis"  (Some len) "Recalculating trees" in
        let trees = Sexpr.map (doit st) run.trees in
        let stored = Sexpr.map (doit st) run.stored_trees in
        Status.finished st;
        { run with trees = trees; stored_trees = stored }
    else
        run

let update_trees_to_data ?(classify=true) force load_data run =
    let classify = (not (Data.has_dynamic run.data)) && classify in
    let data, nodes = 
        if load_data then Node.load_data ~classify run.data 
        else run.data, run.nodes 
    in
    let run = { run with nodes = nodes; data = data } in
    let len = (Sexpr.length run.trees) + (Sexpr.length run.stored_trees) in
    if len > 0 then begin
        let st = Status.create "Diagnosis"  (Some len) "Recalculating trees" in
        let nodes =
            List.fold_left
                (fun acc nd ->
                    let code = Node.taxon_code nd in
                    All_sets.IntegerMap.add code nd acc)
                All_sets.IntegerMap.empty run.nodes
        in
        let replacer nd = 
            let code = Node.taxon_code nd in
            try All_sets.IntegerMap.find code nodes with
            | Not_found as err-> 
                    Status.user_message Status.Error 
                    ("Couldn't find code " ^ string_of_int code);
                    raise err
        in
        let doit replacer tree = 
            let are_leaves_different =
                force || (0 <> 
                All_sets.IntegerMap.fold (fun code node acc ->
                    if acc <> 0 then acc
                    else
                        try
                            let n = Ptree.get_node_data code tree in
                            Node.compare node n
                        with Not_found -> 1)
                    nodes 0)
            in
            let ach = Status.get_achieved st in
            if are_leaves_different then
                let tree = { tree with Ptree.node_data = nodes; data = run.data; } in
                let res = CT.transform_tree replacer tree in
                let () = Status.full_report ~adv:(ach + 1) st in
                res
            else 
                let () = Status.full_report ~adv:(ach + 1) st in
                tree
        in
        let trees = Sexpr.map (doit replacer) run.trees in
        let stored = Sexpr.map (doit replacer) run.stored_trees in
        Status.finished st;
        { run with trees = trees; stored_trees = stored }
    end else begin
        run
    end

let process_transform (run : r) (meth : Methods.transform) =
    match meth with
    | `OriginCost float ->
          { run with
                trees =
                  Sexpr.map (Ptree.set_origin_cost float) run.trees }
    (* Since these transforms are tree-dependent, we do not update trees to
       data, and continue on with the diagnosis *)
    | #Methods.tree_transform as meth ->
        let trees, data, nodes =
            CT.transform_nodes_trees run.trees run.data run.nodes [meth]
        in
        {run with trees = trees; nodes = nodes; data = data;}
    | #Methods.char_transform as meth ->
        let data, nodes =
            CT.transform_nodes run.trees run.data run.nodes [meth] 
        in
        update_trees_to_data false false {run with nodes = nodes; data = data}
    | #Methods.terminal_transform as meth ->
            let data, htbl = Data.randomize_taxon_codes meth run.data in
            let data, nodes = Node.load_data data in
            let trees = 
                Sexpr.map 
                (fun x ->
                    { (Ptree.empty data) with 
                        Ptree.tree = Tree.replace_codes 
                        (fun x -> try Hashtbl.find htbl x with _ -> x)
                        x.Ptree.tree;
                })
                run.trees 
            in
            update_trees_to_data true false
            { run with nodes = nodes; data = data; trees = trees }

let load_data (meth : Methods.input) data nodes =
    let prealigned_files = ref [] in
    let rec reader annotated is_prealigned data (meth : Methods.simple_input) = 
        match meth with
        | `Poyfile files ->
            failwith "TODO"
        | `AutoDetect files ->
            let files = explode_filenames files in
            if is_prealigned then
                prealigned_files := files :: !prealigned_files;
            List.fold_left
                (Data.guess_class_and_add_file annotated is_prealigned)
                data files
        | `PartitionedFile files 
        | `Nucleotides files  as meth ->
            let mode = match meth with
                | `PartitionedFile _ -> `Partitioned Data.Clip
                | `Nucleotides _ -> `DO
            in
            let files = explode_filenames files in
            let data =
                List.fold_left
                    (fun acc x -> Data.add_file acc [Data.Characters] x)
                    data files
            in
            if is_prealigned then
                prealigned_files := files :: !prealigned_files;
            List.fold_left
                (fun d f ->
                    Data.process_molecular_file Data.default_tcm
                            Cost_matrix.Two_D.default Cost_matrix.Two_D.default
                            Cost_matrix.Three_D.default annotated
                            Alphabet.nucleotides mode is_prealigned `Seq d f)
                data files
        (** read chromosome data from files each chromosome is presented simply
            as a long plain nucleotide sequences *)
        | `Chromosome files ->
            let files = explode_filenames files in
            let data =
                List.fold_left
                    (fun acc x -> Data.add_file acc [Data.Characters] x)
                    data files
            in
            List.fold_left
                (fun d f ->
                    Data.process_molecular_file Data.default_tcm
                            Cost_matrix.Two_D.default Cost_matrix.Two_D.default
                            Cost_matrix.Three_D.default annotated
                            Alphabet.nucleotides `DO false `Chromosome d f)
                data files
        (** read genome data from files each genome is presented as a sequence
            of chromosomes separated by @ signs *)
        | `Genome files ->
            let files = explode_filenames files in
            let data =
                List.fold_left
                    (fun acc x -> Data.add_file acc [Data.Characters] x)
                    data files
            in
            List.fold_left
                (fun d f ->
                    Data.process_molecular_file Data.default_tcm
                            Cost_matrix.Two_D.default Cost_matrix.Two_D.default
                            Cost_matrix.Three_D.default annotated
                            Alphabet.nucleotides `DO false `Genome d f)
                data files
        | `Aminoacids (files,read_options) ->
            let files = explode_filenames files in
            let data =
                List.fold_left
                    (fun acc x -> Data.add_file acc [Data.Characters] x)
                    data files
            in
            let is_prealigned = List.mem (`Prealigned) read_options in
            if is_prealigned then
                prealigned_files := files :: !prealigned_files;
            let init3D = List.mem (`Init3D true) read_options in
            let alpha =
                if init3D then Alphabet.aminoacids_use_3d
                          else Alphabet.aminoacids
            in
            let dynastate,default_mode =
                if is_prealigned then `SeqPrealigned,`GeneralNonAdd 
                                 else `CustomAlphabet,`DO
            in
            List.fold_left
                (fun d f ->
                    Data.process_molecular_file Data.default_tcm
                        Cost_matrix.Two_D.default_aminoacids Cost_matrix.Two_D.default_aminoacids
                        (Lazy.force Cost_matrix.Three_D.default_aminoacids)
                        annotated alpha default_mode is_prealigned dynastate d f)
                data files
        | `GeneralAlphabetSeq (f, alph, read_options) ->
            let data = Data.add_file data [Data.Characters] f in
            let orientation = List.mem (`Orientation false) read_options in
            let init3D = List.mem (`Init3D true) read_options in
            let is_prealigned = List.mem (`Prealigned) read_options in
            let tie_breaker_first = List.mem (`TieBreaker `First) read_options in
            let tie_breaker_last = List.mem (`TieBreaker `Last) read_options in
            let tie_breaker_random = List.mem (`TieBreaker `Keep_Random) read_options in
            let tb = (*we use First as custom alphabet's tie breaker*)
                if tie_breaker_first then `First
                else if tie_breaker_last then `Last
                else if tie_breaker_random then `Keep_Random
                else `First
            in
            let level = 2 in
            let respect_case = true in
            let alphabet, (twod_full,twod_original,matrix), threed =
                Alphabet.of_file alph orientation init3D level respect_case tb
            in
            if is_prealigned then
                prealigned_files := [f] :: !prealigned_files;
            let dynastate,default_mode =
                if is_prealigned then `SeqPrealigned,`GeneralNonAdd 
                                 else `CustomAlphabet,`DO
            in
            let tcmfile = FileStream.filename alph in
            Data.process_molecular_file ~respect_case:respect_case
                    (Data.Input_file (tcmfile,matrix)) twod_full twod_original
                    threed annotated alphabet default_mode is_prealigned dynastate data f
        (** read breakinv data from files each breakinv is presented as a
            sequence of general alphabets *)
        | `Breakinv (seq, alph, read_options) ->
            let data = Data.add_file data [Data.Characters] seq in
            let orientation = not (List.mem (`Orientation false) read_options) in
            let init3D = (List.mem (`Init3D true) read_options) in
            let data = Data.add_file data [Data.Characters] seq in
            let respect_case = true in
            let alphabet, (twod,twod_ori,matrix), threed =
                Alphabet.of_file alph orientation init3D 0 respect_case `Keep_Random
            and tcmfile = FileStream.filename alph in
            Data.process_molecular_file (Data.Input_file (tcmfile,matrix)) twod
                    twod_ori threed annotated alphabet `DO is_prealigned
                    `Breakinv data seq
        | `ComplexTerminals files ->
            List.fold_left Data.process_complex_terminals data
                           (explode_filenames files)

    and annotated_reader data (meth : Methods.input) =
        match meth with
        | `AnnotatedFiles files ->
            List.fold_left (reader true false) data files
        | #Methods.simple_input as meth -> 
            reader false false data meth
        | `Prealigned (meth, tcm, gap_opening) ->
            prealigned_files := [];
            let data = reader false true data meth in (* read data as dynamic *)
            let data = Data.categorize (Data.remove_taxa_to_ignore data) in
            let files = List.flatten !prealigned_files in
            let chars =
                (* list of characters to filter after conversion to static *)
                let chars = List.rev_map (function `Local x | `Remote x -> (x ^ ".*")) files in
                `Names (true,chars)
            in
            prealigned_files := [];
            let data = match tcm with
                | `Assign_Transformation_Cost_Matrix file ->
                    Data.assign_tcm_to_characters_from_file data chars (Some file) 
                | `Create_Transformation_Cost_Matrix (trans, gaps) ->
                    Data.assign_transformation_gaps data chars trans gaps 
            in
            let data =
                if gap_opening > 0 then
                    Data.assign_affine_gap_cost data chars
                        (Cost_matrix.Affine gap_opening)
                else data
            in
            Data.prealigned_characters ImpliedAlignment.analyze_tcm data chars 
    in
    let data = annotated_reader data meth in   
    let data = Data.categorize (Data.remove_taxa_to_ignore data) in 
    Node.load_data data

type script = Methods.script

let process_input run (meth : Methods.input) =
    let d, nodes = load_data meth run.data run.nodes in
    let run = { run with data = d; nodes = nodes } in
    (* check whether this read any trees *)
    let run = update_trees_to_data false false run in
    let () = List.iter (Data.verify_trees run.data) run.data.Data.trees in
    let trees = List.map (fun (x, _, _) -> x) run.data.Data.trees in
    let trees = Build.prebuilt trees (run.data, run.nodes) in
    let total_trees = Sexpr.union trees run.trees in
    let d = { d with Data.trees = [] } in
    { run with trees = total_trees; data = d }

let temporary_transform run meth =
    let run1 = process_transform run meth in
    match meth with
    | `OriginCost float ->
          let old_origin_cost =
              try let tree = Sexpr.first run.trees in
                  tree.Ptree.origin_cost
              with _ -> 0. in
          run1,
          `UntransformFn (fun run ->
                              { run with
                                    trees =
                                      Sexpr.map
                                          (Ptree.set_origin_cost
                                               old_origin_cost)
                                          run.trees })
    | #Methods.tree_transform
    | #Methods.char_transform     -> run1, `UntransformChar
    | #Methods.terminal_transform -> run1, `UntransformChar

(** [temporary_transforms run meths] transforms [run] as specified by the list
    of transformations [meths].  It returns the transformed [run] and a list of
    functions to undo the transform. *)
let temporary_transforms meths run =
    let replace_chars myrun =
        let trees =
            Sexpr.map_status "Untransforming Characters"
                (fun o_tree ->
                    let dest =
                        Data.convert_static_to_dynamic_branches
                             ~src:o_tree.Ptree.data ~dest:run.data in
                    let data =
                        Data.sync_static_to_dynamic_model
                                 ~src:o_tree.Ptree.data ~dest in
                    CT.replace_nodes_data data o_tree)
                myrun.trees
        in
        { run with trees = trees; }
    in
    match meths with
    | [] -> run, [(fun a -> a)]
    | meths ->
          let rec r ?(chartransform=false) run meths untrans =
              match meths with
              | [] -> run, untrans
              | meth :: meths ->
                    let run, untranstype =
                        temporary_transform run meth in
                    match untranstype with
                    | `UntransformFn fn ->
                          r ~chartransform run meths (fn :: untrans)
                    | `UntransformChar ->
                          let untrans =
                              if chartransform
                              then untrans
                              else replace_chars :: untrans in
                          let chartransform = true in
                          r ~chartransform run meths untrans
          in
          r run meths []


(** [output_clade_file data filename counter (tree, cost)] writes a
    HENNIG86-format file to disk describing every clade found in the tree.  The
    filename derives from the prefix [filename], followed by the integer
    returned by [counter ()], followed by ".hen".  [data] provides the names of
    the OTUs. *)
let output_clade_file data fn counter tree =
    let cost = Ptree.get_cost `Adjusted tree in
    let all_otus = Ptree.get_all_leaves_ids tree in
    (* Code borrowed from Support.bremer_to_sexpr *)
    let fn = fn ^ string_of_int (counter ()) ^ ".hen" in
    let file = open_out fn in
    output_string file
        ("xread\n'Clades file for tree with cost "
         ^ string_of_float cost ^ "'\n");
    let visit (Tree.Edge (h, n)) =
        let otus = Ptree.get_leaves n tree in
        let n_otus = List.length otus in
        if n_otus = 1
        then None
        else
            let otu_set =
                List.fold_left
                    (fun set elt ->
                        All_sets.Integers.add
                            (Node.taxon_code elt) set)
                    All_sets.Integers.empty otus
            in
            Some (otu_set)
    in
    let list =
        let visit edge list = match visit edge with
            | None -> (Tree.Continue, list)
            | Some s -> (Tree.Continue, s :: list) in
        let list_of_handle h acc =
            Ptree.pre_order_edge_visit visit h tree acc in
        let handles = tree.Ptree.tree.Tree.handles in
        All_sets.Integers.fold list_of_handle handles []
    in
    output_string file
        (string_of_int (List.length list) ^ " "
         ^ string_of_int (List.length all_otus) ^ "\n");
    List.iter
        (fun otu ->
             output_string file
                (let code = Node.taxon_code (Ptree.get_node_data otu tree) in
                 try Data.code_taxon code data
                 with | Not_found -> string_of_int code);
             output_string file "\n";
             List.iter
                 (fun set ->
                      output_string file
                          (if All_sets.Integers.mem otu set then "1" else "0"))
                 list;
             output_string file "\n")
        all_otus;
    output_string file ";\ncc-.;\nproc /;\n";
    close_out file


let runtime_store rediagnose run meth =
    let store name run clas = match clas with
        | `Data      -> Hashtbl.replace run.data_store name (run.data, run.nodes)
        | `Trees     -> Hashtbl.replace run.tree_store name run.trees
        | `Bremer    -> Hashtbl.replace run.bremer_store name run.bremer_support
        | `Jackknife -> Hashtbl.replace run.jackknife_store name run.jackknife_support
        | `Bootstrap -> Hashtbl.replace run.bootstrap_store name run.bootstrap_support
    in
    let remove name run clas = match clas with
        | `Data      -> Hashtbl.remove run.data_store name 
        | `Trees     -> Hashtbl.remove run.tree_store name
        | `Bremer    -> Hashtbl.remove run.bremer_store name
        | `Jackknife -> Hashtbl.remove run.jackknife_store name
        | `Bootstrap -> Hashtbl.remove run.bootstrap_store name
    in
    let set name (changed, run) clas =
        let find hstb = Hashtbl.find hstb name in
        match clas with
        | `Data -> 
                let data, nodes = find run.data_store in
                let changed = changed || not (Data.compare run.data data) in
                changed, { run with data = data; nodes = nodes}
        | `Trees ->
                true, { run with trees = find run.tree_store }
        | `Bremer ->
                true, { run with bremer_support = find run.bremer_store }
        | `Jackknife ->
                true, { run with jackknife_support = find run.jackknife_store }
        | `Bootstrap ->
                true, { run with bootstrap_support = find run.bootstrap_store }
    in
    let do_rediagnose = function
        | `Data
        | `Trees -> true
        | _ -> false
    in
    match meth with
    | `Store (clas, name) -> 
            List.iter (store name run) clas ; run
    | `Set (clas, name) ->
            let changed, nrun = 
                try List.fold_left (set name) (false, run) clas with
                | Not_found -> 
                    Status.user_message Status.Error 
                        ("The@ state@ of@ search@ " ^ name ^
                         "@ has@ not@ been@ defined.@ @ will@ continue@ " ^
                         "with@ the@ current@ state.");
                        false, run
            in
            if changed && List.exists do_rediagnose clas then
                rediagnose nrun
            else nrun
    | `Discard (clas, name) -> 
            List.iter (remove name run) clas ; run
    | `Keep_only (num, meth) ->
          Status.user_message Status.Information "Ignoring@ method:@ only@ keeping@ best";
          let trees =
              List.fast_sort (fun a b -> compare (Ptree.get_cost `Adjusted a)
                                  (Ptree.get_cost `Adjusted b))
                  (Sexpr.to_list run.trees) in
          let rec select i list =
              if i >= num
              then []
              else match list with
              | [] -> []
              | l :: ls -> l :: (select (succ i) ls) 
          in
          { run with trees = Sexpr.of_list (select 0 trees) }

let create_hstb () = Hashtbl.create 13 

let empty () = {
    description = None;
    trees = `Empty; 
    data = Data.empty ();
    nodes = [];
    bremer_support  = `Empty;
    jackknife_support = None;
    bootstrap_support = None;
    runtime_store = create_hstb ();
    data_store = create_hstb ();
    bremer_store = create_hstb ();
    jackknife_store = create_hstb ();
    bootstrap_store = create_hstb ();
    tree_store = create_hstb ();
    queue = Sampler.create ();
    stored_trees = `Empty;
    original_trees = `Empty;
    search_results = empty_search_results;
}

let process_random_seed_set run v =
    let sv = string_of_int v in
    let msg = "@[Setting@ random@ seed@ value@ to@ " ^ sv ^ "@]"in
    Status.user_message Status.Information msg;
    Sadman.start "random_seed" [("process", "0"); ("value", sv)];
    Sadman.finish [];
    Random.init v; 
    run

let process_sets run name ident =
    match ident with
    | `Codon ident ->
        let msg = "@[Creating@ alias@ and@ character@ set@ for@ codons@]" in
        Status.user_message Status.Information msg;
        { run with data = Data.make_codon_partitions false run.data name ident }
    | `Chars ident ->
        let msg = "@[Creating@ alias@ and@ character@ set@ "^name^"@]"in
        Status.user_message Status.Information msg;
        { run with data = Data.make_set_partitions false run.data name ident }

let do_recovery run =
    let trees = Sexpr.to_list run.trees in
    let rec stackadder () = 
        if Stack.is_empty run.queue.Sampler.stack then ()
        else 
            let _ = 
                Queue.push (Stack.pop run.queue.Sampler.stack) 
                run.queue.Sampler.queue 
            in
            stackadder ()
    in
    let rec adder res = 
        if Queue.is_empty run.queue.Sampler.queue then Sexpr.of_list res
        else adder (let trees = List.map (fun (a, _, _) -> 
            let nt = { (Ptree.empty run.data) with Ptree.tree = a } in
            let nt = List.fold_left (fun acc x ->
                Ptree.add_node_data (Node.taxon_code x) x acc) nt run.nodes
            in
            CT.replace_nodes run.nodes nt)
            (Queue.pop run.queue.Sampler.queue) in trees @ res)
    in
    stackadder ();
    { run with trees = adder trees }
                     
let process_characters_handling (run : r) meth = 
    let data, do_nodes = match meth with
        | `RenameCharacters syns ->
                let process a ((_, y) as b) =
                    try Data.process_rename_characters a b with
                    | Data.Illegal_argument ->
                            let msg = "The@ name " ^ StatusCommon.escape y ^ "@ is@ already@ used" ^
                            "@ to@ describe@ a@ character@ in@ this@ dataset." ^
                            "@ I@ will@ ignore@ it@ as@ a@ synonym." in
                            Status.user_message Status.Error msg;
                            a
                in
                List.fold_left process run.data syns, false
        | `AnalyzeOnlyCharacterFiles (dont_complement, files) ->
                Data.process_analyze_only_characters_file true dont_complement run.data files, 
                true
        | `AnalyzeOnlyCharacters c ->
                let c = Data.get_chars_codes_comp run.data c in
                let c = Data.complement_characters run.data (`Some c) in
                Data.process_ignore_characters true run.data (c :> Data.characters),
                true
    in
    if do_nodes then
        let data, nodes = Node.load_data data in
        { run with nodes = nodes; data = data }
    else  { run with data = data }

let ( --> ) a b = b a 

let process_taxon_filter (run : r) meth =
    if 0 < Sexpr.length run.trees then begin
        Status.user_message Status.Error
            ("Selecting@ a@ subset@ of@ taxa@ while@ trees@ " ^
             "are@ in@ memory@ is@ not@ permitted.");
        failwith "Illegal request of select"
    end else
        let data = match meth with
            | `IgnoreTaxaFiles files ->
                files --> List.fold_left Data.process_ignore_file run.data
                      --> Data.remove_taxa_to_ignore
            | `IgnoreTaxa taxa ->
                taxa --> List.fold_left Data.process_ignore_taxon run.data
                     --> Data.remove_taxa_to_ignore
            | `SynonymsFile files ->
                List.fold_left Data.add_synonyms_file run.data files
            | `Synonyms taxa ->
                List.fold_left Data.add_synonym run.data taxa
            | `AnalyzeOnlyFiles (dont_complement, files) ->
                files --> Data.process_analyze_only_file dont_complement run.data
                      --> Data.remove_taxa_to_ignore
            | `AnalyzeOnly c ->
                run.data --> Data.process_analyze_only_taxa c
                         --> Data.remove_taxa_to_ignore
            | `Random _ as c ->
                run.data --> Data.process_analyze_only_taxa c
                         --> Data.remove_taxa_to_ignore
        in
        let data = Data.categorize data in
        let data, nodes = Node.load_data data in
        let () =
            Printf.ksprintf
                (Status.user_message Status.Status)
                ("Selected %d of %d terminals")
                (List.length nodes)  (run.data.Data.number_of_taxa)
        in
        { run with nodes = nodes; data = data }

let sort_trees trees =
    let t = Sexpr.to_list trees in
    let comparison a b =
        compare (Ptree.get_cost `Adjusted a) (Ptree.get_cost `Adjusted b)
    in
    List.sort comparison t

let rec process_tree_handling run meth =
    let rec get_first_n n lst =
        if n < 1 then []
        else match lst with
            | h :: t -> h :: (get_first_n (n - 1) t)
            | [] -> []
    in
    let get_optimal = function
        | (h :: _) as lst ->
            let cost = Ptree.get_cost `Adjusted h in
            List.filter (fun x -> cost = Ptree.get_cost `Adjusted x) lst
        | [] -> []
    in
    let trees = match meth with
        | `BestN (Some n) ->
                get_first_n n (sort_trees run.trees)
        | `BestN None ->
                let trees = sort_trees run.trees in
                get_optimal trees
        | `BestWithin th ->
                begin match sort_trees run.trees with
                    | h :: t ->
                        let cost = th +. (Ptree.get_cost `Adjusted h) in
                        let t =
                            List.filter (fun x -> (Ptree.get_cost `Adjusted x) < cost) t
                        in
                        (h :: t)
                    | [] -> []
                end
        | `RandomTrees n ->
                let lst = Sexpr.to_list run.trees in
                let arr = Array.of_list lst in
                Array_ops.randomize arr;
                let lst = Array.to_list arr in
                get_first_n n lst
        | `Unique ->
                let sexpr = Sexpr.to_list run.trees in
                let res =
                    try TS.get_unique sexpr with
                    | Not_found as err ->
                        Status.user_message Status.Error "This is the path";
                        raise err
                in
                res
    in
    let trees = Sexpr.of_list trees in
    { run with trees = trees }

exception Error_in_Script of (exn * r)

let check_ft_queue run =
    if Stack.is_empty run.queue.Sampler.stack then ()
    else begin
        while not (Stack.is_empty run.queue.Sampler.stack) do
            Queue.push (Stack.pop run.queue.Sampler.stack)
            run.queue.Sampler.queue;
        done;
    end

let explode_trees run =
    Sexpr.map (fun x -> { run with trees = `Single x }) run.trees

let is_forest = function
    | `LocalOptimum x -> x.Methods.oo

let build_has_exact = build_has `Exact

let only_multistatic meth =
    List.fold_left (fun acc x ->
        acc & (match x with
        | `Static_Aprox _
        | `Automatic_Static_Aprox _ -> false
        | _ -> true)) true meth

let warn_if_no_trees_in_memory trees = 
    let items = Sexpr.length trees in
    if items = 0 then
        Status.user_message Status.Warning
            ("There@ are@ no@ active@ trees@ in@ memory!")

let get_trees_for_support support_class run =
    let do_support support_set x = match support_set with
        | None -> `Empty
        | Some (iterations, fs) ->
            begin match x with
                | `Individual ->
                    Sexpr.map 
                        (fun tree ->
                            Ptree.supports 
                                (fun x -> Data.code_taxon x run.data)
                                0.0 (float_of_int iterations) tree.Ptree.tree fs)
                        run.trees
                | `InputFile x ->
                    let trees = Tree.Parse.of_file (`Local x) in
                    Sexpr.of_list
                        (List.map 
                            (fun tree ->
                                Ptree.support_of_input
                                    (fun x -> Data.code_taxon x run.data)
                                    0.0 (float_of_int iterations) tree run.data fs)
                            trees)
                | `Consensus ->
                    `Single 
                        (Ptree.preprocessed_consensus 
                            (fun code -> Data.code_taxon code run.data) 
                            ((float_of_int iterations) /. 2.0)
                            iterations
                            fs)
            end
    in
    match support_class with
    | `Bremer (Some input_files) ->
            S.bremer_of_input_file_but_trust_input_cost
                (match run.data.Data.root_at with
                    | Some x -> x
                    | None ->
                        Status.user_message Status.Error
                            ("I@ require@ a@ specified@ root.@ See@ the@ "^
                             "documentation@ for@ the@ set(root:ID)@ command.");
                            raise Not_found)
                (fun x -> Data.code_taxon x run.data)
                run.data
                input_files
                run.trees,
            "Bremer"
    | `Bremer None ->
            Sexpr.map (S.support_to_string_tree run.data) run.bremer_support,
            "Bremer"
    | `Jackknife x ->
            let res = do_support run.jackknife_support x in
            res, "Jackknife"
    | `Bootstrap x ->
            let res = do_support run.bootstrap_support x in
            res, "Bootstrap"

(** PAUP-like output for site likelihood *)
let report_lk_sites ft tree chars : unit =
    let tree_num = ref 0 in
    let header = [| "Tree"; "-lnL"; "Site"; "Weight"; "-lnL"; |] in
    let modify_array (cost,array) =
        let second = [| string_of_int !tree_num; string_of_float cost; ""; ""; ""; |] in
        incr tree_num;
        Array.init
            ((Array.length array)+2)
            (fun i ->
                if i = 0 then header else
                if i = 1 then second
                else begin
                    let a,b = array.(i-2) in
                    [| ""; "";string_of_int (i-2); string_of_float a; string_of_float b; |]
                end)
    in
    match Ptree.get_roots tree with
    | [x] ->
        let node_root = match x.Ptree.root_median with
            | Some (_,n) -> n
            | None       -> assert false
        in
        List.iter (fun x -> let x = modify_array x in ft x)
                  (Node.get_lk_sites node_root chars)
    |  _  ->
        Status.user_message Status.Error
            ("I@ am@ unsure@ what@ to@ do@ about@ reporting@ site@ likelihood@ "
            ^"on@ a@ disjoint@ tree.@ I@ am@ skipping@ reporting@ its@ data")


let rec handle_support_output run meth = match meth with
    | `Supports (support_class, filename) ->
        begin match support_class with
            | Some support_class ->
                let trees, title = get_trees_for_support support_class run in
                (* If we have bremer supports, we print them *)
                let fo = Status.Output (filename, false, []) in
                Status.user_message fo ("@[<v>@[" ^ title ^ "@ Supports:@]@,");
                let output tree =
                    Status.user_message fo ("@[Support@ tree:@]@,@[");
                    Status.user_message fo (AsciiTree.for_formatter false true false tree);
                    Status.user_message fo "@]";
                in
                Sexpr.leaf_iter output trees;
                Status.user_message fo "@,%!@]";
            | None ->
                let a = Some (`Bremer None)
                and b = Some (`Jackknife `Individual)
                and c = Some (`Bootstrap `Individual) in
                handle_support_output run (`Supports (a, filename));
                handle_support_output run (`Supports (b, filename));
                handle_support_output run (`Supports (c, filename));
        end
    | `GraphicSupports (support_class, filename) ->
        begin match support_class with
            | Some support_class ->
                let trees, title = get_trees_for_support support_class run in
                let trees = Sexpr.to_list trees in
                let trees = Array.of_list trees in
                if 0 = Array.length trees then ()
                else begin
                    let trees = Array.map (fun x -> 0.0, x) trees in
                    match filename with
                        | Some f -> GraphicsPs.display title f trees
                        | None   ->
                            let fo = Status.Output (filename,false, []) in
                            Array.iter
                                (fun (_, x) ->
                                    let r = AsciiTree.to_string ~sep:2 ~bd:6 true x in
                                    Status.user_message fo ("@[@[<v>@[" ^ title ^ "@ Support@ tree@]@,@[");
                                    Status.user_message fo r;
                                    Status.user_message fo "@,@]@]@]%!")
                                trees;
                end
            | None ->
                let a = Some (`Bremer None)
                and b = Some (`Jackknife `Individual)
                and c = Some (`Bootstrap `Individual) in
                handle_support_output run (`GraphicSupports (a, filename));
                handle_support_output run (`GraphicSupports (b, filename));
                handle_support_output run (`GraphicSupports (c, filename));
        end

let update_mergingscript folder mergingscript run tmp =
    let tmp = { run with trees = Sexpr.union tmp.trees
    run.stored_trees } in
    let tmp = List.fold_left folder tmp mergingscript in
    let tmp = { tmp with trees = `Empty; stored_trees = tmp.trees } in
    tmp

let on_each_tree folder set_data dosomething mergingscript run tree =
    let tmp = { run with trees = `Empty } in
    let tmp = folder tmp set_data in
    let tmp = { tmp with trees = `Single tree } in
    let tmp = List.fold_left folder tmp dosomething in
    update_mergingscript folder mergingscript run tmp

let emit_identifier =
    let identifier = ref (-1) in
    fun () -> 
        incr identifier; 
        "__poy_id_" ^ string_of_int !identifier


    let extract_tree tree = tree.Ptree.tree

    (*
    let force tree = 
        let force_root x =
            { x with Ptree.root_median =
                match x.Ptree.root_median with
                | None -> None
                | Some (a, b) -> Some (a, Node.force b) }
        in
        { tree with Ptree.node_data = 
            All_sets.IntegerMap.map Node.force tree.Ptree.node_data;
            edge_data = Tree.EdgeMap.map Edge.force tree.Ptree.edge_data;
            component_root = All_sets.IntegerMap.map force_root
            tree.Ptree.component_root}
        *)

    let encode_trees print_msg run =
        (*
        run.trees, run.stored_trees
        *)
        print_msg "Will encode the trees, first test";
            let () =
                let _ = try TS.get_unique (Sexpr.to_list run.trees) with
                    | err -> 
                            Status.user_message Status.Error "12";
                            raise err
                in
                ()
            in
            let () =
                let _ = try TS.get_unique (Sexpr.to_list run.stored_trees) with
                    | err -> 
                            Status.user_message Status.Error "11";
                            raise err
                in
                ()
            in
        print_msg "The test passed";
        (Sexpr.map extract_tree run.trees), (Sexpr.map extract_tree run.stored_trees)

    let toptree data tree = 
        { (Ptree.empty data) with Ptree.tree = tree } 

    let update_codes run tc tree = 
        let code_ref = ref run.data.Data.number_of_taxa in
        let my_hsh = Hashtbl.create 91 in
        let update_item a = 
            if All_sets.IntegerMap.mem a tc then
                Data.taxon_code (All_sets.IntegerMap.find a tc) run.data
            else if Hashtbl.mem my_hsh a then
                Hashtbl.find my_hsh a
            else 
                let _ = incr code_ref in
                let newa = !code_ref in
                let _ = Hashtbl.add my_hsh a newa in
                newa
        in
        let update_node code node acc = 
            let node = 
                match node with
                | Tree.Interior (a, b, c, d) ->
                        let a = update_item a
                        and b = update_item b
                        and c = update_item c
                        and d = update_item d in
                        Tree.Interior (a, b, c, d)
                | Tree.Leaf (a, b) ->
                        let a = update_item a
                        and b = update_item b in
                        Tree.Leaf (a, b)
                | Tree.Single a ->
                        let a = update_item a in
                        Tree.Single a
            in
            All_sets.IntegerMap.add (update_item code) node acc
        in
        let update_edge (Tree.Edge (a, b)) acc =
            Tree.EdgeSet.add (Tree.Edge (update_item a, update_item b)) acc
        in
        let u_topo = All_sets.IntegerMap.fold update_node tree.Tree.u_topo
        All_sets.IntegerMap.empty
        and d_edges = 
            Tree.EdgeSet.fold update_edge tree.Tree.d_edges
            Tree.EdgeSet.empty
        and handles = 
            All_sets.Integers.fold 
            (fun x acc -> All_sets.Integers.add (update_item x) acc) 
            tree.Tree.handles All_sets.Integers.empty
        in
        { 
            Tree.tree_name = None;
                u_topo = u_topo; 
                d_edges = d_edges;
                handles = handles;
                avail_ids = [];
                new_ids = run.data.Data.number_of_taxa + 1; 
        }

    let decode_trees print_msg (trees, stored_trees) run =
        (*
        { run with trees = Sexpr.union trees run.trees; stored_trees =
            Sexpr.union stored_trees run.stored_trees }
        let trees = Sexpr.map (update_codes run tc) trees in  
        *)
        try
            print_msg "Will attempt to decode trees";
            let my_trees = run.trees 
            and my_stored = run.stored_trees in
            let nrun = 
                let nt = Sexpr.map (toptree run.data) trees 
                and nst = Sexpr.map (toptree run.data) stored_trees in
                update_trees_to_data true false { run with trees = nt; stored_trees = nst }
            in
            let () =
                let _ = try TS.get_unique (Sexpr.to_list my_trees) with
                    | Not_found as err -> 
                            Status.user_message Status.Error "9";
                            raise err
                in
                ()
            in
            let () =
                let _ = try TS.get_unique (Sexpr.to_list my_stored) with
                    | Not_found as err -> 
                            Status.user_message Status.Error "10";
                            raise err
                in
                ()
            in
            let () =
                let _ = try TS.get_unique (Sexpr.to_list nrun.stored_trees) with
                    | Not_found as err -> 
                            Status.user_message Status.Error "8";
                            raise err
                in
                ()
            in
            let () =
                let _ = try TS.get_unique (Sexpr.to_list nrun.trees) with
                    | Not_found as err -> 
                            Status.user_message Status.Error "7";
                            raise err
                in
                ()
            in
            print_msg "Did the tree decoding";
            { run with trees = Sexpr.union my_trees nrun.trees; stored_trees =
                Sexpr.union my_stored nrun.stored_trees }
        with
        | Not_found as err -> 
                print_msg "This occurs during the decoding";
                raise err

    let encode_jackknife run = 
        (run.jackknife_support), (run.data.Data.taxon_codes)

    let to_my_code tc run code =
        Data.taxon_code (All_sets.IntegerMap.find code tc) run.data

    let decode_jackknife set run = 
        let set, tc = set in
        match set with
        | Some (int, set) ->
            let to_my_code code acc = 
                All_sets.Integers.add (to_my_code tc run code) acc
            in
            Some (int,
            Tree.CladeFPMap.fold (fun clades count acc ->
                let clades = 
                    All_sets.Integers.fold to_my_code clades All_sets.Integers.empty
                in
                Tree.CladeFPMap.add clades count acc) set Tree.CladeFPMap.empty)
        | None -> None

    let encode_bootstrap run = run.bootstrap_support, run.data.Data.taxon_codes

    let decode_bootstrap set run = 
        decode_jackknife set run

    let encode_bremer run = run.bremer_support, run.data.Data.taxon_codes

    let decode_bremer set run = set


IFDEF USEPARALLEL THEN
    let filter_my_trees run =
        let rank = (Mpi.comm_rank Mpi.comm_world) 
        and total = (Mpi.comm_size Mpi.comm_world) in
        let _, trees =
            Sexpr.fold_left (fun (pos, trees) tree ->
                if rank = pos mod total then (pos + 1), (tree :: trees)
                else (pos + 1), trees) (0, []) run.trees
        in
        { run with trees = Sexpr.of_list trees }
END

    let add_something get adder run jck =
        let njck =
            match get run, jck with
            | x, None
            | None, x -> x
            | Some (ont, ocnts), Some (nt, cnts) -> 
                    let new_total = ont + nt 
                    and new_clades = 
                        Tree.CladeFPMap.fold (fun clade cnt res ->
                            if Tree.CladeFPMap.mem clade res then
                                let res_cnt = Tree.CladeFPMap.find clade res in
                                Tree.CladeFPMap.add clade (res_cnt + cnt) res
                            else Tree.CladeFPMap.add clade cnt res)
                        cnts
                        ocnts 
                    in
                    Some (new_total, new_clades)
        in
        adder run njck

    let add_jackknifes run bts =
        add_something (fun r -> r.jackknife_support) 
        (fun r x -> { r with jackknife_support = x })
        run bts

    let add_boostraps run bts =
        add_something (fun r -> r.bootstrap_support) 
        (fun r x -> { r with bootstrap_support = x })
        run bts

    let add_bremers run (bms, tc) =
        (* Both should have the same number of trees *)
        assert ((Sexpr.length run.bremer_support) = 
            (Sexpr.length bms));
        let prev = Sexpr.to_list run.bremer_support
        and newb = Sexpr.to_list bms in
        let rec join2 t1 t2 =
            match t1, t2 with
            | Methods.Leaf i, Methods.Leaf j when i = to_my_code tc run j -> Methods.Leaf i
            | Methods.Node (w1, t1l, t1r), Methods.Node (w2, t2l, t2r) ->
                    Methods.Node (min w1 w2, join2 t1l t2l, join2 t1r t2r)
            | _ -> failwith "join_support_trees/join2: incompatible trees"
        in
        let res = List.map2 join2 prev newb in
        { run with bremer_support = Sexpr.of_list res }

IFDEF USEPARALLEL THEN
    let () = 
        let my_rank = Mpi.comm_rank Mpi.comm_world in
        let vbst = Mpi.broadcast Methods.Low 0 Mpi.comm_world in
        match my_rank with
        | 0 -> ()
        | _ -> 
                let printer_function t m =
                    let vbst = 
                        match t with
                        | Status.Output _ -> Methods.High
                        | _ -> vbst
                    in
                    match vbst with
                    | Methods.High -> 
                            Mpi.send (t, m) 
                            0
                            Methods.io
                            Mpi.comm_world
                    | _ -> ()
                in
                Status.is_parallel my_rank (Some printer_function)

    let debug_parallel = false
    let print_msg msg = 
        if debug_parallel then begin
            let my_rank = Mpi.comm_rank Mpi.comm_world in
            print_endline (string_of_int my_rank ^ ":" ^ msg);
            flush stdout
        end else ()
END

let automated_search folder max_time min_time max_memory min_hits target_cost visited user_constraint run =
    let module FPSet = Set.Make (Ptree.Fingerprint) in
    let has_dynamic = Data.has_dynamic run.data in
    let timer = Timer.start () in
    let search_results = ref empty_search_results in
    let get_memory () = (Gc.stat ()).Gc.heap_words in
    let command_processor = PoyCommand.of_parsed false in
    let trees = ref run.trees in
    let exhausted = ref FPSet.empty in
    let ratchets = ref 0 in
    let best_cost = ref (match target_cost with | None -> max_float | Some x -> x)
    and hits = ref 0 
    and memory_change = ref (get_memory ())
    and memory_limit = ref false 
    and run = ref { run with trees = `Empty } in
    let st = Status.create "Automated Search"  None "" in 
    let remaining_time () =
        max_time -. (Timer.wall timer)
    in
    let update_information step =
        begin match step with
        | `Initial run ->
                (* We only update the number of hits and the memory change *)
                let cost = Ptree.get_cost `Adjusted (Sexpr.first run.trees) in
                search_results := add_results cost !search_results;
                if cost < !best_cost then begin 
                    best_cost := cost;
                    hits := 1;
                end else if cost = !best_cost then incr hits
                else ();
                let current_memory = get_memory () in
                let memory_delta = current_memory - !memory_change in
                if current_memory + memory_delta > max_memory then 
                    memory_limit := true
                else ();
                memory_change := current_memory;
        | `Others (before, run) ->
                let cost, hit = 
                    let get_cost run =
                        Sexpr.fold_left (fun (cost, hits) x ->
                            let nc = Ptree.get_cost `Adjusted x in
                            if nc = cost then (cost, hits + 1) 
                            else if nc < cost then (nc, 1)
                            else (cost, hits)) (max_float, 0) run.trees
                    in
                    let cost_b, hit_b = get_cost before 
                    and cost_a, hit_a = get_cost run in
                    if cost_b < cost_a then cost_b, 0
                    else
                        cost_a, (if cost_b = cost_a then hit_a - hit_b else hit_a)
                in
                search_results := add_results cost !search_results;
                best_cost := cost;
                hits := !hits + hit;
        end;
        Status.full_report ~msg:("Best tree: " ^ string_of_float
        !best_cost ^ "; Time left: " ^ Timer.to_string (remaining_time ())
        ^"; Hits: " ^ string_of_int !hits) st;
    in
    let select_if_necessary () =
        if !memory_limit then
            let len = Sexpr.length !trees in
            let com = 
                command_processor (CPOY select (best:[max (len / 2) 1])) 
            in
            let () = 
                Gc.full_major ();
                memory_limit := false;
                let nrun = folder { !run with trees = !trees } com in
                run := { nrun with trees = `Empty };
                trees := nrun.trees;
            in 
            ()
        else ()
    in
    let did_initial = ref false in
    let stop_if_necessary state =
        let time = Timer.wall timer in
        if time >= max_time then raise Exit;
        let () =
            IFDEF USEPARALLEL THEN
                ()
            ELSE
                let time = Timer.wall timer in
                if !hits >= min_hits && (time >= min_time) then
                    raise Exit
                else if !hits >= min_hits && (min_time = max_time) then 
                    raise Exit
                else ()
            END
        in
        let fraction =
            match state with
            | `Initial when not !did_initial -> 0.66
            | `Perturb -> 
                    did_initial := true;
                    0.66
            | _ -> 1.0
        in
        if time >= (fraction *. max_time) then raise Exit 
        else 
            match state with
            | `Initial -> if !hits >= (min_hits / 2) then raise Exit
            | _ -> ()
    in
    let get_cost nrun =
        Ptree.get_cost `Adjusted (Sexpr.first nrun.trees)
    in
    let iterations_counter = ref 0 in
    let fuse_counter = ref 0 in
    let collect_results not_final =
IFDEF USEPARALLEL THEN
        let my_trees = !run.trees in
        let run = 
            if not_final then
                folder !run [`BestN None; `Barrier]
            else folder !run [`Barrier]
        in
        let pot_mem =
            try
                ((get_memory ()) * (Mpi.comm_size Mpi.comm_world)) 
            with
            | _ -> max_int
        in
            let arr = Mpi.allgather 
                (!iterations_counter, !ratchets, !fuse_counter, !hits, !best_cost, pot_mem)
                Mpi.comm_world in
            let iterations_counter, ratchet_counter, fuse_counter, hits, best, mem = 
                Array.fold_left (fun (aic, ars, afs, ahts, abest, ax) (ic, rs, fs, hts, best, x) -> 
                    aic + ic, ars + rs, afs + fs,
                        (if best < abest then hts
                        else if best = abest then hts + ahts
                        else ahts), min best abest, max ax x) 
                    (0, 0, 0, 0, max_float, 0) arr in
            let run =
                if not_final then
                    folder run [`GatherTrees ([`BestN None], [])]
                else if mem < max_memory then
                    folder run [`GatherTrees ([`Unique], [])]
                else
                    let arr = Mpi.allgather (Sexpr.length run.trees) Mpi.comm_world in
                    let mmax = Array.fold_left max 0 arr in
                    let mmax = max 1 (mmax / 2) in
                    let run = folder run [`Unique; `BestN (Some mmax)] in
                    folder run [`GatherTrees ([`Unique; `BestN (Some mmax)],[])]
            in
            let run =
                if not_final then
                    { run with trees = Sexpr.union run.trees my_trees }
                else
                    run
            in
            let search_results =
                !search_results
                    --> incr_search_results `Builds iterations_counter
                    --> incr_search_results `Fuse fuse_counter
                    --> incr_search_results `Ratchet ratchet_counter
            in
            { run with search_results = search_results },
            iterations_counter, fuse_counter, hits, best
ELSE
        let search_results =
            !search_results
                --> incr_search_results `Builds !iterations_counter
                --> incr_search_results `Fuse !fuse_counter
                --> incr_search_results `Ratchet !ratchets
        in
        { !run with search_results = search_results }, 
        !iterations_counter, !fuse_counter, !hits, !best_cost
END
    in
    let exec r c = folder r (command_processor c) in
    let add_visited acc = 
        match visited with
        | Some (Some file) -> (APOY visited:[file]) :: acc
        | Some None -> (APOY visited) :: acc
        | None -> acc
    in
    let add_constraint acc =
        match user_constraint with
        | None -> acc
        | Some file -> 
                Printf.printf "Adding constraint\n%!";
                (APOY constraint_p:(file:[file])) :: acc
    in
    let get_top_n n =
        let trees = sort_trees !trees in
        let rec get_n n tl acc =
            if n > 0 then 
                match tl with
                | h :: t -> get_n (n - 1) t (h :: acc)
                | [] -> acc
            else acc
        in
        get_n n trees []
    in
    let add_constraint_iterations condition trees nrun acc =
        match user_constraint with
        | None when condition ->
                let trees = get_top_n trees in
                let trees = Sexpr.of_list trees in
                let sts = 
                    TreeSearch.sets
                    (`Partition [])
                    nrun.data
                    (Sexpr.union trees nrun.trees)
                in
                (APOY sets:[sts]) :: acc
        | _ -> add_constraint acc
    in
    let add_all acc = 
        if has_dynamic then (APOY all) :: acc
        else acc
    in
    let add_randomized acc =
        if Random.bool () then ((APOY randomized) :: acc) else acc
    in
    try
    while true do
    let search_iteration_status = Status.create "RAS + TBR" None "" in
    Status.full_report search_iteration_status;
    try
        (* We run this command in chunks of 10 trees to try to maintain a 
           set of trees and run fusing on them often but not too much *)
        for i = 1 to 10 do
            incr iterations_counter;
            Status.message search_iteration_status ("Searching on tree number " ^ string_of_int !iterations_counter);
            Status.full_report search_iteration_status;
            let nrun =
                let initial = match user_constraint with
                    | None      -> [APOY trees:1]
                    | Some file -> (`Constraint (Some file)) :: [APOY trees:1]
                in
                exec !run (CPOY build {initial})
            in
            let build_cost = get_cost nrun in
            let prev_time = Timer.wall timer in
            let nrun =
                let command =
                    let initial =
                        [APOY tbr; APOY timeout:[`Dynamic remaining_time]]
                    in
                    CPOY swap { initial --> add_constraint --> add_visited }
                in
                exec nrun command
            in
            let search_time = (Timer.wall timer) -. prev_time in
            let do_perturb = build_cost = get_cost nrun in
            let prev_meth = !Methods.cost in
            let nrun, do_perturb =
                if do_perturb && (0.3 < Random.float 1.) then
                    nrun, do_perturb
                else
                    if not has_dynamic then nrun, true
                    else
                        let initial =
                            let time () =
                                min (remaining_time ()) (search_time /. 2.)
                            in
                            [APOY timeout:[`Dynamic time]; APOY tbr]
                        in
                        try 
                            let args = 
                                if (0.5 > Random.float 1.) &&
                                    ((prev_meth = `Normal_plus_Vitamines) || (prev_meth = `Normal)) then
                                        Methods.cost := `Exhaustive_Weak;
                                initial
                                    --> add_constraint_iterations (!iterations_counter >= 4) 4 nrun
                                    --> add_visited
                                    --> add_all
                            in
                            let nrun = exec nrun (CPOY swap {args}) in
                            nrun, false
                        with
                        | err -> Methods.cost := prev_meth; raise err
            in
            Methods.cost := prev_meth;
            trees := Sexpr.union nrun.trees !trees;
            update_information (`Initial nrun);
            if do_perturb || 0.5 < Random.float 1.0 then begin
                let comparison a b = 
                    let gc = Ptree.get_cost `Adjusted in
                    if (i mod 2) = 0 then
                        compare (gc a) (gc b)
                    else
                        let cost_comparison = 
                            if (gc a) < (gc b) then
                                if 0.66 < Random.float 1. then (-1)
                                else 1
                            else if (gc a) > (gc b) then
                                if 0.66 < Random.float 1. then 1
                                else (-1)
                            else 0
                    in
                    match cost_comparison with
                    | 0 -> if Random.bool () then 1 else -1
                    | x -> x
                in
                if 0. < remaining_time () then incr ratchets;
                (* We need to pick one tree from all those in memory *)
                let nrun = 
                    let arr = !trees --> Sexpr.to_list --> Array.of_list in
                    let pos =
                        if (i mod 2) = 0 && 1 < Array.length arr then
                            (* One of the top 5 *)
                            Random.int (min 5 (Array.length arr))
                        else 0
                    in
                    Array_ops.randomize arr;
                    Array.stable_sort comparison arr;
                    { nrun with trees = `Single arr.(pos) }
                in
                let time () =
                    let r = remaining_time () in
                    if has_dynamic then r
                    else min r (search_time /. 2.)
                in
                let swap_args = 
                    [APOY tbr; APOY timeout:[`Dynamic time]] --> add_constraint
                in
                let nrun = 
                    exec nrun 
                        (if has_dynamic then
                            (CPOY 
                                perturb (iterations:4, transform (tcm:(1,1), static_approx), timeout:[`Dynamic remaining_time], swap { swap_args }))
                        else
                            (CPOY 
                                perturb (iterations:4, timeout:[`Dynamic remaining_time], swap { swap_args }) swap (timeout:[`Dynamic remaining_time], randomized)))
                in
                trees := Sexpr.union nrun.trees !trees;
                update_information (`Initial nrun);
            end;
            select_if_necessary ();
            stop_if_necessary `Perturb;
        done;
        raise Exit
    with
    | Exit ->
        Status.finished search_iteration_status;
        let do_exhaustive_round r =
            if not has_dynamic then r
            else begin
                let r' = r in
                let compare_trees a b =
                    let cost = Ptree.get_cost `Adjusted in
                    compare (cost a) (cost b)
                in
                let prev = !Methods.cost in
                match prev with
                | `Exhaustive_Weak | `Exhaustive_Strong | `Iterative _ -> r
                | `Normal | `Normal_plus_Vitamines ->
                    try
                        let potential, not_potential = 
                            Sexpr.split
                                (fun x ->
                                    not (FPSet.mem (Ptree.Fingerprint.fingerprint x) !exhausted))
                                r.trees 
                        in
                        match List.sort compare_trees (Sexpr.to_list potential) with
                        | [] -> r
                        | h :: t ->
                            let r = { r with trees = `Single h } in
                            Methods.cost := `Exhaustive_Weak;
                            let r = 
                                let args = 
                                    [APOY timeout:[`Dynamic remaining_time]] -->
                                        add_visited --> add_constraint
                                in
                                exec r (CPOY swap {args})
                            in
                            let h = Sexpr.first r.trees in
                            exhausted := FPSet.add (Ptree.Fingerprint.fingerprint h) !exhausted;
                            Methods.cost := prev;
                            let r =
                                { r with trees = Sexpr.union r.trees (Sexpr.union not_potential (Sexpr.of_list t)) }
                            in
                            update_information (`Others (r', r));
                            r
                    with | err ->
                        Methods.cost := prev;
                        raise err
            end
        in
        let fuse_iteration_status = Status.create "Fusing Trees" None "" in
        try
            run := { !run with trees = !trees };
            let len = Sexpr.length !trees in
            trees := `Empty;
            for i = 0 to (2 * len) do
                if 0. < remaining_time () then begin 
                    incr fuse_counter;
                    Status.message fuse_iteration_status ("Generation " ^ string_of_int !fuse_counter);
                    update_information (`Others (!run, !run));
                end;
                Status.full_report fuse_iteration_status;
                stop_if_necessary `Fuse;
                let fus = 
                    let swap_args = 
                        [APOY tbr; APOY timeout:[`Dynamic remaining_time]] -->
                            add_randomized --> add_visited --> add_constraint
                    in
                    CPOY fuse (iterations:1, swap { swap_args })
                in
                let r = exec !run fus in 
                update_information (`Others (!run, r));
                run := r;
                stop_if_necessary `Other;
            done;
            let r = folder !run (command_processor (CPOY select (unique))) in
            let r = do_exhaustive_round r in
            trees := r.trees;
            run := { r with trees = `Empty };
            Status.finished fuse_iteration_status;
        with | Exit ->
            Status.finished fuse_iteration_status;
            raise Exit
    done;
    !run
    with | Exit -> 
        let r, iterations_counter, fuse_counter, hits, _ = collect_results false in
        Status.finished st;
        Status.user_message Status.Information 
            ("The search evaluated " ^ string_of_int iterations_counter ^ 
            " independent repetitions with ratchet and fusing " ^
            "for " ^ string_of_int fuse_counter ^ " generations. " ^
            "The shortest tree was found " ^ string_of_int hits ^ 
            " times.");
        folder r [`Unique]


let rec process_application run item = 
    let run = reroot_at_outgroup run in
    match item with
    | `Algn_Newkk | `Algn_Normal as meth when meth <> !Methods.algn_mode ->
            Methods.algn_mode := meth;
            process_application run `ReDiagnoseTrees
    | `Algn_Newkk | `Algn_Normal
    | `Interactive -> run
    | `Optimization opt_mode ->
        let old_opt_mode = !Methods.opt_mode in
        let () = Methods.set_opt_mode opt_mode in
        let c  = Numerical.compare_opt_mode old_opt_mode opt_mode in
        if c = 0 || c = 1
            then run
            else process_application run `ReDiagnoseTrees

    | `Normal | `Normal_plus_Vitamines | `Exhaustive_Weak 
    | `Exhaustive_Strong | `Iterative _ as meth -> 
        if !Methods.cost <> meth then
            match !Methods.cost with
            | `Normal_plus_Vitamines | `Normal | `Exhaustive_Strong | `Exhaustive_Weak ->
                begin match meth with
                    | `Iterative _ ->
                        let () = Methods.cost := meth in
                        process_application run `ReDiagnoseTrees
                    | _ -> 
                        let () = Methods.cost := meth in
                        run
                end
            | `Iterative _ -> 
                let () = Methods.cost := meth in
                process_application run `ReDiagnoseTrees
        else begin
            run
        end
    | `Exit -> exit 0
    | `Version ->
            Status.user_message Status.Information Version.string;
            run
    | `ChangeWDir dir ->
            let dir = Str.global_replace (Str.regexp "\\\\ ") " " dir in
            Sys.chdir dir;
            run
    | `PrintWDir ->
            Status.user_message Status.Information 
            ("The current working directory is " ^ (StatusCommon.escape
            (Sys.getcwd ())));
            run
    | `Recover -> do_recovery run
    | `ClearRecovered -> Queue.clear run.queue.Sampler.queue; run
    | `ClearMemory items -> 
            let has_item item = List.exists (fun x -> x = item) items in
            if has_item `Matrices then Matrix.flush ();
            Gc.full_major (); 
            run
    | `TimerInterval x -> Tabus.timer_interval := x; run
    | `HistorySize len -> Status.resize_history len; run
    | `Logfile file -> StatusCommon.set_information_output file; run
    | `Redraw -> Status.redraw_screen (); run
    | `SetSeed v -> process_random_seed_set run v
    | `Alias (n,x) -> process_sets run n x
    | `ReDiagnose -> update_trees_to_data true true run
    | `ReDiagnoseTrees -> rediagnose_trees false run
    | `Help item -> HelpIndex.help item; run
    | `Wipe -> empty ()
    | `Echo (s, c) ->
            let s = StatusCommon.escape s in
            let c = match c with
                | `Information -> Status.Information
                | `Error -> Status.Error
                | `Output str -> begin match str with
                    | None -> Status.Output (None, false, [])
                    | Some str -> Status.Output (Some str, false, [])
                end
            in
            Status.user_message c (s ^ "@\n%!");
            run
     | `Graph (filename, collapse) ->
             let run = reroot_at_outgroup run in
             let trees = Sexpr.to_list run.trees in
             let trees = Array.of_list trees in
             if 0 = Array.length trees then run
             else begin
                let trees =
                    Array.map
                        (fun x ->
                            let cost = (Ptree.get_cost `Adjusted x) in
                            cost, (TS.build_forest_as_tree collapse x ""))
                        trees
                in
                match filename with 
                    | Some filename ->
                        GraphicsPs.display "" filename trees; run
                    | None ->
                        Status.user_message Status.Information 
                            ("@[Interactive@ graphics@ are@ not@ supported@ "
                            ^"in@ this@ compiled@ version@ of@ POY.@ Here@ is@ "
                            ^"the@ ascii@ art@ though:@]");
                        process_application run (`Ascii (None, collapse))
            end
     | `Ascii (filename, collapse) ->
            let trees =
                Sexpr.map
                    (fun x ->
                        let cost = int_of_float (Ptree.get_cost `Adjusted x) in
                        let str = string_of_int cost in
                        cost, TS.build_forest_as_tree collapse x str)
                    run.trees
            in
            Sexpr.leaf_iter
                (fun (cost, x) ->
                    let r = AsciiTree.to_string ~sep:2 ~bd:2 false x in
                    Status.user_message (Status.Output (filename,false, [])) 
                        ("@[@[<v>@[Tree@ with@ cost@ " ^ string_of_int cost ^ "@]@,@[");
                    Status.user_message (Status.Output (filename,false, [])) r;
                    Status.user_message (Status.Output (filename, false, [])) "@]@]@]%!";)
                trees;
            run
     | `Memory filename ->
             let memory_usage = report_memory () in
             Status.user_message (Status.Output (filename, false,[])) memory_usage;
             run
     | `KML (plugin, csv, output_file) ->
            let plugin = match plugin with
                | None -> "default"
                | Some x -> x
            in
            let csv = match csv with
                | `Local f -> f
                | `Remote f -> failwith "The csv must be available locally"
            in
            Kml.KFile.kml ~plugin "POY analysis" output_file run.data csv run.trees;
            run
     | `InspectFile str ->
            try let (desc, _, _) = PoyFile.read_file str in
                let desc = match desc with
                     | None -> "No@ description@ available."
                     | Some d -> d
                in
                Status.user_message Status.Information
                    ("@[<v 2>" ^ (StatusCommon.escape str) ^
                     " is a POY file: @,@[" ^ desc ^ "@]@]");
                run
             with | _ -> 
                Status.user_message Status.Information
                    ("The@ file@ " ^ StatusCommon.escape str ^
                     "@ is@ not@ a@ valid" ^ "@ POY@ fileformat.");
                run


let range_timer = ref (Timer.start ())


let get_character_costs trees = 
    (* We will first define a function to compute the median between a pair of
    * lists of non additive characters as generated by the NonaddCS* modules. *)
    let nonadd_median = 
        (* We don't hold the cost of individual non-additive characters within
        * the tree, therefore, we will have to recompute them completely from
        * scratch in to be able to produce the ci and ri scores. We first define
        * the extra cost for a pair of non-additive characters *)
        let distance a b = if 0 <> a land b then 0. else 1.
        and median a b = 
            if 0 <> a land b then a land b else a lor b in
        (* Now that we can compute something for a pair of characters, we define
        * a function that compares a pair of lists of the same length and
        * produces their cost and median. *)
        let compare (codea, (a, va), costa) (codeb, (_, vb), costb) =
            assert (codea = codeb);
            (codea, (a, median va vb), costa +. costb +. distance va vb)
        in
        (* We are ready to define a function to process a pair of lists and
        * generate a new one, notice that what we hold are lists of lists of
        * characters ...  *)
        List.map2 (fun a b -> List.map2 compare a b)
    in
    (* Now  we define a function that can be used to produce the leafs of the
    * tree that we will be processing to generate the character cost for the
    * non-additive characters. The goal is to have a pair of functions
    * compatible with  the Ptree.post_order_downpass_style function. *)
    let leaf_downpass tree _ l =
        let generate get tolist = 
            let lst = get None (Ptree.get_node_data l tree) in
            let lst = List.map fst lst in 
            List.map tolist lst
        in
        (generate Node.get_nonadd_8 NonaddCS8.to_list) @ 
        (generate Node.get_nonadd_16 NonaddCS16.to_list) @
        (generate Node.get_nonadd_32 NonaddCS32.to_list)
    in
    let interior_downpass _ _ = nonadd_median in
    let nonadditive_characters tree = 
        (* We require that every tree has only one component, so we verify it *)
        let handles = Ptree.get_handles tree in
        if 1 = All_sets.Integers.cardinal handles then 
            let handle = All_sets.Integers.choose handles in
            Ptree.post_order_downpass_style (leaf_downpass tree) 
            interior_downpass handle tree
        else failwith "Each tree must have only one handle for indexes"
    in
    let get_character_cost tree =
        (* We need to collect the cost for every vertex and the roots of the
        * tree, so we define a function that can collect over a given handle *)
        let get_cost_handle handle acc = 
            match (Ptree.get_component_root handle tree).Ptree.root_median with
            | Some ((`Edge (a, b)), root) ->
                    let res = 
                        Tree.post_order_node_with_edge_visit_simple 
                        (fun parent vertex acc ->
                            let characters = 
                                let node = Ptree.get_node_data vertex tree in
                                Node.character_costs (Some parent) node
                            in
                            characters :: acc) 
                        (Tree.Edge (a, b)) tree.Ptree.tree []
                    in
                    ((Node.character_costs None root) :: res) :: acc
            | _ -> acc
        in
        let res = 
            []
            --> All_sets.Integers.fold get_cost_handle tree.Ptree.tree.Tree.handles
            --> List.map (List.map (fun x -> List.partition (function (`Sank, _, _) -> true | _ ->
                false) (List.filter (function (`NonAdd, _, _) -> false | _ -> true) x)))
        in
        let merge_lists res = 
            match res with
            | [] -> []
            | h :: t -> 
                    List.fold_left (List.map2 (fun (a, b, c) (d, e, f) -> 
                    assert (a = d); assert (b = e); a, b, c +. f)) h t
        in
        let res = List.filter (function [] -> false | _ -> true) res in
        let sankoff = 
            let res, _ = List.split (List.map List.hd res) in
            merge_lists res
        and additive = 
            let _, res = List.split (List.flatten res) in
            merge_lists res
        and nonadditive = 
            let res = nonadditive_characters tree in
            List.map (fun (x, _, y) -> (`NonAdd, x, y)) (List.flatten res)
        in
        let res = sankoff @ additive @ nonadditive in
        List.fold_left (fun (total, res) (_, code, cost) ->
            total +. cost, 
            (All_sets.IntegerMap.add code cost res)) 
        (0., All_sets.IntegerMap.empty) res
    in
    let lst = Sexpr.to_list trees in 
    List.map get_character_cost lst

let hennig_extensions = ["ss"; "hennig"; "hen"]
let nexus_extensions = ["nexus"; "nex"]

let check_suffix filename lst = 
    match filename with
    | None -> false
    | Some filename ->
            List.exists (Filename.check_suffix filename) lst

IFDEF USEPARALLEL THEN
(* A function to create a complete mask *)
let complete_mask com_size =
    let rec complete_mask bit com_size =
        if bit > com_size then bit - 1
        else complete_mask (bit lsl 1) com_size
    in
    complete_mask 1 com_size

let world_size = Mpi.comm_size Mpi.comm_world 

let my_rank = Mpi.comm_rank Mpi.comm_world 

(* A function to flip a given bit *)
let mask_bit bit complete_mask =
    if (bit land my_rank) = 0 then my_rank lor bit
    else my_rank land (max_int lxor bit)

let compute_other_rank bit = world_size --> complete_mask --> mask_bit bit

END

    module SymbolTable = struct
        type 'a arguments_pv = 
            [ `Float of string 
            | `Int of string 
            | `String of string
            | `Labeled of (string * 'a arguments_pv)
            | `Lident of string
            | `CommandArg of 'a
            | `List of 'a arguments_pv list ]

        type pc_pv = Analyzer.pc_pv

        exception ExistingSymbol of string

        let table = Hashtbl.create 1667

        let add symbol f requirements = 
            if Hashtbl.mem table symbol then raise (ExistingSymbol symbol)
            else begin
                Hashtbl.add Analyzer.dependency_table symbol requirements;
                Hashtbl.add table symbol f;
            end

        let replace symbol f requirements =
            Hashtbl.replace table symbol f;
            Hashtbl.replace Analyzer.dependency_table symbol requirements

        let execute symbol run arguments =
            try (Hashtbl.find table symbol) run arguments with
            | Not_found -> failwith ("Undefined function " ^ symbol)

        let execute (run : r) (command : Analyzer.pc_pv) =
            match command with
            | `Command (symbol, arguments) -> execute symbol run arguments


        module Registration = struct
            exception IllegalArgument of 
                (string * Analyzer.pc_pv Analyzer.arguments_pv) 

            module type Registrar = sig
                type arguments
                val default_argument : arguments
                val process_argument : 
                    arguments ->
                        Analyzer.pc_pv Analyzer.arguments_pv -> 
                            arguments
                val execute : r -> arguments -> r
                val name : string
                val dependencies : arguments -> Analyzer.dependencies
            end

            module Generator (Implementation : Registrar) = struct
                let process_arguments args =
                    List.fold_left Implementation.process_argument
                    Implementation.default_argument args 
                let to_register_execution run arguments = 
                    Implementation.execute run (process_arguments arguments)

                let to_register_dependencies arguments = 
                    Implementation.dependencies (process_arguments arguments)
                let () = 
                    replace Implementation.name to_register_execution
                    to_register_dependencies
            end

            let build run arguments = run
            let change_wdir run arguments = run
            let clear_memory run arguments = run
            let clear_recovered run arguments = run
            let discard run arguments = run
            let echo run arguments = run
            let exit run arguments = run 
            let fuse run arguments = run
            let help run arguments = run
            let inspect run arguments = run
            let load run arguments = run
            let perturb run arguments = run
            let print_wdir run arguments = run
            let read run arguments = run
            let recover run arguments = run
            let rediagnose run arguments = run
            let redraw run arguments = run
            let rename run arguments = run
            let report run arguments = run
            let run run arguments = run
            let save run arguments = run
            let search run arguments = run
            let select run arguments = run
            let set run arguments = run
            let stdsearch run arguments = run
            let store run arguments = run
            let support run arguments = run
            let swap run arguments = run
            let transform run arguments = run
            let use run arguments = run
            let version run arguments = run
            let wipe run arguments = run
        end
    end

let rec folder (run : r) meth =
    check_ft_queue run;
    let res = match meth with
    (* The following methods are only used by the parallel execution *)
    | `Barrier -> (* Wait for synchronization with every other process *)
        IFDEF USEPARALLEL THEN
            print_msg "Entering barrier";
            let print_io_messages () =
                let gotit, rank, tag = 
                    Mpi.iprobe Mpi.any_source Methods.io Mpi.comm_world
                in
                if gotit then
                    let (t : Status.c), (msg : string) = 
                        Mpi.receive rank tag Mpi.comm_world in
                    Status.user_message t msg
            in
            let request_barrier rank =
                Mpi.send () rank Methods.barrier Mpi.comm_world in
            let send_barrier = request_barrier in
            let wait_for_barrier rank = 
                Mpi.receive rank Methods.barrier Mpi.comm_world in
            let wait_for_request = wait_for_barrier in
            let rec do_barrier first_time bit =
                if my_rank = 0 then print_io_messages ();
                let other_rank = compute_other_rank bit in
                if bit < world_size then
                    if other_rank < world_size then
                        if other_rank < my_rank then
                            let () = wait_for_request other_rank in
                            let () = send_barrier other_rank in
                            Mpi.barrier Mpi.comm_world
                        else
                            if my_rank <> 0 then
                                let () = request_barrier other_rank in
                                let () = wait_for_barrier other_rank in
                                do_barrier true (bit lsl 1)
                            else 
                                let () =
                                    if first_time then request_barrier other_rank
                                in
                                let gotit, _, _ = 
                                    Mpi.iprobe other_rank Methods.barrier
                                    Mpi.comm_world
                                in
                                if gotit then
                                    let () = wait_for_barrier other_rank in
                                    do_barrier true (bit lsl 1)
                                else do_barrier false bit
                    else do_barrier true (bit lsl 1)
                else Mpi.barrier Mpi.comm_world
            in
            do_barrier true 1;
            run
        ELSE
            run 
        END

    | `GatherTrees (joiner, continue) ->
        IFDEF USEPARALLEL THEN
            print_msg "Entering Gather Trees";
            (* We define a function to reduce in parallel the results *)
            let receive_trees run other_rank =
        IFDEF USE_LARGE_MESSAGES THEN
                print_msg "It is time for me to receive the trees!!!!";
                let tbl = Queue.create () in
                let got_size = ref false in
                let trees = ref [||] in
                let stored_trees = ref [||] in
                let get_arr a = 
                    if a = 0 then trees
                    else stored_trees 
                in
                let msg_cntr = ref 0 in
                let expected_messages = ref max_int in
                let process_msg msg = 
                    incr msg_cntr;
                    print_msg ("Received message number " ^ string_of_int
                    !msg_cntr);
                    match msg with
                    | (a, b, `Arrays) ->
                            got_size := true;
                            trees := Array.make a Ptree.empty;
                            stored_trees := Array.make b Ptree.empty;
                            expected_messages := ((a * 4) + (b * 4) + 1);
                            print_msg ("I expect " ^ string_of_int
                            !expected_messages);
                    | x when not !got_size ->
                            Queue.push x tbl;
                    | (a, b, `Topology x) ->
                            let arr = get_arr a in
                            !arr.(b) <- { !arr.(b) with Ptree.tree = x }
                    | (a, b, `Component_Root (root, origin_cost)) ->
                            let arr = get_arr a in
                            !arr.(b) <- { !arr.(b) with Ptree.component_root = root
                            ; origin_cost = origin_cost };
                    | (a, b, `EdgeData x) ->
                            let arr = get_arr a in
                            !arr.(b) <- { !arr.(b) with Ptree.edge_data = x }
                    | (a, b, `NodeData x) ->
                            let arr = get_arr a in
                            !arr.(b) <- { !arr.(b) with Ptree.node_data = x }
                    | (a, b, _) -> failwith "Invalid message"
                in
                while !expected_messages > !msg_cntr do
                    print_msg "I expect a new message";
                    process_msg (Mpi.receive other_rank Methods.do_job
                    Mpi.comm_world)
                done;
                print_msg "I received all the messages";
                while not (Queue.is_empty tbl) do
                    process_msg (Queue.pop tbl);
                done;
                let trees = Sexpr.of_list (Array.to_list !trees) 
                and stored_trees = Sexpr.of_list (Array.to_list !stored_trees) in
                { run with trees = Sexpr.union run.trees trees; stored_trees =
                    Sexpr.union run.stored_trees stored_trees }
        ELSE
                print_msg "I will use cheap messages";
                let dec = Mpi.receive other_rank Methods.do_job Mpi.comm_world
                in
                print_msg ("Received from " ^ string_of_int other_rank);
                decode_trees print_msg dec run
        END
            in
            let send_trees run other_rank =
        IFDEF USE_LARGE_MESSAGES THEN
                let map_root r =
                    let nr = 
                        match r.Ptree.root_median with
                        | None -> None
                        | Some (x, y) -> Some (x, Node.force y)
                    in
                    { r with Ptree.root_median = nr }
                in
                print_msg "I will use long messages";
                let send_msg msg = Mpi.send msg other_rank Methods.do_job
                Mpi.comm_world
                in
                let sexpr_to_arr x = Array.of_list (Sexpr.to_list x) in
                let trees = sexpr_to_arr run.trees 
                and stored_trees = sexpr_to_arr run.stored_trees in
                let send_tree a b tree = 
                    print_msg "Sending topology";
                    send_msg (a, b, `Topology tree.Ptree.tree);
                    print_msg "Sending component";
                    send_msg (a, b, `Component_Root 
                    ((All_sets.IntegerMap.map map_root tree.Ptree.component_root),
                    tree.Ptree.origin_cost));
                    print_msg "Sending edge data";
                    send_msg (a, b, `EdgeData (Tree.EdgeMap.map Edge.force
                    tree.Ptree.edge_data));
                    print_msg "Sending node data";
                    send_msg (a, b, `NodeData (All_sets.IntegerMap.map
                    Node.force tree.Ptree.node_data));
                in
                print_msg "Sending array size";
                send_msg (Array.length trees, Array.length stored_trees,
                `Arrays);
                Array.iteri (send_tree 0) trees;
                Array.iteri (send_tree 1) stored_trees;
                ()
        ELSE
                print_msg "I will use cheap messages";
                let enc = encode_trees print_msg run in
                print_msg ("Sending to " ^ string_of_int other_rank);
                Mpi.send enc other_rank Methods.do_job
                Mpi.comm_world
        END
            in
            let reduce_in_pairs bit run =
                let other_rank = compute_other_rank bit in
                print_msg ("Exchanging between " ^ string_of_int my_rank ^ " and " 
                ^ string_of_int other_rank);
                if other_rank < world_size then
                    if other_rank > my_rank then 
                        (* I am in charge of doing the reduction *)
                        let () = print_msg "I'm in charge of the reduction" in
                        let run = 
                            receive_trees run other_rank 
                            --> (fun x -> List.fold_left folder x joiner)
                        in
                        let () = print_msg "Received the trees" in
                        let () = send_trees run other_rank in
                        let () = print_msg "Sent the selected trees" in
                        run
                    else
                        let () = print_msg "I will wait for reduction" in
                        let () = send_trees run other_rank in
                        let () = print_msg "Sent the trees" in
                        let run = { run with trees = `Empty; stored_trees =
                            `Empty }
                        in
                        let () = print_msg "Witing for the selected trees" in
                        let r = receive_trees run other_rank in
                        let () = print_msg "Received the selected trees" in
                        r
                else List.fold_left folder run joiner
            in
            let rec reduce_them_all bit run =
                if bit < world_size then
                    reduce_them_all (bit lsl 1) (reduce_in_pairs bit run)
                else List.fold_left folder run joiner
            in
            let run = { run with stored_trees = `Empty } in
            let run = reduce_them_all 1 run in
            print_msg "Finished all gather";
            List.fold_left folder run continue
        ELSE 
                run 
        END

    | `GatherJackknife ->
        IFDEF USEPARALLEL THEN
            print_msg "Entering jackknife";
            let res = Mpi.allgather (encode_jackknife run) Mpi.comm_world in
            let run =
                Array.fold_left 
                (fun run set ->
                    let jckn = decode_jackknife set run in
                    add_jackknifes run jckn)
                run
                res
            in
            print_msg "Exiting jackknife";
            run
        ELSE
            run
        END
    | `GatherBremer ->
        IFDEF USEPARALLEL THEN
            print_msg "Entering bremer";
            let res = Mpi.allgather (encode_bremer run) Mpi.comm_world in
            let run = 
                Array.fold_left 
                (fun run set ->
                    let bmr = decode_bremer set run in
                    add_bremers run bmr)
                run
                res
            in
            print_msg "Exiting bootstrap";
            run
        ELSE 
            run
        END

    | `GatherBootstrap ->
        IFDEF USEPARALLEL THEN
            print_msg "Entering bootstrap";
            let res = Mpi.allgather (encode_bootstrap run) Mpi.comm_world in
            let run = 
                Array.fold_left 
                (fun run set ->
                    let bstp = decode_bootstrap set run in
                    add_boostraps run bstp)
                run
                res
            in
            print_msg "Exiting bootstrap";
            run
        ELSE
            run
        END
    | `SelectYourTrees ->
        IFDEF USEPARALLEL THEN
            print_msg "Entering select trees";
            let run = filter_my_trees run in
            print_msg "Exiting select trees";
            run
        ELSE
            run
        END

    | `Skip
    | `Entry -> run
    | `Plugin (name, args) ->
            if Hashtbl.mem plugins name then
                let f = Hashtbl.find plugins name in
                f args run
            else begin
                Status.user_message Status.Error 
                ("There is no command " ^ name);
                failwith ("Illegal command " ^ name)
            end
    | `StoreTrees -> 
            { run with trees = `Empty;
                       stored_trees = run.trees }
    | `UnionStored ->
            { run with trees = Sexpr.union run.stored_trees run.trees;
                       stored_trees = `Empty }
    | `GetStored ->
            { run with trees = run.stored_trees;
                       stored_trees = `Empty }
    | `OnEachTree (dosomething, mergingscript) ->
            let name = emit_identifier () in
            let run = folder run (`Store ([`Data], name)) in
            let run = 
                if not (List.exists (function #Methods.transform -> true | _ -> false) dosomething) then
                    { run with original_trees = run.trees } 
                else run
            in
            let run = 
                let eta = true in
                Sexpr.fold_status
                    "Running pipeline on each tree" ~eta
                    (on_each_tree folder (`Set ([`Data], name)) dosomething mergingscript)
                    run run.trees
            in
            let run = 
                {run with   trees = run.stored_trees;
                            original_trees = `Empty; 
                            stored_trees = `Empty; } 
            in
            let run = folder run (`Discard ([`Data], name)) in
            run
    | `ParallelPipeline (times, todo, composer, continue) ->
            let name = emit_identifier () in
            let run = folder run (`Store ([`Data], name)) in
            let st = Status.create "Running Pipeline" (Some times) "times" in
            let run = ref run in
            let for_each = todo @ composer in
            let timer = Timer.start () in
            for adv = 1 to times do
                run := folder !run (`Set ([`Data], name));
                run := List.fold_left folder !run for_each;
                let msg = Timer.status_msg (Timer.wall timer) adv times in
                Status.full_report ~adv ~msg st;
            done;
            run := folder !run (`Discard ([`Data], name));
            Status.finished st;
            List.fold_left folder !run continue
    (* The following methods are user friendly *)
    | #Methods.tree_handling as meth ->
            warn_if_no_trees_in_memory run.trees;
            process_tree_handling run meth
    | #Methods.characters_handling as meth ->
            update_trees_to_data false true (process_characters_handling run meth)
    | #Methods.taxa_handling as meth ->
            process_taxon_filter run meth 
    | #Methods.application as meth ->
            process_application run meth
    | #Methods.input as meth ->
            process_input run meth 
    | #Methods.transform as meth ->
            process_transform run meth
    | #Methods.build as meth ->
        begin match run.nodes with
        | [] ->
            Status.user_message Status.Warning "There@ is@ no@ data@ in@ memory!";
            run
        | _  ->
            let build_initial = Build.build_initial_trees in
            begin match MainBuild.get_transformations meth with
                | [] ->
                    let trees = build_initial run.trees run.data run.nodes meth in
                    { run with trees = trees; }
                | trans ->
                    let run,untransforms = temporary_transforms trans run in
                    let tree = build_initial run.trees run.data run.nodes meth in
                    let run = { run with trees = tree } in
                    List.fold_left (fun r f -> f r) run untransforms
            end
        end
    | #Methods.local_optimum as meth ->
            warn_if_no_trees_in_memory run.trees;
            let sets =
                let `LocalOptimum m = meth in
                TreeSearch.sets m.Methods.tabu_join run.data run.trees 
            in
            let do_search run = match is_forest meth with
                | Some cost ->
                    PTS.forest_search run.data run.queue cost meth run.trees
                | None ->
                    PTS.find_local_optimum run.data run.queue run.trees sets meth
            in
            begin match TreeSearch.get_transformations meth with
                | [] ->
                    let trees = do_search run in
                    { run with trees = trees }
                | trans when only_multistatic trans ->
                    let (run, untransforms) = temporary_transforms trans run in
                    let trees = do_search run in
                    let run = { run with trees = trees } in
                    let run = List.fold_left (fun r f -> f r) run untransforms in
                    { run with trees = run.trees }
                | trans ->
                    let runs = explode_trees run in
                    let runs = 
                        Sexpr.map_status
                            "Transforming each tree independently"
                            ~eta:true (temporary_transforms trans) 
                            runs
                    in
                    let run_and_untransform (run, untransforms) =
                        let trees = do_search run in
                        let run = { run with trees = trees } in
                        let run = 
                            List.fold_left (fun r f -> f r) run untransforms
                        in
                        Sexpr.first run.trees
                    in
                    let trees = Sexpr.map run_and_untransform runs in
                    { run with trees = trees }
            end
    | `StandardSearch (max_time, min_time, min_hits, max_memory, target_cost, visited, user_constraint) ->
            let get_default x def = match x with
                | None -> def
                | Some x -> x
            in
            let max_time = get_default max_time 3600. in
            let min_time = get_default min_time max_time in
            let max_mem  = get_default max_memory (2 * 1000 * 1000 * (1000 / (Sys.word_size / 8))) in
            automated_search (List.fold_left folder) max_time min_time
                        max_mem (get_default min_hits max_int) target_cost
                        visited user_constraint run

    | #Methods.perturb_method as meth ->
            warn_if_no_trees_in_memory run.trees;
            { run with trees = CT.perturbe run.data run.trees meth }
    | `Fusing ((_, _, _, _, x, _) as params) ->
            warn_if_no_trees_in_memory run.trees;
            begin
                try { run with trees = PTS.fusing run.data run.queue run.trees params }
                with | Failure "Tree fusing: must have at least two trees" -> run
            end
    | `Bootstrap (it, a, b, c) -> 
            let meth = `Bootstrap (it, a, b, (run.data.Data.root_at)) in
            let run = reroot_at_outgroup run in
            { run with bootstrap_support = 
                Some (it, S.support run.trees run.nodes meth run.data run.queue) }
    | `Jackknife (a, it, b, c, _) ->
            let meth = `Jackknife (a, it, b, c, (run.data.Data.root_at)) in
            let run = reroot_at_outgroup run in
            { run with jackknife_support = 
                Some (it, S.support run.trees run.nodes meth run.data run.queue) }
    | `Bremer (local_optimum, build, my_rank, modul) ->
            assert (my_rank < modul);
            warn_if_no_trees_in_memory run.trees;
            let run = reroot_at_outgroup run in
            { run with bremer_support = 
                S.bremer_support run.trees my_rank modul run.nodes run.trees local_optimum build 
                run.data run.queue }
    | #Methods.escape_local as meth ->
            warn_if_no_trees_in_memory run.trees;
            let (`PerturbateNSearch (tr, _, search_meth, _, _)) = meth in
            let choose_best trees = 
                (* A function to select the best between the
                resulting trees and the previously existing trees *)
                let initialtrees = Sexpr.to_list run.trees 
                and othertrees = Sexpr.to_list trees in
                let module TreeSet = Set.Make (Ptree.Fingerprint) in
                let no_need, othertrees = 
                    let initialtrees = 
                        let add_elements acc x =
                            TreeSet.add (Ptree.Fingerprint.fingerprint x) acc
                        in
                        let a = 
                            Sexpr.fold_left add_elements TreeSet.empty run.trees
                        in
                        let a = 
                            Sexpr.fold_left add_elements a run.stored_trees in
                        Sexpr.fold_left add_elements a run.original_trees
                    in
                    List.partition (fun x ->
                        let x = Ptree.Fingerprint.fingerprint x in
                        TreeSet.mem x initialtrees) othertrees
                in
                let othertrees = 
                    let _, acc =
                        List.fold_left (fun ((set, acc) as pair) x ->
                            let fp = Ptree.Fingerprint.fingerprint x in
                            if TreeSet.mem fp set then pair
                            else TreeSet.add fp set, x :: acc) 
                        (TreeSet.empty, []) othertrees
                    in
                    Sexpr.of_list acc
                in
                if 0 <> Sexpr.length othertrees then
                    let newres = 
                        folder { run with trees = othertrees } 
                        (search_meth :> script)
                    in
                    let ntrees = List.length initialtrees in
                    let all_trees = 
                        Sexpr.of_list (initialtrees @ no_need @ (Sexpr.to_list newres.trees)) 
                    in
                    folder {run with trees = all_trees } (`BestN (Some ntrees))
                else run
            in
            (match tr with
            | [] ->
                let trees = 
                    CT.escape_local run.data run.queue run.trees meth 
                in
                choose_best trees
            | trans when only_multistatic trans ->
                let (run, untransforms) = temporary_transforms trans run in
                let trees = CT.escape_local run.data run.queue run.trees
                meth in
                let run = { run with trees = trees } in
                let run = List.fold_left (fun r f -> f r) run untransforms in
                choose_best run.trees
            | trans ->
                let run_and_untransform (run, untransforms) =
                    let trees = 
                        CT.escape_local run.data run.queue run.trees meth 
                    in
                    let run = { run with trees = trees } in
                    let run = 
                        List.fold_left (fun r f -> f r) run untransforms
                    in
                    let run = folder run (`BestN (Some 1)) in
                    Sexpr.first run.trees
                in
                let runs = explode_trees run in
                let trees =
                    Sexpr.map_status 
                        "Transforming and searching each tree independently"
                        ~eta:true
                        (fun x ->
                            run_and_untransform (temporary_transforms trans x))
                        runs
                in
                choose_best trees)
    | #Methods.runtime_store as meth -> 
            runtime_store (update_trees_to_data false false) run meth 
    | `Repeat (n, comm) ->
            let res = ref run in
            for i = 1 to n do
            res := (List.fold_left folder !res comm);
            done;
            !res
    | `ReadScript files ->
            let file_folder run item =
                try folder run item
                with | err ->
                    let msg = StatusCommon.escape (Printexc.to_string err) in
                    Status.user_message Status.Error msg;
                    raise (Error_in_Script (err, run))
            in
            let script = PoyCommand.read_script_files true 
                (List.map (fun x -> `Filename x) files) in
            let script = Sexpr.of_list script in
            Sexpr.fold_status "Running commands" ~eta:true file_folder run script
    | #Methods.report as meth ->
            (* Update the trees to reflect the rooting we want *)
            let run = reroot_at_outgroup run in 
            match meth with
            | `SequenceStats (filename, ch) ->
                let arr =
                    let all_of_them = Data.sequence_statistics ch run.data in
                    let arr =
                        Array.of_list 
                            (List.map
                                (fun (name, stats) ->
                                    let cnt = float_of_int stats.Data.sequences in
                                    let len_avg = (float_of_int stats.Data.sum_lengths) /. cnt in
                                    let dis_avg =
                                        let numr = stats.Data.sum_distances
                                        and denm = ((cnt *. cnt) /. 2.) -. (cnt /. 2.) in
                                        numr /. denm
                                    in
                                    [|  name;
                                        string_of_int stats.Data.max_length; 
                                        string_of_int stats.Data.min_length;
                                        string_of_float len_avg; 
                                        string_of_float stats.Data.max_distance; 
                                        string_of_float stats.Data.min_distance; 
                                        string_of_float dis_avg |])
                                all_of_them)
                    in
                    Array.init 
                        (1 + Array.length arr)
                        (function
                            | 0 ->
                                [|"Character"; "Max Length"; "Min Length"; "Average Length";
                                  "Maximum Distance"; "Minimum Distance"; "Average Distance"|]
                            | n -> arr.(n - 1))
                and fo = Status.Output (filename, false, []) in
                Status.user_message fo 
                "@{<b>Sequence Statistics:@}@[<v 2>@,";
                Status.output_table fo arr;
                Status.user_message fo "@]\n%!";
                run
            | `Ci (filename, ch) | `Ri (filename, ch) as meth ->
                    let realch = match ch with
                        | Some x -> x
                        | None -> `All
                    in
                    let print_table title table = 
                        let fo = Status.Output (filename, false, []) in
                        Status.user_message fo 
                            ("@{<b>" ^ title ^ " Statistics:@}@[<v 2>@,");
                        Status.output_table fo table;
                        Status.user_message fo "@]\n%!";
                    in
                    let min_list =
                        Data.apply_on_static 
                            AddCS.min_possible_cost
                            NonaddCS8.min_possible_cost 
                            SankCS.min_possible_cost
                            (fun _ _ -> 0.0)
                            (fun _ -> 0.0)
                            realch
                            run.data
                    in
                    let add_list lst = 
                        List.fold_left (fun acc (_, cost) -> cost +. acc) 0. lst
                    in
                    let get_tree_cost_and_apply f =
                        let trees = Sexpr.to_list run.trees in
                        List.map 
                        (fun x -> 
                            let x = Ptree.get_cost `Adjusted x in
                            [|string_of_float x; f x|]) trees 
                    in
                    let table, title = match meth with
                        | `Ci (filename, ch) ->
                            let ci ci_val x = string_of_float (100. *. ci_val /. x) in
                            begin match ch with
                                | None ->
                                    (* We are dealing with all the characters, so
                                    we print the tree statistics *)
                                    let ci_val = add_list min_list in
                                    let trees = get_tree_cost_and_apply (ci ci_val) in
                                    ([|"Tree Cost"; "CI"|] :: trees), "CI"
                                | Some ch -> 
                                    (* Otherwise we print each individual character score *)
                                    let trees = get_character_costs run.trees in
                                    let trees = 
                                        List.map
                                            (fun (tree_cost, tree_chars_costs) ->
                                                List.map
                                                    (fun (code, ci_val) ->
                                                        let name = Data.code_character code run.data in
                                                        try let tree_char_cost = All_sets.IntegerMap.find code tree_chars_costs in
                                                            let res = ci ci_val tree_char_cost in
                                                            [|name; res|]
                                                        with | Not_found ->
                                                            if ci_val = 0.  then
                                                                [|name; "uninformative"|]
                                                            else
                                                                [|name; "100."|])
                                                    min_list, tree_cost)
                                            trees
                                    in
                                    List.flatten
                                        (List.map
                                            (fun (lst, tree_cost) ->
                                                let n = [|"Character of tree with cost " ^ string_of_float tree_cost; "ci"|] in
                                                n :: lst)
                                            trees), "ci"
                            end
                        | `Ri (filename, ch) ->
                            begin
                                let max_list = 
                                    Data.apply_on_static 
                                        (fun _ -> 0.) NonaddCS8.max_possible_cost
                                        SankCS.max_possible_cost (fun _ _ -> 0.)
                                        (fun _ -> 0.0)
                                        realch run.data
                                in
                                let ri max_cost min_cost x =
                                    if max_cost = min_cost then 
                                        "uninformative" 
                                    else
                                        string_of_float (100. *. (max_cost -. x) /. (max_cost -. min_cost))
                                in
                                match ch with
                                | None ->
                                    let res1 = add_list min_list
                                    and res2 = add_list max_list in
                                    let trees = 
                                        get_tree_cost_and_apply (ri res2 res1)
                                    in
                                    ([|"Tree Cost"; "RI"|] :: trees), "RI"
                                | Some ch ->
                                    let trees = let trees = get_character_costs run.trees in
                                    List.map 
                                        (fun (tree_cost, tree_chars_costs) ->
                                            List.map2
                                                (fun (code, min_val) (code, max_val) ->
                                                    let name =
                                                        Data.code_character code run.data
                                                    in
                                                    try let tree_char_cost =
                                                            All_sets.IntegerMap.find code tree_chars_costs
                                                        in
                                                        let res = ri max_val min_val tree_char_cost in
                                                        [|name; res|]
                                                    with | Not_found ->
                                                        if min_val = max_val then
                                                            [|name; "uninformative"|]
                                                        else [|name; "100."|])
                                                min_list max_list, tree_cost)
                                        trees
                                    in
                                    List.flatten
                                        (List.map
                                            (fun (lst, tree_cost) ->
                                                let n = [|"Character of tree with cost " ^ string_of_float tree_cost; "ri"|] in
                                                n :: lst)
                                            trees),"ri"
                            end
                        in
                        let table = 
                            Array.of_list 
                            (List.filter (fun arr -> arr.(1) <> "uninformative")
                            table)
                        in
                        print_table title table;
                        run
            | `CompareSequences (filename, complement, ch1, ch2) ->
                let all_of_them = Data.compare_pairs ch1 ch2 complement run.data in
                List.iter
                    (fun (n1, n2, c) ->
                        Status.user_message (Status.Output (filename, false, []))
                            (Printf.sprintf "@[%s %s %f@]@\n"
                                (StatusCommon.escape n1) (StatusCommon.escape n2) c))
                    all_of_them;
                run
            | `ExplainScript (script, filename) ->
                let script = PoyCommand.of_file false script in
                Analyzer.explain_tree filename script;
                run
            | `Pairwise (filename,chars) ->
                failwith "NOT DONE"
            | `LKSites (filename,chars) ->
                let ft = Status.output_table (Status.Output (filename, false, [])) in
                let items = Sexpr.length run.trees in
                let chars = Array.of_list (Data.get_chars_codes_comp run.data chars) in
                if items = 0 then begin
                    Status.user_message Status.Warning
                        ("There@ are@ no@ active@ trees@ in@ memory!");
                    run
                end else begin
                    Sexpr.leaf_iter (fun x -> report_lk_sites ft x chars)
                                    run.trees;
                    run
                end
            | `Model (filename,chars) ->
                let fo = Status.user_message (Status.Output (filename, false, [])) in
                let ft = Some (Status.output_table (Status.Output (filename, false, []))) in
                begin match (Sexpr.to_list run.trees) with
                    | [] -> 
                        List.iter
                            (fun xs ->
                                let chars =
                                    Data.get_code_from_characters_restricted
                                                `Likelihood run.data (`Some xs)
                                in
                                match chars with
                                | []    -> ()
                                | chars ->
                                    let model  = Data.get_likelihood_model run.data chars
                                    and name   = Data.get_character_set_name run.data chars
                                    and ntaxa  = run.data.Data.number_of_taxa in
                                    let cname = match name with | Some name -> name | None -> "" in
                                    if cname <> "" then
                                        fo ("@[<hov 0>Set Name: "^cname^"@]@\n");
                                    fo ("@[<hov 0>Number of taxa: "^string_of_int ntaxa^"@]@\n");
                                    MlModel.output_model fo ft `Hennig model None)
                            (Data.categorize_likelihood_chars_by_model_comp run.data chars)
                    | trees ->
                        List.iter
                            (fun t ->
                                let tname = match t.Ptree.tree.Tree.tree_name with
                                    | Some tname -> tname
                                    | None -> ""
                                and cats =
                                    Data.categorize_likelihood_chars_by_model_comp
                                                             t.Ptree.data chars
                                in
                                List.iter
                                    (fun xs ->
                                        let model = Data.get_likelihood_model t.Ptree.data xs
                                        and cost  = TreeOps.total_cost t `Adjusted (Some xs)
                                        and length= TreeOps.tree_size t (Some xs)
                                        (* and prior = TreeOps.prior_cost t (Some xs) *)
                                        and cname = Data.get_character_set_name t.Ptree.data xs in
                                        let cname = match cname with | Some cname -> cname | None -> ""
                                        and ntaxa = t.Ptree.data.Data.number_of_taxa in
                                        if tname <> "" then
                                            fo ("@[<hov 0>Tree Name: "^tname^"@]@\n");
                                        if cname <> "" then
                                            fo ("@[<hov 0>Set Name: "^cname^"@]@\n");
                                        fo ("@[<hov 0>Number of taxa: "^string_of_int ntaxa^"@]@\n");
                                        fo ("@[<hov 0>Tree Size: "^string_of_float length^"@]@\n");
                                        fo ("@[<hov 0>Log-Likelihood: "^string_of_float (~-.cost)^"@]@\n");
                                        (* fo ("@[<hov 0>Log-Prior: "^string_of_float (~-.prior)^"@]@\n");*)
                                        MlModel.output_model fo ft `Hennig model None;
                                        fo "@\n@\n")
                                    cats)
                            (trees)
                end;
                run
            | `Script (filename,script) -> 
                run
            | `Nexus filename ->
                let ic = [ `Branches; `NexusStyle; `Collapse false; ] in
                let parsed_trees,labeling = PTS.process_trees ic run.trees in
                CompOut.to_nexus run.data parsed_trees labeling ic filename;
                run
            | `FasWinClad filename ->
                CompOut.to_faswincladfile run.data filename;
                run
            | `Consensus (filename, v) ->
                PTS.output_consensus run.data run.trees filename v false;
                run
            | `GraphicConsensus (filename, v) ->
                PTS.output_consensus run.data run.trees filename v true;
                run
            | `History size ->
                Status.resize_history size;
                run
            | `Dataset filename ->
                if check_suffix filename nexus_extensions then
                    folder run (`Nexus filename)
                else if check_suffix filename hennig_extensions then
                    folder run (`FasWinClad filename)
                else begin
                   let fmt = (Data.to_formatter [] run.data) in
                    PoyFormaters.data_to_status filename fmt;
                    (* Flush the formatter *)
                    Status.user_message (Status.Output (filename, false, [])) "%!"; 
                    run
                end
            | `Xslt (file, style) ->
                let () =
IFDEF USE_XSLT THEN
                    let filename, chout = Filename.open_temp_file "results" ".xml" in
                    close_out chout;
                    Status.user_message Status.Information 
                            ("Generating xml file in " ^ StatusCommon.escape filename);
                    let ofilename = Some filename in
                    let fmt = Data.to_formatter [] run.data in
                    let trs =
                        Sexpr.map
                            (fun tr ->
                                let run =
                                    update_trees_to_data ~classify:false
                                            false true {run with trees = `Single tr}
                                in
                                match run.trees with
                                | `Single tr -> TreeOps.to_formatter `Normal [] tr
                                | _          -> assert false)
                            run.trees
                    in
                    StatusCommon.Files.set_margin ofilename 0;
                    Status.user_message (Status.Output (ofilename, false, [])) " <Diagnosis>@\n";
                    PoyFormaters.data_to_status ofilename fmt;
                    Sexpr.leaf_iter (PoyFormaters.trees_to_formater ofilename []) trs;
                    Status.user_message (Status.Output (ofilename, false, [])) " </Diagnosis>@\n%!";
                    Xslt.process filename style file
ELSE
                    Status.user_message Status.Error 
                    "This version of POY was not compiled with XSLT support."
END
                in
                run
            | `Diagnosis (diag_report_type,filename) ->
                let trees = Sexpr.map (TreeOps.to_formatter diag_report_type []) run.trees  in
                Status.user_message (Status.Output (filename, false, [])) "@[";
                Sexpr.leaf_iter (PoyFormaters.trees_to_formater filename []) trees;
                Status.user_message (Status.Output (filename, false, [])) "@]%!";
                run
            | `GraphicDiagnosis (diag_report_type,filename) ->
                let trees = Sexpr.map (TreeOps.to_formatter diag_report_type []) run.trees  in 
                GraphicsPs.display_diagnosis "Diagnosis" filename trees;
                run
            | `TimeDelta (title, filename) ->
                let prev_time = !range_timer in
                range_timer := Timer.start ();
                let total_time = Timer.get_user prev_time in
                let title = StatusCommon.escape title in
                Status.user_message
                    (Status.Output (filename, false, []))
                    (Printf.sprintf "@[%s : %.3f s@]@,%!" title total_time);
                run
            | `MstR filename ->
                Build.report_mst run.data run.nodes filename;
                run
            | `TreeCosts filename ->
                let res = Buffer.create 97 in
                Sexpr.leaf_iter (fun x ->
                    let cost = Ptree.get_cost `Adjusted x in
                    let str = string_of_float cost in
                    Buffer.add_string res str;
                    Buffer.add_string res ":") 
                run.trees;
                Buffer.add_string res "\n%!";
                Status.user_message (Status.Output (filename, false, [])) (Buffer.contents res);
                run
            | `TreesStats filename ->
                let fo = Status.Output (filename, false, []) in
                let htbl = Hashtbl.create 19 in
                let adder tree = 
                    let cost = Ptree.get_cost `Adjusted tree in
                    if Hashtbl.mem htbl cost then begin
                        let c = Hashtbl.find htbl cost in
                        Hashtbl.remove htbl cost;
                        Hashtbl.add htbl cost (c + 1)
                    end else Hashtbl.add htbl cost 1
                in
                Sexpr.leaf_iter adder run.trees;
                let arr = 
                    Array.make_matrix (1 + Hashtbl.length htbl) 2 (`Int 0) in
                let folder a b cnt =
                    arr.(cnt).(0) <- `Float a;
                    arr.(cnt).(1) <- `Int b;
                    cnt + 1
                in
                let comparison a b = match a.(0), b.(0) with
                    | `Int a, `Int b -> a - b
                    | `Float a, `Float b -> compare a b
                    | `Int _, _ -> -1
                    | _, _ -> 1
                in
                ignore( Hashtbl.fold folder htbl 1 );
                Array.sort comparison arr;
                let arr =
                    Array.map
                        (Array.map (function
                                    | `Float x -> string_of_float x
                                    | `Int x -> string_of_int x))
                        arr
                in
                arr.(0).(0) <- "Tree length   ";
                arr.(0).(1) <- "Number of hits";
                Status.user_message fo "@{<b>Trees Found:@}@[<v 2>@,";
                Status.output_table fo arr;
                Status.user_message fo "@]\n%!";
                run
            | `SearchStats filename ->
                let fo = Status.Output (filename, false, []) in
                let costs =
                    All_sets.FloatMap.fold
                        (fun a b acc -> [|string_of_float a; string_of_int b|] :: acc)
                        run.search_results.tree_costs_found
                        []
                in
                let costs =
                    List.sort
                        (fun a b -> match a, b with
                            | [|a; _|], [|b; _|] ->
                                compare (float_of_string a) (float_of_string b)
                            | _ -> assert false)
                        costs
                in
                let costs = 
                    [|"# of Builds + TBR "; string_of_int run.search_results.total_builds|] ::
                    [|"# of Fuse         "; string_of_int run.search_results.total_fuse|] ::
                    [|"# of Ratchets     "; string_of_int run.search_results.total_ratchet|] ::
                    [|"                  "; "              "|] ::
                    [|"Tree length       "; "Number of hits"|] :: costs 
                in
                Status.user_message fo "@{<b>Search Stats:@}@[<v 2>@,";
                Status.output_table fo (Array.of_list costs);
                Status.user_message fo "@]\n%!";
                run
            | `Trees (ic, filename) ->
                let rec remove_style acc = function
                    | [] -> acc
                    | `NexusStyle :: tl 
                    | `HennigStyle :: tl-> remove_style acc tl
                    | hd :: tl -> remove_style (hd::acc) tl
                in
                let ic = 
                    if check_suffix filename nexus_extensions then
                        `NexusStyle :: (remove_style [] ic)
                    else if check_suffix filename hennig_extensions then
                        `HennigStyle :: (remove_style [] ic)
                    else ic
                in
                PTS.report_trees ic filename run.trees;
                run
            | `CrossReferences (chars, filename) ->
                Data.report_taxon_file_cross_reference chars run.data filename;
                run
            | `KolmoMachine filename ->
                let data =
                    Data.report_kolmogorov_machine filename run.data
                in
                { run with data = data }
            | `TerminalsFiles filename ->
                Data.report_terminals_files filename
                        run.data.Data.taxon_files run.data.Data.ignore_taxa_set;
                run
            | `Nodes filename ->
                let fo = Status.Output (filename, false, []) in
                let nodes = List.map Node.to_string run.nodes in
                List.iter (Status.user_message fo) nodes;
                Status.user_message fo "%!";
                run
            | `Supports _ 
            | `GraphicSupports _ as meth ->
                    handle_support_output run meth;
                    run
            | `Clades fn ->
                let counter = ref (-1) in
                Sexpr.leaf_iter
                    (output_clade_file run.data fn
                        (fun () -> incr counter; !counter))
                    run.trees;
                Status.user_message (Status.Output (None, false,[])) "%!";
                run
            | #Methods.diagnosis as meth ->
                warn_if_no_trees_in_memory run.trees;
                let () = 
                    Sexpr.leaf_iter (fun x -> D.diagnosis run.data x meth) run.trees 
                in
                run
            | `Save (fn, comment) ->
                begin try
                    let srun = {run with runtime_store = create_hstb ();
                                         data_store = create_hstb (); }
                    in
                    PoyFile.store_file (comment, srun, !Methods.cost) fn;
                    run
                with | err ->
                    let msg = 
                        "Could@ not@ open@ the@ file@ " ^ fn ^ 
                        "@ due@ to@ a@ system@ error.@ The@ error@ message@ is@ : " ^ 
                        StatusCommon.escape (Printexc.to_string err)
                    in
                    Status.user_message Status.Error msg;
                    run
                end
            | `Load fn ->
                begin try
                    let (comment, run, cost_mode)  = PoyFile.read_file fn in
                    Methods.cost := cost_mode;
                    let descr = match comment with
                        | None   -> "No description available."
                        | Some d -> d
                    in
                    Status.user_message Status.Information 
                        ("@[<v 2>Loading file " ^ StatusCommon.escape fn 
                        ^ "@,@[" ^ StatusCommon.escape descr ^ "@]@]");
                    run
                with | PoyFile.InvalidFile ->
                    let msg =
                        "The@ file@ " ^ StatusCommon.escape fn ^ 
                        " is@ not a@ valid" ^ "@ POY@ fileformat."
                    in
                    Status.user_message Status.Error msg;
                    run
                end
            | `Root where ->
                let set_root_at run code_opt =
                    let tres = 
                        Sexpr.map
                            (fun t -> 
                                {t with
                                    Ptree.data = {t.Ptree.data with 
                                                    Data.root_at = where;}})
                            run.trees
                    in
                    { run with  data = { run.data with Data.root_at = where };
                                trees = tres; }
                in
                (* we test that the where root is a taxon; why not an internal
                   node? The tree search can be derived from a search, thus
                   internal nodes are not consistent in each tree; except for
                   leaves. *)
                begin match where with
                    | None   -> set_root_at run where
                    | Some x ->
                        if All_sets.IntegerMap.mem x run.data.Data.taxon_codes then
                            set_root_at run where
                        else
                            let msg =
                                "Terminal@ "^StatusCommon.escape (string_of_int x)^
                                "@ not@ found.@ To@ set@ the@ root@ I@ must@ have"^
                                "@ loaded@ some@ data@ for@ it.@ The@ assigned"^
                                "@ root@ will@ remain@ unchanged."
                            in
                            let () = Status.user_message Status.Error msg in
                            run
                end
            | `RootName name ->
                if All_sets.StringMap.mem name run.data.Data.taxon_names then
                    folder run (`Root (Some (Data.taxon_code name run.data)))
                else
                    let msg =
                        "Terminal@ " ^ StatusCommon.escape name ^ "@ not@ "^
                        "found.@ To@ set@ the@ root@ I@ must@ have@ loaded"^
                        "@ some@ data@ for@ it.@ The@ assigned@ root@ will"^
                        "@ remain@ unchanged."
                    in
                    let () = Status.user_message Status.Error msg in
                    run
    in
    StatusCommon.Files.flush (); (* flush all open files *)
    res
              

let deal_with_error output_file run tmp err =
    Status.user_message  Status.Error "Dumping poy computation";
    flush stdout;
    let ch = open_out_bin output_file in
    Marshal.to_channel ch (run, tmp) [];
    close_out ch;
    let msg = StatusCommon.escape (Printexc.to_string err) in
    Status.user_message Status.Error msg;
    raise err

let run ?(folder=folder) ?(output_file=Arguments.default_dump_file) ?(start=(empty ())) lst =
    let rec continue run tmp =
        let my_folder printit run h =
            if ndebug_no_catch then begin
                let run = folder run h in
                if printit then
                    Status.user_message
                        (Status.SearchReport)
                        (SearchInformation.show_information
                            (Some run.trees) (Some run.data) None None);
                run
            end else try 
                let run = folder run h in
                Status.user_message
                    (Status.SearchReport)
                    (SearchInformation.show_information
                        (Some run.trees) (Some run.data) None None);
                run
                with
                    | Error_in_Script (err, run) ->
                        deal_with_error output_file run tmp err
                    | err ->
                        deal_with_error output_file run tmp err
        in
        match tmp with
        | h::t -> continue (my_folder (t = []) run h) t
        | []   -> run
    in
    continue start lst 

let get_dump ?(file=Arguments.default_dump_file) () =
    let ch = open_in_bin file in
    let (r, lst) = Marshal.from_channel ch in
    ((r : r), (lst : script list))

let restart ?(file=Arguments.default_dump_file) () =
    let ch = open_in_bin file in
    let (r, lst) = Marshal.from_channel ch in
    run ~start:r lst

let console_run_val = ref (empty ())

let parsed_run lst = 
    let res = run ~start:!console_run_val lst in
    console_run_val := res

let console_run str = 
    let todo = PoyCommand.of_string true str in
    let res = run ~start:!console_run_val todo in
    console_run_val := res

let channel_run ch = 
    let todo = PoyCommand.of_channel true ch in
    let res = run ~start:!console_run_val todo in
    console_run_val := res

let get_console_run () = !console_run_val

let set_console_run r = console_run_val := r


    module PhyloTree  = struct
        module PtreeSearch = Ptree.Search (Node) (Edge) (TreeOps)
        type phylogeny = (a, b) Ptree.p_tree
        let get_cost x = Ptree.get_cost `Adjusted x
        let fold_edges f acc tree = List.fold_left f acc (Ptree.get_edges_tree tree)
        let fold_nodes f acc tree = List.fold_left f acc (Ptree.get_nodes tree)
        let fold_vertices f acc tree = List.fold_left f acc (Ptree.get_node_ids tree)
        let add_node_data = Ptree.add_node_data
        let get_node_data = Ptree.get_node_data
        let add_edge_data = Ptree.add_edge_data
        let get_edge_data = Ptree.get_edge_data
        let get_parent = Ptree.get_parent
        let get_neighs node tree = 
            match Ptree.get_node node tree with
            | Tree.Interior (_, b, c, d) -> [b; c; d]
            | Tree.Leaf (_, b) -> [b]
            | Tree.Single _ -> []
        let join x y tree = TreeOps.join_fn None [] x y tree 

        let get_root_costs tree = TreeOps.root_costs tree

        let break x tree = 
            let breakage = TreeOps.break_fn None x tree in
            let tree = 
                TreeOps.incremental_uppass breakage.Ptree.ptree
                breakage.Ptree.incremental 
            in
            tree, breakage.Ptree.tree_delta
        let reroot edge tree = 
            let a, b = TreeOps.reroot_fn None true edge tree in
            let a = TreeOps.incremental_uppass a b in
            a
        let downpass = TreeOps.downpass
        let uppass = TreeOps.uppass
        let branch_table = TreeOps.branch_table

        let of_string str data nodes =
            let tree = List.map (fun x -> None,x) (Tree.Parse.of_string str) in
            Sexpr.to_list (Build.prebuilt tree (data, nodes))

        let to_string collapse tree data = 
            let cost = string_of_float (Ptree.get_cost `Adjusted tree) in
            let res = 
                PtreeSearch.build_forest_with_names_n_costs collapse tree cost false None
            in
            List.map (AsciiTree.for_formatter false true true) res 

        let of_file file data nodes =
            let trees = List.map (fun x -> None,x) (Tree.Parse.of_file (`Local file)) in
            Sexpr.to_list (Build.prebuilt trees (data, nodes))

        let of_nodes data nodes = 
            let codes = List.map Node.taxon_code nodes in
            let codes = List.map string_of_int codes in
            let codes = List.map (fun x -> "(" ^ x ^ ")") codes in
            match of_string (String.concat " " codes) data nodes with
            | [h] -> h
            | [] -> failwith "No nodes"
            | _ -> assert false
        class execute_f f = object inherit [a, b] Sampler.do_nothing
            method evaluate tree = f tree
        end

        let local_neigh neigh f data tree =
            let sampler = new execute_f f in
            let res = 
                let optarg = 
                    let neigh = `ChainNeighborhoods neigh in
                    { PoyCommand.swap_default with Methods.ss = neigh }
                in
                PTS.find_local_optimum 
                ~base_sampler:sampler data (Sampler.create ()) 
                (`Single tree)
                (Lazy.lazy_from_val (All_sets.IntSet.empty)) 
                (`LocalOptimum optarg)
            in
            Sexpr.to_list res

        let spr = local_neigh `Spr
        let tbr = local_neigh `Tbr

        let build data nodes = 
            Sexpr.to_list 
                (Build.build_initial_trees `Empty data nodes 
                    (`Build (1, (`Wagner_Rnd (1, 0.,`Last, [], `UnionBased None)), [],(`MaxCount 20,`JoinDelta))))
    end


    module Runtime = struct
        type phylogeny = (a, b) Ptree.p_tree
        let get_cost cmp () = 
            let run = get_console_run () in
            Sexpr.fold_left (fun acc x ->
                match acc with
                | None -> Some (PhyloTree.get_cost x)
                | Some y -> Some (cmp y (PhyloTree.get_cost x))) 
            None run.trees
        let min_cost () = get_cost min ()
        let max_cost () = get_cost max ()
        let trees () =
            let run = get_console_run () in
            Sexpr.to_list run.trees
        let set_trees trees =
            let trees = Sexpr.of_list trees in
            let r = get_console_run () in
            set_console_run { r with trees = trees }
        let all_costs () = 
            let lst = List.rev_map PhyloTree.get_cost (trees ()) in
            List.sort compare lst

        let data () = 
            let run = get_console_run () in
            run.data

        let nodes () = 
            let run = get_console_run () in
            run.nodes

        let to_string bool =
            let run = get_console_run () in
            let trees = trees () in
            List.map (fun x -> PhyloTree.to_string bool x run.data) trees

        let of_string x = 
            let run = get_console_run () in
            let trees = PhyloTree.of_string x run.data run.nodes in
            let trees = Sexpr.of_list trees in
            console_run_val := { run with trees = trees }
    end

    module Node = Node
end

module FILES = struct

    let explode str = Parser.Wildcard.explode_filenames [(`Local str)]

    let do_ch f str = str, (f str)

    let open_all_in str = List.map (do_ch open_in) (explode str)

    let open_all_out str = List.map (do_ch open_out) (explode str)

    let run_n_close file f =
        let ch = open_in file in
        try
            let res = f ch in
            close_in ch;
            res
        with
        | err ->
                close_in ch;
                raise err
end

module DNA = struct

    let alph = Alphabet.nucleotides

    module CM = struct
        type cm2 = Cost_matrix.Two_D.m
        let of_list l = Cost_matrix.Two_D.of_list ~use_comb:true l 31
        let of_array arr = 
            let lst = Array.map (Array.to_list) arr in
            let lst = Array.to_list lst in
            of_list lst

        let of_sub_indel s i = 
            Cost_matrix.Two_D.of_transformations_and_gaps true 5 s i 31

        let of_sub_indel_affine a b c = 
            let mt,_ = of_sub_indel a b in
            Cost_matrix.Two_D.set_affine mt (Cost_matrix.Affine c);
            mt

        let of_file file = 
            let ch = new FileStream.file_reader (`Local file) in
            try 
                let res, _, _  = 
                    Cost_matrix.Two_D.of_channel 
                        ~orientation:false ~use_comb:true 31 ch
                in
                ch#close_in;
                res
            with 
            | err -> 
                    ch#close_in;
                    raise err
        let all_ones () = of_sub_indel 1 1

        let char_to_base a = 
            Alphabet.match_base (Char.escaped a) Alphabet.nucleotides 

        let median a b cm =
            let a = char_to_base a
            and b = char_to_base b in
            let m = Cost_matrix.Two_D.median a b cm in
            (Alphabet.find_code m Alphabet.nucleotides).[0]

        let cost a b cm =
            let a = char_to_base a
            and b = char_to_base b in
            Cost_matrix.Two_D.cost a b cm
    end

    module Seq = struct
        type s = Sequence.s
        type base = char
        let gap = '_'
        let of_string str = 
            Sequence.of_string str alph
        let to_string s =
            Sequence.to_string s alph
        let length = Sequence.length

        let check_bounds s pos =
            if pos >= 0 && pos < Sequence.length s then ()
            else failwith "index out of bounds"


        let get a pos =
            check_bounds a pos;
            (Alphabet.find_code (Sequence.get a pos) alph).[0]

        let set a b pos =
            check_bounds a pos;
            let pb = Sequence.get a pos 
            and code = Alphabet.match_base (Char.escaped b) alph in
            if pb = code then a
            else 
                let ns = Sequence.clone a in
                let _ = Sequence.set ns pos code in
                ns
        let prepend s c = 
            let code = Alphabet.match_base (Char.escaped c) alph 
            and cap = Sequence.capacity s in
            let new_cap =
                if cap > Sequence.length s then cap
                else 2 * cap
            in
            let ns = Sequence.create new_cap in
            let _ = 
                for i = cap - 1 downto 0 do
                    Sequence.prepend ns (Sequence.get s i);
                done
            in
            Sequence.prepend ns code;
            ns

        let merge a b = 
            let la = length a 
            and lb = length b in
            let ns = Sequence.create (la + lb) in
            for i = la - 1 downto 0 do
                Sequence.prepend ns (Sequence.get a i);
            done;
            for i = lb - 1 downto 0 do
                Sequence.prepend ns (Sequence.get b i);
            done;
            ns

        let delete a pos =
            check_bounds a pos;
            let la = Sequence.length a in
            let ns = Sequence.create la in
            for i = la - 1 downto 0 do
                if i <> pos then Sequence.prepend ns (Sequence.get a i)
                else ()
            done;
            ns

        let slice s a b =
            check_bounds s a;
            check_bounds s b;
            if a <= b then
                let ns = Sequence.create (1 + b - a) in
                let _ =
                    for i = b downto a do
                        Sequence.prepend ns (Sequence.get s i);
                    done
                in
                ns
            else failwith "Illegal slice"

        let reverse a = Sequence.safe_reverse a

        let complement a = Sequence.complement Alphabet.nucleotides a
    end

    module Fasta = struct
        type seqs = (string * Sequence.s) list
        type multi_seqs = seqs list

        let of_channel ?(respect_case = false) prealigned ch = 
            let filter (lst, txn) =
                let lst = List.flatten (List.flatten lst) in
                match lst with
                | [] -> false
                | _ -> true
            in
            let converter (lst, txn) =
                let lst = List.flatten (List.flatten lst) in
                match lst with
                | [seq] -> txn, seq
                | _ -> failwith "Illegal FASTA format"
            in
            let alph = 
                if prealigned then 
                    FileContents.Prealigned_Alphabet (Alphabet.nucleotides)
                else FileContents.Nucleic_Acids
            in
            let res = Fasta.of_channel ~respect_case:respect_case alph ch in
            List.map converter (List.filter filter res)

        let multi_of_channel ch = 
            let rec merger name lst1 lst2 =
                match lst1, lst2 with
                | h1 :: t1, h2 :: t2 -> 
                        ((name, h2) :: h1) :: (merger name t1 t2)
                | [], _ :: _ -> 
                        List.map (fun x -> [name, x]) lst2
                | [], [] -> []
                | _ -> failwith "Illegal multi-sequence file."
            in
            let converter acc (lst, txn) =
                List.fold_left (merger txn) acc (List.flatten lst)
            in
            let res = Fasta.of_channel FileContents.Nucleic_Acids ch in
            List.fold_left converter [] res

        let to_channel ch seqs =
            let converter (lst, txn) = (txn, lst) in
            Fasta.to_channel ch (List.map converter seqs) alph

        let of_file prealigned str = 
            FILES.run_n_close str (of_channel prealigned)

        let multi_of_file str =
            FILES.run_n_close str multi_of_channel

        let to_file str seqs = 
            let ch = open_out str in
            try 
                to_channel ch seqs;
                close_out ch;
            with
            | err ->
                    close_out ch;
                    raise err

        let random_sample seqs name_f samples samplesize =
            let seqs = Array.of_list seqs in
            if Array.length seqs < samplesize then 
                failwith "Not enough sequences"
            else
                for i = 1 to samples do
                    let res = ref [] in
                    Array_ops.randomize seqs;
                    for i = samplesize - 1 downto 0 do
                        res := seqs.(i) :: !res;
                    done;
                    to_file (name_f ()) !res;
                done

        let multi_sample seqs which name_f samples samplesize =
            let seqs = List.nth seqs which in
            random_sample seqs name_f samples samplesize


        let print_sequence s =
            Printf.printf "%s\n%!" (Seq.to_string s)

        let _ = Callback.register "print failed alignment" print_sequence
    end

    module Generic = struct
        let molecular file =
            Fasta.of_channel false (Parser.Files.molecular_to_fasta (`Local file))
    end

    module Align = struct
        let algn s1 s2 cm =
            let s1', s2', c = 
                Sequence.Align.align_2 ~first_gap:true s1 s2 cm Matrix.default
            in
            let median = Sequence.median_2 s1' s2' cm in
            s1', s2', c, median

        let gen_algn_all f seqs cm =
            List.map (fun (t1, s1) ->
                List.map (fun (t2, s2) -> 
                    let s1', s2', c, al = f s1 s2 cm in
                    (t1, s1'), (t2, s2'), c, al) seqs) seqs

        let algn_and_print s1 s2 cm =
            let s1', s2', c, m = algn s1 s2 cm in
            (Seq.to_string s1'), (Seq.to_string s2'), c, (Seq.to_string m)

        let algn_all = gen_algn_all algn 

        let algn_all_and_print = gen_algn_all algn_and_print

    end

end
