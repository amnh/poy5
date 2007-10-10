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

let () = SadmanOutput.register "Graphviz" "$Revision: 1644 $"

type agraph
type anode
type aedge
type gvc
type t

(*list of functions from ocaml to c
*                 returns          input
*                 ->------        ->-----
* gvcontext       gvc              none
* gvParseArgs     none             gvc, argc, argv 
* agopen          Agraph           string, int 
* agsubg          Agraph           Agraph, string
* agraphattr      none             Agraph, string, string
* agnode          Agnode           Agraph, string
* agedge          Agedge           Agraph, Agnode, Agnode
* agnodeattr      node             Agraph, string, string
* agset           none             Anode, string, string
* gvLayoutJobs
* gvRenderJobs
* gvFreeLayout
* agclose
*)

external c_gvContext : unit -> gvc = "graphviz_CAML_gvContext"

external c_gvParseArgs : gvc ->  string  -> string -> unit =
    "graphviz_CAML_gvParseArgs" 

external c_agopen : string -> int -> agraph = "graphviz_CAML_agopen"

external c_agsubg : agraph -> string -> agraph = "graphviz_CAML_agsubg"

external c_agraphattr : agraph -> string -> string -> unit =
    "graphviz_CAML_agraphattr"

external c_agnode : agraph -> string -> anode = "graphviz_CAML_agnode"                  
external c_agedge : agraph -> anode -> anode -> unit = "graphviz_CAML_agedge"                  
external c_agnodeattr : agraph -> string-> string -> unit  =
    "graphviz_CAML_agnodeattr"                  

external c_agedgeattr : agraph -> string-> string -> unit  =
    "graphviz_CAML_agedgeattr"                 

external c_agset : anode -> string-> string -> unit  =
    "graphviz_CAML_agset"                 

external c_gvLayoutJobs : gvc -> agraph-> unit  =
    "graphviz_CAML_gvLayoutJobs"                

external c_gvRenderJobs : gvc -> agraph-> unit  =
    "graphviz_CAML_gvRenderJobs"                

external c_gvFreeLayout : gvc -> agraph-> unit  =
    "graphviz_CAML_gvFreeLayout" 

external c_agclose : agraph-> unit  = "graphviz_CAML_agclose" 

external c_gvFreeContext : gvc -> unit  = "graphviz_CAML_gvFreeContext" 

let get_number x =
    let len = String.length x in
    let y = String.sub x 3 (len-3) in
    int_of_string y

let tree_to_nodes_leafs_edges t = 
    let node_generator = 
        let name = ref 0 in
        fun () -> 
            let res = "HTU" ^ (string_of_int !name) in
            incr name;
            res
    in      
    let leaf_generator = 
        let count = ref 0 in
        fun () -> 
            let res = "OTU" ^ (string_of_int !count) in
            incr count;
            res
    in 
    let tbl = Hashtbl.create 100 in
    let rec builder (leafs, edges, nodes) t =
        match t with
        | Parser.Tree.Leaf name -> 
                let code = leaf_generator () in
                Hashtbl.add tbl t code;
                ((code, name) :: leafs), edges, nodes
        | Parser.Tree.Node (others, _) ->
                let leafs, edges, nodes = 
                    List.fold_left builder (leafs, edges, nodes) others 
                in
                let code = node_generator () in
                Hashtbl.add tbl t code;
                let connect a b = 
                    a, Hashtbl.find tbl b
                in
                let new_edges = 
                    List.fold_left begin fun y x -> 
                        (connect code x) :: y end edges others
                in
                leafs, new_edges, (code :: nodes)
    in
    builder ([], [], []) t

                      

let rec count_nodes_leaves t num_nodes num_leaves=
   match t with
    | Parser.Tree.Node (y, _) -> 
            incr num_nodes;
            List.iter (fun x -> count_nodes_leaves x num_nodes num_leaves) y;
    | Parser.Tree.Leaf y -> incr num_leaves



(** 
*   t = Parser.tree format
*   output_format = -Txxx where xxx is the type of file produced - can be
*   png, bmp, gif, ...
*   output_file = -o<fileName>    , example -otree.png
*   graph_attr = (string * string) list of  attributes
*   node_attr = (string * string) list of  attributes
**)    
let make_graph t output_format output_file graph_attr node_attr leaf_attr = 
    let num_nodes = ref 0 and
    num_leaves = ref 0 in
    count_nodes_leaves t num_nodes num_leaves;
    let gcontext = c_gvContext () in
    c_gvParseArgs gcontext output_format output_file; 
    let graph = c_agopen "g" 1 in
    let subGraph = c_agsubg graph "taxon" in
    let leaves, edges, nodes = tree_to_nodes_leafs_edges t in

    let htu0 = c_agnode graph "HTU0" in
    let graph_nodes = Array.make !num_nodes htu0 and
    graph_leaves = Array.make !num_leaves htu0 in
    List.iter begin fun code ->
        let index = get_number code in
          graph_nodes.(index) <- c_agnode graph code;
        end nodes;
    List.iter begin fun (code, name) ->
        let index = get_number code in
          graph_leaves.(index) <- c_agnode subGraph name;
     end leaves;
    List.iter 
        begin fun (a, b) ->
            let indexA = get_number a and
            indexB = get_number b in
            let firstB = String.get b 0 in
            if firstB = 'O' then
                c_agedge graph graph_nodes.(indexA)
                graph_leaves.(indexB)
            else
                c_agedge graph graph_nodes.(indexA)
                graph_nodes.(indexB);
        end edges;
    List.iter begin fun (attrName, attrValue ) ->
          c_agraphattr graph attrName attrValue;
     end graph_attr;
    List.iter begin fun (attrName, attrValue ) ->
          c_agnodeattr graph attrName attrValue;
     end node_attr;
    List.iter begin fun (attrName, attrValue ) ->
    for i = 0 to !num_leaves - 1 do
          c_agset graph_leaves.(i) attrName attrValue;
    done;
     end leaf_attr;
    c_gvLayoutJobs gcontext graph;
    c_gvRenderJobs gcontext graph;
    c_gvFreeLayout gcontext graph;
    c_agclose graph;
    c_gvFreeContext gcontext


