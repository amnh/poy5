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

let () = SadmanOutput.register "GraphicsPs" "$Revision: 1616 $"

module Old : GraphTree.GRAPHICS_TYPE = struct
    type color = Graphps.color
    let black = Graphps.black
    let close_graph = Graphps.close_graph
    let draw_string = Graphps.draw_string
    let foreground = Graphps.foreground
    let lineto = Graphps.lineto
    let moveto = Graphps.moveto
    let polyline lst = 
        List.iter (fun (a, b) ->
            lineto a b;
            moveto a b) lst
    let open_graph = Graphps.open_graph
    let plot = Graphps.plot
    let red = Graphps.red
    let set_color = Graphps.set_color
    let size_x = Graphps.size_x
    let size_y = Graphps.size_y
    let text_size = Graphps.text_size
    let display () = ()
    let add_page () = ()
    let open_file str = Graphps.open_ps str
end

module Pdf : GraphTree.GRAPHICS_TYPE = struct
    type color = Graphicpdf.color
    let black = Graphicpdf.black
    let close_graph = Graphicpdf.close_graph
    let draw_string = Graphicpdf.draw_string
    let foreground = Graphicpdf.foreground
    let lineto = Graphicpdf.lineto
    let moveto = Graphicpdf.moveto
    let polyline = Graphicpdf.polyline
    let open_graph = Graphicpdf.open_graph
    let plot = Graphicpdf.plot
    let red = Graphicpdf.red
    let set_color = Graphicpdf.set_color
    let size_x = Graphicpdf.size_x
    let size_y = Graphicpdf.size_y
    let text_size = Graphicpdf.text_size
    let display = Graphicpdf.display
    let add_page = Graphicpdf.add_page
    let open_file = Graphicpdf.open_file
end

module Ps = Pdf

module GraphTreePs = GraphTree.Make (Ps)

let d = ref 0
let max_depth = ref 0
let num_leaves = ref 0
let longest_name = ref 0

let reset () = 
    d := 0;
    max_depth := 0;
    num_leaves := 0;
    longest_name := 0


let draw_file title t =
    let str = String.make !longest_name 'm' in
    let size1 = GraphTree.leaf_distance * (!num_leaves + 4) in
    let size2 = (GraphTree.depth_distance * (!max_depth + 1)) + 
        (let x, _ = (Ps.text_size str) in x) in
    let size = (string_of_int size2) ^ "x" ^ (string_of_int size1) in
    Ps.open_graph size;
    let x, y = Ps.text_size title in
    Ps.moveto 10 (size1 - (8 * y));
    Ps.draw_string title;
    Ps.moveto 0 0;
    GraphTreePs.draw t

let display title filename all_trees = 
    (* A function to concatenate files, appending the contents of inch
    * to outch *)
    Array.iter (fun (_, t) ->
        d := 0;
        let max_leaves = !num_leaves in
        num_leaves := 0;
        GraphTreePs.calc_depth_leaves t d max_depth num_leaves longest_name;
        if !num_leaves < max_leaves then num_leaves := max_leaves;
    )
    all_trees;
    (* Print each tree in a temporary file *)
    let filename = 
        (try Filename.chop_extension filename with
        | _ -> filename) ^ ".pdf"
    in
    Ps.open_file filename;
    let len = Array.length all_trees in
    Array.iteri (fun pos (cost, t) ->
        let title = 
            if "" = title then Printf.sprintf "Tree Cost: %.3f" cost
            else title
        in
        draw_file title t;
        if pos < len - 1 then Ps.add_page ()) all_trees;
    Ps.close_graph ();
    reset ()
