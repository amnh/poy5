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

module F : GraphTree.GRAPHICS_TYPE = struct
    type color = Graphps.color
    let black = Graphps.black
    let close_graph = Graphps.close_graph
    let draw_string = Graphps.draw_string
    let foreground = Graphps.foreground
    let lineto = Graphps.lineto
    let moveto = Graphps.moveto
    let open_graph = Graphps.open_graph
    let plot = Graphps.plot
    let red = Graphps.red
    let set_color = Graphps.set_color
    let size_x = Graphps.size_x
    let size_y = Graphps.size_y
    let text_size = Graphps.text_size
    let display () = ()

end

module GraphTreePs = GraphTree.Make (F)

let d = ref 0
let max_depth = ref 0
let num_leaves = ref 0
let longest_name = ref 0

let reset () = 
    d := 0;
    max_depth := 0;
    num_leaves := 0;
    longest_name := 0

let draw_file title ?(filename = "poy_tree.ps") t =
    Graphps.open_ps filename;
    let str = String.make !longest_name 'm' in
    let size1 = GraphTree.leaf_distance * (!num_leaves + 4) in
    let size2 = (GraphTree.depth_distance * (!max_depth + 1)) + 
        (let x, _ = (F.text_size str) in x) in
    let size = (string_of_int size2) ^ "x" ^ (string_of_int size1) in
    Graphps.open_graph size;
    let x, y = F.text_size title in
    F.moveto 10 (size1 - (2 * y));
    F.draw_string title;
    F.moveto 0 0;
    GraphTreePs.draw t;
    Graphps.close_graph ()

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
    let cat inch outch =
        try
            while true do
                output_char outch (input_char inch);
            done
        with
        | End_of_file -> ()
    in
    (* Print each tree in a temporary file *)
    let files = Array.map (fun (cost, t) ->
        let title = 
            if "" = title then Printf.sprintf "Tree Cost: %.3f" cost
            else title
        in
        let filename = Filename.temp_file "poy" ".ps" in
        let _ = draw_file title ~filename:filename t in
        filename) all_trees
    in
    (* Create the output tree and append the .ps extension *)
    let filename = 
        (try (Filename.chop_extension filename) with
        | Invalid_argument _ -> filename) ^ ".ps"
    in
    (* Open the final resulting file and concatenate all the temporary 
    * files *)
    let outch = open_out filename in
    Array.iter (fun x ->
        let inch = open_in x in
        cat inch outch;
        close_in inch;) files;
    close_out outch;
    reset ()


    
