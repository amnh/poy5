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

let () = SadmanOutput.register "GraphTree" "$Revision"

let leaf_distance = 20
let depth_distance = 20

module type GRAPHICS_TYPE = sig
    type color = int
    val black : color
    val close_graph : unit -> unit
    val draw_string : string -> unit
    val foreground : color 
    val lineto : int -> int -> unit 
    val moveto : int -> int -> unit 
    val open_graph : string -> unit
    val plot : int -> int -> unit 
    val red : color
    val set_color : color -> unit 
    val size_x : unit -> int
    val size_y : unit -> int
    val text_size : string -> int * int
    val display : unit -> unit
end


exception Wrong_format of string


let get_code ch =
    let codeLine = input_line ch in
    if Str.string_match (Str.regexp " *\\([0-9]+\\)") codeLine 0 then
        (Str.matched_group 1 codeLine)
    else
        ""


module Make = functor (G : GRAPHICS_TYPE) -> struct

    (** [draw_edges from_point to_point] draws an edge from the from_point which is
    * an (int, int) to the to_point *)
    let draw_edges from_point to_point =
        let x1, y1 = from_point and
        x2, y2 = to_point in
        G.moveto x1 y1;
        G.lineto x1 y2;
        G.moveto x1 y2;
        G.lineto x2 y2;
        ()


    (** [average pointList] takes a list of points and calculates the average y
    * position, it also calculates the minimun x position *)    
    let average pointList =
        let total = ref 0 and
        min_x_point = ref 100000 and
        n = List.length pointList in 
        for i = 0 to n - 1 do
            let x, y = List.nth pointList i in
            total := !total + y;
            if x < !min_x_point then
                min_x_point := x;
        done;
        (min_x_point , (!total / n) )


    (** [calc_depth_leaves tree depth max_depth num_leaves longest_name]
     *  calculates 1) the maximum depth of the tree (max_depth variable)
     *             2) the number of leaves (num_leaves variable)
     *             3) the longest taxon name (longest_name variable)
     *)
    let rec calc_depth_leaves t depth max_depth num_leaves longest_name =
        match t with
        | Parser.Tree.Node (y, _) -> 
                incr depth;
                for i = 0 to (List.length y) - 1 do
                    calc_depth_leaves (List.nth y i) depth max_depth num_leaves
                    longest_name;
                done;
                decr depth;
        | Parser.Tree.Leaf y -> 
                incr depth;
                if !depth > !max_depth then max_depth := !depth;
                decr depth;
                let strLength = String.length y in
                if strLength > !longest_name then
                    longest_name := strLength;
                incr num_leaves

    (** [draw_tree t] takes a t which is a Parser.Tree.t type and draws a tree 
    *   also have optional arguments title, size and leafColor - The graph
    *   should have already been open.
    *   title - takes a string 
    *   size - takes a string of format " 800x600" - notice must start with space
    *   leafColor - takes a G color constant use as G.red 
    *   choices are:
    *   black, red , blue, green, yellow, cyan, magenta
    *   *)    
    let draw ?(size="") ?(leafColor=G.black) t =
        let t = AsciiTree.sort_tree t in
       let d = ref 0 and
       max_depth = ref 0 and
       num_leaves = ref 0 and
       longest_name = ref 0 in
       calc_depth_leaves t d max_depth num_leaves longest_name;
       let fontWidth, fontHeight = G.text_size "A" in
       let deltaY = leaf_distance and
       deltaX = depth_distance in
       let total_height = leaf_distance * !num_leaves in
       let update_counter =
            let counter = ref 0 in
            fun () -> 
                incr counter;
                !counter
       in        
       let rec coord depth deltaX deltaY tree =
           match tree with
           | Parser.Tree.Leaf name ->
                   let row = update_counter () in
                   let y = total_height - ((row - 1) * deltaY)  and
                   x = (!max_depth * deltaX + 10) in
                   G.moveto (x + fontWidth) (y - fontHeight / 2);
                   let prev_color = G.foreground in
                   G.set_color leafColor;
                   G.draw_string name;
                   G.set_color prev_color;
                   (x, y)
           | Parser.Tree.Node (children, name) ->
                   let coord_children =
                       List.map (coord (depth+1) deltaX deltaY) children
                   in
                   let xTemp , avg = average coord_children in
                   let x = !xTemp - deltaX in
                   let prev_color = G.foreground in
                   G.set_color leafColor;
                   G.plot x avg;
                   let name_len = String.length name in
                   let name = 
                       if name_len = 0 then name
                       else if name.[name_len - 1] = '.' then 
                           String.sub name 0 (name_len - 1) 
                       else name
                   in
                   G.moveto (x - (fontWidth * (name_len)) - 1)
                   (avg + fontHeight/2);
                   G.draw_string name;
                   G.set_color prev_color;
                   List.iter (draw_edges (x, avg)) coord_children;
                   G.plot x avg;
                   (x, avg)
       in
       let _ = coord 1 deltaX deltaY t in ()


       let disp_tree str tree =
           try
               try draw ~size:" 800x800" tree with 
               | _ -> begin
                   if str <> ""
                   then print_endline str;
                   AsciiTree.draw false stdout tree
               end
           with
           | e -> ()

end

