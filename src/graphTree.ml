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
let depth_distance = 30

module type GRAPHICS_TYPE = sig
    type display
    type color = int
    val black : color
    val close_graph : display -> unit
    val draw_string : 
        ?tag:string -> ?link:string -> display -> string -> display
    val foreground : color 
    val lineto : ?tag:string -> ?link:string -> display -> int -> int -> display
    val moveto : display -> int -> int -> display
    val polyline : display -> (int * int * string option * string option) list -> display
    val open_file : string -> display
    val open_graph : display -> string -> display
    val plot : display -> int -> int -> display
    val red : color
    val set_color : display -> color -> display
    val size_x : display -> int
    val size_y : display -> int
    val text_size : string -> int * int
    val display : display -> unit
    val add_page : display -> display
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
    let draw_edges display from_point to_point =
        let x1, y1 = from_point and
        x2, y2 = to_point in
        let display = G.moveto display x1 y1 in
        G.polyline display [(x1, y1, None, None); (x1, y2, None, None); (x2, y2,
        None, None)]


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
                List.iter (fun t ->
                    calc_depth_leaves t depth max_depth num_leaves
                    longest_name) y;
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
    let draw ?(size="") ?(leafColor=G.black) display t =
        let t = AsciiTree.sort_tree t in
        let display = ref display in
       let d = ref 0 and
       max_depth = ref 0 and
       num_leaves = ref 0 and
       longest_name = ref 0 in
       calc_depth_leaves t d max_depth num_leaves longest_name;
       let fontWidth, fontHeight = G.text_size "m" in
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
                   display := 
                       G.moveto !display 
                       (x + fontWidth) (y - fontHeight / 2);
                   let prev_color = G.foreground in
                   display := G.set_color !display leafColor;
                   let tag = name in
                   display := G.draw_string ~tag !display name;
                   display := G.set_color !display prev_color;
                   (x, y)
           | Parser.Tree.Node (children, name) ->
                   let coord_children =
                       List.map (coord (depth+1) deltaX deltaY) children
                   in
                   let xTemp , avg = average coord_children in
                   let x = !xTemp - deltaX in
                   let prev_color = G.foreground in
                   display := G.set_color !display leafColor;
                   display := G.plot !display x avg;
                   let name_len = String.length name in
                   let name = 
                       if name_len = 0 then name
                       else if name.[name_len - 1] = '.' then 
                           String.sub name 0 (name_len - 1) 
                       else name
                   in
                   display := 
                       G.moveto !display (x - (fontWidth * (name_len)) - 3)
                   (avg + fontHeight/2);
                   display := G.draw_string !display name;
                   display := G.set_color !display prev_color;
                   let () =
                       match coord_children with
                       | [_] | [] -> failwith "I need at least two children"
                       | (a, b) :: tl ->
                               match List.rev tl with 
                               | (c, d) :: tl ->
                                       display := 
                                           G.polyline 
                                           !display
                                           [(a, b, None, None); (x, b, None,
                                           None); (x, d, None, None); (c, d,
                                           None, None)];
                                       List.iter (fun (c, d) -> 
                                           display := 
                                               G.polyline 
                                               !display
                                               [(x, d, None, None); 
                                               (c, d, None, None)]) tl
                               | _ -> assert false 
                   in
                   display := G.plot !display x avg;
                   (x, avg)
       in
       let _ = coord 1 deltaX deltaY t in 
       !display

       let rec get_children tree =
           match tree with
           | ("Forest", _, `Set lst) 
           | ("Tree", _, `Set lst) ->
                   List.fold_left (fun acc sexpr ->
                       Sexpr.fold_left (fun acc ((name, _, _) as t) ->
                           if name = "Tree" then t :: acc
                           else acc) acc sexpr) [] lst
           | ("Tree", _, `Single item) 
           | ("Forest", _, `Single item) -> get_children item
           | ("Tree", _, `Delayed item) 
           | ("Forest", _, `Delayed item) -> 
                   get_children (Sexpr.first (item ()))
           | (x, _, _) -> []

       let get_name (_, x, _) =
           Tags.value_to_string (List.assoc "Name" x) 

       let get_node (_, _, x) =
           let res = 
               match x with
               | `Set lst ->
                       List.fold_left (fun acc sexpr ->
                           Sexpr.fold_left (fun acc ((name, _, _) as t) ->
                               if name = "Node" then Some t 
                               else acc) acc sexpr) None lst
               | `Delayed f ->
                       Sexpr.fold_left (fun acc ((name, _, _) as t) ->
                           if name = "Node" then Some t 
                           else acc) None (f ())
               | `Single x -> Some x
               | _ -> assert false
           in
           match res with
           | None -> assert false 
           | Some x -> x

       let rec calc_depth_leaves_diag t depth max_depth num_leaves longest_name =
           match get_children t with
           | [] ->
                   let y = 
                       let t = get_node t in
                       get_name t 
                    in
                   incr depth;
                   if !depth > !max_depth then max_depth := !depth;
                   decr depth;
                   let strLength = String.length y in
                   if strLength > !longest_name then
                       longest_name := strLength;
                   incr num_leaves
           | y ->
                   incr depth;
                   List.iter (fun t ->
                       calc_depth_leaves_diag t depth max_depth num_leaves
                       longest_name) y;
                   decr depth

    (** draw_diagnosis *)    
    let draw_diagnosis ?(prefix="") ?(size="") ?(leafColor=G.black) 
    display (t : Tags.output) =
        let nodes = ref [] in
        let t =
            match t with
            | ("Forest", _, `Set [t]) -> Sexpr.first t
            | _ -> assert false
        in
       let display = ref display in
       let d = ref 0 and
       max_depth = ref 0 and
       num_leaves = ref 0 and
       longest_name = ref 0 in
       calc_depth_leaves_diag t d max_depth num_leaves longest_name;
       let fontWidth, fontHeight = G.text_size "m" in
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
           let node = get_node tree in
           nodes := node :: !nodes;
           let name = get_name node in
           match get_children tree with
           | [] ->
                   let row = update_counter () in
                   let y = total_height - ((row - 1) * deltaY)  and
                   x = (!max_depth * deltaX + 10) in
                   display := 
                       G.moveto !display 
                       (x + fontWidth) (y - fontHeight / 2);
                   let prev_color = G.foreground in
                   display := G.set_color !display leafColor;
                   let link = prefix ^ ":" ^ name in
                   display := G.draw_string ~link !display name;
                   display := G.set_color !display prev_color;
                   (x, y)
           | children ->
                   let coord_children =
                       List.map (coord (depth+1) deltaX deltaY) children
                   in
                   let xTemp , avg = average coord_children in
                   let x = !xTemp - deltaX in
                   let prev_color = G.foreground in
                   display := G.set_color !display leafColor;
                   display := G.plot !display x avg;
                   let name_len = String.length name in
                   let link = prefix ^ ":" ^ name in
                   display := 
                       G.moveto !display (x - (fontWidth * (name_len)) - 3)
                   (avg + fontHeight/2);
                   display := G.draw_string ~link !display name;
                   display := G.set_color !display prev_color;
                   let () =
                       match coord_children with
                       | [_] | [] -> failwith "I need at least two children"
                       | (a, b) :: tl ->
                               match List.rev tl with 
                               | (c, d) :: tl ->
                                       display := 
                                           G.polyline 
                                           !display
                                           [(a, b, None, None); (x, b, None,
                                           None); (x, d, None, None); (c, d,
                                           None, None)];
                                       List.iter (fun (c, d) -> 
                                           display := 
                                               G.polyline 
                                               !display
                                               [(x, d, None, None); 
                                               (c, d, None, None)]) tl
                               | _ -> assert false 
                   in
                   display := G.plot !display x avg;
                   (x, avg)
       in
       let _ = coord 1 deltaX deltaY t in 
       !display, !nodes

       let disp_tree display str tree =
           try
               try draw ~size:" 800x800" display tree with 
               | _ -> begin
                   if str <> ""
                   then print_endline str;
                   AsciiTree.draw false stdout tree;
                   display
               end
           with
           | e -> display

end

