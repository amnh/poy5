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

module Pdf : GraphTree.GRAPHICS_TYPE with type display = Graphicpdf.display = struct
    type display = Graphicpdf.display
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


let aux_draw_file to_draw display title t =
    let str = String.make !longest_name 'm' in
    let size1 = GraphTree.leaf_distance * (!num_leaves + 4) in
    let size2 = (GraphTree.depth_distance * (!max_depth + 1)) + 
        (let x, _ = (Ps.text_size str) in x) in
    let size = (string_of_int size2) ^ "x" ^ (string_of_int size1) in
    let display = Ps.open_graph display size in
    let x, y = Ps.text_size title in
    let display = Ps.moveto display 10 (size1 - (8 * y)) in
    let display = Ps.draw_string display title in
    let display = Ps.moveto display 0 0 in
    to_draw display t

let draw_file display title t =
    aux_draw_file GraphTreePs.draw display title t

let draw_diagnosis prefix display title t =
    aux_draw_file (GraphTreePs.draw_diagnosis ~prefix) display title t

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
    let display = ref (Ps.open_file filename) in
    let len = Array.length all_trees in
    Array.iteri (fun pos (cost, t) ->
        let title = 
            if "" = title then Printf.sprintf "Tree Cost: %.3f" cost
            else title
        in
        display := draw_file !display title t;
        if pos < len - 1 then 
            display := Ps.add_page !display) all_trees;
    Ps.close_graph !display;
    reset ()

let translate curr_width height width acc (x, y) =
    height := min !height (!height +. y);
    curr_width := !curr_width +. x;
    width := max !width !curr_width;
    Pdfpages.Op_cm (Transform.matrix_of_transform 
    [Transform.Translate (x,y)]) :: acc

let display_table translate has_header hspan vspan operators table = 
    (* The table must contain strings inside *)
    let height = Array.length table in
    let width = 
        if height > 0 then Array.length table.(0)
        else 0
    in
    let print_row acc arr =
        let acc =
            Array.fold_left (fun acc x ->
                translate 
                (Pdfpages.Op_ET :: Pdfpages.Op_Tj x :: 
                Pdfpages.Op_Tf ("/F0", 10.) ::
                    Pdfpages.Op_BT :: acc) (hspan, 0.)) acc arr
        in
        translate acc ((-. 1.) *. (float_of_int width) *. hspan,
        (-.1.) *. vspan)
    in
    let print_table acc arr = Array.fold_left print_row acc arr in
    print_table operators table


let display_node prefix display node =
    let node_name = ref None in
    let display = Graphicpdf.add_page display in
    let height = ref (-. 20.)
    and width = ref 20.
    and curr_width = ref 20. in
    let add_name translation acc name =
        let x, _ = Graphicpdf.text_size name in
        width := max !width (!curr_width +. (float_of_int x));
        Pdfpages.Op_ET :: Pdfpages.Op_Tj name :: Pdfpages.Op_Tf ("/F0", 10.) ::
            Pdfpages.Op_BT :: translate curr_width height width acc translation
    in
    let one_attribute (name, contents) =
        if name = "Name" && None = !node_name then
            node_name := Some (Tags.value_to_string contents)
        else ();
        name ^ " = \"" ^ Tags.value_to_string contents ^ "\" "
    in
    let make_attributes attributes = 
        String.concat " " (List.map one_attribute attributes)
    in
    let add_contents = add_name (0., 0.) in
    let rec aux_print_it_all acc (name, attributes, children) =
        let acc = add_name (10., -10.) acc name in
        let acc = translate curr_width height width acc (10., -.10.) in
        let acc = add_name (0., 0.) acc (make_attributes attributes)  in
        let acc = translate curr_width height width acc (0., -.10.) in
        let acc =
            match children with
            | #Tags.unstructured as v -> add_contents acc (Tags.value_to_string v)
            | #Tags.structured as v ->
                    (match v with
                    | `Empty -> acc
                    | `Single x -> aux_print_it_all acc x
                    | `Set x -> List.fold_left (fun acc x -> 
                            Sexpr.fold_left (aux_print_it_all) acc x) acc x
                    | `Delayed x ->
                            Sexpr.fold_left (aux_print_it_all) acc (x ()))
            | `CDATA _ -> failwith "GraphPS can't handle CDATA yet"
        in
        translate curr_width height width acc (-.20., 0.)
    in
    let trans = (0., (-1.) *. (!height +. 20.) ) in
    let ops = aux_print_it_all [] node in
    let ops = translate curr_width height width (List.rev ops) trans in
    let display = 
        { display with Graphicpdf.page_text = ops; max_x = !width; 
        max_y = (!height); current_x = !width; current_y = !height} in
    match !node_name with
    | None -> display 
    | Some name ->
            Graphicpdf.add_reference display (0., (-.1.) *. !height) (prefix ^ ":"
            ^ name)


let display_diagnosis title filename all_trees = 
    (* A function to concatenate files, appending the contents of inch
    * to outch *)
    let all_trees = Array.of_list (Sexpr.to_list all_trees) in
    Array.iter (fun t ->
        d := 0;
        let max_leaves = !num_leaves in
        num_leaves := 0;
        GraphTreePs.calc_depth_leaves_diag t d max_depth num_leaves longest_name;
        if !num_leaves < max_leaves then num_leaves := max_leaves;
    )
    all_trees;
    (* Print each tree in a temporary file *)
    let filename = 
        (try Filename.chop_extension filename with
        | _ -> filename) ^ ".pdf"
    in
    let display = ref (Ps.open_file filename) in
    let len = Array.length all_trees in
    Array.iteri (fun pos t ->
        let prefix = "tree" ^ (string_of_int pos) in
        let new_display, nodes = draw_diagnosis prefix !display "" t in
        display := List.fold_left (display_node prefix) new_display nodes;
        if pos < len - 1 then 
            display := Ps.add_page !display) all_trees;
    Ps.close_graph !display;
    reset ()
