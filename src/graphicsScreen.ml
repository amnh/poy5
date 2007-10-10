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

let () = SadmanOutput.register "GraphicsScreen" "$Revision: 1616 $"

module F : GraphTree.GRAPHICS_TYPE = struct
    type color = Graphics.color
    let black = Graphics.black
    let close_graph = Graphics.close_graph
    let draw_string = Graphics.draw_string
    let foreground = Graphics.foreground
    let lineto = Graphics.lineto
    let moveto = Graphics.moveto

    let open_graph string = 
        Graphics.open_graph "";
        Graphics.set_window_title string

    let plot = Graphics.plot
    let red = Graphics.red
    let set_color = Graphics.set_color
    let size_x = Graphics.size_x
    let size_y = Graphics.size_y
    let text_size = Graphics.text_size
    let display ()= ()

end

module GraphTreeScreen = GraphTree.Make (F)

let draw_screen ?(title="") t = 
    Graphics.set_window_title title;
    GraphTreeScreen.draw t

module GraphOfXML = struct
    let gather_info info_list_additive info_list_nonadditive ch =
        try
            let min = ref "" 
            and code = ref ""
            and cost = ref "" 
            and non_add_code = ref ""
            and elt = ref "" 
            and line = ref "" in
            line := input_line ch ;
            while not (Str.string_match (Str.regexp " *</final>") !line 0) do
                if  Str.string_match (Str.regexp 
                "<[a-zA-Z]* code=\"\\([0-9]+\\)\" cost=\"\\([0-9]+\\)") 
                !line 0 then
                    begin
                        code := Str.matched_group 1 !line;
                        cost := Str.matched_group 2 !line;
                    end
                else if  Str.string_match (Str.regexp 
                "<[a-zA-Z:]* code=\"\\([0-9]+\\)\"><elt>\\([0-9; ]+\\)") 
                !line 0 then
                    begin
                        non_add_code := Str.matched_group 1 !line;
                        elt := Str.matched_group 2 !line;
                        info_list_nonadditive := (!non_add_code, !elt) ::
                            !info_list_nonadditive;
                    end
                else if Str.string_match 
                (Str.regexp " *<min>\\([0-9]+\\)") !line 0 then
                    min := Str.matched_group 1 !line
                else if  Str.string_match 
                (Str.regexp " *<max>\\([0-9]+\\)") !line 0 then
                    info_list_additive := (!code, !cost, !min, 
                    (Str.matched_group 1 !line)):: !info_list_additive;
                line := input_line ch;
            done;
            List.rev !info_list_additive, List.rev !info_list_nonadditive;
         with
         |End_of_file -> List.rev !info_list_additive, List.rev
         !info_list_nonadditive

    let find_closest_point mouse_x mouse_y pointlist = 
        let min_distance = ref 0 
        and min_code = ref "" in
        let x, y, code = List.nth !pointlist 0 in
        min_code := code;
        min_distance := abs (mouse_x - x) + abs (mouse_y - y);
        for i = 0 to List.length !pointlist - 1 do
            let x, y, code = List.nth !pointlist i in
            let distance = abs (mouse_x - x) + abs (mouse_y - y) in
            if distance < !min_distance then
                begin
                   min_distance := distance;
                   min_code := code;
                end
        done;
        !min_code
                    
    let print_info_verbose (code, cost, min, max) =
        print_string ("code = " ^ code);
        print_string (" cost = " ^ cost);
        print_string (" min = " ^ min);
        print_string (" max = " ^ max);
        print_newline ()

    let print_info (code, cost, min, max) =
        let semicolon = "; " 
        and comma = ", " in
        print_string ("(" ^ code ^ comma ^ cost ^ 
        comma ^ min ^ comma ^ max ^ ")" ^ semicolon )

    let print_info_nonadd_verbose (code, elt) =
        print_string ("code = " ^ code);
        if String.contains elt ';' then
            print_string (" elt = " ^  "[" ^ elt ^ "]" )
        else
            print_string (" elt = " ^ elt);
        print_newline ()

    let print_info_nonadd (code, elt) =
        let semicolon = "; " 
        and comma = ", " in
        if String.contains elt ';' then
            print_string ("(" ^ code ^ comma ^ "[" ^ elt ^ "]" ^ ")" ^ semicolon )
        else
            print_string ("(" ^ code ^ comma ^ elt ^ ")" ^ semicolon )

    (** [xml_calc_depth_leaves tree depth max_depth num_leaves longest_name]
     *  calculates 1) the maximum depth of the tree (max_depth variable)
     *             2) the number of leaves (num_leaves variable)
     *             3) the longest taxon name (longest_name variable)
     *)
    let rec xml_calc_depth_leaves t depth max_depth num_leaves longest_name =
       match t with
        | Parser.Tree.Node (y, _) -> 
                incr depth;
                for i = 0 to (List.length y) - 1 do
                    xml_calc_depth_leaves (List.nth y i) depth max_depth num_leaves
                    longest_name;
                done;
                decr depth;
                
        | Parser.Tree.Leaf y -> 
                incr depth;
                if !depth > !max_depth then
                    max_depth := !depth;
                decr depth;
                let strLength = String.length y in
                if strLength > !longest_name then
                  longest_name := strLength;
                incr num_leaves

    let rec read_xml ch is_node nodeLst code info_table_add info_table_nonadd =
         try
                 let line = input_line ch in
                 if Str.string_match (Str.regexp " *<htu>") line 0 then
                         read_xml ch true nodeLst code info_table_add
                         info_table_nonadd
                 else if Str.string_match (Str.regexp " *<otu>") line 0 then
                         read_xml ch false nodeLst code info_table_add
                         info_table_nonadd
                 else if Str.string_match (Str.regexp " *</htu>") line 0 then
                         match nodeLst with
                         | hd :: Parser.Tree.Node (children, x):: tl  ->
                                 let newNode = Parser.Tree.Node (hd::children, x) in
                                    read_xml ch is_node (newNode :: tl) code
                                    info_table_add info_table_nonadd
                        | hd :: [] -> (hd, info_table_add, info_table_nonadd);
                        | _ -> raise (GraphTree.Wrong_format "End of Node but no Node in list")
                 else if Str.string_match 
                 (Str.regexp " *<code>") line 0 then
                     let code = GraphTree.get_code ch in
                     if is_node then
                         let newNode = Parser.Tree.Node ([], code) in 
                         read_xml ch is_node (newNode :: nodeLst) code
                         info_table_add info_table_nonadd
                     else
                         match nodeLst with
                         |Parser.Tree.Node (lst, x)::tl ->
                                 let newNode = Parser.Tree.Node 
                                 ((Parser.Tree.Leaf code)::lst, x) in
                                 read_xml ch is_node (newNode::tl) code
                                 info_table_add info_table_nonadd
                         | _ -> raise (GraphTree.Wrong_format "Adding a Leaf ")
                 else if Str.string_match (Str.regexp " *<preliminary>") line 0 then
                     let info_list_add, info_list_nonadd = 
                         gather_info (ref []) (ref []) ch in
                     (*print_endline ("for code: " ^ code);
                     List.iter print_info info_list; *)
                     Hashtbl.add info_table_add code info_list_add;
                     Hashtbl.add info_table_nonadd code info_list_nonadd;
                     read_xml ch is_node nodeLst code info_table_add 
                     info_table_nonadd;
                 else 
                     read_xml ch is_node nodeLst code info_table_add
                     info_table_nonadd;
         with
         |End_of_file ->
                 match nodeLst with
                 | hd :: _ -> (hd, info_table_add, info_table_nonadd)
                 | _ -> raise (GraphTree.Wrong_format "At end of file ")

    (** [draw_screen_from_xml ch] takes a channel to an xml file and draws a tree 
    *   also have optional arguments title, size and leafColor -
    *   title - takes a string 
    *   size - takes a string of format " 800x600" - notice must start with space
    *   leafColor - takes a Graphics color constant use as Graphics.red 
    *   choices are:
    *   black, red , blue, green, yellow, cyan, magenta *)    
    let draw_screen_from_xml ?(title="POY Tree Drawer") ?(size="") 
    ?(leafColor=Graphics.black) ch =
       let d = ref 0 and
       max_depth = ref 0 and
       num_leaves = ref 0 and
       longest_name = ref 0 and
       htbl = Hashtbl.create 30 and
       nonadd_htbl = Hashtbl.create 30 in
       let tree, info_table_add, info_table_nonadd  = read_xml ch true [] "" 
       htbl nonadd_htbl in
       xml_calc_depth_leaves tree d max_depth num_leaves longest_name;
       Graphics.open_graph size;
       Graphics.set_window_title title;
       let fontWidth, fontHeight = Graphics.text_size "A" in
       let sizeY = Graphics.size_y () and
       sizeX = Graphics.size_x () - (fontWidth * 2 * !longest_name) - 20 in
       let deltaY = sizeY / (!num_leaves + 1) and
       deltaX = sizeX  / !max_depth in
       let update_counter =
            let counter = ref 0 in
            fun () -> 
                incr counter;
                !counter
       in        
       let rec coord point_list deltaX deltaY tree =
           match tree with
           | Parser.Tree.Leaf name ->
                   let row = update_counter () in
                   let y = row * deltaY  and
                   x = (!max_depth * deltaX + 10) in
                   Graphics.moveto (x + fontWidth) (y-fontHeight/2);
                   let prev_color = Graphics.foreground in
                   Graphics.set_color leafColor;
                   Graphics.draw_string name;
                   Graphics.set_color prev_color;
                   point_list := (x, y, name)::!point_list;
                   (x, y)
           | Parser.Tree.Node (children, name) ->
                   let coord_children =
                       List.map (coord point_list deltaX deltaY) children
                   in
                   let xTemp , avg = GraphTreeScreen.average coord_children in
                   let x = !xTemp - deltaX in
                   List.iter (GraphTreeScreen.draw_edges (x, avg)) coord_children;
                   let prev_color = Graphics.foreground in
                   Graphics.set_color leafColor;
                   Graphics.plot x avg;
                   Graphics.moveto (x - (fontWidth * ((String.length name) + 1)))
                   (avg + fontHeight/2);
                   Graphics.draw_string name;
                   Graphics.set_color prev_color;
                   point_list := (x, avg, name) :: !point_list;
                   (x, avg)
       in
       let point_list = ref [] in
       let _ = coord point_list deltaX deltaY tree in
       (*List.iter print_points !point_list; *)
           let keep_window_open = ref true in
           while !keep_window_open do
               try
                   let _ = Graphics.wait_next_event
                   [Graphics.Key_pressed; Graphics.Button_down] in
                   if (Graphics.button_down ()) then
                       begin
                           let x, y = Graphics.mouse_pos () in
                           (*print_endline ("X is " ^ (string_of_int x));
                           print_endline ("Y is " ^ (string_of_int y)); *)
                           let code = find_closest_point x y point_list in
                           let code_info_list = Hashtbl.find info_table_add code in
                           print_endline ("\nInformation on code " ^ code);
                           List.iter print_info_verbose code_info_list;
                           List.iter print_info code_info_list;
                           let code_info_list = Hashtbl.find info_table_nonadd 
                           code in
                           List.iter print_info_nonadd_verbose code_info_list;
                           List.iter print_info_nonadd code_info_list;
                           print_newline ();
                       end 
                   else
                       keep_window_open := false;
                with
                   | Graphics.Graphic_failure _ -> 
                           keep_window_open := false;
                           Graphics.close_graph ()
          done
end

exception GetOutOfTheLoop 


let display all_trees = 
    let pos = ref 0 
    and max = (Array.length all_trees) - 1 in
    let do_print p =
        Graphics.clear_graph ();
        if p < 0 then pos := 0
        else if p > max then pos := max
        else pos := p;
        let cost , t = all_trees.(!pos) in
        draw_screen ~title:(string_of_int cost) t;
    in
    Graphics.open_graph "";
    do_print 0;
    try while true do
        try 
            (match Graphics.wait_next_event 
            [Graphics.Button_down; Graphics.Key_pressed] 
       with
        | { Graphics.keypressed = true; key = v } -> 
                if v = 'j' then 
                    if !pos = 0 then ()
                    else do_print (!pos - 1)
                else if v = 'k' then 
                    if !pos = max then ()
                    else do_print (!pos + 1)
                else if v = 'q' then 
                    raise GetOutOfTheLoop
                else ()
        | _ -> 
                Printf.printf "There is an issue here?";
                do_print !pos)
        with 
        | Graphics.Graphic_failure _ ->
                Graphics.open_graph "";
                do_print !pos
    done
    with
    | GetOutOfTheLoop -> Graphics.close_graph ()

