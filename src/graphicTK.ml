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

let () = SadmanOutput.register "GraphicTK" "$Revision: 1616 $"

module Graphictktree : GraphTree.GRAPHICS_TYPE = struct 

    type color=int 
    type my_color = Tk.color

    let  color_conversion x = 
        match  x with 
        | 0 -> `Black
        | 1 -> `Red
        | 2 -> `Blue
        | 3 -> `Green
        | 4 -> `Yellow              
        | _ -> failwith " Unsupported color" 

    (* foreground : color
    Default background and foreground colors (usually, either black foreground 
    on a white background or white foreground on a black background).
    Graphics.clear_graph fills the screen with the background color. The initial
    drawing color is foreground.*)


    let foreground = 0    
    let black = 0
    let red = 1
    let close_graph = Tk.closeTk

    (* set_color : color -> unit
    Set the current drawing color*)

    let current_color = ref `Black
    let set_color (color : int) = current_color := (color_conversion color)



    let can = ref None;;               

    (** [open_graph tex]
    Show the graphics window or switch the screen to graphic mode. The
    graphics window is cleared and the current point is set to (0, 0).
    The string argument is used to pass optional information on the
    desired graphics mode, the graphics window size, and so on. Its
    interpretation is implementation-dependent. If the empty string is
    given, a sensible default is selected.*)


    let open_graph _ =
        let top = Tk.openTk ()  in
        let max_x = 1000 (*info.screenwidth  top*) in 
        let max_y = 1000 (*Winfo.screenheight top*) in          
        let base = 
            Frame.create ~width:max_x ~height:max_y
            ~takefocus:true ~visual:`Best  
        top in Tk.pack [base]; 
        let tmp = 
            Canvas.create ~width: max_x ~height:max_y
            ~borderwidth:2  base
        in               
        Tk.pack [tmp];
        can := Some tmp


    let size_x () = 1000 

    let size_y () = 1000


    let current_point = ref (0, 0)         

    let moveto x y = current_point := (x, y)  

    (** lineto : int -> int -> unit
    Draw a line with endpoints the current point and the given point,
    and move the current point to the given point.*)


    let lineto x y  = 
        match !can with
        | None -> failwith "The canvas is None"
        | Some can ->
                let _ = 
                    Canvas.create_line ~xys:[!current_point; (x, y)]
                    ~fill:!current_color ~width:1 can  
                in 
                moveto x y

    let text_size str = 
        match !can with
        | None -> failwith "The canvas is None"
        | Some can ->              
                let text = Canvas.create_text ~x:0 ~y:0 ~text:str can in
                let (a, b, c, d) = Canvas.bbox can [text]  in
                Canvas.delete can [text];
                c - a, d - b



             (** [draw_string str]
             Draw a character or a character string with lower left corner at current position. After 
             drawing, the current position is set to the lower right corner of the text drawn. *) 
             let draw_string str = 
                   match !can with
                    | None -> failwith "The canvas is None"
                    | Some can ->
                        let w, h = text_size str in 
                        let x, y = !current_point in 
                        let _ = 
                            Canvas.create_text
                            ~x:(x + w/2)  ~y:(y + h/2)
                            ~text: str 
                            ~fill:!current_color can 
                        in        
                        ()  
 
                    
            (* val plot : int -> int -> unit
            Plot the given point with the current drawing color. *)
             let plot x y =
                 match !can with
                 | None -> failwith "The canvas is None"
                 | Some can ->
                         let _ = 
                             Canvas.create_oval ~x1: x ~y1: y ~x2:x ~y2:y                         
                             ~fill:!current_color can 
                         in 
                         ()            


             let display () =    
                 let _= Tk.mainLoop () in 
                 ()             

end


module GraphicsTK = GraphTree.Make (Graphictktree)


let display all_trees = 
    let iterator tree = 
        Graphictktree.open_graph "";
        let _, tree = tree in  
        GraphicsTK.draw tree
    in
    Array.iter iterator all_trees; 
    Graphictktree.display ();
    Graphictktree.close_graph ()          












