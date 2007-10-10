(***********************************************************************)
(*                                                                     *)
(*                           Objective Caml                            *)
(*                                                                     *)
(*            Pierre Weis, projet Cristal, INRIA Rocquencourt          *)
(*                                                                     *)
(*  Copyright 2000 Institut National de Recherche en Informatique et   *)
(*  en Automatique.  All rights reserved.  This file is distributed    *)
(*  under the terms of the GNU Library General Public License.         *)
(*                                                                     *)
(***********************************************************************)

(* $Id: graphps.ml 783 2006-05-15 20:53:49Z andres $ *)

open Printf;;

let graphps_version = "1.2";;

exception Graphic_failure of string;;

type command =
   | Open_graph of int * int
   | Close_graph

   | Set_rgb_color of int * int * int
   | Plot of int * int
   | Moveto of int * int
   | Lineto of int * int
   | Rmoveto of int * int
   | Rlineto of int * int

   | Curveto of (int * int) * (int * int) * (int * int)

   | Draw_rect of int * int * int * int
   | Fill_rect of int * int * int * int

   | Draw_poly of (int * int) array
   | Fill_poly of (int * int) array

   | Draw_arc of int * int * int * int * int * int
   | Fill_arc of int * int * int * int * int * int

   | Draw_ellipse of int * int * int * int
   | Fill_ellipse of int * int * int * int

   | Draw_circle of int * int * int
   | Fill_circle of int * int * int

   | Set_line_width of int

   | Draw_char of char
   | Draw_string of string
   | Set_font of string * int;;

let eps_mode = ref false;;
let ps_out_channel = ref stdout;;

let open_ps s = ps_out_channel := open_out s; eps_mode := false;;
let open_eps s = open_ps s; eps_mode := true;;

let current_ps_program = ref [];;

let window_opened = ref false;;

let opened () = !window_opened;;

let not_opened () = failwith "The graphic window is closed";;

let add_command c =
  if opened () then current_ps_program := c :: !current_ps_program else
  not_opened ();;

let size_x_ref = ref 640 and size_y_ref = ref 640;;

let size_x () = !size_x_ref;;
let size_y () = !size_y_ref;;

let parse_geometry s =
  let lim = String.length s in
  let rec find_x i =
   if i >= lim then raise Not_found else
   if s.[i] = 'x' then i else find_x (i + 1) in
  let rec find_digit i =
   if i >= lim then raise Not_found else
   match s.[i] with
   | '0' .. '9' -> i
   | _ -> find_digit (i + 1) in
  try
   let ix = find_x 0 in
   let dx = find_digit 0 in
   let sx = String.sub s dx (ix - dx) in
   let dy = find_digit ix in
   let sy = String.sub s dy (lim - dy) in
   int_of_string sx, int_of_string sy
  with
  | Not_found | Failure _ -> (640, 640);;

let open_graph s =
 let x, y = parse_geometry s in
 size_x_ref := x; size_y_ref := y;
 window_opened := true;
 add_command (Open_graph (x, y));;

let clear_ps_program p =
  let rec clear accu = function
  | [] -> accu
  | c :: cs -> if is_drawing c then clear accu cs else clear (c :: accu) cs
  and is_drawing = function
  | Set_font (_, _) |
    Set_line_width _ |
    Set_rgb_color (_, _, _) -> false
  | Open_graph (_, _) | Close_graph |
    Plot (_, _) | Moveto (_, _) | Lineto (_, _) |
    Rmoveto (_, _) | Rlineto (_, _) |
    Curveto (_, _, _) |
    Draw_rect (_, _, _, _) | Fill_rect (_, _, _, _) |
    Draw_poly (_) | Fill_poly (_) |
    Draw_arc (_, _, _, _, _, _) | Fill_arc (_, _, _, _, _, _) |
    Draw_ellipse (_, _, _, _) | Fill_ellipse (_, _, _, _) |
    Draw_circle (_, _, _) | Fill_circle (_, _, _) |
    Draw_char _ | Draw_string _ -> true in
  clear [Open_graph (size_x (), size_y ())] p;;

let clear_graph () =
  current_ps_program := clear_ps_program (List.rev !current_ps_program);;

(* Colors *)

type color = int;;

let rgb r g b = (r lsl 16) + (g lsl 8) + b;;

let rgb_of_color c =
 let r = (c lsr 16) land 0xFF in
 let g = (c lsr 8) land 0xFF in
 let b = c land 0xFF in
 r, g, b;;

let black   = 0x000000
and white   = 0xFFFFFF
and red     = 0xFF0000
and green   = 0x00FF00
and blue    = 0x0000FF
and yellow  = 0xFFFF00
and cyan    = 0x00FFFF
and magenta = 0xFF00FF;;

let background = white
and foreground = black;;

let current_rgb = ref min_int;;

let set_color c =
 let c = if c < 0 then background else c in
 if c <> !current_rgb then begin
   let r, g, b = rgb_of_color c in
   current_rgb := c;
   add_command (Set_rgb_color (r, g, b)) end;;

(* Drawing *)

let plot x y = add_command (Plot (x, y));;

let plots points =
  for i = 0 to Array.length points - 1 do
    let (x, y) = points.(i) in
    plot x y;
  done
;;

let curr_x = ref 0;;
let curr_y = ref 0;;

let current_x () = !curr_x;;
let current_y () = !curr_y;;

let current_point () = !curr_x, !curr_y;;

let set_point x y = curr_x := x; curr_y := y;;

let moveto x y =
 set_point x y;
 add_command (Moveto (x, y));;

let lineto x y =
 set_point x y;
 add_command (Lineto (x, y));;

let rlineto dx dy = 
 set_point (current_x () + dx) (current_y () + dy);
 add_command (Rlineto (dx, dy));;

let rmoveto dx dy =
 set_point (current_x () + dx) (current_y () + dy);
 add_command (Rmoveto (dx, dy));;

let curveto (x1, y1 as b) (x2, y2 as c) (x3, y3 as d) =
 add_command (Curveto (b, c, d));
 set_point x3 y3;;

let draw_arc x y rx ry a1 a2 = add_command (Draw_arc (x, y, rx, ry, a1, a2));;
let draw_ellipse x y rx ry = add_command (Draw_ellipse (x, y, rx, ry));;
let draw_circle x y r = add_command (Draw_circle (x, y, r));;

let set_line_width w = add_command (Set_line_width w);;

let draw_rect x y w h = add_command (Draw_rect (x, y, w, h));;
let fill_rect x y w h = add_command (Fill_rect (x, y, w, h));;

let fill_poly v = add_command (Fill_poly v);;
let draw_poly v = add_command (Draw_poly v);;
let draw_poly_line =
  let draw points =
    if Array.length points > 0 then begin
      let (savex, savey) = current_point () in
      moveto (fst points.(0)) (snd points.(0));
      for i = 1 to Array.length points - 1 do
        let (x, y) = points.(i) in
        lineto x y;
      done;
      moveto savex savey;
    end in
  draw;;

let draw_segments segs =
  let (savex, savey) = current_point () in
  for i = 0 to Array.length segs - 1 do
    let (x1, y1, x2, y2) = segs.(i) in
    moveto x1 y1;
    lineto x2 y2;
  done;
  moveto savex savey;;

let fill_arc x y rx ry a1 a2 = add_command (Fill_arc (x, y, rx, ry, a1, a2));;

let fill_ellipse x y rx ry = add_command (Fill_ellipse (x, y, rx, ry));;
let fill_circle x y r = add_command (Fill_circle (x, y, r));;

(* Text *)
let default_font = "Helvetica-Bold";;
let default_font_size = 10;;

let x_text_size = ref default_font_size;;
let y_text_size = ref default_font_size;;

let draw_char c =
 let x = current_x ()
 and y = current_y () in
 set_point (x + !x_text_size) y;
 add_command (Draw_char c);;

let draw_string s =
 let x = current_x ()
 and y = current_y () in
 set_point (x + (String.length s * !x_text_size)) y;
 add_command (Draw_string s);;

(* The font machinery *)
let current_font = ref "";;

let fonts = Hashtbl.create 11;;

let add_font f =
 if not (Hashtbl.mem fonts f) then Hashtbl.add fonts f f;;

let parse_font f =
  match f with
  | "Times" | "Times-Bold" | "Times-Italic"
  | "Times-BoldItalic" | "Times-Roman"
  | "Courier" | "Courier-Bold"
  | "Courier-Italic" | "Courier-BoldItalic"
  | "Helvetica" | "Helvetica-Bold"
  | "Helvetica-Italic" | "Helvetica-BoldItalic" -> f
  | _ ->
      prerr_endline
        (Printf.sprintf
           "Warning: cannot set font to %s.\n\
            Choose Times, Courier, or Helvetica (see the documentation)"
            f);
      failwith "set_font";;

let change_font f =
 let f = parse_font f in
 add_font f;
 current_font := f;;

let set_font f =
  if f <> !current_font then begin
    change_font f;
    let sz = !x_text_size in
    add_command (Set_font (f, sz))
  end;;

let change_text_size i =
  x_text_size := i;
  y_text_size := i;;

let set_text_size i =
 if i <> !x_text_size then begin
   change_text_size i;
   add_command (Set_font (!current_font, i))
 end;;

(* Init fonts for this module: PostScript font initialization is done
   via the prelude. *)
let init_font () =
 change_text_size default_font_size;
 change_font default_font;;

init_font ();;

(* Computing text sizes is highly non trivial:
   strictly speaking we should launch a PostScript interpreter to compute it!
   More or less, we define text_size s as
   ((String.length s * !x_text_size + 1) / 2,
   (!y_text_size + 1) / 2);;
   In the following we perform a desesperate attempt to do a better job ...
*)
let text_size s =
  let size_x_char = function
    | 'i' | 'j' | 'l' | 't' | '\'' | ' ' -> max 1 ((3 * !x_text_size + 2) / 4)
    | 'A' .. 'Z' | 'm' | 'n' | 'g' | 'd' -> (5 * !x_text_size + 2) / 4
    | _ -> !x_text_size in
  let l = String.length s in
  let rec loop accu i =
   if i >= l then accu else loop (accu + size_x_char s.[i]) (i + 1) in
  (loop 1 0) / 2, (!y_text_size + 1) / 2;;

(* Images *)

type image;;

let transp = -1;;

(* Post Script functions embedded into the output *)
let ps_defs = "\
/m { moveto } def
/rm { rmoveto } def
/l { lineto } def
/c { currentpoint } def
/cp { currentlinewidth currentpoint } def
/ms { c stroke m } def
/dl { l ms } def
/rl { rlineto ms } def
/p { newpath m c lineto 1 setlinewidth stroke m setlinewidth } def
/pr  %% define the path for a rectangle w h x y
    { newpath m
      dup %% w h h
      0 exch rlineto %% rlineto 0 h -- w h
      exch 0 rlineto %% rlineto w 0 -- h
      0 exch sub 0 exch rlineto %% rlineto 0 -h --
      closepath } def
/dr %% w h x y
    { pr stroke m } def
/fr { pr fill m } def
/fc { 0 360 newpath arc fill m } def
/dc { 0 360 newpath arc stroke m } def
/fa { gsave translate scale newpath 0 0 m arc closepath fill grestore } def
/da { savematrix gsave translate scale newpath arc
      restorematrix stroke grestore } def
/fe { gsave translate scale newpath 0 0 1 0 360 arc fill grestore } def
/de { savematrix gsave translate scale newpath 0 0 1 0 360 arc
      restorematrix stroke grestore } def
/pc { newpath m } def
/dp { closepath stroke m } def
/fp { closepath fill m } def
/ct { curveto c stroke m } def
/kcolor { 1 255 div } def
/color { kcolor mul } def
/stmp { /tmp exch def } def
/srgb {% r g b -
 color stmp color exch color exch tmp setrgbcolor } def
/t { show } def
/savedmatrix [0 0 0 0 0 0] def
/savematrix { savedmatrix currentmatrix pop } def
/restorematrix { savedmatrix setmatrix } def
%% ISO fonts
/ISOLatin1Encoding [
 /.notdef /.notdef /.notdef /.notdef /.notdef /.notdef /.notdef /.notdef
 /.notdef /.notdef /.notdef /.notdef /.notdef /.notdef /.notdef /.notdef
 /.notdef /.notdef /.notdef /.notdef /.notdef /.notdef /.notdef /.notdef
 /.notdef /.notdef /.notdef /.notdef /.notdef /.notdef /.notdef /.notdef
 /space /exclam /quotedbl /numbersign /dollar /percent /ampersand /quoteright
 /parenleft /parenright /asterisk /plus /comma /minus /period /slash
 /zero /one /two /three /four /five /six /seven
 /eight /nine /colon /semicolon /less /equal /greater /question
 /at /A /B /C /D /E /F /G /H /I /J /K /L /M /N /O
 /P /Q /R /S /T /U /V /W /X /Y /Z /bracketleft /backslash /bracketright
                                                       /asciicircum /underscore
 /quoteleft /a /b /c /d /e /f /g /h /i /j /k /l /m /n /o
 /p /q /r /s /t /u /v /w /x /y /z /braceleft /bar /braceright /asciitilde
                                                                       /.notdef
 /.notdef /.notdef /.notdef /.notdef /.notdef /.notdef /.notdef /.notdef
 /.notdef /.notdef /.notdef /.notdef /.notdef /.notdef /.notdef /.notdef
 /.notdef /.notdef /.notdef /.notdef /.notdef /.notdef /.notdef /.notdef
 /.notdef /.notdef /.notdef /.notdef /.notdef /.notdef /.notdef /.notdef
 /space /exclamdown /cent /sterling /currency /yen /brokenbar /section
 /dieresis /copyright /ordfeminine /guillemotleft /logicalnot /hyphen
                                                            /registered /macron
 /degree /plusminus /twosuperior /threesuperior /acute /mu /paragraph
                                                                /periodcentered
 /cedilla /onesuperior /ordmasculine /guillemotright /onequarter /onehalf
                                                   /threequarters /questiondown
 /Agrave /Aacute /Acircumflex /Atilde /Adieresis /Aring /AE /Ccedilla
 /Egrave /Eacute /Ecircumflex /Edieresis /Igrave /Iacute /Icircumflex
                                                         /Idieresis
 /Eth /Ntilde /Ograve /Oacute /Ocircumflex /Otilde /Odieresis /multiply
 /Oslash /Ugrave /Uacute /Ucircumflex /Udieresis /Yacute /Thorn /germandbls
 /agrave /aacute /acircumflex /atilde /adieresis /aring /ae /ccedilla
 /egrave /eacute /ecircumflex /edieresis /igrave /iacute /icircumflex
                                                         /idieresis
 /eth /ntilde /ograve /oacute /ocircumflex /otilde /odieresis /divide
 /oslash /ugrave /uacute /ucircumflex /udieresis /yacute /thorn /ydieresis
] def
%% usage: isoname oldname makeisofont -
/makeisofont {
  dup findfont length dict dup begin
    exch findfont {
      exch dup /FID eq { pop pop } { exch def } ifelse
    } forall
    /Encoding ISOLatin1Encoding def
  end
  definefont pop
} def
";;

let make_iso_font f = sprintf
 "/ISO%s /%s makeisofont" f f;;

let hashtbl_assocs t =
 let res = ref [] in
 Hashtbl.iter (fun k v -> res := (k, v) :: !res) t;
 !res;;

let hashtbl_vals t =
 let res = ref [] in
 Hashtbl.iter (fun k v -> res := v :: !res) t;
 !res;;

let hashtbl_map f t =
 let res = ref [] in
 Hashtbl.iter (fun k v -> res := f k v :: !res) t;
 !res;;

let make_iso_fonts () =
 let font_list = hashtbl_vals fonts in
 String.concat "\n" (List.map make_iso_font font_list);;

let set_default_font default_font default_font_size = sprintf
 "/ISO%s findfont %i scalefont setfont\n" default_font default_font_size;;

let creator size_x size_y = sprintf
  "%%!PS-Adobe-2.0 EPSF-2.0\n\
   %%%%Creator: GraphPs %s\n\
   %%%%BoundingBox: 0 0 %i %i\n" graphps_version size_x size_y;;

let document_fonts () =
  let font_list = hashtbl_vals fonts in
  sprintf
   "%%%%DocumentFonts: %s\n\
    %%%%EndComments\n" (String.concat " " font_list);;

let header size_x size_y =
  sprintf "%s%s" (creator size_x size_y) (document_fonts ());;

let prelude size_x size_y default_font default_font_size =
  sprintf
   "%s\
    %%BeginProcSet\n\
    gsave\n\
    %s\
    0 0 m\n\
    1 setlinewidth\n\
    1 setlinecap\n\
    2 setlinejoin\n\
    10 setmiterlimit\n\
    %s\n\
    %s\
    %%EndProcSet\n" (header size_x size_y)
        ps_defs
        (make_iso_fonts ())
        (set_default_font default_font default_font_size);;

(* Could be "showpage", if one wants automatic display of the current page. *)
let postlude () =
 if !eps_mode then "grestore\n" else "grestore\nshowpage\n%%end\n";;

let escape_char b = function
  | '(' | ')' | '\\' as c -> Buffer.add_char b '\\'; Buffer.add_char b c
  | c -> Buffer.add_char b c;;

let escape_char_for_ps c =
 let b = Buffer.create 3 in
 escape_char b c;
 Buffer.contents b;;

let escape_string_for_ps s =
 let l = String.length s in
 let b = Buffer.create l in
 for i = 0 to l - 1 do escape_char b s.[i] done;
 Buffer.contents b;;

type filling = Fill | Draw;;

let print_poly oc filling v =
  if Array.length v > 0 then begin
   let x, y = v.(0) in
   fprintf oc "c %i %i pc\n" x y;
   for i = 1 to Array.length v - 1 do
    let x, y = v.(i) in
    fprintf oc "%i %i l\n" x y
   done;
   if filling = Draw then fprintf oc "dp\n" else fprintf oc "fp\n"
  end;;

let to_ps oc = function
  | Open_graph (x, y) ->
     fprintf oc "%s" (prelude x y default_font default_font_size);
     set_color background;
     fill_rect 0 0 x y;
     set_color foreground
  | Close_graph ->
     fprintf oc "%s\n" (postlude ())
  | Set_rgb_color (r, g, b) ->
     fprintf oc "%i %i %i srgb\n" r g b
  | Plot (x, y) ->
     fprintf oc "cp %i %i p\n" x y
  | Moveto (x, y) ->
     fprintf oc "%i %i m\n" x y
  | Lineto (x, y) ->
     fprintf oc "%i %i dl\n" x y
  | Rmoveto (x, y) ->
     fprintf oc "%i %i rm\n" x y
  | Rlineto (x, y) ->
     fprintf oc "%i %i rl\n" x y
  | Curveto ((x1, y1), (x2, y2), (x3, y3)) ->
     fprintf oc "%i %i %i %i %i %i ct\n" x1 y1 x2 y2 x3 y3
  | Set_line_width w ->
     fprintf oc "%i setlinewidth\n" w
  | Draw_char c ->
     fprintf oc "(%s) show\n" (escape_char_for_ps c)
  | Draw_string s ->
     fprintf oc "(%s) show\n" (escape_string_for_ps s)
  | Set_font (f, size) ->
     fprintf oc "/ISO%s findfont %i scalefont setfont\n" f size
  | Draw_rect (x, y, w, h) ->
     fprintf oc "c %i %i %i %i dr\n" w h x y
  | Fill_rect (x, y, w, h) ->
     fprintf oc "c %i %i %i %i fr\n" w h x y
  | Fill_poly v ->
     print_poly oc Fill v
  | Draw_poly v ->
     print_poly oc Draw v
  | Draw_circle (x, y, r) ->
     fprintf oc "c %i %i %i dc\n" x y r
  | Fill_circle (x, y, r) ->
     fprintf oc "c %i %i %i fc\n" x y r
  | Draw_ellipse (x, y, rx, ry) ->
     fprintf oc "%i %i %i %i de\n" rx ry x y
  | Fill_ellipse (x, y, rx, ry) ->
     fprintf oc "%i %i %i %i fe\n" rx ry x y
  | Draw_arc (x, y, rx, ry, a1, a2) ->
     fprintf oc "0 0 1 %i %i %i %i %i %i da\n" a1 a2 rx ry x y
  | Fill_arc (x, y, rx, ry, a1, a2) ->
     fprintf oc "0 0 1 %i %i %i %i %i %i fa\n" a1 a2 rx ry x y
;;

let rev_opt l =
  let same_com com c =
    match com, c with
    | Set_rgb_color (r, g, b), Set_rgb_color (_, _, _) -> true
    | Set_line_width w, Set_line_width _ -> true
    | Set_font (f, size), Set_font (_, _) -> true
    | _, _ -> false in

  let rec opt_rec com accu l =
    match l with
    | [] -> com :: accu
    | c :: l ->
    match c with
    | Set_rgb_color (r, g, b) as c ->
       if same_com com c then opt_rec com accu l else opt_rec c (com :: accu) l
    | Moveto (x, y) ->
       begin match com with
       | Rmoveto (_, _) | Moveto (_, _) -> opt_rec com accu l
       | _ -> opt_rec c (com :: accu) l
       end
    | Rmoveto (x, y) ->
       begin match com with
       | Rmoveto (x0, y0) -> opt_rec (Rmoveto (x + x0, y + y0)) accu l
       | Moveto (x0, y0) -> opt_rec com accu l
       | _ -> opt_rec c (com :: accu) l
       end
    | Set_line_width w ->
       if same_com com c then opt_rec com accu l else opt_rec c (com :: accu) l
    | Set_font (f, size) ->
       if same_com com c then opt_rec com accu l else opt_rec c (com :: accu) l
    | Open_graph (_, _) | Close_graph
    | Plot (_, _)
    | Lineto (_, _) | Rlineto (_, _) | Curveto (_, _, _)
    | Draw_char _ | Draw_string _
    | Draw_rect (_, _, _, _)
    | Fill_rect (_, _, _, _) | Fill_poly _ | Draw_poly _
    | Draw_circle (_, _, _) | Fill_circle (_, _, _)
    | Draw_ellipse (_, _, _, _) | Fill_ellipse (_, _, _, _)
    | Draw_arc (_, _, _, _, _, _) | Fill_arc (_, _, _, _, _, _) ->
       opt_rec c (com :: accu) l in
  match l with
  | [] -> l
  | c :: l -> opt_rec c [] l;;

let close_graph () =
  if opened () then begin
    let oc = !ps_out_channel in
    add_command Close_graph;
    let commands = !current_ps_program in
    List.iter (to_ps oc) (rev_opt commands);
    current_ps_program := [];
    flush oc;
    if oc != stdout then close_out oc;
    window_opened := false
  end else not_opened ();;

at_exit 
    (fun () -> 
        try close_graph () with
        | Failure _ -> ());;
