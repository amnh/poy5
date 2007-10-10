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

let () = SadmanOutput.register "NcursesML" "$Revision: 1644 $"

(** OCaml interface to NCURSES library *)

(** On Not Working when Run Through MPIEXEC
 *  a Short Narrative by Illya Bomash
 * 
 *  Once upon a midnight dreary, MPICH decided to change its mpirun script to a
 *  compiled program called mpiexec.  Unfortunately, this compiled program does
 *  strange things that render our trivial implementation of NCURSES incorrect.
 * 
 *  Symptoms:  number of rows and columns is wrong; things get printed in a very
 *  strange manner.
 * 
 *  It looks like all environmental variables (including $TERM) are preserved
 *  when running through mpiexec (cf: "mpiexec env").  It also looks like the
 *  terminal settings aren't changing, and, in particular, they report the
 *  correct number of rows and columns (see "mpiexec stty -a").
 * 
 *  There's got to be a specific way in which NCURSES reads its information ...
 *  that should be the biggest clue about how to fix this. 
 * 
 *  Other programs displaying this problem:  aptitude *)


type window

external init : unit -> window = "ncurs_CAML_init"

external finalize : unit -> unit = "ncurs_CAML_endwin"

let stdscr = init ()

let redraw = ref (fun () -> ())

let do_redraw () = !redraw ()

let set_redraw_function f = redraw := f

external redrawwin : window -> unit = "ncurs_CAML_redrawwin"

external window_create : window -> int -> int -> int -> int -> window = 
    "ncurs_CAML_create"

external getyx : window -> (int * int) = "ncurs_CAML_getyx"

external getparyx : window -> (int * int) = "ncurs_CAML_getparyx"

external getmaxyx : window -> (int * int) = "ncurs_CAML_getmaxyx"

external getbegyx : window -> (int * int) = "ncurs_CAML_getbegyx"

external box : window -> int -> int -> unit = "ncurs_CAML_box"

external c_wrefresh : window -> unit = "ncurs_CAML_wrefresh"

external scroll : window -> int -> unit = "ncurs_CAML_scroll"

external c_refresh : unit -> unit = "ncurs_CAML_refresh"

let refresh ?w () =
    match w with
    | Some w -> c_wrefresh w
    | None -> c_refresh ()

external echo : bool -> unit = "ncurs_CAML_echo"

external raw : unit -> unit = "ncurs_CAML_raw"

external noraw : unit -> unit = "ncurs_CAML_noraw"

external cbreak : unit -> unit = "ncurs_CAML_cbreak"

external c_keypad : window -> bool -> unit = "ncurs_CAML_keypad"

external do_update : unit -> unit = "ncurs_CAML_doupdate"

external wnoutrefresh : window -> unit = "ncurs_CAML_wnoutrefresh"

let keypad ?(w=stdscr) s =
    c_keypad w s

external meta : window -> bool -> unit = "ncurs_CAML_meta"

external addstr : string -> unit = "ncurs_CAML_addstr"

external waddstr : window -> string -> unit = "ncurs_CAML_waddstr"

external waddch : window -> int ->  unit = "ncurs_CAML_waddch"

let printw (v : ('a, unit, string) format) =
    let str = Printf.sprintf v in
    addstr str;
    str

let wprintw w v = 
    let str = Printf.sprintf v in
    waddstr w str;
    str

external mvaddstr : int -> int -> string -> unit = "ncurs_CAML_mvaddstr"

let mvprintw y x v =
    let str = Printf.sprintf v in
    mvaddstr y x str;
    str

external mvwaddstr : window -> int -> int -> string -> unit =
    "ncurs_CAML_mvwaddstr"

external scrollok : window -> bool -> unit = "ncurs_CAML_scrollok"

external idlok : window -> bool -> unit = "ncurs_CAML_idlok"

let mvwprintw w y x v =
    let str = Printf.sprintf v in
    mvwaddstr w y x str;
    str

external getch : unit -> int = "ncurs_CAML_getch"

external wgetch : window -> int = "ncurs_CAML_wgetch"

type attribut =
    | STANDOUT 
    | UNDERLINE
    | REVERSE
    | BLINK
    | BOLD
    | PROTECT
    | INVIS
    | ALTCHARSET

external attr_standout : window -> bool -> unit = "ncurs_CAML_A_STANDOUT"

external attr_underline : window -> bool -> unit = "ncurs_CAML_A_UNDERLINE"

external attr_reverse : window -> bool -> unit = "ncurs_CAML_A_REVERSE"

external attr_blink : window -> bool -> unit = "ncurs_CAML_A_BLINK"

external attr_bold : window -> bool -> unit = "ncurs_CAML_A_BOLD"

external attr_protect : window -> bool -> unit = "ncurs_CAML_A_PROTECT"

external attr_invis : window -> bool -> unit = "ncurs_CAML_A_INVIS"

external attr_altcharset : window -> bool -> unit = "ncurs_CAML_A_ALTCHARSET"

let attr_single_set clas =
    match clas with
    | STANDOUT -> attr_standout 
    | UNDERLINE -> attr_underline
    | REVERSE -> attr_reverse
    | BLINK -> attr_blink
    | BOLD -> attr_bold
    | PROTECT -> attr_protect
    | INVIS -> attr_invis
    | ALTCHARSET -> attr_altcharset

let attr_set_unset bol ?(window=stdscr) params =
    List.iter (fun x -> attr_single_set x window bol) params

let attr_set = attr_set_unset true

let attr_unset = attr_set_unset false

external columns : unit -> int = "ncurs_CAML_columns"

external lines : unit -> int = "ncurs_CAML_lines"

external is_term_resized : int -> int -> bool = "ncurs_CAML_istermresized"

external wresize : window -> int -> int -> unit = "ncurs_CAML_wresize"

external mvwin : window -> int -> int -> unit = "ncurs_CAML_mvwin"

external wmove : window -> int -> int -> unit = "ncurs_CAML_wmove"

external delete: window -> unit = "ncurs_CAML_delete"

external wdelch: window -> unit = "ncurs_CAML_wdelch"

let unbox window =
    box window (Char.code ' ') (Char.code ' ')

external nocbreak : unit -> unit = "ncurs_CAML_nocbreak"

external mvwhline : window -> int -> int -> int -> int -> unit =
    "ncurs_CAML_mvwhline"

external has_colors : unit -> bool = "ncurs_CAML_has_colors"

external start_color : unit -> unit = "ncurs_CAML_start_color"

external wdeleteln : window -> unit = "ncurs_CAML_wdeleteln"

external isprint : int -> bool = "ncurs_CAML_isprint"

external keyup : unit -> int = "ncurs_CAML_keyup"
external keydown : unit -> int = "ncurs_CAML_keydown"
external keyleft : unit -> int = "ncurs_CAML_keyleft"
external keyright : unit -> int = "ncurs_CAML_keyright"
external keybackspace : unit -> int = "ncurs_CAML_keybackspace"
external keydc : unit -> int = "ncurs_CAML_keydc"
external keynpage : unit -> int = "ncurs_CAML_keynpage"
external keyppage : unit -> int = "ncurs_CAML_keyppage"
external resize : unit -> int = "ncurs_CAML_resize"
external error : unit -> int = "ncurs_CAML_error"


let next_color_pair = ref 1
type color =
    | Red
    | Green
    | Yellow
    | Blue
    | Cyan
    | Magenta
    | White
    | Black
external init_pair : int -> color -> color -> unit = "ncurs_CAML_init_pair"
module CompareColorPair = struct
    type t = color * color
    let compare (a, b) (c, d) =
        let r = Pervasives.compare a c in
        if r = 0 then Pervasives.compare b d
        else r
end
module ColorMap = Map.Make (CompareColorPair)
let color_pairs = ref ColorMap.empty

let get_color fg bg =
    try ColorMap.find (fg, bg) !color_pairs
    with Not_found ->
        let code = !next_color_pair in
        incr next_color_pair;
        init_pair code fg bg;
        color_pairs := ColorMap.add (fg, bg) code !color_pairs;
        code

external use_color : int -> window -> bool -> unit = "ncurs_CAML_use_color"

let red_in_black () = ()

let use_red w = use_color (get_color Red Black) w
let use_green w = use_color (get_color Green Black) w
let use_blue w = use_color (get_color Blue Black) w
let use_cyan w = use_color (get_color Cyan Black) w
let use_yellow w = use_color (get_color Yellow Black) w
let use_magenta w = use_color (get_color Magenta Black) w
let use_white w = use_color (get_color White Black) w
let use_red_in_whit w = use_color (get_color Red White) w
let use_green_in_whit w = use_color (get_color Green White) w
let use_blue_in_whit w = use_color (get_color Blue White) w
let use_cyan_in_whit w = use_color (get_color Cyan White) w
let use_yellow_in_whit w = use_color (get_color Yellow White) w
let use_magenta_in_whit w = use_color (get_color Magenta White) w
let use_black_in_whit w = use_color (get_color Black White) w  

external winsch : window -> int -> unit = "ncurs_CAML_winsch"
external halfdelay : int -> unit = "ncurs_CAML_halfdelay"
external sigwinch_supported : unit -> bool = "ncurs_CAML_sigwinch_supported"
external sigwinch : unit -> int = "ncurs_CAML_sigwinch"
external sigwinchhandler : int -> unit = "ncurs_CAML_sigwinch_handler"
external resize_handler : unit -> unit = "ncurs_CAML_resize_handler"

external wclear : window -> unit = "ncurs_CAML_wclear"
external touchwin : window -> unit = "ncurs_CAML_touchwin"
external clearok : window -> bool -> unit = "ncurs_CAML_clearok"
external werase :  window -> unit = "ncurs_CAML_werase"
