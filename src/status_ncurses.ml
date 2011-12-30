(* POY 5.0 Alpha. A phylogenetic analysis program using Dynamic Homologies.   *)
(* Copyright (C) 2011 Andrés Varón, Lin Hong, Nicholas Lucaroni, Ward Wheeler,*)
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

let () = SadmanOutput.register "" "$Revision: 2592 $"

let () = SadmanOutput.register "Status_ncurses" "$Revision: 2592 $"

type tab_state = Begin | First | Continue

(** [ndebug_keycode]: if true, then the values of unknown keycodes are printed
    to the input window. *)
let ndebug_keycode = true

let get_next_code =
    let cnt = ref (-1) in
    fun () -> incr cnt; !cnt

type status = {
    line : int ref;
    formatter : StatusCommon.Format.formatter;
    window : NcursesML.window;
    code : int;
    max : int option;
    name : string;
    achieved : int ref;
    last_message : string ref;
    cleanup : (unit -> unit);
}

exception Illegal_update


type formatter_output = StatusCommon.formatter_output

type c = 
    | SearchReport 
    | Status 
    | Warning
    | Error 
    | Information 
    | Output of (string option * bool * formatter_output list)  (* Filename and close it or not *)

let are_we_parallel = ref false

let slaves_deal_in_this_way : (c -> string -> unit) ref = ref (fun _ _ -> ())

(* If running in parallel, wht is my rank? *)
let my_rank = ref (-1)

let is_parallel x f = 
    are_we_parallel := true;
    my_rank := x;
    match f with
    | None -> ()
    | Some f -> 
            slaves_deal_in_this_way := f

(* [rank c] sets the stored rank of the calling process to [c] *)
let rank c = my_rank := c

(* The verbosity that will be used in a particular run *)
let verbosity : ([ `None | `Low | `Medium | `High ] ref) = ref `Low

(* [set_verbosity v] sets the verbosity of the calling process to [v] *)
let set_verbosity v = verbosity := v

(* [get_verbosity ()] gets the  current setting for the verbosity level *)
let get_verbosity () = !verbosity

let columns =  
    ref (NcursesML.columns ())

let set_columns v = columns := max 75 v

let lines = 
    ref (NcursesML.lines ())

let set_lines v = lines := max 20 v

let parstatus_height () = 7
let parconsole_height () = 10
let parconsole_width () = (!columns / 2) - 1
let parinformation_height () =
    !lines - (1 + parstatus_height ()) - (1 + parconsole_height ())

let parinformation_width () = !columns

let parinformation_pos1 () = 0
let parinformation_pos2 () = 0

let make_parinformation () =
    (NcursesML.window_create NcursesML.stdscr (parinformation_height ())
    (parinformation_width ()) (parinformation_pos1 ()) (parinformation_pos2 ()))

let parinformation = ref (make_parinformation ())

let make_parconsole () =
    (NcursesML.window_create NcursesML.stdscr (parconsole_height ())
    (parconsole_width ()) (parinformation_height ()) (parinformation_pos2 ()))

let parconsole = ref (make_parconsole ())

let parsearchstatus_width () = !columns - (parconsole_width ()) - 1

let make_parsearchstatus () = 
    (NcursesML.window_create NcursesML.stdscr (parconsole_height ())
    (parsearchstatus_width ()) (parinformation_height ()) 
    ((parconsole_width ()) + 1))

let parsearchstatus = ref (make_parsearchstatus ())

let parstatus_pos1 () = (parinformation_height ()) + (parconsole_height ()) + 1

let make_parstatus () =
    (NcursesML.window_create NcursesML.stdscr (parstatus_height ())
    !columns (parstatus_pos1 ()) 0)

let parstatus = ref (make_parstatus ())

let make_information () =
    (NcursesML.window_create !parinformation ((parinformation_height ()) - 2) 
    (!columns - 2) 1 1)

let information = ref (make_information ())

let make_console () = 
    (NcursesML.window_create !parconsole ((parconsole_height ()) - 2) 
    ((parconsole_width ()) - 2) 1 1)

let console = ref (make_console ())

let make_searchstatus () =
    (NcursesML.window_create !parsearchstatus ((parconsole_height ()) - 2)
    ((parconsole_width ()) - 2) 1 1)

let searchstatus = ref (make_searchstatus ())

let status_stack = ref []

let make_status () = 
    ref (NcursesML.window_create !parstatus 5 (!columns - 2) 1 1)

let status = make_status ()

let error = information
let output = information


let emit_tag window str bool =
    match str with
    | "b" -> NcursesML.attr_bold window bool
    | "u" -> NcursesML.attr_underline window bool
    | "c:red" -> NcursesML.use_red window bool
    | "c:green" -> NcursesML.use_green window bool
    | "c:blue" -> NcursesML.use_blue window bool
    | "c:cyan" -> NcursesML.use_cyan window bool
    | "c:yellow" -> NcursesML.use_yellow window bool
    | "c:magenta" -> NcursesML.use_magenta window bool
    | "c:white" -> NcursesML.use_white window bool
    | "c:red_whit" -> NcursesML.use_red_in_whit window bool
    | "c:green_whit" -> NcursesML.use_green_in_whit window bool
    | "c:blue_whit" -> NcursesML.use_blue_in_whit window bool
    | "c:cyan_whit" -> NcursesML.use_cyan_in_whit window bool
    | "c:yellow_whit" -> NcursesML.use_yellow_in_whit window bool
    | "c:magenta_whit" -> NcursesML.use_magenta_in_whit window bool
    | "c:black_whit" -> NcursesML.use_black_in_whit window bool
    | _ -> ()


(* Scrollback buffer *)
type buffer_item =
    | SBstring of string
    | SBtag of (string * bool)
type sbbuffer =
        { mutable sblen : int;          (** max lines in sb buffer *)
          mutable sbhead : int;         (** index of latest line posted *)
          mutable sbbuf : buffer_item list array;
          mutable sbline : buffer_item list; (** current line so far *)
          sbwindow : NcursesML.window ref;       (** window we're working on *)
          mutable sbscrolled : int option;   (** has the user scrolled? *)
          mutable sbwrapped : bool;       (** have we wrapped once? *)
        }

(** [print_line w line] prints a line of NCURSES text (strings and tags) to the
    given window *)
let rec print_line w = function
    | [] -> ()
    | item :: items ->
          begin
              match item with
              | SBstring str ->
                    NcursesML.waddstr w str
              | SBtag (str, b) ->
                    emit_tag w str b
          end;
          print_line w items

let make_sbbuffer len window =
    { sblen = len; sbhead = 0;
      sbline = [];
      sbbuf = Array.make len [];
      sbwindow = window;
      sbscrolled = None;
      sbwrapped = false;
    }


let sbnorm sb i = (i + sb.sblen) mod sb.sblen
(* TODO: concatenate strings together to reduce memory... *)
let sbadd str sb = sb.sbline <- SBstring str :: sb.sbline
let sbflush sb =
    sb.sbbuf.(sb.sbhead) <- List.rev sb.sbline;
    sb.sbhead <- (succ sb.sbhead) mod sb.sblen;
    if sb.sbhead = 0
    then sb.sbwrapped <- true;
    sb.sbline <- []
let sb_maybeflush sb =
    match sb.sbline with
    | [] -> ()
    | _ -> sbflush sb

let get_visible_rows w =
    let (pagerows, _) = NcursesML.getmaxyx w in
    pagerows - 1

let sb_scrollto sb startrow =
(*     NcursesML.waddstr console (string_of_int startrow); *)
    let pagerows = get_visible_rows !(sb.sbwindow) in

    for i = 0 to pagerows - 1 do
        let row = sbnorm sb (startrow + i) in
(*         NcursesML.waddstr sb.sbwindow (string_of_int i); *)
(*         NcursesML.waddstr sb.sbwindow "/"; *)
(*         NcursesML.waddstr sb.sbwindow (string_of_int pagerows); *)
        print_line !(sb.sbwindow) sb.sbbuf.(row);
        NcursesML.waddstr !(sb.sbwindow) "\n";
    done;
    NcursesML.refresh ~w:!(sb.sbwindow) ();

    sb.sbscrolled <- Some startrow

let sb_scrollby sb =
    let pagerows = get_visible_rows !(sb.sbwindow) in
    pagerows - 3

let sb_back_pageful ?scrollby sb =
    sb_maybeflush sb;
    let pagerows = get_visible_rows !(sb.sbwindow) in
    let scroll = match sb.sbscrolled with
    | None -> sb.sbhead - pagerows
    | Some s -> s in
    let scrollby = match scrollby with
    | None -> sb_scrollby sb
    | Some s -> s in
    let startrow = sbnorm sb (scroll - scrollby) in
    (* check for underwrap *)
    let head =
        if sb.sbwrapped
        then sb.sbhead + 1
        else 0 in
    let scrolldiff = sbnorm sb (head - startrow) in
    let startrow =
        if scrolldiff < pagerows
        then head
        else startrow in
    sb_scrollto sb startrow

let sb_forward_pageful ?scrollby sb =
    sb_maybeflush sb;
    let pagerows = get_visible_rows !(sb.sbwindow) in
    (* previous location *)
    let scroll = match sb.sbscrolled with
    | None -> sb.sbhead - pagerows
    | Some s -> s in
    (* how much to scroll *)
    let scrollby = match scrollby with
    | None -> sb_scrollby sb
    | Some s -> s in

    let scrollpos = sbnorm sb (scroll + scrollby) in

    let head' = sb.sbhead in
    let head' =
        if head' >= scrollpos
        then head'
        else head' + sb.sblen in
    let scrollpos =
        if (scrollpos + pagerows) > head'
        then sbnorm sb (head' - pagerows)
        else scrollpos in
    sb_scrollto sb scrollpos

let sb_reset sb =
    match sb.sbscrolled with
    | None -> ()
    | Some _ ->
          let (pagerows, _) = NcursesML.getmaxyx !(sb.sbwindow) in
          let scroll = sb.sbhead - pagerows in
          sb_scrollto sb scroll;
          sb.sbscrolled <- None

let re_tagon = Str.regexp "^tag-on:\\(.*\\)$"
let re_tagoff = Str.regexp "^tag-off:\\(.*\\)$"

(* Some functions to do things in the status stack *)

let formatter_functions ?window ?(sb=None) ty =
    let do_nothing w = fun () -> NcursesML.waddstr w "\n" in
    let w, newlines =
        match ty with
        | Status -> 
              let mw = 
                  match window with 
                  | Some v -> v
                  | None -> failwith "Status window"
              in
              mw, (fun () -> NcursesML.wmove mw 0 0)
        | SearchReport -> 
                !searchstatus, do_nothing !searchstatus
        | _ -> !information, do_nothing !information
    in
    let output_string str n it =
        let str = String.sub str n it in
        if Str.string_match re_tagon str 0
        then emit_tag w (Str.matched_group 1 str) true
        else if Str.string_match re_tagoff str 0
        then emit_tag w (Str.matched_group 1 str) false
        else NcursesML.waddstr w str in
    (* First the print_string function *)
    (match sb with
     | None ->
           output_string
     | Some sb ->
           (fun str n it ->
                let str = String.sub str n it in
                if str = "\n"
                then (sbflush sb; newlines ())
                else begin
                let belt =
                    if Str.string_match re_tagon str 0
                    then (let str = Str.matched_group 1 str in
                          emit_tag w str true;
                          SBtag (str, true))
                    else if Str.string_match re_tagoff str 0
                    then (let str = Str.matched_group 1 str in
                          emit_tag w (str) false;
                          SBtag (str, false))
                    else (NcursesML.waddstr w str;
                          SBstring str) in
                sb.sbline <- belt :: sb.sbline
                end
           )),
    (* Now a function to flush *)
    (fun () -> NcursesML.refresh ~w:w ()),
    (* Now a function to print newlines *)
    (match sb with
     | None -> newlines
     | Some sb ->
           (fun () -> sbflush sb;
                newlines ())),
    (* Now a function to print spaces *)
    (match sb with
     | None ->
           (fun n ->
                let str = String.make n ' ' in
                NcursesML.waddstr w str)
     | Some sb ->
           (fun n ->
                let str = String.make n ' ' in
                sbadd str sb;
                NcursesML.waddstr w str)),
    w

let tag_functions window =
    let nullunit _ = () in
    let nullstring _ = "" in
    match window with
    | None ->
          { StatusCommon.Format.mark_open_tag = nullstring;
            StatusCommon.Format.mark_close_tag = nullstring;
            StatusCommon.Format.print_open_tag = nullunit;
            StatusCommon.Format.print_close_tag = nullunit; }
    | Some window ->
          let attrs bool str =
              emit_tag window str bool;
              if bool
              then "tag-on:" ^ str
              else "tag-off:" ^ str
          in
          { StatusCommon.Format.mark_open_tag = attrs true;
            StatusCommon.Format.mark_close_tag = attrs false;
            StatusCommon.Format.print_open_tag = nullunit;
            StatusCommon.Format.print_close_tag = nullunit; }

let create_formatter ?win ?sb w = 
    let a, b, c, d, w = 
        match win with
        | None -> formatter_functions ~sb w 
        | Some v -> formatter_functions ~window:v  ~sb w
    in
    let f = StatusCommon.Format.make_formatter a b in
    StatusCommon.Format.pp_set_all_formatter_output_functions f ~out:a ~flush:b ~newline:c
    ~spaces:d;
    StatusCommon.Format.pp_set_tags f true;
    StatusCommon.Format.pp_set_formatter_tag_functions f (tag_functions (Some w));
    StatusCommon.Format.pp_set_margin f (!columns - 3);
    f, w

let info_scrollback = make_sbbuffer 400 information

let make_info_fmt () = 
    let info_fmt, _ = create_formatter ~sb:info_scrollback Information in
    info_fmt

let info_fmt = ref (make_info_fmt ())

let make_searchstatus_fmt () =
    let searchstatus_fmt, _ = create_formatter SearchReport in
    searchstatus_fmt

let searchstatus_fmt = ref (make_searchstatus_fmt ())

let create_fresh_stack line name max lm = 
    let new_window = 
        NcursesML.window_create !parstatus 1 (!columns - 2) line 1
    in
    NcursesML.wdeleteln new_window;
    NcursesML.refresh ~w:new_window ();
    let ff, new_window = 
        create_formatter ~win:new_window Status
    in
    let res = {
        line = ref line;
        formatter = ff;
        window = new_window;
        code = get_next_code ();
        max = max;
        name = name;
        achieved = ref 0;
        last_message = ref lm;
        cleanup = fun () -> NcursesML.wdeleteln new_window;
    } 
    in
    res

let find_status st = 
    List.exists (fun x -> x.code = st.code) !status_stack

let move_window res line =
    NcursesML.mvwin res.window line 1;
    res.line := line

let clear_status_subwindows () =
    List.iter (fun hd -> 
        NcursesML.werase hd.window;
        NcursesML.refresh ~w:hd.window ();
        NcursesML.delete hd.window) !status_stack;
    status_stack := []

let eliminate_status_subwindow status = 
    let rec remove_and_rebuild = function
        | [hd] ->
                if status.code = hd.code then begin
                    NcursesML.werase hd.window;
                    NcursesML.refresh ~w:hd.window ();
                    NcursesML.delete hd.window;
                    [], 1
                end else failwith "No such status in the stack"
        | h :: ((s :: t) as tl) ->
                if status.code = h.code then begin
                    NcursesML.werase h.window;
                    NcursesML.refresh ~w:h.window ();
                    NcursesML.delete h.window;
                    tl, s.code
                end else begin 
                    let others, next = remove_and_rebuild tl in
                    move_window h next;
                    (h :: others), next + 1
                end;
        | [] -> failwith "No such status in the stack"
    in
    let res, l = remove_and_rebuild !status_stack in
    status_stack := res;
    NcursesML.mvwhline !parstatus l 1 (Char.code ' ') (!columns - 2);
    NcursesML.refresh ()

let use_colors = ref false

let formatter_of_type = function
    | Status -> failwith "Invalid type of message"
    | Error 
    | Warning
    | Information 
    | Output (None, _, _) -> !info_fmt, false, (fun () -> ())
    | Output ((Some filename), do_close, fo_ls) ->
          StatusCommon.Files.openf filename fo_ls, do_close, 
          (StatusCommon.Files.closef filename) 
    | SearchReport -> 
            let closer = fun () ->
                NcursesML.wmove !searchstatus 0 0;
                NcursesML.refresh ~w:!searchstatus ()
            in
            NcursesML.werase !searchstatus;
            !searchstatus_fmt, true, closer

let last_search_report_message = 
    ref ("" : ('a, StatusCommon.Format.formatter, unit) format)

let print ty t =
    if not !are_we_parallel || 0 = !my_rank then begin
        begin
            match ty with
            | SearchReport -> last_search_report_message := t
            | _ -> ()
        end;
        let f, do_close, closer = formatter_of_type ty in
        let t =
            match ty with
            | Error -> ("@[<v 4>@{<c:red>@{<b>Error: @}@}@,@[" ^^ t ^^ "@]@]@.")
            | Warning -> ("@[<v 4>@{<c:red>@{<b>Warning: @}@}@,@[" ^^ t ^^ "@]@]@.")
            | _ -> t
        in
        StatusCommon.Format.fprintf f t;
        if do_close then closer ();
        let () = 
            match ty, StatusCommon.information_redirected () with
            | (Output (None, _, opt)), Some filename ->
                    let f = StatusCommon.Files.openf filename opt in
                    StatusCommon.Format.fprintf f t;
            | Information, Some filename
            | Warning, Some filename 
            | Error, Some filename ->
                    let f = StatusCommon.Files.openf filename [] in
                    StatusCommon.Format.fprintf f t;
            | _ -> ()
        in
        ()
    end else !slaves_deal_in_this_way ty (string_of_format t)

(* Now we will fill the requirements of the status interface *)
let send_output _ = ()

let type_io msg rank t =
    match t with
    | SearchReport | Information -> `Information (msg, rank)
    | Error -> `Error (msg, rank)
    | Warning -> `Warning (msg, rank)
    | Status -> `Status (msg, rank)
    | Output _ -> `Output (msg, rank)

let type_string = function
    | SearchReport | Information -> "@[Information@ :@ ", stderr
    | Error -> "@[Error@ :@ ", stderr
    | Warning -> "@[Warning@ :@ ", stderr
    | Status -> "@[Status@ :@ ", stderr
    | Output _ -> "", stdout

let achieved st v =
    st.achieved := v

let get_achieved st = !(st.achieved)

let message st m = st.last_message := m

let update_status st = 
    let output_to_formatter ft = 
        match st.max, !(st.achieved) with
        | None, 0 ->
                StatusCommon.Format.fprintf ft
                "@[%s\t@;@[%s@]@ @]@." st.name !(st.last_message);
        | None, n ->
                StatusCommon.Format.fprintf ft
                "@[%s\t%d\t@;@[%s@]@ @]@." st.name !(st.achieved) 
                !(st.last_message);
        | Some max, _ ->
                StatusCommon.Format.fprintf ft
                "@[%s\t%d of %d@;@[%s@]@ @]@." st.name !(st.achieved) 
                max !(st.last_message)
    in
    output_to_formatter st.formatter;
    match StatusCommon.information_redirected () with
    | Some filename ->
            let f = StatusCommon.Files.openf filename [] in
            output_to_formatter f
    | None -> ()

let report st =
    match !verbosity with
    | `None -> ()
    |    _  ->
        assert (find_status st);
        update_status st

let user_message t v =
    match !verbosity,t with
    |    _ , Error    -> print t (StatusCommon.string_to_format (v ^ "@."))
    |    _ , Output _ -> print t (StatusCommon.string_to_format v)
    | `None, _        -> ()
    |    _ , _        -> print t (StatusCommon.string_to_format (v ^ "@."))

let create name max lm =
    match !status_stack with
    | [] -> (* We create a fresh new one *)
        user_message Information ("Starting " ^ name);
        let res = create_fresh_stack 1 name max lm in
        update_status res;
        status_stack := [res];
        res
    | hd :: _ ->
        let res = create_fresh_stack (!(hd.line) + 1) name max lm in
        status_stack := res :: !status_stack;
        res


let do_output_table t v =
    if not !are_we_parallel || 0 = !my_rank then begin
        let a, b, c = formatter_of_type t in
        StatusCommon.Tables.output a b c v;
        let () = 
            match StatusCommon.information_redirected () with
            | None -> ()
            | Some filename ->
                    let f = StatusCommon.Files.openf filename [] in
                    StatusCommon.Tables.output f false (fun () -> ()) v
        in
        ()
    end else ()


let output_table t v =
    match t with
    | Status -> ()
    | Output _ -> do_output_table t v
    | _ ->
            do_output_table t v;
            print t (StatusCommon.string_to_format "@.")

let finished st =
    match !verbosity with
    | `None -> ()
    | _ -> 
        if 1 = !(st.line) then
            user_message Information ("Finished " ^ st.name);
        assert (find_status st);
        st.last_message := "Finished\n";
        st.achieved := 
            begin match st.max with
            | None -> !(st.achieved)
            | Some v -> v
            end;
        if !use_colors then begin
            NcursesML.use_red st.window true;
            update_status st;
            NcursesML.use_red st.window false;
        end else 
            update_status st;
        eliminate_status_subwindow st

let full_report ?msg ?adv st =
    begin match msg with
    | None -> ()
    | Some v -> 
            st.cleanup ();
            st.last_message := v
    end;
    begin match adv with
    | None -> ()
    | Some v -> st.achieved := v
    end;
    report st

(* Some functions for input *)
exception End_of_line 

(* Store a command history *)
let command_history = ref Dequeue.empty
let store_command str =
    command_history := Dequeue.push_front !command_history str

let list_of_string str =
    let len = String.length str in
    let rec r i acc =
        if i = len
        then acc
        else r (succ i) (str.[i] :: acc)
    in r 0 []
let string_of_char res =
    String.make 1 res
    (*let c = (Char.escaped res) in c*)
let list_of_string str = List.map string_of_char (list_of_string str)

let error_location f t =
    (* We add the number of characters of the prompt *)
    let f = f + 5
    and t = t + 5 in
    let str = String.make (t + 1) ' ' in
    for i = f to t - 1 do
        str.[i] <- '^';
    done;
    str.[t] <- '\n';
    NcursesML.waddstr !console str

let keyup = NcursesML.keyup ()
let keydown = NcursesML.keydown ()
let keytab = 9
let keybackspace = NcursesML.keybackspace ()
let keyppage = NcursesML.keyppage ()
let keynpage = NcursesML.keynpage ()
let keyCp = 16
let keyCn = 14
let keyright = NcursesML.keyright ()
let keyleft = NcursesML.keyleft ()
let resize = NcursesML.resize ()
let err = NcursesML.error ()

let print_prompt () =
    NcursesML.use_cyan !console true;
    NcursesML.attr_underline !console true;
    NcursesML.waddstr !console "poy";
    NcursesML.attr_underline !console false;
    NcursesML.waddstr !console "> ";
    NcursesML.use_cyan !console false;
    NcursesML.refresh ~w:!console ();
    5                               (* length of prompt *)

let draw_screen_items () =
    if not !are_we_parallel || 0 = !my_rank then begin
        NcursesML.halfdelay 1;
        NcursesML.keypad ~w:!console true;
        NcursesML.echo false;
        NcursesML.box !parinformation 0 0;
        NcursesML.box !parstatus 0 0;
        NcursesML.box !parsearchstatus 0 0;
        NcursesML.box !parconsole 0 0;
        NcursesML.mvwaddstr !parinformation 0 0 " POY Output  ";
        NcursesML.mvwaddstr !parstatus 0 0 " Current Job  ";
        NcursesML.mvwaddstr !parconsole 0 0 " Interactive Console  ";
        NcursesML.mvwaddstr !parsearchstatus 0 0 " State of Stored Search  ";
        NcursesML.idlok !information true;
        NcursesML.scrollok !information true;
        NcursesML.scrollok !console true;
        use_colors := NcursesML.has_colors ();
        if !use_colors then begin 
            NcursesML.start_color ();
            NcursesML.red_in_black ();
        end;
        NcursesML.refresh ~w:!parinformation ();
        NcursesML.refresh ~w:!parstatus ();
        NcursesML.refresh ~w:!parsearchstatus ();
        NcursesML.refresh ~w:!parconsole ();
        NcursesML.refresh ~w:!information ();
        NcursesML.refresh ~w:!status ();
        NcursesML.refresh ~w:!searchstatus ();
        NcursesML.refresh ~w:!console ();
    end else ()

let sigwinchhandler x =
    NcursesML.resize_handler ();
    let c = NcursesML.columns () 
    and l = NcursesML.lines () in
    let changed =  c <> !columns || l <> !lines in
    set_columns c;
    set_lines l;
    changed

let redraw_screen () =
    ignore (sigwinchhandler ());
    NcursesML.cbreak ();
    NcursesML.finalize ();
    NcursesML.delete !parinformation;
    NcursesML.delete !parconsole;
    NcursesML.delete !parsearchstatus;
    NcursesML.delete !parstatus;
    NcursesML.delete !information;
    NcursesML.delete !console;
    NcursesML.delete !searchstatus;
    NcursesML.werase NcursesML.stdscr;
    NcursesML.refresh ~w:NcursesML.stdscr ();
    parinformation := make_parinformation ();
    parconsole := make_parconsole ();
    parsearchstatus := make_parsearchstatus ();
    parstatus := make_parstatus ();
    information := make_information ();
    console := make_console ();
    searchstatus := make_searchstatus ();
    info_fmt := make_info_fmt ();
    searchstatus_fmt := make_searchstatus_fmt ();
    draw_screen_items ();
    sb_back_pageful info_scrollback;
    sb_forward_pageful info_scrollback;
    print SearchReport !last_search_report_message;
    NcursesML.refresh ~w:!searchstatus ();
    NcursesML.do_update ()

let using_interface = ref false

let is_interactive () = !using_interface

let ignore_keys = [ 330 (* delete key *);
                    383 (* shift+delete key *) ]

let main_loop f =
    using_interface := true;
    (* output a newline *)
    print (Output (None, false, [])) (StatusCommon.string_to_format "@\n");

    (* keep command history *)
    let cmds = ref (Dequeue.push_back !command_history "") in

    let hist_up () =
        if Dequeue.is_empty !cmds
        then ""
        else begin
            let cmd, newcmds = Dequeue.pop_front !cmds in
            let newcmds = Dequeue.push_back newcmds cmd in
            cmds := newcmds;
            cmd
        end in
    let hist_down () =
        if Dequeue.is_empty !cmds
        then ""
        else begin
            let cmd, newcmds = Dequeue.pop_back !cmds in
            let newcmds = Dequeue.push_front newcmds cmd in
            cmds := newcmds;
            cmd
        end in

    let bschar = Printf.sprintf "%c" (Char.chr 8) in
    let backspace ?(n=1) () =
        for i = 1 to n do
            let (y, x) = NcursesML.getyx !console in
            if x = 0 then 
                NcursesML.wmove !console (y - 1) ((parconsole_width ()) - 3)
            else
                NcursesML.waddstr !console bschar;
            NcursesML.wdelch !console;
            NcursesML.refresh ~w:(!console) ();
        done;
        NcursesML.refresh ~w:!console ()
    in
    let delete ?(n=1) () =
        let (y,x) = NcursesML.getyx !console in
        NcursesML.wmove !console y (x+1);
        NcursesML.wdelch !console;
        NcursesML.wmove !console y x;
        NcursesML.refresh ~w:!console ()
    in

    let fuse accbc accac =
        List.fold_left (fun acc x -> x :: acc) accbc accac
    in

    let do_keyleft bc ac = match bc with
        | h :: t ->
            let (y, x) = NcursesML.getyx !console in
            if x = 0
                then NcursesML.wmove !console (y - 1) ((parconsole_width ()) - 3)
                else NcursesML.wmove !console y (x - 1);
            NcursesML.refresh ~w:!console ();
            t, (h :: ac)
        | [] -> bc, ac
    in

    let do_keyright bc ac = match ac with
        | h :: t ->
            let (y, x) = NcursesML.getyx !console in
            if x = ((parconsole_width ()) - 3)
                then NcursesML.wmove !console (y + 1) 0
                else NcursesML.wmove !console y (x + 1);
            NcursesML.refresh ~w:!console ();
            (h :: bc),t
        | [] -> bc, ac
    in

    let did_tab = ref Begin in

    let update_tab s =
        did_tab := 
            if s then
                match !did_tab with
                | Begin -> First
                | First -> Continue
                | Continue -> Continue
            else begin Begin end
    in

    let last_string_to_try_in_tab = ref "" in

    let check_limit_of_command_autocompletion before =
        let cnt = 
            List.fold_left (fun acc ch ->
                if ch = "\"" then acc + 1
                else acc) 0 before 
        in
        (* If we have an even number, it's a command *)
        if 1 land cnt = 0 then 
            ((fun x -> 
                match x.[0] with
                | 'a' .. 'z' | '_' -> false
                | _ -> true), StatusCommon.Files.Command)
        else 
            ((fun x ->
                match x.[0] with
                | '"' -> true
                | _ -> false), StatusCommon.Files.Filename)
    in
    let rec grab_filename_prefix find_limit before accac = 
        match before with
        | h :: rest when find_limit h -> "", before
        | h :: t ->
            backspace ();
            let (y, x) = NcursesML.getyx !console in
            List.iter (NcursesML.waddstr !console) accac;
            NcursesML.waddstr !console " ";
            NcursesML.wmove !console y x;
            NcursesML.refresh ~w:!console ();
            let a, b = grab_filename_prefix find_limit t accac in
            a ^ h, b
        | [] -> "", []

    in
    let prepend_string list to_add accac = 
        let list = ref list in
        String.iter
            (fun x ->
                let c = String.make 1 x in
                NcursesML.waddstr !console c;
                let (y, x) = NcursesML.getyx !console in
                List.iter (NcursesML.waddstr !console) accac;
                NcursesML.wmove !console y x;
                NcursesML.refresh ~w:!console ();
                list :=  c :: !list)
            to_add;
        NcursesML.refresh ~w:!console ();
        !list 
    in

    let prepend_to_before find_limit new_before completed_name to_complete accac =
        let _, pre = grab_filename_prefix find_limit new_before accac in
        let to_add = if completed_name = "" then to_complete else completed_name in
        prepend_string pre to_add accac
    in

    let do_keytab before after =
        let find_limit, t = check_limit_of_command_autocompletion before in
        let str, new_before = grab_filename_prefix find_limit before after in
        let () = match !did_tab with
            | Begin -> failwith "How is this possible?"
            | First -> last_string_to_try_in_tab := str
            | Continue -> ()
        in
        let str = StatusCommon.Files.complete_filename t !last_string_to_try_in_tab in
        prepend_to_before find_limit new_before str !last_string_to_try_in_tab after, after
    in

    let rec read_line accbc accac =
        NcursesML.halfdelay 1;
        let res = NcursesML.wgetch !console in
        let () =
            if res = keytab then
                update_tab true
            else if res <> err then 
                update_tab false
        in

        (* ignore keys *)
        if List.mem res ignore_keys then begin
            read_line accbc accac

        (* Enter key *)
        end else if res = err then begin
            NcursesML.do_update ();
            if (sigwinchhandler ()) then begin
                redraw_screen ();
                ignore (print_prompt ())
            end;
            read_line accbc accac

        (* tab key *)
        end else if res = keytab then begin
            let accbc, accac = do_keytab accbc accac in
            read_line accbc accac

        end else if res = 10 then begin
            sb_reset info_scrollback;
            let acc = fuse accbc accac in
            let str = List.fold_left (fun acc x -> x ^ acc) "" acc in
            let len = List.length accac in
            let (y, x) = NcursesML.getyx !console in
            NcursesML.wmove !console y (x + len);
            NcursesML.waddstr !console ("\n");
            store_command str;
            f str

        (* foreward delete -- does not work *)
        end else if res = 330 || res = 383 then begin
            delete ();
            read_line accbc accac

        (* Backspace *)
        end else if res = keybackspace || res = 127 || res = 8 then begin
            (if [] <> accbc then backspace ());
            let accbc = (match accbc with _ :: t -> t | [] -> []) in
            let (y, x) = NcursesML.getyx !console in
            List.iter (NcursesML.waddstr !console) accac;
            NcursesML.waddstr !console " ";
            NcursesML.wmove !console y x;
            NcursesML.refresh ~w:!console ();
            read_line accbc accac

        (* left key manipulates cursor position *)
        end else if res = keyleft then begin
            let accbc, accac = do_keyleft accbc accac in
            read_line accbc accac

        (* right key manipulates cursor position *)
        end else if res = keyright then begin
            let accbc, accac = do_keyright accbc accac in
            read_line accbc accac

        (* Up in command history *)
        end else if res = keyCp then begin
            backspace ~n:(List.length accbc) ();
            let hist = hist_up () in
            NcursesML.waddstr !console hist;
            read_line (list_of_string hist) []

        (* Down in command history *)
        end else if res = keyCn then begin
            backspace ~n:(List.length accbc) ();
            let hist = hist_down () in
            NcursesML.waddstr !console hist;
            read_line (list_of_string hist) []

        (* Scroll up one line *)
        end else if res = keyup then begin
            sb_back_pageful ~scrollby:1 info_scrollback;
            read_line accbc accac

        (* Scroll down one line *)
        end else if res = keydown then begin
            sb_forward_pageful ~scrollby:1 info_scrollback;
            read_line accbc accac

        (* Page down *)
        end else if res = keynpage then begin
            sb_forward_pageful info_scrollback;
            read_line accbc accac

        (* Page up *)
        end else if res = keyppage then begin
            sb_back_pageful info_scrollback;
            read_line accbc accac

        (* resize window *)
        end else if res = resize then begin
            redraw_screen ();
            ignore (print_prompt ());
            read_line accbc accac

        (* Normal input *)
        end else if NcursesML.isprint res then begin
            try let c = if res = 92 then "\\" else String.make 1 (Char.chr res) in
                begin match accac with
                    | []  -> NcursesML.waddstr !console c;
                    | _::t-> NcursesML.waddstr !console c;
                             let (y, x) = NcursesML.getyx !console in
                             List.iter (NcursesML.waddstr !console) accac;
                             NcursesML.wmove !console y x;
                end;
                NcursesML.refresh ~w:!console ();
                read_line (c :: accbc) accac
            with | err -> read_line accbc accac
        end else begin
            if ndebug_keycode then
                NcursesML.waddstr !console (string_of_int res);
            read_line accbc accac
        end
    in
    let () = f "" in
    while true do
        ignore (print_prompt ());
        cmds := (Dequeue.push_back !command_history "");
        read_line [] [];
    done

let init () = 
    if not !are_we_parallel || 0 = !my_rank then begin
        at_exit NcursesML.finalize;
        draw_screen_items ();
        redraw_screen ();
    end else ()

let resize_sbbuffer len sbb = 
    if len > 0 then begin
        let oldlen = sbb.sblen 
        and oldarr = sbb.sbbuf in
        let dif = 
            if oldlen < len then 
                (fun x -> if x < oldlen then oldarr.(x) else [])
            else 
                (fun x -> oldarr.(x + oldlen - len))
        in
        sbb.sblen <- len;
        sbb.sbbuf <- Array.init len dif
    end else begin
        let msg = "Sorry,@ but@ setting@ the@ history@ to@ " ^ 
        string_of_int len ^ "@ is@ not@ allowed@ in@ POY.@ History@ values@ \
        below@ 1@ are@ illegal." in
        user_message Error msg
    end

let resize_history len = 
    resize_sbbuffer len info_scrollback

let () = Callback.register "ncursesml_redraw" redraw_screen
