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

exception Illegal_update

type formatter_output = StatusCommon.formatter_output

let default_font = "courier 11"

type c = SearchReport | Status | Error | Information 
         | Output of (string option * bool * formatter_output list)

let _ = GtkMain.Main.init ()

(* The functions we will call every time a menu item is called *)
let open_f = ref (fun () -> ())
let run_f = ref (fun () -> ())
let help_f = ref (fun () -> ())

(* The main window *)
let window = GWindow.window ~width:500 ~height:300 ~title:"POY 4" ()

let mainbox = GPack.paned `VERTICAL ~packing:window#add ()
let lower = GPack.vbox ~packing:mainbox#pack2 ()
let mainbox = 
    let mainbox = GPack.vbox ~packing:mainbox#pack1 () in
    (* Time to build the menus *)
    let file_entries = [
        `I ("Open", fun () -> (!open_f) ());
        `I ("Run", fun () -> (!run_f) ());
        `S;
        `I ("Quit", GMain.Main.quit)]

    and help_entries = [ `I ("Help", fun () -> (!help_f) ()) ] in
    let menubar = 
        GMenu.menu_bar ~packing:(mainbox#pack ~fill:false
        ~expand:false) ()
    in
    menubar#misc#modify_font_by_name default_font;
    let create_menu label = 
        let item = GMenu.menu_item ~label ~packing:menubar#append () in
        GMenu.menu ~packing:item#set_submenu ()
    in
    let _ =
        let menu = create_menu "File" in
        GToolbox.build_menu menu ~entries:file_entries
    in
    let _ =
        let menu = create_menu "Help" in
        GToolbox.build_menu menu ~entries:help_entries
    in
    let label = GMisc.label ~markup:"<tt>Output</tt>" ~packing:(mainbox#pack
    ~fill:false ~expand:false) () in
    label#set_use_markup true;
    mainbox

let median = GPack.hbox ~homogeneous:false ~packing:lower#add ()
let median_right = GPack.vbox ~homogeneous:false ~packing:median#add ()
let median_left = GPack.vbox ~homogeneous:false ~packing:median#add ()

let label_console = 
    let label = 
        GMisc.label ~markup:"<tt>Interactive Console</tt>" ~packing:(median_right#pack
        ~fill:false ~expand:false) () 
    in
    label#set_use_markup true;
    label

let label_status = 
    let label = 
        GMisc.label ~markup:"<tt>Search Status</tt>" ~packing:(median_left#pack ~fill:false
        ~expand:false) () 
    in
    label#set_use_markup true;
    label

(* All the crap for the list of status *)
let cols = new GTree.column_list
let status_column = cols#add Gobject.Data.string
let model = GTree.list_store cols 
let list_view = 
    let view = GTree.view ~height:100 ~model ~packing:(lower#pack ~fill:true
    ~expand:false) () in
    let col = GTree.view_column ~renderer:(GTree.cell_renderer_text [],
    ["text", status_column]) () in
    let _ = view#append_column col in
    let _ = view#misc#modify_font_by_name default_font in
    view


class status header maximum suffix = 
    let to_string maximum achieved =
        match maximum, achieved with
        | None, 0 ->
                Format.sprintf "@[%s\t@;@[%s@]@ @]" header suffix
        | None, n ->
                Format.sprintf "@[%s\t%d\t@;@[%s@]@ @]" header 
                achieved suffix
        | Some max, _ ->
                Format.sprintf "@[%s\t%d of %d@;@[%s@]@ @]" header 
                achieved max suffix
    in
    object

    val mutable achieved = 0
    val mutable suffix = suffix

    val iter = 
        let iter = model#append () in
        let string = to_string maximum 0 in
        let _ = model#set ~row:iter ~column:status_column string in
        list_view#misc#draw None;
        iter

    method advanced = achieved <- succ achieved

    method set_advanced v = achieved <- v

    method set_message msg = suffix <- msg

    method get_achieved = achieved

    method print = 
        let output_to_formatter formatter =
            match maximum, achieved with
            | None, 0 ->
                    Format.fprintf formatter "@[%s\t@;@[%s@]@ @]@." header 
                    suffix
            | None, n ->
                    Format.fprintf formatter
                    "@[%s\t%d\t@;@[%s@]@ @]@." header achieved
                    suffix
            | Some max, _ ->
                    Format.fprintf formatter
                    "@[%s\t%d of %d@;@[%s@]@ @]@." header achieved
                    max suffix
        in
        let string = to_string maximum achieved in
        let _ = model#set ~column:status_column ~row:iter string in
        list_view#misc#draw None;
        match StatusCommon.information_redirected () with
        | Some filename ->
                let f = StatusCommon.Files.openf filename [] in
                output_to_formatter f
        | None -> ()

    method destroy () = 
        let _ = model#remove iter in ()

end

let tag_table = (* Now we specify out tag table *)
    let make_tag is_color name properties = 
        let tag = GText.tag ~name:name () in
        let p = 
            if is_color then 
                match Str.split (Str.regexp ":") name with
                | ["<c"; color] -> (`FOREGROUND (String.sub color 0 ((String.length
                color) - 1))) :: properties 
                | _ -> failwith "Illegal color"
            else properties
        in
        tag#set_properties p;
        tag#as_tag
    in
    let spec = window#misc#colormap in
    let gc = Gdk.Color.alloc spec (`RGB (0, 25700, 0))
    and yc = Gdk.Color.alloc spec (`RGB (65535, 35980, 0))
    and cc = Gdk.Color.alloc spec (`RGB (0, 35723, 35723))
    and mc = Gdk.Color.alloc spec (`RGB (59881, 38550, 31354)) in
    let bold = make_tag false "<b>" [`FONT "courier bold 12"]
    and terminal_type = make_tag false "<tt>" [`FONT default_font]
    and underlined = make_tag false "<u>" [`UNDERLINE `SINGLE]
    and red = make_tag false "<c:red>" [`FOREGROUND "red"]
    and green = make_tag false "<c:green>" [`FOREGROUND_GDK gc]
    and blue = make_tag false "<c:blue>" [`FOREGROUND "blue"]
    and yellow = make_tag false "<c:yellow>" [`FOREGROUND_GDK yc]
    and magenta = make_tag false "<c:magenta>" [`FOREGROUND_GDK mc]
    and white = make_tag false "<c:white>" [`FOREGROUND "white"]
    and cyan = make_tag false "<c:cyan>" [`FOREGROUND_GDK cc]
    and cyan2 = make_tag false "<c:cyan_whit>" [`FOREGROUND_GDK cc]
    and red2 = make_tag false "<c:red_whit>" [`FOREGROUND "red"]
    and green2 = make_tag false "<c:green_whit>" [`FOREGROUND_GDK gc]
    and blue2 = make_tag false "<c:blue_whit>" [`FOREGROUND "blue"]
    and yellow2 = make_tag false "<c:yellow_whit>" [`FOREGROUND_GDK yc]
    and magenta2 = make_tag false "<c:magenta_whit>" [`FOREGROUND_GDK mc]
    and black2 = make_tag false "<c:black_whit>" [] in
    let table = GText.tag_table () in
    List.iter table#add [bold; underlined; red; green; blue; yellow; magenta;
    white; cyan; red2; green2; blue2; cyan2; yellow2; magenta2; black2;
    terminal_type];
    table


class output_field ~packing is_interactive = 
    let b = GText.buffer ~tag_table () in
    let s = 
        GBin.scrolled_window ~vpolicy:`AUTOMATIC ~hpolicy:`AUTOMATIC 
        ~packing ~border_width:2 ~shadow_type:`ETCHED_IN ~show:true () 
    in
    let f = (* We have to generate the formatter for the output field *)
        let tags_names = ref [] in
        let tags_handler =
            let open_tag tag = 
                tags_names := tag :: !tags_names
            and close_tag tag = 
                match !tags_names with
                | h :: t ->
                        if h = tag then 
                            tags_names := t
                        else failwith ("Unmatched tag " ^ tag)
                | [] -> ()
            in
            let attrs prefix tag = 
                if tag = "" then ""
                else begin
                    let tag = "<" ^ tag ^ ">" in
                    (match prefix with
                    | `Open -> open_tag tag
                    | `Close -> close_tag tag);
                    ""
                end
            in
            { Format.mark_open_tag = attrs `Open;
            Format.mark_close_tag = attrs `Close;
            Format.print_open_tag = (fun _ -> ());
            Format.print_close_tag = (fun _ -> ()); }
        in
        let out string be l =
            let str = String.sub string be l in
            b#insert ~tag_names:("<tt>" :: !tags_names) str 
        in
        let fmt = Format.make_formatter out (fun () -> ()) in
        Format.pp_set_tags fmt true;
        Format.pp_set_formatter_tag_functions fmt tags_handler;
        fmt
    in
    let margin_size formatter textbox =
            let alc = textbox#misc#allocation in
            let cols = alc.Gtk.width / 8 in
            Format.pp_set_margin formatter cols;
    in
    object

        val scrollwin = s
        val formatter = f

        val textbox = 
            let box =
                GText.view ~buffer:b ~editable:is_interactive
                ~cursor_visible:is_interactive ~packing:s#add ~show:true ()
            in
            box#set_wrap_mode `CHAR;
            margin_size f box;
            box#misc#modify_font_by_name default_font;
            box

        method set_margin = margin_size formatter textbox

        method scrollwin = scrollwin

        method destroy () =
            textbox#destroy ();
            scrollwin#destroy ()

        method formatter = formatter

        method last_line offset = 
            let stop = b#get_iter_at_mark `INSERT  in
            let line = stop#line in
            let start =
                let start = stop#copy in
                let start = start#set_line (line - 1) in
                start#forward_chars offset
            in
            b#get_text ~start ~stop ()
            
        method textbox = textbox

        method delete =
            let a, b = textbox#buffer#bounds in
            textbox#buffer#delete ~start:a ~stop:b;

        method move_to_screen : unit =
            (textbox#scroll_mark_onscreen `INSERT)

        method underline a c =
            let stop = b#get_iter_at_mark `INSERT in
            let line = stop#line in
            let start = stop#copy in
            let start = start#set_line (line - 1) in
            let stop = start#copy in
            let start = start#forward_chars (4 + a) in
            let stop = stop#forward_chars (4 + c) in
            let text = b#get_text ~start ~stop () in
            b#delete start stop;
            b#insert ~iter:start ~tag_names:["<u>"] text
end

class standard ~packing is_interactive = 
    object (self) inherit (output_field ~packing is_interactive) as super
        method print (string : (unit, Format.formatter, unit) format) : unit = 
            Format.fprintf formatter string;
            self#textbox#misc#draw None
end


(* Each of the boxes that we will be using  *)
let output = new standard ~packing:(mainbox#pack ~fill:true ~expand:true) false
let console = new standard ~packing:(median_right#pack ~fill:true ~expand:true) true
let current_search = 
    new standard ~packing:(median_left#pack ~fill:true ~expand:true) false

(*
let status = new output_field ~packing:bottom#add false
*)

let are_we_parallel = ref false
let my_rank = ref 0
let verbosity : ([ `Low | `Medium | `High ] ref) = ref `Low
let get_verbosity () = !verbosity
let set_verbosity x = verbosity := x

let create header maximum suffix = 
    new status header maximum suffix

let achieved status v = status#set_advanced v

let get_achieved status = status#get_achieved

let slaves_deal_in_this_way : (c -> string -> unit) ref = ref (fun _ _ -> ())

let message status string = status#set_message string

let report status = status#print

let finished status = status#destroy ()

let full_report ?msg ?adv status =
    let _ =
        match msg with
        | None -> ()
        | Some msg -> status#set_message msg
    in
    let _ =
        match adv with
        | None -> ()
        | Some adv -> status#set_advanced adv
    in
    status#print

let is_parallel x = 
    are_we_parallel := true;
    match x with
    | None -> ()
    | Some f -> slaves_deal_in_this_way := f

let rank x = my_rank := x

let init () = ()

let error_location a b = 
    console#underline a b

let user_message ty t =
    if not !are_we_parallel || 0 = !my_rank then begin
        let t = StatusCommon.string_to_format t in
        match ty with
        | Output ((Some filename), do_close, fo_ls) ->
                let formatter = 
                    StatusCommon.Files.openf filename fo_ls
                in
                Format.fprintf formatter t;
                if do_close then StatusCommon.Files.closef filename ();
        | Status
        | SearchReport -> 
                current_search#delete;
                current_search#print (t ^^ "%!");
        | Output (None, _, _) -> 
                output#print t;
                output#move_to_screen
        | Information ->
                output#print (t ^^  "@\n%!");
                output#move_to_screen
        | Error -> 
                output#print ("@[<v 4>@{<c:red>@{<b>Error: @}@}@,@[" ^^ t ^^
                "@]@]@.%!");
                output#move_to_screen
    end else !slaves_deal_in_this_way ty t

let do_output_table t v =
    if not !are_we_parallel || 0 = !my_rank then 
        let a, b, c =
            match t with
            | Output ((Some filename), do_close, fo_ls) ->
                    (StatusCommon.Files.openf filename fo_ls),
                    do_close,
                    (StatusCommon.Files.closef filename)
            | Error
            | Information
            | Output (None, _, _) -> 
                    output#formatter,
                    false,
                    (fun () -> ())
            | Status
            | SearchReport -> failwith "Huh?"
        in
        let _ = StatusCommon.Tables.output a b c v in
        let _ = 
            match StatusCommon.information_redirected () with
            | None -> ()
            | Some filename ->
                    let f = StatusCommon.Files.openf filename [] in
                    StatusCommon.Tables.output f false (fun () -> ()) v
        in
        ()
    else ()

let output_table t v =
    do_output_table t v;
    let _ =
        match t with
        | Output _ -> ()
        | _ -> output#print "@." 
    in
    output#move_to_screen

let redraw_screen () = ()

let resize_history _ = ()

let clear_status_subwindows () = ()

let using_interface = ref false

let is_interactive () = !using_interface

let pressed_enter f key =
    match GdkEvent.Key.keyval key with
    | 65293 -> 
            model#clear ();
            let last_line = console#last_line 4 in
            f last_line;
            console#print "@{<c:cyan>poy>@}%!";
            true
    | int -> 
            false

let send_output _ = ()

let type_io msg rank t =
    match t with
    | SearchReport | Information -> `Information (msg, rank)
    | Error -> `Error (msg, rank)
    | Status -> `Status (msg, rank)
    | Output _ -> `Output (msg, rank)

let create_function f prependword =
    (fun () ->
        let filew = GWindow.file_selection ~title:"Open File" () in
        filew#set_select_multiple true;
        let _ = filew#connect#destroy ~callback:(fun () -> ()) in
        let _ = 
            filew#ok_button#connect#clicked 
            ~callback:(fun () -> 
                let list =  filew#get_selections in
                match list with
                | [] -> filew#destroy ()
                | _ ->
                        let filename = 
                            prependword ^ 
                            " (\"" ^ (String.concat "\", \"" list) ^
                            "\")" 
                        in
                        filew#destroy ();
                        f filename)
        in
        let _ = 
            filew#cancel_button#connect#clicked ~callback:filew#destroy 
        in
        filew#show ())

let main_loop f =
    using_interface := true;
    console#print "@{<c:cyan>poy>@}%!";
    (* We will create the file opening and handling functions *)
    open_f := create_function f "read";
    run_f := create_function f "run";
    help_f := (fun () -> f "help ()");
    f "";
    if !my_rank = 0 then
        let _ = 
            console#textbox#event#connect#key_release ~callback:(pressed_enter
            f)
        and _ =
            window#event#connect#configure ~callback:(fun _ ->
                console#set_margin;
                current_search#set_margin;
                output#set_margin; 
                false)
        and _ =
            window#connect#destroy ~callback:GMain.Main.quit
        in
        window#show ();
        GMain.Main.main ()
    else 
        while true do
            f ""
        done
