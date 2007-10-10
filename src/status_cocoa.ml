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

external cocoa_output : string -> unit = "cocoa_CAML_output"
external cocoa_status : string -> unit = "cocoa_CAML_status"

type formatter_output = StatusCommon.formatter_output

type c = SearchReport | Status | Error | Information 
         | Output of (string option * bool * formatter_output list)

class status print_function header maximum suffix = 
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

    val iter = ()

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
        let _ = 
            let string = to_string maximum achieved in
            print_function string
        in
        match StatusCommon.information_redirected () with
        | Some filename ->
                let f = StatusCommon.Files.openf filename [] in
                output_to_formatter f
        | None -> ()

    method destroy () = ()

end

    (*
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
    *)

class output_field print_function = 
    let f = (* We have to generate the formatter for the output field *)
        (*
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
        *)
        let out string be l =
            let str = String.sub string be l in
            print_function str 
        in
        let fmt = Format.make_formatter out (fun () -> ()) in
        (*
        Format.pp_set_tags fmt true;
        Format.pp_set_formatter_tag_functions fmt tags_handler;
        *)
        fmt
    in
    object

        val formatter = f

        method set_margin = ()

        method destroy () = ()

        method formatter = formatter

        method move_to_screen : unit = ()

        method underline (a : int) (c : int) = ()
end

class standard print_function =
    object (self) inherit (output_field print_function) as super
        method print (string : (unit, Format.formatter, unit) format) : unit = 
            Format.fprintf formatter string
end


(* Each of the boxes that we will be using  *)
let output = new standard cocoa_output
let console = new standard cocoa_status
let current_search = output

(*
let status = new output_field ~packing:bottom#add false
*)

let are_we_parallel = ref false
let my_rank = ref 0
let verbosity : ([ `Low | `Medium | `High ] ref) = ref `Low
let get_verbosity () = !verbosity
let set_verbosity x = verbosity := x

let create header maximum suffix = 
    new status cocoa_output header maximum suffix

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

let send_output _ = ()

let type_io msg rank t =
    match t with
    | SearchReport | Information -> `Information (msg, rank)
    | Error -> `Error (msg, rank)
    | Status -> `Status (msg, rank)
    | Output _ -> `Output (msg, rank)

let run_function : (string -> unit) ref = ref (fun _ -> ())

let run_me str = (!run_function) str

let _ =
    Callback.register "cocoa_interpreter" run_me

external cocoa_main : unit -> unit = "cocoa_CAML_main"
let main_loop f =
    using_interface := true;
    run_function := f;
    (* We will create the file opening and handling functions *)
    if !my_rank = 0 then
        cocoa_main ()
    else 
        while true do
            f ""
        done
