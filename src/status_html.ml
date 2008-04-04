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



type c = SearchReport | Status | Warning | Error | Information 
         | Output of (string option * bool * formatter_output list)

let are_we_parallel = ref false

let my_rank = ref 0

let verbosity : ([ `Low | `Medium | `High ] ref) = ref `Low

let slaves_deal_in_this_way : (c -> string -> unit) ref = ref (fun _ _ -> ())

let output_header ch = 
    output_string ch "<html><head><style type \"text/css\">
    .red {color:red}
    .green {color:green}
    .blue {color:blue}
    .purple {color:purple}
    .maroon {color:maroon}
    .white {color:white}
    .teal {color:teal}
    .black {color:black}</style></head>
        <body>
        "

let create_channel str = 
    let x = ref (open_out str) in
    output_header !x;
    x

let in_channel = ref stdin
let out_output = 
    output_header stdout;
    ref stdout
let out_status = 
    output_header stderr;
    ref stderr
let out_current = ref stderr

let close_output = ref ("" ^^ "")

let html_formatter channel = 
    let tags_to_html_tag oc tag = 
        let open_font color = "<span class=\"" ^ color ^ "\">"
        and close_font = "</span>" in
        match tag, oc with
        | "b", `Open -> "<b>"
        | "u", `Open -> "<u>"
        | "c:red_whit", `Open
        | "c:red", `Open -> open_font "red"
        | "c:green_whit", `Open 
        | "c:green", `Open -> open_font "green"
        | "c:blue_whit", `Open
        | "c:blue", `Open -> open_font "blue"
        | "c:yellow_whit", `Open 
        | "c:yellow", `Open -> open_font "purple"
        | "c:magenta_whit", `Open
        | "c:magenta", `Open -> open_font "maroon"
        | "c:white", `Open -> open_font "white"
        | "c:cyan_whit", `Open 
        | "c:cyan", `Open -> open_font "teal"
        | "c:black_whit", `Open -> open_font "black"
        | "b", `Close -> "</b>"
        | "u", `Close -> "</u>"
        | "c:red", `Close 
        | "c:green", `Close 
        | "c:blue", `Close 
        | "c:yellow", `Close 
        | "c:magenta", `Close 
        | "c:white", `Close 
        | "c:cyan", `Close 
        | "c:cyan_whit", `Close 
        | "c:red_whit", `Close 
        | "c:green_whit", `Close 
        | "c:blue_whit", `Close 
        | "c:yellow_whit", `Close 
        | "c:magenta_whit", `Close 
        | "c:black_whit", `Close -> close_font
        | "", _ -> ""
        | a, _ -> failwith ("Unknown tag " ^ a)
    in
    let tags_handler =
        { Format.mark_open_tag = tags_to_html_tag `Open;
        Format.mark_close_tag = tags_to_html_tag `Close;
        Format.print_open_tag = (fun _ -> ());
        Format.print_close_tag = (fun _ -> ()); }
    in
    let out string be l =
        let max = be + l - 1 in
        for i = be to max do
            output_char !channel string.[i];
        done;
        flush !channel;
    in
    let newline () = output_string !channel "<br>" in
    let flush () = flush !channel 
    and indent len = 
        for i = 1 to len do
            output_string !channel "&nbsp;" 
        done
    in
    let fmt = Format.make_formatter out flush  in
    Format.pp_set_all_formatter_output_functions fmt out flush newline indent;
    Format.pp_set_tags fmt true;
    Format.pp_set_formatter_tag_functions fmt tags_handler;
    Format.pp_set_margin fmt 80;
    fmt

let status_formatter = html_formatter out_status

class status header maximum suffix = 
    let to_string maximum achieved =
        match maximum, achieved with
        | None, 0 ->
                Format.sprintf "<p>@[%s\t@;@[%s@]@ @]</p>" header suffix
        | None, n ->
                Format.sprintf "<p>@[%s\t%d\t@;@[%s@]@ @]</p>" header 
                achieved suffix
        | Some max, _ ->
                Format.sprintf "<p>@[%s\t%d of %d@;@[%s@]@ @]</p>" header 
                achieved max suffix
    in
    object

    val mutable achieved = 0
    val mutable suffix = suffix

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
        Format.fprintf status_formatter "%s" string;
        match StatusCommon.information_redirected () with
        | Some filename ->
                let f = StatusCommon.Files.openf filename [] in
                output_to_formatter f
        | None -> ()

    method destroy () = ()

end

class output_field channel = 
    let f = html_formatter channel in
    object

        val formatter = f

        method set_margin v = 
            Format.pp_set_margin formatter v

        method destroy () = ()

        method formatter = formatter

        method delete = ()

    end

class standard channel =
    object (self) inherit (output_field channel) as super
        method print (string : (unit, Format.formatter, unit) format) : unit = 
            Format.fprintf formatter string
end


(* Each of the boxes that we will be using  *)
let output = new standard out_output
let current_search = new standard out_current

let get_verbosity () = !verbosity
let set_verbosity x = verbosity := x

let create header maximum suffix = 
    new status header maximum suffix

let achieved status v = status#set_advanced v

let get_achieved status = status#get_achieved

let message status string = 
    if (not !are_we_parallel) || (0 = !my_rank) then
        status#set_message string
    else !slaves_deal_in_this_way Status string

let report status = 
    if (not !are_we_parallel) || (0 = !my_rank) then
        status#print
    else ()

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
    if (not !are_we_parallel) || (0 = !my_rank) then
        status#print
    else ()

let is_parallel rank x = 
    are_we_parallel := true;
    my_rank := rank;
    match x with
    | None -> ()
    | Some f -> 
            slaves_deal_in_this_way := f

let rank x = my_rank := x

let init () = ()

let error_location a b =  ()

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
                current_search#print ("@[<v>" ^^ t ^^ "@]@\n%!")
        | Output (None, _, _) -> 
                if !close_output <> "" then ()
                else begin
                    output#print "<pre>";
                    close_output := "</pre>"
                end;
                output#print t
        | Information ->
                output#print !close_output;
                close_output := "";
                output#print ("@[<v>" ^^ t ^^  "@]@\n%!")
        | Warning ->
                output#print !close_output;
                close_output := "";
                output#print ("@[<v 4>@{<c:red>@{<b>Warning: @}@}@,@[" ^^ t ^^
                "@]@]@.%!")
        | Error -> 
                output#print !close_output;
                close_output := "";
                output#print ("@[<v 4>@{<c:red>@{<b>Error: @}@}@,@[" ^^ t ^^
                "@]@]@.%!")
    end else !slaves_deal_in_this_way ty t

let do_output_table t v =
    if 0 = !my_rank then 
        let a, b, c =
            match t with
            | Output ((Some filename), do_close, fo_ls) ->
                    (StatusCommon.Files.openf filename fo_ls),
                    do_close,
                    (StatusCommon.Files.closef filename)
            | Warning
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
    if (not !are_we_parallel) || (0 = !my_rank) then begin
        do_output_table t v;
        let _ =
            match t with
            | Output _ -> ()
            | _ -> output#print "@." 
        in
        ()
    end else ()

let redraw_screen () = ()

let resize_history _ = ()

let clear_status_subwindows () = ()

let using_interface = ref false

let is_interactive () = !using_interface

let send_output _ = ()

let type_io msg rank t =
    match t with
    | SearchReport | Information -> `Information (msg, rank)
    | Warning -> `Warning (msg, rank)
    | Error -> `Error (msg, rank)
    | Status -> `Status (msg, rank)
    | Output _ -> `Output (msg, rank)

let main_loop f =
    f "";
    while true do
        let str = input_line !in_channel in
        f str
    done
