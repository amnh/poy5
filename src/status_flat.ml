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

let () = SadmanOutput.register "" "$Revision: 2554 $"

exception Illegal_update
let () = SadmanOutput.register "Status_flat" "$Revision: 2554 $"

let _ = Format.pp_set_margin Format.std_formatter 78

(* So this is my small attempt to have the status lines working for the new POY,
* we _need_ this! *)

type formatter_output = StatusCommon.formatter_output

type c = SearchReport | Status | Error | Information 
         | Output of (string option * bool * formatter_output list)

let verbosity : [`Low | `Medium | `High ] ref = ref `High

let set_verbosity x = verbosity := x

let get_verbosity () = !verbosity

let stdferr = Format.formatter_of_out_channel stderr
let stdfout = Format.formatter_of_out_channel stdout
let _ = Format.pp_set_margin stdfout 0

let build_prefix_suffix title =
    ("@[" ^ title ^ "@ :@ @["), "@]@]@\n%!"

let type_string = function
    | SearchReport
    | Information ->
            let prefix, suffix = build_prefix_suffix "Information" in
            prefix, stdferr, (fun () -> ()), suffix
    | Error ->
            let prefix, suffix = build_prefix_suffix "Error" in
            prefix, stdferr, (fun () -> ()), suffix
    | Status -> 
            let prefix, suffix = build_prefix_suffix "Status" in
            prefix, stdferr, (fun () -> ()), suffix
    | Output (file, _, fo_ls) -> 
            let ch, f = 
                match file with
                | Some f -> 
                      StatusCommon.Files.openf f fo_ls, (fun () -> ())
                | None -> stdfout, fun () -> ()
            in
            "", ch, f, ""

let parallel = ref false

let using_interface = ref false

let is_interactive () = !using_interface

IFDEF USEREADLINE THEN
external get_line : unit -> string = "rl_CAML_gets"
ELSE
let get_line () =
    print_string "poy> ";
    flush stdout;
    let res = input_line stdin in
    print_newline ();
    res
END

let main_loop f = 
    using_interface := true;
    let _  = f "" in
    while true do
        flush stdout;
        flush stderr;
        let str = get_line () in
        f str
    done

let to_do_if_parallel = ref (fun t m ->
    let msg, ch, f, append = type_string t in
    let msg = (StatusCommon.string_to_format ((msg ^ m ^ append))) in
    Format.fprintf ch msg;
    begin 
        match t with
        | Output _ -> ()
        | _ -> Format.pp_print_flush ch ()
    end;
    let _ = 
        match StatusCommon.information_redirected () with
        | Some filename ->
              let f = StatusCommon.Files.openf filename [] in
              Format.fprintf f msg; 
        | None -> ()
    in
    f ())

let is_parallel what_to_do =
    parallel := true;
    match what_to_do with 
    | Some what_to_do -> 
            to_do_if_parallel := what_to_do
    | None -> ()

let my_rank = ref 0 

let rank i = my_rank := i

let type_io msg rank t =
    match t with
    | SearchReport
    | Information -> `Information (msg, rank)
    | Error -> `Error (msg, rank)
    | Status -> `Status (msg, rank)
    | Output _ -> `Output (msg, rank)

type status = {
    code : int;
    max : int;
    name : string;
    mutable achieved : int;
    mutable last_message : string;
}

let channels = ref [stdout]

let send_output lst = 
    channels := lst

let code =
    let counter = ref 0 in
    fun () ->
        incr counter;
        !counter

let create name max message = 
    let max = 
        match max with
        | None -> 0
        | Some v -> v 
    in
    {
        code = code ();
        achieved = 0;
        max = max;
        name = name;
        last_message = message;
    }

let achieved st v =
    st.achieved <- v

let get_achieved st = st.achieved

let message st m =
    st.last_message <- m

let single_channel_report st ch =
    let str = 
        if st.max > 0 then
            Printf.sprintf "@[%s@ :@ @[%d@ of@ %d@ -- %s@]@]" st.name 
            st.achieved st.max st.last_message
        else 
            Printf.sprintf "@[%s@ :@ @[%d@ %s@]@]" st.name st.achieved 
            st.last_message
    in
    !to_do_if_parallel Status str

let report st =
    List.iter (single_channel_report st) !channels

let single_channel_finished st =
    let s = 
        Printf.sprintf "@[%s@ Finished@]" st.name 
    in
    !to_do_if_parallel Status s

let finished st =
    single_channel_finished st

let single_user_message c msg =
    let str = Printf.sprintf "@[%s@]" msg in
    !to_do_if_parallel c str

let full_report ?msg ?adv st =
    begin match msg with
    | Some v -> message st v
    | None -> ();
    end;
    begin match adv with
    | Some v -> achieved st v
    | None -> ()
    end;
    report st

let user_message c msg =
    single_user_message c msg

let init () = ()

let error_location a b = ()

let do_output_table t v =
    match t with
    | Output (Some f, do_close, fo_ls) ->
          let a = StatusCommon.Files.openf f fo_ls in
         StatusCommon.Tables.output a do_close (StatusCommon.Files.closef f) v 
    | _ ->
            let f = stdfout in
            let do_nothing = fun () -> () in
            StatusCommon.Tables.output f false do_nothing v;
            let _ = 
                match StatusCommon.information_redirected () with
                | None -> ()
                | Some filename ->
                      let f = StatusCommon.Files.openf filename [] in
                      StatusCommon.Tables.output f false do_nothing v 
            in
            ()

let output_table t v =
    match t with
    | Status -> ()
    | Output _ -> do_output_table t v
    | _ -> 
            do_output_table t v;
            Format.fprintf stdferr (StatusCommon.string_to_format "@.")
            

let resize_history _ = ()
let redraw_screen () = ()
let clear_status_subwindows () = ()
