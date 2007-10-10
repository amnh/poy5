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

(* $Id: sadmanOutput.ml 1800 2007-05-08 01:02:29Z andres $ *)
(* Created Mon Jan 23 11:47:07 2006 (Illya Bomash) *)

(** sadmanOutput.ml *)

let do_sadman = ref true

let my_name = Unix.gethostname ()

let my_pid = Unix.getpid ()

let my_data_time = Unix.time ()

let my_id = 
    my_name ^ ";" ^ string_of_int my_pid ^ ";" ^ string_of_float my_data_time

let output_ch = ref stdout
let remote = ref false
let file_open = ref false

let establish_connection ?(port=4567) () = 
    let addr = Unix.inet_addr_of_string "172.16.63.146" in
    let _, o = Unix.open_connection (Unix.ADDR_INET (addr, port)) in
    output_ch := o;
    output_value !output_ch (1, my_id, "");
    remote := true


let output data =
    if !do_sadman then begin
        if !remote then begin
            output_value !output_ch (2, my_id, data);
            flush !output_ch
        end
        else begin
            if not !file_open then begin
                output_ch := open_out "output.xml";
                file_open := true
            end;
            output_string !output_ch data;
            flush !output_ch
        end
    end else ()


let release_connection () =
    if !remote then
        output_value !output_ch (3, my_id, "");
    if !do_sadman then
        close_out !output_ch

let start name props =
    output ("<operation name=\"" ^ name
            ^ "\">\n");
    output "<features>\n";
    List.iter (fun (k, v) ->
                   output ("    <feature name=\""
                           ^ k ^ "\" value=\""
                           ^ v ^ "\" />\n"))
        props;
    output "</features>\n"

let finish props =
    output "<results>\n";
    List.iter (fun (k, v) ->
                   output ("    <result name=\""
                           ^ k ^ "\" value=\""
                           ^ v ^ "\" />\n"))
        props;
    output "</results>\n";
    output ("</operation>\n")

let _ = 
    at_exit release_connection




(** {2 Compile time information management} *)

let module_versions = ref []

let register (m : string) (v : string) =
    module_versions := (m, v) :: !module_versions

let get_module_versions () = !module_versions

let runtime_properties =
    [("os-type", Sys.os_type);
     ("compile.ccopt", CompileFlags.ccopt);
     ("compile.ocamlc", CompileFlags.ocamlc);
     ("compile.ocamlopt", CompileFlags.ocamlopt);
    ]


let () = register "SadmanOutput" "$Revision: 1800 $"

