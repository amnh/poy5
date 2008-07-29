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

(** Module for handling POY's advanced, adaptive help system *)

let () = SadmanOutput.register "HelpIndex" "$Revision: 1644 $"

let index = Help.index

let find_with_regexp reg (a, b) = 
    try 
        let _ = Str.search_forward reg a 0 in
        true
    with
    | _ -> 
        try 
            let _ = Str.search_forward reg b 0 in
            true
        with
        | _ -> false

let output_help_item (name, description) =
    Status.user_message Status.Information description

let rec help item =
    match item with
    | None ->
            help (Some ".")
    | Some "index"
    | Some "all" ->
          List.iter output_help_item index
    | Some it ->
          try
              let item = List.assoc it index in
              output_help_item (it, item)
          with
          | Not_found ->
                  let regex = Str.regexp it in
                  match List.filter (find_with_regexp regex) index with
                  | [] -> 
                          let msg = 
                              "Could not find help file \"" ^ it ^ "\""
                          in
                          Status.user_message Status.Information msg
                  | lst -> 
                          List.iter output_help_item lst

let help_if_exists it =
    try
        if List.exists (fun (x, _) -> x = it) index then
            Status.user_message Status.Information 
            ("You can find information using the command 'help (" ^ it ^ ")'")
    with
    | Not_found -> ()
