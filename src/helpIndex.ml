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

(** Module for handling POY's advanced, adaptive help system *)

let () = SadmanOutput.register "HelpIndex" "$Revision: 1644 $"

let index = Help.index

let find_with_regexp reg (a, b) = 
    try ignore( Str.search_forward reg a 0 ); true
    with | _ -> 
        try ignore( Str.search_forward reg b 0 ); true
        with | _ -> false

let output_help_item description =
    Status.user_message Status.Information description

let output_help_head description =
    let description = "@[<v 2>@{<c:cyan>"^description^"@}@]" in
    Status.user_message Status.Information description

let rec help = function
    | None       ->
        help (Some "index")
    | Some "index" ->
        let build x y = x ^ "@ @[<v 2>@{<c:cyan>" ^ y ^ "@}@]@," in
        let str = "@[@{<b>Help Topics :@}@]@," in
        let str = List.fold_left (fun y (x,_) -> build y x) str index in
        Status.user_message Status.Information (str^"@,@,")
    | Some "all"   ->
        List.iter (fun x -> output_help_item (snd x)) index
    | Some it    ->
        try output_help_item (List.assoc it index)
        with | Not_found ->
            let regex = Str.regexp it in
            begin match List.filter (find_with_regexp regex) index with
                | []  -> output_help_item ("Could not find help file \""^it^"\"")
                | lst -> List.iter (fun x -> output_help_item (snd x)) lst
            end

let help_if_exists it =
    try if List.exists (fun (x, _) -> x = it) index then
        output_help_item ("You can find information using the command 'help ("^it^")'")
    with | Not_found -> ()
