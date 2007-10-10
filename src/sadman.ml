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

(* $Id: sadman.ml 1644 2007-02-14 19:05:47Z andres $ *)
let () = SadmanOutput.register "Sadman" "$Revision: 1644 $"


(** sadman.ml *)

type prop = (string * string)

(* [print_param channel element-name (prop name, prop value)] *)
let print_param ch str (name, value) =
    output_string ch ("    <" ^ str ^ " name=\"" ^ name
                      ^ "\" value=\"" ^ value ^ "\" />\n")

(* [print_block block_name item_name channel list-of-properties] *)
let print_block head item ch list =
    output_string ch ("<" ^ head ^ ">\n");
    List.iter (print_param ch item) list;
    output_string ch ("</" ^ head ^ ">\n")

let print_features = print_block "Features" "Feature"
let print_results = print_block "Results" "Result"

let prefix str =
    List.map (fun (a, b) -> (str ^ "." ^ a, b))

type 'a op = unit -> (prop list * prop list * 'a)

let wrap prop fn () =
    let res, a = fn () in
    (prop, res, a)

let run oc op =
    let features, results, a = op () in
    print_features oc features;
    print_results oc results;
    a

let start = SadmanOutput.start
let finish = SadmanOutput.finish
