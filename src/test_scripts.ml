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

(* A module to handle all the test scripts *)

module Nodes = AllDirNode.AllDirF
module Edges = Edge.LazyEdge
module TreeOps = AllDirChar.F
module CharOps = AllDirChar.CharScripting

module Parsimony = 
    Scripting.Make 
    (Nodes)
    (Edges) 
    (TreeOps)
    (CharOps)

let succeed () =
    prerr_string "PASSED";
    exit 0

let failed () =
    prerr_endline "FAILED";
    exit 1

let test_tree_cost comparator result run = 
    assert (1 = Sexpr.length run.Scripting.trees);
    let nc = Ptree.get_cost `Adjusted (Sexpr.first run.Scripting.trees) in
    comparator nc result

let test_script verifier =
    try 
        Parsimony.channel_run stdin;
        if verifier (Parsimony.get_console_run ()) then succeed ()
        else failed ()
    with
    | err -> failed ()

let test_with_costs mode cost =  
    let verifier = 
        match mode with
        | `Equal -> (fun a b -> a = b)
        | `Less -> (fun a b -> a < b)
    in
    test_script (test_tree_cost verifier cost)

let test_for_termination () = 
    test_script (fun _ -> true)

let tree_cost : float option ref = ref None
let tree_cost_less : float option ref = ref None

let parse_list = [
    ("-c", Arg.Float (fun str -> tree_cost := Some str), 
    "Verify that the cost of the tree at the end equals the specified value");
    ("-cl", Arg.Float (fun str -> tree_cost_less := Some str),
    "Verify that the cost of the tree at the end is less than the specified \ 
    value")
]

let usage = 
    "poy_test [OPTIONS] filename ..."

let anon_fun _ = ()

let _ = 
    Arg.parse parse_list anon_fun usage

let _ = 
    match !tree_cost, !tree_cost_less with
    | None, None -> test_for_termination ()
    | _, Some v -> test_with_costs `Less v
    | Some v, _ -> test_with_costs `Equal v
