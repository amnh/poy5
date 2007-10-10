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

type channels = 
    | StandardError
    | StandardOutput
    | File of string
    | FutureFiles

(* Each command can depend on some source of information, and can change some
* source too. The following are those classes of information, which yield a
* dependency graph between commands. *)
type dependency_class = 
    | Channel of channels     (* Input  and output files *)
    | Data                  (* Characters and terminals *)
    | Trees                 (* Trees in memory *)
    | JackBoot              (* Clade sets for jackknife *)
    | Bremer                (* Bremer support values *)
    | EntryPoint            (* The beginning of the program *)

type exploders = 
    (* If the command is present at the end of a pipeline, it can be composed
    * with other composable functions that follow, and placed inside the
    * pipeline, keeping the invariant that is intended at the very end.
    * Composable function allow us to convert a sequence of operations on trees
    * to be performed sequentially, one tree at a time.*)
    | Composable        
    (* If the command is present, then there is no posibility to build a
    * pipeline. *)
    | NonComposable
    (* This command is nice, we can run it repeated times inside a pipeline and
    * it won't change the results of the analysis. *)
    | Linnearizable
    (* A command that can start or continue a pipeline, or that can be run in
    * parallel yielding correct results *)
    | Parallelizable
    (* A command is invariant if it can be safely moved outside the pipeline,
    * provided that all the Linnearizable functions are copied behind it too. *)
    | Invariant 
    (* The `Exit command. Though not composable, it does not destroy composition
    * and pipeline generation. *)
    | ExitPoint
    | InteractivePoint


type 'a concurrent_continue = (dependency_class list * int list * 'a) list

type tree_dependencies = (dependency_class list * int list) list ref
(* The thread list this node depends on and what are those items it would depend
* on each one of them *)

type parallel_v = {
    todo_p : tree;      (* Script to be run in parallel *)
    composer : tree;    (* Script to compose the commands results *)
    next : tree;        (* Script to continue with after composing them *)
    order_p : int;      (* The position in the script of this operation *)
    weight_p : float;     (* The weight in the priority of this operation *)
}
and 'a conc_tree_v = {
    run : Methods.script; (* The script command to be run before the children *)
    exploders : exploders;(* What freedom do we have with this node *)
    children : tree list; (* The list of trees that can be run concurrently
                            after this node *)
    dep_counter : int ref; (* The number of nodes this vertex depends on, that
                            remain to be executed *)
    order_c : int;        (* The posision in the script of this operation *)
    thread : int list;    (* A unique code of the thread of this node *)
    unique : 'a;          (* Type specific to Concurrent or Tree classes *)
    weight_c : float;     (* The weight in the priority of this operation *)
}
and tree =
    | Parallel of parallel_v
    | Concurrent of (tree concurrent_continue) conc_tree_v
    | Tree of tree_dependencies conc_tree_v
    | InteractiveState of int

val analyze : Methods.script list -> Methods.script list 
val almost_analyze : Methods.script list -> tree 
val script_to_string : Methods.script -> string
val explain_tree : string option -> Methods.script list -> unit
val maketree : Methods.script list -> tree list
val parallel_analysis : int -> int -> Methods.script list -> Methods.script list 
val simplify_store_set : Methods.script list -> Methods.script list
