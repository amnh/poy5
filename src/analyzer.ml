(* POY 5.0. A phylogenetic analysis program using Dynamic Homologies.         *)
(* Copyright (C) 2013 Andrés Varón, Lin Hong, Nicholas Lucaroni, Ward Wheeler,*)
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

let () = SadmanOutput.register "Analyzer" "$Revision: 3482 $"

let debug = false

let (-->) a b = b a

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

let thread_table = Hashtbl.create 13

let get_dependencies list =
    let has_it = Hashtbl.create 13 in
    let not_has_it x = not (Hashtbl.mem has_it x) 
    and add_it x = Hashtbl.add has_it x 1 in
    let dep_relations_to_storage_class (a, parent) =
        let dependency_to_storage_class x = 
            let res = 
                match x with
                | Channel _ -> []
                | Data when not_has_it `Data -> [`Data]
                | Trees when not_has_it `Trees -> [`Trees]
                | JackBoot when not_has_it `Jackknife -> [`Jackknife; `Bootstrap]
                | Bremer when not_has_it `Bremer -> [`Bremer]
                | EntryPoint -> []
                | _ -> [] 
            in
            List.iter add_it res;
            res
        in
        `Set ((List.flatten (List.map dependency_to_storage_class a)), 
            try Hashtbl.find thread_table parent with
            | Not_found as err ->
                    Printf.printf "Couldn't find %s\n%!"
                    (String.concat "," (List.map string_of_int parent));
                    raise err)
    in
    List.map dep_relations_to_storage_class list

(* The most common dependency classes *)
let data = [Data]
let trees = [Trees]
let datantrees = [Data; Trees]
let output = [Channel StandardOutput]
let input files = List.map (fun x -> Channel (File x)) files
let outputf files = List.map (fun x -> Channel (File x)) files
let all ifiles ofiles = 
    JackBoot :: Data :: Trees :: Channel StandardOutput :: ((input ifiles) @
    (outputf ofiles))
let input_files = ref []
let output_files = ref []

let filename_to_list filename =
    match filename with
    | None -> [Channel StandardOutput]
    | Some x -> 
            output_files := x :: !output_files;
            [Channel (File x)]

(* Each command has certain properties by themselves, this type describes the
* properties of each command *)
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

let is_tree_dependent = function
    | `MultiStatic_Aprox _
    | `Static_Aprox _
    | `Automatic_Sequence_Partition _
    | `Automatic_Static_Aprox _ -> true
    | _ -> false

type 'a arguments_pv = 
    [ `Float of string 
    | `Int of string 
    | `String of string
    | `Labeled of (string * 'a arguments_pv)
    | `Lident of string
    | `CommandArg of 'a
    | `List of 'a arguments_pv list ]

type pc_pv = [ `Command of string * (pc_pv arguments_pv list) ]

type dependencies = 
    (dependency_class list * dependency_class list * pc_pv * exploders)

let dependency_table : (string, (pc_pv arguments_pv list -> dependencies)) Hashtbl.t = Hashtbl.create 1667

(* Takes a command, and returns a quadruple containing
* the dependency classes it depends on , the dependency classes it affects,  the
* command itself, and the exploder class the command belongs to *)
let dependency_relations (init : Methods.script) = 
    match init with
    | #Methods.tree_handling as meth ->
            let res = 
                match meth with
                | `RandomTrees _ 
                | `BestN _ 
                | `BestWithin _ 
                | `Unique -> [([Trees], [Trees], init, Composable)]
            in
            res
    | #Methods.characters_handling as meth ->
            let res = 
                match meth with
                | `RenameCharacters _
                | `AnalyzeOnlyCharacters _
                | `AnalyzeOnlyCharacterFiles _ -> 
                        [([Data; Trees], [Data; Trees], init, Linnearizable)]
            in
            res
    | #Methods.taxa_handling as meth ->
            let res = 
                match meth with
                | `SynonymsFile _
                | `Synonyms _
                | `AnalyzeOnlyFiles _
                | `AnalyzeOnly _ -> [([Data; Trees], [Data; Trees], init, Linnearizable)]
            in
            res
    | #Methods.application as meth ->
            let res = 
                match meth with
                | `Version -> [(output, output, init, Invariant)]
                | `Interactive
                | `Exit -> 
                        let all = all !input_files !output_files in
                        [(all, all, init, ExitPoint)]
                | `ChangeWDir _ ->
                        let files = outputf !output_files in
                        [(Data :: files, [Data], init, NonComposable)]
                | `PrintWDir ->
                        let files = filename_to_list None in
                        [(outputf !output_files, files, init, NonComposable)]
                | `Memory filename ->
                        let files = filename_to_list filename in
                        [([Data; Trees; JackBoot; Bremer], files, init,
                        NonComposable)]
                | `TimerInterval _
                | `Parmap _
                | `HistorySize _
                | `Redraw -> [([Data], [Data], init, Linnearizable)]
                | `Echo _ 
                | `Help _ -> 
                        let output_files = all !input_files !output_files in
                        [(output_files, output_files, init, NonComposable)]
                | `Logfile _ ->
                        let output_files = 
                            Data :: Trees :: JackBoot :: Bremer :: 
                                outputf !output_files
                        in
                        [(output_files, output_files, init, NonComposable)]
                | `Skip -> [([], [], init, Composable)]
                | `SetSeed _
                | `Alias _
                | `ClearMemory _
                | `Recover
                | `ClearRecovered ->
                        [([Data; Trees; JackBoot; Bremer], [Data; Trees; JackBoot; Bremer], init, NonComposable)]
                | `Wipe ->
                        [([EntryPoint], [Data; Trees; JackBoot; Bremer], init, NonComposable)]
                | `Algn_Newkk (*not sure here*)
                | `Algn_Normal
                | `Exhaustive_Weak
                | `Exhaustive_Strong
                | `Iterative _
                | `Normal_plus_Vitamines
                | `Optimization _
                | `Normal ->
                        [([Data; Trees], [Trees], init, NonComposable)]
                | `ReDiagnose ->
                        [([Data; Trees], [Trees], init, NonComposable)]
                | `ReDiagnoseTrees ->
                        [([Data; Trees], [Trees], init, NonComposable)]
                | `KML (_, _, filename) ->
                        let output = filename_to_list (Some filename) in
                        [([Trees; Data] @ output, output, init, NonComposable)]
                | `Graph (filename, _)
                | `Ascii (filename, _) ->
                        let output = filename_to_list filename in
                        [([Trees; Data] @ output, output, init, NonComposable)]
                | `InspectFile filename ->
                        [((input [filename]) @ output, output, init, NonComposable)]
            in
            res
    | #Methods.input as meth ->
            let rec processor (meth : [< Methods.input]) = 
                match meth with
                | `Prealigned (input, _, _) ->
                        let data, data1, compos =
                            match ((processor (input  :> Methods.input))) with
                            | [(data, data1, _, compo)] -> data, data1, compo
                            | _ -> assert false
                        in
                        [data, data1, init, compos]
                | `AnnotatedFiles files ->
                        (* TODO: Fix this analyzer step, this will be a hard
                        * stop *)
                        [(data, data, init, NonComposable)]
                | `Poyfile files
                | `AutoDetect files
                | `Nucleotides files
                | `PartitionedFile files
                | `Aminoacids (files,_)
                | `Chromosome files
                | `Genome files
                | `GeneralAlphabetSeq (files, _, _)
                | `Breakinv (files, _, _)
                | `ComplexTerminals files ->
                    let files = List.map FileStream.filename files in
                    [ (Data :: Trees :: input files, 
                      [Data; Trees; Bremer; JackBoot],
                      (meth :> Methods.script),
                      NonComposable) ]
            in
            processor meth
    | #Methods.transform as meth ->
            let res = 
                match meth with
                | `Median_Solver _
                | `Seq_to_Chrom _
                | `Custom_to_Breakinv _
                | `Annchrom_to_Breakinv _
                | `Change_Dyn_Pam _
                | `Breakinv_to_Custom _
                | `Seq_to_Kolmogorov _
                | `Fixed_States _
                | `Partitioned _
                | `Direct_Optimization _
                | `Prioritize
                | `Assign_Level _
                | `ReWeight _
                | `WeightFactor _
                | `Assign_Transformation_Cost_Matrix _
                | `OriginCost _
                | `Create_Transformation_Cost_Matrix _
                | `Assign_Affine_Gap_Cost _
                | `Assign_Tail_Cost _
                | `Prealigned_Transform _
                | `EstLikelihood _ 
                | `UseLikelihood _
                | `UseParsimony _ 
                | `Assign_Prep_Cost _ ->
                        [([Data], [Data; Trees; JackBoot; Bremer], init,
                        Linnearizable)]
                | `RandomizedTerminals 
                | `AlphabeticTerminals 
                | `MultiStatic_Aprox _
                | `Chrom_to_Seq _
                | `Static_Aprox _
                | `Search_Based _
                | `Automatic_Static_Aprox _
                | `Automatic_Sequence_Partition _ ->
                        [([Data; Trees], [Data; Trees; JackBoot; Bremer], init, 
                        NonComposable)]
            in
            res
    | #Methods.build as meth ->
            let res = 
                match meth with
                | `Branch_and_Bound _
                | `Nj
                | `Prebuilt _ ->
                        [([Data], [Data; Trees], init, Linnearizable)]
                | `Build (_, _, lst, _) when List.exists is_tree_dependent lst ->
                        [([Data; Trees], [Data; Trees], init, Parallelizable)]
                | `Build (_, (`Constraint (_, _, Some _, _)), _, _) ->
                        [([Data; Trees], [Trees], init, Parallelizable)]
                | `Build (_, (`Constraint (_, _, None, _)), _, _) ->
                        [([Data; Trees], [Trees], init, NonComposable)]
                | `Build _
                | `Build_Random _ ->
                        [([Data], [Trees], init, Parallelizable)]
            in
            res
    | #Methods.local_optimum as meth ->
            let res = 
                let `LocalOptimum (tmp) = meth in
                match tmp.Methods.tabu_join with
                | `Partition x
                when not (List.exists (function `ConstraintFile _ -> true | _ ->
                    false) x) ->
                        [([Trees], [Trees], init, NonComposable)] 
                | _ ->
                        [([Trees], [Trees], init, Parallelizable)]
            in
            res
    | (`StandardSearch _) ->
            [([Trees;Data], [Trees], init, NonComposable)]
    | #Methods.perturb_method as meth ->
            let res = 
                match meth with
                | `Ratchet _
                | `Resample _
                | `UnResample _
                | `UnRatchet 
                | `UnFixImpliedAlignments
                | `FixImpliedAlignments _ ->
                        [([Trees], [Trees], init, Linnearizable)]
            in
            res
    | `Fusing (_, _, _, _, x, _) -> [(trees, trees, init, NonComposable)]
    | `Bootstrap (it, _, _, _) 
    | `Jackknife (_, it, _, _, _) ->
            [([Data], [JackBoot], init, NonComposable)]
    | `Bremer (local_optimum, build, _, _) ->
            [([Data; Trees], [Bremer], init, NonComposable)]
    | #Methods.escape_local ->
            [([Trees], [Trees], init, Parallelizable)]
    | #Methods.runtime_store as meth -> 
            let res = 
                let all = all !input_files !output_files in
                match meth with
                | `Store _ -> [(trees, all, init, NonComposable)]
                | `Set _ -> [(trees, all, init, NonComposable)]
                | `Discard _ -> [(trees, all, init, NonComposable)]
                | `Keep_only _ -> [(trees, all, init, NonComposable)]
            in
            res
    | `ReadScript _ ->
            let all = all !input_files !output_files in
            [(all, all, init, NonComposable)]
    | `Repeat (n, comm) -> 
            (* This is special, as implies concurrency; 
            * we will treat it like nothing for now though, 
            * as nobody knows about it *)
            [(trees, data, init, NonComposable)]
    | #Methods.report as meth ->
            let res, files, isload = 
                match meth with
                | `DebugData ->
                    let fn = filename_to_list None in
                    [(Data :: fn, fn, init, Invariant)], None, false
                | `ExplainScript (_, filename)
                | `SequenceStats (filename, _)
                | `Ci (filename, _)
                | `Ri (filename, _)
                | `CompareSequences (filename, _, _, _)
                | `FasWinClad filename
                | `Nexus filename
                | `Model (filename,_)
                | `LKSites (filename,_)
                | `Script (filename,_)
                | `Dataset filename
                | `Nodes filename
                | `TerminalsFiles filename
                | `RobinsonFoulds filename
                | `CrossReferences (_, filename) ->
                        let fn = filename_to_list filename in
                        [(Data :: fn, fn, init, Invariant)], filename, false
                | `Xslt (filename, _) ->
                        let filename = Some filename in
                        let fn = filename_to_list filename in
                        [(data @ fn, fn, init, Invariant)], filename, false
                | `Topo_Selection (filename,_) ->
                    let fn = filename_to_list filename in
                    [(Data :: fn, fn, init, Invariant)], filename, false
                | `GraphicSupports (suppoutput, filename)
                | `Supports (suppoutput, filename) ->
                        let fn = filename_to_list filename in
                        let res = match suppoutput with
                            | None              ->
                                [([JackBoot; Bremer; Trees; Data] @ fn, fn, init, NonComposable)]
                            | Some `Jackknife _
                            | Some `Bootstrap _ ->
                                [([JackBoot; Trees; Data] @ fn, fn, init, Linnearizable)]
                            | Some (`Bremer _)  ->
                                [([Bremer; Trees; Data] @ fn, fn, init, Linnearizable)]
                        in
                        res, filename, false
                | `Consensus (filename, _)
                | `GraphicConsensus (filename, _) ->
                        let fn = filename_to_list filename in
                        [(datantrees @ fn, fn, init, NonComposable)], filename, false
                | `GraphicDiagnosis (_,filename) ->
                        Printf.printf "report graphic diagnosis with filename %s\n%!" filename;
                        let fn = filename_to_list (Some filename) in
                        [(datantrees @ fn, fn, init, NonComposable)], 
                        Some filename, false
                | `Diagnosis (_,filename) 
                | `AllRootsCost filename
                | `Trees (_, filename)
                | `Implied_Alignment (filename, _, _) ->
                        let fn = filename_to_list filename in
                        [(datantrees @ fn, fn, init, Linnearizable)], filename, false
                | `MstR filename
                | `TimeDelta (_, filename)
                | `TreeCosts filename
                | `KolmoMachine filename
                | `SearchStats filename
                | `TreesStats filename ->
                        let fn = filename_to_list filename in
                        [(datantrees @fn, fn, init, NonComposable)], filename, false
                | `Clades filename ->
                        let fn = filename_to_list (Some filename) in
                        [(datantrees @ fn, fn, init, Linnearizable)], Some filename, false
                | `History _ -> [([], [], init, Linnearizable)], None, false
                | `Save (filename, _) ->
                        let fn = filename_to_list (Some filename) in
                        [([Data; Trees; JackBoot; Bremer] @ fn, fn, init, NonComposable)], Some filename, false
                | `Load filename ->
                        let fn = filename_to_list (Some filename) in
                        [(fn, [Data; Trees; JackBoot; Bremer], init, NonComposable)], Some filename, true
                | `TimerInterval _
                | `Parmap _
                | `RootName _
                | `Root _ -> 
                        [([Data; Trees], all !input_files !output_files, init, Linnearizable)], None, false
            in
            if not isload then
                match files with
                | None -> ()
                | Some x ->
                        if not (List.exists (fun y -> y = x) !output_files) then
                            output_files := x :: !output_files
                        else ()
            else  ();
            res
    | `Plugin _ -> 
            let all = all [] [] in
            [(all, all, init, NonComposable)]
    | `StoreTrees
    | `UnionStored
    | `OnEachTree _
    | `Skip
    | `Entry
    | `ParallelPipeline _
    | `Barrier 
    | `GatherTrees _
    | `GatherJackknife 
    | `GatherBremer 
    | `SelectYourTrees 
    | `GatherBootstrap 
    | `GetStored -> 
            (* These are produced by the analyzer itself, so they can't occur in
            * a script and make no sense by themselves *)
            []

let compare a b = 
    match a, b with
    | EntryPoint, EntryPoint 
    | Trees, Trees
    | JackBoot, JackBoot
    | Bremer, Bremer
    | Data, Data -> 0
    | Channel a, Channel b -> compare a b
    | Bremer, _ -> -1
    | _, Bremer -> 1
    | Data, _ -> -1
    | _, Data -> 1
    | Channel _, _ -> -1
    | _, Channel _ -> 1
    | JackBoot, _ -> -1
    | _, JackBoot -> 1
    | EntryPoint, _ -> -1
    | _, EntryPoint -> 1

let flatnsort list = List.sort compare list

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

let getchannels clas (chan, _) =
    match chan, clas with
    | StandardError, StandardError
    | StandardOutput, StandardOutput -> true
    | (File str), (File str2) when str = str2 -> true
    | _ -> false

let what_I_affect ((entry, data, trees, jackboot, bremer, channel), acc) item =
    match item with
    | Channel clas -> 
            let my_children, independent =
                List.partition (getchannels clas) channel
            in
            let _, my_children = List.split my_children in
            (entry, data, trees, jackboot, bremer, independent), my_children ::
                acc
    | Data ->
            (entry, [], trees, jackboot, bremer, channel), (data :: acc)
    | Trees ->
            (entry, data, [], jackboot, bremer, channel) , trees :: acc
    | JackBoot ->
            (entry, data, trees, [], bremer, channel), jackboot :: acc
    | Bremer ->
            (entry, data, trees, jackboot, [], channel), bremer :: acc
    | EntryPoint ->
            ([], data, trees, jackboot, bremer, channel), entry :: acc

let what_I_depend_on myt (entry, data, trees, jackboot, bremer, channel) item =
    match item with
    | JackBoot -> 
            (entry, data, trees, myt :: jackboot, bremer, channel)
    | Bremer -> 
            (entry, data, trees, jackboot, myt :: bremer, channel)
    | Trees -> 
            (entry, data, myt :: trees, jackboot, bremer, channel)
    | Data -> 
            (entry, myt :: data, trees, jackboot, bremer, channel)
    | Channel clas -> 
            (entry, data, trees, jackboot, bremer, (clas, myt) :: channel)
    | EntryPoint ->
            (myt :: entry, data, trees, jackboot, bremer, channel)

let int_part_of = function
    | Parallel x -> x.order_p
    | InteractiveState x -> x
    | Concurrent x -> x.order_c
    | Tree x -> x.order_c

let get_weight = function
    | InteractiveState _ -> 0.0
    | Parallel x -> x.weight_p
    | Concurrent x -> x.weight_c
    | Tree x -> x.weight_c

let sort_list list = 
    List.sort (fun a b -> 
        let a = int_part_of a
        and b = int_part_of b in a - b) list

let sort_list2 list = 
    List.sort (fun (_, _, a) (_, _, b) -> 
        let a = int_part_of a
        and b = int_part_of b in a - b) list

(* We specify a partial order for the parallelizable operations *)
type 'a builds = [ `Build of 'a ]
type searches = [ Methods.local_optimum | Methods.escape_local ] 
type 'a all_methods = [ 'a builds | searches ]

let is_pair_of_parallels_composable a b =
    match a, b with
    | #searches, #builds 
    | #builds, #builds -> false
    | #searches, #searches 
    | #builds, #searches -> true

(* Find those paths below the tree that can be composed together. Return the tip
* of what is after those functions that can be composed, and the composed
* functions themselves *)
let rec composers cur_tip tree =
    match tree with
    | Tree x ->
            let generate_composition cur_tip =
                let vtx, tip =
                    match x.children with
                    | [chld] -> 
                            (match composers cur_tip chld with
                            | None, tip -> [], tip
                            | Some subtree, tip -> [subtree], tip)
                            | _ -> [], tree
                in
                Some (Tree { x with children = vtx }), tip
            in
            (match x.exploders with
            | Parallelizable ->
                    let new_tip = 
                        match x.run with
                        | #all_methods as x -> x 
                        | _ -> assert false
                    in
                    if is_pair_of_parallels_composable cur_tip new_tip then
                        generate_composition new_tip
                    else (Some ((InteractiveState 0)), tree)
            | Composable -> generate_composition cur_tip
            | _ -> None, tree)
    | tree -> None, tree

let rec get_tip next_command =
    (* A function that creates a tip. It will attempt to compose the components
    * of the tip *)
    let generate_this_tip compose =
        let rec attempt_to_extend_compose next_command =
            match next_command with
            | Parallel x -> (InteractiveState (-1)), next_command
            | Concurrent x -> (InteractiveState (-1)), next_command
            | InteractiveState _ -> next_command, next_command
            | Tree x ->
                    match x.exploders, x.children with
                    | Composable, lst ->
                            let lst = List.map attempt_to_extend_compose lst in
                            let a, b = List.split lst in
                            let a, b =
                                match a, b with
                                | [], [] -> (InteractiveState (-1)),
                                InteractiveState (-1)
                                | ha :: ta, hb :: tb ->
                                        if List.for_all (function
                                            InteractiveState _ -> true
                                        | _ -> false) ta then
                                            ha,
                                            Tree { x with run = `Skip; thread =
                                                [-1]; children = tb}
                                        else (InteractiveState (-1)),
                                        next_command
                                | _ -> assert false
                            in
                            Tree { x with children = [a]}, b
                    | _ -> (InteractiveState (-1)), next_command
        in
        let composer, tip = 
            if compose then 
                attempt_to_extend_compose next_command
            else (InteractiveState (-1)), next_command
        in
        (InteractiveState (-1)), composer, next_command
    in
    let recurse_to_find_tip x =
        match x.children with
        | [y] -> (* Nice, we have only one child *)
                let todo, composer, next = get_tip y in
                Tree { x with children = [todo] }, composer,
                next
        | _ -> generate_this_tip false
    in
    match next_command with
    | Parallel _
    | Concurrent _
    | InteractiveState _ -> generate_this_tip false
    | Tree x ->
            match x.exploders with
            | NonComposable -> generate_this_tip false
            | InteractivePoint 
            | Composable -> generate_this_tip true
            | ExitPoint -> generate_this_tip false
            | Linnearizable 
            | Invariant -> recurse_to_find_tip x
            | Parallelizable -> 
                    (* Hum we can parallelize this step, can we? The important
                    * thing is that if we get a tree, we should output excatly
                    * one tree again, and that we don't require all the trees
                    * from the previous command in one single group to perform
                    * the operations *)
                    match x.run with
                    | `PerturbateNSearch _ -> recurse_to_find_tip x
                    | `Build (_, y, _, _) ->
                            (match y with
                            | `Constraint (_, _, None, _) -> 
                                    (* Hum, we need to collect all of the trees
                                    * before *)
                                    generate_this_tip false
                            | _ -> recurse_to_find_tip x)
                    | `LocalOptimum (l_opt) ->
                            if List.exists (fun x -> x = `KeepBestTrees)
                            l_opt.Methods.samples then
                                generate_this_tip false
                            else 
                                (match l_opt.Methods.tabu_join with
                                | `Partition y -> 
                                        if List.exists (function `ConstraintFile
                                        _ -> true | _ -> false) y then
                                            recurse_to_find_tip x
                                        else generate_this_tip true
                                | _ -> recurse_to_find_tip x)
                    | _ -> assert false

(* Convert any parallelizable section of the tree into it. *)
let rec explode_tree tree = 
    match tree with
    | Tree x ->
            let do_not_parallelize_me_but_recurse () =
                Tree { x with children = List.map explode_tree x.children }
            in
            let parallelize_only_me () =
                let children = List.map explode_tree x.children in
                Parallel {
                    todo_p = Tree { x with children = [InteractiveState 0] };
                    composer = InteractiveState 0;
                    next = Tree { x with run = `Skip; thread = [-1]; children = children };
                    order_p = x.order_c;
                    weight_p = x.weight_c;
                }
            in
            let parallelize_this_command () =
                match x.children with
                | [] | _ :: _ :: _ -> parallelize_only_me ()
                | [y] -> 
                        match get_tip y with
                        | ((InteractiveState _), (InteractiveState _), _) -> 
                                parallelize_only_me ()
                        | (todo, composers, tip) ->
                                let todo = Tree { x with children = [todo] } in
                                let v = { 
                                    todo_p = todo; 
                                    composer = composers; 
                                    next = explode_tree tip;
                                    order_p = x.order_c;
                                    weight_p = x.weight_c;
                                } 
                                in
                                Parallel v
            in
            (match x.run with
            | `PerturbateNSearch _ -> parallelize_this_command ()
            | `Build (_, x, _, _) ->
                    (match x with
                    | `Constraint (_, _, None, _) ->
                            do_not_parallelize_me_but_recurse ()
                    | _ -> parallelize_this_command ())
            | `LocalOptimum (l_opt) ->
                    if List.exists (fun x -> x = `KeepBestTrees) l_opt.Methods.samples then
                        do_not_parallelize_me_but_recurse ()
                    else 
                        (match l_opt.Methods.tabu_join with
                        | `Partition x ->
                                if List.exists 
                                (function `ConstraintFile _ -> true 
                                    | _ -> false) x then
                                    parallelize_this_command ()
                                else do_not_parallelize_me_but_recurse ()
                        | _ -> parallelize_this_command ())
            | _ -> Tree {x with children = List.map explode_tree x.children})
    | Concurrent x ->
        let children = List.map explode_tree x.children
        and unique = 
            List.map (fun (a, b, c) -> (a, b, explode_tree c)) x.unique
        in
        Concurrent { x with children = children; unique = unique }
    | tree -> tree


    (*
let explode_tree (tree, unresolved) =
    match unresolved with
    | [] -> explode_tree tree
    | _ -> failwith "Something unresolved?"
    *)

(* Build a dependency graph by prepending the queadruple to the accumulator (the
* last pair) *)
let prepend_tree counter (depends, affects, script, exploders) (acc, roots) =
    decr counter;
    let depends = flatnsort depends
    and affects = flatnsort affects in
    let res, my_children = List.fold_left what_I_affect (acc, []) affects in
    let tree =
        match my_children with
        | [] -> InteractiveState !counter
        | _ ->
                let my_children = List.flatten my_children in
                List.iter (function
                    | Tree x -> incr x.dep_counter
                    | _ -> ()) my_children;
                let v = {
                    run = script;
                    exploders = exploders;
                    children = my_children;
                    dep_counter = ref 0;
                    order_c = !counter;
                    thread = [];
                    unique = ref [];
                    weight_c = 0.0 }
                in
                Tree v
    in
    match depends with
    | [] -> (acc, tree :: roots)
    | depends ->
            List.fold_left (what_I_depend_on tree) res depends, roots

let are_all_my_dependencies_within_my_thread mythread dependencies =
    let is_prefix a b =
        let rec is_prefix a b =
            match a, b with
            | [], _ -> true
            | ha :: ta, hb :: tb ->
                    (ha = hb) && is_prefix ta tb
            | _, [] -> failwith "Is this possible?"
        in
        (* Unfortunately the thread code comes inverted *)
        let a = List.rev a
        and b = List.rev b in
        is_prefix a b
    in
    List.fold_left (fun acc (_, x) -> acc && is_prefix x mythread) true dependencies

let extract_shared_prefix threads = 
    let get_shared (depa, a) (depb, b) =
        let rec get_shared a b =
            match a, b with
            | _, []
            | [], _ -> []
            | ha :: ta, hb :: tb ->
                    if ha = hb then ha :: get_shared ta tb
                    else []
        in
        let b = List.rev b in
        (depa @ depb), get_shared a b 
    in
    match threads with
    | (dep, h) :: t ->
            List.fold_left get_shared (dep, List.rev h) t
    | [] -> [], []

let rec is_same_list a b =
    match a, b with
    | ha :: ta, hb :: tb ->
            ha = hb && is_same_list ta tb
    | [], [] -> true
    | _, _ -> false

(* Eliminate duplications, by traversing the tree and decreasing each dependency
* counter, when it's 0, we have reached the location where a particular command
* should be called to be correct *)
let thread_id = ref 0 (* Identifiers for threads *)

let shared_dependencies a b = 
    List.filter (fun x -> List.exists (fun y -> x = y) b) a

let prepend_to_unique item pair = item.unique := pair :: !(item.unique)

let rec remove_duplications prevthread tree : 
    tree * ((dependency_class list * int list * tree) list) =
    match tree with
    | Tree x ->
            (incr thread_id;
            let my_thread = prevthread in
            let children, unresolved, yet_to_resolve = 
                let children = 
                    List.filter (function
                        | Tree y -> 
                                let ndep =
                                    try
                                        let (ndep, _, _, _) = 
                                            List.hd 
                                            (dependency_relations y.run)
                                        and (_, dep, _, _) = 
                                            List.hd 
                                            (dependency_relations x.run) 
                                        in
                                        shared_dependencies ndep dep 
                                    with 
                                    | _ -> []
                                in
                                decr y.dep_counter;
                                prepend_to_unique y (ndep, my_thread);
                                !(y.dep_counter) = 0
                        | _ -> true)
                    x.children
                in
                let children, unresolved = 
                    List.partition
                    (function
                        | Tree y ->
                                are_all_my_dependencies_within_my_thread 
                                my_thread
                                !(y.unique)
                        | _ -> true) 
                    children
                in
                let unresolved = List.map 
                    (function 
                        | (Tree y) as tree ->
                                let a, b = extract_shared_prefix !(y.unique) in
                                a, b, tree
                        | _ -> failwith "Impossible") unresolved
                in
                let children = 
                    List.map (remove_duplications my_thread) children 
                in
                let children, yet_to_resolve = List.split children in
                children, unresolved, yet_to_resolve
            in
            let new_tree = 
                Tree { x with children = children; thread = my_thread }
            in
            match List.flatten yet_to_resolve with
            | [] ->
                    new_tree, unresolved
            | yet_to_resolve ->
                    let my_thread = List.rev my_thread in 
                    match 
                        List.partition 
                        (fun (_, x, _) -> is_same_list my_thread x) 
                        yet_to_resolve
                    with
                    | [], yet_to_resolve ->
                            new_tree, (unresolved @ yet_to_resolve)
                    | mine, yet_to_resolve ->
                            let v = {
                                run = x.run;
                                exploders = x.exploders;
                                children = children;
                                dep_counter = x.dep_counter;
                                order_c = x.order_c;
                                thread = List.rev my_thread;
                                unique = mine;
                                weight_c = x.weight_c;
                            } in
                            (Concurrent v), (unresolved @ yet_to_resolve))
    | x -> x, []

let is_something something tree = 
    match tree with
    | Concurrent x -> x.exploders = something
    | Tree x -> x.exploders = something
    | _ -> true

let single_queue item = 
    let q = Queue.create () in
    Queue.add item q;
    q

let emit_name =
    let identifier = ref (-1) in
    fun thread -> 
        incr identifier; 
        let id = string_of_int !identifier ^ "__poy_analyzer_" in
        Hashtbl.add thread_table thread id;
        id

let all_dependencies = [`Jackknife; `Bootstrap; `Bremer; `Trees; `Data]

let merge_contents a b =
    List.fold_left (fun acc x ->
        if List.exists (fun y -> y = x) b then
            acc
        else x :: acc) b a

let is_poy_internal_name name =
    Str.string_match (Str.regexp "^[0-9]*__poy") name 0 

let simplify_store_set script = 
    (* We first have a function that checks if a store is used, and 
    * what elements on it are used, filters out useless stores, and 
    * for those left only the used fields remain. *)
    let rec simplify_one item (result, used_items, modified) =
        match item with
        | `Store (original_contents, name) ->
                if is_poy_internal_name name then
                    if All_sets.StringMap.mem name used_items then
                        let used_contents = 
                            All_sets.StringMap.find name used_items
                        in
                        if used_contents = original_contents then
                            (item :: result), used_items, modified
                        else
                            (((`Store (used_contents, name)) :: result),
                            used_items, true)
                    else (result, used_items, true)
                else (item :: result), used_items, modified
        | `Set ([], name) -> result, used_items, modified
        | `Set (original_contents, name) ->
                if is_poy_internal_name name then
                    let used_contents = 
                        if All_sets.StringMap.mem name used_items then
                                merge_contents original_contents 
                                (All_sets.StringMap.find name used_items)
                        else original_contents
                    in
                    let used_items = 
                        All_sets.StringMap.add name used_contents
                        used_items
                    in
                    (item :: result), used_items, modified
                else (item :: result), used_items, modified
        | `Repeat (n, script) ->
                let (script, used_items, modified) =
                    List.fold_right simplify_one script ([], used_items,
                    modified)
                in
                (`Repeat (n, script)) :: result, used_items, modified
        | `GatherTrees (script1, script2) ->
                let (script1, used_items, modified) = 
                    List.fold_right simplify_one script1 
                    ([], used_items, modified) 
                in
                let (script2, used_items, modified) = 
                    List.fold_right simplify_one script2 
                    ([], used_items, modified)
                in
                `GatherTrees (script1, script2) :: result, used_items, 
                modified
        | `OnEachTree (script1, script2) ->
                let (script1, used_items, modified) = 
                    List.fold_right simplify_one script1 
                    ([], used_items, modified) 
                in
                let (script2, used_items, modified) = 
                    List.fold_right simplify_one script2 
                    ([], used_items, modified)
                in
                `OnEachTree (script1, script2) :: result, used_items, 
                modified
        | `ParallelPipeline (c, l1, l2, l3) ->
                let stores = 
                    List.find_all (function `Store _ -> true | _ -> false) l2
                in
                let l3, used_items, modified =
                    List.fold_right simplify_one l3
                    ([], used_items, modified)
                in
                let l2, used_items, modified =
                    List.fold_right simplify_one l2
                    ([], used_items, modified)
                in
                let l1, used_items, modified =
                    List.fold_right simplify_one l1
                    ([], used_items, modified)
                in
                let l3 =
                    match l3 with
                    | `GetStored :: `Set _ :: t -> `GetStored :: (stores @ t)
                    | `GetStored :: t -> `GetStored :: (stores @ t)
                    | t -> stores @ t
                and l2 = List.filter (function `Set _ -> false | _ -> true) l2 
                and l1 = 
                    match List.rev l1 with
                    | `Store _ :: t -> List.rev t
                    | t -> List.rev t
                in
                `ParallelPipeline (c, l1, l2, l3) :: result, used_items,
                modified
        | item -> (item :: result), used_items, modified
    in
    let list_minus a b =
        let res = 
            List.filter (fun x -> not (List.exists (fun y -> x = y) b)) a
        in
        (List.length a <> List.length res), res
    in
    let has_state item lst = List.exists (fun x -> x = item) lst in
    let remove_state item lst = List.filter (fun x -> x <> item) lst in
    let get_state_from_tuple (treestate, datastate) (items, name) =
        let treestate = if has_state `Trees items then Some name else treestate
        and datastate = if has_state `Data items then Some name else datastate in 
        treestate, datastate
    in
    let is_filtered_and_modify items element name current_name =
        let tmp_current_name =
            match current_name with
            | Some x -> x
            | None -> "__poy4_bleh"
        in
        let has_it = has_state element items in
        if has_it && name = tmp_current_name then
            remove_state element items, true, 
            (if has_it then Some name else current_name)
        else items, false, (if has_it then Some name else current_name)
    in
    let stored_items = ref All_sets.Strings.empty in
    let rec remove_useless ((treestate, datastate) as acc) (list :
        Methods.script list) =
        match list with
        | [] -> [], acc, false
        | h :: t ->
                let (res : Methods.script), acc, is_modified = 
                    match (h : Methods.script) with
                    | `ParallelPipeline (c, a, b, d) ->
                            let (ar : Methods.script list), acc, am = remove_useless acc a in
                            let br, acc, bm = remove_useless acc b in
                            let ar, acc, am = remove_useless acc a in
                            let br, acc, bm = remove_useless acc b in
                            let dr, acc, dm = remove_useless acc d in
                            (`ParallelPipeline (c, ar, br, dr)), acc, 
                            (am || bm || dm)
                    | `OnEachTree (a, b) ->
                            let ar, acc, am = remove_useless acc a in
                            let br, acc, bm = remove_useless acc b in
                            (`OnEachTree (ar, br)), acc, am || bm
                    | `Repeat (n, a) ->
                            let (ar : Methods.script list), acc, am = remove_useless acc a in
                            ((`Repeat (n, ar)) : Methods.script), acc, am
                    | `GatherTrees (a, b) ->
                            let ar, acc, am = remove_useless acc a in
                            let br, acc, bm = remove_useless acc b in
                            (`GatherTrees (ar, br)), acc, am || bm
                    | `Store ((_, name) as arg) ->
                            stored_items := All_sets.Strings.add name
                            !stored_items;
                            let treestate, datastate = 
                                get_state_from_tuple acc arg
                            in
                            h, (treestate, datastate), false
                    | `Set (items, name) ->
                            if is_poy_internal_name name then
                                if All_sets.Strings.mem name !stored_items then
                                    let items, modified, datastate = 
                                        is_filtered_and_modify items `Data name datastate
                                    in
                                    let items, modified2, treestate = 
                                        is_filtered_and_modify items `Trees name
                                        treestate
                                    in
                                    (`Set (items, name)), (treestate, datastate), 
                                    modified2 || modified
                                else (`Set ([], name)), (treestate, datastate), true
                            else (`Set (items, name)), (treestate, datastate),
                            false
                    | x ->
                            match dependency_relations x with 
                            | [(_, b, _, _)] ->
                                    let treestate =
                                        if has_state Trees b then None
                                        else treestate
                                    and datastate = 
                                        if has_state Data b then None 
                                        else datastate 
                                    in
                                    x, (treestate, datastate), false
                            | [] -> x, acc, false
                            | _ -> failwith "Multiples?"
                in
                let tail, acc, is_modified2 = remove_useless acc t in
                res :: tail, acc, is_modified2 || is_modified
    in
    let rec remove_duplications list =
        match list with
        | [] -> false, []
        | `ParallelPipeline (c, a, b, d) :: t ->
                let at, a = remove_duplications a in
                let bt, b = remove_duplications b in
                let dt, d = remove_duplications d in
                let tt, t = remove_duplications t in
                at || bt || dt || tt, ((`ParallelPipeline (c, a, b, d)) :: t)
        | `OnEachTree (a, b) :: t ->
                let at, a = remove_duplications a in
                let bt, b = remove_duplications b in
                let tt, t = remove_duplications t in
                at ||  bt || tt, ((`OnEachTree (a, b)) :: t)
        | `Repeat (n, a) :: t ->
                let at, a = remove_duplications a in
                let tt, t = remove_duplications t in
                at || tt, ((`Repeat (n, a)) :: t)
        | `GatherTrees (a, b) :: t ->
                let at, a = remove_duplications a in
                let bt, b = remove_duplications b in
                let tt, t = remove_duplications t in
                at || bt || tt, (`GatherTrees (a, b)) :: t
        | h :: t ->
                match remove_duplications t with
                | changed, [] -> changed, [h]
                | changed, (i :: tt) ->
                        match h, i with
                        | (`Set (ch, nameh)), (`Set (ci, namei)) ->
                                if nameh = namei then
                                    true, `Set ((merge_contents ch ci), nameh) :: tt
                                else 
                                    let has_shared, ch = list_minus ch ci in
                                    (changed || has_shared), (`Set (ch, nameh))
                                    :: i :: tt
                        | (`Store (_, nameh)), (`Set (_, namei)) when
                            nameh = namei -> true, h :: tt
                        | _ -> changed, h :: i :: tt
    in
    let rec add_discard item (res, used) =
        match item with
        | `Set (_, name) when is_poy_internal_name name ->
                if All_sets.Strings.mem name used then
                    item :: res, used
                else 
                    let used = All_sets.Strings.add name used in
                    item :: `Discard (all_dependencies, name) :: res, used
        | _ -> item :: res, used
    in
    let rec simplify_script script =
        let result, _, mod1 = 
            List.fold_right simplify_one script ([],
            All_sets.StringMap.empty, false) 
        in
        let mod2, result = remove_duplications result in
        let result, _, mod3 = remove_useless (None, None) result in
        if mod1 || mod2 || mod3 then begin
            stored_items := All_sets.Strings.empty;
            simplify_script result
        end else 
            let res, _ =
                List.fold_right add_discard result ([], All_sets.Strings.empty)
            in
            res
    in
    simplify_script script

let remove_set lst =
    match lst with
    | (`Set _) :: t -> t
    | x -> x

let remove_trees_of_set = function
    | `Set (items, name) ->
            (`Set ((List.filter (fun x -> x <> `Trees) items), name))
    | x -> x

let remove_all_trees_from_set = List.map remove_trees_of_set

let remove_trees_from_set lst = 
    match lst with
    | h :: t -> (remove_trees_of_set h) :: t
    | [] -> []

(* The main function to convert a script into a linnear list of operations *)
let rec linearize2 queue acc = 
    while not (Queue.is_empty queue) do
        let tree = Queue.pop queue in
        match tree with
        | Concurrent x ->
                let my_name = emit_name x.thread in
                acc := (`Store (all_dependencies, my_name)) :: x.run :: !acc;
                let continue = x.unique in
                let children = x.children in
                let myqueue = Queue.create () in
                List.iter (fun x -> Queue.add x myqueue) children;
                List.iter (fun (_, _, x) -> Queue.add x queue) continue;
                linearize2 myqueue acc;
        | Tree x -> 
                let my_name = emit_name x.thread in
                let deps = 
                    try get_dependencies !(x.unique) with 
                    | Not_found as err ->
                        Printf.printf "My f thread is %s\n%!" 
                        (String.concat ", " (List.map string_of_int x.thread));
                        raise err
                in
                acc := (`Store (all_dependencies, my_name)) :: x.run :: (deps @ !acc);
                List.iter (fun y -> linearize2 (single_queue y) acc)
                          (sort_list x.children)
        | InteractiveState _ -> ()
        | Parallel y ->
                match y.todo_p with
                | Tree x ->
                        let my_name = emit_name x.thread in
                        let deps = 
                            try get_dependencies !(x.unique) with 
                            | Not_found as err ->
                                Printf.printf "My thread is %s\n%!" 
                                (String.concat ", " (List.map string_of_int x.thread));
                                raise err
                        in
                    begin match x.run with
                        | `Build (total, b, c, d) ->
                                let item = 
                                    Tree { x with run = `Build (1, b, c, d) } 
                                in
                                let todol = ref []
                                and composerl = ref []
                                and nextl = ref [] in
                                linearize2 (single_queue item) todol;
                                linearize2 (single_queue y.composer) composerl;
                                linearize2 (single_queue y.next) nextl;
                                let todol =
                                    deps @ 
                                    remove_trees_from_set (List.rev !todol)
                                and nextl = remove_trees_from_set (List.rev !nextl) 
                                and composerl = 
                                    remove_trees_from_set (List.rev
                                    (`StoreTrees :: !composerl))
                                in
                                acc := 
                                    `Store (all_dependencies, my_name) ::
                                        `ParallelPipeline (total, todol, 
                                        `UnionStored :: composerl, 
                                            `GetStored :: nextl) :: (!acc);
                        | `LocalOptimum _
                        | `PerturbateNSearch _->
                                let composerl = ref []
                                and nextl = ref []
                                and iteml = ref [] in
                                linearize2 (single_queue y.todo_p) iteml;
                                linearize2 (single_queue y.composer) composerl;
                                linearize2 (single_queue y.next) nextl;
                                let deps = (remove_all_trees_from_set (List.rev deps)) in
                                let iteml = deps @ remove_trees_from_set (List.rev !iteml)
                                and nextl = List.rev (remove_trees_from_set (List.rev !nextl)) in
                                acc := nextl @
                                    `Store (all_dependencies, my_name) ::
                                    (`OnEachTree (iteml, List.rev !composerl)) 
                                        :: (!acc)
                        | _ -> assert false
                    end
                | _ -> assert false
    done

let linearize tree _ =
    let acc = ref [] in
    linearize2 (single_queue tree) acc;
    List.rev !acc

let is_load x = 
    let check_it x =
        match x.run with
        | `Load _ 
        | `Save _ -> true
        | _ -> false
    in
    match x with
    | InteractiveState _ 
    | Parallel _ -> false
    | Concurrent x -> check_it x
    | Tree x -> check_it x

let get_order = function
    | Parallel _ -> 0
    | InteractiveState ord -> ord
    | Concurrent x -> x.order_c
    | Tree x -> x.order_c

let entry res = 
    let v = {
        run = `Entry;
        exploders = NonComposable;
        children = res;
        dep_counter = ref 0;
        order_c = 0;
        thread = [];
        unique = ref [];
        weight_c = 0.0
    } in
    Tree v

let make_root_trees script =
    if debug then
    Printf.printf "make_root_trees,script lst len=%d\n%!"
    (List.length script);
    thread_id := 0;
    script 
    --> List.map dependency_relations
    --> List.flatten
    --> fun x ->
        List.fold_right (prepend_tree (ref 0)) 
        (([], 
        [EntryPoint; Data; Trees; JackBoot; Bremer], 
        `Entry, NonComposable) :: x) 
        (([], [], [], [], [], []), []) 

let post_process_trees res =
    res 
    --> List.map (remove_duplications []) 
    --> List.map 
        (fun (a, b) ->
            a :: (List.map (fun (_, x, y) -> 
                match x with
                | [] -> y
                | _ -> failwith "What?") b))
    --> List.flatten
    --> List.map explode_tree

let maketree script =
    let _, res = make_root_trees script in 
    post_process_trees res

let almost_analyze script = 
    let prepare_if_no_exit script =
        if (List.exists (fun x -> x = `Exit) script) then script
        else script @ [`Interactive]
    in
    let (_, _, _, _, _, channels), res = 
        script --> prepare_if_no_exit --> make_root_trees
    in
    let _, channels_subtrees = List.split channels in
    let channels_subtrees = 
        List.filter is_load channels_subtrees
    in
    entry (post_process_trees (res @ channels_subtrees))

(* Take a valid script and make it a memory efficient script, if I can *)
let analyze (script : Methods.script list) : Methods.script list = 
    (*
    let prepare_if_no_exit script =
        if (List.exists (fun x -> x = `Exit) script) then script
        else script 
    in
    *)
    if debug then
    Printf.printf "Analyzer.analyze1,input listlen=%d\n%!"
    (List.length script);
    let (_, _, _, _, _, channels), res = 
        script --> make_root_trees
    in
    let _, channels_subtrees = List.split channels in
    let channels_subtrees = 
        List.filter is_load channels_subtrees
    in
    let res = entry (post_process_trees (res @ channels_subtrees)) in
    let r = simplify_store_set (linearize res []) in 
    Hashtbl.clear thread_table;
    r

let break_in_independent_sections x =
    let a, b =
        List.fold_right (fun x (cur, acc) ->
            match x with
            | `Wipe | `Set _ -> ([], ((x :: cur) :: acc)) 
            | _ -> (x :: cur), acc) x ([], [])
    in
    a :: b

let rec correct_parallel_pipelines_with_internal_transform (res : Methods.script list)= 
    match res with
    | (`ParallelPipeline (times, ((h :: t) as acc), comp, continue)) as hd :: tl ->
            let tl = correct_parallel_pipelines_with_internal_transform tl in
            let continue = correct_parallel_pipelines_with_internal_transform continue in
            (match h with
            | `Build (a, b, cost_calculation, i) ->
                    (match cost_calculation with
                    | [] -> 
                            (`ParallelPipeline (times, acc, comp, continue)) :: tl
                    | _ ->
                        let name_post = emit_name [] 
                        and name_pre = emit_name [] in
                        let dtlst = [`Data] in
                        let make_discard name = `Discard (dtlst, name) in
                        let make_set name = `Set (dtlst, name) in
                        let make_store name = `Store (dtlst, name) in
                        let continue = 
                            (make_discard name_pre) :: (make_discard name_post)
                            :: continue
                        in
                        let acc = (make_set name_post) :: (`Build (a, b, [],i)) ::
                            (make_set name_pre) :: t
                        in
                        (make_store name_pre) :: 
                        ((cost_calculation :> Methods.script list) @ 
                            ((make_store name_post) ::
                                (`ParallelPipeline (times, acc, comp, continue))
                                :: tl)))
            | `Branch_and_Bound ((_, _, _, _, cost_calculation),_)
            | `LocalOptimum { Methods.cc = cost_calculation } 
            | `PerturbateNSearch (_, _, `LocalOptimum { Methods.cc =
                cost_calculation } , _, _) ->
                    if List.exists is_tree_dependent cost_calculation then
                        let name = emit_name [] in
                        (`Store ([`Trees], name)) ::
                        `ParallelPipeline (times, (`Set ([`Trees], name)) ::
                            acc, comp, (`Discard ([`Trees], name)) ::
                                continue) :: tl
                    else 
                        (`ParallelPipeline (times, acc, comp, continue)) :: tl
            | _ -> hd :: tl)
    | hd :: tl -> hd :: (correct_parallel_pipelines_with_internal_transform tl)
    | [] -> []


let analyze script = 
    if debug then
    Printf.printf "Analyzer.analyze2 input lst len = %d\n%!"
    (List.length script);
    let scripts = break_in_independent_sections script in 
    if debug then
    Printf.printf "after break_in_independent_sections,len=%d\n%!"
    (List.length scripts);
    let res = List.flatten (List.map analyze scripts) in
    if debug then
    Printf.printf "before parallel,lstlen=%d\n%!"
    (List.length res);
    let res =
    correct_parallel_pipelines_with_internal_transform res
    in
    if debug then
    Printf.printf "Analyzer.analyze2 reslst len = %d\n%!"
    (List.length res);
    res

let script_to_string (init : Methods.script) =
    match init with
    | #Methods.tree_handling as meth ->
            let res = 
                match meth with
                | `RandomTrees _ ->
                        "@[select random trees@]"
                | `BestN None ->
                        "@[select the optimal trees@]"
                | `BestN (Some x) ->
                        "@[select the best " ^ string_of_int x ^ "@ trees@]"
                | `BestWithin float ->
                        "@[select the best within " ^ string_of_float float ^ 
                        " percent of the minimum@]"
                | `Unique -> "@[eliminate repeated trees@]"
            in
            res
    | #Methods.characters_handling as meth ->
            let res = 
                match meth with
                | `RenameCharacters _ ->
                        "@[rename the characters according to your list@]"
                | `AnalyzeOnlyCharacters _ 
                | `AnalyzeOnlyCharacterFiles _ -> 
                        "@[analyze only those characters you selected in your list@]"
            in
            res
    | #Methods.taxa_handling as meth ->
            let res = 
                match meth with
                | `SynonymsFile _ 
                | `Synonyms _ ->
                        "@[rename the terminals that you specified in your list@]"
                | `AnalyzeOnlyFiles _
                | `AnalyzeOnly _ -> 
                        "@[analyze the terminals that you specified in your list@]"
            in
            res
    | #Methods.application as meth ->
            let res = 
                match meth with
                | `Version -> 
                        "@[output the version@]"
                | `Exit -> 
                        "@[close POY@]"
                | `Interactive -> 
                        "@[wait for the user to issue some other command@]"
                | `ChangeWDir _ ->
                        "@[change my working directory@]"
                | `PrintWDir ->
                        "@[print my working directory@]"
                | `Memory _ ->
                        "@[print my memory statistics@]"
                | `KML _ ->
                        "@[produce a KML file@]"
                | `TimerInterval _ ->
                        "@[change the timer interval@]"
                | `Parmap _ ->
                        "@[change the number of cores for parmap@]"
                | `HistorySize _ ->
                        "@[change my history size@]"
                | `Redraw -> 
                        "@[redraw the screen@]"
                | `Echo _  ->
                        "@[print a user defined message@]"
                | `Help _ -> 
                        "@[print some help@]"
                | `Logfile _ ->
                        "@[change my log file@]"
                | `Alias _ ->
                        "@[name a set of characters@]"
                | `SetSeed _ ->
                        "@[change the random number generator's seed@]"
                | `ClearMemory _ ->
                        "@[cleanup the memory@]"
                | `Recover ->
                        "@[recover trees from a swap@]"
                | `ClearRecovered ->
                        "@[eliminate the trees I recovered in a swap@]"
                | `Wipe ->
                        "@[get rid of all trees and data@]"
                | `Exhaustive_Strong
                | `Exhaustive_Weak -> 
                        "@[set the cost calculation to Exhaustive Weak DO@]"
                | `Iterative _ ->
                        "@[set the cost calculation to iterative@]"
                | `Normal_plus_Vitamines ->
                        "@[set the cost calculation to normal+ DO@]"
                | `Normal ->
                        "@[set the cost calculation to normal DO@]"
                | `Optimization `Coarse _ ->
                        "@[set the optimization level to coarse@]"
                | `Optimization `None ->
                        "@[set the optimization level to off@]"
                | `Optimization `Exhaustive _ ->
                        "@[set the optimization level to exhaustive@]"
                | `Optimization `Exhaustive_dyn _ ->
                        "@[set the optimization level to dynamic exhaustive@]"
                | `Optimization `Custom _ ->
                        "@[set the optimization level to custom@]"
                | `ReDiagnose ->
                        "@[rediagnose the trees@]"
                | `ReDiagnoseTrees ->
                        "@[rediagnose the trees preserving model in likelihood@]"
                | `Graph (_, _)
                | `Ascii (_, _) ->
                        "@[output the trees in memory@]"
                | `InspectFile _ ->
                        "@[printout the metadata of a poy file@]"
                | `Algn_Newkk -> "@[do alignment with low memory mode@]"
                | `Algn_Normal -> "@[do alignment with normal mode@]"

            in
            res
    | #Methods.input -> "@[read an input file@]"
    | #Methods.transform as meth ->
            let res = 
                match meth with
                | `Median_Solver _
                | `Seq_to_Chrom _
                | `Custom_to_Breakinv _
                | `Annchrom_to_Breakinv _
                | `Change_Dyn_Pam _
                | `Chrom_to_Seq _
                | `Breakinv_to_Custom _
                | `Seq_to_Kolmogorov _
                | `Fixed_States _
                | `Partitioned _
                | `Direct_Optimization _
                | `Prioritize
                | `Assign_Level _
                | `ReWeight _
                | `WeightFactor _
                | `Assign_Transformation_Cost_Matrix _
                | `OriginCost _
                | `Create_Transformation_Cost_Matrix _
                | `Assign_Affine_Gap_Cost _
                | `Assign_Tail_Cost _
                | `Assign_Prep_Cost _ 
                | `RandomizedTerminals 
                | `AlphabeticTerminals
                | `MultiStatic_Aprox _
                | `Static_Aprox _
                | `Search_Based _
                | `Prealigned_Transform _
                | `EstLikelihood _
                | `UseParsimony _ 
                | `UseLikelihood _
                | `Automatic_Static_Aprox _
                | `Automatic_Sequence_Partition _ ->
                        "@[transform some characters@]"
            in
            res
    | #Methods.build as meth ->
            let res = 
                match meth with
                | `Branch_and_Bound _ ->
                        "@[Build trees using branch and bound@]"
                | `Nj -> "@[build trees using neighbor joining@]"
                | `Prebuilt _ ->
                        "@[load the trees from a file@]"
                | `Build _
                | `Build_Random _ ->
                        "@[build some trees from scratch@]"
            in
            res
    | #Methods.local_optimum ->
            "@[swap the trees in memory@]"
    | `StandardSearch _ ->
            "@[execute an automated search@]"
    | #Methods.perturb_method as meth ->
            let res = 
                match meth with
                | `Ratchet _ -> "@[do a ratchet@]"
                | `Resample _ -> "@[resample the characters@]"
                | `UnResample _
                | `UnRatchet 
                | `UnFixImpliedAlignments -> ""
                | `FixImpliedAlignments _ ->
                        "@[some obscure preturbation?"
            in
            res
    | `Fusing (_, _, _, _, x, _) -> 
            "@[fuse the trees I have in memory@]"
    | `Bootstrap (it, _, _, _) ->
            "@[calculate the bootstrap clades@]"
    | `Jackknife (_, it, _, _, _) ->
            "@[calculate the jackknife clades@]"
    | `Bremer (local_optimum, build, _, _) ->
            "@[calculate the bremer support values@]"
    | #Methods.escape_local ->
            "@[try to escape the local optimum@]"
(*            `PerturbateNSearch of (transform list * perturb_method **)
(*            local_optimum * int) *)
    | #Methods.runtime_store as meth -> 
            let res = 
                match meth with
                | `Store _ -> "@[store my state in a variable@]"
                | `Set _ -> "@[use a stored state as my current state@]"
                | `Discard _ -> "@[delete a stored state@]"
                | `Keep_only _ ->  ""
            in
            res
    | `ReadScript _ ->
            "@[run the script you gave me@]"
    | `Repeat (n, comm) -> 
            (* This is special, as implies concurrency; we will treat it like nothing for
            now though, as nobody knows about it *)
            ""
    | #Methods.report as meth ->
            let res = 
                match meth with
                | `ExplainScript _ ->
                        "@[explain you a script@]"
                | `SequenceStats _ ->
                        "@[report the sequence statistics@]"
                | `Ci _ -> "@[report ci@]"
                | `Ri _ -> "@[report ri@]"
                | `CompareSequences _ ->
                        "@[report the average distance between pairs of \
                        sequences@]"
                | `FasWinClad _ -> 
                        "@[report the phastwinclad file@]"
                | `Nexus _ -> 
                        "@[report the nexus file@]"
                | `DebugData -> 
                        "@[report global data to screen (debug)@]"
                | `LKSites _ -> 
                        "@[report the site loglikelihood of trees@]"
                | `Model _ -> 
                        "@[report the likelihood model@]"
                | `Script _ -> 
                        "@[report the script being run@]"
                | `Dataset _ ->
                        "@[report the current data under analysis@]"
                | `Xslt _ ->
                        "@[report the current data and trees under analysis@]"
                | `Nodes _ ->
                        "@[report all unadjusted data in tree@]"
                | `Topo_Selection _ -> 
                        "@[report the tree selection p-values@]"
                | `TerminalsFiles _ ->
                        "@[report the terminals per file@]"
                | `RobinsonFoulds _ ->
                        "@[report robinson foulds tree distance matrix@]"
                | `CrossReferences (_, _) ->
                        "@[report the cross references@]"
                | `GraphicSupports (_, _)
                | `Supports (_, _) ->
                        "@[report the support values@]"
                | `Consensus (_, _)
                | `GraphicConsensus (_, _) ->
                        "@[report the consensus@]"
                | `GraphicDiagnosis (reporttype,filename) ->
                        "@[report in " ^ filename ^ " the diagnosis@]"
                | `Diagnosis (reporttype,filename) ->
                        (match filename with
                        | None -> "@[report the diagnosis@]"
                        | Some x ->
                                "@[report in " ^ x ^ " the diagnosis@]")
                | `MstR _ ->
                        "@[report the minimum spanning tree@]"
                | `AllRootsCost filename ->
                        (match filename with
                        | None -> "@[report the minimum length of each rooting@]"
                        | Some x -> 
                                "@[report in " ^ x ^ 
                                " the minimum length of each rooting@]")
                | `Trees (_, filename) ->
                        (match filename with
                        | None -> "@[report the trees in memory@]"
                        | Some x -> 
                                "@[report in " ^ x ^ " the trees in memory@]")
                | `Implied_Alignment _ ->
                        "@[report an implied alignment@]"
                | `TimeDelta _ ->
                        "@[report the time delta@]"
                | `KolmoMachine _ ->
                        "@[report the complexity of the machine loaded@]"
                | `TreeCosts _ ->
                        "@[report the cost of the trees@]"
                | `SearchStats _ ->
                        "@[report the search results@]"
                | `TreesStats _ ->
                        "@[report the tree statistics@]"
                | `Clades _ ->
                        "@[create a clades file@]"
                | `History _ -> ""
                | `Save (_, _) ->
                        "@[save my memory in a file@]"
                | `Load _ ->
                        "@[load my memory from a file@]"
                | `RootName _
                | `Root _ -> 
                        "@[change the root@]"
            in
            res
    | `Plugin _ -> 
            "@[execute a plugin function@]"
    | `Entry -> "@[beginning of the program@]"
    | `Barrier -> "@[Wait for other processors to reach this point@]"
    | `GatherTrees _ -> "@[Exchange trees with other processes@]"
    | `GatherJackknife -> "@[Exchange Jackknife values with other processes@]"
    | `GatherBremer  -> "@[Exchange Bremer values with other processes@]"
    | `GatherBootstrap -> "@[Exchange Bootstrap values with other processes@]" 
    | `SelectYourTrees -> "@[Each processor selects the trees it will work on@]"
    | `Skip
    | `StoreTrees
    | `UnionStored
    | `OnEachTree _
    | `ParallelPipeline _
    | `GetStored -> 
            (* These are produced by the analyzer itself, so they can't occur in
            * a script and make no sense by themselves *)
            ""


let colors fo item = 
    let nice () = fo "@{<c:green>"
    and tough () = fo "@{<c:red>"
    and cool () = fo "@{<c:cyan>"
    and ok () = fo "@{<c:white>" in
    match item with
    | Composable -> nice ()
    | InteractivePoint
    | NonComposable -> tough ()
    | Parallelizable -> cool ()
    | Linnearizable 
    | Invariant 
    | ExitPoint -> ok ()

let outputlist fo lst = 
    let lst = List.sort (fun a b -> a - b) lst in
    List.iter (fun x -> fo (string_of_int x); fo ".") lst 

let outputdep fo dep lst = 
    fo " [";
    List.iter (function
        | Channel StandardError -> fo "stderr; "
        | Channel StandardOutput -> fo "stdout; "
        | Channel (File string) -> fo ("file \"" ^ string ^ "\"; ")
        | Channel FutureFiles -> ()
        | Data -> fo "data; "
        | Trees -> fo "trees; "
        | JackBoot -> fo "jackboot; "
        | Bremer -> fo "bremer; "
        | EntryPoint -> fo "entrypoint; ") dep;
    fo "] from ";
    outputlist fo lst

let rec explain_tree ?(deps=true) fo colors tree =
    let deps = false in
    match tree with
    | InteractiveState _ -> ()
    | Parallel x ->
            fo ("@[in parallel:@]@,@[<v 2>@,@[<v>");
            explain_tree fo colors x.todo_p;
            fo "@]@,@[<v 2>@[while keeping the following invariant:@]@,@[<v>";
            explain_tree fo colors x.composer;
            fo "@]@]@]@,";
            explain_tree fo colors x.next
    | Concurrent x ->
            fo (script_to_string x.run);
            fo ("@}@,@[<v 2>@,@[<v>@[I@ will@ calculate@ the@ following@ in@ " ^
            "separate@ processors@ (if available)@]@,@,");
            let processor = ref 1 in
            List.iter (fun tree ->
                fo 
                ("processor group " ^ 
                string_of_int !processor ^ ":@,@[<v 2>@,");
                incr processor;
                explain_tree fo colors tree;
                fo "@]@,@,";) (sort_list x.children);
            fo ("@,@]@[<v>@[Then@ @]");
            List.iter (fun (dep, code, tree) ->
                explain_tree ~deps:true fo colors tree;
                fo "@,@,";) (sort_list2 x.unique);
            fo  "----------@,";
            fo "@]@]@,"
    | Tree x ->
            let filter_unique acc (lst1, lst2) = 
                let lst1 = List.sort compare lst1
                and lst2 = List.sort (fun a b -> a - b) lst2 in
                if List.exists (fun (x, y) -> x = lst1 && y = lst2) acc then acc
                else (lst1, lst2) :: acc
            in
            match x.run, x.children with
            | `StoreTrees, [child]
            | `UnionStored, [child]
            | `OnEachTree _, [child]
            | `GetStored, [child] -> explain_tree fo colors child;
            | _, children ->
                    colors x.exploders;
                    (* outputlist fo x.thread;*)
                    if deps then begin
                        let dep = List.fold_left filter_unique [] !(x.unique) in
                        fo "with {@[<v 2>@,@[<v>";
                        List.iter (fun (dep, code) ->
                            outputdep fo dep code; fo "@,") dep;
                        fo "@]@]}, the following@,";
                    end;
                    fo (script_to_string x.run);
                    fo "@}@,";
                    match children with
                    | [] -> ()
                    | [chld] -> explain_tree fo colors chld;
                    | children ->
                            fo ("@}@,@[<v 2>@[<v>@[I@ will@ " ^
                            "calculate@ the@ following@ in@ " ^
                            "separate@ processors@ (if available)@]@,@,");
                            let processor = ref 1 in
                            List.iter (fun tree ->
                                fo ("processor group " ^ 
                                string_of_int !processor ^ ":@,@[<v 2>@,");
                                incr processor;
                                explain_tree fo colors tree;
                                fo "@]@,@,";) (sort_list children);
                                fo "@,@]@]@,"

let remove_first tree = 
    match tree with
    | Tree x -> x.children
    | tree -> [tree]

let explain_tree filename script =
    let output = Status.user_message (Status.Output (filename, false, [])) in
    let buffer = Buffer.create 1000 in
    let fo = Buffer.add_string buffer in
    let colors = colors fo in
    fo "@[<v>";
    let _ = (* Now we print the explanation *)
        script
        --> maketree 
        --> List.iter 
            (fun t -> explain_tree fo colors t; fo "@,@,") 
    in
    fo "@]%!";
    output (Buffer.contents buffer)

let my_part mine n a = 
    let rest = a - ((a / n) * n) in
    if mine >= rest then 
        a / n
    else 1 + (a / n)

let is_master_only (init : Methods.script) = match init with
    | `Algn_Newkk
    | `Algn_Normal
    | `Barrier 
    | `GatherTrees _
    | `GatherJackknife 
    | `GatherBremer 
    | `SelectYourTrees 
    | `GatherBootstrap 
    | #Methods.tree_handling 
    | #Methods.characters_handling 
    | `Exit 
    | `Interactive
    | `ChangeWDir _ 
    | `Normal_plus_Vitamines
    | `Optimization _
    | `Normal
    | `Exhaustive_Strong
    | `Exhaustive_Weak
    | `Iterative _
    | `ReDiagnose
    | `ReDiagnoseTrees
    | `ClearMemory _ 
    | `Recover 
    | `ClearRecovered 
    | `Wipe 
    | `Load _ 
    | `RootName _
    | `Root _
    | `Entry 
    | `Skip
    | `StoreTrees
    | `UnionStored 
    | `OnEachTree _
    | `ParallelPipeline _ 
    | `GetStored 
    | #Methods.input 
    | #Methods.transform
    | #Methods.build 
    | #Methods.local_optimum 
    | `StandardSearch _
    | #Methods.perturb_method 
    | #Methods.escape_local
    | `Fusing _
    | `Bootstrap _
    | `Jackknife _
    | `Bremer _ 
    | #Methods.runtime_store 
    | `ReadScript _ 
    | `Repeat _ 
    | `TimerInterval _
    | `Parmap _
    | #Methods.taxa_handling -> false
    | `Graph _
    | `Ascii _
    | `InspectFile _
    | `SequenceStats _
    | `Ci _
    | `Ri _
    | `CompareSequences _
    | `FasWinClad _
    | `Nexus _
    | `Model _
    | `LKSites _
    | `Script _
    | `ExplainScript _
    | `PrintWDir
    | `Memory _
    | `KML _
    | `HistorySize _
    | `Redraw
    | `Echo _
    | `Help _
    | `Logfile _
    | `SetSeed _
    | `Alias _
    | `Dataset _
    | `Xslt _
    | `Nodes _
    | `TerminalsFiles _
    | `CrossReferences _
    | `RobinsonFoulds _
    | `GraphicSupports _
    | `Topo_Selection _
    | `Supports _
    | `Consensus _ 
    | `GraphicConsensus _
    | `GraphicDiagnosis _
    | `Diagnosis _
    | `AllRootsCost _
    | `Trees _
    | `Implied_Alignment _
    | `TimeDelta _
    | `KolmoMachine _
    | `TreeCosts _
    | `TreesStats _
    | `SearchStats _
    | `Clades _
    | `Save _
    | `MstR _
    | `DebugData
    | `Plugin _ 
    | `Version -> true

    (** TODO CHECK ALL SCRIPTS FOR LOCAL/REMOTE FILE ISSUES *)
let rec make_remote_files (init : Methods.script) =
    let mr = function
        | `Local x | `Remote x -> `Remote x
    in
    let mrl = List.map mr in
    let handle_subtypes filter items =
        let res = List.map (fun x -> 
            make_remote_files (x :> Methods.script)) items
        in
        List.map filter res
    in
    let handle_simple_input_list files =
        handle_subtypes 
        (function #Methods.simple_input as x -> x | _ -> failwith "Huh?")
        files
    in
    let handle_subtype filter item = 
        let res = make_remote_files (item :> Methods.script) in
        filter res
    in
    let handle_lo item =
        handle_subtype
        (function #Methods.local_optimum as x -> x | _ -> failwith "Huh?")
        item
    and handle_bu item =
        handle_subtype
        (function #Methods.build as x -> x | _ -> failwith "Huh?")
        item
    and recursive = List.map make_remote_files in
    match init with
    | `ParallelPipeline (a, b, c, d) -> `ParallelPipeline (a, recursive b, recursive c, recursive d) 
    | `OnEachTree (a, b) -> `OnEachTree (recursive a, recursive b)
    | `Repeat (a, b) -> `Repeat (a, recursive b)
    | `GatherTrees (a, b) -> `GatherTrees (recursive a, recursive b)
    | `AnalyzeOnlyCharacterFiles (b, files) -> `AnalyzeOnlyCharacterFiles (b, mrl files)
    | `Poyfile files -> `Poyfile (mrl files)
    | `AutoDetect files -> `AutoDetect (mrl files)
    | `PartitionedFile files -> `PartitionedFile (mrl files)
    | `Nucleotides files -> `Nucleotides (mrl files)
    | `Aminoacids (files,opt) -> `Aminoacids (mrl files,opt)
    | `GeneralAlphabetSeq (a, b, c) -> `GeneralAlphabetSeq (mrl a, mr b, c)
    | `Breakinv (a, b, c) -> `Breakinv (mrl a, mr b, c)
    | `Chromosome files -> `Chromosome (mrl files)
    | `Genome files -> `Genome (mrl files)
    | `Prealigned (meth, cost, gap_opening) ->
            let cost = match cost with
                | `Assign_Transformation_Cost_Matrix (file,level) -> 
                    `Assign_Transformation_Cost_Matrix ((mr file),level)
                | `Create_Transformation_Cost_Matrix _ -> cost
            in
            let meth = match meth with
                | `Poyfile files -> `Poyfile (mrl files)
                | `AutoDetect files -> `AutoDetect (mrl files)
                | `Nucleotides files -> `Nucleotides (mrl files)
                | `PartitionedFile files -> `PartitionedFile (mrl files)
                | `Aminoacids (files,opt) -> `Aminoacids (mrl files,opt)
                | `GeneralAlphabetSeq (a, b, c) -> `GeneralAlphabetSeq (mrl a, mr b, c)
                | `Breakinv (a, b, c) -> `Breakinv (mrl a, mr b, c)
                | `Chromosome files -> `Chromosome (mrl files)
                | `Genome files -> `Genome (mrl files)
                | `ComplexTerminals files -> `ComplexTerminals (mrl files)
            in
            `Prealigned (meth, cost, gap_opening)
    | `AnnotatedFiles files ->
            let files = handle_simple_input_list files in
            `AnnotatedFiles files
    | `Assign_Transformation_Cost_Matrix (a, b) ->
            let a = match a with
                | Some (a,l) -> Some ((mr a),l)
                | None -> None
            in
            `Assign_Transformation_Cost_Matrix (a, b)
    | `Assign_Tail_Cost ((`File a), b) ->
            `Assign_Tail_Cost ((`File (mr a)), b)
    | `Assign_Prep_Cost ((`File a), b) ->
            `Assign_Prep_Cost ((`File (mr a)), b)
    | `LocalOptimum (l_opt) ->
        begin match l_opt.Methods.tabu_join with
            | `Partition files ->
                let files = 
                    List.map 
                        (function 
                            | `ConstraintFile file -> 
                                `ConstraintFile (mr file)
                            | x -> x)
                        files 
                in
                let tabu = `Partition files in
                 `LocalOptimum { l_opt with Methods.tabu_join = tabu }
            | _ ->  init
        end
    | `PerturbateNSearch (tl, pm, lo, v, timer) ->
            let tl = 
                handle_subtypes 
                (function #Methods.transform as x -> x | _ -> failwith "Huh?")
                tl
            and lo = handle_lo lo in
            `PerturbateNSearch (tl, pm, lo, v, timer)
    | `Fusing (a, b, c, d, lo, f) ->
            let lo = handle_lo lo in
            `Fusing (a, b, c, d, lo, f)
    | `Bootstrap (a, lo, bu, c) ->
            let lo = handle_lo lo
            and bu = handle_bu bu in
            `Bootstrap (a, lo, bu, c)
    | `Jackknife (a, b, lo, bu, c) ->
            let lo = handle_lo lo
            and bu = handle_bu bu in
            `Jackknife (a, b, lo, bu, c)
    | `Bremer (lo, bu, a, b) ->
            let lo = handle_lo lo
            and bu = handle_bu bu in
            `Bremer (lo, bu, a, b)
    | `SynonymsFile files -> `SynonymsFile (mrl files)
    | `AnalyzeOnlyFiles (a, b) -> `AnalyzeOnlyFiles (a, mrl b)
    | x -> x

let rec cleanup_non_master script =
    List.fold_right (fun x acc ->
        if is_master_only x then acc
        else 
            (* Before prepending, we cleanup the inner script if present *)
            (match x with
            | `Repeat (a, b) ->
                    `Repeat (a, cleanup_non_master b)
            | `OnEachTree (a, b) ->
                    `OnEachTree (a, cleanup_non_master b)
            | `ParallelPipeline (a, b, c, d) ->
                    `ParallelPipeline (a, b, c, cleanup_non_master d)
            | `GatherTrees (a, b) ->
                    `GatherTrees (a, cleanup_non_master b)
            | x -> x) :: acc) script []

let rec parallel_analysis mine n (script : Methods.script list) = 
    (* The script must have been analyzed previously *)
    let res =
        List.fold_left (fun script meth ->
            let meth = 
                match meth with
                | #Methods.input ->
                        [meth; `Barrier]
                | `ParallelPipeline (a, b, c, d) ->
                        let d = parallel_analysis mine n d in
                        let d =
                            match d with
                            | `GetStored :: t ->
                                    [`Barrier; `GetStored; `GatherTrees (c, d)]
                            | _ -> failwith "No GetStored?"
                        in
                        let a = my_part mine n a in
                        [`ParallelPipeline(a, b, c, d)]
                | (`OnEachTree (_, m)) ->
                        let m = `UnionStored :: (m @ [`StoreTrees]) in
                        [`GatherTrees (m, [`GetStored]); `Barrier; meth; 
                        `SelectYourTrees] 
                | `Jackknife (a, it, b, c, d) ->
                        let it = my_part mine n it in
                        [`GatherJackknife; `Jackknife (a, it, b, c, d)]
                | `Bootstrap (it, b, c, d) ->
                        let it = my_part mine n it in
                        [`GatherBootstrap; `Bootstrap (it, b, c, d)]
                | `Bremer (a, b, c, d) ->
                        [`GatherBremer; `Bremer (a, b, mine, n)]
                | `Fusing (iterations, max_trees, a, b, c, d) ->
                        let merger = 
                            match max_trees with
                            | None -> [`UnionStored; `StoreTrees]
                            | Some x ->  [`UnionStored; `BestN max_trees; 
                            `StoreTrees]
                        in
                        let meth =
                            match iterations with
                            | None -> meth
                            | Some x -> 
                                    `Fusing 
                                    (Some (my_part mine n x), max_trees, 
                                    a, b, c, d)
                        in
                        [`GatherTrees (merger, [`GetStored]); meth; 
                        `SelectYourTrees ]
                 | meth -> [meth]
            in
            meth @ script) [] script
    in
    let res = List.rev res in
    let res = List.map make_remote_files res in
    if mine <> 0 then cleanup_non_master res
    else res

