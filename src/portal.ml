(* Parameters for the search *)
let input_files = ref []
let terminals = ref None
let synonyms = ref None
let out_group_file = ref None
let randomseed = ref None
let constraintfile = ref None
let substitutions = ref 1
let support = ref ""
let indels = ref 2
let gapopening = ref 0
let minutes = ref (71. *. 60.)
let max_memory = ref 2.
let output_trees = "trees.result"
let tmp_trees = "bremertrees"
let output_alignments = "alignments.result"
let output_analyzed_data = "data.result"
let output_support = "support.result"
let output_graph_support = "graphsupport.result"
let pseudoreplicates = ref None
let remove = ref None

let () = 
    let assgn x y = x := y in
    let assgn_opt x y = x := Some y in
    let params = 
        Arg.align
        [
        ("-constraint", Arg.String (assgn_opt constraintfile), 
        "FILE File to be used as constraint for the tree search.");
        ("-gapopening", Arg.Int (assgn gapopening), "INT Gap opening parameter. Default: 0.");
        ("-indels", Arg.Int (assgn indels), "INT Indel parameter. Default: 1.");
        ("-input", Arg.String (fun x -> input_files := x :: !input_files), "FILE An input file, in any format");
        ("-minutes", Arg.Float (assgn minutes), 
        "FLOAT Total number of minutes spent in the search. Default: 71 hours.");
        ("-outgroup", Arg.String (assgn_opt out_group_file), 
        "FILE Name of the terminal that should be used as outgroup");
        ("-pseudoreplicates", Arg.Int (assgn_opt pseudoreplicates), "INT Number of pseudoreplicates to be used");
        ("-randomseed", Arg.Int (assgn_opt randomseed), "INT Random number generator seed to be used.");
        ("-remove", Arg.Float (assgn_opt remove), "FLOAT Percentage of characters to be removed");
        ("-substitutions", Arg.Int (assgn substitutions), "INT Substitution parameter. By Default: 2.");
        ("-support", Arg.String (assgn support), 
        "[jackknife|bootstrap|bremer] Support values to be computed. Valid options are jackknife, bootstrap, and bremer.");
        ("-synonyms", Arg.String (assgn_opt synonyms), "FILE File containing the list of synonyms to be used.");
        ("-terminals", Arg.String (assgn_opt terminals), "FILE Terminals file to be used.");
    ] in
    let usage = "portal [OPTIONS]" in
    try
        Arg.parse_argv Phylo.args params 
        (fun x -> raise (Arg.Bad ("Illegal anonymous argument: " ^ x ^". Maybe you want to use the -input argument?"))) usage
    with
    | Arg.Help _ -> 
            Arg.usage params usage; 
            exit 0
    | Arg.Bad str ->
            prerr_string ("Bad argument: " ^ str);
            Arg.usage params usage;
            exit 1

let () = 
    match !input_files with
    | [] -> exit 0
    | _ ->
            let () =
                match !randomseed with
                | None -> ()
                | Some x -> POY set (seed:[x])
            in
            let () =
                match !synonyms with
                | None -> ()
                | Some x -> POY rename (terminals, [x])
            in
            List.iter (fun x -> POY read ([x])) !input_files;
            let () =
                POY 
                transform (tcm:([!substitutions], [!indels]), 
                gap_opening:[!gapopening])
            in
            let () =
                match !out_group_file with
                | None -> ()
                | Some x -> (POY set (root:[x]))
            in
            let swap = 
                match !constraintfile with
                | None -> (APOY swap ())
                | Some x -> (APOY swap (constraint_p:(file:[x])))
            in
            let trees = Phylo.Runtime.trees () in
            let () =
                match !support with
                | "" -> 
                        let searchargs =
                            let args = [APOY hits:100; 
                                APOY max_time:0:0:[!minutes]; 
                                APOY memory:gb:[!max_memory]] 
                            in
                            match !constraintfile with
                            | None -> args
                            | Some x -> (APOY constraint_s:[x]) :: args
                        in
                        (* We want to run a standard search *)
                        POY
                            search {searchargs}
                            select ()
                            report ([output_trees], trees:(total))
                            report ([output_alignments], ia)
                | "bremer" when [] <> Phylo.Runtime.trees () ->
                        POY 
                            search (visited:[tmp_trees], hits:10,
                            max_time:0:0:[(!minutes /. 2.)])
                            select (best:0);
                        Phylo.Runtime.set_trees trees;
                        POY
                            report ([output_support], supports:bremer:[tmp_trees])
                            report ([output_graph_support],
                                graphsupports:bremer:[tmp_trees])
                | "jackknife" when [] <> Phylo.Runtime.trees () ->
                        let jackknife = 
                            match !remove, !pseudoreplicates with
                            | None, None -> APOY jackknife
                            | Some v, None -> APOY jackknife:(remove:[v])
                            | None, Some x -> APOY jackknife:(resample:[x])
                            | Some v, Some x -> 
                                    APOY jackknife:(resample:[x], remove:[v])
                        in
                        POY
                            calculate_support ({swap}, {jackknife})
                            report ([output_support], supports:jackknife)
                            report ([output_graph_support], graphsupports:jackknife)
                | "bootstrap" when [] <> Phylo.Runtime.trees () -> 
                        let bootstrap =
                            match !pseudoreplicates with
                            | None -> APOY bootstrap
                            | Some x -> APOY bootstrap:[x]
                        in
                        POY
                            calculate_support ({swap}, {bootstrap})
                            report ([output_support], supports:bootstrap)
                            report ([output_graph_support],
                                graphsupports:bootstrap)
                | "bremer" | "bootstrap" | "jackknife" -> 
                        prerr_string 
                        "In the portal application, computing support values \
                        need trees in one of your input files to report the \
                        values on their branches";
                        prerr_newline ();
                        exit 1
                | x -> 
                        prerr_string "Unknown support value :";
                        prerr_string x;
                        prerr_newline ();
                        exit 1
            in
            POY
                report ([output_analyzed_data], data, cross_references)
                exit ()
