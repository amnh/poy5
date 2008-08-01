type label = String of (string * int) | Integer of int

let assign_labels tree =
    let cnt = ref (-1) in
    let rec assign_labels tree =
        incr cnt;
        let mine = !cnt in
        match tree with
        | Parser.Tree.Leaf x -> 
                Parser.Tree.Leaf (String (x, mine))
        | Parser.Tree.Node (chld, _) ->
                let chld = List.map assign_labels chld in
                Parser.Tree.Node (chld, Integer mine)
    in
    assign_labels tree

let merge_paths left right code =
    let left_pass_back = Sexpr.map (Sexpr.union code) left
    and right_pass_back = Sexpr.map (Sexpr.union code) right in
    let all_connected = Sexpr.all_to_all Sexpr.union left_pass_back right in
    let all_to_pass_back = Sexpr.union left_pass_back right_pass_back in
    all_connected, all_to_pass_back

let get_code = function String (_, mine) | Integer mine -> mine

let produce_paths tree =
    let rec generate_all_paths tree =
        match tree with
        | Parser.Tree.Leaf x ->
                let res = Sexpr.singleton (Sexpr.singleton (get_code x)) in
                `Empty, res
        | Parser.Tree.Node ([ch1; ch2], x) ->
                let mine = 
                    let mine = get_code x in
                    if mine <> 0 then Sexpr.singleton mine
                    else `Empty
                in
                let ch1_ac, ch1_apb = generate_all_paths ch1 in
                let ch2_ac, ch2_apb = generate_all_paths ch2 in
                let all_connected = Sexpr.union ch1_ac ch2_ac in
                let new_all_connected, all_to_pass_back = 
                    merge_paths ch1_apb ch2_apb mine
                in
                Sexpr.union all_connected new_all_connected, all_to_pass_back
        | Parser.Tree.Node _ -> assert false
    in
    let paths, _ = generate_all_paths tree in
    Sexpr.to_list paths 

let list_of_leaves tree = 
    let rec list_of_leaves acc tree =
        match tree with
        | Parser.Tree.Leaf x -> x :: acc
        | Parser.Tree.Node (lst, _) ->
                List.fold_left list_of_leaves acc lst 
    in
    List.map (function String x -> x | Integer _ -> assert false)
    (list_of_leaves [] tree)

let make_constraint leaves distance set =
    let leaves =
        Sexpr.to_list (Sexpr.filter (fun x -> All_sets.Integers.mem x leaves)
        set)
    in
    match leaves with
    | [x; y] ->
            let di = distance x y in
            (di, set)
    | _ -> assert false

let produce_LP_instance tcm file fasta = 
    match Parser.Tree.of_file (`Local file) with
    | [[tree]] ->
            let tree = assign_labels tree in
            let leaves = list_of_leaves tree in
            let costs, leaves =  (* Compute all the costs *)
                let seqs = Scripting.DNA.Fasta.of_file false fasta in
                let tcm = Scripting.DNA.CM.of_file tcm in
                let tbl = Hashtbl.create 1667 in
                let cleanup ((a, _), (b, _), c, _) = 
                    Hashtbl.add tbl (List.assoc a leaves, List.assoc b leaves) c
                in
                List.iter (List.iter cleanup) 
                (Scripting.DNA.Align.algn_all seqs tcm);
                tbl, 
                List.fold_left (fun acc x -> All_sets.Integers.add (snd x) acc)
                All_sets.Integers.empty leaves
            in
            let total_variables = 
                ((All_sets.Integers.cardinal leaves) * 2) - 2
            in
            let minimization_problem = Array.make total_variables 1. in
            let constraints = 
                let paths = produce_paths tree in
                List.map (make_constraint leaves (fun a b -> Hashtbl.find costs (a,
                b))) paths
            in
            let a, b = 
                List.fold_left (fun (stacc, cacc) (c, sexp) ->
                    let cacc = (float_of_int c, infinity) :: cacc in
                    let arr = Array.make total_variables 0. in
                    Sexpr.leaf_iter (fun code ->
                        arr.(code - 1) <- 1.) sexp;
                    (arr :: stacc), cacc) ([], []) constraints
            in
            minimization_problem, Array.of_list a, Array.of_list b,
            Array.make total_variables (0., infinity)
    | _ -> assert false

let process_problem tcm treefile fasta = 
    let a, b, c, d = produce_LP_instance tcm treefile fasta in
    let lp = Glpk.make_problem Glpk.Minimize a b c d in
    Glpk.scale_problem lp;
    Glpk.use_presolver lp true;
    Glpk.simplex lp;
    Printf.printf "Z:%g\n%!" (Glpk.get_obj_val lp)

let () = process_problem Sys.argv.(1) Sys.argv.(2) Sys.argv.(3)
