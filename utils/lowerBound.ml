type label = String of (string * int) | Integer of int

let get_code = function String (_, mine) | Integer mine -> mine

let assign_labels tree =
    let cnt = ref (-1) in
    let rec assign_labels tree =
        match tree with
        | Parser.Tree.Leaf x -> 
                incr cnt;
                Parser.Tree.Leaf (String (x, !cnt))
        | Parser.Tree.Node (chld, _) ->
                let chld = List.map assign_labels chld in
                incr cnt;
                Parser.Tree.Node (chld, Integer !cnt)
    in
    assign_labels tree

let merge_paths left right code =
    let left_pass_back = Sexpr.map (Sexpr.union code) left
    and right_pass_back = Sexpr.map (Sexpr.union code) right in
    let all_connected = Sexpr.all_to_all Sexpr.union left_pass_back right in
    let all_to_pass_back = Sexpr.union left_pass_back right_pass_back in
    all_connected, all_to_pass_back

let produce_paths tree =
    let rec generate_all_paths tree =
        match tree with
        | Parser.Tree.Leaf x ->
                let mine = get_code x in 
                let res = Sexpr.singleton (Sexpr.singleton mine) in
                `Empty, res
        | Parser.Tree.Node ([ch1; ch2], x) ->
                let mine = 
                    let mine = get_code x in
                    Sexpr.singleton mine
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
    match tree with
    | Parser.Tree.Leaf _ -> []
    | Parser.Tree.Node ([ch1; ch2], _) ->
            let ch1_ac, ch1_apb = generate_all_paths ch1
            and ch2_ac, ch2_apb = generate_all_paths ch2 in
            let union = Sexpr.all_to_all Sexpr.union ch1_apb ch2_apb in
            Sexpr.to_list (Sexpr.union ch1_ac (Sexpr.union union ch2_ac))
    | _ -> assert false

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

let produce_LP_instance gap_opening tcm file fasta synonyms = 
    match Parser.Tree.of_file (`Local file) with
    | [[tree]] ->
            let costs, leaves, tree, rtree =  (* Compute all the costs *)
                let seqs = 
                    let seqs = Scripting.DNA.Fasta.of_file false fasta in
                    List.map (fun (a, b) ->
                        (if Hashtbl.mem synonyms a then Hashtbl.find synonyms a 
                        else a), b) seqs
                in
                let rtree, seqs = 
                    let res = ref All_sets.StringMap.empty in
                    let exists = List.map fst seqs in
                    let cnt = ref (-1) in
                    let res2 =
                    Parser.Tree.cleanup (fun name ->
                        let () =
                            try
                                let seq = List.assoc name seqs in
                                res := All_sets.StringMap.add name seq !res;
                            with
                            | Not_found -> ()
                        in
                        let res = (name = "") || not (List.exists (( = ) name) exists) in
                        if res then res
                        else 
                            let () = incr cnt in
                            res
                        ) tree
                    in
                    match res2 with
                    | None -> assert false
                    | Some res2 ->
                            res2, All_sets.StringMap.fold (fun a b acc -> (a, b) :: acc)
                            !res []
                in
                let tree as rtree = assign_labels rtree in
                let leaves, tree = list_of_leaves tree, tree in
                let tcm = 
                    let tcm = Scripting.DNA.CM.of_file tcm in
                    Cost_matrix.Two_D.set_affine tcm (Cost_matrix.Affine
                    gap_opening);
                    tcm
                in
                let tbl = Hashtbl.create 1667 in
                let cleanup ((a, _), (b, _), c, _) = 
                    Printf.printf "%s\t%s\t%d\n%!" a b c;
                    Hashtbl.add tbl (List.assoc a leaves, List.assoc b leaves) c
                in
                List.iter (List.iter cleanup) 
                (Scripting.DNA.Align.algn_all seqs tcm);
                tbl, 
                List.fold_left (fun acc x -> All_sets.Integers.add (snd x) acc)
                All_sets.Integers.empty leaves, tree, rtree
            in
            let total_variables = 
                ((All_sets.Integers.cardinal leaves) * 2) - 3
            in
            let minimization_problem = Array.make total_variables 1. in
            let constraints = 
                let paths = produce_paths tree in
                List.map (make_constraint leaves (fun a b -> Hashtbl.find costs (a,
                b))) paths
            in
            let a, b = 
                List.fold_left (fun (stacc, cacc) (c, sexp) ->
                    let cacc = (float_of_int c, float_of_int max_int) :: cacc in
                    let arr = Array.make total_variables 0. in
                    Sexpr.leaf_iter (fun code ->
                        if code = total_variables then ()
                        else arr.(code) <- 1.) sexp;
                    (arr :: stacc), cacc) ([], []) constraints
            in
            minimization_problem, Array.of_list a, Array.of_list b,
            Array.make total_variables (0., float_of_int max_int), rtree
    | _ -> assert false

let process_problem gap_opening tcm treefile synfile fasta = 
    let synonyms = 
        let ch = open_in synfile in
        let res = Parser.Dictionary.of_channel ch in
        close_in ch;
        res
    in
    let gap_opening = int_of_string gap_opening in
    let a, b, c, d, tree = produce_LP_instance gap_opening tcm treefile fasta synonyms in
    let lp = Glpk.make_problem Glpk.Minimize a b c d in
    Glpk.scale_problem lp;
    Glpk.use_presolver lp true;
    Glpk.simplex lp;
    let res = Glpk.get_obj_val lp in
    let variables = Glpk.get_col_primals lp in
    let tree = Parser.Tree.map (fun x ->
        let name, x =
            match x with
            | String x -> x
            | Integer x -> "", x
        in
        try name ^ ":" ^ string_of_float variables.(x) with
        | _ -> name) tree in
    res, tree

let () = 
    match Array.to_list Sys.argv with
    | _ :: a :: b :: c :: d :: tl -> 
            let res = 
                List.fold_left (fun acc file ->
                    let res, tree = process_problem a b c d file in
                    Printf.printf "Solution to this problem is : %f\n%!" res;
                    AsciiTree.draw_parenthesis stdout tree;
                    acc +. res) 0. tl
            in
            print_float res;
            print_newline ()
    | name :: _ ->
            Printf.fprintf stderr
            "Usage: %s gap_opening tcm_file tree_file synonyms [fasta_file] ...\n%!" name;
            exit 1
    | [] -> assert false
