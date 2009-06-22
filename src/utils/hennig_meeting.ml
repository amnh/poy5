let make_parameter_list min max step =
    let rec aux acc cur =
        if cur < min then acc
        else aux (cur :: acc) (cur -. step)
    in
    aux [] max


let parameters = make_parameter_list 0.01 0.3 0.02

let indel_subst = 
    [0.1, 0.1, 0.8;
    0.2, 0.2, 0.6;
    0.3, 0.3, 0.4;
    0.33, 0.33, 0.33;
    0.4, 0.4, 0.2]

let branches = [ "branch_0.05"; "branch_0.1"; "branch_0.2"; "branch_0.3"]
let max_lengths = 
    [ "max_len_1"; "max_len_2";"max_len_4"; "max_len_5"; 
    "max_len_10"; "max_len_15"]
let taxa = ["50_taxa"]
let bases = ["300_bases"]

let process_parameter parameter = 
    List.iter (fun (insertion, deletion, substitution) ->
        POY transform (kolmogorov:(event:[parameter], indelsub:([insertion],
            [deletion], [substitution])));
        match Phylo.Runtime.min_cost () with
        | None -> assert false
        | Some cost -> Printf.printf "%f\t%!" cost) 
indel_subst


let process_tree dir file_number =
    Printf.printf "%s%d\t" dir file_number;
    let file_number = string_of_int file_number in
    let fasta_file = "simmulation_" ^ file_number ^ ".fasta" in
    let tree_file = "simmulation_" ^ file_number ^ ".tree.poy" in
    POY wipe ()
        read ([fasta_file], [tree_file]);
    List.iter process_parameter parameters;
    Printf.printf "\n%!"

let execute_in_directory prev_string function_to_continue directories =
    List.iter (fun directory ->
        let dir = prev_string ^ directory ^ "\t" in
        POY cd ([directory]);
        function_to_continue dir;
        POY cd ("../")) directories

let () = 
    for i = 5 to 100 do
        execute_in_directory "" (fun dir ->
            execute_in_directory dir (fun dir ->
                execute_in_directory dir (fun dir ->
                    execute_in_directory dir (fun dir ->
                        process_tree dir i;
                        ) bases
                    ) taxa
                ) max_lengths
            ) branches;
    done;
