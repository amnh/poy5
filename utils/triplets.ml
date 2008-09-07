(* A program to compute a set of triplets distances (sequences). *)
let arg = Mpi.init Sys.argv

let make_triplets arr =
    let res = ref [] in
    let len = Array.length arr in
    for i = 0 to len - 3 do
        res := (arr.(i), arr.(i + 1), arr.(i + 2)) :: !res;
    done;
    !res

let all_triplets arr = 
    let res = ref [] in
    let len = Array.length arr in
    for i = 0 to len - 3 do
        for j = i + 1 to len - 2 do
            for k = j + 1 to len - 1 do
                res := (arr.(i), arr.(j), arr.(k)) :: !res;
            done;
        done;
    done;
    !res

let compute_triplet toseq tcm tcm3 ((a, b, c) as triple) =
    let a_s = toseq a
    and b_s = toseq b 
    and c_s = toseq c in
    (triple, 
    (snd (Sequence.Align.align_3_powell_inter a_s b_s c_s tcm tcm3)))

let rank = Mpi.comm_rank Mpi.comm_world

let compute_triplets file sub indel affine =
    let seqs = 
        let seqs = 
            if rank > 0 then []
            else
                Scripting.DNA.Fasta.of_file false file in
        Mpi.broadcast seqs 0 Mpi.comm_world 
    in
    let index = Hashtbl.create 1667 in
    List.iter (fun (taxon, seq) -> Hashtbl.add index taxon seq) seqs;
    let tcm = Scripting.DNA.CM.of_sub_indel_affine sub indel affine in
    let tcm3 = Cost_matrix.Three_D.of_two_dim tcm in
    let terminals_arr = Array.of_list (List.map fst seqs) in
    let compute_triplet = compute_triplet (Hashtbl.find index) in
    Array_ops.randomize terminals_arr;
    let triplets = make_triplets terminals_arr in
    let triplets = Array.of_list triplets in
    let triplet = Mpi.scatter triplets 0 Mpi.comm_world in
    let triplet = compute_triplet tcm tcm3 triplet in
    let triplets = Mpi.gather triplet 0 Mpi.comm_world in
    if 0 = rank then
        Array.iter (fun ((a, b, c), dist) ->
            Printf.printf "%s\t%s\t%s\t%d\n%!" a b c dist)
        triplets
    else ()

let gap_opening = ref 0
let indel = ref 1
let substitution = ref 1
let input_file = ref ""

let assign x y = x := y

let parse_list = [
    ("-gap-opening", (Arg.Int (assign gap_opening)), "Gap opening");
    ("-indel", (Arg.Int (assign indel)), "Indel");
    ("-substitution", (Arg.Int (assign substitution)), "Substitution");
]

let anon_fun x = assign input_file x

let () = Arg.parse_argv arg  parse_list anon_fun "triplets [OPTIONS] fasta_file.fas"

let () = 
    match !input_file with
    | "" -> exit 1
    | _ -> compute_triplets !input_file !substitution !indel !gap_opening
