let lexer = Genlex.make_lexer [] 

let rec process_file = parser
    | [< x = vertex; children = children_parser; 
        'Genlex.String dna; files = process_file >] -> 
            (x, children, dna) :: files 
    | [< >] -> []
and vertex = parser
    | [< 'Genlex.Ident name >] -> name
    | [< 'Genlex.Int int >] -> string_of_int int
and children_parser = parser
    | [< first_child = vertex; second_child = vertex >] ->
            Some (first_child, second_child)
    | [< >] -> None

let process_file filename =
    let ch = open_in filename in
    let stream = lexer (Stream.of_channel ch) in
    let contents = process_file stream in
    let () = close_in ch in
    contents

let kuse_affine = ref false
let event_probability = ref 0.5
let indel_probability = ref 0.66
let extension_probability = ref 0.5

let log2 x = 
    ((log x) /. (log 2.))

let compute_score max_sequence_length (a, b, _, _) =
    let a = Scripting.DNA.Seq.to_string a 
    and b = Scripting.DNA.Seq.to_string b in
    assert ((String.length a) = (String.length b));
    let previous_was_gap = ref None in
    let len = String.length a in 
    let score = ref 0. in
    let update_score x = score := x +. !score in
    for i = 0 to len - 1 do
        if a.[i] = b.[i] then begin
            let () =
                match !previous_was_gap with
                | Some _ -> 
                        update_score (-. (log2 (1. -. !extension_probability)));
                        update_score (log2 !extension_probability);
                | None -> ()
            in
            update_score (-. (log2 (1. -. !event_probability)));
            previous_was_gap := None;
        end else if (a.[i] = '-') then begin
            score := 
                !score +.
                    (if not !kuse_affine then
                        (-. (log2 !event_probability)) +.
                            (-. (log2 (!indel_probability /. 2.)))
                    else
                        (if (Some "a") = !previous_was_gap then 
                            (-. (log2 !extension_probability))
                        else 
                            (if (Some "b") = !previous_was_gap then
                                (-. (log2 (1. -. !extension_probability)))
                                -.
                                (-. (log2 !extension_probability))
                            else 0.) +.
                                (-. (log2 !event_probability)) +.
                                (-. (log2 (!indel_probability /. 2.)))));
            if !kuse_affine then previous_was_gap := Some "a";
        end else if (b.[i] = '-') then begin
            score := 
                !score +. 2. +.
                    (if not !kuse_affine then
                        (-. (log2 !event_probability)) +.
                            (-. (log2 (!indel_probability /. 2.)))
                    else
                        (if (Some "b") = !previous_was_gap then 
                            (-. (log2 !extension_probability))
                        else 
                            (if (Some "a") = !previous_was_gap then
                                (-. (log2 (1. -. !extension_probability)))
                                -.
                                (-. (log2 !extension_probability))
                            else 0.) +.
                                (-. (log2 !event_probability)) +.
                                (-. (log2 (!indel_probability /. 2.)))));
        end else begin
            score := 
                (if ((Some "a") = !previous_was_gap) || ((Some "b") =
                    !previous_was_gap) then
                    (-. (log2 (1. -. !extension_probability))) -.
                    (-. (log2 !extension_probability))
                else 0.) +.
                    !score +. (-. (log2 !event_probability)) +. 2. +.
                        (-. (log2 (1. -. !indel_probability)));
            previous_was_gap := None;
        end;
    done;
    update_score
        (if ((Some "a") = !previous_was_gap) || ((Some "b") =
            !previous_was_gap) then
                (-. (log2 (1. -. !extension_probability)))
        else 0.);
    !score

let root_complexity sequence =
    let len = Sequence.length sequence in
    let coded = float_of_int (len * 2) in
    let lencode = 2. *. (log2 (float_of_int len)) in
    coded +. lencode

let compute_score filename cm =
    let file_contents = process_file filename in
    let table = Hashtbl.create 97 in
    let max_sequence_length = 
        List.fold_left (fun acc (_, _, string) -> 
            acc + (String.length string)) 0 file_contents
    in
    List.iter (fun (name, _, sequence) -> 
        Hashtbl.add table name (Scripting.DNA.Seq.of_string sequence))
    file_contents;
    let root_complexity = root_complexity (Hashtbl.find table "root") in
    List.fold_left (fun acc (name, children, _) ->
        match children with
        | None -> acc
        | Some (child1, child2) ->
                let vertex = Hashtbl.find table name 
                and child1 = Hashtbl.find table child1
                and child2 = Hashtbl.find table child2 in
                let alignment = Scripting.DNA.Align.algn vertex child1 cm in
                let score1 = compute_score max_sequence_length alignment in
                let alignment = Scripting.DNA.Align.algn vertex child2 cm in
                let score2 = compute_score max_sequence_length alignment in
                score1 +. score2 +. acc) root_complexity file_contents

let substitution = ref 1
let indel = ref 1
let gap_opening = ref 0 
let filename = ref ""

let set_affine x = 
    kuse_affine := true;
    extension_probability := x 

let parse_list = [
    ("-gap-opening", (Arg.Set_int gap_opening), "Gap opening. Default: 0");
    ("-indel", (Arg.Set_int indel), "Indel. Default: 1");
    ("-substitution", (Arg.Set_int substitution), "Substitution. Default: 1");
    ("-file", (Arg.Set_string filename), "File to process.");
    ("-kuse_affine", (Arg.Float set_affine), 
        "Use affine gap cost, with the given extension probability. Default: \
        no affine.");
    ("-event_probability", (Arg.Set_float event_probability), 
        "Probability of transformation. Default: 0.5");
    ("-indel_probability", (Arg.Set_float indel_probability),
        "Probability of an indel (as opposed to a substitution). Default: \
            0.66");
]

let process_file () =
    Arg.parse parse_list (fun str -> failwith ("Illegal argument: " ^ str)) 
        "kolmo_diagnoser [OPTIONS] -file FILENAME";
    match !filename with
    | "" -> ()
    | filename ->
        let cm =
            Scripting.DNA.CM.of_sub_indel_affine !substitution !indel 
            !gap_opening
        in
        print_float (compute_score filename cm)

let () = process_file ()
