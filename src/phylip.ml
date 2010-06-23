let convert_ file =
    let taxa_char_lengths = Str.regexp "\\(^[0-9]+\\)[\t ]+\\([0-9]+\\)"
    and name_sequence = Str.regexp " *\\([a-zA-Z0-9]+\\)[\t ]+\\([a-zA-Z?-]+\\)+"
    and continued_seq = Str.regexp " *\\([a-ZA-Z?-]+\\)"
    and reset_order = Str.regexp "[\t ]*" in

    let ch, file = FileStream.channel_n_filename file
    and index = ref 0 in

    let ntaxa,nchar = 
        let line = input_line ch in
        if Str.string_match taxa_char_lengths line 0 then 
                (Str.matched_group 1 line), (Str.matched_group 2 line)
        else
            failwith "Error parsing Phylip format; number of taxa and characters."
    in
    let num_taxa_int = int_of_string ntaxa
    and num_char_int = int_of_string nchar in

    let sequence_array = Array.init num_taxa_int (fun _ -> Buffer.create num_char_int)
    and taxa_array = Array.make num_taxa_int ""
    and line_number = ref 0 in
    try
        while true do
            let line = input_line ch in
            let () = incr line_number in
            if Str.string_match name_sequence line 0 then
                begin
                    taxa_array.(!index) <- (Str.matched_group 1 line);
                    let added_data = (Str.matched_group 2 line) in
                    Buffer.add_string sequence_array.(!index) added_data;
                    index := (!index + 1) mod num_taxa_int;
                end
            else if Str.string_match continued_seq line 0 then
                begin
                    let added_data = (Str.matched_group 1 line) in
                    Buffer.add_string sequence_array.(!index) added_data;
                    index := (!index + 1) mod num_taxa_int;
                end
            else if Str.string_match reset_order line 0 then
                index := 0
            else
                Printf.kprintf failwith "I Cannot parse the Phylip file line %d." !line_number
        done;
        raise End_of_file
    with | End_of_file -> taxa_array, sequence_array

let match_alphabets observed = 
    let is_dna =
        try All_sets.Strings.fold
                (fun e d -> let _ = Alphabet.match_base e Alphabet.dna in d)
                observed true
        with Alphabet.Illegal_Character _ -> false
    and is_prot = 
        try All_sets.Strings.fold
                (fun e p -> let _ = Alphabet.match_base e Alphabet.aminoacids in p)
                observed true
        with Alphabet.Illegal_Character _ -> false
    and is_nuc = 
        try All_sets.Strings.fold
                (fun e n -> let _ = Alphabet.match_base e Alphabet.nucleotides in n)
                observed true
        with Alphabet.Illegal_Character _ -> false
    in
    match is_dna,is_nuc,is_prot with
        | true,_,_ -> Alphabet.dna
        | _,true,_ -> Alphabet.nucleotides
        | _,_,true -> Alphabet.aminoacids
        | false,false,false ->
                failwith "I cannot detect the Alphabet in the Phylip file"

let matrix_iter f matrix =
    let row = Array.length matrix in 
    assert(row > 0);
    let col = Array.length matrix.(0) in
    for i = 0 to row-1 do for j = 0 to col-1 do
        matrix.(i).(j) <- f i j
    done done

let list_of_packed d =
    let rec loop_ c i d = match d land 1 with
        | 0 when d = 0 -> c
        | 0  -> loop_ c (i+1) (d lsr 1)
        | 1  -> loop_ (i::c) (i+1) (d lsr 1)
        | _  -> failwith "Phylip.list_of_packed"
    in loop_ [] 0 d
let (-->) a b = b a

let of_file (file : FileStream.f) = 
    let tax_array, seq_array = convert_ file
    and file = match file with | `Local x | `Remote x -> x in

    let final_taxa_array = Array.map (fun x -> Some x) tax_array
    and nchars = Buffer.length seq_array.(0) 
    and ntaxas = Array.length tax_array in
    let alphabet =
        let observed = 
            Array.fold_left
                (fun acc x ->
                    assert( nchars = (Buffer.length x));
                    let acc2 = ref acc in
                    for i=0 to nchars-1 do
                        let value = String.make 1 (Buffer.nth x i) in
                        if All_sets.Strings.mem value !acc2 then ()
                        else acc2 := All_sets.Strings.add value !acc2
                    done;
                    !acc2)
                (All_sets.Strings.empty)
                seq_array
        in
        match_alphabets observed
    in
    let final_seq_matrix =
        let seq_matrix = Array.make_matrix ntaxas nchars None in
        matrix_iter
            (fun i j ->
                let values = 
                    Buffer.nth (seq_array.(i)) j
                        --> String.make 1 
                        --> (fun x -> Alphabet.match_base x alphabet)
                        --> list_of_packed
                in
                Some (`List values) )
            seq_matrix;
        seq_matrix
    and get_observed seq i = 
        let rec add_elms full = function
            | x::xs when List.exists (fun y -> x = y) full -> add_elms full xs
            | x::xs -> add_elms (x::full) xs
            | [] -> full
        in
        Array.fold_left
            (fun acc x -> 
                let values =
                    match x.(i) with 
                    | Some (`List x) -> x
                    | _ -> failwith "impossible"
                in
                add_elms acc values)
            [] seq
    in
    Printf.printf "Making array of %d characters\n%!" nchars;
    let final_chars_array = 
        Array.init nchars 
            (fun i-> {Nexus.File.st_filesource = file;
                        st_name = file ^ ":" ^ (string_of_int i);
                        st_alph = alphabet;
                        st_observed = get_observed final_seq_matrix i;
                        st_labels = [];
                        st_weight = 1.0;
                        st_type = Nexus.File.STUnordered;
                        st_equivalents = [(Alphabet.gap_repr,[])];
                        st_missing = "?";
                        st_matchstate = None;
                        st_gap = Alphabet.gap_repr;
                        st_eliminate = false;
                        st_case = true;
                        st_used_observed = None;
                        st_observed_used = None; })
    in
    { Nexus.File.empty_parsed () with
      Nexus.File.taxa = final_taxa_array;
               matrix = final_seq_matrix; 
           characters = final_chars_array; }
