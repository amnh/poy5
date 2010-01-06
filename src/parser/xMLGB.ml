let convert_to_fasta file =
    let ch, file = FileStream.channel_n_filename file and
    accession_str = ref "" and
    found_gi_id = ref 0 and
    outName, chout = Filename.open_temp_file "fasta" ".tmp" in
    try
        while true do
            let line = input_line ch in
            if Str.string_match (Str.regexp 
            " *\\(<Textseq-id_accession>+\\)\\([a-zA-Z0-9 ,]+\\)")
            line 0 then
                accession_str :=  (Str.matched_group 2 line) ^ "." 
            else if Str.string_match (Str.regexp 
            " *\\(<Textseq-id_version>+\\)\\([0-9]+\\)")
            line 0 then
                accession_str := !accession_str ^
                 (Str.matched_group 2 line) ^ "|" 
            else if (Str.string_match (Str.regexp 
            " *\\(<Seq-id_gi>+\\)\\([0-9]+\\)")
            line 0) && !found_gi_id == 0 then
                begin
                    output_string chout 
                    ( ">gi|" ^ (Str.matched_group 2 line) ^ "|" ^ 
                    !accession_str);
                    found_gi_id := 1;
                end
            else if Str.string_match (Str.regexp 
            " *\\(<Seqdesc_title>+\\)\\([a-zA-Z0-9 ,]+\\)") line 0 then
                output_string chout (Str.matched_group 2 line) 
            else if Str.string_match (Str.regexp 
            " *\\(<IUPACna>+\\)\\([a-zA-Z]+\\)") line 0 then
                begin
                    output_string chout 
                    ( "\n" ^ (Str.matched_group 2 line) ^ "\n\n");
                end
            else if Str.string_match (Str.regexp 
            " *<Seq-entry>") line 0 then
                    found_gi_id := 0;
        done;
        open_in outName;
    with
    | End_of_file -> (open_in outName) 

