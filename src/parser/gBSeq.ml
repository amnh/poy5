let convert_to_fasta file =
    let ch, file = FileStream.channel_n_filename file and
    definition_str = ref "" and
    accession_str = ref "" and
    outName, chout = Filename.open_temp_file "fasta" ".tmp" in
    try
        while true do
            let line = input_line ch in
            if Str.string_match (Str.regexp 
            " *\\(<GBSeq_accession-version>+\\)\\([a-zA-Z0-9. ,]+\\)") 
            line 0 then
                accession_str := 
                ((Str.matched_group 2 line) ^ "|"  ^ !definition_str)
            else if Str.string_match (Str.regexp 
            " *\\(<GBSeqid>gi|+\\)\\([a-zA-Z0-9 ,]+\\)") line 0 then
                output_string chout (">gi|" ^ (Str.matched_group 2 line) ^ 
                "|" ^ !accession_str)
            else if Str.string_match (Str.regexp 
            " *\\(<GBSeq_definition>+\\)\\([a-zA-Z0-9 ,]+\\)") line 0 then
                definition_str := (Str.matched_group 2 line) ^ "\n"
            else if Str.string_match (Str.regexp 
            " *\\(<GBSeq_sequence>+\\)\\([a-zA-Z]+\\)") line 0 then
                output_string chout ( (Str.matched_group 2 line) ^ "\n\n");
        done;
        open_in outName;
    with
    | End_of_file -> (open_in outName) 
