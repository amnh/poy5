let rec read_sequence ch out_seq =
     let line = input_line ch in
     if (String.contains line '}') then 
         begin
             out_seq := !out_seq ^ line;
             ()
         end
     else
         begin
             out_seq := !out_seq ^ line;
             read_sequence ch out_seq
         end

let rec read_title ch out_seq =
     let line = input_line ch in
     if Str.string_match (Str.regexp 
            " *create-date") line 0 then ()
     else
         begin
             out_seq := !out_seq ^ line;
             read_title ch out_seq
         end

let convert_to_fasta file =
    let ch, file = FileStream.channel_n_filename file and
    accession_str = ref "" and
    found_gi_id = ref 0 and
    outName, chout = Filename.open_temp_file "fasta" ".tmp" in
    try
        while true do
            let line = input_line ch in
            if Str.string_match (Str.regexp 
            " *\\(accession+ *\\) \\([a-zA-Z0-9\"]+\\)") line 0 then
                begin
                    let accession_quotes = Str.matched_group 2 line in
                    accession_str := Str.global_replace (Str.regexp "[\"]+")
                        "" accession_quotes
                end 
            else if (Str.string_match 
            (Str.regexp " *\\(gi+ *\\) \\([0-9]+\\)") line 0) 
            && !found_gi_id == 0 then
                begin
                    output_string chout 
                    (">gi|" ^ (Str.matched_group 2 line) ^ "|" 
                    ^ !accession_str ^ "|" );
                    found_gi_id := 1;
                end
            else if Str.string_match (Str.regexp 
            " *\\(title+\\) \\([a-zA-Z0-9\" ,.]+\\)") line 0 then
                begin
                    let out_seq = ref "" in
                    out_seq := (Str.matched_group 2 line);
                    read_title ch out_seq;
                    let out_list = Str.split (Str.regexp "[\"]") 
                    !out_seq in 
                    output_string chout ((List.nth out_list 0) ^ "\n");
                end 
            else if (Str.string_match (Str.regexp " *seq-data") ) line 0 
                then
                    begin
                        let out_seq = ref "" in
                        read_sequence ch out_seq;
                        let out_list = Str.split (Str.regexp "[']") 
                        !out_seq in 
                        output_string chout (List.nth out_list 1);
                        output_string chout "\n\n";
                    end
            else if (Str.string_match (Str.regexp " *seq {") ) line 0 then
                found_gi_id := 0
        done;
        open_in outName;
    with
    | End_of_file -> open_in outName 

