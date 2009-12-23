let to_fasta file =
    let ch, file = FileStream.channel_n_filename file and
    start_string = ref ">" and
    outName, chout = Filename.open_temp_file "fasta" ".tmp" in
    try
        while true do
            let line = input_line ch in
            if Str.string_match (Str.regexp 
            "\\(^[^ ]+\\)") line 0 then 
                output_string chout 
                (!start_string ^ (Str.matched_group 1 line) ^ "\n")
            else
                if Str.string_match 
                (Str.regexp " *\\([0-9]+\\) *\\([a-zA-Z ]+\\)+") line 0 then
                    let temp = Str.matched_group 2 line in
                    let temp2 = Str.global_replace (Str.regexp "[ ]+")
                    "" temp in
                    output_string chout ((String.uppercase temp2) ^ "\n");
                    start_string := "\n>"
        done;
        close_out chout;
        outName
    with
    | End_of_file -> 
            close_out chout;
            outName

let convert_to_fasta file =
    let file = to_fasta file in
    open_in file

