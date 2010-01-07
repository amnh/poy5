let convert_to_fasta ?filename file =
    let ch, file = FileStream.channel_n_filename file 
    and outName, chout = 
    match filename with
    | None -> Filename.open_temp_file "fasta" ".tmp"
    | Some f -> f, open_out f
in
    try
        while true do
            let line = input_line ch in
            (*
            if Str.string_match (Str.regexp 
            " *\\(VERSION+ *\\) \\([a-zA-Z0-9.]+\\) *GI:\\([0-9]+\\)") 
            line 0 then 
                output_string chout 
                (">gi|" ^ (Str.matched_group 3 line) ^ "|" ^
                (Str.matched_group 2 line) ^ "|" ^ !definition_str)
            *)
            if Str.string_match (Str.regexp 
            " *\\(ORGANISM+ *\\) \\([a-zA-Z0-9 ,]+\\)") line 0 then begin
                let str = Str.matched_group 2 line in
                let str = Str.global_replace (Str.regexp " ") "_" str in
                output_string chout ">";
                output_string chout str;
                output_string chout "\n";
            end else
                if Str.string_match 
                (Str.regexp " *\\([0-9]+\\) \\([a-zA-Z ]+\\)+") line 0 then
                    let temp = Str.matched_group 2 line in
                    let temp2 = Str.global_replace (Str.regexp "[ ]+")
                    "" temp in
                    output_string chout (String.uppercase temp2) ;
                    output_string chout "\n";
        done;
        open_in outName;
    with
    | End_of_file -> (open_in outName) 

