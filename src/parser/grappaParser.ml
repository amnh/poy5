    
    (* Given a text in_channel ch, make it a long string where the newlines
    * are replaced with single spaces. *)
    let stringify_channel ch = 
        let res = ref "" in
        try 
            while true do
                res := !res ^ " " ^ (input_line ch);
            done;
            ""
        with
        | End_of_file -> !res
   
    let rec make_list_of_genes input_list output_list =
        match input_list with
        | [] -> List.rev output_list
        | head :: tail -> 
                let blank_space = Str.regexp " +" in
                let name = List.hd (Str.split blank_space head ) in
                let genes = Str.string_after head (String.length name) in
                let genelist = Str.split blank_space genes in
                let thegenes = Array.of_list 
                   (List.map int_of_string genelist) in
                make_list_of_genes tail (thegenes::output_list)
 
    let of_channel ch =
        let line = stringify_channel ch
        and reg = Str.regexp " *>" in 
        if Str.string_match reg line 0 then
            begin
                let name_gene_list = Str.split reg line in
                let res = make_list_of_genes name_gene_list [] in 
                let res_array = Array.of_list res in
                Array.map (fun x -> FileContents.Genes x) res_array
            end
        else
            raise (E.Unsupported_file_format "Not a proper Grappa file") 
            
