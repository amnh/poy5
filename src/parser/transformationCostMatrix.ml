let of_channel = Cost_matrix.Two_D.of_channel

let of_file ?(use_comb = true) file all_elements =
    let ch = FileStream.Pervasives.open_in file in
    let res = of_channel ~use_comb all_elements ch in
    ch#close_in;
    res

let of_channel_nocomb = Cost_matrix.Two_D.of_channel_nocomb

let fm_of_file file =
    (* explode a string around a character;filtering empty results *)
    let explode str ch =
        let rec expl s i l =
            if String.contains_from s i ch then
                let spac = String.index_from s i ch in
                let word = String.sub s i (spac-i) in
                expl s (spac+1) (word::l)
            else
                let final = String.sub s i ((String.length s)-i) in
                final::l in

        List.filter(fun x-> if x = "" then false else true) 
                    (List.rev (expl str 0 [])) in

    (* read a channel line by line and applying f into a list *)
    let rec read_loop f chan =
        try 
            let line = FileStream.Pervasives.input_line chan in
            (List.map (float_of_string) (f line ' ') ) :: read_loop f chan
        with e -> [] in

    let f = FileStream.Pervasives.open_in file in
    let mat = read_loop (explode) f in
    let _ = FileStream.Pervasives.close_in f in 
    mat

let of_list = Cost_matrix.Two_D.of_list

