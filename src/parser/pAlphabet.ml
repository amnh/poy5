let of_file fn orientation init3D =
    let file = FileStream.Pervasives.open_in fn in
    let alph = FileStream.Pervasives.input_line file in
    let default_gap = Alphabet.gap_repr in
    let elts = ((Str.split (Str.regexp " +") alph) @ [default_gap]) in
    let alph = Alphabet.of_string ~orientation:orientation
    elts default_gap None in
    let size = Alphabet.size alph in
    let level =  2 in (* we set level = 2 by default *)
    let alph, do_comb = 
        if orientation then alph, false
        else 
                Alphabet.explote alph level size, true
    in
    let tcm = 
        try
            let all_elements = -1 (* we don't allow ambiguities here *) in
            if do_comb then
                Cost_matrix.Two_D.of_channel 
                ~orientation:orientation ~level:level all_elements file 
            else
                Cost_matrix.Two_D.of_channel_nocomb
                ~orientation all_elements file
        with
        | Failure "No Alphabet" -> 
                assert false 
                (*
                let size = Alphabet.size alph in
                Status.user_message Status.Information ("I'm following this
                path with size " ^ string_of_int size ^ "!");
                Cost_matrix.Two_D.of_transformations_and_gaps false size 1 1
                *)
    in
    let tcm3 = 
        match init3D with
        | true -> Cost_matrix.Three_D.of_two_dim tcm
        | false  ->  Cost_matrix.Three_D.default 
    in 
    file#close_in;
    alph, tcm, tcm3
