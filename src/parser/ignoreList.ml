let of_channel ch = 
    let input_handler = FileStream.stream_reader ch in
    let res = ref [] in
    try
        while true do
            let input = input_handler#read_line in
            res := input :: !res;
        done;
        List.rev !res
    with
    | _ -> List.rev !res
