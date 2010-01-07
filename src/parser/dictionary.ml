let gen_of_channel adder ch = 
    let input_handler = FileStream.stream_reader ch in
    let rec reader counter = 
        try
            let line = input_handler#read_line in
            match Str.split (Str.regexp " +") line with
            | [] | [_] -> 
                    let msg = 
                        ("Line " ^ string_of_int counter ^ ": " ^ line)
                    in
                    raise (E.Illegal_dictionary msg)
            | hd :: tl ->
                    List.iter (fun x -> adder x hd) tl;
                    reader (counter + 1)
        with
        | _ -> ()
    in
    reader 1

let of_channel ch =
    let table = Hashtbl.create 1667 in
    gen_of_channel (Hashtbl.add table) ch;
    table

let of_channel_assoc ch =
    let table = ref [] in
    gen_of_channel (fun a b -> table := (a, b) :: !table) ch;
    List.rev !table

