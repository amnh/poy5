let anywhere_match regex line =
    try
        let _ = Str.search_forward regex line 0 in
        true
    with
    | _ -> false

let explode_filenames files = 
    let explode_filename file = 
        let file = FileStream.filename file in
        (* Warning: For some strange reason, if we run file_exists in
        * windows, then we can not catch the interrupt signal in the 
        * application. That suck for end users (ctr-C doesn't work). *)
        if Sys.os_type <> "Win32" && Sys.file_exists file then [file]
        else
            let ch = 
                let file = 
                    if Sys.os_type = "Win32" then file
                    else if Sys.os_type = "Unix" then
                        Str.global_replace (Str.regexp "\\\\ ") "\\ " file
                    else file
                in
                let line = 
                    match Sys.os_type with
                    | "Win32" -> ("dir /B \"" ^ file ^ "\" 2> NUL")
                    | _ -> "ls -1 " ^ file ^ " 2> /dev/null"
                in
                Unix.open_process_in line 
                in
            let res = IgnoreList.of_channel ch in
            close_in ch;
            match res with
            | [] -> 
                    let msg = "@[No@ file@ matching@ @{<b>" ^ StatusCommon.escape file ^ 
                    "@}@ found.@]" in
                    Status.user_message Status.Error msg;
                    failwith "File not found"
            | _ -> res
        in
        List.flatten (List.map explode_filename files)

