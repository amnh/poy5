let do_text = ref false

let () = 
    let options = 
        ["-text", Arg.Set do_text, "Output the version in plain text"] 
    in
    Arg.parse options (fun _ -> ()) ""


let () = 
    if !do_text then
        print_string Version.version_string
    else begin
        Printf.printf "\\newcommand{\\smallbuildnumber}{%s}"
        Version.small_version_string;
        Printf.printf "\\newcommand{\\buildnumber}{%s}" Version.version_string;
    end
