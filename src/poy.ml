(* POY 5.0 Alpha. A phylogenetic analysis program using Dynamic Homologies.   *)
(* Copyright (C) 2011 Andrés Varón, Lin Hong, Nicholas Lucaroni, Ward Wheeler,*)
(* and the American Museum of Natural History.                                *)
(*                                                                            *)
(* This program is free software; you can redistribute it and/or modify       *)
(* it under the terms of the GNU General Public License as published by       *)
(* the Free Software Foundation; either version 2 of the License, or          *)
(* (at your option) any later version.                                        *)
(*                                                                            *)
(* This program is distributed in the hope that it will be useful,            *)
(* but WITHOUT ANY WARRANTY; without even the implied warranty of             *)
(* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the              *)
(* GNU General Public License for more details.                               *)
(*                                                                            *)
(* You should have received a copy of the GNU General Public License          *)
(* along with this program; if not, write to the Free Software                *)
(* Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301   *)
(* USA                                                                        *)

let () = SadmanOutput.register "Main" "$Revision: 2814 $"

let seed = truncate (Unix.time ())

let hostname = Unix.gethostname ()

let debug_pass_errors = false

let master =
    IFDEF USEPARALLEL THEN
        let my_rank = Mpi.comm_rank Mpi.comm_world in
        if my_rank <> 0 then SadmanOutput.do_sadman := false; 
        my_rank = 0 
    ELSE
        true
    END

let is_running_alone = ref false

let () =
    let () = Status.init () in
    try
        Arg.parse_argv Phylo.args Arguments.parse_list 
                        (Arguments.anon_fun `Filename) Arguments.usage
    with
        | Arg.Help _ ->
            Arg.usage Arguments.parse_list Arguments.usage;
            exit 0
        | Arg.Bad x ->
            prerr_string ("Bad argument: " ^ x);
            Arg.usage Arguments.parse_list Arguments.usage;
            exit 1

(** Drop the Start-Up Message and instructions to console **)
let () =
    let out = Status.user_message Status.Information in
    let rephrase str = Str.global_replace (Str.regexp " +") "@ " str in
    out Version.string;
    let () =
        IFDEF USEPARALLEL THEN
            let tsize = Mpi.comm_size Mpi.comm_world in
            if master && tsize > 1 then
                out (Printf.sprintf "Running@ in@ parallel@ with@ %d@ processes"
                                    (Mpi.comm_size Mpi.comm_world))
            else if tsize = 1 then
                out "Running@ sequentially."
            else ();
            let arr = Array.init tsize (fun x -> seed + x) in
            let seed = Mpi.scatter_int arr 0 Mpi.comm_world in
            ignore (Phylo.process_random_seed_set (Phylo.empty ()) seed)
        ELSE
            ignore (Phylo.process_random_seed_set (Phylo.empty ()) seed)
        END
    in
    out "";
    out "";
    out "";
    IFDEF USENCURSES THEN
        out (rephrase "@[Type commands in the middle window, titled Interactive \
        Console.@\nJob status will appear below, and output will appear here.@\nA \
        summary of POY's current state will appear to the right of the \
        console.@\nFor help, type @{<u>help()@}.@\n@\nEnjoy!@]")
    ELSE
        out (rephrase "@[For help, type @{<u>help()@}.@\n@\nEnjoy!@]")
    END


let regen_and_save filename codestring regenfn =
    let data = regenfn () in
    let outfile = open_out_bin filename in
    Marshal.to_channel outfile (codestring, data) [];
    close_out outfile;
    data

let unmarshal_or_regen filename codestring regenfn =
    try
        let infile = open_in_bin filename in
        let (code, a) = Marshal.from_channel infile in
        close_in infile;
        if code = codestring
        then a
        else regen_and_save filename codestring regenfn
    with Sys_error _ ->
        regen_and_save filename codestring regenfn

(* Useful regex for matching tokens in syntax errors *)
let token_regexp = Str.regexp "\\[\\([a-zA-Z_]*\\)\\]"

let rec all_matches ?(group=0) ?(start=0) regexp string =
    let matcher =
        if group = 0
            then Str.matched_string
            else Str.matched_group group 
    in
    try
        let m = Str.search_forward regexp string start in
        let m_string = matcher string in
        m_string :: (all_matches ~start:(m + 1) ~group regexp string)
    with _ -> []

let unknown_command_regexp = Str.regexp "^[ \t]*\\([a-zA-Z_]*\\)"

(* First load the input *)
let script =
    try let lst = (List.rev !(Arguments.input)) in
        Some (PoyCommand.read_script_files false lst)
    with | err ->
        let msg = StatusCommon.escape (Printexc.to_string err) in
        Status.user_message Status.Error msg;
        if !(Arguments.just_exit) then
            exit 1
        else None

(*** Begin output. *)
(* let () = SadmanOutput.establish_connection () *)
let () = SadmanOutput.output
    "<?xml version=\"1.0\"?><!-- -*- mode: xml; mode: auto-revert; -*- -->\n"
let () = SadmanOutput.output "<sadman version=\"1.0\">\n"
let output_status = Status.user_message Status.Status
let timer = Timer.start ()

let timestring =
    let time = Unix.time () in
    let tm = Unix.gmtime time in
    Printf.sprintf "%d-%02d-%02d %02d:%02d:%02d"
        (tm.Unix.tm_year + 1900) (tm.Unix.tm_mon + 1) tm.Unix.tm_mday
        tm.Unix.tm_hour tm.Unix.tm_min tm.Unix.tm_sec

let () = 
    Sadman.start "POY"
        ((Sadman.prefix "module-version" (SadmanOutput.get_module_versions ()))
         @ (SadmanOutput.runtime_properties)
         @ [("hostname", hostname); ("timestamp", timestring)])

let safe_exit () =
    let time = Timer.get_user timer in
    let () = Sadman.finish [("exit-code", "0");
                            ("execution-time", string_of_float time)] in
    SadmanOutput.output "</sadman>\n";
    IFDEF USEPARALLEL THEN
        Mpi.finalize ()
    ELSE
        ()
    END

let () = at_exit safe_exit



let () = 
    let initial_script = ref script in
    let script = ref (Phylo.empty ()) in
    let input = ref "" in
    let proc_command str =
        let () = Sys.catch_break true in
        try
            let command = 
                if master then
                    let comm =
                        match !initial_script with
                        | Some comm -> 
                                initial_script := None;
                                comm
                        | None -> []
                    in
                    if 0 = String.length str then comm
                    else begin
                        input := str;
                        comm @ (PoyCommand.of_string false str)
                    end
                else []
            in
            let command =
                IFDEF USEPARALLEL THEN
                    let command = Analyzer.analyze command in
                    let command = Mpi.broadcast command 0 Mpi.comm_world in
                    let size = Mpi.comm_size Mpi.comm_world in
                    Analyzer.parallel_analysis (Mpi.comm_rank Mpi.comm_world) 
                                               size command
                ELSE
                    match command with
                    | [_] -> command
                    | x -> Analyzer.analyze command
                END
            in
            let res = 
                Phylo.run
                    ~output_file:(!(Arguments.dump_file)) ~start:!script command 
            in
            script := res;
            if !Arguments.only_run_argument_script then exit 1
            else ()
        with
        | Camlp4.PreCast.Loc.Exc_located (_, PoyCommand.Exit) -> exit 0
        | Camlp4.PreCast.Loc.Exc_located (a, Stream.Error err) ->
            let beg = Camlp4.PreCast.Loc.start_off a
            and en = Camlp4.PreCast.Loc.stop_off a in
            let is_unknown = "illegal begin of expr" = err in
            Printf.ksprintf
                (Status.user_message Status.Error)
                ("@[<v 4>Command error between characters @{<b>%d@}"^^
                 " and @{<b>%d@}:@,@[%s@]@]\n")
                beg en (if is_unknown then "Unknown command" else err);
            Status.error_location beg en;
            let elements =
                if is_unknown
                then all_matches ~group:1 unknown_command_regexp !input
                else all_matches ~group:1 token_regexp err
            in
            List.iter HelpIndex.help_if_exists elements;
            if !(Arguments.just_exit) || not (Status.is_interactive ()) 
                then exit 1
                else ()
        | Sys.Break -> 
            Status.clear_status_subwindows ();
            Status.user_message Status.Error "Interrupted";
            ()
        | Methods.TimedOut ->
            Status.clear_status_subwindows ();
            let msg = "Search timed out" in
            Status.user_message Status.Information msg
        | err when (((not !(Arguments.just_exit)) && 
                (Status.is_interactive ())) && (not debug_pass_errors)) ->
            Status.clear_status_subwindows ();
            let msg = StatusCommon.escape (Printexc.to_string err) in
            Status.user_message Status.Error msg
        | err when !Arguments.just_exit && (not debug_pass_errors) ->
            let msg = StatusCommon.escape (Printexc.to_string err) in
            Status.user_message Status.Error msg;
            raise err
    in
IFDEF USEPARALLEL THEN
    if master then
        Status.main_loop proc_command
    else 
        while true do
            proc_command ""
        done
ELSE
    Status.main_loop proc_command
END
