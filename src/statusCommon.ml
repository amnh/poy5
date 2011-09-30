(* POY 4.0 Beta. A phylogenetic analysis program using Dynamic Homologies.    *)
(* Copyright (C) 2007  Andrés Varón, Le Sy Vinh, Illya Bomash, Ward Wheeler,  *)
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

let () = SadmanOutput.register "StatusCommon" "$Revision: 2871 $"

(* The common files for all the status interfaces. *)

external string_to_format : string -> ('a, 'b, 'c) format = "%identity"

type formatter_output = | Margin of int | Compress

let escape str =
    let str = Str.global_replace (Str.regexp "%") "%%" str in
    Str.global_replace (Str.regexp "@") "@@" str

module CommandCompletion = struct
    let commands = [
        "_attemptsdistr";
        "_breakvsjoin";
        "_cost";
        "_distances";
        "_mst";
        "_random";
        "_rootuniondistr";
        "_unionstats";
        "albert";
        "all";
        "all_roots";
        "alphabetic_terminals";
        "alternate";
        "aminoacids";
        "annchrom_to_breakinv";
        "annealing";
        "annotated";
        "approximate";
        "around";
        "as_is";
        "asciitrees" ;
        "at_random";
        "auto_sequence_partition";
        "auto_static_approx";
        "best";
        "better";
        "bfs";
        "bootstrap";
        "branch_and_bound";
        "breakinv";
        "breakinv_to_custom";
        "bremer";
        "build";
        "calculate_support";
        "cd";
        "characters";
        "chrom_breakpoint";
        "chrom_hom";
        "chrom_indel";
        "chrom_to_seq";
        "chromosome";
        "circular";
        "clades";
        "clear_memory";
        "clear_recovered";
        "codes";
        "collapse";
        "compare";
        "complex";
        "consensus";
        "constraint";
        "cross_references";
        "current_neighborhood";
        "custom_alphabet";
        "custom_to_breakinv";
        "data";
        "depth";
        "diagnosis";
        "direct_optimization";
        "distance";
        "drifting";
        "dynamic";
        "dynamic_pam";
        "echo";
        "end";
        "error";
        "estimate";
        "exact";
        "exhaustive_do";
        "exit";
        "false";
        "file";
        "files";
        "first";
        "fixed_states";
        "forest";
        "fuse";
        "gamma";
        "gap_opening";
        "gb";
        "genome";
        "graphconsensus";
        "graphsupports";
        "graphtrees";
        "help";
        "hennig";
        "history";
        "hits";
        "ia";
        "implied_alignments";
        "individual";
        "info";
        "init3D";
        "inspect";
        "iterations";
        "iterative";
        "jackknife";
        "keep";
        "last";
        "level";
        "likelihood";
        "load";
        "locus_breakpoint";
        "locus_indel";
        "locus_inversion";
        "log";
        "lookahead";
        "m";
        "margin";
        "max_time";
        "mb";
        "med_approx";
        "median";
        "median_solver";
        "memory";
        "model";
        "min_loci_len";
        "min_rearrangement_len";
        "min_time";
        "missing";
        "multi_static_approx";
        "names";
        "new";
        "newick";
        "nodes";
        "nolog";
        "nomargin";
        "none";
        "normal_do";
        "normal_do_plus";
        "not";
        "nucleotides";
        "of_file";
        "once";
        "optimal";
        "orientation";
        "origin_cost";
        "output";
        "perturb";
        "phastwinclad";
        "prealigned";
        "prioritize";
        "pwd";
        "quit" ;
        "random";
        "randomize_terminals";
        "randomized";
        "ratchet";
        "read";
        "recover";
        "rediagnose";
        "redraw";
        "remove";
        "rename";
        "repeat";
        "replace";
        "report";
        "resample";
        "root";
        "run";
        "s";
        "save";
        "script_analysis";
        "search";
        "search_based";
        "searchsstats";
        "sectorial";
        "seed";
        "seed_length";
        "select";
        "seq_stats";
        "seq_to_chrom";
        "sequence_partition";
        "set";
        "spr";
        "static";
        "static_approx";
        "store";
        "supports";
        "swap";
        "swap_med";
        "synonymize";
        "target_cost";
        "tbr";
        "tcm";
        "td";
        "terminals";
        "theta";
        "threshold";
        "ti";
        "timedprint";
        "timeout";
        "timer";
        "total";
        "trailing_deletion";
        "trailing_insertion";
        "trajectory";
        "transform";
        "trees";
        "treestats";
        "true";
        "uniform";
        "unique";
        "use";
        "version";
        "visited";
        "weight";
        "weightfactor";
        "weighting";
        "wipe";
        "within"
    ]
end

module Format = Format

module Files = struct
    
    type auto_complete = Command | Filename

    (* A module to handle filenames *)
    let hd_previous_results = ref []
    let tl_previous_results = ref []
    let prefix = ref ""

    let is_prefix x y =
        (String.length x <= String.length y) &&
        (x = (String.sub y 0 (String.length x)))

    let dirname str = 
        Filename.dirname str

    let last_basedir = ref ""

    let update_hd_tl_for_prefix_command str =
        last_basedir := "";
        let prefixes = List.filter (is_prefix str) (
        CommandCompletion.commands) in
        let lst = List.sort String.compare prefixes in
        tl_previous_results := lst 

    let update_hd_tl_for_prefix str =
        let dir = dirname str in
        last_basedir := dir;
        let filename = Filename.basename str in
        let prefixes = List.filter (is_prefix filename) (Array.to_list
        (Sys.readdir dir)) in
        let lst = List.sort String.compare prefixes in
        tl_previous_results := lst

    let rec find_next_match str = 
        match !tl_previous_results with
        | h :: t -> 
                hd_previous_results := h :: !hd_previous_results;
                tl_previous_results := t;
                h
        | [] -> 
                tl_previous_results := List.rev !hd_previous_results;
                hd_previous_results := [];
                ""

    (* An implementation of complete_filename in OCaml that matches the behavior
    * of readline's function. In this way we simplify possible interoperability
    * between interfaces. *)
    let complete_filename kind str state =
        let f = 
            match kind with
            | Command -> update_hd_tl_for_prefix_command
            | Filename -> update_hd_tl_for_prefix
        in
        if state = 0 then f str;
        find_next_match str

    let counter = ref 0
    let last_filename = ref ""

    (* If the [str] is the same as the previous request, the state remains 1
    * otherwise 0 *)
    let complete_filename kind str = 
        let comparator, sep_string = 
            if Sys.os_type = "Win32" then '\\', "\\" else '/', "/"
        in
        if str = !last_filename then 
            counter := 1
        else begin 
            counter := 0;
            last_basedir := "";
            hd_previous_results := [];
            tl_previous_results := [];
        end;
        match kind with
        | Filename ->
                last_filename := str;
                let complete_filename = complete_filename kind str !counter in
                if complete_filename = "" then complete_filename
                else if !last_basedir.[(String.length !last_basedir) - 1] <>
                comparator then
                    !prefix ^ !last_basedir ^ sep_string ^ complete_filename 
                else 
                    !prefix ^ !last_basedir ^ complete_filename 
        | Command ->
                last_filename := str;
                complete_filename kind str !counter


    let opened_files = Hashtbl.create 7

    let is_open file = Hashtbl.mem opened_files file

    let close_all_opened_files () =
        let closer _ (ch, _, close) =
            close ();
            match ch with
            | `NotCompressed ch -> close_out ch
            | `Zlib ch -> Gz.close_out ch
        in
        Hashtbl.iter closer opened_files 

    let assign_formatter_output f fo_ls = 
        List.iter (fun fo -> 
                       match fo with
                       | Margin m -> 
                               if m = Format.pp_get_margin f () then ()
                               else begin
                                   Format.pp_print_flush f ();
                                   Format.pp_set_margin f m
                               end
                       | Compress -> ()
                  ) fo_ls

    let get_margin filename =         
        match filename with
        | None -> Format.get_margin ()
        | Some filename ->
              try  
                  let _, fo, _ =  Hashtbl.find opened_files filename in   
                  Format.pp_get_margin  fo ()   
              with  Not_found -> Format.get_margin () 



    let openf ?(mode = `Append) name fo_ls = 
        if Hashtbl.mem opened_files name then 
            let _, f, _ = Hashtbl.find opened_files name in 
            let _ = assign_formatter_output f fo_ls in
            f
        else 
            let is_compressed = 
                List.exists (function Compress -> true | _ ->
                false) fo_ls 
            in
            (let file_options = 
                if is_compressed then
                    [Pervasives.Open_wronly; Pervasives.Open_trunc;
                    Pervasives.Open_creat; Pervasives.Open_binary]
                else
                match mode with
                | `Append ->
                        [Pervasives.Open_wronly; Pervasives.Open_append;
                        Pervasives.Open_creat; Pervasives.Open_text]
                | `New ->
                        [Pervasives.Open_wronly; Pervasives.Open_trunc;
                        Pervasives.Open_creat; Pervasives.Open_text]
            in
            let ch = 
                if not is_compressed then
                    `NotCompressed (open_out_gen file_options 0o644 name)
                else `Zlib (Gz.open_out name)
            in
            let f, before_closing = 
                match ch with
                | `NotCompressed ch ->
                    Format.formatter_of_out_channel ch, fun () -> ()
                | `Zlib ch ->
                    let flush () = Gz.flush ch in
                    let out string initial length =
                        for i = initial to initial + length - 1 do
                            Gz.output_char ch string.[i];
                        done
                    in
                    Format.make_formatter out flush,
                    fun () ->  Gz.flush ~flush:Gz.Finish_flush ch
            in
            assign_formatter_output f fo_ls;
            Hashtbl.add opened_files name (ch, f, before_closing);
            f)

    let set_margin filename margin = 
        match filename with
        | None -> ()
        | Some filename ->
                let _ = openf filename [Margin margin] in
                ()


    let flush () = 
        Hashtbl.iter (fun _ (ch, f, _) -> 
            Format.pp_print_flush f ();
            match ch with
            | `NotCompressed ch -> flush ch
            | `Zlib ch -> Gz.flush ~flush:Gz.Full_flush ch) opened_files 


    let closef name () =
        if Hashtbl.mem opened_files name then
            (let ch, f, close = Hashtbl.find opened_files name in
            Format.pp_print_flush f ();
            let () = 
                match ch with
                | `NotCompressed ch -> 
                        Pervasives.flush ch;
                        close ();
                        Pervasives.close_out ch;
                | `Zlib ch -> 
                        Gz.flush ~flush:Gz.Full_flush ch;
                        close ();
                        Gz.close_out ch;
            in
            Hashtbl.remove opened_files name)
        else ()

    let rec channel name =
        if Hashtbl.mem opened_files name then
            let ch, _, _ = Hashtbl.find opened_files name in
            ch
        else 
            let _ = openf name [] in
            channel name

    let _ = 
        at_exit close_all_opened_files

end

module Tables = struct

    let generate_string row_num item x = 
        let tag = 
            if 0 = item mod 2 then
                if 0 = row_num mod 2 then
                    "@{"
                else "@{<c:black_whit>"
            else 
                if 0 = row_num mod 2 then
                    "@{<c:red>"
                else "@{<c:red_whit>"
        in
        let x = if x = "" then " " else x in
        string_to_format (tag ^ x ^ "@}")

    let output_row row_num f row before after =
        Array.iteri (fun col x ->
            before ();
            let string = 
                if 0 = row_num mod 2 then string_to_format "@{" 
                else string_to_format "@{<c:black_whit>"
            in
            Format.fprintf f string;
            (Format.fprintf f (generate_string row_num col x) : unit);
            Format.fprintf f "@}";
            after ()) row

    (* Formatting and outputing tables *)
    let output f do_close closer v = 
        (* We need to set the tabs first *)
        let widths = Array.create (Array.length v.(0)) 0 in
        Array.iter (Array.iteri (fun p x -> 
            widths.(p) <- max (String.length x) widths.(p))) v;
        let first_row = Array.map (fun x -> String.make (x + 1) ' ') widths in
        Format.pp_open_tbox f ();
        output_row 0 f first_row (Format.pp_set_tab f) (fun () -> ());
        Format.pp_force_newline f ();
        Array.iteri (fun item x -> 
            output_row item f x (Format.pp_print_tab f) (fun () -> ());
        ) v;
        Format.pp_close_tbox f ();
        if do_close then closer ();

end

let information_output = ref None

let set_information_output (filename : string option) =
    information_output := filename

let redirect_information () = 
    match !information_output with 
    | None -> false
    | Some _ -> true

let information_redirected () = !information_output

let redirect_filename () =
    match !information_output with
    | None -> failwith "No redirection"
    | Some f -> f

let output_status_to_formatter formatter maximum achieved header suffix =
    let output_to_formatter formatter =
        match maximum, achieved with
        | None, 0 ->
                Format.fprintf formatter "@[%s\t@;@[%s@]@ @]@." header 
                suffix
        | None, n ->
                Format.fprintf formatter
                "@[%s\t%d\t@;@[%s@]@ @]@." header achieved
                suffix
        | Some max, _ ->
                Format.fprintf formatter
                "@[%s\t%d of %d@;@[%s@]@ @]@." header achieved
                max suffix
    in
    output_to_formatter formatter;
    match information_redirected () with
    | Some filename ->
            let f = Files.openf filename [] in
            output_to_formatter f
    | None -> ()

let process_parallel_messages printer =
IFDEF USEPARALLEL THEN
        let max = Mpi.comm_size Mpi.comm_world in
        let rec check_for_message cnt =
            if 0 = Mpi.comm_rank Mpi.comm_world && cnt < max then
                let gotit, rank, tag = 
                    Mpi.iprobe Mpi.any_source 3 Mpi.comm_world
                in
                if not gotit then ()
                else
                    let t, (msg : string) = 
                        Mpi.receive rank tag Mpi.comm_world in
                    let () = printer t msg in
                    check_for_message (cnt + 1)
            else ()
        in
        check_for_message 0
ELSE
        ()
END
