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

(* Helper functions *)
let truere = Str.regexp ".*true.*"
let get_option lst str =
    let _, str = 
        List.find (fun (a, _) -> Str.string_match (Str.regexp a) str 0) lst
    in
    str

let get_graphics str =
    get_option [(".*tk", "tk"); (".*ocaml", "ocaml"); (".*", "none")] str

let get_interface str =
    get_option [(".*gtk", "gtk2"); (".*ncurses", "ncurses"); (".*", "readline")] str

let is_true str = if Str.string_match truere str 0 then  "on" else "off"
let rephrase str = Str.global_replace (Str.regexp " +") "@ " str

(* The Version Values *)
let name = "Pelletier"
let string = "@[@[Welcome to @{<b>@{<c:blue>P@}O@{<c:red>Y@}@} 4.0 alpha,@ " ^
            "build@ " ^ BuildNumber.build ^ "@ \"@{<u>" ^ name ^ "@}\"!@]@." ^
                     rephrase ("@[compiled" ^ CompileFlags.time
                      ^ "with parallel " ^ is_true CompileFlags.str_parallel
                      ^ ", interface " ^ get_interface CompileFlags.str_interface
                      ^ ", graphics " ^ get_graphics CompileFlags.str_graphics ^
                      "@]@,@[" ^
                     "POY version 4.0 Beta, Copyright (C) 2007  Andres Varon, Le Sy Vinh, Illya Bomash, Ward Wheeler, and the American Museum of Natural History.  POY 4.0 comes with ABSOLUTELY NO WARRANTY; This is free software, and you are welcome to redistribute it under the GNU General Public License Version 2, June 1991.@]@]")
