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

let () = SadmanOutput.register "Lang" "$Revision: 1644 $"

(* All the functions to handle properly a POY language file *)

(* The parsed tree type *)
type tree = 
    | Empty 
    | Function 
    | Import 
    | Character 
    | Prob 
    | Alph

module Character = struct
    exception Empty_function
    exception Undefined_function of string
    exception Inconsistent_function_parameter
    exception Inconsistent_function_type
    exception Illegal_argument
    exception Illegal_parameter

    (* The only primitive types supported by the U language *)
    type utype = | Base | Seq | Int | Bool | Unkn


    (* A character definition specification statistics, this is for type
    * checking and counting occurrences of a function *)
    type chardef = {
        mutable functions : string list;
        mutable function_param : (string, utype list) Hashtbl.t;
        mutable function_types : (string, utype) Hashtbl.t;
        mutable function_count : (string, int) Hashtbl.t;
    }

    (* The initial internal state of the counter *)
    let _internal_state = ref { 
        functions = [];
        function_param = Hashtbl.create 10;
        function_types = Hashtbl.create 10;
        function_count = Hashtbl.create 10;
    }

    (* Returns the complete internal state *)
    let get () = !_internal_state

    (* Adds a function definition to the current state of the parser *)
    let rec addf = function
        | hd :: tl -> 
                !_internal_state.functions <- hd :: !_internal_state.functions;
                let tmp = List.map (fun x -> Unkn) tl in
                Hashtbl.replace !_internal_state.function_param hd tmp;
                Hashtbl.replace !_internal_state.function_types hd Unkn;
                Hashtbl.replace !_internal_state.function_count hd 0;
                List.iter (fun x -> addf [x]) tl;
        | _ -> raise Empty_function

    (* Checks if a function has already been defined *)
    let exists f = List.exists (fun x -> x = f) !_internal_state.functions 

    let getf f =
        try
            let a = Hashtbl.find !_internal_state.function_types f
            and b = Hashtbl.find !_internal_state.function_count f
            and c = Hashtbl.find !_internal_state.function_param f in
            c, a, b
        with
        | Not_found -> raise (Undefined_function f)

    let setft f t =
        try 
            match Hashtbl.find !_internal_state.function_types f with
            | Unkn -> Hashtbl.replace !_internal_state.function_types f t
            | pre when pre = t -> ()
            | _ -> raise Inconsistent_function_type
        with
        | Not_found -> raise (Undefined_function f)

    let setpt f it t =
        try
            let lst = Hashtbl.find !_internal_state.function_param f in
            let rec set_param counter = function
                | Unkn :: tl when counter = it -> t :: tl
                | (hd :: tl) as lst when ((counter = it) && (hd = t)) -> lst
                | hd :: tl when (counter < it) -> 
                        hd :: set_param (counter + 1) tl
                | [] -> raise Illegal_parameter
                | _ -> raise Inconsistent_function_parameter
            in
            let res = set_param 0 lst in
            Hashtbl.replace !_internal_state.function_param f res;
        with
        | Not_found -> raise (Undefined_function f)

    let getpt f it =
        try
            let lst = Hashtbl.find !_internal_state.function_param f in
            List.nth lst it
        with
        | Not_found -> raise (Undefined_function f)
        | Failure "nth" -> raise Illegal_argument

    let incr_count f =
        try
            let res = Hashtbl.find !_internal_state.function_count f in
            Hashtbl.replace !_internal_state.function_count f (res + 1);
        with
        | Not_found -> raise (Undefined_function f)

    let get_count f =
        try Hashtbl.find !_internal_state.function_count f
        with
        | Not_found -> raise (Undefined_function f)

    let init () = 
        addf [":"; "__ulang_hd_pre"; "__ulang_tl_pre"];
        addf ["pre"; "__ulang_pred"];
        addf ["suc"; "__ulang_suc"];
        addf ["hd"; "__ulang_hd"];
        addf ["tl"; "__ulang_tl"];
        setft ":" Seq;
        setft "pre" Int;
        setft "suc" Int;
        setft "hd" Base;
        setft "tl" Seq;
        setpt ":" 0 Base;
        setpt ":" 1 Seq;
        setpt "pre" 0 Int;
        setpt "suc" 0 Int;
        setpt "hd" 0 Seq;
        setpt "tl" 0 Seq

    (* If a new character is to be defined, this is the function to remove all
    * the crap from inside it *)
    let clean () = 
        _internal_state :=  { 
            functions = [];
            function_param = Hashtbl.create 10;
            function_types = Hashtbl.create 10; 
            function_count = Hashtbl.create 10; };
        init ()

    let elements c = 
        c.functions, c.function_param, c.function_types, c.function_count

        
end
