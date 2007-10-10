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

let () = SadmanOutput.register "AsciiTree"  "$Revision"

let ndebug = true

type xml_tree = Leaf of string | Node of string * xml_tree list

let print_points (x, y, code) =
    print_string "x = ";
    print_int x;
    print_string " y = ";
    print_int y;
    print_string (" code = " ^ code);
    print_newline ()

let sort_tree tree = 
    let rec aux_sort_tree tree =
        match tree with
        | Parser.Tree.Leaf _ -> tree, 0, 1
        | Parser.Tree.Node (chld, data) ->
                let chld = List.map aux_sort_tree chld in
                let chld = List.sort (fun (_, x, a) (_, y, b) -> 
                    match y - x with
                    | 0 -> b - a
                    | n -> n) chld in
                let max_depth, total_children, chld =
                    List.fold_left (fun (dep, ch, items) (x, y, z) ->
                        (max dep y), (ch + z), (x :: items)) (0, 0, []) chld
                in
                Parser.Tree.Node (chld, data), max_depth, total_children
    in
    let tree, _, _ = aux_sort_tree tree in
    tree

(** [draw b a] outputs a (crappy) ascii tree [a] in the channel [b]. [sep]
* establishes the number of lines separating each leaf and [bd] the number of
* columns between vertices in the tree. *)
let to_matrix ?(sep=4) ?(bd=4) include_interior t =
    (* Calculate the depth, number of leafs and the lenght of the longest string
    * of the leaves. *)
    let t = sort_tree t in
    let rec depth_and_leafs t = 
        match t with
        | Parser.Tree.Leaf str -> 1, 1, String.length str
        | Parser.Tree.Node (c, _) -> 
                let res = List.map depth_and_leafs c in
                let for_fold = fun (u, v, w) (x, y, z) ->
                    max u x, v + y, max w z
                in
                let depth, leafs, str = List.fold_left for_fold (0, 0, 0) res in
                depth + 1, leafs, str
    in
    (* Creates a matrix of characters with the necessary size to draw the tree
    * *)
    let create_matrix depth leafs strlen =
        if not ndebug then begin
            print_endline "Creating matrix:";
            print_endline ("Depth : " ^ string_of_int depth);
            print_endline ("Leafs : " ^ string_of_int leafs);
            print_endline ("Strlen : " ^ string_of_int strlen);
        end;
        let height = leafs + ((leafs - 1) * sep)
        and width = depth + ((depth - 1) * bd) + strlen in
        let matrix = Array.make_matrix height width ' ' in
        matrix, height, width
    in
    (* A function that fills in the contents of the tree [t] in the [matrix].
    * The total [width], [height] and [strlen] of the tree must be known (and
    * are calculated by the [depth_and_leafs] function). *)
    let set_values t matrix width height strlen =
        (* A horizontal line going from [min] to [max] in the [y] coordinate *)
        let fill_horizontal min max y =
            assert (min >= 0);
            assert (
                if min <= max then true
                else begin
                    Status.user_message Status.Error 
                    ("Max is " ^ string_of_int max ^ " and min is " ^
                    string_of_int min);
                    false
                end);
            assert (max < Array.length matrix.(0));
            assert (y < Array.length matrix);
            for i = min to max do
                try matrix.(y).(i) <- '-' 
                with e -> raise e
            done
        (* Same as previous but for a vertical line *)
        and fill_vertical min max x =
            assert (min >= 0);
            assert (min <= max);
            assert (max < Array.length matrix);
            assert (x < Array.length matrix.(0));
            for i = min to max do
                try matrix.(i).(x) <- '|'
                with e -> raise e
            done
        (* Write a string [str] in the [y] coordinate. Strings are written in
        * the right margin of the matrix. *)
        and fill_string str x y =
            assert (y < Array.length matrix);
            let len = String.length str in
            let width = Array.length (matrix.(0)) in
            for i = 0 to (len - 1) do
                try matrix.(y).(x + i + 1) <- str.[i];
                with e -> 
                    print_endline ("y :" ^ string_of_int y);
                    print_endline ("width : " ^ string_of_int width);
                    print_endline ("strlen : " ^ string_of_int strlen);
                    print_endline ("i : " ^ string_of_int i);
                    raise e
            done
        in
        (* Recursively process the tree [t], starting in the [x], [y] coordinate
        * as it's top left corner in [matrix]. *)
        let rec filler t x y =
            match t with
            | Parser.Tree.Leaf str ->
                    let width = Array.length (matrix.(0)) in
                    fill_string str (width - (strlen + 3)) y;
                    fill_horizontal x (width - (strlen + 3)) y;
                    y, y, y
            | Parser.Tree.Node (chld, str) ->
                    let process_children = 
                        fun (m, y, _) t -> 
                            let n, max, median = filler t (x + bd) y in
                            (min m n), max + sep, median
                    in
                    let min, max, last_median = 
                        List.fold_left process_children (max_int, y, y) chld 
                    in
                    let max = max - sep in
                    let median = (min + max) / 2 in
                    if str <> "forest" then begin
                        if include_interior then fill_string str (x + 1) (median + 1);
                        fill_horizontal x (x + bd) median;
                        fill_vertical min last_median (x + bd);
                    end else ();
                    median, max, median
        in
        let _ = filler t 0 0 in
        ()
    in
    let depth, leafs, strlen = depth_and_leafs t in
    let matrix, height, width = create_matrix depth leafs strlen in
    set_values t matrix width height strlen;
    (* Now print it out filtering out empty lines *)
    matrix

let draw ?(sep = 4) ?(bd=4) include_interior ch t =
    let matrix = to_matrix ~sep:sep ~bd:bd include_interior t in
    Array.iter (fun x -> 
        if Array.fold_left (fun ac y -> ac || (y != ' ')) false x then begin
            Array.iter (output_char ch) x; 
            output_string ch "\n";
        end else ()) matrix;
    print_newline ()

let to_string ?(sep = 4) ?(bd = 4) include_interior t = 
    let t = sort_tree t in
    let matrix = to_matrix ~sep:sep ~bd:bd include_interior t in
    let buffer = Buffer.create (Array.length matrix.(0)) in
    Buffer.add_string buffer "@[<v>";
    let char_adder pos c =
        if c = '|' then begin
            begin match pos mod 5 with
            | 0 -> Buffer.add_string buffer "@{<c:red>"
            | 1 -> Buffer.add_string buffer "@{<c:green>"
            | 2 -> Buffer.add_string buffer "@{"
            | 3 -> Buffer.add_string buffer "@{<c:cyan>"
            | 4 -> Buffer.add_string buffer "@{<c:yellow>"
            | _ -> failwith "impossible"
            end;
            Buffer.add_char buffer c;
            Buffer.add_string buffer "@}";
        end else Buffer.add_char buffer c
    in
    Array.iter (fun x ->
        if Array.fold_left (fun ac y -> ac || (y != ' ')) false x then begin
            Buffer.add_string buffer "@[<h>";
            Array.iteri char_adder x;
            Buffer.add_string buffer "@]%!";
            Buffer.add_string buffer "@,";
        end else ()) matrix;
    Buffer.add_string buffer "@]%!";
    Buffer.contents buffer

(** Outputs the tree [t] in channel [ch] using parenthetical notation. *)
let draw_parenthesis ch t = 
    let t = sort_tree t in
    let rec printer ch t = 
        match t with
        | Parser.Tree.Leaf str ->
                output_string ch str;
                output_string ch " ";
        | Parser.Tree.Node (chld, str) ->
                output_string ch "( ";
                List.iter (printer ch) chld;
                output_string ch ") ";
                output_string ch "[";
                output_string ch str;
                output_string ch "] "

    in
    printer ch t;
    output_string ch "\n" 

let for_formatter ?(separator = " ") newick leafsonly t =
    let t = sort_tree t in
    let separator = 
        match newick with 
        | true -> ","
        | false -> separator
    in 

    let rec generator () =
        let sep = ref "" in
        fun acc t ->
            match t with
            | Parser.Tree.Leaf str ->
                  let res = acc ^ !sep ^ str in
                  sep := "@," ^ separator;
                  res
            | Parser.Tree.Node (chld, str) ->
                  let acc = acc ^ !sep ^ "(" in

                  let acc = List.fold_left (generator ()) acc chld in
                  let acc = acc ^ ")" in
                  sep := "@," ^ separator;
                  if str = "" || leafsonly then acc ^ "@," 
                  else acc ^ "[" ^ str ^ "]@," 
    in
    generator () "" t
