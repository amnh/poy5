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

let () = SadmanOutput.register "Character2" "$Revision: 1644 $"

open Bigarray

exception Illegal_argument of string


let assign_codes start morphological molecular = 
    let tmp_start = ref start in
    let assign_to_morphological arr = 
        let len = Array.length arr in 
        let new_arr = Array.make len 0 in
        for i = 0 to len - 1 do
            new_arr.(i) <- !tmp_start;
            incr tmp_start;
        done;
        new_arr
    in
    let assign_to_molecular it = 
        let res = !tmp_start in
        incr tmp_start;
        res
    in
    let new_morpho = List.map assign_to_morphological morphological in
    let new_molecu = List.map assign_to_molecular molecular in
    new_morpho, new_molecu, !tmp_start

module type GENERAL_SET = sig
    (* All character sets use the following interface. Better implementations
    * can be used for each type of character. *)
    type c
    (* The type of a character *)
    type s 
    (* The set of characters type *)
    type item_id
    (* An item id in the set *)
    type item_code = int
    (* The unique code of a character in the set. Each character holds a unique
    * code (the item_code), which identifies it for specific operations. *)
    type data_element = Changed | Code | Total | Weight | Cost | Other of int
    (* It's the actual type of a single character. This provides a generalized,
    * abstract way to define the characters within a set *)
    val split : s -> s * s
    (* Divides a set of size n in two subsets of size \lfloor n/2 \rfloor and
    * \lceil n/2 \rceil (written a la LaTeX). *)
    val join : s -> s -> s
    (* Take two sets and join them in a single set. The join operation checks if
    * both sets have a disjoint set of item_codes; if not, an Illegal_argument
    * exception is raised. *)
    val set_all : data_element -> c -> s -> unit
    (* Sets a data element in all the characters in the set to some value *)
    val set_val : item_id -> data_element -> c -> s -> unit
    (* Sets a specific character, and a specific data element to
    * some specific value. *)
    val get_val : item_id -> data_element -> s -> c
    (* Returns the value of a given data_element of a specific character *)
    val copy_from : data_element -> data_element -> s -> unit
    (* Copies the contents of a data_element a to another data_element b in all
    * the characters of the set. *)
    val create : int -> int -> s
    (* Creates a new character set with x data elements and capacity for y
    * characters *)
    val data_elements : int
    (* Number of data_elements added by default in the set (standard
    * data_elements as defined in the following constants *)
    val cost : int
    (* Total cost of the character set in a tree *)
    val weight : int
    val total : int
    val code : int
    val changed : int
    val duplicate : s -> s
    (* Creates a new fresh copy of the set *)
    val clone : s -> s
    (* Synonym of duplicate *)
end

module type CHARACTER_SET = sig
    type s
    val downpass : s -> s -> s -> int
    val uppass : s -> s -> s -> s -> unit
end

(* Though not constrained, should follow GENERAL_SET. The Integer_Matrix module
* is used by the Additive and NonAdditive character sets. For any special
* data elements necessary in this matrix by any of those modules, increase the
* size in the required number of rows and start using from rows data_elements to
* the total size - 1. This module isn't visible from outside this library, so
* feel free to modify it, knowing that the changes are affecting the Additive
* and NonAdditive following modules. *)
module Integer_Matrix = struct 
    type s = (int,  Bigarray.int_elt,  Bigarray.c_layout) Array2.t  
    (* The set type *)
    type c = int
    (* The character element in the set type *)
    type item_id = int * int
    (* A coordinate in the matrix where the element is located *)
    type item_code = int
    (* The internal code of the character in the set *)
    type data_element = Changed | Code | Total | Weight | Cost | Other of int

    (* An element in the set *)

    (* ******* SOME CONSTANTS IN THE MODULE ******* *)
    (* Though the cost and the weight are elements that belong to the tree
    * directly, I will keep them mixed in the character definition to reduce the
    * total number of passing arguments later on. It will also simplify some
    * access operation later *)
    let data_elements = 5
    (* Total number of elements obligatory in an Integer_Matrix *)
    let changed = data_elements - 1
    (* Flag marking if an element has changed or not *)
    let code = changed - 1
    (* Code of the character in the character index *)
    let total = code - 1
    (* Total cost of the character evolution in the subtree rooted by the node
    * that contains the set. If there is no node or no tree, the total cost is 
    * 0, otherwise the total cost is the sum of the costs of the character in
    * the nodes it contains plus the characters cost in it's direct children. *)
    let weight = total - 1
    (* The weight of the characters in the set. Ideally, all the characters in a
    * set have the same weight, but this isn't necessarily true so far. *)
    let cost = weight - 1
    (* The evolutionary cost of the current character going from the node that
    * contains the set to it's two children. *)
    (* ******* END OF CONSTANTS IN THE MODULE ******* *)

    (* General error raiser for the module. Most of the possible errors that
    * could happen are illegal rows. *)
    let illegal_row x y =
        let msg = "An illegal row was attempted to access, using \
        coordinates " ^ string_of_int x ^ ", " ^ string_of_int y in
        raise (Illegal_argument msg)

    (* Takes a data element and outputs the row number where it occurs. *)
    let data_element_2_row_int = function
        | Changed -> changed
        | Code -> code
        | Total -> total
        | Weight -> weight
        | Cost -> cost 
        | Other a -> a

    (* Copies the contents of the row x to the row y in the matrix mtrx *)
    let copy_from mtrx x y =
        let x = data_element_2_row_int x
        and y = data_element_2_row_int y 
        and len = Array2.dim2 mtrx in
        try
            for i = 0 to len - 1 do
                mtrx.{y, i} <- mtrx.{x, i};
            done;
            ()
        with
        | _ -> illegal_row x y 

    (* Sets the value of position x, y in the matrix mtrx to v *)
    let set_val y x v mtrx = 
        let x = data_element_2_row_int x in
        try
            mtrx.{x, y} <- v
        with
        | _ -> illegal_row x y

    (* Retrieves the value in position x y of the matrix mtrx *)
    let get_val y x mtrx = 
        let x = data_element_2_row_int x in
        try
            mtrx.{x, y}
        with
        | _ -> illegal_row x y

    (* Sets the value of all the elements in the row x of matrix mtrx to v *)
    let set_all x v mtrx =
        let x = data_element_2_row_int x 
        and len = Array2.dim2 mtrx in
        try
            for i = len - 1 downto 0 do
                mtrx.{x, i} <- v;
            done;
            ()
        with
        | _ -> illegal_row x 0

    (* Splits the matrix mtrx in halves of size a b where a = b +- 1 if the
    * total size is odd, otherwise a = b 
    * creating an extra column - because for some reason don't use column 1*)
    let split mtrx =
        let cols = Array2.dim2 mtrx 
        and rows = Array2.dim1 mtrx in
        let fh = (cols - 1) / 2 in
        let sh = (cols - 1) - fh in
        let nf = 
            Array2.create Bigarray.int Bigarray.c_layout rows (fh + 1)
        and ns = 
            Array2.create Bigarray.int Bigarray.c_layout rows (sh + 1) 
        in
        for i = 0 to rows - 1 do
            for j = 1 to fh do
                try
                    nf.{i, j} <- mtrx.{i, j};
                with
                | _ -> illegal_row i j
            done;
        done;
        for i = 0 to rows - 1 do
            for j = 1 to sh do
                try
                    ns.{i, j} <- mtrx.{i, j + fh};
                with
                | _ -> illegal_row i j
            done;
        done;
        nf, ns

    let shared_item_codes a b = 
        let alen = Array2.dim2 a 
        and blen = Array2.dim2 b 
        and res = ref false 
        and i = ref 0 
        and j = ref 0 in
        while (not !res) && (!i < alen) do
            if (!j = blen) then begin
                i := !i + 1;
                j := 0;
            end else begin
                let av = get_val !i Code a
                and bv = get_val !j Code b in
                res := av = bv;
                j := !j + 1;
            end
        done;
        !res

    let join mtrx1 mtrx2 =
        let len_1 = Array2.dim2 mtrx1
        and len_2 = Array2.dim2 mtrx2 
        and rows = Array2.dim1 mtrx2  in
        if shared_item_codes mtrx1 mtrx2 then begin
            raise (Illegal_argument "There are shared item_codes between \
            elements in two sets to be joined.")
        end;
        let nl = len_1 + len_2 - 1 in (* get rid of the standard first blank *)
        let nch = Array2.create Bigarray.int Bigarray.c_layout rows nl in
        for i = 0 to rows - 1 do
            for j = 1 to len_1 - 1 do
                try
                    nch.{i, j} <- mtrx1.{i, j};
                with
                | _ -> illegal_row i j
            done;
        done;
        for i = 0 to rows - 1 do
            for j = len_1 to nl - 1 do
                try
                    nch.{i, j} <- mtrx2.{i, j + 1 - len_1};
                with
                | _ -> illegal_row i j
            done;
        done;
        nch

    let create a b = Array2.create Bigarray.int Bigarray.c_layout a b

    let get_changed it ch_s =
        let res = 
            try
                ch_s.{changed, it + 1} 
            with
            | _ -> illegal_row changed (it + 1)
        in
        match res with
        | 0 -> false
        | 1 -> true
        | e -> raise (Illegal_argument ("There is a strange value in the \
        changed flag. The unexpected value is: " ^ string_of_int e ^" in the \
        character " ^ string_of_int it ^ "."))

    let get_index_of_code cd mtx =
        let width = Array2.dim2 mtx in
        let rec find it =
            if width != it then begin
                if mtx.{code, it} = cd then it
                else find (it + 1)
            end else begin 
                let msg = "The character with code " ^ string_of_int cd ^ 
                " is not a member of the set." in
                raise (Illegal_argument msg)
            end
        in
        find 1

end

module NonAdditive = struct
    (* A set of non-additive characters *)
    (*module NonAdditive : CHARACTER_SET = struct *) 

    type c = int
    (* A single character in the set. Try to avoid using any individual
    * operation though. *)
    type s = Integer_Matrix.s
    (* Some data elements *)

    let rows = Integer_Matrix.data_elements + 2
    (* The total number of rows needed in non additive characters *)

    let prel_row = Integer_Matrix.data_elements
    let preli_row_data_element = Integer_Matrix.Other prel_row
    (* The preliminary state row *)

    let final_row = Integer_Matrix.data_elements + 1
    let final_row_data_element = Integer_Matrix.Other final_row
    (* The final state row *)

    let set_preliminary it ch ch_s = 
        Integer_Matrix.set_val it preli_row_data_element ch ch_s
    (* Set the preliminary cost of the character in position it to the value ch
    * in the character set ch_s. *)

    let get_preliminary it ch_s = 
        Integer_Matrix.get_val it preli_row_data_element ch_s 
    (* Gets the preliminary value of the character in position it of the
    * character set ch_s *)

    let set_final it ch ch_s = 
        Integer_Matrix.set_val it final_row_data_element ch ch_s 
    (* Set the final value of the character in position it of the character set
    * ch_s to ch *)

    let get_final it ch_s = 
        Integer_Matrix.get_val it final_row_data_element ch_s 
    (* Gets the final value of the character in position it of the
    * character set ch_s *)

    let set_cost it cost ch_s = 
        Integer_Matrix.set_val it Integer_Matrix.Cost cost ch_s 
    (* Set the cost of the character in position it of the character set ch_s to
    * ch *)

    let set_all_cost cost ch_s =
        Integer_Matrix.set_all Integer_Matrix.Cost cost ch_s
    (* Set the cost of all the items in the character set ch_s to cost *)

    let get_cost it ch_s = 
        Integer_Matrix.get_val it Integer_Matrix.Cost ch_s
    (* Gets the cost of the character in position it of the character set ch_s
    * *)

    let set_weight it w ch_s = 
        Integer_Matrix.set_val it Integer_Matrix.Weight w ch_s
    (* Sets the weight of the character in position it of the character set ch_s
    * to ch *)

    let get_weight it ch_s = 
        Integer_Matrix.get_val it Integer_Matrix.Weight ch_s
    (* Gets the weight of the character in position it of the character set
    * ch_s *)

    let set_total it ch ch_s = 
        Integer_Matrix.set_val it Integer_Matrix.Total ch ch_s
    (* Sets the total of the character in position it of the character set ch_s
    * to ch *)

    let get_total it ch_s = 
        Integer_Matrix.get_val it Integer_Matrix.Total ch_s 
    (* Gets the total of the character in position it of the character set ch_s
    * *)

    let set_char_code it ch ch_s = 
        Integer_Matrix.set_val it Integer_Matrix.Code ch ch_s 
    (* Sets the character code of the character in position it of the character
    * set ch_s to ch *)

    let get_char_code it ch_s = 
        Integer_Matrix.get_val it Integer_Matrix.Code ch_s 
    (* Gets the character code of the character in position it of the character
    * set ch_s *)

    let set_changed it ch_s = 
        Integer_Matrix.set_val it Integer_Matrix.Changed 1 ch_s 
    (* Sets the changed flag of the character in position it of the character
    * set ch_s *)

    let set_no_changed it ch_s = 
        Integer_Matrix.set_val it Integer_Matrix.Changed 0 ch_s 
    (* Unsets the changed flag of the character in position it of the character
    * set ch_s *)

    let get_changed it ch_s = 
        match Integer_Matrix.get_val it Integer_Matrix.Changed ch_s with
        | 0 -> false
        | _ -> true
    (* Gets the changed flag of the character in position it of the character
    * set ch_s *)

    let get_index_of_code cd mtx = 
        Integer_Matrix.get_index_of_code cd mtx
    (* Finds the index in the character set mtx of the character with code cd.
    * If the code isn't found, raise an Illegal_argument exception. *)

    external c_downpass : s  -> s -> s -> int -> int = "noaddset_CAML_downpass"

    external c_uppass : s -> s -> s -> s -> int -> unit = "noaddset_CAML_uppass"

    let downpass c1 c2 p = c_downpass c1 c2 p (Array2.dim2 c1)

    let uppass c1 c2 p c = c_uppass c1 c2 p c (Array2.dim2 c1)

    let set_all_weights v ch_s = 
        Integer_Matrix.set_all Integer_Matrix.Weight v ch_s

    let remove_weights ch_s = set_all_weights 0 ch_s 

    external c_set_pre_to_fin : Integer_Matrix.s -> int -> unit = 
        "noaddset_CAML_set_pre_to_fin"

    let preliminary_is_final cs =
        c_set_pre_to_fin cs (Array2.dim2 cs) 

    let join a b = Integer_Matrix.join a b

    let split b = Integer_Matrix.split b

    let create specs codes characters =
        let len = Array.length specs in
        let matrix = Integer_Matrix.create rows (len + 1) in
        for i = len - 1 downto 0 do
            let the_observed_set = 
                Parser.Hennig.Encoding.get_set specs.(i) 
            and weight = Parser.Hennig.Encoding.get_weight specs.(i) in
            let poss_char = All_sets.Integers.cardinal the_observed_set in
            if poss_char < 32 then begin
                match characters.(i) with
                | Parser.Unordered_Character v -> 
                        set_final i v matrix;
                        set_preliminary i v matrix;
                        set_char_code i codes.(i) matrix;
                        set_weight i weight matrix;
                | _ -> 
                        raise (Illegal_argument "Unexpected error.")
            end else begin
                let msg = "Non additive characters with more than 32 \
                states are not supported." in
                raise (Illegal_argument msg)
            end;
        done;
        matrix

end

(*module Additive : CHARACTER_SET = struct*)
module Additive = struct
    type c = int * int
    (* A single character in the set. Try to avoid using any individual
    * operation though. *)
    type s = Integer_Matrix.s

    (* Some data elements *)
    let rows = Integer_Matrix.data_elements + 4
    let prel_row_low = Integer_Matrix.data_elements
    let prel_row_up = Integer_Matrix.data_elements + 1
    let final_row_low = Integer_Matrix.data_elements + 2
    let final_row_up = Integer_Matrix.data_elements + 3

    let prel_row_low_data_element = Integer_Matrix.Other prel_row_low
    let prel_row_up_data_element = Integer_Matrix.Other prel_row_up
    let final_row_low_data_element = Integer_Matrix.Other final_row_low
    let final_row_up_data_element = Integer_Matrix.Other final_row_up

    let set_preliminary it (a, b) ch_s = 
        Integer_Matrix.set_val it prel_row_low_data_element a ch_s;
        Integer_Matrix.set_val it prel_row_up_data_element b ch_s

    let get_preliminary it ch_s = 
        let a = Integer_Matrix.get_val it prel_row_low_data_element ch_s
        and b = Integer_Matrix.get_val it prel_row_up_data_element ch_s in
        a, b

    let set_final it (a, b) ch_s = 
        Integer_Matrix.set_val it final_row_low_data_element a ch_s;
        Integer_Matrix.set_val it final_row_up_data_element b ch_s

    let get_final it ch_s = 
        let a = Integer_Matrix.get_val it final_row_low_data_element ch_s
        and b = Integer_Matrix.get_val it final_row_up_data_element ch_s in
        a, b

    let set_cost it ch ch_s = 
        Integer_Matrix.set_val it Integer_Matrix.Cost ch ch_s

    let get_cost it ch_s = 
        Integer_Matrix.get_val it Integer_Matrix.Cost ch_s

    let set_weight it ch ch_s = 
        Integer_Matrix.set_val it Integer_Matrix.Weight ch ch_s

    let get_weight it ch_s = 
        Integer_Matrix.get_val it Integer_Matrix.Weight ch_s 

    let set_total it ch ch_s = 
        Integer_Matrix.set_val it Integer_Matrix.Total ch ch_s

    let get_total it ch_s = 
        Integer_Matrix.get_val it Integer_Matrix.Total ch_s 

    let set_code it ch ch_s = 
        Integer_Matrix.set_val it Integer_Matrix.Code ch ch_s

    let get_code it ch_s = 
        Integer_Matrix.get_val it Integer_Matrix.Code ch_s

    let set_changed it ch_s = 
        Integer_Matrix.set_val it Integer_Matrix.Changed 1 ch_s

    let set_no_changed it ch_s = 
        Integer_Matrix.set_val it Integer_Matrix.Changed 0 ch_s

    let get_changed it ch_s = 
        match Integer_Matrix.get_val it Integer_Matrix.Changed ch_s with
        | 0 -> false
        | _ -> true

    let get_index_of_code cs cd = Integer_Matrix.get_index_of_code cs cd

    external c_downpass : s  -> s -> s -> int -> unit = "add_CAML_downpass"

    external c_uppass : s -> s -> s -> s -> int -> unit = "add_CAML_uppass"

    let downpass c1 c2 p = c_downpass c1 c2 p (Array2.dim2 c1)

    let uppass c1 c2 p c = 
        c_uppass c1 c2 p c (Array2.dim2 c1)
    
    external c_pre_to_fin : s -> int -> unit = "add_CAML_set_pre_to_fin"

    let preliminary_is_final cs = c_pre_to_fin cs (Array2.dim2 cs)

    let set_all_weights v ch_s = 
        Integer_Matrix.set_all Integer_Matrix.Weight v ch_s

    let remove_weights ch_s = set_all_weights 0 ch_s


    let of_array ch_arr =
        let len = Array.length ch_arr in
        let new_arr =
            Array2.create int c_layout rows (len + 1) in
        let build_single it (min, max, weight, code) =
            set_preliminary it (min, max) new_arr;
            set_code it code new_arr;
            set_weight it weight new_arr
        in
        Array.iteri (build_single) ch_arr

    let join a b = Integer_Matrix.join a b

    let split a = Integer_Matrix.split a  

    let create specs codes characters =
        let len = Array.length specs in
        let matrix = Integer_Matrix.create rows (len + 1) in
        for i = len - 1 downto 0 do
            let the_observed_set = 
                Parser.Hennig.Encoding.get_set (specs.(i)) 
            and weight = Parser.Hennig.Encoding.get_weight specs.(i) in
            let poss_char = All_sets.Integers.cardinal the_observed_set in
            if poss_char < 32 then begin
                match characters.(i) with
                | Parser.Ordered_Character v -> 
                        set_final i v matrix;
                        set_preliminary i v matrix;
                        set_code i codes.(i) matrix;
                        set_weight i weight matrix;
                | _ -> 
                        raise (Illegal_argument "Unexpected error.")
            end else begin
                let msg = "Non additive characters with more than 32 \
                states are not supported." in
                raise (Illegal_argument msg)
            end;
        done;
        matrix

end 

(** Operations on sets of sankoff characters. *)
module Sankoff_Precomputed = struct

    type c = int
    (* Some data elements *)

    type s = Integer_Matrix.s
    (* A single character in the set. Try to avoid using any individual
    * operation though. *)

    let rows = Integer_Matrix.data_elements + 2
    (* The total number of rows needed in non additive characters *)

    let prel_row = Integer_Matrix.data_elements + 1
    (* The row for preliminary states *)
    let final_row = Integer_Matrix.data_elements
    (* The row for final states *)

    external c_downpass : 
        s -> s -> s -> Cost_matrix.Three_D.m -> int -> int -> unit =
        "sankoff_CAML_downpass"

    let downpass ch1 ch2 par cm a_sz = c_downpass ch1 ch2 par cm a_sz

end

module Sankoff = struct

    (** A character state *)
    type c = int array

    (** The actual type of the set *)
    type s = c array

end

module Chromosome = struct

    type c 
    (* A single character in the set. Try to avoid using any individual
    * operation though. *)
    
    type s = c array
    (* Some data elements *)
    (* The constants cost, weight etc must match values in chromosome.h for
    * corresponding COST, WEIGHT etc. *)
    let cost = 0
    let weight = 1
    let total = 2
    let code = 3
    let changed = 4
    let none = 0
    let inversion = 1
    let breakpoint = 2
    let invandbreak = 3

    external c_read_input : in_channel -> int -> int -> s = "read_input2_CAML"

    (* reads input file and returns a tuple of number of genomes (or taxon) and
     * number of genes *) 
    external c_calc_num_genomes_genes : in_channel -> ( int * int ) = 
    "calc_num_genomes_genes2_CAML"


    (* takes a genome array and prints it *)
    external c_print_gene_matrix : s -> unit = 
        "print_gene_matrix2_CAML"

    (* takes a genome array and prints it *)
    external c_print_gene_array : c -> unit = 
    "print_gene_array2_CAML"
    
    external c_set_preliminary : int -> c ->  s -> unit = 
        "chromosome_CAML_set_preliminary"

    external c_get_preliminary : int ->  s -> c = 
        "chromosome_CAML_get_preliminary"
 
    external c_set_final : int -> c ->  s -> unit = 
        "chromosome_CAML_set_final"
 
    external c_get_final : int -> s -> c = 
        "chromosome_CAML_get_final"

    (* Set the preliminary code of the character in position it to the value ch
    * in the character set ch_s. *)
    let set_preliminary it ch ch_s = 
        c_set_preliminary it ch ch_s

    (* Gets the preliminary value of the character in position it of the
    * character set ch_s *)
    let get_preliminary it ch_s = 
        c_get_preliminary it ch_s

    (* Set the final value of the character in position it of the character set
    * ch_s to ch *)
    let set_final it ch ch_s = 
        c_set_final it ch ch_s

    (* Gets the final value of the character in position it of the
    * character set ch_s *)
    let get_final it ch_s = 
        c_get_final it ch_s
    
    external c_set_pre_to_fin : s -> unit = 
        "chromosome_CAML_set_pre_to_fin"

    let preliminary_is_final ch_s =
        c_set_pre_to_fin ch_s  
    
    external c_set_val : int -> int -> int -> s -> unit = 
        "chromosome_CAML_set_val"
    
    external c_get_val : int -> int -> s -> int = "chromosome_CAML_get_val"
    
    (* Set the cost of the character in position it of the character set ch_s to
    * nval *)
    let set_cost it nval ch_s =
        c_set_val cost it nval ch_s 

    external c_set_all_cost: int -> s -> unit =
        "chromosome_CAML_set_all_cost"

    (* Set the cost of all the items in the character set ch_s to cost *)
    let set_all_cost cst ch_s =
        c_set_all_cost cst ch_s

    (* Gets the cost of the character in position it of the character set ch_s*)
    let get_cost it ch_s = 
        c_get_val cost it ch_s

    
    (* Sets the weight of the character in position it of the character set ch_s
    * to nval *)
    let set_weight it nval ch_s = 
        c_set_val weight it nval ch_s 

    (* Gets the weight of the character in position it of the character set
    * ch_s *)
    let get_weight it ch_s = 
        c_get_val weight it ch_s

    (* Sets the total of the character in position it of the character set ch_s
    * to nval *)
    let set_total it nval ch_s = 
        c_set_val total it nval ch_s 

    (* Gets the total of the character in position it of the character set ch_s
    * *)
    let get_total it ch_s = 
        c_get_val total it ch_s

    (* Sets the character code of the character in position it of the character
    * set ch_s to nval *)
    let set_char_code it nval ch_s = 
        c_set_val code it nval ch_s

    (* Gets the character code of the character in position it of the character
    * set ch_s *)
    let get_char_code it ch_s = 
        c_get_val code it ch_s

    (* Sets the changed flag of the character in position it of the character
    * set ch_s *)
    let set_changed it ch_s = 
        c_set_val changed it 1 ch_s

    (* Unsets the changed flag of the character in position it of the character
    * set ch_s *)
    let set_no_changed it ch_s = 
        c_set_val changed it 0 ch_s

    (* Gets the changed flag of the character in position it of the character
    * set ch_s *)
    let get_changed it ch_s = 
        match c_get_val changed it ch_s with 
        | 0 -> false
        | _ -> true

    external c_get_index_of_code: int -> s -> int =
        "chromosome_CAML_get_index_of_code"
    
    (* Finds the index in the character set mtx of the character with code cd.
    * If the code isn't found, raise an Illegal_argument exception. *)
    let get_index_of_code cd mtx = 
        c_get_index_of_code cd mtx

    external c_set_all_weights: int -> s -> unit =
        "chromosome_CAML_set_all_weights"
    
    let set_all_weights v ch_s = 
        c_set_all_weights v ch_s

    let remove_weights ch_s = c_set_all_weights 0 ch_s 


    external c_downpass : s  -> s -> s -> int -> int = 
        "chromosome_CAML_downpass"

    external c_uppass : s -> s -> s -> s -> int -> unit = 
        "chromosome_CAML_uppass"

    let downpass ch1 ch2 p meth = c_downpass ch1 ch2 p meth
    
    let uppass ch1 ch2 gp p meth = c_uppass ch1 ch2 gp p meth

    external c_join : s -> s -> s = "chromosome_CAML_join"
    
    let join a b = c_join a b

    external c_split : s -> s * s = "chromosome_CAML_split"
    
    let split b = c_split b
    
    external c_create_empty_chrom : int -> int -> s = 
        "chromosome_CAML_create_empty"

    external c_set_pre_and_final : int -> int array -> s -> unit =
        "chromosome_CAML_set_pre_and_final"

    let get_numgenes characters =
        match characters.(0) with
        | Parser.Genes g -> Array.length g
        | _ -> raise (Illegal_argument "Unexpected error.")
        
    let create specs codes characters =
        let numgenes = get_numgenes characters   
        and numgenomes = Array.length codes in
        let matrix = c_create_empty_chrom numgenomes numgenes in
        for i = numgenomes - 1 downto 0 do
            match characters.(i) with
            | Parser.Genes v ->
                    c_set_pre_and_final i v matrix;
                    set_char_code i codes.(i) matrix;
            | _ -> 
                    raise (Illegal_argument "Unexpected error.")
        done;
        matrix
end
 
(************************************************************************)
(* Replacing missing character stuff                                    *)
let num_types = 3 

type types = | NonAdditive | Additive | Chromosome

type characters = 
    | NonAdditive_Char of NonAdditive.c | Additive_Char of Additive.c |
    Chromosome_Char of Chromosome.c

let same_type firstch secondch =
    match firstch, secondch with
    | NonAdditive_Char _ , NonAdditive_Char _ -> true
    | Additive_Char _ , Additive_Char _ -> true
    | Chromosome_Char _ , Chromosome_Char _ -> true
    | _, _ -> false
    
let get_type_index = function
            | NonAdditive  -> 1
            | Additive  -> 2
            | Chromosome  -> 4

let get_types c =
    match c with
    | NonAdditive_Char _ -> [NonAdditive]
    | Additive_Char _  -> [Additive]
    | Chromosome_Char _ -> [Chromosome]


let create_types_array somechar =
    ref (Array.make num_types somechar)

let initial () = 
    ref (Array.make num_types (NonAdditive_Char 1))
