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

(* $Id: nonaddCS.ml 2871 2008-05-23 17:48:34Z andres $ *)
let () = SadmanOutput.register "NonaddCS.nonadd_v" "$Revision: 2871 $"


(** char_nonadd_c.ml implements sets of equally-weighted non-additive characters
    in C.  These sets are immutable but can share data through reference
    counting on the C side.  C is used mainly to ensure fast median-taking;
    medians are made with vectorizable bitwise operators.

    This implementation can handle sets of size up to the word size in bits,
    i.e., 32 characters on 32-bit machines, 64 characters on 64-bit machines.
    Using n-1 characters (31, 63) ensures complete interoperability with OCaml;
    using 32 or 64 allows you to use all the standard functions, but some of the
    raw operations included in char_nonadd_c.ml will fail.  However, these
    should only be used for debugging in the first place.

    These sets fit the Character.CharacterSet and the Character.CHARACTER
    signatures.  They are implemented as custom blocks in C, allowing for
    serialization, deserialization, and comparison.  The functions included here
    mainly come from the interfaces listed above;  any additional functions are
    explained in the .mli file, or are commented below. *)

(* Debugging mode for this module *)
let ndebug = false

(* This is a custom block *)
type ct
type t = {
    codes : int array;
    data : ct;
}
type cu = ct

(* The custom block for unions of characters *)
type u = t

(* The block contains integers.  However, note that using elements explicitly
   from OCaml restricts you to cardinalities of 31 or 63.  These pairs should
   consist of (code, data). *)
type e = int * int

external register_unmarshal : unit -> unit =
        "char_nonadd_CAML_register_unmarshal"

let _ = register_unmarshal ()

external cardinal : ct -> int = "char_nonadd_CAML_cardinal"
external union : cu -> ct -> cu -> cu -> unit = "char_nonadd_CAML_basic_union"

external union_states :  cu -> cu -> ct -> unit = "char_nonadd_CAML_basic_union_par"

let to_union x = x

external code : ct -> int = "char_nonadd_CAML_code"
let elt_code (code, _) = code
let block_size = NONADDSIZE

(* This is the default cost of a new element.  Unfortunately, we are not always
   given a cost when creating an element, so this is necessary. *)
let def_cost = 1.

external make_new : int -> int -> ct = "char_nonadd_CAML_make_new"
external make_new_unsafe : int -> int -> ct = "char_nonadd_CAML_make_new_unsafe"
external set_elt : ct -> int -> int -> unit = "char_nonadd_CAML_set_elt"
external set_elt_bit : ct -> int -> int -> unit = "char_nonadd_CAML_set_elt_bit"
external set_elt_code : ct -> int -> int -> unit
    = "char_nonadd_CAML_set_elt_code"
external basic_median : ct -> ct -> ct -> unit = "char_nonadd_CAML_basic_median"
(** [median_3 parent node child1 child2] *)
external median_3 : ct -> ct -> ct -> ct -> ct -> unit = "char_nonadd_CAML_median_3"
let median_3 a b c d =
    let new_item = make_new_unsafe (cardinal a.data) (code a.data) in
    median_3 new_item a.data b.data c.data d.data;
    { codes = a.codes; data = new_item}

external reroot_median : ct -> ct -> ct -> unit = "char_nonadd_CAML_reroot_median"
let reroot_median a b =
    let new_item = make_new_unsafe (cardinal a.data) (code a.data) in
    reroot_median new_item a.data b.data;
    { codes = a.codes; data = new_item}

(*** random characters: *)
(* how many in set, value range, code *)
type gen = int * int * int

(* random number between a and b *)
let rand_btwn a b = (Random.int (b - a)) + a

let rand_gen () = (rand_btwn 100 1000, rand_btwn 5 8, Random.int 1073741823)
let make_rand (num, outof, code) =
    let c = make_new num code in
    for i = 0 to num - 1 do
        set_elt c i (int_of_float (2. ** (float_of_int (Random.int outof))))
    done;
    c


(*** character and character set functions: *)
let median _ a b = 
    let new_item = make_new_unsafe (cardinal a.data) (code a.data) in
    basic_median new_item a.data b.data;
    { codes = a.codes; data = new_item }

external distance_ext : ct -> ct -> int = "char_nonadd_CAML_distance"
let distance a b = float_of_int (distance_ext a.data b.data)
external distance_union_ext : cu -> cu -> int = "char_nonadd_CAML_distance"
let distance_union a b = float_of_int (distance_union_ext a.data b.data)
external median_cost : ct -> float = "char_nonadd_CAML_median_cost"

let median_cost x = median_cost x.data

let union a b c =
    assert (cardinal a.data = cardinal b.data);
    assert (cardinal b.data = cardinal c.data);
    let new_item = make_new_unsafe (cardinal a.data) (code a.data) in
    union new_item a.data b.data c.data;
    { codes = a.codes; data = new_item} 

let union_states a b = 
    assert (cardinal a.data = cardinal b.data);
    let new_item = make_new_unsafe (cardinal a.data) (code a.data) in
    union_states new_item a.data b.data;
    { codes =a.codes; data =new_item }

external poly_items : cu -> int -> int = "char_nonadd_CAML_poly_items"
external cardinal_union : cu -> int = "char_nonadd_CAML_cardinal"

let cardinal_union x = cardinal_union x.data

let poly_saturation u pol =
    let c = cardinal_union u
    and p = poly_items u.data pol in
    (float_of_int p) /. (float_of_int c)

external to_int : ct -> int -> int = "char_nonadd_CAML_to_int"
external elt_to_list : ct -> int -> int list = "char_nonadd_CAML_elt_to_list"

external dist_2_ext : ct -> ct -> ct -> int = "char_nonadd_CAML_dist_2"
let dist_2 a b c = float_of_int (dist_2_ext a.data b.data c.data)
external distance_list : ct -> ct -> (int * float) list
    = "char_nonadd_CAML_distance_list"
external to_list : ct -> (int * e * float) list = "char_nonadd_CAML_to_list"


let to_list x = 
    let len = Array.length x.codes in
    let res = ref [] in
    for i = len - 1 downto 0 do
        res := (x.codes.(i), (x.codes.(i), to_int x.data i), 0.) :: !res;
    done;
    !res

(* it's easier to calculate the length of a list on the Caml side *)
external of_list_helper : ct -> (int * e * float) list -> int -> unit
    = "char_nonadd_CAML_of_list_helper"

let sizes = ref All_sets.Integers.empty

let of_list_helper a b =
    let len = List.length a in
    let new_item = make_new_unsafe len b in
    of_list_helper new_item a b;
    new_item

let of_list l =
    let lst = List.stable_sort (fun (c1, _, _) (c2, _, _) -> compare c1 c2) l in
    let codes = Array.of_list (List.map (fun (a, _, _) -> a) lst) in
    { codes = codes; data = of_list_helper lst (List.length l) }

(* We don't store costs in our data structure, and for many set operations, we
   can use the associative list functions if we omit the cost from the triples
   returned by to_list *)
let to_simple_list a =
    List.map (fun (code, elt, _) -> (code, elt)) (to_list a)

let of_simple_list l =
    of_list (List.map (fun (code, elt) -> (code, elt, def_cost)) l)

let parse =
    let make_set (elts, code) =
        let nelts = Array.length elts in
        let true_nelts = 
            Array.fold_left (fun acc x ->
                match x with
                | FileContents.Unordered_Character _, _ -> acc + 1
                | FileContents.Inactive_Character, _ -> acc
                | _ -> 
                    raise
                  (Invalid_argument "Nonadditive characters must be unordered"))
            0 elts 
        in
        let set = make_new_unsafe true_nelts code in
        let rec filler item true_item =
            if item = nelts then ()
            else 
                let (elt, eltcode) = elts.(item) in
                match elt with
                | FileContents.Unordered_Character (elt, _) -> 
                        set_elt_code set true_item eltcode;
                        set_elt set true_item elt;
                        filler (item + 1) (true_item + 1)
                | FileContents.Inactive_Character -> 
                        filler (item + 1) true_item
                | _ -> failwith "Impossible"
        in
        filler 0 0;
        set
    in
    List.map make_set

let list_to_string l =
    "{" ^ String.concat " " (List.map string_of_int l) ^ "}"
(* utility function: make a list from [from] to [tto], and then do a List.map on it
   (in one step) *)
let make_list_iter map from tto =
    let rec mli i =
        if i > tto
        then []
        else (map i) :: (mli (i + 1))
    in
    mli from

let to_string a =
    let strlist =
        make_list_iter
            (fun i -> list_to_string (List.rev (elt_to_list a.data i)))
            0 ((cardinal a.data) - 1) in
    "[" ^ String.concat " " strlist ^ "]"

let e_to_list (_,elt) =
    let rec elt_iter acc elt cntr =
        if elt = 0 then acc 
        else if (elt land 1) <> 0 then 
            elt_iter (cntr :: acc) (elt lsr 1) (cntr + 1)
        else elt_iter acc (elt lsr 1)  (cntr + 1)
    in
    elt_iter [] elt 0

let to_formatter acc attrs c parent d : Xml.xml Sexpr.t list =
    let module T = Xml.Characters in
    let output_character (code, e) (_, cost) =
        let to_sexp = fun x -> 
            let v= `String (Data.to_human_readable d code x) in
            (PXML -[T.value] [v] --)
        in
        let contents = `Set (List.map to_sexp (e_to_list e)) in 
        (PXML -[T.nonadditive]
                    (* Attributes *)
                    ([T.name]=[`String (Data.code_character code d)]) 
                    ([T.cost]=[`Float cost])
                    ([T.definite]=[`Bool  (cost > 0.0)])
                    ([attrs])
                    (* contents *)
                    { contents } --)
    in
    (* TODO CHANGE HTE CODES HERE *)
    let c_ls = to_simple_list c in 
    let dist_list =  
        match parent with 
        | None -> Array.to_list (Array.make (List.length c_ls) (0, 0.) ) 
        | Some parent -> distance_list c.data parent.data
    in 
    let rec produce lst1 lst2 =
        match lst1, lst2 with
        | h1 :: t1, h2 :: t2 -> 
                (output_character h1 h2) :: (produce t1 t2)
        | [], [] -> acc
        | _ -> assert false
    in
    produce c_ls dist_list
 

let empty = make_new 0 0

(*** High-level set operations: *)
(* These typically use to_list and of_list;  they are less efficient than pure C
   implementations, but the hope is that they aren't used often. *)

let add (code, elt) cost set =
    of_list ((code, (code, elt), def_cost) :: (to_list set))

let del code set =
    let list = to_list set in
    let list = List.filter (fun (c, _, _) -> c <> code) list in
    of_list list

(* We do not yet implement colors as components of these sets. *)
let colors _ = []
let codes set =
    List.map (fun (c, _, _) -> c) (to_list set)
let costs set =
    List.map (fun (c, _, cos) -> c, cos) (to_list set)
let get_elt_withcode code set =
    try
        let (_, elt, _) =
            List.find (fun (c, _, _) -> c = code) (to_list set) in
        Some elt
    with
    | Not_found -> None



let rec compare_list compare l1 l2 =
    match (l1, l2) with
    | ([], []) -> 0
    | (x :: xs, y :: ys) ->
          let cmp = compare x y in
          if cmp <> 0
          then cmp
          else compare_list compare xs ys
    | ([], _) -> -1
    | (_, []) -> 1

let compare_union = compare
(* We can use the comparison function encoded in the custom block *)
let compare_data = compare
let compare_codes a b = compare_list compare (codes a) (codes b)



(* This is a silly helper function that we use later to get the list version of
   several characters at once. *)
let pair_map f (a, b) = (f a, f b)

let substitute a b =
    let alist, blist = pair_map to_simple_list (a,b) in
    let new_blist =
        List.map
            (fun (bcode, belt) ->
                 try
                     (bcode, List.assoc bcode alist)
                 with
                 | Not_found -> (bcode, belt))
            blist
    in
    of_simple_list new_blist

let merge a b =
    let alist, blist = pair_map to_simple_list (a,b) in
    let new_alist =
        List.fold_left
            (fun alist (bcode, belt) ->
                 if List.mem_assoc bcode alist
                 then alist
                 else (bcode, belt) :: alist)
            alist
            blist in
    of_simple_list new_alist

let minus a b =
    let alist, blist = pair_map to_simple_list (a,b) in
    let new_alist =
        List.fold_left
            (fun alist (acode, aelt) ->
                 if List.mem_assoc acode blist
                 then alist
                 else (acode, aelt) :: alist)
            []
            alist in
    of_simple_list new_alist

let random f a =
    of_list (List.filter (fun _ -> f ()) (to_list a))

let fold f s a =
    List.fold_left
        (fun accum (code, elt, cost) ->
             f elt accum)
        s
        (to_list a)

let filter f cf a =
    of_list (List.filter f (to_list a))

let f_codes a codes =
    filter (fun (code, _, _) ->
                All_sets.Integers.mem code codes) 
            (fun code -> All_sets.Integers.mem code codes) 
            a
let f_codes_comp a codes =
    filter (fun (code, _, _) ->
                not (All_sets.Integers.mem code codes)) 
            (fun code -> All_sets.Integers.mem code codes) 
    a

let f_colors a col = empty

let iter f a =
    List.iter (fun (code, elt, _) -> f elt code) (to_list a)

let map f a =
    of_list (List.map (fun (code, elt, cost) ->
                           (code, f elt code, cost)) (to_list a))

let is_empty a =
    cardinal a = 0

let of_parser data codes (elts, code) n =
    let make_set elts =
        let nelts = Array.length elts in
        let true_nelts = Array.length elts in
        let set = make_new true_nelts code in
        let rec filler item =
            if item = nelts then ()
            else 
                let (elt, eltcode) = elts.(item) in
                let observed = 
                    match Hashtbl.find data.Data.character_specs eltcode with
                    | Data.Static enc -> enc.Nexus.Parsed.st_observed 
                    | _ -> assert false
                in
                let observed_arr = Array.of_list observed in
                let matcher x =
                    let max = Array.length observed_arr in
                    let rec aux_matcher pos x =
                        if pos = max then failwith "Not found?" 
                        else if observed_arr.(pos) = x then 1 lsl pos
                        else aux_matcher (pos + 1) x
                    in
                    if max > 0 then aux_matcher 0 x
                    else 0
                in
                let elt =
                    let elt = 
                        match elt with
                        | None -> observed
                        | Some (`List x) -> x
                        | Some (`Bits x) -> BitSet.to_list x
                    in
                    List.fold_left (fun acc item -> acc lor (matcher item)) 0 
                    elt
                in
                let () = set_elt set item elt in
                filler (item + 1) 
        in
        filler 0;
        { codes = codes; data = set }
    in
    (make_set elts, code)

let is_potentially_informative elts =
    let states = 
        let make_list_set x = 
            match x with
            | None -> None
            | Some x -> 
                    let x = Nexus.Parsed.static_state_to_list x in
                    Some 
                    (List.fold_left 
                    (fun acc x -> All_sets.Integers.add x acc) 
                    All_sets.Integers.empty x)
        in
        List.fold_left 
        (fun acc x -> 
            match acc with
            | None -> make_list_set x 
            | Some acc ->
                    match make_list_set x with
                    | None -> Some acc
                    | Some x -> Some (All_sets.Integers.inter x acc))
        None elts
    in
    match states with
    | None -> false 
    | Some states ->
            0 = (All_sets.Integers.cardinal states)

let cardinal a = cardinal a.data
(*
let min_cost elts =
    let n_elts = List.length elts in
    let states, graph = build_graph elts in
    let rec compute_cost states remaining_graph =
        if has_single_elts remaining_graph then
            max_int
        else 
            match states with
            | h :: t ->
                min (1 + compute_cost t (pick_state h remaining_graph))
                (1 + compute_cost t (do_not_pick h remaining_graph))
            | [] -> 0
    in
    compute_cost states graph
*)

let extract_elements_present elts = 
    List.fold_right (fun x acc -> match x with None -> acc | Some h -> 
        (Nexus.Parsed.static_state_to_list h) :: acc)
    elts []

let max_possible_cost elts =
    let elts = extract_elements_present elts in
    let counter = Hashtbl.create 97 in
    List.iter (List.iter (fun x -> 
        try 
            let cnt = Hashtbl.find counter x in
            Hashtbl.replace counter x (cnt + 1) 
        with
        | Not_found -> Hashtbl.add counter x 1)) elts;
    let code, max = Hashtbl.fold (fun a b ((x, y) as acc) ->
        if y > b then acc else (a, b)) counter (0, 0) in
    if max = 0 then 0.
    else 
        List.fold_left (fun c a -> 
            if List.mem code a then c else c +. 1.) 0. elts

let min_possible_cost (elts : Nexus.Parsed.static_state list) =
    let elts = extract_elements_present elts in
    let states = 
        let stateset = 
            List.fold_left (fun acc y -> 
                List.fold_left (fun acc x -> All_sets.Integers.add x acc) acc y)
            All_sets.Integers.empty
            elts
        in
        All_sets.Integers.elements stateset
    in
    (* We define a function to filter the elements that are not covered by a
    * particular element *)
    let filter c left = 
        List.filter (fun x -> not (All_sets.Integers.mem c x)) left
    in
    (* A recursive function to solve the problem *)
    let rec optimal_cost codes left = 
        match codes, left with
        | _, [] -> 0
        | [], _ -> (max_int / 2)
        | (h :: t), left ->
                let left' = filter h left in
                min (1 + optimal_cost t left') (optimal_cost t left)
    in
    let res = 
        optimal_cost states (List.map (fun x -> List.fold_left (fun acc x -> 
            All_sets.Integers.add x acc) All_sets.Integers.empty x) elts)
    in
    float_of_int (res - 1)

(* module Test = struct *)
(*     external median_3_debug : t -> t -> t -> t -> t = *)
(*         "char_nonadd_CAML_median_3_debug" *)
(*     let of_simple_list l = *)
(*         of_simple_list (List.map (fun (a, b) -> (a, (a, b))) l) *)
(*     let asdf = *)
(*         if ndebug *)
(*         then begin *)
(*             let a =  of_simple_list [(1, 1); (2, 1); (3, 2)] in *)
(*             let n =  of_simple_list [(1, 3); (2, 1); (3, 3)] in *)
(*             let d1 = of_simple_list [(1, 1); (2, 3); (3, 1)] in *)
(*             let d2 = of_simple_list [(1, 2); (2, 1); (3, 2)] in *)
(*             print_endline (to_string (median_3_debug a n d1 d2)) *)
(*         end *)
(*         else () *)
(* end *)
