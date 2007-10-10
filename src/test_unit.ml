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

let () = SadmanOutput.register "Test_unit" "$Revision: 1644 $"

(* $Id: test_unit.ml 1644 2007-02-14 19:05:47Z andres $ *)

(** test_unit.ml *)

(** Apply a function [f] a given number of times ([n]) *)
let repeat f n = for i = 1 to n do f () done

(** [assert_bool string bool] prints the message in [string] and then a result
    code corresponding to [bool] *)
let assert_bool str b =
    print_endline
        (str ^ ": " ^ (if b then "yes!" else "FAILED!!"))

(** [list_collapse_pairs f list] takes a function which collapses pairs of
    elements into one element, and it applies this function sequentially on
    every pair of elements in [list].  The return result is the folded list.  If
    there is an odd number of elements in [list], the last element will be
    untouched, and will be in the result list. *)
let rec list_collapse_pairs f list =
    match list with
    | a :: b :: list -> (f a b) :: (list_collapse_pairs f list)
    | a -> a

(** This function calls [list_collapse_pairs] recursively until only one element
    is left in the list. *)
let rec collapse_to_one f list =
    match list with
    | [] -> raise (Invalid_argument "Empty list")
    | [a] -> a
    | list -> collapse_to_one f (list_collapse_pairs f list)


let randomize_list l =
    let a = Array.of_list l in
    let len = Array.length a in
    for i = 0 to len - 1 do
        let j = Random.int len in
        if i <> j
        then begin
            let temp = a.(i) in
            a.(i) <- a.(j);
            a.(j) <- temp
        end
    done;
    Array.to_list a

module Test_character(Char : Character.CHARACTER) = struct

    let maker gen =
        let gen =
            match gen with
            | None -> Char.rand_gen ()
            | Some g -> g
        in
        fun () -> Char.make_rand gen

    let assert_bool_1 s b c =
        assert_bool s b;
        if not b
        then print_endline (Char.to_string c)

    let assert_bool_2 s b c1 c2 =
        assert_bool s b;
        if not b
        then (print_endline ("Char. 1:" ^ Char.to_string c1);
              print_endline ("Char. 2:" ^ Char.to_string c2))
            
    let test_median_noparent gen =
        let make = maker gen in
        let c1 = make () in
        let c2 = make () in
        let median = 
            Char.median None c1 c2 in
        let dist1 = Char.distance c1 c2 in
        let dist2 = Char.distance c1 median in
        let dist3 = Char.distance median c2 in
        (* positive *)
        assert_bool_2 "Distance non-negative (c1, c2)" (dist1 >= 0.) c1 c2;
        assert_bool_2 "Distance non-negative (c1, median)" (dist2 >= 0.) c1 median;
        assert_bool_2 "Distance non-negative (median, c2)" (dist3 >= 0.) median c2;
        (* median should be closer *)
        assert_bool "Median closer to c1" (dist2 <= dist1);
        assert_bool "Median closer to c2" (dist3 <= dist1);
        ()

    let print_median_of c1 c2 =
        let median =
            Char.median None c1 c2 in
        
        print_string "Character 1: ";
        print_string (Char.to_string c1);
        print_newline ();
        print_string "Character 2: ";
        print_string (Char.to_string c2);
        print_newline ();
        print_string "Median: ";
        print_string (Char.to_string median);
        print_newline ();
        print_newline ();
        median

    let print_median gen =
        let make = maker gen in
        let c1 = make () in
        let c2 = make () in
        print_median_of c1 c2

    let print_median_4 gen =
        let make = maker gen in
        let c1 = make () in
        let c2 = make () in
        let c3 = make () in
        let c4 = make () in
        let med1 = print_median_of c1 c2 in
        let med2 = print_median_of c3 c4 in
        print_median_of med1 med2

    let print_equal_diff gen =
        let make = maker gen in
        let c1 = make () in
        let c2 = make () in
        let c3 = Char.make_rand (Char.rand_gen ()) in
        print_endline "A character and itself:";
        assert_bool_1 "...compare codes is zero" (0 = Char.compare_codes c1 c1) c1;
        assert_bool_1 "...compare data is zero" (0 = Char.compare_data c1 c1) c1;
        assert_bool_1 "...distance is non-negative" (0. <= Char.distance c1 c1) c1;
        print_endline "Compared to other from same maker:";
        assert_bool_2 "...compare_codes is zero" (0 = Char.compare_codes c1 c2) c1 c2;
        assert_bool_2 "...compare_data is nonzero" (0 <> Char.compare_data c1 c2) c1 c2;
        assert_bool_2 "...distance is positive" (0. < Char.distance c1 c2) c1 c2;
        print_endline "Compared to other from other maker:";
        assert_bool_2 "...compare_codes is nonzero" (0 <> Char.compare_codes c1 c3) c1 c3;
        assert_bool_2 "...compare_data is nonzero" (0 <> Char.compare_data c1
                                                        c3) c1 c3

    let print_marshal gen =
        let make = maker gen in
        let c1 = make () in
        let str = Marshal.to_string c1 [] in
        try
            let c2 = Marshal.from_string str 0 in
            assert_bool_2 "Marshalling is consistent" (0 = Char.compare_data c1 c2)
                c1 c2
        with
        | Failure "Memory error." ->
              assert_bool_1 "Marshalling is consistent (no mem error)" false c1

    let test_all str gen =
        print_endline ("*** Testing character functions of " ^ str);
        print_equal_diff gen;
(*         print_marshal gen; *)
        ignore (print_median_4 gen)



    module OrderedChar = struct
        type t = Char.t
        let compare = Char.compare_data
    end

    module CharSet = Set.Make(OrderedChar)

    let take_n_medians n gen =
        let p_new = 0.1 in
        let make = maker gen in

        let rand_elt set =
            if Random.float 1. < p_new
            then make ()
            else
                let (elt0, elts) =
                    match CharSet.elements set with
                    | a :: b -> (a, b)
                    | _ -> failwith "Nothing in the set" in
                let (elt, _) = List.fold_left
                    (fun (elt, num) newelt ->
                            if Random.int num = 0
                            then (newelt, num + 1)
                            else (elt, num + 1))
                    (elt0, 2) elts in
                elt
        in
        let rand_elt_not set elt = rand_elt (CharSet.remove elt set) in
        let initset = CharSet.empty in
        let initset' = CharSet.add (make ()) initset in
        let initset'' = CharSet.add (make ()) initset' in
        let set = ref initset'' in

        let setsize = ref 2 in
        let set_max_size = 100 in

        for i = 0 to n - 1 do
            let c1 = rand_elt !set in
            let c2 = rand_elt_not !set c1 in
            let median = Char.median None c1 c2 in
            let set' = CharSet.add median !set in
            set := set';
            incr setsize;

            if !setsize > set_max_size
            then begin
                
                let f = let prev = ref false in fun _ -> prev := not !prev;
                    !prev in
                set := CharSet.filter f !set;
                setsize := !setsize / 2
            end
(*             print_endline (Char.to_string median) *)
        done
        
end

module Test_characterset(Char : Character.CharacterSet) = struct
    module AsChar = Test_character(Char)

    let maker gen =
        let gen =
            match gen with
            | None -> Char.rand_gen ()
            | Some g -> g
        in
        fun () -> Char.make_rand gen

    let assert_bool_1 s b c =
        assert_bool s b;
        if not b
        then print_endline (Char.to_string c)

    let assert_bool_2 s b c1 c2 =
        assert_bool s b;
        if not b
        then (print_endline ("Char. 1:" ^ Char.to_string c1);
              print_endline ("Char. 2:" ^ Char.to_string c2))

    let get_rand_elt c =
        let codes = Char.codes c in
        let size = List.length codes in
        let n = Random.int size in
        let rand_code = List.nth codes n in
        match Char.get_elt_withcode rand_code c with
        | Some e -> e
        | None -> failwith "Code inconsistency"

    let test_empty_cardinality_add_del gen =
        let c0 = Char.empty in
        assert_bool_1 "Empty set cardinality zero" (0 = Char.cardinal c0) c0;
        assert_bool_1 "Empty set equal to self" (Char.empty = c0) c0;
        assert_bool_1 "Empty set is empty" (Char.is_empty c0) c0;
        let relt = get_rand_elt (maker gen ()) in
        let relt_code = Char.elt_code relt in
        let c1 = Char.add relt 1. c0 in
        assert_bool_1 "|Empty plus one| = 1" (1 = Char.cardinal c1) c1;
        assert_bool_1 "...has correct code"
            (match Char.get_elt_withcode relt_code c1 with
             | Some _ -> true
             | None -> false) c1;
        assert_bool_1 "...is not empty" (not (Char.is_empty c1)) c1;
        let c2 = Char.del relt_code c1 in
        assert_bool_2 "Deletion brings us back to empty" (0 = Char.compare_data c0 c2) c0 c2;
        assert_bool_1 "...is empty" (Char.is_empty c2) c2

    let test_to_of_list gen =
        let make = maker gen in
        let c0 = make () in
        let list = Char.to_list c0 in
        let c1 = Char.of_list list in
        assert_bool_2 "to_list and of_list consistency" (0 = Char.compare_data c0 c1) c0 c1

    let test_distance_list_consistent gen =
        let make = maker gen in
        let c0 = make () in
        let c1 = make () in
        let dist = Char.distance c0 c1 in
        let distlist = Char.distance_list c0 c1 in
        let distlistsum =
            List.fold_left (fun sum (_, dist) -> dist +. sum) 0. distlist in
        assert_bool_2 "Distance list is consistent with distance"
            (distlistsum = dist) c0 c1

    let test_all str gen =
        AsChar.test_all str gen;
        print_endline ("*** Testing character set functions of " ^ str);
        test_empty_cardinality_add_del gen;
        test_to_of_list gen;
        test_distance_list_consistent gen

    let take_n_medians = AsChar.take_n_medians

        (* Still to check: substitute, merge, minus, random, fold, f_codes,
           iter, map *)
end

module Test_Nonadditive = Test_character(Char_nonadd)
module Test_Additive = Test_character(Char_add)
module Test_c8 = Test_characterset(NonaddCS8)
module Test_c16 = Test_characterset(NonaddCS16)
module Test_c32 = Test_characterset(NonaddCS32)
module Test_sankoff = Test_character(SankCS)

(* module Set_of_Nonadditive = Charset.Of(NonaddCS8) *)
(* module Set_of_Additive = Charset.Of(Char_add) *)
(* module Set_of_Set_of_Additive = Charset.Of(Set_of_Additive) *)
(* module Set_of_Set_of_Set_of_Additive = Charset.Of(Set_of_Set_of_Additive) *)

(* module Test_NSet = Test_characterset(Set_of_Nonadditive) *)
(* module Test_ASet = Test_characterset(Set_of_Additive) *)
(* module Test_Set_3 = Test_characterset(Set_of_Set_of_Set_of_Additive) *)

let _ = Test_c8.test_all "nonadd8" (Some (400, 8, 1001))
let _ = Test_c16.test_all "nonadd16" (Some (400, 16, 1003))
let _ = Test_c32.test_all "nonadd32" (Some (400, 31, 1005))

(* let _ = Test_c.take_n_medians 100000 None *)
(* let _ = repeat (fun () -> Test_c.test_all None *)
let _ = Test_Additive.test_all "additive" None
(* let _ = Test_NSet.test_all None *)
let _ = Test_sankoff.test_all "sankoff" None
(* let _ = Test_Set_3.test_all None *)















(* let _ = Asynch.register_functions () *)
(* let _ = *)
(*     try *)
(*         Comm.init Sys.argv Comm.Fault_Tolerant Comm.Asynch *)
(*     with *)
(*     | e -> raise e *)


(* let _ = Random.init 23 *)
(* (\* let _ = Test_NSet.print_median_4 (Some (Set_of_Nonadditive.rand_gen_parallel ())) *\) *)
(* let _ = Test_ASet.print_median_4 (Some (Set_of_Additive.rand_gen_parallel ())) *)
(* (\* let _ = Test_Set_3.print_median_4 (Some (Set_of_Set_of_Set_of_Additive.rand_gen_parallel ())) *\) *)
