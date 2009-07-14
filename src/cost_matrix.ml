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

let () = SadmanOutput.register "Cost_matrix" "$Revision: 2871 $"


exception Illegal_Cm_Format;;
exception Map_Not_Found;;
exception Illegal_List_length;;

type gap_opening = int;;
type cost_model =
    | No_Alignment
    | Linnear
    | Affine of gap_opening 

external init : unit -> unit = "cm_CAML_initialize"

let debug = false
let uselevel = true (* don't forget to set "uselevel" inside alphabet.ml and
sequence.ml *)

(* We override the value of max_int to avoid clashes in 64 bits environments 
and our internal C representation of a matrix with 32 bit integers. *)
let max_int = (Int32.to_int Int32.max_int) lsr 1

let _ =
    init ()

module Two_D = struct
    type m
    external set_gap : m -> int -> unit = "cm_CAML_set_gap"
    (* add following for level *)
    external set_ori_a_sz : m -> int -> unit = "cm_CAML_set_ori_a_sz"
    external get_ori_a_sz : m -> int = "cm_CAML_get_ori_a_sz"
    external set_map_sz : m -> int -> unit = "cm_CAML_set_map_sz"
    external get_map_sz : m -> int = "cm_CAML_get_map_sz"
    external ini_combmap : m -> int -> unit = "cm_CAML_ini_combmap"
    external ini_comb2list: m -> int -> unit = "cm_CAML_ini_comb2list"
    external set_combmap: int -> int -> int -> m -> unit = "cm_CAML_set_combmap"
    external get_combmap: int -> int -> m -> int = "cm_CAML_get_combmap"
(*
    external print_matrix: int -> int -> m -> unit = "cm_CAML_print_matrix"
*)
(*
    external set_comb2list: int -> int -> int -> m -> unit = "cm_CAML_set_comb2list"
    external get_comblist : int -> int -> m -> int = "cm_CAML_get_comblist"
*)  
    (*external get_comblist: int -> m -> 
    (int, Bigarray.int32_elt, Bigarray.c_layout) Bigarray.Array1.t =
        "cm_CAML_get_comblist"*)
    external set_level : m -> int -> unit = "cm_CAML_set_level"
    external get_level : m -> int = "cm_CAML_get_level"
    (* add above for level *)
    external set_all_elements : m -> int -> unit = "cm_CAML_set_all_elements"
    external get_all_elements : m -> int = "cm_CAML_get_all_elements"
    external alphabet_size : m -> int = "cm_CAML_get_a_sz"
    external lcm : m -> int = "cm_CAML_get_lcm"
    external set_lcm : m -> int -> unit = "cm_CAML_set_lcm"
    external set_alphabet_size : int -> m -> unit = "cm_CAML_set_a_sz"
    external gap : m -> int = "cm_CAML_get_gap"
    external combine : m -> int = "cm_CAML_get_combinations"
    external c_affine : m -> int = "cm_CAML_get_affine"
    external c_set_aff : m -> int -> int -> unit = "cm_CAML_set_affine"
    external gap_opening : m -> int = "cm_CAML_get_gap_opening"
    external cost : int -> int -> m -> int = "cm_CAML_get_cost"
    external worst_cost : int -> int -> m -> int = "cm_CAML_get_worst"
    external median : int -> int -> m -> int = "cm_CAML_get_median"
    external set_cost : int -> int -> m -> int -> unit = "cm_CAML_set_cost"
    external set_worst : int -> int -> m -> int -> unit = "cm_CAML_set_worst"
    external set_median : 
        int -> int -> m -> int -> unit = "cm_CAML_set_median"
    external set_metric : m -> unit = "cm_CAML_set_is_metric"
    external c_get_is_metric : m -> int = "cm_CAML_get_is_metric"
    let is_metric m =
        1 = c_get_is_metric m
    external set_prepend : int -> int -> m -> unit = "cm_CAML_set_prepend"
    external set_tail : int -> int -> m -> unit = "cm_CAML_set_tail"
    external get_tail : int -> m -> int = "cm_CAML_get_tail"
    external get_prepend : int -> m -> int = "cm_CAML_get_prepend"
    external create : int -> bool -> int -> int -> int -> int -> int -> m = 
        "cm_CAML_create_bytecode" "cm_CAML_create"
    external clone : m -> m = "cm_CAML_clone"

    let use_combinations = true
    and use_cost_model = Linnear
    and use_gap_opening  = 0


    let print_intlist lst =
        Printf.printf "[ %!"; 
        List.iter (Printf.printf "%d,%!") lst;
        Printf.printf " ] \n%!";
    ;;

    let get_pure_cost_mat cost_mat = 
        let size = alphabet_size cost_mat in      
        let pure_cost_mat = 
            Array.make_matrix (size + 1) (size + 1) (-1)
        in 
        for c1 = 1 to size do
            for c2 = 1 to size do 
                pure_cost_mat.(c1).(c2) <- cost c1 c2 cost_mat
            done 
        done;
        pure_cost_mat;;

    let rec make_file_string ch str =
        try
            let l = ch#read_line in
            make_file_string ch (str ^ " " ^  l);
        with
        | End_of_file -> 
                str;;

    let match_num = 
        Str.regexp " *\\([1-90\\-]+\\) *\\(.*\\)";;

    let rec load_all_integers str l =
        match Str.string_match match_num str 0 with
        | true ->
                let v = Str.matched_group 1 str
                and r = Str.matched_group 2 str in
                let v = Pervasives.int_of_string v in
                load_all_integers r (v :: l);
        | false ->
                List.rev l;;

    let load_file_as_list ch =
        let str = make_file_string ch "" in
        load_all_integers str [];;

    (* Calculating the number of elements in the alphabet. *)
    let calculate_alphabet_size l =
        let s = List.length l in Printf.printf "calculate_alphabet_size = %d\n" s;
        let a_sz = int_of_float (sqrt (float_of_int s)) in
        if (a_sz * a_sz != s) then raise Illegal_Cm_Format
        else a_sz;;

    let affine m = 
        match c_affine m with
        | 0 -> Linnear
        | 1 -> Affine (gap_opening m)
        | _ -> No_Alignment

    let cost_mode_to_int m =
        match m with
        | Linnear -> 0;
        | Affine _ -> 1;
        | No_Alignment -> 2;;

    let cost_position a b a_sz =
        (a lsl a_sz) + b;;

    let median_position a b a_sz =
        cost_position a b a_sz;;

    let output ch m =
        let a_sz = alphabet_size m in
        Pervasives.output_string ch (string_of_int a_sz);
        Pervasives.output_string ch "\n";
        for i = 1 to a_sz do
            for j = 1 to a_sz do
                let v = cost i j m in
                Pervasives.output_string ch (string_of_int v);
                Pervasives.output_string ch " ";
            done;
            Pervasives.output_string ch "\n";
        done;;

    let rec store_input_list_in_cost_matrix_some_combinations m l elt1 elt2 a_sz =
        if(  (elt1 > a_sz) || (elt2 > a_sz) ) then begin
            if( elt1 <= a_sz ) then
            store_input_list_in_cost_matrix_some_combinations m l (elt1+1) 1
            a_sz
            else ()
        end else begin
            match l with
            | h :: t -> 
                if debug then
                    Printf.printf "Store the cost %d %d is %d\n" elt1 elt2 h;
                set_cost elt1 elt2 m h;
                store_input_list_in_cost_matrix_some_combinations m t
                elt1 (elt2 + 1) a_sz;
            | [] -> raise Illegal_Cm_Format;
        end


    let rec store_input_list_in_cost_matrix_all_combinations m l el1 el2 a_sz
    =
        if ((el1 > a_sz) || (el2 > a_sz)) then begin
            if (el2 > a_sz) && (el1 <= a_sz) then 
                store_input_list_in_cost_matrix_all_combinations m l (el1 + 1) 1 a_sz 
            else ()
        end else begin
            match l with
            | h :: t -> 
                    let elt1 = 1 lsl (el1 - 1)
                    and elt2 = 1 lsl (el2 - 1) in
                    if debug then 
                            Printf.printf "Storing the cost %d %d is %d\n" elt1
                            elt2 h;
                    set_cost elt1 elt2 m h;
                    store_input_list_in_cost_matrix_all_combinations m t el1 (el2 + 1) a_sz;
            | [] -> raise Illegal_Cm_Format;
        end

    let rec store_input_list_in_cost_matrix_no_comb m l el1 el2 a_sz
    all_elements =
        if ((el1 > a_sz) || (el2 > a_sz)) then 
            if (el1 <= a_sz) then
                store_input_list_in_cost_matrix_no_comb m l (el1 + 1) 1 a_sz
                all_elements
            else ()
        else 
            match l with
            | h :: t ->
                    let h =
                        if el1 = all_elements || el2 = all_elements then 0
                        else h
                    in
                    if debug then
                        let m = "Setting " ^ string_of_int el1 ^ " and " ^
                        string_of_int el2 ^ " to " ^ string_of_int h in
                        Status.user_message Status.Error m
                    else ();
                    set_cost el1 el2 m h;
                    set_worst el1 el2 m h;
                    store_input_list_in_cost_matrix_no_comb m t el1 (el2 + 1)
                    a_sz all_elements
            | [] -> raise Illegal_Cm_Format

    let store_input_list_in_cost_matrix use_comb m l a_sz all_elements =
        if use_comb then 
            if uselevel then
                store_input_list_in_cost_matrix_some_combinations m l 1 1 a_sz
            else
            store_input_list_in_cost_matrix_all_combinations m l 1 1 a_sz
        else 
            store_input_list_in_cost_matrix_no_comb m l 1 1 a_sz all_elements

    let rec split_integer_in_list_of_bits v bit a_sz l = 
        if bit == a_sz then l
        else begin
            if ((1 lsl bit) land v) != 0 then
                split_integer_in_list_of_bits 
                    v (bit + 1) a_sz ((1 lsl bit) :: l)
            else
                split_integer_in_list_of_bits v (bit + 1) a_sz l
        end;;

    let list_of_bits n max =
        split_integer_in_list_of_bits n 0 max []

    let factorial n =
        let rec func t acc = 
           if (t>1) then   func (t-1) t*acc
           else acc
        in
        func n 1

    let p_m_n m n =
            let rec func t acc =
                if ( t < n ) then func (t+1) (m-t)*acc
                else acc
            in
            func 0 1    

   let calc_number_of_combinations_by_level a_sz level =
        let sum = ref 0 in
        for i = 1 to level do
            let numerator = p_m_n a_sz i in
            let denominator = factorial i in
            let c_m_n = (numerator / denominator) in
            sum := !sum + c_m_n;
        done;
        let num_of_combinations = !sum in
        num_of_combinations

   let calc_num_of_comb_with_gap ori_a_sz level =
        let numerator = p_m_n (ori_a_sz-1) (level-1) in
        let denominator = factorial (level-1) in
        (numerator / denominator) 

    let combcode_includes_gap combcode ori_a_sz level =
        let totalnum = calc_number_of_combinations_by_level ori_a_sz level in
        let num_comb_with_gap = calc_num_of_comb_with_gap ori_a_sz level in
        if ( (totalnum - combcode) < num_comb_with_gap ) then true
        else false
    
    let gap_filter_for_combcode combcode level ori_a_sz=
        let num_comb_with_gap = calc_num_of_comb_with_gap ori_a_sz level in
        let total_num = calc_number_of_combinations_by_level ori_a_sz level in
        let num_comb_without_gap = total_num - num_comb_with_gap in
        if ( (total_num-combcode)< num_comb_with_gap ) then
            (combcode - num_comb_without_gap)
        else
            combcode
     

    (* a,b,i belongs to alpha list, if newcost = cost(a,i)+cost(i,b) is less
    * than oldcost = cost(a,b),
    * set best to i. if they are equal, add i to the "best set" of a and b by
    * "lor" 
    * *)
    let rec test_combinations a l a_sz m best c w =
        let gap = gap m (* gap code... where is the initialization? *)
        and go = gap_opening m in
        match l with
        | b :: t ->
                for i = 0 to a_sz - 1 do
                    let v = 1 lsl i in 
                    (* change: v = i+1 *)
                    let goa =
                        if 1 = c_affine m && v = gap && 
                            (0 <> (a land gap)) && (0 <> (b land gap)) then 
                                go
                        else 0
                    in
                    let fh = cost a v m
                    and sh = cost v b m in
                    let tc = fh + sh + goa in
                    if (tc < !c) then begin
                        c := tc;
                        best := v;
                    end else if (tc == !c) then 
                        begin
                            best := v lor !best;         
                        end
                done; 
                let cab = cost a b m in
                if (cab > !w) then w := cab;
                test_combinations a t a_sz m best c w;
        | [] -> ();;
 
    let clear_duplication_in_list oldlist =
        let oldarray = Array.of_list oldlist in
        let len = Array.length oldarray in
        let newlist = ref [oldarray.(0)] in
        for i = 1 to (len-1) do
            if ( oldarray.(i) <> oldarray.(i-1) ) then 
                newlist := (oldarray.(i))::(!newlist) 
            else ()
        done;
        List.rev (!newlist)
    ;;

    
    let rec comblist_to_combcode lst m =
        let func a b = 
            if(a==b) then a
            else get_combmap a b m
        in
        match lst with
        | h::t -> List.fold_right func lst (List.hd (List.rev lst))
        | _ -> raise(Illegal_List_length)
    ;;
    

(*
    let combcode_to_comblist code m =
        (*Printf.printf "code to list:%!";*)
        let ori_a_sz = get_ori_a_sz m in
        let map_sz = get_map_sz m in
        let rec fillinlist code m acclist =
            if (code<=ori_a_sz) then [code] @ acclist
            else begin
                let code1 = get_comblist code 0 m 
                and code2 = get_comblist code 1 m in
                (*Printf.printf "get code=%d,code1=%d,code2=%d;\n%!"code code1
                code2;*)
               (*if(code==20) then
                    begin
                        Printf.printf "20=[%d,%d]" code1 code2;
                        print_newline();
                        for i = 8 to 15 do
                            let i1 = get_comblist i 0 m and
                            i2 = get_comblist i 1 m in 
                            Printf.printf "%d:%d.%d; %!" i i1 i2;
                        done;
                    end
                else ();*)
                assert(code1<=map_sz); assert(code2<=map_sz);
                let addlist1 = fillinlist code1 m acclist in
                let addlist2 = fillinlist code2 m acclist in
                clear_duplication_in_list addlist1 @ addlist2 @ acclist
            end
        in
        let comblist = fillinlist code m [] in
        let newcomblist = List.sort compare comblist in
        (*Printf.printf "end of code to list\n%!";*)
        clear_duplication_in_list newcomblist
    ;;
*)
  
    let build_maps_comb_and_codelist m ori_a_sz level =
        let comb_to_list = ref All_sets.IntegerMap.empty in
        let list_to_comb = ref All_sets.IntegerListMap.empty in
        let init_arr = Array.init ori_a_sz (fun i -> i+1 ) in
        let init_list = Array.to_list init_arr in      
        let rec all_combinations lst=
            match lst with
            | h :: t -> 
                      let res = all_combinations t in
                      res @ (List.map (fun x -> h :: x) res)
            | [] -> [[]]
        in
        let all_comb_list = 
            match all_combinations init_list with
            | [] :: ((_ :: _) as r) -> r
            | _ -> assert false
        in
        let comb_list_by_level = List.filter 
            (fun x -> ( (List.length x)<=level )) all_comb_list in
        let all_comb_list = comb_list_by_level in 
        let count = ref (ori_a_sz+1) in
        let create_maps lst =
            match lst with
            | [code] -> 
                    let realcode = ori_a_sz+1-code in
                    (*Printf.printf "set_combmap:%d,%d\n%!" realcode realcode;*)
                    
                    set_combmap realcode realcode realcode m;
              (*
                    set_comb2list realcode realcode realcode m;
              *)    
                    comb_to_list := 
                        All_sets.IntegerMap.add (realcode) [realcode] (!comb_to_list);
                    list_to_comb := 
                        All_sets.IntegerListMap.add [realcode] realcode (!list_to_comb);
            | lst ->
                    let codelist =
                        List.fold_left 
                        (fun acc_codelist code ->  
                            let realcode = 
                                if code<=ori_a_sz then (ori_a_sz+1-code)
                                else code
                            in
                            realcode::acc_codelist) [] lst 
                    in
                    let codelist = List.sort compare codelist in
                    let revcodelist =  List.rev(codelist) in
                    (*Printf.printf "[ "; List.iter (Printf.printf "%d,") revcodelist; 
                    Printf.printf " ];";*)
                   
                     let combcode_tl =
                        match revcodelist with
                        | h::t ->
                                (*Printf.printf "get combcode:[ "; List.iter (Printf.printf "%d,") t; 
                                Printf.printf " ];";*)
                                comblist_to_combcode t m
                        | _ -> raise(Illegal_List_length)
                    in
                  
                    set_combmap (List.hd revcodelist) combcode_tl (!count) m;
                    set_combmap combcode_tl (List.hd revcodelist) (!count) m;
                   (*
                    set_comb2list (List.hd revcodelist) combcode_tl (!count) m;
                   *)
                    (*let testres = get_combcode revcodelist m in
                    Printf.printf "get %d,%d=%d\n" (List.hd revcodelist)
                    combcode_tl testres;*)
                    comb_to_list :=
                    All_sets.IntegerMap.add (!count) codelist (!comb_to_list);
                    list_to_comb :=
                    All_sets.IntegerListMap.add codelist (!count) (!list_to_comb);
                    incr count;
        in
        let _ = List.map create_maps all_comb_list in
        !comb_to_list, !list_to_comb
    ;;

    
    let rec browse_combinations a l a_sz m best c w comb_to_list list_to_comb level =
        let gap = 
            if uselevel then a_sz
            else gap m in
        let go = gap_opening m in
        let numerator = p_m_n (a_sz-1) (level-1) in
        let denominator = factorial (level-1) in
        let comb_withgap = (numerator / denominator) in
        let comb_num = calc_number_of_combinations_by_level a_sz level in
       (* Printf.printf "\n browse_combinations combwithgap = %d; %!"
       * comb_withgap;*)
        match l with
        | b :: t ->
                for i = 0 to a_sz - 1 do
                    let v = i+1 in 
                    let goa =
                        if 1 = c_affine m && v = gap && 
                            ((comb_num-a)<comb_withgap) && ((comb_num-b)<comb_withgap) then 
                                go
                        else 0
                    in
                    let fh = cost a v m
                    and sh = cost v b m in
                    let tc = fh + sh + goa in
                   (* Printf.printf " a,b,v; c,tc= %d,%d,%d; %d,%d; " a b v (!c)
                    tc;*)
                    if (tc < !c) then begin
                        c := tc;
                        best := v;
                        (*Printf.printf "  v is better, best=%d; " v;*)
                    end else if (tc == !c) then 
                        begin
                            let oldcomblist = 
                                try  All_sets.IntegerMap.find (!best) comb_to_list 
                                with
                                | Not_found ->  raise Map_Not_Found
                            in  
                            let oldlistlen = List.length oldcomblist in
                            if (oldlistlen < level) then
                                begin
                                    let newcomblist = 
                                        try All_sets.IntegerMap.find v comb_to_list
                                        with
                                        | Not_found -> raise Map_Not_Found
                                    in
                                    (*let newcomblist = combcode_to_comblist v m
                                    * in*)
                                    let newlistlen = List.length newcomblist in
                                    if ((newlistlen+oldlistlen)<=level) then
                                        begin
                                            let newcomblist = List.sort compare
                                                      (newcomblist@oldcomblist) in
                                            let newlist =
                                                clear_duplication_in_list
                                                newcomblist in
                                            let newcomb =
                                                try
                                                All_sets.IntegerListMap.find newlist list_to_comb 
                                                with
                                                |Not_found ->  raise Map_Not_Found
                                            in
                                            (*Printf.printf "newcomb=%d,%!"
                                            * newcomb; *)
                                           (*
                                           let newcomb = comblist_to_combcode newlist m in
                                            *)
                                            (*Printf.printf "%d; %!" newcomb;*)
                                            (*print_intlist newlist;*)
                                            best := newcomb;
                                            (*Printf.printf " combine, best=%d;
                                            * " (!best);*)
                                         end  
                                    else ()
                                end
                            else()
                        end
                    else ()
                done;
                let cab = cost a b m in
                if (cab > !w) then w := cab;
                 browse_combinations a t a_sz m best c w comb_to_list list_to_comb level;
        | [] -> ();;

    let rec test_combinations_by_level li lj a_sz m best c w comb_to_list
    list_to_comb level = 
        match li with
        | h :: t ->
                browse_combinations h lj a_sz m best c w comb_to_list
                list_to_comb level;
                test_combinations_by_level t lj a_sz m best c w comb_to_list
                list_to_comb level
        | [] -> ();;

    let rec test_all_combinations li lj a_sz m best c w = 
        match li with
        | h :: t ->
                test_combinations h lj a_sz m best c w ;
                test_all_combinations t lj a_sz m best c w  
        | [] -> ();;

    let calc_number_of_combinations a_sz = 
        (1 lsl a_sz) - 1;;

    let cleanup m = 
        let gap = gap m in
        match combine m, affine m with
        | 0, _
        | _, Linnear
        | _, No_Alignment -> fun item -> item
        | _ ->
                (fun item ->
                    if gap <> item && (0 <> (gap land item)) then
                        gap
                    else item)

    let xor a b =
        ((a || b) && (not (a && b)));;

    let fill_best_cost_and_median_for_all_combinations_bitwise m a_sz =
        let find_best_combination get_cost l1 l2 =
            let aux_find_best_combination m i (best, med, worst) j =
                let mb1 = get_cost i j in
                let best, median = 
                    if mb1 < best then mb1, i lor j
                    else if mb1 = best then mb1, i lor j lor med
                    else best, med
                in
                best, median, max mb1 worst
            in
            let process l1 l2 acc = 
                List.fold_left (fun acc i ->
                List.fold_left (aux_find_best_combination m i) acc l2)
                acc l1
            in
            process l1 l2 (process l2 l1 (max_int, 0, 0))
        in
        let find_worst_combination get_cost l1 l2 =
            let _, _, w = find_best_combination get_cost l1 l2 in
            w
        in
        let number_of_combinations = calc_number_of_combinations a_sz in
        for i = 1 to number_of_combinations do
            let li = split_integer_in_list_of_bits i 0 a_sz [] in
            for j = 1 to number_of_combinations do 
                let lj = split_integer_in_list_of_bits j 0 a_sz [] in
                if 1 = List.length li && 1 = List.length lj then begin
                    set_median i j m (i lor j);
                    set_worst i j m (cost i j m);
                end else if 0 <> (i land j) then begin
                    let shared = i land j 
                    and worst = 
                        find_worst_combination (fun a b -> cost a b m) 
                        li lj 
                    in
                    set_median i j m shared;
                    set_cost i j m 0;
                    set_worst i j m worst;
                end else begin
                    let cost_f = fun a b -> cost a b m in
                    let best, median, worst = 
                        find_best_combination cost_f li lj 
                    in
                    if debug then
                        Printf.printf "The cost between %d and %d is %d\n" i j best;
                    set_median i j m median;
                    set_cost i j m best;
                    set_worst i j m worst;
                end
            done;
        done
   
    let fill_best_cost_and_median_for_some_combinations m a_sz level  = 
        let get_cost = cost in
        let cleanup = cleanup m in
        (*let number_of_combinations = calc_number_of_combinations a_sz in*)
        let num_of_comb = calc_number_of_combinations_by_level a_sz level in
        Printf.printf "fill for some combinations: num_of_comb =%d\n%!" num_of_comb;
        let (comb_to_list,list_to_comb) =  build_maps_comb_and_codelist m a_sz level in
        for i = 1 to num_of_comb do
            let keylist = [i] in
            let li = 
                try All_sets.IntegerMap.find i comb_to_list with
                | Not_found -> raise Map_Not_Found
            in
            (*let li = combcode_to_comblist i m in*)
            for j = 1 to num_of_comb do 
                let keylist = clear_duplication_in_list(List.sort compare (j::keylist))in
                let lj = 
                    try All_sets.IntegerMap.find j comb_to_list with
                    |Not_found -> raise Map_Not_Found
                in
                (*let lj = combcode_to_comblist j m in*)
                let best = ref 0
                and cost = ref max_int 
                and worst = ref 0 in
                (*Printf.printf "\n lj = "; List.iter (Printf.printf "%d.") lj;*)
                test_combinations_by_level li lj a_sz m best cost worst comb_to_list
                list_to_comb level;
                let _ = 
                    match li, lj with
                    | [_], [_] ->
                            let comb_i_j =
                                try  All_sets.IntegerListMap.find keylist list_to_comb
                                with
                                |Not_found ->  raise Map_Not_Found  
                            in
                            (*Printf.printf "comb_i_j =%d; " comb_i_j;*)
                            (*let comb_i_j = comblist_to_combcode keylist m in*)
                            set_median i j m ( comb_i_j );
                            if debug then
                                Status.user_message Status.Error
                                ("Set the cost beteen " ^ string_of_int i ^ " and " ^
                                string_of_int j ^ " to " ^ string_of_int (get_cost (i)
                                (j) m));
                    | _, _ -> 
                            if debug then
                                Status.user_message Status.Error
                                ("Set the cost beteen " ^ string_of_int
                                i ^ " and " ^ string_of_int (j) ^ " to " ^ string_of_int !cost);
                            set_cost (i) (j) m !cost;
                            set_median (i) (j) m (cleanup !best);
                in
                set_worst (i) (j) m !worst;
            done;
        done;;

    let fill_best_cost_and_median_for_all_combinations m a_sz =
        let get_cost = cost in
        let cleanup = cleanup m in
        let number_of_combinations = calc_number_of_combinations a_sz in
        if debug then
        Printf.printf "alpha size = %d, num of combinations = %d \n%!" a_sz number_of_combinations;
        for i = 1 to number_of_combinations do
            let li = split_integer_in_list_of_bits i 0 a_sz [] in
            for j = 1 to number_of_combinations do 
                let lj = split_integer_in_list_of_bits j 0 a_sz []
                and best = ref 0
                and cost = ref max_int 
                and worst = ref 0 in
                test_all_combinations li lj a_sz m best cost worst;
                let _ = 
                    match li, lj with
                    | [_], [_] ->
                            set_median i j m (i lor j);
                            if debug then
                                Status.user_message Status.Error
                                ("Setting the cost beteen " ^ string_of_int i ^ " and " ^
                                string_of_int j ^ " to " ^ string_of_int (get_cost i
                                j m));
                    | _, _ -> 
                            if debug then
                                Status.user_message Status.Error
                                ("Setting the cost beteen " ^ string_of_int i ^ " and " ^
                                string_of_int j ^ " to " ^ string_of_int !cost);
                            set_cost i j m !cost;
                            set_median i j m (cleanup !best);
                in
                set_worst i j m !worst;
            done;
        done;;

    let fill_medians m a_sz =
        let matrix = Array.make_matrix (a_sz + 1) (a_sz + 1) [] in
        let cleanup = cleanup m in
        for i = 1 to a_sz do
            for j = 1 to a_sz do
                let best = ref max_int 
                and res = ref [] in
                if i = get_all_elements m then begin
                    best := 0;
                    res := [j];
                end else if j = get_all_elements m then begin
                    best := 0;
                    res := [i];
                end else
                    for k = 1 to a_sz do
                        if k = get_all_elements m then ()
                        else
                            let c = (cost k j m) + (cost i k m) in
                            if c < !best then begin
                                best := c;
                                res := [k];
                            end else if c = !best then
                                res := k :: !res
                            else ()
                    done;
                matrix.(i).(j) <- !res;
            done;
        done;
        let matrix = 
            Array.map (fun arr ->
                Array.map (fun x -> 
                    let len = List.length x in
                    if len = 0 then 0
                    else List.nth x (Random.int len)) 
                arr) 
            matrix
        in
        Array.iteri (fun a arr ->
            Array.iteri (fun b median ->
                set_median a b m (cleanup median)) arr) matrix

    let fill_default_prepend_tail m = 
        let a_sz = alphabet_size m 
        and gap = gap m 
        in
        for i = 1 to a_sz do
            set_tail i (cost i gap m) m;
            set_prepend i (cost gap i m) m;
        done

    let set_affine m model = 
        let asz = alphabet_size m in
        let _ = 
            match model with
            | No_Alignment -> c_set_aff m 2 0 
            | Linnear -> c_set_aff m 0 0 
            | Affine go -> c_set_aff m 1 go
        in
        if 1 = combine m then 
            fill_best_cost_and_median_for_all_combinations_bitwise m (lcm m)
        else fill_medians m asz;
        fill_default_prepend_tail m

    let fill_tail_or_prepend f g cost m = 
        (* The alphabet size has to match the current capacity of the matrix *)
        if alphabet_size m < Array.length cost then 
            failwith "The size of the tail or prepend matrix is incorrect"
        else ();
        match combine m with
        | 0 ->
                f 0 0 m;
                for i = 0 to (Array.length cost) - 1 do
                    f (i + 1) cost.(i) m;
                done;
        | _ ->
                (* Check that we are within the boundary of the lcm *)
                if lcm m < Array.length cost then
                    failwith 
                    "The size of the tail or prepend matrix is incorrect"
                else ();
                f 0 0 m;
                (* Initially we store the values we have on the array *)
                for i = 0 to (Array.length cost) - 1 do
                    let pos = 1 lsl i in
                    f pos cost.(i) m;
                done;
                (* Now we go over the array and store each element on it *)
                let cur_min = ref max_int in
                let find_minimum it = 
                    let cc = g it m in
                    if cc < !cur_min then cur_min := cc
                    else ()
                in
                let a_sz = alphabet_size m in
                for i = 1 to a_sz do
                    let lst = split_integer_in_list_of_bits i 0 a_sz [] in
                    List.iter find_minimum lst;
                    f i !cur_min m;
                    cur_min := max_int;
                done

    let fill_tail cost m =
        fill_tail_or_prepend set_tail get_tail cost m

    let fill_prepend cost m =
        fill_tail_or_prepend set_prepend get_prepend cost m

    let is_positive l = 
        let res = ref true 
        and w = Array.length l in
        for i = 0 to w - 1 do
            for j = 0 to w - 1 do
                res := !res && (l.(i).(j) >= 0);
            done;
        done;
        !res

    let list_to_matrix l w =
        let arr = Array.make_matrix w w 0 in
        let _ = 
            List.fold_left (fun (row, col) cost ->
                arr.(row).(col) <- cost;
                if col = w - 1 then
                    row + 1, 0
                else row, col + 1) (0, 0) l
        in
        arr

    let is_symmetric l =
        let w = Array.length l 
        and res = ref true in
        for i = 0 to w - 1 do
            for j = i to w - 1 do
                res := !res && (l.(i).(j) = l.(j).(i))
            done;
        done;
        !res

    let is_identity l =
        let h = Array.length l 
        and res = ref true in
        for i = 0 to h - 1 do
            res := !res && (0 = l.(i).(i));
        done;
        !res

    let is_triangle_inequality l =
        let h = Array.length l 
        and res = ref true in
        for i = 0 to h - 1 do
            for j = 0 to h - 1 do
                for k = 0 to h - 1 do
                    res := !res && (l.(i).(j) <= l.(i).(k) + l.(k).(j))
                done;
            done;
        done;
        !res

    let input_is_metric l w = 
        let arr = list_to_matrix l w in
        is_positive arr && 
        is_symmetric arr && 
        is_triangle_inequality arr &&
        is_identity arr

    let fill_cost_matrix ?(use_comb=true) l a_sz all_elements =   
        let level = 2 in
        let num_comb = calc_number_of_combinations_by_level a_sz level in
        let m = 
            create a_sz use_comb (cost_mode_to_int use_cost_model) 
            use_gap_opening all_elements level num_comb
        in
        if uselevel then
            begin
            (*  set_ori_a_sz m a_sz; 
                set_level m level;
                set_gap m a_sz;
               set_map_sz m num_comb;
               ini_combmap m num_comb;
                ini_comb2list m num_comb;
            *)
                (*Printf.printf "check combmap ini:\n%!";
                print_combmap num_comb num_comb m;*)
            end
        else ();
        let all_elements =
           if uselevel then num_comb 
           else all_elements
        in
        Printf.printf "cost_matrix.ml fill_cost_matrix, a_sz/all_elements = %d/%d%!\n"
        a_sz all_elements;
        store_input_list_in_cost_matrix use_comb m l a_sz all_elements;
        if use_comb then
            if input_is_metric l a_sz then
                let () = set_metric m in
                if uselevel then
                fill_best_cost_and_median_for_some_combinations m a_sz level 
                else
                fill_best_cost_and_median_for_all_combinations m a_sz 
            else
                let _ = 
                    Status.user_message Status.Warning 
                    "You@ are@ loading@ a@ non-metric@ TCM"
                in
                fill_best_cost_and_median_for_all_combinations_bitwise m a_sz
        else fill_medians m a_sz;
        let _= Printf.printf "end of fill cost\n%!" in
        
        fill_default_prepend_tail m;
        m

    let of_channel ?(orientation=false) ?(use_comb = true) all_elements ch =
        Printf.printf "of_channel\n";
        match load_file_as_list ch with
        | [] -> failwith "No Alphabet"
        | l ->
                let w = calculate_alphabet_size l in
                let w, l =
                    if w = all_elements && (not use_comb) then
                        (* We must add a placeholder for the all elements item
                        * *)
                        let rec add_every cnt lst =
                            match lst with
                            | [] -> []
                            | h :: t -> 
                                    if 0 = cnt mod w then 
                                        1 :: h :: (add_every (cnt + 1) t)
                                    else h :: (add_every (cnt + 1) t)
                        in
                        w + 1, add_every 1 l
                    else w, l
                in
                let m = 
                    match orientation with 
                    | false ->
                          fill_cost_matrix ~use_comb:use_comb l w all_elements 
                    | true ->
                          let l_arr = Array.of_list l in          
                          let w2 = w * 2 - 1 in  
                          let l2 = ref [] in 
                          for code1 = 1 to w2 do 
                              for code2 = 1 to w2 do 
                                  let index = 
                                      ((((code1 + 1) / 2) - 1)  * w)  + 
                                      (((code2 + 1) / 2) - 1) in 
                                  l2 := l_arr.(index)::!l2 
                              done;  
                          done; 
                          let l2 = List.rev !l2 in 
                          fill_cost_matrix ~use_comb:use_comb l2 w2 all_elements
                in  
                m;;


    let of_list ?(use_comb=true) l all_elements =
        (* This function assumes that the list is a square matrix, list of
        * lists, all of the same size *)
        let w = List.length l in
        let l = List.flatten l in
        fill_cost_matrix ~use_comb:use_comb l w all_elements 

    let of_channel_nocomb ?(orientation=false) ch =
        Printf.printf "cost_matrix.ml of_channel_nocomb\n %!";
        of_channel ~orientation:orientation ~use_comb:false ch

    let of_list_nocomb l all_elements =
        (* This function assumes that the list is a square matrix, list of
        * lists, all of the same size *)
        let w = List.length l in
        let l = List.flatten l in
        fill_cost_matrix ~use_comb:false l w all_elements

    let default = 
        of_list [
            [0;1;1;1;2];
            [1;0;1;1;2];
            [1;1;0;1;2];
            [1;1;1;0;2];
            [2;2;2;2;0];
        ] 31

    let default_nucleotides = default

    let default_aminoacids =
        of_list_nocomb [
            [0; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 2];
            [1; 0; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 2];
            [1; 1; 0; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 2];
            [1; 1; 1; 0; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 2];
            [1; 1; 1; 1; 0; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 2];
            [1; 1; 1; 1; 1; 0; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 2];
            [1; 1; 1; 1; 1; 1; 0; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 2];
            [1; 1; 1; 1; 1; 1; 1; 0; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 2];
            [1; 1; 1; 1; 1; 1; 1; 1; 0; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 2];
            [1; 1; 1; 1; 1; 1; 1; 1; 1; 0; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 2];
            [1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 0; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 2];
            [1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 0; 1; 1; 1; 1; 1; 1; 1; 1; 1; 2];
            [1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 0; 1; 1; 1; 1; 1; 1; 1; 1; 2];
            [1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 0; 1; 1; 1; 1; 1; 1; 1; 2];
            [1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 0; 1; 1; 1; 1; 1; 1; 2];
            [1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 0; 1; 1; 1; 1; 1; 2];
            [1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 0; 1; 1; 1; 1; 2];
            [1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 0; 1; 1; 1; 2];
            [1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 0; 1; 1; 2];
            [1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 0; 1; 2];
            [1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 0; 2];
            [2; 2; 2; 2; 2; 2; 2; 2; 2; 2; 2; 2; 2; 2; 2; 2; 2; 2; 2; 2; 2; 0];
        ] 21

    let of_transformations_and_gaps use_combinations alph_size trans gaps all_elements =
        let list_with_zero_in_position pos =
            Array.to_list
            (Array.init alph_size (fun x ->
                if x = pos then 0
                else if x = alph_size - 1 then gaps
                else if pos = alph_size - 1 then gaps
                else trans))
        in
        of_list ~use_comb:use_combinations (Array.to_list 
        (Array.init alph_size (fun x ->
            list_with_zero_in_position x))) all_elements

    let perturbe cm sev prob =
        let cm = clone cm in
        let lim = if (0 != combine cm) then lcm cm else alphabet_size cm in
        for i = 1 to lim do
            for j = 1 to lim do
                let b = Random.int 101 in
                if b <= prob then set_cost i j cm (sev * cost i j cm);
            done;
        done;
        if (0 != combine cm) then
            fill_best_cost_and_median_for_all_combinations_bitwise cm (lcm cm);
        fill_default_prepend_tail cm;
        cm;;

    let copy_matrix src tgt src_sz tgt_sz = 
        let w = calc_number_of_combinations src_sz in
        for i = 1 to w do
            for j = 1 to w do
                set_cost i j tgt (cost i j src);
                set_median i j tgt (median i j src);
            done;
        done;;
    let get_closest cm a b =
        assert(a>0); assert(b>0);

           (*     Printf.printf "check maps:\n%!";
                for i = 8 to 36 do
                    let reslist = combcode_to_comblist i cm in
                    Printf.printf "i=%d,%!" i;
                    (*let rescode = comblist_to_combcode ([i]@[i+1]) m in*)
                    print_intlist reslist;
                done;
           *)
        (*let _ = Printf.printf "get_closest: a,b=%d,%d; %!" a b in*)
        if 0 = combine cm then b (* We only have one element in b *)
        else
            (* get rid of "gap" in b is 'x or z or...gap' *)
            let gap = gap cm in
            let b =
                if uselevel then
                    let level = get_level cm in
                     let ori_a_sz = get_ori_a_sz cm in 
                    if ( a=gap || b=gap ) then b
                    else if ( combcode_includes_gap a ori_a_sz level) && (
                        combcode_includes_gap b ori_a_sz level) then gap
                    else 
                        gap_filter_for_combcode b level ori_a_sz
                else
                    if a = gap || b = gap then b
                    else if (0 <> a land gap) && (0 <> b land gap) then gap
                    else b land (lnot gap)
            in
            let listofcode =
                if uselevel then 
                        let level = get_level cm in
                        let ori_a_sz = get_ori_a_sz cm in 
                        let (comb_to_list,_) =  
                            build_maps_comb_and_codelist cm ori_a_sz level in 
                        try
                        All_sets.IntegerMap.find b comb_to_list with
                        | Not_found -> raise (Map_Not_Found)
                    (*Printf.printf "check list= %!"; print_intlist checklist;*)
                    (*combcode_to_comblist b cm*) 
                else
                    list_of_bits b (alphabet_size cm) 
            in
            match listofcode(*list_of_bits b (alphabet_size cm)*) with
            | [] -> failwith "~ No bits on?"
            | bits ->
                    let bits = List.rev bits in
                    let best, _ = 
                        List.fold_left (fun ((_, cur_cost) as acc) x ->
                        let nc = cost a x cm in
                        if nc < cur_cost then (x, nc)
                        else acc) (a, max_int) bits
                    in
                    best
end;;

module Three_D = struct
    type m
    external set_gap : m -> int -> unit = "cm_CAML_set_gap_3d"
    external c_set_aff :
        m -> int -> int -> unit = "cm_CAML_set_affine_3d"
    external set_all_elements : m -> int -> unit = "cm_CAML_set_all_elements_3d"
    external get_all_elements : m -> int = "cm_CAML_get_all_elements_3d"
    external alphabet_size : m -> int = "cm_CAML_get_a_sz_3d"
    external set_alphabet_size : int -> m -> unit = "cm_CAML_set_a_sz_3d"
    external lcm : m -> int = "cm_CAML_get_lcm_3d"
    external set_lcm : m -> int -> unit = "cm_CAML_set_lcm_3d"
    external gap : m -> int = "cm_CAML_get_gap_3d"
    external combine : m -> int = "cm_CAML_get_combinations_3d"
    external c_affine : m -> int = "cm_CAML_get_affine_3d"
    external gap_opening : m -> int = "cm_CAML_get_gap_opening_3d"
    external cost : int -> int -> int -> m -> int = "cm_CAML_get_cost_3d"
    external median : int -> int -> int -> m -> int = "cm_CAML_get_median_3d"
    external set_cost :
        int -> int -> int -> m -> int -> unit = "cm_CAML_set_cost_3d"
    external set_median : 
        int -> int -> int -> m -> int -> unit = "cm_CAML_set_median_3d"
    external clone : m -> m = "cm_CAML_clone_3d"
    (* [create a_sz comb aff go ge dim] creates a new cost matrix with
    alphabet size a_sz, considering all the possible combinations
    of members of the alphabet (total efective size 2^a_sz - 1),
    if comb is true with affine gap cost model if aff is not 0, in
    which case will use gap opening cost go and extension gap ge. The
    gap is represented as the next available integer depending on
    the a_sz. IF dim is true the matrix will be three dimensional,
    otherwise it will be two dimensional. *)
    external create : 
        int -> bool -> int -> int -> int -> int -> m = "cm_CAML_create_3d_bc" "cm_CAML_create_3d"
    external make_three_dim : Two_D.m -> m = "cm_CAML_clone_to_3d"

    let use_combinations = true
    and use_cost_model = Linnear
    and use_gap_opening  = 0

    let set_affine m model = 
        match model with
        | No_Alignment -> c_set_aff m 2 0
        | Linnear -> 
                c_set_aff m 0 0 
        | Affine go -> 
                c_set_aff m 1 go


    let affine m = 
        c_affine m;;

    let int_to_affine m = 
        match c_affine m with
        | 0 -> Linnear;
        | 1 -> Affine (gap_opening m)
        | _ -> No_Alignment;;

    let cost_position a b c a_sz =
        (((a lsl a_sz) + b) lsl a_sz) + c;;

    let median_position a b c a_sz = 
        cost_position a b c a_sz;;

    let output ch m = 
        let a_sz = lcm m in
        let w = (1 lsl a_sz) - 1 in
        for i = 1 to w do
            for j = 1 to w do
                for k = 1 to w do
                    let v = cost i j k m in
                    Pervasives.output_string ch (string_of_int i);
                    Pervasives.output_char ch ':';
                    Pervasives.output_string ch (string_of_int j);
                    Pervasives.output_char ch ':';
                    Pervasives.output_string ch (string_of_int k);
                    Pervasives.output_char ch ':';
                    Pervasives.output_string ch (string_of_int v);
                    Pervasives.output_string ch "\n";
                done;
            done;
        done;;

    let of_two_dim_no_comb nm mat = 
        let alph = (alphabet_size nm) - 1 in
        for i = 1 to alph do
            for j = 1 to alph do
                for k = 1 to alph do
                    for l = 1 to alph do
                        let cost = 
                            (Two_D.cost l i mat) + (Two_D.cost l j mat) +
                            (Two_D.cost l k mat) 
                        and old_cost = cost i j k nm in
                        if (cost < old_cost) then begin
                            set_cost i j k nm cost;
                            set_median i j k nm l;
                        end;
                    done;
                done;
            done;
        done;;

    let of_two_dim_comb nm m = 
        let alph = Two_D.alphabet_size m
        and gap = Two_D.gap m 
        and lcm = Two_D.lcm m in
        let max = 1 lsl lcm in
        let rec pick_bit cur item =
            assert (
                if cur < max then true
                else 
                    let _ = Printf.printf 
                    "Max is %d and lcm is %d while cur is %d\n%!" max
                        lcm cur in
                    false);
            if 0 <> cur land item then cur
            else pick_bit (cur lsl 1) item
        in
        let is_metric = Two_D.is_metric m in
        let has_shared_bit x y =
            if 0 <> x land y then 1 else 0 
        in
        for i = 1 to alph do
            for j = 1 to alph do
                for k = 1 to alph do
                    set_cost i j k nm max_int;
                    for l = 0 to lcm - 1 do
                        let inter = 1 lsl l in
                        let cost = 
                            if is_metric || 
                            (2 <= ((has_shared_bit i inter) + 
                            (has_shared_bit j inter) + (has_shared_bit k
                            inter))) || (inter <> gap) then
                                (Two_D.cost inter i m) + 
                                (Two_D.cost inter j m) +
                                (Two_D.cost inter k m) 
                            else max_int
                        and old_cost = cost i j k nm in
                        if (cost < old_cost) then begin
                            set_cost i j k nm cost;
                            set_median i j k nm inter;
                        end else if (cost == old_cost) then begin
                            let v = median i j k nm in
                            set_median i j k nm (v lor inter);
                        end;
                    done;
                    (* We pick only one median ammong all options *)
                    set_median i j k nm (pick_bit 1 (median i j k nm));
                done;
            done;
        done;;

    let of_two_dim m = 
        let nm = make_three_dim m in
        match combine nm with
        | 0 -> 
                of_two_dim_no_comb nm m;
                nm;
        | _ -> 
                of_two_dim_comb nm m;
                nm;;

    let perturbe cm sev prob =
        let cm = clone cm in
        let lim = if (0 != combine cm) then lcm cm else alphabet_size cm in
        for i = 1 to lim do
            for j = 1 to lim do
                for k = 1 to lim do
                    let b = Random.int 101 in
                    if b >= prob then set_cost i j k cm (sev * cost i j k cm);
                done;
            done;
        done;
        cm

    let default = of_two_dim Two_D.default
    let default_nucleotides = default
    (* This is lazy because the POY initialization would take too long *)
    let default_aminoacids = 
        Lazy.lazy_from_fun (
            fun () -> 
                let st = Status.create 
                    "Calculating@ Three@ Dimansional@ Aminoacid@ Cost@ Matrix"
                    None ""
                in
                let res = of_two_dim Two_D.default_aminoacids in
                Status.finished st;
                res
            )

end;;
