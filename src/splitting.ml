let set_cost matrix i j cost = 
    try
        let table = Hashtbl.find matrix i in
        Hashtbl.replace table j cost
    with
    | Not_found ->
            let table = Hashtbl.create 59 in
            Hashtbl.replace table j cost;
            Hashtbl.add matrix i table

let get_cost matrix i j =
    Hashtbl.find (Hashtbl.find matrix i) j

let has_cost matrix i j =
    try Hashtbl.mem (Hashtbl.find matrix i) j with
    | Not_found -> false

let ( --> ) a b = b a

let lazy_alignment cm a b = 
    let lena = Sequence.length a 
    and lenb = Sequence.length b 
    and gap = Cost_matrix.Two_D.gap cm in
    let cost a b = Cost_matrix.Two_D.cost a b cm in
    let matrix = Hashtbl.create lena in
    let set_cost = set_cost matrix
    and get_cost = get_cost matrix
    and has_cost = has_cost matrix in
    let table_to_fill = Hashtbl.create lena in
    let rec fill_alignment current_cost current_to_fill next_row =
        match current_to_fill with
        | (i, j) :: tl ->
                if lena > i && lenb > j && not (has_cost i j) then begin
                    let a = Sequence.get a i
                    and b = Sequence.get b j in
                    let cost = 
                        min
                        (min 
                        (try ((get_cost (i - 1) j) + (cost gap a))
                        with Not_found -> max_int / 3)
                        (try ((get_cost i (j - 1)) + (cost gap b))
                        with Not_found -> max_int / 3)
                        )
                        ((try (get_cost (i - 1) (j - 1))
                        with Not_found -> max_int / 3) + (cost a b))
                    in
                    set_cost i j cost;
                    if i = lena - 1 && j = lenb - 1 then cost
                    else if cost = current_cost then 
                        let add_if u v x =
                            if not (has_cost u v) then (u, v) :: x 
                            else x
                        in
                        let next_row =
                            next_row --> add_if (i + 1) j -->
                                add_if i (j + 1) --> add_if (i + 1) (j + 1)
                        in
                        fill_alignment 
                        current_cost tl
                        next_row
                    else begin
                        if not (has_cost (i + 1) j) then
                            Hashtbl.add table_to_fill cost (i + 1, j);
                        if not (has_cost i (j + 1)) then
                            Hashtbl.add table_to_fill cost (i, j + 1);
                        if not (has_cost (i + 1) (j + 1)) then
                            Hashtbl.add table_to_fill cost (i + 1, j + 1);
                        fill_alignment current_cost tl next_row
                    end;
                end else fill_alignment current_cost tl next_row
        | [] -> 
                match next_row with
                | [] ->
                        let rec find_next cost = 
                            if Hashtbl.mem table_to_fill cost then
                                let list = Hashtbl.find_all table_to_fill cost in
                                let () = Hashtbl.remove table_to_fill cost in
                                fill_alignment cost (List.rev list) []
                            else find_next (cost + 1)
                        in
                        find_next (current_cost + 1)
                | _ -> fill_alignment current_cost (List.rev next_row) []
    in
    set_cost 0 0 0;
    for i = 1 to lena - 1 do
        let a = Sequence.get a i in
        set_cost i 0 
        ((get_cost (i - 1) 0) + (cost gap a))
    done;
    for j = 1 to lenb - 1 do
        let b = Sequence.get b j in
        set_cost 0 j
        ((get_cost 0 (j - 1)) + (cost gap b))
    done;
    let _ = fill_alignment 0 [(1,1)] [] in
    (*
    Printf.printf "The cost is %d\n%!" res;
    let res = 
        Array.init lena 
        (fun a ->
            Array.init lenb
            (fun b ->
                if has_cost a b then
                    let () = incr cnt in
                    get_cost a b
                else (max_int / 3)))
    in
    Printf.printf "Filled %d of %d\n%!" !cnt (lena * lenb);
    *)
    matrix


let align cm a b =
    let lena = Sequence.length a 
    and lenb = Sequence.length b 
    and gap = Cost_matrix.Two_D.gap cm in
    let cost a b = Cost_matrix.Two_D.cost a b cm in
    let matrix = Array.make_matrix lena lenb 0 in
    for i = 1 to lena - 1 do
        let a = Sequence.get a i in
        matrix.(i).(0) <- matrix.(i - 1).(0) + cost gap a;
    done;
    for j = 1 to lenb - 1 do
        let b = Sequence.get b j in
        matrix.(0).(j) <- matrix.(0).(j - 1) + cost gap b;
    done;
    for i = 1 to lena - 1 do
        for j = 1 to lenb - 1 do
            let a = Sequence.get a i
            and b = Sequence.get b j in
            matrix.(i).(j) <- 
                min
                    (min 
                        (matrix.(i - 1).(j) + (cost gap a))
                        (matrix.(i).(j - 1) + (cost gap b)))
                    (matrix.(i - 1).(j - 1) + (cost a b));
        done;
    done;
    matrix


let align cm a b = 
    (*
    let lena = Sequence.length a 
    and lenb = Sequence.length b in
    let m1 = align cm a b in
    *)
    let m2 = lazy_alignment cm a b in
    (*
    assert (m1.(lena - 1).(lenb - 1) = m2.(lena - 1).(lenb - 1));
    *)
    m2


type routes = 
    | Empty 
    | Single of (routes * int * routes * int * int * int) (* The number of paths *)
    | Multiple of (routes list * int * int * int) (* The number of paths *)

let produce_matrices cm a b =
    let m1 = align cm a b 
    and m2 = 
        let reverse a = 
            let gap = Cost_matrix.Two_D.gap cm in
            let lena = (Sequence.length a) - 1 in
            Sequence.init (fun pos ->
                if pos = 0 then gap
                else Sequence.get a (lena + 1 - pos)) (lena + 1)
        in
        let a = reverse a 
        and b = reverse b in
        align cm a b 
    in
    let lena = (Sequence.length a - 1) 
    and lenb = (Sequence.length b - 1) in
    let optimal = get_cost m1 lena lenb in
    let table = Hashtbl.create (lena + 1) in
    Hashtbl.iter (fun apos bpart ->
        Hashtbl.iter (fun bpos cost ->
            let cost = 
                try
                (if apos = 0 || bpos = 0 then 1
                else
                    let a = Sequence.get a apos
                    and b = Sequence.get b bpos in
                    optimal -
                    ((Cost_matrix.Two_D.cost a b cm) +
                    (get_cost m1 (apos - 1) (bpos - 1)) + 
                    (get_cost m2 (lena - apos) (lenb - bpos))))
                with Not_found -> -10
            in
            if cost = 0 then set_cost table apos bpos (cost, Empty)) bpart) m1;
    table
(*
    Array.init (lena + 1)
    (fun apos ->
        Array.init (lenb + 1)
        (fun bpos ->
            (* Tuples consisting of the cost, and the list of intermediate
            * elements that are being mapped, for now, as we have no
            * intermediate elements, the list is empty *)
            (if apos = 0 || bpos = 0 then min_int
            else
            let a = Sequence.get a apos
            and b = Sequence.get b bpos in
            optimal -
            (Cost_matrix.Two_D.cost a b cm +
                m1.(apos - 1).(bpos - 1) + 
            m2.(lena - apos).(lenb - bpos))), Empty))
*)

let produce_matrices cm a b =
    let _ = Sequence.Align.cost_2 a b cm Matrix.default in
    let lena = (Sequence.length a) - 1
    and lenb = (Sequence.length b) - 1 in
    let table = Hashtbl.create (lena + 1) in
    let width = min lena lenb
    and height = max lena lenb in
    let set_cost i j cost = 
        if lena >= lenb then set_cost table i j cost
        else set_cost table j i cost
    in
    let add = All_sets.Integers.add in
    let rec do_row row mymap column nextmap =
        let dir = Matrix.get_dir Matrix.default (width + 1) (height + 1) row
        column in
        let nextmap =
            if 0 <> (1 land dir) then begin 
                set_cost row column (0, Empty);
                add (column - 1) nextmap
            end else nextmap 
        in
        let nextmap =
            if 0 <> (4 land dir) then add column nextmap else nextmap
        in
        if 0 <> (2 land dir) then
            if column = 0 || All_sets.Integers.mem (column - 1) mymap then
                nextmap
            else do_row row mymap (column - 1) nextmap 
        else nextmap
    in
    let rec do_matrix row row_columns =
        let next_row =
            All_sets.Integers.fold (do_row row row_columns) 
            row_columns All_sets.Integers.empty 
        in
        if row = 0 then ()
        else do_matrix (row - 1) next_row
    in
    do_matrix height (All_sets.Integers.singleton width);
    table


let get_route_weight = function
    | Empty -> 0
    | Single (_, _, _, x, _, _)
    | Multiple (_, x, _, _) -> x

let get_route_max_weight = function
    | Empty -> 0
    | Single (_, _, _, _, _, x)
    | Multiple (_, _, _, x) -> x

let get_route_paths = function 
    | Empty -> 1
    | Single (_, _, _, _, x, _)
    | Multiple (_, _, x, _) -> x

let invert len table = 
    let res = Hashtbl.create len in
    Hashtbl.iter (fun a tbl ->
        Hashtbl.iter (fun b v ->
            set_cost res b a v) tbl) table;
    res

type table = (int, (int, (int * routes)) Hashtbl.t) Hashtbl.t

let infer_transitive cm a b c (mab : table) (mbc : table) : table = 
    let mba : table = invert (Sequence.length b) mab in
    let res = Hashtbl.create (Sequence.length a) in
    Hashtbl.iter (fun b_pos atbl ->
        if Hashtbl.mem mbc b_pos then
            let ctbl = Hashtbl.find mbc b_pos in
            Hashtbl.iter (fun a_pos (_, left_route) ->
                let to_remove = ref [] in
                let res_a_table = Hashtbl.create 59 in
                Hashtbl.add res a_pos res_a_table;
                Hashtbl.iter (fun c_pos (_, right_route) ->
                    if has_cost res a_pos c_pos then to_remove := c_pos ::
                        !to_remove
                    else
                    let total = 
                        (get_route_weight left_route) +
                        (get_route_weight right_route) + 
                        b_pos
                    in
                    let max_total =
                        (get_route_max_weight left_route) +
                        (get_route_max_weight right_route) + 
                        b_pos
                    in
                    let paths = 
                        get_route_paths left_route *
                        get_route_paths right_route
                    in
                    if paths > 1 then ()
                    else
                    let new_route = 
                        Single (left_route, b_pos, right_route, total,
                        paths, max_total)
                    in
                    Hashtbl.replace res_a_table c_pos (0, new_route)) ctbl;
                    List.iter (Hashtbl.remove res_a_table) !to_remove)
            atbl) mba;
        res
(*

let infer_transitive cm a b c mab mbc =
    let lena = Sequence.length a in
    let res = Hashtbl.create lena in
    Hashtbl.iter (fun a btbl ->
        and a_pos = ref [] in
        Hashtbl.iter (fun pos (v, inter) ->
                a_pos := (pos, inter) :: !a_pos) btbl;
        if !a_score > 0 then ()
        else begin
            List.iter (fun (b_pos, left_route) ->
                try
                let tbl = Hashtbl.find mbc b_pos in
                Hashtbl.iter (fun c (cost, right_route) ->
                    let total = 
                        (get_route_weight left_route) +
                        (get_route_weight right_route) + 
                        b_pos
                    in
                    let max_total =
                        (get_route_max_weight left_route) +
                        (get_route_max_weight right_route) + 
                        b_pos
                    in
                    let paths = 
                        get_route_paths left_route *
                        get_route_paths right_route
                    in
                    let new_route = 
                        Single (left_route, b_pos, right_route, total,
                        paths, max_total)
                    in
                    try 
                        let (_, x) = get_cost res a c in
                        let w = get_route_weight x in
                        let mw = get_route_max_weight x in
                        let np = get_route_paths x in
                        set_cost res 
                        a c
                        (0, Multiple (new_route :: [x], min w total, np
                        + paths, max mw max_total))
                    with
                    | Not_found -> 
                            set_cost res a c (cost, new_route)) 
                tbl with Not_found -> ()) !a_pos
        end) mab;
    res
*)
let check_transitivity cm a b c =
    let mab = produce_matrices cm a b 
    and mbc = produce_matrices cm b c in
    infer_transitive cm a b c mab mbc 


let correlation threshold a b =
    let rows = Array.length a in
    if rows > 0 then begin
        let res = ref 0 in
        let columns = Array.length a.(0) in
        for i = rows - 1 downto 1 do
            for j = columns - 1 downto 1 do
                if threshold < (fst a.(i).(j)) && threshold < (fst b.(i).(j)) then
                    let r = (fst a.(i).(j)) - (fst b.(i).(j)) in
                    res := (r * r) + !res;
            done;
        done;
        Some !res
    end else None

let zeroes trans trans' real = 
    let res = ref [] in
    Hashtbl.iter (fun i tbl ->
        Hashtbl.iter (fun j tbl ->
            try
                let x, routes = get_cost trans i j in
                let y, _ = get_cost trans' i j in
                let cot = fst (get_cost real i j) in
                if x = cot && y = x then
                    res := (i, j, routes) :: !res;
            with
            | Not_found -> ()) tbl) trans;
    List.sort (fun (a, b, _) (c, d, _) -> 
        if a - c = 0 then b - d
        else a - c) !res

let positions cm sa sb sc a b =
    let trans = infer_transitive cm sa sb sc a b in
    let real = produce_matrices cm sa sc in
    zeroes trans real

let rec sum_smallest route =
    match route with
    | Empty -> []
    | Single (left, center, right, total, _, _) ->
            let solutionleft = sum_smallest left 
            and solutionright = sum_smallest right in
            solutionleft @ [center] @ solutionright
    | Multiple (lst, total, _, _) -> 
            let solution = List.find (function
                | Empty -> total = 0
                | Multiple (_, x, _, _)
                | Single (_, _, _, x, _, _) -> x = total) lst
            in
            sum_smallest solution

let rec sum_largest route =
    match route with
    | Empty -> []
    | Single (left, center, right, _, _, total) ->
            let solutionleft = sum_largest left 
            and solutionright = sum_largest right in
            solutionleft @ [center] @ solutionright
    | Multiple (lst, _, _, total) -> 
            let solution = List.find (function
                | Empty -> total = 0
                | Multiple (_, _, _, x)
                | Single (_, _, _, _, _, x) -> x = total) lst
            in
            sum_largest solution


let rec gen_smallest comparison extreme ((best, bestroute) as acc) route =
    match route with
    | Empty -> acc
    | Single (left, center, right, _, _, _) ->
            if comparison center best then 
                center,
                (snd (gen_smallest comparison extreme (extreme, []) left)) 
                @ [center] @ 
                (snd (gen_smallest comparison extreme (extreme, []) right))
            else acc
    | Multiple (lst, _, _, _) -> List.fold_left (gen_smallest comparison extreme) acc lst
    
let smallest route = snd (gen_smallest ( < ) max_int (max_int, []) route)

let largest route = snd (gen_smallest ( > ) min_int (min_int, []) route)

let rec comparator cmp a b =
    match a, b with
    | ha :: ta, hb :: tb -> cmp ha hb && (comparator cmp ta tb)
    | [], [] -> true
    | _ -> false

let select_paths lst =
    let is_better (a_fst, a_lst, a_int) (b_fst, b_lst, b_int) =
        let compare a b = a >= b in
        compare a_fst b_fst && compare a_lst b_lst && 
        comparator compare a_int b_int
    in
    let newpos (a, _, _) (b, _, _) = a <> b in
    let rec aux acc fst lst =
        match lst with
        | [] -> (fst :: acc)
        | h :: t ->
                if newpos fst h then 
                    aux (fst :: acc) h t
                else if is_better fst h then
                    aux acc h t
                else  aux acc fst t
    in
    match lst with
    | h :: t -> aux [] h t
    | [] -> []

let split_paths lst = 
    let is_at_1 (a, b, lst) (x, y, mst) =
        let compare a b = a = b + 1 in
        compare a x && compare b y && comparator compare lst mst
    in
    let rec aux block acc cur lst =
        match lst with
        | [] -> List.rev_map List.rev ((cur :: block) :: acc)
        | h :: t ->
                if is_at_1 h cur then
                    aux (cur :: block) acc h t
                else aux [] ((cur :: block) :: acc) h t
    in
    match lst with
    | h :: t -> aux [] [] h t
    | [] -> []

let do_splits minlen a b lst seqs = 
    let lst = split_paths lst in
    let lst = List.filter (fun x ->  
        let len = List.length x in
        minlen <= len) lst in
    let () = 
        match lst with
        | [x] -> 
                let lena = float_of_int (Sequence.length a)
                and lenb = float_of_int (Sequence.length b)
                and lenx = float_of_int (List.length x) in
                let lena, lenb = (max lena lenb), (min lena lenb) in
                if lena <= (1.1 *. lenb) && (lenx *. 1.1) >= lena then
                    Status.user_message
                    Status.Information
                    ("This sequence can be analyzed prealigned. For this " ^
                    "reason I will not break it, instead, use the " ^
                    "transform (prealigned) command. See the program " ^
                    "documentation.")
        | _ -> ()
    in
    let lst = List.map (fun x -> 
        let (a, b, c) = List.hd x in
        a :: (c @ [b])) lst 
    in
    let invert lst =
        match lst with
        | h :: _ -> 
                let res = List.map (fun _ -> []) h in
                List.fold_right (List.map2 (fun x y -> x :: y)) lst res
        | [] -> []
    in
    let lst = invert lst in
    let rec make_pairs lst =
        match lst with
        | a :: ((b :: _) as tl) ->
                (a, b) :: make_pairs tl 
        | [a] -> [(a, a)]
        | [] -> []
    in
    List.map make_pairs lst

let of_list cm original = 
    let rec transitives (fst, prev, preex) lst =
        match lst with
        | a :: tl ->
                let nm = produce_matrices cm prev a in
                transitives (fst, a, infer_transitive cm fst prev a preex nm) tl
        | [] -> (fst, prev, preex)
    in
    let a, b, original = 
        match original with
        | f :: s :: original ->
                let original = Array.of_list original in
                Array_ops.randomize original;
                let o = Array.to_list original in
                Array_ops.randomize original;
                let o' = Array.to_list original in
                let a = f :: (List.rev (s :: o)) in
                a, f :: (List.rev (s :: o')), a

        | x -> x, x, x
    in
    let a = List.map snd a in
    let b = List.map snd b in
    match a, b with
    | fa :: sa :: ta, fb :: sb :: tb ->
            transitives (fa, sa, produce_matrices cm fa sa) ta, 
            transitives (fb, sb, produce_matrices cm fb sb) tb, 
            original
    | _ -> assert false

let at_most_2_options (a, b, route) = 2 >= get_route_paths route

let partition cm lst =
    let (a, b, trans), (_, _, trans'), original = of_list cm lst in
    let trans = zeroes trans trans' (produce_matrices cm a b)  in
    let splits = 
        let trans = List.filter at_most_2_options trans in
        let len = (List.length trans) / 2 in
        let cnt = ref 0 in
        do_splits 4 a b 
        (select_paths (List.rev_map (fun (a,b,c) -> 
            incr cnt;
            a, b,
            if !cnt > len then sum_smallest c else sum_largest c) trans))
        original
    in
    if 0 = List.length splits then []
    else begin
        assert (List.length original = List.length splits);
        List.map2 (fun (name, seq) splits -> (name, splits))
        original splits 
    end

