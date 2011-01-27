let failwithf format = Printf.ksprintf failwith format
let (-->) b a = a b

type tree = string option * Tree.Parse.tree_types

let read_tree filename = 
    List.map (fun ts -> List.map (fun t -> (Some filename,t)) ts) 
             (Tree.Parse.of_channel (open_in filename))

let read_tree_nexus filename : (string option * Tree.Parse.tree_types) list =
    let trees : Nexus.P.tree list = 
        filename --> open_in
                 --> Lexing.from_channel
                 --> Nexus.Grammar.trees Nexus.Lexer.tree_tokens
    in
    let parsed = 
        List.map (Nexus.File.generate_parser_friendly [] [| None |]) trees
    in
    List.map
        (function
            | Some "",[x]
            | None,[x]   -> (Some filename,x)
            | Some n,[x] -> (Some (filename^":"^n),x)
            | _ -> failwith "nope")
        parsed

let load_trees is_nexus lst = 
    let res = if is_nexus then List.map read_tree_nexus lst
              else List.flatten (List.map read_tree lst) in
    List.map (fun (f,t) -> (f,Tree.Parse.cannonic_order t)) (List.flatten res)
 
let get_name x = match fst x with | Some x -> x | None -> ""

let create_taxon_table trees = 
    let res = Hashtbl.create 31
    and num = ref 0 in
    List.iter 
        (fun tree -> 
            let _ = 
                Tree.Parse.map_tree
                    (function 
                     | "" -> "" 
                     | x  when Hashtbl.mem res x -> x
                     | x  -> Hashtbl.add res x !num; incr num; x)
                    (snd tree)
            in () )
        trees;
    res

let set_of_tbl table= 
    Hashtbl.fold (fun name code set -> All_sets.Integers.add code set)
                 (table)
                 (All_sets.Integers.empty)

let print_set all set = 
    let opp = All_sets.Integers.diff all set in
    Printf.printf "(";
    All_sets.Integers.iter (fun e -> Printf.printf "%d " e) set;
    Printf.printf ") / (";
    All_sets.Integers.iter (fun e -> Printf.printf "%d " e) opp;
    Printf.printf ")"

let print_tbl tbl = 
    Printf.printf "Taxa Association Table:\n";
    Hashtbl.iter
        (fun k v -> Printf.printf "\t %d -- %s\n" v k)
        (tbl);
    print_newline ()

let compare_set all a b =
    let other_b = All_sets.Integers.diff all b in
    let equality_one = All_sets.Integers.compare a b
    and equality_two = All_sets.Integers.compare a other_b in
    (equality_one = 0) || (equality_two = 0)

let unique all =
    let rec unique acc = function
        | [] -> acc
        | hd::tl ->
            if List.exists (compare_set all hd) acc 
                then unique acc tl
                else unique (hd::acc) tl
    in
    unique []

let create_partitions all tbl tree : All_sets.Integers.t list =
    let rec union_lst acc = function
        | []     -> acc
        | hd::tl -> union_lst (All_sets.Integers.union acc hd) tl
    and traversal f = function
        | Tree.Parse.Leafp d ->
            [(All_sets.Integers.singleton (Hashtbl.find tbl (f d)))]
        | Tree.Parse.Nodep (childs,data) ->
            let subtree_sets = List.concat (List.map (traversal f) childs) in
            (union_lst All_sets.Integers.empty subtree_sets) :: (subtree_sets)
    and skip_root_traversal f = function
        | Tree.Parse.Nodep (childs,data) ->
            List.flatten (List.map (traversal f) childs)
        | Tree.Parse.Leafp d ->
            [(All_sets.Integers.singleton (Hashtbl.find tbl (f d)))]
    in
    let res = match snd tree with
        | Tree.Parse.Annotated (t,_) -> skip_root_traversal (fun x -> x)     t
        | Tree.Parse.Flat t          -> skip_root_traversal (fun x -> x)     t
        | Tree.Parse.Branches t      -> skip_root_traversal (fun x -> fst x) t
        | Tree.Parse.Characters t    -> skip_root_traversal (fun x -> fst x) t
    in
    unique all res
 
let compare_partitions verbose taxon_tbl all (t1,p1) (t2,p2) =
    let remove_partition p ps = List.filter (fun x -> not (compare_set all p x)) ps in
    let p_diff21 = List.fold_left (fun acc x -> remove_partition x acc) p2 p1 in
    let p_diff12 = List.fold_left (fun acc x -> remove_partition x acc) p1 p2 in
    if verbose then begin
        Printf.printf "Partition differences in %s and %s\n%!" (get_name t1) (get_name t2);
        List.iter (fun x -> Printf.printf "%s -- " (get_name t1);
                            print_set all x;
                            print_newline ())
                  (p_diff21);
        List.iter (fun x -> Printf.printf "%s -- " (get_name t2);
                            print_set all x;
                            print_newline ())
                  (p_diff12);
        print_newline ()
    end;
    (List.length p_diff21) + (List.length p_diff12)

let print_matrix mat labels =
    for i = 0 to ((Array.length mat)-1) do
        Printf.printf "%s -\t" (match fst (fst labels.(i)) with | Some x -> x | None -> "");
        for j = 0 to ((Array.length mat.(i))-1) do
            Printf.printf "%d%! " mat.(i).(j);
        done;
        print_newline ()
    done;
    ()

let init v ts =
    let tbl = create_taxon_table ts in
    let set = set_of_tbl tbl in
    let num = List.length ts in
    let tps = Array.of_list (List.map (fun t -> (t,create_partitions set tbl t)) ts) in
    let mat = Array.make_matrix num num 0 in
    if v then print_tbl tbl;
    for i = 0 to (num-1) do
        for j = i+1 to (num-1) do
            let dist = compare_partitions v tbl set tps.(i) tps.(j) in
            mat.(i).(j) <- dist;
            mat.(j).(i) <- dist;
        done;
    done;
    print_matrix mat tps

let cmd_line () = match Array.to_list Sys.argv with
    | []     -> failwith  "no arguments?"
    | hd::[] -> failwithf "usage: %s <treefile1> <treefile2> ..." hd
    | _::tl  ->
        let v, tl = 
            if (List.mem "-v" tl) || ( List.mem "--verbose" tl) then
                let tl_minus = 
                    List.filter (fun x -> match x with 
                                    | "-v" | "--verbose" -> false 
                                    | _ -> true) tl
                in
                true,tl_minus
            else 
                false,tl
        in
        begin match load_trees true tl with
            | []             -> failwith "No trees loaded; require at least two."
            | hd::[]         -> failwith "One tree loaded; require at least two."
            | (_::_ as all)  -> init v all
        end

let () = if !Sys.interactive then () else cmd_line ()
