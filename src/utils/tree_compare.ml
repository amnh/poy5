let failwithf format = Printf.ksprintf failwith format
let (-->) b a = a b

let print_set set = 
    Printf.printf "(";
    All_sets.Integers.iter (fun e -> Printf.printf "%d " e) set;
    Printf.printf ")"

let read_tree filename = Tree.Parse.of_channel (open_in filename)

let load_trees lst = 
    lst --> List.map read_tree
        --> List.flatten --> List.flatten
        --> List.map Tree.Parse.cannonic_order
 
let create_partitions tbl tree =
    let rec union_lst acc = function
        | []     -> acc
        | (hd,_)::tl -> union_lst (All_sets.Integers.union acc hd) tl
    and traversal f b = function
        | Tree.Parse.Leafp d -> 
            [(All_sets.Integers.singleton (Hashtbl.find tbl (f d)),b d)]
        | Tree.Parse.Nodep (childs,data) ->
            let subtree_sets = List.concat (List.map (traversal f b) childs) in
            (union_lst All_sets.Integers.empty subtree_sets,b data) :: (subtree_sets)
    in
    match tree with
    | Tree.Parse.Annotated (t,a) ->
            traversal (fun x -> x) (fun _ -> None) t
    | Tree.Parse.Branches t ->
            traversal (fun x -> fst x) (fun (_,b) -> b) t
    | Tree.Parse.Flat t ->
            traversal (fun x -> x) (fun _ -> None) t
    | Tree.Parse.Characters t ->
            traversal (fun x -> fst x) (fun _ -> None) t
 
let compare_branches table partitions_one partitions_two =
    let modified = ref 0 in
    List.iter
        (fun (set,b) ->
            let (_,y) =
                List.find (fun (x,_) -> 0 = (All_sets.Integers.compare set x))
                          partitions_two
            in
            let b_len = match b with | Some x -> x | None -> 0.0 
            and y_len = match y with | Some x -> x | None -> 0.0 in
            if (abs_float (b_len -. y_len)) >= 1.0e-2 then begin
                    Printf.printf "\t%s\t%s  (%f)\t" 
                        (match b with | Some x -> string_of_float x | None -> "none\t")
                        (match y with | Some x -> string_of_float x | None -> "none\t") 
                        (abs_float (b_len -. y_len));
                    print_set set;
                    print_newline ();
                    incr modified
            end)
        partitions_one;
    if !modified = 0 then Printf.printf "\tSame branch lengths\n%!"

let compare_trees taxon_tbl t1 t2 = 
    let partitions_one = create_partitions taxon_tbl t1
    and partitions_two = create_partitions taxon_tbl t2 in
    let s,a = 
        List.fold_left
            (fun (same,all) (elm,_) ->
                try let _ = 
                        List.find (fun (x,_) -> 0 = (All_sets.Integers.compare elm x))
                                  partitions_two
                    in 
                    (same+1,all+1)
                with | Not_found ->
                    Printf.printf "\t"; print_set elm; print_newline ();
                    (same,all+1))
            (0, 0)
            partitions_one
    in

    if s = a then begin
        Printf.printf "\tSame topology\n%!";
        compare_branches taxon_tbl partitions_one partitions_two
    end else
        Printf.printf "\tOf %d partitions, %d are the same\n%!" a s

let create_taxon_table tree = 
    let res = Hashtbl.create 31
    and num = ref 0 in
    let _ = 
        Tree.Parse.map_tree
            (function 
             | "" -> "" 
             | x  -> Hashtbl.add res x !num; incr num; x)
            tree
    in
    res

let print_table tbl = 
    Hashtbl.iter (fun a b -> Printf.printf "\t%d -- %s\n%!" b a) tbl

let () = match Array.to_list Sys.argv with
    | []     -> failwith  "no arguments?"
    | hd::[] -> failwithf "usage: %s <treefile1> <treefile2> ..." hd
    | _::tl ->
        begin match load_trees tl with
            | []        -> failwith "No trees loaded; require at least two"
            | hd::[]    -> failwith "One tree loaded; require at least two"
            | (hd::rest as all) ->
                let table = create_taxon_table hd in
                Printf.printf "Table of taxon codes:\n";
                print_table table;
                Printf.printf "Comparing first tree with %d trees\n" (List.length rest);
                List.iter (compare_trees table hd) rest
        end


