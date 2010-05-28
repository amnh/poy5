let failwithf format = Printf.ksprintf failwith format
let (-->) b a = a b

let print_set set = 
    Printf.printf "(";
    All_sets.Integers.iter (fun e -> Printf.printf "%d " e) set;
    Printf.printf ")"

let read_tree filename = 
    List.map (fun ts -> List.map (fun t -> (filename,t)) ts) 
             (Tree.Parse.of_channel (open_in filename))

let load_trees lst = 
    lst --> List.map read_tree
        --> List.flatten --> List.flatten
        --> List.map (fun (f,t) -> (f,Tree.Parse.cannonic_order t))
 
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
 
(* they must be equal, or intersection = null *)
let compare_partitions all a b =
    let other_b = All_sets.Integers.inter all b in
    let equality_one = All_sets.Integers.compare a b
    and equality_two = All_sets.Integers.compare a other_b in
    (equality_one = 0) || (equality_two = 0)

let compare_branches table set_of_all partitions_one partitions_two =
    let modified = ref 0 in
    List.iter
        (fun (set,b) ->
            let (_,y) =
                List.find (fun (x,_) -> compare_partitions set_of_all set x)
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

let compare_trees verbose taxon_tbl set_of_all t1 t2 = 
    let partitions_one = create_partitions taxon_tbl t1
    and partitions_two = create_partitions taxon_tbl t2 in
    let s,a = 
        List.fold_left
            (fun (same,all) (elm,_) ->
                try let _ = 
                        List.find (fun (x,_) -> compare_partitions set_of_all elm x)
                                  partitions_two
                    in 
                    (same+1,all+1)
                with | Not_found ->
                    if verbose then begin
                        Printf.printf "\t"; print_set elm; print_newline ()
                    end;
                    (same,all+1))
            (0, 0)
            partitions_one
    in

    if verbose then begin
        if s = a then begin
            Printf.printf "\tSame topology\n%!";
            compare_branches taxon_tbl set_of_all partitions_one partitions_two
        end else
            Printf.printf "\tOf %d partitions, %d are the same\n%!" a s
    end;
    a - s


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

let matrix_destructive_map upper matrix func = 
    for i = 0 to (Array.length matrix) - 1 do
        let j_start = if upper then i else 0 in
        for j = j_start to (Array.length matrix.(i) - 1) do
            let res = func i j in
            if upper then begin
                matrix.(i).(j) <- res; matrix.(j).(i) <- res
            end else begin
                matrix.(i).(j) <- res
            end
        done;
    done;
    ()

let print_matrix row_leads matrix = 
    for i = 0 to (Array.length matrix) - 1 do
        Printf.printf "%s--\t" (fst row_leads.(i));
        for j = 0 to (Array.length matrix.(i) -1) do
            Printf.printf "%d\t" matrix.(i).(j)
        done;
        print_newline ()
    done

let compare_matrix verbose all = 
    let table = create_taxon_table (snd (List.hd all)) in
    let set_of_all = 
        Hashtbl.fold (fun name code set -> All_sets.Integers.add code set)
                     table (All_sets.Integers.empty)
    in
    if verbose then begin
        Printf.printf "Table of taxon codes:\n";
        print_table table
    end;
    let tree_assoc = Array.of_list all in
    let length = Array.length tree_assoc in
    let costs_matrix = Array.make_matrix length length 0 in
    let () = 
        matrix_destructive_map true costs_matrix
            (fun i j -> 
                compare_trees verbose table set_of_all
                    (snd tree_assoc.(i)) (snd tree_assoc.(j)))
    in
    print_matrix tree_assoc costs_matrix;
    ()

let () = match Array.to_list Sys.argv with
    | []     -> failwith  "no arguments?"
    | hd::[] -> failwithf "usage: %s <treefile1> <treefile2> ..." hd
    | _::tl ->
        let verbose,tl = 
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
        begin match load_trees tl with
            | []        -> failwith "No trees loaded; require at least two"
            | hd::[]    -> failwith "One tree loaded; require at least two"
            | (hd::rest as all) -> compare_matrix verbose all
        end


