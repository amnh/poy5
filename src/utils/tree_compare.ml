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

let create_names tbl1 tbl2 =
    let i = ref 0
    and incr x = let res = !x in incr x; res in
    (fun name ->
        if Hashtbl.mem tbl1 name then
            Hashtbl.find tbl1 name
        else begin
            let res = incr i in
            let () = Hashtbl.add tbl1 name res
            and () = Hashtbl.add tbl2 res name in
            res
        end)

and create_tree_names () =
    let names = Hashtbl.create 16 in
    (fun str_opt -> match str_opt with
        | Some x ->
            if Hashtbl.mem names x
                then Hashtbl.replace names x ((Hashtbl.find names x)+1)
                else Hashtbl.add names x 0;
            Some (x^(string_of_int (Hashtbl.find names x)))
        | None -> assert false)

let load_trees is_nexus lst : Tree.u_tree list =
    let names,codes = Hashtbl.create 16, Hashtbl.create 16 in
    let naming_func = create_names names codes
    and naming_tree = create_tree_names () in
    let res = if is_nexus then List.map read_tree_nexus lst
              else List.flatten (List.map read_tree lst) in
    res --> List.flatten
        --> List.map (fun (f,t) -> (f,Tree.Parse.cannonic_order t))
        --> List.map (fun (n,t) ->
                Tree.Parse.convert_to ((naming_tree n),[t]) naming_func)
 
let get_name x = match x.Tree.tree_name with | Some x -> x | None -> ""

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

let print_matrix mat t =
    for i = 0 to ((Array.length mat)-1) do
        Printf.printf
            "%s -\t" (get_name t.(i));
        for j = 0 to ((Array.length mat.(i))-1) do
            Printf.printf "%d%! " mat.(i).(j);
        done;
        print_newline ()
    done;
    ()

let init v ts =
    let num = List.length ts in
    let mat = Array.make_matrix num num 0 in
    let tps = Array.of_list ts in
    for i = 0 to (num-1) do
        for j = i+1 to (num-1) do
            let dist = Tree.robinson_foulds tps.(i) tps.(j) in
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
            | []            -> failwith "No trees loaded; require at least two."
            | hd::[]        -> failwith "One tree loaded; require at least two."
            | (_::_ as all) -> init v all
        end

let () = if !Sys.interactive then () else cmd_line ()
