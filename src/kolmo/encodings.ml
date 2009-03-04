type bit = Zero | One
type encoding = bit list
type natural = int

let to_nat (x : int)  =
    if x < 0 then failwith "Illegal argument"
    else (x : natural)

let of_nat x = x

(* From page 13 l(x) *)
let l (x : natural) = 
    1 + (truncate ((log (float_of_int x)) /. (log 2.)))

let to_bit = function
    | 0 -> Zero
    | 1 -> One
    | _ -> assert false

(* The binary representation of an integer *)
let binary x = 
    if x = 0 then []
    else 
        let rec aux_binary x = 
            if x = 0 then []
            else (to_bit  (x land 1)) :: (aux_binary (x lsr 1)) 
        in
        List.rev (aux_binary x)

(* From page 13, E_i(x) *)
let rec e i x =
    match i with
    | 0 -> 
            let rec make_list len acc = 
                if len = 0 then acc
                else make_list (len - 1) (One :: acc)
            in
            make_list x [Zero] 
    | i -> (e (i - 1) (l x)) @ (binary x)

let e_2 = e 2

let hat = e 1

let nat_of_binary x =
    let aux_nat_of_binary acc = function
        | One -> ((acc lsl 1) lor 1) 
        | Zero -> (acc lsl 1) 
    in
    List.fold_left aux_nat_of_binary 0 x

let decode_list_of_ones lst = 
    let rec aux x = function 
        | Zero :: tl -> x, tl
        | One :: lst -> aux (x + 1) lst
        | [] -> failwith "Illegal encoding"
    in
    aux 0 lst

let extract_first x lst = 
    let rec aux h cnt lst =
        if cnt = 0 then List.rev h, lst
        else 
            match lst with
            | a :: t -> aux (a :: h) (cnt - 1) t
            | [] -> failwith "Illegal encoding"
    in
    aux [] x lst

let decode bin =
    let rec length_decode acc left =
        match left with
        | [] -> acc
        | left ->
                let number, rest_of_list = extract_first acc left in
                let acc = (nat_of_binary number) in
                length_decode acc rest_of_list
    in
    let x, rest_of_list = decode_list_of_ones bin in
    length_decode x rest_of_list

let huffman_tree lst = 
    let rec sort_n_merge lst =
        match List.sort (fun (_, (a : float)) (_, (b : float)) -> 
            compare a b) lst 
        with
        | [(x, _)] ->  x
        | (x, a) :: (y, b) :: tl ->
                let l = 
                    (((Parser.Tree.Node ([x; y], None)), a +. b) :: tl)
                in
                sort_n_merge l
        | [] -> failwith "Empty alphabet"
    in
    let r = List.map (fun (x, y) -> (Parser.Tree.Leaf (Some x)), y) lst in
    sort_n_merge r

let huffman lst =
    let tree = huffman_tree lst in
    let encoded_table = 
        let hshtbl = Hashtbl.create 97 in
        let rec generate_table acc tree =
            match tree with
            | Parser.Tree.Leaf (Some x) ->
                    Hashtbl.add hshtbl x (List.rev acc)
            | Parser.Tree.Node ([a; b], None) ->
                    generate_table (Zero :: acc) a;
                    generate_table (One :: acc) b;
            | Parser.Tree.Node (_, _)
            | Parser.Tree.Leaf None -> assert false
        in
        generate_table [] tree;
        hshtbl
    in
    (fun list_to_decode ->
        let rec aux_decoder (decoded, tree_left) item =
            match tree_left with
            | Parser.Tree.Leaf (Some x) -> 
                    aux_decoder ((x :: decoded), tree) item
            | Parser.Tree.Node ([a; b], None) ->
                    (match item with
                    | Zero -> (decoded, a)
                    | One -> (decoded, b))
            | _ -> failwith "Illegal message"
        in
        match List.fold_left aux_decoder ([], tree) list_to_decode with
        | lst, Parser.Tree.Leaf (Some x) -> List.rev (x :: lst)
        | _ -> failwith "Illegal message"),
    (fun list_to_encode ->
        let res = List.map (Hashtbl.find encoded_table) list_to_encode in
        List.flatten res)

let geometric min max p =
    let q = 1. -. p in
    let rec generator cur mlt acc =
        if cur > max then acc
        else 
            let nmlt = q *. mlt in
            if cur < min then generator (cur + 1) nmlt acc
            else generator (cur + 1) nmlt ((cur - 1, mlt) :: acc)
    in
    List.rev (generator 1 p [])
