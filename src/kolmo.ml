(* A series of functions that reflect the concepts and function definitions from 
* An Introduction to Kolmogorov Complexity and its Applications by Li and
* Vitanyi *)

module type E = sig
    (** A bit representation *)
    type bit = Zero | One
    (** Encoding a natural number *)
    type encoding = bit list 
    (* The representation of a natural number *)
    type natural

    (** Type Conversion *)
    val to_nat : int -> natural
    val of_nat : natural -> int

    (** The length of the encoding of a natural number *)
    val l : natural -> natural

    (** THe binary encoding of a natural number *)
    val binary : natural -> encoding

    (** Some prefix free encodings *)
    val e : natural -> natural -> encoding
    val e_2 : natural -> encoding
    val hat : natural -> encoding

    (** Decoding the previous prefix free codifications *)
    val decode : encoding -> natural
    (* A function to generate the huffman code functions for lists of elements
    * of type 'a *)

    val huffman : ('a * float) list -> 
        ((encoding -> 'a list) * ('a list -> encoding))

    (* The tree representation of a huffman code *)
    val huffman_tree : ('a * float) list -> 'a option Parser.Tree.t

end 

module Encodings : E = struct
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


end

let ( --> ) a b = b a

module SK = struct

    open Parser.Tree

    type 'a t = Leaf of 'a | Node of ('a t list)
    type primitives = S | K | Label of string
    type expression = primitives t
    type sk = String of string | Processed of expression

    exception Illegal_Expression of string Parser.Tree.t list

    let universe = Hashtbl.create 97 

    let of_string string =
        let no_forest = function
            | [x] -> 
                    let rec reverse tree = 
                        match tree with
                        | Parser.Tree.Node (s, "") ->
                                Node (List.rev_map reverse s)
                        | Parser.Tree.Leaf "S" -> Leaf S
                        | Parser.Tree.Leaf "K" -> Leaf K
                        | Parser.Tree.Leaf x -> 
                                if Hashtbl.mem universe x then
                                    Hashtbl.find universe x 
                                else Leaf (Label x)
                        | x -> raise (Illegal_Expression [x])
                    in
                    reverse x
            | x -> raise (Illegal_Expression x)
        in
        let res = 
            string --> Parser.Tree.of_string --> List.map no_forest --> List.hd
        in
        Processed res

    let to_string tree =
        match tree with
        | String x -> x
        | Processed tree ->
            let tree = 
                let to_string = function
                    | S -> "S"
                    | K -> "K"
                    | (Label x) -> x
                in
                let rec string_tree x =
                    match x with
                    | Node lst ->
                            Parser.Tree.Node 
                            (List.map string_tree lst, "")
                    | Leaf x -> Parser.Tree.Leaf (to_string x)
                in
                string_tree tree
            in
            let str = ref "" in
            let my_printer x = str := !str ^ x in
            AsciiTree.draw_parenthesis false my_printer tree;
            !str

    let sk_define name tree = 
        match tree with
        | Processed x -> Hashtbl.replace universe name x 
        | String x ->  
                match of_string x with
                | Processed x -> Hashtbl.replace universe name x
                | String _ -> assert false

    let rec simplify tree = 
        match tree with
        | Leaf _ -> [tree]
        | Node lst ->
                match lst with
                | [(Leaf S); _]
                | [(Leaf K); _]
                | [(Leaf _)] 
                | [(Leaf S); _; _] -> lst
                | (Leaf K) :: a :: b :: t -> 
                        let a = simplify a in
                        [Node (a @ t)]
                | (Leaf S) :: a :: b :: c :: t ->
                        let a = simplify a in
                        [Node ((a @ [c] @ [Node [b; c]] @ t))]
                | h :: t -> 
                        let h = simplify h in
                        [Node (h @ t)]
                | [] -> raise (Illegal_Expression [])


    let rec reduce tree =
        match tree with
        | String tree -> reduce (of_string tree)
        | Processed tree ->
            match simplify tree with
            | [ntree] -> 
                    if ntree = tree then Processed tree 
                    else reduce (Processed ntree)
            | [] -> raise (Illegal_Expression [])
            | lst -> 
                    let ntree = Node lst in
                    if ntree = tree then Processed tree 
                    else reduce (Processed ntree)

    let evaluate x = 
        x --> of_string --> reduce 

    let test lst = 
        reduce 
        (Processed (Node (List.map (function Processed x -> x | String x ->
            (match of_string x with
            | Processed x -> x
            | String _ -> assert false)) lst)))

    let s_encode tree = 
        (* We define a function to encode an SK expression in a list of bits *)
        let rec to_bit_list acc tree =
            match tree with
            | Node lst ->
                    (match lst with
                    | [] -> assert false
                    | [x] -> to_bit_list acc x
                    | [x; y] ->
                            let acc = to_bit_list acc y in
                            let acc = to_bit_list acc x in
                            Encodings.One :: acc
                    | lst ->
                            let rec to_binary_list lst = 
                                match lst with
                                | [h; t] -> Node lst
                                | f :: s :: t -> 
                                        to_binary_list 
                                        ((Node [f; s]) :: t)
                                | _ -> assert false
                            in
                            to_bit_list acc (to_binary_list lst))
            | Leaf x ->
                    match x with
                    | S -> Encodings.Zero :: Encodings.Zero :: acc
                    | K -> Encodings.Zero :: Encodings.One :: acc
                    | Label _ -> failwith "I can't encode labels"
        in
        match tree with
        | Processed tree -> to_bit_list [] tree
        | String tree ->
                match tree --> of_string --> reduce with
                | Processed tree -> to_bit_list [] tree
                | String _ -> assert false

        let rec s_decode lst = 
            match lst with
            | Encodings.One :: t ->
                    let a, t = s_decode t in
                    let b, t = s_decode t in
                    "( " ^ a ^ b ^ ") ", t
            | Encodings.Zero :: Encodings.Zero :: t -> "S ", t
            | Encodings.Zero :: Encodings.One :: t -> "K ", t
            | _ -> failwith "Illegal encoding"

        let create pattern label =
            (* A function to create combinator that produce the desired pattern
            * with the label as argument, in other words, extracts the label
            * that we want to have inside the pattern as an argument *)
            let rec contains_label tree =
                match tree with
                | Leaf (Label x) -> x = label
                | Leaf K
                | Leaf S -> false
                | Node lst -> List.exists contains_label lst
            in
            let k = Leaf K 
            and s = Leaf S 
            and split lst = 
                match List.rev lst with
                | h :: t -> Node (List.rev t), h 
                | [] -> assert false
            in
            let rec extract tree =
                if not (contains_label tree) then
                    Node [k; tree]
                else 
                    match tree with
                    | Leaf (Label x) -> 
                            assert (x = label);
                            Hashtbl.find universe "I"
                    | Leaf x ->
                            Node [k; tree]
                    | Node [x] ->
                            extract x
                    | Node [x; Leaf (Label l)] when 
                        (l = label) && (not (contains_label x)) ->
                            x
                    | Node [x; y] ->
                            Node [s; extract x; extract y]
                    | Node lst ->
                            let a, b = split lst in
                            extract (Node [a; b])
            in
            match pattern with
            | Processed tree -> Processed (extract tree)
            | String tree -> 
                    match of_string tree with
                    | Processed tree -> Processed (extract tree)
                    | String _ -> assert false

        let create pattern labels = 
            List.fold_left create pattern (List.rev labels)

        let () =
            (* The standard list of items that are defined *)
            let predefined = 
                [
                    ("True", "( K )");
                    ("False", "( S K )");
                    ("1", "(K)");
                    ("0", "(S K)");
                    ("I", "(S K K)");
                    ("Lambda", "(K S)");
                    ("Pair", 
                    "(S (S (K S)(S (K K)(S (K S)(S (S (K S)(S K)) K))))(K K))");
                    ("EmptyList", "(Pair Lambda Lambda)");
                    ("Not", "(Pair (S K) K)");
                    ("Hd", "(S (S K K) ((S (S K) ( K K) )))");
                    ("Tl", "(S (S K K) ((S (S K) ( K (S K)))))");
                    ("Predecessor", "(Tl)"); 
                    ("Successor", "(Pair K)"); 
                    ("Apply", "(S (S K) (S K K))");
                    ("NotZero", "(S (S Hd (K 1)) (K 1))"); 
                    ("Fixedpoint", "(S S K (S (K (S S(S(S S K))))K))");
                ]
            in
            List.iter (fun (a, b) -> sk_define a (String b)) predefined;
            let generated = [
                ("Right", "(a (b c))", ["a"; "b"; "c"]);
                ("GenerateRecursive", "(test max (R (next max) (update \
                acc)) acc)", ["test"; "next"; "update"; "R"; "max"; "acc"]);
                ("GenerateRecursive2", "(test max (R (next max) (update \
                max acc)) acc)", ["test"; "next"; "update"; "R"; "max"; "acc"]);
                (* Now we can define the addition, substraction, and the
                * multiplication functions. *)
                ("Add", "(Fixedpoint (GenerateRecursive NotZero 
                Predecessor Successor))", []);
                ("Substract", "(Fixedpoint (GenerateRecursive NotZero 
                Predecessor Predecessor) a b)", ["b"; "a"]);
                ("Multiply", "(Fixedpoint (GenerateRecursive 
                NotZero Predecessor (Add a)) b EmptyList)", ["a"; "b"]);
                ("ProcessBinary", "((Hd a) (Successor EmptyList) EmptyList)", 
                ["a"]);
                ("AppendBinary", "(Add (ProcessBinary a) (Multiply (Successor
                (Successor EmptyList)) b))", ["a"; "b"]);
                ("Decode", "(Fixedpoint (GenerateRecursive2 NotZero
                Predecessor AppendBinary) a EmptyList)", ["a"]);
                ("GenerateInsertion", "(NotZero pos (Pair (Hd seq) (R
                base (Predecessor pos) (Tl seq))) (Pair base seq))", 
                ["R"; "base"; "pos"; "seq"]);
                ("Insert", "(Fixedpoint GenerateInsertion)", []);
                ("GenerateDeletion", "(NotZero pos (Pair (Hd seq) (R base
                (Predecessor pos) (Tl seq))) seq)", ["R"; "pos"; "seq"]);
                ("Delete", "(Fixedpoint GenerateDeletion)", []);
                ("GenerateSubstitution", "(NotZero pos (Pair (Hd seq) (R base
                (Predecessor pos) (Tl seq))) (Pair base (Tl seq)))", ["R";
                "pos"; "base"; "seq"]);
                ("Substitution", "(Fixedpoint GenerateSubstitution)", []);
                ]
            in
            List.iter (fun (a, b, c) -> 
                let b = create (String b) c in
                sk_define a b) generated

end

module SK_f = struct
    (* We define a module for the simplest of all lambda calculus models, the S
    * - K model of computation. We implement it to try new combinators for this
    * lambda calculus *)

    (* We happily define the basic primitives, this is so beautiful! *)
    let s a b c = a c (b c)
    let k a b = a

    let falso a b = k a b
    let verda a b = s k a b

end
