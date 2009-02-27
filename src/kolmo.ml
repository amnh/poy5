(* A series of functions that reflect the concepts and function definitions from 
* An Introduction to Kolmogorov Complexity and its Applications by Li and
* Vitanyi *)

let debug_sk = ref false 


module type E = sig
    (** A bit representation *)
    type bit = Zero | One

    (** Encoding a natural number *)
    type encoding = bit list 

    (** The representation of a natural number *)
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

    val geometric : int -> int -> float -> (int * float) list
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

end

let ( --> ) a b = b a

module S_K = struct
    type primitives = [ `S | `K | `Label of string | `Node of primitives list |
        `Debugger of string | `Lazy of primitives Lazy.t]

    exception Illegal_Expression of string Parser.Tree.t list

    let universe = Hashtbl.create 97 

    let of_string string =
        let no_forest = function
            | [x] -> 
                    let rec reverse tree = 
                        match tree with
                        | Parser.Tree.Node (s, "") ->
                                `Node (List.rev_map reverse s)
                        | Parser.Tree.Leaf "S" -> `S
                        | Parser.Tree.Leaf "K" -> `K
                        | Parser.Tree.Leaf x -> `Label x
                        | x -> raise (Illegal_Expression [x])
                    in
                    reverse x
            | x -> raise (Illegal_Expression x)
        in
        string --> Parser.Tree.of_string --> List.map no_forest --> List.hd


    let rec expand ?(except=[]) tree =
        match tree with
        | `Node lst -> 
                `Node (List.map (expand ~except) lst)
        | `Label x ->
                if (List.exists (fun y -> x = y) except) ||
                    (not (Hashtbl.mem universe x)) then tree
                else expand (Hashtbl.find universe x)
        | x -> x

    let to_string tree =
        let buf = Buffer.create 1000 in
        let rec string_tree x =
            match x with
            | `Node lst ->
                    Buffer.add_char buf '(';
                    List.iter string_tree lst;
                    Buffer.add_char buf ')';
            | `Lazy x -> string_tree (Lazy.force_val x)
            | `S | `K | `Label _ | `Debugger _ as x ->
                    let string = 
                        match x with
                        | `S -> "S"
                        | `K -> "K"
                        | `Debugger x -> "Dbg: " ^ x
                        | `Label x -> x
                    in
                    Buffer.add_string buf string;
        in
        string_tree tree;
        Buffer.contents buf

    let rec sk_define name tree = 
        Hashtbl.replace universe name (expand tree) 

    let debug = ref false

    let rec simplify tree = 
        match tree with
        | `Node lst ->
                (match lst with
                | [`S; _]
                | [`K; _]
                | [ `S ] | [ `K ]
                | [(`S); _; _] -> lst
                | (`K) :: a :: b :: t -> 
                        let a = simplify a in
                        simplify (`Node (a @ t))
                | (`S) :: a :: b :: c :: t ->
                        let a = simplify a in
                        let c = 
                            `Lazy (Lazy.lazy_from_fun 
                            (fun () -> `Node (simplify c))) 
                        in
                        simplify (`Node (a @ [c] @ [`Node [b; c]] @ t))
                | `Label _ :: t -> lst
                | `Debugger str :: t ->
                        if !debug then print_endline str; 
                        simplify (`Node t)
                | h :: t ->
                        let h = simplify h in
                        simplify (`Node (h @ t))
                | [] -> raise (Illegal_Expression []))
        | `S | `K | `Label _ -> [tree]
        | `Lazy x -> [Lazy.force_val x]
        | `Debugger str -> 
                if !debug then print_endline str; 
                []


    let rec reduce tree = 
        match simplify tree with
        | [x] -> x
        | x -> `Node x
    (*
        if !debug_sk then print_endline (to_string tree);
        match simplify tree with
        | [ntree], true -> reduce ntree
        | [ntree], false -> ntree
        | [], _ -> raise (Illegal_Expression [])
        | lst, false -> `Node lst
        | lst, true -> 
                match lst with
                | [`S; _]
                | [`K; _]
                | [`S; _; _] -> `Node lst
                | _ -> reduce (`Node lst)
    *)


    let evaluate x = 
        x --> of_string --> reduce 

    let eval x = x --> expand --> reduce 

    let test lst = 
        reduce (`Node lst)

    let s_encode tree = 
        (* We define a function to encode an SK expression in a list of bits *)
        let rec to_bit_list ?(no_print=false) (tree : primitives) =
            match tree with
            | `Debugger _ -> []
            | `Lazy x -> to_bit_list (Lazy.force_val x)
            | `S -> Encodings.Zero :: Encodings.Zero :: []
            | `K -> Encodings.Zero :: Encodings.One :: []
            | `Label _ -> failwith "I can't encode labels"
            | `Node lst ->
                    (match lst with
                    | [] -> assert false
                    | [x] -> to_bit_list x
                    | [x; y] ->
                            let y = to_bit_list y in
                            let x = to_bit_list x in
                            Encodings.One :: (x @ y)
                    | lst ->
                            let rec to_binary_list lst = 
                                match lst with
                                | [h; t] -> `Node lst
                                | f :: s :: t -> 
                                        to_binary_list 
                                        ((`Node [f; s]) :: t)
                                | _ -> assert false
                            in
                            to_bit_list ~no_print:true (to_binary_list lst))
        in
        to_bit_list tree


        let rec s_decode lst = 
            let rec aux_s_decode lst = 
                match lst with
                | Encodings.One :: t ->
                        let a, t = aux_s_decode t in
                        let b, t = aux_s_decode t in
                        `Node [a; b], t
                | Encodings.Zero :: Encodings.Zero :: t -> `S, t
                | Encodings.Zero :: Encodings.One :: t -> `K, t
                | _ -> failwith "Illegal encoding"
            in
            let a, t = aux_s_decode lst in
            match t with
            | [] -> a
            | t -> (`Node [a; s_decode t])

        let create pattern label =
            (* A function to create combinator that produce the desired pattern
            * with the label as argument, in other words, extracts the label
            * that we want to have inside the pattern as an argument *)
            let rec contains_label tree =
                match tree with
                | `Label x -> x = label
                | `Lazy x -> contains_label (Lazy.force_val x)
                | `Debugger _ 
                | `K
                | `S -> false
                | `Node lst -> List.exists contains_label lst
            in
            let k = `K 
            and s = `S 
            and split lst = 
                match List.rev lst with
                | h :: t -> `Node (List.rev t), h 
                | [] -> assert false
            in
            let rec extract tree =
                if not (contains_label tree) then
                    `Node [k; tree]
                else 
                    match tree with
                    | `Label x -> 
                            assert (x = label);
                            Hashtbl.find universe "I"
                    | `S | `K ->
                            `Node [k; tree]
                    | `Debugger _ -> tree
                    | `Lazy x -> extract (Lazy.force_val x)
                    | `Node [x] ->
                            extract x
                    | `Node [x; `Label l] when 
                        (l = label) && (not (contains_label x)) ->
                            x
                    | `Node [x; y] ->
                            `Node [s; extract x; extract y]
                    | `Node lst ->
                            let a, b = split lst in
                            extract (`Node [a; b])
            in
            extract pattern

        let sk_define_interpreted name args tree =
            let tree = expand ~except:args tree in
            let tree = List.fold_left create tree (List.rev args) in
            Hashtbl.replace universe name tree

        let create pattern labels = 
            List.fold_left create pattern (List.rev labels)

        let def a b c = sk_define a (create c b)

        let () =
            (* Many convenient function definitions uisng kolmoExtensions to
            * parse them. *)
            sk_define "I" (SK (S K K));
            def "True" [] (SK K);
            def "False" [] (SK (S K));
            def "O" [] (SK (K));
            def "Z" [] (SK (S K));
            def "I" [] (SK (S K K));
            def "Lambda" [] (SK (K S));
            def "Pair" [] 
                (SK (S (S (K S)(S (K K)(S (K S)(S (S (K S)(S K)) K))))(K K)));
            def "EmptyList" [] (SK (Pair Lambda Lambda));
            def "Not" [] (SK (Pair (S K) K));
            def "Hd" [] (SK (S (S K K) ((S (S K) ( K K) ))));
            def "Tl" [] (SK (S (S K K) ((S (S K) ( K (S K))))));
            def "Predecessor" [] (SK (Tl)); 
            def "Successor" [] (SK (Pair K)); 
            def "Apply" [] (SK (S (S K) (S K K)));
            def "NotZero" [] (SK (S (S Hd (K O)) (K O))); 
            def "Fixedpoint" [] (SK (S S K (S (K (S S(S(S S K))))K)));
            def "Right" (LS a b c) (SK (a (b c)));
            def "GenerateRecursive" 
                (LS test next update R max acc) 
                (SK (test max (R (next max) (update acc)) acc));
            def "GenerateRecursive2" 
                (LS test next update R max acc)
                (SK (test max (R (next max) (update max acc)) acc));
            def "Add" [] 
                (SK (Fixedpoint (GenerateRecursive NotZero Predecessor 
                    Successor)));
            def "Substract" (LS b a) 
                (SK (Fixedpoint (GenerateRecursive NotZero Predecessor 
                    Predecessor) a b));
            def "Multiply" (LS a b) 
                (SK (Fixedpoint (GenerateRecursive 
                    NotZero Predecessor (Add a)) b EmptyList));
            def "ProcessBinary" (LS a) 
                (SK ((Hd a) (Successor EmptyList) EmptyList));
            def "AppendBinary" (LS a b) 
                (SK (Add (ProcessBinary a) (Multiply (Successor
                    (Successor EmptyList)) b)));
            def "DecodeInt" (LS a) 
                (SK (Fixedpoint (GenerateRecursive2 NotZero Predecessor 
                    AppendBinary) a EmptyList));
            def "GeneratePrepend" (LS R pref suf) 
                (SK (NotZero pref (Pair (Hd pref) (R (Tl pref) suf)) suf));
            def "Prepend" (LS R pos ins seq) 
                (SK (Fixedpoint GeneratePrepend));
            def "GenerateAffineInsert" (LS R pos ins seq) 
                (SK (NotZero pos (Pair (Hd seq) (R pos ins
                    (Tl seq))) (Prepend ins seq)));
            def "AffineInsert" (LS R len seq) 
                (SK (Fixedpoint GenerateAffineInsert));
            def "GenerateCut" (LS R len seq) 
                (SK (NotZero len (R (Predecessor len) (Tl seq)) seq));
            def "Cut" [] (SK (Fixedpoint GenerateCut));
            def "GenerateAffineDelete" (LS R pos len seq) 
                (SK (NotZero pos (Pair (Hd seq) (R (Predecessor pos) len
                    (Tl seq))) (Cut len seq)));
            def "AffineDelete" [] 
                (SK (Fixedpoint GenerateAffineDelete));
            def "GenerateCountK" (LS R acc current) 
                (SK (current (R (Successor acc)) acc));
            def "CountK" [] (SK (Fixedpoint GenerateCountK EmptyList))
end

module PM = struct 
    (** This is potentially garbage, but we leave it because we might need some
    * of the effort placed here somewhere else *)
    (* A phylogenetic analysis machine written for the S_K computer 
    * This machine uses the most simply defined functions, so we will start by
    * defining some basic functionality *)
    let u_Tl = SK (S K)
    let u_Hd = SK K

    exception Unsupported

    module type List = sig
        type ocaml_repr
        type ocaml_list
        val empty_list : S_K.primitives
        val is_not_empty : S_K.primitives -> S_K.primitives
        val ml_zero_signal : S_K.primitives -> bool
        val head : S_K.primitives -> S_K.primitives
        val tail : S_K.primitives -> S_K.primitives
        val prepend : S_K.primitives -> S_K.primitives -> S_K.primitives
        val alphabet : (S_K.primitives * ocaml_repr) list
        val to_ocaml : S_K.primitives -> ocaml_repr
        val of_ocaml : ocaml_repr -> S_K.primitives
        val to_list : S_K.primitives -> ocaml_list
        val of_list : ocaml_list -> S_K.primitives
        val complement : S_K.primitives -> S_K.primitives
    end

    module type Tree = sig
        type ocaml_tree
        val leaf : S_K.primitives -> S_K.primitives
        val is_leaf : S_K.primitives -> S_K.primitives
        val ml_leaf_signal : S_K.primitives -> bool
        val contents : S_K.primitives
        val left : S_K.primitives
        val right : S_K.primitives
        val join : S_K.primitives -> S_K.primitives -> S_K.primitives -> S_K.primitives
        val to_tree : S_K.primitives -> ocaml_tree
        val of_tree : ocaml_tree -> S_K.primitives
    end

    module K_S : List with type ocaml_repr = [ `K | `SK] 
    with type ocaml_list = [`K | `SK] list = struct
            type ocaml_repr = [ `K | `SK ]
            type ocaml_list = ocaml_repr list
            let empty_list = SK (Pair (K S) K)
            let is_not_empty x = SK ([x] [u_Hd] K K)
            let ml_zero_signal x = (SK (K S)) = S_K.eval (SK ([x] [u_Hd]))
            let head x = SK ([x] [u_Hd])
            let tail x = SK ([x] [u_Tl])
            let prepend a b = SK (Pair [a] [b])
            let alphabet = [`K, `K; (SK (S K)), `SK]
            let inv_alph = List.map (fun (a, b) -> (b, a)) alphabet
            let to_ocaml a = 
                let a = S_K.eval a in
                List.assoc a alphabet
            let of_ocaml a = List.assoc a inv_alph
            let to_list lst = 
                let rec to_list lst =
                    if not (ml_zero_signal lst) then 
                        let next = to_ocaml (S_K.eval (head lst)) in
                        next :: (to_list (S_K.eval (tail lst)))
                    else []
                in
                to_list (S_K.eval lst) 
            let rec of_list lst = 
                List.fold_right (fun x acc ->
                    prepend (of_ocaml x) acc) lst empty_list
            let complement x = SK ([x] (S K) K)
    end

    type 'a btree = Node of ('a * 'a btree * 'a btree) | Leaf of 'a

    module Tree : Tree with 
        type ocaml_tree = S_K.primitives btree = struct

        type ocaml_tree = S_K.primitives btree

        let leaf x = SK (Pair K [x])

        let is_leaf x = SK ([x] [u_Hd])

        let ml_leaf_signal x = (SK K) = S_K.eval (SK ([x] [u_Hd]))

        let contents = 
            S_K.create 
            (SK ((tree [u_Hd]) (tree [u_Tl]) (tree [u_Tl][u_Hd])))
            (LS tree)

        let left =
            S_K.create
            (SK (tree [u_Tl] [u_Tl] [u_Hd]))
            (LS tree)

        let right =
            S_K.create
            (SK (tree [u_Tl] [u_Tl] [u_Tl]))
            (LS tree)

        let join a b c =
            SK (Pair (S K) (Pair [c] (Pair [a] [b])))

        let rec to_tree x =
            let contents = (S_K.eval (SK ([contents] [x]))) in
            if ml_leaf_signal x then Leaf contents
            else 
                let left = S_K.eval (SK ([left] [x]))
                and right = S_K.eval (SK ([right] [x])) in
                Node (contents, to_tree left, to_tree right)

        let rec of_tree tree =
            match tree with
            | Leaf x -> leaf x
            | Node (x, a, b) -> join (of_tree a) (of_tree b) x
    end

    module type Integer = sig
        val zero : S_K.primitives
        val not_zero : S_K.primitives -> S_K.primitives
        val ml_zero_signal : S_K.primitives -> bool
        val predecessor : S_K.primitives -> S_K.primitives
        val successor : S_K.primitives -> S_K.primitives
        val to_int : S_K.primitives -> int
        val of_int : int -> S_K.primitives
    end

    module ChurchIntegers : Integer = struct
        (* A representation of integers using Church's original list
        * representation with as many K as there can be *)
        let zero = S_K.eval (SK (Pair (K S) K))
        let not_zero x = S_K.eval (SK ([x] [u_Hd] K K))
        let ml_zero_signal x = 
           (SK (S K)) = S_K.eval (SK ([x] [u_Hd] K K))
        let predecessor x = SK ([x] [u_Tl])
        let successor x = SK (Pair K [x])
        let to_int a = 
            let rec aux_to_int acc a =
                if ml_zero_signal a then acc
                else aux_to_int (acc + 1) (predecessor a)
            in
            aux_to_int 0 a
        let of_int a = 
            let rec aux_of_int acc a =
                if a = 0 then acc
                else aux_of_int (successor acc) (a - 1)
            in
            S_K.eval (aux_of_int zero a)
    end

    module Dna : List with type ocaml_repr = string with type ocaml_list =
        string list = struct
        (* A representation of DNA sequences *)
        type ocaml_repr = string
        type ocaml_list = string list
        let empty_list = S_K.eval (SK (Pair (Pair (K S)  K) K))
        let is_not_empty x = (SK ([x] [u_Hd] [u_Hd] K K))
        let ml_zero_signal x = 
           (SK (K S)) = S_K.eval (SK ([x] [u_Hd] [u_Hd]))
        let head x = SK ([x] [u_Hd])
        let tail x = SK ([x] [u_Tl])
        let prepend x y = SK (Pair [x] [y])
        let alphabet = 
            List.map (fun (a, b) -> S_K.eval a, b)
            [(SK (Pair K K), "A"); (SK (Pair K (S K)), "C"); 
            (SK (Pair (S K) K), "G"); (SK (Pair (S K) (S K)), "T")]
        let inv_alph = List.map (fun (a, b) -> (b, a)) alphabet
        let to_ocaml a = 
            match S_K.eval (SK (Hd [a])), S_K.eval (SK (Tl [a])) with
            | `K, `K -> "A"
            | `K, `Node [`S; `K] -> "C"
            | `Node [`S; `K], `K -> "G"
            | `Node [`S; `K], `Node [`S; `K] -> "T"
            | _ -> raise Not_found

        let of_ocaml a = List.assoc a inv_alph
        let complement x =
            let negate n = S_K.eval (SK ([x] [n] (S K) K)) in
            let a = negate u_Hd
            and b = negate u_Tl in
            S_K.eval (SK (Pair [a] [b]))
        let to_list lst = 
            let rec to_list lst =
                if not (ml_zero_signal (S_K.eval lst)) then 
                    let next = to_ocaml (S_K.eval (head lst)) in
                    next :: (to_list (tail lst))
                else []
            in
            to_list (S_K.eval lst) 
        let rec of_list lst = 
            List.fold_right (fun x acc ->
                prepend (of_ocaml x) acc) lst empty_list
    end

    module ComposableBase = struct
        let decode_function =
            let tmp_df =
                S_K.create 
                (SK ([Tree.is_leaf (SK encoder_in_use)]
                    ([Tree.contents] encoder_in_use acc starter
                    R original_f_encoder lst)
                    (R acc starter original_f_encoder ((lst [u_Hd]) ([Tree.left]
                    encoder_in_use) ([Tree.right] encoder_in_use)) (lst
                    [u_Tl]))))
                (LS R acc starter original_f_encoder encoder_in_use lst)
            in
            S_K.create (SK (Fixedpoint [tmp_df])) []

        (* starter R decode_function original_f_encoder acc lst 
        * decode_function R acc starter original_f_encoder encoder_in_use lst
        * insert R acc starter decode_function original_f_encoder lst
        * delete R acc starter decode_function original_f_encoder lst
        * substitute R acc starter decode_function original_f_encoder lst *)
        let starter =
            let tmp_starter =
                S_K.create 
                (SK (((lst [u_Hd]) K K) 
                    (decode_function acc R original_f_encoder original_f_encoder
                    lst)
                    acc))
                (LS R decode_function original_f_encoder acc lst)
            in
            S_K.create (SK (Fixedpoint [tmp_starter] [decode_function])) []


    end

    module type IntegerRepresentation = sig
        type extras
        (* A module type to convert from a given type to ChurchIntegers which
        * are easier to deal with *)
        val of_int : extras -> int -> S_K.primitives
        val to_int : extras -> S_K.primitives -> int
        val to_church : extras -> S_K.primitives
    end

    module ComposableAtomicSE = struct
        (* insert R pos base acc starter decode_function original_f_encoder lst
        * *)
        let aux_insert =
            let tmp_i =
                S_K.create 
                (SK ([ChurchIntegers.not_zero (SK pos)]
                    [Dna.prepend (Dna.head (SK acc)) 
                        (SK (R ([ChurchIntegers.predecessor (SK pos) ]) base 
                        [Dna.tail (SK acc)] starter
                        decode_function original_f_encoder lst))]
                    ([Dna.prepend (SK base) (SK (starter
                    decode_function original_f_encoder acc lst))])))
                (LS R pos base acc starter decode_function original_f_encoder
                lst)
            in
            S_K.create (SK (Fixedpoint [tmp_i])) []

        let aux_delete = 
            let tmp_d =
                S_K.create
                (SK ([ChurchIntegers.not_zero (SK pos)]
                    [Dna.prepend (SK [Dna.head (SK acc)]) (SK (R
                    ([ChurchIntegers.predecessor (SK pos)]) [Dna.tail (SK acc)] starter
                    decode_function original_f_encoder lst))]
                    (starter decode_function
                    original_f_encoder [Dna.tail (SK acc)] lst)))
                (LS R pos acc starter decode_function original_f_encoder
                lst)
            in
            S_K.create (SK (Fixedpoint [tmp_d])) []

        let aux_substitute =
            let tmp_s =
                S_K.create 
                (SK ([ChurchIntegers.not_zero (SK pos)]
                    [Dna.prepend (Dna.head (SK acc)) (SK (R
                    ([ChurchIntegers.predecessor (SK pos) ]) base [Dna.tail (SK
                    acc)] starter
                    decode_function original_f_encoder lst))]
                    ([Dna.prepend (SK base) (SK (starter
                    decode_function original_f_encoder [Dna.tail (SK acc)] lst))])))
                (LS R pos base acc starter decode_function original_f_encoder
                lst)
            in
            S_K.create (SK (Fixedpoint [tmp_s])) []

        let insert to_church = 
            S_K.create 
            (SK ([aux_insert] ([to_church] (lst [u_Hd])) (lst [u_Tl] [u_Hd]) acc starter
            decode_function original_f_encoder (lst [u_Tl] [u_Tl])))
            (LS acc starter decode_function original_f_encoder lst)

        let delete to_church =
            S_K.create 
            (SK ([aux_delete] ([to_church] (lst [u_Hd])) acc starter
            decode_function original_f_encoder (lst [u_Tl])))
            (LS acc starter decode_function original_f_encoder lst)

        let substitute to_church = 
            S_K.create 
            (SK ([aux_substitute] ([to_church] (lst [u_Hd])) (lst [u_Tl] [u_Hd]) acc starter
            decode_function original_f_encoder (lst [u_Tl] [u_Tl])))
            (LS acc starter decode_function original_f_encoder lst)

        let editions of_int a b =
            let len = String.length a in
            if len <> String.length b then failwith "Sequences not aligned"
            else
                let rec comparator pos delta =
                    if pos = len then (SK (Pair (K S) K))
                    else
                        let np = pos + 1 in
                        if a.[pos] = b.[pos] then comparator np (delta + 1) 
                        else if b.[pos] = '-' then 
                            (SK (Pair (S K) (Pair (S K) (Pair [of_int delta]
                            [comparator np 0]))))
                        else if a.[pos] = '-' then 
                            (SK (Pair (S K) (Pair K (Pair [of_int delta]
                            (Pair [Dna.of_ocaml (String.make 1 b.[pos])] 
                            [comparator np 0])))))
                        else 
                            (SK (Pair K (Pair [of_int delta]
                            (Pair [Dna.of_ocaml (String.make 1 b.[pos])] 
                            [comparator np 0]))))
                in
                comparator 0 0

    end

    module Chromosome : List with type ocaml_repr = Dna.ocaml_list with type
    ocaml_list = Dna.ocaml_list list = struct

        type ocaml_repr = string list
        type ocaml_list = string list list 

        let empty_list = SK (Pair [Dna.empty_list] K)

        let is_not_empty x = Dna.is_not_empty (SK ([x] [u_Hd]))

        let ml_zero_signal x = Dna.ml_zero_signal (SK ([x] [u_Hd]))

        let head x = SK ([x] [u_Hd])
        let tail x = SK ([x] [u_Tl])
        let prepend x y = SK (Pair [x] [y])
        let alphabet = []
        let to_ocaml item = Dna.to_list item
        let of_ocaml lst = Dna.of_list lst
        let rec to_list x = 
            if ml_zero_signal x then []
            else (to_ocaml (head x)) :: (to_list (tail x))
        let of_list x =
            List.fold_right (fun x acc ->
                SK (Pair [of_ocaml x] [acc])) x empty_list
        let complement _ = raise Unsupported
    end

    module Debugger (R : List) = struct
        let rec string_list (lst : S_K.primitives) =
            if R.ml_zero_signal lst then 
                (R.to_ocaml (S_K.eval (R.head lst))) :: 
                    string_list (S_K.eval (R.tail lst))
            else []
    end

    module type SE = sig 
        (* A sequence edition module that provides the functions required to
        * edit sequences of some predefined alphabet using the integer and
        * sequence representation of the programmer choice *)

        (** {2 Sequence Edition Functions}
         *
         * A sequence edition function is a machine in the SK language that
         * supports operations in sequences of a given alphabet. Different 
         * restrictions for representation and degree of functionality can 
         * be specified for each module implementing [SE].  *)

        (** [insert] inserts subsequences in another sequence. *)
        val insert : S_K.primitives

        (** [delete] deletes subsequences in another sequence. *)
        val delete : S_K.primitives

        (** [substitute] substitutes subseuqneces in another sequence. *)
        val substitute : S_K.primitives

        (** {2 Ocaml Interface}
         *
         * The following types and functions simplify the generation of SK
         * machines for particular [insert], [delete], and [substitute]
         * operations using the OCaml types directly. *)

        type insert_generator 
        type delete_generator 
        type substitute_generator
        val insert_generator : insert_generator
        val delete_generator : delete_generator
        val substitute_generator : substitute_generator

    end

    module AtomicSE (I : Integer) (S : List) : SE 
        with type insert_generator = int -> S.ocaml_repr -> S_K.primitives
        with type delete_generator = int -> S_K.primitives
        with type substitute_generator = int -> S.ocaml_repr -> S_K.primitives 
        = struct
        (* Atomic operations on sequences. The integer representation is that
        * provided by I, and the sequence representation is that provided by S. 
        * In AtomicSE, all the operations are atomic, in other words,
        * insertions, deletions, and substitutions occur in a single alphabet
        * element at a time basis. In other words, each occurs independently
        * along the sequence, one base at a time. *)

        (* A simplified function to generate all the atomic operations with
        * standard recursion until a zero condition is true *)
        let until_zero_recursions lst pos if_zero if_not_zero =
            let simple_recursion =
                S_K.create (SK 
                    ([I.not_zero (SK pos)]
                        [if_not_zero] 
                        [if_zero])) 
                lst
            in
            S_K.expand (SK (Fixedpoint [simple_recursion]))

        (** [insert p b s] inserts the element [b] belonging to
        * the alphabet of [S] in the position [p] as represented by the module
        * [I] in the sequence [s] as represented by the module [S]. The
        * positions is counted from the 0 index from left to right.*)
        let insert = 
            until_zero_recursions
            (LS R pos base seq)
            (SK pos)
            (S.prepend (SK base) (SK seq))
            (S.prepend (S.head (SK seq))
            (SK (R [I.predecessor (SK pos)] base [S.tail (SK seq)])))

        (** [delete p s] deletes the element located in position
        * [p] as represented by [I] of the sequence [s] as represented by [S].
        * The position is counted from the 0 index from left to right. *)
        let delete = 
            until_zero_recursions
            (LS R pos seq)
            (SK pos)
            (S.tail (SK seq))
            (S.prepend 
                (S.head (SK seq)) 
                (SK (R [I.predecessor (SK pos)] [S.tail (SK seq)])))

        (** [substitute p b s] substitutes the element located
        * in position [p] as represented by [I] with the element [b] from the 
        * alphabet supported by [S] in the sequence [s] as represented by [S].
        * The position is counted from the 0 index from left to right. *)
        let substitute =
            until_zero_recursions
            (LS R pos base seq)
            (SK pos)
            (S.prepend (SK base) (S.tail (SK seq)))
            (S.prepend (S.head (SK seq))
                (SK (R [I.predecessor (SK pos)] base [S.tail (SK seq)])))

        (** [insert_generator p b] takes the position [p] and element [b] to be 
        * inserted, and produces a machine that only requires the sequence on 
        * which to be inserted. *)
        type insert_generator = int -> S.ocaml_repr -> S_K.primitives

        (** [delete_generator p] takes the position [p] and produces a machine
        * that only requires the sequence on which that position should be
        * deleted. *)
        type delete_generator = int -> S_K.primitives

        (** [substitute_generator p b] takes the position [p] and element [b]
         * and generates a machine that only requires the sequence on which the
         * substitution will be performed. *)
        type substitute_generator = int -> S.ocaml_repr -> S_K.primitives

        let insert_generator pos base = 
            SK ([insert] [I.of_int pos] [S.of_ocaml base])

        let delete_generator pos =
            SK ([delete] [I.of_int pos])

        let substitute_generator pos base =
            SK ([substitute] [I.of_int pos] [S.of_ocaml base])
    end

    module AffineSE (I : Integer) (S : List) : SE 
        with type insert_generator = int -> S_K.primitives
        with type delete_generator = int -> int ->  S_K.primitives
        with type substitute_generator = int -> S_K.primitives
        = struct

        (** [insert p pref s] inserts the sequence [pref] as represented in [S]
         * in the position [p] as represented in [I] in the sequence [s] as
         * represented in [S]. *)
        let insert = 
            (* Instead of passing a single base, we will pass a
            complete sequence to this function *)
            let tmp = 
                S_K.create 
                (SK ([I.not_zero (SK pos)]
                    (* If we are not in position 0 we continue *)
                    [S.prepend (S.head (SK seq)) 
                        (SK (R [I.predecessor (SK pos)] prefix 
                        [S.tail (SK seq)]))]
                    (* We have reached the position, we better merge now *)
                    ([S.is_not_empty (SK prefix)]
                    (* If we have some base left in the sequence that we will 
                    * insert *)
                        [S.prepend (S.head (SK prefix)) 
                            (SK (R pos [S.tail (SK prefix)] seq))]
                            (* Otherwise we have finished inserting *)
                            seq)))
                (LS R pos prefix seq)
            in
            S_K.expand (SK (Fixedpoint [tmp]))

        (** [delete p len s] deletes a section of length [len] starting in
        * position [p] as represented in [I] in the sequence [s] as
         * represented in [S]. *)
        let delete =
            let tmp =
                S_K.create
                (SK ([I.not_zero (SK pos)] 
                    (* If we haven't reached the deletion pos *)
                    [S.prepend (S.head (SK seq))
                        (SK (R [I.predecessor (SK pos)] len 
                            [S.tail (SK seq)]))]
                    (* If we have reached the deletion pos *)
                    ([I.not_zero (SK len)]
                        (* We still have something more to del *)
                        (R pos [I.predecessor (SK len)] [S.tail (SK seq)])
                        (* We are done with the deletion *)
                        seq)))
                (LS R pos len seq)
            in
            S_K.expand (SK (Fixedpoint [tmp]))

        (** [substitute p pref s] substitutes the sequence [pref] as represented in [S]
         * in the position [p] as represented in [I] in the sequence [s] as
         * represented in [S]. *)
        let substitute =
            let tmp =
                S_K.create
                (SK ([I.not_zero (SK pos)] 
                (* If we haven't reached the deletion position *)
                    [S.prepend (S.head (SK seq)) 
                    (SK (R [I.predecessor (SK pos)] repl [S.tail (SK seq)]))]
                    (* We have reached the deletion position *)
                    ([S.is_not_empty (SK repl)] 
                        (* Still some more to substitute *)
                        [S.prepend (S.head (SK repl))
                        (SK (R pos [S.tail (SK repl)] [S.tail (SK seq)]))]
                        (* Nothing else to substitute *)
                        seq)))
                (LS R pos repl seq)
            in
            S_K.expand (SK (Fixedpoint [tmp]))

        (** [insert_generator pos] generates an insertion machine in position
        * [pos] requiring only as input the sequence to insert and the sequence
        * where the insertion is to be performed *)
        type insert_generator = int -> S_K.primitives 

        (** [delete_generator pos len] generates a deletion machine in position
        * [pos] with length [len] and requires only as input the sequence where
        * the deletion is to be performed *)
        type delete_generator = int -> int -> S_K.primitives

        (** [substitute_generator pos] does the same as [insert_generator] but
        * for a substitution event. *)
        type substitute_generator = int -> S_K.primitives 

        let insert_generator pos = 
            SK ([insert] [I.of_int pos])

        let delete_generator pos len =
            SK ([delete] [I.of_int pos] [I.of_int len])

        let substitute_generator pos =
            SK ([substitute] [I.of_int pos])
    end


    module LogInt : IntegerRepresentation with type extras = unit = struct
        type extras = unit

        let of_int () int =
            let rec encode_in_log int acc =
                if int = 0 then acc
                else 
                    encode_in_log (int lsr 1)
                    (SK (Pair [if 1 = (1 land int) then (SK K) else (SK (S K))]
                    [acc]))
            in
            if int < 0 then failwith "Illegal argument"
            else encode_in_log int (SK EmptyList)

        let to_int () x =
            let lst = K_S.to_list x in
            List.fold_left (fun acc item ->
                (2 * acc) +
                (match item with
                | `K -> 1 
                | `SK -> 0)) 0 lst

        let to_church () = SK DecodeInt 

    end

    module EncodedInteger : IntegerRepresentation with 
    type extras = int option Parser.Tree.t = struct
        type extras = int option Parser.Tree.t 

        let rec encoder tree =
            match tree with
            | Parser.Tree.Leaf (Some x) -> Tree.leaf (LogInt.of_int () x)
            | Parser.Tree.Node ([l;r ], None) ->
                    Tree.join (encoder l) (encoder r) (SK K)
            | _ -> failwith "Illegal tree"

        let decoder =
            let tmp =
                S_K.create 
                (SK ([Tree.is_leaf (SK encoder)] 
                ([Tree.contents] encoder)
                (R ([K_S.head (SK lst)] ([Tree.left] encoder) ([Tree.right]
                encoder)) [K_S.tail (SK lst)])))
                (LS R encoder lst)
            in
            S_K.create (SK (Fixedpoint [tmp])) []



        let of_int tree int =
            let rec aux_int tree =
                match tree with
                | Parser.Tree.Leaf (Some x) ->
                        if x = int then Some (SK K)
                        else None
                | Parser.Tree.Node ([l; r], None) ->
                        (match aux_int l with
                        | None ->
                                (match aux_int r with
                                | None -> None
                                | Some x -> Some (SK (Pair (S K) [x])))
                        | Some x -> Some (SK (Pair K [x])))
                | _ -> assert false
            in
            match aux_int tree with
            | None -> prerr_int int; prerr_newline (); assert false
            | Some x -> x

        let rec to_int tree int =
            match tree with
            | Parser.Tree.Leaf (Some x) -> x
            | Parser.Tree.Node ([l; r], None) ->
                    (match S_K.eval (SK (Hd [int])) with
                    | `K -> to_int l (SK (Tl [int])) 
                    | `Node [`S; `K] -> to_int r (SK (Tl [int]))
                    | _ -> assert false)
            | _ -> assert false

        let to_church tree =
            let enc = encoder tree in
            let log = LogInt.to_church () in
            S_K.create (SK ([log] ([decoder] [enc] int))) (LS int)

    end


    module FixedMax : IntegerRepresentation with type extras = int = struct

        type extras = int 

        let of_int max int =
            let rec encode_in_log max int acc =
                if max = 0 then acc
                else 
                    encode_in_log (max lsr 1) (int lsr 1) 
                    (SK (Pair [if 1 = (1 land int) then (SK K) else (SK (S K))]
            [acc]))
            in
            if int < 0 || max < 0 || max < int then failwith "Illegal argument"
            else encode_in_log max int (SK K)

        let to_church max =
            let rec aux_len max acc =
                if max = 0 then acc 
                else aux_len (max lsr 1) (acc + 1) 
            in
            if max < 0 then failwith "Illegal argument" 
            else
                let counter = ChurchIntegers.of_int (aux_len max 0) in
                let tmp_rec = 
                    S_K.create
                    (SK ([ChurchIntegers.not_zero (SK cnt)]
                        (R [ChurchIntegers.predecessor (SK cnt)] 
                        ((int [u_Hd]) 
                            (Add [ChurchIntegers.of_int 1] (Multiply
                            [ChurchIntegers.of_int 2] acc)) 
                            (Multiply [ChurchIntegers.of_int 2] acc))
                            (int [u_Tl]))
                        acc))
                    (LS R cnt acc int)
                in
                S_K.create 
                (SK (Fixedpoint [tmp_rec] [counter] [ChurchIntegers.zero])) 
                []

        let to_int max int = 
            ChurchIntegers.to_int 
            (S_K.eval (SK ([to_church max] [int])))

    end

    module FixedRange : IntegerRepresentation with type extras = (int * int) =
        struct
            type min = int
            type max = int
            type extras = (min * max)

            let of_int (min, max) int =
                if min > max || int < min || int > max then 
                    failwith "Illegal argument"
                else
                    FixedMax.of_int (max - min) int

            let to_int (min, max) int =
                min + (FixedMax.to_int (max - min) int)

            let to_church (min, max) =
                S_K.create
                (SK (Add [ChurchIntegers.of_int min] ([FixedMax.to_church (max -
                min)] pos))) (LS pos)
    end

    let convert_to_log_position x = 
        S_K.expand (S_K.create (SK ([x] (DecodeInt pos))) (LS pos))

    (* We are almost there, let's make a function that encodes an integer in
    * it's logarithmic binary representation *)
    let encode_in_log int =
        let rec encode_in_log int acc =
            if int = 0 then acc
            else 
                encode_in_log (int lsr 1)
                (SK (Pair [if 1 = (1 land int) then (SK K) else (SK (S K))]
                [acc]))
        in
        if int < 0 then failwith "Illegal argument"
        else encode_in_log int (SK EmptyList)

    module AtomicSE_LogInts (S : SE) : SE
        with type insert_generator = int -> S_K.primitives
        with type delete_generator = int -> S_K.primitives
        with type substitute_generator = int -> S_K.primitives 
        = struct

        let generic f = S_K.create (SK ([f] ([LogInt.to_church ()] pos))) (LS pos)
        let insert = generic S.insert
        let delete = generic S.delete
        let substitute = generic S.substitute

        type insert_generator = int -> S_K.primitives 
        type delete_generator = int -> S_K.primitives
        type substitute_generator = int -> S_K.primitives 

        let insert_generator pos = 
            SK ([insert] [LogInt.of_int () pos])

        let delete_generator pos =
            SK ([delete] [LogInt.of_int () pos])

        let substitute_generator pos =
            SK ([substitute] [LogInt.of_int () pos])
    end

    module ComposableSE (L : List) (I : IntegerRepresentation) = struct
        let starter =
            let tmp =
                S_K.create
                (SK (((lst [u_Hd]) K K) (decoder acc R lst) (acc)))
                (LS acc lst decoder)
            in
            S_K.create (SK (Fixedpoint [tmp]))

        let decoder =
            let tmp =
                S_K.create
                (SK ([Tree.is_leaf (SK tree_in_use)] 
                    (([Tree.contents] tree_in_use) 
                    (R starter original_tree original_tree)
                        starter lst)
                    (R starter original_tree 
                        ((lst [u_Hd]) ([Tree.left] tree_in_use) 
                            ([Tree.right] tree_in_use)) (lst [u_Tl]) acc)))
                (LS starter original_tree tree_in_use lst acc)
            in
            S_K.create (SK (Fixedpoint [tmp])) []
    end

    module AffineSE_LogInts (S : SE) : SE 
        with type insert_generator = int -> S_K.primitives
        with type delete_generator = int -> int -> S_K.primitives
        with type substitute_generator = int -> S_K.primitives 
        = struct
        let generic f = S_K.create (SK ([f] ([LogInt.to_church ()] pos))) (LS pos)
        let insert = generic S.insert
        let delete = 
            S_K.create (SK ([S.delete] ([LogInt.to_church ()] pos)
            ([LogInt.to_church ()] len)))
            (LS pos len)
        let substitute = generic S.substitute

        type insert_generator = int -> S_K.primitives 
        type delete_generator = int -> int -> S_K.primitives
        type substitute_generator = int -> S_K.primitives 

        let insert_generator pos = 
            SK ([insert] [LogInt.of_int () pos])

        let delete_generator pos len =
            SK ([delete] [LogInt.of_int () pos] 
            [encode_in_log len])

        let substitute_generator pos =
            SK ([substitute] [LogInt.of_int () pos])
    end

    module type HO = sig
        include SE
        val translocate : S_K.primitives
        val invert : S_K.primitives
        val duplicate : S_K.primitives
        val tandem : S_K.primitives
    end

    module type GO = sig
        include HO
        val chr_translocation : S_K.primitives
        val chr_crossover : S_K.primitives
    end

    module HighOrder (I : Integer) (S : List) : HO = 
        struct
        module Aff = AffineSE (I) (S)

        include Aff 

        let extract_subsequence =
            let tmp =
                S_K.create (SK ([I.not_zero (SK start)]
                    (R [I.predecessor (SK start)] len [S.tail (SK seq)])
                    ([I.not_zero (SK len)]
                        ([S.prepend (S.head (SK seq)) 
                            (SK (R [I.zero] [I.predecessor (SK len)] 
                            [S.tail (SK seq)]))])
                        [S.empty_list])))
                (LS R start len seq)
            in
            S_K.create (SK (Fixedpoint [tmp])) []

        let invert_seq = 
            let tmp =
                S_K.create 
                (SK ([S.is_not_empty (SK seq)]
                    (R ([S.prepend (S.complement (S.head (SK seq))) (SK acc)])
                    ([S.tail (SK seq)]))
                    acc))
                (LS R acc seq)
            in
            S_K.create (SK (Fixedpoint [tmp])) []

        let translocate =
            S_K.create
            (SK ([Aff.insert] ins_pos ([extract_subsequence] tran_pos len seq) 
            ([Aff.delete] tran_pos len seq)))
            (LS ins_pos tran_pos len seq)

        let duplicate = 
            S_K.create
            (SK ([Aff.insert] ins_pos ([extract_subsequence] tran_pos len seq)
            seq))
            (LS ins_pos tran_pos len seq)

        let invert =
            S_K.create 
            (SK ([Aff.insert] pos ([invert_seq] [S.empty_list] ([extract_subsequence]
            pos len seq)) ([Aff.delete] pos len seq)))
            (LS pos len seq)

        let tandem =
            let tmp =
                S_K.create
                (SK ([I.not_zero (SK rep)]
                    ([duplicate] pos pos len (R [I.predecessor (SK rep)] pos len seq ))
                    (seq)))
                (LS R rep pos len seq)
            in
            S_K.create (SK (Fixedpoint [tmp])) []
    end

(*
    (* A encoded is made of pairs in the following way:
        * (Pair K (Pair encoded encoded))
        * (Pair (K S) Result)
        * Where left (K) or right (S K) of the right element represent the
        * actual path to be followed from the list. *)
    let decoder simplified =
        let next_encoder = (SK ((encoder [u_Tl]) (lst [u_Hd]))) in
        let decoder_recursion =
            S_K.create 
            (SK ([next_encoder] [u_Hd] K K (R (lst [u_Tl]) [next_encoder]) 
            [if simplified then
                (SK ([next_encoder] [u_Tl]))
            else
                (SK (Pair ([next_encoder] [u_Tl]) (lst [u_Tl])))]))
            (LS R lst encoder)
        in
        S_K.expand (SK (Fixedpoint [decoder_recursion]))

    let rec make_encoder tree leafs =
        match tree with
        | `Node [a; b] -> 
                (SK (Pair K (Pair [make_encoder a leafs] [make_encoder b
                leafs])))
        | `Label x -> 
                (SK (Pair (K S) [try (List.assoc x leafs) with Not_found ->
                    tree]))
        | `S -> (SK (Pair (K S) S))
        | `K -> (SK (Pair (K S) K))
        | `Node _ -> failwith "Illegal encoder: each node must have two leafs"

    (* An example of encoder and decoder usage *)
    let encoder_example = make_encoder (SK (substitute (delete insert))) []

    let result_example = 
        let decoder = S_K.expand (decoder false) in
        S_K.reduce (S_K.expand (SK (([decoder] [(SKLS (S K) (S K))] [encoder_example])
        [u_Hd])))

    (* Each function is specified the number of items that need to be read. The
    * accumulator for all the calculations is always the last item to be added.
    * As the number of arguments is so small, we are better off specifying them
    * in the simplest way, not using the logarithmic, but the standard notation.
    * 
    * So our representation consists of a list of items that are read one after
    * the other and are fed to the encoder and the sequence representation.
    * *)

    (* The following concept:
        * fun f_apply decoder acc list =
            * fun apply acc list =
                * if list is empty then acc
                * else apply (acc (hd list)) (tl list) 
            * in
    *         ((apply (decoder (hd list)) (tl list)) acc)
    *  In order to use it, we need to have a tree (made of pairs) of the following 
    *  form:
        *  (encoded_function, list_of_arguments)
        *  Where the list of arguments is
        *  (K, (Argument, list_of_arguments))
        *  ((K S), (K S)) if we have reach the end of the arguments.
    *         *)
    let f_apply =
        let apply =
            let tmp =
                S_K.create
                (SK ((lst [u_Hd] K K) (R (acc (lst [u_Tl] [u_Hd])) (lst [u_Tl]
                [u_Tl])) acc))
                (LS R acc lst)
            in
            S_K.create (SK (Fixedpoint [tmp])) []
        in
        S_K.create (SK (([apply] ([decoder false] (lst [u_Hd]) encoder  [u_Hd]) (lst [u_Tl])) acc)) 
        (LS acc encoder lst)

    (* An example using it: 
    let () = 
        let r =  
            S_K.eval (SK ([f_apply] [SKLS [adenine]] [make_encoder (SK (Delete
            (Insert Substitute))) ["Delete", delete; "Insert", insert; "Substitute",
            substitute]] [SKLS [SKLS K] K EmptyList])) in
        ()
    *)

    let swap = S_K.create (SK (a b)) (LS b a)

    (* Now we can define the edition of sequences using a composition of
    * editions. In other words, we don't need to specify the absolute position
    * of a base for a deletion, but we can specify the relative position from
    * the previous edition. This could produce a gain if the editions are dense
    * enough. To do this, we will define another function similar to the
    * [f_apply] that we just defined. *)

    let apply_list_of_functions_on_accumulator =
        let recursion = 
            S_K.create
            (SK ((todo [u_Hd] K K)  (* We are not done yet *)
                ([decoder true] (todo [u_Tl] [u_Hd]) encoder R acc encoder 
                (todo [u_Tl] [u_Tl]))
                (acc)))
            (LS R acc encoder todo)
        in
        S_K.expand (SK (Fixedpoint [recursion]))

    let insert_composable = 
        let insert_composable =
            until_zero_recursions
            (LS R Continue pos base seq encoder todo)
            (SK pos)
            (SK (Pair base (Continue seq encoder todo)))
            (SK (Pair (seq [u_Hd]) (R Continue (pos [u_Tl]) base (seq [u_Tl])
            encoder todo)))
        in
        let pos = SK (DecodeInt (todo [u_Hd]))
        and base = SK (todo [u_Tl] [u_Hd]) 
        and rest = SK (todo [u_Tl] [u_Tl]) in
        S_K.create
        (SK ([insert_composable] Continue [pos] [base] seq encoder [rest]))
        (LS Continue seq encoder todo)

    let delete_composable =
        let delete_composable = 
            until_zero_recursions
            (LS R Continue pos seq encoder todo)
            (SK pos)
            (SK (Continue (seq [u_Tl]) encoder todo))
            (SK (Pair (seq [u_Hd]) (R Continue (pos [u_Tl]) (seq [u_Tl]) encoder
            todo)))
        in
        let pos = SK (DecodeInt (todo [u_Hd]))
        and rest = SK (todo [u_Tl]) in
        S_K.create
        (SK ([delete_composable] Continue [pos] seq encoder [rest])) 
        (LS Continue seq encoder todo)

    let substitute_composable =
        let substitute_composable = 
            until_zero_recursions
            (LS R Continue pos base seq encoder todo)
            (SK pos)
            (SK (Pair base (Continue (seq [u_Tl]) encoder todo)))
            (SK (Pair (seq [u_Hd]) (R Continue (pos [u_Tl]) base (seq [u_Tl])
            encoder todo)))
        in
        let pos = SK (DecodeInt (todo [u_Hd]))
        and base = SK (todo [u_Tl] [u_Hd]) 
        and rest = SK (todo [u_Tl] [u_Tl]) in
        S_K.create
        (SK ([substitute_composable] Continue [pos] [base] seq encoder [rest]))
        (LS Continue seq encoder todo)

    let is_not_empty_sequence x = SK ([x] [u_Hd] [u_Hd] K K)

    let empty_sequence = SK (Pair EmptyList EmptyList)

    (* A function to merge a pair of sequences. The only condition is that the
    * empty sequence is represented by a pair of `Label "EmptySequence". *)
    let merge_sequences =
        let merge_recursive =
            S_K.create
            (SK ([is_not_empty_sequence (SK prefix)]
            (* We have a base *)
            (Pair (prefix [u_Hd]) (R (prefix [u_Tl]) sufix))
            sufix))
            (LS R prefix sufix) 
        in
        S_K.create (SK (Fixedpoint [merge_recursive])) []

    let inversion_without_complement = (* A standard inversion function *)
        let inversion_recursive =
            S_K.create
            (SK ([is_not_empty_sequence (SK original)]
            (* We still have to continue inverting *)
            (R (Pair (original [u_Hd]) result) (original [u_Tl]))
            original))
            (LS R result original)
        in
        S_K.create (SK (Fixedpoint [inversion_recursive])) []

    let inversion_with_complement =
        let complement_base x = 
            SK (Pair (Not ([x] [u_Hd])) (Not ([x] [u_Tl])))
        in
        let inversion_recursive =
            S_K.create
            (SK ([is_not_empty_sequence (SK original)]
            (* We still have to continue inverting *)
            (R (Pair [complement_base (SK (original [u_Hd]))] result) 
            (original [u_Tl])) original))
            (LS R result original)
        in
        S_K.create (SK (Fixedpoint [inversion_recursive])) []

    let extract_subsequence =
        let extract_recursive =
            S_K.create
            (SK ((pos [u_Hd] K K) (* We are not in the 0 position of the
            extraction *)
                (R (pos [u_Tl]) len (seq [u_Tl]))
                ((len [u_Hd] K K) (* We have more things to append *)
                    (Pair (seq [u_Hd]) (R pos (len [u_Tl]) (seq [u_Tl])))
                    [empty_sequence])))
            (LS R pos len seq)
        in
        S_K.create (SK (Fixedpoint [extract_recursive])) []

    let duplicate = 
        S_K.create
        (SK ([affine_insert] pos ([extract_subsequence] beg len seq) seq))
        (LS pos beg len seq)

    let translocate =
        S_K.create
        (SK ([affine_insert] pos ([extract_subsequence] beg len seq)
        ([affine_delete] pos len seq)))
        (LS pos beg len seq)

    (* Now take a list of elements that can be passed to [lst] of f_apply and
    * apply on each one of the [f_apply]. This way we can perform a row of edits
    * in an accumulator that is passed around constantly. *)
    let f_edit =
        let tmp = S_K.create
            (SK ((lst [u_Hd] K K) 
            (R encoder (lst [u_Tl] [u_Tl]) ([f_apply] acc encoder (lst [u_Tl]
            [u_Hd]))) acc)) 
            (LS R encoder lst acc)
        in
        S_K.create (SK (Fixedpoint [tmp])) []

    (* We are almost there, let's make a function that encodes an integer in
    * it's logarithmic binary representation *)
    let encode_in_log int =
        let rec encode_in_log int acc =
            if int = 0 then acc
            else 
                encode_in_log (int lsr 1)
                (SK (Pair [if 1 = (1 land int) then (SK K) else (SK (S K))]
                [acc]))
        in
        if int < 0 then failwith "Illegal argument"
        else encode_in_log int (SK EmptyList)

    (* Finally we can do things with sequences *)
    let pm_encodings = 
        make_encoder 
        (SK (substitute (delete insert)))
        ["substitute", substitute; "delete", delete; "insert", insert]

    (* An example using the previous composable edition functions  *)
    let () =
        let do_test = false in
        if do_test then
            let encoder = 
               make_encoder (SK (substitute (insert delete)))
               ["substitute", substitute_composable; "delete", delete_composable;
               "insert", insert_composable] 
            in
            let todo = 
                SKLS K [SKLS K] [encode_in_log 0] [timine] K [SKLS (S K) (S K)]
                [encode_in_log 1]
            in
            let seq = SKLS [adenine] [guanine] [guanine] [guanine] in
            let res = 
                S_K.eval 
                (SK ([apply_list_of_functions_on_accumulator] 
                [seq] [encoder] [todo]))
            in 
            print_endline (S_K.to_string res)
        else ()

    let create_sequence a =
        let pair = `Label "Pair" in
        let merge a b = `Node [pair; a; b] in
        let rec prepend pos acc =
            if pos < 0 then acc
            else 
                let npos = pos - 1 in
                match a.[pos] with
                | 'A' | 'a' -> prepend npos (merge adenine acc)
                | 'C' | 'c' -> prepend npos (merge citosine acc)
                | 'G' | 'g' -> prepend npos (merge guanine acc)
                | 'T' | 't' -> prepend npos (merge timine acc)
                | _ -> prepend npos acc
        in
        prepend ((String.length a) - 1) (SK (Pair EmptyList EmptyList))

    let generate_editions_machine_from_alignment a b =
        (* We want to produce a machine that given the sequence [a] produces as
        * output the sequence [b]. *)
        assert ((String.length a) = (String.length b));
        let max = String.length a in
        let base = function
            | 'A' | 'a' -> adenine
            | 'C' | 'c' -> citosine
            | 'G' | 'g' -> guanine
            | 'T' | 'U' | 't' | 'u' -> timine
            | '_' -> `Label "Emptylist"
            | _ -> failwith "Illegal DNA sequence"
        in
        let rec compare_and_generate_edition delta position =
            if max = position then []
            else
                let npos = position + 1 in
                if a.[position] = b.[position] then 
                    compare_and_generate_edition (delta + 1) npos
                else 
                    let deltap = encode_in_log delta 
                    and basep = base b.[position] in
                    if a.[position] = '_' then
                        (SK K) :: (SKLS (S K) K) :: deltap :: basep ::
                            compare_and_generate_edition 0 npos
                    else if b.[position] = '_' then
                        (SK K) :: (SKLS (S K) (S K)) :: deltap ::
                            compare_and_generate_edition 0 npos
                    else 
                        (SK K) :: (SKLS K) :: deltap :: basep ::
                            compare_and_generate_edition 0 npos
        in
        let pair = `Label "Pair" in
        let rec produce_list_expression lst = 
            match lst with
            | h :: t -> `Node [pair; h; produce_list_expression t]
            | [] -> `Label "EmptyList"
        in
        produce_list_expression (compare_and_generate_edition 0 0) 

    let generate_alignment_machine a b =
        let edition_machine = generate_editions_machine_from_alignment a b 
        and accumulator = create_sequence a 
        and encoder = 
           make_encoder (SK (substitute (insert delete)))
           ["substitute", substitute_composable; "delete", delete_composable;
           "insert", insert_composable] 
        in
        S_K.expand (SK ([apply_list_of_functions_on_accumulator] [accumulator]
        [encoder] [edition_machine]))

    let generate_acc_and_edition_machine a b =
        let edition_machine = generate_editions_machine_from_alignment a b 
        and accumulator = create_sequence a in
        S_K.expand edition_machine, S_K.expand accumulator


    (* A tree is represented the same as before:
        * An internal node
        * (K ((left_tree, left_editions) (right_tree, right_editions))) 
        * A leaf
        * (EmptyList) *)
    let process_tree =
        let recursive_processor =
            let left_editions = SK (tree [u_Tl] [u_Hd] [u_Tl])
            and right_editions = SK (tree [u_Tl] [u_Tl] [u_Tl])
            and left_tree = SK (tree [u_Tl] [u_Hd] [u_Hd])
            and right_tree = SK (Tree [u_Tl] [u_Tl] [u_Hd]) in
            S_K.create 
            (SK ((tree [u_Hd] K K) 
                (* We are in an interior node *)
                (Pair K (Pair (R processor (processor acc encoder [left_editions])
                encoder [left_tree]) (R processor (processor acc encoder [right_editions])
                encoder [right_tree])))
                (* No we are in a leaf *)
                (Pair (K S) acc)))
            (LS processor acc encoder tree)
        in
        S_K.create (SK (Fixedpoint [recursive_processor])) []

    (* Now we should be able to take a POY diagnosis and generate the complete
    * sequences, right? Let's give it a shot! *)
*)

end

module Compiler = struct

    type 'a kolmo_function =
        [ `Module of (name * 'a kolmo_function list)
        | `LetVal of (name * arguments * 'a definition)
        | `RecVal of (name * arguments * 'a definition) 
        | `Let of (name * arguments * 'a definition)
        | `Rec of (name * arguments * 'a definition) ]
    and name = string
    and arguments = string list 
    and 'a definition =
        [ `Value of 'a values
        | `IfElse of ('a condition * 'a definition * 'a definition)
        | `Letin of ('a kolmo_function * 'a definition) ]
    and 'a condition = [ `Condition of 'a definition ]
    and 'a values =
        [ `Apply of (string list * string * ('a definition list))
        | `Integer of string 
        | `Expr of 'a ]

    type sk_function = S_K.primitives kolmo_function

    type extended_sk_function = 
        [ S_K.primitives | `Identity | `Code of int ] kolmo_function

    type environment_item = 
        | Argument of string 
        | S_K of S_K.primitives
        | SKModule of environment_item All_sets.StringMap.t

    type decoder_environment_item = 
        | DArgument of string
        | DS_K of (S_K.primitives * int)
        | NS_K of S_K.primitives
        | DSKModule of decoder_environment_item All_sets.StringMap.t

    let fixedpoint = (SK (S S K (S (K (S S(S(S S K))))K)))

    let add_function env isrecursive name arguments definition =
        match env with
        | h :: t ->
                let f =
                    if not isrecursive then
                        S_K.create definition arguments 
                    else 
                        let tmp = S_K.create definition (name :: arguments) in
                        S_K.create (SK ([fixedpoint] [tmp])) []
                in
                (All_sets.StringMap.add name (S_K f) h) :: t
        | _ -> assert false

    let add_function_decoder counter isval env isrecursive name arguments 
    definition =
        match env with
        | h :: t ->
                let arguments = "__decoder" :: arguments in
                let f =
                    if not isrecursive then
                        S_K.create definition arguments 
                    else 
                        let tmp = S_K.create definition 
                            (name :: arguments) in
                        S_K.create (SK ([fixedpoint] [tmp])) []
                in
                let toadd, counter = 
                    if not isval then NS_K f, counter
                    else DS_K (f, counter), counter + 1
                in
                ((All_sets.StringMap.add name toadd h) :: t), counter
        | _ -> assert false

    let add_args f environment arg =
        match environment with
        | h :: t ->
                (All_sets.StringMap.add arg (f arg) h) :: t
        | _ -> assert false

    let add_args_decoder environment arg =
        add_args (fun x -> DArgument arg) environment arg

    let add_args environment arg = 
        add_args (fun x -> Argument arg) environment arg

    let rec of_int x =
        let m_true = (SK K) in 
        let m_false = (SK (S K)) in
        let m_pair = S_K.create (SK (c a b)) (LS a b c) in
        let m_church_zero = (SK ([m_pair] [m_false] K)) in
        let s_church_successor x = (SK ([m_pair] [m_true] [x])) in
        if x < 0 then failwith "Invalid argument" 
        else
            if x = 0 then m_church_zero
            else s_church_successor (of_int (x - 1))

    let rec find_in_environment s environment = 
        match environment with
        | [] -> raise Not_found
        | h :: t ->
                try All_sets.StringMap.find s h with
                | Not_found -> find_in_environment s t

    type defined = Always | AtTimes | Never

    let rec create_sk must_be_defined arguments environment definition =
        match definition with
        | `Value v -> 
                (match v with
                | `Apply (modules, s, args) -> 
                        (try 
                            let args = 
                                List.map 
                                (create_sk must_be_defined arguments environment)
                                args
                            in
                            let modenv = 
                                List.fold_left (fun env name ->
                                    match find_in_environment name env with
                                    | SKModule env -> env :: environment
                                    | _ -> assert false) environment modules
                            in
                            match must_be_defined with
                            | Never -> `Label s
                            | Always | AtTimes ->
                                    match find_in_environment s modenv with
                                    | Argument s -> 
                                            (match args with
                                            | [] -> `Label s
                                            | _ -> `Node ((`Label s) :: args))
                                    | S_K s -> 
                                            (match args with
                                            | [] -> s
                                            | _ -> `Node (s :: args))
                                    | SKModule _ -> raise Not_found
                        with Not_found ->
                            if must_be_defined = Always then
                                let name =
                                    match modules with
                                    | [] -> s
                                    | _ -> (String.concat "." modules) ^ "." ^ s
                                in
                                failwith ("Unbound value " ^ name)
                            else `Label s)
                | `Integer s ->
                        of_int (int_of_string s)
                | `Expr x -> x)
        | `IfElse (`Condition cond, a, b) -> 
                let cond = create_sk must_be_defined arguments environment cond
                and a = create_sk must_be_defined arguments environment a
                and b = create_sk must_be_defined arguments environment b in
                (SK ([cond] [a] [b]))
        | `Letin (new_fun, def) -> 
                let env = compile must_be_defined environment new_fun in
                create_sk must_be_defined arguments env def
    and compile must_be_defined environment (spec : S_K.primitives kolmo_function)  =
        match spec with
        | `Module (modname, lst) ->
                let modenv = 
                    List.fold_left (compile must_be_defined) (All_sets.StringMap.empty ::
                        environment) lst
                in
                (match modenv with
                | f :: s :: t ->
                        let s = 
                            All_sets.StringMap.add modname (SKModule f) s
                        in
                        s :: t
                | _ -> assert false)
        | `LetVal _ | `RecVal _ | `Let _ | `Rec _ as spec ->
                let (name, arguments, definition), isrecursive =
                    match spec with
                    | `LetVal x 
                    | `Let x -> x, false
                    | `RecVal x 
                    | `Rec x -> x, true
                in
                let definition = 
                    let new_env = 
                        let arguments = 
                            if isrecursive then name :: arguments
                            else arguments
                        in
                        List.fold_left add_args environment arguments in
                    create_sk must_be_defined arguments new_env definition 
                in
                add_function environment isrecursive name arguments definition

    let environment = ref [All_sets.StringMap.empty]

    type encoded = 
        | Encoded of int
        | UseEncoded
        | Inline
        | Module of encoded All_sets.StringMap.t

    let add_to_environment name item env =
        match env with
        | h :: t ->
                if All_sets.StringMap.mem name h then
                    failwith "Name already exists"
                else (All_sets.StringMap.add name item h) :: t
        | [] -> assert false

    let is_name_encoded path name environment =
        try
            let env = 
                match path with
                | h :: t ->
                        let env = 
                            List.find (All_sets.StringMap.mem h) environment 
                        in
                        List.fold_left (fun acc x -> 
                            match All_sets.StringMap.find x acc with
                            | Module x -> x | _ -> assert false) env t 
                | [] -> 
                        List.find (All_sets.StringMap.mem name) environment
            in
            match All_sets.StringMap.find name env with
            | Encoded _ | UseEncoded -> true
            | Inline -> false
            | Module _ -> assert false
        with
        | Not_found -> false

    let rec has_encoded arguments definition environment counter =
        match definition with
        | `IfElse (`Condition cond, a, b) ->
                has_encoded arguments cond environment counter ||
                has_encoded arguments a environment counter ||
                has_encoded arguments b environment counter
        | `Letin (new_fun, def) ->
                let env, counter = classify (environment, counter) new_fun in
                has_encoded arguments def env counter
        | `Value v -> 
                match v with
                | `Integer _ 
                | `Expr _ -> false
                | `Apply (path, name, args) ->
                        is_name_encoded path name environment ||
                        List.exists (fun x -> 
                            has_encoded arguments x environment counter) args

    and classify (environment, counter) (f : sk_function) = 
        match f with
        | `Module (name, fs) ->
                let res, counter = 
                    List.fold_left classify 
                    (All_sets.StringMap.empty :: environment, counter)
                    fs
                in
                (match res with
                | h :: _ ->
                        (add_to_environment name (Module h) environment),
                        counter
                | [] -> assert false)
        | `LetVal (name, arguments, definition) 
        | `RecVal (name, arguments, definition) ->
                (add_to_environment name (Encoded counter) environment), 
                counter + 1
        | `Rec (name, arguments, definition)
        | `Let (name, arguments, definition) as x ->
                let arguments =
                    match x with
                    | `Rec _ -> name :: arguments
                    | `Let _ -> arguments
                in
                let to_add =
                    if has_encoded arguments definition environment counter then 
                        UseEncoded
                    else Inline
                in
                (add_to_environment name to_add environment), counter

    let identity = "__identity"

    let rec create_sk_decoder counter must_be_defined arguments environment 
    definition =
        match definition with
        | `Value v -> 
                (match v with
                | `Apply (modules, s, args) -> 
                        (try 
                            let args = 
                                List.map 
                                (fun x ->
                                    fst (create_sk_decoder counter 
                                    must_be_defined arguments environment x)) 
                                args
                            in
                            let modenv = 
                                List.fold_left (fun env name ->
                                    match find_in_environment name env with
                                    | DSKModule env -> env :: environment
                                    | _ -> assert false) environment modules
                            in
                            match must_be_defined with
                            | Never -> `Label s, counter
                            | Always | AtTimes ->
                                    (match find_in_environment s modenv with
                                    | DArgument s -> 
                                            (match args with
                                            | [] -> `Label s
                                            | _ -> `Node ((`Label s) :: args))
                                    | NS_K s -> 
                                            (match args with
                                            | [] -> s
                                            | _ -> `Node (s :: args))
                                    | DS_K (_, s) ->
                                            `Node 
                                            ([`Label "__decoder"; 
                                                `Label (string_of_int s)] @ 
                                                args @ [`Label identity])
                                    | DSKModule _ -> raise Not_found), counter
                        with Not_found ->
                            if must_be_defined = Always then
                                let name =
                                    match modules with
                                    | [] -> s
                                    | _ -> (String.concat "." modules) ^ "." ^ s
                                in
                                failwith ("Unbound value " ^ name)
                            else `Label s, counter)
                | `Integer s ->
                        of_int (int_of_string s), counter
                | `Expr x -> x, counter)
        | `IfElse (`Condition cond, a, b) -> 
                let caller counter x =
                    create_sk_decoder counter must_be_defined arguments
                    environment x
                in
                let cond, counter = caller counter cond in
                let a, counter = caller counter a in
                let b, counter = caller counter b in
                (SK ([cond] [a] [b])), counter
        | `Letin (new_fun, def) -> 
                let env, counter = 
                    compile_decoder counter must_be_defined environment new_fun
                in
                create_sk_decoder counter must_be_defined arguments env def

    and compile_decoder counter must_be_defined environment (spec : sk_function)  =
        match spec with
        | `Module (modname, lst) ->
                let modenv, counter = 
                    List.fold_left (fun (env, counter) f -> 
                        compile_decoder counter must_be_defined env f) 
                    ((All_sets.StringMap.empty :: environment), counter) lst
                in
                (match modenv with
                | f :: s :: t ->
                        let s = 
                            All_sets.StringMap.add modname (DSKModule f) s
                        in
                        (s :: t), counter
                | _ -> assert false)
        | `LetVal _ | `RecVal _ | `Let _ | `Rec _ as spec ->
                let (name, arguments, definition), isrecursive, isval =
                    match spec with
                    | `LetVal x -> x, false, true
                    | `Let x -> x, false, false
                    | `RecVal x -> x, true, true
                    | `Rec x -> x, true, false
                in
                let definition, counter = 
                    let new_env = 
                        let arguments = 
                            if isrecursive then name :: "__decoder" :: arguments
                            else "__decoder" :: arguments
                        in
                        List.fold_left add_args_decoder environment arguments 
                    in
                    create_sk_decoder counter must_be_defined arguments new_env 
                    definition 
                in
                add_function_decoder counter isval environment isrecursive 
                name arguments definition

    let decoder = ref []
    let final_code = ref All_sets.IntegerMap.empty

    let add n code cnt =
        if All_sets.IntegerMap.mem code cnt then 
            All_sets.IntegerMap.add code 
            ((All_sets.IntegerMap.find code cnt) + n) cnt
        else All_sets.IntegerMap.add code n cnt

    let rec count_occurrences initial_frequencies identity_code env res =
        let rec internal_counter cnt x =
            match x with
            | `S | `K | `Debugger _ | `Lazy _ -> cnt
            | `Label x -> 
                    if x = identity then 
                        add 1 identity_code cnt
                    else 
                        let code = int_of_string x in 
                        add 1 code cnt
            | `Node lst -> List.fold_left internal_counter cnt lst
        in
        All_sets.StringMap.fold (fun name x res ->
            match x with
            | DArgument _ -> res
            | DSKModule m -> count_occurrences initial_frequencies identity_code m res
            | NS_K prim -> internal_counter res prim 
            | DS_K (prim, code) ->
                        let initial = 
                            try List.assoc name initial_frequencies with
                            | Not_found -> 1
                        in
                        internal_counter (add initial code res) prim) env res

    let rec replace_labels ?(prefix="") replacer env =
        let rec internal_replacer (x : S_K.primitives) : S_K.primitives list = 
            match x with
            | `S | `K | `Debugger _ | `Lazy _ -> [x]
            | `Label x -> replacer x
            | `Node lst ->
                    let lst = List.map internal_replacer lst in
                    let lst = List.flatten lst in
                    [`Node lst]
        in
        All_sets.StringMap.fold (fun name x res ->
            match x with
            | DArgument _ -> res
            | DSKModule m -> 
                    let prefix = 
                        if prefix = "" then name
                        else prefix ^ "." ^ name 
                    in
                    All_sets.StringMap.add
                    name (SKModule (replace_labels ~prefix replacer m)) res
            | DS_K (prim, c) ->
                    let f = 
                        match internal_replacer prim with
                        | [x] -> x
                        | x -> (`Node x)
                    in
                    let () = 
                        let name = 
                            if prefix = "" then name
                            else prefix ^ "." ^ name 
                        in
                        Printf.printf "Adding %s\n%!" name;
                        decoder := (name, c, f) :: !decoder;
                    in
                    All_sets.StringMap.add name (S_K f) res
            | NS_K prim -> 
                    All_sets.StringMap.add
                    name
                    (match internal_replacer prim with
                    | [x] -> S_K x
                    | x -> S_K (`Node x)) res)
        env All_sets.StringMap.empty

    let get_count a =
        match a with
        | Parser.Tree.Leaf (_, a)
        | Parser.Tree.Node (_, (_, a)) -> a

    module OrderedPair = struct
        type t = (int * int) Parser.Tree.t
        let compare a b = (get_count a) - (get_count b)
    end

    module MergeHeap = Heap.Make (OrderedPair)

    let getem heap =
        MergeHeap.findMin heap, MergeHeap.deleteMin heap

    let rec make_tree heap = 
        let first, heap = getem heap in
        if MergeHeap.is_empty heap then first
        else
            let second, heap = getem heap in
            let contents = 
                let count = (get_count first) + (get_count second) in
                Parser.Tree.Node ([first; second], (~-1, count))
            in
            make_tree (MergeHeap.insert contents heap)

    let rec create_codes acc trail tree =
        match tree with
        | Parser.Tree.Node ([l; r], _) ->
                let acc = create_codes acc (`S :: trail) l in
                create_codes acc (`K :: trail) r
        | Parser.Tree.Leaf (code, _) ->
                All_sets.IntegerMap.add code (List.rev trail) acc
        | _ -> assert false


    let assign_code_based_on_frequencies initial_frequencies functions env = 
        let count = 
            count_occurrences initial_frequencies functions env 
            All_sets.IntegerMap.empty
        in
        let heap = 
            All_sets.IntegerMap.fold (fun a b acc ->
            MergeHeap.insert (Parser.Tree.Leaf (a, b)) acc) 
            count MergeHeap.empty 
        in
        let tree = make_tree heap in
        let trails = create_codes All_sets.IntegerMap.empty [] tree in
        final_code := trails;
        fun x -> 
            let code =
                if x = identity then functions
                else int_of_string x
            in
            All_sets.IntegerMap.find code trails

    let add_identity identity_code env =
        let env = List.fold_left (compile Always) env (OCAMLSK let __identity a b = b) in
        match env with
        | [] -> assert false
        | e :: _ ->
                match All_sets.StringMap.find "__identity" e with
                | S_K x ->
                        decoder := ("__identity", identity_code, x) :: !decoder;
                        env
                | _ -> assert false

    let compile_decoder specs initial_frequencies =
        let env, functions = 
            List.fold_left (fun (env, counter) ->
                compile_decoder counter Always env) 
                ([All_sets.StringMap.empty], 0) specs
        in
        match env with
        | [env] -> 
                decoder := [];
                let final_code = 
                    assign_code_based_on_frequencies initial_frequencies functions env in
                let env = replace_labels final_code env in
                environment := add_identity functions [env]
        | _ -> assert false

    let decoder () = !decoder
    let final_code () = !final_code

    let compile spec =
        environment := List.fold_left (compile Always) !environment spec

    let clear () = environment := [All_sets.StringMap.empty]

    let evaluate def =
        create_sk AtTimes [] !environment def

    let get acc name =
        match acc with
        | S_K _ 
        | Argument _ -> assert false
        | SKModule environment -> 
                All_sets.StringMap.find name environment 

    let get name = 
        let path = Str.split (Str.regexp "\\.") name in
        match !environment with
        | [environment] ->
                (match List.fold_left get (SKModule environment) path with
                | S_K x -> x
                | _ -> failwith ("Illegal expression " ^ name))
        | _ -> assert false

    let tree_of_decoder () =
        let initial_decoder = decoder () 
        and final_code = final_code () in
        let decoder = List.map (fun (_, code, def) -> 
            All_sets.IntegerMap.find code final_code, def) initial_decoder
        in
        clear ();
        compile 
            (OCAMLSK 
                let m_true = [SK K] 
                let m_false = [SK (S K)]
                let m_and y x = if x then y else x                                  
                let m_or x y = if x then x else y                          
                let m_not x = if x then m_false else m_true
                let first = m_true                                                    
                let second = m_false
                let m_pair a b c =   c a b
                let pair = m_pair

                module Stream = struct
                    (* We don't use the continuation in this version *)
                    let to_bool continuation x = x [SK S] [SK K] [SK K] m_not m_true
                end
                module Decoder = struct

                    let rec generic_decoder decoder tree next =
                        let side =
                            if (Stream.to_bool next) then first tree
                            else second tree
                        in
                        let tip = second side in
                        if (Stream.to_bool (first side)) then 
                            (generic_decoder decoder tip)
                        else (tip (generic_decoder decoder decoder))
                                
                    let leaf x = pair [SK K] x

                    let node x y = pair (pair [SK S] x) (pair [SK S] y)
                end);
        let rec make_tree lst =
            match lst with
            | [(_, def)] -> 
                    compile (OCAMLSK let tmp = Decoder.leaf [def]);
                    get "tmp"
            | [] -> assert false
            | lst ->
                    let l, r = List.partition (function (`S :: t, _) -> true
                        | `K :: t, _ -> false | _ -> assert false) lst in
                    let reduce (x, y) = (List.tl x), y in
                    let l = make_tree (List.map reduce l)
                    and r = make_tree (List.map reduce r) in
                    compile (OCAMLSK let tmp = Decoder.node [r] [l]);
                    get "tmp"
        in
        make_tree decoder, initial_decoder, final_code

end

module Primitives = struct
    (** Primitive operations in an SK machine. In the following functions, any
    * one starting with an m_ is intended to be a macro that can be directly
    * embedded in an SK expression, while s_ functions are ocaml functions that  *)

    (** {1 Atomic Computation Types}
     * 
     * In order to perform all the computations in an S K machine we have to be able
     * to delay some operations and group arguments together. In order to do this we
     * must be able to define some basic types to allow complex things to happen,
     * like - for instance - if-then statements. The basic types that we are
     * interested on are booleans, pairs, lists, and integers. Using those basic
     * units we can compute almost anything.  
     * {2 Basic Types} 
     * {3 Booleans}
     *
     * A boolean can be represented using [K] for true and [S K] for false. This
     * has some nice properties. *)
    let m_true = (SK K)
    let m_false = (SK (S K))
    (** Basically observe that  
     * K A B -> A, ie true A B -> A, while S K A B -> K B (A * B) -> B ie. 
     * false A B -> B. This is exactly what we need. *)

    (** With this ready, we can now define the not, and, and or functions *)
    let s_and x y = (SK ([x] [y] [m_false]))
    let s_or x y = (SK ([x] [m_true] [y]))
    let s_not x = (SK ([x] (S K) K))

    (** And their corresponding macros for S_K *)
    let m_and = S_K.create (SK (x y [m_false])) (LS y x)
    let m_or = S_K.create (SK (x [m_true] y)) (LS x y)
    let m_not = S_K.create (SK (x (S K) K)) (LS x)
     
    (** We want to be able to
     * have a representation that is more concise if necessary. In other words,
     * observe that the representation of false is slightly longer than that of
     * true, and for compression purposes it would be desirable not to set a special
     * burden in the false case. For this we will also define K and S and true and
     * false using the special convertion functions that follow. *)
    let s_to_bool x = (SK ([x] S K K [m_not] [m_true])) 
    let m_to_bool = S_K.create (SK (x S K K [m_not] [m_true])) (LS x)

    (** {2 Pair and List Representation}
    *
    * If we have boolean representations, we can now pair representations by using
    * the true and false to get the first and last elements of the pair. *)
    let m_pair = S_K.create (SK (c a b)) (LS a b c)
    let m_first = m_true
    let m_second = m_false

    (** In the same way we can define a list as just a sequence of pairs. *)
    let m_head = m_first
    let m_tail = m_second

    (** {2 Recursion} 
     *
     * The classic fixedpoint theorem. *)
    let fixedpoint = (SK (S S K (S (K (S S(S(S S K))))K)))

    (** We create a convenience function to define recursive functions *)
    let rec_create name f ls = 
        let tmp = S_K.create f (name :: ls) in
        S_K.create (SK ([fixedpoint] [tmp])) []

    (** {2 Integer Representation}
    *
    * An integer has multiple possible representations. We will do the simplest and
    * then define special functions for more compact ones.
    *
    * {3 Church Integers} 
    * This is the simplest of all, consisting of a row of [i] [K]'s closed by an [S]
    * to represent the integer [i]. The decoder for a church integer produces a list
    * of integers *)
    let m_church_zero = (SK ([m_pair] [m_false] K))
    let s_church_successor x = (SK ([m_pair] [m_true] [x]))
    let s_church_predecessor x = (SK ([x] [m_second]))
    let s_church_not_zero x = (SK ([x] [m_head]))
    let m_church_successor = S_K.create (SK ([m_pair] [m_true] x)) (LS x)
    let m_church_predecessor = S_K.create (SK (x [m_second])) (LS x)
    let m_church_not_zero = S_K.create (SK (x [m_head])) (LS x)

    (** {4 Conversion Functions}
     * Conversion functions between [OCaml] and [S K] church representations *)
    let rec of_int x =
        if x < 0 then failwith "Invalid argument" 
        else
            if x = 0 then m_church_zero
            else s_church_successor (of_int (x - 1))

    let rec to_int x =
        if (SK (S K)) = (S_K.eval (SK ([x] [m_first]))) then 0
        else 1 + (to_int (S_K.eval (SK ([x] [m_second]))))

    (** {4 Basic arithmetic operations} *)
    let m_church_add = 
        rec_create "add" (SK (([s_church_not_zero  (SK a)])
        (add ([s_church_predecessor (SK a)]) ([s_church_successor (SK b)]))
        b)) (LS a b)

    let m_church_substract =
        rec_create "substract" (SK (([s_church_not_zero (SK b)])
            (substract ([s_church_predecessor (SK a)]) ([s_church_predecessor (SK b)]))
            a)) (LS a b)

    let m_church_multiply =
        rec_create "multiply" (SK (([s_church_not_zero (SK b)])
            ([m_church_add] a (multiply a ([s_church_predecessor (SK b)])))
            [m_church_zero])) (LS a b)

    module type List = sig
        type ocaml_repr
        type ocaml_list
        val empty_list : S_K.primitives
        val is_not_empty : S_K.primitives -> S_K.primitives
        val head : S_K.primitives -> S_K.primitives
        val tail : S_K.primitives -> S_K.primitives
        val prepend : S_K.primitives -> S_K.primitives -> S_K.primitives
        val ml_zero_signal : S_K.primitives -> bool
        val alphabet : (S_K.primitives * ocaml_repr) list
        val to_ocaml : S_K.primitives -> ocaml_repr
        val of_ocaml : ocaml_repr -> S_K.primitives
        val to_list : S_K.primitives -> ocaml_list
        val of_list : ocaml_list -> S_K.primitives
    end

    module type Integer = sig
        type extras
        (* A module type to convert from a given type to ChurchIntegers which
        * are easier to deal with *)
        val of_int : extras -> int -> S_K.primitives
        val to_int : extras -> S_K.primitives -> int
        val to_church : extras -> S_K.primitives
    end

    module ChurchIntegers : Integer with type extras = unit = struct
        (* A representation of integers using Church's original list
        * representation with as many K as there can be *)
        type extras = unit 

        let zero = m_church_zero
        let not_zero x = m_church_not_zero
        let ml_zero_signal x = 
           (SK (S K)) = S_K.eval (not_zero x)

        let to_int () = to_int

        let of_int () = of_int

        let to_church () = SK I
    end

    module K_S : List with type ocaml_repr = [ `K | `SK] 
    with type ocaml_list = [`K | `SK] list = struct
        let u_Hd = m_head
        let u_Tl = m_tail
            type ocaml_repr = [ `K | `SK ]
            type ocaml_list = ocaml_repr list
            let empty_list = SK (Pair (K S) K)
            let is_not_empty x = SK ([x] [u_Hd] K K)
            let ml_zero_signal x = (SK (K S)) = S_K.eval (SK ([x] [u_Hd]))
            let head x = SK ([x] [u_Hd])
            let tail x = SK ([x] [u_Tl])
            let prepend a b = SK (Pair [a] [b])
            let alphabet = [`K, `K; (SK (S K)), `SK]
            let inv_alph = List.map (fun (a, b) -> (b, a)) alphabet
            let to_ocaml a = 
                let a = S_K.eval a in
                List.assoc a alphabet
            let of_ocaml a = List.assoc a inv_alph
            let to_list lst = 
                let rec to_list lst =
                    if not (ml_zero_signal lst) then 
                        let next = to_ocaml (S_K.eval (head lst)) in
                        next :: (to_list (S_K.eval (tail lst)))
                    else []
                in
                to_list (S_K.eval lst) 
            let rec of_list lst = 
                List.fold_right (fun x acc ->
                    prepend (of_ocaml x) acc) lst empty_list
            let complement x = SK ([x] (S K) K)
    end

    module Dna : List with type ocaml_repr = string with type ocaml_list =
        string list = struct
        (* A representation of DNA sequences *)
        type ocaml_repr = string
        type ocaml_list = string list
        let empty_list = S_K.eval (SK (Pair (Pair (K S)  K) K))
        let is_not_empty x = (SK ([x] [m_head] [m_head] K K))
        let ml_zero_signal x = 
           (SK (K S)) = S_K.eval (SK ([x] [m_head] [m_head]))
        let head x = SK ([x] [m_head])
        let tail x = SK ([x] [m_tail])
        let prepend x y = SK (Pair [x] [y])
        let alphabet = 
            List.map (fun (a, b) -> S_K.eval a, b)
            [(SK (Pair K K), "A"); (SK (Pair K (S K)), "C"); 
            (SK (Pair (S K) K), "G"); (SK (Pair (S K) (S K)), "T")]
        let inv_alph = List.map (fun (a, b) -> (b, a)) alphabet
        let to_ocaml a = 
            match S_K.eval (SK (Hd [a])), S_K.eval (SK (Tl [a])) with
            | `K, `K -> "A"
            | `K, `Node [`S; `K] -> "C"
            | `Node [`S; `K], `K -> "G"
            | `Node [`S; `K], `Node [`S; `K] -> "T"
            | _ -> raise Not_found

        let of_ocaml a = List.assoc a inv_alph
        let complement x =
            let negate n = S_K.eval (SK ([x] [n] (S K) K)) in
            let a = negate m_head
            and b = negate m_tail in
            S_K.eval (SK (Pair [a] [b]))
        let to_list lst = 
            let rec to_list lst =
                if not (ml_zero_signal (S_K.eval lst)) then 
                    let next = to_ocaml (S_K.eval (head lst)) in
                    next :: (to_list (tail lst))
                else []
            in
            to_list (S_K.eval lst) 
        let rec of_list lst = 
            List.fold_right (fun x acc ->
                prepend (of_ocaml x) acc) lst empty_list
    end
    module LogInt : Integer with type extras = unit = struct
        type extras = unit

        let of_int () int =
            let rec encode_in_log int acc =
                if int = 0 then acc
                else 
                    encode_in_log (int lsr 1)
                    (SK (Pair [if 1 = (1 land int) then (SK K) else (SK (S K))]
                    [acc]))
            in
            if int < 0 then failwith "Illegal argument"
            else encode_in_log int (SK EmptyList)

        let to_int () x =
            let lst = K_S.to_list x in
            List.fold_left (fun acc item ->
                (2 * acc) +
                (match item with
                | `K -> 1 
                | `SK -> 0)) 0 lst

        let to_church () = SK DecodeInt 

    end


    module FixedMax : Integer with type extras = int = struct

        type extras = int 

        let of_int max int =
            let rec encode_in_log max int acc =
                if max = 0 then acc
                else 
                    encode_in_log (max lsr 1) (int lsr 1) 
                    (SK (Pair [if 1 = (1 land int) then (SK K) else (SK (S K))]
            [acc]))
            in
            if int < 0 || max < 0 || max < int then failwith "Illegal argument"
            else encode_in_log max int (SK K)

        let to_church max =
            let rec aux_len max acc =
                if max = 0 then acc 
                else aux_len (max lsr 1) (acc + 1) 
            in
            if max < 0 then failwith "Illegal argument" 
            else
                let counter = ChurchIntegers.of_int () (aux_len max 0) in
                let tmp_rec = 
                    S_K.create
                    (SK ([s_church_not_zero (SK cnt)]
                        (R [s_church_predecessor (SK cnt)] 
                        ((int [m_head]) 
                            (Add [ChurchIntegers.of_int () 1] (Multiply
                            [ChurchIntegers.of_int () 2] acc)) 
                            (Multiply [ChurchIntegers.of_int () 2] acc))
                            (int [m_tail]))
                        acc))
                    (LS R cnt acc int)
                in
                S_K.create 
                (SK (Fixedpoint [tmp_rec] [counter] [m_church_zero])) 
                []

        let to_int max int = 
            ChurchIntegers.to_int ()
            (S_K.eval (SK ([to_church max] [int])))

    end

    module FixedRange : Integer with type extras = (int * int) =
        struct
            type min = int
            type max = int
            type extras = (min * max)

            let of_int (min, max) int =
                if min > max || int < min || int > max then 
                    failwith "Illegal argument"
                else
                    FixedMax.of_int (max - min) int

            let to_int (min, max) int =
                min + (FixedMax.to_int (max - min) int)

            let to_church (min, max) =
                S_K.create
                (SK (Add [ChurchIntegers.of_int () min] ([FixedMax.to_church (max -
                min)] pos))) (LS pos)
        end
end

module IndPM = struct
    (** An implementation of a standard machine for phylogenetic computations using
    * compression and an S K machine *)

    open Primitives
    (** {2 Mutually Independent Phylogenetic Operations On DNA Sequences} *)
    let m_decode_base = 
        S_K.create (SK (Callback ([m_pair] [s_to_bool (SK a)]
        [s_to_bool (SK b)])))
        (LS Callback a b)

    (** Test the previous function 
        * Expected results K, K, (S K), (S K), K, (S K) 
        S_K.eval (SK ([m_decode_base] I K K [m_head]));;
        S_K.eval (SK ([m_decode_base] I K K [m_tail]));;
        S_K.eval (SK ([m_decode_base] I S S [m_tail]));;
        S_K.eval (SK ([m_decode_base] I S S [m_head]));;
        S_K.eval (SK ([m_decode_base] I K S [m_head]));;
        S_K.eval (SK ([m_decode_base] I K S [m_tail]));;
    *)

    let m_decode_church_integer =
        let tmp = 
            S_K.create 
            (SK ([(s_to_bool (SK next))]
                (decoder [s_church_successor (SK acc)] Callback)
                (Callback acc)))
            (LS decoder acc Callback next)
        in
        S_K.create (SK ([fixedpoint] [tmp] [m_church_zero])) []

    let m_decode_binary_integer =
        let tmp =
            S_K.create 
            (SK ([s_church_not_zero (SK [s_church_predecessor (SK cnt)])]
                ([s_to_bool (SK next)] 
                    (decoder add ([m_church_successor] (add acc acc))
                     Callback [s_church_predecessor (SK cnt)])
                    (decoder add (add acc acc) Callback [s_church_predecessor (SK cnt)] 
                    ))
                (Callback ([s_to_bool (SK next)] ([m_church_successor] (add acc
                acc)) (add acc acc)))))
            (LS decoder add acc Callback cnt next)
        in
        S_K.create (SK ([fixedpoint] [tmp] [m_church_add] [m_church_zero])) []

    let m_decode_binary_integer =
        S_K.create
        (SK ( [m_decode_church_integer] ([m_decode_binary_integer] Callback)))
        (LS Callback)

    (* Testing the previous functions 
    to_int (S_K.eval (SK ([m_decode_binary_integer] K K K K K S K S K K S)));;
    *)

    let decode_function =
        let tmp =
            S_K.create
            (SK ([s_to_bool (SK next)]
                ([s_to_bool (SK (encoded [m_tail] [m_head] [m_head]))]
                    (decoder (encoded [m_tail] [m_head]) Callback)
                    (Callback (encoded [m_tail] [m_head] [m_tail])))
                ([s_to_bool (SK (encoded [m_tail] [m_tail] [m_head]))]
                    (decoder (encoded [m_tail] [m_tail]) Callback)
                    (Callback (encoded [m_tail] [m_tail] [m_tail])))))
            (LS decoder encoded Callback next)
        in
        S_K.create (SK ([fixedpoint] [tmp])) []

    (* Testing the decode_function
    let l1 = (SK ([m_pair] S a))
    let l2 = (SK ([m_pair] S b))
    let l3 = (SK ([m_pair] S c))
    let join a b = (SK ([m_pair] K ([m_pair] [a] [b])))
    let encoded = join (join l2 l3) l1;;
    (S_K.to_string (S_K.eval (SK ([decode_function] [encoded] I S))));;
    *)

    let m_prepender = (* A function to prepend a base to an accumulator *)
        S_K.create (SK (Callback ([m_pair] x acc))) (LS Callback acc x)

    (* Testing m_prepender 
    * This function should output 0 
    to_int (S_K.eval (SK ([m_prepender] [m_church_predecessor] [m_church_zero] K)));;
    * This function should output 2
    to_int (S_K.eval (SK ([m_prepender] [m_church_successor] [m_church_zero] K)));;
    *)

    let m_decode_sequence =
        let shufle_last_two = S_K.create (SK (a e f b d c)) (LS a e f b c d) in
        let tmp =
            S_K.create 
            (SK ([s_church_not_zero (s_church_predecessor (SK cnt))] 
                (m_decode_base
                    (m_prepender
                        ([shufle_last_two] decoder m_decode_base m_prepender
                        Callback 
                            [s_church_predecessor (SK cnt)]) acc))
                (m_decode_base 
                    (m_prepender Callback acc))))
            (LS decoder m_decode_base m_prepender Callback acc cnt)
        in
        let decoder =
            S_K.create (SK ([fixedpoint] [tmp] [m_decode_base] [m_prepender] 
            Callback ([m_pair] ([m_pair] (K S) K) K))) (LS Callback)
        in
        S_K.create (SK ( [m_decode_binary_integer] ([decoder] (Callback stack
        results accumulator))))
        (LS Callback stack results accumulator editedsequence)
    (* An example of the use of this function is available in the module
    * interface. *)

    let merger = S_K.create (SK (a (b c))) (LS a b c)

    let m_insert_delta : S_K.primitives = 
        let tmp =
            S_K.create
            (SK ([s_church_not_zero (s_church_predecessor (SK pos))]
                (insert 
                    Callback
                    ([m_pair] (seq [m_head]) acc)
                    (seq [m_tail])
                    [s_church_predecessor (SK pos)]
                    base)
                (Callback ([m_pair] base acc) seq)))
            (LS insert Callback acc seq pos base)
        in
        let insert_f = S_K.create (SK ([fixedpoint] [tmp])) [] in
        S_K.create (SK ([m_decode_binary_integer] ([merger] [m_decode_base]
            ([insert_f] (Callback stack results) acc seq)))) (LS Callback stack
            results acc seq)

    (**  Testing m_insert_delta
    let fourth = S_K.create (SK d) (LS a b c d);;
    let seq = S_K.eval (SK ([m_decode_sequence] [fourth] stack results
        accumulator editedsequence K K S K K K K K K S K));;
    let empty_sequence = (SK ([Primitives.m_pair] ([Primitives.m_pair] (K S) K) K))
    let res = S_K.eval (SK ([m_insert_delta] Callback stack results
        [empty_sequence] [seq] K K S K S K S));; 
    let tot = S_K.eval (SK ([res] [m_head] [m_tail]));;
    let tot = S_K.eval (SK ([res] [m_tail] [m_tail] [m_head] [m_tail]));;
    *)

    let m_delete_delta =
        let tmp =
            S_K.create
            (SK ([s_church_not_zero (s_church_predecessor (SK pos))]
                (delete 
                    Callback
                    ([m_pair] (seq [m_head]) acc)
                    (seq [m_tail])
                    [s_church_predecessor (SK pos)])
                (Callback acc (seq [m_tail]))))
            (LS delete Callback acc seq pos)
        in
        let delete_f = S_K.create (SK ([fixedpoint] [tmp])) [] in
        S_K.create (SK ([m_decode_binary_integer] 
        ([delete_f] (Callback stack results) acc seq))) (LS Callback stack
        results acc seq)

    (**  Testing m_delete_delta
    let res = S_K.eval (SK ([m_delete_delta] K ([m_pair] K K) B K S K));; 
    let tot = S_K.eval (SK ([res] [m_head] [m_tail]));;
    *)

    let m_substi_delta = 
        let tmp =
            S_K.create
            (SK ([s_church_not_zero (s_church_predecessor (SK pos))]
                (substi 
                    Callback
                    ([m_pair] (seq [m_head]) acc)
                    (seq [m_tail])
                    [s_church_predecessor (SK pos)]
                    base)
                (Callback 
                    ([m_pair] base acc) 
                    (seq [m_tail]))))
            (LS substi Callback acc seq pos base)
        in
        let substi_f = S_K.create (SK ([fixedpoint] [tmp])) [] in
        S_K.create (SK ([m_decode_binary_integer] ([merger] [m_decode_base]
            ([substi_f] (Callback stack results) acc seq)))) (LS Callback stack
            results acc seq)


    (** Testing substitution
    let res = S_K.eval (SK ([m_substi_delta] K ([m_pair] A K) ([m_decode_sequence] I K K S K K K K K K S K) K K S K S K S));;
    let tot = S_K.eval (SK ([res] [m_head] [m_tail]));;
    let tot = S_K.eval (SK ([res] [m_head] [m_head]));;
    let tot = S_K.eval (SK ([res] [m_tail] [m_tail] [m_head] [m_tail]));;
    *)

    let rec make_encoder tree leafs =
        match tree with
        | `Node [a; b] -> 
                (SK (Pair K (Pair [make_encoder a leafs] [make_encoder b
                leafs])))
        | `Label x -> 
                (SK (Pair S [try (List.assoc x leafs) with Not_found ->
                    tree]))
        | `S -> (SK (Pair S S))
        | `K -> (SK (Pair S K))
        | `Lazy x -> make_encoder (Lazy.force_val x) leafs
        | `Debugger _ 
        | `Node _ -> failwith "Illegal encoder: each node must have two leafs
    and no debugging info"

    let m_decode_function =
        let tmp = 
            S_K.create 
            (SK (([s_to_bool (SK (encoder [m_tail] [s_to_bool (SK next)]
            [m_head]))])
                ( decoder original_encoder 
                    (encoder [m_tail] [s_to_bool (SK next)]) 
                    stack results acc seq)
                ((encoder [m_tail] [s_to_bool (SK next)] [m_tail]) 
                    ( decoder original_encoder original_encoder) 
                    stack results acc seq)))
            (LS decoder original_encoder encoder stack results acc seq next)
        in
        S_K.create (SK (  [fixedpoint] [tmp])) []

    (* Testing make_encoder and m_decode_function. The result of this function
        * call should be "Andres", as the [Callback] is selected by the [a]
        * function.
    let [a; b; c] = get_lst res;;
    let a = S_K.create (SK (a)) (LS a b c);;
    let enc = make_encoder (SK (A B)) [("A", a)];;
    let remove_first = S_K.create (SK e) (LS a b c d e);;
    S_K.to_string (S_K.eval (SK ([make_list] A B C D E)));;
    let enc = make_encoder (SK (A B)) ["A", remove_first];;
    let res = S_K.eval (SK ([m_decode_function] [enc] [enc] stack results acc
    seq K));;
    S_K.to_string res;; 
    *)

    let m_invert_full_sequence =
        let tmp =
            S_K.create
            (SK (((seq [m_head] [m_head]) K K)
                (R  ([m_pair] (seq [m_head]) acc) (seq [m_tail]))
                acc))
            (LS R acc seq)
        in
        S_K.create (SK ([fixedpoint] [tmp] ([m_pair] ([m_pair] (K S) K) K))) []

    (** Testing m_invert_full_sequence. The following should return the inverted
    * sequence, that is, the input sequence in the very same order. 
    let a = S_K.eval (SK ([m_decode_sequence] I K K S K K K S K K S K));;
    let b = S_K.eval (SK ([m_invert_full_sequence] I [a] ([m_pair] ([m_pair] (K
    S) K) K)));;
    S_K.eval (SK ([b] [m_head] [m_tail]));;
    S_K.eval (SK ([b] [m_head] [m_head]));;
    S_K.eval (SK ([b] [m_tail] [m_head] [m_tail]));;
    S_K.eval (SK ([b] [m_tail] [m_head] [m_head]));;
    S_K.eval (SK ([b] [m_tail] [m_tail] [m_head] [m_tail]));;
    S_K.eval (SK ([b] [m_tail] [m_tail] [m_head] [m_head]));;
    S_K.eval (SK ([b] [m_tail] [m_tail] [m_tail] [m_head] [m_head]));; 
    *)

    let m_copy_rest =
        let tmp =
            S_K.create
            (SK ((seq [m_head] [m_head] K K) 
            (R ([m_pair] (seq [m_head]) acc) (seq [m_tail]))
            acc))
            (LS R acc seq)
        in
        S_K.create (SK ([fixedpoint] [tmp])) []

    (** Testing m_copy_rest. The following should return the inverted input
    * sequence 
    let a = S_K.eval (SK ([m_decode_sequence] I K K S K K K S K K S K));;
    let b = S_K.eval (SK ([m_copy_rest] I ([m_pair] ([m_pair] (K
    S) K) K) [a] ));;
    S_K.eval (SK ([b] [m_head] [m_tail]));;
    S_K.eval (SK ([b] [m_head] [m_head]));;
    S_K.eval (SK ([b] [m_tail] [m_head] [m_tail]));;
    S_K.eval (SK ([b] [m_tail] [m_head] [m_head]));;
    S_K.eval (SK ([b] [m_tail] [m_tail] [m_head] [m_tail]));;
    S_K.eval (SK ([b] [m_tail] [m_tail] [m_head] [m_head]));;
    S_K.eval (SK ([b] [m_tail] [m_tail] [m_tail] [m_head] [m_head]));;  *)

    let m_end_branch =
        let empty_sequence = (SK ([m_pair] ([m_pair] (K S) K) K)) in
        let end_branch =
            S_K.create 
            (SK ([s_to_bool (SK next)]
            (* We have to add this item twice to the stack list and callback *)
                ( Callback ([m_pair] (copy_rest_n_invert acc seq) stack) results 
                    [empty_sequence] (copy_rest_n_invert acc seq) next2)
            (* We actually reached a leaf, so we don't add it to the stack but to
            * the results and continue *)
            (([s_to_bool (SK next2)])
                (* We have reached the end, time to output the result *)
                ( [m_pair] (copy_rest_n_invert acc seq) results)
                (* We have not reached the end let's keep going *)
                (  Callback (stack [m_tail]) ([m_pair] (copy_rest_n_invert acc seq)
                results) [empty_sequence] (stack [m_head])))))
            (LS copy_rest_n_invert Callback stack results acc seq next next2)
        in
        let copy_rest_n_invert = 
            S_K.create (SK ( [m_invert_full_sequence] ([m_copy_rest] acc seq)))
            (LS acc seq)
        in
        S_K.create (SK ( [end_branch] [copy_rest_n_invert])) 
        []

    (* Test copy_rest_n_invert:
    let res = S_K.eval (SK ([m_end_branch] [first] stack results
    [empty_sequence] [seq] S K));;
    let [a; b; c] = get_lst res;;
    let [x; y] = get_lst c;;
    let a = S_K.eval (SK ([res] [m_tail]));;
    let res = S_K.eval (SK ([copy_rest_n_invert] [empty_sequence] [seq]));;
    let res = S_K.eval (SK ([m_invert_full_sequence]  [seq]));;
        S_K.eval (SK ([a] [m_head] [m_tail]));;
        S_K.eval (SK ([a] [m_head] [m_head]));;
        S_K.eval (SK ([a] [m_tail] [m_head] [m_tail]));;
        S_K.eval (SK ([a] [m_tail] [m_head] [m_head]));;
        S_K.eval (SK ([a] [m_tail] [m_tail] [m_head] [m_tail]));;
        S_K.eval (SK ([a] [m_tail] [m_tail] [m_head] [m_head]));;
        S_K.eval (SK ([seq] [m_head] [m_tail]));;
        S_K.eval (SK ([seq] [m_head] [m_head]));;
        S_K.eval (SK ([seq] [m_tail] [m_head] [m_tail]));;
        S_K.eval (SK ([seq] [m_tail] [m_head] [m_head]));;
        S_K.eval (SK ([seq] [m_tail] [m_tail] [m_head] [m_tail]));;
        S_K.eval (SK ([seq] [m_tail] [m_tail] [m_head] [m_head]));;
        S_K.eval (SK ([res] [m_head] [m_tail]));;
        S_K.eval (SK ([res] [m_head] [m_head]));;
        S_K.eval (SK ([res] [m_tail] [m_head] [m_tail]));;
        S_K.eval (SK ([res] [m_tail] [m_head] [m_head]));;
        S_K.eval (SK ([res] [m_tail] [m_tail] [m_head] [m_tail]));;
        S_K.eval (SK ([res] [m_tail] [m_tail] [m_head] [m_head]));;
        S_K.eval (SK ([res] [m_tail] [m_tail] [m_tail] [m_head]));; 
    *)

end
