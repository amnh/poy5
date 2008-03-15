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

module S_K = struct
    type primitives = [ `S | `K | `Label of string | `Node of primitives list ]

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
        let tree = 
            let rec string_tree x =
                match x with
                | `Node lst ->
                        Parser.Tree.Node 
                        (List.map string_tree lst, "")
                | `S | `K | `Label _ as x -> 
                        let string = 
                            match x with
                            | `S -> "S"
                            | `K -> "K"
                            | `Label x -> x
                        in
                        Parser.Tree.Leaf string
            in
            string_tree tree
        in
        let str = ref "" in
        let my_printer x = str := !str ^ x in
        AsciiTree.draw_parenthesis false my_printer tree;
        !str

    let rec sk_define name tree = 
        Hashtbl.replace universe name (expand tree) 

    let rec simplify tree = 
        match tree with
        | `Node lst ->
                (match lst with
                | [`S; _]
                | [`K; _]
                | [`Label _ ]  | [ `S ] | [ `K ]
                | [(`S); _; _] -> lst
                | (`K) :: a :: b :: t -> 
                        let a = simplify a in
                        [`Node (a @ t)]
                | (`S) :: a :: b :: c :: t ->
                        let a = simplify a in
                        [`Node ((a @ [c] @ [`Node [b; c]] @ t))]
                | h :: t -> 
                        let h = simplify h in
                        [`Node (h @ t)]
                | [] -> raise (Illegal_Expression []))
        | `S | `K | `Label _ -> [tree]


    let rec reduce tree =
        match simplify tree with
        | [ntree] -> 
                if ntree = tree then tree 
                else reduce ntree
        | [] -> raise (Illegal_Expression [])
        | lst -> 
                let ntree = (`Node lst) in
                if ntree = tree then tree 
                else reduce tree

    let evaluate x = 
        x --> of_string --> reduce 

    let eval x = x --> expand --> reduce 

    let test lst = 
        reduce (`Node lst)

    let s_encode tree = 
        (* We define a function to encode an SK expression in a list of bits *)
        let rec to_bit_list ?(no_print=false) (tree : primitives) =
            match tree with
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
                if not (ml_zero_signal lst) then 
                    let next = to_ocaml (S_K.eval (head lst)) in
                    next :: (to_list (S_K.eval (tail lst)))
                else []
            in
            to_list (S_K.eval lst) 
        let rec of_list lst = 
            List.fold_right (fun x acc ->
                prepend (of_ocaml x) acc) lst empty_list
    end

    module Debugger (R : List) = struct
        let rec string_list (lst : S_K.primitives) =
            if R.ml_zero_signal lst then 
                (R.to_ocaml (S_K.eval (R.head lst))) :: 
                    string_list (S_K.eval (R.tail lst))
            else []
    end

    module type SE = sig 
        val insert : S_K.primitives
        val delete : S_K.primitives
        val substitute : S_K.primitives
        type insert_generator 
        type delete_generator 
        type substitute_generator
        val insert_generator : insert_generator
        val delete_generator : delete_generator
        val substitute_generator : substitute_generator
        (*
        val functions_provided : (string * string * S_K.primitives) list
        val machine_description : string
        *)
    end

    module AtomicSE (I : Integer) (S : List) : SE 
        with type insert_generator = int -> S.ocaml_repr -> S_K.primitives
        with type delete_generator = int -> S_K.primitives
        with type substitute_generator = int -> S.ocaml_repr -> S_K.primitives 
        = struct
        let until_zero_recursions lst pos if_zero if_not_zero =
            let simple_recursion =
                S_K.create (SK 
                    ([I.not_zero (SK pos)]
                        [if_not_zero] 
                        [if_zero])) 
                lst
            in
            S_K.expand (SK (Fixedpoint [simple_recursion]))
        (* These functions accept atomic operations *)
        let insert = 
            until_zero_recursions
            (LS R pos base seq)
            (SK pos)
            (S.prepend (SK base) (SK seq))
            (S.prepend (S.head (SK seq))
            (SK (R [I.predecessor (SK pos)] base [S.tail (SK seq)])))

        let delete = 
            until_zero_recursions
            (LS R pos seq)
            (SK pos)
            (S.tail (SK seq))
            (S.prepend 
                (S.head (SK seq)) 
                (SK (R [I.predecessor (SK pos)] [S.tail (SK seq)])))

        let substitute =
            until_zero_recursions
            (LS R pos base seq)
            (SK pos)
            (S.prepend (SK base) (S.tail (SK seq)))
            (S.prepend (S.head (SK seq))
                (SK (R [I.predecessor (SK pos)] base [S.tail (SK seq)])))

        type insert_generator = int -> S.ocaml_repr -> S_K.primitives
        type delete_generator = int -> S_K.primitives
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
        (* The following are atomic affine editions *)
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

        type insert_generator = int -> S_K.primitives 
        type delete_generator = int -> int -> S_K.primitives
        type substitute_generator = int -> S_K.primitives 

        let insert_generator pos = 
            SK ([insert] [I.of_int pos])

        let delete_generator pos len =
            SK ([delete] [I.of_int pos] [I.of_int len])

        let substitute_generator pos =
            SK ([substitute] [I.of_int pos])
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

        let generic f = S_K.create (SK ([f] (DecodeInt pos))) (LS pos)
        let insert = generic S.insert
        let delete = generic S.delete
        let substitute = generic S.substitute

        type insert_generator = int -> S_K.primitives 
        type delete_generator = int -> S_K.primitives
        type substitute_generator = int -> S_K.primitives 

        let insert_generator pos = 
            SK ([insert] [encode_in_log pos])

        let delete_generator pos =
            SK ([delete] [encode_in_log pos])

        let substitute_generator pos =
            SK ([substitute] [encode_in_log pos])
    end

    module AffineSE_LogInts (S : SE) : SE 
        with type insert_generator = int -> S_K.primitives
        with type delete_generator = int -> int -> S_K.primitives
        with type substitute_generator = int -> S_K.primitives 
        = struct
        let generic f = S_K.create (SK ([f] (DecodeInt pos))) (LS pos)
        let insert = generic S.insert
        let delete = 
            S_K.create (SK ([S.delete] (DecodeInt pos) (DecodeInt len)))
            (LS pos len)
        let substitute = generic S.substitute

        type insert_generator = int -> S_K.primitives 
        type delete_generator = int -> int -> S_K.primitives
        type substitute_generator = int -> S_K.primitives 

        let insert_generator pos = 
            SK ([insert] [encode_in_log pos])

        let delete_generator pos len =
            SK ([delete] [encode_in_log pos] 
            [encode_in_log len])

        let substitute_generator pos =
            SK ([substitute] [encode_in_log pos])
    end

    module type HO = sig
        include SE
        val translocate : S_K.primitives
        val invert : S_K.primitives
        val duplicate : S_K.primitives
        val tandem : S_K.primitives
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
    module Log = struct
        (* Now we convert the previous three functions to accept an encoded integer
        * in logarithmic space *)
        let convert_to_log_position x = 
            S_K.expand (S_K.create (SK ([x] (DecodeInt pos))) (LS pos))
        let insert = convert_to_log_position insert
        let delete = convert_to_log_position delete
        let substitute = convert_to_log_position substitute
    end

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
