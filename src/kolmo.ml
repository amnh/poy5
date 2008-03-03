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

    let until_zero_recursions lst pos if_zero if_not_zero =
        let simple_recursion =
            S_K.create (SK ([pos] [u_Hd] K K [if_not_zero] [if_zero])) 
            lst
        in
        S_K.expand (SK (Fixedpoint [simple_recursion]))

    (* These functions accept atomic operations *)
    let insert = 
        until_zero_recursions 
        (LS R pos base seq) 
        (SK pos) 
        (SK (Pair base seq))
        (SK (Pair (seq [u_Hd]) (R (pos [u_Tl]) base (seq [u_Tl]))))

    let delete = 
        until_zero_recursions
        (LS R pos seq)
        (SK pos)
        (SK (seq [u_Tl]))
        (SK (Pair (seq [u_Hd]) (R (pos [u_Tl]) (seq [u_Tl]))))

    let substitute =
        until_zero_recursions
        (LS R pos base seq)
        (SK pos)
        (SK (Pair base (seq [u_Tl])))
        (SK (Pair (seq [u_Hd]) (R (pos [u_Tl]) base (seq [u_Tl]))))

    (* Now we convert the previous three functions to accept an encoded integer
    * in logarithmic space *)
    let convert_to_log_position x = 
        S_K.expand (S_K.create (SK ([x] (DecodeInt pos))) (LS pos))

    let insert = convert_to_log_position insert
    let delete = convert_to_log_position delete
    let substitute = convert_to_log_position substitute

    (* The following are atomic affine editions *)
    let affine_insert = 
        (* Instead of passing a single base, we will pass a
        complete sequence to this function *)
        let tmp = 
            S_K.create 
            (SK ((pos [u_Hd] K K) (* If we are not in position 0 we continue *)
                (Pair (seq [u_Hd]) (R (pos [u_Tl]) prefix (seq [u_Tl])))
                (* We have reached the position, we better merge now *)
                ((prefix [u_Hd] K K) 
                (* If we have some base left in the sequence that we will insert
                * *)
                    (Pair (prefix [u_Tl] [u_Hd]) (R pos (prefix [u_Tl] [u_Tl])
                    seq))
                    (* Otherwise we have finished inserting *)
                    seq)))
            (LS R pos prefix seq)
        in
        S_K.expand (SK (Fixedpoint [tmp]))

    let affine_delete =
        let tmp =
            S_K.create
            (SK ((pos [u_Hd] K K) (* If we haven't reached the deletion pos *)
                (R (pos [u_Tl]) len seq)
                (* If we have reached the deletion pos *)
                ((len [u_Hd] K K) (* We still have something more to del *)
                    (R pos (len [u_Tl]) (seq [u_Tl]))
                    (* We are done with the deletion *)
                    seq)))
            (LS R pos len seq)
        in
        S_K.expand (SK (Fixedpoint [tmp]))

    let affine_substitution =
        let tmp =
            S_K.create
            (SK ((pos [u_Hd] K K) (* If we haven't reached the deletion position *)
                (R (pos [u_Tl]) repl seq)
                (* We have reached the deletion position *)
                ((repl [u_Hd] K K) (* Still some more to substitute *)
                    (Pair (repl [u_Tl] [u_Hd]) (R pos (repl [u_Tl] [u_Tl]) (seq
                    [u_Tl])))
                    (* Nothing else to substitute *)
                    seq)))
            (LS R pos repl seq)
        in
        S_K.expand (SK (Fixedpoint [tmp]))

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

    let adenine = SK (Pair K K)
    let citosine = SK (Pair K (S K))
    let guanine = SK (Pair (S K) K)
    let timine = SK (Pair (S K) (S K))

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
        prepend ((String.length a) - 1) (`Label "EmptyList")

    let generate_machine_from_alignment a b =
        (* We want to produce a machine that given the sequence [a] produces as
        * output the sequence [b]. *)
        assert ((String.length a) = (String.length b));
        let max = String.length a in
        let base = function
            | 'A' -> adenine
            | 'C' -> citosine
            | 'G' -> guanine
            | 'T' | 'U' -> timine
            | '_' -> `Label "Emptylist"
            | _ -> failwith "Illegal DNA sequence"
        in
        let rec compare_and_generate_edition delta position =
            if max = position then []
            else if a.[position] = b.[position] then 
                compare_and_generate_edition (delta + 1) (position + 1)
            else 
                let npos = position + 1 
                and deltap = encode_in_log delta 
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

end
