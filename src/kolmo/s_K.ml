let (-->) a b = b a

type primitives = [ `S | `K | `Label of string | `Node of primitives list |
    `Debugger of string | `Lazy of primitives Lazy.t]

exception Illegal_Expression of string Tree.Parse.t list

let universe = Hashtbl.create 97 

let of_string string =
    let no_forest = function
        | [x] -> 
                let rec reverse tree = 
                    match tree with
                    | Tree.Parse.Nodep (s, "") ->
                            `Node (List.rev_map reverse s)
                    | Tree.Parse.Leafp "S" -> `S
                    | Tree.Parse.Leafp "K" -> `K
                    | Tree.Parse.Leafp x -> `Label x
                    | x -> raise (Illegal_Expression [x])
                in 
                reverse (Tree.Parse.strip_tree x)
        | x -> raise (Illegal_Expression (List.map Tree.Parse.strip_tree x))
    in
    string 
        --> Tree.Parse.of_string
        --> List.map no_forest
        --> List.hd


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
