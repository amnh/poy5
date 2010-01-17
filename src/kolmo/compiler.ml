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

type compiler = {
    environment : environment_item All_sets.StringMap.t list;
    decoder :
        (All_sets.StringMap.key * All_sets.IntegerMap.key * S_K.primitives) list;
    final_code : S_K.primitives list All_sets.IntegerMap.t;
}

let compiler = { environment = [All_sets.StringMap.empty]; 
    decoder = []; final_code = All_sets.IntegerMap.empty }


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

let add n code cnt =
    if All_sets.IntegerMap.mem code cnt then 
        All_sets.IntegerMap.add code 
        ((All_sets.IntegerMap.find code cnt) + n) cnt
    else All_sets.IntegerMap.add code n cnt

let rec count_occurrences ?(prefix="") initial_frequencies identity_code env 
res =
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
        | DSKModule m -> 
                let prefix = 
                    if prefix ^ name = "" then assert false
                    else prefix ^ name ^ "." 
                in
                count_occurrences ~prefix initial_frequencies 
                identity_code m res
        | NS_K prim -> internal_counter res prim 
        | DS_K (prim, code) ->
                    let initial = 
                        let name = prefix ^ name in
                        try List.assoc name initial_frequencies with
                        | Not_found -> 1
                    in
                    internal_counter (add initial code res) prim) env res

let rec replace_labels ?(prefix="") replacer env compiler =
    let rec internal_replacer (x : S_K.primitives) : S_K.primitives list = 
        match x with
        | `S | `K | `Debugger _ | `Lazy _ -> [x]
        | `Label x -> replacer x
        | `Node lst ->
                let lst = List.map internal_replacer lst in
                let lst = List.flatten lst in
                [`Node lst]
    in
    All_sets.StringMap.fold (fun name x ((res, compiler) as acc) ->
        match x with
        | DArgument _ -> acc
        | DSKModule m -> 
                let prefix = 
                    if prefix = "" then name
                    else prefix ^ "." ^ name 
                in
                let labels, compiler = replace_labels ~prefix replacer m compiler in
                All_sets.StringMap.add name (SKModule labels) res, compiler
        | DS_K (prim, c) ->
                let f = 
                    match internal_replacer prim with
                    | [x] -> x
                    | x -> (`Node x)
                in
                let compiler = 
                    let name = 
                        if prefix = "" then name
                        else prefix ^ "." ^ name 
                    in
                    { compiler with decoder = (name, c, f) :: compiler.decoder }
                in
                All_sets.StringMap.add name (S_K f) res, compiler
        | NS_K prim -> 
                All_sets.StringMap.add
                name
                (match internal_replacer prim with
                | [x] -> S_K x
                | x -> S_K (`Node x)) res, compiler)
    env (All_sets.StringMap.empty, compiler)

let get_count a =
    match a with
    | Tree.Parse.Leafp (_, a)
    | Tree.Parse.Nodep (_, (_, a)) -> a

module OrderedPair = struct
    type t = (int * int) Tree.Parse.t
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
            Tree.Parse.Nodep ([first; second], (~-1, count))
        in
        make_tree (MergeHeap.insert contents heap)

let rec create_codes acc trail tree =
    match tree with
    | Tree.Parse.Nodep ([l; r], _) ->
            let acc = create_codes acc (`S :: trail) l in
            create_codes acc (`K :: trail) r
    | Tree.Parse.Leafp (code, _) ->
            All_sets.IntegerMap.add code (List.rev trail) acc
    | _ -> assert false


let assign_code_based_on_frequencies initial_frequencies functions env compiler = 
    let count = 
        count_occurrences initial_frequencies functions env 
        All_sets.IntegerMap.empty
    in
    let heap = 
        All_sets.IntegerMap.fold (fun a b acc ->
        MergeHeap.insert (Tree.Parse.Leafp (a, b)) acc) 
        count MergeHeap.empty 
    in
    let tree = make_tree heap in
    let trails = create_codes All_sets.IntegerMap.empty [] tree in
    { compiler with final_code = trails },
    fun x -> 
        let code =
            if x = identity then functions
            else int_of_string x
        in
        All_sets.IntegerMap.find code trails

let add_identity identity_code compiler =
    let environment = 
        List.fold_left (compile Always) compiler.environment 
        (OCAMLSK let __identity a b = b) in
    match environment with
    | [] -> assert false
    | e :: _ ->
            match All_sets.StringMap.find "__identity" e with
            | S_K x ->
                    { compiler with 
                        decoder = ("__identity", identity_code, x) :: 
                            compiler.decoder } 
            | _ -> assert false

let compile_decoder specs initial_frequencies compiler =
    let env, functions = 
        List.fold_left (fun (env, counter) ->
            compile_decoder counter Always env) 
            ([All_sets.StringMap.empty], 0) specs
    in
    match env with
    | [env] -> 
            let compiler = { compiler with decoder = [] } in
            let compiler, final_code = 
                assign_code_based_on_frequencies initial_frequencies 
                functions env compiler
            in
            let env, compiler = replace_labels final_code env compiler in
            let compiler = { compiler with environment = [env] } in
            add_identity functions compiler
    | _ -> assert false

let decoder compiler = compiler.decoder
let final_code compiler = compiler.final_code

let compile spec compiler =
    { compiler with environment = 
        List.fold_left (compile Always) compiler.environment spec }
 
let clear compiler = 
    { compiler with environment = [All_sets.StringMap.empty] }

let evaluate def compiler =
    create_sk AtTimes [] compiler.environment def

let get acc name =
    match acc with
    | S_K _ 
    | Argument _ -> assert false
    | SKModule environment -> 
            All_sets.StringMap.find name environment 

let get compiler name = 
    let path = Str.split (Str.regexp "\\.") name in
    match compiler.environment with
    | [environment] ->
            (match List.fold_left get (SKModule environment) path with
            | S_K x -> x
            | _ -> failwith ("Illegal expression " ^ name))
    | _ -> assert false

let tree_of_decoder compiler =
    let initial_decoder = decoder compiler
    and final_code = final_code compiler in
    let decoder = List.map (fun (_, code, def) -> 
        All_sets.IntegerMap.find code final_code, def) initial_decoder
    in
    let compiler = 
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
                    let to_bool continuation x = 
                        x [SK S] [SK K] [SK K] m_not m_true
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
                end) (clear compiler)
    in
    let rec make_tree compiler lst =
        match lst with
        | [(_, def)] -> 
                let compiler = 
                    compile (OCAMLSK let tmp = Decoder.leaf [def]) compiler
                in
                get compiler "tmp", compiler
        | [] -> assert false
        | lst ->
                let l, r = List.partition (function (`S :: t, _) -> true
                    | `K :: t, _ -> false | _ -> assert false) lst in
                let reduce (x, y) = (List.tl x), y in
                let l, compiler = make_tree compiler (List.map reduce l) in
                let r, compiler = make_tree compiler (List.map reduce r) in
                let compiler = 
                    compile (OCAMLSK let tmp = Decoder.node [r] [l]) 
                    compiler 
                in
                get compiler "tmp", compiler
    in
    fst (make_tree compiler decoder), initial_decoder, final_code


let uniform_integer compiler decoder integer =
    let r = get compiler decoder in
    let rec prepend items acc = 
        if items = 0 then acc
        else prepend (items - 1) (`S :: acc)
    in
    let rec generate_list cnt acc integer = 
        if integer = 0 then 
            prepend cnt acc
        else 
            generate_list (cnt + 1)
            ((if 1 = (1 land integer) then `K else `S) :: acc)
            (integer lsr 1)
    in
    `Node (r ::
        (if integer = 0 then [`S; `S] 
            else generate_list 0 [] integer))

let (-->) a b = b a

let complexity compiler funct =
    funct --> get compiler --> S_K.s_encode --> List.length
