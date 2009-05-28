type definition = S_K.primitives Compiler.kolmo_function list

module Basic = struct
    let booleans = (OCAMLSK 
        let m_true = [SK K] 
        let m_false = [SK (S K)])

    let logic = (OCAMLSK
        let m_and y x = if x then y else x                                  
        let m_or x y = if x then x else y                          
        let m_not x = if x then m_false else m_true)

    let tuples = (OCAMLSK
        let first = m_true                                                    
        let second = m_false
        let m_pair a b c =   c a b
        let pair = m_pair)

    let lists = (OCAMLSK
        module List = struct
            let head = first
            let tail = second
            let prepend = pair
        end)

    let church_integers = (OCAMLSK
        module Church = struct
            let zero = pair m_false [SK K]
            let successor x = pair m_true x                                   
            let predecessor x = x second                                
            let not_zero x = x first
            let rec equal a b =
                if (m_and (not_zero a) (not_zero b)) then
                    (equal (predecessor a) (predecessor b))
                else
                    (if (m_not (m_or (not_zero a) (not_zero b))) then
                        m_true
                    else m_false)
            let rec gt a b =
                if (m_and (not_zero a) (not_zero b)) then
                    (gt (predecessor a) (predecessor b))
                else 
                    (m_and (not_zero a) (m_not (not_zero b)))
            let rec lt a b =
                if (m_and (not_zero a) (not_zero b)) then
                    (lt (predecessor a) (predecessor b))
                else (m_and (m_not (not_zero a)) (not_zero a))
            let rec add x y = 
                if (not_zero x) then (add (predecessor x) (successor y))
                else y
            let rec substract x y =
                if (not_zero y) then (substract (predecessor x) (predecessor y))
                else x
            let rec multiply x y =
                if (not_zero y) then 
                    (add x (multiply x (predecessor y)))
                else 0

            let rec log2 acc cur x =
                if (m_or (gt cur x) (equal cur x)) then acc
                else (log2 (successor acc) (add cur cur) x)

            let log2 x = log2 0 1 x

        end)

    let stream = (OCAMLSK
        module Stream = struct
            let to_bool x = x [SK S] [SK K] [SK K] m_not m_true

            (* We define a vanilla function that decodes the church integer in
            * the following stream *)
            let rec decode_church_integer acc next = 
                if (to_bool next) then 
                    (decode_church_integer (Church.successor acc))
                else acc
        end)
    let decoder = (OCAMLSK
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

            let rec integer_decoder tree continuation next =
                let side =
                    if (Stream.to_bool next) then first tree
                    else second tree
                in
                let tip = second side in
                if (Stream.to_bool (first side)) then 
                    (integer_decoder continuation tip)
                else (continuation tip)
                        
            let leaf x = pair [SK K] x

            let node x y = pair (pair [SK S] x) (pair [SK S] y)
        end)



    let integer_decoder = (OCAMLSK
        module IntegerDecoder = struct
            (* A function to compute log2 x *)
            let church_stream continuation = 
                let rec _aux_church_stream continuation acc next = 
                    if (Stream.to_bool next) then
                        (_aux_church_stream continuation (Church.successor acc))
                    else continuation acc
                in
                _aux_church_stream continuation 0


            let uniform next =
                let rec _aux_uniform acc bits next =
                    let nacc = 
                        Church.add (Church.add acc acc) 
                        (if (Stream.to_bool next) then 1 else 0)
                    in
                    if (Church.equal bits 1) then nacc
                    else (_aux_uniform nacc (Church.predecessor bits))
                in
                church_stream (_aux_uniform 0) next

            let uniform_max = 
                let rec _aux_uniform_max acc max next =
                    let nacc = 
                        Church.add (Church.add acc acc) 
                        (if (Stream.to_bool next) then 1 else 0) 
                    in
                    if (Church.equal max 1) then nacc
                    else (_aux_uniform_max nacc (Church.predecessor max))
                in
                _aux_uniform_max 0

            let rec uniform_min_max min max next =
                Church.add (uniform_max (Church.substract max min)) min

            let normal mean max threshold =
                let rec pectinated acc mean high low threshold =
                    if Church.lt high threshold then acc
                    else 
                        (Decoder.node (Decoder.leaf low) (Decoder.node
                        (Decoder.leaf high) threshold))
                in
                let low = Church.substract mean (Church.substract max mean) in
                let acc = Decoder.node (Decoder.leaf low) (Decoder.leaf mean) in
                let res = 
                    pectinated acc (Church.predecessor max) 
                    (Church.successor low) threshold
                in
                let rec merger left right a b =
                    if Stack.is_empty left then
                        if Stack.is_empty right then
                            if Stack.is_empty a then
                                (Decoder.integer_decoder (Stack.pop b))
                            else 
                                if Stack.is_empty b then
                                    (Decoder.integer_decoder (Stack.pop a))
                                else (merger a b Stack.empty Stack.empty)
                        else 
                            (merger (Stack.push (Stack.pop right) b) a
                            Stack.empty Stack.empty)
                    else 
                        if Stack.is_empty right then
                            (merger (Stack.push (Stack.pop left) b) a
                            Stack.empty Stack.empty)
                        else 
                            (merger (Stack.rest left) (Stack.rest right) 
                            (Stack.push (Decoder.node (Stack.pop left)
                            (Stack.pop right)) b) a)
                in
                let rec startup left right leftcnt rightcnt leftf rightf =
                    if Church.equal leftcnt mean then 
                        (merger left right Stack.empty Stack.empty)
                    else 
                        (startup right 
                        (Stack.push (Decoder.leaf leftcnt) left)
                        rightcnt (leftf leftcnt) rightf leftf)
                in
                let left = Stack.push (Decoder.leaf mean) Stack.empty in
                let right = Stack.push acc Stack.empty in
                let leftcnt = 
                    Church.substract mean (Church.substract threshold mean) 
                in
                startup left right leftcnt threshold Church.successor
                Church.predecessor

        end)

    let stack = (OCAMLSK
        module Stack = struct
            let empty = pair [SK K] [SK K]
            let push x stack  = pair [SK S] (pair x stack)
            let pop stack = first (second stack)
            let rest stack = second stack
            let is_empty stack = Stream.to_bool (first stack)
        end)


    let dna = (OCAMLSK
        module Dna = struct
            let empty = pair [SK K] [SK K]
            let prepend x stack  = pair [SK S] (pair x stack)
            let rest stack = second stack
            let is_empty stack = Stream.to_bool (first stack)
            let first stack = first (second stack)
            let rec insert sequence base position =
                if Church.not_zero position then 
                    (prepend (first sequence) (insert (rest sequence) base
                    (Church.predecessor position)))
                else (prepend base sequence)
        end)


    let huffman_decoder = (OCAMLSK
        (* We will define the structure of a generic Huffman like decoder
        * for streams of S and K *)
        let leaf x = pair [SK K] x
        let is_leaf x = first x
        let is_internal x = m_not (List.is_leaf x)

        let rec basic_decoder tree next =
            if (List.is_leaf tree) then
                (let a = second tree in
                (a next))
            else
                let b = Stream.to_bool next in
                (basic_decoder (b (second tree)))

        let rec basic_decoder_cb tree acc next = 
            if (is_leaf tree) then
                let a = second tree in
                a acc next
            else 
                let a = Stream.to_bool next in
                basic_decoder_cb (a (second tree)) acc)

    let load compiler = 
        let compile a b = Compiler.compile b a in
        List.fold_left compile compiler
        [booleans; logic; tuples; lists; church_integers; stream;
        stack; decoder; integer_decoder; dna; huffman_decoder]
end 

module Continuation = struct

    let church_integers = (OCAMLSK
        module ChurchC = struct
            let successor continuation x =
                continuation (Church.successor x)
            let predecessor continuation x =
                continuation (Church.predecessor x)
            let add continuation x y =
                continuation (Church.add x y)
            let substract continuation x y =
                continuation (Church.substract x y)
            let multiply continuation x y =
                continuation (Church.multiply x y)
            let log2 continuation x =
                continuation (Church.log2 x)
        end)

    let integer_decoder = (OCAMLSK
        module IntegerDecoderC = struct
            let rec _aux_uniform continuation acc bits next =
                let nacc = 
                    Church.add (Church.add acc acc) 
                    (if (Stream.to_bool next) then 1 else 0) 
                in
                if (Church.gt bits 1) then
                    (_aux_uniform continuation nacc (Church.predecessor bits))
                else (continuation nacc)

            let uniform continuation =
                IntegerDecoder.church_stream (_aux_uniform continuation 0)

            let rec _aux_uniform_max continuation acc max next =
                let nacc = 
                    Church.add (Church.add acc acc) 
                    (if (Stream.to_bool next) then 1 else 0) 
                in
                if (Church.equal max 1) then (continuation nacc)
                else (_aux_uniform_max continuation nacc (Church.predecessor max))

            let uniform_max continuation = _aux_uniform_max continuation 0

            let rec uniform_min_max continuation min max next =
                uniform_max (ChurchC.add continuation min) 
                (Church.substract max min)
        end)


    let stream = (OCAMLSK
        module StreamC = struct
            let rec _aux_decode_church_integer continuation acc next =
                if (Stream.to_bool next) then
                    _aux_decode_church_integer continuation
                    (Church.successor acc)
                else (continuation acc)
            let decode_church_integer continuation = 
                _aux_decode_church_integer continuation 0
        end)

    let dna = (OCAMLSK
        module DnaC = struct
            (* All of the functions that we now apply must take a callback
            * function that will be applied over the accumulator. *)
            let decode_base callback a b =
                callback (pair (Stream.to_bool a) (Stream.to_bool b))
        end)




    let tree = (OCAMLSK
        module Tree = struct

            let what_to_do_in_node process_branch stack results seq next =
                if (Stream.to_bool next) then
                    let next_call = Stack.pop stack in
                    (next_call (Stack.rest stack) (Stack.push seq results))
                else (process_branch stack results seq)

        end)


    end
