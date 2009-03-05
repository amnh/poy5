let pervasives = 
    OCAMLSK 
    let m_true = [SK K] 
    let m_false = [SK (S K)]
    let m_and y x = if x then y else x                                  
    let m_or x y = if x then x else y                          
    let m_not x = if x then m_false else m_true
    let first = m_true                                                    
    let second = m_false
    let m_pair a b c =   c a b
    let pair = m_pair

    module List = struct
        let head = first
        let tail = second
        let prepend = pair
    end

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

    end

    module Stream = struct
        let to_bool x = x [SK S] [SK K] [SK K] m_not m_true

        (* We define a vanilla function that decodes the church integer in
        * the following stream *)
        let rec decode_church_integer acc next = 
            if (to_bool next) then 
                (decode_church_integer (Church.successor acc))
            else acc
    end

    module Stack = struct
        let empty = pair [SK K] [SK K]
        let push x stack  = pair [SK S] (pair x stack)
        let pop stack = first (second stack)
        let rest stack = second stack
        let is_empty stack = Stream.to_bool (first stack)
    end

    module Dna = struct
        let empty = pair [SK K] [SK K]
        let prepend x stack  = pair [SK S] (pair x stack)
        let rest stack = second stack
        let is_empty stack = Stream.to_bool (first stack)
        let first stack = first (second stack)
    end

    let rec insert sequence base position =
        if Church.not_zero position then 
            (Dna.prepend (Dna.first sequence) (insert (Dna.rest sequence) base
            (Church.predecessor position)))
        else (Dna.prepend base sequence)

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
            basic_decoder_cb (a (second tree)) acc

    module IntegerDecoder = struct
        (* A function to compute log2 x *)
        let rec _aux_church_stream continuation acc next = 
            if (Stream.to_bool next) then
                (_aux_church_stream continuation (Church.successor acc))
            else continuation acc

        let church_stream continuation = _aux_church_stream continuation 0

        let rec _aux_uniform acc bits next =
            let nacc = 
                Church.add (Church.add acc acc) 
                (if (Stream.to_bool next) then 1 else 0)
            in
            if (Church.equal bits 1) then nacc
            else (_aux_uniform nacc (Church.predecessor bits))

        let uniform next =
            church_stream (_aux_uniform 0) next

        let rec _aux_uniform_max acc max next =
            let nacc = 
                Church.add (Church.add acc acc) (if (Stream.to_bool next) then 1 else 0) 
            in
            if (Church.equal max 1) then nacc
            else (_aux_uniform_max nacc (Church.predecessor max))

        let uniform_max = _aux_uniform_max 0

        let rec uniform_min_max min max next =
            Church.add (uniform_max (Church.substract max min)) min
    end

    module Continuation = struct
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
        end

        module StreamC = struct
            let rec _aux_decode_church_integer continuation acc next =
                if (Stream.to_bool next) then
                    _aux_decode_church_integer continuation
                    (Church.successor acc)
                else (continuation acc)
            let decode_church_integer continuation = 
                _aux_decode_church_integer continuation 0

        end

        module Dna = struct
            (* All of the functions that we now apply must take a callback
            * function that will be applied over the accumulator. *)
            let decode_base callback a b =
                callback (pair (Stream.to_bool a) (Stream.to_bool b))
        end

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

        end



        module Tree = struct

            let what_to_do_in_node process_branch stack results seq next =
                if (Stream.to_bool next) then
                    let next_call = Stack.pop stack in
                    (next_call (Stack.rest stack) (Stack.push seq results))
                else (process_branch stack results seq)

        end


    end

    module PhyloTree = struct

        module Stream = struct
            (* We don't use the continuation in this version *)
            let to_bool continuation x = x [SK S] [SK K] [SK K] m_not m_true
        end

        module Stack = struct
            let empty = pair [SK K] [SK K]
            let push continuation x stack  = pair [SK S] (pair x stack)
            let pop continuation stack = first (second stack)
            let rest continuation stack = second stack
            let is_empty continuation stack = Stream.to_bool (first stack)
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
        end

        module Dna = struct
            let empty = Stack.empty
            let prepend = Stack.push 
            let prepend_rev continuation seq x = 
                continuation (Stack.push continuation x seq)
            let head = Stack.pop 
            let tail = Stack.rest 
            let is_empty = Stack.is_empty 

            let rec _aux_length acc seq = 
                if (Stack.is_empty seq) then acc
                else _aux_length (Church.successor acc) seq

            let length seq = _aux_length 0 seq

            let decode_base continuation a b =
                continuation (pair (Stream.to_bool a) (Stream.to_bool b))

            let _aux_append_base continuation seq = 
                decode_base (prepend_rev continuation seq)

            let rec _aux_decode_sequence continuation cnt seq =
                if (Church.not_zero cnt) then
                    (_aux_append_base (_aux_decode_sequence continuation
                    (Church.predecessor cnt)))
                else (continuation seq)

            let _aux_decode_sequence continuation cnt =
                _aux_decode_sequence continuation cnt empty

            let decode_sequence continuation =
                Continuation.IntegerDecoderC.uniform 
                (_aux_decode_sequence continuation)
        end

        module Independent = struct

            let rec _aux_do_insert sequence offset base =
                if (Church.not_zero offset) then
                    (Dna.prepend (Dna.head sequence) 
                        (_aux_do_insert (Dna.tail sequence) (Church.predecessor
                        offset) base))
                else (Dna.prepend base sequence)

            let _aux_prepend_to_stack continuation stack sequence offset base =
                continuation 
                (Stack.push (_aux_do_insert sequence offset base) stack)

            let _aux_prepend_base continuation stack sequence offset = 
                Dna.decode_base
                (_aux_prepend_to_stack continuation stack sequence offset)

            let _aux_prepend_offset continuation stack sequence =
                Continuation.IntegerDecoder.uniform 
                (_aux_prepend_base continuation stack sequence) 

            let insert continuation results stack =
                _aux_prepend_offset (continuation results) (Stack.rest stack)
                (Stack.pop stack)

            let rec _aux_do_substitute sequence offset base =
                if (Church.not_zero offset) then
                    (Dna.prepend (Dna.head sequence) 
                        (_aux_do_insert (Dna.tail sequence) (Church.predecessor
                        offset) base))
                else (Dna.prepend base (Dna.tail sequence))

            let _aux_prepend_to_results continuation stack sequence offset base =
                continuation 
                (Stack.push (_aux_do_substitute sequence offset base)
                stack)

            let _aux_prepend_base continuation stack sequence offset = 
                Dna.decode_base
                (_aux_prepend_to_results continuation stack sequence offset)

            let _aux_prepend_offset continuation stack sequence =
                Continuation.IntegerDecoder.uniform 
                (_aux_prepend_base continuation stack sequence) 

            let substitute continuation results stack =
                _aux_prepend_offset (continuation results) (Stack.rest stack)
                (Stack.pop stack)

            let rec _aux_do_delete sequence offset =
                if (Church.not_zero offset) then
                    (Dna.prepend (Dna.head sequence) 
                        (_aux_do_insert (Dna.tail sequence) (Church.predecessor
                        offset)))
                else (Dna.tail sequence)

            let _aux_prepend_to_results continuation stack sequence offset =
                continuation 
                (Stack.push (_aux_do_delete sequence offset) stack)

            let _aux_prepend_offset continuation stack sequence =
                Continuation.IntegerDecoder.uniform 
                (_aux_prepend_to_results continuation stack sequence) 

            let delete continuation results stack =
                _aux_prepend_offset (continuation results) (Stack.rest stack)
                (Stack.pop stack)

        end

        module Branch = struct
            let leaf continuation results stack =
                continuation (Stack.push (Stack.pop stack) results) 
                (Stack.rest stack)

            let interior continuation results stack =
                continuation results (Stack.push (Stack.pop stack) stack)

            let ended continuation results stack = results

        end

        module Tree = struct
            let _aux_root continuation results stack seq =
                continuation results (Stack.push seq stack)

            let root continuation results stack = 
                Dna.decode_sequence (_aux_root continuation results stack)
        end


    end
