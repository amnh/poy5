let simple_indelsubstitution compiler =
    Kolmo.Compiler.compile_decoder (OCAMLSK
        let m_true = [SK K]
        let m_false = [SK (S K)]
        let m_and y x = if x then y else x
        let m_or x y = if x then x else y
        let m_not x = if x then m_false else m_true

        let first = m_true
        let second = m_false
        val pair a b c = c a b

        module Stream = struct
            (* We don't use the continuation in this version *)
            let to_bool continuation x = x [SK S] [SK K] [SK K] m_not m_true
        end

        module Stack = struct
            val push x stack = pair [SK S] (pair x stack)
            val empty = pair [SK K] [SK K]
            val pop stack = first (second stack)
            val rest stack = second stack
            val is_empty stack = Stream.to_bool (first stack)
        end
        module Church = struct
            val successor x = pair m_true x                                   
            val predecessor x = x second                                
            val not_zero x = x first
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

            let log2 x = 
                let rec log2 acc cur x =
                    if (m_or (gt cur x) (equal cur x)) then acc
                    else (log2 (successor acc) (add cur cur) x)
                in
                log2 0 1 x

        end

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

        end

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

        module IntegerDecoderC = struct

            val uniform continuation =
                let rec _aux_uniform continuation acc bits next =
                    let nacc = 
                        Church.add (Church.add acc acc) 
                        (if (Stream.to_bool next) then 1 else 0) 
                    in
                    if (Church.gt bits 1) then
                        (_aux_uniform continuation nacc (Church.predecessor bits))
                    else (continuation nacc)
                in
                IntegerDecoder.church_stream (_aux_uniform continuation 0)

        end

        module Dna = struct
            (* First the functions to decode a base and a sequence *)
            val decode_base continuation a b =
                continuation (pair (Stream.to_bool a) (Stream.to_bool b))

            let decode_sequence continuation =
                let _add_sequence continuation len =
                    let rec _aux_decode_sequence continuation len seq =
                        let prepend_base continuation seq = 
                            let prepend_rev continuation seq base =
                                continuation (Stack.push base seq)
                            in
                            decode_base (prepend_rev continuation seq)
                        in
                        if (Church.not_zero len) then
                            (prepend_base (_aux_decode_sequence continuation len) seq)
                        else (continuation seq)
                    in
                    _aux_decode_sequence continuation len Stack.empty
                in
                IntegerDecoderC.uniform (_add_sequence continuation Stack.empty)

            (* Now we write the functions for insertions *)


            val insert continuation seq =
                let do_insertion continuation sequence position base =
                    let rec do_insertion sequence position base =
                        if (Church.not_zero position) then
                            (Stack.push (Stack.pop sequence) (do_insertion (Stack.rest
                            sequence) (Church.predecessor position) base))
                        else (Stack.push base sequence)
                    in
                    continuation (do_insertion sequence position base)
                in
                IntegerDecoderC.uniform (decode_base (do_insertion continuation seq))

            val delete continuation seq =
                let do_deletion continuation sequence position =
                    let rec do_deletion sequence position =
                        if (Church.not_zero position) then
                            (Stack.push (Stack.pop sequence) (Stack.rest sequence))
                        else (Stack.rest sequence)
                    in
                    continuation (do_deletion sequence position)
                in
                IntegerDecoderC.uniform (do_deletion continuation seq)


            val substitute continuation seq = 
                let do_substitution continuation sequence position base =
                    let rec do_substitution sequence position base =
                        if (Church.not_zero position) then
                            (Stack.push (Stack.pop sequence) 
                            (do_substitution (Stack.rest sequence) 
                                (Church.predecessor position) base))
                        else (Stack.push base (Stack.rest sequence))
                    in
                    continuation (do_substitution sequence position base)
                in
                IntegerDecoderC.uniform (decode_base (do_substitution continuation
                seq))
        end

        module Branch = struct
            val leaf continuation results stack =
                continuation 
                (Stack.push (Stack.pop stack) results) 
                (Stack.rest stack)

            val interior continuation results stack =
                continuation results (Stack.push (Stack.pop stack) stack)

            val ended continuation results stack = results
        end

        module Tree = struct

            val root continuation results = 
                let _aux_root continuation results stack seq =
                    continuation results (Stack.push seq stack)
                in
                Dna.decode_sequence (_aux_root continuation results Stack.empty)

        end

        module Phylogeny = struct
            val start continuation =
                Tree.root continuation Stack.empty
        end) ["Dna.insert", 300; "Dna.delete", 300; "Dna.substitute", 300;
        "Branch.leaf", 50; "Branch.interior", 49; "Branch.ended", 50; 
        "Tree.root", 2] compiler



let apply_model compiler = simple_indelsubstitution compiler
