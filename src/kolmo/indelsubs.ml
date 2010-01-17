let compiler = 
    Compiler.compile (OCAMLSK (** Boolean operations *)

(* [m_true] is defined as follows K x y -> x *)
let m_true = [SK K]

(* [m_false] is the opposite to true S K x y -> K y (x y) -> y *)
let m_false = [SK (S K)]

(** [m_or x y] is the boolean [x] OR [y]. The argument order in the function
* reduces the size. *)
let m_or x y = if x then x else y

let m_and x y = if x then y else x

(** [m_not x] is the negation of [x]. *)
let m_not x = if x then m_false else m_true

(** [pair a b c] is used to represent pairs of elements [a] and [b] by holding
* their computation. The idea is as follows, [z <- pair a b], then we can use
* [m_true] and [m_false] to extract either [a] or [b] by applying it to z, as
* follows [z m_true] -> [m_true a b] -> [a], and [z m_false] -> [m_false a b] ->
* [b]. *)
let pair a b c = c a b

(* From the description of [pair], it is easy to see that [first] and [second]
 * are [m_true] and [m_false] respectively. *)
let first = m_true
let second = m_false

(** Church integers representation. The number n is represented as a list of n
* trues. The list is formed by pairs. *)
module Church = struct

    (** We represent zero as a pair of false and true *)
    let zero = pair m_false m_true
    let successor x = pair m_true x                                   
    let predecessor x = x second                                
    let not_zero x = x first
    let is_zero x = m_not (not_zero x)

    let rec equal a b =
        if not_zero a then
            if not_zero b then equal (predecessor a) (predecessor b)
            else m_false
        else is_zero b

    let rec gt a b =
        if not_zero a then
            if not_zero b then gt (predecessor a) (predecessor b)
            else m_true
        else m_false

    let rec gt a b =
        if not_zero a then
            if not_zero b then gt (predecessor a) (predecessor b)
            else m_true
        else m_false

    let lt a b = gt b a

    let rec add x y = 
        if not_zero x then add (predecessor x) (successor y)
        else y

    let rec substract x y =
        if not_zero y then substract (predecessor x) (predecessor y)
        else x

    let rec multiply x y =
        if not_zero y then 
            add x (multiply x (predecessor y))
        else 0

    let log2 x = 
        let rec log2 acc cur x =
            if m_or (gt cur x) (equal cur x) then acc
            else log2 (successor acc) (add cur cur) x
        in
        log2 1 1 x

end

module Stream = struct
    (* [to_bool] converts a single S into SK and K into K. *)
    let to_bool x = x [SK S] [SK K] [SK K] m_not m_true

    (* [to_bool_c] is the same as [to_bool], but does this procedure using
    * continuation passing style, to continue processing the stream with the
    * output. *)
    let to_bool_c continuation x = 
        continuation (x [SK S] [SK K] [SK K] m_not m_true)
end

module Stack = struct
    let push x stack = pair [SK S] (pair x stack)
    let empty = pair [SK K] [SK K]
    let pop stack = stack second first
    let rest stack = stack second
    let is_empty stack = Stream.to_bool (stack first)

    let rec inv_merge a b = 
        if is_empty a then b
        else inv_merge (rest a) (push (pop a) b)

    let rec merge a b =
        if is_empty a then b
        else push (pop a) (merge (rest a) b)

    let rec length stack = 
        if is_empty stack then Church.zero
        else Church.successor (length (rest stack))
end

module IntegerDecoder = struct
    let church_stream continuation =
        let rec _church_stream continuation acc next =
            if Stream.to_bool next then 
                continuation acc next
            else
                _church_stream continuation (Church.successor acc)
        in
        _church_stream continuation 0

    let apply x y = x y

    let identity x = x

    let rec _uniform_max continuation acc bits next =
        let nacc = 
            apply
            (if (Stream.to_bool next) then Church.successor else identity)
            (Church.add acc acc) 
        in
        let next_to_execute newbits = 
            if Church.not_zero newbits then continuation nacc
            else _uniform_max continuation nacc newbits
        in
        next_to_execute (Church.predecessor bits)

end

module IndelSub = struct

    let prepend c x y = c (Stack.push x y)


let jc_decode_base substitution prepend seq previous a b =
    substitution
    (prepend 
        (pair (Stream.to_bool a) (Stream.to_bool b)))
        (Stack.rest seq)

let k2p_decode_base next prepend seq previous is_transition =
    if Stream.to_bool is_transition then 
        next
            (prepend 
                (pair (first previous) 
                    (m_not (second previous))))
            (Stack.rest seq)
    else 
        let decode_other d =
            next
                (prepend 
                    (pair (m_not (first previous)) 
                        (Stream.to_bool d)))
                (Stack.rest seq)
        in
        decode_other

let gtr_decode_base next prepend seq previous =
        let decode_another opt1 opt2 d = 
            next
                (prepend 
                    (if Stream.to_bool d then opt1 
                    else opt2))
                (Stack.rest seq)
        in
        let continue base = next (prepend base) in
        let decode_base a b c =
            if a then
                if b then
                    (* We are dealing with Adenine *)
                    if Stream.to_bool c then
                        continue
                            (pair [SK (S K)] [SK K]) (* Guanine *)
                    else 
                        decode_another 
                            (pair [SK (S K)] [SK (S K)]) (* Timine *)
                            (pair [SK K] [SK (S K)]) (* Citosine *)
                else 
                    (* We are dealing with Citosine *)
                    if Stream.to_bool c then
                        continue 
                            (pair [SK (S K)] [SK (S K)])  (* Timine *)
                    else 
                        decode_another 
                            (pair [SK K] [SK K]) (* Adenine *)
                            (pair [SK (S K)] [SK K]) (* Guanine *)
            else 
                if b then 
                    (* We are dealing with Guanine *)
                    if Stream.to_bool c then
                        continue
                            (pair [SK K] [SK K]) (* Adenine *)
                    else 
                        decode_another
                            (pair [SK K] [SK (S K)]) (* Guanine *)
                            (pair [SK (S K)] [SK (S K)]) (* Timine *)
                else 
                    (* We are dealing with Timine *)
                    if Stream.to_bool c then
                        continue
                            (pair [SK K] [SK (S K)]) (* Citosine *)
                    else 
                        decode_another
                            (pair [SK K] [SK K]) (* Adenine *)
                            (pair [SK (S K)] [SK K]) (* Guanine *)
        in
        decode_base (first previous) (second previous)


    let rec substitution decode_base continuation seq = 
        if Stack.is_empty seq then
            let do_substitution something_happened =
                if Stream.to_bool something_happened then
                    (decode_base substitution (prepend continuation) 
                        seq (Stack.pop seq))
                else 
                    substitution decode_base
                        (prepend continuation (Stack.pop seq))
                            (Stack.rest seq)
            in
            do_substitution 
        else continuation seq

let jc_substitution = substitution jc_decode_base

let k2p_substitution = substitution k2p_decode_base

let gtr_substitution = substitution gtr_decode_base


    let atomic_indel next continuation seq is_deletion = 
        if Stream.to_bool is_deletion then
            next continuation (Stack.rest seq) 
        else (* Is insertion *)
            let decode_base a b = 
                next 
                    (prepend 
                        (pair (Stream.to_bool a) (Stream.to_bool b)))
                    (Stack.rest seq)
            in
            decode_base

    let affine_indel substitution next continuation seq is_deletion = 
        if Stream.to_bool is_deletion then
            let rec extend_deletion seq should_continue =
                if Stream.to_bool should_continue then
                    extend_deletion (Stack.rest seq)
                else next continuation (Stack.rest seq)
            in
            extend_deletion
        else
            let rec extend_insertion continuation seq should_continue =
                if Stream.to_bool should_continue then
                    substitution extend_insertion (prepend continuation) 
                        (Stack.push [SK K] seq) 
                        [SK K]
                else 
                    substitution
                        next 
                        (prepend continuation)
                        (Stack.push [SK K] seq) 
                        [SK K] 
            in
            extend_insertion continuation seq [SK K]


    let rec indelsub indel substitution continuation seq =
        if Stack.is_empty seq then
            let check_insertion do_insertion = 
                if Stream.to_bool do_insertion then
                    let curry1 a b = a b in
                    indel curry1 continuation seq [SK S]
                else continuation seq
            in
            check_insertion
        else
            let do_operation nothing_happened = 
                if Stream.to_bool nothing_happened then
                    indelsub indel substitution
                        (prepend continuation (Stack.pop seq))
                            (Stack.rest seq)
                else 
                    let do_indelsub is_substitution =
                        if Stream.to_bool is_substitution then
                            substitution (indelsub indel substitution) 
                                (prepend continuation) seq (Stack.pop seq)
                        else (* We are doing a deletion or a insertion *)
                            indel (indelsub indel substitution) continuation seq
                    in
                    do_indelsub
            in
            do_operation 

    let jc_atomic_indelsub = indelsub atomic_indel jc_decode_base 

    let k2p_atomic_indelsub = indelsub atomic_indel k2p_decode_base 

    let gtr_atomic_indelsub = indelsub atomic_indel gtr_decode_base

    let affine_indelsub decode_base = indelsub (affine_indel decode_base) decode_base
    let jc_affine_indelsub = affine_indelsub jc_decode_base

    let k2p_affine_indelsub = 
        indelsub (affine_indelsub jc_decode_base) k2p_decode_base

    let gtr_affine_indelsub = 
        indelsub (affine_indelsub jc_decode_base) gtr_decode_base

end) Compiler.compiler

let jc_main = Compiler.get compiler "IndelSub.jc_substitution"
let k2p_main = Compiler.get compiler "IndelSub.k2p_substitution"
let gtr_main = Compiler.get compiler "IndelSub.gtr_substitution"

let jc_indel_main = Compiler.get compiler "IndelSub.jc_atomic_indelsub"
let jc_indel_len = Compiler.complexity compiler "IndelSub.jc_atomic_indelsub" 

let k2p_indel_main = Compiler.get compiler "IndelSub.k2p_atomic_indelsub"
let k2p_indel_len = Compiler.complexity compiler "IndelSub.k2p_atomic_indelsub" 

let gtr_indel_main = Compiler.get compiler "IndelSub.gtr_atomic_indelsub"
let gtr_indel_len = Compiler.complexity compiler "IndelSub.gtr_atomic_indelsub" 

let jc_aff_main = Compiler.get compiler "IndelSub.jc_affine_indelsub"
let jc_aff_len = Compiler.complexity compiler "IndelSub.jc_affine_indelsub" 

let k2p_aff_main = Compiler.get compiler "IndelSub.k2p_affine_indelsub"
let k2p_aff_len = Compiler.complexity compiler "IndelSub.k2p_affine_indelsub" 

let gtr_aff_main = Compiler.get compiler "IndelSub.gtr_affine_indelsub"
let gtr_aff_len = Compiler.complexity compiler "IndelSub.gtr_affine_indelsub" 

let k2p_len = Compiler.complexity compiler "IndelSub.k2p_substitution"
let gtr_len = Compiler.complexity compiler "IndelSub.gtr_substitution"
let jc_len = Compiler.complexity compiler "IndelSub.jc_substitution"
