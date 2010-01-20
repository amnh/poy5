let mymodule = 
    (OCAMLSK (** Boolean operations *)

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

module Grimm = struct
    let is_genome_limit base = base first

    let rec invert_and_merge stack1 stack2 = 
        if Stack.is_empty stack1 then stack2
        else 
            let a_head = Stack.pop stack1 in
            let a_sign = m_not (a_head second first) in
            let a_base = a_sign second second in
            let new_a_head = pair m_false (pair a_sign a_base) in
            Stack.push new_a_head (invert_and_merge (Stack.rest stack1) stack2)


    let rec invert genome stack offset1 offset2 =
        if Church.not_zero offset1 then
            Stack.push (Stack.pop genome) 
            (invert (Stack.rest genome) stack (Church.predecessor offset1)
            offset2) 
        else 
            if Church.not_zero offset2 then
                invert 
                (Stack.rest genome)
                (Stack.push (Stack.pop genome) stack) 
                offset1 
                (Church.predecessor offset2)
            else
                invert_and_merge stack genome

    let rec beginning_of_offset genome beginning position offset1 offset2 =
        if Church.not_zero offset1 then 
            beginning_of_offset (Stack.rest genome) 
            (if is_genome_limit (Stack.pop genome) then position else beginning)
            (Church.successor position) (Church.predecessor offset1) offset2
        else if Church.not_zero offset2 then
            beginning_of_offset (Stack.rest genome) 
            (if is_genome_limit (Stack.pop genome) then position else beginning)
            (Church.successor position) offset1 (Church.predecessor offset2)
        else beginning

    let rec end_of_offset genome position offset1 = 
        let reached_end_of_chromosome= 
            m_and (Church.is_zero offset1) (is_genome_limit (Stack.pop genome)) 
        in
        if reached_end_of_chromosome then position
        else
            end_of_offset (Stack.rest genome) 
            (Church.successor position) (Church.predecessor offset1)
            

    let do_grimm invert genome offset1 offset2 =
        let chromosome_of_1 = beginning_of_offset genome 0 0 offset1 0 in
        let chromosome_of_2 = beginning_of_offset genome 0 0 offset1 offset2 in
        if Church.equal chromosome_of_1 chromosome_of_2 then
            invert genome Stack.empty offset1 offset2
        else 
            let end_of_chromosome_1 = end_of_offset genome 0 offset1 in
            let choose_translocation is_simple_translocation =
                let genome = 
                    invert genome Stack.empty end_of_chromosome_1 chromosome_of_2
                in
                if Stream.to_bool is_simple_translocation then 
                    invert genome Stack.empty offset1 offset2
                else
                    let genome = 
                        invert genome Stack.empty chromosome_of_1
                        end_of_chromosome_1 
                    in
                    let genome = invert genome Stack.empty 
                        (Church.add 
                            chromosome_of_1 
                            (Church.substract end_of_chromosome_1 offset1))
                        offset2 
                    in
                    invert genome Stack.empty
                    chromosome_of_1 
                        (Church.add 
                            (Church.substract end_of_chromosome_1 offset1)
                            (Church.substract offset2 chromosome_of_2))
            in
            choose_translocation


    let rec my_process _uniform_max continuation genome length do_nothing =
        if do_nothing then continuation genome length
        else 
            let decode_operation offset1 offset2 =
                my_process _uniform_max continuation 
                (do_grimm invert genome Stack.empty offset1 offset2) length
            in
            let decode_offset2 offset1 =
                _uniform_max (decode_operation offset1) 0
                length
            in
            let decode_offset1 =
                _uniform_max decode_offset2 0 length
            in
            decode_offset1

    let grimm continuation genome length do_nothing =
        my_process IntegerDecoder._uniform_max continuation genome length
        do_nothing
end) 

let compiler =
    Compiler.compile mymodule Compiler.compiler
let main = Compiler.get compiler "Grimm.grimm"
let len = Compiler.complexity compiler "Grimm.grimm"
