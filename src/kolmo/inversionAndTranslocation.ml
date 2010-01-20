let mymodule =
    (OCAMLSK (** Boolean operations *)

(* [m_true] is defined as follows K x y -> x *)
let m_true = [SK K]

(* [m_false] is the opposite to true S K x y -> K y (x y) -> y *)
let m_false = [SK (S K)]

(** [m_not x] is the negation of [x]. *)
let m_not x = if x then m_false else m_true

(** [pair a b c] is used to represent pairs of elements [a] and [b] by holding
* their computation. The idea is as follows, [z <- pair a b], then we can use
* [m_true] and [m_false] to extract either [a] or [b] by applying it to z, as
* follows [z m_true] -> [m_true a b] -> [a], and [z m_false] -> [m_false a b] ->
* [b]. *)
val pair a b c = c a b

(* From the description of [pair], it is easy to see that [first] and [second]
 * are [m_true] and [m_false] respectively. *)
let first = m_true
let second = m_false

module Stream = struct
    (* [to_bool] converts a single S into SK and K into K. *)
    let to_bool x = x [SK S] [SK K] [SK K] m_not m_true

end

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
        else equal b a

    let rec add x y = 
        if not_zero x then add (predecessor x) (successor y)
        else y

end

module Stack = struct
    val push x stack = pair [SK S] (pair x stack)
    val empty = pair [SK K] [SK K]
    val pop stack = stack second first
    val rest stack = stack second
    val is_empty stack = Stream.to_bool (stack first)

    val rec inv_merge a b = 
        if is_empty a then b
        else inv_merge (rest a) (push (pop a) b)

    let rec merge a b =
        if is_empty a then b
        else push (pop a) (merge (rest a) b)
end

module IntegerDecoder = struct

    let rec _uniform_max add continuation acc bits next =
        let nacc = 
            add (add acc acc) 
            (if Stream.to_bool next then 1 else 0)
        in
        if Church.not_zero (Church.predecessor bits) then 
            _uniform_max add continuation nacc (Church.predecessor bits)
        else continuation nacc

end

module IT = struct

    let rec invert_and_merge stack1 stack2 = 
        if Stack.is_empty stack1 then stack2
        else 
            let a_head = Stack.pop stack1 in
            let a_sign = m_not (a_head second first) in
            let a_base = a_head second second in
            let new_a_head = pair m_false (pair a_sign a_base) in
            Stack.push new_a_head (invert_and_merge (Stack.rest stack1) stack2)


    let rec do_inversion_translocation genome stack offset1 offset2 do_inversion offset3 =
        if Church.not_zero offset1 then
            Stack.push (Stack.pop genome) 
            (do_inversion_translocation (Stack.rest genome) stack (Church.predecessor offset1)
            offset2 do_inversion offset3) 
        else 
            if Church.not_zero offset2 then
                do_inversion_translocation 
                (Stack.rest genome)
                (Stack.push (Stack.pop genome) stack) 
                offset1 
                (Church.predecessor offset2)
                do_inversion
                offset3
            else
                if Stream.to_bool do_inversion then
                    (* We ignore the contents of offset3 *)
                    invert_and_merge stack genome
                else 
                    if Church.not_zero offset3 then
                        Stack.push (Stack.pop genome)
                        (do_inversion_translocation (Stack.rest genome) stack offset1
                        offset2 do_inversion (Church.predecessor offset3))
                    else
                        Stack.inv_merge stack genome

    let rec my_process process _uniform_max continuation genome length do_nothing =
        if do_nothing then continuation genome length
        else 
            let decode_offset3 offset1 offset2 do_inversion offset3 =
                my_process process _uniform_max continuation 
                (process genome Stack.empty offset1 offset2 do_inversion offset3)
                length
            in
            let decode_operation offset1 offset2 do_inversion =
                if Stream.to_bool do_inversion then 
                    my_process process _uniform_max continuation 
                    (process genome Stack.empty offset1 offset2
                    do_inversion [SK K]) length
                else 
                    _uniform_max (decode_offset3 offset1 offset2
                    do_inversion) 0 length 
            in
            let decode_offset2 offset1 =
                _uniform_max (decode_operation offset1) 0
                length
            in
            let decode_offset1 =
                _uniform_max decode_offset2 0 length
            in
            decode_offset1

    let my_process continuation genome length do_nothing =
        my_process do_inversion_translocation 
        (IntegerDecoder._uniform_max Church.add) 
        continuation genome length do_nothing

end)

let compiler = Compiler.compile mymodule Compiler.compiler
let main = Compiler.get compiler "IT.my_process" 
let len = Compiler.complexity compiler "IT.my_process"
