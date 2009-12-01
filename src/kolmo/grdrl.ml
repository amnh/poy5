let compiler = Compiler.compile (OCAMLSK (** Boolean operations *)

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



(* When we are processing a stream of data, we don't want to represent true and
* false in an asymetric way. Instead, we should be able to convert K and S in
* true and false, as if they represented 1 and 0 respectively. *)
module Stream = struct
    (* [to_bool] converts a single S into SK and K into K. *)
    let to_bool x = x [SK S] [SK K] [SK K] m_not m_true
end

module Sequence = struct
    val prepend x stack = pair [SK S] (pair x stack)
    val empty = pair [SK K] [SK K]
    val head stack = stack second first
    val tail stack = stack second
    val is_empty stack = Stream.to_bool (stack first)
end

(** [m_or x y] is the boolean [x] OR [y]. The argument order in the function
* reduces the size. *)
let m_or x y = if x then x else y

module Church = struct

    (** We represent zero as a pair of false and true *)
    let zero = pair m_false m_true
    val successor x = pair m_true x                                   
    val predecessor x = x second                                
    val not_zero x = x first
    val is_zero x = m_not (not_zero x)

    val rec equal a b =
        if not_zero a then
            if not_zero b then equal (predecessor a) (predecessor b)
            else m_false
        else equal b a

    val rec gt a b =
        if not_zero a then
            if not_zero b then gt (predecessor a) (predecessor b)
            else m_true
        else m_false

    val rec gt a b =
        if not_zero a then
            if not_zero b then gt (predecessor a) (predecessor b)
            else m_true
        else m_false

    let lt a b = gt b a

    val rec add x y = 
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

module GTDRL = struct

    let prepend c x y = c (Sequence.prepend x y)
    let identity x = x

    let curry1 f x = f x

    let rec process fh sh seq position len =
        if Church.not_zero position then
            process (prepend fh (Sequence.head seq)) sh 
            (Sequence.tail seq) (Church.predecessor position)
            len
        else
            if Church.not_zero len then
                let process_next item =
                    curry1
                    (if Stream.to_bool item then 
                        process (prepend fh (Sequence.head seq)) sh 
                        (Sequence.tail seq) 
                    else 
                        process fh (prepend sh (Sequence.head seq)) 
                        (Sequence.tail seq)) (Church.predecessor len)
                    
                in
                process_next
            else sh (fh seq)

    let rec my_process continuation loglen seq another_tdrl = 
        if Stream.to_bool another_tdrl then 
            let decode_len position = 
                IntegerDecoder.church_stream 
                    (process identity (my_process continuation loglen) seq position)
            in
            IntegerDecoder._uniform_max decode_len 0 loglen
        else continuation seq

end) Compiler.compiler

let main = Compiler.get compiler "GTDRL.my_process" 
let len = Compiler.complexity compiler "GTDRL.my_process" 
