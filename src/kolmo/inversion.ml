let mymodule = 
    (OCAMLSK (** Boolean operations *)

(* [m_true] is defined as follows K x y -> x *)
let m_true = [SK K]

(* [m_false] is the opposite to true S K x y -> K y (x y) -> y *)
let m_false = [SK (S K)]

(** [m_and x y] is the boolean [x] AND [y]. The argument order in
* the function reduces the size *)
let m_and y x = if x then y else x

(** [m_or x y] is the boolean [x] OR [y]. The argument order in the function
* reduces the size. *)
let m_or x y = if x then x else y

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

(** Church integers representation. The number n is represented as a list of n
* trues. The list is formed by pairs. *)
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
        else is_zero b


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

module Stream = struct
    (* [to_bool] converts a single S into SK and K into K. *)
    val to_bool x = x [SK S] [SK K] [SK K] m_not m_true

    (* [to_bool_c] is the same as [to_bool], but does this procedure using
    * continuation passing style, to continue processing the stream with the
    * output. *)
    let to_bool_c continuation x = 
        continuation (x [SK S] [SK K] [SK K] m_not m_true)
end

module Stack = struct
    val push x stack = pair [SK S] (pair x stack)
    val empty = pair [SK K] [SK K]
    val pop stack = stack second first
    val rest stack = stack second
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

    val rec _uniform_max continuation acc bits next =
        let nacc = 
            Church.add (Church.add acc acc) 
            (if (Stream.to_bool next) then 1 else 0)
        in
        if Church.equal bits 1 then continuation nacc
        else _uniform_max continuation nacc (Church.predecessor bits)

    val uniform_max continuation = _uniform_max continuation 0 

    let uniform continuation = church_stream (uniform_max continuation)

    let uniform_min_max continuation min max =
        let my_continuation decoded_integer = 
            continuation (Church.add min decoded_integer)
        in
        _uniform_max my_continuation max

end
module Inversions = struct
    let prepend continuation x y = continuation (Stack.push x y)
    let rec interval continuation min max cnt = 
        if Church.is_zero cnt then continuation min max
        else 
            let dif = Church.substract max min in
            if Church.gt cnt dif then
                interval continuation (Church.successor min) max
                    (Church.substract cnt dif)
            else interval continuation min (Church.predecessor max)
                    (Church.predecessor cnt)


    let rec do_invert continuation max seq =
        if Church.not_zero max then
            let mycontinuation tmp = 
                Stack.push (Stack.pop seq) (continuation tmp)
            in
            do_invert mycontinuation (Church.predecessor max) (Stack.rest seq)
        else 
            continuation seq

    let identity x = x

    val rec invert continuation seq min max =
        if Church.is_zero min then continuation (do_invert identity max seq)
        else 
            let mycontinuation tmp =
                continuation (Stack.push (Stack.pop seq) tmp)
            in
            invert mycontinuation (Stack.rest seq) (Church.predecessor min)
            (Church.predecessor max)

    val rec my_process continuation length sequence do_nothing =
        if Stream.to_bool do_nothing then
            continuation length sequence
        else
            let internal_wrapper bits_to_read =
                IntegerDecoder.uniform_max 
                (interval (invert (my_process continuation length) sequence)
                Church.zero bits_to_read) bits_to_read
            in
            internal_wrapper (Church.add length length)
end) 

let compiler =
    Compiler.compile mymodule Compiler.compiler

let decoder = 
    Compiler.compile_decoder mymodule ["my_process", 100] Compiler.compiler

let len = 
    let res, _, _ = Compiler.tree_of_decoder decoder in
    List.length (S_K.s_encode res)
