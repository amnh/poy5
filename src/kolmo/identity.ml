let mymodule = 
    (OCAMLSK (** Boolean operations *)

(* [m_true] is defined as follows K x y -> x *)
let m_true = [SK K]

(* [m_false] is the opposite to true S K x y -> K y (x y) -> y *)
let m_false = [SK (S K)]

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
let pair a b c = c a b

(* From the description of [pair], it is easy to see that [first] and [second]
 * are [m_true] and [m_false] respectively. *)
let first = m_true
let second = m_false

module Stream = struct
    (* [to_bool] converts a single S into SK and K into K. *)
    let to_bool x = x [SK S] [SK K] [SK K] m_not m_true

    (* [to_bool_c] is the same as [to_bool], but does this procedure using
    * continuation passing style, to continue processing the stream with the
    * output. *)
    let to_bool_c continuation x = 
        continuation (x [SK S] [SK K] [SK K] m_not m_true)
end
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

    let rec _uniform_max continuation acc bits next =
        let nacc = 
            Church.add (Church.add acc acc) 
            (if (Stream.to_bool next) then 1 else 0)
        in
        if Church.equal bits 1 then continuation nacc
        else _uniform_max continuation nacc (Church.predecessor bits)

    let uniform_max continuation = _uniform_max continuation 0 

    let uniform continuation = church_stream (uniform_max continuation)

    let uniform_min_max continuation min max =
        let my_continuation decoded_integer = 
            continuation (Church.add min decoded_integer)
        in
        _uniform_max my_continuation max

end
module Stack = struct
    val push x stack = pair [SK S] (pair x stack)
    val empty = pair [SK K] [SK K]
    val pop stack = stack second first
    val rest stack = stack second
    val is_empty stack = Stream.to_bool (stack first)
    val rec merge a b = 
        if is_empty a then b
        else merge (rest a) (push (pop a) b)
    val rec length stack = 
        if is_empty stack then Church.zero
        else Church.successor (length (rest stack))
end

let rec aux_identity n stack =
    if Church.not_zero n then 
        aux_identity (Church.predecessor n) (Stack.push n stack)
    else stack

let chromosome_identity mechanism =
    let create_initial n = 
        mechanism (aux_identity n Stack.empty)
    in
    IntegerDecoder.uniform create_initial

let rec aux_signed_identity n stack = 
    if Church.not_zero n then 
        aux_signed_identity (Church.predecessor n) 
        (Stack.push (pair m_true n) stack)
    else stack

let chromosome_signed_identity mechanism = 
    let create_initial n =
        mechanism (aux_signed_identity n Stack.empty)
    in
    IntegerDecoder.uniform create_initial
    )

let compiler = Compiler.compile mymodule Compiler.compiler
let chromosome_identity = Compiler.get compiler "chromosome_identity"
let chromosome_signed_identity = 
    Compiler.get compiler "chromosome_signed_identity"
let len_chromosome_identity = 
    Compiler.complexity compiler "chromosome_identity"
let len_chromosome_signed_identity =
    Compiler.complexity compiler "chromosome_signed_identity"
