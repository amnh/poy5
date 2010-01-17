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



module Sequence = struct
    val prepend x stack = pair [SK S] (pair x stack)
    val empty = pair [SK K] [SK K]
    val head stack = stack second first
    val tail stack = stack second
    val is_empty stack = Stream.to_bool (stack first)
end

(** Tandem duplication random loss *)
module Tandem = struct

    let rec do_duplicate toduplicate seq offset1 offset2 =
        if Church.not_zero offset1 then
            Stack.push (Stack.pop seq)
            (do_duplicate toduplicate (Stack.rest seq) offset1 
                (Church.predecessor offset2))
        else if Church.not_zero offset2 then
            Stack.push (Stack.pop seq)
            (do_duplicate (Stack.push (Stack.pop seq) toduplicate)
                (Stack.rest seq) offset1 (Church.predecessor offset2))
        else Stack.inv_merge toduplicate seq

    let decode_positions_and_duplicate continuation decoder bits seq =
        let apply_offsets offset1 offset2 = 
            continuation (do_duplicate Stack.empty seq offset1 offset2)
        in
        let decode_offset2 offset1 =
            decoder (apply_offsets offset1) 0 bits 
        in
        let decode_offset1 =
            decoder decode_offset2 0 bits
        in
        decode_offset1

    let rec duplicate continuation bits seq another_duplication = 
        if Stream.to_bool another_duplication then 
            decode_positions_and_duplicate (duplicate continuation bits) 
            IntegerDecoder._uniform_max bits seq
        else continuation bits seq

end) Compiler.compiler
let main = Compiler.get compiler "Tandem.duplicate"
let len = Compiler.complexity compiler "Tandem.duplicate" 
