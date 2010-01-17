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

(** Tandem duplication random loss *)
module TDRL = struct

    let prepend c x y = c (Sequence.prepend x y)
    let identity x = x

    let rec process fh sh seq =
        if Sequence.is_empty seq then sh (fh (Sequence.empty))
        else
            let process_next item =
                (if (Stream.to_bool item) then 
                    (process (prepend fh (Sequence.head seq)) sh 
                    (Sequence.tail seq))
                else 
                    (process fh (prepend sh (Sequence.head seq)) 
                    (Sequence.tail seq)))
                
            in
            process_next

    let rec my_process continuation seq another_tdrl = 
        if Stream.to_bool another_tdrl then 
            process identity (my_process continuation) seq
        else continuation seq

end) Compiler.compiler
let main = Compiler.get compiler "TDRL.my_process" 
let len = Compiler.complexity compiler "TDRL.my_process" 
