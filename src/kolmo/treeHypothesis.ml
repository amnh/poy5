let mymodule = 
    (OCAMLSK 

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

end
module Tree = struct

    let rec process_tree mechanism stack results sequence next_node_is_leaf =
        let continuation sequence =
            if Stream.to_bool next_node_is_leaf then
                if Stack.is_empty stack then (Stack.push sequence results)
                else
                    process_tree mechanism (Stack.rest stack) results
                        (Stack.pop stack)
            else
                process_tree mechanism (Stack.push sequence stack) results
                sequence
        in
        mechanism continuation

    let start_processing mechanism sequence =
        process_tree mechanism (Stack.push sequence Stack.empty)
        Stack.empty sequence

end)

let compiler = Compiler.compile mymodule Compiler.compiler
let mainf = "Tree.start_processing"
let main = Compiler.get compiler mainf
let len = Compiler.complexity compiler mainf
