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

module Stream = struct
    (* [to_bool] converts a single S into SK and K into K. *)
    let to_bool x = x [SK S] [SK K] [SK K] m_not m_true
end

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

    val rec length stack = 
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

    let rec _uniform_max add continuation acc bits next =
        let nacc = 
            add (add acc acc) 
            (if Stream.to_bool next then 1 else 0)
        in
        if Church.not_zero (Church.predecessor bits) then 
            _uniform_max add continuation nacc (Church.predecessor bits)
        else continuation nacc

end

module DCJ = struct

    let rec invert_and_merge stack1 stack2 = 
        if Stack.is_empty stack1 then stack2
        else 
            let a_head = Stack.pop stack1 in
            let a_sign = m_not (a_head second first) in
            let a_base = a_head second second in
            let new_a_head = pair m_false (pair a_sign a_base) in
            Stack.push new_a_head (invert_and_merge (Stack.rest stack1) stack2)


    let is_chromosome_limit marker = marker first 

    let rec are_in_same_chromosome offset1 offset2 genome = 
        let first_base = Stack.pop genome in
        let genome = Stack.rest genome in
        if Church.not_zero offset1 then 
            let offset1 = Church.predecessor offset1 in
            are_in_same_chromosome
            offset1 offset2
            genome
        else 
            if Church.not_zero offset2 then
                if is_chromosome_limit first_base then 
                    m_false
                else
                    are_in_same_chromosome offset1 
                    (Church.predecessor offset2)
                    genome
            else m_true

    let apply1 f x = f x
    let marker x = pair m_true x
    let circular_mark = marker m_true
    let linear_mark = marker m_false

    let invert_in_place inv_merge stack1 stack2 genome =
        invert_and_merge stack1
            (inv_merge stack2 genome)

    let make_circular_in_place inv_merge stack1 stack2 genome =
        inv_merge stack2
            (Stack.push circular_mark 
                (inv_merge stack1 genome))


    let rec in_one_chromosome function_to_apply stack1 stack2 offset1 offset2 genome =
        let first_base = Stack.pop genome in
        let genome = Stack.rest genome in
        if Church.not_zero offset1 then 
            let offset1 = Church.predecessor offset1 in
            Stack.push first_base 
            (in_one_chromosome function_to_apply stack1 stack2 offset1 offset2 genome)
        else 
            if Church.not_zero offset2 then
                in_one_chromosome 
                function_to_apply
                (Stack.push first_base stack1)
                stack2
                offset1
                (Church.predecessor offset2)
                genome
            else 
                if is_chromosome_limit first_base then 
                    function_to_apply Stack.inv_merge stack1 stack2 
                    (Stack.push first_base genome)
                else
                    (in_one_chromosome function_to_apply stack1 
                    (Stack.push first_base stack2) 
                    offset1 offset2 genome)

    let in_one_chromosome f x offset1 offset2 offset3 genome = 
        in_one_chromosome f x x x offset1 offset2 offset3 genome

    let rec second_cut inv_merge invert_and_merge is_circular first_half second_half first_cut position2
    invert genome =
        let first_element = Stack.pop genome in
        let genome = Stack.rest genome in
        if Church.not_zero position2 then
            if is_chromosome_limit first_element then
                Stack.push (marker is_circular)
                (inv_merge first_half
                (second_cut inv_merge invert_and_merge (first_element second) 
                Stack.empty Stack.empty first_cut 
                (Church.predecessor position2) invert genome))
            else 
                second_cut inv_merge invert_and_merge is_circular 
                (Stack.push first_element first_half)
                second_half first_cut (Church.predecessor position2) invert
                genome
        else 
            (* Wow we have reached the breaking point! At last! *)
            if first_cut first then
                (* The first one was circular, this is easy *)
                invert_and_merge first_half
                (invert_and_merge (first_cut second) 
                (Stack.push first_element genome))
            else 
                (* The first cut was linear, then handling depends on the kind
                * of this one *)
                if is_circular then
                    if is_chromosome_limit first_element then
                        Stack.push linear_mark
                        (inv_merge
                        (first_cut second first)
                        (apply1
                        (if invert then invert_and_merge first_half 
                        else inv_merge second_half)
                        (apply1
                        (if invert then invert_and_merge second_half
                        else inv_merge first_half)
                        (inv_merge (first_cut second second)
                            (Stack.push first_element genome)))))
                    else 
                        second_cut inv_merge invert_and_merge is_circular first_half 
                        (Stack.push first_element second_half)
                        first_cut
                        position2 
                        invert genome
                else
                    Stack.push
                    linear_mark
                    (inv_merge first_half
                    (apply1
                    (if invert then 
                        inv_merge (first_cut second second)
                    else invert_and_merge (first_cut second first))
                    (Stack.push linear_mark
                    (apply1
                    (if invert then 
                        inv_merge (first_cut second first) 
                    else 
                        invert_and_merge (first_cut second second))
                    (Stack.push
                    first_element
                    genome)))))

    let rec first_cut invert_and_merge is_circular first_half second_half position1 
    position2 invert genome =
        let position2 = Church.predecessor position2 in
        let first_element = Stack.pop genome in
        let genome = Stack.rest genome in
        if Church.not_zero position1 then
            if is_chromosome_limit first_element then
                Stack.push (marker is_circular)
                (invert_and_merge first_half 
                    (let is_circular = first_element second in
                    first_cut invert_and_merge is_circular Stack.empty Stack.empty 
                    (Church.predecessor position1)
                    position2
                    invert 
                    genome))
            else 
                first_cut invert_and_merge is_circular (Stack.push first_element first_half) 
                second_half (Church.predecessor position1) position2 invert 
                genome
        else 
            if is_chromosome_limit first_element then
                second_cut Stack.inv_merge invert_and_merge (first_element second)
                Stack.empty Stack.empty
                (pair is_circular
                    (if is_circular then
                        Stack.inv_merge second_half 
                        (invert_and_merge first_half Stack.empty)
                    else pair first_half second_half))
                position2
                invert
                genome
            else 
                first_cut invert_and_merge is_circular first_half
                (Stack.push first_element second_half)
                position1 position2 invert genome

    let selector continuation genome position1 position2 invert = 
        if are_in_same_chromosome position1 position2 genome then
            continuation
                (apply1
                (in_one_chromosome 
                (if invert then invert_in_place else make_circular_in_place)
                Stack.empty)
                position1 position2 genome)
        else
            let head = Stack.pop genome in
            let genome = Stack.rest genome in
            let is_circular = head second in
            continuation
                (first_cut invert_and_merge is_circular Stack.empty Stack.empty 
                position1 position2 invert genome)

    let rec my_processor uniform_max continuation length genome do_nothing =
        if do_nothing then
            continuation genome
        else
            let decode_invert continuation position1 position2 invert = 
                my_processor uniform_max continuation length
                (selector (my_processor uniform_max continuation length) 
                position1 position2 (Stream.to_bool invert) )
            in
            let decode_second_position position1 =
                uniform_max (decode_invert position1) 0 length
            in
            let decode_first_position =
                uniform_max decode_second_position 0 length 
            in
            decode_first_position

    let my_processor continuation length genome do_nothing = 
        my_processor (IntegerDecoder._uniform_max Church.add) 
        continuation length genome
        do_nothing

end)

let compiler = Compiler.compile mymodule Compiler.compiler
let main = Compiler.get compiler "DCJ.my_processor"
let len = Compiler.complexity compiler "DCJ.my_processor" 
