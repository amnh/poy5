let mymodule = (OCAMLSK
let m_true = [SK K]
let m_false = [SK (S K)]
let m_and y x = if x then y else x
let m_or x y = if x then x else y
let m_not x = if x then m_false else m_true

let first = m_true
let second = m_false
val pair a b c = c a b

module Stream = struct
    (* We don't use the continuation in this version *)
    let to_bool continuation x = x [SK S] [SK K] [SK K] m_not m_true
end

module Stack = struct
    val push x stack = pair [SK S] (pair x stack)
    val empty = pair [SK K] [SK K]
    val pop stack = first (second stack)
    val rest stack = second stack
    val is_empty stack = Stream.to_bool (first stack)
end


module Church = struct
    let zero = pair m_false [SK K]
    val successor x = pair m_true x                                   
    val predecessor x = x second                                
    val not_zero x = x first
    let rec equal a b =
        if (m_and (not_zero a) (not_zero b)) then
            (equal (predecessor a) (predecessor b))
        else
            (if (m_not (m_or (not_zero a) (not_zero b))) then
                m_true
            else m_false)
    let rec gt a b =
        if (m_and (not_zero a) (not_zero b)) then
            (gt (predecessor a) (predecessor b))
        else 
            (m_and (not_zero a) (m_not (not_zero b)))
    let rec lt a b =
        if (m_and (not_zero a) (not_zero b)) then
            (lt (predecessor a) (predecessor b))
        else (m_and (m_not (not_zero a)) (not_zero a))
    let rec add x y = 
        if (not_zero x) then (add (predecessor x) (successor y))
        else y
    let rec substract x y =
        if (not_zero y) then (substract (predecessor x) (predecessor y))
        else x
    let rec multiply x y =
        if (not_zero y) then 
            (add x (multiply x (predecessor y)))
        else 0

    let log2 x = 
        let rec log2 acc cur x =
            if (m_or (gt cur x) (equal cur x)) then acc
            else (log2 (successor acc) (add cur cur) x)
        in
        log2 0 1 x

end

module Array = struct
    let rec create n =
        if Church.not_zero n then 
            Stack.push Stack.empty (create (Church.predecessor n))
        else Stack.empty

    let rec insert position element arr =
        if Church.not_zero position then 
            Stack.push (Stack.pop arr) 
                (insert (Church.predecessor position) element (Stack.rest arr))
        else 
            Stack.push (Stack.push element (Stack.pop arr)) (Stack.rest arr)
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

module Network = struct

    let rec branch_processor _uniform_max max_levels continuation mechanism levels parent seq =
        let add_sequence_in_position pos seq = 
            continuation (Array.insert pos seq levels)
        in
        let process_position pos = 
            mechanism (add_sequence_in_position pos) parent seq 
        in
        _uniform_max process_position 0 max_levels

    let rec process_level _uniform_max process_network mechanism max_levels parents current levels =
        if Stack.is_empty current then process_network levels
        else 
            let next = Stack.pop current in
            branch_processor max_levels
                (process_level _uniform_max process_network 
                mechanism max_levels parents (Stack.rest current)) 
                mechanism parents next

    let rec process_network _uniform_max mechanism max_levels current_level parents levels =
        if Church.not_zero current_level then 
            process_level _uniform_max
                (process_network _uniform_max mechanism max_levels 
                    (Church.predecessor current_level) (Stack.pop levels)) 
                mechanism max_levels parents (Stack.pop levels) (Stack.rest levels)
        else 
            Stack.pop levels

    let network _uniform_max mechanism root =
        let process_levels levels = 
            process_network _uniform_max mechanism levels levels Stack.empty 
            (Array.insert 0 root (Array.insert 0 root (Array.create levels)))
        in
        IntegerDecoder.church_stream (_uniform_max process_levels 0)

    let network mechanism root = 
        network IntegerDecoder._uniform_max mechanism root

end) 


let compiler =
    Compiler.compile mymodule Compiler.compiler
let main = Compiler.get compiler "Network.network"
let len = Compiler.complexity compiler "Network.network"
