OCAMLSK 
    let m_true = [SK K] 
    let m_false = [SK (S K)]
    let m_and y x = if x then y else x                                  
    let m_or x y = if x then x else y                          
    let m_not x = if x then m_false else m_true
    let first = m_true                                                    
    let second = m_false
    let head = first
    let tail = second
    let m_pair a b c =   c a b
    let pair = m_pair
    let zero = pair m_false [SK K]
    let successor x = pair m_true x                                   
    let predecessor x = x second                                
    let not_zero x = x head    
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
    let to_bool x = x [SK S] [SK K] [SK K] m_not m_true

    (* All of the functions that we now apply must take a callback
    * function that will be applied over the accumulator. *)
    let decode_base callback a b =
        callback (pair (to_bool a) (to_bool b))

    (* We define a vanilla function that decodes the church integer in
    * the following stream *)
    let rec decode_church_integer acc next = 
        if (to_bool next) then 
            (decode_church_integer (successor acc))
        else acc

    (* Now we define a function that decodes church integers with
    * continuation passing style *)
    let rec decode_church_integer_cb callback acc next = 
        if (to_bool next) then
            (decode_church_integer_cb callback (successor acc))
        else (callback acc)

    (* We will define the structure of a generic Huffman like decoder
    * for streams of S and K *)
    let leaf x = pair [SK K] x
    let is_leaf x = first x
    let is_internal x = m_not (is_leaf x)
    let rec basic_decoder tree next =
        if (is_leaf tree) then
            (let a = second tree in
            (a next))
        else 
            (basic_decoder ((to_bool next) (second tree)));;

    let rec basic_decoder_cb tree acc next = 
        if (is_leaf tree) then
            (let a = second tree in
            a acc next)
        else (basic_decoder_cb ((to_bool next) (second tree)) acc)
