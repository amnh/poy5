let res =
    Compiler.compile 
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


(* When we are processing a stream of data, we don't want to represent true and
* false in an asymetric way. Instead, we should be able to convert K and S in
* true and false, as if they represented 1 and 0 respectively. *)
module Stream = struct
    (* [to_bool] converts a single S into SK and K into K. *)
    let to_bool x = x [SK S] [SK K] [SK K] m_not m_true

    (* [to_bool_c] is the same as [to_bool], but does this procedure using
    * continuation passing style, to continue processing the stream with the
    * output. *)
    let to_bool_c continuation x = 
        continuation (x [SK S] [SK K] [SK K] m_not m_true)
end


(** A module to represent and extract functions encoded in some Huffman code *)
module Decoder = struct
    (** A Huffman code uses a tree representation, we stick to it, by using
    * pairs for the nodes, and the trees themselves.  *)

    (** [leaf x] creates a tree that contains only [x] as a leaf. A leaf is a
    * pair in which the first element is K, and the second element is the
    * contents [x]. *)
    let leaf x = pair [SK K] x

    (** [node x y] joins the two decoder trees [x] and [y] in one node. To
    * differenciate it from a leaf, we associate the tag [S] as its first
    * element. *)
    let node x y = pair [SK S] (pair x y)

    (** A function that uses the boolean [next] to decide what side of a pair 
    * should be the argument to call the continuation. *)
    let choose_side continuation tree next = continuation (next tree)

    (* [generic_decoder decoder tree] takes a stream of S and K's and consumes
    * them until reaching a tree leaf. At this point it calls the function
    * contained in the tip, and passes it the decoder as the first argument, to
    * be called again as a continuation once the function is finished with the
    * computation. *)
    let rec generic_decoder decoder tree = 
        let contents_of_node = tree second in
        if Stream.to_bool (tree first) then 
            contents_of_node (generic_decoder decoder decoder)
        else 
            Stream.to_bool_c 
                (choose_side (generic_decoder decoder) contents_of_node) 
end

(** Church integers representation. The number n is represented as a list of n
* trues. The list is formed by pairs. *)
module Church = struct

    (** We represent zero as a pair of false and true *)
    let zero = pair m_false m_true
    val successor x = pair m_true x                                   
    val predecessor x = x second                                
    val not_zero x = x first
    val is_zero x = m_not (not_zero x)

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

(** Equal probabilities *)
module Dna = struct
    val decode_base continuation a b =
        continuation (pair (Stream.to_bool a) (Stream.to_bool b))

    val rec decode_sequence continuation accumulator length =
        if Church.not_zero length then 
            let prepend base = 
                decode_sequence continuation 
                    (Stack.push base accumulator) 
                    (Church.predecessor length) 
            in
            decode_base prepend 
        else continuation accumulator
end

(** Substitutions only *)
module JC69 = struct
    val rec substitute continuation acc seq =
        if Stack.is_empty seq then continuation acc
        else 
            let rest = Stack.rest seq in
            let process_next_element do_substitution =
                if Stream.to_bool do_substitution then 
                    let do_substitute base =
                        substitute continuation (Stack.push base acc) 
                        rest
                    in
                    Dna.decode_base do_substitute
                else 
                    substitute continuation 
                    (Stack.push (Stack.pop seq) acc) rest
            in
            process_next_element

end

module K2P = struct
    let adenine =  pair m_true m_true
    let guanine = pair m_true m_false
    let cytosine = pair m_false m_true
    let thymine = pair m_false m_false

    let decode_dna continuation base next =
        (* The base is what has been observed in the dna strand and is going to
        * be substitued *)
        let first_half_of_base = base first in
        if Stream.to_bool next then 
            (* We are dealing with a transition, which is more likely, so just
            * flip the second part of the input base *)
            continuation (pair first_half_of_base (m_not (base second)))
        else  
            (* We are doing a transversion, so we filp the first bit, and grab
            the next item in the input stream to decide which one *)
            let decide_transversion next = 
                continuation (pair (m_not first_half_of_base) 
                (Stream.to_bool next))
            in
            decide_transversion

    val rec substitute continuation acc seq =
        if Stack.is_empty seq then continuation seq
        else 
            let rest = Stack.rest seq in
            let process_next_element do_nothing =
                let current_base = Stack.pop seq in
                if Stream.to_bool do_nothing then 
                    substitute continuation (Stack.push current_base acc) rest
                else 
                    let do_substitute new_base =
                        substitute continuation (Stack.push new_base acc) rest
                    in
                    decode_dna do_substitute current_base
            in
            process_next_element

end

module GTR = struct
    let adenine =  pair m_true m_true
    let guanine = pair m_true m_false
    let cytosine = pair m_false m_true
    let thymine = pair m_false m_false
    
    let decode_dna continuation base next = 
        (* We will simply assume some matching in the order, just to compute 
        * the complexity of each model *)
        if (base first) then (* Adenine or Guanine *)
            if (base second) then (* Adenine *)
                if Stream.to_bool next then  (* To Guanine *)
                    continuation guanine
                else (* Based on the next one we decide what the output is *)
                    let next continuation item =
                        continuation (pair m_false (Stream.to_bool item))
                    in
                    next continuation
            else (* Guanine *)
                if Stream.to_bool next then  (* To Adenine *)
                    continuation adenine
                else 
                    let next continuation item =
                        continuation (pair m_false (Stream.to_bool item))
                    in
                    next continuation
        else
            if (base second) then (* Cytosine *)
                if Stream.to_bool next then  (* To Thymine *)
                    continuation thymine
                else 
                    let next continuation item =
                        continuation (pair m_true (Stream.to_bool item))
                    in
                    next continuation
            else (* Thymine *)
                if Stream.to_bool next then  (* To Cytosine *)
                    continuation cytosine
                else 
                    let next continuation item =
                        continuation (pair m_true (Stream.to_bool item))
                    in
                    next continuation

    val rec substitute continuation acc seq =
        if Stack.is_empty seq then continuation seq
        else 
            let rest = Stack.rest seq in
            let process_next_element what_to_do =
                let current_base = Stack.pop seq in
                if Stream.to_bool what_to_do then (* We see K, no subst *)
                    substitute continuation (Stack.push current_base acc) rest
                else 
                    let do_substitute new_base =
                        substitute continuation (Stack.push new_base acc) rest
                    in
                    decode_dna do_substitute current_base
            in
            process_next_element


end

(** Indels only, same probability of deletion or insertion for everyone  *)
module IndelsOnly = struct
    val insert continuation acc seq =
        let do_insert base =
            continuation (Stack.push base acc) seq
        in
        Dna.decode_base do_insert

    val delete continuation acc seq = 
        continuation acc (Stack.rest seq)


    val rec process_options continuation acc seq no_edition =
        if Stream.to_bool no_edition then (* We have to process an indel *)
            if Stack.is_empty seq then continuation acc
            else
                process_options continuation (Stack.push (Stack.pop seq) acc)
                (Stack.rest seq)
        else 
            if Stack.is_empty seq then 
                insert (process_options continuation) acc seq
            else
                let decide_what_to_do do_insert = 
                    if Stream.to_bool do_insert then 
                        insert (process_options continuation) acc seq
                    else 
                        process_options continuation acc (Stack.rest seq)
                in
                decide_what_to_do
end


(** Indels and Substitutions *)
module IndelsSubs_JC = struct

    let insert continuation acc seq base =
        continuation (Stack.push base acc) seq

    (* We will assume here that nothing happening is more likely still *)
    val rec process_options continuation acc seq no_edition =
        if Stream.to_bool no_edition then 
            if Stack.is_empty seq then continuation acc
            else 
                process_options continuation (Stack.push (Stack.pop seq) acc)
                (Stack.rest seq)
        else 
            (* We have one of the three pals occurring *)
            let substitution_or_indel is_substitution =
                if Stream.to_bool is_substitution then 
                    (* All the bases have the same probability *)
                    let do_substitution base = 
                        process_options continuation (Stack.push base acc) 
                        (Stack.rest seq)
                    in
                    Dna.decode_base do_substitution 
                else (* We are doing an indel, now to decide what is it *)
                    if Stack.is_empty seq then 
                        (* We are surely doing an insertion *)
                        Dna.decode_base (insert (process_options continuation)
                        acc seq) 
                    else
                        let insertion_or_deletion is_insertion =
                            if Stream.to_bool is_insertion then 
                                Dna.decode_base 
                                (insert (process_options continuation) acc seq)
                            else 
                                process_options continuation acc 
                                (Stack.rest seq)
                        in
                        insertion_or_deletion 
            in
            substitution_or_indel
end

module IndelsSubs_K2P = struct

    val rec process_options continuation acc seq no_edition = 
        if Stream.to_bool no_edition then 
            if Stack.is_empty seq then continuation acc
            else 
                process_options continuation (Stack.push (Stack.pop seq) acc)
                (Stack.rest seq)
        else 
            (* We have one of the three pals occurring *)
            let substitution_or_indel is_substitution =
                if Stream.to_bool is_substitution then 
                    (* All the bases have the same probability *)
                    let do_substitution base = 
                        process_options continuation (Stack.push base acc) 
                        (Stack.rest seq)
                    in
                    K2P.decode_dna do_substitution (Stack.pop seq)
                else (* We are doing an indel, now to decide what is it *)
                    if Stack.is_empty seq then 
                        (* We are surely doing an insertion *)
                        Dna.decode_base 
                        (IndelsSubs_JC.insert 
                            (process_options continuation) acc seq) 
                    else
                        let insertion_or_deletion is_insertion =
                            if Stream.to_bool is_insertion then 
                                Dna.decode_base 
                                (IndelsSubs_JC.insert 
                                    (process_options continuation) acc seq)
                            else 
                                process_options continuation acc (Stack.rest seq)
                        in
                        insertion_or_deletion 
            in
            substitution_or_indel

end

module IndelsSubs_GTR = struct
    let insert continuation acc seq base =
        continuation (Stack.push base acc) seq

    val rec process_options continuation acc seq do_something = 
        if Stream.to_bool do_something then 
            (* We have one of the three pals occurring *)
            let substitution_or_indel is_substitution =
                if Stream.to_bool is_substitution then 
                    (* All the bases have the same probability *)
                    let do_substitution base = 
                        process_options continuation (Stack.push base acc) 
                        (Stack.rest seq)
                    in
                    GTR.decode_dna do_substitution (Stack.pop seq)
                else (* We are doing an indel, now to decide what is it *)
                    if Stack.is_empty seq then 
                        (* We are surely doing an insertion *)
                        Dna.decode_base (insert (process_options continuation)
                        acc seq) 
                    else
                        let insertion_or_deletion is_insertion =
                            if Stream.to_bool is_insertion then 
                                Dna.decode_base (insert (process_options
                                continuation) acc seq)
                            else 
                                process_options continuation acc 
                                (Stack.rest seq)
                        in
                        insertion_or_deletion 
            in
            substitution_or_indel
        else 
            if Stack.is_empty seq then continuation acc
            else 
                process_options continuation 
                (Stack.push (Stack.pop seq) acc) (Stack.rest seq)
end

    let apply2 a b c = a b c

(** Affine Indels *)
module AffIndelsSubs_GTR = struct


    let rec insertion continuation acc seq =
        let make_choice base should_continue =
            apply2 
                (if Stream.to_bool should_continue then
                    insertion continuation a b
                else continuation a b )
                (Stack.push base acc) seq
        in
        Dna.decode_base make_choice 

    let rec deletion continuation acc seq should_continue =
        if Stack.is_empty seq then continuation acc seq
        else
            if Stream.to_bool should_continue then
                deletion continuation acc (Stack.rest seq)
            else continuation acc (Stack.rest seq)

    val rec process_options continuation acc seq no_edition = 
        if Stream.to_bool no_edition then 
            if Stack.is_empty seq then continuation acc
            else 
                process_options continuation (Stack.push (Stack.pop seq) acc)
                (Stack.rest seq)
        else 
            (* We have one of the three pals occurring *)
            let substitution_or_indel is_substitution =
                if Stream.to_bool is_substitution then 
                    (* All the bases have the same probability *)
                    let do_substitution base = 
                        process_options continuation (Stack.push base acc) 
                        (Stack.rest seq)
                    in
                    GTR.decode_dna do_substitution (Stack.pop seq)
                else (* We are doing an indel, now to decide what is it *)
                    if Stack.is_empty seq then 
                        (* We are surely doing an insertion *)
                        insertion (process_options continuation) acc seq
                    else
                        let insertion_or_deletion is_insertion =
                            apply2 (if Stream.to_bool is_insertion then 
                                insertion (process_options continuation)
                            else 
                                deletion (process_options continuation)) acc seq
                        in
                        insertion_or_deletion 
            in
            substitution_or_indel

end


(** Substitutions and Affine Indels *)
module AffineIndel_AffineSubs_GTR = struct

    let rec insertion continuation acc seq =
        let make_choice base should_continue =
            apply2
            (if Stream.to_bool should_continue then
                insertion continuation 
            else continuation) (Stack.push base acc) seq
        in
        Dna.decode_base make_choice 

    let rec deletion continuation acc seq should_continue =
        if Stack.is_empty seq then continuation acc seq
        else
            apply2 (if Stream.to_bool should_continue then
                deletion continuation 
            else continuation) acc (Stack.rest seq)

    let rec substitution continuation acc seq = 
        let do_substitution base = 
            let acc = Stack.push base acc in
            let seq = Stack.rest seq in
            let should_continue substitutions_continue =
                apply2 (if Stream.to_bool substitutions_continue then
                    substitution continuation
                else continuation) acc seq
            in
            should_continue 
        in
        GTR.decode_dna do_substitution (Stack.pop seq)


    let apply3 a b c d = a b c d

    val rec process_options continuation acc seq no_edition = 
        if Stream.to_bool no_edition then 
            if Stack.is_empty seq then continuation acc
            else 
                process_options continuation (Stack.push (Stack.pop seq) acc)
                (Stack.rest seq)
        else 
            (* We have one of the three pals occurring *)
            let substitution_or_indel is_substitution =
                if Stream.to_bool is_substitution then 
                    (* All the bases have the same probability *)
                    substitution (process_options continuation) acc seq
                else (* We are doing an indel, now to decide what is it *)
                    if Stack.is_empty seq then 
                        (* We are surely doing an insertion *)
                        insertion (process_options continuation) acc seq
                    else
                        let insertion_or_deletion is_insertion =
                            apply3
                            (if Stream.to_bool is_insertion then insertion 
                             else deletion)
                            (process_options continuation) acc seq
                        in
                        insertion_or_deletion 
            in
            substitution_or_indel

end

(** Inversion only *)
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
end

module SignedInversions = struct
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


    let negate seq = 
        pair (m_not (seq first)) (seq second)

    let rec do_invert continuation max seq =
        if Church.not_zero max then
            let mycontinuation tmp = 
                Stack.push (negate (Stack.pop seq)) (continuation tmp)
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
end

module DCJ = struct

    let is_circular x = x second
    let is_linnear x = m_not (x second)

    let invert a =
        let rec _invert x acc =
            if Stack.is_empty x then acc
            else 
                _invert (Stack.rest x)
                    (Stack.push (SignedInversions.negate (Stack.pop x)) acc)
        in
        _invert a Stack.empty

    let join_two_separate_linnear continuation chromosmes upper1 lower1 
    upper2 lower2 kind_of_join =
        let process_pair result =
            continuation (Stack.push (pair (m_false (result first))) 
                (Stack.push (pair (m_false (result second))) chromosmes))
        in
        process_pair
        (if Stream.to_bool kind_of_join then 
            (* We are joininng upper1 with upper2 *)
            (pair (prepend upper1 (invert upper2)) 
                (prepend (invert lower1) lower2))
        else (* we re joining upper1 with lower2 *)
            (pair 
                (prepend upper1 lower2) 
                (prepend upper2 lower1)))

    let join_two_separate_circular continuation chromosomes chrom1 chrom2 
    kind_of_join =
        continuation
            (Stack.push (pair m_true 
                (if Stream.to_bool kind_of_join then prepend chrom1 chrom2 
                else prepend chrom1 (invert chrom2)))
             chromosmes)

    let join_one_circular_in_two continuation chromosmes half1 half2 
    kind_of_join =
        continuation
            (if Stream.to_bool kind_of_join then
                Stack.push 
                (pair m_true (prepend (invert half1) half2))
                chromosmes
            else 
                Stack.push 
                (pair m_true half1)
                (Stack.push (pair m_true half2) chromosmes))

    let join_circular_and_linnear chromosomes upper1 lower1 circular 
    kind_of_join =
        continuation
            (Stack.push (pair m_false 
                (if Stream.to_bool kind_of_join then 
                    prepend upper1 (prepend circular lower1)
                else prepend upper1 (prepend (invert circular) lower1)))
            chromosomes)

    let join_three_linnear continuation chromosomes upper1 center lower1 
    kind_of_join =
        continuation
            (if Stream.to_bool kind_of_join then 
                Stack.push (pair m_false 
                (prepend upper1 (prepend (invert center) lower1))) chromosomes
            else 
                Stack.push (pair m_true center) 
                (Stack.push (pair m_false upper1 lower1) chromosomes))

    let join_of_circular_tips a b = 
        prepend (prepend a Stack.empty) (prepend b Stack.empty)

    let apply1 a b = a b 

    let do_second_cut continuation is_first_call head prev pos 
        is_circular_of_first head_of_first prev_of_first tail_of_first 
        is_circular tail rest =
            if Church.is_zero pos then 
                (* We are done! First merge the other chromosomes in one big
                * stack *)
                let new_chromosome = 
                    Stack.merge prev_of_first (Stack.merge prev rest) 
                in
                if is_first_call then 
                    if is_circular_first then 
                        join_one_circular_in_two 
                        continuation
                        new_chromosome 
                        (join_of_circular_tips head_of_first tail)
                        head
                    else 
                        join_three_linnear 
                        continuation
                        new_chromosome 
                        head_of_first
                        head 
                        tail
                else
                    if is_circular_first then
                        if is_circular then
                            join_two_separate_circular
                            continuation
                            new_chromosome
                            (join_of_circular_tips head_of_first tail_of_first)
                            (join_of_circular_tips head tail)
                        else 
                            join_circular_and_linnear
                            continuation
                            new_chromosome
                            head
                            tail 
                            (join_of_circular_tips head_of_first tail_of_first)
                    else 
                        if is_circular then
                            join_circular_and_linnear
                            continuation
                            new_chromosome
                            head_of_first
                            tail_of_first
                            (join_of_circular_tips head tail)
                        else 
                            join_two_separate_linnear
                            continuation
                            new_chromosome
                            head
                            tail
                            head_of_first
                            tail_of_first
            else 
                let new_pos = Church.predecessor pos in
                if Stack.is_empty tail then 
                    let new_prev = 
                        if is_first_call then prev
                        else (Stack.push (pair is_circular head) prev)
                    in
                    let new_tail = apply1 (Stack.pop rest) second in
                    let new_is_circular = apply1 (Stack.pop rest) first in
                    let new_rest = Stack.rest rest in
                    let new_head = Stack.empty in
                    do_second_cut continuation m_false 
                    new_head new_prev new_pos 
                    is_circular_of_first head_of_first prev_of_first 
                    tail_of_first new_is_circular new_tail new_rest
                else
                    let new_head = Stack.push (Stack.pop tail) in
                    let new_tail = Stack.rest tail in
                    do_second_cut continuation is_first_call 
                    new_head prev new_pos is_circular_of_first 
                    head_of_first prev_of_first tail_of_first
                    is_circular new_tail rest

    let do_initial_cut continuation pos kind head tail prev rest = 
        if Church.is_zero pos then 
            continuation kind head prev tail kind tail rest 
        else
            if Stack.is_empty tail then
                do_initial_cut continuation (Church.predecessor pos) 
                    apply1 (apply1 (Stack.pop rest) first) Stack.empty 
                    apply1 (apply1 (Stack.pop rest) second) 
                    (Stack.push this_chrom prev)
                    (Stack.rest rest)
            else 
                do_initial_cut continuation (Church.predecessor pos) kind 
                    (Stack.push (Stack.pop tail) head) (Stack.rest tail) prev
                    rest

    let process_dcj continuation log_number_of_genes genome stop =
        if Stream.to_bool stop then continuation genome
        else
            let second_cut pos =
                do_second_cut continuation m_true Stack.empty Stack.empty pos
            in
            let first_cut first_cut_pos second_cut_pos = 
                do_initial_cut 
                (second_cut second_cut_pos)
                pos
                (apply1 (Stack.pop genome) first)
                Stack.empty
                (apply1 (Stack.pop genome) second) 
                Stack.empty
                (Stack.rest genome)
            in
            IntegerDecoder.uniform_max 
            (IntegerDecoder.uniform_max first_cut_pos log_number_of_genes) 
            log_number_of_genes

end

(** Tandem duplication random loss *)
module TDRL = struct

    let prepend c x y = c (Stack.push x y)
    let identity x = x

    let rec process fh sh seq =
        if Stack.is_empty seq then sh (fh (Stack.empty))
        else
            let process_next item =
                (if (Stream.to_bool item) then 
                    (process (prepend fh (Stack.pop seq)) sh 
                    (Stack.rest seq))
                else 
                    (process fh (prepend sh (Stack.pop seq)) 
                    (Stack.rest seq)))
                
            in
            process_next

    val rec my_process continuation seq another_tdrl = 
        if Stream.to_bool another_tdrl then 
            process identity (my_process continuation) seq
        else continuation seq

end)
