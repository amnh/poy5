(** An implementation of a standard machine for phylogenetic computations using
* compression and an S K machine *)

open Primitives
(** {2 Mutually Independent Phylogenetic Operations On DNA Sequences} *)
let m_decode_base = 
    S_K.create (SK (Callback ([m_pair] [s_to_bool (SK a)]
    [s_to_bool (SK b)])))
    (LS Callback a b)

(** Test the previous function 
    * Expected results K, K, (S K), (S K), K, (S K) 
    S_K.eval (SK ([m_decode_base] I K K [m_head]));;
    S_K.eval (SK ([m_decode_base] I K K [m_tail]));;
    S_K.eval (SK ([m_decode_base] I S S [m_tail]));;
    S_K.eval (SK ([m_decode_base] I S S [m_head]));;
    S_K.eval (SK ([m_decode_base] I K S [m_head]));;
    S_K.eval (SK ([m_decode_base] I K S [m_tail]));;
*)

let m_decode_church_integer =
    let tmp = 
        S_K.create 
        (SK ([(s_to_bool (SK next))]
            (decoder [s_church_successor (SK acc)] Callback)
            (Callback acc)))
        (LS decoder acc Callback next)
    in
    S_K.create (SK ([fixedpoint] [tmp] [m_church_zero])) []

let m_decode_binary_integer =
    let tmp =
        S_K.create 
        (SK ([s_church_not_zero (SK [s_church_predecessor (SK cnt)])]
            ([s_to_bool (SK next)] 
                (decoder add ([m_church_successor] (add acc acc))
                 Callback [s_church_predecessor (SK cnt)])
                (decoder add (add acc acc) Callback [s_church_predecessor (SK cnt)] 
                ))
            (Callback ([s_to_bool (SK next)] ([m_church_successor] (add acc
            acc)) (add acc acc)))))
        (LS decoder add acc Callback cnt next)
    in
    S_K.create (SK ([fixedpoint] [tmp] [m_church_add] [m_church_zero])) []

let m_decode_binary_integer =
    S_K.create
    (SK ( [m_decode_church_integer] ([m_decode_binary_integer] Callback)))
    (LS Callback)

(* Testing the previous functions 
to_int (S_K.eval (SK ([m_decode_binary_integer] K K K K K S K S K K S)));;
*)

let decode_function =
    let tmp =
        S_K.create
        (SK ([s_to_bool (SK next)]
            ([s_to_bool (SK (encoded [m_tail] [m_head] [m_head]))]
                (decoder (encoded [m_tail] [m_head]) Callback)
                (Callback (encoded [m_tail] [m_head] [m_tail])))
            ([s_to_bool (SK (encoded [m_tail] [m_tail] [m_head]))]
                (decoder (encoded [m_tail] [m_tail]) Callback)
                (Callback (encoded [m_tail] [m_tail] [m_tail])))))
        (LS decoder encoded Callback next)
    in
    S_K.create (SK ([fixedpoint] [tmp])) []

(* Testing the decode_function
let l1 = (SK ([m_pair] S a))
let l2 = (SK ([m_pair] S b))
let l3 = (SK ([m_pair] S c))
let join a b = (SK ([m_pair] K ([m_pair] [a] [b])))
let encoded = join (join l2 l3) l1;;
(S_K.to_string (S_K.eval (SK ([decode_function] [encoded] I S))));;
*)

let m_prepender = (* A function to prepend a base to an accumulator *)
    S_K.create (SK (Callback ([m_pair] x acc))) (LS Callback acc x)

(* Testing m_prepender 
* This function should output 0 
to_int (S_K.eval (SK ([m_prepender] [m_church_predecessor] [m_church_zero] K)));;
* This function should output 2
to_int (S_K.eval (SK ([m_prepender] [m_church_successor] [m_church_zero] K)));;
*)

let m_decode_sequence =
    let shufle_last_two = S_K.create (SK (a e f b d c)) (LS a e f b c d) in
    let tmp =
        S_K.create 
        (SK ([s_church_not_zero (s_church_predecessor (SK cnt))] 
            (m_decode_base
                (m_prepender
                    ([shufle_last_two] decoder m_decode_base m_prepender
                    Callback 
                        [s_church_predecessor (SK cnt)]) acc))
            (m_decode_base 
                (m_prepender Callback acc))))
        (LS decoder m_decode_base m_prepender Callback acc cnt)
    in
    let decoder =
        S_K.create (SK ([fixedpoint] [tmp] [m_decode_base] [m_prepender] 
        Callback ([m_pair] ([m_pair] (K S) K) K))) (LS Callback)
    in
    S_K.create (SK ( [m_decode_binary_integer] ([decoder] (Callback stack
    results accumulator))))
    (LS Callback stack results accumulator editedsequence)
(* An example of the use of this function is available in the module
* interface. *)

let merger = S_K.create (SK (a (b c))) (LS a b c)

let m_insert_delta : S_K.primitives = 
    let tmp =
        S_K.create
        (SK ([s_church_not_zero (s_church_predecessor (SK pos))]
            (insert 
                Callback
                ([m_pair] (seq [m_head]) acc)
                (seq [m_tail])
                [s_church_predecessor (SK pos)]
                base)
            (Callback ([m_pair] base acc) seq)))
        (LS insert Callback acc seq pos base)
    in
    let insert_f = S_K.create (SK ([fixedpoint] [tmp])) [] in
    S_K.create (SK ([m_decode_binary_integer] ([merger] [m_decode_base]
        ([insert_f] (Callback stack results) acc seq)))) (LS Callback stack
        results acc seq)

(**  Testing m_insert_delta
let fourth = S_K.create (SK d) (LS a b c d);;
let seq = S_K.eval (SK ([m_decode_sequence] [fourth] stack results
    accumulator editedsequence K K S K K K K K K S K));;
let empty_sequence = (SK ([Primitives.m_pair] ([Primitives.m_pair] (K S) K) K))
let res = S_K.eval (SK ([m_insert_delta] Callback stack results
    [empty_sequence] [seq] K K S K S K S));; 
let tot = S_K.eval (SK ([res] [m_head] [m_tail]));;
let tot = S_K.eval (SK ([res] [m_tail] [m_tail] [m_head] [m_tail]));;
*)

let m_delete_delta =
    let tmp =
        S_K.create
        (SK ([s_church_not_zero (s_church_predecessor (SK pos))]
            (delete 
                Callback
                ([m_pair] (seq [m_head]) acc)
                (seq [m_tail])
                [s_church_predecessor (SK pos)])
            (Callback acc (seq [m_tail]))))
        (LS delete Callback acc seq pos)
    in
    let delete_f = S_K.create (SK ([fixedpoint] [tmp])) [] in
    S_K.create (SK ([m_decode_binary_integer] 
    ([delete_f] (Callback stack results) acc seq))) (LS Callback stack
    results acc seq)

(**  Testing m_delete_delta
let res = S_K.eval (SK ([m_delete_delta] K ([m_pair] K K) B K S K));; 
let tot = S_K.eval (SK ([res] [m_head] [m_tail]));;
*)

let m_substi_delta = 
    let tmp =
        S_K.create
        (SK ([s_church_not_zero (s_church_predecessor (SK pos))]
            (substi 
                Callback
                ([m_pair] (seq [m_head]) acc)
                (seq [m_tail])
                [s_church_predecessor (SK pos)]
                base)
            (Callback 
                ([m_pair] base acc) 
                (seq [m_tail]))))
        (LS substi Callback acc seq pos base)
    in
    let substi_f = S_K.create (SK ([fixedpoint] [tmp])) [] in
    S_K.create (SK ([m_decode_binary_integer] ([merger] [m_decode_base]
        ([substi_f] (Callback stack results) acc seq)))) (LS Callback stack
        results acc seq)


(** Testing substitution
let res = S_K.eval (SK ([m_substi_delta] K ([m_pair] A K) ([m_decode_sequence] I K K S K K K K K K S K) K K S K S K S));;
let tot = S_K.eval (SK ([res] [m_head] [m_tail]));;
let tot = S_K.eval (SK ([res] [m_head] [m_head]));;
let tot = S_K.eval (SK ([res] [m_tail] [m_tail] [m_head] [m_tail]));;
*)

let rec make_encoder tree leafs =
    match tree with
    | `Node [a; b] -> 
            (SK (Pair K (Pair [make_encoder a leafs] [make_encoder b
            leafs])))
    | `Label x -> 
            (SK (Pair S [try (List.assoc x leafs) with Not_found ->
                tree]))
    | `S -> (SK (Pair S S))
    | `K -> (SK (Pair S K))
    | `Lazy x -> make_encoder (Lazy.force_val x) leafs
    | `Debugger _ 
    | `Node _ -> failwith "Illegal encoder: each node must have two leafs
and no debugging info"

let m_decode_function =
    let tmp = 
        S_K.create 
        (SK (([s_to_bool (SK (encoder [m_tail] [s_to_bool (SK next)]
        [m_head]))])
            ( decoder original_encoder 
                (encoder [m_tail] [s_to_bool (SK next)]) 
                stack results acc seq)
            ((encoder [m_tail] [s_to_bool (SK next)] [m_tail]) 
                ( decoder original_encoder original_encoder) 
                stack results acc seq)))
        (LS decoder original_encoder encoder stack results acc seq next)
    in
    S_K.create (SK (  [fixedpoint] [tmp])) []

(* Testing make_encoder and m_decode_function. The result of this function
    * call should be "Andres", as the [Callback] is selected by the [a]
    * function.
let [a; b; c] = get_lst res;;
let a = S_K.create (SK (a)) (LS a b c);;
let enc = make_encoder (SK (A B)) [("A", a)];;
let remove_first = S_K.create (SK e) (LS a b c d e);;
S_K.to_string (S_K.eval (SK ([make_list] A B C D E)));;
let enc = make_encoder (SK (A B)) ["A", remove_first];;
let res = S_K.eval (SK ([m_decode_function] [enc] [enc] stack results acc
seq K));;
S_K.to_string res;; 
*)

let m_invert_full_sequence =
    let tmp =
        S_K.create
        (SK (((seq [m_head] [m_head]) K K)
            (R  ([m_pair] (seq [m_head]) acc) (seq [m_tail]))
            acc))
        (LS R acc seq)
    in
    S_K.create (SK ([fixedpoint] [tmp] ([m_pair] ([m_pair] (K S) K) K))) []

(** Testing m_invert_full_sequence. The following should return the inverted
* sequence, that is, the input sequence in the very same order. 
let a = S_K.eval (SK ([m_decode_sequence] I K K S K K K S K K S K));;
let b = S_K.eval (SK ([m_invert_full_sequence] I [a] ([m_pair] ([m_pair] (K
S) K) K)));;
S_K.eval (SK ([b] [m_head] [m_tail]));;
S_K.eval (SK ([b] [m_head] [m_head]));;
S_K.eval (SK ([b] [m_tail] [m_head] [m_tail]));;
S_K.eval (SK ([b] [m_tail] [m_head] [m_head]));;
S_K.eval (SK ([b] [m_tail] [m_tail] [m_head] [m_tail]));;
S_K.eval (SK ([b] [m_tail] [m_tail] [m_head] [m_head]));;
S_K.eval (SK ([b] [m_tail] [m_tail] [m_tail] [m_head] [m_head]));; 
*)

let m_copy_rest =
    let tmp =
        S_K.create
        (SK ((seq [m_head] [m_head] K K) 
        (R ([m_pair] (seq [m_head]) acc) (seq [m_tail]))
        acc))
        (LS R acc seq)
    in
    S_K.create (SK ([fixedpoint] [tmp])) []

(** Testing m_copy_rest. The following should return the inverted input
* sequence 
let a = S_K.eval (SK ([m_decode_sequence] I K K S K K K S K K S K));;
let b = S_K.eval (SK ([m_copy_rest] I ([m_pair] ([m_pair] (K
S) K) K) [a] ));;
S_K.eval (SK ([b] [m_head] [m_tail]));;
S_K.eval (SK ([b] [m_head] [m_head]));;
S_K.eval (SK ([b] [m_tail] [m_head] [m_tail]));;
S_K.eval (SK ([b] [m_tail] [m_head] [m_head]));;
S_K.eval (SK ([b] [m_tail] [m_tail] [m_head] [m_tail]));;
S_K.eval (SK ([b] [m_tail] [m_tail] [m_head] [m_head]));;
S_K.eval (SK ([b] [m_tail] [m_tail] [m_tail] [m_head] [m_head]));;  *)

let m_end_branch =
    let empty_sequence = (SK ([m_pair] ([m_pair] (K S) K) K)) in
    let end_branch =
        S_K.create 
        (SK ([s_to_bool (SK next)]
        (* We have to add this item twice to the stack list and callback *)
            ( Callback ([m_pair] (copy_rest_n_invert acc seq) stack) results 
                [empty_sequence] (copy_rest_n_invert acc seq) next2)
        (* We actually reached a leaf, so we don't add it to the stack but to
        * the results and continue *)
        (([s_to_bool (SK next2)])
            (* We have reached the end, time to output the result *)
            ( [m_pair] (copy_rest_n_invert acc seq) results)
            (* We have not reached the end let's keep going *)
            (  Callback (stack [m_tail]) ([m_pair] (copy_rest_n_invert acc seq)
            results) [empty_sequence] (stack [m_head])))))
        (LS copy_rest_n_invert Callback stack results acc seq next next2)
    in
    let copy_rest_n_invert = 
        S_K.create (SK ( [m_invert_full_sequence] ([m_copy_rest] acc seq)))
        (LS acc seq)
    in
    S_K.create (SK ( [end_branch] [copy_rest_n_invert])) 
    []

(* Test copy_rest_n_invert:
let res = S_K.eval (SK ([m_end_branch] [first] stack results
[empty_sequence] [seq] S K));;
let [a; b; c] = get_lst res;;
let [x; y] = get_lst c;;
let a = S_K.eval (SK ([res] [m_tail]));;
let res = S_K.eval (SK ([copy_rest_n_invert] [empty_sequence] [seq]));;
let res = S_K.eval (SK ([m_invert_full_sequence]  [seq]));;
    S_K.eval (SK ([a] [m_head] [m_tail]));;
    S_K.eval (SK ([a] [m_head] [m_head]));;
    S_K.eval (SK ([a] [m_tail] [m_head] [m_tail]));;
    S_K.eval (SK ([a] [m_tail] [m_head] [m_head]));;
    S_K.eval (SK ([a] [m_tail] [m_tail] [m_head] [m_tail]));;
    S_K.eval (SK ([a] [m_tail] [m_tail] [m_head] [m_head]));;
    S_K.eval (SK ([seq] [m_head] [m_tail]));;
    S_K.eval (SK ([seq] [m_head] [m_head]));;
    S_K.eval (SK ([seq] [m_tail] [m_head] [m_tail]));;
    S_K.eval (SK ([seq] [m_tail] [m_head] [m_head]));;
    S_K.eval (SK ([seq] [m_tail] [m_tail] [m_head] [m_tail]));;
    S_K.eval (SK ([seq] [m_tail] [m_tail] [m_head] [m_head]));;
    S_K.eval (SK ([res] [m_head] [m_tail]));;
    S_K.eval (SK ([res] [m_head] [m_head]));;
    S_K.eval (SK ([res] [m_tail] [m_head] [m_tail]));;
    S_K.eval (SK ([res] [m_tail] [m_head] [m_head]));;
    S_K.eval (SK ([res] [m_tail] [m_tail] [m_head] [m_tail]));;
    S_K.eval (SK ([res] [m_tail] [m_tail] [m_head] [m_head]));;
    S_K.eval (SK ([res] [m_tail] [m_tail] [m_tail] [m_head]));; 
*)
