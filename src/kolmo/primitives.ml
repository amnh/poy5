(** Primitive operations in an SK machine. In the following functions, any
* one starting with an m_ is intended to be a macro that can be directly
* embedded in an SK expression, while s_ functions are ocaml functions that  *)

(** {1 Atomic Computation Types}
 * 
 * In order to perform all the computations in an S K machine we have to be able
 * to delay some operations and group arguments together. In order to do this we
 * must be able to define some basic types to allow complex things to happen,
 * like - for instance - if-then statements. The basic types that we are
 * interested on are booleans, pairs, lists, and integers. Using those basic
 * units we can compute almost anything.  
 * {2 Basic Types} 
 * {3 Booleans}
 *
 * A boolean can be represented using [K] for true and [S K] for false. This
 * has some nice properties. *)
let m_true = (SK K)
let m_false = (SK (S K))
(** Basically observe that  
 * K A B -> A, ie true A B -> A, while S K A B -> K B (A * B) -> B ie. 
 * false A B -> B. This is exactly what we need. *)

(** With this ready, we can now define the not, and, and or functions *)
let s_and x y = (SK ([x] [y] [m_false]))
let s_or x y = (SK ([x] [m_true] [y]))
let s_not x = (SK ([x] (S K) K))

(** And their corresponding macros for S_K *)
let m_and = S_K.create (SK (x y [m_false])) (LS y x)
let m_or = S_K.create (SK (x [m_true] y)) (LS x y)
let m_not = S_K.create (SK (x (S K) K)) (LS x)
 
(** We want to be able to
 * have a representation that is more concise if necessary. In other words,
 * observe that the representation of false is slightly longer than that of
 * true, and for compression purposes it would be desirable not to set a special
 * burden in the false case. For this we will also define K and S and true and
 * false using the special convertion functions that follow. *)
let s_to_bool x = (SK ([x] S K K [m_not] [m_true])) 
let m_to_bool = S_K.create (SK (x S K K [m_not] [m_true])) (LS x)

(** {2 Pair and List Representation}
*
* If we have boolean representations, we can now pair representations by using
* the true and false to get the first and last elements of the pair. *)
let m_pair = S_K.create (SK (c a b)) (LS a b c)
let m_first = m_true
let m_second = m_false

(** In the same way we can define a list as just a sequence of pairs. *)
let m_head = m_first
let m_tail = m_second

(** {2 Recursion} 
 *
 * The classic fixedpoint theorem. *)
let fixedpoint = (SK (S S K (S (K (S S(S(S S K))))K)))

(** We create a convenience function to define recursive functions *)
let rec_create name f ls = 
    let tmp = S_K.create f (name :: ls) in
    S_K.create (SK ([fixedpoint] [tmp])) []

(** {2 Integer Representation}
*
* An integer has multiple possible representations. We will do the simplest and
* then define special functions for more compact ones.
*
* {3 Church Integers} 
* This is the simplest of all, consisting of a row of [i] [K]'s closed by an [S]
* to represent the integer [i]. The decoder for a church integer produces a list
* of integers *)
let m_church_zero = (SK ([m_pair] [m_false] K))
let s_church_successor x = (SK ([m_pair] [m_true] [x]))
let s_church_predecessor x = (SK ([x] [m_second]))
let s_church_not_zero x = (SK ([x] [m_head]))
let m_church_successor = S_K.create (SK ([m_pair] [m_true] x)) (LS x)
let m_church_predecessor = S_K.create (SK (x [m_second])) (LS x)
let m_church_not_zero = S_K.create (SK (x [m_head])) (LS x)

(** {4 Conversion Functions}
 * Conversion functions between [OCaml] and [S K] church representations *)
let rec of_int x =
    if x < 0 then failwith "Invalid argument" 
    else
        if x = 0 then m_church_zero
        else s_church_successor (of_int (x - 1))

let rec to_int x =
    if (SK (S K)) = (S_K.eval (SK ([x] [m_first]))) then 0
    else 1 + (to_int (S_K.eval (SK ([x] [m_second]))))

(** {4 Basic arithmetic operations} *)
let m_church_add = 
    rec_create "add" (SK (([s_church_not_zero  (SK a)])
    (add ([s_church_predecessor (SK a)]) ([s_church_successor (SK b)]))
    b)) (LS a b)

let m_church_substract =
    rec_create "substract" (SK (([s_church_not_zero (SK b)])
        (substract ([s_church_predecessor (SK a)]) ([s_church_predecessor (SK b)]))
        a)) (LS a b)

let m_church_multiply =
    rec_create "multiply" (SK (([s_church_not_zero (SK b)])
        ([m_church_add] a (multiply a ([s_church_predecessor (SK b)])))
        [m_church_zero])) (LS a b)

module type List = sig
    type ocaml_repr
    type ocaml_list
    val empty_list : S_K.primitives
    val is_not_empty : S_K.primitives -> S_K.primitives
    val head : S_K.primitives -> S_K.primitives
    val tail : S_K.primitives -> S_K.primitives
    val prepend : S_K.primitives -> S_K.primitives -> S_K.primitives
    val ml_zero_signal : S_K.primitives -> bool
    val alphabet : (S_K.primitives * ocaml_repr) list
    val to_ocaml : S_K.primitives -> ocaml_repr
    val of_ocaml : ocaml_repr -> S_K.primitives
    val to_list : S_K.primitives -> ocaml_list
    val of_list : ocaml_list -> S_K.primitives
end

module type Integer = sig
    type extras
    (* A module type to convert from a given type to ChurchIntegers which
    * are easier to deal with *)
    val of_int : extras -> int -> S_K.primitives
    val to_int : extras -> S_K.primitives -> int
    val to_church : extras -> S_K.primitives
end

module ChurchIntegers : Integer with type extras = unit = struct
    (* A representation of integers using Church's original list
    * representation with as many K as there can be *)
    type extras = unit 

    let zero = m_church_zero
    let not_zero x = m_church_not_zero
    let ml_zero_signal x = 
       (SK (S K)) = S_K.eval (not_zero x)

    let to_int () = to_int

    let of_int () = of_int

    let to_church () = SK I
end

module K_S : List with type ocaml_repr = [ `K | `SK] 
with type ocaml_list = [`K | `SK] list = struct
    let u_Hd = m_head
    let u_Tl = m_tail
        type ocaml_repr = [ `K | `SK ]
        type ocaml_list = ocaml_repr list
        let empty_list = SK (Pair (K S) K)
        let is_not_empty x = SK ([x] [u_Hd] K K)
        let ml_zero_signal x = (SK (K S)) = S_K.eval (SK ([x] [u_Hd]))
        let head x = SK ([x] [u_Hd])
        let tail x = SK ([x] [u_Tl])
        let prepend a b = SK (Pair [a] [b])
        let alphabet = [`K, `K; (SK (S K)), `SK]
        let inv_alph = List.map (fun (a, b) -> (b, a)) alphabet
        let to_ocaml a = 
            let a = S_K.eval a in
            List.assoc a alphabet
        let of_ocaml a = List.assoc a inv_alph
        let to_list lst = 
            let rec to_list lst =
                if not (ml_zero_signal lst) then 
                    let next = to_ocaml (S_K.eval (head lst)) in
                    next :: (to_list (S_K.eval (tail lst)))
                else []
            in
            to_list (S_K.eval lst) 
        let rec of_list lst = 
            List.fold_right (fun x acc ->
                prepend (of_ocaml x) acc) lst empty_list
        let complement x = SK ([x] (S K) K)
end

module Dna : List with type ocaml_repr = string with type ocaml_list =
    string list = struct
    (* A representation of DNA sequences *)
    type ocaml_repr = string
    type ocaml_list = string list
    let empty_list = S_K.eval (SK (Pair (Pair (K S)  K) K))
    let is_not_empty x = (SK ([x] [m_head] [m_head] K K))
    let ml_zero_signal x = 
       (SK (K S)) = S_K.eval (SK ([x] [m_head] [m_head]))
    let head x = SK ([x] [m_head])
    let tail x = SK ([x] [m_tail])
    let prepend x y = SK (Pair [x] [y])
    let alphabet = 
        List.map (fun (a, b) -> S_K.eval a, b)
        [(SK (Pair K K), "A"); (SK (Pair K (S K)), "C"); 
        (SK (Pair (S K) K), "G"); (SK (Pair (S K) (S K)), "T")]
    let inv_alph = List.map (fun (a, b) -> (b, a)) alphabet
    let to_ocaml a = 
        match S_K.eval (SK (Hd [a])), S_K.eval (SK (Tl [a])) with
        | `K, `K -> "A"
        | `K, `Node [`S; `K] -> "C"
        | `Node [`S; `K], `K -> "G"
        | `Node [`S; `K], `Node [`S; `K] -> "T"
        | _ -> raise Not_found

    let of_ocaml a = List.assoc a inv_alph
    let complement x =
        let negate n = S_K.eval (SK ([x] [n] (S K) K)) in
        let a = negate m_head
        and b = negate m_tail in
        S_K.eval (SK (Pair [a] [b]))
    let to_list lst = 
        let rec to_list lst =
            if not (ml_zero_signal (S_K.eval lst)) then 
                let next = to_ocaml (S_K.eval (head lst)) in
                next :: (to_list (tail lst))
            else []
        in
        to_list (S_K.eval lst) 
    let rec of_list lst = 
        List.fold_right (fun x acc ->
            prepend (of_ocaml x) acc) lst empty_list
end
module LogInt : Integer with type extras = unit = struct
    type extras = unit

    let of_int () int =
        let rec encode_in_log int acc =
            if int = 0 then acc
            else 
                encode_in_log (int lsr 1)
                (SK (Pair [if 1 = (1 land int) then (SK K) else (SK (S K))]
                [acc]))
        in
        if int < 0 then failwith "Illegal argument"
        else encode_in_log int (SK EmptyList)

    let to_int () x =
        let lst = K_S.to_list x in
        List.fold_left (fun acc item ->
            (2 * acc) +
            (match item with
            | `K -> 1 
            | `SK -> 0)) 0 lst

    let to_church () = SK DecodeInt 

end


module FixedMax : Integer with type extras = int = struct

    type extras = int 

    let of_int max int =
        let rec encode_in_log max int acc =
            if max = 0 then acc
            else 
                encode_in_log (max lsr 1) (int lsr 1) 
                (SK (Pair [if 1 = (1 land int) then (SK K) else (SK (S K))]
        [acc]))
        in
        if int < 0 || max < 0 || max < int then failwith "Illegal argument"
        else encode_in_log max int (SK K)

    let to_church max =
        let rec aux_len max acc =
            if max = 0 then acc 
            else aux_len (max lsr 1) (acc + 1) 
        in
        if max < 0 then failwith "Illegal argument" 
        else
            let counter = ChurchIntegers.of_int () (aux_len max 0) in
            let tmp_rec = 
                S_K.create
                (SK ([s_church_not_zero (SK cnt)]
                    (R [s_church_predecessor (SK cnt)] 
                    ((int [m_head]) 
                        (Add [ChurchIntegers.of_int () 1] (Multiply
                        [ChurchIntegers.of_int () 2] acc)) 
                        (Multiply [ChurchIntegers.of_int () 2] acc))
                        (int [m_tail]))
                    acc))
                (LS R cnt acc int)
            in
            S_K.create 
            (SK (Fixedpoint [tmp_rec] [counter] [m_church_zero])) 
            []

    let to_int max int = 
        ChurchIntegers.to_int ()
        (S_K.eval (SK ([to_church max] [int])))

end

module FixedRange : Integer with type extras = (int * int) =
    struct
        type min = int
        type max = int
        type extras = (min * max)

        let of_int (min, max) int =
            if min > max || int < min || int > max then 
                failwith "Illegal argument"
            else
                FixedMax.of_int (max - min) int

        let to_int (min, max) int =
            min + (FixedMax.to_int (max - min) int)

        let to_church (min, max) =
            S_K.create
            (SK (Add [ChurchIntegers.of_int () min] ([FixedMax.to_church (max -
            min)] pos))) (LS pos)
    end
