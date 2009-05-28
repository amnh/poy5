
module Encodings :
  sig
    type bit = Zero | One
    type encoding = bit list
    type natural
    val to_nat : int -> natural
    val of_nat : natural -> int
    val l : natural -> natural
    val binary : natural -> encoding
    val e : natural -> natural -> encoding
    val e_2 : natural -> encoding
    val hat : natural -> encoding
    val decode : encoding -> natural
    val huffman :
      ('a * float) list -> (encoding -> 'a list) * ('a list -> encoding)
    val huffman_tree : ('a * float) list -> 'a option Parser.Tree.t
    val geometric : int -> int -> float -> (int * float) list
  end

module S_K :
  sig
    exception Illegal_Expression of string Parser.Tree.t list
    type primitives = [`S | `K | `Label of string | `Node of primitives list |
    `Debugger of string | `Lazy of primitives Lazy.t]
    val debug : bool ref
    val universe : (string, primitives) Hashtbl.t
    val of_string : string -> primitives
    val expand : 
        ?except:string list -> primitives -> primitives
    val to_string : primitives -> string
    val sk_define : string -> primitives -> unit
    val simplify : primitives -> primitives list
    val reduce : primitives -> primitives
    val evaluate : string -> primitives
    val eval : primitives -> primitives
    val test : primitives list -> primitives
    val def : string -> string list -> primitives -> unit
    val s_encode : primitives -> Encodings.bit list
    val s_decode : Encodings.bit list -> primitives
    val sk_define_interpreted : string -> string list -> primitives -> unit
    val create : primitives -> string list -> primitives
  end

module Compiler : sig

    type compiler 
    val compiler : compiler
    type 'a kolmo_function =
        [ `Module of (name * ('a kolmo_function list))
        | `LetVal of (name * arguments * 'a definition)
        | `RecVal of (name * arguments * 'a definition) 
        | `Let of (name * arguments * 'a definition)
        | `Rec of (name * arguments * 'a definition) ]
    and name = string
    and arguments = string list 
    and 'a definition =
        [ `Value of 'a values
        | `IfElse of ('a condition * 'a definition * 'a definition)
        | `Letin of ('a kolmo_function * 'a definition) ]
    and 'a condition = [ `Condition of 'a definition ]
    and 'a values =
        [ `Apply of (string list * string * ('a definition list))
        | `Integer of string 
        | `Expr of 'a ]

    type sk_function = S_K.primitives kolmo_function
    val compile_decoder : S_K.primitives kolmo_function list ->  
        (string * int) list -> compiler -> compiler
    val compile : S_K.primitives kolmo_function list -> compiler -> compiler
    val clear : compiler -> compiler
    val evaluate : S_K.primitives definition -> compiler -> S_K.primitives
    val get : compiler -> string -> S_K.primitives
    val decoder : compiler -> (string * int * S_K.primitives) list
    val final_code : compiler -> S_K.primitives list All_sets.IntegerMap.t
    val tree_of_decoder : compiler -> 
        S_K.primitives * (string * int * S_K.primitives) list * 
        S_K.primitives list All_sets.IntegerMap.t
end

(** Primitive operations in an SK machine *)
module Primitives : sig
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
     * {2 Booleans}
     *
     * A boolean can be represented using [K] for true and [S K] for false. This
     * has some nice properties: observe that  
     * K A B -> A, ie true A B -> A, while S K A B -> K B (A * B) -> B ie. 
     * false A B -> B. This is exactly what we need. *)

    val m_true : S_K.primitives
    val m_false : S_K.primitives
    val s_and : S_K.primitives -> S_K.primitives -> S_K.primitives
    val s_or : S_K.primitives -> S_K.primitives -> S_K.primitives
    val s_not : S_K.primitives -> S_K.primitives
    val m_and : S_K.primitives
    val m_or : S_K.primitives
    val m_not : S_K.primitives

    (** {2 Short Booleans} 
     *
     * We want to be able to
     * have a representation that is more concise if necessary. In other words,
     * observe that the representation of false is slightly longer than that of
     * true, and for compression purposes it would be desirable not to set a special
     * burden in the false case. For this we will also define K and S and true and
     * false using the special convertion functions that follow. *)

    val s_to_bool : S_K.primitives -> S_K.primitives
    val m_to_bool : S_K.primitives

    (** {2 Pair and List Representation}
    *
    * If we have boolean representations, we can now pair representations by using
    * the true and false to get the first and last elements of the pair. *)
    val m_pair : S_K.primitives
    val m_first : S_K.primitives
    val m_second : S_K.primitives

    (** In the same way we can define a list as just a sequence of pairs. *)

    val m_head : S_K.primitives
    val m_tail : S_K.primitives

    (** {2 Recursion} *)

    (** The classic fixedpoint theorem. *)
    val fixedpoint : S_K.primitives

    (** [rec_create name def lst] is the same as [S_K.create (SK ([fixedpoint]
     * [S_K.create def (name :: lst)]))]. It is A convenience function to 
     * define recursive SK functions.*)
    val rec_create : string -> S_K.primitives -> string list -> S_K.primitives

    (* {3 Church Integers} 
    * This is the simplest of all, consisting of a row of [i] [K]'s closed by an [S]
    * to represent the integer [i]. The decoder for a church integer produces a list
    * of integers *)

    val m_church_zero : S_K.primitives
    val m_church_successor : S_K.primitives
    val m_church_predecessor : S_K.primitives
    val m_church_not_zero : S_K.primitives
    val s_church_successor : S_K.primitives -> S_K.primitives
    val s_church_predecessor : S_K.primitives -> S_K.primitives
    val s_church_not_zero : S_K.primitives -> S_K.primitives
    val m_church_add : S_K.primitives
    val m_church_substract : S_K.primitives
    val m_church_multiply : S_K.primitives


    (** {2 Utility Modules}
     * The following set of modules provide simplified utilities for conversion
     * between OCaml types and S K repesentations widely used. *)

    (** {3 Abstract Signatures} *)

    (** Standard representation of Lists. *)
    module type List = sig
        (** A list could be represented in different ways depending on its
        * contents. This signature specifies the functions that need to be
        * provided to have full support for a particular kind of list. *)

        (** The ocaml representation of a single element in the list. *)
        type ocaml_repr

        (** The ocaml representation of the list of elements in ocaml_repr. *)
        type ocaml_list

        (** The representation of an empty list *)
        val empty_list : S_K.primitives

        (** SK boolean to check whether or not a list is empty. The resulting
        * function, when reduced, should produce [SK K] or [SK (S K)]. *)
        val is_not_empty : S_K.primitives -> S_K.primitives

        (** Lists's head (first element), and tail (the rest excepting the last
        * element). These functions do not check if the list is the empty_list
        * or not, as exceptions can't be handled by [S K]. *)
        val head : S_K.primitives -> S_K.primitives
        val tail : S_K.primitives -> S_K.primitives

        (** [prepend x y] produces the SK representation of [x :: y]. *)
        val prepend : S_K.primitives -> S_K.primitives -> S_K.primitives

        (** {2 OCaml conversions} *)

        (** [ml_zero_signal x] is true if [S_K.eval (SK ([is_not_empty x])) = SK
        * (S K)], otherwise its true. *)
        val ml_zero_signal : S_K.primitives -> bool

        (** The list of elements accepted inside the list type represented by
        * the functionality of the module. *)
        val alphabet : (S_K.primitives * ocaml_repr) list

        (** [to_ocaml x] returns the ocaml representation of the SK expression
        * belonging to the [alphabet]. *)
        val to_ocaml : S_K.primitives -> ocaml_repr

        (** [of_ocaml x] is symmetric to [to_ocaml]. *)
        val of_ocaml : ocaml_repr -> S_K.primitives

        (** [to_list x] returns the ocaml representation of the SK list [x]
         * containing elements of the [alphabet]. *)
        val to_list : S_K.primitives -> ocaml_list

        (** [of_list x] is symmetric to [to_list]. *)
        val of_list : ocaml_list -> S_K.primitives

    end

    (** Standard representation of integers. *)
    module type Integer = sig
        type extras
        (* A module type to convert from a given type to ChurchIntegers and
        * other representations which are easier to deal with *)

        val of_int : extras -> int -> S_K.primitives
        val to_int : extras -> S_K.primitives -> int
        val to_church : extras -> S_K.primitives
    end

    (** {3 Concrete Lists and Integers} *)

    (** An implementation of a list for DNA sequences *)
    module Dna : List with type ocaml_repr = string with type ocaml_list =
        string list

    (** An implementation of a list holding only [(SK K)] and [SK (S K)]. *)
    module K_S : List with type ocaml_repr = [`K | `SK ] with type ocaml_list =
        [`K | `SK ] list

    (** An implementation of [Integer] using church's representation, the same
    * as that in the [Primitives] module. *)
    module ChurchIntegers : Integer with type extras = unit


    (** An implementation of [Integer] for logarithmic integers. It is optimal
    * for unbounded integers with equal probabilities. *)
    module LogInt : Integer with type extras = unit

    module FixedMax : Integer with type extras = int

    module FixedRange : Integer with type extras = (int * int)
end

(** A module with an SK machine for Independent Phylogenetic events
including only atomic substitutions, insertions, and deletions.
*)
module IndPM : sig

    (** [m_decode_sequence Callback stack results acc edited] reads the stream
    * of a sequence [s] and as a result executes [Callback stack results acc s].
    * The format of the encoded sequence is a logarithmic integer [l]
    * representing the length of the sequence to be decoded, followed by
    * [l] bases one after the other. The resulting sequence is closed by the
    * empty sequence represented by [([m_pair] ([m_pair] (K S) K) K)]. The input
    * sequence [edited] is discarded in the [Callback] *)
    val m_decode_sequence : S_K.primitives
    (** The following is an example of the use of [m_decode_sequence] in an SK
    * machine:
        {[
        let fourth = S_K.create (SK d) (LS a b c d);;
        let seq = S_K.eval (SK ([m_decode_sequence] [fourth] stack results
        accumulator editedsequence K K S K K K K K K S K));;
        S_K.eval (SK ([res] [m_head] [m_tail]));;
        S_K.eval (SK ([res] [m_head] [m_head]));;
        S_K.eval (SK ([res] [m_tail] [m_head] [m_tail]));;
        S_K.eval (SK ([res] [m_tail] [m_head] [m_head]));;
        S_K.eval (SK ([res] [m_tail] [m_tail] [m_head] [m_tail]));;
        S_K.eval (SK ([res] [m_tail] [m_tail] [m_head] [m_head]));;
        S_K.eval (SK ([res] [m_tail] [m_tail] [m_tail] [m_head]));;]}
    *)

    (** [m_insert_delta Callback stack results acc edited] reads a logarithmic
    * integer [l] and a base [b] to be inserted in the sequence [edited]. The
    * function moves [l - 1] bases in [edited] to [acc] and adds [b] to [acc] to
    * produce a new longer [acc'] and the new [edited']. The result is [Callback
    * stack result acc' edited']. *)
    val m_insert_delta : S_K.primitives
    (**  An example of usage of [m_insert_delta]:
        {[
    let fourth = S_K.create (SK d) (LS a b c d);;
    let seq = S_K.eval (SK ([m_decode_sequence] [fourth] stack results
        accumulator editedsequence K K S K K K K K K S K));;
    let empty_sequence = 
        (SK ([Primitives.m_pair] ([Primitives.m_pair] (K S) K) K))
    let res = S_K.eval (SK ([m_insert_delta] Callback stack results
        [empty_sequence] [seq] K K S K S K S));;]} *)

    (** [m_delete_delta Callback stack results acc edited] reads a logarithmic
    * integer [l] marking the position in the sequence [edited] to be deleted. 
    * The function moves [l - 1] bases in [edited] to [acc] and removes the
    * remaining head of [edited] to
    * produce a new longer [acc'] and the new [edited']. The result is [Callback
    * stack result acc' edited']. *)
    val m_delete_delta : S_K.primitives
    (**  An example of usage of [m_delete_delta]:
        {[
    let fourth = S_K.create (SK d) (LS a b c d);;
    let seq = S_K.eval (SK ([m_decode_sequence] [fourth] stack results
        accumulator editedsequence K K S K K K K K K S K));;
    let empty_sequence = 
        (SK ([Primitives.m_pair] ([Primitives.m_pair] (K S) K) K))
    let res = S_K.eval (SK ([m_delete_delta] Callback stack results
        [empty_sequence] [seq] K K S K S));;]} *)

    (** [m_substi_delta Callback stack results acc edited] reads a logarithmic
    * integer [l] and a base [b] to substitute the corresponding base in
    * position [l] of [edited] with [b]. The
    * function moves [l - 1] bases in [edited] to [acc], adds [b] to [acc], and
    * removes the remaining head in [edited] to
    * produce a new longer [acc'] and the new [edited']. The result is [Callback
    * stack result acc' edited']. *)
    val m_substi_delta : S_K.primitives
    (** The usage of the function is the same as [m_insert_delta]. *)

    (** [make_encoder tree replacement] reads a standard SK expression
    * representing a tree with labels as leaves, and produces an encoding tree
    * suitable for {!Kolmo.IndPM.m_decode_function} where the leaves are the corresponding
    * primitives for the leaves stored in the [replacement] association list. *)
    val make_encoder : 
        S_K.primitives -> (string * S_K.primitives) list -> S_K.primitives
    (** An example of usage of [make_encoder]:
        {[
        let a = S_K.create (SK (a)) (LS a b c);;
        let enc = make_encoder (SK (A B)) [("A", a)];;]}
        The resulting encoder holds a leaf with the function [a] represented by
        [K] and another leaf with the label "B" represented by [S]. *)

    (** [m_decode_function o_encoder encoder stack results acc seq] decodes the
    * incoming stream of bits, decoding the contents of [encoder]. If the
    * [encoder] reaches a leaf [l], the result is the application of the
    * function in [l] as in [l (m_decode_function o_encoder o_encoder) stack
    * results acc seq]. The encoder contains functions to be called
    * sequentially, upon termination, each function is expected to call
    * [(m_decode_function o_encoder o_encoder)] together with the current
    * contents of [stack], [results], [acc], and [seq]. Some functions can,
    * however, continue a different path of execution, for example, when
    * reaching the end of a computation, as in the end of a tree with
    * {!Kolmo.IndPM.m_end_branch}. *)
    val m_decode_function : S_K.primitives
    (** An example of usage is:
        * {[
    let [a; b; c] = get_lst res;;
    let a = S_K.create (SK (a)) (LS a b c);;
    let enc = make_encoder (SK (A B)) [("A", a)];;
    let remove_first = S_K.create (SK e) (LS a b c d e);;
    S_K.to_string (S_K.eval (SK ([make_list] A B C D E)));;
    let enc = make_encoder (SK (A B)) ["A", remove_first];;
    let res = S_K.eval (SK ([m_decode_function] [enc] [enc] stack results acc
    seq K));;
    S_K.to_string res;;]}. *)

    (** [m_end_branch Callback stack results acc seq n1 n2] reads the [n1] and
    * [n2] elements and decides whether or not it is necessary to continue
    * processing more branches in a tree or we have reached a leaf in a tree
    * ([n2] is necessary only if [n1 = S]. Before continuing, the function
    * copies the rest of [seq] in [acc] ([seq] and [acc] are sequences) to
    * generate [acc'], and inverts [acc'] to produce [acc'']. 
    * If [n1 = K] then the branch is internal, and the result of the function is
    * [Callback (acc '' :: stack) results empty_sequence acc''
    * n2] otherwise, if [n2 = S] we have reached the end of a leaf and the
    * result is [Callback stack (acc'' :: results) empty_sequence
    * empty_sequence], but if [n2 = K], then we have finished processing a tree,
    * and the result [acc'' :: results] is the output. This marks the end of
    * processing a phylogenetic tree. *)
    val m_end_branch : S_K.primitives

    (** {2 Using IndPM, an example}
    {[
    open Kolmo;;
    open Primitives;;
    open IndPM;;
    let encoded = 
        make_encoder 
        (SK ((ins del) (subs (endbranch decodesequence))))
        ["ins", m_insert_delta; "del", m_delete_delta; "subs", m_substi_delta; 
            "endbranch", m_end_branch; "decodesequence", m_decode_sequence];;
    let empty_sequence = (SK ([m_pair] ([m_pair] (K S) K) K)) ;;
    let res = S_K.eval (SK
        ([m_decode_function] [encoded] [encoded] A results [empty_sequence]
        [empty_sequence] S S S K K S K K S K K K S S K K K S K S S S S K S
        K));;
    let res = S_K.eval (SK ([res] [m_head]));;
    PM.Dna.to_list res;;]}
    *)

end

module Align : sig
    type matrix = {
        event_cost : float;
        first_event_cost : float;
        matrix : Cost_matrix.Two_D.m;
    }

    val align : Sequence.s -> Sequence.s -> matrix -> float * Sequence.s *
    Sequence.s * Sequence.s
end
