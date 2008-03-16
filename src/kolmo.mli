module type E =
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
  end
module Encodings : E

module S_K :
  sig
    exception Illegal_Expression of string Parser.Tree.t list
    type primitives = [`S | `K | `Label of string | `Node of primitives list ]
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

module PM : sig
    val encode_in_log : int -> S_K.primitives
    module type List = sig
        type ocaml_repr
        type ocaml_list
        val empty_list : S_K.primitives
        val is_not_empty : S_K.primitives -> S_K.primitives
        val ml_zero_signal : S_K.primitives -> bool
        val head : S_K.primitives -> S_K.primitives
        val tail : S_K.primitives -> S_K.primitives
        val prepend : S_K.primitives -> S_K.primitives -> S_K.primitives
        val alphabet : (S_K.primitives * ocaml_repr) list
        val to_ocaml : S_K.primitives -> ocaml_repr
        val of_ocaml : ocaml_repr -> S_K.primitives
        val to_list : S_K.primitives -> ocaml_list
        val of_list : ocaml_list -> S_K.primitives
        val complement : S_K.primitives -> S_K.primitives
    end

    module type Integer = sig
        val zero : S_K.primitives
        val not_zero : S_K.primitives -> S_K.primitives
        val ml_zero_signal : S_K.primitives -> bool
        val predecessor : S_K.primitives -> S_K.primitives
        val successor : S_K.primitives -> S_K.primitives
        val to_int : S_K.primitives -> int
        val of_int : int -> S_K.primitives
    end

    module ChurchIntegers : Integer

    module Dna : List with type ocaml_repr = string with type ocaml_list =
        string list

    module Chromosome : List with type ocaml_repr = Dna.ocaml_list with type
    ocaml_list = Dna.ocaml_list list

    module type SE = sig 
        val insert : S_K.primitives
        val delete : S_K.primitives
        val substitute : S_K.primitives
        type insert_generator 
        type delete_generator 
        type substitute_generator
        val insert_generator : insert_generator
        val delete_generator : delete_generator
        val substitute_generator : substitute_generator
    end

    module AtomicSE : functor (I : Integer) ->
        functor (S : List) -> SE 
        with type insert_generator = int -> S.ocaml_repr -> S_K.primitives
        with type delete_generator = int -> S_K.primitives
        with type substitute_generator = int -> S.ocaml_repr -> S_K.primitives

    module AffineSE : functor (I : Integer) ->
        functor (S : List) -> SE
        with type insert_generator = int -> S_K.primitives
        with type delete_generator = int -> int -> S_K.primitives
        with type substitute_generator = int -> S_K.primitives

    module AtomicSE_LogInts : functor (S : SE) -> SE
        with type insert_generator = int -> S_K.primitives
        with type delete_generator = int -> S_K.primitives
        with type substitute_generator = int -> S_K.primitives

    module AffineSE_LogInts : functor (S : SE) -> SE
        with type insert_generator = int -> S_K.primitives
        with type delete_generator = int -> int -> S_K.primitives
        with type substitute_generator = int -> S_K.primitives 

    module type HO = sig
        include SE
        val translocate : S_K.primitives
        val invert : S_K.primitives
        val duplicate : S_K.primitives
        val tandem : S_K.primitives
    end

    module HighOrder : functor (I : Integer) -> 
        functor (S : List) -> HO

    (*
    val encoder_example : S_K.primitives
    val decoder : bool -> S_K.primitives
    val result_example : S_K.primitives
    val make_encoder : 
        S_K.primitives -> (string * S_K.primitives) list -> S_K.primitives
    val f_apply : S_K.primitives
    val f_edit : S_K.primitives
    val pm_encodings : S_K.primitives
    val encode_in_log : int -> S_K.primitives

    val pm_encodings : S_K.primitives
    val substitute : S_K.primitives
    val insert : S_K.primitives
    val delete : S_K.primitives
    val adenine : S_K.primitives
    val citosine : S_K.primitives
    val guanine : S_K.primitives
    val timine : S_K.primitives
    val apply_list_of_functions_on_accumulator : S_K.primitives
    val insert_composable : S_K.primitives
    val delete_composable : S_K.primitives
    val substitute_composable : S_K.primitives
    val create_sequence : string -> S_K.primitives
    val generate_alignment_machine : string -> string -> S_K.primitives
    val generate_acc_and_edition_machine : 
        string -> string -> S_K.primitives * S_K.primitives
    val process_tree : S_K.primitives
    *)
end
