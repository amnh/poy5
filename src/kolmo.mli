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
    val encoder_example : S_K.primitives
    val decoder : S_K.primitives
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
end
