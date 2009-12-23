(** A bit representation *)
type bit = Zero | One

(** Encoding a natural number *)
type encoding = bit list 

(** The representation of a natural number *)
type natural

(** Type Conversion *)
val to_nat : int -> natural
val of_nat : natural -> int

(** The length of the encoding of a natural number *)
val l : natural -> natural

(** THe binary encoding of a natural number *)
val binary : natural -> encoding

(** Some prefix free encodings *)
val e : natural -> natural -> encoding
val e_2 : natural -> encoding
val hat : natural -> encoding

(** Decoding the previous prefix free codifications *)
val decode : encoding -> natural
(* A function to generate the huffman code functions for lists of elements
* of type 'a *)

val huffman : ('a * float) list -> 
    ((encoding -> 'a list) * ('a list -> encoding))

(* The tree representation of a huffman code *)
val huffman_tree : ('a * float) list -> 'a option Tree.Parse.t

val geometric : int -> int -> float -> (int * float) list
