(** Type of data contained in a file (either Fasta, POY or Hennig86
 * formats *)
type t = 
    | Nucleic_Acids 
    | Proteins 
    | Prealigned_Alphabet of Alphabet.a
    | AlphSeq of Alphabet.a
    | Inactive_Character 
    | Ordered_Character of (int * int * bool)
    | Unordered_Character of (int * bool)
    | Sankoff_Character of (int list * bool)
    | Genes of int array


