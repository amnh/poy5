type taxon = string

module type S = sig

    val of_channel : FileContents.t -> in_channel -> (Sequence.s list list list * string) list
    (** [of_channel x y] takes an input channel y with information in ascii
    * format of a molecular data of type t and outputs the list of sequences and
    * taxon names contained in the file. If the function finds an unexpected
    * character or an illegal condition for a specific file format, an
    * Illegal_molecular_format or Unsupported_file_format exception is raised. *)

    val to_channel : 
        out_channel -> (Sequence.s * string) list -> Alphabet.a -> unit
    (** [to_channel x y z] writes in the output channel x the sequences and taxon
    * list y using the alphabet z for the sequences in y. If the function finds
    * an illegal element in any of the sequences for the alphabet z, raises an
    * Alphabet.Illegal_Code exception. There is no guarantee on the state of the
    * output file if the exception is raised.*)

    val of_file : FileContents.t -> FileStream.f -> (Sequence.s list list list * string) list
end

type fl = {
    filename : string;
    taxon : string;
    sequence : string;
    character : string;
    line : int;
}

exception Illegal_molecular_format of fl


include S
(** A parser implementation for the Fasta file format *)

val of_string : FileContents.t -> string -> (Sequence.s list list list * string) list 

val process_sequence : 
    bool -> (char Stream.t -> int list -> int -> int list * int) ->
        Alphabet.a -> string list -> Sequence.s

