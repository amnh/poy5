(* POY 4.0 Beta. A phylogenetic analysis program using Dynamic Homologies.    *)
(* Copyright (C) 2007  Andrés Varón, Le Sy Vinh, Illya Bomash, Ward Wheeler,  *)
(* and the American Museum of Natural History.                                *)
(*                                                                            *)
(* This program is free software; you can redistribute it and/or modify       *)
(* it under the terms of the GNU General Public License as published by       *)
(* the Free Software Foundation; either version 2 of the License, or          *)
(* (at your option) any later version.                                        *)
(*                                                                            *)
(* This program is distributed in the hope that it will be useful,            *)
(* but WITHOUT ANY WARRANTY; without even the implied warranty of             *)
(* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the              *)
(* GNU General Public License for more details.                               *)
(*                                                                            *)
(* You should have received a copy of the GNU General Public License          *)
(* along with this program; if not, write to the Free Software                *)
(* Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301   *)
(* USA                                                                        *)

(** Parsing files for phylogenetic analysis. *)

(** A file location specification. This type is used in the error reporting for
 * molecular files, such as Fasta and POY formats. *)
type fl = {
    filename : string;  (** The filename of the error *)
    taxon : string;     (** Name of the taxon containing the error *)
    sequence : string;  (** Sequence where the error was detected *)
    character : string; (** Character where the error was detected *)
    line : int;         (** Line number of the error occurrence *)
}
    
(** Illegal format for Hennig86 files. The starting position of the error is
* reported in the string, but note that the string will contain the complete
* file starting there. *)
exception Illegal_hennig86_format of string

(** Illegal format for DNA and Protein sequences. The error location is reported
 * in the fl dataype. *)
exception Illegal_molecular_format of fl

(** In the current version of POY only DNA sequences are supported. In the near
 * future, Protein sequences should supported, but in the meanwhile any attempt
 * to analyze this type will raise this exception. If necessary, other file
 * formats that should but are not supported yet can raise this exception. *)
exception Unsupported_file_format of string

(** The file type as detected or according to the user suplied format. *)
type ft = 
    | Is_Hennig     (** A Hennig86 file format *)
    | Is_Dpread     (** A POY dpread file format *)
    | Is_Clustal    (** A Clustal file format *)
    | Is_Fasta      (** A Fasta file format *)
    | Is_Poy        (** A Poy file format *)
    | Is_Genome     (** A Genome (Grappa) file format *)
    | Is_ASN1       (** An ASN1 file format *)
    | Is_Genbank    (** A Genbank file format *)
    | Is_INSDSeq    (** A INSDSeq XML file format *)
    | Is_GBSeq      (** A GBSeq XML file format *)
    | Is_TinySeq    (** A TinySeq XML file format *)
    | Is_XML        (** A XML file format *)
    | Is_Nexus      (** A Nexus file format *)
    | Is_NewSeq     
    | Is_Dictionary
    | Is_Fixed_States_Dictionary
    | Is_Phylip     (** PHYLIP file format *)     
    | Is_Unknown    (** Unknown file format *)
    | Is_Trees      (** Trees file in parenthetical format *)
    | Is_ComplexTerminals (** Define complex terminals *)

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

(** [is_unknown t] tells observers whether this present character was listed as
    unknown by the user *)
val is_unknown : t -> bool

(** A taxon name ... *)
type taxon = string

type filename = [ `Local of string | `Remote of string ]

(** [test_file a] tries to guess the appropriate format of a file with path [a].
* If the file doesn't exists raises a Sys_error of string exception, otherwise it
 * will return the guessed type (or Is_Unknown) *)
val test_file : filename -> ft


(** [lor_list_withhash l hash] returns the logical or of the hash values of all
    of the values in [l] *)
val lor_list_withhash : int list -> (int, int) Hashtbl.t option -> int


module Tree : sig
    (** Parser for tree in format (a (b c)) *)

    type 'a t = Leaf of 'a | Node of 'a t list * 'a

    (** [of_string x] given a string x of a tree in the form (a (b c)), returns
    * its representation as an internal t type. If an error occurs, an
    * Illegal_tree_format error is raised. *)
    val of_string : string -> string t list list 

    (** [of_channel x] creates a list of all the trees contained in the input
    * channel x in ascii format. Each tree should be in a single line. In case
    * of error the function raises Illegal_molecular_format. *)
    val of_channel : in_channel -> string t list list 

    (** [of_file] is a shortcut of [of_channel], when the channel is an opened
    * file. *)
    val of_file : filename ->  string t list list

    (** [of_string_annotated] is simmilar to [of_string] excepting that the
    * returned list includes the associated information contained in square
    * brackets together with each tree, instead of ignoring it as [of_string]
    * does. *)
    val of_string_annotated : string -> (string t * string) list list

    (** [of_channel_annotated] is to [of_channel] as [of_string_annotated] is to
    * [of_string]. *)
    val of_channel_annotated : in_channel -> (string t * string) list list

    (** [of_file_annotated] is to [of_file] as [of_string_annotated] is to
    * [of_string]. *)
    val of_file_annotated : filename -> (string t * string) list list

    (** [stream_of_file f] produces a function that returns on each call one of
    * the trees in the input file [f]. If no more trees are found, an
    * End_of_file exception is raised. The function _requires_ that the trees be
    * separated with associated information, semicolons, or stars.*)
    val stream_of_file : filename -> (unit -> (string t * string))
    val cannonic_order : string t -> string t

end


(** Hennig file format parser *)
module OldHennig : sig

    (** Encoding specification handler. *)
    module Encoding : sig 
        (** An encoding specification *)
        type s 

        (** Creates the default specification *)
        val default : 'a -> s

        (** Get the minimum value observed in the states of a character. *)
        val get_min : s -> int

        (** Get the maximum value observed in the states of a character. *)
        val get_max : s -> int 

        (** Set the minimum value observed in the states of a character. *)
        val set_min : s -> int -> s

        (** Set the maximum value observed in the states of a character. *)
        val set_max : s -> int -> s

        (** Get the set of observed states in a character. *)
        val get_set : s -> All_sets.Integers.t

        val set_set : s -> All_sets.Integers.t -> s

        val set_tcm : s -> int array array -> s

        (** Get the assigned weight to a character *)
        val get_weight : s -> int

        val set_weight : s -> int -> s

        (** True of the character is active, otherwise false *)
        val is_active : s -> bool

        (** True if the character is ordered (additive), false otherwise (non
        * additive) *)
        val is_ordered : s -> bool
        
        (** True if the character is sankoff *)
        val is_sankoff : s -> bool
        
        val get_observed_used : s -> (int, int) Hashtbl.t option

        val get_used_observed : s -> (int, int) Hashtbl.t option

        val to_string : s -> string
         
        val has_states : int -> int -> s -> bool

        val get_tcm : s -> int array array

        val dna_encoding : s

        val gap_encoding : int -> s

        val set_unordered : s -> s

        val set_sankoff : s -> int array array -> s
    end

    (** [of_channel x] takes as input a channel_in of text containing the
    * Hennig86 file to be read, and parses it, to return a tuple [(a, b)], where
    * [a] is the encoding specifications of the contents of the character array 
    * (the mapping used to store the information in the character array so it 
    * can be decoded later), and  [b] is the list of all taxa. Each taxon in 
    * this list is a tuple (c, d), where [c] is the array of its characters in 
    * the Hennig86 file, and [d] is the name of the taxon on it. *)
    val of_channel : 
        in_channel -> Encoding.s array * (t array * string) list * string Tree.t
        list list

    val of_file :
        filename -> Encoding.s array * (t array * string) list * string Tree.t
        list list
    val convert_to_Phylip_format_file : filename -> string -> unit

    (** Splits the set of characters of a taxa in additive and non additive in a
    * tuple. *)
    val split_ordered : Encoding.s array -> (t array * string) list -> 
        (Encoding.s array * (t array * string) list) * 
        (Encoding.s array * (t array * string) list)

    val merger : (Encoding.s array * (t array * string) list) list -> 
        (int * int) array * (t array * int) list 

    val print_character_specs : Pervasives.out_channel -> unit
    
    (** [character_states_minmax min max] returns the character specification
    * codes  that have the state codes between [min] and [max] inclusive. *)
    val character_states_minmax : int -> int -> int list

    (** [character_states_number min max] returns the character specification
    * codes that have cardinality of the states set between [min] and [max]
    * inclusive. *)
    val character_states_number : int -> int -> int list

    (** [character_additive ()] returns the list of character specification
    * codes that have additive specification. *)
    val character_additive : unit -> int list

    (** [character_nonadditive ()] returns the list of character specification
    * codes that have nonadditive specification. *)
    val character_nonadditive : unit -> int list

    (** [character_active ()] returns the list of character specification
    * codes that are active. *)
    val character_active : unit -> int list

    (** [character_inactive ()] returns the list of character specification
    * codes that are inactive. *)
    val character_inactive : unit -> int list

    (** [character_complement lst] returns the list of character specification
    * codes that are not in the list of character specification codes [lst]. *)
    val character_complement : int list -> int list 
    
    (** [character_additive_minmax min max] returns the character specification
    * codes that are additive and have the state codes between [min] and [max] 
    * inclusive. *)
    val character_additive_minmax : int -> int -> int list
   
    (** [character_nonadditive_minmax min max] returns the character 
    * specification codes that are nonadditive and have the state codes between 
    * [min] and [max] inclusive. *)
    val character_nonadditive_minmax : int -> int -> int list
    
    (** [character_active_minmax min max] returns the character specification
    * codes that are active and have the state codes between [min] and [max] 
    * inclusive. *)
    val character_active_minmax : int -> int -> int list
    
    (** [character_inactive_minmax min max] returns the character specification
    * codes that are inactive and have the state codes between [min] and [max] 
    * inclusive. *)
    val character_inactive_minmax : int -> int -> int list

    (** [general_character_filter f] returns the list of character specification
    * codes with specification [s] such that [f s] is [true]. *)
    val general_character_filter : (Encoding.s -> bool) -> int list

    (** [filter_matrix l m] selects those characters with character specification
    * code in the list [l] from the character matrix [m]. *)
    val filter_matrix : int list -> (int * int) array * (t array * int) list ->
        (int * int) array * (t array * int) list 

    val clear_characters : unit -> unit

    val clear_taxa : unit -> unit

    val flatten_matrix : (int * int) array * (t array * int) list -> 
        ((t * int) array * int) list

    val code_character : int -> Encoding.s

    val character_code : Encoding.s -> int * int

    val code_character_hom : int -> Encoding.s
    
    val code_taxon : int -> string

    val taxon_code : string -> int

    (** [categorize_chars (array, taxon_code)] returns a 5-tuple of [(sub_array,
     * taxon_code)] values, where the sub-arrays are, in order:
     * - additive characters
     * - nonadditive characters with <= 8 states
     * - nonadditive chars with 9 <= states <= 16
     * - nonadditive chars with 17 <= states <= 32
     * - nonadditive chars with >=33 states
     *)
    val categorize_chars : ((t * int) array * int) ->
        (((t * int) array * int) * ((t * int) array * int) *
             ((t * int) array * int) * ((t * int) array * int) *
             ((t * int) array * int))

end

module type MOLECULAR = sig

    val of_channel : t -> in_channel -> (Sequence.s list list list * taxon) list
    (** [of_channel x y] takes an input channel y with information in ascii
    * format of a molecular data of type t and outputs the list of sequences and
    * taxon names contained in the file. If the function finds an unexpected
    * character or an illegal condition for a specific file format, an
    * Illegal_molecular_format or Unsupported_file_format exception is raised. *)

    val to_channel : 
        out_channel -> (Sequence.s * taxon) list -> Alphabet.a -> unit
    (** [to_channel x y z] writes in the output channel x the sequences and taxon
    * list y using the alphabet z for the sequences in y. If the function finds
    * an illegal element in any of the sequences for the alphabet z, raises an
    * Alphabet.Illegal_Code exception. There is no guarantee on the state of the
    * output file if the exception is raised.*)

    val of_file : t -> filename -> (Sequence.s list list list * taxon) list
end

module Poy : MOLECULAR
(** A parser implementation for the Poy file format *)

val alphabet_of_t : t -> Alphabet.a
(** [alphabet_of_t t] returns the alphabet associated with this character
    type *)

module Fasta : MOLECULAR
(** A parser implementation for the Fasta file format *)

module ASN1 : sig
(** A parser implementation for the ASN.1 file format - just convert to 
* Fasta Format*)
   val convert_to_fasta: filename -> in_channel
   
end 

module Genbank : sig
(** A parser implementation for the Genbank file format - just convert to 
* Fasta Format*)
    val convert_to_fasta: ?filename:string -> filename -> in_channel
   
end 


module INSDSeq : sig
(** A parser implementation for the INSDSeq XML file format - just convert to 
* Fasta Format*)
   val convert_to_fasta: filename -> in_channel
   
end 

module GBSeq : sig
(** A parser implementation for the GBSeq XML file format - just convert to 
* Fasta Format*)
   val convert_to_fasta: filename -> in_channel
   
end  

module TinySeq : sig
(** A parser implementation for the TinySeq XML file format - just convert to 
* Fasta Format*)
   val convert_to_fasta: filename -> in_channel
   
end  

module XML : sig
(** A parser implementation for the XML file format - just convert to 
* Fasta Format*)
   val convert_to_fasta: filename -> in_channel
   
end   

module NewSeq : sig
(** A parser implementation for the file format that looks like the following
*   Coleosporium_asterum_AB011052
*      1     AAAGAUUAAG CCAUGCAUGU CUAAGUAUAA ACAAUUAUAC AGUGAAACUG CGAAUGGCUC
*     61     AUUAAAUCAG UUAUAGUUUA UUUGAUGAAA CCUUACUACA UGGGAUAACU GUGGUAAUUC
*)
   val convert_to_fasta: filename -> in_channel
   val to_fasta : filename -> string
   
end   

module Phylip : sig
(** A parser implementation for the Phylip file format - just convert to 
* Hennig Format*)
   val convert_to_hennig: filename -> in_channel
   
end 

module GrappaParser : sig
    val of_channel : in_channel -> t array 

end

module TransformationCostMatrix : sig

    val of_channel : 
        ?orientation:bool -> ?use_comb:bool -> FileStream.greader -> Cost_matrix.Two_D.m

    val of_channel_nocomb: ?orientation:bool -> FileStream.greader -> Cost_matrix.Two_D.m

    val of_list : ?use_comb:bool -> int list list -> Cost_matrix.Two_D.m 

    val of_file : ?use_comb:bool -> filename -> Cost_matrix.Two_D.m

end

module Dictionary : sig
    val of_channel : in_channel -> (string, string) Hashtbl.t
end

module FixedStatesDict : sig

    val of_channel : t -> in_channel -> Sequence.s list

    val create_dp_read_file : 
        ?filename:string -> string -> (Sequence.s * string) list -> 
            Sequence.s list -> Cost_matrix.Two_D.m -> string

end

module IgnoreList : sig
    val of_channel : in_channel -> string list
end

val print_error_message : fl -> unit

val molecular_to_fasta : filename -> in_channel

module SC : sig
    type st_type = 
        | STOrdered
        | STUnordered  
        | STSankoff of int array array (* If Sankoff, the cost matrix to use *)


    type static_spec = {
        st_filesource : string; (* The file that contained the character
        originally *)
        st_name : string;       (* The name assigned to the character *)
        st_alph : Alphabet.a;   (* The set of potential character symbols *)
        st_observed : int list; (* The set of observed states *)
        st_labels : string list;(* The labels assigned to the states *)
        st_weight : float;      (* The character weight *)
        st_type : st_type;      (* The type of character *)
        st_equivalents : (string * string list) list;
                                (* Things that are the same in the input *)
        st_missing : string;       (* The character that represents missing data *)
        st_matchstate : string option; 
            (* The chaaracter that marks the same state as teh first taxon *)
        st_gap : string;            (* The gap representation *)
        st_eliminate : bool;       (* Wether or not the user wants to get rid of it *)
        st_case : bool;       (* Wether or not the user wants be case sensistive *)
        st_used_observed : (int, int) Hashtbl.t option;
        st_observed_used : (int, int) Hashtbl.t option;
    }

    type static_state = int list option

    (** [spec_of_alph alphabet missing gap] generates a specification that can read
    * the elements in the [alphabet] using when the matrix represents [missing]
    * data and [gaps] as specified. *)
    val spec_of_alph : string list -> string -> string -> string -> static_spec

    (** [to_string v] outputs a string representation of the static_specification
    * [v] *)
    val to_string : static_spec -> string

    (** [to_formatter v] generates a standard [Tags.output] representation of
    * [to_formatter]. *)
    val to_formatter : static_spec -> Tags.output

    (** [file_output] is the standard output generated by a Static Homology file
    * parser *)
    type file_output = 
    (string option array * static_spec array * static_state array array * 
    string Tree.t list list * 
    ((Alphabet.a * (Sequence.s list list list * taxon) list) list))

    (** [of_channel style ch filename] generates a parsed output from the input
    * file contained in the [ch] channel, using the appropriate [style] to parse
    * either Hennig or Nexus files. The characters are assigned [st_filesource]
    * [filename]. *)
    val of_channel : [`Hennig | `Nexus ] -> in_channel -> string -> file_output

    (** [fill_observed parsed] takes a parsed output and updates whatever
    * internal information is needed in a post parsing step (for now only the observed
    * states in the [static_spec] field) *)
    val fill_observed : file_output -> unit

    (** [of_old_spec filename alph old position] converts an old static homology
    * specification to the new style. *)
    val of_old_spec : 
        string -> Alphabet.a option -> OldHennig.Encoding.s -> int -> 
            static_spec

    (** Symmetric to the previous one, but for the observed state of a taxon *)
    val of_old_atom : static_spec -> OldHennig.Encoding.s -> t -> static_state

    (** [of_old_parser filename alphabets old_parsed] converts the [old_parsed]
     * style of static homology parsed file to the new style, with (optional)
     * assignment of [alphabets], and is assigned [st_filesource] [filename]. *)
    val of_old_parser : 
        string ->
        Alphabet.a array option ->
        OldHennig.Encoding.s array * (t array * string) list * string Tree.t list list ->
            file_output
end

module PAlphabet : sig
    val of_file : filename -> bool -> bool -> Alphabet.a * Cost_matrix.Two_D.m *
        Cost_matrix.Three_D.m
    end

module SetGroups : sig
    type set_type =
            | Group 
            | Any of float
            | Tree of float             (* origin/loss cost *)
    type 'a ct =
            | Elt of 'a
            | Set of 'a * set_type * 'a ct list
    type t = string ct

    val of_file : filename -> t list
    val to_string : t -> string
    val unify : t -> t -> t             (* can fail *)
    val unify_list : t list -> t
    val coerce_to_type : t -> t -> t
    end


