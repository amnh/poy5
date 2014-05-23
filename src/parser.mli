(* POY 5.1.1. A phylogenetic analysis program using Dynamic Homologies.       *)
(* Copyright (C) 2014 Andrés Varón, Lin Hong, Nicholas Lucaroni, Ward Wheeler,*)
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

module E : sig
        
    (** Illegal format for Hennig86 files. The starting position of the error is
    * reported in the string, but note that the string will contain the complete
    * file starting there. *)
    exception Illegal_hennig86_format of string

    (** In the current version of POY only DNA sequences are supported. In the near
     * future, Protein sequences should supported, but in the meanwhile any attempt
     * to analyze this type will raise this exception. If necessary, other file
     * formats that should but are not supported yet can raise this exception. *)
    exception Unsupported_file_format of string

    (** A taxon name ... *)
    type taxon = string

end

module Files : sig
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

    (** [is_unknown t] tells observers whether this present character was listed as
        unknown by the user *)
    val is_unknown : FileContents.t -> bool

    (** [test_file a] tries to guess the appropriate format of a file with path [a].
    * If the file doesn't exists raises a Sys_error of string exception, otherwise it
     * will return the guessed type (or Is_Unknown) *)
    val test_file : FileStream.f -> ft

    val molecular_to_fasta : FileStream.f -> in_channel

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

        val set_likelihood_model : MlModel.model option -> s -> s

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

        val rearr_encoding : int -> s

        val set_unordered : s -> s

        val set_sankoff : s -> int array array -> s

        val print : s -> unit 
    end

    (** [of_channel x] takes as input a channel_in of text containing the
    * Hennig86 file to be read, and parses it, to return a tuple [(a, b)], where
    * [a] is the encoding specifications of the contents of the character array 
    * (the mapping used to store the information in the character array so it 
    * can be decoded later), and  [b] is the list of all taxa. Each taxon in 
    * this list is a tuple (c, d), where [c] is the array of its characters in 
    * the Hennig86 file, and [d] is the name of the taxon on it. *)
    val of_channel : 
        in_channel -> Encoding.s array * (FileContents.t array * string) list *
            (string option * Tree.Parse.tree_types list) list

    val of_file :
        FileStream.f -> Encoding.s array * (FileContents.t array * string) list * 
            (string option * Tree.Parse.tree_types list) list

    (** Splits the set of characters of a taxa in additive and non additive in a
    * tuple. *)
    val split_ordered : Encoding.s array -> (FileContents.t array * string) list -> 
        (Encoding.s array * (FileContents.t array * string) list) * 
        (Encoding.s array * (FileContents.t array * string) list)

    val merger : (Encoding.s array * (FileContents.t array * string) list) list -> 
        (int * int) array * (FileContents.t array * int) list 

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
    val filter_matrix : int list -> (int * int) array * (FileContents.t array * int) list ->
        (int * int) array * (FileContents.t array * int) list 

    val clear_characters : unit -> unit

    val clear_taxa : unit -> unit

    val flatten_matrix : (int * int) array * (FileContents.t array * int) list -> 
        ((FileContents.t * int) array * int) list

    val code_character : int -> Encoding.s

    val process_matrix : 
        string -> bool -> int -> int -> (string * int list list) list

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
     * - nonadditive chars with >=33 states *)
    val categorize_chars : ((FileContents.t * int) array * int) ->
        (((FileContents.t * int) array * int) * ((FileContents.t * int) array * int) *
             ((FileContents.t * int) array * int) * ((FileContents.t * int) array * int) *
             ((FileContents.t * int) array * int))

    val generate_alphabet : Alphabet.a option -> Encoding.s -> Alphabet.a

    (** [of_old_spec filename alph old position] converts an old static homology
    * specification to the new style. *)
    val to_new_spec : 
        ?separator:string ->
        string -> Alphabet.a -> Encoding.s -> int -> 
           Nexus.File.static_spec

    (** Symmetric to the previous one, but for the observed state of a taxon *)
    val to_new_atom :
        (FileContents.t, Nexus.File.static_state) Hashtbl.t ->
            Nexus.File.static_spec -> Alphabet.a option ->
                FileContents.t -> Nexus.File.static_state

    (** [of_old_parser filename alphabets old_parsed] converts the [old_parsed]
     * style of static homology parsed file to the new style, with (optional)
     * assignment of [alphabets], and is assigned [st_filesource] [filename]. *)
    val to_new_parser : 
        ?separator:string ->
        string ->
        Alphabet.a array option ->
        Encoding.s array * (FileContents.t array * string) list * 
        (string option * Tree.Parse.tree_types list) list -> Nexus.File.nexus
end

module Asn1 : sig
(** A parser implementation for the ASN.1 file format - just convert to 
* Fasta Format*)
   val convert_to_fasta: FileStream.f -> in_channel
   
end 

module Genbank : sig
(** A parser implementation for the Genbank file format - just convert to 
* Fasta Format*)
    val convert_to_fasta: ?filename:string -> FileStream.f -> in_channel
   
end 

module INSDSeq : sig
(** A parser implementation for the INSDSeq XML file format - just convert to 
* Fasta Format*)
   val convert_to_fasta: FileStream.f -> in_channel
   
end 

module GBSeq : sig
(** A parser implementation for the GBSeq XML file format - just convert to 
* Fasta Format*)
   val convert_to_fasta: FileStream.f -> in_channel
   
end  

module TinySeq : sig
(** A parser implementation for the TinySeq XML file format - just convert to 
* Fasta Format*)
   val convert_to_fasta: FileStream.f -> in_channel
   
end  

module XMLGB : sig
(** A parser implementation for the XML file format - just convert to 
* Fasta Format*)
   val convert_to_fasta: FileStream.f -> in_channel
   
end   

module NewSeq : sig
(** A parser implementation for the file format that looks like the following
*   Coleosporium_asterum_AB011052
*      1     AAAGAUUAAG CCAUGCAUGU CUAAGUAUAA ACAAUUAUAC AGUGAAACUG CGAAUGGCUC
*     61     AUUAAAUCAG UUAUAGUUUA UUUGAUGAAA CCUUACUACA UGGGAUAACU GUGGUAAUUC
*)
   val convert_to_fasta: FileStream.f -> in_channel
   val to_fasta : FileStream.f -> string
   
end   

module Dictionary : sig
    val of_channel : in_channel -> (string, string) Hashtbl.t
    val of_channel_assoc : in_channel -> (string * string) list
end

module FixedStatesDict : sig

    val of_channel : FileContents.t -> in_channel -> Sequence.s list

    val create_dp_read_file : 
        ?filename:string -> string -> (Sequence.s * string) list -> 
            Sequence.s list -> Cost_matrix.Two_D.m -> string

end

module IgnoreList : sig
    val of_channel : in_channel -> string list
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

    val of_file : FileStream.f -> t list
    val to_string : t -> string
    val unify : t -> t -> t             (* can fail *)
    val unify_list : t list -> t
    val coerce_to_type : t -> t -> t
    end


module Wildcard : sig
    val explode_filenames : FileStream.f list -> string list
end
