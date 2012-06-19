(* POY 5.0 Alpha. A phylogenetic analysis program using Dynamic Homologies.   *)
(* Copyright (C) 2011 Andrés Varón, Lin Hong, Nicholas Lucaroni, Ward Wheeler,*)
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

exception Illegal_argument

(** This type is needed for the dynamic homologies datasets. Contains all the
 * valid options to perform a dynamic homology analysis. *)
type dynhom_opts = 
    | Tcm of string     (** A transformation cost matrix to be used *)

(* The valid contents of an input file. These are not all of them, add others if
* needed, these are all that we need for now. *)
type contents = 
    | Characters  (** Terminal characters *)
    | CostMatrix  (** A transformation cost matrix *)
    | Trees       (** Trees *)

type parsed_trees = ((string option * Tree.Parse.tree_types list) * string * int)

type dyna_state_t = [
|`SeqPrealigned
(** A short sequence, no rearrangements are allowed*)
| `Seq
| `Ml

(** A long sequence, genes inside are broken down
 * automatically, rearrangements are allowed*)
| `Chromosome

| `Genome 

(** A list of shorted sequences 
 * annotated by piles, rearrangements are allowed *)
| `Annotated

(** A sequence of gene names, rearrangements are allowed *)
| `Breakinv 

| `CustomAlphabet]

type polymorphism_t = Methods.polymorphism_arg

(* none of the following should be in this module ... *)
type annotate_tool_t =
    [ `Default of int * int * int | `Mauve of float * float * float * float ]
type median_solver_t =
    [ `Albert
    | `BBTSP
    | `COALESTSP
    | `ChainedLK
    | `MGR
    | `Siepel
    | `SimpleLK
    | `Vinh ]
type re_meth_t = [ `Locus_Breakpoint of int | `Locus_Inversion of int ]
type dyna_pam_t = {
  median_solver : median_solver_t option;
  annotate_tool : annotate_tool_t option;
  re_meth : re_meth_t option;
  circular : int option;
  locus_indel_cost : (int * int) option;
  chrom_indel_cost : (int * int) option;
  chrom_hom : int option;
  translocation : int option;
  keep_median : int option;
  swap_med : int option;
  approx : bool option;
  symmetric : bool option;
  max_3d_len : int option;
  max_kept_wag : int option;
  mode : [ `Chromosome | `Genome | `Breakinv ] option;
}
val dyna_pam_default : dyna_pam_t
(* done *)

val use_mauve_annotator : dyna_pam_t -> bool 

type clip = Clip | NoClip

type tcm_definition = 
    | Substitution_Indel of (int * int)
    | Input_file of (string * (int list list))
    | Substitution_Indel_GapOpening of (int * int * int)
    | Input_file_GapOpening of (string * (int list list) * int)
    | Level of (tcm_definition * int)

val default_tcm : tcm_definition



type dyna_initial_assgn = [ 
    | `Partitioned of clip
    | `AutoPartitioned of (clip * int * (int,  ((int * int) list)) Hashtbl.t)
    | `GeneralNonAdd
    | `DO ]

type dynamic_hom_spec = {
    filename : string;
    fs : string;
    tcm : tcm_definition;
    initial_assignment : dyna_initial_assgn;
    tcm2d : Cost_matrix.Two_D.m;
    tcm3d : Cost_matrix.Three_D.m; 
    lk_model : MlModel.model option;
    alph : Alphabet.a;
    state : dyna_state_t;
    pam : dyna_pam_t;
    weight : float;
    (** choose a way to deal with polymorphism. during transformation to
    * fixed_state, we resolve polymorphism. here are 3 ways of doing it: 
        * 1. Do_All: do the full "get_closest" thing, which might take a long time
        * 2. Pick_One: just pick one.
        * 3. Do_Nothing: do nothing, leave the input sequence as it is.*)
    polymorphism : polymorphism_t;
}

type fixed_state_spec =
    {   costs   : float array array;
        seqs    : Sequence.s array;
        codes   : (int, int) Hashtbl.t;
        opt_bls : float array array option;
        original_dynspec : dynamic_hom_spec;
    }

val get_weight_from_fs_spec : fixed_state_spec -> float

type distr =
    | MaxLength of int
                        (* Any of the distributions with a maximum length *)

type affine_f = {
    selfp : float; (* The natural logarithm of the probability *)
    distr : distr; (* The probability distribution *)
}

type model = 
    | InDels of (float * float)
    | InDelSub of (float * float * float)
    | Subs of float
    | AffInDelSub of (affine_f * affine_f * float)
    | AffInDelAffSub of (affine_f * affine_f * affine_f)

type arr = (float, Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array1.t

type basic_kolmo_spec = {
    event_prob : float;         (** Probability of an event *)
    tm : float list list;     (* The transformation cost matrix *)
    be : arr;         (* Base encodings *)
    simplebe : float array;              (** The original encoding cost *)
    simplebed : float array;            (** The deletion cost in the tips *)
    ins : float; (* Insertion encoding cost *)
    del : float; (* Deletion encoding cost *)
    sub : float; (* Substitution encoding cost *)
    ins_opening : float;        (** Insertion opening *)
    del_opening : float;        (** Deletion opening *)
    sub_opening : float;        (** Substitution opening *)
    root_cost : float;          (** Extra cost incurred by each root *)
    branch_cost : float;        (** Extra cost incurred by each branch *)
    leaf_cost : float;          (** Extra cost incurred by marking a leaf *)
    end_cost : float;           (** Extra cost of ending the compuations and
    producing the output *)
    mo : model; (* The model *)
}

type aux_kolmo_spec = {
    funset : (Methods.indel_prob option) * (Methods.substitution_prob option);
    wordset : int;
    intset : int;
    kolmo_spec : basic_kolmo_spec;
}


type kolmo_spec = {
    dhs : dynamic_hom_spec;
    ks : aux_kolmo_spec;
}

type static_hom_spec =  
    | NexusFile of Nexus.File.static_spec 
    | FixedStates of fixed_state_spec 

(** A character specification. Contains the information regarding the
 * characteristics of a particular character *)
type specs = 
    (** Static homology characters, includes the encoding specs from the
    * parser *)
    | Static of static_hom_spec (* Nexus.File.static_spec or fixed states*)
    (** A dynamic homology based character type, with three parameters, the
    * file name containing the set of sequences, the filename of the valid 
    * fixed states that can be used for that set of sequences, and the file
    * containing the transformation cost matrix to be used to perform their
     * alignments.  Also, we store the alphabet used. *)
    | Dynamic of dynamic_hom_spec
    | Set 
    | Kolmogorov of kolmo_spec


(** [specified] tells us whether or not the user left out the given character
    in the current taxon *)
type specified = [ `Specified | `Unknown ]


type bool_characters = Methods.characters (*[
    | `All
    | `Some of (bool * int list)
    | `Names of (bool * string list)
    | `Random of float
    | `CharSet of (bool * string list)
    | `AllStatic
    | `AllDynamic
    | `Missing of (bool * int)
    | `Range of ( bool * string * int * int)
]*)



type characters = [
    | `All
    | `Some of int list 
    | `Names of string list
    | `Random of float
    | `CharSet of string list
    | `AllStatic
    | `AllDynamic
    | `Missing of (bool * int)
    | `Range of (string * int * int)
]

val transform_range_to_codes : string -> int -> int -> [> `Names of string list]

type 'a seq_t = {
    seq : 'a;
    delimiter : int list;
    code : int;
}

type 'a dyna_data = {
    (** typically 'a is a sequence, the integer followed is the code of 'a *)
    seq_arr : ('a seq_t) array;
}

type 'a fixedstates_data = {
        states : int list;
            dynamic_data : 'a dyna_data;
}

(* A valid loaded character from a file *)
type cs_d = 
    (* A dynamic homology, containing its code, the sequence, the
    * transformation cost matrix and its three dimensional transformation
    * cost matrix *)
    |  Dyna of (int * Sequence.s dyna_data)
    (* fixed states , is static, but still carry everything from dynamic side*)
    | FS of int 
    (* A static homology character, containing its code, and the character
    * itself. If None, means missing data. *)
    | Stat of (int * Nexus.File.static_state)

type cs = cs_d * specified

(** A transformation cost matrix for Sankoff characters *)
type tcm = int array array

module OutputInformation : sig
    type treelengths_information = [
        | `Minimum
        | `Maximum
        | `Summary
    ]

    type file_information = [
        | `Filename 
    ]

    type character_information = [ `Type ]

    type search_script_information = [
        | `Done
        | `ToDo
        | `Filename
        | `Current
    ]

    type t = [
        | `CostMode
        | `TreeInformation of [ treelengths_information | `Number ] list
        | `TaxonInformation
        | `CharacterInformation of character_information list
        | `FileInformation of file_information list
        | `SearchInformation of search_script_information list
        | `Timer
    ]
end

type alph = 
    | Nucleotides
    | Aminoacids
    | GeneralAlphabet of 
        (string * Cost_matrix.Two_D.m * Cost_matrix.Three_D.m * Alphabet.a)


type name = string 

type kolmo_range = (float * float)

type kolmo_parameter_pairs = (string * float)

type kolmo_options = 
    [ `EProbability of  kolmo_parameter_pairs list
    | `FProbability of (string * kolmo_parameter_pairs list) ]

type kolmogorov_modules = 
    [ `SK of Kolmo.Compiler.sk_function list
    | `Alphabet of (name * (string list) * kolmo_options option)
    | `Character of (name * Kolmo.Compiler.sk_function list * kolmo_options option)
    | `WordSet of (name * name * kolmo_range * kolmo_options option)
    | `IntSet of (name * kolmo_range * kolmo_options option) ]

type d = {
    (* The number of terminals currently loaded *)
    number_of_taxa : int;
    (** The pairs of synonyms for the loaded taxa *)
    synonyms : string All_sets.StringMap.t;
    do_fs : bool;
    (** The current fixed states list to be used in whatever sequence is
    * loaded next in the dataset *)
    current_fs : Sequence.s list;
    (** The name of the fixed states file that contained the contents of
    * current_fs. *)
    current_fs_file : string;
    (** A function to generate codes for character *)
    character_code_gen : int ref;
    (* A function to generate codes for segments in chromosomes *)
    seg_code_gen : int ref;
    (* Set of files where each taxon appears *)
    taxon_files : All_sets.Strings.t All_sets.StringMap.t;
    (** A map between taxon names and their assigned codes *)
    taxon_names : int All_sets.StringMap.t;
    (** A map between the taxon codes and their corresponding names *)
    taxon_codes : string All_sets.IntegerMap.t;
    (* A mapping of the dynamic data to static character codes *)
    dynamic_static_codes : (int list) All_sets.IntegerMap.t;
    static_dynamic_codes : (int) All_sets.IntSetMap.t;
    (** A map of each taxon code and their corresponding character list *)
    taxon_characters : (int, (int, cs) Hashtbl.t) Hashtbl.t;
    (* A map of extra taxon for fixed_state median nodes*)
    searchbase_files : All_sets.Strings.t All_sets.StringMap.t;
    searchbase_names : int All_sets.StringMap.t;
    searchbase_codes : string All_sets.IntegerMap.t;
    searchbase_characters : (int, (int, cs) Hashtbl.t) Hashtbl.t;
    (**)
    (* A map between the character names and their corresponding codes *)
    character_names : (string, int) Hashtbl.t;
    (* A map between the character codes and their corresponding names *)
    character_codes : (int, string) Hashtbl.t;
    (* character set names to the list of charcter names *)
    character_sets : (string, string list) Hashtbl.t;
    (* A map between the character name and it's set *)
    character_nsets : (string, string) Hashtbl.t;
    (* A map between the character codes and their specifications *)
    character_specs : (int, specs) Hashtbl.t;
    (** The set of taxa to be ignored in the analysis *)
    ignore_taxa_set : All_sets.Strings.t;
    (* The set of taxa to be ignored in the analysis *)
    ignore_character_set : string list;
    trees : parsed_trees list;
    branches : (string,((string,float) Hashtbl.t) All_sets.IntSetMap.t) Hashtbl.t option;
    iterate_branches : bool;
    non_additive_1 : int list;
    non_additive_8 : int list;
    non_additive_16 : int list;
    non_additive_32 : int list;
    non_additive_33 : int list;
    additive : int list;
    sankoff : int list list;
    fixed_states : int list;
    dynamics : int list;
    kolmogorov : int list;
    static_ml : int list;
    complex_schema : Parser.SetGroups.t list;
    (** Tree for how to arrange taxa into complex terminals *)
    files : (string * contents list) list;
    machine : Kolmo.Compiler.compiler;
    search_information : OutputInformation.t list;

    (** At what taxon to root output trees *)
    root_at : int option;
}

(** [empty ()] creates a fresh empty dataset. *)
val empty : unit -> d

val duplicate : d -> d

val compare : d -> d -> bool

(** [add_character_spec s c x] returns fresh [d] generated by adding the
 * specification [s] to [x] with code [c]. *)
val add_character_spec : specs -> int -> d -> d

(** [to_channel ch b] outputs the specification of characters and loaded
 * characters of [b] in a human readable manner in channel [ch]. *)
val to_channel : out_channel -> d -> unit

(** [code_taxon c d] retrieves the name associated with a taxon code [c] in the
 * data [d]. *)
val code_taxon : int -> d -> string

(** [taxon_code n d] finds the code assigned to taxon [n] in the dataset 
 * [d]. *)
val taxon_code : string -> d -> int

val character_code : string -> d -> int

val code_character : int -> d -> string

val get_tcm : int -> d -> tcm

val get_weight : int -> d -> float

val get_weights : d -> (int * float) list

val process_parsed_sequences :
     bool -> float -> tcm_definition -> Cost_matrix.Two_D.m -> Cost_matrix.Three_D.m ->
        dyna_initial_assgn -> bool -> Alphabet.a -> string -> dyna_state_t ->
            d -> (Sequence.s list list list * Parser.E.taxon) list -> MlModel.model option
                -> dyna_pam_t option -> d

val process_molecular_file : ?respect_case:bool -> tcm_definition -> Cost_matrix.Two_D.m ->
    Cost_matrix.Three_D.m -> bool -> Alphabet.a -> dyna_initial_assgn-> bool -> 
        dyna_state_t -> d -> FileStream.f -> d

val add_static_file : ?report:bool -> [`Hennig | `Nexus] -> d -> FileStream.f -> d

val process_trees : d -> FileStream.f -> d

val process_fixed_states : d -> FileStream.f option -> d

val add_synonyms_file : d -> FileStream.f -> d

val add_synonym : d -> (string * string) -> d

val process_ignore_file : d -> FileStream.f -> d

val process_ignore_taxon : d -> string -> d

val process_analyze_only_file : bool -> d -> FileStream.f list -> d

val number_of_taxa : d -> int

val process_analyze_only_taxa : 
    [`Random of float | `CharSet of (bool * string list) 
    | `Names of (bool * string list) | `Missing of (bool * int) ] -> d -> d


(* Functions to manipulate and determine character sets and demarcation of data
 * sets by type, model, ... *)
val get_set_of_character : d -> int -> string option

val categorize : d -> d

val categorize_likelihood_chars_by_model : bool_characters -> d -> int list list

val categorize_sets : d -> int list list

val categorize_characters : d -> characters -> int list list

val categorize_characters_comp : d -> bool_characters -> int list list

val make_set_partitions : bool -> d -> string -> Methods.characters -> d

val make_codon_partitions : bool -> d -> string -> Methods.characters -> d

val remove_taxa_to_ignore : d -> d

val get_sequence_tcm : int -> d -> Cost_matrix.Two_D.m

val get_tcm2d : d -> int -> Cost_matrix.Two_D.m 
val get_tcm3d : d -> int -> Cost_matrix.Three_D.m
val get_tcmfile : d -> int -> tcm_definition

val get_sequence_alphabet : int -> d -> Alphabet.a

val add_file : d -> contents list -> FileStream.f -> d

val get_taxa : d -> string list

val add_static_parsed_file : d -> string -> Nexus.File.nexus -> d

val add_multiple_static_parsed_file :
    d -> (int option * (string * Nexus.File.nexus)) list -> d * int list

val characters_to_formatter :  d -> Xml.xml

val character_spec_to_formatter : specs -> Xml.xml

val ignored_taxa_to_formatter : d -> Xml.xml

val files_to_formatter : d -> Xml.xml

val taxon_names_to_formatter : d -> Xml.xml

val synonyms_to_formatter : d -> Xml.xml

val to_formatter : Xml.attributes -> d -> Xml.xml 

type classes = 
    [ `Fixedstates |`Dynamic |  `NonAdditive | `StaticLikelihood | `DynamicLikelihood | `Likelihood
    | `Additive | `Sankoff | `Kolmogorov | `AllStatic | `AllDynamic ] 

val get_code_from_characters_restricted : classes -> d -> characters -> int list

val get_code_from_characters_restricted_comp : classes -> d -> bool_characters -> int list 

val filter_non_static_approx_characters : ?comp:bool -> d -> int list -> int list

val can_all_chars_do_static_approx : d -> int list -> bool

val transform_dynamic : Methods.dynamic_char_transform -> d -> d

val transform_chrom_to_rearranged_seq :
  d -> Methods.dynamic_char_transform -> 'c -> Methods.implied_alignment list -> d

val print : d -> unit

val get_chars_codes : d -> characters -> int list

val get_chars_codes_comp : d -> bool_characters -> int list

val process_ignore_characters : bool -> d -> characters -> d

val process_ignore_characters_file : bool -> d -> FileStream.f -> d

val complement_characters : d -> characters -> [ `All | `Some of int list | `Names of string list ]

val complement_taxa : d -> int list -> int list

val process_analyze_only_characters_file : bool -> bool -> d -> FileStream.f list -> d

val process_rename_characters : d -> (string * string) -> d

val assign_transformation_gaps :
    d -> bool_characters -> int -> int -> d

val assign_level : d -> bool_characters -> Methods.keep_method -> int -> d 
    
val assign_affine_gap_cost : 
    d -> bool_characters -> Cost_matrix.cost_model -> d

val assign_tail : 
    d -> bool_characters -> [ `File of FileStream.f | `Array of int array ] -> d

val assign_prepend : 
    d -> bool_characters -> [ `File of FileStream.f | `Array of int array ] -> d

val assign_tcm_to_characters_from_file :
    d -> bool_characters -> (FileStream.f * (int*Methods.keep_method) option) option -> d

val add_search_base_from_file : 
    d -> bool_characters -> (string * string) list -> d

val process_complex_terminals : d -> FileStream.f -> d

val get_alphabet : d -> int -> Alphabet.a

val get_pam : d -> int -> dyna_pam_t

val get_character_state : d -> int -> dyna_state_t

val process_taxon_code : d -> All_sets.StringMap.key -> string -> d * int
 
val set_dyna_data : 'a seq_t array -> 'a dyna_data

val set_fs_data : 'a seq_t array -> int list -> 'a fixedstates_data

val get_recost : dyna_pam_t -> int

val get_locus_indel_cost : dyna_pam_t -> int * int

val get_character_set_name : d -> int list -> string option

val modified_characters : d -> d -> int

val get_likelihood_model : d -> int list -> MlModel.model

val apply_likelihood_model_on_chars : d -> int list -> MlModel.model -> d

val verify_alphabet : d -> int list -> Methods.ml_alphabet -> int * Alphabet.a

val update_priors : d -> int list -> bool -> d

val apply_heuristic_cost_model : ?cost_model:Methods.ml_costfn -> d -> d option

val set_parsimony  : d -> Methods.characters -> d 

val set_likelihood : d -> Methods.ml_spec    -> d 

val get_likelihood_model : d -> int list -> MlModel.model 

val report_taxon_file_cross_reference : 
    bool_characters option -> d -> string option -> unit

val report_terminals_files : 
    string option -> All_sets.Strings.t All_sets.StringMap.t -> All_sets.Strings.t ->
        unit

val dyna_pam_default : dyna_pam_t

val get_empty_seq: Alphabet.a -> Sequence.s seq_t

val find_max_seq_id : d -> int

val flush : d -> unit

val kolmo_round_factor : float

val transform_weight : 
    [ `ReWeight of (bool_characters * float) | `WeightFactor of (bool_characters * float) ] -> d -> d

val file_exists : d -> FileStream.f -> bool

val make_fixed_states : string option -> bool_characters -> polymorphism_t option -> d -> d

val make_direct_optimization : bool_characters -> d -> d

val make_partitioned : [`Clip | `NoClip] -> bool_characters -> d -> d

val has_dynamic : d -> bool

val has_likelihood: d -> bool

val type_of_dynamic_likelihood: d -> Methods.ml_costfn option

(* functions on translating data from static to dynamic *)
val convert_dynamic_to_static_branches : src:d -> dest:d -> d

val convert_static_to_dynamic_branches : src:d -> dest:d -> d

val sync_dynamic_to_static_model : src:d -> dest:d -> d

val sync_static_to_dynamic_model : src:d -> dest:d -> d

val remove_absent_present_encodings : ?ignore_data:bool -> d -> bool_characters -> d * int list

val randomize_taxon_codes : Methods.terminal_transform -> d -> d * (int, int) Hashtbl.t

module Sample : sig
    val generate : d -> [ `Bootstrap | `Jackknife of float ] -> d
end

type tcm_class =
    [ `AllOne of int
    | `AllOneGapSame of (int * int)
    | `AffinePartition of (int * int * int)
    | `AllSankoff of (string -> int) option]

val prealigned_characters :   
    (Cost_matrix.Two_D.m -> MlModel.model option -> Alphabet.a ->
              tcm_class * ([> `Exists ] -> int -> FileContents.t list -> 
                  FileContents.t list) *
                 (int -> (Alphabet.a * Parser.OldHennig.Encoding.s) list ->
                     (Alphabet.a * Parser.OldHennig.Encoding.s) list)) ->
                           d -> bool_characters -> d

(** [compare_all_pairs a b c d] compare for each taxon the characters with code
* [a] and [b], making the reverse complement of [b] if [c] is true, as stored in
* the data structure [d], returning for each character an optional tuple holding
* the number of taxa with comparable characters, and their overal comparison
* index. If none of the taxa hold the required characters, [None] is returned. *)
val compare_all_pairs : int -> int -> bool -> d -> (int * float) option

(** [compare_pairs a b c d] calculates the compare_all_pairs for every pair of
* characters in [a] and [b], calculating the reverse complement of [b] if [c] is
* true, as stored in the structure [d]. The output is a list holding for each
* pair of character names, the comparison index (see [compare_all_pairs]). *)
val compare_pairs : bool_characters -> bool_characters -> bool -> d -> 
    (string * string * float) list

type sequence_statistics = {
    max_length : int;
    min_length : int;
    sum_lengths : int;
    sequences : int;
    max_distance : float;
    min_distance : float;
    sum_distances : float; }

(** [sequence_statistics ch d] returns a list containing pairs consisting 
* of [(n, (w, x, y, z, d, e, f))], where [n] is the name of each character included in
* [ch], [w] is the maximum sequence length, [x] is the minimum sequence length,
* [y] is the sum of the length of all taxa for the character [n], [z] is the
* total number of characters included (taxa containing it), [d] is the maximum
* ammong the all pairs distances, [e] is the minimum ammong the all pairs
* distances, and [f] is the sum of all the distances. In total (z ^ 2 / 2 - z)
* distances are computed. *)
val sequence_statistics : bool_characters -> d -> (string * sequence_statistics) list

val to_human_readable : d -> int -> int -> string

(** [apply_boolean nadd add d c] returns the result of applying [nadd] or [add]
 * to the list of states observed for the character [c]  in the data [d]. [nadd]
 * is applied iff [c] is nonadditive, and [add] iff [c] is additive. For any
 * other character the function returns true.*)
val apply_boolean : 
    (Nexus.File.static_state list -> bool) -> 
        (Nexus.File.static_state list -> bool) -> d -> int -> bool

(** [get_model c d] returns the model of the ML character with code [c] in data
* [d]. If the character is not of type ML, it will raise an exception. *)
val get_model : int -> d -> MlModel.model
val get_model_opt : d -> int -> MlModel.model option

(** [min_max_possible_cost a b c d e] applies the functions [a], [b] and [c] in
 * the ordered, unordered, and sankoff characters respectively listed in [d], 
 * of all the terminals stored in [e], and returns the result per character in a
 * list of tuples holding the character code and the result. *)
val apply_on_static :
    (Nexus.File.static_state list -> float) -> (Nexus.File.static_state list -> 
        float) -> 
        (int array array -> Nexus.File.static_state list -> float) -> 
            (MlModel.model -> Nexus.File.static_state list -> float) -> bool_characters ->
            d -> (int * float) list

val repack_codes : d -> d

val verify_trees : d -> parsed_trees -> unit

val guess_class_and_add_file : bool -> bool -> d -> FileStream.f -> d

val report_kolmogorov_machine : string option -> d -> d
(*
val is_fs : d -> int -> bool *)
