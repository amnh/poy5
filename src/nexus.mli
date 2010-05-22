module P : sig
    type datatype = 
         DStandard | Dna | Rna | Nucleotide | Protein | Continuous 

    type item = 
         Min | Max | Median | Average | Variance | Stderror | SampleSize |
        States 

    type statesformat =
         StatesPresent | Individuals | Count | Frequency 

    type triangleformat =  Lower | Upper | Both 

    type format_options =  
        | Datatype of datatype  
        | RespectCase
        | FMissing of string
        | Gap of string
        | Symbols of string  (** still need to parse the symbols *)
        | Equate of string   (** still need to parse the tuples *)
        | MatchChar of string
        | Labels of bool
        | Transpose
        | Interleave
        | Items of item
        | StatesFormat of statesformat
        | Triangle of triangleformat
        | Tokens of bool 

    type charset = 
        | Range of (string * string option * int)
        | Single of string
        | Name of string
        | CharSet of string

    type char_data = {
        char_taxon_dimensions : string option;
        char_char_dimensions : string;
        char_format : format_options list;
        char_eliminate : charset option;
        char_taxlabels : string list;
        char_statelabels : (string * string * string list) list;
        char_charlabels : string list;
        char_charstates : (string * string list) list;
        chars : string;
    }

    type unalg_data = {
        unal_taxon_dimensions : string option;
        unal_format : format_options list;
        unal : string;
    }

    type standard_item = 
         Code of (string * charset list) | IName of (string * charset list) 

    type standard_list = 
        | STDVector of string list | STDStandard of charset list

    type set_type = 
        | Standard of standard_item list | Vector of string list 

    type set_pair =  
        | TaxonSet of charset list
        | CharacterSet of charset list
        | StateSet of charset list
        | TreeSet of charset list
        | CharPartition of set_type
        | TaxPartition of set_type
        | TreePartition of set_type

    type source =  Inline | File | Resource 

    type pictureformat =  Pict | Tiff | Eps | Jpeg | Gif 

    type pictureencoding =  None | UUEncode | BinHex 
    type polytcount =  MinSteps | MaxSteps 
    type gapmode =  Missing | NewState 

    type user_type =  StepMatrix of (string * string list) | CSTree of string 

    type assumption_set = (bool * string * bool * set_type)

    type assumption_items = 
        | Options of (string option * polytcount * gapmode)
        | UserType of (string * user_type)
        | TypeDef of assumption_set
        | WeightDef of assumption_set
        | ExcludeSet of (bool * string * standard_list)
        | AncestralDef of assumption_set

    type likelihood_model = 
        | Model of string
        | Variation of string
        | Variation_Sites of string
        | Variation_Alpha of string
        | Variation_Invar of string
        | Given_Priors of (string * float) list
        | Other_Priors of string 
        | Chars of charset list
        | Parameters of float list
        | GapMode of bool
        | Files of string

    type poy_data =          (* trees , characters, (nodes , length) *)
        | CharacterBranch of string list * charset list * (string * float) list
        | Likelihood of likelihood_model list
        | Tcm of (bool * string * standard_item list)
        | GapOpening of (bool * string * standard_item list)
        | DynamicWeight of (bool * string * standard_item list)

    type block = 
        | Taxa of (string * string list) 
        | Characters of char_data 
        | Distances of ((bool * string * string) option * format_options list * string list * string)
        | Ignore of string
        | Unaligned of unalg_data
        | Trees of (string * string) list * string list 
        | Notes of ((set_pair list * source * string) option * (set_pair list *
        pictureformat option * pictureencoding option * source * string) option) 
        | Assumptions of assumption_items list 
        | Error of string
        | Sets of (string * set_pair) list
        | Poy of poy_data list

    type tree_i = 
        | Leaf of (string * (float option * string option))
        | Node of (tree_i list * string option * (float option * string option))


    type tree = string * tree_i

    val print_error : (string -> unit) ref
end

module Grammar : sig
    type token
    val tree : (Lexing.lexbuf -> token) -> Lexing.lexbuf -> P.tree
    val header : (Lexing.lexbuf -> token) -> Lexing.lexbuf -> unit
    val block : (Lexing.lexbuf -> token) -> Lexing.lexbuf -> P.block
    val symbol_pair : 
        (Lexing.lexbuf -> token) -> Lexing.lexbuf -> (string * string list) 
    val symbol_list : 
        (Lexing.lexbuf -> token) -> Lexing.lexbuf -> string list
end

module Lexer : sig
    exception Eof
    val token : Lexing.lexbuf -> Grammar.token
    val tree_tokens : Lexing.lexbuf -> Grammar.token
end

module File : sig
    type st_type = 
        | STOrdered
        | STUnordered  
        | STSankoff of int array array   (* If Sankoff, the cost matrix to use *)
        | STLikelihood of MlModel.model  (* The ML model to use *)

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

    type static_state = [ `Bits of BitSet.t | `List of int list ]  option

    type taxon = string

    type nexus = {
        char_cntr : int ref;
        taxa : string option array;
        characters : static_spec array;
        matrix : static_state array array;
        csets : (string, P.charset list) Hashtbl.t;
        unaligned : 
            (float * int option * (string * int array array) option * 
             Alphabet.a * MlModel.model option * 
             (Sequence.s list list list * taxon) list)
            list;
        trees : (string option * Tree.Parse.tree_types list) list;
        branches : (string, (string, (string , float) Hashtbl.t) Hashtbl.t) Hashtbl.t;
        assumptions : (string, string array * float array array) Hashtbl.t;
    }

    val get_character_names : 
            static_spec array -> 
                (string, P.charset list) Hashtbl.t -> P.charset -> string list

    val empty_parsed : unit -> nexus

    val static_state_to_list : 
        [ `Bits of BitSet.t | `List of int list ] -> int list
    (** [spec_of_alph alphabet missing gap] generates a specification that can read
    * the elements in the [alphabet] using when the matrix represents [missing]
    * data and [gaps] as specified. *)
    val spec_of_alph : Alphabet.a -> string -> string -> static_spec

    (** [to_string v] outputs a string representation of the static_specification
    * [v] *)
    val to_string : static_spec -> string

    (** [to_formatter v] generates a standard [Xml.xml] representation of
    * [to_formatter]. *)
    val to_formatter : static_spec -> Xml.xml

    val make_symbol_alphabet : string -> string list ->  (string * string list) list
    -> P.format_options list -> Alphabet.a * (string * string list) list 

    val process_matrix : 
       bool ->
       [ `Hennig | `Nexus | `None ] ->
       static_state array array ->
       string option array ->
       static_spec array ->
       (string -> int) ->
       (int -> int -> static_state -> unit) -> string -> unit

    val find_taxon : string option array -> string -> int

    val of_channel : in_channel -> string -> nexus

    val generate_alphabet : string list -> string -> Alphabet.a

    val spec_of_alph : Alphabet.a -> string -> string -> static_spec
    val fill_observed : 
       static_spec array ->
       [< `Bits of BitSet.t | `List of All_sets.Integers.elt list ]
       option array array -> unit

end
