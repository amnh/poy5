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
    Range of (string * string option) | Single of string | Name of string

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

type set_pair =  
    | TaxonSet of charset 
    | CharacterSet of charset 
    | StateSet of charset 
    | TreeSet of charset 

type source =  Inline | File | Resource 

type pictureformat =  Pict | Tiff | Eps | Jpeg | Gif 

type pictureencoding =  None | UUEncode | BinHex 
type polytcount =  MinSteps | MaxSteps 
type gapmode =  Missing | NewState 

type user_type =  StepMatrix of (string * string list) | CSTree of string 

type standard_item = 
     Code of (string * charset list) | IName of (string * charset list) 
type set_type =  Standard of standard_item list | Vector of string list 

type assumption_set = (bool * string * bool * set_type)

type assumption_items = 
     Options of (string option * polytcount * gapmode)
    | UserType of (string * user_type)
    | TypeDef of assumption_set
    | WeightDef of assumption_set
    | ExcludeSet of assumption_set
    | AncestralDef of assumption_set

type block = 
     Taxa of (string * string list) 
    | Characters of char_data 
    | Distances of ((bool * string * string) option * format_options list * string
    list * string)
    | Ignore of string
    | Unaligned of unalg_data
    | Trees of ((string * string) list * string list) 
    | Notes of ((set_pair list * source * string) option * (set_pair list *
    pictureformat option * pictureencoding option * source * string) option) 
    | Assumptions of assumption_items list 
    | Error of string

type tree = 
    | Leaf of (string * float option)
    | Node of (tree list * string option * float option)
