(* POY 5.0. A phylogenetic analysis program using Dynamic Homologies.         *)
(* Copyright (C) 2013 Andrés Varón, Lin Hong, Nicholas Lucaroni, Ward Wheeler,*)
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
    | Range of (int * int option * int)
    | Single of int
    | Name of string

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
    | Cost_Mode of string
    | Chars of charset list
    | Parameters of float list
    | Gap_Mode of (string * float option)
    | Files of string

type character_data =
    | Tree_Names of string list
    | Set_Names of charset list
    | Labeling of (string * float) list 

type annot_data =
    | Annot_Quality of float
    | Annot_Min of int
    | Annot_Max of int
    | Annot_Min_Percent of float
    | Annot_Max_Percent of float
    | Annot_Coverage of float
    | Annot_Type of [`Mauve | `Default]
    | Annot_Rearrangement of int

type chrom_data =
    | Chrom_Solver of string
    | Chrom_Locus_Indel of int * float
    | Chrom_Locus_Breakpoint of int
    | Chrom_Locus_Inversion of int
    | Chrom_Approx of bool
    | Chrom_Median of int
    | Chrom_Symmetric of bool
    | Chrom_Annotations of annot_data list

type genome_data = 
    | Genome_Median of int
    | Genome_Indel of int * float
    | Genome_Circular of bool
    | Genome_Breakpoint of int
    | Genome_Distance of float

type poy_data =  
    | Chrom of chrom_data list * charset list
    | Genome of genome_data list * charset list
    | BreakInv of chrom_data list * charset list
    | CharacterBranch of character_data list
    | Likelihood of likelihood_model list
    | Tcm of (bool * string * standard_item list)
    | GapOpening of (bool * string * standard_item list)
    | DynamicWeight of (bool * string * standard_item list)
    | Level of (bool * string * standard_item list)

type block = 
    | Taxa of (string * string list) 
    | Characters of char_data 
    | Distances of ((bool * string * string) option * format_options list * string list * string)
    | Ignore of string
    | Unaligned of unalg_data
    | Trees of (string * string) list * string list 
    | Notes of 
        ((set_pair list * source * string) option * (set_pair list *
            pictureformat option * pictureencoding option * source * string) option) 
    | Assumptions of assumption_items list 
    | Error of string
    | Sets of (string * set_pair) list
    | Poy of poy_data list

type tree_i = 
    | Leaf of (string * (float option * string option))
    | Node of (tree_i list * string option * (float option * string option))

type tree = string * tree_i

let print_error = ref prerr_string

