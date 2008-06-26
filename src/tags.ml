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

let () = SadmanOutput.register "Tags" "$Revision: 2554 $"

type tag = string
type value = tag
type attribute = tag * value
type attributes = attribute list
type output = 
    (tag * attributes * [`String of string | `Structured of output Sexpr.t])

let make tag atv out = (tag, atv, out)

let remove_non_alpha_numeric str =
    Str.global_replace (Str.regexp "[^0-9a-zA-Z]") "_" str

let to_xml fo item =
    let output_attrs (fo : string -> unit) (a, b) =
        fo (remove_non_alpha_numeric a);
        fo "=\"";
        fo b;
        fo "\" "
    in

    let rec to_xml fo ((tag, attributes, contents) : output) : unit = 
        fo " <";
        fo (remove_non_alpha_numeric tag);
        fo " ";
        List.iter (output_attrs fo) attributes;
        fo ">@\n%!";
        begin match contents with
        | `String x -> 
              fo x; 
              fo "@\n%!"
        | `Structured x ->
                Sexpr.leaf_iter (to_xml fo) x;
        end;
        fo ( " </" ^ (remove_non_alpha_numeric tag) ^ ">@\n%!")
    in
    fo "";
    to_xml fo item;
    fo "%!"

module Alphabet = struct
    let element = "Element"
    let value = "Value"
    let code = "Code"
end

module Characters = struct
    let suffix = " Character"
    let character = "Character"
    (* The characters themselves *)
    let additive = "Additive" ^ suffix
    let nonadditive = "Non Additive" ^ suffix
    let molecular = "Molecular" ^ suffix
    let sankoff = "Sankoff" ^ suffix
    let likelihood = "Likelihood" ^ suffix
    let set = "Set" ^ suffix
    let kolmogorov = "Kolmogorov" ^ suffix

    (* Their attributes *)
    let name = "Name"    
    let cost = "Cost"
    let recost = "Rearrangment Cost"
    let definite = "Definite"
    let weight = "Weight"
    let cclass = "Class"
    let alphabet = "Alphabet"
    let words = "Words"
    let ints = "Integers"
    let chars = "Character Functions"

    let prior = "Prior"
    let p_mat = "Probability_Matrix"
    let model = "Model"
    let mle = "Likelihood"
    let p_vec = "Likelihood_Vector"
    let characters = "Characters"

    (* Their values *)
    let min = "Min"
    let max = "Max"
    let value = "Value"
    let sequence = "Sequence"
    let annchrom = "Annotated chromosome"
    let chromosome = "Chromosome"
    let genome = "Genome"
    let breakinv = "Breakinv"

    (* The cost matrices and calculation methods *)
    let fixed_states = "Fixed States"
    let tcm = "Transformation Cost Matrix"
    let gap_opening = "Gap Opening"
    let state = "State"
    let states = "States"
    let ref_code = "ReferenceCode"
    let chrom_map = "Pairwise_Alignment_Map"

    (* The labels of the states of the characters *)
    let label = "Label"
    let labels = "Labels"
    let item = "Item"
    let observed = "Observed"
    let source = "File Source"
    let missing_symbol = "Missing Symbol"
    let matchstate_symbol = "Match Symbol"
    let gap_symbol = "Gap Symbol"
    let ignore = "Ignore"
    let case = "Case Sentive"  

    let equivalencies = "Equivalencies"
    let equivalent = "Equivalency"
    let from = "From"
    let towards = "To"

    (* All the chromosomal information character contents *)
    let clas = "class"
    let sequence = "sequence"
    let chromosome = "chromosome"
    let genome = "genome"
    let annotated = "annotated"
    let breakinv = "breakinv"
    let seed_len = "seed_len"
    let re_meth = "re_meth"
    let circular = "circular"
    let locus_indel_cost = "locus_indel_cost"
    let chrom_indel_cost = "chrom_indel_cost"
    let chrom_hom = "chrom_hom"
    let chrom_breakpoint = "chrom_breakpoint"
    let sig_block_len = "min_loci_len"
    let rearranged_len = "rearranged_len"
    let keep_median = "keep_median"
    let swap_med = "swap_med"
    let approx = "approx"
    let symmetric = "symmetric"
    let max_3d_len = "max_3d_len"
end

module Nodes = struct
    let node = "Node"
    let preliminary = "Preliminary"
    let single = "Single"
    let final = "Final"
    let cost = Characters.cost
    let recost = Characters.recost
    let node_cost = "Node cost"
    let name = Characters.name
    let nce = "Number of Children Edges"
    let notu = "Number of OTUS"
    let child1_name = "Child1's name"
    let child2_name = "Child2's name"
    let time = "Time"
end

module Trees = struct
    let forest = "Forest"
    let tree = "Tree"
    let cost = Nodes.cost    
    let recost = Nodes.recost
end

module Data = struct
    let cost_matrix = "Cost Matrix"
    let gap_opening = Characters.gap_opening
    let trees = "Trees"
    let data = "Data"
    let synonyms = "Synonyms"
    let synonym = "Synonym"
    let value =  Characters.value
    let code = "Code"
    let name = Nodes.name
    let taxon = "Taxon"
    let taxa = "Taxa"
    let file_contents = "File Contents"
    let filename = "File Name"
    let file = "File"
    let files = "Files"
    let ignored_taxa = "Ignored Taxa"
    let ignored_characters = "Ignored Characters"
    let characters = "Characters"
    let modeltype = "Model_Type"
end

(* The Kolmogorov complexity characters specifications *)
module KolSpecs = struct
    let spec_index = "Specifications Index"
    let char_index = "Character Index"
    let set_spec = "Set Specification"
    let char_spec = "Character Specification"
    let spec_class = "Class"
    let alphabet = "Alphabet"
    let integers = "Integers"
    let words = "Words"
    let spec_name = "Name"
    let char_name = spec_name
    let fun_name = spec_name
    let char_fun = "Function"

    let alph_element = "Alphabet Element"
    let value = "Value"
    let prob = "Nats"
    let int_set = "Integer Set"
    let word_set = "Word Set"
    let min = "Min"
    let max = "Max"

    let model = "model"
    let indelsonly =  "Insertions and Deletions"
    let indelsub =  "Insertions, Deletions, and Substitutions"
    let subsonly =  "Substitutions"
    let affineindelssub = "Affine insertions and deletions, and atomic substitutions"
    let affineindelsaffsub = "Affine insertions, deletions, and substitutions"
end


module GenomeMap = struct
    let genome = "GenomeMap"
    let chrom = "ChromosomeMap"
    let seg = "SegmentMap"

    let ref_code = "ReferenceCode"
    let a_ref_code = "AncestorReferenceCode"
    let d_ref_code = "DescendantReferenceCode"


    let a_chrom_id = "AncestorChromosomeID"
    let d_chrom_id = "DescendantChromosomeID"

    let seq_order = "SequenceOrder"
    let a_seq_order = "AncestorSequenceOrder"
    let d_seq_order = "DescendantSequenceOrder"

    let start_seg = "StartPosition"
    let a_start_seg = "AncestorStartPosition"
    let d_start_seg = "DescendantStartPosition"

    let end_seg = "EndPosition"
    let a_end_seg = "AncestorEndPosition"
    let d_end_seg = "DescendantEndPosition"

    let dir_seg = "Direction"        
    let a_dir_seg = "AncestorDirection"        
    let d_dir_seg = "DescendantDirection"        
end
