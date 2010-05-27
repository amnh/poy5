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

let () = SadmanOutput.register "Xml" "$Revision: 2554 $"

type tag = string

type unstructured = [ `Bool of bool 
    | `IntFloatTuple of (int * float) 
    | `IntTuple of (int * int)
    | `FloatFloatTuple of (float * float) 
    | `String of string 
    | `Int of int 
    | `Float of float 
    | `Fun of (unit -> string) ]

type 'b structured =
    [ `Delayed of (unit -> 'b Sexpr.t)
    | `Empty
    | `Set of 'b Sexpr.t list
    | `Single of 'b ]

type 'a structured_xml =
    [ 'a structured 
    | `CDATA of [ unstructured | 'a structured] ]

type attribute = tag * unstructured
type attributes = attribute list
type 'a contents = [ unstructured | 'a structured_xml]
type xml = (tag * attributes * xml contents)

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

    let priors = "Priors"
    let llike = "LnLikelihood"
    let model = "Model"
    let sites = "Sites"
    let characters = "Characters"
    let alpha = "Alpha"
    let invar = "Invarient"
    let vector = "Vector"
    let categories = "Rate Classes"
    let gapascharacter = "Gap Character"

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
    let initial_assignment = "Initial Assignment"
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
    let max_kept_wag = "max_kept_wag"
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
    let min_time = "Minimum Code Child Time"
    let oth_time = "Other Child Time"
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
    let param x = "Parameter_"^(string_of_int x)
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

let attribute name ((_, atts, _) : xml) =
    List.assoc name atts

let tag ((t, _, _) : xml) = t

let contents ((_, _, c) : xml) = c

let attributes ((_, a, _) : xml) = a

let rec children search_tag ((a, b, chld) : xml) = 
    match chld with
    | `Delayed f -> 
            children search_tag (a, b, (f () :> xml contents))
    | #unstructured 
    | `CDATA _ -> `Empty
    | #Sexpr.t as chld ->
            Sexpr.filter (fun x -> search_tag = tag x) chld

let value ((_, _, c) : xml) =
    match c with
    | #unstructured as v -> v
    | _ -> failwith "Not a value"

let eagerly_compute x = 
    match x with
    | #Sexpr.t as x -> x
    | `Delayed f -> f ()

let structured (_, _, c) =
    match c with
    | #structured as x -> x
    | _ -> failwith "Not structured"


let coherce x = (x :> xml contents)

let make tag atv out = (tag, atv, out)

let remove_non_alpha_numeric str =
    Str.global_replace (Str.regexp " ") "_"
    (Str.global_replace (Str.regexp "[^0-9a-zA-Z]") "_" str)

let get_float unstr = 
    match unstr with
    | `Float x -> x
    | `Int x -> x
    | _ -> failwith "get float/int only, in xml.ml"

let value_to_string = function
    | `Bool x -> string_of_bool x
    | `IntTuple (a, b) -> string_of_int a ^ "," ^ string_of_int b
    | `IntFloatTuple (a, b) -> string_of_int a ^ "," ^ string_of_float b
    | `FloatFloatTuple (a, b) -> string_of_float a ^ "," ^ string_of_float b
    | `String x -> x
    | `Int x -> string_of_int x
    | `Float x -> string_of_float x
    | `Fun x -> x ()

let print_string ch = function
    | `Bool x -> Printf.fprintf  ch "%B" x
    | `IntTuple (a, b) -> Printf.fprintf ch "%d,%d" a b
    | `IntFloatTuple (a, b) -> Printf.fprintf ch "%d,%f" a b
    | `FloatFloatTuple (a, b) -> Printf.fprintf ch "%f,%f" a b
    | `String x -> Printf.fprintf ch "%s" x
    | `Int x -> Printf.fprintf ch "%d" x
    | `Float x -> Printf.fprintf ch "%s" (string_of_float x)
    | `Fun x -> Printf.fprintf ch "%s" (x ())


let to_file ch (item : xml) =
    let fo = output_string ch in
    let output_attrs (fo : string -> unit) (a, b) =
        fo " ";
        fo (remove_non_alpha_numeric a);
        fo "=\"";
        print_string ch b;
        fo "\""
    in
    let stack = Stack.create () in
    let simplify (contents : xml contents) = 
        match contents with
        | #structured as x -> 
                let x = eagerly_compute x in
                `Structured (Sexpr.to_list x)
        | #unstructured as x -> x 
        | `CDATA (#unstructured as x) -> `CDATA (Some x)
        | `CDATA (#structured as x) ->
                let x = eagerly_compute x in
                `CDATA (Some (`Structured (Sexpr.to_list x)))
    in
    let process_queue () =
        while not (Stack.is_empty stack) do
            let close_tag tag =
                match tag with
                | None -> ()
                | Some tag -> fo "</"; fo tag; fo ">\n"
            in
            let (tag, contents) = Stack.pop stack in
            match contents with
            | #unstructured as x -> 
                    print_string ch x;
                    close_tag tag;
            | `CDATA None -> fo "]]>";
            | `CDATA (Some (#unstructured as v)) ->
                    fo "<![CDATA[";
                    Stack.push (None, `CDATA None) stack;
                    Stack.push (None, v) stack;
            | `CDATA (Some ((`Structured _) as t)) ->
                    Stack.push (tag, `Structured []) stack;
                    fo "<![CDATA[";
                    Stack.push (None, `CDATA None) stack;
                    Stack.push (None, t) stack;
            | `Structured []
            | `Empty -> close_tag tag
            | `Structured ((ntag, attributes, contents) :: t) ->
                    Stack.push (tag, `Structured t) stack;
                    let ntag = remove_non_alpha_numeric ntag in
                    fo "<"; fo ntag; 
                    List.iter (output_attrs fo) attributes;
                    fo ">";
                    let contents = simplify contents in
                    Stack.push ((Some ntag), contents) stack
        done;
    in
    fo "";
    let tag, attributes, contents = item in
    let tag = remove_non_alpha_numeric tag in
    fo " <"; fo tag; fo " ";
    List.iter (output_attrs fo) attributes;
    fo ">\n";
    let contents = simplify contents in
    Stack.push ((Some tag), contents) stack;
    process_queue ()

