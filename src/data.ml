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

include Dyn_pam (* TODO: this needs to be seperated properly. *)

open StdLabels

exception Nam_Collision of (CharacSpec.t * string)
exception Illegal_argument

type filename = string 

module FullTupleMap = All_sets.FullTupleMap
module IntMap = All_sets.IntegerMap
    
let ( --> ) a b = b a
let output_error = Status.user_message Status.Error
let output_warning = Status.user_message Status.Warning
let output_info = Status.user_message Status.Information
let failwithf format = Printf.ksprintf failwith format

let output_errorf format = 
    Printf.ksprintf (Status.user_message Status.Error) format
let output_warningf format = 
    Printf.ksprintf (Status.user_message Status.Warning) format
let output_infof format = 
    Printf.ksprintf (Status.user_message Status.Information) format

let debug_kolmo = false
let debug_level = false
let debug_search_base = false
let debug_parsed_seq = false 

type dynhom_opts = 
    | Tcm of string     (* A transformation cost matrix to be used *)

(** The valid types of contents of a file *)
type contents = Characters | CostMatrix | Trees 

(* (name * tree ) * file * num *)
type parsed_trees = ((string option * Tree.Parse.tree_types list) * string * int)

type dyna_state_t = [
    | `SeqPrealigned
    (** A short sequence, no rearrangements are allowed*)
    | `Seq
    (** same as above but for likelihood characters *)
    | `Ml
    (** A long sequence, genes inside are broken automatically, rearrangements are allowed*)
    | `Chromosome
    (** A set of chromosomes *)
    | `Genome
    (** A list of shorted sequences annotated by piles, rearrangements are allowed *)
    | `Annotated
    (** A sequence of gene names, rearrangements are allowed *)
    | `Breakinv ]

type polymorphism_t = Methods.polymorphism_arg

type clip = Clip | NoClip

type fixed_state =
    {   costs   : float array array;
        seqs    : Sequence.s array;
        codes   : (int, int) Hashtbl.t;
        opt_bls : float array array option;
    }

type dyna_initial_assgn =
    [ `Partitioned of clip
    | `AutoPartitioned of clip * int * (int,  ((int * int) list)) Hashtbl.t
    | `GeneralNonAdd
    | `DO
    | `FS of fixed_state ]

type tcm_definition =
    | Substitution_Indel of (int * int)
    | Input_file of (string * (int list list))
    | Substitution_Indel_GapOpening of (int * int * int)
    | Input_file_GapOpening of (string * (int list list) * int)
    | Level of (tcm_definition * int)

let default_tcm = Substitution_Indel (1,1)

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

type distr = | MaxLength of int (* Any of the distributions with a maximum length *)

type affine_f = {
    selfp : float;
    distr : distr;
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
    tm : float list list;       (** The transformation cost matrix *)
    be : arr;                   (** Base encodings, and insertion in tips *)
    simplebe : float array;     (** The original encoding cost *)
    simplebed : float array;    (** The deletion cost in the tips *)
    ins : float;                (** Insertion encoding cost *)
    del : float;                (** Deletion encoding cost *)
    sub : float;                (** Substitution encoding cost *)
    ins_opening : float;        (** Insertion opening *)
    del_opening : float;        (** Deletion opening *)
    sub_opening : float;        (** Substitution opening *)
    root_cost : float;          (** Extra cost incurred by each root *)
    branch_cost : float;        (** Extra cost incurred by each branch *)
    leaf_cost : float;          (** Extra cost incurred by marking a leaf *)
    end_cost : float;           (** Extra cost of ending the compuations and producing the output *)
    mo : model;                 (** The model *)
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

(** A character specification. Contains the information regarding the
    characteristics of a particular character *)
type specs = 
    (** Static homology characters, includes the encoding specs from the
        parser and the name of the file it comes from (the name includes the
        column number). *)
    | Static of Nexus.File.static_spec

    (** A dynamic homology based character type, with three parameters, the
        file name containing the set of sequences, the filename of the valid 
        fixed states that can be used for that set of sequences, and the file
        containing the transformation cost matrix to be used to perform their
        alignments.  Also, we store the alphabet used. *)
    | Dynamic of dynamic_hom_spec

    (** A set of characters that could be of several types. Probably some
        information is needed, but in the current implementation there is no
        need. As we get to more complex character types, this will change. *)
    | Set 

    (** The characters that should be analyzed using Kolmogorov complexity
        as optimality criterion. The only parameters really needed are the name
        of the character, the specification of the alphabet the character uses,
        the wordset the character uses, the integer set for functions defined in
        the character definition itself. Note that the words stored must also be
        specified in some file. *)
    | Kolmogorov of kolmo_spec

(** [specified] tells us whether or not the user left out the given character
    in the current taxon *)
type specified = [ | `Specified | `Unknown ]

(** typically 'a is a sequence, the integer followed is the code of 'a *)
type 'a seq_t = {
    seq : 'a;
    delimiter : int list;
    code : int;
}

(** wrap the seq_t structure to an array for partitioned data *)
type 'a dyna_data = {
    seq_arr : ('a seq_t) array;
}

(** A valid loaded character from a file *)
type cs_d = 
    (** A dynamic homology, containing its code, the sequence, the
        transformation cost matrix and its three dimensional transformation
        cost matrix *)
    | Dyna of (int * Sequence.s dyna_data)
    (** A static homology character, containing its code, and the character
        itself *)
    | Stat of (int * Nexus.File.static_state)

type cs = cs_d * specified

(* A transformation cost matrix for Sankoff characters *)
type tcm = int array array

module OutputInformation = struct
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
    (* The pairs of synonyms for the loaded taxa *)
    synonyms : string All_sets.StringMap.t;
    (* Whether or not fixed states is to be used in the following molecular
    * character *)
    do_fs : bool;
    (* The current fixed states list to be used in whatever sequence is
    * loaded next in the dataset *)
    current_fs : Sequence.s list;
    (* The name of the fixed states file that contained the contents of
    * current_fs. *)
    current_fs_file : string;
    (* A function to generate codes for character *)
    character_code_gen : int ref;
    (* A function to generate codes for segments in chromosomes *)
    seg_code_gen : int ref;
    (* Set of files where each taxon appears *)
    taxon_files : All_sets.Strings.t All_sets.StringMap.t;
    (* A map between taxon names and their assigned codes *)
    taxon_names : int All_sets.StringMap.t;
    (* A map between the taxon codes and their corresponding names *)
    taxon_codes : string All_sets.IntegerMap.t;
    (* A mapping of the dynamic data to static character codes *)
    dynamic_static_codes : (int list) All_sets.IntegerMap.t;
    static_dynamic_codes : (int) All_sets.IntSetMap.t;
    (* A map of each taxon code and their corresponding character list *)
    taxon_characters : (int, (int, cs) Hashtbl.t) Hashtbl.t;
    (* A map of extra taxon for median nodes*)
    searchbase_files : All_sets.Strings.t All_sets.StringMap.t;
    searchbase_names : int All_sets.StringMap.t;
    searchbase_codes : string All_sets.IntegerMap.t;
    searchbase_characters : (int, (int, cs) Hashtbl.t) Hashtbl.t;
    (**)
    (* A map between the character names and their corresponding codes *)
    character_names : (string, int) Hashtbl.t;
    (* A map between the character codes and their corresponding names *)
    character_codes : (int, string) Hashtbl.t;
    (* A map between the character set name and it's contents *)
    character_sets : (string, string list) Hashtbl.t;
    (* A map between the character name and it's set *)
    character_nsets : (string, string) Hashtbl.t;
    (* A map between the character codes and their specifications *)
    character_specs : (int, specs) Hashtbl.t;
    (* The set of taxa to be ignored in the analysis *)
    ignore_taxa_set : All_sets.Strings.t;
    (* The set of taxa to be ignored in the analysis *)
    ignore_character_set : string list;
    (* The set of loaded trees *)
    trees : parsed_trees list;
    (* branch lengths read from nexus file: tree -> node -> char *)
    (* node is defined by partition of leaves in the tree        *)
    branches : (string, ((string, float) Hashtbl.t) All_sets.IntSetMap.t) Hashtbl.t option;
    (* determine if we should iterate the branches *)
    iterate_branches : bool;
    (* The set of codes that belong to the class of Non additive with up to 1
    * states (useless!) *)
    non_additive_1 : int list;
    (* The set of codes that belong to the class of Non additive with up to 8
    * states *)
    non_additive_8 : int list;
    (* The set of codes that belong to the class of Non additive with up to 16
    * states *)
    non_additive_16 : int list;
    (* The set of codes that belong to the class of Non additive with up to 32
    * states *)
    non_additive_32 : int list;
    (* The set of codes that belong to the class of Non additive with up to 33
    * states *)
    non_additive_33 : int list;
    (* The set of codes that belong to the class of additive characters *)
    additive : int list;
    (* The set of codes that belong to the class of sankoff characters *)
    sankoff : int list list;
    (* The set of codes that belong to the class of sequence characters *)
    dynamics : int list;
    (* The set of codes that belong to the class of Kolmogorov characters *)
    kolmogorov : int list;
    (* The set of codes that belong to the class of Likelihood characters *)
    static_ml : int list;
    (** Tree for how to arrange taxa into complex terminals *)
    complex_schema : Parser.SetGroups.t list;
    (** The association list of files and kind of data they could hold *)
    files : (string * contents list) list;
    (* An index of characters defined using the MDL principle *)
    machine : Kolmo.Compiler.compiler;
    (* The search information to be presented to the user *)
    search_information : OutputInformation.t list;
    (* At what taxon to root output trees *)
    root_at : int option;
}

(** [compare d1 d2] Compare two data; we only return a boolean since ordering
    would be completely undefined and chaotic. The comparison ignores things
    that are not relevant to the DIAGNOSIS of the tree. This is to avoid
    rediagnosising a tree, even though it has inconsequential details changed in
    the parallel pipeline (scripting module). *)
let compare data1 data2 : bool =
    (data1.number_of_taxa = data2.number_of_taxa)
        && (data1.synonyms = data2.synonyms)
        && (data1.do_fs = data2.do_fs)
        && (data1.current_fs = data2.current_fs)
        && (data1.current_fs_file = data2.current_fs_file)
        && (data1.taxon_files = data2.taxon_files)
        && (data1.taxon_names = data2.taxon_names)
        && (data1.taxon_codes = data2.taxon_codes)
        && (data1.taxon_characters = data2.taxon_characters)
        && (data1.searchbase_files = data2.searchbase_files)
        && (data1.searchbase_names = data2.searchbase_names)
        && (data1.searchbase_codes = data2.searchbase_codes)
        && (data1.searchbase_characters = data2.searchbase_characters)
        && (data1.character_names = data2.character_names)
        && (data1.character_codes = data2.character_codes)
        && (data1.character_sets = data2.character_sets)
        && (data1.character_nsets = data2.character_nsets)
        && (data1.character_specs = data2.character_specs)
        && (data1.ignore_taxa_set = data2.ignore_taxa_set)
        && (data1.ignore_character_set = data2.ignore_character_set)
        && (data1.iterate_branches = data2.iterate_branches)
        && (data1.complex_schema = data2.complex_schema)
        && (data1.files = data2.files)
        && (data1.machine = data2.machine)

type bool_characters = Methods.characters (* [
    | `All
    | `Some of (bool * int list)
    | `Names of (bool * string list)
    | `CharSet of (bool * string list)
    | `Random of float
    | `AllStatic
    | `AllDynamic
    | `Missing of (bool * int)
    | `Range of (bool * string * int * int)
] *)

type characters = [
    | `All
    | `Some of int list 
    | `CharSet of string list
    | `Names of string list
    | `Random of float
    | `AllStatic
    | `AllDynamic
    | `Missing of (bool * int)
    | `Range of (string * int * int)
]

let transform_range_to_codes file x y =
    assert( x < y );
    let rec loop_ acc x y =
        if x > y then acc
        else loop_ ((file^":"^(string_of_int x))::acc) (x+1) y
    in
    `Names (loop_ [] x y)

let create_ht () = Hashtbl.create 1667 

(* Empty data! *)
let empty () = 
    {
        number_of_taxa = 0;
        synonyms = All_sets.StringMap.empty;
        do_fs = false;
        current_fs = [];
        current_fs_file = "";
        character_code_gen = ref 0;
        seg_code_gen = ref 0;
        taxon_names = All_sets.StringMap.empty;
        taxon_files = All_sets.StringMap.empty;
        taxon_codes = All_sets.IntegerMap.empty;
        dynamic_static_codes = All_sets.IntegerMap.empty;
        static_dynamic_codes = All_sets.IntSetMap.empty;
        taxon_characters = create_ht ();
        searchbase_names = All_sets.StringMap.empty;
        searchbase_files = All_sets.StringMap.empty;
        searchbase_codes = All_sets.IntegerMap.empty;
        searchbase_characters = create_ht ();
        character_names = create_ht ();
        character_codes = create_ht ();
        character_specs = create_ht ();
        character_sets = create_ht ();
        character_nsets = create_ht ();
        ignore_taxa_set = All_sets.Strings.empty;
        ignore_character_set = [];
        trees = [];
        branches = None;
        iterate_branches = true;
        non_additive_1 = [];
        non_additive_8 = [];
        non_additive_16 = [];
        non_additive_32 = [];
        non_additive_33 = [];
        additive = [];
        sankoff = [];
        dynamics = [];
        kolmogorov = [];
        static_ml = [];
        files = [];
        machine = Kolmo.Compiler.compiler;
        search_information = [`TreeInformation [`Summary]; `CostMode];
        root_at = None;
        complex_schema = [];
}


let copy_taxon_characters tc = 
    let new_tc = create_ht () in
    Hashtbl.iter (fun code othertbl ->
        Hashtbl.replace new_tc code (Hashtbl.copy othertbl)) tc;
    new_tc

let copy_branches = function
    | Some branches ->
        let tbl = Hashtbl.create (Hashtbl.length branches) in
        Hashtbl.iter
            (fun k v_map ->
                let n_map = All_sets.IntSetMap.map (fun v -> Hashtbl.copy v) v_map in
                Hashtbl.add tbl k n_map) 
            branches;
        Some tbl
    | None   -> None


let get_set_of_character data char : string option =
    let char_name = Hashtbl.find data.character_codes char in
    let all = Hashtbl.find_all data.character_nsets char_name in
    match all with
    | []  -> None
    | [x] -> assert( List.mem char_name (Hashtbl.find data.character_sets x) );
             Some x
    |  _  -> failwithf "Character %d/%s is associated with multiple sets" char char_name


let duplicate data = 
    { data with
        character_code_gen = ref !(data.character_code_gen);
        taxon_characters = copy_taxon_characters data.taxon_characters;
        searchbase_characters = copy_taxon_characters
        data.searchbase_characters;
        character_names = Hashtbl.copy data.character_names;
        character_codes = Hashtbl.copy data.character_codes;
        character_specs = Hashtbl.copy data.character_specs;
        character_sets  = Hashtbl.copy data.character_sets;
        character_nsets = Hashtbl.copy data.character_nsets;
        branches = copy_branches data.branches;
    }

(* recreate static_dynamic codes from dynamic_static map *)
let reverse_dynamic_static_codes map = 
    All_sets.IntegerMap.fold
        (fun i lst acc ->
            let key = List.fold_right ~f:All_sets.Integers.add
                                      ~init:All_sets.Integers.empty lst in
            All_sets.IntSetMap.add key i acc)
        map 
        All_sets.IntSetMap.empty

(* [convert_dynamic_to_static_branches src dest] Use the static_dynamic_codes
 * map from destination, and the branches structure from source. This behavior
 * is used after an implied alignment applies the dynamic_static_codes map to
 * the newly created data (dest), and needs to copy over and update the
 * branches. *)
let convert_dynamic_to_static_branches ~src ~dest = 
    let add_data desttbl char_name dat =
        let c_code =Hashtbl.find src.character_names char_name in
        if IntMap.mem c_code dest.dynamic_static_codes then
            List.iter
                (fun c ->
                    let n = Hashtbl.find dest.character_codes c in
                    Hashtbl.add desttbl n dat)
                (IntMap.find c_code dest.dynamic_static_codes)
        else ()
    in
    match src.branches with
    | Some branches ->
        let copy_tree = Hashtbl.create 10 in
        Hashtbl.iter
            (fun tree map ->
                let results =
                    All_sets.IntSetMap.fold
                        (fun intset ntbl acc -> 
                            let copy_node = Hashtbl.create 10 in
                            Hashtbl.iter (add_data copy_node) ntbl;
                            All_sets.IntSetMap.add intset copy_node acc)
                        map
                        All_sets.IntSetMap.empty
                in
                Hashtbl.add copy_tree tree results)
            branches;
        { dest with branches = Some copy_tree; }
    | None -> dest

(* the converse of the above function *)
let convert_static_to_dynamic_branches ~src ~dest = 
    let add_data ctbl ntbl = 
        All_sets.IntSetMap.iter
            (fun set new_code ->
                let new_name = Hashtbl.find dest.character_codes new_code
                and old_name,rep_code = 
                    let rep_code = All_sets.Integers.max_elt set in
                    Hashtbl.find src.character_codes rep_code,rep_code
                in
                let lengths = 
                    try Hashtbl.find ctbl old_name
                    with | Not_found ->
                        failwithf "Could not find old_name:%s; new_code:%d new_name:%s; rep_code:%d"
                                    old_name new_code new_name rep_code
                in
                Hashtbl.add ntbl new_name lengths)
            src.static_dynamic_codes;
        ntbl
    in
    match src.branches with
    | Some branches ->
        let copy_tree = Hashtbl.create 10 in
        Hashtbl.iter
            (fun tree map -> 
                let results = 
                    All_sets.IntSetMap.fold
                        (fun partition lengths intset -> 
                            let ntbl = add_data lengths (Hashtbl.create 17) in
                            All_sets.IntSetMap.add partition ntbl intset)
                        map
                        All_sets.IntSetMap.empty
                in
                Hashtbl.add copy_tree tree results)
            branches;
        {dest with branches = Some copy_tree; }
    | None -> dest

let set_dyna_data seq_arr  = {seq_arr = seq_arr}

(** [get_recost pams] returns the rearrangement cost in [pams] *)
let get_recost user_pams = 
    match user_pams.re_meth with
    | None -> failwith "The rearrangement cost is not specified"
    | Some (`Locus_Breakpoint c)
    | Some (`Locus_Inversion c) -> c

(** [get_locus_indel_cost user_pams] returns the locus indel cost in [pams] *)
let get_locus_indel_cost user_pams = 
    match user_pams.locus_indel_cost with
    | None -> failwith "The locus indel cost is not specified"
    | Some c -> c

let get_character_set_name data codes : string option = match codes with
    | [] -> failwith "No characters specified"
    | code::codes -> 
        let char = Hashtbl.find data.character_codes code in
        if Hashtbl.mem data.character_nsets char then
            Some (Hashtbl.find data.character_nsets char)
        else
            None

let modified_characters data_one data_two : int =
    let data_one = data_one.character_specs
    and data_two = data_two.character_specs in
    let data_one,data_two =
        if (Hashtbl.length data_one) < (Hashtbl.length data_two)
            then data_two,data_one
            else data_one,data_two
    and modified = ref 0 in
    Hashtbl.iter
        (fun k v1 ->
            try let v2 = Hashtbl.find data_two k in
                if v2 <> v1 then incr modified else ()
            with | Not_found -> let () = incr modified in ())
        data_one;
    !modified

let get_likelihood_model data chars = 
    let get_model x = match Hashtbl.find data.character_specs x with
        | Static dat ->
            begin match dat.Nexus.File.st_type with
                | Nexus.File.STLikelihood model -> model
                | _ -> failwith "unsupported static character"
            end
        | Dynamic s when s.state = `Ml ->
            begin match s.lk_model with
                | Some x -> x
                | None   -> failwith "inconsistent dynamic likelihood state"
            end
        | _ -> failwith "unsupported characters"
    in
    match List.map get_model chars with
    | h::t ->
        if List.fold_left ~f:(fun acc x -> acc && (x = h)) ~init:true t then
            h
        else 
            failwith "Inconsistent Model over characters"
    | []   ->
        failwith "No Characters found"

let is_fs data code =
    match Hashtbl.find data.character_specs code with
    | Kolmogorov _ | Set | Static _ -> false
    | Dynamic x -> match x.initial_assignment with
        | `FS _ -> true
        | `Partitioned _
        | `AutoPartitioned _
        | `DO 
        | `GeneralNonAdd -> false

let get_empty_seq alph = 
    let seq = Sequence.create 1 in
    let seq = Sequence.prepend_char seq (Alphabet.get_gap alph) in
    { seq = seq; delimiter = []; code = -1; }


let print (data : d) = 
    let print_taxon (key : int) (ch_ls : (int, cs) Hashtbl.t) = 
        let len = Hashtbl.fold (fun _ _ acc -> acc + 1) ch_ls 0 in
        let  taxa_name = All_sets.IntegerMap.find key data.taxon_codes in  
        Printf.fprintf stdout "Taxon: %s, number chars: %i\n" taxa_name len;
        Hashtbl.iter
            (fun _ ch ->
                 match ch with 
                 | Dyna (code, dyna_data), _ ->         
                       let  char_name = Hashtbl.find data.character_codes code in 
                       Printf.fprintf stdout "%s -> " char_name; 
                       Array.iter (fun seq -> 
                                       Printf.fprintf stdout "%d:" seq.code;
                                       Sequence.print stdout seq.seq Alphabet.nucleotides;
                                       Printf.fprintf stdout " | ")
                                  dyna_data.seq_arr;
                 | Stat (code, None), _ ->
                       let a = match Hashtbl.find data.character_specs code with
                            | Static x -> x.Nexus.File.st_alph
                            | _ -> failwith "Nope"
                       in
                       Printf.fprintf stdout "[%d]%s |" code
                            (Alphabet.match_code (Alphabet.get_gap a) a)
                 | Stat (code, (Some stuff)), _ ->
                       let a =
                           try match Hashtbl.find data.character_specs code with
                                | Static x -> x.Nexus.File.st_alph
                                | _ -> failwith "Nope"
                            with | _ -> failwithf "Couldn't find %d in specs" code
                       in
                       begin match Nexus.File.static_state_to_list stuff with 
                            | []  -> 
                                Printf.fprintf stdout "[%d]%s |" code
                                        (Alphabet.match_code (Alphabet.get_gap a) a)
                            | [x] -> 
                                Printf.fprintf stdout "[%d]%s |" code
                                        (Alphabet.match_code x a)
                            | xs  ->
                                Printf.fprintf stdout "[%d](" code;
                                List.iter 
                                    ~f:(fun x -> Printf.fprintf stdout "%s" 
                                                (Alphabet.match_code x a))
                                    xs;
                                Printf.fprintf stdout ")"
                       end)
            ch_ls;
            print_newline ()
    and print_branches tree_name set = 
        Printf.printf "Tree: [%s]\n" tree_name;
        All_sets.IntSetMap.iter
            (fun intset code_tbl ->
                Printf.printf "\t[";
                All_sets.Integers.iter (Printf.printf "%d, ") intset;
                Printf.printf "]:";
                Hashtbl.iter 
                    (fun n f -> Printf.printf "\t\t%s -- %f\n" n f)
                    code_tbl)
            set
   and print_models chars = 
        let get_function code = 
            match Hashtbl.find data.character_specs code with
            | Static spec ->
                begin match spec.Nexus.File.st_type with
                | Nexus.File.STLikelihood x -> x.MlModel.spec
                | _ -> assert false
                end
            | Dynamic ({state = s} as x) when s = `Ml ->
                begin match x.lk_model with
                | Some m -> m.MlModel.spec
                | _ -> assert false
                end
            | _ -> assert false
        in
        List.iter
            (fun x -> 
                let m = get_likelihood_model data x in
                Printf.printf "%d:" (List.hd x);
                MlModel.output_model (print_string) `Nexus m None;
                print_newline ())
            (MlModel.categorize_by_model get_function chars)
    and print_specs (code : int) (spec : specs) = 
        let name = Hashtbl.find data.character_codes code in 
        Printf.fprintf stdout "Key: %i, name: %s " code name; 
        (match spec with 
        | Dynamic dspec ->
              (match dspec.state with 
               | `SeqPrealigned -> Printf.fprintf stdout "Seq Prealigned"
               | `Seq -> Printf.fprintf stdout "Seq"
               | `Ml -> Printf.fprintf stdout "Dynamic ML"
               | `Breakinv -> Printf.fprintf stdout "Breakinv"
               | `Chromosome -> Printf.fprintf stdout "Chromosome"
               | `Genome -> Printf.fprintf stdout "Genome"
               | `Annotated -> Printf.fprintf stdout "Annotated")
        | Static {Nexus.File.st_type=st_type} ->
              (match st_type with
               | Nexus.File.STOrdered -> Printf.fprintf stdout "Ordered"
               | Nexus.File.STUnordered -> Printf.fprintf stdout "Unordered"
               | Nexus.File.STLikelihood m -> Printf.fprintf stdout "Static ML"
               | Nexus.File.STSankoff m -> Printf.fprintf stdout "Sankoff")
        | _ -> Printf.fprintf stdout "Not Dynamic");
        print_newline ()
    and print_csets name chars = 
        Printf.printf "%s -- (" name;
        List.iter (fun x -> Printf.printf "%s, " x) chars;
        Printf.printf ")";
        print_newline ()
    and print_intset_codes intset i =
        Printf.printf "(";
        All_sets.Integers.iter (Printf.printf "%d, ") intset;
        Printf.printf ") : %d\n%!" i;
    and print_int_intlist i list =
        Printf.printf "%d: " i;
        List.iter (Printf.printf "%d, ") list;
        print_newline ()
    and print_character_codes key bind =
        Printf.printf "[%d,%s],%!" key bind;
    in
    Printf.fprintf stdout "Number of sequences: %i\n" (List.length data.dynamics);
    List.iter (Printf.fprintf stdout "%i ") data.dynamics; print_newline ();
    Printf.printf "\n check character_codes: \n%!";
    Hashtbl.iter print_character_codes data.character_codes;
    Printf.printf "\n check taxon_characters:\n %!";
    Hashtbl.iter print_taxon data.taxon_characters;
    Printf.printf "\n check search_base_characters:\n %!";
    Hashtbl.iter print_taxon data.searchbase_characters;
    Printf.printf "\n check character_specs:\n%!";
    Hashtbl.iter print_specs data.character_specs;
    Printf.printf "\n check character_sets:\n%!";
    Hashtbl.iter print_csets data.character_sets;
    Printf.printf "\n check models:\n%!";
    print_models (data.dynamics @ data.static_ml);
    let () =
        Printf.printf "\n check branches:\n%!";
        match data.branches with
        | Some d -> Hashtbl.iter print_branches d
        | None   -> Printf.printf "None%!\n"
    in
    Printf.printf "\n Check Dynamic->Static Codes\n%!";
    IntMap.iter print_int_intlist data.dynamic_static_codes;
    Printf.printf "\n Check Static->Dynamic Codes\n%!";
    All_sets.IntSetMap.iter print_intset_codes data.static_dynamic_codes;
    print_newline ()


let get_weight c data = 
    match Hashtbl.find data.character_specs c with
    | Dynamic spec -> spec.weight
    | Static spec -> spec.Nexus.File.st_weight
    | _ -> 1.0

let get_weights data =
    Hashtbl.fold (fun x _ acc -> (x, get_weight x data) :: acc) 
    data.character_specs []

(* Returns a fresh object with the added synonym from [a] to [b] to [data]. *)
let rec aux_add_synonym stack data (a, b) =
    (* We need to check if [a] already has an assigned synonym *)
    if (a = b) then begin
        let msg = String.concat " to " (a :: (List.rev (b :: stack))) in
        output_warningf "I'm@ ignoring@ the@ self-synonym@ %s." msg;
        data
    end else if All_sets.StringMap.mem a data.synonyms then begin
        (* Hum ... so we do have a synonym, we need now to check if [(a, b)]
         * is consistent with whatever is stored currently in [data],
         * otherwise this is an exception *)
        let cur = All_sets.StringMap.find a data.synonyms in
        if cur <> b then begin
            output_errorf "@[Overriding@ synonym:@ %s@ to@ %s@ will@ now@ map@ to@ %s@]" a cur b;
            if All_sets.StringMap.mem b data.synonyms then 
                aux_add_synonym (b :: stack) data 
                (a, All_sets.StringMap.find b data.synonyms)
            else
                let ns = All_sets.StringMap.add a b data.synonyms in
                { data with synonyms = ns }
        end else data
    end else if All_sets.StringMap.mem b data.synonyms then begin
        aux_add_synonym (b :: stack) data (a, All_sets.StringMap.find b data.synonyms)
    end else begin 
        (* [a] doesn't have a synonym, so we simply add it *)
        let ns = All_sets.StringMap.add a b data.synonyms in
        { data with synonyms = ns }
    end

let add_synonym = aux_add_synonym []

let add_synonyms_file data file = 
    try
        let ch, file = FileStream.channel_n_filename file in
        let syns = Parser.Dictionary.of_channel_assoc ch in
        close_in ch;
        let len = List.length syns in
        output_infof "@[The@ file@ %s@ contains@ %d@ synonyms.@]"
                     (StatusCommon.escape file) len;
        List.fold_left ~f:add_synonym ~init:data syns
    with
    | (Sys_error err) as exn ->
        output_errorf
            ("@[Couldn't@ open@ file@ %s@ to@ load@ the@ list@ of@ "^^
             "synonyms.@ @ The@ system@ error@ message@ is@ %s.@]")
            (FileStream.filename file)
            (err);
        raise exn

let warn_if_repeated_and_choose_uniquely list str file =
    let repeated, selected = 
        List.fold_left ~f:(fun (repeated, seen) item ->
            if All_sets.Strings.mem item seen then
                All_sets.Strings.add item repeated, seen
            else repeated, All_sets.Strings.add item seen) 
        ~init:(All_sets.Strings.empty, All_sets.Strings.empty) list
    in
    let total = All_sets.Strings.cardinal repeated in
    if total > 0 then
        if total = 1 then
            let item = All_sets.Strings.choose repeated in
            output_errorf
                "%s@ is@ duplicated@ in@ the@ %s@ %s"
                (StatusCommon.escape item)
                (StatusCommon.escape str)
                (StatusCommon.escape file)
        else begin
            let message, _ = 
                All_sets.Strings.fold 
                    (fun item (str, cnt) ->
                        if cnt = 0 then str ^ item, 1
                        else str ^ ",@ " ^ item, 1) 
                    (repeated)
                    (("The@ following@ items@ are@ duplicated@ in@ " ^
                      "@ the@ " ^ StatusCommon.escape str ^ "@ " ^ StatusCommon.escape file
                     ^":@ "), 0)
            in
            Status.user_message Status.Warning message
        end
    else ();
    All_sets.Strings.fold (fun x acc -> x :: acc) selected []


let trim taxon =
    let rec non_empty_position x n fmod res =
        if x = n then (res x)
        else 
            match taxon.[x] with
            | ' ' | '\t' -> non_empty_position (fmod x) n fmod res
            | _ -> x
    in
    let len = String.length taxon in
    if 0 < len then 
        let start = non_empty_position 0 len succ pred
        and final = non_empty_position (len - 1) (-1) pred succ in
        String.sub taxon start (final - start + 1) 
    else taxon

(* convert branch lengths with subsets of leaves *)
let branches_to_map data initial_table branch_table trees = 
    let complete_set =
        All_sets.StringMap.fold
            (fun k v acc -> All_sets.Integers.add v acc)
            data.taxon_names
            (All_sets.Integers.empty)
    and new_tree_table = match initial_table with
        | Some x -> x
        | None -> Hashtbl.create 13
    and found = ref false in
    let create_complement to_remove =
        All_sets.Integers.fold
            (All_sets.Integers.remove)
            to_remove
            complete_set
    and create_all_char_table dist = 
        let new_internal = Hashtbl.create 1229 in
        Hashtbl.iter
            (fun _ v -> Hashtbl.add new_internal v dist)
            data.character_codes;
        new_internal
    (* add two table lengths, to resolve rooted trees *)
    and add_two_lengths tbl1 tbl2 =
        Hashtbl.iter
            (fun k dist1 ->
                let dist2 =
                    try Hashtbl.find tbl2 k
                    with | Not_found -> failwith "Sets don't match"
                in
                Hashtbl.replace tbl1 k (dist1 +. dist2))
            tbl1
    in
    (* assoc partition with a value, if it exists, then add them, as it is
     * due to a rooted tree *)
    let add_single (t_name:string) (set:All_sets.Integers.t)
                   (chartbl_opt: (string,float) Hashtbl.t option)
                   (partition_map: 'a All_sets.IntSetMap.t): 'a All_sets.IntSetMap.t =
        match chartbl_opt with
        | Some n ->
            let comp = create_complement set in
            if (All_sets.IntSetMap.mem set partition_map)
                or (All_sets.IntSetMap.mem comp partition_map)
                then begin
                    let t = All_sets.IntSetMap.find set partition_map in
                    add_two_lengths t n;
                    partition_map -->
                        All_sets.IntSetMap.add comp t -->
                        All_sets.IntSetMap.add set t
                end else begin
                    partition_map -->
                        All_sets.IntSetMap.add comp n -->
                        All_sets.IntSetMap.add set n
                end
        | None -> partition_map
    in
    let rec continually_add t_name f_ext mapp t_node =
        match t_node with
        | Tree.Parse.Leafp (x,dat) ->
            let single_set =
                try
                    let t_code = All_sets.StringMap.find x data.taxon_names in
                    All_sets.Integers.add t_code (All_sets.Integers.empty)
                with Not_found -> 
                    output_errorf "While parsing tree, I cannot find taxon name %s in read data." x;
                    raise Not_found
            in
            (single_set,add_single t_name single_set (f_ext dat) mapp)
        | Tree.Parse.Nodep (xs,(_,dat)) ->
            let children_set,children_map =
                List.fold_left
                    ~f:(fun (oths,othm) nd ->
                         let (one_set,one_map) = continually_add t_name f_ext othm nd in
                         (All_sets.Integers.union oths one_set, one_map))
                    ~init:(All_sets.Integers.empty,mapp)
                    xs
            in
            children_set,add_single t_name children_set (f_ext dat) children_map
    in 
    let add_single_tree (t_name_opt,t) =
        let tree_name =
            match t_name_opt with | Some n -> String.uppercase n | None -> ""
        in
        let rec process_tree acc_map x = match x with
            | Tree.Parse.Branches t ->
                found := true;
                let _,map = continually_add
                    tree_name
                    (fun x -> match (x:float option) with
                        | Some x -> Some (create_all_char_table x)
                        | None -> None)
                    acc_map
                    t
                in
                map
            | Tree.Parse.Characters t ->
                found := true;
                let old_table = match branch_table with
                    | None -> failwith "Cannot find defined table"
                    | Some branch_table -> 
                        try Hashtbl.find branch_table tree_name
                        with | Not_found ->
                            failwithf "Cannot find pre-defined tree, %s." tree_name
                in
                let _,map =
                    continually_add
                        tree_name
                        (fun x -> match (x:string option) with
                            | Some nd -> Some (Hashtbl.find old_table (String.uppercase nd))
                            | None -> None)
                        acc_map
                        t
                in 
                map
            | Tree.Parse.Annotated (t,_) -> process_tree acc_map t
            | Tree.Parse.Flat _ -> 
                acc_map
        in
        let tree_map =
            List.fold_left ~f:process_tree ~init:(All_sets.IntSetMap.empty) t
        in
        Hashtbl.add new_tree_table tree_name tree_map
    in
    List.iter add_single_tree trees;
    new_tree_table,!found


let verify_trees data (((name,tree), file, position) : parsed_trees) =
    let esc_file = StatusCommon.escape file in
    let rec leaves acc tree = 
        let rec aux_leaves f acc subtree = match subtree with
            | Tree.Parse.Nodep (c, _) ->
                List.fold_left ~f:(aux_leaves f) ~init:acc c
            | Tree.Parse.Leafp x -> (f x) :: acc
        in
        match tree with
            | Tree.Parse.Flat t          -> aux_leaves (fun x -> x) acc t
            | Tree.Parse.Characters t    -> aux_leaves (fun x -> fst x) acc t
            | Tree.Parse.Branches t      -> aux_leaves (fun x -> fst x) acc t
            | Tree.Parse.Annotated (t,_) -> leaves acc t
    in
    let rec stop_if_not_all_terminals_in_tree map taxon =
        let taxon = trim taxon in
        if All_sets.StringMap.mem taxon data.synonyms then
            stop_if_not_all_terminals_in_tree map
            (All_sets.StringMap.find taxon data.synonyms)
        else 
            if All_sets.StringMap.mem taxon map then begin
                All_sets.StringMap.remove taxon map
            end else 
                let msg = 
                    ("input@ tree@ " ^ string_of_int position ^ 
                    (if "" <> file then "@ of@ file@ " ^ esc_file else "") ^
                    "@ has@ the@ terminal@ "
                    ^ StatusCommon.escape taxon ^ 
                    "@ and@ there@ is@ no@ data@ loaded@ for@ it")
                in
                let () = Status.user_message Status.Error msg in
                failwith "Data not found"
                    
    in
    let leafs = List.fold_left ~f:leaves ~init:[] tree in
    ignore
        (warn_if_repeated_and_choose_uniquely leafs 
            ("input@ tree@ " ^ string_of_int position ^ "@ of@ file@ ") file);
    let res = 
        List.fold_left ~f:(stop_if_not_all_terminals_in_tree)
                       ~init:data.taxon_names
                        leafs
    in
    if All_sets.StringMap.is_empty res then ()
    else 
        let taxa = 
            (String.concat ", "
            (List.map StatusCommon.escape 
                (All_sets.StringMap.fold (fun a _ acc -> a :: acc) res [])))
        in
        let file_string =
            if "" <> file then "@ of@ file@ " ^ esc_file else ""
        in
        let msg = 
            "The@ following@ terminals@ do@ not@ appear@ in@ the@ input@ "
            ^ "tree@ " ^ string_of_int position ^ file_string ^ "@ :@ " ^ taxa ^
            ".@ Beware@ that@ this@ tree@ will@ be@ incompatible@ with@ any@ " ^
            "@ other@ trees@ built@ by@ POY,@ as@ some@ terminals@ appearing@ "
            ^ "in@ the@ new@ trees@ will@ be@ missing@ on@ this@ one@ and@ " ^ 
            "could@ cause@ errors." 
        in
        Status.user_message Status.Warning msg

let process_trees data file =
    try
        let ch, file = FileStream.channel_n_filename file in
        let trees = Tree.Parse.of_channel ch in
        let () = close_in ch in
        let len = List.length trees in
        let msg = 
            "@[The@ file@ " ^ StatusCommon.escape file ^ 
            "@ contains@ " ^ string_of_int len ^ "@ trees.@]"
        in
        let cnt = ref 0 in
        let trees = List.map ~f:(fun x -> incr cnt; (None,x), file, !cnt) trees in
        let branches, found =
            branches_to_map data None None (List.map (fun (x,_,_) -> x) trees)
        in
        let branches = if found then Some branches else None in
        { data with trees = data.trees @ trees;
                    branches = branches;
                    iterate_branches = (not found);}
    with
    | Sys_error err ->
        let file = FileStream.filename file in
        let msg = "@[Couldn't@ open@ file@ " ^ file ^ "@ to@ load@ the@ " ^
                  "trees.@ @ The@ system@ error@ message@ is@ " ^ err ^ ".@]" in
        output_error msg;
        data

let process_fixed_states data = function
    | Some file ->
        begin try
            let ch, file = FileStream.channel_n_filename file in
            let fs = 
                Parser.FixedStatesDict.of_channel 
                FileContents.Nucleic_Acids ch 
            in
            close_in ch;
            let len = List.length fs in
            let msg =
                Printf.sprintf "@[The@ file@ %s@ defines@ %d@ valid@ states.@]"
                                (StatusCommon.escape file) len
            in
            Status.user_message Status.Information msg;
            { 
                data with do_fs = true;
                current_fs = fs;
                current_fs_file = file; 
            }
        with | Sys_error err ->
            let file = FileStream.filename file in
            let msg =
                Printf.sprintf
                    ("@[Couldn't@ open@ file@ %s@ to@ load@ the@ fixed@ "^^
                     "states.@ The@ system@ error@ message@ is@ %s.@]")
                    file err
            in
            output_error msg;
            data
        end
    | None ->
        { data with do_fs = true; }

let find_code_for_root_if_removed data =
    (* We want to test, if a terminal that is currently root is removed, we choose
    * the one with the lowest terminal code, but if that's not the root, we
    * continue like nothing *)
    match data.root_at with
    | Some c -> 
        if Hashtbl.mem data.taxon_characters c
        then 
            data
        else begin
            let nc = 
                Hashtbl.fold 
                    (fun c _ acc -> match acc with
                        | None -> Some c
                        | Some accc ->
                            if c < accc then Some c
                            else acc)
                    data.taxon_characters None
            in
            { data with root_at = nc }
        end
    | None -> data

let pick_code_for_root code data = match data.root_at with
    | None ->
        { data with root_at = Some code }
    | Some previous -> 
        if Hashtbl.mem data.taxon_characters previous then
            data 
        else 
            { data with root_at = Some code }

(** Returns a [d] with the following condition: if [taxon] is present in [data]
 * then [data] is returned, otherwise, [data] with the added name and assigned code
 * is returned. In addition, if there is no preferred taxon assigned to the
 * data structure, we set whatever is being load, that way we ensure that if no
 * root is defined, the first loaded taxon will be it by default. *)
let rec process_taxon_code data taxon file =
    (* Check first if the taxon is in the current list of taxa *)
    let taxon = trim taxon in
    if All_sets.StringMap.mem taxon data.taxon_names then
        (*already exists, we must have a root already! *)
        let new_taxon_file =
            let set = All_sets.StringMap.find taxon data.taxon_files in
            let set = All_sets.Strings.add file set in
            All_sets.StringMap.add taxon set data.taxon_files
        in
        { data with taxon_files = new_taxon_file },
        All_sets.StringMap.find taxon data.taxon_names
    else if All_sets.StringMap.mem taxon data.synonyms then
        (* Is a synonym *)
        process_taxon_code data
        (All_sets.StringMap.find taxon data.synonyms) file
    else begin
        (* It is new, so lets assign it a code and add it *)
        let code = data.number_of_taxa + 1 in
        let taxon_names =
            All_sets.StringMap.add taxon code data.taxon_names
        in
        let taxon_codes =
            All_sets.IntegerMap.add code taxon data.taxon_codes
        in
        let taxon_files =
            All_sets.StringMap.add taxon (All_sets.Strings.singleton file)
            data.taxon_files
        in
        { data with
            number_of_taxa = data.number_of_taxa + 1;
            taxon_names = taxon_names;
            taxon_codes = taxon_codes;
            taxon_files = taxon_files }, code
    end

let rec process_searchbase_code data taxon file =
    let taxon = trim taxon in
    if debug_search_base then 
        Printf.printf "process_searchbase_code with taxon=%s\n%!" taxon;
    if All_sets.StringMap.mem taxon data.taxon_names then 
        (*already exists, we must have a root already! *)
        (*let new_taxon_file = 
            let set = All_sets.StringMap.find taxon data.taxon_files in
            let set = All_sets.Strings.add file set in
            All_sets.StringMap.add taxon set data.taxon_files 
        in*)
        let new_searchbase_file =
        if All_sets.StringMap.mem taxon data.searchbase_files then   
            let set = All_sets.StringMap.find taxon data.searchbase_files in
            let set = All_sets.Strings.add file set in
            All_sets.StringMap.add taxon set data.searchbase_files
        else
            All_sets.StringMap.add taxon (All_sets.Strings.singleton file)
            data.searchbase_files
        in
        { data with searchbase_files = new_searchbase_file }, 
        All_sets.StringMap.find taxon data.taxon_names
    else if All_sets.StringMap.mem taxon data.synonyms then
        (* Is a synonym *)
        process_searchbase_code data 
        (All_sets.StringMap.find taxon data.synonyms) file
    else begin
        let code = data.number_of_taxa + 1 in
        let searchbase_names = 
            All_sets.StringMap.add taxon code data.searchbase_names
        in
        let searchbase_codes = 
            All_sets.IntegerMap.add code taxon data.searchbase_codes
        in
        let searchbase_files =
            All_sets.StringMap.add taxon (All_sets.Strings.singleton file)
            data.searchbase_files 
        in
        let data = pick_code_for_root code data in
        { data with
            number_of_taxa = data.number_of_taxa + 1;
            searchbase_names = searchbase_names; searchbase_codes = searchbase_codes;
            searchbase_files = searchbase_files },
        code
    end


let get_taxon_characters data tcode =
    try Hashtbl.find data.taxon_characters tcode with 
    | _ -> create_ht ()

let get_searchbase_characters data tcode =
    try Hashtbl.find data.searchbase_characters tcode with 
    | _ -> create_ht ()

(* Changes in place *)
let add_static_character_spec data (code, spec) =
    if not spec.Nexus.File.st_eliminate then begin
        Hashtbl.replace data.character_specs code (Static spec);
        Hashtbl.replace data.character_names spec.Nexus.File.st_name code;
        Hashtbl.replace data.character_codes code spec.Nexus.File.st_name;
    end else ()

let report_static_input file f_out =
    let characters = Array.length f_out.Nexus.File.characters 
    and unaligned = List.length (f_out.Nexus.File.unaligned)
    and tree = List.length f_out.Nexus.File.trees 
    and taxa = Array.length f_out.Nexus.File.taxa in
    let msg =
        "@[The@ file@ " ^ StatusCommon.escape file ^ "@ defines@ " ^
        string_of_int characters ^ "@ static@ homology@ characters,@ " ^
        string_of_int unaligned ^ " unaligned@ sequences,@ and@ " ^
        string_of_int tree ^ "@ trees,@ containing@ " ^ 
        string_of_int taxa ^ "@ taxa.@]"
    in
    Status.user_message Status.Information msg

let repack_codes data = 
    (* We verify that all the terminals have consecutive codes, and that all the
    * characters have consecutive codes. Finally we restore the code generation
    * reference to the necessary code *)
    (* We first collect the data about all the taxa that we have and build a map
    * of codes, then we do the same about all the characters that we have *)
    let taxa_used_codes = ref []
    and char_used_codes = ref [] in
    All_sets.StringMap.iter (fun name code ->
        taxa_used_codes := code :: !taxa_used_codes) data.taxon_names;
    Hashtbl.iter (fun _ code ->
        char_used_codes := code :: !char_used_codes) data.character_names;
    (* We now sort them *)
    let taxa_used_codes = List.sort ( - ) !taxa_used_codes
    and char_used_codes = List.sort ( - ) !char_used_codes in
    let rec available_codes_generator acc inverted cnt to_check =
        match to_check with
        | h :: t ->
                if h = cnt then 
                    available_codes_generator acc (max h inverted)
                    (cnt + 1) t
                else if cnt < h then 
                    available_codes_generator (cnt :: acc) inverted 
                    (cnt + 1) to_check
                else assert false
        | [] -> (acc, inverted, cnt)
    in
    let available_taxa_codes, greatest_taxa_code, taxa_counter = 
        match taxa_used_codes with
        | h :: t ->
                available_codes_generator [] 0 h taxa_used_codes
        | [] -> [], 0, 0
    in
    let available_char_codes, greatest_char_code, character_counter  = 
        match char_used_codes with
        | h :: t -> available_codes_generator [] 0 h char_used_codes
        | [] -> [], 0, 0
    in
    let rec process_available recode_function check data used_code available_codes =
        match available_codes with
        | av :: ta ->
                if check data used_code then
                    let uc = used_code in
                    let data = recode_function data uc av in
                    process_available recode_function check data (uc - 1) ta
                else 
                    process_available recode_function check data (used_code - 1)
                    available_codes
        | [] -> data
    in
    let recode_taxon data used_code available_code =
        (* Fix taxon names *)
        let taxon_codes, name =
            let name = All_sets.IntegerMap.find used_code data.taxon_codes in
            let taxon_codes = All_sets.IntegerMap.remove used_code
            data.taxon_codes in
            All_sets.IntegerMap.add available_code name taxon_codes, name
        in
        let taxon_names = 
            let taxon_names = All_sets.StringMap.remove name data.taxon_names
            in
            All_sets.StringMap.add name available_code taxon_names
        in
        let () =
            try 
                let taxon_characters = Hashtbl.find data.taxon_characters used_code in
                Hashtbl.remove data.taxon_characters used_code;
                Hashtbl.replace data.taxon_characters available_code taxon_characters;
            with
            | Not_found -> ()
        in
        let root_at = match data.root_at with
            | None -> None
            | Some x when x = used_code -> 
                    Some available_code
            | x -> x
        in
        { data with taxon_codes = taxon_codes; 
                    taxon_names = taxon_names;
                    root_at = root_at }
    in
    let recode_character data used_code available_code =
        let () = 
            try let name = Hashtbl.find data.character_codes used_code in
                Hashtbl.remove data.character_names name;
                Hashtbl.remove data.character_codes used_code;
                Hashtbl.replace data.character_codes available_code name;
                Hashtbl.replace data.character_names name available_code;
            with | Not_found -> ()
        in
        let () = 
            try let spec = Hashtbl.find data.character_specs used_code in
                Hashtbl.replace data.character_specs available_code spec;
                Hashtbl.remove data.character_specs used_code;
            with | Not_found -> ()
        in
        Hashtbl.iter 
            (fun code chars ->
                try let char = Hashtbl.find chars used_code in
                    Hashtbl.remove chars used_code;
                    Hashtbl.replace chars available_code 
                    (match char with
                    | (Dyna (_, d)), x -> 
                            ((Dyna (available_code, d)), x)
                    | (Stat (_, d)), x ->
                            ((Stat (available_code, d)), x));
                with | Not_found -> ())
            data.taxon_characters;
        let dyn_stat_codes =
            All_sets.IntegerMap.map
                (fun lst ->
                    let tl = List.filter (fun x -> not (x = used_code)) lst in
                    available_code :: tl)
                data.dynamic_static_codes
        in
        {data with dynamic_static_codes = dyn_stat_codes;
                   static_dynamic_codes = reverse_dynamic_static_codes dyn_stat_codes; }
    in
    let data = 
        let check data code = 
            All_sets.IntegerMap.mem code data.taxon_codes
        in
        process_available recode_taxon check data greatest_taxa_code
                            available_taxa_codes 
    in
    let check data code = Hashtbl.mem data.character_codes code in
    process_available recode_character check data greatest_char_code 
                        available_char_codes



let add_character data file tcmfile tcm tcm3 default_mode lk_model alphabet dyna_state dyna_pam weight = 
        
        incr data.character_code_gen;
        let chcode = !(data.character_code_gen) in
        let dspec = {
            filename = file;
            fs = data.current_fs_file;
            tcm = tcmfile;
            initial_assignment = default_mode;
            tcm2d = tcm;
            tcm3d = tcm3;
            lk_model = lk_model;
            alph = alphabet;
            state = dyna_state;
            pam = dyna_pam;
            weight = weight; 
            polymorphism = `Do_All; }
        in
        if debug_parsed_seq then
        Printf.printf "add character with character code=%d,file=%s\n%!" chcode file;
        Hashtbl.replace data.character_specs chcode (Dynamic dspec);
        Hashtbl.replace data.character_names file chcode;
        Hashtbl.replace data.character_codes chcode file;
        chcode, data

let process_parsed_breakinv data res
original_filename file tcmfile tcm tcm3 default_mode lk_model alphabet dyna_state dyna_pam weight =
    if debug_parsed_seq then Printf.printf "process parsed breakinv\n%!";
    let get_multichrom_and_delimiter file_taxon_chrom_loci_frag_seq =
        (*breakinv could be multi-chromosome. get the seq and its length out =>
        * [
             [(chrom1_seq of taxon1,seqlen);(chrom2_seq of taxon1,seqlen);....];
             [(chrom1_seq of taxon2,seqlen);(chrom2_seq of taxon2,seqlen);....];
             ....
          ]
        *)
        List.map (fun (taxon_chrom_loci_frag_seq,taxon) ->
            let extract_seq_and_deli (accseq,accdeli) chrom_loci_frag_seq =
                match chrom_loci_frag_seq with 
                | [[seq]] -> 
                        let seq = (*ignore first gap*)
                        Sequence.sub seq 1 ((Sequence.length seq)-1) in
                        (seq::accseq,(Sequence.length seq)::accdeli)
                | _ -> failwith ("we are doing multi-chromosome for breakinv, there\
                should not be any loci('|') or fragment('#') delimiters")
            in
            let seqlst,delilst = List.fold_left ~f:extract_seq_and_deli ~init:([],[]) taxon_chrom_loci_frag_seq in
            (Sequence.concat (List.rev seqlst),List.rev delilst,taxon)
        ) file_taxon_chrom_loci_frag_seq
    in
    let parsed_bkinv_seq = get_multichrom_and_delimiter res in
    if debug_parsed_seq then begin
        Printf.printf "parsed bkinv seq = \n%!";
    List.iter(fun (seq,delilst,taxon) ->
        Sequence.printseqcode seq;
        Printf.printf " deli= %!";
        Utl.printIntArr (Array.of_list delilst);
    ) parsed_bkinv_seq;
    end;
    let process_a_taxon (data,chcode) (singleseq,delilst,taxon) =
        let data, tcode =
            process_taxon_code data taxon original_filename
        in
        let tl = get_taxon_characters data tcode in
        let seqa = 
            {seq = singleseq; 
            delimiter=
                if (List.length delilst)>1 then delilst 
                else []
            ; 
            code= -1} in
        let dyna_data = { seq_arr = [|seqa|]; } in
        let spc = `Specified in
        if debug_parsed_seq then  
            Printf.printf "replace Data.taxon_characters with tcode=%d chcode=%d\n%!" tcode chcode;
        Hashtbl.replace tl chcode (Dyna (chcode,dyna_data),spc);
        Hashtbl.replace data.taxon_characters tcode tl;
        data,chcode
    in
    let chcode,newdata = add_character data file 
    tcmfile tcm tcm3 default_mode lk_model alphabet dyna_state dyna_pam weight in
    let resdata,_ = List.fold_left ~f:process_a_taxon ~init:(newdata,chcode)
    parsed_bkinv_seq in
    resdata
    



let process_parsed_genome data res
original_filename file tcmfile tcm tcm3 default_mode lk_model alphabet dyna_state dyna_pam weight =
    let get_multichromseq file_taxon_chrom_loci_frag_seq =
        (*for genome data, each taxon is a genome. genome might be
        * multi-chromosome, or not, this function works for both.*)
        List.map (fun (taxon_chrom_loci_frag_seq,b) ->
            List.map (fun chrom_loci_frag_seq ->
                match chrom_loci_frag_seq with 
                | [[seq]] -> (seq,b)
                | _ -> failwith ("we are doing multi-chromosome here, there\
                should not be any loci('|') or fragment('#') delimiters")
            ) taxon_chrom_loci_frag_seq
        ) file_taxon_chrom_loci_frag_seq
    in
    let parsed_genome_seq = get_multichromseq res in
    let num_genome = List.length parsed_genome_seq in
    let num_chrom = List.length (List.hd parsed_genome_seq) in
    if debug_parsed_seq then  
        Printf.printf "process_genome num_genome=%d,num_chrom=%d\n%!" num_genome num_chrom;
    let max_chrom_len = ref 0 in
    let genome_arr =
        Array.init num_genome
            (fun ti ->
                Array.init num_chrom
                    (fun ci ->
                        let chrom, name = List.nth (List.nth parsed_genome_seq ti) ci in
                        max_chrom_len := max !max_chrom_len (Sequence.length chrom);
                        (chrom, name)))
    in
    let chcode, data =  add_character data file tcmfile tcm tcm3 default_mode lk_model alphabet dyna_state dyna_pam weight in
    Array.fold_left
        ~f:(fun data chrom_arr ->
            let chrom_ls =
                Array.fold_left
                    ~f:(fun chrom_ls (chrom, _) ->
                        let chrom_len = Sequence.length chrom in
                        match chrom_len with
                        | 1 -> chrom_ls
                        | _ ->
                            let clean_chrom = Sequence.sub chrom 1 (chrom_len - 1) in
                            incr data.seg_code_gen;
                            let code = data.seg_code_gen in
                            if debug_parsed_seq then begin
                                Printf.printf "add to chrom list: code=%d,seq=%!" !code;
                                Sequence.printseqcode clean_chrom;
                            end;
                            {seq = clean_chrom; delimiter = []; code = !code}::chrom_ls)
                    ~init:[]
                    chrom_arr
            in
            let genome_data = {seq_arr = Array.of_list (List.rev chrom_ls)} in
            let _, taxon_name = chrom_arr.(0) in
            let data, tcode =
                process_taxon_code data taxon_name original_filename
            in
            let tl = get_taxon_characters data tcode in
            if debug_parsed_seq then 
                Printf.printf "replace taxon_characters ,chcode=%d with new genome_data\n%!" chcode;
            Hashtbl.replace tl chcode (Dyna (chcode, genome_data), `Specified);
            Hashtbl.replace data.taxon_characters tcode tl;
            data)
        ~init:data
        genome_arr
        
(* normal sequence, sequence could be divided by fragment "#" *)
let process_parsed_normal_sequence data res 
original_filename tcmfile tcm tcm3 default_mode lk_model alphabet 
dyna_state dyna_pam weight prealigned domerge =
    if debug_parsed_seq then Printf.printf "process normal sequence\n%!";
    let get_multi_segment_seq file_taxon_chrom_loci_frag_seq =
        (* input seq list list list list become a matrix looks like this:
        * [ 
            [frag1 in taxon 1;frag2 in taxon 1;...];
            [frag1 in taxon 2;frag2 in taxon 2;...];
            ...
           ]
        *)
        let max_num_fragment = ref 0 in
        let file_taxon_frag_seq = Array.of_list (List.map (fun (taxon_chrom_loci_frag_seq,t) ->
            (*we are only expecting one chromosome for each taxon*)
            (*we are only expecting one loci for each chromosome*)
            (*each loci is a list of fragment, each fragment is a sequence*)
            let seq_t_arr = match taxon_chrom_loci_frag_seq with 
                | [[loci_frag_seq]] -> (*this is a sequence list*) 
                       if !max_num_fragment < (List.length loci_frag_seq) then
                           begin
                           if !max_num_fragment=0 then
                                max_num_fragment := (List.length loci_frag_seq)
                           else
                               failwith ("we have an inconsistent number of \
                               fragments in sequence from input file:"^original_filename);
                           end
                        else ();
                        Array.of_list ( List.map (fun fragseq -> (fragseq,t) )
                        loci_frag_seq )
                | _ -> failwith ("we are working on normal sequence here, there\
                should not be any loci('|') or chromosome ('@') delimiters")
            in
            seq_t_arr
        ) file_taxon_chrom_loci_frag_seq ) in
        (*tranpose the matrix above, make it looks like this
        * [
            [seq belong to frag1 in taxon1;seq belong to frag1 in taxon2;...];
            [seq belong to frag2 in taxon1;seq belong to frag3 in taxon2;...];
            ...
          ]
        *)
        let num_taxon = Array.length file_taxon_frag_seq in
        let max_num_fragment = !max_num_fragment in
        if debug_parsed_seq then  
            Printf.printf "num_taxon=%d,max_num_fragment=%d\n%!" num_taxon max_num_fragment;
        let frag_taxon_seq = ref [] in
        for j = max_num_fragment - 1 downto 0 do
            let acc = ref [] in
            for i = num_taxon - 1 downto 0 do
                match j < Array.length file_taxon_frag_seq.(i) with
                | true ->
                        (*just add seq in taxon.i belong to fragment.j*) 
                        let seq,t =  file_taxon_frag_seq.(i).(j) in
                        acc := ([|seq|],t) :: !acc
                | false ->
                        (*we have missing data, create empty seq with a gap*)
                        let _,taxon = file_taxon_frag_seq.(i).(0) in
                        let seq = Sequence.create 1 in
                        let seq = 
                          Sequence.prepend_char seq (Alphabet.get_gap alphabet)
                        in
                        acc := ([|seq|],taxon) :: !acc
            done;
            frag_taxon_seq := !acc :: !frag_taxon_seq ;
        done;
        file_taxon_frag_seq,!frag_taxon_seq
    in
    let file_taxon_frag_seq,file_frag_taxon_seq = get_multi_segment_seq res in
    (* Now a function to process one taxon at a time to be
     * folded over the taxa collected in the [file]. *)
    let process_a_taxon (data,chcode) (seqarr, taxon)  = 
        (* Here is where, at the parser level, the fixed
         * states support should be added *)
        let data, tcode =
            process_taxon_code data taxon original_filename
        in
        let tl = get_taxon_characters data tcode in
        let seqa =
            let makeone seqa = {seq=seqa; delimiter=[]; code = -1} in
            match dyna_state with
            | `Ml  when not prealigned -> Array.map makeone seqarr 
            | `Seq when not prealigned -> Array.map makeone seqarr
            | `SeqPrealigned           -> Array.map makeone seqarr
            | `Chromosome | `Genome
            | `Annotated  | `Breakinv 
            | `Seq | `Ml (* when prealigned *) ->
                Array.map (fun x -> makeone (Sequence.del_first_char x)) seqarr
        in
        let dyna_data = { seq_arr = seqa; } in 
        let () = 
            let spc =
                let maxlen =
                    Array.fold_left ~f:(fun acc x ->
                        max acc (Sequence.length x.seq)) ~init:0 seqa
                in
                if 2 <= maxlen 
                    then `Specified else 
                if prealigned && maxlen <> 0
                    then `Specified 
                    else `Unknown
            in
            if debug_parsed_seq then begin
                Printf.printf "replace taxon_characters.taxon code=%d,character code=%d with %!"
            tcode chcode;
            Array.iter (fun seq -> Sequence.printseqcode seq;) seqarr;
            end;
            Hashtbl.replace tl chcode (Dyna (chcode, dyna_data),spc) 
        in
        Hashtbl.replace data.taxon_characters tcode tl;
        data,chcode
    in
    let filename_with_charactorNO = 
        let c = ref (-1) in
        ref (fun () -> incr c; original_filename ^ ":" ^ string_of_int !c)
    in
    (*work on the taxon seq list belong to a fragment *)
    let fold_over_a_taxon_lst olddata lst  =
        (*create new character for each fragment*)
        let chcode, newdata = add_character olddata (!(filename_with_charactorNO) ())
        tcmfile tcm tcm3 default_mode lk_model alphabet dyna_state dyna_pam weight in
        let resdata,_ = List.fold_left ~f:process_a_taxon ~init:(newdata,chcode) lst 
        in
        resdata
    in
    if domerge then
        let file_taxon_seq = Array.map (fun taxon_frag_seq ->
            let _,taxon = taxon_frag_seq.(0) in
            let accseq = ref [] in
            Array.iter (fun (seq,t) ->
                accseq := seq::(!accseq)
            ) taxon_frag_seq;
            Array.of_list (List.rev !accseq),taxon
        ) file_taxon_frag_seq in
        fold_over_a_taxon_lst data (Array.to_list file_taxon_seq)
    else 
        List.fold_left ~f:fold_over_a_taxon_lst ~init:data file_frag_taxon_seq
    
    

let process_annotated_chrom data res 
original_filename file tcmfile tcm tcm3 default_mode lk_model alphabet dyna_state dyna_pam weight =
    let get_annotated_seq file_taxon_chrom_loci_frag_seq =
        (*each file has a list of taxons*)
        List.map (fun (taxon_chrom_loci_frag_seq,b) ->
        (*for annotated chromosome data, each taxon is a chromosome*)
        match taxon_chrom_loci_frag_seq with 
        | [chrom_loci_frag_seq] ->
                (* chromosome might be annotated by '|'*)
                List.map (fun loci_frag_seq ->
                    (*we are not expecting multi fragment here*)
                    match loci_frag_seq with
                    | [fragseq] -> (fragseq,b)
                    | _ -> failwith "we are doing annotated-chromosome here, there\
                should not be any fragment('#') delimiters"
                ) chrom_loci_frag_seq
        | _ -> failwith "we are doing annotated-chromosome here, there\
                should not be any chromosome('@') delimiters"
        ) file_taxon_chrom_loci_frag_seq
    in
    let res = get_annotated_seq res in
    let arr = Array.of_list res in
    let num_taxa = Array.length arr in
    let mat = Array.map Array.of_list arr in
    let num_loci = Array.length mat.(0) in
    if debug_parsed_seq then Printf.printf "process_annotated_chrom,num_loci=%d,num_taxa=%d\n%!"
    num_loci num_taxa;
    let chcode, data = add_character data file tcmfile tcm tcm3 default_mode lk_model alphabet dyna_state dyna_pam weight in
    data.seg_code_gen := 0;
    let rec add_taxon data t = match t >= num_taxa with
        | true -> data
        | false ->
            let _, taxon = mat.(t).(0) in
            let data, tcode =
                process_taxon_code data taxon original_filename
            in
            let rec get_annchrom rev_seq_ls l =
                if l = num_loci then
                    List.rev rev_seq_ls
                else begin
                    let seq, _ = mat.(t).(l) in
                    let len = Sequence.length seq in
                    match len > 1 with
                    | false -> get_annchrom rev_seq_ls (l + 1)
                    | true ->
                        incr data.seg_code_gen;
                        let code = !(data.seg_code_gen) in
                        (** this code is for the negative state *)
                        incr data.seg_code_gen;
                        let clean_seq = Sequence.sub seq 1 (len - 1) in
                        let rev_seq_ls =
                            {seq = clean_seq; delimiter = []; code = code}::rev_seq_ls
                        in
                        if debug_parsed_seq then begin
                        Printf.printf "add seq to res_seq_lst with code=%d,seq = %!"
                        code; Sequence.printseqcode clean_seq;
                        end;
                        get_annchrom rev_seq_ls (l + 1)
                end
            in
            let seq_arr = Array.of_list (get_annchrom [] 0) in
            let chrom_data = {seq_arr =  seq_arr} in
            let tl = get_taxon_characters data tcode in
            if debug_parsed_seq then 
                Printf.printf "replace taxon_characters tcode=%d with new chrom_data (chcode = %d)\n%!" tcode chcode;
            Hashtbl.replace tl chcode (Dyna (chcode, chrom_data), `Specified);
            Hashtbl.replace data.taxon_characters tcode tl;
            add_taxon data (t+1)
    in
    add_taxon data 0
    
    

let process_parsed_sequences prealigned weight tcmfile tcm tcm3 default_mode
                             annotated alphabet file dyna_state data res
                             lk_model dyna_pam =
    let data = duplicate data in
    (* in a file , we have a list of taxon followed by '>taxonname'
     ['>'taxon['@'chrom['|'loci['#'frag'#']'|']'@']]
    * *)
    (*let x=ref 0 and y = ref 0 and z = ref 0 and w = ref 0 in
    Printf.printf "Data.process_parsed_sequences\n%!";
    List.map (fun ((taxon_chrom_loci_frag_seq:Sequence.s list list list),b) ->
        Printf.printf "taxon.%d { \n%!" !x; 
        x:=!x+1; y :=0;
        (*each taxon has a list of chromosome,seperated by '@'*)
        List.map( fun (chrom_loci_frag_seq:Sequence.s list list) ->
            Printf.printf "chrom.%d @ \n%!" !y; 
            y:=!y+1; z := 0;
            (*each chromosome has a list of locus, seperated by '|',[ loci0; loci1; ...,]*)
            List.map ( fun (loci_frag_seq:Sequence.s list) ->
                Printf.printf "loci.%d | \n%!" !z;
                z := !z+1; w := 0;
                (*each loci has a list of fragments, seperated by '#',each fragment is a sequence*)
                List.map (fun (frag_seq:Sequence.s) ->
                    Printf.printf "freg.%d # %!" !w;
                    w := !w +1;
                    Sequence.printseqcode frag_seq;
                    Printf.printf " # %!";
                ) loci_frag_seq;
                Printf.printf "| \n%!";
            ) chrom_loci_frag_seq;
            Printf.printf "@ \n%!";
        ) taxon_chrom_loci_frag_seq;
        Printf.printf "} \n%!";
    ) res;*)
    let original_filename = file in
    let locus_name = 
        let c = ref (-1) in
        ref (fun () -> incr c; file ^ ":" ^ string_of_int !c)
    in
    let dyna_pam = match dyna_pam with
        | Some x -> x
        | None   -> dyna_pam_default
    in
    let dyna_state : dyna_state_t =
        match annotated,dyna_pam.Dyn_pam.mode with
        | true, _    -> `Annotated
        | false,None -> dyna_state
        |  _, Some `Chromosome -> `Chromosome
        |  _, Some `Genome -> `Genome
        |  _, Some `Breakinv -> `Breakinv
    in

    let file = 
        if  annotated || (dyna_state = `Chromosome) then 
            original_filename
        else match default_mode with
        | `FS _ | `DO | `GeneralNonAdd |`AutoPartitioned _ -> (!locus_name) () 
        | `Partitioned _ -> original_filename
    in
    let data = 
        if annotated then 
            process_annotated_chrom data res original_filename file 
            tcmfile tcm tcm3 default_mode lk_model alphabet dyna_state dyna_pam weight
        else if dyna_state = `Genome then 
            process_parsed_genome data res original_filename file 
            tcmfile tcm tcm3 default_mode lk_model alphabet dyna_state dyna_pam weight
        else if dyna_state = `Breakinv then
            process_parsed_breakinv data res original_filename file 
            tcmfile tcm tcm3 default_mode lk_model alphabet dyna_state dyna_pam weight
        else if `DO = default_mode || `GeneralNonAdd = default_mode then
            process_parsed_normal_sequence data res original_filename  
            tcmfile tcm tcm3 default_mode lk_model alphabet dyna_state dyna_pam weight prealigned false
        else 
            process_parsed_normal_sequence data res original_filename 
            tcmfile tcm tcm3 default_mode lk_model alphabet dyna_state dyna_pam weight prealigned true
    in 
    data

(* convert Nexus.File.file_output to Data.d *)
let gen_add_static_parsed_file do_duplicate data file file_out =
    let data = if do_duplicate then duplicate data else data in
    (* A function to report the taxa loading *)
    let len_taxa = Array.length file_out.Nexus.File.taxa in
    let st = Status.create "Loading Characters" (Some len_taxa) "taxa" in
    (* [codes] contain the character speecification for the sequences in
    * this file *)
    let codes = 
        let builder x = 
            incr data.character_code_gen;
            (!(data.character_code_gen)), x
        in
        Array.map ~f:builder file_out.Nexus.File.characters
    in
    let new_codes = Array.map fst codes in
    Status.full_report ~msg:"Storing the character specifications" st;
    (* Now we add the codes to the data *)
    Array.iter ~f:(add_static_character_spec data) codes;
    (* And now a function to add one taxon at a time to the data *)
    let data =
        Array.fold_left
            ~f:(fun data name -> match name with
                | Some name -> fst (process_taxon_code data name file)
                | None      -> data)
            ~init:data
            file_out.Nexus.File.taxa
    in
    for row = 0 to len_taxa - 1 do
        (* We ignore the data because we already processed the names *)
        match file_out.Nexus.File.taxa.(row) with
        | None -> ()
        | Some tname ->
                let _, tcode = process_taxon_code data tname file in
                let tl = get_taxon_characters data tcode in
                let add_character column it =
                    let chcode, spec = codes.(column) in
                    if not spec.Nexus.File.st_eliminate then begin
                        let specified = 
                            match it with
                            | Some _ -> `Specified
                            | None -> `Unknown
                        in
                        Hashtbl.replace tl chcode ((Stat (chcode, it)), specified);
                    end;
                    (column + 1)
                in
                if Array.length (file_out.Nexus.File.matrix) > 0 then
                    ignore
                        (Array.fold_left ~f:add_character 
                                         ~init:0
                                         file_out.Nexus.File.matrix.(row));
                Hashtbl.replace data.taxon_characters tcode tl;
                let did = Status.get_achieved st in
                Status.full_report ~adv:(did + 1) st;
    done;
    (* We add the trees *)
    let data = 
        let cnt = ref 0 in
        let trees = 
            List.rev_map
                ~f:(fun x -> (x, file, (incr cnt; !cnt)))
                file_out.Nexus.File.trees
        in
        { data with trees = data.trees @ trees }
    in
    (* combine the branch lengths of the root if it's a rooted tree *)
    let () =
        let unroot_branch_lengths table tree node1 node2 =
            try
                let t_tbl = Hashtbl.find table (String.uppercase tree) in
                let n_1 = Hashtbl.find t_tbl (String.uppercase node1)
                and n_2 = Hashtbl.find t_tbl (String.uppercase node2) in
                Hashtbl.iter 
                    (fun c_name length ->
                        let n_length = length +. (Hashtbl.find n_2 c_name) in
                        let () = Hashtbl.replace n_2 c_name n_length
                        and () = Hashtbl.replace n_1 c_name n_length in
                        () )
                    n_1
            with | Not_found -> failwith ("Cannot find tree "^tree^
                                            " or nodes ("^node1^","^node2^")")
        in
        let get_stuff = function | Tree.Parse.Nodep (_,d) | Tree.Parse.Leafp d -> snd d in
        let combine_on_one ((name,trees),_,_) =
            let rec deal_with_tree name t = match t with
                | Tree.Parse.Flat t -> ()
                | Tree.Parse.Annotated (t,_) -> deal_with_tree name t
                | Tree.Parse.Branches t -> ()
                | Tree.Parse.Characters t -> 
                    begin match t with
                        | Tree.Parse.Leafp _ -> ()
                        | Tree.Parse.Nodep (lst,_) -> 
                            begin match lst with
                                | [l1;l2] ->
                                    begin match (get_stuff l1),(get_stuff l2) with
                                         | Some x,Some y ->
                                            unroot_branch_lengths 
                                                file_out.Nexus.File.branches name x y
                                        | _ -> ()
                                    end
                                | _ -> ()
                            end
                    end
            in
            match name with
            | Some name -> List.iter (deal_with_tree name) trees
            | None      -> ()
        in
        List.iter (combine_on_one) data.trees
    in
    (* Now time to add the molecular sequences *)
    let data = 
        let single_sequence_adder data (u : Nexus.File.unaligned) =
            let alph = u.Nexus.File.u_alph in
            let size = Alphabet.distinct_size (Alphabet.to_sequential alph) in
            let all_elements =
                if alph = Alphabet.nucleotides then
                    31
                else if Alphabet.is_aminoacids alph then
                    21
                else if Alphabet.is_aminoacids alph then
                    21
                else 
                    ~-1
            in
            let seq = u.Nexus.File.u_data in
            match u.Nexus.File.u_model with
            | Some _ ->
                process_parsed_sequences 
                    false u.Nexus.File.u_weight default_tcm
                    Cost_matrix.Two_D.default Cost_matrix.Three_D.default
                    `DO false u.Nexus.File.u_alph file `Ml data seq
                    u.Nexus.File.u_model u.Nexus.File.u_pam
            | None ->
                let tcm,name = match u.Nexus.File.u_tcm with
                    | None ->
                        Cost_matrix.Two_D.of_transformations_and_gaps 
                                (size < 7) size 1 2 all_elements,
                        default_tcm
                    | Some (name,matrix) ->
                        let size = Array.length matrix in
                        assert( size > 0 );
                        let lst = Array.map Array.to_list matrix in
                        let lst = Array.to_list lst in
                        let use_comb = size < 7 in
                        Cost_matrix.Two_D.of_list ~use_comb lst 
                                    (if use_comb then (1 lsl size) - 1 else size),
                        (Input_file (name,lst))
                in
                let tcm,name = match u.Nexus.File.u_opening with
                    | None -> tcm, name
                    | Some v ->
                        if v=0 then tcm,name
                        else begin
                            Cost_matrix.Two_D.set_affine tcm (Cost_matrix.Affine v);
                            let name = match name with
                                | Substitution_Indel (a,b) -> 
                                    Substitution_Indel_GapOpening (a,b,v)
                                | Input_file (name,lst) ->
                                    Input_file_GapOpening (name,lst,v)
                                | _ -> assert( false )
                            in
                            tcm, name
                        end
                in
                let name = match u.Nexus.File.u_level with
                    | Some level -> Level (name, level)
                    | None       -> name
                in
                let tcm3d = Cost_matrix.Three_D.of_two_dim tcm in
                process_parsed_sequences false u.Nexus.File.u_weight name tcm
                                        tcm3d `DO false u.Nexus.File.u_alph file
                                        `Seq data seq None u.Nexus.File.u_pam
        in
        List.fold_left ~f:(single_sequence_adder) 
                       ~init:data
                        file_out.Nexus.File.unaligned
    in
    let cnsets,csets = (* create `character <--> setname` tables *)
        let character_nsets = create_ht () in
        let character_sets = create_ht () in
        Hashtbl.iter 
            (fun set_name set_specs ->
                List.iter
                    (fun set_spec ->
                        let cnames = 
                            Nexus.File.get_character_names 
                                file_out.Nexus.File.characters
                                file_out.Nexus.File.csets set_spec
                        in
                        if Hashtbl.mem character_sets set_name then
                            Hashtbl.replace character_sets set_name
                                (cnames@(Hashtbl.find character_sets set_name))
                        else
                            Hashtbl.add character_sets set_name cnames;
                        List.iter (fun name -> Hashtbl.add character_nsets name set_name)
                                  cnames)
                    set_specs)
            file_out.Nexus.File.csets;
        (character_nsets,character_sets)
    in
    (* create set for branches *)
    let tbl,found =
        branches_to_map data None (Some file_out.Nexus.File.branches)
                        file_out.Nexus.File.trees
    in
    let tbl = if found then Some tbl else None in
    let d = 
        {data with 
            character_sets = csets;
            character_nsets = cnsets;
            branches = tbl;
            iterate_branches = (not found); }
    in
    Status.finished st;
    new_codes, d

let add_static_parsed_file (data:d) (file:string) (triple:Nexus.File.nexus) : d =
    let _,d = gen_add_static_parsed_file true data file triple in d

let print_parsed_data lst = 
    let print_single (iopt,(str,nexus)) =
        Printf.printf "%s -- %s"
            (match iopt with | Some x -> string_of_int x | None -> "none") str;
        Printf.printf "Taxa: ";
        Array.iter (fun x -> Printf.printf "%s, "
                            (match x with | Some x -> x | None -> "none") )
                  nexus.Nexus.File.taxa;
        print_newline ();
        Printf.printf "Data: %d x %d" (Array.length nexus.Nexus.File.matrix)
                                      (Array.length nexus.Nexus.File.matrix.(0))
    in
    List.iter print_single lst


let add_multiple_static_parsed_file data list =
(*    print_parsed_data list;*)
    let data = duplicate data in
    let data,map = 
        List.fold_left 
            ~f:(fun (data,map) (code,(file, triple)) ->
                    let ncodes,data = gen_add_static_parsed_file false data file triple in
                    match code with
                    | Some code -> 
                        (data,All_sets.IntegerMap.add code (Array.to_list ncodes) map)
                    | None -> (data,map))
            ~init:(data,data.dynamic_static_codes)
            list
    in
    { data with dynamic_static_codes = map;
                static_dynamic_codes = reverse_dynamic_static_codes map; }


let add_static_file ?(report = true) style data (file : FileStream.f) = 
    try
        let ch, file = FileStream.channel_n_filename file in
        let r = 
            match style with
            | `Hennig -> Hennig.File.of_channel ch file 
            | `Nexus  -> Nexus.File.of_channel ch file
        in
        if report then
            report_static_input file r;
            close_in ch;
            add_static_parsed_file data file r
    with | Sys_error err ->
        let file = FileStream.filename file in
        let msg =
            Printf.sprintf
               ("Couldn't@ open@ file@ %s@ to@ load@ the@ data.@ The@ "^^
                "system@ error@ message@ is@ %s.")
               file err
        in
        output_error msg;
        data


let dna_lexer = Alphabet.Lexer.make_lexer false false Alphabet.nucleotides

let check_if_taxa_are_ok file taxa = 
    let the_great_majority_is_acgt (lst, _) = 
        let base = ref 0
        and others = ref 0 in
        List.iter (function 1 | 2 | 4 | 8 -> incr base | _ -> incr others) lst;
        !base > (3 * !others)
    in
    let has_spaces x = 
        FileStream.has_condition FileStream.is_taxon_delimiter (trim x)
    in
    let has_unacceptable x =
        FileStream.has_condition FileStream.is_unacceptable_in_taxon_name x
    in
    (* We want to check if the names are unique and if there is a suspicious
    * name of a taxon *)
    let _, second = List.fold_left ~f:(fun (acc, is_ok) x ->
        let msg =
            ("There@ is@ a@ taxon@ name@ that@ has@ "
            ^ "illegal@ characters@ on@ it@ ([]();,: ).@ This@ leaves@ the@ "
            ^ "generated@ trees@ "
            ^ "unreadable!.@ If you@ want@ to@ continue,@ that's@ "
            ^ "your@ call...@ the@ file@ is@ " ^ StatusCommon.escape file 
            ^ "@ and the@ taxon@ is@ " ^ StatusCommon.escape x)
        in
        (* We'll see, if we can parse it with the dna parser, we will tell the
        * user that the name of the taxon looks suspicious. *)
        let is_ok =
            let name = Stream.of_string x in
            try 
                let l = dna_lexer name [] 0 in
                if the_great_majority_is_acgt l then
                    Status.user_message Status.Warning 
                    ("There@ is@ a@ taxon@ name@ that@ is@ "
                    ^ "suspiciously@ simmilar@ to@ a@ DNA@ sequence@ in@ the@ "
                    ^ "file@ " ^ StatusCommon.escape file ^ 
                    "!.@ The@ taxon@ is@ " ^ x)
                else if has_spaces x then
                    Status.user_message Status.Warning msg
                else ();
                if has_unacceptable x then begin
                    Status.user_message Status.Error
                    ("Illegal@ terminal@ name:@ the@ characters@ @@ and %% may " ^
                    "cause@ POY@ to@ crash.@ Please@ remove@ them@ from@ your@ "
                    ^ "terminal@ names.");
                    false
                end else is_ok
            with
            | _ -> 
                    if has_spaces x then
                        Status.user_message Status.Error msg
                    else ();
                    if has_unacceptable x then begin
                        Status.user_message Status.Error
                        ("Illegal@ terminal@ name:@ the@ characters@ @@ and %% may " ^
                        "cause@ POY@ to@ crash.@ Please@ remove@ them@ from@ your@ "
                        ^ "terminal@ names.");
                        false
                    end else is_ok
        in
        if All_sets.Strings.mem x acc then begin
            Status.user_message Status.Warning 
            ("There@ is@ a@ taxon@ duplicated@ in@ the@ "
            ^ "file@ " ^ StatusCommon.escape file ^ "!.@ The@ duplicated@ taxon@ is@ " 
            ^ StatusCommon.escape x);
            acc, is_ok
        end else (All_sets.Strings.add x acc, is_ok)) ~init:(All_sets.Strings.empty,
        true) taxa
    in
    second

let print_error_message fl =
    let msg = "Unexpected@ character@ in@ file@ " ^ fl.Fasta.filename ^ 
    "@ in@ taxon@ " ^ StatusCommon.escape fl.Fasta.taxon ^ ".@ The@ character@ '" ^ 
    StatusCommon.escape fl.Fasta.character ^ "' " ^ 
    "is@ illegal@ in@ this@ file@ format." in
    Status.user_message Status.Error msg


let aux_process_molecular_file ?(respect_case = false) tcmfile tcm tcm3 alphabet processor builder dyna_state data file =
    begin try
        let ch = Parser.Files.molecular_to_fasta file in
        let res = 
            try Fasta.of_channel ~respect_case:respect_case (builder alphabet) ch with
            | Fasta.Illegal_molecular_format fl ->
                    let file = FileStream.filename file in
                    let fl = { fl with Fasta.filename = file } in
                    print_error_message fl;
                    raise (Fasta.Illegal_molecular_format fl)
        in
        let res = List.filter (function [[]], _ | [], _ -> false | _ -> true) res
        in
        let () = (* Output a message with the contents of the file *)
            let num_taxa = List.length res in
            let taxa_contents = 
                let file = FileStream.filename file in
                "@[The@ file@ " ^ StatusCommon.escape file ^ 
                "@ contains@ sequences@ of@ " ^
                string_of_int num_taxa ^ "@ taxa" 
            and sequence_contents = 
                if 0 = num_taxa then ""
                else begin
                    let lst, _ = List.hd res in
                    let add acc x = acc + (List.length (List.flatten x)) in
                    let len = List.fold_left ~f:add ~init:0 lst in
                    ",@ each@ sequence@ holding@ " ^
                    string_of_int len ^ "@ " ^ (if len > 1 then
                        "fragments.@]" else "fragment.@]@.")
                end
            in
            Status.user_message Status.Information (taxa_contents ^ sequence_contents);
        in
        close_in ch;
        if check_if_taxa_are_ok (FileStream.filename file) 
            (let _, names = List.split res in names) then
                processor alphabet res
        else begin
            Status.user_message Status.Error 
            ("Ignoring@ the@ file@ " ^ StatusCommon.escape (FileStream.filename
            file) ^ "@ as@ it@ contains@ illegal@ characters@ in@ the@ taxon@ "
            ^ "names@ (Andres@ is@ sure@ that@ POY@ will@ crash@ when@ you@ " ^
            "attempt@ to@ get@ the@ results).");
            data
        end
    with
    | Sys_error err ->
            let file = FileStream.filename file in
            let msg = "Couldn't@ open@ file@ " ^ file ^ "@ to@ load@ the@ " ^
                "dna@ sequences@ file.@ @ The@ system@ error@ message@ is@ "
                ^ err ^
                "." in
            output_error msg;
            data
    end

let process_molecular_file ?(respect_case = false) tcmfile tcm tcm3 annotated alphabet
                            mode is_prealigned dyna_state data file =
    let data =
        aux_process_molecular_file ~respect_case:respect_case
            tcmfile tcm tcm3 alphabet
            (fun alph parsed -> 
                process_parsed_sequences is_prealigned 1.0 
                tcmfile tcm tcm3 mode annotated 
                alph (FileStream.filename file) dyna_state data parsed None None)
            (fun x -> 
                if not is_prealigned then FileContents.AlphSeq x
                else FileContents.Prealigned_Alphabet x)
            dyna_state data file
    in
    data

let process_ignore_taxon data taxon =
    let res = All_sets.Strings.add taxon data.ignore_taxa_set in
    let taxon_names, taxon_codes =
        let code = All_sets.StringMap.find taxon data.taxon_names in
        let () =
            try Hashtbl.remove data.taxon_characters code with
            | Not_found -> ()
        in
        All_sets.StringMap.remove taxon data.taxon_names,
        All_sets.IntegerMap.remove code data.taxon_codes 
    in
    { data with ignore_taxa_set = res; taxon_codes = taxon_codes; taxon_names =
                    taxon_names}


let process_ignore_file data file = 
    try
        let ch, file = FileStream.channel_n_filename file in
        let taxa = Parser.IgnoreList.of_channel ch in
        let data = List.fold_left ~f:process_ignore_taxon ~init:data taxa in
        data
    with
    | Sys_error err ->
            let file = FileStream.filename file in
            let msg = 
                "Couldn't open file " ^ StatusCommon.escape file ^ 
                " to load the ignore taxa. The system error message is " 
                ^ err 
            in
            Status.user_message Status.Error msg;
            data

let complement_taxa data taxa = 
    let remover acc taxon = All_sets.StringMap.remove taxon acc in
    let the_taxa = List.fold_left ~f:remover ~init:data.taxon_names taxa in
    let elements = 
        All_sets.StringMap.fold (fun name _ lst -> name :: lst) the_taxa []
    in
    elements

let code_taxon code data = 
    All_sets.IntegerMap.find code data.taxon_codes

let report included excluded =
    let len1 = List.length included
    and len2 = List.length excluded in
    let total = max len1 len2 in
    let arr = Array.make_matrix (total + 1) 2 "" in
    let add row position item =
        arr.(position).(row) <- StatusCommon.escape item;
        position + 1
    in
    let total_included = List.fold_left ~f:(add 0) ~init:1 included
    and total_excluded = List.fold_left ~f:(add 1) ~init:1 excluded in
    arr.(0).(0) <- "Included";
    arr.(0).(1) <- "Excluded";
    Status.user_message Status.Information 
    "@[<v>@[Selected@ Terminals:@]@,";
    Status.output_table Status.Information arr;
    Status.user_message Status.Information 
    "@]@,%!";
    Status.user_message Status.Information 
    ("@[Total@ included:@ " ^ string_of_int (total_included - 1) ^ "@]@,");
    Status.user_message Status.Information 
    ("@[Total@ excluded:@ " ^ string_of_int (total_excluded - 1) ^ "@]@,")

let get_all_taxon_active_codes data = 
    Hashtbl.fold (fun code _ acc -> code :: acc) data.taxon_characters [] 

let absolute_of_percentage n precentage =
    truncate ((precentage *. (float_of_int n)) /. 100.)

let rec process_analyze_only_taxa meth data = 
    match meth with
    | `Random fraction ->
            let all_codes = Array.of_list (get_all_taxon_active_codes data) in
            Array_ops.randomize all_codes;
            if (fraction > 100. || fraction < 0.) then
                failwith "Illegal fraction";
            let n = absolute_of_percentage (Array.length all_codes) fraction in
            let taxa = 
                Array.to_list (Array.map (fun x -> code_taxon x data) 
                (Array.sub all_codes 0 n)) 
            in
            let res = List.fold_left ~f:process_ignore_taxon ~init:data taxa in
            res
    | `Names (dont_complement, taxa) ->
            let taxa = 
                warn_if_repeated_and_choose_uniquely taxa 
                "selected@ names" ""
            in
            let taxa = 
                if not dont_complement then taxa
                else complement_taxa data taxa 
            in
            let res = List.fold_left ~f:process_ignore_taxon ~init:data taxa in
            res
    | `CharSet (dont_complement, name_lst) ->
            let names = List.flatten
                (List.map (fun x -> try Hashtbl.find data.character_sets x
                                    with | Not_found -> [] ) name_lst)
            in
            process_analyze_only_taxa (`Names (dont_complement,names)) data

    | `Missing (dont_complement, fraction) ->
            Status.user_message Status.Information "Here";
            let fraction = (float_of_int fraction) /. 100. in
            let characters = 
                float_of_int 
                (Hashtbl.fold (fun _ _ acc -> acc + 1)
                data.character_codes 0)
            in
            let process_taxon txn_lst =
                Hashtbl.fold
                    (fun _ (_, item) spec -> match item with
                        | `Specified -> spec + 1
                        | `Unknown -> spec)
                    txn_lst 0
            in
            let add_or_remove_taxon code taxon_cs (included, excluded) =
                let spec = process_taxon taxon_cs in
                if fraction <= (float_of_int spec) /. characters then
                    included, (code_taxon code data) :: excluded
                else (code_taxon code data) :: included, excluded
            in
            let included, excluded =
                let included, excluded = 
                    Hashtbl.fold add_or_remove_taxon 
                    data.taxon_characters ([], [])
                in
                if dont_complement then included, excluded
                else excluded, included
            in
            report included excluded;
            process_analyze_only_taxa (`Names (false, excluded)) data

let process_analyze_only_file dont_complement data files =
    let data = duplicate data in
    try
        let appender acc file = 
            try 
                let ch, file = FileStream.channel_n_filename file in
                let taxa = Parser.IgnoreList.of_channel ch in
                let taxa = 
                    warn_if_repeated_and_choose_uniquely taxa "terminals file" file
                in
                close_in ch;
                taxa @ acc
                with
                | Sys_error err ->
                        let file = FileStream.filename file in
                        let msg = 
                            "Couldn't open file " ^ file ^ " to load the " ^
                            "terminals file. The system error message is " ^ 
                            err 
                        in
                        failwith msg
        in
        let taxa = List.fold_left ~f:appender ~init:[] files in
        let taxa = List.map trim taxa in
        let ignored, taxa = 
            if dont_complement then complement_taxa data taxa, taxa
            else taxa, complement_taxa data taxa
        in
        report taxa ignored;
        let res = List.fold_left ~f:process_ignore_taxon ~init:data ignored in
        res
    with
    | Failure msg ->
            Status.user_message Status.Error msg;
            data

let process_ignore_taxon data taxon = 
    (process_ignore_taxon data taxon)


let remove_taxa_to_ignore data = 
    let data = duplicate data in
    let process_data taxon data =
        try
            let tcode = 
                All_sets.StringMap.find taxon data.taxon_names 
            in
            Hashtbl.remove data.taxon_characters tcode;
            data
        with
        | _ -> data
    in
    let data = All_sets.Strings.fold process_data data.ignore_taxa_set data in
    let data = find_code_for_root_if_removed data in
    data

let report_terminals_files filename taxon_files ignored_taxa =
    let files = All_sets.StringMap.empty 
    and fo = Status.user_message (Status.Output (filename, false, [])) in
    let make_per_file_set taxon files acc =
        let is_ignored = All_sets.Strings.mem taxon ignored_taxa in
        let process_all_files file acc =
            if All_sets.StringMap.mem file acc then
                let (included, excluded) = All_sets.StringMap.find file acc in
                if is_ignored then
                    All_sets.StringMap.add file (included, (taxon :: excluded)) 
                    acc
                else 
                    All_sets.StringMap.add file ((taxon :: included), excluded)
                    acc
            else 
                if is_ignored then
                    All_sets.StringMap.add file ([], [taxon]) acc
                else All_sets.StringMap.add file ([taxon], []) acc
        in
        All_sets.Strings.fold process_all_files files acc
    and print_file file_name (included, excluded) =
        fo ("@,@[<v 2>@{<u>Input File:@} " ^ StatusCommon.escape file_name 
        ^ "@,@[<v 0>");
        let output_list str lst = 
            let total = string_of_int (List.length lst) in
            fo ("@,@[Terminals " ^ StatusCommon.escape str ^ 
            " (" ^ total ^ ")@]@[<v 2>@,@[<v 0>");
            List.iter (fun x -> fo (StatusCommon.escape x); fo "@,") lst;
            fo "@]@]";
        in
        output_list "Included" included;
        output_list "Excluded" excluded;
        fo "@]@]@,"
    in
    let files = All_sets.StringMap.fold make_per_file_set taxon_files files in
    fo "@[<v 2>@{<b>Terminals Files@}:@,@[<v 0>";
    All_sets.StringMap.iter print_file files;
    fo "@]@]%!"

let rec tcm_definition_to_string = function
    | Substitution_Indel_GapOpening (a, b, _)
    | Substitution_Indel (a, b) -> 
            "tcm:(" ^ string_of_int a ^ "," ^ string_of_int b ^ ")"
    | Input_file_GapOpening (name, _, _) 
    | Input_file (name, _) -> "'" ^ name ^ "'"
    | Level (def,_) -> tcm_definition_to_string def

let rec gap_opening_to_string = function
    | Input_file _ 
    | Substitution_Indel _ ->  "0"
    | Substitution_Indel_GapOpening (_, _, v)
    | Input_file_GapOpening (_, _, v)  -> string_of_int v
    | Level (def,_) -> tcm_definition_to_string def

(* This function will output in channel [ch] what was recorded in the
* [data]. This will have to be changed to XML. *)
let to_channel ch data = 
    let print_synonyms a b = 
        output_string ch a;
        output_string ch " -> ";
        output_string ch b;
        output_string ch "\n"
    in
    let print_taxa a _ =
        output_string ch a;
        output_string ch "\n"
    in
    let print_ignored_taxa a = 
        output_string ch a;
        output_string ch "\n"
    in
    let print_characters _ = function
        | Static a ->
                let str = Nexus.File.to_string a in
                output_string ch str;
                output_string ch "\n"
        | Dynamic dspec ->
                output_string ch dspec.filename;
                output_string ch ", ";
                output_string ch dspec.fs;
                output_string ch ", ";
                output_string ch (tcm_definition_to_string dspec.tcm);
                output_string ch "\n"
        | _ -> ()
    in
    output_string ch "List of synonyms: \n";
    All_sets.StringMap.iter print_synonyms data.synonyms;
    output_string ch "List of taxa: \n";
    All_sets.StringMap.iter print_taxa data.taxon_names;
    output_string ch "List of ignored taxa: \n";
    All_sets.Strings.iter print_ignored_taxa data.ignore_taxa_set;
    output_string ch "List of loaded characters: \n";
    Hashtbl.iter print_characters data.character_specs 


let add_character_spec spec code data =
    let data = duplicate data in
    (* We have to test whether or not it's a member already *)
    if not (Hashtbl.mem data.character_specs code) then begin
        Hashtbl.replace data.character_specs code spec 
    end else raise Illegal_argument;
    data

let rec taxon_code name data =
    try All_sets.StringMap.find name data.taxon_names with
    | Not_found -> taxon_code (All_sets.StringMap.find name data.synonyms) data

let get_tcm code data = 
    match Hashtbl.find data.character_specs code with
    | Static spec -> 
        begin match spec.Nexus.File.st_type with
            | Nexus.File.STSankoff x -> x
            | _ -> failwith "Unexpected"
        end
    | _ -> failwith "Unexpected"



(** funtion to see if a cost matrix matches the cost matrix associated with 
 *  the list first_lst.
 *  Parameters are
 *  code_cm is a cost matrix (int array array) 
 *  first_lst is an int list
 *  data is a Parser.Data.d type and this contains the encoding specs which
 *  contain the cost matrices.
*)
let match_cost_matrix code_cm first_lst data =
    match first_lst with
    | hd :: _ -> 
            let hd_cm = get_tcm hd data in
            if code_cm = hd_cm then true
            else false
    | [] -> false

(** [classify code data g] - is used for Sankoff characters to produce an
 *  int list list. Every code that has the same cost matrix is placed
 *  in the same int list.  
 *  This function is placed here rather than in parser.ml because is needs
 *  to use Data.get_tcm function and if this is used in parser.ml
 *  then there is a circular dependency between poyParser.ml and parser.ml
 *  Parameters are
 *  code is an int  
 *  data is a Parser.Data.d type and this contains the encoding specs which
 *  contain the cost matrices.
 *  g is the int list list being built
 *)
let classify code data = 
    let code_cm = get_tcm code data in
    let rec builder lst = 
        match lst with
        | [] -> [[code]]
        | hd :: tl -> 
                if match_cost_matrix code_cm hd data then
                    (code :: hd) :: tl
                else hd :: builder tl
    in
    { data with sankoff = builder data.sankoff }

let categorize data =
    (* We must return 7 lists of integers containing the codes for each class of
     * character *)

    (* We recategorize the data, so we must clear any already-loaded
       data *)
    let data = { data with
                     non_additive_1 = [];
                     non_additive_8 = [];
                     non_additive_16 = [];
                     non_additive_32 = [];
                     non_additive_33 = [];
                     additive = [];
                     sankoff = [];
                     dynamics = [];
                     kolmogorov = [];
                     static_ml = [];
               } in                         
    let data = repack_codes data in
    let categorizer code spec data =
        match spec with
        | Static enc -> (* Process static characters *)
                let observed = List.length enc.Nexus.File.st_observed in
                let between x y = observed >= x && observed <= y in
                (match enc.Nexus.File.st_type with
                | Nexus.File.STSankoff _ -> classify code data
                | Nexus.File.STOrdered when observed > 1 ->
                      { data with additive = code :: data.additive }
                | Nexus.File.STOrdered -> data
                | Nexus.File.STUnordered ->
                        if between 0 1 then
                            { data with non_additive_1 = code ::
                                data.non_additive_1 }
                        else if between 2 8 then
                            { data with non_additive_8 = code ::
                                data.non_additive_8 }
                        else if between 9 16 then
                            { data with non_additive_16 = code :: 
                                data.non_additive_16 }
                        else if between 17 32 then
                            { data with non_additive_32 = code :: 
                                data.non_additive_32 }
                        else if observed > 32 then
                            { data with non_additive_33 = code :: 
                                data.non_additive_33 }
                        else data
                | Nexus.File.STLikelihood _ ->
                        { data with static_ml = code :: data.static_ml })
        | Dynamic _ -> 
                { data with dynamics = code :: data.dynamics }
        | Set ->  data
        | Kolmogorov _ -> { data with kolmogorov = code :: data.kolmogorov }
    in
    let res = Hashtbl.fold categorizer data.character_specs data in
    res


let character_code name data = 
    Hashtbl.find data.character_names name

let code_character code data =
    Hashtbl.find data.character_codes code

let get_sequence_tcm seqcode data = 
    try match Hashtbl.find data.character_specs seqcode with
        | Dynamic dspec    -> dspec.tcm2d
        | Kolmogorov dspec -> dspec.dhs.tcm2d
        | _ -> failwith "Data.get_sequence_tcm: Not a dynamic character"
    with | Not_found as err ->
        let name = code_character seqcode data in
        let msg = "Could not find the code " ^ string_of_int seqcode ^ 
                  " with name " ^ StatusCommon.escape name in
        Status.user_message Status.Error msg;
        raise err

let get_sequence_alphabet seqcode data = 
    let chars = data.character_specs in
    try
        match Hashtbl.find chars seqcode with
        | Dynamic dspec ->
                dspec.alph
        | Kolmogorov dspec -> dspec.dhs.alph
        | _ -> failwith "Data.get_sequence_alphabet: Not a dynamic character"
    with
    | err ->
            let name = code_character seqcode data in
            let msg = "Could not find the code " ^ string_of_int seqcode ^ 
            " with name " ^ StatusCommon.escape name in
            Status.user_message Status.Error msg;
            raise err


let get_files data = List.rev data.files

let add_file data contents file = 
    let file = (FileStream.filename file), contents in
    { data with files = file :: data.files }

let get_taxa data = 
    All_sets.IntegerMap.fold (fun _ name acc -> name :: acc) data.taxon_codes []

let get_used_observed code d =
    let a = 
        match Hashtbl.find d.character_specs code with
        | Static a -> a
        | _ -> failwith "Data.get_used_observed expects a static character"
    in
    match a.Nexus.File.st_used_observed with
    | Some x -> x
    | None -> failwith "Data.get_used_observed"

let make_value_formatter x = (PXML -[Xml.Data.value][`String x]--)

let synonyms_to_formatter d : Xml.xml =
    let module T = Xml.Data in
    let synonym t1 t2 acc : Xml.xml Sexpr.t list =
        (PXML -[T.synonym]
            { make_value_formatter t1 }
            { make_value_formatter t2 } --) :: acc
    in
    let res = All_sets.StringMap.fold synonym d.synonyms [] in
    (RXML -[T.synonyms] { set res }--)

let taxon_names_to_formatter d : Xml.xml =
    let module T = Xml.Data in
    let taxon_name code name acc =
        if All_sets.Strings.mem name d.ignore_taxa_set then acc
        else
            (PXML -[T.taxon] 
                ([T.name] = [`String name])
                ([T.code] = [`Int code]) --) :: acc
    in
    let res = All_sets.IntegerMap.fold taxon_name d.taxon_codes [] in
    (RXML -[T.taxa] { set res } --)

let files_to_formatter d : Xml.xml =
    let module T = Xml.Data in
    let file (f, contents) : Xml.xml Sexpr.t =
        let aux item =
            let tmp = 
                `String 
                    (match item with
                    | Characters -> Xml.Data.characters
                    | CostMatrix -> Xml.Data.cost_matrix
                    | Trees -> Xml.Data.trees)
            in
            (PXML -[T.file_contents] [tmp] --)
        in
        (PXML -[T.file] 
            ([T.filename] = [`String f]) { set List.map aux contents } --)
    in
    (RXML -[T.files] { set List.map file d.files } --)

let ignored_taxa_to_formatter d : Xml.xml = 
    let module T = Xml.Data in
    let taxon x acc = (make_value_formatter x) :: acc in
    let res = All_sets.Strings.fold taxon d.ignore_taxa_set [] in
    (RXML -[T.ignored_taxa] { set res } --)

let ignored_characters_to_formatter d : Xml.xml = 
    let module T = Xml.Data in
    let res = List.map make_value_formatter d.ignore_character_set in
    (RXML -[T.ignored_characters] { set res } --)

let states_set_to_formatter enc : Xml.xml = 
    let module T = Xml.Characters in
    let set = Parser.OldHennig.Encoding.get_set enc in
    let add item acc = (PXML -[T.state] ([T.value] = [`Int item]) --) :: acc in
    (RXML -[T.states] { set All_sets.Integers.fold add set [] } --)

let pam_spec_to_formatter (state : dyna_state_t) pam =
    let module T = Xml.Characters in
    let option_to_string contents = function
        | Some x -> contents x
        | None -> assert false
    in
    let handle_bool = option_to_string (fun x -> `Bool x) 
    and handle_int = option_to_string (fun x -> `Int x)
    and handle_re_meth x = 
        let conversion =
            (function `Locus_Breakpoint x | `Locus_Inversion x -> `Int x)
        in
        match x with
        | Some x -> conversion x
        | None -> assert false
    in
    match (state : dyna_state_t) with
    | `Seq |`SeqPrealigned -> [T.clas, `String T.sequence]
    | others -> 
            let clas = 
                match others with
                | `Chromosome -> `String T.chromosome
                | `Genome -> `String T.genome
                | `Annotated -> `String T.annotated
                | `Breakinv -> `String T.breakinv
                | _ -> assert false
            in

            let deref ptr = 
                match ptr with
                | Some content -> content 
                | None -> failwith "It is a null pointer" 
            in      

            let locus_indel_o, locus_indel_e = deref pam.locus_indel_cost in
            let locus_indel_e = float locus_indel_e /. 100.0 in
            let locus_indel_str = `IntFloatTuple (locus_indel_o, locus_indel_e)
            in

            let chrom_indel_o, chrom_indel_e = deref pam.chrom_indel_cost in
            let chrom_indel_e = float chrom_indel_e /. 100.0 in
            let chrom_indel_str = `IntFloatTuple (chrom_indel_o, chrom_indel_e) 
            in

            (AXML
            ([T.clas] = [clas]) 
            (*([T.seed_len] = [handle_int pam.seed_len]) *)
            ([T.re_meth] = [handle_re_meth pam.re_meth])
            ([T.circular] = [handle_int pam.circular])
            ([T.locus_indel_cost] = [locus_indel_str])
            ([T.chrom_indel_cost] = [chrom_indel_str])
            ([T.chrom_hom] = [handle_int pam.chrom_hom])
            ([T.translocation] = [handle_int pam.translocation])
            (*([T.sig_block_len] = [handle_int pam.sig_block_len])
            ([T.rearranged_len] = [handle_int pam.rearranged_len])*)
            ([T.keep_median] = [handle_int pam.keep_median])
            ([T.swap_med] = [handle_int pam.swap_med])
            ([T.approx] = [handle_bool pam.approx])
            ([T.symmetric] =  [handle_bool pam.symmetric]))

let character_spec_to_formatter enc : Xml.xml =
    let module T = Xml.Characters in
    match enc with
    | Kolmogorov d ->
            (RXML -[T.kolmogorov]
                ([T.name] = [`String d.dhs.filename])
                (* TODO *)
                ([T.chars] = [`String ""])
                ([T.words] = [int d.ks.wordset])
                ([T.ints] = [int d.ks.intset]) --)
    | Static enc -> Nexus.File.to_formatter enc    
    | Dynamic dspec ->
            let initial =
                match dspec.initial_assignment with
                | `AutoPartitioned (_, size, _) -> 
                        `String 
                        ("AutoPartition of size " ^ string_of_int size ^ 
                        " with DO")
                | `Partitioned _ -> 
                        `String "User provided partition with DO"
                | `DO -> `String "Direct Optimization"
                | `GeneralNonAdd -> `String "Prealigned sequence."
                | `FS _ -> `String "Fixed States"
            in
            (RXML -[T.molecular]
                ([T.name] = [`String dspec.filename])
                ([T.initial_assignment] = [initial])
                ([T.tcm] = [`String (tcm_definition_to_string dspec.tcm)])
                ([T.gap_opening] = [`String (gap_opening_to_string dspec.tcm)])
                ([Xml.Characters.weight] = [`Float dspec.weight])
                ([pam_spec_to_formatter dspec.state dspec.pam])
                { single Alphabet.to_formatter dspec.alph } --)
    | Set -> failwith "TODO Set in Data.character_spec_to_formatter"

let characters_to_formatter d : Xml.xml =
    let module T = Xml.Data in
    let create code name acc =
        let enc = 
            try Hashtbl.find d.character_specs code with
            | Not_found as err -> 
                    prerr_int code;
                    prerr_newline ();
                    raise err 
        in
        (PXML -[T.characters] 
            ([T.name] = [`String name])
            ("code" = [`Int code])
            { single character_spec_to_formatter enc } --) :: acc
    in
    (RXML -[T.characters] { set Hashtbl.fold create d.character_codes [] } --)

let to_formatter attr d : Xml.xml =
    (RXML -[Xml.Data.data] ([attr]) 
        { single taxon_names_to_formatter d }
        { single synonyms_to_formatter d }
        { single ignored_taxa_to_formatter d }
        { single characters_to_formatter d }
        { single ignored_characters_to_formatter d }
        { single files_to_formatter d } 
        --)

(** transform dyna_pam_ls which is a list of dynamic parameters
* taken from poyCommand into dyna_pam which is structured as a record*)
let set_dyna_pam dyna_pam_ls old_dynpam =
    List.fold_left 
    ~f:(fun dyna_pam pam ->
        match pam with
        | `Align_Meth v -> {dyna_pam with align_meth = Some v}
        | `Median_Solver c -> {dyna_pam with median_solver = Some c}
        | `Annotate_Tool c -> {dyna_pam with annotate_tool = Some c}
        | `Locus_Inversion c -> {dyna_pam with re_meth = Some (`Locus_Inversion c)}
        | `Locus_Breakpoint c -> {dyna_pam with re_meth = Some (`Locus_Breakpoint c)}
        | `Translocation c -> {dyna_pam with translocation = Some c}
        | `Circular c -> 
              if c then {dyna_pam with circular = Some 1}
              else {dyna_pam with circular = Some 0}
        | `Locus_Indel_Cost c -> {dyna_pam with locus_indel_cost = Some c}
        | `Chrom_Indel_Cost c -> {dyna_pam with chrom_indel_cost = Some c}
        | `Chrom_Hom c -> {dyna_pam with chrom_hom = Some c}
        | `Keep_Median c -> 
                {dyna_pam with keep_median = Some c}
        | `SwapMed c -> {dyna_pam with swap_med = Some c}    
        | `Approx c  -> {dyna_pam with approx = Some c} 
        | `Symmetric c  -> {dyna_pam with symmetric = Some c}
        | `Max_3D_Len c -> {dyna_pam with max_3d_len = Some c} 
        | `Max_kept_wag c -> {dyna_pam with max_kept_wag = Some c}) 
    ~init:old_dynpam dyna_pam_ls

let use_mauve_annotator user_pam =
    match user_pam.annotate_tool with
    | Some (`Mauve (_,_,_,_))  -> true
    | Some (`Default (_,_,_)) -> false
    | None -> false

let get_dynas data dyna_code = 
    Hashtbl.fold
        (fun taxa_code ch_ls dyna_ls -> 
             let acs = 
                 try
                     match Hashtbl.find ch_ls dyna_code with
                     | ((Dyna _), _) as res -> Some res
                     | _ -> None
                 with  
                 | Not_found -> None
             in   
             match acs with 
             | Some acs ->
                   let dyna =   
                       match acs with   
                       | Dyna (code, dyna), _ -> dyna
                       | _, _ -> failwith "get dynamics"                          
                   in   
                   dyna :: dyna_ls  
             | None -> dyna_ls
        ) data.taxon_characters [] 

(** transform all sequences whose codes are on the code_ls into kolmogorov *)    
let transform_seqs_to_kolmogorov d codes newspec =
    let d = duplicate d in
    let set = 
        List.fold_left 
        ~f:(fun acc x -> All_sets.Integers.add x acc) 
        ~init:All_sets.Integers.empty
        codes 
    in
    let transformer code = function
        | Dynamic _ when (All_sets.Integers.mem code set) -> 
                Hashtbl.replace d.character_specs code (Kolmogorov newspec)
        | x -> ()
    in
    Hashtbl.iter transformer d.character_specs;
    d

(** [create_alpha_c2_breakinvs data chcode] transforms all 
* annotated chromosomes whose codes are on the code list [chcode]
* into breakinvs characters *)    
let create_alpha_c2_breakinvs (data : d) chcode =  

    let spec = Hashtbl.find data.character_specs chcode in  
    let c2, alpha,dynpam = match spec with 
        | Dynamic dspec -> dspec.tcm2d, dspec.alph, dspec.pam
        | _ -> failwith "Transfrom_annchroms_to_breakinvs: Not Dynamic" 
    in  
    let use_ukk =
        match dynpam.align_meth with
        | Some `NewKK   -> true
        | Some `Default -> false
        | None          -> false
    in
    let chrom_ls = get_dynas data chcode in 
        
    let max_code = List.fold_left  
        ~f:(fun max_code chrom -> 
                Array.fold_left ~f:(fun max_code seq -> max max_code seq.code 
                                   ) ~init:max_code chrom.seq_arr 
           ) ~init:0 chrom_ls 
    in  

    let gen_alpha = ref [] in 
    for code = 1 to max_code + 1 do 
        let char =  
            match code mod 2  with 
            | 1 -> string_of_int ( (code + 1)  / 2)  
            | 0 -> Alphabet.elt_complement ^ ( string_of_int (code / 2) ) 
            | _ -> failwith "compiler is error" 
        in  
        gen_alpha := (char, code, None)::!gen_alpha;         
    done;  

    let gen_gap_code = max_code + 2 in  
    gen_alpha := (Alphabet.gap_repr, gen_gap_code, None)::!gen_alpha;  
    gen_alpha := (Alphabet.elt_complement ^ Alphabet.gap_repr, (gen_gap_code + 1), None)::!gen_alpha;

    let gen_com_code = gen_gap_code + 2 in  
    gen_alpha := ("*", gen_com_code, None)::!gen_alpha;  
    gen_alpha := ("~*", (gen_com_code + 1), None)::!gen_alpha;  

    let max_code = gen_com_code + 1 in  
    let gen_alpha = 
        Alphabet.list_to_a (List.rev !gen_alpha) Alphabet.gap_repr (Some "*")
        Alphabet.Sequential
    in 
    (** Finish creating alphabet*)      


    let all_seq_arr = List.fold_left  
        ~f:(fun all_seq_arr chrom -> Array.append all_seq_arr chrom.seq_arr 
           ) ~init:[||] chrom_ls  
    in  

    let gen_cost_mat = Array.make_matrix max_code max_code max_int in  

    let num_seq = Array.length all_seq_arr in  
    for idx1 = 0 to num_seq - 2 do 
        for idx2 = idx1 + 1 to num_seq - 1 do 

            let seq1 = all_seq_arr.(idx1).seq in 
            let code1 = all_seq_arr.(idx1).code in 
            
            let seq2 = all_seq_arr.(idx2).seq in 
            let code2 = all_seq_arr.(idx2).code in 
                
            let _, _, cost =
                if use_ukk then
                Sequence.NewkkAlign.align_2 ~first_gap:false 
                seq1 seq2 c2 Sequence.NewkkAlign.default_ukkm
                else
                Sequence.Align.align_2 ~first_gap:false seq1 seq2 c2 Matrix.default
            in 
            gen_cost_mat.(code1).(code2) <- cost;
            gen_cost_mat.(code1).(code2 + 1) <- cost; 
            gen_cost_mat.(code1 + 1).(code2) <- cost; 
            gen_cost_mat.(code1 + 1).(code2 + 1) <- cost; 
            gen_cost_mat.(code2).(code1) <- cost; 
            gen_cost_mat.(code2).(code1 + 1) <- cost; 
            gen_cost_mat.(code2 + 1).(code1) <- cost; 
            gen_cost_mat.(code2 + 1).(code1 + 1) <- cost;                                 
        done; 
    done; 
    let gap = Alphabet.get_gap alpha in
    Array.iter  
        (fun chrom ->  
             let seq = chrom.seq in   
             let code = chrom.code in   
             let gap_seq = Sequence.create 1 in  
             let gap_seq = Sequence.prepend_char gap_seq gap in
             let _, _, cost = 
                 if use_ukk then 
                 Sequence.NewkkAlign.align_2 ~first_gap:false seq
                 gap_seq c2 Sequence.NewkkAlign.default_ukkm
                 else
                 Sequence.Align.align_2 ~first_gap:false seq
                 gap_seq c2 Matrix.default 
             in  
             gen_cost_mat.(code).(gen_gap_code) <- cost; 
             gen_cost_mat.(code + 1).(gen_gap_code) <- cost; 
             gen_cost_mat.(gen_gap_code).(code) <- cost; 
             gen_cost_mat.(gen_gap_code).(code + 1) <- cost;                 
        ) all_seq_arr;  
    gen_cost_mat.(gen_gap_code).(gen_gap_code) <- 0; 
    let gen_cost_ls = List.tl (Array.to_list gen_cost_mat) in 
    let gen_cost_ls = List.map  
        (fun cost_arr -> List.tl (Array.to_list cost_arr) ) gen_cost_ls 
    in  

    let gen_cost_mat = 
        Cost_matrix.Two_D.of_list ~use_comb:false gen_cost_ls 
        (if alpha = Alphabet.nucleotides then 31 
              else if (Alphabet.is_aminoacids alpha) 
              then 21 else (-1))
    in  

    gen_alpha, gen_cost_mat 

module Kolmogorov = struct

    let error_msg str =
        Status.user_message Status.Error 
        ("You@ have@ not@ defined@ yet@ the@ " ^ str)

    let extend_encodings enc = 
        let get_minimum x =
            let min_v = ref max_float in
            Array.iteri (fun pos v ->
                let mask = 1 lsl pos in
                if 0 <> x land mask then
                    min_v := min !min_v v
                else ()) enc;
            !min_v
        in
        let arr = Array.init 32 get_minimum in
        Bigarray.Array1.of_array Bigarray.float64 Bigarray.c_layout arr

    let calculate data indels substitution max_word_len max_indel_length
    event_prob =
        (* We begin by compiling the machine that has been loaded so far *)
        let find_function = 
            let (_, codes, repr) = Kolmo.Compiler.tree_of_decoder data.machine in
            fun str ->
                    let (_, code, _) = 
                        List.find ~f:(fun (name, _, _) -> name = str) codes 
                    in
                    let res =
                        float_of_int (List.length (All_sets.IntegerMap.find code
                        repr))
                    in
                    res
        in
        let rec log2 x acc = 
            if x > 0 then 
                log2 (x lsr 1) (acc + 1)
            else float_of_int acc
        in
        let bits_of_prob prob =
            ~-. ((log prob) /. (log 2.))
        in
        (* We can only handle for now homogeneous positions *)
        let root_cost = find_function "Tree.root"
        and branch_cost = find_function "Branch.interior"
        and leaf_cost = find_function "Branch.leaf"
        and end_cost = (find_function "Branch.ended") +. (find_function
        "Phylogeny.start") in
        (* We first look for the pairwise sequence alignment parameters *)
        let m_pairwise_algn = 
            (* The model must at least have the insertion and deletion functions
            * *)
            match indels, substitution with
            | Some (insertion_deletion), None ->
                    let ins, del =
                        match insertion_deletion with
                        | `Encoding (insertion, deletion) ->
                                find_function insertion, find_function deletion
                        | `Probs (insertion, deletion) ->
                                bits_of_prob insertion, bits_of_prob deletion
                    in
                    InDels (del, ins)
            | None, Some substitution ->
                    let sub = 
                        match substitution with
                        | `Encoding substitution ->
                                find_function substitution
                        | `Probs substitution ->
                                bits_of_prob substitution
                    in
                    Subs sub
            | Some (insertion_deletion), Some substitution ->
                    let sub, ins, del =
                        match insertion_deletion, substitution with
                        | `Encoding (insertion, deletion), `Encoding substitution ->
                                let sub = find_function substitution
                                and ins = find_function insertion
                                and del = find_function deletion in
                                sub, ins, del
                        | `Probs (insertion, deletion), `Probs substitution ->
                                bits_of_prob insertion,
                                bits_of_prob deletion,
                                bits_of_prob substitution
                        | _ -> assert false
                    in
                    if max_indel_length = 1 then
                        InDelSub (del, ins, sub)
                    else if max_indel_length > 1 then
                        let del = 
                            { selfp = del; distr = MaxLength max_indel_length } in
                        let ins = { del with selfp = ins } in
                        AffInDelSub (del, ins, sub)
                    else Subs sub
            | None, None -> assert false
        in
        let a = 2.
        and c = 2.
        and g = 2.
        and t = 2. in
        (* Now we need to find the chromosomal alignment parameters *)
        let calculate_encodings ins del =
            let enc v = v 
            and delhead v = v +. del in
            let simpleenc = [|enc a ; enc c; enc g; enc t; 0.0|] in
            let enc = extend_encodings simpleenc 
            and simplebend = 
                [| delhead a; delhead c; delhead g; delhead t; 0.0|] 
            in
            enc, simpleenc, simplebend
        in
        let kolmo_pairwise_align_specs =
            match m_pairwise_algn with
            | InDels (del, ins) ->
                    (* This is the expected list of functions in the case of 
                    * MDL0 *)
                    (* All the substitutions have a ridiculously high cost, 
                    * aligning has cost 0 *)
                    (* and pos = CharacSpec.length model "Pos" *)
                    let enc, simpleenc, simplebend = 
                        calculate_encodings ins del 
                    in
                    let a_ins = ins +. a
                    and c_ins = ins +. c
                    and g_ins = ins +. g
                    and t_ins = ins +. t in
                    let a_del = del 
                    and c_del = del 
                    and g_del = del 
                    and t_del = del in 
                    let sub = 1000. in
                    let tcm_matrix = 
                        [
                            [0.0; sub; sub; sub; a_del;];
                            [sub; 0.0; sub; sub; c_del;];
                            [sub; sub; 0.0; sub; g_del;];
                            [sub; sub; sub; 0.0; t_del;];
                            [a_ins; c_ins; g_ins; t_ins; 0.0;];
                        ]
                    in
                    { tm = tcm_matrix; be = enc; simplebe = simpleenc; 
                    simplebed = simplebend; ins = ins; del = del; sub = sub; 
                    ins_opening = 0.; del_opening = 0.; sub_opening = 0.; 
                    mo = m_pairwise_algn; root_cost = root_cost; branch_cost =
                        branch_cost; leaf_cost = leaf_cost; end_cost = end_cost; 
                        event_prob = event_prob; }
            | InDelSub (del, ins, sub) ->
                    (* This is the expected list of functions in the case of 
                    * MDL1 *)
                    (* and pos = CharacSpec.length model "Pos" *)
                    let enc, simpleenc, simplebend = 
                        calculate_encodings ins del 
                    in
                    let a_ins = ins +. a
                    and c_ins = ins +. c
                    and g_ins = ins +. g
                    and t_ins = ins +. t in
                    let a_del = del 
                    and c_del = del 
                    and g_del = del 
                    and t_del = del in
                    let a_sub = sub +. a
                    and c_sub = sub +. c
                    and g_sub = sub +. g
                    and t_sub = sub +. t in
                    let tcm_matrix = 
                        [
                            [0.0; a_sub; a_sub; a_sub; a_del];
                            [c_sub; 0.0; c_sub; c_sub; c_del];
                            [g_sub; g_sub; 0.0; g_sub; g_del];
                            [t_sub; t_sub; t_sub; 0.0; t_del];
                            [a_ins; c_ins; g_ins; t_ins; 0.0];
                        ]
                    in
                    if debug_kolmo then begin
                        let print_float v str = 
                            Status.user_message Status.Information (str ^ ": " ^
                            string_of_float v)
                        in
                        print_float a_ins "A Insertion";
                        print_float c_ins "C Insertion";
                        print_float g_ins "G Insertion";
                        print_float t_ins "T Insertion";
                        print_float a_del "A Deletion";
                        print_float c_del "C Deletion";
                        print_float g_del "G Deletion";
                        print_float t_del "T Deletion";
                        print_float a_sub "A Substitution";
                        print_float c_sub "C Substitution";
                        print_float g_sub "G Substitution";
                        print_float t_sub "T Substitution";
                        print_float a "A Encoding";
                        print_float c "C Encoding";
                        print_float g "G Encoding";
                        print_float t "T Encoding";
                    end;
                    { tm = tcm_matrix; be = enc; simplebe = simpleenc; 
                    simplebed = simplebend; ins = ins; del = del; sub = sub; 
                    ins_opening = 0.; del_opening = 0.; sub_opening = 0.;
                    mo = m_pairwise_algn; root_cost = root_cost; branch_cost =
                        branch_cost; leaf_cost = leaf_cost; end_cost = end_cost;
                    event_prob = event_prob;}
            | Subs sub ->
                    (* This is the expected list of functions in the case of 
                    * MDL2 *)
                    let indels = 1000. in
                    let enc, simpleenc, simplebend = 
                        calculate_encodings indels indels 
                    in
                    let a_sub = sub +. a
                    and c_sub = sub +. c
                    and g_sub = sub +. g
                    and t_sub = sub +. t in
                    let tcm_matrix = 
                        [
                            [0.0; a_sub; a_sub; a_sub; indels];
                            [c_sub; 0.0; c_sub; c_sub; indels];
                            [g_sub; g_sub; 0.0; g_sub; indels];
                            [t_sub; t_sub; t_sub; 0.0; indels];
                            [indels; indels; indels; indels; 0.0];
                        ]
                    in
                    { tm = tcm_matrix; be = enc; simplebe = simpleenc; 
                    simplebed = simplebend; ins = indels; del = indels; 
                    sub = sub; ins_opening = 0.; del_opening = 0.; 
                    sub_opening = 0.; mo = m_pairwise_algn; root_cost = root_cost; 
                    branch_cost = branch_cost; leaf_cost = leaf_cost; 
                    end_cost = end_cost; event_prob = event_prob; }
            | AffInDelSub (del, ins, sub) -> 
                    let calculate_encodings ins del =
                        let enc v = v 
                        and delhead v = v +. del.selfp in
                        let simpleenc = [|enc a ; enc c; enc g; enc t; 0.0|] in
                        let enc = extend_encodings simpleenc 
                        and simplebend = 
                            [| delhead a; delhead c; delhead g; delhead t; 0.0|]
                        in
                        enc, simpleenc, simplebend
                    in
                    let calculate_extension x = 2. in
                    let calculate_extra x = 
                        match x.distr with
                        | MaxLength y -> (log2 y 1) +. x.selfp +. 2. 
                    in
                    let enc, simpleenc, simplebend = 
                        calculate_encodings ins del 
                    in
                    let insop = calculate_extra ins
                    and delop = calculate_extra del 
                    and inex = calculate_extension ins
                    and delex = calculate_extension del in
                    let a_ins = a +. inex
                    and c_ins = c +. inex
                    and g_ins = g +. inex
                    and t_ins = t +. inex in
                    (* We will assign the cost of the base because we readjust the
                    * cost later in the Kolmogorov class of characters. However, by
                    * keeping it, we ensure that we will try to minimize the length
                    * of the resulting sequence. *)
                    let a_del = a +. delex
                    and c_del = c +. delex
                    and g_del = g +. delex
                    and t_del = t +. delex in
                    let a_sub = sub +. a
                    and c_sub = sub +. c
                    and g_sub = sub +. g
                    and t_sub = sub +. t in
                    if debug_kolmo then begin
                        let print_float v str = 
                            Status.user_message Status.Information (str ^ ": " ^
                            string_of_float v)
                        in
                        print_float a_ins "A Insertion";
                        print_float c_ins "C Insertion";
                        print_float g_ins "G Insertion";
                        print_float t_ins "T Insertion";
                        print_float a_del "A Deletion";
                        print_float c_del "C Deletion";
                        print_float g_del "G Deletion";
                        print_float t_del "T Deletion";
                        print_float a_sub "A Substitution";
                        print_float c_sub "C Substitution";
                        print_float g_sub "G Substitution";
                        print_float t_sub "T Substitution";
                        print_float a "A Encoding";
                        print_float c "C Encoding";
                        print_float g "G Encoding";
                        print_float t "T Encoding";
                        print_float insop "Insertion Opening";
                        print_float delop "Deletion Opening";
                    end;
                    let tcm_matrix = 
                        [
                            [0.0; a_sub; a_sub; a_sub; a_del];
                            [c_sub; 0.0; c_sub; c_sub; c_del];
                            [g_sub; g_sub; 0.0; g_sub; g_del];
                            [t_sub; t_sub; t_sub; 0.0; t_del];
                            [a_ins; c_ins; g_ins; t_ins; 0.0];
                        ]
                    in
                    { tm = tcm_matrix; be = enc; simplebe = simpleenc; 
                    simplebed = simplebend; ins = ins.selfp;
                    del = del.selfp; sub = sub; ins_opening = insop;
                    del_opening = delop; sub_opening = 0.;
                    mo = m_pairwise_algn; root_cost = root_cost; 
                    branch_cost = branch_cost; leaf_cost = leaf_cost; 
                    end_cost = end_cost; event_prob = event_prob;}
            | AffInDelAffSub _ -> failwith "Programming it"
        in
        kolmo_pairwise_align_specs, [], data

    let rec calculate_precision v p =
        if v < 1.0 then calculate_precision (v *. 10.0) (p *. 10.0)
        else p
end

let convert_dyna_taxon_data data (ch_ls : (int, cs) Hashtbl.t)
    tran_code_ls transform =
    let add chcode (ch_cs : cs) (acc_cs_ls : (int, cs) Hashtbl.t)
        : (int, cs) Hashtbl.t =
        match ch_cs with 
        | ((Dyna (chcode, dyna_data)) as ch, specified) ->  
              let new_ch = 
                  match (List.mem chcode tran_code_ls) with
                  | true -> 
                        let seq =  dyna_data.seq_arr.(0).seq in
                        let seq_code = dyna_data.seq_arr.(0).code  in
                        let specs = Hashtbl.find data.character_specs chcode in 
                        let specs = 
                            match specs with 
                            | Dynamic d -> d 
                            | Kolmogorov d -> d.dhs 
                            | _ -> failwith "Data.convert_dyna_taxon_data 1"
                        in 
                        let gap = Alphabet.get_gap specs.alph in 
                        let seq = 
                            match transform with 
                            | `Seq_to_Chrom _ | `Custom_to_Breakinv _ -> 
                                  (match (Sequence.get seq 0) = gap with
                                   | true -> Sequence.del_first_char seq  
                                   | false -> seq)  

                            | `Chrom_to_Seq _ | `Breakinv_to_Custom _ ->
                                    (match (Sequence.get seq 0) = gap with
                                    | true -> seq
                                    | false -> Sequence.prepend_char seq gap)

                            | `Annchrom_to_Breakinv _ ->
                                    let len = Array.length dyna_data.seq_arr in 
                                    Sequence.init 
                                    (fun idx -> dyna_data.seq_arr.(idx).code)
                                    len

                            | `Seq_to_Kolmogorov _
                            | `Change_Dyn_Pam _ -> seq 
                        in 
                        let seq_arr = 
                            match transform with 
                            | `Seq_to_Kolmogorov _
                            | `Change_Dyn_Pam _ -> dyna_data.seq_arr
                            | _ -> [|{seq = seq; delimiter = []; code = seq_code}|]
                        in 
                        let dyna_data = {seq_arr = seq_arr} in
                        Dyna (chcode, dyna_data)
                  | false -> ch
              in
              Hashtbl.replace acc_cs_ls chcode (new_ch, specified);
              acc_cs_ls
        | _ -> 
                Hashtbl.replace acc_cs_ls chcode ch_cs;
                acc_cs_ls
    in  
    Hashtbl.fold add ch_ls (create_ht ())

let get_state default = function
    | `Seq_to_Chrom _ -> `Chromosome
    | `Custom_to_Breakinv _ 
    | `Annchrom_to_Breakinv _ -> `Breakinv
    | `Chrom_to_Seq _ | `Breakinv_to_Custom _ -> `Seq
    | `Change_Dyn_Pam _ -> default
    | `Seq_to_Kolmogorov _ -> failwith "Illegal Data.get_state argument"


(** [get_sequences code data] outputs a stack containing all the sequences of the
* character [code] stored in [data]. The empty sequences (according to
* [Sequence.is_empy]) are not included in the stack. If the input code does not
* correspond to a sequence character, or the sequence contains more than one
* fragment, then it outputs an empty stack. *)
let get_sequences code data = 
    let alpha = get_sequence_alphabet code data in
    let gap = Alphabet.get_gap alpha in
    let seqs = Stack.create () in
    let process_taxon a b = 
        match Hashtbl.find b code with
        | (Stat _), _ -> ()
        | (Dyna (_, d)), _ ->
                match d.seq_arr with
                | [|dv|] ->
                        if not (Sequence.is_empty dv.seq gap) then 
                            Stack.push dv.seq seqs
                | _ -> ()
    in
    Hashtbl.iter process_taxon data.taxon_characters;
    seqs


(** [all_pairs_alignments seqs cm] outputs a stack containing all the pairwise
* comparisons between the sequences in the stack [seqs] employing the distance 
* function specified by [cm]. The output consists of the aligned versions of
* every sequence in [seqs].  *)
let all_pairs_alignments seqs cm =
    let seqs = Stack.copy seqs in
    let res = Stack.create () in
    while not (Stack.is_empty seqs) do
        let a = Stack.pop seqs in
        Stack.iter (fun b ->
            let alignment = Sequence.Align.align_2 a b cm Matrix.default in
            Stack.push alignment res) seqs;
    done;
    res


(** [alignments_of_code code data] is the same as [all_pairs_alignments]
 * excepting that the input requests the result for a particular character
 * [code] stored in [data]. The function generates an empty stack if
 * [get_sequences code data] would output an empty stack. *)
let alignments_of_code code data = 
    let seqs = get_sequences code data in
    let cm = get_sequence_tcm code data in
    let pairs = all_pairs_alignments seqs cm in
    pairs

type sequence_statistics = {
    max_length : int;
    min_length : int;
    sum_lengths : int;
    sequences : int;
    max_distance : float;
    min_distance : float;
    sum_distances : float; }

let statistics_of_alignments seqs pairs =
    let cnt = Stack.length seqs
    and d_min = ref max_int
    and d_max = ref 0
    and d_sum = ref 0 
    and s_max = ref 0
    and s_min = ref max_int 
    and s_sum = ref 0 in
    (* Gathed the distance data *)
    Stack.iter
        (fun (_, _, cost) ->
            d_min := min !d_min cost;
            d_max := max !d_max cost;
            d_sum := !d_sum + cost;)
        pairs;
    (* Gather the sequence data *)
    Stack.iter
        (fun seq ->
            let len = Sequence.length seq in
            s_max := max !s_max len;
            s_min := min !s_min len;
            s_sum := !s_sum + len)
        seqs;
    { max_length = !s_max;
        min_length = !s_min;
        sum_lengths = !s_sum;
        sequences = cnt;
        max_distance = float_of_int !d_max;
        min_distance = float_of_int !d_min;
        sum_distances = float_of_int !d_sum }

let sequence_code_statistics data code =
    let seqs = get_sequences code data 
    and pairs = alignments_of_code code data
    and name = code_character code data in
    name, statistics_of_alignments seqs pairs

let kolmo_round_factor = 100.

let event_frequency algnments = 
    let samples = ref 0
    and best_cost = ref (float_of_int max_int)
    and adder = ref 0. in
    let add x = 
        incr samples;
        adder := x +. !adder 
    in
    let update_item (a, b, cost) =
        let cost = float_of_int cost in
        if cost < !best_cost  then begin
            samples := 0;
            adder := 0.;
            best_cost := cost;
            let len = Sequence.length a in
            assert (len = Sequence.length b);
            let distance = ref 1 in
            for i = 0 to len - 1 do
                let a = Sequence.get a i
                and b = Sequence.get b i in
                if 0 <> (a land b) then incr distance
                else begin
                    add (1. /. (float_of_int !distance));
                    distance := 1
                end
            done;
        end else ()
    in
    Stack.iter update_item algnments;
    let probs = !adder /. (float_of_int !samples) in
    probs


let convert_dyna_spec data chcode spec transform_meth =  
    match spec with   
    | Kolmogorov _ 
    | Dynamic _ ->
            let dspec = 
                match spec with
                | Kolmogorov x -> x.dhs
                | Dynamic x -> x
                | _ -> failwith "Impossible"
            in
        let () =
            (* First check if the transformation is legal  *)
            match dspec.state, transform_meth with 
            | `Seq, `Seq_to_Chrom _
            | `Seq, `Custom_to_Breakinv _
            | `Annotated, `Annchrom_to_Breakinv _
            | `Chromosome, `Chrom_to_Seq  _
            | `Breakinv, `Breakinv_to_Custom _
            | `Seq, `Seq_to_Kolmogorov _
            | _, `Change_Dyn_Pam _ -> ()
            | _, _ -> failwith "Illegal character transformation requested"
        in
        begin match transform_meth with (*
        TODO ANDRES
        | `Seq_to_Kolmogorov model ->
                let a, b = 
                    match model with
                    | `AtomicIndel (_, Some (insert, delete, substitute)) 
                    | `AffineIndel (_, Some (insert, delete, substitute)) ->
                            (Some (`Probs (insert, delete))),
                                Some (`Probs substitute)
                    | `AtomicIndel (_, None) 
                    | `AffineIndel (_, None) ->
                            (Some (`Encoding ("Dna.insert", "Dna.delete"))), 
                            Some (`Encoding "Dna.substitute") 
                in
                let seqs = get_sequences chcode data 
                and algnments = alignments_of_code chcode data in
                let stats = statistics_of_alignments seqs algnments in
                let kolmo_event_probability = event_frequency algnments in
                let kolmo_max_length = 
                    stats.max_length + (stats.max_length / 2)
                in
                let c, d, e =
                    match model with
                    | `AtomicIndel (user_kolmo_event_probability, _)
                    | `AffineIndel (user_kolmo_event_probability, _) -> 
                            kolmo_max_length, 1, 
                                match user_kolmo_event_probability with
                                | Some x -> x
                                | None -> kolmo_event_probability
                in
                let data = 
                    { data with 
                        machine = Kolmo.SimpleIndels.apply_model
                        kolmo_max_length Kolmo.Compiler.compiler }
                in
                let kolmospec, dyn_spec_options, data = 
                    Kolmogorov.calculate data a b c d e in
                let tcm = 
                    (* Round the tcm *)
                    List.map ~f:(List.map 
                        ~f:(fun x -> truncate (x *.  kolmo_round_factor))) 
                    kolmospec.tm
                in
                let ks = { 
                    funset = a, b;
                    wordset = c;
                    intset = d;
                    kolmo_spec = kolmospec;
                }
                and dspec = 
                    let tcm = Cost_matrix.Two_D.of_list tcm 31 in 
                    let bed = 
                        Array.map (fun x -> truncate (x *.
                        kolmo_round_factor)) kolmospec.simplebed
                    and be = 
                        Array.map (fun x -> truncate (x *.
                        kolmo_round_factor)) kolmospec.simplebe
                    in
                    Cost_matrix.Two_D.fill_tail bed tcm;
                    Cost_matrix.Two_D.fill_prepend be tcm;
                    if (kolmospec.ins_opening <> 0.) then
                        Cost_matrix.Two_D.set_affine tcm (Cost_matrix.Affine
                        (truncate (kolmospec.ins_opening *.
                        kolmo_round_factor)))
                    else ();
                    let old_dynpam = dspec.pam in
                    let pam = set_dyna_pam dyn_spec_options old_dynpam in
                    { dspec with 
                    tcm = "Kolmogorov"; 
                    tcm2d = tcm;
                    pam = pam;} 
                in
                Kolmogorov { dhs = dspec; ks = ks }, data *)
        | transform_meth ->
            let (old_dynpam:dyna_pam_t) = dspec.pam in
            let (al, c2), pam = 
                (* Now we can transform *)
                match transform_meth with 
                | `Seq_to_Kolmogorov _ -> failwith "Impossible"
                | `Annchrom_to_Breakinv pam_ls -> 
                        create_alpha_c2_breakinvs data chcode,
                        set_dyna_pam pam_ls old_dynpam

                | `Seq_to_Chrom pam_ls
                | `Custom_to_Breakinv pam_ls
                | `Change_Dyn_Pam pam_ls
                | `Chrom_to_Seq pam_ls 
                | `Breakinv_to_Custom pam_ls ->
                        let pam = set_dyna_pam pam_ls old_dynpam in
                        (dspec.alph, dspec.tcm2d), pam
            in
            let state = get_state dspec.state transform_meth in
            let pam_state = match state with
                | `Chromosome -> Some `Chromosome
                | `Genome     -> Some `Genome
                | `Breakinv   -> Some `Breakinv
                | _           -> None
            in
            let new_spec = 
                Dynamic { dspec with 
                          alph = al; 
                          tcm2d = c2; 
                          pam = { pam with mode = pam_state } }
            in
            new_spec,data
        end
    | _ -> failwith "Convert_dyna_spec: Not a dynamic character" 

let check_fraction fraction =
    if (100. < fraction || fraction < 0.) then
        failwith "Illegal fraction"
    else ()

let select_random_sublist fraction lst =
    let n = absolute_of_percentage (List.length lst) fraction in
    let arr = Array.of_list lst in
    Array_ops.randomize arr;
    Array.to_list (Array.sub arr 0 n )

type classes = 
    [ `Dynamic |  `NonAdditive | `StaticLikelihood | `DynamicLikelihood | `Likelihood
    | `Additive | `Sankoff | `Kolmogorov | `AllStatic | `AllDynamic ] 

let has_likelihood d = 
    match d.static_ml with
    | [] -> false
    | _  -> true

let type_of_dynamic_likelihood d = match d.dynamics with
    | []      -> None
    | hd :: _ ->
        match Hashtbl.find d.character_specs hd with
        | Dynamic spec when spec.state = `Ml ->
            begin match spec.lk_model with
                | Some s -> Some s.MlModel.spec.MlModel.cost_fn 
                | None   -> assert false
            end
        | _ -> None

let rec get_code_from_name data name_ls = 
  let code_ls = List.fold_right 
      ~f:(fun name acc -> 
              try 
                  let code = Hashtbl.find_all data.character_names name in
                  code @ acc
              with 
              | Not_found -> 
                      (* We will try to use a regular expression to match the
                      * requested item *)
                    let nname = Str.regexp name in
                    match Hashtbl.fold (fun stored_name code acc -> 
                        if Str.string_match nname stored_name 0 then code :: acc
                        else acc) data.character_names []
                    with
                    | [] -> 
                            failwith 
                            ("The@ character@ " ^ name ^ 
                            ",@ target@ of@ the@ transformation@ does@ not@ " ^
                            "exist.");
                    | r -> r @ acc
         ) name_ls ~init:[] 
  in
  code_ls

and get_code_with_missing dont_complement data fraction = 
    let taxa = 
        let tmp = 
            All_sets.StringMap.fold (fun _ _ acc -> acc + 1)
            data.taxon_names 0
        in
        float_of_int tmp
    in
    let count_occurrences data =
        let extract_info data =
            Hashtbl.fold (fun code _ acc -> 
                All_sets.IntegerMap.add code 0 acc) 
            data.character_codes All_sets.IntegerMap.empty
        in
        let add_counter _ x counter =
            match x with
            | Dyna (y, _), `Specified 
            | Stat (y, _), `Specified -> 
                    let cnt = All_sets.IntegerMap.find y counter in
                    All_sets.IntegerMap.add y (cnt + 1) counter
            | _, `Unknown -> counter
        in
        let add_taxon_to_character_counters _ taxon_cs_lst counters = 
            Hashtbl.fold add_counter taxon_cs_lst counters
        in
        Hashtbl.fold add_taxon_to_character_counters
        data.taxon_characters (extract_info data)
    in
    let fraction = (float_of_int fraction) /. 100. in
    let codes = 
        All_sets.IntegerMap.fold (fun code count acc ->
        if fraction <= ((float_of_int count) /. taxa) then code :: acc
        else acc) (count_occurrences data) []
    in
    if dont_complement then codes
    else 
        match complement_characters data (`Some codes ) with
        | `Some x -> x
        | _ -> failwith "Data.get_code_with_missing"
(**Give a list of characters, return their codes*)    
and get_code_from_characters_restricted kind (data : d) (chs : characters) = 
    let kind_lst = match kind with
        | `Dynamic ->
            data.dynamics
        | `DynamicLikelihood ->
            List.filter
                (fun x -> match Hashtbl.find data.character_specs x with
                    | Dynamic x when x.state = `Ml -> true
                    | _ -> false)
                data.dynamics
        | `NonAdditive ->
                data.non_additive_1 @ data.non_additive_8 @
                data.non_additive_16 @ data.non_additive_32 @
                data.non_additive_33
        | `Additive -> data.additive
        | `Sankoff -> List.flatten data.sankoff
        | `Kolmogorov -> data.kolmogorov
        | `AllDynamic -> data.kolmogorov @ data.dynamics
        | `StaticLikelihood -> data.static_ml
        | `Likelihood -> 
            (get_code_from_characters_restricted `StaticLikelihood data chs) @
            (get_code_from_characters_restricted `DynamicLikelihood data chs)
        | `AllStatic -> 
                data.non_additive_1 @ data.non_additive_8 @
                data.non_additive_16 @ data.non_additive_32 @
                data.additive @ data.non_additive_33 @
                data.static_ml @ (List.flatten data.sankoff)
    in
    let rec items chs = match chs with
        | `Some code_ls -> code_ls 
        | `Range (file, x,y) -> items (transform_range_to_codes file x y)
        | `Names name_ls -> get_code_from_name data name_ls
        | `Random fraction ->
                check_fraction fraction;
                select_random_sublist fraction kind_lst
        | `All -> kind_lst
        | `AllDynamic | `AllStatic as m -> 
                get_code_from_characters_restricted m data `All
        | `Missing (dont_complement, fraction) ->
                get_code_with_missing dont_complement data fraction
        | `CharSet str_lst ->
            let names =
                List.flatten 
                    (List.map 
                        ( fun x -> try Hashtbl.find data.character_sets x
                                   with | Not_found -> [])
                        str_lst)
            in
            items (`Names names)
    in
    List.filter (fun x -> List.exists (fun y -> y = x) kind_lst) (items chs)

and get_all_codes data =
    Hashtbl.fold (fun c _ acc -> c :: acc) data.character_codes  []

and get_chars_codes data = function
    | `All -> get_all_codes data 
    | `Random fraction ->
            check_fraction fraction;
            select_random_sublist fraction (get_all_codes data)
    | `Some codes -> codes
    | `Range (file,x,y) -> get_chars_codes data (transform_range_to_codes file x y)
    | `Names names ->
            let names = 
                warn_if_repeated_and_choose_uniquely names 
                "selected@ characters@ " ""
            in
            let get_code acc name =
                try 
                    (Hashtbl.find data.character_names name) :: acc 
                with
                | Not_found as err ->
                        let nname = Str.regexp name in
                        match
                            Hashtbl.fold (fun item code acc ->
                                if Str.string_match nname item 0 then 
                                    code :: acc
                                else acc)
                            data.character_names acc
                        with
                        | [] -> 
                                Status.user_message Status.Error
                                ("Could@ not@ find@ any@ character@ matching@
                                the@ expression@ " ^ StatusCommon.escape name);
                                raise err
                        | r -> r 
            in
            List.fold_left ~f:get_code ~init:[] names
    | `AllStatic | `AllDynamic as m -> 
            get_code_from_characters_restricted m data `All
    | `Missing (dont_complement, fraction) ->
            get_code_with_missing dont_complement data fraction
    | `CharSet sets ->
            let names =
                List.flatten
                    (List.map
                        (fun x -> Hashtbl.find data.character_sets (String.uppercase x))
                        sets)
            in
            get_chars_codes data (`Names names)

and complement_characters data characters = 
    let codes = get_chars_codes data characters in
    let res = Hashtbl.fold (fun x _ acc -> 
        if List.exists (fun y -> x = y) codes then acc
        else x :: acc) data.character_codes [] in
    `Some res

let get_all_codes data =
    Hashtbl.fold (fun c _ acc -> c :: acc) data.character_codes  []

let get_code_from_characters_restricted_comp kind d ch =
    let dont_complement, chars = match ch with
        | `Some (dont_complement, x) -> dont_complement, `Some x
        | `Names (dont_complement, x) -> dont_complement, `Names x
        | `CharSet (dont_complement, x) -> dont_complement, `CharSet x
        | `Range (dont_complement, file, x, y) -> dont_complement, transform_range_to_codes file x y
        | `Random _ | `Missing _ | `All | `AllDynamic | `AllStatic as x -> true, x
    in
    let chars = get_code_from_characters_restricted kind d chars in
    if dont_complement then chars
    else 
        match complement_characters d (`Some chars) with
        | `Some x -> x
        | _ -> failwith "Impossible?"


let get_chars_codes_comp data ch =
    let dont_complement, ch = match ch with
        | `Some (x, y) -> x, `Some y
        | `Names (x, y) -> x, `Names y 
        | `CharSet (x, y) -> x, `CharSet y
        | `Range (x,file,a,b) -> x, transform_range_to_codes file a b
        | `Random _ | `Missing _ | `All | `AllStatic | `AllDynamic as x -> true, x
    in
    let codes = get_chars_codes data ch in
    if dont_complement then codes 
    else 
        match complement_characters data (`Some codes) with
        | `Some x -> x
        | _ -> failwith "Impossible?"


(* non/functional creation of sets *)
let make_set_partitions (functional:bool) (data:d) (name:string) (ccodes:Methods.characters) = 
    let name = String.uppercase name in
    let sets,nsets = 
        if functional 
            then Hashtbl.copy data.character_sets, Hashtbl.copy data.character_nsets
            else data.character_sets, data.character_nsets
    in
    let ncodes =
        List.map
            (fun x -> 
                try Hashtbl.find data.character_codes x 
                with | Not_found ->
                    failwithf "Cannot find %d in character_codes" x)
            (get_chars_codes_comp data ccodes)
    in
    let () = 
        if Hashtbl.mem sets name then
            let old_data = Hashtbl.find sets name in
            Hashtbl.replace sets name (ncodes @ old_data)
        else
            Hashtbl.add sets name ncodes
    in
    let () = List.iter (fun c -> Hashtbl.replace nsets c name) ncodes in
    if functional then
        { data with character_sets = sets; character_nsets = nsets; }
    else
        data


let categorize_static_likelihood_by_model chars data =
    let get_spec i = 
        let model = 
            match Hashtbl.find data.character_specs i with
                | Static spec ->
                    begin match spec.Nexus.File.st_type with
                        | Nexus.File.STLikelihood model -> model
                        | _ -> assert false
                    end
                | _ -> assert false
        in
        model.MlModel.spec
    in
    `Some (get_chars_codes_comp data chars)
        --> get_code_from_characters_restricted `StaticLikelihood data
        --> MlModel.categorize_by_model get_spec


let categorize_sets data : int list list =
    let rec inner_find i = function
        | is::_ when List.mem i is -> true
        |  _::iss -> inner_find i iss
        | []      -> false
    in
    let named_sets = 
        List.map
            ~f:(fun xs ->
                    List.map 
                        ~f:(fun x -> Hashtbl.find data.character_names x)
                        xs)
            (Hashtbl.fold (fun k v acc -> v::acc) data.character_sets [])
    in
    (* add un-named sets from characters *)
    let fcodes = 
        List.filter
            (fun x -> not (inner_find x named_sets))
            (get_chars_codes_comp data `All)
    in
    match fcodes with
    | [] -> named_sets
    | fc -> fc :: named_sets


let categorize_characters_comp data chars = match chars with
    | `All -> categorize_sets data
    | othr ->
        let found set1 set2 = match set1 with
            | x::xs when List.mem x set2 ->
                assert( List.fold_left ~f:(fun acc x -> acc & (List.mem x set2)) ~init:true xs);
                true
            | _::xs -> 
                assert( not (List.fold_left ~f:(fun acc x -> acc & (List.mem x set2)) ~init:true xs));
                false
            | [] -> 
                false
        in
        let selected = get_chars_codes_comp data chars in
        let all_sets = categorize_sets data in
        List.fold_left
            ~f:(fun acc xs ->
                if found xs selected then xs :: acc else acc)
            ~init:[]
            all_sets

let categorize_characters data chars = match chars with
    | `All -> categorize_sets data
    | othr ->
        let found set1 set2 = match set1 with
            | x::xs when List.mem x set2 ->
                assert( List.fold_left ~f:(fun acc x -> acc & (List.mem x set2)) ~init:true xs);
                true
            | _::xs -> 
                assert( not (List.fold_left ~f:(fun acc x -> acc & (List.mem x set2)) ~init:true xs));
                false
            | [] -> 
                false
        in
        let selected = get_chars_codes data chars in
        let all_sets = categorize_sets data in
        List.fold_left
            ~f:(fun acc xs ->
                if found xs selected then xs :: acc else acc)
            ~init:[]
            all_sets


let make_codon_partitions functional data name ccodes =
    let process_set name set nset codes = 
        let concat = if name = "" then "" else ":" in
        let name1 = String.uppercase (name^concat^"codon1")
        and name2 = String.uppercase (name^concat^"codon2")
        and name3 = String.uppercase (name^concat^"codon3") in
        let rec process_three acc1 acc2 acc3 = function
            | [] -> acc1, acc2, acc3
            | one::two::three::xss ->
                let () = Hashtbl.replace nset one name1
                and () = Hashtbl.replace nset two name2
                and () = Hashtbl.replace nset three name3 in
                process_three (one::acc1) (two::acc2) (three::acc3) xss
            (* non-divisible situations *)
            | _::_::[]
            |    _::[] ->
                failwithf ("Data in %s has %d characters; not divisible by "^^
                           "three and transformable to codons") name (List.length codes)
        in
        let acc1,acc2,acc3 = process_three [] [] [] codes in
        let () = Hashtbl.remove set name
        and () = Hashtbl.add set name1 acc1
        and () = Hashtbl.add set name2 acc2
        and () = Hashtbl.add set name3 acc3 in
        ()
    in
    let sets,nsets = 
        if functional 
            then Hashtbl.copy data.character_sets, Hashtbl.copy data.character_nsets
            else data.character_sets, data.character_nsets
    in
    let ncodes =
        List.map (fun x -> Hashtbl.find data.character_codes x) 
                 (get_chars_codes_comp data ccodes)
    in
    let () = process_set name sets nsets ncodes in
    { data with
        character_sets = sets;
        character_nsets = nsets; }

let get_tcm2d data c =
    match Hashtbl.find data.character_specs c with
    | Dynamic dspec -> dspec.tcm2d
    | Kolmogorov dspec -> dspec.dhs.tcm2d
    | _ -> failwith "Data.get_tcm2d"

let get_tcm3d data c =
    match Hashtbl.find data.character_specs c  with
    | Dynamic dspec -> dspec.tcm3d
    | Kolmogorov dspec -> dspec.dhs.tcm3d
    | _ -> failwith "Data.get_tcm3d"

let get_tcmfile data c = match Hashtbl.find data.character_specs c  with
    | Dynamic dspec    -> dspec.tcm
    | Kolmogorov dspec -> dspec.dhs.tcm
    | Static dspec     -> failwith "Data.get_tcmfile; Cannot transform static data"
    | Set              -> failwith "Data.get_tcmfile; Cannot transform set data"

let get_model_opt data c = 
    match Hashtbl.find data.character_specs c with
    | Dynamic dspec -> dspec.lk_model
    | Static sspec ->
        begin match sspec.Nexus.File.st_type with
            | Nexus.File.STLikelihood model -> Some model
            | _ -> None
        end
    | _ -> None


let get_alphabet data c =
    try match Hashtbl.find data.character_specs c  with
        | Dynamic dspec    -> dspec.alph
        | Kolmogorov dspec -> dspec.dhs.alph
        | Static  sspec    -> sspec.Nexus.File.st_alph
        | Set              -> failwithf "Data.get_alphabet: Finding %d alphabet in Set" c
    with  Not_found        -> failwithf "Data.get_alphabet: Couldn't find %d in character specs" c


let verify_alphabet data chars =
    match List.map (get_alphabet data) chars with
    | h :: t ->
        ignore 
            (List.fold_left
                ~f:(fun acc x -> 
                        if x = h then acc 
                        else begin
                            Alphabet.print x;
                            Alphabet.print h;
                            failwith "The alphabet of the characters is different"
                        end)
                ~init:true t);
        let alph = Alphabet.to_sequential h in
        (Alphabet.size alph), alph
    | [] -> failwith "No alphabet to verify?"

(* [independent c d] make each character in the characterset independent. Used
 * so each character can have a different model
let independent chars data = 
    let data_loop tbl num =
        (* chain to work through to change the code... *)
        (* specs -> Static (static_spec) -> st_type -> STLikelihood (ml_model) -> set_code *)
        let c_spec = Hashtbl.find tbl num in
        let n_ = match c_spec with
            | Static ss -> Static (Nexus.File.change_ml_code (next_likelihood_set ()) ss)
            | _ -> failwith "I only do this for likelihood characters" in
        Hashtbl.replace tbl num n_
    in
    List.iter (data_loop data.character_specs) 
              (get_code_from_characters_restricted_comp (`Likelihood) data chars);
    data
*)
(*
let failwithf format = Printf.ksprintf failwith format
 [make_char_sets sets d] makes a character set with a name based on a list of
 * ranges or single values
let make_char_sets sets d =
    let rec perform_on_each fn x = 
        let rec on_range l h = 
            if l = h then [fn h]
            else (fn l) :: (on_range (l+1) h) in
        match x with
        | (`Single x) :: tl -> (fn x) :: (perform_on_each fn tl)
        | (`Range (x,y)) :: tl -> assert ( x < y );
            (on_range x y) :: (perform_on_each fn tl)
        | [] -> []
    in
    let make_list_of_set_and_add setname setcode stuffs =
        perform_on_each
            (fun x -> 
                let c_name =  try Hashtbl.find d.character_codes x
                    with | Not_found -> failwithf "Cannot find character %d" x
                in
                let tbl = Hashtbl.find data.character_specs x in
                Hashtbl.replace tbl x setcode;
                c_name
            ) stuffs
    in
    List.iter 
        (fun x -> match x with
         | (name,stuffs)::tl ->
                let scode =
                    try Hashtbl.find d.character_codes name
                    with | Not_found ->
                        let c = next_likelihood_set () in
                        Hashtbl.add d.character_codes name c;
                        c
                in
                Hashtbl.add d.character_sets
                            name
                            (make_list_of_set_and_add name scode stuffs)
         | [] -> ()
        ) sets;
    d
*)

(** [compute_priors data chars] computes the observed frequencies for all the
 * elements in the alphabet shared by the characters *)
let compute_priors data chars u_gap = 
    let debug_priors = false in
    let size, alph = verify_alphabet data chars in
    let size = if u_gap then size else size-1 in
    let priors = Array.make size 0.0 in
    (* A function that takes a list of states and add the appropriate value to
    * each of the components of the priors *)
    let inverse = 1. /. (float_of_int size) in
    let longest = ref 0 and lengths = ref [] in
    let gap_char = Alphabet.get_gap alph in (* 4, usually *)
    let counter = ref 0 and gap_counter = ref 0 in
    (* A function to add the frequencies in all the taxa from the characters
    * specified in the list. *)
    let taxon_adder tax taxon_chars =
        let total = ref 0 in
        let adder char_code =
            let (cs, _) = Hashtbl.find taxon_chars char_code in
            match cs with
            | Dyna (_, dyna_data ) ->
                Array.iter
                    (fun x ->
                        total := (Sequence.length x.seq) - 1 + !total;
                        counter := (Sequence.length x.seq) - 1 + !counter;
                        for i = 1 (* skip initial gap *) to (Sequence.length x.seq) - 1 do
                            let lst = BitSet.Int.list_of_packed (Sequence.get x.seq i) in
                            if List.exists (fun x -> x = gap_char) lst && not u_gap || (lst = []) then
                                assert false
                            else
                                let inv = 1.0 /. (float_of_int (List.length lst)) in
                                List.iter (fun x -> priors.(x) <- priors.(x) +. inv) lst
                        done)
                    dyna_data.seq_arr
            | Stat (_,s) -> 
                Nexus.File.compute_static_priors alph u_gap (priors,counter,gap_counter) inverse s
        in
        List.iter adder chars;
        longest := max !total !longest;
        lengths := !total :: !lengths
    in
    Hashtbl.iter taxon_adder data.taxon_characters;
    let counter = float_of_int !counter
    and gcounter = float_of_int !gap_counter
    and gap_contribution = (float_of_int !gap_counter) /. (float_of_int size) in
    if debug_priors then begin
        Printf.printf "Computed Priors of %.1f char + %.1f gaps: " counter gcounter;
        Array.iter (Printf.printf "|%f") priors;
        Printf.printf "|]\n%!";
    end;
    let final_priors = 
        if u_gap then begin
            let total_added_gaps = 
                float_of_int
                    (List.fold_left ~f:(fun acc x -> (!longest - x) + acc) ~init:0 !lengths)
            in
            if debug_priors then Printf.printf "Total added gaps = %f\n%!" total_added_gaps;
            priors.(gap_char) <- priors.(gap_char) +. total_added_gaps;
            let counter = counter -. gcounter +. total_added_gaps;
            and weight  = gcounter /. (float_of_int size) in
            Array.map (fun x -> (x -. weight) /. counter) priors
        end else begin
            Array.map (fun x ->(x -. gap_contribution) /. counter) priors
        end
    in
    if debug_priors then begin
        let sum = Array.fold_left ~f:(fun a x -> a +. x) ~init:0.0 final_priors in
        Printf.printf "Final Priors (%f): [" sum;
        Array.iter (Printf.printf "|%f") final_priors;
        Printf.printf "|]\n%!";
    end;
    final_priors


(* this is the only function that applies a likelihood model to a set of
 * characters. it ensures that the dynamic characters have gap as an additional
 * state *)
let apply_likelihood_model_on_char_table replace data table codes model = 
    let model_opt = ref None
    and model_enc = Nexus.File.STLikelihood model
    and is_likelihood = function Nexus.File.STLikelihood _ -> true | _ -> false
    and priors_ () = compute_priors data codes true in
    if replace then begin
        (* replace only likelihood characters in set of codes *)
        List.iter
            (fun code -> match Hashtbl.find table code with
                | Static x when (is_likelihood x.Nexus.File.st_type) -> 
                    Hashtbl.replace table code
                        (Static {x with Nexus.File.st_type = model_enc; })
                | Dynamic ({state = state} as x) when state = `Ml ->
                    let () = match !model_opt,model.MlModel.spec.MlModel.use_gap with
                        | None,`Independent | None,`Coupled _ ->
                            model_opt := Some model
                        | None,`Missing -> 
                            model_opt := Some (MlModel.add_gap_to_model priors_ model)
                        | Some _,_ -> ()
                    in
                    Hashtbl.replace table code
                        (Dynamic {x with state = `Ml;
                                     lk_model = !model_opt;})
                | _ -> ())
            codes
    end else begin
        (* replace all characters in the set of codes *)
        List.iter
            (fun code -> match Hashtbl.find table code with
            | Static x ->
                Hashtbl.replace table code
                    (Static {x with Nexus.File.st_type = model_enc; })
            | Dynamic x ->
                let () = match !model_opt,model.MlModel.spec.MlModel.use_gap with
                    | None,`Independent | None,`Coupled _ ->
                        model_opt := Some model
                    | None, `Missing -> 
                        model_opt := Some (MlModel.add_gap_to_model priors_ model)
                    | Some _,_ -> ()
                in
                Hashtbl.replace table code
                    (Dynamic {x with state = `Ml;
                                     lk_model = !model_opt;})
            | _ -> failwith "Illegal conversion")
        codes
    end

let update_priors data charcodes use_gap = 
    let data = {data with character_specs = Hashtbl.copy data.character_specs; } in
    let current = get_likelihood_model data charcodes in
    match current.MlModel.spec.MlModel.base_priors with
    | MlModel.Estimated _  -> 
        let new_priors = compute_priors data charcodes use_gap
        and model_s = ref None and model_d = ref None in
        List.iter
            (fun code -> 
                match Hashtbl.find data.character_specs code,!model_s,!model_d with
                | Static x, Some model,_ ->
                    Hashtbl.replace data.character_specs code
                        (Static {x with Nexus.File.st_type = model})
                | Static x, None, _ -> 
                    let get_lk_model y = match y.Nexus.File.st_type with 
                        | Nexus.File.STLikelihood x -> x
                        | _ -> assert false
                    in
                    let model =
                        Nexus.File.STLikelihood
                            (MlModel.replace_priors (get_lk_model x) new_priors)
                    in
                    model_s := Some model;
                    Hashtbl.replace data.character_specs code
                        (Static {x with Nexus.File.st_type = model; })
                | Dynamic x,_,((Some model) as m) ->
                    Hashtbl.replace data.character_specs code
                        (Dynamic {x with state = `Ml; lk_model = m;})
                | Dynamic x,_,None ->
                    let get_lk_model = function
                        | Some m -> m
                        | None -> assert false
                    in
                    let model = MlModel.replace_priors (get_lk_model x.lk_model) new_priors in
                    model_d := Some model;
                    Hashtbl.replace data.character_specs code
                        (Dynamic {x with state = `Ml; lk_model = !model_d;})
                | (Kolmogorov _ | Set) , _ , _ -> assert false)
            charcodes;
        data
    | MlModel.Equal 
    | MlModel.Given _ -> data


let apply_likelihood_model_on_chars data char_codes (model:MlModel.model) = 
    let new_specs = Hashtbl.copy data.character_specs in
    apply_likelihood_model_on_char_table false data new_specs char_codes model;
    { data with character_specs = new_specs } --> categorize


(* remove the absent/present columns in the parsed data, and modify the previous
 * column (the column the active/present column is adding data to), and modify
 * it to a gap for likelihood. This is only for likelihood characters, unless
 * ignore is set to true, in which case it will remove them with impunity! The
 * active/present columns in other characters are used for an additional cost to
 * that column, combined with the weight parameter. *)
let remove_absent_present_encodings ?(ignore_data=false) data chars =
    let is_likelihood = function
        | Nexus.File.STLikelihood _ -> (true || ignore_data)
        | _ -> false
    (* transform a character to a gap *)
    and apply_gap_to_cs taxon_code char_code spec state = match spec,state with
        | Static sspec, (Stat (code,data),specified) ->
            let gap = Alphabet.get_gap sspec.Nexus.File.st_alph in
            (Stat (code, Some (`List [gap])),specified)
        | _, _ -> failwith "Data.remove_absent_present_encodings failure"
    in
    (* test if we are using STATIC likelihood somewhere *)
    let test_using_likelihood data =
        let ret = ref false in
        Hashtbl.iter
            (fun x spec -> match spec with
                | Static {Nexus.File.st_type = st_type} ->
                    ret := !ret || (is_likelihood st_type)
                | _ -> ())
            data.character_specs;
        (!ret || ignore_data)
    (* find if we are using the absent/present alphabet *)
    and is_present_absent a =
        try let () = ignore (Alphabet.match_base "present" a) 
            and () = ignore (Alphabet.match_base "absent" a) in
            true
        with | _ -> false
    (* apply the absent/present encoding column to the previous column *)
    and is_present tcode code state spec = match spec,state with
        | Static sspec, (Stat (code,Some data),_) ->
            let data_list = Nexus.File.static_state_to_list data
            and present = Alphabet.match_base "present" sspec.Nexus.File.st_alph in
            List.mem present data_list
        | Static sspec, (Stat (code,None),_) ->
            ignore (Alphabet.match_base "present" sspec.Nexus.File.st_alph);
            true
        | _,_ -> failwith "Incorrect Match1"
    and is_absent tcode code state spec = match spec,state with
        | Static sspec, (Stat (code,Some data),_) ->
            let data_list = Nexus.File.static_state_to_list data
            and absent = Alphabet.match_base "absent" sspec.Nexus.File.st_alph in
            List.mem absent data_list
        | Static sspec, (Stat (code,None),_) ->
            (* do the lookup anyway to ensure it exists *)
            ignore (Alphabet.match_base "absent" sspec.Nexus.File.st_alph);
            false
        | _,_ -> failwith "Incorrect Match2"
    and branch_remove branch name = match branch with
        | Some branch ->
            Hashtbl.iter
                (fun name v_set ->
                    All_sets.IntSetMap.iter 
                        (fun _ tbl -> Hashtbl.remove tbl name) (v_set) )
                branch
        | None -> ()
    in
    let apply_absent_encoding branch specs chars names codes code = 
        Hashtbl.iter
            (fun t_code t_state_tbl ->
                let apstate = Hashtbl.find t_state_tbl code
                and apspec  = Hashtbl.find specs code in
                if is_present t_code code apstate apspec then begin
                    let lkstate = Hashtbl.find t_state_tbl (code-1) 
                    and lkspec  = Hashtbl.find specs (code-1) in
                    let nstate = apply_gap_to_cs t_code (code-1) lkspec lkstate in
                    Hashtbl.replace t_state_tbl (code-1) nstate;
                    Hashtbl.remove t_state_tbl code;
                end else begin
                    assert(is_absent t_code code apstate apspec);
                    Hashtbl.remove t_state_tbl code
                end)
            chars;
        let name = Hashtbl.find codes code in
        Hashtbl.remove names name;
        Hashtbl.remove codes code;
        Hashtbl.remove specs code;
        branch_remove branch name
    in
    if test_using_likelihood data then begin
        let copy_spec   = Hashtbl.copy data.character_specs
        and copy_names  = Hashtbl.copy data.character_names
        and copy_codes  = Hashtbl.copy data.character_codes
        and copy_branch = copy_branches data.branches 
        and copy_char   =
            let n = Hashtbl.create (Hashtbl.length data.taxon_characters) in
            Hashtbl.iter (fun k v -> Hashtbl.add n k (Hashtbl.copy v))
                         data.taxon_characters;
            n
        in
        let map,chars =
            List.fold_left
                ~f:(fun (acc,chars) k -> match Hashtbl.find data.character_specs k with
                    | Static {Nexus.File.st_alph = st_alph;}
                        when is_present_absent st_alph ->
                            apply_absent_encoding copy_branch copy_spec copy_char copy_names copy_codes k;
                            IntMap.map
                                (fun x -> List.filter (fun y -> not (y = k)) x)
                                acc, k::chars
                    | Static x ->
                            acc, chars
                    | (Dynamic _ | Set | Kolmogorov _) -> assert false)
                ~init:(data.dynamic_static_codes,[])
                (get_code_from_characters_restricted_comp `AllStatic data chars)
        in
        { data with character_specs      = copy_spec;
                    taxon_characters     = copy_char;
                    character_names      = copy_names;
                    character_codes      = copy_codes;
                    dynamic_static_codes = map;
                    branches             = copy_branch;
                    static_dynamic_codes = reverse_dynamic_static_codes map; },
        chars
    end else begin
        data,[]
    end



(** [set_parsimony lk chars] transforms the characters specified in [chars] to
 * the likelihood model specified in [lk] *)
let set_parsimony data chars =
IFDEF USE_LIKELIHOOD THEN
    (* apply a type to a set of STATIC characters in data *)
    let transform_chars data chars =
        let new_specs = Hashtbl.copy data.character_specs in
        List.iter
            (fun code -> match Hashtbl.find new_specs code with
                | Static x ->
                    Hashtbl.replace new_specs code
                        (Static { x with Nexus.File.st_type = Nexus.File.STUnordered})
                | Dynamic x ->
                    Hashtbl.replace new_specs code
                        (Dynamic {x with lk_model = None; state = `Seq;})
                | Set ->
                    failwith "Cannot transform Set characters to parsimony")
            chars;
        new_specs
    in
    let chars = get_chars_codes_comp data chars in
    match chars with
        | [] -> data (* no characters selected; do not duplicate data *)
        | chars ->
            let new_specs = transform_chars data chars in
            categorize { data with character_specs = new_specs }
ELSE
    data
END


(** [set_likelihood lk chars] transforms the characters specified in [chars] to
* the likelihood model specified in [lk] *)
let set_likelihood data (((chars,_,_,_,_,use_gap) as m_spec):Methods.ml_spec) =
IFDEF USE_LIKELIHOOD THEN
    let u_gap = match use_gap with 
        | `Independent | `Coupled _ -> true | `Missing -> false
    in
    let transform_char_set data chars =
        match chars with
        | [] -> data
        | x::xs ->
            let dynamic = match Hashtbl.find data.character_specs x with
                | Dynamic _ -> true
                | _ -> false
            in
            (* We get the characters and filter them out to have only static types *)
            let model =
                let compute_priors () = compute_priors data chars u_gap in
                let alph_size,alph = verify_alphabet data chars in
                let lk_spec = MlModel.convert_methods_spec alph_size compute_priors m_spec in
                let lk_spec = 
                    if dynamic then MlModel.remove_gamma_from_spec lk_spec 
                               else lk_spec 
                in
                MlModel.create alph lk_spec
            in
            apply_likelihood_model_on_chars data chars model
    in
    let d,removed = remove_absent_present_encodings ~ignore_data:true data chars in
    let all_chars = categorize_characters_comp d chars in
    List.fold_left ~f:(transform_char_set) ~init:d all_chars

ELSE
    failwith "Likelihood not enabled: download different binary or contact mailing list"

END

let get_tran_code_meth data meth = 
    let tran_code_ls, meth =
        let a, b = 
            match meth with
            | `Seq_to_Chrom (a, b) ->
                    a, `Seq_to_Chrom b
            | `Custom_to_Breakinv (a, b) ->
                    a, `Custom_to_Breakinv b
            | `Annchrom_to_Breakinv (a, b) ->
                    a, `Annchrom_to_Breakinv b
            | `Change_Dyn_Pam (a, b) ->
                    a, `Change_Dyn_Pam b
            | `Chrom_to_Seq (a, b) ->
                    a, `Chrom_to_Seq b
            | `Breakinv_to_Custom (a, b) ->
                    a, `Breakinv_to_Custom b
            | `Seq_to_Kolmogorov (a, b) ->
                    a, `Seq_to_Kolmogorov b
        in
        let dont_complement, codes =
            match a with
            | `Some (dont_complement, codes) ->
                    dont_complement, `Some codes
            | `Range (dont_complement, file, x,y) ->
                    dont_complement, transform_range_to_codes file x y
            | `Names (dont_complement, names) ->
                    dont_complement, `Names names
            | `CharSet (dont_complement, names) ->
                    dont_complement, `CharSet names
            | `Random _ | `Missing _ | `All | `AllDynamic | `AllStatic as x -> true, x
        in
        let codes = get_code_from_characters_restricted `AllDynamic data codes in
        let codes = 
            if dont_complement then codes
            else 
                match complement_characters data (`Some codes) with
                | `Some x -> x
                | _ -> failwith "Impossible?"
        in
        codes, b
    in
    tran_code_ls, meth

(** transform all sequences whose codes are on the code_ls into chroms *)    
let transform_dynamic (meth: Methods.dynamic_char_transform) data =
    let tran_code_ls, meth = get_tran_code_meth data meth in 
    let data = ref (duplicate data) in
    Hashtbl.iter
    (fun code spec ->
         if List.mem code tran_code_ls then begin
             let dyn, d = convert_dyna_spec !data code spec meth in
             Hashtbl.replace !data.character_specs code dyn;
             data := d;
        end else ()) !data.character_specs;
    let new_taxon_chs = 
        let new_tbl = create_ht () in
        Hashtbl.iter 
            (fun code ch_ls -> 
                let new_ls = convert_dyna_taxon_data !data ch_ls tran_code_ls meth in
                Hashtbl.replace new_tbl code new_ls) 
            !data.taxon_characters;
        new_tbl
    in 
    {!data with taxon_characters = new_taxon_chs}

let intmap_filter f y =
    All_sets.IntegerMap.fold (fun a b acc ->
        if f a b then All_sets.IntegerMap.add a b acc
        else acc) y All_sets.IntegerMap.empty

let hashtbl_filter f y =
    Hashtbl.iter (fun a b ->
        if not (f a b) then Hashtbl.remove y a
        else ()) y;
    y

let process_ignore_character report data code_set =
    let data = duplicate data in
    let compare x = not (All_sets.Integers.mem x code_set) in
    let compare1 code _ = not (All_sets.Integers.mem code code_set) in
    try
        (* Remove each character in the code set *)
        let new_cign =
            All_sets.Integers.fold
                (fun code new_cign ->
                    let name = Hashtbl.find data.character_codes code in
                    Hashtbl.remove data.character_names name;
                    Hashtbl.remove data.character_codes code;
                    Hashtbl.remove data.character_specs code;
                    name :: new_cign)
                code_set 
                data.ignore_character_set
        in
        let non_additive_1 = List.filter compare data.non_additive_1 
        and non_additive_8 = List.filter compare data.non_additive_8
        and non_additive_16 = List.filter compare data.non_additive_16
        and non_additive_32 = List.filter compare data.non_additive_32
        and non_additive_33 = List.filter compare data.non_additive_33 
        and additive = List.filter compare data.additive 
        and static_ml = List.filter compare data.static_ml
        and sankoff = List.map (fun x -> List.filter compare x) data.sankoff 
        and dynamics = List.filter compare data.dynamics in
        Hashtbl.iter 
            (fun code lst ->
                Hashtbl.replace data.taxon_characters code (hashtbl_filter compare1 lst))
            data.taxon_characters;
        let sankoff = List.filter (function [] -> false | _ -> true) sankoff in
        { data with
            ignore_character_set = new_cign;
            non_additive_1 = non_additive_1;
            non_additive_8 = non_additive_8;
            non_additive_16 = non_additive_16;
            non_additive_32 = non_additive_32;
            non_additive_33 = non_additive_33;
            additive = additive;
            sankoff = sankoff;
            dynamics = dynamics;
            static_ml = static_ml;
        }
    with | Not_found -> 
        let msg = 
            "Could not find a character. " ^
            "I will ignore it and continue without processing it." 
        in
        Status.user_message Status.Error msg;
        data

let process_ignore_characters report data characters = 
    let codes = get_chars_codes data characters in
    let codes = 
        List.fold_left  ~f:(fun acc x -> All_sets.Integers.add x acc)
                        ~init:All_sets.Integers.empty codes
    in
    process_ignore_character report data codes

let process_analyze_only_characters_file report dont_complement data files =
    let codes = 
        List.fold_left 
            ~f:(fun acc x ->
                let ch, x = FileStream.channel_n_filename x in
                let items = Parser.IgnoreList.of_channel ch in
                close_in ch;
                let items = 
                    warn_if_repeated_and_choose_uniquely items "characters@ file@ " x
                in
                let codes = get_chars_codes data (`Names items) in
                acc @ codes)
            ~init:[]
            files
    in
    let items = 
        if dont_complement then  (`Some codes)
        else complement_characters data (`Some codes) 
    in
    process_ignore_characters report data items

let process_ignore_characters_file report data file =
    try
        let ch, file = FileStream.channel_n_filename file in
        let items = Parser.IgnoreList.of_channel ch in
        let items = List.map trim items in
        close_in ch;
        process_ignore_characters report data (`Names items)
    with | err ->
        let file = FileStream.filename file in
        let msg = 
            Printf.sprintf
                ("Error while attempting to read the ignore characters file %s."
                    ^^" I will continue without processing it.")
                (StatusCommon.escape file)
        in
        Status.user_message Status.Error msg;
        data

let replace_name_in_spec name = function
    | Static e -> Static { e with Nexus.File.st_name = name }
    | Dynamic dspec -> Dynamic { dspec with filename = name }
    | Kolmogorov d ->
            Kolmogorov { d with dhs = { d.dhs with filename = name } }
    | Set -> Set


let process_rename_characters data (a, b) = 
    let data = duplicate data in
    if Hashtbl.mem data.character_names b then
        raise Illegal_argument
    else begin
        let code = character_code a data in
        let spec = Hashtbl.find data.character_specs code in
        Hashtbl.replace data.character_names b code;
        Hashtbl.replace data.character_codes code b;
        Hashtbl.replace data.character_specs code (replace_name_in_spec b spec);
        data
    end

exception Invalid_Character of int

(** transform all sequences whose codes are on the code_ls into chroms 
 * each ia for one character*)    
let transform_chrom_to_rearranged_seq data meth tran_code_ls 
        (ia_ls : Methods.implied_alignment list) = 
    let data = duplicate data in
    let tran_code_ls, _ = get_tran_code_meth data meth in 
    let num_ia = ref 0 in  
    let t_ch_ia_map = List.fold_left 
    ~f:(fun (t_ch_ia_map) ia  ->    
        List.fold_left 
        ~f:(fun t_ch_ia_map ((handle_ia, _), _) ->
            List.fold_left 
            ~f:(fun t_ch_ia_map (t_id, t_ia) ->
                List.fold_left 
                ~f:(fun t_ch_ia_map ch_set_ia -> 
                    IntMap.fold  
                        (fun chcode t_ch_ia t_ch_ia_map ->                   
                            num_ia := Array.length t_ch_ia;
                            FullTupleMap.add (t_id,chcode) t_ch_ia 
                            t_ch_ia_map
                        ) ch_set_ia t_ch_ia_map
                    ) ~init:t_ch_ia_map t_ia
                ) ~init:t_ch_ia_map handle_ia
            ) ~init:t_ch_ia_map ia
        ) ~init:FullTupleMap.empty ia_ls  
    in  
    let data = List.fold_left 
        ~f:(fun data char_code ->
            let char_name = Hashtbl.find data.character_codes char_code in
            let tcm = get_tcm2d data char_code 
            and tcm3d = get_tcm3d data char_code 
            and alph = get_alphabet data char_code
            and tcmfile = get_tcmfile data char_code in
            let data = 
                process_ignore_characters true data (`Some [char_code]) 
            in
            let seqs = 
                IntMap.fold 
                    (fun (t_code : int) (t_name : string) seqs ->
                        try
                            let ia_arr = FullTupleMap.find (t_code, char_code) t_ch_ia_map in
                            let ia_arr = 
                                Array.map 
                                    (function `DO ia | `Last ia | `First ia ->
                                        let seq = Sequence.of_code_arr ia Alphabet.gap in 
                                        Sequence.prepend_char seq Alphabet.gap)
                                    ia_arr
                            in 
                            ([[Array.to_list ia_arr]], t_name) :: seqs
                        with Not_found -> seqs) 
                    data.taxon_codes []
            in                   
            process_parsed_sequences false 1.0 tcmfile tcm tcm3d `DO false 
                                     alph char_name `Seq data seqs None None)
        ~init:data tran_code_ls 
    in
    categorize data


let get_sequences_of_code data code =
    let process_taxon tcode chars acc =
        try let (tc, _) = Hashtbl.find chars code in
            match tc with
            | Dyna (_, c) -> (tcode, (c.seq_arr.(0)).seq)  :: acc
            | _ -> failwith "Impossible?"
        with | Not_found -> acc
    in
    (Hashtbl.fold process_taxon 
    data.taxon_characters [])

let auto_partition mode data code =
    let mode = match mode with
        | `Clip -> Clip
        | `NoClip -> NoClip
    in
    let dhs = match Hashtbl.find data.character_specs code with
        | Dynamic dhs -> dhs
        | _ -> assert false
    in
    let taxa = 
        let gap = Alphabet.get_gap dhs.alph in
        let taxa = get_sequences_of_code data code in
        List.filter (fun (_, a) -> not (Sequence.is_empty a gap)) taxa 
    in
    let partitions = Splitting.partition dhs.tcm2d taxa in
    match partitions with
    | (_, ((_ :: _) as h)) :: _ ->  
        let res = Hashtbl.create 1667 in
        let size = 1 + (List.length h) in
        List.iter 
            (fun (code, part) -> Hashtbl.replace res code part)
            partitions;
        Hashtbl.replace data.character_specs code 
            (Dynamic { dhs with 
                initial_assignment = `AutoPartitioned (mode, size, res)})
    | _ -> 
        Status.user_message Status.Information "There are no potential partitions"


let compute_fixed_states filename data code polymph =
    let debug = false and debug2 = false in
    if debug then Printf.printf "Data.compute_fixed_states, code=%d \n%!" code;
    let dhs = match Hashtbl.find data.character_specs code with
        | Dynamic dhs -> dhs
        | _ -> assert false
    in
    let annotate_with_mauve = match dhs.pam.annotate_tool with
        | Some (`Mauve _) -> true
        | Some (`Default _) -> false
        | None -> false
    and using_likelihood = match dhs.lk_model with
        | Some _ -> true
        | None   -> false
    in
    let align_with_newkk = match dhs.pam.align_meth with
        | Some `NewKK   -> true
        | None          -> false
        | Some `Default -> false
    in
    let taxon_sequences = Hashtbl.create 1667 in
    let sequences_taxon = Hashtbl.create 1667 in
    let states = ref 0 in
    let process_taxon tcode chars acc =
        try let (tc, _) = Hashtbl.find chars code in
            match tc with
            | Dyna (_, c) -> (tcode, (c.seq_arr.(0)).seq)  :: acc
            | _ -> failwith "Impossible?"
        with | Not_found -> acc
    in
    let taxa_arr = 
        Array.of_list
            (Hashtbl.fold process_taxon data.taxon_characters [])
    in
    let searchbase_arr =
        Array.of_list
            (Hashtbl.fold process_taxon data.searchbase_characters [])
    in
    let taxalen = (Array.length taxa_arr) 
    and searchbaselen = (Array.length searchbase_arr) in
    let total_arr = Array.append taxa_arr searchbase_arr in
    let taxa = Array.map fst total_arr in
    let initial_sequences = Array.map snd total_arr in
    if debug then Printf.printf "total size = %d + %d = %d\n%!" 
    taxalen searchbaselen  (Array.length total_arr);
    (* find all single assignments between two sequences; these become the
       states that can be placed on the internal nodes of the tree.
       this process resolves polymorphism of input sequences.
       NOTE: this process cost n-square(n=number of taxons),it might take long time.*)
    let () = 
        for x = 0 to (Array.length taxa) - 1 do
            for y = 0 to (Array.length taxa) - 1 do
                let a,b = 
                match polymph with
                | `Do_All -> if debug then Printf.printf "polymorphism=do all\n%!";
                let yclose, _ = 
                    if align_with_newkk then
                        Sequence.NewkkAlign.closest initial_sequences.(x)
                        initial_sequences.(y) dhs.tcm2d Sequence.NewkkAlign.default_ukkm
                    else
                    Sequence.Align.closest 
                        initial_sequences.(x) initial_sequences.(y) dhs.tcm2d Matrix.default 
                in
                let xclose, cost =
                    if align_with_newkk then
                        Sequence.NewkkAlign.closest yclose initial_sequences.(x)
                        dhs.tcm2d Sequence.NewkkAlign.default_ukkm
                    else
                        Sequence.Align.closest yclose initial_sequences.(x) dhs.tcm2d Matrix.default
                in
                xclose,yclose
                | `Pick_One ->
                        if debug then Printf.printf "polymorphism=pick_one\n%!";
                        Sequence.to_single_seq initial_sequences.(x) dhs.tcm2d,
                        Sequence.to_single_seq initial_sequences.(y) dhs.tcm2d
                | `Do_Nothing ->
                        if debug then Printf.printf "polymorphism=do_nothing\n%!";
                        initial_sequences.(x),initial_sequences.(y)
                in
                (*we add more states through extra sequence file -- check out
                * search_base. these extra states are for internal nodes only, 
                * we cannot assign them to leaf nodes. that's why we don't
                * add them to hashtbl taxon_sequences.*)
                if x<taxalen && y<taxalen then begin
                Hashtbl.replace taxon_sequences taxa.(y) b;
                Hashtbl.replace taxon_sequences taxa.(x) a;
                end;
                if not (Hashtbl.mem sequences_taxon b) then begin
                    Hashtbl.replace sequences_taxon b !states;
                    incr states;
                end;
                if not (Hashtbl.mem sequences_taxon a) then begin
                    Hashtbl.replace sequences_taxon a !states;
                    incr states;
                end;
            done;
        done;
    in
    let states = !states in
    if debug then Printf.printf "states num = %d\n%!" states;
    let sequences = Array.init states (fun _ -> Sequence.create 1) in
    let distances = Array.make_matrix states states 0. in
    let branches = ref None in
    Hashtbl.iter (fun seq pos -> sequences.(pos) <- seq) sequences_taxon;
    (* Fill the costs for all pairs of the single assignment sequences *)
    for x = 0 to states - 1 do
        for y = x + 1 to states - 1 do
            if debug then begin
                        Printf.printf "work on seqx,seqy=\n%!";
                        if debug2 then begin
                            Sequence.printseqcode sequences.(x);
                            Sequence.printseqcode  sequences.(y);
                        end;
                    end;
            let cost =
                if annotate_with_mauve then
                    let min_lcb_ratio,min_cover_ratio,min_lcb_len,max_lcb_len = 
                        match dhs.pam.annotate_tool with
                        | Some (`Mauve (a,b,c,d)) -> a,b,c,d
                        | _ -> assert(false)                        
                    in
                    let l_i_c = match dhs.pam.locus_indel_cost with
                        | Some cost -> cost
                        | None -> (10,100)
                    in
                    (*why delete first char?*)
                    let seqx,seqy = Sequence.del_first_char sequences.(x),
                    Sequence.del_first_char sequences.(y) in 
                    (*we don't call AliMap.create_general_ali_mauve directly here,
                    * because parameter of that function is with low level module
                    * type, like "Block.pairChromPam_t". Data.ml is suppose to
                    * be on top of those modules -- including aliMap. To avoid 
                    * circular denpendency, we call get_matcharr_and_costmatrix and
                    * output2mauvefile seperately.*)
                    let code1_arr,code2_arr,gen_cost_mat,ali_mat,gen_gap_code,
                            edit_cost,indel_cost,full_code_lstlst =
                        Block_mauve.get_matcharr_and_costmatrix seqx seqy
                                min_lcb_ratio min_cover_ratio min_lcb_len
                                max_lcb_len l_i_c dhs.tcm2d align_with_newkk
                    in
                    if debug then begin 
                        Printf.printf "code1/code2 arr from block_mauve:\n%!";
                    Utl.printIntArr code1_arr; 
                    Utl.printIntArr code2_arr; 
                    end;
                    let re_meth = match dhs.pam.re_meth with
                        | Some value -> value
                        | None -> `Locus_Breakpoint 10
                    and circular = match dhs.pam.circular with
                        | Some value -> value 
                        | None -> 0 
                    in                    
                    let cost, rc, alied_gen_seq1, alied_gen_seq2 = 
                    GenAli.create_gen_ali_new code1_arr code2_arr gen_cost_mat 
                    gen_gap_code re_meth circular false in
                    (*remember the editing cost between lcbs is not included in
                    * the gen_cost_mat, we set cost between matching lcb to 0 in
                    * block_mauve, just to make sure they are aligned to each other*)
                    if debug then begin
                        Printf.printf "cost <- cost(%d,rc=%d)+edit_cost(%d),\
                        alied_code1/code2=\n%!" cost rc edit_cost;
                        Utl.printIntArr alied_gen_seq1;
                        Utl.printIntArr alied_gen_seq2;
                    end;
                    let cost = cost + edit_cost in 
                    let xname,yname = string_of_int x,string_of_int y in
                    let fullname = match filename with 
                        | None -> ""
                        | Some fname -> (fname^"_"^xname^"_"^yname^".alignment")
                    in
                    Block_mauve.output2mauvefile fullname cost None 
                        alied_gen_seq1 alied_gen_seq2 full_code_lstlst ali_mat
                        gen_gap_code (Sequence.length seqx) (Sequence.length seqy);
                    float_of_int cost 
                else if using_likelihood then
                    let model =
                        match dhs.lk_model with | Some x -> x | None -> assert false
                    in
                    let bl,cst =
                        FloatSequence.pair_distance
                            (FloatSequence.make_model dhs.alph model)
                            sequences.(x) sequences.(y)
                    in
                    let () = match !branches with
                        | Some z ->
                            z.(x).(y) <- bl; z.(y).(x) <- bl
                        | None ->
                            let z = Array.make_matrix states states 0.0 in
                            z.(x).(y) <- bl; z.(y).(x) <- bl; 
                            branches := Some z
                    in
                    cst
                else
                    let cost = 
                    if align_with_newkk then
                    Sequence.NewkkAlign.cost_2 sequences.(x) sequences.(y)
                    dhs.tcm2d Sequence.NewkkAlign.default_ukkm
                    else
                    Sequence.Align.cost_2 sequences.(x) 
                    sequences.(y) dhs.tcm2d Matrix.default
                    in
                    float_of_int cost
            in
            if debug then 
                Printf.printf "set dis.%d.%d with cost %f\n%!" x y cost;
            distances.(x).(y) <- cost;
            distances.(y).(x) <- cost;
        done;
    done;
    let taxon_codes = Hashtbl.create 97 in
    Hashtbl.iter
        (fun code seq ->
            let tmp = Hashtbl.find sequences_taxon seq in
            Hashtbl.replace taxon_codes code tmp;
            (*
            seq --> Hashtbl.find sequences_taxon
                --> Hashtbl.replace taxon_codes code*)
        )
        taxon_sequences;
    let fs_data = 
        { costs = distances;
           seqs = sequences;
          codes = taxon_codes;
        opt_bls = !branches; }
    in
    let dyn_data =
        { dhs with
            initial_assignment = `FS fs_data;
            state = `Seq;
            polymorphism = polymph }
    in
    Hashtbl.replace data.character_specs code (Dynamic dyn_data)

(* make sure we call Data.categorize (Data.remove_taxa_to_ignore data) before
     calling this function to update tcm. *)
let assign_tcm_to_characters data chars foname tcm newalph =
    (* Get the character codes and filter those that are of the sequence class.
    * This allows simpler specifications by the users, for example, even though
    * morphological characters are loaded, an (all, create_tcm:(1,1)) will
    * operate properly in all the characters that are valid in the context. *)
    let per_size = Hashtbl.create 97 in
    let data = duplicate data in
    (*let chars = get_chars_codes_comp data chars in
    Printf.printf "dynamics list size = %d\n%!" (List.length chars);
     why we only update tcm file for dynamic character type?
    * *)let chars =  
        List.filter (fun x -> (List.exists (fun y -> x = y) data.dynamics))
                    (get_chars_codes_comp data chars)
    in(**)
    let chars_specs =
        List.fold_left 
        ~f:(fun acc x -> 
            let res = Hashtbl.find data.character_specs x in
            let acc = (res, x) :: acc in
            Hashtbl.remove data.character_specs x;
            acc
        ) 
        ~init:[] chars
    in
    let ref_tcmfile = ref None in
    let new_charspecs = 
        List.map 
            (function ((Dynamic dspec), code) ->
                let tcm, tcm3, tcmfile =
                    if Hashtbl.mem per_size dspec.alph then 
                        Hashtbl.find per_size dspec.alph
                    else 
                        let all_elements = 
                            if dspec.alph = Alphabet.nucleotides then 31
                            else if (Alphabet.is_aminoacids dspec.alph) then 21
                            else (-1)
                        in
                        let tcm, tcmfile = tcm all_elements in
                        ref_tcmfile := Some tcmfile;
                        if debug_level then 
                            Printf.printf "assign_tcm_to_characters,calc tcm3d if init3D=true\n%!";
                        let tcm3d =
                            if (Alphabet.use_3d dspec.alph) then
                                Cost_matrix.Three_D.of_two_dim tcm
                            else
                                dspec.tcm3d
                        in
                        tcm, tcm3d, tcmfile
                in
                (Dynamic { dspec with tcm = tcmfile; tcm2d = tcm; tcm3d = tcm3; alph = newalph}), 
                code
                | _, code -> raise (Invalid_Character code))
            chars_specs
    in
    let files = 
        match !ref_tcmfile with
        | Some (Input_file_GapOpening (tcmfile, _, _)) 
        | Some (Input_file (tcmfile, _)) -> 
                if List.exists (function (x, _) -> x = tcmfile) data.files 
                then data.files
                else (tcmfile, [CostMatrix]) :: data.files 
        | _ -> data.files
    in
    List.iter ~f:(fun (spec, code) -> 
                    Hashtbl.replace data.character_specs code spec)
              new_charspecs;
    List.iter ~f:(fun (spec, code) -> 
                    if is_fs data code then begin
                        let polymph = match spec with 
                        | Dynamic dspec ->
                                dspec.polymorphism
                        | _ -> `Do_All
                        in
                        compute_fixed_states None data code polymph;
                    end;
    )  new_charspecs;
    { data with files = files }


let assign_tcm_to_characters_from_file data chars file =
    let alphabet = 
        match get_chars_codes_comp data chars with
        | x::_ -> get_alphabet data x
        | []   -> failwith "No characters selected in transform"
    in
    let is_nucleotides = if alphabet=Alphabet.nucleotides then true else false in
    let is_aminoacids = Alphabet.is_aminoacids alphabet in
    let is_dna = if alphabet = Alphabet.dna then true else false in
    let is_dna_or_ami_or_nucleotides = (is_dna || is_aminoacids || is_nucleotides) in
    let oldlevel,ori_sz = 
        Alphabet.get_level alphabet,  Alphabet.get_ori_size alphabet in
    let tcm,newalph= match file with
        | None ->
            (fun x -> Cost_matrix.Two_D.default, default_tcm),
            alphabet
        | Some (f,level_and_tie_breaker) ->
            let level,tie_breaker,use_comb =
                match level_and_tie_breaker with
                | None -> 
                        if is_dna||is_nucleotides then 0,`Keep_Random,true
                        else if (Alphabet.check_level alphabet) then
                            oldlevel,`Keep_Random,true
                        else if oldlevel=1 then oldlevel,`Keep_Random,false
                        else 0,`Keep_Random,false
                | Some (l,tb) ->
                        if is_dna then 0,tb,true
                        else if l<=1 then l,tb,false
                        else if l>ori_sz then ori_sz,tb,true
                        else l,tb,true in
            (fun x ->
                if debug_level then Printf.printf
                "assign_tcm_to_characters_from_file,ori_sz=%d,oldlevel=%d,\
                newlevel=%d,use_comb=%b,is_dna?%b,is_ami?%b,is_nucl?%b\n%!"
                ori_sz oldlevel level use_comb is_dna is_aminoacids is_nucleotides;
                if debug_level then Alphabet.print alphabet;
                let tcm,mat = 
                        Cost_matrix.Two_D.of_file ~tie_breaker:tie_breaker ~use_comb:use_comb 
                        ~level:level f x is_dna_or_ami_or_nucleotides 
                in
                tcm, Input_file ((FileStream.filename f), mat)),
                Alphabet.set_level alphabet level
    in
    assign_tcm_to_characters data chars None tcm newalph 

let add_search_base_for_one_character_from_file data chars file character_name =
    let chcode = 
        try (Hashtbl.find data.character_names character_name)
        with | Not_found -> failwith "cannot add search base for non existing character"
    in
    if debug_search_base then 
        Printf.printf "Data.add_search_base_from_file %s, chcode=%d\n%!" file chcode;
    let original_filename = file in
    let file = `Local file in 
    let alphabet = get_alphabet data chcode in
    let ch = Parser.Files.molecular_to_fasta file in
    let res = 
        try Fasta.of_channel ~respect_case:true (FileContents.AlphSeq alphabet) ch with
        | Fasta.Illegal_molecular_format fl ->
                let file = FileStream.filename file in
                let fl = { fl with Fasta.filename = file } in
                print_error_message fl;
                raise (Fasta.Illegal_molecular_format fl)
    in
    let res = List.filter (function [[]], _ | [], _ -> false | _ -> true) res
    in
    close_in ch;
    let data = duplicate data in (*really?*)
    let res = 
        let res = List.map (fun (res3, b) -> (List.flatten res3),b) res in 
        let tmp = 
            List.map (fun (res2, b) -> 
                          let res1 = List.flatten res2 in
                          List.map (fun s -> s, b) res1) res
        in
        let arr = 
            let lst = List.map (Array.of_list) tmp in
            Array.of_list lst 
        in
        let num_taxa = Array.length arr in
        let loci =
            Array.fold_left 
                ~f:(fun max_loci loci_arr -> 
                     max max_loci (Array.length loci_arr))
                ~init:0 arr 
        in
        if debug_search_base then Printf.printf "num_taxa=%d,loci=%d\n%!" num_taxa loci;
        let res = ref [] in
        for j = loci - 1 downto 0 do 
            let loci = ref [] in
            for i = num_taxa - 1 downto 0 do
                match j < Array.length arr.(i) with
                | true -> loci := arr.(i).(j) :: !loci
                | false ->
                      let _, taxon = arr.(i).(0) in 
                      let seq = Sequence.create 1 in
                      let seq = 
                          Sequence.prepend_char seq (Alphabet.get_gap alphabet)
                      in
                      loci := (seq, taxon ) :: !loci
            done;
            res := !loci :: !res;
        done; 
        !res
    in
    (*get the old code*)
    (*let chcode = !(data.character_code_gen) in*)
    let locus_name = 
        let c = ref (-1) in
        ref (fun () ->
            incr c;
            original_filename  ^ ":" ^ string_of_int !c)
    in
    let single_loci_processor acc res = 
        let debug_search_base = false in
        if debug_search_base then Printf.printf "single_loci_processor\n%!";
        (*let data = add_searchbase_to_character acc chcode in*) 
        (* Now a function to process one taxon at a time to be 
        * folded over the taxa collected in the [file]. *)
        let process_a_taxon data (seq, taxon) =
            if debug_search_base then Printf.printf "process_a_taxon --> %!";
            (* Here is where, at the parser level, the fixed 
            * states support should be added *)
            let data, tcode = 
                process_searchbase_code data taxon original_filename
                (*process_taxon_code data taxon original_filename*) 
            in
            if debug_search_base then Printf.printf "tcode=%d, %!" tcode;
            if debug_search_base then 
                Array.iter (fun x -> Sequence.printseqcode x) seq;
            let tl = get_searchbase_characters data tcode in
            let dyna_state = `Seq and prealigned = false in
            let seqa = 
                let makeone seqa = {seq = seqa; delimiter = []; code = -1} in
                match dyna_state with 
                | `Ml  when not prealigned -> Array.map makeone seq 
                | `Seq when not prealigned -> Array.map makeone seq
                | _ -> Array.map (fun x -> x --> 
                        Sequence.del_first_char --> makeone) seq 
            in 
            let dyna_data = {seq_arr =  seqa} in 
            let () = 
                Hashtbl.replace tl chcode (Dyna (chcode, dyna_data), `Specified) 
            in
            Hashtbl.replace data.searchbase_characters tcode tl;
            data
        in
        List.fold_left ~f:process_a_taxon ~init:data res
    in
    let individual_fragments x = 
        List.map (List.map ~f:(fun (seq, t) -> [|seq|], t)) x
    in
    let data = 
       (* if annotated then process_annotated_chrom data 
        else if dyna_state = `Genome then process_genome data
        else if `DO = default_mode then *)
            match (res --> individual_fragments) with
            | [x] ->
                    locus_name := (fun () -> original_filename);
                    Printf.printf "only one fragments,call single_loci_processor\n%!";
                    single_loci_processor data x
            | [] -> 
                    Printf.printf "data remain unchanged\n%!";
                    data
            | x -> 
                    Printf.printf "fold on single_loci_processor (len=%d) \n%!" (List.length x);
                    List.fold_left ~f:single_loci_processor ~init:data x
       (* else 
            (res --> merge_fragments -->
            single_loci_processor data) *)
    in 
    
    let sbfile = FileStream.filename file in
    let files = 
            if List.exists (function (x, _) -> x = sbfile) data.files 
                then data.files
            else (sbfile, [Characters]) :: data.files 
    in
    let chars = 
        List.filter (fun x -> (List.exists (fun y -> x = y) data.dynamics))
                    (get_chars_codes_comp data chars)
    in
    let chars_specs =
        List.fold_left 
        ~f:(fun acc x -> 
            let res = Hashtbl.find data.character_specs x in
            let acc = (res, x) :: acc in
            (*Hashtbl.remove data.character_specs x;*)
            acc
        ) 
        ~init:[] chars
    in
    List.iter ~f:(fun (spec, code) -> 
        if debug_search_base then Printf.printf "call compute fixed states with code = %d\n%!" code;
        let polymph = match spec with 
        | Dynamic dspec ->
                dspec.polymorphism
        | _ -> `Do_All in
        compute_fixed_states None data code polymph) chars_specs;
    { data with files = files }


let add_search_base_from_file data chars file_chname_lst = 
    List.fold_left ~f:(fun accdata name ->
        let chname,filename = name in
        let newdata = add_search_base_for_one_character_from_file accdata chars filename
        chname in
        newdata
    ) ~init:data file_chname_lst 


let classify_characters_by_alphabet_size data chars =
    let is_dynamic_character x = 
        (List.exists (fun y -> x = y) data.dynamics)
    in
    let make_tuple_of_character_and_size acc char =
        let size = 
            data 
                --> get_sequence_alphabet char
                --> Alphabet.to_sequential 
                --> Alphabet.distinct_size
        in
        (char, size) :: acc
    in
    let classify_by_size list =
        let sets = 
            List.fold_left ~f:(fun acc (code, size) ->
                if All_sets.IntegerMap.mem size acc then
                    let prev = All_sets.IntegerMap.find size acc in
                    All_sets.IntegerMap.add size (code :: prev) acc
                else
                    All_sets.IntegerMap.add size [code] acc)
            ~init:All_sets.IntegerMap.empty list
        in
        All_sets.IntegerMap.fold (fun a b acc -> (a, `Some (true, b)) :: acc)
        sets []
    in
    chars 
        --> get_chars_codes_comp data 
        --> List.filter ~f:(is_dynamic_character)
        --> List.fold_left ~f:make_tuple_of_character_and_size ~init:[]
        --> classify_by_size

let assign_transformation_gaps data chars transformation gaps =
    let alphabet = match get_chars_codes_comp data chars with
        | x::_ -> get_alphabet data x
        | []   -> failwith "No characters selected in transform"
    in
    let uselevel = Alphabet.check_level alphabet in
    if uselevel then begin
        output_error "We are using level, input tcm file instead";
        data
    end else 
        let alphabet_sizes = classify_characters_by_alphabet_size data chars in
        List.fold_left 
            ~f:(fun data (size, chars) ->
                let tcm = 
                    (fun x -> 
                        Cost_matrix.Two_D.of_transformations_and_gaps 
                                        (size < 7) size transformation gaps x, 
                        (Substitution_Indel (transformation, gaps)))
                in
                assign_tcm_to_characters data chars None tcm alphabet)
            ~init:data alphabet_sizes

let codes_with_same_tcm codes data =
    (* This function assumes that the codes have already been filtered by class *)
    let rec assign_matching acc ((code : int), tcm, alph, tcmfile) =
        match acc with
        | (codes, assgn, alph,tcmfile) :: tl when tcm == assgn ->
                ((code :: codes), assgn, alph, tcmfile) :: tl
        | hd :: tl ->
                hd :: (assign_matching tl (code, tcm, alph, tcmfile))
        | [] -> [([code], tcm, alph, tcmfile)]
    in
    let codes = 
        List.map 
            ~f:(fun x -> 
                    x, get_tcm2d data x, get_alphabet data x,get_tcmfile data x)
            codes
    in
    List.fold_left ~f:assign_matching ~init:[] codes


let assign_level data chars tie_breaker level =
    if debug_level then Printf.printf "Data.assign_level,level=%d\n%!" level;
    let make_level level = function
        | Level (otcm,_) -> Level (otcm,level)
        | x -> Level (x, level)
    in
    let codes = get_chars_codes_comp data chars in
    let codes =
        List.map 
          (fun (a, b, alph, tcmfile) ->
            if debug_level then begin
                Printf.printf "before assigning:%!"; Alphabet.print alph;
            end;
            let name = make_level level tcmfile in
            let b =
                if (Alphabet.dna = alph) || (Alphabet.nucleotides = alph) then
                    b
                else
                    let all_elements = Cost_matrix.Two_D.get_all_elements b in
                    let oldlevel = Cost_matrix.Two_D.get_level b in
                    let ori_sz = Cost_matrix.Two_D.get_ori_a_sz b in
                    let pure_sz = if all_elements>0 then ori_sz-1 else ori_sz in
                    let combnum =
                        if level > 1 then
                            Cost_matrix.Two_D.calc_number_of_combinations_by_level
                            pure_sz level
                        else
                            ori_sz
                    in
                    if debug_level then Printf.printf
                    "Data.assign_level,oldlevel=%d,newlevel=%d\n%!" oldlevel level;
                    if combnum <= 0 then begin
                        output_info ("The alphabet size based on the new level is"^
                                     " too large. I will NOT apply any changes.\n%!");
                        b
                    end else begin
                        let b = Cost_matrix.Two_D.clone b in
                        let b = Cost_matrix.Two_D.create_cm_by_level b level
                        oldlevel all_elements tie_breaker in
                        b
                    end
                in
                (true, a),(fun _ -> b,name),Alphabet.set_level alph level)
            (codes_with_same_tcm codes data)
    in
    List.fold_left 
        ~f:(fun acc (a, tcm, newalph) -> assign_tcm_to_characters acc (`Some a) None tcm newalph)
        ~init:data codes

let rec make_affine cost_model tcmfile = match tcmfile with
    | Substitution_Indel (a, b) ->
            (match cost_model with
            | Cost_matrix.No_Alignment
            | Cost_matrix.Linnear
            | Cost_matrix.Affine 0 -> tcmfile
            | Cost_matrix.Affine x -> Substitution_Indel_GapOpening (a, b, x))
    | Input_file (a, b) ->
            (match cost_model with
            | Cost_matrix.No_Alignment
            | Cost_matrix.Linnear
            | Cost_matrix.Affine 0 -> tcmfile
            | Cost_matrix.Affine x -> Input_file_GapOpening (a, b, x))
    | Substitution_Indel_GapOpening (a, b, _) ->
            (match cost_model with
            | Cost_matrix.No_Alignment
            | Cost_matrix.Linnear
            | Cost_matrix.Affine 0 -> Substitution_Indel (a ,b)
            | Cost_matrix.Affine x -> Substitution_Indel_GapOpening (a, b, x))
    | Input_file_GapOpening (a, b, _) ->
            (match cost_model with
            | Cost_matrix.No_Alignment
            | Cost_matrix.Linnear
            | Cost_matrix.Affine 0 -> tcmfile
            | Cost_matrix.Affine x -> Input_file_GapOpening (a, b, x))
    | Level (inner,n) -> 
            Level (make_affine cost_model inner,n)

let rec assign_affine_gap_cost data chars cost =
    let codes = get_code_from_characters_restricted_comp `AllDynamic data chars in 
    let codes = 
        List.map 
            (fun (a, b, alph, tcmfile) -> 
                let b =
                    if Alphabet.nucleotides = alph then
                        let b = Cost_matrix.Two_D.clone b in
                        let () = Cost_matrix.Two_D.set_affine b cost in
                        b
                    else 
                        b
                in 
                (true, a), (fun _ -> b, make_affine cost tcmfile),alph)
            (codes_with_same_tcm codes data)
    in
    let cost = 
        match cost with
        | Cost_matrix.Linnear
        | Cost_matrix.No_Alignment -> "0"
        | Cost_matrix.Affine x -> string_of_int x
    in
    List.fold_left 
        ~f:(fun acc (a, b, alph) ->
            assign_tcm_to_characters acc (`Some a) (Some cost) b alph)
        ~init:data codes

let rec assign_prep_tail filler data chars filit =
    match filit with
    | `File x ->
            let ch = new FileStream.file_reader x in
            let lst = Cost_matrix.Two_D.load_file_as_list ch in
            ch#close_in;
            let arr = Array.of_list lst in
            assign_prep_tail filler data chars (`Array arr)
    | `Array arr ->
            let codes = 
                get_code_from_characters_restricted_comp `AllDynamic data chars
            in
            let codes = codes_with_same_tcm codes data in
            let codes =
                List.map
                    (fun (a, b, alph,tcmfile) -> 
                        let b = Cost_matrix.Two_D.clone b in
                        let () = filler arr b in
                        (true, a), (fun _ -> b,tcmfile),alph)
                    codes
            in
            List.fold_left 
                ~f:(fun acc (a, b, alph) ->
                    assign_tcm_to_characters acc (`Some a) None b alph) 
                ~init:data codes

let assign_prepend data chars filit =
    assign_prep_tail Cost_matrix.Two_D.fill_prepend data chars filit

let assign_tail data chars filit =
    assign_prep_tail Cost_matrix.Two_D.fill_tail data chars filit

(** [process_complex_terminals data filename] reads a complex terminals file
    ([filename]) and adds the specification to [data].  It will raise
    exceptions if the reading fails. *)
let process_complex_terminals data filename =
    let groups = Parser.SetGroups.of_file filename in
    if groups <> []
    then begin
            let group_type = Parser.SetGroups.unify_list groups in
            let groups =
                List.map (Parser.SetGroups.coerce_to_type group_type)
                    groups in
            (* create ids and types *)
            { data with complex_schema = groups }
        end
    else { data with complex_schema = [] }

let get_character_state data c =
    match Hashtbl.find data.character_specs c  with
    | Dynamic dspec -> dspec.state
    | Kolmogorov dspec -> dspec.dhs.state
    | _ -> failwith "Data.get_character_state"

let get_pam data c =
    match Hashtbl.find data.character_specs c  with
    | Dynamic dspec -> dspec.pam
    | Kolmogorov dspec -> dspec.dhs.pam
    | _ -> failwith "Data.get_pam"

let report_taxon_file_cross_reference chars data filename =
    let files_arr, taxa = 
        match chars with
        | None -> (* The user requested the data per file *)
                let files_arr = 
                    let filtered_files = 
                        let has_char (_, c) =
                            List.exists (function Characters -> true | _ ->
                                false) c
                        in
                        let filtered = List.filter has_char data.files in
                        List.map (fun (a, _) -> a) filtered
                    in
                    Array.of_list ("Terminal" :: filtered_files) in
                let taxa =
                    All_sets.StringMap.fold 
                     (fun taxon files arr ->
                         if All_sets.Strings.mem taxon data.ignore_taxa_set then
                             arr
                         else
                             let new_arr = Array.mapi (fun pos file ->
                                 if pos = 0 then taxon
                                 else if All_sets.Strings.mem file files then "+"
                                 else "-") files_arr
                             in
                             new_arr :: arr) data.taxon_files []
                in
                files_arr, taxa
        | Some chars -> (* The user requested the data for some characters *)
                let is_specified code_taxon code_char data = 
                    (* Establish if a character was parte of the input of a
                    * taxon or not *)
                    let taxon_specs = 
                        Hashtbl.find data.taxon_characters code_taxon
                        
                    in
                    try 
                        match Hashtbl.find taxon_specs code_char with
                        | (_, `Specified) -> true
                        | (_, `Unknown) -> false
                    with
                    | Not_found -> false
                in
                let codes = get_chars_codes_comp data chars in
                let codes = List.sort Pervasives.compare codes in
                let codes_arr = Array.of_list codes 
                and chars_arr = 
                    let name x = StatusCommon.escape (Hashtbl.find
                    data.character_codes x) in
                    Array.of_list ("Terminal" :: List.map name codes)
                in
                let taxa = 
                    All_sets.IntegerMap.fold
                    (fun code taxon acc ->
                        try
                            let new_arr = Array.mapi (fun pos file ->
                                if pos = 0 then taxon
                                else if 
                                    is_specified code (codes_arr.(pos - 1)) data
                                then "+"
                                else "-") chars_arr
                            in
                            new_arr :: acc
                        with Not_found -> acc) data.taxon_codes []
                in
                chars_arr, taxa
    in
    let files_arr = Array.map (fun x ->
        "@{<u>" ^ StatusCommon.escape (Filename.basename x) ^ "@}") files_arr
    in
    let fo = Status.Output (filename, false, []) in
    Status.user_message fo "@[<v 2>@{<b>File References:@}@,@[";
    Status.output_table fo 
    (Array.of_list (files_arr :: (List.rev taxa)));
    Status.user_message fo "@]@]@."



let find_max_seq_id data = 
    let max_seq_id = Hashtbl.fold  
        (fun key cs_ls max_seq_id ->
            Hashtbl.fold
            (fun _ cs max_seq_id -> 
                let csd, _ = cs in 
                match csd with
                | Dyna (_, dyna_data) ->
                        Array.fold_left 
                        ~f:(fun max_seq_id seq -> max max_seq_id seq.code
                ) ~init:max_seq_id dyna_data.seq_arr 
                | _ -> max_seq_id) 
            cs_ls 
            max_seq_id) 
        data.taxon_characters 
        0  
    in 
    max_seq_id + 2



let flush d = 
    Hashtbl.iter (fun _ item ->
        match  item with
        | Dynamic dspec -> ()
        | _ -> ()) d.character_specs 

let set_weight weight spec =
    match spec with
    | Static enc ->
            Static { enc with Nexus.File.st_weight = weight}
    | Dynamic enc ->
            Dynamic ({ enc with weight = weight })
    | Kolmogorov enc ->
            let res = { enc.dhs with weight = weight } in
            Kolmogorov { enc with dhs = res }
    | _ -> spec

let set_weight_factor weight spec =
    match spec with
    | Static enc ->
            Static 
            { enc with Nexus.File.st_weight = enc.Nexus.File.st_weight *. weight }
    | Dynamic enc ->
            Dynamic ({ enc with weight = enc.weight *. weight })
    | Kolmogorov enc ->
            let res = { enc.dhs with weight = enc.dhs.weight *. weight } in
            Kolmogorov { enc with dhs = res }
    | _ -> spec

let aux_transform_weight meth data =
    let f = 
        match meth with
        | `ReWeight (chars, weight) ->
                let chars = 
                    let codes = get_chars_codes_comp data chars in
                    List.fold_left ~f:(fun acc x -> All_sets.Integers.add x acc)
                    ~init:All_sets.Integers.empty codes
                in
                (fun code spec ->
                    if All_sets.Integers.mem code chars then
                        set_weight weight spec
                    else spec)
        | `WeightFactor (chars, weight) ->
                let chars = 
                    let codes = get_chars_codes_comp data chars in
                    List.fold_left ~f:(fun acc x -> All_sets.Integers.add x acc)
                    ~init:All_sets.Integers.empty codes
                in
                (fun code spec ->
                    if All_sets.Integers.mem code chars then
                        set_weight_factor weight spec
                    else spec)
    in
    Hashtbl.iter (fun code char ->
        Hashtbl.replace data.character_specs code (f code char)) 
    data.character_specs

let transform_weight meth data = 
    let data = duplicate data in
    aux_transform_weight meth data;
    data

let file_exists data filename =
    List.exists (fun (x, _) -> x = FileStream.filename filename) data.files

let complement_taxa data taxa = 
    let taxa = 
        List.fold_left ~f:(fun acc x -> All_sets.Integers.add x acc)
        ~init:All_sets.Integers.empty taxa
    in
    All_sets.IntegerMap.fold (fun c _ acc -> 
        if All_sets.Integers.mem c taxa then acc
        else c :: acc) data.taxon_codes []

let make_fixed_states filename chars polymph data =
    let data = duplicate data in
    let polymph =
        match polymph with 
        | None -> `Do_All 
        | Some x -> x 
    in
    let convert_and_process code =
        match Hashtbl.find data.character_specs code with
        | Dynamic dhs -> 
                (match dhs.initial_assignment with
                | `FS _ -> ()
                | `Partitioned _
                | `AutoPartitioned _
                | `DO
                | `GeneralNonAdd -> 
                        compute_fixed_states filename data code polymph)
        | _ -> failwith "How could this happen?"
    in
    let codes = get_code_from_characters_restricted_comp `Dynamic data chars in
    List.iter ~f:convert_and_process codes;
    data

let make_direct_optimization chars data =
    let data = duplicate data in
    let convert_and_process code =
        match Hashtbl.find data.character_specs code with
        | Dynamic dhs ->
                (match dhs.initial_assignment with
                | `AutoPartitioned _
                | `Partitioned _
                | `FS _ -> 
                    Hashtbl.replace data.character_specs code 
                        (Dynamic { dhs with initial_assignment = `DO })
                | `DO | `GeneralNonAdd -> ())
        | _ -> ()
    in
    let codes = get_code_from_characters_restricted_comp `Dynamic data chars in
    List.iter ~f:convert_and_process codes;
    data

let make_partitioned mode chars data =
    let data = duplicate data in
    let convert_and_process code =
        match Hashtbl.find data.character_specs code with
        | Dynamic dhs -> 
                (match dhs.initial_assignment with
                | `FS _ -> ()
                | `AutoPartitioned _
                | `Partitioned _
                | `DO | `GeneralNonAdd -> auto_partition mode data code)
        | _ -> failwith "How could this happen?"
    in
    let codes = get_code_from_characters_restricted_comp `Dynamic data chars in
    List.iter ~f:convert_and_process codes;
    data


let number_of_taxa d = 
    Hashtbl.fold  (fun _ _ num_taxa -> num_taxa + 1) d.taxon_characters 0  

let has_dynamic d = 
    match d.dynamics, d.kolmogorov with
    | [], [] -> false
    | _      -> true


let can_do_static_approx_code d x =
    let appropriate_alphabet_size ds =
        10 > (Alphabet.distinct_size (Alphabet.to_sequential ds.alph))
    in
    match Hashtbl.find d.character_specs x with
        | Dynamic d when appropriate_alphabet_size d ->
            begin match d.state with
                | `Seq      | `Annotated  | `Ml                      -> true
                | `Breakinv | `Chromosome | `Genome | `SeqPrealigned -> false
            end
        (* only dynamics with alphabet < 10 *)
        | Dynamic d     -> false | Static _      -> false
        | Set           -> false | Kolmogorov _  -> false

let can_all_chars_do_static_approx d xs =
    List.fold_left ~f:(fun acc x -> acc && (can_do_static_approx_code d x)) ~init:true xs

let filter_non_static_approx_characters ?(comp=true) d xs =
    if comp then List.filter (can_do_static_approx_code d) xs
            else List.filter (fun x -> not (can_do_static_approx_code d x)) xs


(* [sync_model_branches copy translate src dest] sync the data from the src to
 * destination. Copy defines if the data returned is a copy, or if the
 * destination should be returned. The translate parameter will translate the
 * codes in the branch lengths with the map defined in the dynamic_static_codes
 * map. *)
let sync_dynamic_to_static_model ~src ~dest =
    let char_specs =
        let spec = Hashtbl.copy dest.character_specs in
        All_sets.IntegerMap.iter
            (fun src_c dest_cs -> match get_model_opt src src_c with
                | Some x ->
                    apply_likelihood_model_on_char_table false dest spec dest_cs x;
                    Hashtbl.remove spec src_c
                | None -> ())
            dest.dynamic_static_codes;
        spec
    in
    categorize { dest with character_specs = char_specs; }

(* converse of above function *)
let sync_static_to_dynamic_model ~src ~dest = 
    (* confirm a set is consistent and return model *)
    let get_model_from_set data codes = 
        get_model_opt data (All_sets.Integers.choose codes)
    in
    let char_specs = 
        let spec = Hashtbl.copy dest.character_specs in
        All_sets.IntSetMap.iter
            (fun src_cs dest_c -> 
                match get_model_from_set src src_cs with
                | Some x ->
                    apply_likelihood_model_on_char_table false dest spec [dest_c] x;
                    All_sets.Integers.iter (fun x -> Hashtbl.remove spec x) src_cs
                | None -> ())
            src.static_dynamic_codes;
        spec
    in
    categorize { dest with character_specs = char_specs; }


(** Functions to modify the taxon codes *)
let change_taxon_codes reorder_function data =
    let data = duplicate data in
    (* First we produce a hash table with an randomized reassignment of codes
     * for the taxa in data *)
    let htbl =
        let taxon_codes = 
            All_sets.StringMap.fold (fun name code acc -> 
                if All_sets.Strings.mem name data.ignore_taxa_set then acc
                else code :: acc) data.taxon_names []
        in
        let chars = Array.of_list taxon_codes 
        and chars_org = Array.of_list taxon_codes in
        reorder_function chars;
        let htbl = Hashtbl.create 1667 in
        for i = 0 to (Array.length chars_org) - 1 do
            Hashtbl.replace htbl chars_org.(i) chars.(i);
        done;
        htbl
    in
    (* Now that we have the assignment, we proceed to modify the contents of the
     * new data. *)
    let find htbl code =
        try Hashtbl.find htbl code with
        | Not_found -> code 
    in
    let taxon_names = 
        All_sets.StringMap.fold (fun name old_code acc ->
            All_sets.StringMap.add name (find htbl old_code) acc)
        data.taxon_names All_sets.StringMap.empty
    and taxon_codes =
        All_sets.IntegerMap.fold (fun old_code name acc ->
            All_sets.IntegerMap.add (find htbl old_code) name acc)
        data.taxon_codes All_sets.IntegerMap.empty
    and taxon_characters =
        let res = Hashtbl.create 1667 in
        Hashtbl.iter (fun old_code contents ->
            Hashtbl.replace res (find htbl old_code) contents)
        data.taxon_characters;
        res
    in
    let root = match data.root_at with
        | None -> None
        | Some code -> 
                try Some (Hashtbl.find htbl code) with
                | _ -> None
    in
    { data with taxon_names = taxon_names;
                taxon_codes = taxon_codes;
                taxon_characters = taxon_characters;
                root_at = root },
    htbl

let randomize_taxon_codes meth data = 
    let f = match meth with
        | `RandomizedTerminals -> Array_ops.randomize 
        | `AlphabeticTerminals ->
                Array.stable_sort ~cmp:(fun a b ->
                    String.compare (code_taxon a data) (code_taxon b
                    data)) 
    in
    change_taxon_codes f data

let lexicographic_taxon_codes data = 
    let lexicographic_sort (a : int) b =
        let namea = code_taxon a data
        and nameb = code_taxon b data in
        String.compare namea nameb
    in
    change_taxon_codes (Array.stable_sort ~cmp:lexicographic_sort) data

(* A function to produce the alignment of prealigned data *)
let process_prealigned analyze_tcm data code : (string * Nexus.File.nexus) =
    let alph = get_sequence_alphabet code data in
    let model = get_model_opt data code in
    let gap = Alphabet.get_gap alph in
    let character_name = code_character code data in
    let tcm_case, do_states, do_encoding = 
        let cm = get_sequence_tcm code data in
        analyze_tcm cm model alph
    in
    (* We first define a function that collects all the limits of the gaps so
    * that we can define the gap characters if the gap opening parameter is
    * present. We use bitsets, 1 for starting, 2 for ending (a position could be
    * both start and end). *)
    let mark_gaps sequences =
        if Array.length sequences > 0 then
            let first = sequences.(0) in
            let res = Array.make (Sequence.length first) 0 in
            res.(0) <- 0;
            let handle_start seq pos =
                if Sequence.get seq pos <> gap && 
                Sequence.get seq (pos + 1) = gap then
                    res.(pos + 1) <- res.(pos + 1) lor 1
                else ()
            and handle_end seq pos =
                if Sequence.get seq pos = gap && 
                Sequence.get seq (pos + 1) <> gap then
                    res.(pos) <- res.(pos) lor 1
                else ()
            in
            for i = 0 to (Array.length sequences) - 1 do
                let seq = sequences.(i) in
                handle_start seq 0;
                for j = 1 to (Sequence.length seq) - 2 do
                    if gap = Sequence.get seq j then
                        res.(j) <- res.(j) lor 2;
                    handle_start seq j;
                    handle_end seq j;
                done;
                handle_end seq ((Sequence.length seq) - 2);
            done;
            res
        else 
            [||]
    in
    (* Now we define a function to compute the length of the indel blocks *)
    let rec find_end acc start pos mark =
        if pos = Array.length mark then 
            ((start, (pos - start) + 1) :: acc)
        else if 0 = mark.(pos) then
            find_start ((start, (pos - start)) :: acc) pos mark
        else if 0 <> 1 lor mark.(pos) then
            find_end ((start, (pos - start)) :: acc) pos (pos + 1) mark
        else find_end acc start (pos + 1) mark
    and find_start acc pos mark = 
        if pos = Array.length mark then acc
        else if 0 = mark.(pos) then 
            find_start acc (pos + 1) mark
        else find_end acc pos (pos + 1) mark
    in
    (* A function to compute the cost of an indel block *)
    let compute_cost = 
        match tcm_case with
        | `AllSankoff None
        | `AllOne  _
        | `AllOneGapSame _ -> (fun _ -> 0)
        | `AllSankoff (Some f) -> 
                (fun len -> 
                    f (String.make len 'A'))
        | `AffinePartition (_, gapcost, gapopening) ->
                (fun len -> gapopening + (len * gapcost))
    in
    let encoding len = 
        Alphabet.present_absent, Parser.OldHennig.Encoding.gap_encoding (compute_cost len)
    in
    let make_indel_blocks_encoding lst = 
        let res = List.rev_map (fun (_, x) -> encoding x) lst in
        List.rev res
    in
    let make_indel_blocks lst seq = 
        let res = List.rev_map (fun (start, _) ->
            if gap = Sequence.get seq start then
                FileContents.Unordered_Character (2, false)
            else FileContents.Unordered_Character (1, false)) lst in
        List.rev res
    in
    let compute_blocks_of_indels () =
        let process_taxon a b acc =
            if Hashtbl.mem b code then
                match Hashtbl.find b code with
                | (Stat _), _
                | _, `Unknown -> acc
                | (Dyna (_, d)), `Specified ->
                    begin match d.seq_arr with
                        | [|v|] -> v.seq :: acc
                        | _ -> assert false
                    end
            else 
                acc
        in
        let sequences =
            Hashtbl.fold process_taxon data.taxon_characters [] 
        in
        let sequences = Array.of_list sequences in
        let gaps = mark_gaps sequences in
        let blocks = find_start [] 0 gaps in
        blocks
    in
    let blocks = 
        match tcm_case with
        | `AllSankoff (Some _)
        | `AffinePartition _ -> compute_blocks_of_indels ()
        | `AllOneGapSame _ 
        | `AllOne _
        | `AllSankoff None -> []
    in
    let process_taxon a b ((enc, names, acc) as res) =
        if Hashtbl.mem b code then
            match Hashtbl.find b code with
            | (Stat _), _ -> res
            | _, `Unknown -> res
            | (Dyna (_, d)), _ -> 
                begin match d.seq_arr with
                | [|v|] ->
                    let enc = match enc with
                        | [||] -> 
                            (* We have to generate the encoding *)
                            let initial_acc = match tcm_case with
                                | `AllSankoff (Some _)
                                | `AffinePartition _ -> make_indel_blocks_encoding blocks
                                | _ -> []
                            in
                            let res = 
                                Sequence.fold_right 
                                    (fun acc base -> do_encoding base acc)
                                    initial_acc v.seq
                            in
                            Array.of_list res
                        | _ -> enc
                    and seq =
                        let initial_acc = match tcm_case with
                            | `AllSankoff (Some _)
                            | `AffinePartition _ -> make_indel_blocks blocks v.seq
                            | _ -> []
                        in
                        Sequence.fold_right 
                            (fun acc base ->
                                let base = if base = gap then 0 else base in
                                do_states `Exists base acc)
                            initial_acc v.seq
                    and name = code_taxon a data in
                    let seq = Array.of_list seq in
                    if Array.length seq <> Array.length enc then begin
                        Status.user_message Status.Error
                        ("The@ prealigned@ sequences@ in@ " ^
                        code_character code data ^ "@ do@ not@ have@ \
                        the@ same@ length.@ The@ taxon@ " ^ name ^ "@ \
                        has@ a@ sequence@ of@ length@ " ^ string_of_int 
                        (Array.length seq) ^ "@ while@ the@ expected@ \
                        length@ is@ " ^ string_of_int (Array.length enc));
                        failwith "Illegal prealigned molecular sequences."
                    end;
                    enc, (Some name) :: names, seq :: acc
                | _ -> assert false
            end
        else 
            res
    in
    let enc, names, matrix = 
        Hashtbl.fold process_taxon data.taxon_characters ([||], [], [])
    in
    let newenc = 
        let table = Hashtbl.create 1667 in
        Array.iter 
            (fun ((alph, enc) as r) ->
                if not (Hashtbl.mem table r) then
                    let alph = Alphabet.to_sequential alph in
                    let alph = Parser.OldHennig.generate_alphabet (Some alph) enc in
                    Hashtbl.replace table r alph
                else
                    ())
            enc;
        Array.init 
            (Array.length enc) 
            (fun pos ->
                let alph = Hashtbl.find table enc.(pos) in
                Parser.OldHennig.to_new_spec character_name alph (snd enc.(pos)) pos)
    in
    let matrix = 
        let matrix = Array.of_list matrix in
        let table = Hashtbl.create 1667 in
        Array.init (Array.length matrix) 
            (fun x -> 
                Array.init (Array.length enc)
                    (fun y -> 
                        Parser.OldHennig.to_new_atom table newenc.(y) 
                                            (snd enc.(y)) matrix.(x).(y)))
    in
    let res =
        { Nexus.File.empty_parsed () with
            Nexus.File.taxa = Array.of_list names;
            characters = newenc;
            matrix = matrix;
        }
    in
    Nexus.File.fill_observed res.Nexus.File.characters res.Nexus.File.matrix;
    character_name, res

type tcm_class =
    [ `AllOne of int
    | `AllOneGapSame of (int * int)
    | `AffinePartition of (int * int * int)
    | `AllSankoff of (string -> int) option]

let prealigned_characters analyze_tcm data chars =
    let codes = get_chars_codes_comp data chars in
    let names = List.map (fun x -> code_character x data) codes in
    let res = List.rev_map (process_prealigned analyze_tcm data) codes in
    let d =
        List.fold_left
            ~f:(fun acc (file,triple) ->
                    let _,d = gen_add_static_parsed_file false acc file triple in
                    d)
            ~init:data
            res
    in
    process_ignore_characters false d (`Names names)

let sequence_statistics ch data =
    let codes = get_chars_codes_comp data ch in
    List.map (sequence_code_statistics data) codes

let compare_all_pairs char1 char2 complement data = 
    let alpha1 = get_sequence_alphabet char1 data
    and alpha2 = get_sequence_alphabet char2 data 
    and cm = get_sequence_tcm char1 data in
    let gap = Alphabet.get_gap alpha1 in
    if alpha1 = alpha2 then 
        let process_taxon a b (cnt, res) =
            match Hashtbl.find b char1, Hashtbl.find b char2 with
            | ((Stat _), _), _
            | ((_, `Unknown), _) | (_, (_, `Unknown)) -> (cnt, res)
            | ((Dyna (_, d)), `Specified), ((Dyna (_, e), `Specified)) ->
                    (match d.seq_arr, e.seq_arr with
                    | [|dv|], [|de|] ->
                            if Sequence.is_empty dv.seq gap || 
                            Sequence.is_empty de.seq gap then
                                (cnt, res)
                            else
                                let de = 
                                    if complement then 
                                        Sequence.complement alpha1 de.seq
                                    else de.seq
                                in
                                (cnt + 1), res +. 
                                ((float_of_int (Sequence.Align.cost_2 dv.seq de
                                cm Matrix.default)) /.
                                (float_of_int ((max (Sequence.length dv.seq) (Sequence.length
                                de)))))
                    | _ -> (cnt, res))
            | _ -> assert false
        in
        Some (Hashtbl.fold process_taxon data.taxon_characters (0, 0.))
    else None

let compare_pairs ch1 ch2 complement data =
    let rec process_pair acc a b =
        let namea = code_character a data
        and nameb = code_character b data in
        match compare_all_pairs a b complement data with
        | None -> acc
        | Some (cnt, sum) ->
                (namea, nameb, (sum) /. (float_of_int cnt)) :: acc
    in
    let codes1 = get_chars_codes_comp data ch1
    and codes2 = get_chars_codes_comp data ch2 in
    List.fold_left ~f:(fun acc x ->
        List.fold_left ~f:(fun acc y -> process_pair acc x y) ~init:acc codes1)
    ~init:[] codes2

module Sample = struct
    let characters_to_arr chars =
        Array.of_list (Hashtbl.fold (fun code spec acc ->
            (0, code, spec) :: acc) chars [])

    let bootstrap_spec ?(rand = Random.int) arr m =
        let n = Array.length arr in
        if n > 0 then begin
            for i = 0 to m do
                let p = rand n in
                let (cnt, code, spec) = arr.(p) in
                arr.(p) <- (cnt + 1, code, spec);
            done;
            arr
        end else arr

    (** [jackknife_spec ar m] resamples [Array.length ar - m] characters without
    * replacement uniformly at random *)
    let jackknife_spec ar prob =
        let n = Array.length ar in
        assert (prob <= 1.);
        for i = 0 to n - 1 do
            if prob < (Random.float 1.) then
                let (_, code, spec) = ar.(i) in
                ar.(i) <- (1, code, spec);
        done;
        ar

    let generate data perturb =
        let arr = characters_to_arr data.character_specs in
        let data = duplicate data in
        let arr =
            match perturb with
            | `Bootstrap -> 
                    bootstrap_spec arr (Array.length arr)
            | `Jackknife m ->
                    jackknife_spec arr m 
        in
        Array.iter (fun (cnt, code, spec) ->
            Hashtbl.replace data.character_specs code 
                    (set_weight_factor (float_of_int cnt) spec)) arr;
        data

end

let to_human_readable data code item =
    let spec = 
        match Hashtbl.find data.character_specs code with
        | Static x -> x
        | _ -> assert false
    in
    let item = 
        match spec.Nexus.File.st_type with
        | Nexus.File.STUnordered -> List.nth spec.Nexus.File.st_observed item
        | Nexus.File.STLikelihood _ | Nexus.File.STOrdered | Nexus.File.STSankoff _ -> item
    in
    let name =
        try 
            try List.nth spec.Nexus.File.st_labels item with
            | _ -> Alphabet.match_code item spec.Nexus.File.st_alph 
        with
            | (Alphabet.Illegal_Code n) as err ->
                Alphabet.print spec.Nexus.File.st_alph;
                raise err
    in
    name

let apply_boolean nonadd_f add_f data char = 
    match Hashtbl.find data.character_specs char with
    | Static x ->
            (* We can test only for non-additive
            * characters first *)
            let generate_characters_list () = 
                Hashtbl.fold (fun taxon chars acc ->
                    try 
                        match Hashtbl.find chars char with
                        | (Stat (_, states)), _  -> 
                                states :: acc
                        | _ -> assert false
                    with
                    | Not_found -> None :: acc)
                data.taxon_characters []
            in
            (match x.Nexus.File.st_type with
            | Nexus.File.STUnordered -> 
                    let all_specs = generate_characters_list () in
                    nonadd_f all_specs
            | Nexus.File.STOrdered -> 
                    let all_specs = generate_characters_list () in
                    add_f all_specs
            | _ -> true)
    | _ -> true

let get_model code data =
    match Hashtbl.find data.character_specs code with
    | Static x -> 
        (match x.Nexus.File.st_type with
            | Nexus.File.STLikelihood x -> x
            | _ -> failwith "Data.get_model")
    | _ -> failwith "Data.get_model 2"

(* We define a function that adds the min possible cost to the tree *)
let apply_on_static ordered unordered sankoff likelihood char data =
    let process_code acc code =
        match Hashtbl.find data.character_specs code with
        | Static spec ->
                let specified_static_chars = 
                    Hashtbl.fold (fun a b acc -> 
                        try match Hashtbl.find b code with
                        | (Stat (_, x), `Specified) -> x :: acc
                        | _ -> acc
                        with Not_found -> acc) data.taxon_characters []
                in
                let res = 
                    match spec.Nexus.File.st_type with
                    | Nexus.File.STOrdered -> 
                            ordered specified_static_chars
                    | Nexus.File.STUnordered -> 
                            unordered specified_static_chars
                    | Nexus.File.STSankoff mtx -> 
                            sankoff mtx specified_static_chars 
                    | Nexus.File.STLikelihood model -> 
                            likelihood model specified_static_chars 
                in
                (code, res) :: acc
        | _ -> acc
    in
    let codes = get_chars_codes_comp data char in
    List.fold_left ~f:process_code ~init:[] codes

let guess_class_and_add_file annotated is_prealigned data filename =
    if file_exists data filename then
        let () =
            let filename = FileStream.filename filename in
            let msg = 
                "@[A@ file@ with@ name@ " ^ StatusCommon.escape filename ^ 
                "@ has@ previously@ " 
                ^ "been@ loaded.@ Sorry,@ I@ will@ cowardly@ refuse@ to@ "
                ^ "load@ its@ contents@ again.@ However,@ I@ will@ continue@ "
                ^ "loading@ any@ files@ remaining.@]"
            in
            Status.user_message Status.Error msg
        in
        data
    else
        let file_type_message str = 
            let msg =
                let filename = FileStream.filename filename in
                "@[Reading@ file@ " ^ StatusCommon.escape filename ^ 
                "@ of@ type@ " ^ str ^ "@]@." 
            in
            Status.user_message Status.Information msg
        in
        let add_file contents = add_file data contents filename in
        let res = 
            match Parser.Files.test_file filename with
            | Parser.Files.Is_Poy -> 
                    failwith "TODO Is_poy"
                    (*
                    let data = add_file [] in
                    file_type_message "POY";
                    (match filename with
                    | `Local filename
                    | `Remote filename -> of_file data filename)
                    *)
            | Parser.Files.Is_Clustal| Parser.Files.Is_TinySeq
            | Parser.Files.Is_Fasta| Parser.Files.Is_Genome| Parser.Files.Is_ASN1
            | Parser.Files.Is_Genbank| Parser.Files.Is_INSDSeq| Parser.Files.Is_GBSeq
            | Parser.Files.Is_XML| Parser.Files.Is_NewSeq->
                    let data = add_file [Characters] in
                    file_type_message "input@ sequences";
                    process_molecular_file default_tcm
                                           Cost_matrix.Two_D.default 
                                           Cost_matrix.Three_D.default
                                           annotated Alphabet.nucleotides `DO
                                           is_prealigned `Seq data filename
            | Parser.Files.Is_Phylip->
                    file_type_message "phylip";
                    let fo = Phylip.of_file filename in
                    let fn = FileStream.filename filename in
                    report_static_input fn fo;
                    add_static_parsed_file data fn fo
            | Parser.Files.Is_Hennig-> 
                    let data = add_file [Characters; Trees] in
                    file_type_message "hennig86/Nona";
                    add_static_file `Hennig data filename
            | Parser.Files.Is_Dpread->
                    let data = add_file [Characters; Trees] in
                    file_type_message "dpread file";
                    let parsed = Parser.OldHennig.of_file filename in
                    let fn = FileStream.filename filename in
                    let conv = Parser.OldHennig.to_new_parser fn None parsed in
                    add_static_parsed_file data fn conv
            | Parser.Files.Is_Fixed_States_Dictionary->
                    let data = add_file [] in
                    file_type_message "Fixed@ States@ Dictionary";
                    process_fixed_states data (Some filename)
            | Parser.Files.Is_Dictionary->
                    let data = add_file [Characters] in
                    file_type_message "Synonyms@ Dictionary";
                    add_synonyms_file data filename
            | Parser.Files.Is_Trees->
                    let data = add_file [Trees] in
                    file_type_message "Tree@ List";
                    process_trees data filename
            | Parser.Files.Is_Nexus-> 
                    file_type_message "Nexus@ File";
                    add_static_file `Nexus data filename 
            | Parser.Files.Is_Unknown->
                    let data = add_file [Characters; Trees; CostMatrix] in
                    file_type_message "input@ sequences@ (default)";
                    process_molecular_file default_tcm
                                            Cost_matrix.Two_D.default 
                                            Cost_matrix.Three_D.default
                                            annotated Alphabet.nucleotides
                                            `DO is_prealigned `Seq data filename
            | Parser.Files.Is_ComplexTerminals->
                    let data = add_file [Characters] in
                    file_type_message "Complex@ terminals@ definition@ file";
                    process_complex_terminals data filename
        in
        categorize (remove_taxa_to_ignore res)


let report_kolmogorov_machine file data =
    let fo = Status.Output (file, false, []) in
    let fo = Status.user_message fo in
    let (res, _, _) = Kolmo.Compiler.tree_of_decoder data.machine in
    fo (string_of_int (List.length (Kolmo.S_K.s_encode res)));
    fo "\n%!";
    data
