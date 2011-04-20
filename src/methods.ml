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


let () = SadmanOutput.register "Methods" "$Revision: 2871 $"

exception TimedOut

(** Data *)

let do_job = 11
let process_management = 2
let io = 3
let debugging = 4
let barrier = 5

type cost_modes = [ `Normal | `Normal_plus_Vitamines | `Exhaustive_Weak |
`Exhaustive_Strong | `Iterative of [`ThreeD of int option | `ApproxD of int
option ]  ]
let cost : cost_modes ref = 
    ref `Normal

type filename = [ `Local of string | `Remote of string ]

type support_tree = 
    | Leaf of int
    | Node of float * support_tree * support_tree

type orientation_t = [
| `Orintation
| `NoOrientation
]


type init3D_t = [
| `Init3D
| `NoInit3D
]

type read_option_t = [
| `Init3D of bool
| `Orientation of bool
]

type prealigned_costs = [
    | `Assign_Transformation_Cost_Matrix of (filename * int option) 
    | `Create_Transformation_Cost_Matrix of (int * int) ]

type simple_input = [
    | `Poyfile of filename list
    | `AutoDetect of filename list
    | `Nucleotides of filename list
    | `PartitionedFile of filename list
    | `Aminoacids of filename list
    | `GeneralAlphabetSeq of (filename * filename * read_option_t list)
    | `Breakinv of filename * filename * read_option_t list
    | `Chromosome of filename list
    | `Genome of filename list
    | `ComplexTerminals of filename list
(* sequence file, alphabet file *)
]

type input = [
    | simple_input
    | `Prealigned of (simple_input * prealigned_costs * int)
    | `AnnotatedFiles of simple_input list
]


type information_contained = 
    [ `Nothing | `Cost | `HennigStyle | `Total | `Branches  | `Newick | `Margin of int
    | `NexusStyle | `Collapse of bool ]

type taxon_and_characters = [
    | `Random of float
    | `Names of (bool * string list)
    | `CharSet of (bool * string list)
    | `Missing of bool * int ]

type characters = [
    | `All
    | `Some of (bool * int list)
    | `AllStatic
    | `AllDynamic
    | taxon_and_characters
]

type prep_tail_spec = [
    | `File of filename
    | `Array of int array
]

type level = [
    | `Assign_Level of (int * characters)
]

type transform_cost_matrix = [
    | `Assign_Transformation_Cost_Matrix of ((filename * int option) option * characters)
    | `Create_Transformation_Cost_Matrix of (int * int * characters)
    | `Assign_Affine_Gap_Cost of (int * characters)
    | `Assign_Tail_Cost of (prep_tail_spec * characters)
    | `Assign_Prep_Cost of (prep_tail_spec * characters)
]


type median_solver_chosen = [ `MGR | `SimpleLK | `ChainedLK | `COALESTSP | `BBTSP |
`Siepel | `Albert | `Vinh  ]

type annotate_tool = [ `Mauve of (float*int*float*int) | `Default of (int*int*int) ]

(** parameters used in determining the medians between two chromosomes or genomes *)
type chromosome_pam_t = [
    | `Locus_Inversion of int
    | `Locus_Breakpoint of int
    | `Chrom_Breakpoint of int
    | `Circular of bool
    | `Locus_Indel_Cost of (int * int)
    | `Chrom_Indel_Cost of (int * int)
    | `Chrom_Hom of int
    | `Keep_Median of int
    | `SwapMed of int
    | `Approx of bool 
    | `Symmetric of bool 
    | `Max_3D_Len of int
    | `Max_kept_wag of int
    | `Median_Solver of median_solver_chosen
    | `Annotate_Tool of annotate_tool
]


type indel_prob = 
    [ `Probs  of (float * float) | `Encoding of (string * string) ]

type substitution_prob =
    [ `Probs of float | `Encoding of string ]

type probs =  (indel_prob * substitution_prob)

type kolmo_model = 
    [ `AtomicIndel of ((float option) * ((float * float * float)  option)) 
    | `AffineIndel  of ((float option) * ((float * float * float) option)) ]

type dynamic_char_transform = [
    | `Seq_to_Chrom of (characters * chromosome_pam_t list)
    | `Custom_to_Breakinv of (characters * chromosome_pam_t list)
    | `Annchrom_to_Breakinv of (characters * chromosome_pam_t list)
    | `Change_Dyn_Pam of (characters * chromosome_pam_t list)
    | `Chrom_to_Seq of (characters * chromosome_pam_t list)
    | `Breakinv_to_Custom of (characters * chromosome_pam_t list)
    | `Seq_to_Kolmogorov of (characters * kolmo_model) 
    (*((string * string) option * string
        option * int * int * float))*)
]

    

type terminal_transform = [
    | `RandomizedTerminals
    | `AlphabeticTerminals
]

type ml_substitution = [
    | `JC69
    | `F81
    | `F84 of float list option
    | `HKY85 of float list option
    | `K2P of float list option
    | `TN93 of float list option
    | `F84 of float list option
    | `GTR of float list option
    | `File of string
]


type ml_costfn = [ `MAL     (* maximum average likelihood *)
                 | `MPL     (* most parsimonious likelihood *)
                 | `FLK     (* dynamic alignment with a single matrix *)
                 ] 

type ml_site_variation= [   | `Gamma of int * float option
                            | `Theta of int * (float * float) option ]
type ml_priors = [ `Estimate | `Given of float list | `Equal | `Consistent ]
type ml_gap = [ `Missing | `Independent | `Coupled of float ]
type ml_spec = 
    (characters * ml_costfn * ml_substitution * ml_site_variation option
        * ml_priors * ml_gap)

type char_transform = [
    | dynamic_char_transform
    | `MultiStatic_Aprox of (characters * bool)
    | `Static_Aprox of (characters * bool)
    | `Search_Based of characters
    | `Fixed_States of (characters*(string option))
    | `Partitioned of ([`Clip | `NoClip] * characters)
    | `Direct_Optimization of characters
    | `Automatic_Sequence_Partition of (characters * bool * (int option))
    | `Automatic_Static_Aprox of bool
    | `Prioritize
    | `ReWeight of (characters * float)
    | `WeightFactor of (characters * float)
    | `Prealigned_Transform of characters
    | `UseLikelihood of ml_spec
    | `UseParsimony of characters
    | transform_cost_matrix
    | level 
]

type tree_transform = [ | `EstLikelihood of ml_spec ]

type transform = [
    | char_transform
    | tree_transform
    | terminal_transform
    | `OriginCost of float
]

(* Method employed in a certain cost calculation procedure.
*
* [Exact] calculates the precise cost of each character.
* [Approximate] uses various shortcuts (explained in the code) to estimate the
* cost in a faster way.
* [Static_Aprox] only updates the alignments when a new tree with better cost is
* found.
* [Search_Based] uses a set of possible states as the values in the HTU's. *)
type cost_calculation = [
    | transform
]

type diagnosis = [
    | `AllRootsCost of string option
    | `Implied_Alignment of (string option * characters * bool)]

type summary_class = [ `Individual | `Consensus | `InputFile of string ]

type support_output = [
    | `Bremer of filename list option
    | `Jackknife of summary_class
    | `Bootstrap of summary_class
]

type report = [
    | `KolmoMachine of string option
    | `MstR of string option
    | `Trees of (information_contained list * string option)
    | `TreeCosts of string option
    | `TreesStats of string option
    | `SearchStats of string option
    | `TimeDelta of (string * string option)
    | `Clades of string                 (* file prefix *)
    | `CrossReferences of (characters option * string option)
    | `TerminalsFiles of string option
    | `Supports of (support_output option * string option)
    | `GraphicSupports of (support_output option * string option)
    | `GraphicDiagnosis of string
    | `Dataset of string option
    | `Xslt of (string * string)
    | `Diagnosis of string option
    | `Consensus of (string option * float option)
    | `GraphicConsensus of (string option * float option)
    | `FasWinClad of string option
    | `Nexus of string option
    | `Model of string option
    | `Script of string option * string list
    | `SequenceStats of (string option * characters)
    | `Ci of (string option * characters option)
    | `Ri of (string option * characters option)
    | `CompareSequences of (string option * bool * characters * characters)
    | `ExplainScript of (string * string option)
    | diagnosis
    | `Save of (string * string option)
    | `Load of string
    | `Root of int option
    | `Nodes of string option
    | `RootName of string
]

(* When a maximum number of trees max is to be kept, and the total
* number of trees is higher than that max, the method used to mantain the
* max value.
*
* [Fifo] will simply keep the first trees found.
* [Lifo] will keep the last trees found.
* [Random] will keep a random subset of the trees found. *)
type keep_method = [
    | `First
    | `Last
    | `Keep_Random ]


type store_class = [
    | `Data
    | `Trees
    | `Bremer
    | `Jackknife
    | `Bootstrap ]

type runtime_store = [
    | `Store of (store_class list * string)
    | `Set of (store_class list * string)
    | `Discard of (store_class list * string)
    | `Keep_only of (int * keep_method)
]

type ('c, 'd) character_input_output = [
    | `Characters of 'c
    | `Floats of 'd
]

type ia_seq = [ `DO of int array | `First of int array | `Last of int array ]

(* Note: there are a list of alignments coresponding to a character
 * if the character is a chromosome and broken into diffrent segments *)
type implied_alignment = 
(
    (
        (int * ia_seq array All_sets.IntegerMap.t list) list *
        (int * string * int * [ `Deletion | `Insertion | `Missing] * int Sexpr.t) Sexpr.t list list
    ) * (int * int Sexpr.t) Sexpr.t list list
) list

type ('a, 'b, 'c, 'd) parallel_input = [ 
    | `Trees of 'a 
    | `Data of ('b * int)
    | `DataNTrees of ('b * int * 'a)
    | ('c, 'd) character_input_output
    | `Support of support_tree Sexpr.t
    | `Random_Seed of int
]


type taxa = characters


(** Errors *)
type parallel_special_condition = [ 
    | `Caught of (int * string * exn)   
            (* A non terminating exception that has
            to be informed to the user *)
    | `Uncaught of (int * string * exn) 
            (* A terminating exception that will be 
            * reported for debugging purposes or error causes *)
    | `CleanExit of int
            (* A perfect exit, no need of messaging. If a slave receives this
            * message, it should go down cleanly. *)
    | `Unknown of int
            (* Unknown error, we will simply exit *)
    | `Cleanup (* Deallocate all the memory you can *)
]

type verbosity = Low | Medium | High 

type io = [ 
    | `Status of (string * int)
    | `Error of (string * int)
    | `Information of (string * int)
    | `Output of (string * int)
]

type output_class = [
    | `Information
    | `Error
    | `Output of string option
]

(** Options *)

type tabu_join_strategy = [
    | `UnionBased of int option
    | `AllBased of int option
    | `Partition of [ `Sets of (All_sets.IntSet.t Lazy.t)
                    | `MaxDepth of int 
                    | `ConstraintFile of filename ] list ]

(* defines how branches and the model for likelihood are iterated *)
type tabu_modeli_strategy = [
    | `Threshold of float   (* model iteration based on score improvement *)
    | `MaxCount of int      (* number of search iterations before iterating model *)
    | `Always               (* always iterate the model *)
    | `Both of float * int  (* both threshold and count *)
    | `Neighborhood of float(* same as threshold but weighted around join neighborhood *)
    | `Null ]               (* no iteration of model *)
type tabu_branchi_strategy = [
    | `JoinDelta           (* path along join -> break *)
    | `Neighborhood         (* neighborhood around join point *)
    | `AllBranches          (* iterate all the branches in the tree *)
    | `Null ]               (* no iteration of branches *)
type tabu_iteration_strategy = tabu_modeli_strategy * tabu_branchi_strategy

(* New tree build_method methods.
* 
* [Wagner_Rnd a] creates Wagner (Minimum spanning) trees adding the taxa in
* random order. Note that at most a trees are keep, as different trees could
* have the same cost during the sequential addition.
* [Wagner_Ordered a] creates Wagner (Minimum spanning) trees adding the taxa in
* the order they are found in the input files. Note that at most a trees are
* keep, as different trees could have the same cost during the sequential
* addition.
* [Random a] produces a random trees.
* [Input_file chan] loads the trees stored in the input channel cha. *)
(* So I will modify this set of build methods to wrap them in a nicer manner *)
type build_strategy = 
    int * float * keep_method * cost_calculation list * tabu_join_strategy


type build_method = [
    | `Constraint of (int * float * filename option * cost_calculation list)
    | `Branch_and_Bound  of
        (float option * float option * keep_method * int * cost_calculation list) * tabu_iteration_strategy
    | `Wagner_Rnd of build_strategy
    | `Wagner_Mst of build_strategy
    | `Wagner_Distances of build_strategy
    | `Wagner_Ordered of build_strategy
    | `Build_Random of build_strategy * tabu_iteration_strategy
    | `Nj
    | `Prebuilt of filename 
]
(* because of type constraints the above uses tabu_iteration_strat when in fact
 * the one contained in the Build constructor below will be utilized. *)

type parallelizable_build = build_method

type build = [
    | `Nj
    | `Prebuilt of filename
    | `Build of int * build_method * cost_calculation list * tabu_iteration_strategy
    | `Build_Random of build_strategy * tabu_iteration_strategy
    | `Branch_and_Bound  of
        (float option * float option * keep_method * int * cost_calculation list) * tabu_iteration_strategy
]

(* Optimality criterion employed. Self explanatory. *)
type optimality_criterion = [
    | `Parsimony
    | `Likelihood 
]

(* Neighborhood definition and parameters.
*
* [Spr threshold max method best position] define the neighborhood as all the 
* trees around a
* given tree t produced by the Spr procedure. The method is expected to produce
* a set of trees with cost within the threshold of the local minimum in that
* neighborhood; such set should have at most max trees (unlimited if set to -1),
* kept by using the method defined.
* [Tbr threshold max method best] same as Spr, using Tbr for the neighborhood
* definition instead. 
* [Itself] simply defines the neighborhood of a given tree t as t itself. *)
type neighborhood = [
    | `Spr
    | `Tbr ]

type search_space = [
    | `SingleNeighborhood of neighborhood
    | `ChainNeighborhoods of neighborhood
    | `Alternate of neighborhood * neighborhood
    | `None ]

type tabu_break_strategy = [
    | `Randomized
    | `DistanceSorted of bool
    | `OnlyOnce
]

type tabu_reroot_strategy = [
    | `Bfs of int option
]

type origin_cost = float option

type trajectory_method = [
    | `AllAround of string option
    | `AllThenChoose
    | `BestFirst
    | `PoyDrifting of (float * float)
    | `Annealing of (float * float) ]

type timer = [ `Fixed of float | `Dynamic of (unit -> float) ]

type samples = [
    | `KeepBestTrees 
    | `AllVisited of string option
    | `PrintTrajectory of string option
    | `MaxTreesEvaluated of int
    | `TimeOut of timer 
    | `TimedPrint of (float * string option)
    | `UnionStats of (string option * int)
    | `RootUnionDistr of string option
    | `AttemptsDistr of string option
    | `BreakVsJoin of string option
    | `Likelihood of string option
    | `LikelihoodModel of string option
]

type local_opt = {
    ss : search_space;
    threshold : float;
    num_keep : int;
    keep : keep_method;
    cc : cost_calculation list;
    oo : origin_cost;
    tm : trajectory_method;
    tabu_break   : tabu_break_strategy;
    tabu_join    : tabu_join_strategy;
    tabu_reroot  : tabu_reroot_strategy;
    tabu_iterate : tabu_iteration_strategy;
    samples : samples list;
}

(** [local_optimum] parameters: what to search, threshold, number to keep, keep
    method, parameters for calculating the cost, origin cost for forest search,
    trajectory method *)
type local_optimum = 
       [ `LocalOptimum of local_opt ]

type tree_weights = [
    | `Uniform
]

type fusing_keep_method = [
    | `Best
    | `Better
]

type driven_search = [
    (** [`Fusing (iterations, max_trees, weighting, keep_method, (min, max))] *)
| `Fusing of int option * int option * tree_weights * fusing_keep_method * local_optimum
          * (int * int)
]

(* Method for calculating the support of a given tree. *)
type support_method = [
    | `Bootstrap of int * local_optimum * build * int option 
    | `Jackknife of float * int * local_optimum * build * int option ]

type bremer_support = [
    | `Bremer of (local_optimum * build * int * int)
]


(* Perturbation method. 
*
* [Ratchet p s] is the method for modifying p% of the characters in a given
* tree t, by reweighting their cost with severity s. 
* [No_perturb] does not modify a tree t . *)
type perturb_method = [
    | `Ratchet of float * int 
    | `Resample of [ `Characters of int | `Taxa of int ] 
    | `UnResample of [ `Characters of int | `Taxa of int ] 
    | `UnRatchet 
    | `UnFixImpliedAlignments
    | `FixImpliedAlignments of (characters * bool) ]

type escape_local = [`PerturbateNSearch of (transform list * perturb_method *
local_optimum * int * timer option) ]

(** Method employed in the character weighting. 
*
* [Original_Goloboff a] sets the value of the weighted character parameter as p
* using the original function published by Goloboff, 1994.
* [Product_Goloboff a] uses Original_Goloboff a times the cost of the current
* best tree.
* [Retention index] uses the regular retention index as the weight of each
* character. *)
type character_weight_method = [
    | `Original_Goloboff of int
    | `Product_Goloboff of int
    | `Retention_Index ]
     
(** Method to calculate the HTU states.
*
* [Iterative_Pass a] performs an iterative pass calculation for at most a
* rounds. If should be performed until no changes occur, set a = 0.
* [Fixed_States] calculates the states of each HTU using a list of possible
* sequences.
* [Dynamic_Optimization] makes a full downpass uppass for the statet
* calculation, while performing the sequence alignment. *)
type state_calculation_method = [
    | `Iterative_Pass of int  
    | `Fixed_States of int list
    | `Simple_Dynamic
    | `Goloboff_Approximation
    | `Goloboff_Exact ]

(** Schemes for origin and loss costs *)
type origin_loss_cost = [
| `Flat of int * int ]                 (** Flat cost for each origin and loss *)

(** Methods for complex terminal alignments *)
type complex_term_method = [
| `Strictly_Same  (** Don't allow recombination between elements of a set:  sets
                      are only used for grouping *)
| `Origin_Loss of origin_loss_cost ]   (** Allow simple origin and loss costs *)

(** Method to calculate a consensus tree.
*
* [Majority_Rule a] produces a tree where each branch appears in at least a%
* trees from the set, where 49 < a < 101.
* [Strict_Majority] is equivalent to [Majority_Rule 100]
* [Adams] produces a tree that contains the shared triples in the set of trees.
* *)
type consensus_method = [
    | `Strict_Majority
    | `Majority_Rule of int
    | `Semistrict ]

type consistency_method = 
    | Something;;

(** [make_consis_arr ()] returns an array of the appropriate size to hold the
* consistency values of a node. The size should be the maximum of the
* [consistency_enc] values. *)
let make_consis_arr () =
    Array.make 1 0;;

let consistency_enc = function
    | Something -> 0;;

(** Integer encoding for some calculations in arrays. *)
let support_enc = function
    | `Bremmer _ -> 1
    | `Bootstrap -> 2
    | `Jackknife -> 4
    | `No_support -> raise (Invalid_argument 
        "Encoding for No_support doesn't exist.");;

(** [make_supp_arr ()] is simetric to [make_consis_arr ()] but for the support
* methods. *)
let make_supp_arr () =
    Array.make 4 0;;

let consensus_encoding = function 
    | `Strict_Majority -> 1
    | `Majority_Rule _ -> 2

let state_calculation_enc = function 
    | `Iterative_Pass _ -> 1
    | `Fixed_States _ -> 2
    | `Simple_Dynamic -> 4
    | `Goloboff_Approximation -> 8
    | `Goloboff_Exact -> 16;;

let cost_calculation_enc = function
    | `Exact -> 1
    | `Approximate -> 2
    | `Static_Aprox -> 4
    | `Search_Based -> 8
    | `MultiStatic_Aprox -> 16

let character_weight_enc = function
    | `Original_Goloboff _ -> 1
    | `Product_Goloboff _ -> 2
    | `Retention_Index -> 4;;

let perturb_enc = function 
    | `Ratchet _ -> 1
    | `Char_subset _ -> 2
    | `Taxon_subset _ -> 4
    | `No_perturb -> 8;;

let neighborhood_enc = function
    | `Spr _ -> 1
    | `Tbr _ -> 2
    | `Itself -> 4;;

let keep_enc = function
    | `First -> 1
    | `Last -> 2
    | `Keep_Random -> 4;;

let optimality_criterion_enc = function
    | `Parsimony -> 1
    | `Likelihood -> 2;;

let build_enc = function
    | `Wagner_Rnd _ -> 1
    | `Wagner_Ordered _ -> 2
    | `Build_Random _ -> 4
    | `Input_file  _ -> 8
    | `Wagner_Mst _ -> 16

let default () =
    (`Itself, `Parsimony, `Original_Goloboff 0, `Exact);;

type extract_characters = [
    | `OfNodes of int
]

type 'a character_operations = [
    | `Distance of ('a Sexpr.t * 'a Sexpr.t)
    | `Median of ('a Sexpr.t * 'a Sexpr.t)
]

type clear_item = [
    | `Matrices
    | `SequencePool
]

type 'a plugin_arguments = 
    [ `Empty
    | `Float of float 
    | `Int of int
    | `String of string
    | `Lident of string
    | `Labled of (string * 'a plugin_arguments)
    | `List of 'a plugin_arguments list 
    | `Command of 'a
]

type application = [
    | `Version
    | `Exit
    | `Interactive
    | `ChangeWDir of string
    | `PrintWDir 
    | `Recover
    | `ClearRecovered
    | `Redraw
    | `Help of string option
    | `Wipe
    | `Graph of (string option * bool)
    | `Ascii of (string option * bool)
    | `Memory of string option
    | `KML of (string option * filename * string)
    | `TimerInterval of int
    | `HistorySize of int
    | `Logfile of string option
    | `Normal
    | `Normal_plus_Vitamines
    | `Exhaustive_Weak
    | `Exhaustive_Strong
    | `Iterative of [`ThreeD of int option | `ApproxD of int option ]
    | `ReDiagnose
    | `SetSeed of int
    | `InspectFile of string
    | `ClearMemory of clear_item list
    | `Echo of (string * output_class)
]

type characters_handling = [
    | `RenameCharacters of ((string * string) list)
    | `AnalyzeOnlyCharacters of characters
    | `AnalyzeOnlyCharacterFiles of (bool * filename list)
]

type taxa_handling = [
    | `SynonymsFile of filename list
    | `Synonyms of (string * string) list
    | `AnalyzeOnlyFiles of (bool * filename list)
    | `AnalyzeOnly of taxon_and_characters
]

type tree_handling = [
    | `BestN of int option
    | `BestWithin of float
    | `RandomTrees of int
    | `Unique
]

type script = [ 
    | `Skip
    | `Entry (* Beginning of the program for analyzer *)
    | `StoreTrees
    | `UnionStored
    | `GetStored
    | `ReadScript of string list
    | `Repeat of (int * script list)
    | `OnEachTree of (script list * script list)
    | `ParallelPipeline of (int * script list * script list * script list)
    | `Barrier
    | `GatherTrees of (script list * script list)
    | `GatherJackknife
    | `GatherBootstrap
    | `GatherBremer
    | `SelectYourTrees
    | `StandardSearch of 
        (float option * float option * int option * 
            int option * float option * string option option * string option)
    | `Plugin of (string * script plugin_arguments)
    | input
    | transform
    | build
    | perturb_method 
    | local_optimum
    | support_method
    | bremer_support
    | report
    | runtime_store
    | escape_local
    | application
    | driven_search
    | characters_handling
    | taxa_handling
    | tree_handling
]

type ('a, 'b, 'c, 'd) checkpoints = [
    | `Item of (int * ('a, 'b, 'c, 'd) checkpoints)
    | `ParallelPipeline of (int * script list * script list * script list)
    | ('a Sexpr.t, 'b, 'd Sexpr.t, 'c Sexpr.t) parallel_input
    | `RunForTrees of script list
    | input
    | build
    | build_method
    | local_optimum
    | perturb_method 
    | escape_local
    | extract_characters
    | 'd character_operations
    | support_method
]


(* We are missing some methods to operate directly over characters, for example,
* to be able to calculate things not associated with a tree but rather directly
* over some set of observations (for example a large matrix of distances between
* sequences. *)


