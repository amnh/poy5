(* POY 5.0 Beta. A phylogenetic analysis program using Dynamic Homologies.    *)
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

exception Exit 

let () = SadmanOutput.register "PoyCommand" "$Revision: 3383 $"

let debug = false 

type read_option_t = [
    | `Init3D of bool
    | `Orientation of bool
    | `InputFile of Methods.filename
    | `CostMatrix of Methods.filename
    | `Level of int * Methods.keep_method
    | `Tie_Breaker of Methods.keep_method
    | `Affine of int
]

type otherfiles = [
    | `AutoDetect of string
    | `Nucleotides of string list
    | `PartitionedFile of string list
    | `Aminoacids of (string list * read_option_t list)
    | `GeneralAlphabetSeq of (string list * string * read_option_t list) 
    | `Breakinv of (string list * string * read_option_t list)
    | `Chromosome of string list
    | `Prealigned of (otherfiles * Methods.prealigned_costs * int)
    | `Genome of string list
    | `ComplexTerminals of string list
]

type reada = Methods.input

type keepa = Methods.tree_handling

type old_identifiers = [
    | `All
    | `Names of (bool * string list)
    | `Some of (bool * int list)
    | `CharSet of (bool * string list)
    | `AllDynamic
    | `AllStatic
    | `Missing of bool * int
    | `Random of float
    | `Range of (bool * string * int * int)
]

type identifiers = [
    | old_identifiers
    | `Files of (bool * string list)
]

(*chromosome args is a list of dynamic pam args, not just for chromosome.
* we used to use transform(dyn_pam:(...)), everything for dyamic pam is here, we
* need to split them into their own characters,like genome/chrom/dna/etc *)
type chromosome_args = [
    | `Median_Solver of Methods.median_solver_chosen
    | `Annotate_Tool of Methods.annotate_tool (*annotated tool = mauve or default*)
    | `Locus_Inversion of int (** the cost of a locus inversion operation inside a chromosome *)
    | `Locus_Breakpoint of int (* the cost of a locus breakpoint operation inside a chromosome *)
    | `Circular of bool (** indicate if the chromosome is circular or not *)

    (** [(a, b)] is indel cost of a locus in a chromosome
    * where [a] is opening cost, [b] is extension cost *)
    | `Locus_Indel_Cost of (int * int) 

    (** [(a, b)] is indel cost of a chromosome in a genome 
    * where [a] is opening cost, [b] is extension cost *)
    | `Chrom_Indel_Cost of (int * int)

    (** The maximum cost between two chromosomes
    * at which they are considered as homologous chromosomes *)
    | `Chrom_Hom of int 
    
    (** The cost of a breakpoint happing between two chromosome *)
    | `Translocation of int (* Breakpoint cost between loci of two different chromosomes *)

    (** The maximum number of medians at one node kept during the search*)
    | `Keep_Median of int 

(** number iterations are applied in refining alignments with rearrangements *)
    | `SwapMed of int 

(** approx = true, the median sequence of X and Y is approximated by either X or Y,
* otherwise, calculate as a set of median between X and Y *)
    | `Approx of bool

    (** symmetric = true, calculate the both distances between X to Y
     * and vice versa. Otherwise, only form X to Y *) 
    | `Symmetric of bool 

    (** maximum length of sequences aligned by 3D-alignment *)
    | `Max_3D_Len of int 

    | `Max_kept_wag of int
]

type polymorphism_arg = Methods.polymorphism_arg


type transform_method = [
    | `RandomizedTerminals
    | `AlphabeticTerminals
    | `Prealigned_Transform
    | `EstLikelihood of 
        ( Methods.ml_alphabet * Methods.ml_costfn * Methods.ml_model * 
          Methods.ml_site_variation option * Methods.ml_priors * Methods.ml_gap)
    | `UseParsimony
    | `UseLikelihood of 
        ( Methods.ml_alphabet * Methods.ml_costfn * Methods.ml_model * 
          Methods.ml_site_variation option * Methods.ml_priors * Methods.ml_gap)
    | `Level of (int * Methods.keep_method)
    | `Tcm of (string * (int * Methods.keep_method) option)
    | `Gap of (int * int)
    | `AffGap of int
    | `TailInput of int array
    | `TailFile of string
    | `PrepInput of int array
    | `PrepFile of string
    | `StaticApproximation of bool
    | `MultiStaticApproximation of bool
    | `Automatic_Static_Aprox of bool
    | `ReWeight of float
    | `WeightFactor of float
    | `Automatic_Sequence_Partition of (bool * int option)
    | `Prioritize
    | `SearchBased of (string * string) list
    | `Fixed_States of (string option * polymorphism_arg option) 
    | `Partitioned of [`Clip | `NoClip]
    | `Direct_Optimization
    | `SeqToChrom of chromosome_args list
    | `ChangeDynPam of chromosome_args list
    | `CustomToBreakinv of chromosome_args list
    | `AnnchromToBreakinv of chromosome_args list
    | `ChromToSeq of chromosome_args list
    | `BreakinvToSeq of chromosome_args list
    | `OriginCost of float
(*    | `Seq_to_Kolmogorov of Methods.kolmo_model*)
]

type transforma = (identifiers * transform_method)

type transform = [
    | `Transform of transforma list
]

type cost_calculation = [
    | `Exhaustive_Weak
    | `Exhaustive_Strong
    | `Iterative of [ `ThreeD of int option | `ApproxD of int option ] 
    | `Normal_plus_Vitamines
    | `Normal
]

type alignment_mode = [
    | `Algn_Newkk
    | `Algn_Normal
]

type keep_method = [
    | `Last
    | `First
    | `Keep_Random
]

type thresh_trees = [
    | `Threshold of float
    | `Trees of int
]

type iteration_strategy = [
    | `NullModel
    | `AlwaysModel
    | `ThresholdModel of float
    | `MaxCountModel of int
    | `BothModel of float * int
    | `NullBranches
    | `AllBranches
    | `JoinDeltaBranches
    | `NeighborhoodBranches of int
]

type builda = [
    | thresh_trees
    | `Lookahead of int
    | `Nj
    | `Prebuilt of Methods.filename
    | `Mst
    | `DistancesRnd
    | `Branch_and_Bound of float option
    | `Constraint of string option
    | `Random
    | `RandomTree
    | `Ordered
    | keep_method
    | transform
    | Methods.tabu_join_strategy
    | `IterationB of iteration_strategy list
]

type swap_neighborhood = [
    | `Spr
    | `Tbr
]

type swap_strategy = [
    | `SingleNeighborhood of swap_neighborhood
    | `ChainNeighborhoods of swap_neighborhood
    | `Alternate of (swap_neighborhood * swap_neighborhood)
    | `None
]

type swap_trajectory = [
    | `AllAround of string option
    | `AllThenChoose
    | `BestFirst
    | `PoyDrifting of (float * float)
    | `Annealing of (float * float)
]

type swapa = [
    | `IterationS of iteration_strategy list
    | `Forest of float
    | thresh_trees
    | keep_method
    | swap_strategy
    | transform
    | swap_trajectory
    | Methods.tabu_break_strategy
    | Methods.tabu_join_strategy
    | Methods.tabu_reroot_strategy
    | Methods.samples
]

type swap = [
    | `Swap of swapa list
]

type build = [ `Build of builda list ]

type supporta = [
    | build
    | swap
    | `Bremer
    | `Jackknife of [ `Select of float | `Resample of int ] list
    | `Bootstrap of int option
]

type poy_file_commands = [
    | `Load of string
    | `Save of (string * string option)
    | `InspectFile of string
]

type internal_memory = [
    | `Store of (Methods.store_class list * string)
    | `Use of (Methods.store_class list * string)
    | `Discard of (Methods.store_class list * string)
]

type settings = [
    | `TimerInterval of int
    | `Parmap of int
    | `HistorySize of int
    | `Logfile of string option
    | cost_calculation
    | `SetSeed of int
    | `Root of int option
    | `RootName of string
    | `Alias of string * [ `Codon of old_identifiers | `Chars of old_identifiers ]
    | alignment_mode
    | `Optimization of Numerical.opt_modes
]

type output_class = [
    | `Information
    | `Error
    | `Output of string option
]

type application = [
    | `Version
    | `ChangeWDir of string
    | `PrintWDir
    | `Exit
    | `Recover
    | `ClearRecovered
    | `Echo of (string * output_class list)
    | `Help of string option
    | `Set of settings list
    | `Redraw
    | `Wipe 
    | `ReDiagnose
    | `ReDiagnoseTrees
    | `ClearMemory of Methods.clear_item list
    | `ReadScript of string list
    | poy_file_commands
    | internal_memory
]

type charortax = [
    | `Characters
    | `Taxa
]

type charoper = [
    | `Distance
    | `Median
    | `Taxa of identifiers
    | `Characters of identifiers
]

type reporta = [
    | `File of string
    | `Data
    | `Xslt of (string * string)
    | `KML of (string * string)
    | `Ascii of Methods.report_branch option
    | `Memory
    | `Graph of Methods.report_branch option
    | `Trees of Methods.information_contained list
    | `MstR
    | `TreeCosts
    | `KolmoMachine
    | `TreesStats
    | `SearchStats 
    | `TimeDelta of string
    | `SequenceStats of old_identifiers
    | `Ci of old_identifiers option
    | `Ri of old_identifiers option
    | `CompareSequences of (bool * old_identifiers * old_identifiers)
    | `FasWinClad
    | `Nexus
    | `Model of old_identifiers
    | `LKSites of old_identifiers
    | `DebugData
    | `Pairwise of old_identifiers
    | `Script of string list
    | `ExplainScript of string
    | `Consensus of float option
    | `GraphicConsensus of float option
    | `Clades
    | `CrossReferences of old_identifiers option
    | `RobinsonFoulds
    | `TerminalsFiles
    | `Supports of Methods.support_output option
    | `GraphicSupports of Methods.support_output option
    | `AllRootsCost
    | `Implied_Alignments of identifiers * bool
    | `Topo_Selection of Methods.ml_topo_test list
    | `GraphicDiagnosis of Methods.diagnosis_report_type
    | `Diagnosis of Methods.diagnosis_report_type 
    | `Nodes
]

type perturba = [
    | `Ratchet of (float * int)
    | `Resample of int
    | `Iterations of int
    | `TimeOut of Methods.timer
    | swap
    | transform
]

type selecta = [
    | charortax
    | identifiers
    | Methods.tree_handling
]

type renamea = [
    | charortax
    | `File of string
    | `Syn of (string * string)
]

type fusea = [
    | `Keep of int
    | `Iterations of int
    | `Replace of [`Better | `Best]
    | `Swap of swapa list
    | `Weighting of [`Uniform]
    | `Clades of int * int option
    | `IterationF of iteration_strategy list
]

type searcha = [
    | `Build of bool
    | `Transform of bool
]

type std_searcha = [
    | `MaxRam of int
    | `MinHits of int
    | `MaxTime of float
    | `MinTime of float
    | `Target of float
    | `Visited of string option
    | `ConstraintFile of string 
]

type command = [
    | `Read of reada list 
    | build
    | swap
    | `Fuse of fusea list
    | `Support of supporta list
    | `Calculate of charoper list
    | `Report of reporta list
    | `Plugin of (string * command Methods.plugin_arguments)
    | `Select of selecta list
    | `Rename of renamea list
    | `Search of searcha list
    | `StandardSearch of std_searcha list
    | transform
    | `Perturb of perturba list
    | `Repeat of (int * command list)
    | application
]

    let all_store_types = [ `Bremer; `Jackknife; `Bootstrap; `Data; `Trees ]
    (* Transform *)
    let transform_transform acc (id, x) = match id with
        | `Files _ ->
            let msg = "In@ transformations@ of@ characters,@ specifying@ a@ " ^
                      "character@ set@ using@ a@ file@ is@ not@ allowed.@ I@ will@ " ^
                      "ignore@ this@ command." in
            Status.user_message Status.Error msg;
            acc
        | #Methods.characters as id ->
            begin match x with
                | `RandomizedTerminals -> `RandomizedTerminals :: acc
                | `AlphabeticTerminals -> `AlphabeticTerminals :: acc
                | `Prealigned_Transform -> (`Prealigned_Transform id) :: acc
                | `EstLikelihood (a,b,c,d,e,f) -> (`EstLikelihood (id,a,b,c,d,e,f))::acc
                | `UseParsimony -> (`UseParsimony id) :: acc
                | `UseLikelihood (a,b,c,d,e,f) -> (`UseLikelihood (id,a,b,c,d,e,f))::acc
                | `Level (l,tb) -> (`Assign_Level (l,tb,id))::acc
                | `Tcm (f,l_and_tb) -> (`Assign_Transformation_Cost_Matrix ((Some ((`Local f),l_and_tb)), id)) :: acc
                | `SearchBased f -> (`Search_Based (f,id)) :: acc
                | `Gap (a, b) -> (`Create_Transformation_Cost_Matrix (a, b, id)) :: acc
                | `AffGap c -> (`Assign_Affine_Gap_Cost (c, id)) :: acc
                | `PrepInput x -> (`Assign_Prep_Cost ((`Array x), id)) :: acc
                | `PrepFile x -> (`Assign_Prep_Cost ((`File (`Local x)), id)) :: acc
                | `TailInput x -> (`Assign_Tail_Cost ((`Array x), id)) :: acc
                | `TailFile x -> (`Assign_Tail_Cost ((`File (`Local x)), id)) :: acc
                | `MultiStaticApproximation noninf -> (`MultiStatic_Aprox (id, noninf)) :: acc
                | `StaticApproximation noninf -> (`Static_Aprox (id, noninf)) :: acc
                | `Automatic_Static_Aprox sens -> (`Automatic_Static_Aprox sens) :: acc
                | `ReWeight w -> (`ReWeight (id, w)) :: acc
                | `WeightFactor w -> (`WeightFactor (id, w)) :: acc
                | `Automatic_Sequence_Partition (sens, x) -> (`Automatic_Sequence_Partition (id, sens, x)) :: acc
                | `Prioritize -> `Prioritize :: acc
                | `Fixed_States (filename,polymph) -> (`Fixed_States (id,filename,polymph)) :: acc
                | `Partitioned mode -> (`Partitioned (mode, id)) :: acc
                | `Direct_Optimization -> (`Direct_Optimization id) :: acc
                | `SeqToChrom x -> (`Seq_to_Chrom (id, x)) :: acc
                | `CustomToBreakinv x -> (`Custom_to_Breakinv (id, x)) :: acc
                | `AnnchromToBreakinv x -> (`Annchrom_to_Breakinv (id, x)) :: acc
    (*            | `Seq_to_Kolmogorov x -> (`Seq_to_Kolmogorov (id, x)) :: acc*)
                | `ChangeDynPam x -> (`Change_Dyn_Pam (id, x)) :: acc
                | `ChromToSeq x -> (`Chrom_to_Seq (id, x)) :: acc
                | `BreakinvToSeq x -> (`Breakinv_to_Custom (id, x)) :: acc
                | (`OriginCost _) as id -> id :: acc
            end

    let transform_transform_arguments x =
        List.fold_left transform_transform [] x

    (* Reading files *)
    let modify_acc acc c = function
        | [] -> acc
        | files -> (`Other (files, c)) :: acc

    (*let console_script = ref [] *)
    (*let add_command_to_console_script str : char Stream.t =*)
    (*    console_script := str :: !console_script; str*)

    let iter_default =  `Always, `AllBranches

    (* Building *)
    let build_default_method_args = (1, 0.0, `Last, [], `UnionBased None)
    let build_default_method = `Wagner_Rnd build_default_method_args
    let build_default = (10, build_default_method,[],iter_default)

    let transform_iterations 
        (items : iteration_strategy list) :
            Methods.tabu_modeli_strategy * Methods.tabu_branchi_strategy = 
        List.fold_left
            (fun (a,b) -> function
                | `NullModel        -> (`Null,b)
                | `AlwaysModel      -> (`Always,b)
                | `ThresholdModel x -> (`Threshold x, b)
                | `MaxCountModel x  -> (`MaxCount x, b)
                | `BothModel (x,y)  -> (`Both (x,y), b)
                | `NullBranches     -> (a,`Null)
                | `AllBranches      -> (a,`AllBranches)
                | `JoinDeltaBranches -> (a,`JoinDelta)
                | `NeighborhoodBranches y -> (a,`Neighborhood y))
            iter_default items

    let transform_build
            ((n, (meth:Methods.build_method), 
             (trans:Methods.transform list),
             (iter:Methods.tabu_iteration_strategy)) as acc)
            (builder: builda) = match builder with
        | `Nj -> (n, `Nj, trans,iter)
        | `Prebuilt fn -> (n, (`Prebuilt fn), trans,iter)
        | `RandomTree ->
            begin match meth with
                | `Nj
                | `Prebuilt _ -> acc
                | `Wagner_Ordered x 
                | `Wagner_Distances x 
                | `Wagner_Mst x
                | `Build_Random (x,_)
                | `Wagner_Rnd x -> (n, (`Build_Random (x,iter)), trans, iter)
                | `Constraint _ ->
                    failwith "Constraint tree has already been selected as build method."
                | `Branch_and_Bound _ ->
                    failwith "Branch and bound tree has already been selected as build method."
            end
        | `Constraint file ->
            let file = match file with
                | None -> None
                | Some x -> Some (`Local x)
            in
            begin match meth with
                | `Wagner_Ordered (keep_max, _, keep_method, lst, _) 
                | `Wagner_Distances (keep_max, _, keep_method, lst, _) 
                | `Wagner_Mst (keep_max, _, keep_method, lst, _)
                | `Build_Random ((keep_max, _, keep_method, lst, _),_)
                | `Wagner_Rnd (keep_max, _, keep_method, lst, _) -> 
                        (n, (`Constraint (1, 0.0, file, lst)), trans, iter)
                | `Branch_and_Bound _ ->
                    failwith
                        "Branch and bound has already been selected as build method."
                | `Constraint _ ->
                    failwith "Constraint has already been selected as build method."
                | `Nj -> 
                    failwith "Neighbor joining has already been selected as build method."
                | `Prebuilt _ -> 
                    failwith "Prebuilt has already been selected as build method."
            end
        | `Branch_and_Bound bound ->
            begin match meth with
                | `Wagner_Ordered (keep_max, _, keep_method, lst, _) 
                | `Wagner_Distances (keep_max, _, keep_method, lst, _) 
                | `Wagner_Mst (keep_max, _, keep_method, lst, _)
                | `Build_Random ((keep_max, _, keep_method, lst, _),_)
                | `Wagner_Rnd (keep_max, _, keep_method, lst, _) -> 
                    (n,`Branch_and_Bound ((bound, None, keep_method, keep_max, lst),iter),trans,iter)
                | `Branch_and_Bound (x,_) ->
                    (n, `Branch_and_Bound (x,iter), trans, iter)
                | `Constraint _ ->
                    failwith "Constraint has already been selected as build method."
                | `Nj -> 
                    failwith "Neighbor joining has already been selected as build method."
                | `Prebuilt _ -> 
                    failwith "Prebuilt has already been selected as build method."
            end
        | `DistancesRnd ->
            begin match meth with
                | `Prebuilt _
                | `Wagner_Distances _ -> acc
                | `Wagner_Mst x
                | `Wagner_Rnd x 
                | `Wagner_Ordered x
                | `Build_Random (x,_) -> (n, (`Wagner_Distances x), trans, iter)
                | `Constraint _ ->
                    failwith "Constraint has already been selected as build method."
                | `Branch_and_Bound _ -> 
                    failwith "Branch and bound tree has already been selected as build method."
                | `Nj -> 
                    failwith "Neighbor joining has already been selected as build method."
            end
        | `Mst ->
            begin match meth with
                | `Prebuilt _
                | `Wagner_Mst _ -> acc
                | `Wagner_Distances x
                | `Wagner_Rnd x
                | `Wagner_Ordered x
                | `Build_Random (x,_) -> (n, (`Wagner_Mst x), trans, iter)
                | `Constraint _ ->
                    failwith "Constraint has already been selected as build method."
                | `Branch_and_Bound _ -> 
                    failwith "Branch and bound tree has already been selected as build method."
                | `Nj -> 
                    failwith "Neighbor joining has already been selected as build method."
            end
        | `Random ->
            begin match meth with
                | `Prebuilt _
                | `Wagner_Rnd _ -> acc
                | `Wagner_Distances x
                | `Wagner_Mst x
                | `Wagner_Ordered x
                | `Build_Random (x,_) -> (n, (`Wagner_Rnd x), trans, iter)
                | `Constraint _ ->
                    failwith "Constraint has already been selected as build method."
                | `Branch_and_Bound _ -> 
                    failwith "Branch and bound tree has already been selected as build method."
                | `Nj -> 
                    failwith "Neighbor joining has already been selected as build method."
            end
        | `Ordered ->
            begin match meth with
                | `Prebuilt _
                | `Wagner_Ordered _ -> acc
                | `Wagner_Distances x
                | `Wagner_Mst x
                | `Build_Random (x ,_)
                | `Wagner_Rnd x -> (n, (`Wagner_Ordered x), trans, iter)
                | `Constraint _ ->
                    failwith "Constraint has already been selected as build method."
                | `Branch_and_Bound _ -> 
                    failwith "Branch and bound tree has already been selected as build method."
                | `Nj -> 
                    failwith "Neighbor joining has already been selected as build method."
            end
        | `Threshold x ->
            let converter (a, _, c, d, e) = (a, x, c, d, e) in
            let nmeth = match meth with
                | `Branch_and_Bound _
                | `Constraint _
                | `Build_Random _
                | `Nj
                | `Prebuilt _ -> meth
                | `Wagner_Distances y -> `Wagner_Distances (converter y)
                | `Wagner_Mst y -> `Wagner_Mst (converter y)
                | `Wagner_Rnd y -> `Wagner_Rnd (converter y)
                | `Wagner_Ordered y -> `Wagner_Ordered (converter y)
            in
            n, nmeth, trans, iter
        | `Trees x ->
            (x, meth, trans, iter)
        | `Lookahead x ->
            let converter (_, b, c, d, e) = (x, b, c, d, e) in
            let nmeth = match meth with
                | `Branch_and_Bound _
                | `Constraint _
                | `Build_Random _
                | `Nj
                | `Prebuilt _ -> meth
                | `Wagner_Distances y -> `Wagner_Distances (converter y)
                | `Wagner_Mst y -> `Wagner_Mst (converter y)
                | `Wagner_Rnd y -> `Wagner_Rnd (converter y)
                | `Wagner_Ordered y -> `Wagner_Ordered (converter y)
            in
            n, nmeth, trans, iter
        | `Last
        | `First
        | `Keep_Random as x ->
            let converter (a, b, _, c, d) = (a, b, x, c, d) in
            let nmeth = match meth with
                | `Constraint _
                | `Nj
                | `Prebuilt _ -> meth
                | `Wagner_Distances y -> `Wagner_Distances (converter y)
                | `Wagner_Mst y -> `Wagner_Mst (converter y)
                | `Wagner_Rnd y -> `Wagner_Rnd (converter y)
                | `Wagner_Ordered y -> `Wagner_Ordered (converter y)
                | `Build_Random (y,i) -> `Build_Random ((converter y),i)
                | `Branch_and_Bound ((a,b, _,c,d),i) ->
                        `Branch_and_Bound ((a,b,x,c,d),i)
            in
            n, nmeth, trans, iter
        | `Transform x ->
            let t = transform_transform_arguments x in
            (n, meth, (t @ trans), iter)
        | `IterationB xs ->
            let (m,b) as i = transform_iterations xs in
            n, meth, trans, i
        | #Methods.tabu_join_strategy as tabu ->
            let converter (a, b, c, d, _) = (a, b, c, d, tabu) in
            let nmeth = match meth with
                | `Constraint _
                | `Branch_and_Bound _
                | `Nj
                | `Prebuilt _ -> meth
                | `Wagner_Distances y -> `Wagner_Distances (converter y)
                | `Wagner_Mst y -> `Wagner_Mst (converter y)
                | `Wagner_Rnd y -> `Wagner_Rnd (converter y)
                | `Wagner_Ordered y -> `Wagner_Ordered (converter y)
                | `Build_Random (y,i) -> `Build_Random ((converter y),i)

            in
            n, nmeth, trans, iter

    let transform_build_arguments x
        : int * Methods.build_method * Methods.cost_calculation list * Methods.tabu_iteration_strategy =
        let (x, y, z, i) = List.fold_left transform_build build_default x in
        (x, y, List.rev z, i)

    (* Swapping *)
    let swap_default ={ Methods.ss =  `Alternate (`Spr, `Tbr);
                        Methods.threshold = 0.0;
                        Methods.num_keep = 1;
                        Methods.keep = `Last;
                        Methods.cc = [];
                        Methods.oo = None;
                        Methods.tm = `BestFirst;
                        Methods.tabu_break = `DistanceSorted false;
                        Methods.tabu_join = `UnionBased None;
                        Methods.tabu_reroot = `Bfs None;
                        Methods.tabu_iterate = iter_default;
                        Methods.samples = []; }

    let swap_default_none = { swap_default with Methods.ss = `None }

    let transform_swap l_opt (param : swapa) = match param with
        | `Threshold thres ->
            { l_opt with Methods.threshold=thres }
        | `Trees keep ->
            { l_opt with Methods.num_keep = keep }
        | #Methods.keep_method as keepm ->
            { l_opt with Methods.keep = keepm }
        | #swap_strategy as space ->
            { l_opt with Methods.ss = space }
        | `Transform x ->
            let t = transform_transform_arguments x in 
            let cclist = t @ l_opt.Methods.cc in
            { l_opt with Methods.cc = cclist }
        | `IterationS xs ->
            let i = transform_iterations xs in
            { l_opt with Methods.tabu_iterate = i; }
        | `Forest cost ->
            let origin = Some cost in
            print_endline ("Forest: "^string_of_float cost);
            { l_opt with Methods.oo = origin }
        | #swap_trajectory as traj ->
            { l_opt with Methods.tm = traj }
        | #Methods.tabu_break_strategy as tabu ->
            { l_opt with Methods.tabu_break = tabu }
        | #Methods.tabu_join_strategy as tabu ->
            { l_opt with Methods.tabu_join = tabu }
        | #Methods.tabu_reroot_strategy as tabu ->
            { l_opt with Methods.tabu_reroot = tabu }
        | #Methods.samples as s ->
            { l_opt with Methods.samples = s :: l_opt.Methods.samples }


    let transform_swap_arguments (lst : swapa list) =
        let options = List.fold_left transform_swap swap_default lst in
        `LocalOptimum { options with Methods.cc = List.rev options.Methods.cc } 

    (* Fusing *)
    let rec transform_fuse ?(iterations=None) ?(keep=None) ?(replace=`Better)
            ?(search=`LocalOptimum swap_default_none) ?(weighting=`Uniform) ?(clades=(4,7)) =
        function
        | [] -> `Fusing (iterations, keep, weighting, replace, search, clades)
        | x :: xs ->
            begin match x with
               | `Keep keep ->
                     let keep = Some keep in
                     transform_fuse ~iterations ~keep ~replace ~search ~weighting ~clades xs
               | `Iterations iterations ->
                     let iterations = Some iterations in
                     transform_fuse ~iterations ~keep ~replace ~search ~weighting ~clades xs
               | `Replace replace ->
                    transform_fuse ~iterations ~keep ~replace ~search ~weighting ~clades xs
               | `Weighting weighting ->
                    transform_fuse ~iterations ~keep ~replace ~search ~weighting ~clades xs
               | `Clades (i, Some j) ->
                     let clades = (i, j) in
                     transform_fuse ~iterations ~keep ~replace ~search ~weighting ~clades xs
               | `Clades (i, None) ->
                     let clades = (i, i) in
                     transform_fuse ~iterations ~keep ~replace ~search ~weighting ~clades xs
               | `Swap swapas ->
                     let search = transform_swap_arguments swapas in
                     transform_fuse ~iterations ~keep ~replace ~search ~weighting ~clades xs
               | `IterationF stuff -> 
                    let strat = transform_iterations stuff in
                    let search = match search with
                        | `LocalOptimum x -> `LocalOptimum {x with Methods.tabu_iterate = strat; }
                    in
                    transform_fuse ~iterations ~keep ~replace ~search ~weighting ~clades xs
            end


    (* Perturbing *)
    let pert_def_swap : Methods.local_optimum = `LocalOptimum swap_default
    let pert_def_iter    = 1
    let pert_def_perturb = `Ratchet (0.25, 2)
    let pert_def_trans   = []
    let pert_def_time    = None
    let perturb_default  =
        (pert_def_trans,pert_def_perturb,pert_def_swap,pert_def_iter,pert_def_time)

    let transform_perturb (tr, m, sw, it, timeout) = function
        | `TimeOut x -> (tr, m, sw, it, Some x)
        | `Transform x  ->
            let x = transform_transform_arguments x in
            ((x @ tr), m, sw, it, timeout)
        | `Swap x -> 
            let x = transform_swap_arguments x in
            tr, m, x, it, timeout
        | `Ratchet _ as x -> (tr, x, sw, it, timeout)
        | `Resample x -> tr, `Resample x, sw, it, timeout
        | `Iterations it -> tr, m, sw, it, timeout

    let transform_perturb_arguments x : Methods.script list = 
        let tr, a, b, c, d = List.fold_left transform_perturb perturb_default x in
        `PerturbateNSearch (List.rev tr, a, b, c, d) :: []

    (* Support *)
    let support_default_swap = `LocalOptimum swap_default
    let support_select = 36.0
    let support_resamplings = 5
    let support_default_build = (1, build_default_method, [],iter_default)
    let support_default = 
        `Bremer, (support_select, support_resamplings, support_default_swap, support_default_build)

    let transform_support (meth, (ss, sr, ssw, sb)) = function
        | `Swap s ->
            let v = transform_swap_arguments s in
            (meth, (ss, sr, v, sb))
        | `Build s ->
            let z = transform_build_arguments s in
            (meth, (ss, sr, ssw, z))
        | `Bremer as x -> (x, (ss, sr, ssw, sb))
        | `Bootstrap x ->
            let nx = match x with
                | None -> sr
                | Some v -> v
            in
            (`Bootstrap, (ss, nx, ssw, sb))
        | `Jackknife items ->
            let process_item (a, b, c, sb) = function
                | `Select x   -> (x, b, c, sb)
                | `Resample x -> (a, x, c, sb)
            in
            let r = List.fold_left process_item (ss, sr, ssw, sb) items in
            (`Jackknife, r)

    let rec transform_support_arguments args =
        (* lift iteration to the local_optimum; for ease of use *)
        let lift (`LocalOptimum l_opt) ((_,_,_,i) as all) = 
            `LocalOptimum {l_opt with Methods.tabu_iterate = i},all
        in
        match List.fold_left transform_support support_default args with
        | `Bremer, (_, _, c, d) ->
            let c,d = lift c d in `Bremer (c, (`Build d), 0, 1)
        | `Jackknife, (a, b, c, d) ->
            let c,d = lift c d in `Jackknife (a, b, c, (`Build d), None)
        | `Bootstrap, (_, b, c, d) ->
            let c,d = lift c d in `Bootstrap (b, c, (`Build d), None)
        
    (* Reporting things *)
    let transform_report ((acc : Methods.script list), file) (item : reporta) = 
        match item with
        | `File x -> (acc, Some x)
        | `Data -> (`Dataset file) :: acc, file
        | `Xslt x -> (`Xslt x) :: acc, file
        | `Ascii x -> (`Ascii (file, x)) :: acc, file
        | `KML (plugin, out_file) -> 
            begin match file with
                | None -> failwith "KML generation requires an output file"
                | Some f -> (`KML (Some plugin, `Local out_file, f)) :: acc, file
            end
        | `Memory -> (`Memory file) :: acc, file
        | `Graph x -> 
            begin match acc, file with
                | (_ :: _), None -> ()
                | (_ :: _), Some _ ->
                    let msg = "You@ have@ requested@ to@ output@ a@ postscript@ "^
                        "file@ in@ a@ flat@ text@ file,@ this@ is@ probably@ not@ "^
                        "what@ you@ expect.@ If@ you@ want@ to@ output@ a@ "^
                        "graphical@ version@ of@ your@ tree@ in@ a@ file,@ "^
                        "assign@ a@ different@ filename@ to@ the@ graph@ command@ "^
                        "in@ the@ report@ line."
                    in
                    Status.user_message Status.Warning msg                       
                | _ -> ()
            end;
            (`Graph (file, x)) :: acc, file
        | `Trees lst         -> (`Trees (lst, file)) :: acc, file
        | `MstR              -> (`MstR file) :: acc, file
        | `KolmoMachine      -> (`KolmoMachine (file)) :: acc, file
        | `TreeCosts         -> (`TreeCosts (file)) :: acc, file
        | `TreesStats        -> (`TreesStats (file)) :: acc, file
        | `SearchStats       -> (`SearchStats (file)) :: acc, file
        | `TimeDelta str     -> (`TimeDelta (str, file)) :: acc, file
        | `Consensus v       -> (`Consensus (file, v)) :: acc, file
        | `GraphicConsensus v-> (`GraphicConsensus (file, v)) :: acc, file
        | `SequenceStats c   -> (`SequenceStats (file, c)) :: acc, file
        | `Ci c              -> (`Ci (file, c)) :: acc, file
        | `Ri c              -> (`Ri (file, c)) :: acc, file
        | `CompareSequences (a,b,c)
                             -> (`CompareSequences (file, a, b, c)) :: acc, file
        | `FasWinClad        -> (`FasWinClad (file)) :: acc, file
        | `Nexus             -> (`Nexus (file)) :: acc, file
        | `Pairwise x        -> (`Pairwise (file,x)) :: acc, file
        | `LKSites x         -> (`LKSites (file,x)) :: acc, file
        | `DebugData         -> `DebugData :: acc, file
        | `Topo_Selection x  ->
            (`Topo_Selection (file,x)) :: acc, file
        | `Model x           -> (`Model (file,x)) :: acc, file
        | `Script lst        -> (`Script (file,lst)) :: acc, file
        | `ExplainScript scr -> (`ExplainScript (scr, file)) :: acc, file
        | `Clades ->
            begin match file with
                | None ->
                    let msg = "Sorry,@ I@ need@ a@ filename@ to@ output@ the@ " ^
                        "clades@ file.@ I@ will@ ignore@ your@ clades@ request."
                    in
                    Status.user_message Status.Error msg;
                    acc, file
                | Some f -> (`Clades f) :: acc, file
            end
        | `CrossReferences x -> (`CrossReferences (x, file)) :: acc, file
        | `RobinsonFoulds -> (`RobinsonFoulds file) :: acc, file
        | `TerminalsFiles -> (`TerminalsFiles file) :: acc, file
        | `Supports c -> (`Supports (c, file)) :: acc, file
        | `GraphicSupports c -> (`GraphicSupports (c, file)) :: acc, file
        | `AllRootsCost -> (`AllRootsCost file) :: acc, file
        | `Implied_Alignments (id, include_header) ->
            begin match id with
                | #Methods.characters as id ->
                    (`Implied_Alignment (file, id, include_header)) :: acc, file
                | _ -> acc, file
            end
        | `GraphicDiagnosis c -> 
            begin match file with
                | None -> (`Diagnosis (c,file)) :: acc, file
                | Some file -> (`GraphicDiagnosis (c,file)) :: acc, Some file
            end
        | `Diagnosis c -> (`Diagnosis (c,file)) :: acc, file
        | `Nodes       -> (`Nodes file) :: acc, file

    let transform_report_arguments x = match x with
        | [`File file] ->
                let file = Some file in
                [`Ascii (file, Some `Single); `Diagnosis (`Normal,file); `Trees ([], file)]
        | [] -> [`Ascii (None, Some `Single); `Diagnosis (`Normal,None); `Trees ([], None)]
        | _  -> let def = [], None in
                let x, _ = List.fold_left transform_report def x in
                x

    (* Selecting *)
    let transform_select (choose, (acc : Methods.script list)) = function
        | `Characters
        | `Taxa as x -> (x, acc)
        | (`Random _) | (`Missing _) | (`Names _) as meth ->
            begin match choose with
                | `Taxa       -> (choose, ((`AnalyzeOnly meth) :: acc))
                | `Characters -> (choose, (`AnalyzeOnlyCharacters meth) :: acc)
            end
        | `Files (do_complement, x) ->
            let x = List.map (fun x -> `Local x) x in
            begin match choose with
                | `Taxa ->
                    (choose, ((`AnalyzeOnlyFiles (do_complement, x)) :: acc))
                | `Characters ->
                    (choose, ((`AnalyzeOnlyCharacterFiles (do_complement, x)) :: acc))
            end
        | `BestN _ | `BestWithin _ | `Unique | `RandomTrees _ as x ->
            (choose, (x :: acc))
        | `AllStatic | `AllDynamic  as x ->
            begin match choose with
                | `Taxa ->
                    let msg = "I@ only@ support@ taxa@ selection@ with@ the names@ of" ^
                        "@ the@ taxa@.@ So@ please,@ if@ it@ is@ worth@ " ^
                        "the@ time@ and@ effort,@ place@ the@ feature@ request."
                    in
                    Status.user_message Status.Error msg;
                    (choose, acc)
                | `Characters ->
                    (choose, ((`AnalyzeOnlyCharacters x) :: acc))
            end
        | _ ->
            let msg = "I@ only@ support@ taxa@ selection@ with@ the names@ of" ^
                "@ the@ taxa@.@ So@ please,@ if@ it@ is@ worth@ " ^
                "the@ time@ and@ effort,@ place@ the@ feature@ request."
            in
            Status.user_message Status.Error msg;
            (choose, acc)

    let transform_select_arguments x = match x with
        | [] -> [`BestN None; `Unique]
        | x  -> snd (List.fold_left transform_select (`Taxa, []) x)

    (* Renaming *)
    let transform_rename (on, (files : Methods.filename list), ren, acc) x = 
        let is_empty = match files, ren with
            | [], [] -> true
            | _ -> false
        in
        match x with
        | `Characters ->
            begin match on with
                | `Characters -> (on, files, ren, acc)
                | `Taxa ->
                    if is_empty then (`Characters, files, ren, acc)
                    else
                        let acc = (`RenameCharacters ren) :: acc in
                        (`Characters, [], [], acc)
            end
        | `Taxa ->
            begin match on with
                | `Taxa -> (on, files, ren, acc)
                | `Characters ->
                    if is_empty then (`Taxa, files, ren, acc)
                    else
                        let acc =
                            (`SynonymsFile (List.rev files))::(`Synonyms ren)::acc
                        in
                        (`Taxa, [], [], acc)
            end
        | `File f ->
            begin match on with
                | `Taxa -> (on, (`Local f) :: files, ren, acc)
                | `Characters ->
                    let msg = "We@ only@ suport@ command@ line@ character@ " ^
                        "renaming@ for@ now.@ If@ you@ think@ this@ is@ a@ useful"
                        ^ "@ feature,@ file@ the@ request.@ For@ now,@ the@ file@ "
                        ^ StatusCommon.escape f ^ "@ will@ be@ ignored."
                    in
                    Status.user_message Status.Error msg;
                    (on, files, ren, acc)
            end
        | `Syn s ->
            (on, files, s :: ren, acc)
        
    let transform_rename_arguments x =
        match List.fold_left transform_rename (`Taxa, [], [], []) x with
        | `Characters, [], [], x 
        | `Taxa, [], [], x -> x
        | `Taxa, files, ren, acc ->
                (`SynonymsFile files) :: (`Synonyms ren) ::  acc
        | `Characters, _, ren, acc ->
                (`RenameCharacters ren) :: acc

    let default_search : Methods.script list = 
        let change_transforms x = { swap_default with Methods.cc = x } in
        let s1 = change_transforms [(`Static_Aprox (`All, false))]
        and s2 = change_transforms [(`Automatic_Sequence_Partition (`All, false, None)); (`Automatic_Static_Aprox false)]
        and s3 = change_transforms [`Automatic_Sequence_Partition (`All, false, None)] in
        List.rev (`Build build_default :: `LocalOptimum s1 :: `LocalOptimum s2 :: [`LocalOptimum s3])

    let transform_search items =
        let do_transform =
            List.exists (function `Transform true -> true | _ -> false) items
        and do_build =
            List.fold_left
                (fun acc ction ->
                    acc && (match ction with `Build false -> false | _ -> true))
                true items
        in
        match default_search with
        | [a; b; c; d] ->
            if not do_build && not do_transform then 
                [`LocalOptimum swap_default]
            else if not do_transform then
                [`LocalOptimum swap_default; d]
            else if not do_build then
                [a; b; c]
            else default_search
        | _ -> failwith "Forgot to update the list of options of search?"


    let process_likelihood_commands est lst =
        let process (alph,cost,model,vari,prior,gaps) = function
            | `ML_alph  x -> ( x,   cost, model, vari, prior, gaps)
            | `ML_cost  x -> (alph,  x,   model, vari, prior, gaps)
            | `ML_subst x -> (alph, cost,  x,    vari, prior, gaps)
            | `ML_vars  x -> (alph, cost, model,  x,   prior, gaps)
            | `ML_prior x -> (alph, cost, model, vari,  x,    gaps)
            | `ML_gaps  x -> (alph, cost, model, vari, prior,  x  )
        in
        let (_,_,model,_,_,_) as tuple =
            List.fold_left process MlModel.default_command lst
        in
        match model with
            | #Methods.ml_optimization          -> `EstLikelihood tuple
            | #Methods.ml_meta                  -> `UseLikelihood tuple
            | #Methods.ml_substitution when est -> `EstLikelihood tuple
            | #Methods.ml_substitution          -> `UseLikelihood tuple

    let organize_read_options lst = 
        let rec organize_data (ss,cm,ro) = function
            ( `Level _
            | `Init3D _
            | `Tie_Breaker _
            | `Affine _
            | `Orientation _) as a -> (ss,cm,a::ro)
            | (`CostMatrix x) ->
                begin match cm with
                    | None -> (ss,Some x,ro)
                    | _  -> failwith "I detected multiple cost matrices specified"
                end
            | (`InputFile x) -> (x::ss, cm, ro)
        in
        List.fold_left organize_data ([],None,[]) lst

    let transform_stdsearch items = 
        `StandardSearch 
            (List.fold_left
                (fun (a, e, b, c, d, f, g) x -> match x with
                    | `MaxTime x -> (Some x, e, b, c, d, f, g)
                    | `MinTime x -> (a, Some x, b, c, d, f, g)
                    | `MinHits x -> (a, e, Some x, c, d, f, g)
                    | `MaxRam x  -> (a, e, b, Some x, d, f, g)
                    | `Target x  -> (a, e, b, c, Some x, f, g)
                    | `Visited x -> (a, e, b, c, d, Some x, g)
                    | `ConstraintFile x -> (a, e, b, c, d, f, Some x))
                (None, None, None, None, None, None, None)
                items)


    let rec transform_command (acc : Methods.script list) (meth : command) : Methods.script list =
        match meth with
        | `Plugin (name, commands) ->
            let rec process_contents x = match x with
                | `Lident _ | `Float _ | `Empty | `Int _ | `String _ as x -> x
                | `Labled (name, args) -> `Labled (name, process_contents args)
                | `List lst -> `List (List.map process_contents lst)
                | `Command c -> 
                    let c = transform_command [] c in
                    `List (List.map (fun x -> `Command x) c)
            in
            (`Plugin (name, process_contents commands)) :: acc
        | `Version
        | `Exit
        | `Recover
        | `ClearRecovered
        | `Redraw
        | `Wipe
        | `ReDiagnose
        | `ReDiagnoseTrees
        | `ClearMemory _
        | `ReadScript _
        | `Store _
        | `Discard _
        | `Load _
        | `Save _
        | `InspectFile _
        | `ChangeWDir _
        | `PrintWDir
        | `Help _ as meth -> meth :: acc
        | `Use x -> `Set x :: acc
        | `Echo (a, []) -> (`Echo (a,`Information)) :: acc
        | `Echo (a, y) -> List.fold_left (fun acc x -> `Echo (a, x) :: acc) acc y
        | `Set settings -> ((List.rev settings) :> Methods.script list)@acc
        | `Read meth -> ((List.rev meth) :> Methods.script list) @ acc
        | `Build x -> (`Build (transform_build_arguments x)) :: acc
        | `Swap x -> (transform_swap_arguments x) :: acc
        | `Search args -> (transform_search args) @ acc
        | `StandardSearch args -> (transform_stdsearch args) :: acc
        | `Fuse x -> (transform_fuse x) :: acc
        | `Support x -> (transform_support_arguments x) :: acc
        | `Calculate _ -> acc
        | `Report x -> (transform_report_arguments x) @ acc
        | `Select x -> (transform_select_arguments x) @ acc
        | `Transform x -> (transform_transform_arguments x) @ acc
        | `Perturb x -> (transform_perturb_arguments x) @ acc
        | `Rename x -> (transform_rename_arguments x) @ acc
        | `Repeat (it, com) ->
            let res = List.rev (List.fold_left transform_command [] com) in
            `Repeat (it, res) :: acc

    let transform_all_commands (x : command list) = 
        List.rev (List.fold_left transform_command [] x)

    let to_local x = List.map (fun x -> `Local x) x

    (* The necessary types to produce the tree of the parsed input. *)
    open Camlp4.PreCast
    module Gram = Camlp4.Struct.Grammar.Static.Make (CommandLexer.Lexer)
    let create_expr () = 
        let expr = Gram.Entry.mk "expr" in
        EXTEND Gram
            GLOBAL: expr;
            expr: [ [a = LIST1 command; `EOI -> (a : command list)] ];
            (* Application commands *)
            (* Transforming taxa or characters *)
            transform:
                [
                    [ LIDENT "transform"; left_parenthesis; x = id_opt_lst; right_parenthesis -> x ]
                ];
            id_opt_lst:
                [
                    [ y = LIST1 [left_parenthesis; x = identifiers_opt; right_parenthesis -> x] SEP "," ->
                        `Transform (List.flatten y) ] |
                    [ y = identifiers_opt -> `Transform y ]
                ];
            identifiers_opt:
                [ [id = identifiers; ","; left_parenthesis; 
                    x = LIST0 [ y = transform_method -> y] SEP ","; right_parenthesis ->
                        (List.map (fun trf -> (id,trf)) x) ] |
                  [ x = LIST0 [ y = transform_method -> y] SEP "," ->
                        (List.map (fun trf -> (`All,trf)) x) ]
                ];
            ml_floatlist:
                [
                    [ ":";left_parenthesis; x = LIST1 integer_or_float SEP ",";
                        right_parenthesis -> List.map float_of_string x ]
                ];
            ml_alphabet: 
                [
                    [ LIDENT "max" -> `Max ] |
                    [ LIDENT "min" -> `Min ] |
                    [ x = INT      -> `Int (int_of_string x) ]
                ];
            ml_substitution: 
                [
                    (* Information Criteria Selection *)
                    [ LIDENT "aic"; ":"; x = STRING  -> `AIC  (Some x) ] |
                    [ LIDENT "bic"; ":"; x = STRING  -> `BIC  (Some x) ] |
                    [ LIDENT "aicc"; ":"; x = STRING -> `AICC (Some x) ] |
                    [ LIDENT "aic"  -> `AIC  None ] |
                    [ LIDENT "bic"  -> `BIC  None ] |
                    [ LIDENT "aicc" -> `AICC None ] |
                    (* Meta Likelihood; transforms as a special command *)
                    [ LIDENT "ncm" -> `NCM ] |
                    (* Standard likelihood transformation models *)
                    [ LIDENT "jc69" -> `JC69 ] |
                    [ LIDENT "f81"  -> `F81  ] |
                        (* values of these types get checked later *)
                    [ LIDENT "f84";   x = OPT ml_floatlist ->
                        let x = (function None -> [] | Some x -> x) x in `F84 x  ] |
                    [ LIDENT "k80";   x = OPT ml_floatlist ->
                        let x = (function None -> [] | Some x -> x) x in `K2P x  ] |
                    [ LIDENT "k2p";   x = OPT ml_floatlist ->
                        let x = (function None -> [] | Some x -> x) x in `K2P x  ] |
                    [ LIDENT "hky";   x = OPT ml_floatlist ->
                        let x = (function None -> [] | Some x -> x) x in `HKY85 x] |
                    [ LIDENT "hky85"; x = OPT ml_floatlist ->
                        let x = (function None -> [] | Some x -> x) x in `HKY85 x] |
                    [ LIDENT "tn93";  x = OPT ml_floatlist ->
                        let x = (function None -> [] | Some x -> x) x in `TN93 x ] |
                    [ LIDENT "gtr";   x = OPT ml_floatlist ->
                        let x = (function None -> [] | Some x -> x) x in `GTR x  ] |
                    [ LIDENT "file"; ":"; x = STRING   -> `File x ] |
                    [ LIDENT "custom"; ":"; x = STRING -> `Custom x ]
                ];
            site_properties:
                [
                    [ right_parenthesis -> None ] |
                    [ ","; y = LIST1 integer_or_float SEP ",";right_parenthesis ->
                        Some (List.map float_of_string y) ]
                ];
            ml_rates: 
                [
                    [ LIDENT "gamma";":";left_parenthesis;x = INT; d = site_properties -> 
                        let d = match d with | Some [x] -> Some x | None -> None
                            | _ -> failwith "Improper Gamma Argument" in
                        Some (`Gamma (int_of_string x, d)) ] |
                    [ LIDENT "theta";":";left_parenthesis;x = INT; d = site_properties -> 
                        let d = match d with | Some [x;y] -> Some (x,y) | None -> None
                            | _ -> failwith "Improper Theta Argument" in
                            Some (`Theta (int_of_string x, d)) ] |
                    [ LIDENT "constant" -> None ] |
                    [ LIDENT "none" -> None ]
                ];
            ml_priors:
                [ 
                    [ LIDENT "estimate" -> `Estimate ] |
                    [ LIDENT "equal" -> `Equal ] | 
                    [ LIDENT "given"; ":"; left_parenthesis; 
                        x = LIST1 [ x = FLOAT -> float_of_string x ] SEP ",";
                        right_parenthesis -> `Given x ]
                ];
            ml_gap :
                [ [ LIDENT "missing"   -> `Missing ] |
                  [ LIDENT "character" -> `Independent] |
                  [ LIDENT "coupled"   -> `Coupled 1.0 ] |
                  [ LIDENT "coupled"; ":"; x = integer_or_float -> `Coupled (float_of_string x) ]
                ];
            ml_costfn:
                [
                    [LIDENT "mal" -> `MAL] | [LIDENT "mpl" -> `MPL]
                ];
            partitioned_mode:
                [
                    [ LIDENT "clip" -> `Clip ] | [ LIDENT "noclip" -> `NoClip ]
                ];
            ml_properties:
                [
                    [ x = ml_substitution -> `ML_subst x ] |
                    [ LIDENT "alphabet"; ":"; x = ml_alphabet -> `ML_alph x ] |
                    [ LIDENT "rates";":"; x = ml_rates -> `ML_vars x ] |
                    [ x = ml_costfn -> `ML_cost x ] |
                    [ LIDENT "priors";":"; x = ml_priors -> `ML_prior x ] |
                    [ LIDENT "gap"; ":"; x = ml_gap -> `ML_gaps x ]
                ];
            optional_poly :
                [
                    [","; x = polymorphism_parameter -> x ]
                ];
            polymorphism_parameter:
                [
                    [ LIDENT "ignore_polymorphism" -> `Do_Nothing ] |
                    [ LIDENT "simple_polymorphism" -> `Pick_One] |
                    [ LIDENT "full_polymorphism" -> `Do_All]
                ];
            fixed_states_option :
                [
                   [":"; left_parenthesis; 
                        x = OPT [z=STRING->z]; y = OPT optional_poly;
                   right_parenthesis -> (x,y) ]
                ];
            level_and_tiebreaker :
                [
                    [ left_parenthesis; x = INT; ","; y = keep_method; right_parenthesis
                    ->(int_of_string x,y) ]|
                    [ x=INT -> (int_of_string x,`First) ]
                ];
            transform_method:
                [
                    [ LIDENT "origin_cost"; ":"; x = integer_or_float ->
                        `OriginCost (float_of_string x) ] |
                    [ LIDENT "prioritize" -> `Prioritize ] |
                    [ LIDENT "elikelihood"; ":"; left_parenthesis;
                        lst = LIST1 [x = ml_properties -> x] SEP ","; right_parenthesis ->
                            process_likelihood_commands true lst ] |
                    [ LIDENT "parsimony" -> `UseParsimony ] |
                    [ LIDENT "likelihood"; ":"; left_parenthesis;
                        lst = LIST1 [x = ml_properties -> x] SEP ","; right_parenthesis ->
                            process_likelihood_commands false lst ] |
                    [ LIDENT "prealigned" -> `Prealigned_Transform ] |
                    [ LIDENT "randomize_terminals" -> `RandomizedTerminals ] |
                    [ LIDENT "alphabetic_terminals" -> `AlphabeticTerminals ] |
                    [ LIDENT "level"; ":"; x = level_and_tiebreaker -> `Level x ] |
                    [ LIDENT "tcm"; ":"; left_parenthesis; 
                         x = tcm_arguments; right_parenthesis -> x ] |
                    [ LIDENT "partitioned"; ":"; x = partitioned_mode -> `Partitioned x ] |  
                    [ LIDENT "fixed_states"; x = OPT fixed_states_option -> match x with
                        | Some y -> `Fixed_States y
                        | None -> `Fixed_States (None,None)
                    ] | 
                    [ LIDENT "direct_optimization" -> `Direct_Optimization ] |
                    [ LIDENT "do" -> `Direct_Optimization ] |
                    [ LIDENT "gap_opening"; ":"; x = INT -> `AffGap (int_of_string x) ] |
                    [ LIDENT "trailing_deletion"; ":"; x = STRING -> `TailFile x ] |
                    [ LIDENT "td"; ":"; x = STRING -> `TailFile x ] |
                    [ LIDENT "trailing_deletion"; ":"; left_parenthesis; 
                        x = LIST1 [ x = INT -> x]  SEP ","; right_parenthesis -> 
                            `TailInput (Array.of_list (List.map (int_of_string) x)) ] |
                    [ LIDENT "td"; ":"; left_parenthesis; x = LIST1 [ x = INT -> x] SEP ",";
                        right_parenthesis -> 
                            `TailInput (Array.of_list (List.map (int_of_string) x)) ] |
                    [ LIDENT "trailing_insertion"; ":"; x = STRING -> `PrepFile x ] |
                    [ LIDENT "ti"; ":"; x = STRING -> `PrepFile x ] |
                    [ LIDENT "trailing_insertion"; ":"; left_parenthesis; 
                        x = LIST1 [ x = INT -> x] SEP ","; right_parenthesis -> 
                            `PrepInput (Array.of_list (List.map (int_of_string) x)) ] |
                    [ LIDENT "ti"; ":"; left_parenthesis; x = LIST1 [ x = INT -> x] SEP ",";
                        right_parenthesis -> 
                            `PrepInput (Array.of_list (List.map (int_of_string) x)) ] |
                    [ LIDENT "static_approx"; x = OPT uninformative_characters -> match x with 
                        | None -> `StaticApproximation true 
                        | Some v -> `StaticApproximation v ] |
                    [ LIDENT "multi_static_approx"; x = OPT uninformative_characters -> match x with 
                        | None -> `MultiStaticApproximation true 
                        | Some v -> `MultiStaticApproximation v ] |
                    [ LIDENT "auto_static_approx"; x = OPT optional_boolean -> match x with
                        | None -> `Automatic_Static_Aprox false 
                        | Some x -> `Automatic_Static_Aprox x ] |
                    [ LIDENT "auto_sequence_partition"; x = OPT optional_boolean -> match x with
                        | None -> `Automatic_Sequence_Partition (false, None)
                        | Some x -> `Automatic_Sequence_Partition (x, None) ] |
                    [ LIDENT "sequence_partition"; ":"; x = INT -> 
                        `Automatic_Sequence_Partition (false, Some (int_of_string x)) ] |
                    [ LIDENT "weight"; ":"; x = neg_integer_or_float -> `ReWeight (float_of_string x) ] |
                    [ LIDENT "weightfactor"; ":"; x = neg_integer_or_float -> `WeightFactor (float_of_string x) ] |
                    [ LIDENT "search_based"; ":"; left_parenthesis; file_couple_lst = LIST0
                        [left_parenthesis; chname = STRING; ","; filename = STRING;
                            right_parenthesis -> (chname,filename)] SEP ","; 
                        right_parenthesis -> `SearchBased file_couple_lst ] |
                    [ LIDENT "seq_to_chrom"; ":"; left_parenthesis; x = LIST0
                            [ x = chromosome_argument -> x] SEP ","; right_parenthesis -> `SeqToChrom x ] | 
                    [ LIDENT "annchrom_to_breakinv"; ":"; left_parenthesis; x = LIST0
                            [x = chromosome_argument -> x] SEP ","; right_parenthesis -> `AnnchromToBreakinv x ] |
                    [ LIDENT "custom_alphabet"; ":"; left_parenthesis; x = LIST0 
                            [ x = chromosome_argument -> x] SEP ","; right_parenthesis -> `ChangeDynPam x] |
                    [ LIDENT "breakinv"; ":"; left_parenthesis; x = LIST0 
                            [ x = chromosome_argument -> x] SEP ","; right_parenthesis -> `ChangeDynPam x ] |
                    [ LIDENT "chromosome"; ":"; left_parenthesis; x = LIST0 
                            [ x = chromosome_argument -> x] SEP ","; right_parenthesis -> `ChangeDynPam x ] |
                    [ LIDENT "annotated"; ":"; left_parenthesis; x = LIST0 
                            [ x = chromosome_argument -> x] SEP ","; right_parenthesis -> `ChangeDynPam x ] |
                    [ LIDENT "genome"; ":"; left_parenthesis; x = LIST0 
                            [ x = genome_argument -> x] SEP ","; right_parenthesis -> `ChangeDynPam x ] |
                    [ LIDENT "chrom_to_seq" -> `ChromToSeq [] ] |
                    [ LIDENT "breakinv_to_custom" -> `BreakinvToSeq [] ]
    (*              | [ LIDENT "kolmogorov"; y = OPT optional_kolmogorov_parameters ->*)
    (*                    match y with*)
    (*                    | None -> `Seq_to_Kolmogorov (`AtomicIndel (None, None))*)
    (*                    | Some x -> x ]*)

                ];
            tcm_arguments:
                [
                    [ x = INT; ","; y = INT ->
                        `Gap (int_of_string x, int_of_string y)] |
                    [ x = STRING; level_value = OPT optional_level ->
                        match level_value with 
                            | None   -> `Tcm (x,None)
                            | Some y -> `Tcm (x,Some y)
                    ]
                ];
    (*        optional_kolmogorov_parameters: *)
    (*            [ *)
    (*                [ ":"; left_parenthesis; *)
    (*                    x = LIST0 [ x = kolmogorov_parameters -> x] SEP ","; *)
    (*                    right_parenthesis  -> *)
    (*                        let default = (None, None) in*)
    (*                        let x = *)
    (*                            List.fold_left *)
    (*                                (fun acc x -> match x with *)
    (*                                    | `Event n -> (Some n, snd acc)*)
    (*                                    | `IndelSub n -> (fst acc, Some n))*)
    (*                                default *)
    (*                                x*)
    (*                        in*)
    (*                        `Seq_to_Kolmogorov (`AtomicIndel x) ] *)
    (*            ];*)
    (*        kolmogorov_parameters:*)
    (*            [*)
    (*                [ LIDENT "event"; ":"; x = FLOAT -> *)
    (*                    `Event (float_of_string x) ] | *)
    (*                [ LIDENT "indelsub"; ":"; left_parenthesis; insertion = FLOAT;*)
    (*                    ","; deletion = FLOAT; ","; substitution = FLOAT;*)
    (*                    right_parenthesis -> *)
    (*                        `IndelSub*)
    (*                        (float_of_string insertion,*)
    (*                        float_of_string deletion, *)
    (*                        float_of_string substitution) ]*)
    (*            ];*)
            uninformative_characters:
                [
                    [ ":"; LIDENT "keep" -> false ] |
                    [ ":"; LIDENT "remove" -> true ]
                ];
            median_solvers:
                [
                    [ LIDENT "mgr" -> `MGR   ] |
                    [ LIDENT "simplelk" -> `SimpleLK   ] |
                    [ LIDENT "chainedlk" -> `ChainedLK ] |
                    [ LIDENT "coalestsp" -> `COALESTSP ] |
                    [ LIDENT "bbtsp" -> `BBTSP         ] |
                    [ LIDENT "siepel" -> `Siepel       ] |
                    [ LIDENT "default" -> `Albert      ] |
                    [ LIDENT "caprara" -> `Albert      ] |
                    [ LIDENT "vinh" -> `Vinh           ]
                ];
            annotate_param:
                [
                    [ LIDENT "vinh"; ","; x = INT; ","; y = INT; ","; z = INT ->
                        `Default (int_of_string x,int_of_string y, int_of_string z) ] |
                    [ LIDENT "mauve"; ","; a = FLOAT; ","; b = FLOAT; ","; c = FLOAT; ","; d = FLOAT ->
                        `Mauve (float_of_string a,float_of_string b, float_of_string c, float_of_string d) ]
                ];
            genome_argument:
                [
                    [ LIDENT "translocation"; ":"; c = INT -> 
                          `Translocation (int_of_string c) ]  |
                    [ LIDENT "chrom_indel"; ":"; left_parenthesis; o = INT; 
                    ","; e = integer_or_float; right_parenthesis ->
                          `Chrom_Indel_Cost ( (int_of_string o), 
                          int_of_float ((float_of_string e) *. 100.0) ) ] | 
                    [ LIDENT "chrom_hom"; ":"; c = FLOAT -> 
                          `Chrom_Hom (int_of_float ((float_of_string c) *. 100.)) ] |
                    [ LIDENT "median"; ":"; c = INT -> `Keep_Median (int_of_string c) ] 
                ];
            chromosome_argument:
                [
                    [ LIDENT "circular"; ":"; e = boolean -> `Circular e] |
                    [ LIDENT "median_solver"; ":"; c = median_solvers ->
                        match c with
                        | `MGR -> `Median_Solver `MGR
                        | `SimpleLK -> `Median_Solver `SimpleLK
                        | `ChainedLK -> `Median_Solver `ChainedLK
                        | `COALESTSP -> `Median_Solver `COALESTSP
                        | `BBTSP -> `Median_Solver `BBTSP
                        | `Siepel -> `Median_Solver `Siepel                     
                        | `Albert -> `Median_Solver `Albert 
                        | `Vinh -> `Median_Solver `Vinh
                    ]|
                    [ LIDENT "annotate"; ":"; left_parenthesis;
                       c = annotate_param; right_parenthesis
                       -> `Annotate_Tool c
                    ]|
                    [ LIDENT "locus_inversion"; ":"; c = INT -> 
                          `Locus_Inversion (int_of_string c) ]  |
                    [ LIDENT "locus_breakpoint"; ":"; c = INT -> 
                          `Locus_Breakpoint (int_of_string c) ]  |
                    [ LIDENT "locus_indel"; ":"; left_parenthesis; o = INT; 
                        ","; e = integer_or_float; right_parenthesis ->
                          `Locus_Indel_Cost ( (int_of_string o), 
                          int_of_float ((float_of_string e) *. 100.0) ) ] | 
                    [ LIDENT "swap_med"; ":"; iters = INT -> `SwapMed (int_of_string iters) ] | 
                    [ LIDENT "median"; ":"; c = INT -> `Keep_Median (int_of_string c) ] |
                    [ LIDENT "med_approx"; ":"; ans = boolean -> `Approx ans] |
                    [ LIDENT "symmetric"; ":"; ans = boolean -> `Symmetric ans] |
                    [ LIDENT "max_3d_len"; ":"; l = INT -> 
                          `Max_3D_Len (int_of_string l) ]  |
                    [ LIDENT "max_kept_wag"; ":"; l = INT -> 
                          `Max_kept_wag (int_of_string l) ]  
                ];

            (* Applications *)
            application_command:
                [
                    [ LIDENT "version"; left_parenthesis; right_parenthesis -> `Version ] |
                    [ LIDENT "exit"; left_parenthesis; right_parenthesis -> `Exit ] |
                    [ LIDENT "recover"; left_parenthesis; right_parenthesis -> `Recover ] |
                    [ LIDENT "clear_recovered"; left_parenthesis; right_parenthesis -> 
                        `ClearRecovered ] |
                    [ LIDENT "quit" ; left_parenthesis; right_parenthesis -> `Exit ] |
                    [ LIDENT "echo"; left_parenthesis; a = STRING; OPT ",";
                        x = LIST0 [x = output_class -> x ] SEP ","; right_parenthesis -> `Echo (a, x) ] |
                    [ LIDENT "help"; left_parenthesis; a = OPT string_or_ident; 
                        right_parenthesis -> `Help a ] |
                    [ LIDENT "set"; left_parenthesis; b = LIST0 [x = setting -> x] SEP ","; 
                        right_parenthesis -> `Set b ] |
                    [ LIDENT "redraw"; left_parenthesis; right_parenthesis -> `Redraw ] |
                    [ LIDENT "wipe"; left_parenthesis; right_parenthesis -> `Wipe ] |
                    [ LIDENT "clear_memory"; left_parenthesis; x = LIST0 [x =
                        clear_options -> x]
                        SEP ","; right_parenthesis -> `ClearMemory x ] |
                    [ LIDENT "load"; left_parenthesis; a = STRING; 
                        right_parenthesis -> `Load a ] |
                    [ LIDENT "save"; left_parenthesis; a = STRING; 
                        b = OPT file_comment; right_parenthesis -> `Save (a, b) ] |
                    [ LIDENT "inspect"; left_parenthesis; a = STRING; 
                        right_parenthesis -> `InspectFile a ] |
                    [ LIDENT "store"; left_parenthesis; a = STRING; y = OPT store_class; 
                    right_parenthesis -> 
                            let st =
                                match y with
                                | None -> all_store_types
                                | Some x -> x
                            in
                            `Store (st, a) ] |
                    [ LIDENT "use"; left_parenthesis; a = STRING; y = OPT
                    store_class; right_parenthesis -> 
                            let st = 
                                match y with
                                | None -> all_store_types
                                | Some x -> x 
                            in
                            `Use (st, a) ] |
                    [ LIDENT "rediagnose"; left_parenthesis; 
                        x = OPT clear_diagnosis; right_parenthesis -> 
                            match x with
                            | None  
                            | Some false -> `ReDiagnose
                            | Some true  -> `ReDiagnoseTrees ] |
                    [ LIDENT "run"; left_parenthesis; a = LIST0 [x = STRING -> x] SEP ","; 
                        right_parenthesis -> `ReadScript a ] |
                    [ LIDENT "cd"; left_parenthesis; a = STRING; right_parenthesis ->
                        `ChangeWDir a ] |
                    [ LIDENT "pwd"; left_parenthesis; right_parenthesis -> 
                        `PrintWDir ]
                ];
            clear_diagnosis:
                [ [ LIDENT "clear" -> false ] | [ LIDENT "preserve" -> true ] ];
            store_class:
                [ [ ","; left_parenthesis; x = LIST0 [ x = store_class_list -> x ]
                    SEP ","; right_parenthesis -> x ] ];
            store_class_list:
                [ 
                    [ LIDENT "data" -> `Data ] |
                    [ LIDENT "trees" -> `Trees ] |
                    [ LIDENT "jackknife" -> `Jackknife ] | 
                    [ LIDENT "bremer" -> `Bremer ] |
                    [ LIDENT "bootstrap" -> `Bootstrap ] 
                ];
            clear_options:
                [
                    [ LIDENT "m" -> `Matrices ] |
                    [ LIDENT "s" -> `SequencePool ]
                ];
            file_comment:
                [
                    [ ","; x = STRING -> x ]
                ];
            output_class:
                [
                    [ LIDENT "info" -> `Information ] |
                    [ LIDENT "error" -> `Error ] |
                    [ LIDENT "output"; x = OPT optional_string -> `Output x ]
                ];
            optional_string:
                [ [  ":"; x = STRING -> x ] ];
            optional_level:
                [ [ ","; LIDENT "level"; ":"; y=level_and_tiebreaker -> y ] ];
            setting:
                [
                    [ LIDENT "timer"; ":"; x = INT -> `TimerInterval (int_of_string x) ] |
                    [ LIDENT "history"; ":"; x = INT -> `HistorySize (int_of_string x) ] |
                    [ LIDENT "log"; ":"; x = STRING -> `Logfile (Some x) ] |
                    [ LIDENT "log"; ":"; LIDENT "new"; ":"; x = STRING ->
                        StatusCommon.Files.closef x ();
                        ignore (StatusCommon.Files.openf ~mode:`New x []);
                        `Logfile (Some x) ] |
                    [ LIDENT "nolog" -> `Logfile None ] |
                    [ LIDENT "seed"; ":"; x = neg_integer -> `SetSeed (int_of_string x) ] |
                    [ LIDENT "root"; ":"; x = STRING -> `RootName x ] |
                    [ LIDENT "root"; ":"; x = INT -> `Root (Some (int_of_string x)) ] |
                    [ LIDENT "exhaustive_do"    -> `Exhaustive_Weak ] |
                    [ LIDENT "iterative"; ":"; x = iterative_mode
                                                -> `Iterative x ] |
                    [ LIDENT "normal_do"        -> `Normal ] |
                    [ LIDENT "normal_do_plus"   -> `Normal_plus_Vitamines ] |
                    [ LIDENT "opt"; ":"; x = opt_level   -> `Optimization x ] |
                    [ LIDENT "codon_partition"; ":"; 
                        left_parenthesis; n = STRING; ","; x = old_identifiers; right_parenthesis ->
                            `Alias (n,`Codon x) ] |
                    [ LIDENT "partition"; ":"; 
                        left_parenthesis; n = STRING; ","; x = old_identifiers; right_parenthesis ->
                            `Alias (n,`Chars x) ] |
                    [ LIDENT "parmap"; ":"; x = INT -> `Parmap (int_of_string x) ]
    (*                [ LIDENT "space_saving_alignment" -> `Algn_Newkk ]|*)
    (*                [ LIDENT "normal_alignment" -> `Algn_Normal ]*)
                ];
            neg_integer:
                [
                    [ "-"; x = INT -> "-" ^ x ] |
                    [ x = INT -> x ]
                ];
            opt_level:
                [
                    [ LIDENT "exhaustive";  iterations = OPT optional_integer_or_float ->
                        let iterations = match iterations with
                            | None -> None | Some x -> Some (int_of_string x)
                        in
                        `Exhaustive iterations ] |
                    [ LIDENT "coarse" ;     iterations = OPT optional_integer_or_float ->
                        let iterations = match iterations with
                            | None -> None | Some x -> Some (int_of_string x)
                        in
                        `Coarse iterations ] |
                    [ LIDENT "no_opt" -> `None ]
                ];
            iterative_mode:
                [ 
                    [ LIDENT "exact"; iterations = OPT optional_integer_or_float  ->
                        let iterations = match iterations with 
                            None -> None | Some x -> Some (int_of_string x) in
                        `ThreeD iterations ] |
                    [ LIDENT "approximate"; iterations = OPT optional_integer_or_float -> 
                        let iterations = match iterations with
                            | None -> None | Some x -> Some (int_of_string x) in
                        `ApproxD iterations ]
                ];
            (* Reporting *)
            report:
                [
                    [ LIDENT "report"; 
                        left_parenthesis; a = LIST0 [x = report_argument -> x] SEP ","; 
                        right_parenthesis -> `Report a ]
                ];
            topology_test_params : 
                [   [ LIDENT "sh" -> `SH ] |
                    [ LIDENT "kh" -> `KH ] |
                    [ LIDENT "n"; ":"; x = INT -> `Replicates (int_of_string x) ] |
                    [ LIDENT "k"; ":"; x = INT -> `ScaleFactors (int_of_string x) ] |
                    [ LIDENT "rell" -> `ReplicateOpt (false, false) ] |
                    [ LIDENT "full" -> `ReplicateOpt (true, true) ] |
                    [ LIDENT "part" -> `ReplicateOpt (false, true) ] |
                    [ x = old_identifiers -> `Characters x]
                ];
            topology_tests:
                [ [ left_parenthesis;
                        x = LIST0 [x = topology_test_params -> x] SEP ",";
                    right_parenthesis -> x ]
                ];
            report_argument:
                [
                    [ x = STRING -> `File x ] |
                    [ LIDENT "kml"; ":"; left_parenthesis; 
                        plugin = LIDENT; ","; csv = STRING;
                    right_parenthesis -> `KML (plugin, csv) ] |
                    [ LIDENT "new"; ":"; x = STRING ->
                        StatusCommon.Files.closef x ();
                        ignore (StatusCommon.Files.openf ~mode:`New x []);
                        `File x ] |
                    [ LIDENT "asciitrees" ; y = OPT optional_collapse -> 
                        match y with
                        | Some (`Collapse y) -> `Ascii y
                        | None -> `Ascii (Some `Single) ] |
                    [ LIDENT "memory" -> `Memory ] |
                    [ LIDENT "graphtrees"; y = OPT optional_collapse ->
                        match y with
                            | Some (`Collapse y) -> `Graph y
                            | None -> `Graph (Some `Single) ] |
                    [ LIDENT "trees"; x = OPT tree_information_list ->
                        match x with | Some x -> `Trees x | None -> `Trees [] ] |
                    [ LIDENT "treestats" -> `TreesStats ] |
                    [ LIDENT "searchstats" -> `SearchStats ] |
                    [ LIDENT "treecosts" -> `TreeCosts ] |
    (*                [ LIDENT "kolmo_machine" -> `KolmoMachine ] |*)
                    [ LIDENT "timer"; ":"; x = STRING -> `TimeDelta x ] |
                    [ LIDENT "_mst" -> `MstR ] |
                    [ LIDENT "consensus"; x = OPT optional_integer_or_float ->
                        match x with
                            | None -> `Consensus None
                            | Some x -> `Consensus (Some (float_of_string x)) ] |
                    [ LIDENT "graphconsensus"; x = OPT optional_integer_or_float ->
                        match x with
                            | None -> `GraphicConsensus None
                            | Some x -> `GraphicConsensus (Some (float_of_string x))] |
                    [ LIDENT "clades" -> `Clades ] |
                    [ LIDENT "phastwinclad" -> `FasWinClad ] | 
                    [ LIDENT "nexus" -> `Nexus ] | 
                    [ LIDENT "lkmodel"; ":"; x = old_identifiers -> `Model x ] | 
                    [ LIDENT "lkmodel" -> `Model `All ] | 
                    [ LIDENT "lksites"; ":"; x = old_identifiers -> `LKSites x ] | 
                    [ LIDENT "lksites" -> `LKSites `All ] | 
                    [ LIDENT "debug_data" -> `DebugData ] | 
    (*                [ LIDENT "script" -> `Script (!console_script) ] |*)
    (*                [ LIDENT "pairwise"; ":"; x = old_identifiers -> `Pairwise x] |*)
    (*                [ LIDENT "pairwise" -> `Pairwise `All ] |*)
                    [ LIDENT "seq_stats"; ":"; ch = old_identifiers -> `SequenceStats ch ] |
                    [ LIDENT "seq_stats" -> `SequenceStats `All ] |
                    [ LIDENT "ci"; ":"; ch = old_identifiers -> `Ci (Some ch) ] |
                    [ LIDENT "ci" -> `Ci None ] |
                    [ LIDENT "ri"; ":"; ch = old_identifiers -> `Ri (Some ch) ] |
                    [ LIDENT "ri" -> `Ri None ] |
                    [ LIDENT "compare"; ":"; left_parenthesis; complement = boolean;
                        ","; ch1 = old_identifiers; ","; ch2 = old_identifiers;
                        right_parenthesis -> `CompareSequences (complement, ch1, ch2) ] |
                    [ LIDENT "script_analysis"; ":"; x = STRING -> `ExplainScript x ] |
                    [ LIDENT "supports"; y = OPT opt_support_names -> `Supports y ] |
                    [ LIDENT "topologytest"; ":"; y = topology_tests -> `Topo_Selection y] |
                    [ LIDENT "graphsupports"; y = OPT opt_support_names -> `GraphicSupports y ] |
                    [ LIDENT "diagnosis"; y = OPT opt_report_type -> match y with 
                        | None -> `Diagnosis `Normal
                        | Some x -> `Diagnosis x] |
                    [ LIDENT "graphdiagnosis";  y = OPT opt_report_type -> match y with
                        | None -> `GraphicDiagnosis `Normal
                        | Some x -> `GraphicDiagnosis x ] |
                    [ LIDENT "data" -> `Data ] |
                    [ LIDENT "xslt"; ":"; "("; a = STRING; ","; b = STRING; ")" -> `Xslt (a, b) ] |
                    [ LIDENT "implied_alignments"; ":"; x = identifiers -> `Implied_Alignments (x, true) ] |
                    [ LIDENT "fasta"; ":"; x = identifiers -> `Implied_Alignments (x, false) ] |
                    [ LIDENT "fasta" -> `Implied_Alignments (`All, false) ] |
                    [ LIDENT "all_roots" -> `AllRootsCost ] |
                    [ LIDENT "implied_alignments" -> `Implied_Alignments (`All, true)] |
                    [ LIDENT "ia"; ":"; x = identifiers -> `Implied_Alignments (x, true) ] | 
                    [ LIDENT "ia" -> `Implied_Alignments (`All, true) ] |
                    [ LIDENT "nodes" -> `Nodes ] |
                    [ LIDENT "cross_references"; ":"; x = old_identifiers -> `CrossReferences (Some x) ] |
                    [ LIDENT "terminals" -> `TerminalsFiles ] | 
                    [ LIDENT "cross_references" -> `CrossReferences None ] |
                    [ LIDENT "robinson_foulds" -> `RobinsonFoulds ]
                ];
            (* Perturbation method *)
            perturb:
                [
                    [ LIDENT "perturb"; left_parenthesis; x = LIST0 [x =
                        perturb_argument -> x] SEP ","; 
                        right_parenthesis -> `Perturb x ]
                ];
            perturb_argument:
                [
                    [ x = ratchet -> x ]  |
                    [ x = resample -> x ] |
                    [ x = swap -> (x :> perturba) ] |
                    [ x = transform -> (x :> perturba) ] |
                    [ LIDENT "timeout"; ":"; x = integer_or_float -> `TimeOut
                    (`Fixed (float_of_string x))] |
                    [ LIDENT "iterations"; ":"; x = INT -> `Iterations (int_of_string x) ]
                ];
            ratchet:
                [
                    [ LIDENT "ratchet"; ":"; left_parenthesis; x = FLOAT; ","; y = INT; 
                        right_parenthesis -> `Ratchet (float_of_string x,
                        int_of_string y) ] |
                    [ LIDENT "ratchet" -> pert_def_perturb ]
                ];
            resample:
                [
                    [ LIDENT "resample"; ":"; x = INT -> `Resample (int_of_string x) ]
                ];
            charortax:
                [
                    [ LIDENT "characters" -> `Characters ] |
                    [ LIDENT "terminals" -> `Taxa ]
                ];
            (* Selecting characters or taxa *)
            select:
                [
                    [ LIDENT "select"; left_parenthesis; x = LIST0 [x =
                        select_argument -> x] SEP ","; 
                        right_parenthesis -> `Select x ]
                ];
            select_argument:
                [
                    [ x = identifiers -> (x :> selecta) ] |
                    [ x = charortax -> (x :> selecta) ] |
                    [ x = seltrees -> (x :> selecta) ] |
                    [ x = STRING -> `Files (true, [x]) ] 
                ];
            seltrees:
                [
                    [ LIDENT "optimal" -> `BestN None ] |
                    [ LIDENT "unique" -> `Unique ] |
                    [ LIDENT "best"; ":"; x = INT -> `BestN (Some (int_of_string x)) ] |
                    [ LIDENT "within"; ":"; x = integer_or_float -> 
                        `BestWithin (float_of_string x) ] |
                    [ LIDENT "random"; ":"; x = INT -> `RandomTrees (int_of_string x) ]
                ];
            (* Renaming characters or taxa *)
            rename:
                [
                    [ rename_cmd; left_parenthesis; x = LIST0 [x = rename_argument
                    -> x] SEP ","; 
                        right_parenthesis -> `Rename x ]
                ];
            rename_cmd:
                [
                    [ LIDENT "synonymize" ] |
                    [ LIDENT "rename" ]
                ];
            rename_argument:
                [
                    [ x = charortax -> (x :> renamea) ] |
                    [ x = STRING -> `File x ] | 
                    [ left_parenthesis; a = STRING; ","; b = STRING; right_parenthesis -> 
                        `Syn (a, b) ]
                ];
            (* POY commands *)
            command:
                [
                    [ t = read -> (t :> command) ] |
                    [ t = build -> (t :> command) ] |
                    [ t = swap -> (t :> command) ] |
                    [ t = search -> (t :> command) ] |
                    [ t = calculate_support -> (t :> command) ] |
                    [ t = perturb -> (t :> command) ] | 
                    [ t = transform -> (t :> command) ] |
                    [ t = report -> (t :> command) ] |
                    [ t = select -> (t :> command) ] |
                    [ t = rename -> (t :> command) ] |
                    [ t = application_command -> (t :> command) ] |
                    [ t = fuse -> (t :> command) ] |
                    [ t = loop -> (t :> command) ] |
                    [ t = plugin -> (t :> command) ]
                ];
            plugin: 
                [
                    [ x = LIDENT; left_parenthesis; 
                        y = LIST0 [ x = plugin_args -> x] SEP ",";
                            right_parenthesis -> 
                        match y with
                        | [] -> `Plugin (x, `Empty)
                        | [y] -> `Plugin (x, y) 
                        | y -> `Plugin (x, `List y) ]
                ];
            plugin_args:
                [   
                    [ x = FLOAT -> `Float (float_of_string x) ] |
                    [ x = INT -> `Int (int_of_string x) ] |
                    [ x = STRING -> `String x ] |
                    [ x = LIDENT; ":"; y = plugin_args -> `Labled (x, y)] |
                    [ x = LIDENT -> `Lident x ] | 
                    [ left_parenthesis; x = LIST0 [ x = plugin_args -> x] SEP ",";
                    right_parenthesis -> `List x ]
                ];
            loop:
                [
                    [ LIDENT "repeat"; x = INT; LIDENT "times"; com = LIST0 [x =
                        command -> x]; 
                        LIDENT "end" -> `Repeat (int_of_string x, com) ] 
                ];
            read:
                [
                    [ LIDENT "read"; left_parenthesis; a = LIST0 [x = read_argument
                    -> x] SEP ","; 
                        right_parenthesis -> `Read a ]
                ];
            build:
                [
                    [ LIDENT "build"; left_parenthesis; a = LIST0 [x =
                        build_argument -> x] SEP ","; 
                        right_parenthesis -> `Build a ]
                ];
            swap:
                [
                    [ LIDENT "swap"; left_parenthesis; a = LIST0 [x = swap_argument
                    -> x] SEP ","; 
                        right_parenthesis -> (`Swap a :> swap) ]
                ];
            time:
                [
                    [ days = integer_or_float; ":"; hours = integer_or_float; ":";
                    minutes = integer_or_float ->
                        (int_of_float (((float_of_string days) *. 60. *. 60. *. 24.) +.
                        ((float_of_string hours) *. 60. *. 60.) +.
                        ((float_of_string minutes) *. 60. ))) ] 
                ];
            memory:
                [
                    [ LIDENT "gb"; ":"; x = INT -> ((int_of_string x) * 
                        1000 * 1000 * (1000 / (Sys.word_size / 8))) ] |
                    [ LIDENT "mb"; ":"; x = INT ->((int_of_string x) * 
                        1000 * 1000 / (Sys.word_size / 8)) ]
                ];
            model_iter :
                [
                    [LIDENT "never" -> `NullModel] |
                    [LIDENT "always" -> `AlwaysModel] |
                    [LIDENT "threshold"; ":"; x = FLOAT -> `ThresholdModel (float_of_string x) ]|
                    [LIDENT "max_count"; ":"; x = INT -> `MaxCountModel (int_of_string x) ]|
                    [LIDENT "never" -> `NullModel] |
                    [LIDENT "always" -> `AlwaysModel]
                ];
            branch_iter :
                [
                    [LIDENT "never"        -> `NullBranches] |
                    [LIDENT "all_branches" -> `AllBranches] |
                    [LIDENT "all"          -> `AllBranches] |
                    [LIDENT "join_delta"   -> `JoinDeltaBranches] |
                    [LIDENT "join_region"  -> `NeighborhoodBranches 0] |
                    [LIDENT "join_region"; ":"; x = INT -> `NeighborhoodBranches (int_of_string x) ]
                ];
            iterate_options:
                [
                    [LIDENT "model"; ":";  x = model_iter -> x] |
                    [LIDENT "branch"; ":"; x = branch_iter -> x]
                ];
            iteration_method:
                [
                    [ LIDENT "optimize"; ":";
                        left_parenthesis;
                            a = LIST1 [x = iterate_options -> x] SEP ",";
                        right_parenthesis -> a ]
                ];
            std_search_argument:
                [   
                    [ LIDENT "target_cost"; ":"; x = integer_or_float -> `Target
                    (float_of_string x) ] |
                    [ LIDENT "memory"; ":"; x = memory -> `MaxRam x ] |
                    [ LIDENT "hits"; ":"; x = INT -> `MinHits (int_of_string x) ] |
                    [ LIDENT "max_time"; ":"; x = time -> `MaxTime (float_of_int x)
                    ] |
                    [ LIDENT "visited"; x = OPT string_arg -> `Visited x ] |
                    [ LIDENT "min_time"; ":"; x = time -> 
                        `MinTime (float_of_int x) ] |
                    [ LIDENT "constraint"; ":"; x = STRING -> `ConstraintFile x ]
                ];
            search:
                [
                    [ LIDENT "search"; left_parenthesis; a = LIST0 [ x =
                        std_search_argument -> x ] SEP ",";  right_parenthesis ->
                        `StandardSearch a] |
                    [ LIDENT "_search"; left_parenthesis; a = LIST0 [x =
                        search_argument -> x]  SEP
                    ","; right_parenthesis -> `Search a ]
                ];
            fuse:
                [
                    [ LIDENT "fuse"; left_parenthesis; a = LIST0 [x = fuse_argument
                    -> x] SEP ","; 
                        right_parenthesis -> (`Fuse a) ]
                ];
            fuse_argument:
                [
                    [ LIDENT "keep"; ":"; i = INT -> `Keep (int_of_string i) ]
                |   [ LIDENT "iterations"; ":"; i = INT -> `Iterations (int_of_string i) ]
                |   [ LIDENT "replace"; ":"; r = fuseareplace -> `Replace r]
                |   [ x = swap -> (x :> fusea) ]
                |   [ LIDENT "weighting"; ":"; w = fuseaweighting -> `Weighting w]
                |   [ LIDENT "clades"; ":"; cfrom = INT; cto = OPT fusea_cto ->
                          `Clades (int_of_string cfrom, cto)]
                |   [ a = iteration_method -> `IterationF a]
                ];
            fuseareplace:
                [ [ LIDENT "better" -> `Better ]
                | [ LIDENT "best" -> `Best ] ];
            fuseaweighting:
                [ [ LIDENT "uniform" -> `Uniform ] ];
            fusea_cto:
                [ [ "-"; cto = INT -> int_of_string cto ] ];
            calculate_support:
                [
                    [ LIDENT "calculate_support"; left_parenthesis; a = LIST0
                    [x = support_argument -> x]  SEP ","; 
                    right_parenthesis -> `Support a ]
                ];
            (* Reading a file *)

            prealigned_costs:
                [
                    [ ","; LIDENT "tcm"; ":"; left_parenthesis; x = tcm_arguments; right_parenthesis ->
                        let x = match x with
                            | `Tcm (f,l) -> (`Assign_Transformation_Cost_Matrix (`Local f,l))
                            | `Gap (a,b) -> (`Create_Transformation_Cost_Matrix (a, b))
                            | _          -> assert false
                        in
                        x
                    ]
                ];
            prealigned_level:
                [ 
                    [ ","; LIDENT "level"; ":"; x = level_and_tiebreaker -> `Level x ]
                ];
            read_argument:
                [
                    [ LIDENT "annotated"; ":"; left_parenthesis;
                        a = LIST1 [x = otherfiles -> x] SEP ","; right_parenthesis ->
                            ((`AnnotatedFiles a) :> Methods.input) ] |
                    [ LIDENT "prealigned"; ":"; left_parenthesis; a = otherfiles_pre;
                        b = OPT prealigned_costs; c = OPT prealigned_level;
                        right_parenthesis -> match a with
                            | `GeneralAlphabetSeq (a,_,_) ->
                                begin match b with
                                    | Some (`Assign_Transformation_Cost_Matrix (f,_)) ->
                                        let oth = match c with
                                            | None   -> [`Prealigned]
                                            | Some c -> [`Prealigned;c]
                                        in
                                        ((`GeneralAlphabetSeq (a,f,oth)) :> Methods.input)
                                    | Some ( `Create_Transformation_Cost_Matrix _)
                                    | None -> failwith "I@ require@ an@ explicit@ cost@ matrix@ for@ custom@ alphabet@ characters."
                                end
                            | `Aminoacids (f,o) ->
                                let x : Methods.read_option_t list = match c with
                                    | None   -> o
                                    | Some c -> c::o
                                in
                                ((`Aminoacids (f,x)) :> Methods.input)
                            | _ ->
                                let b = match b with
                                    | Some b -> b
                                    | None   -> `Create_Transformation_Cost_Matrix (1,1)
                                in
                                `Prealigned (a, b, 0) ] |
                    [ x = otherfiles -> (x :> Methods.input) ]
                ];
        otherfiles_pre: (** subset of characters that are prealigned **)
            [   [ f = STRING -> `AutoDetect [`Local f] ] |
                [LIDENT "nucleotide"; ":"; left_parenthesis;
                    a = LIST1 [x = STRING -> x] SEP ","; right_parenthesis ->
                        `Nucleotides (to_local a) ] |
                [ LIDENT "nucleotide"; ":"; left_parenthesis;
                    a = LIST1 [x = STRING -> x] SEP ","; right_parenthesis ->
                        `Nucleotides (to_local a) ] |
                [ LIDENT "custom_alphabet"; ":"; left_parenthesis;
                    a = LIST0 [x = STRING -> `Local x] SEP ","; right_parenthesis ->
                        `GeneralAlphabetSeq (a,`Local "",[]) ] |
                [LIDENT "aminoacid"; ":"; left_parenthesis;
                    a = LIST1 [x = STRING -> x] SEP ","; right_parenthesis ->
                        `Aminoacids (to_local a,[`Prealigned]) ] |
                [ LIDENT "aminoacids"; ":"; left_parenthesis;
                    a = LIST1 [x = STRING -> x] SEP ","; right_parenthesis ->
                        `Aminoacids (to_local a,[`Prealigned]) ]
                ];
        otherfiles:
            [   [ f = STRING -> `AutoDetect [`Local f] ] |
                [ LIDENT "partitioned"; ":"; left_parenthesis;
                    a = LIST1 [x = STRING -> x] SEP ","; right_parenthesis ->
                        `PartitionedFile (to_local a) ] |
                [ LIDENT "nucleotides"; ":"; left_parenthesis;
                    a = LIST1 [x = STRING -> x] SEP ","; right_parenthesis ->
                        `Nucleotides (to_local a) ] |
                [ LIDENT "nucleotide"; ":"; left_parenthesis;
                    a = LIST1 [x = STRING -> x] SEP ","; right_parenthesis ->
                        `Nucleotides (to_local a) ] |
                [ LIDENT "chromosome"; ":"; left_parenthesis;
                    a = LIST1 [x = STRING -> x] SEP ","; right_parenthesis ->
                        `Chromosome (to_local a) ] |
                [ LIDENT "genome"; ":"; left_parenthesis;
                    a = LIST1 [x = STRING -> x] SEP ","; right_parenthesis ->
                        `Genome (to_local a) ] |
                [ LIDENT "aminoacid"; ":"; left_parenthesis; 
                    a = LIST1 [x = read_optiona -> x] SEP ","; right_parenthesis ->
                        let files,options : read_option_t list * read_option_t list = 
                            List.partition (function `InputFile _ -> true | _  -> false) a
                        in
                        let files = 
                            List.map (function `InputFile x -> x | _ -> assert false) files
                        in
                        (`Aminoacids (files,options) :> Methods.simple_input)] |
                [ LIDENT "aminoacids"; ":"; left_parenthesis; 
                    a = LIST1 [x = read_optiona -> x] SEP ","; right_parenthesis ->
                        let files,options : read_option_t list * read_option_t list = 
                            List.partition
                                (function `InputFile _ -> true | _ -> false) a
                        in
                        let files =
                            List.map (function `InputFile x -> x | _ -> assert false) files
                        in
                        (`Aminoacids (files,options) :> Methods.simple_input) ] |
                [ LIDENT "custom_alphabet"; ":"; left_parenthesis;
                    a = LIST0 [x = read_optiona -> x] SEP ","; right_parenthesis ->
                        let files,cmo,ros = organize_read_options a in
                        begin match cmo with
                            | Some x -> `GeneralAlphabetSeq (files, x, ros)
                            | None   -> failwith "I require a cost matrix in the read command"
                        end ] |
                [ LIDENT "breakinv"; ":"; left_parenthesis;
                    a = LIST0 [x = read_optiona -> x] SEP ","; right_parenthesis ->
                        let files,cmo,ros = organize_read_options a in
                        begin match cmo with
                            | Some x -> `Breakinv (files, x, ros)
                            | None   -> failwith "I require a cost matrix in the read command"
                        end ] |
                [ LIDENT "complex"; left_parenthesis;
                    a = LIST1 [x = STRING -> x] SEP ","; right_parenthesis ->
                        `ComplexTerminals (to_local a) ] 
            ];
        read_optiona:
            [
                [ x = STRING -> `InputFile (`Local x) ] |
                [LIDENT "init3D"; ":"; init3D = boolean -> `Init3D init3D] |
                [LIDENT "orientation"; ":"; ori = boolean -> `Orientation ori] |
                [LIDENT "tcm"; ":"; left_parenthesis; cm = STRING; right_parenthesis -> `CostMatrix (`Local cm) ] |
                [LIDENT "level"; ":"; x = level_and_tiebreaker -> `Level x ] |
                [LIDENT "tie_breaker"; ":"; x = keep_method -> `Tie_Breaker x] |
                [LIDENT "gap_opening"; ":"; x = INT -> `Affine (int_of_string x) ]
            ];
        tree_information_list:
            [   
                [ ":"; "("; x = LIST0 [x = tree_information -> x] SEP ","; ")" -> x ]
            ];
        report_branch :
            [
                [LIDENT "min"    -> `Final ] |
                [LIDENT "max"    -> `Max ] |
                [LIDENT "single" -> `Single ] |
                [LIDENT "false"  -> `None ]
            ];
        report_branch_opt :
            [   [ ":"; x = report_branch -> x]  ];
        tree_information:
            [
                [ LIDENT "_cost" -> `Cost ] |
                [ LIDENT "hennig" -> `HennigStyle ] |
                [ LIDENT "nexus" -> `NexusStyle ] |
                [ LIDENT "total" -> `Total ] |
                [ LIDENT "newick" -> `Newick ] |
                [ LIDENT "branches"; p = OPT report_branch_opt ->
                    let x :> Methods.report_branch option = match p with
                        | Some `None   -> None
                        | Some `Final  -> Some `Final
                        | Some `Max    -> Some `Max
                        | Some `Single -> Some `Single
                        | None         -> Some `Single
                    in
                    `Branches x ] |
                [ LIDENT "margin"; ":"; m = INT -> `Margin (int_of_string m) ] |
                [ LIDENT "nomargin" -> `Margin (1000000010 - 1) ] |
                [ x = old_identifiers -> `Chars x ] |
                [ x = collapse -> (x :> Methods.information_contained)  ]

            ];
        collapse:
            [
                [ LIDENT "collapse"; y = OPT report_branch_opt ->
                    let x :> Methods.report_branch option = match y with
                        | Some `None   -> None
                        | Some `Final  -> Some `Final
                        | Some `Max    -> Some `Max
                        | Some `Single -> Some `Single
                        | None         -> Some `Single
                    in
                    `Collapse x]
            ];
        optional_collapse:
            [ 
                [ ":"; x = collapse -> x ]
            ];
        optional_boolean:
            [
                [ ":"; x = boolean -> x ]
            ];
        boolean: 
            [
                [ LIDENT "true" -> true ] |   
                [ LIDENT "false" -> false ]    
            ];
        (* Building a tree *)
        join_method:
            [
                [ LIDENT "sectorial"; x = OPT integer -> 
                    ((`UnionBased x) : Methods.tabu_join_strategy)  ] |
                [ LIDENT "all"; x = OPT integer -> `AllBased x ] |
                [ LIDENT "constraint"; ":"; x = INT ->
                    `Partition [`MaxDepth (int_of_string x)] ] |
                [ LIDENT "constraint"; ":"; left_parenthesis; 
                    x = LIST1 [x = constraint_options -> x] SEP ","; right_parenthesis ->
                    `Partition x ] |
                [ LIDENT "constraint" -> `Partition [] ]
            ];
        build_argument:
            [
                [ x = threshold_and_trees -> (x :> builda) ] |
                [ x = build_method -> (x :> builda) ] |
                [ x = join_method -> (x :> builda) ] |
                [ x = keep_method -> (x :> builda) ] |
                [ a = iteration_method -> `IterationB a] |
                [ x = cost_calculation -> (x :> builda) ] |
                [ LIDENT "lookahead"; ":"; x = INT -> 
                    `Lookahead (int_of_string x) ]
            ];
        threshold_and_trees:
            [
                [ LIDENT "threshold"; ":"; x = integer_or_float -> `Threshold (float_of_string x) ] |
                [ LIDENT "trees"; ":"; x = INT -> `Trees (int_of_string x) ] |
                [ x = INT -> `Trees (int_of_string x) ]
            ];
        build_method:
            [
                [ x = STRING -> `Prebuilt (`Local x) ] |
                [ LIDENT "of_file"; ":"; x = STRING -> `Prebuilt (`Local x) ] |
                [ LIDENT "randomized" -> `Random ] |
                [ LIDENT "random" -> `RandomTree ] |
                [ LIDENT "as_is" -> `Ordered ] |
                [ LIDENT "branch_and_bound"; x = OPT optional_integer_or_float -> 
                    let thresh = 
                        match x with
                        | None ->  None
                        | Some x -> Some (float_of_string x)
                    in
                    `Branch_and_Bound thresh ] |
                [ LIDENT "constraint"; x = OPT optional_string -> `Constraint x ]
                        |
                [ LIDENT "nj" -> `Nj ] | 
                [ LIDENT "_mst" -> `Mst ] |
                [ LIDENT "_distances" -> `DistancesRnd ]
            ];
        keep_method:
            [ 
                [ LIDENT "last" -> `Last ] |
                [ LIDENT "first" -> `First ] | 
                [ LIDENT "at_random" -> `Keep_Random ] 
            ];
        cost_calculation:
            [
                [ x = transform -> (x :> transform) ]
            ];
        (* Swaping *)
        search_argument:
            [ 
                [ LIDENT "build"; x = OPT optional_boolean -> 
                    (match x with
                    | None -> `Build true
                    | Some x -> `Build x) ] |
                [ LIDENT "transform"; x = OPT optional_boolean -> 
                    match x with
                    | None -> `Transform false
                    | Some x -> `Transform x ]
            ];
        swap_argument:
            [
                [ x = sample_method -> (x :> swapa) ] |
                [ x = swap_method -> (x :> swapa) ] |
                [ x = keep_method -> (x :> swapa) ] |
                [ x = threshold_and_trees -> (x :> swapa) ] |
                [ x = transform -> (x :> swapa) ] |
                [ LIDENT "forest"; a = OPT optional_integer_or_float -> 
                    match a with
                    | None -> `Forest 0.
                    | Some a -> `Forest (float_of_string a) ] |   
                [ a = trajectory_method -> a ] |
                [ a = break_method -> a ] |
                [ a = reroot_method -> a ] |
                [ a = join_method -> (a :> swapa) ] |
                [ a = iteration_method -> `IterationS a]
            ];
        trajectory_method:
            [
                [ LIDENT "annealing"; ":"; left_parenthesis; 
                  coeff = integer_or_float;
                  ","; exp = integer_or_float; right_parenthesis ->
                      `Annealing (float_of_string coeff,
                                  float_of_string exp)] |
                [ LIDENT "drifting"; ":"; left_parenthesis;
                    equalprob = integer_or_float; ","; 
                    worstfactor = integer_or_float; right_parenthesis ->
                        `PoyDrifting (float_of_string equalprob, float_of_string worstfactor) ] |
                [ LIDENT "current_neighborhood"; f = OPT string_arg -> 
                    `AllAround f ] |
                [ LIDENT "around" -> `AllThenChoose ]
            ];
        sample_method:
            [
                [ LIDENT "_maxtrees"; x = integer ->
                    `MaxTreesEvaluated x ] |
                [ LIDENT "timeout"; ":"; x = integer_or_float -> 
                    `TimeOut (`Fixed (float_of_string x)) ] |
                [ LIDENT "timedprint"; ":"; left_parenthesis; x = integer_or_float; ","; 
                    y = STRING; right_parenthesis -> 
                        `TimedPrint (float_of_string x, Some y) ] |
                [ LIDENT "visited"; x = OPT string_arg ->
                    `AllVisited x ] |
                [ LIDENT "trajectory"; x = OPT string_arg -> 
                        `PrintTrajectory x ] |
                [ LIDENT "recover" -> `KeepBestTrees ] |
                [ LIDENT "_unionstats"; ":"; left_parenthesis; 
                    x = STRING; ","; y = INT; right_parenthesis ->
                    `UnionStats (Some x, int_of_string y) ] |
                [ LIDENT "_rootuniondistr"; ":"; x = STRING -> 
                    `RootUnionDistr (Some x) ] |
                [ LIDENT "_attemptsdistr"; ":"; x = STRING ->
                    `AttemptsDistr (Some x) ] |
                [ LIDENT "_breakvsjoin"; x = OPT string_arg ->
                    `BreakVsJoin x ] |
                [ LIDENT "_likelihood"; x = OPT string_arg ->
                    `Likelihood x ] |
                [ LIDENT "_model"; x = OPT string_arg ->
                    `LikelihoodModel x ]
            ];
        swap_method:
            [
                [ LIDENT "spr"; y = OPT single_option -> match y with
                    | None -> `ChainNeighborhoods `Spr
                    | Some _ -> `SingleNeighborhood `Spr ] |
                [ LIDENT "tbr"; y = OPT single_option -> match y with
                    | None -> `ChainNeighborhoods `Tbr
                    | Some _ -> `SingleNeighborhood `Tbr ] |
                [ LIDENT "alternate" -> `Alternate (`Spr, `Tbr) ]
            ];
        single_option:
            [
                [ ":"; LIDENT "once" -> 1 ]
            ];
        break_method:
            [
                [ LIDENT "randomized" -> `Randomized ] |
                [ LIDENT "distance"; x = OPT optional_boolean -> match x with
                    | None -> `DistanceSorted false
                    | Some x -> `DistanceSorted x] |
                [ LIDENT "once" -> `OnlyOnce ]
            ];
        constraint_options:
            [
                [ x = INT -> `MaxDepth (int_of_string x) ] |
                [ LIDENT "depth"; ":"; x = INT -> `MaxDepth (int_of_string x) ] |
                [ x = STRING -> `ConstraintFile (`Local x) ] |
                [ LIDENT "file"; ":"; x = STRING -> `ConstraintFile (`Local x) ]
            ];
        reroot_method:
            [
                [ LIDENT "bfs"; x = OPT integer -> `Bfs x ]
            ];
        integer:
            [
                [ ":"; x = INT -> int_of_string x ]
            ];
        string_arg:
            [
                [ ":"; x = STRING -> x ]
            ];
        report_type:
            [
                [ LIDENT "statename_only" -> `StateOnly ]
            ];
        opt_report_type:
            [
                [":"; x = report_type -> x]
            ];
        (* Support values *)
        support_argument:
            [ [ x = build -> (x :> supporta) ] |
                [ x = swap -> (x :> supporta) ] |
                [ x = support_method -> (x :> supporta) ]
            ];
        support_method:
            [
                [ LIDENT "bremer" -> `Bremer ] |
                [ LIDENT "jackknife"; x = OPT list_of_jackknifea -> match x with
                    | None -> `Jackknife []
                    | Some v -> `Jackknife v ] |
                [ LIDENT "bootstrap"; x = OPT integer -> `Bootstrap x ]
            ];
        opt_support_names:
            [ [ ":"; y = support_names -> y ] ];
        support_names:
            [
                [ LIDENT "bremer"; ":"; x = STRING -> `Bremer (Some [(`Local x)]) ] |
                [ LIDENT "bremer"; ":"; left_parenthesis; x = LIST1 [ y = STRING -> y] SEP ",";
                    right_parenthesis ->
                        `Bremer (Some (List.map (fun x -> `Local x) x))] |
                [ LIDENT "bremer" -> `Bremer None ] |
                [ LIDENT "jackknife"; y = OPT [ x = summary_class -> x ] -> match y with
                    | None -> `Jackknife `Individual
                    | Some y -> `Jackknife y] |
                [ LIDENT "bootstrap"; y = OPT [ x = summary_class -> x ] -> match y with
                    | None -> `Bootstrap `Individual
                    | Some y -> `Bootstrap y]
            ];
        summary_class:
            [
                [ ":"; LIDENT "individual" -> `Individual ] |
                [ ":"; LIDENT "consensus" -> `Consensus ] |
                [ ":"; x = STRING -> `InputFile x ]
            ];
        list_of_jackknifea:
            [
                [ ":"; left_parenthesis; x = LIST0 [x = jackknifea -> x] SEP ","; right_parenthesis -> x ]
            ];
        jackknifea:
            [
                [ LIDENT "remove"; ":"; x = integer_or_float -> `Select (float_of_string x) ] |
                [ LIDENT "resample"; ":"; x = INT -> `Resample (int_of_string x) ]
            ];
        (* Shared items *)
        left_parenthesis: [ [ "(" ] ];
        right_parenthesis: [ [ ")" ] ];

        (* How we identify characters; in files *)
        identifiers:
            [
                [ LIDENT "not" ; LIDENT "files"; ":"; 
                    left_parenthesis; x = LIST0 [x = STRING -> x] SEP ","; right_parenthesis ->
                        `Files (false, x) ] |
                [ x = old_identifiers -> (x :> identifiers) ] |
                [ LIDENT "files"; ":"; left_parenthesis; x = LIST0 [x = STRING
                -> x] SEP ","; 
                    right_parenthesis -> `Files (true, x) ] 
            ];
        old_identifiers:
            [
                [ LIDENT "all" -> `All ] |
                [ LIDENT "not"; LIDENT "names"; ":"; left_parenthesis;
                    x = LIST0 [x = STRING -> x] SEP ","; right_parenthesis ->
                        `Names (false, x) ] |
                [ LIDENT "range"; ":"; left_parenthesis;
                    file = STRING; ","; x = INT; ","; y = INT; right_parenthesis ->
                        `Range (true, file, int_of_string x, int_of_string y) ] |
                [ LIDENT "not"; LIDENT "range"; ":"; left_parenthesis;
                    file = STRING; ","; x = INT; ","; y = INT; right_parenthesis ->
                        `Range (false, file, int_of_string x, int_of_string y) ] |
                [ LIDENT "not"; LIDENT "sets"; ":"; left_parenthesis;
                    x = LIST0 [x = STRING -> x] SEP ","; right_parenthesis ->
                        `CharSet (false, x) ] |
                [ LIDENT "not"; LIDENT "codes"; ":"; left_parenthesis; 
                    x = LIST0 [x = INT -> x] SEP ","; right_parenthesis ->
                        `Some (false, List.map int_of_string x) ] |
                [ LIDENT "names"; ":"; left_parenthesis; 
                    x = LIST0 [x = STRING -> x] SEP ","; right_parenthesis ->
                        `Names (true, x) ] |
                [ LIDENT "sets"; ":"; left_parenthesis;
                    x = LIST0 [x = STRING -> x] SEP ","; right_parenthesis ->
                        `CharSet (true, x) ] |
                [ LIDENT "codes"; ":"; left_parenthesis;
                    x = LIST0 [x = INT -> x] SEP ","; right_parenthesis ->
                        `Some (true, List.map int_of_string x) ] |
                [ LIDENT "static" -> `AllStatic ] | 
                [ LIDENT "dynamic" -> `AllDynamic ] |
                [ LIDENT "missing"; ":"; x = INT -> 
                    `Missing (true, 100 - int_of_string x) ] |
                [ LIDENT "not"; LIDENT "missing"; ":"; x = INT ->
                    `Missing (false, 100 - int_of_string x) ] |
                [ LIDENT "_random"; ":"; x = integer_or_float -> 
                    `Random (100. -. (float_of_string x)) ]
            ];
        optional_integer_or_float: 
            [
                [ ":"; x = integer_or_float -> x ]
            ];
        neg_integer_or_float:
            [
                [ "-"; x = integer_or_float -> "-" ^ x ] |
                [ x = integer_or_float -> x ]
            ];
        integer_or_float:
            [
                [ x = INT -> x ] |
                [ x = FLOAT -> x ]
            ];
        string_or_ident:
            [
                [ x = STRING -> x ]
            |   [ x = LIDENT -> x ]
            ];
    END
    ;
    expr

let ( --> ) a b = b a

let rec process_commands optimize command = 
    match command with
    | `ChangeWDir dir -> let dir = simplify_directory dir in Sys.chdir dir; [command]
    | `ReadScript files -> 
            read_script_files false (List.map (fun x -> `Filename x) files)
    | x -> [x]

and read_script_files optimize (files : [`Inlined of string | `Filename of string] list)  = 
    let res = 
        List.map
            (fun f -> 
                let where = match f with
                    | `Filename f -> "file@ @{<b>" ^ f ^ "@}"
                    | `Inlined f -> "inlined@ script"
                in
                try match f with
                    | `Inlined f ->
                        let comm = of_string false f in
                        `Echo ("Running inlined script", `Information) :: comm
                    | `Filename f ->
                        let f = simplify_directory f in
                        let comm = of_file false f in
                        `Echo ("Running file " ^ f, `Information) :: comm
                with 
                    | Loc.Exc_located (a, Stream.Error err) ->
                        let is_unknown = "illegal begin of expr" = err in
                        let msg =
                            Printf.sprintf
                                ("@[<v 4>@[Command@ error@ in@ %s@}@ line@ "^^
                                 "@{<b>%d@}@ between@ characters@ @{<b>%d@}@ "^^
                                 "and @{<b>%d@} :@]@,@[%s@]@]\n")
                                (where) (Loc.start_line a)
                                ((Loc.start_off a) - (Loc.start_bol a))
                                ((Loc.stop_off a) - (Loc.stop_bol a))
                                (if is_unknown then "Unknown command" else err)
                        in
                        Status.user_message Status.Error msg;
                        failwith "Script execution stopped"
                    | err ->
                        Status.user_message Status.Error 
                            ("Error@ while@ processing@ script@ " ^ where);
                        raise err)
            files 
    in
    do_analysis optimize (List.flatten res)

and do_analysis optimize res =
    if debug then Printf.printf "do analysis,optimize=%b\n%!" optimize;
    if optimize then Analyzer.analyze res
    else res

and simplify_directory dir = 
    Str.global_replace (Str.regexp "\\\\ ") " " dir 

and of_parsed optimize lst =
    let cur_directory = Sys.getcwd () in
    let res =
        lst --> transform_all_commands
            --> List.map (process_commands false)
            --> List.flatten
            --> do_analysis optimize
    in
    let cur_directory = simplify_directory cur_directory in
    Sys.chdir cur_directory;
    res

and of_stream optimize str =
    let cur_directory = Sys.getcwd () in
    let expr = create_expr () in
    let res = 
        str (* --> add_command_to_console_script *)
            --> Gram.parse expr (Loc.mk "<stream>")
            --> transform_all_commands
    in
    let res =
        res --> List.map (process_commands false)
            --> List.flatten
            --> do_analysis optimize
    in
    let cur_directory = simplify_directory cur_directory in
    Sys.chdir cur_directory;
    res


and of_channel optimize ch = 
    of_stream optimize (Stream.of_channel ch)

and of_file optimize f =
    if debug then Printf.printf "poyCommand.of_file %s,optimize=%b\n%!" f optimize;
    let ch = open_in f in
    let r = of_channel optimize ch in
    if debug then Printf.printf "end of of_file %s\n%!" f;
    close_in ch;
    r

and of_string optimize str =
    if StatusCommon.redirect_information () then
        let file = StatusCommon.redirect_filename () in
        Status.user_message (Status.Output ((Some file), false, [])) (str ^ "@\n")
    else ();
    of_stream optimize (Stream.of_string str)
