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

(** A dynamic character set implementation. 
* The dynamic character set allows rearrangements *)

exception Illegal_Arguments
let () = SadmanOutput.register "DynamicCS" "$Revision: 3106 $"

module IntMap = All_sets.IntegerMap
module IntSet = All_sets.Integers

exception No_Union

(** A dynamic character type. 'a can be SeqCS, ChromCS, BreakCS, GenomeCS ...*)
type t = 
    | SeqCS  of SeqCS.t
    | MlCS of MlDynamicCS.t
    | BreakinvCS of  BreakinvCS.t 
    | ChromCS of ChromCS.t 
    | AnnchromCS of AnnchromCS.t
    | GenomeCS of GenomeCS.t

type u =
    | U_SeqCS of SeqCS.Union.u
    | U_Others


let failwith_todo f_name = 
    failwith ("Todo: " ^ f_name ^ " dynamicCS")

(** [alpha a] returns the alphabet of dynamic character set [a] *)
let alpha (a : t) = match a with
    | MlCS a -> MlDynamicCS.alph a
    | SeqCS a -> a.SeqCS.alph
    | ChromCS a -> a.ChromCS.alph
    | GenomeCS a -> a.GenomeCS.alph
    | BreakinvCS a -> a.BreakinvCS.alph
    | AnnchromCS a -> a.AnnchromCS.alph

(** [total_cost a] returns the total cost to create
*  dynamic character set [a]*)
let total_cost (a : t) = match a with
    | MlCS a -> MlDynamicCS.total_cost a
    | SeqCS a -> a.SeqCS.total_cost
    | ChromCS a -> a.ChromCS.total_cost
    | GenomeCS a -> a.GenomeCS.total_cost
    | BreakinvCS a -> a.BreakinvCS.total_cost
    | AnnchromCS a -> a.AnnchromCS.total_cost

(** [get_subtree_cost a] returns the cost of subtree *)
let get_subtree_cost (a : t) = match a with
    | MlCS a -> MlDynamicCS.total_cost a (*to do, nic*)
    | SeqCS a -> a.SeqCS.subtree_cost
    | ChromCS a -> a.ChromCS.subtree_cost
    | GenomeCS a -> a.GenomeCS.subtree_cost
    | BreakinvCS a -> a.BreakinvCS.subtree_cost
    | AnnchromCS a -> a.AnnchromCS.subtree_cost


(** [total_recost a] returns the total recost to create
*  dynamic character set [a]*)
let total_recost (a : t) = match a with
    | MlCS _ | SeqCS _ -> 0.
    | ChromCS a -> a.ChromCS.total_recost
    | GenomeCS a -> a.GenomeCS.total_recost
    | BreakinvCS a -> a.BreakinvCS.total_recost
    | AnnchromCS a -> a.AnnchromCS.total_recost

(** [subtree_recost a] returns the total recost of
* subtree whose root contains dynamic character set [a] *)
let subtree_recost (a : t) = match a with
    | MlCS _ | SeqCS _ -> 0.0
    | ChromCS a -> a.ChromCS.subtree_recost
    | GenomeCS a -> a.GenomeCS.subtree_recost
    | BreakinvCS a -> a.BreakinvCS.subtree_recost
    | AnnchromCS a -> a.AnnchromCS.subtree_recost

(** [c2_full a] returns the two dimentional cost matrix
* of dynamic character set [a] *)
let c2_full (a : t) = match a with 
    | MlCS a -> MlDynamicCS.get_cm a
    | SeqCS a -> a.SeqCS.heuristic.SeqCS.c2_full
    | ChromCS a -> a.ChromCS.c2_full
    | GenomeCS a -> a.GenomeCS.c2_full
    | BreakinvCS a -> a.BreakinvCS.c2_full
    | AnnchromCS a -> a.AnnchromCS.c2_full


(** [lk_model a] returns the likelihood model for the dynamic likelihood
    character, or raise a Not_found for other characters *)
let lk_model (a : t) = match a with 
    | MlCS a       -> MlDynamicCS.model a
    | SeqCS _ 
    | ChromCS _
    | GenomeCS _
    | BreakinvCS _
    | AnnchromCS _ -> raise Not_found

(** [chrom_pam a] returns the user-defined chromosome parameters
* of dynamic character set [a] *)
let chrom_pam (a : t) = match a with
    | ChromCS a -> a.ChromCS.chrom_pam
    | GenomeCS a -> a.GenomeCS.chrom_pam
    | AnnchromCS a -> a.AnnchromCS.annchrom_pam
    | BreakinvCS a -> a.BreakinvCS.breakinv_pam
    | _ -> Data.dyna_pam_default 

(** [state a] returns the character type of
* dynamic character set [a] *)
let state (a : t) : Data.dyna_state_t = match a with
    | MlCS a -> `Ml
    | SeqCS a -> `Seq
    | ChromCS a -> `Chromosome
    | GenomeCS a -> `Genome
    | BreakinvCS a -> `Breakinv
    | AnnchromCS a -> `Annotated

(** [code a] returns the code of dynamic character set [a] *)
let code (a : t) = match a with
    | MlCS a -> MlDynamicCS.code a
    | SeqCS a -> a.SeqCS.code
    | ChromCS a -> a.ChromCS.code
    | GenomeCS a -> a.GenomeCS.code
    | BreakinvCS a -> a.BreakinvCS.code
    | AnnchromCS a -> a.AnnchromCS.code

let codes = function
    | MlCS a       -> a.MlDynamicCS.codes
    | SeqCS a      -> a.SeqCS.codes
    | ChromCS a    ->
        Array.of_list (IntMap.fold (fun k _ a -> k::a) a.ChromCS.meds [])
    | GenomeCS a   ->
        Array.of_list (IntMap.fold (fun k _ a -> k::a) a.GenomeCS.meds [])
    | BreakinvCS a ->
        Array.of_list (IntMap.fold (fun k _ a -> k::a) a.BreakinvCS.meds [])
    | AnnchromCS a ->
        Array.of_list (IntMap.fold (fun k _ a -> k::a) a.AnnchromCS.meds [])

let print dyn = match dyn with
    | ChromCS ch -> ChromCS.print ch
    | AnnchromCS ch -> AnnchromCS.print ch 
    | SeqCS ch -> SeqCS.print ch
    | BreakinvCS ch -> BreakinvCS.print ch  
    | _ -> print_endline "Do not print non-chromosome characters"

(** [copy_chrom_map s_ch d_ch] copies the choromosome map
* of dynamic character set [s_ch] to [d_ch] *)
let copy_chrom_map s_ch d_ch = match s_ch, d_ch with
    | ChromCS s_ch, ChromCS d_ch ->
          ChromCS (ChromCS.copy_chrom_map s_ch d_ch)
    | AnnchromCS s_ch, AnnchromCS d_ch ->
          AnnchromCS (AnnchromCS.copy_chrom_map s_ch d_ch)
    | GenomeCS s_ch, GenomeCS d_ch ->
          GenomeCS (GenomeCS.copy_chrom_map s_ch d_ch) 
    | _, _ -> d_ch

    
(** [leaf_sequences a] turns dynamic character set [a] 
* into a set of chromosome arrays *)
let leaf_sequences (a : t) = match a with
    | MlCS a  -> MlDynamicCS.leaf_sequences a
    | SeqCS a ->
            let map = ref IntMap.empty in
            for i = (SeqCS.cardinal a) - 1 downto 0 do
                map :=
                    IntMap.add a.SeqCS.codes.(i)
                        (match a.SeqCS.characters.(i) with
                        | SeqCS.Partitioned x ->
                            Array.map
                                (function
                                    | SeqCS.PartitionedDOS.Last x ->
                                        `Last x.SeqCS.DOS.sequence
                                    | SeqCS.PartitionedDOS.DO x  ->
                                        `DO x.SeqCS.DOS.sequence
                                    | SeqCS.PartitionedDOS.First x ->
                                        `First x.SeqCS.DOS.sequence)
                                x
                        | SeqCS.Heuristic_Selection x -> [|`DO x.SeqCS.DOS.sequence|]
                        | SeqCS.General_Prealigned x  -> [|`DO x.GenNonAdd.seq|])
                    !map
            done;
            !map
    | BreakinvCS a ->
          IntMap.map
              (fun med ->
                   [|`DO (List.hd med.Breakinv.med_ls).BreakinvAli.seq|])
              a.BreakinvCS.meds
    | ChromCS a ->
          IntMap.map
              (fun med  ->
                   [|`DO (List.hd med.Chrom.med_ls).ChromAli.seq|])
              a.ChromCS.meds
    | AnnchromCS a ->
          IntMap.map
              (fun med ->
                   let seqt_arr =
                       (List.hd med.Annchrom.med_ls).AnnchromAli.seq_arr
                   in
                   Array.map (fun seqt -> `DO seqt.AnnchromAli.seq) seqt_arr)
              a.AnnchromCS.meds
    | GenomeCS a ->
          IntMap.map
              (fun med ->
                   let seq_arr = Array.map
                       (fun chromt ->
                            `DO chromt.GenomeAli.seq)
                       (List.hd med.Genome.med_ls).GenomeAli.chrom_arr
                   in
                   seq_arr)
              a.GenomeCS.meds

(** [unions a] returns the union of dynamic character set [a] *)
let unions (a : u) = match a with
    | U_SeqCS (Some a) -> 
            let map = ref IntMap.empty in
            for i = (Array.length a.SeqCS.Union.u_codes) - 1 
            downto 0 do
                map := IntMap.add (a.SeqCS.Union.u_codes.(i))
                    (a.SeqCS.Union.unions.(i)) !map
            done;
            !map
    | _ -> failwith "DynamicCS.unions"

let to_union a = match a with
    | SeqCS a -> U_SeqCS (SeqCS.to_union a)
    | _ -> U_Others

let union a b c = match a, b, c with
    | SeqCS a, U_SeqCS b, U_SeqCS c -> U_SeqCS (SeqCS.Union.union a b c)
    | SeqCS _, _, _ -> failwith "DynamicCS.union"
    | _, U_Others, U_Others -> U_Others
    | _, _, _ -> failwith_todo "union"


let cardinal_union a = match a with
    | U_SeqCS a -> SeqCS.Union.cardinal_union a
    | U_Others -> 0

let poly_saturation x v =
    let polyacc, polylen =
        match x with
        | U_SeqCS x ->
                let card = SeqCS.Union.cardinal_union x in
                let poly = SeqCS.Union.poly_saturation x v in
                (int_of_float (poly *. (float_of_int card))), 
                card
        | U_Others -> 0, 0
    in
    (float_of_int polyacc) /. (float_of_int polylen)

(** [of_array spec genome_arr code taxon num_taxa] 
* creates a dynamic character set from genome array [genome_arr] *)
let of_array spec genome_arr code taxon num_taxa =
    match spec.Data.state with
    | `SeqPrealigned
    | `Ml | `Seq | `CustomAlphabet as meth ->
            let seq_arr = 
                Array.map 
                (fun (genome_data, genome_code) ->
                    let seqs = 
                        Array.map (fun x -> x.Data.seq)
                        genome_data.Data.seq_arr
                    in  
                    (seqs, genome_code)) genome_arr
            in 
            let t = SeqCS.of_array spec seq_arr code taxon in
            begin match meth, spec.Data.lk_model with
                | `SeqPrealigned,_ 
                | `CustomAlphabet,_
                | `Seq,_ -> SeqCS t
                | `Ml,Some m -> MlCS (MlDynamicCS.make (spec.Data.alph) t m)
                | `Ml,None -> assert false
            end
    | `Breakinv ->
            let seq_arr =
                Array.map
                    (fun (genome_data, genome_code) ->
                        assert( 1 = (Array.length genome_data.Data.seq_arr));
                        let first_seq_deli = genome_data.Data.seq_arr.(0) in
                        let first_seq  = first_seq_deli.Data.seq in 
                        let first_deli = first_seq_deli.Data.delimiter in
                        (first_seq, first_deli, genome_code))
                    genome_arr
            in
            let t = BreakinvCS.of_array spec seq_arr code in
            BreakinvCS t
    | `Chromosome ->
            let seq_arr =
                Array.map
                    (fun (genome_data, genome_code) ->
                        assert( 1 = (Array.length genome_data.Data.seq_arr));
                        let first_seqa = genome_data.Data.seq_arr.(0) in
                        let first_seq = first_seqa.Data.seq in 
                        (first_seq, genome_code))
                    genome_arr
            in
            let t = ChromCS.of_array spec seq_arr code taxon num_taxa in
            ChromCS t
    | `Annotated -> 
          let t = AnnchromCS.of_array spec genome_arr code  taxon num_taxa in
          AnnchromCS t 
    | `Genome  ->
            let t = 
                GenomeCS.of_array spec genome_arr code taxon num_taxa 
            in 
            GenomeCS t 


(** [of_list spec genome_arr code taxon num_taxa] 
* creates a dynamic character set from genome list [genome_ls] *)
let of_list spec genome_ls =
    of_array spec (Array.of_list genome_ls) 

(** [median a b] creates the median set between dynamic 
* character sets [a] and [b] *)
let median code a b t1 t2 = match a, b with
    | MlCS a, MlCS b -> MlCS (MlDynamicCS.median code a b t1 t2)
    | SeqCS a, SeqCS b -> SeqCS (SeqCS.median code a b)
    | ChromCS a, ChromCS b -> ChromCS (ChromCS.median2 a b)
    | GenomeCS a, GenomeCS b -> GenomeCS (GenomeCS.median2 a b)
    | BreakinvCS a, BreakinvCS b -> BreakinvCS (BreakinvCS.median2 a b)
    | AnnchromCS a, AnnchromCS b -> AnnchromCS (AnnchromCS.median2 a b)
    | _, _ -> failwith_todo "median"

(** [median_3 p n c1 c2] creates the median among
* three dynamic character sets [p], [c1] and [c2] *)
let median_3 p n c1 c2 = match p, n, c1, c2 with 
    | MlCS p, MlCS n, MlCS c1, MlCS c2 -> 
          MlCS (MlDynamicCS.median_3 p n c1 c2)
    | SeqCS p, SeqCS n, SeqCS c1, SeqCS c2 -> 
          SeqCS (SeqCS.median_3 p n c1 c2)
    | BreakinvCS p, BreakinvCS n, BreakinvCS c1, BreakinvCS c2 -> 
          BreakinvCS (BreakinvCS.median3 p n c1 c2)
    | AnnchromCS p, AnnchromCS n, AnnchromCS c1, AnnchromCS c2 -> 
          AnnchromCS (AnnchromCS.median3 p n c1 c2)
    | ChromCS p, ChromCS n, ChromCS c1, ChromCS c2 -> 
          ChromCS (ChromCS.median3 p n c1 c2)
    | GenomeCS p, GenomeCS n, GenomeCS c1, GenomeCS c2 -> 
          GenomeCS (GenomeCS.median3 p n c1 c2)
    | _, _, _, _ -> failwith_todo "median_3"

(** [cost_for_root root ] return distance between two algned children of
* the root. due to non-zero diagonal this might be different from alignment cost
* of the two children.*)
let extra_cost_for_root (a : t) = match a with
    | MlCS x        -> 0.0
    | SeqCS x       -> SeqCS.extra_cost_for_root x
    | ChromCS x     -> ChromCS.get_extra_cost_for_root x
    | AnnchromCS x  -> AnnchromCS.get_extra_cost_for_root x 
    | GenomeCS x    -> GenomeCS.get_extra_cost_for_root x
    | BreakinvCS x  -> BreakinvCS.get_extra_cost_for_root x


(* Like [distance] but calculates it only 
* if the type of the characters match one of those listed. *)
let distance_of_type t missing_distance a b len =
    let has_t x = List.exists (fun z -> z = x) t in
    let has_seq = has_t `Seq    and has_chrom = has_t `Chrom 
    and has_gen = has_t `Genome and has_break = has_t `Breakinv
    and has_lk = has_t `Ml      and has_ann = has_t `Annchrom in
    match a, b with
    | MlCS a, MlCS b when has_lk -> MlDynamicCS.distance missing_distance a b len
    | SeqCS a, SeqCS b when has_seq -> SeqCS.distance missing_distance a b
    | ChromCS a, ChromCS b when has_chrom -> ChromCS.distance a b  
    | GenomeCS a, GenomeCS b when has_gen -> GenomeCS.distance a b  
    | BreakinvCS a, BreakinvCS b when has_break -> BreakinvCS.distance a b  
    | AnnchromCS a, AnnchromCS b when has_ann -> AnnchromCS.distance a b  
    | _, _ -> 0.0

(** [distance_of_type a b] returns the distance between
* two dynamic character sets [a] and [b] *)
let distance missing_distance a b = match a, b with   
    | MlCS a, MlCS b -> MlDynamicCS.distance missing_distance a b None
    | SeqCS a, SeqCS b -> SeqCS.distance missing_distance a b
    | ChromCS a, ChromCS b -> ChromCS.distance a b  
    | GenomeCS a, GenomeCS b -> GenomeCS.distance a b  
    | BreakinvCS a, BreakinvCS b -> BreakinvCS.distance a b  
    | AnnchromCS a, AnnchromCS b -> AnnchromCS.distance a b  
    | _ , _ -> failwith_todo "distance"  



(** [distance_union a b] returns the union distance between
* two dynamic character sets [a] and [b] *)
let distance_union a b = match a, b with
    | U_SeqCS a, U_SeqCS b -> SeqCS.Union.distance_union a b
    | U_Others, U_Others -> 0.0
    | _, _ -> failwith "DynamicCS.distance_union"


(** [to_string a] returns dynamic character set [a] 
* into the string format *)
let to_string a = match a with 
    | MlCS a -> MlDynamicCS.to_string a
    | SeqCS a -> SeqCS.to_string a
    | BreakinvCS a -> BreakinvCS.to_string a
    | ChromCS a -> ChromCS.to_string a
    | GenomeCS a -> GenomeCS.to_string a
    | AnnchromCS a -> AnnchromCS.to_string a


let name_string a = match a with 
    | MlCS t -> MlDynamicCS.name_string t
    | SeqCS _ -> "Sequence"
    | BreakinvCS _ -> "Break Inversion"
    | ChromCS _ -> "Chromosome"
    | GenomeCS _ -> "Genome"
    | AnnchromCS _ -> "Annotated Chromosome"

(* [dist_2 delta n a b] calculates the cost of joining 
* the node containing  chromosome character set [n]
* between two nodes containing [a] and [b]. 
* [a] must be the parent (ancestor) of [b] *)
let dist_2 delta n a b =
    match n, a, b with 
    | MlCS n, MlCS a, MlCS b -> MlDynamicCS.dist_2 delta n a b
    | SeqCS n, SeqCS a, SeqCS b -> SeqCS.dist_2 delta n a b
    | ChromCS n, ChromCS a, ChromCS b -> ChromCS.dist_2 n a b
    | GenomeCS n, GenomeCS a, GenomeCS b -> GenomeCS.dist_2 n a b
    | BreakinvCS n, BreakinvCS a, BreakinvCS b -> BreakinvCS.dist_2 n a b
    | AnnchromCS n, AnnchromCS a, AnnchromCS b -> AnnchromCS.dist_2 n a b
    | _, _, _ -> failwith_todo "dist_2"



(** [f_codes s c] returns a dynamic character subset 
* of  dynamic character set [s] whose codes are 
* also in  chromosome character set [c] *)
let f_codes s c = 
    match s with 
    | MlCS s -> MlCS (MlDynamicCS.f_codes s c)
    | SeqCS s -> SeqCS (SeqCS.f_codes s c)
    | ChromCS s -> ChromCS (ChromCS.f_codes s c)
    | GenomeCS s -> GenomeCS (GenomeCS.f_codes s c)
    | BreakinvCS s -> BreakinvCS (BreakinvCS.f_codes s c)
    | AnnchromCS s -> AnnchromCS (AnnchromCS.f_codes s c)



(** [f_codes_comp s c] returns a dynamic character subset 
* of  dynamic character set [s] whose codes are NOT
* also in  chromosome character set [c] *)
let f_codes_comp s c = 
    match s with 
    | MlCS s -> MlCS (MlDynamicCS.f_codes_comp s c)
    | SeqCS s -> SeqCS (SeqCS.f_codes_comp s c)
    | ChromCS s -> ChromCS (ChromCS.f_codes_comp s c)
    | GenomeCS s -> GenomeCS (GenomeCS.f_codes_comp s c)
    | BreakinvCS s -> BreakinvCS (BreakinvCS.f_codes_comp s c)
    | AnnchromCS s -> AnnchromCS (AnnchromCS.f_codes_comp s c)


(** [compare_data a b] returns 1 if dynamic character 
* set [a] is the same dynamic character set [b], 
* otherwise (-1) or 1 *)
let compare_data a b =
    match a, b with 
    | SeqCS a, SeqCS b -> SeqCS.compare_data a b
    | ChromCS a, ChromCS b -> ChromCS.compare_data a b
    | GenomeCS a, GenomeCS b -> GenomeCS.compare_data a b
    | BreakinvCS a, BreakinvCS b -> BreakinvCS.compare_data a b
    | AnnchromCS a, AnnchromCS b -> AnnchromCS.compare_data a b
    | _, _ -> failwith_todo "compare_data"

let rec compare_union a b = 
    match a, b with
    | (U_SeqCS ha), (U_SeqCS hb) -> SeqCS.Union.compare_union ha hb 
    | U_Others, U_Others -> 0
    | _, _ -> failwith "DynamicCS.compare_union"


(** [to_formatter ref_codes attr t parent_t] returns
* dynamic character set [t] into Tag.output format *) 
let to_formatter report_type (node_name:string option) ref_codes attr t parent_t bl d : Xml.xml Sexpr.t list =
    match t,parent_t with 
        | MlCS t, None ->
            MlDynamicCS.to_formatter attr t None bl d
        | MlCS t, Some (MlCS parent_t) ->
            MlDynamicCS.to_formatter attr t (Some parent_t) bl d
        | SeqCS t, None ->
            SeqCS.to_formatter report_type attr t None d 
        | SeqCS t, Some (SeqCS parent_t) ->
            SeqCS.to_formatter report_type attr t (Some parent_t) d
        | ChromCS t, None ->
            ChromCS.to_formatter node_name ref_codes attr t None d
        | ChromCS t, Some (ChromCS parent_t) ->
            ChromCS.to_formatter node_name ref_codes attr t (Some parent_t) d
        | AnnchromCS t, None ->
            AnnchromCS.to_formatter ref_codes attr t None d
        | AnnchromCS t, Some (AnnchromCS parent_t) ->
            AnnchromCS.to_formatter ref_codes attr t (Some parent_t) d
        | BreakinvCS t, None ->
            BreakinvCS.to_formatter ref_codes attr t None d
        | BreakinvCS t, Some (BreakinvCS parent_t) ->
            BreakinvCS.to_formatter ref_codes attr t (Some parent_t) d
        | GenomeCS t, None ->
            GenomeCS.to_formatter ref_codes attr t None d
        | GenomeCS t, Some (GenomeCS parent_t) -> 
            GenomeCS.to_formatter ref_codes attr t (Some parent_t) d
        | (GenomeCS _ | AnnchromCS _ | ChromCS _ | BreakinvCS _ | MlCS _ | SeqCS _), _ ->
            failwith "DynamicCS.to_formatter; Inconsistent types"

(** [tabu_distance a_final b_final] returns the 
* tabu distance between dynamic character set [a_final] and [b_final] *)
let tabu_distance a_final b_final = 
    match a_final, b_final with 
    | _, MlCS b_final -> MlDynamicCS.tabu_distance b_final
    | _, SeqCS b_final -> (SeqCS.tabu_distance b_final)
    | ChromCS a_final, ChromCS b_final -> ChromCS.max_distance a_final b_final
    | GenomeCS a_final, GenomeCS b_final -> GenomeCS.max_distance a_final b_final
    | BreakinvCS a_final, BreakinvCS b_final -> BreakinvCS.max_distance a_final b_final
    | AnnchromCS a_final, AnnchromCS b_final -> AnnchromCS.max_distance a_final b_final
    | _, _ -> failwith_todo "tabu_distance"

(** [get_active_ref_code t] returns the set of active codes
* of dynamic character set [t] *)
let get_active_ref_code t = match t with
    | ChromCS t -> ChromCS.get_active_ref_code t
    | AnnchromCS t -> AnnchromCS.get_active_ref_code t
    | BreakinvCS t -> BreakinvCS.get_active_ref_code t
    | GenomeCS t -> GenomeCS.get_active_ref_code t
    | _ -> IntSet.empty, IntSet.empty

let mem c t = match c with
    | None    -> true
    | Some [] -> false
    | Some xs -> List.mem (code t) xs

let cardinal = function
    | MlCS x -> MlDynamicCS.cardinal x
    | SeqCS x -> SeqCS.cardinal x
    | ChromCS x -> ChromCS.cardinal x
    | GenomeCS x -> GenomeCS.cardinal x
    | BreakinvCS x -> BreakinvCS.cardinal x
    | AnnchromCS x -> AnnchromCS.cardinal x

let get_sequence_union code x = match x with
    | U_SeqCS x -> SeqCS.Union.get_sequence_union code x
    | U_Others -> failwith "DynamicCS.get_sequence_union"

let encoding enc x = match x with
    | MlCS x -> MlDynamicCS.encoding enc x
    | SeqCS x -> SeqCS.encoding enc x
    | _ -> failwith "Unsupported DynamicCS.encoding"

(** [classify_transformations] build a map/set of base frequencies and
    transformations that exist across this branch length defined by the two nodes
    passed to it. **)
let classify_transformations leafa nodea leafb nodeb chars acc = match nodea,nodeb with
    | SeqCS a, SeqCS b when mem chars nodea ->
        SeqCS.classify_transformations leafa a leafb b acc
    | SeqCS _, SeqCS _ ->
        acc
    (* We do not support other characters; yet? *)
    | MlCS _, MlCS _             -> failwith "DynamicML cannot classify transformations"
    | BreakinvCS _, BreakinvCS _ -> failwith "BreakInv cannot classify transformations"
    | AnnchromCS _, AnnchromCS _ -> failwith "Annchrom cannot classify transformations"
    | ChromCS    _, ChromCS    _ -> failwith "Chrom cannot classify transformations"
    | GenomeCS   _, GenomeCS   _ -> failwith "Genome cannot classify transformations"
    | (SeqCS _ | MlCS _ | BreakinvCS _ | AnnchromCS _ | ChromCS _ | GenomeCS _), _ -> assert false


let flatten t_lst = match List.hd t_lst with
    | BreakinvCS _ ->
        let bkCS_t_lst =
            List.map
                (fun x -> match x with
                    | BreakinvCS x_bkinvCS -> x_bkinvCS
                    | _ -> failwith "ERROR data type in flatten of dynamicCS.ml")
                t_lst
            in
            BreakinvCS.flatten bkCS_t_lst
    | _ -> failwith "we don't deal with this type of dynmaic data now"


(* return 1 if we are dealing with this kind of dynamic data for multi-chromosome.*)
let is_available in_data = match in_data with
    | BreakinvCS bk_t -> 1
    | SeqCS seq_t -> SeqCS.is_available seq_t
    | _ -> 0


let update_t oldt file_median_seq file_median_chrom_seqdeli =
    let newt = match oldt with
        | BreakinvCS bk_t ->
            let new_bk_t_lst =
                BreakinvCS.update_t bk_t file_median_seq file_median_chrom_seqdeli
            in
            List.map (fun x -> BreakinvCS x ) new_bk_t_lst
        |_ -> failwith "dynaicCS.update_t,we don't update this datatype for multi-chromosome functions under MGR"
    in
    newt
    
let single_to_multi single_t =
    let tlist = match single_t with
        |BreakinvCS bk_t ->
            List.map (fun bkinvCS_t -> BreakinvCS bkinvCS_t)
                     (BreakinvCS.single_to_multi bk_t)
        |_ -> failwith "we only deal with breakinv now"
    in
    tlist

(** [readjust ch1 ch2 par mine] attempts to (heuristically) readjust the character 
* set [mine] to somewhere in between [ch1], [ch2], and [par] (the children and
* parent of [mine] respectively). The function returns a triple [(a, b, c)],
* where [a] is the previous cost of [mine], [b] is the new cost of [c] as [ch1]
* and [ch2] parent, and [c] is the new readjusted [mine]. *)
(** [readjust_3d to_adjust modified ch1 ch2 ch1 ch2 parent mine]
* readjusts the current median [mine] of three medians [ch1],
* [ch2], and [parent] using three dimentional alignments*)
let readjust mode to_adjust modified ch1 ch2 parent mine =
    let no_iterative_other_than_for_seqs = false in
    match ch1, ch2, parent, mine with
    | SeqCS ch1, SeqCS ch2, SeqCS parent, SeqCS mine ->
        (*modified is a list of character id that get changed, new_cost is the cost(ch1,ch2), nc is the new SeqCS.t updated with new characters and new total_cost*)
        let modified, new_cost, new_sumcost, nc  = 
            if ch1.SeqCS.alph = Alphabet.nucleotides then  
                SeqCS.readjust ch1.SeqCS.alph mode to_adjust modified ch1 ch2 parent mine 
            else 
                SeqCS.readjust_custom_alphabet ch1.SeqCS.alph mode to_adjust modified ch1 ch2 parent mine
        in
	    modified,  new_cost, new_sumcost, (SeqCS nc)    
    | _, _, _, mine when no_iterative_other_than_for_seqs ->  
            let prev_cost = total_cost mine in
            modified, prev_cost,get_subtree_cost mine,  mine
    | ChromCS ch1, ChromCS ch2, ChromCS parent, ChromCS mine when ch1.ChromCS.alph =
        Alphabet.nucleotides -> 
            let modified, new_cost,new_sumcost, nc = 
                ChromCS.readjust to_adjust modified ch1 ch2 parent mine in
            modified, new_cost, new_sumcost,(ChromCS nc)
    | AnnchromCS ch1, AnnchromCS ch2, AnnchromCS parent, AnnchromCS mine 
        when ch1.AnnchromCS.alph = Alphabet.nucleotides -> 
          let modified, new_cost, new_sumcost, nc = 
              AnnchromCS.readjust to_adjust modified ch1 ch2 parent mine in
          modified, new_cost, new_sumcost, (AnnchromCS nc)
    | BreakinvCS ch1, BreakinvCS ch2, BreakinvCS parent, BreakinvCS mine ->
          let modified, new_cost,new_sumcost, nc = 
              BreakinvCS.readjust to_adjust modified ch1 ch2 parent mine in
          modified,  new_cost, new_sumcost, (BreakinvCS nc)
    | GenomeCS ch1, GenomeCS ch2, GenomeCS parent, GenomeCS mine ->
          let modified, new_cost, new_sumcost, nc = 
              GenomeCS.readjust to_adjust modified ch1 ch2 parent mine in
          modified, new_cost, new_sumcost, (GenomeCS nc)
    | _, _, _, mine ->  
            let prev_cost = total_cost mine in
            modified,  prev_cost, get_subtree_cost mine, mine

(* readjust the likelihood characters; has the additional branch length
 * arguments included *)
let readjust_lk mode to_adjust modified ch1 ch2 mine t1 t2 =
     match ch1, ch2, mine with
    | MlCS ch1, MlCS ch2, MlCS mine ->
        let m,pc,nc,ts,res = MlDynamicCS.readjust ch1 ch2 mine t1 t2 in
        if not m then (modified,pc,nc,(t1,t2),(MlCS mine))
        else begin
            let x =
                Array.fold_right
                    (fun c s -> All_sets.Integers.add c s) 
                    (MlDynamicCS.get_codes mine)
                    (All_sets.Integers.empty)
            in
            (x, pc, nc, ts, (MlCS res))
        end
     | _,_,_ -> assert false

(* readjust the likelihood characters; has the additional branch length
 * arguments included *)
let readjust_lk3 mode to_adjust modified ch1 ch2 mine par t1 t2 t3 =
    match ch1, ch2, mine, par with
    | MlCS ch1, MlCS ch2, MlCS mine, MlCS par ->
        let m,pcost,cost,ts,res = MlDynamicCS.readjust3_opt mine ch1 ch2 par t1 t2 t3 in
        if not m then (modified, pcost, pcost, (t1,t2,t3), (MlCS mine))
        else begin
            let x =
                Array.fold_right
                    (fun c s -> All_sets.Integers.add c s) 
                    (MlDynamicCS.get_codes mine)
                    (All_sets.Integers.empty)
            in
            (x, pcost, cost, ts, (MlCS res))
        end
     | _,_,_,_ -> assert false


let final_states mine ch1 ch2 par t1 t2 t3 = 
    match ch1,ch2,mine,par with
    | MlCS ch1, MlCS ch2, MlCS mine, MlCS par ->
        let res = MlDynamicCS.readjust3 mine ch1 ch2 par t1 t2 t3 in
        MlCS res
    | _,_,_,_ -> assert false


(** [to_single  ref_codes root alied_map p n] returns a node that contains per character a single state
 * which is closest to [p] among those available in [n].  
 * we don't modify cost of node to cost_to_parent_node anymore. all we do is
 * update sequence assigned to that node.
 * when root is passed, we intend to replace the sequence, leave alied
 * children sequence unchanged.
 * when root is None, we will use mine as root back to to_single function of
 * node.ml. 
 **)
let to_single ref_codes root parent mine time =
    match parent, mine with
    | SeqCS parent, SeqCS mine ->
        let parent, root  = match root with
            | None           -> parent, None
            | Some (SeqCS x) -> x, Some x
            | _              -> assert false
        in
        let prev_cost, new_cost, median = SeqCS.to_single parent mine root in
        prev_cost, new_cost, SeqCS median
    | MlCS parent, MlCS mine ->
        let min_bl = MlStaticCS.minimum_bl () in
        let parent = match root with
            | None          -> parent
            | Some (MlCS x) -> x
            | _             -> assert false
        and time = match time with
            | Some x -> max min_bl x
            (* we are dealing with a single node *)
            | None when MlDynamicCS.compare parent mine -> 0.0
            | None -> assert false
        in
        let p_cost, n_cost, med = MlDynamicCS.to_single parent mine time in
        p_cost, n_cost, MlCS med
    | ChromCS parent, ChromCS mine ->
        let root = match root with
            | Some (ChromCS root) -> Some root
            | None   -> None
            | Some _ -> assert false
        in
        let prev_cost, new_cost, median =
            ChromCS.to_single ref_codes root parent mine
        in
        prev_cost, new_cost, ChromCS median
    | BreakinvCS parent, BreakinvCS mine ->
        let root = match root with
            | Some (BreakinvCS root) -> Some root
            | None   -> None
            | Some _ -> assert false
        in
        let prev_cost, new_cost, median =
            BreakinvCS.to_single ref_codes root parent mine
        in
        prev_cost, new_cost, BreakinvCS median
    | AnnchromCS parent, AnnchromCS mine ->
        let root = match root with
            | Some (AnnchromCS root) -> Some root
            | None   -> None
            | Some _ -> assert false
        in
        let prev_cost, new_cost, median =
            AnnchromCS.to_single ref_codes root parent mine
        in
        prev_cost, new_cost, AnnchromCS median
    | GenomeCS parent, GenomeCS mine ->
        let root = match root with
            | Some (GenomeCS root) -> Some root
            | None   -> None
            | Some _ -> assert false
        in
        let prev_cost, new_cost, median =
            GenomeCS.to_single ref_codes root parent mine
        in
        prev_cost, new_cost, GenomeCS median
    | (SeqCS _ | GenomeCS _ | AnnchromCS _ | MlCS _ | BreakinvCS _ | ChromCS _), _ ->
        assert false

(** [to_single_root ref_codes mine] creates
* the single states for dynamic character set at root [mine] *)
let to_single_root ref_codes mine times = 
    to_single ref_codes (Some mine) mine mine times

