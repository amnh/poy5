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

let () = SadmanOutput.register "AllDirChar" "$Revision: 1616 $"

(** A dynamic character set implementation. 
* The dynamic character set allows rearrangements *)

exception Illegal_Arguments
let () = SadmanOutput.register "DynamicCS" "$Revision: 1006 $"


module IntMap = All_sets.IntegerMap
module IntSet = All_sets.Integers

exception No_Union

(** A dynamic character type. 'a can be SeqCS, ChromCS, BreakCS, GenomeCS ...*)
type t = 
    | SeqCS  of SeqCS.t
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
let alpha (a : t) = 
    match a with 
    | SeqCS a -> a.SeqCS.alph
    | ChromCS a -> a.ChromCS.alph
    | GenomeCS a -> a.GenomeCS.alph
    | BreakinvCS a -> a.BreakinvCS.alph
    | AnnchromCS a -> a.AnnchromCS.alph

(** [total_cost a] returns the total cost to create
*  dynamic character set [a]*)
let total_cost (a : t) = 
    match a with 
    | SeqCS a -> a.SeqCS.total_cost
    | ChromCS a -> a.ChromCS.total_cost
    | GenomeCS a -> a.GenomeCS.total_cost
    | BreakinvCS a -> a.BreakinvCS.total_cost
    | AnnchromCS a -> a.AnnchromCS.total_cost

(** [total_recost a] returns the total recost to create
*  dynamic character set [a]*)
let total_recost (a : t) = 
    match a with 
    | SeqCS a -> 0.
    | ChromCS a -> a.ChromCS.total_recost
    | GenomeCS a -> a.GenomeCS.total_recost
    | BreakinvCS a -> a.BreakinvCS.total_recost
    | AnnchromCS a -> a.AnnchromCS.total_recost

(** [subtree_recost a] returns the total recost of
* subtree whose root contains dynamic character set [a] *)
let subtree_recost (a : t) = 
    match a with 
    | SeqCS a -> 0.0
    | ChromCS a -> a.ChromCS.subtree_recost
    | GenomeCS a -> a.GenomeCS.subtree_recost
    | BreakinvCS a -> a.BreakinvCS.subtree_recost
    | AnnchromCS a -> a.AnnchromCS.subtree_recost

(** [c2 a] returns the two dimentional cost matrix
* of dynamic character set [a] *)
let c2 (a : t) = 
    match a with 
    | SeqCS a -> a.SeqCS.heuristic.SeqCS.c2
    | ChromCS a -> a.ChromCS.c2
    | GenomeCS a -> a.GenomeCS.c2
    | BreakinvCS a -> a.BreakinvCS.c2
    | AnnchromCS a -> a.AnnchromCS.c2

(** [chrom_pam a] returns the user-defined chromosome parameters
* of dynamic character set [a] *)
let chrom_pam (a : t) = 
    match a with 
    | ChromCS a -> a.ChromCS.chrom_pam
    | GenomeCS a -> a.GenomeCS.chrom_pam
    | AnnchromCS a -> a.AnnchromCS.annchrom_pam
    | BreakinvCS a -> a.BreakinvCS.breakinv_pam
    | _ -> Data.dyna_pam_default 

(** [state a] returns the character type of
* dynamic character set [a] *)
let state (a : t) : Data.dyna_state_t= 
    match a with 
    | SeqCS a -> `Seq
    | ChromCS a -> `Chromosome
    | GenomeCS a -> `Genome
    | BreakinvCS a -> `Breakinv
    | AnnchromCS a -> `Annotated

(** [code a] returns the code of dynamic character set [a] *)
let code (a : t) = 
    match a with 
    | SeqCS a -> a.SeqCS.code
    | ChromCS a -> a.ChromCS.code
    | GenomeCS a -> a.GenomeCS.code
    | BreakinvCS a -> a.BreakinvCS.code
    | AnnchromCS a -> a.AnnchromCS.code

let print dyn = 
    match dyn with
    | ChromCS ch ->
          ChromCS.print ch
    | _ -> print_endline "Do not print non-chromosome characters"

(** [copy_chrom_map s_ch d_ch] copies the choromosome map
* of dynamic character set [s_ch] to [d_ch] *)
let copy_chrom_map s_ch d_ch =
    match s_ch, d_ch with
    | ChromCS s_ch, ChromCS d_ch ->
          ChromCS (ChromCS.copy_chrom_map s_ch d_ch)
    | AnnchromCS s_ch, AnnchromCS d_ch ->
          AnnchromCS (AnnchromCS.copy_chrom_map s_ch d_ch)
    | GenomeCS s_ch, GenomeCS d_ch ->
          GenomeCS (GenomeCS.copy_chrom_map s_ch d_ch) 
    | _, _ -> d_ch

    
(** [leaf_sequences a] turns dynamic character set [a] 
* into a set of chromosome arrays *)
let leaf_sequences (a : t) = 
    match a with 
    | SeqCS a -> 
            let map = ref IntMap.empty in
            for i = (SeqCS.cardinal a) - 1 downto 0 do
                map := IntMap.add a.SeqCS.codes.(i)
                    (match a.SeqCS.characters.(i) with
                    | SeqCS.Heuristic_Selection x -> [|x.SeqCS.DOS.sequence|]
                    | SeqCS.Relaxed_Lifted (t, x) -> 
                            let p = SeqCS.RL.find_smallest x in
                            [|t.SeqCS.RL.sequence_table.(p)|]) !map
            done;
            !map
    | BreakinvCS a ->          
          IntMap.map 
              (fun med  -> 
                   [|(List.hd med.Breakinv.med_ls).BreakinvAli.seq|]
              ) a.BreakinvCS.meds


    | ChromCS a ->          
          IntMap.map 
              (fun med  -> 
                   [|(List.hd med.Chrom.med_ls).ChromAli.seq|]
              ) a.ChromCS.meds

    | AnnchromCS a ->
          IntMap.map 
              (fun med  -> 
                   let seqt_arr = 
                       (List.hd med.Annchrom.med_ls).AnnchromAli.seq_arr 
                   in
                   Array.map (fun seqt -> seqt.AnnchromAli.seq) seqt_arr
              ) a.AnnchromCS.meds


    | GenomeCS a ->
          IntMap.map 
              (fun med  -> 
                   let seq_arr = Array.map 
                       (fun chromt -> 
                            chromt.GenomeAli.seq
                       ) (List.hd med.Genome.med_ls).GenomeAli.chrom_arr 
                   in
                   seq_arr
              ) a.GenomeCS.meds

(** [unions a] returns the union of dynamic character set [a] *)
let unions (a : u) = 
    match a with 
    | U_SeqCS (Some a) -> 
            let map = ref IntMap.empty in
            for i = (Array.length a.SeqCS.Union.u_codes) - 1 
            downto 0 do
                map := IntMap.add (a.SeqCS.Union.u_codes.(i))
                    (a.SeqCS.Union.unions.(i)) !map
            done;
            !map
    | _ -> failwith "DynamicCS.unions"



let to_union a = 
    match a with 
    | SeqCS a -> U_SeqCS (SeqCS.to_union a)
    | _ -> U_Others

let union a b c =
    match a, b, c with
    | SeqCS a, U_SeqCS b, U_SeqCS c -> U_SeqCS (SeqCS.Union.union a b c)
    | SeqCS _, _, _ -> failwith "DynamicCS.union"
    | _, U_Others, U_Others -> U_Others
    | _, _, _ -> failwith_todo "union"


let cardinal_union a =
    match a with
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
    | `Seq | `Breakinv | `Chromosome as meth ->
            let seq_arr = 
                Array.map 
                (fun (genome_data, genome_code) ->
                    let first_seq = genome_data.Data.seq_arr.(0).Data.seq in  
                    (first_seq, genome_code)) genome_arr
            in 
            begin match meth with
            | `Seq -> 
                    let t = SeqCS.of_array spec seq_arr code taxon in
                    SeqCS t
            | `Breakinv -> 
                    let t = BreakinvCS.of_array spec seq_arr code in
                    BreakinvCS t
            | `Chromosome  ->
                    let t = 
                        ChromCS.of_array spec seq_arr code taxon num_taxa 
                    in 
                    ChromCS t
            end
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
let median code a b =
    match a, b with 
    | SeqCS a, SeqCS b -> SeqCS (SeqCS.median code a b)
    | ChromCS a, ChromCS b -> ChromCS (ChromCS.median2 a b)
    | GenomeCS a, GenomeCS b -> GenomeCS (GenomeCS.median2 a b)
    | BreakinvCS a, BreakinvCS b -> BreakinvCS (BreakinvCS.median2 a b)
    | AnnchromCS a, AnnchromCS b -> AnnchromCS (AnnchromCS.median2 a b)
    | _, _ -> failwith_todo "median"

(** [median_3 p n c1 c2] creates the median among
* three dynamic character sets [p], [c1] and [c2] *)
let median_3 p n c1 c2 =
    match p, n, c1, c2 with 
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


(* Like [distance] but calculates it only 
* if the type of the characters match one of those listed. *)
let distance_of_type t a b =
    let has_t x = List.exists (fun z -> z = x) t in
    let has_seq = has_t `Seq 
    and has_chrom = has_t `Chrom 
    and has_gen = has_t `Genome
    and has_break = has_t `Breakinv
    and has_ann = has_t `Annchrom in
    match a, b with
    | SeqCS a, SeqCS b when has_seq -> (SeqCS.distance a b)
    | ChromCS a, ChromCS b when has_chrom -> ChromCS.distance a b  
    | GenomeCS a, GenomeCS b when has_gen -> GenomeCS.distance a b  
    | BreakinvCS a, BreakinvCS b when has_break -> BreakinvCS.distance a b  
    | AnnchromCS a, AnnchromCS b when has_ann -> AnnchromCS.distance a b  
    | _, _ -> 0.0

(** [distance_of_type a b] returns the distance between
* two dynamic character sets [a] and [b] *)
let distance a b = 
    match a, b with   
    | SeqCS a, SeqCS b -> (SeqCS.distance a b)
    | ChromCS a, ChromCS b -> ChromCS.distance a b  
    | GenomeCS a, GenomeCS b -> GenomeCS.distance a b  
    | BreakinvCS a, BreakinvCS b -> BreakinvCS.distance a b  
    | AnnchromCS a, AnnchromCS b -> AnnchromCS.distance a b  
    | _, _ -> failwith_todo "distance"  



(** [distance_union a b] returns the union distance between
* two dynamic character sets [a] and [b] *)
let distance_union a b =
    match a, b with
    | U_SeqCS a, U_SeqCS b -> SeqCS.Union.distance_union a b
    | U_Others, U_Others -> 0.0
    | _, _ -> failwith "DynamicCS.distance_union"


(** [to_string a] returns dynamic character set [a] 
* into the string format *)
let to_string a =
    match a with 
    | SeqCS a -> SeqCS.to_string a
    | BreakinvCS a -> BreakinvCS.to_string a
    | ChromCS a -> ChromCS.to_string a
    | GenomeCS a -> GenomeCS.to_string a
    | AnnchromCS a -> AnnchromCS.to_string a


(* [dist_2 delta n a b] calculates the cost of joining 
* the node containing  chromosome character set [n]
* between two nodes containing [a] and [b]. 
* [a] must be the parent (ancestor) of [b] *)
let dist_2 delta n a b =
    match n, a, b with 
    | SeqCS n, SeqCS a, SeqCS b -> (SeqCS.dist_2 delta n a b)
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
let to_formatter ref_codes attr t (parent_t : t option) d : Tags.output list = 
    match t, parent_t with 
    | SeqCS t, _ -> begin
            match parent_t with
            | None -> SeqCS.to_formatter attr t None d 
            | Some (SeqCS parent_t) ->  
                    SeqCS.to_formatter attr t (Some parent_t) d
            | _ -> failwith "to_formatter in dynamicCS"
    end
    | ChromCS t, _  -> begin 
          match parent_t with 
          | None ->  ChromCS.to_formatter ref_codes attr t None d
          | Some (ChromCS parent_t) -> ChromCS.to_formatter ref_codes attr t (Some parent_t) d
          | _ -> failwith "to_formatter in dynamicCS"
      end 
    | AnnchromCS t, _ -> begin
          match parent_t with 
          | None ->  AnnchromCS.to_formatter ref_codes attr t None d
          | Some (AnnchromCS parent_t) -> AnnchromCS.to_formatter ref_codes attr t (Some parent_t) d
          | _ -> failwith "to_formatter in dynamicCS"
      end 
    | BreakinvCS t, __ -> begin
          match parent_t with 
          | None ->  BreakinvCS.to_formatter ref_codes attr t None d
          | Some (BreakinvCS parent_t) ->BreakinvCS.to_formatter ref_codes attr t (Some parent_t) d
          | _ -> failwith "to_formatter in dynamicCS"
      end 
    | GenomeCS t,  _ -> begin
          match parent_t with
          | None -> GenomeCS.to_formatter ref_codes attr t None d
          | Some (GenomeCS parent_t) -> GenomeCS.to_formatter ref_codes attr t (Some parent_t) d
          | _ -> failwith "to_formatter in dynamicCS"
      end 



(** [tabu_distance a_final b_final] returns the 
* tabu distance between dynamic character set [a_final] and [b_final] *)
let tabu_distance a_final b_final = 
    match a_final, b_final with 
    | _, SeqCS b_final -> (SeqCS.tabu_distance b_final)
    | ChromCS a_final, ChromCS b_final -> ChromCS.max_distance a_final b_final
    | GenomeCS a_final, GenomeCS b_final -> GenomeCS.max_distance a_final b_final
    | BreakinvCS a_final, BreakinvCS b_final -> BreakinvCS.max_distance a_final b_final
    | AnnchromCS a_final, AnnchromCS b_final -> AnnchromCS.max_distance a_final b_final
    | _, _ -> failwith_todo "tabu_distance"

(** [get_active_ref_code t] returns the set of active codes
* of dynamic character set [t] *)
let get_active_ref_code t = 
    match t with
    | ChromCS t -> ChromCS.get_active_ref_code t
    | AnnchromCS t -> AnnchromCS.get_active_ref_code t
    | BreakinvCS t -> BreakinvCS.get_active_ref_code t
    | GenomeCS t -> GenomeCS.get_active_ref_code t
    | _ -> IntSet.empty, IntSet.empty


let cardinal x =
    match x with
    | SeqCS x -> SeqCS.cardinal x
    | ChromCS x -> ChromCS.cardinal x
    | GenomeCS x -> GenomeCS.cardinal x
    | BreakinvCS x -> BreakinvCS.cardinal x
    | AnnchromCS x -> AnnchromCS.cardinal x

let get_sequence_union code x = 
    match x with
    | U_SeqCS x -> SeqCS.Union.get_sequence_union code x
    | U_Others -> failwith "DynamicCS.get_sequence_union"

let encoding enc x =
    match x with
    | SeqCS x -> SeqCS.encoding enc x
    | _ -> failwith "Unsupported DynamicCS.encoding"

(* We are turning off iterative for high order characters until the algorithms
* are properly fixed. *)
let no_iterative_other_than_for_seqs = true

(** [readjust ch1 ch2 par mine] attempts to (heuristically) readjust the character 
* set [mine] to somewhere in between [ch1], [ch2], and [par] (the children and
* parent of [mine] respectively). The function returns a triple [(a, b, c)],
* where [a] is the previous cost of [mine], [b] is the new cost of [c] as [ch1]
* and [ch2] parent, and [c] is the new readjusted [mine]. *)
(** [readjust_3d to_adjust modified ch1 ch2 ch1 ch2 parent mine]
* readjusts the current median [mine] of three medians [ch1],
* [ch2], and [parent] using three dimentional alignments*)
let readjust mode to_adjust modified ch1 ch2 parent mine =
    match ch1, ch2, parent, mine with
    | SeqCS ch1, SeqCS ch2, SeqCS parent, SeqCS mine when ch1.SeqCS.alph =
        Alphabet.nucleotides -> 
            let modified, new_cost, nc = 
                SeqCS.readjust mode to_adjust modified ch1 ch2 parent mine in
            let prev_cost = SeqCS.distance ch1 mine +. SeqCS.distance ch2 mine in
            modified, prev_cost, new_cost, (SeqCS nc)

    | _, _, _, mine when no_iterative_other_than_for_seqs ->  
            let prev_cost = total_cost mine in
            modified, prev_cost, prev_cost, mine

    | ChromCS ch1, ChromCS ch2, ChromCS parent, ChromCS mine when ch1.ChromCS.alph =
        Alphabet.nucleotides -> 
            let modified, new_cost, nc = 
                ChromCS.readjust to_adjust modified ch1 ch2 parent mine in
            let prev_cost = ChromCS.distance ch1 mine +. ChromCS.distance ch2 mine in
            modified, prev_cost, new_cost, (ChromCS nc)

    | AnnchromCS ch1, AnnchromCS ch2, AnnchromCS parent, AnnchromCS mine when ch1.AnnchromCS.alph =
          Alphabet.nucleotides -> 
          let modified, new_cost, nc = 
              AnnchromCS.readjust to_adjust modified ch1 ch2 parent mine in
          let prev_cost = AnnchromCS.distance ch1 mine +. AnnchromCS.distance ch2 mine in
          modified, prev_cost, new_cost, (AnnchromCS nc)


    | BreakinvCS ch1, BreakinvCS ch2, BreakinvCS parent, BreakinvCS mine ->
          let modified, new_cost, nc = 
              BreakinvCS.readjust to_adjust modified ch1 ch2 parent mine in
          let prev_cost = BreakinvCS.distance ch1 mine +. BreakinvCS.distance ch2 mine in
          modified, prev_cost, new_cost, (BreakinvCS nc)


    | GenomeCS ch1, GenomeCS ch2, GenomeCS parent, GenomeCS mine ->
          let modified, new_cost, nc = 
              GenomeCS.readjust to_adjust modified ch1 ch2 parent mine in
          let prev_cost = GenomeCS.distance ch1 mine +. GenomeCS.distance ch2 mine in
          modified, prev_cost, new_cost, (GenomeCS nc)
              

    | _, _, _, mine ->  
            let prev_cost = total_cost mine in
            modified, prev_cost, prev_cost, mine

(** [to_single ?is_root pre_ref_code alied_map p n] returns a node that contains per character a single state
 * which is closest to [p] among those available in [n]. Useful for tree length
 * verification. is_root optional paramter indicates that if n is root. The
 * default is false. pre_ref_code contains active codes for chromosome characters. 
 * Inactive codes are eliminated from diagnosis. 
 * If p is the handle, alied_map is the root containing the aligned map between p
 * and n for chromosome stuff, else alied_map is assigned by p *)
let to_single ref_codes root parent mine = 
    match parent, mine with
    | SeqCS parent, SeqCS mine -> 
            let parent = 
                match root with
                | None -> parent
                | Some (SeqCS x) -> x
                | _ -> assert false
            in
            let prev_cost, new_cost, median = SeqCS.to_single parent mine in
            prev_cost, new_cost, SeqCS median
    | ChromCS parent, ChromCS mine -> 
          let root = match root with 
          | Some (ChromCS root) -> Some root
          | _ -> None
          in 
          let prev_cost, new_cost, median = 
              ChromCS.to_single ref_codes root parent mine 
          in
          prev_cost, new_cost, ChromCS median          


    | BreakinvCS parent, BreakinvCS mine ->
          let root = match root with 
          | Some (BreakinvCS root) -> Some root
          | _ -> None
          in 
          let prev_cost, new_cost, median = 
              BreakinvCS.to_single ref_codes root parent mine 
          in
          prev_cost, new_cost, BreakinvCS median          


    | AnnchromCS parent, AnnchromCS mine ->
          let root = match root with 
          | Some (AnnchromCS root) -> Some root
          | _ -> None
          in 
          let prev_cost, new_cost, median = 
              AnnchromCS.to_single ref_codes root parent mine 
          in
          prev_cost, new_cost, AnnchromCS median          

    | GenomeCS parent, GenomeCS mine ->
          let root = match root with 
          | Some (GenomeCS root) -> Some root
          | _ -> None
          in 
          let prev_cost, new_cost, median = 
              GenomeCS.to_single ref_codes root parent mine 
          in
          prev_cost, new_cost, GenomeCS median          
    | _, mine ->

            let cst = total_cost mine in
            cst, cst, mine

(** [to_single_root ref_codes mine] creates
* the single states for dynamic character set at root [mine] *)
let to_single_root ref_codes mine = to_single ref_codes (Some mine) mine mine


module Kolmogorov = struct
    let correct_cost t m = 
        match t with
        | SeqCS t -> SeqCS (SeqCS.Kolmogorov.correct_cost t m)
        | _ -> t
end
