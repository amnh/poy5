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

(** A chromosome character set implementation. 
* The chromosome character set allows rearrangements *)

exception Illegal_Arguments
let () = SadmanOutput.register "ChromCS" "$Revision: 2754 $"

let fprintf = Printf.fprintf

module IntMap = All_sets.IntegerMap
module IntSet = All_sets.Integers

type meds_t = Chrom.meds_t

(** A chromosme character set *)
type t = { 
    meds : meds_t IntMap.t; (** a set of chromosome characters *)
    costs : float IntMap.t; (** cost2 from two children*)
    recosts : float IntMap.t;
    total_cost : float;             (** The total cost of the character set:
        cost from aligning its children *)
    subtree_cost : float;           (** The total_cost of subtree root in this node*)
    total_recost : float;           (** The total recost of the character set *)
    subtree_recost : float;         (** The total subtree recost of the
                                        character set *)
    c2_full : Cost_matrix.Two_D.m;       (** The two dimensional cost matrix to be used in the character set *)
    c2_original : Cost_matrix.Two_D.m;       (** The two dimensional cost matrix to be used in the character set *)
    c3 : Cost_matrix.Three_D.m;     (** The three dimensional cost matrix to be used in the character set *)
    chrom_pam : Data.dyna_pam_t;
    alph : Alphabet.a;              (** The alphabet of the sequence set *)
    code : int;                     (** The set code *)
}

let cardinal x = IntMap.fold (fun _ _ x -> x + 1) x.meds 0


(** [of_array spec arr chcode tcode num_taxa] 
* creates a  chromosome character set from an array of sequences [arr] *)
let of_array spec arr chcode tcode num_taxa = 
    let adder (meds, costs, recosts) (seq, key) = 
        let med = Chrom.init_med  seq spec.Data.tcm2d_full
        spec.Data.tcm2d_original spec.Data.pam tcode num_taxa in 
        (IntMap.add key med meds), 
        (IntMap.add key 0.0 costs), 
        (IntMap.add key 0.0 recosts)
    in
    let meds, costs, recosts = 
        let empty = IntMap.empty in
        Array.fold_left adder (empty, empty, empty) arr 
    in
    {
        meds = meds;
        costs = costs;
        recosts = recosts;
        total_cost = 0.0;
        subtree_cost = 0.0;
        total_recost = 0.0;
        subtree_recost = 0.0;
        c2_full = spec.Data.tcm2d_full;
        c2_original = spec.Data.tcm2d_original;
        c3 = spec.Data.tcm3d;
        alph = spec.Data.alph;
        code = chcode;
        chrom_pam = spec.Data.pam;
    }

(** [of_list spec arr chcode tcode num_taxa] creates 
* a chromosome character set from a list of sequences [lst] *)
let of_list spec lst =
    let arr = Array.of_list lst in
    of_array spec arr 

(** [to_list t] returns a list of chromosome characters 
* from a set of chromosome characters *)
let to_list t =
    IntMap.fold (fun code med acc -> (med, code) :: acc) t.meds []

(** [same_codes a b] returns true if code set of [a] 
* is the same to code set of [b], otherwise false *)
let same_codes a b =
    let checker x _ res = res && (IntMap.mem x b) in
    IntMap.fold checker a true

let print ch =
    IntMap.iter (fun _ med -> Chrom.print med) ch.meds


(** [median2 a b] returns the chromosome median set 
* between  chromosome character sets [a] and [b] *)
let median2 (a : t) (b : t) =
    (* We will use imperative style for this function *)
    let empty = IntMap.empty in

    let median code (meda : meds_t) (medians, costs, recosts, total_cost, total_recost) = 
        let medb : meds_t = IntMap.find code b.meds in
        let medab = Chrom.find_meds2 meda medb in
        
        let new_median = IntMap.add code medab medians 
        and new_costs = 
            IntMap.add code (float_of_int medab.Chrom.total_cost) costs  
        and new_recosts = IntMap.add code (float_of_int medab.Chrom.total_recost) recosts  
        and new_total_cost = total_cost + medab.Chrom.total_cost 
        and new_total_recost = total_recost + medab.Chrom.total_recost 
        in 
        new_median, new_costs, new_recosts, new_total_cost, new_total_recost
    in
    let medab_map, new_costs, new_recosts, total_cost, total_recost =  
        IntMap.fold median a.meds (empty, empty, empty, 0, 0) 
    in
    let subtree_recost = a.subtree_recost +. b.subtree_recost +. (float_of_int total_recost) in 
    let subtree_cost = a.subtree_cost +. b.subtree_cost +. (float_of_int total_cost) in
    { a with meds = medab_map; costs = new_costs; recosts = new_recosts;
          total_cost = float_of_int total_cost; 
          subtree_cost = subtree_cost;
          total_recost = float_of_int total_recost;
          subtree_recost = subtree_recost;
    }
    

(** [median3 p n c1 c2] returns the chromosome median set of
*  chromosome character sets [p], [c1] and [c2] *)
let median3 p n c1 c2 =
(*    print_endline "median3 in ChromCs module"; *)
    let median code  medp res_medians = 
        let med1= IntMap.find code c1.meds in 
        let med2 = IntMap.find code c2.meds in
        
        let medp12 = Chrom.find_meds3 medp med1 med2 in
          IntMap.add code medp12 res_medians 
(*        let med12 = Chrom.find_meds2 med1 med2 p.c2 in 
        IntMap.add code med12 res_medians *)
    in
    let acc = IntMap.empty in
    let medp12_map = IntMap.fold median p.meds acc in
    { n with meds = medp12_map; }


(** [readjust to_adjust modified ch1 ch2 parent mine] returns
* the readjusted median set [mine] of three  chromosome character
* sets [ch1], [ch2] and [parent] *)
let readjust to_adjust modified ch1 ch2 parent mine = 
    let empty = IntMap.empty and
            c2 = parent.c2_full and
            c3 = parent.c3 
    in
    let adjusted code parent_chrom acc =
        let to_adjust =
            match to_adjust with
            | None -> All_sets.Integers.singleton code
            | Some x -> x
        in
        let (modified, res_medians, res_costs, total,totalsum) = acc in
        let my_chrom = IntMap.find code mine.meds
        and ch1_chrom = IntMap.find code ch1.meds
        and ch2_chrom = IntMap.find code ch2.meds in
        if (not (All_sets.Integers.mem code to_adjust)) then 
            let new_costs = IntMap.add code 0. res_costs 
            and new_single = IntMap.add code my_chrom res_medians in
            modified, new_single, new_costs, total, totalsum
        else begin
            let rescost,ressumcost, seqm, changed = 
                Chrom.readjust_3d ch1_chrom ch2_chrom my_chrom
                    c2 c3 parent_chrom
            in
            let new_single = IntMap.add code seqm res_medians
            and new_costs = IntMap.add code (float_of_int rescost) res_costs 
            and new_total = total + rescost in
            let new_totalsum = ressumcost + totalsum in
            let modified = 
                if changed then All_sets.Integers.add code modified
                else modified
            in
            modified, new_single, new_costs, new_total, new_totalsum       
        end 
    in 
    let modified, meds, costs, total_cost, total_sum_cost = 
        IntMap.fold adjusted parent.meds (modified, empty, empty, 0, 0)
    in
    let tc = float_of_int total_cost in
    let tsc = float_of_int total_sum_cost in
    modified,
    tc,
    tsc,
    { mine with meds = meds; costs = costs; total_cost = tc; subtree_cost = tsc; }

(**[get_extra_cost_for_root root] return the extra cost result from non-zero
* diagonal in cost matrix.*)
let get_extra_cost_for_root (a :t) =
    let get_ec code medst acc =
        acc + Chrom.get_extra_cost_for_root medst a.c2_original 
    in
    let is_identity = Cost_matrix.Two_D.is_identity a.c2_original in
    if is_identity then 0. (*0 diagonal in cost matrix*)
    else (*non-0 diagonal in cost matrix*)
        float_of_int (IntMap.fold get_ec a.meds 0) 

(** [distance a b] returns total distance between 
* two  chromosome character sets [a] and [b] *)
let distance (a : t) (b : t)  = 
    let single_distance code meda (acc_cost, acc_recost) =
        let medb = IntMap.find code b.meds in
        let cost, recost = Chrom.cmp_min_pair_cost meda medb in 
        acc_cost + cost, acc_recost + recost
    in
    let cost, _ = IntMap.fold single_distance a.meds (0, 0) in 
    float_of_int cost


(** [max_distance a b] returns total maximum distances
* between two  chromosome character sets [a] and [b]*)
let max_distance (a : t) (b : t)  = 
    let single_distance code meda (acc_cost, acc_recost) =
        let medb = IntMap.find code b.meds in
        let cost, recost = Chrom.cmp_max_pair_cost meda medb in 
        acc_cost + cost, acc_recost + recost
    in
    let cost, _ = IntMap.fold single_distance a.meds (0,0) in 
    float_of_int cost

(** [to_string a] converts chromosome character set [a]
* into string format *)
let to_string a =
    let builder code med acc =
        let code = string_of_int code in 
        let seq_ls = List.map (fun med -> Sequence.to_string med.ChromAli.seq a.alph) 
            med.Chrom.med_ls 
        in 

        let seq = String.concat ":" seq_ls in 
        acc ^ code ^ ": " ^ seq ^ "; "
    in
    IntMap.fold builder a.meds ""


(* [dist_2 n a b] calculates the cost of joining 
* the node containing  chromosome character set [n]
* between two nodes containing [a] and [b]. 
* [a] must be the parent (ancestor) of [b] *)
let dist_2 n a b =
    let cost_calculator code medb (acc_cost, acc_recost) =
        let medn = IntMap.find code n.meds
        and meda = IntMap.find code a.meds 
        and medb = IntMap.find code b.meds in
        
        let medab : meds_t = Chrom.find_meds2 meda medb in 
        
        let cost, recost = Chrom.cmp_min_pair_cost medn medab in 
        acc_cost + cost, acc_recost + recost
    in
    let cost, _ =     IntMap.fold cost_calculator b.meds (0,0) in 
    float_of_int cost


(** [f_codes s c] returns a chromosome 
* character subset of  chromosome 
* character set [s] whose codes are 
* also in  chromosome character set [c] *)
let f_codes s c = 
    let check x = All_sets.Integers.mem x c in
    let adder x y acc = 
        if check x then IntMap.add x y acc 
        else acc
    in
    let n_meds = IntMap.fold adder s.meds IntMap.empty
    and n_costs = IntMap.fold adder s.costs IntMap.empty in
    { s with meds = n_meds; costs = n_costs}

(** [f_codes_comp s c] returns  chromosome character subset
* of  chromosome character set [s] whose codes are NOT
* also in  chromosome character set [c] *)
let f_codes_comp s c = 
    let check x = not (All_sets.Integers.mem x c) in
    let adder x y acc = 
        if check x then IntMap.add x y acc 
        else acc
    in
    let n_meds = IntMap.fold adder s.meds IntMap.empty
    and n_costs = IntMap.fold adder s.costs IntMap.empty in
    { s with meds = n_meds; costs = n_costs }


(** [compare_data a b] compares  chromosome character
* set [a] and [b]. Returns 0 if they are the same, 
* otherwise (-1) or (1) *)
let compare_data a b =
    let comparator code medb acc =
        if acc = 0 then begin
            let meda = IntMap.find code a.meds in
            Chrom.compare meda medb
        end else acc
    in
    IntMap.fold comparator b.meds 0


(** [to_formatter ref_codes attr t parent_t d] returns
* the map between  chromosme character set [t] and 
* its parents [parent_t] in the Tag.output format *)
let to_formatter (node_name:string option) ref_codes attr t (parent_t : t option) d : Xml.xml Sexpr.t list = 
    let _, state = List.hd attr in   
    let output_chrom code med acc =
        let med = 
            try
                List.find (fun med ->                     
                               IntSet.mem med.ChromAli.ref_code  ref_codes
                          ) med.Chrom.med_ls
            with Not_found -> List.hd med.Chrom.med_ls
        in         
        let cost, recost,  map = 
            match parent_t with  
            | None -> 0, 0, None
            | Some parent -> begin 
                  let parent_med = IntMap.find code parent.meds in  
                  let parent_med =  
                      try 
                          List.find 
                              (fun med -> 
                                   IntSet.mem med.ChromAli.ref_code ref_codes 
                              ) parent_med.Chrom.med_ls
                      with Not_found -> List.hd parent_med.Chrom.med_ls
                  in 
                  let cost, recost, map = 
                      match state with
                      | `String "Preliminary" ->
                            ChromAli.create_map parent_med med.ChromAli.ref_code  
                      | `String "Final" ->
                            ChromAli.create_map med parent_med.ChromAli.ref_code   
                      | _ ->      
                            let filename = 
                                match node_name with 
                                | Some x -> 
                                        Some ("mauve_node#"^x^".alignment")
                                | None -> None
                            in
                            let cost, recost, med_ls = 
                                ChromAli.find_med2_ls 
                                med parent_med t.c2_original t.chrom_pam filename
                            in 
                            let med = List.hd med_ls in                             
                            let map = ChromAli.create_single_map med in 
                            cost, recost, map
                  in 
                  cost, recost, Some map
              end 
        in 
        let seq = Sequence.to_formater med.ChromAli.seq t.alph in
        let name = Data.code_character code d in  
        let cost_str = 
            match state with
            | `String "Single" -> `IntTuple (cost, cost)
            | _ -> `IntTuple (0, cost)
        in 
        let module T = Xml.Characters in
        (PXML 
            -[T.chromosome] 
                (* Attributes *)
                ([T.name] = [`String name])
                ([T.cost] = [cost_str])
                ([T.recost] = [`Int recost])
                ([T.definite] = [`Bool (cost > 0)])
                ([T.ref_code] = [`Int med.ChromAli.ref_code])
                ([attr])

                (* Contents, first the optional map, then everything else *)
                { match map with
                | None -> (PXML ---)
                | Some map -> `Single map}

                -[T.sequence] { `String seq } --
            --) :: acc
    in
    IntMap.fold output_chrom t.meds []


(** [to_single ref_codes root single_parent mine] returns
* the single states of  chromosome character set [mine] *) 
let to_single ref_codes (root : t option) single_parent mine = 
    let previous_total_cost = mine.total_cost in 
    let c2_full = mine.c2_full in 
    let c2_original = mine.c2_original in 

    let median code med (acc_meds, acc_costs, acc_recosts, acc_total_cost) =        

        let amed = 
            try
                List.find (fun med -> 
                               IntSet.mem med.ChromAli.ref_code ref_codes
                          ) med.Chrom.med_ls
            with Not_found -> begin
(*                failwith "Not_found: ChromCS.to_single"*)
                List.hd med.Chrom.med_ls
            end 
        in         


        let parent_med = IntMap.find code single_parent.meds in  
        let aparent_med = 
            try
                List.find  
                    (fun med -> 
                         IntSet.mem med.ChromAli.ref_code ref_codes 
                    ) parent_med.Chrom.med_ls
            with Not_found -> begin
                List.hd parent_med.Chrom.med_ls
            end 
        in            

        let cost,  recost, single_seq = 
            match root with
            | None -> 
                  let single_seq = 
                      ChromAli.to_single aparent_med amed.ChromAli.ref_code
                      c2_full  med.Chrom.chrom_pam
                  in 

                  let cost, recost = ChromAli.cmp_cost 
                      ({amed with ChromAli.seq = (Sequence.delete_gap
                      single_seq)} ) aparent_med c2_original
                      med.Chrom.chrom_pam `Chromosome
                  in 
                  cost, recost, single_seq

            | Some root ->              
                  let single_root = ChromAli.to_single_root amed
                      aparent_med.ChromAli.ref_code c2_full 
                  in
                  0, 0, single_root
        in 
        
        let single_med = ChromAli.change_to_single amed single_seq c2_original in                     
        let single_med = {med with Chrom.med_ls = [single_med]} in 
        let new_single = IntMap.add code single_med acc_meds in
        let new_costs = IntMap.add code (float_of_int cost) acc_costs in 
        let new_recosts = IntMap.add code (float_of_int recost) acc_recosts in 

        new_single, new_costs, new_recosts, (acc_total_cost + cost)
    in
    let meds, costs,  recosts, total_cost = 
        match root with
        | Some root ->
              IntMap.fold median root.meds (IntMap.empty, IntMap.empty, IntMap.empty, 0)
        | None ->
              IntMap.fold median mine.meds (IntMap.empty, IntMap.empty, IntMap.empty, 0)
    in 


    previous_total_cost, float_of_int total_cost, 
    {mine with meds = meds; 
         costs = costs;
         recosts = recosts;
         total_cost = float_of_int total_cost}
         

(** [get_active_ref_code t] returns active reference codes
* of  chromosome character set [t]. 
* One active reference code for one  chromosome character *)
let get_active_ref_code t = 
    IntMap.fold 
        (fun _ meds (acc_ref_code, acc_child_ref_code) ->
             let ref_code, child1_ref_code, child2_ref_code = 
                 Chrom.get_active_ref_code meds 
             in
             IntSet.add ref_code acc_ref_code,
             IntSet.add child2_ref_code (IntSet.add child1_ref_code acc_child_ref_code)  
        ) t.meds (IntSet.empty, IntSet.empty)






(** [copy_chrom_map s_ch d_ch] copies the chromosome
* map of  chromosome set [s_ch] to chromosome set [d_ch] *)
let copy_chrom_map s_ch d_ch =
    let copied_meds = IntMap.mapi 
        (fun code ad_ch ->
             let as_ch = IntMap.find code s_ch.meds in 
             let copied_ad = Chrom.copy_chrom_map as_ch ad_ch in 
             copied_ad
        ) d_ch.meds
    in 
    {d_ch with meds = copied_meds}

