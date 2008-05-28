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

let () = SadmanOutput.register "Genome" "$Revision: 975 $"

(** Genome module implements functions to create medians
    between two lists of genomes *)

let fprintf = Printf.fprintf


type med_t = GenomeAli.med_t

(** [meds_t] is a data structure for a list of  medians 
* between two genomes. Rearrangements are allowed *)
type meds_t = {
    med_ls : med_t list; (** genome list *)
    total_cost : int;   (** the total cost to create this genome list *)
    total_recost : int;  (** total recost to create this genome list *)
    c2 : Cost_matrix.Two_D.m; 

    approx_cost_arr : int array; (* An array of genome medians 
                                  * created from this median to other medians.
                                  * This is used for approximation *)
    approx_recost_arr : int array;
    code : int; (** the taxa code containing this genome list *)

    chrom_pam : Data.dyna_pam_t (** user-defined parameters used for creating genome medians *)   
}

let max_taxa_id = ref 0 

(** [init_med genome chrom_pam tcode num_taxa] 
* returns a genome list with only one element
* created from genome [genome] *) 
let init_med (genome : Sequence.s Data.dyna_data) chrom_pam tcode num_taxa
        =  
    let med = GenomeAli.init genome in 
     max_taxa_id := max !max_taxa_id num_taxa;
   
    { 
        med_ls = [med]; 
    
        total_cost = 0; 
        total_recost = 0; 
        c2 = Cost_matrix.Two_D.default; 
 
        approx_cost_arr = (Array.make !max_taxa_id max_int);
        approx_recost_arr = (Array.make !max_taxa_id max_int);
 
       code = tcode;

        chrom_pam = chrom_pam 
    } 


(** [keep chrom_pam med_ls] selects randomly a subset 
* of medians kept to process further *)
let rec keep chrom_pam med_ls = 
    match chrom_pam.Data.keep_median with 
    | None -> med_ls 
    | Some keep_median ->
          if  keep_median >= List.length med_ls then med_ls
          else Utl.get_k_random_elem med_ls keep_median


(** [update_cost_mat meds1 meds2] creates the median
* between genome [meds1] and [meds2]
* and update their approximate cost, recost
* arrays if they are not yet computed *)
let update_cost_mat meds1 meds2 =     
    let code1 = meds1.code and code2 = meds2.code in  
    if meds1.approx_cost_arr.(code2) = max_int then begin 

        let med1 = List.hd meds1.med_ls in  
        let med2 = List.hd meds2.med_ls in  
        let cost, (recost1, recost2) = GenomeAli.cmp_cost med1 med2 meds1.c2 meds1.chrom_pam in  
        meds1.approx_cost_arr.(code2) <- cost; 
        meds2.approx_cost_arr.(code1) <- cost; 

        meds1.approx_recost_arr.(code2) <- (recost1 + recost2); 
        meds2.approx_recost_arr.(code1) <- (recost1 + recost2); 
    end


(** [find_meds2 keep_all_meds meds1 meds2] finds median list
 * between two lists of medians [meds1=(x1,...,xk)] and [meds2=(y1,...,yt)]
 * where xi and yj are medians. For each pair (xi, yj) we have 
 * a list of medians z_ij with the same cost c_ij. 
 * Find z*_ij = minargv(z_ij )(c_ij) *)
let find_meds2 ?(keep_all_meds=false) (meds1 : meds_t) (meds2 : meds_t) =
    let find_exact () = 
        let best_meds = List.fold_left  
            (fun best_meds med1 -> 
                 List.fold_left  
                     (fun best_meds med2 -> 
                          let cost, (recost1, recost2), med_ls =
                              GenomeAli.find_med2_ls med1 med2 meds1.c2
                                  meds1.chrom_pam 
                          in  
                          if cost < best_meds.total_cost then  
                              { best_meds with med_ls = med_ls; 
                                    total_cost = cost; total_recost = (recost1 + recost2)}
                          else best_meds                      
                     ) best_meds meds2.med_ls
            ) {meds1 with med_ls = []; total_cost = max_int} meds1.med_ls
        in 
                                
        match keep_all_meds with 
        | true -> best_meds 
        | false ->             
              let kept_med_ls = keep meds1.chrom_pam best_meds.med_ls in  
              {best_meds with med_ls = kept_med_ls} 
    in 

    match meds1.chrom_pam.Data.approx with 
    | Some approx ->
          if approx then begin 
              update_cost_mat meds1 meds2;
              {meds1 with total_cost = meds1.approx_cost_arr.(meds2.code); 
                   total_recost = meds1.approx_recost_arr.(meds2.code)}
          end else find_exact ()

    | None -> find_exact ()
      


(** [find_meds3 medsp meds1 meds2] creates the median list
 * of three lists of medians [medsp=(x1,...,xk)], [meds1=(y1,...,yt)]
 * and [meds2=(z1,...,zq)] where xi, yj, and zp are medians. 
 * For each triplet (xi, yj, zp) we have 
 * a list of medians w_ijp with the same cost c_ijp. 
 * Find w*ijp = minargv_(w_ijp) (c_ijp) *)
let find_meds3 (medsp: meds_t) (meds1: meds_t) (meds2: meds_t) =    
    let meds1p = find_meds2 meds1 medsp in  
    let meds2p = find_meds2 meds2 medsp in  
    if meds1p.total_cost < meds2p.total_cost then meds1p
    else meds2p
            
       
(** [cmp_min_pair_cost] computes the min median cost
 * between two lists of medians [meds1=(x1,...,xk)] and [meds2=(y1,...,yt)]
 * where xi and yj are medians. For each pair (xi, yj) we have 
 * a list of medians z_ij with the same cost c_ij. 
 * returns c*_ij = min (c_ij) *)
let cmp_min_pair_cost (meds1 : meds_t) (meds2 : meds_t) = 
    let cmp_exact () = 
          List.fold_left 
              (fun (best_cost, recost) med1 ->
                   List.fold_left 
                       (fun (best_cost, recost) med2 ->
                            let acost, (a_recost1, a_recost2) = GenomeAli.cmp_cost med1 med2 meds1.c2 meds1.chrom_pam in 
                            if acost < best_cost then acost, a_recost1 + a_recost2
                            else best_cost, recost
                       ) (best_cost, recost) meds2.med_ls
              ) (max_int, 0) meds1.med_ls
    in 

    match meds1.chrom_pam.Data.approx with 
    | Some approx ->
          if approx then begin 
              update_cost_mat meds1 meds2;
              meds1.approx_cost_arr.(meds2.code), meds1.approx_recost_arr.(meds2.code)
          end else cmp_exact ()
    | None -> cmp_exact ()


(** [cmp_max_pair_cost] computes the max median cost
 * between two lists of medians [meds1=(x1,...,xk)] and [meds2=(y1,...,yt)]
 * where xi and yj are medians. For each pair (xi, yj) we have 
 * a list of medians z_ij with the same cost c_ij. 
 * returns c*_ij = max (c_ij) *)
let cmp_max_pair_cost (meds1 : meds_t) (meds2 : meds_t) = 
    let cmp_exact () = 
          List.fold_left 
              (fun (max_cost, max_recost) med1 ->
                   List.fold_left 
                       (fun (max_cost2, max_recost2) med2 ->
                            let cost, (recost1, recost2) = GenomeAli.cmp_cost med1 med2 meds1.c2 meds1.chrom_pam in 
                            if cost > max_cost2 then cost, recost1 + recost2
                            else max_cost2, max_recost2
                       ) (max_cost, max_recost) meds2.med_ls
              ) (0, 0) meds1.med_ls
    in 

    match meds1.chrom_pam.Data.approx with 
    | Some approx ->
          if approx then begin 
              update_cost_mat meds1 meds2;
              meds1.approx_cost_arr.(meds2.code), meds1.approx_recost_arr.(meds2.code)
          end else cmp_exact ()

    | None -> cmp_exact ()




(** [Compare meds1 meds2] returns 0 if these 
* two lists [meds1] and [meds2] are the same, 
* otherwise (-1) or 1 *)
let compare (meds1 : meds_t) (meds2 : meds_t) = 
    let num_med1 = List.length meds1.med_ls in 
    let num_med2 = List.length meds2.med_ls in 
    if num_med1 != num_med2 then num_med1 - num_med2
    else begin
        let rec check cur_meds1 =
            match cur_meds1 with    
            | [] -> 0 
            | med1::tl ->  
                  if List.exists 
                      (fun med2 -> (GenomeAli.compare med1 med2) = 0
                      ) meds2.med_ls then check tl 
                  else 1 
        in
        check meds1.med_ls              
    end


let clean_median meds1 meds2 = 
    meds1, meds2




(** [get_active_ref_code meds] returns active 
* reference codes of genome medians [meds] *)
let get_active_ref_code meds = 
    let med = List.hd meds.med_ls in
    med.GenomeAli.genome_ref_code, med.GenomeAli.genome_ref_code1, med.GenomeAli.genome_ref_code2



(** [copy_chrom_map s_ch d_ch] copies genome map
* from median [s_ch] into median [d_ch] *)
let copy_chrom_map s_ch d_ch =
    let copied_med_ls = List.map (fun ad_med -> 
                  let as_med = List.find 
                      (fun as_med -> as_med.GenomeAli.genome_ref_code1 = ad_med.GenomeAli.genome_ref_code  || 
                              as_med.GenomeAli.genome_ref_code2 = ad_med.GenomeAli.genome_ref_code
                      ) s_ch.med_ls
                  in 
                  GenomeAli.copy_chrom_map ad_med as_med
             ) d_ch.med_ls
    in 
    {d_ch with med_ls = copied_med_ls}





(** [readjust_3d ch1 ch2 mine c2 c3 parent] readjusts
* the current median [mine] of three medians [ch1],
* [ch2], and [parent] using three dimentional alignments*)
let readjust_3d ch1 ch2 mine c2 c3 parent = 
    let chrom_pam = mine.chrom_pam in 

    let ach1 = List.hd ch1.med_ls in 
    let ach2 = List.hd ch2.med_ls in
    let amine = List.hd mine.med_ls in
    let aparent = List.hd parent.med_ls in
    
    let cost1, _  = cmp_min_pair_cost ch1 mine in
    let cost2, _  = cmp_min_pair_cost ch2 mine in
    let costp, _ = cmp_min_pair_cost parent mine in
    let old_cost = cost1 + cost2 + costp in 


    let cost, adjust_med = GenomeAli.find_med3 ach1 ach2 aparent amine c2 c3
        chrom_pam in 
    let adjust_med = {mine with med_ls = [adjust_med]} in 

(*    fprintf stdout "old_cost: %i, new_cost: %i\n" old_cost cost; flush stdout; *)
    if old_cost <= cost then  old_cost, mine, false
    else cost, adjust_med, true
