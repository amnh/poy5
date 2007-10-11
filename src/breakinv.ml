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
let () = SadmanOutput.register "Breakinv" "$Revision: 911 $"

(** Median module contains functions to create medians
    between two lists of general chracters *)

let fprintf = Printf.fprintf

type breakinv_t = BreakinvAli.breakinv_t

(** A list medians between two sequences of general characters. 
 * Rearrangements are allowed *)
type meds_t = {
    med_ls : breakinv_t list;
    num_med : int;
(** total_cost = editing cost + rearrangement cost *)
    total_cost : int;   
    total_recost : int;
    breakinv_pam : Data.dyna_pam_t;    
    gen_cost_mat : Cost_matrix.Two_D.m;
    pure_gen_cost_mat : int array array;
    alpha : Alphabet.a
}


let init_med (seq : Sequence.s) gen_cost_mat alpha breakinv_pam = 
(*    Alphabet.print alpha;
    Cost_matrix.Two_D.output stdout gen_cost_mat;
*)
    {med_ls = [(BreakinvAli.init seq)];
     num_med = 1;   
     total_cost = 0;  
     total_recost = 0;  
     breakinv_pam =  breakinv_pam; 
     gen_cost_mat = gen_cost_mat;  
     pure_gen_cost_mat = Cost_matrix.Two_D.get_pure_cost_mat gen_cost_mat;  
     alpha = alpha;}  


(** Given a list of medians, determine a subset of medians to be kept to process
    further based on the customs's defined paramaters *)
let rec keep chrom_pam med_ls = 
    match chrom_pam.Data.keep_median with 
    | None -> med_ls 
    | Some keep_median ->
          if  keep_median >= List.length med_ls then med_ls
          else Utl.get_k_random_elem med_ls keep_median


(** Given two lists of medians [meds1=(x1,...,xk)] and [meds2=(y1,...,yt)]
 * where xi and yj are medians. For each pair (xi, yj) we have 
 * a list of medians z_ij with the same cost c_ij. 
 * Find z*_ij = minargv(z_ij )(c_ij) *)
let find_meds2 (meds1 : meds_t) (meds2 : meds_t) = 
    let update med1 med2 (best_meds : meds_t) :
            meds_t = 
        let cost, (recost1, recost2), med_ls =   
            BreakinvAli.find_med2_ls med1 med2 meds1.gen_cost_mat  
                meds1.pure_gen_cost_mat  
                meds1.alpha meds1.breakinv_pam   
        in   

        if cost < best_meds.total_cost then   
            {best_meds with  
                 total_cost = cost;  
                 total_recost = recost1 + recost2;
                 med_ls = med_ls;   
                 num_med = List.length med_ls} 
        else best_meds
    in 
                            
        
    let init_meds : meds_t = {
        med_ls = []; num_med = 0; 
         total_cost = Utl.infinity; 
         total_recost = 0;
         breakinv_pam = meds1.breakinv_pam;
         gen_cost_mat = meds1.gen_cost_mat;
         pure_gen_cost_mat = Cost_matrix.Two_D.get_pure_cost_mat meds1.gen_cost_mat;
         alpha = meds1.alpha} 
    in 

    let best_meds = 
        List.fold_left (fun best_meds1 med1 ->                                
                            List.fold_left (fun best_meds2 med2 ->
                                 update med1 med2 best_meds2 
                            ) best_meds1 meds2.med_ls 
                   ) init_meds meds1.med_ls
    in 
             
    let kept_med_ls = keep meds1.breakinv_pam best_meds.med_ls in
    {best_meds with  
         med_ls = kept_med_ls;         
         num_med = List.length kept_med_ls} 

    
(** Given two lists of medians [meds1=(x1,...,xk)] and [meds2=(y1,...,yt)]
 * where xi and yj are medians. For each pair (xi, yj) we have 
 * a list of medians z_ij with the same cost c_ij. 
 * Find c*_ij = min (c_ij) *)
let cmp_min_pair_cost (meds1 : meds_t) (meds2 : meds_t) =
    let min_cost, min_recost = List.fold_left 
        (fun (min_cost, min_recost) med1 -> 
                List.fold_left 
                    (fun (min_cost2, min_recost2) med2 -> 
                         let cost, (recost1, recost2) = BreakinvAli.cmp_cost med1 med2
                                meds1.gen_cost_mat 
                                meds1.pure_gen_cost_mat 
                                meds1.alpha meds1.breakinv_pam 
                         in  
                         if  min_cost2 > cost then cost, (recost1 + recost2)
                         else min_cost2, min_recost2
                    ) (min_cost, min_recost) meds2.med_ls
        ) (Utl.infinity, 0) meds1.med_ls 
    in 
    min_cost, min_recost


(** Given two lists of medians [meds1=(x1,...,xk)] and [meds2=(y1,...,yt)]
 * where xi and yj are medians. For each pair (xi, yj) we have 
 * a list of medians z_ij with the same cost c_ij. 
 * Find c*_ij = max (c_ij) *)
let cmp_max_pair_cost (meds1 : meds_t) (meds2 : meds_t) =
    let max_cost, max_recost = List.fold_left 
        (fun (max_cost, max_recost) med1 -> 
                List.fold_left 
                    (fun (max_cost2, max_recost2)  med2 -> 
                         let cost, (recost1, recost2) = BreakinvAli.cmp_cost med1 med2
                             meds1.gen_cost_mat  meds1.pure_gen_cost_mat  
                             meds1.alpha meds1.breakinv_pam                                
                         in  
                         if max_cost2 < cost then cost, (recost1 + recost2)
                         else max_cost2, max_recost2

                    ) (max_cost, max_recost) meds2.med_ls 
        ) (0, 0) meds1.med_ls 
    in 
    max_cost, max_recost


(** Given three lists of medians [medsp=(x1,...,xk)], [meds1=(y1,...,yt)]
 * and [meds2=(z1,...,zq)] where xi, yj, and zp are medians. 
 * For each triplet (xi, yj, zp) we have 
 * a list of medians w_ijp with the same cost c_ijp. 
 * Find w*ijp = minargv_(w_ijp) (c_ijp) *)
let find_meds3 (medsp: meds_t) (meds1: meds_t) (meds2: meds_t) =
    let meds1p = find_meds2 meds1 medsp in 
    let meds2p = find_meds2 meds2 medsp in 
    if meds1p.total_cost < meds2p.total_cost then meds1p
    else meds2p
            

let readjust_3d ch1 ch2 mine c2 c3 parent = 
    let seq1 = (List.hd ch1.med_ls).BreakinvAli.seq in
    let seq2 = (List.hd ch2.med_ls).BreakinvAli.seq in
    let seq3 = (List.hd parent.med_ls).BreakinvAli.seq in

    let ali_pam = BreakinvAli.get_breakinv_pam ch1.breakinv_pam in          

    let adjust_seq, cost = GenAli.create_gen_ali3 seq1 seq2 seq3 
        ch1.pure_gen_cost_mat ch1.alpha ali_pam.BreakinvAli.re_meth
        ali_pam.BreakinvAli.swap_med ali_pam.BreakinvAli.circular
    in 
    let amed = List.hd mine.med_ls in
    let adjust_med = {amed with BreakinvAli.seq = adjust_seq} in 
    let adjust_mine = {mine with med_ls = [adjust_med]} in 


    let cost1, _ = cmp_min_pair_cost ch1 mine in 
    let cost2, _ = cmp_min_pair_cost ch2 mine in 
    let cost3, _ = cmp_min_pair_cost parent mine in 
    let old_cost = cost1 + cost2 + cost3 in 
(*    fprintf stdout "Cost: %i  -> adjusted cost: %i\n" old_cost cost;*)
    if old_cost < cost then old_cost, mine, false
    else cost, adjust_mine, 0 <> (compare mine adjust_mine)

       
    

(** Compare two list of medians *)
let compare (meds1 : meds_t) (meds2 : meds_t) = 
    let num_med1 = meds1.num_med in 
    let num_med2 = meds2.num_med in 
    match num_med1 != num_med2 with
    | true -> num_med1 - num_med2  
    | false -> 
        let rec compare_meds med1_ls med2_ls = 
            match med1_ls, med2_ls with 
                | med1::tl1, med2::tl2 ->
                    if compare med1 med2 = 0 then compare_meds tl1 tl2
                    else compare med1 med2
                | _ -> 0
        in  


        compare_meds meds1.med_ls meds2.med_ls

    

let get_active_ref_code meds = 
    let med = List.hd meds.med_ls in
    med.BreakinvAli.ref_code, med.BreakinvAli.ref_code1, med.BreakinvAli.ref_code2

