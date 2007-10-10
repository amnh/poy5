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

let () = SadmanOutput.register "Chrom" "$Revision: 975 $"
(** Median module contains functions to create medians
    between two lists of chromosomes *)

let fprintf = Printf.fprintf


type med_t = ChromAli.med_t

type meds_t = {
    med_ls : med_t list;
    total_cost : int;   
    total_recost : int;
    c2 : Cost_matrix.Two_D.m; 

    approx_cost_arr : int array;
    approx_med_arr : med_t array;
    approx_recost_arr : int array;
    code : int;

    chrom_pam : Data.dyna_pam_t 
}

let max_taxa_id = ref 0 

let init_med (seq : Sequence.s) c2 chrom_pam tcode num_taxa =  

    let seq = 
        if Sequence.is_empty seq (Alphabet.get_gap Alphabet.nucleotides) then
            UtlPoy.get_empty_seq ()
        else seq
    in

    let med = ChromAli.init_med seq in
    max_taxa_id := max !max_taxa_id num_taxa;

    { 
        med_ls = [med]; 
    
        total_cost = 0; 
        total_recost = 0;
        c2 = c2;
        approx_med_arr = (Array.make !max_taxa_id med);
        approx_cost_arr = (Array.make !max_taxa_id max_int);
        approx_recost_arr = (Array.make !max_taxa_id max_int);
        code = tcode;

        chrom_pam = chrom_pam 
    } 

let print meds =
       List.iter (fun med -> fprintf stdout "%i -> %i %i\n " med.ChromAli.ref_code
                      med.ChromAli.ref_code1 med.ChromAli.ref_code2;
                      ChromAli.print_map med.ChromAli.chrom_map;
                 ) meds.med_ls;
    flush stdout


let update_approx_mat meds1 meds2 =     
    let med1 = List.hd meds1.med_ls in  
    let med2 = List.hd meds2.med_ls in  
    let code2 = meds2.code in
    (if meds1.approx_cost_arr.(code2) = max_int then begin 
        let cost, recost, med2_ls = ChromAli.find_med2_ls med1 med2 meds1.c2 meds1.chrom_pam in  
        meds1.approx_med_arr.(code2) <- List.hd med2_ls;
        meds1.approx_cost_arr.(code2) <- cost; 
        meds1.approx_recost_arr.(code2) <- recost; 
     end) 

(** Compute the minimum pair cost between medians in two lists *)
let cmp_min_pair_cost (meds1 : meds_t) (meds2 : meds_t) = 
    let cmp_exact () = 
          List.fold_left 
              (fun (min_cost, min_recost) med1 ->
                   List.fold_left 
                       (fun (min_cost2, min_recost2) med2 ->
                            let cost, recost = ChromAli.cmp_cost med1 med2
                                meds1.c2 meds1.chrom_pam `Chromosome 
                            in 
                            if cost < min_cost2 then cost, recost
                            else min_cost2, min_recost2
                       ) (min_cost, min_recost) meds2.med_ls
              ) (max_int, 0) meds1.med_ls
    in 

    match meds1.chrom_pam.Data.approx with 
    | Some approx ->
          if approx then begin 
              update_approx_mat meds1 meds2;
              meds1.approx_cost_arr.(meds2.code), meds1.approx_recost_arr.(meds2.code)
          end else cmp_exact ()

    | None -> cmp_exact ()

(** Compute the maximum pair cost between medians in two lists *)
let cmp_max_pair_cost (meds1 : meds_t) (meds2 : meds_t) = 
    let cmp_exact () = 
          List.fold_left 
              (fun (max_cost, max_recost) med1 ->
                   List.fold_left 
                       (fun (max_cost2, max_recost2) med2 ->
                            let cost, recost = ChromAli.cmp_cost med1 med2
                                meds1.c2 meds1.chrom_pam `Chromosome
                            in 
                            if cost > max_cost2 then cost, recost
                            else max_cost2, max_recost2
                       ) (max_cost, max_recost) meds2.med_ls
              ) (0, 0) meds1.med_ls
    in 

    match meds1.chrom_pam.Data.approx with 
    | Some approx ->
          if approx then begin 
              update_approx_mat meds1 meds2;
              meds1.approx_cost_arr.(meds2.code), meds1.approx_recost_arr.(meds2.code)
          end else cmp_exact ()

    | None -> cmp_exact ()


(** Find the best median list between two lists of medians *)
let find_meds2 (meds1 : meds_t) (meds2 : meds_t) =
    let find_exact () = 
        let best_meds = List.fold_left  
            (fun best_meds med1 -> 
                 List.fold_left  
                     (fun best_meds med2 -> 
                          let cost, recost, med_ls =
                              ChromAli.find_med2_ls med1 med2 meds1.c2 meds1.chrom_pam in  
    
                        if cost < best_meds.total_cost then 
                            { best_meds with med_ls = med_ls; total_cost = cost;
                                  total_recost = recost}
                        else best_meds                      
                     ) best_meds meds2.med_ls
            ) {meds1 with med_ls = []; total_cost = max_int} meds1.med_ls
        in 
                                
        best_meds
    in 


    match meds1.chrom_pam.Data.approx with 
    | Some approx ->
          if approx then begin 
              update_approx_mat meds1 meds2;

              let med1 = List.hd meds1.med_ls in  
              let med2 = List.hd meds2.med_ls in  
              let code2 = meds2.code in

              let med12 = ChromAli.find_approx_med2 med1 med2
                  meds1.approx_med_arr.(code2) 
              in 
              {meds1 with med_ls = [med12]; 
                   total_cost = meds1.approx_cost_arr.(code2);
                   total_recost = meds1.approx_recost_arr.(code2)}
                   
          end else find_exact ()

    | None -> find_exact ()
      

(** Find the median between three meds lists*)
let find_meds3 (medsp: meds_t) (meds1: meds_t) (meds2: meds_t) =    
    let meds1p = find_meds2 meds1 medsp in  
    let meds2p = find_meds2 meds2 medsp in  
    if meds1p.total_cost < meds2p.total_cost then meds1p
    else meds2p
            

let readjust_3d ch1 ch2 mine c2 c3 parent = 
    let adjust_med = find_meds3 parent ch1 ch2 in 
    let amed = List.hd adjust_med.med_ls in
    let single_seq = UtlPoy.get_single_seq amed.ChromAli.seq c2 in 
    
    
    let adjust_med = {adjust_med with med_ls = [{amed with ChromAli.seq = single_seq}]} in 
    let cost1, _  = cmp_min_pair_cost ch1 adjust_med in
    let cost2, _  = cmp_min_pair_cost ch2 adjust_med in
    let costp, _ = cmp_min_pair_cost parent adjust_med in
    (cost1 + cost2 + costp), adjust_med, 0 <> (compare mine adjust_med)


       


(** Check if median list meds1 is the same median list meds2 *)
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
                      (fun med2 -> (Sequence.compare med1.ChromAli.seq med2.ChromAli.seq) = 0
                      ) meds2.med_ls then check tl 
                  else 1 
        in
        check meds1.med_ls              
    end
    

(** ============================================================== **)
let get_active_ref_code meds = 
(*
    List.iter ChromAli.print  meds.med_ls;
    flush stdout;  
*)
  
    let med = List.hd meds.med_ls in
    med.ChromAli.ref_code, med.ChromAli.ref_code1, med.ChromAli.ref_code2

let copy_chrom_map s_ch d_ch =
    let copied_med_ls = List.map (fun ad_med -> 
                  let as_med = List.find 
                      (fun as_med -> as_med.ChromAli.ref_code1 = ad_med.ChromAli.ref_code  || 
                              as_med.ChromAli.ref_code2 = ad_med.ChromAli.ref_code
                      ) s_ch.med_ls
                  in 
                  ChromAli.copy_chrom_map as_med ad_med
             ) d_ch.med_ls
    in 
    {d_ch with med_ls = copied_med_ls}
