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

(** Breakinv module contains functions to create medians
*    between two lists of breakinv chracters *)

let debug = false 

let fprintf = Printf.fprintf

type breakinv_t = BreakinvAli.breakinv_t

(** [meds_t] is a data structure for a list medians 
* between two breakinv character lists *)
type meds_t = {
    med_ls : breakinv_t list; (** breakinv list *)
    num_med : int; (** number of breakinv characters *)
    total_cost : int;   (** the cost to create this breakinv list *)
    total_recost : int; (** the recost to create this breakinv list *)
    breakinv_pam : Data.dyna_pam_t; (** breakinv paramenters used to create breakinv median *)
    gen_cost_mat : Cost_matrix.Two_D.m; 
    pure_gen_cost_mat : int array array;
    alpha : Alphabet.a 
}


(** [init_med seq gen_cost_mat alpha breakinv_pam] returns
* a breakinv character list with only one element 
* created from a sequence of general character [seq]*)
let init_med (seq : Sequence.s) gen_cost_mat alpha breakinv_pam = 
    {med_ls = [(BreakinvAli.init seq [])];
     num_med = 1;   
     total_cost = 0;  
     total_recost = 0;  
     breakinv_pam =  breakinv_pam; 
     gen_cost_mat = gen_cost_mat;  
     pure_gen_cost_mat = Cost_matrix.Two_D.get_pure_cost_mat gen_cost_mat;  
     alpha = alpha;}  


(** [keep chrom_pam med_ls] returns a sublist of median list
* [med_ls] to be kept to process further based on the customs's defined paramaters *)
let rec keep chrom_pam med_ls = 
    match chrom_pam.Data.keep_median with 
    | None -> med_ls 
    | Some keep_median ->
          if  keep_median >= List.length med_ls then med_ls
          else Utl.get_k_random_elem med_ls keep_median

let transform_matrix (inmat : int array array) invcost =
    let sizei = Array.length inmat 
    and sizej = Array.length inmat.(0) in
    let isodd x = if (x mod 2)<>0 then true else false in
    let iseven x = if (x mod 2)=0 then true else false in
    (* the first line and column of inmat is not being used. 
    *  the last line and column of inmat is for gap, no need 
    *  to modify cost there*)
    Array.mapi (fun i arri ->
        if (i>=1)&&(i<=sizei-2) then  
            Array.mapi ( fun j cost ->
                if (j>=1)&&(j<=sizej-2) then
                    if (isodd i)&&(isodd j) then cost
                    else if (iseven i)&&(iseven j) then cost+2*invcost
                    else cost+invcost
                else cost
            ) arri
        else arri
    ) inmat

(** [find_meds2 meds1 meds2] returns a list of 
* breakinv character medians created for two lists of medians 
* [meds1=(x1,...,xk)] and [meds2=(y1,...,yt)]
* where xi and yj are medians. For each pair (xi, yj) we have 
* a list of medians z_ij with the same cost c_ij. 
* Find z*_ij = minargv(z_ij )(c_ij) *)
let find_meds2 (meds1 : meds_t) (meds2 : meds_t) = 
    let update med1 med2 (best_meds : meds_t) : meds_t =
        let bkpam = meds1.breakinv_pam in
        let pure_gen_cost_mat =
            match bkpam.Data.re_meth with
            | Some re_meth ->
                (match re_meth with
                | `Locus_Breakpoint c -> meds1.pure_gen_cost_mat
                | `Locus_Inversion invc -> 
                    transform_matrix meds1.pure_gen_cost_mat invc
                )
            | None -> meds1.pure_gen_cost_mat
        in
        let cost, (recost1, recost2), med_ls =   
            BreakinvAli.find_med2_ls med1 med2 meds1.gen_cost_mat  
                pure_gen_cost_mat meds1.alpha bkpam   
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
         total_cost = Utl.large_int; 
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

    
(** [cmp_min_pair_cost] returns the minimum cost
* between two lists of medians [meds1=(x1,...,xk)] and [meds2=(y1,...,yt)]
 * where xi and yj are medians. For each pair (xi, yj) we have 
 * a list of medians z_ij with the same cost c_ij. 
 * Find c*_ij = min (c_ij) *)
let cmp_min_pair_cost (meds1 : meds_t) (meds2 : meds_t) =
    let bkpam = meds1.breakinv_pam in
    let pure_gen_cost_mat =
            match bkpam.Data.re_meth with
            | Some re_meth ->
                (match re_meth with
                | `Locus_Breakpoint c -> meds1.pure_gen_cost_mat
                | `Locus_Inversion invc -> 
                    transform_matrix meds1.pure_gen_cost_mat invc
                )
            | None -> meds1.pure_gen_cost_mat
    in
    let min_cost, min_recost = List.fold_left 
        (fun (min_cost, min_recost) med1 -> 
                List.fold_left 
                    (fun (min_cost2, min_recost2) med2 -> 
                         let cost, (recost1, recost2) = BreakinvAli.cmp_cost med1 med2
                                meds1.gen_cost_mat 
                                pure_gen_cost_mat 
                                meds1.alpha bkpam 
                         in  
                         if  min_cost2 > cost then cost, (recost1 + recost2)
                         else min_cost2, min_recost2
                    ) (min_cost, min_recost) meds2.med_ls
        ) (Utl.large_int, 0) meds1.med_ls 
    in 
    min_cost, min_recost


(** [cmp_max_pair_cost] returns the maximum cost
* between two lists of medians [meds1=(x1,...,xk)] and [meds2=(y1,...,yt)]
 * where xi and yj are medians. For each pair (xi, yj) we have 
 * a list of medians z_ij with the same cost c_ij. 
 * Find c*_ij = min (c_ij) *)
let cmp_max_pair_cost (meds1 : meds_t) (meds2 : meds_t) =
    let bkpam = meds1.breakinv_pam in
    let pure_gen_cost_mat =
            match bkpam.Data.re_meth with
            | Some re_meth ->
                (match re_meth with
                | `Locus_Breakpoint c -> meds1.pure_gen_cost_mat
                | `Locus_Inversion invc -> 
                    transform_matrix meds1.pure_gen_cost_mat invc
                )
            | None -> meds1.pure_gen_cost_mat
    in
    let max_cost, max_recost = List.fold_left 
        (fun (max_cost, max_recost) med1 -> 
                List.fold_left 
                    (fun (max_cost2, max_recost2)  med2 -> 
                         let cost, (recost1, recost2) = BreakinvAli.cmp_cost med1 med2
                             meds1.gen_cost_mat pure_gen_cost_mat  
                             meds1.alpha bkpam  in  
                         if max_cost2 < cost then cost, (recost1 + recost2)
                         else max_cost2, max_recost2

                    ) (max_cost, max_recost) meds2.med_ls 
        ) (0, 0) meds1.med_ls 
    in 
    max_cost, max_recost


(** [find_meds3 medsp meds1 meds2] returns 
* the median of three lists of medians [medsp=(x1,...,xk)], [meds1=(y1,...,yt)]
 * and [meds2=(z1,...,zq)] where xi, yj, and zp are medians. 
 * For each triplet (xi, yj, zp) we have 
 * a list of medians w_ijp with the same cost c_ijp. 
 * Find w*ijp = minargv_(w_ijp) (c_ijp) *)
let find_meds3 (medsp: meds_t) (meds1: meds_t) (meds2: meds_t) =
    let meds1p = find_meds2 meds1 medsp in 
    let meds2p = find_meds2 meds2 medsp in 
    if meds1p.total_cost < meds2p.total_cost then meds1p
    else meds2p

            
(** [readjust_3d ch1 ch2 mine c2 c3 parent] readjusts
* the breakinv median [mine] of three breakinv medians 
* [ch1], [ch2] and [parent] *) 
let readjust_3d ch1 ch2 mine c2 c3 parent =
    let ali_pam = BreakinvAli.get_breakinv_pam ch1.breakinv_pam in
    let cost1, _ = cmp_min_pair_cost ch1 mine in 
    let cost2, _ = cmp_min_pair_cost ch2 mine in 
    let cost3, _ = cmp_min_pair_cost parent mine in 
    let old_cost : int = cost1 + cost2 + cost3 in 
    let seq1 = (List.hd ch1.med_ls).BreakinvAli.seq in
    let seq2 = (List.hd ch2.med_ls).BreakinvAli.seq in
    let seq3 = (List.hd parent.med_ls).BreakinvAli.seq in
    let mine_seq = (List.hd mine.med_ls).BreakinvAli.seq in
    let delimiters1 = (List.hd ch1.med_ls).BreakinvAli.delimiter_lst in 
    let delimiters2 = (List.hd ch2.med_ls).BreakinvAli.delimiter_lst in 
    let delimiters3 = (List.hd parent.med_ls).BreakinvAli.delimiter_lst in 
    (* no need to keep 0 in delimiter list*)
    let delimiters1 = List.filter (fun x -> x<>0 ) delimiters1 in
    let delimiters2 = List.filter (fun x -> x<>0 ) delimiters2 in
    let delimiters3 = List.filter (fun x -> x<>0 ) delimiters3 in
    (* debug msg 
    Printf.printf "### Breakinv.readjust_3d compute the old cost first,check seq1/seq2/seq3/mine_seq ###\n%!";
    Sequence.printseqcode seq1; Sequence.printseqcode seq2;
    Sequence.printseqcode seq3; Sequence.printseqcode mine_seq;
    Printf.printf "old~cost1,cost2,cost3=%d,%d,%d;\n%!" cost1 cost2 cost3;
    Printf.printf "deli = [%!";
    List.iter (Printf.printf "%d,") delimiters1;
    Printf.printf "] [%!";
    List.iter (Printf.printf "%d,") delimiters2;
    Printf.printf "] [%!";
    List.iter (Printf.printf "%d,") delimiters3;
    Printf.printf "]\n%!";
    debug msg*)
    let  adjust_seq, adjust_delimiters, new_cost1, new_cost2, new_cost3 =
    if (cost1=0)&&(cost2=0)&&(cost3=0) then
        mine_seq,[],cost1,cost2,cost3
    else
        match ali_pam.BreakinvAli.median_solver with
        | `Vinh ->
            if debug then Printf.printf "Vinh median solver\n%!";
            let adjust_seq, cost1,cost2,cost3 = 
                GenAli.create_gen_ali3 ali_pam.BreakinvAli.kept_wag 
                seq1 seq2 seq3 mine_seq 
                ch1.pure_gen_cost_mat 
                ch1.alpha 
                ali_pam.BreakinvAli.re_meth
                ali_pam.BreakinvAli.swap_med 
                ali_pam.BreakinvAli.circular
                (Alphabet.get_orientation ch1.alpha) 
                ali_pam.BreakinvAli.symmetric 
            in
            if debug then begin
                Printf.printf "median solver: new cost=%d/%d/%d\n%!" cost1 cost2 cost3;
                Printf.printf "adjustseq = %!";
                Sequence.printseqcode adjust_seq;
            end;
            adjust_seq,[],cost1,cost2,cost3
        | _ ->
            if debug then Printf.printf "Grappa median solver\n%!";
            let adjust_seq,delimiter_lst,cost1,cost2,cost3,_,_,_,_,_,_,_,_,_ = 
                GenAli.create_gen_ali3_by_medsov 
                ali_pam.BreakinvAli.median_solver
                ali_pam.BreakinvAli.kept_wag
                seq1 seq2 seq3 
                [delimiters1;delimiters2;delimiters3]
                ch1.gen_cost_mat
                ch1.alpha 
                ali_pam.BreakinvAli.re_meth
                ali_pam.BreakinvAli.swap_med 
                ali_pam.BreakinvAli.circular
                (Alphabet.get_orientation ch1.alpha) 
                ali_pam.BreakinvAli.symmetric 
            in
            adjust_seq,delimiter_lst,cost1,cost2,cost3
    in
    let new_cost = new_cost1+new_cost2+new_cost3 in
    if old_cost <= new_cost then 
        cost3, mine, false
    else 
        let amed = List.hd mine.med_ls in
        let adjust_med_ls = 
        {amed with BreakinvAli.seq = adjust_seq
        } in
        let adjust_med_ls = 
        match adjust_delimiters with
        | [] -> adjust_med_ls
        | _ -> 
            {adjust_med_ls with BreakinvAli.delimiter_lst = adjust_delimiters}
        in
        new_cost3, {mine with med_ls = [adjust_med_ls]}, true
       
    

(** [compare meds1 meds2] returns 0 if breakinv list [meds1]
* is the same as breakinv list [meds2], otherwise (-1) or (1) *)
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

    
(** [get_active_ref_code meds] return active reference codes
* of breakinv medians [meds] *)
let get_active_ref_code meds = 
    let med = List.hd meds.med_ls in
    med.BreakinvAli.ref_code, med.BreakinvAli.ref_code1, med.BreakinvAli.ref_code2

let update_medst old_medst newseq delimiters =
    let new_bkinv_t = 
        BreakinvAli.update_bkinv_t (List.hd old_medst.med_ls) newseq delimiters
    in
    {old_medst with med_ls = [new_bkinv_t] }

let single_to_multi medst = 
    let medls = List.hd medst.med_ls in
    let (medls_lst:breakinv_t list) = BreakinvAli.single_to_multi medls in
    List.map (fun bkinvt_hd ->
        { medst with med_ls = [bkinvt_hd] }
    )medls_lst
        

