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

(** A Breakinv character set implementation. The breakinv
* character set allows rearrangements *)

let () = SadmanOutput.register "BreakinvCS" "$Revision: 2754 $"

exception Illegal_Arguments
let () = SadmanOutput.register "Breakinv Character" "$Revision: 2754 $"

let debug = false

let fprintf = Printf.fprintf

module IntMap = All_sets.IntegerMap
module IntSet = All_sets.Integers
type meds_t = Breakinv.meds_t

(** [t] is data structure to present a set of breakinv characters *)
type t = { 
    meds : meds_t IntMap.t; (** a set of breakinv characters *)
    costs : float IntMap.t;
    recosts : float IntMap.t;
    total_cost : float;    (** The total cost of the character set *)
    total_recost : float;  (** The total cost of the character set *)

    subtree_cost : float;         (** The total subtree cost of the character set *)
    subtree_recost : float;         (** The total subtree recost of the character set *)

    c2_full : Cost_matrix.Two_D.m;       (** The two dimensional cost matrix to be used in the character set *)
    c2_original : Cost_matrix.Two_D.m;       (** The two dimensional cost matrix to be used in the character set *)
    c3 : Cost_matrix.Three_D.m;     (** The three dimensional cost matrix to be used in the character set *)
    alph : Alphabet.a;              (** The alphabet of the sequence set *)
    breakinv_pam : Data.dyna_pam_t; 
    code : int;                     (** The set code. n = this is nth 
    chromosome of multichromosome genes, n = 1,2,3,...... ,
    also the key to hashtbl data.character_codes, the map between the character
    codes and their corresponding names. in data.ml, each chromosome is a character.*) 
}

let print thist =
    IntMap.iter (fun _ med -> Breakinv.print med) thist.meds

let cardinal x = IntMap.fold (fun _ _ acc -> acc + 1) x.meds 0

(** [of_array spec arr code] creates a breakinv set
* from an array of sequences [arr] *)
let of_array spec arr code = 
    let adder (meds, costs, recosts) (seq,delimiter,key) = 
        let med = 
            Breakinv.init_med seq delimiter spec.Data.tcm2d_full
            spec.Data.tcm2d_original spec.Data.alph spec.Data.pam 
        in 
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
        total_recost = 0.0;
        subtree_cost = 0.;
        subtree_recost = 0.;
        c2_full = spec.Data.tcm2d_full;
        c2_original = spec.Data.tcm2d_original;
        c3 = spec.Data.tcm3d;
        alph = spec.Data.alph;
        breakinv_pam = spec.Data.pam;
        code = code;
    }

(** [of_list spec lst code] create a breakinv set
* from a list sequence [lst] *)
let of_list spec lst code = 
    let arr = Array.of_list lst in
    of_array spec arr code

(** [to_list t] returns a list of breakinv characters 
* from a set of breakinv characters *)
let to_list t =
    IntMap.fold (fun code med acc -> (med, code) :: acc) t.meds []

(** [same_codes a b] returns true if code set of [a] 
* is the same to code set of [b], otherwise false *)
let same_codes a b =
    let checker x _ res = res && (IntMap.mem x b) in
    IntMap.fold checker a true

(** [median2 a b] returns the median set 
* between breakinv character sets [a] and [b] *)
let median2 (a : t) (b : t) =
    (* We will use imperative style for this function *)
    let debug = false in
    let empty = IntMap.empty in
    let median code (meda : meds_t) (medians, costs, recosts, total_cost, total_recost) = 
        let medb : meds_t = IntMap.find code b.meds in
        if debug then begin
        Printf.printf "breakinvCS.median2, meda is : %!";
        List.iter (fun x -> Sequence.printseqcode x.BreakinvAli.seq)
        meda.Breakinv.med_ls;
        Printf.printf "medb is : %!";
        List.iter (fun x -> Sequence.printseqcode x.BreakinvAli.seq)
        medb.Breakinv.med_ls;
        end;
        let medab = Breakinv.find_meds2 meda medb in
        if debug then begin
        Printf.printf "medab is : %!";
        List.iter (fun x -> Sequence.printseqcode x.BreakinvAli.seq;)
        medab.Breakinv.med_ls;
        Printf.printf "\n%!";
        end;
        let new_median = IntMap.add code medab medians 
        and new_costs = 
            IntMap.add code (float_of_int medab.Breakinv.total_cost) costs  
        and new_recosts = 
            IntMap.add code (float_of_int medab.Breakinv.total_recost) recosts  
        and new_total_cost = total_cost + medab.Breakinv.total_cost 
        and new_total_recost = total_recost + medab.Breakinv.total_recost in
        new_median, new_costs, new_recosts, new_total_cost, new_total_recost
    in
    let medab_map, new_costs, new_recosts, total_cost, total_recost = 
        IntMap.fold median a.meds (empty, empty, empty, 0, 0)
    in
    let subtree_recost = a.subtree_recost +. b.subtree_recost +. (float_of_int total_recost) in
    let subtree_cost = a.subtree_cost +. b.subtree_cost +. (float_of_int total_recost) in
    { a with meds = medab_map; costs = new_costs; recosts = new_recosts;
          total_cost = float_of_int total_cost; 
          total_recost = float_of_int total_recost;
          subtree_recost = subtree_recost;
          subtree_cost = subtree_cost;
    }
   
(** [median3 p n c1 c2] returns the median set of
* breakinv character sets [p], [c1] and [c2] *)
let median3 p n c1 c2 =
    let median code  medp res_medians = 
        let med1= IntMap.find code c1.meds in 
        let med2 = IntMap.find code c2.meds in
        let medp12 = 
               Breakinv.find_meds3 medp med1 med2 
        in
        (*debug msg
        let tmp_breakinv_t = (List.hd medp12.Breakinv.med_ls) in
        let tmpseq = tmp_breakinv_t.BreakinvAli.seq in
        let tmprecost1 = tmp_breakinv_t.BreakinvAli.recost1 in
        let tmprecost2 = tmp_breakinv_t.BreakinvAli.recost2 in
        Printf.printf "recost1/recost2 = %d/%d, seq = \n%!" tmprecost1 tmprecost2;
        Sequence.printseqcode tmpseq;
        debug msg*)
        IntMap.add code medp12 res_medians
    in 
    let acc = IntMap.empty in
    let medp12_map = IntMap.fold median p.meds acc in
    { n with meds = medp12_map; }

let flatten t_lst = 
    let parent = List.hd t_lst in (*we only have one item in t_lst to start with*)
    parent.code,
    IntMap.fold ( fun kn parent_medst (seq_lstlst:Sequence.s list list) ->
       let medst_lst = List.map (fun x ->  IntMap.find kn x.meds
       ) (List.tl t_lst) in
       let medls_lst = List.map (fun x -> x.Breakinv.med_ls
       ) medst_lst in
       let add_seqlst:Sequence.s list = 
           List.map (fun x -> 
               (List.hd x).BreakinvAli.seq 
               ) medls_lst 
       in
       let parent_seq = (List.hd (parent_medst.Breakinv.med_ls)).BreakinvAli.seq in
       let add_seqlst = parent_seq::add_seqlst in
       match seq_lstlst with
       | [[]] ->(*this is the first median we visit in the IntMap*)
              (*  [[ch1_seq];[ch2_seq];[parent_seq];[mine_seq]] *)
              List.map (fun x -> [x]) add_seqlst
       | _ ->  
           List.map2 
               (fun (addseq:Sequence.s) (seqlst:Sequence.s list) -> 
                addseq :: seqlst  
               ) add_seqlst seq_lstlst;
    ) parent.meds [[]]

(* is the nth of IntMap corresponding to nth of the sequence list? *)
let update_t oldt (file_median_seq:Sequence.s list list) (file_median_chrom_seqdeli: int list list list) =
    let i = ref 0 in
    List.map2 (fun median_seq median_chrom_seqdeli ->
    i := 0;
    let newmedsMap = 
      IntMap.mapi ( fun key old_meds_t ->
          (* Printf.printf "key is %d\n%!" key; when we have only one median of
          * each node, key is 1. what if we have more than one median?*)
          let res = Breakinv.update_medst old_meds_t (List.nth median_seq !i) (List.nth
          median_chrom_seqdeli !i) in
          i := !i +1 ;
          res
      ) oldt.meds
    in
    (*we only update median sequences and code here, I hope it's enough.
    * also I'm not sure about the code. this code used to be the key for hashtbl
    * in data.ml : character_codes. if we merge multi-chromosome of a node in a
    * file into a single chromsome with delimiters, the mapping between old keys
    * and their characters won't be correct.*)
    { oldt with meds = newmedsMap; code = !i }
    ) file_median_seq file_median_chrom_seqdeli

let single_to_multi single_t =
    let maxlen = ref 0 in
    let medianlst_single_mapset = to_list single_t in
    let medianlst_multi_mapset (*(Breakinv.meds_t list, IntMap.key) list *) = 
        List.map (fun (med,code) ->
        let medst_lst:meds_t list = Breakinv.single_to_multi med in
        let newlen = List.length medst_lst in
        if (newlen > (!maxlen)) then maxlen := newlen;
        (medst_lst,code) 
        ) medianlst_single_mapset in
    let maxlen = !maxlen in 
    let empty = IntMap.empty in
    let maplst = Array.to_list (Array.make maxlen empty) in
    let maplst = List.map (fun x -> ref x) maplst in
    let medlen = (List.length medianlst_multi_mapset) in
    for i = 0 to (List.length medianlst_multi_mapset)-1 do
        begin
            let (medst_lst,code) = List.nth medianlst_multi_mapset i in
            let len = List.length medst_lst in
            (* length of medst_lst maybe shorter than maplist,that's why we cannot use
            * map2 here *)
            if debug then
                Printf.printf "bkinvCS.ml single_to_multi,medianNUM=%d,code=%d,maxlen = %d,len=%d\n%!" medlen code maxlen len;
            for j = 0 to len-1 do
                let add_medst = List.nth medst_lst j in
                let add_to_map = List.nth maplst j in
                add_to_map := IntMap.add (j+1) add_medst !add_to_map
            done;
        end
    done;
    let count = ref 0 in
    let t_list = 
        List.map (fun intmap ->
        count := !count + 1;
        let intmap = !intmap in
        { single_t with meds = intmap; code = !count }
    ) maplst
    in
    t_list
    (* wait, what we do with costs?recosts?etc *)



(** [readjust to_adjust modified ch1 ch2 parent mine] returns
* the readjusted median set [mine] of three breakinv character
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
                Breakinv.readjust_3d ch1_chrom ch2_chrom my_chrom
                    c2 c3 parent_chrom
            in
            let new_single = IntMap.add code seqm res_medians
            and new_costs = IntMap.add code (float_of_int rescost) res_costs 
            and new_total = total + rescost in
            let new_totalsum = totalsum + ressumcost in
            let modified = 
                if changed then 
                    All_sets.Integers.add code modified
                else modified
            in
            modified, new_single, new_costs, new_total, new_totalsum        
        end 
    in 
    let modified, meds, costs, total_cost,total_sum_cost = 
        IntMap.fold adjusted parent.meds (modified, empty, empty, 0, 0)
    in
    let tc = float_of_int total_cost in
    let tsc = float_of_int total_sum_cost in
    let subtree_recost = tc +. ch1.subtree_recost +. ch2.subtree_recost in
    let subtree_cost = tsc in
    if debug then
        Printf.printf "ADJUST.... subtree_cost = %f,subtree_recost = %f+%f+%f = %f \n%!" 
        subtree_cost ch1.subtree_recost ch2.subtree_recost tc subtree_recost; 
    modified,
    tc,
    tsc,
    { mine with meds = meds; costs = costs; 
       total_cost = tc; 
       subtree_recost = subtree_recost; 
       subtree_cost = subtree_cost; 
        }

(**[get_extra_cost_for_root root] return the extra cost result from non-zero
* diagonal in cost matrix.*)
let get_extra_cost_for_root (a :t) =
    let get_ec code medst acc =
        acc + Breakinv.get_extra_cost_for_root medst 
    in
    let is_identity = Cost_matrix.Two_D.is_identity a.c2_original in
    if is_identity then 0. (*0 diagonal in cost matrix*)
    else (*non-0 diagonal in cost matrix*)
    float_of_int (IntMap.fold get_ec a.meds 0) 


(** [distance a b] returns total distance between 
* two breakinv character sets [a] and [b] *)
let distance (a : t) (b : t)  = 
    let single_distance code meda (acc_cost, acc_recost) =
        let medb = IntMap.find code b.meds in
        let cost, recost = Breakinv.cmp_min_pair_cost meda medb in 
        acc_cost + cost, acc_recost + recost
    in
    let cost, _ = IntMap.fold single_distance a.meds (0, 0) in 
    float_of_int cost

(** [max_distance a b] returns total maximum distances
* between two breakinv character sets [a] and [b]*)
let max_distance (a : t) (b : t)  = 
    let single_distance code meda (acc_cost, acc_recost) =
        let medb = IntMap.find code b.meds in
        let cost, recost = Breakinv.cmp_max_pair_cost meda medb in 
        acc_cost + cost, acc_recost + recost
    in
    let cost, _ = IntMap.fold single_distance a.meds (0,0) in 
    float_of_int cost

(** [to_string a] converts the breakinv character set [a]
* into string format *)
let to_string a =
    let builder code med acc =
        let code = string_of_int code in 
        let seq_ls = List.map (fun med -> Sequence.to_string med.BreakinvAli.seq a.alph) 
            med.Breakinv.med_ls 
        in 

        let seq = String.concat ":" seq_ls in 
        acc ^ code ^ ": " ^ seq ^ "; "
    in
    IntMap.fold builder a.meds ""

(* [dist_2 n a b] calculates the cost of joining 
* the node containing breakinv character set [n]
* between two nodes containing [a] and [b]. 
* [a] must be the parent (ancestor) of [b] *)
let dist_2 n a b =
    let cost_calculator code medb (acc_cost, acc_recost) =
        let medn = IntMap.find code n.meds
        and meda = IntMap.find code a.meds 
        and medb = IntMap.find code b.meds in
        
        let medab : meds_t = Breakinv.find_meds2 meda medb in 
        
        let cost, recost = Breakinv.cmp_min_pair_cost medn medab in 
        acc_cost + cost, acc_recost + recost
    in
    let cost, _ =     IntMap.fold cost_calculator b.meds (0,0) in 
    float_of_int cost


(** [f_codes s c] returns breakinv character subset
* of breakinv character set [s] whose codes are also in breakinv character set [c] *)
let f_codes s c = 
    let check x = All_sets.Integers.mem x c in
    let adder x y acc = 
        if check x then IntMap.add x y acc 
        else acc
    in
    let n_meds = IntMap.fold adder s.meds IntMap.empty
    and n_costs = IntMap.fold adder s.costs IntMap.empty in
    { s with meds = n_meds; costs = n_costs}

(** [f_codes s c] returns breakinv character subset
* of breakinv character set [s] whose codes are NOT in breakinv character set [c] *)
let f_codes_comp s c = 
    let check x = not (All_sets.Integers.mem x c) in
    let adder x y acc = 
        if check x then IntMap.add x y acc 
        else acc
    in
    let n_meds = IntMap.fold adder s.meds IntMap.empty
    and n_costs = IntMap.fold adder s.costs IntMap.empty in
    { s with meds = n_meds; costs = n_costs }

(** [compare_data a b] compares breakinv character
* set [a] and [b]. Returns 0 if they are the same, 
* otherwise (-1) or (1) *)
let compare_data a b =
    let comparator code medb acc =
        if acc = 0 then begin
            let meda = IntMap.find code a.meds in
            Breakinv.compare meda medb
        end else acc
    in
    IntMap.fold comparator b.meds 0


(** [to_formatter ref_codes attr t parent_t d] returns
* the map between breakinv character set [t] and its parents
* [parent_t] in the Tag.output format *)
let to_formatter ref_codes attr t (parent_t : t option) d : Xml.xml Sexpr.t list =
    let _, state = List.hd attr in 
    let output_breakinv code med acc =
        let med = 
            try
                List.find (fun med -> 
                               IntSet.mem med.BreakinvAli.ref_code  ref_codes
                          ) med.Breakinv.med_ls
            with Not_found -> failwith "Not found med -> to_formatter -> BreakinvCS"
        in         
        let cost, recost =
            match parent_t with  
            | None ->
                    0, 0
            | Some parent -> begin
                  let parent_med_ls = IntMap.find code parent.meds in
                  let breakinv_pam = parent_med_ls.Breakinv.breakinv_pam in
                  let gen_cost_mat_original =
                      parent_med_ls.Breakinv.gen_cost_mat_original in
                  let pure_gen_cost_mat_original =
                      parent_med_ls.Breakinv.pure_gen_cost_mat_original in    
                  let pure_gen_cost_mat_original =
                      match breakinv_pam.Data.re_meth with
                      | Some re_meth ->
                      (match re_meth with
                         | `Locus_Breakpoint c -> pure_gen_cost_mat_original
                         | `Locus_Inversion invc -> 
                            Breakinv.transform_matrix pure_gen_cost_mat_original invc
                      )
                      | None -> pure_gen_cost_mat_original
                  in
                  let alpha = parent_med_ls.Breakinv.alpha in
                  let parent_med = List.find 
                      (fun med -> 
                           IntSet.mem med.BreakinvAli.ref_code ref_codes 
                      ) parent_med_ls.Breakinv.med_ls
                  in            
                  let cost, recost =
                      match state with
                      | `String "Preliminary" ->
                      let cost,(recost1,recost2) = 
                      BreakinvAli.cmp_cost parent_med med 
                      gen_cost_mat_original pure_gen_cost_mat_original alpha breakinv_pam in
                      cost, recost1+recost2
                      (*  BreakinvAli.get_costs parent_med med.BreakinvAli.ref_code*)  
                      | `String "Final" ->
                      let cost,(recost1,recost2) = 
                      BreakinvAli.cmp_cost med parent_med 
                      gen_cost_mat_original pure_gen_cost_mat_original alpha breakinv_pam in
                      cost,recost1+recost2
                      (*BreakinvAli.get_costs med parent_med.BreakinvAli.ref_code*)
                      | `String "Single" ->
                      let cost,(recost1,recost2) = 
                      BreakinvAli.cmp_cost med parent_med 
                      gen_cost_mat_original pure_gen_cost_mat_original alpha breakinv_pam in
                      cost,recost1+recost2
(*
                        let cost = IntMap.find code t.costs in 
                            let recost = IntMap.find code t.recosts in 
                            (int_of_float cost), (int_of_float recost)    *)             
                      | _ ->
                            let cost = IntMap.find code t.costs in 
                            let recost = IntMap.find code t.recosts in 
                            (int_of_float cost), (int_of_float recost)
                  in 
                  cost, recost
              end 
        in 
        let seq = Sequence.to_formater med.BreakinvAli.seq t.alph in
        let module T = Xml.Characters in
        (PXML 
            -[T.breakinv]
                (* Attributes *)
                ([T.name] = [`String (Data.code_character code d)])
                ([T.cost] = [`Int cost])
                ([T.recost] = [`Int recost])
                ([T.definite] = [`Bool (cost > 0)])
                ([attr])

                { `String seq }
            --) :: acc
    in
    IntMap.fold output_breakinv t.meds []



(** [to_single ref_codes root single_parent mine] returns
* the single states of breakinv character set [mine] *) 
let to_single ref_codes (root : t option) single_parent mine = 
    let previous_total_cost = mine.total_cost in 
    let median code med (acc_meds, acc_costs, acc_recosts, acc_total_cost,
    acc_total_recost) =
        let amed = List.hd med.Breakinv.med_ls in
        let parent_med = IntMap.find code single_parent.meds in 
        let aparent_med = List.hd parent_med.Breakinv.med_ls in
        let cost, (recost1, recost2) = 
            match root with 
            | Some root -> 0,(0, 0) 
            | None ->
                 let bkpam = med.Breakinv.breakinv_pam in
                 let pure_gen_cost_mat_original =
                      match bkpam.Data.re_meth with
                      | Some re_meth ->
                      (match re_meth with
                         | `Locus_Breakpoint _ ->
                                 med.Breakinv.pure_gen_cost_mat_original
                         | `Locus_Inversion invc -> 
                            Breakinv.transform_matrix
                            med.Breakinv.pure_gen_cost_mat_original invc
                      )
                      | None -> med.Breakinv.pure_gen_cost_mat_original
                  in
                  BreakinvAli.cmp_cost amed aparent_med 
                      med.Breakinv.gen_cost_mat_original pure_gen_cost_mat_original
                      med.Breakinv.alpha bkpam
                      
        in 
        let single_med = {med with Breakinv.med_ls = [amed]} in 
        let new_single = IntMap.add code single_med acc_meds in
        let new_costs = IntMap.add code (float_of_int cost) acc_costs in 
        let new_recosts = IntMap.add code (float_of_int (recost1 + recost2) ) acc_recosts in
        let new_total_recost = acc_total_recost + recost1 + recost2 in
        new_single, new_costs, new_recosts, (acc_total_cost + cost),
        new_total_recost
    in  
    let meds, costs,  recosts, total_cost, total_recost = 
        match root with
        | Some root ->
              IntMap.fold median root.meds (IntMap.empty, IntMap.empty,
              IntMap.empty, 0, 0)
        | None ->
              IntMap.fold median mine.meds (IntMap.empty, IntMap.empty,
              IntMap.empty, 0, 0)
    in
    let subtree_recost = 
        match root with 
        | Some root -> root.subtree_recost
        | None -> mine.subtree_recost 
    in
    previous_total_cost, float_of_int total_cost, 
    {mine with meds = meds; 
         costs = costs;
         recosts = recosts;
         total_cost = float_of_int total_cost;
         subtree_recost = subtree_recost}
   



(** [get_active_ref_code t] returns active reference codes
* of breakinv character set [t]. One active reference code
* for one breakinv character *)
let get_active_ref_code t = 
    IntMap.fold 
        (fun _ meds (acc_ref_code, acc_child_ref_code) ->
             let ref_code, child1_ref_code, child2_ref_code = 
                 Breakinv.get_active_ref_code meds 
             in
             IntSet.add ref_code acc_ref_code,
             IntSet.add child2_ref_code (IntSet.add child1_ref_code acc_child_ref_code)  
        ) t.meds (IntSet.empty, IntSet.empty)
