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

let () = SadmanOutput.register "UtlGrappa" "$Revision: 2784 $"

(** utlGrappa module provides functions for rearrangement 
* operations such as computing inversion, breakpoint distances 
* between two gene orders arrays *)

let fprintf = Printf.fprintf 

(** [standardize genomeX genomeY] standardizes  
 * gene order [genomeX=(x1, x2,..., xk)] and  one of its 
 * permutations [genomeY=(y1, y2,...,yk)] 
 * into [sta_X] and [sta_Y] such that [genomeX=(1,...,k)]) 
 * For example: [genomeX] = (5, -3, 2) and [genomeY] = (-2, 3, 5)
 * [sta_X] = (1, 2, 3) and [sta_Y] = (-3, -2, 1 *)
let standardize genomeX genomeY =
(*  debug msg
    let print_intarr arr = 
        Printf.printf "[%!";
        Array.iter (Printf.printf "%d,%!") arr;
        Printf.printf "],%!";
    in
    Printf.printf "standardize,input :%! ";
    print_intarr genomeX; print_intarr genomeY; *)
    let max_index = Array.fold_left 
        (fun max_gene gene -> max max_gene (abs gene) ) (-1) genomeX 
    in 
    let num_gene = Array.length genomeX in  
    let sta_genomeX = Array.init num_gene (fun index -> index + 1) in
    let index_arr = Array.make (max_index + 1) 0 in     
    for idx = 0 to num_gene - 1 do
        match genomeX.(idx) > 0 with
        | true -> index_arr.( genomeX.(idx) ) <- (idx + 1)
        | false -> index_arr.( -genomeX.(idx) ) <- -(idx + 1)
    done; 
    let sta_genomeY = Array.init num_gene 
        (fun idx ->
            if ( genomeY.(idx) > max_index )||( genomeY.(idx) < -max_index )
            then
                Array.iteri (fun a arr ->
                    Printf.printf "[%d]=%d,%!" a arr
            ) genomeY
                ;
             match genomeY.(idx) > 0 with 
             | true -> index_arr.(genomeY.(idx)) 
             | false -> -index_arr.(-genomeY.(idx)) 
        )
    in  
  (*  debug msg 
    Printf.printf "end of standardize, output : %!";
    print_intarr sta_genomeX; print_intarr sta_genomeY;
    Printf.printf "\n%!"; *)
    sta_genomeX, sta_genomeY

let standardize3 (genomeX:int array) (genomeY:int array) (genomeZ:int array) = 
    let max_index = Array.fold_left 
        (fun max_gene gene -> max max_gene (abs gene) ) (-1) genomeX 
    in 
    let num_gene = Array.length genomeX in
    let sta_genomeX = Array.init num_gene (fun index -> index + 1) in
    let index_arr = Array.make (max_index + 1) 0 in     
    for idx = 0 to num_gene - 1 do
        match genomeX.(idx) > 0 with
        | true -> index_arr.( genomeX.(idx) ) <- (idx + 1)
        | false -> index_arr.( -genomeX.(idx) ) <- -(idx + 1)
    done; 
    let sta_genomeY = Array.init num_gene 
        (fun idx -> 
             match genomeY.(idx) > 0 with 
             | true -> index_arr.(genomeY.(idx)) 
             | false -> -index_arr.(-genomeY.(idx)) 
        )
    in  
    let sta_genomeZ = Array.init num_gene 
        (fun idx -> 
             match genomeZ.(idx) > 0 with 
             | true -> index_arr.(genomeZ.(idx)) 
             | false -> -index_arr.(-genomeZ.(idx)) 
        )
    in  
    sta_genomeX, sta_genomeY, sta_genomeZ

let de_standardize3 ori_arr standar_arr arrsize =
    assert(arrsize>0);
    let res_arr = Array.init arrsize 
    ( fun idx ->  
        let tmp = standar_arr.(idx) in
        if (tmp>0) then ori_arr.(tmp-1) 
        else - ori_arr.((abs tmp)-1)
    ) in
    res_arr

let to_ori_arr arr =
    Array.init (Array.length arr) (fun index->
        assert((arr.(index)<>0));
        if( arr.(index) mod 2 == 0) then -(arr.(index)/2)
        else (arr.(index)+1)/2 )

let from_ori_arr arr = 
    Array.init (Array.length arr) (fun index->
        assert((arr.(index)<>0));
        if (arr.(index)<0) then (abs arr.(index))*2
        else (arr.(index))*2-1
        )

let cmp_inversion_dis_multichrom (genomeX :int array) (genomeY : int array)
(delimiterX : int array) (delimiterY : int array) num_gen = 
    let set_seq = 1 and set_delimiters = 0 in
    let genomeX,genomeY = to_ori_arr genomeX, to_ori_arr genomeY
    in
    let sta_genomeX, sta_genomeY = standardize genomeX genomeY in
    let num_gen = Array.length genomeX in  
    let genome_arr = Grappa.c_create_empty_genome_arr 2 num_gen in  
    for index = 0 to num_gen - 1 do
        Grappa.c_set set_seq genome_arr 0 index sta_genomeX.(index);   
        Grappa.c_set set_seq genome_arr 1 index sta_genomeY.(index);   
    done;
    (* debug msg 
    let print_intarr arr = 
        Printf.printf "[%!";
        Array.iter (Printf.printf "%d,%!") arr;
        Printf.printf "],%!";
    in
     debug msg *)
    (* debug msg 
     Printf.printf "standardize input seq: %!";
    print_intarr genomeX; print_intarr genomeY; 
    Printf.printf "standardized seq: %!";
    print_intarr sta_genomeX; print_intarr sta_genomeY; 
    Printf.printf "deli array =  %!";
    print_intarr delimiterX; print_intarr delimiterY;
    Printf.printf "\n%!";
     debug msg *)
    let num_deliX = Array.length delimiterX and 
    num_deliY = Array.length delimiterY in
    (* if number of deli is 1, it is a single chromosome, no need to pass delimiters array *)
    for index = 0 to num_deliX - 1 do
        Grappa.c_set set_delimiters genome_arr 0 index delimiterX.(index);   
    done;
    for index = 0 to num_deliY - 1 do
        Grappa.c_set set_delimiters genome_arr 1 index delimiterY.(index);   
    done;
    let g0 = Grappa.c_get_one_genome genome_arr 0 in
    let g1 = Grappa.c_get_one_genome genome_arr 1 in  
    let inv_dis = Grappa.c_cmp_inv_dis g0 g1 num_gen 0 in   
    inv_dis 


(** [cmp_inversion_dis genomeX genomeY circular] computes
 * the inversion distance between two given gene orders 
 * [genomeX=(x1, x2, ... xk)] and [genomeY=(y1, y2,..yk)]. 
 * Note [genomeY] is a permutation of [genomeX].
 * Compute the inversion distance between [genomeX] and [genomeY] using GRAPPA functions
 * For example: [genomeX] = (-6, 1, 5), [genomeY] = (-5, 1, 6) *)
let cmp_inversion_dis (genomeX : int array) (genomeY : int array) circular  =
    let set_seq = 1 in 
    let sta_genomeX, sta_genomeY = standardize genomeX genomeY in
    let num_gen = Array.length genomeX in  
    let genome_arr = Grappa.c_create_empty_genome_arr 2 num_gen in  
    for index = 0 to num_gen - 1 do
        Grappa.c_set set_seq genome_arr 0 index sta_genomeX.(index);   
        Grappa.c_set set_seq genome_arr 1 index sta_genomeY.(index);   
    done;
    (* distance calc for multichromosome should be different than
    * singlechromosome, back to this later
    * for index = 0 to num_gen - 1 do
        Grappa.c_set set_delimiters genome_arr 0 index sta_genomeX.(index);   
    done;
    for index = 0 to num_gen - 1 do
        Grappa.c_set set_delimiters genome_arr 1 index sta_genomeY.(index);   
    done;
    *)
    let g0 = Grappa.c_get_one_genome genome_arr 0 in
    let g1 = Grappa.c_get_one_genome genome_arr 1 in  
    let inv_dis = Grappa.c_cmp_inv_dis g0 g1 num_gen circular in   
    inv_dis 


let find_better_capping (genomeX : int array) (genomeY : int array) (delimiterX:
    int array) (delimiterY :int array) =
    (* debug message 
    let print_intarr arr = 
        Printf.printf "[%!";
        Array.iter (Printf.printf "%d,%!") arr;
        Printf.printf "]\n%!";
    in
     debug message *)
    (*debug message 
    Printf.printf "find better capping ,input seqcodes/delimiters:\n %!";
    print_intarr genomeX; print_intarr genomeY; 
    print_intarr delimiterX; print_intarr delimiterY;
     debug message *)
    let set_seq = 1 and set_delimiters = 0 in
    let genomeX,genomeY = to_ori_arr genomeX, to_ori_arr genomeY
    in
    let ori_genomeX = genomeX in
    let sta_genomeX, sta_genomeY = standardize genomeX genomeY in
    let num_gen = Array.length genomeX in  
    let genome_arr = Grappa.c_create_empty_genome_arr 2 num_gen in  
    for index = 0 to num_gen - 1 do
        Grappa.c_set set_seq genome_arr 0 index sta_genomeX.(index);   
        Grappa.c_set set_seq genome_arr 1 index sta_genomeY.(index);   
    done;
    let num_deliX = Array.length delimiterX 
    and num_deliY = Array.length delimiterY in
    for index = 0 to num_deliX - 1 do
        Grappa.c_set set_delimiters genome_arr 0 index delimiterX.(index);   
    done;
    for index = 0 to num_deliY - 1 do
        Grappa.c_set set_delimiters genome_arr 1 index delimiterY.(index);   
    done;
    let g0 = Grappa.c_get_one_genome genome_arr 0 in
    let g1 = Grappa.c_get_one_genome genome_arr 1 in  
    let new_g1 = Grappa.get_better_capping_genome g0 g1 num_gen in 
    let new_g1_intarr = Grappa.genome_to_gene_intarr new_g1 num_gen in
    let delinum = Grappa.c_get_delimiter_num new_g1 in
    let deli_arr = Grappa.get_delimiter_arr new_g1 delinum in
    let resarr = de_standardize3 ori_genomeX new_g1_intarr num_gen in
    let resarr = from_ori_arr resarr in
    (*debug msg
    Printf.printf "res seq/deli = \n%!"; 
    print_intarr resarr; print_intarr deli_arr;
    debug msg*)
    resarr,deli_arr

let inv_med (medsov : Data.median_solver_t) (genomeX : int array) (genomeY : int
array) (genomeZ : int array) (delimiters_lstlst : int list list) circular =
    let set_seq = 1 and set_delimiters = 0 in
    (* debug message 
    let print_intarr arr = 
        Printf.printf "[%!";
        Array.iter (Printf.printf "%d,%!") arr;
        Printf.printf "],%!";
    in
    Printf.printf "inv_med ,input seqcodes: %!";
    print_intarr genomeX; print_intarr genomeY; print_intarr genomeZ;
    debug message*) 
    let ori_genomeX = genomeX in
    let genomeX, genomeY, genomeZ = standardize3 genomeX genomeY genomeZ in
    let num_gen = Array.length genomeX in 
    let genome_arr = Grappa.c_create_empty_genome_arr 3 num_gen in  
     assert (num_gen>0); assert(num_gen<=2048);
    (* for alert-median3 solver to work , sequence cannot be empty, also there
    * is a MAX_STR_LEN=2048 macro in grappa, if we need to work with longer sequence, 
    * modify the macro in structs.h to accomodate new requirement*)
    for index = 0 to num_gen - 1 do
        Grappa.c_set set_seq genome_arr 0 index genomeX.(index);   
        Grappa.c_set set_seq genome_arr 1 index genomeY.(index);   
        Grappa.c_set set_seq genome_arr 2 index genomeZ.(index);
    done;
    assert( (List.length delimiters_lstlst) = 3 );
    let deli1 = List.hd delimiters_lstlst 
    and deli2 = List.nth delimiters_lstlst 1
    and deli3 = List.nth delimiters_lstlst 2 in
    let deliX = Array.of_list deli1 
    and deliY = Array.of_list deli2
    and deliZ = Array.of_list deli3 in
    let num_deliX = Array.length deliX
    and num_deliY = Array.length deliY
    and num_deliZ = Array.length deliZ in
    (* debug msg  
   (*  Printf.printf "standardize input seq: %!";
    print_intarr genomeX; print_intarr genomeY; print_intarr genomeZ; *)
     Printf.printf "circular = %d, num_deli = %d,%d,%d,deli array = \n%!" 
     circular num_deliX num_deliY num_deliZ;
     if num_deliX >1 then  print_intarr deliX; 
     if num_deliY >1 then  print_intarr deliY; 
     if num_deliZ >1 then  print_intarr deliZ;
     debug msg *)
    (* if number of deli is 1, it is a single chromosome, no need to pass
    * delimiters array *)
    for index = 0 to num_deliX - 1 do
        Grappa.c_set set_delimiters genome_arr 0 index deliX.(index);   
    done;
    for index = 0 to num_deliY - 1 do
        Grappa.c_set set_delimiters genome_arr 1 index deliY.(index);   
    done;
    for index = 0 to num_deliZ - 1 do
        Grappa.c_set set_delimiters genome_arr 2 index deliZ.(index);
    done;
       
    let g0 = Grappa.c_get_one_genome genome_arr 0 in
    let g1 = Grappa.c_get_one_genome genome_arr 1 in 
    let g2 = Grappa.c_get_one_genome genome_arr 2 in
    let gmed =
        match medsov with
        |`Vinh ->
                failwith "Vinh median solver is not in grappa"
        |`Albert ->
                Grappa.get_med3_genome 1 g0 g1 g2 num_gen circular
        |`Siepel ->
              Grappa.get_med3_genome 2 g0 g1 g2 num_gen circular
        |`BBTSP ->
            Grappa.get_med3_genome 3 g0 g1 g2 num_gen circular
        |`COALESTSP ->
             Grappa.get_med3_genome 4 g0 g1 g2 num_gen circular
        |`SimpleLK ->
            Grappa.get_med3_genome 5 g0 g1 g2 num_gen circular
        |`ChainedLK ->
            Grappa.get_med3_genome 6 g0 g1 g2 num_gen circular
        |`MGR ->
            Grappa.get_med3_genome 7 g0 g1 g2 num_gen circular
    in
    let med3arr = Grappa.genome_to_gene_intarr gmed num_gen in
    let delinum = 
        match medsov with
        | `MGR ->
                Grappa.c_get_delimiter_num gmed
        | _ -> 0
    in
    let deli_arr = 
        match medsov with
        | `MGR ->
                Grappa.get_delimiter_arr gmed delinum
        | _ -> [||]
    in
    (* debug msg 
     Printf.printf "output ori seqcode = %!"; print_intarr med3arr;
     print_newline();
     Printf.printf "delinum=%d,delimiters = %!" delinum; 
     print_intarr deli_arr;
     print_newline();
      debug msg*)
    let resarr = de_standardize3 ori_genomeX med3arr num_gen in
    (* debug msg  
     Printf.printf "output seqcode = %!"; print_intarr resarr;
     print_newline();
    debug msg*)
    resarr,deli_arr


(** [cmp_breakpoint_dis genomeX genomeY circular] computes
 * the breakpoint distance between two given gene orders 
 * [genomeX=(x1, x2, ... xk)] and [genomeY=(y1, y2,..yk)]. 
 * Note that orientations are ignored. *)
let cmp_breakpoint_dis (genomeX : int array) (genomeY : int array) circular = 

    let _, sta_genomeY = standardize genomeX genomeY in 
                        
    let sta_genomeY = Array.map (fun gene -> abs gene) sta_genomeY in 

    let num_gene = Array.length sta_genomeY in     
    let dis = ref 0 in 
    for pos = 0 to num_gene - 2 do
        if abs(sta_genomeY.(pos) - sta_genomeY.(pos + 1) ) != 1 then 
            dis := !dis + 1
    done;
    

    match circular with
    | 0 -> !dis 
    | _ -> 
          let cir_dis = abs ( sta_genomeY.(0) - sta_genomeY.(num_gene - 1) ) in 
          if (cir_dis = 1) || (cir_dis = num_gene - 1) then !dis
          else !dis + 1


(** [cmp_oriented_breakpoint_dis genomeX genomeY circular] computes
 * the breakpoint distance between two given gene orders 
 * [genomeX=(x1, x2, ... xk)] and [genomeY=(y1, y2,..yk)]. 
 * Note that orientations are taken into account. 
 * For example: [genomeX] = (-6, 1, 5), [genomeY] = (-5, 1, 6) *)
let cmp_oriented_breakpoint_dis (genomeX : int array) (genomeY : int array)
        circular = 
    let _, sta_genomeY = standardize genomeX genomeY in
                        
    let num_gene = Array.length sta_genomeY in 
    let dis = ref 0 in 
    for pos = 0 to num_gene - 2 do
        if sta_genomeY.(pos) + 1 != sta_genomeY.(pos + 1) then
            dis := !dis + 1
    done;
    

    match circular with
    | 0 -> !dis 
    | _ -> 
          let fi_g = sta_genomeY.(0) 
          and la_g = sta_genomeY.(num_gene -1 ) in 
          if ( (la_g = num_gene) && (fi_g = 1) ) || 
              ( (la_g = -1) && (fi_g = -num_gene) ) || 
              (la_g + 1 = fi_g) then !dis
          else !dis + 1

(** [get_ordered_permutation genomeX] returns 
 * the ordered permutation [genomeY]=(y1,..,yk | y1 < ... < yk) 
 * of [genomeX] =(x1, x2, ... xk), 
 * For example. [genomeX] = (-6, 1, 5), [genomeY] = (1, 5, 6) *)
let get_ordered_permutation genomeX = 
    let max_index = Array.fold_left 
        (fun max_gene gene -> max max_gene (abs gene) ) (-1) genomeX 
    in 

    let num_gene = Array.length genomeX in  
    let index_arr = Array.make (max_index + 1) (-1) in 
    for idx = 0 to num_gene - 1 do
        index_arr.( abs genomeX.(idx) ) <- idx;
    done; 

    let genomeY_ls = 
        Array.fold_right (fun idx genomeY ->
                          match idx with
                          | -1 -> genomeY
                          | _ -> (abs genomeX.(idx))::genomeY) index_arr []
    in 
    Array.of_list genomeY_ls
    

(** [cmp_self_inversion_dis genome circular] 
 * computes the inversion distance between a gene orders [genome=(x1, x2, ... xk)]
 * and [Y=(y1,..,yk | y1 < ... < yk) using GRAPPA functions
 * where [Y] is an ordered permutation of |genome|
 * For example. [genome] = (-6, 1, 5), Y = (1, 5, 6) *)
let cmp_self_inversion_dis (genome : int array)  circular  =
    let ordered_permutation = get_ordered_permutation genome in     
    let dis = cmp_inversion_dis genome ordered_permutation circular in 
    dis


(** [cmp_self_breakpoint_dis genome circular] 
 * computes the breakpoint distance between a gene orders [genome=(x1, x2, ... xk)]
 * and [Y=(y1,..,yk | y1 < ... < yk) using GRAPPA functions
 * where [Y] is an ordered permutation of |genome|
 * For example. [genome] = (-6, 1, 5), Y = (1, 5, 6) *)
let cmp_self_breakpoint_dis (genome : int array)  circular  =
    let ordered_permutation = get_ordered_permutation genome in 
    cmp_breakpoint_dis genome ordered_permutation circular



(** [cmp_self_oriented_breakpoint_dis genome circular] 
 * computes the breakpoint distance between a gene orders [genome=(x1, x2, ... xk)]
 * and [Y=(y1,..,yk | y1 < ... < yk) using GRAPPA functions
 * where [Y] is an ordered permutation of |genome|
 * For example. [genome] = (-6, 1, 5), Y = (1, 5, 6) *)
let cmp_self_oriented_breakpoint_dis (genome : int array)  circular  =
    let ordered_permutation = get_ordered_permutation genome in 
    cmp_oriented_breakpoint_dis genome ordered_permutation circular

