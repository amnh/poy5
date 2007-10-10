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

let () = SadmanOutput.register "UtlGrappa" "$Revision: 2243 $"

(** This module provide functions for rearrangement operations such as computing
    inversion, breakpoint distances between two gene orders arrays *)
let fprintf = Printf.fprintf 

(** This is only for debuging purpose *)
let write_2genome_file (filename : string) 
        (genomeX : int array) (genomeY : int array) = 

    let outfile = open_out filename in 
    fprintf outfile ">S0\n";
    Array.iter (fun gene -> fprintf outfile "%i " gene) genomeX;
    fprintf outfile "\n>S1\n";
    Array.iter (fun gene -> fprintf outfile "%i " gene) genomeY;
    fprintf outfile "\n";
    close_out outfile
    

(** Given gene order [X=(x1, x2,..., xk)] and  one of its permutations [Y=(y1, y2,...,yk)] 
 * Map [X] and [Y] into [sta_X] and [sta_Y] such that [X=(1,...,k)]) 
 * For example: X = (5, -3, 2) and Y = (-2, 3, 5)
 * sta_X = (1, 2, 3) and sta_Y = (-3, -2, 1 *)
let standardize genomeX genomeY = 
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
(*
    Utl.printIntArr genomeX;
    Utl.printIntArr genomeY;
    Utl.printIntArr index_arr;
*)
    let sta_genomeY = Array.init num_gene 
        (fun idx -> 
             match genomeY.(idx) > 0 with 
             | true -> index_arr.(genomeY.(idx)) 
             | false -> -index_arr.(-genomeY.(idx)) 
        )
    in  
    sta_genomeX, sta_genomeY



(** Given two gene orders [X=(x1, x2, ... xk)] and [Y=(y1, y2,..yk)]. 
 * Note [Y] is a permutation of [X].
 * Compute the inversion distance between X and Y using GRAPPA functions
 * For example:
 * 1. X = (-6, 1, 5), Y = (-5, 1, 6) *)
let cmp_inversion_dis (genomeX : int array) (genomeY : int array) circular  =
    let sta_genomeX, sta_genomeY = standardize genomeX genomeY in
                        
    let num_gen = Array.length genomeX in  
    let genome_arr = Grappa.c_create_empty_genome_arr 2 num_gen in  
    for index = 0 to num_gen - 1 do
        Grappa.c_set genome_arr 0 index sta_genomeX.(index);   
        Grappa.c_set genome_arr 1 index sta_genomeY.(index);   
        
    done;

    let g0 = Grappa.c_get_one_genome genome_arr 0 in
    let g1 = Grappa.c_get_one_genome genome_arr 1 in  

    let inv_dis = Grappa.c_cmp_inv_dis g0 g1 num_gen circular in   
    inv_dis 
 


(** Given two gene orders [X=(x1, x2, ... xk)] and [Y=(y1, y2,..yk)]. 
 * Note [Y] is a permutation of [X].
 * Compute the breakpoint distance between X and Y (orientations are ignored).
 * For example: X = (-6, 1, 5), Y = (-5, 1, 6) *)
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


(** Given two gene orders [X=(x1, x2, ... xk)] and [Y=(y1, y2,..yk)]. 
 * Note [Y] is a permutation of [X].
 * Compute the breakpoint distance between X and Y (orientations are taken into account).
 * For example: X = (-6, 1, 5), Y = (-5, 1, 6) *)
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

(** Given a gene orders [X=(x1, x2, ... xk)], return 
 * the ordered permutation [Y=(y1,..,yk | y1 < ... < yk| of |X|
 * For example. X = (-6, 1, 5), Y = (1, 5, 6) *)
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
    

(** Given a gene orders [X=(x1, x2, ... xk)], compute 
 * Compute the inversion distance between [X] and 
 * [Y=(y1,..,yk | y1 < ... < yk| using GRAPPA functions
 * where [Y] is an ordered permutation of |X|
 * For example. X = (-6, 1, 5), Y = (1, 5, 6) *)
let cmp_self_inversion_dis (genome : int array)  circular  =
    let ordered_permutation = get_ordered_permutation genome in     
    let dis = cmp_inversion_dis genome ordered_permutation circular in 
    dis



(** Given a gene orders [X=(x1, x2, ... xk)], compute 
 * Compute the breakpoint distance between [X] and 
 * [Y=(y1,..,yk | y1 < ... < yk| using GRAPPA functions
 * where [Y] is an ordered permutation of |X|
 * For example. X = (-6, 1, 5), Y = (1, 5, 6) *)
let cmp_self_breakpoint_dis (genome : int array)  circular  =
    let ordered_permutation = get_ordered_permutation genome in 
    cmp_breakpoint_dis genome ordered_permutation circular



(** Given a gene orders [X=(x1, x2, ... xk)], compute 
 * Compute the breakpoint distance between [X] and 
 * [Y=(y1,..,yk | y1 < ... < yk| using GRAPPA functions
 * where [Y] is an ordered permutation of |X|
 * For example. X = (-6, 1, 5), Y = (1, 5, 6) *)
let cmp_self_oriented_breakpoint_dis (genome : int array)  circular  =
    let ordered_permutation = get_ordered_permutation genome in 
    cmp_oriented_breakpoint_dis genome ordered_permutation circular

