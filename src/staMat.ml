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

let () = SadmanOutput.register "StaMat" "$Revision: 1644 $"
(** Statically local alignment:
 * This module concerns about the statistical significance of local
 * alignments. The local alignments which are statistically significant are used
 * to constrcut the map between two
 * chromosomes 
 *)

(** the default parameter in equation 4 of Karlin and Altschul, 1990 *)
let k_val = ref 0.335141
let lambda_val = ref 1.097076
let sta_sig_prob = 5.0

let get_score = Cost_matrix.Two_D.cost
let get_float_score code1 code2  score_mat  = 
    float (Cost_matrix.Two_D.cost code1 code2 score_mat )

let find_lambda (score_mat : Cost_matrix.Two_D.m) (code_frq_arr : float array) =
    let lower_lambda = ref 0.0 in 
    let upper_lambda = ref 1.0 in 
    
    let size : int = Cost_matrix.Two_D.alphabet_size score_mat in 
    let compute_lambda (a_lambda : float) = 
        let sum = ref 0.0 in
        for code1 = 1 to size - 1 do
            for code2 = 1 to size - 1 do
                sum := !sum +. code_frq_arr.(code1) *. 
                    code_frq_arr.(code2) *. exp(a_lambda *. 
                    (get_float_score code1 code2 score_mat));
            done
        done;
        
        !sum
    in
    
    while compute_lambda (!upper_lambda) < 1.0 do
        upper_lambda := !upper_lambda *. 2.0
    done;
    
    let lambda_sensitivity = 0.0001 in
    while (!upper_lambda -. !lower_lambda > lambda_sensitivity) do
        let mid_lambda = (!lower_lambda +. !upper_lambda) /. 2.0 in 
        if compute_lambda (mid_lambda) > 1.0 then upper_lambda := mid_lambda
        else lower_lambda := mid_lambda
    done;

    (!lower_lambda +. !upper_lambda) /. 2.0 
(*******************************************************************)
    
(** Find the value of K in Equation 5 and the 
 * Appendix of Karlin and Altschul, 1990 *)
let find_K (lambda : float) (score_mat : Cost_matrix.Two_D.m) 
        (code_frq_arr : float array) = 

    let size = Cost_matrix.Two_D.alphabet_size score_mat in 

    let find_C () =
        let min_score  = ref (get_score 1 1 score_mat) in 
        let max_score  = ref (get_score 1 1 score_mat) in 
        for code1 = 1 to size - 1 do
            for code2 = 1 to size - 1 do
                let score = get_score code1 code2 score_mat in 
                min_score := min !min_score score;
                max_score := max !max_score score;
            done
        done;
        
        
        let converged_k = 20 in 
        let min_S_k = !min_score * converged_k in
        let max_S_k = !max_score * converged_k in 

        (* We have to shift because the starting point of array is 0*)
        let shift = abs(min_S_k) in 
        let pre_prob_S_k_arr = ref (Array.create (abs(min_S_k) + 
                                      abs(max_S_k) + 1) 0.0) in 
        !pre_prob_S_k_arr.(0 + shift) <- 1.0;
        
        
        let compute_numerator () = 
            let numerator = ref 0.0 in
            for k = 1 to converged_k do
                let prob_S_k_arr = Array.create (abs(min_S_k) + 
                                                     abs(max_S_k) + 1) 0.0 in 
                for value = min_S_k + abs(!min_score) to
                    max_S_k - abs(!max_score) do

                    for code1 = 1 to size - 1 do
                        for code2 = 1 to size -1 do
                            let score = get_score code1 code2 score_mat in 
                            let prob = code_frq_arr.(code1) *. 
                                code_frq_arr.(code2) in 
                            let s_k = value + score in 
                            prob_S_k_arr.(s_k  + shift) <- 
                                prob_S_k_arr.(s_k + shift) +. 
                                !pre_prob_S_k_arr.(value + shift) *. prob
                        done (* end of code2 *)
                    done (* end of code1 *)         
                done;
            
                let sum = ref 0.0 in
                for s_k = min_S_k to max_S_k do 
                    if (s_k < 0) then sum := !sum +. exp(lambda *. 
                                     (float s_k)) *. prob_S_k_arr.(s_k + shift)
                    else sum := !sum +. prob_S_k_arr.(s_k + shift)
                done;
            
                numerator := !numerator +. (!sum /. (float k));
                
                pre_prob_S_k_arr := prob_S_k_arr;
            done; (* End of for k *)
            
            exp( -2.0 *. !numerator)
        in (* End of compute numerator*)
        
        
        let compute_denominator () = 
            let denominator = ref 0.0 in 
            for code1 = 1 to size - 1 do
                for code2 = 1 to size - 1 do
                    let s_1 = get_float_score code1 code2 score_mat in
                    let prob = code_frq_arr.(code1) *. code_frq_arr.(code2) in 
                    denominator := !denominator +. prob *. 
                        ( s_1 *. exp (lambda *. s_1));
                done (* end of code2 *)         
            done; (* end of code1 *)        
            
            !denominator *. lambda
        in  (* End of compute denominator*) 
        
        let numerator = compute_numerator () in
        let denominator = compute_denominator () in 


        
        let c = numerator /. denominator in 
        c
    in (* End of find_C*)
    
    let find_delta = 
        let rec find_largest_common_divisor () = 
            (* Find the greatest common divisor between two integers a and b*)
            let rec find_gcd (a : int) (b : int) = 
                let r = a mod b in
                    if r = 0 then b
                    else find_gcd b r
            in 
            
            let gcd = ref (get_score 1 1 score_mat) in
            for code1 = 1 to size - 1 do
                for code2 = 1 to size - 1 do
                    gcd := find_gcd !gcd (get_score code1 code2 score_mat)                  
                done
            done;
            
            !gcd
        in (* End of find_largest_common_divisor_code*)
        
        let delta = find_largest_common_divisor () in 
        delta
    in  (* End of find_delta*)
    
    let delta = find_delta in       
(*  fprintf stdout "delta: %i \n" delta; flush stdout; *)
    
    let c = find_C () in
(*  fprintf stdout "C value: %f \n" c; flush stdout; *)
    let k =  c *. (lambda *. (float delta) /. (1.0 -. 
           exp (-.lambda *. (float delta))) ) 
    in     
    k
    
    
(*******************************************************************)
let init (score_mat : Cost_matrix.Two_D.m) (code_frq_arr : float array) = 
    lambda_val  := find_lambda score_mat code_frq_arr;
(*  fprintf stdout "lambda: %f \n" !lambda_val; flush stdout; *)
    k_val := find_K (!lambda_val) score_mat code_frq_arr
(*  fprintf stdout "k: %f \n" !k_val; flush stdout *)
(*******************************************************************)


let create_default_score () : Cost_matrix.Two_D.m = 
    let score_mat = Cost_matrix.Two_D.clone (Cost_matrix.Two_D.default) in 
    let size = Cost_matrix.Two_D.alphabet_size score_mat in 
    for code1 = 1 to size - 1 do
        for code2 = 1 to size - 1 do
            match code1 = code2 with
                | true -> Cost_matrix.Two_D.set_cost code1 code2 score_mat 1
                | false -> Cost_matrix.Two_D.set_cost code1 code2 
                      score_mat (-1)
        done
    done;
    
    score_mat
(*======================================================*)

let cmp_code_frq (size : int) (seq1 : Sequence.s) (seq2 : Sequence.s) = 
    let code_frq_arr = Array.make (size + 1) 0.0 in 
    let count seq = 
        for index = 0 to (Sequence.length seq) - 1 do
            let code = Sequence.get seq index in 
            code_frq_arr.(code) <- code_frq_arr.(code) +. 1.0
        done
    in
    count seq1;
    count seq2;
    let num_code = (Sequence.length seq1) + (Sequence.length seq2) in 
    Array.iteri (fun code frq -> 
                     code_frq_arr.(code) <- frq /. (float num_code) 
                ) code_frq_arr;
    
    code_frq_arr
    
(*==================================================*)  


(** Compute the parameter lambda_val and k_val based on two sequences seq1 and
    seq2 and default cost matrix in Cost_matrix.Two_D *)
let init_default (seq1 : Sequence.s) (seq2 : Sequence.s) = 
    let score_mat : Cost_matrix.Two_D.m = create_default_score () in 
    let size = Cost_matrix.Two_D.alphabet_size score_mat in 
    let code_frq_arr = cmp_code_frq size seq1 seq2 in 
    
    lambda_val := find_lambda score_mat code_frq_arr;
    k_val := find_K !lambda_val score_mat code_frq_arr
    

(** Compute the minimum score of a statically significant local alignment. This
    score depends on the lengths of two sequences, cost matrix and base frequency*)    
let cmp_min_score (seq1_len : int) (seq2_len : int) = 
    let x = -.log( (sta_sig_prob /. 100.0) /. !k_val) /. !lambda_val in 

    let min_score = ( log( float seq1_len) +. log( float seq2_len) ) /. 
        !lambda_val +. x in 
    int_of_float min_score
