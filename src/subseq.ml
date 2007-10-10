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

let () = SadmanOutput.register "Subseq" "$Revision: 1644 $"
(** The data and functions applied to a subsequence in a chromosome *)

type type_t = Alied | Deleted | Both


let fprintf = Printf.fprintf

type subseq_t = {
    mutable id : int;
    mutable sta : int;
    mutable en : int;
    (* One subsequence may belong to different blocks due to the duplications *)
    mutable block_id_ls : int list; 
}

let print ?(channel = stdout) (subseq : subseq_t) = 
    fprintf channel "id: %i --> (%i, %i -> %i), block_ls_id: " subseq.id 
        subseq.sta subseq.en (subseq.en - subseq.sta + 1);
    List.iter (fprintf channel "%i  ") subseq.block_id_ls;
    fprintf channel "\n";
    flush channel
    
let get_sum (subseq_ls : subseq_t list) = 
    List.fold_left (fun sum subseq -> 
                        match subseq.block_id_ls with
                        | [] -> sum
                        | _ -> sum + (subseq.en - subseq.sta + 1)
                    ) 0 subseq_ls
                    

let delete_block (sq : subseq_t) (block_id : int) = 
    sq.block_id_ls <- List.filter (fun id -> id != block_id) sq.block_id_ls
    
let is_equal (s1 : subseq_t) (s2 : subseq_t) = 
    (s1.sta = s2.sta) && (s1.en = s2.en)
    

let get_intersection s1 e1 s2 e2 = 
    (max s1 s2), (min e1 e2)


(** Compute the cost of a subseq when aligning with all gaps *)
let cmp_del_cost del_subseq seq cost_mat = 
    let cost_model = Cost_matrix.Two_D.affine cost_mat in 
    let gap_opening_cost = 
        match cost_model with
        | Cost_matrix.Affine (gap_opening_cost) -> gap_opening_cost
        | _ -> 0 
    in
    let gap_code = Cost_matrix.Two_D.gap cost_mat in
    let rec travel pos cost = 
        if pos > del_subseq.en then cost
        else begin
            let code = Sequence.get seq pos in 
            travel (pos + 1) (cost + Cost_matrix.Two_D.cost gap_code code cost_mat)                      
        end
    in
    travel del_subseq.sta gap_opening_cost


let is_free sq = (List.length sq.block_id_ls) = 0

let get_subseq seq ss  =
    Sequence.sub seq ss.sta (ss.en - ss.sta + 1)

let get_len ss = ss.en - ss.sta + 1
