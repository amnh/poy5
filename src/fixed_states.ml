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

(** A Sequence Character Set implementation *)
exception Illegal_Arguments
let () = SadmanOutput.register "Fixed_states" "$Revision: 2831 $"

type elt = SankCS.elt (*cside*)

type t = SankCS.t (*cside*)

type t_w_seqtbl = {
    states_and_distbl : t;
    sequence_table: Sequence.s array; (*we carry original sequence before
                                        transform(fixed_state) with us*)
    alph : Alphabet.a;
}

let cardinal t = 1 (*we have only 1 elt for each t*)

let f_codes x codes = (*this function is for filtering out elt from elt_array with some code, but
                        we only have one elt here, nothing should be done anyway*)
    x

let f_codes_comp x codes = x (*similar to f_codes*)

(*t = 1,2,3.  1:get states,2:get best states of leftchild,3:get best states of rightchild*)
let get_states x t = 
    SankCS.get_states x.states_and_distbl t

let get_earray x = 
    SankCS.get_earray x.states_and_distbl

let get_min_states x =
    let debug = false in
    let states = get_states x 1 in
    let earray = get_earray x in
    if debug then begin 
        Printf.printf "get min states, state arr = %!";
        Utl.printIntArr states;
        Printf.printf " e array = %!";
        Utl.printIntArr earray;
    end;
    let idx = ref (-1) in
    let _ , bestidx =
    Array.fold_left
        (fun (min_extra_cost,pos) extra_cost ->
            idx := !idx + 1;
            if extra_cost<min_extra_cost
                then extra_cost,!idx
                else min_extra_cost,pos)
        (Utl.large_int,0)
        earray
    in
    states.(bestidx),bestidx

let of_array fs taxon code = 
    let debug = false in
    if debug then
        Printf.printf "Fixed_states.off_array,taxon=%d,code:%d\n%!" taxon code;
    let distbl = Utl.float_to_int_mat fs.Data.costs in
    let seqtbl = fs.Data.seqs in
    let len = Array.length fs.Data.seqs in
    let states =
        let seqs = try Hashtbl.find_all fs.Data.codes taxon
                   with | Not_found -> []
        in
        let is_empty = seqs = [] in
        (*infinity here is not infinity on the c side, we pass (-1) instead*)
        Array.init len
            (fun x ->
                if is_empty || List.exists (fun y -> y = x) seqs
                    then Int32.of_int 0
                    else Int32.of_int (-1))
    in
    let ecode_arr = [|Int32.of_int code|] in
    let new_eltarr =
        SankCS.create_eltarr taxon code len ecode_arr [|states|] distbl true
    in
    {
        states_and_distbl = new_eltarr;
        sequence_table = seqtbl; 
        alph = fs.Data.original_dynspec.Data.alph;
    }


let to_formatter report_type attr tws d =
    match report_type with 
    | `Normal ->
    SankCS.to_formatter_with_seq true tws.sequence_table tws.alph attr tws.states_and_distbl None d
    | `StateOnly ->
    SankCS.to_formatter_with_seq false tws.sequence_table tws.alph attr tws.states_and_distbl None d
 
(*do nothing*)
let to_single parent child =  child, 0

let median code a b =
    let debug = false in
    if debug then Printf.printf "Fixed_states.median \n%!";
    let med,cost = SankCS.median code a.states_and_distbl
    b.states_and_distbl  in
    {states_and_distbl = med; sequence_table = a.sequence_table; alph = a.alph },
    cost

    
let median_3 p n c1 c2 = 
    let debug = false in
    if debug then Printf.printf "Fixed_states.median_3\n%!";
    let med3_eltarr = SankCS.median_3 p.states_and_distbl n.states_and_distbl
    c1.states_and_distbl c2.states_and_distbl in
    {n with states_and_distbl = med3_eltarr}


let distance a b =
    let debug = false in
    if debug then Printf.printf "Fixed_states.distance\n%!";
    SankCS.distance a.states_and_distbl b.states_and_distbl  

let dist_2 n a b =
    let cost = SankCS.dist_2 n.states_and_distbl a.states_and_distbl b.states_and_distbl in
    cost


let to_string x =
    let b_states = get_states x 1 in
    let b_leftchild = get_states x 2 in
    let b_rightchild = get_states x 3 in
    let res = 
        Array_ops.map_3 (fun a b c ->
            Printf.sprintf "Cost: %d - left: %d - right: %d" a b
            c) b_states b_leftchild b_rightchild 
    in
    String.concat "--" (Array.to_list res)





let compare_data a b = SankCS.compare_data a.states_and_distbl b.states_and_distbl


