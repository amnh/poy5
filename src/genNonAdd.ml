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


(*this module is for prealigned custom alphabet datatype.
* call it like this
* read(custom_alphabet:("data file","cost matrix file",prealigned))
* *)

let debug = true

let fprintf = Printf.fprintf

(*cost tuple, like the one in seqCS.ml*)
type cost_tuple = 
{
    min : float;
    max : float;
}

type gnonadd_sequence = {
    seq : Sequence.s;
    costs : cost_tuple;
}

let init_gnonadd_t in_seq = 
    { seq = in_seq; costs = { min = 0.0; max = 0.0 }; }

let get_max_cost ct = ct.max

let get_min_cost ct = ct.min

(*[make_cost tmpcost] make cost_tuple with tmpcost -- to do: what is maxcost?*)
let make_cost tmpcost = 
    let tmpcost = float_of_int tmpcost in
    {min = tmpcost; max= tmpcost}

(*[get_cost] and [get_median] are two functions from cost_matrx*)
let get_cost = Cost_matrix.Two_D.cost
let get_median = Cost_matrix.Two_D.median

let get_cost3d = Cost_matrix.Three_D.cost
let get_median3d = Cost_matrix.Three_D.median



(**[get_distance gnoadd1 gnoadd2 cost_mat] return distance between two sequence.
* also return a median sequence.
* remember that these sequence are prealigned. just add cost from each column.*)
let get_distance_3d gnoadd1 gnoadd2 gnoadd3 cost_mat = 
    let size = Cost_matrix.Three_D.alphabet_size cost_mat in
    let seq1,seq2,seq3 = gnoadd1.seq,gnoadd2.seq,gnoadd3.seq in
    if ( (Sequence.length seq1)<>(Sequence.length seq2) )||( (Sequence.length
    seq3)<>(Sequence.length seq2) ) then begin
        Status.user_message Status.Error
        (
            "The@ prealigned@ sequences@ do@ not@ have@ the@ same@ length." ^
            string_of_int(Sequence.length seq1) ^ " !=@ " ^ string_of_int(Sequence.length seq2)
            ^ " or " ^ string_of_int(Sequence.length seq3) ^ " !=@ " ^ string_of_int(Sequence.length seq2)
        );
        failwith "Illegal prealigned molecular sequences."
    end;
    let arr1,arr2,arr3 = Sequence.to_array seq1, Sequence.to_array seq2,
    Sequence.to_array seq3 in
    if debug then begin
       Printf.printf "genNonAdd.distance, alphabet size = %d,arr1/arr2/arr3=%!" size;
       Utl.printIntArr arr1;
       Utl.printIntArr arr2;
       Utl.printIntArr arr3;
    end;
    let rescost = ref 0 in
    let medarr = 
    Array_ops.map_3 (fun code1 code2 code3 ->
        let c = get_cost3d code1 code2 code3 cost_mat in
        rescost := !rescost + c;
        get_median3d code1 code2 code3 cost_mat
    ) arr1 arr2 arr3
    in
    if debug then begin
        Printf.printf "distance = %d, med arr = %!" !rescost;
        Utl.printIntArr medarr;
    end;
    !rescost, Sequence.of_array medarr



(**[get_distance gnoadd1 gnoadd2 cost_mat] return distance between two sequence.
* also return a median sequence.
* remember that these sequence are prealigned. just add cost from each column.*)
let get_distance gnoadd1 gnoadd2 cost_mat = 
    let size = Cost_matrix.Two_D.alphabet_size cost_mat in
    let seq1,seq2 = gnoadd1.seq,gnoadd2.seq in
    if (Sequence.length seq1)<>(Sequence.length seq2) then begin
        Status.user_message Status.Error
        (
            "The@ prealigned@ sequences@ do@ not@ have@ the@ same@ length." ^
            string_of_int(Sequence.length seq1) ^ " !=@ " ^ string_of_int(Sequence.length seq2)
        );
        failwith "Illegal prealigned molecular sequences."
    end;
    let arr1,arr2 = Sequence.to_array seq1, Sequence.to_array seq2 in
    if debug then begin
       Printf.printf "genNonAdd.distance, alphabet size = %d,arr1/arr2=%!" size;
       Utl.printIntArr arr1;
       Utl.printIntArr arr2;
    end;
    let rescost = ref 0 in
    let medarr = 
    Array_ops.map_2 (fun code1 code2 ->
        let c = get_cost code1 code2 cost_mat in
        rescost := !rescost + c;
        get_median code1 code2 cost_mat
    ) arr1 arr2
    in
    if debug then begin
        Printf.printf "distance = %d, med arr = %!" !rescost;
        Utl.printIntArr medarr;
    end;
    !rescost, Sequence.of_array medarr

(** [distance gnoadd1 gnoadd2 cost_mat] just return the distance between two
* general nonaddictive sequence. ignore the median seq.*)
let distance gnoadd1 gnoadd2 cost_mat =
    if debug then Printf.printf "genNonAdd.distance, call get_distance\n%!";
    let dis,_ = get_distance gnoadd1 gnoadd2 cost_mat in
    dis

(** [median cost_mat a b] return median of two general nonaddictive sequence*)
let median cost_mat a b =
    if debug then Printf.printf "genNonAdd.median, call get_distance\n%!";
    let dis, medseq = get_distance a b cost_mat in
    {seq = medseq; costs = make_cost dis },dis

let median_3 cost_mat parent mine child1 child2 =
    let dis, medseq = get_distance_3d parent child1 child2 cost_mat in
    if debug then Printf.printf "genNonAdd.median_3, res cost = %d\n%!" dis;
    { seq = medseq; costs = make_cost dis }
    (* fake median3
    let medpc1,dis1 = median cost_mat parent child1 
    and medpc2,dis2 = median cost_mat parent child2 in
    if debug then Printf.printf "median_3, cost_p_c1 = %d, cost_p_c2 = %d,%!"
    dis1 dis2;
    let newmedseq,newdis = 
    if dis1<dis2 then medpc1.seq,dis1
    else medpc2.seq,dis2
    in
    if debug then Printf.printf "res cost = %d\n%!" newdis;
    { seq = newmedseq; costs = make_cost newdis } 
    *)


(*get closest code for code2 based on code1*)
let get_closest_code alph cost_mat code1 code2 =
    let comblst2 = Alphabet.find_codelist code2 alph in
    let comblst1 = Alphabet.find_codelist code1 alph in
    let mincost,bestc2 = List.fold_left (fun (acc_cost,acc_c2) c2 ->
        let lowest_cost = 
        List.fold_left (fun lower_cost c1 ->
            let thiscost = get_cost c1 c2 cost_mat in
            if thiscost < lower_cost then thiscost
            else lower_cost
        ) Utl.large_int comblst1 in
        if lowest_cost < acc_cost then lowest_cost,c2 
        else acc_cost,acc_c2
    ) (Utl.large_int,(-1)) comblst2 in
    if debug then Printf.printf "get_closest_code, mincost=%d, bestc2 = %d\n%!"
    mincost bestc2;
    bestc2

(**[to_single alph cost_mat parent mine] return single assignment of mine based
* on parent. return cost between parent and new single*)
let to_single alph cost_mat parent mine = 
    if debug then Printf.printf "genNonAdd.to_single,%!";
    let arrp = Sequence.to_array parent.seq
    and arrm = Sequence.to_array mine.seq in
    let arrclosest = Array_ops.map_2 (fun codep codem ->
        get_closest_code alph cost_mat codep codem
    ) arrp arrm in
    let newsingle = { mine with seq = Sequence.of_array arrclosest } in
    let cst = distance parent newsingle cost_mat in
    if debug then Printf.printf "cst = %d\n%!" cst;
    newsingle, cst


let compare a b =
    Sequence.compare a.seq b.seq
