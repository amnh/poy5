(* POY 5.0 Beta. A phylogenetic analysis program using Dynamic Homologies.    *)
(* Copyright (C) 2013 Andrés Varón, Lin Hong, Nicholas Lucaroni, Ward Wheeler,*)
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


(*this module is for prealigned custom alphabet datatype. *)

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

let (-->) a b = b a

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

let get_worst_cost = Cost_matrix.Two_D.worst_cost

let get_median = Cost_matrix.Two_D.median

let get_cost3d = Cost_matrix.Three_D.cost

let get_median3d = Cost_matrix.Three_D.median

let is_empty seq cst = Sequence.is_empty seq (Cost_matrix.Two_D.gap cst)


(**[get_distance gnoadd1 gnoadd2 cost_mat] return distance between two sequence.
* also return a median sequence.
* remember that these sequence are prealigned. just add cost from each column.*)
let get_distance gnoadd1 gnoadd2 cost_mat =
    let seq1,seq2 = gnoadd1.seq,gnoadd2.seq in
    if is_empty seq1 cost_mat then
        0,seq2
    else if is_empty seq2 cost_mat then
        0,seq1
    else begin
        assert ((Sequence.length seq1) = (Sequence.length seq2) );
        let arr1,arr2 = Sequence.to_array seq1, Sequence.to_array seq2 in
        let rescost = ref 0 in
        let medarr = match Cost_matrix.Two_D.get_cost_model cost_mat with
            | Cost_matrix.No_Alignment
            | Cost_matrix.Linnear ->
                Array_ops.map_2
                    (fun code1 code2 ->
                        let c = get_cost code1 code2 cost_mat in
                        rescost := !rescost + c;
                        get_median code1 code2 cost_mat)
                    arr1 arr2
            | Cost_matrix.Affine _ -> assert false
        in
        !rescost, Sequence.of_array medarr
    end


(**[get_distance gnoadd1 gnoadd2 cost_mat] return distance between two sequence.
* also return a median sequence.
* remember that these sequence are prealigned. just add cost from each column.*)
let get_distance_3d gnoadd1 gnoadd2 gnoadd3 cost_mat cost_mat2 =
    let seq1,seq2,seq3 = gnoadd1.seq,gnoadd2.seq,gnoadd3.seq in
    match (is_empty seq1 cost_mat2),(is_empty seq2 cost_mat2),(is_empty seq3 cost_mat2) with
    | true, true, _    -> 0, seq3
    | true, _, true    -> 0, seq2
    | _, true, true    -> 0, seq1
    | true,false,false -> get_distance gnoadd2 gnoadd3 cost_mat2
    | false,true,false -> get_distance gnoadd1 gnoadd3 cost_mat2
    | false,false,true -> get_distance gnoadd1 gnoadd2 cost_mat2
    | false,false,false->
        assert( ((Sequence.length seq1) = (Sequence.length seq2))
             && ((Sequence.length seq3) = (Sequence.length seq2)) );
        let arr1 = Sequence.to_array seq1 in
        let arr2 = Sequence.to_array seq2 in
        let arr3 = Sequence.to_array seq3 in
        let rescost = ref 0 in
        let medarr =
            Array_ops.map_3
                (fun code1 code2 code3 ->
                    let c = get_cost3d code1 code2 code3 cost_mat in
                    rescost := !rescost + c;
                    get_median3d code1 code2 code3 cost_mat)
                arr1 arr2 arr3
        in
        !rescost, Sequence.of_array medarr


(** [distance gnoadd1 gnoadd2 cost_mat] just return the distance between two
* general nonaddictive sequence. ignore the median seq.*)
let distance gnoadd1 gnoadd2 cost_mat =
    let dis,_ = get_distance gnoadd1 gnoadd2 cost_mat in
    dis


let get_max_distance gnoadd1 gnoadd2 cost_mat = 
    let seq1,seq2 = gnoadd1.seq,gnoadd2.seq in
    if (is_empty seq1 cost_mat) then
        0
    else if (is_empty seq2 cost_mat) then
        0
    else begin
        assert( (Sequence.length seq1) = (Sequence.length seq2) );
        let arr1,arr2 = Sequence.to_array seq1, Sequence.to_array seq2 in
        Array_ops.fold_right_2
            (fun acc code1 code2 -> acc + (get_worst_cost code1 code2 cost_mat))
            0 arr1 arr2
    end


let max_distance gnoadd1 gnoadd2 cost_mat =
    get_max_distance gnoadd1 gnoadd2 cost_mat


(** [median cost_mat a b] return median of two general nonaddictive sequence*)
let median cost_mat a b =
    let dis, medseq = get_distance a b cost_mat in
    {seq = medseq; costs = make_cost dis },dis

(*when there is no 3d matrix, what do we do for median3 function? 
*  we align both child and parent, then return the better one *)
let median_3_fake cost_mat parent mine child1 child2 =
    let medpc1,dist1 = median cost_mat parent child1 in
    let medpc2,dist2 = median cost_mat parent child2 in
    if dist1<dist2 then
        medpc1
    else
        medpc2
    
(* [median_3] determine the cost of aligning the three sequences *)
let median_3 cost_mat3 cost_mat2 parent mine child1 child2 =
    let dis, medseq = get_distance_3d parent child1 child2 cost_mat3 cost_mat2 in
    { seq = medseq; costs = make_cost dis }


(* get closest code for code2 based on code1 *)
let get_closest_code alph cost_mat code1 code2 =
    let comblst2 = Alphabet.find_codelist code2 alph in
    let comblst1 = Alphabet.find_codelist code1 alph in
    let mincost,bestc2 =
        List.fold_left
            (fun (acc_cost,acc_c2) c2 ->
                let lowest_cost = 
                    List.fold_left
                        (fun lower_cost c1 ->
                            let thiscost = get_cost c1 c2 cost_mat in
                            if thiscost < lower_cost then thiscost
                            else lower_cost)
                        Utl.large_int comblst1
                in
                if lowest_cost < acc_cost then lowest_cost,c2 
                else acc_cost,acc_c2)
            (Utl.large_int,(-1))
            comblst2
    in
    bestc2

(**[to_single alph cost_mat parent mine] return single assignment of mine based
* on parent. return cost between parent and new single*)
let to_single alph cost_mat parent mine =
    match is_empty parent.seq cost_mat, is_empty mine.seq cost_mat with
    | true, _ -> mine,0
    | _,true  -> parent,0
    | false,false ->
        let arrp = Sequence.to_array parent.seq
        and arrm = Sequence.to_array mine.seq in
        let arrclosest =
            Array_ops.map_2
                (fun codep codem -> get_closest_code alph cost_mat codep codem)
                arrp arrm
        in
        let newsingle = { mine with seq = Sequence.of_array arrclosest } in
        let cst = distance parent newsingle cost_mat in
        newsingle, cst

let compare a b =
    Sequence.compare a.seq b.seq
