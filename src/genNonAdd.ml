(* POY 5.1.1. A phylogenetic analysis program using Dynamic Homologies.       *)
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
    weights : float array;
}

let (-->) a b = b a

let init_gnonadd_t in_seq weights =
  let weights = match weights with
    | None   -> Array.make (Sequence.length in_seq) 1.0
    | Some x -> x
  in
  { seq = in_seq; costs = { min = 0.0; max = 0.0 }; weights = weights; }

let get_max_cost ct = ct.max

let get_min_cost ct = ct.min

(*[make_cost tmpcost] make cost_tuple with tmpcost -- to do: what is maxcost?*)
let make_cost tmpcost = 
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
        0.0,seq2
    else if is_empty seq2 cost_mat then
        0.0,seq1
    else begin
        assert ((Sequence.length seq1) = (Sequence.length seq2) );
        let rescost = ref 0.0 in
        let medarr = match Cost_matrix.Two_D.get_cost_model cost_mat with
            | Cost_matrix.No_Alignment
            | Cost_matrix.Linnear ->
                Array.mapi
                    (fun i weight ->
                        let code1 = Sequence.get seq1 i
                        and code2 = Sequence.get seq2 i in
                        let c = float_of_int (get_cost code1 code2 cost_mat) in
                        rescost := !rescost +. (c *. weight);
                        get_median code1 code2 cost_mat)
                    gnoadd2.weights
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
    | true, true, _    -> 0.0, seq3
    | true, _, true    -> 0.0, seq2
    | _, true, true    -> 0.0, seq1
    | true,false,false -> get_distance gnoadd2 gnoadd3 cost_mat2
    | false,true,false -> get_distance gnoadd1 gnoadd3 cost_mat2
    | false,false,true -> get_distance gnoadd1 gnoadd2 cost_mat2
    | false,false,false->
        assert( ((Sequence.length seq1) = (Sequence.length seq2))
             && ((Sequence.length seq3) = (Sequence.length seq2)) );
        let rescost = ref 0.0 in
        let medarr =
            Array.mapi
                (fun i weight ->
                    let code1 = Sequence.get seq1 i
                    and code2 = Sequence.get seq2 i
                    and code3 = Sequence.get seq3 i in
                    let c = float_of_int (get_cost3d code1 code2 code3 cost_mat) in
                    rescost := !rescost +. (c *.weight);
                    get_median3d code1 code2 code3 cost_mat)
                gnoadd1.weights
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
        0.0
    else if (is_empty seq2 cost_mat) then
        0.0
    else begin
        assert((Sequence.length seq1) = (Sequence.length seq2));
        Sequence.foldi_2
            (fun acc i code1 code2 ->
                let ncost = float_of_int (get_worst_cost code1 code2 cost_mat) in
                acc +. (ncost *. gnoadd1.weights.(i)))
            0.0 seq1 seq2
    end


let max_distance gnoadd1 gnoadd2 cost_mat =
    get_max_distance gnoadd1 gnoadd2 cost_mat


(** [median cost_mat a b] return median of two general nonaddictive sequence*)
let median cost_mat a b =
    let dis, medseq = get_distance a b cost_mat in
    {seq = medseq; costs = make_cost dis; weights = a.weights},dis


(* [median_3] determine the cost of aligning the three sequences *)
let median_3 cost_mat3 cost_mat2 parent mine child1 child2 =
    let dis, medseq = get_distance_3d parent child1 child2 cost_mat3 cost_mat2 in
    { seq = medseq; costs = make_cost dis; weights = child1.weights; }


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
    | true, _ -> mine,0.0
    | _,true  -> parent,0.0
    | false,false ->
        let arrclosest =
            assert( (Sequence.length parent.seq) = (Array.length parent.weights) );
            Array.mapi
                (fun i weight ->
                    let codep = Sequence.get parent.seq i
                    and codem = Sequence.get mine.seq i in
                    get_closest_code alph cost_mat codep codem)
                parent.weights
        in
        let newsingle = { mine with seq = Sequence.of_array arrclosest } in
        let cst = distance parent newsingle cost_mat in
        newsingle, cst

(*when there is no 3d matrix, what do we do for median3 function? 
*  we align both child and parent, then return the better one *)
let median_3_fake alph m s1 s0 s2 s3 =
    let algn    s1 s2 = median  m s1 s2
    and closest s1 s2 = to_single alph m s1 s2
    and is_empty s1 m = is_empty s1.seq m in
    let make_center_with_missing a b : float * (gnonadd_sequence * float) * float =
        let sab, cab = algn a b in
        cab,(sab,0.0),cab
    and make_center s1 s2 s3 : float * (gnonadd_sequence * float) * float =
            (* first median  *)
        let s12, c12 = algn s1 s2 
        and s23, c23 = algn s2 s3 
        and s13, c13 = algn s1 s3 in
            (* second median *)
        let s123, c123 = algn s12 s3
        and s231, c231 = algn s23 s1
        and s132, c132 = algn s13 s2 in
            (* sum costs *)
        let c123 = c123 +. c12
        and c231 = c231 +. c23
        and c132 = c132 +. c13 in
            (* determine best... *)
        if c123 <= c231 then
            if c123 <= c132
                then c123, closest s3 s12, c123
                else c132, closest s2 s13, c123
        else if c231 < c132
            then c231, closest s1 s23, c123
            else c132, closest s2 s13, c123
    in
    let cst, (s, _), previous =
        match (is_empty s1 m), (is_empty s2 m), (is_empty s3 m) with
        | true, true, true  -> 0.0, (s1,0.0), 0.0
        | false, true, true -> 0.0, (s1,0.0), 0.0
        | true, false, true -> 0.0, (s2,0.0), 0.0
        | true, true, false -> 0.0, (s3,0.0), 0.0
        | false,false,true  -> make_center_with_missing s1 s2
        | false,true,false  -> make_center_with_missing s1 s3
        | true,false,false  -> make_center_with_missing s2 s3
        | false,false,false -> make_center s1 s2 s3
    in
  s,cst

let compare a b =
    Sequence.compare a.seq b.seq
