
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

(** [init_gnonadd_t inseq] create a new gnoadd with input seq*)
val init_gnonadd_t : Sequence.s -> gnonadd_sequence

(**[to_single alph cost_mat parent mine] return single assignment of mine based
* on parent. return cost between parent and new single*)
val to_single : Alphabet.a -> Cost_matrix.Two_D.m -> gnonadd_sequence ->
    gnonadd_sequence -> gnonadd_sequence * int

(**[distance gnoadd1 gnoadd2 cost_mat] return distance between two sequence.*)
val distance : gnonadd_sequence -> gnonadd_sequence -> Cost_matrix.Two_D.m -> int

(** [median cost_mat a b] return median of two general nonaddictive sequence*)
val median : Cost_matrix.Two_D.m -> gnonadd_sequence -> gnonadd_sequence ->
    gnonadd_sequence * int

(** [median_3 cost_mat parent mine child1 child2] return a new median for mine
* with info from parent,child1 and child2 *)
val median_3 : Cost_matrix.Two_D.m -> gnonadd_sequence -> gnonadd_sequence ->
    gnonadd_sequence -> gnonadd_sequence -> gnonadd_sequence

(** [compare gnoadd1 gnoadd2] compare the sequence of two input.*)
val compare : gnonadd_sequence -> gnonadd_sequence -> int

(** [get_max_cost cost_tuple] return the max of cost_tuple*)
val get_max_cost : cost_tuple -> float
