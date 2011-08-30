

type cost_tuple = 
{
    min : float;
    max : float;
}

type gnonadd_sequence = {
    seq : Sequence.s;
    costs : cost_tuple;
}


val init_gnonadd_t : Sequence.s -> gnonadd_sequence

val to_single : Alphabet.a -> Cost_matrix.Two_D.m -> gnonadd_sequence ->
    gnonadd_sequence -> gnonadd_sequence * int

val distance : gnonadd_sequence -> gnonadd_sequence -> Cost_matrix.Two_D.m -> int

val median : Cost_matrix.Two_D.m -> gnonadd_sequence -> gnonadd_sequence ->
    gnonadd_sequence * int

val median_3 : Cost_matrix.Two_D.m -> gnonadd_sequence -> gnonadd_sequence ->
    gnonadd_sequence -> gnonadd_sequence -> gnonadd_sequence

val compare : gnonadd_sequence -> gnonadd_sequence -> int

val get_max_cost : cost_tuple -> float
