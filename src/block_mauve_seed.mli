(*this is the seed table we are going to use.*)    
val palindromic_spaced_seed_tbl : (int,int list list) Hashtbl.t
   
val fill_in_hmatrix : Cost_matrix.Two_D.m -> unit

val get_proper_seedlen : float -> int

val return_a_seedNO : int -> (int array) ref -> unit

val get_a_seedNO : (int array) ref -> int

val radix_sort : (int*int*int*int array*int) array -> (int*int*int*int array*int) array

val extend_seq_in_both_dir : (int*int*int) array -> int -> int array array -> (int*int*int) array * int

val get_score_from_2seq : int array -> int array -> int -> int
