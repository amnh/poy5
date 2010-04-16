type genome
type genome_arr
val c_get_num_genome : genome_arr -> int 
val c_get_num_gene : genome_arr -> int 
val c_print_genome : genome -> int -> unit 
val c_print_genome_arr : genome_arr -> int -> int -> unit
val c_get_one_genome : genome_arr -> int -> genome
val c_cmp_inv_dis : genome -> genome -> int -> int -> int
val inversion_distance : genome -> genome -> int -> bool -> int
val c_create_empty_genome_arr : int -> int -> genome_arr
val c_set : int -> genome_arr -> int -> int -> int -> unit 
val inversions : genome -> genome -> int -> int -> (int * int) list 
val genomes : int array array -> genome_arr
val genes : int array array -> genome array
val get_med3_genome : int -> genome -> genome -> genome -> int -> int -> genome
val get_med3_arr : genome -> int -> int array
val get_delimiter_arr : genome -> int -> int array
val c_get_delimiter_num : genome ->  int
