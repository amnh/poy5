type genome
type genome_arr
external c_get_num_genome : genome_arr -> int = "get_num_genome_CAML"
external c_get_num_gene : genome_arr -> int = "get_num_gene_CAML"
external c_print_genome : genome -> int -> unit = "print_genome_CAML"
external c_print_genome_arr : genome_arr -> int -> int -> unit
  = "print_genome_arr_CAML"
external c_get_one_genome : genome_arr -> int -> genome
  = "get_one_genome_CAML"
external c_cmp_inv_dis : genome -> genome -> int -> int -> int
  = "cmp_inv_dis_CAML"
external c_create_empty_genome_arr : int -> int -> genome_arr
  = "create_empty_genome_arr_CAML"
external c_set : genome_arr -> int -> int -> int -> unit = "set_CAML"
