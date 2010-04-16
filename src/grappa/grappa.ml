

(* read file of format:
    * >S1
    * 1 2 3 4 5
    * >S2
    * 1 -3 -2 4 5 
    * ...
    * S1 and S2 are the taxon names always start with a >
    * numbers are the traits  *)
type genome
type genome_arr


external c_get_num_genome : genome_arr -> int = "grappa_CAML_get_num_genome"
external c_get_num_gene : genome_arr -> int = "grappa_CAML_get_num_gene"


(* takes a genome and the number of genes and print it *)
external c_print_genome : genome -> int -> unit = "grappa_CAML_print_genome"

(* takes a genome and the number of genes and print it *)
external c_print_genome_arr : genome_arr -> int -> int -> unit = "grappa_CAML_print_genome_arr"


(*give index in array and returns the one genome at that index *)
external c_get_one_genome : genome_arr -> int -> genome = "grappa_CAML_get_one_genome"

(*  for the distance returned *)
external c_cmp_inv_dis : genome -> genome -> int -> int ->  int 
=  "grappa_CAML_cmp_inv_dis"

(* for the median3, call albert-median resolver*)
(* genomeX, genomeY, genomeZ, number of genes in each genome, circular/or not,
* output the median of X, Y and Z *)
external c_inv_med : int -> genome -> genome -> genome -> int -> int ->
    genome_arr
  (*
    (int32, Bigarray.int32_elt, Bigarray.c_layout) Bigarray.Array1.t 
  *)
   = "grappa_CAML_inv_med_bytecode" "grappa_CAML_inv_med"

external c_get_gene_bigarr : genome -> int ->
    (int32, Bigarray.int32_elt, Bigarray.c_layout) Bigarray.Array1.t 
    = "grappa_CAML_get_gene_bigarr"

external c_get_delimiter_num : genome -> int =
    "grappa_CAML_get_delimiter_num"

external c_get_delimiter_bigarr : genome -> int ->
    (int32, Bigarray.int32_elt, Bigarray.c_layout) Bigarray.Array1.t 
    = "grappa_CAML_get_delimiter_bigarr"

external c_init : int -> unit = "grappa_CAML_initialize" 

let _ =
    c_init(64)

let bigarr_to_intarr bigarr = 
    let len = Bigarray.Array1.dim bigarr in
    let intarr = Array.init len ( 
        fun index ->
           Int32.to_int (bigarr.{index})
    ) in
    intarr

let get_med3_arr (g1 : genome) num_gene = 
    let med3_bigarr =  c_get_gene_bigarr g1 num_gene in
    bigarr_to_intarr med3_bigarr

let get_delimiter_arr (g1 : genome) num_delimiter = 
    let deli_bigarr =  c_get_delimiter_bigarr g1 num_delimiter in
    bigarr_to_intarr deli_bigarr

let get_med3_genome (medsov : int) (g0 : genome) (g1 : genome) (g2 : genome) num_gen circular =
    let output_genome_arr = c_inv_med medsov g0 g1 g2 num_gen circular in
    c_get_one_genome output_genome_arr 0 
  (*  let intarr_med3 = bigarr_to_intarr bigarr_med3 in
    intarr_med3 *)

let inversion_distance a b c d = 
    c_cmp_inv_dis a b c (if d then 1 else 0)

(* See [inversions] *)
external c_inversions : 
    genome -> genome -> int -> int -> (int * int) list  = 
        "grappa_CAML_inversions"

(** [inversions g1 g2 ngenes dist] produce a list of tuples, with a sequence of
* inversions that convert [g2] into [g1], with [ngenes] genes, and at inversion
* distance [dist]. *)
let inversions a b c d = List.rev (c_inversions a b c d)

external c_create_empty_genome_arr : int -> int -> genome_arr = "grappa_CAML_create_empty_genome_arr"

external c_set : int -> genome_arr -> int -> int -> int -> unit = "grappa_CAML_set"

let genomes arr = 
    (* Create an array of genomes *)
    let total = Array.length arr in
    if total = 0 then c_create_empty_genome_arr 0 0
    else 
        let genes = Array.length arr.(0) in
        let res = c_create_empty_genome_arr total genes in
        let () = 
            Array.iteri (fun a chrom ->
                Array.iteri (c_set 1 res a) chrom) arr
        in
        res

let genes arr = 
    let genomes = genomes arr in
    Array.init (Array.length arr) (c_get_one_genome genomes)
