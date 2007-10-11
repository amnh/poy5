val infinity : int
val max_seq_len : int
val fprintf : out_channel -> ('a, out_channel, unit) format -> 'a
val gen_chrom_ref_code : int ref
val gen_seq_ref_code : int ref
val gen_genome_ref_code : int ref
val get_new_chrom_ref_code : unit -> int
val get_new_genome_ref_code : unit -> int
val get_new_seq_ref_code : unit -> int
val deref : 'a option -> 'a
val is_null : 'a option -> bool
val compare_non_dec_list : int list -> int list -> bool
val get_sum_arr : int array -> int -> int -> int
val invert_arr : 'a array -> unit
val invert_subarr : 'a array -> int -> int -> unit
val invert_direction_subarr : int array -> int -> int -> unit
val binary_search : int array -> int -> int
val find_index : 'a array -> 'b -> ('b -> 'a -> int) -> int
val insert : 'a array -> int -> 'a -> 'a array
val move_forward : 'a array -> int -> int -> int -> 'a array
val swap_item : int -> int -> 'a array -> 'a array
val printIntArr : int array -> unit
val printIntMat : int array array -> unit
val create_ls : int -> 'a -> 'a list
module IntSet :
  sig
    type elt = int
    type t = All_sets.Integers.t
    val empty : t
    val is_empty : t -> bool
    val mem : elt -> t -> bool
    val add : elt -> t -> t
    val singleton : elt -> t
    val remove : elt -> t -> t
    val union : t -> t -> t
    val inter : t -> t -> t
    val diff : t -> t -> t
    val compare : t -> t -> int
    val equal : t -> t -> bool
    val subset : t -> t -> bool
    val iter : (elt -> unit) -> t -> unit
    val fold : (elt -> 'a -> 'a) -> t -> 'a -> 'a
    val for_all : (elt -> bool) -> t -> bool
    val exists : (elt -> bool) -> t -> bool
    val filter : (elt -> bool) -> t -> t
    val partition : (elt -> bool) -> t -> t * t
    val cardinal : t -> int
    val elements : t -> elt list
    val min_elt : t -> elt
    val max_elt : t -> elt
    val choose : t -> elt
    val split : elt -> t -> t * bool * t
  end
val remove_nth : ?acc:'a list -> 'a list -> int -> 'a * 'a list
val insert_arr : int array -> int array -> int -> int array
val pairwisep : ('a -> 'a -> bool) -> 'a list -> bool
val filterArr : 'a array -> ('a -> bool) -> 'a array
val get_k_random_elem : 'a list -> int -> 'a list
val equalArr : 'a array -> 'b array -> ('a -> 'b -> int) -> bool
val filterArray : ('a -> bool) -> 'a array -> 'a array
val break_array : 'a array -> (int * int) list -> 'a array list
val printIntSet : IntSet.t -> unit
val get_dir : [> `Negative | `Positive ] -> string
val max_arr : 'a array -> 'a
val min_arr : 'a array -> 'a
val get_common :
  'a array -> 'a array -> ('a -> 'a -> int) -> 'a array * 'a array

