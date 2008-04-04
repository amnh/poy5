type pairChromPam_t = ChromPam.chromPairAliPam_t
type seed_t = {
  mutable id : int;
  mutable sta1 : int;
  mutable sta2 : int;
  mutable len : int;
  mutable score : int;
  mutable is_chained : bool;
}
type sufTree_t = SufTree.sufTree_t
type incList_t = seed_t IncList.incList_t
val num_located_cell : int
val num_used_cell : int ref
val located_cell_arr : incList_t array ref
val fprintf : out_channel -> ('a, out_channel, unit) format -> 'a
val deref : 'a option -> 'a
val copy : seed_t -> seed_t
val swap : seed_t -> unit
val create_empty_seed : unit -> seed_t
val get_dum_first_seed : int -> int -> seed_t
val get_dum_last_seed : int -> int -> seed_t
val set_seed : seed_t -> int -> int -> int -> int -> bool -> unit
val create_new_seed : int -> int -> int -> int -> seed_t
val cmp_dis1 : seed_t -> seed_t -> int
val cmp_dis2 : seed_t -> seed_t -> int
val cmp_dis : seed_t -> seed_t -> int
val cmp_dia_dis : seed_t -> seed_t -> int
val invert : int -> int -> seed_t -> unit
val print : seed_t -> unit
val compare_end1 : seed_t -> seed_t -> int
val compare_start1 : seed_t -> seed_t -> int
val cmp_connecting_cost : seed_t -> seed_t -> pairChromPam_t -> int
val create_local_connection :
  seed_t list -> pairChromPam_t -> int array * int array * seed_t array
val locate_mem : unit -> unit
val create_cell : int -> int -> int -> bool -> incList_t
val connect :
  int ->
  int list ->
  incList_t -> pairChromPam_t -> int array -> int array -> seed_t list
val determine_pos_seed :
  sufTree_t -> Sequence.s -> pairChromPam_t -> seed_t list
val determine_neg_seed :
  sufTree_t -> Sequence.s -> pairChromPam_t -> seed_t list
val print_sta : seed_t list -> unit
val determine_seed :
  Sequence.s ->
  Sequence.s ->
  pairChromPam_t -> ChromPam.direction_t -> seed_t list * seed_t list
val create_alied_seed :
  seed_t ->
  Sequence.s ->
  Sequence.s -> Cost_matrix.Two_D.m -> Sequence.s * Sequence.s * int
val create_alied_subseq :
  seed_t ->
  seed_t ->
  Sequence.s ->
  Sequence.s -> Cost_matrix.Two_D.m -> Sequence.s * Sequence.s * int
