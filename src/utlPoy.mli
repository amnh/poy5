val fprintf : out_channel -> ('a, out_channel, unit) format -> 'a
val is_existed_code : int -> Sequence.s -> bool
val is_existed_char : string -> Sequence.s -> bool
val printDNA : Sequence.s -> unit
val create_gap_seq : ?gap:int -> int -> Sequence.s
val cmp_num_all_DNA : Sequence.s -> int
val cmp_num_not_gap : Sequence.s -> int
val cmp_gap_cost : int * int -> Sequence.s -> int
val cmp_ali_cost :
  Sequence.s -> Sequence.s -> [> `Positive ] -> Cost_matrix.Two_D.m -> int
val get_empty_seq : unit -> Sequence.s
val subseq : Sequence.s -> int -> int -> Sequence.s
val align2 :
  Sequence.s ->
  Sequence.s -> Cost_matrix.Two_D.m -> Sequence.s * Sequence.s * int * int
val align3 :
  Sequence.s ->
  Sequence.s ->
  Sequence.s ->
  Cost_matrix.Three_D.m -> Sequence.s * Sequence.s * Sequence.s * int * int
val closest_alied_seq :
  Sequence.s -> Sequence.s -> Cost_matrix.Two_D.m -> Sequence.s * int
val concat : Sequence.s list -> Sequence.s
val create_subalign2 :
  Sequence.s ->
  Sequence.s ->
  Cost_matrix.Two_D.m ->
  int -> int -> int -> int -> Sequence.s * Sequence.s * int
val dna_gap : int
val get_num_base : Sequence.s -> int
val delete_gap : ?gap_code:int -> Sequence.s -> Sequence.s
val create_median_gap :
  Sequence.s ->
  ?start_pos:int -> ?end_pos:int -> Cost_matrix.Two_D.m -> Sequence.s
val create_median_seq :
  ?approx:[< `BothSeq | `First | `Second > `BothSeq ] ->
  Sequence.s -> Sequence.s -> Cost_matrix.Two_D.m -> Sequence.s * int
val create_median_deled_seq : Sequence.s -> 'a -> Sequence.s
val create_median :
  ?approx:[< `BothSeq | `First | `Second > `BothSeq ] ->
  Sequence.s ->
  Sequence.s ->
  ?s1:int ->
  ?e1:int ->
  ?s2:int ->
  ?e2:int ->
  Cost_matrix.Two_D.m -> Sequence.s * Sequence.s * Sequence.s * int
val check_repeated_char : Sequence.s -> Alphabet.a -> unit
val create_general_ali :
  int array ->
  int array -> int -> Cost_matrix.Two_D.m -> int array * int array * int
val map : (int -> int) -> Sequence.s -> Sequence.s
val of_array : int array -> Sequence.s
val test_general_ali : unit -> unit
val get_single_seq : Sequence.s -> Cost_matrix.Two_D.m -> Sequence.s
val cmp_locus_indel_cost :
  Sequence.s -> Cost_matrix.Two_D.m -> int * int -> int
