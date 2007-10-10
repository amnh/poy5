exception Illegal_Arguments
module IntMap :
  sig
    type key = int
    type 'a t = 'a All_sets.IntegerMap.t
    val empty : 'a t
    val is_empty : 'a t -> bool
    val add : key -> 'a -> 'a t -> 'a t
    val find : key -> 'a t -> 'a
    val remove : key -> 'a t -> 'a t
    val mem : key -> 'a t -> bool
    val iter : (key -> 'a -> unit) -> 'a t -> unit
    val map : ('a -> 'b) -> 'a t -> 'b t
    val mapi : (key -> 'a -> 'b) -> 'a t -> 'b t
    val fold : (key -> 'a -> 'b -> 'b) -> 'a t -> 'b -> 'b
    val compare : ('a -> 'a -> int) -> 'a t -> 'a t -> int
    val equal : ('a -> 'a -> bool) -> 'a t -> 'a t -> bool
  end
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
exception No_Union
type t =
    SeqCS of SeqCS.t
  | BreakinvCS of BreakinvCS.t
  | ChromCS of ChromCS.t
  | AnnchromCS of AnnchromCS.t
  | GenomeCS of GenomeCS.t
type u = U_SeqCS of SeqCS.u | U_Others
val failwith_todo : string -> 'a
val total_cost : t -> float
val total_recost : t -> float
val subtree_recost : t -> float
val alpha : t -> Alphabet.a
val c2 : t -> Cost_matrix.Two_D.m
val chrom_pam : t -> Data.dyna_pam_t
val state : t -> Data.dyna_state_t
val code : t -> int
val leaf_sequences : t -> Sequence.s array SeqCS.Codes.t
val unions : u -> Sequence.Unions.u SeqCS.Codes.t
val reprioritize : t -> t -> t
val prioritize : t -> t
val to_union : t -> u
val union : t -> u -> u -> u
val cardinal_union : u -> int
val poly_saturation : u -> int -> float
val of_array :
    Data.dynamic_hom_spec ->
        (Sequence.s Data.dyna_data * AnnchromCS.IntMap.key) array ->
            int -> int -> int -> t
val of_list :
    Data.dynamic_hom_spec ->
        (Sequence.s Data.dyna_data * AnnchromCS.IntMap.key) list ->
            int -> int -> int -> t
val median : t -> t -> t

(** [to_single ?is_root pre_ref_code alied_map p n] returns a node that contains per character a single state
 * which is closest to [p] among those available in [n]. Useful for tree length
 * verification. is_root optional paramter indicates that if n is root. The
 * default is false. pre_ref_code contains active codes for chromosome characters. 
 * Inactive codes are eliminated from diagnosis. 
 * If p is the handle, alied_map is the root containing the aligned map between p
 * and n for chromosome stuff, else alied_map is assigned by p *)

val to_single : ChromCS.IntSet.t -> t option -> t -> t -> float * float * t

(** [to_single_root n] is the same as [to_single n n]. *)
val to_single_root : ChromCS.IntSet.t -> t -> float * float * t

(** [readjust ch1 ch2 par mine] attempts to (heuristically) readjust the character 
* set [mine] to somewhere in between [ch1], [ch2], and [par] (the children and
* parent of [mine] respectively). The function returns a triple [(a, b, c)],
* where [a] is the previous cost of [mine], [b] is the new cost of [c] as [ch1]
* and [ch2] parent, and [c] is the new readjusted [mine]. *)
val readjust : All_sets.Integers.t option -> All_sets.Integers.t -> t -> t -> t -> t -> 
    All_sets.Integers.t * float * float * t 
val median_3 : t -> t -> t -> t -> t
(* Like [distance] but calculates it only if the type of the characters match
* one of those listed. *)
val distance_of_type : 
    [`Seq | `Breakinv | `Chrom | `Annchrom | `Genome] list  -> t -> t -> float
val distance : t -> t -> float
val distance_union : u -> u -> float
val to_string : t -> string
val dist_2 : float -> t -> t -> t -> float
val f_codes : t -> All_sets.Integers.t -> t
val f_codes_comp : t -> All_sets.Integers.t -> t
val compare_data : t -> t -> int
val compare_union : u -> u -> int
val to_formatter :
  ChromCS.IntSet.t ->
  Tags.attribute list -> t -> t option -> Data.d -> Tags.output list
val tabu_distance : t -> t -> float
val get_active_ref_code : t -> ChromCS.IntSet.t * ChromCS.IntSet.t
val cardinal : t -> int
val get_sequence_union : SeqCS.Codes.key -> u -> Sequence.Unions.u option
val encoding :
  (float, Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array1.t ->
  t -> float
module Kolmogorov : sig val correct_cost : t -> Data.kolmo_spec -> t end

val print : t -> unit
val copy_chrom_map : t -> t -> t
