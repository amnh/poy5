(* POY 4.0 Beta. A phylogenetic analysis program using Dynamic Homologies.    *)
(* Copyright (C) 2007  Andrés Varón, Le Sy Vinh, Illya Bomash, Ward Wheeler,  *)
(* and the American Museum of Natural History.                                *)
(*                                                                            *)
(* This program is free software; you can redistribute it and/or modify       *)
(* it under the terms of the GNU General Public License as published by       *)
(* the Free Software Foundation; either version 2 of the License, or          *)
(* (at your option) any later version.                                        *)
(*                                                                            *)
(* This program is distributed in the hope that it will be useful,            *)
(* but WITHOUT ANY WARRANTY; without even the implied warranty of             *)
(* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the              *)
(* GNU General Public License for more details.                               *)
(*                                                                            *)
(* You should have received a copy of the GNU General Public License          *)
(* along with this program; if not, write to the Free Software                *)
(* Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301   *)
(* USA                                                                        *)

(** A dynamic character set implementation. 
* The dynamic character set allows rearrangements *)

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


(** A dynamic character type. 'a can be SeqCS, ChromCS, BreakCS, GenomeCS ...*)
type t = 
    | SeqCS  of SeqCS.t
    | MlCS of MlDynamicCS.t
    | BreakinvCS of  BreakinvCS.t 
    | ChromCS of ChromCS.t 
    | AnnchromCS of AnnchromCS.t
    | GenomeCS of GenomeCS.t


type u = U_SeqCS of SeqCS.Union.u | U_Others
val failwith_todo : string -> 'a

(** [total_cost a] returns the total cost to create
*  dynamic character set [a]*)
val total_cost : t -> float

(** [total_recost a] returns the total recost to create
*  dynamic character set [a]*)
val total_recost : t -> float

(** [subtree_recost a] returns the total recost of
* subtree whose root contains dynamic character set [a] *)
val subtree_recost : t -> float

(** [alpha a] returns the alphabet of dynamic character set [a] *)
val alpha : t -> Alphabet.a

(** [c2 a] returns the two dimentional cost matrix
* of dynamic character set [a] *)
val c2 : t -> Cost_matrix.Two_D.m

(** [chrom_pam a] returns the user-defined chromosome parameters
* of dynamic character set [a] *)
val chrom_pam : t -> Data.dyna_pam_t

(** [state a] returns the character type of
* dynamic character set [a] *)
val state : t -> Data.dyna_state_t

(** [code a] returns the code of dynamic character set [a] *)
val code : t -> int

(** [leaf_sequences a] turns dynamic character set [a] 
* into a set of chromosome arrays *)
val leaf_sequences : t -> 
    [ `DO of Sequence.s | `First of Sequence.s | `Last of Sequence.s ] 
    array SeqCS.Codes.t

(** [unions a] returns the union of dynamic character set [a] *)
val unions : u -> SeqCS.union_element option SeqCS.Codes.t

val to_union : t -> u
val union : t -> u -> u -> u
val cardinal_union : u -> int
val poly_saturation : u -> int -> float

(** combine a dynamic and static character set *)
val combine : t -> MlStaticCS.t -> t

(** [of_array spec genome_arr code taxon num_taxa] 
* creates a dynamic character set from genome array [genome_arr] *)
val of_array :
    Data.dynamic_hom_spec ->
        (Sequence.s Data.dyna_data * AnnchromCS.IntMap.key) array ->
            int -> int -> int -> t

(** [of_list spec genome_arr code taxon num_taxa] 
* creates a dynamic character set from genome list [genome_ls] *)
val of_list :
    Data.dynamic_hom_spec ->
        (Sequence.s Data.dyna_data * AnnchromCS.IntMap.key) list ->
            int -> int -> int -> t

(** [median a b] creates the median set between dynamic 
* character sets [a] and [b] *)
val median : int -> t -> t -> float option -> float option -> t

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
val readjust : [`ThreeD of int option | `ApproxD of int option ] -> All_sets.Integers.t option -> All_sets.Integers.t -> t -> t -> t -> t -> 
    All_sets.Integers.t * float * float * t 

(** [median_3 p n c1 c2] creates the median among
* three dynamic character sets [p], [c1] and [c2] *)
val median_3 : t -> t -> t -> t -> t

(* Like [distance] but calculates it only 
* if the type of the characters match one of those listed. *)
val distance_of_type : 
    [`Seq | `Ml | `Breakinv | `Chrom | `Annchrom | `Genome] list  -> float -> t -> t -> float

(** [distance_of_type a b] returns the distance between
* two dynamic character sets [a] and [b] *)
val distance : float -> t -> t -> float

(** [distance_union a b] returns the union distance between
* two dynamic character sets [a] and [b] *)
val distance_union : u -> u -> float

(** [to_string a] returns dynamic character set [a] 
* into the string format *)
val to_string : t -> string

(* name of the character type *)
val name_string : t -> string

(* [dist_2 delta n a b] calculates the cost of joining 
* the node containing  chromosome character set [n]
* between two nodes containing [a] and [b]. 
* [a] must be the parent (ancestor) of [b] *)
val dist_2 : float -> t -> t -> t -> float

(** [f_codes s c] returns a dynamic character subset 
* of  dynamic character set [s] whose codes are 
* also in  chromosome character set [c] *)
val f_codes : t -> All_sets.Integers.t -> t

(** [f_codes_comp s c] returns a dynamic character subset 
* of  dynamic character set [s] whose codes are NOT
* also in  chromosome character set [c] *)
val f_codes_comp : t -> All_sets.Integers.t -> t

(** [compare_data a b] returns 1 if dynamic character 
* set [a] is the same dynamic character set [b], 
* otherwise (-1) or 1 *)
val compare_data : t -> t -> int

val compare_union : u -> u -> int

(** [to_formatter ref_codes attr t parent_t] returns
* dynamic character set [t] into Tag.output format *) 
val to_formatter :
  ChromCS.IntSet.t ->
  Xml.attribute list -> t -> t option -> Data.d -> Xml.xml Sexpr.t list

(** [tabu_distance a_final b_final] returns the 
* tabu distance between dynamic character set [a_final] and [b_final] *)
val tabu_distance : t -> t -> float

(** [get_active_ref_code t] returns the set of active codes
* of dynamic character set [t] *)
val get_active_ref_code : t -> ChromCS.IntSet.t * ChromCS.IntSet.t

val cardinal : t -> int
val get_sequence_union : SeqCS.Codes.key -> u -> SeqCS.union_element option

val encoding :
  (float, Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array1.t ->
  t -> float

val print : t -> unit

(** [copy_chrom_map s_ch d_ch] copies the choromosome map
* of dynamic character set [s_ch] to [d_ch] *)
val copy_chrom_map : t -> t -> t
