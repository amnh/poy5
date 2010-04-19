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

type t  
IFDEF USE_LIKELIHOOD THEN
type s  (* abstract type: contains matrix of character codes *)

(** two functions to convert from double*,int* to Array3,Array2 *)
external s_bigarray: 
    s -> (float,Bigarray.float64_elt,Bigarray.c_layout) Bigarray.Array3.t *
    ((int32,Bigarray.int32_elt,Bigarray.c_layout) Bigarray.Array1.t) option =
    "likelihood_CAML_StoBigarray"
external bigarray_s: 
    (float,Bigarray.float64_elt,Bigarray.c_layout) Bigarray.Array3.t ->
    ((int32,Bigarray.int32_elt,Bigarray.c_layout) Bigarray.Array1.t) option -> s =
    "likelihood_CAML_BigarraytoS"

(** [median_*** pi U D [Ui] at bt achars bchars rates probs] -> cchars
 * calculate the median of two nodes with character sets [achars] and [bchars]
 * and branch lenghts [at] and [bt], respectively. The probability matrix for
 * the transition is constructed via [U] [D] and, possibly, [Ui] (GTR only). *)
external median_gtr: FMatrix.m ->
    (float,Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array2.t -> 
    (float,Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array2.t ->
    (float,Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array2.t ->
     float -> float -> s -> s -> 
    (float,Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array1.t ->
    (float,Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array1.t -> s = 
     "likelihood_CAML_median_gtr" "likelihood_CAML_median_wrapped_gtr" 
external median_sym: FMatrix.m ->
    (float,Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array2.t -> 
    (float,Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array2.t ->
     float -> float -> s-> s -> 
    (float,Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array1.t ->
    (float,Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array1.t -> s = 
     "likelihood_CAML_median_sym" "likelihood_CAML_median_wrapped_sym" 
(** [readjust_*** U D [Ui] a b c at bt %inv rates probs priors ll] -> at*bt*ll
 * readjusts the branch lengths for the children of [c], [a] and [b], with
 * branch lengths [at] and [bt], respectively, using priors given by [priors],
 * and the gamma rates given by [rates]. The old likelihood is given in
 * [ll], and the new one is returned. [c] is updated to the new vector based on
 * the new branch length, returned with the new [ll]. *)
external readjust_sym: FMatrix.m ->
    (float,Bigarray.float64_elt,Bigarray.c_layout) Bigarray.Array2.t ->
    (float,Bigarray.float64_elt,Bigarray.c_layout) Bigarray.Array2.t ->
    s -> s -> s -> float -> float -> float ->
    (float,Bigarray.float64_elt,Bigarray.c_layout) Bigarray.Array1.t ->
    (float,Bigarray.float64_elt,Bigarray.c_layout) Bigarray.Array1.t ->
    (float,Bigarray.float64_elt,Bigarray.c_layout) Bigarray.Array1.t ->
    (float,Bigarray.float64_elt,Bigarray.c_layout) Bigarray.Array1.t ->
    float -> float*float =
        "likelihood_CAML_readjust_sym" "likelihood_CAML_readjust_sym_wrapped"
external readjust_gtr: FMatrix.m ->
    (float,Bigarray.float64_elt,Bigarray.c_layout) Bigarray.Array2.t ->
    (float,Bigarray.float64_elt,Bigarray.c_layout) Bigarray.Array2.t ->
    (float,Bigarray.float64_elt,Bigarray.c_layout) Bigarray.Array2.t ->
    s -> s -> s -> float -> float -> float ->
    (float,Bigarray.float64_elt,Bigarray.c_layout) Bigarray.Array1.t ->
    (float,Bigarray.float64_elt,Bigarray.c_layout) Bigarray.Array1.t ->
    (float,Bigarray.float64_elt,Bigarray.c_layout) Bigarray.Array1.t ->
    (float,Bigarray.float64_elt,Bigarray.c_layout) Bigarray.Array1.t ->
    float -> float*float =
        "likelihood_CAML_readjust_gtr" "likelihood_CAML_readjust_gtr_wrapped"

(** [loglikelihood s w pi prob %invar] -> float : calculates the mle of a character set *) 
external loglikelihood: 
    s -> (float,Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array1.t ->
         (float,Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array1.t ->
         (float,Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array1.t ->
        float -> float = "likelihood_CAML_loglikelihood"
(** [filter s as] -> s
 * filters s with indexes of as and returns new character set, *)
external filter: s -> int array -> s = "likelihood_CAML_filter"
(** [compare_chars s s] -> int
 * compares two sets, returns -1,0,1 depending on which set is larger *)
external compare_chars: s -> s -> int = "likelihood_CAML_compare"

(* report the minimum branch length allowed *)
external minimum_bl: unit -> float = "likelihood_CAML_minimum_bl"

(** calculates the proporation of similarity between two sequences *)
external proportion: s -> s -> float = "likelihood_CAML_proportion"

(** [estimate_time a b ] -> time
* estimates the time between two nodes using proporation (above) *)
val estimate_time : t -> t -> float * float

(** [gc_alloc_max] -> how many nodes to alloc before a GC is triggered *)
external gc_alloc_max : int -> unit = "likelihood_GC_custom_max"

(** [register] -> ()
 * register the likelihood operations for the garbage collection deserialization *)
external register : unit -> unit = "likelihood_CAML_register"

(** A string representation of the character set, used only for debugging purposes *)
val to_string : t -> string

(** [median a b at bt]
* computes the median between [a] and [b] with branch lengths [at] and [bt] *)
val median : t -> t -> float -> float -> int -> int -> t

(** The cost, likelihood, of the median [t]. wrapper around internally
* saved value *)
val median_cost : t -> float

(** [median_3 pn nn c1 c2]
* computes whatever heuristic values will be assigned to [Node.final] in
* a vertex, for the vertex [nn] with parent [pn] and children [c1] and [c2]. *)
val median_3 : t -> t -> t -> t -> t


val print_barray1 : (float,Bigarray.float64_elt,Bigarray.c_layout) Bigarray.Array1.t -> unit
val print_barray2 : (float,Bigarray.float64_elt,Bigarray.c_layout) Bigarray.Array2.t -> unit
val print_barray3 : (float,Bigarray.float64_elt,Bigarray.c_layout) Bigarray.Array3.t -> unit

(** [reroot_median a b] 
* computes the median that should be assigned to the root of a tree as the
* median between [a] and [b] with total branch length between the two nodes,
* [at]+[bt] in an unrooted or rooted tree *)
val reroot_median : t -> t -> float -> float -> t

(** [dist_2 a b c at bt ct]  
* calculates the Likelihood of connecting the root of a subtree [a] in between
* the pair of neighbor vertices [b] and [c]. This is used for fast evaluation
* during SPR and TBR. *)
val dist_2 : t -> t -> t -> float -> float -> float -> float option -> float

(** [f_codes x c]
* creates a new character set where all characters with code appearing in [c]
* have been filtered out. *)
val f_codes : t -> All_sets.Integers.t -> t

(** [f_codes_comp x c]
* creates a new character set where all characters with code NOT appearing in
* [c] have been filtered out (the complement of [f_codes]).*)
val f_codes_comp : t -> All_sets.Integers.t -> t

(** [cardinal x] returns the cardinality of the character set [x].*)
val cardinal : t -> int

(** [union prev ch1 ch2] *)
val union: t -> t -> t -> t

(** [compare_data a b]
* is a total ordering of the character sets, where [compare a b < 0]
* iff [a < b], [compare a b = 0] iff [a = b], otherwise [compare a b > 0]. *)
val compare_data : t -> t -> int

val compare : t -> t -> int

(** [readjust check has_changed c1 c2 mine t1 t2 ]
 *       -> modified set * old_mle * new_mle * (new_branch_lengths) * new node
 *
 *   Check/xopt and has_changed/x are not used, as they determine iteration on a
 * subset of characters, which doesn't make sense when branch lengths affect all
 * the characters in the set. Seperation of the characters should be done in
 * different complete character sets at the node level. 
 *   The function returns the set of modified characters (all of them), the old
 * likelihood, the new likelihood, the new branch lengths and an updated t
 * with new likelihood_vector and same model. *)
val readjust : All_sets.Integers.t option -> All_sets.Integers.t ->  
    t -> t -> t -> float -> float ->
    All_sets.Integers.t * float * float * (float*float) * t


(** [of_parser spec characters] creates a character set with specification
* [spec] and characters defined in the array [characters], where each element is
* a tuple [(states, code)], where [states] is the list of states observed in the
* characters with code [code]. If [states = None] then the character is missing
* (should be treated as if [states] held all the possible states for the
* character).*)
(*  [spec] -> [weights] -> [characters: ([states]*[code]) array] -> t *)
val of_parser : Nexus.File.static_spec -> float array -> ((Nexus.File.static_state * int) array) -> t

(* The extra cost incurred by the root of the tree. *)
val root_cost : t -> float

(* [distance a b at bt] computes the -log likelihood of [b] given [a] with
* branch lengths [at] and [bt]. *)
val distance : t -> t -> float -> float -> float

(* to be able to see the results on each vertex of the tree. *)
val to_formatter : Xml.attributes -> t -> float option * float option -> 
                        Data.d -> Xml.xml Sexpr.t list

val extract_states : t -> (float * int * MlModel.chars) list 



val get_codes : t -> int array

ELSE

val likelihood_error : string

END
