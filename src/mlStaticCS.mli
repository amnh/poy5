type t

val to_string : t -> string

(** [median prev a b] computes the median between [a] and [b] with [prev] as an 
* optional value of a previous median computed between [a] and [b] being [prev].
* [prev] can be used for heuristic purposes.
* *)
val median : t option -> t -> t -> t

(** The cost of the median [t]. This is not the overall subtree cost of the tree
* rooted by [t]. *)
val median_cost : t -> float

(** [median_3 pn nn c1 c2] computes whatever heuristic values will be assigned
* to [Node.final] in a vertex, for the vertex [nn] with parent [pn] and children
* [c1] and [c2]. *)
val median_3 : t -> t -> t -> t -> t

(** [reroot_median a b] computes the median that should be assigned to the root
* of a tree as the median between [a] and [b]. *)
val reroot_median : t -> t -> t

(** [dist_2 a b c]  estimates a lower bound of the additional cost of connecting
* the root of a subtree [a] in between the pair of neighbor vertices [b] and 
* [c]. This is used for fast evaluation during SPR and TBR. *)
val dist_2 : t -> t -> t -> float

(** [f_codes x c] creates a new character set where all characters with code
* appearing in [c] have been filtered out. *)
val f_codes : t -> All_sets.Integers.t -> t

(** [f_codes_comp x c] creates a new character set where all characters with
* code NOT appearing in [c] have been filtered out (the complement of
* [f_codes]).*)
val f_codes_comp : t -> All_sets.Integers.t -> t

(** [cardinal x] returns the cardinality of the character set [x].*)
val cardinal : t -> int

(** [compare_data a b] is a total ordering of the character sets, where 
* [compare a b < 0] iff [a < b], [compare a b = 0] iff [a = b], otherwise
* [compare a b > 0]. *)
val compare_data : t -> t -> int

(** [readjust check has_changed time c1 c2 par mine] readjusts the edge time for 
* some characters in the vertex [mine] with parent [par] and children [c1] and
* [c2], for an edge with overall [time]. If [check] is 
 * [None] then the function attempts to readjust the values in all the
 * characters of [mine] otherwise, only the codes included in [check] are
 * attempted to readjust. [has_changed] is the accumulator of characters that
 * have been adjusted in other character sets, and any character that is
 * effectively readjusted in [mine] should be added to this accumulator in the
 * output quintuple. The function outputs [(acc, prev_cost, cost, length, adjusted)]
 * where [acc] is the new accumulator from the [has_changed] set of codes,
 * [prev_cost] is the previous cost of the edge connecting [mine] and [par],
 * [cost] is the new cost of the edge connecting [mine] and [par] after the
 * readjustement, [length] is the new edge length (time), and [adjusted] is 
 * the newly adjusted vertex [mine]. *)
val readjust : All_sets.Integers.t option -> All_sets.Integers.t -> float -> t -> t -> t -> t ->
    All_sets.Integers.t * float * float * float * t

(** [of_parser spec characters] creates a character set with specification
* [spec] and characters defined in the array [characters], where each element is
* a tuple [(states, code)], where [states] is the list of states observed in the
* characters with code [code]. If [states = None] then the character is missing
* (should be treated as if [states] held all the possible states for the
* character).*)
val of_parser : Parser.SC.static_spec -> ((int list option * int) array) -> t

(* The extra cost incurred by the root of the tree. *)
val root_cost : t -> float

(* [distance a b] computes (in ML), the -log likelihood of [b] given [a]. *)
val distance : t -> t -> float
(** Non urgent functions, after finishing with the previous functions, 
* add the to_formatter to be able to see the results on each vertex of the tree.
val to_formatter : Tags.attributes -> t -> t option -> Data.d -> Tags.output
list
*)


