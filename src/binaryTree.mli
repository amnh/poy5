
type ('k,'a) b_tree = 
    | Empty
    | Leaf of 'k * 'a 
    | Node of 'k * ('k,'a) b_tree * ('k,'a) b_tree

val just_a_leaf: ('k,'a) b_tree -> bool

val just_an_empty_leaf: ('k,'a) b_tree -> bool

val create_btree: 'k -> 'a -> ('k,'a) b_tree

val add_to_btree: 'k -> 'a -> ('k,'a) b_tree -> ('k -> unit)
-> ('a -> unit) -> ('a -> 'a -> bool) -> ('k,'a) b_tree * bool

val search_in_btree: 'k -> ('k,'a) b_tree -> ('k -> unit) -> ('a -> unit) -> 'a

val remove_from_btree: 'k -> ('k,'a) b_tree -> ('k -> unit) -> ('a -> unit) -> 
    ('k,'a) b_tree

val iter_b_tree: ('k,'a) b_tree -> ('a -> unit) -> ('k -> unit) -> bool -> unit

val print_b_tree: ('k,'a) b_tree -> ('a -> unit) -> ('k -> unit) -> unit
