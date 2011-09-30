(** binary tree type. could be empty leafnode, or a leaf, or a tree.
*   k is the type of middle-key, a is the type of leaf.
*   middle key is bigger than keys in its left subtree*)
type ('k,'a) b_tree = 
    | Empty
    | Leaf of 'k * 'a 
    | Node of 'k * ('k,'a) b_tree * ('k,'a) b_tree

(** [just_a_leaf btree] returns true if the tree is a single non-empty leafnode*)
val just_a_leaf: ('k,'a) b_tree -> bool

(** [just_an_empty_leaf btree] returns true if the tree is empty*)
val just_an_empty_leaf: ('k,'a) b_tree -> bool

(** [create_btree key leafnode] create a single leafnode, return it as the tree*)
val create_btree: 'k -> 'a -> ('k,'a) b_tree

(** [add_to_btree key leafnode btree print_f_for_key compare_f_for_leafnode] add
* leafnode with key to btree, return the new btree and a bool sign, true if the
* leafnode is new to the tree, or if we already have the same leafnode with same 
* key, false if the same leafnode with different key is in the tree. usually
* when it's false, something is wrong.*)
val add_to_btree: 'k -> 'a -> ('k,'a) b_tree -> ('k -> unit)
-> ('a -> unit) -> ('a -> 'a -> bool) -> ('k,'a) b_tree * bool

(** [search_in_btree key btree print_f_for_key print_f_for_node] return the
* leafnode with key. print_f_for_key and print_f_for_node are the debug function
* to printout middle-key and leafnode. *)
val search_in_btree: 'k -> ('k,'a) b_tree -> ('k -> unit) -> ('a -> unit) -> 
    'a * ('k,'a) b_tree

(** [remove_from_btree key btree print_f_for_key print_f_for_node] remove the
* node with key, return the new btree. print_f_for_key and print_f_for_node are the debug function to printout middle-key and leafnode.*)
val remove_from_btree: 'k -> ('k,'a) b_tree -> ('k -> unit) -> ('a -> unit) -> 
    ('k,'a) b_tree

(**[iter_b_tree btree f_for_node print_f_for_key debug] iter the binary
* tree, with function f_for_node and pinrt_f_for_key, f_for_node work on each
* node of the tree, don't need to be debug function. but f_for_key debug is the bool option.*)
val iter_b_tree: ('k,'a) b_tree -> ('a -> unit) -> ('k -> unit) -> bool -> unit

(** [print_b_tree btree print_f_for_node print_f_for_key] prints out a btree, it calls iter_b_tree with function print_f_for_node print_f_for_key.*)
val print_b_tree: ('k,'a) b_tree -> ('a -> unit) -> ('k -> unit) -> unit
