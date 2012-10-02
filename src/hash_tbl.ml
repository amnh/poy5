(* POY 5.0 Alpha. A phylogenetic analysis program using Dynamic Homologies.   *)
(* Copyright (C) 2011 Andrés Varón, Lin Hong, Nicholas Lucaroni, Ward Wheeler,*)
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

let () = SadmanOutput.register "Hash_tbl" "$Revision: 2684 $"

(** The hash table record that contains the id generated by h2, the
 * number of times the partition occurs so far and the number of elements
 * in the partition. *)
type hash_tbl_record = { id : int ; 
                         (** identity generated by h2 *)
                         count : int ;
                         (** number of times this bi-partition occurs. *)
                         size : int ;
                         (** number of elements in the bi-partition. *) 
                         index : int
                         (** if the split was of size 1, then
                            this field holds the index of the leaf node.
                            since all legal indexes start at 0, this field
                            needs to be initialized to -1. *) }
                         
(** The array to hold the list of records. *)
type hash_tbl = { records : hash_tbl_record list array ;
                  (** the records of the hash table. *)
                  h1_params : int array ;
                  (** the a_i's for h1. *)
                  h2_params : int array ;
                  (** the a_i's for h2. *) 
                  m1 : int ;
                  (** prime number >= num_trees * num_nodes_per_tree *)
                  m2 : int 
                  (** prime number atleast twice as large as m1 *) }

(** An empty record. *)
let empty_record = { id = 0 ; 
                    count = 0 ; 
                    size = 0 ; 
                    index = -1 } 

(** [make h1_p h2_p m1 m2]
    @param h1_p the parameters for the universal hash function h1.
    @param h2_p the parameters for the universal hash function h2.
    @param m1 a prime number greater than the number of trees *
        number_of_nodes_per_tree.
    @param m2 a prime number greater than m1. 
    @return a newly created hash table with n empty records and the given
        parameters. *)
let make h1_p h2_p m1 m2 = 
    (* make sure the parameters are consistent. *)
    assert((m1 > 0) && (m2 > m1)) ;
    assert((Array.length h1_p) = (Array.length h2_p)) ;
    for i = 0 to (Array.length h1_p) - 1 do
        assert((h1_p.(i) >= 0) && (h1_p.(i) < m1)) ;
        assert((h2_p.(i) >= 0) && (h2_p.(i) < m2)) ;
    done ;
    (* Make a hash table with empty records. *)
    { records = Array.make m1 [] ;
      h1_params = h1_p ;
      h2_params = h2_p ;
      m1 = m1 ;
      m2 = m2 }

(** {1 Interface} *)

(** [Interface] provides a more reasonable interface to the functions in this
    module. *)
module Interface = struct
    (** [fp] is an individual fingerprint *)
    type fp = (int * int)
            
    (** [t] is an object that responds to queries for the hash values on a tree
    *)
    type t = (int, int * int * int) Hashtbl.t * hash_tbl

    (** [def_safety] is the default safety level (see paper, above) *)
    let def_safety = 500

    (** [make_params ?safety_factor ntrees nodes_per_tree] makes an
        initialization parameter for a certain number of trees and nodes per
        tree.  These parameters can be used in [make_wparams], below.  The
        advantage (compared to using [make]) is being able to make several hash
        tables with the same parameters, so that they end up compatible.
        Otherwise, fingerprints mean different things in different hash
        tables. *)
    let make_params ?(safety_factor = def_safety) ntrees nodes_per_tree =
        (* m1 is prime > ntrees * nodes_per_tree *)
        let m1 = Primes.Probable.next_prime_gt (ntrees * nodes_per_tree) in
        (* m2 is prime > m1 *)
        let m2 = Primes.Probable.next_prime_gt (m1 * safety_factor) in

        let num_leaves = (nodes_per_tree + 1) / 2 in

        let h1_p = Array.make num_leaves 0 in
        let h2_p = Array.make num_leaves 0 in
        for i = 0 to num_leaves - 1 do
            h1_p.(i) <- (Random.int m1) ;
            h2_p.(i) <- (Random.int m2)
        done;

        (m1, m2, h1_p, h2_p, nodes_per_tree)

    (** [make_wparams (m1, m2, h1_p, h2_p, nodes_per_tree)] makes a hash table
        using the parameters made above. *)
    let make_wparams (m1, m2, h1_p, h2_p, nodes_per_tree) =
        let hash_tbl = make h1_p h2_p m1 m2 in
        let ht = Hashtbl.create nodes_per_tree in
        ((ht, hash_tbl) : t)

    (** [make] is a combination of [make_params] and [make_wparams], above. *)
    let make ?(safety_factor = def_safety) ntrees nodes_per_tree =
        (* m1 is prime > ntrees * nodes_per_tree *)
        let m1 = Primes.Probable.next_prime_gt (ntrees * nodes_per_tree) in
        (* m2 is prime > m1 *)
        let m2 = Primes.Probable.next_prime_gt (m1 * safety_factor) in

        let num_leaves = (nodes_per_tree + 1) / 2 in

        let h1_p = Array.make num_leaves 0 in
        let h2_p = Array.make num_leaves 0 in
        for i = 0 to num_leaves - 1 do
            h1_p.(i) <- (Random.int m1) ;
            h2_p.(i) <- (Random.int m2)
        done;

        let hash_tbl = make h1_p h2_p m1 m2 in
        let ht = Hashtbl.create nodes_per_tree in
        ((ht, hash_tbl) : t)

    (** [query t nodeid] returns the hash of the clade below [nodeid] from [t] *)        
    let query (t : t) nodeid =
        let (node_hashes, hash_tbl) = t in

        let (h1, h2, _) = Hashtbl.find node_hashes nodeid in
        ((h1, h2) : fp)
end
