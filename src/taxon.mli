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


module Taxon : sig
    (* Basic structure for the taxon functionality *)
    type t = 
        {
            code: int;  
            (* Internal unique code used in all tree operations. *)
            name: string;  
            (* Name of the taxon as provided by the user. *)
            characters: Character2.characters array ref;
            (* Reference to the characters of the taxon. *)
        };;

    (* [create_characters t c] Creates the character arrays required for
    * the taxon to work in a tree operation. *)
    val create_characters : t ref -> int -> 
        int list array ref * Character2.characters array ref array ref * 
        int array array ref * Perturb.p array ref * int ref * 
        int ref * int list ref

    (* [get_name t] returns the name of taxon referenced by [t] *)
    val get_name : t ref -> string

    (* [get_code t] returns the internal code assigned to taxon referenced by 
    * [t] *)
    val get_code : t ref -> int
    
    val get_code_no_ref : t -> int

    (* [get_code_str t] returns the string version of the numberic code of the 
    * taxon referenced by [t] *)
    val get_code_str : t ref -> string

    (* [get_characters t] returns a reference to the array of characters 
    * assigned to the taxon referenced by [t] *)
    val get_characters : t ref -> Character2.characters array ref

    (* [empty_taxon ()] returns a reference to a fresh type taxon. *)
    val empty_taxon : unit -> t ref

    (* [new_htu ()] is a synonim of [empty_taxon ()] *)
    val new_htu : unit -> t ref

    (* Produces an array with all the otu's index values *)
 (*   val all_index_array : unit -> int array *)

    val output : out_channel -> t -> unit 

end


module Taxa : sig
    (* A set of taxons *)

    type s 
    (* A set of taxa type *)

    val create : unit -> s
    (* [create ()] creates a new set of taxa *)

    val add : Taxon.t -> s -> unit
    (* [add x y] adds the taxon x to the taxa set y *)

    val remove : Taxon.t -> s -> unit
    (* Removes a taxon from the taxa set *)

    val find : 'a -> ('a, exn) Hashtbl.t -> exn
    (*val find : int -> s -> Taxon.t *)
    (* [find x y] returns the taxon with code x in the set y. If the taxon is
    * not found raises a Not_found exception. *)

    val to_list : s -> Taxon.t list 
    (* [to_list x] converts a taxa set in a list of taxons *)

end

