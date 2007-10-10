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

let () = SadmanOutput.register "Taxon" "$Revision: 1644 $"


module Taxon = struct

    type t = 
        {
            code: int; 
            name: string; 
            characters: Character2.characters array ref
        };;

    (*  Produces a ref to a list containing all integers from it to s excepting 
    *  c.*)
    let rec create_index_avail c s it l =
        if not (it = s) then begin
                if not (c = it) then create_index_avail c s (it + 1) (it :: l)
                else create_index_avail c s (it + 1) l
        end else ref l;;

    (*  Produce an array ref of size s with all integers from 0 to s *)
    let create_active_characters s =
        let tmp = Array.create s (Perturb.create_perturb ())  
        and n = s - 1 in
        for i = 0 to n do
            Perturb.set_position (tmp.(i)) i;
        done;
        ref tmp;;

    (*  Appends the element it to the list contained in the item t of  the array
    *  referenced by arr. *)
    let append_to_type_list arr it t =
        !arr.(t) <- it :: !arr.(t);;

    (* Creates a list of indexes in the character array c where each item in the
    * array corresponds to the respective get_type_index of the type. The 
    * produced array index the different types found in the array of 
    * character c for faster access to the different types. *)
    let create_types_index s c =
        let n = s - 1 
        and tmp = Character2.create_types_array [] in
        for i = 0 to n do
            let t = Character2.get_types !c.(i) in
            let t = List.map (Character2.get_type_index) t in
            List.iter (append_to_type_list tmp i) t;
        done;
        tmp;;

    let create_characters t (s : int) =
        let size = Array.length !(!t.characters) in
        let types_index = create_types_index size !t.characters
        and characters = ref (Array.make s !t.characters)
        and costs = ref (Array.make_matrix s s 0)
        and active_characters = create_active_characters size 
        and number_active = ref size
        and number_characters = ref size
        and index_avail = create_index_avail !t.code size 0 [] in
        types_index, characters, costs, active_characters, number_active,
        number_characters, index_avail;;
        
    let make c n =
        {code = c; name = n; characters = Character2.initial ()};;

    let get_name t =
        !t.name;;

    let get_code t =
        !t.code;;

    let get_code_no_ref t =
        t.code;;

    let get_code_str t =
        Pervasives.string_of_int (get_code t);;

    let get_characters t =
        !t.characters;;

    let empty_taxon () =
        ref {code = 0; name = ""; characters = ref [||]};;

    let new_htu () = empty_taxon ();;

    let output o t =
        ();;

end


module Taxa = struct

    type s = (int, Taxon.t) Hashtbl.t
    
    let create () = Hashtbl.create 100 

    let add x y = 
        if (Hashtbl.mem y x.Taxon.code) then () 
        else Hashtbl.add y x.Taxon.code x

    let remove x y = 
        if (Hashtbl.mem y x.Taxon.code) then Hashtbl.remove y
        x.Taxon.code 
        else ()

    let find x y =
        if (Hashtbl.mem y x) then Hashtbl.find y x 
        else Not_found

(*    let to_list x =
        Hashtbl.fold (fun x y z -> z :: x) [] x *)
    let aux_to_list x y z =
        y :: z

    let to_list x =
        let res = ref [] in
        Hashtbl.fold aux_to_list x !res

    let of_list x =
        let res = create () in
        List.iter (fun x -> add x res) x;
        res

end

