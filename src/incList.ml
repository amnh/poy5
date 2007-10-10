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

let () = SadmanOutput.register "IncList" "$Revision: 736 $"
(* increasing list *)
type 'a incList_t = { 
    mutable data : 'a;
    mutable next : 'a incList_t ref option;
}

let len (incList : 'a incList_t) = 
    let rec cmp_len (cur_incList : 'a incList_t) (num_cell: int) = 
        match cur_incList.next with
            | None -> num_cell 
            | Some tail -> cmp_len !tail (num_cell + 1)
    in
    
    cmp_len incList 0
    
    
let rec iter (incList : 'a incList_t) (app_fun : 'a -> unit) = 
    app_fun incList.data;
    match incList.next with
        | None -> ()
        | Some tail -> iter !tail app_fun
    

let rec insert (incList : 'a incList_t) (cmp_fun : 'a -> 'a  -> int) 
        (new_cell : 'a incList_t) =   

    match incList.next with
        | None -> 
            incList.next <- Some (ref new_cell);
            new_cell.next <- None;
        | Some next_cell ->         
            if cmp_fun !next_cell.data new_cell.data >= 0 then 
            begin
                incList.next <- Some (ref new_cell);
                new_cell.next <- Some next_cell
            end 
            else insert !next_cell cmp_fun new_cell
            
            
let rec check_increasing (incList : 'a incList_t) (cmp_fun : 'a -> 'a  -> int) =
    match incList.next with
        | None -> true
        | Some next_cell ->
            if cmp_fun incList.data !next_cell.data > 0 then false
            else check_increasing !next_cell cmp_fun
