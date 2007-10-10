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

let () = SadmanOutput.register "SparceMatrix" "$Revision: 1616 $"

(* SparceMatrix provides a simple matrix for which some
columns and rows can be added and removed every now and then. Instead of
recreating the whole matrix, we simple update a map containing the information
(remove the row and remove the column) with a log n complexity type, where n is
the size of the longest side matrix. Retrieving and adding an element is also
log n *)
type m = {
    items : All_sets.Integers.t;
    matrix : (float Lazy.t All_sets.IntegerMap.t) All_sets.IntegerMap.t;
}

let empty = 
    {
        items = All_sets.Integers.empty;
        matrix = All_sets.IntegerMap.empty;
    }

let unity matrix = 
    All_sets.IntegerMap.fold 
    (fun code map acc ->
        acc && All_sets.Integers.mem code matrix.items && 
        (All_sets.IntegerMap.fold 
        (fun code _ acc -> acc && All_sets.Integers.mem code matrix.items)
        map
        true)) matrix.matrix true

let find x y m =
    let x = min x y 
    and y = max x y in
    let res = 
        try All_sets.IntegerMap.find y (All_sets.IntegerMap.find x m.matrix)
        with
        | Not_found ->
            All_sets.IntegerMap.find x (All_sets.IntegerMap.find y m.matrix)
    in
    Lazy.force_val res


let minimum v m =
    let res = 
        All_sets.Integers.fold (fun x acc ->
            if x = v then acc
            else 
                match acc with
                | None -> Some (find x v m)
                | Some c ->
                        let tmp = find x v m in
                        if tmp < c then Some tmp
                        else acc)
        m.items
        None
    in
    match res with
    | None -> failwith "Empty matrix?"
    | Some v -> v

let row x m = 
    All_sets.IntegerMap.map
        Lazy.force_val 
        (All_sets.IntegerMap.find x m.matrix)

let column x m =
    All_sets.IntegerMap.fold 
    (fun y set acc ->
        if All_sets.IntegerMap.mem x set then
            All_sets.IntegerMap.add y 
            (Lazy.force_val (All_sets.IntegerMap.find x set))
            acc
        else acc)
    m.matrix
    All_sets.IntegerMap.empty

let add x f m =
    let items = All_sets.Integers.add x m.items in
    let row, matrix = 
        All_sets.Integers.fold 
        (fun y (row, matrix) ->
            let dir1 = Lazy.lazy_from_fun (fun () -> f x y) in
            (All_sets.IntegerMap.add y dir1 row),
            (
                All_sets.IntegerMap.add x
                (All_sets.IntegerMap.add y dir1
                (try All_sets.IntegerMap.find y matrix with
                | Not_found -> All_sets.IntegerMap.empty))
                matrix
            )
        )
        items
        (All_sets.IntegerMap.empty, m.matrix)
    in
    { 
        items = items;
        matrix = All_sets.IntegerMap.add x row matrix;
    }

let remove x m =
    assert (unity m);
    let items = All_sets.Integers.remove x m.items 
    and matrix = All_sets.IntegerMap.remove x m.matrix in
    let matrix =
        All_sets.IntegerMap.fold 
        (fun y z acc ->
            let newr = All_sets.IntegerMap.remove x z in
            All_sets.IntegerMap.add y newr acc)
        matrix
        matrix
    in
    let res = { items = items; matrix = matrix } in
    assert (unity res);
    res

let mem code m = All_sets.Integers.mem code m.items

let rows m = m.items

let columns = rows
