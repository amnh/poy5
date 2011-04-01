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

let () = SadmanOutput.register "BinaryTree" "$Revision: 1751 $"


type ('k,'a) b_tree =
    | Empty
    | Leaf of 'k * 'a 
    | Node of 'k * ('k,'a) b_tree * ('k,'a) b_tree


let create_btree leafkey leafdata =  Leaf (leafkey,leafdata)

let just_a_leaf bt = 
    match bt with 
    | Leaf (leafkey,leafmum) -> true
    | _ -> false

let just_an_empty_leaf bt =
    match bt with 
    | Empty -> true
    | _ -> false

(*add new mum_node to b_tree, if same mum already exist, replace it with new one
* if other mum with same match sequence already exist, do nothing, return false*)
let add_to_btree seedseq newmum old_bt printkey_f printnode_f compare_node_f =
    let debug = false in
    if debug then 
        begin Printf.printf "add to btree\n%!"; printnode_f newmum; end;
    let sign = ref true in
    let rec insert_node sub_bt key mum =
        match sub_bt with
        | Node (middle_key, left_bt, right_bt) ->
                let insert_dir = compare middle_key key in
                if insert_dir=1 then (*middle_key>key*)
                    let new_left_bt = insert_node left_bt key mum in
                    Node (middle_key, new_left_bt, right_bt)
                else 
                    let new_right_bt = insert_node right_bt key mum in
                    Node (middle_key,
                    left_bt,
                    new_right_bt)
        | Leaf (leafkey, leafmum) ->
                let insert_dir = compare key leafkey in
                if debug then begin
                    printkey_f key;
                    printkey_f leafkey;
                    Printf.printf " insert_dir=%d,%!" 
                    insert_dir;
                    printnode_f leafmum;
                end;
                if insert_dir=1 then (*key>leafkey*)
                    Node (key, 
                          Leaf (leafkey, leafmum),
                          Leaf (key,mum))
                else if insert_dir=(-1) then
                    Node (leafkey, 
                          Leaf (key,mum),
                          Leaf (leafkey, leafmum))
                else 
                    if (compare_node_f newmum leafmum) then  
                        Leaf (leafkey,newmum)
                    else begin
                        sign := false;
                        Leaf (leafkey,leafmum)
                    end
    in
    let res_bt = insert_node old_bt seedseq newmum in
    res_bt,!sign

let search_in_btree key bt printkey_f printnode_f = 
    let rec search_node sub_bt key =
        match sub_bt with
        | Node (middle_key, left_bt, right_bt) ->
                let search_dir = compare middle_key key in
                if search_dir=1 then(*key<middle_key*)
                    search_node left_bt key
                else 
                    search_node right_bt key
        | Leaf (leafkey, leafmum) ->
                if (compare leafkey key)=0 then leafmum
                else begin
                        Printf.printf "we are searching for :%!";
                        printkey_f key;
                        Printf.printf "search failed, we reach\n:%!";
                        printnode_f leafmum;
                    failwith "search reach a dead end in btree";
                end
    in
    search_node bt key

let remove_from_btree seedseq bt printkey_f printnode_f = 
    let rec remove_node sub_bt key = 
        match sub_bt with 
        | Node (middle_key, left_bt, right_bt) ->
                let remove_dir = compare middle_key key in
                if remove_dir=1 then (*key<middle_key*)
                    let new_left_child = remove_node left_bt key in
                    (match new_left_child with
                    | Empty ->  right_bt
                    | x ->  
                            Node (middle_key,x,right_bt))
                else 
                    let new_right_child = remove_node right_bt key in
                    (match new_right_child with
                    | Empty -> left_bt
                    | x -> 
                            Node (middle_key,left_bt, x))
        | Leaf (leafkey, leafmum) ->
                if (compare leafkey key)=0 then Empty
                else begin
                    Printf.printf "we are searching for :%!";
                        printkey_f key;
                        Printf.printf "search failed, we reach\n:%!";
                        printnode_f leafmum;
                    failwith "remove node reach a dead end in btree"
                end
    in
    remove_node bt seedseq 
    
let iter_b_tree bt node_f printkey_f debug =
    let rec iter_sub_tree sub_bt =
        match sub_bt with
        | Node (middle_key, left_bt, right_bt) ->
                if debug then begin
                    Printf.printf "{ mkey=%!";
                    printkey_f middle_key;
                    Printf.printf "\n ( %!";
                end;
                iter_sub_tree left_bt;
                iter_sub_tree right_bt;
                if debug then Printf.printf " ) }\n%!";
        | Leaf (leafkey, leafmum) -> 
                node_f leafmum 
    in
    iter_sub_tree bt


let print_b_tree bt printnode_f printkey_f =
    Printf.printf "btree is : \n%!";
    iter_b_tree bt printnode_f printkey_f true


