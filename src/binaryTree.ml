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
        | Empty -> create_btree key mum
    in
    let res_bt = insert_node old_bt seedseq newmum in
    res_bt,!sign
(*
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
        | Empty -> failwith "search in empty btree"
    in
    search_node bt key
*)

let rotate_child parent_bt leftchild printkey_f =
    let debug = false in
    match parent_bt with
    | Node (middle_key, left_bt, right_bt) ->
            (
                if debug then Printf.printf "rotate child,is_left_child=%b,parent_key=%!" leftchild;
                if debug then printkey_f middle_key;
            if leftchild then
                match left_bt with
                | Node (key, l_bt, r_bt) ->
                    let new_right_bt = Node (middle_key,r_bt,right_bt) in
                    Node(key,l_bt,new_right_bt)
                | Leaf _ -> if debug then Printf.printf "leafnode, do nothing\n"; 
                Node (middle_key, left_bt, right_bt)
            else 
                match right_bt with
                | Node (key, l_bt, r_bt) ->
                    let new_left_bt = Node (middle_key,left_bt,l_bt) in
                    Node(key,new_left_bt,r_bt) 
                | Leaf _ -> if debug then Printf.printf "leafnode, do nothing\n"; 
                Node (middle_key, left_bt, right_bt)
            )
    | _ -> failwith "parent node cannot be a leaf1"

let rotate_twice grand_parent_bt p_g_left(*parent is left child of grandparent*)
me_p_left  (*me is left child of parent*) printkey_f =
    let debug = false in
    match grand_parent_bt with
    | Node (g_key,g_l_bt,g_r_bt) ->
    (
        if debug then Printf.printf "rotate twice, p_g_left=%b,me_p_left=%b,grand parent key="
        p_g_left me_p_left;
        if debug then printkey_f g_key;
        match p_g_left,me_p_left with
        | true,true 
        | false,false ->
            let new_bt = rotate_child grand_parent_bt p_g_left printkey_f in
            (match new_bt with
            | Node (k,l,r) -> rotate_child new_bt me_p_left printkey_f
            | _ -> failwith "new bt cannot be a leaf1")
        | true,false 
        | false,true -> 
            let new_bt =
                match p_g_left with
                | true ->
                        let new_l_bt = rotate_child g_l_bt me_p_left printkey_f in
                        Node (g_key,new_l_bt,g_r_bt) 
                | false ->
                        let new_r_bt = rotate_child g_r_bt me_p_left printkey_f in
                        Node (g_key,g_l_bt,new_r_bt)
            in
            (match new_bt with
            | Node (k,l,r) -> rotate_child new_bt p_g_left printkey_f
            | _ -> failwith "new bt cannot be a leaf2")
                
    )
    | _ -> failwith "grand parent cannot be a leaf2"

(*0: leaf, 1: left, -1:right, 3: root*)
let search_in_btree key bt printkey_f printnode_f = (*mumseq is the key inside b_tree*)
    let debug = false in
    let rec search_node dir sub_bt key =
        match sub_bt with
        | Node (middle_key, left_bt, right_bt) ->
                let search_dir = compare middle_key key in
                let search_dir = 
                    if search_dir=0 then -1 else search_dir in
                if debug then Printf.printf "search node, thisdir=%d,search_dir=%d,thiskey=\n%!"
                dir search_dir;
                if debug then printkey_f middle_key;
                (*update left or right child with search res*)
                let resleaf,current_bt,p_g_left,me_p_left =
                    if search_dir=1 then(*key<middle_key*)
                        let res,new_l_bt,pgl,mepl = 
                            search_node search_dir left_bt key
                        in
                        res,Node (middle_key,new_l_bt,right_bt),
                        pgl,mepl
                    else 
                        let res,new_r_bt,pgl,mepl =
                        search_node search_dir right_bt key
                        in
                        res,Node (middle_key,left_bt,new_r_bt),
                        pgl,mepl
                in
                let new_bt,new_p_g_left,new_me_p_left =
                match p_g_left,me_p_left with
                | 0,0 -> current_bt,0,dir (*back from a leaf child,continue*)
                | 0,y -> current_bt,dir,y (*back from a parent of a leaf,continue*)
                | x,y ->
                        let pgl =
                            if x=1 then true else false
                        and mepl =
                            if y=1 then true else false in
                        rotate_twice current_bt pgl mepl printkey_f,
                        0,dir
                in
                if debug then 
                    Printf.printf "return newbt,%d,%d\n%!" new_p_g_left new_me_p_left;
                resleaf,new_bt,new_p_g_left,new_me_p_left
        | Leaf (leafkey, leafmum) ->
                if (compare leafkey key)=0 then
                    leafmum,Leaf(leafkey,leafmum),0,0
                else begin
                        Printf.printf "we are searching for :%!";
                        printkey_f key;
                        Printf.printf "search failed, we reach\n:%!";
                        printnode_f leafmum;
                    failwith "search reach a dead end in btree";
                end
    in
    let res,newbt,_,_ = search_node 3 bt key in
    res,newbt


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
        | Empty -> failwith "remove from empty btree"
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
        | Empty -> failwith "iter in empty btree"
    in
    iter_sub_tree bt


let print_b_tree bt printnode_f printkey_f =
    Printf.printf "btree is : \n%!";
    iter_b_tree bt printnode_f printkey_f true


