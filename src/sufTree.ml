(* POY 5.0 Beta. A phylogenetic analysis program using Dynamic Homologies.    *)
(* Copyright (C) 2013 Andrés Varón, Lin Hong, Nicholas Lucaroni, Ward Wheeler,*)
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

let () = SadmanOutput.register "SufTree" "$Revision: 3160 $"

(** sufTree module implements structures and
* functions of the suffix tree. The implementations
* follow the suffix tree algorithm in Algorithms on Strings, Trees, 
* and Sequences: Computer science and Computational Biology by Dan Gusfield *)

let deref = Utl.deref
open SufNode;;


(**  Rule1:  the character X is added at the end of a leaf
 *   Rule21: create a new internal node, the X is a leaf of the internal node
 *   Rule22: create a new leaf corresponding to the X from an existed internal one
 *   Rule3:  The [...X] has been already in the tree. Do nothing
 *   Rule31: X before the stop_node
 *   Rule32: X end after the stop_node *)
 
type rule = Rule1 | Rule21 | Rule22 | Rule31 | Rule32;;
let max_located_node = 100000
let cur_phase_ptr = ref 0


type sufTree_t = {
    mutable code_arr : int array;
    mutable len : int;

    mutable root : sufNode_t;

    mutable num_located_node : int;
    mutable num_used_node : int;
    mutable located_node_arr : sufNode_t ref array;

    mutable gre_node_id : int;
}


let create_fresh_suf_tree () = {
    code_arr    = [||];
    len = 0;
    root = SufNode.create_fresh_node ();

    num_located_node = 0;
    num_used_node = 0;
    located_node_arr = [||];

    gre_node_id = -1;
}

    
(* Check if the suf tree was constructed correctly *)
let check_tree (tree : sufTree_t) = 
(* travel down from the cur_node for a number of rem_label_length character *)
    for start_pos = 0 to tree.len - 1 do    
        let is_match = SufNode.travel tree.root tree.code_arr 
            start_pos tree.len 
        in  
        if is_match = false then 
        begin
            print_endline (string_of_int start_pos);
            failwith "WRROOOOOOOOOOOOOOOONNNNNNNNNNNNNGGGGGGGGGGGGGGGGGGGGGGG"
        end
    done;       
    print_endline "CORRRRRRRRRRREEEEETLY RECONSTRUCTED"
    
        

(* Create the first node of the suf tree *)
let do_first_phase (tree : sufTree_t) =
    let root = tree.root in 
    root.id <- (tree.len + 1);
    tree.gre_node_id <- (tree.len + 1);
    cur_phase_ptr := 0;
        
    let first_node = !(tree.located_node_arr.(0)) in
    tree.num_used_node <- 1;

    first_node.id <- 0;
    SufNode.set_leaf_label_pos first_node 0 cur_phase_ptr;
    
    first_node.father_ptr <- Some (ref root);
    first_node.suf_link_ptr <- Some (ref root);
        
    root.first_child_ptr <- Some (ref first_node);
    root.suf_link_ptr <- Some (ref root);
    root.father_ptr <- Some (ref root)

    
let locate_mem (tree : sufTree_t) =
    if (tree.num_used_node = tree.num_located_node) then 
    begin
        tree.num_used_node <- 0;
        tree.num_located_node <- (min max_located_node tree.len * 3);
        tree.located_node_arr <- Array.init (tree.num_located_node) 
            (fun _ -> ref (SufNode.create_fresh_node ()) );
    end

        
(* This is the rule22, create a new child (leaf) for the node. 
   This leaf becomes the first child *)
let do_rule22 (tree : sufTree_t) (extension: int) (node : sufNode_t) 
        (ext_code : int)  =   

    let new_leaf = !(tree.located_node_arr.(tree.num_used_node)) in 
    tree.num_used_node <- tree.num_used_node + 1;
    locate_mem tree;

(*  let new_leaf = SufNode.create_fresh_leaf () in *)
    new_leaf.id <- extension; 
    new_leaf.start_pos <- !cur_phase_ptr;  
    new_leaf.end_pos_ptr <- cur_phase_ptr; 
    new_leaf.father_ptr <- Some (ref node); 

    new_leaf.right_sibling_ptr <- node.first_child_ptr;
    node.first_child_ptr <- Some (ref new_leaf)



(* Do rule21, create a new internal node *)
let do_rule21 (tree : sufTree_t) (extension : int) (stop_node : sufNode_t) 
        stop_len : rule * sufNode_t option * sufNode_t option  * sufNode_t = 
        
    let start_label_pos = stop_node.start_pos in 
    let father_stop_node = !(deref stop_node.father_ptr) in         


    let new_inter_node = !(tree.located_node_arr.(tree.num_used_node)) in
    tree.num_used_node <- tree.num_used_node + 1;
    locate_mem tree;    
    

    let new_leaf = !(tree.located_node_arr.(tree.num_used_node)) in
    tree.num_used_node <- tree.num_used_node + 1;
    locate_mem tree;    
    
    new_leaf.id <- extension; 
    new_leaf.start_pos <- !cur_phase_ptr;  
    new_leaf.end_pos_ptr <- cur_phase_ptr; 
    new_leaf.father_ptr <- Some (ref new_inter_node);
                                                 
    tree.gre_node_id  <- tree.gre_node_id + 1;  
    new_inter_node.id  <- tree.gre_node_id;
    new_inter_node.start_pos <-  start_label_pos;
    new_inter_node.end_pos_ptr <- ref (start_label_pos + stop_len - 1);
    new_inter_node.father_ptr <- Some (ref father_stop_node);
    

    (match !(deref father_stop_node.first_child_ptr).id = stop_node.id with 
        | true -> (* the stop_node is the first node *)
            father_stop_node.first_child_ptr <- Some (ref new_inter_node)
        | false ->
            let left_sibling = get_left_sibling 
                !(deref father_stop_node.first_child_ptr) stop_node 
            in 
            left_sibling.right_sibling_ptr <- Some (ref new_inter_node)
    );
    new_inter_node.first_child_ptr <- Some (ref stop_node);
    new_inter_node.right_sibling_ptr <- stop_node.right_sibling_ptr;

    stop_node.right_sibling_ptr <- Some (ref new_leaf);
    stop_node.father_ptr <- Some (ref new_inter_node);
    stop_node.start_pos <- (start_label_pos + stop_len);
    
    
    Rule21, (Some new_inter_node), (Some new_inter_node), stop_node
    


(* create the implicit suf tree for [..phase]. We will add the
    the code_arr[phase] to the tree from the stop_pos
    which is on edge of stop_node *)
let do_explicit_extension (tree : sufTree_t) (extension : int) 
        (stop_node : sufNode_t) (stop_len : int) 
        : rule  * sufNode_t option  * sufNode_t option  * sufNode_t = 
(* The stop_pos is within the edge *)               
    match (SufNode.get_label_len stop_node) > stop_len with  
        | true ->       
            if ( tree.code_arr.(stop_node.start_pos + stop_len) = 
                    tree.code_arr.(!cur_phase_ptr) ) then 
                Rule31, None, None, stop_node
            else do_rule21 tree extension stop_node stop_len
                
        | false -> (*The stop_pos is at the end of the stop_node *)             
            if SufNode.is_leaf stop_node then 
                Rule1, None, None, stop_node
            else begin                
                let child = SufNode.find_child stop_node
                    tree.code_arr.(!cur_phase_ptr) tree.code_arr in 
                match child with
                | None -> do_rule22 tree extension stop_node 
                          tree.code_arr.(!cur_phase_ptr);
                    Rule22, Some stop_node, None, stop_node
                | Some a_child -> Rule32, Some stop_node, None, a_child
            end
    
        
        
(* travel down from the cur_node for a number of rem_label_length character *)
let rec travel (tree : sufTree_t) (node : sufNode_t) (rem_label_len : int) =
    let start_code = tree.code_arr.(!cur_phase_ptr - rem_label_len) in
    let child = SufNode.find_child node start_code tree.code_arr in 
    match child with 
        | Some a_child ->
              let child_label_len = SufNode.get_label_len a_child in 
              if rem_label_len <= child_label_len then 
                  a_child, rem_label_len
              else
                  travel tree a_child (rem_label_len - child_label_len)

        | None-> failwith "Can not travel down the suf tree"

    
(* Create the implicit suf tree T_{phase} *)
let rec do_extension_rec (tree : sufTree_t) (extension : int) 
        (stop_node : sufNode_t) (stop_len : int) 
        (waited_suf_link_node : sufNode_t option) : sufNode_t * int * int =  

    let father_stop_node = !(deref stop_node.father_ptr) in
    let father_suf_link = !(deref father_stop_node.suf_link_ptr) in 
    let ext_rule, inter_node, new_inter, new_stop_node =
        do_explicit_extension tree extension stop_node stop_len in

    match ext_rule with
        | Rule31 -> new_stop_node, (stop_len + 1), (extension - 1)
        | Rule32 -> 
            (match waited_suf_link_node with
                | None -> ()
                | Some a_waited_inter_node -> 
                      a_waited_inter_node.suf_link_ptr <- 
                          Some (ref (deref inter_node))
            );
            new_stop_node, 1, (extension - 1)
        | _ -> 
            let new_stop_node, new_stop_len = 
                if SufNode.is_root father_stop_node then 
                    match stop_len with
                        | 0 | 1 -> tree.root, 0 
                        | _ -> travel tree tree.root (stop_len - 1) 
                else travel tree father_suf_link stop_len 
            in

            (if (ext_rule =  Rule21) || (ext_rule = Rule22) then
                match waited_suf_link_node with
                    | None -> ()
                    | Some a_waited_inter_node -> 
                          a_waited_inter_node.suf_link_ptr <- 
                              Some (ref (deref inter_node))
            );
                
            if extension = !cur_phase_ptr then tree.root, 0, extension
            else do_extension_rec tree (extension + 1) 
                new_stop_node new_stop_len new_inter
            
    
(* Create the suf tree using the code_arr data and Ukkonen's algorithm
   as described by Dan Gusfield *)
let create_suf_tree (code_arr : int array) = 
    
(*    Gc.set { (Gc.get()) with Gc.space_overhead = 500};*)
    let tree = create_fresh_suf_tree () in  
    tree.code_arr <- code_arr;
    tree.len <- Array.length code_arr;  
    locate_mem tree;

(*   first create the implicit suf tree *)  
    do_first_phase tree;    
    let last_extension = ref 0 in
    let stop_node = ref tree.root in 
    let stop_label_len = ref 0 in 
    
    
    for phase = 1 to tree.len - 1 do
        cur_phase_ptr := phase;
        (if (phase mod 1000000 = 0) then 
            print_endline (string_of_int phase) );
        
        let tmp_stop_node, tmp_stop_label_len, tmp_last_extension = 
            do_extension_rec tree (!last_extension + 1) 
                !stop_node !stop_label_len None in
                
        stop_node := tmp_stop_node;
        stop_label_len := tmp_stop_label_len;
        last_extension := tmp_last_extension;        

(*        print_newline ();
          SufNode.print_node tree.root tree.code_arr; *)

    done;

(*   Create the suf_tree from the implicit tree adding a special character 
     which does not appear whereelse in the sequence *)      
    let special_code = max_int in 
    tree.code_arr <- Array.append tree.code_arr [|special_code|];
    tree.len <- tree.len + 1;

    cur_phase_ptr := tree.len - 1;
    let _ = do_extension_rec tree (!last_extension + 1) !stop_node 
        !stop_label_len None in
    
(*    print_endline ("suffix tree len: " ^ (string_of_int tree.len)); *)
    tree
    
(**************************************************************)    


    
let create_k_mer_tree (code_arr : int array) (k : int) = 
    let tree = create_suf_tree code_arr in 
    SufNode.create_depth tree.root tree.root;   
    SufNode.create_leaf_ls tree.root k;
    tree.located_node_arr <- [||];
(*    print_endline "Finish creating k_mer_tree";*)
    tree
        

    
(* find the stop node for this k_mer *) 
let find_stop_node (tree : sufTree_t) (k_mer : int list) =
    let gap_code = Alphabet.get_gap Alphabet.nucleotides in
    let gapless_k_mer = List.map (fun ch -> ch land (gap_code - 1) ) k_mer in 

    let node, rest_k_mer = SufNode.find_stop_node ~exact:false tree.root
        tree.code_arr gapless_k_mer 
    in

    let rest_len = List.length rest_k_mer in 	
    let rest_k_mer, _ = List.fold_right 
        (fun ch (nascent_k_mer, nascent_len) -> 
             if nascent_len = rest_len then nascent_k_mer, nascent_len
             else (ch::nascent_k_mer), (nascent_len + 1) 
        ) k_mer ([], 0)

    in 
    node, rest_k_mer
