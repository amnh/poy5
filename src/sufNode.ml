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

let () = SadmanOutput.register "SufNode" "$Revision: 1644 $"
(** The implementation of suffix node in suffix tree *)

let fprintf = Printf.fprintf
let deref = Utl.deref


type sufNode_t = {
    mutable id                  : int;
    mutable father_ptr          : sufNode_t ref option;
    mutable suf_link_ptr        : sufNode_t ref option;

    mutable first_child_ptr     : sufNode_t ref option;
    mutable right_sibling_ptr   : sufNode_t ref option;
    mutable depth               : int;
    
    mutable leaf_ls             : int list;
    mutable start_pos           : int;
    mutable end_pos_ptr         : int ref
}


let create_fresh_node () = {
    id = -1; 
    father_ptr = None; 
    suf_link_ptr = None;
    first_child_ptr = None;
    right_sibling_ptr = None; 

    depth = 0;  

    leaf_ls = [];
    start_pos = -1;
(* to make sure that the label_len = end_pos - start_pos + 1 = 0 *)
    end_pos_ptr = ref (-2)  
}

let get_label_len (node : sufNode_t) = !(node.end_pos_ptr) - node.start_pos + 1

let set_leaf_label_pos (node : sufNode_t) (new_start_pos : int) 
        (new_end_pos_ptr : int ref) =
    node.start_pos <- new_start_pos;
    node.end_pos_ptr <- new_end_pos_ptr
        
let set_inter_label_pos (node : sufNode_t) (new_start_pos : int) 
        (new_end_pos : int) =
    node.start_pos <- new_start_pos;
    node.end_pos_ptr := new_end_pos

let is_root (node : sufNode_t) = 
    node.id = !(deref node.father_ptr).id
            
let is_leaf (node : sufNode_t)  : bool = node.first_child_ptr = None                        

    
let rec create_depth (father_node : sufNode_t) (node : sufNode_t) = 
    node.depth <- father_node.depth + !(node.end_pos_ptr) - node.start_pos + 1;
    let rec check_child (child : sufNode_t ref option) =    
        match child with
            | None -> ()
            | Some a_child_ptr -> 
                create_depth node !a_child_ptr;
                check_child !a_child_ptr.right_sibling_ptr; 
    in
    check_child node.first_child_ptr


let cmp_num_child (node : sufNode_t) = 
    let rec check_child (child : sufNode_t ref option) (num_child : int) =  
        match child with
        | None -> num_child
        | Some a_child_ptr -> 
              check_child !a_child_ptr.right_sibling_ptr (num_child + 1);   
    in
    check_child node.first_child_ptr 0

let rec create_leaf_ls (node : sufNode_t) (min_depth : int) =           
    if is_leaf node then 
        node.leaf_ls <- [node.id]
    else begin
        let rec check_child (child : sufNode_t ref option) =
            match child with
                | None -> ()
                | Some a_child_ptr -> 
                    create_leaf_ls !a_child_ptr min_depth;
                    (if node.depth >= min_depth then 
                        node.leaf_ls <- List.merge compare 
                            node.leaf_ls !a_child_ptr.leaf_ls);
                    check_child !a_child_ptr.right_sibling_ptr
    
        in

        check_child node.first_child_ptr;
        
        if node.depth >= min_depth then 
            node.first_child_ptr <- None
        else ()
    end     



let find_child ?(exact=true) (node : sufNode_t) looking_code code_arr =
    let rec check_child (child : sufNode_t ref option) =
        match child with
        | None -> None 
        | Some a_child_ptr ->  
              let code = code_arr.(!a_child_ptr.start_pos) in 
              match (looking_code = code)  || ( (looking_code land code > 0) && (exact = false) ) with
              | true -> Some !a_child_ptr    
              | false -> check_child !a_child_ptr.right_sibling_ptr  

    
    in
    check_child node.first_child_ptr    

let rec get_left_sibling (node : sufNode_t) (looking_node  : sufNode_t) = 
    let right_sibling = !(deref node.right_sibling_ptr) in
    match right_sibling.id = looking_node.id with
        | true -> node
        | false -> get_left_sibling right_sibling looking_node
    


let rec travel (node : sufNode_t) (code_arr : int array) (cur_pos : int)
        (stop_pos : int) =
    if cur_pos = stop_pos then true
    else begin
        match find_child node code_arr.(cur_pos) code_arr with
        | None ->  failwith "Can not go down to any child"  
        | Some a_child ->  
              let rec check offset =  
                  match code_arr.(cur_pos + offset) =
                  code_arr.(a_child.start_pos + offset) with  
                  | false -> false  
                  | true ->            
                        if offset = !(a_child.end_pos_ptr) - a_child.start_pos then true  
                        else check (offset + 1)  
              in  
                
              match check 0 with 
              | false -> false     
              | true -> travel  a_child code_arr (cur_pos +
                  !(a_child.end_pos_ptr) - a_child.start_pos + 1)stop_pos 


    end
    
    
(* print the implicit suf tree *)
let rec print_node  (node : sufNode_t) (code_arr : int array) = 
    if is_leaf node = false then 
    begin
        let suf_link_id = 
            match node.suf_link_ptr with 
                | None -> -1
                | Some a_link_ptr -> !a_link_ptr.id
        in

        fprintf stdout "Id: %i, depth: %i, suf_link: %i -> \n" 
            node.id node.depth suf_link_id;

        let rec check_child (child : sufNode_t ref option) =
            match child with
            | None -> print_newline ()
            | Some a_child_ptr -> 
                check_child !a_child_ptr.right_sibling_ptr;
                print_node !a_child_ptr code_arr
        in
        check_child node.first_child_ptr            
    end
(********************************************************************)          

        
let check_node ?(exact=true) (node : sufNode_t) (code_arr : int array) (cur_mer : int list) =
    let edge_len = !(node.end_pos_ptr) - node.start_pos + 1 in 
    let rec check offset cur_mer =
        if offset = edge_len then 
            true, cur_mer 
        else
            match cur_mer with
            | [] -> true, [] 
            | head::tail ->  
                  let code = code_arr.(node.start_pos + offset) in
                  match (code = head) || ( (code land head > 0) && (exact = false) ) with  
                  | false -> false, cur_mer 
                  | true -> check (offset + 1) tail           
    in
    check 0 cur_mer
    

(* find the stop node for this k_mer *) 
let rec find_stop_node ?(exact=true) (node : sufNode_t) (code_arr : int array) 
        (cur_mer : int list) : sufNode_t * int list =

    if is_leaf node then node, cur_mer
    else
        match cur_mer with
        | [] -> node, []  
        | head_code :: _ -> 
              match find_child ~exact:exact node head_code code_arr with  
              | None -> node, cur_mer  
              | Some a_child ->   
                    let is_identical, res_mer =   
                        check_node ~exact:exact a_child code_arr cur_mer in   
                    match is_identical with  
                    | false -> node, cur_mer  
                    | true -> find_stop_node ~exact:exact a_child code_arr res_mer 


(* find the stop node for this k_mer *) 
let rec find_leaves ?(exact=true) node code_arr (cur_mer : int list) =
    match cur_mer with 
    | [] -> node.leaf_ls 
    | code :: _ -> begin 
          let leaf_ls = ref [] in  
          let rec check_child (child : sufNode_t ref option) = 
              match child with 
              | None -> ()
              | Some chi ->   
                    let chi_code = code_arr.(!chi.start_pos) in  
                    if (code = chi_code)  || ((code land chi_code > 0) && (exact = false)) then begin                            
                        let is_identical, res_mer =  
                            check_node ~exact:exact !chi code_arr cur_mer  
                        in    
                            
                        if is_identical then begin 
                            let chi_leaf_ls = find_leaves ~exact:exact !chi code_arr res_mer in
                            leaf_ls := chi_leaf_ls @ !leaf_ls
                        end;   
                    end;                               
                    check_child !chi.right_sibling_ptr      
          in 
          check_child node.first_child_ptr;     
          !leaf_ls 
      end  

