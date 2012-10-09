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

let () = SadmanOutput.register "AllDirNode" "$Revision: 2717 $"

let eager        = false
let uppass_debug = false

type exclude = Node.exclude

type direction = (int * int) option

type 'a my_lazy = 
    | Lazy of 'a Lazy.t
    | Eager of 'a

type a_node = Node.Standard.n my_lazy

let error_user_message format = Printf.ksprintf (Status.user_message Status.Error) format
let info_user_message format = Printf.ksprintf (Status.user_message Status.Information) format
let failwithf format = Printf.ksprintf failwith format

let pp_list to_str chan v = List.iter (fun x -> output_string chan (to_str x)) v
let pp_opt_list chan v =
        pp_list (function | Some x -> string_of_float x | None -> "none")
                chan v

(* to print lists through Printf *)
let pp_list depth to_string chan lst = 
	output_string chan "[ ";
	List.iter
		(fun x -> 
			output_string chan (to_string x);
			output_char chan ' ';)
		lst;
	output_char chan ']'

let pp_opt_lst depth to_string = 
    pp_list depth (fun x -> match x with | Some x -> to_string x | None -> "None")
(* ------ *)

let force_val = function
    | Lazy x -> Lazy.force_val x
    | Eager x -> x

let force_opt = function
    | Some x -> Some (force_val x)
    | None -> None 

let lazy_from_fun x =
    if eager then
        let l = Lazy (Lazy.lazy_from_fun x) in
        Eager (force_val l)
    else
        Lazy (Lazy.lazy_from_fun x)
    
let lazy_from_val x = Eager x

let force x = Eager (force_val x)

type node_dir = {
    lazy_node : a_node;
    dir : direction;
    code : int;
}

type node_data = {
    unadjusted : node_dir list; (** The standard downpass node *)
    adjusted : node_dir option;  (** adjusted data calculated in uppass *)
}

let get_adjusted_nodedata data msg =
    match data.adjusted with
    | None -> failwith msg
    | Some x -> x

let print_node_data ndata print_unadjusted=
    Printf.printf "allDirNode, node_data.adjusted:\n%!";
    let () = match ndata.adjusted with
        | Some nodedir ->
            let dir = nodedir.dir in
            let () = match dir with 
                | Some (x,y) -> Printf.printf "dir=(%d,%d),%!" x y
                | None -> Printf.printf "no dir,%!"
            in
            let anode = nodedir.lazy_node in
            Node.print (force_val anode);
        | None -> Printf.printf "no adjusted data\n%!"
    in
    if print_unadjusted then begin
        Printf.printf "\n node_data.unadjusted:\n%!";
        List.iter 
            (fun nodedir ->
                let dir = nodedir.dir in
                let () = match dir with 
                    | Some (x,y) -> Printf.printf "dir=(%d,%d),%!" x y
                    | None -> Printf.printf "no dir,%!"
                in
                let anode = nodedir.lazy_node in
                Node.print (force_val anode);)
            ndata.unadjusted;
        print_newline();
    end

let to_n node = Eager node

let has_code code n =
    match n.dir with
    | None -> true
    | Some (a, b) ->
            a <> code && b <> code

let yes_with code n =
    assert (0 <> List.length n);
    match List.filter (fun x -> not (has_code code x)) n with
    | [x; y] -> x, y
    | x ->  raise Not_found

let not_with code n =
    assert (0 <> List.length n);
    match List.filter (has_code code) n with
    | [x] -> x
    | x   -> raise Not_found

(* grabs a node in the direction with children c1 and c2 --Internal Nodes only *)
let with_both c1 c2 n = 
    match List.filter (fun x -> match x.dir with
                        | Some (xa,xb) -> (xa=c1 && xb=c2) || (xa=c2 && xb=c1)
                        | None -> false ) n with
    | [x] -> x
    |  x  -> raise Not_found

let either_with c n = (* grabs a node that has this branch c--n, in n *)
    assert( 0 <> List.length n);
    match List.filter (fun x -> not (has_code c x)) n with
        | hd::tl -> hd
        | [] -> (match n with
                 | [x] -> x
                 |  _  -> raise Not_found)

module OneDirF : 
    NodeSig.S with type e = exclude with type n = a_node with type other_n =
        Node.Standard.n with type nad8 = Node.Standard.nad8 = struct

    type n = a_node
    type other_n = Node.Standard.n

    let to_other x = force_val x

    type e = exclude

    let recode f n = 
        lazy_from_fun (fun () -> Node.Standard.recode f (force_val n))

    let compare a b = 
        Node.Standard.compare (force_val a) (force_val b)

    let load_data ?(is_fixedstates=false) ?(silent=true) ?(classify=true) data = 
        let data, nodes = Node.Standard.load_data ~is_fixedstates ~classify data in
        data, List.map to_n nodes

    let fix_preliminary x = x

    let apply_f_on_lazy f a b = 
        let a = force_val a
        and b = force_val b in
        f a b

    let distance ?(para=None) ?(parb=None) missing_distance a b =
        apply_f_on_lazy (Node.Standard.distance missing_distance) a b

    let character_costs x n = 
        Node.Standard.character_costs x (force_val n)

    let get_nonadd_8 x n =
        Node.Standard.get_nonadd_8 x (force_val n)

    let get_nonadd_16 x n =
        Node.Standard.get_nonadd_16 x (force_val n)

    let get_nonadd_32 x n =
        Node.Standard.get_nonadd_32 x (force_val n)

    let get_addgen x n =
        Node.Standard.get_addgen x (force_val n)

    let get_addvec x n =
        Node.Standard.get_addvec x (force_val n)

    let get_sank x n =
        Node.Standard.get_sank x (force_val n)

    let get_fixedstates x n =
        Node.Standard.get_fixedstates x (force_val n)

    let get_dynamic x n =
        Node.Standard.get_dynamic x (force_val n)

    let get_mlstatic x n =
        Node.Standard.get_mlstatic x (force_val n)

    let median ?branches my_code old a b = 
        lazy_from_fun 
        (fun () ->
            apply_f_on_lazy
                (Node.Standard.median ?branches my_code None)
                a b)

    let edge_iterator par_opt mine ch1 ch2 = 
        let n_mine = Node.Standard.edge_iterator (force_opt par_opt)
            (force_val mine) (force_val ch1) (force_val ch2)
        in (to_n n_mine)

    let readjust mode set_2adjust a b c mine =
        let node,set = Node.Standard.readjust mode set_2adjust
                            (force_val a) (force_val b) (force_val c) (force_val mine)
        in
        (to_n node,set)

    let final_states x par cur a b = 
        lazy_from_fun 
        (fun () -> apply_f_on_lazy
            (Node.Standard.final_states x (force_val par) (force_val cur)) a b)

    let apply_time root x y =
        lazy_from_fun 
            (fun () -> 
                let x = force_val x in
                let y = force_val y in
                Node.Standard.apply_time root x y)
            
    let get_times_between a b = match b with
        | Some b -> Node.Standard.get_times_between (force_val a) (Some (force_val b))
        | None   -> Node.Standard.get_times_between (force_val a) None

    let uppass_heuristic pcode ptime mine a b = mine

    let extract_states alph d a c0 b = Node.Standard.extract_states alph d None c0 (force_val b)

    let to_string v = Node.Standard.to_string (force_val v)

    let apply_single_f_on_lazy f a = f (force_val a)

    let total_cost x n =
        apply_single_f_on_lazy (Node.Standard.total_cost x) n

    let min_prior c n =
        apply_single_f_on_lazy (Node.Standard.min_prior c) n

    let node_cost par n = 
        apply_single_f_on_lazy (Node.Standard.node_cost par) n

    let update_leaf x = x

    let taxon_code n = 
        apply_single_f_on_lazy Node.Standard.taxon_code n

    let union_distance _ _ = 0.0

    let is_collapsable clas a b = 
        apply_f_on_lazy (Node.Standard.is_collapsable clas) a b

    let classify_data l1 n1 l2 n2 chars acc =
        Node.Standard.classify_data l1 (force_val n1) l2 (force_val n2) chars acc

    let to_xml _ _ _ = ()

    let num_height c n = 
        apply_single_f_on_lazy (Node.Standard.num_height c) n

    let num_otus c n =
        apply_single_f_on_lazy (Node.Standard.num_otus c) n

    let get_sequences c n = 
        apply_single_f_on_lazy (Node.Standard.get_sequences c) n

    let get_dynamic_preliminary code n = 
        apply_single_f_on_lazy (Node.Standard.get_dynamic_preliminary code) n

    let get_dynamic_adjusted code n =
        apply_single_f_on_lazy (Node.Standard.get_dynamic_adjusted code) n

    let edge_distance a b = 
        apply_f_on_lazy Node.Standard.edge_distance a b

    let support_chars starting code n =
        apply_single_f_on_lazy (Node.Standard.support_chars starting code) n

    let n_chars ?(acc=0) n =
        apply_single_f_on_lazy (Node.Standard.n_chars ~acc) n

    let prioritize x = x

    let reprioritize _ x = x

    let f_codes ints n = 
        lazy_from_val (apply_single_f_on_lazy (Node.Standard.f_codes ints) n)

    let min_child_code code n = 
        apply_single_f_on_lazy (Node.Standard.min_child_code code) n

    type nad8 = Node.Standard.nad8

    let new_characters = Node.Standard.new_characters

    let build_node static_characters chars n = 
        lazy_from_val 
        (apply_single_f_on_lazy 
        (Node.Standard.build_node static_characters chars) 
        n)

    let set_exclude_info e n = 
        lazy_from_val 
        (apply_single_f_on_lazy (Node.Standard.set_exclude_info e) n)

    let excludes_median code a b = 
        apply_f_on_lazy (Node.Standard.excludes_median code) a b

    let has_excluded = Node.Standard.has_excluded

    module T = struct

        let add_exclude set = 
            fun x ->
                lazy_from_val (Node.Standard.T.add_exclude set 
                (force_val x))

        let remove_exclude = 
            fun x ->
                lazy_from_val (Node.Standard.T.remove_exclude 
                (force_val x))
    end

    module Union = struct
        type u = Node.Standard.Union.u my_lazy

        let union code a b c =
            let a = force_val a 
            and b = force_val b 
            and c = force_val c in
            lazy_from_val (Node.Standard.Union.union code a b c)

        let union_one which code a b =
            let a = force_val a
            and b = force_val b in
            lazy_from_val (which code a b)

        let union_preliminary a b c =
            union_one Node.Standard.Union.union_preliminary a b c

        let union_final a b c =
            union_one Node.Standard.Union.union_final a b c

        let leaf taxon_code code a = 
            lazy_from_val 
                (Node.Standard.Union.leaf taxon_code code (force_val a))

        let distance a b = 
            Node.Standard.Union.distance 
            (force_val a)
            (force_val b)

        let saturation a = 
            Node.Standard.Union.saturation (force_val a)

        let distance_node code n u =
            Node.Standard.Union.distance_node code
            (force_val n)
            (force_val u)

        let compare a b =
            Node.Standard.Union.compare (force_val a) (force_val b)

        let get_sequence a b c = 
            Node.Standard.Union.get_sequence a b (force_val c)
    end

    let for_support a b c d : n list =
        let nb = 
            List.map (fun (a, (b : n)) ->
                a, force_val b) b
        in
        let res = Node.Standard.for_support a nb c d in
        List.map lazy_from_val res

    let root_cost a = 
        Node.Standard.root_cost (force_val a)

    let extra_cost_from_root a =
        Node.extra_cost_from_root (force_val a)

    let tree_cost a b =
        (total_cost a b) +. (root_cost b)

    let to_single root a b c d set =
        let root' = force_opt root
        and b' = (force_val b)
        and d' = (force_val d) in
        lazy_from_val (Node.Standard.to_single root' a b' c d' set)

    let force x = force x
end

let q_print n = 
    let adjusted_data_lst = match n.adjusted with
        | None -> []
        | Some x -> [x]
    in
    let taxon_code n = match n.unadjusted @ adjusted_data_lst with
        | h :: _ -> h.code
        | []     -> failwith "AllDirNode.taxon_code"
    in
    Printf.printf "Node %d has adjusted|unadjusted:\t" (taxon_code n);
    List.iter (fun x -> match x.dir with
                    | Some (a,b) -> Printf.printf "(%d,%d) " a b
                    | None       -> Printf.printf "none ")
              adjusted_data_lst;
    print_string " | ";
    List.iter (fun x -> match x.dir with
                    | Some (a,b) -> Printf.printf "(%d,%d) " a b
                    | None       -> Printf.printf "none ")
              n.unadjusted;
    Printf.printf "\n\tUnadjusted Data:\n";
    List.iter (fun x ->
                print_string "\t\t";
                Node.print_times (force_val x.lazy_node))
              n.unadjusted;
    Printf.printf "\n\tAdjusted Data:\n";
    List.iter (fun x ->
                print_string "\t\t";
                Node.print_times (force_val x.lazy_node))
              adjusted_data_lst;
    print_newline ()


let using_static_likelihood a = match a.unadjusted with
    | x::_ -> Node.using_likelihood `Static (force_val x.lazy_node)
    | []   -> (* no data; impossible *) assert( false );

module AllDirF : NodeSig.S with type e = exclude
                           with type n = node_data
                           with type other_n = Node.Standard.n
                           with type nad8 = Node.Standard.nad8 =
struct

    type n = node_data
    type other_n = Node.Standard.n

    let to_other x = match x.unadjusted with
        | [x] -> OneDirF.to_other x.lazy_node
        | _ -> failwith "illegal argument"

    type e = exclude 

    let to_n_nodir node = 
        let node_dir = {
            lazy_node = to_n node;
            dir = None;
            code = node.Node.taxon_code; }
        in
        { unadjusted = [node_dir]; adjusted = (Some node_dir)}

    let force x =
        let force_item lst = 
            List.map (fun x -> { x with lazy_node = OneDirF.force x.lazy_node }) lst
        in
        let force_item2 item =
            match item with 
            | None -> None
            | Some x -> Some { x with lazy_node = OneDirF.force x.lazy_node }
        in
        { unadjusted = force_item x.unadjusted; 
          adjusted = force_item2 x.adjusted }

    let get_something f code n =
        let node = match code with
            | Some code -> not_with code n
            | None -> match n with
                | [x] -> x
                | [] -> failwith "AllDirNode.num_height 1"
                | _ -> failwith "AllDirNode.num_height 2"
        in
        f code node.lazy_node

    (* replace a direction in a set *)
    let replace_dir node nodes =
        let found = ref false in
        let nodes = 
            List.map
                (fun ns -> match ns.dir,node.dir with
                    | None, None ->
                        assert( not !found ); found := true;
                        node
                    | Some (a,b), Some (c,d) when a = c && b = d -> 
                        assert( not !found ); found := true;
                        node
                    | Some (a,b), Some (c,d) when a = d && b = c ->
                        assert( not !found ); found := true;
                        node
                    | _ -> ns)
                nodes
        in
        assert( !found );
        nodes

    let compare a b =
        let rec dir_compare lst1 lst2 =
            match lst1, lst2 with
            | h1 :: t1, h2 :: t2 ->
                    if 0 = compare h1.dir h2.dir then
                        if 0 = compare h1.code h2.code then
                            let r =
                                OneDirF.compare h1.lazy_node h2.lazy_node
                            in
                            if 0 = r then dir_compare t1 t2
                            else r
                        else compare h1.code h2.code
                    else compare h1.dir h2.dir
            | [], [] -> 0
            | [], _ -> -1
            | _, [] -> 1
        in
        dir_compare a.unadjusted b.unadjusted

    let recode_anode f n = 
        { 
            lazy_node = OneDirF.recode f n.lazy_node;
            code = f n.code;
            dir = 
                match n.dir with
                | None        -> None
                | Some (a, b) -> Some (f a, f b);
        }

    let taxon_code n = 
        match n.unadjusted with
        | h :: _ -> (** All the elements in the list have the same code *)
                h.code
        | [] -> failwith "AllDirNode.taxon_code"

    let recode f n = { 
        unadjusted = List.map (recode_anode f) n.unadjusted;
        adjusted = match n.adjusted with
        | None -> None
        | Some x -> Some (recode_anode f x) 
    }

    let load_data ?(is_fixedstates=false) ?(silent=true) ?(classify=true) data = 
        let data, nodes = Node.Standard.load_data ~is_fixedstates ~classify data in
        data, List.map to_n_nodir nodes

    let fix_preliminary x = x

(* TODO: nobody calls this now -- move functions from allDirChar here *)
    let to_single root a b c d set =
        let b',d' = match a,c with
            | None,None ->
                begin
                    let one = match b.adjusted with
                        | None -> failwith "allDirNode,to_single,no adjusted data"
                        | Some x -> x 
                    and two = match d.adjusted with
                        | None -> failwith "allDirNode,to_single,no adjusted data"
                        | Some y -> y 
                    in one,two
                end
            | Some x,Some y -> not_with x b.unadjusted, not_with y d.unadjusted
            | _,_ -> failwithf "AllDirNode.to_single: Ambiguous directions"
        in

        let root = match root with 
            | Some r -> Some ((List.hd r.unadjusted).lazy_node)
            | None -> None
        in
        let lazy_node = OneDirF.to_single root a (b'.lazy_node) c (d'.lazy_node) set in
        let node = { d' with lazy_node = lazy_node } in
        { unadjusted = [node]; adjusted = Some node }

    let distance ?(para=None) ?(parb=None) missing_distance a b =
        let one_f = OneDirF.distance missing_distance in
        match a.unadjusted, b.unadjusted with
        | [a], [b] -> one_f a.lazy_node b.lazy_node
        | a, [b] -> 
                (match parb with
                | None -> failwith "AlldNode.distance 1"
                | Some para -> 
                        let a = not_with para a in
                        one_f a.lazy_node b.lazy_node)
        | [a], b -> 
                (match para with
                | None -> failwith "AlldNode.distance 2"
                | Some parb -> 
                        let b = not_with parb b in
                        one_f a.lazy_node b.lazy_node)
        | a, b -> 
                (match para, parb with
                | Some para, Some parb -> 
                        let a = not_with para a 
                        and b = not_with parb b in
                        one_f a.lazy_node b.lazy_node
                | _ -> failwith "AlldNode.distance 3")

    let min_child_code code n = 
        get_something OneDirF.min_child_code code n.unadjusted

    let apply_on_one_direction f x n =
        match x, n.unadjusted with
        | None, [n] -> 
                f x n.lazy_node
        | Some cx, n ->
                let n = not_with cx n in
                f x n.lazy_node
        | _ -> assert false

    let character_costs = apply_on_one_direction OneDirF.character_costs 

    let get_nonadd_8  = apply_on_one_direction OneDirF.get_nonadd_8
    let get_nonadd_16 = apply_on_one_direction OneDirF.get_nonadd_16
    let get_nonadd_32 = apply_on_one_direction OneDirF.get_nonadd_32
    let get_addgen    = apply_on_one_direction OneDirF.get_addgen
    let get_addvec    = apply_on_one_direction OneDirF.get_addvec
    let get_sank      = apply_on_one_direction OneDirF.get_sank
    let get_fixedstates = apply_on_one_direction OneDirF.get_fixedstates
    let get_dynamic   = apply_on_one_direction OneDirF.get_dynamic
    let get_mlstatic  = apply_on_one_direction OneDirF.get_mlstatic

    (* likelihood only function to optimize the branch length *)
    let edge_iterator par mine ch1 ch2 =
        let get_dir p n = (not_with (taxon_code p) n.unadjusted).lazy_node in
        let atom, btom, parofm = match par with
            | Some x -> get_dir mine ch1, get_dir mine ch2, Some (get_dir mine x)
            | None   -> get_dir ch2 ch1, get_dir ch1 ch2, None
        and ab_to_m = 
            (with_both (taxon_code ch1) (taxon_code ch2) mine.unadjusted)
        in
        let node_dir  = 
            { ab_to_m with
                lazy_node = OneDirF.edge_iterator parofm ab_to_m.lazy_node atom btom
            }
        in
        { unadjusted = [node_dir]; adjusted = None; }

    (* calculate the median between a and b. old can be used as a heuristic,
     * branches are the supplied branch lengths of the children a and b, *)
    let median ?branches my_code old a b =
        let na, nb,code = match my_code with
            | Some code ->
                let in_a = not_with code a.unadjusted
                and in_b = not_with code b.unadjusted in
                in_a, in_b, code
            | None ->
                let in_a =
                    try not_with (taxon_code b) a.unadjusted
                    with Not_found -> q_print a; q_print b;
                        failwithf "Direction in A to %d; for root" (taxon_code b)
                and in_b =
                    try not_with (taxon_code a) b.unadjusted
                    with Not_found -> q_print a; q_print b;
                        failwithf "Direction in B to %d; for root" (taxon_code a)
                in
                in_a, in_b, -1
        in
        let old = match old with
            | Some oldness ->
                begin try 
                    let n = with_both na.code nb.code oldness.unadjusted in
                    Some (n.lazy_node)
                 with | Not_found -> None end
            | None -> None
        in
        let lazy_node =
                OneDirF.median ?branches my_code old na.lazy_node nb.lazy_node
        in
        let node = {
            lazy_node = lazy_node;
            dir = Some (na.code, nb.code);
            code = code;
        } in
        { unadjusted = [node]; adjusted = None; }

    (* The final_states is calculated with the unadjusted component of the vertex *)
    let final_states grandcode par cur a b =
        let get_desired_dir par x = not_with (taxon_code par) x.unadjusted in
        let na = get_desired_dir cur a
        and nb = get_desired_dir cur b 
        and ncur = get_desired_dir par cur 
        and npar = match grandcode with
            | Some x -> not_with x par.unadjusted
            | None   ->
                begin match par.unadjusted with
                    | [par] -> par
                    | _     -> failwith "AllDirNode.final_states"
                end
        in
        let node = {
            lazy_node = 
                OneDirF.final_states None 
                        npar.lazy_node ncur.lazy_node na.lazy_node nb.lazy_node;
            dir = Some (na.code, nb.code);
            code = taxon_code cur;
        } in
        match cur.unadjusted with
        | [_] -> { cur with unadjusted = [node] }
        | _   -> let x, y = yes_with (taxon_code par) cur.unadjusted in
                 { cur with unadjusted = [x; y; node] }

    (** [get_times_between child par] 
     *
     * returns the branch length between the child and parent, contained in
     * parent, although, it shouldn't matter as long as the directions are
     * available *)
    let get_times_between (child:n) (par: n option) =
        match par with
        | None ->
            begin match child.unadjusted with
                | [x] -> OneDirF.get_times_between x.lazy_node None
                | _   -> failwith "A direction should be specified"
            end
        | Some par ->
            (* Use the leaf if available *)
            begin match child.unadjusted, par.unadjusted with
                | _,[x] -> OneDirF.get_times_between x.lazy_node None
                | [x],_ -> OneDirF.get_times_between x.lazy_node None
                |  _, _ ->
                    try let child = not_with (taxon_code par) child.unadjusted
                        and par = either_with (taxon_code child) par.unadjusted in
                        OneDirF.get_times_between child.lazy_node (Some par.lazy_node)
                    with | _ ->
                        let child = either_with (taxon_code child) par.unadjusted
                        and par   = not_with (taxon_code par) child.unadjusted in
                        OneDirF.get_times_between par.lazy_node (Some child.lazy_node)
            end


    (** [extract_states par child] extract the states of child toward par *)
    let extract_states alph data parc codes mine = 
        let n = match parc with
            | None -> begin match mine.unadjusted with
                | [x] -> x
                |  _  -> failwith "AllDirNode.extract_states; No Direction" end
            | Some x -> not_with x mine.unadjusted
        in
        OneDirF.extract_states alph data None codes n.lazy_node

    (** [apply_time child parent] applies time from parent into child --used on
     * leaves when the uppass_heuristic doesn't normally run over internal nodes
     *
     * first argument is for when time is from root median *)
    let apply_time root child parent =
        if uppass_debug then
            info_user_message "Applying time from %d to %d"
                    (taxon_code parent) (taxon_code child);
        let pc = taxon_code parent and mc = taxon_code child in
        let time_M2P =
            lazy_from_fun 
                (fun () -> 
                    if root then
                        let timedat = (either_with mc parent.unadjusted).lazy_node in
                        Node.get_times_between (force_val timedat) None
                    else
                        let timedat = (either_with mc parent.unadjusted).lazy_node in
                        Node.get_times_between (force_val timedat) 
                                (Some (min_child_code (Some pc) child)) )
        and data_m2p = match child.unadjusted with
            | [x] -> x.lazy_node
            |  _  -> failwith "Apply time is used on leaves."
        in
        let node = lazy_from_fun
            (fun () -> 
                Node.replace_parent_time (force_val data_m2p) (force_val time_M2P))
        in
        let node = {
            lazy_node = node;
            dir = None;
            code = taxon_code child; }
        in
        { unadjusted = [node]; adjusted = Some node; }
     

    (** [uppass_heuristic tbl curr ch1 ch2 par r]
     * [cur], at this point will have one direction ([ch1],[ch2]), after this call
     * it will fill in the other directions, lazily, with [par]. This process will
     * update the nodes, the times (likelihood), et cetera, with lazy values.
     *
     * If a node [has] more directions then it should have all three and is
     * based around the root (since that time data has been filled in prior),
     * the new nodes should be the same (with an assertion). [r] will hold the
     * root median between m and p, and is used in this case.
     *
     * time_M2A [B]
     *      __  //
     *    [A]==[M] < data_m2a 
     *          \\
     *          [P]
    **)
    let uppass_heuristic p_data p_time m a b =
        if uppass_debug then 
            info_user_message "Performing uppass Heuristic on %d with ch1=%d,ch2=%d,p=%d"
                (taxon_code m) (taxon_code a) (taxon_code b) (taxon_code p_data);
        (* open all the taxon codes; they get used often for directions *)
        let mc = taxon_code m and ac = taxon_code a
        and bc = taxon_code b and pc = taxon_code p_data
        and get_dir parc x = (not_with parc x.unadjusted).lazy_node in
        let data_m2p = 
            try match m.unadjusted with
            (* then it hasn't been resolved by an earlier uppass on root *)
            | [x] -> 
                if (match x.dir with
                    | None -> false
                    | Some (xa,xb) -> (xa=ac && xb=bc) || (xa=bc && xb=ac))
                  then x.lazy_node
                  else raise Not_found
            (* ...was resolved, so get the direction *)
            |  _  -> get_dir pc m 
            with | Not_found -> 
                q_print m;
                failwithf "Data from %d toward parent %d is missing.\n" mc pc
        and data_p2m = match p_data.unadjusted with
            (* then it HASN'T been resolved by an earlier uppass, and is truely
             * the parent with a calculated/lazy median *)
            | [x] -> x.lazy_node
            (* ...was resolved, so get one of the directions that has m as child *)
            |  _  -> get_dir mc p_data
        (* both of these should exist *)
        and data_b2m = get_dir mc b
        and data_a2m = get_dir mc a in
        (* get the times in all directions --AB have them in M and M has it in P *)
        let smc = Some mc in
        let time_M2A = 
            lazy_from_fun (fun () ->
                Node.get_times_between 
                    (force_val data_m2p) (Some (min_child_code smc a)) )
        and time_M2B = 
            lazy_from_fun (fun () ->
                Node.get_times_between 
                    (force_val data_m2p) (Some (min_child_code smc b)) )
        and time_M2P =
            lazy_from_fun (fun () -> match p_time with
                | Some timedat ->
                    let timedat = force_val (either_with mc timedat.unadjusted).lazy_node in
                    Node.get_times_between timedat None
                | None -> 
                    let timedat = force_val (either_with mc p_data.unadjusted).lazy_node in
                    Node.get_times_between timedat (Some (min_child_code (Some pc) m)))
        in
        (* Printf.printf "(%d--%d):%a\t(%d--%d):%a\t(%d--%d):%a\n%!"
           mc ac pp_opt_list (force_val time_M2A) mc bc pp_opt_list
           (force_val time_M2B) mc pc pp_opt_list (force_val time_M2P);*)
        (* call medians with times supplied *)
        let node_A = lazy_from_fun
            (fun () ->
                let data_m2p = force_val data_m2p in
                let data_p2m = force_val data_p2m in
                let time_M2P = force_val time_M2P in
                let time_M2B = force_val time_M2B in
                let time_M2A = force_val time_M2A in
                let data_b2m = force_val data_b2m in
                Node.median_w_times (Some mc) (Some (data_m2p)) data_p2m
                                    data_b2m time_M2P time_M2B (Some time_M2A))
        and node_B = lazy_from_fun
            (fun () -> 
                let data_m2p = force_val data_m2p in
                let data_p2m = force_val data_p2m in
                let time_M2P = force_val time_M2P in
                let time_M2B = force_val time_M2B in
                let time_M2A = force_val time_M2A in
                let data_a2m = force_val data_a2m in
                Node.median_w_times (Some mc) (Some (data_m2p)) data_p2m
                                    data_a2m time_M2P time_M2A (Some time_M2B))
        and node_C = lazy_from_fun
            (fun () -> 
                Node.replace_parent_time (force_val data_m2p) (force_val time_M2P))
        in
        (* build the nodes *)
        let allDir = 
            let dir_A= { code= mc; lazy_node= node_A; dir= Some(bc,pc); }
            and dir_B= { code= mc; lazy_node= node_B; dir= Some(ac,pc); }
            and dir_C= { code= mc; lazy_node= node_C; dir= Some(ac,bc); } in
            [ dir_A ; dir_B ; dir_C ]
        in
        let res = { unadjusted = allDir; adjusted = None } in
        if uppass_debug then
	    begin
            info_user_message
        "End of Performing uppass Heuristic on %d with ch1=%d,ch2=%d,p=%d"
        (taxon_code m) (taxon_code a) (taxon_code b) (taxon_code p_data);
        (*print_node_data a true;
        print_node_data b true;*)
        (*info_user_message "return node#.%d with node data:" (taxon_code m);
		print_node_data res true;
	*)  
	  end;
        res

    (** adjust the branches in the tree, including branch lengths, uses
        adjusted from a single assignment **)
    let readjust mode to_adjust ch1 ch2 par mine =
        let errmsg =  "allDirChar,readjust,no adj-data" in
        let c1_adj = (get_adjusted_nodedata ch1 errmsg).lazy_node
        and c2_adj = (get_adjusted_nodedata ch2 errmsg).lazy_node
        and pr_adj = (get_adjusted_nodedata par errmsg).lazy_node
        and my_adj = (get_adjusted_nodedata mine errmsg).lazy_node in
        let a1,modi = OneDirF.readjust mode to_adjust c1_adj c2_adj pr_adj my_adj in
        let node_dir =
            {
                lazy_node = a1;
                dir = Some ((taxon_code ch1), (taxon_code ch2));
                code = taxon_code mine;
            } 
        in
        let node =
            { adjusted = Some node_dir;
              unadjusted = replace_dir node_dir mine.unadjusted; }
        in
        (node,modi)

    let to_string nodes =
        let adj_data = get_adjusted_nodedata nodes "allDirNode,to_string no adj-data" in
        let res =
            List.map 
                (fun x -> OneDirF.to_string x.lazy_node)
                nodes.unadjusted
            @
            [OneDirF.to_string adj_data.lazy_node]
        in
        String.concat "\n" res

    let total_cost par n = match par with
        | Some code ->
            OneDirF.total_cost par (not_with code n.unadjusted).lazy_node
        | None -> match n.unadjusted with
            | [x] -> OneDirF.total_cost par x.lazy_node
            | []  -> failwith "The emtpy median? AllDirNode.total_cost"
            | _   -> failwith "AllDirNode.total_cost"

    let min_prior chars n =
        let lst = match n.adjusted with | None -> n.unadjusted | Some x -> [x] in
        let lst = List.map (fun x -> OneDirF.min_prior chars x.lazy_node) lst in
        List.fold_left (fun acc x -> min acc x) infinity lst

    let node_cost par n = match par with
        | Some code ->
                OneDirF.node_cost par (not_with code n.unadjusted).lazy_node
        | None -> match n.unadjusted with
            | [x] -> OneDirF.node_cost par x.lazy_node
            | [] -> failwith "The emtpy median? AllDirNode.node_cost"
            | _ -> failwith "AllDirNode.node_cost"

    let update_leaf x = x

    let union_distance _ _ = 0.0

    let rec is_collapsable clas a b =
        let acode = taxon_code a
        and bcode = taxon_code b in
        match clas with
        | `Static ->
                let da = not_with bcode a.unadjusted
                and db = not_with acode b.unadjusted in
                OneDirF.is_collapsable `Static da.lazy_node db.lazy_node
        | `Dynamic ->
                let errmsg = "allDirNode,is_collapsable,no adj-data" in
                let da = get_adjusted_nodedata a errmsg 
                and db = get_adjusted_nodedata b errmsg in
                OneDirF.is_collapsable `Dynamic da.lazy_node db.lazy_node
        | `Any ->
                let s = try is_collapsable `Static a b  with | _ -> false
                and d = try is_collapsable `Dynamic a b with | _ -> false in
                (s && d)

    let classify_data l1 n1 l2 n2 chars acc =
        match n1.adjusted, n2.adjusted with
        | Some x, Some y ->
            OneDirF.classify_data l1 x.lazy_node l2 y.lazy_node chars acc
        | (Some _ | None), _ ->
            failwith "Cannot process non-final data in classification"

    let to_xml _ _ _ = ()


    let run_any f n = match n with
        | h :: _ -> lazy_from_val (f h.lazy_node)
        | [] -> failwith "AllDirNode.run_any"

    let run_all f n =
        let processor x = { x with lazy_node = f x.lazy_node } in
        let adj_data = match n.adjusted with
            | None   -> None
            | Some x -> Some (processor x)
        in
        { unadjusted = List.map processor n.unadjusted;
            adjusted = adj_data }

    let num_height code n = 
        get_something OneDirF.num_height code n.unadjusted

    let num_otus code n =
        get_something OneDirF.num_otus code n.unadjusted

    let get_sequences code n = 
        let adj_data = get_adjusted_nodedata n "allDirNode,get_sequences,no adj-data" in
        get_something OneDirF.get_sequences code [adj_data]

    let get_dynamic_preliminary code n =
        get_something OneDirF.get_dynamic_preliminary code n.unadjusted

    let get_dynamic_adjusted code n =
        let adj_data = get_adjusted_nodedata n "allDirNode,get_dynamic_sequences,no adj-data" in
        get_something OneDirF.get_dynamic_preliminary code [adj_data]

    let edge_distance a b =
        let acode = taxon_code a
        and bcode = taxon_code b in
        let da = not_with acode b.unadjusted
        and db = not_with bcode a.unadjusted in
        OneDirF.edge_distance da.lazy_node db.lazy_node

    let support_chars starting code n =
        get_something (OneDirF.support_chars starting) code n.unadjusted

    let n_chars ?(acc=0) n =
        force_val (run_any (OneDirF.n_chars ~acc) n.unadjusted)

    let prioritize x = x

    let reprioritize _ x = x

    let get_others err code n = match code with
        | Some code -> Some (yes_with code n)
        | None -> None

    let get_node err code n = match code with
        | Some code -> not_with code n
        | None -> match n with 
            | [x] -> x
            | _ -> failwith err

    let f_codes ints n = 
        run_all (OneDirF.f_codes ints) n

    type nad8 = OneDirF.nad8

    let new_characters = Node.Standard.new_characters

    let build_node static_characters chars node =
        let processor x = 
            let res = OneDirF.build_node static_characters chars x.lazy_node in
            { x with lazy_node = res }
        in
        let adj_data = get_adjusted_nodedata node "allDirNode,build_node,no adj-data" in
        let uadj = List.map processor node.unadjusted
        and adj = processor adj_data in
        { unadjusted = uadj; adjusted = Some adj; }

    let set_exclude_info e x =
        let processor x = 
            let res = OneDirF.set_exclude_info e x.lazy_node in
            { x with lazy_node = res }
        in
        let adj_data = get_adjusted_nodedata x "allDirNode,set_exclude_info,no adj-data" in
        let uadj = List.map processor x.unadjusted
        and adj = processor adj_data in
        { unadjusted = uadj; adjusted = Some adj }

    let excludes_median code a b = 
        let err = "AllDirNode.excludes_median" in
        let nodea, nodeb = get_node err code a.unadjusted,
                           get_node err code b.unadjusted in
        OneDirF.excludes_median code nodea.lazy_node nodeb.lazy_node

    let has_excluded = Node.Standard.has_excluded

    module T = struct

        let add_exclude set n =
            let processor x = 
                { x with lazy_node = OneDirF.T.add_exclude set x.lazy_node}
            in
            let adj_data = get_adjusted_nodedata n "allDirNode,add_exclude,no adj-data" in
            let uadj = List.map processor n.unadjusted
            and adj = processor adj_data in
            { unadjusted = uadj; adjusted = Some adj }

        let remove_exclude n = 
            let processor x = 
                { x with lazy_node = OneDirF.T.remove_exclude x.lazy_node}
            in
            let adj_data = get_adjusted_nodedata n "allDirNode,remove_exclude,no adj-data" in
            let uadj = List.map processor n.unadjusted
            and adj = processor adj_data in
            { unadjusted = uadj; adjusted =Some adj }
    end

    module Union = struct

        type u = OneDirF.Union.u
        let union code n a b =
            let node = 
                get_node "AllDirNode.AllDirF.Union.union" code n.unadjusted 
            in
            OneDirF.Union.union code node.lazy_node a b

        let union_final code a b = 
            let b = 
                get_node "AllDirNode.AllDirF.Union.union_final" code b.unadjusted 
            in
            OneDirF.Union.union_final code a b.lazy_node

        let union_preliminary code a b = 
            let b = 
                get_node "AllDirNode.AllDirF.Union.union_preliminary" code b.unadjusted 
            in
            OneDirF.Union.union_preliminary code a b.lazy_node

        let leaf taxon_code code n =
            let x = get_node "AllDirNode.AllDirF.Union.leaf" code n.unadjusted in
            OneDirF.Union.leaf taxon_code code x.lazy_node

        let distance a b = 
            let debug = false in
            if debug then
                Printf.printf "AllDirNode.AllDirF.Union.distance,%!";
            OneDirF.Union.distance a b

        let saturation = OneDirF.Union.saturation

        let distance_node code n u =
            let node = 
                get_node "AllDirNode.AllDirF.Union.distance_node" code n.unadjusted
            in
            OneDirF.Union.distance_node code node.lazy_node u

        let compare = OneDirF.Union.compare

        let get_sequence = OneDirF.Union.get_sequence

    end

    let for_support a b c d : n list=
        let nb = 
            List.map (fun (a, (b : n)) ->
                match b.unadjusted with
                | [b] -> a, force_val b.lazy_node
                | _ -> failwith "AllDirNode.for_support") b 
        in
        let res = Node.Standard.for_support a nb c d in
        let rec merge a b =
            match a, b with
            | (_, [ha]) :: ta, hb :: tb ->
                    [{ ha with lazy_node = lazy_from_val hb }] ::
                        merge ta tb
            | [], [] -> []
            | _, _ -> failwith "AllDirNode.for_support 2"
        in
        let b = List.map (fun (a, x) -> a, x.unadjusted) b in
        let res = merge b res in
        List.map (fun x -> { unadjusted = x; adjusted = None; }) res

    let root_cost a = match a.adjusted with
        | None   -> assert( using_static_likelihood a ); 0.0
        | Some x -> OneDirF.root_cost x.lazy_node

    let extra_cost_from_root a = match a.adjusted with
        | None -> assert( using_static_likelihood a ); 0.0
        | Some x ->  OneDirF.extra_cost_from_root x.lazy_node
    
    let tree_cost a b = 
        (total_cost a b) +. (root_cost b)

end

type 'a node_hybrid = {
    st : Node.Standard.n option;
    dy : 'a;
}

module HybridF = struct
    let get_dynamic x = x.dy
end

let create_root_from_child_branch child parent : AllDirF.n =
    let get_dir p n = not_with (AllDirF.taxon_code p) n.unadjusted in
    let child_w_time = (get_dir parent child).lazy_node
    and parent_w_out = (get_dir child parent).lazy_node in
    let lnode = 
        lazy_from_fun
            (fun () ->
                Node.median_of_child_branch None 
                    (force_val child_w_time) (force_val parent_w_out))
    in
    let dir = (AllDirF.taxon_code parent, AllDirF.taxon_code child) in
    let node = { lazy_node = lnode; dir = Some dir; code = -1; } in
    { unadjusted = [node]; adjusted = Some node; }

let create_root ?branches a aa ab b ba bb opt =
    let middle = match opt with
        | Some x -> x
        | None   -> AllDirF.median ?branches None None a b
    in
    let a_final = match aa,ab with
        | Some aa,Some ab -> AllDirF.uppass_heuristic b (Some middle) a aa ab
        | None, None      -> AllDirF.apply_time true a middle
        | _, _            -> failwith "Failed with children of A @ Create Roots"
    and b_final = match ba,bb with
        | Some ba,Some bb -> AllDirF.uppass_heuristic a (Some middle) b ba bb
        | None, None      -> AllDirF.apply_time true b middle
        | _, _            -> failwith "Failed with children of B @ Create Roots"
    and l_middle = match middle.unadjusted with
        | [x] ->
            assert( match x.dir with | None -> false
                    | Some (xa,xb) ->
                        (xa = (AllDirF.taxon_code a) && xb = (AllDirF.taxon_code b)) ||
                        (xa = (AllDirF.taxon_code b) && xb = (AllDirF.taxon_code a)));
            x.lazy_node
        |  _  -> failwith "root median has multiple directions?"
    in
    a_final,b_final,l_middle

let create_root_w_times left right =
    let get_dir parc x = force_val (not_with parc x.unadjusted).lazy_node
    and get_a_dir child x = match x.unadjusted with
        | [x] -> assert( match x.dir with | None -> true | _ -> false); 
                 force_val x.lazy_node
        |  x  -> force_val (either_with child x).lazy_node
    in
    lazy_from_fun
        (fun () ->
            let l_code = AllDirF.taxon_code left
            and r_code = AllDirF.taxon_code right in
            let left2right = get_dir r_code left
            and right2left = get_dir l_code right in
            let in_l_time =
                Node.get_times_between (get_a_dir r_code left)
                    (Some (AllDirF.min_child_code (Some l_code) right))
            and in_r_time = 
                Node.get_times_between (get_a_dir l_code right)
                    (Some (AllDirF.min_child_code (Some r_code) left))
            in
(*            Printf.printf "Edge (%d,%d) : %a -- %a\n%!" l_code r_code*)
(*                          pp_opt_list in_l_time pp_opt_list in_r_time;*)
            Node.median_w_times 
                None None left2right right2left in_l_time in_r_time None)

