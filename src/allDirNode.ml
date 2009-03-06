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

let () = SadmanOutput.register "AllDirNode" "$Revision: 1616 $"

let eager = false
let uppass_debug = false

type exclude = Node.exclude

type direction = (int * int) option

type 'a my_lazy = 
    | Lazy of 'a Lazy.t
    | Eager of 'a

type a_node = Node.node_data my_lazy

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
    adjusted : node_dir list;   (** An adjuted node value calculated after the
    downpass *)
}

let to_n node = Eager node

let failwithf format = Printf.ksprintf failwith format
let info_user_message format = 
    Printf.ksprintf (Status.user_message Status.Information) format
let error_user_message format = 
    Printf.ksprintf (Status.user_message Status.Error) format

let has_code code n =
    match n.dir with
    | None -> true
    | Some (a, b) ->
            a <> code && b <> code
(*
let get_code n =
    match n with
    | x :: _ ->
            let code = x.code in
            assert (List.for_all (fun x -> x.code = code) n);
            x.code
    | [] -> failwith "AllDirNode.get_code"

let print_pairs x = 
    let msg = 
        match x.dir with
        | None -> "No direction!"
        | Some (a, b) -> Printf.ksprintf (fun x -> x) "%d, %d" a b
    in
    Status.user_message Status.Error msg
*)

let yes_with code n =
    assert (0 <> List.length n);
    match List.filter (fun x -> not (has_code code x)) n with
    | [x; y] -> x, y
    | x -> raise Not_found

let not_with code n =
    assert (0 <> List.length n);
    match List.filter (has_code code) n with
    | [x] -> x
    | x -> raise Not_found

(* grabs a node in the direction with children c1 and c2 --Internal Nodes only *)
let with_both c1 c2 n = 
    assert( 3 = List.length n );
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

    let load_data ?(silent=true) ?(classify=true) data = 
        let data, nodes = Node.Standard.load_data ~classify data in
        data, List.map to_n nodes

    let fix_preliminary x = x

    let apply_f_on_lazy f a b = 
        let a = force_val a
        and b = force_val b in
        f a b

    let distance ?(para=None) ?(parb=None) missing_distance a b =
        apply_f_on_lazy (Node.Standard.distance missing_distance) a b

    let distance_of_type t missing_dist p_code m a b = 
        Node.Standard.distance_of_type t missing_dist p_code
                (force_val m) (force_val a) (force_val b)

    let character_costs x n = 
        Node.Standard.character_costs x (force_val n)

    let get_nonadd_8 x n =
        Node.Standard.get_nonadd_8 x (force_val n)

    let get_nonadd_16 x n =
        Node.Standard.get_nonadd_16 x (force_val n)

    let get_nonadd_32 x n =
        Node.Standard.get_nonadd_32 x (force_val n)

    let get_add x n =
        Node.Standard.get_add x (force_val n)

    let get_sank x n =
        Node.Standard.get_sank x (force_val n)

    let get_dynamic x n =
        Node.Standard.get_dynamic x (force_val n)

    let get_mlstatic x n =
        Node.Standard.get_mlstatic x (force_val n)

    let dump_node f x y = Node.Standard.dump_node f (force_val x) (force_val y)

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

    let readjust gp mode set_2adjust a b c mine =
        let node,set = Node.Standard.readjust gp mode set_2adjust
                            (force_val a) (force_val b) (force_val c) (force_val mine)
        in
        (to_n node,set)

    let final_states x par cur a b = 
        lazy_from_fun 
        (fun () -> apply_f_on_lazy
            (Node.Standard.final_states x (force_val par) (force_val cur))
            a
            b
        )

    let apply_time root x y =
        lazy_from_fun (fun () -> apply_f_on_lazy (Node.Standard.apply_time root) x y)

    let estimate_time x y =
        Node.Standard.estimate_time (force_val x) (force_val y)

    let uppass_heuristic pcode mine a b par = mine

    let to_string v = Node.Standard.to_string (force_val v)

    let apply_single_f_on_lazy f a = 
        f (force_val a)

    let total_cost x n =
        apply_single_f_on_lazy (Node.Standard.total_cost x) n

    let node_cost par n = 
        apply_single_f_on_lazy (Node.Standard.node_cost par) n

    let update_leaf x = x

    let taxon_code n = 
        apply_single_f_on_lazy Node.Standard.taxon_code n

    let union_distance _ _ = 0.0

    let is_collapsable clas a b = 
        apply_f_on_lazy (Node.Standard.is_collapsable clas) a b

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

    let tree_cost x n =
        (total_cost x n) +. (root_cost n)

    let to_single root a b c d set =
        let root' = force_opt root
        and b' = (force_val b)
        and d' = (force_val d) in
        lazy_from_val (Node.Standard.to_single root' a b' c d' set)

    let force x = force x
end

let q_print n = 
    let taxon_code n = match n.unadjusted with
        | h :: _ -> h.code
        | [] -> failwith "AllDirNode.taxon_code"
    in
    Printf.printf "Node %d has:\t" (taxon_code n);
    List.iter (fun x -> match x.dir with 
                    | Some (a,b) -> Printf.printf "(%d,%d) " a b
                    | None -> Printf.printf "none ") n.adjusted;
    print_string " | ";
    List.iter (fun x -> match x.dir with 
                    | Some (a,b) -> Printf.printf "(%d,%d) " a b
                    | None -> Printf.printf "none ") n.unadjusted;
    print_newline ()

module AllDirF : NodeSig.S with type e = exclude with type n = node_data with
type other_n = Node.Standard.n with
type nad8 = Node.Standard.nad8 = struct

    type n = node_data
    type other_n = Node.Standard.n

    let to_other x = 
        match x.unadjusted with
        | [x] -> OneDirF.to_other x.lazy_node
        | _ -> failwith "illegal argument"

    type e = exclude 

    let to_n_nodir node = 
        let node_dir = {
            lazy_node = to_n node;
            dir = None;
            code = node.Node.taxon_code;
        }
        in
        { unadjusted = [node_dir]; adjusted = [node_dir]}


    let force x =
        let force_item lst = 
            List.map (fun x -> { x with lazy_node = OneDirF.force x.lazy_node }) lst
        in
        { unadjusted = force_item x.unadjusted; adjusted = force_item x.adjusted }

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
                | None -> Printf.printf "\t Making None with %d -> %d\n%!" n.code (f n.code);None
                | Some (a, b) -> Some (f a, f b);
        }

    let recode f n = { 
        unadjusted = List.map (recode_anode f) n.unadjusted;
        adjusted = List.map (recode_anode f) n.adjusted;
    }

    let load_data ?(silent=true) ?(classify=true) data = 
        let data, nodes = Node.Standard.load_data ~classify data in
        data, List.map to_n_nodir nodes

    let fix_preliminary x = x

    let to_single root a b c d set =
        let b',d' = match a,c with
            | None,None -> (match b.adjusted,d.adjusted with
                    | [x],[y] -> x,y
                    | x,y -> failwithf "AllDirNode.to_single: no GC but %d %d direction(s)"
                                        (List.length x) (List.length y)
                    )
            | Some x,Some y -> not_with x b.unadjusted,not_with y d.unadjusted
            | _,_ -> failwithf "AllDirNode.to_single: Ambiguous directions"
        in

        let root = match root with 
            | Some r -> Some ((List.hd r.unadjusted).lazy_node)
            | None -> None
        in  
        let lazy_node = OneDirF.to_single root a (b'.lazy_node) c (d'.lazy_node) set in
        let node = [{ d' with lazy_node = lazy_node }] in
        { unadjusted = node; adjusted = node }

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

    let apply_on_one_direction f x n =
        match x, n.unadjusted with
        | None, [n] -> 
                f x n.lazy_node
        | Some cx, n ->
                let n = not_with cx n in
                f x n.lazy_node
        | _ -> assert false

    let character_costs = apply_on_one_direction OneDirF.character_costs 

    let get_nonadd_8 = apply_on_one_direction OneDirF.get_nonadd_8
    let get_nonadd_16 = apply_on_one_direction OneDirF.get_nonadd_16
    let get_nonadd_32 = apply_on_one_direction OneDirF.get_nonadd_32
    let get_add = apply_on_one_direction OneDirF.get_add
    let get_sank = apply_on_one_direction OneDirF.get_sank
    let get_dynamic = apply_on_one_direction OneDirF.get_dynamic
    let get_mlstatic = apply_on_one_direction OneDirF.get_mlstatic

    let taxon_code n = 
        match n.unadjusted with
        | h :: _ -> (** All the elements in the list have the same code *)
                h.code
        | [] -> failwith "AllDirNode.taxon_code"

    let dump_node f child par = 
        OneDirF.dump_node f
        (not_with (taxon_code par) child.unadjusted).lazy_node
        (either_with (taxon_code child) par.unadjusted).lazy_node

    let distance_of_type t missing_dist p_code mine ch1 ch2 =
        let get_dir p n = (not_with p n.adjusted).lazy_node in
        OneDirF.distance_of_type
                t missing_dist p_code 
                (get_dir p_code mine)
                (get_dir (taxon_code mine) ch1)
                (get_dir (taxon_code mine) ch2)

    let edge_iterator par mine ch1 ch2 = 
        let get_dir p n = (not_with (taxon_code p) n.adjusted).lazy_node in
        let atom, btom, parofm = match par with
                    | Some x -> 
                        get_dir mine ch1, get_dir mine ch2, Some (get_dir mine x)
                    | None -> 
                        get_dir ch2 ch1, get_dir ch1 ch2, None
        and ab_to_m = (with_both (taxon_code ch1) (taxon_code ch2) mine.adjusted) in
        
        let node_dir  = 
            { ab_to_m with
                lazy_node = OneDirF.edge_iterator parofm ab_to_m.lazy_node atom btom
            } in
        { unadjusted = [node_dir]; adjusted = [node_dir]; }

    (* adjust the branches in the tree, including branch lengths *)
    let readjust gp_opt mode to_adjust ch1 ch2 par mine =
        (* in [n], we want the direction toward [p], the parent *)
        let get_dir p_code n = (not_with (taxon_code p_code) n.adjusted).lazy_node in
        let mine2par = match gp_opt with
            | Some gpc -> (not_with gpc par.adjusted).lazy_node
            | None -> 
                (match par.adjusted with
                    | [x] -> x.lazy_node
                    |  _  -> failwithf "AllDirNode.readjust got no GP but with %d directions."
                                                (List.length par.adjusted)
                )
        in

        let a1,modified = OneDirF.readjust gp_opt mode to_adjust (get_dir mine ch1) 
                                (get_dir mine ch2) mine2par (get_dir par mine)
        in

        let node_dir = {
                lazy_node = a1;
                dir = Some( (taxon_code ch1),(taxon_code ch2));
                code = taxon_code mine; } 
        in
        let node = { mine with adjusted=[node_dir]; } in
        (node,modified)

    (* calculate the median between a and b. old can be used as a heuristic,
     * branches are the supplied branch lengths of the children a and b, *)
    let median ?branches my_code old a b =
        let na, nb,code = 
            match my_code with
            | Some code ->
                    not_with code a.unadjusted, not_with code b.unadjusted,code
            | None ->
                    not_with (taxon_code b) a.unadjusted,
                    not_with (taxon_code a) b.unadjusted, -1
        in
        let old = match old with
            | Some oldness -> 
                (try
                    let n = with_both na.code nb.code oldness.unadjusted in
                    Some (n.lazy_node)
                 with
                  | Not_found -> None 
                )
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
        { unadjusted = [node]; adjusted = [node] }

    (* The final_states is calculated with the unadjusted component of the vertex *)
    let final_states grandcode par cur a b =
        let get_desired_dir par x = not_with (taxon_code par) x.unadjusted in
        let na = get_desired_dir cur a
        and nb = get_desired_dir cur b 
        and ncur = get_desired_dir par cur 
        and npar = 
            match grandcode with
            | Some x -> not_with x par.unadjusted
            | None ->
                    match par.unadjusted with
                    | [par] -> par
                    | _ -> failwith "AllDirNode.final_states"
        in
        let node = {
            lazy_node = OneDirF.final_states None 
                        npar.lazy_node ncur.lazy_node na.lazy_node nb.lazy_node;
            dir = Some (na.code, nb.code);
            code = taxon_code cur;
        } in
        match cur.unadjusted with
        | [_] -> { cur with unadjusted = [node] }
        | _ ->  let x, y = yes_with (taxon_code par) cur.unadjusted in
                { cur with unadjusted = [x; y; node] }

    let estimate_time left right = 
        let get_dir p c = (not_with (taxon_code p) c.unadjusted).lazy_node in
        OneDirF.estimate_time (get_dir left right) (get_dir right left)

    (** [apply_time child parent] applies time from parent into child --used on
     * leaves when the uppass_heuristic doesn't normally run over internal nodes
     *
     * first argument is for when time is from root median *)
    let apply_time root child parent =
        if uppass_debug then
            Printf.printf "Applying time from %d to %d"
                    (taxon_code parent) (taxon_code child);
        let has_one code x = match x.dir with
            | None -> true (* cannot be a leaf as well, no time! use median *)
            | Some (a,b) -> a = code || b = code
        in

        let c_data = match child.adjusted with
            | [x] -> assert (match x.dir with | None -> true | Some _ -> false);
                     x.lazy_node
            |  _  -> failwith "AllDirNode.apply_time only apply_time on leaves"
        and p_data = match parent.adjusted with
            | [x] -> assert(has_one (taxon_code child) x);
                     x.lazy_node
            |  _  -> let x,y = yes_with (taxon_code child) parent.adjusted in
                     x.lazy_node
        in

        let node = [{
            lazy_node = OneDirF.apply_time root c_data p_data;
            dir = None;
            code = taxon_code child;
        }] in
        { unadjusted = node; adjusted = node; }

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
     *    [A]==[M] < node_A
     *          \\
     *          [P]
    **)
    let uppass_heuristic p_code m a b p =
        if uppass_debug then
            Printf.printf "Performing uppass Heuristic on %d with (%d,%d) and %d/%d\n%!"
                (taxon_code m) (taxon_code a) (taxon_code b) (taxon_code p) p_code;
        let has_one code x = match x.dir with
            | None -> true
            | Some (a,b) -> a = code || b = code
        in

        let mc = taxon_code m and ac = taxon_code a
        and bc = taxon_code b 
        and get_dir parc x =
            try
                force_val (not_with parc x.adjusted).lazy_node
            with | Not_found -> 
                Printf.printf "\n-----\nCannot find direction %d in %d\n\t%!" 
                                parc (taxon_code x);
                q_print x;
                Printf.printf "-----\n\n%!";
                assert(false)
        in

        let data_m2p = match m.adjusted with
            (* then it hasn't been resolved by an earlier uppass on root *)
            | [x] -> 
                assert ( match x.dir with
                         | None -> false
                         | Some (xa,xb) -> (xa=ac && xb=bc) || (xa=bc && xb=ac)
                       );
                force_val x.lazy_node
            (* ...was resolved, so get the direction *)
            |  _  -> get_dir p_code m

        and data_p2m = match p.adjusted with
            (* then it HASN'T been resolved by an earlier uppass, and is truely
             * the parent with a calculated/lazy median *)
            | [x] -> if has_one mc x then
                        force_val x.lazy_node
                     else begin
                        Printf.printf "Cannot find %d in %d(%d)\n%!"
                            mc (taxon_code p) p_code;
                        q_print p;
                        q_print m;
                        print_newline ();
                        assert (false)
                     end
            (* ...was resolved, as above *)
            |  _  -> get_dir mc p

        (* both of these should exist *)
        and data_b2m = get_dir mc b
        and data_a2m = get_dir mc a in

        (* get the times in all directions --AB have them in M and M has it in P *)
        let time_M2A = Node.get_times_between data_m2p (Some (get_dir mc a))
        and time_M2B = Node.get_times_between data_m2p (Some (get_dir mc b))
        and time_M2P = match m.adjusted with
            | [x] -> Node.get_times_between data_p2m (Some data_m2p)
            |  _  -> Node.get_times_between data_p2m None
        in

        (* call medians with times supplied *)
        let node_A = lazy_from_fun
                (fun () -> Node.median_w_times 
                                (Some mc) (Some data_m2p) data_p2m
                                    data_b2m time_M2P time_M2B )
        and node_B = lazy_from_fun
                (fun () -> Node.median_w_times 
                                (Some mc) (Some data_m2p) data_p2m
                                    data_a2m time_M2P time_M2A ) in

        let dir_A= { code= mc; lazy_node= node_A; dir= Some(bc,p_code); }
        and dir_B= { code= mc; lazy_node= node_B; dir= Some(ac,p_code); }
        and dir_C= { code= mc; lazy_node= lazy_from_val data_m2p; dir= Some(ac,bc); } in
        let allDir = [ dir_A ; dir_B ; dir_C ] in
        { unadjusted = allDir; adjusted = allDir }

    let to_string nodes =
        let res =
            List.map (fun x -> OneDirF.to_string x.lazy_node)
            nodes.unadjusted
        in
        String.concat "\n" res

    let total_cost par n =
        match par with
        | Some code ->
                OneDirF.total_cost par (not_with code n.unadjusted).lazy_node
        | None ->
                match n.unadjusted with
                | [x] -> OneDirF.total_cost par x.lazy_node
                | [] -> failwith "The emtpy median? AllDirNode.total_cost"
                | _ -> failwith "AllDirNode.total_cost"

    let node_cost par n =
        match par with
        | Some code ->
                OneDirF.node_cost par (not_with code n.unadjusted).lazy_node
        | None ->
                match n.unadjusted with
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
                let da = 
                    match a.adjusted with
                    | [x] -> x
                    | _ -> failwith "AllDirNode.is_collapsable 1"
                and db = 
                    match b.adjusted with
                    | [x] -> x
                    | _ -> failwith "AllDirNode.is_collapsable 1"
                in
                OneDirF.is_collapsable `Dynamic da.lazy_node db.lazy_node
        | `Any ->
                (is_collapsable `Static a b) && (is_collapsable `Dynamic a b)

    let to_xml _ _ _ = ()

    let run_any f n =
        match n with
        | h :: _ -> lazy_from_val (f h.lazy_node)
        | [] -> failwith "AllDirNode.run_any"

    let run_all f n =
        let processor x = 
            let res = f x.lazy_node in
            { x with lazy_node = res }
        in
        { unadjusted = List.map processor n.unadjusted;
        adjusted = List.map processor n.adjusted }

    let get_something f code n =
        let node = 
            match code with
            | Some code ->
                    not_with code n
            | None ->
                    match n with
                    | [x] -> x
                    | [] -> failwith "AllDirNode.num_height 1"
                    | _ -> failwith "AllDirNode.num_height 2"
        in
        f code node.lazy_node

    let num_height code n = 
        get_something OneDirF.num_height code n.unadjusted

    let num_otus code n =
        get_something OneDirF.num_otus code n.unadjusted

    let get_sequences code n = 
        get_something OneDirF.get_sequences code n.adjusted

    let get_dynamic_preliminary code n =
        get_something OneDirF.get_dynamic_preliminary code n.unadjusted

    let get_dynamic_adjusted code n =
        get_something OneDirF.get_dynamic_preliminary code n.adjusted

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

    let get_others err code n =
            match code with
            | Some code -> Some (yes_with code n)
            | None -> None

    let get_node err code n =
            match code with
            | Some code -> not_with code n
            | None ->
                    match n with 
                    | [x] -> x
                    | _ -> failwith err

    let f_codes ints n = 
        run_all (OneDirF.f_codes ints) n

    let min_child_code code n = 
        get_something OneDirF.min_child_code code n.unadjusted

    type nad8 = OneDirF.nad8

    let new_characters = Node.Standard.new_characters

    let build_node static_characters chars node =
        let processor x = 
            let res = OneDirF.build_node static_characters chars x.lazy_node in
            { x with lazy_node = res }
        in
        let uadj = List.map processor node.unadjusted
        and adj = List.map processor node.adjusted in
        { unadjusted = uadj; adjusted = adj; }

    let set_exclude_info e x =
        let processor x = 
            let res = OneDirF.set_exclude_info e x.lazy_node in
            { x with lazy_node = res }
        in
        let uadj = List.map processor x.unadjusted
        and adj = List.map processor x.adjusted in
        { unadjusted = uadj; adjusted = adj }

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
            let uadj = List.map processor n.unadjusted
            and adj = List.map processor n.adjusted in
            { unadjusted = uadj; adjusted = adj }
    end

    module Union = struct

        type u = OneDirF.Union.u
        let union code n a b =
            let node = get_node "AllDirNode.AllDirF.Union.union" code
            n.unadjusted in
            OneDirF.Union.union code node.lazy_node a b

        let union_final code a b = 
            let b = get_node "AllDirNode.AllDirF.Union.union_final" code
            b.unadjusted in
            OneDirF.Union.union_final code a b.lazy_node

        let union_preliminary code a b = 
            let b = 
                get_node "AllDirNode.AllDirF.Union.union_preliminary" code
                b.unadjusted
            in
            OneDirF.Union.union_preliminary code a b.lazy_node

        let leaf taxon_code code n =
            let x = get_node "AllDirNode.AllDirF.Union.leaf" code n.unadjusted in
            OneDirF.Union.leaf taxon_code code x.lazy_node

        let distance a b = OneDirF.Union.distance a b

        let saturation = OneDirF.Union.saturation

        let distance_node code n u =
            let node = 
                get_node "AllDirNode.AllDirF.Union.distance_node" code
                n.unadjusted
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
        List.map (fun x -> { unadjusted = x; adjusted = x }) res

    let root_cost a = 
        let lst =
            match !Methods.cost with
            | `Iterative _ -> a.adjusted
            | _ -> a.unadjusted
        in
        match lst with
        | [a] -> OneDirF.root_cost a.lazy_node
        | _ -> failwith "AllDirNode.root_cost"

    let tree_cost a b = (total_cost a b) +. (root_cost b)

end

type 'a node_hybrid = {
    st : Node.Standard.n option;
    dy : 'a;
}

module HybridF = struct
    let get_dynamic x = x.dy
end

let create_root ?branches a aa ab b ba bb opt =
    let middle = match opt with
        | Some x -> x 
        | None -> AllDirF.median ?branches None None a b
    in
    let a_final = match aa,ab with
        | Some aa,Some ab ->
                let a = AllDirF.uppass_heuristic (AllDirF.taxon_code b) a aa ab middle in
                assert( (List.length a.adjusted) = 3);
                a
        | None,None ->
                AllDirF.apply_time true a middle
        | _,_ -> failwith "Failed with children of A @ Create Roots"
    and b_final = match ba,bb with
        | Some ba,Some bb -> 
                let b = AllDirF.uppass_heuristic (AllDirF.taxon_code a) b ba bb middle in
                assert( (List.length b.adjusted) = 3);
                b
        | None,None -> 
                AllDirF.apply_time true b middle
        | _,_ -> failwith "Failed with children of B @ Create Roots"
    and l_middle = match middle.adjusted with
        | [x] -> 
            assert( match x.dir with | None -> false
                    | Some (xa,xb) -> 
                        (xa = (AllDirF.taxon_code a) && xb = (AllDirF.taxon_code b)) ||
                        (xa = (AllDirF.taxon_code b) && xb = (AllDirF.taxon_code a))
                  );                            
            x.lazy_node
        |  _  -> failwith "median has multiple directions?"
    in
    a_final,b_final,l_middle

let create_root_w_times (adjusted:bool) left right =
    let get_dir parc x =
        try
            let refresh_from = if adjusted then x.adjusted else x.unadjusted in
            force_val (not_with parc refresh_from).lazy_node
        with | Not_found -> failwithf "rootw/times: Cannot find %d in %d"
                                            parc (AllDirF.taxon_code x)

    and get_a_dir child x = 
        let refresh_from = if adjusted then x.adjusted else x.unadjusted in
        match refresh_from with (* one direction now means it's a leaf *)
        | [x] -> assert( match x.dir with | None -> true | _ -> false); 
                 force_val x.lazy_node
        |  _  -> let x,y = yes_with child refresh_from in 
                 force_val x.lazy_node
    in
        
    let left2right = get_dir (AllDirF.taxon_code right) left
    and right2left = get_dir (AllDirF.taxon_code left) right in

    let in_l_time = Node.get_times_between 
                    (get_a_dir (AllDirF.taxon_code right) left)
                    (Some right2left)
    and in_r_time = Node.get_times_between 
                    (get_a_dir (AllDirF.taxon_code left) right)
                    (Some left2right) in
    lazy_from_fun
        (fun () -> 
            Node.median_w_times 
                None None left2right right2left in_l_time in_r_time)

(** [verify_branch_lengths a b c m] verify an internal node [m] *)
let verify_branch_lengths a b c m : bool =
    let taxon_code = AllDirF.taxon_code in
    let get_dir parc x =
        try
            force_val (not_with parc x.adjusted).lazy_node
        with | Not_found ->
            q_print x;
            failwithf "Cannot find %d in %d\n%!" parc (taxon_code x)
    in

    (* verify x with internal node *)
    let verify internal x =
        let ic,xc = taxon_code internal, taxon_code x in
        let m_x1,m_x2 = yes_with xc internal.adjusted in
        match List.filter (fun y -> match y.dir with 
                                    | None -> true
                                    | Some (a,b) -> a = ic || b = ic)
                          x.adjusted with
            | [a;b] -> (* x is also an internal node *)
                (Node.verify_time (force_val m_x1.lazy_node) (get_dir xc internal)
                                  (force_val a.lazy_node) (Some (get_dir ic x))) &&
                (Node.verify_time (force_val m_x2.lazy_node) (get_dir xc internal)
                                  (force_val b.lazy_node) (Some (get_dir ic x)))
            | [a]   ->
                assert (match a.dir with | None -> true | _ -> false );
                (Node.verify_time (force_val m_x1.lazy_node) (get_dir xc internal)
                                  (force_val a.lazy_node) None) &&
                (Node.verify_time (force_val m_x2.lazy_node) (get_dir xc internal)
                                  (force_val a.lazy_node) None)
            |  _    -> Printf.printf "%d doesn't have %d at all!\n%!" xc ic;
                       q_print x;
                       q_print internal;
                       failwith "Verification Failed"
    in
    (* verify each direction *)
    let adjusted_passed = 
        List.fold_right (fun x acc -> acc && (verify m x)) [a;b;c] true
    in
    adjusted_passed
