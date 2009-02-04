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

let eager_lazy = false

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
    if eager_lazy then
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

let yes_with code n =
    assert (0 <> List.length n);
    match List.filter (fun x -> not (has_code code x)) n with
    | [x; y] -> x, y
    | x -> List.iter print_pairs x;
            failwithf "AllDirNode.yes_with %d of %d has %d matching values."
                        code (get_code n) (List.length x)
let not_with code n =
    assert (0 <> List.length n);
    match List.filter (has_code code) n with
    | [x] -> x
    | x -> List.iter print_pairs x;
            failwithf "AllDirNode.not_with %d of %d has %d matching values."
                        code (get_code n) (List.length x)

let with_both c1 c2 n = (* grabs a node in the direction not in c1 or c2 *)
    assert( 0 <> List.length n );
    let lil_filter c = List.filter (fun x -> not (has_code c x)) in
    match lil_filter c1 (lil_filter c2 n) with
    | [x] -> x
    |  x  -> raise Not_found

let either_with c n = (* grabs a node that has this branch c--n, in n *)
    assert( 0 <> List.length n);
    match List.filter (fun x -> not (has_code c x)) n with
        | hd::tl -> hd
        | [] -> (match n with
            | [x] ->
                Printf.printf "Singular with (%s)\n%!"
                        (match x.dir with
                        | None -> "None"
                        | Some (a, b) -> Printf.ksprintf (fun x -> x) "%d, %d" a b);
                    x
            |  _  -> List.iter print_pairs n;
            failwithf ("AllDirNode.either_with %d of %d has 0 matching"^^
                        " values in %d total.") c (get_code n) (List.length n))

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

    let apply_time x y =
        lazy_from_fun (fun () -> apply_f_on_lazy Node.Standard.apply_time x y)

    let estimate_time x y =
        Node.Standard.estimate_time (force_val x) (force_val y)

    let uppass_heuristic ?branches = final_states

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

    let to_single root a b c d set =
        let root' = force_opt root
        and b' = (force_val b)
        and d' = (force_val d) in
        lazy_from_val (Node.Standard.to_single root' a b' c d' set)

    let force x = force x
end

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
        let node = { d' with lazy_node = lazy_node } in
        { d with unadjusted = [node] }

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

    let edge_iterator par mine ch1 ch2 = 
        let get_dir p n = (not_with (taxon_code p) n.unadjusted).lazy_node in
        let atom, btom, parofm = match par with
                    | Some x -> 
                        get_dir mine ch1, get_dir mine ch2, Some (get_dir mine x)
                    | None -> 
                        get_dir ch2 ch1, get_dir ch1 ch2, None
        and ab_to_m = (with_both (taxon_code ch1) (taxon_code ch2) mine.unadjusted) in
        
        let node_dir  = 
            { ab_to_m with
                lazy_node = OneDirF.edge_iterator parofm ab_to_m.lazy_node atom btom
            } in
        { unadjusted = [node_dir]; adjusted = [node_dir]; }

    let readjust gp_opt mode to_adjust ch1 ch2 par mine =
        (* in [n], we want the direction toward [p], the parent *)
        let get_dir p_code n = (not_with (taxon_code p_code) n.unadjusted).lazy_node in
        let mine2par = match gp_opt with
            | Some gpc -> (not_with gpc par.unadjusted).lazy_node
            | None -> (match par.adjusted with
                        | [x] -> x.lazy_node
                        |  _  -> failwithf "AllDirNode.readjust got no GP but with %d directions."
                                                (List.length par.unadjusted)
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

    let median ?branches my_code old a b =
        let na, nb,code = 
            match my_code with
            | Some code ->
                    not_with code a.unadjusted, not_with code b.unadjusted,code
            | None ->
                    not_with (taxon_code b) a.adjusted,
                    not_with (taxon_code a) b.adjusted, -1
        in
        let old = match old with
            | Some oldness -> 
                (try let n = with_both na.code nb.code oldness.unadjusted in
                     Some (n.lazy_node) with
                 | Not_found -> None )
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
     * leaves when the uppass_heuristic doesn't normally run over *)
    let apply_time child parent =
        let get_dir p xun = (not_with (taxon_code p) xun).lazy_node in
        let has_with code n = match n.dir with
            | None -> true 
            | Some (a,b) -> a = code || b = code in
        let node = match
            List.filter (has_with (taxon_code child)) parent.unadjusted with
            | [hd] 
            | [hd;_] -> OneDirF.apply_time (get_dir parent child.unadjusted)
                                           hd.lazy_node
            | [] -> 
                OneDirF.apply_time (get_dir parent child.unadjusted)
                                   (get_dir child parent.unadjusted)
            | _  -> failwith "AllDirNode.apply_time: Too many directions to find time"
        in
        let node = [{
            lazy_node= node;
            dir = None;
            code = taxon_code child;
        }] in
        { unadjusted = node; adjusted = node; }

    (** [uppass_heuristic tbl curr ch1 ch2 par]
     * [cur], at this point will have one direction ([ch1],[ch2]), after this call
     * it will fill in the other directions, lazily, with [par]. This process will
     * update the nodes, the times (likelihood), et cetera, with lazy values.
     *            cb'
     *   ca'       |
     *    \      [cB]---cb'
     *     \     //
     *    [cA]==[M] < dir_A
     *     /     \\
     *    /      [pC]--gp?
     *   ca'       |
     *            gp?
    **)
    let uppass_heuristic ?branches _ m a b c =
        let get_dir parc x = force_val (not_with parc x.unadjusted).lazy_node in
        let ac,bc,cc,mc = taxon_code a, taxon_code b, taxon_code c, taxon_code m in
        let data_c2m = get_dir cc m in
        let time_M2C = 
            match List.filter (fun x -> not (has_code mc x)) c.unadjusted with
                | h :: _ ->
                    Node.get_times_between (force_val h.lazy_node) data_c2m
                | [] -> match c.unadjusted with
                    | [single] ->
                        (match branches with
                        | Some tbl -> 
                            Node.get_times_between_tbl 
                                        (Hashtbl.find tbl cc)
                                        (force_val single.lazy_node)
                        | None ->
                            Node.get_times_between (get_dir mc c) data_c2m
                        )
                    | _ -> failwith "AllDirNode: uppass_heuristic"
        and nd_c = get_dir mc c
        and time_M2A,time_M2B = (* ignore table since its used for root *)
                Node.get_times_between data_c2m (get_dir mc a),
                Node.get_times_between data_c2m (get_dir mc b)
        in

        let node_A = lazy_from_fun
                (fun () -> Node.median_w_times 
                                (Some mc) (Some data_c2m) nd_c
                                    (get_dir mc b) time_M2C time_M2B )
        and node_B = lazy_from_fun
                (fun () -> Node.median_w_times 
                                (Some mc) (Some data_c2m) nd_c
                                    (get_dir mc a) time_M2C time_M2A )
        in

        let dir_A= { code= mc; lazy_node= node_A; dir= Some(bc,cc); }
        and dir_B= { code= mc; lazy_node= node_B; dir= Some(ac,cc); }
        and dir_C= { code= mc; lazy_node= lazy_from_val data_c2m; dir= Some(ac,bc); } in

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
end

type 'a node_hybrid = {
    st : Node.Standard.n option;
    dy : 'a;
}

module HybridF = struct
    let get_dynamic x = x.dy
end

(** [verify_branch_lengths a b c m] verify an internal node [m] *)
let verify_branch_lengths a b c m : bool =
    info_user_message "Verifying Branch Lengths";
    let taxon_code = AllDirF.taxon_code in
    let get_dir parc x = force_val (not_with parc x.unadjusted).lazy_node in

    (* verify x with internal *)
    let verify internal x =
        let ic = taxon_code internal and xc = taxon_code x in
        let m_1,m_2 = yes_with xc internal.unadjusted
        and m_oth = get_dir xc internal
        and x_s =
            let x_oth = get_dir ic x in
            match List.filter (fun y -> not (has_code ic y))
                              x.unadjusted with
            | [one;two] -> [(one,Some x_oth);(two,Some x_oth)]
            | [] -> [(not_with ic x.unadjusted,None)]
            | err -> failwithf "Verify -- %d directions found" (List.length err)
        and string_pairs r = match r.dir with
            | None -> "None"
            | Some (a,b) -> Printf.ksprintf (fun x -> x) "%d,%d" a b
        in

        List.fold_right
            (fun (x,x_oth) ac1 ->
                List.fold_right
                    (fun m ac2 ->
                        let res = Node.verify_time (force_val m.lazy_node) m_oth
                                                   (force_val x.lazy_node) x_oth
                        in
                        if res then
                            info_user_message "Verified (%s) and (%s) between %d and %d"
                                (string_pairs m) (string_pairs x) m.code x.code
                        else
                            error_user_message "Failed (%s) and (%s) between %d and %d"
                                (string_pairs m) (string_pairs x) m.code x.code;
                        res && ac2
                    ) [m_1;m_2] ac1
            ) x_s true
    in
    (* verify each direction *)
    List.fold_right (fun x acc -> acc && (verify m x)) [a;b;c] true
