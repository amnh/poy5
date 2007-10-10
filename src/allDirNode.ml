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

type exclude = Node.exclude

type direction = (int * int) option

type a_node = Node.node_data Lazy.t

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

let my_own_code () = incr Data.median_code_count; !Data.median_code_count

let to_n node = Lazy.lazy_from_val node

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
        | Some (a, b) ->
                string_of_int a ^ ", " ^ string_of_int b
    in
    Status.user_message Status.Error msg

let yes_with code n =
    assert (0 <> List.length n);
    match List.filter (fun x -> not (has_code code x)) n with
    | [x; y] -> x, y
    | x -> 
            List.iter print_pairs x;
            failwith 
            ("AllDirNode.not_with " ^ string_of_int code ^ 
            " of " ^ string_of_int (get_code x) ^ " has " ^
            string_of_int (List.length x) ^ " matching values.")

let not_with code n =
    assert (0 <> List.length n);
    match List.filter (has_code code) n with
    | [x] -> x
    | x -> 
            List.iter print_pairs x;
            failwith 
            ("AllDirNode.not_with " ^ string_of_int code ^ 
            " of " ^ string_of_int (get_code x) ^ " has " ^
            string_of_int (List.length x) ^ " matching values.")

module OneDirF : 
    NodeSig.S with type e = exclude with type n = a_node with type
    nad8 = Node.Standard.nad8 = struct

    type n = a_node
    type e = exclude

    let load_data ?taxa ?codes ?(classify=true) data = 
        let data, nodes = 
            match taxa, codes with
            | None, None -> Node.Standard.load_data ~classify data 
            | Some v, None -> Node.Standard.load_data ~taxa:v ~classify data
            | None, Some v -> Node.Standard.load_data ~codes:v ~classify data
            | Some v, Some w -> 
                    Node.Standard.load_data ~taxa:v ~codes:w ~classify data
        in 
        data, List.map to_n nodes

    let fix_preliminary x = x

    let apply_f_on_lazy f a b = 
        let a = Lazy.force_val a
        and b = Lazy.force_val b in
        f a b

    let distance ?(para=None) ?(parb=None) a b =
        apply_f_on_lazy Node.Standard.distance a b

    let median code my_code old a b = 
        Lazy.lazy_from_fun 
        (fun () -> apply_f_on_lazy 
        (Node.Standard.median None my_code None) a b)

    let median_3 x par cur a b = 
        Lazy.lazy_from_val
        (apply_f_on_lazy
        (Node.Standard.median_3 x (Lazy.force_val par) (Lazy.force_val cur))
        a b)

    let to_string v = Node.Standard.to_string (Lazy.force_val v)

    let apply_single_f_on_lazy f a = 
        f (Lazy.force_val a)

    let total_cost x n =
        apply_single_f_on_lazy (Node.Standard.total_cost x) n

    let node_cost par n = 
        apply_single_f_on_lazy (Node.Standard.node_cost par) n

    let update_leaf x = x

    let taxon_code n = 
        apply_single_f_on_lazy Node.Standard.taxon_code n

    let union_distance _ _ = 0.0

    let is_collapsable a b = 
        apply_f_on_lazy Node.Standard.is_collapsable a b

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
        Lazy.lazy_from_val (apply_single_f_on_lazy (Node.Standard.f_codes ints) n)

    let min_child_code code n = 
        apply_single_f_on_lazy (Node.Standard.min_child_code code) n

    type nad8 = Node.Standard.nad8

    let new_characters = Node.Standard.new_characters

    let build_node static_characters chars n = 
        Lazy.lazy_from_val 
        (apply_single_f_on_lazy 
        (Node.Standard.build_node static_characters chars) 
        n)

    let set_exclude_info e n = 
        Lazy.lazy_from_val 
        (apply_single_f_on_lazy (Node.Standard.set_exclude_info e) n)

    let excludes_median code a b = 
        apply_f_on_lazy (Node.Standard.excludes_median code) a b

    let has_excluded = Node.Standard.has_excluded

    module T = struct

        let add_exclude set = 
            fun x ->
                Lazy.lazy_from_val (Node.Standard.T.add_exclude set 
                (Lazy.force_val x))
    end

    module Union = struct
        type u = Node.Standard.Union.u Lazy.t

        let union code a b c =
            let a = Lazy.force_val a 
            and b = Lazy.force_val b 
            and c = Lazy.force_val c in
            Lazy.lazy_from_val (Node.Standard.Union.union code a b c)

        let union_one which code a b =
            let a = Lazy.force_val a
            and b = Lazy.force_val b in
            Lazy.lazy_from_val (which code a b)

        let union_preliminary a b c =
            union_one Node.Standard.Union.union_preliminary a b c

        let union_final a b c =
            union_one Node.Standard.Union.union_final a b c

        let leaf taxon_code code a = 
            Lazy.lazy_from_val 
                (Node.Standard.Union.leaf taxon_code code (Lazy.force_val a))

        let distance a b = 
            Node.Standard.Union.distance 
            (Lazy.force_val a)
            (Lazy.force_val b)

        let saturation a = 
            Node.Standard.Union.saturation (Lazy.force_val a)

        let distance_node code n u =
            Node.Standard.Union.distance_node code
            (Lazy.force_val n)
            (Lazy.force_val u)

        let compare a b =
            Node.Standard.Union.compare (Lazy.force_val a) (Lazy.force_val b)

        let get_sequence a b c = 
            Node.Standard.Union.get_sequence a b (Lazy.force_val c)
    end

    let for_support a b c d : n list =
        let nb = 
            List.map (fun (a, (b : n)) ->
                a, Lazy.force_val b) b
        in
        let res = Node.Standard.for_support a nb c d in
        List.map Lazy.lazy_from_val res


    let root_cost a = 
        Node.Standard.root_cost (Lazy.force_val a)

    let to_single root a b c d =       
        let root' = match root with
        | Some root -> Some (Lazy.force_val root) 
        | None -> None
        in
        let b' = (Lazy.force_val b) in
        let d' = (Lazy.force_val d) in
        lazy (Node.Standard.to_single root' a b' c d')
end

module AllDirF : NodeSig.S with type e = exclude with type n = node_data with
type nad8 = Node.Standard.nad8 = struct

    type n = node_data
    type e = exclude 

    let to_n node = 
        let node_dir = {
            lazy_node = to_n node;
            dir = None;
            code = node.Node.taxon_code;
        }
        in
        { unadjusted = [node_dir]; adjusted = [node_dir]}

    let load_data ?taxa ?codes ?(classify=true) data = 
        let data, nodes = 
            match taxa, codes with
            | None, None -> Node.Standard.load_data ~classify data 
            | Some v, None -> Node.Standard.load_data ~taxa:v ~classify data
            | None, Some v -> Node.Standard.load_data ~codes:v ~classify data
            | Some v, Some w -> 
                    Node.Standard.load_data ~taxa:v ~codes:w ~classify data
        in 
        data, List.map to_n nodes

    let fix_preliminary x = x

    let to_single root a b c d =
        let get_code x = 
            match x with
            | Some x -> x
            | None -> failwith "huh?"
        in
        let a' = get_code a
        and c' = get_code c in
        let b' = not_with a' b.unadjusted
        and d' = not_with c' d.unadjusted in

        let root = match root with 
        | Some root ->
              let root' = List.hd root.unadjusted in
              Some root'.lazy_node
        | None -> None
        in  

        let lazy_node = OneDirF.to_single root a b'.lazy_node c d'.lazy_node in
        let node = { d' with lazy_node = lazy_node } in
        { d with unadjusted = [node] }

    let distance ?(para=None) ?(parb=None) a b =
        match a.unadjusted, b.unadjusted with
        | [a], [b] -> OneDirF.distance a.lazy_node b.lazy_node
        | a, [b] -> 
                (match parb with
                | None -> failwith "AlldNode.distance 1"
                | Some para -> 
                        let a = not_with para a in
                        OneDirF.distance a.lazy_node b.lazy_node)
        | [a], b -> 
                (match para with
                | None -> failwith "AlldNode.distance 2"
                | Some parb -> 
                        let b = not_with parb b in
                        OneDirF.distance a.lazy_node b.lazy_node)
        | a, b -> 
                (match para, parb with
                | Some para, Some parb -> 
                        let a = not_with para a 
                        and b = not_with parb b in
                        OneDirF.distance a.lazy_node b.lazy_node
                | _ -> failwith "AlldNode.distance 3")

    let median code my_code old a b =
        let my_code =
            match my_code with
            | Some code -> code
            | None -> my_own_code ()
        in
        let na, nb = 
            match code with
            | Some code ->
                    not_with code a.unadjusted, not_with code b.unadjusted
            | None ->
                    match a.unadjusted, b.unadjusted with
                    | [a], [b] -> a, b
                    | _, _ -> failwith "AlldNode.median"
        in
        let lazy_node = 
            OneDirF.median None (Some my_code) None na.lazy_node nb.lazy_node
        in
        let node = {
            lazy_node = lazy_node;
            dir = Some (na.code, nb.code); 
            code = my_code;
        }
        in
        { unadjusted = [node]; adjusted = [node] }

    let taxon_code n = 
        match n.unadjusted with
        | h :: _ -> (** All the elements in the list have the same code *)
                h.code
        | [] -> failwith "AllDirNode.taxon_code"

    (* The median_3 is calculated with the unadjusted component of the vertex *)
    let median_3 grandcode par cur a b =
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
                    | _ -> failwith "AllDirNode.median_3"
        in
        let node = {
            lazy_node = 
                OneDirF.median_3 None npar.lazy_node ncur.lazy_node
                na.lazy_node nb.lazy_node;
            dir = Some (na.code, nb.code);
            code = taxon_code cur;
        }
        in
        match cur.unadjusted with
        | [_] -> { cur with unadjusted = [node] }
        | _ ->
                let x, y = yes_with (taxon_code par) cur.unadjusted in
                { cur with unadjusted = [x; y; node] }

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

    let is_collapsable a b =
        let acode = taxon_code a
        and bcode = taxon_code b in
        let da = not_with bcode a.unadjusted
        and db = not_with acode b.unadjusted in
        OneDirF.is_collapsable da.lazy_node db.lazy_node

    let to_xml _ _ _ = ()

    let run_any f n =
        match n with
        | h :: _ ->
                Lazy.lazy_from_val (f h.lazy_node)
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
        get_something OneDirF.get_sequences code n.unadjusted

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
        Lazy.force_val (run_any (OneDirF.n_chars ~acc) n.unadjusted)

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

    let get_nodes code a b =
        match code with
        | Some code ->
                not_with code a.unadjusted, not_with code b.unadjusted 
        | None ->
                match a.unadjusted, b.unadjusted with
                | [a], [b] -> a, b
                | _, _ -> failwith "AllDirNode.excludes_median"

    let excludes_median code a b = 
        let nodea, nodeb = get_nodes code a b in
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
                | [b] -> a, Lazy.force_val b.lazy_node
                | _ -> failwith "AllDirNode.for_support") b 
        in
        let res = Node.Standard.for_support a nb c d in
        let rec merge a b =
            match a, b with
            | (_, [ha]) :: ta, hb :: tb ->
                    [{ ha with lazy_node = Lazy.lazy_from_val hb }] ::
                        merge ta tb
            | [], [] -> []
            | _, _ -> failwith "AllDirNode.for_support 2"
        in
        let b = List.map (fun (a, x) -> a, x.unadjusted) b in
        let res = merge b res in
        List.map (fun x -> { unadjusted = x; adjusted = x }) res

    let root_cost a = 
        match a.unadjusted with
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
