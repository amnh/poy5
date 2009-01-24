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

(* The Sankoff characters have some particular properties in the tree cost
 * calculation.  Each HTU has a set of possible states, in principle, any allowed
 * state in the character could be a state in the HTU. As the tree evaluation moves
 * down in the downpass, it continues to test all the possible combinations, so
 * potentially, down in the tree, a suboptimal local solution becomes optimal. For
 * this reason, it is necessary to keep the cost of a particular character if
 * assigned on each node, and this cost corresponds to the transformation between
 * it and the optimal elements in the subtree that it roots. In reality we never
 * handle unrooted trees for this kind of operations (remember the tree module has
 * a handle for "Unrooted" trees, meaning that we can safely keep this meaning
 * properly. *)
let () = SadmanOutput.register "SankCS" "$Revision: 2871 $"


let infinity = max_int / 4
let is_infinity x = x >= infinity

type cost = int

let cost_less a b = 
    if is_infinity a then false
    else if is_infinity b then true
    else a < b

let cost_min a b =
    if cost_less a b
    then a else b

let cost_plus a b =
    if is_infinity a || is_infinity b then infinity
    else a + b

let cost_minus a b =
    if is_infinity a then
        let () = assert (not (is_infinity b)) in
        max_int
    else if is_infinity b then
        let () = assert (not (is_infinity a)) in
        0
    else a - b

let float_of_cost x = 
    if is_infinity x then Pervasives.infinity else float_of_int x

let string_of_cost x = 
    if is_infinity x then "inf" else string_of_int x

let ( +$ ) a b = cost_plus a b
let ( -$ ) a b = cost_minus a b
let ( <$ ) a b = cost_less a b

type cm = int array array               (* never have infinity *)

type gen = int * int * cm

let color = Character.Black

type elt = {
    ecode : int;
    
    (* The array of states, each index in the array corresponds to a single
     * state and the integer in the array corresponds to the cost of the that
     * state given the neighbors of the character in the node that contains it in
     * the tree.
     *)
    s : cost array;

    (* The following values come from Goloboff 1998: Tree Searches under Sankoff
     * Parsimony. *)
    
    (* This is a downpass-calculated intermediate value *)
    beta : cost array;

    (* e.(s) is the preliminary added cost for using state s
     * (uppass-calculated) *)
    e : cost array;

    (* M.(s) is a cached value; it depends on E and beta, so we reset it during
     * the downpass. *)
    mutable m : cost array option;

    (* Scratch area for the diagnosis, it is safe to do funcky things with it
    * because we only use it during the diagnosis *)
    best_states : int list;

}
(* The Sankoff character type *)
type t = { 
    (* The code the character belongs to. Every pair of characters with the same
     * code correspond to homologous units *)
    code : int;
    (* The transformation cost matrix to be used for this character (and its
     * homologous characters) *)
    tcm : cm;
    (* The color of the character being used *)
    color : Character.c;

    (* A bunch of individual characters *)
    elts : elt array;
}

(* A default string representation as a list with the minimal states (those with
* minimal cost in the state array). *)
let elt_to_string a = 
(*     let min = Array.fold_left cost_min Infinity a.s in *)
    let sep = ref "" in
    let res, _ = Array.fold_left 
            begin fun (b, pos) a -> 
                let next =
                    if a = infinity
                    then b
                    else
                        let next =
                    b ^ !sep
                    ^ string_of_int pos ^ "="
                    ^ string_of_cost a
                        in
                        sep := ";";
                        next
                in
(*                     if a = min then begin *)
(*                         let res = b ^ !sep ^ string_of_int pos in *)
(*                         sep := ";"; *)
(*                         res *)
(*                     end *)
(*                     else b *)
(*                 in *)
                next, pos + 1
            end ("{",0) a.s in
    res ^ "}"

let to_string s =
    let list = Array.to_list s.elts in
    let strings = List.map elt_to_string list in
    "[" ^ String.concat "; " strings ^ "]"

let to_list s = 
    let list = Array.map (fun x -> x.ecode, 
        Array.map (fun x ->
            if is_infinity x then max_int
            else x) x.s) s.elts in
    Array.to_list list

let elt_to_full_string a =
    let string_of_costarray a =
        let ar =
            Array.mapi (fun i c -> (string_of_int i ^ "=" ^ string_of_cost c)) a
        in String.concat ";" (Array.to_list ar) in
    let states = string_of_costarray a.s in
    let beta = string_of_costarray a.beta in
    let e = string_of_costarray a.e in
    "states=[" ^ states ^ "]; beta=[" ^ beta ^ "]; e=[" ^ e ^ "]"
let to_full_string s =
    let list = Array.to_list s.elts in
    let strings = List.map elt_to_full_string list in
    "[" ^ String.concat "; " strings ^ "]"

let to_string = to_full_string
    

let assert_ninf a x y=
    assert (
        if is_infinity a then begin
         print_string (elt_to_full_string x);
            print_string (elt_to_full_string y);
            false
        end else true);
    a



(* An empty character. This is used only for the parsing process. *)
let empty =
    { elts = [||];
      code = (-1);
      color = Character.Black;
      tcm = [|[||]|];
    }

let empty_elt =
    { s = [||];
      beta = [||];
      e = [||];
      m = None;
      ecode = (-1);
      best_states = [];
    }

let codes {elts=elts} =
    let len = Array.length elts in
    let acc = ref [] in
    for i = len - 1 downto 0 do
        acc := (elts.(i).ecode) :: !acc
    done;
    !acc

(* The cost of a transformation occurring in the following order i -> j + i -> k
 * as defined in the transformation cost matrix par. *)
let median_cost par i j k = 
    (par.(i).(j) + par.(i).(k))

let get_min a = Array.fold_left cost_min infinity a


(* The code of the character a *)
let code a = a.code

(* Code of the element *)
let ecode e = e.ecode

(* The transformation cost matrix associated with character a *)
let tcm a = a.tcm

let set_tcm cm t  = {t with tcm = cm;}

let print_tcm tcm =
    Array.iter (fun row ->
                    Array.iter (fun v ->
                                    print_string (string_of_int v ^ " ")) row;
                    print_newline ())
        tcm

(* How many states there are *)
let nstates a = Array.length a.tcm

(* We'll be using this a lot ...
 * (make it monomorphic; more efficient?) *)
let store_min (a : cost) b = if a <$ !b then b := a

let is_inf x = is_infinity x


(* Even OTUs must have [beta] values, so we may as well store [e], also. *)
let canonize tcm a =
    let states = Array.length a.s in
    let minval = Array.fold_left cost_min infinity a.s in
    let e = Array.init states
        (fun i ->
            if is_infinity a.s.(i) then infinity
            else a.s.(i) -$ minval) in
    let beta = Array.init states
        (fun s ->
             let best = ref infinity in
             for x = states - 1 downto 0 do
                 store_min
                     (if is_inf e.(x) then e.(x)
                      else e.(x) +$ (tcm.(s).(x)))
                     best
             done;
             !best) in
    { a with e = e; beta = beta; }

(* A random character generator specification *)
let rand_gen () = 
    let states = 4 in
    (* Random.int only accepts arguments up to 2^30 ... *)
    let code = Random.int 1073741823
    and tcm = Array.make_matrix 4 4 ( 1) in
    let num = Random.int 15 in 
    num, code, states, tcm

(* The random state generator, using the random character generator created by
 * rand_gen () *)
(* let make_rand (num, code, states, tcm) = *)
(*     let states = Array.init states (fun _ -> Cost (1 + (Random.int states))) in *)
(*     canonize *)
(*         { empty with s = states; code = code; color = Character.Black; tcm = tcm } *)

let get_minstates {e=a} =
    let (_, list) = 
        Array.fold_right (fun v (num, lst) ->
                              if v = 0
                              then (num - 1, num :: lst)
                              else (num - 1, lst)) a (Array.length a - 1, [])
    in list


(* Calculates the median between to characters a and b. a and b have to be
* homologous (same code), and must also have the same transformation cost matrix
* associated. If they are homologous, it must also be the case that they hold
* the same number of valid states. *)
let elt_median tcm a b = 
    assert ((Array.length a.s) = (Array.length b.s));
    assert (a.ecode = b.ecode);
    let states = Array.length a.s in

    let min_cost = ref infinity in
    (* calculate the preliminary costs, and store our minimum value *)
    let states_init i = 
        let best = ref infinity in
        for j = states - 1 downto 0 do
            for k = states - 1 downto 0 do
                let combination_cost = median_cost tcm i j k in
                let combination_cost = combination_cost +$ a.s.(j) +$ b.s.(k) in
                store_min combination_cost best
            done;
        done;
        store_min !best min_cost;
        assert_ninf !best a b
    in
    let c = Array.init states states_init in

    (* preliminary added cost *)
    let e = Array.init states
        (fun s ->
             c.(s) -$ !min_cost) in

    (* beta value: see Goloboff 1998 *)
    let beta_init s =
        let best = ref infinity in
        for x = states - 1 downto 0 do
            store_min ((tcm.(s).(x)) +$ e.(x)) best
        done;
        !best
    in
    let beta = Array.init states beta_init in
    
    { a with s = c; beta = beta; e = e; m = None; }

let median _ a b =
    assert (a.code = b.code);
    assert (a.tcm = b.tcm);
    { a with elts = Array.mapi (fun i -> elt_median a.tcm a.elts.(i)) b.elts }

(* Calculates the distance between two characters a and b *)
(* Note that we only calculate the _added_ distance *)
let elt_distance tcm a b =
    let median = elt_median tcm a b in
    let get_cost m =
        float_of_cost (assert_ninf (Array.fold_left cost_min infinity m.s) a b) in
    let med_cost = get_cost median in
    let a_cost = get_cost a in
    let b_cost = get_cost b in
    med_cost -. a_cost -. b_cost

let distance a b =
    let tcm = a.tcm in
    let acc = ref 0. in
    for i = Array.length a.elts - 1 downto 0 do
        acc := !acc +. elt_distance tcm a.elts.(i) b.elts.(i)
    done;
    !acc

(* Compares two characters a and b. Note that this comparison is used basically
 * in the sets, and therefore, comparing the codes is enough (all the characters in
 * a set have different codes). *)
let compare_codes a b = 
    compare a.code b.code

(* Compares the data in the two characters a and b.  This is used to find cases
 * where a graph update does not change values far away in the graph, so that
 * propagation can be stopped. *)
let compare_data a b =
    compare a b

(* A dummy function to be improved later using the official parser *)
let parse a = [("", empty) ]
let init2 len1 len2 fn =
    Array.init len1
        (fun i ->
             let fn = fn i in
             Array.init len2 (fun j -> fn j))


(* This algorithm is taken from Goloboff 1998.  Calculate D', then the final E
 * value from that. *)
let elt_median_3 tcm a n l r =          (* ancestor, node, left, right *)
    let states = Array.length n.s in
    (* NOTE that this is REVERSED: d'[i][s], not d'[s][i] *)
    let init_d' i =
        let best = ref infinity in
        for x = states - 1 downto 0 do
            store_min ((tcm.(i).(x)) +$ l.beta.(x) +$ r.beta.(x)) best
        done;
        fun s ->
            a.e.(i) +$ (tcm.(i).(s)) +$ l.beta.(s) +$ r.beta.(s) -$ !best
    in
    let d' = init2 states states init_d' in
    let init_e s =
        let best = ref infinity in
        for x = states - 1 downto 0 do
            store_min (d'.(x).(s)) best
        done;
        !best
    in
    let mybest =
        let abest = 
            match a.best_states with
            | [] -> 
                    (* We pick the smallest of the parent *)
                    let cost = ref infinity in
                    let res = ref [] in
                    Array.iteri (fun pos poscost ->
                        if !cost = poscost then res := pos :: !res
                        else if !cost < poscost then begin 
                            cost := poscost;
                            res := [pos];
                        end;) a.s;
                    !res
            | x -> x
        in
        let lst = ref [] in
        let best_cost = ref infinity in
        Array.iteri (fun pos cost -> 
            let mycost = List.fold_left (fun best parentpos ->
                min best (cost + tcm.(pos).(parentpos))) infinity abest
            in
            if mycost = !best_cost then 
                lst := pos :: !lst
            else if mycost < !best_cost then begin
                best_cost := mycost;
                lst := [pos]
            end else ()) n.s;
        !lst
    in
    let e = Array.init states init_e in
    { n with e = e; m = None; best_states = mybest }

let median_3 a n l r =
    let tcm = n.tcm in
    let elts = Array.init (Array.length n.elts)
        (fun i ->
             elt_median_3 tcm a.elts.(i) n.elts.(i) l.elts.(i) r.elts.(i))
    in { n with elts = elts }

let elt_exists_shared_state r a =
    let states = Array.length r.s in
    (* a state is optimal if e.(i) = 0 *)
    let rec exists n =
        if n = states
        then false
        else if r.e.(n) = 0 && a.e.(n) = 0
        then true
        else exists (n + 1)
    in exists 0

let elt_dist_2 tcm r a d =
    let states = Array.length r.s in

    (* We first check whether there are shared states between r and a, or
       between r and d.  If so, we return a delta of 0. *)
    if elt_exists_shared_state r a || elt_exists_shared_state r d
    then 0.
    else begin
        
        (* We need the array M to find the delta.  We calculate this the first
           time, then cache it.  This is safe because we create a new record each time we
           do a downpass or uppass. *)
        let m =
            match d.m with
            | Some m -> m
            | None -> begin
                  let init_d'' i =
                      let best = ref infinity in
                      for x = states - 1 downto 0 do
                          store_min ((tcm.(i).(x)) +$ d.beta.(x)) best
                      done;
                      let e = a.e.(i) in
                      fun s -> e +$ (tcm.(i).(s)) +$ d.beta.(s) -$ !best
                  in
                  let d'' = init2 states states init_d'' in

                  let init_m s =
                      let best = ref infinity in
                      for x = states - 1 downto 0 do
                          store_min d''.(x).(s) best
                      done;
                      !best
                  in
                  let m = Array.init states init_m in
                  d.m <- Some m;
                  m
              end
        in

        (* Find the best value *)
        let best = ref infinity in
        for x = states - 1 downto 0 do
            let m = m.(x) in
            let tcm = tcm.(x) in
            let e = r.e in
            for y = states - 1 downto 0 do
                store_min (m +$ (tcm.(y)) +$ e.(y)) best
            done
        done;
        float_of_cost !best
    end

let dist_2 r a d =
    let acc = ref 0. in
    for i = (Array.length r.elts) - 1 downto 0 do
        acc := !acc +. elt_dist_2 r.tcm r.elts.(i) a.elts.(i) d.elts.(i)
    done;
    !acc

let elt_to_formatter attr d tcm elt elt_parent : Xml.xml Sexpr.t =
    let module T = Xml.Characters in
    match attr with
    | [_, `String x] -> 
            let cost, lst =
                if x = Xml.Nodes.preliminary then
                    (* In this case we just pick the smallest of all mine *)
                    let lst = ref []
                    and cost = ref infinity in
                    let () =
                        Array.iteri (fun pos mycost ->
                            if mycost < !cost then begin
                                cost := mycost;
                                lst := [pos];
                            end else if mycost = !cost then lst := pos :: !lst
                            else ()) elt.s
                    in
                    !cost, !lst
                else if x = Xml.Nodes.final then 
                    let of_parent = elt_parent.best_states in
                    let (_, cost, lst) = Array.fold_left (fun ((pos, min, minlist)  as acc) x ->
                        let of_parent =
                            match of_parent with
                            | [] -> [pos]
                            | _ -> of_parent
                        in
                        if is_inf x then (pos + 1, min, minlist)
                        else
                            let dists = List.map (fun y -> 
                                tcm.(pos).(y) + x) of_parent in
                            let didit = ref false in
                            let newpos = pos + 1 in
                            List.fold_left (fun (_, min, minlist) x ->
                            if x < min then begin
                                didit := true;
                                (newpos, x, [pos])
                            end else if x = min && not !didit then begin 
                                didit := true;
                                (newpos, x, (pos :: minlist))
                            end else (newpos, min, minlist)) acc dists)  (0, max_int, []) elt.s
                    in
                    cost, lst
                else assert false
            in
            let create x = 
                (PXML -[T.value] 
                    [ `String (Data.to_human_readable d elt.ecode x) ] -- ) 
            in
            (PXML 
                -[T.sankoff] 
                    (* Attributes *)
                    ([T.name] = [`String (Data.code_character elt.ecode d)])
                    ([T.cost] = [`Int cost])
                    ([T.definite] = [`Bool (cost > 0)])
                    ([attr])

                    (* Contents *)
                    { `Set (List.map create lst) } 
                --)
    | _ -> assert false


let to_formatter attr a (parent : t option) d : Xml.xml Sexpr.t list =
    let items = Array.to_list a.elts in
    let items_parent = match parent with 
    | Some parent -> Array.to_list parent.elts 
    | None -> items   
    in 
    let tcm = a.tcm  in
    List.map2 (elt_to_formatter attr d tcm) items  items_parent

let make_onestate code tcm state =
    { empty with
        code = code;
        tcm = tcm;
        elts = [|canonize tcm { empty_elt with
                     ecode = code;
                     s = Array.init (Array.length tcm)
                         (fun x ->
                              if x = state
                              then 0
                              else infinity)}|];
    }

let make_randstate code tcm =
    make_onestate code tcm (Random.int (Array.length tcm))

let make_n_randstate scode codefn n tcm =
    let states = Array.length tcm in
    { empty with
          code = scode;
          tcm = tcm;
          elts =
            Array.init n
                (fun _ ->
                     let state = Random.int states in
                     canonize tcm
                         { empty_elt with
                               ecode = codefn ();
                               s = Array.init states
                                 (fun x ->
                                      if x = state
                                      then 0
                                      else infinity) }
                )
    }

let make_sank code tcm eltlist =
    let nstates = Array.length tcm in
    { empty with
          code = code;
          tcm = tcm;
          elts =
            let list = List.map
                (fun (ecode, states) ->
                     canonize tcm
                         { empty_elt with
                               ecode = ecode;
                               s = Array.init nstates
                                 (fun x ->
                                      if List.mem x states
                                      then 0
                                      else infinity) } )
                eltlist
            in
            Array.of_list list
    }
                                          
                               

let make_random_tcm ?(max=7) len =
    let res = init2 len len (fun _ _ -> 0) in
    for i = 0 to len - 1 do
        for j = i + 1 to len - 1 do
            let n = 1 + Random.int (max - 1) in
            res.(i).(j) <- n;
            res.(j).(i) <- n
        done
    done;
    res

let of_parser tcm (arr, taxcode) mycode =
    let nstates = Array.length tcm in
    let all_states = Array.to_list (Array.init nstates (fun x -> x)) in
    let make_elt (elt, ecode) =
        let states = 
            match elt with
            | Some (`List states) -> states
            | Some (`Bits states) -> BitSet.to_list states
            | None -> all_states
        in
        assert (List.fold_left (fun acc x -> acc && x < nstates) true states);
        canonize tcm
            { empty_elt with
                  ecode = ecode;
                  s = Array.init nstates
                    (fun i ->
                         if List.mem i states
                         then 0
                         else infinity);}
    in

    let elts = Array.map make_elt arr in
    ({ empty with
           elts = elts;
           tcm = tcm;
           code = mycode;
     }, taxcode)

let reroot_elt tcm old p q =
    (* See Goloboff 1998, p234 (p6) *)
    let states = Array.length p.s in
    let tempjk = Array.init states
        (fun j ->
             Array.init states
                 (fun k ->
                      let min = ref max_int in
                      for x = 0 to states - 1 do
                          let c = tcm.(x).(j) + tcm.(x).(k) in
                          if c < !min
                          then min := c
                      done;
                      p.e.(j) +$ q.e.(k) -$ (!min))) in
    let alpha5 = Array.init states
        (fun i ->
             Array.init states
                 (fun j ->
                      Array.init states
                          (fun k ->
                               (tcm.(i).(j) + tcm.(i).(k))
                               +$ tempjk.(j).(k)))) in
    let newe = Array.init states
        (fun i ->
             let best = ref infinity in
             for x = 0 to states - 1 do
                 for y = 0 to states - 1 do
                     store_min alpha5.(i).(x).(y) best
                 done
             done;
             !best) in
    { old with e = newe }

let reroot old p q =
    let tcm = old.tcm in
    let newelts = Array.init (Array.length old.elts)
        (fun i ->
             reroot_elt tcm old.elts.(i) p.elts.(i) q.elts.(i)) in
    { old with elts = newelts }

let get_elt_array {elts=e} = e

let elt_filter f elt = 
    Array_ops.filter f elt

let filter f t =
    { t with elts = elt_filter f t.elts }

let f_codes t codes = 
    let check x = All_sets.Integers.mem x.ecode codes in
    filter check t

let f_codes_comp t codes = 
    let check x = not (All_sets.Integers.mem x.ecode codes) in
    filter check t

let cardinal t = Array.length t.elts


(* We will program branch and bound only for sankoff characters *)
let large_number = max_int / 8 

let make_leaf states x = 
    let arr = Array.make states large_number in
    List.iter (fun x -> arr.(x) <- 0) x;
    Parser.Tree.Leaf arr

let rec aux_bb bound lst mtx tree = 
    let states = Array.length mtx in
    let calc_node left right = 
        let get_arr left = 
            match left with
            | Parser.Tree.Leaf x
            | Parser.Tree.Node (_, x) -> x 
        in
        let left = get_arr left 
        and right = get_arr right in
        Array.init states (fun pos ->
            let res = ref large_number in
            for i = 0 to states - 1 do
                for j = 0 to states - 1 do
                    let cost =  left.(i) + right.(j) + mtx.(i).(pos) + mtx.(pos).(j) in
                    if cost < !res then res := cost
                done;
            done;
            !res) 
    in
    let join a b = Parser.Tree.Node ([a; b], calc_node a b) in
    let get_cost tree =
        let arr = 
            match tree with
            | Parser.Tree.Node (_, arr) | Parser.Tree.Leaf arr -> arr 
        in
        Array.fold_left min large_number arr 
    in
    let rec cleanup_list (best_cost, acc) lst =
        match lst with
        | h :: t -> 
                let cost = get_cost h in
                if best_cost > cost then 
                    cleanup_list (best_cost, (h :: acc)) t
                else cleanup_list (best_cost, acc) t
        | [] -> best_cost, acc
    in
    let rec append_everywhere terminal tree =
        match tree with
        | Parser.Tree.Leaf x -> [join terminal tree]
        | Parser.Tree.Node ([a; b], _) ->
                (join b (join a terminal)) :: (join a (join b terminal)) ::
                    (List.map (join a) (append_everywhere terminal b)) @
                    (List.map (join b) (append_everywhere terminal a))
        | Parser.Tree.Node _ -> assert false
    in
    match lst with
    | [] -> 
            get_cost tree
    | h :: t ->
            let res = append_everywhere (make_leaf states h) tree in
            let _, res = cleanup_list (bound, []) res in
            match res with
            | [] -> bound
            | trees -> List.fold_left min bound (List.map (aux_bb bound t mtx) trees)
                 
let maximum_states_for_warning = 8

let bb mtx characters =
    let states = Array.length mtx in
    if states > 8 then
        Status.user_message Status.Error 
        "Computing CI and RI for Sankoff characters takes exponential time! you
        have more than 8 states for at least one character, so you should expect
        a long execution time for this calculation ... good luck!";
    match characters with
    | h :: t -> 
            aux_bb large_number t mtx (make_leaf states h)
    | [] -> 0


(* Sankoff characters can't have an exact minimum possible cost (the problem is
* NP-Hard by itself), therefore, we will simply define a function that returns
* all possible sets of assignments, and that way we will be able to search for
* the best cost _found_ during a search. The function assumes that the cost
* matrix is metric, otherwise the method that this function is used for will be
* ... ehm ... broken. *)
let get_all_possible_assignments (elts : int list option list) = 
    let rec aux_all_possiblities left sets = 
        match left with
        | [] -> sets
        | (Some []) :: t 
        | None :: t -> aux_all_possiblities t sets
        | (Some lst) :: t ->
                let sets = 
                    List.flatten 
                    (List.map (fun x -> List.map (fun y -> 
                        if All_sets.Integers.mem x y then y 
                        else All_sets.Integers.add x y) sets) lst)
                in
                aux_all_possiblities t sets
    in
    let rec remove_supersets lst acc =
        match lst with
        | h :: t -> 
                let res = List.filter (fun x ->
                    not (All_sets.Integers.subset h x)) t in
                remove_supersets res (h :: acc) 
        | [] -> acc
    in
    let x = 
        let x = aux_all_possiblities elts [All_sets.Integers.empty] in
        let x = List.sort (fun x y -> 
            (All_sets.Integers.cardinal x) - (All_sets.Integers.cardinal y)) 
            x 
        in
        remove_supersets x []
    in
    List.map All_sets.Integers.elements x

let min_possible_cost mtx (elts : Parser.SC.static_state list) = 
    let all_possible = 
        let rec filter_none acc lst =
            match lst with
            | None :: t -> filter_none acc t
            | (Some x) :: t -> 
                    filter_none ((Parser.SC.static_state_to_list x) :: acc) t
            | [] -> acc
        in
        filter_none [] elts
    in
    let res = bb mtx all_possible in
    float_of_int res

let max_possible_cost mtx elts = 
    let mtx = Array.map (Array.map (fun x -> x * (-1))) mtx in
    let res = min_possible_cost mtx elts in
    ((-.1.) *. res)
