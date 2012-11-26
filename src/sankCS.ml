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
let () = SadmanOutput.register "SankCS" "$Revision: 2831 $"

let debug = false

external register_eltarr : unit -> unit = "sankoff_CAML_register_eltarr"
external register_elt : unit -> unit = "sankoff_CAML_register_elt"
let () = register_eltarr ()
let () = register_elt ()

let error_user_message format = Printf.ksprintf (Status.user_message Status.Error) format
    
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
(*
type elt in sankoff.c*)
type elt
(*
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
*)
type t
(*
(* The Sankoff character type *)
type t = { 
    (* The code the character belongs to. Every pair of characters with the same
     * code correspond to homologous units *)
    code : int;
    (* The transformation cost matrix to be used for this character (and its
     * homologous characters) *)
    tcm : cm;

    (* A bunch of individual characters *)
    elts : elt array;
}
*)

(* A default string representation as a list with the minimal states (those with
* minimal cost in the state array). *)

(*function for fixed states*)
external get_best_child_state : t -> int -> int = "sankoff_CAML_get_best_child_state"

external set_gc_alloc_max : int -> unit = "sankoff_GC_custom_max"

external get_ecode : elt -> int = "sankoff_CAML_get_ecode"

external get_taxon_code : t -> int = "sankoff_CAML_get_taxon_code"

external get_num_elts : t -> int = "sankoff_CAML_get_num_elts" 

external get_elt : t -> int -> elt = "sankoff_CAML_get_elt"

external get_tcm : t ->
    (int32,Bigarray.int32_elt,Bigarray.c_layout) Bigarray.Array2.t =
        "sankoff_CAML_get_tcm"

external get_earray_cside : elt -> 
    (int32,Bigarray.int32_elt,Bigarray.c_layout) Bigarray.Array1.t = "sankoff_CAML_get_e_array"

external get_extra_cost_for_root : t -> int = "sankoff_CAML_get_extra_cost_for_root"


let to_string s = " to do "

external get_states_cside : elt -> int ->
    (int32,Bigarray.int32_elt,Bigarray.c_layout) Bigarray.Array1.t
    = "sankoff_CAML_get_states"
    
let to_list s =
    let num_elts = get_num_elts s in
    let res = Array.init num_elts (fun eltNO ->
        let thiselt = get_elt s eltNO in
        let states_bigarr = get_states_cside thiselt 1 in
        let states = Array.init (Bigarray.Array1.dim states_bigarr) 
        (fun x -> 
            Int32.to_int (Bigarray.Array1.get states_bigarr x)
        ) in
        get_ecode thiselt,
        states
    ) in
    Array.to_list res

(*function for fixed states, where we have only one elt for each t*)
(*1=states,2=leftstates,3=rightstates*)
let get_states s this_or_left_or_right =
    let num_elts = get_num_elts s in
    assert(num_elts==1);
    let thiselt = get_elt s 0 in
    let states_bigarr = get_states_cside thiselt this_or_left_or_right in
    let states = Array.init (Bigarray.Array1.dim states_bigarr) 
    (fun x -> 
        Int32.to_int (Bigarray.Array1.get states_bigarr x)
    )  in
    states

(*function for fixed states, where we have only one elt for each t*)
let get_earray s =
    let num_elts = get_num_elts s in
    assert(num_elts==1);
    let thiselt = get_elt s 0 in
    let e_bigarr = get_earray_cside thiselt in
    Array.init
        (Bigarray.Array1.dim e_bigarr)
        (fun x -> Int32.to_int (Bigarray.Array1.get e_bigarr x))


(*function for fixed states, where we have only one elt for each t*)
let get_earray_camlside thiselt = 
    let e_bigarr = get_earray_cside thiselt in
    let earr = Array.init (Bigarray.Array1.dim e_bigarr) 
    (fun x -> 
        Int32.to_int (Bigarray.Array1.get e_bigarr x)
    )  in
    earr



let assert_ninf a x y=
    assert (
        if is_infinity a then begin
         (*print_string (elt_to_full_string x);
            print_string (elt_to_full_string y);*)
            false
        end else true);
    a


(* The cost of a transformation occurring in the following order i -> j + i -> k
 * as defined in the transformation cost matrix par. *)
let median_cost par i j k = 
    (par.(i).(j) + par.(i).(k))

let get_min a = Array.fold_left cost_min infinity a


external get_code : t -> int = "sankoff_CAML_get_code"


let print_tcm tcm =
    Array.iter (fun row ->
                    Array.iter (fun v ->
                                    print_string (string_of_int v ^ " ")) row;
                    print_newline ())
        tcm



(* We'll be using this a lot ...
 * (make it monomorphic; more efficient?) *)
let store_min (a : cost) b = if a <$ !b then b := a

let is_inf x = is_infinity x




(* Even OTUs must have [beta] values, so we may as well store [e], also. *)

(* A random character generator specification *)
let rand_gen () = 
    let states = 4 in
    (* Random.int only accepts arguments up to 2^30 ... *)
    let code = Random.int 1073741823
    and tcm = Array.make_matrix 4 4 ( 1) in
    let num = Random.int 15 in 
    num, code, states, tcm




external median_cside : int -> t -> t -> t = "sankoff_CAML_median"

external get_sumcost : t -> int = "sankoff_CAML_get_sumcost"


let median median_node_code a b =
    let debug = false in
    if debug then begin
        let tca = get_taxon_code a in
        let tcb = get_taxon_code b in
        Printf.printf "SankCS.median,median_node_code=%d,taxon code a/b = %d,%d\n%!" median_node_code
        tca tcb;
    end;
    let med = median_cside median_node_code a b in
    let sumcost = get_sumcost med in
    if (sumcost<0) then 
    error_user_message "subtree cost exceed max_int of current system.";
    let cost = float_of_int sumcost in
    if debug then Printf.printf "return median with cost = %f\n%!" cost;
    med,cost
    




external distance_cside : t -> t -> int = "sankoff_CAML_distance"
(** [distance a b] return the sankoff distance. Note that it calls
* [elt_distance], which will call [elt_median]. if you need median and distance,
* don't call two functions seperately, use [distance_and_median] instead*)
let distance a b =
    let debug = false in
    if debug then Printf.printf "SankCS.distance\n%!";
    let dis = distance_cside a b in
    float_of_cost dis
    



(* Compares the data in the two characters a and b.  This is used to find cases
 * where a graph update does not change values far away in the graph, so that
 * propagation can be stopped. *)
external compare_cside : t -> t -> int = "sankoff_CAML_compare_eltarr"

let compare_data a b =
    compare_cside a b
    (*compare a b*)



let init2 len1 len2 fn =
    Array.init len1
        (fun i ->
             let fn = fn i in
             Array.init len2 (fun j -> fn j))



external median_3_cside : t -> t -> t -> t -> t = "sankoff_CAML_median_3"

let median_3 a n l r =
    let debug = false in
    let tcn = get_taxon_code n in
    let tcr = get_taxon_code r in
    if debug then begin
        let tcl = get_taxon_code l in
        let tca = get_taxon_code a in
        Printf.printf "median 3, taxon code, nodeN=%d =?= nodeR=%d, nodeL=%d, \
        nodeA=%d\n%!" tcn tcr tcl tca;
    end;
    if tcn=tcr then(*leafnode*)
        n
    else
        median_3_cside a n l r


let get_min_cost_between_same_states same_states cm = 
    let best = ref infinity in
    List.iter (fun state ->
        store_min cm.(state).(state) best;
    ) same_states;
    if debug then begin
       Printf.printf "get_min_cost_between_same_states:%!";
       Utl.printIntList same_states; 
       Printf.printf "min_cost = %d\n%!" !best;
    end;
    !best



external dist_2_cside : t -> t -> t -> int = "sankoff_CAML_dist_2"

let dist_2 r a d =
    let dis = dist_2_cside d a r in
    float_of_cost dis
   

let elt_to_formatter attr d tcm idx elt elt_parent : Xml.xml Sexpr.t =
    let module T = Xml.Characters in
    match attr with
    | [_, `String x] -> 
            let cost, lst =
                let states_bigarr = get_states_cside elt 1 in
                let states = Array.init (Bigarray.Array1.dim states_bigarr) 
                (fun x -> 
                Int32.to_int (Bigarray.Array1.get states_bigarr x)
                ) in
                let idx = ref (-1) in
                let best_states_idxlst,bestcost = Array.fold_left (fun (acc,best) x ->
                    idx := !idx + 1;
                    if best>x then ([!idx],x)
                    else if best=x then ((!idx)::acc, best)
                    else (acc,best)
                ) ([],states.(0))  states in
                bestcost, best_states_idxlst
            in
            let ecode = get_ecode elt in
            let create x = 
                (PXML -[T.value] 
                    [ `String (Data.to_human_readable d ecode x) ] -- ) 
            in
            (PXML 
                -[T.sankoff] 
                    (* Attributes *)
                    ([T.name] = [`String (Data.code_character ecode d)])
                    ([T.cost] = [`Int cost])
                    ([T.definite] = [`Bool (cost > 0)])
                    ([attr])
                    (* Contents *)
                    { `Set (List.map create lst) } 
                --)
    | _ -> assert false

let elt_to_formatter_with_seq print_seq seq_array seq_alph attr d elt elt_parent : Xml.xml Sexpr.t =
    let module T = Xml.Characters in
    match attr with
    | [_, `String x] -> 
            let cost, lst, seqlst =
                let earray = get_earray_camlside elt in
                let idx = ref (-1) in
                let best_states_idxlst,best_seq_lst,bestcost = 
                    Array.fold_left (fun (idxacc,seqacc,best) x ->
                    idx := !idx + 1;
                    if best>x then ([!idx],[seq_array.(!idx)],x)
                    else if best=x then ((!idx)::idxacc,(seq_array.(!idx))::seqacc, best)
                    else (idxacc,seqacc,best)
                ) ([],[],earray.(0)) earray in
                bestcost, best_states_idxlst, best_seq_lst
            in
            let ecode = get_ecode elt in
            let seq_string_lst = 
                List.map (fun x -> 
                    let x = Sequence.del_first_char x in
                    Sequence.to_formater x seq_alph) seqlst
            in
            let create x y = 
                if print_seq then
                (PXML -[T.value] 
                    [ `String ("state:"^string_of_int(x)^",seq="^y) ] -- ) 
                else
                    (PXML -[T.value] 
                    [ `String ("state:"^string_of_int(x)) ] -- )
            in
            (PXML 
                -[T.sankoff] 
                    (* Attributes *)
                    ([T.name] = [`String (Data.code_character ecode d)])
                    ([T.cost] = [`Int cost])
                    ([T.definite] = [`Bool (cost > 0)])
                    ([attr])
                    (* Contents *)
                    { `Set (List.map2 create lst seq_string_lst) } 
                --)
    | _ -> assert false

let get_elts_from_t thist =
    let num_elts = get_num_elts thist in
    Array.init num_elts (fun i -> get_elt thist i)
    

let to_formatter attr a (parent : t option) d : Xml.xml Sexpr.t list =
    let items = Array.to_list (get_elts_from_t a) in 
    let items_parent = match parent with 
    | Some parent -> Array.to_list (get_elts_from_t parent)  
    | None -> items   
    in 
    let tcm = get_tcm a in 
    let idx = ref 0 in
    List.map2 (fun this_elt parent_elt -> 
        let res = elt_to_formatter attr d tcm !idx this_elt parent_elt in
        idx := !idx + 1;
        res) items  items_parent

(*to_formatter function for fixed states, pass the original sequence and
* alphabet in case we want to print them out*)
let to_formatter_with_seq print_seq seq_arr seq_alph attr a (parent : t option) d : Xml.xml Sexpr.t list =
    let items = Array.to_list (get_elts_from_t a) in 
    assert((List.length items)=1);(*we have only one elt per t for fixed_states*)
    let items_parent = match parent with 
    | Some parent -> Array.to_list (get_elts_from_t parent)  
    | None -> items   
    in 
    let idx = ref 0 in
    List.map2 (fun this_elt parent_elt -> 
        let res = elt_to_formatter_with_seq print_seq seq_arr seq_alph attr d this_elt parent_elt in
        idx := !idx + 1;
        res) items  items_parent

let is_identity tcm =
    let sizex = Array.length tcm in
    let res = ref true in
    for i = 0 to sizex-1 do
        if tcm.(i).(i)<>0 then res:=false;
    done;
    !res

external create_eltarr_cside : int -> int -> int -> int ->
    (int32, Bigarray.int32_elt, Bigarray.c_layout) Bigarray.Array1.t -> 
    (int32, Bigarray.int32_elt, Bigarray.c_layout) Bigarray.Array2.t -> 
    (int32, Bigarray.int32_elt, Bigarray.c_layout) Bigarray.Array2.t -> 
        t = "sankoff_CAML_create_eltarr_bytecode" "sankoff_CAML_create_eltarr"


let create_eltarr taxcode mycode nstates ecode_arr state_arrarr tcm isidentity =
    let tcm_int32 = Utl.int_to_int32_mat tcm in
    let tcm_bigarr 
    = Bigarray.Array2.of_array Bigarray.int32 Bigarray.c_layout tcm_int32 in
    let state_bigarr = Bigarray.Array2.of_array Bigarray.int32 Bigarray.c_layout
    state_arrarr in
    let ecode_bigarr = Bigarray.Array1.of_array Bigarray.int32 Bigarray.c_layout
    ecode_arr in
    let isidentity = if isidentity then 1 else 0 in (*we pass int instead of bool*)
    create_eltarr_cside isidentity taxcode mycode nstates ecode_bigarr state_bigarr tcm_bigarr

(*create sankoff chr from input file*)
let of_parser tcm (arr, taxcode) mycode =
    let debug = false in
    let iside = is_identity tcm in
    if debug then Printf.printf "SankCS.of_parser,taxcode=%d,mycode=%d,is identity=%b\n%!" taxcode mycode iside;
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
        (*infinity here is not infinity on the c side, we pass (-1) instead*)
        let state_arr =  Array.init nstates (fun i -> 
            if List.mem i states then Int32.of_int 0 
            else Int32.of_int (-1) (*infinity*)
        ) in
        Int32.of_int ecode, state_arr
    in
    let elts = Array.map make_elt arr in
    let eltlst = Array.to_list elts in 
    let ecode_lst,state_arrlst = List.split eltlst in
    let ecode_arr,state_arrarr = Array.of_list ecode_lst, Array.of_list
    state_arrlst in
    create_eltarr taxcode mycode nstates ecode_arr state_arrarr tcm iside,
    taxcode
    

external filter_character : t -> 
    (int32,Bigarray.int32_elt,Bigarray.c_layout) Bigarray.Array1.t -> int -> t = 
        "sankoff_CAML_filter_character"

let f_codes t codes =
    let ecodelst = All_sets.Integers.elements codes in
    let ecodearr = Array.of_list ecodelst in
    let ecodearr = Array.map (fun x -> Int32.of_int x ) ecodearr in
    let ecode_bigarr = Bigarray.Array1.of_array Bigarray.int32 Bigarray.c_layout
    ecodearr
    in
    filter_character t ecode_bigarr 0
    

let f_codes_comp t codes = 
    let ecodelst = All_sets.Integers.elements codes in
    let ecodearr = Array.of_list ecodelst in
    let ecodearr = Array.map (fun x -> Int32.of_int x ) ecodearr in
    let ecode_bigarr = Bigarray.Array1.of_array Bigarray.int32 Bigarray.c_layout
    ecodearr
    in
    filter_character t ecode_bigarr 1
    

let cardinal t = get_num_elts t (*Array.length t.elts*)


(* We will program branch and bound only for sankoff characters *)
let large_number = max_int / 8 

let make_leaf states x = 
    let arr = Array.make states large_number in
    List.iter (fun x -> arr.(x) <- 0) x;
    Tree.Parse.Leafp arr

let rec aux_bb bound lst mtx tree = 
    let states = Array.length mtx in
    let calc_node left right = 
        let get_arr left = 
            match left with
            | Tree.Parse.Leafp x
            | Tree.Parse.Nodep (_, x) -> x 
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
    let join a b = Tree.Parse.Nodep ([a; b], calc_node a b) in
    let get_cost tree =
        let arr = 
            match tree with
            | Tree.Parse.Nodep (_, arr) | Tree.Parse.Leafp arr -> arr 
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
        | Tree.Parse.Leafp x -> [join terminal tree]
        | Tree.Parse.Nodep ([a; b], _) ->
                (join b (join a terminal)) :: (join a (join b terminal)) ::
                    (List.map (join a) (append_everywhere terminal b)) @
                    (List.map (join b) (append_everywhere terminal a))
        | Tree.Parse.Nodep _ -> assert false
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



let min_possible_cost mtx (elts : Nexus.File.static_state list) = 
    let all_possible = 
        let rec filter_none acc lst =
            match lst with
            | None :: t -> filter_none acc t
            | (Some x) :: t -> 
                    filter_none ((Nexus.File.static_state_to_list x) :: acc) t
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
