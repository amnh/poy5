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
let () = SadmanOutput.register "SankCS" "$Revision: 2871 $"

let debug = false

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

external get_e_array : elt -> 
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


(*let elt_to_full_string a =
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
*)  

let assert_ninf a x y=
    assert (
        if is_infinity a then begin
         (*print_string (elt_to_full_string x);
            print_string (elt_to_full_string y);*)
            false
        end else true);
    a



(* An empty character. This is used only for the parsing process. 
let empty =
    { elts = [||];
      code = (-1);
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
*)

(* nobody calls this
let codes {elts=elts} =
    let len = Array.length elts in
    let acc = ref [] in
    for i = len - 1 downto 0 do
        acc := (elts.(i).ecode) :: !acc
    done;
    !acc
*)

(* The cost of a transformation occurring in the following order i -> j + i -> k
 * as defined in the transformation cost matrix par. *)
let median_cost par i j k = 
    (par.(i).(j) + par.(i).(k))

let get_min a = Array.fold_left cost_min infinity a


external get_code : t -> int = "sankoff_CAML_get_code"
(* The code of the character a 
let code a = a.code *)

(*nobody calls this
let mem cs t = match cs with
    | None    -> true
    | Some [] -> false
    | Some xs -> List.mem (code t) xs
*)

(* Code of the element 
let ecode e = e.ecode*)

(* The transformation cost matrix associated with character a --nobody calls
* this 
let tcm a = a.tcm*)

(*let set_tcm cm t  = {t with tcm = cm;} nobody calls this*)

let print_tcm tcm =
    Array.iter (fun row ->
                    Array.iter (fun v ->
                                    print_string (string_of_int v ^ " ")) row;
                    print_newline ())
        tcm

(* How many states there are-- nobody calls this
let nstates a = Array.length a.tcm*)

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

(* The random state generator, using the random character generator created by
 * rand_gen () *)
(* let make_rand (num, code, states, tcm) = *)
(*     let states = Array.init states (fun _ -> Cost (1 + (Random.int states))) in *)
(*     canonize *)
(*         { empty with s = states; code = code; tcm = tcm } *)

(*let get_minstates {e=a} =
    let (_, list) = 
        Array.fold_right (fun v (num, lst) ->
                              if v = 0
                              then (num - 1, num :: lst)
                              else (num - 1, lst)) a (Array.length a - 1, [])
    in list
no one calls this *)



external median_cside : int -> t -> t -> t = "sankoff_CAML_median"

external get_sumcost : t -> int = "sankoff_CAML_get_sumcost"


let median median_node_code a b =
    let debug = false in
    if debug then begin
        let tca = get_taxon_code a in
        let tcb = get_taxon_code b in
        Printf.printf "SankCS.median,median_node_code=%d,taxon code a/b = %d,%d\n%!" median_node_code
        tca tcb;
        (*let alst = to_list a in
        let blst = to_list b in
        List.iter (fun (ec, states_arr) -> 
            Printf.printf "nodeA=%d, ecode=%d, " tca ec; 
            Utl.printIntArr states_arr;
        ) alst;
        List.iter (fun (ec, states_arr) -> 
            Printf.printf "nodeB=%d, ecode=%d, " tcb ec; 
            Utl.printIntArr states_arr;
        ) blst;*)
    end;
    let med = median_cside median_node_code a b in
    let sumcost = get_sumcost med in
    if (sumcost<0) then 
    error_user_message "subtree cost exceed max_int of current system.";
    let cost = float_of_int sumcost in
    if debug then Printf.printf "return median with cost = %f\n%!" cost;
    med,cost
    (*
    assert (a.code = b.code);
    assert (a.tcm = b.tcm);
    if debug then Printf.printf "median\n a = %s\n b = %s\n%!" (to_string a) (to_string b);
    let res = { a with elts = Array.mapi (fun i -> elt_median a.tcm a.elts.(i)) b.elts }
    in
    if debug then Printf.printf "median <- %s\n%!" (to_string res);
    res*)



(* move to cside 
* Calculates the distance between two characters a and b *)
(* Note that we only calculate the _added_ distance 
let elt_distance tcm a b =
    let median = elt_median tcm a b in
    let get_cost m =
        float_of_cost (assert_ninf (Array.fold_left cost_min infinity m.s) a b) in
    let med_cost = get_cost median in
    let a_cost = get_cost a in
    let b_cost = get_cost b in
    (*this is wrong*)
    med_cost -. a_cost -. b_cost
*)
(* to do
let elt_distance_and_median tcm a b = 
    let median = elt_median tcm a b in
    let get_cost m =
        float_of_cost (assert_ninf (Array.fold_left cost_min infinity m.s) a b) in
    let med_cost = get_cost median in
    let a_cost = get_cost a in
    let b_cost = get_cost b in
    med_cost -. a_cost -. b_cost, 
    median
*)

external distance_cside : t -> t -> int = "sankoff_CAML_distance"
(** [distance a b] return the sankoff distance. Note that it calls
* [elt_distance], which will call [elt_median]. if you need median and distance,
* don't call two functions seperately, use [distance_and_median] instead*)
let distance a b =
    let debug = false in
    if debug then Printf.printf "SankCS.distance\n%!";
    let dis = distance_cside a b in
    float_of_cost dis
    (*let tcm = a.tcm in
    let acc = ref 0. in
    for i = Array.length a.elts - 1 downto 0 do
        acc := !acc +. elt_distance tcm a.elts.(i) b.elts.(i)
    done;
    if debug then Printf.printf "distance = %f\n%!" !acc;
    !acc
*)

(* to do : get the distance function right, then back to this
external distance_and_median_cside : t -> t -> int * t = ""   

let distance_and_median a b =
    let dis,med = distance_and_median_cside a b in
    float_of_cost dis, med
*)
(*    if debug then Printf.printf "SankCS.distance_and_median\n%!";
    assert (a.code = b.code);
    assert (a.tcm = b.tcm);
    let tcm = a.tcm in
    let acc = ref 0. in
    let med = 
        { a with elts = 
          Array_ops.map_2 (fun aelt belt -> 
            let disi,medi = elt_distance_and_median tcm aelt belt in
            acc := !acc +. disi;
            medi 
            ) a.elts b.elts
        }
    in
    if debug then Printf.printf "distance = %f,median = %s\n%!" 
    !acc (to_string med);
    !acc, med*)


(* Compares two characters a and b. Note that this comparison is used basically
 * in the sets, and therefore, comparing the codes is enough (all the characters in
 * a set have different codes). 
let compare_codes a b = 
    compare a.code b.code
no one calls this*)

(* Compares the data in the two characters a and b.  This is used to find cases
 * where a graph update does not change values far away in the graph, so that
 * propagation can be stopped. *)
external compare_cside : t -> t -> int = "sankoff_CAML_compare_eltarr"

let compare_data a b =
    compare_cside a b
    (*compare a b*)

(* A dummy function to be improved later using the official parser 
let parse a = [("", empty) ]
*)

let init2 len1 len2 fn =
    Array.init len1
        (fun i ->
             let fn = fn i in
             Array.init len2 (fun j -> fn j))


(* move to cside
* This algorithm is taken from Goloboff 1998.  Calculate D', then the final E
 * value from that.
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
move this to c side*)


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
        (*let alst = to_list a in
        let blst = to_list n in
        List.iter (fun (ec, states_arr) -> 
            Printf.printf "nodeA=%d, ecode=%d, " tca ec; 
            Utl.printIntArr states_arr;
        ) alst;
        List.iter (fun (ec, states_arr) -> 
            Printf.printf "nodeN=%d, ecode=%d, " tcn ec; 
            Utl.printIntArr states_arr;
        ) blst; *)
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

(*
let elt_return_shared_states r a =
    let states = Array.length r.s in
    let acc = ref [] in
    (* a state is optimal if e.(i) = 0 *)
    let rec exists n =
        if n = states then !acc
        else if r.e.(n) = 0 && a.e.(n) = 0
        then begin
            acc := n :: (!acc);
            exists (n + 1)
        end
        else exists (n + 1)
    in exists 0
*)

(* not in use any more
let elt_exists_shared_state r a =
    let states = Array.length r.s in
    (* a state is optimal if e.(i) = 0 *)
    let rec exists n =
        if n = states
        then false
        else if r.e.(n) = 0 && a.e.(n) = 0
        then true
        else exists (n + 1)
    in exists 0 *)


(* move to cside
let elt_dist_2 tcm r a d =
    let debug = false in
    if debug then Printf.printf "elt_dist_2, r.elts=%s\n a.elts=%s \n d.elts=%s \n%!"
    (elt_to_string r) (elt_to_string a) (elt_to_string d);
    let states = Array.length r.s in
    (* We first check whether there are shared states between r and a, or
       between r and d.  If so, we return a delta of 0. *)
    let shared_states_ra = elt_return_shared_states r a 
    and shared_states_rd = elt_return_shared_states r d in
    (*if elt_exists_shared_state r a || elt_exists_shared_state r d
    then 0. this won't be true if cost between same state is non-zero*)
    if (List.length shared_states_rd)>0 || (List.length shared_states_ra)>0 then
        float_of_cost (2*(min (get_min_cost_between_same_states shared_states_ra tcm)
        (get_min_cost_between_same_states shared_states_rd tcm)))
    else begin
        (* We need the array M to find the delta.  We calculate this the first
           time, then cache it.  This is safe because we create a new record each time we do a downpass or uppass. *)
        let m =
            match d.m with
            | Some m -> m
            | None -> begin
                  let init_d'' i =
                      let best = ref infinity in
                      let bests = ref (-1) in
                      for x = states - 1 downto 0 do
                          let oldbest = !best in
                          store_min ((tcm.(i).(x)) +$ d.beta.(x)) best;
                          if !best<oldbest then
                              bests := x;
                      done;
if debug then 
Printf.printf "when nodeA take statei.%d, min(tcm.i.x + nodeD.beta.x) = %d with \
nodeM take states=%d\n%!" i !best !bests;
                      let e = a.e.(i) in
                      fun s -> 
if debug then
Printf.printf "when nodeA take statei.%d,nodeM take states.%d,\
E(i,A(d))=%d + tcm.i.s=%d + nodeD.beta.s=%d - best(%d)\n%!" 
      i s e tcm.(i).(s) d.beta.(s) !best;
                          e +$ (tcm.(i).(s)) +$ d.beta.(s) -$ !best
                  in
                  let d'' = init2 states states init_d'' in
                  let init_m s =
                      let best = ref infinity in
                      let besti = ref [] in
                      for x = states - 1 downto 0 do
                          let oldbest = !best in
                          store_min d''.(x).(s) best;
                          if !best<=oldbest then
                              besti := x::(!besti);
                      done;
                      (*when the best assignment to nodeM(s) is included in
                  * nodeA, an extra cost between same states must be added *)
                      if (List.mem s !besti) then
                          best := !best + tcm.(s).(s);
                      if debug then Printf.printf 
"when nodeM take state %d,M <- %d, best statei for nodeA is %!" s !best;
if debug then Utl.printIntList (!besti);
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
if debug then Printf.printf "when nodeM take state.%d, nodeR take state.y=%d,\
M(%d) + tcm(%d) + E(%d)\n%!" x y m tcm.(y) e.(y);
                store_min (m +$ (tcm.(y)) +$ e.(y)) best
            done
        done;
        float_of_cost !best
    end
*)

external dist_2_cside : t -> t -> t -> int = "sankoff_CAML_dist_2"

let dist_2 r a d =
    let dis = dist_2_cside d a r in
    float_of_cost dis
    (*let debug = false in
    if debug then Printf.printf "dist_2 on r = %s\n a = %s\n d = %s\n%!"
    (to_string r) (to_string a) (to_string d);
    let acc = ref 0. in
    for i = (Array.length r.elts) - 1 downto 0 do
        let disti = elt_dist_2 r.tcm r.elts.(i) a.elts.(i) d.elts.(i) in
        if debug then Printf.printf "i=%d, acc += %f\n%!" i disti;
        acc := !acc +. disti;
    done;
    !acc
*)


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



let get_elts_from_t thist =
    let num_elts = get_num_elts thist in
    Array.init num_elts (fun i -> get_elt thist i)
    

let to_formatter attr a (parent : t option) d : Xml.xml Sexpr.t list =
    let items = Array.to_list (get_elts_from_t a) in (*Array.to_list a.elts in*)
    let items_parent = match parent with 
    | Some parent -> Array.to_list (get_elts_from_t parent) (*parent.elts*) 
    | None -> items   
    in 
    let tcm = get_tcm a in (*a.tcm  in*)
    let idx = ref 0 in
    List.map2 (fun this_elt parent_elt -> 
        let res = elt_to_formatter attr d tcm !idx this_elt parent_elt in
        idx := !idx + 1;
        res) items  items_parent

(* no one calls these 
let make_onestate code tcm state =
    { code = code;
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
    { code = scode;
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
    { code = code;
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
*)

(* move to cside
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
*)


external create_eltarr_cside : int -> int -> int ->
    (int32, Bigarray.int32_elt, Bigarray.c_layout) Bigarray.Array1.t -> 
    (int32, Bigarray.int32_elt, Bigarray.c_layout) Bigarray.Array2.t -> 
    (int32, Bigarray.int32_elt, Bigarray.c_layout) Bigarray.Array2.t -> 
        t = "sankoff_CAML_create_eltarr_bytecode" "sankoff_CAML_create_eltarr"

let create_eltarr taxcode mycode nstates ecode_arr state_arrarr tcm =
    let tcm_int32 = Utl.int_to_int32_mat tcm in
    let tcm_bigarr 
    = Bigarray.Array2.of_array Bigarray.int32 Bigarray.c_layout tcm_int32 in
    let state_bigarr = Bigarray.Array2.of_array Bigarray.int32 Bigarray.c_layout
    state_arrarr in
    let ecode_bigarr = Bigarray.Array1.of_array Bigarray.int32 Bigarray.c_layout
    ecode_arr in
    create_eltarr_cside taxcode mycode nstates ecode_bigarr state_bigarr tcm_bigarr

(*create sankoff chr from input file*)
let of_parser tcm (arr, taxcode) mycode =
    let debug = false in
    if debug then Printf.printf "SankCS.of_parser,taxcode=%d,mycode=%d\n%!" taxcode mycode;
    (*let tcm_int32 = Utl.int_to_int32_mat tcm in
    let tcm_int32 = Utl.int_to_int32_mat tcm in
    let tcm_bigarr 
    = Bigarray.Array2.of_array Bigarray.int32 Bigarray.c_layout tcm_int32 in*)
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
    create_eltarr taxcode mycode nstates ecode_arr state_arrarr tcm,
    (*
    let states_bigarr = Bigarray.Array2.of_array Bigarray.int32 Bigarray.c_layout
    state_arrarr in
    let ecode_bigarr = Bigarray.Array1.of_array Bigarray.int32 Bigarray.c_layout
    ecode_arr in
    create_eltarr taxcode mycode nstates ecode_bigarr states_bigarr
    tcm_bigarr,*)
    taxcode
    

(* no one calls this
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
no one calls this *)

(*
let get_elt_array {elts=e} = e*)

(*  no longer in use
let elt_filter f elt = 
    Array_ops.filter f elt

let filter f t =
    { t with elts = elt_filter f t.elts }
no longer in use*)

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
    (*
    let check x = All_sets.Integers.mem x.ecode codes in
    filter check t*)

let f_codes_comp t codes = 
    let ecodelst = All_sets.Integers.elements codes in
    let ecodearr = Array.of_list ecodelst in
    let ecodearr = Array.map (fun x -> Int32.of_int x ) ecodearr in
    let ecode_bigarr = Bigarray.Array1.of_array Bigarray.int32 Bigarray.c_layout
    ecodearr
    in
    filter_character t ecode_bigarr 1
    (*
    let check x = not (All_sets.Integers.mem x.ecode codes) in
    filter check t
    *)

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


(* Sankoff characters can't have an exact minimum possible cost (the problem is
* NP-Hard by itself), therefore, we will simply define a function that returns
* all possible sets of assignments, and that way we will be able to search for
* the best cost _found_ during a search. The function assumes that the cost
* matrix is metric, otherwise the method that this function is used for will be
* ... ehm ... broken. 
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
no one calls this *)

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
