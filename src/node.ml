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

let () = SadmanOutput.register "Node" "$Revision: 2871 $"
let infinity = float_of_int max_int

let debug = false
let debug_exclude = false
let debug_sets = false
let debug_set_cost = false
let odebug = Status.user_message Status.Information

let likelihood_error = 
    "Likelihood not enabled: download different binary or contact mailing list" 

module IntSet = All_sets.Integers

exception Illegal_argument of string
type 'a r = {
    preliminary : 'a;
    final : 'a;
    cost : float;
    sum_cost : float;
    weight : float;
    time : float;
}

(** Schemes for origin and loss costs *)
type origin_loss_cost = [
| `Flat of int * int ]                 (** Flat cost for each origin and loss *)

(** Methods for complex terminal alignments.  These variants also contain the
    extra median information appropriate for that set type. *)
type complex_term_method = [
| `Strictly_Same  (** Don't allow recombination between elements of a set:  sets
                      are only used for grouping *)

| `Any_Of of ((int * int * (int * int) list) * float)
(** Allow choosing of any of the elements to use for the median.  The median
    information is:
    - ID of "left" node from which this median was made
    - ID of "right" node from which this median was made
    - List with the same length as the set;  the contents are the indices of the
    element from which that median was made.
*)

(* | `Tree of int * int * float * tree option *)
(** Allow simple origin and loss costs *)
]

let same_ct_method a b =
    match a, b with
    | `Strictly_Same, `Strictly_Same -> true
    | `Any_Of _, `Any_Of _ -> true
    | `Strictly_Same, _ | `Any_Of _, _ -> false

let ct_same_methods = Utl.pairwisep same_ct_method


type 'a css = {
    sid : int;
    set : 'a list;
    smethod : complex_term_method;
}
let empty_css =
    { sid = -1;
      set = [];
      smethod = `Strictly_Same;
    }

(** [ct_verify sa] is a helper debugging function *)
let ct_verify sa =
    assert (sa.sid >= 0);
    match sa.smethod with
    | `Strictly_Same -> ()
    | `Any_Of ((l, r, list), _) ->
          if l <> 0
          then assert (l <> r);
          assert (List.length list = List.length sa.set)

IFDEF USE_LIKELIHOOD THEN
    type ml_rep = MlStaticCS.t r
ELSE
    type ml_rep = unit
END

type cs = 
    | Nonadd8 of NonaddCS8.t r 
    | Nonadd16 of NonaddCS16.t r 
    | Nonadd32 of NonaddCS32.t r 
    | Add of AddCS.t r
    | Sank of SankCS.t r
    | Dynamic of DynamicCS.t r
    | Kolmo of KolmoCS.t r
    | Set of cs css r
    | StaticMl of ml_rep

let rec to_string_ch ch1 =
    match ch1 with
    | Nonadd8 a ->
          ("na8: " ^ NonaddCS8.to_string a.final)
    | Nonadd16 a ->
          ("na16: " ^ NonaddCS16.to_string a.final)
    | Nonadd32 a ->
          ("na32: " ^ NonaddCS32.to_string a.final)
    | Add a ->
          AddCS.to_string a.final
    | Sank a ->
          ("sank: " ^ SankCS.to_string a.final)
    | Dynamic a ->
          "dynamic: " ^ DynamicCS.to_string a.final
    | Set a ->
          let sub = List.map to_string_ch a.final.set in
          let stype = match a.final.smethod with
          | `Strictly_Same -> "same"
          | `Any_Of _ -> "any-of" in
          "set(" ^ stype ^ "): [" ^ (String.concat "; " sub) ^ "]"
    | Kolmo a ->
            ("kolmo: " ^ KolmoCS.to_string a.final)
    | StaticMl a ->
        IFDEF USE_LIKELIHOOD THEN
            ("static ML: " ^ MlStaticCS.to_string a.preliminary)
        ELSE
            failwith likelihood_error
        END
let extract_cost = function
    | Nonadd8 v -> v.cost
    | Nonadd16 v -> v.cost
    | Nonadd32 v -> v.cost
    | Add v -> v.cost
    | Sank v -> v.cost
    | Dynamic v -> v.cost
    | Set v -> v.cost
    | Kolmo v -> v.cost
    | StaticMl v -> 
        IFDEF USE_LIKELIHOOD THEN
            v.cost
        ELSE
            failwith likelihood_error
        END

(** Helper function for recursing into sets *)
let setrec a fn = { a with set = List.map fn a.set }

(** [set_update_cost set] calculates the new cost of the set based on its
    contents, and updates the node's cost field.  Only applicable for
    [`Strictly_Same] sets. *)
let set_update_cost = function
    | Set a ->
          (match a.preliminary.smethod with
           | `Strictly_Same ->
                 let cost = 
                     List.fold_left
                         (fun acc m -> acc +. extract_cost m)
                         0. a.preliminary.set in
                 if debug_set_cost
                 then odebug ("`Strictly_Same cost: " ^ string_of_float cost);
                 Set { a with cost = cost }
           | `Any_Of _ -> assert false
          )
    | _ -> assert false

type exclude = ([`Excluded | `NotExcluded | `Either] * int * int * int) list

type node_data = 
    { 
        characters : cs list;
        total_cost : float;
        node_cost : float;
        taxon_code : int;
        min_child_code : int;
        num_child_edges : int;
        num_height : int;
        num_otus : int;                 (** How many OTUs are a child of this
                                            node.  No longer necessary? *)
        exclude_sets : All_sets.Integers.t list;
        exclude_info : exclude;
        cost_mode : [ `Likelihood | `Parsimony ];
        (** This allows us to count how many taxa from a set are children of the
            given node *)
    }

let recode f n = 
    { n with taxon_code = f n.taxon_code }

(* Static characters only *)
let character_costs node = 
    let rec get_individual_cost acc lst = 
        match lst with
        | h :: t -> 
                let process_nonadd f x = 
                    List.map (fun (a, _, b) -> (`NonAdd, a, b)) (f x.preliminary) in
                let new_costs = 
                    match h with
                    | Nonadd8 x -> 
                            process_nonadd NonaddCS8.to_list x
                    | Nonadd16 x -> process_nonadd NonaddCS16.to_list x
                    | Nonadd32 x -> process_nonadd NonaddCS32.to_list x
                    | Add x -> 
                            let lst = AddCS.to_list_with_cost x.preliminary in
                            List.map (fun (_, _, code, cost) -> `Add, code, cost) lst
                    | Sank x -> 
                            let lst = SankCS.to_list x.preliminary in
                            List.map (fun (a, b) ->
                                `Sank, a, float_of_int (Array.fold_left min max_int b))
                            lst
                    | _ -> []
                in
                get_individual_cost (new_costs @ acc) t
        | [] -> acc
    in
    (get_individual_cost [] node.characters)

let empty cost_mode = {
    characters = [];
    total_cost = 0.0;
    node_cost = 0.0;
    taxon_code = -1;
    min_child_code = 0;
    num_child_edges = 0;
    num_height = 0;
    num_otus = 0;
    exclude_sets = [];
    exclude_info = [];
    cost_mode = cost_mode;
}

let print nd = 
    List.iter 
        (fun a_cs ->
             match a_cs with
             | Dynamic a_dyn ->
                   print_endline "Preliminary state";
                   DynamicCS.print a_dyn.preliminary;
                   print_endline "Final state";
                   DynamicCS.print a_dyn.final
             | _ -> print_endline "Do not print non-dynamic characters"
        ) nd.characters
    
let to_string {characters=chs; total_cost=cost; taxon_code=tax_code} =
    ("[[NODE tax_code=" ^ string_of_int tax_code
     ^ " cost=" ^ string_of_float cost
     ^ " elts: "
     ^ (String.concat "; " (List.map to_string_ch chs))
     ^ "]]")

(*** Helper functions for node data ***)
let cs_prelim_to_final a =
    {a with final = a.preliminary}

let rec prelim_to_final =
    function
        | StaticMl a -> 
            IFDEF USE_LIKELIHOOD THEN        
                StaticMl (cs_prelim_to_final a)
            ELSE
                failwith likelihood_error
            END
        | Nonadd8 a -> Nonadd8 (cs_prelim_to_final a)
        | Nonadd16 a -> Nonadd16 (cs_prelim_to_final a)
        | Nonadd32 a -> Nonadd32 (cs_prelim_to_final a)
        | Add a -> Add (cs_prelim_to_final a)
        | Sank a -> Sank (cs_prelim_to_final a)
        | Dynamic a -> Dynamic (cs_prelim_to_final a)
        | Kolmo a -> Kolmo (cs_prelim_to_final a)
        | Set a ->
              let r = setrec a.preliminary prelim_to_final in
              Set { a with preliminary = r; final = r; }

let all_prelim_to_final ({characters = chars} as node) =
    {node with
         characters = List.map prelim_to_final chars}

let float_close ?(epsilon=0.001) a b =
    let diff = a -. b in
    (abs_float diff) < epsilon

let rec cs_median code anode bnode prev a b = 
    match a, b with 
    | StaticMl ca, StaticMl cb ->
        IFDEF USE_LIKELIHOOD THEN
            assert (ca.weight = cb.weight);
            let catime,cbtime = 
                if code < 0 then (* it's a median *)
                    (
                (* TODO:: during set(iterative), this doens't hold sometimes! *)
                    (* assert( ca.time = cb.time ); *)
                    (* 
                        let () = Printf.printf "t1: %f,\tt2: %f\n" ca.time cb.time
                    *)
                    let half_time = ca.time /. 2.0 in
                    half_time,half_time
                    )
                else
                    ca.time,cb.time
            in
            let median = MlStaticCS.median ca.preliminary cb.preliminary catime cbtime in
            let n_cost = MlStaticCS.root_cost median in
            let res =
                { ca with 
                    preliminary = median;
                    final = median;
                    cost = n_cost *. ca.weight; 
                    sum_cost = n_cost;
                } 
            in
            StaticMl res
        ELSE
            failwith likelihood_error
        END
    | Nonadd8 ca, Nonadd8 cb ->
            assert (ca.weight = cb.weight);
            let prev =
                match prev with
                | None -> None
                | Some (Nonadd8 prev) -> Some (prev.preliminary)
                | _ -> raise (Illegal_argument "cs_median")
            in
            let median = NonaddCS8.median prev ca.preliminary cb.preliminary in
            let cost = NonaddCS8.median_cost median in
            let res = { ca with preliminary = median; 
                            final = median; 
                            cost = ca.weight *. cost;
                            sum_cost = ca.sum_cost +. cb.sum_cost +. cost;
                      } in
            Nonadd8 res
    | Nonadd16 ca, Nonadd16 cb -> 
            assert (ca.weight = cb.weight);
          let prev =
              match prev with
              | None -> None
              | Some (Nonadd16 prev) -> Some (prev.preliminary)
              | _ -> raise (Illegal_argument "cs_median");
          in
          let median = NonaddCS16.median prev ca.preliminary cb.preliminary in
          let cost = NonaddCS16.median_cost median in
          let res = { ca with preliminary = median; 
                      final = median; 
                      sum_cost = ca.sum_cost +. cb.sum_cost +. cost;
                      cost = ca.weight *. cost } in
          Nonadd16 res
    | Nonadd32 ca, Nonadd32 cb -> 
          assert (ca.weight = cb.weight);
          let prev = 
              match prev with
              | None -> None
              | Some (Nonadd32 prev) -> Some (prev.preliminary)
              | _ -> raise (Illegal_argument "cs_median")
          in
          let median = NonaddCS32.median prev ca.preliminary cb.preliminary in
          let cost = NonaddCS32.median_cost median in
          let res = { ca with preliminary = median; 
                      final = median; 
                      sum_cost = ca.sum_cost +. cb.sum_cost +. cost;
                      cost = ca.weight *. cost } in
          Nonadd32 res
    | Add ca, Add cb -> 
          assert (float_close ca.weight cb.weight);
          let old = match prev with
          | Some Add old -> Some old.preliminary
          | None -> None
          | _ -> assert false in
          let median = AddCS.median old ca.preliminary cb.preliminary in
          let cost = AddCS.median_cost median in
          let res = { ca with
              preliminary = median; 
              final = median; 
              sum_cost = ca.sum_cost +. cb.sum_cost +. cost;
              cost = ca.weight *. cost; } 
          in
          Add res
    | Sank ca, Sank cb ->
            assert (ca.weight = cb.weight);
            let prev = 
                match prev with
                | None -> None
                | Some (Sank prev) -> Some (prev.preliminary)
                | _ -> raise (Illegal_argument "cs_median") 
            in
            let median = SankCS.median prev ca.preliminary cb.preliminary in
            let cost = SankCS.distance ca.preliminary cb.preliminary in
            let res = 
                { 
                    ca with preliminary = median;
                    final = median;
                    sum_cost = ca.sum_cost +. cb.sum_cost +. cost;
                    cost = ca.weight *. cost }
            in
            Sank res
    | Dynamic ca, Dynamic cb ->
            assert (ca.weight = cb.weight);
            let ca, cb =
                if anode.min_child_code < bnode.min_child_code then ca, cb
                else cb, ca
            in
            let median = DynamicCS.median code ca.preliminary cb.preliminary in
            let total_cost = DynamicCS.total_cost median in 
            let res = 
                { ca with 
                    preliminary = median;
                    final = median;
                    cost = ca.weight *. total_cost;
                    sum_cost = ca.sum_cost +. cb.sum_cost +. total_cost;
                } 
            in
            Dynamic res
    | Kolmo ca, Kolmo cb ->
            assert (ca.weight = cb.weight);
            let ca, cb =
                if anode.min_child_code < bnode.min_child_code then ca, cb
                else cb, ca
            in
            let median = KolmoCS.median code ca.preliminary cb.preliminary in
            let total_cost = KolmoCS.total_cost median in 
            let res = 
                { ca with 
                    preliminary = median;
                    final = median;
                    cost = ca.weight *. total_cost;
                    sum_cost = ca.sum_cost +. cb.sum_cost +. total_cost;
                } 
            in
            Kolmo res
    | Set ca, Set cb ->
(*           assert (ca.sid = cb.sid); *)
          assert (same_ct_method ca.preliminary.smethod cb.preliminary.smethod);
          assert (ca.weight = cb.weight);
          let l1 = ca.preliminary.set in
          let l2 = cb.preliminary.set in
          begin
              match ca.preliminary.smethod with
              | `Strictly_Same ->
                    (* just recursively apply it to our children, in order *)
                    let res =
                        List.map2
                            (cs_median code anode bnode None) l1 l2 in
                    let result = { cb.preliminary with set = res } in
                    set_update_cost
                        (Set { cb with
                                  preliminary = result;
                                  final = result; })
              | `Any_Of (_, v) ->
                    (* take all the possible medians, keep those of minimum cost
                     * *)
                    let min_cost = ref (2. *. v) in
                    let update_cost m = min_cost := min m !min_cost in
                   let _, medians =
                        List.fold_left
                            (fun (i, list) l1i ->
                                 let _, list = List.fold_left
                                 (fun (j, list) l2i ->
                                      let median = cs_median code anode
                                          bnode None l1i l2i in
                                      let cost = extract_cost median in
                                      update_cost cost;
                                      (succ j,
                                       (cost, median, i, j) :: list))
                                 (0, list)
                                 l2 in
                                 (succ i, list))
                            (0, [])
                            l1 in
                    let medians =
                        List.filter
                            (fun (cost, _, _, _) -> cost = !min_cost)
                            medians in
                    (* Make the proper sets *)
                    let median =
                        { cb.preliminary with
                              set = List.map
                                (fun (_, median, _, _) -> median) medians;
                              smethod =
                                `Any_Of ((anode.taxon_code, bnode.taxon_code,
                                    (List.map (fun (_, _, i, j) -> (i, j))
                                         medians)), v);
                        } in
                    (* Update the cost and return *)
                    let res = Set { cb with
                                        preliminary = median;
                                        final = median;
                                        cost = ca.weight *. 
                                        (if max_float = !min_cost then 0. else
                                            !min_cost);
                                        sum_cost = ca.sum_cost +. cb.sum_cost
                            +. !min_cost;
                                  } in
                    res
          end
    | Nonadd8 _, _ | Nonadd16 _, _| Nonadd32 _, _ | Add _, _ | Sank _, _ 
    | Dynamic _, _ |  Set _, _ | Kolmo _, _ | StaticMl _, _ -> 
            raise (Illegal_argument "cs_median")

let map2 f a b =
    let rec mapper a b acc =
        match a, b with
        | (ha :: ta), (hb :: tb) ->
                mapper ta tb ((f ha hb) :: acc)
        | [], [] -> List.rev acc
        | _ -> raise (Illegal_argument "map2")
    in
    mapper a b []

let map3 f a b c =
    let rec mapper a b c acc =
        match a, b, c with
        | (ha :: ta), (hb :: tb), (hc :: tc) ->
              mapper ta tb tc ((f ha hb hc) :: acc)
        | [], [], [] -> List.rev acc
        | _ -> raise (Illegal_argument "map3")
    in
    mapper a b c []

let map4 f a b c d =
    let rec mapper a b c d acc =
        match a, b, c, d with
        | (ha :: ta), (hb :: tb), (hc :: tc), (hd :: td) ->
                mapper ta tb tc td ((f ha hb hc hd) :: acc)
        | [], [], [], [] -> List.rev acc
        | _ -> 
                raise (Illegal_argument "map4")
    in
    mapper a b c d []

(** [edge iterator gp rent c1 c2]
 * iterates the branch lengths in the character set that are ML characters. *) 
let edge_iterator (gp:node_data option) (c0:node_data) (c1:node_data) (c2:node_data) = 
    (* preliminary: check if we are adjusting a handle *)
    let root_e = match gp with 
        | Some x -> false
        | None -> true
    in
    (* accumulators [pa],[aa],[ba] to iterate the branch length for each element in
     * the character set [p] as parent, and [a], [b] as children. *)
    let rec ei_map p a b pa aa ba = match p,a,b with
        | (StaticMl pml)::ptl, (StaticMl aml)::atl,(StaticMl bml)::btl ->
        IFDEF USE_LIKELIHOOD THEN
            let modf = ref All_sets.Integers.empty in
            (* if root edge then split the branch in half *)
            let amltime,bmltime = if root_e
                then aml.time /. 2.0, bml.time /. 2.0 
                else aml.time,bml.time in

            let mine,pcost,cost,(t1,t2),res = 
                MlStaticCS.readjust None !modf aml.preliminary bml.preliminary
                    pml.preliminary amltime bmltime pml.time in

            (* pull edges back together *)
            let t1,t2 = if root_e then t1 +. t2,t1 +. t2 else t1,t2 in

            let first   = StaticMl { aml with time = t1; }
            and second  = StaticMl { bml with time = t2; }
            and mine    = StaticMl { pml with preliminary=res;final=res;
                                              cost=cost;sum_cost=cost; } in
            ei_map ptl atl btl (mine::pa) (first::aa) (second::ba)
        ELSE 
            let first = StaticMl aml and second = StaticMl bml and mine = StaticMl pml in

            ei_map ptl atl btl (mine::pa) (first::aa) (second::bb)
        END
        (* ignore non-likelihood characters *)
        (* | Nonadd8 | Nonadd16 | Nonadd32 | Add | Sank | Dynamic | Kolmo | Set *)
        | pml::ptl,aml::atl,bml::btl -> ei_map ptl atl btl (pml::pa) (aml::aa) (bml::ba)
        | [],[],[] -> (pa,aa,ba)
        | _ -> failwith "Number of characters is inconsistent"
    in
    let mine,ch1,ch2 = ei_map c0.characters c1.characters c2.characters [] [] [] in
    let mine_cost = List.fold_left 
            (fun x y -> match y with
                | StaticMl a ->
                    IFDEF USE_LIKELIHOOD THEN
                        x +.a.cost
                    ELSE
                        x
                    END
                | _ -> x) 0.0 mine
    in (
        {c0 with characters = mine;
                 total_cost = mine_cost;
                 node_cost  = mine_cost;},
        {c1 with characters = ch1;}, (*childrens cost doesn't change *)
        {c2 with characters = ch2;}
       )

(** [apply_time node1 node2] - new node_data with time associated. Each
 * node should be in the opposite direction, and ca1 should have the proper time
 * data to apply to the ca2 node. Optionally apply a function to the two times.
 * This is currently useful for dividing the times by 1/2 or adding them
 * together after adjusting the root of three dir tree. *)
let apply_time (ca1:node_data) (ca2:node_data) : node_data = 
    let rec at_map todata fromdata acc = match todata,fromdata with
        | (StaticMl curr)::tl,(StaticMl from)::t2 ->
            IFDEF USE_LIKELIHOOD THEN
                let mine = StaticMl { curr with time = from.time } in
                at_map tl t2 (mine::acc) 
            ELSE
                at_map tl t2 (StaticMl curr::acc)
            END
        | curr::tl,from::t2 -> at_map tl t2 (curr::acc)
        | [],[] -> acc
        | _ -> failwith "apply_time: Characters of inconsistent lengths"
    in
    {ca2 with characters = at_map ca1.characters ca2.characters [] }

let rec cs_median_3 pn nn c1n c2n p n c1 c2 =
    match p, n, c1, c2 with
    | StaticMl cp, StaticMl cn, StaticMl cc1, StaticMl cc2 -> (*
            let m = MlStaticCS.median_3 cp.final cn.preliminary cc1.preliminary 
            cc2.preliminary in
            StaticMl { cn with final = m } *)
        n
    | Nonadd8 cp, Nonadd8 cn, Nonadd8 cc1, Nonadd8 cc2 -> 
          let m = NonaddCS8.median_3 cp.final cn.preliminary cc1.preliminary
              cc2.preliminary in
          Nonadd8 { cn with final = m }
    | Nonadd16 cp, Nonadd16 cn, Nonadd16 cc1, Nonadd16 cc2 -> 
          let m = NonaddCS16.median_3 cp.final cn.preliminary cc1.preliminary
              cc2.preliminary in
          Nonadd16 { cn with final = m }
    | Nonadd32 cp, Nonadd32 cn, Nonadd32 cc1, Nonadd32 cc2 -> 
          let m = NonaddCS32.median_3 cp.final cn.preliminary cc1.preliminary
              cc2.preliminary in
          Nonadd32 { cn with final = m }
    | Add cp, Add cn, Add cc1, Add cc2 -> 
          let m = AddCS.median_3 cp.final cn.preliminary cc1.preliminary
              cc2.preliminary in
          Add { cn with final = m }
    | Sank cp, Sank cn, Sank cc1, Sank cc2 ->
          let m = SankCS.median_3 cp.final cn.preliminary cc1.preliminary
              cc2.preliminary in
          Sank { cn with final = m }
    | Dynamic cp, Dynamic cn, Dynamic cc1, Dynamic cc2 ->
          let m = 
              DynamicCS.median_3 cp.final cn.preliminary cc1.preliminary
              cc2.preliminary
          in
          Dynamic { cn with final = m }
    | Kolmo cp, Kolmo cn, Kolmo cc1, Kolmo cc2 ->
          let m = 
              KolmoCS.median_3 cp.final cn.preliminary cc1.preliminary
              cc2.preliminary
          in
          Kolmo { cn with final = m }
    | Set cp, Set cn, Set cc1, Set cc2 ->
          (match cn.preliminary.smethod with
           | `Strictly_Same ->
                 (* Apply in order to our children *)
                 let res = map4 (cs_median_3 pn nn c1n c2n)
                     cp.preliminary.set cn.preliminary.set
                     cc1.preliminary.set cc2.preliminary.set in
                 (* Note: if anything else updates the cost, the cost should
                  * also be recalculated here.  (But this shouldn't happen in
                  * the 3-median...) *)
                 let res = { cn.preliminary with set = res } in
                 Set { cn with preliminary = res; final = res; }
           | `Any_Of ((froml, fromr, medians), v) ->
                 (* Get the "from" value for a node *)
                 let get_from {smethod=m} = match m with
                 | `Any_Of (l, _) -> l
                 | _ -> assert false in (* Sets should be compatible *)

                 (* Get the from list specifically *)
                 let get_from_list n =
                     let (_, _, l) = get_from n in l in
                 (* Am I the left node? *)
                     let (parent_left, _, _) = get_from cp.final in
                 let is_left =
                     nn.taxon_code = parent_left in
                     let (_, parent_right, _) = get_from cp.final in
                 let is_right =
                     nn.taxon_code = parent_right in
                 assert (is_left || is_right);

                 (* Make a list: parent median and which of my children it came
                    from *)
                 let set = map2
                     (fun median (l, r) ->
                          if is_left
                          then (median, l)
                          else (median, r))
                     cp.final.set (get_from_list cp.final) in
                 let set, from =
                     match medians with
                     (* we're a leaf node (OTU); this is a special case for
                        3-median *)
                     | [] ->
                           (List.map
                                (fun (parent, my_index) ->
                                     let median = List.nth cn.preliminary.set  my_index in
                                     let left, right = median, median in
                                     cs_median_3 pn nn c1n c2n
                                         parent median left right)
                                set,
                            [])
                     (* real 3-median *)
                     | _ ->
                           List.split (List.map
                               (fun (parent, my_index) ->
                                    let median = List.nth cn.preliminary.set  my_index in
                                    let (l, r) = List.nth medians my_index      in
                                    let left =   List.nth cc1.preliminary.set l in
                                    let right =  List.nth cc2.preliminary.set r in
                                    (cs_median_3 pn nn c1n c2n
                                         parent median left right,
                                     (l, r))
                               )
                               set) in
                 let res = { cn.preliminary with
                                 set = set;
                                 smethod = `Any_Of ((froml, fromr, from), v);
                           } in
                 Set { cn with final = res }
          )
    | Nonadd8 _, _, _, _ | Nonadd16 _, _, _, _ | Nonadd32 _, _, _, _ 
    | Add _, _, _, _ | Sank _, _, _, _ | Dynamic _, _, _, _  
    | Set _, _, _, _ | Kolmo _, _, _, _ | StaticMl _, _, _ ,_ ->
          raise (Illegal_argument "cs_median_3")

let new_node_stats a b =
    let num_child_edges =
        a.num_child_edges + b.num_child_edges + 2 in
    let num_height =
        (max a.num_height b.num_height) + 1 in
    (num_child_edges, num_height)
    

let excludes_median a b =
    let exclude_median = function
        | `Excluded, `Excluded -> 0, `Excluded
        | `NotExcluded, `NotExcluded -> 0, `NotExcluded
        | `Either, `Either -> 0, `Either
        | `Either, a -> 0, a
        | a, `Either -> 0, a
        | `Excluded, `NotExcluded
        | `NotExcluded, `Excluded -> 1, `Either in

    List.map2
        (fun (s1, c1, card, count1) (s2, c2, _, count2) ->
             let cost, med = exclude_median (s1, s2) in
             (med, c1 + c2 + cost, card, count1 + count2))
        a.exclude_info
        b.exclude_info

let has_excluded list =
    let excluded (_, cost, size, count) =
        cost = 1 && count = size in
    List.exists excluded list

let get_characters_cost chars =
    List.fold_left (fun a b -> a +. (extract_cost b)) 0. chars 

let total_cost _ a =
    match a.cost_mode with
    | `Likelihood -> a.node_cost
    | `Parsimony -> a.total_cost

let root_cost root =
    let adder acc character = 
        match character with
        | Kolmo v -> acc +. KolmoCS.root_cost v.preliminary
        | _ -> acc
    in
    List.fold_left adder (total_cost root root) root.characters

let get_characters_of_type map t node =
    List.map map 
    (List.filter (function
        | Nonadd8 _ when t = `NonAdd8 -> true
        | Nonadd16 _ when t = `NonAdd8 -> true
        | Nonadd32 _ when t = `NonAdd8 -> true
        | Add _ when t = `Add -> true
        | Sank _ when t = `Sank -> true
        | StaticMl _ when t = `StaticML -> true
        | Dynamic _ when t = `Dynamic -> true
        | Set _ when t = `Set -> true
        | _ -> false) node.characters)

let get_nonadd_8 = 
    get_characters_of_type 
    (function Nonadd8 x -> x.preliminary, x.final | _ -> assert false) 
    `NonAdd8 
let get_nonadd_16 = 
    get_characters_of_type 
    (function Nonadd16 x -> x.preliminary, x.final | _ -> assert false)
    `NonAdd16
let get_nonadd_32 = 
    get_characters_of_type 
    (function Nonadd32 x -> x.preliminary, x.final | _ -> assert false)
    `NonAdd32 
let get_add = 
    get_characters_of_type 
    (function Add x -> x.preliminary, x.final | _ -> assert false)
    `Add
let get_sank = 
    get_characters_of_type 
    (function Sank x -> x.preliminary, x.final | _ -> assert false)
    `Sank
let get_dynamic = 
    get_characters_of_type 
    (function Dynamic x -> x.preliminary, x.final | _ -> assert false)
    `Dynamic

IFDEF USE_LIKELIHOOD THEN
    let get_mlstatic = 
        get_characters_of_type 
            (function StaticMl x -> x.preliminary, x.final | _ -> assert false)
            `StaticML
ELSE
    let get_mlstatic x = []
END

let get_set = 
    get_characters_of_type 
    (function Set x -> x.preliminary, x.final | _ -> assert false)
    `Set

let median_counter = ref (-1)

let median code old a b =
    (* the code is negative if we are calculating on an edge *)
    let code = 
        match code with
        | Some code -> code
        | None -> 
                decr median_counter;
                !median_counter
    in
    (*
    Printf.printf "Code of median is %d with children %d and %d\n%!" code 
    a.taxon_code b.taxon_code;
    *)
    let new_characters =
        match old with
        | None -> 
                map2 (cs_median code a b None) 
                a.characters b.characters
        | Some c ->
                    map3 (fun x -> cs_median code a b (Some x))
                    c.characters
                    a.characters
                    b.characters
    and children_cost = a.total_cost +. b.total_cost in
    let node_cost = get_characters_cost new_characters in
    let total_cost = 
        assert (a.cost_mode = b.cost_mode);
        match a.cost_mode with
        | `Likelihood -> node_cost
        | `Parsimony -> node_cost +. children_cost in
    let num_child_edges, num_height = new_node_stats a b in
    let exclude_info = excludes_median a b in
    let excluded = has_excluded exclude_info in
    { 
        characters = new_characters;
        total_cost =
            if excluded
            then infinity
            else total_cost;
        node_cost = node_cost;
        taxon_code = code;
        min_child_code = min a.min_child_code b.min_child_code;
        num_child_edges = num_child_edges;
        num_height = num_height;
        num_otus = a.num_otus + b.num_otus;
        exclude_sets = a.exclude_sets;
        exclude_info = exclude_info;
        cost_mode = a.cost_mode;
    }

(* let median_final old a b = *)
(*     let old = match old with *)
(*     | None -> None *)
(*     | Some a -> Some (all_prelim_to_final a) in *)
(*     median old (all_prelim_to_final a) (all_prelim_to_final b) *)

let median_3 p n c1 c2 = 
    let new_characters = 
        map4 (cs_median_3 p n c1 c2) 
        p.characters n.characters c1.characters c2.characters
    in
    { n with characters = new_characters }

let median_no_cost median =
    { median with total_cost = 0.0; }
    
let median_self_cost a = a.node_cost

let update_leaf n =
    let exclude_info =
        List.map
            (fun set ->
                 let card = All_sets.Integers.cardinal set in
                 if All_sets.Integers.mem n.taxon_code set
                 then (`Excluded, 0, card, 1)
                 else (`NotExcluded, 0, card, 0))
            n.exclude_sets in
    { n with
          total_cost = 0.;
          exclude_info = exclude_info;
    }

let node_height {num_height = h} = h
let node_child_edges {num_child_edges = c} = c

let get_code {taxon_code=taxcode} = taxcode

let get_cost_mode a = a.cost_mode

let compare_data_final {characters=chs1} {characters=chs2} =
    let rec compare_two ch1 ch2 =
        match ch1, ch2 with
        | Nonadd8 a, Nonadd8 b ->
              NonaddCS8.compare_data a.final b.final
        | Nonadd16 a, Nonadd16 b ->
              NonaddCS16.compare_data a.final b.final
        | Nonadd32 a, Nonadd32 b ->
              NonaddCS32.compare_data a.final b.final
        | Add a, Add b ->
              AddCS.compare_data a.final b.final
        | Sank a, Sank b ->
              SankCS.compare_data a.final b.final
        | Dynamic a, Dynamic b ->
              DynamicCS.compare_data a.final b.final
        | Kolmo a, Kolmo b ->
              KolmoCS.compare_data a.final b.final
        | StaticMl a, StaticMl b ->
            (* preliminary == final *)
            IFDEF USE_LIKELIHOOD THEN
                MlStaticCS.compare_data a.preliminary b.preliminary
            ELSE
                failwith likelihood_error  
            END
        | Set { final = { set = a } }, Set { final = { set = b } } ->
              (* Temporary solution.  The problem is: our `Any_Of nodes need to
                 know from which nodes they were made as a median.  Even if two
                 sets are equal, this is information that we have to update.
                 As a result, claiming the equality of two Set types breaks
                 their processing.  Returning -1 forces their re-update on the
                 uppass. *)
              -1
(*               let cmp = compare (List.length a) (List.length b) in *)
(*               if cmp <> 0 *)
(*               then cmp *)
(*               else compare_lists a b *)
        | _ -> failwith "Incompatible characters (1)"
    and compare_lists chs1 chs2 =
        match chs1, chs2 with
        | ch1 :: chs1, ch2 :: chs2 ->
              let c = compare_two ch1 ch2 in
              if c <> 0 then c
              else compare_lists chs1 chs2
        | [], [] -> 0
        | _ -> failwith "Incompatible characters (2)" in
    compare_lists chs1 chs2


let compare_data_preliminary {characters=chs1} {characters=chs2} =
    let rec compare_two ?(fail=true) ch1 ch2 =
        match ch1, ch2 with
        | Nonadd8 a, Nonadd8 b ->
              NonaddCS8.compare_data a.preliminary b.preliminary
        | Nonadd16 a, Nonadd16 b ->
              NonaddCS16.compare_data a.preliminary b.preliminary
        | Nonadd32 a, Nonadd32 b ->
              NonaddCS32.compare_data a.preliminary b.preliminary
        | Add a, Add b ->
              AddCS.compare_data a.preliminary b.preliminary
        | Sank a, Sank b ->
              SankCS.compare_data a.preliminary b.preliminary
        | Dynamic a, Dynamic b ->
              DynamicCS.compare_data a.preliminary b.preliminary
        | Kolmo a, Kolmo b ->
              KolmoCS.compare_data a.preliminary b.preliminary
        | StaticMl a, StaticMl b ->
            IFDEF USE_LIKELIHOOD THEN
                MlStaticCS.compare_data a.preliminary b.preliminary
            ELSE
                failwith likelihood_error
            END
        | Set { preliminary = { set = a } }, Set {preliminary = { set = b }} ->
              (* See above... *)
              -1
(*               compare_lists ~fail:false a b *)
        | _ -> failwith "Incompatible characters (3)"
    and compare_lists ?(fail=true) chs1 chs2 =
        match chs1, chs2 with
        | ch1 :: chs1, ch2 :: chs2 ->
              let c = compare_two ~fail ch1 ch2 in
              if c <> 0 then c
              else compare_lists ~fail chs1 chs2
        | [], [] -> 0
        | [], _ ->
              if fail
              then failwith "Incompatible characters (4)"
              else -1
        | _, [] ->
              if fail
              then failwith "Incompatible characters (4)"
              else 1
    in
    compare_lists chs1 chs2

let edge_distance clas nodea nodeb =
    (* This function assumes that nodea and nodeb come from a valid tree and
    * nodea is the parent of nodeb *)
    let rec distance_two ch1 ch2 =
        match ch1, ch2 with
        | Nonadd8 a, Nonadd8 b ->
                (match clas with
                | `Static | `Any -> 
                        a.weight *. NonaddCS8.distance a.final b.final
                | `Dynamic -> 0.)
        | Nonadd16 a, Nonadd16 b ->
                (match clas with
                | `Static | `Any -> 
                        a.weight *. NonaddCS16.distance a.final b.final
                | `Dynamic -> 0.)
        | Nonadd32 a, Nonadd32 b ->
                (match clas with
                | `Static | `Any -> 
                        a.weight *. NonaddCS32.distance a.final b.final
                | `Dynamic -> 0.)
        | Add a, Add b ->
                (match clas with
                | `Static | `Any -> 
                        a.weight *. AddCS.distance a.final b.final
                | `Dynamic -> 0.)
        | Sank a, Sank b ->
                (match clas with
                | `Static | `Any -> 
                        a.weight *. SankCS.distance a.final b.final
                | `Dynamic -> 0.)
        | Dynamic a, Dynamic b ->
                (match clas with
                | `Dynamic | `Any -> 
                        (* Observe that we REQUIRE the single assignment for
                        * this collpse to be correct. *)
                        let d = 
                            DynamicCS.distance a.preliminary b.preliminary 
                        in
                        a.weight *. d
                | `Static -> 0.)
        | Kolmo a, Kolmo b ->
              a.weight *. KolmoCS.tabu_distance a.final b.final
        | StaticMl a, StaticMl b ->
            IFDEF USE_LIKELIHOOD THEN
                a.weight *. (MlStaticCS.distance a.final b.final a.time b.time)
            ELSE
                failwith likelihood_error
            END
        | Set a, Set b ->
              (match a.final.smethod with
               | `Strictly_Same ->
                     distance_lists a.final.set b.final.set 0.
               | `Any_Of _ ->
                     (* unf. we just take the full median and check the distance
                     *)
                     decr median_counter;
                     let m = cs_median !median_counter nodea nodeb None ch1 ch2 in
                     extract_cost m)
        | _ -> failwith "Incompatible characters (5)"
    and distance_lists chs1 chs2 acc =
        match chs1, chs2 with
        | ch1 :: chs1, ch2 :: chs2 ->
              distance_lists chs1 chs2 (acc +. distance_two ch1 ch2)
        | [], [] -> acc
        | _ -> failwith "Incompatible characters (6)" in
    distance_lists nodea.characters nodeb.characters 0.

let has_to_single : [ `Add | `Annchrom | `Breakinv | `Chrom | `Genome | `Kolmo
| `Nonadd | `Sank | `Seq | `StaticMl ] list = [`Seq ; `Chrom; `Annchrom; `Breakinv]

let distance_of_type ?(para=None) ?(parb=None) t
    ({characters=chs1} as nodea) ({characters=chs2} as nodeb) =
    let has_t x = List.exists (fun z -> z = x) t
    and filter_dynamic res x = 
        match x with
        | `Seq | `Breakinv | `Chrom | `Annchrom | `Genome as x -> x :: res
        | _ -> res 
    in
    let has_nonadd = has_t `Nonadd
    and has_add = has_t `Add
    and has_sank = has_t `Sank
    and has_staticml = has_t `StaticMl
    and dy_t = List.fold_left filter_dynamic [] t 
    and has_kolmo = has_t `Kolmo in
    let rec distance_two ch1 ch2 =
        match ch1, ch2 with
        | Nonadd8 a, Nonadd8 b when has_nonadd ->
              a.weight *. NonaddCS8.distance a.final b.final
        | Nonadd16 a, Nonadd16 b when has_nonadd ->
              a.weight *. NonaddCS16.distance a.final b.final
        | Nonadd32 a, Nonadd32 b when has_nonadd ->
              a.weight *. NonaddCS32.distance a.final b.final
        | Add a, Add b when has_add ->
              a.weight *. AddCS.distance a.final b.final
        | Sank a, Sank b when has_sank ->
              a.weight *. SankCS.distance a.final b.final
        | Dynamic a, Dynamic b ->
              a.weight *. DynamicCS.distance_of_type dy_t a.final b.final
        | Kolmo a, Kolmo b when has_kolmo ->
              a.weight *. KolmoCS.distance a.final b.final
        | StaticMl a, StaticMl b when has_staticml ->
            IFDEF USE_LIKELIHOOD THEN
                0.0
                (* a.weight *. MlStaticCS.distance a.final b.final a.time b.time *)
            ELSE
                failwith likelihood_error
            END
        | Set a, Set b ->
              (match a.final.smethod with
               | `Strictly_Same ->
                     distance_lists a.final.set b.final.set 0.
               | `Any_Of _ ->
                     (* unf. we just take the full median and check the distance
                     *)
                     decr median_counter;
                     let m = cs_median !median_counter nodea nodeb None ch1 ch2 in
                     extract_cost m)
        | _ -> 0.0
    and distance_lists chs1 chs2 acc =
        match chs1, chs2 with
        | ch1 :: chs1, ch2 :: chs2 ->
              distance_lists chs1 chs2 (acc +. distance_two ch1 ch2)
        | [], [] -> acc
        | _ -> failwith "Incompatible characters (6)" in
    distance_lists chs1 chs2 0.


let distance ?(para=None) ?(parb=None) 
    ({characters=chs1} as nodea) ({characters=chs2} as nodeb) =
    let rec distance_two ch1 ch2 =
        match ch1, ch2 with
        | Nonadd8 a, Nonadd8 b ->
              a.weight *. NonaddCS8.distance a.final b.final
        | Nonadd16 a, Nonadd16 b ->
              a.weight *. NonaddCS16.distance a.final b.final
        | Nonadd32 a, Nonadd32 b ->
              a.weight *. NonaddCS32.distance a.final b.final
        | Add a, Add b ->
              a.weight *. AddCS.distance a.final b.final
        | Sank a, Sank b ->
              a.weight *. SankCS.distance a.final b.final
        | Dynamic a, Dynamic b ->
              a.weight *. DynamicCS.distance a.final b.final
        | Kolmo a, Kolmo b ->
              a.weight *. KolmoCS.distance a.final b.final
        | StaticMl a, StaticMl b ->
            IFDEF USE_LIKELIHOOD THEN
                a.weight *. MlStaticCS.distance a.final b.final a.time b.time
            ELSE
                failwith likelihood_error
            END
        | Set a, Set b ->
              (match a.final.smethod with
               | `Strictly_Same ->
                     distance_lists a.final.set b.final.set 0.
               | `Any_Of _ ->
                     (* unf. we just take the full median and check the distance
                     *)
                     decr median_counter;
                     let m = cs_median !median_counter nodea nodeb None ch1 ch2 in
                     extract_cost m)
        | _ -> failwith "Incompatible characters (5)"
    and distance_lists chs1 chs2 acc =
        match chs1, chs2 with
        | ch1 :: chs1, ch2 :: chs2 ->
              distance_lists chs1 chs2 (acc +. distance_two ch1 ch2)
        | [], [] -> acc
        | _ -> failwith "Incompatible characters (6)" in
    distance_lists chs1 chs2 0.

(* Calculates the cost of joining the node [n] between [a] and [b] in a tree *)
(* [a] must be the parent (ancestor) of [b] *)
let dist_2 minimum_delta n a b =
    let rec ch_dist delta_left n a' b' =
        match n, a', b' with
        | StaticMl n, StaticMl a, StaticMl b ->
            IFDEF USE_LIKELIHOOD THEN
                n.weight *. MlStaticCS.dist_2 n.final a.final b.final n.time
                a.time b.time None
            ELSE
                failwith likelihood_error
            END
        | Nonadd8 n, Nonadd8 a, Nonadd8 b ->
              n.weight *. NonaddCS8.dist_2 n.final a.final b.final
        | Nonadd16 n, Nonadd16 a, Nonadd16 b ->
              n.weight *. NonaddCS16.dist_2 n.final a.final b.final
        | Nonadd32 n, Nonadd32 a, Nonadd32 b ->
              n.weight *. NonaddCS32.dist_2 n.final a.final b.final
        | Add n, Add a, Add b ->
              n.weight *. AddCS.distance_2 n.final a.final b.final
        | Sank n, Sank a, Sank b ->
              n.weight *. SankCS.dist_2 n.final a.final b.final
        | Dynamic nc, Dynamic ac, Dynamic bc ->
                let ac, bc = 
                    if a.min_child_code < b.min_child_code then ac, bc
                    else bc, ac
                in
                let d  = (DynamicCS.dist_2 delta_left nc.final ac.final bc.final) in
                nc.weight *. d
        | Kolmo nc, Kolmo ac, Kolmo bc ->
                let ac, bc = 
                    if a.min_child_code < b.min_child_code then ac, bc
                    else bc, ac
                in
                nc.weight *. (KolmoCS.dist_2 delta_left nc.final ac.final bc.final)
        | Set sn, Set sa, Set sb ->
              let b_code = b.taxon_code in
              (match sn.final.smethod with
               | `Strictly_Same ->
                     (* simply map over all characters and sum *)
                     (* r delta_left total_delta list_n list_a list_b *)
                     let rec r dl td ln la lb = match ln, la, lb with
                     | n :: ln, a :: la, b :: lb ->
                           let d = ch_dist dl n a b in
                           if d > dl
                           then d
                           else 
                               let td = td +. d in
                               let dl = dl -. d in
                               r dl td ln la lb
                     | [], [], [] -> td
                     | _ -> assert false in
                     r delta_left 0.
                         sn.final.set
                         sa.final.set
                         sb.final.set
               | `Any_Of ((_, _, nlist), _) ->
                     (* find the combination with the least cost *)
                     (* (since all final states can be in a min-cost tree) *)
                     let cmin = ref infinity in
                     let record n = if n < !cmin then cmin := n in
                     let aleft, aright, afromlist =
                         match sa.final.smethod with
                         | `Any_Of ((l, r, list), _) -> l, r, list
                         | _ -> assert false in
                     let alist = sa.final.set in
                     let is_left = b_code = aleft in
                     let is_right = b_code = aright in
                     assert (is_left || is_right);
                     let blist =
                         let fn (l, r) =
                             if is_left
                             then l
                             else r in
                         List.map fn afromlist in
                     let blist =
                         List.map (List.nth sb.preliminary.set) blist in
                       (* TODO: I am not updating this delta_left properly? *)
                     List.iter
                         (fun n -> List.iter2
                              (fun a b ->
                                   record (ch_dist delta_left n a b))
                              alist
                              blist)
                         sn.final.set;
                     !cmin)
        | Nonadd8 _, _, _ | Nonadd16 _, _, _ | Nonadd32 _, _, _
        | Add _, _, _ | Sank _, _, _ | Dynamic _, _, _ 
        | Set _, _, _ | Kolmo _, _, _ | StaticMl _, _, _ ->
                (* These are explicitly left so that modifying code is easier,
                 * the compiler guides the changes *)
              raise (Illegal_argument "dist_2")
    in
    let rec chars acc lst =
        if acc > minimum_delta then acc
        else
            match lst with
            | n :: ncs, a :: acs, b :: bcs ->
                    let max_delta = minimum_delta -. acc in
                  chars (acc +. ch_dist max_delta n a b) (ncs, acs, bcs)
            | [], [], [] -> acc
            | _ -> raise (Illegal_argument "dist_2_chars")
    in
    let excludes = excludes_median n a in
    if has_excluded excludes
    then infinity
    else chars 0. (n.characters, a.characters, b.characters)

let extract_stat = function
    | (Data.Stat (a, b), _) -> (b, a)
    | _ -> raise (Illegal_argument "extract_stat")

let extract_dynamic data dyna tcode = 
    match dyna with 
    | (Data.Dyna (dyna_code, chrom_data), _) -> 
          let specs = Hashtbl.find data.Data.character_specs dyna_code in  
          let specs, weight =
              match specs with 
              | Data.Dynamic d -> d, d.Data.weight
              | _ -> failwith "extract dynamic"
          in 
          let dyna =  
              let num_taxa = Data.number_of_taxa data in
              DynamicCS.of_array specs [|(chrom_data, dyna_code)|] 
              dyna_code tcode num_taxa 
          in 
          { preliminary = dyna; 
            final = dyna; 
            cost = 0.; 
            weight = weight;
            sum_cost = 0.;
            time = 0.
          } 
    | Data.Stat (code, _), _ ->
          raise (Illegal_argument ("Stat" ^ (string_of_int code))) 

let extract_kolmo data kolmo tcode = 
    match kolmo with
    | (Data.Dyna (chcode, chrom_data), _) -> 
          let specs = Hashtbl.find data.Data.character_specs chcode in  
          let kspec, weight =
              match specs with 
              | Data.Kolmogorov d -> d, d.Data.dhs.Data.weight
              | _ -> failwith "extract dynamic"
          in 
          let dyna =  
              let num_taxa = Data.number_of_taxa data in
              KolmoCS.of_array kspec [|(chrom_data, chcode)|] chcode tcode 
              num_taxa
          in 
          { preliminary = dyna; 
            final = dyna; 
            cost = 0.; 
            weight = weight;
            sum_cost = 0.;
            time = 0.;
          } 
    | Data.Stat (code, _), _ ->
          raise (Illegal_argument ("Stat" ^ (string_of_int code))) 
              

type ms = All_sets.Integers.t

module OrderedLists = struct
    type t = float * (BitSet.t list)
    let rec single_compare a b =
        match a, b with
        | ha :: ta, hb :: tb -> 
                let res = ha - hb in
                if res = 0 then single_compare ta tb
                else res
        | [], [] -> 0
        | [], _ -> -1
        | _, [] -> 1


    let rec list_compare a b = 
        match a, b with
        | ha :: ta, hb :: tb ->
                let res = BitSet.compare ha hb in
                if 0 = res then compare ta tb
                else res
        | [], [] -> 0
        | [], _ -> -1
        | _, [] -> 1

    let compare (aw, a) (bw, b) =
        match Pervasives.compare aw bw  with
        | 0 -> list_compare a b
        | x -> x
end

module SetLists = Map.Make (OrderedLists)

let debug_profile_memory = false

let current_snapshot x = 
    if debug_profile_memory then MemProfiler.current_snapshot x
    else ()

let to_bitset size lst = 
    let set  = BitSet.create (size * (List.length lst)) in
    let cnt = ref (-1) in
    List.iter (fun lst ->
        incr cnt;
        List.iter (fun x ->
            BitSet.set set (x + (!cnt * 32))) lst) lst;
    set

let bitset_table  = Hashtbl.create 1667 

let collapse size characters all_static =
    (*
    let simplify item = 
        let hash = Hashtbl.create 255 
        and cur_counter = ref 1 in
        let simplify item = 
            let rec simplify item =
                if Hashtbl.mem hash item then Hashtbl.find hash item 
                else if Hashtbl.mem hash !cur_counter then begin
                    incr cur_counter;
                    simplify item
                end else begin
                    Hashtbl.add hash !cur_counter item;
                    Hashtbl.add hash item !cur_counter;
                    incr cur_counter;
                    (!cur_counter) - 1
                end
            in
            List.map simplify item
        in
        List.map simplify item
    in
    *)
    let process_all_lists lists = 
        current_snapshot "This is before lists";
        let lists = 
            List.fold_left (fun acc x ->
                let (code, weight, lst) = characters x in
                let lst = 
                    List.map (function
                    | `List x -> 
                            if Hashtbl.mem bitset_table x then 
                                Hashtbl.find bitset_table x
                            else
                                let set = BitSet.create 31 in
                                let () = List.iter (BitSet.set set) x in
                                let () = Hashtbl.add bitset_table x set in
                                set
                    | `Bits x -> x) lst 
                in
                (*
                let lst = simplify lst in
                let lst = to_bitset size lst in
                *)
                if SetLists.mem (weight, lst) acc then
                    let code, nweight = SetLists.find (weight, lst) acc in
                    SetLists.add (weight, lst) (code, (weight +. nweight)) acc
                else SetLists.add (weight, lst) (code, weight) acc) SetLists.empty 
                all_static
        in
        current_snapshot "This is before elements";
        SetLists.fold (fun _ it acc -> it :: acc)
        lists []
    in
    let res = process_all_lists characters in
    current_snapshot "This is after lists";
    res

let classify size chars data =
    let all_static = 
        Hashtbl.fold (fun code spec acc ->
            match spec with
            | Data.Static spec ->
                    (match spec.Parser.SC.st_type with
                    | Parser.SC.STUnordered ->
                        (code, spec) :: acc
                    | _ -> acc)
            | _ -> acc) data.Data.character_specs []
    in
    let taxa (code, spec) = 
        let weight, observed = 
            match spec.Parser.SC.st_type with
            | Parser.SC.STUnordered ->
                    spec.Parser.SC.st_weight, spec.Parser.SC.st_observed
            | _ -> assert false
        in
        Hashtbl.fold (fun _ taxon_chars acc ->
            let lst =
                try
                    match Hashtbl.find taxon_chars code with
                    | (Data.Stat (c, (Some v)), `Specified) -> (c, weight, v)
                    | (Data.Stat (c, v), _) -> (c, weight, `List observed)
                    | _ -> failwith "Impossible 2?"
                with
                | Not_found -> (code, weight, `List observed)
            in
            lst :: acc) data.Data.taxon_characters []
    in
    current_snapshot "Generating characters";
    let characters = 
        let add_taxon_to_accumulator acc (_, _, v) = v :: acc in
        let reshape chars (a, b, _) = (a, b, chars) in
        let chars item =
            match taxa item with
            | ((_, _, x) as h) :: t -> 
                    let chars = List.fold_left add_taxon_to_accumulator [x] t in
                    reshape chars h
            | [] -> failwith "Nothing?"
        in
        collapse size chars all_static
    in
    current_snapshot "Final fold characters";
    let r = List.fold_left (fun acc (a, b) ->
        All_sets.IntegerMap.add a b acc) All_sets.IntegerMap.empty
        characters
    in
    current_snapshot "Done fold";
    r

let classify size doit chars data =
    if doit then Some (classify size chars data)
    else None

let generate_taxon do_classify (laddcode : ms) (lnadd8code : ms) 
    (lnadd16code : ms) (lnadd32code : ms) (lnadd33code : ms) lsankcode dynamics 
    kolmogorov static_ml data cost_mode =
        let add_character =  Data.add_character_spec 
        and set = Data.Set 
        and data = ref data in
        let character_code_gen () =
            incr (!data).Data.character_code_gen;
            !((!data).Data.character_code_gen)
        in
        let cg () = 
            let code = character_code_gen () in
            data := add_character set code !data;
            code
        in
        let group_in_weights weights codes = 
            let get_weight c = 
                match weights with
                | None -> Data.get_weight c !data
                | Some v ->
                        let a = All_sets.IntegerMap.find c v in
                        a
            in
            let table = Hashtbl.create 1667 in
            let weights = 
                All_sets.Integers.fold 
                    (fun x acc -> 
                        try 
                            let w = get_weight x in
                            Hashtbl.add table w x; 
                            All_sets.Floats.add w acc
                        with
                        | Not_found -> acc) 
                    codes All_sets.Floats.empty 
            in
            let res = All_sets.Floats.fold 
                (fun x acc -> 
                    let lst = Hashtbl.find_all table x in (x, lst) :: acc) 
                weights []
            in
            List.fold_left
                (fun acc (w, lst) -> 
                    (character_code_gen (), 
                    (List.map (fun x -> w, x) lst)) :: acc)
                [] res
        in
        let group_ml_by_codes lst = 
            let hstbl = Hashtbl.create 97 in
            let set_codes = 
                List.fold_left (fun acc code ->
                    match Hashtbl.find (!data).Data.character_specs code with
                    | Data.Static spec -> 
                            let model = 
                                match spec.Parser.SC.st_type with
                                | Parser.SC.STLikelihood x -> x
                                | _ -> assert false
                            in
                            Hashtbl.add hstbl model.Parser.SC.set_code code;
                            let group = model.Parser.SC.set_code in
                            if All_sets.Integers.mem group acc then
                                acc
                            else All_sets.Integers.add group acc
                    | _ -> assert false)
                All_sets.Integers.empty lst
            in
            List.map (Hashtbl.find_all hstbl) 
            (All_sets.Integers.elements set_codes)
        in
        let nadd8weights = classify 8 do_classify lnadd8code !data
        and nadd16weights = classify 16 do_classify lnadd16code !data
        and nadd32weights = classify 32 do_classify lnadd32code !data in
        let laddcode = group_in_weights None laddcode
        and lnadd8code = group_in_weights nadd8weights lnadd8code
        and lnadd16code = group_in_weights nadd16weights lnadd16code
        and lnadd32code = group_in_weights nadd32weights lnadd32code
        and lstaticmlcode = group_ml_by_codes static_ml
        and lsankcode = List.map (fun x -> cg (), x) lsankcode in
        let add_codes ((_, x) as y) = y, Array.map snd (Array.of_list x) in
        let laddcode = List.map add_codes laddcode 
        and lnadd8code = List.map add_codes lnadd8code
        and lnadd16code = List.map add_codes lnadd16code
        and lnadd32code = List.map add_codes lnadd32code in
        (* We need ways of making empty characters when a character is
           unspecified *)
        let get_static_encoding code =
            let specs =
                Hashtbl.find !data.Data.character_specs code in
            match specs with
            | Data.Static encoding -> encoding
            | Data.Dynamic _
            | Data.Kolmogorov _
            | Data.Set -> failwith "get_static_encoding" in
        let module Enc = Parser.OldHennig.Encoding in
        let gen_add code =
            let enc = get_static_encoding code in
            (Data.Stat (code, Some (`List enc.Parser.SC.st_observed)), 
            `Unknown) 
        in
        let gen_nadd code =
            let enc = get_static_encoding code in
            (Data.Stat (code, Some (`List enc.Parser.SC.st_observed)), `Unknown) 
        in
        let gen_dynamic code =
            let alph = Data.get_alphabet !data code in
            (* print_endline ("adding sequence with code " ^ string_of_int code); *)
            let empty_seq = Data.get_empty_seq alph in
            let chrom_data = Data.set_dyna_data [|empty_seq|] in 
            (Data.Dyna (code, chrom_data), `Unknown)
        in
        let gen_sank code =
            (* print_endline ("adding sankoff with code " ^ string_of_int code);
            * *)
            let specs = Hashtbl.find !data.Data.character_specs code in
            let states = 
                match specs with
                | Data.Static enc -> `List enc.Parser.SC.st_observed
                | _ -> assert false 
            in
            (Data.Stat (code, Some states), `Unknown)
        in
        current_snapshot "Done";
        !data, 
        fun tcode acc ->
            current_snapshot "Generating taxon";
            let tcharacters = Hashtbl.find !data.Data.taxon_characters tcode in
            let get_character_with_code_n_weight gen_new (w, acc, cnt) (weight, code) = 
                try 
                    weight, (Hashtbl.find tcharacters code) :: acc, cnt + 1
                with
                | Not_found -> weight, (gen_new code) :: acc, cnt
            in
            let get_character_with_code gen_new acc code = 
                try 
                    (Hashtbl.find tcharacters code) :: acc
                with
                | Not_found -> (gen_new code) :: acc
            in
            let addmapper gen_new ((x, y), arr) =
                let a, b, cnt = 
                    List.fold_left (get_character_with_code_n_weight gen_new) 
                    (1., [], 0) y
                in
                x, (a, b), arr
            in
            let ladd_chars    = List.map (addmapper gen_add)  laddcode    
            and lnadd8_chars  = 
                List.map (addmapper gen_nadd) lnadd8code
            and lnadd16_chars = 
                List.map (addmapper gen_nadd) lnadd16code
            and lnadd32_chars = 
                List.map (addmapper gen_nadd) lnadd32code
            and lnadd33_chars = []
            and ldynamic_chars = 
                List.fold_left (get_character_with_code gen_dynamic) [] dynamics
            and lkolmo_chars = 
                List.fold_left (get_character_with_code gen_dynamic) [] kolmogorov
            and lsank_chars =
                List.map 
                (fun (x, y) -> x, List.fold_left (get_character_with_code gen_sank) 
                [] y) 
                lsankcode
            in
            let result = { 
                characters = []; 
                total_cost = 0.;
                node_cost = 0.;
                taxon_code = tcode;
                min_child_code = tcode;
                num_child_edges = 0;
                num_height = 0;
                num_otus = 1;
                exclude_info = [];
                exclude_sets = [];
                cost_mode = cost_mode;
            } in
            let add_characters of_parser buildme result (code, (w, char)) =
                let v = List.map extract_stat char in
                let arr = Array.of_list v in
                let c, _ = of_parser (arr, tcode) code in
                let c = buildme c w in
                { result with characters = c :: result.characters }
            in
            let make_with_w c w =
                { preliminary = c; final = c; cost = 0.; sum_cost = 0.; weight =
                    w; time = 0. }
            in
            let result = 
                List.fold_left 
                (fun acc (a, b, c) -> 
                    add_characters (NonaddCS8.of_parser !data c)
                    (fun c w -> Nonadd8 (make_with_w c w)) acc (a, b))
                result lnadd8_chars
            in
            let result =
                List.fold_left 
                (fun acc (a, b, c) -> 
                    add_characters (NonaddCS16.of_parser !data c)
                    (fun c w -> Nonadd16 (make_with_w c w)) acc (a, b))
                result lnadd16_chars
            in
            let result =
                List.fold_left 
                (fun acc (a, b, c) -> add_characters (NonaddCS32.of_parser !data c)
                (fun c w -> Nonadd32 (make_with_w c w)) acc (a, b))
                result lnadd32_chars
            in
            let result = 
                match lnadd33_chars with
                | _ -> result
            in
            let result = 
                List.fold_left 
                (fun acc (a, b, _) -> add_characters (AddCS.of_parser !data)
                (fun c w -> Add (make_with_w c w)) acc (a, b))
                result ladd_chars
            in
            let result = 
                match ldynamic_chars with
                | [] -> result
                | _ ->
                      let c = 
                          List.map 
                          (fun dyna -> extract_dynamic !data dyna tcode) 
                          ldynamic_chars 
                      in
                      let c : cs list = List.map (fun c -> Dynamic c) c in 
                      { result with characters = c @ result.characters } 
            in
            let result = 
                match lkolmo_chars with
                | [] -> result
                | _ ->
                        let c = 
                            List.map 
                            (fun kolm -> extract_kolmo !data kolm tcode)
                            lkolmo_chars
                        in
                        let c = List.map (fun c -> Kolmo c) c in
                        { result with characters = c @ result.characters }
            in
            let result = 
                let single_lsank_chars_process result (code, lst) =
                    match lst with
                    | [] -> result
                    | _ ->
                            let v = List.map extract_stat lst in
                            let tcm =
                                match v with
                                | (_, code) :: _ -> Data.get_tcm code !data
                                | _ -> failwith "This is impossible"
                            in
                            let arr= Array.of_list v in
                            let c, _ = SankCS.of_parser tcm (arr, tcode) code in
                            let c = Sank { preliminary = c; final = c; cost = 0.;
                                          sum_cost = 0.;
                            weight = 1.; time = 0. } in
                            { result with characters = c :: result.characters }
                in
                List.fold_left single_lsank_chars_process result lsank_chars
            in
            let result = 
                IFDEF USE_LIKELIHOOD THEN
                let single_ml_group result lst =
                    match lst with
                    | [] -> result
                    | h :: t -> (* We do have some characters to add *)
                            let spec = 
                                match Hashtbl.find
                                    (!data).Data.character_specs h with
                                | Data.Static x -> x 
                                | _ -> assert false
                            in
                            let lst = List.map (Hashtbl.find tcharacters) lst in
                            let v = List.map extract_stat lst in
                            let c = 
                                MlStaticCS.of_parser spec (Array.of_list v)
                            in
                            let c = StaticMl { preliminary = c; final = c; cost = 0.;
                                          sum_cost = 0.;
                            weight = 1.; time = 1.5 } in (* TODO: estimate *)
                            { result with characters = c :: result.characters }
                in
                List.fold_left single_ml_group result lstaticmlcode
                ELSE
                result
                END
            in
            let () = current_snapshot "Finished taxon" in
            result :: acc

let node_contents_compare a b =         (* sets? *)
    match a, b with
    | Nonadd8 _ , Nonadd8 _
    | Nonadd16 _, Nonadd16 _
    | Nonadd32 _, Nonadd32 _
    | Add _     , Add _
    | Sank _    , Sank _
    | Dynamic _     , Dynamic _
    | Kolmo _     , Kolmo _
    | StaticMl _, StaticMl _
    | Set _     , Set _      -> 0

    | Nonadd8 _ , _          -> (-1)
    | _         , Nonadd8 _  -> ( 1)
    | Nonadd16 _, _          -> (-1)
    | _         , Nonadd16 _ -> ( 1)
    | Nonadd32 _, _          -> (-1)
    | _         , Nonadd32 _ -> ( 1)
    | Add _     , _          -> (-1)
    | _         , Add _      -> ( 1)
    | Sank _    , _          -> (-1)
    | _         , Sank _     -> ( 1)
    | Dynamic _     , _          -> (-1)
    | _         , Dynamic _      -> ( 1)
    | Kolmo _, _ -> (-1)
    | _, Kolmo _ -> (1)
    | StaticMl _, _ -> (-1)
    | _, StaticMl _ -> (1)


(** [structure_into_sets data nodes] reads the complex terminal structure
    stored in [data.complex_schema] and combines the parsed nodes in [nodes]
    into the proper structure.  Note that this function modifies [data]. *)
let structure_into_sets data (nodes : node_data list) =
    let data' = ref data in
    let add_taxon str =
        let newdata, code = Data.process_taxon_code !data' str "" in
        data' := newdata;
        code in
    let new_sid () =
        let sid = ref 0 in
        fun () -> (incr sid; !sid) 
    in
    let node_to_cs { characters = characters; taxon_code = taxon_code } =
        let set = 
            { sid = taxon_code;
              set = characters;
              smethod = `Strictly_Same; } in
        Set { preliminary = set;
              final = set;
              cost = 0.;
              sum_cost = 0.;
              weight = 1.; time = 0.} 
    in
    let cs_list_to_set new_sid meth css : cs =
        let sid = new_sid () in
        let set =
            { sid = sid;
              set = css;
              smethod = meth; } in
        Set { preliminary = set;
              final = set;
              cost = 0.;
              sum_cost = 0.;
              weight = 1.; time = 0.} 
    in
    let rec group nid new_sid = function
        | Parser.SetGroups.Elt taxon ->
              (try
                   let taxon_id = All_sets.StringMap.find taxon 
                       data.Data.taxon_names in 
                   (try 
                        (node_to_cs 
                             (List.find (fun node -> node.taxon_code = taxon_id) nodes), 
                         taxon) 
                    with Not_found -> 
                        assert false) (* node id should exist in list of nodes *) 
               with Not_found ->
                   failwith ("Complex terminals: Taxon named \"" ^ taxon ^ "\" not found!"))
        | Parser.SetGroups.Set (name, meth, list) ->
              let meth = match meth with
              | Parser.SetGroups.Group -> `Strictly_Same
              | Parser.SetGroups.Any f -> `Any_Of ((nid, nid, []), f)
              | Parser.SetGroups.Tree f -> assert false (* `Tree (nid, nid, f, None) *) in
              let css, _ = List.split (List.map (group nid new_sid) list) in
              (cs_list_to_set new_sid meth css, name) 
    in
    let cost_mode = 
        match nodes with
        | h :: _ -> h.cost_mode
        | [] -> `Parsimony
    in
    let eg_node =
        { 
            characters = [];
            total_cost = 0.;
            node_cost = 0.;
            taxon_code = -1;
            min_child_code = -1;
            num_child_edges = 0;
            num_height = 0;
            num_otus = 1;
            exclude_sets = [];
            exclude_info = [];
            cost_mode = cost_mode;
        } 
    in
    match data.Data.complex_schema with
    | [] -> 
            nodes, data
    | schema ->
    let nodes = List.map
        (function
             | (Parser.SetGroups.Set (name, meth, list)) as t ->
             let taxon_code = add_taxon name in
             let cs, name = group taxon_code (new_sid ()) t in
             (* add taxon name *)
             { eg_node with
                   characters = [cs];
                   taxon_code = taxon_code;
                   min_child_code = taxon_code; }
             | Parser.SetGroups.Elt taxon ->
                   let taxon_id = All_sets.StringMap.find taxon
                       data.Data.taxon_names in
                   List.find (fun node -> node.taxon_code = taxon_id)
                       nodes
        )
        schema 
    in
    nodes, !data'

let load_data ?(silent=true) ?(classify=true) data = 
    (* Not only we make the list a set, we filter those characters that have
    * weight 0. *)
    current_snapshot "Node.load_data start";
    let make_set_of_list lst =
        List.fold_left (fun acc x -> 
            if 0. = Data.get_weight x data then acc
            else All_sets.Integers.add x acc) 
        All_sets.Integers.empty lst
    in
    let is_mem =
        (* We check for informative characters among all the terminals in data
        * *)
        if classify then 
            (* We need to verify if the character is
            potentially informative or not *)
            (fun char -> 
                Data.apply_boolean
                NonaddCS8.is_potentially_informative 
                AddCS.is_potentially_informative data char
                && 0. <> Data.get_weight char data)
        else (fun _ -> true)
    in
    let data, generate_taxon = 
        current_snapshot "start nonadd sets";
        let n8 = List.filter is_mem data.Data.non_additive_8
        and n16 = List.filter is_mem data.Data.non_additive_16
        and n32 = List.filter is_mem data.Data.non_additive_32
        and n33 = List.filter is_mem data.Data.non_additive_33
        and add = List.filter is_mem data.Data.additive
        and sank = List.map (List.filter is_mem) data.Data.sankoff
        and kolmogorov = List.filter is_mem data.Data.kolmogorov
        and static_ml = List.filter is_mem data.Data.static_ml
        and dynamics = List.filter is_mem data.Data.dynamics in
        current_snapshot "start nonadd set2";
        let n8 = make_set_of_list n8
        and n16 = make_set_of_list n16
        and n32 = make_set_of_list n32
        and n33 = make_set_of_list n33
        and add = make_set_of_list add in
        let cost_mode = 
            match static_ml with
            | [] -> `Parsimony
            | _ -> `Likelihood
        in
        current_snapshot "end nonadd set2";
        let r = 
            generate_taxon classify add n8 n16 n32 n33 sank dynamics kolmogorov 
            static_ml data cost_mode
        in
        current_snapshot "end generate taxon";
        r
    in
    let nodes = 
        let ntaxa = 
            All_sets.IntegerMap.fold (fun _ _ acc -> acc + 1)
            data.Data.taxon_codes 0 
        in
        let st, finalize = 
            if not silent then
                let status = 
                    Status.create "Loading terminals" (Some ntaxa)
                    "terminals loaded" 
                in
                let cnt = ref 0 in
                (fun () -> 
                    incr cnt; 
                    Status.achieved status !cnt;
                    Status.full_report status),
                (fun () -> Status.finished status)
            else (fun () -> ()), (fun () -> ())
        in
        let res = 
            All_sets.IntegerMap.fold (fun x _ acc -> 
                try 
                    let res = generate_taxon x acc in
                    st ();
                    res
                with 
                Not_found -> acc)
            data.Data.taxon_codes []
        in
        finalize ();
        res
    in
    let nodes = 
        let sorted x = List.stable_sort node_contents_compare x in
        List.map (fun x -> { x with characters = sorted x.characters }) nodes
    in
    let nodes, data = structure_into_sets data nodes in
    current_snapshot "Node.load_data end";
    data, nodes

(* OUTPUT TO XML *)
let generate_print_endline ch =
    fun str ->
        output_string ch str;
        output_string ch "\n";
        flush ch

let generate_print_string ch =
    fun str -> output_string ch str; flush ch

let generate_print_int ch =
    fun i -> output_string ch (string_of_int i); flush ch

let generate_print_float ch =
    fun f -> (output_string ch (string_of_float f)); flush ch

let print_item printf before after between =
    printf before;
    between ();
    printf after

let output_taxon_code ch i = 
    let print_endline = generate_print_endline ch in
    let between () = print_endline (string_of_int i) in
    print_item print_endline "<code>" "</code>" between

let output_total_cost ch i =
    let print_endline = generate_print_endline ch in
    let between () = print_endline (string_of_float i) in
    print_item print_endline "<cost>" "</cost>" between

let pre = [ (Tags.Characters.cclass, Tags.Nodes.preliminary) ]
let fin = [ (Tags.Characters.cclass, Tags.Nodes.final) ]
let sing = [ (Tags.Characters.cclass, Tags.Nodes.single) ]

let rec cs_to_single (pre_ref_code, fi_ref_code) (root : cs option) parent_cs mine : cs =
    match parent_cs, mine with
    | Dynamic parent, Dynamic mine ->
            (* Do we need this only for dynamic characters? I will first get it
            * going here only *)
          let root_pre = match root with
          | Some (Dynamic root) -> Some root.preliminary
          | _ -> None
          in 
          let prev_cost, cost, res = 
              DynamicCS.to_single pre_ref_code 
                    root_pre parent.preliminary mine.preliminary 
            in
          Dynamic {preliminary = res; final = res; 
                   cost = (mine.weight *.  cost);
                   sum_cost = (mine.weight *. cost);
                   weight = mine.weight; time = 0.}
              
    | _ -> mine

let to_single (pre_ref_codes, fi_ref_codes) root parent mine = 
    (*
    Printf.printf "Assigning single to %d\n%!" mine.taxon_code;
    *)
    match root with
    | Some root ->
          let root_char_opt = List.map (fun c -> Some c) root.characters in 
          { mine with 
                characters = map3 (cs_to_single (pre_ref_codes, fi_ref_codes) ) 
                  root_char_opt parent.characters mine.characters }
    | None ->
          { mine with 
                characters = map2 (cs_to_single (pre_ref_codes, fi_ref_codes) None ) 
                  parent.characters mine.characters }

let readjust mode to_adjust ch1 ch2 parent mine = 
    let ch1, ch2 =
        if ch1.min_child_code < ch2.min_child_code then ch1, ch2
        else ch2, ch1
    in
    let modified = ref All_sets.Integers.empty in
    let adjusted_likelihood = ref false in
    let cs_readjust c1 c2 parent mine =
        match c1, c2, parent, mine with 
        | StaticMl c1, StaticMl c2, StaticMl parent, StaticMl mine -> 
            IFDEF USE_LIKELIHOOD THEN
                adjusted_likelihood := true;
                let m, prev_cost, cost, (t1,t2), res = 
                    MlStaticCS.readjust to_adjust !modified c1.preliminary
                    c2.preliminary mine.preliminary c1.time c2.time mine.time
                in

                modified := m;
                let cost = mine.weight *. cost in

                (( StaticMl {mine with 
                        preliminary=res;final=res;cost=cost;sum_cost=cost;},
                StaticMl {c1 with time =t1;}), StaticMl { c2 with time = t2; })
            ELSE
                failwith likelihood_error
            END
        | ((Dynamic c1) as c1'), ((Dynamic c2) as c2'), Dynamic parent, Dynamic mine ->
                let m, prev_cost, cost, res =
                    DynamicCS.readjust mode to_adjust !modified c1.preliminary c2.preliminary
                    parent.preliminary mine.preliminary
                in
                modified := m;
                let cost = mine.weight *. cost in
                (Dynamic { preliminary = res; final = res; 
                cost = cost;
                sum_cost = c1.sum_cost +. c2.sum_cost +. cost;
                weight = mine.weight; time=0.0}, c1'), c2'
        | _ -> (mine, c1), c2
    in
    if mine.total_cost = infinity then (`Mine mine), !modified
    else
        let characters = 
            map4 cs_readjust ch1.characters ch2.characters
            parent.characters mine.characters
        in
        let characters, b, c = 
            let a, c = List.split characters in
            let a, b = List.split a in
            a, b, c 
        in
        let children_cost = ch1.total_cost +. ch2.total_cost 
        and node_cost = get_characters_cost characters in
        let total_cost =
            match mine.cost_mode with
            | `Likelihood -> node_cost 
            | `Parsimony  -> node_cost +. children_cost in
        
        let res = 
            { mine with characters = characters; total_cost = total_cost; 
            node_cost = node_cost }
        in
        if !adjusted_likelihood then
            let ch1 = { ch1 with characters = b }
            and ch2 = { ch2 with characters = c } in
            `MineNChildren (res, ch1, ch2), !modified
        else (`Mine res), !modified

let to_single_root (pre_ref_codes, fi_ref_codes) mine = 
    to_single (pre_ref_codes, fi_ref_codes) (Some mine) mine mine

(** [get_active_ref_code node_data] returns codes of
* all active chromosomes which are used for single state process
* on the subtree rooted at [node_data] *)
let get_active_ref_code node_data = 
    List.fold_left 
        (fun (acc_pre, acc_pre_child, acc_fi, acc_fi_child) cs -> 
             match cs with                   
             | Dynamic d -> 
                   let pre, pre_child =
                       DynamicCS.get_active_ref_code d.preliminary 
                   in

                   let fi, fi_child =
                       DynamicCS.get_active_ref_code d.final 
                   in

                   IntSet.union acc_pre pre, IntSet.union acc_pre_child pre_child,
                   IntSet.union acc_fi fi, IntSet.union acc_fi_child fi_child

             | _ -> acc_pre, acc_pre_child, acc_fi, acc_fi_child

        ) (IntSet.empty, IntSet.empty, IntSet.empty, IntSet.empty) node_data.characters 

let rec cs_to_formatter (pre_ref_codes, fi_ref_codes) d 
    (cs , cs_single) (parent_cs : (cs * cs) option) : Tags.output list = 
    match cs,  parent_cs, cs_single with
    | Nonadd8 cs, _, _ -> begin
          match parent_cs with 
          | None ->
                NonaddCS8.to_formatter pre cs.preliminary None d
                @ NonaddCS8.to_formatter fin cs.final None  d
          | Some ((Nonadd8 parent_cs), _) ->
                NonaddCS8.to_formatter pre cs.preliminary (Some parent_cs.preliminary) d
            @ NonaddCS8.to_formatter fin cs.final (Some parent_cs.final) d  
            | _ -> failwith "Fucking up with Nonadd at cs_to_formatter in node.ml"
      end 
    | Nonadd16 cs, _, _ -> begin 
          match parent_cs with  
          | None ->
                NonaddCS16.to_formatter pre cs.preliminary None d
                @ NonaddCS16.to_formatter fin cs.final None  d
          | Some ((Nonadd16 parent_cs), _) ->
                NonaddCS16.to_formatter pre cs.preliminary (Some parent_cs.preliminary) d
            @ NonaddCS16.to_formatter fin cs.final (Some parent_cs.final) d  
          | _ -> failwith "Fucking up with Nonadd at cs_to_formatter in node.ml" 
      end  
    | Nonadd32 cs, _, _ -> begin
          match parent_cs with  
          | None -> 
                NonaddCS32.to_formatter pre cs.preliminary None d
                @ NonaddCS32.to_formatter fin cs.final None  d 
          | Some ((Nonadd32 parent_cs), _) -> 
                NonaddCS32.to_formatter pre cs.preliminary (Some parent_cs.preliminary) d
                @ NonaddCS32.to_formatter fin cs.final (Some parent_cs.final) d  
          | _ -> failwith "Fucking up with Nonadd at cs_to_formatter in node.ml" 
      end   
    | Add cs, _, _-> begin
          match parent_cs with 
          | None ->
                AddCS.to_formatter pre cs.preliminary None d
                @ AddCS.to_formatter fin cs.final None  d
          | Some ((Add parent_cs), _) ->
                AddCS.to_formatter pre cs.preliminary (Some parent_cs.preliminary) d
                @ AddCS.to_formatter fin cs.final (Some parent_cs.final) d  
          | _ -> failwith "Fucking up with Add at cs_to_formatter in node.ml" 
      end 
    | Sank cs, _, _-> begin
          match parent_cs with 
          | None ->
                SankCS.to_formatter pre cs.preliminary None d
                @ SankCS.to_formatter fin cs.final None  d
          | Some ((Sank parent_cs), _) ->
                SankCS.to_formatter pre cs.preliminary (Some parent_cs.preliminary) d
                @ SankCS.to_formatter fin cs.final (Some parent_cs.final) d  
          | _ -> failwith "Fucking up with Sank at cs_to_formatter in node.ml" 
      end  
    | Dynamic cs, _, Dynamic cs_single -> begin
          match parent_cs with 
          | None ->
                DynamicCS.to_formatter pre_ref_codes pre cs.preliminary None d
                @ DynamicCS.to_formatter fi_ref_codes fin  cs.final None  d
                @ 
                (DynamicCS.to_formatter pre_ref_codes sing 
                cs_single.preliminary None d)

          | Some ((Dynamic parent_cs), (Dynamic parent_cs_single)) ->
                (DynamicCS.to_formatter pre_ref_codes pre cs.preliminary
                (Some parent_cs.preliminary) d) 
                @ 
                (DynamicCS.to_formatter
                fi_ref_codes fin cs.final (Some parent_cs.final)
                d) 
                @ 
                (DynamicCS.to_formatter pre_ref_codes sing 
                cs_single.preliminary (Some parent_cs_single.preliminary) d)
          | _ -> failwith "Fucking up with Dynamic at cs_to_formatter in node.ml"
      end 
    | Kolmo x, _, _ -> 
          KolmoCS.to_formatter pre_ref_codes pre x.preliminary d @
              KolmoCS.to_formatter fi_ref_codes  fin x.final d
    | Set x, _, Set x_single ->
          let attributes =
              [(Tags.Characters.name, (Data.code_character x.final.sid d))] in
          let sub a = 
              (* SET BUG!!!! I'm pasing here `Left, but I believe this is an
              * error, though I can't see where ... *)
              List.map (fun a -> `Single a) 
              (cs_to_formatter (pre_ref_codes, fi_ref_codes) d a None) in
          let sub : (Tags.output Sexpr.t list) = List.map2 (fun a b -> `Set
          (sub (a, b))) x.final.set x_single.final.set in
          let cont = `Structured (`Set sub) in
          [(Tags.Characters.set, attributes, cont)]
    | StaticMl cs,_, _ -> 
        IFDEF USE_LIKELIHOOD THEN
        begin
          match parent_cs with 
          | None ->
                MlStaticCS.to_formatter pre cs.preliminary cs.time None d
          | Some ((StaticMl parent_cs), _) ->
                MlStaticCS.to_formatter pre cs.preliminary cs.time (Some parent_cs.preliminary) d
          | _ -> failwith "Fucking up with StaticMl at cs_to_formatter in node.ml" 
        end
        ELSE
            failwith likelihood_error
        END
    | _ -> assert false

(* Compute total recost of the NODE, NOT THE SUBTREE*)
let cmp_node_recost node_data = 
    List.fold_left (fun recost cs -> 
                        match cs with
                        | Dynamic dyn -> recost +. DynamicCS.total_recost dyn.preliminary
                        | _ -> recost
                   ) 0. node_data.characters

(* Compute total recost of the subtree rooted by this NODE *)
let cmp_subtree_recost node_data = 
    List.fold_left 
        (fun subtree_recost cs -> 
             match cs with 
             | Dynamic dyn -> subtree_recost +. DynamicCS.subtree_recost dyn.preliminary
             | _ -> subtree_recost
        ) 0. node_data.characters


let to_formatter_single (pre_ref_codes, fi_ref_codes) 
    (acc : Tags.attributes) 
    d (node_data, node_single) (node_id : int) parent_data : Tags.output =
    let get_node_name id = 
        try Data.code_taxon id d 
        with | Not_found -> string_of_int id 
    in
    let node_name = get_node_name node_id in 
    let attr = 
        (Tags.Nodes.cost, string_of_float 0.) :: 
        (Tags.Nodes.recost, string_of_float 0.) :: 
        (Tags.Nodes.node_cost, string_of_float node_data.node_cost) ::
        (Tags.Nodes.name, node_name) ::
        (Tags.Nodes.child1_name, "") ::
        (Tags.Nodes.child2_name, "") ::
        (Tags.Nodes.nce, string_of_int node_data.num_child_edges) ::
        (Tags.Nodes.notu, string_of_int node_data.num_otus) :: acc
    in
    let get_parent item = 
        match parent_data with
        | None -> None
        | Some (a, b) -> 
                (Some ((List.nth a.characters item), (List.nth b.characters
                item)))
    in
    let children, _ = List.fold_left 
        (fun (acc, item) cs ->
            let parent_data = get_parent item in
             let res = 
                 cs_to_formatter 
                 (pre_ref_codes, fi_ref_codes) d cs parent_data 
            in 
             let res = List.map (fun ch -> `Single ch) res in 
             res @ acc, item + 1
        ) ([], 0) 
        (List.map2 (fun a b -> a, b) 
        node_data.characters node_single.characters)
    in
    (Tags.Nodes.node, attr, `Structured (`Set children))

(** [copy_chrom_map source des] copies the chromosome map
* which creates chromosome [source] to chromosome map which creates chromosome [des] *)
let copy_chrom_map source des =
    let s_ch_arr = Array.of_list source.characters in 
    let d_ch_arr = Array.of_list des.characters in 
    let d_ch_arr = 
        Array.mapi 
            (fun idx d_ch ->
                 match s_ch_arr.(idx), d_ch with 
                 | Dynamic s_ch, Dynamic d_ch ->
                       let pre  = DynamicCS.copy_chrom_map s_ch.preliminary
                           d_ch.preliminary in 
                       let fi = DynamicCS.copy_chrom_map s_ch.final
                           d_ch.final in

                       Dynamic {d_ch with preliminary = pre; final = fi}
                 | _, _ -> d_ch
            ) d_ch_arr
    in 
    {des with characters = Array.to_list d_ch_arr}

let to_formatter_subtree (pre_ref_codes, fi_ref_codes)
        acc d (node_data, node_single) node_id (child1_id,  child1_node_data)
        (child2_id,  child2_node_data) (parent_node_data_opt : (node_data * node_data) option) : Tags.output =

    let get_node_name id =  
        try Data.code_taxon id d with 
        | Not_found -> string_of_int id 
    in
    let child1_name = get_node_name child1_id in 
    let child2_name = get_node_name child2_id in 

    let node_name =
        match parent_node_data_opt with  
        | None -> "root"
        | _ -> get_node_name node_id
    in 
    let child1_recost = cmp_subtree_recost child1_node_data in 
    let child2_recost = cmp_subtree_recost child2_node_data in 
    let subtree_recost = 
        child1_recost +. 
        child2_recost +. 
        (cmp_node_recost node_data) in 
    let attr = 
        (Tags.Nodes.cost, 
         (match parent_node_data_opt with 
          | None -> "0." 
          | _ -> (string_of_float node_data.total_cost))):: 
        (Tags.Nodes.recost, string_of_float subtree_recost) :: 
        (Tags.Nodes.node_cost, string_of_float node_data.node_cost) ::
        (Tags.Nodes.name, node_name) ::
        (Tags.Nodes.child1_name, child1_name) ::        
        (Tags.Nodes.child2_name, child2_name) ::        
        (Tags.Nodes.nce, string_of_int node_data.num_child_edges) ::
        (Tags.Nodes.notu, string_of_int node_data.num_otus) :: acc
    in
    let children, _ = List.fold_left 
        (fun (acc, idx) cs ->
             let parent_cs = 
                 match parent_node_data_opt with 
                 | Some (parent_node_data, single_parent_node_data) -> 
                         let pn = List.nth parent_node_data.characters idx 
                         and pnsingle = List.nth
                             single_parent_node_data.characters idx
                         in 
                         Some (pn, pnsingle)
                 | _ -> None
             in  
             let res = cs_to_formatter (pre_ref_codes, fi_ref_codes) d cs parent_cs in 

             let res = List.map (fun ch -> `Single ch) res in 
             (res @ acc), (idx + 1)) ([], 0) (List.map2 (fun a b -> a, b) node_data.characters 
             node_single.characters)
    in
    (Tags.Nodes.node, attr, `Structured (`Set children))

let to_xml data ch node = ()

type 'a converter = 'a r -> cs list -> cs list

let listify a acc = a :: acc
let lnon8 cs acc = listify (Nonadd8 cs) acc
let lnon16 cs acc = listify (Nonadd16 cs) acc
let lnon32 cs acc = listify (Nonadd32 cs) acc
let ladd cs acc = listify (Add cs) acc
let lsank cs acc = listify (Sank cs) acc
let ldynamic cs acc = listify (Dynamic cs) acc
let lkolmo cs acc = listify (Kolmo cs) acc
let lstaticml cs acc = listify (StaticMl cs) acc

let rec convert_data
        ?(tnon8=lnon8)
        ?(tnon16=lnon16)
        ?(tnon32=lnon32)
        ?(tadd=ladd)
        ?(tsank=lsank)
        ?(tdynamic=ldynamic)
        ?(tkolmo=lkolmo)
        ?(tstaticml=lstaticml)
        node =
    let rec conv cs acc = match cs with
             | Nonadd8 cs -> tnon8 cs acc
             | Nonadd16 cs -> tnon16 cs acc
             | Nonadd32 cs -> tnon32 cs acc
             | Add cs -> tadd cs acc
             | Sank cs -> tsank cs acc
             | Dynamic cs -> tdynamic cs acc
             | Kolmo cs -> tkolmo cs acc
             | StaticMl cs -> tstaticml cs acc
             | Set set ->
                   let res = List.fold_right conv set.preliminary.set [] in
                   let res = { set.preliminary with set = res } in
                   (Set { set with
                              preliminary = res;
                              final = res; }) :: acc 
    in
    let characters' =
        List.fold_right
            conv
            node.characters [] in
    { node with characters = characters' }

let of_data data =
    load_data data

let filter_characters f n =
    { n with characters = f n.characters }

let do_filter cardinal f c codes =
    if 0 = cardinal c.preliminary then c
    else
        { c with preliminary = f c.preliminary codes; final = f c.final codes }

let rec filter_character_codes (codes : All_sets.Integers.t) item = 
    match item with
    | StaticMl c ->
        IFDEF USE_LIKELIHOOD THEN
            StaticMl (do_filter MlStaticCS.cardinal MlStaticCS.f_codes c codes)
        ELSE
            failwith likelihood_error
        END
    | Nonadd8 c ->
          Nonadd8 (do_filter NonaddCS8.cardinal NonaddCS8.f_codes c codes)
    | Nonadd16 c ->
          Nonadd16 (do_filter NonaddCS16.cardinal NonaddCS16.f_codes c codes)
    | Nonadd32 c ->
          Nonadd32 (do_filter NonaddCS32.cardinal NonaddCS32.f_codes c codes)
    | Add c ->
          Add (do_filter AddCS.cardinal AddCS.f_codes c codes)
    | Sank c -> 
          Sank (do_filter SankCS.cardinal SankCS.f_codes c codes)
    | Dynamic c ->
          Dynamic (do_filter DynamicCS.cardinal DynamicCS.f_codes c codes)
    | Kolmo c ->
            Kolmo (do_filter KolmoCS.cardinal KolmoCS.f_codes c codes)
    | Set s ->
          if All_sets.Integers.mem s.final.sid codes
          then Set s
          else begin
                  let res = List.map (filter_character_codes codes)
                      s.preliminary.set in
                  let res = { s.preliminary with set = res } in
                  Set { s with preliminary = res; final = res; }
              end

let rec filter_character_codes_complement codes = function
    | StaticMl c ->
        IFDEF USE_LIKELIHOOD THEN
          Some  (StaticMl (do_filter MlStaticCS.cardinal MlStaticCS.f_codes_comp c codes))
        ELSE
          failwith likelihood_error
        END
    | Nonadd8 c ->
          Some (Nonadd8 (do_filter NonaddCS8.cardinal NonaddCS8.f_codes_comp c codes))
    | Nonadd16 c ->
          Some (Nonadd16 (do_filter NonaddCS16.cardinal NonaddCS16.f_codes_comp c codes))
    | Nonadd32 c ->
          Some (Nonadd32 (do_filter NonaddCS32.cardinal NonaddCS32.f_codes_comp c codes))
    | Add c ->
          Some (Add (do_filter AddCS.cardinal AddCS.f_codes_comp c codes))
    | Sank c -> 
          Some (Sank (do_filter SankCS.cardinal SankCS.f_codes_comp c codes))
    | Dynamic c ->
          Some (Dynamic (do_filter DynamicCS.cardinal DynamicCS.f_codes_comp c codes))
    | Kolmo c ->
          Some (Kolmo (do_filter KolmoCS.cardinal KolmoCS.f_codes_comp c codes))
    | Set s ->
          if All_sets.Integers.mem s.final.sid codes
          then None
          else begin
                  let res = List.map (filter_character_codes_complement codes)
                      s.preliminary.set in
                  let res = List.fold_left
                      (fun list elt ->
                           match elt with
                           | None -> list
                           | Some e -> e :: list)
                      [] res in
                  let res = List.rev res in
                  let res = { s.preliminary with set = res } in
                  Some (Set { s with preliminary = res; final = res; })
              end

let rec mmap fn list =
    match list with
    | [] -> []
    | i :: is ->
          match fn i with
          | Some v -> v :: mmap fn is
          | None -> mmap fn is

let f_codes codes n =
    let codes = 
        List.fold_left 
        (fun acc x -> All_sets.Integers.add x acc)
        All_sets.Integers.empty
        codes
    in
    let chars = List.map (filter_character_codes codes) n.characters in
    { n with characters = chars }

let f_codes_comp codes n =
    let codes = 
        List.fold_left 
        (fun acc x -> All_sets.Integers.add x acc)
        All_sets.Integers.empty
        codes
    in
    { n with characters =
            mmap (filter_character_codes_complement codes) n.characters }


let f_characters n = 
    (* We are required, for this operation, to have only one kind of character
    * filtered out in the node n *)
    match n.characters with
    | [item] -> item
    | _ -> 
            failwith "To extract a list of characters, a single type of \
            character is acceptable in the set"

let add_exclude set n =
    let excluded = All_sets.Integers.mem n.taxon_code set in
    let card = All_sets.Integers.cardinal set in
    { n with
          exclude_sets =
            set :: n.exclude_sets;
          exclude_info =
            (if excluded
             then (`Excluded, 0, card, 1)
             else (`NotExcluded, 0, card, 0)) :: n.exclude_info;
    }


let print_node_data (data : node_data) = 
    Printf.fprintf stdout "Taxon code: %i, number characters: %i\n" data.taxon_code (List.length data.characters);
    let dynamic_ls = List.filter (fun ch -> 
                                  match ch with 
                                  | Dynamic _ -> true
                                  | _ -> false ) data.characters 
    in 
    Printf.fprintf stdout "Number dynamics: %i, " (List.length dynamic_ls); 

    print_newline ()


let get_sequences _ node =
    let all_seq = 
        List.fold_left 
        (fun acc x -> 
            match x with Dynamic x -> 
                (match x.preliminary with
                | DynamicCS.SeqCS x -> x :: acc 
                | _ -> acc)
            | _ -> acc)
        []
        node.characters
    in
    let all_seq = List.map SeqCS.explode all_seq in
    List.flatten all_seq


(*============================================================*)
(*============================================================*)

let fprintf = Printf.fprintf 

let prioritize node = node

let reprioritize a b = b

let is_collapsable clas a b = 
    0.0 = edge_distance clas a b

let rec internal_n_chars acc (chars : cs list) = 
    match chars with
    | [] -> acc
    | c :: cs -> 
            let count = 
                match c with
                | StaticMl r -> 
                    IFDEF USE_LIKELIHOOD THEN
                        MlStaticCS.cardinal r.preliminary
                    ELSE
                        failwith likelihood_error
                    END
                | Nonadd8 r -> NonaddCS8.cardinal r.preliminary
                | Nonadd16 r -> NonaddCS16.cardinal r.preliminary
                | Nonadd32 r -> NonaddCS32.cardinal r.preliminary
                | Add r -> AddCS.cardinal r.preliminary
                | Sank r -> SankCS.cardinal r.preliminary
                | Dynamic r -> DynamicCS.cardinal r.preliminary
                | Kolmo r -> KolmoCS.cardinal r.preliminary
                | Set s -> internal_n_chars acc s.preliminary.set
            in
            internal_n_chars (acc + count) cs
(** [n_chars chars] returns the total number of characters in the list of
    characters [chars] *)
let rec n_chars ?(acc=0) n =
    internal_n_chars acc n.characters

module T = struct
    let add_exclude = add_exclude
end

module Union = struct
    let debug = false

    type 'a ru = {
        ch : 'a;
        u_weight : float;
    }

IFDEF USE_LIKELIHOOD THEN
    type ml_repu = MlStaticCS.t ru
ELSE
    type ml_repu = unit
END

    type c = 
        | Nonadd8U of NonaddCS8.u ru
        | Nonadd16U of NonaddCS16.u ru
        | Nonadd32U of NonaddCS32.u ru
        | AddU of AddCS.t ru
        | SankU of SankCS.t ru
        | DynamicU of DynamicCS.u ru
        | KolmoU of KolmoCS.u ru
        | StaticMlU of ml_repu

    (* In this first implementation, sets can not use the union heuristics, this is
    * work to be done *)

    type u = {
        charactersu : c list;
        min_child_codeu : int;
    }

    type p = Preliminary | Final

    let get_sequence _ code c = 
        let res = 
            List.fold_left 
            (fun acc c ->
                match acc with
                | Some _ -> acc
                | None -> 
                        match c with
                        | Nonadd8U _
                        | Nonadd16U _
                        | Nonadd32U _
                        | AddU _ 
                        | StaticMlU _
                        | SankU _ -> None
                        | KolmoU x 
                        | DynamicU x -> DynamicCS.get_sequence_union code x.ch)
            None
            c.charactersu
        in
        match res with
        | Some v -> v
        | None -> 
                failwith 
                ("Node.get_sequence: could not find code") 
                (* Warning: This exception is catched in charTransform for proper 
                * behavior of the program *)

    let saturation x =
        let single_saturation (polyacc, lenacc) = function
            | Nonadd8U x -> 
                    let card = float_of_int (NonaddCS8.cardinal_union x.ch) in
                    (NonaddCS8.poly_saturation x.ch 1) *. card, card +. lenacc
            | Nonadd16U x -> 
                    let card = float_of_int (NonaddCS16.cardinal_union x.ch) in
                    (NonaddCS16.poly_saturation x.ch 1)*. card, card +. lenacc
            | Nonadd32U x -> 
                    let card = float_of_int (NonaddCS32.cardinal_union x.ch) in
                    (NonaddCS32.poly_saturation x.ch 1) *. card, card +. lenacc
            | StaticMlU x ->
                    failwith "TODO Node.Union.poly_saturation"
            | AddU x -> 
                    failwith "TODO Node.Union.poly_saturation"
                    (*
                    let card = float_of_int (AddCS.cardinal x) in
                    AddCS.poly_saturation x *. card, card +. lenacc
                    *)
            | SankU x -> 
                    failwith "TODO Node.Union.poly_saturation"
                    (*
                    let card = float_of_int (SankCS.cardinal x) in
                    SankCS.poly_saturation x *. card, card +. lenacc
                    *)
            | KolmoU x
            | DynamicU x ->
                    let card = float_of_int (DynamicCS.cardinal_union x.ch) in
                    (DynamicCS.poly_saturation x.ch 1) *. card, card +. lenacc
        in
        let a, b = 
            List.fold_left single_saturation (0.0, 0.0) x.charactersu 
        in
        a /. b

    let union _ cs u1 u2 =
        let u1, u2 = 
            if u1.min_child_codeu < u2.min_child_codeu then u1, u2
            else u2, u1
        in
        if debug then 
            Printf.printf "I am calculating the union of %d and %d.\n%!"
            u1.min_child_codeu u2.min_child_codeu;
        let single_union cs u1 u2 =
            match cs, u1, u2 with
            | Nonadd8 cs, Nonadd8U u1, Nonadd8U u2 ->
                    let r = {
                        ch = (NonaddCS8.union cs.preliminary u1.ch u2.ch);
                        u_weight = u1.u_weight;
                    }
                    in
                    Nonadd8U r
            | Nonadd16 cs, Nonadd16U u1, Nonadd16U u2 ->
                    let r = {
                        ch = (NonaddCS16.union cs.preliminary u1.ch u2.ch);
                        u_weight = u1.u_weight;
                    }
                    in
                    Nonadd16U r
            | Nonadd32 cs, Nonadd32U u1, Nonadd32U u2 ->
                    let r = {
                        ch = (NonaddCS32.union cs.preliminary u1.ch u2.ch);
                        u_weight = u1.u_weight;
                    }
                    in
                    Nonadd32U r
            | Add cs, AddU u1, AddU u2 -> 
                    AddU { ch = cs.preliminary; u_weight = u1.u_weight }
            | Sank cs, SankU u1, SankU u2 ->
                    SankU { ch = cs.preliminary; u_weight = u1.u_weight }
            | Dynamic cs, DynamicU u1, DynamicU u2 ->
                    let r = {
                        ch = (DynamicCS.union cs.preliminary u1.ch u2.ch);
                        u_weight = u1.u_weight;
                    }
                    in
                    DynamicU r
            | Kolmo cs, KolmoU u1, KolmoU u2 ->
                    let r = {
                        ch = (KolmoCS.union cs.preliminary u1.ch u2.ch);
                        u_weight = u1.u_weight;
                    }
                    in
                    KolmoU r
            | StaticMl cs, StaticMlU u1, StaticMlU u2 ->
                IFDEF USE_LIKELIHOOD THEN
                    StaticMlU { ch = cs.preliminary; u_weight = u1.u_weight; }
                ELSE
                    failwith likelihood_error
                END
    
            | _ -> failwith "Node.Union.union"
        in
        let rec map3 c d e =
            match c, d, e with
            | hc :: tc, hd :: td, he :: te ->
                    (single_union hc hd he) :: (map3 tc td te)
            | [], [], [] -> []
            | _ -> failwith "Node.Union.map4"
        in
        let c = map3 cs.characters u1.charactersu u2.charactersu in
        { charactersu = c; 
        min_child_codeu = min u1.min_child_codeu u2.min_child_codeu }

    let union_preliminary _ u par =
        let single_union u par =
            match u, par with
            | Nonadd8U u, Nonadd8 par ->
                    Nonadd8U { ch = NonaddCS8.union_states u.ch par.preliminary; 
                    u_weight = u.u_weight; }
            | Nonadd16U u, Nonadd16 par ->
                    Nonadd16U { ch = NonaddCS16.union_states u.ch
                    par.preliminary;
                    u_weight = u.u_weight; }
            | Nonadd32U u, Nonadd32 par ->
                    Nonadd32U { ch = NonaddCS32.union_states u.ch par.preliminary;
                    u_weight = u.u_weight; }
            | x, _ -> x
        in
        let rec map2 c d = 
            match c, d with
            | hc :: tc, hd :: td -> (single_union hc hd) :: map2 tc td
            | [], [] -> []
            | _, _ -> failwith "Node.Union.union_preliminary"
        in
        { u with charactersu = map2 u.charactersu par.characters }

    let union_final _ u par =
        let single_union u par =
            match u, par with
            | Nonadd8U u, Nonadd8 par ->
                    Nonadd8U { ch = NonaddCS8.union_states u.ch par.final;
                    u_weight = u.u_weight; }
            | Nonadd16U u, Nonadd16 par ->
                    Nonadd16U { ch = NonaddCS16.union_states u.ch par.final;
                    u_weight = u.u_weight }
            | Nonadd32U u, Nonadd32 par ->
                    Nonadd32U { ch = NonaddCS32.union_states u.ch par.final;
                    u_weight = u.u_weight }
            | x, _ -> x
        in
        let rec map2 c d = 
            match c, d with
            | hc :: tc, hd :: td -> (single_union hc hd) :: map2 tc td
            | [], [] -> []
            | _, _ -> failwith "Node.Union.union_final"
        in
        { u with charactersu = map2 u.charactersu par.characters }

    let leaf taxon_code _ c =
        let create_union f x = 
            { ch = f x.preliminary; u_weight = x.weight }
        in
        let single_leaf acc char = 
            match char with
            | Nonadd8 c ->
                    (Nonadd8U (create_union NonaddCS8.to_union c)) :: acc
            | Nonadd16 c ->
                    (Nonadd16U (create_union NonaddCS16.to_union c)) :: acc
            | Nonadd32 c ->
                    (Nonadd32U (create_union NonaddCS32.to_union c)) :: acc
            | Dynamic c ->
                    (DynamicU (create_union DynamicCS.to_union c)) :: acc 
            | Add c ->
                    (AddU { ch = c.preliminary; u_weight = c.weight }) :: acc
                    (* TODO
                    Add (AddCS.to_union c.Node.preliminary)
                    *)
            | Sank c ->
                    (SankU { ch = c.preliminary; u_weight = c.weight}) :: acc
                    (* TODO
                    Sank (SankCS.to_union c.Node.preliminary)
                    *)
            | Kolmo c ->
                    (KolmoU (create_union KolmoCS.to_union c)) :: acc
            | StaticMl c ->
                IFDEF USE_LIKELIHOOD THEN
                    (StaticMlU { ch = c.preliminary; u_weight = c.weight }) :: acc
                ELSE
                    failwith likelihood_error
                END
            | Set _ -> failwith "Node.Union.leaf TODO"
        in
        let nc = List.fold_left single_leaf [] c.characters in
        let nc = List.rev nc in
        let code = 
            match taxon_code with
            | None -> c.taxon_code
            | Some v -> v
        in
        { charactersu = nc; min_child_codeu = code }

    let distance a b = 
        let rec distance acc a b = 
            match a, b with
            | (Nonadd8U a) :: at, (Nonadd8U b) :: bt ->
                    distance (acc +. (a.u_weight *. (NonaddCS8.distance_union
                    a.ch b.ch))) at bt
            | (Nonadd16U a) :: at, (Nonadd16U b) :: bt ->
                    distance (acc +. (a.u_weight *. (NonaddCS16.distance_union
                    a.ch b.ch))) at bt
            | (Nonadd32U a) :: at, (Nonadd32U b) :: bt ->
                    distance (acc +. (a.u_weight *. (NonaddCS32.distance_union
                    a.ch b.ch))) at bt
            | (StaticMlU a) :: at, (StaticMlU b) :: bt ->
                    distance acc at bt
                    (* TODO
                    MlStaticCS.distance_union a b, at, bt
                    *)
            | (AddU a) :: at, (AddU b) :: bt ->
                    distance acc at bt
                    (* TODO
                    AddCS.distance a b, at, bt
                    *)
            | (SankU a) :: at, (SankU b) :: bt ->
                    distance acc at bt
                    (* TODO
                    SankCS.distance a b, at, ct
                    *)
            | (KolmoU a) :: at, (KolmoU b) :: bt ->
                    distance (acc +. (a.u_weight *. (KolmoCS.distance_union a.ch
                    b.ch))) at bt
            | (DynamicU a) :: at, (DynamicU b) :: bt ->
                    distance (acc +. (a.u_weight *. (DynamicCS.distance_union
                    a.ch b.ch))) at bt
            | [], [] -> acc
            | (Nonadd8U _) :: _, _
            | (Nonadd16U _) :: _, _
            | (Nonadd32U _) :: _, _
            | (AddU _) :: _, _
            | (SankU _) :: _, _
            | (DynamicU _) :: _, _
            | (KolmoU _) :: _, _ 
            | (StaticMlU _) :: _, _ 
            | [], _ -> failwith "Node.Union.distance TODO"
        in
        distance 0.0 a.charactersu b.charactersu

    let distance_node code a b =
        let a = leaf None code a in
        distance a b

    let single_compare a b =
        match a, b with
        | (Nonadd8U a), (Nonadd8U b) ->
                NonaddCS8.compare_union a.ch b.ch
        | (Nonadd16U a), (Nonadd16U b) ->
                NonaddCS16.compare_union a.ch b.ch
        | (Nonadd32U a), (Nonadd32U b) ->
                NonaddCS32.compare_union a.ch b.ch
        | (StaticMlU a), (StaticMlU b) ->
            IFDEF USE_LIKELIHOOD THEN
                MlStaticCS.compare_data a.ch b.ch
            ELSE
                failwith likelihood_error
            END
        | (AddU a), (AddU b) ->
                AddCS.compare_data a.ch b.ch
                (* TODO
                AddCS.distance a b, at, bt
                        *)
        | (SankU a), (SankU b) ->
                SankCS.compare_data a.ch b.ch
                (* TODO
                0.0, at, bt
                SankCS.distance a b, at, ct
                *)
        | (KolmoU a), (KolmoU b)
        | (DynamicU a), (DynamicU b) ->
                DynamicCS.compare_union a.ch b.ch
        | _ -> failwith "Node.Union.distance TODO"

    let rec compare a b =
        match a, b with
        | [], [] -> 0
        | [], _ -> -1
        | _, [] -> 1
        | ha :: ta, hb :: tb ->
                match single_compare ha hb with
                | 0 -> compare ta tb
                | x -> x

    let compare a b = compare a.charactersu b.charactersu
end 
(* Given a map of characters, their codes and a node, create a new node that
* removes the set of characters and adds the characters stored in the map
* for it *)
let build_node static_characters chars node = 
    let new_char = 
        All_sets.IntegerMap.find node.taxon_code static_characters 
    in
    let new_node = f_codes_comp chars node in
    { new_node with characters = (Nonadd8 new_char) ::
        new_node.characters }

let support_chars starting _ n = 
            let filter_map fn lst =
                let lst = List.map fn lst in
                let lst = List.filter (function Some _ -> true | _ -> false) lst in
                let lst = List.map (function Some a -> a | _ -> assert false) lst in
                lst in
            let char_list =
                let extract n =
                    match n with
                    | Nonadd8 r ->
                          let code = List.hd (NonaddCS8.codes r.final) in
                          if code >= starting
                          then Some r
                          else None
                    | _ -> None
                in
                filter_map extract n.characters
            in
            List.map
                (fun char ->
                     (* clade not present iff cost > 1 *)
                     char.sum_cost,
                     List.hd (NonaddCS8.codes char.final) - starting)
                char_list

let get_dynamic_preliminary data =
    let data = data.characters in
    let data =
        List.filter (function Dynamic _ | Kolmo _ -> true | _ -> false) data
    in

    List.map 
        (function 
            | Dynamic d -> d.preliminary 
            | Kolmo d -> KolmoCS.get_dynamic_preliminary d.preliminary
            | _ -> failwith "Impossible")
        data

let new_characters character_code acc taxa =
    let max_code = ref 0 in
    (* A function that takes a sequence code and a sequence and 
    * it to the current set of elements for a character creation *)
    let add_sequence_to_character_set _ sequence (cur_code, acc) =
        let add_base (counter, acc) item =
            max_code := counter;
            (counter + 1, ((counter, (counter, item), 0.0) :: acc)) 
        in
        Array.fold_left add_base (cur_code, acc) sequence
    in
    (* A function that will handle a taxon and a list of sets of sequences *)
    let build_taxon current_code acc (taxcode, sequences) =
        let handle acc sequence = 
            All_sets.IntegerMap.fold add_sequence_to_character_set sequence acc
        in
        let _, lst = List.fold_left handle (current_code, []) sequences in
        let chars = NonaddCS8.of_list lst in
        let chars = 
            { preliminary = chars; final = chars; cost = 0.0;
              sum_cost = 0.0;
            weight = 1.0 ; time = 0.}
        in
        All_sets.IntegerMap.add taxcode chars acc
    in
    (* Add each taxon to a map where each taxon is assigned a nonadditive
    * character set *)
    List.fold_left (build_taxon character_code) acc taxa

let for_support starting leaves leaves_id nodes =
    let cost_mode = 
        match leaves with
        | (_, h) :: t ->
                let mode = h.cost_mode in
                assert (List.for_all (fun (_, x) -> x.cost_mode = mode) t);
                mode
        | [] -> `Parsimony
    in
    let leaves = List.map (fun (x, y) -> (x, [], y)) leaves in
    let load_clade node_data node_id =
        List.map
            (fun (leaf_id, charlist, leaf_data) ->
                 let char = NonaddCS8.of_list [(node_id + starting,
                                                (node_id + starting,
                                                 if List.mem leaf_id
                                                     leaves_id
                                                 then 2
                                                 else 1), 0.)] in
                 (leaf_id, char :: charlist, leaf_data))
            node_data in
    let node_data = List.fold_left load_clade leaves nodes in
    let node_data = List.map
        (fun (id, charlist, leaf_data) ->
             let cslist =
                 List.fold_left
                     (fun l c ->
                          Nonadd8
                              { preliminary = c;
                                final = c;
                                cost = 0.;
                                sum_cost = 0.;
                                weight = 0.;
                                time = 0.;
                              } :: l)
                     leaf_data.characters charlist 
            in
             { characters = cslist;
               total_cost = 0.;
               node_cost = 0.;
               taxon_code = id;
               min_child_code = id;
               num_child_edges = 0;
               num_height = 0;
               num_otus = 1;
               exclude_sets = [];
               exclude_info = [];
               cost_mode = cost_mode;
             })
        node_data 
    in
    node_data

let compare_uppass = compare_data_final
let compare_downpass = compare_data_preliminary

let set_node_cost a b = { b with node_cost = a }

let extract_time nd = 
    IFDEF USE_LIKELIHOOD THEN
        let rec t_filter chs = match chs with
            | [] -> []
            | StaticMl a :: ttail -> a.time :: t_filter ttail
            | _ :: ttail -> t_filter ttail in
        let lst = nd.characters in t_filter lst 
    ELSE
        []
    END

module Standard : 
    NodeSig.S with type e = exclude and type n = node_data and type other_n =
        node_data =
        struct
        type e = exclude
        type n = node_data
        type other_n = n
        let to_other x = x
        type nad8 = NonaddCS8.t r
        let recode f n = recode f n
        let fix_preliminary = all_prelim_to_final
        let distance = distance
        let set_exclude_info a b = { b with exclude_info = a }
        let excludes_median _ = excludes_median
        let character_costs _ = character_costs
        let get_characters _ = get_characters_of_type
        let median _ = median
        let median_3 _ = median_3
        let to_string = to_string
        let total_cost = total_cost
        let node_cost _ a = a.node_cost
        let update_leaf = update_leaf
        let exclude_info _ x = x.exclude_info
        let has_excluded = has_excluded
        let taxon_code x = x.taxon_code
        (* TODO This function must be removed *)
        let union_distance _ _ = 0.
        let is_collapsable = is_collapsable
        let to_xml = to_xml
        let median_self_cost = median_self_cost
        let num_height _ x = x.num_height
        let num_otus _ x = x.num_otus
        let get_sequences = get_sequences
        let build_node = build_node
        let get_dynamic_preliminary _ = get_dynamic_preliminary
        let get_dynamic_adjusted _ _ = 
            failwith "No single assignment available"
        let edge_distance = edge_distance `Any
        let support_chars = support_chars
        let load_data = load_data
        let n_chars = n_chars
        let f_codes = f_codes
        let new_characters = new_characters
        let min_child_code _ x = x.min_child_code 
        let prioritize = prioritize
        let reprioritize = reprioritize
        let extract_time _  = extract_time
        let apply_time = apply_time
        let edge_iterator = edge_iterator
        module T = T
        module Union = Union
        let for_support = for_support
        let root_cost = root_cost
        let to_single root  _ a _ b = to_single (IntSet.empty, IntSet.empty) root a b
        let get_nonadd_8 _ = get_nonadd_8 
        let get_nonadd_16 _ = get_nonadd_16 
        let get_nonadd_32 _ = get_nonadd_32 
        let get_add _ = get_add
        let get_sank _ = get_sank
        let get_dynamic _ = get_dynamic
        let get_mlstatic _ = get_mlstatic
        let compare a b = 
            let rec aux_cmt a b =
                match a, b with
                | (Nonadd8 ha) :: ta, (Nonadd8 hb) :: tb ->
                        let cmp = compare ha.preliminary hb.preliminary in
                        if cmp = 0 then 
                            aux_cmt ta tb
                        else cmp
                | (Nonadd8 ha) :: ta, _ -> -1
                | (Nonadd16 ha) :: ta, (Nonadd16 hb) :: tb ->
                        let cmp = compare ha.preliminary hb.preliminary in
                        if cmp = 0 then 
                            aux_cmt ta tb
                        else cmp
                | (Nonadd16 ha) :: ta, _ -> -1
                | (Nonadd32 ha) :: ta, (Nonadd32 hb) :: tb ->
                        let cmp = compare ha.preliminary hb.preliminary in
                        if cmp = 0 then 
                            aux_cmt ta tb
                        else cmp
                | (Nonadd32 ha) :: ta, _ -> -1
                | (Add ha) :: ta, (Add hb) :: tb ->
                        let cmp = compare ha.preliminary hb.preliminary in
                        if cmp = 0 then 
                            aux_cmt ta tb
                        else cmp
                | (Add ha) :: ta, _ -> -1
                | (Sank ha) :: ta, (Sank hb) :: tb ->
                        let cmp = compare ha.preliminary hb.preliminary in
                        if cmp = 0 then 
                            aux_cmt ta tb
                        else cmp
                | (Sank ha) :: ta, _ -> -1
                | (Dynamic ha) :: ta, (Dynamic hb) :: tb ->
                        let cmp = compare ha.preliminary hb.preliminary in
                        if cmp = 0 then 
                            aux_cmt ta tb
                        else cmp
                | (Dynamic ha) :: ta, _ -> -1
                | (Kolmo ha) :: ta, (Kolmo hb) :: tb ->
                        let cmp = compare ha.preliminary hb.preliminary in
                        if cmp = 0 then 
                            aux_cmt ta tb
                        else cmp
                | (Kolmo ha) :: ta, _ -> -1
                | (StaticMl ha) :: ta, (StaticMl hb) :: tb ->
                        IFDEF USE_LIKELIHOOD THEN
                        let cmp = compare ha.preliminary hb.preliminary in
                        if cmp = 0 then 
                            aux_cmt ta tb
                        else cmp
                        ELSE 
                        0
                        END
                | (StaticMl ha) :: ta, _ -> -1
                | (Set ha) :: ta, (Set hb) :: tb ->
                        let cmp = compare ha.preliminary hb.preliminary in
                        if cmp = 0 then 
                            aux_cmt ta tb
                        else cmp
                | (Set ha) :: ta, _ -> 1
                | [], [] -> 0
                | [], _ -> -1
            in
            aux_cmt a.characters b.characters

        let force x = x
end 

let merge a b =
    { a with
        characters = a.characters @ b.characters;
        total_cost = a.total_cost +. b.total_cost;
        node_cost = a.node_cost +. b.node_cost;
    }

let total_cost_of_type t n =
    let rec total_cost_cs acc item = 
        match item, t with
        | Nonadd8 x, `Nonadd -> acc +. (x.sum_cost *. x.weight)
        | Nonadd16 x, `Nonadd -> acc +. (x.sum_cost *. x.weight)
        | Nonadd32 x, `Nonadd -> acc +. (x.sum_cost *. x.weight)
        | Add x, `Add -> acc +. (x.sum_cost *. x.weight)
        | Sank x, `Sank -> acc +. (x.sum_cost *. x.weight)
        | Set x, t -> List.fold_left total_cost_cs acc x.preliminary.set
        | StaticMl x, `StaticMl ->
           IFDEF USE_LIKELIHOOD THEN
            acc +. (x.cost)
           ELSE
            failwith likelihood_error
           END
        | Dynamic x, t -> 
                (match x.preliminary, t with
                | DynamicCS.SeqCS _, `Seq -> acc +. (x.sum_cost *. x.weight)
                | DynamicCS.BreakinvCS _, `Breakinv -> 
                        acc +. (x.sum_cost *.  x.weight)
                | DynamicCS.ChromCS _, `Chrom -> 
                        acc +. (x.sum_cost *. x.weight)
                | DynamicCS.AnnchromCS _, `Annchrom -> 
                        acc +. (x.sum_cost *. x.weight)
                | DynamicCS.GenomeCS _, `Genome -> 
                        acc +. (x.sum_cost *. x.weight)
                | _ -> acc)
        | _ -> acc
    in
    List.fold_left total_cost_cs 0.0 n.characters
