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

let () = SadmanOutput.register "Node" "$Revision: 2773 $"
let infinity = float_of_int max_int

open Numerical.FPInfix

let debug           = false
let debug_bl        = false
let debug_exclude   = false
let debug_sets      = false
let debug_set_cost  = false
let debug_treebuild = false
let debug_tosingle  = false
let debug_formatter = false
let debug_distance  = false
let debug_cs_median = false

let likelihood_error = 
    "Likelihood not enabled: download different binary or contact mailing list" 

let (-->) b a = a b

let odebug = Status.user_message Status.Information
let info_user_message format = 
    Printf.ksprintf (Status.user_message Status.Information) format
let failwithf format = Printf.ksprintf (failwith) format

let pp_fopt chan fopt = match fopt with
    | Some f -> output_string chan (string_of_float f)
    | None   -> output_string chan "none"

module IntSet = All_sets.Integers

exception Illegal_argument of string

type to_single = 
    [ `Add | `Annchrom | `Breakinv | `Chrom | `Genome | `Kolmo | `Nonadd
    | `FixedStates | `Sank | `Seq | `StaticMl | `Ml ]

type 'a r = {
    preliminary : 'a;
    final : 'a;
    cost : float; (*cost of this charactor on this node,* by weight*)
    sum_cost : float; (*cost of subtree root on this node, * by weight*)
    weight : float;
    time : float option * float option * float option;
}

(** Schemes for origin and loss costs *)
type origin_loss_cost = [ | `Flat of int * int ]
(** Flat cost for each origin and loss *)

(** Methods for complex terminal alignments.  These variants also contain the
    extra median information appropriate for that set type. *)
type complex_term_method = [
    | `Strictly_Same
    (** Don't allow recombination between elements of a set:  sets are only used for grouping *)
    | `Any_Of of ((int * int * (int * int) list) * float)
    (** Allow choosing of any of the elements to use for the median.  The median
        information is:
        - ID of "left" node from which this median was made
        - ID of "right" node from which this median was made
        - List with the same length as the set;  the contents are the indices of the
          element from which that median was made. *)
]

let same_ct_method a b = match a, b with
    | `Strictly_Same, `Strictly_Same -> true
    | `Any_Of _, `Any_Of _ -> true
    | `Strictly_Same, _ 
    | `Any_Of _, _ -> false

let ct_same_methods = Utl.pairwisep same_ct_method

type 'a css = {
    sid : int;
    set : 'a list;
    smethod : complex_term_method;
}
let empty_css = { 
    sid = -1;
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
    | AddVec of AddCS.Vector.t r
    | AddGen of AddCS.General.t r
    | Sank of SankCS.t r
    | FixedStates of Fixed_states.t_w_seqtbl r
    | Dynamic of DynamicCS.t r
    | Kolmo of KolmoCS.t r
    | Set of cs css r
    | StaticMl of ml_rep

let cs_string = function 
    | Nonadd8 _ -> "nonadd8"
    | Nonadd16 _ -> "nonadd16"
    | Nonadd32 _ -> "nonadd32"
    | AddVec _ -> "add"
    | AddGen _ -> "add"
    | Sank _ -> "sank"
    | FixedStates _ -> "fixed_states"
    | Dynamic x -> DynamicCS.name_string x.preliminary
    | Set _ -> "set"
    | Kolmo _ -> "kolmo"
    | StaticMl _ -> "staticml"

let rec to_string_ch ch1 = match ch1 with
    | Nonadd8 a ->
            ("na8: " ^ NonaddCS8.to_string a.final)
    | Nonadd16 a ->
            ("na16: " ^ NonaddCS16.to_string a.final)
    | Nonadd32 a ->
            ("na32: " ^ NonaddCS32.to_string a.final)
    | AddVec a ->
            AddCS.Vector.to_string a.final
    | AddGen a ->
            AddCS.General.to_string a.final
    | Sank a ->
            ("sank: " ^ SankCS.to_string a.final)
    | Dynamic a ->
            "dynamic: " ^ DynamicCS.to_string a.final
    | FixedStates a -> ("fixed states: " ^ Fixed_states.to_string a.final)
    | Set a ->
            let sub = List.map to_string_ch a.final.set in
            let stype = match a.final.smethod with
                | `Strictly_Same -> "same"
                | `Any_Of _ -> "any-of"
            in
            "set(" ^ stype ^ "): [" ^ (String.concat "; " sub) ^ "]"
    | Kolmo a ->
            ("kolmo: " ^ KolmoCS.to_string a.final)
    | StaticMl a ->
        IFDEF USE_LIKELIHOOD THEN
            ("static ML: " ^ MlStaticCS.to_string a.preliminary)
        ELSE
            failwith MlStaticCS.likelihood_error
        END

let extract_cost = function
    | Nonadd8 v -> v.cost
    | Nonadd16 v -> v.cost
    | Nonadd32 v -> v.cost
    | AddVec v -> v.cost
    | AddGen v -> v.cost
    | Sank v -> v.cost
    | FixedStates v -> v.cost
    | Dynamic v -> v.cost
    | Set v -> v.cost
    | Kolmo v -> v.cost
    | StaticMl v -> 
        IFDEF USE_LIKELIHOOD THEN
            v.cost
        ELSE
            failwith MlStaticCS.likelihood_error
        END

let get_characters_cost chars =
    List.fold_left (fun a b -> a +. (extract_cost b)) 0. chars

(** Helper function for recursing into sets *)
let setrec a fn = { a with set = List.map fn a.set }

(** [set_update_cost set] calculates the new cost of the set based on its
    contents, and updates the node's cost field.  Only applicable for
    [`Strictly_Same] sets. *)
let set_update_cost = function
    | Set a ->
        begin match a.preliminary.smethod with
           | `Strictly_Same ->
                let cost =
                    List.fold_left
                        (fun acc m -> acc +. extract_cost m)
                        0. a.preliminary.set in
                if debug_set_cost
                    then odebug ("`Strictly_Same cost: " ^ string_of_float cost);
                Set { a with cost = cost }, cost
            | `Any_Of _ -> assert false
        end
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
        cost_mode : [ `Likelihood | `NotLikelihood | `SumLikelihood  ];
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


let print_times n =
    let get_times = function 
        | Nonadd8 x  -> x.time
        | Nonadd16 x -> x.time
        | Nonadd32 x -> x.time
        | AddVec x   -> x.time
        | AddGen x   -> x.time
        | Sank x     -> x.time
        | Dynamic x  -> x.time
        | FixedStates x -> x.time
        | Set x      -> x.time
        | Kolmo x    -> x.time
        | StaticMl x ->
            IFDEF USE_LIKELIHOOD THEN
                x.time
            ELSE
                failwith MlStaticCS.likelihood_error
            END
    in
    List.iter
        (fun ncs ->
            let one,two,thr = get_times ncs in
            Printf.printf "%d -- (%s,%s|%s)\n%!" n.taxon_code
                (match one with | Some x -> string_of_float x | None -> "none")
                (match two with | Some x -> string_of_float x | None -> "none")
                (match thr with | Some x -> string_of_float x | None -> "none") )
        n.characters

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

let map5 f a b c d e=
    let rec mapper a b c d e acc =
        match a, b, c, d, e with
        | (ha :: ta), (hb :: tb), (hc :: tc), (hd :: td),(he :: te) ->
                mapper ta tb tc td te ((f ha hb hc hd he) :: acc)
        | [], [], [], [], [] -> List.rev acc
        | _ -> 
                raise (Illegal_argument "map5")
    in
    mapper a b c d e []

let mapN f in_lstlst =
    let len = 
         List.fold_left
            (fun last_lst_len this_lst -> 
                let this_lst_len = List.length this_lst in
                (* make sure each list of input listlist has the same length*)
                assert(this_lst_len = last_lst_len); 
                this_lst_len)
            (List.length (List.hd in_lstlst))
            (List.tl in_lstlst)
    in
    let res_lst = ref [] in
    for i = len-1 downto 0 do
        let f_input = List.map (fun lst -> List.nth lst i) in_lstlst in
        let res = f f_input in
        res_lst := res :: (!res_lst);
    done;
    !res_lst


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
                | Nonadd8 x  -> process_nonadd NonaddCS8.to_list x
                | Nonadd16 x -> process_nonadd NonaddCS16.to_list x
                | Nonadd32 x -> process_nonadd NonaddCS32.to_list x
                | AddVec x ->
                    let lst = AddCS.Vector.to_list_with_cost x.preliminary in
                    List.map (fun (_, _, code, cost) -> `Add, code, cost) lst
                | AddGen x ->
                    let lst = AddCS.General.to_list_with_cost x.preliminary in
                    List.map (fun (_, _, code, cost) -> `Add, code, cost) lst
                | Sank x ->
                    let lst = SankCS.to_list x.preliminary in
                    List.map
                        (fun (a, b) ->
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
    Printf.printf "{ taxon_code = %d, total_cost = %f, node_cost = %f \n%!"
    nd.taxon_code nd.total_cost nd.node_cost;
    let count = ref 0 in
    List.iter 
        (fun a_cs ->
            count := !count + 1;
            match a_cs with
            | Dynamic a_dyn ->
                Printf.printf "character#.%d, cost = %f, sum_cost = %f \n%!" 
                !count a_dyn.cost a_dyn.sum_cost;
                print_endline "Preliminary state";
                DynamicCS.print a_dyn.preliminary;
                print_endline "Final state";
                DynamicCS.print a_dyn.final;
            | _ -> print_endline "Do not print non-dynamic characters")
        nd.characters;
    Printf.printf " }\n%!"
    
let to_string {characters=chs; total_cost=cost; taxon_code=tax_code} =
    ("[[NODE tax_code=" ^ string_of_int tax_code
     ^ " total cost=" ^ string_of_float cost
     ^ " elts: "
     ^ (String.concat "; " (List.map to_string_ch chs))
     ^ "]]")

(*** Helper functions for node data ***)
let cs_prelim_to_final a =
    {a with final = a.preliminary}


let total_cost _ a = a.total_cost

let rec prelim_to_final = function
    | StaticMl a -> 
        IFDEF USE_LIKELIHOOD THEN
            StaticMl (cs_prelim_to_final a)
        ELSE
            failwith MlStaticCS.likelihood_error
        END
    | Nonadd8 a -> Nonadd8 (cs_prelim_to_final a)
    | Nonadd16 a -> Nonadd16 (cs_prelim_to_final a)
    | Nonadd32 a -> Nonadd32 (cs_prelim_to_final a)
    | AddVec a -> AddVec (cs_prelim_to_final a)
    | AddGen a -> AddGen (cs_prelim_to_final a)
    | Sank a -> Sank (cs_prelim_to_final a)
    | Dynamic a -> Dynamic (cs_prelim_to_final a)
    | FixedStates a -> FixedStates (cs_prelim_to_final a)
    | Kolmo a -> Kolmo (cs_prelim_to_final a)
    | Set a ->
          let r = setrec a.preliminary prelim_to_final in
          Set { a with preliminary = r; final = r; }

let all_prelim_to_final ({characters = chars} as node) =
    {node with
         characters = List.map prelim_to_final chars}

let using_likelihood types x =
    let x,y = 
        IFDEF USE_LIKELIHOOD THEN
            List.fold_left
                (fun ((static,dyn) as res) x -> match x with 
                    | StaticMl _ -> (true,dyn)
                    | Dynamic x ->
                        begin match x.preliminary with
                        | DynamicCS.MlCS _ -> (static,true)
                        | _ -> res
                    end
                    | _          -> res)
                (false,false) x.characters
        ELSE
            false,false
        END
    in
    match types with
        | `Static   -> x
        | `Dynamic  -> y
        | `Either   -> x || y

(* extracts data on whole character set *)
let extract_states alph data in_codes node =
    let alph = Alphabet.to_sequential alph in
    let include_cs = match in_codes with
        | Some in_codes ->
            (fun cs_codes ->    
                let r = ref false in
                for i = 0 to (Array.length in_codes) - 1 do begin
                    if !r then ()
                    else begin
                        for j = 0 to (Array.length cs_codes) - 1 do
                            r := !r || (in_codes.(i) = cs_codes.(j));
                        done;
                    end
                end done;
                !r)
        | None -> (fun _ -> true) (* all data *)
    and norm c = 
        List.map
            (fun x -> let str = Data.to_human_readable data c x in
                      Alphabet.find_base str alph)
    in
    let extract_states_cs = function
        | Nonadd8 a ->
            if not (include_cs a.final.NonaddCS8.codes) then []
            else List.map
                    (fun (c,e) -> a.weight, c, `List (norm c (NonaddCS8.e_to_list e)) )
                    (NonaddCS8.to_simple_list a.final)
        | Nonadd16 a ->
            if not (include_cs a.final.NonaddCS16.codes) then []
            else List.map
                    (fun (c,e) -> a.weight, c, `List (norm c (NonaddCS16.e_to_list e) ))
                    (NonaddCS16.to_simple_list a.final)
        | Nonadd32 a ->
            if not (include_cs a.final.NonaddCS32.codes) then []
            else List.map
                    (fun (c,e) -> a.weight, c, `List (norm c (NonaddCS32.e_to_list e)))
                    (NonaddCS32.to_simple_list a.final)
        | StaticMl a ->
            IFDEF USE_LIKELIHOOD THEN
                MlStaticCS.extract_states a.final
            ELSE
                failwith "Unsupported character type"
            END
        | _ -> failwith "Unsupported character type"
    in
    List.flatten (List.map extract_states_cs node.characters)


(*[cs_update_cost_only mine ch1 ch2] when any of the children has a different
* sum_cost, but the assignment for that children(from median function) remain
* the same, here we only need to update sum_cost of mine*)
let cs_update_cost_only mine ch1 ch2 = 
    let debug = false in
    match mine,ch1,ch2 with
    | StaticMl m, StaticMl c1, StaticMl c2 ->
        IFDEF USE_LIKELIHOOD THEN
            let sumcost = m.cost in
            StaticMl { m with sum_cost = sumcost; }, sumcost
        ELSE
            failwith MlStaticCS.likelihood_error
        END
    | Dynamic m, Dynamic c1, Dynamic c2 ->
            let sumcost = c1.sum_cost +. c2.sum_cost +. m.cost in
            if debug then Printf.printf "update_cost_only, sum_cost <- %f\n%!" sumcost;
            Dynamic { m with sum_cost = sumcost;}, sumcost
    | AddVec m, AddVec c1, AddVec c2 -> 
            let sumcost = c1.sum_cost +. c2.sum_cost +. m.cost in
            AddVec { m with sum_cost = sumcost; }, sumcost
    | AddGen m, AddGen c1, AddGen c2 -> 
            let sumcost = c1.sum_cost +. c2.sum_cost +. m.cost in
            AddGen { m with sum_cost = sumcost; }, sumcost
    | Kolmo m, Kolmo c1, Kolmo c2 ->  
            let sumcost = c1.sum_cost +. c2.sum_cost +. m.cost in
            Kolmo { m with sum_cost = sumcost; }, sumcost
    | Sank m, Sank c1, Sank c2 ->
            Sank m, m.sum_cost
    | FixedStates m, FixedStates c1, FixedStates c2 ->
            FixedStates m, m.sum_cost
    | Nonadd8 m, Nonadd8 c1, Nonadd8 c2 ->  
            let sumcost = c1.sum_cost +. c2.sum_cost +. m.cost in
            Nonadd8 { m with sum_cost = sumcost }, sumcost
    | Nonadd16 m, Nonadd16 c1, Nonadd16 c2 ->  
            let sumcost = c1.sum_cost +. c2.sum_cost +. m.cost in
            Nonadd16 { m with sum_cost = sumcost }, sumcost
    | Nonadd32 m, Nonadd32 c1, Nonadd32 c2 ->  
            let sumcost = c1.sum_cost +. c2.sum_cost +. m.cost in
            Nonadd32 { m with sum_cost = sumcost }, sumcost
    | _, _ , _ ->
            failwith "cs_update_cost_only, wrong character type mix"

(* calculate the median between two nodes *)
let rec cs_median code anode bnode prev t1 t2 a b =
    if debug_cs_median then Printf.printf "node.cs_median, on node %d and %d:\n%!"
	anode.taxon_code bnode.taxon_code;
    match a, b with
    | StaticMl ca, StaticMl cb ->
        IFDEF USE_LIKELIHOOD THEN
            assert (ca.weight = cb.weight);
            let min_bl = MlStaticCS.minimum_bl () in
            let t1, t2 = match t1,t2 with
                | Some (t1), Some (t2) when code <= 0 ->
                    ( max min_bl (t1/.2.0), max min_bl (t2/.2.0) )
                (* this situation happens on edge data for roots *)
                | (Some t1, None)
                | (None, Some t1) when code <= 0 ->
                    ( max min_bl (t1/.2.0), max min_bl (t1/.2.0) )
                | (Some t1, None)
                | (None, Some t1) ->
                    ( max min_bl (t1/.2.0), max min_bl (t1/.2.0) )
                | Some (t1), Some (t2) -> 
                    (max min_bl t1, max min_bl t2)
                | None, None ->
                    let t1,t2 = MlStaticCS.estimate_time ca.preliminary cb.preliminary in
                    if debug_bl then Printf.printf "estimating BL: %f, %f\n%!" t1 t2;
                    t1,t2
            in 
            let model = MlStaticCS.get_model ca.preliminary in
            let median = begin match model .MlModel.spec.MlModel.cost_fn with
                | `MPL when code > 0 -> 
                    MlStaticCS.median2 ca.preliminary cb.preliminary
                                       t1 t2 anode.taxon_code bnode.taxon_code
                | `MAL -> 
                    MlStaticCS.median2 ca.preliminary cb.preliminary
                                       t1 t2 anode.taxon_code bnode.taxon_code
                | `MPL ->
                    MlStaticCS.median1 ca.preliminary cb.preliminary (t1+.t2)
                end
            in
            let n_cost = MlStaticCS.root_cost median in
            if debug then 
                info_user_message
                    "Calculating %d with %f(%d) and %f(%d) = %f (static)"
                    code t1 anode.taxon_code t2 bnode.taxon_code n_cost;
            let t1,t2 = 
                if anode.min_child_code < bnode.min_child_code then (t1, t2)
                else (t2, t1)
            in
            let n_cost  =  n_cost *. ca.weight in
            let res =
                {
                    preliminary = median;
                    final = median;
                    cost = n_cost;
                    sum_cost = n_cost;
                    time = Some t1, Some t2, None;
                    weight = ca.weight;
                }
            in
            StaticMl res, n_cost
        ELSE
            failwith MlStaticCS.likelihood_error
        END
    | Nonadd8 ca, Nonadd8 cb ->
            assert (ca.weight = cb.weight);
            let prev = match prev with
                | None -> None
                | Some (Nonadd8 prev) -> Some (prev.preliminary)
                | _ -> raise (Illegal_argument "cs_median")
            in
            let median = NonaddCS8.median prev ca.preliminary cb.preliminary in
            let cost = NonaddCS8.median_cost median in
            let cost = cost *. ca.weight in
            let sumcost = ca.sum_cost +. cb.sum_cost +. cost in
            if debug_cs_median then 
                Printf.printf "node.ml cs_median Nonadd8 weight=%f, cost <- %f,\
                sumcost<-%f(l)+%f(r)+cost=%f\n%!" ca.weight cost ca.sum_cost
                cb.sum_cost sumcost;
            let res = { ca with preliminary = median; 
                            final = median; 
                            cost = cost;
                            sum_cost = sumcost;
                      } in
            Nonadd8 res, sumcost
    | Nonadd16 ca, Nonadd16 cb -> 
            assert (ca.weight = cb.weight);
            let prev = match prev with
                | None -> None
                | Some (Nonadd16 prev) -> Some (prev.preliminary)
                | _ -> raise (Illegal_argument "cs_median");
            in
            let median = NonaddCS16.median prev ca.preliminary cb.preliminary in
            let cost = NonaddCS16.median_cost median in
            let cost = cost *. ca.weight in
            let sumcost = ca.sum_cost +. cb.sum_cost +. cost in
            let res =
                { ca with
                    preliminary = median; 
                    final = median; 
                    sum_cost = sumcost;
                    cost = cost; }
            in
            Nonadd16 res, sumcost
    | Nonadd32 ca, Nonadd32 cb -> 
            assert (ca.weight = cb.weight);
            let prev = match prev with
                | None -> None
                | Some (Nonadd32 prev) -> Some (prev.preliminary)
                | _ -> raise (Illegal_argument "cs_median")
            in
            let median = NonaddCS32.median prev ca.preliminary cb.preliminary in
            let cost = NonaddCS32.median_cost median in
            let cost = cost *. ca.weight in
            let sumcost = ca.sum_cost +. cb.sum_cost +. cost in
            let res =
                { ca with
                    preliminary = median; 
                    final = median; 
                    sum_cost = sumcost;
                    cost = cost }
            in
            Nonadd32 res, sumcost
    | AddVec ca, AddVec cb -> 
            assert (ca.weight = cb.weight);
            let old = match prev with
                | Some AddVec old -> Some old.preliminary
                | None -> None
                | _ -> assert false in
            let median = AddCS.Vector.median old ca.preliminary cb.preliminary in
            let cost = AddCS.Vector.median_cost median in
            let cost = cost *. ca.weight in
            let sumcost = ca.sum_cost +. cb.sum_cost +. cost in
            let res =
                { ca with
                    preliminary = median; 
                    final = median; 
                    sum_cost = sumcost;
                    cost = cost; } 
            in
            AddVec res, sumcost
    | AddGen ca, AddGen cb -> 
            assert (ca.weight = cb.weight);
            let old = match prev with
                | Some AddGen old -> Some old.preliminary
                | None -> None
                | _ -> assert false in
            let median = AddCS.General.median old ca.preliminary cb.preliminary in
            let cost = AddCS.General.median_cost median in
            let cost = cost *. ca.weight in
            let sumcost = ca.sum_cost +. cb.sum_cost +. cost in
            if debug_cs_median then Printf.printf "node.cs_median,\
            AddGen,weight=%f, cost = %f, sumcost = %f+%f+cost = %f\n%!" 
            ca.weight cost ca.sum_cost cb.sum_cost
            sumcost;
            let res =
                { ca with
                    preliminary = median; 
                    final = median; 
                    sum_cost = sumcost;
                    cost = cost; } 
            in
            AddGen res, sumcost
    | Sank ca, Sank cb ->
            assert (ca.weight = cb.weight);
            let median,cost = SankCS.median code ca.preliminary cb.preliminary in
            let cost = cost *. ca.weight in
            let res = 
                { ca with
                    preliminary = median;
                    final = median;
                    sum_cost = cost;
                    cost = cost }
            in
            if debug_cs_median then Printf.printf "node.cs_median
            sankoff,weight=%f, cost=%f=sum cost\n %!" ca.weight cost;
            Sank res, cost
    | FixedStates ca, FixedStates cb ->
            let median,cost = Fixed_states.median code ca.preliminary cb.preliminary in
            let cost = cost *. ca.weight in
            let res = { 
                    ca with preliminary = median;
                    final = median;
                    sum_cost = cost;
                    cost = cost }
            in
            FixedStates res, cost
    | Dynamic ca, Dynamic cb ->
            assert (ca.weight = cb.weight);
            let ca, cb =
                if anode.min_child_code < bnode.min_child_code then ca, cb
                else cb, ca
            in
            let median, times = match ca.preliminary, cb.preliminary with
                | DynamicCS.MlCS ca_pre, DynamicCS.MlCS cb_pre ->
                    IFDEF USE_LIKELIHOOD THEN
                        let min_bl = MlStaticCS.minimum_bl () in
                        let t1, t2 = match t1,t2 with
                            | Some (t1), Some (t2) when code <= 0 -> max min_bl t1, 0.0
                            | Some (t1), Some (t2) -> max min_bl t1, max min_bl t2
                            | (Some t1, None)
                            | (None, Some t1) when code <= 0 -> max min_bl t1, 0.0
                            | _ when code <= 0 -> 
                                let t1,t2 = MlDynamicCS.estimate_time ca_pre cb_pre in
                                if debug_bl then Printf.printf "estimating BL: %f, %f\n%!" t1 t2;
                                (t1+.t2), 0.0
                            | _ ->
                                let t1,t2 = MlDynamicCS.estimate_time ca_pre cb_pre in
                                if debug_bl then Printf.printf "estimating BL: %f, %f\n%!" t1 t2;
                                t1,t2
                        in
                        let median =
                            MlDynamicCS.median code ca_pre cb_pre (Some t1) (Some t2)
                        in
                        let res = DynamicCS.MlCS median in
                        let cst = DynamicCS.total_cost res in
                        if debug then 
                            info_user_message
                                "Calculating %d with %f(%d) and %f(%d) = %f (dyn)"
                                code t1 anode.taxon_code t2 bnode.taxon_code cst;
                        let t1,t2 = 
                            if anode.min_child_code < bnode.min_child_code 
                            then t1,t2 else t2,t1
                        in
                        (res,(Some t1,Some t2,None))
                    ELSE
                        failwith MlStaticCS.likelihood_error
                    END
                | ca_prelim, cb_prelim ->
                    DynamicCS.median code ca.preliminary cb.preliminary None None,ca.time
            in
            let total_cost = ca.weight *. (DynamicCS.total_cost median) in
            let sumcost = ca.sum_cost +. cb.sum_cost +. total_cost in
            if debug_cs_median then
                info_user_message "Return Median with costs: Cost:%f,Sum_cost:%f(ca=%f,cb=%f)" 
                total_cost sumcost ca.sum_cost cb.sum_cost;
            let res = 
                { ca with 
                    preliminary = median;
                    final = median;
                    cost = total_cost;
                    sum_cost = sumcost;
                    time = times;
                } 
            in
            if debug_cs_median then Printf.printf "end of node.cs_median, dynamic, \
            weight=%f,cost=%f,sum_cost=%f\n%!" ca.weight total_cost sumcost;
            Dynamic res, sumcost
    | Kolmo ca, Kolmo cb ->
            assert (ca.weight = cb.weight);
            let ca, cb =
                if anode.min_child_code < bnode.min_child_code then ca, cb
                else cb, ca
            in
            let median = KolmoCS.median code ca.preliminary cb.preliminary in
            let total_cost = KolmoCS.total_cost median in
            let total_cost = ca.weight *. total_cost in
            let sumcost = ca.sum_cost +. cb.sum_cost +. total_cost in
            let res = 
                { ca with 
                    preliminary = median;
                    final = median;
                    cost = total_cost;
                    sum_cost = sumcost;
                } 
            in
            Kolmo res, sumcost
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
                    let cs_lst,_ (*ignore cost_list*)= List.split (
                        List.map2 (cs_median code anode bnode None t1 t2) l1 l2 ) 
                    in
                    let result = { cb.preliminary with set = cs_lst } in
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
                                      let median,_(*ignore sumcost*) =  
                                           cs_median code anode bnode
                                          None t1 t2 l1i l2i
                                      in
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
                            medians
                    in
                    (* Make the proper sets *)
                    let median =
                        { cb.preliminary with
                              set = List.map
                                (fun (_, median, _, _) -> median) medians;
                              smethod =
                                `Any_Of ((anode.taxon_code, bnode.taxon_code,
                                    (List.map (fun (_, _, i, j) -> (i, j)) medians)), v);
                        }
                    in
                    (* Update the cost and return *)
                    let cost = ca.weight *.  
                    (if max_float = !min_cost then 0. else !min_cost) in
                    let sumcost = ca.sum_cost +. cb.sum_cost +. cost in
                    let res =
                        Set { cb with
                                preliminary = median;
                                final = median;
                                cost = cost;
                                sum_cost = sumcost;
                            }
                    in
                    res, sumcost
          end
    | Nonadd8 _, _ | Nonadd16 _, _| Nonadd32 _, _ | AddVec _, _ | Sank _, _  
    | FixedStates _, _ | Dynamic _, _ |  Set _, _ | Kolmo _, _ | StaticMl _, _ | AddGen _,_ -> 
            raise (Illegal_argument "cs_median")

(** [edge iterator gp rent c1 c2]
 * iterates the branch lengths in the character set that are ML characters. *) 
let edge_iterator (gp:node_data option) (c0:node_data) (c1:node_data) (c2:node_data) = 
    IFDEF USE_LIKELIHOOD THEN
    (* preliminary: check if we are adjusting a handle/root *)
    let root_e = match gp with 
        | Some x -> false
        | None -> true
    and c1,c2 =
        if c1.min_child_code < c2.min_child_code then c1,c2 else c2,c1
    in
    (* accumulators [pa],[aa],[ba] to iterate the branch length for each element in
     * the character set [p] as parent, and [a], [b] as children. *)
    let rec ei_map p a b pa = match p,a,b with
        | (StaticMl pm)::ptl, (StaticMl am)::atl,(StaticMl bm)::btl ->
            let model = MlStaticCS.get_model am.preliminary in
            let min_bl = MlStaticCS.minimum_bl () in
            let mine,sumcost = begin match  model.MlModel.spec.MlModel.cost_fn with
                (* MAL has full C functionality *)
                | `MAL -> 
                    let modf = ref All_sets.Integers.empty in
                    let t1,t2,t3_opt = match pm.time with 
                        | Some x, Some y, z -> ( max min_bl x, max min_bl y, z)
                        | None, None, None -> 
                            let (x,y) = MlStaticCS.estimate_time am.preliminary bm.preliminary in
                            (x, y, None )
                        | _ -> failwith "something happened terribly wrong"
                    in
                    let t1,t2 = if root_e then (t2 +. t1),0.0 else t1,t2 in
                    let _,pcost,cost,(t1,t2),res = 
                        MlStaticCS.readjust None !modf am.preliminary bm.preliminary
                                            pm.preliminary t1 t2 
                    in
                    let cost = cost *. pm.weight in
                    { pm with  preliminary = res; final = res;
                               cost = cost; sum_cost = cost;
                               time = Some t1, Some t2, t3_opt; },
                    cost
                (* calculate the median1 for MPL with OCAML Brents *)
                | `MPL ->
                    let calculate_single t = 
                        let m = MlStaticCS.median1 am.preliminary bm.preliminary t in
                        (m,MlStaticCS.root_cost m)
                    and t1,t2 = match pm.time with 
                        | Some x,Some y, _ -> ( max min_bl x, max min_bl y )
                        | None, None, _ -> 
                            let (x,y) = MlStaticCS.estimate_time am.preliminary bm.preliminary in
                            if debug_bl then Printf.printf "estimating BL: %f, %f\n%!" x y;
                            (x, y)
                        | _ -> failwith "something happened terribly wrong"
                    in
                    let fstart = pm.preliminary,MlStaticCS.root_cost pm.preliminary in
                    let (v,(dv,fv)) =
                        Numerical.brents_method calculate_single (t1+.t2,fstart)
                    in
                    let cost = fv *. bm.weight in
                    { pm with preliminary = dv; final = dv;
                              cost = cost; sum_cost = cost;
                              time = Some v, Some 0.0, None; },
                    cost
                end
            in
            ei_map ptl atl btl ((StaticMl mine, sumcost)::pa)
        (* dynamic characters *)
        | ((Dynamic pm) as pml)::ptl, (Dynamic am)::atl,(Dynamic bm)::btl ->
            begin match am.preliminary, bm.preliminary with
                | DynamicCS.MlCS c1_pre, DynamicCS.MlCS c2_pre ->
                    let min_bl = MlStaticCS.minimum_bl () in
                    let modf = ref All_sets.Integers.empty in
                    let t1,t2,t3_opt = match pm.time with 
                        | Some x,Some y, z -> ( max min_bl x, max min_bl y,z )
                        | None, None, z -> 
                            let (x,y) = MlDynamicCS.estimate_time c1_pre c2_pre in
                            if debug_bl then Printf.printf "estimating BL: %f, %f\n%!" x y;
                            ( x, y, z )
                        | _ -> failwith "something happened terribly wrong"
                    in
                    if debug then
                        info_user_message "Starting Cost: %f\t%f\t%f%!" t1 t2 pm.cost; 
                    let ot1,ot2 = if root_e then (t2 +. t1),0.0 else t1,t2 in
                    let m,pcost,cost,(t1,t2),res = 
                        DynamicCS.readjust_lk (`ThreeD None) None !modf
                            am.preliminary bm.preliminary pm.preliminary ot1 ot2
                    in
                    let cost = cost *. am.weight in
                    if debug then
                        info_user_message "Ending Cost: %f\t%f\t%f%!" t1 t2 cost;
                    let sumcost = cost +. bm.sum_cost +. am.sum_cost in
                    let mine =  
                        { pm with  preliminary = res; final = res;
                                   cost = cost;
                                   sum_cost = sumcost;
                                   time = Some t1, Some t2, t3_opt; }
                    in
                    ei_map ptl atl btl ((Dynamic mine, sumcost)::pa)
                | _,_ ->
                    let old_sumcost = pm.sum_cost in
                    ei_map ptl atl btl ((pml,old_sumcost)::pa)
            end
        (* ignore non-likelihood characters *)
        (* | Nonadd8 | Nonadd16 | Nonadd32 | AddVec | AddGen | Sank | FixedStates |Dynamic | Kolmo | Set *)
        | ((Nonadd8 pm) as pml)::ptl, aml::atl, bml::btl ->
                ei_map ptl atl btl ((pml,pm.sum_cost)::pa)
        | ((Nonadd16 pm) as pml)::ptl, aml::atl, bml::btl ->
                ei_map ptl atl btl ((pml,pm.sum_cost)::pa)
        | ((Nonadd32 pm) as pml)::ptl, aml::atl, bml::btl ->
                ei_map ptl atl btl ((pml,pm.sum_cost)::pa)
        | ((AddVec pm) as pml)::ptl, aml::atl, bml::btl ->
                ei_map ptl atl btl ((pml,pm.sum_cost)::pa)
        | ((AddGen pm) as pml)::ptl, aml::atl, bml::btl ->
                ei_map ptl atl btl ((pml,pm.sum_cost)::pa)
        | ((Sank pm) as pml)::ptl, aml::atl, bml::btl ->
                ei_map ptl atl btl ((pml,pm.sum_cost)::pa)
        | ((FixedStates pm) as pml)::ptl, aml::atl, bml::btl ->
                ei_map ptl atl btl ((pml,pm.sum_cost)::pa)
        | ((Dynamic pm) as pml)::ptl, aml::atl, bml::btl ->
                ei_map ptl atl btl ((pml,pm.sum_cost)::pa)
        | ((Kolmo pm) as pml)::ptl, aml::atl, bml::btl ->
                ei_map ptl atl btl ((pml,pm.sum_cost)::pa)
        | ((Set pm) as pml)::ptl, aml::atl, bml::btl ->
                ei_map ptl atl btl ((pml,pm.sum_cost)::pa)
        | [],[],[] -> pa
        | _ -> failwith "Number of characters is inconsistent"
    in
    let mine, sumcost_list =
        List.split ( List.rev (ei_map c0.characters c1.characters c2.characters
        []) )
    in
    let mine_cost = get_characters_cost mine in
    let total_cost = List.fold_left (fun acc x -> acc +. x ) 0. sumcost_list in 
    { c0 with characters = mine;
              total_cost = total_cost;
              node_cost  = mine_cost;
    }
    ELSE
        c0
    END

let apply_time root child parent =
    IFDEF USE_LIKELIHOOD THEN
        let fst (a,_,_) = a and snd (_,a,_) = a and thr (_,_,a) = a in
        let p = if child.min_child_code = parent.min_child_code then fst else snd in
        let p_opt to_string = function | None -> "None" | Some x -> to_string x in
        let apply_it p pndtime = 
            let res = if root then begin
                        match fst pndtime, snd pndtime with
                        | Some x, Some y -> 
                            if debug then
                                info_user_message "Adding to %f -- %f%!" x y;
                            Some (x +. y)
                        | _,_ -> failwith "Incorrect branch lengths in a root node"
                      end else p pndtime
            in
            if debug then
                info_user_message "Applying %s (root:%b) to %d -- %d%!"
                                  (p_opt string_of_float res) root
                                  (child.taxon_code) (parent.taxon_code);
            res
        and replace_third (a,b,_) x = (a,b,x) in
        let rec apply_times ch par =
            match ch,par with
            | StaticMl cnd,StaticMl pnd ->
                let time = replace_third cnd.time (apply_it p pnd.time) in
                if debug then begin
                    info_user_message "%d has %s | %s | %s" (parent.taxon_code)
                                  (p_opt string_of_float (fst time))
                                  (p_opt string_of_float (snd time))
                                  (p_opt string_of_float (thr time));
                    info_user_message "%d has %s | %s | %s" (child.taxon_code)
                                  (p_opt string_of_float (fst time))
                                  (p_opt string_of_float (snd time))
                                  (p_opt string_of_float (thr time))
                end;
                StaticMl { cnd with time = time }
            | Dynamic cnd, Dynamic pnd -> 
                begin match cnd.preliminary,pnd.preliminary with
                | DynamicCS.MlCS _,DynamicCS.MlCS _ ->
                    if debug then begin
                        info_user_message "%d has %s | %s | %s" (parent.taxon_code)
                                      (p_opt string_of_float (fst pnd.time))
                                      (p_opt string_of_float (snd pnd.time))
                                      (p_opt string_of_float (thr pnd.time));
                        info_user_message "%d has %s | %s | %s" (child.taxon_code)
                                      (p_opt string_of_float (fst cnd.time))
                                      (p_opt string_of_float (snd cnd.time))
                                      (p_opt string_of_float (thr cnd.time));
                    end;
                    let time = replace_third cnd.time (apply_it p pnd.time) in
                    Dynamic { cnd with time = time }
                | _,_ -> ch
                end
            | _ -> ch
        in
        { child with characters =
                 List.map2 apply_times child.characters parent.characters }
    ELSE
        child
    END

let compare_opt a b f1 f2 one two = match one,two with
    | Some x,Some y ->
            info_user_message "Checking times: %f(%d)(%s) == %f(%d)(%s)??"
                    x a.taxon_code (f1 ("fst","snd")) y b.taxon_code 
                    (f2 ("fst","snd")) ; x = y
    | None, None ->
            info_user_message "Checking times: None(%d)(fst) == None(%d)(fst)??"
                    a.taxon_code b.taxon_code; true
    | Some x,None ->
            info_user_message "Checking times: %f(%d)(%s) == None(%d)(fst)??"
                    x a.taxon_code (f1 ("fst","snd")) b.taxon_code; false
    | None,Some y ->
            info_user_message "Checking times: None(%d)(fst) == %f(%d)(%s)??"
                    a.taxon_code y b.taxon_code (f2 ("fst","snd")); false

let rec cs_final_states pn nn c1n c2n p n c1 c2 =
    match p, n, c1, c2 with
    | StaticMl cp, StaticMl cn, StaticMl cc1, StaticMl cc2 -> 
        n
    | Nonadd8 cp, Nonadd8 cn, Nonadd8 cc1, Nonadd8 cc2 -> 
        let m = NonaddCS8.median_3 cp.final cn.preliminary cc1.preliminary cc2.preliminary in
        Nonadd8 { cn with final = m }
    | Nonadd16 cp, Nonadd16 cn, Nonadd16 cc1, Nonadd16 cc2 -> 
        let m = NonaddCS16.median_3 cp.final cn.preliminary cc1.preliminary cc2.preliminary in
        Nonadd16 { cn with final = m }
    | Nonadd32 cp, Nonadd32 cn, Nonadd32 cc1, Nonadd32 cc2 -> 
        let m = NonaddCS32.median_3 cp.final cn.preliminary cc1.preliminary cc2.preliminary in
        Nonadd32 { cn with final = m }
    | AddVec cp, AddVec cn, AddVec cc1, AddVec cc2 -> 
        let m = AddCS.Vector.median_3 cp.final cn.preliminary cc1.preliminary cc2.preliminary in
        AddVec { cn with final = m }
    | AddGen cp, AddGen cn, AddGen cc1, AddGen cc2 -> 
        let m = AddCS.General.median_3 cp.final cn.preliminary cc1.preliminary cc2.preliminary in
        AddGen { cn with final = m }
    | Sank cp, Sank cn, Sank cc1, Sank cc2 ->
        let m = SankCS.median_3 cp.final cn.preliminary cc1.preliminary cc2.preliminary in
        Sank { cn with final = m }
    | FixedStates cp, FixedStates cn, FixedStates cc1, FixedStates cc2 ->
        let m = Fixed_states.median_3 cp.final cn.preliminary cc1.preliminary cc2.preliminary in
        FixedStates { cn with final = m }
    | Dynamic cp, Dynamic cn, Dynamic cc1, Dynamic cc2 ->
        IFDEF USE_LIKELIHOOD THEN
            begin match cp.final with
            | DynamicCS.MlCS _ ->
                begin match cn.time with
                    | Some x,Some y,Some z ->
                        let res = DynamicCS.final_states cn.preliminary
                                    cc1.preliminary cc2.preliminary cp.final x y z
                        in
                        Dynamic { cn with final = res; }
                    | None  ,None  ,Some z -> n (* leaf node *)
                    | _  -> print_times nn; assert false
                end
            | _ -> 
                let m = DynamicCS.median_3 cp.final cn.preliminary cc1.preliminary cc2.preliminary in
                Dynamic { cn with final = m }
            end
        ELSE
            let m = DynamicCS.median_3 cp.final cn.preliminary cc1.preliminary cc2.preliminary in
            Dynamic { cn with final = m }
        END
    | Kolmo cp, Kolmo cn, Kolmo cc1, Kolmo cc2 ->
        let m = KolmoCS.median_3 cp.final cn.preliminary cc1.preliminary cc2.preliminary in
        Kolmo { cn with final = m }
    | Set cp, Set cn, Set cc1, Set cc2 ->
        begin match cn.preliminary.smethod with
            | `Strictly_Same ->
                (* Apply in order to our children *)
                let res =
                    map4 (cs_final_states pn nn c1n c2n)
                         cp.preliminary.set cn.preliminary.set
                         cc1.preliminary.set cc2.preliminary.set
                in
                (* Note: if anything else updates the cost, the cost should
                 * also be recalculated here.  (But this shouldn't happen in
                 * the 3-median...) *)
                let res = { cn.preliminary with set = res } in
                Set { cn with preliminary = res; final = res; }
            | `Any_Of ((froml, fromr, medians), v) ->
                (* Get the "from" value for a node *)
                let get_from {smethod=m} = match m with
                    | `Any_Of (l, _) -> l
                    | _ -> assert false
                in
                 (* Get the from list specifically *)
                let get_from_list n =
                     let (_, _, l) = get_from n in l in
                 (* Am I the left node? *)
                let (parent_left, _, _) = get_from cp.final in
                let is_left = nn.taxon_code = parent_left in
                let (_, parent_right, _) = get_from cp.final in
                let is_right = nn.taxon_code = parent_right in
                assert (is_left || is_right);
                (* Make a list: parent median and which of my children it came from *)
                let set =
                    map2 (fun median (l, r) -> if is_left then (median, l) else (median, r))
                         cp.final.set
                         (get_from_list cp.final)
                in
                let set, from = match medians with
                    (* we're a leaf node (OTU); this is a special case for 3-median *)
                    | [] ->
                        (List.map
                            (fun (parent, my_index) ->
                                let median = List.nth cn.preliminary.set  my_index in
                                let left, right = median, median in
                                cs_final_states pn nn c1n c2n parent median left right)
                            set,
                        [])
                     (* real 3-median *)
                    | _ ->
                        List.split
                            (List.map
                                (fun (parent, my_index) ->
                                    let median = List.nth cn.preliminary.set  my_index in
                                    let (l, r) = List.nth medians my_index      in
                                    let left =   List.nth cc1.preliminary.set l in
                                    let right =  List.nth cc2.preliminary.set r in
                                    (cs_final_states pn nn c1n c2n parent median left right, (l, r)))
                                set)
                in
                let res =
                    { cn.preliminary with
                        set = set;
                        smethod = `Any_Of ((froml, fromr, from), v); }
                in
                Set { cn with final = res }
        end
    | Nonadd8 _, _, _, _ | Nonadd16 _, _, _, _ | Nonadd32 _, _, _, _ 
    | AddGen _, _, _, _  | Sank _, _, _, _     | FixedStates _, _, _, _ | Dynamic _, _, _, _  
    | AddVec _, _, _, _  | Set _, _, _, _      | Kolmo _, _, _, _ 
    | StaticMl _, _, _ ,_ ->
        raise (Illegal_argument "cs_final_states")

let new_node_stats a b =
    let num_child_edges = a.num_child_edges + b.num_child_edges + 2 in
    let num_height = (max a.num_height b.num_height) + 1 in
    (num_child_edges, num_height)
    
let excludes_median a b =
    let exclude_median = function
        | `Excluded, `Excluded -> 0, `Excluded
        | `NotExcluded, `NotExcluded -> 0, `NotExcluded
        | `Either, `Either -> 0, `Either
        | `Either, a -> 0, a
        | a, `Either -> 0, a
        | `Excluded, `NotExcluded
        | `NotExcluded, `Excluded -> 1, `Either
    in
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

let root_cost root =
    let adder acc character = 
        match character with
        | Kolmo v -> 
                acc +. KolmoCS.root_cost v.preliminary root.num_otus
                root.num_child_edges
        | _ -> acc
    in
    List.fold_left adder 0. root.characters

let get_characters_of_type map t node =
    List.map map 
    (List.filter (function
        | Nonadd8 _ when t = `NonAdd8 -> true
        | Nonadd16 _ when t = `NonAdd8 -> true
        | Nonadd32 _ when t = `NonAdd8 -> true
        | AddGen _ when t = `Add -> true
        | AddVec _ when t = `Add -> true
        | Sank _ when t = `Sank -> Printf.printf "get_characters_of_type,is Sank\n%!"; true
        | FixedStates _ when t = `FixedState -> Printf.printf "get_characters_of_type,is FixedState\n%!"; true
        | StaticMl _ when t = `StaticML -> true
        | Dynamic _ when t = `Dynamic -> true
        | Set _ when t = `Set -> true
        | Dynamic _ when t = `Ml -> true (* TODO *)
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
let get_addgen = 
    get_characters_of_type 
    (function AddGen x -> x.preliminary, x.final | _ -> assert false)
    `AddGen
let get_addvec = 
    get_characters_of_type 
    (function AddVec x -> x.preliminary, x.final | _ -> assert false)
    `AddVec
let get_sank = 
    get_characters_of_type 
    (function Sank x -> x.preliminary, x.final | _ -> assert false)
    `Sank
let get_fixedstates = 
    get_characters_of_type 
    (function FixedStates x -> x.preliminary, x.final | _ -> assert false)
    `FixedStates
let get_dynamic = 
    get_characters_of_type 
    (function Dynamic x -> x.preliminary, x.final | _ -> assert false)
    `Dynamic

IFDEF USE_LIKELIHOOD THEN
    let get_mlstatic = 
        get_characters_of_type 
            (function StaticMl x -> x.preliminary, x.final | _ -> assert false)
            `StaticML

    and get_mldynamic = 
        get_characters_of_type 
            (function StaticMl x -> x.preliminary, x.final | _ -> assert false)
            `Ml
ELSE
    let get_mlstatic x = []
    and get_mldynamic x = []
END

let get_set = 
    get_characters_of_type 
    (function Set x -> x.preliminary, x.final | _ -> assert false)
    `Set

let median_counter = ref (-1)

(* Convert a hashtable representation to list representation of branch lengths;
   The branch length information has character codes as hash keys that we use to
   associate to the character sets. *)
let convert_2_lst chars tbl : float option list =
    List.map
        (function
            | Dynamic z -> IFDEF USE_LIKELIHOOD THEN
                begin match z.preliminary with
                | DynamicCS.MlCS zdat ->
                    begin match tbl with
                     | Some x ->
                        begin try
                            let x = Hashtbl.find x chars.taxon_code in
                            let codes = MlDynamicCS.get_codes zdat in
                            Some (Hashtbl.find x codes.(0))
                        with | Not_found -> None end
                     | None -> None
                    end
                | _ -> None
                end
                ELSE
                    None
                END
            | StaticMl z -> 
                IFDEF USE_LIKELIHOOD THEN
                  begin match tbl with
                    | Some x ->
                        begin try
                            let x = Hashtbl.find x chars.taxon_code in
                            let codes = MlStaticCS.get_codes z.preliminary in
                            Some (Hashtbl.find x codes.(0))
                        with | Not_found -> None end
                    | None -> None
                  end
                ELSE
                    None
                END
            | _ -> None)
        chars.characters


(** [update_cost_only mine child1 child2] update mine with new sum_cost,
* calculated by new sum_cost of child1 or/and child2*)
let update_cost_only mine child1 child2 = 
    let debug = false in
    let new_exclude_info = 
        excludes_median child1 child2
    in
    let new_excluded = has_excluded new_exclude_info in
    if new_excluded then (*only update total_cost to inf*)
        { mine with total_cost = infinity; exclude_info = new_exclude_info; }
    else (*update sum_cost of each character, then the total_cost*)
        let characters_with_new_cost, sumcost_list = List.split (
            map3 (fun a b c -> cs_update_cost_only a b c) mine.characters child1.characters child2.characters )
        in
        let total_cost = List.fold_left (fun acc x -> acc +. x ) 0. sumcost_list in
        let _ = 
            if debug then Printf.printf 
            "update total_cost=%f only to node#.%d (old cost = %f, ch1(%d,%f),ch2(%d,%f)\n%!"
            total_cost mine.taxon_code mine.total_cost child1.taxon_code child1.total_cost
            child2.taxon_code child2.total_cost;
        in
        {mine with 
            characters = characters_with_new_cost; 
            total_cost = total_cost;
            exclude_info = new_exclude_info;
        }


let median ?branches code old a b =
    (* the code is negative if we are calculating on an edge *)
    let code = match code with
        | Some code -> code
        | None -> decr median_counter; !median_counter
    in
    if  debug_treebuild then info_user_message "node.ml median,nodea:%d(%f),nodeb:%d(%f),nodeab:%d"
    a.taxon_code a.total_cost b.taxon_code b.total_cost code;
    (* if code>0 then assert(false); *)
    let brancha = convert_2_lst a branches
    and branchb = convert_2_lst b branches in
    let new_characters,sumcost_list =
        List.split (
        match old with
        | None -> 
            if debug_treebuild then info_user_message "node.ml median begin of map4";
            map4 (cs_median code a b None) brancha branchb 
                 a.characters b.characters
        | Some c ->
            if debug_treebuild then info_user_message "node.ml median begin of map5";
            map5 (fun x -> cs_median code a b (Some x)) c.characters
                 brancha branchb a.characters b.characters
        )
    in
    let total_cost = List.fold_left (fun acc x -> acc +. x ) 0. sumcost_list in 
    let node_cost  = get_characters_cost new_characters in
    if debug_treebuild then
         info_user_message "end of mapx in node.ml,return nodedata with node_cost=%f,total_cost=%f" 
	node_cost total_cost;
    let num_child_edges, num_height = new_node_stats a b in
    let exclude_info = excludes_median a b in
    let excluded = has_excluded exclude_info in
    let results = 
        { 
            characters = new_characters;
            total_cost = if excluded then infinity else total_cost;
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
    in
    results

(** [get_times_between nd child_code] returns a list of times between [nd] and
 * the child node. data must be contained already in [nd] *)
let get_times_between_dir (parent) (nd:node_data) (child_code : int option) =
    let func =
IFDEF USE_LIKELIHOOD THEN
        let fst (a,_,_) = a and snd (_,a,_) = a and thr (_,_,a) = a in
        let f = match child_code with
            | _ when parent   -> thr
            | Some child_code ->
                let f = if nd.min_child_code = child_code then fst else snd in
                (* If only one direction is set; we can safely assume that we
                 * are dealing with a leaf. This avoids the unrulely circumstance
                 * that times could be located in either the first or second
                 * element; now it is always the first. *)
                (function | (Some _) as x, None, _ -> x
                          | x -> f x)
            | None ->
                (function | Some x,Some y, _ -> Some (x+.y)
                          | _ -> None)
        in
        function | StaticMl z -> f z.time
                 | Dynamic z ->
                    begin match z.preliminary with
                        | DynamicCS.MlCS _ -> f z.time
                        | _ -> None
                    end
                 | _ -> None
ELSE
        fun _ -> None
END
    in
    List.map func nd.characters

let get_times_between n c = get_times_between_dir false n c


let get_times_between_plus_codes (child:node_data) (parent:node_data option) =
(*    Printf.printf "Finding time between: %s -- %s : %d : "*)
(*        (string_of_int child.taxon_code)*)
(*        (match parent with | Some x -> string_of_int x.taxon_code | None -> "none")*)
(*        (List.length child.characters);*)
    let null = ([||],None) in
    let func =
IFDEF USE_LIKELIHOOD THEN
        let fst (a,_,_) = a and snd (_,a,_) = a and thr (_,_,a) = a in
        let f = match parent with
            | None     -> thr
            | Some par ->
                if par.min_child_code = child.min_child_code 
                    then fst
                    else snd
        in
        function
            | StaticMl z ->
                (*Printf.printf " : %a,%a,%a\n" pp_fopt (fst z.time) pp_fopt
                                  (snd z.time) pp_fopt (thr z.time);    *)
                MlStaticCS.get_codes z.preliminary, f z.time
            | Dynamic z ->
                (*Printf.printf " : %a,%a,%a\n" pp_fopt (fst z.time) pp_fopt
                                  (snd z.time) pp_fopt (thr z.time);    *)
                begin match z.preliminary with
                    | DynamicCS.MlCS x -> MlDynamicCS.get_codes x, f z.time
                    | _ -> null
                end
            | _ -> null
ELSE
        fun _ -> null
END
    in match parent with
        | Some par -> List.map func par.characters
        | None     -> List.map func child.characters

(** replace the third component of the time with the value passed *)
let replace_parent_time node time =
    IFDEF USE_LIKELIHOOD THEN
        let replace_third (a,b,_) x = (a,b,x) in
        let cs_replace node time = match node with
            | StaticMl a -> 
                StaticMl { a with time = replace_third a.time time; }
            | Dynamic a ->
                Dynamic { a with time = replace_third a.time time; }
            | x -> x (* we don't use time in other characters *)
        in
        { node with
            characters = List.map2 cs_replace node.characters time; }
    ELSE
        node
    END

(** [median_w_times gp nd1 nd2 t1 t2]
 * uses time data from two correct nodes, [time_1] and [time_2] for the
 * calculation of the median between [nd1] and [nd2]. **)
let median_w_times code prev nd_1 nd_2 times_1 times_2 times_3 =
    let code = match code with
        | Some code -> code
        | None      -> decr median_counter;
                       !median_counter
    in
    if debug_treebuild then
        Printf.printf "node.ml median_w_times, nd1 = %d,nd2 = %d,nd12=%d\n%!"
        nd_1.taxon_code nd_2.taxon_code code;
    let new_characters, sumcost_list = 
        List.split (
        match prev with
        | Some prev ->
            map5 (fun x -> cs_median code nd_1 nd_2 (Some x))
                 prev.characters times_1 times_2 nd_1.characters nd_2.characters
        | None ->
            map4 (cs_median code nd_1 nd_2 None)
                 times_1 times_2 nd_1.characters nd_2.characters
        )
    in
    let node_cost = get_characters_cost new_characters in
    let total_cost = List.fold_left (fun acc x -> acc +. x ) 0. sumcost_list in
    let num_child_edges, num_height = new_node_stats nd_1 nd_2
    and exclude_info = excludes_median nd_1 nd_2 in
    let excluded = has_excluded exclude_info in
    if debug_treebuild then
         Printf.printf "end of median_w_times, node.ml.tottal_cost=%f\n\n%!" total_cost;
    let node = 
      { characters = new_characters;
        total_cost = if excluded then infinity else total_cost;
        node_cost = node_cost;
        taxon_code = code;
        min_child_code = min nd_1.min_child_code nd_2.min_child_code;
        num_child_edges = num_child_edges;
        num_height = num_height;
        num_otus = nd_1.num_otus + nd_2.num_otus;
        exclude_sets = nd_1.exclude_sets;
        exclude_info = exclude_info;
        cost_mode = nd_1.cost_mode;
      }
    in
    match times_3 with
        | Some times_3 -> replace_parent_time node times_3
        | None         -> node

let median_of_child_branch code child parent =
    let times = get_times_between_dir true child None in
    let otime = List.map (fun _ -> None) times in
    median_w_times code None child parent times otime None

let final_states p n c1 c2 =
    let debug = false in
    if debug then begin
        Printf.printf "node.ml final_states, parent is:\n%!";
        print p;
        Printf.printf "mine is :\n%!";
        print n;
        Printf.printf "child 1 is :\n%!";
        print c1;
        Printf.printf "child 2 is :\n%!";
        print c2;
        Printf.printf "call median3 to assign final states\n%!";
    end;
    let new_characters =
        map4 (cs_final_states p n c1 c2)
        p.characters n.characters c1.characters c2.characters
    in { n with characters = new_characters }

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

let prior c n = 
IFDEF USE_LIKELIHOOD THEN
    let priors = 
(*        match c with*)
(*        | None ->*)
(*            (fun acc -> function*)
(*                | Dynamic a ->*)
(*                    begin match a.preliminary with*)
(*                    | DynamicCS.MlCS x -> MlDynamicCS.prior x*)
(*                    | _                -> acc*)
(*                    end*)
(*                | _ -> acc)*)
(*        | Some c ->*)
            (fun acc -> function
                | Dynamic a when DynamicCS.mem c a.final ->
                    begin match a.preliminary with
                    | DynamicCS.MlCS x -> MlDynamicCS.prior x
                    | _                -> acc
                    end
                | _ -> acc)
                
    in
    List.fold_left (priors) 0.0 n.characters
ELSE
    0.0
END


let get_cost_mode a = a.cost_mode

let classify_data leafa dataa leafb datab chars acc =
    let rec classify_two acc ch1 ch2 = match ch1, ch2 with
        | Dynamic a, Dynamic b ->
              DynamicCS.classify_transformations leafa a.final leafb b.final chars acc
        | (Kolmo _ | StaticMl _ | Nonadd8 _ | Nonadd16 _ | Nonadd32 _
          | AddGen _ | AddVec _ | FixedStates _ | Sank _ | Set _ | Dynamic _ ),_ ->
            assert false
    in
    List.fold_left2 classify_two acc dataa.characters datab.characters

let compare_data_final {characters=chs1} {characters=chs2} =
    let rec compare_two ch1 ch2 =
        match ch1, ch2 with
        | Nonadd8 a, Nonadd8 b   -> NonaddCS8.compare_data a.final b.final
        | Nonadd16 a, Nonadd16 b -> NonaddCS16.compare_data a.final b.final
        | Nonadd32 a, Nonadd32 b -> NonaddCS32.compare_data a.final b.final
        | AddVec a, AddVec b     -> AddCS.Vector.compare_data a.final b.final
        | AddGen a, AddGen b     -> AddCS.General.compare_data a.final b.final
        | Sank a, Sank b         -> SankCS.compare_data a.final b.final
        | FixedStates a, FixedStates b -> Fixed_states.compare_data a.final b.final
        | Dynamic a, Dynamic b   -> DynamicCS.compare_data a.final b.final
        | Kolmo a, Kolmo b       -> KolmoCS.compare_data a.final b.final
        | StaticMl a, StaticMl b ->
            IFDEF USE_LIKELIHOOD THEN
                MlStaticCS.compare_data a.final b.final
            ELSE
                failwith MlStaticCS.likelihood_error  
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
        | Nonadd8 a, Nonadd8 b   -> NonaddCS8.compare_data a.preliminary b.preliminary
        | Nonadd16 a, Nonadd16 b -> NonaddCS16.compare_data a.preliminary b.preliminary
        | Nonadd32 a, Nonadd32 b -> NonaddCS32.compare_data a.preliminary b.preliminary
        | AddGen a, AddGen b     -> AddCS.General.compare_data a.preliminary b.preliminary
        | AddVec a, AddVec b     -> AddCS.Vector.compare_data a.preliminary b.preliminary
        | Sank a, Sank b         -> SankCS.compare_data a.preliminary b.preliminary
        | FixedStates a, FixedStates b -> Fixed_states.compare_data a.preliminary b.preliminary
        | Dynamic a, Dynamic b   -> DynamicCS.compare_data a.preliminary b.preliminary
        | Kolmo a, Kolmo b       -> KolmoCS.compare_data a.preliminary b.preliminary
        | StaticMl a, StaticMl b ->
            IFDEF USE_LIKELIHOOD THEN
                MlStaticCS.compare_data a.preliminary b.preliminary
            ELSE
                failwith MlStaticCS.likelihood_error
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

(* This function assumes that nodea and nodeb come from a valid tree and
 * nodea is the parent of nodeb.  *)
let edge_distance clas nodea nodeb =
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
        | AddGen a, AddGen b ->
                (match clas with
                | `Static | `Any -> 
                        a.weight *. AddCS.General.distance a.final b.final
                | `Dynamic -> 0.)
        | AddVec a, AddVec b ->
                (match clas with
                | `Static | `Any -> 
                        a.weight *. AddCS.Vector.distance a.final b.final
                | `Dynamic -> 0.)
        | Sank a, Sank b ->
                (match clas with
                | `Static | `Any -> 
                        a.weight *. SankCS.distance a.final b.final
                | `Dynamic -> 0.)
        | FixedStates a, FixedStates b ->
                (match clas with
                | `Static | `Any -> 
                        a.weight *. Fixed_states.distance a.final b.final
                | `Dynamic -> 0.)
        | Dynamic a, Dynamic b ->
                (match clas with
                | `Dynamic | `Any -> 
                        (* Observe that we REQUIRE the single assignment for
                        * this collapse to be correct. *)
                        let d = 
                            DynamicCS.distance 0. a.preliminary b.preliminary 
                        in
                        a.weight *. d
                | `Static -> 0.)
        | Kolmo a, Kolmo b ->
              a.weight *. KolmoCS.tabu_distance a.final b.final
        | StaticMl a, StaticMl b ->
            IFDEF USE_LIKELIHOOD THEN
                let x , _(*ignore sumcost*) = cs_median 0 nodea nodeb None None None ch1 ch2 in
                match x with | StaticMl x -> 0.0 *. x.cost | _ -> assert false
            ELSE
                failwith MlStaticCS.likelihood_error
            END
        | Set a, Set b ->
              (match a.final.smethod with
               | `Strictly_Same ->
                     distance_lists a.final.set b.final.set 0.
               | `Any_Of _ ->
                     (* unf. we just take the full median and check the distance
                     *)
                     decr median_counter;
                     let m,_ (*ignore sumcost*) = cs_median !median_counter nodea nodeb None None None ch1 ch2 in
                     extract_cost m)
        | _ -> failwith "Incompatible characters (5)"
    and distance_lists (chs1:cs list) chs2 acc =
        match chs1, chs2 with
        | ch1 :: chs1, ch2 :: chs2 ->
              distance_lists chs1 chs2 (acc +. distance_two ch1 ch2)
        | [], [] -> acc
        | _ -> failwith "Incompatible characters (6)" in
    distance_lists nodea.characters nodeb.characters 0.

let all_types = 
    [ `Add ; `Annchrom ; `Breakinv ; `Chrom ; `Genome ; `Kolmo ; `Nonadd ;
      `Sank; `Seq ; `StaticMl; `Ml; `FixedStates ]

let type_string = function
    | `Add -> "additive"
    | `Annchrom -> "annotated chromosomes"
    | `Breakinv -> "break inversion"
    | `Chrom -> "chromosome"
    | `Genome -> "genome"
    | `Kolmo -> "kolmo"
    | `Nonadd -> "nonadditive"
    | `Sank -> "sankoff"
    | `FixedStates -> "fixed states"
    | `Seq -> "sequence"
    | `StaticMl -> "static likelihood"
    | `Ml -> "dynamic likelihood"

let has_to_single : to_single list = [`Seq ; `Chrom; `Annchrom; `Breakinv;]

module ToSingleModule = Set.Make (struct type t = to_single let compare = compare end)

let all_characters_single =
    List.fold_left (fun acc x -> ToSingleModule.add x acc) 
                   (ToSingleModule.empty)
                   (all_types)

let not_to_single =
    ToSingleModule.elements 
        (List.fold_left (fun acc x -> ToSingleModule.remove x acc)
                        (all_characters_single) 
                        (has_to_single) )

(**[distance_of_type] return the distance of two node : nodea and nodeb. it's
* being called by [chekc_cost] of allDirChar.ml *)
let distance_of_type ?branches ?(para=None) ?(parb=None) t missing_distance
    ({characters=chs1} as nodea) ({characters=chs2} as nodeb) =
    if debug_distance then
        Printf.printf "\n Node.distance_of_type on node#.%d and node#.%d -> %!" nodea.taxon_code nodeb.taxon_code;
    let has_t x = List.exists (fun z -> z = x) t
    and filter_dynamic res x = match x with
        | `Seq    | `Breakinv
        | `Chrom  | `Annchrom
        | `Genome | `Ml as x -> x :: res
        | _ -> res
    in
    let has_nonadd = has_t `Nonadd
    and has_add = has_t `Add
    and has_sank = has_t `Sank
    and has_fixedstate = has_t `FixedStates 
    and dy_t = List.fold_left filter_dynamic [] t 
    and has_kolmo = has_t `Kolmo in
    let rec distance_two ch1 ch2 (bl: float option) : float =
        match ch1, ch2 with
        | Nonadd8 a, Nonadd8 b when has_nonadd ->
            a.weight *. NonaddCS8.distance a.final b.final
        | Nonadd16 a, Nonadd16 b when has_nonadd ->
            a.weight *. NonaddCS16.distance a.final b.final
        | Nonadd32 a, Nonadd32 b when has_nonadd ->
            a.weight *. NonaddCS32.distance a.final b.final
        | AddVec a, AddVec b when has_add ->
            a.weight *. AddCS.Vector.distance a.final b.final
        | AddGen a, AddGen b when has_add ->
            a.weight *. AddCS.General.distance a.final b.final
        | Sank a, Sank b when has_sank ->
            a.weight *. SankCS.distance a.final b.final
        | FixedStates a, FixedStates b when has_fixedstate ->
            a.weight *. Fixed_states.distance a.final b.final
        | Dynamic a, Dynamic b ->
            a.weight *. DynamicCS.distance_of_type dy_t missing_distance a.final b.final bl
        | Kolmo a, Kolmo b when has_kolmo ->
            a.weight *. KolmoCS.distance a.final b.final
        | Set a, Set b ->
            begin match a.final.smethod with
                | `Strictly_Same ->
                    let bls = Array.to_list (Array.create (List.length b.final.set) bl) in
                    distance_lists a.final.set b.final.set bls 0.
               | `Any_Of _ ->
                     (* unf. we just take the full median and check the distance *)
                     decr median_counter;
                     let m, _ (*ignore sumcost*) = cs_median !median_counter nodea nodeb None None None ch1 ch2 in 
                     extract_cost m
            end
        | StaticMl a, StaticMl b ->
            IFDEF USE_LIKELIHOOD THEN
                let x, _ (*ignore sumcost*) = cs_median 0 nodea nodeb None bl (Some 0.0) ch1 ch2 in
                match x with | StaticMl x -> a.weight *. (x.cost -. (a.cost +. b.cost))
                             | _ -> assert false
            ELSE
                failwith MlStaticCS.likelihood_error
            END
        | _ -> 0.0
    and distance_lists chs1 chs2 branchs acc =
        match chs1, chs2, branchs with
        | ch1 :: chs1, ch2 :: chs2, b :: bs ->
              distance_lists chs1 chs2 bs (acc +. distance_two ch1 ch2 b)
        | [], [], [] -> acc
        | _ -> failwith "Incompatible characters (6)"
    in
    distance_lists chs1 chs2 (convert_2_lst nodea branches) 0.


let distance ?(para=None) ?(parb=None) missing_distance
    ({characters=chs1} as nodea) ({characters=chs2} as nodeb) =
        if debug_distance then 
            Printf.printf "\n Node.distance on node#.%d and node#.%d -> %!" nodea.taxon_code nodeb.taxon_code;
    let rec distance_two ch1 ch2 =
        match ch1, ch2 with
        | Nonadd8 a, Nonadd8 b ->
              a.weight *. NonaddCS8.distance a.final b.final
        | Nonadd16 a, Nonadd16 b ->
              a.weight *. NonaddCS16.distance a.final b.final
        | Nonadd32 a, Nonadd32 b ->
              a.weight *. NonaddCS32.distance a.final b.final
        | AddVec a, AddVec b ->
              a.weight *. AddCS.Vector.distance a.final b.final
        | AddGen a, AddGen b ->
              a.weight *. AddCS.General.distance a.final b.final
        | Sank a, Sank b ->
              a.weight *. SankCS.distance a.final b.final
        | FixedStates a, FixedStates b ->
              a.weight *. Fixed_states.distance a.final b.final
        | Dynamic a, Dynamic b ->
              a.weight *. DynamicCS.distance missing_distance a.final b.final
        | Kolmo a, Kolmo b ->
              a.weight *. KolmoCS.distance a.final b.final
        | StaticMl a, StaticMl b ->
            IFDEF USE_LIKELIHOOD THEN
                let x, _ (*ignore sumcost*) = cs_median 0 nodea nodeb None None None ch1 ch2 in
                match x with | StaticMl x -> a.weight *. (x.cost -. (a.cost +. b.cost))
                             | _ -> assert false
            ELSE
                failwith MlStaticCS.likelihood_error
            END
        | Set a, Set b ->
              (match a.final.smethod with
               | `Strictly_Same ->
                     distance_lists a.final.set b.final.set 0.
               | `Any_Of _ -> (* TODO:: check this is correct *)
                     (* unf. we just take the full median and check the distance *)
                     decr median_counter;
                     let m, _ (*ignore sumcost*) = cs_median !median_counter nodea nodeb None None None ch1 ch2 in
                     extract_cost m)
        | _ -> failwith "Incompatible characters (5)"
    and distance_lists chs1 chs2 acc =
        match chs1, chs2 with
        | ch1 :: chs1, ch2 :: chs2 ->
              distance_lists chs1 chs2 (acc +. distance_two ch1 ch2)
        | [], [] -> acc
        | _ -> failwith "Incompatible characters (6)"
    in
    let res = distance_lists chs1 chs2 0. in
    if debug_distance then Printf.printf "resdis=%f\n%!" res;
    res

(* Calculates the cost of joining the node [n] between [a] and [b] in a tree *)
(* [a] must be the parent (ancestor) of [b] *)
let dist_2 minimum_delta n a b =
    let debug = false in
    if debug then Printf.printf "node.dist_2, join node#.%d and node#.%d together\n%!"
    a.taxon_code b.taxon_code;
    let rec ch_dist delta_left n' a' b' =
        match n', a', b' with
        | StaticMl nn, StaticMl aa, StaticMl bb ->
               (**   a
                * c_/|    ~I am using the cs_median function to calculate
                *   \|     exactly. This shouldn't be done, but this function
                *    x     does not get called at all under likelihood, at
                *   / \    least to the best of my knowledge, since a grep only
                *  n   b   shows results in chartree.ml --one direction tree.
               **)
            IFDEF USE_LIKELIHOOD THEN
                let x', _ (*ignore sumcost*)= cs_median 0 n b None None None n' b' in
                (* the min_code is the only thing used so this is sufficient *)
                let x = if n.min_child_code < b.min_child_code then n else b in
                let c,_ = cs_median (-1) x a None None None a' x' in
                match c with | StaticMl c -> nn.weight *. c.cost | _ -> assert false
            ELSE
                failwith MlStaticCS.likelihood_error
            END
        | Nonadd8 n, Nonadd8 a, Nonadd8 b ->
                let dist = NonaddCS8.dist_2 n.final a.final b.final in
                let res = n.weight *. dist in
                if debug then Printf.printf "Nonadd8, res=%f*%f=%f; %!" dist n.weight res;
                res
        | Nonadd16 n, Nonadd16 a, Nonadd16 b ->
              n.weight *. NonaddCS16.dist_2 n.final a.final b.final
        | Nonadd32 n, Nonadd32 a, Nonadd32 b ->
              n.weight *. NonaddCS32.dist_2 n.final a.final b.final
        | AddVec n, AddVec a, AddVec b ->
              n.weight *. AddCS.Vector.distance_2 n.final a.final b.final
        | AddGen n, AddGen a, AddGen b ->
                let dist = AddCS.General.distance_2 n.final a.final b.final in
                let res = n.weight *. dist in
                if debug then Printf.printf "AddGen, res=%f*%f=%f; %!" dist n.weight res;
                res
        | Sank n, Sank a, Sank b ->
              n.weight *. SankCS.dist_2 n.final a.final b.final
        | FixedStates n, FixedStates a, FixedStates b ->
              n.weight *. Fixed_states.dist_2 n.final a.final b.final
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
        | AddGen _, _, _ | Sank _, _, _ | FixedStates _, _, _ | Dynamic _, _, _ | AddVec _, _, _
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
                let add_dist = ch_dist max_delta n a b in
                  chars (acc +. add_dist) (ncs, acs, bcs)
            | [], [], [] -> acc
            | _ -> raise (Illegal_argument "dist_2_chars")
    in
    let excludes = excludes_median n a in
    if has_excluded excludes
    then infinity
    else 
        let res = chars 0. (n.characters, a.characters, b.characters) in
        if debug then Printf.printf "return join cost = %f\n%!" res;
        res

let extract_stat = function
    | (Data.Stat (a, b), _) -> (b, a)
    | _ -> raise (Illegal_argument "extract_stat")

(*
    let tcharacters = Hashtbl.find !data.Data.taxon_characters tcode in
            let chfilenames = !data.Data.character_codes in
    Hashtbl.find tcharacters code,Hashtbl.find chfilenames code
*)
let extract_fixedstates data fs tcode = (*tcode = taxon code*) 
    match fs with
    | Data.FS code, _  -> (*code is charactor code belong to the taxon*)
            let specs = Hashtbl.find data.Data.character_specs code in  
            let specs, weight =
                match specs with 
                | Data.Static (Data.FixedStates fs_spec) ->  fs_spec, Data.get_weight_from_fs_spec fs_spec
                | _ -> failwith "extract fiexed states"
            in 
            (* val off_array : Data.fixed_state_spec -> int -> t_w_seqtbl *)
            let tws =  Fixed_states.off_array specs tcode in 
            {
                preliminary = tws;
                final = tws; 
                cost = 0.; 
                weight = weight;
                sum_cost = 0.;
                time = None,None,None;
            }
    | Data.Dyna code, _ -> 
                        failwith("we have cs=dyna instead of fs")
    | _ -> raise (Illegal_argument "extract fixedstates")

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
            time = None,None,None;
          } 
    | Data.Stat (code, _), _ ->
          raise (Illegal_argument ("Stat" ^ (string_of_int code))) 
    | Data.FS code, _  ->
          raise (Illegal_argument ("FS " ^ (string_of_int code))) 

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
              let seq = Array.map (fun x -> x.Data.seq) chrom_data.Data.seq_arr
              in
              KolmoCS.of_array kspec [|(seq, chcode)|] chcode num_taxa
          in 
          { preliminary = dyna; 
            final = dyna; 
            cost = 0.; 
            weight = weight;
            sum_cost = 0.;
            time = None,None,None;
          } 
    | Data.Stat (code, _), _ ->
          raise (Illegal_argument ("Stat" ^ (string_of_int code))) 
    | Data.FS code, _  ->
          raise (Illegal_argument ("FS " ^ (string_of_int code)))
              

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
    let process_all_lists lists = (* list of each column *)
        current_snapshot "This is before lists";
        let lists = 
            (* fold over each column into SetList, increasing weight if it is
             * allready a member *)
            List.fold_left (fun acc x ->
                let (code, weight, lst) = characters x in
                let lst = 
                    (* transform each column to bitset *)
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

(* Classify into discreet sections returning a map of the last character with
 * the identifying property (exact column in this case), and a float representing
 * the weight imposed on that character.
     * size  --of alphabet
     * chars --list of column codes to compress
     * data  --Data.d *)
let classify size chars data =
    let all_static = 
        Hashtbl.fold
            (fun code spec acc -> match spec with
                | Data.Static (Data.NexusFile nf) ->
                    (* categorize likelihood and unordered characters *)
                    begin match nf.Nexus.File.st_type with
                        | Nexus.File.STUnordered
                        | Nexus.File.STLikelihood _ -> (code, nf) :: acc
                        | _ -> acc
                    end
                | Data.Static (Data.FixedStates fs) -> acc (*(code, spec) :: acc*)
                | _ -> acc)
            data.Data.character_specs
            []
    in
    (* transform each static column into a list with current weight. ~NL *)
    let taxa (code,spec) = 
        let weight,observed = match spec.Nexus.File.st_type with
            | Nexus.File.STUnordered ->
                spec.Nexus.File.st_weight, spec.Nexus.File.st_observed
            | Nexus.File.STLikelihood m ->
                spec.Nexus.File.st_weight, spec.Nexus.File.st_observed
            | _ -> assert false
        in
        Hashtbl.fold 
            (fun _ taxon_chars acc ->
                let lst =
                    try match Hashtbl.find taxon_chars code with
                        | (Data.Stat (c, (Some v)), `Specified) -> (c, weight, v)
                        | (Data.Stat (c, v), _) -> (c, weight, `List observed)
                        | _ -> failwith "Impossible 2?"
                    with | Not_found -> (code, weight, `List observed)
                in
                lst :: acc)
            data.Data.taxon_characters
            []
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
            | [] -> failwith "Nothing?" (* must be of least length one, *)
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

let generate_taxon do_classify laddgencode laddveccode lnadd8code lnadd16code
                   lnadd32code lnadd33code lsankcode dynamics fixedstates kolmogorov
                   static_ml data cost_mode =
        let calc_total treesnum directions nodenum = 
            let res = treesnum * ( directions *(nodenum-1) + nodenum )
            in
            res
        in
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
            (* get weight of character; map has Data.weight included *)
            let get_weight c = match weights with
                | None   -> Data.get_weight c !data
                | Some v -> All_sets.IntegerMap.find c v
            in
            let table = Hashtbl.create 1667 in
            let weights =
                All_sets.Integers.fold
                    (fun x acc ->
                        try let w = get_weight x in
                            if Hashtbl.mem table w then
                                let lst = Hashtbl.find table w in
                                Hashtbl.replace table w (x :: lst)
                            else
                                Hashtbl.add table w [x];
                            All_sets.Floats.add w acc
                        with | Not_found -> acc)
                    codes All_sets.Floats.empty
            in
            let res = All_sets.Floats.fold
                (fun x acc ->
                    let lst = Hashtbl.find table x in
                    (x, lst) :: acc)
                weights []
            in
            List.fold_left
                (fun acc (w, lst) ->
                    (character_code_gen (),
                    (List.map (fun x -> w, x) lst)) :: acc)
                [] res
        in
        let group_ml_by_model lst =
            let get_function code =
                match Hashtbl.find (!data).Data.character_specs code with
                | Data.Static (Data.NexusFile spec) ->
                    begin match spec.Nexus.File.st_type with
                        | Nexus.File.STLikelihood x -> x.MlModel.spec
                        | _ -> assert false
                    end
                | _ -> assert false
            in
            MlModel.categorize_by_model get_function lst
        and group_by_sets lst =
            let curr = Hashtbl.create 1667 in
            let sets = List.fold_left (* list of all the set names *)
                (fun acc code ->
                    try let name = Hashtbl.find (!data).Data.character_codes code in
                        let set = Hashtbl.find (!data).Data.character_nsets name in
                        Hashtbl.add curr set code;
                        if List.mem set acc then
                            acc
                        else
                            set::acc
                    with | Not_found ->
                        Hashtbl.add curr "" code; acc)
                [""] lst
            in
            List.map (Hashtbl.find_all curr) sets
        in
        let nadd8weights = classify 8 do_classify lnadd8code !data
        and nadd16weights = classify 16 do_classify lnadd16code !data
        and nadd32weights = classify 32 do_classify lnadd32code !data in
        let laddveccode = group_in_weights None laddveccode
        and laddgencode = group_in_weights None laddgencode
        and lnadd8code = group_in_weights nadd8weights lnadd8code
        and lnadd16code = group_in_weights nadd16weights lnadd16code
        and lnadd32code = group_in_weights nadd32weights lnadd32code
        and lstaticmlcode = 
            let () = 
                IFDEF USE_LIKELIHOOD THEN
                    (* set garbage collector frequency; multiplier = 4 *)
                    MlStaticCS.gc_alloc_max 
                        (calc_total 4 3 (!data).Data.number_of_taxa)
                ELSE
                    ()
                END
            in
            (* apply classify on a list of integers to classify *)
            let set_of_list =
                List.fold_left (fun a v -> All_sets.Integers.add v a)
                                All_sets.Integers.empty
            in
            (* create map of weight classes *)
            let lk_classify_weights = function
                | (x::xs) as all ->
                    (* already characterized by group and model, so = alph *)
                    let alph_len =
                         let x = match Hashtbl.find (!data).Data.character_specs x with
                            | Data.Static (Data.NexusFile x) -> x 
                            | _ -> assert false
                        in
                        match x.Nexus.File.st_type with
                        | Nexus.File.STLikelihood s -> MlModel.get_alphabet s
                        | _ -> assert false
                    in
                    classify alph_len (MlStaticCS.compress && do_classify) all !data
                | [] -> assert false
            (* convert characters and group them by weight *)
            and lk_group_weights weights chars = 
                chars --> set_of_list --> group_in_weights weights
            in
            (* group sets, model then compress columns *)
            static_ml
                --> group_by_sets
                --> List.map group_ml_by_model
                --> List.flatten
                --> List.map (fun x -> lk_group_weights (lk_classify_weights x) x)
        and lsankcode =
            let () = 
                 SankCS.set_gc_alloc_max 
                        (calc_total 4 1 (!data).Data.number_of_taxa)
            in
            List.map (fun x -> cg (), x) lsankcode in
        let add_codes ((_, x) as y) = 
            y, Array.map snd (Array.of_list (List.rev x)) in
        let laddveccode = List.map add_codes laddveccode 
        and laddgencode = List.map add_codes laddgencode 
        and lnadd8code = List.map add_codes lnadd8code
        and lnadd16code = List.map add_codes lnadd16code
        and lnadd32code = List.map add_codes lnadd32code in
        (* We need ways of making empty characters when a character is
           unspecified *)
        let get_static_encoding code =
            let specs = Hashtbl.find !data.Data.character_specs code in
            match specs with
            | Data.Static (Data.NexusFile encoding) -> encoding
            | Data.Static _ 
            | Data.Dynamic _
            | Data.Kolmogorov _
            | Data.Set -> failwith "get_static_encoding"
        in
        let module Enc = Parser.OldHennig.Encoding in
        let gen_add code =
            let enc = get_static_encoding code in
            (Data.Stat (code, Some (`List enc.Nexus.File.st_observed)), `Unknown) 
        in
        let gen_nadd code =
            let enc = get_static_encoding code in
            (Data.Stat (code, Some (`List enc.Nexus.File.st_observed)), `Unknown) 
        in
        let gen_dynamic code =
            let alph = Data.get_alphabet !data code in
            (* print_endline ("adding sequence with code " ^ string_of_int code); *)
            let empty_seq = Data.get_empty_seq alph in
            let chrom_data = Data.set_dyna_data [|empty_seq|] in 
            (Data.Dyna (code, chrom_data), `Unknown)
        in
        let gen_fixedstates code = 
            (*looks like we just need the code, then later we can get
            * fixedstates's spec from data.Data.character_specs with that code.
            * everything -- including the original sequence is in there.*)
            (*to do : what do we need here?*)
            (*let alph = Data.get_alphabet !data code in
            let empty_seq = Data.get_empty_seq alph in
            let specs = Hashtbl.find !data.Data.character_specs code in
            let states = 
                match specs with
                | Data.Static x -> 
                        (match x with  (*to do: what should we put here?*)
                        | Data.FixedStates enc -> []
                        | _ -> failwith "gen_sank is not for fixedstates" )
                | _ -> assert false 
            in
            let fs_data = Data.set_fs_data [|empty_seq|] [] in*)
            (Data.FS code, `Unknown)
        in
        let gen_sank code =
            (* print_endline ("adding sankoff with code " ^ string_of_int code);
            * *)
            let specs = Hashtbl.find !data.Data.character_specs code in
            let states = 
                match specs with
                | Data.Static x -> 
                        (match x with 
                        | Data.NexusFile enc -> `List enc.Nexus.File.st_observed
                        | _ -> failwith "gen_sank is not for fixedstates" )
                | _ -> assert false 
            in
            (Data.Stat (code, Some states), `Unknown)
        in
        current_snapshot "Done";
        !data, 
        fun tcode acc ->
            current_snapshot "Generating taxon";
            let debug = false in
            let tcharacters = Hashtbl.find !data.Data.taxon_characters tcode in
            let chfilenames = !data.Data.character_codes in
            if debug then Printf.printf "\n Generating Taxon %d Has Characters: \n%!" tcode;
(*            Hashtbl.iter (fun k _ -> Printf.printf "%d, " k) tcharacters;*)
(*            print_newline ();*)
            let get_character_with_code_n_weight gen_new (w, acc, cnt) (weight, code) = 
                try 
                    weight, (Hashtbl.find tcharacters code) :: acc, cnt + 1
                with
                | Not_found -> weight, (gen_new code) :: acc, cnt
            in
            let get_character_with_code gen_new acc code = 
                if debug then begin
                    Printf.printf "get char with code=%d\n%!" code;
                end;
                try 
                    (Hashtbl.find tcharacters code,Hashtbl.find chfilenames code) :: acc
                with
                | Not_found ->
                        if debug then Printf.printf "not found, gen new;\n%!";
                        (gen_new code,"") :: acc
            in
            let addmapper gen_new ((x, y), arr) =
                let a, b, cnt = 
                    List.fold_left (get_character_with_code_n_weight gen_new) 
                    (1., [], 0) y
                in
                x, (a, b), arr
            in
            let ladd_gen_chars= List.map (addmapper gen_add)  laddgencode
            and ladd_vec_chars= List.map (addmapper gen_add)  laddveccode
            and lnadd8_chars  = List.map (addmapper gen_nadd) lnadd8code
            and lnadd16_chars = List.map (addmapper gen_nadd) lnadd16code
            and lnadd32_chars = List.map (addmapper gen_nadd) lnadd32code
            and lnadd33_chars = []
            and ldynamic_chars = 
                List.fold_left (get_character_with_code gen_dynamic) [] dynamics
            and lfixedstates_chars =
                List.fold_left (get_character_with_code gen_fixedstates) [] fixedstates
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
                    w; time = None,None,None }
            in
            let result = (* NONADD8 *)
                List.fold_left 
                (fun acc (a, b, c) -> 
                    add_characters (NonaddCS8.of_parser !data c)
                    (fun c w -> Nonadd8 (make_with_w c w)) acc (a, b))
                result lnadd8_chars
            in
            let result =(* NONADD16 *)
                List.fold_left 
                (fun acc (a, b, c) -> 
                    add_characters (NonaddCS16.of_parser !data c)
                    (fun c w -> Nonadd16 (make_with_w c w)) acc (a, b))
                result lnadd16_chars
            in
            let result = (* NONADD32 *)
                List.fold_left 
                (fun acc (a, b, c) -> add_characters (NonaddCS32.of_parser !data c)
                (fun c w -> Nonadd32 (make_with_w c w)) acc (a, b))
                result lnadd32_chars
            in
            let result = 
                match lnadd33_chars with
                | _ -> result
            in
            let result =  (* ADDITIVE General *)
                List.fold_left 
                (fun acc (a, b, _) -> add_characters (AddCS.General.of_parser !data)
                (fun c w -> AddGen (make_with_w c w)) acc (a, b))
                result ladd_gen_chars
            in
            let result =  (* ADDITIVE Vectorized *)
                List.fold_left 
                (fun acc (a, b, _) -> add_characters (AddCS.Vector.of_parser !data)
                (fun c w -> AddVec (make_with_w c w)) acc (a, b))
                result ladd_vec_chars
            in
            let result = (* DYNAMIC *)
                match ldynamic_chars with
                | [] -> result
                | _ ->
                      let c = 
                          List.map 
                          (fun (dyna,fname) -> extract_dynamic !data dyna tcode) 
                          ldynamic_chars 
                      in
                      let c : cs list = List.map (fun c -> Dynamic c) c in
                      { result with characters = c @ result.characters } 
            in
            let result = (* FIXED STATES *)
                match lfixedstates_chars with 
                | [] -> result 
                | _ ->
                      let c = 
                          List.map 
                          (fun (fs,fname) -> extract_fixedstates !data fs tcode) 
                          lfixedstates_chars 
                      in
                      let c : cs list = List.map (fun c -> FixedStates c) c in
                      { result with characters = c @ result.characters }
            in
            let result = (* SANK *)
                let single_lsank_chars_process result (code, lst) =
                    match lst with
                    | [] -> result
                    | _ ->
                            let v = List.map (fun (sankcs,fname) -> extract_stat
                            sankcs ) lst in
                            let tcm =
                                match v with
                                | (_, code) :: _ -> Data.get_tcm code !data
                                | _ -> failwith "This is impossible"
                            in
                            let arr= Array.of_list v in
                            let c, _ = SankCS.of_parser tcm (arr, tcode) code in
                            let c = Sank { preliminary = c; final = c; cost = 0.;
                                          sum_cost = 0.;
                            weight = 1.; time = None,None,None; } in
                            { result with characters = c :: result.characters }
                in
                List.fold_left single_lsank_chars_process result lsank_chars
            in
            let result = (* KOLMO *)
                match lkolmo_chars with
                | [] -> result
                | _ ->
                        let c = 
                            List.map 
                            (fun (kolm,fname) -> extract_kolmo !data kolm tcode)
                            lkolmo_chars
                        in
                        let total_cost = 
                            List.fold_left (fun acc c -> acc +.  c.sum_cost)
                            0. c 
                        in
                        let c = List.map (fun c -> Kolmo c) c in
                        { result with characters = c @ result.characters; 
                        total_cost = total_cost +. result.total_cost;
                        node_cost = total_cost +. result.node_cost;}
            in
            let result = (* ML *)
                let single_ml_group =
                  IFDEF USE_LIKELIHOOD THEN
                    let seperate_data dat =
                        let pairs = dat --> List.map snd --> List.flatten in
                        let fsts = List.map fst pairs and snds = List.map snd pairs in
                        fsts,snds
                    in
                    fun result -> function
                        | [] -> result
                        | all_data ->
                            let ws,cs = seperate_data all_data in
                            let spec = 
                                match Hashtbl.find (!data).Data.character_specs
                                                    (List.hd cs) with
                                | Data.Static (Data.NexusFile x) -> x 
                                | _ -> assert false
                            in
                            let c = 
                                cs --> List.map 
                                        (fun x -> 
                                            try Hashtbl.find tcharacters x
                                            with | Not_found -> 
                                                failwithf "Couldn't find char %d" x)
                                   --> List.map extract_stat (* bitset * code *)
                                   --> Array.of_list
                                   --> MlStaticCS.of_parser spec (Array.of_list ws)
                            in
                            let cost = MlStaticCS.root_cost c in
                            let c = 
                                StaticMl { preliminary = c;
                                           final = c;
                                           cost = cost;
                                           sum_cost = cost;
                                           weight = 1.;
                                           time = None,None,None; }
                            in
                            { result with characters = c :: result.characters;
                                          total_cost = result.total_cost +. cost; }
                  ELSE
                    fun result _ -> result
                  END
                in
                List.fold_left single_ml_group result lstaticmlcode
            in
            let () = current_snapshot "Finished taxon" in
            result :: acc

let node_contents_compare a b = match a, b with
    | Nonadd8 _ , Nonadd8 _
    | Nonadd16 _, Nonadd16 _
    | Nonadd32 _, Nonadd32 _
    | AddGen _  , AddGen _
    | AddVec _  , AddVec _
    | Sank _    , Sank _
    | FixedStates _ , FixedStates _ 
    | Dynamic _ , Dynamic _
    | Kolmo _   , Kolmo _
    | StaticMl _, StaticMl _
    | Set _     , Set _      ->  0 (* Sets correct? *)
    | Nonadd8 _ , _          -> -1
    | _         , Nonadd8 _  ->  1
    | Nonadd16 _, _          -> -1
    | _         , Nonadd16 _ ->  1
    | Nonadd32 _, _          -> -1
    | _         , Nonadd32 _ ->  1
    | AddGen _  , _          -> -1
    | _         , AddGen _   ->  1
    | AddVec _  , _          -> -1
    | _         , AddVec _   ->  1
    | Sank _    , _          -> -1
    | _         , Sank _     ->  1
    | FixedStates _ , _      -> -1
    | _         , FixedStates _ -> 1
    | Dynamic _ , _          -> -1
    | _         , Dynamic _  ->  1
    | Kolmo _   , _          -> -1
    | _         , Kolmo _    ->  1
    | StaticMl _, _          -> -1
    | _         , StaticMl _ ->  1


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
              weight = 1.; time = None,None,None } 
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
              weight = 1.; time = None,None,None } 
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
    let cost_mode = match nodes with
        | h :: tl ->
            List.iter (fun x -> assert(h.cost_mode = x.cost_mode)) tl;
            h.cost_mode
        | [] -> `NotLikelihood
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

let is_available in_data =
    let res = 
        match (List.hd in_data) with
        |  Dynamic x_dycs ->                
                DynamicCS.is_available x_dycs.preliminary
        | _ -> 0
    in
    res

let flatten cs_lst=
    let dycs_lst = List.map (fun x -> match x with
    | Dynamic x_dycs -> x_dycs.preliminary
    | _ -> failwith ("we only deal with DynamicCS now")
    ) cs_lst
    in
    DynamicCS.flatten dycs_lst 

let flatten_cslist (characters_lst: cs list list) data =
    let debug = false in
    if debug then Printf.printf "node.flatten_cslst -> %!";
    let (seq_lstlstlst:(int * Sequence.s list list ) list)
    = mapN flatten characters_lst in 
    (* for BreakinvCS, seq_lstlstlst is like following ( other dynamicCS data
    * types I know in are simpler than this):
      {
          (
            [seq of 1th median of 1th chromosome of 1th node;
             seq of 2th median of 1th chromosome of 1th node;
             .... ];
            [seq of 1th median of 1th chromosome of 2th node;
             ...... 2th ..........1th ......................;
             .....];
            ......
          )
          (
            [seq of 1th median of 2th chromosome of 1th node;
             seq of 2th median of 2th chromosome of 1th node;
             .... ];
            [seq of 1th median of 2th chromosome of 2th node;
             ...... 2th ..........2th ......................;
             .....];
            ......
          )
          ......
          ......
          (
              [seq of 1th median of nth chromosome of 1th node;
               seq of 2th ....................................;
               ..... ];
              [seq of 1th median of nth chromosome of 1th node;
               seq of 2th median of ..........................
               ..... ];
               .......
          )
      }
    *)
    if debug then Printf.printf "before flatten_cslist seq_lstlstlst = \n%!";
    let getfilename filename_w_taxoncode =
        (*if we have more than one chromsome in a taxon of a file
        * fiel1.fas
        * >t1
        * a b c @ d e
        * >t2
        * a c b @ e d
        * each chromosome become a character, with name filename^":"^"which
        * taxon is it in"
        * ( see data.ml
        *   let locus_name = 
            let c = ref (-1) in
            ref (fun () -> incr c; file ^ ":" ^ string_of_int !c)
            in)
        * for the example above, we have
        *         chmosome1 chmosome2
        * file1   file1:0   file1:1
        *
        * this function get the filename out of string "filename:int"
        * *)
        let idx = ref 0 and len = ref 0 in
        String.iter (fun chr -> 
            if chr=':' then len := !idx
            else idx := !idx + 1;
        ) filename_w_taxoncode;
        if !len=0 then (*just one chromosome in that taxon of that file, there is no ":"*)
            filename_w_taxoncode
        else 
            String.sub filename_w_taxoncode 0 !len
    in
    let tbl = Hashtbl.create 50 in (*filename -> chrom_node_med lstlstlst*)
    List.iter ( fun (character_code,thischrom) ->
        let filename_w_taxoncode = 
            try Hashtbl.find data.Data.character_codes character_code 
            with Not_found -> failwith "could not find some charater in data.charcter_codes" 
        in
        let filename = getfilename filename_w_taxoncode in
        if filename="" then failwith "node.ml getfilename return empty string";
        if debug then Printf.printf "add chrom_node_med to tbl with filename = %s\n%!" filename;
        if Hashtbl.mem tbl filename then
            let oldrecord = Hashtbl.find tbl filename in
            Hashtbl.replace tbl filename (oldrecord@[thischrom])
        else
            Hashtbl.add tbl filename [thischrom]
    ) seq_lstlstlst;
    let file_chrom_node_med = Hashtbl.fold (fun filename chrom_node_med acc ->
        acc@[(filename,chrom_node_med)]
    ) tbl [] in
    (*[tranpose take a matrix, return its transposed form]
     * { 
         * a1,a2,a3,....,an;
         * b1,b2,b3,....,bn;
         * ......
         * x1,x2,x3,....,xn;
         * ....
         * }
         become
        { a1,b1,....,x1,...;
          a2,b2,....,x2,...;
          .....
          an,bn,....,xn,...;}*)
    let transpose in_lstlst =
        let len = List.length (List.hd in_lstlst) in
        let out_lstlst = ref [] in
        for i = 0 to len -1 do
            out_lstlst := (ref [])::(!out_lstlst)
        done;
        let out_lstlst = !out_lstlst in
        List.iter (fun in_lst ->
            let len = List.length in_lst in
            for i = 0 to len-1 do
                let out_i = List.nth out_lstlst i in
                out_i := (List.nth in_lst i)::(!out_i);  (*reversed?*)
            done;
        ) in_lstlst;
        let out_lstlst = List.map (fun out_i -> !out_i) out_lstlst in
        out_lstlst
    in
    let generate_delimiters (seqlstlst: Sequence.s list list) = 
        List.map (fun seqlst ->
            List.map (fun seq -> Sequence.length seq) seqlst 
        ) seqlstlst 
    in
    let file_node_chrom_med = 
        List.map (fun (filename,chrom_node_med) -> filename,transpose chrom_node_med
        ) file_chrom_node_med in
    let file_node_median_chrom = 
        List.map ( fun (filename,node_chrom_med) ->
            List.map (fun chrom_med -> transpose chrom_med
            ) node_chrom_med 
    ) file_node_chrom_med in 
    let node_file_median_chrom_seq = transpose file_node_median_chrom in 
    let node_file_median_chrom_seqdeli = 
        List.map (fun file_median_chrom ->
            List.map generate_delimiters file_median_chrom
    ) node_file_median_chrom_seq
    in
    node_file_median_chrom_seq, node_file_median_chrom_seqdeli

let single_to_multi_chromosome single_cs =
    let dynCS_t_list:DynamicCS.t list =
        match single_cs with
            | Dynamic cs ->
                    DynamicCS.single_to_multi cs.preliminary
            | _ -> 
                failwith ("single-chromosome to multi-chromosome :\ 
                    we only deal with DynamicCS now")
    in
    let dynCS_t_r_list = 
        List.map (fun dynCS_t -> 
            match single_cs with 
            |Dynamic s_cs -> {s_cs with preliminary = dynCS_t; final = dynCS_t}
            | _ -> failwith ("single to multichromosome : we only deal with DynamicCS now")
            )dynCS_t_list 
    in
    List.map (fun dynCS_t_r -> Dynamic dynCS_t_r ) dynCS_t_r_list 

let multi_to_single_chromosome node_data file_median_seq file_median_chrom_seqdeli (*newseq delimiters*) =
    let old_characters = node_data.characters in
    let new_characters = 
    match (List.hd old_characters) with
    | Dynamic cs ->
            let new_preliminary_lst = 
                DynamicCS.update_t cs.preliminary file_median_seq file_median_chrom_seqdeli in
            List.map (fun x -> 
                Dynamic {cs with preliminary = x; final=x;}
            ) new_preliminary_lst
    | _ -> 
    failwith ("multichromosome to singlechromosome : we only deal with DynamicCS now")
    in
    { node_data with characters = new_characters }   

let transform_multi_chromosome ( nodes : node_data list ) data =
    let available = 
        is_available (List.hd nodes).characters 
    in
    if (available=1) then begin
    let characters_lst: cs list list = List.map (fun x -> x.characters) nodes in
    let (node_file_median_chrom_seq:Sequence.s list list list list),
    (node_file_median_chrom_seqdeli: int list list list list) = 
        flatten_cslist characters_lst data
    in
        (* now the seq_lstlstlst is like this:
          file 1:
            {
            (
            [seq of 1th median of 1th chromosome of 1th node;
             seq of 1th median of 2th chromosome of 1th node;
             .... ];
            [seq of 2th median of 1th chromosome of 1th node;
             ...... 2th ..........2th ......................;
             .....];
            ......
            )
            (
            [seq of 1th median of 1th chromosome of 2th node;
             seq of 1th median of 2th chromosome of 2th node;
             .... ];
            [seq of 2th median of 1th chromosome of 2th node;
             ...... 2th ..........2th ......................;
             .....];
            ......
            )
            .......
            ( ... )
          };
          file 2:
            {
                ......
            }
        * *)
        (*debug msg 
        Printf.printf "check delimiters list after flatten =>\n%!";
        List.iter (fun intlstlstlst ->
        Printf.printf "node{\n%!";
        List.iter(fun intlstlst -> 
            Printf.printf "file(%!";
            List.iter (fun intlst -> 
                Printf.printf "median[%!";
                List.iter (Printf.printf "%d,%!") intlst; 
                Printf.printf "]\n%!"; )intlstlst;  
            Printf.printf ")\n%!";
        ) intlstlstlst;
        Printf.printf "}\n%!";
        )node_file_median_chrom_seqdeli;
        Printf.printf "check seq list =>\n%!";
        List.iter (fun seq_lstlstlst ->
            Printf.printf "node{\n%!";
            List.iter (fun seqlstlst ->
            Printf.printf "file(\n%!";
            List.iter (fun seqlst -> 
                Printf.printf "median[%!";   List.iter (Sequence.printseqcode) seqlst; 
                Printf.printf "]\n%!") seqlstlst; 
            Printf.printf ")\n%!";
            )seq_lstlstlst;
            Printf.printf "}\n%!";
        ) node_file_median_chrom_seq;
         debug msg*)
        let node_file_median_seq = 
            List.map ( fun file_median_chrom_seq ->
                List.map (fun median_chrom_seq ->
                    List.map (fun chrom_seq -> Sequence.concat chrom_seq)
                    median_chrom_seq;
                ) file_median_chrom_seq;
            ) node_file_median_chrom_seq;(* file_seq_lstlstlst;*)
        in
        (* note: after concat, nodelst_medlst_seq is
        * { ( sequence of 1th median); ( sequence of 2th median ); .... }
        * deli_lstlst is
        * { 
            (delimiter lst of seqeunce of 1th median);
            (delimiter lst of sequence of 2th median);
            ......
        * *)
        (* debug msg
        Printf.printf "after concat: %!";
        List.iter (fun file_median_seq ->
            Printf.printf "node{\n%!";
            List.iter (fun seqlst ->
                Printf.printf "file (\n%!";
                List.iter (Sequence.printseqcode) seqlst;
                Printf.printf ")\n%!";
            ) file_median_seq;
            Printf.printf "}\n%!";
        ) node_file_median_seq;
         debug msg*)
        let new_nodedata_lst = 
        map3 (fun old_nodedata file_median_seq file_median_chrom_seqdeli -> 
            multi_to_single_chromosome old_nodedata file_median_seq file_median_chrom_seqdeli)
         nodes 
         node_file_median_seq 
         node_file_median_chrom_seqdeli 
        in
        new_nodedata_lst
    end
    else 
        nodes
    
let load_data ?(is_fixedstates=false) ?(silent=true) ?(classify=true) data =
    (* classify -- Not only we make the list of characters into sets of
       characers, but we also filter those characters that have weight 0. *)
    let classify = false in
    current_snapshot "Node.load_data start";
    let classify = (not (Data.has_dynamic data)) && classify in
    let make_set_of_list lst =
        List.fold_left 
            (fun acc x -> 
                if 0. = Data.get_weight x data then acc
                else All_sets.Integers.add x acc)
            (All_sets.Integers.empty)
            (lst)
    in
    let is_mem =
        (* We check for informative characters among all the terminals in data *)
        if classify then 
            (fun char ->
                Data.apply_boolean NonaddCS8.is_potentially_informative
                                   AddCS.is_potentially_informative data char
                && 
                    0. <> Data.get_weight char data)
        else (fun _ -> true)
    in
    (*let sign_dyna = List.filter is_mem data.Data.dynamics in*)
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
        and fixedstates = List.filter is_mem data.Data.fixed_states 
        and dynamics = List.filter is_mem data.Data.dynamics in
        let has_dynamic_mpl,has_dynamic_mal =
            List.fold_left
                (fun ((mpl,mal) as acc) code ->
                    match Hashtbl.find data.Data.character_specs code with
                    | Data.Dynamic ({Data.state = state} as s) when state = `Ml  ->
                        begin match s.Data.lk_model with
                            | Some m when m.MlModel.spec.MlModel.cost_fn = `MPL -> true,mal
                            | Some m when m.MlModel.spec.MlModel.cost_fn = `MAL ->  mpl,true
                            | _ -> assert false (* above pattern should be exhaustive *)
                        end
                    | _ -> acc)
                (false,false)
                (dynamics)
        in
        current_snapshot "start nonadd set2";
        let addvec,addgen = match add with
            | [] -> [],[]
            |  _ -> AddCS.split_vectorized_characters data add
        in
        let n8 = make_set_of_list n8
        and n16 = make_set_of_list n16
        and n32 = make_set_of_list n32
        and n33 = make_set_of_list n33
        and addvec = make_set_of_list addvec
        and addgen = make_set_of_list addgen in
        let cost_mode = match static_ml with
            | _::_                    -> `Likelihood
            | [] when has_dynamic_mpl -> `SumLikelihood
            | [] when has_dynamic_mal -> `Likelihood
            | _                       -> `NotLikelihood
        in
        current_snapshot "end nonadd set2";
        let r =
            generate_taxon classify addgen addvec n8 n16 n32 n33 sank dynamics
                           fixedstates kolmogorov static_ml data cost_mode
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
                    Status.create "Loading terminals" (Some ntaxa) "terminals loaded"
                in
                let cnt = ref 0 in
                (fun () ->
                    incr cnt;
                    Status.achieved status !cnt;
                    Status.full_report status),
                (fun () -> Status.finished status)
            else 
                (fun () -> ()),
                (fun () -> ())
        in
        let res =
            All_sets.IntegerMap.fold
                (fun x _ acc ->
                    let res = generate_taxon x acc in st (); res)
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

let pre = [ (Xml.Characters.cclass, `String Xml.Nodes.preliminary) ]
let fin = [ (Xml.Characters.cclass, `String Xml.Nodes.final) ]
let sing = [ (Xml.Characters.cclass, `String Xml.Nodes.single) ]

let to_single (pre_ref_codes, fi_ref_codes) combine_bl root parent mine =
    if debug_tosingle then begin
        Printf.printf "Node.ml to_single on parent and mine:\n%!";
        print parent;
        print mine;
    end;
    let rec cs_to_single (pre_ref_code, fi_ref_code) (root : cs option) parent_cs minet : cs =
        match parent_cs, minet with
            | Dynamic parentt, Dynamic minet ->
                let fst (a,_,_) = a and snd (_,a,_) = a in
                let root_pre,bl = match root with
                    | Some (Dynamic root) ->
                        begin match root.time with
                            | None,None,None  -> Some root.preliminary, None
                            | Some x,Some y,_ -> Some root.preliminary, Some (x +. y)
                            | _               -> assert( false )
                        end
                    | None ->
                        if combine_bl then begin match parentt.time with
                            | Some x,Some y, _   -> None, Some (x +. y)
                            | None, None, Some z -> None, Some z
                            | None, None, None   -> None, None
                            | _, _ , _           ->
                                 let (one,two,thr) = parentt.time in
                                 failwithf "Inconsistent branches for combine: %s, %s | %s"
                                    (match one with | Some x -> string_of_float x | None -> "none")
                                    (match two with | Some x -> string_of_float x | None -> "none")
                                    (match thr with | Some x -> string_of_float x | None -> "none")
                        end else if mine.min_child_code = parent.min_child_code then
                            None, fst parentt.time
                        else
                            None, snd parentt.time
                    | Some _ -> failwith "Inconsistent root passed to to_single"
                in
                let _, cost, res =
                    DynamicCS.to_single pre_ref_code
                        root_pre parentt.preliminary minet.preliminary bl
                in
                let cost = minet.weight *. cost in
                if debug_tosingle then Printf.printf "cs_to_single on mine=%d,parent=%d cost<-%f(old:%f),sum_cost <- %f(old:%f)\n%!"
                mine.taxon_code parent.taxon_code cost minet.cost cost minet.sum_cost;
                Dynamic {
                            preliminary = res; final = res;
                            (*we should NOT replace cost&sum_cost with cost to
                            * its parent,they are node cost and tree cost, and
                            * they should remain that way. if we really need a
                            * cost to parent, maybe we can add cost_to_parent to the
                            * data-structure.*)
                            cost = minet.cost;
                            sum_cost = minet.sum_cost;
                            weight = minet.weight;
                            time = minet.time;
                        }
                | _ -> match root with
                | Some mine -> mine
                | None -> minet
    in
    match root with
    | Some root ->
        let root_char_opt = List.map (fun c -> Some c) root.characters in
        let chars = map3 (cs_to_single (pre_ref_codes, fi_ref_codes) )
                         root_char_opt parent.characters mine.characters in
        { root with characters = chars; }
    | None ->
        let chars = map2 (cs_to_single (pre_ref_codes, fi_ref_codes) None ) 
                        parent.characters mine.characters in
        { mine with characters = chars; }

let readjust mode to_adjust ch1 ch2 parent mine = 
    let debug = false and debug2 = false in
    if debug then 
        Printf.printf "\n Node.ml readjust on mine:%d, parent:%d,ch1:%d, ch2:%d\n%!"
        mine.taxon_code parent.taxon_code ch1.taxon_code ch2.taxon_code;
    let ch1, ch2 =
        if ch1.min_child_code < ch2.min_child_code then
            ch1, ch2
        else ch2, ch1
    in
    let modified = ref All_sets.Integers.empty in
    let cs_readjust c1 c2 parent mine =
        match c1, c2, parent, mine with 
        | StaticMl c1, StaticMl c2, StaticMl parent, StaticMl mine -> 
            IFDEF USE_LIKELIHOOD THEN
                let min_bl = MlStaticCS.minimum_bl () in
                let t1,t2,t3_opt = match mine.time with 
                    | Some x,Some y,par -> max min_bl x,max min_bl y,par
                    | _ -> 
                        let (x,y) = MlStaticCS.estimate_time c1.preliminary c2.preliminary in
                        if debug_bl then Printf.printf "estimating BL: %f, %f\n%!" x y;
                        max min_bl x,max min_bl y, None
                in
                let m, prev_cost, cost, (t1,t2), res = 
                    MlStaticCS.readjust to_adjust !modified c1.preliminary
                                        c2.preliminary mine.preliminary t1 t2
                in
                modified := m;
                let cost = mine.weight *. cost in
                StaticMl 
                    { mine with 
                        preliminary=res;final=res;
                        cost=cost;sum_cost=cost;
                        time = Some t1, Some t2,t3_opt;
                    },
                cost
            ELSE
                failwith MlStaticCS.likelihood_error
            END
        | Dynamic c1, Dynamic c2, Dynamic parent, Dynamic mine ->
            begin match c1.preliminary,c2.preliminary with
                | DynamicCS.MlCS c1_pre, DynamicCS.MlCS c2_pre ->
                  IFDEF USE_LIKELIHOOD THEN
                    let ot1,ot2,ot3 = match mine.time with 
                        | Some x,Some y,Some z -> x, y, z
                        | _  -> assert false
                    in
                    let m,p_cost,n_cost,(t1,t2,t3),res =
                        DynamicCS.readjust_lk3 mode to_adjust !modified
                                mine.preliminary c1.preliminary c2.preliminary
                                parent.preliminary ot1 ot2 ot3
                    in
                    modified := m;
                    let cost = mine.weight *. n_cost in
                    let sumcost = cost +. c1.sum_cost +. c2.sum_cost in
(*                    Printf.printf "Optimized Cost %f(%f,%f,%f) --> %f(%f,%f,%f)\n%!"*)
(*                                  (p_cost *. mine.weight) ot1 ot2 ot3 cost t1 t2 t3;*)
                    Dynamic
                        { mine with
                            preliminary = res;
                            final = res;
                            cost = cost; 
                            sum_cost = sumcost;
                            time = Some t1, Some t2, Some t3;
                        },
                    sumcost
                  ELSE
                    failwith MlStaticCS.likelihood_error
                  END
                | _ -> 
                    let m, cost, sumcost, res = 
                        DynamicCS.readjust mode to_adjust !modified c1.preliminary
                                c2.preliminary parent.preliminary mine.preliminary
                    in
                    modified := m;
                    let cost = mine.weight *. cost in
                    let sumcost = mine.weight *. sumcost in
                    if debug2 then begin
                        Printf.printf " update mine with sum_cost\
                        <-- c1.sum_cost(%f)+c2.sum_cost(%f)+new node cost(%f) = %f,\
                        old sum_cost = %f, old node cost = %f;\n%!" 
                        c1.sum_cost c2.sum_cost cost sumcost mine.sum_cost mine.cost;
                        (*DynamicCS.print res;*)
                    end;
                    (*bug to fix: even when both the seq assignment and cost for 
                    * mine remain the same, the sum_cost could change. 
                    * in this case, modified is empty.
                    * in allDirChar.[adjust_node], we won't udpate node_data with
                    * the new seq & cost set if modified is empty.
                    *)
		            if (IntSet.is_empty !modified)&&((cost<>mine.cost)||(sumcost<>mine.sum_cost)) then 
				    failwith "node.ml readjust function, nothing changed from \
                    lower function, but different node cost or subtree cost";	
                    let res = 
                    Dynamic
                        { mine with
                            preliminary = res; 
                            final = res; 
                            cost = cost;
                            sum_cost = sumcost;
                            time=None,None,None;
                        }
                    in
                    if debug2 then begin
                        Printf.printf "end of node.readjust,return new med\n%!";
                    end;
                    res, sumcost
            end
        | Nonadd8 c1, Nonadd8 c2, Nonadd8 parent, Nonadd8 mine -> Nonadd8 mine, mine.sum_cost
        | Nonadd16 c1, Nonadd16 c2, Nonadd16 parent, Nonadd16 mine -> Nonadd16 mine, mine.sum_cost
        | Nonadd32 c1, Nonadd32 c2, Nonadd32 parent, Nonadd32 mine -> Nonadd32 mine, mine.sum_cost
        | AddVec c1, AddVec c2, AddVec parent, AddVec mine -> AddVec mine, mine.sum_cost
        | AddGen c1, AddGen c2, AddGen parent, AddGen mine ->  AddGen mine, mine.sum_cost
        | Sank c1, Sank c2, Sank parent, Sank mine -> Sank mine, mine.sum_cost
        | FixedStates c1, FixedStates c2, FixedStates parent, FixedStates mine -> FixedStates mine, mine.sum_cost
        | Kolmo c1, Kolmo c2, Kolmo parent, Kolmo mine -> Kolmo mine, mine.sum_cost
        | Set c1, Set c2, Set parent, Set mine -> Set mine, mine.sum_cost
        | _ -> failwith "Wrong type in matching, node.ml readjust"
    in
    if mine.total_cost = infinity then 
        mine, !modified
    else
        let _ = 
            if debug2 then Printf.printf "map4 on ch1,ch2,parent and mine's characters (len = %d)\n%!" 
            (List.length ch1.characters) 
        in
        let characters, sumcost_list = 
            List.split ( map4 cs_readjust ch1.characters ch2.characters
                             parent.characters mine.characters )
        in
        let node_cost = get_characters_cost characters in
        let total_cost = List.fold_left (fun acc x -> acc +. x ) 0. sumcost_list in 
	if debug then
		Printf.printf "end of Node.readjust, return mine with total_cost=%f(old \
        total_cost=%f),node_cost=%f(old node_cost=%f)\n\n%!" total_cost
        mine.total_cost node_cost mine.node_cost;
        let res = 
            { mine with characters = characters; 
                        total_cost = total_cost; 
                         node_cost = node_cost; }
        in
        res, !modified

let to_single_root (pre_ref_codes, fi_ref_codes) mine =
    to_single (pre_ref_codes, fi_ref_codes) true (Some mine) mine mine

(** [get_active_ref_code node_data] returns codes of
* all active chromosomes which are used for single state process
* on the subtree rooted at [node_data] *)
let get_active_ref_code node_data = 
    List.fold_left 
        (fun (acc_pre, acc_pre_child, acc_fi, acc_fi_child) cs -> 
             match cs with                   
             | Dynamic d -> 
                   let pre, pre_child = DynamicCS.get_active_ref_code d.preliminary in
                   let fi, fi_child = DynamicCS.get_active_ref_code d.final in
                   IntSet.union acc_pre pre, IntSet.union acc_pre_child pre_child,
                   IntSet.union acc_fi fi, IntSet.union acc_fi_child fi_child
             | _ -> acc_pre, acc_pre_child, acc_fi, acc_fi_child) 
        (IntSet.empty, IntSet.empty, IntSet.empty, IntSet.empty) 
        (node_data.characters)

let rec cs_to_formatter report_type node_name (pre_ref_codes,fi_ref_codes) d (cs,cs_single) parent_cs = 
    match cs, cs_single with
    | Nonadd8 cs, Nonadd8 _ ->
        begin match parent_cs with
        | None ->
            NonaddCS8.to_formatter 
                (NonaddCS8.to_formatter [] fin cs.final None d)
                pre cs.preliminary None d
        | Some ((Nonadd8 parent_cs), _) ->
            NonaddCS8.to_formatter
                (NonaddCS8.to_formatter [] fin cs.final (Some parent_cs.final) d)
                pre cs.preliminary (Some parent_cs.preliminary) d
        | _ -> assert false
        end
    | Nonadd16 cs, Nonadd16 _ ->
        begin match parent_cs with
        | None ->
            NonaddCS16.to_formatter 
                (NonaddCS16.to_formatter [] fin cs.final None d)
                pre cs.preliminary None d
        | Some ((Nonadd16 parent_cs), _) ->
            NonaddCS16.to_formatter
                (NonaddCS16.to_formatter [] fin cs.final (Some parent_cs.final) d)
                pre cs.preliminary (Some parent_cs.preliminary) d
        | _ -> assert false
        end
    | Nonadd32 cs, Nonadd32 _ ->
        begin match parent_cs with
        | None ->
            NonaddCS32.to_formatter 
                (NonaddCS32.to_formatter [] fin cs.final None d)
                pre cs.preliminary None d
        | Some ((Nonadd32 parent_cs), _) ->
            NonaddCS32.to_formatter
                (NonaddCS32.to_formatter [] fin cs.final (Some parent_cs.final) d)
                pre cs.preliminary (Some parent_cs.preliminary) d
        | _ -> assert false
        end
    | AddGen cs, AddGen _ ->
        begin match parent_cs with
            | None ->
                AddCS.General.to_formatter pre cs.preliminary None d
                @ AddCS.General.to_formatter fin cs.final None  d
            | Some ((AddGen parent_cs), _) ->
                AddCS.General.to_formatter pre cs.preliminary (Some parent_cs.preliminary) d
                @ AddCS.General.to_formatter fin cs.final (Some parent_cs.final) d
            | _ -> assert false
        end

    | AddVec cs, AddVec _ ->
        begin match parent_cs with
            | None ->
                AddCS.Vector.to_formatter pre cs.preliminary None d
                @ AddCS.Vector.to_formatter fin cs.final None  d
            | Some ((AddVec parent_cs), _) ->
                AddCS.Vector.to_formatter pre cs.preliminary (Some parent_cs.preliminary) d
                @ AddCS.Vector.to_formatter fin cs.final (Some parent_cs.final) d
            | _ -> assert false
        end
    | Sank cs, Sank _ ->
        begin match parent_cs with
            | None ->
                SankCS.to_formatter fin cs.final None d
            | Some ((Sank parent_cs), _) ->
                SankCS.to_formatter fin cs.final (Some parent_cs.final) d
            | _ -> assert false
        end
    | FixedStates cs, FixedStates _ ->
            Fixed_states.to_formatter report_type fin cs.final d
    | Dynamic cs, Dynamic cs_single ->
        let time = match cs.time with | (a,b,_) -> (a,b) in
        begin match parent_cs with
            | None ->
                DynamicCS.to_formatter report_type node_name pre_ref_codes pre cs.preliminary None time d
              @ DynamicCS.to_formatter report_type node_name fi_ref_codes fin cs.final None time d
              @ DynamicCS.to_formatter report_type node_name pre_ref_codes sing cs_single.preliminary None time d
            | Some ((Dynamic parent_cs), (Dynamic parent_cs_single)) ->
                DynamicCS.to_formatter report_type node_name pre_ref_codes pre cs.preliminary
                    (Some parent_cs.preliminary) time d
              @ DynamicCS.to_formatter report_type node_name fi_ref_codes fin cs.final
                    (Some parent_cs.final) time d
              @ DynamicCS.to_formatter report_type node_name pre_ref_codes sing cs_single.preliminary
                    (Some parent_cs_single.preliminary) time d
            | _ -> assert false
        end
    | Kolmo x, Kolmo x_single ->
        KolmoCS.to_formatter report_type pre_ref_codes pre x.preliminary d @
            KolmoCS.to_formatter report_type fi_ref_codes  fin x.final d @
            KolmoCS.to_formatter report_type pre_ref_codes sing x_single.preliminary d
    | Set x, Set x_single ->
        let attributes =
              [(Xml.Characters.name, `String (Data.code_character x.final.sid d))] in
        let sub a =
            (* we are getting rid of type Set in POY6 *)
            (cs_to_formatter report_type None (pre_ref_codes, fi_ref_codes) d a None)
        in
        let sub : (Xml.xml Sexpr.t list) = 
            List.map2 (fun a b -> `Set (sub (a, b))) 
                      x.final.set x_single.final.set 
        in
        [`Single (Xml.Characters.set, attributes, `Set sub)]
    | StaticMl cs, StaticMl _ ->
        IFDEF USE_LIKELIHOOD THEN
            let time = match cs.time with | (a,b,_) -> (a,b) in
            MlStaticCS.to_formatter pre cs.preliminary time d
        ELSE
            failwith MlStaticCS.likelihood_error
        END
    | (FixedStates _ |StaticMl _ |Kolmo _ |Dynamic _ |Sank _ |AddGen _ | AddVec _
        |Nonadd32 _ |Nonadd16 _ |Nonadd8 _ |Set _), _ -> assert false

(* Compute total recost of the NODE, NOT THE SUBTREE*)
let cmp_node_recost node_data = 
    List.fold_left 
        (fun recost cs -> match cs with
            | Dynamic dyn -> recost +. DynamicCS.total_recost dyn.preliminary
            | _ -> recost) 
        0.0 node_data.characters

(* Compute total recost of the subtree rooted by this NODE *)
let cmp_subtree_recost node_data =
    if debug_formatter then
        Printf.printf "node.ml cmp_subtree_recost:%!";
    List.fold_left 
        (fun subtree_recost cs -> match cs with 
            | Dynamic dyn ->
                let res = subtree_recost +. (DynamicCS.subtree_recost dyn.preliminary) in
                if debug_formatter then
                    Printf.printf "res = %f + %f \n%!"  subtree_recost
                        (DynamicCS.subtree_recost dyn.preliminary);
                res
            | _ -> subtree_recost) 
        0.0
        node_data.characters

let to_formatter_single report_type (pre_ref_codes,fi_ref_codes) acc d (node_data,node_single) node_id parent_data =
    let get_node_name id =
        try Data.code_taxon id d
        with | Not_found -> string_of_int id 
    in
    let node_name = get_node_name node_id in
    let get_parent item = match parent_data with
        | None -> None
        | Some (a, b) ->
            (Some ((List.nth a.characters item), (List.nth b.characters item)))
    in
    if debug_formatter then begin
        Printf.printf "to_formatter_single on node %s\n%!" node_name;
        print node_data;
        print node_single;
    end;
    let lazy_children () =
        `Set (fst
            (List.fold_left
                (fun (acc, item) cs ->
                    let parent_data = get_parent item in
                    if debug_formatter then
                        Printf.printf "Delayed function of to_formatter_single on node %s\n%!" node_name;
                    let res =
                         cs_to_formatter report_type (Some node_name) (pre_ref_codes, fi_ref_codes) d cs parent_data
                    in
                    res @ acc, item + 1)
                ([], 0)
                (List.map2 (fun a b -> a, b) node_data.characters node_single.characters)))
    in
    let module T = Xml.Nodes in
    (RXML
        -[T.node]
            ([T.cost] = 0.)
            ([T.recost] = 0.)
            ([T.node_cost] = [`Float node_data.node_cost])
            ([T.name] = [`String node_name])
            ([T.child1_name] = "")
            ([T.child2_name] = "")
            ([T.nce] = [`Int node_data.num_child_edges])
            ([T.notu] = [`Int node_data.num_otus])
            ([acc])
            { `Delayed lazy_children } --)

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

let to_formatter_subtree diag_report_type (pre_ref_codes, fi_ref_codes)
        acc d (node_data, node_single) node_id (child1_id,  child1_node_data)
        (child2_id,  child2_node_data) (parent_node_data_opt : (node_data *
        node_data) option) : Xml.xml =
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
    if debug_formatter then begin
        Printf.printf "node.ml to_formatter_subtree, node name : %s, child1 is %s, child2 is %s\n%!" node_name child1_name child2_name;
    print node_data;
    print node_single;
    end;
    let child1_recost = cmp_subtree_recost child1_node_data in 
    let child2_recost = cmp_subtree_recost child2_node_data in 
    let node_recost = (cmp_node_recost node_data) in
    let subtree_recost = 
        child1_recost +. 
        child2_recost +. 
        node_recost in 
    if debug_formatter then
        Printf.printf "%f + %f + %f => subtree_recost=%f \n%!" 
    child1_recost child2_recost node_recost subtree_recost;
    let module T = Xml.Nodes in
    let attr : Xml.attributes = 
        (AXML
            ([T.cost] =
                 [match parent_node_data_opt with 
                  | None -> `Float 0. (*that's why in diagnosis, the "Cost"
                  under "root" is 0, while "Rearrangement cost" is not*)
                  | _ -> `Float node_data.total_cost])
            ([T.recost] = [`Float subtree_recost])
            ([T.node_cost] = [`Float subtree_recost])
            ([T.name] = [`String node_name])
            ([T.child1_name] = [`String child1_name])
            ([T.child2_name] = [`String child2_name])
            ([T.nce] = [`Int node_data.num_child_edges])
            ([T.notu] = [`Int node_data.num_otus])
            ([acc]))
    in
    let children : Xml.xml Xml.contents =
        `Delayed (fun () ->
        `Set (fst (List.fold_left 
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
             if debug_formatter then 
                Printf.printf "Delayed function in to_formatter_subtree on node.%s"
             node_name;
             let res = cs_to_formatter diag_report_type None (pre_ref_codes, fi_ref_codes) d cs parent_cs in 
             (res @ acc), (idx + 1)) ([], 0) 
             (List.map2 (fun a b -> a, b) node_data.characters 
             node_single.characters))))
    in
    (RXML -[T.node] ([attr]) { children } --)

let to_xml data ch node = ()

type 'a converter = 'a r -> cs list -> cs list

let listify a acc = a :: acc
let lnon8 cs acc = listify (Nonadd8 cs) acc
let lnon16 cs acc = listify (Nonadd16 cs) acc
let lnon32 cs acc = listify (Nonadd32 cs) acc
let laddvec cs acc = listify (AddVec cs) acc
let laddgen cs acc = listify (AddGen cs) acc
let lsank cs acc = listify (Sank cs) acc
let lfixedstates cs acc = listify (FixedStates cs) acc
let ldynamic cs acc = listify (Dynamic cs) acc
let lkolmo cs acc = listify (Kolmo cs) acc
let lstaticml cs acc = listify (StaticMl cs) acc

let rec convert_data
        ?(tnon8=lnon8)
        ?(tnon16=lnon16)
        ?(tnon32=lnon32)
        ?(taddvec=laddvec)
        ?(taddgen=laddgen)
        ?(tsank=lsank)
        ?(tfixedstates=lfixedstates)
        ?(tdynamic=ldynamic)
        ?(tkolmo=lkolmo)
        ?(tstaticml=lstaticml)
        node =
    let rec conv cs acc = match cs with
             | Nonadd8 cs -> tnon8 cs acc
             | Nonadd16 cs -> tnon16 cs acc
             | Nonadd32 cs -> tnon32 cs acc
             | AddGen cs -> taddgen cs acc
             | AddVec cs -> taddvec cs acc
             | Sank cs -> tsank cs acc
             | Dynamic cs -> tdynamic cs acc
             | FixedStates cs -> tfixedstates cs acc
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


let cardinal item = match item with
    | StaticMl c ->
        IFDEF USE_LIKELIHOOD THEN
            MlStaticCS.cardinal c.preliminary
        ELSE
            failwith MlStaticCS.likelihood_error
        END
    | Nonadd8 c  -> NonaddCS8.cardinal c.preliminary
    | Nonadd16 c -> NonaddCS16.cardinal c.preliminary
    | Nonadd32 c -> NonaddCS32.cardinal c.preliminary
    | AddVec c   -> AddCS.Vector.cardinal c.preliminary
    | AddGen c   -> AddCS.General.cardinal c.preliminary
    | FixedStates c -> Fixed_states.cardinal c.preliminary   
    | Sank c     -> SankCS.cardinal c.preliminary
    | Dynamic c  -> DynamicCS.cardinal c.preliminary
    | Kolmo c    -> KolmoCS.cardinal c.preliminary
    | Set s      -> 0


let rec filter_character_codes (codes : All_sets.Integers.t) item = match item with
    | StaticMl c ->
        IFDEF USE_LIKELIHOOD THEN
            StaticMl (do_filter MlStaticCS.cardinal MlStaticCS.f_codes c codes)
        ELSE
            failwith MlStaticCS.likelihood_error
        END
    | Nonadd8 c ->
          Nonadd8 (do_filter NonaddCS8.cardinal NonaddCS8.f_codes c codes)
    | Nonadd16 c ->
          Nonadd16 (do_filter NonaddCS16.cardinal NonaddCS16.f_codes c codes)
    | Nonadd32 c ->
          Nonadd32 (do_filter NonaddCS32.cardinal NonaddCS32.f_codes c codes)
    | AddGen c ->
          AddGen (do_filter AddCS.General.cardinal AddCS.General.f_codes c codes)
    | AddVec c ->
          AddVec (do_filter AddCS.Vector.cardinal AddCS.Vector.f_codes c codes)
    | Sank c -> 
          Sank (do_filter SankCS.cardinal SankCS.f_codes c codes)
    | FixedStates c ->  
          FixedStates (do_filter Fixed_states.cardinal Fixed_states.f_codes c codes)
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
          failwith MlStaticCS.likelihood_error
        END
    | Nonadd8 c ->
          Some (Nonadd8 (do_filter NonaddCS8.cardinal NonaddCS8.f_codes_comp c codes))
    | Nonadd16 c ->
          Some (Nonadd16 (do_filter NonaddCS16.cardinal NonaddCS16.f_codes_comp c codes))
    | Nonadd32 c ->
          Some (Nonadd32 (do_filter NonaddCS32.cardinal NonaddCS32.f_codes_comp c codes))
    | AddGen c ->
          Some (AddGen (do_filter AddCS.General.cardinal AddCS.General.f_codes_comp c codes))
    | AddVec c ->
          Some (AddVec (do_filter AddCS.Vector.cardinal AddCS.Vector.f_codes_comp c codes))
    | Sank c -> 
          Some (Sank (do_filter SankCS.cardinal SankCS.f_codes_comp c codes))
    | FixedStates c ->
        Some (FixedStates (do_filter Fixed_states.cardinal Fixed_states.f_codes_comp c codes))
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

let pp_int_set chan set =
    All_sets.Integers.iter (fun i -> Printf.fprintf chan "%d, " i) set

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
    let chars = mmap (filter_character_codes_complement codes) n.characters in
    { n with characters = chars }


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

let remove_exclude n = { n with exclude_sets = []; exclude_info = [] }


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
                        failwith MlStaticCS.likelihood_error
                    END
                | Nonadd8 r  -> NonaddCS8.cardinal r.preliminary
                | Nonadd16 r -> NonaddCS16.cardinal r.preliminary
                | Nonadd32 r -> NonaddCS32.cardinal r.preliminary
                | AddGen r   -> AddCS.General.cardinal r.preliminary
                | AddVec r   -> AddCS.Vector.cardinal r.preliminary
                | Sank r     -> SankCS.cardinal r.preliminary
                | Dynamic r  -> DynamicCS.cardinal r.preliminary
                | FixedStates r -> Fixed_states.cardinal r.preliminary
                | Kolmo r    -> KolmoCS.cardinal r.preliminary
                | Set s      -> internal_n_chars acc s.preliminary.set
            in
            internal_n_chars (acc + count) cs
(** [n_chars chars] returns the total number of characters in the list of
    characters [chars] *)
let rec n_chars ?(acc=0) n =
    internal_n_chars acc n.characters

module T = struct
    let add_exclude = add_exclude
    let remove_exclude = remove_exclude
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
        | AddGenU of AddCS.General.t ru
        | AddVecU of AddCS.Vector.t ru
        | SankU of SankCS.t ru
        | DynamicU of DynamicCS.u ru
        | FixedStatesU of Fixed_states.t_w_seqtbl ru
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
                        | AddVecU _ 
                        | AddGenU _ 
                        | StaticMlU _
                        | KolmoU _
                        | FixedStatesU _ 
                        | SankU _ -> None
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
            | AddVecU x -> 
                    failwith "TODO Node.Union.poly_saturation"
            | AddGenU x -> 
                    failwith "TODO Node.Union.poly_saturation"
                    (*
                    let card = float_of_int (AddCS.cardinal x) in
                    AddCS.poly_saturation x *. card, card +. lenacc
                    *)
            | KolmoU _ 
            | FixedStatesU _ 
            | SankU _ -> 
                    failwith "TODO Node.Union.poly_saturation"
                    (*
                    let card = float_of_int (SankCS.cardinal x) in
                    SankCS.poly_saturation x *. card, card +. lenacc
                    *)
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
            | AddVec cs, AddVecU u1, AddVecU u2 -> 
                    AddVecU { ch = cs.preliminary; u_weight = u1.u_weight }
            | AddGen cs, AddGenU u1, AddGenU u2 -> 
                    AddGenU { ch = cs.preliminary; u_weight = u1.u_weight }
            | Sank cs, SankU u1, SankU u2 ->
                    SankU { ch = cs.preliminary; u_weight = u1.u_weight }
            | FixedStates cs, FixedStatesU u1, FixedStatesU u2 ->
                    FixedStatesU { ch = cs.preliminary; u_weight = u1.u_weight }
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
                    failwith MlStaticCS.likelihood_error
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
        let debug = false in
        let code = 
            match taxon_code with
            | None -> c.taxon_code
            | Some v -> v
        in
        if debug then Printf.printf "node.Union.leaf create union for Node:%d\n%!" code;
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
            | FixedStates c ->(* TODO *)
                    (FixedStatesU { ch = c.preliminary; u_weight = c.weight}) :: acc
            | AddVec c ->
                    (AddVecU { ch = c.preliminary; u_weight = c.weight }) :: acc
            | AddGen c ->
                    (AddGenU { ch = c.preliminary; u_weight = c.weight }) :: acc
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
                    (* Do nothing special *)
                    (StaticMlU { ch = c.preliminary; u_weight = c.weight }) :: acc
                ELSE
                    failwith MlStaticCS.likelihood_error
                END
            | Set _ -> failwith "Node.Union.leaf TODO"
        in
        let nc = List.fold_left single_leaf [] c.characters in
        let nc = List.rev nc in
        if debug then Printf.printf "end of node.Union.leaf \n%!";
        { charactersu = nc; min_child_codeu = code }

    let distance a b =
        let debug = false in
        if debug then  
            Printf.printf "node.Union.distance on Anode:%d and Bnode:%d --> %!"
            a.min_child_codeu b.min_child_codeu;
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
                    (* TODO MlStaticCS.distance_union a b, at, bt *)
            | (AddGenU a) :: at, (AddGenU b) :: bt ->
                    distance acc at bt
            | (AddVecU a) :: at, (AddVecU b) :: bt ->
                    distance acc at bt
                    (* TODO AddCS.distance a b, at, bt *)
            | (FixedStatesU a) :: at, (FixedStatesU b) :: bt ->
                    distance acc at bt
                    (* TODO Fixed_states.distance a b, at, ct *)
            | (SankU a) :: at, (SankU b) :: bt ->
                    distance acc at bt
                    (* TODO SankCS.distance a b, at, ct *)
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
            | (AddGenU _) :: _, _
            | (AddVecU _) :: _, _
            | (SankU _) :: _, _
            | (FixedStatesU _) :: _, _
            | (DynamicU _) :: _, _
            | (KolmoU _) :: _, _ 
            | (StaticMlU _) :: _, _ 
            | [], _ -> failwith "Node.Union.distance TODO"
        in
        let res = distance 0.0 a.charactersu b.charactersu in
        res

    let distance_node code a b =
        let debug = false in
        if debug then Printf.printf "node.Union.distance_node,%!";
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
                failwith MlStaticCS.likelihood_error
            END
        | (AddGenU a), (AddGenU b) ->
                AddCS.General.compare_data a.ch b.ch
        | (AddVecU a), (AddVecU b) ->
                AddCS.Vector.compare_data a.ch b.ch
        | (SankU a), (SankU b) ->
                SankCS.compare_data a.ch b.ch
        | (KolmoU a), (KolmoU b) ->
                0
        | (DynamicU a), (DynamicU b) ->
                DynamicCS.compare_union a.ch b.ch
        | _ -> failwith "Node.Union.distance TODO"

    let rec compare a b = match a, b with
        | [], [] -> 0
        | [], _  -> -1
        | _, []  -> 1
        | ha::ta,hb::tb ->
            begin match single_compare ha hb with
                | 0 -> compare ta tb
                | x -> x
            end

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
            weight = 1.0 ; time = None,None,None }
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
        | [] -> `NotLikelihood
    in
    let leaves = List.map (fun (x, y) -> (x, [], y)) leaves in
    let load_clade node_data node_id =
        List.map
            (fun (leaf_id, charlist, leaf_data) ->
                let char =
                    let leaf = if List.mem leaf_id leaves_id then 2 else 1 in
                    NonaddCS8.of_list
                        [(node_id + starting, (node_id + starting,leaf), 0.)]
                in
                (leaf_id, char :: charlist, leaf_data))
            node_data
    in
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
                                time = None,None,None;
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

(*when we have non-zero diagonal, cost from aligning two children of the root will
* be higher than distance of two aligned children of the root. this function
    * return the difference between them.*)
let extra_cost_from_root n treecost =
    let debug = false and debug2 = false in
    if debug then begin
        Printf.printf "node.ml extra cost from root, treecost=%f, %!" treecost;
	    if debug2=true then begin
        Printf.printf "root nodedata:%!";
        print n;
        end;
    end;
    if treecost = 0. then 
        (*when we build wagner tree, we add nodes one by one. so there will be
        * time when the tree has only one node, which is a leafnode, tree cost
        * should be 0 at this point.
        * extra_cost_from_root will get distance of root's alied children,
        * since they carry same sequence as the leaf itself,
        * [get_extra_cost_for_root] might return some cost. *)
        0.
    else 
        let acc_cost_cs acc item =
             match item with 
            | Sank x -> 
                    let ec = SankCS.get_extra_cost_for_root x.preliminary in
                    if debug2 then Printf.printf "sankCS,acc(%f) += %d\n%!" acc ec;
                    acc +. (float_of_int ec)
            | Dynamic x ->
                    let disc = DynamicCS.extra_cost_for_root x.preliminary  in
                    if debug2 then Printf.printf "DynamicCS,acc(%f) += %f\n%!" acc disc;
                    acc +. disc
            | _ -> 0.0
        in
        let cost = List.fold_left acc_cost_cs 0.0 n.characters in
        if debug then Printf.printf "return %f\n%!" cost;
        (*I know for dynamic character we can just return distance cost as new cost for root, 
        * there is no need to get cost difference between align cost and distance cost. but for
        * sankoff charactor type, there is no algned children. *)
        cost
    



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
        let median = median
        let apply_time = apply_time
        let extract_states a d _ c n = extract_states a d c n
        let min_prior = prior
        let get_times_between = get_times_between_plus_codes 
        let final_states _ = final_states
        let uppass_heuristic pcode ptime mine a b = mine
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
        let classify_data = classify_data
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
        let edge_iterator = edge_iterator
        let readjust = readjust
        module T = T
        module Union = Union
        let for_support = for_support
        let root_cost = root_cost
        let extra_cost_from_root = extra_cost_from_root
        let tree_cost a b =
		let tc = (root_cost b) +. (total_cost a b) in
		let ec = (extra_cost_from_root b tc) in
		tc -. ec
        let to_single root _ a _ b sets =
            let combine = match root with
                | Some _ -> true
                | None   -> false
            in
            to_single sets combine root a b
        let get_nonadd_8 _ = get_nonadd_8 
        let get_nonadd_16 _ = get_nonadd_16 
        let get_nonadd_32 _ = get_nonadd_32 
        let get_addgen _ = get_addgen
        let get_addvec _ = get_addvec
        let get_sank _ = get_sank
        let get_dynamic _ = get_dynamic
        let get_fixedstates _ = get_fixedstates
        let get_mlstatic _ = get_mlstatic
        let compare a b = 
            let rec aux_cmt a b =
                match a, b with
                | (Nonadd8 ha) :: ta, (Nonadd8 hb) :: tb ->
                        let x = compare ha.weight hb.weight in
                        if x = 0 then
                            let cmp = compare ha.preliminary hb.preliminary in
                            if cmp = 0 then 
                                aux_cmt ta tb
                            else cmp
                        else x
                | (Nonadd8 ha) :: ta, _ -> -1
                | (Nonadd16 ha) :: ta, (Nonadd16 hb) :: tb ->
                        let x = compare ha.weight hb.weight in
                        if x = 0 then
                            let cmp = compare ha.preliminary hb.preliminary in
                            if cmp = 0 then 
                                aux_cmt ta tb
                            else cmp
                        else x
                | (Nonadd16 ha) :: ta, _ -> -1
                | (Nonadd32 ha) :: ta, (Nonadd32 hb) :: tb ->
                        let x = compare ha.weight hb.weight in
                        if x = 0 then
                            let cmp = compare ha.preliminary hb.preliminary in
                            if cmp = 0 then 
                                aux_cmt ta tb
                            else cmp
                        else x
                | (Nonadd32 ha) :: ta, _ -> -1
                | (AddGen ha) :: ta, (AddGen hb) :: tb ->
                        let x = compare ha.weight hb.weight in
                        if x = 0 then
                            let cmp = compare ha.preliminary hb.preliminary in
                            if cmp = 0 then 
                                aux_cmt ta tb
                            else cmp
                        else x
                | (AddGen ha) :: ta, _ -> -1
                | (AddVec ha) :: ta, (AddVec hb) :: tb ->
                        let x = compare ha.weight hb.weight in
                        if x = 0 then
                            let cmp = compare ha.preliminary hb.preliminary in
                            if cmp = 0 then 
                                aux_cmt ta tb
                            else cmp
                        else x
                | (AddVec ha) :: ta, _ -> -1
                | (Sank ha) :: ta, (Sank hb) :: tb ->
                        let x = compare ha.weight hb.weight in
                        if x = 0 then
                            let cmp = SankCS.compare_data ha.preliminary hb.preliminary in
                            if cmp = 0 then 
                                aux_cmt ta tb
                            else cmp
                        else x
                | (Sank ha) :: ta, _ -> -1
                | (FixedStates ha) :: ta, (FixedStates hb) :: tb ->
                            let cmp = Fixed_states.compare_data ha.preliminary hb.preliminary in
                            if cmp = 0 then 
                                aux_cmt ta tb
                            else cmp
                | (FixedStates ha) :: ta, _ -> -1
                | (Dynamic ha) :: ta, (Dynamic hb) :: tb ->
                        let x = compare ha.weight hb.weight in
                        if x = 0 then
                            let cmp = compare ha.preliminary hb.preliminary in
                            if cmp = 0 then 
                                aux_cmt ta tb
                            else cmp
                        else x
                | (Dynamic ha) :: ta, _ -> -1
                | (Kolmo ha) :: ta, (Kolmo hb) :: tb ->
                        let x = compare ha.weight hb.weight in
                        if x = 0 then
                            let cmp = compare ha.preliminary hb.preliminary in
                            if cmp = 0 then 
                                aux_cmt ta tb
                            else cmp
                        else x
                | (Kolmo ha) :: ta, _ -> -1
                | (StaticMl ha) :: ta, (StaticMl hb) :: tb ->
                    IFDEF USE_LIKELIHOOD THEN
                        let x = compare ha.weight hb.weight in
                        if x = 0 then
                            (* specialized compare --only compares general model *)
                            let cmp = MlStaticCS.compare ha.preliminary hb.preliminary in
                            if cmp = 0 then 
                                aux_cmt ta tb
                            else cmp
                        else x
                    ELSE 
                        0
                    END
                | (StaticMl ha) :: ta, _ -> -1
                | (Set ha) :: ta, (Set hb) :: tb ->
                        let x = compare ha.weight hb.weight in
                        if x = 0 then
                            let cmp = compare ha.preliminary hb.preliminary in
                            if cmp = 0 then 
                                aux_cmt ta tb
                            else cmp
                        else x 
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

(*this function [total_cost_of_type] is only being called by allDirChar.ml [check_cost] for static
* charactors. so why we need to match Dynamic/etc/ here?*)
let total_cost_of_type t n =
    let debug = false in
    let rec total_cost_cs acc item =
        let single = match item, t with
            | Nonadd8 x, `Nonadd -> 
                    if debug then Printf.printf
                    "total_cost_of_type,Nonadd8,%f \n%!" x.sum_cost; 
            x.sum_cost 
            | Nonadd16 x, `Nonadd -> 
                    if debug then Printf.printf
                    "total_cost_of_type,Nonadd16,%f\n%!" x.sum_cost;
                    x.sum_cost
            | Nonadd32 x, `Nonadd -> 
                    if debug then Printf.printf
                    "total_cost_of_type,Nonadd32,%f \n%!" x.sum_cost;
                    x.sum_cost
            | AddGen x, `Add -> 
                    if debug then Printf.printf
                    "total_cost_of_type,AddGen,%f \n%!" x.sum_cost;
                    x.sum_cost 
            | AddVec x, `Add -> 
                    if debug then Printf.printf
                    "total_cost_of_type,AddVec,%f\n%!" x.sum_cost;
                    x.sum_cost 
            | Sank x, `Sank ->
                    let ec = SankCS.get_extra_cost_for_root x.preliminary in
                    (x.sum_cost -. (float_of_int ec) *. x.weight) 
            | StaticMl x, `StaticMl ->
                IFDEF USE_LIKELIHOOD THEN
                    x.cost
                ELSE
                    0.0
                END
            | Set x, t -> List.fold_left total_cost_cs acc x.preliminary.set
            | FixedStates x, `FixedStates ->
                    x.sum_cost
            | Dynamic x, t ->
                    begin match x.preliminary, t with
                    | DynamicCS.MlCS _, `Ml when n.cost_mode = `Likelihood ->
                        x.cost 
                    | DynamicCS.MlCS _, `Ml when n.cost_mode = `SumLikelihood ->
                        x.sum_cost 
                    | DynamicCS.MlCS _, `Ml -> assert false
                    | DynamicCS.SeqCS _, `Seq ->
                        x.sum_cost 
                    | DynamicCS.BreakinvCS _, `Breakinv ->
                        x.sum_cost 
                    | DynamicCS.ChromCS _, `Chrom ->
                        x.sum_cost
                    | DynamicCS.AnnchromCS _, `Annchrom ->
                        x.sum_cost 
                    | DynamicCS.GenomeCS _, `Genome ->
                        x.sum_cost 
                    | _ -> 0.0
                end
            | Kolmo x, `Kolmo -> x.sum_cost 
            | _,_ ->  0.0
        in
        if debug then
            info_user_message "acc=%f, %s contributed additional cost %f for %s" 
                              acc (cs_string item) (single) (type_string t) ;
        acc +. single
    in
    List.fold_left total_cost_cs 0.0 n.characters



