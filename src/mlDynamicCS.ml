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

let () = SadmanOutput.register "MlDynamicCS" "$Revision: 2650 $"

(*---- non-external helper functions/settings *)
open Numerical.FPInfix
open FloatSequence  (* the modules inside are descriptive enough where opening
                       this module is not determental to readability. *)

let debug      = false   (* show debug information during median calculations *)
let debug_est  = false   (* output for optimization of data *)
let verify     = false   (* Verify pure ocaml and C implementations of MPL *)
let pure_ocaml = false   (* Use a pure ocaml implementation of MPL *)

let (-->) a b = b a

let failwith_todo func =
    failwith ("Your lazy developer needs to write: MlDynamicCS."^func)

let debug_printf = 
    if debug then (fun format -> print_string format)
             else (fun _ -> ())

IFDEF USE_LIKELIHOOD THEN

(*---- type for the module *)
type 'a align = { ss : 'a array; }

(* Define how we should align and create the medians; 
 *   MPLAlign    - Two transformation matrices are created for each sequence,
 *                 these measure the probability of each transition, and when
 *                 find the minimum value (X->I)(Y->I) for the assignment I.
 *                 This is only turned on when 'pure_ocaml' above is set to
 *                 true, otherwise, CMPL is used.
 *   CMPLAlign   - is an C implementation of MPL. (Default).
 *   Verify      - Checks MPL and CMPL against each other for verification. This
 *                 is turned on above when 'verify' is set.  *)
type r = | MPLAlign  of MPLAlign.s align
         | CMPLAlign of CMPLAlign.s align
         | Verify    of (CMPLAlign.s * MPLAlign.s) align
    (** Define characters for maximum parsimonious likelihood or ancestral
        likelihood. These characters find the best assignment of a transformation
        from the two children. If X and Y are the states of the children, $P(a,b|t)$
        is the probability of state b over branch length t, then MPL is defined as,
            $\max_{\alpha = \{ACTG-\}} P(X,\alpha|t_1) * P(Y,\alpha|t_2)$ *)

type t = { model  : dyn_model;
            data  : r;
            codes : int array;
            code  : int;
            times : float * float; (* branch lengths of composition *)
            cost  : float; }

(*---- basic function set for gathering data of the module *)

(* although we have a cost matrix for the characters, this method is used to
   ensure that we do not use it without a branch length present *)
let get_cm _ = assert false

let print_combined_cm matrix = 
    for i = 0 to (Array.length matrix) -1 do
        for j = 0 to (Array.length matrix.(i))-1 do
            Printf.printf "[%f|%d]\t" (fst matrix.(i).(j)) (snd matrix.(i).(j));
        done;
        Printf.printf "\n%!";
    done;
    ()

(* compare the cost matrix and assignments generated from FloatSequence.A.get_cm *)
let compare_combined_cm r1 r2 : bool =
    let ret = ref true in
    for i = 0 to (Array.length r1) -1 do
        for j = 0 to (Array.length r1.(i))-1 do
            let t1 = ((fst r1.(i).(j)) = (fst r2.(i).(j)))
            and t2 = ((snd r1.(i).(j)) = (snd r2.(i).(j))) in
            ret := !ret & t1 & t2;
        done;
    done;
    !ret

let compare a b = 
    let compare_array2 f a b = 
        if (Array.length a) = (Array.length b) then begin
            let res = ref true in
            for i = 0 to (Array.length a)-1 do
                res := !res && (f a.(i) b.(i))
            done;
            !res
        end else begin
            false
        end
    in
    let model_compare = MlModel.compare a.model.static b.model.static in
    let data_compare = match a.data,b.data with
        |  MPLAlign ad,  MPLAlign bd -> compare_array2   MPLAlign.compare ad.ss bd.ss
        | CMPLAlign ad, CMPLAlign bd -> compare_array2  CMPLAlign.compare ad.ss bd.ss
        |    Verify ad,    Verify bd ->
            (compare_array2 (fun (x,_) (y,_) -> CMPLAlign.compare x y) ad.ss bd.ss) &&
            (compare_array2 (fun (_,x) (_,y) ->  MPLAlign.compare x y) ad.ss bd.ss)
        | (MPLAlign _ | Verify _ | CMPLAlign _) , _ -> false
    in
    0 = model_compare && data_compare

let pp_seq model chan seq = Sequence.print chan seq model.alph

let to_string t = failwith_todo "to_string"

let alph t = t.model.alph

let code t = t.code

let get_codes t = t.codes

let cardinal t = Array.length t.codes

let encoding e t = match t.data with
    | MPLAlign r ->
        Array.fold_left
            (fun acc x ->
                Sequence.encoding (e) (MPLAlign.seq_of_s x))
            (0.0)
            (r.ss)
    | CMPLAlign r ->
        Array.fold_left
            (fun acc x ->
                Sequence.encoding (e) (CMPLAlign.seq_of_s x))
            (0.0)
            (r.ss)
    | Verify r ->
        Array.fold_left
            (fun acc (cmpl,mpl) ->
                let one = Sequence.encoding (e) (MPLAlign.seq_of_s mpl)
                and oth = Sequence.encoding (e) (CMPLAlign.seq_of_s cmpl) in
                assert( one =. oth );
                one)
            (0.0)
            (r.ss)

let name_string t = match t.data with
    | MPLAlign _  -> "Pure OCaml Maximum Parsimonious Dynamic Likelihood"
    | CMPLAlign _ -> "Maximum Parsimonious Dynamic Likelihood" 
    | Verify    _ -> "Verifying Maximum Parsimonious Dynamic Likelihood" 

let total_cost t = t.cost

let model t = t.model

let static_model t = t.model.static

let s_of_seq seq : Sequence.s array =
    seq --> SeqCS.get_sequences --> Array.to_list --> Array.concat

let leaf_sequences t = match t.data with
    | MPLAlign r ->
        Array_ops.fold_right_2
            (fun acc x y ->
                let y = [| `DO (MPLAlign.seq_of_s y) |] in
                All_sets.IntegerMap.add x y acc)
            (All_sets.IntegerMap.empty)
            (t.codes)
            (r.ss)
    | CMPLAlign r ->
        Array_ops.fold_right_2
            (fun acc x y ->
                let y = [| `DO (CMPLAlign.seq_of_s y) |] in
                All_sets.IntegerMap.add x y acc)
            (All_sets.IntegerMap.empty)
            (t.codes)
            (r.ss)
    | Verify r ->
        Array_ops.fold_right_2
            (fun acc x (y1,y2) ->
                let y = 
                    assert(0 = Sequence.compare (CMPLAlign.seq_of_s y1) (MPLAlign.seq_of_s y2) );
                    [| `DO (CMPLAlign.seq_of_s y1) |] 
                in
                All_sets.IntegerMap.add x y acc)
            (All_sets.IntegerMap.empty)
            (t.codes)
            (r.ss)


(*---- to formatter, and printing functions *)
let to_formatter attr mine par_opt (t1,t2) d : Xml.xml Sexpr.t list = 
    let str_time = function | Some x -> `Float x | None -> `String "None" in
    let alphabet= alph mine in
    let model_d = MlModel.to_formatter (static_model mine) in
    let seq_f () = match mine.data with
        | MPLAlign r ->
            let seq_lst =
                Array.fold_left
                    (fun acc s ->
                        s --> MPLAlign.seq_of_s
                          --> Sequence.del_first_char
                          --> fun x -> (Sequence.to_formater x alphabet)::acc)
                    []
                    r.ss
            in
            String.concat "#" seq_lst
        | CMPLAlign r ->
            let seq_lst =
                Array.fold_left
                    (fun acc s ->
                        s --> CMPLAlign.seq_of_s
                          --> Sequence.del_first_char
                          --> fun x -> (Sequence.to_formater x alphabet)::acc)
                    []
                    r.ss
            in
            String.concat "#" seq_lst
        | Verify r -> 
            let seq_lst =
                Array.fold_left
                    (fun acc s ->
                        s --> (fun x -> x --> snd --> MPLAlign.seq_of_s)
                          --> Sequence.del_first_char
                          --> fun x -> (Sequence.to_formater x alphabet)::acc)
                    []
                    r.ss
            in
            String.concat "#" seq_lst
    in
    let seq_d = 
        (PXML
            -[Xml.Characters.characters]
                ([Xml.Characters.name] = [`String (Data.code_character mine.code d)])
                { `Fun seq_f }
            --)
    in
    (PXML
        -[Xml.Characters.dlikelihood]
            ([Xml.Characters.cost] = [`Float mine.cost])
            ([Xml.Nodes.min_time] = [str_time t1])
            ([Xml.Nodes.oth_time] = [str_time t2])
            ([attr])
            { `Set (model_d @ [seq_d]) }
        --) :: []


(*---- special functions for likelihood in a dynamic context *)
let min_bl = MlStaticCS.minimum_bl ()
let estimate_time a b =
(*    let is_transition one two =*)
(*        let atran =*)
(*            try match Alphabet.complement one (alph b) with*)
(*                | None   -> false*)
(*                | Some x -> true*)
(*            with | _     -> false*)
(*        and btran =*)
(*            try match Alphabet.complement two (alph a) with*)
(*                | None   -> false*)
(*                | Some x -> true*)
(*            with | _     -> false*)
(*        in*)
(*        atran || btran*)
(*    in*)
    let rec count_bits n =
        if n = 0 then 0 else 1 + (count_bits (n land (n-1)))
    in
    let same = ref 0.0 (* the number of aligned columns *)
    and tran = ref 0.0 (* the number of transitions aligned *)
    and gaps = ref 0.0 (* the number of characters aligned with gaps *)
    and othr = ref 0.0 (* other aligned characters --transversions *)
    and incr r s = r := !r +. s in
    let gchar = Alphabet.get_gap (alph a) in
    let align_s_distance align_a align_b =
        for n = 0 to (Sequence.length align_a) - 1 do
            let achar = Sequence.get align_a n
            and bchar = Sequence.get align_b n in
            (* a_tot and b_tot are the number of set bits in each character, the
             * total of set columns in each become the total minus the aligned
             * states (s_tot). *)
            let a_tot = count_bits achar
            and b_tot = count_bits bchar
            and s_tot = count_bits (achar land bchar land (lnot gchar)) in
            let t_tot = float_of_int ((a_tot + b_tot) - s_tot) in
            (* incr gap if necessary *)
            if ((gchar land achar) > 0) then incr gaps (1.0/.t_tot);
            if ((gchar land achar) > 0) then incr gaps (1.0/.t_tot);
            incr same ((float_of_int s_tot) /. t_tot);
            let m_tot =
                let m = if ((gchar land achar) > 0) then 1 else 0 in
                let m = if ((gchar land bchar) > 0) then m+1 else m in
                let m = t_tot -. (float_of_int (s_tot+m)) in
                if m < 0.0 then 0.0 else m
            in
            incr othr (m_tot /. t_tot)
        done;
        ()
    in
    let () = match a.data,b.data with
        | MPLAlign ar, MPLAlign br ->
            assert( (Array.length ar.ss) = (Array.length br.ss) );
            for i = 0 to (Array.length ar.ss)-1 do
                let aligna,alignb,_ = 
                    Sequence.Align.align_2 
                        (MPLAlign.seq_of_s ar.ss.(i))
                        (MPLAlign.seq_of_s br.ss.(i))
                        (Cost_matrix.Two_D.default) (Matrix.default)
                in
                align_s_distance aligna alignb
            done
        | CMPLAlign ar, CMPLAlign br ->
            assert( (Array.length ar.ss) = (Array.length br.ss) );
            for i = 0 to (Array.length ar.ss)-1 do
                let aligna,alignb,_ = 
                    Sequence.Align.align_2 
                        (CMPLAlign.seq_of_s ar.ss.(i))
                        (CMPLAlign.seq_of_s br.ss.(i))
                        (Cost_matrix.Two_D.default) (Matrix.default)
                in
                align_s_distance aligna alignb
            done
        | Verify ar, Verify br ->
            assert( (Array.length ar.ss) = (Array.length br.ss) );
            for i = 0 to (Array.length ar.ss)-1 do
                let aligna,alignb,_ = 
                    Sequence.Align.align_2 
                        (CMPLAlign.seq_of_s (fst ar.ss.(i)))
                        (CMPLAlign.seq_of_s (fst br.ss.(i)))
                        (Cost_matrix.Two_D.default) (Matrix.default)
                and v_aligna,v_alignb,_ = 
                    Sequence.Align.align_2 
                        (MPLAlign.seq_of_s (snd ar.ss.(i)))
                        (MPLAlign.seq_of_s (snd br.ss.(i)))
                        (Cost_matrix.Two_D.default) (Matrix.default)
                in
                assert( 0 = Sequence.compare aligna v_aligna );
                assert( 0 = Sequence.compare alignb v_alignb );
                align_s_distance aligna alignb
            done
        (* must be consistent *)
        | (MPLAlign _ | Verify _ | CMPLAlign _), _ -> assert false
    in
    let total = (!same +. !tran +. !gaps +. !othr) in
    (* below, similar to mlStatic estimate_time function *)
    let proportion = (!same) /. (total) in
    let p = match (1.0 -. proportion) with
        | x when x < 0.75 -> x
        | x -> 0.70
    in
    let nt2 = ~-. 0.75 *. (log (1.0 -. (p *. 4.0 /. 3.0))) in
    let nt = if nt2 <= min_bl then min_bl else nt2 /. 2.0 in
    if debug then
        Printf.printf
            ("Aligned      : %f\nTransitions  : %f\nGaps         : %f\n"^^
             "Transversions: %f\nTotals       : %f\nProportion   : %f\n\n"^^
             "Distance     : %f\n\n%!")
            !same !tran !gaps !othr total proportion nt2;
    (nt,nt)


let prior a = 
    let m = static_model a in
    let priors = m.MlModel.pi_0 in
    let seq_data = match a.data with
        | MPLAlign a -> Array.map MPLAlign.seq_of_s a.ss
        | CMPLAlign a -> Array.map CMPLAlign.seq_of_s a.ss
        | Verify a -> Array.map (fun (x,_) -> CMPLAlign.seq_of_s x) a.ss
    in
    let max_prior_of_seq (acc:float) sequence =
        Sequence.foldi
            (fun acc pos i -> 
                if pos = 0 then acc (* skip gap opening *)
                else begin
                    let c_pi = 
                        List.fold_left
                            (fun acc i -> max acc priors.{i})
                            0.0 (BitSet.Int.list_of_packed i)
                    in
                    (~-. (log c_pi)) +. acc
                end)
            acc
            sequence
    in
    Array.fold_left
        (fun acc i -> max_prior_of_seq acc i) 0.0 seq_data


let remove_ambiguities dyn =
    let gap = Alphabet.get_gap (alph dyn) in
    let rec remove_low_order_bits x = match x land (x-1) with
        | 0 -> x
        | x -> remove_low_order_bits x
    and remove_low_order_bits_from_seq seq =
        for n = 0 to (Sequence.length seq) - 1 do
            let sstate = Sequence.get seq n in
            let ret = match sstate with
                | x when x = gap -> sstate
                | x -> remove_low_order_bits x
            in
            Sequence.set seq n ret
        done
    in
    if (static_model dyn).MlModel.spec.MlModel.cost_fn = `MPL && false then
        All_sets.IntegerMap.iter
            (fun k v ->
                Array.iter
                    (fun v -> match v with
                        | `DO v | `Last v | `First v ->
                            remove_low_order_bits_from_seq v)
                    v)
            (leaf_sequences dyn);
    dyn

(*---- median functions *)
let median code a b t1 t2 = 
    assert( 0 = MlModel.compare (static_model a) (static_model b) );
    match a.data,b.data with
    | CMPLAlign ar, CMPLAlign br -> 
        let bla,blb = match t1,t2 with
            | Some t1, Some t2 -> t1, t2
            | _,_ -> failwith "branches not specified by caller"
        in
        let cost = ref 0.0 in
        let meds = 
            Array_ops.map_2
                (fun sa sb -> 
                    let mem = CMPLAlign.get_mem sa sb in
                    if debug then CMPLAlign.clear_mem mem;
                    let cst,med = CMPLAlign.median_2_cost sa sb a.model bla blb mem in
                    cost := !cost +. cst;
                    med)
                (ar.ss)
                (br.ss)
        in
        { a with cost = !cost; data = CMPLAlign { ss = meds }; times = bla,blb; }

    | Verify ar, Verify br -> 
        let bla,blb = match t1,t2 with
            | Some t1, Some t2 -> t1, t2
            | _,_ -> failwith "branches not specified by caller"
        in
        let cost = ref 0.0 in
        let meds = 
            Array_ops.map_2
                (fun (csa,osa) (csb,osb) ->
                    let cmem = CMPLAlign.get_mem csa csb and omem = MPLAlign.get_mem osa osb in
                    if debug then begin CMPLAlign.clear_mem cmem;MPLAlign.clear_mem omem end;
                    let ce1,ce2,ccst,cmed = CMPLAlign.gen_all_2 csa csb a.model bla blb cmem in
                    let oe1,oe2,ocst,omed = MPLAlign.gen_all_2 osa osb a.model bla blb omem in
                    let () = 
                        let fm1 = MPLAlign.get_cm a.model bla blb 
                        and fm2 = CMPLAlign.get_cm a.model bla blb in
                        if not (compare_combined_cm fm1 fm2) then begin
                            print_combined_cm fm1;
                            print_newline ();
                            print_combined_cm fm2;
                            assert false
                        end else if not (ccst =. ocst) then begin
                            Printf.printf "%f =/= %f\n" ocst ccst;
                            Printf.printf "O1:%a\nC1:%a\n" (pp_seq a.model) (MPLAlign.seq_of_s oe1) 
                                                         (pp_seq a.model) (CMPLAlign.seq_of_s ce1);
                            Printf.printf "O2:%a\nC2:%a\n" (pp_seq a.model) (MPLAlign.seq_of_s oe2) 
                                                         (pp_seq a.model) (CMPLAlign.seq_of_s ce2);
                            Printf.printf "OM:%a\nCM:%a\n" (pp_seq a.model) (MPLAlign.seq_of_s omed) 
                                                         (pp_seq a.model) (CMPLAlign.seq_of_s cmed);
                            MPLAlign.print_mem omem;
                            assert false
                        end else if 0 != Sequence.compare (CMPLAlign.seq_of_s cmed) (MPLAlign.seq_of_s omed) then begin
                            Printf.printf "%f =/= %f\n" ocst ccst;
                            Printf.printf "O1:%a\nC1:%a\n" (pp_seq a.model) (MPLAlign.seq_of_s oe1) 
                                                         (pp_seq a.model) (CMPLAlign.seq_of_s ce1);
                            Printf.printf "O2:%a\nC2:%a\n" (pp_seq a.model) (MPLAlign.seq_of_s oe2) 
                                                         (pp_seq a.model) (CMPLAlign.seq_of_s ce2);
                            Printf.printf "OM:%a\nCM:%a\n" (pp_seq a.model) (MPLAlign.seq_of_s omed) 
                                                         (pp_seq a.model) (CMPLAlign.seq_of_s cmed);
                            MPLAlign.print_mem omem;
                            assert false
                        end else if 0 != Sequence.compare (CMPLAlign.seq_of_s ce1) (MPLAlign.seq_of_s oe1) then begin
                            Printf.printf "%f =/= %f\n" ocst ccst;
                            Printf.printf "O1:%a\nC1:%a\n" (pp_seq a.model) (MPLAlign.seq_of_s oe1) 
                                                         (pp_seq a.model) (CMPLAlign.seq_of_s ce1);
                            Printf.printf "O2:%a\nC2:%a\n" (pp_seq a.model) (MPLAlign.seq_of_s oe2) 
                                                         (pp_seq a.model) (CMPLAlign.seq_of_s ce2);
                            Printf.printf "OM:%a\nCM:%a\n" (pp_seq a.model) (MPLAlign.seq_of_s omed) 
                                                         (pp_seq a.model) (CMPLAlign.seq_of_s cmed);
                            MPLAlign.print_mem omem;
                            assert false
                        end else if 0 != Sequence.compare (CMPLAlign.seq_of_s ce2) (MPLAlign.seq_of_s oe2) then begin
                            Printf.printf "%f =/= %f\n" ocst ccst;
                            Printf.printf "O1:%a\nC1:%a\n" (pp_seq a.model) (MPLAlign.seq_of_s oe1) 
                                                         (pp_seq a.model) (CMPLAlign.seq_of_s ce1);
                            Printf.printf "O2:%a\nC2:%a\n" (pp_seq a.model) (MPLAlign.seq_of_s oe2) 
                                                         (pp_seq a.model) (CMPLAlign.seq_of_s ce2);
                            Printf.printf "OM:%a\nCM:%a\n" (pp_seq a.model) (MPLAlign.seq_of_s omed) 
                                                         (pp_seq a.model) (CMPLAlign.seq_of_s cmed);
                            MPLAlign.print_mem omem;
                            assert false
                        end
                    in
                    cost := !cost +. ccst;
                    cmed,omed)
                (ar.ss)
                (br.ss)
        in
        { a with cost = !cost; data = Verify { ss = meds }; times = bla,blb; }

    | MPLAlign ar, MPLAlign br -> 
        let bla,blb = match t1,t2 with
            | Some t1, Some t2 -> t1, t2
            | _,_ -> failwith "branches not specified by caller"
        in
        let cost = ref 0.0 in
        let meds = 
            Array_ops.map_2
                (fun sa sb -> 
                    let mem = MPLAlign.get_mem sa sb in
                    if debug then MPLAlign.clear_mem mem;
                    let cst,med = MPLAlign.median_2_cost sa sb a.model bla blb mem in
                    cost := !cost +. cst;
                    med)
                (ar.ss)
                (br.ss)
        in
        { a with cost = !cost; data = MPLAlign { ss = meds }; times = bla,blb; }
    | (MPLAlign _ | Verify _ | CMPLAlign _), _ -> assert false


let readjust c1 c2 mine t1 t2 = 
    let internal_loop (t:float) : t * float =
        let d = median (code mine) c1 c2 (Some t) (Some t2) in
        (d,total_cost d)
    in
    let nt,(nmine,ncost) = 
        Numerical.brents_method internal_loop (t1,(mine,total_cost mine))
    in
    if debug_est then
        Printf.printf "Optimized Branch: %f -(%f->%f)-> %f\n%!" 
                      (total_cost mine) t1 nt ncost;
    let modified = not ((total_cost mine) =. ncost) in
    modified,total_cost mine,ncost,(nt,t2),nmine

let readjust3 mine c1 c2 par t1 t2 t3 =
    assert( 0 = MlModel.compare (static_model mine) (static_model c2) );
    assert( 0 = MlModel.compare (static_model mine) (static_model c1) );
    assert( 0 = MlModel.compare (static_model mine) (static_model par));
    let node = match c1.data, c2.data, par.data with
        | CMPLAlign c1, CMPLAlign c2, CMPLAlign p ->
            let ns = 
                Array_ops.map_3
                    (fun x y z ->
                        let c,n,_ = CMPLAlign.readjust x y z mine.model t1 t2 t3 in n)
                    c1.ss c2.ss p.ss
            in
            CMPLAlign { ss = ns; }
        | MPLAlign c1, MPLAlign c2, MPLAlign p ->
            let ns = 
                Array_ops.map_3
                    (fun x y z ->
                        let c,n,_ = MPLAlign.readjust x y z mine.model t1 t2 t3 in n)
                    c1.ss c2.ss p.ss
            in
            MPLAlign { ss = ns; }
        | Verify c1, Verify c2, Verify p                          -> assert false
        | (MPLAlign _ | CMPLAlign _ | Verify _), _, _ -> assert false
    in
    {mine with data = node;}



let readjust3_opt mine c1 c2 par t1 t2 t3 =
    let single_result =
        (function
            | [| t1; t2; t3 |] ->
                let cst = ref 0.0 in
                let nds = match c1.data, c2.data, par.data with
                    | CMPLAlign c1, CMPLAlign c2, CMPLAlign p ->
                        let ns = 
                            Array_ops.map_3
                                (fun x y z ->
                                    let c,n,_ = CMPLAlign.readjust x y z mine.model t1 t2 t3 in
                                    cst := !cst +. c;
                                    n)
                                c1.ss c2.ss p.ss
                        in
                        CMPLAlign { ss = ns; }

                    | Verify c1, Verify c2, Verify p -> (* TODO *) assert false

                    | MPLAlign c1, MPLAlign c2, MPLAlign p ->
                        let ns = 
                            Array_ops.map_3
                                (fun x y z ->
                                    let c,n,_ = MPLAlign.readjust x y z mine.model t1 t2 t3 in
                                    cst := !cst +. c;
                                    n)
                                c1.ss c2.ss p.ss
                        in
                        MPLAlign { ss = ns; }
                    | (MPLAlign _ | CMPLAlign _ | Verify _), _, _ -> assert false
                in
                (nds,!cst)
            | _ -> assert false)
    in
    let (times,(node,score)),pscore =
        let init_vec = [| t1;t2;t3 |] in
        let (_,pscore) as init_res = single_result init_vec in
        Numerical.run_method
            (Numerical.default_numerical_optimization_strategy !Methods.opt_mode 3)
            single_result (init_vec,init_res),
        pscore
    in
    true, pscore, score, (times.(0),times.(1),times.(2)), {mine with data = node;}

let median_i code a b (t1:float) (t2:float) : t * float * float =
    let m = median code a b (Some t1) (Some t2) in
    let _,_,_,(t1,t2),m = readjust a b m t1 t2 in
    m,t1,t2

and median_3 p n c1 c2 = assert false (** TODO **)
(*    assert( 0 = MlModel.compare n.model c1.model);*)
(*    assert( 0 = MlModel.compare n.model c2.model);*)
(*    assert( 0 = MlModel.compare n.model p.model );*)
(*    { *)
(*        seq = SeqCS.median_3 p.seq n.seq c1.seq c2.seq;*)
(*      model = n.model;*)
(*         ia = None; *)
(*    }*)

(*---- distance functions; all for Sequence characters. *)
let distance missing_distance a b blen =
    let dist1,dist2 = match blen with
        | None ->
            let d1,d2 = estimate_time a b in
            (Some d1), (Some d2)
        | Some _ -> blen, Some 0.0
    in
    total_cost (median (-1) a b dist1 dist2)

let dist_2 delta n a b =
    let dist1,dist2 = estimate_time a b in
    let tmp = median (-1) a b (Some dist1) (Some dist2) in
    let dist1,dist2 = estimate_time n tmp in
    total_cost (median (-1) a b (Some dist1) (Some dist2))
            
let tabu_distance a = match a.data with
    | MPLAlign _    -> failwith_todo "tabu_distance (MPLAlign)"
    | CMPLAlign _   -> failwith_todo "tabu_distance (CMPLAlign)" 
    | Verify _      -> failwith_todo "tabu_distance (Verify)" 

(*---- filter functions *)
let array_filter f a b = 
    let len = Array.length a in
    assert( len = Array.length b );
    let rec array_filter_pair acc1 acc2 i =
        if i = len then (List.rev acc1, List.rev acc2)
        else if not (f a.(i)) then array_filter_pair acc1 acc2 (i+1)
        else array_filter_pair (a.(i)::acc1) (b.(i)::acc2) (i+1)
    in
    let one,two = array_filter_pair [] [] 0 in
    Array.of_list one, Array.of_list two

let f_codes s c = match s.data with
    | MPLAlign r ->
        let codes,seqs = 
            array_filter (fun x -> All_sets.Integers.mem x c) s.codes r.ss
        in
        { s with codes = codes; data = MPLAlign { ss = seqs; } }
    | Verify r ->
        let codes,seqs = 
            array_filter (fun x -> All_sets.Integers.mem x c) s.codes r.ss
        in
        { s with codes = codes; data = Verify { ss = seqs; } }
    | CMPLAlign r ->
        let codes,seqs = 
            array_filter (fun x -> All_sets.Integers.mem x c) s.codes r.ss
        in
        { s with codes = codes; data = CMPLAlign { ss = seqs; } }


let f_codes_comp s c = match s.data with
    | MPLAlign r ->
        let codes,seqs = 
            array_filter (fun x -> not (All_sets.Integers.mem x c)) s.codes r.ss
        in
        { s with codes = codes; data = MPLAlign { ss = seqs; } }
    | Verify r ->
        let codes,seqs = 
            array_filter (fun x -> not (All_sets.Integers.mem x c)) s.codes r.ss
        in
        { s with codes = codes; data = Verify { ss = seqs; } }
    | CMPLAlign r ->
        let codes,seqs = 
            array_filter (fun x -> not (All_sets.Integers.mem x c)) s.codes r.ss
        in
        { s with codes = codes; data = CMPLAlign { ss = seqs; } }

(*---- make an initial leaf node; still requires IA. *)
let make a s m = 
    let r = match m.MlModel.spec.MlModel.cost_fn with
        | `MPL ->
            if verify then begin
                let data =
                    Array.map
                        (fun s -> CMPLAlign.s_of_seq s,MPLAlign.s_of_seq s)
                        (s_of_seq s)
                in
                Verify { ss = data; }
            end else if pure_ocaml then begin
                let data = Array.map MPLAlign.s_of_seq (s_of_seq s) in
                MPLAlign { ss = data; }
            end else begin (* DEFAULT *)
                let data = Array.map CMPLAlign.s_of_seq (s_of_seq s) in
                CMPLAlign { ss = data; }
            end
        | `MAL ->
            failwith "Dynamic Maximum Average Likelihood is not implemented to diagnose trees."
    in
    let data = 
        {    data = r;
            model = { static = m; alph = a; };
             cost = 0.0; (* fill this in with call to prior *)
             code = s.SeqCS.code;
            times = (0.0,0.0);
            codes = s.SeqCS.codes; }
    in
    { data with cost = prior data; }

let to_single parent mine t = 
    let pcost = total_cost mine in
    match parent.data,mine.data with    
    | CMPLAlign ps, CMPLAlign ms -> 
        let score = ref 0.0 in
        let n_data = 
            Array_ops.map_2 
                (fun p m ->
                    let mem = CMPLAlign.get_mem p m in
                    let r,s = CMPLAlign.closest ~p ~m mine.model t mem in
                    score := s +. !score; r)
                ps.ss
                ms.ss
        in
        if debug_est then
            Printf.printf "MlDynamicCS.to_single: t:%f / c:%f\n\n%!" t !score;
        pcost, !score, { mine with data = CMPLAlign { ss = n_data }; }

    | Verify ps, Verify ms -> 
        let score = ref 0.0 in
        let n_data = 
            Array_ops.map_2 
                (fun (cp,op) (cm,om) ->
                    let omem = MPLAlign.get_mem op om 
                    and cmem = CMPLAlign.get_mem cp cm in
                    let c_r,cs = CMPLAlign.closest ~p:cp ~m:cm mine.model t cmem
                    and o_r,os = MPLAlign.closest ~p:op ~m:om mine.model t omem in
                    if not ((cs =. os) || ( 0 = Sequence.compare (CMPLAlign.seq_of_s c_r) (MPLAlign.seq_of_s o_r))) then begin
                        Printf.printf "C:\nP:%a\nM:%a\n%f\n\n%!" (pp_seq mine.model) (CMPLAlign.seq_of_s cp)
                                                        (pp_seq mine.model) (CMPLAlign.seq_of_s c_r) cs;
                        Printf.printf "O:\nP:%a\nM:%a\n%f\n\n%!" (pp_seq mine.model) (MPLAlign.seq_of_s op)
                                                        (pp_seq mine.model) (MPLAlign.seq_of_s o_r) os;

                        assert false;
                    end;
                    score := cs +. !score; c_r,o_r)
                ps.ss
                ms.ss
        in
        if debug_est then
            Printf.printf "MlDynamicCS.to_single: t:%f / c:%f\n\n%!" t !score;
        pcost, !score, { mine with data = Verify { ss = n_data }; }

    | MPLAlign ps, MPLAlign ms -> 
        let score = ref 0.0 in
        let n_data = 
            Array_ops.map_2 
                (fun p m ->
                    let mem = MPLAlign.get_mem p m in
                    let r,s = MPLAlign.closest ~p ~m mine.model t mem in
                    score := s +. !score; r)
                ps.ss
                ms.ss
        in
        if debug_est then
            Printf.printf "MlDynamicCS.to_single: t:%f / c:%f\n\n%!" t !score;
        pcost, !score, { mine with data = MPLAlign { ss = n_data }; }
    (* although weak, this is the only solution *)
    | (MPLAlign _ | Verify _ | CMPLAlign _), _ -> assert false

ELSE

(* empty; required functions when likelihood is not enabled *)

type r = unit

type t = { model  : dyn_model;
            data  : r;
            codes : int array;
            code  : int;
            times : float * float; (* branch lengths of composition *)
            cost  : float; }

let alph _ = failwith MlStaticCS.likelihood_error
let total_cost _ = failwith MlStaticCS.likelihood_error
let get_cm _ = failwith MlStaticCS.likelihood_error
let model _ = failwith MlStaticCS.likelihood_error
let code _ = failwith MlStaticCS.likelihood_error
let get_codes _ = failwith MlStaticCS.likelihood_error
let make _ _ _ = failwith MlStaticCS.likelihood_error
let estimate_time _ _ = failwith MlStaticCS.likelihood_error
let median _ _ _ _ _ = failwith MlStaticCS.likelihood_error
let median_i _ _ _ _ _ = failwith MlStaticCS.likelihood_error
let median_3 _ _ _ _ = failwith MlStaticCS.likelihood_error
let readjust _ _ _ _ _ = failwith MlStaticCS.likelihood_error
let readjust3 _ _ _ _ _ _ _ = failwith MlStaticCS.likelihood_error
let distance _ _ _ = failwith MlStaticCS.likelihood_error
let to_string _ = failwith MlStaticCS.likelihood_error
let compare _ _ = failwith MlStaticCS.likelihood_error
let name_string _ = failwith MlStaticCS.likelihood_error
let dist_2 _ _ _ _ = failwith MlStaticCS.likelihood_error
let f_codes _ _ = failwith MlStaticCS.likelihood_error
let f_codes_comp _ _ = failwith MlStaticCS.likelihood_error
let tabu_distance _ = failwith MlStaticCS.likelihood_error
let cardinal _ = failwith MlStaticCS.likelihood_error
let encoding _ _ = failwith MlStaticCS.likelihood_error
let leaf_sequences _ = failwith MlStaticCS.likelihood_error
let to_formatter _ _ _ _ _ = failwith MlStaticCS.likelihood_error
let prior _ = failwith MlStaticCS.likelihood_error
let to_single _ _ _ = failwith MlStaticCS.likelihood_error
let readjust3_opt _ _ _ _ _ _ _ = failwith MlStaticCS.likelihood_error

ENDIF
