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

let () = SadmanOutput.register "MlDynamicCS" "$Id$"

(*---- non-external helper functions/settings *)
open Numerical.FPInfix
open FloatSequence
let debug     = false
let debug_est = false
let verify    = false
let (-->) a b = b a
let failwith_todo func = failwith ("Your lazy developer needs to write: MlDynamicCS."^func)
let debug_printf = 
    if debug then (fun format -> print_string format)
             else (fun _ -> ())

IFDEF USE_LIKELIHOOD THEN

(*---- type for the module *)
type 'a align = { ss : 'a array; }

(* Define how we should align and create the medians; 
 *   FPAlign     - Same as the above, but using new routines to align by
 *                 floating point values.
 *   MPLAlign    - Two transformation matrices are created for each sequence,
 *                 these measure the probability of each transition, and when
 *                 find the minimum value (X->I)(Y->I) for the assignment I. *)
type r = | FPAlign  of FloatAlign.s align
         | MPLAlign of MPLAlign.s align

type t = { model  : dyn_model;
            data  : r;
            codes : int array;
            code  : int;
            cost  : float; }

(*---- basic function set for gathering data of the module *)
let get_cm t = match t.data with
    | MPLAlign _
    | FPAlign _     -> assert false

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
    match a.data,b.data with
    | FPAlign ad, FPAlign bd -> 
        (0 = MlModel.compare a.model.static
                             b.model.static) &&
        (compare_array2 FloatAlign.compare ad.ss bd.ss)
    | MPLAlign ad, MPLAlign bd -> 
        (0 = MlModel.compare a.model.static
                             b.model.static) &&
        (compare_array2 MPLAlign.compare ad.ss bd.ss)
    | (MPLAlign _ | FPAlign _) , _ -> false

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
    | FPAlign r -> 
        Array.fold_left
            (fun acc x -> 
                Sequence.encoding (e) (FloatAlign.seq_of_s x))
            (0.0)
            (r.ss)

let name_string t = match t.data with
    | MPLAlign _ -> "Maximum Parsimonious Dynamic Likelihood"
    | FPAlign _  -> "Floating Point Ancesteral Dynamic Likelihood"

let total_cost t = t.cost

let model t = t.model

let static_model t = t.model.static

let s_of_seq seq : Sequence.s array =
    seq --> SeqCS.get_sequences --> Array.to_list --> Array.concat

let leaf_sequences t = match t.data with
    | FPAlign r ->
        Array_ops.fold_right_2
            (fun acc x y ->
                let y = [| `DO (FloatAlign.seq_of_s y) |] in
                All_sets.IntegerMap.add x y acc)
            (All_sets.IntegerMap.empty)
            (t.codes)
            (r.ss)
    | MPLAlign r ->
        Array_ops.fold_right_2
            (fun acc x y ->
                let y = [| `DO (MPLAlign.seq_of_s y) |] in
                All_sets.IntegerMap.add x y acc)
            (All_sets.IntegerMap.empty)
            (t.codes)
            (r.ss)

(*---- to formatter, and printing functions *)
let to_formatter attr mine par_opt (t1,t2) d : Xml.xml Sexpr.t list = 
    let str_time = function | Some x -> `Float x | None -> `String "None" in
    match mine.data with
    | FPAlign r     ->
        let model_d = MlModel.to_formatter (static_model mine) in
        let alphabet= alph mine in
        let seq_f () =
            let seq_lst = 
                Array.fold_left
                    (fun acc s ->
                        s --> FloatAlign.seq_of_s
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

    | MPLAlign r    ->
        let model_d = MlModel.to_formatter (static_model mine) in
        let alphabet= alph mine in
        let seq_f () =
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
        | FPAlign ar, FPAlign br ->
            assert( (Array.length ar.ss) = (Array.length br.ss) );
            for i = 0 to (Array.length ar.ss)-1 do
                let aligna,alignb,_ = 
                    Sequence.Align.align_2 
                        (FloatAlign.seq_of_s ar.ss.(i))
                        (FloatAlign.seq_of_s br.ss.(i))
                        (Cost_matrix.Two_D.default) (Matrix.default)
                in
                align_s_distance aligna alignb
            done
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
        (* must be consistent *)
        | (FPAlign _ | MPLAlign _ ), _ -> assert false
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

(* Verify Cost; only used for FPAlign *)
let verify_fp_cost cst1 sa sb model bla blb mem = 
    if verify then begin
        let cst2 = FloatAlign.full_cost_2 cst1 sa sb model bla blb mem in
        if cst1 =. cst2 then true
        else if not debug then false 
        else begin 
            Printf.printf "%f --- %f\nFull Align:\n" cst1 cst2;
            FloatAlign.print_cm model bla;
            FloatAlign.print_cm model blb;
            FloatAlign.print_mem mem;
            FloatAlign.clear_mem mem;
            Printf.printf "Ukkonen Align:\n";
            ignore (FloatAlign.cost_2 sa sb model bla blb mem);
            FloatAlign.print_mem mem;
            false
        end 
    end else begin
        true
    end

let verify_mpl_cost cst1 sa sb model bla blb mem = 
    if verify then begin
        let cst2 = MPLAlign.full_cost_2 cst1 sa sb model bla blb mem in
        if cst1 =. cst2 then true
        else if not debug then false 
        else begin 
            Printf.printf "%f --- %f\nFull Align:\n" cst1 cst2;
            MPLAlign.print_cm model bla;
            MPLAlign.print_cm model blb;
            MPLAlign.print_mem mem;
            MPLAlign.clear_mem mem;
            Printf.printf "Ukkonen Align:\n";
            ignore (MPLAlign.cost_2 sa sb model bla blb mem);
            MPLAlign.print_mem mem;
            false
        end 
    end else begin
        true
    end

let prior a = 
    let m = static_model a in
    let priors = m.MlModel.pi_0 in
    match a.data with
    | FPAlign a ->
        let avg_prior_of_seq (acc:float) sequence =
            Sequence.foldi
                (fun acc pos i -> 
                    if pos = 0 then acc (* skip gap opening *)
                    else begin
                        let c_pi,len = 
                            List.fold_left
                                (fun (acc,n) i -> acc +. priors.{i},n+1)
                                (0.0,0)
                                (BitSet.Int.list_of_packed i)
                        in
                        (~-. ((log c_pi) -. (log (float len)))) +. acc
                    end)
                acc
                sequence
        in
        Array.fold_left
            (fun acc i -> avg_prior_of_seq acc (FloatAlign.seq_of_s i))
            (0.0) 
            (a.ss)
    | MPLAlign a ->
        let max_prior_of_seq (acc:float) sequence =
            Sequence.foldi
                (fun acc pos i -> 
                    if pos = 0 then acc (* skip gap opening *)
                    else begin
                        let c_pi = 
                            List.fold_left
                                (fun acc i -> max acc priors.{i})
                                (0.0)
                                (BitSet.Int.list_of_packed i)
                        in
                        (~-. (log c_pi)) +. acc
                    end)
                (acc)
                (sequence)
        in
        Array.fold_left
            (fun acc i -> max_prior_of_seq acc (MPLAlign.seq_of_s i))
            (0.0) 
            (a.ss)

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
                    if debug then assert( verify_mpl_cost cst sa sb a.model bla blb mem );
                    cost := !cost +. cst;
                    med)
                (ar.ss)
                (br.ss)
        in
        { a with cost = !cost; data = MPLAlign { ss = meds } }
    | FPAlign ar, FPAlign br -> 
        let bla,blb = match t1,t2 with
            | Some t1, Some t2 -> t1, t2
            | _,_ -> failwith "branches not specified by caller"
        in
        let cost = ref 0.0 in
        let meds = 
            Array_ops.map_2
                (fun sa sb -> 
                    let mem = FloatAlign.get_mem sa sb in
                    if debug then FloatAlign.clear_mem mem;
                    let cst,med = FloatAlign.median_2_cost sa sb a.model bla blb mem in
                    if debug then assert( verify_fp_cost cst sa sb a.model bla blb mem );
                    cost := !cost +. cst;
                    med)
                (ar.ss)
                (br.ss)
        in
        { a with cost = !cost; data = FPAlign { ss = meds; } } 
    | (FPAlign _ | MPLAlign _ ), _ -> assert false

let readjust c1 c2 mine t1 t2 = 
    let internal_loop (t:float) : t * float =
        let d = median (code mine) c1 c2 (Some t) (Some t2) in
        (d,total_cost d)
    in
    let nt,(nmine,ncost) = 
        Numerical.brents_method (t1,(mine,total_cost mine)) internal_loop
    in
    if debug_est then
        Printf.printf "Optimized Branch: %f -(%f->%f)-> %f\n%!" 
                      (total_cost mine) t1 nt ncost;
    let modified = not ((total_cost mine) =. ncost) in
    modified,total_cost mine,ncost,(nt,t2),nmine

let readjust3 mine c1 c2 par t1 t2 t3 =
    let single_result =
        (function
            | [| t1; t2; t3 |] ->
                let cst = ref 0.0 in
                let nds = match c1.data, c2.data, par.data with
                    | FPAlign c1, FPAlign c2, FPAlign p ->
                        let ns = 
                            Array_ops.map_3
                                (fun x y z ->
                                    let c,n,_ = FloatAlign.readjust x y z mine.model t1 t2 t3 in 
                                    cst := !cst +. c;
                                    n)
                                c1.ss c2.ss p.ss
                        in
                        FPAlign { ss = ns; }

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
                    | (MPLAlign _ | FPAlign _), _, _ -> assert false
                in
                (nds,!cst)
            | _ -> assert false)
    in
    let (times,(node,score)),pscore =
        let init_vec = [| t1;t2;t3 |] in
        let (_,pscore) as init_res = single_result init_vec in
        Numerical.bfgs_method
            ~epsilon:1e-5 single_result init_vec init_res, pscore
    in
    true, pscore, score, (times.(0),times.(1),times.(2)), {mine with data = node;}

let median_i code a b (t1:float) (t2:float) : t * float * float =
    match a.data, b.data with
    | FPAlign _, FPAlign _ 
    | MPLAlign _, MPLAlign _ -> 
        median code a b (Some t1) (Some t2),t1,t2
    | (FPAlign _ | MPLAlign _ ), _ -> assert false

and median_3 p n c1 c2 = n
(*    assert( 0 = MlModel.compare n.model c1.model);*)
(*    assert( 0 = MlModel.compare n.model c2.model);*)
(*    assert( 0 = MlModel.compare n.model p.model );*)
(*    { *)
(*        seq = SeqCS.median_3 p.seq n.seq c1.seq c2.seq;*)
(*      model = n.model;*)
(*         ia = None; *)
(*    }*)

(*---- distance functions; all for Sequence characters. *)
let distance missing_distance a b blen = match a.data, b.data with
    | MPLAlign _ , MPLAlign _
    | FPAlign _ , FPAlign _ ->
        let dist1,dist2 = match blen with
            | None ->
                let d1,d2 = estimate_time a b in
                (Some d1), (Some d2)
            | Some _ -> blen, Some 0.0
        in
        total_cost (median (-1) a b dist1 dist2)
    | (FPAlign _ | MPLAlign _ ), _ -> assert false

let dist_2 delta n a b = match a.data,b.data,n.data with
    | MPLAlign _, MPLAlign _, MPLAlign _
    | FPAlign _, FPAlign _, FPAlign _ -> 
        let dist1,dist2 = estimate_time a b in
        let tmp = median (-1) a b (Some dist1) (Some dist2) in
        let dist1,dist2 = estimate_time n tmp in
        total_cost (median (-1) a b (Some dist1) (Some dist2))
    | (FPAlign _ | MPLAlign _ ), _, _ -> assert false
            
let tabu_distance a = match a.data with
    | FPAlign a -> failwith_todo "tabu_distance (FPAlign)"
    | MPLAlign a -> failwith_todo "tabu_distance (FPAlign)"

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
    | FPAlign r ->
        let codes,seqs = 
            array_filter (fun x -> All_sets.Integers.mem x c) s.codes r.ss
        in
        { s with codes = codes; data = FPAlign { ss = seqs; } }

let f_codes_comp s c = match s.data with
    | MPLAlign r ->
        let codes,seqs = 
            array_filter (fun x -> not (All_sets.Integers.mem x c)) s.codes r.ss
        in
        { s with codes = codes; data = MPLAlign { ss = seqs; } }
    | FPAlign r ->
        let codes,seqs = 
            array_filter (fun x -> not (All_sets.Integers.mem x c)) s.codes r.ss
        in
        { s with codes = codes; data = FPAlign { ss = seqs; } }


(*---- make an initial leaf node; still requires IA. *)
let make a s m = 
    let r = match m.MlModel.spec.MlModel.cost_fn with
        | `FLK ->
            let data = Array.map FloatAlign.s_of_seq (s_of_seq s) in
            FPAlign { ss = data; }
        | `MPL ->
            let data = Array.map MPLAlign.s_of_seq (s_of_seq s) in
            MPLAlign { ss = data; }
        | `MAL ->
            failwith "Dynamic Maximum Average Likelihood is not implemented to diagnose trees."
    in
    let data = 
        {    data = r;
            model = { static = m; alph = a; };
             cost = 0.0; (* fill this in with call to prior *)
             code = s.SeqCS.code;
            codes = s.SeqCS.codes; }
    in
    { data with cost = prior data; }

let to_single parent mine t = 
    let pcost = total_cost mine in
    match parent.data,mine.data with    
    | FPAlign ps, FPAlign ms ->
        let score = ref 0.0 in
        let n_data = 
            Array_ops.map_2
                (fun p m -> 
                    let mem = FloatAlign.get_mem p m in
                    let r,s = FloatAlign.closest ~p ~m mine.model t mem in
                    score := s +. !score; r)
                ps.ss
                ms.ss
        in
        pcost, !score, {mine with data = FPAlign {ss = n_data}; }
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
    | (FPAlign _ | MPLAlign _ ), _ -> assert false

ELSE

(* empty; required functions when likelihood is not enabled *)

type 'a align = { ss : 'a array; }

type r = | FPAlign      of FloatAlign.s align
         | MPLAlign     of MPLAlign.s align

type t = { model  : dyn_model;
            data  : r;
            codes : int array;
            code  : int;
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

ENDIF

