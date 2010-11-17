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


IFDEF USE_LIKELIHOOD THEN
(*---- type for the module *)

(* The process of constructing this record is as follows,
    -- Downpass/Uppass to determine seq alignments
    -- Implied Alignment of sequence data
    -- Downpass/Uppass of Implied Alignment data using likelihood criteria
    -- Combine two versions of data to get likelihood score from IA. *)
type t_integerized =
    {     seq : SeqCS.t; 
       imodel : MlModel.model;
           ia : MlStaticCS.t option; }

(* We ignore the implied alignment and only discuss the cost on the alignment
 * via the product of the entire sequence states, across all edges *)
type t_fpalign = 
    {     ss : FloatSequence.FloatAlign.s array;
      fmodel : MlModel.model;
       codes : int array;
        code : int;
        cost : float; }

(* Define how we should align and create the medians; 
 *   Integerized - the original method; uses Seq alignment methods and creates
 *                 an integerized cost matrix from the floating-point values.
 *   FPAlign     - Same as the above, but using new routines to align by
 *                 floating point values.
 *   FixedStates - ....
 *   Likelihood  - ....
*)
type t = | Integerized of t_integerized
         | FPAlign     of t_fpalign
(*         | IFPAlign    of t_integerized * t_fpalign*)
(*         | FixedStates*)
(*         | DynamicMPL*)
(*         | DynamicMAL*)

(*---- non-external helper functions/settings *)
let debug     = false
let debug_est = false
let verify    = false
let (-->) a b = b a
let (=.) ?(epsilon=10e-6) a b = (abs_float (a-.b)) < epsilon
let failwith_todo func = failwith ("Your lazy developer needs to write: MlDynamicCS."^func)
let debug_printf = 
    if debug then (fun format -> print_string format)
             else (fun _ -> ())

(*---- basic function set for gathering data of the module *)
let get_cm t = match t with
    | Integerized t -> t.seq.SeqCS.heuristic.SeqCS.c2
    | FPAlign t     -> failwith_todo "get_cm (FPAlign)"

let to_string t = failwith_todo "to_string"

let alph t = match t with
    | FPAlign t     -> t.fmodel.MlModel.alph
    | Integerized t -> t.imodel.MlModel.alph

let code t = match t with
    | FPAlign t     -> t.code
    | Integerized t -> t.seq.SeqCS.code

let cardinal t = match t with
    | FPAlign t     -> Array.length t.codes
    | Integerized t -> SeqCS.cardinal t.seq

let encoding e t = match t with
    | FPAlign t -> 
        Array.fold_left
            (fun acc x -> 
                Sequence.encoding (e) (FloatSequence.FloatAlign.seq_of_s x))
            (0.0)
            (t.ss)
    | Integerized t -> SeqCS.encoding e t.seq

let name_string t = match t with
    | FPAlign _     -> "Floating Point Ancesteral Dynamic Likelihood"
    | Integerized t -> 
        begin match t.ia with
            | Some _ -> "Integerized Ancesteral Dynamic Likelihood (with ia)"
            | None   -> "Integerized Ancesteral Dynamic Likelihood (w/out ia)"
        end

let total_cost t = match t with
    | FPAlign t     -> t.cost
    | Integerized t -> 
        begin match t.ia with
            | Some x -> MlStaticCS.median_cost x
            | None   -> t.seq.SeqCS.total_cost
        end

let get_codes t = match t with
    | Integerized t -> t.seq.SeqCS.codes
    | FPAlign t     -> t.codes

let combine t (b:MlStaticCS.t) = match t with
    | Integerized t -> Integerized {t with ia = Some b}
    | FPAlign c     ->
        debug_printf "MlDynamicCS.combine; Ignoring Combine of with an IA";
        t

let model t = match t with
    | Integerized t -> t.imodel
    | FPAlign t     -> t.fmodel

let s_of_seq seq : FloatSequence.FloatAlign.s array = 
    seq --> SeqCS.get_sequences
        --> Array.to_list
        --> Array.concat
        --> Array.map FloatSequence.FloatAlign.s_of_seq

let leaf_sequences t = match t with
    | Integerized t ->
        let a = t.seq in
        let map = ref All_sets.IntegerMap.empty in
        for i = (SeqCS.cardinal a) - 1 downto 0 do
            map := All_sets.IntegerMap.add a.SeqCS.codes.(i)
                (match a.SeqCS.characters.(i) with
                | SeqCS.Partitioned x ->
                    Array.map (function 
                        | SeqCS.PartitionedDOS.Last x ->
                                `Last x.SeqCS.DOS.sequence
                        | SeqCS.PartitionedDOS.DO x  ->
                                `DO x.SeqCS.DOS.sequence
                        | SeqCS.PartitionedDOS.First x ->
                                `First x.SeqCS.DOS.sequence) x
                | SeqCS.Heuristic_Selection x -> [|`DO x.SeqCS.DOS.sequence|]
                | SeqCS.Relaxed_Lifted (t, x) -> 
                    let p = SeqCS.RL.find_smallest x in
                    [|`DO t.SeqCS.RL.sequence_table.(p)|])
                !map
        done;
        !map
    | FPAlign t ->
        Array_ops.fold_right_2
            (fun acc x y ->
                let y = [| `DO (FloatSequence.FloatAlign.seq_of_s y) |] in
                All_sets.IntegerMap.add x y acc)
            (All_sets.IntegerMap.empty)
            (t.codes)
            (t.ss)

(*---- to formatter, and printing functions *)
let to_formatter attr t par_opt bl d = match t with
    | FPAlign t     -> failwith_todo "to_formatter"
    | Integerized t -> 
        begin match t.ia with
            | Some t -> MlStaticCS.to_formatter attr t bl d
            | None   -> failwith "MlDynamicCS.to_formatter does not have an implied alignment"
        end

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
    let gchar = Alphabet.get_gap (Alphabet.explote (alph a) 1 0) in
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
    let () = match a,b with
        | Integerized at,Integerized bt -> 
            let aligned_triples = SeqCS.align_2 at.seq bt.seq in
            assert( (alph a) = (alph b) );
            Array.iter 
                (fun triple_array ->
                    Array.iter 
                        (fun (align_a,align_b,_) -> align_s_distance align_a align_b)
                        (triple_array) )
                aligned_triples;
        | FPAlign at, FPAlign bt ->
            assert( (Array.length at.ss) = (Array.length bt.ss) );
            for i = 0 to (Array.length at.ss)-1 do
                let aligna,alignb,_ = 
                    Sequence.Align.align_2 
                        (FloatSequence.FloatAlign.seq_of_s at.ss.(i))
                        (FloatSequence.FloatAlign.seq_of_s bt.ss.(i))
                        (Cost_matrix.Two_D.default) (Matrix.default)
                in
                align_s_distance aligna alignb
            done
        (* must be consistent *)
        | _ , _ -> failwith "MlDynamicCS.estimate_time; inconsistent procedures"
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
let verify_cost cst1 sa sb model bl mem = 
    if verify then begin
        let cst2 = FloatSequence.FloatAlign.verify_cost_2 cst1 sa sb model bl mem in
        if cst1 =. cst2 then true
        else if not debug then false 
        else begin 
            Printf.printf "%f --- %f\nFull Align:\n" cst1 cst2;
            FloatSequence.FloatAlign.print_cm model bl;
            FloatSequence.FloatAlign.print_mem mem;
            FloatSequence.FloatAlign.clear_mem mem;
            Printf.printf "Ukkonen Align:\n";
            ignore (FloatSequence.FloatAlign.cost_2 ~debug:true sa sb model bl mem);
            FloatSequence.FloatAlign.print_mem mem;
            false
        end 
    end else begin
        true
    end


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
    if (model dyn).MlModel.spec.MlModel.cost_fn = `MPL && false then
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
let median code a b t1 t2 = match a,b with
    | Integerized a, Integerized b -> 
        assert( 0 = MlModel.compare a.imodel b.imodel );
        let heur_cm = match t1,t2 with
            | Some t1,Some t2 -> 
                SeqCS.make_default_heuristic 
                    ~c3:(Cost_matrix.Three_D.default)
                        (MlModel.model_to_cm a.imodel (t1+.t2))
            | _,_ -> failwith "branches not specified by caller"
        in
        let aseq = { a.seq with SeqCS.heuristic = heur_cm; }
        and bseq = { b.seq with SeqCS.heuristic = heur_cm; } in
        let medn = { a with seq = SeqCS.median code aseq bseq; ia = None; } in
        remove_ambiguities (Integerized medn)
    | FPAlign a, FPAlign b -> 
        let bl = match t1,t2 with
            | Some t1, Some t2 -> t1 +. t2
            | _,_ -> failwith "branches not specified by caller"
        in
        let cost = ref 0.0 in
        let meds = 
            Array_ops.map_2
                (fun sa sb -> 
                    let mem = FloatSequence.FloatAlign.get_mem sa sb in
                    FloatSequence.FloatAlign.clear_mem mem;
                    let cst,med = FloatSequence.FloatAlign.median_2_cost sa sb a.fmodel bl mem in
                    assert( verify_cost cst sa sb a.fmodel bl mem );
                    if debug then begin
                        Printf.printf "\nA:%d:\t" (a.code); FloatSequence.FloatAlign.print_s sa;
                        Printf.printf "\nB:%d:\t" (a.code); FloatSequence.FloatAlign.print_s sb;
                        Printf.printf "\nC:%d:\t" (a.code); FloatSequence.FloatAlign.print_s med;
                        Printf.printf "\nCost: %f\tTime: %f\n%!" cst bl
                    end;
                    cost := !cost +. cst;
                    med)
                (a.ss)
                (b.ss)
        in
        FPAlign { a with cost = !cost; ss = meds }

    | _ , _ -> assert false

let readjust c1 c2 mine t1 t2 = 
    let internal_loop (t:float) : t * float =
        let d = median (code mine) c1 c2 (Some t) (Some t2) in
        (d,total_cost d)
    in
    let nt,(nmine,ncost) = 
        Numerical.brents_method (t1,(mine,total_cost mine)) internal_loop
    in
    (not (nt = t1)),total_cost mine,ncost,(nt,t2),nmine

let median_i code a b (t1:float) (t2:float) : t * float * float =
    let start = t1+.t2 in
    let _,_,ncost,(nt,_),nmine = 
        readjust a b (median code a b (Some t1) (Some t2)) (t1+.t2) 0.0
    in
    let half_time = nt /. 2.0 in
    if debug_est then
        Printf.printf "Optimized Branch of %d from %f --> %f: %f\n%!" 
                        code start nt ncost;
    nmine, half_time, nt -. half_time

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
let distance missing_distance a b = match a,b with
    | Integerized a,Integerized b -> SeqCS.distance missing_distance a.seq b.seq
    | FPAlign _ , FPAlign _       ->
        let dist1,dist2 = estimate_time a b in
        total_cost (median (-1) a b (Some dist1) (Some dist2))
    | _ , _ -> assert false

let dist_2 delta n a b = match a,b,n with
    | Integerized a,Integerized b,Integerized n -> SeqCS.dist_2 delta n.seq a.seq b.seq
    | FPAlign at, FPAlign bt,FPAlign nt -> 
        let dist1,dist2 = estimate_time a b in
        let tmp = median (-1) a b (Some dist1) (Some dist2) in
        let dist1,dist2 = estimate_time n tmp in
        total_cost (median (-1) a b (Some dist1) (Some dist2))
    | _ , _ , _ -> assert false
            
let tabu_distance a = match a with
    | Integerized a -> SeqCS.tabu_distance a.seq
    | FPAlign a -> failwith_todo "tabu_distance (FPAlign)"

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

let f_codes s c = match s with
    | Integerized s -> 
        Integerized { s with seq = SeqCS.f_codes s.seq c; ia = None; }
    | FPAlign s ->
        let codes,seqs = 
            array_filter (fun x -> All_sets.Integers.mem x c) s.codes s.ss
        in
        FPAlign { s with ss = seqs; codes = codes; }

let f_codes_comp s c = match s with
    | Integerized s ->
        Integerized { s with seq = SeqCS.f_codes_comp s.seq c; ia = None; }
    | FPAlign s ->
        let codes,seqs =
            array_filter (fun x -> not (All_sets.Integers.mem x c)) s.codes s.ss
        in
        FPAlign { s with ss = seqs; codes = codes; }


(*---- make an initial leaf node; still requires IA. *)
let make s m = match m.MlModel.spec.MlModel.cost_fn with
    | `ILK -> Integerized { seq = s; imodel= m; ia = None; }
    | `FLK ->
        let data = s_of_seq s in
        assert( (Array.length data) = (Array.length s.SeqCS.codes) );
        FPAlign {
                    ss = data;
                 codes = s.SeqCS.codes;
                fmodel = m;
                  code = s.SeqCS.code;
                  cost = 0.0;   }

    | _ -> failwith "not done yet"

let prior a = 
    let prior_of_seq alph priors (acc:float) sequence =
        Sequence.fold_right
            (fun acc i -> 
                let states = MlModel.list_of_packed i in
                let length = float_of_int (List.length states) in
                let c_pi = 
                    List.fold_left
                        (fun acc i -> acc +. (priors.{i} /. length))
                        (0.0)
                        (states)
                in
                c_pi *. acc)
            (acc)
            (sequence)
    in
    match a with
    | Integerized a ->
        let f = prior_of_seq a.imodel.MlModel.alph a.imodel.MlModel.pi_0 in
        Array.fold_left
            (fun acc data -> match data with
                | SeqCS.Heuristic_Selection x -> f acc x.SeqCS.DOS.sequence
                | SeqCS.Partitioned x ->
                    Array.fold_left
                        (fun acc -> function
                            | SeqCS.PartitionedDOS.DO x
                            | SeqCS.PartitionedDOS.First x
                            | SeqCS.PartitionedDOS.Last x ->
                                f acc x.SeqCS.DOS.sequence)
                        (acc) (x)
                | SeqCS.Relaxed_Lifted (x,_) ->
                    Array.fold_left f acc x.SeqCS.RL.sequence_table)
            (0.0) (a.seq.SeqCS.characters);
    | FPAlign a ->
        Array.fold_left
            (fun acc i ->
                prior_of_seq a.fmodel.MlModel.alph a.fmodel.MlModel.pi_0
                             acc (FloatSequence.FloatAlign.seq_of_s i))
            (0.0) 
            (a.ss)

ELSE

(* empty; required functions when likelihood is not enabled *)
type t = unit

let alph _ = failwith MlStaticCS.likelihood_error
let total_cost _ = failwith MlStaticCS.likelihood_error
let get_cm _ = failwith MlStaticCS.likelihood_error
let code _ = failwith MlStaticCS.likelihood_error
let get_codes _ = failwith MlStaticCS.likelihood_error
let combine _ _ = failwith MlStaticCS.likelihood_error
let make _ _ = failwith MlStaticCS.likelihood_error
let median _ _ _ _ _ = failwith MlStaticCS.likelihood_error
let median_i _ _ _ _ _ = failwith MlStaticCS.likelihood_error
let median_3 _ _ _ _ = failwith MlStaticCS.likelihood_error
let readjust _ _ _ _ _ = failwith MlStaticCS.likelihood_error
let distance _ _ _ = failwith MlStaticCS.likelihood_error
let to_string _ = failwith MlStaticCS.likelihood_error
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

ENDIF
