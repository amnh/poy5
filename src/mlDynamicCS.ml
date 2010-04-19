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
*)
type t =
    {     seq : SeqCS.t; 
        model : MlModel.model;
           ia : MlStaticCS.t option; }


(*---- non-external helper functions *)
let debug = true
let failwith_todo func = 
    failwith ("Your lazy developer needs to write: MlDynamicCS."^func)
let debug_printf format = if debug then Printf.printf format else ()


(*---- basic function set for gathering data of the module *)
let alph t        = t.seq.SeqCS.alph
let get_cm t      = t.seq.SeqCS.heuristic.SeqCS.c2
let code t        = t.seq.SeqCS.code
let to_string t   = failwith_todo "to_string"
let cardinal t    = SeqCS.cardinal t.seq
let encoding e t  = SeqCS.encoding e t.seq
let name_string t = match t.ia with 
                  | Some _ -> "Dynamic Likelihood (with ia)"
                  | None   -> "Dynamic Likelihood (w/out ia)"

let total_cost t = match t.ia with
    | Some x -> MlStaticCS.median_cost x
(*    | None when debug ->
        Status.user_message Status.Warning 
            "MlDynamicCS.total_cost; no implied alignment found";
        t.seq.SeqCS.total_cost *)
    | None   -> t.seq.SeqCS.total_cost

let get_codes t  = t.seq.SeqCS.codes

let combine t (b:MlStaticCS.t) = {t with ia = Some b}

let leaf_sequences {seq = a} = 
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
                    [|`DO t.SeqCS.RL.sequence_table.(p)|]) !map
    done;
    !map


(*---- to formatter, and printing functions *)
let to_formatter attr t par_opt bl d = match t.ia with
    | Some t -> MlStaticCS.to_formatter attr t bl d
    | None   -> failwith "MlDynamicCS.to_formatter does not have an implied alignment"


(*---- special functions for likelihood in a dynamic context *)
let min_bl = MlStaticCS.minimum_bl ()
let estimate_time a b = 
    (* set up the IA functor *)
    let is_transition one two = 
        let atran =
            try match Alphabet.complement one b.model.MlModel.alph with
                | None -> false
                | Some x -> true
            with | _ -> false
        and btran = 
            try match Alphabet.complement two a.model.MlModel.alph with
                | None -> false
                | Some x -> true
            with | _ -> false
        in
        atran || btran
    in
    let aligned_triples = SeqCS.align_2 a.seq b.seq in
    assert( (a.model.MlModel.alph) = (b.model.MlModel.alph) );
    let gchar = Alphabet.get_gap a.model.MlModel.alph in
    let same = ref 0 (* the number of aligned columns *)
    and tran = ref 0 (* the number of transitions aligned *)
    and gaps = ref 0 (* the number of characters aligned with gaps *)
    and othr = ref 0 (* other aligned characters --transversions *) in

    Array.iter (fun triple_array -> 
    Array.iter (fun (align_a,align_b,_) -> 
            for n = 0 to (Sequence.length align_a) - 1 do
                let achar = Sequence.get align_a n 
                and bchar = Sequence.get align_b n in
                if achar = bchar then incr same
                else if achar = gchar then incr gaps
                else if bchar = gchar then incr gaps
                else if is_transition achar bchar then incr tran
                else incr othr
            done;)
        triple_array)
        aligned_triples;

    let total = (!same + !tran + !gaps + !othr) in
    (* below, similar to mlStatic estimate_time function *)
    let proportion = (float_of_int !same) /. (float_of_int total) in
    let p = match (1.0 -. proportion) with
        | x when x < 0.75 -> x
        | x -> 0.70
    in
    let nt2 = ~-. 0.75 *. (log (1.0 -. (p *. 4.0 /. 3.0))) in
    let nt = if nt2 <= min_bl then min_bl else nt2 /. 2.0 in
(*    Printf.printf*)
(*        ("Aligned      : %d\nTransitions  : %d\nAligned Gaps : %d\n"^^*)
(*         "Transversions: %d\nTotals       : %d\nProportion   : %f\n\n"^^*)
(*         "Distance     : %f\n\n%!")*)
(*        !same !tran !gaps !othr total proportion nt2;*)
    (nt,nt)


(*---- median functions *)
let median code a b t1 t2 =
    assert( 0 = MlModel.compare a.model b.model);
    let heur_cm = match t1,t2 with
        | Some t1,Some t2 -> 
            SeqCS.make_default_heuristic (MlModel.model_to_cm a.model (t1+.t2))
        | _,_ -> failwith "branches not specified by caller"
    in
    let aseq = {a.seq with SeqCS.heuristic = heur_cm; }
    and bseq = {b.seq with SeqCS.heuristic = heur_cm; } in
    {
          seq = SeqCS.median code aseq bseq;
        model = a.model;
           ia = None; 
    }

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
let distance missing_distance a b = SeqCS.distance missing_distance a.seq b.seq
let dist_2 delta n a b            = SeqCS.dist_2 delta n.seq a.seq b.seq
let tabu_distance a               = SeqCS.tabu_distance a.seq


(*---- filter functions *)
let f_codes s c = 
    { seq = SeqCS.f_codes s.seq c; model = s.model; ia = None; }
let f_codes_comp s c = 
    { seq = SeqCS.f_codes_comp s.seq c; model = s.model; ia = None; }


(*---- make an initial leaf node; still requires IA. *)
let make s m = { seq = s; model = m; ia = None };

ELSE
type t = unit

let alph _          = failwith MlStaticCS.likelihood_error
let total_cost _    = failwith MlStaticCS.likelihood_error
let get_cm _        = failwith MlStaticCS.likelihood_error
let code _          = failwith MlStaticCS.likelihood_error
let combine _ _     = failwith MlStaticCS.likelihood_error
let make _ _        = failwith MlStaticCS.likelihood_error
let median _ _ _ _ _= failwith MlStaticCS.likelihood_error
let median_3 _ _ _ _= failwith MlStaticCS.likelihood_error
let distance _ _ _  = failwith MlStaticCS.likelihood_error
let to_string _     = failwith MlStaticCS.likelihood_error
let name_string _   = failwith MlStaticCS.likelihood_error
let dist_2 _ _ _ _  = failwith MlStaticCS.likelihood_error
let f_codes _ _     = failwith MlStaticCS.likelihood_error
let f_codes_comp _ _= failwith MlStaticCS.likelihood_error
let tabu_distance _ = failwith MlStaticCS.likelihood_error
let cardinal _      = failwith MlStaticCS.likelihood_error
let encoding _ _    = failwith MlStaticCS.likelihood_error
let leaf_sequences _= failwith MlStaticCS.likelihood_error
let to_formatter _ _ _ _ _ = failwith MlStaticCS.likelihood_error
ENDIF
