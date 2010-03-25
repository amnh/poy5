
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
let alph t       = t.seq.SeqCS.alph
let get_cm t     = t.seq.SeqCS.heuristic.SeqCS.c2
let code t       = t.seq.SeqCS.code
let to_string t  = failwith_todo "to_string"
let cardinal t   = SeqCS.cardinal t.seq
let encoding e t = SeqCS.encoding e t.seq

let total_cost t = match t.ia with
    | Some x -> MlStaticCS.median_cost x
    | None   ->
(*        Printf.printf "Stealing Cost from Dynamic\n%!";*)
        t.seq.SeqCS.total_cost (* we need this for before IA is done *)

let get_codes t  = t.seq.SeqCS.codes

let combine t (b:MlStaticCS.t) = {t with ia = Some b}

(*---- to formatter, and printing functions *)
let to_formatter attr t par_opt d = failwith_todo "to_formatter"


(*---- special functions for likelihood in a dynamic context *)
let estimate_time a b = 0.05,0.05


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

and median_3 p n c1 c2 = 
    assert( 0 = MlModel.compare n.model c1.model);
    assert( 0 = MlModel.compare n.model c2.model);
    assert( 0 = MlModel.compare n.model p.model);
    { 
        seq = SeqCS.median_3 p.seq n.seq c1.seq c2.seq;
      model = n.model;
         ia = None; 
    }


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
