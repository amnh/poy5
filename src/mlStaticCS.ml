type mlchar = {
    ccode : int;
    p_v : float array;  (* probability vector *)
}
type cm = {
    p: float -> float array array;  (* probability matrix;arg- branch length *)
    pi_0: float array;              (* [prior] probability *)
}
type t = {
    code: int;             (* identifying code *)
    mle: float;            (* MLE of this char set *)
    model: cm;             (* model: probability vector * prior *)
    chars: mlchar array;   (* char set *)
    (* MLE is eventually stored in cost in r *)
}

(* empty types *)
let emp_c = { ccode = 0; p_v = [||] }
let emp_m = { p = (fun x -> [|[||]|]); pi_0 = [||]; }
let emp_t = { code = 0; mle = 0.0; chars = [| emp_c |]; model=emp_m; }

(* ZIP o MAP in 'one' step *)
let array_zipmap f a1 a2 = 
    assert (Array.length a1 == Array.length a2);
    let x = Array.make (Array.length a1) (f a1.(0) a2.(0)) in
        for i = 1 to ((Array.length a1)-1) do
            x.(i) <- f a1.(i) a2.(i)
        done; x

(* recall: sum of product of two vectors *)
let dot_product v1 p1 = 
    let r = array_zipmap ( *. ) v1 p1 in
    Array.fold_right (+.) r 0.0

(* log liklihood: float->(float array)->mlchar *)
let mle a b = log (dot_product a b.p_v)

(* (float->float)->mlchar->mlchar *)
let median_char p_1 p_2 a b = 
    let pa = a.p_v and pb = b.p_v in
    (** median_element a b p1 p2 x                      *)
    (* calculates one element of the probability vector *)
    (* denoted by 'x'; a/b are vectors of the nodes, &  *)
    (* p1 p2 are the probability matrices for each node *)
    let median_element a b p1 p2 x = 
        let x1 = Array.get p1 x in
        let x2 = Array.get p2 x in
            let right_sum = dot_product a x1 in
            let lefts_sum = dot_product b x2 in
                lefts_sum *. right_sum
    in
        let curried_c = median_element pa pb p_1 p_2 in
        let basis = 
            if Array.length a.p_v == 4 then 
                [| 0;1;2;3 |] 
            else
                (* TODO: linear algebra instead *)
                Array.init (Array.length a.p_v) (fun x -> x) 
        in
        { ccode = 0; p_v = Array.map curried_c basis; }


(* empty argument is 'previous' *)
let median _ a b t1 t2 = 
    let p_1 = a.model.p t1 and
        p_2 = b.model.p t2 in
    (* for computing new char prob vectors *)
    let medcc = median_char p_1 p_2 in
    let c_map = (array_zipmap medcc a.chars b.chars) in
    (* for computing new MLE *)
    let n_mle = Array.map (mle a.model.pi_0) c_map in
    {   a with
        chars = c_map;
        mle = Array.fold_right (+.) n_mle 0.0 ;
        model = a.model;
    }

(* find branch length between two nodes *)
let find_t a b = 0.1

let root_cost t = t.mle
let distance a b t1 t2 = 
    let t = median a a b t1 t2 in
    t.mle

(* (a (b)) -> (a (x (b n))) *)
let dist_2 n a b nt at bt xt= 
    let x = median a a b at bt in
    let tt = match xt with
            | None -> find_t x n
            | Some y -> y in
    let y = median x x n tt nt in
    y.mle

let reroot_median a b at bt = median a a b at bt 
let median_cost ta = ta.mle
let median_3 _ x _ _ =  x

let to_string a = "MLStaticCS"
let cardinal ta = Array.length ta.chars

(* filtering functions *)
let filter f t =
    let filter_e f e = Array_ops.filter f e in
    {t with chars = filter_e f t.chars }
let f_codes_comp t codes = 
    let check mlc = not (All_sets.Integers.mem mlc.ccode codes) in 
    filter check t
let f_codes t codes = 
    let check mlc = All_sets.Integers.mem mlc.ccode codes in
    filter check t

let compare_data a b = compare a b

let f_abs x = 
    if x < 0.0 then -.x else x

(* roam by step to find best result *)
let rec g_roam c1 c2 t1 t2 mle_ step epsilon =
    let best = 
        let left = (distance c1 c2 (t1-.step) t2,(t1-.step),t2) in
        let right= (distance c1 c2 (t1+.step) t2,(t1+.step),t2) in
             match left, right with
            |((l,_,_),(r,_,_)) -> if l > r then left else right in
    let best = 
        let up = (distance c1 c2 t1 (t2+.step),t1,(t2+.step)) in
             match best,up with
            |((b,_,_),(u,_,_)) -> if b > u then best else up in
    let best =
        let down = (distance c1 c2 t1 (t2-.step), t1, (t2-.step)) in
             match best,down with
            |((b,_,_),(d,_,_)) -> if b > d then best else down in
     match best with
    |(b_t,b_t1,b_t2) ->
            if b_t < mle_ || epsilon >= f_abs(b_t-.mle_) then (mle_,t1,t2) else
                g_roam c1 c2 b_t1 b_t2 b_t step epsilon

(* All_sets.Integers.t option ->  to_adjust
 * All_sets.Integers.t ->         !modified
 * float ->                       mine.time
 * t -> t -> t -> t ->            c1,c2,parent,mine
 * float -> float ->              current time/branch length for t1 and t2
 * (All_sets.Integers.t * float * float * float * t) *)
let readjust _ x _ c1 c2 _ mine t1 t2 = (x,mine.mle,mine.mle,(t1,t2),mine) 

(*    let data = g_roam c1 c2 t1 t2 mine.mle 0.1 0.0001 in
    match data with
    | (bmle,bt1,bt2) when bt1==t1 && bt2==t2 -> (* no change *)
            (x,mine.mle,mle,(bt1,bt2),mine)
    | (bmle,bt1,bt2) -> 
            (x,mine.mle,mle,(bt1,bt2), {mine with mle = bmle}
*)
    
    
(* Parser.SC.static_spec -> ((int list option * int) array) -> t *)
let of_parser spec characters = 

    (* for looping over characters array *)
    let loop_ (states,code) =
        match states with (* mlchar: ccode/p_v *)
        | None -> { ccode = code; p_v = Array.make 4 1.0; }
        | Some s ->
            let pv = Array.map (fun x -> float x) (Array.of_list s) in
            { ccode = code; p_v = pv; }
    in
    {
        code  = 0;
        mle   = 0.0;
        model = emp_m;
        chars = Array.map loop_ characters; 
    }

(* ------------------------------------------------ *)
(* models *)

(* val jc69 *)
let jc69 d = 
    let p_0 = 0.25 +. 0.75 *.  exp(-4.*. d /. 3.) in
    let p_1 = 1.0 -. p_0 in
     [| [|p_0;p_1;p_1;p_1|];
        [|p_1;p_0;p_1;p_1|];
        [|p_1;p_1;p_0;p_1|];
        [|p_1;p_1;p_1;p_0|];
     |]

(* val k80 : alpha/2*beta = 1             *) 
(* k = alpha/beta and d= (alpha+2*beta)/t *)
let k80 alpha beta t =
    (* assert (alpha /. (2.0 *. beta) == 1.0); *)
    let term_1 = 0.25 *. exp(-4.*.beta*.t) in
    let term_2 = 0.5 *. exp(-2. *. t *. (alpha +. beta) ) in
    
    let p_0 = 0.25 +. term_1 +. term_2 in
    let p_1 = 0.25 +. term_1 -. term_2 in
    let p_2 = 0.25 -. term_1 in

    [| [|p_0;p_1;p_2;p_2|];
       [|p_1;p_0;p_2;p_2|];
       [|p_2;p_2;p_0;p_1|];
       [|p_2;p_2;p_1;p_0|]
    |]
(* val k80: since r=a/2b=1, a = 2b, b=a/2 *) 
(* r = transition/transversion ratio      *)
let k80_r r t = 
    let a = r /. (r +. 1.0) in
    let b = 1.0 /. (2.0 *. (r +. 1.0)) in
    k80 a b t

(* val F81      -pi_ACTG     --TN93 with 1=k1=k2                *)
(* val F84      -pi_ACTGYR,k --TN93 with k1=(1+k/Y), k2=(1+k/R) *)
(* val HKY85    -pi_ACTG,k   --TN93 with k2=k1                  *)
(* val TN93     -pi_ACTG,k1,k2                                  *)
(* val REV      -pi_ACTG,a,b,c,d,e                              *)


(* ------------------------------- *)
(* manual test of single character *)
let ch_a = {ccode=0;p_v=[|0.0;0.0;1.0;0.0|]; }
let ch_c = {ccode=0;p_v=[|0.0;1.0;0.0;0.0|]; }
let ch_t = {ccode=0;p_v=[|1.0;0.0;0.0;0.0|]; }
let ch_g = {ccode=0;p_v=[|0.0;0.0;0.0;1.0|]; }
let s_model = { p = k80 0.5 0.25; pi_0 = [|0.25;0.25;0.25;0.25|]; }

let l_1 = {code=0;mle=0.0;model=s_model;chars=[|ch_c|]; }
let l_2 = {code=1;mle=0.0;model=s_model;chars=[|ch_t|]; }
let l_3 = {code=2;mle=0.0;model=s_model;chars=[|ch_c|]; }
let l_4 = {code=3;mle=0.0;model=s_model;chars=[|ch_c|]; }
let l_5 = {code=4;mle=0.0;model=s_model;chars=[|ch_a|]; }

let l_7 = median l_1 l_1 l_2 0.2 0.2
let l_6 = median l_7 l_7 l_5 0.1 0.2
let l_8 = median l_4 l_4 l_3 0.2 0.2
let l_0 = median l_6 l_6 l_8 0.1 0.1
(* correct when...
val l_0 : t =
  {code = 0; mle = -7.58140757255770126;
   model = {p = <fun>; pi_0 = [|0.25; 0.25; 0.25; 0.25|]};
   chars =
    [|{ccode = 0;
       p_v =
        [|0.000112371114942437172; 0.00183822623771382442;
          7.51370241327383381e-05; 1.36379295632049643e-05|]}|]}
*)

(* ------------------------------------------------------------ *)
(* manual test with two same characters (double score as above) *)
let l_1 = {code=0;mle=0.0;model=s_model;chars=[|ch_c;ch_c|]; }
let l_2 = {code=1;mle=0.0;model=s_model;chars=[|ch_t;ch_t|]; }
let l_3 = {code=2;mle=0.0;model=s_model;chars=[|ch_c;ch_c|]; }
let l_4 = {code=3;mle=0.0;model=s_model;chars=[|ch_c;ch_c|]; }
let l_5 = {code=4;mle=0.0;model=s_model;chars=[|ch_a;ch_a|]; }

let l_7 = median l_1 l_1 l_2 0.2 0.2
let l_6 = median l_7 l_7 l_5 0.1 0.2
let l_8 = median l_4 l_4 l_3 0.2 0.2
let l_0 = median l_6 l_6 l_8 0.1 0.1
(* correct when...
  {code = 0; mle = -15.1628151451154025;
   model = {p = <fun>; pi_0 = [|0.25; 0.25; 0.25; 0.25|]};
   chars =
    [|{ccode = 0;
       p_v =
        [|0.000112371114942437172; 0.00183822623771382442;
          7.51370241327383381e-05; 1.36379295632049643e-05|]};
      {ccode = 0;
       p_v =
        [|0.000112371114942437172; 0.00183822623771382442;
          7.51370241327383381e-05; 1.36379295632049643e-05|]}|]}
*)
(* ------------------------------------------------------------ *)
