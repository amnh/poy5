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
let () = SadmanOutput.register "MlStaticCS" "$Revision %r $"

let debug = false 

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

(* val jc69 *)
let jc69 mu t = 
    let p_1 = 0.25 -. 0.25 *. exp(-.t*.mu) in
    let p_0 = 1.0 -. (3.0 *. p_1 ) in

        (* A   C   G   T *)
     [| [|p_0;p_1;p_1;p_1|];
        [|p_1;p_0;p_1;p_1|];
        [|p_1;p_1;p_0;p_1|];
        [|p_1;p_1;p_1;p_0|];
     |]

(* val k80: since r=a/2b=1, a = 2b, b=a/2 *) 
(* r = transition/transversion ratio      *)
let k2p k t = 
    let term1 = 0.25 *. exp ( -4.*.t/.(k+.2.) ) in
    let term2 = 0.5 *. exp (-2.*.t*.(k+.1.)/.(k+.2.)) in
    
    let p_0 = 0.25 +. term1 +. term2 and
        p_1 = 0.25 +. term1 -. term2 and
        p_2 = 0.25 -. term1 in

       (* A   C   G   T *)
    [| [|p_0;p_2;p_1;p_2|];
       [|p_2;p_0;p_2;p_1|];
       [|p_1;p_2;p_0;p_2|];
       [|p_2;p_1;p_2;p_0|];
    |]

(* printing commands *)
let print_float x = Printf.printf "%2.10f\t" x
let print_array xs = Array.iter print_float xs; print_newline ()

(* ZIP o MAP *)
let array_zipmap f a1 a2 =
    let x = Array.make (Array.length a1) (f a1.(0) a2.(0)) in
        for i = 1 to ((Array.length a1)-1) do
            x.(i) <- f a1.(i) a2.(i)
        done; x

(* recall: sum of product of two vectors *)
let dot_product v1 v2 = 
    let r = array_zipmap ( *. ) v1 v2 in
    let res = Array.fold_right (+.) r 0.0 in
    res

(* negative log liklihood *)
let mle a model = 
    let res = -.log (dot_product a model.p_v) in
    res

(* calculate a new nodes mle and prob_vectors *)
let median_char p_1 p_2 a b = 
    let pa = a.p_v and pb = b.p_v in
    (* calculates one element of the probability vector *)
    let median_element a b p1 p2 x = 
        let x1 = Array.get p1 x in
        let x2 = Array.get p2 x in
            let right_sum = dot_product a x1 in
            let lefts_sum = dot_product b x2 in
                lefts_sum *. right_sum
    in
        let curried_c = median_element pa pb p_1 p_2 in
        let npv = Array.init (Array.length a.p_v) curried_c in
        { a with p_v = npv; }

(* empty argument is 'previous' *)
let median _ a_node b_node t1 t2 =
    let p_1 = a_node.model.p t1 and
        p_2 = b_node.model.p t2 in
    (* for computing new char prob vectors *)
    let medcc = median_char p_1 p_2 in
    let c_map = (array_zipmap medcc a_node.chars b_node.chars) in
    (* for computing new MLE *)
    let n_mle = Array.map (mle a_node.model.pi_0) c_map in
    let n_mle = Array.fold_right (+.) n_mle 0.0 in
    { a_node with
        chars = c_map;
        mle = n_mle;
    }

let find_t a b = 0.0

let root_cost t = t.mle

let distance a_node b_node t1 t2 = 
    let t = median a_node a_node b_node t1 t2 in
    t.mle

(* insert a node between two *)
(* (c, b) -> (c,(a),b)  *)
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

let f_abs x = if x < 0.0 then -.x else x
(* roam by step to find best result *)
let rec g_roam c1 c2 t1 t2 mle_ step epsilon =
    let best = 
        let left = if (t1-.step) < 0.0 then (mle_,t1,t2) else
            (distance c1 c2 (t1-.step) t2,(t1-.step),t2) in
        let right= if (t1+.step) < 0.0 then (mle_,t1,t2) else
            (distance c1 c2 (t1+.step) t2,(t1+.step),t2) in
        match left, right with
            |((l,_,_),(r,_,_)) -> if l < r then left else right in
    let best = 
        let up = if (t2+.step) < 0.0 then (mle_,t1,t2) else
            (distance c1 c2 t1 (t2+.step),t1,(t2+.step)) in
             match best,up with
            |((b,_,_),(u,_,_)) -> if b < u then best else up in
    let best =
        let down = if (t2-.step) < 0.0 then (mle_,t1,t2) else
            (distance c1 c2 t1 (t2-.step), t1, (t2-.step)) in
             match best,down with
            |((b,_,_),(d,_,_)) -> if b < d then best else down in
     match best with
    |(b_t,b_t1,b_t2) ->
        if b_t >= mle_ || epsilon >= f_abs(b_t-.mle_) 
        then (mle_,t1,t2)
        else g_roam c1 c2 b_t1 b_t2 b_t step epsilon

(* readjust the branch lengths to create better mle score *)
let readjust xopt x child1 child2 mine c_t1 c_t2 mine_t =
    let mmle = mine.mle in

    let data = g_roam child1 child2 c_t1 c_t2 mmle 0.01 0.001 in (* TODO: constants *)
    match data with
    | (bmle,bt1,bt2) when bt1==c_t1 && bt2==c_t2 ->
        (x,mine.mle,bmle,(bt1,bt2),mine)
    | (bmle,bt1,bt2) ->
        (*  Printf.printf  "%f  -->  %f  -->  %f\n" mine.mle mmle bmle; *)
        (* add all changed characters to set *)
        let x = Array.fold_right (* TODO: bottleneck *)
            (fun c s -> All_sets.Integers.add c.ccode s)
            mine.chars x in
        (x,mine.mle,bmle,(bt1,bt2), {mine with mle = bmle} )

(* of_parser stuff *)
let rec list_of n x =
    match n with
    | 0 -> []
    | e -> x :: (list_of (n-1) x)
let rec set_in num lst =
    match lst with
    | hd :: tl ->
       (match num with
        | 0 -> 1.0::tl
        | n -> hd:: set_in (n-1) tl)
    | [] ->
        match num with
        | 0 -> [1.0]
        | e -> (list_of (e) 0.0) @ [1.0]
let rec sublist l a b =
    match l with
    | hd::tl ->
       (match a,b with
        | 0,0 -> []
        | 0,_ -> hd::sublist tl 0 (b-1)
        | _,_ -> sublist tl (a-1) b)
    | _ -> []

(* Parser.SC.static_spec -> ((int list option * int) array) -> t *)
let of_parser spec characters = 
    let a_size = Alphabet.size(Alphabet.to_sequential(spec.Parser.SC.st_alph))-1 in
    let a_gap = Alphabet.get_gap( spec.Parser.SC.st_alph ) in

    let model = spec.Parser.SC.st_type in
    let model =
        match model with
        | Parser.SC.STLikelihood m -> m 
        | _ -> failwith "wrong model" in 

    let priors =
        match model.Parser.SC.base_priors with
        | Parser.SC.Estimated p
        | Parser.SC.Given p -> p in

    let p_funk =
        match model.Parser.SC.substitution with
        | Parser.SC.Constant lambda -> jc69 lambda
        | Parser.SC.K2P kappa -> k2p kappa
(*      | Parser.SC.F81 pi_s -> tn93 pi_s 1.0 1.0       *)
(*      | Parser.SC.F84 pi_se k -> 
                let k1 = 1.0+.k/.Y and k2 = 1.0+.k/.R in
                tn93 pi_s k1 k2                         *)
(*      | Parser.SC.HKY85 pi_s k                        *)
(*      | Parser.SC.TN93 pi_s k1 k2 -> tn93 pis k1 k2   *)
(*      | Parser.SC.REV ...                             *)
        | _ -> failwith ("Likelihood model not available")
    in

    (* loop for each character *)
    let loop_ (states,code) =
        match states with 
        | None -> { ccode = code;
                    p_v = Array.make a_size 1.0; } 
        | Some s when List.hd s = a_gap ->
                  { ccode = code;
                    p_v = Array.make a_size 1.0; }
        | Some s ->
                let pl = List.fold_right set_in s (list_of a_size 0.0) in
                let pv = Array.of_list (sublist pl 0 a_size) in
                { ccode = code; p_v = pv; }
    in
    {   code = model.Parser.SC.set_code;
         mle = 0.0;
       model = { pi_0 = Array.sub priors 0 a_size; p = p_funk };
       chars = Array.map loop_ characters; }

(* filtering functions *)
let filter f t =
    let filter_e f e = Array_ops.filter f e in
    { t with chars = filter_e f t.chars }
let f_codes_comp t codes = 
    let check mlc = not (All_sets.Integers.mem mlc.ccode codes) in 
    filter check t
let f_codes t codes = 
    let check mlc = All_sets.Integers.mem mlc.ccode codes in
    filter check t
let compare_data a b = compare a b

(* Tags.attributes ->t ->t option ->Data.d *)
let to_formatter attr mine minet _ data :Tags.output list =
    let i_pow b n = 
            let fb = float_of_int b and fn = float_of_int n in
            int_of_float (fb**fn) in

    let pack lst = `Structured (`Set (lst)) in
    (* list->map->map(i) for arrays to eventually pack *)
    let _ray f v = pack (Array.to_list 
                        (Array.map (fun x -> `Single x) (Array.map f v))) in
    let _rayi f v = pack (Array.to_list 
                        (Array.map (fun x -> `Single x) (Array.mapi f v))) in

    (** pack individual elements into an ID
     *   cs is the nucleotide id
     *   p is the value 
     **)
    let f_element cs p :Tags.output =
        let alpha = Alphabet.nucleotides in 
        (Alphabet.find_code (i_pow 2 cs) alpha, [], `String (string_of_float p)) in

    (** pack an array into an ID 
     *   cs is the nucleotide id
     *   p is an array
     **)
    let f_array cs (p:float array) :Tags.output =
        let alpha = Alphabet.nucleotides in
        (Alphabet.find_code (i_pow 2 cs) alpha, [], (_rayi f_element p) ) in

    let f_pmat m :Tags.output = (Tags.Characters.p_mat,[],(_rayi f_array m)) in
    let f_char c :Tags.output = 
        let a_char = (Tags.Data.code, string_of_int c.ccode) in
        (Tags.Characters.character, [a_char], (_rayi f_element c.p_v)) in
    let f_prior cm :Tags.output =
            (Tags.Characters.prior, [], (_rayi f_element cm.pi_0)) in 
    let f_chars cs :Tags.output = 
            (Tags.Characters.characters, [], (_ray f_char cs)) in
    let attrib :Tags.attribute list = 
                (Tags.Data.code,string_of_int mine.code) :: 
                (Tags.Nodes.time, string_of_float minet) ::
                (Tags.Characters.mle, string_of_float mine.mle) :: attr in
    let con = `Structured 
                (`Set [
                    `Single (Tags.Characters.model,[], `Structured
                        (`Set [
                            `Single (f_prior mine.model);
                            `Single (f_pmat (mine.model.p minet))]));
                    `Single (Tags.Data.characters,[],
                        `Structured 
                            (`Single (f_chars mine.chars)))] ) in

    [(Tags.Characters.likelihood, attrib, con)]
(* -> Tags.output list *)

(* -------------------------------------------------------------------
(* manual test of single character *)
let ch_a = {ccode=0;p_v=[|1.0;0.0;0.0;0.0|]; }
let ch_c = {ccode=0;p_v=[|0.0;1.0;0.0;0.0|]; }
let ch_t = {ccode=0;p_v=[|0.0;0.0;0.0;1.0|]; }
let ch_g = {ccode=0;p_v=[|0.0;0.0;1.0;0.0|]; }
let s_model = { p = k80_r 2.0; pi_0 = [|0.25;0.25;0.25;0.25|]; }

let l_1 = {code=0;mle=0.0;model=s_model;chars=[|ch_c|]; }
let l_2 = {code=1;mle=0.0;model=s_model;chars=[|ch_t|]; }
let l_3 = {code=2;mle=0.0;model=s_model;chars=[|ch_a|]; }

let l_7 = median l_1 l_1 l_2 0.2 0.2
let l_6 = median l_1 l_7 l_3 0.1 0.2
let l_8 = median l_1 l_1 l_1 0.2 0.2
let l_0 = median l_1 l_6 l_8 0.1 0.1
let nl_8= g_roam l_1 l_3 0.2 0.2 1.7544841622 0.1 0.1

(* correct when...
val l_0 : t =
  {code = 0; mle = 7.58140757255770126;
   model = {p = <fun>; pi_0 = [|0.25; 0.25; 0.25; 0.25|]};
   chars =
    [|{ccode = 0;
       p_v =
        [|0.000112371114942437172; 0.00183822623771382442;
          7.51370241327383381e-05; 1.36379295632049643e-05|]}|]}
*)

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

(*
PHYML and DNAML input:
  RATIO::
    PHYML:2.0, DNAML:4.0
  TREE::
    (((Gamma:0.0001,Epsilon:0.0001):0.0001,Delta:0.0001):0.0002,Beta:0.0001,Alpha:0.0001);
  SEQ::
    ----------
    5   2 
    Alpha     TT
    Beta      TT
    Gamma     AC
    Delta     AC
    Epsilon   AA
    ----------
  Note:: running poy with this data set produces this as it's first tree as
    well when all branch lengths are set to 0.0001. The 0.0002 in tree above
    is because the input is unrooted, so branch lengths are added.
*)

let a = {ccode=0;p_v=[|1.0;0.0;0.0;0.0|]; }
let c = {ccode=0;p_v=[|0.0;1.0;0.0;0.0|]; }
let t = {ccode=0;p_v=[|0.0;0.0;0.0;1.0|]; }
let g = {ccode=0;p_v=[|0.0;0.0;1.0;0.0|]; }
let s_model = { p = k80_r 2.0; pi_0 = [|0.25;0.25;0.25;0.25|]; }

let alpha = {code=0;mle=0.0;model=s_model;chars = [| t;t |]; }
let beta = {code=1;mle=0.0;model=s_model;chars =  [| t;t |]; }
let gamma = {code=2;mle=0.0;model=s_model;chars = [| a;c |]; }
let delta = {code=3;mle=0.0;model=s_model;chars = [| a;c |]; }
let epsilon = {code=4;mle=0.0;model=s_model;chars=[| a;a |]; }

let l_3 = median alpha epsilon gamma 0.0001 0.0001
let l_2 = median alpha delta l_3 0.0001 0.0001
let l_1 = median alpha alpha beta 0.0001 0.0001
let l_0 = median alpha l_2 l_1 0.0001 0.0001

*)
