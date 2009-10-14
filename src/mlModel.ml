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

let () = SadmanOutput.register "MlMatrix" "$Revision"

let epsilon = 0.000001
let (=.) a b = abs_float (a-.b) < epsilon

type site_var = 
    (* #categories, alpha, beta, %invar *)
    | Gamma of int * float * float
    | Theta of int * float * float * float
    | Constant 

type subst_model =
    | JC69  of float
    | F81   of float
    | K2P   of float * float
    | F84   of float * float
    | HKY85 of float * float
    | TN93  of float * float * float
    | GTR   of float array
    | File  of float array array 

type priors = 
    (* diferented to determine if we should iterate *)
    | Estimated of float array
    | Given     of float array

type spec = {
    substitution : subst_model;
    site_variation : site_var option;
    base_priors : priors;
    use_gap : bool;
}

type model = {
    name: string;
    pi_0: (float, Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array1.t;
    alpha : float option; invar : float option; sites: int;
    rate: (float, Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array1.t;
    prob: (float, Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array1.t;
    u: (float, Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array2.t;
    d: (float, Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array2.t;
    ui:(float, Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array2.t option; 
}

(* EXTERNAL FUNCTIONS -- maintained in likelihood.c *)

(* diagonalize a symmetric or gtr matrix, WARNING: modifies passed matrices *)
external diagonalize_gtr: (* U D Ui *) FMatrix.m ->
    (float,Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array2.t ->
    (float,Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array2.t ->
    (float,Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array2.t ->
    unit = "likelihood_CAML_diagonalize_gtr"
external diagonalize_sym: (* U D *) FMatrix.m ->
    (float,Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array2.t ->
    (float,Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array2.t ->
    unit = "likelihood_CAML_diagonalize_sym"

(* compose matrices -- for testing purposes, as this composition is
 * usually done on the C side exclusively *)
external compose_gtr: (* U D Ui t *) FMatrix.m ->
    (float,Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array2.t -> 
    (float,Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array2.t ->
    (float,Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array2.t -> float
    -> (float,Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array2.t =
        "likelihood_CAML_compose_gtr"
external compose_sym: (* U D *) FMatrix.m ->
    (float,Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array2.t ->
    (float,Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array2.t -> float
    -> (float,Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array2.t =
        "likelihood_CAML_compose_sym"

(* calculates the gamma rates for specific alpha, beta and #classes *)
external gamma_rates: float -> float -> int ->
    (float,Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array1.t = 
        "gamma_CAML_rates"

(* LOCAL FUNCTIONS -- *)

(* divide a matrix by the mean rate so it will equal 1 *)
let m_meanrate srm pi_ =
    let mr = ref 0.0 and a_size = Bigarray.Array2.dim1 srm in
    for i = 0 to (a_size-1) do
        mr := !mr +. (~-.(srm.{i,i}) *. pi_.(i));
    done;
    for i = 0 to (a_size-1) do
        for j = 0 to (a_size-1) do
            srm.{i,j} <- srm.{i,j} /. !mr;
        done;
    done

(* val jc69 :: ANY ALPHABET size *)
let m_jc69 pi_ mu a_size =
    let srm = Bigarray.Array2.create Bigarray.float64 Bigarray.c_layout a_size a_size in
    Bigarray.Array2.fill srm mu;
    let diag = -. mu *. float (a_size-1) in
    for i = 0 to (a_size-1) do
        srm.{i,i} <- diag
    done;
    m_meanrate srm pi_;
    srm

(* val k2p :: only 4 or 5 characters *)
let m_k2p pi_ alpha beta a_size =
    let srm = Bigarray.Array2.create Bigarray.float64 Bigarray.c_layout a_size a_size in
    Bigarray.Array2.fill srm beta;
    (* modify transition elements to alpha *)
    srm.{1, 3} <- alpha; srm.{3, 1} <- alpha;
    srm.{2, 0} <- alpha; srm.{0, 2} <- alpha;
    (* set up the diagonal elements *)
    let diag = if a_size = 4 then
                    -. alpha -. beta -. beta
               else -. alpha -. beta *. 3.0  in
    for i = 0 to (a_size-1) do
        srm.{i,i} <- diag
    done;
    m_meanrate srm pi_;
    srm

(* val tn93 :: only 4 or 5 characters *)
let m_tn93 pi_ alpha beta gamma a_size =
    let srm = Bigarray.Array2.create Bigarray.float64 Bigarray.c_layout a_size a_size in
    Bigarray.Array2.fill srm gamma;
    srm.{0,2} <- alpha; srm.{1,3} <- beta; (* ACGT -- R=AG -- Y=CT *)
    srm.{2,0} <- alpha; srm.{3,1} <- beta; (* 0123 -- R=02 -- Y=13 *)
    for i = 0 to (a_size-1) do
        for j = 0 to (a_size-1) do
            srm.{i,j} <- srm.{i,j} *. pi_.(j);
        done;
    done;
    (* normalize diagonal so row sums to 0 *)
    for i = 0 to (a_size-1) do
        let diag = ref 0.0 in
        for j = 0 to (a_size-1) do
            if (i <> j) then diag := !diag +. srm.{i,j};
        done;
        srm.{i,i} <- -. !diag;
    done;
    m_meanrate srm pi_;
    srm

let m_tn93_ratio pi_ kappa1 kappa2 a_size = 
    let beta = (pi_.(0) *. pi_.(2) *. kappa1) +. (pi_.(1) *. pi_.(3) *. kappa2) +.
                ((pi_.(0) +. pi_.(2)) *. (pi_.(1)+.pi_.(3)) ) in
    let beta = 1.0 /. (2.0 *. beta) in
    let alpha1 = kappa1 *. beta and alpha2 = kappa2 *. beta in
    m_tn93 pi_ alpha1 alpha2 beta a_size

(* val f81 :: ANY ALPHABET size *)
let m_f81 pi_ lambda a_size =
    let srm = Bigarray.Array2.create Bigarray.float64 Bigarray.c_layout a_size a_size in
    for i = 0 to (a_size-1) do
        for j = 0 to (a_size-1) do
            srm.{i,j} <- 
                if i = j then
                     -. lambda *. (1.0 -. pi_.(i))
                else 
                    pi_.(j)*.lambda;
        done;
    done;
    m_meanrate srm pi_;
    srm

(* val hky85_ratio :: only 4 or 5 characters *)
let m_hky85_ratio pi_ kappa a_size =
    let beta = (pi_.(0) *. pi_.(2) *. kappa) +. (pi_.(1) *. pi_.(3) *. kappa) +.
                ((pi_.(0) +. pi_.(2)) *. (pi_.(1)+.pi_.(3)) ) in
    let beta = 1.0 /. (2.0 *. beta) in
    let alpha = kappa *. beta in
    m_tn93 pi_ alpha alpha beta a_size

(* val hky85 :: only 4 or 5 characters *)
let m_hky85 pi_ alpha beta a_size = m_tn93 pi_ alpha alpha beta a_size

(* val f84 :: only 4 or 5 characters *)
let m_f84 pi_ gamma kappa a_size =
    let y = pi_.(1) +. pi_.(3) in (* Y = C + T *)
    let r = pi_.(0) +. pi_.(2) in (* R = A + G *)
    let alpha = (1.0+.kappa/.r) *. gamma in
    let beta = (1.0+.kappa/.y) *. gamma in
    m_tn93 pi_ alpha beta gamma a_size 

(* val gtr :: ANY ALPHABET size
   pi_ == n; co_ == ((1+n)*n)/2
   form of lower packed storage mode, excluding diagonal, *)
let m_gtr pi_ co_ a_size =
    let n = ref 0 in (* array index *)
    let srm = Bigarray.Array2.create Bigarray.float64 Bigarray.c_layout a_size a_size in
    for i = 0 to (a_size-1) do
        for j = (i+1) to (a_size-1) do
            srm.{i,j} <- co_.(!n) *. pi_.(j);
            srm.{j,i} <- co_.(!n) *. pi_.(i);
            incr n;
        done;
    done;
    (* normalize diagonal so row sums to 0 *)
    for i = 0 to (a_size-1) do begin
        let diag = ref 0.0 in
        for j = 0 to (a_size-1) do
            if (i <> j) then diag := !diag +. srm.{i,j};
        done;
        srm.{i,i} <- -. !diag;
    end; done;
    (* divide through by mean-rate *)
    m_meanrate srm pi_;
    srm

(* val m_file :: any alphabet size -- recomputes diagonal and divides by meanrate *)
let m_file pi_ f_rr a_size =
    assert(a_size = Array.length f_rr);
    let srm = Bigarray.Array2.create Bigarray.float64 Bigarray.c_layout a_size a_size in
    for r = 0 to (a_size-1) do
        assert(a_size = Array.length f_rr.(r));
        let diag = ref 0.0 in
        for c = 0 to (a_size-1) do
            if (c <> r) then 
                (diag := !diag +. f_rr.(r).(c); srm.{r,c} <- f_rr.(r).(c);)
        done;
        srm.{r,r} <- ~-. !diag;
    done;
    m_meanrate srm pi_; (* TODO: should this be done if given? *)
    srm

(* test substitution matrix by diagonalizing, and composing with [t], uses GTR *)
let compose_model sub_mat t = 
    let a_size = Bigarray.Array2.dim1 sub_mat in
    let (u_,d_,ui_) = 
        let n_d = Bigarray.Array2.create Bigarray.float64 Bigarray.c_layout a_size a_size
        and n_ui = Bigarray.Array2.create Bigarray.float64 Bigarray.c_layout a_size a_size in
        Bigarray.Array2.fill n_d 0.0;
        Bigarray.Array2.fill n_ui 0.0;
        diagonalize_gtr FMatrix.scratch_space sub_mat n_d n_ui;
        (sub_mat, n_d, n_ui)
    in
    compose_gtr FMatrix.scratch_space u_ d_ ui_ t

(* create a model based on a specification and an alphabet *)
let create alph lk_spec = 
    let (a_size,a_gap) = match lk_spec.use_gap with
        | true -> Alphabet.size alph, (-1)
        | false -> (Alphabet.size alph) - 1, Alphabet.get_gap alph
    in
    (* set up all the probability and rates *)
    let variation,probabilities,alpha,invar,sites =
        match lk_spec.site_variation with
        | None ->
            Bigarray.Array1.of_array Bigarray.float64 Bigarray.c_layout ([| 1.0 |]),
            Bigarray.Array1.of_array Bigarray.float64 Bigarray.c_layout ([| 1.0 |]),
            None,None,1
        | Some a -> begin match a with
            | Constant ->  (* same as NONE *)
                 Bigarray.Array1.of_array Bigarray.float64 Bigarray.c_layout ([| 1.0 |]),
                 Bigarray.Array1.of_array Bigarray.float64 Bigarray.c_layout ([| 1.0 |]),
                 None,None,1
            | Gamma (x,y,z) -> (* SITES,ALPHA,BETA *)
                let p = Bigarray.Array1.create Bigarray.float64 Bigarray.c_layout x in
                Bigarray.Array1.fill p (1.0 /. (float_of_int x));
                gamma_rates y z x,p,Some y,None,x
            | Theta (w,x,y,z) -> (* GAMMA_SITES,ALPHA,BETA,PERCENT_INVAR *)
                if w < 1 then
                    failwith "Number of rate categories must be >= 1 (if invar w/out gamma)"
                else if w = 1 then begin
                    Bigarray.Array1.of_array Bigarray.float64 Bigarray.c_layout ([| 1.0 |]),
                    Bigarray.Array1.of_array Bigarray.float64 Bigarray.c_layout ([| 1.0 |]),
                    None,Some z,1
                end else begin
                    let p = Bigarray.Array1.create Bigarray.float64 Bigarray.c_layout w in
                    Bigarray.Array1.fill p (1.0 /. (float_of_int w));
                    let r = gamma_rates x y w in
                    r,p,Some x,Some z,w
                end
        end
    in
    (* check the rates so SUM(r_k*p_k) == 1 and SUM(p_k) == 1 |p| == |r| *)
    assert( (Bigarray.Array1.dim probabilities) = (Bigarray.Array1.dim variation) );
    let rsps = ref 0.0 and ps = ref 0.0 in
    for i = 0 to (Bigarray.Array1.dim probabilities) - 1 do
        rsps := !rsps +. ((Bigarray.Array1.get probabilities i) *. (Bigarray.Array1.get variation i));
        ps := !ps +. Bigarray.Array1.get probabilities i;
    done;
    assert( !rsps =. 1.0 && !ps =. 1.0 );
    (* extract the prior probability *)
    let priors = match lk_spec.base_priors with
        | Estimated p
        | Given p -> assert(a_size = Array.length p); p
    in
    (*  get the substitution rate matrix and set sym variable and to_formatter vars *)
    let sym,mname,sub_mat = match lk_spec.substitution with
        | JC69 lambda -> 
                true,  "JC69",  m_jc69 priors lambda a_size
        | K2P (alpha,beta) ->
                true,  "K2P",   m_k2p priors alpha beta a_size
        | F81 lambda ->
                false, "F81",   m_f81 priors lambda a_size
        | F84 (kappa,beta) ->
                false, "F84",   m_f84 priors kappa beta a_size
        | HKY85 (alpha,beta) ->
                false, "HKY85", m_hky85_ratio priors (alpha /. beta) a_size
        | TN93 (alpha1,alpha2,beta) -> 
                false, "TN93",  m_tn93 priors alpha1 alpha2 beta a_size
        | GTR coeff ->
                false, "GTR",   m_gtr priors coeff a_size
        | File matrix ->
                false, "File",  m_file priors matrix a_size
    in
    (* diagonalize to get factored probability matrix --submat is destroyed
     * but can be reconstructed through the model parameters *)
    let (u_,d_,ui_) = 
        let n_d = Bigarray.Array2.create Bigarray.float64 Bigarray.c_layout a_size a_size in
        Bigarray.Array2.fill n_d 0.0;
        if sym then begin
            diagonalize_sym FMatrix.scratch_space sub_mat n_d;
            (sub_mat, n_d, None)
        end else begin
            let n_ui = Bigarray.Array2.create
                            Bigarray.float64 Bigarray.c_layout a_size a_size in
            Bigarray.Array2.fill n_ui 0.0;
            diagonalize_gtr FMatrix.scratch_space sub_mat n_d n_ui;
            (sub_mat, n_d, Some n_ui)
        end
    in
    { rate = variation;
      prob = probabilities;
      name = mname;
     alpha = alpha;
     invar = invar;
     sites = sites;
      pi_0 = Bigarray.Array1.of_array Bigarray.float64 Bigarray.c_layout priors;
         u = u_; d = d_;ui = ui_ }


(* some testing functions and for pure ocaml median functions *)
let compose model t = match model.ui with
    | Some ui -> compose_gtr FMatrix.scratch_space model.u model.d ui t
    | None    -> compose_sym FMatrix.scratch_space model.u model.d t

