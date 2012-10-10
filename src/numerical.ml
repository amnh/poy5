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

let () = SadmanOutput.register "Numerical" "$Revision: 2727 $"

let (-->) b a = a b

let debug = false

let coarse_debug = false

let failwithf format = Printf.ksprintf failwith format

let warning_message format = Printf.ksprintf (Status.user_message Status.Warning) format

let debug_printf format = 
    Printf.ksprintf (fun x -> if debug then print_string x; flush stdout) format

and cdebug_printf format = 
    Printf.ksprintf (fun x -> if coarse_debug then print_string x; flush stdout) format

and pp_farray xs =
    (Array.fold_left (fun acc x -> acc^"| "^(string_of_float x)^" ") "[" xs)^" |]"

and pp_iarray xs =
    (Array.fold_left (fun acc x -> acc^"| "^(string_of_int x)^" ") "[" xs)^" |]"


(** {6 Constants} *)

let tolerance = 1e-6
let epsilon   = 1e-10
let minimum   = tolerance *. 2.0


(** {6 Floating Point Functions} *)

(** Check if a value is zero (takes into account negative; and subnormal *)
let is_zero x = match classify_float x with
        | FP_zero | FP_subnormal -> true
        | FP_infinite | FP_nan | FP_normal -> false

(** Check if a value is Not A Number *)
and is_nan x = match classify_float x with
        | FP_zero | FP_subnormal
        | FP_infinite | FP_normal -> false
        | FP_nan -> true

(** Check if a value is +/- infinity *)
and is_inf x = match classify_float x with
        | FP_infinite -> true
        | FP_nan | FP_zero | FP_subnormal | FP_normal -> false


(** {6 Special functions} *)

external gamma : float -> float = "gamma_CAML_gamma"

external lngamma : float -> float = "gamma_CAML_lngamma"

external gamma_rates: float -> float -> int ->
    (float,Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array1.t = 
        "gamma_CAML_rates"

external rand_normal : float -> float -> float = "gamma_CAML_randnormal"

external rand_exp    : float -> float = "gamma_CAML_randexp"

external rand_gamma  : float -> float -> float = "gamma_CAML_randgamma"

(** {6 Statistics Object} *)

(** An object that keeps track of running variance, mean, and standard deviation
    of a simulation. **)
class running_stats = object(self)

    val mutable   k = 0
    val mutable m_k = 0.0
    val mutable s_k = 0.0

    val mutable max_x = 0.0
    val mutable min_x = max_float

    method push x_k =
        k <- k+1;
        max_x <- max max_x x_k;
        min_x <- min min_x x_k;
        if k = 1 then begin
            m_k <- x_k;
            s_k <- 0.0;
        end else begin
            let m_k1 = m_k +. ((x_k -. m_k) /. (float_of_int k)) in
            let s_k1 = s_k +. ((x_k -. m_k1) *. (x_k -. m_k)) in
            s_k <- s_k1;
            m_k <- m_k1;
        end

    method iter ()     = k
    method min ()      = min_x
    method max ()      = max_x
    method variance () = if k > 1 then s_k /. (float_of_int (k-1)) else 0.0
    method mean ()     = if k = 0 then 0.0 else m_k
    method std_dev ()  = sqrt (self#variance ())

end


(** {6 Infix Module} *)
module type I =
    sig
        val set_eps : float -> unit
        val set_ops : int   -> unit
        val reset   : unit  -> unit

        val (=.) : float -> float -> bool
        val (>.) : float -> float -> bool
        val (<.) : float -> float -> bool
    end

module FPInfix : I =
    struct

        let default = tolerance

        let warn =
            "Numerical.Infix; Tolerance for floating point equality is being"
            ^" increased. This could be due to the length of your sequences or"
            ^" any situation where a large number of floating point processes"
            ^" build up error."

        let l_eps     = ref default

        let reset ()  = l_eps := default

        let set_eps n_eps =
            if n_eps > default then begin
                Status.user_message Status.Warning warn;
                l_eps := n_eps
            end else begin
                l_eps := default
            end

        let set_ops i = 
            set_eps (Pervasives.epsilon_float *. (float_of_int i))

        let (=.) a b = (abs_float (a-.b)) < !l_eps 
            (* (a = b) || (* basic case; will short-circuit *)
            (if a *. b = 0.0
                (* one is zero, thus relative error is misleading *)
                then (a -. b) < (!l_eps *. !l_eps)
                (* Here it is and use relative error *)
                else ((a -. b) /. ((abs_float a) +. (abs_float b))) < !l_eps) *)
            
        let (>.) a b = (a > b)

        let (<.) a b = (a < b)

    end

open FPInfix


(** {6 Types} *)

(** Simplex strategy defines how to expand, contract, and reflect a simplex *)
type simplex_strategy = 
    {   alpha : float;  (** The Reflection factor *)
         beta : float;  (** The Contraction factor *)
        gamma : float;  (** The Expansion factor *)
        delta : float;  (** The Shrinkage (Massive Contraction) factor *)
    } 

(** The Subplex strategy contains a simplex strategy and added features *)
type subplex_strategy = 
    { simplex : simplex_strategy;   (** How to perform the simplex *)
          psi : float;              (** The Simplex Reduction coefficient *)
        omega : float;              (** The Step Reduction coefficient *)
        nsmin : int;                (** Minimum subspace dimension; or 2 *)
        nsmax : int;                (** Maximum subspace dimension; or 5 *)
    }

(** Define an optimization strategy. This will call the appropriate optimization
    routine with the specified parameters. A wrapper around the four methods we
    have to optimize multi-dimensional functions. *)
type optimization_strategy = 
    {   routine : [ `Simplex of simplex_strategy option
                  | `Subplex of subplex_strategy option
                  | `BFGS    of float option
                  | `Brent_Multi ];
        (* Below are optional to use the algorithm default *)
        max_iter : int option;
        tol      : float option; 
    }

(** Define levels of optimization for the chooser *)
type opt_modes =
    [ `None
    | `Coarse of int option
    | `Exhaustive of int option
    | `Custom of optimization_strategy list ]

(** Here we define a simplex as a collection of points with attached data -the
    tree or some added information carried through the computation. The array
    should be the size of the dimension plus one. *)
type 'a simplex = (float array * ('a * float)) array

(** Set the brents tolerance for the C-Brents method in optimization branches *)
external set_tol_c_brents_method : float -> unit = "likelihood_CAML_set_tol"

(** Return the default tolerance set for the opt_mode **)
let get_tol = function
    | `None         -> max_float
    | `Exhaustive _ -> tolerance
    | `Coarse _     -> sqrt tolerance
    | `Custom _     -> assert false

(** get the brents tolerance for the C-Brents method in optimization branches *)
external get_tol_c_brents_method : unit -> float = "likelihood_CAML_get_tol"

(** Default Simplex Strategy defined by NMS *)
let default_simplex =
    { alpha = 1.0; beta = 0.5; gamma = 2.0; delta = 0.5; }

(** Default Simplex Strategy defined by FSAoNA paper *)
let default_subplex =
    { simplex = default_simplex;
      omega = 0.1; psi = 0.25; nsmin = 2; nsmax = 5; }

(** Define the default strategy from a specific routine *)
let default_strategy ?tol ?iter routine =
    let routine = match routine with
        | `Simplex      -> `Simplex None
        | `Subplex      -> `Subplex None
        | `Brent_Multi  -> `Brent_Multi
        | `BFGS         -> `BFGS None
    in
    { routine = routine; max_iter = iter; tol = tol; }


(** {6 Numerical Functions for Optimization Routines *)

(** find the derivative of a single variable function *)
let derivative_at_x ?(epsilon=epsilon) f x fx =
    let _,f_new = f (x +. epsilon) in
    (f_new -. fx) /. epsilon

(** find the magnitude of a vector x_array *)
let magnitude x_array = 
    sqrt (Array.fold_left (fun acc x -> acc +. (x *. x)) 0.00 x_array)

(** find the gradient of a multi-variant function at a point x_array *)
let gradient_at_x ?(epsilon=epsilon) f_ x_array f_array : float array =
    let i_replace i x v = let y = Array.copy x in Array.set y i v; y in
    Array.mapi
        (fun i i_val ->
            derivative_at_x
                    ~epsilon
                    (fun x ->
                        let newvec = i_replace i x_array x in
                        let newlk = f_ newvec in
                        debug_printf "\t[%s] -- %f (d:%F)\n" (pp_farray newvec)
                                     (snd newlk) ((snd newlk)-.f_array);
                        newlk)
                    i_val
                    f_array)
        x_array

(** dot product of two arrays *)
let dot_product x_array y_array = 
    let n = Array.length x_array and r = ref 0.0 in
    assert (n = Array.length y_array);
    for i = 0 to n-1 do
        r := !r +. (x_array.(i) *. y_array.(i));
    done;
    !r

(** map a matrix with a function *)
let matrix_map f mat = 
    let n1 = Array.length mat and n2 = Array.length mat.(0) in
    let output_matrix = Array.create_matrix n1 n2 0.0 in
    for i = 0 to n1 - 1 do for j = 0 to n2 -1 do 
        output_matrix.(i).(j) <- f i j mat.(i).(j);
    done; done;
    output_matrix

(** Calculates the infinity-norm in L^p space. The Maximum Norm. *)
let inf_norm_vec x =
    Array.fold_left (max) (~-.max_float) x

(** Calculates the 2-norm in L^p space. The Euclidean Norm. *)
let two_norm_vec x =
    sqrt (Array.fold_left (fun acc x -> acc +. (x *. x)) 0.0 x)

(** Calculates the 1-norm in L^p space. The Taxicab Norm. *)
let one_norm_vec x =
    Array.fold_left (fun acc x -> acc +. (abs_float x)) 0.0 x

(** subtract one vector from another *)
let sub_vec x y = 
    Array.mapi (fun i _ -> x.(i) -. y.(i)) x

(** Add two vectors together; Imperative style *)
let add_veci x y = 
    for i = 0 to (Array.length x)-1 do
        x.(i) <- x.(i) +. y.(i);
    done;
    ()


(** {6 Numerical Optimization Functions *)

(** Uses a combination of golden section searches and parabolic fits to find the
    optimal value of a function of one variable. **)
let brents_method ?(max_iter=200) ?(v_min=minimum) ?(v_max=300.0)
                  ?(tol=tolerance) ?(epsilon=epsilon) f ((v_orig,f_orig) as orig) =
    debug_printf "Starting Brents Method max_iter=%d, tol=%f, epsilon=%f\n%!" max_iter tol epsilon;
  (*-- ensure value falls between range; if using one *)
    let minmax value = max (min v_max value) v_min in
  (*-- constant for the golden ratio *)
    let golden = 0.3819660 in
  (*-- approximation of equality; based on optional argument above (over-ride
    local fuzzy equality to allow more control over the particular run). *)
    FPInfix.set_eps epsilon;
  (*-- a function that copies the sign of the second argument to the first *)
    let sign a b = if b > 0.0 then abs_float a else ~-. (abs_float a) in
  (*-- auxillary functions to bracket a region *)
    let rec create_initial_three_and_bracket (o,(_,fo) as o') =
        let rec create_scaled (v,_) s = 
            let vs = minmax (v *. s) in let fvs = f vs in
            debug_printf "Calculated [%f,%f] in Brent\n%!" vs (snd fvs);
            (vs,fvs)
        and push_left a b c = (create_scaled a 0.5),a,b
        and push_right a b c = b,c,(create_scaled c 2.0)
        and bracket_region ((l,(_,fl)) as low) ((m,(_,fm)) as med) ((h,(_,fh)) as hi) =
            if l =. h then (low,med,hi)       (* converged *)
            else if fl =. fm && fm =. fh then (* converged; flat  *)
                begin
                debug_printf ("Cannot bracket flat region for brents method;"^^
                              "[%f,%f] [%f,%f] [%f,%f]\n%!") l fl m fm h fh;
                (low,med,hi)
                end
            else if fl <= fm && fm <= fh then (* increasing *)
                let a,b,c = push_left low med hi in bracket_region a b c
            else if fl >= fm && fm >= fh then (* decreasing *)
                let a,b,c = push_right low med hi in bracket_region a b c
            else if fm <= fl && fm <= fh then (* a bracket! *) (low,med,hi)
            else begin (* bracketed a maximum... wut? *)
                (* let us do something gracefully, push ourselves to the minimum
                 * on the left or right, and continue with the algorithm;
                 * priority pushes ourselves to smaller branches. *)
                warning_message "Numerical.brent; Cannot bracket a region.";
                if fl <= fh 
                    then let a,b,c = push_left  low med hi in bracket_region a b c
                    else let a,b,c = push_right low med hi in bracket_region a b c
            end
        in
        debug_printf "Trying to bracket around %f,%f\n%!" o fo;
        bracket_region (create_scaled o' 0.2) o' (create_scaled o' 2.0)
  (*-- brents method as in Numerical Recipe in C; 10.2 *)
    and brent ((x,(_,fx)) as x') ((w,(_,fw)) as w') ( v') a b d e iters pu =
        let (v,(_,fv))  = v' in
        debug_printf "Iteration %d, bracketing (%f,%f) with: %f,%f,%f\n%!" 
                        iters a b x w v;
        let xm = (a +. b) /. 2.0
        and tol1 = tol *. (abs_float x) +. epsilon in
        (* check ending conditions *)
        if iters > max_iter then begin 
            warning_message "Numerical.brent; hit max number of iterations.";
            x'
        end else if (abs_float (x-.xm)) <= ((2.0 *. tol) -. (b -. a) *. 0.5) then x'
        else begin
            let d,e =
                if (abs_float e) > tol1 then begin
                    (* calculate the abscissa *)
                    let r = (x -. w) *. (fx -. fv) 
                    and q = (x -. v) *. (fx -. fw) in
                    let p = ((x -. v) *. q) -. ((x-.w) *. r) in
                    let q = 2.0 *. (q -. r) in
                    let p = if q > 0.0 then ~-. p else p in
                    let q = abs_float q in
                    (* the acceptability of the parabolic fit? *)
                    if (abs_float p) >= (abs_float (0.5 *. q *. e))
                        || p <= q *. (a -. x) || p >= q *. (b -. x) then
                        (* do a golden section instead of parabolic fit *)
                        let e = if x >= xm then a-.x else b -. x in
                        let d = golden *. e in
                        debug_printf "\tDoing a golden section search1: %f//%f\n" d e;
                        d,e
                    else begin
                        (* take the parabolic step *)
                        let d_new = p /. q in
                        let u = x +. d_new in
                        let d,e = 
                            if (u -. a) < (tol1 *. 2.0) || (b -. u) < (tol1 *. 2.0) 
                                then sign tol1 (xm -. x),d
                                else d_new,d
                        in
                        debug_printf "\tDoing a parabolic fit: %f//%f\n" d e;
                        d,e
                    end
                end else begin
                    let e = if x >= xm then a -. x else b -. x in
                    let d = golden *. e in
                    debug_printf "\tDoing a golden section search2: %f//%f\n" d e;
                    d,e
                end
            in
            (* the ONLY function evalution for each iteration *)
            let u =
                if (abs_float d) >= tol1
                    then minmax (x+.d)
                    else minmax (x+.(sign tol1 d))
            in
            let fu = f u in
            let u',fu = (u,fu), snd fu in
            debug_printf "\tCalculated [%f,%.11f] in Brent\n%!" u fu;
            (* what to do with results for next iteration *)
            if pu = u then
                begin if fx < fu then x' else u' end
            else if fu <= fx then begin
                let a,b = if u >= x then x,b else a,x in
                brent u' x' w' a b d e (iters+1) u
            end else begin
                let a,b = if u < x then u,b else a,u in
                if fu <= fw || w =. x then
                    brent x' u' w' a b d e (iters+1) u
                else 
                    brent x' w' u' a b d e (iters+1) u
            end
        end
    in
    let (lv,_),m,(hv,_) = create_initial_three_and_bracket orig in
    cdebug_printf "Bracketed by Brents Method : (%f,%f,%f)\n%!" lv (fst m) hv;
    let (b,(_,fb)) as res = brent m m m lv hv 0.0 0.0 0 (fst m) in
    cdebug_printf "Iterated Brents Method from (%f,%f) to (%f,%f)\n%!"
                    v_orig (snd f_orig) b fb;
    FPInfix.reset ();
    res


(** Meta function above; we sequentially modify each variable ONCE; RAxML *)
let brents_method_multi ?(max_iter=200) ?(v_min=minimum) ?(v_max=300.0)
                        ?(tol=tolerance) ?(epsilon=epsilon) f orig =
    let rec do_single i ((array,x) as data) =
        debug_printf "Optimizing %d in %s\n%!" i (pp_farray array);
        if i < (Array.length array) then
            let (v,fv) =
                brents_method ~max_iter ~v_min ~v_max ~tol ~epsilon
                          (update_single i array) (array.(i),x)
            in
            let () = Array.set array i v in
            do_single (i+1) (array,fv)
        else
            data
    and update_single i array v =
        let array = Array.copy array in
        let () = Array.set array i v in
        f array
    in
    do_single 0 orig


(** line search along a specified direction; Numerical Recipes in C : 9.7 *)
let line_search ?(epsilon=tolerance) f point fpoint gradient maxstep direction =
    (* as in the previous function; we over-ride the local equality function for
       better control over the estimation process *)
    FPInfix.set_eps epsilon;
    let get_score x = snd x in
    (* set up some globals for the function to avoid tons of arguments *)
    let n = Array.length point and origfpoint = get_score fpoint in
    (* scale direction so, |pstep| <= maxstep *)
    let setup_function point direction gradient = 
        let direction = (* ||dir|| <= maxstep *)
            let magstep = magnitude direction in
            if magstep > maxstep
                then Array.map (fun x -> x *. maxstep /. magstep) direction 
                else direction
        in
        let slope = dot_product direction gradient
        and minstep = 
            let r = ref 0.0 in
            Array.iteri 
                (fun i x -> 
                    let tmp = (abs_float direction.(i)) /. (max (abs_float x) 1.0) in
                    if tmp > !r then r := tmp)
                point;
            epsilon /. (!r)
        and step = 1.0 in
        direction, slope, minstep, step
    (* find a new point and *)
    and next_step step prevstep slope newfpoint prevfpoint = 
        let newstep =
            if step =. 1.0 then
                ~-. slope /. (2.0 *. (newfpoint -. origfpoint -.  slope))
            else begin
                let tstep =
                    let rhs1 = newfpoint -. origfpoint -. (step *. slope)
                    and rhs2 = prevfpoint -. origfpoint -. (prevstep *. slope) in
                    let rhs1divstepstep = rhs1 /. (step *. step)
                    and rhs2divpsteppstep = rhs2 /. (prevstep *. prevstep) in
                    let a =  rhs1divstepstep -. rhs2divpsteppstep
                    and b = (step *. rhs2divpsteppstep) -. (prevstep *. rhs1divstepstep) in
                    let a = a /. (step -. prevstep) and b = b /. (step -. prevstep) in
                    if a =. 0.0 then ~-. slope /. (2.0 *. b)
                    else begin
                        let disc = (b *. b) -. (3.0 *. a *. slope) in
                        if disc < 0.0 then 0.5 *. step
                        else if b <= 0.0 then (~-. b +. (sqrt disc)) /. (3.0 *. a)
                        else ~-. slope /. (b +. (sqrt disc))
                    end
                in
                min tstep (0.5 *. step)
            end
        in
        max newstep (0.1 *. step)
    in
    (* main algorithm -- first instance sets up some variables *)
    let rec main_ prevfpoint slope direction step prevstep minstep = 
        if step < minstep then
            (point,fpoint,true) (* -~verify convergence~- *)
        else begin
            let newpoint = Array.init n (fun i -> abs_float (point.(i) +. (step *. direction.(i)))) in
            let newfpoint = f newpoint in
            if (get_score newfpoint) <= origfpoint then begin
                debug_printf "\t\t%f--Accepting %f\n" prevstep (get_score newfpoint);
                (newpoint,newfpoint,false)
            end else begin
                debug_printf "\t\t%f--Rejecting %f\n" step (get_score newfpoint);
                let newstep = next_step step prevstep slope (get_score newfpoint) prevfpoint in
                main_ (get_score newfpoint) slope direction newstep step minstep
            end
        end
    in
    (* initialize and run... *)
    let direction, slope, minstep, step = setup_function point direction gradient in
    (* this could happen if the delta for gradient is huge (ie, errors in rediagnose) 
     * or some major instability in the tree/algorithm. The function will continue, 
     * but this warning message should report that the results are questionable. *)
    if (abs_float slope) > 1_000_000.0 then
        warning_message "Numerical.linesearch; Very large slope in optimization function.";
    cdebug_printf "\tInitial LineSearch: %f, slope: %f, step: %f, direction: %s\n"
                 origfpoint slope maxstep (pp_farray direction);
    let results = main_ origfpoint slope direction step step minstep in
    FPInfix.reset ();
    results


(** BFGS Algorithm; Gradient Search Function; Numerical Recipes in C : 10.7 *)
let bfgs_method ?(max_iter=200) ?(epsilon=epsilon) ?(max_step=10.0) ?(tol=tolerance) f (p,fp) =
   let n = Array.length p and get_score x = snd x in
    (* test convergence of a point --that it equals the direction, essentially *)
    let converged_l direction test_array =
        let test = ref 0.0 in
        Array.iteri
            (fun i x->
                let temp = max (abs_float x) 1.0 in
                test := max ((abs_float direction.(i)) /. temp) !test)
            test_array;
        (!test < (epsilon *. 4.0))
    (* Test tolerance for zeroing the gradient *)
    and converged_g fp gradient test_array =
        let test = ref 0.0
        and denom = max fp 1.0 in
        Array.iteri
            (fun i x ->
                let temp = (max (abs_float x) 1.0) /. denom in
                let temp = temp *. (abs_float gradient.(i)) in
                if temp > !test then test := temp)
            test_array;
        (!test < tol)
    (* Setup initial hessian (identity), initial gradiant vector, maximum step and direction *)
    and setup_function f_array x_array fx_array =
        let hessian =
            let h = Array.make_matrix n n 0.0 in
            for i = 0 to n-1 do h.(i).(i)  <- 1.0 done;
            h
        and x_grad = gradient_at_x ~epsilon f_array x_array fx_array in
        let dir = Array.map (fun x -> ~-. x) x_grad
        and mxstep = max_step *. (max (magnitude x_array) (float_of_int n)) in
        hessian, x_grad, mxstep, dir in
    (* Do a line search step --return new p, new fp, new dir, if converged *)
    let line_searcher f p fp gradient step dir =
        let np,nfp,_ = line_search f p fp gradient step dir in
        let dir = Array.init n (fun i -> np.(i) -. p.(i) ) in
        np, nfp, dir, (converged_l dir np)
    (* update gradient --ret new gradient, difference of gradients,
     * difference of gradient times hessian matrix, if converged *)
    and gradient_update hessian ograd f p fp = 
        let ngrad = gradient_at_x ~epsilon f p fp in
        let dgrad = Array.init n (fun i -> ngrad.(i) -. ograd.(i)) in
        let hgrad = 
            Array.init n
                (fun i -> 
                    let res = ref 0.0 in
                    for j = 0 to n-1 do
                        res := !res +. (hessian.(i).(j) *. dgrad.(j));
                    done; !res)
        in
        ngrad, dgrad, hgrad, (converged_g fp ngrad p)
    (* Update the hessian matrix --skips the update if fac not sufficiently positive *)
    and bfgs_update_matrix dgrad hgrad direc hessian =
        let fac = dot_product dgrad direc
        and fae = dot_product dgrad hgrad
        and sumdgr = Array.fold_left (fun a x -> (x *. x) +. a) 0.0 dgrad
        and sumdir = Array.fold_left (fun a x -> (x *. x) +. a) 0.0 direc in
        if (fac *. fac) <= (epsilon *. sumdgr *. sumdir) then
            hessian
        else begin
            let fac = 1.0 /. fac and fad = 1.0 /. fae in
            let dgrad = 
                Array.init n
                    (fun i -> (fac *. direc.(i)) -. (fad *. hgrad.(i)))
            in
            matrix_map
                (fun i j x -> x +. (fac *. direc.(i) *. direc.(j))
                                -. (fad *. hgrad.(i) *. hgrad.(j))
                                +. (fae *. dgrad.(i) *. dgrad.(j)) )
                hessian
        end
    (* calculate the new direction for the line search by the hessian matrix *)
    and calculate_direction (hessian:float array array) (gradient:float array) : float array =
        Array.init n
            (fun i -> 
                let acc = ref 0.0 in
                for j = 0 to n-1 do
                    acc := !acc -. (hessian.(i).(j) *. gradient.(j))
                done;
                !acc) in
    (* main loop of algorithm *)
    let iter = ref 0 in
    let rec main_loop hessian f p fp step direction grad =
        incr iter;
        let np, nfp, direction, c = line_searcher f p fp grad step direction in
        if c then np,nfp 
        else if (!iter > max_iter) then begin
            warning_message "Numerical.bfgs; hit max number of iterations.";
            np,nfp
        end else begin
            let grad, dgrad, hgrad, c = gradient_update hessian grad f np (get_score nfp) in
            debug_printf "\tNext Gradient: %s\n%!" (pp_farray grad);
            if c then begin
                np,nfp
            end else begin
                let hessian = bfgs_update_matrix dgrad hgrad direction hessian in
                let direction = calculate_direction hessian grad in
                debug_printf "\tNew Direction: %s\n%!" (pp_farray direction);
                main_loop hessian f np nfp step direction grad
            end
        end in
    (* initiate algorithm *)
    let hessian, pgrad, mxstep, dir = setup_function f p (get_score fp) in
    cdebug_printf "Initial Gradient: %s\n%!" (pp_farray pgrad);
    cdebug_printf "Initial Direction: %s\n%!" (pp_farray dir);
    let pf,fpf = main_loop hessian f p fp mxstep dir pgrad in
    cdebug_printf "Performed BFGS:\n\t(%s,%f)\n\t\t--[%d]-->\n\t(%s,%f)\n%!"
            (pp_farray p) (get_score fp) !iter (pp_farray pf) (get_score fpf);
    (pf,fpf)


(** {6 Simplex / Subplex specific functions **)

(** Takes a function [f] for an array [ray] ,and a subpsace [sub], like (0,2),
    that replaces the elements of the new array into [ray] according to the
    subspace association. We keep this function functional by copying the array
    each time, although it may not be necessary. *)
let function_of_subspace that_shit ray sub_assoc =
    (fun sub ->
        let cray = Array.copy ray in
        let () = Array.iteri (fun i x -> cray.(x) <- sub.(i) ) sub_assoc in
        that_shit cray)

(** Return the vector that represents the subspace *)
let make_subspace_vector subs x =
    Array.init (Array.length subs) (fun i -> x.( subs.(i) ))

(** Replace elements of a subspace vector into the main vector *)
let replace_subspace_vector sub nx x =
    Array.iteri (fun i _ -> x.( sub.(i) ) <- nx.(i)) sub

(** A subplex routine to find the optimal step-size for the next iteration of
    the algorithm. The process is outlined in 5.3.2. Each step, the vector is
    re-scaled in proportion to how much progress was made previously. If little
    progress is made then step is reduced and if onsiderable progress is made
    then it is increased. Lower and Upper Bounds in the strategy ensure we do
    not do anything rash. *)
let find_stepsize strat nsubs x dx steps =
    let n = Array.length x
    and sign x = if x >= 0.0 then 1.0 else -1.0
    and minmax x lo hi = 
        let lo = min lo hi and hi = max lo hi in
        max lo (min x hi)
    in
    (* find the scale factor for new step *)
    let stepscale =
        if nsubs > 1 then begin
            let stepscale = (one_norm_vec dx) /. (one_norm_vec steps) in
            minmax stepscale (1.0/.strat.omega) strat.omega
        end else begin
            strat.psi
        end
    in
    (* scale step vector by stepscale *)
    let nstep = Array.map (fun x -> x *. stepscale) steps in
    (* orient step in the proper direction *)
    for i = 0 to n-1 do
        if dx.(i) = 0.0
            then nstep.(i) <- ~-. (nstep.(i))
            else nstep.(i) <- (sign dx.(i)) *. (nstep.(i))
    done;
    nstep

(** get the worst, second worst, and best element of the simplex. return the
    index so we know to replace the proper point. *)
let get_simplex_hsl (simplex: 'a simplex) : int * int * int =
    let get_cost (_,(_,x)) = x in
    Array.sort (fun x y -> compare (get_cost y) (get_cost x)) simplex;
    (0, 1, ((Array.length simplex)-1))

(** Determine the subspaces by a randomization of the vector, and randomizing
    the size of the subspaces; conditions will match the strategy *)
let rand_subspace strat vec : int array list =
    let randomize ar = 
        let l = Array.length ar - 1 in
        for i = 0 to l do
            let rnd = i + Random.int (l - i + 1) in
            let tmp = ar.(i) in
            ar.(i) <- ar.(rnd);
            ar.(rnd) <- tmp;
        done;
        ar
    in
    let n = Array.length vec in
    let rec take acc n lst = match lst with
        | _ when n = 0 -> List.rev acc,lst
        | []           -> assert false
        | x::lst       -> take (x::acc) (n-1) lst
    in
    let rec continually_take taken acc lst =
        if taken = n then
            List.rev acc
        else if (n - taken) < 2 * strat.nsmin then
            continually_take (n) (lst::acc) []
        else begin
            let lo = min strat.nsmax (n-taken) in
            let r = strat.nsmin + (Random.int  (lo - strat.nsmin)) in
            if n-(r+taken) < strat.nsmin
                then continually_take taken acc lst
                else begin
                    let f,l = take [] r lst in
                    continually_take (taken+r) (f::acc) l
                end
        end
    in
    vec --> Array.mapi (fun i _ -> i)
        --> randomize
        --> Array.to_list
        --> continually_take 0 []
        --> List.map (Array.of_list)

(** A subplex routine that splits up an delta array into subspaces that match
    the criteria found in the subplex paper section 5.3.3 *)
let find_subspace strat vec =
    let n = Array.length vec in
    (* resolve a specific value of k *)
    let rec resolve_k k lvec : float =
        let rec left (acc:float) (i:int) lvec : float  = match lvec with
            | []            -> acc /. (float_of_int k)
            | xs when i = k -> right (acc /. (float_of_int k)) 0.0 xs
            | (_,x)::tl     -> left (acc +. (abs_float x)) (i+1) tl
        and right leftval acc = function
            | []       -> leftval -. (acc /. (float_of_int (n - k)))
            | (_,x)::t -> right leftval (acc +. (abs_float x)) t
        in
        left 0.0 0 lvec
    and apply_ks nleft k lvec : (int * float) list =
        if nleft < k then []
                 else (k,resolve_k k lvec)::(apply_ks nleft (k+1) lvec)
    and take acc n lst = match lst with
        | _ when n = 0 -> List.rev acc,lst
        | []           -> assert false
        | x::lst       -> take (x::acc) (n-1) lst
    in
    (* return a list of lengths for subspaces *)
    let rec partition acc nsubs nused nleft lvec : int list =
        let sorted_ks = 
            List.sort (fun (_,kv1) (_,kv2) -> compare kv2 kv1) (apply_ks nleft 1 lvec)
        in
        List.iter (fun (k,f) -> debug_printf "%d - %f\n%!" k f) sorted_ks;
        let constraint_1 k = (* fall in appropriate range *)
            (strat.nsmin <= k) && (k <= strat.nsmax)
        and constraint_2 k = (* can be partitioned further *)
            let r =
                strat.nsmin * (int_of_float
                    (ceil ((float_of_int (n-nused-k)) /. (float_of_int strat.nsmax)))) 
            in
            r <= (n-nused-k)
        in
        let rec get_next_k = function
            | (k,_)::_ when (constraint_1 k) && (constraint_2 k) ->
                if (nused+k) = n then 
                    List.rev (k::acc)
                else begin
                    debug_printf "Found K %d\n%!" k;
                    partition (k::acc) (nsubs+1) (nused+k) 
                              (nleft-k) (snd (take [] k lvec))
                end
            | _::tl -> get_next_k tl
            | []    -> assert false
        in
        get_next_k sorted_ks
    (* take a list of size of subspaces and  build association vectors for
        subspaces *)
    and build_association_arrays lvec = function
        | []    -> []
        | h::tl ->
            let this,oth = take [] h lvec in
            this::(build_association_arrays oth tl)
    in
    let lvec = 
        vec --> Array.mapi (fun i x -> (i,x))
            --> Array.to_list
            --> List.sort (fun (_,x) (_,y) -> compare (abs_float y) (abs_float x))
    in
    lvec --> partition [] 0 0 n
         --> build_association_arrays lvec
         --> List.map (List.map fst)
         --> List.map (Array.of_list)
 
(** Define how termination of the algorithm should be done. This is outlined in
    the paper, section 5.3.4, This test, because of a noisy function, uses the
    distance between the vertices of the simplex to see if the function has
    converged. *)
let subplex_termination strat tol dx x stp =
(*    let numr = max (inf_norm_vec dx) ((inf_norm_vec stp) *. strat.psi)*)
(*    and denm = max (inf_norm_vec x) 1.0 in*)
    let ret = ref false in
    Array.iteri
        (fun i _ ->
            let numr = max dx.(i) (abs_float (stp.(i) *. strat.psi))
            and denm = max (abs_float x.(i)) 1.0 in
            ret := !ret || ((numr /. denm) > tol);
            debug_printf "\tTermination: (%f/%f) = %f <= %f?\n%!"
                            numr denm (numr /. denm) tol)
        x;
    not (!ret)

(** General Simplex termination; this is done through the standard deviation of
    the simplex. *)
let simplex_termination_stddev tol simplex = 
    let n_plus_one = float_of_int (Array.length simplex) in
    let mean = 
        let s = Array.fold_left (fun acc (_,(_,x)) -> acc +. x) 0.0 simplex in
        s /. n_plus_one
    in
    let std_dev =
        Array.fold_left
            (fun acc (_,(_,x)) -> let y = x -. mean in acc +. (y *. y)) 0.0 simplex
    in
    let std_dev = (sqrt (std_dev /. (n_plus_one))) in
    std_dev < tol

(** Simplex termination test defined by Gill, Murray and Wright. This method
    looks to see that the points are in a stationary position. This is more
    appropriate for optimizing smooth functions. It can also be used for noisy
    functions, but would be equivlent to convergance of the simplex. *)
let simplex_termination_stationary tol simplex = 
    let high,low = 
        Array.fold_left
            (fun (hi,lo) (_,(_,x)) -> (max hi x),(min lo x))
            (~-.max_float, max_float)
            (simplex)
    in
    ((high-.low) /. (1.0 +. (abs_float low))) < tol

(** Calculate the centroid of a simplex. Defined by the mean of the value,
    excluding the highest (worst) point. *)
let centroid (simplex:'a simplex) h_i =
    let n = Array.length (fst simplex.(0)) in
    let c_array = Array.make n 0.0 in
    Array.iteri
        (fun s_i (r,_) ->
            if s_i = h_i then () else add_veci c_array r)
        simplex;
    for i = 0 to n-1 do
        c_array.(i) <- c_array.(i) /. (float_of_int n);
    done;
    c_array

(** Create a simplex point. All of the operations on a simplex are linear
    equations, and can be generalized to this function with different
    coefficients passed to it.
        n = x + coef (x - y)
    in each particular situation, the following
        reflection  - x = centroid, y = high point, coef =  alpha
        contraction - x = centroid, y = high point, coef = -beta
        expansion   - x = centroid, y = reflection, coef = -gamma
        shrink      - x = high point, y = all,      coef = -delta *)
let create_new_point f t strategy xvec yvec : float array * ('a * float) =
    let coef = match t with
        | `Reflection   -> strategy.alpha
        | `Contraction  -> ~-. (strategy.beta)
        | `Expansion    -> ~-. (strategy.gamma)
        | `Shrink       -> ~-. (strategy.delta)
    in
    let ret = Array.copy xvec in
    for i = 0 to (Array.length xvec) - 1 do
        ret.(i) <- xvec.(i) +. (coef *. (xvec.(i) -. yvec.(i)))
    done;
    let fret = f ret in
    let () = match t with
        | `Shrink       ->
            debug_printf "Shrinking %s -> %s -> %s:%f\n%!"
                (pp_farray xvec) (pp_farray yvec) (pp_farray ret) (snd fret);
        | `Reflection   ->
            debug_printf "Reflection %s -> %s -> %s:%f\n%!"
                (pp_farray xvec) (pp_farray yvec) (pp_farray ret) (snd fret);
        | `Contraction  ->
            debug_printf "Contraction %s -> %s -> %s:%f\n%!"
                (pp_farray xvec) (pp_farray yvec) (pp_farray ret) (snd fret);
        | `Expansion    ->
            debug_printf "Expansion %s -> %s -> %s:%f\n%!"
                (pp_farray xvec) (pp_farray yvec) (pp_farray ret) (snd fret);
    in
    ret,fret

(** Set up the initial simplex by randomly selecting points *)
let random_simplex f (p,fp) step =
    Array.init ((Array.length p)+1)
        (fun i -> 
            let p = Array.init (Array.length p) (fun _ -> Random.float 1.0) in
            p,f p)

(** Set up the initial simplex from a point and a step size *)
let initial_simplex f (p,fp) (step : float array option) =
    let step : float array = match step with
        | Some step -> step
        | None -> Array.init (Array.length p) (fun _ -> Random.float 5.0)
    in
    let simplex =
        Array.init
            ((Array.length p)+1)
            (fun i ->
                if i = 0 then (p,fp)
                         else let x = Array.copy p in
                              x.(i-1) <- x.(i-1) +. step.(i-1);
                              x, f x)
    in
    let get_cost (_,(_,x)) = x in
    Array.sort (fun x y -> compare (get_cost y) (get_cost x)) simplex;
    simplex

(* shrink involves modifying each point except the best *)
let shrink_simplex simplex f strategy i_l =
    let replace_simplex = Array.set in
    for i = 0 to (Array.length simplex)-1 do
        if i_l = i then ()
        else
            let s_i = create_new_point f `Shrink strategy
                                (fst simplex.(i_l)) (fst simplex.(i)) in
            replace_simplex simplex i s_i
    done;
    ()

(** Verify that the strategy elements are consistent with their intended
    transformation. Thus, avoids negative and fractional values in some. *)
let verify_strategy strat =
    (strat.alpha > 0.0) &&
    (strat.beta > 0.0 && strat.beta < 1.0) &&
    (strat.gamma > 0.0 && strat.gamma > strat.alpha) &&
    (strat.delta > 0.0 && strat.delta < 1.0)


(** {6 Main Simplex Algorithm} *)

(** The simplex uses an n+1 dimensional figure to move, like an amoeba, around
    the surface. The extermities of the object are evaluated, where the maximal
    values are modified, while the others are progressed with improvements. The
    four moves that can be done on the simplex are, reflection, expansion,
    contraction, and shrinkage. The degree to which these are done is modified
    by the strategy used. *)
let simplex_method ?(termination_test=simplex_termination_stddev) ?(tol=tolerance)
                   ?(simplex_strategy=default_simplex) ?(max_iter=100) ?(step=None)
                    f (p,fp) = 
    (* wrap function to keep track of the number of evaluations *)
    let i = ref 0 in
    let f = (fun x -> incr i; f x) in
    let strategy = simplex_strategy in
    assert( verify_strategy strategy );
    (* set up some alias functions to make the algorithm more readable. *)
    let get_cost (_,(_,x)) = x in
    let replace_simplex = Array.set in
    let rec simplex_loop step simplex =
        let i_h,_,i_l = get_simplex_hsl simplex in
        let s_c = centroid simplex i_h in
        (* first do a reflection *)
        let r = create_new_point f `Reflection strategy s_c (fst simplex.(i_h)) in
        if (get_cost r) < (get_cost simplex.(i_l)) then begin
            (* since it's so good, do an expansion *)
            let e = create_new_point f `Expansion strategy s_c (fst r) in
            if (get_cost e) < (get_cost simplex.(i_l))
                then replace_simplex simplex i_h e
                else replace_simplex simplex i_h r
        end else begin
        (* do a contraction of the simplex instead *)
            let c =
                (* contract from the worst point *)
                if (get_cost simplex.(i_h)) < (get_cost r)
                    then create_new_point f `Contraction strategy s_c (fst simplex.(i_h))
                    else create_new_point f `Contraction strategy s_c (fst r)
            in
            if (get_cost c) < (min (get_cost r) (get_cost simplex.(i_h)))
                (* successful contraction *)
                then replace_simplex simplex i_h c
                (* contraction failed; shrink --massive contraction *)
                else shrink_simplex simplex f strategy i_l
        end;
        if ((not (termination_test tol simplex)) || !i < max_iter)
            then simplex
            else simplex_loop (step+1) simplex
    in
    (* setup the initial simplex *)
    let simplex = simplex_loop 0 (initial_simplex f (p,fp) step) in
    let _,_,best = get_simplex_hsl simplex in
    cdebug_printf "\tSimplex Found : %s -- %f\n%!"
                (pp_farray (fst simplex.(best))) (snd (snd simplex.(best)));
    simplex.( best )


(** {6 Main Subplex Algorithm} *)

(** The Subplex method is a generalization to a number of algorithms, paramount
    the Nelder-Mead simplex method, with alternating variables, and Nelder-Mead
    Simplex with restart. The advantages are outlined in the previously
    mentioned paper, in section 5.3.6. *)
let subplex_method ?(subplex_strategy=default_subplex) ?(tol=tolerance) ?(max_iter=50) f (p,fp) =
    let i = ref 0 in
    let rec subplex_loop step subs ((x,fx) as xfx) dx =
        incr i;
        let step = find_stepsize subplex_strategy (List.length subs) x dx step in
(*        let subs = rand_subspace strategy dx in*)
        let subs = find_subspace subplex_strategy dx in
        debug_printf "Iteration %d\n%!" !i;
        debug_printf "\tStep : %s\n%!" (pp_farray step);
        debug_printf "\tBest : %s\n%!" (pp_farray x);
        debug_printf "\tSubs : %d\n%!" (List.length subs);
        List.iter (fun x -> debug_printf "\t\t%s\n%!" (pp_iarray x)) subs;
        let (nx,nfx) as nxnfx =
            let simplex_strategy = subplex_strategy.simplex in
            List.fold_left
                (fun (x,fx) sub ->
                    debug_printf "\tModifying %s\n%!" (pp_iarray sub);
                    let sub_vec = make_subspace_vector sub x in
                    let step    = Some (make_subspace_vector sub step) in
                    let (nx,nfx) =
                        simplex_method ~simplex_strategy ~step
                                    (function_of_subspace f x sub) (sub_vec,fx)
                    in
                    replace_subspace_vector sub nx x;
                    debug_printf "\t\tAV : %s -- %f\n%!" (pp_farray x) (snd nfx);
                    (x,nfx))
                xfx
                subs
        in
        debug_printf "\tSimplexAV found %s -- %f\n%!" (pp_farray nx) (snd nfx);
        let dx = sub_vec x nx in
        debug_printf "\tdx : %s\n%!" (pp_farray dx);
        if (subplex_termination subplex_strategy tol dx nx step) || (!i > max_iter)
            then nxnfx
            else subplex_loop step subs nxnfx dx
    in
    let dx = Array.make (Array.length p) 0.0 in
    let ((p,(_,fp)) as pdfp) =
        subplex_loop
            (find_stepsize subplex_strategy 1 p dx (Array.make (Array.length p) 1.0))
            [ (Array.init (Array.length p) (fun x -> x)) ]
            (p,fp)
            (dx)
    in
    cdebug_printf "Subplex Found %s -- %f\n%!" (pp_farray p) fp;
    pdfp

(** Determine the numerical optimization strategy from the methods cost mode *)
let default_numerical_optimization_strategy o p =
    let tol = Some (get_tol o) in
    let meth = if p < 3 then `Brent_Multi else `BFGS in
    match o with
    | `None         -> []
    | `Coarse     _ -> (default_strategy ?tol meth) :: []
    | `Exhaustive _ -> (default_strategy meth) :: []
    | `Custom     s -> s

(** Determine the default number of branch and model optimizations to do. These
    are related to the current optimization mode stored in Methods (usually). *)
let default_number_of_passes = function
    | `None                -> 0
    | `Coarse (Some i)
    | `Exhaustive (Some i) -> i
    | `Coarse None         -> 1
    | `Custom _
    | `Exhaustive None     -> max_int

(** Run an optimization strategy; call the proper algorithm w/ convergence
    properties; a list of the appropriate ones are included here *)
let run_method opts f pfp =
    List.fold_left
        (fun pfp opt ->
            let max_iter = opt.max_iter and tol = opt.tol in
            match opt.routine with
            | `Simplex simplex_strategy -> simplex_method ?tol ?max_iter ?simplex_strategy f pfp
            | `Subplex subplex_strategy -> subplex_method ?tol ?max_iter ?subplex_strategy f pfp
            | `BFGS max_step            -> bfgs_method ?max_iter ?tol ?max_step f pfp
            | `Brent_Multi              -> brents_method_multi ?max_iter ?tol f pfp)
        pfp
        opts

