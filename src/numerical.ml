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

let tolerance = 1e-6
let epsilon   = 1e-10
let minimum   = tolerance *. 2.0

let (=.) a b = abs_float (a-.b) < tolerance
let (-->) b a = a b

let debug = false
let failwithf format = Printf.ksprintf failwith format
let warning_message format = Printf.ksprintf (Status.user_message Status.Warning) format

let is_zero x = match classify_float x with
        | FP_zero | FP_subnormal -> true
        | FP_infinite | FP_nan | FP_normal -> false
and is_nan x = match classify_float x with
        | FP_zero | FP_subnormal
        | FP_infinite | FP_normal -> false
        | FP_nan -> true
and is_inf x = match classify_float x with
        | FP_infinite -> true
        | FP_nan | FP_zero | FP_subnormal | FP_normal -> false

let debug_printf msg format = 
    Printf.ksprintf (fun x -> if debug then print_string x; flush stdout) msg format

and pp_farray xs =
    (Array.fold_left (fun acc x -> acc^"| "^(string_of_float x)^" ") "[" xs)^" |]"

(* calculates the gamma rates for specific alpha, beta and #classes *)
external gamma_rates: float -> float -> int ->
    (float,Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array1.t = 
        "gamma_CAML_rates"

(** GENERAL BRENTS METHOD **)
let brents_method ?(max_iter=100) ?(v_min=minimum) ?(v_max=300.0)
                  ?(tol=tolerance) ?(epsilon=epsilon) ((v_orig,f_orig) as orig) f =
    debug_printf "Starting Brents Method max_iter=%d, tol=%f, epsilon=%f\n%!" max_iter tol epsilon;
  (*-- ensure value falls between range; if using one *)
    let minmax value = max (min v_max value) v_min in
  (*-- constant for the golden ratio *)
    let golden = 0.3819660 in
  (*-- approximation of equality; based on optional argument above *)
    let (=.) a b = (abs_float (a -. b)) < epsilon in
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
    and brent ((x,(_,fx)) as x') ((w,(_,fw)) as w') ((v,(_,fv)) as v') a b d e iters =
        debug_printf "Iteration %d, bracketing (%f,%f) with: %f,%f,%f\n%!" 
                        iters a b x w v;
        let xm = (a +. b) /. 2.0
        and tol1 = tol *. (abs_float x) +. epsilon in
        (* check ending conditions *)
        if iters > max_iter then begin 
            warning_message "Numerical.brent; hit max number of iterations.";
            x'
        end else if (abs_float (x-.xm)) <= ((2.0 *. tol) -. (b -. a) /. 2.0) then x'
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
                        let d,e = golden *. e, e in
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
                    let d,e = golden *. e, e in
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
            if fu <= fx then begin
                let a,b = if u >= x then x,b else a,x in
                brent u' x' w' a b d e (iters+1)
            end else begin
                let a,b = if u < x then u,b else a,u in
                if fu <= fw || w =. x then
                    brent x' u' w' a b d e (iters+1)
                else 
                    brent x' w' u' a b d e (iters+1)
            end
        end
    in
    let (lv,_),m,(hv,_) = create_initial_three_and_bracket orig in
    let (b,(_,fb)) as res = brent m m m lv hv 0.0 0.0 0 in
    debug_printf "Iterated Brents Method from (%f,%f) to (%f,%f)\n%!"
                    v_orig (snd f_orig) b fb;
    res

(* find the derivative of a single variable function *)
let derivative_at_x ?(epsilon=epsilon) f x fx =
    let _,f_new = f (x +. epsilon) in
    (f_new -. fx) /. epsilon
(* find the magnitude of a vector x_array *)
let magnitude x_array = sqrt (Array.fold_left (fun acc x -> acc +. (x *. x)) 0.00 x_array)
(* find the gradient of a multi-variant function at a point x_array *)
let gradient_at_x ?(epsilon=epsilon) f_ x_array f_array : float array = 
    let i_replace i x v = let y = Array.copy x in Array.set y i v; y in
    Array.mapi 
        (fun i i_val ->
            derivative_at_x ~epsilon 
                            (fun x -> 
                                let newvec = i_replace i x_array x in
                                let newlk = f_ newvec in
                                debug_printf "\t[%s] -- %f\n" (pp_farray newvec) (snd newlk);
                                newlk)
                            i_val
                            f_array)
        x_array
(* dot product of two arrays *)
let dot_product x_array y_array = 
    let n = Array.length x_array and r = ref 0.0 in
    assert (n = Array.length y_array);
    for i = 0 to n-1 do
        r := !r +. (x_array.(i) *. y_array.(i));
    done;
    !r
(* map a matrix with a function *)
let matrix_map f mat = 
    let n1 = Array.length mat and n2 = Array.length mat.(0) in
    let output_matrix = Array.create_matrix n1 n2 0.0 in
    for i = 0 to n1 - 1 do for j = 0 to n2 -1 do 
        output_matrix.(i).(j) <- f i j mat.(i).(j);
    done; done;
    output_matrix

(* line search along a specified direction *)
(* Numerical Recipes in C : 9.7            *)
let line_search ?(epsilon=tolerance) f point fpoint gradient maxstep direction =
    let (=.) a b = epsilon > (abs_float (a -. b)) and get_score x = snd x in
    (* set up some globals for the function to avoid tons of arguments *)
    let n = Array.length point and origfpoint = get_score fpoint in
    (* scale direction so, |pstep| <= maxstep *)
    let setup_function point direction gradient = 
        let direction = (* ||dir|| <= maxstep *)
            let magstep = magnitude direction in
            if magstep > maxstep then Array.map (fun x -> x *. maxstep /. magstep) direction 
            else direction in
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
                debug_printf "\t\t%f--Rejecting %f\n" prevstep (get_score newfpoint);
                let newstep = next_step step prevstep slope (get_score newfpoint) prevfpoint in
                main_ (get_score newfpoint) slope direction newstep step minstep
            end
        end
    in
    (* initialize and run... *)
    let direction, slope, minstep, step = setup_function point direction gradient in
    debug_printf "\tInitial LineSearch: %f, slope: %f\n" origfpoint slope;
    (* this could happen if the delta for gradient is huge (ie, errors in rediagnose) 
     * or some major instability in the tree/algorithm. The function will continue, 
     * but this warning message should report that the results are * questionable   *)
    if (abs_float slope) > 100000.0 then 
        warning_message "Numerical.linesearch; Very large slope in optimization function.";
    main_ origfpoint slope direction step step minstep 

(** BFGS Algorithm                   **)
(** Numerical Recipes in C : 10.7    **)
let bfgs_method ?(max_iter=200) ?(epsilon=epsilon) ?(mx_step=10.0) ?(g_tol=tolerance) f p fp =
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
        (!test < g_tol)
    (* Setup initial hessian (identity), initial gradiant vector, maximum step and direction *)
    and setup_function f_array x_array fx_array =
        let hessian =
            let h = Array.make_matrix n n 0.0 in
            for i = 0 to n-1 do h.(i).(i)  <- 1.0 done;
            h
        and x_grad = gradient_at_x f_array x_array fx_array in
        let dir = Array.init n (fun i -> ~-. (x_grad.(i)) )
        and mxstep = mx_step *. (max (magnitude x_array) (float_of_int n)) in
        hessian, x_grad, mxstep, dir in
    (* Do a line search step --return new p, new fp, new dir, if converged *)
    let line_searcher f p fp gradient step dir =
        let np,nfp,_ = line_search f p fp gradient step dir in
        let dir = Array.init n (fun i -> np.(i) -. p.(i) ) in
        np, nfp, dir, (converged_l dir np)
    (* update gradient --ret new gradient, difference of gradients,
     * difference of gradient times hessian matrix, if converged *)
    and gradient_update hessian ograd f p fp = 
        let ngrad = gradient_at_x f p fp in
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
    debug_printf "Initial Gradient: %s\n%!" (pp_farray pgrad);
    debug_printf "Initial Direction: %s\n%!" (pp_farray dir);
    let pf,fpf = main_loop hessian f p fp mxstep dir pgrad in
    debug_printf "Performed BFGS:\n\t(%s,%f)\n\t\t--[%d]-->\n\t(%s,%f)\n%!" (pp_farray p)
                    (get_score fp) !iter (pp_farray pf) (get_score fpf);
    (pf,fpf)

