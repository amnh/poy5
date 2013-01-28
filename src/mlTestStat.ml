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

let () = SadmanOutput.register "MlTestStat" "$Revision: 3049 $"

let (-->) a b = b a
let (-|>) a b = let () = b a in a

let debug_cdf = false
and debug_boot= false
and debug_sh  = false
and debug_au  = false
and debug     = false

let failwithf format = Printf.ksprintf failwith format

let info_user_message format =
    Printf.ksprintf (Status.user_message Status.Information) format

let warn_user_message format =
    Printf.ksprintf (Status.user_message Status.Warning) format

let error_user_message format =
    Printf.ksprintf (Status.user_message Status.Error) format

let outputt table =
    Status.output_table Status.Information table


module type S = sig

    type a 
    type b

    type tree = (a,b) Ptree.p_tree
    type replicate = { m:bool; b:bool }

    val rell : replicate
    val full : replicate

    val can_perform_stat_tests : tree -> bool
    val analyze_tree : tree -> (float * float) * (float * float)
    val analyze_pair : tree -> tree -> float * float

    val bootstrap_data : ?n:int -> replicate -> float array -> float array
    val get_cdf : float array -> float array

    val kh : ?n:int -> ?p_star:float -> ?rep:replicate -> tree -> tree -> ()
    val sh : ?n:int -> ?p_star:float -> ?rep:replicate -> tree list -> ()
    val au : ?n:int -> ?p_star:float -> ?rep:replicate -> tree list -> ()

    val analyze : tree list -> unit
end

module Make (Node : NodeSig.S with type other_n = Node.Standard.n)
            (Edge : Edge.EdgeSig with type n = Node.n) 
            (TreeOps : Ptree.Tree_Operations with type a = Node.n with type b = Edge.e) =
struct

    type a = Node.n
    type b = Edge.e
    type tree = (a, b) Ptree.p_tree 

    type replicate = { m:bool; b:bool }

    let rell = {m=false;b=false;}
    let full = {m=true; b=true; }

    let can_perform_stat_tests tree = match (Ptree.get_roots tree) with
        | [x] ->
            let root = match x.Ptree.root_median with
                | Some (`Edge _,n)   -> n
                | Some (`Single _,n) -> n
                | None               -> assert false
            in
            root --> (fun x -> List.hd x.AllDirNode.unadjusted)
                 --> (fun x -> AllDirNode.force_val x.AllDirNode.lazy_node)
                 --> (fun x -> List.hd x.Node.characters)
                 --> (fun x -> match x with
                        | Node.StaticMl _ -> true
                        | _               -> false)
        | []  -> false
        |  _  -> false


    (* return the singular characters mlstatic data of a tree *)
    let get_ml_root_data tree =
        let root = match (List.hd (Ptree.get_roots tree)).Ptree.root_median with
            | Some (`Edge _,n)   -> n
            | Some (`Single _,n) -> n
            | None               -> assert false
        in
        root --> (fun x -> List.hd x.AllDirNode.unadjusted)
             --> (fun x -> AllDirNode.force_val x.AllDirNode.lazy_node)
             --> (fun x -> List.hd x.Node.characters)
             --> (fun x -> match x with
                    | Node.StaticMl x -> x.Node.preliminary
                    | _ -> assert false)

    (* general helper function to return cost; all these methods assume no
      negation of loglikelihood, so we transform to that standard temporarily *)
    let get_ml_cost x =
        ~-. (MlStaticCS.median_cost x)

    (* output the rell and tree variance and mean *)
    let analyze_tree tree =
        let root  = get_ml_root_data tree in
        let sm,sv = MlStaticCS.variance root in
        let rm,rv = MlStaticCS.rell_bootstrap root 10 in
        (rm,rv),(sm,sv)

    (* calculate the pair; and return string to format output later *)
    and analyze_pair r1 r2 =
        let r1 = get_ml_root_data r1 in
        let r2 = get_ml_root_data r2 in
        MlStaticCS.variance_ratio r1 r2

    (* print matrix of mean/variance information for a list of trees *)
    let analyze x =
        let n  = Array.length x in
        let m  =
            Array.init (n+1)
                (fun i ->
                    if i = 0 then
                        Array.init (n+1) 
                            (fun x -> if x = 0 then "" else string_of_int (x-1))
                    else
                        Array.init (n+1)
                            (fun j ->
                                if j = 0 then string_of_int (i-1)
                                         else analyze_pair x.(i-1) x.(j-1)))
        in
        outputt m

(*
    let rell_bootstrap_data =
        let rec idx_of_cdf (cdf:float array) (g:float) (il:int) (ih:int) : int =
            let ig = (il + ih) / 2 in
            if debug_cdf then
                Printf.printf "IDX OF CDF: %d(%f) < (%d/%f?)%f < %d(%f)\n%!"
                              il cdf.(il) ig cdf.(ig) g ih cdf.(ih);
            assert( cdf.(ih) > g );
            assert( cdf.(il) < g );
            if (cdf.(ig) < g) && (g < cdf.(ig+1))
                then ig
            else if cdf.(ig) > g
                then idx_of_cdf cdf g il ig
                else idx_of_cdf cdf g ig ih
        and bootstrap_data ?n cdf =
            let n = match n with
                | None   -> int_of_float (cdf.((Array.length cdf)-1))
                | Some x -> x
            in
            let boot_data = Array.make (Array.length cdf) 0.0 in
            for i = 0 to n do
                let r = (Random.float ((cdf.((Array.length cdf)-1))-.1.0))+.1.0 in
                let i = idx_of_cdf cdf r 0 ((Array.length cdf)-1) in
                if debug_cdf then Printf.printf "IDX OF CDF: %d = %f\n%!" i r;
                boot_data.(i) <- boot_data.(i) +. 1.0;
            done;
            boot_data
        in
        bootstrap_data

    and get_cdf t_lks =
        let cdf = Array.make (Array.length t_lks) 0.0 in
        let sum = ref 0.0 in
        for i = 0 to (Array.length t_lks)-1 do
            sum := !sum +. (fst t_lks.(i));
            cdf.(i) <- !sum;
        done;
        if debug_cdf then begin
            Printf.printf "CDF: ";
            Array.iter (fun x -> Printf.printf "%f, " x) cdf;
        end;
        cdf

    and cost_of_rell_bootstrap t_lks boot_weights : float =
        assert( (Array.length t_lks) = (Array.length boot_weights) );
        let t_cost = ref 0.0 in
        for i = 0 to (Array.length boot_weights)-1 do
            t_cost := !t_cost +. (t_lks.(i) *. boot_weights.(i));
        done;
        ~-. !t_cost

    let rell_kh ?(alpha=0.05) t1 t2 n =
        (** Calculate the cost of each tree from rell re-sampling *)
        let get_cost_of_bootstrap t1 t2 boot_weights : float =
            let t1_cost = get_cost_of_bootstrap t1 boot_weights
            and t2_cost = get_cost_of_bootstrap t2 boot_weights in
            t1_cost -. t2_cost
        (** Calculate the CDF and Site likelihoods for both trees; compression is
            data dependent, and are the same in each tree, assertion follows. *)
        and get_cdf_lks t1_root t2_root =
            let t1_lks = MlStaticCS.site_likelihood t1_root
            and t2_lks = MlStaticCS.site_likelihood t2_root in
            let cdf = get_cdf t1_lks in
            let t1_lks = Array.map snd t1_lks and t2_lks = Array.map snd t2_lks in
            cdf, t1_lks, t2_lks
        in
        (** Main Components of the algorithm.. *)
        (* 0. initial variables for computation *)
        let t1_root = get_ml_root_data t1
        and t2_root = get_ml_root_data t2 in
        let cdf, t1_lks, t2_lks = get_cdf_lks t1_root t2_root in
        (* 1. determine the test statistic *)
        let d_0 = (get_ml_cost t1_root) -. (get_ml_cost t2_root) in
        (* 2. test statistic of bootstrap replicates *)
        let d_i =
            Array.init n
                (fun _ -> get_cost_of_bootstrap t1_lks t2_lks (bootstrap_data cdf))
        in
        if debug_kh then begin
            Printf.printf "BOOTSTRAP STATISTICS: %f" d_0;
            Array.iter (fun x -> Printf.printf ", %f" x) d_i;
            print_newline ()
        end;
        (* 3. centered test-statistic for bootstrap replicates *)
        let dc_i =
            let avg = (Array.fold_left (fun a x -> a+.x) 0.0 d_i) /. (float_of_int n) in
            (Array.map (fun x -> x -. avg) d_i)
        in
        if debug_kh then begin
            Printf.printf "CENTERED BOOTSTRAP STATISTICS: %f" d_0;
            Array.iter (fun x -> Printf.printf ", %f" x) dc_i;
            print_newline ()
        end;
        (* 4. Does it pass confidence test of alpha value? *)
        ()
            

    let rell_sh ?(replicates=10000) ?(p_star:0.05) ts =
        (** STEP 0: Setup Structures and variables necessary *)
        let ts =
            ts --> List.map get_ml_root_data
               --> Array.of_list
               -|> Array.sort
                    (fun x y -> Pervasives.compare (get_ml_cost y) (get_ml_cost x))
        in
        let l = Array.map (fun x -> x, MlStaticCS.site_likelihood x) ts in
        let m = Array.length l in
        assert( m > 1 );
        let cdf = get_cdf (snd l.(0)) in
        let l = Array.map (fun (a,b) -> a,(Array.map snd b)) l in
        if debug_sh then begin
            Printf.printf "Initial Costs\n\t%!";
            Array.iter (fun (x,_) -> Printf.printf "%f, " (get_ml_cost x)) l;
            print_newline ()
        end;
        (** STEP 1: Calculate test statistic *)
        let t =
            Array.map (fun (x,_) -> (get_ml_cost (fst l.(0))) -. (get_ml_cost x)) l
        in
        if debug_sh then begin
            Printf.printf "Test Statistics\n\t%!";
            Array.iter (fun x -> Printf.printf "%f, " x) t;
            print_newline ()
        end;
        (** STEP 2: Generate bootstrap replicates *)
        let b = 
            let res = Array.make_matrix m replicates 0.0 in
            for i = 0 to replicates-1 do
                let weights = bootstrap_data cdf in
                for alpha = 0 to m-1 do
                    res.(alpha).(i) <- get_cost_of_bootstrap (snd l.(alpha)) weights
                done;
            done;
            res
        in
        if debug_sh then begin
            Printf.printf "B Matrix - Bootstrap Replicates\n%!";
            Array.iteri
                (fun alpha ib ->
                    Printf.printf "%d\t" alpha;
                    Array.iter (fun x -> Printf.printf "%f\t" x) ib;
                    let avg =
                        (Array.fold_left (fun x h -> x+.h) 0.0 ib)
                            /. (float_of_int replicates)
                    in
                    Printf.printf "\t= %f\n%!" avg)
                b
        end;
        (** STEP 3: Centering *)
        let r =
            Array.init m
                (fun alpha ->
                    let avg =
                        let s = Array.fold_left (fun x h -> x +. h) 0.0 b.(alpha) in
                        s /. (float_of_int n)
                    in
                    Array.map (fun x -> x -. avg) b.(alpha))
        in
        if debug_sh then begin
            Printf.printf "R Matrix - Centered Bootstrap Replicates\n%!";
            Array.iteri
                (fun alpha ir ->
                    Printf.printf "%d\t" alpha;
                    Array.iter (fun x -> Printf.printf "%f\t" x) ir;
                    print_newline ())
                r
        end;
        (** STEP 4: Create Replicate of test-statistic (t) *)
        let s =
            let res = Array.make_matrix m n 0.0 in
            for i = 0 to replicates-1 do
                let alpha_hat =
                    let best = ref 0 in
                    for j = 1 to m-1 do
                        if r.(!best).(i) < r.(j).(i) then best := j
                    done;
                    !best
                in
                for alpha = 0 to m-1 do
                    res.(alpha).(i) <- r.(alpha_hat).(i) -. r.(alpha).(i)
                done;
            done;
            res
        in
        if debug_sh then begin
            Printf.printf "S Matrix - Replicate Test Statistic\n%!";
            Array.iteri
                (fun alpha is ->
                    Printf.printf "%d\t" alpha;
                    Array.iter (fun x -> Printf.printf "%f\t" x) is;
                    print_newline ())
                s
        end;
        (** STEP 5: Calculate P-Values for each test-statistic *)
        let p =
            Array.init m
                (fun alpha ->
                    let p_alpha = ref 0 in
                    for i = 0 to replicates-1 do
                        if s.(alpha).(i) > t.(alpha) then incr p_alpha
                    done;
                    (float_of_int !p_alpha) /. (float_of_int replicates))
        in
        if debug_sh then begin
            Printf.printf "P-Values\n%!";
            Array.iteri
                (fun i _ ->
                    Printf.printf "%d\t%f\t%f\t%f\n%!"
                                  i (get_ml_cost (fst l.(i))) t.(i) p.(i))
                p;
        end;
        (** STEP 6: Return set such that, P_star < P_i *)
        Array.mapi (fun i _ -> ts.(i),t.(i),p.(i)) p
        

    let au ts k b =
        (** STEP 0 : Basic setup for every one of these algorithms *)
        let ts =
            ts --> List.map get_ml_root_data
               --> Array.of_list
               -|> Array.sort
                    (fun x y -> Pervasives.compare (get_ml_cost y) (get_ml_cost x))
        in
        let l = Array.map (fun x -> x, MlStaticCS.site_likelihood x) ts in
        let m = Array.length l in
        assert( m > 1 );
        let cdf = get_cdf (snd l.(0)) in
        let n   = cdf.((Array.length cdf)-1) in
        let l = Array.map (fun (a,b) -> a,(Array.map snd b)) l in
        if debug_au then begin
            Printf.printf "Initial Costs\n\t%!";
            Array.iter (fun (x,_) -> Printf.printf "%f, " (get_ml_cost x)) l;
            print_newline ()
        end;
        (** STEP 1 : Define r_k and B_k; the scale factor is uniformly distributed
            around mean 1.0,  and B are the number of replicates, equal for each
            scale factor in the algorithm. *)
        let r = Array.init k (fun i -> 0.5 +. ((float_of_int i) /. (float_of_int k))) in
        let b = Array.init k (fun _ -> b) in
        (** STEP 2.1 : Generate the bootstrap replicates of len N' = r_k*N. *)
        let y =
            Array.init k    (* k * b * t *)
                (fun i ->
                    let n' : int = int_of_float (r.(i) *. n) in
                    let scale : float = n /. (float_of_int n') in
                    Array.init b.(i)
                        (fun _ ->
                            let w = bootstrap_data ~n:n' cdf in
                            Array.map
                                (fun (_,t) -> (get_cost_of_bootstrap t w) *. scale) l))
        in
        if debug_au then begin
            for i = 0 to k-1 do
                Printf.printf "\n\n%d -- %f\n\tTs\\Bs\t" i r.(i);
                for j = 0 to (b.(i)-1) do
                    Printf.printf " % 9d  " j;
                done;
                for t = 0 to m-1 do
                    Printf.printf "\n\t%d\t" t;
                    for j = 0 to (b.(i)-1) do
                        Printf.printf "% 9.4f  " y.(i).(j).(t);
                    done;
                done;
            done;
            print_newline ()
        end;
        (** STEP 2.2 : calculate the BP-values; BP_0 = #{Y(r_k) \mem H_0}/b_k *)
        let bp =
            Array.init k
                (fun i ->
                    let cnt = ref 0 in
                    for b_i = 0 to (b.(i)-1) do
                        let best = ref true in
                        for j = 1 to m-1 do
                            best := !best && (y.(i).(b_i).(0) > y.(i).(b_i).(j));
                        done;
                        if !best then incr cnt
                    done;
                    (float_of_int !cnt) /. (float_of_int b.(i)))
        in
        if debug_au then begin
            Printf.printf "BP-Values\n%!";
            Array.iteri (fun i r_i -> Printf.printf "%d\t%f\t%f\n%!" i r_i bp.(i)) r
        end;
        (** STEP 3 : estimate d and c from weighted least squares by minimizing RSS *)
        let d,c =
            let opt_function ray =
                assert( 2 = (Array.length ray));
                let d,c = ray.(0), ray.(1) in
                let v_inv k =
                    let norm_inv = Numerical.qnorm bp.(k) in
                    let numr = bp.(k) *. (1.0 -. bp.(k))
                    and deno = Numerical.dnorm (norm_inv**2.0) in
                    (deno *. (float_of_int b.(k))) /. numr
                in
                let sum = ref 0.0 in
                for i = 0 to k-1 do
                    let rss_i = (d *. (sqrt r.(i))) +. c /. (sqrt r.(i)) in
                    let rss_i = rss_i -. (Numerical.qnorm (1.0 -. bp.(i))) in
                    let v_inv_i = v_inv i in
                    sum := !sum +. v_inv_i +. rss_i**2.0;
                done;
                (),!sum
            in
            (* the starting position below seems to be popular with the ladies *)
            let i = [| 1.0; 0.0 |],opt_function [| 1.0; 0.0 |] in
            (* below; tested brent_multi but couldn't bracket region *)
            let a,((),_) = Numerical.bfgs_method opt_function i in
            a.(0),a.(1)
        in
        if debug_au then
            Printf.printf "D:%f\nC:%f\n" d c;
        (** STEP 4 : calculate p-values, AU = 1 - \Phi(d-c) **)
        let p = Numerical.pnorm (d -. c) in
        if debug_au then
            Printf.printf "P-Value:%f\n" p;
        p

*)
end

(** Here are some tests; for these to work, we need to relax the dependencies
    and remove them from 'scripting.ml'
module IC = Make (AllDirNode.AllDirF) (Edge.LazyEdge) (AllDirChar.F)

let testing_0 () =
    POY run("test.script")
    
let testing_0_1 () =
    POY run("test2.script")

let testing_1 () =
    IC.analyze ()

let testing_2 () =
    let ts = Array.of_list (Phylo.Runtime.trees ()) in
    let res= Array.make_matrix (Array.length ts) (Array.length ts) () in
    for i = 0 to (Array.length ts)-1 do
        for j = 0 to (Array.length ts)-1 do
            if i = j then ()
                     else res.(i).(j) <- IC.rell_kh ts.(i) ts.(j) 100
        done;
    done;
    ()

let testing_3 n =
    ignore (IC.sh (Phylo.Runtime.trees ()) n)

let testing_4 r n =
    ignore(IC.au (Phylo.Runtime.trees ()) r n)
*)
