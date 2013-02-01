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

let () = SadmanOutput.register "MlTestStat" "$Revision: 3056 $"

let (-->) a b = b a
let (-|>) a b = let () = b a in a

let debug_cdf = false
and debug_boot= false

and debug_kh  = true
and debug_sh  = true
and debug_au  = true

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
    type wtree = { t : tree; slk : float array; root : MlStaticCS.t; }

    val rell : replicate
    val full : replicate

    val create_wrapped_tree : Data.bool_characters -> tree -> wtree
    val can_perform_stat_tests : Data.bool_characters -> tree -> bool

    val bootstrap_weights : ?n:int -> float array -> float array
    val replicate_cost : replicate -> wtree -> float array -> float
    val get_cdf : wtree -> float array

    val kh : ?n:int -> ?p:float -> ?rep:replicate -> ?chars:Data.bool_characters -> tree -> tree -> unit
    val sh : ?n:int -> ?p:float -> ?rep:replicate -> ?chars:Data.bool_characters -> tree list -> unit
    val au : ?n:int -> ?rep:replicate -> ?k:int -> ?chars:Data.bool_characters -> tree list -> unit

end

module Make (NodeF : NodeSig.S with type other_n = Node.Standard.n)
            (Edge : Edge.EdgeSig with type n = NodeF.n) 
            (TreeOps : Ptree.Tree_Operations with type a = NodeF.n with type b = Edge.e) =
struct

    type a = NodeF.n
    type b = Edge.e
    type tree = (a, b) Ptree.p_tree 

    type replicate = { m:bool; b:bool }

    type wtree = { t : tree; slk : float array; root : MlStaticCS.t; }

    (* RELL method; called Relative Estimated Log-Likelihood *)
    let rell = {m=false;b=false;}

    (* Full optimization; model and branches *)
    let full = {m=true; b=true; }

    (* Return the number of replicates to perform; this is based on the
       replicate optimization level; 10000 with full optimization would take an
       extremely long time; these numbers are based on papers on these subjects *)
    let get_default_reps n rep = match n with
        | Some x -> x
        | None when rep.m || rep.b -> 50
        | None -> 10_000

    (* from a list of characters, select the ONE data-set applicable; raise
       Not_found if the dataset is not contained in the list. opps! *)
    let get_matching_dataset t char roots : MlStaticCS.t =
        let chars = Data.get_chars_codes_comp (Ptree.get_data t) char in
        let (r,_) =
            List.find
                (fun (r,_) ->
                    let codes = MlStaticCS.get_codes r in
                    assert( Array.fold_left (fun a x -> a && (List.mem x chars)) true codes);
                    List.mem codes.(0) chars)
                roots
        in
        r

    (* Return if we can perform these tests on the tree; these requrements are
       that it is mlstatic data, of a complete (non-disjoint) tree (one root). *)
    let can_perform_stat_tests chars tree : bool =
        match Ptree.get_roots tree with
        | [x] ->
            begin match x.Ptree.root_median with
                | Some (`Edge _,n)
                | Some (`Single _,n) ->
                    let ml_chars = NodeF.get_mlstatic None n in
                    begin
                        try ignore (get_matching_dataset tree chars ml_chars);
                            true
                        with | Not_found -> false
                    end
                | None -> false
            end
        | _ -> false
    
    (* return the singular characters mlstatic data of a tree *)
    let get_ml_root_data chars tree : MlStaticCS.t =
        assert( can_perform_stat_tests chars tree );
        match (List.hd (Ptree.get_roots tree)).Ptree.root_median with
        | Some (`Edge _,n)
            (* preliminary = final in mlstatic data; fst or snd is fine below *)
        | Some (`Single _,n) ->
            get_matching_dataset tree chars (NodeF.get_mlstatic None n)
        | None               -> assert false

    (* create a structure with the necessary parts for test statistics; reduces
       cost of queries to this data, which would happen regularly.*)
    let create_wrapped_tree chars tree : wtree =
        let t_root = get_ml_root_data chars tree in
        let t_slk  = Array.map (snd) (MlStaticCS.site_likelihood t_root) in
        { t = tree; slk = t_slk; root=t_root; }

    (* general helper function to return cost; all these methods assume no
      negation of loglikelihood, so we transform to that standard temporarily *)
    let get_ml_cost x =
        ~-. (MlStaticCS.median_cost x.root)

    (* return a set of weights for replicated data *)
    let bootstrap_weights =
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

    (* return the cumulative distribution function; paired with a bisect to
       lookup randomly uniform data that was compressed and have weights
       associated to each column in the data-set. *)
    let get_cdf wt =
        let t_lks = MlStaticCS.site_likelihood wt.root in
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

    (* return the cost of a replicate; apply the replicate optimization level *)
    let replicate_cost rep t w =
        let cost_of_rell_bootstrap t_lks boot_weights : float =
            assert( (Array.length t_lks) = (Array.length boot_weights) );
            let t_cost = ref 0.0 in
            for i = 0 to (Array.length boot_weights)-1 do
                t_cost := !t_cost +. (t_lks.(i) *. boot_weights.(i));
            done;
            ~-. !t_cost
        and cost_of_full_bootstrap model branch tree boot_weights : float =
            assert false
        in
        match rep.m,rep.b with
        | false, false -> cost_of_rell_bootstrap t.slk w
        | _    , _     -> cost_of_full_bootstrap rep.m rep.b t.t w


    (* KH test of a priori trees *)
    let kh ?n ?(p=0.05) ?(rep=rell) ?(chars=`All) t1 t2 =
        let n = get_default_reps n rep in
        (** Calculate the cost of each tree from rell re-sampling *)
        let get_cost_of_bootstrap t1 t2 boot_weights : float =
            let t1_cost = replicate_cost rep t1 boot_weights
            and t2_cost = replicate_cost rep t2 boot_weights in
            t1_cost -. t2_cost
        in
        (** Main Components of the algorithm.. *)
        (* 0. initial variables for computation *)
        let t1  = create_wrapped_tree chars t1
        and t2  = create_wrapped_tree chars t2 in
        let cdf = get_cdf t1 in
        (* 1. determine the test statistic *)
        let d_0 = (get_ml_cost t1) -. (get_ml_cost t2) in
        (* 2. test statistic of bootstrap replicates *)
        let d_i =
            Array.init n
                (fun _ -> get_cost_of_bootstrap t1 t2 (bootstrap_weights cdf))
        in
        if debug_kh then begin
            Printf.printf "BOOTSTRAP STATISTICS: [%f]" d_0;
            Array.iter (fun x -> Printf.printf ", %f" x) d_i;
            print_newline ()
        end;
        (* 3. centered test-statistic for bootstrap replicates *)
        let dc_i =
            let f_n = float_of_int n in
            let avg = (Array.fold_left (fun a x -> a+.x) 0.0 d_i) /. f_n in
            (Array.map (fun x -> x -. avg) d_i)
        in
        if debug_kh then begin
            Printf.printf "CENTERED BOOTSTRAP STATISTICS: %f" d_0;
            Array.iter (fun x -> Printf.printf ", %f" x) dc_i;
            print_newline ()
        end;
        (* 4. Does it pass confidence test of alpha value? *)
        ()
            

    let sh ?n ?(p=0.05) ?(rep=rell) ?(chars=`All) ts =
        (** STEP 0: Setup Structures and variables necessary *)
        let replicates = get_default_reps n rep in
        let ts =
            ts --> List.map (create_wrapped_tree chars)
               --> Array.of_list
               -|> Array.sort
                    (fun x y -> Pervasives.compare (get_ml_cost y) (get_ml_cost x))
        in
        let m = Array.length ts in
        assert( m > 1 );
        let cdf = get_cdf ts.(0) in
        if debug_sh then begin
            Printf.printf "Initial Costs\n\t%!";
            Array.iter (fun x -> Printf.printf "%f, " (get_ml_cost x)) ts;
            print_newline ()
        end;
        (** STEP 1: Calculate test statistic *)
        let t =
            Array.map (fun x -> (get_ml_cost ts.(0)) -. (get_ml_cost x)) ts
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
                let weights = bootstrap_weights cdf in
                for alpha = 0 to m-1 do
                    res.(alpha).(i) <- replicate_cost rep ts.(alpha) weights
                done;
            done;
            res
        in
        if debug_sh then begin
            Printf.printf "B Matrix - Bootstrap Replicates\n%!";
            let f_n = float_of_int replicates in
            Array.iteri
                (fun alpha ib ->
                    Printf.printf "%d\t" alpha;
                    Array.iter (fun x -> Printf.printf "%f\t" x) ib;
                    let avg = (Array.fold_left (fun x h -> x+.h) 0.0 ib) /. f_n in
                    Printf.printf "\t= %f\n%!" avg)
                b
        end;
        (** STEP 3: Centering *)
        let r =
            let f_n = float_of_int replicates in
            Array.init m
                (fun alpha ->
                    let avg = (Array.fold_left (fun x h -> x +. h) 0.0 b.(alpha)) /. f_n in
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
            let res = Array.make_matrix m replicates 0.0 in
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
                                  i (get_ml_cost ts.(i)) t.(i) p.(i))
                p;
        end;
        (** STEP 6: Return set such that, P_star < P_i *)
(*        Array.mapi (fun i _ -> ts.(i),t.(i),p.(i)) p*)
        ()


    let au ?n ?(rep=rell) ?(k=5) ?(chars=`All) ts =
        (** STEP 0 : Basic setup for every one of these algorithms *)
        let ts =
            ts --> List.map (create_wrapped_tree chars)
               --> Array.of_list
               -|> Array.sort
                    (fun x y -> Pervasives.compare (get_ml_cost y) (get_ml_cost x))
        in
        let b = get_default_reps n rep in
        let m = Array.length ts in
        assert( m > 1 );
        let cdf = get_cdf ts.(0) in
        let n   = cdf.((Array.length cdf)-1) in
        if debug_au then begin
            Printf.printf "Initial Costs\n\t%!";
            Array.iter (fun x -> Printf.printf "%f, " (get_ml_cost x)) ts;
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
                            let w = bootstrap_weights ~n:n' cdf in
                            Array.map
                                (fun t -> (replicate_cost rep t w) *. scale) ts))
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
            let i = [| 1.0; 1.0 |],opt_function [| 1.0; 1.0 |] in
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
        ()

end


(** This is a little test application for the module. Uncomment and compile, the
    camlp4 tags need to be set (modify _tags file and add, 
        "mlTestStat.ml" : pp(camlp4orf), use_camlp4o, use_extensions
    usage: ./mlTestStat.native <LOAD SCRIPT> <STAT TYPE> <REPLICATES>    *)

module MLTest = Make (AllDirNode.AllDirF) (Edge.LazyEdge) (AllDirChar.F)
let test file s_type n =
    let phylo_to_mltest : (Phylo.a, Phylo.b) Ptree.p_tree -> (MLTest.a, MLTest.b) Ptree.p_tree = Obj.magic in
    Status.set_verbosity `None;
    let ()    = (POY run ([file])) in
    let ts = List.map phylo_to_mltest (Phylo.Runtime.trees ()) in
    let () = match s_type,ts with
        | "kh",x::y::_ -> MLTest.kh ~n x y
        | "sh", ts     -> MLTest.sh ~n ts
        | "au", ts     -> MLTest.au ~n ts
    in
    ()
let () =
    test Sys.argv.(1) Sys.argv.(2) (int_of_string Sys.argv.(3))
