(* POY 5.0 Beta. A phylogenetic analysis program using Dynamic Homologies.    *)
(* Copyright (C) 2013 Andrés Varón, Lin Hong, Nicholas Lucaroni, Ward Wheeler,*)
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

let () = SadmanOutput.register "MlTestStat" "$Revision: 3395 $"

let (-->) a b = b a
let (-|>) a b = let () = b a in a

let debug_cdf = false
and debug_boot= false

and debug_kh  = false
and debug_sh  = false

let failwithf format =
    Printf.ksprintf failwith format

let info_user_message format =
    Printf.ksprintf (Status.user_message Status.Information) format

let warn_user_message format =
    Printf.ksprintf (Status.user_message Status.Warning) format

let error_user_message format =
    Printf.ksprintf (Status.user_message Status.Error) format

let outputt table =
    Status.output_table Status.Information table


(** {6 Outside Types --irrespective of the functorization. This is allow modules
       to process arguments as soon as they come in *)

(** Define a type to determine the way to optimize the tree to obtain
    likelihood scores for replicates of the tree. *)
type replicate = { m:bool; b:bool }

(* RELL method; called Relative Estimated Log-Likelihood *)
let rell = {m=false;b=false;}

(* Full optimization; model and branches *)
let full = {m=true; b=true; }

(* Full optimization; model and branches *)
let part = {m=false; b=true; }

(* Return the number of replicates to perform; this is based on the
   replicate optimization level; 10000 with full optimization would take an
   extremely long time; these numbers are based on papers on these subjects *)
let get_default_reps n rep = match n with
    | Some x -> x
    | None when rep.m || rep.b -> 50
    | None -> 10_000

let process_methods_arguments args =
    let folder (t,c,n,r) = function
        | `SH | `KH as t -> (Some t,c,n,r)
        | `Characters c -> (t,c,n,r)
        | `Replicates n -> (t,c,Some n,r)
        | `ReplicateOpt (mt,bt) -> (t,c,n,{m=mt;b=bt})
    in
    let (t,c,n,r) = List.fold_left folder (None,`All,None,rell) args in
    match t with
    | Some t -> t,c,n,r
    | None   ->
        Status.user_message Status.Error
            "No@ Topology@ Selection@ Method@ specified.@ Please@ include@ sh@ or@ kh.";
        raise Not_found

module type S = sig

    type a 
    type b

    exception Incorrect_Data

    type tree = (a,b) Ptree.p_tree
    type wtree =
        { t : tree;
          slk : (int * float * float) array;
          root : MlStaticCS.t;
          chars : Data.bool_characters; }

    val create_wrapped_tree : Data.bool_characters -> tree -> wtree
    val can_perform_stat_tests : Data.bool_characters -> tree -> bool

    val bootstrap_weights : ?n:int -> (int * float) array -> (int * float) array
    val replicate_cost : replicate -> wtree -> (int * float) array -> float
    val get_cdf : wtree -> (int * float) array

    val kh : ?n:int -> ?p:float -> ?rep:replicate -> ?chars:Data.bool_characters -> tree -> tree -> unit
    val sh : ?n:int -> ?p:float -> ?rep:replicate -> ?chars:Data.bool_characters -> tree list -> unit

end

module Make (NodeF : NodeSig.S with type other_n = Node.Standard.n)
            (Edge : Edge.EdgeSig with type n = NodeF.n) 
            (TreeOps : Ptree.Tree_Operations with type a = NodeF.n with type b = Edge.e) =
struct

    type a = NodeF.n
    type b = Edge.e
    type tree = (a, b) Ptree.p_tree 

    exception Incorrect_Data

    type wtree =
        { t : tree; slk : (int * float * float) array; root : MlStaticCS.t; chars : Data.bool_characters; }

    IFDEF USE_LIKELIHOOD THEN
        let get_ml_codes = MlStaticCS.get_codes
        let get_ml_cost x = ~-. (MlStaticCS.median_cost x.root)
        let ml_site_likelihood = MlStaticCS.site_likelihood
    ELSE
        let get_ml_codes _ = assert false
        let get_ml_cost _  = assert false
        let ml_site_likelihood _ = assert false
    END

    (* from a list of characters, select the ONE data-set applicable; raise
       Not_found if the dataset is not contained in the list. opps! *)
    let get_matching_dataset t char roots : MlStaticCS.t =
        let chars = Data.get_chars_codes_comp (Ptree.get_data t) char in
        let (r,_) =
            List.find
                (fun (r,_) ->
                    let codes = get_ml_codes r in
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
        let t_slk  = ml_site_likelihood t_root in
        { t = tree; slk = t_slk; root=t_root; chars = chars; }

    (* return a set of weights for replicated data *)
    let bootstrap_weights : ?n:int -> (int * float) array -> (int * float) array =
        let rec idx_of_cdf cdf (g:float) (il:int) (ih:int) : int =
            let ig = (il + ih) / 2 in
            if debug_cdf then
                Printf.printf "IDX OF CDF: %d(%f) < (%d/%f?)%f < %d(%f)\n%!"
                              il (snd cdf.(il)) ig (snd cdf.(ig)) g ih (snd cdf.(ih));
            assert( (snd cdf.(ih)) > g );
            assert( (snd cdf.(il)) < g );
            if ((snd cdf.(ig)) < g) && (g < (snd cdf.(ig+1)))
                then ig
            else if (snd cdf.(ig)) > g
                then idx_of_cdf cdf g il ig
                else idx_of_cdf cdf g ig ih
        and bootstrap_data ?n cdf =
            let n = match n with
                | None   -> int_of_float (snd cdf.((Array.length cdf)-1))
                | Some x -> x
            in
            let boot_data : (int * float) array=
                Array.init (Array.length cdf) (fun i -> fst cdf.(i),0.0) in
            let incr_tuple (a,b) = (a, b +. 1.0) in

            let upper = (snd (cdf.((Array.length cdf)-1))-.1.0) in
            for i = 0 to n do
                let r = (Random.float upper)+.1.0 in
                let i = idx_of_cdf cdf r 0 ((Array.length cdf)-1) in
                if debug_cdf then Printf.printf "IDX OF CDF: %d = %f\n%!" i r;
                boot_data.(i) <- incr_tuple boot_data.(i);
            done;
            boot_data
        in
        bootstrap_data

    (* return the cumulative distribution function; paired with a bisect to
       lookup randomly uniform data that was compressed and have weights
       associated to each column in the data-set. *)
    let get_cdf wt =
        let t_lks = wt.slk in
        let cdf = Array.make (Array.length t_lks) (0,0.0) in
        let sum = ref 0.0 in
        for i = 0 to (Array.length t_lks)-1 do
            let code,weight,_ = t_lks.(i) in
            sum := !sum +. weight;
            cdf.(i) <- (code,!sum);
        done;
        if debug_cdf then begin
            Printf.printf "CDF: ";
            Array.iter (fun (c,x) -> Printf.printf "%d|%f, " c x) cdf;
        end;
        cdf

    (* apply weights to a tree data-set with restricted costs; *)
    let reweight_data (d:Data.d) (chars:Data.bool_characters) weights =
        Array.fold_left
            (fun data (code,weight) -> Data.set_weight code weight data)
            d
            weights

    (* return the cost of a replicate; apply the replicate optimization level *)
    let replicate_cost rep t w =
        let cost_of_rell_bootstrap t_lks boot_weights : float =
            assert( (Array.length t_lks) = (Array.length boot_weights) );
            let t_cost = ref 0.0 in
            for i = 0 to (Array.length boot_weights)-1 do
                let c1,_,x = t_lks.(i) and c2,w = boot_weights.(i) in
                assert( c1 = c2 );
                t_cost := !t_cost +. (x *. w);
            done;
            ~-. !t_cost
        and cost_of_full_bootstrap model branch tree boot_weights : float =
            let data,nodes =
                boot_weights
                    --> reweight_data (Ptree.get_data tree.t) tree.chars
                    --> NodeF.load_data
            in
            let node_data : a All_sets.IntegerMap.t =
                List.fold_left
                    (fun acc x -> All_sets.IntegerMap.add (NodeF.taxon_code x) x acc)
                    All_sets.IntegerMap.empty
                    nodes
            in
            { tree.t with
                Ptree.node_data = node_data; Ptree.data = data; }
                --> TreeOps.downpass
                --> TreeOps.uppass
                --> (fun t -> ~-. (TreeOps.total_cost t `Adjusted None))
        in
        match rep.m,rep.b with
        | false, false -> cost_of_rell_bootstrap t.slk w
        | _    , _     -> cost_of_full_bootstrap rep.m rep.b t w


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
        let dc_i = let f_n = float_of_int n in
            let avg = (Array.fold_left (fun a x -> a+.x) 0.0 d_i) /. f_n in
            (Array.map (fun x -> x -. avg) d_i)
        in
        if debug_kh then begin
            Printf.printf "CENTERED BOOTSTRAP STATISTICS: %f" d_0;
            Array.iter (fun x -> Printf.printf ", %f" x) dc_i;
            print_newline ()
        end;
        (** 4: Output information *)
        let p_value =
            let p = 
                Array.fold_left
                    (fun acc x -> if x > d_0 then succ acc else acc) 0 dc_i
            in
            (float_of_int p) /. (float_of_int n)
        in
        info_user_message "@[<4>@[KH@ test@ p-value@ from@ %d@ replicates@]" n;
        let matrix =
            [| [| "MLE Tree:"; Printf.sprintf "%f" (get_ml_cost t1); |];
               [| "Tree:"; Printf.sprintf "%f" (get_ml_cost t2); |];
               [| "p-value:"; Printf.sprintf "%f" p_value; |] |]
        in
        outputt matrix;
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
        let alpha_mle = 0 in
        let m = Array.length ts in
        assert( m > 1 );
        let cdf = get_cdf ts.(alpha_mle) in
        if debug_sh then begin
            Printf.printf "Initial Costs\n\t%!";
            Array.iter (fun x -> Printf.printf "%f, " (get_ml_cost x)) ts;
            print_newline ()
        end;
        (** STEP 1: Calculate test statistic *)
        let t =
            Array.map (fun x -> (get_ml_cost ts.(alpha_mle)) -. (get_ml_cost x)) ts
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
                    (* ???: for KH we may replace this with alpha_mle *)
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
        (** STEP 6: Output information *)
        let matrix =
            Array.init ((Array.length p)+1)
                (fun i ->
                    if i = 0 then
                        [| "Tree"; "loglk"; "p-value" |]
                    else
                        let i = i-1 in
                        [| Printf.sprintf "%d" i;
                           Printf.sprintf "%f" (get_ml_cost ts.(i));
                           Printf.sprintf "%f" p.(i); |])
        in
        info_user_message
            "@[<4>@[SH@ test@ P-Values@ from@ %d@ replicates@]" replicates;
        outputt matrix;
        info_user_message "@]";
        ()

end
