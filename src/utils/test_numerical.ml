(* POY 5.1.1. A phylogenetic analysis program using Dynamic Homologies.       *)
(* Copyright (C) 2014 Andrés Varón, Lin Hong, Nicholas Lucaroni, Ward Wheeler,*)
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

let () = Random.self_init ()

module IntMap = All_sets.IntegerMap

let (-->) a b = b a

module Test = struct

    (** Quadratic function; simple test to ensure algorithms are working
        properly. The minimum is x= (0,0...0) *)
    let quadratic_function p =
        Array.fold_left (fun acc x -> acc +. (x *. x)) 0.0 p

    (** Quadratic function with stochasitc coefficients; same solution as above,
        the randomness ensures undiffernentiable, is x= (0,0...0). The default
        is a uniform distribution, others can be passed. *)
    let quadratic_function_stoc ?(dist=(fun () -> Random.float 1.0)) p =
        Array.fold_left (fun acc x -> acc +. (dist ()) *. (x *. x)) 0.0 p

    (** Mutli-dimensional generalisation of the rosenbrock function; may have
        multiple minimums based on the number of dimensions (N=3 has one, N=4 has 2) *)
    let rosenbrock_function p =
        let total = ref 0.0 in
        for i = 0 to (Array.length p) - 2 do
            let two_d = (1.0 -. p.(i))**2.0 +.
                        (100.0 *. (p.(i+1) -. p.(i)**2.0)**2.0)
            in
            total := !total +. two_d;
        done;
        !total

    (** Stochastic version of the rosenbrock function, passes into a distribution of
        a single parameter. The default is a uniform distribution. This makes it
        impossible to find an optimal point with a gradient method. The default
        is a uniform distribution, others can be passed. *)
    let rosenbrock_function_stoc ?(dist=(fun () -> Random.float 1.0)) p =
        let total = ref 0.0 in
        for i = 0 to (Array.length p) - 2 do
            let two_d = (1.0 -. p.(i))**2.0 +.
                        (100.0 *. (dist ()) *. (p.(i+1) -. p.(i)**2.0)**2.0)
            in
            total := !total +. two_d;
        done;
        !total

    (** The rastrigein function is a non-convext function for performance
        testings. global minimum at  {0,0...}, and range between +/- 5.12. *)
    let rastigin_function p =
        let total = ref 0.0 and a = 10.0 and pi = acos (-1.0) in
        for i = 0 to (Array.length p)-1 do
            let p_i = min (max ~-.5.12 p.(i)) 5.12 in
            total := !total +. (p_i *. p_i) -. a *. (cos (2.0 *. pi *. p_i))
        done;
        !total +. (a *. (float_of_int (Array.length p)))

    (** Schwefel's function had a global minium that is geometrically distant
        from the other nearest local minimum. This tests convergance of the wrong
        direction of the optimization routine *)
    let schwefel_function p =
        Array.fold_left
            (fun acc x -> acc +. (~-. x *. (sin (sqrt (abs_float x))))) 0.0 p

    (** Sum of different powers is a uni-modal function *)
    let sum_powers p =
        let summ = ref 0.0 in
        for i = 0 to (Array.length p)-1 do
            summ := !summ +. (abs_float p.(i))**(float_of_int (i+1));
        done;
        !summ

    (** Ackley function; global minima at f( x* ) = 0. *)
    let ackley_function p =
        let n = Array.length p and pi = acos (-1.0) in
        let left_sum = ref 0.0 and right_sum = ref 0.0 in
        for i = 0 to n-1 do
            left_sum := !left_sum +. (p.(i) ** 2.0);
            right_sum := !right_sum +. (cos (2.0 *. pi *. p.(i)));
        done;
        let left = ~-. 0.2 *. (sqrt (!left_sum /. (float_of_int n)))
        and right = !right_sum /. (float_of_int n) in
        ~-. 20.0 *. (exp left) -. (exp right) +. 20.0 +. (exp 1.0)
    
    (** Michalewiczs function is a multi-modal test function with n! local
        optima, m, defines the steepness (set to 10.0), points outside an area
        give no information, while a few steep cliffs direct to minima. *)
    let michalewicz_function p =
        let sum = ref 0.0 and pi = acos (-1.0) in
        for i = 0 to (Array.length p) - 1 do
            let part = (sin ((float_of_int i)*.p.(i)*.p.(i) /. pi))**(20.0) in
            sum := !sum +. (sin p.(i)) *. part;
        done;
        !sum

    (** The Hummelblau function is a multi-modal function with four optimal points.
        often used in testing numerical functions,
        f( 3.000, 2.000) = 0.0     f(-2.805, 3.131) = 0.0
        f(-3.779,-3.283) = 0.0     f( 3.584,-1.848) = 0.0 *)
    let hummelblau_function p = 
        assert( (Array.length p) = 2 );
        ((p.(0)**2.0) +. p.(1) -. 11.0)**2.0 +.  ((p.(1)**2.0) +. p.(0) -. 7.0)**2.0

    (** The Booth function has several local minima, but a global minima at 1,3. *)
    let booth_function p = 
        assert( (Array.length p) = 2 );
        let left = ((p.(0) +. p.(1) +. p.(1) -. 7.0)**2.0)
        and rght = ((p.(0) +. p.(0) +. p.(1) -. 5.0)**2.0) in
        left +. rght

    (** Diagnose a tree based on an array of rate parameters that get normalized
        before the tree is diagnosed. files is a list of files that will load a
        tree into memory and do an initial transform to a likelihood model that
        this function should test against. *)
    let poy_diagnose_likelihood files =
        let print_model model = 
            MlModel.output_model (output_string stdout) None `Nexus model None in
        let update_model_to_vector tree chars =
            let model = Data.get_likelihood_model tree.Ptree.data chars in
            let func = match MlModel.get_update_function_for_model model with
                | Some f -> f | None   -> assert false
            in
            (fun vec ->
                let model = func model vec in
                let data,nodes =
                    model
                        --> Data.apply_likelihood_model_on_chars tree.Ptree.data chars
                        --> Data.categorize
                        --> AllDirNode.AllDirF.load_data
                in
                let nodes =
                    List.fold_left
                        (fun acc x -> IntMap.add (AllDirNode.AllDirF.taxon_code x) x acc)
                        IntMap.empty
                        nodes
                in
                let tree = 
                    { tree with Ptree.data = data; Ptree.node_data = nodes; }
                in
                model,(tree --> Phylo.PhyloTree.downpass
                            --> Phylo.PhyloTree.uppass
                            --> Phylo.PhyloTree.get_cost))
        and get_model_dimension tree chars =
            chars --> Data.get_likelihood_model tree.Ptree.data
                  --> MlModel.get_current_parameters_for_model
                  --> Array.length
        in
        Status.set_verbosity `None;
        let () = (POY run ([files])) in
        let t = List.hd (Phylo.Runtime.trees ()) in
        let chars = 
            match Data.categorize_likelihood_chars_by_model t.Ptree.data `All with
            | [x] -> x
            |  _  -> failwith "This application only supports one character set"
        in
        let f = update_model_to_vector t chars in
        f, print_model, get_model_dimension t chars
end

let test_battery channel f n p r () =
    Status.set_verbosity `None; (* suppress warning messages *)
    (* for calculating the mean, standard deviation of algorithm results *)
    let sim_obj = new Numerical.running_stats and bfg_obj = new Numerical.running_stats
    and sub_obj = new Numerical.running_stats and bnt_obj = new Numerical.running_stats in
    (* for calculating avg, and worst/best timings for algorithm *)
    let sim_tme = new Numerical.running_stats and bfg_tme = new Numerical.running_stats
    and sub_tme = new Numerical.running_stats and bnt_tme = new Numerical.running_stats in
    let f t = let t = ref t in
              (fun p -> incr t; let _,fp = f p in float_of_int !t,fp)
    in
    Printf.fprintf channel "I       Simplex         Subplex         Brents          BFGS       \n%!";
    Printf.fprintf channel "-------------------------------------------------------------------\n%!";
    for i = 0 to r do
        (* Randomly selection an initial vector *)
        let   init = Array.init n (fun _ -> Random.float 1.0) in
        let f_init = f 0 init in
        (* Run algorithms, send each cost to the running_statistics module *)
        let _,(sim_time, sim_cost) = Numerical.simplex_method      (f 1) (init,(f_init)) in
        let _,(sub_time, sub_cost) = Numerical.subplex_method      (f 1) (init,(f_init)) in
        let _,(bnt_time, bnt_cost) = Numerical.brents_method_multi (f 1) (init,(f_init)) in
        let _,(bfg_time, bfg_cost) = Numerical.bfgs_method         (f 1) (init,(f_init)) in
        Printf.fprintf channel "%d\t%f\t%f\t%f\t%f\n%!" i sim_cost sub_cost bnt_cost bfg_cost;
        sim_obj#push sim_cost; bfg_obj#push bfg_cost; sub_obj#push sub_cost; bnt_obj#push bnt_cost;
        sim_tme#push sim_time; bfg_tme#push bfg_time; sub_tme#push sub_time; bnt_tme#push bnt_time;
    done;
    Printf.fprintf channel "-------------------------------------------------------------------\n%!";
    Printf.fprintf channel "Cost    Simplex         Subplex         Brents          BFGS       \n%!";
    Printf.fprintf channel "-------------------------------------------------------------------\n%!";
    Printf.fprintf channel "Mean\t%f\t%f\t%f\t%f\t\n%!"
        (sim_obj#mean ()) (sub_obj#mean ()) (bnt_obj#mean ()) (bfg_obj#mean ());
    Printf.fprintf channel "StdDev\t%f\t%f\t%f\t%f\n%!"
        (sim_obj#std_dev ()) (sub_obj#std_dev ()) (bnt_obj#std_dev ()) (bfg_obj#std_dev ());
    Printf.fprintf channel "Min\t%f\t%f\t%f\t%f\n%!"
        (sim_obj#min ()) (sub_obj#min ()) (bnt_obj#min ()) (bfg_obj#min ());
    Printf.fprintf channel "Max\t%f\t%f\t%f\t%f\n%!"
        (sim_obj#max ()) (sub_obj#max ()) (bnt_obj#max ()) (bfg_obj#max ());
    Printf.fprintf channel "-------------------------------------------------------------------\n%!";
    Printf.fprintf channel "Time    Simplex         Subplex         Brents          BFGS       \n%!";
    Printf.fprintf channel "-------------------------------------------------------------------\n%!";
    Printf.fprintf channel "Mean\t%f\t%f\t%f\t%f\t\n%!"
        (sim_tme#mean ()) (sub_tme#mean ()) (bnt_tme#mean ()) (bfg_tme#mean ());
    Printf.fprintf channel "StdDev\t%f\t%f\t%f\t%f\n%!"
        (sim_tme#std_dev ()) (sub_tme#std_dev ()) (bnt_tme#std_dev ()) (bfg_tme#std_dev ());
    Printf.fprintf channel "Min\t%f\t%f\t%f\t%f\n%!"
        (sim_tme#min ()) (sub_tme#min ()) (bnt_tme#min ()) (bfg_tme#min ());
    Printf.fprintf channel "Max\t%f\t%f\t%f\t%f\n%!"
        (sim_tme#max ()) (sub_tme#max ()) (bnt_tme#max ()) (bfg_tme#max ());
    Printf.fprintf channel "-------------------------------------------------------------------\n%!";
    ()

let main_num f d i channel =
    let dist () = Random.float 1.0 in
    let unit_wrapper f = (fun x -> ((), f x)) and unit_printer () = () in
    let f = match f with
        | 0 when d = 2 -> unit_wrapper Test.hummelblau_function
        | 1 when d = 2 -> unit_wrapper Test.booth_function
        | 2            -> unit_wrapper Test.ackley_function
        | 3            -> unit_wrapper Test.rosenbrock_function
        | 4            -> unit_wrapper (Test.rosenbrock_function_stoc ~dist)
        | 5            -> unit_wrapper Test.quadratic_function
        | 6            -> unit_wrapper (Test.quadratic_function_stoc ~dist)
        | 7            -> unit_wrapper Test.rastigin_function
        | 8            -> unit_wrapper Test.schwefel_function
        | 9            -> unit_wrapper Test.sum_powers
        | 10           -> unit_wrapper Test.michalewicz_function
        | _            -> assert false
    in
    test_battery channel f d unit_printer i ()

let main () =
    try let f_num = int_of_string Sys.argv.(1) in
        let i = int_of_string Sys.argv.(3) in
        let channel =
            try match Sys.argv.(4) with
                | "stdout" -> stdout
                | "stderr" -> stderr
                | filename -> open_out filename
            with | _       -> stdout
        in
        main_num f_num (int_of_string Sys.argv.(2)) i channel
    with | Failure "int_of_string" ->
        let f,p,d = Test.poy_diagnose_likelihood Sys.argv.(1) in
        let channel =
            try match Sys.argv.(3) with
                | "stdout" -> stdout
                | "stderr" -> stderr
                | filename -> open_out filename
            with | _       -> stdout
        in
        test_battery channel f d p (int_of_string Sys.argv.(2)) ()

let usage = (Sys.argv.(0))^" <FILE|FUNC DIMS> <ITER> [OUTF]"

let desc = "  FILE\t- A POY script that loads a tree and transforms to likelihood\n"
         ^ "  FUNC\t- A function number 0-10\n"
         ^ "  DIMS\t- Number of dimensions when using functions\n"
         ^ "  ITER\t- Size of sample; number of times to run optimization\n"
         ^ "  OUTF\t- (optional, default to screen). Output location a file, or stderr or stdout\n"

let () =
    try
        if ("--help" = Sys.argv.(1)) || ("-h" = Sys.argv.(1))
            then Printf.printf "%s\n%s\n%!" usage desc
            else main ()
    with _ ->
        Printf.printf "%s\n%s\n%!" usage desc

