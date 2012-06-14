
let () = Random.self_init ()

module Test = struct

    (** The Hummelblau function is a multi-modal function with four optimal points.
        often used in testing numerical functions,
        f( 3.000, 2.000) = 0.0     f(-2.805, 3.131) = 0.0
        f(-3.779,-3.283) = 0.0     f( 3.584,-1.848) = 0.0 *)
    let hummelblau_function p = 
        assert( (Array.length p) = 2 );
        ((p.(0)**2.0) +. p.(1) -. 11.0)**2.0 +.  ((p.(1)**2.0) +. p.(0) -. 7.0)**2.0

    (** The Booth function has several local minima, but a global minima at 1,3.  *)
    let booth_function p = 
        assert( (Array.length p) = 2 );
        ((p.(0) +. p.(1) +. p.(1) -. 7.0)**2.0) +.  ((p.(0) +. p.(0) +. p.(1) -. 5.0)**2.0)

    (** Dixon & Price Function; global minima at f( x* ) = 0. *)
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


    (** The rosenbrock function is a non-convex function used to test numerical
        functions. The optimial point is at (1,1) *)
    let rosenbrock_function p =
        assert( (Array.length p) = 2 );
        (1.0 -. p.(0))**2.0 +. (100.0 *. (p.(1) -. p.(0)**2.0)**2.0)

    (** Mutli-dimensional generalisation of the rosenbrock function; may have
        multiple minimums based on the number of dimensions (N=3 has one, N=4 has 2) *)
    let rosenbrock_function_multi p =
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
        impossible to find an optimal point with a gradient method. *)
    let rosenbrock_function_stoc ?(dist=(fun () -> Random.float 1.0)) p =
        let total = ref 0.0 in
        for i = 0 to (Array.length p) - 2 do
            let two_d = (1.0 -. p.(i))**2.0 +.
                        (100.0 *. (dist ()) *. (p.(i+1) -. p.(i)**2.0)**2.0)
            in
            total := !total +. two_d;
        done;
        !total

    (** Quadratic function; simple test to ensure algorithms are working
        properly. The minimum is x= (0,0...0) *)
    let quadratic_function p =
        Array.fold_left (fun acc x -> acc +. (x *. x)) 0.0 p

    (** Quadratic function with stochasitc coefficients; same solution as above,
        the randomness ensures undiffernentiaable, is x= (0,0...0) *)
    let quadratic_function_stoc ?(dist=(fun () -> Random.float 1.0)) p =
        Array.fold_left (fun acc x -> acc +. (dist ()) *. (x *. x)) 0.0 p

    (** The rastrigein function is a non-convext function for performance
        testings. global minimum at  {0,0...}, and range between +/- 5.12. *)
    let rastigin_function p =
        let total = ref 0.0 and a = 10.0 and pi = acos (-1.0) in
        for i = 0 to (Array.length p)-1 do
            let p_i = min (max ~-.5.12 p.(i)) 5.12 in
            total := !total +. (p_i *. p_i) -. a *. (cos (2.0 *. pi *. p_i))
        done;
        !total +. (a *. (float_of_int (Array.length p)))
end

(** An object that keeps track of running variance, mean, and standard deviation
    of a simulation. **)
class running_stats = object(self)

    val mutable   k = 0
    val mutable m_k = 0.0
    val mutable s_k = 0.0

    val mutable max_x = 0.0
    val mutable min_x = max_float

    method push x_k (*t_k*) =
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

let test_battery channel f n r () =
    Status.set_verbosity `None; (* suppress warning messages *)
    let f = (fun p -> let fp = f p in ((), fp)) in
    let sim_obj = new running_stats and bfgs_obj  = new running_stats
    and sub_obj = new running_stats and brent_obj = new running_stats in
    Printf.fprintf channel "Simplex         Subplex         Brents          BFGS    \n%!";
    Printf.fprintf channel "--------------------------------------------------------\n%!";
    for i = 0 to r do
        (* Randomly selection an initial vector *)
        let init = Array.init n (fun _ -> Random.float 1.0) in
        (* Run algorithms, send each cost to the running_statistics module *)
        let _,(_,  sim_cost) = Numerical.simplex_method init f init (f init) in
        let _,(_,  sub_cost) = Numerical.subplex_method f init (f init) in
        let _,(_,brent_cost) = Numerical.brents_method_multi (init,(f init)) f in
        let _,(_, bfgs_cost) = Numerical.bfgs_method f init (f init) in
        Printf.fprintf channel "%f\t%f\t%f\t%f\n%!"
                       sim_cost sub_cost brent_cost bfgs_cost;
        sim_obj#push sim_cost; bfgs_obj#push bfgs_cost;
        sub_obj#push sub_cost; brent_obj#push brent_cost;
    done;
    Printf.fprintf channel "--------------------------------------------------------\n%!";
    Printf.fprintf channel "Simplex         Subplex         Brents          BFGS    \n%!";
    Printf.fprintf channel "--------------------------------------------------------\n%!";
    Printf.fprintf channel "%f\t%f\t%f\t%f\n%!"
        (sim_obj#mean ()) (sub_obj#mean ()) (brent_obj#mean ()) (bfgs_obj#mean ());
    Printf.fprintf channel "%f\t%f\t%f\t%f\n%!"
        (sim_obj#std_dev ()) (sub_obj#std_dev ()) (brent_obj#std_dev ()) (bfgs_obj#std_dev ());
    Printf.fprintf channel "%f\t%f\t%f\t%f\n%!"
        (sim_obj#min ()) (sub_obj#min ()) (brent_obj#min ()) (bfgs_obj#min ());
    Printf.fprintf channel "%f\t%f\t%f\t%f\n%!"
        (sim_obj#max ()) (sub_obj#max ()) (brent_obj#max ()) (bfgs_obj#max ());
    Printf.fprintf channel "--------------------------------------------------------\n%!";
    ()

let main f_num d_num i_num channel =
    let dist () = Random.float 1.0 in
    let f = match f_num with
        | 0 when d_num = 2 -> Test.hummelblau_function
        | 1 when d_num = 2 -> Test.booth_function
        | 2                -> Test.ackley_function
        | 3 when d_num = 2 -> Test.rosenbrock_function
        | 3                -> Test.rosenbrock_function_multi
        | 4                -> Test.rosenbrock_function_stoc ~dist
        | 5                -> Test.quadratic_function
        | 6                -> Test.quadratic_function_stoc ~dist
        | 7                -> Test.rastigin_function
        | _                -> failwith "Usage: ./test_numerical [0-7] D I [F]"
    in
    test_battery channel f d_num i_num ()

let () =
    let channel = 
        try match Sys.argv.(4) with
            | "stdout" -> stdout
            | "stderr" -> stderr
            | filename -> open_out filename
        with | _       -> stdout
    in
    main (int_of_string Sys.argv.(1)) (int_of_string Sys.argv.(2)) (int_of_string Sys.argv.(3)) channel
