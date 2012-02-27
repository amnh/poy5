(* Program to generate MAP tcm matrices *)
(*Need to do:
       Add ratio param for indel prior and transitions
       Clean up main code--functionalize matrix adds
       output files based on stem
       read configuration data from file
*)
let (-->) a b = b a
let ba_of_array1 x = Bigarray.Array1.of_array Bigarray.float64 Bigarray.c_layout x
and ba_of_array2 x = Bigarray.Array2.of_array Bigarray.float64 Bigarray.c_layout x

let choose_2 in_val =
let out_val = (((in_val * in_val) - in_val) / 2) in
out_val
;;

let set_up_params () =
  Printf.printf "Setting values \n";
  let alpha_bet_size = 5 in
  let 
     priors = Array.create alpha_bet_size 1.0 and
     transitions = Array.create (choose_2 alpha_bet_size) 1.0 and
     invariant_min = 0.0 and
     invariant_max = 1.0 and
     alpha_min = 0.0 and
     alpha_max = 50.0 and
     edge_length = 10.0 and
     num_rate_classes_min = 1 and
     num_rate_classes_max = 7 and
     iterations = 100000000 and
     epsilon = 0.000001 and
     prior_dist = "Dirichlet" and  (*Fixed or Dirichlet*)
     transition_dist = "Dirichlet" and   (*Fixed or Dirichlet*)
     invariant_dist = "Uniform" and   (*Fixed or Uniform*)
     alpha_dist = "Uniform" and   (*Fixed or Uniform*)
     num_rate_classes_dist = "Uniform" and   (*Fixed or Uniform*)
     edge_time_dist = "Exponential" and   (*Exponential or Uniform*)
     file_stub_name = "MAPA_"
  in
  (alpha_bet_size, priors, transitions, invariant_min, invariant_max, alpha_min, 
       alpha_max, num_rate_classes_min, num_rate_classes_max, iterations, epsilon, edge_length,
       prior_dist, transition_dist, invariant_dist, alpha_dist, num_rate_classes_dist, 
       edge_time_dist, file_stub_name)
;;

let print_vals alphabet_size prior_array transition_array invariant_min invariant_max 
  alpha_min alpha_max num_rate_classes_min num_rate_classes_max max_iterations epsilon edge_length 
  prior_dist transition_dist invariant_dist alpha_dist num_rate_classes_dist edge_time_dist 
  file_stub_name =
Printf.printf "Alphabet size = \t\t\t%d\n" alphabet_size;
Printf.printf "Prior parameters = %s\t\t" prior_dist;
  for  i = 0 to alphabet_size - 1 do
    Printf.printf "%f " prior_array.(i);
  done;
Printf.printf "\n";
Printf.printf "Transition parameters = %s\t" transition_dist;
  for  i = 0 to (choose_2 alphabet_size) - 1 do
    Printf.printf "%f " transition_array.(i);
  done;
Printf.printf "\n";
Printf.printf "Invariant sites = %s\t\t[%f, %f]\n" invariant_dist invariant_min invariant_max;
Printf.printf "Alpha rate param = %s\t\t[%f, %f]\n" alpha_dist alpha_min alpha_max;
Printf.printf "Number rate classes = %s\t\t[%d, %d]\n" num_rate_classes_dist num_rate_classes_min num_rate_classes_max;
Printf.printf "Edge time = %s\t\t\t%f\n" edge_time_dist edge_length;
Printf.printf "Maximum iterations = \t\t\t%d\n" max_iterations;
Printf.printf "Epsilon = \t\t\t\t%f\n" epsilon;
Printf.printf "File stub = %s\n" file_stub_name;
;;

let _ = Random.self_init ();;

let exponential_variant alpha =
  let rand_float = Random.float 1.0 in
  let exp_rand_float = (-1.0 *. (log rand_float) /. alpha) in
  exp_rand_float
;;
let normal_variant mu sigma =
  let u1 = Random.float 1.0 and
  u2 = Random.float 1.0 in
  let r = sqrt( -2.0 *. (log u1)) and
  theta = (2.0 *. 3.14159265358979323846264338327 *. u2) in
  let return = mu +. (sigma *. r *. (sin theta)) in
  return
;;

let rec gamma_variant alpha beta = 
  let return_val = ref 0.0 in
  if (not (alpha < 1.0)) then 
    begin
       let d = alpha -. (1.0 /. 3.0) in
       let c = 1.0 /. (sqrt (9.0 *. d)) in
       let condition = ref true in
        while !condition do
           let x = ref (normal_variant 0.0 1.0) in
           let v = ref (1.0 +. (c *. !x)) in
           while (!v <= 0.0) do
             x := normal_variant 0.0 1.0;
             v := 1.0 +. (c *. !x);
           done; 
           v := !v *. !v *. !v;
           let u = Random.float 1.0 and
           xsq = !x *. !x in
           if ((u < (1.0 -. 0.0331 *. xsq *. xsq)) or ((log u) < (0.5 *. xsq +. (d *. (1.0 -. !v +. (log !v)))))) then 
             begin
                return_val := beta *. d *. !v;
                condition:= false;
             end;
        done;
    end
  else 
    begin
      let g = gamma_variant (alpha +. 1.0) 1.0 in
      let w = Random.float 1.0 in
      return_val := beta *. g *. (w ** (1.0 /. alpha)); 
    end; 
!return_val
;;

(*Generates Dirichlet form a series of Gammas (alpha_i, 1.0) *)
let dirichlet_variant param_array =
  let n = Array.length param_array in
  let y_array = Array.create n 0.0 and
  x_array =  Array.create n 0.0 and
  y_sum = ref 0.0 in
  for i = 0 to n - 1 do
    y_array.(i) <- gamma_variant param_array.(i) 1.0;
    y_sum := !y_sum +. y_array.(i);
  done;
  for i = 0 to n - 1 do
    x_array.(i) <- y_array.(i) /. !y_sum;
  done;
x_array
;;

(* Add float elements of b to a
   assumes both have same dimensions*)
let matrix_increment a b =
  for i = 0 to (Array.length a) - 1 do
     for j = 0 to (Array.length a.(i)) - 1 do
        a.(i).(j) <- a.(i).(j) +.  b.(i).(j);
     done;
  done;
;;

(*Divide elements in matrix by float*)
let matrix_divide a b =
  for i = 0 to (Array.length a) - 1 do
     for j = 0 to (Array.length a.(i)) - 1 do
        a.(i).(j) <- a.(i).(j) /.  b;
     done;
  done;
;;

(* multiply bigarray b by c then add to a--all float*)
let matrix_increment_multiply a b c =
  for i = 0 to (Array.length a) - 1 do
     for j = 0 to (Array.length a.(i)) - 1 do
        a.(i).(j) <- a.(i).(j) +.  (b.{i,j} *. c);
     done;
  done;
;;

let () =
Printf.printf "Beginning run...\n";

(*Get parameter values*)

let (alphabet_size, prior_array, transition_array, invariant_min, invariant_max, 
      alpha_min, alpha_max, num_rate_classes_min, num_rate_classes_max, max_iterations, epsilon, edge_length,
      prior_dist, transition_dist,invariant_dist, alpha_dist, num_rate_classes_dist, 
      edge_time_dist, file_stub_name) = set_up_params ()
      in
print_vals alphabet_size prior_array transition_array invariant_min invariant_max alpha_min 
  alpha_max num_rate_classes_min num_rate_classes_max max_iterations epsilon edge_length prior_dist 
  transition_dist invariant_dist alpha_dist num_rate_classes_dist edge_time_dist file_stub_name;

(*Final matrix to hold all iterations*)

let final_matrix = Array.make_matrix alphabet_size alphabet_size 0.0 in
   (*prior_array.(alphabet_size - 1) <- prior_array.(alphabet_size - 1) /. 2.0;*)
   (*transition_array.(3) <- 0.1;
   transition_array.(6) <- 0.1;
   transition_array.(8) <- 0.1;
   transition_array.(9) <- 0.1;*)

for iteration = 0 to max_iterations -1 do
 (*Draw one replicate for all parameters*)
  let invariant_draw = (invariant_min +. (Random.float (invariant_max  -. invariant_min))) and
    alpha_draw = (alpha_min +. (Random.float (alpha_max  -. alpha_min))) and
    edge_draw = (Random.float edge_length) and
    rate_classes_draw = (num_rate_classes_min + (Random.int (1 + num_rate_classes_max  - num_rate_classes_min))) and
    prior_draw = dirichlet_variant prior_array and
    transition_draw = dirichlet_variant transition_array
  in
  (*Get a single matrix with drawn values*)
  let temp_matrix = Array.make_matrix alphabet_size alphabet_size 0.0 and
  rates = Numerical.gamma_rates alpha_draw alpha_draw rate_classes_draw in
  (*Rate classes and Invariant*)
  if (rate_classes_draw < 2) then rates.{0} <- 1.0;
   for i = 0 to rate_classes_draw - 1 do
     let subst_matrix = 
       try MlModel.m_gtr (ba_of_array1 prior_draw) transition_draw alphabet_size None 
       with | _ -> failwith "Matrix size and Priors are inconsistent"
     in
     let t_matrix = MlModel.compose_model subst_matrix (edge_draw *. rates.{i}) in
       matrix_increment_multiply temp_matrix t_matrix  ((1.0 -. invariant_draw) /. (float_of_int rate_classes_draw));
   done;
   for i = 0 to alphabet_size - 1 do
      temp_matrix.(i).(i) <- temp_matrix.(i).(i) +. invariant_draw;
   done;
   matrix_increment final_matrix temp_matrix;
  done;    
  matrix_divide final_matrix (float_of_int max_iterations);

  (*Output result matrices*)

  Printf.printf "Final matrix\n";
        for i = 0 to alphabet_size -1  do
          for j = 0 to alphabet_size -1  do
            Printf.printf "%f " final_matrix.(i).(j);
          done;
          Printf.printf "\n";
        done;
        Printf.printf "\n";

  Printf.printf "-log Final matrix\n";
        for i = 0 to alphabet_size -1  do
          for j = 0 to alphabet_size -1  do
            Printf.printf "%f " (-1.0 *. (log final_matrix.(i).(j)));
          done;
          Printf.printf "\n";
        done;
        Printf.printf "\n";

  let precision = 100000.0 in
  Printf.printf "Integerized -log Final matrix\n";
        for i = 0 to alphabet_size -1  do
          for j = 0 to alphabet_size -1  do
            Printf.printf "%d " (int_of_float (floor (precision *. -1.0 *. (log final_matrix.(i).(j)))));
          done;
          Printf.printf "\n";
        done;
        Printf.printf "\n";

;;
