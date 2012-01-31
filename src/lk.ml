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

let (-->) a b = b a
let ba_of_array1 x = Bigarray.Array1.of_array Bigarray.float64 Bigarray.c_layout x
and ba_of_array2 x = Bigarray.Array2.of_array Bigarray.float64 Bigarray.c_layout x

let print_barray1 a =
    for i = 0 to (Bigarray.Array1.dim a)-1 do
        Printf.printf "%2.10f\t" a.{i};
    done; Printf.printf "\n"; ()

and print_barray2 a =
    for i = 0 to (Bigarray.Array2.dim1 a)-1 do
        for j = 0 to (Bigarray.Array2.dim2 a)-1 do
            Printf.printf "%2.10f\t" a.{i,j};
        done; Printf.printf "\n";
    done; Printf.printf "\n"; ()

let adjust_sub_matrix subst_matrix priors alpha_size = ()
(*   let column_total = Array.create alpha_size 0.0 in
	for i = 0 to alpha_size - 1 do
           for j = 0 to alpha_size - 1 do
               subst_matrix.{i,j} <- subst_matrix.{i,j} /. priors.(j);
               if (i <> j) then column_total.(j) <- column_total.(j) +. subst_matrix.{i,j};
           done;
        done;
        for i = 0 to alpha_size - 1 do
           subst_matrix.{i,i} <- -1.0 *. column_total.(i);
        done;
*)
;;

(*mutiply by priors to make symmetrical and adjust such that rows sum to 1*)
let make_symmetrical subst_matrix priors alpha_size =
    let row_total = Array.create alpha_size 0.0 in
    for i = 0 to alpha_size - 1 do
        for j = 0 to alpha_size - 1 do
            if (i <> j) then subst_matrix.(i).(j) <- subst_matrix.(i).(j) *.
            priors.(i) *. (float_of_int alpha_size); 
            row_total.(i) <- row_total.(i) +. subst_matrix.(i).(j);
            Printf.printf "%f " subst_matrix.(i).(j);
        done;
        Printf.printf "\n" ;
     done;
    Printf.printf "\n";

     (*for i = 0 to alpha_size - 1 do
        for j = 0 to alpha_size - 1 do
            subst_matrix.(i).(j) <- subst_matrix.(i).(j) /. row_total.(i);
        done;
     done;
     for i = 0 to alpha_size - 1 do
        for j = 0 to alpha_size - 1 do
            Printf.printf "%f " subst_matrix.(i).(j);
        done;
        Printf.printf " total %f \n" row_total.(i);
     done;
 Printf.printf "\n";*)

 
;;

let main () = 
    if not (3 = Array.length Sys.argv) && not (2 = Array.length Sys.argv) && not
    (6 = Array.length Sys.argv) && not (8 = Array.length Sys.argv) && not (9 = Array.length Sys.argv) then
        failwith "Usage: ./lk <matrixfile> <t> <step> <distribution> <dist param> <rate classes> <alpha> <theta>"
    else begin
        let data = FileStream.read_floatmatrix Sys.argv.(1)
        and branch_length_max = try float_of_string Sys.argv.(2) 
                            with | _ -> ~-.1.0
        and steps = try int_of_string Sys.argv.(3) 
                            with | _ -> ~-1 
        and distro = try Sys.argv.(4)
                            with | _ -> "no_distro"
        and distro_param = try float_of_string Sys.argv.(5)
                            with | _ -> ~-.1.0
        and nclasses = try int_of_string Sys.argv.(6)
                            with | _ -> ~-1
        and alpha = try float_of_string Sys.argv.(7)
                            with | _ -> ~-.1.0
        and theta = try float_of_string Sys.argv.(8)
                            with | _ -> ~-.0.0
        in
        (* process read_floatmatrix output *)
        let priors,coefficients = match data with
            | hd::co::[] ->
                let pi = Array.of_list hd and co = Array.of_list co in
                (pi,co)
            | _ -> failwith "Matrix file in an incorrect format."
        in
        (* create the matrix; and multiply through by priors *)
        let alpha_size = Array.length priors in
        let integrated_matrix = Array.make_matrix alpha_size alpha_size 0.0
	and final_matrix = Array.make_matrix alpha_size alpha_size 0.0
        in
	(*MlModel.gamma_rates <alpha:float> <beta:float> <categories:int>*)
        let rates = Numerical.gamma_rates alpha alpha nclasses in
        for i = 0 to nclasses -1  do
          Printf.printf "rate %d = %f\n" i rates.{i};
        done;
        for i = 0 to alpha_size -1  do
          Printf.printf "prior %d = %f\n" i priors.(i);
        done;
        let subst_matrix = 
            try MlModel.m_gtr (ba_of_array1 priors) coefficients alpha_size None 
            with | _ -> failwith "Matrix size and Priors are inconsistent"
        in
        (* Adjust subst_matrix with priors for symmetry*)
        adjust_sub_matrix subst_matrix priors alpha_size;
        
        (*let final_matrix = MlModel.compose_model subst_matrix
        * branch_length_max in*)
        (*print_barray2 subst_matrix;*)

        let interval = branch_length_max /. (float_of_int steps)  
        and bl = ref 0.0 in
        for iterations = 0 to steps - 1 do
            if (nclasses < 1) then 
                begin
                  let s_m = MlModel.m_gtr (ba_of_array1 priors) coefficients alpha_size None in
		  (* Adjust subst_matrix with priors for symmetry*)
                  adjust_sub_matrix s_m priors alpha_size;
                  let t_matrix = MlModel.compose_model s_m !bl in
                  for i = 0 to alpha_size - 1 do
                      for j = 0 to alpha_size - 1 do
                         final_matrix.(i).(j) <- (1.0 -. theta) *. t_matrix.{i,j};
                         if (i=j) then  final_matrix.(i).(j) <-  final_matrix.(i).(j) +. theta;
                      done;
                   done;
                 end
            else
               begin
                 let temp_matrix = Array.make_matrix alpha_size alpha_size 0.0 in
                 for k = 0 to nclasses -1 do
                     let s_m = MlModel.m_gtr (ba_of_array1 priors) coefficients alpha_size None in
                      (* Adjust s_m with priors for symmetry*)
		      adjust_sub_matrix s_m priors alpha_size;
		  let t_matrix = MlModel.compose_model s_m (rates.{k} *. !bl) in
                     (*Printf.printf "rate %d %f bl %f\n" k rates.{k} (rates.{k} *. !bl);*)
                     for i = 0 to alpha_size - 1 do
                       for j = 0 to alpha_size - 1 do
                        temp_matrix.(i).(j) <- temp_matrix.(i).(j) +. ((t_matrix.{i,j} /. (float_of_int nclasses)) *. (1.0 -. theta));
                        (* Printf.printf "fm %f tm %f \n"  temp_matrix.(i).(j) (t_matrix.{i,j} /. (float_of_int nclasses));*)
                       done;
                     done;
                  done;
                  for i = 0 to alpha_size - 1 do
                     for j = 0 to alpha_size - 1 do
                         final_matrix.(i).(j) <- temp_matrix.(i).(j);
			 if (i=j) then final_matrix.(i).(j) <- final_matrix.(i).(j) +. theta;
                     done;
                  done;
               end;
            (*branch_length_max in *)
            let branch_prob = ref 1.0 in
            if ( distro = "uniform") then branch_prob := ( interval /. branch_length_max)
            else if (distro = "exponential") then  branch_prob := (interval *.
                distro_param *. exp(-1.0 *. (!bl +. (interval /. 2.0)) *. distro_param))
            else Printf.printf "Unkown distribution\n";

            for i = 0 to alpha_size - 1 do
                for j = 0 to alpha_size - 1 do
                    integrated_matrix.(i).(j) <- integrated_matrix.(i).(j) +.
                    (final_matrix.(i).(j)  *. !branch_prob);
                done;
            done;
            bl := !bl +. interval; 
            (*Printf.printf "bl %f bp %f int %f \n" !bl !branch_prob interval;*)
        done;
        Printf.printf "Rate matrix multiplied priors w/ mean rate = 1\n";
        print_barray2 subst_matrix;
        if branch_length_max < 0.0 then Printf.printf "Instantanious Rate Matrix:\n"
        else Printf.printf "Transition Probability Matrix for t = %f\n"
        branch_length_max;
        let final_matrix2 = MlModel.compose_model subst_matrix branch_length_max in
        print_barray2 final_matrix2;

	(*Make symmetrical by normalizing by priors*)
	 (*for i = 0 to alpha_size - 1 do
            for j = 0 to alpha_size - 1 do
                 integrated_matrix.(i).(j) <- integrated_matrix.(i).(j) /. priors.(j);
                 done;
            done;*)

        (*Print integrated matrix*)
        Printf.printf "Integrated rate matrix 0 to %f by %d steps of %f gamma classes %d alpha-beta %f theta %f\n" !bl
        steps interval nclasses alpha theta;
        for i = 0 to alpha_size - 1 do
            for j = 0 to alpha_size - 1 do
                 Printf.printf "%f " integrated_matrix.(i).(j);
                 done;
                 Printf.printf "\n";
            done;
            Printf.printf "\n";

        (*Print log weight matrix*)
        Printf.printf "Natural Log weight matrix\n";
        (*make_symmetrical integrated_matrix priors alpha_size;*) 
        for i = 0 to alpha_size - 1 do
            for j = 0 to alpha_size - 1 do
                if (i <> j) then Printf.printf "%f " (-1.0*. (log (integrated_matrix.(i).(j))
                 +. (log (priors.( i)))))
                else  Printf.printf "%f " (-1.0*. (log
                (integrated_matrix.(i).(j))));
                 done;
                 Printf.printf "\n";
            done;
            Printf.printf "\n";

        Printf.printf "Integerized Natural Log weight matrix\n";
        let precision = 10000.0 in
        (*make_symmetrical integrated_matrix priors alpha_size;*) 
        for i = 0 to alpha_size - 1 do
            for j = 0 to alpha_size - 1 do
                if (i <> j) then Printf.printf "%d " (int_of_float (floor((-1.0
                *. precision *. (log (integrated_matrix.(i).(j))
                 +. (log (priors.( i))))))))
                else  Printf.printf "%d " (int_of_float (floor((-1.0 *.
                precision *. (log
                (integrated_matrix.(i).(j)))))));
                 done;
                 Printf.printf "\n";
            done;
            Printf.printf "\n";


   end

let () = main ()
