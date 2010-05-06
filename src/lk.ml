let main () = 
    if not (3 = Array.length Sys.argv) && not (2 = Array.length Sys.argv) then
        failwith "Usage: ./lk <matrixfile> <t>"
    else begin
        let data = FileStream.read_floatmatrix Sys.argv.(1)
        and branch_length = try float_of_string Sys.argv.(2) 
                            with | _ -> ~-.1.0
        in
        (* process read_floatmatrix output *)
        let priors,rate_matrix = match data with
            | hd::([])::tl | hd::tl ->
                let pi = Array.of_list hd
                and ms = Array.of_list (List.map Array.of_list tl) in
                (pi,ms)
            | _ -> failwith "Matrix file in an incorrect format."
        in
        (* create the matrix; and multiply through by priors *)
        let subst_matrix = 
            let alpha_size = Array.length priors in
            let subst = MlModel.m_file priors rate_matrix alpha_size in
            for i = 0 to alpha_size-1 do
                for j = (i+1) to alpha_size - 1 do
                    subst.{i,j} <- subst.{i,j} *. priors.(j);
                    subst.{j,i} <- subst.{j,i} *. priors.(i);
                done;
            done;
            MlModel.m_meanrate subst priors;
            subst
        in
        Printf.printf "Rate matrix multiplied priors w/ mean rate = 1\n";
        MlStaticCS.print_barray2 subst_matrix;

        let final_matrix = MlModel.compose_model subst_matrix branch_length in

        if branch_length < 0.0 then Printf.printf "Instantanious Rate Matrix:\n"
        else Printf.printf "Transition Probability Matrix for t = %f\n" branch_length;
        MlStaticCS.print_barray2 final_matrix
    end

let () = main ()
