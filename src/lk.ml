
let (-->) a b = b a

let main () = 
    if not (3 = Array.length Sys.argv) && not (2 = Array.length Sys.argv) then
        failwith "Usage: ./lk <matrixfile> <t>"
    else begin
        let data = FileStream.read_floatmatrix Sys.argv.(1)
        and branch_length = try float_of_string Sys.argv.(2) 
                            with | _ -> ~-.1.0
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
        let subst_matrix = 
            try MlModel.m_gtr priors coefficients alpha_size
            with | _ -> failwith "Matrix size and Priors are inconsistent"
        in
        let final_matrix = MlModel.compose_model subst_matrix branch_length in

        Printf.printf "Rate matrix multiplied priors w/ mean rate = 1\n";
        MlStaticCS.print_barray2 subst_matrix;
        if branch_length < 0.0 then Printf.printf "Instantanious Rate Matrix:\n"
        else Printf.printf "Transition Probability Matrix for t = %f\n" branch_length;
        MlStaticCS.print_barray2 final_matrix
    end

let () = main ()
