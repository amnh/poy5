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
