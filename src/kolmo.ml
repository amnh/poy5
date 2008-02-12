(* A series of functions that reflect the concepts and function definitions from 
* An Introduction to Kolmogorov Complexity and its Applications by Li and
* Vitanyi *)

module Encodings = struct
    (* From page 13 l(x) *)
    let l x = int_of_float ((log (float_of_int (1 + x))) /. (log 2.)) 

    (* The binary representation of an integer *)
    let rec binary x = 
        if x < 0 then failwith "Illegal argument"
        else
            if x = 0 then ""
            else (binary (x lsr 1)) ^ (if (x land 1) = 1 then "1" else "0")

    (* From page 13, E_i(x) *)
    let rec e i x =
        match i with
        | 0 -> (String.make x '1') ^ "0"
        | i -> (e (i - 1) (l x)) ^ (binary x)

    let e_2 = e 2
end
