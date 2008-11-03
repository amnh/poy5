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

let () = SadmanOutput.register "Primes" "$Revision: 1644 $"

module Probable = struct

    let calculate_r n =
        let m = n -  1 in
        let rec shifter n cnt =
            if 0 <> n land 1 then n, cnt
            else shifter (n lsr 1) (cnt + 1)
        in
        shifter m 0

    let modpower a r n =
        let a = Int64.of_int a
        and r = Int64.of_int r 
        and n = Int64.of_int n in
        let rec builder acc r =
            if r = 1L then Int64.rem (Int64.mul a acc) n
            else builder (Int64.rem (Int64.mul a acc) n) (Int64.pred r)
        in
        Int64.to_int (builder 1L r)

    let rec iterative_test i xim1 n nm1 t =
        if i > t then xim1 = 1
        else
            let xi = modpower xim1 2 n in
            if xi = 1 && xim1 <> 1 && xim1 <> nm1 then false
            else iterative_test (i + 1) xi n nm1 t

    let karp_rabi n =
        let a = 2 + (Random.int (n - 2)) in
        let r, s = calculate_r n in
        let res = modpower a r n in
        if res <> 1 && res <> n - 1 then 
            iterative_test 1 res n (n - 1) s
        else true

    let rec is_prime ?(tests = 20) n =
        (* Check by hand the cases smaller than 10 to avoid an endless loop in
        * the Karp-Rabi test *)
        if n < 10 then 
            n = 2 || n = 3 || n = 5 || n = 7 
        else if 0 = n mod 2 then false
        else
            (* The Karp-Rabi test 
            let rec test_primality it max =
                if it > max then true
                else 
                    if 0 = n mod it then false
                    else test_primality (it + 2) max
            in
            test_primality 3 (int_of_float (sqrt (float_of_int n)))
            *)
            let rec tester it =
                if it = 0 then true
                else karp_rabi n && (tester (it - 1))
            in 
            tester tests

    let rec next_prime_gt m =
        if is_prime m then m
        else  next_prime_gt (m + 1)

    let next_prime_lt m =
        if m < 2 then failwith "Primes.next_prime_lt" 
        else
            let rec next_prime_lt m =
                if is_prime m then m
                else next_prime_lt (m - 1)
            in
            next_prime_lt m

end
