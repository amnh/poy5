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

let () = SadmanOutput.register "FingerPrint" "$Revision: 1644 $"

exception Illegal_argument
exception Not_match

let sigma_p a b p f = (f a b) mod p

let plus a b p = sigma_p a b p (+)

let minus a b p = sigma_p a b p (-)

let times a b p = sigma_p a b p ( * )

module type SERIALIZABLE = sig
    type ser
    type v
    val equal : v -> v -> bool
    val get_pos : ser -> int -> v
    val to_int : v -> int
    val of_int : int -> v
    val length : ser -> int
    val sub : ser -> int -> int -> ser
end

module Make = functor (Elt : SERIALIZABLE) -> struct

        let fp t q i n asz =
            assert (Elt.length t >= i + n - 1);
            let rec builder acc pos max =
                if max > pos then begin
                    let it = Elt.to_int (Elt.get_pos t pos) in
                    builder (plus (times acc asz q) it q) (pos + 1) max
                end else acc
            in
            builder 0 i (i + n)

        let complete_fp a q asz = fp a q 0 (Elt.length a) asz

        let pow v n q = 
            let rec builder acc n =
                if n <> 0 then builder ((v * acc) mod q) (n - 1)
                else acc
            in
            builder 1 n

        let res a b = 
            let tmp = a / b in
            a - (tmp * b)

        let fp_incremental t q i n cur asz = 
            let new_n = i + n - 1 in
            assert (Elt.length t > new_n);
            (* First remove what we get from the previous value *)
            let pos = Elt.to_int (Elt.get_pos t i) in
            let removed = 
                (q + (minus cur (times pos (pow asz (n - 1) q) q) q)) mod q in
            let shifted = times removed asz q 
            and new_pos = Elt.to_int (Elt.get_pos t (new_n + 1)) in
            plus shifted new_pos q

        let fp_match t p asz q =
            let n = Elt.length t 
            and m = Elt.length p in
            let p_fp = complete_fp p q asz in
            let rec compare st it = 
                if it = m then true
                else if Elt.equal (Elt.get_pos t st) (Elt.get_pos p it) then 
                    compare (st + 1) (it + 1) 
                else false
            in
            let rec finder cur_pos cur_fp =
                if cur_fp = p_fp then begin
                    if compare cur_pos 0 then cur_pos
                    else if cur_pos < (n - m) then begin
                        let next_fp = 
                            fp_incremental t q cur_pos m cur_fp asz in
                        finder (cur_pos + 1) next_fp
                    end else raise Not_match
                end else if cur_pos < (n - m) then begin
                    let next_fp = fp_incremental t q cur_pos m cur_fp asz in
                    finder (cur_pos + 1) next_fp
                end else raise Not_match
            in
            let initial_fp = fp t q 0 m asz in
            finder 0 initial_fp

        let do_match t p asz = 
            let n = Elt.length t 
            and m = Elt.length p in
            let q = Primes.Probable.next_prime_lt (n * m * m) in
            fp_match t p asz q

            let matcher asz q = 
                fun t p -> fp_match t p asz q

            let fp_incr m asz p = 
                fun t it cur -> fp_incremental t p it m cur asz
    end

module MakeMult = functor (Elt : SERIALIZABLE) -> struct
    module EltFp = Make (Elt)
    
    let fp s pns i n a =
        List.map (fun x -> EltFp.fp s x i n a) pns

    let complete_fp s pns n a =
        List.map (fun x -> EltFp.complete_fp s x a) pns

    let fp_incremental s pns i n curr a =
        let rec calculator pns curr acc =
            match pns, curr with
            | hp :: tp, hc :: tc ->
                    let new_it = EltFp.fp_incremental s hp i n hc a in
                    calculator tp tc (new_it :: acc)
            | [], [] -> List.rev acc
            | _, _ -> failwith "Lists of different length"
        in
        calculator pns curr []
end

(* Some generics that we will be using *)
module SerGenArray = struct
    type v = int
    type ser = v array
    let equal a b = 0 = compare a b 
    let get_pos a p = a.(p)
    let to_int v = v 
    let of_int v = v 
    let length = Array.length
    let sub a b c =  try Array.sub a b c with | _ -> raise Illegal_argument
end

module SerString = struct
    type v = char
    type ser = string
    let equal a b = a = b 
    let get_pos a p = a.[p]
    let to_int = Char.code
    let of_int = Char.chr
    let length = String.length
    let sub a b c = try String.sub a b c with | _ -> raise Illegal_argument
end

(* A matcher for arrays of integers *)
module IntArrayFp = Make (SerGenArray)

(* A matcher for strings *)
module StringFp = Make (SerString)
