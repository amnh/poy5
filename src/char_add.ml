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

let () = SadmanOutput.register "Char_add" "$Revision: 1644 $"

(* $Id: char_add.ml 1644 2007-02-14 19:05:47Z andres $ *)

(* Each character is (from, to) *)
type e = int * int

type t = {
    tlen : int;
    (* costs? *)
    tweight : int array ref;
    tcode : int array;
    tmin : int array;
    tmax : int array;
    tcost : int;
    tmy_code : int;
}

type gen = int * int * int array * int
let randint a b = (Random.int (b - a)) + a
let rand_gen () =
    let min_randlen = 6 in
    let max_randlen = 12 in
    let min_randwidth = 12 in
    let max_randwidth = 40 in


    let len = randint min_randlen max_randlen in
    let width = randint min_randwidth max_randwidth in
    let weights = Array.make len 1 in
    (len, width, weights, Character.new_code ())

let make_rand (len, width, weights, code) =
    let v = Array.init len (fun _ -> randint 0 (width - 1)) in
    {tlen = len;
     tweight = ref weights;
     tcode = Array.init len (fun a -> a);
     tmin = v;
     tmax = v;
     tcost = 0;
     tmy_code = code;
    }

let color = Character.Brown
let code {tmy_code=c} = c

let median _ a b =
    (* We currently ignore the previous median.  We could use it as a reference
       in making the new one. *)

    assert (!(a.tweight) = !(b.tweight));

    let cost = ref 0 in
    let medmin = Array.make a.tlen 0 in
    let medmax = Array.make a.tlen 0 in

    for i = 0 to a.tlen - 1 do
        let mmin = max a.tmin.(i) b.tmin.(i) in
        let mmax = min a.tmax.(i) b.tmax.(i) in
        if mmin < mmax
        then begin
            medmin.(i) <- mmin;
            medmax.(i) <- mmax;
        end
        else begin
            medmin.(i) <- mmax;
            medmax.(i) <- mmin;
            cost := !cost + ((!(a.tweight)).(i) * (mmin - mmax));
        end
    done;
    {a with
         tmin = medmin;
         tmax = medmax;
         tcost = !cost;
    }
    
let median_3 _ _ _ _ =
    failwith "3-medians are unsupported for additive characters"

let distance a b =
    assert (a.tweight = b.tweight);

    let cost = ref 0 in
    let len = a.tlen in

    for i = 0 to len - 1 do
        let mmin = max a.tmin.(i) b.tmin.(i) in
        let mmax = min a.tmax.(i) b.tmax.(i) in
        if mmin > mmax
        then
            cost := !cost + ((!(a.tweight)).(i) * (mmin - mmax));
    done;
    (float_of_int !cost)

let parse i = []
let to_string t =
    let tostring i =
        let min = t.tmin.(i) in
        let max = t.tmax.(i) in
        if min = max
        then string_of_int min
        else string_of_int min ^ "-" ^ string_of_int max
    in
    let string = ref "" in
    for i = 0 to t.tlen - 1 do
        string :=
            (if !string = ""
             then "("
             else !string ^ " ") ^ tostring i
    done;
    !string ^ ")"

let (|>) a b =
    if a = 0
    then b
    else a

let compare_codes {tmy_code = c1} {tmy_code = c2} =
    compare c1 c2

let compare_data
        {tmy_code=c1;tmin=min1;tmax=max1}
        {tmy_code=c2;tmin=min2;tmax=max2} =
    (compare c1 c2) |> (compare min1 min2) |> (compare max1 max2)
