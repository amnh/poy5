(* POY 5.1.1. A phylogenetic analysis program using Dynamic Homologies.       *)
(* Copyright (C) 2014 Andrés Varón, Lin Hong, Nicholas Lucaroni, Ward Wheeler,*)
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

(* $Id: array_ops.ml 2871 2008-05-23 17:48:34Z andres $ *)
let () = SadmanOutput.register "Array_ops" "$Revision: 3649 $"


exception Empty

let use_array_append_from_caml3 = true

(** [array_append_caml3 arr1 arr2] return Array.append arr1 arr2.
* this is a walk around of seg fault problem causing by Array.append of Ocaml4.0.X. 
* here we create our version of Array.append, does the job just like append function from Ocaml3.12.X. *)
let array_append a1 a2 =
  if use_array_append_from_caml3 then begin
      let l1 = Array.length a1 and l2 = Array.length a2 in
      if l1 = 0 && l2 = 0 then [||] else begin
          (*Array.get calls safe_get from array.c, this is not what they did, they call unsafe_get in Ocaml3.12.X, but there is
          * no entrance for unsafe_get in Array.mli*)
        let r = Array.make (l1 + l2) (Array.get (if l1 > 0 then a1 else a2) 0) in
        for i = 0 to l1 - 1 do Array.set r i (Array.get a1 i) done;
        for i = 0 to l2 - 1 do Array.set r (i + l1) (Array.get a2 i) done;
        r
      end;
  end
  else 
      Array.append a1 a2



let rec _calculate_size a len lst sum = match len with 
    | (-1) -> sum, lst
    | _ -> let s = Array.length (a.(len)) in
            _calculate_size a (len - 1) (s :: lst) (sum + s)

let rec _transfer_cont a tgt lst cura curtgt =
    match lst with
    | [] -> tgt
    | 0 :: t -> _transfer_cont a tgt t (cura + 1) curtgt
    | h :: t ->
            let new_curtgt = curtgt + h in
            for i = curtgt to new_curtgt - 1 do
                tgt.(i) <- a.(cura).(i - curtgt);
            done;
            _transfer_cont a tgt t (cura + 1) new_curtgt

(* Returns the first element in any of the non empty arrays in the array array
* a. *)
let rec _get_first a lst count =
    match lst with
    | 0 :: t -> _get_first a t (count + 1)
    | _ :: t -> a.(count).(0)
    | [] -> raise (Empty)

(* Note that tgt should have the same size as the sum of all the
 arrays pointed by a. The lst, cura and curtgt are variables for tail
 recursion. *)
let flatten_array a =
    let len = Array.length a in
    let size, lst = _calculate_size a (len - 1) [] 0 in
    match size with
    | 0 -> [||]
    | _ ->
            let sample = _get_first a lst 0 in
            let tmp = Array.make size sample in
            _transfer_cont a tmp lst 0 0

let rec _rec_fold_min f v a l it max =
    if max = it then (Array.of_list l)
    else begin
        let res, cost = f a.(it) in
        if cost > v then _rec_fold_min f v a l (it + 1) max
        else if cost = v then _rec_fold_min f v a (res :: l) (it + 1) max
        else _rec_fold_min f cost a [res] (it + 1) max
    end

let rec _rec_fold_max f v a l it max =
    if max = it then (Array.of_list l)
    else begin
        let res, cost = f a.(it) in
        if cost < v then _rec_fold_max f v a l (it + 1) max
        else if cost = v then _rec_fold_max f v a (res :: l) (it + 1) max
        else _rec_fold_max f cost a [res] (it + 1) max
    end

let int_f func c elem = 
    let r = func elem in
    r, (c r)

let fold_min (f : 'a -> 'b) (c : 'b -> int) a =
    _rec_fold_min (int_f f c) max_int a [] 0 (Array.length a)

let fold_min_p f c a min =
    _rec_fold_min (int_f f c) min a [] 0 (Array.length a)


let fold_max f c a =
    _rec_fold_max (int_f f c) min_int a [] 0 (Array.length a)

let fold_max_p f c a max =
    _rec_fold_max (int_f f c) max a [] 0 (Array.length a)

let map_ip f a =
    let len = Array.length a in
    if 0 < len then begin
        for i = 0 to len - 1 do
            a.(i) <- f (a.(i));
        done;
    end else ()

let mapi_ip f a =
    let len = Array.length a in
    if 0 < len then begin
        for i = 0 to len - 1 do
            a.(i) <- f i (a.(i));
        done;
    end else ()

let randomize ar = (* Fisher-Yates *)
    let l = Array.length ar - 1 in
    for i = 0 to l do
        let rnd = i + Random.int (l - i + 1) in
        let tmp = ar.(i) in
        ar.(i) <- ar.(rnd);
        ar.(rnd) <- tmp;
    done;;

let mem ar a =
    let rec short_circuit i =
        if i = (Array.length ar) 
            then false 
        else if ar.(i) = a 
            then true 
            else short_circuit (i+1)
    in
    short_circuit 0

let filter f arr =
    let res = Array.fold_right (fun x acc -> 
        if (f x) then x :: acc else acc) arr [] in
    Array.of_list res

let map_2 f a b = 
    assert (Array.length a = Array.length b);
    Array.init (Array.length a) (fun x -> f a.(x) b.(x)) 

let map_3 f a b c = 
    assert (Array.length a = Array.length b);
    assert (Array.length c = Array.length b);
    Array.init (Array.length a) (fun x -> f a.(x) b.(x) c.(x)) 

let map_4 f a b c d = 
    assert (Array.length a = Array.length b);
    assert (Array.length c = Array.length b);
    assert (Array.length d = Array.length b);
    Array.init (Array.length a) (fun x -> f a.(x) b.(x) c.(x) d.(x)) 

let map_5 f a b c d e = 
    assert (Array.length a = Array.length b);
    assert (Array.length c = Array.length b);
    assert (Array.length d = Array.length b);
    assert (Array.length e = Array.length b);
    Array.init (Array.length a) (fun x -> f a.(x) b.(x) c.(x) d.(x) e.(x)) 

let fold_right_2 f acc a b =
    let acc = ref acc in
    assert (Array.length a = Array.length b);
    for i = (Array.length a) - 1 downto 0 do
        acc := f !acc a.(i) b.(i);
    done;
    !acc

let fold_right_3 f acc a b c =
    let acc = ref acc in
    assert (Array.length a = Array.length b);
    assert (Array.length c = Array.length b);
    for i = (Array.length a) - 1 downto 0 do
        acc := f !acc a.(i) b.(i) c.(i);
    done;
    !acc

let fold_right_4 f acc a b c d =
    let acc = ref acc in
    assert (Array.length a = Array.length b);
    assert (Array.length c = Array.length b);
    assert (Array.length d = Array.length b);
    for i = (Array.length a) - 1 downto 0 do
        acc := f !acc a.(i) b.(i) c.(i) d.(i);
    done;
    !acc

let split size arr =
    let len = Array.length arr in
    let remainder = len - ((len / size) * size) in
    let fraction = len / size in
    Array.init size (fun pos ->
        if remainder = 0 then
            Array.sub arr (pos * fraction) fraction
        else if remainder > pos then
            Array.sub arr ((pos * fraction) + pos)
            (fraction + 1)
        else
            Array.sub arr ((pos * fraction) + remainder)
            fraction)

let fold_righti f init ray =
    let acc = ref init in
    for i = (Array.length ray) - 1 downto 0 do
        acc := f i ray.(i) !acc
    done;
    !acc

let is_identical2 arr1 arr2 =
    let len1 = Array.length arr1 and len2 = Array.length arr2 in
    if (len1<>len2) then
        0
    else begin
        let sign = ref 1 in 
        let idx = ref (len1-1) in
        while (!sign = 1)&&( !idx > 0) do
            sign := 
                if (arr1.(!idx) = arr2.(!idx)) then 1
                else 0 ;
            idx := !idx -1;
        done;
        !sign
    end

let is_identical3 arr1 arr2 arr3 =
    let len1 = Array.length arr1 and len2 = Array.length arr2 
    and len3 = Array.length arr3 in
    if ((len1<>len2)||(len2<>len3)) then begin
        0
    end else begin
        let sign = ref 1 in 
        let idx = ref (len1-1) in
        while (!sign = 1)&&( !idx > 0) do
            sign := 
                if (arr1.(!idx) = arr2.(!idx))&& (arr3.(!idx) = arr2.(!idx)) then 1
                else 0 ;
            idx := !idx -1 ;
        done;
        !sign
    end


IFDEF USE_PARMAP THEN
    let fill_symmetric_square_matrix ?(status=None) f src des =
        let states = Array.length des in
        assert( (states > 0) && (states = (Array.length des.(0))) );
        let xy_of_i i =
            let rec xy_of_i s y i =
                if i > s
                    then xy_of_i (s+1) (y+1) (i-(s+1))
                    else (i,y)
            in
            xy_of_i 0 0 i
        in
        let tmp_array = Array.make ((states * (states+1))/2) des.(0).(0) in
        let tmp_array =
            Parmap.array_float_parmapi
                (fun i _ -> let x,y = xy_of_i i in f x y src.(x) src.(y))
                tmp_array
        in
        Array.iteri
            (fun i c -> let x,y = xy_of_i i in des.(x).(y) <- c; des.(y).(x) <- c)
            tmp_array;
        ()

ELSE
    let fill_symmetric_square_matrix ?(status=None) f src des =
        let report stride =
            match status with
            | Some x -> 
                let status = Status.create x (Some 100) "percent complete"
                and each = (1.0 /. (float_of_int stride)) *.  100.0 in
                let curr = ref 0.0 in
                (fun incr ->
                    curr := !curr +. (each *. (float_of_int incr));
                    Status.achieved status (int_of_float !curr);
                    Status.full_report status)
            | None ->
                (fun _ -> ())
        in
        let states = Array.length des in
        assert( states = (Array.length des.(0)) );
        assert( states = (Array.length src) );
        let report = report ((states-1)*(states-1)) in
        for x = 0 to states-1 do
            for y = x to states-1 do
                let res = f x y src.(x) src.(y) in
                des.(x).(y) <- res;
                des.(y).(x) <- res;
            done;
            report (2*(states-2-x));
        done;
        ()
END
