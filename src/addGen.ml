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

(** {6 This is an implementation of the add.c in OCaml, and also not constrained
       to a the size of a char. This allows for general characters and
       functorization of the entire Additive module *)


(** {6 Types *)

type ct =
    {   min_state  : int array;
        max_state  : int array;
        cost       : int array;
        total      : float;
    }

(** {6 General Usability Functions *)

(** Comparison function; return first hit of difference **)
let compare_data x y =
    if (Array.length x.min_state) = (Array.length y.min_state) then begin
        let ret = ref 0 and i = ref 0 in
        while (!ret <> 0) && (!i < Array.length x.min_state) do begin
            ret := x.min_state.(!i) - y.min_state.(!i);
            incr i;
        end done;
        i := 0;
        while !ret <> 0 && !i < Array.length x.max_state do begin
            ret := x.max_state.(!i) - y.max_state.(!i);
            incr i;
        end done;
        !ret
    end else begin
        (Array.length x.min_state) - (Array.length y.min_state)
    end


(** Copy the contents of [x] into [y] **)
let copy x y =
    for i = 0 to (Array.length x.min_state) -1 do 
        y.min_state.(i) <- x.min_state.(i);
        y.max_state.(i) <- x.max_state.(i);
             y.cost.(i) <- x.cost.(i);
    done;
    ()

(** Create an exact copy of [x] **)
let clone x =
    {  min_state = Array.copy x.min_state;
       max_state = Array.copy x.max_state;
            cost = Array.copy x.cost;
           total = x.total; }

(** Convert data to string for debugging; (min,max,cost) **)
let to_string ct = 
    Array.fold_left
        (fun x y -> x ^ y) 
        ""
        (Array.mapi
            (fun i _ ->
                Printf.sprintf "(%d, %d, %d)"
                    ct.min_state.(i) ct.max_state.(i) ct.cost.(i))
            ct.min_state)


let print ct =
    Printf.printf "Min:\t";
    for i = 0 to (Array.length ct.min_state)-1 do
        Printf.printf "%d\t" ct.min_state.(i);
    done;
    Printf.printf "\nMax:\t";
    for i = 0 to (Array.length ct.min_state)-1 do
        Printf.printf "%d\t" ct.max_state.(i);
    done;
    Printf.printf "\nCst:\t";
    for i = 0 to (Array.length ct.min_state)-1 do
        Printf.printf "%d\t" ct.cost.(i);
    done;
    Printf.printf "\n%!";
    ()



(** {6 Accessor/Setter functions **)

(** Create a new type from array of mins and maxes *)
let create mins maxs =
    assert( (Array.length mins ) = (Array.length maxs) );
    { min_state = mins; 
      max_state = maxs;
           cost = Array.create (Array.length mins) 0;
          total = 0.0; }

(** Set the value at [i] in [x] to ([x_min],[x_max]) *)
let pos_set_state x i x_min x_max =
    assert( i < (Array.length x.min_state) );
    x.min_state.(i) <- x_min;
    x.max_state.(i) <- x_max;
    ()

(** get max value at [i] in [x] **)
let pos_get_max x i = x.max_state.(i)

(** get min value at [i] in [x] **)
let pos_get_min x i = x.min_state.(i)

(** get cost value at [i] in [x] **)
let pos_get_cost x i = float_of_int x.cost.(i)

(** return the total cost of the node **)
let median_cost x = x.total


(** {6 Median functions. *)

(** Find the median of a particular character; we modify res and return the cost
    for that characters; thus this can be used imperitively or functionally
    abstracted. **)
let rec add_item a_min b_min a_max b_max res i =
    if a_min <= b_min then begin
        if b_min <= a_max then begin (* intersection *)
            if b_max > a_max then begin (* not contained *)
                res.max_state.(i)   <- a_max;
                res.min_state.(i)   <- b_min;
                res.cost.(i)        <- 0;
            end else begin  (* full containment *)
                res.max_state.(i)   <- b_max;
                res.min_state.(i)   <- b_min;
                res.cost.(i)        <- 0;
            end
        end else begin (* no intersection *)
            res.max_state.(i) <- b_min;
            res.min_state.(i) <- a_max;
            res.cost.(i)      <- (b_min - a_max);
        end;
        res.cost.(i)
    end else begin
        add_item b_min a_min b_max a_max res i 
    end

(** General median function for characters; non-functional style. *)
let mediani a b r =
    let total_cost = ref 0 in
    for i = 0 to (Array.length a.min_state) - 1 do 
        let add_cost = 
            add_item a.min_state.(i) b.min_state.(i) a.max_state.(i) 
                     b.max_state.(i) r i
        in
        total_cost := !total_cost + add_cost;
    done;
    float_of_int !total_cost

(** General median function for characrers; abstracted non-functional for
    functional useage *)
let median a b = 
    let r = clone a in
    let m = mediani a b r in
    { r with total = m }

(** To make the interface consistant; although it should be noted that the total
    cost of the node will be incorrect. *)
let mediani a b r = ignore (mediani a b r)

(** median of three *)
let median_3i p n c1 c2 r =
    (* check if x is contained in y *)
    let rec is_contained xmin xmax ymin ymax =
        (xmin >= ymin) && (xmax <= ymax)
    (* add containment to result; x is contained in y *)
    and add_contained xmin xmax i =
        r.max_state.(i) <- xmax;
        r.min_state.(i) <- xmin
    (* find the full union; largest interval *)
    and full_union xmin ymin xmax ymax =
        if xmin <= ymin then
            if xmax <= ymax then (xmin,ymax) else (xmin,xmax)
        else
            if xmax <= ymax then (ymin,ymax) else (ymin,xmax)
    (* Check if we are dealing with an intersection *)
    and is_intersection xmin xmax ymin ymax =
        if xmin <= ymin then
            not (xmax >= ymin)
            else is_intersection ymin ymax xmin xmax
    (* apply the intersection to the result *)
    and add_intersection xmin xmax ymin ymax i =
        if xmin <= ymin then begin
            if ymax >= xmax then begin
                r.min_state.(i) <- ymin;
                r.max_state.(i) <- xmax
            end else begin
                r.min_state.(i) <- ymin;
                r.max_state.(i) <- ymax
            end
        end else begin
            add_intersection ymin ymax xmin xmax i
        end
    (* add non-intersecting/containing union to result *)
    and add_union xmin xmax ymin ymax i =
        if xmax <= ymax then begin
            r.min_state.(i) <- xmax;
            r.max_state.(i) <- ymin
        end else begin
            r.min_state.(i) <- ymin;
            r.max_state.(i) <- xmax
        end
    in
    for i = 0 to (Array.length c1.max_state)-1 do
        if not (is_contained p.min_state.(i) p.max_state.(i)
                             n.min_state.(i) n.max_state.(i)) then
            let xmin,xmax = full_union c1.min_state.(i) c2.min_state.(i)
                                       c1.max_state.(i) c2.max_state.(i) in
            if not (is_intersection xmin xmax p.min_state.(i) p.max_state.(i)) then
                if p.min_state.(i) < xmin then
                    add_union p.min_state.(i) p.max_state.(i) xmin xmax i
                else
                    add_union xmin xmax p.min_state.(i) p.max_state.(i) i
            else
                add_intersection xmin xmax p.min_state.(i) p.max_state.(i) i
        else
            add_contained p.min_state.(i) p.max_state.(i) i
    done;
    ()

(* median of three; result should replace [n] in the returned data *)
let median_3 p n c1 c2 =
    let r = clone n in
    let () = median_3i p n c1 c2 r in
    r

(* Basic union function for calculating a median; full union medians *)
let full_union a b r =
    for i = 0 to (Array.length a.min_state)-1 do
        let rmin,rmax =
            if a.min_state.(i) <= b.min_state.(i) then
                if a.max_state.(i) <= b.max_state.(i)
                    then (a.min_state.(i),b.max_state.(i))
                    else (a.min_state.(i),a.max_state.(i))
            else
                if a.max_state.(i) <= b.max_state.(i)
                    then (b.min_state.(i),b.max_state.(i))
                    else (b.min_state.(i),a.max_state.(i))
        in
        r.min_state.(i) <- rmin;
        r.max_state.(i) <- rmax
    done;
    ()


(** {6 Distance functions. *)

let set_minimum d s =
    let total_cost = ref 0 in
    for i = 0 to (Array.length d.min_state) - 1 do
        if d.cost.(i) > s.cost.(i) then d.cost.(i) <- s.cost.(i);
        total_cost := !total_cost + d.cost.(i)
    done;
    {d with total = float_of_int !total_cost; }

let distance a b = 
    (median a b).total

let distance_2 n a b = 
    (set_minimum (median n a) (median n b)).total

let distance_median a b = 
    let m = median a b in m.total, m

