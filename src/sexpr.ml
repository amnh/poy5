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

let () = SadmanOutput.register "Sexpr" "$Revision: 1810 $"

type 'a t = [ `Empty | `Set of 'a t list | `Single of 'a ]

let rec fold_left f acc = function
    | `Empty -> acc
    | `Single tree -> f acc tree
    | `Set lst -> List.fold_left (fold_left f) acc lst

let rec fold_right f acc = function
    | `Empty -> acc
    | `Single tree -> f acc tree
    | `Set lst -> List.fold_right (fun x y -> fold_right f y x) lst acc 

let rec map f = function
    | `Empty -> `Empty
    | `Single tree -> `Single (f tree)
    | `Set lst -> `Set (List.map (map f) lst)


let rec map_fold_left f acc = function
    | `Empty -> (`Empty : 'a t), acc
    | `Single tree -> 
          let acc, tree = f acc tree in 
          (`Single tree : 'a t), acc
    | `Set lst ->
          let acc, lst = 
              List.fold_right (fun item (acc, lst) ->
                                   let tree, acc = map_fold_left f acc item in
                                   acc, (tree :: lst)
                              )  lst (acc, [])
          in
          (`Set lst : 'a t), acc


let cardinal t = fold_right (fun acc _ -> succ acc) 0 t

let map_feedback cb f sexpr =
    let count = ref 0 in
    let f t =
        incr count;
        cb !count;
        f t in
    map f sexpr

let fold_feedback cb f init sexpr =
    let count = ref 0 in
    let f t =
        incr count;
        cb !count;
        f t in
    fold_left f init sexpr

let make_single a = `Single a 
let singleton = make_single
let of_list = function
    | [] -> `Empty
    | [a] -> make_single a
    | lst -> `Set (List.map make_single lst)

let to_list t = 
    List.rev (fold_left (fun acc x -> x :: acc) [] t)

let rec first = function
    | `Single x -> x
    | `Empty -> failwith "Sexpr.first"
    | `Set lst -> 
            match lst with
            | h :: t ->
                    begin try first h with
                    | Failure "Sexpr.first" -> first (`Set t)
                    end
            | [] -> failwith "Sexpr.first"

let rec length = function
    | `Empty -> 0
    | `Single _ -> 1
    | `Set lst ->
            List.fold_left (fun acc n -> acc + length n) 0 lst

let rec all_to_all status f a b =
    match a, b with
    | `Empty, _ | _, `Empty -> 
            `Empty
    | `Single x, `Single y -> 
            let ach = Status.get_achieved status in
            Status.full_report ~adv:(ach + 1) status;
            `Single (f x y)
    | `Single x, `Set y ->
            let mapper = fun x -> all_to_all status f a x in
            `Set (List.map mapper y)
    | `Set x, `Single y ->
            let mapper = fun x -> all_to_all status f x b in
            `Set (List.map mapper x)
    | `Set x, `Set y ->
            let folder = fun x acc ->
                (all_to_all status f x b :: acc)
            in
            `Set (List.fold_right folder x [])

let all_to_all f a b =
    let status = 
        let lena = length a
        and lenb = length b in
        Status.create "All pairs" (Some (lena * lenb)) ""
    in
    let res = all_to_all status f a b in
    Status.finished status;
    res

let shallow_all_to_all f a b =
    match a, b with
    | `Empty, _ | _, `Empty -> `Empty
    | `Single _, `Single _
    | `Single _, `Set _
    | `Set _, `Single _ -> `Single (f a b)
    | `Set x, `Set _ ->
            `Set (List.map (fun c -> (`Single (f c b) : 'b t)) x)

let rec rev = function
    | `Empty | `Single _ as x -> x
    | `Set x ->
            `Set (List.rev (List.map rev x))

let rec flatten = function
    | `Empty -> `Empty
    | `Single it -> it
    | `Set lst ->
            `Set (List.map flatten lst)


let rec iter f = function
    | `Empty 
    | `Single _ as m -> f m
    | `Set lst as m -> 
            List.iter (iter f) lst;
            f m

let rec leaf_iter f = function
    | `Empty -> ()
    | `Single l -> f l
    | `Set lst ->
            List.iter (leaf_iter f) lst

let choose_random s =
    let chosen = ref None in
    let item = ref 1. in
    leaf_iter (fun i ->
                   if Random.float 1. <= (1. /. !item)
                   then chosen := Some i;
                   item := !item +. 1.) s;
    !chosen

let rec to_array (x : 'a t) =
    match x with
    | `Set lst ->
            let comparator (x : 'a t)= 
                match x with
                | `Single _ -> true
                | `Empty | `Set _ -> false
            in
            if List.for_all comparator lst then
                let res = 
                    Array.of_list (List.map (function `Single x -> x | _ -> 
                        failwith "Error") lst)
                in
                `Single res
            else `Set (List.map to_array lst)
    | `Empty | `Single _ -> failwith "Invalid input"

type 'a searcher = Continue of int | Found of 'a

let rec aux_nth n t =
    match t, n with
    | _, Found _ | `Empty, Continue _ -> n
    | `Single x, Continue 0 -> Found x
    | `Single _, Continue n -> Continue (n - 1)
    | `Set lst, _ ->
            let finder acc x =
                match acc with
                | Continue _ -> aux_nth acc x
                | Found _ -> acc
            in
            List.fold_left finder n lst

let nth n t =
    match aux_nth (Continue n) t with
    | Found x -> x
    | Continue v -> failwith "Sexpr.nth"

let rec filter (f : ('b -> bool)) (m : 'b t) : 'b t = 
    match m with
    | `Single x ->
            if f x then `Single x
            else `Empty
    | `Set lst ->
            let x = List.map (filter f) lst in
            (match List.filter (function `Empty -> false | _ -> true) x with
            | [] -> `Empty
            | [x] -> x
            | x -> `Set x)
    | `Empty -> `Empty

let rec combine x =
    match x with
    | `Empty, `Empty -> `Empty
    | (`Single x), y -> `Single (x, y)
    | (`Set lsta), (`Set lstb) ->
            (try `Set (List.map combine (List.combine lsta lstb))
            with _ -> failwith "Sexpr.combine")
    | _, _ -> failwith "Sexpr.combine"

let rec full_combine = function
    | `Empty, `Empty -> `Empty
    | `Single x, `Single y -> `Single (x, y)
    | `Set x, `Set y ->
            (try `Set (List.map full_combine (List.combine x y))
            with _ -> failwith "Sexpr.full_combine")
    | _, _ -> failwith "Sexpr.full_combine"

let rec map_insexpr (f : ('a -> 'b t list)) (x : 'a t) : 'b t = 
    match x with
    | `Empty -> `Empty
    | `Single x -> 
            (match f x with
            | [] -> `Empty
            | [x] -> x
            | x -> `Set x)
    | `Set x ->
            `Set (List.map (map_insexpr f) x)

let status_msg elapsed_wall adv n =
    let present remain msg = 
        (Printf.sprintf "%.0f" remain) ^ " " ^ msg
    in
    let tpi = elapsed_wall /. (float_of_int adv) in
    let remain = tpi *. (float_of_int (n - adv)) in
    let remain = 
        let mine_remain v factor = 
            let mine = floor (v /. factor) in
            let remain = v -. (mine *. factor) in
            mine, remain
        in
        let rec full_calculation remain =
            if (remain) < 60. then present remain "s"
            else if (remain) < 3600. then (* Show minutes *)
                let mine, remain = mine_remain remain 60. in
                (present mine "m, ") ^ full_calculation remain
            else if (remain) < 3600. *. 24. then (* Show days *)
                let mine, remain = mine_remain remain 3600. in
                (present mine "h, ") ^ full_calculation remain
            else 
                let mine, remain = mine_remain remain (3600. *. 24.) in
                (present mine "d, ") ^ full_calculation remain
        in
        full_calculation remain
    in
    "Estimated finish in " ^ remain

let fold_status str ?(eta=true) fn init sexpr =
    let n = cardinal sexpr in
    if n > 1
    then begin
        let status = Status.create str (Some n) "" in
        let time = Timer.start () in
        let feedback adv =
            let elapsed_wall = Timer.wall time in
            let msg =
                match adv - 1 with
                | 0 -> ""
                | adv ->
                      status_msg elapsed_wall adv n
            in
            if eta
            then Status.full_report ~adv ~msg status
            else Status.full_report ~adv status
        in
        let res =
            fold_feedback
                feedback
                fn init sexpr in
        Status.finished status;
        res
    end else fold_left fn init sexpr
    

let map_status str ?(eta=true) fn sexpr =
    let n = cardinal sexpr in
    if n > 1
    then begin
        let status = Status.create str (Some n) "" in
        let time = Timer.start () in
        let feedback adv =
            let elapsed_wall = Timer.wall time in
            let msg =
                match adv - 1 with
                | 0 -> ""
                | adv ->
                      status_msg elapsed_wall adv n
            in
            if eta
            then Status.full_report ~adv ~msg status
            else Status.full_report ~adv status
        in
        let res =
            map_feedback
                feedback
                fn sexpr in
        Status.finished status;
        res
    end else map fn sexpr

(** [compose_status string ?eta fn n start] actually has nothing to do with
    sexprs, but it's grouped here with other progress-indicating functions.
    [compose_status] calculates [fn]^[n]([start]), that is, [fn] applied to
    [start] [n] times.  However, we also pass [fn] an extra parameter, the
    first one, which is the number of the iteration, starting with 0. *)
let compose_status string ?(eta=true) fn n start =
    if n > 1
    then begin
            let status = Status.create string (Some n) "" in
            let time = Timer.start () in
            let feedback adv =
                let elapsed_wall = Timer.wall time in
                let msg =
                    match adv - 1 with
                    | 0 -> ""
                    | adv ->
                          status_msg elapsed_wall adv n
                in
                if eta
                then Status.full_report ~adv ~msg status
                else Status.full_report ~adv status
            in
            let rec f i start =
                feedback (n - i + 1);
                if i = 0
                then start
                else f (pred i) (fn (n - i) start) in
            let res = f n start in
            Status.finished status;
            res
        end
    else fn 0 start

let union a b =
    `Set [a; b]
