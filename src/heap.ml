(* From Okasaki's book *)

module type S = sig
    
    type contents
    type heap
    exception Empty 
    val empty : heap    
    val is_empty : heap -> bool
    val insert : contents -> heap -> heap
    val findMin : heap -> contents
    val deleteMin : heap -> heap

    val iter : (contents -> unit) -> heap  -> unit
    val fold : ('a -> contents -> 'a) -> 'a -> heap -> 'a
    val elements : heap -> contents list 
end

module Naked (O : Set.OrderedType) = struct

    type contents = O.t
    type tree = Node of (int * O.t * tree list)
    type heap = tree list

    exception Empty

    let empty = []
    let is_empty = function [] -> true | _ -> false

    let rank (Node (r, _, _)) = r
    let root (Node (_, r, _)) = r
    let link ((Node (r, x1, c1)) as t1) ((Node (_, x2, c2)) as t2) =
        if 1 > O.compare x1 x2 then
            Node (r + 1, x1, t2 :: c1)
        else Node (r + 1, x2, t1 :: c2)

    let rec insTree t1 ts =
        match ts with
        | [] -> [t1]
        | t2 :: rest ->
                if rank t1 < rank t2 then t1 :: ts
                else insTree (link t1 t2) rest

    let insert x ts =
        insTree (Node (0, x, [])) ts

    let rec merge ts1 ts2 =
        match ts1, ts2 with
        | [], x
        | x, [] -> x
        | (t1 :: s1), (t2 :: s2) ->
                if rank t1 < rank t2 then t1 :: (merge s1 ts2)
                else if rank t2 < rank t1 then t2 :: (merge ts1 s2)
                else insTree (link t1 t2) (merge s1 s2)

    let rec findMin = function
        | [] -> raise Empty
        | [t] -> root t
        | t :: ts ->
                let x = root t 
                and y = findMin ts in
                if 1 > O.compare x y then x else y

    let deleteMin = function
        | [] -> raise Empty
        | ts ->
                let rec getMin = function
                    | [t] -> (t, [])
                    | t :: ts ->
                            let (t', ts') = getMin ts in
                            if 1 > O.compare (root t) (root t') then
                                (t, ts)
                            else (t', t :: ts')
                    | [] -> assert false
                in
                let ((Node (_, x, ts1)), ts2) = getMin ts in
                merge (List.rev ts1) ts2

    let rec iter f heap =
        List.iter (function (Node (_, x, lst)) ->
            f x;
            iter f lst) heap

    let rec fold f acc heap =
        List.fold_left (fun acc (Node (_, x, lst)) ->
            fold f (f acc x) lst) acc heap

    let elements heap = 
        fold (fun acc x -> x :: acc) [] heap
end


module Make (O : Set.OrderedType) : S with type contents = O.t = Naked (O)
