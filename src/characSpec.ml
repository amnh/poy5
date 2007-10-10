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

let () = SadmanOutput.register "CharacSpec" "$Revision: 1644 $"

exception Illegal_Assigned_Probabilities of float
exception Undefined_Variable of string
exception Undefined_Function of string
exception Illegal_Combination
exception Precission_Failure

type c = 
    | Base of string 
    | Seq of string
    | Int of string
    | Prep 
    | Hd 
    | Tl 
    | Pre 
    | Suc

type t = {
    functions: string list;
    variables: (string * string * c * SpecIndex.s) list;
    counter : int All_sets.StringMap.t;
    probabilities: (string * float) list;
    decoders_length : float;
}

type p = Counter | Assigned of (string * float) list

let exists s n = 
    List.exists (fun x -> x = n) s.functions

let find_variable s n =
    try List.find (fun (x, _, _, _) -> x = n) s.variables with
    | Not_found -> raise (Undefined_Variable n)

let builtin = ["prep"; "hd"; "tl"; "pre"; "suc"]

let non_standard t = 
    List.filter (fun x -> not (List.exists (fun y -> y = x) builtin)) t.functions

let get_string = function
    | Base s | Seq s | Int s -> s
    | Prep -> "prep"
    | Hd -> "hd"
    | Tl -> "tl"
    | Pre -> "pre"
    | Suc -> "suc"

let empty =
    let the_counter =
        List.fold_left (fun acc x ->
            All_sets.StringMap.add x 0 acc) 
        All_sets.StringMap.empty builtin
    in
    { 
        functions = []; 
        variables = []; 
        counter = the_counter;
        probabilities = []; 
        decoders_length = 0.0;
    }

let add s n p i =
    try
        let i = SpecIndex.to_list i in
        let parameter_processor = fun (name, nametype) ->
            match SpecIndex.find ~index:i ~name:nametype with
            | SpecIndex.Alph _ as x -> (name, nametype, Base nametype, x)
            | SpecIndex.Int _ as x -> (name, nametype, Int nametype, x)
            | SpecIndex.Word _ as x -> (name, nametype, Seq nametype, x)
        in
        let np = List.map parameter_processor p in
        let tmp = 
            { s with
                functions = n :: s.functions;
                variables = np;
            }
        in
        let tmp =
            List.fold_left (fun acc (name, nametype) ->
                if List.exists (fun x -> x = nametype) acc.functions then 
                    acc
                else
                    let nf = nametype :: acc.functions 
                    and nc = All_sets.StringMap.add nametype 0 acc.counter in 
                    { acc with functions = nf; counter = nc; }) 
            tmp p
        in
        { tmp with counter = All_sets.StringMap.add n 0 tmp.counter }
    with
    | e ->
            prerr_string "CharacSpec.add\n";
            raise e

let rec count s n v =
    try if List.exists (fun (x, _, _, _) -> x = n) s.variables then begin
        let (_, x, _, _) = List.find (fun (x, _, _, _) -> x = n) s.variables in
        count s x v
    end else begin 
            let ct = All_sets.StringMap.find n s.counter in
            { s with counter = All_sets.StringMap.add n (ct + v) s.counter }
    end with
    | Not_found -> 
            raise (Undefined_Function n)

let to_list x = All_sets.StringMap.fold (fun a b c -> (a, b) :: c) x.counter []

let merge x y =
    let rec builder x y acc =
    match x, y with
    | (nx, cx) :: tx, (ny, cy) :: ty when nx = ny ->
            builder tx ty ((nx, cx + cy) :: acc)
    | [], [] -> acc
    | _, _ -> raise Illegal_Combination
    in 
    builder x y []

let prepend lst str = 
    let prepender (x, y) = (str^x, y) in
    List.map prepender lst

let join x y nx ny idx = 
    let lx = to_list x
    and ly = to_list y in
    let is_builtin = fun (f, c) -> List.exists (fun x -> f = x) builtin
    in
    let blx, nblx = List.partition is_builtin lx
    and bly, nbly = List.partition is_builtin ly in
    let bl = merge blx bly 
    and nblx = prepend nblx nx
    and nbly = prepend nbly ny in
    let add_element acc (x, y) = count (add acc x [] idx) x y in
    let res = List.fold_left add_element empty bl in
    let res = List.fold_left add_element res nblx in
    let res = List.fold_left add_element res nbly in
    res

let prob s p =
    let cp = 
        match p with
        | Counter -> []
        | Assigned lst -> 
                List.fold_left (fun res ((x, _) as z) ->
                    z :: (List.filter (fun (y, _) -> x <> y) res)) 
                s.probabilities lst 
    in
    { s with probabilities = cp }

let estimate_prob s =
    (* A function that takes lists a and b returns a - b *)
    let difference = 
        fun a b ->
            List.fold_left (fun z (x, _) -> List.filter (fun y -> x <> y) z) a b
    in
    (* The summation of the current set of assigned probabilities *)
    let sum = 
        List.fold_left (fun z (x, y) -> y +. z) 0.0 s.probabilities 
    in
    if sum < 1.0 then begin
        let r = 1.0 -. sum in
        let all_functions, _ = List.split (to_list s) in
        let not_included = difference all_functions s.probabilities in
        let find_count a = 
            try All_sets.StringMap.find a s.counter with
            | Not_found -> assert (false)
        in
        let counts = List.map find_count not_included in
        let total_count = List.fold_left (fun x y -> x + y) 0 counts in
        let fraction = r /. (float_of_int total_count) in
        let res = List.map (fun x -> (float_of_int x) *. fraction) counts in
        let res = List.combine not_included res in 
        prob s (Assigned res)
    end else if sum = 1.0 then begin 
        (* If the summation is 1.0 then all the counted functions must have
        * also an assigned probability set *)
        let all_functions, _ = List.split (to_list s) in
        match difference all_functions s.probabilities with
        | [] -> s
        | _ -> raise (Illegal_Assigned_Probabilities sum)
    end else raise (Illegal_Assigned_Probabilities sum)

let length s n =
    if List.exists (fun (x, _) -> x = n) s.probabilities then begin
        let _, p = List.find (fun (x, _) -> x = n) s.probabilities in
        -. (log p)
    end else 
        raise (Undefined_Function n)

let add_decoder s v =
    { s with decoders_length = s.decoders_length +. v }

let get_decoder s = s.decoders_length

let k c = 
    let summation name count accu =
        let len = length c name in
        match (nan <> len), count with
        | _, 0 | false, _ -> accu
        | _, count -> (len *. (float_of_int count)) +. accu 
    in
    All_sets.StringMap.fold summation c.counter c.decoders_length

let to_formatter items =
    let mapper (name, spec) =
        `Single (Tags.KolSpecs.char_spec, 
        [(Tags.KolSpecs.char_name, name); (Tags.KolSpecs.prob, string_of_float
        (k spec))], 
        `Structured 
        (`Set 
            (List.map (fun (func, prob) ->
                `Single (Tags.KolSpecs.char_fun, [(Tags.KolSpecs.fun_name,
                func); (Tags.KolSpecs.prob, string_of_float prob)], `Structured
                `Empty)) spec.probabilities)))
    in
    let res = 
        match items with
        | [] -> `Empty
        | items -> `Set (List.map mapper items)
    in
    (Tags.KolSpecs.char_index, [], `Structured res)
