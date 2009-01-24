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

let () = SadmanOutput.register "SpecIndex" "$Revision: 1644 $"

open StdLabels

exception Not_Defined of string
exception Illegal_Element of string
exception Name_Used 
exception Illegal_Combination

type s = Alph of AlphSpec.t | Int of IntSpec.t | Word of WordSpec.t
type e = InInt of int | InAlph of string | InWord of int
type t = (string * s) list

let empty () = []

let add ~index:i ~name:n ~spec:s =
    match List.exists (fun (x, _) -> x = n) i with
    | false -> (n, s) :: i
    | true -> raise Name_Used

let remove ~index:i ~name:s =
    List.filter (fun (x, _) -> x <> s) i

let length ~index:i ~spec:s ~elem:e =
    match s, e with
    | Alph s, InAlph x -> AlphSpec.length ~alph:s ~elem:x
    | Int s, InInt x -> IntSpec.length s x
    | Word s, InWord x -> WordSpec.length s x
    | _, _ -> raise Illegal_Combination

let find ~index:i ~name:n =
    try
        let (x, y) = List.find ~f:(fun (x, y) -> x = n) i in
        y
    with
    | _ -> raise (Not_Defined n)

let to_list x = x

let name ~index:s ~item:it =
    try
        let (res, _) = List.find ~f:(fun (_, x) -> x = it) s in
        res
    with
    | e ->
            print_endline "This item wasn't found.";
            raise e

let k c =
    let process accu it =
        let res =
            match it with
            | _, Alph a -> AlphSpec.decoder a
            | _, Int i -> IntSpec.decoder i
            | _, Word w -> 0.0 (* For now we will leave this as zero *)
        in 
        res +. accu
    in
    List.fold_left ~f:process ~init:0.0 c


type ff = Xml.xml Xml.contents 
let to_formatter s : Xml.xml =
    let mapper (name, item) =
        let name = (Xml.KolSpecs.spec_name, `String name) in
        match item with
        | Alph t -> 
                let clas = (Xml.KolSpecs.spec_class, `String Xml.KolSpecs.alphabet) in
                ((`Single (Xml.KolSpecs.set_spec, [clas; name], 
                ((AlphSpec.to_formatter t) :> ff))) : Xml.xml Sexpr.t)
        | Int t ->
                let clas = (Xml.KolSpecs.spec_class, `String Xml.KolSpecs.integers) in
                ((`Single (Xml.KolSpecs.set_spec, [clas; name], 
                ((IntSpec.to_formatter t) :> ff))) : Xml.xml Sexpr.t)
        | Word t ->
                let clas = (Xml.KolSpecs.spec_class, `String Xml.KolSpecs.words) in
                ((`Single (Xml.KolSpecs.set_spec, [clas; name], 
                ((WordSpec.to_formatter t) :> ff))) : Xml.xml Sexpr.t)
    in
    let items : Xml.xml Sexpr.t list = List.map mapper s in
    (Xml.KolSpecs.spec_index, [], (`Set items))
