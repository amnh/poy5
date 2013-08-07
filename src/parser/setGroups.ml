(* POY 5.0. A phylogenetic analysis program using Dynamic Homologies.         *)
(* Copyright (C) 2011  Andrés Varón, Lin Hong, Nicholas Lucaroni, Ward Wheeler*)
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

type set_type =
        | Group 
        | Any of float
        | Tree of float             (* origin/loss cost *)
type 'a ct =                        (* this is like a sexpr *)
        | Elt of 'a
        | Set of 'a * set_type * 'a ct list
type t = string ct

(* What constitutes a taxon name *)
let taxon_name = FileStream.is_not FileStream.is_ws_nl

let set_type_to_string = function
    | Group -> "G:"
    | Any f -> "A(" ^ string_of_float f ^ "):"
    | Tree f -> "T(" ^ string_of_float f ^ "): "
let rec to_string (groups : t) : string =
    match groups with
    | Elt taxon -> taxon
    | Set (name, t, list) ->
          "{"
          ^ set_type_to_string t
          ^ (String.concat " " (List.map to_string list)) ^ "}"

(** [of_file fn] reads a CT definition from a file *)
let of_file fn : t list =
    let r = new FileStream.file_reader fn in
    let line = r#read_line in
    if not (Wildcard.anywhere_match (Str.regexp "COMPLEX") line) then begin
        Status.user_message Status.Error 
        ("Illegal@ Complex@ Terminal@ File:@ The@ file@ doesn't@ start@ " ^
        "with@ a@ line@ containing@ the@ COMPLEX@ keyword.");
        failwith "Illegal Complex file format";
    end else
        let rec read (curr : t list) : t list =
            let read_taxon () =
                r#read_while taxon_name in
            r#skip_ws_nl;
            if r#match_prefix "}"
            then List.rev curr
            else
            try
                let taxon = read_taxon () in
                if taxon = ""
                then failwith "unexpected in group parser 000001";

                (* check whether this is a group *)
                r#skip_ws_nl;
                if r#match_prefix "{"
                then begin
                        (* try matching a specification prefix *)
                        let set_type =
                            if r#match_prefix "G:"
                            then Group
                            else if r#match_prefix "A("
                            then begin
                                    let cost = r#read_float in
                                    if not (r#match_prefix "):")
                                    then failwith "Expected close parenthesis and colon in any cost spec.";
                                    Any cost
                            end else if r#match_prefix "T("
                            then begin
                                    let cost = r#read_float in
                                    if not (r#match_prefix "):")
                                    then failwith "Expected close parenthesis and colon in tree cost spec.";
                                    Tree cost
                                end
                            else Any 0. in
                        let set_elts = read [] in
                        let set = Set (taxon, set_type, set_elts) in
                        read (set :: curr)
                    end
                else read (Elt taxon :: curr)
            with End_of_file -> List.rev curr
        in 
        read []

let unify_names n1 n2 =
    if n1 = n2
    then n1
    else n1 ^ " " ^ n2

(** [unify a b] unifies two set specification types, or raises an exception
    if they cannot be unified.  We do this because an [`Any_Of] set could be matched
    with another type, in which case that second object should be implicitly
    converted to an [`Any_Of] *)
let rec unify a b = match a, b with
| Elt _, Elt _ -> a
| Set (n1, s1, c1), Set (n2, s2, c2) ->
      let n = unify_names n1 n2 in
      begin
          match s1, s2 with
          | Group, Group ->
                if List.length c1 = List.length c2
                then begin
                        let types = List.map2 unify c1 c2 in
                        Set (n, Group, types)
                    end
                else failwith "unify: mismatched cardinality of grouping sets"
          | Any v1, Any v2 ->
                (* We need to unify a type for all of the elements *)
                if v1 <> v2 
                then failwith "unify: trees with different origin/loss costs"
                else
                    let utype = unify_list (c1 @ c2) in
                    Set (n, Any v1, [utype])
          | Tree f1, Tree f2 ->
                (* Verify that the floats are the same *)
                if f1 <> f2
                then failwith "unify: trees with different origin/loss costs"
                else begin
                        let utype = unify_list (c1 @ c2) in
                        Set (n, Tree f1, [utype])
                    end
          (* If Any is set against another type, we insert an implicit Any *)
          | Any v, _ ->
                let utype = unify_list (b :: c1) in
                Set (n, Any v, [utype])
          | _, Any v ->
                let utype = unify_list (a :: c2) in
                Set (n, Any v, [utype])
          | _ ->
                failwith "unify: incompatible set types"
      end
| Set (sn, Any j, sc), Elt eltn
| Elt eltn, Set (sn, Any j, sc) ->
      Set (sn, Any j, [Elt eltn])
| _ -> failwith "unify: cannot unify element with set"
and unify_list = function
    | [] -> failwith "unify_list: empty list"
    | l :: ls -> List.fold_left unify l ls

(** [coerce_to_type utype set] puts a type into canonical form *)
let rec coerce_to_type utype set = match utype, set with
| Elt _, Elt _ -> set
| Set (_, Group, subtypes), Set (n, Group, subelts) ->
      Set (n, Group, List.map2 coerce_to_type subtypes subelts)
| Set (_, Any v, subtype), Set (n, Any _, subelts) ->
      let subtype = match subtype with
      | [s] -> s
      | _ -> assert false in
      Set (n, Any v, List.map (coerce_to_type subtype) subelts)
| Set (_, Tree f, subtype), Set (n, Tree f', subelts) ->
      let subtype = match subtype with
      | [s] -> s
      | _ -> assert false in
      assert(f = f');
      Set (n, Tree f, List.map (coerce_to_type subtype) subelts)
| Set (n, Any v, subtype), elt ->
      let subtype = match subtype with
      | [s] -> s
      | _ -> assert false in
      Set (n, Any v, [coerce_to_type subtype elt])
| _ -> assert false
