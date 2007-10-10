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

let () = SadmanOutput.register "CharacterScripting" "$Revision: 1914 $"

type characters = [
    | `All
    | `Some of (bool * int list)
    | `Names of (bool * string list)
    | `Random of float
    | `AllDynamic
    | `AllStatic
    | `Missing of (bool * int)
]
type std_cs = 
    | Nonadd8 of NonaddCS8.t           (** A set of non additive characters with
                                           at most 8 states *)
    | Nonadd16 of NonaddCS16.t         (** A set of non additive characters with
                                           at most 16 states *)
    | Nonadd32 of NonaddCS32.t         (** A set of non additive characters with
                                           at most 32 states *)
    | Add of AddCS.t                    (** A set of additive characters *)
    | Sank of SankCS.t                  (** A set of sankoff characters *)
    | Dynamic of DynamicCS.t                    (** A set of dynamics *)
    | Set of std_cs list                    (** A set of characters *)


type taxa = characters

module type S = sig
    type cs
    type n
    (* This type must be consistant with (cs, float) Methods.character_input_output
    * *)
    type character_input_output = [ 
        | `Characters of cs Sexpr.t
        | `Floats of float Sexpr.t
    ]

    val distance : cs -> cs -> float

    val median : cs -> cs -> cs

    val scriptchar_operations :
        n list ->
        Data.d ->
        [ `Distance of (characters * characters) * characters
        | `Median of (characters * characters) * characters ] -> character_input_output

    val filter_char_operations : n list -> Data.d -> (taxa * taxa) ->
        characters -> (cs Sexpr.t * cs Sexpr.t) list

    val extract_character : int -> n -> cs 

    val character_operations : [ `Distance of (cs Sexpr.t * cs Sexpr.t) | `Median of (cs
        Sexpr.t * cs Sexpr.t) ] -> [ `Characters of cs Sexpr.t | `Floats of float Sexpr.t ]
end

module Standard = struct
type cs = std_cs
type n = Node.node_data

type character_input_output = [ 
    | `Characters of cs Sexpr.t
    | `Floats of float Sexpr.t
]

let distance a b =
    (* TODO: Sets *)
    match a, b with
    | Nonadd8 s, Nonadd8 t -> 
          NonaddCS8.distance s t
    | Nonadd16 s, Nonadd16 t ->
          NonaddCS16.distance s t
    | Nonadd32 s, Nonadd32 t ->
          NonaddCS32.distance s t
    | Add s, Add t ->
          AddCS.distance s t
    | Sank s, Sank t ->
          SankCS.distance s t
    | Dynamic s, Dynamic t ->
          let d =  DynamicCS.distance s t in
          d
    | Set s, Set t -> assert false      (* TODO: Set *)
    | Dynamic _, _ | Sank _, _ | Add _, _ 
    | Nonadd32 _, _ | Nonadd16 _, _ | Nonadd8 _, _
    | Set _, _ ->
          failwith "Invalid combination"

let median a b =
    match a, b with
    | Nonadd8 s, Nonadd8 t -> 
            Nonadd8 (NonaddCS8.median None s t)
    | Nonadd16 s, Nonadd16 t ->
            Nonadd16 (NonaddCS16.median None s t)
    | Nonadd32 s, Nonadd32 t ->
            Nonadd32 (NonaddCS32.median None s t)
    | Add s, Add t ->
          Add (AddCS.median None s t)
    | Sank s, Sank t ->
            Sank (SankCS.median None s t)
    | Dynamic s, Dynamic t ->
            Dynamic (DynamicCS.median s t)

    | Set s, Set t -> assert false      (* TODO: Set *)
    | Dynamic _, _ | Sank _, _ | Add _, _ 
    | Nonadd32 _, _ | Nonadd16 _, _ | Nonadd8 _, _
    | Set _, _ ->
            failwith "Invalid combination"

let character_operations = function
    | `Distance (a, b) ->
            let res = 
                Sexpr.all_to_all distance a b
            in
            `Floats res
    | `Median (a, b) -> 
        let res =
            Sexpr.all_to_all median a b
        in
        `Characters res

let rec filter_nodes data a nodes =
    match a with
    | `All -> nodes
    | `Some (dont_complement, codes) ->
            let codes = 
                if dont_complement then codes
                else Data.complement_taxa data codes
            in
            List.filter (fun x -> List.mem x.Node.taxon_code codes) nodes
    | `Names (dont_complement, names) ->
            let codes = Data.get_chars_codes data (`Names names) in
            filter_nodes data (`Some (dont_complement, codes)) nodes
    | _ -> failwith "Unsupported CharacterScripting.filter_nodes"

let rec extract_character (code : int) (n : Node.node_data) : cs =
    match (Node.Standard.f_codes [code] n).Node.characters with
    | Node.Nonadd8 x :: _ -> Nonadd8 x.Node.preliminary
    | Node.Nonadd16 x :: _ -> Nonadd16 x.Node.preliminary
    | Node.Nonadd32 x :: _ -> Nonadd32 x.Node.preliminary
    | Node.Add x :: _ -> Add x.Node.preliminary
    | Node.Sank x :: _ -> Sank x.Node.preliminary
    | Node.Dynamic x :: _ -> Dynamic x.Node.preliminary
    | Node.Set x :: _ -> assert false        (* is this right?? TODO *)
    | _ -> failwith "No characters in CharacterScripting.extract_character"

let extract_character_name data name n =
    let code = Data.character_code name data in
    extract_character code n

let rec filter_char_operations nodes data (a, b) c :
    (cs Sexpr.t * cs Sexpr.t) list =
    match c with
    | `Some (dont_complement, codes) ->
            (* [build_of_code code] creates a character set Sexpr.t holding only
            * the characters with code [code]. *)
            let codes =
                if dont_complement then codes
                else 
                    match Data.complement_characters data (`Some codes) with
                    | `Some codes -> codes
                    | _ -> failwith "IMpossible?"
            in
            let build_of_code code =
                let nodesa = filter_nodes data a nodes in
                let nodesb = filter_nodes data b nodes in
                let builder x : cs Sexpr.t =
                    ((`Single (extract_character code x)) : cs Sexpr.t) 
                in
                let charsa = `Set (List.map builder nodesa) 
                and charsb = `Set (List.map builder nodesb) in
                ((charsa, charsb) : (cs Sexpr.t * cs Sexpr.t))
            in
            List.map build_of_code codes
    | `Names (dont_complement, names) ->
            let codes = Data.get_chars_codes data (`Names names) in
            filter_char_operations nodes data (a, b) (`Some (dont_complement,
            codes))
    | `All -> 
            let codes = 
                Hashtbl.fold (fun code _ acc -> code :: acc) 
                data.Data.character_codes []
            in
            filter_char_operations nodes data (a, b) (`Some (true, codes))
    | _ -> failwith "Unsupported CharacterScripting.filter_char_operations"


let scriptchar_operations nodes data meth : character_input_output = 
    match meth with
    | `Distance (tx, code) ->
            let tuples_of_characters =
                filter_char_operations nodes data tx code 
            in
            let distance = fun (charsa, charsb) ->
                Sexpr.all_to_all distance charsa charsb 
            in
            `Floats 
            ((`Set (List.map distance tuples_of_characters) : float Sexpr.t))
    | `Median (tx, code) ->
            let tuples_of_characters = filter_char_operations nodes data tx code in
            let median = fun (charsa, charsb) -> 
                Sexpr.all_to_all median charsa charsb 
            in
            `Characters 
            ((`Set (List.map median tuples_of_characters)) : cs Sexpr.t)
end

