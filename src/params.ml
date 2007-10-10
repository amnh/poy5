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

let () = SadmanOutput.register "Params" "$Revision: 1644 $"

exception Empty_Parameters

type p = 
    { 
        n : Methods.neighborhood; 
        oc : Methods.optimality_criterion;
        cwm : Methods.character_weight_method;
        cc : Methods.cost_calculation;
    }

type t = 
    | Empty 
    | Non_Empty of p

let empty = Empty

let create n oc cwm cc =
    Non_Empty {n = n; oc = oc; cwm = cwm; cc = cc}

let get_neighborhood = function
    | Non_Empty p -> p.n
    | Empty -> raise Empty_Parameters

let get_optimality_criterion = function
    | Non_Empty p -> p.oc
    | Empty -> raise Empty_Parameters

let get_character_weight_method = function
    | Non_Empty p -> p.cwm
    | Empty -> raise Empty_Parameters

let get_cost_calculation = function
    | Non_Empty p -> p.cc
    | Empty -> raise Empty_Parameters

let set_neigborhood n = function
    | Non_Empty p -> create n p.oc p.cwm p.cc
    | Empty -> raise Empty_Parameters

let set_optimality_criterion n = function
    | Non_Empty p -> create p.n n p.cwm p.cc
    | Empty -> raise Empty_Parameters

let set_character_weight_method n = function
    | Non_Empty p -> create p.n p.oc n p.cc
    | Empty -> raise Empty_Parameters

let set_cost_calculation n = function
    | Non_Empty p -> create p.n p.oc p.cwm n
    | Empty -> raise Empty_Parameters

let to_tuple = function
    | Non_Empty v -> (v.n, v.oc, v.cwm, v.cc)
    | Empty -> raise Empty_Parameters

let defaults = ref empty

let set_defaults a b c d = 
    defaults := create a b c d

let default = !defaults
