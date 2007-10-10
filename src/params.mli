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

type t 
exception Empty_Parameters
val get_neighborhood : t -> Methods.neighborhood
val get_optimality_criterion : t -> Methods.optimality_criterion
val get_character_weight_method : t -> Methods.character_weight_method
val get_cost_calculation : t -> Methods.cost_calculation
val set_neigborhood : Methods.neighborhood -> t -> t
val set_optimality_criterion : Methods.optimality_criterion -> t -> t
val set_character_weight_method : 
    Methods.character_weight_method -> t -> t
val set_cost_calculation : Methods.cost_calculation -> t -> t
val to_tuple : t -> 
    Methods.neighborhood * Methods.optimality_criterion * 
    Methods.character_weight_method * Methods.cost_calculation
val create : 
    Methods.neighborhood -> Methods.optimality_criterion ->  
        Methods.character_weight_method -> Methods.cost_calculation -> t
val empty : t
val set_defaults : 
    Methods.neighborhood -> Methods.optimality_criterion ->  
        Methods.character_weight_method -> Methods.cost_calculation -> unit
