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


(** {6 Types *)



(** {6 Information metrics for Model selection *)

val aic : model -> int -> int -> float -> float
(** [aic n k l_max] Calculate the Akaike Information Criterion *)

val aicc : model -> int -> int -> float -> float
(** [aic n k l_max] Calculate the AICc; as n->inf this is equal to aic, but has
    a higher parameter penalty for smaller n *)

val bic : model -> int -> int -> float -> float
(** [bic n k l_max] Calculate the Bayesian Information Criterion *)

val hqic : model -> int -> int -> float -> float
(** [hqic is defined as, HQIC =2 ln(ln(n)) * k - 2 * ln(L_max), where n is the number
    of observations, k are the number of parameters, and L_max is the
    log-likelihood of the model under the data. *)


