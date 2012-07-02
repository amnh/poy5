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

let () = SadmanOutput.register "ModelSelection" "$Revision"

let (-->) b a = a b
let failwithf format = Printf.ksprintf failwith format

(** {6 Types *)

(** This is the same type as Data.d.branches *)
type branch_lengths = 
    (string, ((string, float) Hashtbl.t) All_sets.IntSetMap.t) Hashtbl.t

(** Defines the details required to diagnose a tree; report diagnosis; and
    select a model based on their AIC, BIC, LRT information *)
type model_selection =
    {   tree : Ptree.t;
          ml : (float * MlModel.model * branch_lengths) list; }


(** {6 General Information Metrics for Model Selection *)

(** [aic] is defined as, AIC = (2*n / (n-k-1)) * k - 2 (ln L_max), where n is
    the number of observations, k are the number of parameters, and L_max is the
    log-likelihood of the model under the data. *)
let aic k n l_max =
    let k = float_of_int k and n = float_of_int n in
    (2.0 *. n *. k /. (n -. k -. 1.0)) -. 2.0 *. l_max

(** [aicc] is defined as AICc = AIC + (2 k (k -1)) / (n-k-1); this criteria is the
    same as AIC as n -> inf (since n is in the denominator), as is a stricter
    parameter penalized version AIC. *)
let aicc k n l_max =
    let aic = aic m k n l_max and n = float_of_int n and k = float_of_int k in
    let c   = (2.0 *. k *. ( k +. 1.0)) /. (n -. k -. 1.0) in
    aic +. c

(** [bic is defined as, BIC = ln(n) * k - 2 * ln(L_max), where n is the number
    of observations, k are the number of parameters, and L_max is the
    log-likelihood of the model under the data. *)
let bic k n l_max =
    let k = float_of_int k and n = float_of_int n in
    ((log n) *. k) -. 2.0 *. l_max

(** [hqic is defined as, HQIC =2 ln(ln(n)) * k - 2 * ln(L_max), where n is the number
    of observations, k are the number of parameters, and L_max is the
    log-likelihood of the model under the data. *)
let hqic k n l_max =
    let k = float_of_int k and n = float_of_int n in
    (2.0 *. (log (log n)) *. k) -. 2.0 *. l_max



