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

(** [Metric] currently serves two purposes:  It implements various functions for
    / over norm and metric spaces, and it also interfaces with specific
    characters (so far just Sankoff) to embed them into L1, which can be
    simulated with additive characters. *)

(** {2 Norm and metric signatures} *)

(** Norm space over a field.  In theory, this should have an addition operation
    and a zero, but we use subtraction for turning norm spaces into metric
    spaces, so that's what we do here. *)
module type Norm =
sig
    type elt                            (** Element of a norm space *)
    type space                          (** The type of the space itself *)
    val sub : space -> elt -> elt -> elt (** Binary subtraction operator *)
    val norm : space -> elt -> float    (** Unary norm *)
end

(** Generic metric space.  We only provide a distance function. *)
module type Metric =
sig
    type elt                            (** Element of a norm space *)
    type space                          (** The type of the space itself *)
    val dist : space -> elt -> elt -> float (** Binary distance *)
end

(** A functor to create the metric of a norm space *)
module MetricOfNorm :
    functor (N : Norm) ->
sig
    type elt = N.elt
    type space = N.space
    val dist : N.space -> N.elt -> N.elt -> float
end

(** {2 Specific norms and metrics} *)

(** [DiscreteMetric] implements a simple discrete metric space with specified
    costs between pairs of elements *)
module DiscreteMetric :
sig
    type elt = int                      (** Elements are integer indices *)
    type space = int * float array array
    (** We store the space as [(size, cost_array)], where [size] is the number
        of elements *)

    val dist : space -> elt -> elt -> float
    (** [dist space a b] will return the distance from [a] to [b] *)
    val make_space : int -> float array array -> space
    (** [make_space n_elts dist_array] makes a new discrete metric space *)
    val make_elts : int -> int list
    (** [make_elts i] just makes a list [\[0; 1; ...; i-1;\]] *)
end

(** [LpdNorm] defines a d-dimensional l-family norm space {% ($\ell_p^d$) %}. *)
module LpdNorm :
sig
    type space = float * int            (** Simply store p and d *)
    type elt = float array              (** d-length array of values *)
    val add : space -> elt -> elt -> elt
    (** [add space a b] provides addition *)
    val sub : space -> elt -> elt -> elt
    (** [sub space a b] provides subtraction *)
    val zero : space -> elt
    (** [zero space] provides the zero element *)
    val norm : space -> elt -> float
    (** [norm space elt] calculates the norm of the element *)
    val make_space : float -> int -> space
    (** [make_space p d] creates an lp metric of dimension d *)
    val make_elt : int -> (int -> float) -> elt
    (** Make a d-dimensional element given a float-generating function *)    
end

(** [LpdMetric] is the metric of [LpdNorm] *)
module LpdMetric :
sig
    type elt = LpdNorm.elt
    type space = LpdNorm.space
    val dist : LpdNorm.space -> LpdNorm.elt -> LpdNorm.elt -> float
end

(** The [MetricUtil] functor creates useful utility functions over a metric.
    Currently, it implements minimum distance to a set and a Bourgain embedding
    in L1. *)
module MetricUtil :
    functor (M : Metric) ->
sig
    val min_dist : M.space -> M.elt -> M.elt list -> float
    (** [min_dist space elt elts] gives you the minimum distance from [elt] to a
        member of [elts] *)
    val embed_Bourgain :
        M.space ->
        M.elt list ->
        float -> (unit -> bool) -> LpdMetric.space * LpdMetric.elt list
        (** [embed space elts c bool_gen] embeds a list of elements from a space
            [space] into l{_1}{^d}.  See Bourgain 1985.
            @param c This is a constant argument for the embedding
            @param bool_gen This should be a function to generate random boolean
            values; Random.bool works fine. *)
end

(** This is {!Metric.MetricUtil} instantiated for discrete metric spaces *)
module DiscreteUtil :
sig
    val min_dist :
        DiscreteMetric.space ->
        DiscreteMetric.elt -> DiscreteMetric.elt list -> float
    val embed_Bourgain :
        DiscreteMetric.space ->
        DiscreteMetric.elt list ->
        float -> (unit -> bool) -> LpdMetric.space * LpdMetric.elt list
end


(** {2 Conversions and distortions} *)

(** a function that returns the distance between two elements *)
type 'a distfn = 'a -> 'a -> float

(** a set (list) of elements with their associated distance function *)
type 'a eltset = 'a list * 'a distfn

val distortion :
    'a eltset -> 'b eltset -> float
(** [distortion (as, adist) (bs, bdist)] calculates the distortion between two
    metrics.  The lists [as] and [bs] must contain the same A-embedded and
    B-embedded elements in the same order. *)
val convert_lpd_Bourgain :
    float array array -> float -> LpdMetric.space * LpdMetric.elt list
(** [convert_lpd_Bourgain array c] returns the Lpd embedding of distance matrix
    [array] using parameter [c] *)
val pairwise_dist : 'a eltset -> float array array
(** [pairwise_dist (elts, distfn)] returns the distance matrix for the given
    elements *)
val make_discrete : float array array -> DiscreteMetric.elt eltset
(** make a discrete metric from a float matrix *)
val make_lpd :
    LpdNorm.space * LpdNorm.elt list -> LpdNorm.elt eltset
(** make an {!Metric.eltset} out of an Lpd space and element list *)
val avg_n_trials : ('a -> float) -> int -> 'a -> float
(** [avg_n_trials fn num input] returns the average result of [fn] when called
    with [input] [num] times *)
val verify_triangle : float array array -> unit
(** [verify_triangle] verifies that the triangle inequality holds in a distance
    matrix *)
val make_float_list : float -> float -> float -> float list
(** [make_float_list from to step] makes a list of floats from [from] to [to]
    with increment [step] *)

(** {3 Sankoff conversion} *)

val sankoff_convert_fn :
    ?c:float ->
    int array array -> SankCS.t Node.r -> Node.cs list -> Node.cs list
(** [sankoff_convert_fn ?(c=0.95) tcm node] returns a function [conv c list]
    which will take a Sankoff character, convert it to additive characters, and
    return the passed list with these additive characters prepended. *)

val sankoff_null_convert : 'a -> 'b -> 'c
(** a "conversion function" which raises a failure *)

val remember_sankoff_convert :
    ?c:float ->
    SankCS.t ->
    (SankCS.t Node.r -> Node.cs list -> Node.cs list) ->
    SankCS.t Node.r -> Node.cs list -> Node.cs list
(** This function allows us to chain together a single Sankoff conversion
    function which will work for all Sankoff characters present in an input.
    This function takes a Sankoff character [sc] and an existing Sankoff conversion
    function [conv] (such as {!Metric.sankoff_null_convert} above) and returns a
    function which will convert Sankoffs just like [conv] but will also convert
    characters homologous [sc] (checked using [sc]'s code). *)
