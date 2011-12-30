(* POY 5.0 Alpha. A phylogenetic analysis program using Dynamic Homologies.   *)
(* Copyright (C) 2011  Andrés Varón, Le Sy Vinh, Illya Bomash, Ward Wheeler,  *)
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

(* Debug variables/ combinators *)
let (-->) a b = b a

let debug_mem = true
let debug_aln = false
let debug_cst = false

let use_cost_fn = false (* do not memoize cost function  *)
let use_align   = false (* use full alignment matrix     *)
let failwithf format = Printf.ksprintf (failwith) format

(* a wrapper for printf to allow us to direct output for testing *)
let printf format = Printf.printf format

let to_char_list string =
    let rec to_ i len acc =
        if i = len then List.rev acc
                   else to_ (i+1) len ((String.get string i)::acc)
    in
    to_ 0 (String.length string) []

(* general remove gap function for a sequence *)
let remove_gaps gap seq =
    let remove_gaps seq base =
        if base = gap
            then seq
            else let () = Sequence.prepend seq base in seq
    in
    let res = Sequence.fold_right (remove_gaps)
                                  (Sequence.create (Sequence.length seq))
                                  (seq)
    in
    Sequence.prepend res gap;
    res

let sequence_of_string ?(filter_gap=true) alph seq =
    seq --> to_char_list
        --> List.map (String.make 1)
        --> List.map (fun x -> Alphabet.match_base x alph)
        --> Array.of_list
        --> Sequence.of_array
        --> (fun x ->
                if filter_gap
                    then remove_gaps (Alphabet.get_gap alph) x
                    else x)

let pp_ilst chan lst =
    List.iter (fun x -> Printf.fprintf chan "%d| " x) lst

and print_barray2 a =
    for i = 0 to (Bigarray.Array2.dim1 a)-1 do
        for j = 0 to (Bigarray.Array2.dim2 a)-1 do
            printf "%2.10f\t" a.{i,j};
        done; printf "\n";
    done; printf "\n"; ()

type dyn_model = { static : MlModel.model; alph : Alphabet.a; }

let pp_seq model chan seq = Sequence.print chan seq model.alph

let pp_aseq alph chan seq = Sequence.print chan seq alph

let cost_fn m = m.static.MlModel.spec.MlModel.cost_fn

let make_model alph model = { static = model; alph = alph; }

let spec_model alph spec = make_model alph (MlModel.create alph spec)

open Numerical.FPInfix (* fuzzy comparison functions: =., <., >. *)

(* minimum of three with annotations; this will work for all methods *)
let min3_ a at b bt c ct =
    let min_v,min_t =
        if a =. b then
            min a b, [at;bt]
        else if a <  b then
            a, [at] 
        else
            b, [bt]
    in
    if min_v =. c then
        min min_v c, ct::min_t
    else if c < min_v then
        c, [ct]
    else 
        min_v, min_t
 

(* Sequence alignment module for two or three sequences. *)
module type A = sig

    type floatmem
    type s

    (* memory functions *)
    val get_mem     : s -> s -> floatmem
    val create_mem  : int -> int -> floatmem
    val clear_mem   : floatmem -> unit
    (* converting our datatypes to an external type *)
    val s_of_seq    : Sequence.s -> s
    val seq_of_s    : s -> Sequence.s
    val compare     : s -> s -> bool
    (* functions for testing externally *)
    val print_mem   : floatmem -> unit
    val print_s     : s -> Alphabet.a -> unit
    val print_raw   : s -> unit
    val print_cm    : dyn_model -> float -> unit

    (* cost matrix/function; maps states to their cost and possible assignments *)
    val get_cf : dyn_model -> float -> float -> (int -> int -> float * int)
    val get_cm : dyn_model -> float -> float -> (float * int) array array
    val cost   : dyn_model -> float -> float -> (int -> int -> float * int)
    val aln_cost_2      : s -> s -> dyn_model -> float -> float

    (* 2d operations *)
    val cost_2          : ?deltaw:int -> s -> s -> dyn_model -> float -> float -> floatmem -> float
    val full_cost_2     : float -> s -> s -> dyn_model -> float -> float -> floatmem -> float
    val create_edited_2 : s -> s -> dyn_model -> float -> float -> floatmem -> s * s
    val align_2         : ?first_gap:bool -> s -> s -> dyn_model -> float -> float -> floatmem -> s * s * float
    val clip_align_2    : ?first_gap:bool -> Sequence.Clip.s -> Sequence.Clip.s -> dyn_model -> float -> float
                            -> Sequence.Clip.s * Sequence.Clip.s * float * int * Sequence.Clip.s * Sequence.Clip.s
    val median_2        : s -> s -> dyn_model -> float -> float -> floatmem -> s
    val median_2_cost   : s -> s -> dyn_model -> float -> float -> floatmem -> float * s
    val full_median_2   : s -> s -> dyn_model -> float -> float -> floatmem -> s
    val gen_all_2       : s -> s -> dyn_model -> float -> float -> floatmem -> s * s * float * s
    val backtrace       : ?filter_gap:bool -> dyn_model -> floatmem -> s -> s -> float -> float -> s

    (* uppass operations *)
    val closest     : p:s -> m:s -> dyn_model -> float -> floatmem -> s * float
    val get_closest : dyn_model -> float -> i:int -> p:int -> m:int -> int * float

    (* (pseudo) 3d operations *)
    val readjust : s -> s -> s -> dyn_model -> float -> float -> float -> float * s * bool
    val optimize : s -> s -> dyn_model -> float -> floatmem -> float * float

end

IFDEF USE_LIKELIHOOD THEN

module CMPLAlign : A = struct

    type s        = Sequence.s
    type floatmem = (FMatrix.m * Matrix.m)

    (* We define all of the external functions, and compose them properly in
       to the interface below *)

    external s_of_seq : Sequence.s -> s = "%identity"
    external seq_of_s : s -> Sequence.s = "%identity"

    external fm_compose :
        FMatrix.m -> 
        (float,Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array2.t ->
        (float,Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array2.t ->
        ((float,Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array2.t) option ->
        float -> float -> int ->
            (int32,Bigarray.int32_elt, Bigarray.c_layout) Bigarray.Array2.t *
            (float,Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array2.t 
        = "fm_CAML_compose_wrapper" "fm_CAML_compose"

    external full_align :
        FMatrix.m -> Matrix.m -> 
        (float,Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array2.t ->
        (float,Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array2.t ->
        ((float,Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array2.t) option ->
        float -> float -> s -> s -> int ->
            float
        = "falign_CAML_align_2_wrapper" "falign_CAML_align_2"

    external nukk_align :
        FMatrix.m -> Matrix.m -> 
        (float,Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array2.t ->
        (float,Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array2.t ->
        ((float,Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array2.t) option ->
        float -> float -> s -> s -> int ->
            float
        = "falign_CAML_nukk_2_wrapper" "falign_CAML_nukk_2"

    external full_backtrace :
        FMatrix.m -> Matrix.m -> 
        (float,Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array2.t ->
        (float,Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array2.t ->
        ((float,Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array2.t) option ->
        s -> s -> s -> s -> s -> int ->
            unit 
        = "falign_CAML_backtrace_wrapper" "falign_CAML_backtrace"

    external cost_2 :
        FMatrix.m ->
        (float,Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array2.t ->
        (float,Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array2.t ->
        ((float,Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array2.t) option ->
        float -> s -> s ->
            float
        = "falign_CAML_cost_2_wrapper" "falign_CAML_cost_2"

    external closest :
        FMatrix.m ->
        (float,Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array2.t ->
        (float,Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array2.t ->
        ((float,Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array2.t) option ->
        float -> s -> s -> s -> int ->
            float 
        = "falign_CAML_closest_wrapper" "falign_CAML_closest"

    external median_backtrace :
        FMatrix.m -> Matrix.m -> 
        (float,Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array2.t ->
        (float,Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array2.t ->
        ((float,Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array2.t) option ->
        s -> s -> s -> int ->
            unit 
        = "falign_CAML_median_wrapper" "falign_CAML_median"

    external get_closest :
        FMatrix.m ->
        (float,Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array2.t ->
        (float,Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array2.t ->
        ((float,Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array2.t) option ->
        int -> int -> float -> int * float
        = "fm_CAML_get_closest_wrapper" "fm_CAML_get_closest"

    external print_alignment_matrix : 
        FMatrix.m -> Matrix.m -> int -> int -> unit =
            "falign_CAML_print_alignments"


    let create_s ch1 ch2 =
        let len = 
            (Sequence.length (seq_of_s ch1)) + (Sequence.length (seq_of_s ch2))
        in
        s_of_seq (Sequence.create len)


    let get_closest model t ~i ~p ~m =
        get_closest FMatrix.scratch_space model.static.MlModel.u
                    model.static.MlModel.d model.static.MlModel.ui p m t


    let full_backtrace ?(filter_gap=true) a b m ta tb (fmat,mat) =
        let a,ta, b,tb, switch =
            if Sequence.length a > Sequence.length b
                then a,ta, b,tb, false
                else b,tb, a,ta, true
        in
        let gap = Alphabet.get_gap m.static.MlModel.alph in
        let med= create_s a b and ea= create_s a b and eb= create_s a b in
        let () = full_backtrace fmat mat m.static.MlModel.u m.static.MlModel.d
                                m.static.MlModel.ui a b ea eb med gap in
        let ea,eb,med = 
            if filter_gap
                then ea, eb, remove_gaps (Alphabet.get_gap m.alph) med
                else ea, eb, med
        in
        if switch then eb,ea,med else ea,eb,med


    let median_backtrace ?(filter_gap=true) a b m ta tb (fmat,mat) =
        let a,ta, b,tb =
            if Sequence.length a > Sequence.length b
                then a,ta, b,tb
                else b,tb, a,ta
        in
        let med= create_s a b in
        let gap = Alphabet.get_gap m.static.MlModel.alph in
        let () = median_backtrace fmat mat m.static.MlModel.u m.static.MlModel.d
                                  m.static.MlModel.ui a b med gap in
        if filter_gap
            then remove_gaps (Alphabet.get_gap m.alph) med
            else med


    let nukk_align fmat mat u d ui ta tb a b g =
        let a,ta,b,tb =
            if Sequence.length a > Sequence.length b 
                then a,ta, b,tb
                else b,tb, a,ta
        in
        if use_align then
            full_align fmat mat u d ui ta tb a b g
        else
            nukk_align fmat mat u d ui ta tb a b g

    let full_align fmat mat u d ui ta tb a b g =
        let a,ta,b,tb =
            if Sequence.length a > Sequence.length b 
                then a,ta, b,tb
                else b,tb, a,ta
        in
        full_align fmat mat u d ui ta tb a b g

    (* --- implement the A module *)
    let compare a b      = 0 = Sequence.compare (seq_of_s a) (seq_of_s b)

    let print_mem (fm,m) = FMatrix.print fm; Matrix.print_2d m 0 0

    let create_mem _ _   = FMatrix.scratch_space, Matrix.default

    let get_mem _ _      = FMatrix.scratch_space, Matrix.default

    let clear_mem (fm,m) = FMatrix.freeall fm; Matrix.clear_direction m

    let print_cm m t     = 
        let mat = MlModel.compose m.static t in
        for i = 0 to (Bigarray.Array2.dim1 mat) - 1 do
            for j = 0 to (Bigarray.Array2.dim2 mat) - 1 do
                printf " [%f] " (~-. (log mat.{i,j}));
            done;
            printf "\n%!";
        done;
        ()

    let print_s (seq:s) (a) =
        for i = 0 to (Sequence.length seq)-1 do
            try printf "%s" (Alphabet.match_code (Sequence.get seq i) a);
            with (Alphabet.Illegal_Code _ ) -> printf "%d" (Sequence.get seq i);
        done;
        printf "\n%!";
        ()

    let print_raw seq =
        for i = 0 to (Sequence.length seq)-1 do
            printf "|%d" (Sequence.get seq i);
        done;
        printf "|%!";
        ()
    
    let get_cm model t1 t2 =
        let assgn,csts = 
            fm_compose FMatrix.scratch_space model.static.MlModel.u 
                       model.static.MlModel.d model.static.MlModel.ui t1 t2 1
        in
        let a = Bigarray.Array2.dim1 csts and b = Bigarray.Array2.dim2 csts in
        let r = Array.make_matrix a b (0.0,0) in
        for i = 0 to a-1 do for j = 0 to b-1 do
            r.(i).(j) <- (csts.{i,j}, Int32.to_int assgn.{i,j});
        done; done; 
        r

    let get_cf model t1 t2 =
        let assgn,csts =
            fm_compose FMatrix.scratch_space model.static.MlModel.u 
                       model.static.MlModel.d model.static.MlModel.ui t1 t2 1
        in
        (fun i j -> csts.{i,j}, Int32.to_int assgn.{i,j})

    let aln_cost_2 x y model t =
        let fmem,_ = get_mem x y in
        cost_2 fmem model.static.MlModel.u model.static.MlModel.d
               model.static.MlModel.ui t x y

    let cost_2 ?deltaw x y m tx ty (fmat,mat) = 
        let gap = Alphabet.get_gap m.static.MlModel.alph in
        nukk_align fmat mat m.static.MlModel.u m.static.MlModel.d
                        m.static.MlModel.ui tx ty x y gap

    let median_2_cost x y m tx ty ((fmat,mat) as mem) =
        let cost = cost_2 x y m tx ty mem in
        let med  = median_backtrace x y m tx ty mem in
        cost,med

    let median_2 x y m tx ty mem = 
        median_2_cost x y m tx ty mem --> snd

    let cost = get_cf

    let full_cost_2 _ x y m tx ty (fmat,mat) =
        let gap = Alphabet.get_gap m.static.MlModel.alph in
        nukk_align fmat mat m.static.MlModel.u m.static.MlModel.d
                        m.static.MlModel.ui tx ty x y gap

    let create_edited_2 x y m tx ty mem : s * s =
        let () = ignore (cost_2 x y m tx ty mem) in
        let ex,ey,_ = full_backtrace x y m tx ty mem in
        ex,ey

    let clip_align_2 ?first_gap _ _ _ _ _ = failwith "not implemented"

    let backtrace ?(filter_gap=true) m mem x y tx ty =
        median_backtrace ~filter_gap x y m tx ty mem

    let align_2 ?first_gap x y m tx ty ((fmat,mat) as mem) =
        let cost= cost_2 x y m tx ty mem in
        let ex,ey,_ = full_backtrace x y m tx ty mem in
        ex,ey,cost

    let optimize s1 s2 model t mem =
        let update t = (),cost_2 s1 s2 model t 0.0 mem in
        let (t,((),c)) = Numerical.brents_method (t,update t) update in
        (t,c)

    let full_median_2 x y m tx ty mem =
        median_2_cost x y m tx ty mem --> snd

    let gen_all_2 x y m tx ty ((fmat,mat) as mem) =
        let cost= cost_2 x y m tx ty mem in
        let ex,ey,med = full_backtrace x y m tx ty mem in
        ex,ey,cost,med

    let closest ~p ~m model t ((fmat,_) as mem) =
        let gap = Alphabet.get_gap model.alph in
        let (s_new,c) as res = 
            if Sequence.is_empty m gap then begin
                m, 0.0
            end else if 0 = Sequence.compare p m then begin
                let mask = lnot gap in
                let m = Sequence.mapi (fun x p -> if p > 0 then x land mask else x) m in
                let cost = aln_cost_2 p m model t in
                remove_gaps gap m, cost
            end else begin
                let paln, maln, cst = align_2 p m model 0.0 t mem in
                let m = create_s paln maln in
                assert( Sequence.length paln = Sequence.length maln );
                let cost = 
                    closest fmat model.static.MlModel.u  model.static.MlModel.d
                            model.static.MlModel.ui t paln maln m gap in
                remove_gaps gap m, cost
            end
        in
        res

    let readjust s1 s2 s3 model t1 t2 t3 =
        let algn s1 s2 t1 t2 = median_2_cost s1 s2 model t1 t2 (get_mem s1 s2) in
        if debug_aln then begin
            Printf.printf "A : %f - %a\n" t1 (pp_seq model) s1;
            Printf.printf "B : %f - %a\n" t2 (pp_seq model) s2;
            Printf.printf "C : %f - %a\n" t3 (pp_seq model) s3
        end;
        let make_center s1 s2 s3=
            (* first median  *)
            let c12, s12 = algn s1 s2 t1 t2
            and c23, s23 = algn s2 s3 t2 t3
            and c13, s13 = algn s1 s3 t1 t2 in
            (* second median *)
            let c123, s123 = algn s12 s3 0.0 t3
            and c231, s231 = algn s23 s1 0.0 t1
            and c132, s132 = algn s13 s2 0.0 t2 in
            (* sum costs *)
            let c123 = c123 +. c12
            and c231 = c231 +. c23
            and c132 = c132 +. c13 in
            if debug_aln then
                printf "Cost123: %f\tCost231: %f\tCost132: %f\n" c123 c231 c132;
            (* determine best... *)
            if c123 <= c231 then
                if c123 <= c132 then
                    false, c123, closest s3 s12 model t3 (get_mem s3 s12), c123
                else
                    true, c132, closest s2 s13 model t2 (get_mem s2 s13), c123
            else if c231 < c132 then
                true, c231, closest s1 s23 model t1 (get_mem s1 s23), c123
            else
                true, c132, closest s2 s13 model t2 (get_mem s2 s13), c123
        in
        let has_to_print, cst, (s, _), previous = make_center s1 s2 s3 in
        cst, s, has_to_print

end


module FloatAlign : A = struct

    type s        = Sequence.s
    type dir      = Root | Align | Delete | Insert
    type floatmem = ( float * ( int * dir list ) ) array array

    external s_of_seq : Sequence.s -> s = "%identity"
    external seq_of_s : s -> Sequence.s = "%identity"

    let compare a b =
        0 = Sequence.compare (seq_of_s a) (seq_of_s b)

    let choose_dir dirs =
             if List.mem Delete dirs then Delete
        else if List.mem Align  dirs then Align
        else if List.mem Insert dirs then Insert
        else match dirs with | [Root] -> Root | _ -> assert false

    let int_of_dir dirs =
        List.fold_left
            (fun x y -> match y with
                | Delete -> 4 lor x
                | Align  -> 1 lor x
                | Insert -> 2 lor x
                | Root   -> 0 lor x)
            0
            dirs

    let print_mem (mem:floatmem) =
        printf "O:\n";
        for i = 0 to (Array.length mem)-1 do
            for j = 0 to (Array.length mem.(i))-1 do
                let cost,(indel,dirs) = mem.(i).(j) in
                printf "| %d %08.5f[%02d]  " (int_of_dir dirs) (abs_float cost) indel;
            done;
            printf "|\n%!";
        done;
        print_newline ()

    let print_raw seq =
        for i = 0 to (Sequence.length seq)-1 do
            printf "|%d" (Sequence.get seq i);
        done;
        printf "|%!";
        ()

    let print_s (seq:s) (a) =
        for i = 0 to (Sequence.length seq)-1 do
            try printf "%s" (Alphabet.match_code (Sequence.get seq i) a);
            with (Alphabet.Illegal_Code _ ) -> printf "%d" (Sequence.get seq i);
        done;
        printf "\n%!";
        ()

    (* memory processing functions *)
    let mem = ref None

    let create_mem a b =
        let a = max a b in
        let res = Array.make_matrix a a (~-.1.0,(0,[Root])) in
        mem := Some res;
        res

    and clear_mem mem =
        Array.iteri
            (fun i v ->
                Array.iteri
                    (fun j _ -> mem.(i).(j) <- (~-.1.0,(0,[Root])))
                    v)
            mem

    let get_mem a b = match !mem with
        | None -> create_mem ((Sequence.length a)) ((Sequence.length b))
        | Some x ->
            if (Sequence.length a) > (Array.length x) then
                create_mem ((Sequence.length a)) ((Sequence.length b))
            else if (Sequence.length b) > (Array.length x.(0)) then
                create_mem ((Sequence.length a)) ((Sequence.length b))
            else
                x

    let create_cost_matrix model t =
        let mat = MlModel.compose model.static t in
        for i = 0 to (Bigarray.Array2.dim1 mat) - 1 do
            for j = 0 to (Bigarray.Array2.dim2 mat) - 1 do
                mat.{i,j} <- ~-. (log mat.{i,j})
            done;
        done;
        mat

    let print_cm m t = print_barray2 (create_cost_matrix m t)

    (* [create_align_cost_fn m t] Compose the model [m] into a cost matrix with
     * branch length [t] *)
    let create_align_cost_fn m tx ty =
        let cm = create_cost_matrix m (tx+.ty) in
        fun x_i y_i ->
            let cst =
                List.fold_left
                    (fun acc x ->
                        List.fold_left (fun acc y -> min acc cm.{x,y})
                                       (acc)
                                       (BitSet.Int.list_of_packed y_i))
                    (infinity) 
                    (BitSet.Int.list_of_packed x_i)
            in
            if debug_cst then
                printf "Cost: %d(%a) ->%f<- %d(%a)\n%!"
                    x_i pp_ilst (BitSet.Int.list_of_packed x_i) cst
                    y_i pp_ilst (BitSet.Int.list_of_packed y_i);
            cst

    let get_cm model tx ty =
        let f = create_align_cost_fn model tx ty in
        let s = Alphabet.distinct_size model.alph in
        let m = Array.make_matrix s s (~-.1.0,0) in
        for i = 0 to (Array.length m)-1 do
            for j = i to (Array.length m.(i))-1 do
                let cost = f (i+1) (j+1) in
                m.(i).(j) <- cost, (i+1) lor (j+1);
                m.(j).(i) <- cost, (i+1) lor (j+1);
            done;
        done;
        m

    let get_cf m t1 t2 =
        let f = create_align_cost_fn m t1 t2 in
        (fun a b -> f a b, a lor b)

    (* p is single; m is not; find in m the least cost to p *)
    let get_closest model t =
        let cost_matrix = create_cost_matrix model t
        and gap = Alphabet.get_gap model.alph in
        (fun ~i ~p ~m ->
            let p = if p = 0 then gap else p in
            let m = if m = 0 then gap else m in
            (* now choose the best from parent; since we MUST choose one, we do
               not need to worry about floating point issues. *)
            let state,cst =
                List.fold_left
                    (fun acc m ->
                        List.fold_left
                            (fun ((assgn,cst) as acc) p ->
                                let ncst = cost_matrix.{p,m} in
                                if ncst < cst then (m,ncst) else acc)
                            (acc)
                            (BitSet.Int.list_of_packed p))
                    ((~-1),infinity)
                    (BitSet.Int.list_of_packed m)
            in
            let res = 1 lsl state in
(*            Printf.printf "O:%d\t P:%d M:%d --> %d : %f\n" i p m res cst;*)
(*            if debug_cst || ~-1 = state then*)
(*                printf "%d -- p:%02d(%a) m:%02d(%a)\t-(%f/%f)->%02d(%02d)\n%!"*)
(*                              i p pp_ilst (BitSet.Int.list_of_packed p) m*)
(*                              pp_ilst (BitSet.Int.list_of_packed m) cst t state res;*)
            assert( state <> ~-1 );
            res,cst)

    let cost m tx ty =
        if use_cost_fn then
            get_cf m tx ty
        else
            let cm = get_cm m tx ty in
            (fun a b -> cm.(a-1).(b-1))

    (* [alignments mem x y m] builds the alignment of [x] and [y] from the
     * alignment cost matrix, [mem]; the only purpose for the model [m] is to
     * obtain the alphabet, and gap character. *)
    let alignments mem x y m =
        let get_direction i j = mem.(i).(j) --> snd --> snd --> choose_dir
        and gap = Alphabet.get_gap m.alph in
        let rec build_alignments one two i j = match get_direction i j with
            | Align  -> build_alignments ((Sequence.get x i)::one) ((Sequence.get y j)::two) (i-1) (j-1)
            | Delete -> build_alignments ((Sequence.get x i)::one) (gap::two) (i-1) (j)
            | Insert -> build_alignments (gap::one) ((Sequence.get y j)::two) (i) (j-1)
            | Root   -> Sequence.of_array (Array.of_list (gap::one)),
                        Sequence.of_array (Array.of_list (gap::two))
        in
        build_alignments [] [] ((Sequence.length x)-1) ((Sequence.length y)-1)

    (* [backtrace men x y] find the median of x and y *)
    let backtrace ?(filter_gap=true) model mem x y _ _ =
        let get x i = Sequence.get x i
        and union x y i j = (Sequence.get x i) lor (Sequence.get y j) in
        let get_direction i j = mem.(i).(j) --> snd --> snd --> choose_dir in
        let rec build_median acc i j = match get_direction i j with
            | Align  -> build_median ((union x y i j)::acc) (i-1) (j-1)
            | Delete -> build_median ((get x i)::acc) (i-1) j
            | Insert -> build_median ((get y j)::acc) i (j-1)
            | Root   -> Sequence.of_array (Array.of_list ((get x 0)::acc))
        in
        let m = build_median [] ((Sequence.length x)-1) ((Sequence.length y)-1) in
        if filter_gap
            then remove_gaps (Alphabet.get_gap model.alph) m
            else m

    (* [align_2 mem x y m t] Align the sequence [x] and [y] with the cost
     * matrix from [m] and branch length [t], and completely fills the matrix of
     * [mem]. Return minimum cost; the alignment can be recovered from [mem] *)
    let align_2 (mem:floatmem) (x:s) (y:s) (m:dyn_model) (tx:float) (ty:float) =
        let xlen = Sequence.length x and ylen = Sequence.length y in
        Numerical.FPInfix.set_ops (2 * (xlen + ylen));
        let cost = let cf = cost m tx ty in (fun a b -> fst (cf a b) ) in
        let gap = Alphabet.get_gap m.alph in
        let get a b = Sequence.get a b in
        let get_cost i j =
            if i = 0 && j = 0 then begin
                (0.0,[Root])
            end else if i = 0 then begin
                (fst mem.(i).(j-1) +. (cost gap (get y j))),[Insert]
            end else if j = 0 then begin
                (fst (mem.(i-1).(j))) +. (cost (get x i) gap),[Delete]
            end else begin
                min3_ ((fst mem.(i-1).(j))   +. (cost (get x i) gap)) Delete
                      ((fst mem.(i).(j-1))   +. (cost gap (get y j))) Insert
                      ((fst mem.(i-1).(j-1)) +. (cost (get x i) (get y j))) Align
            end
        in
        for i = 0 to xlen - 1 do
            for j = 0 to ylen - 1 do
                let c,ds = get_cost i j in
                mem.(i).(j) <- c,(0,ds)
            done;
        done;
        Numerical.FPInfix.reset ();
        if debug_aln then begin
            let ex,ey = alignments mem x y m in
            printf "C1: "; print_s ex m.alph;
            printf "C2: "; print_s ey m.alph;
            printf "%f ME: " (fst mem.(xlen-1).(ylen-1));
            print_s (backtrace ~filter_gap:false m mem x y tx ty) m.alph;
        end;
        fst mem.(xlen-1).(ylen-1)

    (* [ukkonen_align_2 uk_min uk_max mem x y m t] Align the sequences [x] and
     * [y] using the a cost matrix from the mlmodel [m] and branch length [t],
     * and fills in the matrix [mem]. Returns the cost; alignments can be
     * recovered from the matrix therafter. uk_min and uk_max define the
     * appropriate bounds for using the Ukkonen barrier. *)
    let ukkonen_align_2 ?(uk_min=0) ?(uk_max=max_int) (mem:floatmem) x y m tx ty =
        (* these assertions should be taken care of by other functions. This
         * avoids extra logic within this function, and possible errors *)
        let lenX = Sequence.length x and lenY = Sequence.length y in
        assert( lenX <= lenY );
        Numerical.FPInfix.set_ops (2 * (lenX + lenY));
        if (Array.length mem) < lenX then
            failwithf "FloatSequence.FloatAlign; not enough memory in scratch space; Has %d, requires %d"
                        (Array.length mem) lenX;
        if (Array.length mem.(0)) < lenY then
            failwithf "FloatSequence.FloatAlign; not enough memory in scratch space; Has %d, requires %d"
                        (Array.length mem.(0)) lenY;

        (* obtain the gap representation for this model *)
        let gap = Alphabet.get_gap m.alph in

        (* Obtain the cost matrix for the alignment *)
        let cost_fn = let cm = cost m tx ty in (fun a b -> fst (cm a b)) in

        (* A general function to calculate the barrier of k; this is the length
         * of the horizontal and vertical bands that build the diagonal strip. *)
        let barrier k = (k - (lenY - lenX)) / 2 in

        (* function to wrap up the cost_fn and align cost. *)
        let get_cost i j addX addY : float * int =
            let icost = cost_fn addX addY in
            let ocost = fst mem.(i).(j)
            and indel = fst (snd mem.(i).(j)) in
            if ocost < 0.0 then begin
                print_mem mem;
                failwithf "Illegal access of (%d,%d) = %f" i j ocost
            end;
            (icost+.ocost),indel
        in

        (* update a cell in the matrix by ALL its neighbors; This should only be
         * used to calculate the cost of a cell when all the neighbors exist. *)
        let update_all i j =
            let aln,at = get_cost (i-1) (j-1) (Sequence.get x i) (Sequence.get y j)
            and del,dt = get_cost (i-1) (j)   (gap) (Sequence.get y j)
            and ins,it = get_cost (i)   (j-1) (Sequence.get x i) (gap) in
            (* modify the indel/edit count *)
            let at = if (Sequence.get x i) = (Sequence.get y j) then at else 1+at
            and it = it+1 and dt = dt+1 in
            (* the minimum cost with minimum indel, adjusted with additional
             * indel if necessary. *)
            let m : float * (int * dir list) =
                List.fold_left
                    (fun ((cmin,(imin,dmin)) as amin) (ccur,(icur,dcur)) ->
                        if cmin < ccur then amin
                        else if cmin > ccur then (ccur,(icur,[dcur]))
                        else begin
                            if imin < icur then amin
                            else if imin > icur then (ccur,(icur,[dcur]))
                            else (ccur,(icur,dcur::dmin))
                        end)
                    (aln,(at,[Align]))
                    [(del,(dt,Delete)); (ins,(it,Insert))]
            in
            mem.(i).(j) <- m

        (* Same as above, but will not look at the node to the right (j-1) *)
        and update_row i j =
            let aln,at = get_cost (i-1) (j-1) (Sequence.get x i) (Sequence.get y j)
            and del,dt = get_cost (i-1) (j) gap (Sequence.get y j) in
            let at = if (Sequence.get x i) = (Sequence.get y j) then at else 1+at
            and dt = dt+1 in
            let (cost,indel) as m =
                if del < aln then del,(dt,[Delete])
                else if aln < del then aln,(at,[Align])
                else begin
                    if dt < at then del,(dt,[Delete])
                    else if dt > at then aln,(at,[Align])
                    else aln,(at,[Delete;Align])
                end
            in
            mem.(i).(j) <- m
  
        (* Same as above, but will not look at the node above (i-1) *)
        and update_col i j =
            let aln,at = get_cost (i-1) (j-1) (Sequence.get x i) (Sequence.get y j)
            and ins,it = get_cost (i) (j-1)   (Sequence.get x i) gap in
            let at = if (Sequence.get x i) = (Sequence.get y j) then at else 1+at
            and it = it+1 in
            let (cost,indel) as m =
                if ins < aln then ins,(it,[Insert])
                else if aln < ins then aln,(at,[Align])
                else begin
                    if it < at then ins,(it,[Insert])
                    else if it > at then aln,(at,[Align])
                    else aln,(at,[Insert;Align])
                end
            in
            mem.(i).(j) <- m
        in

        (* Update the Ukkonen barrier cells; The new area is bounded by two
         * strips along the edges of the diagonal band, perpendicular to each
         * other. We move along the barrier, down, and across until a node does
         * not change it's value and stop. *)
        let rec update_matrix (ok:int) (nk:int): unit =
            let ob = barrier ok and nb = barrier nk in
            (* move down each column and update until nothing changes *)
            (* move across each row and update until nothing changes *)
            let run_row i j_min j_max =
                let rec run_row i j =
                    if j >= lenY then
                        ()
                    else if j = j_max then
                        begin update_col i j end
                    else 
                        begin update_all i j; run_row (i) (j+1) end
                in
                run_row i j_min
            in
            for i = 1 to (lenX-1) do
                (* .___._______.___.
                 * |___|_______|___| ; update new sections : right by column
                 *  new   old   new  ;                     : left by row *)
                let old_j_max = i+ob+(lenY-lenX)
                and new_j_max = i+nb+(lenY-lenX)
                and new_j_min = i - nb in
(*                Printf.printf "O:Updating Strip: (%d,%d) -> (%d,%d)\n" new_j_min i new_j_max i;*)
                begin match i - ob with
                    | old_j_min when old_j_min <= 1 ->
                        run_row i (old_j_max) (new_j_max)
                    | old_j_min when new_j_min < 1 ->
                        run_row i 1 new_j_max
                    | old_j_min ->
                        update_row i new_j_min;
                        run_row i (new_j_min+1) new_j_max
                end
            done
        (* If dolphins are so smart, why do they live in Igloos? *)
        and initial_matrix () =
            mem.(0).(0) <- (0.0,(0,[Root]));
            for j = 1 to (lenY-1) do
                let cost,_ = get_cost 0 (j-1) gap (Sequence.get y j) in
                mem.(0).(j) <- cost,(j,[Insert]);
            done;
            for i = 1 to (lenX-1) do
                let cost,_ = get_cost (i-1) 0 (Sequence.get x i) gap in
                mem.(i).(0) <- cost,(i,[Delete]);
            done;
            build_strip (max uk_min ((lenY - lenX) + 1))

        (* build a single strip/band in matrix *)
        and build_strip k =
            let b = barrier k in
            let p_max = ref 0 in
            for i = 1 to (lenX-1) do
                let j_min = max 1 (i - b - 1)
                and j_max = min (lenY-1) (i+b+(lenY-lenX)) in
                Printf.printf "O:Building Strip: (%d,%d) -> (%d,%d)\n" j_min i j_max i;
                if j_min = 1 
                    then update_all i 1
                    else update_row i j_min;
                for j = j_min+1 to j_max-1 do
                    update_all i j
                done;
                if !p_max = (lenX-1)
                    then update_all i j_max
                    else update_col i j_max;
                p_max := j_max;
            done;
            update k

        (* this is to update k and matrix until ending condition *)
        and update k =
            let mat_k = fst (snd (mem.(lenX-1).(lenY-1))) in
            if (k <= mat_k) && (k < uk_max) then begin
                update_matrix k (k*2);
                update (k*2)
            end else begin
                ()
            end
        in
        if not use_align then begin
            initial_matrix ();
            fst (mem.(lenX-1).(lenY-1))
        end else begin
            align_2 mem x y m tx ty
        end

    let aln_cost_2 p m model t : float =
        let cost_matrix = create_cost_matrix model t in
        let cost_single p m =
            List.fold_left
                (fun acc m ->
                    List.fold_left
                        (fun cst p -> min cst cost_matrix.{p,m})
                        (acc)
                        (BitSet.Int.list_of_packed p))
                (infinity)
                (BitSet.Int.list_of_packed m)
        in
        assert( (Sequence.length p) = (Sequence.length m) );
        Sequence.foldi
            (fun acc pos _ ->
                if pos = 0 then acc (* skip gap opening *)
                else
                    let p_i = Sequence.get p pos
                    and m_i = Sequence.get m pos in
                    acc +. cost_single p_i m_i)
            (0.0)
            (m)

(* Functions that implement the module Align *)
    let cost_2 ?deltaw s1 s2 model t1 t2 mem : float =
        if debug_mem then clear_mem mem;
        let s1,s2,t1,t2 =
            if Sequence.length s1 <= Sequence.length s2
                then s1,s2,t1,t2 else s2,s1,t2,t1
        in
        match deltaw with
        | Some uk_max -> ukkonen_align_2 ~uk_max mem s1 s2 model t1 t2
        | None        -> ukkonen_align_2 mem s1 s2 model t1 t2

    let full_median_2 s1 s2 model t1 t2 mem : s =
        if debug_mem then clear_mem mem;
        let s1,s2,t1,t2 =
            if Sequence.length s1 <= Sequence.length s2
                then s1,s2,t1,t2 else s2,s1,t2,t1
        in
        ignore (align_2 mem s1 s2 model t1 t2);
        backtrace model mem s1 s2 t1 t2

    let full_cost_2 cost1 s1 s2 model t1 t2 mem : float =
        if debug_mem then clear_mem mem;
        let s1,s2,t1,t2 =
            if Sequence.length s1 <= Sequence.length s2
                then s1,s2,t1,t2 else s2,s1,t2,t1
        in
        align_2 mem s1 s2 model t1 t2

    let create_edited_2 s1 s2 model t1 t2 mem : s * s =
        if debug_mem then clear_mem mem;
        let s1,s2,t1,t2 =
            if Sequence.length s1 <= Sequence.length s2
                then s1,s2,t1,t2 else s2,s1,t2,t1
        in
        ignore (ukkonen_align_2 mem s1 s2 model t1 t2);
        alignments mem s1 s2 model

    let align_2 ?first_gap s1 s2 model t1 t2 mem : s*s*float =
        if debug_mem then clear_mem mem;
        if Sequence.length s1 <= Sequence.length s2 then
            let cost = ukkonen_align_2 mem s1 s2 model t1 t2 in
            let se1,se2 = alignments mem s1 s2 model in
            se1,se2,cost
        else
            let cost = ukkonen_align_2 mem s2 s1 model t2 t1 in
            let se2,se1 = alignments mem s2 s1 model in
            se1,se2,cost

    let clip_align_2 ?first_gap cs1 cs2 model t1 t2 =
        let s2 = s_of_seq (Sequence.Clip.extract_s cs2)
        and s1 = s_of_seq (Sequence.Clip.extract_s cs1) in
        let mem = get_mem s1 s2 in
        let s1,s2,cost = align_2 s1 s2 model t1 t2 mem in
        let s1 = `DO (seq_of_s s1) and s2 = `DO (seq_of_s s2) in
        s1,s2,cost,0,s1,s2

    let gen_all_2 s1 s2 model t1 t2 mem =
        if debug_mem then clear_mem mem;
        let gen_all_2 s1 s2 t1 t2 =
            let cost  = ukkonen_align_2 mem s1 s2 model t1 t2 in
            let med   = backtrace model mem s1 s2 t1 t2 in
            let s1,s2 = alignments mem s1 s2 model in
            s1,s2,cost,med
        in
        if Sequence.length s1 <= Sequence.length s2 then
            gen_all_2 s1 s2 t1 t2
        else
            let s2,s1,cost,m = gen_all_2 s2 s1 t2 t1 in
            s1,s2,cost,m

    let algn s1 s2 model t1 t2 mem : float * s =
        if debug_mem then clear_mem mem;
        let s1,s2,t1,t2 =
            if Sequence.length s1 <= Sequence.length s2
                then s1,s2,t1,t2
                else s2,s1,t2,t1
        in
        let c = ukkonen_align_2 mem s1 s2 model t1 t2 in
        c,backtrace model mem s1 s2 t1 t2

    let median_2 s1 s2 model t1 t2 mem : s = snd (algn s1 s2 model t1 t2 mem)

    let median_2_cost s1 s2 model t1 t2 mem : float * s = algn s1 s2 model t1 t2 mem

    let closest ~p ~m model t mem : s * float =
        let gap = Alphabet.get_gap model.alph in
        if debug_aln then begin
            printf "P: ";print_s p model.alph;
            printf "M: ";print_s m model.alph;
        end;
        let mask_gaps seq gap =
            let mask = lnot gap in
            Sequence.mapi (fun x p -> if p > 0 then x land mask else x) seq
        and get_closest par : int -> int -> int =
            let gc = get_closest model t in
            (fun om pos ->
                let p = Sequence.get par pos in
                let m =
                    if om = gap || p = gap then om
                    else if (0 <> p land gap) && (0 <> om land gap) then gap
                    else om land (lnot gap)
                in
                fst (gc ~i:pos ~p ~m))
        in
        let (new_m,cst) as res =
            if Sequence.is_empty m gap then
                m, 0.0
            else if 0 = Sequence.compare p m then
                let masked = mask_gaps m gap in
                Sequence.mapi (get_closest masked) masked --> remove_gaps gap, 0.0
            else
                let paln, maln, cst = align_2 p m model (t/.2.0) (t/.2.0) mem in
                assert( Sequence.length paln = Sequence.length maln );
                Sequence.mapi (get_closest paln) maln --> remove_gaps gap, cst
        in
        if debug_aln then begin
            printf " -%f(%d)-> " cst gap; print_s new_m model.alph; print_newline ();
        end;
        res


    (* pseudo 3d operation to find a median of three *)
    let readjust s1 s2 s3 model t1 t2 t3 =
        let algn s1 s2 t1 t2 = algn s1 s2 model t1 t2 (get_mem s1 s2) in
        let make_center s1 s2 s3=
            (* first median  *)
            let c12, s12 = algn s1 s2 t1 t2
            and c23, s23 = algn s2 s3 t2 t3
            and c13, s13 = algn s1 s3 t1 t2 in
            (* second median *)
            let c123, s123 = algn s12 s3 0.0 t3
            and c231, s231 = algn s23 s1 0.0 t1
            and c132, s132 = algn s13 s2 0.0 t2 in
            (* sum costs *)
            let c123 = c123 +. c12
            and c231 = c231 +. c23
            and c132 = c132 +. c13 in
            if debug_aln then
                printf "Cost123: %f\tCost231: %f\tCost132: %f\n" c123 c231 c132;
            (* determine best... *)
            if c123 <= c231 then
                if c123 <= c132 then
                    false, c123, closest s3 s12 model t3 (get_mem s3 s12), c123
                else
                    true, c132, closest s2 s13 model t2 (get_mem s2 s13), c123
            else if c231 < c132 then
                true, c231, closest s1 s23 model t1 (get_mem s1 s23), c123
            else
                true, c132, closest s2 s13 model t2 (get_mem s2 s13), c123
        in
        let has_to_print, cst, (s, _), previous = make_center s1 s2 s3 in
        cst, s, has_to_print

    let optimize s1 s2 model t mem =
        let update t = (),cost_2 s1 s2 model t 0.0 mem in
        let (t,((),c)) = Numerical.brents_method (t,update t) update in
        (t,c)

end


module MPLAlign : A = struct

    type s   = Sequence.s
    type dir = Root | Align  of int
                    | Delete of int
                    | Insert of int

    type floatmem = ( float * ( int * dir list ) ) array array

    external s_of_seq : Sequence.s -> s = "%identity"
    external seq_of_s : s -> Sequence.s = "%identity"

    let compare a b =
        0 = Sequence.compare (seq_of_s a) (seq_of_s b)

    (** choose a standard direction for consistent alignments *)
    let choose_dir dirs =
        let is_del = function | Delete _ -> true | _ -> false
        and is_aln = function | Align  _ -> true | _ -> false
        and is_ins = function | Insert _ -> true | _ -> false  in
             if List.exists is_del dirs then List.find is_del dirs
        else if List.exists is_aln dirs then List.find is_aln dirs
        else if List.exists is_ins dirs then List.find is_ins dirs
        else match dirs with | [Root] -> Root | _ -> assert false

    let int_of_dir dirs =
        List.fold_left
            (fun x y -> match y with
                | Delete _ -> 4 lor x
                | Align  _ -> 1 lor x
                | Insert _ -> 2 lor x
                | Root     -> 0 lor x)
            0
            dirs

    let get_assignment = function
        | Root     -> 0
        | Align a | Delete a | Insert a -> a

    let print_mem (mem:floatmem) =
        printf "O:\n";
        for i = 0 to (Array.length mem)-1 do
            for j = 0 to (Array.length mem.(i))-1 do
                let cost,(indel,dirs) = mem.(i).(j) in
                printf "| %d %8.4f[%02d] " (int_of_dir dirs) (abs_float cost) indel;
            done;
            printf "|\n%!";
        done;
        printf "\n%!"

    let print_raw seq =
        for i = 0 to (Sequence.length seq)-1 do
            printf "|%d" (Sequence.get seq i);
        done;
        printf "|\n%!";
        ()

    let print_s (seq:s) (a) =
        for i = 0 to (Sequence.length seq)-1 do
            try printf "%s" (Alphabet.match_code (Sequence.get seq i) a);
            with (Alphabet.Illegal_Code _) -> printf "%d" (Sequence.get seq i);
        done;
        printf "\n%!";
        ()

    let print_cm m t =
        let mat = MlModel.compose m.static t in
        for i = 0 to (Bigarray.Array2.dim1 mat) - 1 do
            for j = 0 to (Bigarray.Array2.dim2 mat) - 1 do
                printf " [%f] " (~-. (log mat.{i,j}));
            done;
            printf "\n%!";
        done;
        ()

    (* memory processing functions *)
    let mem = ref None

    let create_mem a b =
        let a = max a b in
        let res = Array.make_matrix a a (~-.1.0,(0,[Root])) in
        mem := Some res;
        res

    and clear_mem mem =
        Array.iteri
            (fun i v ->
                Array.iteri
                    (fun j _ -> mem.(i).(j) <- (~-.1.0,(0,[Root]))) v)
            mem

    let get_mem a b = match !mem with
        | None -> create_mem ((Sequence.length a)) ((Sequence.length b))
        | Some x ->
            if (Sequence.length a) > (Array.length x) then
                create_mem ((Sequence.length a)) ((Sequence.length b))
            else if (Sequence.length b) > (Array.length x.(0)) then
                create_mem ((Sequence.length a)) ((Sequence.length b))
            else
                x

    (* [create_mpl_cost_fn m t1 t2]
     * In this cost function, an analogue to MPL, we assign the node to the
     * minimum cost of transforming each child to that node, across their
     * respective branch lengths. *)
    let create_mpl_cost_fn m t1 t2 =
        let cost1,cost2 =
            let mat1 = MlModel.compose m.static t1
            and mat2 = MlModel.compose m.static t2 in
            assert( (Bigarray.Array2.dim1 mat1) = (Bigarray.Array2.dim1 mat2) );
            assert( (Bigarray.Array2.dim2 mat1) = (Bigarray.Array2.dim2 mat2) );
            for i = 0 to (Bigarray.Array2.dim1 mat1) - 1 do
                for j = 0 to (Bigarray.Array2.dim2 mat1) - 1 do
                    mat1.{i,j} <- ~-. (log mat1.{i,j});
                    mat2.{i,j} <- ~-. (log mat2.{i,j})
                done;
            done;
            mat1,mat2
        in
        (* find the cost of median state [me] from [xe] and [ye] *)
        (* return function to calculate costs of polymorphisms *)
        fun x_i y_i ->
            let med_cost xe ye me =
                try cost1.{xe,me} +. cost2.{ye,me}
                with | _ ->
                    failwithf "Error looking up %d/%d -> %d; size/gap: %d/%d from (%d/%d)"
                               xe ye me (Bigarray.Array2.dim1 cost1)
                               (Alphabet.get_gap m.alph) x_i y_i;
            in
            let xs = BitSet.Int.list_of_packed x_i
            and ys = BitSet.Int.list_of_packed y_i in

            let cst,states =
                List.fold_left
                    (fun acc x ->
                        List.fold_left
                            (fun acc y ->
                                let c_min = ref acc in
                                for i = 0 to (m.static.MlModel.alph_s)-1 do
                                    let c = med_cost x y i in
                                    if (fst !c_min) =. c then
                                        c_min := (min c (fst !c_min), i::(snd !c_min))
                                    else if c < (fst !c_min) then
                                        c_min := (c,[i])
                                done;
                                !c_min)
                            acc
                            ys)
                    (infinity,[])
                    xs
            in
            if debug_cst || (0 = List.length states) then begin
                let packed = BitSet.Int.packed_of_list states in
                printf "\tBest Cost: %d(%a) ->%d(%a)<- %d(%a) = %f\n%!"
                    x_i pp_ilst xs packed pp_ilst states y_i pp_ilst ys cst;
            end;
            match classify_float cst, states with
            | _, [] ->
                failwithf "%d -(%f)->?<-(%f)- %d: has no median assignment" x_i t1 t2 y_i
            | (FP_infinite | FP_nan), _ ->
                failwithf "%d -(%f)->inf<-(%f)- %d: has infinite cost" x_i t1 t2 y_i
            | _, _  -> cst, BitSet.Int.packed_of_list states

    let get_cm model tx ty =
        let f = create_mpl_cost_fn model tx ty in
        let s = Alphabet.distinct_size model.alph in
        let m = Array.make_matrix s s (~-.1.0,0) in
        for i = 0 to (Array.length m)-1 do
            for j = 0 to (Array.length m.(i))-1 do
                let cost_state = f (i+1) (j+1) in
                m.(i).(j) <- cost_state;
            done;
        done;
        m

    let get_cf m t1 t2 = create_mpl_cost_fn m t1 t2

    (* p is single; m is not; find in m the least coast to p *)
    let get_closest model t =
        let cost_matrix =
            let mat = MlModel.compose model.static t in
            for i = 0 to (Bigarray.Array2.dim1 mat) - 1 do
                for j = 0 to (Bigarray.Array2.dim2 mat) - 1 do
                    mat.{i,j} <- ~-. (log mat.{i,j});
                done;
            done;
            mat
        and gap = Alphabet.get_gap model.alph in
        (* find the cost of median state [me] from [xe] and [ye] *)
        (fun ~i ~p ~m ->
            let p = if p = 0 then gap else p in
            let m = if m = 0 then gap else m in
            (* now choose the best from parent; since we MUST choose one, we do
               not need to worry about floating point issues. *)
            let state,cst =
                List.fold_left
                    (fun acc m ->
                        List.fold_left
                            (fun ((assgn,cst) as acc) p ->
                                let ncst = cost_matrix.{p,m} in
                                if ncst < cst then (m,ncst) else acc)
                            (acc)
                            (BitSet.Int.list_of_packed p))
                    ((~-1),infinity)
                    (BitSet.Int.list_of_packed m)
            in
            let res = 1 lsl state in
(*            Printf.printf "O:%d\t P:%d M:%d --> %d : %f\n" i p m res cst;*)
(*            if debug_cst || ~-1 = state then*)
(*                printf "%d -- p:%02d(%a) m:%02d(%a)\t-(%f)->%02d(%02d)\n%!"*)
(*                              i p pp_ilst (BitSet.Int.list_of_packed p) m*)
(*                              pp_ilst (BitSet.Int.list_of_packed m) cst state res;*)
            assert( state <> ~-1 );
            res,cst)


    let cost m tx ty =
        if use_cost_fn then
            let cf = create_mpl_cost_fn m tx ty in
            (fun a b ->
                let temp = cf a b in
(*                Printf.printf "\t\tcost[%d][%d] = %f\n" a b (fst temp);*)
                temp)
        else
            let cm = get_cm m tx ty in
            (fun a b ->
                let temp = cm.(a-1).(b-1) in
(*                Printf.printf "\t\tcost[%d][%d] = %f\n" a b (fst temp);*)
                temp)
               


    (* [alignments mem x y m] builds the alignment of [x] and [y] from the
     * alignment cost matrix, [mem]; the only purpose for the model [m] is to
     * obtain the alphabet, and gap character. *)
    let alignments mem x y m =
        let get_direction i j = mem.(i).(j) --> snd --> snd --> choose_dir
        and gap = Alphabet.get_gap m.alph in
        let rec build_alignments one two i j = match get_direction i j with
            | Align  _ -> build_alignments ((Sequence.get x i)::one) ((Sequence.get y j)::two) (i-1) (j-1)
            | Insert _ -> build_alignments (gap::one) ((Sequence.get y j)::two) (i) (j-1)
            | Delete _ -> build_alignments ((Sequence.get x i)::one) (gap::two) (i-1) (j)
            | Root     -> Sequence.of_array (Array.of_list (gap::one)),
                          Sequence.of_array (Array.of_list (gap::two))
        in
        build_alignments [] [] ((Sequence.length x)-1) ((Sequence.length y)-1)

    (* [backtrace men x y] find the median of x and y *)
    let backtrace ?(filter_gap=true) model mem x y _ _ =
        let get x i = Sequence.get x i in
        let get_direction i j = mem.(i).(j) --> snd --> snd --> choose_dir in
        let rec build_median acc i j = match get_direction i j with
            | Align  s -> build_median (s::acc) (i-1) (j-1)
            | Delete s -> build_median (s::acc) (i-1) j
            | Insert s -> build_median (s::acc) i (j-1)
            | Root     -> Sequence.of_array (Array.of_list ((get x 0)::acc))
        in
        let m = build_median [] ((Sequence.length x)-1) ((Sequence.length y)-1) in
        if filter_gap
            then remove_gaps (Alphabet.get_gap model.alph) m
            else m

    (* [align_2 mem x y m t] Align the sequence [x] and [y] with the cost
     * matrix from [m] and branch length [t], and completely fills the matrix of
     * [mem]. Return minimum cost; the alignment can be recovered from [mem] *)
    let align_2 (mem:floatmem) (x:s) (y:s) (m:dyn_model) (tx:float) (ty:float) =
        let xlen = Sequence.length x and ylen = Sequence.length y in
        Numerical.FPInfix.set_ops (3 * (xlen + ylen));
        let cost = cost m tx ty in
        let gap = Alphabet.get_gap m.alph in
        let get a b = Sequence.get a b in
        let get_cost i j =
            if i = 0 && j = 0 then begin
                (0.0,[Root])
            end else if i = 0 then begin
                let cst,s = cost gap (get y j) in
                (fst mem.(i).(j-1) +. cst),[Insert s]
            end else if j = 0 then begin
                let cst,s = cost (get x i) gap in
                (fst mem.(i-1).(j) +. cst),[Delete s]
            end else begin
                let dcst,sd = cost (get x i) gap
                and icst,si = cost gap (get y j)
                and acst,sa = cost (get x i) (get y j) in
                let dcst = (fst mem.(i-1).(j))   +. dcst
                and icst = (fst mem.(i).(j-1))   +. icst
                and acst = (fst mem.(i-1).(j-1)) +. acst in
                let (temp_min,temp_dirs) as temp =
                    min3_ dcst (Delete sd) icst (Insert si) acst (Align sa) in
(*                Printf.printf "(%d,%d)[00]: A:%f\tI:%f\tD:%f --> %f[%d]\n"*)
(*                        i j acst icst dcst temp_min (int_of_dir temp_dirs);*)
                temp
            end
        in
        for i = 0 to xlen - 1 do
            for j = 0 to ylen - 1 do
                let c,ds = get_cost i j in
                mem.(i).(j) <- c,(0,ds)
            done;
        done;
        Numerical.FPInfix.reset ();
        if debug_aln then begin
            let ex,ey = alignments mem x y m in
            printf "%f C1: " tx; print_s ex m.alph;
            printf "%f C2: " ty; print_s ey m.alph;
            printf "%f ME: " (fst mem.(xlen-1).(ylen-1));
            print_s (backtrace ~filter_gap:false m mem x y tx ty) m.alph;
        end;
        fst mem.(xlen-1).(ylen-1)

    (* [ukkonen_align_2 uk_min uk_max mem x y m t] Align the sequences [x] and
     * [y] using the a cost matrix from the mlmodel [m] and branch length [t],
     * and fills in the matrix [mem]. Returns the cost; alignments can be
     * recovered from the matrix therafter. uk_min and uk_max define the
     * appropriate bounds for using the Ukkonen barrier. *)
    let ukkonen_align_2 ?(uk_min=0) ?(uk_max=max_int) (mem:floatmem) x y m tx ty =
        (* these assertions should be taken care of by other functions. This
         * avoids extra logic within this function, and possible errors *)
        let lenX = Sequence.length x and lenY = Sequence.length y in
        assert( lenX <= lenY );
        Numerical.FPInfix.set_ops (3 * (lenX + lenY));
        if (Array.length mem) < lenX then
            failwithf "FloatSequence.FloatAlign; not enough memory in scratch space; Has %d, requires %d"
                        (Array.length mem) lenX;
        if (Array.length mem.(0)) < lenY then
            failwithf "FloatSequence.FloatAlign; not enough memory in scratch space; Has %d, requires %d"
                        (Array.length mem.(0)) lenY;

        (* obtain the gap representation for this model *)
        let gap = Alphabet.get_gap m.alph in

        (* Obtain the cost matrix for the alignment *)
        let cost_fn = cost m tx ty in

        (* A general function to calculate the barrier of k; this is the length
         * of the horizontal and vertical bands that build the diagonal strip. *)
        let barrier k = (k - (lenY - lenX)) / 2 in

        (* function to wrap up the cost_fn and align cost. *)
        let get_cost i j addX addY : float * int * int =
            let cost,states = cost_fn addX addY in
            let ocost = fst mem.(i).(j)
            and indel = fst (snd mem.(i).(j)) in
            if ocost < 0.0 then begin
                failwithf "Illegal access of (%d,%d) = %f" i j ocost
            end;
            (cost+.ocost),indel,states
        in

        (* update a cell in the matrix by ALL its neighbors; This should only be
         * used to calculate the cost of a cell when all the neighbors exist. *)
        let update_all i j =
(*            Printf.printf "\tO:A:%f\tI:%f\tD:%f\n%!"*)
(*                (fst mem.(i-1).(j-1)) (fst mem.(i).(j-1)) (fst mem.(i-1).(j));*)
            let aln,at,sa = get_cost (i-1) (j-1) (Sequence.get x i) (Sequence.get y j)
            and del,dt,sd = get_cost (i-1) (j)   (Sequence.get x i) gap
            and ins,it,si = get_cost (i) (j-1) gap (Sequence.get y j) in
            (* modify the indel/edit count *)
            let at = if (Sequence.get x i) = (Sequence.get y j) then at else 1+at
            and it = it+1
            and dt = dt+1 in
            (* the minimum cost with minimum indel, adjusted with additional
             * indel if necessary. *)
            let (cost,indel) as m : float * (int * dir list) =
                List.fold_left
                    (fun ((cmin,(imin,dmin)) as amin) (ccur,(icur,dcur)) ->
                        if cmin < ccur then amin
                        else if cmin > ccur then (ccur,(icur,[dcur]))
                        else begin
                            if imin < icur then amin
                            else if imin > icur then (ccur,(icur,[dcur]))
                            else (ccur,(icur,dcur::dmin))
                        end)
                    (aln,(at,[Align sa]))
                    [(del,(dt,Delete sd)); (ins,(it,Insert si))]
            in
(*            Printf.printf "\tO:(%d,%d): A:%f\tI:%f\tD:%f --> %f[%d]\n"*)
(*                          j i aln ins del cost (int_of_dir (snd indel));*)
            mem.(i).(j) <- m

        (* Same as above, but will not look at the node to the right (j-1) *)
        and update_row i j =
            let aln,at,sa = get_cost (i-1) (j-1) (Sequence.get x i) (Sequence.get y j)
            and del,dt,sd = get_cost (i-1) (j)   (Sequence.get x i) gap in
            let at = if (Sequence.get x i) = (Sequence.get y j) then at else 1+at
            and dt = dt+1 in
            let (cost,indel) as m =
                if del < aln then del,(dt,[Delete sd])
                else if aln < del then aln,(at,[Align sa])
                else begin
                    if dt < at then del,(dt,[Delete sd])
                    else if dt > at then aln,(at,[Align sa])
                    else aln,(at,[Delete sd;Align sa])
                end
            in
(*            Printf.printf "\tO:(%d,%d): A:%f\tI:%f\tD:%f --> %f[%d]\n"*)
(*                          j i aln 0.0 del cost (int_of_dir (snd indel));*)
            mem.(i).(j) <- m
  
        (* Same as above, but will not look at the node above (i-1) *)
        and update_col i j =
            let aln,at,sa = get_cost (i-1) (j-1) (Sequence.get x i) (Sequence.get y j)
            and ins,it,si = get_cost (i) (j-1) gap (Sequence.get y j) in
            let at = if (Sequence.get x i) = (Sequence.get y j) then at else 1+at
            and it = it+1 in
            let (cost,indel) as m =
                if ins < aln then ins,(it,[Insert si])
                else if aln < ins then aln,(at,[Align sa])
                else begin
                    if it < at then ins,(it,[Insert si])
                    else if it > at then aln,(at,[Align sa])
                    else aln,(at,[Insert si;Align sa])
                end
            in
(*            Printf.printf "\tO:(%d,%d): A:%f\tI:%f\tD:%f --> %f[%d]\n"*)
(*                          j i aln ins 0.0 cost (int_of_dir (snd indel));*)
            mem.(i).(j) <- m
        in

        (* Update the Ukkonen barrier cells; The new area is bounded by two
         * strips along the edges of the diagonal band, perpendicular to each
         * other. We move along the barrier, down, and across until a node does
         * not change it's value and stop. *)
        let rec update_matrix (ok:int) (nk:int): unit =
            let ob = barrier ok and nb = barrier nk in
            (* move down each column and update until nothing changes *)
            (* move across each row and update until nothing changes *)
            let run_row i j_min j_max =
                let rec run_row i j =
                    if j >= lenY then
                        ()
                    else if j = j_max then
                        begin update_col i j end
                    else 
                        begin update_all i j; run_row (i) (j+1) end
                in
                run_row i j_min
            in
            for i = 1 to (lenX-1) do
                (* .___._______.___.
                 * |___|_______|___| ; update new sections : right by column
                 *  new   old   new  ;                     : left by row *)
                let old_j_max = i+ob+(lenY-lenX)
                and new_j_max = i+nb+(lenY-lenX)
                and new_j_min = i - nb in
(*                Printf.printf "O:Updating Strip: (%d,%d) -> (%d,%d)\n" new_j_min i new_j_max i;*)
                begin match i - ob with
                    | old_j_min when old_j_min <= 1 ->
                        run_row i (old_j_max) (new_j_max)
                    | old_j_min when new_j_min < 1 ->
                        run_row i 1 new_j_max
                    | old_j_min ->
                        update_row i new_j_min;
                        run_row i (new_j_min+1) new_j_max
                end
            done
        (* If dolphins are so smart, why do they live in Igloos? *)
        and initial_matrix () =
            mem.(0).(0) <- (0.0,(0,[Root]));
            for j = 1 to (lenY-1) do
                let cost,_,s = get_cost 0 (j-1) gap (Sequence.get y j) in
                mem.(0).(j) <- cost,(j,[Insert s]);
            done;
            for i = 1 to (lenX-1) do
                let cost,_,s = get_cost (i-1) 0 (Sequence.get x i) gap in
                mem.(i).(0) <- cost,(i,[Delete s]);
            done;
            build_strip (max uk_min ((lenY - lenX) + 1))

        (* build a single strip/band in matrix *)
        and build_strip k =
            let b = barrier k in
            let p_max = ref 0 in
            for i = 1 to (lenX-1) do
                let j_min = max 1 (i - b - 1)
                and j_max = min (lenY-1) (i+b+(lenY-lenX)) in
(*                Printf.printf "O:Building Strip: (%d,%d) -> (%d,%d)\n" j_min i j_max i;*)
                if j_min = 1
                    then update_all i 1
                    else update_row i j_min;
                for j = j_min+1 to j_max-1 do
                    update_all i j
                done;
                if !p_max = (lenY-1)
                    then update_all i j_max
                    else update_col i j_max;
                p_max := j_max;
            done;
(*            Printf.printf "Finished Building Strip k = %d\n" k;*)
            update k

        (* this is to update k and matrix until ending condition *)
        and update k =
            let mat_k = fst (snd (mem.(lenX-1).(lenY-1))) in
            if (k <= mat_k) && (k < uk_max) then begin
(*                Printf.printf "O:Increasing K to %d\n" (k*2);*)
                update_matrix k (k*2);
                update (k*2)
            end else begin
(*                Printf.printf "O:Finished Alignment in k:%d = %f\n"*)
(*                                k (fst mem.(lenX-1).(lenY-1));*)
                ()
            end
        in
        if not use_align then begin
            initial_matrix ();
            fst (mem.(lenX-1).(lenY-1))
        end else begin
            align_2 mem x y m tx ty
        end


    let aln_cost_2 p m model t : float =
        let cost_matrix =
            let mat = MlModel.compose model.static t in
            for i = 0 to (Bigarray.Array2.dim1 mat) - 1 do
                for j = 0 to (Bigarray.Array2.dim2 mat) - 1 do
                    mat.{i,j} <- ~-. (log mat.{i,j});
                done;
            done;
            mat
        in
        let cost_single p m =
            List.fold_left
                (fun acc m ->
                    List.fold_left
                        (fun cst p -> min cst cost_matrix.{p,m})
                        (acc)
                        (BitSet.Int.list_of_packed p))
                (infinity)
                (BitSet.Int.list_of_packed m)
        in
        assert( (Sequence.length p) = (Sequence.length m) );
        Sequence.foldi
            (fun acc pos _ ->
                if pos = 0 then acc (* skip gap opening *)
                else
                    let cst = cost_single (Sequence.get p pos) (Sequence.get m pos) in
                    acc +. cst)
            (0.0)
            (m)

(* Functions that implement the module Align *)
    let cost_2 ?deltaw s1 s2 model t1 t2 mem : float =
        if debug_mem then clear_mem mem;
        let s1,s2,t1,t2 =
            if Sequence.length s1 <= Sequence.length s2
                then s1,s2,t1,t2 else s2,s1,t2,t1
        in
        match deltaw with
        | Some uk_max -> ukkonen_align_2 ~uk_max mem s1 s2 model t1 t2
        | None        -> ukkonen_align_2 mem s1 s2 model t1 t2

    let full_median_2 s1 s2 model t1 t2 mem : s =
        if debug_mem then clear_mem mem;
        let s1,s2,t1,t2 =
            if Sequence.length s1 <= Sequence.length s2
                then s1,s2,t1,t2 else s2,s1,t2,t1
        in
        ignore (align_2 mem s1 s2 model t1 t2);
        backtrace model mem s1 s2 t1 t2

    let full_cost_2 cost1 s1 s2 model t1 t2 mem : float =
        if debug_mem then clear_mem mem;
        let s1,s2,t1,t2 =
            if Sequence.length s1 <= Sequence.length s2
                then s1,s2,t1,t2 else s2,s1,t2,t1
        in
        align_2 mem s1 s2 model t1 t2

    let create_edited_2 s1 s2 model t1 t2 mem : s * s =
        if debug_mem then clear_mem mem;
        let s1,s2,t1,t2 =
            if Sequence.length s1 <= Sequence.length s2
                then s1,s2,t1,t2 else s2,s1,t2,t1
        in
        ignore (ukkonen_align_2 mem s1 s2 model t1 t2);
        alignments mem s1 s2 model

    let align_2 ?first_gap s1 s2 model t1 t2 mem : s*s*float =
        if debug_mem then clear_mem mem;
        if Sequence.length s1 <= Sequence.length s2 then
            let cost = ukkonen_align_2 mem s1 s2 model t1 t2 in
            let se1,se2 = alignments mem s1 s2 model in
            se1,se2,cost
        else
            let cost = ukkonen_align_2 mem s2 s1 model t2 t1 in
            let se2,se1 = alignments mem s2 s1 model in
            se1,se2,cost

    let clip_align_2 ?first_gap cs1 cs2 model t1 t2 =
        let s2 = s_of_seq (Sequence.Clip.extract_s cs2)
        and s1 = s_of_seq (Sequence.Clip.extract_s cs1) in
        let mem = get_mem s1 s2 in
        let s1,s2,cost = align_2 s1 s2 model t1 t2 mem in
        let s1 = `DO (seq_of_s s1) and s2 = `DO (seq_of_s s2) in
        s1,s2,cost,0,s1,s2

    let gen_all_2 s1 s2 model t1 t2 mem =
        if debug_mem then clear_mem mem;
        let gen_all_2 s1 s2 t1 t2 =
            let cost  = ukkonen_align_2 mem s1 s2 model t1 t2 in
            let med   = backtrace model mem s1 s2 t1 t2 in
            let s1,s2 = alignments mem s1 s2 model in
            s1,s2,cost,med
        in
        if Sequence.length s1 <= Sequence.length s2 then
            gen_all_2 s1 s2 t1 t2
        else
            let s2,s1,cost,m = gen_all_2 s2 s1 t2 t1 in
            s1,s2,cost,m

    let algn s1 s2 model t1 t2 mem : float * s =
        if debug_mem then clear_mem mem;
        let s1,s2,t1,t2 =
            if Sequence.length s1 <= Sequence.length s2
                then s1,s2,t1,t2
                else s2,s1,t2,t1
        in
        let c = ukkonen_align_2 mem s1 s2 model t1 t2 in
        let b = backtrace model mem s1 s2 t1 t2 in
        c,b

    let median_2 s1 s2 model t1 t2 mem : s = snd (algn s1 s2 model t1 t2 mem)

    let median_2_cost s1 s2 model t1 t2 mem : float * s = algn s1 s2 model t1 t2 mem

    let closest ~p ~m model t mem : s * float =
        let gap = Alphabet.get_gap model.alph in
        let mask_gaps seq gap =
            let mask = lnot gap in
            Sequence.mapi (fun x p -> if p > 0 then x land mask else x) seq
        and get_closest par : int -> int -> int =
            let gc = get_closest model t in
            (fun om pos ->
                let p = Sequence.get par pos in
                let m =
                    if om = gap || p = gap then om
                    else if (0 <> p land gap) && (0 <> om land gap) then gap
                    else om land (lnot gap)
                in
                fst (gc ~i:pos ~p ~m))
        in
        let (s_new,c) as res =
            if Sequence.is_empty m gap then
                m, 0.0
            else if 0 = Sequence.compare p m then
                let masked = mask_gaps m gap in
                let m = Sequence.mapi (get_closest masked) masked in
                let cost = aln_cost_2 p m model t in
                remove_gaps gap m, cost
            else
                let paln, maln, cst = align_2 p m model 0.0 t mem in
                assert( Sequence.length paln = Sequence.length maln );
                let m = Sequence.mapi (get_closest paln) maln in
                remove_gaps gap m, aln_cost_2 paln m model t
        in
        if debug_aln then begin
            printf "\nP: ";print_s p model.alph;
            printf "\nM: ";print_s m model.alph;
            printf " -%f(%f)-> " c t; print_s s_new model.alph;
            print_newline ();
        end;
        res


    let readjust s1 s2 s3 model t1 t2 t3 =
        let algn s1 s2 t1 t2 = algn s1 s2 model t1 t2 (get_mem s1 s2) in
        let make_center s1 s2 s3 =
            (* first median  *)
            let c12, s12 = algn s1 s2 t1 t2
            and c23, s23 = algn s2 s3 t2 t3
            and c13, s13 = algn s1 s3 t1 t2 in
            (* second median *)
            let c123, s123 = algn s12 s3 0.0 t3
            and c231, s231 = algn s23 s1 0.0 t1
            and c132, s132 = algn s13 s2 0.0 t2 in
            (* sum costs *)
            let c123 = c123 +. c12
            and c231 = c231 +. c23
            and c132 = c132 +. c13 in
            if debug_aln then
                printf "Cost123: %f\tCost231: %f\tCost132: %f\n" c123 c231 c132;
            (* determine best... *)
            if c123 <= c231 then
                if c123 <= c132 then
                    false, c123, closest s3 s12 model t3 (get_mem s3 s12), c123
                else
                    true, c132, closest s2 s13 model t2 (get_mem s2 s13), c123
            else if c231 < c132 then
                true, c231, closest s1 s23 model t1 (get_mem s1 s23), c123
            else
                true, c132, closest s2 s13 model t2 (get_mem s2 s13), c123
        in
        let has_to_print, cst, (s, _), previous = make_center s1 s2 s3 in
        cst, s, has_to_print

    let optimize s1 s2 model t mem =
        let update t = (),cost_2 s1 s2 model t 0.0 mem in
        let (t,((),c)) = Numerical.brents_method (t,update t) update in
        (t,c)
end

module MALAlign : A = struct

    type s   = Sequence.s

    type floatmem = float array array

    external s_of_seq : Sequence.s -> s = "%identity"
    external seq_of_s : s -> Sequence.s = "%identity"

    let compare a b =
        0 = Sequence.compare (seq_of_s a) (seq_of_s b)

    let print_mem (mem:floatmem) =
        for i = 0 to (Array.length mem)-1 do
            for j = 0 to (Array.length mem.(i))-1 do
                printf "| %08.5f " (abs_float mem.(i).(j));
            done;
            printf "|\n%!";
        done;
        print_newline ()

    (* memory processing functions *)
    let mem = ref None

    let create_mem a b =
        let a = max a b in
        let res = Array.make_matrix a a (~-.1.0) in
        mem := Some res;
        res

    and clear_mem mem =
        Array.iteri
            (fun i v ->
                Array.iteri
                    (fun j _ -> mem.(i).(j) <- (~-.1.0)) v)
            mem

    let get_mem a b = match !mem with
        | None -> create_mem ((Sequence.length a)) ((Sequence.length b))
        | Some x ->
            if (Sequence.length a) > (Array.length x) then
                create_mem ((Sequence.length a)) ((Sequence.length b))
            else if (Sequence.length b) > (Array.length x.(0)) then
                create_mem ((Sequence.length a)) ((Sequence.length b))
            else
                x

    let create_cost_matrix model t =
        let mat = MlModel.compose model.static t in
        for i = 0 to (Bigarray.Array2.dim1 mat) - 1 do
            for j = 0 to (Bigarray.Array2.dim2 mat) - 1 do
                mat.{i,j} <- ~-. (log mat.{i,j})
            done;
        done;
        mat

    let print_cm m t = print_barray2 (create_cost_matrix m t)

    (* to determine the cost of MAL static likelihood *)
    let sum3 a b c = min (min a b) c

    (* to determine the cost of total likelihood *)
    let sum3 a b c = ~-. (log ((exp (~-.a)) +. (exp (~-.b)) +. (exp (~-.c))))

    (** [create_align_cost_fn m t] Compose the model [m] into a cost matrix with
        branch length [t]. We assign the cost of aligning two bases as (for dna),
            C(x,y) = \prod_{m=\{ACTG-\} P(x,m|t_x) * P(y,m|t_y) * \pi_m *)
    let create_align_cost_fn m tx ty =
        let cm1 = MlModel.compose m.static tx in
        let cm2 = MlModel.compose m.static ty in
        let num = Alphabet.size m.static.MlModel.alph in
        fun x_i y_i ->
            let fn (cst,n) x y =
                let cst = ref 0.0 in
                for i = 0 to num-1 do
                    cst := !cst +. cm1.{x,i} *. cm2.{y,i} *. m.static.MlModel.pi_0.{i};
                done;
                !cst,n+1
            and xs = BitSet.Int.list_of_packed x_i
            and ys = BitSet.Int.list_of_packed y_i in
            let cst,n =
                List.fold_left
                    (fun acc x ->
                        List.fold_left (fun acc y -> fn acc x y) acc ys)
                    (0.0,0) (xs)
            in
            if debug_cst then
                printf "Cost: %d(%a) ->%f/%f(%d)<- %d(%a)\n%!"
                    x_i pp_ilst (BitSet.Int.list_of_packed x_i) cst (log cst) n
                    y_i pp_ilst (BitSet.Int.list_of_packed y_i);
            ~-. (log cst)

    let memoize_cost_fn m tx ty =
        let f = create_align_cost_fn m tx ty in
        let s = Alphabet.distinct_size m.alph in
        let m = Array.make_matrix s s ~-.1.0 in
        for i = 0 to (Array.length m)-1 do
            for j = i to (Array.length m.(i))-1 do
                let cost =
                    try f (i+1) (j+1)
                    with | _ -> failwithf "Error in calculating %d -- %d" (i+1) (j+1)
                in
                m.(i).(j) <- cost;
                m.(j).(i) <- cost;
            done;
        done;
        m

    let basic_cost m tx ty =
        if use_cost_fn then
            let cf = create_align_cost_fn m tx ty in
            (fun a b -> cf a b)
        else
            let cm = memoize_cost_fn m tx ty in
            (fun a b -> cm.(a-1).(b-1) )

    let dyn_2 (mem:floatmem) (x:s) (y:s) (m:dyn_model) (tx:float) (ty:float) =
        let cost = basic_cost m tx ty in
        let gap     = Alphabet.get_gap m.alph in
        let get a b = Sequence.get a b in
        let get_cost i j =
            if i = 0 && j = 0 then begin
                0.0
            end else if i = 0 then begin
                mem.(i).(j-1) +. (cost gap (get y j))
            end else if j = 0 then begin
                mem.(i-1).(j) +. (cost (get x i) gap)
            end else begin
                sum3 (mem.(i-1).(j)   +. (cost (get x i) (gap    )))
                     (mem.(i).(j-1)   +. (cost (gap    ) (get y j)))
                     (mem.(i-1).(j-1) +. (cost (get x i) (get y j)))
            end
        in
        let xlen = Sequence.length x and ylen = Sequence.length y in
        for i = 0 to xlen - 1 do
            for j = 0 to ylen - 1 do
                mem.(i).(j) <- get_cost i j;
            done;
        done;
        mem.(xlen-1).(ylen-1)

    (* only cost function defined for module *)
    let cost_2 ?deltaw s1 s2 model t1 t2 mem : float =
        if debug_mem then clear_mem mem;
        dyn_2 mem s1 s2 model t1 t2

    let optimize s1 s2 model t mem =
        let update t = (),cost_2 s1 s2 model t 0.0 mem in (* don't worry; pully holds! *)
        let (t,((),c)) = Numerical.brents_method (t,update t) update in
        (t,c)

    (* functions not implemented *)
    let print_s _ _                 = failwith "not implemented"
    let print_raw _                 = failwith "not implemented"
    let aln_cost_2 _ _ _ _          = failwith "not implemented"
    let median_2 _ _ _ _ _ _        = failwith "not implemented"
    let median_2_cost _ _ _ _ _ _   = failwith "not implemented"
    let full_median_2 _ _ _ _ _ _   = failwith "not implemented"
    let gen_all_2 _ _ _ _ _ _       = failwith "not implemented"
    let closest ~p ~m _ _ _         = failwith "not implemented"
    let get_closest _ _ ~i ~p ~m    = failwith "not implemented"
    let readjust _ _ _ _ _ _ _      = failwith "not implemented"
    let get_cm _ _ _                = failwith "not implemented"
    let cost _ _ _                  = failwith "not implemented"
    let get_cf _ _ _                = failwith "not implemented"
    let full_cost_2 _ _ _ _ _ _     = failwith "not implemented"
    let create_edited_2 _ _ _ _ _ _ = failwith "not implemented"
    let clip_align_2 ?first_gap _ _ _ _ _ = failwith "not implemented"
    let align_2 ?first_gap _ _ _ _ _ _ = failwith "not implemented"
    let backtrace ?filter_gap _ _ _ _ _ _ = failwith "not implemented"

end

ELSE

module Empty : A = struct

    type s        = unit
    type floatmem = unit

    let s_of_seq _                  = failwith "not implemented"
    let seq_of_s _                  = failwith "not implemented"
    let compare _ _                 = failwith "not implemented"
    let print_mem _                 = failwith "not implemented"
    let create_mem _ _              = failwith "not implemented"
    let clear_mem _                 = failwith "not implemented"
    let get_mem _ _                 = failwith "not implemented"
    let print_cm _ _                = failwith "not implemented"
    let cost_2 ?deltaw _ _ _ _ _    = failwith "not implemented"
    let optimize _ _ _ _            = failwith "not implemented"
    let print_s _ _                 = failwith "not implemented"
    let print_raw _                 = failwith "not implemented"
    let aln_cost_2 _ _ _ _          = failwith "not implemented"
    let median_2 _ _ _ _ _ _        = failwith "not implemented"
    let median_2_cost _ _ _ _ _ _   = failwith "not implemented"
    let full_median_2 _ _ _ _ _ _   = failwith "not implemented"
    let gen_all_2 _ _ _ _ _ _       = failwith "not implemented"
    let closest ~p ~m _ _ _         = failwith "not implemented"
    let get_closest _ _ ~i ~p ~m    = failwith "not implemented"
    let readjust _ _ _ _ _ _ _      = failwith "not implemented"
    let get_cm _ _ _                = failwith "not implemented"
    let cost _ _ _                  = failwith "not implemented"
    let get_cf _ _ _                = failwith "not implemented"
    let full_cost_2 _ _ _ _ _ _     = failwith "not implemented"
    let create_edited_2 _ _ _ _ _ _ = failwith "not implemented"
    let clip_align_2 ?first_gap _ _ _ _ _ = failwith "not implemented"
    let align_2 ?first_gap _ _ _ _ _ _ = failwith "not implemented"
    let backtrace ?filter_gap _ _ _ _ _ _ = failwith "not implemented"

end

module MPLAlign = Empty
module FloatAlign = Empty
module MALAlign = Empty
module CMPLAlign = Empty

END

let pair_distance model sq1 sq2 = match cost_fn model with
    | `FLK ->
        let s1 = FloatAlign.s_of_seq sq1 and s2 = FloatAlign.s_of_seq sq2 in
        FloatAlign.optimize s1 s2 model 0.1 (FloatAlign.get_mem s1 s2)
    | `MAL ->
        let s1 = MALAlign.s_of_seq sq1 and s2 = MALAlign.s_of_seq sq2 in
        MALAlign.optimize s1 s2 model 0.1 (MALAlign.get_mem s1 s2)
    | `MPL ->
        let s1 = MPLAlign.s_of_seq sq1 and s2 = MPLAlign.s_of_seq sq2 in
        MPLAlign.optimize s1 s2 model 0.1 (MPLAlign.get_mem s1 s2)

(* a simple function to test the scores of each of the methods above *)
let test_all alignments channel seq1 seq2 bl1 bl2 model =
    let seq1,seq2,bl1,bl2 =
        if Sequence.length seq1 <= Sequence.length seq2
            then seq1,seq2,bl1,bl2
            else seq2,seq1,bl2,bl1
    and filter_gap = false in (* for backtrace code *)
    let flk_cost,flk_opt_cost,flk_seq1,flk_seq2,flk_median =
        let faln1     = FloatAlign.s_of_seq seq1
        and faln2     = FloatAlign.s_of_seq seq2 in
        let mem       = FloatAlign.get_mem faln1 faln2 in
        let a1,a2,cst = FloatAlign.align_2 faln1 faln2 model bl1 bl2 mem in
        let med       = FloatAlign.backtrace ~filter_gap model mem faln1 faln2 bl1 bl2 in
        let _,ocst    = FloatAlign.optimize faln1 faln2 model (bl1+.bl2) mem in
        cst, ocst, (FloatAlign.seq_of_s a1), (FloatAlign.seq_of_s a2), (FloatAlign.seq_of_s med)
    and mpl_cost,mpl_opt_cost,mpl_seq1,mpl_seq2,mpl_median =
        let dmpl1     = MPLAlign.s_of_seq seq1
        and dmpl2     = MPLAlign.s_of_seq seq2 in
        let mem       = MPLAlign.get_mem dmpl1 dmpl2 in
        let a1,a2,cst = MPLAlign.align_2 dmpl1 dmpl2 model bl1 bl2 mem in
        let med       = MPLAlign.backtrace ~filter_gap model mem dmpl1 dmpl2 bl1 bl2 in
        MPLAlign.print_mem mem;
        let _,ocst    = MPLAlign.optimize dmpl1 dmpl2 model (bl1+.bl2) mem in
        cst, ocst, (MPLAlign.seq_of_s a1), (MPLAlign.seq_of_s a2), (MPLAlign.seq_of_s med)
    and mal_cost, mal_opt_cost =
        let dmal1    = MALAlign.s_of_seq seq1
        and dmal2    = MALAlign.s_of_seq seq2 in
        let mem      = MALAlign.get_mem dmal1 dmal2 in
        let cst      = MALAlign.cost_2 dmal1 dmal2 model bl1 bl2 mem
        and obl,ocst = MALAlign.optimize dmal1 dmal2 model (bl1+.bl2) mem in
        cst,ocst
    and cmpl_cost,cmpl_opt_cost,cmpl_seq1,cmpl_seq2,cmpl_median =
        let dmpl1     = CMPLAlign.s_of_seq seq1
        and dmpl2     = CMPLAlign.s_of_seq seq2 in
        let mem       = CMPLAlign.get_mem dmpl1 dmpl2 in
        let a1,a2,cst = CMPLAlign.align_2 dmpl1 dmpl2 model bl1 bl2 mem in
        let med       = CMPLAlign.backtrace ~filter_gap model mem dmpl1 dmpl2 bl1 bl2 in
        let _,ocst    = CMPLAlign.optimize dmpl1 dmpl2 model (bl1+.bl2) mem in
        cst, ocst, (CMPLAlign.seq_of_s a1), (CMPLAlign.seq_of_s a2), (CMPLAlign.seq_of_s med)
    in
    if alignments then begin
        Printf.fprintf channel "%a  %a %a\n%!"
                        (pp_seq model) mpl_seq1 (pp_seq model) flk_seq1 (pp_seq model) cmpl_seq1;
        Printf.fprintf channel "%a  %a %a\n%!"
                        (pp_seq model) mpl_seq2 (pp_seq model) flk_seq2 (pp_seq model) cmpl_seq2;
    end;
    Printf.fprintf channel "%f %f %f %f %f %f %f %f  %a  %a %a\n%!"
                   mal_cost     mpl_cost     flk_cost cmpl_cost
                   mal_opt_cost mpl_opt_cost flk_opt_cost cmpl_opt_cost
                   (pp_seq model) mpl_median (pp_seq model) flk_median (pp_seq model) cmpl_median

