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

(* Debug variables/ combinators *)
let (-->) a b = b a
let debug_mem = false
let debug_aln = false
let use_cost_fn = true
let use_align = true
let failwithf format = Printf.ksprintf (failwith) format

(* a wrapper for printf to allow us to direct output for testing *)
let printf format = Printf.printf format
(*let printf format =*)
(*    Printf.ksprintf (Status.user_message Status.Information) format*)

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

let sequence_of_string seq alph = 
    seq --> to_char_list
        --> List.map (String.make 1)
        --> List.map (fun x -> Alphabet.match_base x alph)
        --> Array.of_list
        --> Sequence.of_array
        --> remove_gaps (Alphabet.get_gap alph)

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

let cost_fn m = m.static.MlModel.spec.MlModel.cost_fn

open Numerical.FPInfix (* fuzzy comparison functions: =., <., >. *)

(* minimum of three with annotations; this will work for all methods *)
let min3_ a at b bt c ct =
    if a >. b then begin
        if b >. c then (c,[ct]) 
        else if b =. c then ((min b c),[bt;ct]) 
        else (b,[bt])
    end else if a =. b then begin
        if b =. c then ((min a (min b c)),[at;bt;ct]) 
        else if a >. c then (c,[ct]) 
        else ((min a b),[at;bt])
    end else begin
        if c >. a then (a,[at]) 
        else if c =. a then ((min a c),[at;ct]) 
        else (c,[ct])
    end

(* Sequence alignment module for two or three sequences. *)
module type A = sig

    type floatmem
    type s

    (* auxiliary/helper functions *)
    val get_mem     : s -> s -> floatmem
    val create_mem  : int -> int -> floatmem
    val clear_mem   : floatmem -> unit
    (* converting our datatypes to an external type *)
    val s_of_seq    : Sequence.s -> s 
    val seq_of_s    : s -> Sequence.s
    (* functions for testing externally *)
    val print_mem   : floatmem -> unit
    val print_s     : s -> Alphabet.a -> unit
    val print_raw   : s -> unit
    val print_cm    : dyn_model -> float -> unit

    (* cost matrix/function; maps states to their cost and possible assignments *)
    val get_cf : dyn_model -> float -> float -> (int -> int -> float * int)
    val get_cm : dyn_model -> float -> float -> (float * int) array array
    val cost   : dyn_model -> float -> float -> (int -> int -> float * int)

    (* 2d operations *)
    val cost_2          : ?deltaw:int -> s -> s -> dyn_model -> float -> float -> floatmem -> float
    val aln_cost_2      : s -> s -> dyn_model -> float -> float
    val full_cost_2     : float -> s -> s -> dyn_model -> float -> float -> floatmem -> float
    val create_edited_2 : s -> s -> dyn_model -> float -> float -> floatmem -> s * s
    val align_2         : ?first_gap:bool -> s -> s -> dyn_model -> float -> float -> floatmem -> s * s * float
    val clip_align_2    : ?first_gap:bool -> Sequence.Clip.s -> Sequence.Clip.s -> dyn_model -> float -> float
                            -> Sequence.Clip.s * Sequence.Clip.s * float * int * Sequence.Clip.s * Sequence.Clip.s
    val median_2        : s -> s -> dyn_model -> float -> float -> floatmem -> s
    val median_2_cost   : s -> s -> dyn_model -> float -> float -> floatmem -> float * s
    val full_median_2   : s -> s -> dyn_model -> float -> float -> floatmem -> s
    val gen_all_2       : s -> s -> dyn_model -> float -> float -> floatmem -> s * s * float * s
    val backtrace       : ?filter_gap:bool -> dyn_model -> floatmem -> s -> s -> s

    (* uppass operations *)
    val closest     : p:s -> m:s -> dyn_model -> float -> floatmem -> s * float
    val get_closest : dyn_model -> float -> i:int -> p:int -> m:int -> int * float

    (* (pseudo) 3d operations *)
    val readjust : s -> s -> s -> dyn_model -> float -> float -> float -> floatmem -> float * s * bool
    val optimize : s -> s -> dyn_model -> float -> floatmem -> float * float

end 


module FloatAlign : A = struct

    type s        = Sequence.s
    type dir      = Root | Align | Delete | Insert
    type floatmem = ( float * ( int * dir list ) ) array array

    external s_of_seq : Sequence.s -> s = "%identity"
    external seq_of_s : s -> Sequence.s = "%identity"

    let s_of_seq = s_of_seq
    let seq_of_s = seq_of_s

    let choose_dir dirs = 
             if List.mem Delete dirs then Delete
        else if List.mem Align  dirs then Align
        else if List.mem Insert dirs then Insert
        else match dirs with | [Root] -> Root | _ -> assert false

    let print_mem (mem:floatmem) = 
        let dirToString = function
            | Root   -> "X" | Align  -> "\\"
            | Delete -> "|" | Insert -> "-"
        in
        for i = 0 to (Array.length mem)-1 do
            for j = 0 to (Array.length mem.(i))-1 do
                let cost,(indel,dirs) = mem.(i).(j) in
                printf "|  %s %02.3f[%d]  " 
                    (dirToString (choose_dir dirs)) (abs_float cost) indel;
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
        let cost_matrix = create_cost_matrix m (tx+.ty) in
        fun x_i y_i ->
            let fn (cst,n) x y =
                try (cst +. cost_matrix.{x,y},n+1)
                with | _ -> failwithf "Cannot find from cost matrix: %d, %d" x y
            and xs = BitSet.Int.list_of_packed x_i 
            and ys = BitSet.Int.list_of_packed y_i in
            let cst,n =
                List.fold_left
                    (fun acc x ->
                        List.fold_left (fun acc y -> fn acc x y) acc ys)
                    (0.0,0) (xs)
            in
            if debug_aln then
                printf "Cost: %d(%a) ->%f/%d<- %d(%a)\n%!" 
                    x_i pp_ilst (BitSet.Int.list_of_packed x_i) cst n
                    y_i pp_ilst (BitSet.Int.list_of_packed y_i);
            cst /. (float_of_int n)

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
            (* now choose the best from parent *)
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
            if debug_aln || ~-1 = state then
                printf "%d -- p:%02d(%a) m:%02d(%a)\t-(%f/%f)->%02d(%02d)\n%!"
                              i p pp_ilst (BitSet.Int.list_of_packed p) m
                              pp_ilst (BitSet.Int.list_of_packed m) cst t state res;
            assert( state <> ~-1 );
            res,cst)

    let cost m tx ty = 
        if use_cost_fn then
            get_cf m tx ty
        else
            let cm = get_cm m tx ty in 
            (fun a b -> cm.(a-1).(b-1))

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
        fst mem.(xlen-1).(ylen-1)


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
            | Root   -> Sequence.of_array (Array.of_list ((Sequence.get x 0)::one)),
                        Sequence.of_array (Array.of_list ((Sequence.get y 0)::two))
        in
        build_alignments [] [] ((Sequence.length x)-1) ((Sequence.length y)-1)

    (* [backtrace men x y] find the median of x and y *)
    let backtrace ?(filter_gap=true) model mem x y = 
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
            (* move down each column and update until nothing changes *)
            let run_col i j i_max =
                let rec run_col i j i_max = 
                    if i < i_max then begin
                        update_all i j;
                        run_col (i+1) (j) i_max
                    end else if i = i_max then begin
                        update_row i j
                    end
                in
                update_col i j;
                run_col (i+1) j i_max
            (* move across each row and update until nothing changes *)
            and run_row i j j_max = 
                let rec run_row i j j_max =
                    if j < j_max then begin
                        update_all i j;
                        run_row (i) (j+1) j_max
                    end else if j = j_max then begin
                        update_col i j
                    end
                in
                update_row i j;
                run_row i (j+1) j_max
            in
            (* If dolphins are so smart, why do they live in Igloos? *)
            let ob = barrier ok and nb = barrier nk in
            for i = 1 to (lenX-1) do
                (* ___ _______ ___
                 * |___|_______|___| ; update new sections : right by column
                 *  new   old   new  ;                     : left by row *)
                let old_j_max = min (lenY-1) (i+ob+(lenY-lenX))
                and new_j_max = min (lenY-1) (i+nb+(lenY-lenX))
                and new_j_min = max 1 (i - nb) in
                for j = (old_j_max+1) to new_j_max do
                    run_col i j (i+ob)
                done;
                run_row i new_j_min new_j_max;
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
            for i = 1 to (lenX-1) do
                let j_min = max 1 (i - b - 1)
                and j_max = min (lenY-1) (i+b+(lenY-lenX)) in
                update_row i j_min;
                for j = j_min+1 to j_max-1 do
                    update_all i j
                done;
                update_col i j_max;
            done;
            update k

        (* this is to update k and matrix until ending condition *)
        and update k = 
            let mat_k = fst (snd (mem.(lenX-1).(lenY-1))) in
            if debug_mem then print_mem mem;
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
            let t_cost,t_count = 
                List.fold_left
                    (fun acc m ->
                        List.fold_left
                            (fun (cst,cnt) p -> (cst +. cost_matrix.{p,m}, (cnt+1)))
                            (acc)
                            (BitSet.Int.list_of_packed p))
                    (0.0,0)
                    (BitSet.Int.list_of_packed m)
            in
            assert( t_count >= 1);
            if t_count > 1 then t_cost /. (float_of_int t_count)
                           else t_cost
        in
        assert( (Sequence.length p) = (Sequence.length m) );
        let r = 
            Sequence.foldi
                (fun acc pos _ -> 
                    let p_i = Sequence.get p pos 
                    and m_i = Sequence.get m pos in
                    acc +. cost_single p_i m_i)
                (0.0)
                (m)
        in
        Numerical.FPInfix.reset ();
        r


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
        backtrace model mem s1 s2

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
            let med   = backtrace model mem s1 s2 in
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
        c,backtrace model mem s1 s2

    let median_2 s1 s2 model t1 t2 mem : s = snd (algn s1 s2 model t1 t2 mem)

    let median_2_cost s1 s2 model t1 t2 mem : float * s = algn s1 s2 model t1 t2 mem

    let closest ~p ~m model t mem : s * float =
        let alph= Alphabet.explote model.alph 1 0 in
        let gap = Alphabet.get_gap alph in
        if debug_aln then begin
            printf "\nP: ";print_s p alph;
            printf "\nM: ";print_raw m; (* raw; SM is not single *)
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
            printf " -%f(%d)-> " cst gap; print_s new_m alph; print_newline ();
        end;
        res


    (* pseudo 3d operation to find a median of three *)
    let readjust s1 s2 s3 model t1 t2 t3 mem =
        let algn s1 s2 t1 t2 : float * s = algn s1 s2 model t1 t2 mem in
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
                    false, c123, closest s3 s12 model t3 mem, c123
                else 
                    true, c132, closest s2 s13 model t2 mem, c123
            else if c231 < c132 then
                true, c231, closest s1 s23 model t1 mem, c123 
            else 
                true, c132, closest s2 s13 model t2 mem, c123 
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

    (** choose a standard direction for consistent alignments *)
    let choose_dir dirs = 
        let is_del = function | Delete _ -> true | _ -> false
        and is_aln = function | Align  _ -> true | _ -> false
        and is_ins = function | Insert _ -> true | _ -> false  in
             if List.exists is_del dirs then List.find is_del dirs
        else if List.exists is_aln dirs then List.find is_aln dirs
        else if List.exists is_ins dirs then List.find is_ins dirs
        else match dirs with | [Root] -> Root | _ -> assert false
 
    let print_mem (mem:floatmem) = 
        let dir_string = function
            | Root     -> "X" | Align _  -> "\\"
            | Delete _ -> "|" | Insert _ -> "-"
        and get_assignment = function
            | Root     -> 0 
            | Align a | Delete a | Insert a -> a
        in
        for i = 0 to (Array.length mem)-1 do
            for j = 0 to (Array.length mem.(i))-1 do
                let cost,(indel,dirs) = mem.(i).(j) in
                printf "|  %s%d %02.3f[%d]  " 
                    (dir_string (choose_dir dirs)) (List.length dirs) (abs_float cost) 
                    (get_assignment (choose_dir dirs));
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
        MlModel.check_metricity m.static t1 t2 cost1 cost2;
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
                                        c_min := (c, i::(snd !c_min))
                                    else if c < (fst !c_min) then
                                        c_min := (c,[i])
                                done;
                                !c_min)
                            acc
                            ys)
                    (infinity,[])
                    xs
            in
            if debug_aln || (0 = List.length states) then begin
                let packed = BitSet.Int.packed_of_list states in
                printf "Cost: %d(%a) ->%d(%a)<- %d(%a) = %f\n%!"
                    x_i pp_ilst xs packed pp_ilst states y_i pp_ilst ys cst;
            end;
            match classify_float cst, states with
            | _, [] ->
                failwithf "%d -(%f-%f)- %d: has no median" x_i t1 t2 y_i
            | FP_infinite, _ | FP_nan, _ ->
                failwithf "%d -- %d: returned infinite cost (%f,%f)" x_i y_i t1 t2
            | _, _  -> cst, BitSet.Int.packed_of_list states

    let get_cm model tx ty =
        let f = create_mpl_cost_fn model tx ty in
        let s = Alphabet.distinct_size model.alph in
        let m = Array.make_matrix s s (~-.1.0,0) in
        for i = 0 to (Array.length m)-1 do
            for j = i to (Array.length m.(i))-1 do
                let cost_state = f (i+1) (j+1) in
                m.(i).(j) <- cost_state;
                m.(j).(i) <- cost_state;
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
            if debug_aln || ~-1 = state then
                printf "%d -- p:%02d(%a) m:%02d(%a)\t-(%f)->%02d(%02d)\n%!"
                              i p pp_ilst (BitSet.Int.list_of_packed p) m 
                              pp_ilst (BitSet.Int.list_of_packed m) cst state res;
            assert( state <> ~-1 );
            res,cst)

    let cost m tx ty = 
        if use_cost_fn then
            let cf = create_mpl_cost_fn m tx ty in
            (fun a b -> cf a b)
        else
            let cm = get_cm m tx ty in 
            (fun a b -> cm.(a-1).(b-1) )

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
                min3_ ((fst mem.(i-1).(j))   +. dcst) (Delete sd)
                      ((fst mem.(i).(j-1))   +. icst) (Insert si)
                      ((fst mem.(i-1).(j-1)) +. acst) (Align  sa)
            end
        in
        for i = 0 to xlen - 1 do
            for j = 0 to ylen - 1 do
                let c,ds = get_cost i j in
                mem.(i).(j) <- c,(0,ds)
            done;
        done;
        Numerical.FPInfix.reset ();
        fst mem.(xlen-1).(ylen-1)


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
            | Root     -> Sequence.of_array (Array.of_list ((Sequence.get x 0)::one)),
                          Sequence.of_array (Array.of_list ((Sequence.get y 0)::two))
        in
        build_alignments [] [] ((Sequence.length x)-1) ((Sequence.length y)-1)

    (* [backtrace men x y] find the median of x and y *)
    let backtrace ?(filter_gap=true) model mem x y = 
        let get x i = Sequence.get x i in
        let get_direction i j = mem.(i).(j) --> snd --> snd --> choose_dir in
        let rec build_median acc i j = match get_direction i j with
            | Align  s -> build_median (s::acc) (i-1) (j-1)
            | Delete s -> build_median (s::acc) (i-1) j
            | Insert s -> build_median (s::acc) i (j-1)
            | Root     -> Sequence.of_array (Array.of_list ((get x 0)::acc))
        in
        let m = build_median [] ((Sequence.length x)-1) ((Sequence.length y)-1) in
        let s = 
            if filter_gap
                then remove_gaps (Alphabet.get_gap model.alph) m
                else m
        in
        s

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
            let aln,at,sa = get_cost (i-1) (j-1) (Sequence.get x i) (Sequence.get y j)
            and del,dt,sd = get_cost (i-1) (j)   (gap) (Sequence.get y j)
            and ins,it,si = get_cost (i)   (j-1) (Sequence.get x i) (gap) in
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
                    (aln,(at,[Align sa]))
                    [(del,(dt,Delete sd)); (ins,(it,Insert si))]
            in
            mem.(i).(j) <- m

        (* Same as above, but will not look at the node to the right (j-1) *)
        and update_row i j = 
            let aln,at,sa = get_cost (i-1) (j-1) (Sequence.get x i) (Sequence.get y j)
            and del,dt,sd = get_cost (i-1) (j) gap (Sequence.get y j) in
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
            mem.(i).(j) <- m
   
        (* Same as above, but will not look at the node above (i-1) *)
        and update_col i j =
            let aln,at,sa = get_cost (i-1) (j-1) (Sequence.get x i) (Sequence.get y j)
            and ins,it,si = get_cost (i) (j-1)   (Sequence.get x i) gap in
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
            mem.(i).(j) <- m
        in

        (* Update the Ukkonen barrier cells; The new area is bounded by two
         * strips along the edges of the diagonal band, perpendicular to each
         * other. We move along the barrier, down, and across until a node does
         * not change it's value and stop. *)
        let rec update_matrix (ok:int) (nk:int): unit =
            (* move down each column and update until nothing changes *)
            let run_col i j i_max =
                let rec run_col i j i_max = 
                    if i < i_max then begin
                        update_all i j;
                        run_col (i+1) (j) i_max
                    end else if i = i_max then begin
                        update_row i j
                    end
                in
                update_col i j;
                run_col (i+1) j i_max
            (* move across each row and update until nothing changes *)
            and run_row i j j_max = 
                let rec run_row i j j_max =
                    if j < j_max then begin
                        update_all i j;
                        run_row (i) (j+1) j_max
                    end else if j = j_max then begin
                        update_col i j
                    end
                in
                update_row i j;
                run_row i (j+1) j_max
            in
            (* If dolphins are so smart, why do they live in Igloos? *)
            let ob = barrier ok and nb = barrier nk in
            for i = 1 to (lenX-1) do
                (* ___ _______ ___
                 * |___|_______|___| ; update new sections : right by column
                 *  new   old   new  ;                     : left by row *)
                let old_j_max = min (lenY-1) (i+ob+(lenY-lenX))
                and new_j_max = min (lenY-1) (i+nb+(lenY-lenX))
                and new_j_min = max 1 (i - nb) in
                for j = (old_j_max+1) to new_j_max do
                    run_col i j (i+ob)
                done;
                run_row i new_j_min new_j_max;
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
            for i = 1 to (lenX-1) do
                let j_min = max 1 (i - b - 1)
                and j_max = min (lenY-1) (i+b+(lenY-lenX)) in
                update_row i j_min;
                for j = j_min+1 to j_max-1 do
                    update_all i j
                done;
                update_col i j_max;
            done;
            update k

        (* this is to update k and matrix until ending condition *)
        and update k = 
            let mat_k = fst (snd (mem.(lenX-1).(lenY-1))) in
            if debug_mem then print_mem mem;
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
            let t_cost,t_count = 
                List.fold_left
                    (fun acc m ->
                        List.fold_left
                            (fun (cst,cnt) p -> (cst +. cost_matrix.{p,m}, (cnt+1)))
                            (acc)
                            (BitSet.Int.list_of_packed p))
                    (0.0,0)
                    (BitSet.Int.list_of_packed m)
            in
            assert( t_count >= 1 );
            if t_count > 1 then t_cost /. (float_of_int t_count)
                           else t_cost
        in
        assert( (Sequence.length p) = (Sequence.length m) );
        let r = 
            Sequence.foldi
                (fun acc pos _ -> 
                    let p_i = Sequence.get p pos 
                    and m_i = Sequence.get m pos in
                    acc +. cost_single p_i m_i)
                (0.0)
                (m)
        in
        Numerical.FPInfix.reset ();
        r

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
        backtrace model mem s1 s2

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
            let med   = backtrace model mem s1 s2 in
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
        let b = backtrace model mem s1 s2 in
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
                remove_gaps gap m, aln_cost_2 p m model t
            else
                let paln, maln, cst = align_2 p m model 0.0 t mem in
                assert( Sequence.length paln = Sequence.length maln );
                let m = Sequence.mapi (get_closest paln) maln in
                remove_gaps gap m, aln_cost_2 paln m model t
        in
        if debug_aln then begin
            printf "\nP: ";print_raw p;
            printf "\nM: ";print_raw m; (* raw; SM is not single *)
            printf " -%f-> " c; print_s s_new model.alph;
            print_newline ();
        end;
        res

    (* requires not implemented functions *)
    let readjust s1 s2 s3 model t1 t2 t3 mem =
        let algn s1 s2 t1 t2 : float * s = algn s1 s2 model t1 t2 mem in
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
                    false, c123, closest s3 s12 model t3 mem, c123
                else 
                    true, c132, closest s2 s13 model t2 mem, c123
            else if c231 < c132 then
                true, c231, closest s1 s23 model t1 mem, c123 
            else 
                true, c132, closest s2 s13 model t2 mem, c123 
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

    let print_mem (mem:floatmem) = 
        for i = 0 to (Array.length mem)-1 do
            for j = 0 to (Array.length mem.(i))-1 do
                printf "| %02.3f " (abs_float mem.(i).(j));
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
            if debug_aln then
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
    let readjust _ _ _ _ _ _ _ _    = failwith "not implemented"
    let get_cm _ _ _                = failwith "not implemented"
    let cost _ _ _                  = failwith "not implemented"
    let get_cf _ _ _                = failwith "not implemented"
    let full_cost_2 _ _ _ _ _ _     = failwith "not implemented"
    let create_edited_2 _ _ _ _ _ _ = failwith "not implemented"
    let clip_align_2 ?first_gap _ _ _ _ _ = failwith "not implemented"
    let align_2 ?first_gap _ _ _ _ _ _ = failwith "not implemented"
    let backtrace ?filter_gap _ _ _ _ = failwith "not implemented" 

end

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
        let med       = FloatAlign.backtrace ~filter_gap model mem faln1 faln2 in
        let _,ocst    = FloatAlign.optimize faln1 faln2 model (bl1+.bl2) mem in
        cst, ocst, (FloatAlign.seq_of_s a1), (FloatAlign.seq_of_s a2), (FloatAlign.seq_of_s med)
    and mpl_cost,mpl_opt_cost,mpl_seq1,mpl_seq2,mpl_median =
        let dmpl1     = MPLAlign.s_of_seq seq1
        and dmpl2     = MPLAlign.s_of_seq seq2 in
        let mem       = MPLAlign.get_mem dmpl1 dmpl2 in
        let a1,a2,cst = MPLAlign.align_2 dmpl1 dmpl2 model bl1 bl2 mem in
        let med       = MPLAlign.backtrace ~filter_gap model mem dmpl1 dmpl2 in
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
    in
    if alignments then begin
        Printf.fprintf channel "%a  %a\n%!" 
                        (pp_seq model) mpl_seq1 (pp_seq model) flk_seq1;
        Printf.fprintf channel "%a  %a\n%!"
                        (pp_seq model) mpl_seq2 (pp_seq model) flk_seq2;
    end;
    Printf.fprintf channel "%f %f %f %f %f %f  %a  %a\n%!" 
                   mal_cost     mpl_cost     flk_cost
                   mal_opt_cost mpl_opt_cost flk_opt_cost
                   (pp_seq model) mpl_median (pp_seq model) flk_median;
