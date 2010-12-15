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

(* Debug variables/ combinators *)
let (-->) a b = b a
let debug_mem = false
let debug_aln = false
let failwithf format = Printf.ksprintf (failwith) format

let to_char_list string = 
    let rec to_ i len acc = 
        if i = len then List.rev acc
                   else to_ (i+1) len ((String.get string i)::acc)
    in
    to_ 0 (String.length string) []

let sequence_of_string seq alph = 
    seq --> to_char_list
        --> List.map (String.make 1)
        --> List.map (fun x -> Alphabet.match_base x alph)
        --> Array.of_list
        --> Sequence.of_array

let pp_ilst chan lst =
    List.iter (fun x -> Printf.fprintf chan "%d| " x) lst

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
    val print_cm    : MlModel.model -> float -> unit

    (* cost matrix; maps states to their cost and possible states *)
    val get_cm : MlModel.model -> float -> float -> (int -> int -> float * int)

    (* 2d operations *)
    val cost_2          : ?deltaw:int -> s -> s -> MlModel.model -> float -> float -> floatmem -> float
    val verify_cost_2   : float -> s -> s -> MlModel.model -> float -> float -> floatmem -> float
    val c_cost_2        : s -> s -> MlModel.model -> float -> float -> floatmem -> int -> float
    val create_edited_2 : s -> s -> MlModel.model -> float -> float -> floatmem -> s * s
    val align_2         : ?first_gap:bool -> s -> s -> MlModel.model -> float -> float -> floatmem -> s * s * float
    val median_2        : s -> s -> MlModel.model -> float -> float -> floatmem -> s
    val median_2_cost   : s -> s -> MlModel.model -> float -> float -> floatmem -> float * s
    val full_median_2   : s -> s -> MlModel.model -> float -> float -> floatmem -> s
    val gen_all_2       : s -> s -> MlModel.model -> float -> float -> floatmem -> s * s * float * s

    (* uppass operations *)
    val closest  : p:s -> m:s -> MlModel.model -> float -> floatmem -> s * float

    (* (pseudo) 3d operations *)
    val readjust : s -> s -> s -> MlModel.model -> float -> float -> float -> floatmem -> float * s * bool

end 


module FloatAlign : A = struct

    type s        = Sequence.s
    type dir      = Root | Align | Delete | Insert
    type floatmem = ( float * ( int * dir list ) ) array array

    external s_of_seq : Sequence.s -> s = "%identity"
    external seq_of_s : s -> Sequence.s = "%identity"

    let s_of_seq = s_of_seq
    let seq_of_s = seq_of_s

    let print_mem (mem:floatmem) = 
        let dirToString = function
            | Root   -> "X" | Align  -> "\\"
            | Delete -> "|" | Insert -> "-"
        in
        for i = 0 to (Array.length mem)-1 do
            for j = 0 to (Array.length mem.(i))-1 do
                let cost,(indel,dirs) = mem.(i).(j) in
                Printf.printf "|  %s %02.3f[%d]  " 
                    (dirToString (List.hd dirs)) (abs_float cost) indel;
            done;
            Printf.printf "|\n%!";
        done;
        print_newline ()

    let print_raw seq =
        for i = 0 to (Sequence.length seq)-1 do
            Printf.printf "|%d" (Sequence.get seq i);
        done;
        Printf.printf "|%!";
        ()

    let print_s (seq:s) (a) =
        for i = 0 to (Sequence.length seq)-1 do
            try Printf.printf "%s" (Alphabet.match_code (Sequence.get seq i) a);
            with (Alphabet.Illegal_Code _ ) -> Printf.printf "%d" (Sequence.get seq i);
        done;
        ()

    let print_cm m t = 
        let mat = MlModel.compose m t in
        for i = 0 to (Bigarray.Array2.dim1 mat) - 1 do
            for j = 0 to (Bigarray.Array2.dim2 mat) - 1 do
                Printf.printf " [%f] " (~-. (log mat.{i,j}));
            done;
            print_newline ();
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

    (* minimum of three *)
    let min3 a b c = min a (min b c)

    (* minimum of three with annotations *)
    let min3_ a at b bt c ct : float * dir list =
        if a > b then begin
            if b > c then (c,[ct]) else if b = c then (b,[bt;ct]) else (b,[bt])
        end else if a = b then begin
            if b = c then (a,[at;bt;ct]) else if a > c then (c,[ct]) else (a,[at;bt])
        end else begin
            if c > a then (a,[at]) else if c = a then (c,[at;ct]) else (c,[ct])
        end

    let create_cost_matrix model t = 
        let mat = MlModel.compose model t in
        for i = 0 to (Bigarray.Array2.dim1 mat) - 1 do
            for j = 0 to (Bigarray.Array2.dim2 mat) - 1 do
                mat.{i,j} <- ~-. (log mat.{i,j})
            done;
        done;
        mat


    (* [create_align_cost_fn m t] Compose the model [m] into a cost matrix with
     * branch length [t] *)
    let create_align_cost_fn m tx ty =
        let cost_matrix = create_cost_matrix m (tx+.ty) in
        fun x_i y_i ->
            let fn (cst,n) x y =
                try (cst +. cost_matrix.{x,y},n+1)
                with | _ -> failwithf "Cannot find from cost matrix: %d, %d" x y
            and xs = MlModel.list_of_packed x_i 
            and ys = MlModel.list_of_packed y_i in
            let cst,n =
                List.fold_left
                    (fun acc x ->
                        List.fold_left (fun acc y -> fn acc x y) acc ys)
                    (0.0,0) (xs)
            in
            if debug_aln then
                Printf.printf "Cost: %d(%a) ->%f/%d<- %d(%a)\n%!" 
                    x_i pp_ilst (MlModel.list_of_packed x_i) cst n
                    y_i pp_ilst (MlModel.list_of_packed y_i);
            cst /. (float_of_int n)

    let get_cm m t1 t2 = (fun a b -> create_align_cost_fn m t1 t2 a b, a land b)

    (* p is single; m is not; find in m the least coast to p *)
    let get_closest gap model t =
        let cost_matrix = create_cost_matrix model t in
        (fun ~i ~p ~m -> 
            (* Determine if this is a gap; or remove that state otherwise *)
            let nm =
                if m = gap || p = gap then m
                else if (0 <> p land gap) && (0 <> m land gap) then gap
                else m land (lnot gap)
            in
            assert(  p > 0 );
            assert( nm > 0 );
            (* now choose the best from parent *)
            let state,cst = 
                List.fold_left
                    (fun acc m ->
                        List.fold_left
                            (fun ((assgn,cst) as acc) p ->
                                let ncst = cost_matrix.{p,m} in
(*                                Printf.printf "\t Cost: %d -- %d = %f; old: %f\n" p m ncst cst;*)
                                if ncst < cst then (m,ncst) else acc)
                            (acc)
                            (MlModel.list_of_packed p))
                    ((~-1),infinity)
                    (MlModel.list_of_packed nm)
            in
            let res = 1 lsl state in
             if debug_aln then
                Printf.printf "%d -- p:%02d(%a) m:%02d/%02d(%a)\t-(%f)->%02d(%02d)\n%!"
                              i p pp_ilst (MlModel.list_of_packed p) m nm 
                              pp_ilst (MlModel.list_of_packed nm) cst state res;
            assert( state <> ~-1 );
            res)


    (* [align_2 mem x y m t] Align the sequence [x] and [y] with the cost
     * matrix from [m] and branch length [t], and completely fills the matrix of
     * [mem]. Return minimum cost; the alignment can be recovered from [mem] *)
    let align_2 (mem:floatmem) (x:s) (y:s) (m:MlModel.model) (tx:float) (ty:float) = 
        let cost = create_align_cost_fn m tx ty in
        let gap = Alphabet.get_gap (Alphabet.explote m.MlModel.alph 1 0) in
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
        let xlen = Sequence.length x and ylen = Sequence.length y in
        for i = 0 to xlen - 1 do
            for j = 0 to ylen - 1 do
                let c,ds = get_cost i j in
                mem.(i).(j) <- c,(0,ds)
            done;
        done;
        fst mem.(xlen-1).(ylen-1)


    (* [alignments mem x y m] builds the alignment of [x] and [y] from the
     * alignment cost matrix, [mem]; the only purpose for the model [m] is to
     * obtain the alphabet, and gap character. *)
    let alignments mem x y m =
        let get_direction i j = mem.(i).(j) --> snd --> snd --> List.hd
        and gap = Alphabet.get_gap (Alphabet.explote m.MlModel.alph 1 0) in
        let rec build_alignments one two i j = match get_direction i j with
            | Align  -> build_alignments ((Sequence.get x i)::one) ((Sequence.get y j)::two) (i-1) (j-1)
            | Delete -> build_alignments ((Sequence.get x i)::one) (gap::two) (i-1) (j)
            | Insert -> build_alignments (gap::one) ((Sequence.get y j)::two) (i) (j-1) 
            | Root   -> Sequence.of_array (Array.of_list ((Sequence.get x 0)::one)),
                        Sequence.of_array (Array.of_list ((Sequence.get y 0)::two))
        in
        build_alignments [] [] ((Sequence.length x)-1) ((Sequence.length y)-1)

    (* [backtrace men x y] find the median of x and y *)
    let backtrace mem x y = 
        let get x i = Sequence.get x i
        and union x y i j = (Sequence.get x i) lor (Sequence.get y j) in
        let get_direction i j = mem.(i).(j) --> snd --> snd --> List.hd in
        let rec build_median acc i j = match get_direction i j with
            | Align  -> build_median ((union x y i j)::acc) (i-1) (j-1)
            | Delete -> build_median ((get x i)::acc) (i-1) j
            | Insert -> build_median ((get y j)::acc) i (j-1)
            | Root   -> Sequence.of_array (Array.of_list ((get x 0)::acc))
        in
        build_median [] ((Sequence.length x)-1) ((Sequence.length y)-1)

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
        if (Array.length mem) < lenX then
            failwithf "FloatSequence.FloatAlign; not enough memory in scratch space; Has %d, requires %d"
                        (Array.length mem) lenX;
        if (Array.length mem.(0)) < lenY then
            failwithf "FloatSequence.FloatAlign; not enough memory in scratch space; Has %d, requires %d"
                        (Array.length mem.(0)) lenY;

        (* obtain the gap representation for this model *)
        let gap = Alphabet.get_gap (Alphabet.explote m.MlModel.alph 1 0) in

        (* Obtain the cost matrix for the alignment *)
        let cost_fn = create_align_cost_fn m tx ty in

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
        initial_matrix ();
        fst (mem.(lenX-1).(lenY-1))


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

    let c_cost_2 s1 s2 model t1 t2 mem delta : float = failwith "FloatAlign.c_cost_2 not implemented"

    let full_median_2 s1 s2 model t1 t2 mem : s =
        if debug_mem then clear_mem mem;
        let s1,s2,t1,t2 = 
            if Sequence.length s1 <= Sequence.length s2 
                then s1,s2,t1,t2 else s2,s1,t2,t1
        in
        ignore (ukkonen_align_2 mem s1 s2 model t1 t2);
        backtrace mem s1 s2

    let verify_cost_2 cost1 s1 s2 model t1 t2 mem : float =
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

    let gen_all_2 s1 s2 model t1 t2 mem =
        if debug_mem then clear_mem mem;
        let gen_all_2 s1 s2 t1 t2 = 
            let cost  = ukkonen_align_2 mem s1 s2 model t1 t2 in
            let med   = backtrace mem s1 s2 in
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
        c,backtrace mem s1 s2

    let median_2 s1 s2 model t1 t2 mem : s = snd (algn s1 s2 model t1 t2 mem)

    let median_2_cost s1 s2 model t1 t2 mem : float * s = algn s1 s2 model t1 t2 mem

    let closest ~p ~m model t mem : s * float =
        let alph= Alphabet.explote model.MlModel.alph 1 0 in
        let gap = Alphabet.get_gap alph in
        if debug_aln then begin
            Printf.printf "\nP: ";print_s p alph;
            Printf.printf "\nM: ";print_raw m; (* raw; SM is not single *)
        end;
        let remove_gaps seq =
            let remove_gaps seq base = 
                if base <> gap then 
                    let () = Sequence.prepend seq base in seq
                else seq
            in
            let res = Sequence.fold_right (remove_gaps)
                                          (Sequence.create (Sequence.length seq)) (seq)
            in
            Sequence.prepend res gap;
            res
        and mask_gaps seq gap =
            let mask = lnot gap in
            Sequence.mapi (fun x p -> if p > 0 then x land mask else x) seq
        and get_closest par : int -> int -> int =
            let gc = get_closest gap model t in
            (fun m pos -> gc ~i:pos ~p:(Sequence.get par pos) ~m)
        in
        let (new_m,cst) as res =
            if Sequence.is_empty m gap then
                m, 0.0
            else if 0 = Sequence.compare p m then
                let masked = mask_gaps m gap in
                Sequence.mapi (get_closest masked) masked --> remove_gaps, 0.0
            else
                let paln, maln, cst = align_2 p m model (t/.2.0) (t/.2.0) mem in
                assert( Sequence.length paln = Sequence.length maln );
                Sequence.mapi (get_closest paln) maln --> remove_gaps, cst
        in
        if debug_aln then begin
            Printf.printf " -%f(%d)-> " cst gap; print_s new_m alph; print_newline ();
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
                Printf.printf "Cost123: %f\tCost231: %f\tCost132: %f\n" c123 c231 c132;
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

end


module MPLAlign : A = struct

    type s   = Sequence.s
    type dir = Root | Align  of int
                    | Delete of int
                    | Insert of int

    type floatmem = ( float * ( int * dir list ) ) array array

    external s_of_seq : Sequence.s -> s = "%identity"
    external seq_of_s : s -> Sequence.s = "%identity"

    let print_mem (mem:floatmem) = 
        let dirToString = function
            | Root     -> "X" | Align _  -> "\\"
            | Delete _ -> "|" | Insert _ -> "-"
        in
        for i = 0 to (Array.length mem)-1 do
            for j = 0 to (Array.length mem.(i))-1 do
                let cost,(indel,dirs) = mem.(i).(j) in
                Printf.printf "|  %s %02.3f[%d]  " 
                    (dirToString (List.hd dirs)) (abs_float cost) indel;
            done;
            Printf.printf "|\n%!";
        done;
        print_newline ()

    let print_raw seq =
        for i = 0 to (Sequence.length seq)-1 do
            Printf.printf "|%d" (Sequence.get seq i);
        done;
        Printf.printf "|\n%!";
        ()

    let print_s (seq:s) (a) =
        for i = 0 to (Sequence.length seq)-1 do
            try Printf.printf "%s" (Alphabet.match_code (Sequence.get seq i) a);
            with (Alphabet.Illegal_Code _) -> Printf.printf "%d" (Sequence.get seq i);
        done;
        ()

    let print_cm m t = 
        let mat = MlModel.compose m t in
        for i = 0 to (Bigarray.Array2.dim1 mat) - 1 do
            for j = 0 to (Bigarray.Array2.dim2 mat) - 1 do
                Printf.printf " [%f] " (~-. (log mat.{i,j}));
            done;
            print_newline ();
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

    (* minimum of three *)
    let min3 a b c = min a (min b c)

    (* minimum of three with annotations *)
    let min3_ a at b bt c ct : float * dir list =
        if a > b then begin
            if b > c then (c,[ct]) else if b = c then (b,[bt;ct]) else (b,[bt])
        end else if a = b then begin
            if b = c then (a,[at;bt;ct]) else if a > c then (c,[ct]) else (a,[at;bt])
        end else begin
            if c > a then (a,[at]) else if c = a then (c,[at;ct]) else (c,[ct])
        end


    (* [create_mpl_cost_fn m t1 t2]
     * In this cost function, an analogue to MPL, we assign the node to the
     * minimum cost of transforming each child to that node, across their
     * respective branch lengths. *)
    let create_mpl_cost_fn ?(epsilon=1e-5) m t1 t2 =
        let (=.) a b = (abs_float (a-.b)) < epsilon in
        let cost1,cost2 = 
            let mat1 = MlModel.compose m t1 and mat2 = MlModel.compose m t2 in
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
        let med_cost xe ye me = cost1.{xe,me} +. cost2.{ye,me} in
        (* return function to calculate costs of polymorphisms *)
        fun x_i y_i ->
            let xs = MlModel.list_of_packed x_i 
            and ys = MlModel.list_of_packed y_i in
            let cst,states = 
                List.fold_left 
                    (fun acc x ->
                        List.fold_left 
                            (fun acc y ->
                                let c_min = ref acc in
                                for i = 0 to (m.MlModel.alph_s)-1 do
                                    let c = med_cost x y i in
                                    if c < (fst !c_min) then 
                                        c_min := (c,[i])
                                    else if (fst !c_min) =. c then 
                                       c_min := (c, i::(snd !c_min))
                                done;
                                !c_min)
                            acc 
                            ys)
                    (infinity,[])
                    xs
            in
            match classify_float cst with
            | FP_infinite | FP_nan -> failwith "returned infinite cost"
            | _                    -> cst, MlModel.packed_of_list states

    let get_cm m t1 t2 = (fun a b -> create_mpl_cost_fn m t1 t2 a b)

    let single_cost_fn m t = 
        let mat = 
            let mat = MlModel.compose m t in
            for i = 0 to (Bigarray.Array2.dim1 mat) - 1 do
                for j = 0 to (Bigarray.Array2.dim2 mat) - 1 do
                    mat.{i,j} <- ~-. (log mat.{i,j});
                done;
            done;
            mat
        in
        fun x y ->
            assert( 0 = ( x land (x-1)));
            let x = List.hd (MlModel.list_of_packed x) in
            List.fold_left
                (fun ((_,cst) as acc) yi ->
                    let ncst = mat.{yi,x} in
                    if ncst < cst then (x,ncst) else acc)
                (x,max_float)
                (MlModel.list_of_packed y)


    (* p is single; m is not; find in m the least coast to p *)
    let get_closest i gap cst_fn ~p ~m : int =
        let m = 
            if m = gap || p = gap then m
            else if (0 <> p land gap) && (0 <> m land gap) then gap
            else m land (lnot gap)
        in
        assert( m > 0 ); 
        let state,_ = 
            List.fold_left
                (fun ((assn,cst) as acc) x -> 
                    let ncst = cst_fn p x in
                    if ncst < cst then (x,ncst) else acc)
                (~-1,infinity)
                (MlModel.list_of_packed m)
        in
        assert( state <> ~-1 );
        let res = 1 lsl state in
        if debug_aln then
            Printf.printf "%d -- p:%02d m:%02d\t%a\t-(%d)->%02d\n" 
                          i p m pp_ilst (MlModel.list_of_packed m) state res;
        res


    (* [align_2 mem x y m t] Align the sequence [x] and [y] with the cost
     * matrix from [m] and branch length [t], and completely fills the matrix of
     * [mem]. Return minimum cost; the alignment can be recovered from [mem] *)
    let align_2 (mem:floatmem) (x:s) (y:s) (m:MlModel.model) (tx:float) (ty:float) = 
        let cost = create_mpl_cost_fn m tx ty in
        let gap = Alphabet.get_gap (Alphabet.explote m.MlModel.alph 1 0) in
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
        let xlen = Sequence.length x and ylen = Sequence.length y in
        for i = 0 to xlen - 1 do
            for j = 0 to ylen - 1 do
                let c,ds = get_cost i j in
                mem.(i).(j) <- c,(0,ds)
            done;
        done;
        fst mem.(xlen-1).(ylen-1)


    (* [alignments mem x y m] builds the alignment of [x] and [y] from the
     * alignment cost matrix, [mem]; the only purpose for the model [m] is to
     * obtain the alphabet, and gap character. *)
    let alignments mem x y m =
        let get_direction i j = mem.(i).(j) --> snd --> snd --> List.hd
        and gap = Alphabet.get_gap (Alphabet.explote m.MlModel.alph 1 0) in
        let rec build_alignments one two i j = match get_direction i j with
            | Align  _ -> build_alignments ((Sequence.get x i)::one) ((Sequence.get y j)::two) (i-1) (j-1)
            | Insert _ -> build_alignments (gap::one) ((Sequence.get y j)::two) (i) (j-1) 
            | Delete _ -> build_alignments ((Sequence.get x i)::one) (gap::two) (i-1) (j)
            | Root     -> Sequence.of_array (Array.of_list ((Sequence.get x 0)::one)),
                          Sequence.of_array (Array.of_list ((Sequence.get y 0)::two))
        in
        build_alignments [] [] ((Sequence.length x)-1) ((Sequence.length y)-1)

    (* [backtrace men x y] find the median of x and y *)
    let backtrace mem x y = 
        let get x i = Sequence.get x i in
        let get_direction i j = mem.(i).(j) --> snd --> snd --> List.hd in
        let rec build_median acc i j = match get_direction i j with
            | Align  s -> build_median (s::acc) (i-1) (j-1)
            | Delete s -> build_median (s::acc) (i-1) j
            | Insert s -> build_median (s::acc) i (j-1)
            | Root     -> Sequence.of_array (Array.of_list ((get x 0)::acc))
        in
        build_median [] ((Sequence.length x)-1) ((Sequence.length y)-1)

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
        if (Array.length mem) < lenX then
            failwithf "FloatSequence.FloatAlign; not enough memory in scratch space; Has %d, requires %d"
                        (Array.length mem) lenX;
        if (Array.length mem.(0)) < lenY then
            failwithf "FloatSequence.FloatAlign; not enough memory in scratch space; Has %d, requires %d"
                        (Array.length mem.(0)) lenY;

        (* obtain the gap representation for this model *)
        let gap = Alphabet.get_gap (Alphabet.explote m.MlModel.alph 1 0) in

        (* Obtain the cost matrix for the alignment *)
        let cost_fn = create_mpl_cost_fn m tx ty in

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
        initial_matrix ();
        fst (mem.(lenX-1).(lenY-1))



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

    let c_cost_2 s1 s2 model t1 t2 mem delta : float = failwith "FloatAlign.c_cost_2 not implemented"

    let full_median_2 s1 s2 model t1 t2 mem : s =
        if debug_mem then clear_mem mem;
        let s1,s2,t1,t2 = 
            if Sequence.length s1 <= Sequence.length s2 
                then s1,s2,t1,t2 else s2,s1,t2,t1
        in
        ignore (align_2 mem s1 s2 model t1 t2);
        backtrace mem s1 s2

    let verify_cost_2 cost1 s1 s2 model t1 t2 mem : float =
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

    let gen_all_2 s1 s2 model t1 t2 mem =
        if debug_mem then clear_mem mem;
        let gen_all_2 s1 s2 t1 t2 = 
            let cost  = ukkonen_align_2 mem s1 s2 model t1 t2 in
            let med   = backtrace mem s1 s2 in
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
        c,backtrace mem s1 s2

    let median_2 s1 s2 model t1 t2 mem : s = snd (algn s1 s2 model t1 t2 mem)

    let median_2_cost s1 s2 model t1 t2 mem : float * s = algn s1 s2 model t1 t2 mem

    let closest ~p ~m model t mem : s * float =
        let alph= Alphabet.explote model.MlModel.alph 1 0 in
        let gap = Alphabet.get_gap alph in
        let remove_gaps seq =
            let remove_gaps seq base = 
                if base <> gap then 
                    let () = Sequence.prepend seq base in seq
                else seq
            in
            let res = 
                Sequence.fold_right (remove_gaps)
                                    (Sequence.create (Sequence.length seq)) 
                                    (seq)
            in
            Sequence.prepend res gap;
            res
        and mask_gaps seq gap =
            let mask = lnot gap in
            Sequence.mapi (fun x p -> if p > 0 then x land mask else x) seq
        and get_closest par : int -> int -> int =
            let cst i j = 
                fst (create_mpl_cost_fn model (t/.2.0) (t/.2.0) i j) in
            (fun m pos -> 
                get_closest pos gap cst ~p:(Sequence.get par pos) ~m)
        in
        let (s_new,c) as res =
            if Sequence.is_empty m gap then
                m, 0.0
            else if 0 = Sequence.compare p m then
                let masked = mask_gaps m gap in
                Sequence.mapi (get_closest masked) masked --> remove_gaps, 0.0
            else
                let paln, maln, cst = align_2 p m model (t/.2.0) (t/.2.0) mem in
                assert( Sequence.length paln = Sequence.length maln );
                Sequence.mapi (get_closest paln) maln --> remove_gaps, cst
        in
        if debug_aln then begin
            Printf.printf "\nP: ";print_s p alph;
            Printf.printf "\nM: ";print_raw m; (* raw; SM is not single *)
                Printf.printf " -%f-> " c; print_s s_new alph;
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
                Printf.printf "Cost123: %f\tCost231: %f\tCost132: %f\n" c123 c231 c132;
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

end

