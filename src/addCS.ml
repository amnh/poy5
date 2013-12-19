(* POY 5.1.1. A phylogenetic analysis program using Dynamic Homologies.       *)
(* Copyright (C) 2013 Andrés Varón, Lin Hong, Nicholas Lucaroni, Ward Wheeler,*)
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

let () = SadmanOutput.register "AddCS" "$Revision: 3641 $"

let tests = false (* run additive self tests *)

let ( --> ) a b = b a

module type AdditiveInterface = 
  sig
    type ct

    val create : int array -> int array -> ct
    val copy : ct -> ct -> unit
    val clone : ct -> ct
    val compare_data : ct -> ct -> int

    val median : ct -> ct -> ct
    val distance : ct -> ct -> float
    val distance_2 : ct -> ct -> ct -> float
    val distance_median : ct -> ct -> float * ct
    val median_cost : ct -> float
    val median_3 : ct -> ct -> ct -> ct -> ct
    val full_unioni : ct -> ct -> ct -> unit
    val full_union : ct -> ct -> ct
    val mediani : ct -> ct -> ct -> unit

    val pos_set_state : ct -> int -> int -> int -> unit
    val pos_get_max : ct -> int -> int
    val pos_get_min : ct -> int -> int
    val pos_get_cost : ct -> int -> float
    val to_string : ct -> string
    val cardinal : ct -> int

    val print : ct -> unit
  end

module type Make = sig
    exception Exists
    exception Illegal_Arguments
    exception Duplicated
    exception Illegal_State
    exception Not_Found

    type t
    type c = int * int * int
    val of_array : c array -> int -> t
    val of_list : c list -> int -> t
    val to_list : t -> c list
    val to_list_with_cost : t -> (int * int * int * float) list
    val copy : t -> t -> unit 
    val clone : t -> t 
    val median : t option -> t -> t -> t
    val reroot_median : t -> t -> t
    val median_3 : t -> t -> t -> t -> t 
    val distance : t -> t -> float 
    val distance_2 : t -> t -> t -> float 
    val dist_2 : t -> t -> t -> float 
    val distance_median : t -> t -> float * t 
    val median_cost : t -> float 
    val compare_codes : t -> t -> int 
    val compare_data : t -> t -> int
    val set_state : t -> int -> int -> int -> t
    val get_max : t -> int -> int 
    val get_min : t -> int -> int 
    val get_cost : t -> int -> float
    val find_pos : t -> int -> int
    val codes : t -> int list
    val code_mem : int list option -> t -> bool
    val char_mem : int list option -> t -> bool
    val get_length : t -> int
    val cardinal : t -> int
    val deep_cardinal : t -> int
    val get_state : t -> int -> c
    val elt_code : c -> int
    val set_code : t -> int
    val code : t -> int
    val get_code : t -> int -> int
    val elements : t -> c list
    val map : (c -> c) -> t -> t
    val empty : int -> t
    val is_empty : t -> bool
    val mem : c -> t -> bool
    val code_exists : int -> t -> bool
    val add : c -> t -> t
    val singleton : c -> int -> t
    val remove : c -> t -> t
    val union : t -> t -> t  
    val inter : t -> t -> int -> t
    val diff : t -> t -> int -> t
    val subset : t -> t -> bool
    val equal : t -> t -> bool
    val iter : (c -> unit) -> t -> unit
    val fold : ('a -> c -> 'a) -> t -> 'a -> 'a
    val for_all : (c -> bool) -> t -> bool
    val exists : (c -> bool) -> t -> bool
    val filter : (c -> bool) -> t -> t
    val f_codes : t -> All_sets.Integers.t -> t
    val f_codes_comp : t -> All_sets.Integers.t -> t
    val partition : (c -> bool) -> t -> t * t
    val min_elt : t -> c
    val max_elt : t -> c
    val choose : t -> c
    val full_union : t -> t -> t -> unit 
    val split : c -> t -> t * bool * t
    val of_parser : Data.d -> ((Nexus.File.static_state * int) array * int) -> int -> t * int
    val to_string : t -> string 
    val state_to_xml : Pervasives.out_channel -> t -> Data.d -> unit
    val to_formatter : Xml.attributes -> t -> t option -> Data.d -> Xml.xml Sexpr.t list

    module Imperative : sig
        type it 
        type ic
        val create : t -> it
        val of_array : ic array -> int -> it
        val of_list : ic list -> int -> it
        val to_list : it -> ic list
        val copy : it -> it -> unit
        val clone : it -> it
        val median : it -> it -> it -> unit
        val distance : it -> it -> float 
        val median_cost : it -> float
        val compare : it -> it -> int
        val set_state : it -> ic -> unit
        val get_max : it -> int -> int
        val get_min : it -> int -> int
        val get_cost : it -> int -> float 
        val cardinal : it -> int
        val get_state : it -> int -> ic
        val elt_code : ic -> int
        val get_set_code : it -> int
    end

end



module Make (Add : AdditiveInterface) : Make =
  struct

    type ct = Add.ct

    (* Internal only exceptions *)
    exception Success
    exception Failed

    (* Visible exception *)
    exception Exists
    exception Illegal_Arguments
    exception Duplicated
    exception Illegal_State
    exception Not_Found

    type t = {
        codes : (int, Bigarray.int_elt, Bigarray.c_layout) Bigarray.Array1.t;
        characters : ct;
        scode : int;
    }

    type c = (int * int * int)

    let set_code s = s.scode

    let code = set_code

    let char_mem cs t = match cs with
        | None    -> true
        | Some [] -> false
        | Some xs ->
            let x = ref false in
            for i = 0 to (Bigarray.Array1.dim t.codes)-1 do
                x := List.mem t.codes.{i} xs;
            done;
            !x

    let code_mem cs t = match cs with
        | None    -> true
        | Some [] -> false
        | Some xs -> List.mem (code t) xs

    let codes {codes=codes} =
        let len = Bigarray.Array1.dim codes in
        let acc = ref [] in
        for i = len - 1 downto 0 do
            acc := (Bigarray.Array1.get codes i) :: !acc
        done;
        !acc

    let of_array arr k = 
        (* A function to check that there are no duplications *)
        let check_contents =
            let codes = ref All_sets.Integers.empty in
            fun (a, b, c) ->
                if All_sets.Integers.mem c !codes then raise Duplicated
                else begin
                    codes := All_sets.Integers.add c !codes;
                    if b < a then 
                        raise Illegal_State
                    else ()
                end
        in
        (* A comparison function to be able to sort the arrays *)
        let compare (_, _, a) (_, _, b) = a - b in
        (* The arrays to be used in the final construction *)
        let len = Array.length arr in
        (* Store the proper information on each array. *)
        Array.iter check_contents arr;
        Array.sort compare arr;
        let minarr = Array.init len (fun x -> let (a, _, _) = arr.(x) in a)
        and maxarr = Array.init len (fun x -> let (_, a, _) = arr.(x) in a)
        and codes = Array.init len (fun x -> let (_, _, a) = arr.(x) in a) in
        let codes = Bigarray.Array1.of_array Bigarray.int Bigarray.c_layout codes in
        { characters = Add.create minarr maxarr; codes = codes; scode = k }

    let find_pos t it = 
        let len = Bigarray.Array1.dim t.codes in
        let rec finder it pos max =
            if (pos < max) then begin
                if (t.codes.{pos} = it) then pos
                else finder it (pos + 1) max 
            end else failwith "Not_found"
        in
        finder it 0 len

    let median _ a b =
        let c = Add.median a.characters b.characters in
        { characters = c; codes = a.codes; scode = a.scode }

    let reroot_median = median None

    let distance a b = 
        Add.distance a.characters b.characters

    let distance_2 a b c =
        Add.distance_2 a.characters b.characters c.characters

    let dist_2 = distance_2

    let distance_median a b =
        assert (a.scode = b.scode);
        let (x, y) = Add.distance_median a.characters b.characters in
        (x, { characters = y; codes = a.codes; scode = a.scode })

    let median_cost c = 
        Add.median_cost c.characters

    let compare_codes a b =
        let rec comparator a b it max =
            if it < max then begin
                match a.{it} - b.{it} with
                | 0 -> comparator a b (it + 1) max
                | n -> n
            end else 0
        in
        let lena = Bigarray.Array1.dim a.codes in
        let lenb = Bigarray.Array1.dim b.codes in
        match comparator a.codes b.codes 0 (min lena lenb) with
        | 0 -> lena - lenb
        | n -> n

    let compare_data a b =
        Add.compare_data a.characters b.characters

    let get_length a =
        Bigarray.Array1.dim a.codes

    let cardinal = get_length

    let deep_cardinal = cardinal

    let copy a b = (* The two sets must have the same cardinality *)
        if cardinal a = cardinal b then 
            Add.copy a.characters b.characters
        else begin
            print_string "AddCS.copy\n";
            raise Illegal_Arguments
        end

    let clone a =
        let res = Add.clone a.characters 
        and len = Bigarray.Array1.dim a.codes in
        let codes = 
            Bigarray.Array1.create Bigarray.int Bigarray.c_layout len
        in
        for i = 0 to len - 1 do
            codes.{i} <- a.codes.{i};
        done;
        { characters = res; codes = codes; scode = a.scode }

    let set_state a x y c = 
        let resa = clone a in
        let pos = find_pos resa c in
        Add.pos_set_state resa.characters x y pos;
        resa

    let get_max a c =
        Add.pos_get_max a.characters c

    let get_min a c =
        Add.pos_get_min a.characters c

    let get_cost a c =
        Add.pos_get_cost a.characters c

    let get_code a p = a.codes.{p}

    let get_triple t pos = get_min t pos, get_max t pos, get_code t pos

    let get_state = get_triple

    let to_list_with_cost t = 
        let len = get_length t in
        let rec build i acc =
            if i < 0 then
                acc
            else begin
                let r = get_min t i, get_max t i, get_code t i, get_cost t i in
                build (i-1) (r :: acc)
            end
        in
        build (len - 1) []

    let to_list t = 
        let len = get_length t in
        let rec build it acc = 
            if it < 0 then acc
            else begin 
                try let res = get_min t it, get_max t it, get_code t it in
                    build (it - 1) (res :: acc)
                with | Failure "Not_found" -> build (it - 1) acc

            end
        in
        build (len - 1) []

    let of_list l k = 
        let len = List.length l in
        let min = Array.make len 0
        and max = Array.make len 0
        and code = Array.make len 0 in
        let () =
            ignore
                (List.fold_left
                    (fun x (a, b, c) ->
                        min.(x) <- a;
                        max.(x) <- b;
                        code.(x) <- c;
                        x + 1)
                    0 l)
        in
        {
            characters = Add.create min max;
            codes = Bigarray.Array1.of_array Bigarray.int Bigarray.c_layout code;
            scode = k;
        }

    let map f t =
        let lst = to_list t in
        let res = List.map f lst in
        of_list res (set_code t)

    let fold f i t = 
        let lst = to_list t in
        List.fold_left f i lst


    let empty c = 
        { 
            characters = Add.create [||] [||]; 
            codes = Bigarray.Array1.of_array Bigarray.int Bigarray.c_layout [||];
            scode = c;
        }

    let is_empty t = 0 = cardinal t

    let mem (min, max, code) t = 
        (* This function could be log n time but I will do the simplest
        * implementation here *)
        try
            let len = get_length t in
            for i = len - 1 downto 0 do
                if code = get_code t i && min = get_min t i && max = get_max t i then 
                    raise Success
                else ();
            done;
            false
        with
        | Success -> true

    let code_exists code t = 
        try 
            let len = get_length t in
            for i = len - 1 downto 0 do
                if code = get_code t i then raise Success
                else ()
            done;
            false
        with 
        | Success -> true

    let add ((_, _, code) as it) t = 
        if code_exists code t then raise Exists
        else of_list (it :: to_list t) (set_code t)

    let singleton it k = 
        of_list [it] k

    let remove (a, b, c) t = 
        let lst = to_list t in
        let lst = List.filter begin fun (x, y, z) -> 
            not (z == c && a = x && y == b) end lst in
        of_list lst (set_code t)

    let union a b = 
        (* This will be truly slow! *)
        let lst = to_list a in
        List.fold_left (fun x y -> add y (remove y x)) b lst

    let inter a b k =
        let lst = to_list b in
        List.fold_left begin fun x y -> 
            if mem y a then add y x
            else x
        end (empty k) lst

    let diff a b k =
        let lst = to_list a in
        let lst = List.fold_left begin fun x y ->
            if not (mem y b) then y :: x
            else x
        end [] lst in
        of_list (List.rev lst) k

    let equal a b = 0 = compare a b

    let subset a b =
        let lena = cardinal a 
        and lenb = cardinal b in
        if lena <= lenb then begin
            try
                for i = lena - 1 downto 0 do
                    if not (mem (get_triple a i) b) then raise Failed
                    else ()
                done;
                true
            with
            | Failed -> false
        end else false

    let iter f t = 
        let lst = to_list t in
        List.iter f lst

    let fold f t init = 
        let lst = to_list t in
        List.fold_left f init lst

    let for_all f t =
        let lst = to_list t in
        List.fold_left (fun x y -> x && f y) true lst

    let exists f t =
        List.exists f (to_list t)

    let filter f t =
        of_list (List.filter f (to_list t)) (set_code t)

    let f_codes t codes =
        let check (_, _, c) = All_sets.Integers.mem c codes in
        filter check t

    let f_codes_comp t codes =
        let check (_, _, c) = not (All_sets.Integers.mem c codes) in
        filter check t

    let partition f a =
        let lst = to_list a in
        let t, f = List.fold_left 
            (fun (tr, fl) x -> if f x then (x :: tr, fl) else (tr, x :: fl)) ([], []) lst
        in
        of_list t (set_code a), of_list f (set_code a)

    let elements = to_list

    let min_elt x = 
        if 0 = cardinal x then raise Not_Found
        else get_triple x 0

    let max_elt x =
        let len = cardinal x in
        if len != 0 then get_triple x (len - 1) 
        else raise Not_Found

    let choose t =
        min_elt t

    let split it t =
        let comparer (less, is, greater) a =
            let c = Pervasives.compare a it in
            if c = 0 then (less, true, greater)
            else if c < 0 then (a :: less, is, greater)
            else (less, is, a :: greater)
        in
        let less, b, greater = 
            List.fold_left comparer ([], false, []) (to_list t) 
        in
        of_list less (set_code t), b , of_list greater (set_code t)

    let elt_code (_, _, c) = c

    let median_3 p n c1 c2 = 
        let res = 
            Add.median_3 p.characters n.characters c1.characters c2.characters 
        in
        { n with characters = res }

    let full_union a b c =
        Add.full_unioni a.characters b.characters c.characters

    let state_to_xml ch c (_ : Data.d) =
        let print_character (min, max, code, cost) =
            let beg = "<characterAdditive code=\"" ^ string_of_int code ^ 
            "\" cost=\"" ^ string_of_float cost ^ "\">\n" in
            output_string ch beg;
            output_string ch ("<min>" ^ string_of_int min ^ "</min>\n");
            output_string ch ("<max>" ^ string_of_int max ^ "</max>\n");
            output_string ch "</characterAdditive>\n"
        in
        let c = to_list_with_cost c in
        List.iter print_character c

    let to_formatter attr c parent d : Xml.xml Sexpr.t list =
        let module T = Xml.Characters in
        let c_ls = to_list c in
        let c_parent_ls = match parent with
            | Some parent -> to_list parent
            | None -> c_ls
        in
        let idx = ref 0 in
        let output_character (min, max, code, cost) =
            let cost = distance (singleton (List.nth c_ls !idx) 0)
                                (singleton (List.nth c_parent_ls !idx) 0)
            in
            incr idx;
            (PXML
                -[T.additive]
                    (* Attributes *)
                    ([T.name] = [`String (Data.code_character code d)])
                    ([T.cost] = [`Float cost])
                    ([T.definite] = [`Bool (cost > 0.0)])
                    ([attr])

                    (*Contents *)
                    -[T.min] { `String (Data.to_human_readable d code min) } --
                    -[T.max] { `String (Data.to_human_readable d code max) } --
                --)
        in
        let c = to_list_with_cost c in
        List.map output_character c


    (** Now a purely imperative version of this library *)
    module Imperative =
      struct
        type it = t
        type ic = { mutable min : int; mutable max : int; mutable code : int }

        let of_array arr k =
            let len = Array.length arr in
            let min = Array.init len (fun i -> arr.(i).min)
            and max = Array.init len (fun i -> arr.(i).max)
            and code = Array.init len (fun i -> arr.(i).code) in
            {   characters = Add.create min max; 
                codes = Bigarray.Array1.of_array Bigarray.int Bigarray.c_layout code;
                scode = k }
        
        let create a = clone a

        let of_list lst k =
            of_array (Array.of_list lst) k

        let to_list r = 
            let create (a, b, c) =
                { 
                    min = a;
                    max = b;
                    code = c;
                }
            in
            List.map create (to_list r)

        let copy = copy

        let clone = clone

        let median a b c =
            Add.mediani a.characters b.characters c.characters

        let distance = distance

        let median_cost = median_cost

        let compare = compare

        let set_state a k =
            let pos = find_pos a k.code in
            Add.pos_set_state a.characters k.min k.max pos

        let get_max = get_max

        let get_min = get_min

        let get_cost = get_cost

        let cardinal = cardinal

        let get_state a b =
            let (a, b ,c) = get_state a b in
            { min = a; max = b; code = c }

        let elt_code x = x.code

        let get_set_code = set_code

        let median_3 p n c1 c2 = Add.median_3 p n c1 c2

    end


    let of_parser data (it, taxon) code =
        let first lst = match lst with
            | h::_ -> h
            | _    -> assert false
        in
        let rec last lst = match lst with
            | [h]   -> h
            | _::tl -> last tl
            | []    -> assert false
        in
        let check_type_and_val acc = function
            | Some v, code ->
                let v = match v with
                    | `List x -> x
                    | `Bits x -> BitSet.to_list x
                in
                let v = List.sort compare v in
                (first v, last v, code) :: acc
            | None, code ->
                begin match Hashtbl.find data.Data.character_specs code with
                    | Data.Static (Data.NexusFile enc) ->
                            (first enc.Nexus.File.st_observed,
                                    last enc.Nexus.File.st_observed,code) :: acc
                    | _ -> assert false
                end
        in
        let arr = 
            Array.of_list (List.rev (Array.fold_left check_type_and_val [] it))
        in
        of_array arr code, taxon

    let to_string a = 
        Add.to_string a.characters

end

(* Battery of tests for the library *)
module Test = struct
 
    type ct = AddVec.ct * AddGen.ct

    let failwithf format = Printf.ksprintf (failwith) format

    let compare_vector_general vt gt =
        if (AddGen.cardinal gt) = (AddVec.cardinal vt) then begin
            let same = ref true in
            for i = 0 to (AddGen.cardinal gt)-1 do
                let this_char =
                    let max_match = 
                        if (AddVec.pos_get_max vt i) = 255 then
                            let max_gen = AddGen.pos_get_max gt i in
                            (max_gen = 255) || (max_gen = max_int)
                        else
                            (AddVec.pos_get_max vt i) = (AddGen.pos_get_max gt i)
                    in
                    let min_match = (AddGen.pos_get_min gt i) = (AddVec.pos_get_min vt i)
                    and cst_match = (AddGen.pos_get_cost gt i) = (AddVec.pos_get_cost vt i) in
                    max_match && min_match && cst_match
                in
                same := !same && this_char;
            done;
            !same
        end else begin
            false
        end

    let create x y : ct =
        assert(Array.fold_left (fun acc x -> acc && (x <= 255) && (x >= 0)) true x);
        assert(Array.fold_left (fun acc y -> acc && (y <= 255) && (y >= 0)) true y);
        let v = AddVec.create x y
        and g = AddGen.create x y in
        assert( compare_vector_general v g );
        (v,g)

    let copy ((x1,x2):ct) ((y1,y2):ct) : unit =
        AddGen.copy x2 y2;
        AddVec.copy x1 y1;
        assert( compare_vector_general y1 y2 );
        ()

    let clone (xv,xg) =
        let g = AddGen.clone xg in
        let v = AddVec.clone xv in
        assert( compare_vector_general v g );
        (v,g)

    let compare_data (xv,xg) (yv,yg) =
        let g_comp = AddGen.compare_data xg yg in
        let v_comp = AddVec.compare_data xv yv in
        assert( compare_vector_general xv xg );
        assert( compare_vector_general yv yg );
        assert( g_comp = v_comp );
        g_comp

    let median (xv,xg) (yv,yg) =
        let v = AddVec.median xv yv
        and g = AddGen.median xg yg in
        assert( compare_vector_general v g );
        (v,g)

    let distance (xv,xg) (yv,yg) =
        let dg = AddGen.distance xg yg 
        and dv = AddVec.distance xv yv in
        assert( dg = dv );
        dg

    let distance_2 (xv,xg) (yv,yg) (zv,zg) =
        let dv = AddVec.distance_2 xv yv zv
        and dg = AddGen.distance_2 xg yg zg in
        assert( dg = dv );
        dg

    let distance_median (xv,xg) (yv,yg) =
        let vd,vm = AddVec.distance_median xv yv
        and gd,gm = AddGen.distance_median xg yg in
        assert( compare_vector_general vm gm );
        assert( vd = gd );
        gd,(vm,gm)

    let median_cost (xv,xg) =
        let v = AddVec.median_cost xv
        and g = AddGen.median_cost xg in
        assert( v = g );
        v

    let median_3 (wv,wg) (xv,xg) (yv,yg) (zv,zg) =
        let dg = AddGen.median_3 wg xg yg zg
        and dv = AddVec.median_3 wv xv yv zv in
        assert (compare_vector_general dv dg);
        assert( compare_vector_general dv dg );
        (dv,dg)

    let full_unioni (xv,xg) (yv,yg) (zv,zg) =
        let () = AddVec.full_unioni xv yv zv
        and () = AddGen.full_unioni xg yg zg in
        assert( compare_vector_general zv zg);
        ()

    let full_union (xv,xg) (yv,yg) =
        let uv = AddVec.full_union xv yv
        and ug = AddGen.full_union xg yg in
        assert( compare_vector_general uv ug);
        (uv,ug)
   
    let mediani (xv,xg) (yv,yg) (zv,zg) =
        let () = AddVec.mediani xv yv zv
        and () = AddGen.mediani xg yg zg in
        assert( compare_vector_general zv zg);
        ()

    let pos_get_max (v,g) i =
        let m = AddVec.pos_get_max v i in
        assert(m = (AddGen.pos_get_max g i));
        m

    let pos_get_min (v,g) i =
        let m = AddVec.pos_get_min v i in
        assert(m = (AddGen.pos_get_min g i));
        m

    let pos_get_cost (v,g) i =
        let m = AddVec.pos_get_cost v i in
        assert(m = (AddGen.pos_get_cost g i));
        m

    let cardinal (v,g) =
        let x = AddVec.cardinal v in
        assert (x = AddGen.cardinal g);
        x

    let pos_set_state ((v,g) as t) i smin smax =
        let () = AddVec.pos_set_state v i smin smax in
        let () = AddGen.pos_set_state g i smin smax in
        ignore (pos_get_min t i);
        ignore (pos_get_max t i);
        ()

    let print _ = ()
        let to_string _ = ""

end 

module Vector  = Make (AddVec)
(*module Vector  = Make (AddGen)*)
module General = Make (AddGen)


(* We assume the characters are all static continuous/additive characters. We
   want to find the subset that is vectorizable; range is < sizeof(char). This
   is previously worked out in st_normal. If used, we know that the characters
   have been normalized to fit into the alloted vectorized space. *)
let split_vectorized_characters data codes =
    let is_character_vectorizable code = 
        let spec = match Hashtbl.find data.Data.character_specs code with
            | Data.Static (Data.NexusFile spec) -> spec
            | Data.Static _ | Data.Dynamic _ | Data.Set | Data.Kolmogorov _ -> assert false
        in
        match spec.Nexus.File.st_normal with
        | Some _ -> true
        | None   -> false
    in
    match codes with
    | [] -> [],[]
    | _  ->
        let vectorized,general = List.partition is_character_vectorizable codes in
        Printf.ksprintf (Status.user_message Status.Information)
                        ("@[Vectorized:@[%d of %d characters.@]@]")
                        (List.length vectorized) (List.length codes);
        vectorized,general


let is_potentially_informative elts = 
    let intersection a b = match a with
        | None -> None
        | Some (x, y) -> match b with
            | None  -> a
            | Some lst ->
                match Nexus.File.static_state_to_list lst with
                | [] -> a
                | lst ->
                    let b = List.fold_left min max_int lst
                    and c = List.fold_left max 0 lst in
                    if (b <= x && x <= c) then
                        Some (x, (min y c))
                    else if (x <= b && b <= y) then
                        Some (b, (min y c))
                    else 
                        None
    in
    match List.fold_left intersection (Some (0, max_int)) elts with
    | None -> true
    | Some _ -> false

let min_possible_cost (elts : Nexus.File.static_state list ) =
    let get_last lst = 
        assert (lst <> []);
        List.hd (List.rev lst) 
    in
    let elts = NonaddCS8.extract_elements_present elts in
    let elts = List.map (List.sort ( - )) elts in
    let codes = List.fold_left (fun acc lst ->
        try
            let fst, lst = match lst with
                | [h]  -> h, h
                | h::t -> h, get_last t
                | []   -> raise Exit
            in
            acc --> All_sets.Integers.add fst --> All_sets.Integers.add lst
        with
        | Exit -> acc) All_sets.Integers.empty elts
    in
    let elts = 
        List.map 
            (function [x] -> All_sets.Integers.singleton x
                | [] -> All_sets.Integers.empty
                | h :: t ->
                    let t = get_last t in
                    All_sets.Integers.fold
                        (fun a acc ->
                            if a <= t && a >= h then 
                                All_sets.Integers.add a acc
                            else acc)
                        codes All_sets.Integers.empty)
            elts
    in
    let codes = List.sort ( - ) (All_sets.Integers.elements codes) in
    let filter c left =
        (List.filter (fun x ->not (All_sets.Integers.mem c x)) left)
    in
    let rec optimal_cost codes left = match codes, left with
        | _, [] -> None, 0
        | [], _ -> 
            (* We assign a big value that won't wrap
            * to the negatives if we add something to it *)
            None, (max_int / 2) 
        | (h :: t), left ->
            let left' = filter h left in
            let leftl', leftc' = optimal_cost t left'
            and (_, leftc) as second = optimal_cost t left in
            let (_, leftc') as first =
                match leftl' with
                | None -> (Some h), leftc'
                | Some x -> (Some h), leftc' + (x - h)
            in
            (* We have to prefer the solution that contains the smallest, so
            * in case of equality, we choose leftc' *)
            if leftc' <= leftc then first
            else second
    in
    let _, cost = optimal_cost codes elts in
    float_of_int cost

