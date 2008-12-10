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

(* $Id: charset.ml 1644 2007-02-14 19:05:47Z andres $ *)
let () = SadmanOutput.register "Charset" "$Revision: 1644 $"


type parscheme =
    | ParCodewise
    | Nopar

module IntSet = All_sets.Integers

exception Unmatched_code of int

let _charset_num = ref 0 
let get_par_id () =
    incr _charset_num;
    !_charset_num

module Of(Elt : Character.CHARACTER) = struct

    (* Type of elements in our set *)
    type e = Elt.t

    (* Type of elements in the actual Set we use *)
    type mye = {elt : Elt.t;
                cost : float;
                code : int;
               }

    module CharSetComp = struct
        type t = mye
        let compare i j =
            compare i.code j.code
    end
    module CharSet = Set.Make(CharSetComp)
    module Elt = Elt

    (* Replace a value in a char set *)
    let set_replace elt set =
        CharSet.add elt (CharSet.remove elt set)
        
    (* This is our type *)
    type t = {set : CharSet.t;
              heur : Character.h;
              changeset : IntSet.t;
              parstyle : parscheme;

              (* We need this to keep track of changes *)
              from1 : t option;
              from2 : t option;
              from3 : t option;
              parent : t option;
              set_code : int;
             }


    let color = Character.Blue
    let mycode = Character.new_code ()
    let code {set_code=c} = c
    let set_code t code = {t with set_code = code}
    let elt_code e = Elt.code e

    (* We need a default heuristic *)
    let empty = {
        set=CharSet.empty;
        heur=Character.Exact;
        changeset=IntSet.empty;
        parstyle=Nopar;

        from1=None;
        from2=None;
        from3=None;
        parent=None;
        set_code=0;
    }

    let make_empty h c p () =
        {empty with
             heur = h;
             set_code = c;
             parstyle = p;
        }
            
    let make_empty_like s =
        {s with
             set = CharSet.empty;
             changeset = IntSet.empty;
             from1 = None;
             from2 = None;
             from3 = None;
             parent = None;
        }

    let add (elt : e) cost ({set=set} as t) =
        {t with set=CharSet.add {elt = elt; cost = cost; code = elt_code elt} set}
    let del delcode ({set=set} as t) =
        {t with set=CharSet.filter (fun {code=code} -> code <> delcode) set}
            
    let cardinal {set=set} = CharSet.cardinal set
    let deep_cardinal {set=set} =
        CharSet.fold (fun e n -> Elt.deep_cardinal e.elt + n) set 0

    (* Type to make new random character sets *)
(*     type gen = int * Elt.gen * Character.h * parscheme * int *)

(*     let rand_gen () = *)
(*         let maxelts = 20 in *)
(*         let minelts = 10 in *)
(*         let elt_gen = Elt.rand_gen () in *)
(*         ((Random.int (maxelts - minelts)) + minelts, elt_gen, *)
(*          Character.Exact, Nopar, Character.new_code ()) *)

(*     let rand_gen_parallel () = *)
(*         let maxelts = 20 in *)
(*         let minelts = 10 in *)
(*         let elt_gen = Elt.rand_gen () in *)
(*         ((Random.int (maxelts - minelts)) + minelts, elt_gen, *)
(*          Character.Exact, ParCodewise, Character.new_code ()) *)

(*     let make_rand (nelts, gen, h, p, c) = *)
(*         let set = ref (make_empty h c p ()) in *)
(*         for i = 0 to nelts do *)
(*             set := add (Elt.make_rand gen) 1. !set *)
(*         done; *)
(*         !set *)
        
    (* Note: colors are static to a module;  so this function doesn't really
       make sense, the way things are currently implemented *)
(*     let colors {set=set} = *)
(*         List.map (fun {elt=elt; code=code; cost=cost} -> (code, Elt.color)) *)
(*             (CharSet.elements set) *)
    let codes {set=set} =
        List.map (fun {elt=elt; code=code; cost=cost} -> code)
            (CharSet.elements set)
    let costs {set=set} =
        List.map (fun {elt=elt; code=code; cost=cost} -> (code, cost))
            (CharSet.elements set)
    let median_cost t =
        List.fold_left (fun s (_, f) -> s +. f) 0. (costs t)

    let merge {set=set1} ({set=set2} as t) =
        {t with set=CharSet.union set1 set2}
    let minus {set=set1} ({set=set2} as t) =
        {t with set=CharSet.diff set1 set2}
    let random rand ({set=set} as t) =
        {t with set=CharSet.filter (fun _ -> rand ()) set}

    (* a folding function useful for us, internally *)
    let myfold f start {set=set} =
        CharSet.fold f set start

    (* helper fn: return the element with a given code, if it exists *)
    let get_elt_withcode code set =
        myfold (fun {elt=nelt; code=ncode} elt -> if code=ncode then Some nelt else elt)
            None set
    let get_eltrec_withcode code set =
        myfold (fun ({elt=nelt; code=ncode} as erec) elt -> if code=ncode then Some erec else elt)
            None set
    (* TODO: get rid of this; elt lists are ordered *)
    let get_code_pairs :
            t -> t -> (mye * mye) list * mye list * mye list =
        let insert set1 set2 (paired, list1, list2) code =
            let elts = (get_eltrec_withcode code set1,
                        get_eltrec_withcode code set2) in
            match elts with
            | (Some elt1, Some elt2) ->
                  ((elt1, elt2) :: paired, list1, list2)
            | (Some elt1, None) ->
                  (paired, elt1 :: list1, list2)
            | (None, Some elt2) ->
                  (paired, list1, elt2 :: list2)
            | (None, None) ->
                  assert ("This should never happen!  Where'd we get the code?"
                          <> "");
                  (paired, list1, list2)
        in
        fun set1 set2 ->
            let codes = List.merge compare (codes set1) (codes set2) in
            List.fold_left (insert set1 set2) ([],[],[]) codes
            

    let update code elt set =
        match get_eltrec_withcode code set with
        | None ->                (* This shouldn't actually happen.  But, OK. *)
              {set with
                   (* TODO: default cost *)
                   set = CharSet.add {elt=elt;cost=0.;code=code} set.set;
                   changeset = IntSet.add code set.changeset;
              }
        | Some oldrec ->
              if (Elt.compare_data elt oldrec.elt) = 0
              then                      (* Do we want to replace it anyway? *)
                  set
              else
                  {set with
                       set = set_replace {oldrec with elt = elt} set.set;
                       changeset = IntSet.add code set.changeset;
                  }

    let substitute ({set=set1} as a) (t2) =
        let newset = CharSet.fold
            (fun ({code=code} as orig_elt) newset ->
                 match get_eltrec_withcode code t2 with
                 | None -> CharSet.add orig_elt newset
                 | Some new_elt -> CharSet.add new_elt newset)
            set1 CharSet.empty in
        {a with set=newset}
        

    let to_list {set=set} =
        List.map (fun {elt=elt; code=code; cost=cost} -> (code, elt, cost))
            (CharSet.elements set)
    let of_list list =
        List.fold_left
            (fun set (code, elt, cost) -> add elt cost set)
            empty
            list

    let set_heu h t = {t with heur = h}
    let get_heu {heur = h} = h

    let set_par p t = {t with parstyle = p}
    let get_par {parstyle = p} = p

    (* utility functions *)
    let set_changeset set t = {t with changeset=set}
    let get_changeset {changeset=set} = set

    let fold f =
        myfold (fun {elt=elt} -> f elt)
    let filter f ({set=set} as t) =
        {t with set=CharSet.filter f set}
    let f_codes t codes =
        filter (fun {code=code} -> List.mem code codes) t
(*     let f_colors t color = *)
(*         filter (fun {elt=elt} -> color = Elt.color) t *)
    let iter f {set=set} =
        CharSet.iter (fun {elt=elt; code=code; cost=cost} -> f elt code) set
    let map f ({set=set} as t) =
        {t with set =
                CharSet.fold
                    (fun e s ->
                         CharSet.add {e with elt=f e.elt e.code} s)
                    set CharSet.empty}

    let is_empty {set=set} = CharSet.is_empty set

    let derived_copy set from1 from2 from3 =
        {set with
             changeset = IntSet.empty;
             from1 = from1;
             from2 = from2;
             from3 = from3;
             parent = Some set;
        }

    (* no real parsing yet *)
    let parse _ = []
    let to_string t =
        "(" ^
            (myfold (fun {elt=elt; code=code; cost=cost} str ->
                         (Elt.to_string elt) ^
                             if str = ""
                             then ")"
                             else " " ^ str)
                 "" t)

    let to_formatter attr c d : Tags.xml list =
        []                              (* TODO: output *)

    let par_id = get_par_id ()
    let assert_par_id pid =
        assert (pid = par_id)
    let marshal (s : t) = (par_id, Marshal.to_string s [])
    let unmarshal (pid, s) =
        assert (pid = par_id);
        ((Marshal.from_string s 0) : t)

    let make_parfunc fn =
        (par_id, Marshal.to_string fn [Marshal.Closures])
    let get_parfunc (pid, mfn) =
        assert (pid = par_id);
        (Marshal.from_string mfn 0)

    let dispatch_gather parid ms1 ms2 code fntag =
        assert (parid = par_id);
        let s1 = unmarshal ms1 in
        let s2 = unmarshal ms2 in
        let fn : int -> Elt.t -> Elt.t -> 'a
            = (get_parfunc fntag) in
        let result =
            match
                get_elt_withcode code s1, get_elt_withcode code s2
            with
            | Some e1, Some e2 ->
                  fn code e1 e2
            | _ -> raise (Unmatched_code code)
        in
        result

    let pair_gather fn s1 s2 icodes parscheme init incfn =
        let icodes = if icodes = [] then codes s1 else icodes in
(*         match parscheme with *)
(*         | ParCodewise -> *)
(*               let queue = Asynch.queue_list icodes in *)
(*               let tag = Register.set_pair_sum_tag in *)
(*               let fntag = make_parfunc fn in *)

(*               Asynch.run_jobs_fold *)
(*                   (ref queue) tag *)
(*                   (\* Message that we send *\) *)
(*                   (fun code -> *)
(*                        let subset1 = f_codes s1 [code] in *)
(*                        let subset2 = f_codes s2 [code] in *)
(*                        (make_parfunc dispatch_gather, *)
(*                         marshal subset1, *)
(*                         marshal subset2, *)
(*                         code, *)
(*                         fntag)) *)
(*                   init *)
(*                   (\* Received message should be (code, res) *\) *)
(*                   (fun (code, result) acc -> incfn result acc) *)
(*         | Nopar -> *)
              List.fold_left
                  (fun acc code ->
                       let result =
                           match
                               get_elt_withcode code s1, get_elt_withcode code s2
                           with
                           | Some e1, Some e2 ->
                                 fn code e1 e2
                           | _ -> raise (Unmatched_code code)
                       in
                       (incfn result acc))
                  init icodes
    let tri_gather fn s1 s2 s3 icodes init incfn =
        let icodes = if icodes = [] then codes s1 else icodes in
        List.fold_left
            (fun acc code ->
                 let result =
                     match
                         get_elt_withcode code s1,
                         get_elt_withcode code s2,
                         get_elt_withcode code s3
                     with
                     | Some e1, Some e2, Some e3 ->
                           fn code e1 e2 e3
                     | _ -> raise (Unmatched_code code)
                 in
                 (incfn result acc))
            init icodes

    (* This is a "maybe-do" operator, similar to one of Haskell's monads.  If
       the first argument is None, we return None;  if it is Some a, we apply f
       to a (which must return some sort of option type itself).  This idiom for
       chaining together commands comes from Xavier Leroy (see
       http://groups.google.com/group/fa.caml/msg/f194a3240d240e71 ) *)
    let (++?) x f =
        match x with
        | None -> None
        | Some v -> f v
    let (++?>) x f =
        match x with
        | None -> None
        | Some v -> Some (f v)

    let dispatch_fold parid mp ms1 ms2 code fntag =
        assert (parid = par_id);
        let parset = mp ++?> unmarshal in
        let s1 = unmarshal ms1 in
        let s2 = unmarshal ms2 in
        let fn : int -> Elt.t option -> Elt.t -> Elt.t -> Elt.t
            = (get_parfunc fntag) in
        let result =
            match
                (parset ++? (get_elt_withcode code),
                 get_elt_withcode code s1,
                 get_elt_withcode code s2)
            with
            | parelt, Some e1, Some e2 ->
                  fn code parelt e1 e2
            | _ -> raise (Unmatched_code code)
        in
        result

    let pair_fold fn parset s1 s2 icodes parscheme into =
        let icodes = if icodes = [] then codes s1 else icodes in
(*         match parscheme with *)
(*         | ParCodewise -> *)
(*               let queue = Asynch.queue_list icodes in *)
(*               let tag = Register.set_pair_fold_tag in *)
(*               let fntag = make_parfunc fn in *)

(*               Asynch.run_jobs_fold *)
(*                   (ref queue) tag *)
(*                   (\* Message that we send *\) *)
(*                   (fun code -> *)
(*                        let subset1 = f_codes s1 [code] in *)
(*                        let subset2 = f_codes s2 [code] in *)
(*                        (make_parfunc dispatch_fold, *)
(*                         parset *)
(*                         ++?> (fun set -> f_codes set [code]) *)
(*                         ++?> marshal, *)
(*                         marshal subset1, *)
(*                         marshal subset2, *)
(*                         code, *)
(*                         fntag)) *)
(*                   into *)
(*                   (\* Received message should be (code, elt) *\) *)
(*                   (fun (code, result) acc -> *)
(*                        update code result acc) *)
(*         | Nopar -> *)
        List.fold_left
            (fun acc code ->
                 let result =
                     match
                         (parset ++? (get_elt_withcode code),
                          get_elt_withcode code s1,
                          get_elt_withcode code s2)
                     with
                     | parelt, Some e1, Some e2 ->
                           fn code parelt e1 e2
                     | _ -> raise (Unmatched_code code)
                 in
                 update code result acc)
            into icodes
    let tri_fold fn parset s1 s2 s3 icodes into =
        let icodes = if icodes = [] then codes s1 else icodes in
        List.fold_left
            (fun acc code ->
                 let result =
                     match
                         (parset ++? (get_elt_withcode code),
                          get_elt_withcode code s1,
                          get_elt_withcode code s2,
                          get_elt_withcode code s3)
                     with
                     | parelt, Some e1, Some e2, Some e3 ->
                           fn code parelt e1 e2 e3
                     | _ -> raise (Unmatched_code code)
                 in
                 update code result acc)
            into icodes
    let quad_fold fn parset s1 s2 s3 s4 icodes into =
        let icodes = if icodes = [] then codes s1 else icodes in
        List.fold_left
            (fun acc code ->
                 let result =
                     match
                         (parset ++? (get_elt_withcode code),
                          get_elt_withcode code s1,
                          get_elt_withcode code s2,
                          get_elt_withcode code s3,
                          get_elt_withcode code s4)
                     with
                     | parelt, Some e1, Some e2, Some e3, Some e4 ->
                           fn code parelt e1 e2 e3 e4
                     | _ -> raise (Unmatched_code code)
                 in
                 update code result acc)
            into icodes
                  
    let dist_help _ = Elt.distance
    let distance t1 t2 =
        pair_gather dist_help t1 t2 [] t1.parstyle 0. (+.)

    let dist_list_help code a b = (code, Elt.distance a b)
    let distance_list t1 t2 =
        pair_gather
            dist_list_help
            t1 t2 [] t1.parstyle [] (fun a b -> a :: b)

    let median_help _ p a b =
        Elt.median p a b
    let median prev s1 s2 =
        let codes =
            match prev with
            | None -> []
            | Some old ->
                   (* Find what's changed in i and j (check if they're an
                      updated version of what we saw earlier) *)
                  let get_changelist oldt newt =
                      if oldt == newt.parent       (* We're functional!
                                                      This works. *)
                      then IntSet.elements (newt.changeset)
                      else codes newt in
                  let (oldi, oldj) = old.from1, old.from2 in
                  let changelist_i = get_changelist oldi s1 in
                  let changelist_j = get_changelist oldj s2 in
                  let changelist =
                      List.merge compare changelist_i changelist_j in
                  changelist
        in
        let into =
            match prev with
            | Some p -> derived_copy p (Some s1) (Some s2) None
            | None -> make_empty_like s1 in
        pair_fold median_help prev s1 s2 codes s1.parstyle into

    let reroot_median s1 s2 =
        let into = make_empty_like s1 in
        pair_fold (fun _ _ -> Elt.reroot_median) None s1 s2 [] () into

    let median_3_help _ _ a b c d =
        Elt.median_3 a b c d
    let median_3 s1 s2 s3 s4 =
        let into = make_empty_like s1 in
        quad_fold median_3_help None s1 s2 s3 s4 [] into

    let dist_2 s1 s2 s3 =
        tri_gather (fun _ -> Elt.dist_2) s1 s2 s3 [] (0.) ( +. )
                  
    let rec compare_list compare l1 l2 =
        match (l1, l2) with
        | ([], []) -> 0
        | (x :: xs, y :: ys) ->
              let cmp = compare x y in
              if cmp <> 0
              then cmp
              else compare_list compare xs ys
        | ([], _) -> -1
        | (_, []) -> 1

    let compare_codes {set=set1} {set=set2} = CharSet.compare set1 set2
    let compare_data {set=set1} {set=set2} =
        let cmp = CharSet.compare set1 set2 in
        if cmp <> 0
        then cmp
        else compare_list
            (fun a b -> Elt.compare_data a.elt b.elt)
            (CharSet.elements set1)
            (CharSet.elements set2)

                  
end
