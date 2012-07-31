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

let () = SadmanOutput.register "Alphabet" "$Revision: 2871 $"


(* $Id: alphabet.ml 2871 2008-05-23 17:48:34Z andres $ *)

exception Illegal_Character of string
exception Illegal_Code of int
exception Illegal_List of int list

let debug = false
(*
let uselevel = false (* Don't forget to set "uselevel" inside cost_matrix.ml and
sequence.ml*) 
*)
type kind = 
    | Sequential
    | Simple_Bit_Flags
    | Extended_Bit_Flags
    | Continuous
    | Combination_By_Level

type a = {
    comb_to_list : int list All_sets.IntegerMap.t;
    list_to_comb : int All_sets.IntegerListMap.t;
    string_to_code : int All_sets.StringMap.t;
    code_to_string : string All_sets.IntegerMap.t;
    complement : int option All_sets.IntegerMap.t;
    level : int;
    ori_size : int;
    gap : int;
    all : int option;
    size : int;
    kind : kind;
    orientation : bool;
    threeD : bool;
}

let is_combination_by_level alpha =
    match alpha.kind with
    | Combination_By_Level -> true 
    | _ -> false

let continuous =
    {
        comb_to_list   = All_sets.IntegerMap.empty;
        list_to_comb   = All_sets.IntegerListMap.empty;
        string_to_code = All_sets.StringMap.empty;
        code_to_string = All_sets.IntegerMap.empty;
        complement     = All_sets.IntegerMap.empty;
        orientation    = false;
        level    = 0;
        ori_size = 0;
        gap      = -1;
        all      = None;
        size     = 0;
        kind     = Continuous;
        threeD   = false;
    }

let print alpha = 
    Printf.printf "Alphabet.print alpha, ";
    let () = match alpha.kind with
        | Simple_Bit_Flags -> Printf.printf "alph kind:Simple_Bit_Flags\n%!"
        | Sequential -> Printf.printf "alph kind:Sequential\n%!"
        | Extended_Bit_Flags -> Printf.printf "alph kind:Extended_Bit_Flags\n%!"
        | Continuous -> Printf.printf "alph kind:Continuous\n%!"
        | Combination_By_Level -> Printf.printf "alph kind:Combination_By_Level\n%!"
    in
    Printf.printf
        "level=%d,ori_size=%d,size=%d,orientation=%b,threeD=%b, gap=%d,%!"
        alpha.level alpha.ori_size  alpha.size alpha.orientation alpha.threeD alpha.gap;
    let () = match alpha.all with
        | None   -> Printf.printf "no all elements\n%!"
        | Some x -> Printf.printf "all elements=%d\n%!" x
    in
    All_sets.IntegerMap.iter 
        (fun code char -> Printf.fprintf stdout "%i %s\n" code char)
        alpha.code_to_string;
    Printf.printf "end of print alph\n%!"

let to_formatter alph : Xml.xml =
    let element_to_tags string code acc =
        (`Single (Xml.Alphabet.element, [(Xml.Alphabet.value, `String string);
        (Xml.Alphabet.code, `Int code)], `Empty)) :: acc
    in
    let res = 
        `Set (All_sets.StringMap.fold element_to_tags alph.string_to_code [])
    in
    Xml.Characters.alphabet, [], res


(* The alphabet type *)

(* Some constant alphabets  *)

(* Nucleic Acids. Each is assigned a different bit in an integer for better
* representation of different combinations of bases (e.g. adenine or citosine =
* 3 = 1 + 2 = adenine lor citosine). *)
let adenine = 1
let citosine = 2
let guanine = 4
let timine = 8 
let gap = 16
let gap_repr = "-"
let elt_complement = "~"
let uracile = timine

(* Amino Acids. Each is assigned a unique number *)
let alanine = 1        
let arginine = 2       
let asparagine = 3     
let aspartic = 4       
let cysteine = 5       
let glutamine = 6      
let glutamic = 7       
let glycine = 8        
let histidine = 9      
let isoleucine = 10    
let leucine = 11       
let lysine = 12        
let methionine = 13    
let phenylalanine = 14 
let proline = 15       
let serine = 16        
let threonine = 17     
let tryptophan = 18    
let tyrosine = 19      
let valine = 20        
let unspecified = 21   
let all_aminoacids = 21
let aa_gap = 22        


let check_level alph =
    let level = alph.level in
    let size = alph.ori_size in 
    if (level>1)&&(level<=size) then true
    else false

(*[list_to_a] init alphabet with input lst, note that level is set to 0, size
* and ori_size are set to length of input list. *)
let list_to_a ?(respect_case = false) ?(orientation=false) ?(init3D=false) lst gap all kind =
    let a_size = List.length lst in
    let debug = false in
    if debug then 
        Printf.printf "Alphabet.list_to_a, sz=%d,init3D=%b, case sensitive=%b\n%!"
                        a_size init3D respect_case;
    let add (s2c, c2s, cmp, cnt) (a, b, c) =
        if debug then Printf.printf "add %s,%d to s2c; %!" a b;
        let uppa = if respect_case then (String.uppercase a) else a in
        All_sets.StringMap.add uppa b s2c,
        (if All_sets.IntegerMap.mem b c2s then c2s else All_sets.IntegerMap.add b a c2s),
        (All_sets.IntegerMap.add b c cmp),
        cnt + 1
    in
    let empty = 
        All_sets.StringMap.empty, 
        All_sets.IntegerMap.empty,
        All_sets.IntegerMap.empty, 0 
    in
    let s2c, c2s, cmp, cnt = List.fold_left add empty lst in
    let gap_code = 
        let gap = if respect_case then String.uppercase gap else gap in
        try All_sets.StringMap.find gap s2c with
        | Not_found as err ->
            List.iter (fun (x, _, _) -> Status.user_message Status.Error x) lst;
            Status.user_message Status.Error ("could not find the gap " ^ gap);
            raise err
    and all_code = match all with
        | Some all ->
            let all = 
                All_sets.StringMap.find 
                    (if respect_case then String.uppercase all else all)
                    s2c
            in
            Some all
        | None -> None
    in
    let (comb_to_list,list_to_comb) =
        All_sets.StringMap.fold
            (fun _ (v:int) (c2l,l2c) ->
                let vec = match kind with
                    | Sequential           -> [v]
                    | Simple_Bit_Flags     -> BitSet.Int.list_of_packed v
                    | Extended_Bit_Flags   -> BitSet.Int.list_of_packed v
                    | Continuous           -> assert false (* singleton *)
                    | Combination_By_Level -> assert false (* requires other function *)
                in
                All_sets.IntegerMap.add v vec c2l,
                All_sets.IntegerListMap.add vec v l2c)
            (s2c)
            (All_sets.IntegerMap.empty,All_sets.IntegerListMap.empty)
    in
    if debug then Printf.printf "end of list_to_a\n%!";
    { comb_to_list = comb_to_list;
      list_to_comb = list_to_comb;
      level = 0; 
      string_to_code = s2c;
      code_to_string = c2s;
      gap = gap_code;
      all = all_code;
      size = a_size;
      ori_size = a_size;
      kind = kind;
      complement = cmp;
      orientation = orientation;
      threeD = init3D; }

(* used to calculate costs of gaps in static characters / implied alignments *)
let present_absent =
    list_to_a [("present", 1, None); ("absent", 2, None) ]
              "absent" None Sequential

(* used in IO for implied alignments *)
let present_absent_io string =
    list_to_a [(string, 1, None); ("-", 2, None)]
              "-" None Sequential

(* The alphabet limited to the four bases *)
let dna =
    let all = adenine lor citosine lor guanine lor timine lor gap in
    list_to_a 
    [
        ("A", adenine, Some timine);
        ("C", citosine, Some guanine);
        ("G", guanine, Some citosine);
        ("T", timine, Some adenine);
        (gap_repr, gap, Some all);
        ("X", all, Some all)
    ] gap_repr (Some "X") Simple_Bit_Flags 

(* The alphabet of accepted IUPAC codes (up to N), and other codes used in the
* POY file format (_ up to |). *)
let nucleotides =
    let all = gap lor timine lor guanine lor adenine lor citosine in
    list_to_a
    [ 
        ("A", adenine, Some timine); 
        ("C", citosine, Some guanine); 
        ("G", guanine, Some citosine);
        ("T", timine, Some adenine); 
        ("U", uracile, Some adenine);
        ("M", adenine lor citosine, Some (timine lor guanine)); 
        ("R", adenine lor guanine, Some (timine lor citosine)); 
        ("W", adenine lor timine, Some (timine lor adenine)); 
        ("S", citosine lor guanine, Some (guanine lor citosine)); 
        ("Y", citosine lor timine, Some (guanine lor adenine)); 
        ("K", guanine lor timine, Some (citosine lor adenine)); 
        ("V", adenine lor citosine lor guanine, 
        Some (timine lor guanine lor citosine)); 
        ("H", adenine lor citosine lor timine, 
        Some (timine lor guanine lor adenine)); 
        ("D", adenine lor guanine lor timine,
        Some (timine lor citosine lor adenine)); 
        ("B", citosine lor guanine lor timine,
        Some (guanine lor citosine lor adenine)); 
        ("N", adenine lor citosine lor timine lor guanine,
        Some (timine lor guanine lor adenine lor citosine)); 
        ("X", adenine lor citosine lor timine lor guanine, 
        Some (timine lor guanine lor adenine lor citosine)); 
        (gap_repr, gap, Some gap); 
        ("1", 17, Some (all land (lnot 17)));
        ("2", 18, Some (all land (lnot 18)));
        ("3", 19, Some (all land (lnot 19)));
        ("4", 20, Some (all land (lnot 20)));
        ("5", 21, Some (all land (lnot 21)));
        ("6", 22, Some (all land (lnot 22)));
        ("7", 23, Some (all land (lnot 23)));
        ("8", 24, Some (all land (lnot 24)));
        ("9", 25, Some (all land (lnot 25)));
        ("0", 26, Some (all land (lnot 26)));
        ("!", 27, Some (all land (lnot 27)));
        ("^", 28, Some (all land (lnot 28)));
        ("$", 29, Some (all land (lnot 29)));
        ("#", 30, Some (all land (lnot 30)));
        ("*", 31, Some (all));
        ("?", 31, Some (all));
    ] gap_repr (Some "*") Extended_Bit_Flags 

        

(* The list of aminoacids *)
let aminoacids_char_list = 
    [
        ("A", alanine, None); 
        ("R", arginine, None); 
        ("N", asparagine, None); 
        ("D", aspartic, None); 
        ("C", cysteine, None); 
        ("Q", glutamine, None); 
        ("E", glutamic, None); 
        ("G", glycine, None); 
        ("H", histidine, None); 
        ("I", isoleucine, None); 
        ("L", leucine, None); 
        ("K", lysine, None); 
        ("M", methionine, None); 
        ("F", phenylalanine, None); 
        ("P", proline, None); 
        ("S", serine, None); 
        ("T", threonine, None); 
        ("W", tryptophan, None); 
        ("Y", tyrosine, None); 
        ("V", valine, None); 
        ("X", all_aminoacids, None); 
        (gap_repr, aa_gap, None);
    ]

let aminoacids =
    list_to_a aminoacids_char_list gap_repr (Some "X") Sequential 

let aminoacids_use_3d =
    list_to_a ~init3D:true aminoacids_char_list gap_repr (Some "X") Sequential 

(** [is_aminoacids] return true is alph is aminoacids, it might be with
    any level value, it doesn't matter*)
let is_aminoacids alph = 
    (alph.code_to_string = aminoacids.code_to_string)

let find_codelist comb alpha = match alpha.kind with
    | Continuous -> assert false
    | (Simple_Bit_Flags | Extended_Bit_Flags | Sequential |Combination_By_Level) ->
        try All_sets.IntegerMap.find comb alpha.comb_to_list 
        with | Not_found -> raise (Illegal_Code comb)


let find_comb codelist alpha = match alpha.kind with
    | Continuous -> assert false
    | (Simple_Bit_Flags | Extended_Bit_Flags | Sequential |Combination_By_Level) ->
        try All_sets.IntegerListMap.find codelist alpha.list_to_comb
        with | Not_found -> raise (Illegal_List codelist)

let match_base x alph = match alph.kind with
    | Continuous when x = "?" || x = "-" -> alph.gap
    | Continuous -> int_of_string x
    | (Simple_Bit_Flags | Extended_Bit_Flags | Sequential |Combination_By_Level) ->
        try All_sets.StringMap.find x alph.string_to_code 
        with | Not_found ->
            try All_sets.StringMap.find (String.uppercase x) alph.string_to_code 
            with | Not_found -> raise (Illegal_Character x)

let find_base = match_base

let match_code x alph = match alph.kind with
    | Continuous when x = alph.gap -> "?"
    | Continuous -> string_of_int x
    | (Simple_Bit_Flags | Extended_Bit_Flags | Sequential | Combination_By_Level) ->
        try All_sets.IntegerMap.find x alph.code_to_string
        with | Not_found -> raise (Illegal_Code x)

let find_code = match_code

let zero_indexed a =
    try ignore (match_code 0 a); true
    with Illegal_Code _ -> false

let set_ori_size alph size =
    {alph with ori_size = size}


(*[of_string x gap all] expand input string x with '~' sign if
    * orientation is true. then call [list_to_a] for a plain alphabet*)
let of_string ?(respect_case = false) ?(orientation = false) ?(init3D=false) x gap all =
    let debug = false in
    let osize = (List.length x)  in
    if debug then
        Printf.printf "Alphabet.of_string,orientation=%b,init3D=%b,original size=%d, call list_to_a\n%!" 
        orientation init3D osize;
    let rec builder alph counter = function
        | h :: t -> 
              if orientation then 
                  builder ((h, counter, None) :: ("~" ^ h, counter + 1, None) :: alph)
                      (counter + 2) t
              else 
                  builder ((h, counter, None):: alph) (counter + 1)  t

        | [] -> List.rev alph
    in
    let res = builder [] 1 x in
    (*call list_to_a with the expanded list and kind = Sequential*)
    let alpha = list_to_a ~respect_case:respect_case ~orientation:orientation ~init3D:init3D res gap all Sequential in
    if debug then print alpha;
    if debug then 
        Printf.printf "end of of_string,return alpha with ori_size=%d\n%!" osize;
    { alpha with ori_size = osize }

let size a = a.size
let get_orientation a = a.orientation
let get_ori_size a = a.ori_size

let use_3d a = a.threeD

let rnd a = 
    fun () ->
        let it = Random.int a.size in
        let finder a _ (cnt, res) =
            if cnt = it then (cnt + 1, Some a)
            else (cnt + 1, res)
        in
        match All_sets.IntegerMap.fold finder a.code_to_string (0, None) with
        | _, Some x -> x
        | _, None -> failwith "Alphabet.rnd"


let get_all a = a.all

let get_gap a = 
    a.gap

let get_missing _ = "?"

let get_level a = a.level

let set_level a newlevel = { a with level = newlevel }

let set_size a newsize = { a with size = newsize } 

(*[to_list a] return a list of pair (string,code),like (a,1), ordered by code*)
let to_list a =
    let res =
        All_sets.IntegerMap.fold (
            fun b a acc -> 
                (a, b) :: acc
            ) a.code_to_string 
        []
    in
    List.sort (fun (_, a) (_, b) -> a - b) res

module Lexer = struct
    (** A module to make a stream processor for a given alphabet type *)
    module OrderedChar = struct
        type t = char
        let compare (a : char) (b : char) = Pervasives.compare a b
    end

    module CM = Map.Make (OrderedChar)

    type p = Code of int | Unfinished of p CM.t

    let internal_lexer respect_case a =
        if debug then Printf.printf "alphabet.Lexer.internal_lexer,respect_case = %b\n%!" respect_case;
        let rec add_stream stream code acc =
            match acc with
            | Code x -> failwith "This alphabet is not prefix free1"
            | Unfinished set ->
                    try
                        let c = Stream.next stream in
                        if CM.mem c set then 
                            let nacc = CM.find c set in
                            match add_stream stream code nacc with
                            | Code x ->
                                    Printf.printf "we already have %c -> code=%d,%!" 
                                    c x;
                                    failwith "This alphabet is not prefix free2"
                            | res -> Unfinished (CM.add c res set)
                        else 
                            let nacc = Unfinished CM.empty in
                            let res = add_stream stream code nacc in
                            let nacc = CM.add c res set in
                            Unfinished nacc
                    with
                    | Stream.Failure ->
                            Code code
        in
        let lst = All_sets.StringMap.fold (fun a b acc ->
            let a = if not respect_case then String.uppercase a else a in
            (Stream.of_string a, b) :: acc) a.string_to_code []
        in
        List.fold_left (fun acc (a, b) ->
            add_stream a b acc) (Unfinished CM.empty) lst

    let rec single_processor issue_warnings respect_case stream acc = 
        function
        | Code x -> 
                x :: acc
        | Unfinished x ->
                let c = 
                    let c = Stream.next stream in
                    if respect_case then c
                    else Char.uppercase c 
                in
                try single_processor issue_warnings respect_case stream acc 
                (CM.find c x) with
                | Not_found as err ->
                    if issue_warnings then begin
                        Status.user_message Status.Error 
                            ("I@ could@ not@ find@ the@ character@ " ^ 
                             (StatusCommon.escape (String.make 1 c)) ^ 
                             "@ in@ position@ " ^ string_of_int (Stream.count stream));
                        Status.user_message Status.Error
                            ("I@ found@ an@ illegal@ character@ in@ " ^
                             "the@ " ^ "last@ file@ I@ was@ reading.");
                    end else 
                        ();
                    raise err

    let make_simplified_lexer style respect_case issue_warnings a =
        let lexer = internal_lexer respect_case a in
        fun stream ->
            let rec ignore_comment () =
                try
                    while true do
                        match Stream.peek stream with
                        | None -> failwith "Unterminated comment in Nexus file"
                        | Some ']' ->
                                Stream.junk stream;
                                raise Exit
                        | Some '[' ->
                                Stream.junk stream;
                                ignore_comment ()
                        | _ -> ()
                    done
                with
                | Exit -> ()
            in
            let rec full_processor acc = 
                match Stream.peek stream with
                | None -> failwith "Illegal non closed polymorphism"
                | Some v ->
                        match v with
                        | '\001' .. '\032' ->
                                Stream.junk stream;
                                full_processor acc
                        | '}' | ')' when style = `Nexus ->
                                Stream.junk stream;
                                acc
                        | '{' | '(' when style = `Nexus ->
                                failwith "Illegal nested polymorphism"
                        | ']' when style = `Hennig ->
                                Stream.junk stream;
                                acc
                        | '[' when style = `Hennig ->
                                failwith "Illegal nested polymorphism"
                        | _ ->
                                let res = 
                                    single_processor issue_warnings respect_case stream acc lexer
                                in
                                full_processor res
            in
            let rec processor_driver () =
                match Stream.peek stream with
                | None -> raise Exit
                | Some v ->
                        match v with
                        | '\001' .. '\032' ->
                                Stream.junk stream;
                                processor_driver ();
                        | '{' | '(' | '[' when style = `Nexus ->
                                Stream.junk stream;
                                if v <> '[' then
                                    full_processor []
                                else begin 
                                    ignore_comment ();
                                    processor_driver ()
                                end
                        | '[' when style = `Hennig ->
                                Stream.junk stream;
                                full_processor []
                        | _ -> single_processor issue_warnings respect_case stream [] lexer
            in
            processor_driver ()

    let make_lexer issue_warnings respect_case a = 
        let lexer = internal_lexer respect_case a in
        fun stream lst len ->
            let rec full_processor acc cnt =
                match Stream.peek stream with
                | Some v ->
                        begin match v with
                        | '\001' .. '\032' ->
                                Stream.junk stream;
                                full_processor acc cnt
                        | _ ->
                                let res = single_processor issue_warnings
                                respect_case stream acc lexer in
(*                                Printf.fprintf stdout " "; *)
                                full_processor res (cnt + 1)
                        end
                | None -> acc, cnt
            in
            let res, cnt = full_processor lst len in
            res, cnt
end


let kind alpha = alpha.kind

let simplify alph =
    if debug then Printf.printf "Alphabet.simplify\n%!";
     (* for full combination,We need to extract those numbers that only have one bit on 
     * for level combination, we need to get codes without combination*)
    match alph.kind with
    | Simple_Bit_Flags
    | Continuous
    | Sequential -> alph
    | Combination_By_Level ->
        if debug then Printf.printf "Combination_By_Level,call of_string\n%!";
        let ori_a_size = get_ori_size alph in
        let ori_alst:string list  = 
            let res = ref [] in
            for i = 0 to (ori_a_size - 1) do
                let pos = (i+1) in
                let str = 
                    try All_sets.IntegerMap.find pos alph.code_to_string 
                    with
                    | Not_found -> raise (Illegal_Code pos)
                in
                res := (!res)@[str]
            done;
            (!res)
        in
        assert(alph.orientation=false);
        assert((get_all alph)=None);
        (*of_string will give us a alphabet with kind=Sequential*)
        of_string ~orientation:false ~init3D:(use_3d alph) ori_alst gap_repr None 
    | Extended_Bit_Flags ->
            if debug then Printf.printf "Extended_Bit_Flags,call list_to_a\n%!";
            let gap = get_gap alph 
            and all = match get_all alph with
                | Some all -> all 
                | None -> (-1) (*failwith "Impossible"*)
            in
            (* Printf.printf "gap=%d,all=%d,level=%d,ori_size=%d\n%!"
             *               gap all alph.level alph.ori_size;          *)
            let has_one_bit_or_all v =
                if (check_level alph) then begin
                    (*if v<= gap then true
                    else false
                    *) true
                end
                else begin
                    if v = all then true
                    else
                        let rec has_only_one_bit_on v =
                            if v = 1 then true
                            else if 0 <> (1 land v) then false
                            else has_only_one_bit_on (v lsr 1)
                        in
                        has_only_one_bit_on v
                    ;
                end
            in
            let add_those_who_have_it v name acc =
                if has_one_bit_or_all v
                    then (name, v, None) :: acc
                    else acc
            in
            let list =
                All_sets.IntegerMap.fold add_those_who_have_it alph.code_to_string []
            in
            let allcode = match all with
                | (-1) -> None
                | _    -> Some (try find_code all alph with _ -> "*")
            in
            list_to_a ~orientation:alph.orientation list
                      (find_code gap alph) allcode  Simple_Bit_Flags


let rec to_sequential alph =
    let uselevel = check_level alph in
    match alph.kind with
    | Continuous | Sequential -> 
            alph
    | Extended_Bit_Flags ->  to_sequential (simplify alph) 
    | Combination_By_Level -> simplify alph 
    | Simple_Bit_Flags -> 
            (* We only really need to handle this case *)
            let all_code = 
                match get_all alph with
                | None -> max_int
                | Some x -> x 
            in
            let bit_to_code_position x = 
                let rec aux cnt x =
                    if x = 0 then (cnt - 1)
                    else aux (cnt + 1) (x lsr 1)
                in
                aux 0 x
            in
            let new_string_to_code =
                if uselevel then alph.string_to_code
                else
                All_sets.StringMap.fold (fun a b acc ->
                    if b = all_code then acc
                    else All_sets.StringMap.add a (bit_to_code_position b) acc) 
                alph.string_to_code All_sets.StringMap.empty
            in
            let new_code_to_string =
                if uselevel then alph.code_to_string
                else
                All_sets.IntegerMap.fold (fun a b acc ->
                    if a = all_code then acc
                    else All_sets.IntegerMap.add (bit_to_code_position a) b acc) 
                alph.code_to_string All_sets.IntegerMap.empty
            in
            let new_complement =
                All_sets.IntegerMap.fold (fun a b acc ->
                    if a = all_code then acc
                        else 
                            let to_add =
                                match b with
                                | None -> None
                                | Some x -> Some (bit_to_code_position x)
                            in
                            All_sets.IntegerMap.add (bit_to_code_position a)
                            to_add acc) alph.complement All_sets.IntegerMap.empty
            in
            let res = 
                { 
                    comb_to_list = alph.comb_to_list; 
                    list_to_comb = alph.list_to_comb;
                    level = alph.level;
                    ori_size = alph.ori_size;
                    string_to_code = new_string_to_code;
                    code_to_string = new_code_to_string;
                    complement = new_complement;
                    gap = All_sets.StringMap.find gap_repr new_string_to_code;
                    all = None;
                    size = 
                    All_sets.StringMap.fold (fun _ _ acc -> acc + 1)
                    new_string_to_code 0;
                    kind = Sequential;
                    orientation = alph.orientation;
                    threeD = alph.threeD }
            in
            res

(*[explote alph level ori_sz]. this is a recursive function, pass the alphabet with tag 'Sequential',
* this function will replace the code with bitwise code if it's doing full
* combination, if it's doing level combination, does nothing. 
* Then call itself again with tag 'Simple_Bit_Flags', expand the
* alphabet with possible combinations, return the alphabet with finish tag 
* 'Combination_By_Level' or 'Extended_Bit_Flags', end the function.*)
let rec explote alph level ori_sz=
    (* function "check_level" is not ready to use here, for the alphbet is not set up yet*)
    let uselevel =
        if (level>=1)&&(level<=ori_sz) then true
        else false
    in
    match alph.kind with
    | Continuous | Extended_Bit_Flags | Combination_By_Level ->  alph 
    | Simple_Bit_Flags -> 
            (* for full combination, We have each element as one bit, we have now to extend it into
            * all possible combinations, for combination by level, we use
            * continuous number for each element, extend the list just like full combination code*)
            if debug then Printf.printf
            "Alphabet.explote,Simple_Bit_Flags->Extended_Bit_Flags or Combination_By_Level\n%!";
            (* we do the List.rev here because we want the gap to be the first element in following combination calculation *)
            let list = 
                if uselevel then
                    List.rev(to_list alph)
                else (to_list alph)
            in
            let all_combinations =  
                (* sanity check *)
                assert (0 <> List.length list);
                let list = 
                    match alph.all with
                    | None -> list
                    | Some code -> 
                            List.filter (fun (_, b) -> b <> code) list
                in
                assert (0 <> List.length list);
                (* add all possible combinations to the alphabet *)
                let rec all_combinations lst =
                    match lst with
                    | h :: t -> 
                            let res = all_combinations t in
                            (*when doing level combination, get rid of any
                            * combination that has more elements than level*)
                            let newres = List.filter 
                            (fun x -> ( (List.length x)<level )) res in
                            (*Printf.printf "res = { %!";
                            List.iter (fun xlst -> 
                                Printf.printf "[%!";
                                List.iter (fun (x,y) -> Printf.printf "%s,%d;" x y) xlst;
                                Printf.printf "],%!";
                            ) newres;
                            Printf.printf "]\n%!"; *)
                            res @ (List.map (fun x -> h :: x) newres)
                    | [] -> [[]]
                in
                match all_combinations list with
                | [] :: ((_ :: _) as r) -> r
                | _ -> assert false
            in
            let new_comb_to_list = ref All_sets.IntegerMap.empty in
            let new_list_to_comb = ref All_sets.IntegerListMap.empty in
            let a_size = List.length list in
            let count = ref (a_size + 1 ) in
            let new_alphabet = 
                let merge_combination lst = 
                    (*merge a list of states into a string [a;b;c] -> "abc". then add it to the two maps*)
                    match lst with
                    | [(item, code)] -> 
                            if debug then 
                                Printf.printf "add [%s , %d] \n" item code;
                        if uselevel then begin
                            new_comb_to_list := 
                                All_sets.IntegerMap.add code [code] (!new_comb_to_list);
                            new_list_to_comb :=
                                All_sets.IntegerListMap.add [code] code (!new_list_to_comb);  
                        end else ();
                       (item, code, None)
                    | lst ->
                        assert (lst <> []);
                        if uselevel then
                        begin
                            let item, codelist =
                                List.fold_left 
                                (fun (acc_item, acc_codelist) (item, code) ->
                                     item ^ acc_item, code::acc_codelist 
                                ) ("", []) lst 
                            in
                            let item = "["^item^"]" in
                            let codelist = List.sort compare codelist in
                            if debug then 
                                Printf.printf "add  %s  ,  %d  \n" (item) (!count);
                            new_comb_to_list :=
                            All_sets.IntegerMap.add (!count) codelist (!new_comb_to_list);
                            new_list_to_comb :=
                            All_sets.IntegerListMap.add codelist (!count) (!new_list_to_comb);
                            incr count;
                            (item, (!count-1), None)
                            end
                        else begin
                            let item, code =
                            List.fold_left 
                            (fun (acc_item, acc_code) (item, code) ->
                                acc_item ^ item, (acc_code lor code))
                            ("[", 0) lst
                            in
                            if debug then 
                                Printf.printf "add  %s  ,  %d  \n" (item^"]") (code);
                            (item ^ "]", code, None)
                         end
                in
                List.map merge_combination all_combinations
            in
            let all_repr, _, _ = 
                match new_alphabet with
                | h :: t ->
                    List.fold_left 
                    (fun ((_, code, _) as acc) ((_, codet, _) as item) ->
                        if codet <= code then acc else item) h t
                | [] -> assert false
            in
            (*combination by level and full combination should have different tags.*)
            let res_kind = 
                if uselevel then Combination_By_Level 
                else Extended_Bit_Flags in
            let return_alpha = list_to_a ~orientation:alph.orientation new_alphabet
            gap_repr (Some all_repr) res_kind in
            let return_alpha =
            if uselevel then
                { return_alpha with 
                  comb_to_list = (!new_comb_to_list);  list_to_comb = (!new_list_to_comb);
                  level = level;
                  ori_size = ori_sz;
                  all = None 
                }
            else { return_alpha with ori_size = ori_sz }
            in
            return_alpha
    | Sequential ->
            if debug then Printf.printf "Alphabet.explote,Sequential->Simple_Bit_Flags\n%!";
            if uselevel then (*this part is not necessary for level combiantion?*)
                begin
                let new_alph_list =
                    let old_alph_list = to_list alph in
                    let count = ref 1 in
                    (*for level combination, code are continues*)
                    List.map (fun (str,_) ->
                        let res = !count in
                        count := !count + 1;
                        str,res, None) old_alph_list
                in
                (*call list_to_a with expanded (string,code,None) list, like (a,1,None),(b,2,None),(c,3,None),..... *)
                let res = list_to_a ~orientation:alph.orientation new_alph_list gap_repr None Simple_Bit_Flags in
                (*call itself again with new tag*)
                explote res level ori_sz 
                end
            else begin
                (*change code from (1,2,3,4...) to (1,2,4,8,....) for full combination*)
                let new_alphabet =
                    let list = to_list alph in
                    let bit_position = ref 1 in
                    (*for full combination, code are set for bitwise operation
                    * later, that is 2^n for each state in original alphabet*)
                    List.map (fun (str, _) ->
                        let res = !bit_position in
                        bit_position := !bit_position lsl 1;
                        str, res, None) list
                in
                (*call list_to_a with expanded (string,code,None) list, like (a,1,None),(c,2,None),(t,4,None),..... *)
                let res = list_to_a ~orientation:alph.orientation new_alphabet gap_repr None Simple_Bit_Flags in
                (*call itself again with new tag*)
                explote res level ori_sz
        end

let n_to_the_powers_of_m n m =
    let rec func t acc =
        if (t>=1) then func (t-1) n*acc
        else acc
    in
    func m 1

let create_alph_by_level alph level oldlevel =
    let debug = false in
    let ori_a_size = get_ori_size alph in
    if debug then 
        Printf.printf "Alphabet.create new alph: ori_size=%d,old_level/newlevel=%d/%d\n%!"
                        ori_a_size oldlevel level; 
    let orientation= get_orientation alph in
    let get_ori_alst  = 
        let res = ref [] in
        for i = 0 to (ori_a_size - 1) do
            let pos = (i+1) in
            let str = 
                try All_sets.IntegerMap.find pos alph.code_to_string 
                with
                | Not_found -> raise (Illegal_Code pos)
            in
            res := (!res)@[str]
        done;
        (!res)
    in
    let (ori_alst:string list) =   get_ori_alst  in
    let default_gap = gap_repr in
    let orialph = of_string ~orientation:orientation ori_alst default_gap None in
    let newalph = 
        if (level > ori_a_size) then 
            explote orialph ori_a_size ori_a_size 
        else if (level <1 ) then
            orialph
        else (*we need to do this even level=1*)
            explote orialph level ori_a_size
    in
    if debug then begin
        Printf.printf "end of Alphabet.create_alph_by_level, print newalph\n%!";
        print newalph;
    end;
    newalph
            
let distinct_size alph =
    All_sets.IntegerMap.fold (fun _ _ acc -> acc + 1) alph.code_to_string 0

let complement c alph =
    try 
    All_sets.IntegerMap.find c alph.complement
    with | Not_found ->
        Printf.printf  "cannot find complement of %d in alphabet" c;
        assert(false)

let complement2 c alph =
    let res = 
    try 
    All_sets.IntegerMap.find c alph.complement
    with | Not_found ->
        Printf.printf  "cannot find complement of %d in alphabet" c;
        assert(false)
    in
     match res with 
    | Some y -> y
    | None -> failwith "we are not expecting NONE complement"


(*NOTE: rev_comp_lst just give us the complement seq, not reverse, use List.rev to do that*)
let rev_comp_lst seqlst alph =
    List.map (fun x ->  complement2 x alph) seqlst

(*rev_comp_arr give us the reverse complement of seqarr*)
let rev_comp_arr seqarr alph = 
    let size = Array.length seqarr in
    Array.mapi (fun idx x -> complement2 seqarr.(size-idx-1) alph) seqarr 

(*[of_file] is only called by scripting.ml, when reading in custom alphabet and breakinversion file*)
let of_file fn orientation init3D level respect_case tie_breaker =
    let debug = false in
    if debug then Printf.printf
    "Alphabet.of_file,orientation=%b,level=%d,init3D=%b,respect_case=%b\n%!"
    orientation level init3D respect_case;
    let file = FileStream.Pervasives.open_in fn in
    let alph = FileStream.Pervasives.input_line file in
    let default_gap = gap_repr in
    let elts = ((Str.split (Str.regexp " +") alph) @ [default_gap]) in
    (*all_element is None for custom alphabet and breakinv*)
    let alph = of_string ~respect_case:respect_case ~orientation:orientation ~init3D:init3D
    elts default_gap None in 
    (*at this point, alph is a plain alphabet with only orientation but no
    * combination, size set to original size, level set to 0, two maps set to empty, *)
    let size = size alph in(*size here is still the original size without combination*)
    let alph, do_comb = 
        (*when orientation is true, it's breakinv, we don't do combination/level*)
        if orientation then alph, false 
        else explote alph level size, true
    in
    (*get cost matrix out of input cm file*)
    let tcm, matrix = 
        let all_elements = -1 (* breakinv and customalphabet don't have all_elements *) in
        if do_comb then
            Cost_matrix.Two_D.of_channel 
                ~tie_breaker:tie_breaker ~orientation:orientation ~level:level all_elements file 
        else
            Cost_matrix.Two_D.of_channel_nocomb
                ~orientation all_elements file
    in
    let tcm3 = match init3D with
        | true -> Cost_matrix.Three_D.of_two_dim tcm
        | false  ->  Cost_matrix.Three_D.default_nucleotides 
    in 
    file#close_in;
    if debug then Printf.printf "end of of_file, return alph and tcm\n%!";
    if debug then print alph;
    alph, (tcm,matrix), tcm3

(*    code_to_string : string All_sets.IntegerMap.t;    *)
(* vim: set et sw=4 tw=80: *)
