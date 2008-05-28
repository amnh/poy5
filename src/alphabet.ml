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

let () = SadmanOutput.register "Alphabet" "$Revision: 2871 $"


(* $Id: alphabet.ml 2871 2008-05-23 17:48:34Z andres $ *)

exception Illegal_Character of string
exception Illegal_Code of int

type kind = 
    | Sequential
    | Simple_Bit_Flags
    | Extended_Bit_Flags

type a = {
    string_to_code : int All_sets.StringMap.t;
    code_to_string : string All_sets.IntegerMap.t;
    complement : int option All_sets.IntegerMap.t;
    gap : int;
    all : int option;
    size : int;
    kind : kind;
    orientation : bool;
}

let print alpha = 
    All_sets.IntegerMap.iter (fun code char -> 
                       Printf.fprintf stdout "%6i %s\n" code char)
    alpha.code_to_string;
    print_newline ()

let to_formatter alph : Tags.output =
    let element_to_tags string code acc =
        let code = string_of_int code in
        (`Single (Tags.Alphabet.element, [(Tags.Alphabet.value, string);
        (Tags.Alphabet.code, code)], `Structured `Empty)) :: acc
    in
    let res = 
        `Set (All_sets.StringMap.fold element_to_tags alph.string_to_code [])
    in
    Tags.Characters.alphabet, [], `Structured res


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

let list_to_a ?(orientation=false)  lst gap all kind = 
    let add (s2c, c2s, cmp, cnt) (a, b, c) =
        All_sets.StringMap.add a b s2c,
        (if All_sets.IntegerMap.mem b c2s then c2s
        else All_sets.IntegerMap.add b a c2s),
        (All_sets.IntegerMap.add b c cmp),
        cnt + 1
    in
    let empty = All_sets.StringMap.empty, All_sets.IntegerMap.empty,
    All_sets.IntegerMap.empty, 0 in
    let s2c, c2s, cmp, cnt = List.fold_left add empty lst in
    let gap_code = 
        try All_sets.StringMap.find gap s2c with
        | Not_found as err ->
                List.iter (fun (x, _, _) -> 
                    Status.user_message Status.Error x) lst;
                Status.user_message Status.Error
                ("could not find the gap " ^ gap);
                raise err
    and all_code = 
        match all with
        | Some all -> Some (All_sets.StringMap.find all s2c)
        | None -> None
    in
    { string_to_code = s2c; code_to_string = c2s; gap = gap_code; all = all_code;
      size = cnt; kind = kind; complement = cmp; orientation = orientation }

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
let aminoacids =
    list_to_a
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
    ] gap_repr (Some "X") Sequential 

let match_base x alph =
    try
        let x = String.uppercase x in
        All_sets.StringMap.find x alph.string_to_code 
    with
    | Not_found -> raise (Illegal_Character x)



let find_base = match_base

let match_code x alph =
    try All_sets.IntegerMap.find x alph.code_to_string with
    | Not_found -> raise (Illegal_Code x)

let find_code = match_code

let of_string ?(orientation = false) x gap all =
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
    let alpha = list_to_a ~orientation:orientation res gap all Sequential in
    alpha

let size a = a.size
let get_orientation a = a.orientation

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

let get_gap a = a.gap

let to_list a =
    let res =
        All_sets.StringMap.fold (fun a b acc -> (a, b) :: acc) a.string_to_code 
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
        let rec add_stream stream code acc =
            match acc with
            | Code x -> failwith "This alphabet is not prefix free"
            | Unfinished set ->
                    try
                        let c = Stream.next stream in
                        if CM.mem c set then 
                            let nacc = CM.find c set in
                            match add_stream stream code nacc with
                            | Code _ -> failwith "This alphabet is not prefix free"
                            | res -> Unfinished (CM.add c res set)
                        else 
                            let nacc = Unfinished CM.empty in
                            let res = add_stream stream code nacc in
                            let nacc = CM.add c res set in
                            Unfinished nacc
                    with
                    | Stream.Failure -> Code code
        in
        let lst = All_sets.StringMap.fold (fun a b acc ->
            let a = if not respect_case then String.uppercase a else a in
            (Stream.of_string a, b) :: acc) a.string_to_code []
        in
        List.fold_left (fun acc (a, b) ->
            add_stream a b acc) (Unfinished CM.empty) lst

    let rec single_processor issue_warnings respect_case stream acc = function
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
                            "@ in@ position@ " ^
                            string_of_int (Stream.count stream));
                            Status.user_message Status.Error
                            ("I@ found@ an@ illegal@ character@ in@ " ^
                            "the@ " ^ "last@ file@ I@ was@ reading.");
                        end else ();
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

    let make_lexer issue_warnings a = 
        let lexer = internal_lexer false a in
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
                                false stream acc lexer in
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
    match alph.kind with
    | Simple_Bit_Flags
    | Sequential -> alph
    | Extended_Bit_Flags ->
            (* We need to extract those numbers that only have one bit on *)
            let gap = get_gap alph
            and all = 
                let all = get_all alph in
                match all with
                | Some all -> all 
                | None -> failwith "Impossible"
            in
            let has_one_bit_or_all v =
                if v = all then true
                else
                    let rec has_only_one_bit_on v =
                        if v = 1 then true
                        else if 0 <> (1 land v) then false
                        else has_only_one_bit_on (v lsr 1)
                    in
                    has_only_one_bit_on v
            in
            let add_those_who_have_it v name acc =
                if has_one_bit_or_all v then (name, v, None) :: acc
                else acc
            in
            let list = 
                All_sets.IntegerMap.fold add_those_who_have_it 
                alph.code_to_string []
            in
            list_to_a ~orientation:alph.orientation list (find_code gap alph) 
            (Some (try find_code all alph with _ -> "*"))  Simple_Bit_Flags 

let rec to_sequential alph =
    match alph.kind with
    | Sequential -> 
            alph
    | Extended_Bit_Flags -> 
            to_sequential (simplify alph) 
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
                All_sets.StringMap.fold (fun a b acc ->
                    if b = all_code then acc
                    else All_sets.StringMap.add a (bit_to_code_position b) acc) 
                alph.string_to_code All_sets.StringMap.empty
            in
            let res = 
                { string_to_code = new_string_to_code;
                code_to_string =
                    All_sets.IntegerMap.fold (fun a b acc ->
                        if a = all_code then acc
                        else All_sets.IntegerMap.add (bit_to_code_position a) b acc) 
                    alph.code_to_string All_sets.IntegerMap.empty;
                complement =
                    All_sets.IntegerMap.fold (fun a b acc ->
                        if a = all_code then acc
                        else 
                            let to_add =
                                match b with
                                | None -> None
                                | Some x -> Some (bit_to_code_position x)
                            in
                            All_sets.IntegerMap.add (bit_to_code_position a)
                            to_add acc) alph.complement
                            All_sets.IntegerMap.empty;
                gap = All_sets.StringMap.find gap_repr new_string_to_code;
                all = None;
                size = 
                    All_sets.StringMap.fold (fun _ _ acc -> acc + 1)
                    new_string_to_code 0;
                    kind = Sequential;
                    orientation = alph.orientation }
            in
            res

let rec explote alph =
    match alph.kind with
    | Extended_Bit_Flags -> (* We already have what we want *) alph 
    | Simple_Bit_Flags -> 
            (* We have each element as one bit, we have now to extend it into
            * all possible combinations *)
            let all_combinations =  
                let list = to_list alph in
                assert (0 <> List.length list);
                let list = 
                    match alph.all with
                    | None -> list
                    | Some code -> 
                            List.filter (fun (_, b) -> b <> code) list
                in
                assert (0 <> List.length list);
                let rec all_combinations lst =
                    match lst with
                    | h :: t -> 
                            let res = all_combinations t in
                            res @ (List.map (fun x -> h :: x) res)
                    | [] -> [[]]
                in
                match all_combinations list with
                | [] :: ((_ :: _) as r) -> r
                | _ -> assert false
            in
            let new_alphabet = 
                let merge_combination lst =
                    match lst with
                    | [(item, code)] -> (item, code, None)
                    | lst ->
                            let item, code =
                                List.fold_left 
                                (fun (acc_item, acc_code) (item, code) ->
                                    acc_item ^ item, (acc_code lor code))
                                ("[", 0) lst
                            in
                            (item ^ "]", code, None)
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
            list_to_a ~orientation:alph.orientation new_alphabet gap_repr (Some all_repr) Extended_Bit_Flags 
    | Sequential ->
            let new_alphabet =
                let list = to_list alph in
                let bit_position = ref 1 in
                List.map (fun (str, _) ->
                    let res = !bit_position in
                    bit_position := !bit_position lsl 1;
                    str, res, None) list
            in
            let res = list_to_a ~orientation:alph.orientation new_alphabet gap_repr None Simple_Bit_Flags in
            explote res

let distinct_size alph =
    All_sets.IntegerMap.fold (fun _ _ acc -> acc + 1) alph.code_to_string 0

let complement c alph = 
    All_sets.IntegerMap.find c alph.complement
(*    code_to_string : string All_sets.IntegerMap.t;    *)
(* vim: set et sw=4 tw=80: *)
