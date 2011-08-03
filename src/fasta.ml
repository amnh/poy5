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

type taxon = string

let debug = false

let (-->) b a = a b

let true_trim str = 
    Str.global_replace (Str.regexp "\\(^ +\\)\\|\\( +$\\)") "" str

module type S = sig
    (** Interface for molecular kind of data parsers *)

    (** [of_channel x y] takes an input channel [y] with information in ascii
    * format of a molecular data of type [x] and outputs the list of sequences 
    * and
    * taxon names contained in the file. If the function finds an unexpected
    * character or an illegal condition for a specific file format, an
    * Illegal_molecular_format or Unsupported_file_format exception is raised. 
    * *)
    val of_channel : ?respect_case:bool -> FileContents.t -> in_channel -> 
        (Sequence.s list list list * taxon) list

    (** [to_channel x y z] writes in the output channel [x] the sequences and 
    * taxon list [y] using the alphabet [z] for the sequences in [y]. If the 
    * function finds
    * an illegal element in any of the sequences for the alphabet [z], raises an
    * Alphabet.Illegal_Code exception. There is no guarantee on the state of the
    * output file if the exception is raised. *)
    val to_channel : 
        out_channel -> (Sequence.s * taxon) list -> Alphabet.a -> unit

    val of_file : ?respect_case:bool -> FileContents.t -> FileStream.f -> 
        (Sequence.s list list list * taxon) list 
end


let alphabet_of_t t =
    match t with
    | FileContents.Nucleic_Acids -> Alphabet.nucleotides
    | FileContents.Proteins -> Alphabet.aminoacids
    | FileContents.Prealigned_Alphabet a
    | FileContents.AlphSeq a -> a
    | _ -> failwith "alphabet_of_t"

let empty_line = Str.regexp "^[ \t]*$"

type read = Eof | Eot | Read of string

let get_line () =
    fun ch ->
    try 
        let line = FileStream.Pervasives.input_line ch in
        if Str.string_match empty_line line 0 then Eot else Read line
    with
    | _ -> Eof

let is_taxon line = line.[0] = '>'
exception Finished

let split_subsequences pattern sequence =
    let divider = Str.regexp pattern in        
    let builder item (prev, acc)  =
        let items = Str.split_delim divider item in 
        let len = List.length items in 
        match items with
        | [h] -> (h :: prev), acc
        | [] -> prev, acc
        | hd :: tl ->
              let acc, _ = 
                  List.fold_right 
                      (fun item (acc, idx) ->
                           if idx = len - 1 then  
                               (item::prev)::acc, idx - 1  
                           else if idx = 0 then acc, idx - 1 
                           else ([item]::acc), idx - 1 
                      ) items (acc, len - 1) 
              in  
              [hd], acc 
    in  
    let lst =  
        match List.fold_right builder sequence ([], []) with 
        | [], tl -> tl 
        | hd, tl -> hd :: tl  
    in  
    lst

let split_frag = split_subsequences "#"
let split_loci = split_subsequences "|"
let split_chromosomes = split_subsequences "@"

type fl = {
    filename : string;
    taxon : string;
    sequence : string;
    character : string;
    line : int;
}

exception Illegal_molecular_format of fl


let doprepend alph seq v = 
    try Sequence.prepend seq (Alphabet.match_base v alph) with
    | Alphabet.Illegal_Character c -> 
            if v <> " " then
                let fl = 
                    { filename = ""; taxon = ""; sequence = ""; character = c; 
                    line = 0 }
                in
                raise (Illegal_molecular_format fl)
            else ()

let process_sequence remove_gaps lexer alph lst =
    let lst = List.rev lst in
    let gap = Alphabet.get_gap alph in
    let seq, length =
        List.fold_right (fun x (lst, cnt) ->
                let stream = Stream.of_string x in
                lexer stream lst cnt) lst ([], 0)
    in
    (* Watch it, the sequence list is in reversed order *)
    let s = Sequence.create (length + 1) in
    List.iter (fun x -> 
        if x <> gap || not remove_gaps then
            Sequence.prepend s x
        else ()) seq;
    Sequence.prepend s gap;
    s

let process_taxon_name comment_str name =
    (* build a regexp to remove from 'comment_str' to end of the line *)
    let regexp = Str.regexp ((Str.quote comment_str)^".*$") in
    true_trim (Str.replace_first regexp "" name)

let process_file_imp ?(respect_case=false) remove_gaps ch alph  =
    if debug then Printf.printf "process_file_imp,respect_case=%b\n%!" respect_case;
    (* Apparenty ocaml is not making the next tail-recursive function a
     * loop? for some reason I get stack overflow, so I will write it using a
     * loop and see how it goes. *)
    let res = ref [] 
    and sequence = ref []
    and lexer = Alphabet.Lexer.make_lexer true respect_case alph 
    and taxon = ref "" in
    let get_line = get_line () in
    let last_line_empty = ref true in
    try
        while true do
            match get_line ch with
            | Read line ->
                if is_taxon line || !last_line_empty then begin
                    try let tmps =  
                            !sequence
                                --> List.rev
                                --> split_chromosomes
                                --> List.map split_loci
                                --> List.map (List.map split_frag)
                                --> List.map (List.map (List.map (process_sequence remove_gaps lexer alph)))
                        in 
                        let result = tmps, (process_taxon_name "$" !taxon) in
                        res := result :: !res;
                        taxon := line;
                        if (!taxon).[0] = '>' then
                            taxon := String.sub line 1 ((String.length line) - 1);
                        sequence := [];
                        last_line_empty := false;
                    with | Illegal_molecular_format fl ->
                        let fl = { fl with taxon = line } in
                        raise (Illegal_molecular_format fl)
                end else (sequence := line :: !sequence;)
            | Eot -> last_line_empty := true
            | Eof -> 
                    let tmps = 
                        !sequence 
                            --> List.rev
                            --> split_chromosomes
                            --> List.map split_loci
                            --> List.map (List.map split_frag)
                            --> List.map (List.map (List.map (process_sequence remove_gaps lexer alph)))
                    in
                    let result = tmps, (process_taxon_name "$" !taxon) in
                    res := result :: !res;
                    raise Finished
        done;
        []
    with
    | Finished -> 
            !res

exception Unsupported_file_format of string

let of_channel_obj ?(respect_case=false) t ch  =
    let res = 
        match t with
        | FileContents.Nucleic_Acids 
        | FileContents.AlphSeq _ 
        | FileContents.Proteins -> process_file_imp ~respect_case:respect_case true ch (alphabet_of_t t) 
        | FileContents.Prealigned_Alphabet _ ->
                process_file_imp  ~respect_case:respect_case false ch (alphabet_of_t t)
        | _ -> 
                let msg = 
                    "Unexpected error 01. Contact the AMNH programming team." 
                in
                raise (Unsupported_file_format msg)
    in
    List.fold_left (fun acc ((_, b) as item) -> 
        if b <> "" then item :: acc
        else acc) [] res

let of_channel ?(respect_case=false) t ch  =
    of_channel_obj ~respect_case:respect_case t (FileStream.stream_reader ch) 

let of_string ?(respect_case=false) t str =
    of_channel_obj ~respect_case:respect_case t (new FileStream.string_reader str) 

let output ch alp (seq, name) = 
    let str = Sequence.to_string seq alp in
    output_string ch (">" ^ name ^ "\n");
    output_string ch (str ^ "\n")

let to_channel ch l alp =
    List.iter (output ch alp) l

let of_file ?(respect_case=false) t f =
    try 
        let ch = FileStream.Pervasives.open_in f in
        let res = of_channel_obj ~respect_case:respect_case t ch in
        FileStream.Pervasives.close_in ch;
        res
    with
    | Illegal_molecular_format fl ->
            let f = FileStream.filename f in
            let fl = { fl with filename = f } in
            raise (Illegal_molecular_format fl)
