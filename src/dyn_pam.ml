(* POY 4.0 Beta. A phylogenetic analysis program using Dynamic Homologies.    *\
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
\* USA                                                                        *)

type annotate_tool_t = 
    [ `Mauve of (float*float*float*float) | `Default of (int*int*int) ]

type median_solver_t =
    [ `MGR | `Vinh | `Albert | `Siepel | `BBTSP | `COALESTSP | `ChainedLK | `SimpleLK ]

type align_meth_t =
    [ `NewKK | `Default ]

type re_meth_t = 
    [ `Locus_Breakpoint of int | `Locus_Inversion of int ]

type dyna_pam_t = {

    align_meth : align_meth_t option; 

    median_solver: median_solver_t option; 

    annotate_tool: annotate_tool_t option;

    (** Cost parameters of rearrangement function which is either
    * breakpoint distance or inversion distance *)
    re_meth : re_meth_t option;

    (** Circular = 0 means linear chromosome, otherwise circular chromosome *)
    circular : int option;

    (** [(a, b)] is indel cost of a locus in a chromosome
    * where [a] is opening cost, [b] is extension cost *)
    locus_indel_cost : (int * int) option;

    (** [(a, b)] is indel cost of a chromosome in a genome 
    * where [a] is opening cost, [b] is extension cost *)
    chrom_indel_cost : (int * int) option;

    (** The maximum cost between two chromosomes
    * at which they are considered as homologous chromosomes *)
    chrom_hom : int option;

    (** The cost of a breakpoint happing between two chromosome *)
    chrom_breakpoint : int option;
    
    (** The maximum number of medians at one node kept during the search*)
    keep_median : int option;

    (** number iterations are applied in refining alignments with rearrangements *)
    swap_med : int option; 

    (** approx = true, the median sequence of X and Y is approximated by either X or Y,
    * otherwise, calculate as a set of median between X and Y *)
    approx : bool option;

    (** symmetric = true, calculate the both distances between X to Y
     * and vice versa. Otherwise, only form X to Y *) 
    symmetric : bool option;

    (** maximum length of sequences aligned by 3D-alignment *)
    max_3d_len : int option;

    (** maximum of Wagner-based alignments are kept during 
    * the pairwise alignment with rearrangements *)
    max_kept_wag: int option; 

    (** what variables are enabled in this structure; since things are
    * demarcated properly **)
    mode : [`Chromosome | `Genome | `Breakinv ] option
}

(** [dyna_pam_default] assigns default values for parameters 
* used to create the median between two chromosomes or genomes *)
let dyna_pam_default = {
    align_meth = Some `Default;
    median_solver = Some `Albert;
    annotate_tool = Some (`Default (100,100,9));
    re_meth = Some (`Locus_Breakpoint 10);
    circular = Some 0;
    locus_indel_cost = Some (10, 100);
    chrom_indel_cost = Some (10, 100);
    chrom_hom = Some 200;
    chrom_breakpoint = Some 100;
    keep_median = Some 1;
    swap_med = Some 1;
    approx = Some false;
    symmetric = Some true;
    max_3d_len = Some 100;
    max_kept_wag = Some 1;
    mode = None;
}

let rec meth_to_string = function
    | Some `MGR       -> "mgr"
    | Some `Vinh      -> "vinh"
    | Some `Albert    -> "albert"
    | Some `Siepel    -> "siepel"
    | Some `BBTSP     -> "bbtsp"
    | Some `COALESTSP -> "coalestsp"
    | Some `ChainedLK -> "chainedlk"
    | Some `SimpleLK  -> "simplelk"
    | None            -> assert false

let annotate_to_nexus fo annot : unit = match annot with
    | Some `Default (x,y,z) ->
        fo ("@[Annotations@, ");
        fo ("@[Model = Default;@]@, ");
        fo ("@[Min = "^(string_of_int x)^";@]@, ");
        fo ("@[Max = "^(string_of_int y)^";@]@, ");
        fo ("@[Rearrangement = "^(string_of_int y)^";@]@, ");
        fo ("@,;@,@]")
    | Some `Mauve (w,x,y,z) ->
        fo ("@[Annotations@, ");
        fo ("@[Model = Mauve;@]@, ");
        fo ("@[Quality = "^(string_of_float w)^";@]@, ");
        fo ("@[Coverge = "^(string_of_float x)^";@]@, ");
        fo ("@[Min = "^(string_of_float z)^";@]@, ");
        fo ("@[Max = "^(string_of_float y)^";@]@, ");
        fo ("@,;@,@]")
    | None ->
        assert false

let print_characters fo chars = 
    match chars with
    | []    -> ()
    | [x]   -> fo "@,@]:";
               fo x
    | x::xs -> fo "@,@]:";
               fo x;
               List.iter (fun x -> Printf.printf ", %s" x) xs

let to_nexus fo pam chars : unit =
    let rec get_some = function | Some x -> x | None -> assert false in
    (* we assume nothing is NONE *)
    match pam.mode with
    | Some `Chromosome ->
        fo "@[Chromosome@, ";
        fo ("@[solver = "^meth_to_string pam.median_solver^";@]@, ");
        fo ("@[median = "^string_of_int (get_some pam.keep_median)^";@]@, ");
        fo (if get_some pam.symmetric then "@[symmetric;@]@, " else "");
        fo (match pam.re_meth with
            | Some `Locus_Breakpoint x -> "@[Breakpoint = "^(string_of_int x)^";@]@, "
            | Some `Locus_Inversion x  -> "@[Inversion = "^(string_of_int x)^";@]@, "
            | None -> assert false);
        fo (if get_some pam.approx then "@[approximate;@]@, " else "");
        annotate_to_nexus fo pam.annotate_tool;
        fo (match pam.locus_indel_cost with
            | Some (x,y) ->
                let y = (float_of_int y) /. 100.0 in
                "@[indel = "^(string_of_int x)^","^(string_of_float y)^";@]@, "
            | None -> assert false);
        print_characters fo chars;
        fo "@,;@]@."

    | Some `Breakinv ->
        fo "@[BreakInversion@, ";
        fo ("@[solver = "^meth_to_string pam.median_solver^";@]@, ");
        print_characters fo chars;
        fo "@,;@]@."

    | Some `Genome ->
        fo "@[Genome@, ";
        fo ("@[median = "^string_of_int (get_some pam.keep_median)^";@]@, ");
        fo (if 1 = (get_some pam.circular) then "@[circular;@]@, " else "");
        fo ("@[breakpoint = "^string_of_int (get_some pam.chrom_breakpoint)^";@]@, ");
        fo ("@[distance = "^string_of_float ((float_of_int (get_some pam.chrom_hom))/. 100.0) ^";@]@, ");
        fo (match pam.locus_indel_cost with
            | Some (x,y) ->
                let y = (float_of_int y) /. 100.0 in
                "@[indel = "^(string_of_int x)^","^(string_of_float y)^";@]@, "
            | None -> assert false);
        print_characters fo chars;
        fo "@,;@]@."

    | None -> ()
