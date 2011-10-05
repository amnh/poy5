
type annotate_tool_t = 
    [ `Mauve of (float*int*float*int) | `Default of (int*int*int) ]

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
    max_kept_wag : int option; 
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
}


