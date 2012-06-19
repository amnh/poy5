

type elt = SankCS.elt (*cside*)
type t = SankCS.t (*cside*)

type t_w_seqtbl = {
    states_and_distbl : t;
    sequence_table: Sequence.s array; (*we carry original sequence before
    transform(fixed_state) with us*)
    alph : Alphabet.a;
}

val cardinal : t_w_seqtbl -> int

val f_codes : t_w_seqtbl -> All_sets.Integers.t -> t_w_seqtbl

val f_codes_comp : t_w_seqtbl -> All_sets.Integers.t -> t_w_seqtbl

val get_states : t_w_seqtbl -> int -> int array

val get_earray : t_w_seqtbl -> int array

val get_min_states : t_w_seqtbl -> int * int

val off_array : Data.fixed_state_spec -> int -> t_w_seqtbl 

val distance : t_w_seqtbl -> t_w_seqtbl -> float

val median : int -> t_w_seqtbl -> t_w_seqtbl -> t_w_seqtbl * float

val median_3 : t_w_seqtbl ->  t_w_seqtbl ->  t_w_seqtbl ->  t_w_seqtbl -> t_w_seqtbl 

val dist_2 : t_w_seqtbl -> t_w_seqtbl -> t_w_seqtbl -> float 

val to_single : t_w_seqtbl -> t_w_seqtbl -> t_w_seqtbl * int

val to_string : t_w_seqtbl -> string

val compare_data : t_w_seqtbl -> t_w_seqtbl -> int 

val to_formatter :
    Methods.diagnosis_report_type -> Xml.attributes -> t_w_seqtbl -> Data.d 
    -> Xml.xml Sexpr.t list 


