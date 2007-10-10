val fprintf : out_channel -> ('a, out_channel, unit) format -> 'a
type breakinv_t = BreakinvAli.breakinv_t
type meds_t = {
  med_ls : breakinv_t list;
  num_med : int;
  total_cost : int;
  total_recost : int;
  breakinv_pam : Data.dyna_pam_t;
  gen_cost_mat : Cost_matrix.Two_D.m;
  pure_gen_cost_mat : int array array;
  alpha : Alphabet.a;
}
val init_med :
  Sequence.s ->
  Cost_matrix.Two_D.m -> Alphabet.a -> Data.dyna_pam_t -> meds_t
val keep : Data.dyna_pam_t -> 'a list -> 'a list
val find_meds2 : meds_t -> meds_t -> meds_t
val find_meds3 : meds_t -> meds_t -> meds_t -> meds_t
val cmp_min_pair_cost : meds_t -> meds_t -> int * int
val cmp_max_pair_cost : meds_t -> meds_t -> int * int
val compare : meds_t -> meds_t -> int
val get_active_ref_code : meds_t -> int * int * int
val readjust_3d : meds_t -> meds_t -> meds_t -> 'a -> 'b -> meds_t -> int * meds_t * bool

