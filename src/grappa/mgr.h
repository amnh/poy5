
#include "mgrstructs.h"
#include "mgr_scenario.h"

G_struct Genomes;
G_struct Genomes_copy;
G_struct Preancestors;

struct mgr_genome_struct *mgr_genome_list;
struct mgr_genome_struct *mgr_genome_list_copy;
struct mgr_genome_struct median_genome;

mgr_distmem_t distmem_mgrmed;

cbounds_t  * cb_max_chromo_size;// mgr_genome_ops.c, "find_max_chromo_size"
//cbounds_t  * cb_carry_on_reag;//mgr_genome_ops.c, "carry_on_reag" 
//cbounds_t  * cb_print_one_reag;//mgr_genome_ops.c, "print_one_reag"
//cbounds_t  * cb_build_list_reag; //mgr_my_testrev.c, "build_list_reag" 
//cbounds_t  * cb_compute_Rscore; //mgr_list_ops, "compute_Rscore"

list_listreag list_carry_best_reag;
list_listreag list_carry_best_reag2;

list_listreag list_do_k_reag;
//list_listreag list_find_k_reag2;

int carry_best_reag(G_struct *Genomes , int *nbreag, int nb_spec,
					int *spec_left,
					int reversals,      // TRUE if we want to include reversals
					int transloc, 	// TRUE if we want to include translocations
					int fusion,
					int fission, // TRUE if we want to include fusions and fissions
					mgr_distmem_t *distmem, treemem_t *treemem,
					int *reduction_type,
					int *last_genomenb, int heuristic,
					int verbose, int pair1, int pair2, // usually both equal to -1 but if we want to merge two nodes,
					// we use pair1 pair2 to indicate which node need to be merged
					struct mgr_genome_struct *old_genome);


void solve_with_good_reag(G_struct *Genomes, int *nbreag,
						  int nb_spec, int *spec_left, 
						  mgr_distmem_t *distmem, treemem_t *treemem, 
						  int heuristic, int pair1, int pair2, int verbose) ;

void mgr_ini_mem (int num_genes);

int mgr_invdist (int * g1, int * g2, int num_genes, int * deli1, int * deli2, int num_deli1, int num_deli2);

void mgr_med(int * g1, int * g2,int * g3, int * deli1, int * deli2, int * deli3, int num_deli1, int num_deli2, int num_deli3, int SIZE_ALPHA, int CIRCULAR, struct genome_struct * g_med);

int solve_triplet (G_struct *Genomes, int *nbreag, int *perfect_triplet, int reduction_type, int nb_spec, int *spec_left, int depth, mgr_distmem_t *distmem, treemem_t *treemem,
G_struct *Preancestors,int heuristic, int condensing, int verbose) ;
