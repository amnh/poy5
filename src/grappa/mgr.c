/* head file from grappa */
#include "structs.h"
/* head files from mgr */

#include "mgr.h"
//#include "mgrstructs.h"
#include "mgr_genome_print.h"
#include "mgr_genome_ops.h"
#include "mgr_list_ops.h"
#include "mgr_my_testrev.h"
#include "mgr_genome_ops.h"
//#include "mgr_write_data.h"



//[1.3.1] from rearrange.c
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
					struct mgr_genome_struct *old_genome) 
{

    struct mgr_genome_struct *genome_list = Genomes->genome_list;
    list_reag *the_list = NULL, *temp_list = NULL;
    int score_list, total_dist, found_something = FALSE;
	int did_merge;
/* let's consider heuristic issue later   		
	if (heuristic != 0) {
		if (heuristic == 1) {
			return carry_best_reag_h1(Genomes, nbreag, nb_spec, spec_left, reversals, transloc, fusion, fission,
										distmem, treemem, reduction_type, last_genomenb,
										verbose, pair1, pair2, old_genome);
		}
		else {
			if (heuristic == 2) {
				return carry_best_reag_h2(Genomes, nbreag, nb_spec, spec_left, reversals, transloc, fusion, fission,
										  distmem, treemem, reduction_type, last_genomenb,
										  verbose, pair1, pair2, old_genome);
			}
		}
	}
	else
    */
    {
    
		// Regular case without heuristic
	    // first compute current total distance
		total_dist = find_total_dist(Genomes, nb_spec, distmem);
		if (total_dist == 0)
			return FALSE;
	
        // build the list of good rearrangements
		score_list = build_list_reag(&the_list, Genomes, nb_spec, *spec_left,
									 reversals, transloc, fusion, fission,
									 distmem, FALSE, // FALSE because we want the list, not just it's size
									 FALSE, *reduction_type, -1, // we look for all rearrangements
									 *last_genomenb, pair1, pair2, verbose);

		if (verbose) {
			fprintf(stdout,"\ntotal distance: %d and score_list %d (-%d), spec_left %d\n",
					total_dist, score_list, *reduction_type, *spec_left);
			fflush(stdout);
		}


		if (verbose) {
			//print_list_reag(the_list, nb_spec);

			print_dist_mat(Genomes, nb_spec);
			fflush(stdout);
		}

		if (the_list != NULL) { // we did find a good rearrangement
	
			int max_reduction=-1, new_reduction;
			int max_index = -1, the_index=1;
			
			temp_list=the_list;
			while (temp_list!=NULL) {
				list_reag *new_list = NULL;
					

				// carry-on the rearrangement with chromo flip if necessary
				copy_genes(genome_list[temp_list->spec].genes, old_genome->genes, genome_list[temp_list->spec].num_g);
				carry_on_reag(temp_list, Genomes, nb_spec, 
							  spec_left, &did_merge, distmem);
				

				new_reduction = build_list_reag(&new_list, Genomes, nb_spec, *spec_left,
												TRUE, (genome_list[0].num_chr>0)?TRUE:FALSE,
												(genome_list[0].num_chr>0)?TRUE:FALSE,(genome_list[0].num_chr>0)?TRUE:FALSE,
													/* because we want to know the total number of
													* potential rearragements */
												distmem, TRUE, // TRUE because we only want size
												FALSE, *reduction_type, -1, temp_list->spec, pair1, pair2, verbose);
				
				if (new_reduction>max_reduction) {
					//this one is better
					max_reduction = new_reduction;
					max_index = the_index;
				}
				
				// in any case undo last reversals to check for others...
				copy_genes(old_genome->genes, Genomes->genome_list[temp_list->spec].genes,  Genomes->genome_list[temp_list->spec].num_g);
				undo_reag(temp_list, Genomes, nb_spec, spec_left, did_merge, distmem);

				the_index++;
				temp_list=temp_list->next;
			}
				
			if (max_index!=-1) { // we found the best rearrangement
				int i;
				
				temp_list = the_list;
            // we go back to the best reversal we found
				for (i=1;i<max_index;i++) {
					temp_list=(*temp_list).next;
				}
/*	
 *	print_one_reag is in genome_ops.c and it's long....
            	if (verbose) {
					print_one_reag(temp_list, Genomes, nb_spec, old_genome);
				}
*/			
				carry_on_reag(temp_list, Genomes, nb_spec,							
							  spec_left, &did_merge, distmem);
/*				
				if (verbose) {
					special_print_genomes2(stdout, Genomes, nb_spec);
				}
*/			
				nbreag[temp_list->spec]++;
				*last_genomenb = temp_list->spec;
				found_something = TRUE ;
				
				//potentially update max_chromo_size
				find_max_chromo_size(Genomes, nb_spec);
				
				if (temp_list->merge_with != -1) {// there was a merge and we need to update tree
					update_for_merge(Genomes, nbreag, distmem, treemem, temp_list->spec,
									 temp_list->merge_with, nb_spec, *spec_left, verbose);
					if (*reduction_type>*spec_left-1) {
			        // there was a merge, need to update
						*reduction_type = *spec_left-1;
					}
				}
			}
		}
	}
	
	erase_list_reag(the_list);
	
	return found_something;
}

//[1.3.2.1] from rearrange.c
int find_k_reag(list_reag **a_list, G_struct *Genomes,
				int depth, int nb_spec, int *spec_left, int reduction_type,
				int *last_genomenb,mgr_distmem_t *distmem, int *list_score,
				int reversal, int transloc, int fusion, int fission,
				int pair1, int pair2, int verbose) // usually both equal to -1 but use them if you want to merge two nodes
{
	struct mgr_genome_struct old_genome;
	struct mgr_genome_struct *genome_list = Genomes->genome_list;
	list_reag *the_list=NULL, *best_reag=NULL;
	list_reag *temp_list;
	int total_distance, next_distance, next_best_distance = MAXINT, 
		new_distance, min_distance = MAXINT;
	int best_list_score=0;
	int did_merge, potential_merge;
	int i, j, k;
  
	total_distance = find_total_dist(Genomes, nb_spec, distmem) ;

#if 0
	fprintf(stdout,"in find %d reag (reduction %d)\n", depth, reduction_type);
#endif
	
	potential_merge = FALSE;
	for (j=0;j<nb_spec;j++) {
		for (k=0;k<nb_spec;k++) {
			if (Genomes->dist_mat[j*nb_spec+k] == 1) {
				potential_merge = TRUE;
			}
		}
	}
	
	if (total_distance>0) {

		// allocate memory for the genomes we use to try various rearragements
		alloc_simple_genome(genome_list[0].num_g, &old_genome);
		
		*list_score = build_list_reag(&the_list, Genomes, nb_spec, *spec_left,
									  reversal, transloc, fusion, fission, 
									  distmem, FALSE, // FALSE because we want the list, not just it's size	
									  FALSE, reduction_type, -1, *last_genomenb, pair1, pair2, verbose);
		
	
#if 0
		if (depth == 2) {
			fprintf(stdout, "list score is %d\n", *list_score);
			print_list_reag(the_list, nb_spec);
		}
#endif
		
    if (the_list != NULL) { // we did find at least one rearrangement
	
		temp_list=the_list;
		while (temp_list!=NULL) {
			list_reag *new_list = NULL;
	
#if 0
			if (depth==2) {
				print_one_reag(temp_list, Genomes, nb_spec, &old_genome);
			}
#endif
		
#if 0
			fprintf(stdout,"checking in %d type %d (%d %d) (%d %d)\n",
					(*temp_list).spec, (*temp_list).optype,
					genome_list[temp_list->spec].genes[(*temp_list).start], 
					genome_list[temp_list->spec].genes[(*temp_list).end],  
					(*temp_list).sc1, (*temp_list).sc2);
#endif


			
			// carry-on the rearrangement
			copy_genes(genome_list[temp_list->spec].genes, old_genome.genes, genome_list[temp_list->spec].num_g);
			/*if (depth==2) {
				old_dist_mouse = (Genomes->dist_mat)[(temp_list->spec)*nb_spec + 2]; // TEMPORARY: old distance with mouse
			}*/
			carry_on_reag(temp_list, Genomes, nb_spec, spec_left, &did_merge, distmem);
	
			next_distance = find_total_dist(Genomes, nb_spec, distmem) ;
			if (depth>1) { //new distance is the distance possible after the current move
				new_distance = find_k_reag(&new_list, Genomes, depth-1,
										   nb_spec, spec_left, reduction_type,
										   last_genomenb,
										   distmem, list_score,
										   reversal, transloc, fusion, fission,
										   pair1, pair2, verbose);
			}
			else {
				new_distance = next_distance;
			}
	
#if 0
			fprintf(stdout, "next_distance %d new_distance %d, did_merge %d\n", next_distance, new_distance, did_merge);
#endif

			if (new_distance<min_distance || 
				((new_distance==min_distance) && ((next_distance<next_best_distance) || (*list_score>best_list_score)))) { 
				// distance is better or the first step is better than previous best

#if 0
					fprintf(stdout, "this is a new best (min_dist %d, list_score %d)\n", new_distance, *list_score);
#endif

				if (min_distance == MAXINT) {
					// it's the first one
					best_reag = list_reag_insert(best_reag, temp_list->spec, temp_list->optype, 
												 temp_list->start, temp_list->end,
												 temp_list->sc1, temp_list->sc2, temp_list->merge_with,
												 temp_list->reduction, temp_list->changes, nb_spec);
				}
				else { // we erase previous best
					best_reag->spec = temp_list->spec;
					best_reag->optype = temp_list->optype;
					best_reag->start = temp_list->start;
					best_reag->end = temp_list->end;
					best_reag->sc1 = temp_list->sc1;
					best_reag->sc2 = temp_list->sc2;
					best_reag->merge_with = temp_list->merge_with;
					best_reag->reduction = temp_list->reduction;
					best_reag->Rscore = temp_list->Rscore;
					for (i=0; i<nb_spec; i++) {
						best_reag->changes[i] = temp_list->changes[i];
					}
	      
					erase_list_reag(best_reag->next);
	      
				}
				best_reag->next = new_list; // because best_reag has the best rearr followed by list of potential rearr.
				
				next_best_distance = next_distance;
				min_distance = new_distance;
				best_list_score = *list_score;
			}
			else {
				erase_list_reag(new_list);
			}
			//}
	  
			// in any case undo last reversals to check for others...
			copy_genes(old_genome.genes, Genomes->genome_list[temp_list->spec].genes, Genomes->genome_list[temp_list->spec].num_g);
			undo_reag(temp_list, Genomes, nb_spec, spec_left, did_merge, distmem);
			
			temp_list=temp_list->next;	
			if (best_reag != NULL && reduction_type == 0 && new_distance<total_distance &&
				potential_merge == FALSE && depth == 1) {
				temp_list = NULL; // shortcut, at depth one, take any rearrangement reducing the distance
			}
		}
		
    }
    erase_list_reag(the_list);	
    
    *a_list = best_reag;

    free_simple_genome(&old_genome);
    
	} else {
		min_distance = 0;
		*a_list = NULL;
	}
  
#if 0
fprintf(stdout,"leaving find %d reag (min_distance %d)\n", depth, min_distance);
#endif
  
fflush(stdout);
  
return min_distance;
}

//[1.3.2] from rearrange.c
int do_k_reag(G_struct *Genomes,
			  int depth, int *nbreag,
			  int nb_spec, int *spec_left,					
			  mgr_distmem_t *distmem, treemem_t *treemem, 
			  int *last_genomenb, int verbose,
			  int pair1, int pair2, //usually both equal to -1 but if we want to merge, use them...
			  struct mgr_genome_struct *old_genome) 
{
  list_reag *the_list = NULL;
  int total, new_total;
  int found_something = FALSE;
  int did_merge;
  int num_chromosomes = Genomes->genome_list[0].num_chr;
  
  
  if (verbose) {
	  fprintf(stdout,"in do %d reag (spec_left %d)\n", depth, *spec_left);
  }
  total = find_total_dist(Genomes, nb_spec, distmem) ;
  
  
  if (total>0) {
    int list_score=0;

    // look for rearrangements (rev/transloc only) reducing the distance by at least 1 

	  new_total = find_k_reag(&the_list, Genomes, depth,		 
							  nb_spec, spec_left, 1, last_genomenb,
							  distmem, &list_score,
							  TRUE, (num_chromosomes>0)?TRUE:FALSE, FALSE, FALSE,
							  pair1, pair2, verbose);

    // look for any rearrangements reducing the distance by at least 1
	  if (total-new_total<3 && num_chromosomes>0) {  // not good enough
		  erase_list_reag(the_list);
		  if (verbose)
			  fprintf(stdout, "look including fusion/fission\n");
		  new_total = find_k_reag(&the_list, Genomes, depth,		 
								  nb_spec, spec_left, 1, last_genomenb,
								  distmem, &list_score,
								  TRUE, (num_chromosomes>0)?TRUE:FALSE, 
								  (num_chromosomes>0)?TRUE:FALSE, (num_chromosomes>0)?TRUE:FALSE,
								  pair1, pair2, verbose);
	  } 

	  
    
    // otherwise consider any rearrangements not making situation worst
    if (the_list==NULL) {
      new_total = find_k_reag(&the_list, Genomes, depth,
							  nb_spec, spec_left, 0, last_genomenb,
							  distmem, &list_score,
							  TRUE, (num_chromosomes>0)?TRUE:FALSE, 
							  (num_chromosomes>0)?TRUE:FALSE, 
							  (num_chromosomes>0)?TRUE:FALSE,
							  pair1, pair2, verbose);
	}
	
	  
#if 0
	  if (verbose) {
		  print_list_reag(the_list, nb_spec);
		  fprintf(stdout,"total was %d and now is %d\n", total, new_total);
	  }
#endif

//	if (total-new_total>=0) {
    if (total-new_total>0) {
      
		if (the_list!=NULL) {
			if (verbose)
				print_one_reag(the_list, Genomes, nb_spec, old_genome);
			
			
			carry_on_reag(the_list, Genomes, nb_spec,  
						  spec_left, &did_merge, distmem);

			if (verbose) {
				special_print_genomes2(stdout, Genomes, nb_spec);
			}
			
			if (the_list->merge_with != -1) /* there was a merge and we need to update tree */
				update_for_merge(Genomes, nbreag, distmem, treemem, the_list->spec, the_list->merge_with,
								 nb_spec, *spec_left, verbose);
			
			nbreag[the_list->spec]++;
			*last_genomenb = the_list->spec; 
			found_something = TRUE ;
		}
	}    
	
	  
erase_list_reag(the_list);
	}
	
#ifdef DEBUGM
  fprintf(stdout,"leaving do_k_reag\n");	
#endif
	
  return found_something;
}



//[1.3]from rearrange.c
void solve_with_good_reag(G_struct *Genomes, int *nbreag,
						  int nb_spec, int *spec_left, 
						  mgr_distmem_t *distmem, treemem_t *treemem, 
						  int heuristic, int pair1, int pair2, int verbose) 
{

  int found_something = TRUE;
  int last_genomenb = 1;
  int reduction_type;
  struct mgr_genome_struct old_genome;

  // allocate memory for the genomes we use to try various rearragements
  alloc_simple_genome(Genomes->genome_list[0].num_g, &old_genome);

  while (found_something == TRUE &&
		 find_total_dist(Genomes, nb_spec, distmem)!=0) {
    
    reduction_type = *spec_left-1;
    
    // start with reversals
	if (verbose) {
		fprintf(stdout,"look for a reversal\n");		
		fflush(stdout);
	}

    found_something = carry_best_reag(Genomes , nbreag, nb_spec,
									  spec_left,
									  TRUE, 	// TRUE because we want reversals
									  FALSE, 	// no translocations
									  FALSE,
									  FALSE, 	// FALSE because we don't want to include fusions and fissions
									  distmem, treemem, &reduction_type, &last_genomenb,
									  heuristic, verbose,
									  pair1, pair2,
									  &old_genome);

    // otherwise look for transloc
    if (found_something == FALSE && Genomes->genome_list[0].num_chr>0) {
		if (verbose) {
			fprintf(stdout,"look for a transloc\n");		
			fflush(stdout);
		}

		found_something = carry_best_reag(Genomes , nbreag, nb_spec,
										  spec_left,
										  FALSE, 	// FALSE because we don't want reversals
										  TRUE, 	// FALSE because we don't want translocations
										  FALSE,
										  FALSE, 	// TRUE because we want to include fusions and fissions
										  distmem, treemem, &reduction_type, &last_genomenb, 
										  heuristic, verbose,
										  pair1, pair2,
										  &old_genome);
    }

    
    // otherwise look for fusion/fission
    if (found_something == FALSE && Genomes->genome_list[0].num_chr>0) {
		
		if (verbose) {
			fprintf(stdout,"look for a fusion/fission\n");		
			fflush(stdout);
		}
		
		found_something = carry_best_reag(Genomes , nbreag, nb_spec,
										  spec_left,
										  FALSE, 	// FALSE because we don't want reversals
										  FALSE, 	// FALSE because we don't want translocations
										  TRUE,
										  TRUE, 	// TRUE because we want to include fusions and fissions
										  distmem, treemem, &reduction_type, &last_genomenb, heuristic, verbose,
										  pair1, pair2,
										  &old_genome);
    }
		 
    while (found_something == FALSE && reduction_type>0) {

      // can't find a good rearrangement, need to relax condition
		reduction_type--;
      
       
		if (verbose) {
			fprintf(stdout,"look for a reduction of %d\n", reduction_type);		
			fflush(stdout);
		}
      
          
		if (reduction_type==0) {
   
      	  // do a depth 2 search to find next rearrangement
			if (do_k_reag(Genomes, 2, nbreag, 
						  nb_spec, spec_left, distmem, treemem, 
						  &last_genomenb, verbose, pair1, pair2,
						  &old_genome) == TRUE) {	    
				found_something = TRUE;
			} else {
				fprintf(stdout, "couldn't find depth 2 rearrangement in solve_with_good_reag\n");
	      //exit(-1);
			}
		} else {
	  
      // start with reversals
			if (verbose) {
				fprintf(stdout,"look for a reversal\n");		
				fflush(stdout);
			}

			found_something = carry_best_reag(Genomes , nbreag, nb_spec,
											  spec_left,
											  TRUE, 	// TRUE because we want reversals
											  FALSE, 	// no translocations
											  FALSE,
											  FALSE, 	// FALSE because we don't want to include fusions and fissions
											  distmem, treemem, &reduction_type, &last_genomenb, heuristic, verbose,
											  pair1, pair2,
											  &old_genome);
			
      // otherwise look for transloc
			if (found_something == FALSE && Genomes->genome_list[0].num_chr>0) {
				if (verbose) {
					fprintf(stdout,"look for a transloc\n");		
					fflush(stdout);
				}

				found_something = carry_best_reag(Genomes , nbreag, nb_spec,
												  spec_left,
												  FALSE, 	// FALSE because we don't want reversals
												  TRUE, 	// FALSE because we don't want translocations
												  FALSE,
												  FALSE, 	// TRUE because we want to include fusions and fissions
												  distmem, treemem, &reduction_type, &last_genomenb, 
												  heuristic, verbose,
												  pair1, pair2,
												  &old_genome);
			}

    
      // otherwise look for fusion/fission
			if (found_something == FALSE && Genomes->genome_list[0].num_chr>0) {
				
				if (verbose) {
					fprintf(stdout,"look for a fusion/fission\n");		
					fflush(stdout);
				}

				found_something = carry_best_reag(Genomes , nbreag, nb_spec,
												  spec_left,
												  FALSE, 	// FALSE because we don't want reversals
												  FALSE, 	// FALSE because we don't want translocations
												  TRUE,
												  TRUE, 	// TRUE because we want to include fusions and fissions
												  distmem, treemem, &reduction_type, &last_genomenb, 
												  heuristic, verbose,
												  pair1, pair2,
												  &old_genome);
			}
		}  
    }
    
  }
  
  
  free(old_genome.genes);

}

/* initialize the treemem structure that contains the phylogeny */
void initialize_treemem(treemem_t *treemem, int nb_spec) {
	 
	 int i;
	 // the nodes already in the tree, 1 implies in the tree, 
	 //								   0 implies not in tree
	 treemem->tree_array = (int *)malloc(nb_spec * sizeof(int));
	 for (i=0;i<nb_spec;i++)
		  treemem->tree_array[i]=0;
	 treemem->tree_size = 0;
	 treemem->new_node = nb_spec;
	 treemem->merge_edge_list = NULL;
	 treemem->the_edge_list = NULL;
}

/* free the memory of the phylogeny */
void free_treemem(treemem_t *treemem) {
	 
	 free(treemem->tree_array);
	 erase_list_edge(treemem->merge_edge_list);
	 erase_list_edge(treemem->the_edge_list);
}

void free_mem_4_mgr()
{
    free_G_struct(&Genomes,4); 
    free_G_struct(&Genomes_copy,3); 
    mcdist_freemem(&MGR_DISTMEM);
}

void mgr_ini_mem (int num_genes, int num_chromosomes)
{
    int num_genomes = 3; //median of 3 genomes
    mcdist_allocmem(num_genes, num_chromosomes, &MGR_DISTMEM);
    mgr_genome_list_copy =
    (struct mgr_genome_struct *) malloc(num_genomes*sizeof(struct mgr_genome_struct));
    int i;
    for (i=0 ; i<num_genomes; i++) {
    (mgr_genome_list_copy)[i].gnamePtr = (char *) malloc(MAX_STR_LEN*sizeof(char));
    (mgr_genome_list_copy)[i].genes = (int*)malloc(num_genes*sizeof(int));
    }
    mgr_genome_list = (struct mgr_genome_struct *) malloc(2*(num_genomes-1)*sizeof(struct mgr_genome_struct));
    for (i=0 ; i<2*(num_genomes-1); i++) {
    (mgr_genome_list)[i].gnamePtr = (char *) malloc(MAX_STR_LEN*sizeof(char));
    (mgr_genome_list)[i].genes = (int*)malloc(num_genes*sizeof(int));
    }

    int nb_spec = num_genomes;
    Genomes.label = (int *)malloc(nb_spec * sizeof(int));
    Genomes.same_as = (int *)malloc(nb_spec * sizeof(int));
    Genomes.nb_chromo = (int *)malloc(nb_spec * sizeof(int));
    Genomes.dist_mat = (int *)malloc(nb_spec * nb_spec * sizeof(int));

    Genomes_copy.label = (int *)malloc(nb_spec * sizeof(int));
    Genomes_copy.same_as = (int *)malloc(nb_spec * sizeof(int));
    Genomes_copy.nb_chromo = (int *)malloc(nb_spec * sizeof(int));
    Genomes_copy.dist_mat = (int *)malloc(nb_spec * nb_spec * sizeof(int));

}

void mgr_med (int * g1, int * g2,int * g3, int SIZE_ALPHA, int CIRCULAR, struct genome_struct * g_med)
{
    int i,j;
    int NumGenomes = 3;
    int circular = CIRCULAR;
    int nb_spec=3;
    int num_genes;
    int size_alpha = SIZE_ALPHA;
    int num_chromosomes;
    int genome_type=GLINEAR;
    int condensing = FALSE;
    int optimize = FALSE;
	int alternatives = FALSE;
    int verbose = FALSE;
    int depth=2;/* depth of the search when running out of good rearrangements */
    int heuristic = 0;
    // -H: heuristic to speed up triplet resolution
    //    1 -> only look at reversals initially and picks first good one
    //    2 -> only look at reversals initially, take shortest one

    /* Data structure for the phylogeny */
    treemem_t treemem;	/* the tree itself */
    int spec_left;		/* number of genomes left to be put in the tree */
    int *nbreag = NULL;/* number of rearrangements that were carried in each genome */
    int total_weight;	/* total weight of the tree */
    double avg_nb_rev=0;
    
    /* consider unichromosome first
    if (genome_type == GMULTI) { // multichromosomal distance
    // we've assume that we have the maximum number of chromosomes
    // in order not to reassign memory
       num_genes = 3*size_alpha;
       num_chromosomes = size_alpha;
    }
    else */ 
    { // unichromosomal distance
                num_genes = size_alpha;
                num_chromosomes = 0;
    }

    
   for (j=0;j<size_alpha;j++)
   {
    (mgr_genome_list_copy)[0].genes[j] = g1[j]; 
    (mgr_genome_list_copy)[1].genes[j] = g2[j]; 
    (mgr_genome_list_copy)[2].genes[j] = g3[j];
   }
  
    // finish initialization of the input genomes, including actual nb_chromo
    init_G_struct
        (&Genomes_copy, mgr_genome_list_copy, nb_spec, num_genes, num_chromosomes);

  
    // assign memory for genomes assuming each marker is his own chromosome (worst case)
    initialize_genome_list(&mgr_genome_list, 2*(nb_spec-1), num_genes, num_chromosomes);
    // initialize structure for Genomes which includes true nb_chromo for each genome
    init_G_struct(&Genomes, mgr_genome_list, 2*(nb_spec-1), num_genes, num_chromosomes);
    // we also keep the genome_list_copy as a reference
    copy_first_n(&Genomes_copy, &Genomes, nb_spec);
    for (i=0; i<nb_spec; i++)
        strcpy((mgr_genome_list)[i].gnamePtr, (mgr_genome_list_copy)[i].gnamePtr);

    // corresponds to the true size of the alphabet
    size_alpha = num_genes - 2*num_chromosomes;


    if (condensing == TRUE) {
        condense_genomes(&Genomes, nb_spec, verbose);
		//uncondense_genomes(&Genomes, nb_spec);
		//special_print_genomes2(stdout, &Genomes, nb_spec);
	}

    compute_dist_mat(&Genomes, nb_spec, &MGR_DISTMEM);

     // allocate memory for the tree/phylogeny that we want to build
    initialize_treemem(&treemem, nb_spec);

    // the number of species left to place on the tree
	spec_left=nb_spec;
			
    // the number of transformations we did on each genome
    nbreag = (int *)malloc(nb_spec * sizeof(int));
    for (i=0;i<nb_spec;i++)
          nbreag[i]=0;

    // check to see if two genomes are identical before we begin
	check_for_merge(&Genomes, nbreag, nb_spec, &spec_left,
					&MGR_DISTMEM, &treemem, verbose);

    if (spec_left == 3) {
        int perfect_triplet;
        if (verbose) {
        solve_triplet(&Genomes, nbreag, &perfect_triplet, 2,
                          nb_spec, &spec_left,
                          depth, &MGR_DISTMEM, &treemem, &Preancestors, 
                          heuristic, condensing, verbose);
        } else {
            solve_triplet(&Genomes, nbreag, &perfect_triplet, 2,
                          nb_spec, &spec_left,
                          depth, &MGR_DISTMEM, &treemem, NULL, 
                          heuristic, condensing, verbose);
        }
    }
/* debug msg
    struct mgr_genome_struct *tmplist= Genomes_copy.genome_list;
    tmplist = Genomes.genome_list;
    for(ii=0;ii<4;ii++)
    {
                fprintf(stdout,"ii=%d : [",ii);
                int * point = tmplist->genes;
                for(jj=0;jj<size_alpha;jj++)
                    fprintf(stdout,"%d,",point[jj]);
                fprintf(stdout,"]\n");
                fflush(stdout);
                tmplist ++;
    }
debug msg */
    /* do we need this ?
  // only need to do this if we haven't been able to solve the problem
    // using only good reversals
    if (find_total_dist(&Genomes, nb_spec, &MGR_DISTMEM) !=0) {
        int did_grow = TRUE;
        G_struct Star;
        struct genome_struct *star_list;
        if (verbose) {
            fprintf(stdout, "\nStart building phylogeny...\n");
            fflush(stdout);
        }
        if (condensing == TRUE) {
            condense_genomes(&Genomes, nb_spec, verbose);
            //uncondense_genomes(&Genomes, nb_spec);
        }
        // memory and initialization needed in start_tree and growing_tree
        initialize_genome_list(&star_list, 3, num_genes, num_chromosomes);
        init_G_struct(&Star, star_list, 3, num_genes, num_chromosomes);
        start_tree(&Genomes, &Star, nbreag, nb_spec, &spec_left,
                   depth, &MGR_DISTMEM, &treemem, heuristic, condensing, verbose) ;
        // at each iteration, adds a new genome to the tree
        while (treemem.tree_size < nb_spec && did_grow==1) {
            if (verbose) {
                fprintf(stdout, "Building phylogeny, %d %s remaining to place\n",
                        nb_spec-treemem.tree_size, (nb_spec-treemem.tree_size>1)?"genomes":"genome");
                fflush(stdout);
            }
            growing_tree(&Genomes, &Star, nbreag, nb_spec, &spec_left,
                         depth, &MGR_DISTMEM, &treemem, &did_grow, heuristic, condensing, verbose) ;
            if (verbose) {
                special_print_genomes2(stdout, &Genomes, 2*(nb_spec-1));
            }
        }
        free_G_struct(&Star, 3);
        // these good rearrangements weren't counted in the edges yet
        add_nbrev_edges(treemem.the_edge_list, nbreag, nb_spec);
    }// if (find_total_dist(&Genomes, nb_spec, &MGR_DISTMEM) !=0)
    else { // when the tree was solve only using good rearrangements
        // need to finish the list with the last 3 genomes
        fprintf(stdout,"finish edge list\n");
        fflush(stdout);
        finish_edge_list(&Genomes, nbreag, nb_spec, &treemem);
    }// if (find_total_dist(&Genomes, nb_spec, &MGR_DISTMEM) !=0)
    */
    // merge the list of merged that occured while doing good rearrangements
    // with the list of edges that occured when building the phylogeny
   //  merge_two_list(&treemem); 
    // uncondense
    if (condensing == TRUE) {
        if (verbose) {
            special_print_genomes2(stdout, &Genomes,2*(nb_spec-1));
        }												
        uncondense_genomes(&Genomes, 2*(nb_spec-1), verbose);	
    }
        
    copy_genes((Genomes.genome_list)->genes,g_med->genes,size_alpha);

    // free memory for the number of rearrangements
    free(nbreag);
    // free memory for phylogeny
    free_treemem(&treemem);
   
    // free preancestors
    if (verbose) {
        free_G_struct(&Preancestors, nb_spec);
    }
    
    

    
}


//[1]from rearrange.c
int solve_triplet (G_struct *Genomes, int *nbreag, int *perfect_triplet, int reduction_type, int nb_spec, int *spec_left, int depth, mgr_distmem_t *distmem, treemem_t *treemem,G_struct *Preancestors,int heuristic, int condensing, int verbose) 
{
	
	int current_total_dist, i, total_nb_reag;
	int last_genomenb = 1; // actually this initialization should probably be done only once (before first call)
	struct mgr_genome_struct old_genome;
	int found_something, first_pass;
	
	// allocate memory for the genomes we use to try various rearragements
	alloc_simple_genome(Genomes->genome_list[0].num_g, &old_genome);
		
	
	current_total_dist = find_total_dist(Genomes, nb_spec, distmem);
	// fprintf(stdout, "current_total: %d\n", current_total_dist);

    //add the condensing later
	if (condensing == TRUE && current_total_dist != 0)
        condense_genomes(Genomes, nb_spec, verbose);

	found_something = first_pass = TRUE;
	while (current_total_dist!=0 && found_something == TRUE) {
	  	/* start with all good rearrangements */
		solve_with_good_reag(Genomes, nbreag, 
							 nb_spec, spec_left,
							 distmem, treemem, heuristic, -1, -1, verbose);
		current_total_dist = find_total_dist(Genomes, nb_spec, distmem);
		if (first_pass == TRUE) {
			if (condensing == TRUE && current_total_dist != 0) 
				condense_genomes(Genomes, nb_spec, verbose);
			if (Preancestors!=NULL) {
     	    // keep track of preancestors
				copy_first_n(Genomes, Preancestors,  nb_spec);
				for (i=0; i<nb_spec; i++)
					sprintf((Preancestors->genome_list[i]).gnamePtr, "%s*", (Genomes->genome_list[i]).gnamePtr);
			}
			first_pass = FALSE;
		}
		if (current_total_dist != 0) {
			*perfect_triplet = FALSE;
        	  // do a depth k search to find next rearrangement
			if (do_k_reag(Genomes, depth, nbreag,
						  nb_spec, spec_left, distmem, treemem, 
						  &last_genomenb, verbose, -1, -1,
						  &old_genome) == TRUE) 
			{ 
			//we found one restart normally
				found_something = TRUE;
			} else {
				if (depth < 3) {
					fprintf(stdout,"we are increasing depth of search to %d\n", ++depth);
					if (do_k_reag(Genomes, depth, nbreag, 
								  nb_spec, spec_left, distmem, treemem, 
								  &last_genomenb, verbose, -1, -1,
								  &old_genome) == TRUE) { //restart
						
						found_something = TRUE;
					}
					else {
						found_something = FALSE;
						fprintf(stdout,"couldn't solve this triplet\n");
						exit(-1);
					}
				}
			}
		}
	}
	if (condensing == TRUE) {
		uncondense_genomes(Genomes, nb_spec, verbose);
	}
	total_nb_reag = 0; 
	for (i=0;i<nb_spec;i++)
		total_nb_reag += nbreag[i];
	free(old_genome.genes);
	return total_nb_reag;	
}


