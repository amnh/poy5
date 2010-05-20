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
#include <caml/fail.h>
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
if (verbose) {
    fprintf(stdout,"mgr.c carry best reag,\n "); fflush(stdout); }
    struct mgr_genome_struct *genome_list = Genomes->genome_list;
    //list_reag *the_list = NULL;
    int score_list, total_dist, found_something = FALSE;
	int did_merge;
    list_reag *temp_list = NULL;
    
    list_listreag * the_list = &list_carry_best_reag;
    list_listreag * new_list = &list_carry_best_reag2;
    
    if((the_list==NULL)||(the_list->reaglist == NULL))
    {
		if (verbose) {
            fprintf(stdout, "init list_carry_best_reag \n"); fflush(stdout); }
       init_list_listreag_memory(the_list, nb_spec, 100);
    } else {}
    if((new_list==NULL)||(new_list->reaglist == NULL))
    {
		if (verbose) {
         fprintf(stdout, "init list_carry_best_reag2 \n"); fflush(stdout); }
       init_list_listreag_memory(new_list, nb_spec, 100);
    } else {}

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
		if (total_dist == 0)	return FALSE;
		if (verbose) {
           fprintf(stdout,"call build_list_reag 1\n"); fflush(stdout); }
        // build the list of good rearrangements
        the_list->list_size = 0;
		score_list = build_list_reag(the_list, Genomes, nb_spec, *spec_left,
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
		if (the_list->list_size >0 ) {// we did find a good rearrangement
            //(the_list != NULL) { // we did find a good rearrangement
			int max_reduction=-1, new_reduction;
			int max_index = -1, the_index=1;
            int idx=0;
			//temp_list=the_list; 
			//while (temp_list!=NULL) {
		   if (verbose) {
            fprintf(stdout,"result of build_list_reag,list_size=%d,start for loop \n",
                    the_list->list_size); fflush(stdout);}
            for(idx=0;idx<the_list->list_size;idx++) {
                temp_list = &(the_list->reaglist[idx]);
				// list_reag *new_list = NULL; 
                //this is a bug, memory belongs to new_list is never freed.
				// carry-on the rearrangement with chromo flip if necessary
				copy_genes(genome_list[temp_list->spec].genes, old_genome->genes, genome_list[temp_list->spec].num_g);
				carry_on_reag(temp_list, Genomes, nb_spec, 
							  spec_left, &did_merge, distmem);
				new_reduction = build_list_reag(new_list, Genomes, nb_spec, *spec_left,
												TRUE, (genome_list[0].num_chr>0)?TRUE:FALSE,
												(genome_list[0].num_chr>0)?TRUE:FALSE,
                                                (genome_list[0].num_chr>0)?TRUE:FALSE,
												/* because we want to know the total number of
												* potential rearragements */
												distmem, TRUE, // TRUE because we only want size
												FALSE, *reduction_type, -1, temp_list->spec,
                                                pair1, pair2, verbose);
				if (new_reduction>max_reduction) {
					//this one is better
					max_reduction = new_reduction;
					max_index = the_index;
				}
				// in any case undo last reversals to check for others...
				copy_genes(old_genome->genes, Genomes->genome_list[temp_list->spec].genes, 
                        Genomes->genome_list[temp_list->spec].num_g);
				undo_reag(temp_list, Genomes, nb_spec, spec_left, did_merge, distmem);
				the_index++;
				//temp_list=temp_list->next;
			}
			if (max_index!=-1) { // we found the best rearrangement
				int i;
			    //temp_list = the_list;
                // we go back to the best reversal we found
                temp_list = &(the_list->reaglist[max_index-1]);
			    /*	for (i=1;i<max_index;i++) { temp_list=(*temp_list).next;} */
				carry_on_reag(temp_list, Genomes, nb_spec,							
							  spec_left, &did_merge, distmem);
				nbreag[temp_list->spec]++;
				*last_genomenb = temp_list->spec;
				found_something = TRUE ;
				//potentially update max_chromo_size
				find_max_chromo_size(Genomes, nb_spec, cb_max_chromo_size);
				if (temp_list->merge_with != -1) {// there was a merge and we need to update tree
					update_for_merge(Genomes, nbreag, distmem, treemem, temp_list->spec,
									 temp_list->merge_with, nb_spec, *spec_left, verbose,
                                     cb_max_chromo_size);
					if (*reduction_type>*spec_left-1) {
			        // there was a merge, need to update
						*reduction_type = *spec_left-1;
					}
				}
			}
		}
	}
//	erase_list_reag(the_list);
	if (verbose) {
    fprintf(stdout, "end of carry best reag\n"); fflush(stdout); }
	return found_something;
}

//[1.3.2.1] from rearrange.c 
//NOTE: this function is recursive
int find_k_reag(list_reag **a_list, G_struct *Genomes,
				int depth, int nb_spec, int *spec_left, int reduction_type,
				int *last_genomenb, mgr_distmem_t *distmem, int *list_score,
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
	fprintf(outfile,"in find %d reag (reduction %d)\n", depth, reduction_type);
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
		
        list_listreag reaglist;
        list_listreag * the_list = NULL;
        the_list = &reaglist;
        if((the_list==NULL)||(the_list->reaglist == NULL))
        {
		if (verbose) {
            fprintf(stdout, "init listreag in find_k_reag \n"); fflush(stdout); }
            init_list_listreag_memory(the_list, nb_spec, 100);
        } else {}
		*list_score = build_list_reag(the_list, Genomes, nb_spec, *spec_left,
									  reversal, transloc, fusion, fission, 
									  distmem, FALSE, // FALSE because we want the list, not just it's size	
									  FALSE, reduction_type, -1, *last_genomenb, pair1, pair2, verbose);
		
	
#if 0
		if (depth == 2) {
			fprintf(outfile, "list score is %d\n", *list_score);
			print_list_reag(the_list, nb_spec);
		}
#endif
		
    //if (the_list != NULL) { // we did find at least one rearrangement
    if (the_list->list_size > 0) { // we did find at least one rearrangement
	   //temp_list=the_list;
        int idx = 0;   
		if (verbose) {
        fprintf(stdout,"the_list->list_size = %d, start while loop:\n",the_list->list_size);} 
        fflush(stdout);
		//while (temp_list!=NULL) {
		while (idx<the_list->list_size) {
		    if (verbose) {
            fprintf(stdout, "idx=%d,",idx);fflush(stdout); }
            temp_list = &(the_list->reaglist[idx]);
			list_reag *new_list = NULL;
	
#if 0
			if (depth==2) {
				print_one_reag(temp_list, Genomes, nb_spec, &old_genome);
			}
#endif
		
#if 0
			fprintf(outfile,"checking in %d type %d (%d %d) (%d %d)\n",
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
			fprintf(outfile, "next_distance %d new_distance %d, did_merge %d\n", next_distance, new_distance, did_merge);
#endif
			if (new_distance<min_distance || 
				((new_distance==min_distance) && ((next_distance<next_best_distance) || (*list_score>best_list_score)))) { 
				// distance is better or the first step is better than previous best
#if 0
					fprintf(outfile, "this is a new best (min_dist %d, list_score %d)\n", 
                            new_distance, *list_score);
#endif
				if (min_distance == MAXINT) {
					// it's the first one
					best_reag = list_reag_insert
                        (best_reag, temp_list->spec, temp_list->optype, 
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
			}// if (new_distance<min_distance....
			else {
				erase_list_reag(new_list);
			}//end of if (new_distance<min_distance.....
			// in any case undo last reversals to check for others...
			copy_genes(old_genome.genes, Genomes->genome_list[temp_list->spec].genes, Genomes->genome_list[temp_list->spec].num_g);
			undo_reag(temp_list, Genomes, nb_spec, spec_left, did_merge, distmem);
			
			idx ++;//temp_list=temp_list->next;	
			
            if (best_reag != NULL && reduction_type == 0 && new_distance<total_distance &&
				potential_merge == FALSE && depth == 1) {
				idx=the_list->list_size;//temp_list = NULL; 
                // shortcut, at depth one, take any rearrangement reducing the distance
			}
		}//end of 	while (idx<the_list->list_size) ...
		if (verbose) fprintf(stdout, "end of while \n");
    }//end of if (the_list->list_size > 0)
     free_list_listreag_memory(the_list);//erase_list_reag(the_list);	
    
    *a_list = best_reag;

    free_simple_genome(&old_genome);
    
	}//	if (total_distance>0)  
    else {
		min_distance = 0;
		*a_list = NULL;
	}//  end of if (total_distance>0) 
  
#if 0
fprintf(outfile,"leaving find %d reag (min_distance %d)\n", depth, min_distance);
#endif
  
//fflush(outfile);
  
return min_distance;
}


int find_k_reag_withmem(list_listreag * the_list, //list_reag **a_list,
                G_struct *Genomes,
				int depth, int nb_spec, int *spec_left, int reduction_type,
				int *last_genomenb,mgr_distmem_t *distmem, int *list_score,
				int reversal, int transloc, int fusion, int fission,
				int pair1, int pair2, int verbose) // usually both equal to -1 but use them if you want to merge two nodes
{
    //remember where we start at the list_listreag
    int start_list_size = the_list->list_size; 
    int copysize = 0;
	struct mgr_genome_struct old_genome;
	struct mgr_genome_struct *genome_list = Genomes->genome_list;
    list_reag *best_reag=NULL;
	list_reag *temp_list;
	int total_distance, next_distance, next_best_distance = MAXINT, 
		new_distance, min_distance = MAXINT;
	int best_list_score=0;
	int did_merge, potential_merge;
	int i, j, k;
	total_distance = find_total_dist(Genomes, nb_spec, distmem) ;
	
	if (verbose) {
    fprintf(stdout,"\nin find %d reag (reduction %d),list start at:%d\n",
            depth, reduction_type,start_list_size);
    fflush(stdout); }
	
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
		*list_score = build_list_reag(the_list,Genomes, nb_spec, *spec_left,
									  reversal, transloc, fusion, fission, 
									  distmem, FALSE, // FALSE because we want the list, not just it's size	
									  FALSE, reduction_type, -1, *last_genomenb, pair1, pair2, verbose);
    int idx_after_build_list_reag = the_list->list_size;
	if (verbose) {
    fprintf(stdout, "list size after build_list_reag: %d\n", idx_after_build_list_reag ); 
    fflush(stdout); }
    if (idx_after_build_list_reag > start_list_size) { 
        // we did find at least one rearrangement
       int idx = start_list_size;//0;
       int bestreag_idx = the_list->list_size;
       while(idx<idx_after_build_list_reag) {
           temp_list = &(the_list->reaglist[idx]);
			// carry-on the rearrangement
			copy_genes(genome_list[temp_list->spec].genes, old_genome.genes, genome_list[temp_list->spec].num_g);
			//if (depth==2) {
		//		old_dist_mouse = (Genomes->dist_mat)[(temp_list->spec)*nb_spec + 2]; // TEMPORARY: old distance with mouse
		//	}
			carry_on_reag(temp_list, Genomes, nb_spec, spec_left, &did_merge, distmem);
			next_distance = find_total_dist(Genomes, nb_spec, distmem) ;
            the_list->list_size = bestreag_idx+1;
			if (depth>1) { //new distance is the distance possible after the current move
                if (verbose) fprintf(stdout, "depth>1, start dept=%d recursive \n",depth-1);
				new_distance = find_k_reag_withmem(the_list, Genomes, depth-1,
										   nb_spec, spec_left, reduction_type,
										   last_genomenb,
										   distmem, list_score,
										   reversal, transloc, fusion, fission,
										   pair1, pair2, verbose);
		      if (verbose) 
                fprintf(stdout, "list size after dept = %d recursive: %d\n",
                        depth-1,the_list->list_size);
			}
			else { 
                new_distance = next_distance;	}
			if (new_distance<min_distance || 
				((new_distance==min_distance) && ((next_distance<next_best_distance) || (*list_score>best_list_score)))) {
				// distance is better or the first step is better than previous best
                if (verbose) 
				fprintf(stdout, "find a new best (min_dist %d, list_score %d), list_size set back to %d\n", new_distance, *list_score, bestreag_idx);
                the_list->list_size = bestreag_idx;
				if (min_distance == MAXINT) {
                    if (verbose) {
                    fprintf(stdout, "first best_reag, add to list\n");fflush(stdout); }
					// it's the first one
					list_reag_insert2(the_list, temp_list->spec, temp_list->optype, 
												 temp_list->start, temp_list->end,
												 temp_list->sc1, temp_list->sc2, temp_list->merge_with,
												 temp_list->reduction, temp_list->changes, nb_spec);
                }
				else { // we erase previous best
                    if (verbose) {
                    fprintf(stdout, "erase previous best_reag,"); fflush(stdout); }
                    list_reag_insert2(the_list, temp_list->spec, temp_list->optype, 
												 temp_list->start, temp_list->end,
												 temp_list->sc1, temp_list->sc2, temp_list->merge_with,
												 temp_list->reduction, temp_list->changes, nb_spec);

				}
				next_best_distance = next_distance;
				min_distance = new_distance;
				best_list_score = *list_score;
			}
			else {
                the_list->list_size =  bestreag_idx;
            }//end of if (new_distance<min_distance...
			// in any case undo last reversals to check for others...
			copy_genes(old_genome.genes, Genomes->genome_list[temp_list->spec].genes, Genomes->genome_list[temp_list->spec].num_g);
			undo_reag(temp_list, Genomes, nb_spec, spec_left, did_merge, distmem);
            idx ++;
			if (the_list->list_size > bestreag_idx 
                    && reduction_type == 0 && new_distance<total_distance &&
				potential_merge == FALSE && depth == 1) {
                // shortcut, at depth one, take any rearrangement reducing the distance
		        if (verbose) {
                fprintf(stdout, "depth=1, take any reag reduce the dis\n");
                fflush(stdout);}
                idx=idx_after_build_list_reag; //cut the while loop
			}
		    if (verbose) 
                fprintf(stdout, "end of while, list_size = %d\n",the_list->list_size);
           }//end of while(idx<idx_after_build_list_reag) 	
        }//end of if(idx_after_build_list_reag > start_list_size)  
       ///what to do with min_distance? = total_distance?
        free_simple_genome(&old_genome);
	}//if (total_distance>0) 
    else {
		min_distance = 0;
	}//end of if (total_distance>0)
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
  list_listreag *the_list = NULL;
  //list_reag *the_list = NULL;
  
  list_reag *tmp_reag = NULL;
  int total, new_total;
  int found_something = FALSE;
  int did_merge;
  int num_chromosomes = Genomes->genome_list[0].num_chr;
  
  
  if (verbose) {
	  fprintf(stdout,"in do %d reag (spec_left %d)\n", depth, *spec_left);
  }
  total = find_total_dist(Genomes, nb_spec, distmem) ;
  
  the_list = &list_do_k_reag;
  if((the_list==NULL)||(the_list->reaglist == NULL))
  {
  if (verbose) {
       fprintf(stdout, "init list_do_k_reag \n"); fflush(stdout); }
       init_list_listreag_memory(the_list, nb_spec, 100);
  } else {}

  if (total>0) {
    int list_score=0;

    // look for rearrangements (rev/transloc only) reducing the distance by at least 1 

	  new_total = find_k_reag_withmem(the_list, Genomes, depth,		 
	 // new_total = find_k_reag(&the_list, Genomes, depth,		 
							  nb_spec, spec_left, 1, last_genomenb,
							  distmem, &list_score,
							  TRUE, (num_chromosomes>0)?TRUE:FALSE, FALSE, FALSE,
							  pair1, pair2, verbose);

    // look for any rearrangements reducing the distance by at least 1
	  if (total-new_total<3 && num_chromosomes>0) {  // not good enough
		  //erase_list_reag(the_list);
          the_list->list_size = 0;
		  if (verbose)
			  fprintf(stdout, "look including fusion/fission\n");
		  new_total = find_k_reag_withmem(the_list, 
		  //new_total = find_k_reag(&the_list, 
                                  Genomes, depth,		 
								  nb_spec, spec_left, 1, last_genomenb,
								  distmem, &list_score,
								  TRUE, (num_chromosomes>0)?TRUE:FALSE, 
								  (num_chromosomes>0)?TRUE:FALSE, (num_chromosomes>0)?TRUE:FALSE,
								  pair1, pair2, verbose);
	  } 

	  
    
    // otherwise consider any rearrangements not making situation worst
     if(the_list->list_size==0) {
     //if (the_list==NULL) {
    //  new_total = find_k_reag(&the_list,
      new_total = find_k_reag_withmem(the_list,
              Genomes, depth, nb_spec, spec_left, 0, last_genomenb,
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

    if (total-new_total>0) {
      
		tmp_reag = &(the_list->reaglist[0]); 
		if(the_list->list_size>0) {
        //if (the_list!=NULL) {
			if (verbose)
				print_one_reag(tmp_reag,
                               //the_list, 
                        Genomes, nb_spec, old_genome);
			
			carry_on_reag(tmp_reag,
                          //the_list, 
                    Genomes, nb_spec, spec_left, &did_merge, distmem);

			if (verbose) {
				special_print_genomes2(stdout, Genomes, nb_spec);
			}
			
			//if (the_list->merge_with != -1) /* there was a merge and we need to update tree */
			if (tmp_reag->merge_with != -1) /* there was a merge and we need to update tree */
				update_for_merge(Genomes, nbreag, distmem, treemem, 
                        //the_list->spec, the_list->merge_with,
                        tmp_reag->spec, tmp_reag->merge_with,
								 nb_spec, *spec_left, verbose,cb_max_chromo_size);
			
			nbreag[tmp_reag->spec]++; 
            //nbreag[the_list->spec]++;

			*last_genomenb = tmp_reag->spec;
			//*last_genomenb = the_list->spec; 
			found_something = TRUE ;
		}
	}    
	
	  
//erase_list_reag(the_list);
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
	if (verbose) 
	{
		fprintf(stdout,"solve with good reag, look for a reversal\n");		
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
		if (verbose) 
		{
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
		if (verbose) 
        {
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
		if (verbose)
        {
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
				fprintf(stdout,"look for a reversal 2\n");		
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
					fprintf(stdout,"look for a transloc2\n");		
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
					fprintf(stdout,"look for a fusion/fission2\n");		
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

void free_cap_mem (int nb_spec)
{
    int i;
    if( mgr_genome_list_cap != NULL)
    {
    	for (i=0; i<nb_spec; i++) {
			free(mgr_genome_list_cap[i].gnamePtr);
			free(mgr_genome_list_cap[i].genes);
            free(mgr_genome_list_cap[i].delimiters);
			if (mgr_genome_list_cap[i].alphabet != NULL) {
				if (mgr_genome_list_cap[i].num_chr>0) {
					free_alphabet(mgr_genome_list_cap[i].alphabet, 
                            mgr_genome_list_cap[i].num_g/3);
				}
				else {
					free_alphabet(mgr_genome_list_cap[i].alphabet, 
                            mgr_genome_list_cap[i].num_g);
				}		
			}
		}
    }
    free(mgr_genome_list_cap);
}

void free_mem_4_mgr()
{
    free_G_struct(&Genomes,4); 
    free_G_struct(&Genomes_copy,3);

   free_cap_mem(3);

    free_cbounds(cb_max_chromo_size);

    free(nbreag);

    mcdist_freemem(&distmem_mgrmed);
    mcdist_freemem(&distmem_mgrinvdist);
    mcdist_freemem(&distmem_capgraph);

    free_list_listreag_memory (&list_carry_best_reag);
    free_list_listreag_memory (&list_carry_best_reag2);
    free_list_listreag_memory (&list_do_k_reag);

}

void init_cbounds_memory (int num_genes, int num_chromosomes, cbounds_t * cb)
{
//    cb  = (cbounds_t *) malloc( sizeof(cbounds_t));
    if( (cb == NULL) || ( cb == ( cbounds_t *)NULL))
        failwith ("initialization for cbouns_t failed");
    cb->num_genes = num_genes;
    cb->num_chromosomes = num_chromosomes;
    cb->cNum = (int *) malloc(num_genes*sizeof(int));
    if( (cb->cNum == NULL) || (cb->cNum == (int *)NULL))
        failwith ("initialization for cbouns_t.cNum failed");
	cb->cBound   = (int *) malloc((num_chromosomes+1)*sizeof(int));
    if( (cb->cBound == NULL) || (cb->cBound == (int *)NULL)) 
        failwith ("initialization for cbouns_t.cBound failed");
}

void mgr_ini_mem (int num_genes)
{
    int i;
    int num_genomes = 3; //median of 3 genomes
    int num_chromosomes = num_genes; //consider the worst case, each loci is a chromosome
     
    mcdist_allocmem(3*num_genes, num_chromosomes, &distmem_mgrmed);
    mcdist_allocmem(3*num_genes, num_chromosomes, &distmem_mgrinvdist);
    mcdist_allocmem(3*num_genes, num_chromosomes, &distmem_capgraph);

    nbreag = (int *)malloc(3 * sizeof(int));
   
    cb_max_chromo_size = (cbounds_t *) malloc( sizeof(cbounds_t));
    init_cbounds_memory(3*num_genes, num_genes, cb_max_chromo_size); 
	
    mgr_genome_list_copy =
    (struct mgr_genome_struct *) malloc(num_genomes*sizeof(struct mgr_genome_struct));
    for (i=0 ; i<num_genomes; i++) {
    (mgr_genome_list_copy)[i].gnamePtr = (char *) malloc(MAX_STR_LEN*sizeof(char));
   //if we are dealing with multichromosome, in the worst case, each loci is a singlechromosome ($a$b$c$....), we will need 3 times of original num_genes
    (mgr_genome_list_copy)[i].genes = (int*)malloc(3*num_genes*sizeof(int));
    (mgr_genome_list_copy)[i].delimiters = (int*)malloc(num_genes*sizeof(int));
    (mgr_genome_list_copy)[i].num_delimiters = 0;
     // alphabet is useful only in condense function inside MGR, anyway, we still
   // initialize it here, just in case.
	initialize_alphabet(&(mgr_genome_list_copy[i].alphabet), num_genes, TRUE);
 //   fprintf(stdout, "check alphabet of mgr_genome_list copy, i=%d : \n",i); fflush(stdout);
   // print_alphabet((mgr_genome_list_copy[i].alphabet), num_genes);
    }
    mgr_genome_list = (struct mgr_genome_struct *) malloc(2*(num_genomes-1)*sizeof(struct mgr_genome_struct));
    for (i=0 ; i<2*(num_genomes-1); i++) {
    (mgr_genome_list)[i].gnamePtr = (char *) malloc(MAX_STR_LEN*sizeof(char));
    //if we are dealing with multichromosome, in the worst case, each loci is a singlechromosome ($a$b$c$....), we will need 3 times of original num_genes
    (mgr_genome_list)[i].genes = (int*)malloc(3*num_genes*sizeof(int));
    (mgr_genome_list)[i].delimiters = (int*)malloc(num_genes*sizeof(int));
    (mgr_genome_list)[i].num_delimiters = 0;
   // (mgr_genome_list)[i].alphabet = (a_strip *) NULL; 
   // alphabet is useful only in condense function inside MGR, anyway, we still
   // initialize it here, just in case.
	initialize_alphabet(&(mgr_genome_list[i].alphabet), num_genes, TRUE);
  //  fprintf(stdout, "check alphabet in mgr_genome_list, i=%d : \n",i); fflush(stdout);
  //  print_alphabet((mgr_genome_list[i].alphabet), num_genes);
    }

    int nb_spec = num_genomes;
    Genomes.label = (int *)malloc(2*(nb_spec-1) * sizeof(int));
    Genomes.same_as = (int *)malloc(2*(nb_spec-1) * sizeof(int));
    Genomes.nb_chromo = (int *)malloc(2*(nb_spec-1) * sizeof(int));
    Genomes.dist_mat = (int *)malloc(2*(nb_spec-1) * 2*(nb_spec-1) * sizeof(int));

    Genomes_copy.label = (int *)malloc(nb_spec * sizeof(int));
    Genomes_copy.same_as = (int *)malloc(nb_spec * sizeof(int));
    Genomes_copy.nb_chromo = (int *)malloc(nb_spec * sizeof(int));
    Genomes_copy.dist_mat = (int *)malloc(nb_spec * nb_spec * sizeof(int));

    mgr_genome_list_cap = 
     (struct mgr_genome_struct *) malloc(3*sizeof(struct mgr_genome_struct));
    for (i=0 ; i<3; i++) {
    (mgr_genome_list_cap)[i].gnamePtr = (char *) malloc(MAX_STR_LEN*sizeof(char));
    (mgr_genome_list_cap)[i].genes = (int*)malloc(3*num_genes*sizeof(int));
    (mgr_genome_list_cap)[i].delimiters = (int*)malloc(num_genes*sizeof(int));
    (mgr_genome_list_cap)[i].num_delimiters = 0;
    (mgr_genome_list_cap)[i].alphabet = (a_strip *) NULL; 
    }


    

}

void clear_up_genome (struct mgr_genome_struct * in_genome)
{

    int i=0; int size;
    int * gene = in_genome->genes;
    size = in_genome->num_g; 
    for(i=0;i<size;i++)
    {
        gene[i]=0;
    }
}
//input : genes array, delimiters array
//create a single-chromosome genome that mimic the multi-chromosome genome,
//store it in mgr_genome_list_copy.
void multi_to_single ( int * out_genes, int * in_genes, int * in_delimiters, int size_gene, int size_delimiters, int max_size_deli)
{
  int cap = size_gene+1;
  int i,j;
  int outpos = 0, inpos = 0;
//  fprintf(stdout, "start ..... size_deli=%d \n",size_delimiters);
  int check_sum=0;
  for (i=0;i<size_delimiters;i++)
      check_sum+= in_delimiters[i];
  if(check_sum != size_gene)
  {
      fprintf(stderr,"check_sum = %d != size_gene = %d\n",
              check_sum, size_gene);
      fflush(stdout);
      failwith("ERROR: in multi_to_single, mgr.c, nonsence delimiters array ");
  }
  for (i=0;i<size_delimiters;i++)
  {
      int count = in_delimiters[i];
      out_genes[outpos] = cap;
      cap++; outpos++;
      for(j=inpos;j<inpos+count;outpos++,j++)
      {
          out_genes[outpos] = in_genes[j];
      }
      inpos = inpos + count;
      out_genes[outpos] = cap; 
      cap ++; outpos ++;
  } 
  if(max_size_deli>size_delimiters)
  {
      int diff = max_size_deli - size_delimiters;
      for(i=outpos;i<diff*2+outpos;i++)
      {
            out_genes[i] = cap; 
            cap++;
      }
  }
}

int single_to_multi ( int * in_gene, int * out_gene, int * out_delimiters, int gene_size, int delimiters_size)
{
    assert(in_gene != (int*)NULL); assert( out_gene != (int*)NULL); assert(out_delimiters!= (int*)NULL);
    int i=0,j=0,k=0,count=0;
    int alpha_size = gene_size - delimiters_size*2;
    for (i=0;i<gene_size;i++)
    {
        if( (in_gene[i]<=alpha_size)&&(in_gene[i]>=(-alpha_size)) )
        {
            out_gene[j]=in_gene[i];  
            j++; count++; 
        }
        else
        {
            if(count!=0)
            {
                out_delimiters[k] = count;
                k++; count = 0;
            }
        }
    }
    return k;

}




//distance with mgr distance function.
int mgr_invdist (int * g1, int * g2, int num_genes, int * deli1, int * deli2, int num_deli1, int num_deli2)
{
   // mgr_distmem_t dist_mem;//move this out later
    mgr_distmem_t * dist_mem = &distmem_mgrinvdist;
    struct mgr_genome_struct * mgr_genome_pair;
    int max_num_deli = 0;
    int res = 0;
    if(num_deli1>num_deli2) max_num_deli = num_deli1;
    else max_num_deli = num_deli2;
    /* debug msg 
    int x=0;
    fprintf(stdout,"mgr_invdist, g1 = {");
    for(x=0;x<num_genes;x++)
        fprintf(stdout,"%d,",g1[x]);
    fprintf(stdout,"} g2 = {");
    for(x=0;x<num_genes;x++)
    {
        if (g2[x]==0) failwith("g2=0");
        fprintf(stdout,"%d,",g2[x]);
    }
    fprintf(stdout,"}\n deli1 = [");
    for(x=0;x<num_deli1;x++)
        fprintf(stdout,"%d,",deli1[x]);
    fprintf(stdout,"] deli2 = [");
    for(x=0;x<num_deli2;x++)
        fprintf(stdout,"%d,",deli2[x]);
    fprintf(stdout,"] \n"); fflush(stdout);
     debug msg */
    mgr_genome_pair = 
     (struct mgr_genome_struct *) malloc(2*sizeof(struct mgr_genome_struct));
    int i;
    for (i=0 ; i<2; i++) {
    (mgr_genome_pair)[i].gnamePtr = (char *) malloc(MAX_STR_LEN*sizeof(char));
    (mgr_genome_pair)[i].genes = (int*)malloc(3*num_genes*sizeof(int));
    (mgr_genome_pair)[i].delimiters = (int*)malloc(num_genes*sizeof(int));
    (mgr_genome_pair)[i].num_delimiters = 0;
    (mgr_genome_pair)[i].alphabet = (a_strip *) NULL; 
    }
    multi_to_single (mgr_genome_pair[0].genes, g1, deli1, num_genes, num_deli1, max_num_deli);
    multi_to_single (mgr_genome_pair[1].genes, g2, deli2, num_genes, num_deli2, max_num_deli);
    num_genes = num_genes + max_num_deli * 2 ;
   // not sure, do we need this? 
   // add_null_chromos(mgr_genome_list_copy, nb_spec, &num_genes, &num_chromosomes);
  // fprintf(stdout, "alloc mem for mcdist, %d,%d\n",num_genes,max_num_deli);
  //  mcdist_allocmem(num_genes, max_num_deli, &dist_mem);
    res = mcdist_noncircular
        (&(mgr_genome_pair[0]), &(mgr_genome_pair[1]), num_genes, max_num_deli, dist_mem, NULL);

    
  //  mcdist_freemem(&dist_mem);
    free_genome_list(mgr_genome_pair,2);
    /*debug msg
     fprintf(stdout," res = %d\n",res); fflush(stdout);
     debug msg*/
    return res;
}

/* given a pair of delimeters: in_deli1 and in_deli2, fill out_g2 with a better "capping" of g2 */
void better_capping (int * g1, int * g2, int num_genes, int * in_deli1, int * in_deli2, int num_deli1, int num_deli2, struct genome_struct * out_g2)
{
    int i;
    struct mgr_genome_struct *tmplist = NULL;
    mgr_distmem_t * dist_mem_cap = &distmem_capgraph;
    int max_num_deli = 0;
    if(num_deli1>num_deli2) max_num_deli = num_deli1;
    else max_num_deli = num_deli2;
/*
    struct mgr_genome_struct * mgr_genome_list;
    mgr_genome_list = 
     (struct mgr_genome_struct *) malloc(3*sizeof(struct mgr_genome_struct));
    int i;
    for (i=0 ; i<3; i++) {
    (mgr_genome_list)[i].gnamePtr = (char *) malloc(MAX_STR_LEN*sizeof(char));
    (mgr_genome_list)[i].genes = (int*)malloc(3*num_genes*sizeof(int));
    (mgr_genome_list)[i].delimiters = (int*)malloc(num_genes*sizeof(int));
    (mgr_genome_list)[i].num_delimiters = 0;
    (mgr_genome_list)[i].alphabet = (a_strip *) NULL; 
    }
*/
    multi_to_single (mgr_genome_list_cap[0].genes, g1, in_deli1, num_genes, num_deli1, max_num_deli);
    multi_to_single (mgr_genome_list_cap[1].genes, g2, in_deli2, num_genes, num_deli2, max_num_deli);
    num_genes = num_genes + max_num_deli * 2 ;

    int num_chromosomes = max_num_deli;
    tmplist = mgr_genome_list_cap;
    int dis = mcdist_capgraph 
        (tmplist, tmplist++, num_genes, num_chromosomes, dist_mem_cap,NULL);
    for(i=0;i<num_genes;i++)
    {
        mgr_genome_list_cap[2].genes[i] = dist_mem_cap->cappedp2[i];
    }
    int real_deli_num = single_to_multi ( mgr_genome_list_cap[2].genes, out_g2->genes, out_g2->delimiters, num_genes,max_num_deli);
    out_g2->deli_num =  real_deli_num;

   //  free_genome_list(mgr_genome_list,3);

}

void mgr_med (int * g1, int * g2, int * g3, int * deli1, int * deli2, int * deli3, int num_deli1, int num_deli2, int num_deli3, int SIZE_ALPHA, int CIRCULAR, struct genome_struct * g_med)
{  
    mgr_distmem_t * dist_mem = &distmem_mgrmed;
    int i,j;
    int NumGenomes = 3;
    int circular;// = CIRCULAR;
    int nb_spec=3;
    int num_genes;
    int size_alpha = SIZE_ALPHA;
    int num_chromosomes;
    int genome_type=GLINEAR;
    if (CIRCULAR) genome_type = GCIRCULAR; 
    if((num_deli1>=1)||(num_deli2>=1)||(num_deli3>=1))
    {
        genome_type = GMULTI;
    }
    int condensing = FALSE;
    int optimize = FALSE;
	int alternatives = FALSE;
    int verbose = FALSE;
    int depth=2;/* depth of the search when running out of good rearrangements */
    int heuristic = 0; // -H: heuristic to speed up triplet resolution
    /* Data structure for the phylogeny */
    treemem_t treemem;	/* the tree itself */
    int spec_left;		/* number of genomes left to be put in the tree */
  //  int *nbreag = NULL;/* number of rearrangements that were carried in each genome */
    int total_weight;	/* total weight of the tree */
    double avg_nb_rev=0;
    int max_num_deli = 0;
    if (genome_type == GMULTI) { // multichromosomal distance
    // we've assume that we have the maximum number of chromosomes
    // in order not to reassign memory
       if(num_deli1>num_deli2) max_num_deli = num_deli1;
       else max_num_deli = num_deli2;
       if(max_num_deli<num_deli3) max_num_deli = num_deli3;
       num_genes = size_alpha+2*max_num_deli;
        num_chromosomes = max_num_deli; 
    }
    else  { // unichromosomal distance
        num_genes = size_alpha;  num_chromosomes = 0;
    }
    if (genome_type == GMULTI)
    {
        multi_to_single 
        ( (mgr_genome_list_copy)[0].genes, g1, deli1, size_alpha, num_deli1, max_num_deli);
        multi_to_single 
        ( (mgr_genome_list_copy)[1].genes, g2, deli2, size_alpha, num_deli2, max_num_deli);
        multi_to_single 
        ( (mgr_genome_list_copy)[2].genes, g3, deli3, size_alpha, num_deli3, max_num_deli);
    }
    else {
       for (j=0;j<size_alpha;j++)
       {
        (mgr_genome_list_copy)[0].genes[j] = g1[j]; 
        (mgr_genome_list_copy)[1].genes[j] = g2[j]; 
        (mgr_genome_list_copy)[2].genes[j] = g3[j];
       }
    }
    // add_null_chromos adds null chromo at the end of each genome
    // num_genes will be 3*(original gene size).
    if (num_chromosomes > 0)
        add_null_chromos(mgr_genome_list_copy, nb_spec, &num_genes, &num_chromosomes);
    /* debug msg
    fprintf(stdout,"check mgr gene: \n");
    for(i=0;i<3;i++)
    {
        for(j=0;j<num_genes;j++)
            fprintf(stdout,"%d,",(mgr_genome_list_copy)[i].genes[j]);
        fprintf(stdout,"\n");
    }
    fprintf(stdout,"end of checking\n"); fflush(stdout);
    debug msg*/
    // finish initialization of the input genomes, including actual nb_chromo
  //  num_chromosomes = max_num_deli;
  if(verbose)
  { 
      fprintf(stdout, "\n mgr.c mgr_med , %d,%d,%d\n",num_deli1,num_deli2,num_deli3); fflush(stdout);
      fprintf(stdout,"init G struct, %d,%d,%d\n",nb_spec,num_genes,num_chromosomes); 
      fflush(stdout);
      fprintf(stdout, "check cb: %d\n",cb_max_chromo_size->num_genes); }
    init_G_struct
        (&Genomes_copy, mgr_genome_list_copy, nb_spec, num_genes, num_chromosomes,
         cb_max_chromo_size);
    // assign memory for genomes assuming each marker is his own chromosome (worst case)
    initialize_genome_list(&mgr_genome_list, 2*(nb_spec-1), num_genes, num_chromosomes);
    // initialize structure for Genomes which includes true nb_chromo for each genome
    init_G_struct(&Genomes, mgr_genome_list, 2*(nb_spec-1), num_genes, num_chromosomes,
            cb_max_chromo_size);
    // we also keep the genome_list_copy as a reference
    copy_first_n(&Genomes_copy, &Genomes, nb_spec);
    for (i=0; i<nb_spec; i++)
        strcpy((mgr_genome_list_copy)[i].gnamePtr, (mgr_genome_list)[i].gnamePtr);

    // corresponds to the true size of the alphabet
    size_alpha = num_genes - 2*num_chromosomes;
 //   if (condensing == TRUE) {
//        condense_genomes(&Genomes, nb_spec, verbose);
		//uncondense_genomes(&Genomes, nb_spec);
		//special_print_genomes2(stdout, &Genomes, nb_spec);
//	}

  //  mcdist_allocmem(num_genes, num_chromosomes, &dist_mem);
    compute_dist_mat(&Genomes, nb_spec, dist_mem, cb_max_chromo_size);
    if (verbose){
			fprintf(stdout, "\nInitial pairwise distance matrix:\n"); fflush(stdout);
			print_dist_mat(&Genomes, nb_spec);
	}	

   // struct mgr_genome_struct *tmplist = NULL;//this is for debug msg
        
     // allocate memory for the tree/phylogeny that we want to build
    initialize_treemem(&treemem, nb_spec);

    // the number of species left to place on the tree
	spec_left=nb_spec;
			
    // the number of transformations we did on each genome
  //  nbreag = (int *)malloc(nb_spec * sizeof(int));
    for (i=0;i<nb_spec;i++)
          nbreag[i]=0;

    // check to see if two genomes are identical before we begin
	check_for_merge(&Genomes, nbreag, nb_spec, &spec_left,
					dist_mem, &treemem, verbose, cb_max_chromo_size);
    
    
    if (spec_left == 3) {
        int perfect_triplet;
        if (verbose) {
        solve_triplet(&Genomes, nbreag, &perfect_triplet, 2,
                          nb_spec, &spec_left,
                          depth, dist_mem, &treemem, &Preancestors, 
                          heuristic, condensing, verbose);
        } else {
            solve_triplet(&Genomes, nbreag, &perfect_triplet, 2,
                          nb_spec, &spec_left,
                          depth, dist_mem, &treemem, NULL, 
                          heuristic, condensing, verbose);
        }
    }
/* debug msg 
    int ii,jj;
    //struct mgr_genome_struct *
     //   tmplist= Genomes_copy.genome_list;
    tmplist = Genomes.genome_list;
    for(ii=0;ii<4;ii++)
    {
                fprintf(stdout,"ii=%d : [",ii);
                int * point = tmplist->genes;
                for(jj=0;jj<num_genes;jj++)
                    fprintf(stdout,"%d,",point[jj]);
                fprintf(stdout,"]\n");
                fflush(stdout);
                tmplist ++;
    }
 debug msg */
    /* do we need this ?
  // only need to do this if we haven't been able to solve the problem
    // using only good reversals
    if (find_total_dist(&Genomes, nb_spec, &dist_mem) !=0) {
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
                   depth, &dist_mem, &treemem, heuristic, condensing, verbose) ;
        // at each iteration, adds a new genome to the tree
        while (treemem.tree_size < nb_spec && did_grow==1) {
            if (verbose) {
                fprintf(stdout, "Building phylogeny, %d %s remaining to place\n",
                        nb_spec-treemem.tree_size, (nb_spec-treemem.tree_size>1)?"genomes":"genome");
                fflush(stdout);
            }
            growing_tree(&Genomes, &Star, nbreag, nb_spec, &spec_left,
                         depth, &dist_mem, &treemem, &did_grow, heuristic, condensing, verbose) ;
            if (verbose) {
                special_print_genomes2(stdout, &Genomes, 2*(nb_spec-1));
            }
        }
        free_G_struct(&Star, 3);
        // these good rearrangements weren't counted in the edges yet
        add_nbrev_edges(treemem.the_edge_list, nbreag, nb_spec);
    }// if (find_total_dist(&Genomes, nb_spec, &dist_mem) !=0)
    else { // when the tree was solve only using good rearrangements
        // need to finish the list with the last 3 genomes
        fprintf(stdout,"finish edge list\n");
        fflush(stdout);
        finish_edge_list(&Genomes, nbreag, nb_spec, &treemem);
    }// if (find_total_dist(&Genomes, nb_spec, &dist_mem) !=0)
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
    
  
    if (genome_type == GMULTI)
    {
 /*       fprintf(stdout,"\n TESTING, num_genes=%d,num_chrom=%d\n",
                                    num_genes,num_chromosomes);
       tmplist = Genomes_copy.genome_list;
        int test;
        for(ii=0;ii<4;ii++) {
        test = mcdist_capgraph 
            (Genomes.genome_list, tmplist+i, num_genes, num_chromosomes, dist_mem_cap,NULL);
        fprintf(stdout,"mcdist capgraph = %d\n",test);
        }
        fflush(stdout);
*/
        remove_null_chromos (Genomes.genome_list, nb_spec, &num_genes, &num_chromosomes, size_alpha, max_num_deli);
        int real_deli_num = single_to_multi ( (Genomes.genome_list)->genes, g_med->genes, g_med->delimiters, num_genes,max_num_deli);
        g_med->deli_num =  real_deli_num;//max_num_deli;
    }
    else
    {
        copy_genes((Genomes.genome_list)->genes,g_med->genes,num_genes);
        g_med -> deli_num = 0 ;
    }

    // free memory for the number of rearrangements
    //free(nbreag);
    // free memory for phylogeny
    free_treemem(&treemem);
   
   // mcdist_freemem(&dist_mem);
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


