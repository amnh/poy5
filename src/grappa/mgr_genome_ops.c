/*  genome_ops.c
*
*  SUMMARY:
*
*		All the basic operations associated with genome allocation and
*		manipulation (including condensing, uncondensing).
*
* Copyright (C) 2004-2006 Genome Insitute of Singapore
* Copyright (C) 2000-2002 University of Southern California
*
* See file COPYRIGHT for details. 
*****************************************************************************
* This file is part of MGR.
*
* MGR is free software; you can redistribute it and/or modify
* it under the terms of the GNU General Public License, Version 2,
* dated June 1991, as published by the Free Software Foundation.
*
* MGR is distributed in the hope that it will be useful, but
* WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
* General Public License for more details.
*
* You should have received a copy of the GNU General Public License along
* with this program; if not, write to the Free Software Foundation, Inc.,
* 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
*/

/* Last modified on Wed Aug 2, 2006, by Guillaume Bourque
*/

#include "mgrstructs.h"
#include "mgr_e_malloc.h"
#include "mgr_mcread_input.h"
#include "mgr_write_data.h"
#include "mgr_my_testrev.h"
#include "mgr_list_ops.h"
#include "mgr_genome_print.h"
#include "mgr_genome_ops.h"

void initialize_alphabet(a_strip **alphabet, int size_alpha, int identity) {
	int j;
	(*alphabet) = (a_strip *) e_malloc(size_alpha*sizeof(a_strip), "alphabet");
	for (j=0; j<size_alpha; j++) {
		(*alphabet)[j].size = MAXSTRIP;
		(*alphabet)[j].strip = (int *)e_malloc(MAXSTRIP*sizeof(int), "alphabet");
		if (identity == TRUE) {
			(*alphabet)[j].last = 0;
			(*alphabet)[j].strip[0] = j+1;  // each strip is actually the gene itself
		}
		else {
			(*alphabet)[j].last = -1;
		}
	}
}

void free_alphabet(a_strip *alphabet, int size_alpha) {
	int j;
    if (alphabet != NULL)
    {
        //fprintf(stdout, "in free alphabet %d\n", size_alpha);
        for (j=0; j<size_alpha; j++) {
            //fprintf(stdout, "freeing strip %d of size %d\n", j, alphabet[j].size);
           if(alphabet[j].strip != (int*) NULL)  free(alphabet[j].strip);
        }
        free(alphabet);
    }
}

void check_alpha(a_strip *alphabet, int size_alpha, int ori_size_alpha) {
	int i, j;
	int *array;
	
	array = (int *)e_malloc(ori_size_alpha*sizeof(int), "array");
	for (i=0;i<ori_size_alpha;i++) {
		array[i]=0;
	}
	for (i=0;i<size_alpha;i++) {
		for (j=0;j<=alphabet[i].last;j++) {
			if (abs(alphabet[i].strip[j])>=1 && abs(alphabet[i].strip[j])<=ori_size_alpha) {
				array[abs(alphabet[i].strip[j])-1]++;
			}
			else {
				fprintf(stdout, "alphabet[%d].strip[%d]=%d is not between 1 and %d\n", i, j, alphabet[i].strip[j], ori_size_alpha);
				exit(-1);
			}
		}
	}
	for (i=0;i<ori_size_alpha;i++) {
		if (array[i]!=1) {
			fprintf(stdout, "There are %d copies of gene %d in alphabet\n", array[i], i+1);
			exit(-1);
		}
	}
	

	free(array);
}
	

void copy_alphabet(a_strip *alphabet1, a_strip **alphabet2, int size_alpha1, int size_alpha2) {
	int i, j;
	
//	fprintf(stdout, "in copy_alphabet (%d %d)\n", size_alpha1, size_alpha2);
//	print_alphabet(alphabet1, size_alpha1);
//	print_alphabet(*alphabet2, size_alpha2);
	
	if (size_alpha1 != size_alpha2) {
		free_alphabet(*alphabet2, size_alpha2);
		initialize_alphabet(alphabet2, size_alpha1, TRUE);
	}
			
	for (i=0; i<size_alpha1; i++) {
		if ((*alphabet2)[i].size != alphabet1[i].size) {
		//if (alphabet2[i].size-1 < alphabet1[i].last) {
			free((*alphabet2)[i].strip);
			(*alphabet2)[i].strip = (int *)e_malloc((alphabet1[i].size)*sizeof(int), "alphabet");
			(*alphabet2)[i].size = alphabet1[i].size;
		}
		for (j=0; j<alphabet1[i].size; j++)  {
			(*alphabet2)[i].strip[j] = alphabet1[i].strip[j];
		}
		(*alphabet2)[i].last = alphabet1[i].last;
	}
	
	//print_alphabet(*alphabet2, size_alpha1);

}
	
void copy_genome_struct(struct mgr_genome_struct *genome1, struct mgr_genome_struct *genome2, int need_copy_alphabet) {
	int size_alpha1, size_alpha2;
	
	if (genome1->num_chr>0) {
		size_alpha1 = genome1->num_g/3;
		size_alpha2 = genome2->num_g/3;
	}
	else {
		size_alpha1 = genome1->num_g;
		size_alpha2 = genome2->num_g;
	}
	
	//fprintf(stdout, "in copy struct (%d %d)\n", genome1->num_g, genome2->num_g);
	//fprintf(stdout, "genome2 before\n");
	//special_print_genome2(stdout, genome2, 1);
		
	copy_genes(genome1->genes, genome2->genes, genome1->num_g);

	
	if ((need_copy_alphabet == TRUE)&&(genome1->alphabet != NULL)) {
		copy_alphabet(genome1->alphabet, &(genome2->alphabet), size_alpha1, size_alpha2);
		genome2->num_g = genome1->num_g;
		genome2->num_chr = genome1->num_chr;
		genome2->ori_num_g = genome1->ori_num_g;
		
		//fprintf(stdout, "we've copied the alphabet\n");
	}
	else {
		// we need to free old alphabet if there is one and create a new empty alphabet

		
		if (genome2->alphabet != NULL) {
			//fprintf(stdout, "old alphabet\n");
			//print_alphabet(genome2->alphabet, genome2->num_g);
			free_alphabet(genome2->alphabet, size_alpha2);
		}
		
		initialize_alphabet(&(genome2->alphabet), size_alpha1, TRUE);
		genome2->num_g = genome1->num_g;
		genome2->num_chr = genome1->num_chr;
		genome2->ori_num_g = genome2->num_g; // num_g == ori_num_g
				
		//fprintf(stdout, "we've emptied the alphabet\n");

	}
	//print_alphabet(genome2->alphabet, genome2->num_g);

	//fprintf(stdout, "genome2 after %d\n", genome2->num_g);
	//special_print_genome2(stdout, genome2, 1);

}



void copy_genomes(G_struct *Genomes1, int gindex1,
				  G_struct *Genomes2, int gindex2,
				  int copy_alphabet) 
{
	struct mgr_genome_struct *genome_list1 = Genomes1->genome_list;
	struct mgr_genome_struct *genome_list2 = Genomes2->genome_list;
	
	//fprintf(stdout, "in copy_genomes\n");
	
	copy_genome_struct(&genome_list1[gindex1], &genome_list2[gindex2], copy_alphabet);
	//print_alphabet(genome_list2[gindex2].alphabet, genome_list2[gindex2].num_g);

	
	Genomes2->nb_chromo[gindex2] = Genomes1->nb_chromo[gindex1];
	Genomes2->same_as[gindex2] = Genomes1->same_as[gindex1];
	
	//fprintf(stdout, "two copies?\n");
	//special_print_genome2(stdout, &genome_list1[gindex1], Genomes1->nb_chromo[gindex1]);
	//special_print_genome2(stdout, &genome_list2[gindex2], Genomes2->nb_chromo[gindex2]);
	//print_alphabet(genome_list2[gindex2].alphabet, genome_list2[gindex2].num_g);
	
	//fprintf(stdout, "out copy_genomes\n");

	
}

void copy_first_n(G_struct *Genomes1, G_struct *Genomes2, int nb_spec)
{
    int i;

    for (i=0; i<nb_spec; i++)
        copy_genomes(Genomes1, i, Genomes2, i, TRUE);
	
	Genomes2->max_chromo_size = Genomes1->max_chromo_size;
		
}



// In this procedure gindex1 is now like gindex2 so we modify treemem accordingly
void update_for_merge(G_struct *Genomes, int *nbreag, mgr_distmem_t *distmem,
					  treemem_t *treemem, int gindex1, int gindex2, int nb_spec, 
					  int spec_left, int verbose)
{
	struct mgr_genome_struct *genome_list = Genomes->genome_list;

	if (verbose)
		fprintf(stdout,"genome %d is now like genome %d\n", gindex1, gindex2);
	
	if (spec_left >= 3) {
		treemem->merge_edge_list = list_edge_insert(treemem->merge_edge_list, 
													gindex1, Genomes->label[gindex1], 
													treemem->new_node,treemem->new_node,
													nbreag[gindex1]);
		
		treemem->merge_edge_list = list_edge_insert(treemem-> merge_edge_list, 
													gindex2, Genomes->label[gindex2], 
													treemem->new_node,treemem->new_node,
													nbreag[gindex2]);
	
		nbreag[gindex1]=0;
		nbreag[gindex2]=0;
		Genomes->label[gindex2]=treemem->new_node;
		copy_genome_struct(&genome_list[gindex2], &genome_list[treemem->new_node], TRUE);
		Genomes->nb_chromo[treemem->new_node] = Genomes->nb_chromo[gindex2];
		treemem->new_node++;
	} else {
		// there are only 3 genomes left so if two of them are merged we are done
		int gindex3;
		
		// find the third genome
		gindex3 = 0;
		while ((gindex3 == gindex1 || gindex3 == gindex2 || Genomes->same_as[gindex3]!=-1) && gindex3<nb_spec) {
			gindex3++;
		}
		
		//update nb of rearrangements
		nbreag[gindex3] += Genomes->dist_mat[gindex2*nb_spec + gindex3];
			
		//copy gene order
		copy_genome_struct(&genome_list[gindex2], &genome_list[gindex3], TRUE);
		Genomes->nb_chromo[gindex3] = Genomes->nb_chromo[gindex2];
			
fprintf(stdout,"update for merge\n"); fflush(stdout);
		// should now all be zero
		compute_dist_mat(Genomes, nb_spec, distmem);
		
		// reset to pretend we have 3 final genomes (all identical)
		Genomes->same_as[gindex1] = -1;
		spec_left = 3;
	
	}
}


void check_for_merge(G_struct *Genomes, int *nbreag, 
						 int nb_spec, int *spec_left, 
						 mgr_distmem_t *distmem, treemem_t *treemem, int verbose) {
  //	struct mgr_genome_struct *genome_list = Genomes->genome_list;
	int i, j;
	
	for (i=0; i<nb_spec; i++) {
		if (Genomes->same_as[i] == -1) {
			for (j=i+1; j<nb_spec; j++) {
				if (Genomes->same_as[j] == -1 &&
					Genomes->dist_mat[i*nb_spec + j] == 0) {// Genomes i and j are identical
					
					if (*spec_left>=3) {
//					if (*spec_left>3) {
						Genomes->same_as[i] = j;
						(*spec_left)--;
					}
					update_for_merge(Genomes, nbreag, distmem, treemem, i, j, nb_spec, *spec_left, verbose);
					j=nb_spec;	//to leave the loop
					
				}
			}
		}
	}
	fprintf(stdout,"check_for_merge\n"); fflush(stdout);				
	compute_dist_mat(Genomes, nb_spec, distmem);
					
//	}
}

// Note: initialize_genome_list here no longer alloc memory for mgr_genome_list, instead,
//all mgr_genome_list should be defined in mgr.h and initlized in mgr_ini_mem of mgr.c
// just start with identity genome for all species 
void initialize_genome_list(struct mgr_genome_struct **genome_list, int nb_spec,
                            int num_genes, int num_chromosomes) {
	int i,j;
	int size_alpha;
	int *new_genes;
	
	/* allocate space for the genomes */
//	*genome_list = (struct mgr_genome_struct *) malloc(nb_spec*sizeof(struct mgr_genome_struct));
	
	for (i=0; i<nb_spec; i++) {

		/* allocate space to store the name */
//		((*genome_list)[i]).gnamePtr = (char *) malloc(MAX_STR_LEN*sizeof(char));

		if (i<nb_spec/2 + 1)
			 sprintf(((*genome_list)[i]).gnamePtr, "Genome%d", i+1);
		else {
		    //if (nb_spec!=4)
			 sprintf(((*genome_list)[i]).gnamePtr, "A%d", i+1);
		    // else
			// sprintf(((*genome_list)[i]).gnamePtr, "A");
		}


		/* allocate space to store the genes */
		//new_genes = (int *) malloc(num_genes*sizeof(int));
		new_genes = (*genome_list)[i].genes;
		
		if (num_chromosomes>0) {
			size_alpha = num_genes/3;
			/* initialize as a single identity chromosome with remaining null chromo at the end */
			
			new_genes[0] = size_alpha+1;		/* this is the first cap */
			for (j=1; j<=size_alpha; j++) {
				new_genes[j] = j;
			}
			for (j=size_alpha+1; j<num_genes; j++)
				new_genes[j] = j+1;
		}
		else {
			for (j=0; j<num_genes; j++)
				new_genes[j] = j+1;
		}	
	//	((*genome_list)[i]).genes = new_genes;		
		((*genome_list)[i]).num_g = num_genes;
		((*genome_list)[i]).ori_num_g = num_genes;
		((*genome_list)[i]).num_chr = num_chromosomes;
   //     fprintf(stdout,"check alphabet in initialize_genome_list :\n");
    //    print_alphabet(((*genome_list)[i]).alphabet, size_alpha);
	//	((*genome_list)[i]).alphabet = NULL;
	}
	

}


void free_genome_list(struct mgr_genome_struct *genome_list, int nb_spec) 
{
	int i;


	if (genome_list != NULL) {
		for (i=0; i<nb_spec; i++) {
			free(genome_list[i].gnamePtr);
			free(genome_list[i].genes);
            free(genome_list[i].delimiters);
			if (genome_list[i].alphabet != NULL) {
				if (genome_list[i].num_chr>0) {
					free_alphabet(genome_list[i].alphabet, genome_list[i].num_g/3);
				}
				else {
					free_alphabet(genome_list[i].alphabet, genome_list[i].num_g);
				}		
				//free((genome_list[i].alphabet)->strip);
				//free(genome_list[i].alphabet);
			}
		}
	}
	free(genome_list);
}




void print_G_struct(G_struct *Genomes, int nb_spec) {
	
	struct mgr_genome_struct *genome_list=Genomes->genome_list;
	int i;

	for (i=0; i<nb_spec; i++) {
		//fprintf(stdout, "num_g ori_num_g num_chr: %d %d %d\n", genome_list[i].num_g, genome_list[i].ori_num_g, genome_list[i].num_chr);
		//fprintf(stdout, "Same_as:%d Nb_chromo:%d\n", Genomes->same_as[i], Genomes->nb_chromo[i]);
		print_genome_multichrom_wcaps(&genome_list[i], genome_list[i].num_g, genome_list[i].num_chr);
		if (genome_list[i].num_chr>0) {
			//print_alphabet(genome_list[i].alphabet, genome_list[i].num_g/3);
		}
		else {
			//print_alphabet(genome_list[i].alphabet, genome_list[i].num_g);
		}		
	}
}

//init_G_struct here no longer alloc memory for G_struct, instead, all G_structs
//should be defined in mgr.h and initlized in mgr_ini_mem of mgr.c
void init_G_struct(G_struct *Genomes, struct mgr_genome_struct *genome_list,
                   int nb_spec, int num_genes, int num_chromosomes)
{
    int i, j;
	int size_alpha;
/*
    Genomes->label = (int *)malloc(nb_spec * sizeof(int));
    Genomes->same_as = (int *)malloc(nb_spec * sizeof(int));
    Genomes->nb_chromo = (int *)malloc(nb_spec * sizeof(int));
    Genomes->dist_mat = (int *)malloc(nb_spec * nb_spec * sizeof(int));
*/
    Genomes->genome_list = genome_list;
	if (num_chromosomes>0) {
		size_alpha = num_genes/3;
	}
	else {
		size_alpha = num_genes;
	}
    for (i=0; i<nb_spec; i++) {
        if (num_chromosomes == 0) {
            /* unichromosomal case: count it as 1 chromosome */
            Genomes->nb_chromo[i] = 1;
        }
        else {
            /* find last nonnull chromosome in i */
            Genomes->nb_chromo[i] = num_chromosomes;
            j = num_genes;
            while (Genomes->nb_chromo[i] > 0 && (&Genomes->genome_list[i])->genes[j-2] == j-1) {
                /* ends in two caps, which is our standardized null */
                j -= 2;
                Genomes->nb_chromo[i]--;
            }
        }
        Genomes->label[i] = i;
        Genomes->same_as[i] = -1;
		/* alphabet won't be used unless condensing option is selected */
		genome_list[i].num_g = num_genes;
		genome_list[i].ori_num_g = num_genes;
		genome_list[i].num_chr = num_chromosomes;
	//	move initialize_alphabet to mgr_ini_mem
    //	initialize_alphabet(&(genome_list[i].alphabet), size_alpha, TRUE);
		
        /*genome_list[i].alphabet = (a_strip *) e_malloc(num_genes*sizeof(a_strip), "alphabet");
		
		for (j=0; j<num_genes; j++) {
			(genome_list[i].alphabet[j]).size = MAXSTRIP;
			(genome_list[i].alphabet[j]).strip = (int *)e_malloc(MAXSTRIP*sizeof(int), "alphabet");
			
			(genome_list[i].alphabet[j]).last = 0;
			(genome_list[i].alphabet[j]).strip[0] = j+1;  // each strip is actually the gene itself
		}*/
		
    }
	find_max_chromo_size(Genomes, nb_spec);    
}

void free_G_struct(G_struct *Genomes, int nb_spec)
{

	free(Genomes->label);
	free(Genomes->same_as);
	free(Genomes->nb_chromo);
	free(Genomes->dist_mat);

	free_genome_list(Genomes->genome_list, nb_spec);
}

void print_dist_mat(G_struct *Genomes, int nb_spec) 
{
	int i, j;

	for (i=0; i<nb_spec; i++) {
		for (j=0; j<nb_spec; j++)
			fprintf(stdout, "%d ", Genomes->dist_mat[i*nb_spec + j]);
		fprintf(stdout, "\n");
	}	
}
void compute_dist_mat(G_struct *Genomes, int nb_spec, mgr_distmem_t *distmem) 
{
	struct mgr_genome_struct *genome_list = Genomes->genome_list;
	int *dist_mat = Genomes->dist_mat;
	
	int i, j;
	 	 
	
	for (i=0; i<nb_spec; i++) {
		dist_mat[i*nb_spec + i] = 0;
		for (j=i+1; j<nb_spec; j++) {
						
			if (Genomes->same_as[i]==-1 && Genomes->same_as[j]==-1) {
				dist_mat[i*nb_spec + j] = mcdist_noncircular(&genome_list[i], &genome_list[j], 
															 genome_list[i].num_g, genome_list[i].num_chr,
															 distmem, NULL);
				dist_mat[j*nb_spec + i] = dist_mat[i*nb_spec + j];
			}
			else dist_mat[i*nb_spec + j] = dist_mat[j*nb_spec + i] = 0;
		}
	}
fprintf(stdout,"compute_dist_mat, call find_max_chromo_size\n "); fflush(stdout);	
	// also compute max_chromo_size
	find_max_chromo_size(Genomes, nb_spec);

					 	 
}

void local_distmat_update(int gindex, G_struct *Genomes, int nb_spec, 
						  mgr_distmem_t *distmem)
{
	struct mgr_genome_struct *genome_list = Genomes->genome_list;
	int *dist_mat = Genomes->dist_mat;
	int i;

	for (i=0; i<nb_spec; i++) {
		if ((Genomes->same_as[gindex] == -1) &&
			(Genomes->same_as[i] == -1) && i!=gindex) {
			dist_mat[gindex*nb_spec + i] = mcdist_noncircular(&genome_list[gindex], &genome_list[i],
															  genome_list[gindex].num_g, genome_list[gindex].num_chr, 
															  distmem, NULL);
			dist_mat[i*nb_spec + gindex] = dist_mat[gindex*nb_spec + i];
		}
		else
			dist_mat[gindex*nb_spec + i] = dist_mat[i*nb_spec + gindex] = 0;
	}

}



int find_total_dist(G_struct *Genomes, int nb_spec, mgr_distmem_t *distmem) 
{
	 int i, j, total = 0;
	 
	 //compute_dist_mat(Genomes, nb_spec, distmem);
	 for (i=0; i<nb_spec; i++)
		  for (j=i+1; j<nb_spec; j++)
				total += Genomes->dist_mat[i*nb_spec + j];
					 
	 return total;
}

void remove_null_chromos (struct mgr_genome_struct *genome_list, int nb_spec,
					  int *num_genes, int *num_chromosomes, int alpha_size, 
                      int max_num_deli) 
{
	int i, j, old_num_genes, k;
	int *new_genes;
	//	fprintf(stdout,"num_genes =%d,num_chromosomes=%d,alpha_size=%d,max_num_deli=%d\n",
    //            (*num_genes),(*num_chromosomes),alpha_size,max_num_deli);
    //    fflush(stdout);
	old_num_genes = *num_genes;
	*num_chromosomes = max_num_deli;  
	*num_genes = alpha_size + max_num_deli*2;	
}


// this procedure adds null chromo at the end of each genome to allow for the max
// number of possible fissions
void add_null_chromos(struct mgr_genome_struct *genome_list, int nb_spec,
					  int *num_genes, int *num_chromosomes) 
{
    
	int i, j, old_num_genes;
	int *new_genes;
		
	old_num_genes = *num_genes;
	*num_genes -= 2*(*num_chromosomes); // this is the actual number of genes without caps
	*num_chromosomes = *num_genes;  
	*num_genes = 3*(*num_genes);
		
	for (i=0; i<nb_spec; i++) {
		 		 
		 new_genes = (int *)e_malloc((*num_genes) * sizeof(int), "new_genes");
		 // this is the real genome
		 for (j=0; j<old_num_genes; j++) 
			  new_genes[j] = genome_list[i].genes[j];
		 
		 // these are all the null chromosomes
		 for (j=old_num_genes; j<*num_genes; j++)
			  new_genes[j] = j+1;
		 
		 free(genome_list[i].genes);
		 genome_list[i].genes = new_genes;
		 
//		 print_genome_multichrom(&genome_list[i], *num_genes, *num_chromosomes);
		 
	}
	
}


void find_max_chromo_size(G_struct *Genomes, int nb_spec) {
	int gindex, c1, size, max_size = 0;
    struct mgr_genome_struct *genome_list = Genomes->genome_list;
	
	if (genome_list[0].num_chr > 0) {

		for (gindex = 0; gindex<nb_spec; gindex++) {
            cbounds_t * cb;
            cb = ( cbounds_t *) malloc ( sizeof( cbounds_t));
            fprintf(stdout, "count_cbounds =%d ", count_cbounds); fflush(stdout);
              count_cbounds ++;
	cb->cNum     = (int *) malloc(genome_list[gindex].num_g*sizeof(int));
	cb->cBound   = (int *) malloc((genome_list[gindex].num_chr+1)*sizeof(int));


            fprintf(stdout,"find max chromo size, gindex=%d \n",gindex);
		//	init_cbounds(genome_list[gindex].num_g, genome_list[gindex].num_chr, &genome_list[gindex], &cb);
        	init_cbounds_wmem(genome_list[gindex].num_g, 
                    genome_list[gindex].num_chr,  
                    &genome_list[gindex], cb);

		
			//fprintf(stdout, "In genome %d, sizes:\n", gindex);
			for (c1 = 1; c1 <= Genomes->nb_chromo[gindex]; c1++) {
				size = cb->cBound[c1]-1 - (cb->cBound[c1-1]+1);
				if (size > max_size) {
					max_size = size;
				}
				//fprintf(stdout, "%d (%d)\n", size, max_size);
			}
            //fprintf(stdout, "free cbounds\n"); fflush(stdout);
            free_cbounds(cb);
		}
		Genomes->max_chromo_size = max_size;
	}
	else {
		Genomes->max_chromo_size = genome_list[0].num_g;
	}
}

void carry_on_reag(list_reag *the_list, G_struct *Genomes, int nb_spec,
				   int *spec_left, int *did_merge, mgr_distmem_t *distmem)
{
	struct mgr_genome_struct *genome_list = Genomes->genome_list;
	//	int *dist_mat = Genomes->dist_mat;
	int s1=0, e1=0, s2=0, e2=0;
	cbounds_t cb;        /* structure with chromosome boundaries */
	int gindex;
	
#if 0
	fprintf(stdout, "we are doing\nin %d do %d (%d %d %d %d)\n", the_list->spec, the_list->optype,
      		the_list->start, the_list->end, the_list->sc1, the_list->sc2);
#endif

	gindex = the_list->spec;	
			
	if (genome_list[0].num_chr>0 && (the_list->sc1<0 || the_list->sc2<0)) {
		/* get the chromosome boundaries */
		fprintf(stdout,"carry_on_reag, "); fflush(stdout);
		init_cbounds(genome_list[gindex].num_g, genome_list[gindex].num_chr,
			&genome_list[gindex],
			&cb);
		
		if (the_list->sc1 < 0) { // we need to flip this chromo
			
			s1 = cb.cBound[(-the_list->sc1)-1];
			e1 = cb.cBound[(-the_list->sc1)]-1;
			
			reverse_in_place_genes(genome_list[gindex].genes, s1,e1);
		}
		if (the_list->sc2 < 0) { // we need to flip this chromo
			
			s2 = cb.cBound[(-the_list->sc2)-1];
			e2 = cb.cBound[(-the_list->sc2)]-1;
			
			reverse_in_place_genes(genome_list[gindex].genes, s2,e2);
		}      
	}
	
	// the actual rearrangement
	reverse_in_place_genes(genome_list[gindex].genes, the_list->start,the_list->end); 
	
	if (the_list->optype == CFUS) {
		// because of null bug
		reformat_caps(&genome_list[gindex], genome_list[gindex].num_g, genome_list[gindex].num_chr);
		Genomes->nb_chromo[gindex]-=1;
	}
	else if (the_list->optype == CFIS) 
		Genomes->nb_chromo[gindex]+=1;
    
	// if this rearrangement created a merge
//	if (*spec_left>3 && the_list->merge_with != -1) {
	if (the_list->merge_with != -1) {
		Genomes->same_as[gindex] = the_list->merge_with;
		(*spec_left)--;
		*did_merge = TRUE;
	}
	else {
		*did_merge = FALSE;
	}
	
	
	// update the distance matrix
	local_distmat_update(gindex, Genomes, nb_spec, distmem);

	if (genome_list[gindex].num_chr>0 && (the_list->sc1<0 || the_list->sc2<0)) 
    {
        //fprintf(stdout," free cb\n "); fflush(stdout);
        free_cbounds(&cb);
    }
	

}
	

void undo_reag(list_reag *the_list, G_struct *Genomes, int nb_spec,
			   int *spec_left, 
			   int did_merge, mgr_distmem_t *distmem) 
{
	// if nb_chromo was modified, undo
	if (the_list->optype == CFUS) 
		Genomes->nb_chromo[the_list->spec]+=1;
	else if (the_list->optype == CFIS) 
		Genomes->nb_chromo[the_list->spec]-=1;

	// if there was a merge, undo
	if (did_merge==TRUE) {
		Genomes->same_as[the_list->spec] = -1;
		(*spec_left)++;
	}
	
	// update the distance matrix
	local_distmat_update(the_list->spec, Genomes, nb_spec, distmem);

}

void print_one_reag(list_reag *the_list, G_struct *Genomes,
					int nb_spec, struct mgr_genome_struct *old_genome)
{
	struct mgr_genome_struct *genome_list = Genomes->genome_list;
	// int *dist_mat = Genomes->dist_mat;
	int s1=0, e1=0, s2=0, e2=0;  /* used to flip chromo if necessary */
	cbounds_t cb;        /* structure with chromosome boundaries */
	int gindex;


#if 0
	fprintf(stdout, "we are trying: %d do %d (%d %d %d %d)\n", the_list->spec, the_list->optype,
      		the_list->start, the_list->end, the_list->sc1, the_list->sc2);
#endif


	gindex = the_list->spec;	

	// make a copy before rearrangement
	copy_genes(genome_list[gindex].genes, old_genome->genes, genome_list[gindex].num_g);
			
	if (genome_list[gindex].num_chr>0) {
		/* get the chromosome boundaries */
		fprintf(stdout,"print_one_reag, ");
		init_cbounds(genome_list[gindex].num_g, genome_list[gindex].num_chr,
			&genome_list[gindex],
			&cb);
		
		if (the_list->sc1 < 0) { // we need to flip this chromo
			s1 = cb.cBound[(-the_list->sc1)-1];
			e1 = cb.cBound[(-the_list->sc1)]-1;
			reverse_in_place_genes(genome_list[gindex].genes, s1, e1);
		}
		else {
			s1 = cb.cBound[the_list->sc1-1];
			e1 = cb.cBound[the_list->sc1]-1;
		}
		
		if (the_list->sc2 < 0) { // we need to flip this chromo
			s2 = cb.cBound[(-the_list->sc2)-1];
			e2 = cb.cBound[(-the_list->sc2)]-1;
			reverse_in_place_genes(genome_list[gindex].genes, s2,e2);
		}
		else {
			s2 = cb.cBound[the_list->sc2-1];
			e2 = cb.cBound[the_list->sc2]-1;
		}
	}

	fprintf(stdout, "in %d do ", gindex);
	switch (the_list->optype) {
	case CREV:
		fprintf(stdout, "reversal from %d to %d\n", genome_list[gindex].genes[the_list->start],
		genome_list[gindex].genes[the_list->end]);
		break;
	case CTRA:
		fprintf(stdout, "translocation from %d to %d\n", genome_list[gindex].genes[the_list->start],
		genome_list[gindex].genes[the_list->end]);
		break;
	case CFIS:
		fprintf(stdout, "fission at %d\n", genome_list[gindex].genes[the_list->start]);
		break;
	case CFUS:
		fprintf(stdout, "fusion between ");
		if (the_list->start == s1+1 && the_list->end == s2)
			fprintf(stdout, "%d and %d\n", -genome_list[gindex].genes[s1+1],
			genome_list[gindex].genes[s2+1]);
		else
			fprintf(stdout, "%d and %d\n", genome_list[gindex].genes[e1-1],
			-genome_list[gindex].genes[e2-1]);
		break;
	}
	//fprintf(stdout, "(%d %d -> %d)\n", the_list->start, the_list->end, the_list->end-the_list->start);

	copy_genes(old_genome->genes, genome_list[gindex].genes, genome_list[gindex].num_g);
	
	if (genome_list[gindex].num_chr>0)
		free_cbounds(&cb);

}

// we want to concatenate the gene at the end of the string after checking
// that there is enough space
void add_to_strip(a_strip *alphabet, int strip, int gene) {
	int i;
	int *new_strip;
	
	//fprintf(stdout, "adding %d\n", gene);
	// check that the size of strip is sufficient
	if (alphabet[strip].last == alphabet[strip].size-1) {

		alphabet[strip].size = (alphabet[strip].size)*2;
		//fprintf(stdout, "the strip_size %d: ", alphabet[strip].size);
		new_strip = (int *)e_malloc((alphabet[strip].size)*sizeof(int), "new_strip in add_to_strip");
		for (i=0; i<= alphabet[strip].last; i++) {
			new_strip[i] = alphabet[strip].strip[i];
		}
		free(alphabet[strip].strip);
		alphabet[strip].strip = new_strip;
	} else {
		if (alphabet[strip].last > alphabet[strip].size-1) {
			fprintf(stdout, "in add_to_strip, last > size-1\n");
			exit(-1);
		}
	}	
	
	alphabet[strip].last = alphabet[strip].last+1;
	alphabet[strip].strip[alphabet[strip].last] = gene;

}

// old_alphabet gives the sequence of genes in a strip
// new_alphabet gives strips of those strips
// we want old_alphabet to give the sequence of genes in the superstrips
void merge_alphabet(a_strip **old_alphabet, a_strip *new_alphabet, int old_nb_strips, int nb_super_strips) {
	int strip;
	a_strip *super_strip;
	int i, j, k;
	
	//fprintf(stdout, "we are merging two alphas... %d %d\n", old_nb_strips, nb_super_strips);
	//print_alphabet(*old_alphabet, old_nb_strips);
	//print_alphabet(new_alphabet, nb_super_strips);
	
	initialize_alphabet(&super_strip, nb_super_strips, FALSE);
/*	super_strip = (a_strip *)e_malloc(nb_super_strips*sizeof(a_strip), "all super_strips in merge_alphabet");
	for (i=0; i<nb_super_strips; i++) {
		super_strip[i].size = MAXSTRIP;
		super_strip[i].strip = (int *)e_malloc(MAXSTRIP*sizeof(int), "super_strip[j]");
		super_strip[i].last = -1;
	}	*/
		
	for (i=0; i< nb_super_strips; i++) {
		
		// for each strip in super_strip
		for (j=0; j<=new_alphabet[i].last; j++) {
			
			// add each gene of strip of super_strip
			strip = new_alphabet[i].strip[j];
			//fprintf(stdout, "for strip %d: ", strip);
			if (strip>=0) {
				for (k=0; k<=(*old_alphabet)[strip-1].last; k++) {
					add_to_strip(super_strip, i, (*old_alphabet)[strip-1].strip[k]);
					//fprintf(stdout, "%d ", (*old_alphabet)[strip-1].strip[k]);
				}
			}
			else {
				for (k=(*old_alphabet)[-strip-1].last; k>=0; k--) {
					add_to_strip(super_strip, i, -(*old_alphabet)[-strip-1].strip[k]);
					//fprintf(stdout, "%d ", -(*old_alphabet)[-strip-1].strip[k]);
				}
			}
			//fprintf(stdout, "\n");
		}
		
	}
	//fprintf(stdout, "superstrip\n");
	//print_alphabet(super_strip, nb_super_strips);

	
	free_alphabet(*old_alphabet, old_nb_strips);
	(*old_alphabet) = super_strip;
	
/*	for (i=0; i<nb_super_strips; i++) {
		free(old_alphabet[i].strip);
		old_alphabet[i].strip = super_strip[i].strip;
		old_alphabet[i].size = super_strip[i].size;
		old_alphabet[i].last = super_strip[i].last;
	}
	free(super_strip); */
}
								


void condense_genomes(G_struct *Genomes, int nb_spec, int verbose) {
	
	struct mgr_genome_struct *genome_list = Genomes->genome_list;
	char *breakpoints;    /* bit array of genes that are breakpoints */
	int lowcap, max_pos;
	int gref = 0; // the reference genome is the first one by default
	int inchromo = FALSE;	
	int i, gene, current_strip, nb_strips = 0;
	int *convert;
	int genomenb, new_num_genes, ori_num_chromosomes, new_num_chromosomes;
	int *new_genes, cap, newpos, oldpos;
	int initial_alpha_size;
	a_strip *new_alphabet;
	
	if (verbose) {
		fprintf(stdout, "\nCondensing...\n\n");
	}
	
	i=0;
	while (i<nb_spec && Genomes->same_as[i]!=-1) {
		i++;
	}
	if (i<nb_spec) {
		gref = i;
	}
	else {
		fprintf(stdout, "no non-empty genome in condense\n");
	}
	
	
	//fprintf(stdout, "num_g num_chr: %d %d\n", genome_list[0].num_g, genome_list[0].num_chr);
	for (i=0;i<nb_spec;i++) {
		if (Genomes->same_as[i] == -1) {
			if (genome_list[i].num_g != genome_list[gref].num_g ||
				genome_list[i].ori_num_g != genome_list[gref].ori_num_g ||
				genome_list[i].num_chr != genome_list[gref].num_chr ) {
				fprintf(stdout, "problem...\n");
				fprintf(stdout, "%d %d: %d %d %d %d %d %d\n", gref, i, genome_list[gref].num_g, genome_list[i].num_g,
						genome_list[gref].ori_num_g, genome_list[i].ori_num_g, genome_list[gref].num_chr, genome_list[i].num_chr);
				exit(-1);
			}
		}
	}
	
	// two stages, first find new alphabet
	/* the lowest cap number */
	if (genome_list[gref].num_chr > 0) {
		max_pos = genome_list[gref].num_g/3 + 2*Genomes->nb_chromo[gref];	// to avoid all the nulls at the end		
		lowcap = genome_list[gref].num_g - 2*(genome_list[gref].num_chr) + 1; 
		initial_alpha_size = genome_list[gref].num_g/3;
	}
    else {
		max_pos = genome_list[gref].num_g;
		lowcap = genome_list[gref].num_g + 1;
		initial_alpha_size = genome_list[gref].num_g;
	}
	
	// allocate memory for conversion table
	convert = (int *) e_malloc(lowcap*sizeof(int), "conversion table");
	
	// allocate memory for breakpoints
    init_isBP_mem(genome_list[gref].num_g, genome_list[gref].num_chr, &breakpoints);
	init_isBP_array(genome_list[gref].num_g, genome_list[gref].num_chr, nb_spec, // only consider real genome for finding breakpoints
					Genomes, gref, //use the first genome as the reference
					breakpoints);

	// first step is to identify strips and count them
	//fprintf(stdout, "looking for strips\n");
	for (i=0; i < max_pos; i++) {
		gene = genome_list[gref].genes[i];
		if (abs(gene) >= lowcap) {
			inchromo=!inchromo;
			if (!inchromo && i<max_pos-1) {
				//fprintf(stdout, "\n");
				nb_strips++;
			}
		}
		else {
			if (i>0 && abs(genome_list[gref].genes[i-1]) < lowcap && //not the first gene of chromosome
				isBP(gene, breakpoints)) {
				//fprintf(stdout, "\n");
				nb_strips++;
			}
			//fprintf(stdout, "%d ", gene);
		}
	}
	//fprintf(stdout, "\n");
	nb_strips++;
	
	if (nb_strips < initial_alpha_size) {
		if (verbose) {
			fprintf(stdout, "found %d strips\n", nb_strips);
		}
		
		// generate new_alphabet for these strings
		initialize_alphabet(&new_alphabet, nb_strips, FALSE);
		/*new_alphabet = (a_strip *) e_malloc(nb_strips*sizeof(a_strip), "new alphabet");
		for (j=0; j<nb_strips; j++) {
			new_alphabet[j].size = MAXSTRIP;
			new_alphabet[j].strip = (int *)e_malloc(MAXSTRIP*sizeof(int), "new_alphabet[j]");
			new_alphabet[j].last = -1;
		}*/
		
		current_strip = 0;
		// store new alphabet
		for (i=0; i < max_pos; i++) {
			gene = genome_list[gref].genes[i];
			if (abs(gene) >= lowcap) {
				inchromo=!inchromo;
				if (!inchromo && i<max_pos-1) {
					current_strip++;
				}
			}
			else {
				if (i>0 && abs(genome_list[gref].genes[i-1]) < lowcap && //not the first gene of chromosome
					isBP(gene, breakpoints)) {
					current_strip++;
				}
				if (gene > 0) {
					convert[gene] = current_strip+1;
				}
				else {
					convert[-gene] = -(current_strip+1);
				}
				add_to_strip(new_alphabet, current_strip, gene);
			}
		}
		
		if (verbose) {
			print_alphabet(new_alphabet, nb_strips);
		}
		
		// we want to pad our genomes with null chromosomes at the end to allow for fissions
		if (genome_list[gref].num_chr > 0) {
			new_num_genes = 3*nb_strips;
			new_num_chromosomes = nb_strips;
		}
		else {
			new_num_genes = nb_strips;
			new_num_chromosomes = 0;
		}
		
		for (genomenb=0; genomenb<nb_spec; genomenb++) {
			
			// only condense unmerged genomes
			if (Genomes->same_as[genomenb] == -1) {
				
				cap = nb_strips+1;
			/* allocate space to store the new genes */
				new_genes = (int *) e_malloc(new_num_genes*sizeof(int), "new genes");
			
			// scan old genes and update new genes
				newpos = oldpos = 0;
				while (newpos < new_num_genes) {
					gene = genome_list[genomenb].genes[oldpos];
				
					if (abs(gene) >= lowcap) {
					//we are reading a cap, update new sequence
						new_genes[newpos] = cap;
						cap++;
						newpos++;
					}
					else {
						if (newpos > 0 && abs(convert[abs(gene)]) == abs(new_genes[newpos-1])) {
						// we are another gene from the same strip
							oldpos++;
							continue;
						}
						else {
						// we are starting a new strip
							if (gene > 0) {
								new_genes[newpos] = convert[gene];
							}
							else {
								new_genes[newpos] = -convert[-gene];
							}
							newpos++;
						}
					}
				//fprintf(stdout, "oldpos gene, newpos newgene: %d %d, %d %d\n", oldpos, gene, newpos-1, new_genes[newpos-1]);
					oldpos++;
				}
			
			
				free(genome_list[genomenb].genes);
				genome_list[genomenb].genes = new_genes;
			
			
			// update alphabet
				//if (verbose) {
				//	fprintf(stdout, "genomenb %d (%d)\n", genomenb, gref);
				//}
				merge_alphabet(&(genome_list[genomenb].alphabet), new_alphabet, initial_alpha_size, nb_strips);
				genome_list[genomenb].num_g = new_num_genes;
				genome_list[genomenb].num_chr = new_num_chromosomes;
			//genome_list[genomenb].num_ori remains unchanged
			
			//print_alphabet(genome_list[genomenb].alphabet, nb_strips);
			
		
			//print_genome_multichrom_wcaps(&genome_list[genomenb], new_num_genes, new_num_chromosomes);
			}
		}
		free_alphabet(new_alphabet, nb_strips);
		
		if (verbose) {
			fprintf(stdout, "nb_strips %d\n", nb_strips);
			print_alphabet(genome_list[gref].alphabet, nb_strips);
		}
		if (genome_list[gref].num_chr>0) {
			ori_num_chromosomes = genome_list[gref].ori_num_g/3;
		}
		else {
			ori_num_chromosomes = 0;
		}

	}
	else {
		if (verbose) {
			fprintf(stdout, "No condensing necessary\n");
		}
	}
			
	free(convert);
	clean_isBP_mem(breakpoints);
	
	if (verbose) {
		special_print_genomes2(stdout, Genomes,  nb_spec);
	}
	//fprintf(stdout, "Leaving condensing\n");
	
	//print_G_struct(Genomes, nb_spec);
	
}


void uncondense_genomes(G_struct *Genomes, int nb_spec, int verbose) {
	
	struct mgr_genome_struct *genome_list = Genomes->genome_list;
	int genomenb, oldpos, newpos, lowcap, cap;
	int strip, j;
	int *new_genes;
	int current_num_genes, current_num_chromosomes;
			
	for (genomenb=0; genomenb<nb_spec; genomenb++) {
		
		current_num_genes = genome_list[genomenb].num_g;
		
		if (current_num_genes != genome_list[genomenb].ori_num_g) {
			// this genome is condensed, need to uncondense
			
			//fprintf(stdout, "\n\nUncondensing... %d %d\n", current_num_genes, genome_list[genomenb].ori_num_g);
			
			if (verbose) {
				if (genomenb==0)
					fprintf(stdout, "\n\nUncondensing...\n");
			}
			
			
			if (genome_list[genomenb].num_chr>0) {
				current_num_chromosomes = current_num_genes/3;
				lowcap = current_num_genes - 2*(current_num_chromosomes) + 1;
			}
			else {
				current_num_chromosomes = 0;
				lowcap = current_num_genes + 1;
			}
			
			//print_genome_multichrom_wcaps(&genome_list[genomenb], genome_list[genomenb].num_g, genome_list[genomenb].num_chr);
			//print_alphabet(genome_list[genomenb].alphabet, current_num_genes);
			
			cap = (genome_list[genomenb].ori_num_g/3)+1;
			/* allocate space to store the new genes */
			new_genes = (int *) e_malloc(genome_list[genomenb].ori_num_g*sizeof(int), "new genes");

			
			// scan old genes and update new genes
			newpos = oldpos = 0;
			while (oldpos < current_num_genes) {
				strip = genome_list[genomenb].genes[oldpos];
				//fprintf(stdout, "reading strip %d (%d)\n", strip, lowcap);

				if (abs(strip) >= lowcap) {
					//we are reading a cap, update new sequence
					//fprintf(stdout, "inserting %d\n", cap);
					new_genes[newpos++] = cap++;
				}
				else {
					// we are reading a strip, need to insert all genes of this strip
					if (strip >= 0) {
						for (j=0; j<=genome_list[genomenb].alphabet[strip-1].last; j++) {
							new_genes[newpos++] = genome_list[genomenb].alphabet[strip-1].strip[j];
							//fprintf(stdout, "inserting %d\n", genome_list[genomenb].alphabet[strip-1].strip[j]);
						}
					}
					else {
						for (j=genome_list[genomenb].alphabet[-strip-1].last; j>=0; j--) {
							new_genes[newpos++] = -genome_list[genomenb].alphabet[-strip-1].strip[j];
							//fprintf(stdout, "inserting %d\n",  -genome_list[genomenb].alphabet[-strip-1].strip[j]);
						}
					}
				}
				oldpos++;
			}
			// fill with null chromosomes
			while (newpos < genome_list[genomenb].ori_num_g) {
				//fprintf(stdout, "inserting %d\n", cap);
				new_genes[newpos++] = cap++;
			}
			
			
			free(genome_list[genomenb].genes);
			genome_list[genomenb].genes = new_genes;
			genome_list[genomenb].num_g = genome_list[genomenb].ori_num_g;
			
			// reset alphabet
			//fprintf(stdout, "size alphabet before/after %d -> %d\n", current_num_genes, genome_list[genomenb].ori_num_g);
			if (genome_list[genomenb].num_chr>0) {
				free_alphabet(genome_list[genomenb].alphabet, current_num_genes/3);
				initialize_alphabet(&(genome_list[genomenb].alphabet), genome_list[genomenb].ori_num_g/3, TRUE);
				genome_list[genomenb].num_chr = genome_list[genomenb].num_g/3;
			}
			else {
				free_alphabet(genome_list[genomenb].alphabet, current_num_genes);
				initialize_alphabet(&(genome_list[genomenb].alphabet), genome_list[genomenb].ori_num_g, TRUE);
				genome_list[genomenb].num_chr = 0;
			}				
			/*if (genome_list[genomenb].num_chr > 0) {
				alphabet_size = genome_list[genomenb].ori_num_g/3;
			}
			else {
				alphabet_size = genome_list[genomenb].ori_num_g;
			}
			for (i=0; i< alphabet_size; i++) {
				genome_list[genomenb].alphabet[i].strip[0] = i+1;
				genome_list[genomenb].alphabet[i].last = 0;
			}*/
			
			//print_genome_multichrom_wcaps(&genome_list[genomenb], new_num_genes, new_num_chromosomes);
		}
	}
}
	
