/* list_ops.c
*
*  SUMMARY:
*
*		All the basic operations associated list of both rearrangements
*		and of edges of the tree.
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
#include "mgr_check.h"
#include "mgrstructs.h"
#include "mgr_genome_ops.h"
#include "mgr_my_testrev.h"
#include "mgr_list_ops.h"

list_edge *RemoveItem_edge(list_edge *the_list) {
    list_edge *temp;
    
    temp = the_list->next;
    free(the_list);
    return temp;
}

void erase_list_edge(list_edge *the_list) {
	while (the_list!=NULL) {
		the_list = RemoveItem_edge(the_list);   
	}
}


list_reag *RemoveItem_reag(list_reag *the_list) {
    list_reag *temp;
    
	free(the_list->changes);
    temp = (*the_list).next;
    free(the_list);
    return temp;
}


void erase_list_reag(list_reag *the_list) {
	while (the_list!=NULL) {
		the_list = RemoveItem_reag(the_list);   
	}
}


void print_list_edge(list_edge *the_list) {
	
	fprintf(stdout,"The list of edges:\n\n");
	fprintf(stdout,"node1\t node2\t score\n");
	while(the_list!=NULL) {
		fprintf(stdout,"  %d\t   %d\t   %d\n",
				  the_list->label1, the_list->label2,
 					(*the_list).score);
		
		the_list=(*the_list).next;
	}
	fprintf(stdout,"----------------------\n");
}

void merge_two_list(treemem_t *treemem) {
	 list_edge *temp_edge_list;
	 
	 temp_edge_list = treemem->the_edge_list;
	 while (temp_edge_list->next != NULL)
	 	temp_edge_list = temp_edge_list->next;
	 temp_edge_list->next = treemem->merge_edge_list;
	 
	 treemem->merge_edge_list = NULL;
	 
	 return ;
}




list_edge *list_edge_insert(list_edge *the_list, int node1, int label1,
									 int node2, int label2, int score) {	
    list_edge *temp_list = the_list;
	
    if (the_list!=NULL) {
		temp_list = the_list;
		the_list = (list_edge *) ckalloc (sizeof (list_edge));
		the_list->next = temp_list;
		the_list->node1=node1;
		the_list->label1=label1;
		the_list->node2=node2;
		the_list->label2=label2;
		the_list->score=score;
		
		return the_list;
    }
    else {
		the_list = (list_edge *) ckalloc (sizeof (list_edge));
		the_list->next = NULL;
		the_list->node1=node1;
		the_list->label1=label1;
		the_list->node2=node2;
		the_list->label2=label2;
		the_list->score=score;
		
		return the_list;
    }
}

list_reag *list_reag_insert(list_reag *the_list, int spec, int optype, int start, int end,
							int sc1, int sc2, int merge_with, int reduction, int *changes, int nb_spec) {
	int i;
    list_reag *temp_list = the_list;
		
    if (the_list!=NULL) {
		temp_list = the_list;
		the_list = (list_reag *) ckalloc (sizeof (list_reag));
		the_list->next = temp_list;
		the_list->spec=spec;
		the_list->optype=optype;
		the_list->start=start;
		the_list->end=end;
		the_list->sc1=sc1;
		the_list->sc2=sc2;
		the_list->merge_with=merge_with;
		the_list->reduction = reduction;
		the_list->Rscore = 0;
		the_list->changes = (int *) ckalloc (nb_spec*sizeof(int));
		for (i=0; i<nb_spec; i++) {
			the_list->changes[i] = changes[i];
		}
		
		return the_list;
    }
    else {
		the_list = (list_reag *) ckalloc (sizeof (list_reag));
		the_list->next = NULL;
		the_list->spec=spec;
		the_list->optype=optype;
		the_list->start=start;
		the_list->end=end;
		the_list->sc1=sc1;
		the_list->sc2=sc2;
		the_list->merge_with=merge_with;
		the_list->reduction = reduction;
		the_list->Rscore = 0;
		the_list->changes = (int *) ckalloc (nb_spec*sizeof(int));
		for (i=0; i<nb_spec; i++) {
			the_list->changes[i] = changes[i];
		}
		
		return the_list;
    }
}


void print_list_reag(list_reag *the_list, int nb_spec) {
    int i, j;
    
    i=1;
    while (the_list!=NULL) {
	
		fprintf(stdout, "%d> in %d do %d (%d %d %d %d) (reduction=%d, merge=%d, Rscore=%d) (", i++, 
				the_list->spec, the_list->optype,
				the_list->start, the_list->end, the_list->sc1, the_list->sc2,
				the_list->reduction, the_list->merge_with, the_list->Rscore);
		for (j=0; j<nb_spec; j++) {
			fprintf(stdout, "%d ", the_list->changes[j]);
		}
		fprintf(stdout, ")\n");
	the_list = the_list->next;
    }
}

int rear_length(list_reag *the_list) {
	int length;
	
	if (the_list != NULL) {
		length = the_list->end - the_list->start;
		
		return abs(length);
	}
	else {
		return MAXINT;
	}
}

int total_reduction(list_reag *the_list) {
	if (the_list!=NULL)
		return total_reduction(the_list->next)+the_list->reduction;
	else return 0;
}

int total_edge_score(list_edge *the_list) {
	if (the_list!=NULL)
		return total_edge_score(the_list->next)+the_list->score;
	else return 0;
}
  
  
// adds the number of good reversals we did on each genome initially
void  add_nbrev_edges(list_edge *the_edge_list, int *nbrev, 
							 int nb_spec) {
	list_edge *temp_edge_list;
	
	temp_edge_list = the_edge_list;
	while (temp_edge_list != NULL) {
		if (temp_edge_list->node1 < nb_spec)
			temp_edge_list->score += nbrev[temp_edge_list->node1];
		else
			if (temp_edge_list->node2 < nb_spec)
				temp_edge_list->score += nbrev[temp_edge_list->node2];
			temp_edge_list = temp_edge_list->next;
	}
}

					
void check_list_edge(G_struct *Genomes, list_edge *the_edge_list, mgr_distmem_t *distmem,
					 int verbose) 
{
	struct mgr_genome_struct *genome_list = Genomes->genome_list;	 
	 
	 while (the_edge_list!=NULL) {
		  int edge_dif;
		 
		  edge_dif = mcdist_noncircular(&genome_list[the_edge_list->label1],
										&genome_list[the_edge_list->label2], genome_list[the_edge_list->label1].num_g, 
										genome_list[the_edge_list->label1].num_chr,
										distmem, NULL) - the_edge_list->score;
		  

		  if (edge_dif!=0) {
			  if (verbose) {
				  fprintf(stdout, "edge_dif is %d between %d and %d\n", edge_dif, 
						  the_edge_list->label1, the_edge_list->label2);
			  }
		      if (edge_dif < 0) // we can improve the tree by putting the actual distance
					the_edge_list->score+=edge_dif;
			  else {
			  	 fprintf(stderr, "Problem in list_ops.c in check_list_edge\n");
			  	 exit(-1);
			  }
		  }
		  
				
		  the_edge_list = the_edge_list->next;
		  
	 }
	 
	
	 
	 fprintf(stdout, "\n\n");
}

void compute_Rscore(list_reag *the_list, G_struct *Genomes) {
	struct mgr_genome_struct *genome_list = Genomes->genome_list;
	int *counter, gindex, i;
	cbounds_t cb;        /* structure with chromosome boundaries */
	list_reag *first_reag, *temp_list;
	int s1,e1,s2,e2,position;
	
	if (the_list != NULL) {
	// allocate space for counter
		temp_list = the_list;
		
		counter = (int *) ckalloc(genome_list[0].num_g * sizeof(int));
	
	    // we count the number of occurences of each edge independently for each genome
		while (temp_list != NULL) {
			
			//reset counter for each genome
			for (i=0;i<genome_list[0].num_g;i++) {
				counter[i] = 0;
			}
			
			first_reag = temp_list;
			gindex = temp_list->spec;
			
			if (genome_list[0].num_chr>0)
				init_cbounds(genome_list[gindex].num_g, genome_list[gindex].num_chr,
							 &genome_list[gindex],
							 &cb);
			
			while (temp_list!=NULL && temp_list->spec == gindex) {
				
				// for each rearrangement, we want to update edge counter
				//counter[temp_list->start]++;
				//counter[temp_list->end]++;


				if (temp_list->sc1 < 0) {
						// need to reverse chromosome, broken edge will be after current start
					s1 = cb.cBound[(-temp_list->sc1)-1];
					e1 = cb.cBound[(-temp_list->sc1)]-1;
					position = temp_list->start - s1;
					position = e1 - position;
					counter[position+1]++;
				}
				else {
					s1 = cb.cBound[temp_list->sc1-1];
					e1 = cb.cBound[temp_list->sc1]-1;
					counter[temp_list->start]++;
				}
				if (temp_list->sc2 < 0) {
						// need to reverse chromosome, broken edge will be before current end
					s2 = cb.cBound[(-temp_list->sc2)-1];
					e2 = cb.cBound[(-temp_list->sc2)]-1;
					position = temp_list->end - s2;
					position = e2 - position;
					counter[position]++;
				}
				else {
					s2 = cb.cBound[temp_list->sc2-1];
					e2 = cb.cBound[temp_list->sc2]-1;
					counter[temp_list->end+1]++;
				}
								
				temp_list = temp_list->next;								
				
			}

    		// go back to first rearrangement and update Rscore
			temp_list = first_reag;
			while (temp_list!=NULL && temp_list->spec == gindex) {
				
				if (temp_list->sc1 < 0) {
						// need to reverse chromosome, broken edge will be after current start
					s1 = cb.cBound[(-temp_list->sc1)-1];
					e1 = cb.cBound[(-temp_list->sc1)]-1;
					position = temp_list->start - s1;
					position = e1 - position;
					temp_list->Rscore = counter[position+1]; 
				}
				else {
					s1 = cb.cBound[temp_list->sc1-1];
					e1 = cb.cBound[temp_list->sc1]-1;
					temp_list->Rscore = counter[temp_list->start];
				}
				if (temp_list->sc2 < 0) {
						// need to reverse chromosome, broken edge will be before current end
					s2 = cb.cBound[(-temp_list->sc2)-1];
					e2 = cb.cBound[(-temp_list->sc2)]-1;
					position = temp_list->end - s2;
					position = e2 - position;
					temp_list->Rscore += counter[position];
				}
				else {
					s2 = cb.cBound[temp_list->sc2-1];
					e2 = cb.cBound[temp_list->sc2]-1;
					temp_list->Rscore += counter[temp_list->end+1];
				}
				
				temp_list = temp_list->next;
			}
			
			if (genome_list[0].num_chr>0)
				free_cbounds(&cb);

		}
		
		free(counter);
	}
}
 
