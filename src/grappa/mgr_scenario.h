/* scenario.h
*    Generate sequence of steps of a most parsimonious scenario
*    (original version; some routines still used, but see opt_scenario.{c,h}
	  *    for improved version)
*
* Copyright (C) 2001-2006 The Regents of the University of California
* by Glenn Tesler
*
* See file COPYRIGHT for details.
*****************************************************************************
* This file is part of GRIMM.
*
* GRIMM is free software; you can redistribute it and/or modify
* it under the terms of the GNU General Public License, Version 2,
* dated June 1991, as published by the Free Software Foundation.
*
* GRIMM is distributed in the hope that it will be useful, but
* WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
* General Public License for more details.
*
* You should have received a copy of the GNU General Public License along
* with this program; if not, write to the Free Software Foundation, Inc.,
* 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
*/

/* Last modified on Tue Aug 1, 2006, by Glenn Tesler
*/

/* 13 -2 14 15 1 3 4 5 16 
   13 15 are LTAIL, 14 16 are RTAIL
   1 is LEND, 5 is REND, -2 is LREND
   3 4 are NORMAL
   */
#define T_NORMAL 0
#define T_LTAIL 1
#define T_RTAIL 2
#define T_LEND 3
#define T_REND 4
#define T_LREND 5
#define T_ERROR 6

void print_scenario(struct mgr_genome_struct *g1, struct mgr_genome_struct *g2,
		    int num_genes, int num_chromosomes);

void copy_genes(int *source_perm, int *dest_perm, int num_genes);

void copy_and_reverse_genes(int *source_perm, int *dest_perm,
			    int start, int end,
			    int num_genes);

void reverse_in_place_genes(int *source_perm,
			    int start, int end);

void print_genes(int *source_perm, int num_genes);
void print_genes_nocr(int *source_perm, int num_genes);

int position_type(int pos, int *genes, int num_genes, int cap);

int chromosome_number(int pos, int *genes, int num_genes, int cap);

int is_empty_chromosome(int *genes, int ch_num, int num_genes, int cap);

int is_null_chromosome(int pos, int *genes, int num_genes, int cap);


void alloc_simple_genome(int num_genes,
			 struct mgr_genome_struct *genome);

void free_simple_genome(struct mgr_genome_struct *genome);

void init_scenario_mem(int num_genes, int num_chromosomes, mgr_distmem_t *distmem,
		       struct mgr_genome_struct *current_genome,
		       struct mgr_genome_struct *trial_genome,
		       struct mgr_genome_struct *g1,
		       int pad);

void init_scenario_bp_wmem(int num_genes, int num_chromosomes,
			   /* mgr_distmem_t *distmem, */
			   struct mgr_genome_struct *g2,
			   int pad,
			   int **successors);

void init_scenario_bp(int num_genes, int num_chromosomes,
		      /* mgr_distmem_t *distmem, */
		      struct mgr_genome_struct *g2,
		      int pad,
		      int **successors);

void clean_scenario_bp_mem(int *successors);

void clean_scenario_mem(struct mgr_genome_struct *current_genome,
			struct mgr_genome_struct *trial_genome,
			mgr_distmem_t *distmem);



typedef struct {
  /* these arrays use positions in the signed permutation,
     not the doubled unsigned one that is used
     in distmem->chromNum,chromBd
     Also cBound is indexed differently */

  int *cNum;           /* array mapping gene position to chromosome # */

  int *cBound;         /* array mapping chromosome # to left end of chrom */
  /* chromosome i=1,2,... is in positions
     cBound[i-1], cBound[i-1]+1, ..., cBound[i]-1
     */

     
  
} cbounds_t;


#define SUCCESSOR(k) successors[ (k)>0 ? (k)-1 : -(k)-1+num_genes+2 ]
#define ISUCCESSOR(k) (*successors)[ (k)>0 ? (k)-1 : -(k)-1+num_genes+2 ]
