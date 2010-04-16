/* my_testrev.h
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
#include "mgr_scenario.h"


Breakpoint_array bkarr_build_list_reag;

list_reag list_build_list_reag;


void init_cbounds_wmem(int num_genes, int num_chromosomes,
					   struct mgr_genome_struct *g1,
					   cbounds_t *cb);
void init_cbounds(int num_genes, int num_chromosomes,
				  struct mgr_genome_struct *g1,
				  cbounds_t *cb);

void free_cbounds(cbounds_t *cb);



void init_isBP_mem(int num_genes, int num_chromosomes,
				   char **breakpoints);
void clean_isBP_mem(char *breakpoints);

void init_isBP_array(int num_genes, int num_chromosomes,
					 int num_genomes,
					 G_struct *Genomes,
					 int gindex, /* which genome to compute b.p. against */
					 char *breakpoints);

int isBP(int g, char *breakpoints);




/* private declarations */

#define BP_FIS     " /"
#define BP_BP      " :"
#define BP_ABFUS     " | + |"
#define BP_ABNOFUS   " | . |"
#define BP_BAFUS     " | +"
/* #define BP_BANOFUS   " | ." */
#define BP_BANOFUS   ""


#define CREV 0         /* counts for reversals */
#define CTRA 1         /* counts for translocations */
#define CFIS 2         /* fissions */
#define CFUS 3         /* fusions */
#define CTOT0 4        /* total */
#define CNF0 5         /* number of fields */

#define CDOWN1 0       /* events lowering distance by 1 */
#define CSAME  1       /* events keeping distance the same */
#define CUP1   2       /* events raising distance by 1 */
#define CERROR 3       /* another amount, error */
#define CTOT1  4       /* total */
#define CFIRST 5       /* count # tries till first success */
#define CNF1   6       /* number of fields */

cbounds_t  cb_build_list_reag; //mgr_my_testrev.c, "build_list_reag"

int build_list_reag(list_listreag * in_list,//list_reag **a_list, 
        G_struct *Genomes , int nb_spec, int spec_left,
					int reversals,      // TRUE if we want to include reversals
					int transloc, 	// TRUE if we want to include translocations
					int fusion,
					int fission, // TRUE if we want to include fusions and fissions
					mgr_distmem_t *distmem,
					int only_size,
					int optimize,
					int reduction_type,
					int rev_size,
					int last_genomenb,
					int pair1, int pair2,
					int verbose);



