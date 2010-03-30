/* mcrdist.h
*    Pairwise distance between two multichromosomal genomes.
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


struct mgr_genome_struct trial_genome_build_list_reag;

void buildG_grey(int num_genes, int num_chromosomes, mgr_distmem_t *distmem,
		 graph_t *G);

void repairG_chrom(int num_genes, int num_chromosomes, mgr_distmem_t *distmem,
		   graph_t *G);

int mcdist_noncircular_nomem(struct mgr_genome_struct *g1,
			     struct mgr_genome_struct *g2,
			     int num_genes, int num_chromosomes,
			     graphstats_t *graphstats);

int mcdist_noncircular(struct mgr_genome_struct *g1, struct mgr_genome_struct *g2,
		       int num_genes, int num_chromosomes,
		       mgr_distmem_t *distmem,
		       graphstats_t *graphstats);


int mcdist_capgraph_nomem(struct mgr_genome_struct *g1,
			  struct mgr_genome_struct *g2,
			  int num_genes, int num_chromosomes,
			  graphstats_t *graphstats);

int mcdist_capgraph(struct mgr_genome_struct *g1, struct mgr_genome_struct *g2,
		    int num_genes, int num_chromosomes,
		    mgr_distmem_t *distmem,
		    graphstats_t *graphstats);

int mcdist_capgraph_connections(struct mgr_genome_struct *g1,
				struct mgr_genome_struct *g2,
				int num_genes, int num_chromosomes,
				mgr_distmem_t *distmem,
				graph_t *G, pathcounts_t *pcounts_G,
				graphstats_t *graphstats);
void mcdist_capgraph_finalperms(struct mgr_genome_struct *g1,
			       struct mgr_genome_struct *g2,
			       int num_genes, int num_chromosomes,
			       mgr_distmem_t *distmem,
			       graph_t *G,
			       graphstats_t *graphstats);

#ifdef NULLFIX1
void kludge_fix_nulls(struct mgr_genome_struct **g1, struct mgr_genome_struct **g2,
		      int *num_genes, int *num_chromosomes);
#else
#define kludge_fix_nulls(g1,g2,num_genes,num_chromosomes)
#endif






void mcdist_allocmem(int num_genes, int num_chromosomes,
		     mgr_distmem_t *distmem);
void mcdist_freemem(mgr_distmem_t *distmem);


void setmcdistmatrix(int **distmatrix, struct mgr_genome_struct *genomes,
		     int num_genes, int num_chromosomes, int num_genomes,
		     mgr_distmem_t *distmem);


void mcdist_noncircular_graph_d(struct mgr_genome_struct *g1,
				struct mgr_genome_struct *g2,
				int num_genes, int num_chromosomes,
				mgr_distmem_t *distmem);

void num_breakpoints_mc(int *perm,
			int size,
			int lowcap1,  /* lowest cap in non-doubled perms */
			struct mgr_genome_struct *g1,
			struct mgr_genome_struct *g2,

			/* return values: */
			int *bp_int,   /* # internal breakpoints */
			int *bp_ext    /* # external breakpoints */
			);

void mcdist_partial_G(struct mgr_genome_struct *g1,
		      struct mgr_genome_struct *g2,
		      int num_genes, int num_chromosomes,
		      mgr_distmem_t *distmem,
		      graph_t *G,
		      pathcounts_t *pcounts_G);

