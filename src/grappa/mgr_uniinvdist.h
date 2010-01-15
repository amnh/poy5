/* uniinvdist.h
*    1-line description
*
* Copyright (C) 2001-2006 The Regents of the University of California
* by Glenn Tesler
*
* Contains code from invdist.h in GRAPPA 1.02
* Copyright (C) 2000-2001  The University of New Mexico and
*                          The University of Texas at Austin
* by David A. Bader, Bernard M.E. Moret, Tandy Warnow, Stacia K Wyman, Mi Yan
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

/* Modified from GRAPPA 1.02: invdist.h */

#ifndef UNIINVDIST_H
#define UNIINVDIST_H

#include <sys/time.h>
/*#include "structs.h"*/

/*
extern double time_linear;
#ifdef BH
extern double time_BH;
#endif // BH
*/
/*
void calc_invmatrix(struct mgr_genome_struct *genomes, int num_genes, int num_genomes,
		    mgr_distmem_t *distmem, int CIRCULAR);
*/

#ifdef BH
void calc_invmatrix_BH(struct mgr_genome_struct *genomes, int num_genes, int num_genomes,
		       mgr_distmem_t *distmem, int CIRCULAR);
#endif // BH
/*
void setinvmatrix(int **distmatrix,struct mgr_genome_struct *genomes,
		  int num_genes,int num_genomes,mgr_distmem_t *distmem,int CIRCULAR);
*/
int invdist_noncircular_G(struct mgr_genome_struct *g1,
			  struct mgr_genome_struct *g2,
			  int offset, int num_genes,
			  mgr_distmem_t *distmem,
			  graph_t *G,
			  pathcounts_t *pcounts_G,
			  graphstats_t *graphstats);

int mgr_invdist_noncircular(struct mgr_genome_struct *g1, struct mgr_genome_struct *g2,
			int offset, int num_genes, mgr_distmem_t *distmem);

int invdist_noncircular_v(struct mgr_genome_struct *g1, struct mgr_genome_struct *g2,
			  int offset, int num_genes, mgr_distmem_t *distmem,
			  graphstats_t *graphstats);
/*
int invdist_circular(struct mgr_genome_struct *g1, struct mgr_genome_struct *g2,
		     int num_genes, mgr_distmem_t *distmem);
             */
/*
#ifdef BH
int invdist_noncircular_BH(struct mgr_genome_struct *g1, struct mgr_genome_struct *g2,
			   int offset, int num_genes, mgr_distmem_t *distmem);
int invdist_circular_BH(struct mgr_genome_struct *g1, struct mgr_genome_struct *g2,
			int num_genes, mgr_distmem_t *distmem);
#endif // BH
*/
int mgr_invdist_noncircular_nomem(struct mgr_genome_struct *g1,
			      struct mgr_genome_struct *g2,
			      int offset, int num_genes);
int invdist_noncircular_nomem_v(struct mgr_genome_struct *g1,
				struct mgr_genome_struct *g2,
				int offset, int num_genes,
				graphstats_t *graphstats);
int mgr_invdist_circular_nomem(struct mgr_genome_struct *g1,
			   struct mgr_genome_struct *g2,
			   int num_genes);

int mgr_calculate_offset(struct mgr_genome_struct *g1, struct mgr_genome_struct *g2,
		     int num_genes);

#endif
