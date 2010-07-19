/* uniinvdist.c
*    Unichromosomal reversal distance algorithm.
*
* Copyright (C) 2001-2006 The Regents of the University of California
* by Glenn Tesler
*
* Contains code from invdist.c in GRAPPA 1.02
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

/* Unichromosomal reversal distance algorithm */
/* This implements parts of the externally accessible interface to
   the original invdist.c, but some of the innards have been
   enhanced for the multichromosomal case and moved to other
   files. */

#include <stdio.h>
#include <stdlib.h>
#include "mgrstructs.h"
#include "mgr_uniinvdist.h"
#include "mgr_graph_edit.h"
#include "mgr_graph_components.h"
#include "mgr_mcrdist.h"

/* same as in invdist.c
#ifdef TESTING
double time_linear, time_BH;
#endif
*/



/* same as in invdist.c 

void calc_invmatrix(struct mgr_genome_struct *genomes, int num_genes, int num_genomes,
		     mgr_distmem_t *distmem, int CIRCULAR) {
  int i, j;
  int dist;

  for (i=0 ; i<num_genomes ; i++) {
    for (j=i+1 ; j<num_genomes ; j++) {
      if (CIRCULAR)
	dist = invdist_circular(genomes+i, genomes+j, num_genes, distmem);
      else 
	dist = invdist_noncircular(genomes+i, genomes+j, 0, num_genes, distmem);
    }
  }
    
  return;
}
*/

/* same as in invdist.c 
void setinvmatrix(int **distmatrix,struct mgr_genome_struct *genomes,
		  int num_genes,int num_genomes,mgr_distmem_t *distmem,
		  int CIRCULAR) {
  int i, j;

  for (i=0 ; i<num_genomes ; i++) {
    distmatrix[i][i] = 0;
    for (j=i+1 ; j<num_genomes ; j++) {
      if (CIRCULAR)
	distmatrix[i][j] = distmatrix[j][i] =
		invdist_circular(genomes+i,genomes+j,num_genes,distmem);
      else 
	distmatrix[i][j] = distmatrix[j][i] =
		invdist_noncircular(genomes+i,genomes+j,0,num_genes,distmem);
    }
  }
  return;
}
*/

/* same as invdist.c */
int mgr_calculate_offset(struct mgr_genome_struct *g1, struct mgr_genome_struct *g2,
		     int num_genes) {
  int i;
  int gA;
  int *genes2;

#ifdef DEBUG
  for (i=0 ; i<num_genes ; i++) {
    fprintf(stdout,"calc_offset: g1[%3d]: %3d  g2[%3d]: %3d\n",
	    i, g1->genes[i], i, g2->genes[i]);
  }
#endif
  
  genes2 = g2->genes;

  gA    = g1->genes[0];
  
  for (i=0 ; i<num_genes ; i++) {
    if (genes2[i] == gA)
      return(i);
    if (genes2[i] == -gA)
      return(i + num_genes);
  }

  fprintf(stdout,"ERROR: calc_offset() could not locate gene %d\n", gA);
  
  return(-1);
}


/* same as in invdist.c */
int mgr_num_breakpoints(int *perm, int size) {
  int i, b;
  int pA, pA1;
  int pB, pB1;
  
  b = 0;
  if (perm[1] != 1)
    b++;
  for (i=2 ; i<size-1 ; i+=2) {
    pA  = perm[i];
    pA1 = (pA == size ? 1 : pA + 1);
    pB  = perm[i+1];
    pB1 = (pB == size ? 1 : pB + 1);
    if ((pB != pA1) && (pA != pB1))
      b++;
  }

  return (b);
}


/* return reversal distance
   and the graph and parameters about the graph, so that computations
   may continue
   new to grappa
*/


int invdist_noncircular_G(struct mgr_genome_struct *g1,
			  struct mgr_genome_struct *g2,
			  int offset, int num_genes,
			  mgr_distmem_t *distmem,
			  graph_t *G,
			  pathcounts_t *pcounts_G,
			  graphstats_t *graphstats)
{
  int b, c;
  int reversal_dist;

  hurdlecounts_t   hurdles_G;

#ifdef DEBUG
  fprintf(stdout,"invdist_noncircular: offset: %3d\n", offset);
#endif

  double_perms(g1, g2, offset,
	       num_genes,
	       distmem);
  buildG_grey(num_genes, 0 /* num_chromosomes */,
	      distmem, G);
 
  classify_cycle_path(G, pcounts_G);

  c = pcounts_G->num_cycles;
  form_connected_components(G);
  classify_connected_components(G);
  calc_num_hurdles_and_fortress(G, C_U, &hurdles_G);

  b = mgr_num_breakpoints(distmem->perm, G->size);

  /* reformulation of distance formula, page 215 */
  /* note: prior to page 215, c counts # cycles, including adjacencies.
     Starting middle of page 215, c counts # cycles, excluding adjacencies.
   */
  reversal_dist = b - c + hurdles_G.num_hurdles + hurdles_G.num_fortress;

#if 0  
  if (verbose) {
    fprintf(stdout,"Number of Breakpoints:        \t%d\n",b);
    fprintf(stdout,"Number of Cycles:             \t%d\n",c);
    fprintf(stdout,"Number of Hurdles:            \t%d\n",
	    hurdles_G.num_hurdles);
    fprintf(stdout,"Fortress:                     \t%d\n",
	    hurdles_G.num_fortress);
#ifdef DEBUG
    fprintf(stdout,"REVERSAL DISTANCE:            \t%3d\n",reversal_dist);
#endif
  }
#endif

  if (graphstats) {
    graphstats->br = b;
    graphstats->c4 = c;
    graphstats->h  = hurdles_G.num_hurdles;
    graphstats->f  = hurdles_G.num_fortress;
    graphstats->d  = reversal_dist;

    graphstats->num_components_u = hurdles_G.num_unoriented;
    graphstats->num_components_o = G->num_components - hurdles_G.num_unoriented;
  }
  
  return(reversal_dist);
}



/* unichromosomal reversal distance, linear
   Uses same parameter list as function of same name in GRAPPA

   _v versions: graphstats return structure is included
                (used to be "verbose option is included")
*/
int mgr_invdist_noncircular(struct mgr_genome_struct *g1, struct mgr_genome_struct *g2,
			int offset, int num_genes, mgr_distmem_t *distmem)
{
  graph_t G;
  pathcounts_t pcounts_G;

  /* return the distance, but discard the graph info */
  return invdist_noncircular_G(g1, g2, offset, num_genes, distmem,
			       &G, &pcounts_G, FALSE);
}

int invdist_noncircular_v(struct mgr_genome_struct *g1, struct mgr_genome_struct *g2,
			  int offset, int num_genes, mgr_distmem_t *distmem,
			  graphstats_t *graphstats)
{
  graph_t G;
  pathcounts_t pcounts_G;


  /* return the distance, but discard the graph info */
  return invdist_noncircular_G(g1, g2, offset, num_genes, distmem,
			       &G, &pcounts_G, graphstats);
}

/* unichromosomal reversal distance, linear
   Uses same parameter list as function of same name in GRAPPA

   _v versions: graphstats return structure included
                (used to be "verbose option is included")
*/
/* same as in invdist.c 
int invdist_circular(struct mgr_genome_struct *g1, struct mgr_genome_struct *g2,
		     int num_genes, mgr_distmem_t *distmem) {

  int offset;

  offset = calculate_offset(g1, g2, num_genes);
  
  return invdist_noncircular(g1, g2, offset, num_genes, distmem);
}
*/

// add the memory init/free function for distmem later
int mgr_invdist_noncircular_nomem(struct mgr_genome_struct *g1,
			      struct mgr_genome_struct *g2,
			      int offset, int num_genes) {
  mgr_distmem_t distmem;
  int dist;

  mcdist_allocmem(num_genes,
		  0, /* num_chromosomes */
		  &distmem);

  dist = invdist_noncircular(g1, g2, offset, num_genes, &distmem);

  mcdist_freemem(&distmem);
  
  return(dist);
}

int invdist_noncircular_nomem_v(struct mgr_genome_struct *g1,
				struct mgr_genome_struct *g2,
				int offset, int num_genes,
				graphstats_t *graphstats) {
  mgr_distmem_t distmem;
  int dist;

  mcdist_allocmem(num_genes,
		  0, /* num_chromosomes */
		  &distmem);

  dist = invdist_noncircular_v(g1, g2, offset, num_genes, &distmem,
			       graphstats);

  mcdist_freemem(&distmem);
  
  return(dist);
}

int mgr_invdist_circular_nomem(struct mgr_genome_struct *g1, struct mgr_genome_struct *g2,
			   int num_genes) {
  int offset;

  offset = mgr_calculate_offset(g1, g2, num_genes);
  
  return (mgr_invdist_noncircular_nomem(g1, g2, offset, num_genes));
}
