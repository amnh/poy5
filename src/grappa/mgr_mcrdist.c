/* mcrdist.c
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

/* Multichromosomal rearrangement distance algorithm
* Pevzner book, page 226,
* plus extensive changes described in Tesler JCSS article
* plus changes in Ozery-Flato and Shamir article
*/

/* contains some code fragments from invdist.c */

#include <stdio.h>
#include <stdlib.h>
#include "mgrstructs.h"
#include "mgr_uniinvdist.h"
#include "mgr_mcrdist.h"
#include "mgr_graph_edit.h"
#include "mgr_graph_components.h"
#include "mgr_e_malloc.h"



int get_caps_pgpath(graph_t *G, int capno);


/************************************************************************/

/* p. 226 step 1: Build the graph */


/* construct grey edges */
/* adapted from num_cycles */ 
void buildG_grey(int num_genes, int num_chromosomes, mgr_distmem_t *distmem,
		 graph_t *G) {
  int i, ind, j1, j2;
  int *greyEdges;
  int *invperm;
#if 0
  int *perm1;
#endif
  int *perm2, *perm;

  int size;
  int inchrom;       /* Boolean: TRUE  = inside chromosome
			         FALSE = at boundary
		      */
  int lowcap;        /* smallest doubled gene # corresponding to a cap */
  int cnum;          /* index 1,2,...,N of the current chromosome */
  int *chromNum;     /* chromNum[i]=chromosome # that vertex i is in */
  int *chromBd;      /* chromosome i=1,...,N (in doubled pi) occupies
			vertices chromBd[i],chromBd[i]+1,...,chromBd[i+1]-1
		      */

#ifdef DEBUG
  fprintf(outfile,"buildG_grey: num_genes: %3d  num_chromosomes: %3d\n",
	  num_genes, num_chromosomes);
#endif

  G->distmem = distmem;
  size = G->size = 2*num_genes+2;


  greyEdges  = distmem->greyEdges;
  invperm    = distmem->done;
#if 0
  perm1      = distmem->perm1;          /* doubled (gamma-hat)^(-1) */
#endif
  perm2      = distmem->perm2;          /* doubled pi-hat */
  perm       = distmem->perm;           /* composition */

  chromNum   = distmem->chromNum;
  chromBd    = distmem->chromBd;


#if 0
  /* this is always followed by a call to classify_cycle_path,
     which will set the following anyway: */
  G->components_good = FALSE;           /* info on components not computed */
#endif


#if 0
  /* This is now implicit */
  /* set black -edges */
  for (i=0 ; i<size ; i+=2) {
    blackEdges[i  ] = i+1;
    blackEdges[i+1] = i;
  }
#endif
  
  /* set grey -edges */
  for (i=0 ; i<size ; i++) {
    invperm[perm[i]] = i;
    greyEdges[i] = V_ADJ;
  }
  
  j1 = invperm[1];
  if (j1 != 1)
    greyEdges[0] = j1;
  
  for (i=1 ; i<size-1 ; i+=2) {
    ind = perm[i];
    if (ind < perm[i+1]) {
      j1 = invperm[ind-1];
      j2 = invperm[ind+2];
    }
    else {
      j1 = invperm[ind+1];
      j2 = invperm[ind-2];
    }
    if (j1 != i-1) {
      greyEdges[i] = j1;
    }
    if (j2 != i+2) {
      greyEdges[i+1] = j2;
    }
  }

  j1 = invperm[size-2];
  if (j1 != size-2)
    greyEdges[size-1] = j1;

  /* now delete some edges incident with caps or tails:
     mark
           pi-cap tail tail pi-cap
     regions between chromosomes

     lowest # for cap in unexpanded permutation:
        num_genes - 2*num_chromosomes + 1
     but under doubling it becomes
        2*num_genes - 4*num_chromosomes + 1;

     and mark the # of the chromosome each vertex is in
  */

  /* TODO: separate this into another routine */
  if (num_chromosomes == 0) {
    /* initialize chromosome numbers: it's one huge chromosome */
    for (i=0; i<size; i++)
      chromNum[i] = 0;

    chromBd[0] = 0;
    chromBd[1] = size;

  } else {
    /* multichromosomal */

    lowcap = 2*num_genes - 4*num_chromosomes + 1;

    cnum=0;  /* chromosome counter */

    i=0;
    inchrom = FALSE; /* at start of a chromosome */

    chromNum[i] = cnum;
    greyEdges[i++] = V_TAIL;

    while (i < size-1) {
      if (perm2[i] >= lowcap) {          /* regions between chromosomes */
	if (inchrom == FALSE) {
	  cnum++; /* started a new chromosome */
	  inchrom = TRUE;
	  chromBd[cnum] = i;

	  /* deleting tail is easy */
	  chromNum[i] = cnum; greyEdges[i++] = V_TAIL;

	  /* deleting grey edge leaving tail is longer */
	  chromNum[i] = cnum;
	  j1 = greyEdges[i];

	  /* JCSS - NULL CHROMOSOME BUG FIX */
	  greyEdges[i] = V_PICAP;

	  if (j1 == V_ADJ) {
	    /* adjacency, restore vertex number */
	    if (i%2 == 0) j1=i+1; else j1=i-1;
	  }
	  if (j1>=0) greyEdges[j1] = V_GTAIL;
	  i++;
	} else {
	  /* delete grey edge leaving tail */
	  j1 = greyEdges[i];

	  /* JCSS - NULL CHROMOSOME BUG FIX */
	  greyEdges[i] = V_PICAP;

	  if (j1 == V_ADJ) {
	    /* adjacency, restore vertex number */
	    if (i%2 == 0) j1=i+1; else j1=i-1;
	  }
	  if (j1>=0) greyEdges[j1] = V_GTAIL;

	  chromNum[i] = cnum;
	  i++;

	  /* delete tail */
	  chromNum[i] = cnum; greyEdges[i++] = V_TAIL;
	  inchrom = FALSE;
	}

      } else {                           /* region within a chromosome */
	chromNum[i] = cnum;
#if 0
	/* delete gamma-tails */
	j1 = greyEdges[i];
	if (j1>=0 && perm2[j1] >= lowcap) {
	  greyEdges[i] = V_GTAIL;
	  /* greyEdges[j1] = V_PICAP; */  /* already done */
	}
#endif
	i++;
      }
    }
    chromNum[i] = cnum+1;
    chromBd[cnum+1] = i;
    greyEdges[i++] = V_TAIL;  /* end */
  }


#ifdef DEBUG
#if 0
  for (i=0 ; i<size ; i++) {
    fprintf(outfile,"i: %3d Grey edge: (%3d, %3d)\n", perm[i],
	    perm[i], greyEdges[i] == V_ADJ ? V_ADJ : perm[greyEdges[i]]);
  }
#endif

  fprintf(outfile,"IP: ");
  for (i=0; i<size ; i++)
    fprintf(outfile,"%4d ",
	    invperm[i]);
  fprintf(outfile,"\n");


  fprintf(outfile,"Grey edges:");
#if 0
  fprintf(outfile,"\nP:  ");
  for (i=0; i<size ; i++)
    fprintf(outfile,"%4d ",perm[i]);
#endif

#if 0
  fprintf(outfile,"\nE:  ");
  for (i=0; i<size ; i++)
    fprintf(outfile,"%4d ",
	    greyEdges[i]<0 ? greyEdges[i] : perm[greyEdges[i]]);
#endif

  fprintf(outfile,"\nE0: ");
  for (i=0; i<size ; i++)
    fprintf(outfile,"%4d ",
	    greyEdges[i]);
  fprintf(outfile,"\n");

  fflush(outfile);
#endif

}



/*************************************************************************
 * Build the chromosome number/bdry arrays but not the edges.
 * This is a subset of step 1.  It is only used for special purposes,
 * not in normal graph construction.
 *************************************************************************/
void repairG_chrom(int num_genes, int num_chromosomes, mgr_distmem_t *distmem,
		   graph_t *G) {
  int i, j1;
  int *greyEdges;
#if 0
  int *perm1;
#endif
  int *perm2, *perm;

  int size;
  int inchrom;       /* Boolean: TRUE  = inside chromosome
			         FALSE = at boundary
		      */
  int lowcap;        /* smallest doubled gene # corresponding to a cap */
  int cnum;          /* index 1,2,...,N of the current chromosome */
  int *chromNum;     /* chromNum[i]=chromosome # that vertex i is in */
  int *chromBd;      /* chromosome i=1,...,N (in doubled pi) occupies
			vertices chromBd[i],chromBd[i]+1,...,chromBd[i+1]-1
		      */

#ifdef DEBUG
  fprintf(outfile,"buildG_grey: num_genes: %3d  num_chromosomes: %3d\n",
	  num_genes, num_chromosomes);
#endif

  G->distmem = distmem;
  size = G->size = 2*num_genes+2;


  greyEdges  = distmem->greyEdges;
#if 0
  perm1      = distmem->perm1;          /* doubled (gamma-hat)^(-1) */
#endif
  perm2      = distmem->perm2;          /* doubled pi-hat */
  perm       = distmem->perm;           /* composition */

  chromNum   = distmem->chromNum;
  chromBd    = distmem->chromBd;


  if (num_chromosomes == 0) {
    /* initialize chromosome numbers: it's one huge chromosome */
    for (i=0; i<size; i++)
      chromNum[i] = 0;

    chromBd[0] = 0;
    chromBd[1] = size;

  } else {
    /* multichromosomal */

    lowcap = 2*num_genes - 4*num_chromosomes + 1;

    cnum=0;  /* chromosome counter */

    i=0;
    inchrom = FALSE; /* at start of a chromosome */

    chromNum[i] = cnum;

    while (i < size-1) {
      if (perm2[i] >= lowcap) {          /* regions between chromosomes */
	if (inchrom == FALSE) {
	  cnum++; /* started a new chromosome */
	  inchrom = TRUE;
	  chromBd[cnum] = i;

	  /* deleting tail is easy */
	  chromNum[i] = cnum; i++;

	  /* deleting grey edge leaving tail is longer */
	  chromNum[i] = cnum;
	  j1 = greyEdges[i];

	  if (j1 == V_ADJ) {
	    /* adjacency, restore vertex number */
	    if (i%2 == 0) j1=i+1; else j1=i-1;
	  }
	  i++;
	} else {
	  /* delete grey edge leaving tail */
	  j1 = greyEdges[i];

	  if (j1 == V_ADJ) {
	    /* adjacency, restore vertex number */
	    if (i%2 == 0) j1=i+1; else j1=i-1;
	  }

	  chromNum[i] = cnum;
	  i++;

	  /* delete tail */
	  chromNum[i] = cnum; i++;
	  inchrom = FALSE;
	}

      } else {                           /* region within a chromosome */
	chromNum[i] = cnum;
	i++;
      }
    }
    chromNum[i] = cnum+1;
    chromBd[cnum+1] = i;
  }
}

/**************************************************************************
 * Step 1.5? deal with null chromosomes in g1 (gamma) that cause
 * pi-pi or gamma-gamma edges
 **************************************************************************/

#if 0
void check_null_edges(graph_t *G)
{
  int *greyEdges;
  int size;
  int n_pp, n_gg;
  int i;

  greyEdges  = G->distmem->greyEdges;
  size = G->size;

  n_pp = n_gg = 0;

  for (i=0 ; i<size ; i += 2)
    if (greyEdges[i] == greyEdges[i+1]) {
      if (greyEdges[i] == V_PICAP) n_pp++;
      else if (greyEdges[i] == V_GTAIL) n_gg++;
    }

  fprintf(outfile,"# pi-pi edges = %d, # gamma-gamma edges = %d\n",
	  n_pp,n_gg);
  
}
#endif


#if 0
/* p. 226, Step 2: Close all Pi-Gamma paths in simple components */
void build_Gbar(graph_t *G, pathcounts_t *pcounts)
{
  int i;

  mgr_component_t *components;
  int *cc;
  int *ptype;

  int size = G->size;
  mgr_distmem_t *distmem = G->distmem;

  ptype      = distmem->ptype;
  components = distmem->components;
  cc         = distmem->cc;

#ifdef DEBUG
  fprintf(outfile,"build_Gbar\n");
  fflush(outfile);
#endif


  for (i=0 ; i<size ; i++) {
    /* find Pi-Gamma paths */
    if (ptype[i] != P_PG  &&  ptype[i] != P_GP) continue;

    if (cc[i] < 0  ||  (components[cc[i]].hurdle & SIMPLE) != 0) {
      closepath(G, i, pcounts);        /* close this path */
    }
  }
}
#endif

/* p. 226, Step 2: Close all Pi-Gamma paths in simple components
 * return number of paths closed
 */
int build_Gbar(graph_t *G, pathcounts_t *pcounts, int num_chromosomes)
{
  int v, capno;
  int num_closed = 0;

  mgr_distmem_t *distmem = G->distmem;
  mgr_component_t *components = distmem->components;
  int *cc                 = distmem->cc;


#ifdef DEBUG
  fprintf(outfile,"build_Gbar\n");
  fflush(outfile);
#endif

  for (capno=2*num_chromosomes+1 ; capno>=2 ; capno--) {
    /* if the vertex for this cap position is in a pg-path,
     * move to the rightmost vertex of the pg-path
     */
    v = get_caps_pgpath(G, capno);
    if (v < 0) continue;

    if (cc[v] < 0  ||  (components[cc[v]].hurdle & SEMIKNOT) == 0) {
      closepath(G, v, pcounts);        /* close this path */
      num_closed++;
    }

  }

  return num_closed;
}



#if 0
/* p. 226, Step 3: Close all but one Pi-Gamma path in components
   having more than one such "inside" (within?) them */
void close_allbut1_pgpath(graph_t *G, pathcounts_t *pcounts)
{
  int i,j;

  mgr_component_t *components;
  int *cc;
  int *ptype;

  int size = G->size;
  mgr_distmem_t *distmem = G->distmem;

#ifdef DEBUG
  fprintf(outfile,"close_allbut1_pgpath: %d pgpaths\n",
	  pcounts->num_pg);

  fflush(outfile);
#endif

  ptype      = distmem->ptype;
  components = distmem->components;
  cc         = distmem->cc;

  /* set all components "pending" until their first Pi-Gamma path
     encountered */
  for (j=0 ; j < G->num_components ; j++)
    components[j].scanstate = SCAN_PENDING;


  /* find Pi-Gamma paths */
  for (i=0 ; i<size ; i++) {
    /* skip it if it's not the left vertex of a pi-gamma path */
    if (ptype[i] != P_PG  &&  ptype[i] != P_GP) continue;

    /* If it's the first one encountered in this component,
       it stays open.
       If it's the 2nd or later, close it.
       (Which one to keep open is arbitrary; this is just the easiest choice.)
    */

    if (components[cc[i]].scanstate == SCAN_PENDING) {
      components[cc[i]].scanstate = SCAN_ACTIVE;
    } else {
      closepath(G, i, pcounts);        /* close this path */
    }
  }
}
#endif

/* p. 226, Step 3: Close all but one Pi-Gamma path in components
 * having more than one such within them.
 * Since we closed all paths in simple components, the only
 * components with >1 pg-path are semi-knots with 2 pg-paths, both
 * on the same chromosome.
 * Seek these out and close the pg-path terminating on the righthand cap.
 */
void close_allbut1_pgpath(graph_t *G, pathcounts_t *pcounts,
			  int num_chromosomes)
{
  int v1, v2;
  int chr;
  int compno1;

  mgr_distmem_t *distmem = G->distmem;
  mgr_component_t *components = distmem->components;
  int *cc                 = distmem->cc;


#ifdef DEBUG
  fprintf(outfile,"close_allbut1_pgpath: %d pgpaths\n",
	  pcounts->num_pg);

  fflush(outfile);
#endif

  for (chr=num_chromosomes ; chr>=1 ; chr--) {
    /* pg-path index, if any, for cap at right of chromosome # chr */
    v2 = get_caps_pgpath(G, 2*chr+1);
    if (v2<0) continue;

    /* pg-path index, if any, for cap at left of chromosome # chr */
    v1 = get_caps_pgpath(G, 2*chr);
    if (v1<0) continue;

    /* are they in the same component, and is it a semiknot? */
    compno1 = cc[v1];
    if (compno1 >=0 && compno1 == cc[v2] &&
	(components[compno1].hurdle & SEMIKNOT)) {
      closepath(G, v2, pcounts);       /* close righthand path in semiknot */
    }
  }
}


/****************************************************************************/

int is_interchrom_or_oriented(graph_t *G, int v1, int v2) {
  if ((v1-v2)%2 == 0) return TRUE;   /* oriented */

  return (v1-v2)%2 == 0      /* oriented */
    || G->distmem->chromNum[v1] != G->distmem->chromNum[v2];  /* interchrom */
}



/* p. 226 step 5: join pi-pi and gamma-gamma paths
 * and JCSS - NULL CHROMOSOME BUG FIX: close remaining pi-pi paths
 */

void join_pp_gg_paths(graph_t *G, pathcounts_t *pcounts)
{
  int pp_index=0, gg_index=0;   /* scan for pi-pi and gamma-gamma paths */
  int v1, v2;                   /* vertices to join */
  int v3, v4;
  int *Bcap, *Gcap, *ptype;

  mgr_distmem_t *distmem = G->distmem;
  Bcap       = distmem->Bcap;
  Gcap       = distmem->Gcap;
  ptype      = distmem->ptype;

  v1=v2=v3=v4=0;


#ifdef DEBUG
  fprintf(outfile,"join_pp_gg_paths: num_pp=%d\n",pcounts->num_pp);
  fflush(outfile);
#endif

  while (pcounts->num_gg > 0) {
    while (ptype[pp_index] != P_PP) pp_index++;  /* next pi-pi path */
    while (ptype[gg_index] != P_GG) gg_index++;  /* next gamma-gamma path */

    /* 4 possible edges to test.
     * Need a pair of edges where at least one edge is
     * oriented or interchromosomal.
     * Either use (v1,v2) & (v3,v4) if @ least one is oriented/intrachrom,
     * or else use (v1,v3) & (v2,v4).
     */
    v1 = Bcap[pp_index]; v2 = Bcap[gg_index];
    v3 = Gcap[pp_index]; v4 = Gcap[gg_index];

    if (!is_interchrom_or_oriented(G,v1,v2) &&
	!is_interchrom_or_oriented(G,v3,v4)) {
      v2 = v4; v4 = Bcap[gg_index];
#if DEBUG
      if (!is_interchrom_or_oriented(G,v1,v2) &&
	  !is_interchrom_or_oriented(G,v3,v4)) {
	fprintf(outfile,"ERROR: no interchromosomal or oriented edge between paths %d, %d\n",gg_index,pp_index);
	fflush(outfile);
      }
#endif
    }

    /* join the paths into a cycle */
    join_paths_2_cycle(G, v1,v2, v3,v4, pcounts);

    pp_index++; gg_index++;  /* skip past these now that they're joined */
  }


  /* JCSS - NULL CHROMOSOME BUG FIX
   * close any remaining pi-pi paths */
  while (pcounts->num_pp > 0) {
    while (ptype[pp_index] != P_PP) pp_index++;  /* next pi-pi path */
    closepath(G, pp_index, pcounts);
  }
}



/****************************************************************************/
/*                       NEW version of steps 7-18                          */
/****************************************************************************/


/****************************************************************************
 * improved version of join_close_semiknots:
 * attempt to join pi-gamma paths on same chromosome when possible
 ****************************************************************************/

/**** BEGIN CHANGE ****/

/*
 * Get the index of the path or cycle associated with cap # capno
 * If it is a pi-gamma path, return it, else a code -1
 *
 * On chromosome i=1,2,...,N,
 * the left cap is indexed capno=2i, the right cap is indexed capno=2i+1
 *
 */
int get_caps_pgpath(graph_t *G, int capno)
{
  mgr_distmem_t *distmem = G->distmem;
  int *chromBd       = distmem->chromBd;
  int *cycle         = distmem->labeled;
  int *ptype         = distmem->ptype;

  int v;


  /* set v to the index of the vertex at the left or right end of
   * the chromosome, beyond the tail
   */
  if (capno % 2 == 0) {
    /* left end of chromosome (capno/2)+1 */
    v = chromBd[capno/2] + 2;
  } else {
    /* right end of chromosome ((capno-1)/2)+1 */
    v = chromBd[capno/2 + 1] - 2;
  }

  /* Get the index of the cycle/path containing v */
  v = cycle[v];

  if (v >= 0 && (ptype[v] == P_PG || ptype[v] == P_GP)) {
    return v;
  } else {
    return -1;
  }
}

/* the only paths that are left are exactly one pi-gamma path per semiknot
 * hence, # semiknots = # pi-gamma paths
 * Join them in pairs as appropriate, preferring ones on the same chromosome,
 * and then close the remaining ones, in the manner dictated by the algorithm.
 *
 * s2flag: flag concerning treatment of 2 semiknots
 */

void join_close_semiknots(graph_t *G, pathcounts_t *pcounts,
			  int s2flag
			  )
{
  /* unmatched pi-gamma paths from anywhere */
  int index_pg_u1 = -1;       /* index of unmatched pi-gamma path
			       * -1 means none available right now
			       */
  int index_pg_u2 = -1;


  
  int index_pg1, index_pg2; /* indices of pi-gamma paths on current chrom */

  int join_till;           /* when to stop joining paths */

  int capno;


  /* quick exit if no semiknots */
  if (pcounts->num_pg == 0) return;



  /* While there are more than 2 semiknots, we may pair them up into cycles
   * arbitrarily.
   *
   * For aesthetic reasons (drawing breakpoint graphs),
   * if any chromosome has 2 pi-gamma paths, we join them both together.
   * The remaining pi-gamma paths are paired up consecutively left to right.
   *
   * The last remaining 0, 1, or 2 pi-gamma paths are closed instead of
   * joined, depending on the parameter gr and the parity of the number
   * of pi-gamma paths.
   */

  if (pcounts->num_pg % 2)
    join_till = 1;          /* odd # semi-knots, pair up till end */
  else if (s2flag)
    join_till = 2;          /* even # & has gr. real knot, don't join last 2 */
  else
    join_till = 0;          /* even #, no gr. r.k., pair them all up */


  /* Loop over all chromosomes.
   * If a pair of pi-gamma paths is found on the same chromosome,
   * join them together.
   * If stray paths are found, keep up to two of them and then join them.
   */

  /* cap number
   * chromosome i=1,2,...,N has cap numbers 2i, 2i+1
   */
  capno = 2;

  while (pcounts->num_pg > join_till) {
    /* determine if the cap at each end of current chromosome is in a
     * pi-gamma path
     */

    index_pg1 = get_caps_pgpath(G,capno++);
    index_pg2 = get_caps_pgpath(G,capno++);
    if (index_pg1 >= 0 && index_pg2 >= 0) {
      /* make sure they weren't accidentally the same one */
      if (index_pg1 == index_pg2) {
	index_pg_u2 = index_pg1;
      } else {
	/* bind the paths together into a cycle */
	join_pgpaths_2_cycle(G, index_pg1, index_pg2, pcounts);
      }
    } else if (index_pg1 >= 0) {
      index_pg_u2 = index_pg1;
    } else if (index_pg2 >= 0) {
      index_pg_u2 = index_pg2;
    }

    if (index_pg_u2 >= 0) {
      if (index_pg_u1 >= 0) {
	if (index_pg_u1 != index_pg_u2) {
	  /* bind the paths together into a cycle */
	  join_pgpaths_2_cycle(G, index_pg_u1, index_pg_u2, pcounts);
	  index_pg_u1 = index_pg_u2 = -1;
	} else {
	  /* same as previous path, discard */
	  index_pg_u2 = -1;
	}
      } else {
	/* don't yet have a path to match it with */
	index_pg_u1 = index_pg_u2; index_pg_u2 = -1;
      }
    }
  }

  /* either 2 semiknots left but no greatest real knot,
     OR 1 semiknot left,
     OR no semiknots left.
     In all cases, close all remaining semiknots */
  if (pcounts->num_pg > 0) {
    /* close the stored path, if any */
    if (index_pg_u1 >= 0) {
      closepath(G, index_pg_u1, pcounts);
    }

    /* find the remaining ones and close them */
    while (pcounts->num_pg > 0) {
      index_pg1 = get_caps_pgpath(G,capno++);
      if (index_pg1 >= 0) {
	closepath(G, index_pg1, pcounts);
      }
    }
  }
}




/**** END CHANGE ****/





/****************************************************************************/

/* p. 226, step 19 */
/* As described in the book, it's incomplete.
 * This is replacement for it.
 */


/****************************************************************************
 *  Correct algorithm for forming concatenations, using notions of
 *  "proper flippings" and "proper bondings"
 *
 * find_proper_flipping: flip chromosomes of pi so that all
 * interchromosomal components are oriented
 * This version flips the leftmost chromosome of an interchrom unor comp
 * (and if more than one such comp has same leftmost chrom, just flip it
 * once)
 *
 * find_optimal_bonding: rearranges chromosomes of pi so that they
 * induce optimal concats of pi and gamma that are properly flipped, and
 * have 0 or 1 improper bonds.
 *
 * compute_concatenates: whether pi is optimally bonded or not, read off
 * the concatenates of pi and gamma it induces.  For gamma, there can be
 * multiple segments, and they are put in some order.
 ****************************************************************************/


/* find proper flipping
   This version chooses the left-most chromosome incident on
   an unoriented interchromosomal component, as the one to flip.
   */
void find_proper_flipping(graph_t *G, int num_chromosomes)
{
  int j;
  int ind, chrom;   /* vertex and chromosome */



  mgr_distmem_t *distmem = G->distmem;

  /*int num_components = G->num_components;*/ /* changes during this routine */
  mgr_component_t *components = distmem->components;
  int         *chromNum   = distmem->chromNum;

  /* stack not in use at this point, be careful if this changes */
  /* it's 1-based, with indices 1,2,...,num_chromosomes */
  int         *flipped    = distmem->stack;

  pathcounts_t pcounts_G;



#ifdef DEBUGf
  fprintf(outfile,"find_proper_flipping\n");
  fflush(outfile);
#endif


  form_connected_components(G);
  classify_connected_components(G);

  /* each improper chromosome should be flipped just once */
  for (j=1 ; j <= num_chromosomes; j++)
    flipped[j] = FALSE;


  /* check each component */

  for (j=0 ; j < G->num_components ; j++) {
#ifdef DEBUGf
    fprintf(outfile,"check flipping: Component %d: index=%d, oriented=%d, intrachrom=%d, real=%d, hurdle=oct'%o\n",
	    j,components[j].index,components[j].oriented,
	    components[j].intrachrom,components[j].real,
	    components[j].hurdle);
#endif
    if (components[j].oriented == FALSE &&
	components[j].intrachrom == FALSE) {
      /* leftmost chromosome with this component */
      ind = components[j].index;
      chrom = chromNum[ind];

      /* if already flipped, fine, else flip it */
      if (!flipped[chrom]) {
	flipped[chrom] = TRUE;
#ifdef DEBUGf
	fprintf(outfile,"Flipping chromosome %d\n", chrom);
#endif
	flip_chromosomes_nocycles(G,chrom,chrom);
      }
    }
  }


  /* If any chromosomes were flipped, the cycle & path info needs to
   * be recomputed.
   * If components need recomputation, it will be done later.
   */

  if (!G->components_good) {
    classify_cycle_path(G,&pcounts_G);
  }

#ifdef DEBUGf
  fprintf(outfile,"Properly flipped.\n");
  fflush(outfile);
#endif
}


/* compute the concatenates of pi and gamma induced by the graph */
/* input: dop2=TRUE to do cappedp1, cappedp2, and add in tail edges
               FALSE to do only cappedp1 */
/* output:
 * when dop2=FALSE: only do cappedp1
 * when dop2=TRUE:
 *    G->distmem->cappedp1 & cappedp2 = cappings in orig alphabet
 *    grayEdges updated to include the tail cycles
 */
void compute_concatenates(graph_t *G, int num_genes, int num_chromosomes,
			  int dop2)
{
  int *greyEdges, *perm2, *cappedp1, *cappedp2, *chromBd;
  int *capMark;
  int num_read, next, i1, i2;
  int j;

  int lowcap1;
  int next_freecap, capno;
  int atend;



  mgr_distmem_t *distmem = G->distmem;


  greyEdges  = distmem->greyEdges;
  chromBd    = distmem->chromBd;

  perm2      = distmem->perm2;          /* doubled pi-hat */
  cappedp1   = distmem->cappedp1;       /* new capping of gamma-hat */
  cappedp2   = distmem->cappedp2;       /* new capping of pi-hat */

  capMark    = distmem->capMark;  /* which caps have been used? */


#ifdef DEBUGb
  fprintf(outfile,"compute_concatenates\n");
  /*print_arrays(G);*/
#endif



  /* show haven't read any chromosomes yet */
  lowcap1 = num_genes - 2*num_chromosomes + 1;     /* with single #s */
  for (j=2*num_chromosomes-1; j>=0 ; j--) {
    capMark[j] = 0;   /* cap lowcap+j hasn't been used yet */
  }



  num_read = 0;                         /* how many genes have been read? */

  /* chromosome i=1,2,... has lcap at position 2i-2, rcap at 2i-1 */
  /* the caps in positions 0,1,2,..., next_freecap-1  are known to be used */
  next_freecap = 1;                     

  /* alternately follow    white gap + grey edge + w + g + w + g ... + w */
  next = 1;
  while (num_read < num_genes) {
#if 0
#ifdef DEBUGb
    fprintf(outfile,"nr=%d nx=%d\n",num_read,next);
#endif
#endif

    if (greyEdges[next] == V_TAIL) {
      /* entering the chromosome pointed to by next */
      /* have we either run off the right end, or entered a chrom that's
	 already been used?
	 Either means a segment of gamma has been completed and we
	 need to find a new start. */


      i1 = perm2[next];
      capno = (i1+1)/2 - lowcap1;
      atend = (next == G->size-1) || (capno>=0 && capMark[capno]);

      /* if have reached an end, pick any unused cap as a new start */
      while (atend) {
	next = chromBd[(next_freecap+1)/2 + 1];
	if (next_freecap % 2 == 1)
	  next--;             /* right end of chromosome =
				 one before left end of next chromosome */

	next_freecap++;

	i1 = perm2[next];
	capno = (i1+1)/2 - lowcap1;   /* normalized cap# 0,1,...,2*nchromo-1 */
	atend = capMark[capno];       /* has it been used yet? */

#ifdef DEBUGb
	fprintf(outfile,
		"next_freecap=%d  next=%d  i1=%d  capno=%d  atend=%d\n",
		next_freecap-1, next, i1, capno, atend);
	fflush(outfile);
#endif

      }


      /* mark current cap as used */
      capMark[capno] = 1;
#ifdef DEBUGb
      fprintf(outfile,"Marking capno=%d on entry\n", capno);
      fflush(outfile);
#endif
    }



    /* already got this from above */
    i1 = perm2[next];         /* doubled perm, one end of entry */

    if (next % 2 == 1) next++; else next--;     /* follow white gap */

    i2 = perm2[next];                 /* doubled perm, other end of entry */

    /* to copy in doubled form: append i1, then i2 */
    /* signed form: */
    if (i1<i2)
      cappedp1[num_read] = i2/2;      /* positive entry */
    else
      cappedp1[num_read] = -i1/2;     /* negative entry */
    num_read++;

    /* follow grey edge */
    i1 = greyEdges[next];

    /* if it's a tail, mark that we have read it */
    if (i1 == V_TAIL) {
      capno = (i2+1)/2 - lowcap1;
      capMark[capno] = 1;
#ifdef DEBUGb
      fprintf(outfile,"Marking capno=%d on exit\n", capno);
      fprintf(outfile,"Following tail %d at next=%d\n",(i2+1)/2,next);
      fflush(outfile);
#endif
    }

    if (i1 == V_ADJ || i1 == V_TAIL) {
      /* adjacencies were ignored before, but are needed now */
      if (next % 2 == 0) i1=next+1; else i1=next-1;
    }
    next = i1;
  }

  /* output new capping of pi */
  if (dop2)
    for (j=0 ; j<num_genes ; j++) {
      i1 = perm2[2*j+1]; i2 = perm2[2*j+2];
      if (i1<i2)
	cappedp2[j] = i2/2;
      else
	cappedp2[j] = -i1/2;
    }

#if 0
#ifdef DEBUG
  fprintf(outfile,"\npi-cap:\n");
  for (j=0; j < num_genes ; j++) {
    fprintf(outfile,"%4d ", cappedp2[j]);
  }
  fprintf(outfile,"\n");

  fprintf(outfile,"\ngamma-cap:\n");
  for (j=0; j < num_genes ; j++)
    fprintf(outfile,"%4d ",
	    cappedp1[j]);
  fprintf(outfile,"\n");

  fflush(outfile);
#endif // DEBUG
#endif // 0
}


/* Form an array of the caps of chromosomes of Gamma */
/* Output:
   array i=lowcap,lowcap+1,...,hicap

   G->capMark[i-lowcap] = j
   G->capMark[j-lowcap] = i

   where the chromo of gamma beginning with cap +/- i  ends with cap +/- j
*/
/* TODO: redo with more efficient algorithm that uses
   perm1 (the inverse of permutation 1) array to locate cap locations in
   gamma;
   and uses chromBd to locate cap locations in pi;
   and compare these appropriately.
   */
void compute_gammacaps(graph_t *G, int num_genes, int num_chromosomes)
{
  mgr_distmem_t *distmem = G->distmem;
  int *cappedp1, *capMark;
  int lowcap;
  int i,g;
  int lcap,rcap;

  cappedp1   = distmem->cappedp1;       /* new capping of gamma-hat */
  capMark    = distmem->capMark;        /* array of gamma chromo cap pairs */

  /* lowest cap number in nondoubled system */
  lowcap = num_genes - 2*num_chromosomes + 1;


  /* overkill: get gamma as ordinary perm */
  compute_concatenates(G, num_genes, num_chromosomes, FALSE);

  /* scan gamma left to right for caps */
  lcap = -1;
  for (i=0; i<num_genes; i++) {
    g = cappedp1[i];
    if (g<0) g=-g;
    if (g>=lowcap) {
      if (lcap == -1) {
	/* entering chromosome */
	lcap = g;
      } else {
	/* exiting chromosome */
	rcap = g;
	capMark[lcap-lowcap] = rcap;
	capMark[rcap-lowcap] = lcap;

	lcap = -1;
      }
    }
  }

#ifdef DEBUGb
  fprintf(outfile,"capMark: ");
  for (i=0; i<2*num_chromosomes ; i++)
    fprintf(outfile,"%4d ",
	    capMark[i]);
  fprintf(outfile,"\n");
#endif
}



/* determine a proper flipping with at most one improper bonding */
/* return number of improper bonds */
int find_optimal_bonding(graph_t *G, int num_genes, int num_chromosomes)
{
  int c;
  int *chromBd, *perm2, *capMark;
  int i1, next, cap1, cap2, lowcap;
  int num_badbonds = 0;

  mgr_distmem_t *distmem = G->distmem;

  chromBd    = distmem->chromBd;

  perm2      = distmem->perm2;          /* doubled pi-hat */

  capMark    = distmem->capMark;  /* which caps have been used? */



  /* lowest cap number in nondoubled system */
  lowcap = num_genes - 2*num_chromosomes + 1;


#ifdef DEBUGb
  fprintf(outfile,"find_optimal_bonding\n");
  /*print_arrays(G);*/
#endif




  /* initialize array G->capMark so that

   G->capMark[i-lowcap] = j
   G->capMark[j-lowcap] = i

   where the chromo of gamma beginning with cap +/- i  ends with cap +/- j.

   Then as we bond chromos together, the internal caps are not updated,
   but the terminal ones should satisfy the same relations as above,
   for a string of chromos beginning with +/- i and ending with +/- j.
   */

     

  /* initialize the array */
  compute_gammacaps(G, num_genes, num_chromosomes);

  /* form a proper left flipping of G */
  find_proper_flipping(G, num_chromosomes);

  /* scan chromos right to left, bonding when allowed,
     rearranging when not allowed to bond */
  for (c=num_chromosomes; c>=2; c--) {
    /* bond chromo c to chromo c-1, if permitted */
    next = chromBd[c];
    i1 = perm2[next];
    cap1 = (i1+1)/2;
    cap1 = (cap1>0) ? cap1 : -cap1;

    i1 = perm2[next-1];
    cap2 = (i1+1)/2;
    cap2 = (cap2>0) ? cap2 : -cap2;

#ifdef DEBUGb
    fprintf(outfile,"Examining bond %d:%d between chromos %d and %d\n",
	    cap1,cap2, c-1,c);
#endif

    if (capMark[cap1-lowcap] == cap2) {
      /* bond not permitted! */

#ifdef DEBUGb
      fprintf(outfile,"Bad bond, c=%d, cap2=%d! Avoiding\n", c, cap2);
#endif

#define KEEP_CHROM_IN_ORDER
      if (c>2) {
#ifdef KEEP_CHROM_IN_ORDER
	/* Attempt to just flip chromo c-1 so that chromos stay in order.
	 * This is NOT guaranteed to work; it may be that all intrachrom
	 * components incident on chromo c-1 have leftmost vertex in c-1.
	 * But then, the find_proper_flipping below will just flip c-1 back.
	 * If that happens, we swap two chromos, which _is_ guaranteed to work.
	 */
	flip_chromosomes(G, c-1, c-1, TRUE);
	find_proper_flipping(G, num_chromosomes);
	if (perm2[next-1] == i1) {
	  /* proper flipping put this chromo back as it was; we failed! */
#endif // KEEP_CHROM_IN_ORDER
#ifdef DEBUGb
	  fprintf(outfile,"   flipping 1 chromo unsuccessful, trying flip 2\n");
#endif // DEBUGb
	  /* Swap chromos c-1 and c-2 */
	  flip_chromosomes(G, c-2, c-1, TRUE);

	  /* proper left flipping keeps chromos c,c+1,... intact */
	  find_proper_flipping(G, num_chromosomes);
#ifdef KEEP_CHROM_IN_ORDER
	}
#ifdef DEBUGb
	else
	  fprintf(outfile,"   flipping 1 chromo successful!\n");
#endif // DEBUGb
#endif // KEEP_CHROM_IN_ORDER

	/* now bond is guaranteed */
      } else {
	/* bond between chromo 2 and chromo 1;
	   may not be possible to get proper bond and proper flipping.
	   If not possible, proper flipping wins.
	   */
	flip_chromosomes(G, c-1, c-1, TRUE);
	find_proper_flipping(G, num_chromosomes);
      }

      /* update cap2 */
      i1 = perm2[next-1];
      cap2 = (i1+1)/2;
      cap2 = (cap2>0) ? cap2 : -cap2;

      /* should only happen when c=2, and even then only
       * on rare occasions
       */
      if (capMark[cap1-lowcap] == cap2)
	num_badbonds++;
    }

    /* record new ends of block chromo due to this bond */
    /* dubious when c=2 and could not get proper bond,
       but we abort right after this w/o further using the info */

    cap2 = capMark[cap2-lowcap];
    cap1 = capMark[cap1-lowcap];

    capMark[cap1-lowcap] = cap2;
    capMark[cap2-lowcap] = cap1;
  }

  return num_badbonds;
}


/****************************************************************************
 * Count number of multichromosomal breakpoints
 * Breaks down into internal vs external breakpoints (@ caps)
 ****************************************************************************/
void num_breakpoints_mc(int *perm,
			int size,
			int lowcap1,  /* lowest cap in non-doubled perms */
			struct mgr_genome_struct *g1,
			struct mgr_genome_struct *g2,

			/* return values: */
			int *bp_int,   /* # internal breakpoints */
			int *bp_ext    /* # external breakpoints */
			)
{
  int i, b_i, b_e;
  int pA, pA1;
  int pB, pB1;

  int *g1genes = g1->genes;
  int *g2genes = g2->genes;

  int g2a, g2b, cap2a, cap2b;
  int pA2, pB2;

  int k;

  b_i = b_e = 0;
  if (perm[1] != 1)
    {
      b_e++;
    }
  for (i=2 ; i<size-1 ; i+=2) {
    pA  = perm[i];
    pA1 = (pA == size ? 1 : pA + 1);
    
    pB  = perm[i+1];
    pB1 = (pB == size ? 1 : pB + 1);

    if ((pB != pA1) && (pA != pB1)) {
      /* check if any caps are involved */

      g2a = g2genes[(i-1)/2];     cap2a = (g2a >= lowcap1 || g2a <= -lowcap1);
      g2b = g2genes[i/2];         cap2b = (g2b >= lowcap1 || g2b <= -lowcap1);

      if (cap2a || cap2b) {
	/* caps are involved */
	if (cap2a) {
	  if (cap2b) {
	    /* form (cap1,cap2).
	     * If tail, don't count.
	     * Else, null chromosome;
	     *   should count as long as g1 doesn't have matching null chrom.
	     *   TODO: fix this.  For now, assume only g1 or g2 has null chr.
	     */
	    if ((g2a > 0 && (g2a-lowcap1) % 2 == 0) ||
		(g2a < 0 && ((-g2a)-lowcap1) % 2 == 1)) {
	      b_e++;
	    }
	  } else { /* cap2a && !cap2b */
	    /* form (C,x) */
	    pA2 = ((pB % 2) == 0) ? pB1 : pB-1;

	    k = g1genes[(pA2-1)/2];
	    if (k < lowcap1 && k > -lowcap1) {
	      /* neighbor of 2nd # is not a cap in g1 */
	      b_e++;
	    }
	  }
	} else { /* !cap2a && cap2b */
	  /* form (x,C) */
	  pB2 = ((pA % 2) == 0) ? pA1 : pA-1;
	  k = g1genes[(pB2-1)/2];
	  if (k < lowcap1 && k > -lowcap1) {
	    /* neighbor of 1st # is not a cap in g1 */
	    b_e++;
	  }
	}
      } else {
	/* no caps, internal b.p. */
	b_i++;
      }
    }
  }

  *bp_int = b_i;
  *bp_ext = b_e;

  return;
}



/****************************************************************************/
/****************************************************************************/

/* copy signed permutation into genome structure */
void copy_perm_to_genome(int *perm,
			 struct mgr_genome_struct *g,
			 int num_genes)
{
  int i;

  for (i=0; i<num_genes; i++)
    g->genes[i] = perm[i];

}




/* p. 226 step 20-21: in scenario.c */

/**********************************************************************/




/* if allocation errors occur, it prints messages but
   doesn't exit cleanly */
void mcdist_allocmem(int num_genes, int num_chromosomes,
		     mgr_distmem_t *distmem)
{
  int gsize = 2*num_genes+2;      /* # entries in doubled genome */

  distmem->perm1 = (int *) e_malloc(gsize*sizeof(int), "perm1");

  distmem->perm2 = (int *) e_malloc(gsize*sizeof(int), "perm2");

  distmem->perm = (int *) e_malloc(gsize*sizeof(int), "perm");

  distmem->done = (int *) e_malloc(gsize*sizeof(int), "done");

  distmem->greyEdges = (int *) e_malloc(gsize*sizeof(int), "greyEdges");

  distmem->stack = (int *) e_malloc(gsize*sizeof(int), "stack");

  distmem->oriented = (int *) e_malloc(gsize*sizeof(int), "oriented");

  distmem->cc = (int *) e_malloc(gsize*sizeof(int), "cc");

  distmem->labeled = (int *) e_malloc(gsize*sizeof(int), "labeled");


  /* # components is at most # genes+1, since every pair of verts is
     joined by a black edge */
  distmem->components = (mgr_component_t *)
     e_malloc((num_genes+1)*sizeof(mgr_component_t), "components");

  distmem->Bcap = (int *) e_malloc(gsize*sizeof(int), "Bcap");

  distmem->Gcap = (int *) e_malloc(gsize*sizeof(int), "Gcap");

  distmem->ptype = (int *) e_malloc(gsize*sizeof(int), "ptype");

  distmem->chromNum = (int *) e_malloc(gsize*sizeof(int), "chromNum");

  distmem->chromBd = (int *) e_malloc((num_chromosomes+2)*sizeof(int),
				      "chromBd");
  distmem->capMark = (int *) e_malloc((2*num_chromosomes)*sizeof(int),
				      "capMark");

  distmem->cappedp1 = (int *) e_malloc(num_genes * sizeof(int), "cappedp1");
  
  distmem->cappedp2 = (int *) e_malloc(num_genes * sizeof(int), "cappedp2");
}

void mcdist_freemem(mgr_distmem_t *distmem) {

  free(distmem->cappedp2);
  free(distmem->cappedp1);
  free(distmem->capMark);
  free(distmem->chromBd);
  free(distmem->chromNum);
  free(distmem->ptype);
  free(distmem->Gcap);
  free(distmem->Bcap);

  free(distmem->components);
  free(distmem->labeled);
  free(distmem->cc);
  free(distmem->oriented);
  free(distmem->stack);

  free(distmem->greyEdges);
  free(distmem->done);
  free(distmem->perm1);
  free(distmem->perm2);
  free(distmem->perm);
}


int mcdist_noncircular_nomem(struct mgr_genome_struct *g1,
			     struct mgr_genome_struct *g2,
			     int num_genes, int num_chromosomes,
			     graphstats_t *graphstats)
{
  mgr_distmem_t distmem;
  int dist;


  /* allocate memory */
  mcdist_allocmem(num_genes, num_chromosomes,
		  &distmem);


  /* compute distance */
  dist = mcdist_noncircular(g1, g2,
			    num_genes, num_chromosomes,
			    &distmem,
			    graphstats);

  /* free memory */
  mcdist_freemem(&distmem);

  return(dist);
}



/* Formula on page 224 is computed w/info from Algorithm
   on p. 226, steps 1 and 2

   This "full" version returns info on the graph
   so that we can continue with the remainder of the algorithm.
   Returns G, pcounts_G, s2flag

   if graphstats != NULL: fill graphstats with values of the parameters

   (was: verbose=TRUE prints out info)
*/
int mcdist_noncircular_full(struct mgr_genome_struct *g1,
			    struct mgr_genome_struct *g2,
			    int num_genes, int num_chromosomes,
			    mgr_distmem_t *distmem,

			    graph_t *G,
			    pathcounts_t *pcounts_G,
			    int *s2flag,   /* treatment of 2-semiknots */
			    graphstats_t *graphstats)
{
  /* hurdlecounts_t   hurdles_G; */
  /* hurdlecounts_t   hurdlesIU_G; */
  hurdlecounts_t   hurdlesRU_G;

  /* the parameters in the distance formula, p. 224 */
#if 0
  int param_b, param_c, param_p, param_rr, param_s, param_gr, param_fr;
#endif
  int param_b, param_c, param_p, param_r, param_s, param_gr, param_fr;
  int has_sgrk;
#if 0
  int num_simple;
#endif
  int dist;



  /* Algorithm on page 226, Step 1 */

  /* build the doubled permutations */
  double_perms(g1, g2,
	       0, /* offset */
	       num_genes,
	       distmem);



  /* construct the grey edges, caps, tails (p. 225(d)) */
  buildG_grey(num_genes, num_chromosomes, distmem, G);

  if (graphstats) {
    /* calc # breakpoints before further graph editing performed */
    num_breakpoints_mc(distmem->perm,
		       G->size,
		       num_genes - 2*num_chromosomes + 1,
		       g1, g2,

		       &graphstats->bp_int,
		       &graphstats->bp_ext);
  }


  /* # black edges in the graph */
  param_b = num_genes - num_chromosomes;
#ifdef DEBUG
  fprintf(outfile,"# black edges = %d\n",param_b);
#endif



  classify_cycle_path(G, pcounts_G);

  /* total number of paths and cycles of all kinds */
  param_c =
    pcounts_G->num_adjacencies + pcounts_G->num_cycles
    + pcounts_G->num_pp + pcounts_G->num_gg + pcounts_G->num_pg;
  
  /* number of gamma-gamma paths in G */
  param_p = pcounts_G->num_gg;


  form_connected_components(G);
  classify_connected_components(G);


  /* find hurdles of all kinds */

  /* ordinary hurdles aren't relevant to this algorithm
     The test code worked anyway */
  /*
  calc_num_hurdles_and_fortress(G, C_U, &hurdles_G);
  */

#if 0
  calc_num_hurdles_and_fortress(G, C_IU, &hurdlesIU_G);
  calc_num_hurdles_and_fortress(G, C_RU, &hurdlesRU_G);
  param_s = calc_num_semiknots(G);
#endif

  calc_num_hurdles_and_fortress(G, C_RU, &hurdlesRU_G);
  param_r = hurdlesRU_G.num_hurdles;                /* # real knots in G */
  param_s = calc_num_semirealknots(G,&has_sgrk);    /* # semi-real-knots */


  param_gr = 0;
  if (param_s &&                        /* # semiknots > 0 */
	   hurdlesRU_G.num_great)       /* has greatest real knot */
    param_gr = 1;


  /* Algorithm p. 226, Step 2: build G-bar, and compute its statistics */
  if (build_Gbar(G, pcounts_G, num_chromosomes)) {
    /* only redo classification if paths were closed in building Gbar */

    /* reclassify the components */
    classify_connected_components(G);

    /* compute new statistics */
    calc_num_hurdles_and_fortress(G, C_RU, &hurdlesRU_G);
  }


#if 00
  num_simple = calc_num_simple(G);


  /* Algorithm p. 226, Step 2: build G-bar, and compute its statistics */
  if (num_simple>0) {
    build_Gbar(G, pcounts_G);

    /* reclassify the components */
    classify_connected_components(G);

    /* compute new statistics */
#if 0
    calc_num_hurdles_and_fortress(G, C_IU, &hurdlesIU_G);
#endif
    calc_num_hurdles_and_fortress(G, C_RU, &hurdlesRU_G);
  }
#endif // 00

#if 0
  param_rr = hurdlesRU_G.num_hurdles;   /* # real knots in G-bar */

  /* fr: fortress of real-knots or weak-fortress-of-real-knots */
  param_fr = param_gr = 0;

  if (param_s &&                        /* # semiknots > 0 */
	   hurdlesRU_G.num_great)       /* has greatest real knot */
    param_gr = 1;

  if (hurdlesRU_G.num_fortress)
    param_fr = 1;
  else if (param_gr &&                 /* gr == 1 */
	   (param_rr % 2) == 1 &&      /* odd # real knots */

	   /* and all real knots but the greatest are super-real-knots */
	   (hurdlesRU_G.num_super + 1 ==
	    hurdlesRU_G.num_hurdles))
    param_fr = 1;
#endif

  param_fr = 0;
  if (hurdlesRU_G.num_fortress && !has_sgrk)
    param_fr = 1;
  else if (param_gr &&                 /* gr == 1 */
	   (hurdlesRU_G.num_hurdles % 2) == 1 &&       /* odd # real knots */

	   /* and all real knots but the greatest are super-real-knots */
	   (hurdlesRU_G.num_super + 1 ==
	    hurdlesRU_G.num_hurdles))
    param_fr = 1;



  /* modified distance formula, page 224 */
#if 0
  dist = param_b - param_c + param_p + param_rr
         + (param_s - param_gr + param_fr + 1)/2;
#endif
  dist = param_b - param_c + param_p + param_r
         + (param_s - param_gr + param_fr + 1)/2;


  /* value used later in capping algorithm
   * If reach just 2 semiknots, close them if this flag is true, join
   * them if this flag is false:
   */
  *s2flag = hurdlesRU_G.num_fortress > 0;

#if 0
  *s2flag = param_gr;
#endif

#if 0
  if (verbose) {
#endif
#if 0
    fprintf(outfile,"b=%d c=%d p=%d rr=%d s=%d gr=%d fr=%d     d=%d\n",
	    param_b, param_c, param_p, param_rr, param_s, param_gr, param_fr,
	    dist);
#endif
#if 0
    fprintf(outfile,"Number of black edges in G:   \t%d\n",param_b);
    fprintf(outfile,"Number of cycles and paths:   \t%d\n",param_c);
    fprintf(outfile,"Number of gamma-gamma paths:  \t%d\n",param_p);
    fprintf(outfile,"Number of semi-knots:         \t%d\n",param_s);
    fprintf(outfile,"Number of real knots in G-bar:\t%d\n",param_rr);
    fprintf(outfile,"Parameter gr:                 \t%d\n",param_gr);
    fprintf(outfile,"Parameter fr:                 \t%d\n",param_fr);
    fflush(outfile);
  }
#endif

  /* return the values of the different parameters */
  if (graphstats) {
    graphstats->bl  = param_b;
    graphstats->cp  = param_c;
    graphstats->pgg = param_p;
    graphstats->s   = param_s;
#if 0
    graphstats->rr  = param_rr;
#endif
    graphstats->r   = param_r;
    graphstats->gr  = param_gr;
    graphstats->fr  = param_fr;
    
    graphstats->d   = dist;
  }

  return(dist);
}

/* just do enough of the above alg to produce the graph G(Pi,Gamma) */
void mcdist_noncircular_graph_d(struct mgr_genome_struct *g1,
				struct mgr_genome_struct *g2,
				int num_genes, int num_chromosomes,
				mgr_distmem_t *distmem)
{
  graph_t G0,   *G = &G0;

  /* build the doubled permutations */
  double_perms(g1, g2,
	       0, /* offset */
	       num_genes,
	       distmem);

  /* construct the grey edges, caps, tails (p. 225(d)) */
  buildG_grey(num_genes, num_chromosomes, distmem, G);

  return;
}


/* Do enough of the alg to produce the graph G(Pi,Gamma) and
 * classify the components as IU, RU, or other.
 * This is the amount needed for the unsigned multichromosomal alg.
 */

void mcdist_partial_G(struct mgr_genome_struct *g1,
		      struct mgr_genome_struct *g2,
		      int num_genes, int num_chromosomes,
		      mgr_distmem_t *distmem,
		      graph_t *G,
		      pathcounts_t *pcounts_G)
{
  /* build the doubled permutations */
  double_perms(g1, g2,
	       0, /* offset */
	       num_genes,
	       distmem);

  /* construct the grey edges, caps, tails (p. 225(d)) */
  buildG_grey(num_genes, num_chromosomes, distmem, G);

  /* classify path/cycle types */
  classify_cycle_path(G, pcounts_G);

  /* classify components: IU, RU, other */
  form_connected_components(G);
  classify_connected_components(G);

  /* omit computing actual parameters, or computing hurdle statistics */

  return;
}



/* interface to compute the distance (p. 226, steps 1 & 2)
   that discards the other info not needed for only distance */

int mcdist_noncircular(struct mgr_genome_struct *g1,
		       struct mgr_genome_struct *g2,
		       int num_genes, int num_chromosomes,
		       mgr_distmem_t *distmem,
		       graphstats_t *graphstats)
{
  graph_t G;
  pathcounts_t pcounts_G;
  int s2flag;


  /* allow calling for unichromosomal genomes */
  if (num_chromosomes == 0)
    return invdist_noncircular_v(g1, g2,
				 0, /* offset */
				 num_genes, distmem, graphstats);

  /* multichromosomal */
  return mcdist_noncircular_full(g1, g2,
				 num_genes, num_chromosomes, distmem,
				 &G, &pcounts_G, &s2flag,
				 graphstats);

}


/*
   This does steps 1-18 of page 226.
   It forms the graph from those steps, ending with something like 225(e);
   all the paths have been joined/closed into cycles, but the
   chromos are not yet flipped for the final stage of the capping.

   return value: distance
   also returns the graph in G and the pathcounts in pcounts_G

   does not do the final flippings (steps "18.9" and on)
   or alter the input gene data

*/
int mcdist_capgraph_connections(struct mgr_genome_struct *g1,
				struct mgr_genome_struct *g2,
				int num_genes, int num_chromosomes,
				mgr_distmem_t *distmem,
				graph_t *G, pathcounts_t *pcounts_G,
				graphstats_t *graphstats)
{
  int dist;
  int s2flag;    /* treatment of 2 semiknots */

 
  /* page 226 algorithm, steps 1 & 2 */
  dist = mcdist_noncircular_full(g1, g2,
				 num_genes, num_chromosomes, distmem,
				 G, pcounts_G,
				 &s2flag,
				 graphstats);

  /* step 3: Close all but one Pi-Gamma path in components w >1 such */
  if (pcounts_G->num_pg > 1)
    close_allbut1_pgpath(G, pcounts_G, num_chromosomes);
#if 0
    close_allbut1_pgpath(G, pcounts_G);
#endif

  /* step 5-6: Join pairs of  Pi-Pi & Gamma-Gamma paths */
  join_pp_gg_paths(G, pcounts_G);

  /* steps 4-18 (except 5-6): join/close all Pi-Gamma paths */
  /* join_close_pg_paths(G, pcounts_G); */

  /* new version */
  join_close_semiknots(G, pcounts_G, s2flag);

  return dist;
}


/* return value: distance */
/* also changes the gene data in g1, g2,
 * to be fully capped, properly oriented, concatenated */
int mcdist_capgraph(struct mgr_genome_struct *g1,
		    struct mgr_genome_struct *g2,
		    int num_genes, int num_chromosomes,
		    mgr_distmem_t *distmem,
		    graphstats_t *graphstats)
{
  int dist;
  graph_t G;
  pathcounts_t pcounts_G;



  dist = mcdist_capgraph_connections(g1, g2,
				     num_genes, num_chromosomes,
				     distmem,
				     &G, &pcounts_G,
				     graphstats);

  /* read off the final perms from the graph */
  mcdist_capgraph_finalperms(g1, g2,
			     num_genes, num_chromosomes,
			     distmem,
			     &G,
			     graphstats);


  /* output results by overwriting original genome */
  copy_perm_to_genome(distmem->cappedp1, g1, num_genes);
  copy_perm_to_genome(distmem->cappedp2, g2, num_genes);

  return dist;
}

/* do steps 18.9-19 to compute the final capped permutations */
/* output them in distmem->cappedp1, distmem->cappedp2 */
void mcdist_capgraph_finalperms(struct mgr_genome_struct *g1,
			       struct mgr_genome_struct *g2,
			       int num_genes, int num_chromosomes,
			       mgr_distmem_t *distmem,
			       graph_t *G,
			       graphstats_t *graphstats)
{
  int num_badbonds;

  pathcounts_t pcounts_G;

  classify_cycle_path(G, &pcounts_G);
  form_connected_components(G);
  classify_connected_components(G);

#if 0
  /* step 18.9 (omission): Find a properly-flipped graph (p. 217) */
  find_proper_flipping(G, num_chromosomes);
  /* step 19: find capping of Gamma-hat and Pi-hat given by graph.
     Algorithm in book is incomplete, this routine includes fixes */
  find_induced_capping(G, num_genes, num_chromosomes);
#endif


  /* Replacement for step 19:
     Find optimal concatenation.
     Includes finding proper flipping and proper bonding.
     */
  num_badbonds = find_optimal_bonding(G, num_genes, num_chromosomes);
  compute_concatenates(G, num_genes, num_chromosomes, TRUE);

#if 0    
  if (graphstats && num_badbonds>0) {
    fprintf(outfile,"Number of bad bonds:          \t%d\n",num_badbonds);
    fflush(outfile);
  }
#endif
  if (graphstats) {
    graphstats->badbonds = num_badbonds;
  }

}



int mcdist_capgraph_nomem(struct mgr_genome_struct *g1,
			  struct mgr_genome_struct *g2,
			  int num_genes, int num_chromosomes,
			  graphstats_t *graphstats)
{
  mgr_distmem_t distmem;
  int dist;



  /* allocate memory */
  mcdist_allocmem(num_genes, num_chromosomes,
		  &distmem);


  /* compute capping */
  dist = mcdist_capgraph(g1, g2,
			 num_genes, num_chromosomes,
			 &distmem,
			 graphstats);

  /* free memory */
  mcdist_freemem(&distmem);

  return dist;
}


/* like setinvmatrix, but for multiple chromosomes */
void setmcdistmatrix(int **distmatrix, struct mgr_genome_struct *genomes,
		     int num_genes, int num_chromosomes, int num_genomes,
		     mgr_distmem_t *distmem)
{
  int i, j;

  for (i=0 ; i<num_genomes ; i++) {
    distmatrix[i][i] = 0;
    for (j=i+1 ; j<num_genomes ; j++) {
	distmatrix[i][j] = distmatrix[j][i] =
		mcdist_noncircular(genomes+i,genomes+j,
				   num_genes,num_chromosomes,distmem,
				   (graphstats_t *)NULL); /* not verbose */
    }
  }
  return;
}
