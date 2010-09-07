/* graph_edit.c
*    Build the breakpoint graph, and routines to reverse portions of
*    it and add edges.
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

/* contains code from invdist_noncircular in GRAPPA 1.02: invdist.c */

/* build the basic graph, and utility routines to
   reverse portions of it and add edges */
/* contains some code from invdist_noncircular */


#include <stdio.h>
#include <stdlib.h>
#include "mgrstructs.h"
#include "mgr_graph_edit.h"

#include <caml/fail.h>

/* build the graph G(Pi,Gamma) (page 220; page 224(d)) */

void double_perms(struct mgr_genome_struct *g1, struct mgr_genome_struct *g2,
           int offset,
	   int num_genes,
	   mgr_distmem_t *distmem)
{
  int i, twoi, n;
  int g;
  int *perm1, *perm2, *perm;

#ifdef DEBUG
  fprintf(stdout,"double_perms: offset: %3d  num_genes: %3d\n",
	  offset,num_genes);
#endif

  n = 2*num_genes+2;
  
  perm1      = distmem->perm1;          /* doubled (gamma-hat)^(-1) */
  perm2      = distmem->perm2;          /* doubled pi-hat */
  perm       = distmem->perm;           /* composition */

/*for(i=0;i<=2*num_genes;i++)
{
    perm1[i]=0; perm2[i]=0; perm[i]=0;
}*/

  /* create doubled (gamma-hat)^(-1) */
  for (i=0 ; i<num_genes ; i++) {
    g = 2*g1->genes[i];
    twoi = 2*i+1;
    if (g > 0) {
      perm1[g-1] =  twoi;
      perm1[g  ] =  twoi+1;
    }
    else {
      perm1[-g  ] = twoi;
      perm1[-g-1] = twoi+1;
    }
  }
#ifdef DEBUG
  fprintf(stdout,"G1: ");
  for (i=0 ; i<num_genes ; i++) {
    fprintf(stdout,"%4d ",g1->genes[i]);
  }
  fprintf(stdout,"\n");

  fprintf(stdout,"P1: ");
  for (i=1 ; i<=2*num_genes ; i++) {
    fprintf(stdout,"%4d ",perm1[i]);
  }
  fprintf(stdout,"\n");

  fflush(stdout);
#endif



  /* create doubled pi-hat */

  for (i=0 ; i<num_genes ; i++) {
    if (offset<num_genes)
      g = g2->genes[(offset + i)%num_genes];
    else {
      g = - g2->genes[(offset - i)%num_genes];
    }
    twoi = 2*i+1;
    if (g > 0) {
      perm2[twoi  ] =  2*g-1;
      perm2[twoi+1] =  2*g;
    }
    else {
      perm2[twoi  ] = -2*g;
      perm2[twoi+1] = -2*g-1;
    }
  }

#ifdef DEBUG
  fprintf(stdout,"G2: ");
  for (i=0 ; i<num_genes ; i++) {
    fprintf(stdout,"%4d ",
	    offset < num_genes ?
	    g2->genes[(offset+i)%num_genes] :
	    - g2->genes[(offset-i)%num_genes] );
  }
  fprintf(stdout,"\n");

  fprintf(stdout,"P2: ");
  for (i=1 ; i<=2*num_genes ; i++) {
    fprintf(stdout,"%4d ",perm2[i]);
  }
  fprintf(stdout,"\n");

  fflush(stdout);
fprintf(stdout, "i:  ");
  for (i=0 ; i<2*num_genes ; i++) {
    fprintf(stdout,"%4d ",i);
  }
  fprintf(stdout, "\n");

  fflush(stdout);
#endif


  /* composition */

  perm[0] = 0;
  for (i=1 ; i<n-1 ; i++) {
    perm[i] = perm1[perm2[i]];
  }
  perm[n-1] = n-1;

//#ifdef DEBUG
/*
 if((perm[61]==60)&&(perm[63]==60))
{
  fprintf(stdout," P: ");
  for (i=0 ; i<n ; i++) {
    fprintf(stdout,"[%d]%4d ",i,perm[i]);
  }
  fprintf(stdout,"\n");

  fflush(stdout);

   fprintf(stdout,"G1: ");
  for (i=0 ; i<num_genes ; i++) {
    fprintf(stdout,"[%d]%4d ",i,g1->genes[i]);
  }
  fprintf(stdout,"\n");

  fprintf(stdout,"P1: ");
  for (i=1 ; i<=2*num_genes ; i++) {
    fprintf(stdout,"[%d]%4d ",i,perm1[i]);
  }
  fprintf(stdout,"\n");

  fflush(stdout);

  fprintf(stdout,"G2: ");
  for (i=0 ; i<num_genes ; i++) {
    fprintf(stdout,"[%d]%4d ",i,
	    offset < num_genes ?
	    g2->genes[(offset+i)%num_genes] :
	    - g2->genes[(offset-i)%num_genes] );
  }
  fprintf(stdout,"\n");

  fprintf(stdout,"P2: ");
  for (i=1 ; i<=2*num_genes ; i++) {
    fprintf(stdout,"[%d]%4d ",i,perm2[i]);
  }
  fprintf(stdout,"\n");

  fflush(stdout);
}*/
//#endif

}



/* adapted from 2nd part of num_cycles */
/* builds several arrays:
   cycle = distmem->labeled:  (THIS INCLUDED PATHS TOO!!!!)
      cycle[j] = i   when j is in a cycle/path whose leftmost vertex is i
                 P_NONE
		     when j is in an adjacency or tail (which are ignored)
		     i.e. when greyEdges[i] == V_ADJ or V_TAIL

   several arrays, wasteful of space.  The leftmost vertex "i" in each
   cycle/path will have an entry filled in for ARRAY[i],
   and the other vertices will have some "NOT APPLICABLE" marker filled in

      right[i] = j       IF i is the leftmost vertex in a cycle, set this to
                         j, where j is the rightmost vertex
		 P_NONE  when i is NOT the leftmost vertex in its cycle

      Bcap[i]  = P_NONE  when i is NOT the leftmost vertex in its cycle/path
                 otherwise:
		 P_CYCLE This is a cycle
		 j        j is the last vertex on following the path from i
		          starting with the black edge first
      Gcap[i]  = P_NONE   when i is NOT the leftmost vertex in its cycle/path
                 otherwise:
		 P_CYCLE  This is a cycle
		 j        j is the last vertex on following the path from i
		          starting with the grey edge first,
			  or i if there is no grey edge associated with i

      ptype[i] = P_NONE     not the leftmost vertex of a path/cycle
           if i is the leftmost vertex of a cycle, set ptype[i]
	   to one of these:
                 P_CYCLE    
                 P_PP
                 P_PG
                 P_GP
                 P_GG
		 

   Returns several counts in the fields of pcounts:
      num_adjacencies (# V_ADJ V_ADJ pairs)
      num_cycles  (ignoring V_ADJ's and V_TAIL's)
      num_pp      number of pi-pi paths
      num_gg      number of gamma-gamma paths
      num_pg      number of pi-gamma (or gamma-pi) paths
*/

void classify_cycle_path(graph_t *G,
			 pathcounts_t *pcounts)
{
  int *greyEdges;
  int *done;
  int *cycle;
  int *right;
  int *ptype;
  int *Bcap, *Gcap;

  int i;
  //int j;
  int next;

  int rightmost; /* the rightmost vertex in a cycle/path */
  int Bpos, Gpos;
  int endtag;

  //G = &G_mcdist_noncircular;

  int size = G->size;
  mgr_distmem_t *distmem = G->distmem;

  int n_adjacencies, n_cycles, n_pp, n_pg, n_gg;

  n_adjacencies = n_cycles = n_pp = n_pg = n_gg = 0;

  G -> components_good = FALSE; /* components not yet computed */
/*
  fprintf(stdout, "greyEdges[size=%d] = \n",size);
  for (i=0 ; i<size ; i++)
  {
      fprintf(stdout, "%d,",distmem->greyEdges[i]);
  }
  fprintf(stdout, "\n"); fflush(stdout);
*/
  greyEdges  = distmem->greyEdges;

  done       = distmem->done;
  right      = distmem->oriented;
  cycle      = distmem->labeled;

  Bcap       = distmem->Bcap;
  Gcap       = distmem->Gcap;
  ptype      = distmem->ptype;

  /* Initialization */
  /* Count adjacencies */
  /* Mark adjacencies and tails as done, everything else as not done */
  /* Clear the other arrays */
  for (i=0 ; i<size ; i++) {
    done[i] = 0;
    if (greyEdges[i] == V_ADJ) { 
      n_adjacencies++;
      done[i] = 1;
    } else if (greyEdges[i] == V_TAIL) {
      done[i] = 1;
    };
    cycle[i] = right[i] = Bcap[i] = Gcap[i] = ptype[i] = P_NONE;
  };

  n_adjacencies /=2;


   
#if 0
  /* old code to mark the vertices in a cycle */
  for (i=0 ; i<size ; i++) {
    if (done[i] == 0) {
      cycle[i]=i;
      done[i] = 1;
      next = i;
      do {
	if (next % 2 == 0)
	  next++;
	else
	  next--;
	done[next] = 1;
        cycle[next]=i;
	next = greyEdges[next];
	done[next] = 1;
	cycle[next]=i;
      }
      while(next != i);
      c++;
    }
  }
#endif

  /* Scan left to right for the leftmost vertices of paths/cycles that
     haven't been examined.
     On finding one, traverse the whole path or cycle to
     mark it (by leftmost vertex), to compute its rightmost vertex,
     its caps (if a path), and its classification.
     */

  for (i=0; i<size; i++) {
    if (done[i] == 0) {
    
#if 000
      if (greyEdges[i] < 0) {
	/* It cannot be a cycle, start from the black edge */
	Gpos = i;
	Bpos = 0;
	endtag = 1;
	cycle[i] = i;
	done [i] = 1;
	next = i;
	rightmost = i;
	do {
	  /* follow black edge */
	  if (next % 2 == 0)  next++;	  else   next--;

	  done[next] = 1;
	  cycle[next] = i;
	  rightmost = rightmost >= next ? rightmost : next;

	  /* follow gray edge */
	  if (greyEdges[next] < 0) {
	    endtag = 0;
	    Bpos = next;
	  }
	  else {
	    next = greyEdges[next];
	    done[next] = 1;
	    cycle[next] = i;
	    rightmost  = rightmost >= next ? rightmost : next;
	  }
	}
	while(endtag);

      } else {
#endif // 000


	/* It can be a path or a cycle. We will start from the black edge.
	   If it reaches the beginning, it is a cycle.
	   Otherwise, we will start from the grey edge again. */
	Bpos = Gpos = 0;
	endtag = 1;
	done[i] = 1;
	cycle[i] = i;
	rightmost = i;
	next = i;
	do {
	  /* follow black edge */
	  if (next % 2 == 0)  next++;	  else   next--;
	  done[next] = 1;
	  cycle[next] = i;
	  rightmost = rightmost >= next ? rightmost : next;

	  /* follow gray edge */
	  if (greyEdges[next] < 0) {
	    endtag = 0;
	    Bpos = next;
	  } else {
	    next = greyEdges[next];
	    done[next] = 1;
	    cycle[next] = i;
	    rightmost  = rightmost >= next ? rightmost : next;
//	  if (next % 2 == 0)  next++;	  else   next--;
	  }   
	} while (endtag   &&   next != i);

	if (next == i) {
	  /* it is a cycle */
	  Bpos = Gpos = P_CYCLE;
	} else {
	  /* path, follow other half, starting from the gray edge */
	  next  = i;
	  endtag = 1;
      

      int count = 0;
	  do { 
          count ++;
          if(count >1000)
          {
              fprintf(stdout, "\ngreyEdges[%d]=%d;",next,greyEdges[next]);
              fflush(stdout);
          }
	    /* gray edge */
	    if (greyEdges[next] < 0) {
	      endtag = 0;
	      Gpos = next;
	    } else {
            if(count > 1000)
            {
                fprintf(stdout, "next<--%d,done[%d]<--1,cycle[%d]<--%d;",
                        greyEdges[next],next,next,i);
            }
	      next = greyEdges[next];
	      done[next] = 1;
	      cycle[next] = i;
	      rightmost = rightmost >= next ? rightmost : next;
	      /* black edge */
	      if (next % 2 == 0)  next++;	  else   next--;
          if(count>1000) { fprintf(stdout, "next <-- %d ;",next); fflush(stdout); }
	      done[next] = 1;
	      cycle[next] = i;
	      rightmost = rightmost >= next? rightmost : next;
	    };

      if(count>1005)
      {
            fprintf(stdout, "i=%d,greyEdges[size=%d] = \n",i,size);
              for (i=0 ; i<size ; i++)
           {
          fprintf(stdout, "[%d]:%d,",i,distmem->greyEdges[i]);
          }
      fprintf(stdout, "\n"); fflush(stdout);
      failwith("endless loop");
      }


	  } 
	  while(endtag); 
	};
#if 000
      };
#endif // 000
      right[i] = rightmost;
      Gcap[i] = Gpos;
      Bcap[i] = Bpos;
      if (Gpos>=0 && Bpos>=0) {
	if (greyEdges[Bpos] == V_PICAP) {
	  if (greyEdges[Gpos] == V_PICAP) {
	    n_pp++;
	    ptype[i] = P_PP;
	  } else { /* greyEdges[Gpos] == V_GTAIL */
	    n_pg++;
	    ptype[i] = P_PG;
	  }
	} else { /* greyEdges[Bpos] == V_GTAIL */
	  if (greyEdges[Gpos] == V_PICAP) {
	    n_pg++;
	    ptype[i] = P_GP;
	  } else { /* greyEdges[Gpos] == V_GTAIL */
	    n_gg++;
	    ptype[i] = P_GG;
	  }
	}
#if 0
	if (greyEdges[Gpos]==V_PICAP && greyEdges[Bpos]==V_PICAP) {
	  n_pp++;
	  ptype[i] = P_PP;
	};
	if (greyEdges[Gpos]==V_GTAIL && greyEdges[Bpos]==V_PICAP) {
	  n_pg++;
	  ptype[i] = P_PG;
	};
	if (greyEdges[Gpos]==V_GTAIL && greyEdges[Bpos]==V_GTAIL) {
	  n_gg++;
	  ptype[i] = P_GG;
	};
	if (greyEdges[Gpos]==V_PICAP && greyEdges[Bpos]==V_GTAIL) {
	  n_pg++;
	  ptype[i] = P_GP;
	};
#endif
      } else {
	n_cycles++;
	ptype[i] = P_CYCLE;
      };
    }     

  };

  pcounts->num_adjacencies = n_adjacencies;
  pcounts->num_cycles = n_cycles;
  pcounts->num_pg = n_pg;
  pcounts->num_pp = n_pp;
  pcounts->num_gg = n_gg;

#ifdef DEBUG
  fprintf(stdout,"classify_cycle_path:\n");
  fprintf(stdout,"     num_adjacencies=%d, num_cycles=%d, num_gg=num_pp=%d, num_pg=%d, total=%d\n",
	  n_adjacencies,
	  n_cycles,
	  n_pp,
	  n_pg,
	  n_adjacencies + n_cycles + 2*n_pp + n_pg);

  print_arrays(G);
#endif

}



/* flip chromosome cnum_start,...,cnum_end
   The graph consists of cycles only, no paths.

   Chromosomes occupy vertices vLeft,vLeft+1,...,vRight
   where greyEdges[vLeft...vRight] = V_TAIL <other entries> V_TAIL
   and the other entries are either vertex #'s (>=0) or V_ADJ

   The reversal operation should update the following arrays:
       perm2:     entries in positions vLeft...vRight should be reversed
                  i.e., swap perm2[vLeft+i] & perm2[vRight-i] for
		  i=0,...,(vRight-vLeft)/2
       greyEdges:
                  Need to mirror image the grey edge incidencies within
		  chromosome c.  To do this:
                  1. reverse the entries between vLeft...vRight as in perm2.
		  2. for ALL entries i=0,...,size-1,
		        if vLeft<=greyEdges[i]<=vRight,
		           replace greyEdges[i] by Mirror-greyedges[i]
		           where Mirror=vLeft+vRight
		     Negative entries V_ADJ, V_TAIL may occur, and
		     should not be changed.
		     Other negative entries V_???? should not occur at all.
       cycle, right:
                  update A[i]  (i=0,...,size-1,  A=one of these arrays)
		  when its value is between vLeft & vRight;
		  don't waste time computing it if it's not in that range.
		  1. for i=0,...,size-1
		         done[i] = FALSE if A[i] is between vLeft & vRight,
			           TRUE  otherwise
		  2. Do the original cycle marking algorithm
		     (but also keeping track of what to put in "right")
		  
       Bcap, Gcap, ptype:
                  There are no paths, only cycles.
		  So, for i=vLeft...vRight, set
		  ptype[i]=Bcap[i]=Gcap[i]=P_CYCLE
		             if i is leftmost vertex in cycle after flipping
		                          =P_NONE   otherwise


       all entries (indices 0,...,size-1) whose values
                  are between vLeft & vRight need to be recomputed.



       chromosome boundaries are also moved if flipping consecutive
       chromosomes:

       chromNum:
             the chromosome numbers of each vertex must be mirrored
       chromBd:
             the chromosome boundaries must be mirrored
*/
 
   
/* G = graph,  cnum = chromosome # */
/* Boolean updatecycleinfo: TRUE to update cycle/path info, FALSE to not */

void flip_chromosomes(graph_t *G, int cnum_start, int cnum_end,
		      int updatecycleinfo)
{
  int *greyEdges;
  int *done;
  int *cycle;
  int *right;
  int *ptype;
  int *perm2;
  int *Bcap, *Gcap;
  int *chromBd, *chromNum;

  int vLeft, vRight, Mirror;
  int i,j;
  int temp;
#if 00
  int next,next1,rightmost;
#endif // 00
  int size = G->size;

  mgr_distmem_t *distmem = G->distmem;

  pathcounts_t pcounts_G;


  perm2      = distmem->perm2;          /* doubled pi-hat */
  greyEdges  = distmem->greyEdges;

  done       = distmem->done;
  right      = distmem->oriented;
  cycle      = distmem->labeled;

  Bcap       = distmem->Bcap;
  Gcap       = distmem->Gcap;
  ptype      = distmem->ptype;

  chromBd    = distmem->chromBd;
  chromNum   = distmem->chromNum;

  vLeft      = chromBd[cnum_start];
  vRight     = chromBd[cnum_end+1]-1;
  Mirror     = vLeft+vRight;
#ifdef DEBUG
  fprintf(stdout,"Chromosomes %d-%d go between vertices %d and %d\n",
	  cnum_start, cnum_end, vLeft, vRight);
  print_arrays(G);
#endif // DEBUG

  /* modify perm2 */
  for (i=vLeft, j=vRight;
       i<j;
       i++, j--) {
    temp = perm2[i];
    perm2[i] = perm2[j];
    perm2[j] = temp;
  };

#ifdef DEBUG
  fprintf(stdout,"perm2 done\n");
#endif // DEBUG

#if 0
  for (i=0; i<size; i++) {
    done[i] = 0;
  }; 
#endif

  /* modify greyEdges */
  for (i=vLeft, j=vRight;
       i<j;
       i++, j--) {
    temp = greyEdges[i];
    greyEdges[i] = greyEdges[j];
    greyEdges[j] = temp;
  };

  for (i=0; i<size; i++) {
      if (greyEdges[i]>=vLeft && greyEdges[i]<=vRight) {
	/* destination is inside the range */
	greyEdges[i] = Mirror-greyEdges[i];
      }
  };

#ifdef DEBUG
  fprintf(stdout,"greyEdges done\n");
#endif // DEBUG

  if (updatecycleinfo) {

    /* TODO: fix the disabled "00" code correctly to incrementally
       fix up the various arrays.
       The cycle and right arrays are not being properly updated
       by the code because of other routines that reuse these
       arrays for other purposes... need to track this down.
       Meanwhile, update the info inefficiently by one call:
       */

    classify_cycle_path(G,&pcounts_G);

    /* Since chromosome flipping is only done when capping,
       this inefficiency probably won't matter.
       Ordinary distance computations will never execute this.
       */

#if 00

#if 0
    for (i = 0; i< size; i++)
      {
	if (greyEdges[i]<0) done[i] = 1;
	Gcap[i] = Bcap[i] = ptype[i] = P_NONE;
	cycle[i] = right[i] = P_NONE;
      };
#endif

    for (i=0 ; i<size ; i++) {
      /* only need to update a cycle if either of its ends
	 was in the flipped chromosomes */
      if ((vLeft<=cycle[i] && cycle[i]<=vRight) ||
	  (vLeft<=right[i] && right[i]<=vRight)) {
	done[i] = FALSE;
	Gcap[i] = Bcap[i] = ptype[i] = P_NONE;
	cycle[i] = right[i] = P_NONE;
      } else done[i] = TRUE;
    }

    /* modify cycle/right */
    for (i=0; i<size; i++) {
      if (done[i] == 0) {
	cycle[i] = i;
	done[i] = 1;
	next = i;
	rightmost = i;
	do {
	  /* follow black edge */
	  if (next % 2 == 0)  next++;	  else   next--;

	  if (next > rightmost) rightmost = next;

	  done[next] = 1;
	  cycle[next]=i;

	  /* follow grey edge */
	  next1 = greyEdges[next];
	  if (next1>=0)
	    next = next1;
	  else if (next1 == V_ADJ || next1 == V_TAIL) {
	    if (next % 2 == 0)  next++;	  else   next--;
	  } else { /* PI-CAP or GAMMA-TAIL, except there should only be
		      cycles now */
	    break;
	  }

	  if (next > rightmost) rightmost = next;
	  done[next] = 1;
	  cycle[next]= i;
	}
	while(next != i);
	right[i] = rightmost;
	ptype[i] = Gcap[i] = Bcap[i] = P_CYCLE;
      }
    };
#endif // 00
  }

  /* flipped multiple chromosome, must update chromosome boundaries */
  /* ex. bdries at  ... start:10 16 30 42 end:60  100 (next one) ... */
  /*     goes to    ...       10 50 68 80 94      100 ... */
  if (cnum_start != cnum_end) {

#ifdef DEBUG
    fprintf(stdout,"before reflection: chromNum:\n    ");
    for (i=0 ; i<size ; i++) {
      fprintf(stdout,"%4d ",chromNum[i]);
    }
    fprintf(stdout,"\n");

    fprintf(stdout,"chromBd:  ... ");
    for (i=cnum_start ; i<=cnum_end+1 ; i++) {
      fprintf(stdout,"%d ",chromBd[i]);
    }
    fprintf(stdout," ...\n");
#endif


    /* first, swap the boundaries */
    for (i=cnum_start+1, j=cnum_end;
	 i<j;
	 i++, j--) {
      temp = chromBd[j];
      chromBd[j] = chromBd[i];
      chromBd[i] = temp;
    }

    /* then, mirror the boundaries */
    for (i=cnum_start+1; i<=cnum_end; i++)
      chromBd[i] = Mirror - chromBd[i] + 1;

    /* then, update the chromosome numbers */
    for (i=cnum_start; i<=cnum_end; i++) {
      for (j=chromBd[i]; j<=chromBd[i+1]-1; j++)
	chromNum[j]=i;
    }

#ifdef DEBUG
    fprintf(stdout,"after reflection: chromNum:\n    ");
    for (i=0 ; i<size ; i++) {
      fprintf(stdout,"%4d ",chromNum[i]);
    }
    fprintf(stdout,"\n");

    fprintf(stdout,"chromBd:  ... ");
    for (i=cnum_start ; i<=cnum_end+1 ; i++) {
      fprintf(stdout,"%d ",chromBd[i]);
    }
    fprintf(stdout," ...\n");
#endif
  }

  /* some interchromosomal components split/merged */
  G->components_good = FALSE;

#ifdef DEBUG
  fprintf(stdout,"Arrays after flipping:\n");
  print_arrays(G);
#endif
}

/* flip chromsomes w/o updating path/cycle info */
void flip_chromosomes_nocycles(graph_t *G, int cnum_start, int cnum_end)
{
  flip_chromosomes(G, cnum_start, cnum_end, FALSE);
}

/* flip chromosome with vertex # i */
void flip_chromosome_wvertex(graph_t *G, int i)
{
  int cnum = G->distmem->chromNum[i];
  flip_chromosomes(G, cnum, cnum, TRUE);
}


void test_flipping(graph_t *G, int num_chromosomes)
{
  int j;
  print_arrays(G);
  for (j=1 ; j<=num_chromosomes ;  j++) {
    fprintf(stdout,"Flipping chromosome %d\n",j);
    flip_chromosomes(G,j,j, TRUE);
  }
}

/***************************************************************************/

/************************************************************************/
/* Do surgery on the graph */

/* close the path whose leftmost vertex is   index  */
/* Does NOT update the classification of this component */
void closepath(graph_t *G, int index, pathcounts_t *pcounts)
{
  int i,j;

  int *Bcap, *Gcap;
  int *greyEdges;
  int *ptype;
  int *cycle;

  mgr_distmem_t *distmem = G->distmem;

  Bcap       = distmem->Bcap;
  Gcap       = distmem->Gcap;
  greyEdges  = distmem->greyEdges;
  ptype      = distmem->ptype;
  cycle      = distmem->labeled;



  /* Closing path doesn't alter components. */
  /* It could affect the classification of it being real or not real. */
  /* G->components_good = FALSE; */




  /* get the vertices at the two ends of the path */
  i = Gcap[index];
  j = Bcap[index];

#ifdef DEBUG
  fprintf(stdout,"Closing path: index=%d   edge=(%d,%d)\n",
	  index,i,j);

  fflush(stdout);
#endif

  /* adjust path count */
  if (ptype[index] == P_PG || ptype[index] == P_GP) {
    pcounts->num_pg--;
  } else {
    if (ptype[index] == P_PP)
      pcounts->num_pp--;
    else
      pcounts->num_gg--;         /* should never happen */
  }

  /* add a grey edge */
  if (i == index  &&  j == i+1) { /* forms an adjacency once closed */
    greyEdges[i] = greyEdges[j] = V_ADJ;
    cycle[index] = P_NONE;        /* it's no longer regarded as a cycle */

    /* They are no longer ends of a path */
    Gcap[index] = Bcap[index] = P_NONE;

    /* reclassify the left vertex as NOT starting a path or cycle */
    ptype[index] = P_NONE;

    /* adjust path/cycle counts */
    pcounts->num_adjacencies++;
#if 0
    pcounts->num_pg--;
#endif
  } else {                        /* doesn't form an adjacency */
    /* new edge (i,j) */
    greyEdges[i] = j;
    greyEdges[j] = i;

    /* Now it's a cycle, so it doesn't have path ends */
    Gcap[index] = Bcap[index] = P_CYCLE;

    /* reclassify the path as a cycle */
    ptype[index] = P_CYCLE;

    /* adjust path/cycle counts */
    pcounts->num_cycles++;
#if 0
    pcounts->num_pg--;
#endif
  }
}


/* Add grey edge between vertices i & j in graph G.
   On input, one of i or j is a pi-cap, the other a gamma-tail.
   They must lie on different paths.
   No error checking is done to enforce these conditions.
*/
void addedge(graph_t *G, int i, int j, pathcounts_t *pcounts)
{
  int *greyEdges;
  int *right;
  int *cycle;
  int *Bcap, *Gcap;
  int *ptype;

  int tmp;
  int next;
  int i_is_B, j_is_B;
  int index_i, index_j;
  int new_Bcap, new_Gcap;

  mgr_distmem_t *distmem = G->distmem;

  greyEdges  = distmem->greyEdges;
  right      = distmem->oriented;
  cycle      = distmem->labeled;
  Bcap       = distmem->Bcap;
  Gcap       = distmem->Gcap;
  ptype      = distmem->ptype;


  /***********************************************************************/

  G->components_good = FALSE;          /* component computations invalidated */

  greyEdges[i] = j; greyEdges[j] = i;  /* insert edge */

  /* cycle/path index = leftmost vertex */
  index_i = cycle[i]; index_j = cycle[j];

  /* swap indices & i,j if necessary so that index_i < index_j */
  if (index_i > index_j) {
    tmp = i; i = j; j = tmp;
    tmp = index_i; index_i = index_j; index_j = tmp;
  }

#ifdef DEBUG
  fprintf(stdout,"Joining paths: indices=%d,%d   edge=(%d,%d)\n",
	  index_i,index_j,i,j);

  fflush(stdout);
#endif


  /* everything in path j should now be indexed by index_i instead */
  next = j;
  do {
    cycle[next]=index_i;            /* relabel index of vertex */

    if (next % 2 == 0)  next++;	  else   next--;  /* follow black edge */

    cycle[next]=index_i;            /* relabel index of vertex */
    next = greyEdges[next];         /* follow grey edge */
  } while (next >= 0);  /* terminate on reaching vertex w/o grey edge */


  /* are i,j the B-cap (TRUE) or G-cap (FALSE) of their paths? */
  i_is_B =  Bcap[index_i] == i;
  j_is_B =  Bcap[index_j] == j;

  /* compute the new B-cap and G-cap of the merged paths */
  /* depends on the 4 combinations of which ends of each path i,j are on */
  if (i_is_B) {
    new_Gcap = Gcap[index_i];
    if (j_is_B) {
      new_Bcap = Gcap[index_j];
    } else {
      new_Bcap = Bcap[index_j];
    }
  } else {
    new_Bcap = Bcap[index_i];
    if (j_is_B) {
      new_Gcap = Gcap[index_j];
    } else {
      new_Gcap = Bcap[index_j];
    }
  }
  Bcap[index_i] = new_Bcap; Gcap[index_i] = new_Gcap;
  Bcap[index_j] = P_NONE;   Gcap[index_j] = P_NONE;


  /* compute the new rightmost vertex */
  if (right[index_j] > right[index_i])  right[index_i] = right[index_j];
  right[index_j] = P_NONE;   /* and j no longer has a rightmost vertex */


  /* new path type */
  /* The path types joined are either
       pi-pi    & gamma-gamma  ->  pi-gamma  (or gamma-pi)
   OR  pi-gamma & pi-gamma     ->  pi-gamma  (or gamma-pi)
   */
  /* first update the path counts */
  if (ptype[index_i] == P_PP || ptype[index_i] == P_GG) {
    pcounts->num_pp--; pcounts->num_gg--; pcounts->num_pg++;
  } else {
    pcounts->num_pg--;
  }
  /* and then the path types */
  ptype[index_i] =  (greyEdges[new_Bcap] == V_PICAP) ? P_PG : P_GP;
  ptype[index_j] = P_NONE;
}

/***************************************************************************/

/* Join two paths into a cycle by adding edges (v1,v2) and (v3,v4)
   In each edge, one vertex is a pi-cap and the other a gamma-tail
   and the added edges must be at opposite ends of two paths.
   No error checking is done to enforce these conditions.

   Does NOT update the classification of this component.

   Adds in the edges.
   Updates the cycle array for the upper path.
   Updates the value of right.
   Updates Bcap, Gcap, ptype.
   Updates pcounts.
*/
void join_paths_2_cycle(graph_t *G,
			int v1, int v2, int v3, int v4,
			pathcounts_t *pcounts)
{
  int *greyEdges;
  int *right;
  int *cycle;
  int *Bcap, *Gcap;
  int *ptype;

  int tmp;
  int next;
//  int i_is_B, j_is_B;
  int index_i, index_j;
//  int new_Bcap, new_Gcap;

  mgr_distmem_t *distmem = G->distmem;

  greyEdges  = distmem->greyEdges;
  right      = distmem->oriented;
  cycle      = distmem->labeled;
  Bcap       = distmem->Bcap;
  Gcap       = distmem->Gcap;
  ptype      = distmem->ptype;


  /***********************************************************************/

  G->components_good = FALSE;          /* component computations invalidated */

  /* one path is "i", other path is "j".
     These names are symbollic only.
     */
//  greyEdges[i] = j; greyEdges[j] = i;  /* insert edge */
  greyEdges[v1] = v2; greyEdges[v2] = v1; /* insert first edge */

  /* cycle/path index = leftmost vertex */
//  index_i = cycle[i]; index_j = cycle[j];
  index_i = cycle[v1]; index_j = cycle[v2];

  /* swap indices & v1,v2 if necessary so that index_i < index_j */
  if (index_i > index_j) {
//    tmp = i; i = j; j = tmp;
    tmp = v1; v1 = v2; v2 = v1;
    tmp = index_i; index_i = index_j; index_j = tmp;
  }

#ifdef DEBUG
  fprintf(stdout,"Joining paths: indices=%d,%d   edge=(%d,%d)\n",
	  index_i,index_j,i,j);

  fflush(stdout);
#endif


  /* everything in path j should now be indexed by index_i instead */
//  next = j;
  next = v2;
  do {
    cycle[next]=index_i;            /* relabel index of vertex */

    if (next % 2 == 0)  next++;	  else   next--;  /* follow black edge */

    cycle[next]=index_i;            /* relabel index of vertex */
    next = greyEdges[next];         /* follow grey edge */
  } while (next >= 0);  /* terminate on reaching vertex w/o grey edge */


  greyEdges[v3] = v4; greyEdges[v4] = v3; /* insert second edge */


//  /* are i,j the B-cap (TRUE) or G-cap (FALSE) of their paths? */
//  i_is_B =  Bcap[index_i] == i;
//  j_is_B =  Bcap[index_j] == j;
//
//  /* compute the new B-cap and G-cap of the merged paths */
//  /* depends on the 4 combinations of which ends of each path i,j are on */
//  if (i_is_B) {
//    new_Gcap = Gcap[index_i];
//    if (j_is_B) {
//      new_Bcap = Gcap[index_j];
//    } else {
//      new_Bcap = Bcap[index_j];
//    }
//  } else {
//    new_Bcap = Bcap[index_i];
//    if (j_is_B) {
//      new_Gcap = Gcap[index_j];
//    } else {
//      new_Gcap = Bcap[index_j];
//    }
//  }
//  Bcap[index_i] = new_Bcap; Gcap[index_i] = new_Gcap;
//  Bcap[index_j] = P_NONE;   Gcap[index_j] = P_NONE;

    Bcap[index_i] = Gcap[index_i] = P_CYCLE; /* leftmost vertex in cycle */
    Bcap[index_j] = Gcap[index_j] = P_NONE;  /* not leftmost vertex any more */
    


  /* compute the new rightmost vertex */
  if (right[index_j] > right[index_i])  right[index_i] = right[index_j];
  right[index_j] = P_NONE;   /* and j no longer has a rightmost vertex */


  /* new path type */
  /* The path types joined are either
       pi-pi    & gamma-gamma  ->  cycle
   OR  pi-gamma & pi-gamma     ->  cycle
   */
  /* first update the path counts */
  if (ptype[index_i] == P_PP || ptype[index_i] == P_GG) {
    pcounts->num_pp--;
    pcounts->num_gg--;
    pcounts->num_cycles++;
  } else {
    pcounts->num_pg -= 2;
    pcounts->num_cycles++;
  }
  /* and then the path types */
  ptype[index_i] = P_CYCLE;
  ptype[index_j] = P_NONE;
}





/* i, j are the leftmost vertices of pi-gamma or gamma-pi paths.
   Form a cycle by joining each path's pi-cap with the other's gamma-tail.
   */

void join_pgpaths_2_cycle(graph_t *G,
			  int i, int j,
			  pathcounts_t *pcounts) {
//  int i,j;
//  int index_i;
  int v1i=0, v2i=0, v1j=0, v2j=0;

//  int *greyEdges, *cc;
  int *ptype, *Bcap, *Gcap;
//  component_t *components;

//  int size = G->size;
  mgr_distmem_t *distmem = G->distmem;

  ptype      = distmem->ptype;
//  greyEdges  = distmem->greyEdges;
//  components = distmem->components;
//  cc         = distmem->cc;
  Bcap       = distmem->Bcap;
  Gcap       = distmem->Gcap;

#ifdef DEBUG
  fprintf(stdout,"join_pgpaths_2_cycle: i=%d j=%d\n", i, j);
  fflush(stdout);
#endif


  /* must join pi-cap with gamma-tail */
  /* get it so that v1i, v1j have same type (pi or gamma),
     v2i, v2j have the other type */
  /* if both paths are pi-gamma or both are gamma-pi, use one perm of
     the vertices;
     if one is pi-gamma and other is gamma-pi, use the other perm.
     */
  if (ptype[i] == ptype[j]) {
    v1i = Bcap[i]; v2i = Gcap[i];   /* ends of 1st path */
    v1j = Bcap[j]; v2j = Gcap[j];   /* ends of 2nd path */
  } else {
    v1i = Bcap[i]; v2i = Gcap[i];   /* ends of 1st path */
    v1j = Gcap[j]; v2j = Bcap[j];   /* ends of 2nd path */
  }

  /* now join (v1i,v2j) and (v1j,v2i) */
  join_paths_2_cycle(G, v1i,v2j, v1j,v2i, pcounts);
}


/***************************************************************************
 * Compute distance on cycle starting at i, following black edge, gray edge,
 * etc., till reach j.
 * With paths, may have to try i, gray, black, gray, ..., too.
 * d(i,i) = 0
 * if i and j are in different paths/cycles, returns -1
 ***************************************************************************/

int graphdist(graph_t *G, int i, int j)
{
  mgr_distmem_t *distmem = G->distmem;
  int *cycle         = distmem->labeled;
  int *greyEdges     = distmem->greyEdges;

  int c_i,c_j,d,v;

  c_i = cycle[i];
  c_j = cycle[j];

  /* different paths or cycles */
  if (c_i != c_j) return -1;
  if (i == j) return 0;

  if (c_i < 0) {
    /* it's an adjacency or tail, check specially */
    if (i % 2 == 1) {
      return (i + 1 == j) ? 1 : -1;
    } else {
      return (i - 1 == j) ? 1 : -1;
    }
  }

  /* same path or cycle */
  d = 0;
  v = i;
  while (v >= 0) {
    /* follow black edge */
    if (v % 2 == 0) v++; else v--;
    d++;

    if (v == j) return d;

    /* follow gray edge */
    v = greyEdges[v];
    d++;
    if (v == j) return d;

    if (v == i) {
      fprintf(stderr,"ERROR: vertex j misclassified, i=%d, j=%d, c_i=%d\n", i, j, c_i);
      exit(-1);
    }
  }

  /* if we got here w/o a bug, it means it's a path and
   * that j is on the other side of i.
   */
  d = 0;
  v = i;
  while (v >= 0) {
    /* follow gray edge */
    v = greyEdges[v];
    d++;
    if (v == j) return d;
    if (v < 0) {
      fprintf(stderr,"ERROR: vertices misclassified, i=%d, j=%d, c_i=%d\n", i, j, c_i);
      exit(-1);
    }

    /* follow black edge */
    if (v % 2 == 0) v++; else v--;
    d++;

    if (v == j) return d;
  }

  /* should never get here */

  fprintf(stderr,"ERROR: graphdist, shouldn't be here\n");
  exit(-1);


}




/***************************************************************************/

/* Debugging: print various arrays */

/* 0-offset array */
void print_int_array(char *name, int *p, int size) {
  int i;

  fprintf(stdout,"%s:\n",name);
  for (i=0; i<size; i++)
    fprintf(stdout,"%4d ",p[i]);

  fprintf(stdout,"\n");

}
/* 1-offset array */
void print_int_array1(char *name, int *p, int size) {
  int i;

  fprintf(stdout,"%s:\n     ",name);
  for (i=1; i<size; i++)
    fprintf(stdout,"%4d ",p[i]);

  fprintf(stdout,"\n");

}



void print_arrays(graph_t *G)
{
  int size = G->size;

  mgr_distmem_t *distmem = G->distmem;


  /*done       = distmem->done;*/



  print_int_array("cycle=distmem->labeled", distmem->labeled, size);
  print_int_array("right=distmem->oriented", distmem->oriented, size);
  print_int_array("Bcap", distmem->Bcap, size);
  print_int_array("Gcap", distmem->Gcap, size);
  print_int_array("ptype", distmem->ptype, size);
  print_int_array1("perm2", distmem->perm2, size-1);
  print_int_array("greyEdges (E)", distmem->greyEdges, size);
  print_int_array("cc", distmem->cc, size);

  fflush(stdout);
}

