/* graph_components.c
*    Form components of overlap graph.
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

/* Calculations to form connected components,
classify them by type (oriented? intrachromosomal? real?)
and count and classify hurdles, etc.

Contains code from GRAPPA 1.02: invdist.c, routines
connected_component
num_hurdles_and_fortress
*/


#include <stdio.h>
#include <stdlib.h>
#include "mgrstructs.h"
#include "mgr_graph_components.h"

/* invdist.c connected_component with very minor modifications for tails */
/* sets G->num_components */
/* sets  cc and components  arrays */
void form_connected_components(graph_t *G)
{
  int i;    
  int stack_ptr;
  int *greyEdges;   
  int *cycle,*range,*cc;
  mgr_component_t *components;
  int *stack_root,*stack_range;
  int *parent,*next;
  int right,p;
#ifdef TESTING
  struct timeval tm;
  double time0,time1;
#endif

  int size = G->size;
  mgr_distmem_t *distmem = G->distmem;



  /* if no changes to graph, it's already fine. */
  if (G->components_good) return;


#ifdef TESTING
  gettimeofday(&tm,NULL);   
  time0=(double)tm.tv_sec+(double)tm.tv_usec/1000000.0;
#endif
  
  greyEdges            = distmem->greyEdges;

  next = stack_root    = distmem->stack;

  range = distmem->oriented;

  stack_range = cc     = distmem->cc;
  parent=cycle         = distmem->labeled;
  components           = distmem->components;
  stack_ptr = -1;
  G->num_components = 0;


  /* Use Linear algorithm to compute connected component */
#if 0  /* now it's done in classify_cycle_path */
  for (i=0; i<size ; i++) {
    if (greyEdges[i]==V_ADJ)
      continue;
    range[cycle[i]]=i;
  }
#endif

  for (i=0 ; i<size ; i++) {
    /* it is self loop, tail, or lone black edge; discard it*/
    if (greyEdges[i]==V_ADJ || greyEdges[i]==V_TAIL

      /* bare edge: p. 224 -- get 3 semi-knots w/o this, 1 w/this */
	|| (cycle[i]==i && range[i]==i+1)
	|| (cycle[i]==i-1 && range[i-1]==i))  continue;
    if (parent[i]==i) {
      stack_ptr++;   
      stack_root[stack_ptr]=i;
      stack_range[stack_ptr]=range[i];
    }
    else { /*check the top of stack for intersection*/
      right=i; 
#ifdef DEBUG
#if 0
      fprintf(outfile,"top:%d,rang:%d,i:%d,parent[i]:%d\n",
              stack_root[stack_ptr],stack_range[stack_ptr],i,parent[i]);
      fflush(outfile);
#endif // 0
#endif // DEBUG
      // stack_ptr will go down below 0. I don't think it should
      while (/*(stack_ptr>=0) &&*/ (stack_root[stack_ptr] > parent[i])) {
         // fprintf(stdout,"stack_root[%d]=%d > parent[%d] = %d\n",
           //       stack_ptr,stack_root[stack_ptr],i,parent[i]);
         // fflush(stdout);
        /*union top to the i's connected component*/
        parent[stack_root[stack_ptr]]=parent[i];
        if (right < stack_range[stack_ptr])
          right=stack_range[stack_ptr]; /*extend the active range*/
        stack_ptr--;
      }
      if (stack_range[stack_ptr] < right)
        stack_range[stack_ptr]=right;
      if (stack_range[stack_ptr] <=i) {
        /*the top connected-component is INACTIVE*/
        components[G->num_components].index=stack_root[stack_ptr];
        G->num_components++;
        stack_ptr--;
      }
    }
  }


#ifdef DEBUG
#if 0
  fprintf(outfile,"stack_range (later reused as cc): \n    ");
  for (i=0 ; i<size ; i++) {
    fprintf(outfile,"%4d ",stack_range[i]);
  }
  fprintf(outfile,"\n");
#endif // 0
#endif // DEBUG


  /*turn the forest to set of linked lists whose list head is index of
    component*/
  for (i=0 ; i<size ; i++)
    next[i]=-1;

  for (i=0 ; i<size ; i++) {
#if 0
    if (greyEdges[i]==V_ADJ || greyEdges[i]==V_TAIL)
#endif
    /* it is self loop, tail, or lone black edge; discard it*/
    if (greyEdges[i]==V_ADJ || greyEdges[i]==V_TAIL
	|| (cycle[i]==i && range[i]==i+1)
	|| (cycle[i]==i-1 && range[i-1]==i))
      cc[i]=-1;
    else
      if (i!=parent[i]) {
        /* insert i between parent(i) and next of parent(i) */
        next[i]=next[parent[i]];
        next[parent[i]]=i;
      }
  }


  /*label each node with its root*/
  for (i=0; i< G->num_components; i++) {
    p=components[i].index;
  //  fprintf(stdout,"p=%d,i=%d\n",p,i);
    while(p!=-1) {
    //  fprintf(stdout,"cc[p]=i,");
      cc[p]=i;
    //  fprintf(stdout,"p=next[p]=%d\n",next[p]); fflush(stdout);
      p=next[p];
    }
  }

#ifdef DEBUG
  fprintf(outfile,"cc:\n");
  for (i=0;i<size;i++)
    fprintf(outfile,"%4d ",cc[i]);
  fprintf(outfile,"\n");
#endif

  /* components have been computed */
  G->components_good = TRUE;

  /* TODO: clean up, optimize */
  /* Now that this routine is potentially called multiple times
     as edges are added, the original "cycle" array (which is
     renamed and altered into the "parent" array) must be
     restored */
  for (i=0 ; i<size ; i+=2) {
    cycle[i] = cycle[i+1];
  }

#ifdef DEBUG
  fprintf(outfile,"cc: ");
  for (i=0 ; i<size ; i++) {
    fprintf(outfile,"%4d ",cc[i]);
  }
  fprintf(outfile,"\n");

  fprintf(outfile,"next, stack_root, distmem->stack: \n    ");
  for (i=0 ; i<size ; i++) {
    fprintf(outfile,"%4d ",stack_root[i]);
  }
  fprintf(outfile,"\n");

  fprintf(outfile,"range, distmem->oriented: \n    ");
  for (i=0 ; i<size ; i++) {
    fprintf(outfile,"%4d ",range[i]);
  }
  fprintf(outfile,"\n");

  fprintf(outfile,"parent, cycle, distmem->labeled: \n    ");
  for (i=0 ; i<size ; i++) {
    fprintf(outfile,"%4d ",parent[i]);
  }
  fprintf(outfile,"\n");

  fprintf(outfile,"component indices: ");
  for (i=0 ; i < G->num_components ; i++)
    fprintf(outfile,"%3d ",components[i].index);
  fprintf(outfile,"  (%d total)\n",G->num_components);

  fflush(outfile);
#endif


#ifdef TESTING
  gettimeofday(&tm,NULL);
  time1=(double)tm.tv_sec+(double)tm.tv_usec/1000000.0;
  time_linear+=time1-time0;
#endif
  return;
}


/* Classify components according to whether they are
       oriented

       unoriented
       IU              intrachromosomal unoriented
       RU              real unoriented

   Replacement for code from the first half of num_hurdles_and_fortress

   IU, RU are new
*/

void classify_connected_components(graph_t *G)
{
  int i, j;
  int *cc;           /* cc[i] = component # with vertex i */
  int compno;        /* compno = cc[i] */
  mgr_component_t
    *components,
    *comp = (mgr_component_t *) NULL;

  int pt=0, ptmask;

  int *chromNum;     /* chromNum[i]=chromosome # that vertex i is in */
  int *greyEdges;
  int *cycle;
  /*int *right;*/
  int *ptype;

  int size = G->size;
  int num_components = G->num_components;
  mgr_distmem_t *distmem = G->distmem;

  greyEdges  = distmem->greyEdges;
  chromNum   = distmem->chromNum;
  ptype      = distmem->ptype;


  cycle      = distmem->labeled;
  /*right      = distmem->oriented;*/

  cc         = distmem->cc;
  components = distmem->components;

  if (num_components == 0) {
#ifdef DEBUG
    fprintf(outfile,"No components -- quick exit\n");
#endif
    return;
  }

  /* reset all flags for all components */
  for (j=0 ; j < num_components ; j++) {
    components[j].oriented   = FALSE;   /* until proven */
    components[j].intrachrom = TRUE;    /* until disproven */
    components[j].hurdle     = 0;
  }


  /* Calculate whether each gray edge is oriented or unoriented;
     intrachromosomal or interchromosomal;
     and set flags for its component if it's oriented or interchromosomal
   */

  for (i=0 ; i<size ; i++) {
    j = greyEdges[i];

#if 0  /* moved functionality to form_connected_components */
    if (j < 0) {         /* no grey edge */
      /* if component is just one black edge, consider it oriented */
      if (cc[i]>=0 && cycle[i]==i && right[i]==i+1) {
	components[cc[i]].oriented = TRUE;
      }
      continue;
    }
#endif

    /* Note that deleted/ignored edges (j<0)
       are in fact ignored by this, too */
    if (i < j) {

      /* oriented or not? */
      if ((j - i) % 2 == 0) {
	components[cc[i]].oriented = TRUE;
      }

      /* intrachromosomal or not? */
      if (chromNum[i] != chromNum[j]) {
	components[cc[i]].intrachrom = FALSE;
      }
    }
  }

  /* determine the "real" unoriented components */
  /* Scan the vertices from right to left.
   * When a component number is encountered, set the component ACTIVE.
   * When the root vertex of the component is reached, set the component
   * INACTIVE.
   * If any paths are encountered while a component is active, the
   * component is NOT real.
   * Also note presence of gamma-gamma paths while component is active.
   * TODO: consider speeding it up using a stack instead of scanstate flags.
  */

  for (j=0 ; j < num_components ; j++) {
    components[j].scanstate = SCAN_PENDING;

    if (components[j].oriented || !components[j].intrachrom) {

      /* If it's oriented or INTERchromosomal, it's not real */
      components[j].realf = NOTREAL;

    } else { /* intrachromosomal unoriented */
      components[j].realf = 0;
    }

  }

  /* scan vertices left to right to determine which components are real */
  for (i=size-1 ; i>=0 ; i--) {
#if 0
    if (cc[i]<0) continue;  /* vertex isn't in a component */
#endif


    /* activate the component with vertex i if we just encountered it's
       left end */
    if ((compno=cc[i])>=0) {
      comp = &components[compno];
      if (comp->scanstate == SCAN_PENDING)
	comp->scanstate = SCAN_ACTIVE;
    }


    /* if vertex i is left end of a path, set all active components to
     * not real.
     * Also, if path is gamma-gamma or pi-pi path,
     * mark all active components as containing such.
     */

    pt = ptype[i];
    ptmask = 0;
    if (pt == P_PP || pt == P_GG) {
      ptmask = NOTREAL | HASGG;
    } else if (pt == P_PG || pt == P_GP) {
      ptmask = NOTREAL;
    }

    /*    if (greyEdges[i] == V_PICAP || greyEdges[i] == V_GTAIL) {*/
    if (ptmask) {
      for (j=0 ; j < num_components ; j++) {
	if (components[j].scanstate == SCAN_ACTIVE) {
	  /*	  components[j].scanstate = SCAN_DONE; */
	  /*	  components[j].real = FALSE; */
	  components[j].realf |= ptmask;
	}
      }
    }


    /* if vertex i is at the right end of this component,
       deactivate the component */
    if (compno>=0) {
      if (comp->index == i
	  && comp->scanstate == SCAN_ACTIVE) {
	comp->scanstate = SCAN_DONE;
	/* comp->real = TRUE; */  /* defaulted to this */
      }
    }
  }


#ifdef DEBUG
  for (i=0; i < num_components ; i++) {
   fprintf(outfile,"Component %d: index=%d, oriented=%d, intrachrom=%d, real=%d, hurdle=oct'%o\n",
	    i,components[i].index,components[i].oriented,
	   components[i].intrachrom,
	   /*components[i].real,*/
	   components[i].realf,
	   components[i].hurdle);
  }
  fprintf(outfile,"\n");

  fflush(outfile);
#endif

}


/* adapted from num_hurdles_and_fortress */
/* Add in ability to restrict computation to
          all unoriented components;  RU;  IU
*/
  /* flags       oriented     intrachromosomal    real
     type C_U:   FALSE        -                   -
          C_IU:  FALSE        TRUE                -
          C_RU:  FALSE        TRUE                TRUE
  */

/* a bitmask would be more efficient */
#if 0
#define isctype(comp) \
    (comp_type == C_U ? !(comp.oriented) : \
     comp_type == C_IU ? (!(comp.oriented) && comp.intrachrom) : \
                         (!(comp.oriented) && comp.real))
                         /* comp_type == C_RU */
#endif
#define isctype(comp) (!(comp.oriented) && \
    (comp_type == C_U ? TRUE : \
     comp_type == C_IU ? comp.intrachrom : \
                         !comp.realf))
                         /* comp.real)) */
                         /* comp_type == C_RU */

        
/* sets flags in the components[i].hurdle bitfields.
   Sets all the fields of counts to the number of hurdles,
   superhurdles, greatest hurdles, fortresses
   */
void calc_num_hurdles_and_fortress(graph_t *G,
				   int comp_type, /* U, IU, RU */
				   hurdlecounts_t *counts)
{
  int cIdx;
  int i;
  int *cc;
  mgr_component_t *components;
  int num_unoriented;
  int first_comp, last_comp, num_block;

  /* treat hurdles/knots/real knots by reporting status in different bits */
  int hurdle_shift;
  int m_HURDLE;
  int m_HURDLE_GREATHURDLE;
  int m_HURDLE_MASK;
  int m_not_HURDLE_MASK;
  int m_GREATHURDLE;
  int m_SUPERHURDLE;


  int size = G->size;
  int num_components = G->num_components;
  mgr_distmem_t *distmem = G->distmem;

  /* By default, set number of hurdles and fortresses to 0 */
  counts->num_hurdles = 0;
  counts->num_fortress = 0;
  counts->num_super = 0;
  counts->num_great = 0;
  counts->num_unoriented = 0;


  cc     = distmem->cc;
  components           = distmem->components;


  if (num_components == 0) {
#ifdef DEBUG
    fprintf(outfile,"No components -- quick exit\n");
#endif
    return;
  }


  /* Count components of specified type */
  
  num_unoriented = 0;
  for (i=0 ; i<num_components ; i++) {
    /*   if (components[i].oriented == FALSE) */
    if (isctype(components[i])) {
      num_unoriented++;
    }
  }

  counts->num_unoriented = num_unoriented;

  if (num_unoriented == 0) {
#ifdef DEBUG
    fprintf(outfile,"Quick exit -- no hurdles\n");
#endif
    return;
  }

  hurdle_shift = (comp_type == C_U ? U_SHIFT :
		  comp_type == C_IU ? IU_SHIFT :
		  /* comp_type == C_RU ? */
		  RU_SHIFT);

  m_HURDLE = HURDLE << hurdle_shift;
  m_HURDLE_GREATHURDLE = (HURDLE | GREATHURDLE) << hurdle_shift;
  m_HURDLE_MASK = HURDLE_MASK << hurdle_shift;
  m_not_HURDLE_MASK = ~m_HURDLE_MASK;
  m_GREATHURDLE = GREATHURDLE << hurdle_shift;
  m_SUPERHURDLE = SUPERHURDLE << hurdle_shift;
  

  
  for (i=0 ; i<num_components ; i++) {
    components[i].blocks   = 0;
    components[i].hurdle   &= m_not_HURDLE_MASK;  /* ~(HURDLE_MASK << hurdle_shift) */
    components[i].left     = -1;
    components[i].right    = -1;
  }
  


  /* HURDLES 
     Hurdles are a subset of the nonoriented components. 
     There are two types of hurdles (in the KST-sense):
     "simple" and "superhurdle".
     First, we implicitly eliminate oriented components.
     Second, if a nonoriented component is one contiguous block of
     vertices, it is a hurdle.
     Third, if a hurdle "protects" the same non-hurdle on its left and
     right side, then it is a "superhurdle".
  */

  first_comp = -1;
  last_comp = -1;
  num_block = -1;
  for (i=0 ; i<size ; i++) {
    cIdx = cc[i];
    if (cIdx != -1) {
      /* if (components[cIdx].oriented == FALSE) */
      if (isctype(components[cIdx])) {
	if (cIdx != last_comp) {
	  if (last_comp == -1) {
	    first_comp = cIdx;
	  }
	  else {
	    components[last_comp].right = cIdx;
	    components[cIdx].left       = last_comp;
	  }
	  last_comp = cIdx;
	  num_block++;
	  components[cIdx].blocks++;
	}
      }
    }
  }

#if 0
#ifdef DEBUG
  for (i=0 ; i<num_components ; i++) {
    if (components[i].oriented == FALSE) {
      fprintf(outfile,"nonoriented component %3d: blocks: %3d\n",
	      i,components[i].blocks);
    }
  }
#endif // DEBUG
#endif // 0


  /* only components of comp_type get block counts > 0, so no need to test
   * for comp_type membership
   */

  for (i=0 ; i<num_components ; i++) {
    if ( /* (components[i].oriented == FALSE) && */
	(components[i].blocks == 1)) {
#if 0
#ifdef DEBUG
      fprintf(outfile,"nonoriented component %3d: blocks: %3d\n",
	      i,components[i].blocks);
#endif // DEBUG
#endif // 0
      components[i].hurdle |= m_HURDLE;
      (counts->num_hurdles)++;
    }
  }

  if ((first_comp == last_comp) && (components[first_comp].blocks == 2)) {
    components[first_comp].hurdle |= m_HURDLE_GREATHURDLE;
    (counts->num_hurdles)++;
    (counts->num_great)++;
#ifdef DEBUG
    fprintf(outfile,"Component %d is a greatest hurdle\n",first_comp);
#endif
  }

#ifdef DEBUG
  for (i=0; i <num_components ; i++) {
    fprintf(outfile,"Component %d: index=%d, oriented=%d, intrachrom=%d, real=%d, blocks=%d, hurdle=oct'%o, left=%d, right=%d\n",
	    i,components[i].index,
	    components[i].oriented,components[i].intrachrom,components[i].real,
	    components[i].blocks,components[i].hurdle,
	    components[i].left,components[i].right);
  }
#endif


  /* handle special case of hurdle protecting a nonhurdle strictly to
   * its left or right, rather than above it
   */
  if (first_comp != -1 && first_comp < last_comp) {
    /* hurdle on left, nonhurdle on right */
    if (components[first_comp].right == last_comp &&
	(components[first_comp].hurdle & m_HURDLE) &&
	!(components[last_comp].hurdle & m_HURDLE)) {
      components[first_comp].hurdle |= m_SUPERHURDLE;
      (counts->num_super)++;
    } else
      /* hurdle on right, nonhurdle on left */
      if (components[last_comp].left == first_comp &&
	  (components[last_comp].hurdle & m_HURDLE) &&
	  !(components[first_comp].hurdle & m_HURDLE)) {
	components[last_comp].hurdle |= m_SUPERHURDLE;
	(counts->num_super)++;
    }

  }




#if 0
  /* The early returns were because def of fortress is
   *    odd number of hurdles, all of which are superhurdles
   * and the original code only wanted to detect if fortress or not.
   * We want info on how many superhurdles vs. hurdles, so cannot do
   * early returns.
   */

  if (counts->num_hurdles < 3)
    return;
#endif
  
  /* counts->num_super = 0; */
  for (i=0 ; i<num_components ; i++) {
    if (components[i].hurdle & m_HURDLE_MASK ) {
      if ((components[i].left == components[i].right) &&
	  (components[i].left != -1)) {
	if ((components[components[i].left].blocks == 2) &&
	    ((components[components[i].left].hurdle & m_GREATHURDLE) == 0)) {
#ifdef DEBUG
	  fprintf(outfile,"Component %3d is a superhurdle \n",i);
#endif
	  components[i].hurdle |=   m_SUPERHURDLE;
	  (counts->num_super)++;
	}
	else {
	  /*  return; */
	}
      }
      else {
	/* return; */
      }
    }
  }




#if 0
#ifdef DEBUG
  fprintf(outfile,"Number of superhurdles: %3d\n",counts->num_super);
#endif
#endif

  /* Set num_fortress if there are an odd number of hurdles,
     all of which are superhurdles. */
  if ((counts->num_hurdles == counts->num_super)
      && (counts->num_super % 2 == 1))
    counts->num_fortress = 1;

#ifdef DEBUG
  fprintf(outfile,"%s:  # %s=%d   # super=%d   # great=%d    # fortresses=%d\n",
	  comp_type == C_U ? "U" :
	  comp_type == C_IU ? "IU" : "RU",
	  comp_type == C_U ? "hurdles" :
	  comp_type == C_IU ? "knots" : "real knots",
	  counts->num_hurdles,
	  counts->num_super,
	  counts->num_great,
	  counts->num_fortress);

  fflush(outfile);
#endif


  return;
}



/* calculate number of semi-real-knots (Flato's definition to correct
 * the definition of semi-knots)
 * and mark the corresponding components
 */


int calc_num_semirealknots(graph_t *G,
			   int *has_sgrk)
{


  /* detecting semi-greatest-real-knot (s-g-r-k) is as follows.
   * If there is a greatest-real-knot then there is no s-g-r-k.
   * Else, the pattern of vertices and their components in IU is this,
   * after deleting real oriented components:
   *  u1 u2 ...  x .. x  r1 r2 .... x .. x  v1 v2 ...
   * where r1,r2,... are verts of "real" components (>=1 such),
   * u's, v's are verts of nonreal components (>=0 such)
   * and x's are vertices of the s-g-r-k.
   */

  /*
   * Scan the vertices left to right in a manner similar to computing hurdles.
   * The condition for minimum s-r-k is similar to the one for min hurdles.
   * The condition for s-g-r-k is different.
   * The code for super hurdles and fortresses is not needed.
   *
   * Code for computing s-g-r-k uses
   *
   * p_index: last index of an IU\RU component detected prior to the 1st
   *          real knot
   *          = -2: proved that there is no s-g-r-k
   *          = -1: no candidate yet
   *          >= 0: candidate component
   *
   * p_nblocks <= # blocks in the potential s-g-r-k
   *              but in above example would go 0 -> 1 -> 2
   *
   * num_real: number of real components encountered so far
   */




  int cIdx;
  int i;

  int num_IU;
  int first_comp, last_comp, num_block;

  int p_index, p_nblocks, num_real;

  int size = G->size;
  int num_components = G->num_components;
  mgr_distmem_t *distmem = G->distmem;
  int *cc = distmem->cc;
  mgr_component_t *components = distmem->components;
  mgr_component_t *comp;

  int num_srk = 0;
  int has_grk = 0;

  *has_sgrk = 0;

  if (num_components == 0) {
#ifdef DEBUG
    fprintf(outfile,"No components -- quick exit\n");
#endif
    return 0;
  }


  /* count components in IU with no gg/pp paths inside */
  num_IU = 0;
  for (i=0 ; i<num_components ; i++) {
    if (!components[i].oriented
	&& components[i].intrachrom
	&& components[i].realf == NOTREAL   /* bitmask HASGG false */
	) {
      num_IU++;
    }
  }

  if (num_IU == 0) {
#ifdef DEBUG
    fprintf(outfile,"Quick exit -- no semi-real-knots\n");
#endif
    return 0;
  }

  for (i=0 ; i<num_components ; i++) {
    components[i].blocks   = 0;
    components[i].hurdle   &= ~SEMIKNOT;
    components[i].left     = -1;
    components[i].right    = -1;

    if (components[i].hurdle & (GREATHURDLE<<RU_SHIFT)) {
      has_grk++;
    }
  }


  /* Scan for minimal semi-real-knots.
   * Also get info used for computing semi-greatest-real-knot.
   */

  p_index = has_grk ? -2 : -1;
  p_nblocks = 0;
  num_real = 0;

  first_comp = -1;
  last_comp = -1;
  num_block = -1;
  for (i=0 ; i<size ; i++) {
    cIdx = cc[i];

    /* eliminate adjacencies, tails */
    if (cIdx != -1) {
      comp = &components[cIdx];

      /* only examine components in IU */
      if (!comp->oriented && comp->intrachrom) {

	/* only treat 1st vertex of consecutive verts in a component */
	if (cIdx != last_comp) {
	  if (last_comp == -1) {
	    first_comp = cIdx;
	  }
	  else {
	    components[last_comp].right = cIdx;
	    comp->left                  = last_comp;
	  }
	  last_comp = cIdx;
	  num_block++;
	  comp->blocks++;

	  /* if real, does it affect a potential s-g-r-k?
	   * if not real, is it a potential s-g-r-k?
	   */
	  if (p_index > -2) {
	    if (comp->realf) {
	      /* not real, may be potential s-g-r-k */
	      if (num_real == 0) {
		p_index = cIdx;
		p_nblocks = 1;
	      } else {
		if (p_index == cIdx) {
		  p_nblocks++;
		} else if (p_nblocks == 1) {
		  p_index = -2;
		}
	      }
	    } else {
	      num_real++;

	      /* real, see if influences potential s-g-r-k */
	      if (p_index == -1) {
		/* not within a potential s-g-r-k, hence NO s-g-r-k */
		p_index = -2;
	      } else {
		/* possibly within a potential s-g-r-k */
		if (p_nblocks != 1) {
		  /* not within the region since potential s-g-r-k is
		   * first encountered, so it's not an s-g-r-k
		   */
		  p_index = -2;
		}
	      }
	    }
	  }


	}
      }
    }
  }


  /* count and mark minimal semi-real-knots */

  for (i=0 ; i<num_components ; i++) {
    comp = &components[i];
    if (comp->blocks == 1 &&      /* min hurdle */
	comp->realf == NOTREAL    /* has PG/GP path but no GG/PP path */
	) {
      comp->hurdle |= SEMIKNOT;
      num_srk++;
    }
  }


  /* count and mark semi-greatest-real-knot */
  if (p_index >= 0 && num_real > 0) {
    comp = &components[p_index];
    if (comp->blocks == 2 &&     /* possible max hurdle */
	comp->realf == NOTREAL
	) {
      comp->hurdle |= SEMIKNOT;
      num_srk++;
      ++*has_sgrk;
    }
  }


  return num_srk;
}


/* Calculate number of semi-knots, and mark the components as semi-knots.
 * Assumes that the RU and IU hurdles have already been computed:
 *     calc_num_hurdles_and_fortress(G, C_IU, &hurdlesIU_G);
 *     calc_num_hurdles_and_fortress(G, C_RU, &hurdlesRU_G);
 *     calc_num_semiknots(G);
 *
 * Routine is no longer used; superceded by calc_num_semirealknots
 */
int calc_num_semiknots(graph_t *G)
{
  int i,j;
  mgr_component_t *components;
#if 0
  int checkfor;  /* path type to check for */
#endif

  int *cc;
  int *ptype;
  int *greyEdges;

  int num_semiknots;

  mgr_component_t *comp;  /* the component, if any, with the current vertex */

  int size = G->size;
  int num_components = G->num_components;
  mgr_distmem_t *distmem = G->distmem;


  ptype      = distmem->ptype;
  greyEdges  = distmem->greyEdges;
  components = distmem->components;
  cc         = distmem->cc;

#ifdef DEBUG
  for (i=0; i < num_components ; i++) {
   fprintf(outfile,"Component %d: index=%d, oriented=%d, intrachrom=%d, real=%d, hurdle=oct'%o\n",
	    i,components[i].index,components[i].oriented,
	   components[i].intrachrom,components[i].real,
	   components[i].hurdle);
  }
  fprintf(outfile,"\n");

  fflush(outfile);
#endif


#ifdef DEBUG
  fprintf(outfile,"Semiknots @ indices  ");
#endif



  num_semiknots = 0;

  /* semi-knots are a subset of (knots - real knots) */
  for (j=0 ; j<num_components ; j++) {
    if ((components[j].hurdle &
	 ((HURDLE<<IU_SHIFT) | (HURDLE<<RU_SHIFT)))  /* knot | real-knot */
	==
	 ((HURDLE<<IU_SHIFT))) {               /* knot but not real-knot */
      components[j].hurdle |= SEMIKNOT;
      components[j].scanstate = SCAN_PENDING;  /* value may change */
    } else {
      components[j].hurdle &= ~SEMIKNOT;
      components[j].scanstate = SCAN_DONE;     /* it's not a semi-knot */
    }
  }

  /* Eliminate components that have pi-pi or gamma-gamma vertices inside.
   * No IU components have pi-pi paths inside anyway, so only have to
   * eliminate components with gamma-gamma paths.
   */


  for (i=size-1 ; i>=0 ; i--) {
    comp = (cc[i] < 0 ) ?
      (mgr_component_t *)NULL  :  &components[cc[i]];

#if 0
    if (cc[i]<0) continue;  /* vertex isn't in a component */
#endif


    /* activate the component with vertex i if we just encountered its
       left end */
    if (comp != (mgr_component_t *)NULL &&
	comp->scanstate == SCAN_PENDING)
      comp->scanstate = SCAN_ACTIVE;


    /* instead of seeing if this is a pi-pi or gamma-gamma vertex,
       we see if this is the left end of a pi-pi or gamma-gamma path;

       the p-p/g-g vertex @ end of this path is "inside" a component
       <=> this vertex is "inside",
       because if one is inside and one is outside, the path crosses
       the component & is thus part of the component.
       */

    if (    /* ptype[i] == P_PP || */
	ptype[i] == P_GG) {
#if 0
    /* see if this is a pi-pi or gamma-gamma vertex */

    if (greyEdges[i] == V_PICAP)
      checkfor = P_PP;
    else if (greyEdges[i] == V_GTAIL)
      checkfor = P_GG;
    else
      checkfor = P_NONE;   /* not pi-pi or gamma-gamma */

    if (checkfor != P_NONE &&
	(ptype[components[cc[i]].index] == checkfor) ) {
#endif
      /* is pi-pi or gamma-gamma, hence active components are not semi-knots */

      for (j=0 ; j<num_components ; j++) {
	if (components[j].scanstate == SCAN_ACTIVE) {
	  components[j].scanstate = SCAN_DONE;
	  components[j].hurdle &= ~SEMIKNOT;
	}
      }
      
    }

    /* if vertex i is at the right end of this component,
       deactivate this component and adjust semiknot counter,
       if appropriate */
    
    if (comp != NULL &&
	comp->index == i &&
	comp->scanstate == SCAN_ACTIVE) {

      comp->scanstate = SCAN_DONE;
      /* components[cc[i]].hurdle |= SEMIKNOT; */  /* defaulted to this */

#ifdef DEBUG
      fprintf(outfile,"%d ",i);
#endif
      num_semiknots++;
    }
  }

#ifdef DEBUG
  fprintf(outfile," (%d total)\n",num_semiknots);

  fflush(outfile);
#endif

  return num_semiknots;
}


/* determine which components are "simple" */
/* assumes semi-knot flags have already been set */
/* Now only used for testing and for annotating the graph, no longer
 * used in normal course of algorithm */
int calc_num_simple(graph_t *G)
{
  int i;
  mgr_component_t *components;

  int num_simple = 0;

  int *cc;
  int *ptype;

  int size = G->size;
  int num_components = G->num_components;
  mgr_distmem_t *distmem = G->distmem;


  ptype      = distmem->ptype;
  components = distmem->components;
  cc         = distmem->cc;



  /* simple components are components containing pi-gamma paths but
     which are not semi-knots */

  for (i=0 ; i<num_components ; i++ ) {
    components[i].hurdle &= ~SIMPLE;        /* default: not simple */
  }

  for (i=0 ; i<size ; i++) {
    if (ptype[i] == P_PG || ptype[i] == P_GP) {
      /* found a pi-gamma or gamma-pi path */

      if (cc[i] < 0) {
	/* bare edge, treat as simple */
	num_simple++;
      }
      else if ((components[cc[i]].hurdle & SEMIKNOT) == 0)  /* not semiknot */
      {
	components[cc[i]].hurdle |= SIMPLE;
      }
    }
  }
#ifdef DEBUG
  fprintf(outfile,"Simple components:  ");
#endif

  for (i=0 ; i<num_components ; i++) {
    if ((components[i].hurdle & SIMPLE) != 0) {
      num_simple++;
#ifdef DEBUG
      fprintf(outfile,"%d (index %d)  ",i,components[i].index);
#endif
    }
  }
#ifdef DEBUG
  fprintf(outfile,"[total=%d]\n",num_simple);

  fflush(outfile);
#endif

  return num_simple;
}

