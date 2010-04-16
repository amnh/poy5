/* mgrstructs.h
*
*  SUMMARY:
*
*		Structures specific to MGR program.
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
#ifndef MGRSTRUCTS_H
#define MGRSTRUCTS_H


#include <math.h>
#include <string.h>
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#ifdef DARWIN
#include <getopt.h>
#endif
#include "structs.h"
//#include "mgr_mcstructs.h"
//#include "mgr_mcrdist.h"
//#include "mgr_check.h"
//#include "mgr_scenario.h"

#define MAXINT 999999
#define MAXANCESTOR 3000
#define MAXNEWICK 10000
// the initital maximum size of a strip for alphabet (it is increased if necessary later on)
#define MAXSTRIP 2

/* definitions for types of genomes */
#define GCIRCULAR   0	/* for unichromosomal circular genomes */
#define GLINEAR		1	/* for unichromosomal linear oriented genomes */
#define GMULTI		2	/* for multichromosomal or linear unoriented genomes */

/* structs below are copy from mgr_mcstructs.h(mcstructs.h is the original name)*/


/* greyEdges[i] values for non-edges */
/* Grey edge handling
      It was: greyEdges[i]=j with j>=0:  grey edge from i to j
                                  j=-1:  adjacency, generally ignored
   Multichromosome case:

      tail pi-cap inside-verts pi-cap tail
   and if one of the inside-verts is connected by grey edge to pi-cap,
   the grey edge is removed and this inside-vert becomes a gamma-tail
*/
#define V_ADJ    -1    /* adjacency between 2i & 2i+1 */
#define V_TAIL   -2    /* pad at end of chromosome; delete from graph */
#define V_PICAP  -3    /* connected w/black edge to actual end of chromosome */
#define V_GTAIL  -4    /* needs connection w/grey edge to some pi-cap */

/* cycle and path types */
#define P_NONE    -1   /* vertex is NOT the leftmost one in its path/cycle */
#define P_CYCLE   -2   /* it's a cycle */
#define P_PP      -3   /* Bcap is Pi-cap,      Gcap is Pi-cap */
#define P_PG      -4   /* Bcap is Pi-cap,      Gcap is Gamma-tail */
#define P_GP      -5   /* Bcap is Gamma-tail,  Gcap is Pi-cap */
#define P_GG      -6   /* Bcap is Gamma-tail,  Gcap is Gamma-tail */

/* hurdle flags: */
#define HURDLE          1
#define GREATHURDLE (1<<1)
#define SUPERHURDLE (1<<2)
#define HURDLE_MASK (HURDLE|GREATHURDLE|SUPERHURDLE)
#define SEMIKNOT    (1<<9)
#define SIMPLE      (1<<10)
#define U_SHIFT  0
#define IU_SHIFT 3
#define RU_SHIFT 6

/* values of scanstate: */
#define SCAN_PENDING   0   /* haven't yet activated component */
#define SCAN_ACTIVE    1   /* activated component and have not yet
			      seen enough to determine the final answer */
#define SCAN_DONE      2   /* have final answer for this component */

/* REAL flags: */
#define NOTREAL  1
#define HASGG    (1<<1)



typedef struct {
	int size;
	int last;
	int *strip;
} a_strip;


struct mgr_genome_struct {
	int *genes;
	int genome_num;
	char *encoding;
	char *gnamePtr;       /* Used to copy the gname when adding genomes to the tree */
	int num_g;            /* size of permutation (including caps */
	int num_chr;          /* number of chromosome in this genome */
	int ori_num_g;        /* original number of genes in this genome before alphabet conversion */
	a_strip *alphabet;	/* each gene can actually be a strip, alphabet stores conversion table */
    int *delimiters;
    int num_delimiters;
	
};

/* Structs for invdist.c */
typedef struct {
  int index;    /* index of component's root */
  int oriented; /* Boolean: Is component oriented */
  int intrachrom;
                /* Boolean: Is component intrachromosomal? */

#if 0
  int real;     /* Boolean: Is component "real"? */
#endif
  int realf;    /* Bitmask of REAL properties 
		 * Bit 0: 1=has pi-cap or gamma-tail inside, so not real
		 * Bit 1: 1=has gg path inside
		 */

  int scanstate; /* used while scanning components */

  int blocks;   /* Number of blocks in nonoriented component */
  int hurdle;   /* Bitmask of HURDLE properties.
		 * Set U (unoriented components):
		 * Bit 0 = hurdle
		 * Bit 1 = wrap-around hurdle
		 * Bit 2 = superhurdle

		 * Set IU (intrachromosomal unor. components): bits 3,4,5
		 * Set RU (real intrachrom. unor. components): bits 6,7,8
		 */

  int left;     /* Index of component to the left of my rightmost block */
  int right;    /* Index of component to the right of my rightmost block */
} mgr_component_t;


typedef struct {
#if 0
  int *hammingArr; /* DIST_BP:  (NUMGENES+1)*2 */
#endif
  
  int *perm1;      /* DIST_INV: 2*num_genes + 2 */
  int *perm2;
  int *perm;
  int *done;
  int *greyEdges;

  int *ptype;      /* path/cycle type */
  int *Bcap;       /* vertex @ end of path when start @ left, follow
		      black edge, and continue until reach the end */
  int *Gcap;       /* same, but starting @ left and following grey edge */
  int *chromNum;   /* chromosome # 1,...,N with the vertex
		      0 at start, N+1 at end vertex
		   */


  int *chromBd;   /* chromosome i=1,...,N (in doubled pi) occupies
		     vertices chromBd[i],chromBd[i]+1,...,chromBd[i+1]-1
		     size = num_chromosomes+2
		  */

  int *capMark;   /* normalized caps i=0,1,...,2*num_chromosomes-1

		     in compute_concatenates:
		     Array used to mark that normalized cap i has been
		     used.

		     in find_optimal_bond:
		     Array used to store which caps are at opposite ends
		     of same chromosome strings
		     */

  int *cappedp1, *cappedp2;   /* size = num_genes */
  
  int *stack;      /* DIST_INV: size */
  int *oriented;
  int *cc;
  int *labeled;
  mgr_component_t *components;

#if 0
  UFelem  *uf;
#endif

} mgr_distmem_t;

mgr_distmem_t MGR_DISTMEM;

typedef struct {
  int size;             /* # vertices */
  int num_components;   /* # components */

  int components_good;  /* Boolean:
			 * TRUE if components have been computed/updated
			 * FALSE if not yet computed/updated
			 *
			 * This only describes the accuracy of the
			 * components, not their classsification flags.
			 */


  mgr_distmem_t *distmem;   /* all the structures and arrays defining the graph */
} graph_t;


typedef struct {
  int num_pp;           /* # pi-pi paths */
  int num_gg;           /* # gamma-gamma paths */
  int num_pg;           /* # pi-gamma + gamma-pi paths */
  int num_cycles;       /* # cycles other than adjacencies */
  int num_adjacencies;  /* # adjacencies = 1 black + 1 grey edge */
} pathcounts_t;


typedef struct {
  int num_hurdles;
  int num_great;
  int num_super;
  int num_fortress;

  int num_unoriented;
} hurdlecounts_t;

/* graph parameters for verbose output */
typedef struct {
  /* all */
  int n;
  int d;

  /* unichromosomal */
  int br;                /* "b" = # breakpoints */
  int c4;                /* "c" = # cycles except adjacencies */
  int h;                 /* h = # hurdles */
  int f;                 /* f = # fortresses */

  int num_components_u;  /* number of unoriented components */
  int num_components_o;  /* number of oriented components */
  /* adjacencies are not counted in either of those, but can be computed
   * as   n + 1 - br
   */

  /* multichromosomal */
  int bl;                /* "b" = # black edges */
  int cp;                /* "c" = # cycles and paths, of all lengths */
  int pgg;               /* "p" = # gamma-gamma paths */
  int s;                 /* s = # semi-knots (semi-real-knots) */
#if 0
  int rr,fr,gr;          /* parameters about G-bar */
#endif
  int r;                 /* r = # real-knots */
  int fr,gr;


  int badbonds;          /* number of bad bonds, usu. 0, sometimes 1 */
  int bp_int;            /* number of internal breakpoints */
  int bp_ext;            /* number of external breakpoints */

  /* (un)signed */

} graphstats_t;


/* structs above are copy from mgr_mcstructs.h*/



/* list of edges in the phylogeny */
typedef struct a_list_edge{
  int node1;		/* node is used internally to keep track of genomes */
  int label1;		/* label corresponds to the actual node # for output */
  int node2;		/* second genome involved in this edge */
  int label2;		/* ... */
  int score;		/* number of rearrangements on this edge */
  struct a_list_edge *next;
} list_edge;


/* data structure for the phylogeny */
typedef struct {
    int *tree_array;					/* 1: this genome is the tree
                    					 * 0: this gneome is not yet in the tree */
    int tree_size;					/* nb of genomes already in the tree */
    int new_node;						/* index of the next empty node */
    list_edge *the_edge_list;		/* list of edges of the tree on which we 
                             		 * can add new genomes */
    list_edge *merge_edge_list;	/* list of edges of the tree generated by
                               	 * merges while looking for good rearrangements
											 * we can't add genomes to these edges */
} treemem_t;



/* list of possible rearrangements */
typedef struct a_list_reag{
  int spec;     		/* this rearrangement is for which genome */
  int optype;			/* what type of rearrangement is it */
  int start; 			/* we start the reversal at... */
  int end;				/* we end the reversal at... */
  int sc1;				/* indicates first chromo involved, <0 means it's reverse */
  int sc2;				/* indicates second chromo ... */
  int merge_with;		/* if this rearrangement is carried on to spec it becomes 
                 		 * identical to genome #merge_with */
  int reduction;		/* if this rearrangement is carried on the
                		 * total distance is reduced by ... */
  int *changes;			/* changes associated with this rearrangemnt */
  int Rscore;           /* counts how many times rearrangements involving the same edges
	                     * were seen in the list */
  
  struct a_list_reag *next;
} list_reag;

/* contains all the info for the genomes */
typedef struct {
	struct mgr_genome_struct *genome_list;	/* the actual genomes */
	int *label;								/* label used in treemem for output */
    int *same_as;								/* this genome is now the same as... */
    int *nb_chromo;							/* nb of chromosomes per genome */
    int *dist_mat;							/* pairwise distance matrix */
	int max_chromo_size;                /* the size of the longest chromosome in genome_list */
} G_struct;





extern FILE *ancestor_file;


#endif /* ifndef MGRSTRUCTS_H*/
