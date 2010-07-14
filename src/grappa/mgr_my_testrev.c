/* my_testrev.c
*
*  SUMMARY:
*		
*		Modified version of the testrev.c file from Glenn Tesler.
*
*      Do all possible multichromosomal ops between g1 and g2         
*      Original code from Glenn Tesler, modified so 
*      that it produces a list of good rearrangements      
*
* Copyright (C) 2004-2006 Genome Insitute of Singapore
* Copyright (C) 2000-2002 University of Southern California
*
* Modified from file testrev.c in GRIMM. See file COPYRIGHT for details. 
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

/* #include <stdio.h> */

#include <stdlib.h>
#include "mgrstructs.h"
#include "mgr_mcrdist.h"
#include "mgr_mcread_input.h"
#include "mgr_e_malloc.h"
#include "mgr_write_data.h"
#include "mgr_genome_ops.h"
#include "mgr_list_ops.h"
#include "mgr_my_testrev.h"

#include <caml/fail.h>

typedef struct {
	struct mgr_genome_struct *trial_genome;
	struct mgr_genome_struct *dest_genome;
	int num_genes;
	int num_chromosomes;
	mgr_distmem_t *distmem;
	
//	int counts[CNF0][CNF1]; nobody use this array, just get rid of it.
} tryrevparams_t;

int my_try_operation(int optype, int d, tryrevparams_t *p);

int try_good_operation(int optype, G_struct *Genomes, int nb_spec, int spec_left,
					   int gindex1, tryrevparams_t *p, 
					   int *merge_with, int reduction_type, int *changes,
					   int pair1, int pair2);	
    


/* init chromosome boundary arrays */
/* init_cbounds: allocates memory for ptrs in cb
* init_cbounds_wmem: the pointers in cb already point to
* blocks of memory large enough
*/

void init_cbounds_wmem(int num_genes, int num_chromosomes,
					   struct mgr_genome_struct *g1,
					   cbounds_t *cb) {
	
    if(cb == NULL) fprintf(stderr, "`error cb in init_cbounds_wmem\n");
    if ((cb->cBound == NULL) || (cb->cBound == (int*) NULL)) fprintf(stderr, "`error cBound in init_cbounds_wmem\n");
     if ((cb->cNum == NULL) || (cb->cNum == (int*) NULL)) fprintf(stderr, "`error cNum \n");
fflush(stderr);

    if(cb->num_genes<num_genes) 
        fprintf(stderr, "Warning in init_cbounds_wmem: num_gene we need is %d, bigger than %d\n", num_genes,cb->num_genes);
    if(cb->num_chromosomes < num_chromosomes)
        fprintf(stderr, "Warning in init_cbounds_wmem: num_chromosome=%d>%d\n", num_chromosomes,cb->num_chromosomes);
    fflush(stderr);

	int i;
	int inchrom, cnum;
	int lowcap;
	int *genes, g;
	
	lowcap = num_genes - 2*num_chromosomes + 1;
	cnum = 0;
	inchrom = FALSE;
	genes = g1->genes;
	
	for (i=0 ; i<num_genes ; i++) {
		g = genes[i];
		if (g<0)   g = -g;
		if (g >= lowcap) {
			if (inchrom) {
				/* end of chromosome */
				inchrom = FALSE;
			} else {
				/* start of chromosome */
				inchrom = TRUE;
				cb->cBound[cnum] = i;
				cnum++;
			}
		}
		cb->cNum[i] = cnum;
	}
	
	if (num_chromosomes > 0)
    {
		cb->cBound[cnum] = num_genes;  /* end of last chromosome */
    }
	
#if 0
	if (num_chromosomes == 0) {
		cb->cBound[0] = 0;
		cb->cBound[1] = num_genes;
	} else {
		cb->cBound[cnum] = num_genes;
	}
#endif
	
#if 0
    fprintf(stdout,"\ncb->cNum: ");
	for (i=0 ; i < num_genes ; i++) {
		fprintf(stdout,"%d ",cb->cNum[i]);
	}
	fprintf(stdout,"\n");
	
	fprintf(stdout,"\ncb->cBound: ");
	for (i=0 ; i <= num_chromosomes ; i++) {
		fprintf(stdout,"%d ",cb->cBound[i]);
	}
	fprintf(stdout,"\n");
#endif
	
}

void init_cbounds(int num_genes, int num_chromosomes,
				  struct mgr_genome_struct *g1,
				  cbounds_t *cb) {
    fprintf(stdout,"init_cbounds, num_genes/num_chromo=%d/%d, count_cbounds=%d\n",
            num_genes, num_chromosomes,count_cbounds);
    //fflush(stdout);
//	cb->cNum     = (int *) e_malloc(num_genes*sizeof(int), "cb->cNum");
//	cb->cBound   = (int *) e_malloc((num_chromosomes+1)*sizeof(int),
//									"cb->cBound");
    count_cbounds ++;
	cb->cNum     = (int *) malloc(num_genes*sizeof(int));
	cb->cBound   = (int *) malloc((num_chromosomes+1)*sizeof(int));
	init_cbounds_wmem(num_genes, num_chromosomes, g1, cb);
}

/*
void init_cbounds_memory (int num_genes, int num_chromosomes, cbounds_t * cb)
{
    fprintf(stdout,"init cbounds memory, num_genes=%d\n",num_genes);
    cb  = (cbounds_t *) malloc( sizeof(cbounds_t));
    if( (cb == NULL) || ( cb == ( cbounds_t *)NULL))
        failwith ("initialization for cbouns_t failed");
    cb->num_genes = num_genes;
    cb->num_chromosomes = num_chromosomes;
    cb->cNum = (int *) malloc(num_genes*sizeof(int));
    if( (cb->cNum == NULL) || (cb->cNum == (int *)NULL))
        failwith ("initialization for cbouns_t.cNum failed");
	cb->cBound   = (int *) malloc((num_chromosomes+1)*sizeof(int));
    if( (cb->cBound == NULL) || (cb->cBound == (int *)NULL)) 
        failwith ("initialization for cbouns_t.cBound failed");
}
*/

void free_cbounds(cbounds_t *cb) {
   // fprintf(stdout,"free_cbounds,count_cbounds=%d\n",count_cbounds);
    if(cb!=NULL)
    {
        count_cbounds --;
	    free(cb->cBound);
	    free(cb->cNum);
        free(cb);
    }
    else 
        fprintf(stderr, "Warning: try to free NULL cbounds\n");
}



/*********************************************************/
/* create/maintain isBP "array":                        */
/* If a genome has consecutive genes (a,b) or (-b,-a),  */
/* then the position between the two genes is called:   */
/* "after a", "before b", "after -b", "before -a"       */
/* and isBP(-a)=isBP(b) is boolean valued:              */
/*             =1  if the position is a breakpoint      */
/*             =0  if it is not a breakpoint            */
/*                                                      */
/* Special gene values:                                 */
/*     isBP(num_genes+1) = 1      if end is b.p.        */
/*     isBP(-(num_genes+1)) = 1   if start is b.p.      */
/*                                                      */
/* maintain it as a bit array                           */
/* for genes g=x,   x=1,...,num_genes, use bit 2x-2     */
/* for genes g=-x,  x=1,...,num_genes, use bit 2x-1     */
/********************************************************/

/* allocate memory for breakpoint array.
Return pointer  (char *) *breakpoints
*/

#define BPsize(num_genes)   (num_genes/4 + 1)

void init_isBP_mem(int num_genes, int num_chromosomes,
				   char **breakpoints)
{
	/* # bits required = 2*(num_genes+1), rounded up to a byte */
	/* # bytes = floor((2*(num_genes+1) + 7)/8) = floor(num_genes/4) + 1 */
	*breakpoints = (char *) e_malloc( BPsize(num_genes) * sizeof(char),
									  "breakpoints");
	
}

void clean_isBP_mem(char *breakpoints)
{
	free(breakpoints);
}

/* set or clear whole b.p. array */
/* value = 0 or 1 */
void set_isBP_array_whole(int num_genes, char *breakpoints, int value)
{
	int bp_bytes = BPsize(num_genes);
	int i;
	char setto = value ? 0xff : 0;
	
	for (i=0; i<bp_bytes; i++)
		breakpoints[i] = setto;
}


/* set_isBP_array_1: set bit for gene g in b.p. array
set_isBP_array_0: clear bit for gene g
*/
void set_isBP_array_1(int g, char *breakpoints)
{
	/* for genes g=x,   x=1,...,num_genes, use bit 2x-2     */
	/* for genes g=-x,  x=1,...,num_genes, use bit 2x-1     */
	
	int bitno = (g>0)  ?   (g-1)<<1  :  (((-g) - 1)<<1)  + 1 ;
	
	
	/* break it up into byte number and bit w/in byte */
	int byteno = bitno >> 3;
	bitno = bitno & 7;
	
	
	breakpoints[byteno] |=   1<<bitno;
}


/* breakpoint determination method 2:
Caps are 2-sided b.p.
For adjacencies not involving caps, if (a,b)/(-b,-a) is in all genomes,
there's no b.p. before b or -a, else there is.

In unichromo case, genome starts (ends) w/b.p. iff any genomes differ in
starting (ending) gene
*/
void init_isBP_array_m2(int num_genes, int num_chromosomes,
						int num_genomes,
						G_struct *Genomes,
						int gindex,
						char *breakpoints)
{
	struct mgr_genome_struct *genome_list = Genomes->genome_list;
	
	int lowcap = num_genes - 2*num_chromosomes + 1;
	
	int gindex2;         /* secondary genome */
	int i;               /* loop over gene positions */
	int gene1, gene2;
	
	int *successors;     /* SUCCESSOR array for primary genome */
	
	
	/* determine successor array of primary genome */
	/* The memory allocated for this could be done just once
		if this routine is repeatedly called */
	init_scenario_bp(num_genes, num_chromosomes,
					 &genome_list[gindex],
					 0,             /* no pad, it's already accounted for */
					 &successors);
	
	/* clear all b.p.'s */
	set_isBP_array_whole(num_genes, breakpoints, 0);
	
	
	
	/* if multichrom, the start and end will always be considered b.p.
		These values are never examined, but to be consistent with
		how the positions between internal caps are treated, call them b.p.
		*/
	
	if (num_chromosomes>0) {
		/* b.p. after start position */
		set_isBP_array_1(-(num_genes+1), breakpoints);
		
		/* b.p. before end position */
		set_isBP_array_1(num_genes+1, breakpoints);
	}
	
	
	
	
	/* Form union of b.p. over all pairs (primary genome, other genome) */
	for (gindex2 = 0; gindex2 < num_genomes; gindex2++) {
		/* skip if this is the primary genome */
		if (gindex2 == gindex)
			continue;
		
		// GB: new in MGR_1.35, BP array shouldn't consider merged genomes
		if (Genomes->same_as[gindex2] != -1) {
			//fprintf(stdout, "in BP array, skipping merged genome\n");
			continue;
		}
		
		/* add b.p. for this genome to the breakpoint list */
		/* first check internal b.p., and then check start/end b.p. */
		for (i=0; i<num_genes-1; i++) {
			gene1 = genome_list[gindex2].genes[i];
			gene2 = genome_list[gindex2].genes[i+1];
			
#if 0
			fprintf(stdout,
					"g1=%d, g2=%d, SUCC(g1)=%d, SUCC(-g2)=%d",
					gene1, gene2, SUCCESSOR(-gene1), SUCCESSOR(-gene2));
#endif
			
			/* caps are always b.p.
				adjacencies not involving caps are ignored */
			
			if (SUCCESSOR(gene1) == gene2              /* pair is adjacency */
				&&
				gene1 < lowcap && gene1 > -lowcap      /* gene1 not a cap */
				&&
				gene2 < lowcap && gene2 > -lowcap      /* gene2 not a cap */
				) {
#if 0
				fprintf(stdout,"\n");
#endif
				continue;
			}
#if 0
fprintf(stdout,"   B.P.\n");
#endif


/* add b.p. to list */
set_isBP_array_1(gene2, breakpoints);  /* b.p. before gene2 */
set_isBP_array_1(-gene1, breakpoints); /* b.p. after gene1 */
		}


/* now check start/end b.p. */

if (num_chromosomes > 0) {
	/* multichrom: these values should never be examined.
	To be consistent w/the positions between internal cocaps,
	call it a breakpoint.
	*/
	
	
	/* b.p. before leading cap */
	gene2 = genome_list[gindex2].genes[0];
	set_isBP_array_1(gene2, breakpoints);
	
	/* b.p. after start position: already set */
	/* set_isBP_array_1(-(num_genes+1), breakpoints); */
	
	
	/* b.p. after ending cap */
	gene2 = genome_list[gindex2].genes[num_genes-1];
	set_isBP_array_1(-gene2, breakpoints);
	
	/* b.p. before end position: already set */
	/* set_isBP_array_1(num_genes+1, breakpoints); */
	
} else {
	/* unichrom: b.p. at start (end) if genomes start (end) w/different genes
	*/
	
	/* check start */
	gene1 = genome_list[gindex].genes[0];
	gene2 = genome_list[gindex2].genes[0];
	
	if (gene1 != gene2) {
		/* b.p. before first gene */
		set_isBP_array_1(gene2, breakpoints);
		
		/* b.p. after start position */
		set_isBP_array_1(-(num_genes+1), breakpoints);
	}
	
	
	
	/* check end */
	gene1 = genome_list[gindex].genes[num_genes-1];
	gene2 = genome_list[gindex2].genes[num_genes-1];
	
	if (gene1 != gene2) {
		/* b.p. after last gene */
		set_isBP_array_1(-gene2, breakpoints);
		
		/* b.p. before end position */
		set_isBP_array_1(num_genes+1, breakpoints);
	}
}


	} /* for (gindex2...) */



clean_scenario_bp_mem(successors);
}


/* initialize the breakpoint array for a list of genomes */
void init_isBP_array(int num_genes, int num_chromosomes,
					 int num_genomes,
					 G_struct *Genomes,
					 int gindex, /* which genome to compute b.p. against */
					 char *breakpoints)
{
#if 0
	int i;
#endif
	
	/* pick method _m1, _m2, or _m3 */
	init_isBP_array_m2(num_genes,num_chromosomes,
					   num_genomes, Genomes,
					   gindex,
					   breakpoints);
	
#if 0
	fprintf(stdout,"Breakpoint array: ");
	for (i=0; i<BPsize(num_genes); i++)
		fprintf(stdout," %02x", (unsigned char) breakpoints[i]);
	fprintf(stdout,"\n");
#endif
}

/* check if a position is a breakpoint */
int isBP(int g, char *breakpoints)
{
	/* for genes g=x,   x=1,...,num_genes, use bit 2x-2     */
	/* for genes g=-x,  x=1,...,num_genes, use bit 2x-1     */
	
	int bitno = (g>0)  ?   (g-1)<<1  :  (((-g) - 1)<<1)  + 1 ;
	
	
	/* break it up into byte number and bit w/in byte */
	int byteno = bitno >> 3;
	bitno = bitno & 7;
	
	/* retrieve bit */
	return (breakpoints[byteno] >> bitno) & 1;
}



						
/*****************************************************************************/
/*            Do all possible multichromosomal ops between g1 and g2         */
/*            Original code from Glenn Tesler, modified by Guillaume Bourque */
/*            So that it produces a list of good rearrangements              */
/*****************************************************************************/
/*   	      The variable 'output_coded' was removed from the function      */
/*	      Because the format of the output is not relevant here	     */


int build_list_reag(list_reag **a_list, G_struct *Genomes , int nb_spec, int spec_left,
					int reversals,            // TRUE if we want to include reverals
					int transloc, 	// TRUE if we want to include translocations
					int fusion,
					int fission, 	// TRUE if we want to include fusions and fissions
					mgr_distmem_t *distmem,
					int only_size, 			// don't build the list, just return it's size
					int optimize,   // we are optimizing the median
					int reduction_type,	/* we only accept rearrangements reducing the total distance by at least
					* this amount ( = spec_left-1 for good rearrangements)*/
					int rev_size,  // for heurisitic 5, only look for reversal of this size
					int last_genomenb,/* the last genome in which we did a rearrangement, we want to alternate
					                  how we start the list */
					int pair1, int pair2,  //usually both -1 but use them if we want to inforce a merge
					int verbose)
{
   if(verbose) { fprintf(stdout, "build list reag, num_spec = %d, rev_size = %d, last_genomenb = %d (%d,%d,%d,%d)\n", nb_spec, rev_size, last_genomenb, reversals, transloc, fusion, fission);
   fflush(stdout);}
    struct mgr_genome_struct *genome_list = Genomes->genome_list;
    //	int *dist_mat = Genomes->dist_mat;		  /* GB: the matrix with the pairwise distances */
    struct mgr_genome_struct trial_genome;
    list_reag *the_list = NULL;
    int score_list = 0;
	int *changes; // the changes associated with a particular rearrangement

    /* GB: The following are not necessary since we use genomes directly

        (no need to add caps anymore)
        mgr_distmem_t distmem;
    structmgr_genome_struct current_genome;
    structmgr_genome_struct Dest_genome, *dest_genome=&Dest_genome;
    */

    int optype;

    //	int num_bkpt;

    int foundany_i;      /* has a reversal been found for this start? */
    int i,j;
    /* GB: d is remove because we use dist_mat in Genomes
        int d;*/
    int lowcap;

    int c1,c2;           /* chromosomes being examined */
    int sc1, sc2;        /* signed versions, in case need to flip them */
    int f1;              /* is chromo c1 flipped? */

    int s1=0,e1=0,s2=0,e2=0;  /* start & end of chromosomes being examined */
    int s1i, e1i;        /* start/end of interior, dropping caps if multichrom */
    int s2_orig=0, e2_orig=0;  /* original start/end of c2, in case it's moved */
    int move_c2;         /* boolean, indicating if c2 is moved to abut c1 */
	
	/* GB: new for heuristic 5 and to look for reversals of fixed size */
	int start_j, end_j;
	int first_loop = TRUE;
	

    /* FUNCTIONALITY MOVED to init_isBP_array */
    /* test at breakpoints only, or test at all locations? */
    /*  int testbp = TRUE;  */
    /*  int testbp = FALSE; */

    cbounds_t cb;        /* structure with chromosome boundaries */
//	cb.cNum = NULL;
//	cb.cBound = NULL;

    tryrevparams_t tryrevparams;

    char *breakpoints;    /* bit array of genes that are breakpoints */

    //int *successors;     /* array of successors in final genome */

    // GB not necessary anymore, we garantee that genomes have an extra null chromo
    /* we pad the genomes with a null chromosome to allow fissions */
    //  int num_genes_p = num_genes + 1;                /* # genes, incl. pad */
    //  int num_chromosomes_p = num_chromosomes + 1;    /* # chromos, incl. pad */
    //  int num_chromosomes_nn;                      /* # nonnull chromos in g1 */

    // GB
    int genomenb;
    int gindex1;		// the index of the current g1 in the array
                  //	int gindex2;		// the index of the current g2 in the array
    int merge_with;		// if current rearrangement implies a merge
    int total_reduction; // when we try a rearrangement, this is the
								// total reduction of the distances
	
	cb.cNum = (int *) NULL;
	cb.cBound = (int *) NULL;

	changes = (int *)e_malloc(nb_spec*sizeof(int), "changes in my_tesrev.c");
	
    // GB: this is the equivalent of init_scenario_mem in testrev.c
    
    alloc_simple_genome(genome_list[0].num_g, &trial_genome);

    // GB: allocate memory for breakpoints, the initialization is done for each genome later

    init_isBP_mem(genome_list[0].num_g, genome_list[0].num_chr,
                  &breakpoints);
	
	// computes the closest pair of genomes
//	find_closest_pair(Genomes, &pair1, &pair2, nb_spec);

//	fprintf(stdout, "reduction type: %d %d %d\n", reduction_type, nb_spec, spec_left);
//	special_print_genomes2(stdout, Genomes, nb_spec);
	//fprintf(stdout, "only look for rearrangements in closest pair: %d %d\n", pair1, pair2);

	// GB first compute the pairwise distances
	//compute_dist_mat(Genomes, nb_spec, num_genes, num_chromosomes, distmem);

	while (first_loop == TRUE || (the_list == NULL && rev_size >= 0 && rev_size <= Genomes->max_chromo_size)) {
		//fprintf(stdout, "in while loop, rev_size == %d\n", rev_size);
		first_loop = FALSE; 
		/* in most case we will only do this loop once, except when we are looking for reversals
		of a fixed size and we want to increase that size */
		// GB we now want to loop over all possible g1 and pick g2 to be any other genome
		for (genomenb = 0; genomenb<nb_spec; genomenb++) {
			if (rev_size < 0) {
				gindex1 = (last_genomenb+genomenb)%nb_spec; // we alternate the genome we look at first
															// we actually start from the last_genomenb
															// because it's a FirstInLastOut queue 
			}
			else {
				// if rev_size >= 0 we will exit as soon as we find something good
				// so don't start with last_genomenb
				gindex1 = (last_genomenb+genomenb+1)%nb_spec;
			}
			if (pair1==-1 || gindex1 == pair1 || gindex1 == pair2) {
				if (Genomes->same_as[gindex1] == -1 && (optimize==FALSE || gindex1 == 3)) { // otherwise the genome has been merged already
				
				//   if (Genomes->same_as[gindex1] == -1) { // otherwise the genome has been merged already
				
				
				//fprintf(stdout, "last_genomenb %d so now gindex1 %d\n", last_genomenb, gindex1);
				
				/*		if (gindex1 == 0) // GB the second genome can be anything but gindex1
				gindex2 = 1;
				else gindex2 = 0;
				*/		
				// initialization from old procedure
				s1=e1=s2=e2=0;
				s2_orig=e2_orig=0;
				// GB: Fonctionnality now is isBP
				//compute_successor(num_genes, &genome_list[gindex2], &successors);
				//  reformat_caps(&genome_list[gindex1], num_genes, num_chromosomes);
				//  reformat_caps(&genome_list[gindex2], num_genes, num_chromosomes);
				if (genome_list[gindex1].num_chr > 0) {
					/* get the chromosome boundaries */
					init_cbounds(genome_list[gindex1].num_g,
								 genome_list[gindex1].num_chr,
								 &genome_list[gindex1],
								 &cb);
				}
				/* initialize breakpoint array */
				/* functionality separated out and changed from init_scenario_mem above */
				//		init_isBP_mem(num_genes,num_chromosomes,
				//			&breakpoints);
				init_isBP_array(genome_list[gindex1].num_g,
								genome_list[gindex1].num_chr, nb_spec,
								Genomes, gindex1,
								breakpoints);
				/* the lowest cap number */
				lowcap = genome_list[gindex1].num_g - 2*genome_list[gindex1].num_chr + 1; 
				/* initialize parameters used for every distance trial */
				tryrevparams.trial_genome = &trial_genome;
				// probably don't need next line anymore
				//tryrevparams.dest_genome = &genome_list[gindex2];
				tryrevparams.num_genes = genome_list[gindex1].num_g;
				tryrevparams.num_chromosomes = genome_list[gindex1].num_chr;
				tryrevparams.distmem = distmem;
				/* categorize by operation type */
				/************************************************************************/
				/*               Test all 1 chromo operations                           */
				/************************************************************************/
				/* start and end of the first null chromosome;
				used in fission operations */
				if (genome_list[gindex1].num_chr>0) {
					s2 = cb.cBound[Genomes->nb_chromo[gindex1]+1 -1];
					e2 = cb.cBound[Genomes->nb_chromo[gindex1]+1] -1;
				}
				for (c1 = 1; c1 <= Genomes->nb_chromo[gindex1]; c1++) {
					/* start & end of chromosome */
					if (genome_list[gindex1].num_chr > 0) {
						s1 = cb.cBound[c1-1];
						e1 = cb.cBound[c1]-1;
					} else {
						s1 = 0;
						e1 = genome_list[gindex1].num_g-1;
					}
					/**************/
					/*  fissions  */
					/**************/
					if (fission) {
						if (genome_list[gindex1].num_chr>0) {
							/* for fissions, need to use padded # genes, chromos */
							//tryrevparams.num_genes = num_genes_p;
							//tryrevparams.num_chromosomes = num_chromosomes_p; 
							//tryrevparams.num_genes = num_genes;
							//tryrevparams.num_chromosomes = num_chromosomes; 
							foundany_i = FALSE; 
							for (i=s1+2; i<e1 ; i++)  {
								/* test the fission starting just before position i */
								if (!isBP(genome_list[gindex1].genes[i], breakpoints))
									continue;
								copy_and_reverse_genes(genome_list[gindex1].genes, trial_genome.genes, i,s2, genome_list[gindex1].num_g);
								//i,s2, num_genes_p);
								total_reduction = try_good_operation(CFIS, Genomes,nb_spec,spec_left,gindex1,&tryrevparams, 
																	 &merge_with, reduction_type, changes, pair1, pair2);
								if (total_reduction >= reduction_type) {														score_list+=total_reduction;
								//	fprintf(stdout, "socre_list += %d", score_list);fflush(stdout);
                                    if (!only_size) 
										the_list = list_reag_insert(the_list, gindex1, CFIS, 
																	i, s2, 0, 0, merge_with, total_reduction, changes, nb_spec);
									//print_genome_multichrom(&trial_genome,
									//			      num_genes_p, num_chromosomes_p);
								}
					}   
				}
			}
					
					
					/************************/
					/*  internal reversals  */
					/************************/
					if (reversals) {
						
						/* for reversals, need to use unpadded # genes, chromos */
						//tryrevparams.num_genes = num_genes;
						//tryrevparams.num_chromosomes = num_chromosomes;
						
						/* determine start/end of interior of chromosome, by
						removing caps if multichromosomal */
						if (genome_list[gindex1].num_chr==0) {
							s1i = s1; e1i = e1;
						} else {
							s1i = s1+1; e1i = e1-1;
						}
						
						for (i=s1i ; i <= e1i ; i++) {
							foundany_i = FALSE;
							
							/* only start at breakpoints */
							
							if (!isBP(genome_list[gindex1].genes[i], breakpoints))
								continue;
							
							// for heuristic 5, only look for reversals of size rev_size and stop when you see one
							// we don't even look for longer ones
	
							if (rev_size >= 0) {
								if (i+rev_size <= e1i) {
									// if this reversal is legal
									start_j = end_j = i+rev_size;
								}
								else {
									start_j = 1;
									end_j = 0;
								}
							}
							else {
								start_j = i;
								end_j =e1i;
							}
							
							for (j=start_j ; j<=end_j ; j++) {
								//					for (j=i ; j<=e1i ; j++) {
								
								if (i == s1i && j == e1i && genome_list[gindex1].num_chr>0)
									continue; /* no chromosome flips in multichromosomal mode */
								
								/* only end at breakpoints */
								if (!isBP(-genome_list[gindex1].genes[j], breakpoints))
									continue;
								
								copy_and_reverse_genes(genome_list[gindex1].genes,
													   trial_genome.genes,
													   i,j, genome_list[gindex1].num_g);
								
								total_reduction = try_good_operation(CREV, Genomes, nb_spec, spec_left,gindex1,&tryrevparams, 
																	 &merge_with, reduction_type, changes, pair1, pair2);
								
								if (total_reduction >= reduction_type) {
									score_list+=total_reduction;
									//fprintf(stdout, "socre_list += %d", score_list);
                                    //fflush(stdout);
                                    if (!only_size) 
										the_list = list_reag_insert(the_list, gindex1, CREV,
																	i, j, 0, 0, merge_with, total_reduction, changes, nb_spec);
									
									// for heuristic 5, if rev_size is positive, stop when you a good reversal
									if (rev_size >= 0) {
										// we want to exit
										//j = e1i+1;
										i = e1i+1;
										c1 = Genomes->nb_chromo[gindex1]+1;
										genomenb = nb_spec;
										
									}
									//print_genome_multichrom(&trial_genome, num_genes, num_chromosomes);
								}
								
								
								} /* end for j (ending point within chromosome) */
							}
					
				

						}     /* end for i (starting point within chromosome) */
			
			
					}  /* end for c1 */
  
  
  
  
  /************************************************************************/
  /*               Test all 2 chromo operations                           */
  /************************************************************************/
  
  if (genome_list[gindex1].num_chr>0) {
	  
	  
	  /* for fusions and translocations, need to use unpadded # genes, chromos */
	  //tryrevparams.num_genes = num_genes;
	  //tryrevparams.num_chromosomes = num_chromosomes;
	  
	  for (c1=1 ; c1 < Genomes->nb_chromo[gindex1] ; c1++) {
		  for (f1=0 ; f1<=1 ; f1++) { /* do we flip chromo c1 or not? */
			  /* start & end of chromosome */
			  s1 = cb.cBound[c1-1];
			  e1 = cb.cBound[c1]-1;
			  
			  
			  if (f1) {
				  /* when we flip this chromosome, we will never be looking
				  at any earlier chromosomes again, so there is no
				  need to restore it later */
				  
				  reverse_in_place_genes(genome_list[gindex1].genes, s1,e1);
			  }
			  
			  
			  for (c2=c1+1 ; c2 <= Genomes->nb_chromo[gindex1] ; c2++) {
				  
				  //printf("chromos: %d %d\n", c1, c2);
				  
				  /* identify chromosomes involved */
				  /* - indicates chromo is flipped */
				  sc1 = f1 ? -c1 : c1;
				  sc2 = c2;
				  
				  /* start & end of chromosome, pointing at caps */
				  s2 = cb.cBound[c2-1];
				  e2 = cb.cBound[c2]-1;
				  
				  /* speedup: move chromo c2 to just after c1 */
				  move_c2 = (e1+1 < s2);
				  if (move_c2) {
					  s2_orig = s2;
					  e2_orig = e2;
					  
					  /* from:  xxx 1 A 2 yyy 3 B 4 zzz */
					  /* to:    xxx 1 A 2 -4 -B -3 -yyy zzz */
					  reverse_in_place_genes(genome_list[gindex1].genes,
											 e1+1,e2);
					  
					  /* new location */
					  s2 = e1+1;
					  e2 = s2+(e2_orig-s2_orig);
					  
					  /* indicate chromo2 was flipped */
					  sc2 = -sc2;
				  }
				  
				  
				  
				  
				  
				  
				  /*******************************/
				  /*  fusions & translocations   */
				  /*******************************/
				  
				  if (fusion || transloc) {
					  for (i=s1+1 ; i <= e1 ; i++) {
						  foundany_i = FALSE;
						  
						  /* only start at breakpoints */
						  if (!isBP(genome_list[gindex1].genes[i], breakpoints))
							  continue;
						  
						  for (j=s2 ; j<e2 ; j++) {
							  
							  /* must take a nonnull portion of at least one chromosome;
							  otherwise it's just a cap flip */
							  if (i==e1 && j==s2)
								  continue;
							  
							  /* also disallow exchanging both chromos */
							  if (i==s1+1 && j==e2-1)
								  continue;
							  
							  /* only end at breakpoints */
							  if (!isBP(-genome_list[gindex1].genes[j], breakpoints))
								  continue;
							  
							  
							  copy_and_reverse_genes(genome_list[gindex1].genes,
													 trial_genome.genes,
													 i,j, 
													 genome_list[gindex1].num_g);
							  
							  optype = CTRA;
							  /* due to the null bug, if the above gave a fusion, the
								  caps must be reformatted for things to work correctly */
							  if ((i == s1+1 && j == s2) ||    /* fusion -chrom1 chrom2 */
								  (i == e1   && j == e2-1)) {  /* fusion chrom1 -chrom2 */
								  optype = CFUS;
								  reformat_caps(&trial_genome,
												genome_list[gindex1].num_g,
												genome_list[gindex1].num_chr);
							  }
							  if ((optype == CTRA && transloc) || (optype == CFUS && fusion)) {
								  
								  total_reduction = try_good_operation(optype, Genomes,nb_spec,spec_left,gindex1,&tryrevparams, 
																	   &merge_with, reduction_type, changes, pair1, pair2);
								  
								  if (total_reduction>=reduction_type) {
									  score_list+=total_reduction;
									  //fprintf(stdout, "socre_list += %d", score_list);fflush(stdout);
                                      if (!only_size)  {
										  
										  //print_genome_multichrom(&genome_list[gindex], num_genes, num_chromosomes);
										  //fprintf(stdout, "i %d j %d sc1 %d sc2 %d\n", i, j, sc1, sc2);
										  if (move_c2)
											  the_list = list_reag_insert(the_list, gindex1, optype,
																		  i, j + (s2_orig-s2), sc1, sc2, merge_with, total_reduction, changes, nb_spec);
										  else
											  the_list = list_reag_insert(the_list, gindex1, optype,
																		  i, j, sc1, sc2, merge_with, total_reduction, changes, nb_spec);
										  //print_genome_multichrom(&trial_genome, num_genes, num_chromosomes);
									  }
									  
								  }
							  }
						  } /* end for j (ending point within chromosome c2) */
						  
					  } /* end for i (starting point within chromosome c1) */
				  }
				  
				  /* if moved c2 for speed, move it back */
				  if (move_c2) {
					  /* from:  xxx 1 A 2 -4 -B -3 -yyy zzz */
					  /* to:    xxx 1 A 2 yyy 3 B 4 zzz */
					  reverse_in_place_genes(genome_list[gindex1].genes,
											 e1+1,e2_orig);
					  
				  }
				  
			  } /* end: c2 */
	if (f1) {
		/* acually better to put them back... GB */
		
		reverse_in_place_genes(genome_list[gindex1].genes, s1,e1);
	}
		  } /* end: f1 */
	  } /* end: c1 */
    
	
  } /* end: if (num_chromosomes>0) */
  
  /* clean up memory */
  
  if (genome_list[gindex1].num_chr>0) {
	  
	  free_cbounds(&cb);
	  //  free_simple_genome(dest_genome);
  }
  
		} /* else gindex1 is not same_as -1 */
		
  } /* GB end: gindex1 */

} // if pair1 or pair2

 if (rev_size >= 0) {
	 if (the_list == NULL) {
		 // couldn't find a good reversal of size rev_size, try increasing it
		 rev_size++;
		 //fprintf(stdout, "rev_size is now %d\n", rev_size);
	 }
	 else {
		 if (verbose) {
			 fprintf(stdout, "found a reversal of size %d\n", rev_size);
		 }
	 }
 }

  } /* GB end: while */

  clean_isBP_mem(breakpoints);
  free(trial_genome.genes);
  free(changes);
  
  /********************************************************/
  
  /* print summary */
  //print_tryops_counts(&tryrevparams);

  
  *a_list = the_list;
  if(verbose) {fprintf(stdout, "return score list = %d\n", score_list); fflush(stdout);}
  return score_list;
  
}			   			    
 



/*
int build_list_reag(list_listreag * in_list,//list_reag **a_list,
        G_struct *Genomes , int nb_spec, int spec_left,
					int reversals,            // TRUE if we want to include reverals
					int transloc, 	// TRUE if we want to include translocations
					int fusion,
					int fission, 	// TRUE if we want to include fusions and fissions
					mgr_distmem_t *distmem,
					int only_size, 			// don't build the list, just return it's size
					int optimize,   // we are optimizing the median
					int reduction_type,	// we only accept rearrangements reducing the total distance by at least
					* this amount ( = spec_left-1 for good rearrangements)
					int rev_size,  // for heurisitic 5, only look for reversal of this size
					int last_genomenb,//the last genome in which we did a rearrangement, want to alternate how we start the list 
					int pair1, int pair2,  //usually both -1 but use them if we want to inforce a merge
					int verbose
                    )
{
   // fprintf(stdout, "build_list_reag,last_genomenb=%d,check nb_chromo:%d,%d,%d\n",last_genomenb,Genomes->nb_chromo[0],Genomes->nb_chromo[1],Genomes->nb_chromo[2]);
    int tmp1, tmp2;
    struct mgr_genome_struct *genome_list = Genomes->genome_list;
    //	int *dist_mat = Genomes->dist_mat;		  // GB: the matrix with the pairwise distances 
    struct mgr_genome_struct * trial_genome; // struct mgr_genome_struct trial_genome;
    //list_reag *the_list = NULL;
    //the_list = *a_list; 
    list_listreag * the_list = in_list;
    int score_list = 0;
	int *changes; // the changes associated with a particular rearrangement
    int optype;
    //	int num_bkpt;
    int foundany_i;      // has a reversal been found for this start? 
    int i,j;
    // GB: d is remove because we use dist_mat in Genomes
        int d;
    int lowcap;
    int c1,c2;           // chromosomes being examined 
    int sc1, sc2;        // signed versions, in case need to flip them 
    int f1;              // is chromo c1 flipped? 
    int s1=0,e1=0,s2=0,e2=0;  // start & end of chromosomes being examined 
    int s1i, e1i;        // start/end of interior, dropping caps if multichrom 
    int s2_orig=0, e2_orig=0;  // original start/end of c2, in case it's moved 
    int move_c2;         // boolean, indicating if c2 is moved to abut c1 
	// GB: new for heuristic 5 and to look for reversals of fixed size 
	int start_j, end_j;
	int first_loop = TRUE;
// FUNCTIONALITY MOVED to init_isBP_array 
    // test at breakpoints only, or test at all locations? 
    //  int testbp = TRUE;  
    //  int testbp = FALSE; 
//    cbounds_t  cb;        // structure with chromosome boundaries 
//	cb.cNum = NULL;	cb.cBound = NULL;
    tryrevparams_t tryrevparams;
    char *breakpoints;    // bit array of genes that are breakpoints 
    //int *successors;     // array of successors in final genome 
    // GB not necessary anymore, we garantee that genomes have an extra null chromo
    // we pad the genomes with a null chromosome to allow fissions 
    //  int num_genes_p = num_genes + 1;                // # genes, incl. pad 
    //  int num_chromosomes_p = num_chromosomes + 1;    // # chromos, incl. pad 
    //  int num_chromosomes_nn;                      // # nonnull chromos in g1 
    // GB
    int genomenb;
    int gindex1;		// the index of the current g1 in the array
                  //	int gindex2;		// the index of the current g2 in the array
    int merge_with;		// if current rearrangement implies a merge
    int total_reduction; // when we try a rearrangement, this is the
								// total reduction of the distances

    cbounds_t * cb = &cb_build_list_reag;
    //move this out later
	changes = (int *)e_malloc(nb_spec*sizeof(int), "changes in my_tesrev.c");
  
    // GB: this is the equivalent of init_scenario_mem in testrev.c
    trial_genome = &trial_genome_build_list_reag;
    if( (trial_genome == NULL) ||(trial_genome->genes == NULL) )
    {
        if(verbose) {
        fprintf(stdout, "init memory for trial_genome_build_list_reag,num_g=%d\n",
                genome_list[0].num_g);
        fflush(stdout);}
        alloc_simple_genome(genome_list[0].num_g, trial_genome);
    }
    else {}
    
    // GB: allocate memory for breakpoints, the initialization is done for each genome later
    Breakpoint_array * bkarr = &bkarr_build_list_reag;
    if((bkarr->array == NULL)||(bkarr->num_genes==0))
    {
        if(verbose) fprintf(stdout, "init bkarr_build_list_reag\n");
        int num_g = genome_list[0].num_g;
        bkarr->num_genes = num_g;
        bkarr->array =  (char *) e_malloc((num_g/4 + 1) * sizeof(char),
									  "breakpoints");
    }
    else { //since all gene array have the same size, we don't need to worry about (bkarr-?num_genes =?= num_).  
    }
    breakpoints = (bkarr->array);
   
	// GB first compute the pairwise distances
	//compute_dist_mat(Genomes, nb_spec, num_genes, num_chromosomes, distmem);
 //   fprintf(stdout, "start while loop\n"); fflush(stdout);
	while (first_loop == TRUE || ( (the_list->list_size == 0)//the_list == NULL 
                && rev_size >= 0 && rev_size <= Genomes->max_chromo_size)) {
     //   	fprintf(stdout, "in while loop, rev_size == %d\n", rev_size);
		first_loop = FALSE; 
		// in most case we will only do this loop once, except when we are looking for reversals of a fixed size and we want to increase that size 
		// GB we now want to loop over all possible g1 and pick g2 to be any other genome
		for (genomenb = 0; genomenb<nb_spec; genomenb++) {
			if (rev_size < 0) {
				gindex1 = (last_genomenb+genomenb)%nb_spec; 
                // we alternate the genome we look at first
				// we actually start from the last_genomenb
				// because it's a FirstInLastOut queue
			}
			else {
				// if rev_size >= 0 we will exit as soon as we find something good
				// so don't start with last_genomenb
				gindex1 = (last_genomenb+genomenb+1)%nb_spec;
			}		
			if (pair1==-1 || gindex1 == pair1 || gindex1 == pair2) {
				if (Genomes->same_as[gindex1] == -1 && (optimize==FALSE || gindex1 == 3))
                { 
                    // otherwise the genome has been merged already
				   //fprintf(stdout, "last_genomenb %d so now gindex1 %d\n", 
                   //last_genomenb, gindex1);
				// initialization from old procedure
				s1=e1=s2=e2=0;
				s2_orig=e2_orig=0;
				// GB: Fonctionnality now is isBP
				//compute_successor(num_genes, &genome_list[gindex2], &successors);
				//  reformat_caps(&genome_list[gindex1], num_genes, num_chromosomes);
				//  reformat_caps(&genome_list[gindex2], num_genes, num_chromosomes);
				if (genome_list[gindex1].num_chr > 0) {
					// get the chromosome boundaries 
                    if((cb->cNum==NULL)||(cb->cBound==NULL))
                    {
                        if(verbose) {
                        fprintf(stdout, "init cb_build_list_reag\n"); fflush(stdout); }
                        int num_g    = genome_list[gindex1].num_g;
                        int num_chr  = genome_list[gindex1].num_chr;
                        cb->cNum     = (int *) malloc(num_g*sizeof(int));
                        cb->cBound   = (int *) malloc((num_chr+1)*sizeof(int));
                        cb->num_genes = num_g;
                        cb->num_chromosomes = num_chr;
                    } else {}
					init_cbounds_wmem(genome_list[gindex1].num_g,
								 genome_list[gindex1].num_chr,
								 &genome_list[gindex1],
								 cb);
				}
				// initialize breakpoint array 
				// functionality separated out and changed from init_scenario_mem above 
				//		init_isBP_mem(num_genes,num_chromosomes,
				//			&breakpoints);
				init_isBP_array(genome_list[gindex1].num_g,
								genome_list[gindex1].num_chr, nb_spec,
								Genomes, gindex1,
								breakpoints);
				// the lowest cap number 
				lowcap = genome_list[gindex1].num_g - 2*genome_list[gindex1].num_chr + 1; 
				// initialize parameters used for every distance trial 
				tryrevparams.trial_genome = &trial_genome_build_list_reag;//trial_genome;
                tryrevparams.dest_genome = NULL;
				// probably don't need next line anymore
				//tryrevparams.dest_genome = &genome_list[gindex2];
				tryrevparams.num_genes = genome_list[gindex1].num_g;
				tryrevparams.num_chromosomes = genome_list[gindex1].num_chr;
				tryrevparams.distmem = distmem;
				// categorize by operation type 
				//               Test all 1 chromo operations                           
				// start and end of the first null chromosome;	used in fission operations 
				if (genome_list[gindex1].num_chr>0) {
                    s2 = cb->cBound[Genomes->nb_chromo[gindex1]+1 -1];
					e2 = cb->cBound[Genomes->nb_chromo[gindex1]+1] -1;
				}
             //   fprintf(stdout, "Test all 1 chromo oper, for c1=1 to %d\n",
              //          Genomes->nb_chromo[gindex1]);
                fflush(stdout);
				for (c1 = 1; c1 <= Genomes->nb_chromo[gindex1]; c1++) {
					// start & end of chromosome 
					if (genome_list[gindex1].num_chr > 0) {
						s1 = cb->cBound[c1-1];
						e1 = cb->cBound[c1]-1;
					} else {
						s1 = 0;
						e1 = genome_list[gindex1].num_g-1;
					}
					if (fission) {
						if (genome_list[gindex1].num_chr>0) {
							//tryrevparams.num_genes = num_genes_p;
							//tryrevparams.num_chromosomes = num_chromosomes_p; 
							//tryrevparams.num_genes = num_genes;
							//tryrevparams.num_chromosomes = num_chromosomes; 
							foundany_i = FALSE; 
							for (i=s1+2; i<e1 ; i++)  {
								if (!isBP(genome_list[gindex1].genes[i], breakpoints))
									continue;
								copy_and_reverse_genes
                                    (genome_list[gindex1].genes, trial_genome->genes,
                                     i,s2, genome_list[gindex1].num_g);
                                total_reduction = try_good_operation
                                    (CFIS, Genomes,nb_spec,spec_left,gindex1,&tryrevparams, 
									 &merge_with, reduction_type, changes, pair1, pair2);
								if (total_reduction >= reduction_type) {
									score_list+=total_reduction;
									if (!only_size)
                                    {
                                        list_reag_insert2
                                          (the_list, gindex1, CFIS, 
										 i, s2, 0, 0, merge_with, total_reduction,
                                         changes, nb_spec);
                                   //     fprintf(stdout, "fission, %d>=%d,reaglist size++ = %d\n",total_reduction, reduction_type,the_list->list_size); fflush(stdout); 
                                    }
									//print_genome_multichrom(&trial_genome,
									//	num_genes_p, num_chromosomes_p);
								}
                            }//for (i=s1+2; i<e1 ; i++)   
                        }//if (genome_list[gindex1].num_chr>0)
                    }//if (fission)
					if (reversals) {
						//tryrevparams.num_genes = num_genes;
						//tryrevparams.num_chromosomes = num_chromosomes;
						if (genome_list[gindex1].num_chr==0) {
							s1i = s1; e1i = e1;
						} else {
							s1i = s1+1; e1i = e1-1;
						}
                    //    fprintf(stdout, "for i=%d, i< %d, i++ \n",s1i, e1i); fflush(stdout);
						for (i=s1i ; i <= e1i ; i++) {
							foundany_i = FALSE;
							if (!isBP(genome_list[gindex1].genes[i], breakpoints))
								continue;
							// for heuristic 5, only look for reversals of size rev_size and stop when you see one
							// we don't even look for longer ones
							if (rev_size >= 0) {
								if (i+rev_size <= e1i) { // if this reversal is legal
									start_j = end_j = i+rev_size;
								}
								else {	start_j = 1;  end_j = 0; }
							}
							else { 
                                start_j = i;	end_j =e1i;
                            }
                       //     fprintf(stdout, "for j=%d, j<=%d,j++\n",start_j, end_j);
                        //    fflush(stdout);
							for (j=start_j ; j<=end_j ; j++) {
								if (i == s1i && j == e1i && genome_list[gindex1].num_chr>0)
									continue; 
								if (!isBP(-genome_list[gindex1].genes[j], breakpoints))
									continue;
                                 copy_and_reverse_genes(genome_list[gindex1].genes,
													   trial_genome->genes,
													   i,j, genome_list[gindex1].num_g);
							   total_reduction = try_good_operation
                                    (CREV, Genomes, nb_spec, spec_left,gindex1,&tryrevparams, 
									 &merge_with, reduction_type, changes, pair1, pair2);
								if (total_reduction >= reduction_type) {
									score_list+=total_reduction;
									if (!only_size)
                                    {
                                      list_reag_insert2
                                          (the_list, gindex1, CREV,i, j, 0, 0,
                                             merge_with, total_reduction, changes, nb_spec) ;
                                    //    fprintf(stdout, "interal rev, %d>=%d,reaglist size++ = %d\n",total_reduction, reduction_type, the_list->list_size); fflush(stdout);
                                    }
									// for heuristic 5, if rev_size is positive, stop when you a good reversal
									if (rev_size >= 0) {
										// we want to exit
										//j = e1i+1;
										i = e1i+1;
										c1 = Genomes->nb_chromo[gindex1]+1;
										genomenb = nb_spec;
									}
									//print_genome_multichrom(&trial_genome, num_genes, num_chromosomes);
								}
							} // end for j (ending point within chromosome) 
						}//end of for (i=s1i ; i <= e1i ; i++)
					} //end of if (reversals)     
				}  // end for c1 
                //end of for (c1 = 1; c1 <= Genomes->nb_chromo[gindex1]; c1++) Test all 1 chromo operations
                if (genome_list[gindex1].num_chr>0) {
                  for (c1=1 ; c1 < Genomes->nb_chromo[gindex1] ; c1++) {
                  //    fprintf(stdout, "for f1 = 0 to 1\n"); fflush(stdout);
                      for (f1=0 ; f1<=1 ; f1++) { // do we flip chromo c1 or not? 
                          s1 = cb->cBound[c1-1];
                          e1 = cb->cBound[c1]-1;
                          if (f1) {
                              reverse_in_place_genes(genome_list[gindex1].genes, s1,e1);
                          }
                          for (c2=c1+1 ; c2 <= Genomes->nb_chromo[gindex1] ; c2++) {
                              //printf("chromos: %d %d\n", c1, c2);
                              sc1 = f1 ? -c1 : c1;
                              sc2 = c2;
                              s2 = cb->cBound[c2-1];
                              e2 = cb->cBound[c2]-1;
                              move_c2 = (e1+1 < s2);
                              if (move_c2) {
                                  s2_orig = s2;
                                  e2_orig = e2;
                                  reverse_in_place_genes(genome_list[gindex1].genes,
                                                         e1+1,e2);
                                  s2 = e1+1;
                                  e2 = s2+(e2_orig-s2_orig);
                                  sc2 = -sc2;
                              }
                              if (fusion || transloc) {
                              //    fprintf(stdout, "for i=%d,i<=%d,i++\n", s1+1, e1); fflush(stdout);
                                  for (i=s1+1 ; i <= e1 ; i++) {
                                      foundany_i = FALSE;
                                      if (!isBP(genome_list[gindex1].genes[i], breakpoints))
                                          continue;
                                 //     fprintf(stdout, "for j=%d, j<%d, j++\n",s2,e2); fflush(stdout);
                                      for (j=s2 ; j<e2 ; j++) {
                                          if (i==e1 && j==s2)
                                              continue;
                                          if (i==s1+1 && j==e2-1)
                                              continue;
                                          if (!isBP(-genome_list[gindex1].genes[j], breakpoints))
                                              continue;
                                          copy_and_reverse_genes(genome_list[gindex1].genes,
                                                                 trial_genome->genes,
                                                                 i,j, 
                                                                 genome_list[gindex1].num_g);
                                          optype = CTRA;
                                          if ((i == s1+1 && j == s2) ||   
                                              (i == e1   && j == e2-1)) {  
                                              optype = CFUS;
                                              reformat_caps(trial_genome,
                                                            genome_list[gindex1].num_g,
                                                            genome_list[gindex1].num_chr);
                                          }
                                          if ((optype == CTRA && transloc) || (optype == CFUS && fusion)) {
                                             total_reduction = try_good_operation
                                                  (optype, Genomes,nb_spec,spec_left,gindex1,&tryrevparams,
                                                   &merge_with, reduction_type, changes, pair1, pair2);
                                              if (total_reduction>=reduction_type) {
                                                  score_list+=total_reduction;
                                                  if (!only_size)  {
                                                      //print_genome_multichrom(&genome_list[gindex], num_genes, num_chromosomes);
                                                      //fprintf(stdout, "i %d j %d sc1 %d sc2 %d\n", i, j, sc1, sc2);
                                                      if (move_c2)    
                                                      {
                                                        list_reag_insert2(the_list, gindex1, optype,
                                                                    i, j + (s2_orig-s2), sc1, sc2, merge_with, 
                                                                    total_reduction, changes, nb_spec);
                                                      }
                                                      else
                                                      {
                                                        
                                                          list_reag_insert2(the_list, gindex1, optype,
                                                                                      i, j, sc1, sc2, merge_with,
                                                                                      total_reduction, changes, nb_spec);
                                                          
                                                      }
                                                      //fprintf (stdout, "fusions & translocations \n");
                                                      fflush(stdout);
                                                      //print_genome_multichrom(&trial_genome, num_genes, num_chromosomes);
                                                  }
                                              }
                                          }
                                      } // end for j (ending point within chromosome c2) 
                                  } // end for i (starting point within chromosome c1) 
                              }//end of if (fusion || transloc)
                              // if moved c2 for speed, move it back 
                              if (move_c2) {
                                  reverse_in_place_genes(genome_list[gindex1].genes, e1+1,e2_orig);
                              }
                          }  //end of  for (c2=c1+1 ;....)
                          if (f1) {// acually better to put them back... GB 
                            reverse_in_place_genes(genome_list[gindex1].genes, s1,e1);
                          }
                      } // end: f1 // end of for (f1=0 ; f1<=1 ; f1++) 
                  } // end: c1  // end of for (c1=1 ; c1 < ....)
              } // end: if (num_chromosomes>0) 
              //end of if (genome_list[gindex1].num_chr>0), Test all 2 chromo operations 
              // clean up memory 
              if (genome_list[gindex1].num_chr>0) {
                  //free_cbounds(&cb);
                  //  free_simple_genome(dest_genome);
              }
          } // else gindex1 is not same_as -1 
          //end of if (Genomes->same_as[gindex1] == -1 && (optimize==FALSE || gindex1 == 3))
       }  
       //end of if (pair1==-1 || gindex1 == pair1 || gindex1 == pair2)
     } 
     //end of for (genomenb = 0; genomenb<nb_spec; genomenb++)
     if (rev_size >= 0) {
         if (the_list->list_size == 0) {//if (the_list == NULL) {
             // couldn't find a good reversal of size rev_size, try increasing it
             rev_size++;
             //fprintf(stdout, "rev_size is now %d\n", rev_size);
         }
         else {
         //    if (verbose) 
               {
                 fprintf(stdout, "found a reversal of size %d\n", rev_size);
               }
         }
     }//end of if (rev_size >= 0)
  } // GB end: while 
     //fprintf(stdout, "end of while loop\n"); fflush(stdout);
  //end of while (first_loop == TRUE || (the_list == NULL .....)
//  clean_isBP_mem(breakpoints);
//  free(trial_genome.genes);
  free(changes);
  
  //print_tryops_counts(&tryrevparams);
 // *a_list = the_list;
  return score_list;
}			   			    
 
*/

/* try an operation
   return previous distance - new distance */   
int my_try_operation(int optype,
		  int d,
		  tryrevparams_t *p)
{
  int d2;
  d2 = mcdist_noncircular(p->trial_genome, p->dest_genome,
			  p->num_genes, p->num_chromosomes,
			  p->distmem,
			  NULL);
  return (d-d2);
}

/* try a given operation
 * return the total distance reduction on the tree */
int try_good_operation(int optype, G_struct *Genomes, int nb_spec, int spec_left,
					   int gindex1, tryrevparams_t *p, 
					   int *merge_with, int reduction_type, int *changes, int pair1, int pair2)
{
	int i;
	int dist;
	struct mgr_genome_struct *genome_list = Genomes->genome_list;
	int *dist_mat = Genomes->dist_mat;
	int a_reduction, total_reduction;
//	fprintf(stdout,"try good operation, gindex1=%d,nb_spec=%d; ",gindex1,nb_spec); 
	 *merge_with = -1;
	 total_reduction = 0;
    	 for (i=0; i<nb_spec; i++) {
	 	 if (i!=gindex1 && Genomes->same_as[i]==-1) {
	 	 	 p->dest_genome = &genome_list[i];
			 dist = dist_mat[gindex1*nb_spec+i];
			 //a_reduction = my_try_operation(optype, dist, p);
 			 int d2 = mcdist_noncircular(p->trial_genome, p->dest_genome,
			  p->num_genes, p->num_chromosomes,
			  p->distmem, NULL);
             a_reduction = dist - d2;
			 if (a_reduction != 1) { // this pairwise distance is not reduced 
//				  printf("reduction type %d spec_left-1 %d\n", reduction_type, spec_left-1);
				  
				 if (reduction_type == spec_left-1) {/* if we are looking for rearrangements
				                                      * reducing each pairwise distance (good) */
						
				  return 0; 
				 }
				 if (pair1 != -1 && ((gindex1==pair1 && i==pair2) || (gindex1==pair2 && i==pair1))) {
					 //we want to merge pair1 and pair2, this rearrangement doesn't do that...
					 return 0;
				 }
			 } else { // the pairwise distance is reduced
				 if (dist==1) {
					 /* this operation corresponds to a merge with i
					 * we assume that before this operation no two genomes are the same
					 * which was tested for before */
					 (*merge_with) = i;
				 }
			 }
			 total_reduction+=a_reduction;
			 changes[i] = a_reduction;
	 	 }	 
		 else {
			 changes[i] = 0;
		 }
	 }
	 return total_reduction;
}
	 
