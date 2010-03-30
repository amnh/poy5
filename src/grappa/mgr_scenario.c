/* scenario.c
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

/* Last modified on Wed Aug 2, 2006, by Glenn Tesler
*/

#include <stdio.h>
#include <stdlib.h>
#include "mgrstructs.h"
#include "mgr_uniinvdist.h"
#include "mgr_scenario.h"
#include "mgr_mcrdist.h"
#include "mgr_e_malloc.h"

#include <caml/fail.h>


#define ABS(n) (n >=0 ? n : (0-n))



/* construct optimal reversal scenario */
/* go from genome g1 (=gamma)  to  genome g2  (=pi) */


void init_cbounds(int num_genes, int num_chromosomes,
		  struct mgr_genome_struct *g1,
		  cbounds_t *cb);

void free_cbounds(cbounds_t *cb);


void init_cbounds_memory (int num_genes, int num_chromosomes, cbounds_t * cb);

/* initialize memory
   compute breakpoints

   pad=0: current_genome, trial_genome are the same # genes/chroms as
          given by num_genes, num_chromosomes
   pad=1: add space for 1 null chromosome
*/
void init_scenario_mem(int num_genes, int num_chromosomes, mgr_distmem_t *distmem,
		       struct mgr_genome_struct *current_genome,
		       struct mgr_genome_struct *trial_genome,
		       struct mgr_genome_struct *g1,
		       int pad)
{
  int num_genes_p = num_genes + 2*pad;
  int num_chromosomes_p = num_chromosomes + pad;

  mcdist_allocmem(num_genes_p,num_chromosomes_p,
		  distmem);

  /* allocate space to store some genomes */
  alloc_simple_genome(num_genes_p, current_genome);
  alloc_simple_genome(num_genes_p, trial_genome);


  /* copy g1=gamma into current_genome,
   i.e., copy the genes array */
  copy_genes(g1->genes,current_genome->genes,num_genes);
  if (pad) {
    /* add a null chromosome if neccessary */
    current_genome->genes[num_genes] = num_genes+1;
    current_genome->genes[num_genes+1] = num_genes+2;
  }
}

/* initialize SUCCESSOR array used for computing breakpoints */
/* TODO: could stick space for array into distmem, if it will be
 * called repeatedly */
/* init_scenario_bp: allocates memory, stores ptr in *sucessors
 * init_scenario_bp_wmem: *successors already set
 */
void init_scenario_bp_wmem(int num_genes, int num_chromosomes,
			   /* mgr_distmem_t *distmem, */
			   struct mgr_genome_struct *g2,
			   int pad,
			   int **successors)
{
  int g,h;  /* successive elements when computing breakpoints */
  int i;
  int num_genes_p = num_genes + 2*pad;
  /* int num_chromosomes_p = num_chromosomes + pad; */



#if 0
  fprintf(stdout,"init_scenario_bp: ");
  for (i=0; i<num_genes; i++)
    fprintf(stdout," %d", g2->genes[i]);
  fprintf(stdout,"\n");
#endif



  /* For unsigned permutations, we never have to break "long strips";
     for signed permutations, this means that if   A B   appears in
     the final genome, then if the current genome has  A B  or  -B -A,
     we should not break those up in a reversal.
     So this is the same as avoiding non-breakpoints.

     For   A B C    the successor of B is C, the successor of -B is -A.
     In multichromosomal case, also cannot start/end reversal at start/end of
     perm.

     successors[i-1] is the successor of i            for i=1,...,num_genes+2
     successors[i+num_genes+2-1]  is successor of -i  for i=1,...,num_genes+2
     (num_genes+2 is used because of the possible null pad, whether
      or not the pad option is TRUE)

     */
  //  *successors = (int *) e_malloc(2*(num_genes+2+1)*sizeof(int),
  //				 "successors");

  for (i=0; i<num_genes-1; i++) {
    g = g2->genes[i];
    h = g2->genes[i+1];

    ISUCCESSOR(g) = h;
    ISUCCESSOR(-h) = -g;
  }

  /* deal with null pad */
  g = g2->genes[num_genes-1];
  if (pad) {
    h = num_genes+1;
    ISUCCESSOR(g) = h; ISUCCESSOR(-h) = -g;

    g = h; h++;
    ISUCCESSOR(g) = h; ISUCCESSOR(-h) = -g;

    g = h;
  }
  h = g2->genes[0];

  /* also deal with start and end;
     however, this initialization for start and end could be
     discarded, since it is tested for separately. */
  ISUCCESSOR(g) = num_genes_p+1;
  ISUCCESSOR(-h) = -(num_genes_p+1);


#if 0
  for (i=1; i<=num_genes_p+1; i++)
    fprintf(stdout,"SUCC(%d)=%d\tSUCC(%d)=%d\n",
	    i,ISUCCESSOR(i), -i,ISUCCESSOR(-i));
  fprintf(stdout,"end init_scenario_bp\n");
#endif
}

void init_scenario_bp(int num_genes, int num_chromosomes,
		      /* mgr_distmem_t *distmem, */
		      struct mgr_genome_struct *g2,
		      int pad,
		      int **successors)
{
  *successors = (int *) e_malloc(2*(num_genes+2+1)*sizeof(int),
				 "successors");
  init_scenario_bp_wmem(num_genes, num_chromosomes,
			g2,
			pad,
			successors);
}




void print_scenario(struct mgr_genome_struct *g1, struct mgr_genome_struct *g2,
		    int num_genes, int num_chromosomes) {
  mgr_distmem_t distmem;
  struct mgr_genome_struct current_genome, trial_genome;
  struct mgr_genome_struct *dest_genome = g2;
  int d, d2, i, j, start, end, chr1, chr2, tp1, tp2, empty1, empty2;
  int startvalue, endvalue;
  int t = 1;
  int lowcap;
  int found;

  /* avoid some reversals based on cap properties */
  int i_is_ltail, i_is_lend, i_is_rtail, postype, badrev;

  int *successors;  /* array of successors in final genome */


  char *description = NULL;     /* description of the step */

  init_scenario_mem(num_genes,num_chromosomes,&distmem,
		    &current_genome,&trial_genome,
		    g1,
		    0); /* no pad */
  init_scenario_bp(num_genes,num_chromosomes,
		   /*&distmem,*/
		    g2,
		    0, /* no pad */
		    &successors);




  /**********************************************************************/
  /*                   print an optimal scenario                        */
  /**********************************************************************/

  fprintf(stdout,"Step 0: (Source)\n");
  print_genes(current_genome.genes,num_genes);
    
  /* UNIchromosomal reversal distance computation */
  d = mgr_invdist_noncircular(&current_genome, dest_genome,
			  0, /* offset */
			  num_genes, &distmem);
  /* the lowest cap number */
  lowcap = num_genes - 2*num_chromosomes + 1; 


  /* print out starting message and pi */

  while (d>0) {
    d2 = d;
    found = FALSE; start = end = 0;
    for (i=0 ; i<num_genes  &&  !found ; i++) {

      /* if i points to an ltail, then it can only end at an rtail
       * (chromosome flip)
       */
#if 0
      i_is_ltail =
	position_type(i,current_genome.genes, num_genes, lowcap) == T_LTAIL;
#endif
      postype = position_type(i,current_genome.genes, num_genes, lowcap);
      i_is_ltail =  postype == T_LTAIL;
      i_is_lend  =  postype == T_LEND || postype == T_LREND;
      i_is_rtail =  postype == T_RTAIL;


      /* don't split signed 2-strip at start of reversal */
      if (i==0) {
	/* if genomes start w/same gene, don't reverse it */
	if (current_genome.genes[i] == dest_genome->genes[i]) continue;

      } else {
	/* don't start here if it breaks adjacency */
	if (SUCCESSOR(current_genome.genes[i-1]) == current_genome.genes[i])
	  continue;
      }


  
      for (j=i ; j<num_genes && !found ; j++) {

	/* avoid doing the distance computation if emulating certain
	 * operations: cap exchange & illegal reversal starting/ending @ cap
	 */
	postype = position_type(j, current_genome.genes, num_genes, lowcap);
	badrev = 0;
	switch (postype) {
	case T_RTAIL:
	  badrev = !i_is_ltail;  /* chrom flip OK, illegal reversals BAD */
	  break;
	case T_LTAIL:
	  badrev = i_is_rtail;   /* cap exchange BAD */
	  break;
	case T_LREND:
	case T_REND:
	  badrev = i_is_lend;    /* cap exchange BAD */
	  break;
	default:  /* T_NORMAL, T_LEND, T_ERROR */
          break;
	}
	if (badrev) continue;

#if 0
	if ((position_type(j, current_genome.genes,
			  num_genes, lowcap) == T_RTAIL)
	    != i_is_ltail)
	  continue;
#endif

	/* don't split signed 2-strip at end of reversal */
	if (j==num_genes-1) {
	  /* if genomes end w/same gene, don't reverse it */
	  if (current_genome.genes[j] == dest_genome->genes[j]) continue;

	} else {
	  /* don't end here if it breaks adjacency */
	  if (SUCCESSOR(current_genome.genes[j]) == current_genome.genes[j+1])
	    continue;
	}

	
	copy_and_reverse_genes(current_genome.genes,trial_genome.genes,i,j, num_genes);
	
	d2 = mgr_invdist_noncircular(&trial_genome, dest_genome,
				 0, /* offset */
				 num_genes, &distmem);

 	/* error check:
	   make sure d2 = d, d-1, or d+1
	   if not, output bug report.
	   */
	if ( ((d2 - d) > 1) || ((d2 - d) < -1) ) {
	  fprintf(stdout,"\nERROR: Potential distance caculation error:\n");

	  fprintf(stdout,"Source:\n\t");
	  print_genes(current_genome.genes, num_genes);

	  fprintf(stdout,"Reversal from gene %d to gene %d:\n\t", i, j );
	  print_genes(trial_genome.genes, num_genes);
	};
	if (d2 < d)  {
	  found = TRUE; start  = i; end = j;
	}
      }
    }
    if (!found) {
      /* no proper reversal found, report failure */
      fprintf(stdout, "ERROR: reversals ended prematurely\n");
      break;
    };
    /* so now reversal from position i to j decreased distance */
    d = d2;


    /* TODO: speed these up */
    chr1 = chromosome_number(start, current_genome.genes,num_genes, lowcap);
    chr2 = chromosome_number(end, current_genome.genes,num_genes,lowcap);

    tp1 = position_type(start, current_genome.genes,num_genes, lowcap);
    tp2 = position_type(end, current_genome.genes,num_genes,lowcap);

#if 0
    empty1 = is_empty_chromosome(current_genome.genes,
				 chr1, num_genes, lowcap);
    empty2 = is_empty_chromosome(current_genome.genes,
				 chr2, num_genes, lowcap);
#endif
    empty1 = is_null_chromosome(start,
				current_genome.genes, num_genes, lowcap);
    empty2 = is_null_chromosome(end+1,
				current_genome.genes, num_genes, lowcap);

    startvalue = current_genome.genes[start];
    endvalue = current_genome.genes[end];


    /* print out description of what reversal mimics:
       reversal, fission, fusion, translocation,
       flip (possible but rare)
       cap exchange (impossible)
       illegal operation (impossible)
       */


    if (tp1 == T_LTAIL || tp2 == T_RTAIL) {
      if (tp1 == T_LTAIL  &&  tp2 == T_RTAIL) {
	description = "Trivial chromosome flip";
      } else {
	description = "Illegal reversal";
      }
    } else if ((tp1 == T_RTAIL && tp2 == T_LTAIL) ||
	       ((tp1 == T_LEND || tp1 == T_LREND) &&
		(tp2 == T_REND || tp2 == T_LREND))) {
      description = "Cap exchange";
    } else if (chr1 == chr2) {
      description = "Reversal";
    } else if (empty1 || empty2) {
      description = "Fission";
    } else if ((tp1 == T_RTAIL && (tp2 == T_REND || tp2 == T_LREND)) ||
	       ((tp1 == T_LEND || tp1 == T_LREND) && tp2 == T_LTAIL)) {
      description = "Fusion";
    } else {
      description = "Translocation";
    }


    fprintf(stdout, "Step %d: ", t);
  
    /*
    fprintf(stdout, "Reverse from gene %d [%d] in chromosome %d to gene %d [%d] in chromosome %d: ", start, startvalue, chr1, end, endvalue, chr2);
    */

    fprintf(stdout,
	    "Chrom. %d, gene %d [%d] through chrom. %d, gene %d [%d]: ",
	    chr1, start, startvalue, chr2, end, endvalue);

    fprintf(stdout, "%s", description);

    if (d == 0)
      fprintf(stdout, " (Destination)");

    fprintf(stdout, "\n");

    print_genes(trial_genome.genes,num_genes);


    /* copy the successful trial_perm to current_perm */
    copy_genes(trial_genome.genes,current_genome.genes,num_genes);


    t++;
  };
  
  /*
  fprintf(stdout,"Destination:\n");
  print_genes(dest_genome->genes,num_genes);
  */


  clean_scenario_bp_mem(successors);
  clean_scenario_mem(&current_genome,&trial_genome,&distmem);

}

void clean_scenario_bp_mem(int *successors)
{
  free(successors);
}

void clean_scenario_mem(struct mgr_genome_struct *current_genome,
			struct mgr_genome_struct *trial_genome,
			mgr_distmem_t *distmem)
{
  /* free memory */
  free(trial_genome->genes);
  free(current_genome->genes);

  mcdist_freemem(distmem);
}

void print_genes_nocr(int *source_perm, int num_genes)
{
  int k;
  for (k = 0; k< num_genes;k++) {
    fprintf(stdout," %4d",source_perm[k]);
  };
};

void print_genes(int *source_perm, int num_genes)
{
  print_genes_nocr(source_perm, num_genes);
  fprintf(stdout, "\n");
}

void copy_genes(int *source_perm, int *dest_perm, int num_genes)
{
 int k;
	
	for (k = 0; k < num_genes; k++) {
		dest_perm[k] = source_perm[k];
	};   
}


/* copy from source to dest while doing a reversal */
void copy_and_reverse_genes(int *source_perm, int *dest_perm,
			    int start, int end,
			    int num_genes)
{
  int k;
  for (k = 0; k< start; k++)
    dest_perm[k] = source_perm[k];
  for (k = start; k<= end; k++)
    dest_perm[k] = 0 - source_perm[end+start-k];
  for (k = end+1; k<num_genes; k++)
    dest_perm[k] = source_perm[k];
}



void reverse_in_place_genes(int *source_perm,
			    int start, int end)
{
  int k,l,temp;
  for (k=start,l=end;k<=l;k++,l--)
    {
      temp = source_perm[k];
      source_perm[k] = 0 - source_perm[l];

      /* note: if odd number of genes, the middle one
	 will have k=l, and the result above will simply
	 be clobbered by the result below. */

      source_perm[l] = 0 - temp; 
    };
}

int position_type(int pos, int *genes, int num_genes, int cap)
{
  if (pos<0 || pos >= num_genes) return T_ERROR;
  if (ABS(genes[pos]) < cap) {
    if ((pos>0 && ABS(genes[pos-1]) >= cap) &&
	(pos<num_genes-1 && ABS(genes[pos+1]) >= cap))
      return T_LREND;
    if (pos>0 && ABS(genes[pos-1]) >= cap)
      return T_LEND;
    if (pos<num_genes-1 && ABS(genes[pos+1]) >= cap)
      return T_REND;
    return T_NORMAL;
  };
  if ((genes[pos] > 0 && (genes[pos] - cap)%2 == 0) 
      || (genes[pos]  <0 && (genes[pos] - cap)%2 != 0)) return T_LTAIL;
  return T_RTAIL;
}

int chromosome_number(int pos, int *genes, int num_genes, int cap)
{
  int i;
  int count = 0;
  if (pos<0 || pos >= num_genes) return -1;
  for (i=0;i<=pos;i++)
    if (ABS(genes[i])>=cap) count++;
  return ((count+1)/2);
}

#if 0
/* not robust for error input */
int is_empty_chromosome(int *genes, int ch_num, int num_genes, int cap)
{
  int i; 
  int count = 0;
  int left,right; 
  left = right = 0;
  for (i = 0; i< num_genes; i++) {
    if (ABS(genes[i])>= cap) 
      {
	count++;
	if (count == 2*ch_num -1 ) left = i;
	if (count == 2*ch_num) 
	  {
	    right = i;
	    break;
	  };
      };
  };
  return (left+1==right);
}
#endif

/* is pos at the right cap of a null chromosome? */
int is_null_chromosome(int pos, int *genes, int num_genes, int cap)
{
  return
    position_type(pos, genes, num_genes, cap) == T_RTAIL &&
    position_type(pos-1, genes, num_genes, cap) == T_LTAIL;
}



/****************************************************************************/

/* construct a genome used for computations
   The important part is the gene array */
void alloc_simple_genome(int num_genes,
			 struct mgr_genome_struct *genome)
{
  genome->genome_num = 0;           /* UNUSED */
  genome->encoding   = NULL;        /* UNUSED */
  genome->gnamePtr   = NULL;        /* UNUSED */
  genome->genes      = (int *) e_malloc(num_genes*sizeof(int),
										"genes (alloc_simple_genome)");
  
  genome->num_g = 0;            /* number of genes (strip) in this genome */
  genome->num_chr = 0;          /* number of chromosome in this genome */
  genome->ori_num_g = 0;        /* original number of genes in this genome before alphabet conversion */
  genome->alphabet = NULL;	/* each gene can actually be a strip, alphabet stores conversion table */
  
}

void free_simple_genome(struct mgr_genome_struct *genome)
{
  free(genome->genes);
}

