/* write_data.c
*    Formatted output of genome permutations.
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

#include <stdio.h>
#include <stdlib.h>
#include "mgrstructs.h"

void print_genome_unichrom(struct mgr_genome_struct *g,
			   int num_genes)
{
  int i;

  fprintf(stdout,">%s\n", g->gnamePtr);

  for (i=0 ; i<num_genes ; i++)
    fprintf(stdout,"%d ",g->genes[i]);

  fprintf(stdout,"\n");
  fflush(stdout);
}



/* print genome on 1 line.  Caps suppressed.  Chromosomes delimited by $. */
void print_genome_multichrom(struct mgr_genome_struct *g,
			     int num_genes, int num_chromosomes)
{
  int i;
  int lowcap = num_genes - 2*num_chromosomes + 1;
  int inchrom = FALSE;
  int gene;

  if (g->gnamePtr != (char *)NULL) 
    fprintf(stdout,">%s\n", g->gnamePtr);

  for (i=0 ; i<num_genes ; i++) {
    gene = g->genes[i];
    if (gene >= lowcap || gene <= -lowcap) {
      inchrom = !inchrom;
      if (!inchrom && i < num_genes-1)
	fprintf(stdout, "$(%d) ", i);
    } else
      fprintf(stdout,"%d(%d) ",gene, i);
  }

  fprintf(stdout,"\n");
  fflush(stdout);
}

/* Print genome on multiple lines.  Caps included.
 * Chromosomes delimited by "/"
 */
void print_genome_multichrom_wcaps(struct mgr_genome_struct *g,
				   int num_genes, int num_chromosomes)
{
  int i;
  int lowcap = num_genes - 2*num_chromosomes + 1;
  int inchrom = FALSE;
  int gene;

  if (g->gnamePtr != (char *)NULL) 
    fprintf(stdout,">%s\n", g->gnamePtr);

  for (i=0 ; i<num_genes ; i++) {
    gene = g->genes[i];
    fprintf(stdout,"%d ",gene);

    if (gene >= lowcap || gene <= -lowcap) {
      inchrom = !inchrom;
      if (!inchrom && i < num_genes-1)
	fprintf(stdout, "/\n");
    }
  }

  fprintf(stdout,"/\n");
  fflush(stdout);
}


/* for printing columns of numbers in the range -num_genes..num_genes
 * with one space between them, compute the optimal width for 80 char screen.
 * Either w=3, 4, or 7 (as w+1 divides 80 for these)
 */
int optimal_gene_width(int num_genes)
{
  if (num_genes >= 1000) return 7;
  if (num_genes >= 100)  return 4;
  return 3;
}


/* Print multichromosomal genome.
 * width = -1: compute a width to use for all genes for alignment
 * width = 0: don't attempt to align, print each in min width
 * width > 0: attempt to right align each gene to this width
 *
 * multiline = 0: on one line
 * multiline = 1: one chromosome per line
 *
 * showcaps = 0: Delimit chromos by "$", omit caps
 * showcaps = 1: Delimit chromos by "/", show caps
 * showcaps = 2: Show caps
 *
 * showname = 0: don't show genome name
 * showname = 1: do show it
 */
void print_genome_multichrom_3(struct mgr_genome_struct *g,
			       int num_genes,
			       int num_chromosomes,
			       int width,
			       int multiline,
			       int showcaps,
			       int showname
			       )
{
  int i;
  int lowcap = num_genes - 2*num_chromosomes + 1;
  int inchrom = FALSE;
  int gene;

  if (width == -1) {
    width = optimal_gene_width(num_genes);
  }

  if (showname && g->gnamePtr != (char *)NULL) 
    fprintf(stdout,">%s\n", g->gnamePtr);

  for (i=0 ; i<num_genes ; i++) {
    gene = g->genes[i];

    if (gene >= lowcap || gene <= -lowcap) {
      if (showcaps > 0) {
	if (width == 0)
	  fprintf(stdout, "%d ",gene);
	else
	  fprintf(stdout, "%*d ", width, gene);
      }

      inchrom = !inchrom;
      if (!inchrom && i < num_genes-1) {
	switch(showcaps) {
	case 0: fprintf(stdout, "$"); break;
	case 1: fprintf(stdout, "/"); break;
	default: break;
	}
	if (multiline == 0) {
	  if (showcaps < 2) fprintf(stdout, " ");
	} else {
	  fprintf(stdout, "\n");
	}
      }
    } else {
      if (width == 0)
	fprintf(stdout, "%d ",gene);
      else
	fprintf(stdout, "%*d ", width, gene);
    }
  }

  switch(showcaps) {
  case 0: fprintf(stdout, "$"); break;
  case 1: fprintf(stdout, "/"); break;
  default: break;
  }
  fprintf(stdout,"\n");

  fflush(stdout);
}


