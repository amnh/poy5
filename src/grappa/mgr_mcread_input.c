/* mcread_input.c
*    Read input genomes.
*
* Copyright (C) 2001-2006 The Regents of the University of California
* by Glenn Tesler
*
* Contains code from read_input.h in GRAPPA 1.02
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

/* Replacement for GRAPPA 1.02: read_input.c */

/* New data file format:
a) Allow any whitespace & \n's w/in gene list
b) Allow comment lines beginning with "#"
c) Allow multiple chromosomes.  Symbols ";" or "$" are used to
separate chromosomes.  For neatness, they are often placed at
the end of a line, but the line break is not the chromosome delimeter.

>name1
2 -3 1 13 ; -4 6 -5 -8 9 ;
# comment
10 11 12 ; 7
>name2
6 -1 2 3 5 7 ;
4 12 13 9 -10 -11 8
>name3
1 -10 11 $
8 -9 -10 12 2 $
5 6 -3 4 7 $
13



d) Allow bracketed segments [...]
whose sign is not known;
symbols (){} will also be allowed but ignored
>name1
2 -3 1 13 ; -4 6 -5 -8 9 ;
# comment
10 11 12 ; 7
>name2
6 [ -1 2 ] 3 5 7 ;
4 12 13 [ 9 ] -10 -11 8
>name3
1 -10 11 $
8 [ -9 -10 ] 12 2 $
[ 5 6 ] -3 [ 4 7 ] $
13
*/

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include "mgrstructs.h"
#include "mgr_mcread_input.h"
/*#include "binencode.h"*/
/*#include "condense.h"*/
#ifdef MPBPA
#include "mpi.h"
#endif

#include "mgr_e_malloc.h"

void check_genes(int *genes, int num_genes,
		 int at_genome) {
  int i, g0, g1;

  int *checkArr;

  checkArr = (int *) e_malloc((num_genes+1)*sizeof(int), "checkArr");

  for (i=1 ; i<=num_genes ; i++)
    checkArr[i] = 0;

  for (i=0 ; i<num_genes ; i++) {
    g0 = genes[i];
    g1 = (g0 > 0 ? g0 : -g0);
    if (g1 == 0) {
      fflush(stdout);
      fclose(stdout);
      fprintf(stderr,"ERROR: Reading Input Genome %d\n", at_genome);
      fprintf(stderr,"ERROR: Input gene == 0\n");
      exit(-1);
    }
    if (g1 >num_genes) {
      fflush(stdout);
      fclose(stdout);
      fprintf(stderr,"ERROR: Reading Input Genome %d\n", at_genome);
      fprintf(stderr,"ERROR: Input gene (%d) larger than number of genes %d\n",
	      g0, num_genes);
      exit(-1);
    }
    if (checkArr[g1] != 0) {
      fflush(stdout);
      fclose(stdout);
      fprintf(stderr,"ERROR: Reading Input Genome %d\n", at_genome);
      fprintf(stderr,"ERROR: Duplicate gene (%d)\n", g1);
      exit(-1);
    }

    checkArr[g1] = 1;
  }

  for (i=1 ; i<=num_genes ; i++) {
    if (checkArr[i] != 1) {
      fflush(stdout);
      fclose(stdout);
      fprintf(stderr,"ERROR: Reading Input Genome %d\n", at_genome);
      fprintf(stderr,"ERROR: Missing gene (%d)\n", i);
      exit(-1);
    }
  }

  free(checkArr);
  return;
}


/**********************************************************************
 *                file i/o utility routines                           *
 **********************************************************************/

/* preview the next character */
int peek_char(FILE *fd) {
   int c;
   if ((c=getc(fd)) !=EOF)
      ungetc(c,fd);
   return c;
}

/* discard all input until end of line */
void discard_line(FILE *fd) {
  int c;
  int numchars = 0;

  while ((c=getc(fd)) != EOF && c != '\n') {
    numchars++;
    if (numchars >= MAX_STR_LEN) {
      fflush(stdout);
      fclose(stdout);
      fprintf(stderr,"ERROR: \tBuffer overflow reading input.\n");
      fprintf(stderr,"\tMCDIST is compiled with a maximum ");
      fprintf(stderr,"read buffer length of %d.\n", MAX_STR_LEN);
      fprintf(stderr,"\tPlease edit the header file \"mcstructs.h\", ");
      fprintf(stderr,"double the value\n");
      fprintf(stderr,"\tof MAX_STR_LEN (currently set to %d), ",MAX_STR_LEN);
      fprintf(stderr,"and recompile.\n");
      exit(-1);
    }
  }
}



/* Read data file of gene orders in format: 
   >genome name 1 2 4 3 5 -6 -10 -9 -8 7 
   genes 1..n.  */ 

/* input             = FILE stream for input

on return:
   num_genes         = number of actual genes + 2*number of chromosomes
   num_genomes       = total number of genomes
   num_chromosomes   = # chromosomes per genome
                       (genomes initially having fewer chromosomes are
		       padded with null chromosomes to reach this #)
   *genome_list      = array containing all the genomes
*/

void mcread_data(FILE *input, struct mgr_genome_struct **genome_list,
		 int *num_genes, int *num_genomes,
		 int *num_chromosomes)
{
  int genome_num;
  int i;
  int c;
  char *buf;                /* for genome name */

  int gnum;                 /* # genes in this genome */
  int gnumc;                /* # genes in this chromosome */
  int cnum;                 /* # chromosomes listed in this genome */

  /* counts for all genomes */
  int NumGenomes, NumGenes;
  int NumChromosomes;       /* max # chromosomes per genome */

  int cap;                  /* #'s used to delimit chromosomes */

  NumGenomes = 0;
  NumGenes   = -1;
  NumChromosomes = 1;


  /* pass 1: read all the data, counting the # genomes, genes, chromosomes,
     so that memory can be allocated to store it */
  /* pass 2: read the data, store & validate it */

  /* pass 1 */
  gnum = gnumc = cnum = 0;

  while (1) {
    /* skip white space & get next character */
    while (((c=getc(input)) != EOF) && isspace(c)) ;

    if (c == '#') {                          /* comment */
      discard_line(input);
      continue;
    }

    if (c == '>' || c == EOF) {              /* new genome */
      if (NumGenomes == 1) {                 /* first genome complete */
	NumGenes = gnum;
	NumChromosomes = cnum;
      }

      if (NumGenomes > 0 &&
	  NumGenes != gnum) {   /* do not permit # real genes to vary */
	fflush(stdout);
	fclose(stdout);
	fprintf(stderr,"ERROR: mcread_data() inconsistent num of genes\n");
	fprintf(stderr,"ERROR: Expecting %d genes, but ",
		NumGenes);
	fprintf(stderr,"Genome %d has %d genes\n",
		NumGenomes, gnum);
	exit(-1);
      }

      /* find max # chromosomes */
      if (NumChromosomes < cnum)
	NumChromosomes = cnum;


      if (c == EOF) break;

      discard_line(input);                   /* ignore the genome name */
      NumGenomes++;
      gnum = gnumc = cnum = 0;
      continue;
    }

    if (NumGenomes>0) {
      if (c == ';' || c == '$') {            /* new chromosome */
	gnumc = 0;                           /* no genes in this chromosome */

	continue;
      }

      if (c == '-' || isdigit(c)) {            /* gene */
	gnum++; gnumc++;
	if (gnumc == 1) cnum++;          /* only count non-null chromosomes */

	while ((c = getc(input)) != EOF && isdigit(c))
	  ;
	if (c != EOF) ungetc(c,input);
      }

      continue;
    }

    /* character unknown */
    fflush(stdout); fclose(stdout);
    fprintf(stderr,"ERROR: Improper file format.  Bad character %c at position %d.\n",c,(int)ftell(input));
    exit(-1);
  }


  /* finished pass 1 */

  /* adjust and validate parameters and allocate memory
     before doing pass 2 */

  /*
  fprintf(stderr,"Pass 1: %d genes, %d chromosomes, %d genomes\n",
	  NumGenes, NumChromosomes, NumGenomes);
	  */




  /* Adjust number of genes to include caps */
  NumGenes += 2*NumChromosomes;


  if (NumGenomes > MAX_GENOMES) {
    fprintf(stderr,"ERROR: \tNumber of genomes in the input is %d, but\n",
	    NumGenomes);
    fprintf(stderr,"\tGRAPPA is compiled with a maximum of %d genomes.\n",
	    MAX_GENOMES);
    fprintf(stderr,"\tPlease edit the header file \"mcstructs.h\", ");
    fprintf(stderr,"increase the value\n");
    fprintf(stderr,"\tof MAX_GENOMES (currently set to %d), ",
	    MAX_GENOMES);
    fprintf(stderr,"and recompile.\n");
    exit(-1);
  }

  if (NumGenes > MAX_NUM_GENES) {
    fprintf(stderr,"ERROR: \tNumber of genes in the input is %d, and\n",
	    NumGenes-2*NumChromosomes);
    fprintf(stderr,"ERROR: \tnumber of chromosomes is %d, but\n",
	    NumChromosomes);
    fprintf(stderr,"\tMCDIST is compiled with a maximum of %d genes (=genes+2*chromosomes).\n",
	    MAX_NUM_GENES);
    fprintf(stderr,"\tPlease edit the header file \"mcstructs.h\", ");
    fprintf(stderr,"increase the value\n");
    fprintf(stderr,"\tof MAX_NUM_GENES (currently set to %d), ",
	    MAX_NUM_GENES);
    fprintf(stderr,"and recompile.\n");
    exit(-1);
  }

  rewind(input);

  *genome_list =
    (struct mgr_genome_struct *) e_malloc(NumGenomes*sizeof(struct mgr_genome_struct),
				      "genome_list");

  for (i=0 ; i<NumGenomes ; i++) {
    (*genome_list)[i].gnamePtr = (char *) e_malloc(MAX_STR_LEN*sizeof(char),
						   "gnamePtr");
  }


  /* pass 2: read the data, store & validate it */
  genome_num = -1;
  gnum = gnumc = cnum = 0;
  cap = NumGenes-2*NumChromosomes+1;

  while (1) {
    /* skip white space & get next character */
    while (((c=getc(input)) != EOF) && isspace(c)) ;

    if (c == '#') {                          /* comment */
      discard_line(input);
      continue;
    }

    if (c == '>' || c == EOF) {             /* new genome */
      /* complete & validate genome just read in */
      if (genome_num >= 0) {
	/* Terminate the current chromosome, and pad with
	   null chromosomes if necessary to reach the
	   total number of chromosomes there should be */

	while (cap <= NumGenes) {
	  (*genome_list)[genome_num].genes[gnum++] = cap++;
	}

	/* validate it */
	check_genes((*genome_list)[genome_num].genes, NumGenes,
		    genome_num+1);
      }


      /* finish if EOF, continue if new genome */
      if (c == EOF) break;

      /* advance genome counter */
      genome_num++;

      /* read in the genome name */
      buf = (*genome_list)[genome_num].gnamePtr;
      fgets(buf,
	    MAX_STR_LEN,
	    input);
      /* delete terminal \n */
      while (*buf != '\0') {
	if (*buf == '\n') {
	  *buf = '\0'; break;
	}
	buf++;
      }

      /* index # of genome */
      (*genome_list)[genome_num].genome_num = genome_num+1;

      /* allocate space to store the genes */
      (*genome_list)[genome_num].genes =
	(int *) e_malloc(NumGenes*sizeof(int), "genes");

      gnum = gnumc = cnum = 0;
      cap = NumGenes-2*NumChromosomes+1;
      continue;
    }

    if (genome_num>=0) {
      if (c == ';' || c == '$') {            /* new chromosome */
	gnumc = 0;                           /* no genes in this chromosome */

	continue;
      }

      if (c == '-' || isdigit(c)) {            /* gene */
	gnumc++;
	if (gnumc == 1) {               /* only count non-null chromosomes */
	  /* close previous chromosome, if there was one */
	  if (cnum > 0)
	    (*genome_list)[genome_num].genes[gnum++] = cap++;

	  /* start new chromosome */
	  (*genome_list)[genome_num].genes[gnum++] = cap++;
	  
	  /* count # chromosomes */
	  cnum++;
	}

	/* read in the gene */
	ungetc(c,input);             /* need to put the lead digit back in */
	fscanf(input,"%d",
	       &(*genome_list)[genome_num].genes[gnum++]);
      }

      continue;
    }

    /* character unknown */
    fflush(stdout); fclose(stdout);
    fprintf(stderr,"ERROR: Pass 2, Improper file format.  Bad character %c at position %d.\n",c,(int)ftell(input));
    exit(-1);
  }

  /**********************************************************************/
  /*                     finished reading input                         */
  /**********************************************************************/

  fclose(input);

  *num_genes   = NumGenes;
  *num_genomes = NumGenomes;
  *num_chromosomes = NumChromosomes;

  return;
}


/* delete the caps that mcread_data put in */
void delete_caps(struct mgr_genome_struct *genome,
		 int num_genes, int num_chromosomes)
{
  int source, dest;
  int *genes = genome -> genes;
  int lowcap = num_genes - 2*num_chromosomes + 1;

  dest=0;
  for (source=0; source < num_genes; source++) {
    if (genes[source]<lowcap && genes[source] > -lowcap) {
      genes[dest++] = genes[source];
    }
  }
}


/* verify that user entered caps are valid
   Number duplication/omission was already checked by mcread_data
   This just does parity/sign checks.  It doesn't check that
   multiple genomes are properly co-capped.

   return: 0=invalid, 1=valid
*/
int check_caps(struct mgr_genome_struct *genome,
	       int num_genes, int num_chromosomes)
{
  int source;
  int *genes = genome -> genes;
  int lowcap = num_genes - 2*num_chromosomes + 1;
  int inchrom = 0;

  int s,g;

  for (source=0; source < num_genes; source++) {
    /* pick apart gene # and sign */
    g = genes[source];
    if (g>0)
      s=0;
    else {
      g = -g; s = 1;
    }

    /* is it a cap? */
    if (g < lowcap) {                /* not a cap */
      if (!inchrom) return FALSE;    /* if not in chromosome,
					there MUST be a cap here */

    } else {                         /* is a cap */
      if ((s + g + inchrom + num_genes) % 2 == 0)
	return FALSE;                /* doesn't have correct sign */

      inchrom = !inchrom;
    }
  }

  return TRUE;
}


// hit
/* reformat a genome so that its null chromosomes are at the
   end, and the cap #'s are consecutive */
void reformat_caps(struct mgr_genome_struct *genome,
		   int num_genes, int num_chromosomes)
{
  int source, dest;
  int *genes = genome -> genes;
  int lowcap = num_genes - 2*num_chromosomes + 1;
  int cap = lowcap;

  int inchrom = 0;
  /* 0: awaiting opening cap
     1: read opening cap, but haven't read a gene in the chromosome yet
     2: have read a gene
     */
  

  int g;

  /* doesn't apply in unichromosomal case */
  if (num_chromosomes == 0)
    return;

  for (source=0, dest=0; source < num_genes; source++) {
    g = genes[source];

    /* is it a cap? */
    if (g >= lowcap || g <= -lowcap) {
      switch(inchrom) {
      case 0: /* awaiting opening cap */
	inchrom = 1; break;

      case 1: /* read opening cap, but haven't read a gene yet */
	/* it was a null chromosome */
	break;

      case 2: /* have read a gene */
	/* close chromosome */
	genes[dest++] = cap++;
	inchrom = 0; /* get ready to start next chromosome */
	break;
      }

    } else {       /* not a cap */
      if (inchrom==1) { /* first gene in chromosome */
	inchrom++;                   /* chromosome has a non-cap entry */
	genes[dest++] = cap++;       /* add opening cap */
      }
      genes[dest++] = genes[source]; /* copy the current gene */
    }
  }

  /* put null chromosomes at the end of the genome */
  while (dest < num_genes)
    genes[dest++] = cap++;


#if 0
  fprintf(stdout,"reformat_genes: ");
  for (source = 0; source < num_genes ; source++)
    fprintf(stdout,"%d ",genes[source]);
  fprintf(stdout,"\n");
  fflush(stdout);
#endif
}


