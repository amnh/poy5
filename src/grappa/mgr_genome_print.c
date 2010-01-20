/* genome_print.c
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
//#include "mgrstructs.h" //include this in mgr_genome_print.h
#include "mgr_genome_print.h"

void special_print_genome(FILE *output, struct mgr_genome_struct *g,
						  int actual_nb_chromo) 
{	
    int j, max_pos;	
    int lowcap = g->num_g - 2*g->num_chr + 1;	
    int inchromo = FALSE;	
    int gene;
    
    if (g->gnamePtr != (char *)NULL)
	fprintf(output,">%s\\n", g->gnamePtr);
    
    if (g->num_chr > 0)
	max_pos = g->num_g/3 + 2*actual_nb_chromo;	// to avoid all the nulls at the end		
    else max_pos = g->num_g;
    
    for (j=0; j<max_pos; j++) {
	gene = g->genes[j];
	if (abs(gene) >= lowcap) {
	    inchromo=!inchromo;
	    if (!inchromo && j<max_pos-1)
		fprintf(output, "\\$\\n");
	}
	else
	    fprintf(output, "%d ", gene);
    }
    if (g->num_chr > 0)
		fprintf(output, "\\$ ");
    fprintf(output, "\\n");
    fflush(output);
}

void special_print_genomes(FILE *output, G_struct *Genomes, 
						   int nb_spec) 
{
	int i;
	for (i=0; i<nb_spec; i++) {
		fprintf(output, "\t\t\"");
		special_print_genome(output, &(Genomes->genome_list[i]), 
							 Genomes->nb_chromo[i]);
		fprintf(output, "\",\n");
	}
}

void print_alphabet(a_strip *alphabet, int nb_strips) {
	int i, j;
	
	fprintf(stdout, "alphabet:\n");
	if (alphabet != NULL) {
		for (i=0; i<nb_strips; i++) {
			fprintf(stdout, "%d: ", i+1);
			if (alphabet[i].strip != NULL) {
				for (j=0;j<=alphabet[i].last;j++)
					fprintf(stdout, "%d ", alphabet[i].strip[j]);
				fprintf(stdout, "\n");
			}
		}
	}
}

void special_print_genome2(FILE *output, struct mgr_genome_struct *g,
						   int actual_nb_chromo) 

{	
    int j, max_pos;	
    //int lowcap = num_genes - 2*num_chromosomes + 1;	
	int lowcap = g->num_g - 2*g->num_chr + 1;
    int inchromo = FALSE;	
    int gene;
	int alphabet_size;
	
    if (g->gnamePtr != (char *)NULL)
		fprintf(output,">%s\n", g->gnamePtr);
    
/*    if (num_chromosomes > 0)
		max_pos = num_genes/3 + 2*actual_nb_chromo;	// to avoid all the nulls at the end		
    else max_pos = num_genes;
*/
	if (g->num_chr > 0)
		max_pos = g->num_g/3 + 2*actual_nb_chromo;	// to avoid all the nulls at the end		
    else max_pos = g->num_g;
	
	//printf("num_chromo, num_g, max_pos, actual_nb_chromo: %d %d %d %d\n", num_chromosomes, g->num_g, max_pos, actual_nb_chromo);
	
    
    for (j=0; j<max_pos; j++) {
		gene = g->genes[j];
		if (abs(gene) >= lowcap) {
			inchromo=!inchromo;
			if (!inchromo && j<max_pos-1)
				fprintf(output, "$\n");
		}
		else {
			fprintf(output, "%d ", gene);
		}
	}
    if (g->num_chr > 0)
		fprintf(output, "$ ");
    fprintf(output, "\n");
    fflush(output);
	
	if (g->num_chr > 0)
		alphabet_size = g->num_g/3;
	else
		alphabet_size = g->num_g;
	
	//print_alphabet(g->alphabet, alphabet_size);
}

void special_print_genomes2(FILE *output, G_struct *Genomes, 
							int nb_spec) 
{
    int i;
	for (i=0; i<nb_spec; i++) {
		//fprintf(output, "%d: num_g=%d ori_num_g=%d\n", i, Genomes->genome_list[i].num_g, Genomes->genome_list[i].ori_num_g);
		//fprintf(output, "nb_chromo=%d, same_as=%d\n", Genomes->nb_chromo[i], Genomes->same_as[i]);
		special_print_genome2(output, &(Genomes->genome_list[i]),
							  Genomes->nb_chromo[i]);
		//fprintf(output, "\n");
	}
}

void print_all_genomes(FILE *output, G_struct *Genomes, G_struct *Preancestors,
		       int num_genes, int num_chromosomes, int nb_spec) 
{
    int i;

    // the following assumes that nb_spec == 3
    if (nb_spec != 3) {
        fprintf(stderr, "print_all_genomes in genome_print.c assumes nb_spec=3\n");
        exit(-1);
    }
    
    for (i=0; i<2; i++) {
	special_print_genome2(output, &(Genomes->genome_list[i]),
						  Genomes->nb_chromo[i]);
	fprintf(output, "\n");
	special_print_genome2(output, &(Preancestors->genome_list[i]),
						  Preancestors->nb_chromo[i]);
	fprintf(output, "\n");
    }
    
    // this is where we print the ancestor
    special_print_genome2(output, &(Genomes->genome_list[3]),
						  Genomes->nb_chromo[3]);

    fprintf(output, "\n");

    	special_print_genome2(output, &(Preancestors->genome_list[2]),
							  Preancestors->nb_chromo[2]);
		fprintf(output, "\n");
	
		special_print_genome2(output, &(Genomes->genome_list[2]),
							  Genomes->nb_chromo[2]);
		fprintf(output, "\n");


}



