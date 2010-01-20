/* genome_print.h
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
#ifndef GENOME_PRINT_H
#define GENOME_PRINT_H

#include "mgrstructs.h"

void special_print_genome(FILE *output, struct mgr_genome_struct *genome,
						  int actual_nb_chromo) ;



void special_print_genomes(FILE *output, G_struct *Genomes, int nb_spec) ;

void print_alphabet(a_strip *alphabet, int nb_strips);

void special_print_genome2(FILE *output, struct mgr_genome_struct *g,
						   int actual_nb_chromo) ;

void special_print_genomes2(FILE *output, G_struct *Genomes,
							int nb_spec) ;

void print_all_genomes(FILE *output, G_struct *Genomes, G_struct *Preancestors,
		       int num_genes, int num_chromosomes, int nb_spec) ;
		   
#endif /* DEFINE GENOME_PRINT_H*/
