/* mcread_input.h
*    Read input genomes
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

/* Modifed from GRAPPA 1.02: read_input.h */

#ifndef READ_INPUT_H
#define READ_INPUT_H


void check_genes(int *genes, int num_genes,
		 int at_genome);

void mcread_data(FILE *input, struct mgr_genome_struct **genome_list,
		 int *num_genes, int *num_genomes,
		 int *num_chromosomes);

void delete_caps(struct mgr_genome_struct *genome,
		 int num_genes, int num_chromosomes);

int check_caps(struct mgr_genome_struct *genome,
	       int num_genes, int num_chromosomes);

void reformat_caps(struct mgr_genome_struct *genome,
		   int num_genes, int num_chromosomes);


#endif
