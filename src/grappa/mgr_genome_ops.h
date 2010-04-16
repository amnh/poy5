/* genome_ops.h
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

#ifndef GENOME_OPS_H
#define GENOME_OPS_H

#include "mgrstructs.h"

void initialize_alphabet(a_strip **alphabet, int num_g, int identity);

void copy_alphabet(a_strip *alphabet1, a_strip **alphabet2, int num_g1, int num_g2) ;

void free_alphabet(a_strip *alphabet, int num_g);

void check_alpha(a_strip *alphabet, int num_g, int ori_num_g);

void copy_genomes(G_struct *Genomes1, int gindex1,
				   G_struct *Genomes2, int gindex2,
				   int need_copy_alphabet) ;

void copy_first_n(G_struct *Genomes1, G_struct *Genomes2,
				  int nb_spec) ;
					


void update_for_merge(G_struct *Genomes, int *nbreag, mgr_distmem_t *distmem,
					  treemem_t *treemem, int gindex1, int gindex2, int nb_spec, 
					  int spec_left, int verbose) ;

void check_for_merge(G_struct *Genomes, int *nbreag, 
						 int nb_spec, int *spec_left, 
						 mgr_distmem_t *distmem, treemem_t *treemem, int verbose) ;

void initialize_genome_list(struct mgr_genome_struct **genome_list, int nb_spec,
							int num_genes, int num_chromosomes);


void free_genome_list(struct mgr_genome_struct *genome_list, int nb_spec) ;


void print_G_struct(G_struct *Genomes, int nb_spec);

void init_G_struct(G_struct *Genomes, struct mgr_genome_struct *genome_list,
						 int nb_spec, int num_genes,
						 int num_chromosomes) ;

void free_G_struct(G_struct *Genomes, int nb_spec);



void print_dist_mat(G_struct *Genomes, int nb_spec) ;

void compute_dist_mat(G_struct *Genomes, int nb_spec,
						 mgr_distmem_t *distmem);

void local_distmat_update(int gindex, G_struct *Genomes, int nb_spec, 
						  mgr_distmem_t *distmem);

void print_dist_mat(G_struct *Genomes, int nb_spec) ;

int find_total_dist(G_struct *Genomes, int nb_spec,  
					mgr_distmem_t *distmem) ;
					
void add_null_chromos(struct mgr_genome_struct *genome_list, int nb_spec,
		      int *num_genes, int *num_chromosomes) ;

void remove_null_chromos (struct mgr_genome_struct *genome_list, int nb_spec,
					  int *num_genes, int *num_chromosomes, int alpha_size, 
                      int max_num_deli) ;

void find_max_chromo_size(G_struct *Genomes, int nb_spec);


void carry_on_reag(list_reag *the_list, G_struct *Genomes, int nb_spec,
				  int *spec_left, 
				   int *did_merge, mgr_distmem_t *distmem);


void undo_reag(list_reag *the_list, G_struct *Genomes, int nb_spec,
			  int *spec_left, 
			   int did_merge, mgr_distmem_t *distmem); 

void print_one_reag(list_reag *the_list, G_struct *Genomes,
					int nb_spec, struct mgr_genome_struct *old_genome);

void condense_genomes(G_struct *Genomes, int nb_spec, int verbose);

void uncondense_genomes(G_struct *Genomes, int nb_spec, int verbose);

#endif /* IFNDEF GENOME_OPS_H*/
