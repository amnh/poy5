/* list_ops.h
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


cbounds_t cb_compute_Rscore; //mgr_list_ops, "compute_Rscore"

list_edge *RemoveItem_edge(list_edge *the_list);

void erase_list_edge(list_edge *the_list);

list_reag *RemoveItem_reag(list_reag *the_list);

void erase_list_reag(list_reag *the_list);

void print_list_edge(list_edge *the_list);

void merge_two_list(treemem_t *treemem);

list_edge *list_edge_insert(list_edge *the_list, int node1, int label1,
			    int node2, int label2, int score);

list_reag *list_reag_insert(list_reag *the_list, int spec, int optype, int start, int end,
			    int sc1, int sc2, int merge_with, int reduction, int *changes, int nb_spec);

void print_list_reag(list_reag *the_list, int nb_spec) ;

int rear_length(list_reag *the_list);

int total_reduction(list_reag *the_list) ;

int total_edge_score(list_edge *the_list);

void add_nbrev_edges(list_edge *the_edge_list, int *nbrev, int nb_spec);


void check_list_edge(G_struct *Genomes, list_edge *the_edge_list, mgr_distmem_t *distmem, int verbose) ;


void compute_Rscore(list_reag *the_list, G_struct *Genomes);
