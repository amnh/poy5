/* graph_edit.h
*    Build the breakpoint graph, and routines to reverse portions of
*    it and add edges.
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

void double_perms(struct mgr_genome_struct *g1, struct mgr_genome_struct *g2, int offset, int num_genes, mgr_distmem_t *distmem);

void classify_cycle_path(graph_t *G,
			 pathcounts_t *pcounts);

void flip_chromosomes(graph_t *G, int cnum_start, int cnum_end,
		      int updatecycleinfo);

void flip_chromosomes_nocycles(graph_t *G, int cnum_start, int cnum_end);

void flip_chromosome_wvertex(graph_t *G, int i);

void test_flipping(graph_t *G, int num_chromosomes);

void closepath(graph_t *G, int index, pathcounts_t *pcounts);

void addedge(graph_t *G, int i, int j, pathcounts_t *pcounts);

void join_paths_2_cycle(graph_t *G,
			int v1, int v2, int v3, int v4,
			pathcounts_t *pcounts);


void join_pgpaths_2_cycle(graph_t *G,
			  int i, int j,
			  pathcounts_t *pcounts);

void print_arrays(graph_t *G);

int graphdist(graph_t *G, int i, int j);









