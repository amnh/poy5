/* graph_components.h
*    Form components of overlap graph.
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


void form_connected_components(graph_t *G);

void classify_connected_components(graph_t *G);



/* component types for computing hurdles */
#define C_U      0
#define C_IU     1
#define C_RU     2
void calc_num_hurdles_and_fortress(graph_t *G,
				   int comp_type, /* U, IU, RU */
				   hurdlecounts_t *counts);

int calc_num_semirealknots(graph_t *G, int *has_sgrk);

int calc_num_semiknots(graph_t *G);

int calc_num_simple(graph_t *G);



