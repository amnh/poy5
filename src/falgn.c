/* POY 4.0 Beta. A phylogenetic analysis program using Dynamic Homologies.    */
/* Copyright (C) 2007  Andrés Varón, Le Sy Vinh, Illya Bomash, Ward Wheeler,  */
/* and the American Museum of Natural History.                                */
/*                                                                            */
/* This program is free software; you can redistribute it and/or modify       */
/* it under the terms of the GNU General Public License as published by       */
/* the Free Software Foundation; either version 2 of the License, or          */
/* (at your option) any later version.                                        */
/*                                                                            */
/* This program is distributed in the hope that it will be useful,            */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of             */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the              */
/* GNU General Public License for more details.                               */
/*                                                                            */
/* You should have received a copy of the GNU General Public License          */
/* along with this program; if not, write to the Free Software                */
/* Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301   */
/* USA                                                                        */

#include "floatmatrix.h"
#include "likelihood.h"
#include "seq.h"
#include "falgn.h"

#define dmat DIRECTION_MATRIX

void full_align( seq x, seq y, mat *f_space, dmat *d_space )
{






}


void init_align( seq x, seq y, mat f_space, dmat d_space, mat U, mat D, mat UI, int tx, int ty)
{
    int work;

    /* calculate the necessary space and ensure we have it in matrices */
    work = x.len * y.len + 2 * 


}
