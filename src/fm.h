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

#include <stdio.h>
#include <string.h>
#include <assert.h>         

#include <caml/alloc.h>     //copy_double, et cetera
#include <caml/mlvalues.h>
#include <caml/memory.h>    //caml_param, et cetera
#include <caml/bigarray.h>
#include <caml/custom.h>    //abstract types
#include <caml/intext.h>    //serialization
#include <caml/fail.h>      //failwith('')

#include "floatmatrix.h"
#include "likelihood.h"

/** tie up two matrices; one defines the cost, another the assignment **/
struct fmat {
    int size;       /** size of the alphabet  < {32,64} **/
    int comb;       /** is the matrix that of combinations? 0=no, 1=yes **/
    int metric;     /** does the matrix follow metricity? 0=no, 1=yes **/
    
    double* cost;   /** this matrix holds either the single character
                        transformation, or a transformation of all combinations,
                        depending on the [comb] variable --if it equals 0 then
                        we are dealing with single character transformations **/
    int* assign;    /** bitset matrix for the assignments; size as the cost matrix **/
};
typedef struct fmat fm;

double calculate_cost( int* a, const double* X, const double* Y,
                            const int x, const int y, const int n );

void precalc( fm *FM, const double *A, const double *B, const int n );

