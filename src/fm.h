/* POY 5.0 Alpha. A phylogenetic analysis program using Dynamic Homologies.   */
/* Copyright (C) 2011  Andrés Varón, Lin Hong, Nicholas Lucaroni, Ward Wheeler*/
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

#include "config.h"         //defines, if likelihood, USE_LIKELIHOOD
#ifdef USE_LIKELIHOOD

#include <caml/bigarray.h>

#include "seq.h"
#include "floatmatrix.h"
#include "likelihood.h"

/** How should the direction matrix be assigned? **/
#define CASN int
#define OASN BIGARRAY_INT32

/** tie up two matrices; one defines the cost, another the assignment **/
struct fmat {
    int size;       /** Size of the cost matrix; alph^N or alph **/
    int alph;       /** the alphabet size **/
    int comb;       /** Is the matrix that of combinations? 0=no, 1=yes **/
    int gap;        /** What is the gap state? **/
    double* cost;   /** this matrix holds either the single character
                        transformation, or a transformation of all combinations,
                        depending on the [comb] variable --if it equals 0 then
                        we are dealing with single character transformations **/
    int* cost_asgn; /** bitset matrix for the assignments; size as the cost matrix **/
};
typedef struct fmat fm;

void neg_log_comp( double * mat, const int n, const int m );

double fcost( fm *FCM, const CASN x, const CASN y );
CASN    fasn( fm *FCM, const CASN x, const CASN y );

double calculate_cost( CASN* a, const double* X, const double* Y,
                            const int x, const int y, const int n );

int get_closest( double* cost, const int p, const int m, const double* cm, const int alph );

void precalc( fm *FM, const double *A, const double *B );

#endif /* USE_LIKELIHOOD */
