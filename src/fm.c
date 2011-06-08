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
#include <stdlib.h>         //malloc, calloc, srand, RAND_MAX, posix_memalign
#include <string.h>
#include <assert.h>         
#include <math.h>           //log10,exp

#include "config.h"         //defines, if likelihood, USE_LIKELIHOOD
#ifdef USE_LIKELIHOOD

#include <caml/alloc.h>     //copy_double, et cetera
#include <caml/mlvalues.h>
#include <caml/memory.h>    //caml_param, et cetera
#include <caml/bigarray.h>
#include <caml/custom.h>    //abstract types
#include <caml/intext.h>    //serialization
#include <caml/fail.h>      //failwith('')

#include "fm.h"

/** CONSTANTS **/
#define EPSILON         1e-10
#ifndef INFINITY
    #define INFINITY    1/0
#endif

/** defines for ease of use **/
#define MAX(a,b) (a >= b)?a:b
#define MIN(a,b) (a <= b)?a:b
#define EQUALS(a,b) abs(a-b) < EPSILON
#define ISSET(x,i) 0<((1 << i)&x)

/** For unpacking ML abstract/option types **/
#define FM_val(v) (*((struct matrix**)Data_custom_val(v)))
#define FM_ptr(v)   ((struct matrix**)Data_custom_val(v))
#define Val_none Val_int(0)
#define Some_val(v) Field(v,0)

/* This function is in likelihood.c
#ifdef _WIN32
__inline
#else
inline
#endif
value Val_some( value v )
{   
    CAMLparam1( v );
    CAMLlocal1( some );
    some = caml_alloc(1, 0);
    Store_field( some, 0, v );
    CAMLreturn( some );
}
*/

/** A file to compose and create matrices for the costs of alignment of
 * likelihood characters. The file standards are similar to those in the
 * static likelihood file. (that is variables that start with capital letters
 * are matrices, and lower case are vectors, t is a time/distance, and a/b/c
 * define sequence data for the alignments. **/

/** [neg_log_comp mat n m] Take the negative log of all the components of a
 * matrix [mat] of size [n] and [m] **/
void neg_log_comp( double * mat, const int n, const int m )
{
    int i;
    for( i =0; i < n*m; ++i )
        mat[i] = - log( mat[i] );
}

/** A preliminary function to calculate the cost of an alignment of two
 * characters [x] and [y] with composed P matrices [X] and [Y], respectively,
 * over alphabet of size [n]. The optimal assignment in [a] as a bitset. **/
double calculate_cost( CDIR* a_ptr, const double* X, const double* Y,
                            const int x, const int y, const int n )
{
    double best_cost, temp_cost;
    int i, j, k;
    CDIR a;

    a = 0;
    best_cost = INFINITY;
    for(i = 0; i < n; ++i){
        if( ISSET( x, i ) ){
            for( j = 0; j < n; ++j ){
                if( ISSET( y, j ) ){
                    for( k = 0; k < n; ++k ){
                        temp_cost = X[i*n + k] + Y[j*n + k];
                        if( EQUALS( temp_cost, best_cost ) ){
                            best_cost = MIN(temp_cost, best_cost);
                            a |= (1 << k);
                        } else if( temp_cost < best_cost ){
                            a = (1 << k);
                            best_cost = temp_cost;
                        }
                    }
                }
            }
        }
    }
    *a_ptr = a;
    return best_cost;
}

/** Fill in an expanded cost matrix of all combination of states; requires that
 * fm be allocated to the appropriate dimensions --reflect comb and size parts.
 * The assignment is discovered through an MPL cost formulation (the assignment
 * is not necessarily a decendents state. **/
void precalc( fm *FM, const double *A, const double *B )
{
    int i, j, k, max;
    CDIR temp_best_asgn;
    double temp_best_cost, temp_cost;

    /* Build an NxN matrix; Combine the two matrices to create a cost and
     * assignment between Row and Col from matrix A and B respectively. */
    if( 0 == FM->comb ){
        max = FM->size;
        for(i = 0; i < max; ++i ){
            for(j = 0; j < max; ++j ){
                temp_best_asgn = 0;
                temp_best_cost = INFINITY;
                for(k = 0; k < max; ++k){
                    temp_cost = A[i*max+k] + B[j*max+k];
                    if( EQUALS( temp_cost, temp_best_cost )){
                        temp_best_cost = MIN(temp_cost,temp_best_cost);
                        temp_best_asgn |= (1 << k);
                    } else if( temp_cost < temp_best_cost ){
                        temp_best_cost = temp_cost;
                        temp_best_asgn = (1 << k);
                    }
                }
                FM->cost[i*max+j]   = temp_best_cost;
                FM->assign[i*max+j] = temp_best_asgn;
            }
        }
    /* Build the matrix 2^N * 2^N, to create matrix of all combinations of
     * assignment of Row and Col (again, combining both matrices). **/
    } else {
        max = 1 << (FM->size);
        for( i = 0; i < max; ++i ){
            for( j = 0; j < max; ++j ){
                temp_best_cost = calculate_cost( &temp_best_asgn, A, B, (i+1), (j+1), FM->size );
                FM->cost[i*max+j]   = temp_best_cost;
                FM->assign[i*max+j] = temp_best_asgn;
            }
        }
    }
}

/** Wrapper for creating a pre-calculated cost matrix from instances; return
 * costs and assignments into individual bigarrays, and return the pair. **/
value
fm_CAML_compose( value FM, value U, value D, value Ui, value cta, value ctb, value full)
{
    CAMLparam5( FM, U, D, Ui, cta );
    CAMLxparam2( ctb, full );
    CAMLlocal3( assign, costs, pair );

    double ta, tb, *PA,*PB, *TMP;
    int alph;
    fm *results;
    mat *space;
    long dims[2]; /** both matrices are the same size **/

    ta = Double_val( cta );
    tb = Double_val( ctb );

    results = (fm*) malloc( sizeof(fm) );
    space = FM_val( FM );
    alph = Bigarray_val(U)->dim[0];
    expand_matrix( space, (3*alph*alph) );
    PA = register_section( space, alph*alph, 0 );
    PB = register_section( space, alph*alph, 0 );
    TMP = register_section( space, alph*alph, 0 );

    /** fill in a new fm to hold data for results **/
    results->size = alph;
    results->comb = Int_val(full);
    if( 0 == results->comb ){
        dims[0] = alph; dims[1] = alph;
        results->cost = (double*) malloc( alph * alph * sizeof(double));
        results->assign = (CDIR*) malloc( alph * alph * sizeof(CDIR));
    } else if (1 == results->comb){
        dims[0] = dims[1] = (1 << alph);
        results->cost = (double*) malloc( (1 << alph) * (1 << alph) * sizeof(double) );
        results->assign = (CDIR*) malloc( (1 << alph) * (1 << alph) * sizeof(CDIR) );
    } else {
        failwith( "fm_CAML_compose :: value for determining combinations" );
    }
    
    /** pre-calculate the matrix **/
    if( Val_none == Ui ){
        compose_sym( PA, Data_bigarray_val(U), Data_bigarray_val(D), ta, alph, TMP );
        compose_sym( PB, Data_bigarray_val(U), Data_bigarray_val(D), tb, alph, TMP );
    } else {
        compose_gtr( PA, Data_bigarray_val(U), Data_bigarray_val(D),
                     Data_bigarray_val(Some_val(Ui)), ta, alph, TMP );
        compose_gtr( PB, Data_bigarray_val(U), Data_bigarray_val(D),
                     Data_bigarray_val(Some_val(Ui)), tb, alph, TMP );
    }
    neg_log_comp( PA, alph, alph );
    neg_log_comp( PB, alph, alph );
    precalc( results, PA, PB );

    /** compose into ocaml values and return **/
    assign = alloc_bigarray( ODIR | BIGARRAY_C_LAYOUT, 2, results->assign, dims );
    costs  = alloc_bigarray( BIGARRAY_FLOAT64  | BIGARRAY_C_LAYOUT, 2, results->cost, dims );
    pair = caml_alloc_tuple( 2 );
    Store_field( pair, 0, assign );
    Store_field( pair, 1, costs );
    CAMLreturn( pair );
}

value
fm_CAML_compose_wrapper(value * argv, int argn)
{
    return fm_CAML_compose
            (argv[0], argv[1], argv[2], argv[3], argv[4], argv[5], argv[6]);
}

#endif /* USE_LIKELIHOOD */
