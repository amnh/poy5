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
#include "config.h"         //defines, if likelihood, USE_LIKELIHOOD
#ifdef USE_LIKELIHOOD   
#include <stdio.h>
#include <stdlib.h> //malloc, srand, RAND_MAX
#include <string.h> //memcpy, memset
#include <assert.h>

#include <math.h>   //log,exp,sin
#ifdef __APPLE__
#include <vecLib/cblas.h>  //dgemm
#else
#include <cblas.h>  //dgemm
#endif

#include <caml/mlvalues.h>
#include <caml/memory.h>
#include <caml/bigarray.h>
#include <caml/custom.h>
#include <caml/intext.h>

#include "array_pool.h"
#include "gamma.h" //contains special functions

/**
 * C struct for a character set and likelihood vectors
 */ 
struct ml {
    int stride;     //stride of matrix (ie: size of alphabet)
    int c_len;      //length of character set
    double* lv_s;   //likelihood vectors x character set (length is n*k)
};
typedef struct ml mll;

struct pt {
    double ta;      //branch length of a
    double tb;      //branch length of b
    double ll;      //likelihood
    mll *vs;        //resultant vector
};
typedef struct pt ptr;

/**
 * testing functions...
 */ 
//print formated array
void printarray( const double* a, const int n );
//print formated matrix
void printmatrix( const double* Z, const int n, const int m );
//testing... create random substituation rate matrix
void rand_sub_mat_gtr( double* A, const int n);
//testing... create symmetric sub-rate matrix
void rand_sub_mat_sym( double* A, const int n);
//converts value -> struct ml and prints it's contents
void CAML_debug( value s );

/**
 *  custom functions for garbage collection
 */
void likelihood_CAML_free_likelihood( value v );
int likelihood_CAML_compare( value a,value b );
void likelihood_CAML_serialize(value v, unsigned long* wsize_32, unsigned long* wsize_64);
unsigned long likelihood_CAML_deserialize( void* dst );
value likelihood_CAML_register (value u0);
/**
 * help functions...
 */
// make diagonal matrix from vector
#ifdef _WIN32
__inline void
#else
inline void
#endif
mk_diag(const double* diag, double* M, const int n, const int m);
// applies exp() to a diagonal [n]x[m] matrix, [D], with constant scalor [t]
#ifdef _WIN32
__inline void
#else
inline void
#endif
apply_exp( double* D, const int n, const int m, const double t);

/**
 * diagonalization functions
 */
//reconstructs a diagonalized matrix after applying exp(D*t)
//  P is *output*.. thus: U*exp(D*t)*Ut
void compose_sym( double* P, const double* U, const double* D,
                        const float t, const int n,double* tmp);
//  P is *output*.. thus: U*exp(D*t)*Ui
void compose_gtr( double* P, const double* U, const double* D,
                        const double* Ui, const double t, const int n, double* tmp);
//diagonalize a symmetric [n]x[n] matrix [A] into A*D*At 
//  A is overwritten as eigenvectors
//  D is a diagonal matrix of eigenvalues --assumed clean
int diagonalize_sym( double* A, double* D, const int n );
int diagonalize_gtr( double* A, double* D, double* Ui, const int n );

/**
 * median functions
 */
//calculates the half median 
void median_h( const double* P, const double* l,const int c, double* nl, const int a);
//calculates median of two character sets, [a] and [b] into [c]
void median_charset( const double* Pa, const double* Pb, const struct ml* a,
                    const struct ml* b, struct ml* c, const double p );
//loglikelihood of a median with priors [p] and character set [l]
double loglikelihood( const struct ml *l, const double* p );

/**
 * functions on character sets
 */ 
void filter_s( const struct ml* c1, struct ml* c2, const int* idxs, const int n );
int compare_chars( const struct ml* c1, const struct ml* c2 );

/**
 * readjust functions
 */
void //0 for no change, 1 for change
readjust_walk_gtr(const double* U,const double* D,const double* Ui,
                   const struct ml* a,const struct ml* b,struct ml* b_c,
                    double* ta,double* tb,double* b_mle,const double* rates,
                     const double *p, const int g_n, const double* pi);
void //0 for no change, 1 for change
readjust_walk_sym(const double* U,const double* D,
                   const struct ml* a,const struct ml* b,struct ml* b_c,
                    double* ta,double* tb,double* b_mle,const double* rates,
                     const double *p, const int g_n,const double pi);
#endif /* USE_LIKELIHOOD */
