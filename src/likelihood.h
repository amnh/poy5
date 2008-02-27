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

/**
 * C struct for a character set and likelihood vectors
 */ 
struct ml {
    int stride;       //stride of matrix (ie: size of alphabet)
    int c_len;        //length of character set
    double* lv_s;     //likelihood vectors x character set (length is n*k)
};

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
 * functions for gamma distribution
 */
//gamma function by Lanczos Approx (~10e-15 relative error)
double gamma( const double z );
//gamma density function for: int( pi*f(x|t*r)*gd(r), 0 , inf, r )
double gamma_pdf( const double r, double alpha, double beta );
//calculates the percetage points for k categories
void gamma_pp( double* out, const int k, const double alpha, const double beta);
//confluent hypergeometric (for igamma)
double gamma_M( const double z, const double a);
//generates rates
void gamma_rates(double* rates, const double alpha, const double beta, 
                    const double* cuts, const int k);
//incomplete gamma ratio; based on AS32 by Bhattacharjee (1970)
double gamma_i( const double x, const double a);

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
void median_h( const double* P, const double* l, double* nl, const int a);
//calculates median of two character sets, [a] and [b] into [c]
void median_charset( const double* Pa, const double* Pb, const struct ml* a,
                    const struct ml* b, struct ml* c, const double r );
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
//double readjust(const struct ml* ca, const struct ml* cb, struct ml*cc,
//                     const double ta, const double tb, const double o_mle)
