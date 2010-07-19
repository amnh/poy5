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

//#ifdef __APPLE__
//#include <vecLib/cblas.h>  //dgemm
//#include <vecLib/clapack.h>
//#else
//#include <cblas.h>  //dgemm
//#include <clapack.h>
//#endif

#include <caml/mlvalues.h>
#include <caml/memory.h>
#include <caml/bigarray.h>
#include <caml/custom.h>
#include <caml/intext.h>

/**
 * C struct for a character set and likelihood vectors
 */ 
struct ml {
    int stride;     //stride of matrix (ie: size of alphabet)
    int c_len;      //length of character set
    int rates;      //number of rates being processed
    double* lv_s;   //likelihood vector for characters/rates, length = stride*c_len*rates
    int invar;      //1 if invar, else <= 0, needed for serlization
    int* lv_invar;  //likelihood vector for invariants, length = c_len
};
typedef struct ml mll;

struct pt {
    double time;    //branch length of b
    double ll;      //likelihood
    mll *vs;        //resultant vector
};
typedef struct pt ptr;

#define DECLARE_LAPACK_ROUTINES 1

#ifdef DECLARE_LAPACK_ROUTINES
/** double-precision-real lapack/blas routines **/
// double symmetric eigenvalue/vector routine
int dsyev_( char *jobz, char *uplo, int *n, double *a, int *lda, double *w,
            double *work, int *lwork, int *info);
// double general eigenvalue/vector routine
int dgeev_( char *jobvl, char *jobvr, int *n, double * a, int *lda, double *wr,
            double *wi, double *vl, int *ldvl, double *vr, int *ldvr,
            double *work, int *lwork, int *info);
//double general real LU-factorization 
int dgetrf_(int *m, int *n, double *a, int * lda, int *ipiv, int *info);
//double general real matrix inverse using LU
int dgetri_(int *n, double *a, int *lda, int *ipiv, double *work,
            int *lwork, int *info);
int dgemm_( char *transa, char *transb, const int *m, const int *n, const int *k,
            double *alpha, const double *a, const int *lda, const double *b,
            const int *ldb, double *beta, double *c, const int *ldc);
#endif

/**
 * testing functions...
 */ 
//print formated array
//void printarray( const double* a, const int n );
//print formated matrix
//void printmatrix( const double* Z, const int n, const int m );
//testing... create random substituation rate matrix
//void rand_sub_mat_gtr( double* A, const int n);
//testing... create symmetric sub-rate matrix
//void rand_sub_mat_sym( double* A, const int n);
//converts value -> struct ml and prints it's contents
//void CAML_debug( value s );

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
mk_diag(double* M, const int n);

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
int diagonalize_sym( mat *space, double* A, double* D, int n );
int diagonalize_gtr( mat *space, double* A, double* D, double* Ui, int n );

/**
 * median functions
 */
//calculates the half median 
//void median_h( const double* P, const double* l,const int c, double* nl, const int a);
//calculates median of two character sets, [a] and [b] into [c]
void median_MAL( const double* Pa, const double* Pb, const struct ml* a,
                        const struct ml* b, struct ml* c, const int rate_idx );
void median_MPL( const double* Pa, const double* Pb, const struct ml* a,
                        const struct ml* b, struct ml* c, const int rate_idx );
//loglikelihood of a median with priors [p], probabilities [prob] and percent
//invariant sites, [pinvar] and character set [l]
double logMALlikelihood( const struct ml *l, const double* ws, const double* pi, 
                        const double* prob, const double pinvar );
double logMPLlikelihood( const struct ml *l, const double* ws, const double* pi, 
                        const double* prob );


/**
 * functions on character sets
 */ 
void filter_s( const struct ml* c1, struct ml* c2, const int* idxs, const int n );
int compare_chars( const struct ml* c1, const struct ml* c2 );

#endif /* USE_LIKELIHOOD */
