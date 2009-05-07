/* POY 4.0 Beta. A phylogenetic analysis program using Dynamic Homologies.    *\
(* Copyright (C) 2007  Andrés Varón, Le Sy Vinh, Illya Bomash, Ward Wheeler,  *)
(* and the American Museum of Natural History.                                *)
(*                                                                            *)
(* This program is free software; you can redistribute it and/or modify       *)
(* it under the terms of the GNU General Public License as published by       *)
(* the Free Software Foundation; either version 2 of the License, or          *)
(* (at your option) any later version.                                        *)
(*                                                                            *)
(* This program is distributed in the hope that it will be useful,            *)
(* but WITHOUT ANY WARRANTY; without even the implied warranty of             *)
(* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the              *)
(* GNU General Public License for more details.                               *)
(*                                                                            *)
(* You should have received a copy of the GNU General Public License          *)
(* along with this program; if not, write to the Free Software                *)
(* Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301   *)
\* USA                                                                        */
       
#include <stdio.h>
#include <stdlib.h>         //malloc, calloc, srand, RAND_MAX
#include <string.h>         //memcpy, memset
#include <assert.h>         
#include <time.h>
    
#include "config.h"         //defines, if likelihood, USE_LIKELIHOOD
#ifdef USE_LIKELIHOOD   
#include <math.h>           //log10,exp

//caml specific headers
#include <caml/alloc.h>     //copy_double, et cetera
#include <caml/mlvalues.h>
#include <caml/memory.h>    //caml_param, et cetera
#include <caml/bigarray.h>
#include <caml/custom.h>    //abstract types
#include <caml/intext.h>    //serialization
#include <caml/fail.h>      //failwith('')
#include <caml/callback.h>
// hello worlk
#include "floatmatrix.h"
#include "likelihood.h"     //includes floatmatrix

//decent values (time/accuracy)
#define EPSILON 1e-4    //error for numerical calculations
#define MAX_ITER 10000  //number of iterations for numerical calculations
#define CAML_ALLOC_2 1000000
     
/* for unboxing abstract types */
#define FM_val(v) (*((struct matrix**)Data_custom_val(v)))
#define FM_ptr(v)   ((struct matrix**)Data_custom_val(v))
#define ML_ptr(v)   ((struct ml**)Data_custom_val(v))
#define ML_val(v) (*((struct ml**)Data_custom_val(v)))

/* ~CHECK_MEM       -- verify enough memory */
#define CHECK_MEM(a) if(a==NULL) failwith("I cannot allocate more memory.")
/* ~CHECK_ZEROES    -- verify array is all zeroes with EPSILON error */
#define CHECK_ZEROES(a,n); int Z;for(Z=0;Z<n;Z++){ if(a[Z] > EPSILON || a[Z] < -EPSILON)\
                            failwith("Imaginary eigenvalues"); }
//------------------------------------------------------------------------------
/** 
 * CONVENTIONS:
 *
 * variables that start with a:
 *        ...uppercase letter --> MATRICES
 *        ...lowercase letter --> VECTORS
 *
 * Xt is the transpose of X
 * Xi is the inverse of X
 *
 * D is a diagonal matrix of eigenvalues
 * U are eigenvectors (column major, Ut is rowmajor)
 * pi are the priors as a vector
 *
 * (int) variable:
 *        n is the number of columns/size of alphabet
 *        t is the branch length
 */

//------------------------------------------------------------------------------
/* prints a matrix (block format) */
void printmatrix( const double* Z, const int n, const int m)
{
    int i,j;
    for (i=0; i<m; ++i) {
        putchar('\t');
        for (j=0; j<n; ++j)
            printf("[%6.8f] ", Z[i*n+j]);
        putchar('\n'); 
    }
}
/* prints an array horizontally -- located in gamma.c now
void printarray( const double* a, const int n )
{
    int i;
    for(i=0;i<n;++i)
        printf("[%6.5f] ", a[i]);
    putchar('\n');
}
*/

/* print out a caml value */
void CAML_debug( value s )
{
    mll* a;
    a = ML_val( s );
    printf("DEBUG:\n\tN: %d\n\tL: %d", a->stride, a->c_len );
    printmatrix(a->lv_s,a->stride,a->c_len);
    putchar('\n');
}
/** creates a random substituation rate matrix
 * each row sums to 0 on diagonal elements
 */
void rand_sub_mat_gtr( double* A, const int n)
{
    srand( time(NULL) );
    int i,j;
    for(i=0;i<n;++i){
        double diag = 0;
        for(j=0;j<n;++j){
            if( i == j ){ continue; }
            A[i*n+j] = ((double)rand()/(double)RAND_MAX);
            diag = diag + A[i*n+j];
        }
        A[i*n+i] = -diag;
    }
}
/** transpose of matrix [A] of [n]x[n], must be square since no reallocation
 * happens 
 */
void transpose( double *A, const int n )
{
    double tmp;
    int i,j;
    for(i=0;i<n;++i){
        for(j=i+1;j<n;++j){
            tmp = A[i*n+j];
            A[i*n+j] = A[j*n+i];
            A[j*n+i] = tmp;
        }
    }
}
/** creates a random symmetric sub-rate matrix
 * each row sums to 0 on diagonal elements
 */
void rand_sub_mat_sym( double* A, const int n)
{
    srand( time(NULL) );
    int i,j;
    double temp;
    for(i=0;i<n;++i){
        double diag = 0;
        for(j=0;j<n;++j){
            if( i == j ){ continue; }
            if( i < j ) { diag = diag + A[i*n+j]; continue; }
            temp = ((double)rand()/(double)RAND_MAX);
            A[i*n+j] = temp;
            A[j*n+i] = temp;
            diag = diag + A[i*n+j];
        }
        A[i*n+i] = -diag;
    }
}

//------------------------------------------------------------------------------
// O'CAML interface functions for Garbage Collection + callback for a GC

/* release struct ml->lv_s array to pool */
void likelihood_CAML_free( value v )
{
    mll* val;
    val = ML_val(v);
    free (val->lv_s);
    free (val);
}

/* serialization of struct ml */
void
likelihood_CAML_serialize(value v,unsigned long* wsize_32,unsigned long* wsize_64)
{
    mll* s;
    s = ML_val(v);
    caml_serialize_int_4( s->stride );
    caml_serialize_int_4( s->c_len );
    int i,j;
    j = s->stride * s->c_len;
    for(i=0;i<j;++i)
        caml_serialize_float_8( s->lv_s[i] );
}

/* deserialization of struct ml */  
unsigned long likelihood_CAML_deserialize( void* v )
{
    mll* s;
    s = (mll*) malloc( sizeof(mll) );
    CHECK_MEM(s);
    s->stride = caml_deserialize_sint_4();
    s->c_len  = caml_deserialize_sint_4();
    int i,j=s->stride*s->c_len;
    double* stuff;
    stuff = (double*) malloc( sizeof(double)*j );
    CHECK_MEM(stuff);
    for(i=0;i<j;++i)
        stuff[i] = caml_deserialize_float_8();

    return (sizeof(struct ml **));
}

/* compares two character sets */
int compare_chars( const mll* c1, const mll* c2)
{
    int ret,i,n;
    ret = c1->stride - c2->stride;
    if (ret != 0)
        return ret;
    ret = c1->c_len - c2->c_len;
    if (ret != 0)
        return ret;
    n = c1->c_len * c1->stride;
    for( i = 0; i < n; ++i ){
        ret = c1->lv_s[i] - c2->lv_s[i];
        if (ret != 0)
            return ret;
    }
    return 0;
}
int likelihood_CAML_compare( value c1, value c2 )
{ return compare_chars( ML_val(c1), ML_val(c2) ); }

/* custom garbage collection */
static struct custom_operations likelihood_custom_operations  = {
    "http://www.amnh.org/poy/likelihood/likelihood.0.1", //identifier
    (&likelihood_CAML_free),        //finalize
    (&likelihood_CAML_compare),     //compare
    custom_hash_default,          //hash
    (&likelihood_CAML_serialize),   //serialize
    (&likelihood_CAML_deserialize)  //deserialize
};

/* registration with ocaml GC */
value 
likelihood_CAML_register (value u) {
    CAMLparam1(u);
    register_custom_operations (&likelihood_custom_operations);
    CAMLreturn (Val_unit);
}

/* calls a full garbage collection, a is used to ensure success. */
void likelihood_CAML_gc()
{
    int a; a = 0;
    a = Int_val( caml_callback(*caml_named_value("likelihood_gc_compact"), Val_int(a)) );
    assert( a == 0 );
}

value likelihood_CAML_StoBigarray( value s )
{
    CAMLparam1( s );
    CAMLlocal1( res );
    mll* work;
    long dims[2];
    double *newone;

    work = ML_val( s );
    dims[0] = work->c_len;
    dims[1] = work->stride;
    newone = (double*) malloc( dims[0] * dims[1] * sizeof(double));
    memcpy(newone,work->lv_s, dims[0] * dims[1] * sizeof(double));

    res = alloc_bigarray(BIGARRAY_FLOAT64 | BIGARRAY_C_LAYOUT,2,newone,dims);
    CAMLreturn( res );
}

/**
 *  convert bigarray.array2 to struct ml
 */
value likelihood_CAML_BigarraytoS( value A )
{
    CAMLparam1( A );
    CAMLlocal1( s );
    double *stuff,*o_stuff;
    mll* ret;

    o_stuff = (double*) Data_bigarray_val( A );

    ret = (mll*) malloc( sizeof(mll) );
    CHECK_MEM(ret);
    ret->c_len = Bigarray_val(A)->dim[0];
    ret->stride = Bigarray_val(A)->dim[1];

    s = caml_alloc_custom(&likelihood_custom_operations, (sizeof(mll*)),
            ret->c_len*ret->stride, CAML_ALLOC_2); 
    ML_val(s) = ret;

    stuff = (double*) malloc( ret->c_len * ret->stride * sizeof(double));
    CHECK_MEM(stuff);
    memcpy( stuff, o_stuff, ret->c_len * ret->stride * sizeof(double));
    ret->lv_s = stuff;

    assert( ML_val(s) == ret);
    CAMLreturn( s ); 
}

/**
 * filters an ml struct with all the indexes in the passed array --ordered
 */
value likelihood_CAML_filter(value as, value ibs)
{  
    CAMLparam2( as, ibs );
    CAMLlocal1( package );
    mll *norm, *ret;
    double* area;
    int m,i,j,n,k;

    m = Wosize_val( ibs );
    norm = ML_val( as );
    n = norm->c_len;
    area = (double*) malloc( sizeof(double)* (norm->c_len-m)*norm->stride );
    CHECK_MEM(area);

    for(i=0,j=0;i<n;++i){     
        if( i != Int_val(Field(ibs, j)) ){
            for(k=0;k<norm->stride;k++)
                area[i*norm->stride + k] = norm->lv_s[ i*norm->stride + k]; 
        } else if(j++ == m)
            break;
    }
    for(;i<m;++i)
        area[i] = norm->lv_s[i];

    package = caml_alloc_custom(&likelihood_custom_operations, (sizeof(mll*)),
                                                    (n-m)*norm->stride, CAML_ALLOC_2);
    ret = ML_val( package );
    ret->stride = norm->stride;
    ret->c_len = norm->c_len - m;

    assert( ret == ML_val(package));
    CAMLreturn( package );
}



/**  [mk_diag diag mat n m] ~ makes a diagonal square matrix from a first [n]
 * elements of the matrix. Put first [n] elements of [M] along diagonal of [M]
 */
#ifdef _WIN32
__inline void 
#else
inline void 
#endif
mk_diag(double* M, const int n)
{
    int i;
    for(i = 1; i < n; ++i){
        M[i*n+i] = M[i];
        M[i] = 0;
    }
}

/** [scale_VL VR VL n] scales VL so VLt*VR = I
 *
 * Since LAPACK doesn't guarentee that VL^t * VR = I, where VL are the left
 * eigen vectors and VR are the right eigenvectors, but that VLt * VR = D mod 1,
 * where D is a diagonal matrix. Thus we perform a multiplication, and grab the
 * diagonal elements as scale factors for VL (thus dividing through each row).
 * Conversly this can be done by scaling VR, by scaling though each column 
 * --more of a fortran method though. 
 *
 *  **output** is in VL
void scale_VL( const double *VR, double *VL, const int n)
{   NOT USED CURRENTLY --calculating the inverse is better.
    double *alphas,a,b;
    int i,j; 
    char ntran, _tran;

    ntran = 'N'; _tran = 'T';
    a = 1; b = 0;
    alphas = (double*) malloc(n*n*sizeof(double));
    CHECK_MEM( alphas );

    cblas_dgemm( CblasRowMajor, CblasTrans, CblasNoTrans,
            n, n, n, a, VL, n, VR, n, b, alphas, n );

    for(i = 0;i < n; i++){
        for(j=0;j<n;j++)
            VL[i*n+j] = VL[i*n+j] / alphas[i*n+i];
    }
    free(alphas);
} */

/** [mk_inverse VL D VR] makes the inverse of VL from the diagonalization
 *
 * Since lapack doesn't guarentee that VL = inv(VR) when we have duplicate
 * eigenvalues (this happens less often in GTR, but symmetric matrices have
 * an easy solution, Ut = Ui). Otherwise we just scale Ui's rows (scale_VL).
 *
 * output is in VL.
 */
int
mk_inverse(mat *space,double *VL, const double *D, const double *VR, int n, double *tmp)
{
    //LU factorization and inverse via DGETRF and DGETRI
    int i,lwork,*pivot;
    double work_size, *work;

    lwork = -1;
    pivot = (int*) malloc(n * sizeof(double));
    memcpy(VL, VR, n*n*sizeof(double) ); //TODO: is this a problem?
    dgetrf_(&n, &n, VL, &n, pivot, &i);
    if( 0 == i ){ //hello workd
        dgetri_(&n, VL, &n, pivot, &work_size, &lwork, &i); //optimal work
        if( 0 == i ){
            lwork = (int)work_size;
            work = register_section( space, lwork, 1 );
            dgetri_(&n, VL, &n, pivot, work, &lwork, &i);
            free( work );
        } else if (i < 0) { 
            failwith ( "dgetri_ argument failed." );
        } else {
            failwith ( "dgetri_ matrix is singular and inverse cannot be computed." );
        }
    } else {
        failwith ( "degetri_ unknown error" );
    }
    return i;
}

/**  [apply_exp diag n m t]
 * Multiplies the diagonal of [diag], an [n]x[m] matrix by [t] and applies exp()
 */
#ifdef _WIN32
__inline void 
#else
inline void 
#endif
apply_exp(double* D, const int n, const int m, const double t)
{
    int i,s;
    s = m * n;
    for(i=0;i<s;i=i+n+1 )
        D[i] = exp( D[i] * t );
}

/** [proportion a b] - number of SIMILAR sites **/
double proportion( const mll* a, const mll* b)
{ 
    double prop, maxa, maxb;
    int i /*current char*/  ,j /*current state*/,
        k /*array location*/,t /*total set state in char*/,
        s /* characters that are equal in a set */;
    //length of characters is equal
    assert(a->stride == b->stride);
    assert(a->c_len == b->c_len);

    prop = k = 0;
    //i holds current char, j holds current array location, k holds curr char
    for(i=0;i<a->c_len;i++){
        maxa = maxb = 0;
        s = t = 0;
        for(j=0;j<a->stride;j++){
            if (a->lv_s[k+j] > maxa){ maxa = a-> lv_s[k+j]; }
            if (b->lv_s[k+j] > maxb){ maxb = b-> lv_s[k+j]; }
        }
        for(j=0;j<a->stride;j++){
            if(a->lv_s[k] == maxa && b->lv_s[k] == maxb){ s++; t+=2; }
            else if (a->lv_s[k] == maxa || b->lv_s[k] == maxb){ t++; }
            k++;
        }
        prop += (2*s)/t;
    }
    prop = prop / a->c_len;
    return prop;
}

value likelihood_CAML_proportion( value a, value b )
{
    CAMLparam2(a,b);
    CAMLlocal1( prop );
    double p;
    mll *a_ml, *b_ml;
    a_ml = ML_val( a );
    b_ml = ML_val( b) ;
    p = proportion( a_ml , b_ml );
    prop = caml_copy_double( p );
    CAMLreturn( prop );
}


/**  [diagonalize_sym A D N]
 * Diagonalizes [A] that is [N]x[N] and upper symmetric and places
 * the resultant eigenvalues in the matrix [D] along the diagonal.
 * Overwrites [A] with the eigenvectors (rowmajor).
 *
 * @returns =0 on success
 *          <0 ith argument had an issue
 *          >0 failed in convergence
 */
int diagonalize_sym( mat *space,double* A, double* D, int n)
{
    char jobz = 'V';
    char uplo = 'U';
    int info = 0,lwork=-1;
    double *work,work_size;

    //find the optimal work (call func with -1 in lwork)
    //result ends up in work_size
    dsyev_(&jobz, &uplo, &n, A, &n, D, &work_size, &lwork, &info);
    if( info == 0 ){
        lwork = (int)work_size;
        work = (double*) register_section (space,lwork,1);
        /** dsyev - calculate eigenvalues and vectors
         *          possibly use PDSYGVX or PDSYEVX --less precision)

         * jobz   - V to return eigenvectors (N otherwise)
         * uplo   - Upper/Lower triangular? (U/L)
         * n      - order of the matrix
         * E_vecs - (in/out) dim(lda,n) --note: lda==n
         * n      - lda
         * e_vals - (out) eigenvalues
         * work   - (workspace/out) = dim(n)
         * lwork  - length of work >= (NB+2)*N
         * info   - retval; see comments above */
        dsyev_(&jobz, &uplo, &n, A, &n, D, work, &lwork, &info);
        if( 0 == info ){
            mk_diag(D,n);
        } else if ( info < 0 ){
            failwith ("dsyev_ argument failed");
        } else { 
            failwith ("dsyev_ diagonalization failed to converge. Singular matrix?");
        }
    }
    return info;
}
/* Diagonlize matrix interfaces: symmetric */
void likelihood_CAML_diagonalize_sym( value tmp, value Q, value D)
{
    CAMLparam3( tmp, Q, D );
    int n;
    n = Bigarray_val ( Q )->dim[0];
    diagonalize_sym( FM_val( tmp ), (double*) Data_bigarray_val( Q ),
                        (double*) Data_bigarray_val( D ), n);
    free_all( FM_val(tmp) );
    CAMLreturn0;
}

int diagonalize_gtr(mat *space, double* A, double* D, double* Ui, int n)
{
    char jobv_ = 'V'; //we went to find left and right eigenvectors
    double *wi,*U,*work,work_size;
    int lwork,info;
    wi = register_section( space, n, 0 );
    U  = register_section( space, n*n, 0 );

    //find the optimal work (call func with -1 in lwork)
    lwork = -1;
    //D holds the real values eigen values through the computation
    dgeev_(&jobv_,&jobv_,&n,A,&n,D,wi,U,&n,Ui,&n,&work_size,&lwork,&info);
    if( info == 0 ) {
        lwork = (int)work_size;
        work = register_section( space, lwork, 1);
        /** dgeev   - A * v(j) = lambda(j) * v(j)
         *            u(j)**H * A = lambda(j) * u(j)**H
         *            where:
         *             ` **H is conjugate transpose)
         *             ` u(j) is a left eigenvector
         *             ` v(j) is a right eigenvector
         *             ` lambda(j) is an eigenvalue

         * JOBVL    -calulate left eigenvectors? V/N
         * JOBVR    -calulate right eigenvectors? V/N
         * N        - order(A)
         * A        - matrix of LDAxN           (**MODIFIED**)
         * LDA      - dim(A) 
         * WR       - real parts of eigenvalues (**OUTPUT**) //putting them in D
         * WI       - imaginary parts           (**OUTPUT**) //ignored
         * VL       - left eigenvectors         (**OUTPUT**)
         * LDVL     - dim(VL)
         * VR       - right eigenvectors        (**OUTPUT**)
         * LDVR     - dim(VR)
         * WORK     - workspace
         * LWORK    - dim(work)
         * INFO     - =0 if successful
         *            <0 if the ith argument is illegal
         *            >0 if QR failed; i+1 is first eigenvalue */
        dgeev_(&jobv_,&jobv_,&n,A,&n,D,wi,U,&n,Ui,&n,work,&lwork,&info);
        //imaginary eigenvals should all = 0 --since the Matrix is similar to 
        //some symmetric matrix, thus have same eigenvalues (Keilson 1979)
        CHECK_ZEROES(wi,n);
        if( 0 == info ) {
            mk_diag(D,n);
            mk_inverse(space, U, D, Ui, n, wi); /* wi is extra */
        } else if ( info < 0 ){
            failwith ("dgeev_ argument failed");
        } else {
            failwith ("dgeev_ QR failed, possibly singular matrix?");
        }
    }
    memcpy(A,U,n*n*sizeof(double));
    return info;
}
/* Diagonalize matrix interface: general */
void likelihood_CAML_diagonalize_gtr( value tmp, value Q, value D, value Qi)
{
    CAMLparam4( tmp, Q, D, Qi );
    diagonalize_gtr(
            FM_val( tmp ),
            (double*) Data_bigarray_val( Q ),
            (double*) Data_bigarray_val( D ),
            (double*) Data_bigarray_val( Qi),
            Bigarray_val( Q )->dim[0]);
    free_all( FM_val(tmp) );
    CAMLreturn0;
}

/** [mk_probmat_*** P U D [Ui] t]
 * Finds the probability matrix based on the diagonalized substitution
 * matrix, [U] and [D], with branch length [t]. returns result in [P].
 *
 * UNCHECKED EXCEPTION ::
 * [n]x[n] == dim([D]) == dim([P]) == dim([U]) == dim([Ui]) == dim([TMP])
 */  
void
compose_sym(double* P,const double* U,const double* D,const float t,int n,double *TMP)
{
    char _tran,ntran;
    double alpha,beta;

    alpha = 1; beta = 0; _tran = 'T'; ntran = 'N';
    memcpy(P, D, n*n*sizeof(double) );
    apply_exp(P,n,n,t); //exp(D*t); along diagonal only
    //calculates: C = op(A)*op(B)*a + C*b
    dgemm_(&ntran,&ntran,        //format, op(A), op(B)
            &n, &n, &n, &alpha,  //rows A, cols B, cols A, multiplier of (A*B)
            U, &n,               //MATRIX A, stride for A
            P, &n,               //MATRIX B, stride for B
            &beta, TMP, &n );    //multiplier for C, MATRIX C, stride for C
    //if scalor mult of C == 0, C isn't used to add
    dgemm_(&ntran,&_tran,&n,&n,&n,&alpha,TMP,&n,U,&n,&beta,P,&n);
}
value likelihood_CAML_compose_sym(value tmp,value U, value D, value t)
{   /* only called for testing purposes */
    CAMLparam4( tmp,U,D,t );
    CAMLlocal1( res );
    double *c_P, c_t,*c_D,*c_U,*c_T;
    mat *space;
    int n;
    long dims[2];

    space = FM_val (tmp);
    n = Bigarray_val( U )->dim[0];
    c_t = Double_val( t );
    c_U = (double *) Data_bigarray_val( U );
    c_D = (double *) Data_bigarray_val( D );
    c_P = (double*) malloc( n*n*sizeof(double));
    CHECK_MEM(c_P);
    c_T = register_section( space, n*n, 0 );

    compose_sym(c_P,c_U,c_D,c_t,n,c_T);
    dims[0] = n; dims[1] = n;
    res = alloc_bigarray(BIGARRAY_FLOAT64 | BIGARRAY_C_LAYOUT, 2, c_P, dims);
    free_all( space );
    CAMLreturn ( res );
}

void
compose_gtr(double* P, const double* U, const double* D, const double* Ui, 
        const double t, const int n,double *tmp)
{
    double alpha,beta; char ntran;
    alpha = 1; beta = 0; ntran = 'N';

    memcpy(P,D,n*n*sizeof(double));
    apply_exp(P,n,n,t);
    //S is U*exp(D)...
    dgemm_(&ntran,&ntran,&n,&n,&n,&alpha,Ui,&n,P,&n,&beta,tmp,&n);
    //P becomes U*expD*Ui... done --note: VL = inv(VR)
    dgemm_(&ntran,&ntran,&n,&n,&n,&alpha,tmp,&n,U,&n,&beta,P,&n);
}
value 
likelihood_CAML_compose_gtr(value tmp,value U, value D, value Ui, value t)
{
    CAMLparam5( tmp,U,D,Ui,t );
    CAMLlocal1( res );
    double *c_P, c_t,*c_D,*c_U,*c_Ui,*c_T;
    int n;
    mat *space;
    long dims[2];

    space = FM_val( tmp );
    n = Bigarray_val( U )->dim[0];
    c_t = Double_val( t );
    c_U = (double*) Data_bigarray_val( U );
    c_D = (double*) Data_bigarray_val( D );
    c_Ui= (double*) Data_bigarray_val( Ui);
    c_P = (double*) malloc( n*n*sizeof(double) );
    CHECK_MEM(c_P);
    c_T = register_section( space, n*n, 0 );

    compose_gtr(c_P,c_U,c_D,c_Ui,c_t,n,c_T);
    dims[0] = n; dims[1] = n;
    res = alloc_bigarray(BIGARRAY_FLOAT64 | BIGARRAY_C_LAYOUT, 2, c_P, dims);

    free_all(space);
    CAMLreturn( res );
}

/**  [median_h P l c nl a]
 * Finds half the likelihood vector of child [l] and probability matrix [P] into [nl]
 *
 *  [a] - length of each vector (in practice, the alphabet size) [P] is [a]x[a]
 *  [c] - the start of the character to do work on
 */
void median_h( const double* P, const double* l, const int c, double* nl, const int a)
{
    int i,j; double elm;
    //for each row of P
    for(i=0; i<a; ++i){
        elm = 0;
        for(j=0;j<a;++j)
            elm += (P[(i*a)+j] * l[c+j]);
        nl[i] = elm;
    }
}

/** [loglikelihood ml p]
 * returns loglikelihood of a character set. |p| = l->stride
 */
double loglikelihood( const mll* l, const double* p )
{
    int i, j, nchars, size; 
    double ret,tmp;
    nchars = l->c_len * l->stride;  //number of characters
    size = l->stride;   //length of vectors
    ret=0;tmp=0;

    for(i=0,j=0; i<nchars; ++i,++j){
        if( j == size ){
            ret -= log( tmp );
            j = 0; tmp = 0;
        }
        tmp += l->lv_s[i] * p[j];
    }
    ret -= log (tmp);
    return ret;
}
/* [likelihoood_CAML_loglikelihood s p] wrapper for loglikelihood */
value likelihood_CAML_loglikelihood(value s, value p)
{
    CAMLparam2( s, p );
    CAMLlocal1( mle );
    double *c_p,ll;
    mll* chars;
    c_p = (double*) Data_bigarray_val( p );
    chars = ML_val( s );

    ll = loglikelihood( chars, c_p );
    mle = caml_copy_double( ll );
    CAMLreturn( mle );
}


/** [median_c Pa Pb a b c]
 * Finds the likelihood vector [c] based on vectors of its children [a] and
 * [b] with the probability matrices [Pa] and [Pb], respectively, and
 * probability [p].
 *
 * tmp are size of alphabet.
 */
void
median_charset(const double* Pa,const double* Pb,
        const mll* a,const mll* b, mll* c,const double p,
        double* tmp1, double* tmp2)
{
    int i,j,len;
    assert( a->stride == b->stride );
    assert( b->c_len == a->c_len );

    for(i=0; i < a->c_len; ++i){
        len = i*a->stride; 
        /* find half medians for each character */
        median_h( Pa, a->lv_s, len, tmp1, a->stride );
        median_h( Pb, b->lv_s, len, tmp2, a->stride );
        /* pairwise multiplication */
        for(j=0;j < a->stride;++j)
            c->lv_s[len+j] += (tmp2[j]*tmp1[j])*p;
    }
    c->c_len = a->c_len;
    c->stride = a->stride;
}

value likelihood_CAML_median_wrapped_sym
    (value tmp,value U,value D,value ta,value tb,value ml_a,value ml_b,value rates,value probs)
{
    /* ocaml macros */
    CAMLparam5( tmp,U,D,ta,tb );
    CAMLxparam4( ml_a,ml_b,rates,probs );
    CAMLlocal1( ml_c );
    /* all variables */
    double cta,ctb,*c_U,*c_D,*PA,*PB,*tmp1,*g_rs,*p_rs;
    mll *a,*b,*c;
    mat *space;
    int num_rates,num_probs,i;
    /* probabilities/rates */
    space = FM_val( tmp );
    num_rates = Bigarray_val(rates)->dim[0];
    num_probs = Bigarray_val(probs)->dim[0];
    assert( num_rates == num_probs );
    g_rs = (double*) Data_bigarray_val( rates );
    p_rs = (double*) Data_bigarray_val( probs );
    /* diagonalized transition matrix */
    c_U = (double*) Data_bigarray_val( U );
    c_D = (double*) Data_bigarray_val( D );
    cta = Double_val( ta );
    ctb = Double_val( tb );
    assert( cta > 0 && ctb > 0 );
    /* character sets */
    a = ML_val( ml_a );
    b = ML_val( ml_b );
    assert( a->stride == b->stride );
    assert( a->c_len == b->c_len );
    /* create new character set, c */
    ml_c = caml_alloc_custom(&likelihood_custom_operations,(sizeof(mll*)),a->stride*a->c_len,CAML_ALLOC_2);
    c = (mll*) malloc( sizeof(mll) );
    CHECK_MEM(c);
    ML_val( ml_c ) = c;
    c->stride = a->stride;
    c->c_len = a->c_len;
    c->lv_s = (double*) calloc( c->c_len * c->stride, sizeof(double));
    CHECK_MEM(c->lv_s);
    /* register temp variables */
    PA = (double*) register_section(space, a->stride*a->stride, 1);
    PB = (double*) register_section(space, a->stride*a->stride, 1);
    tmp1 = (double*) register_section(space, a->stride*a->stride, 1);
    /** 
     *  main loop --
     *
     *   ___ num_rates
     *   \
     *    >   pr_i * F(x | r_i,P )
     *   /__ 
     *    i
     *        where:  F(x|r,P) = f_right(x*r|P) * f_left(x*r|P)
     *              f_j(x*r|P) = SUM(P_jk,k:{ACTG})
     *
     *              P    = probability matrix
     *              pr_i = the rate probability
     *              r_i  = the rate of i
     *              x    = the branch length
     *              P_j  = is the jth row of P
     */
    for(i=0;i<num_rates;i++){
        compose_sym( PA, c_U, c_D, cta*g_rs[i], a->stride,tmp1 );
        compose_sym( PB, c_U, c_D, ctb*g_rs[i], b->stride,tmp1 );
        median_charset( PA, PB, a,b,c,p_rs[i],tmp1,&(tmp1[a->stride]));
    }
    /* free up space, return */
    free_all( space );
    assert( c == ML_val(ml_c));
    assert( c->c_len == b->c_len );
    assert( c->stride == a->stride );
    CAMLreturn(ml_c);
}

value likelihood_CAML_median_wrapped_gtr
    (value tmp,value U,value D,value Ui,value ta,value tb,value ml_a,value ml_b,
        value rates,value probs)
{
    /* ocaml macros */
    CAMLparam5( tmp,U,D,Ui,ta );
    CAMLxparam5( tb,ml_a,ml_b,rates,probs );
    CAMLlocal1( ml_c );
    /* declare variables */
    double cta,ctb,*c_U,*c_D,*c_Ui,*PA,*PB,*tmp1,*g_rs,*p_rs;
    mll *a,*b,*c;
    mat *space;
    int num_rates,num_probs,i;
    space = FM_val (tmp);
    /* rates and probability vectors */
    num_rates = Bigarray_val(rates)->dim[0];
    num_probs = Bigarray_val(probs)->dim[0];
    assert( num_rates == num_probs );
    g_rs= (double*) Data_bigarray_val(rates);
    p_rs= (double*) Data_bigarray_val(probs);
    /* diagonalized transition matrix */
    c_U = (double*) Data_bigarray_val( U );
    c_D = (double*) Data_bigarray_val( D );
    c_Ui= (double*) Data_bigarray_val( Ui);
    /* branch lengths */
    cta = Double_val( ta );
    ctb = Double_val( tb );
    assert( cta > 0 && ctb > 0 );
    /* character set */
    a = ML_val( ml_a );
    b = ML_val( ml_b );
    assert(a->stride == b->stride);
    assert(b->c_len == a->c_len);
    /* create new character set */
    ml_c = caml_alloc_custom(&likelihood_custom_operations, (sizeof(mll*)), 
                                                    a->stride*a->c_len, CAML_ALLOC_2);
    c = (mll*) malloc( sizeof(mll) );
    CHECK_MEM(c);
    ML_val( ml_c ) = c;
    c->lv_s = calloc(a->c_len * a->stride, sizeof(double)); 
    /* register temp variables */
    PA = register_section( space, b->stride * b->stride, 1 );
    PB = register_section( space, b->stride * b->stride, 1 );
    tmp1 = register_section( space, b->stride * b->stride, 1 );
    /* main loop, see symmetric for description */
    for(i=0;i<num_rates;i++){
        compose_gtr( PA, c_U, c_D, c_Ui, cta*g_rs[i], a->stride, tmp1);
        compose_gtr( PB, c_U, c_D, c_Ui, ctb*g_rs[i], b->stride, tmp1);
        median_charset( PA, PB, a,b,c, p_rs[i],tmp1,&(tmp1[b->stride]));
    }
    /* free all variables */
    free_all( space );
    assert( c == ML_val(ml_c) );
    CAMLreturn(ml_c);
}

/* [likelihood_CAML_median_sym ,,,] argument wrapper for median_sym */
value likelihood_CAML_median_sym(value * argv, int argn)
{
    return likelihood_CAML_median_wrapped_sym
        ( argv[0],argv[1],argv[2],argv[3],argv[4],argv[5],argv[6],argv[7],argv[8] ); 
}

/* [likelihood_CAML_median_gtr ...] argument wrapper for median_gtr */
value likelihood_CAML_median_gtr(value * argv, int argn)
{
    return likelihood_CAML_median_wrapped_gtr
        ( argv[0],argv[1],argv[2],argv[3],argv[4],argv[5],argv[6],argv[7],argv[8],argv[9] );
}


//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
void single_sym(ptr *simp,double *PA,double *PB,const double *U,const double *D,const mll *a,
    const mll *b,const double time,const double *gam,const double *prob,
    const double *pi,const int g_n,double *tmp)
{
    int i; double half;
    assert( a->stride == b->stride );
    assert( a->c_len == b->c_len );
    half = time / 2;
    mll *c; c = simp->vs;
    simp->time = time;
    memset(c->lv_s,0,a->c_len*a->stride*sizeof(double));
    for(i=0;i<g_n;i++){
        compose_sym( PA, U, D, half*gam[i], a->stride,tmp );
        compose_sym( PB, U, D, half*gam[i], b->stride,tmp );
        median_charset( PA, PB, a,b,c,prob[i],tmp,&(tmp[a->stride]));
    }
    simp->ll = loglikelihood( c,pi );
}
void single_gtr(ptr *simp,double *PA,double *PB,const double *U,const double *D,const double *Ui,
    const mll *a,const mll *b,const double time,const double *gam,const double *prob,
    const double *pi,const int g_n,double *tmp)
{
    int i;
    double half;
    mll *c; c = simp->vs;
    simp->time = time;
    half = time / 2;
    memset(c->lv_s,0,a->c_len*a->stride*sizeof(double));
    for(i=0;i<g_n;i++){
        compose_gtr( PA, U, D, Ui, half*gam[i], a->stride, tmp );
        compose_gtr( PB, U, D, Ui, half*gam[i], b->stride, tmp );
        median_charset( PA, PB, a,b,c,prob[i],tmp,&(tmp[a->stride]));
    }
    simp->ll = loglikelihood( c,pi );
}

#define golden_exterior(a,b) a->time + ((b->time - a->time)*2.0/(sqrt(5.0)-1))

void
readjust_brents_sym(mat *space,const double* U,const double* D,const mll* a,
        const mll* b,mll* c,double* b_ta,double* b_tb,double* b_mle,
        const double* rates, const double *prob, const int g_n,const double* pi)
{
    /* variables */
    ptr *left,*right,*middle,*best,*temp;
    double *PA,*PB,*TMP,ttime1,ttime2,deno,numr;
    int max_iter,size,bracked;
    void *temp_ptr;
    max_iter = MAX_ITER;
    size = a->c_len*a->stride;
    /* register temp space and matrices */
    PA = register_section(space, a->stride*a->stride , 0);
    PB = register_section(space, a->stride*a->stride , 0);
    TMP = register_section(space, a->stride*a->stride , 0);
    /* set up the points */
    best = (ptr*) malloc( sizeof(ptr) );    middle = (ptr*) malloc(sizeof (ptr));
    best->vs = (mll*) malloc(sizeof(mll));  middle->vs = (mll*) malloc(sizeof(mll));
    left = (ptr*) malloc( sizeof (ptr) );   right = (ptr*) malloc(sizeof (ptr));
    temp = (ptr*) malloc( sizeof (ptr) );   temp->vs = (mll*) malloc (sizeof(mll));
    left->vs = (mll*) malloc(sizeof(mll));  right->vs = (mll*) malloc(sizeof(mll));
    /* register space */
    (best->vs)->lv_s = register_section(space,size,0);
    (middle->vs)->lv_s = register_section(space,size,0);
    (left->vs)->lv_s = register_section(space,size,0);
    (right->vs)->lv_s = register_section(space,size,0);
    (temp->vs)->lv_s = register_section(space,size,0);
    /* set initial times and constants */
    best->time = *b_ta + *b_tb;
    middle->time = *b_ta + *b_tb;
    left->time = middle->time / 10;
    right->time = middle->time * 2;
    best->ll = *b_mle;
    middle->ll = *b_mle;
    // set up best
    memcpy( (best->vs)->lv_s, c->lv_s, size * sizeof(double) );
    memcpy( (middle->vs)->lv_s, c->lv_s, size * sizeof(double) );
    // fill in initial data
    single_sym(left,PA,PB,U,D,a,b,left->time,rates,prob,pi,g_n,TMP);
    single_sym(right,PA,PB,U,D,a,b,right->time,rates,prob,pi,g_n,TMP);
    /* bracket minimum */
    bracked = 1; // true, since we are going to expect the best
    while( (max_iter-- != 0) || !(left->ll > middle->ll && middle->ll < right->ll)){
        if (left->ll <= middle->ll && middle->ll <= right->ll){
            // increasing,  move all points down
            temp_ptr = right;
            right = middle;
            middle = left;
            left = temp_ptr;
            ttime1 = golden_exterior( middle, right );
            if (ttime1 < EPSILON) { ttime1 = EPSILON; break; }
            single_sym(left,PA,PB,U,D,a,b,ttime1,rates,prob,pi,g_n,TMP);
        } else if ( left->ll >= middle->ll && middle->ll >= right->ll) {
            // decreasing,  move all points up (middle,right,new)
            temp_ptr = left;
            left = middle;
            middle = right;
            right = temp_ptr;
            ttime1 = golden_exterior( left, middle );
            single_sym(right,PA,PB,U,D,a,b,ttime1,rates,prob,pi,g_n,TMP);
        } else {
            printf("Curvature?:%f/%f -- %f/%f -- %f/%f\n",left->time,
                    left->ll,middle->time,middle->ll,right->time,right->ll);
            bracked = 0;
            if (left->ll < right->ll){ //new left middle
                temp_ptr = right;
                right = middle;
                middle = left;
                left = temp_ptr;
                ttime1 = golden_exterior (middle, right);
                if (ttime1 < EPSILON) { ttime1 = EPSILON; break; }
                single_sym(left,PA,PB,U,D,a,b,ttime1,rates,prob,pi,g_n,TMP);
            } else { //middle right new
                temp_ptr = left;
                left = middle;
                middle = right;
                right = temp_ptr;
                ttime1 = golden_exterior( left, middle );
                single_sym(right,PA,PB,U,D,a,b,ttime1,rates,prob,pi,g_n,TMP);
            }
            break;
        }
    }
    
    /* parabolic interpolation */
    while( bracked && ( ((max_iter--) != 0 ) ||
                      (fabs (middle->time - left->time) <= EPSILON) ||
                      (fabs (right->time - middle->time) <= EPSILON) ) ){
        printf ("Iteration: %d\n",max_iter);
        // calculate absissca
        ttime1 = (middle->time - left->time) * (middle->ll - right->ll);
        ttime2 = (middle->time - right->time) * (middle->ll - left->ll);
        numr = (ttime1 * (middle->time - left->time)) - (ttime2 * (middle->time - right->time));
        deno = (ttime1 - ttime2) * 2.0;
        if (deno < 0.0){ //points are colinear
            numr = -numr;
            deno = -deno;
        }
        if (deno <= EPSILON){ //points are colinear
            break; //TODO : something better then break
        } else {
            ttime1 = middle->time - (numr/deno);
            if (ttime1 <= 0.0)
                ttime1 = EPSILON;
        }
        // ttime1 is the next time to try 
        single_sym(temp,PA,PB,U,D,a,b,ttime1,rates,prob,pi,g_n,TMP);
        // choose best
        if (best->ll > temp->ll){
            best->ll = temp->ll;
            best->time = temp->time;
            memcpy( (best->vs)->lv_s, (temp->vs)->lv_s, size * sizeof(double) );
        }
        // reorder points
        if(left->time < temp->time && temp->time < middle->time){ //left new middle
            temp_ptr = right;
            right = middle;
            middle = temp;
            temp = temp_ptr;
        } else if (temp->time < left->time){ //new left middle
            temp_ptr = right;
            right = middle;
            middle = left;
            left = temp;
            temp = temp_ptr;
        } else if (middle->time < temp->time && temp->time < right->time){ //middle new right
            temp_ptr = left;
            left = middle;
            middle = temp;
            temp = temp_ptr;
        } else if (right->time < temp->time){ //middle right new
            temp_ptr = left;
            left = middle;
            middle = right;
            right = temp;
            temp = temp_ptr;
        }
    }
    memcpy( c->lv_s, (best->vs)->lv_s, size * sizeof(double));

    *b_ta = best->time / 2;
    *b_tb = best->time / 2;
    *b_mle= best->ll;

    free(best->vs); free(middle->vs); free(left->vs); free(right->vs); free(temp->vs);
    free(best); free(middle); free(left); free(right); free(temp);
}

void
readjust_brents_gtr(mat * space,const double* U,const double* D,const double* Ui,
        const mll* a, const mll* b, mll* c,double* b_ta,double* b_tb,double* b_mle,
        const double *rates,const double *prob,const int g_n,const double* pi)
{
}

//-----------------------------------------------------------------------------
/** [..._readjust_gtr... tmp U D Ui A B C ta tb gammas probs pis ll] **/
value
likelihood_CAML_readjust_gtr_wrapped
    (value tmp,value U,value D,value Ui,value A,value B,value C,value ta,value tb,
        value tc, value g,value p,value pis,value ll)
{   
    /* ocaml macros */
    CAMLparam5(U,D,Ui,A,B);
    CAMLxparam5(C,ta,tb,g,p);
    CAMLxparam3(ll,pis,tmp);
    CAMLlocal1( res );
    /* declare variables */
    double cta,ctb,ctc,likelihood;
    /* branch lengths */
    cta = Double_val( ta );
    ctb = Double_val( tb );
    ctc = Double_val( tc );
    likelihood = Double_val( ll );
    /* check rates, probabilities and priors */
    assert( (ML_val(A)->stride) == Bigarray_val(pis)->dim[0] );
    assert( Bigarray_val(g)->dim[0] == Bigarray_val(p)->dim[0] );
    /* readjust */
    //printf("%f\t%f\t%f\n", cta,ctb,likelihood);
    readjust_brents_gtr(
            FM_val(tmp),
            Data_bigarray_val( U ),
            Data_bigarray_val( D ),
            Data_bigarray_val( Ui ),
            ML_val( A ), ML_val( B ), ML_val( C ),
            &cta, &ctb, &likelihood,
            Data_bigarray_val(g),
            Data_bigarray_val(p),
            Bigarray_val( g )->dim[0],
            Data_bigarray_val(pis)
    );
    //printf("%f\t%f\t%f\n", cta,ctb,likelihood);
    /* construct tuple / freespace */
    res = caml_alloc_tuple( 3 );
    Store_field(res, 0, caml_copy_double(cta));
    Store_field(res, 1, caml_copy_double(ctb));
    Store_field(res, 2, caml_copy_double(likelihood));
    free_all( FM_val(tmp) );
    CAMLreturn(res);
}
value
likelihood_CAML_readjust_sym_wrapped
    (value tmp,value U,value D,value A,value B,value C,value ta,value tb,value tc,value g,value p,value pis,value ll)
{
    /* ocaml macros for GC */
    CAMLparam5(U,D,ll,A,B);
    CAMLxparam5(C,ta,tb,g,p);
    CAMLxparam2( pis,tmp );
    CAMLlocal1( res );
    /* variables in function */
    double cta,ctb,ctc,likelihood;
    /* branch lengths // likelihood */
    cta = Double_val( ta );
    ctb = Double_val( tb );
    ctc = Double_val( tc );
    likelihood = Double_val( ll );
    /* check probabilities, rates, priors */
    assert( ML_val(A)->stride == Bigarray_val( pis )->dim[0]);
    assert( Bigarray_val( g )->dim[0] == Bigarray_val( p )->dim[0] );
    /* readjust loop */
    printf("s: %f\t%f\t%f\n", cta,ctb,likelihood);
    readjust_brents_sym(
            FM_val(tmp),
            Data_bigarray_val( U ), Data_bigarray_val( D ),
            ML_val(A), ML_val(B), ML_val(C),
            &cta, &ctb, &likelihood,
            Data_bigarray_val( g ),
            Data_bigarray_val( p ),
            Bigarray_val( g )->dim[0],
            Data_bigarray_val( pis ) );
    printf("e: %f\t%f\t%f\n", cta,ctb,likelihood);
    /* construct tuple / free space */
    res = caml_alloc_tuple( 3 );
    Store_field(res, 0, caml_copy_double (cta));
    Store_field(res, 1, caml_copy_double (ctb));
    Store_field(res, 2, caml_copy_double (likelihood));
    free_all( FM_val(tmp) );
    CAMLreturn(res);
}

/** functions for interface **/
value
likelihood_CAML_readjust_gtr(value * argv, int argn)
{
    return likelihood_CAML_readjust_gtr_wrapped
        (argv[0], argv[1],argv[2],argv[3],argv[4],argv[5],argv[6],argv[7],
         argv[8],argv[9],argv[10],argv[11],argv[12],argv[13]);
}
value
likelihood_CAML_readjust_sym(value * argv, int argn)
{
    return likelihood_CAML_readjust_sym_wrapped
        (argv[0], argv[1],argv[2],argv[3],argv[4],argv[5],argv[6],
         argv[7],argv[8],argv[9],argv[10],argv[11],argv[12]);
}

#endif /* USE_LIKELIHOOD */
