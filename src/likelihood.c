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
 
#define MAX(a,b) (a >= b)?a:b
#define MIN(a,b) (a <= b)?a:b

/* ~CHECK_MEM       -- verify memory was returned properly */
#define CHECK_MEM(a) if(a==NULL) failwith("I cannot allocate more memory.")
/* ~CHECK_ZEROES    -- verify array is all zeroes with EPSILON error */
/*                     used to ensure no imaginary roots when diagonalizing */
#define CHECK_ZEROES(a,n); int Z;for(Z=0;Z<n;Z++){ if(a[Z] > EPSILON || a[Z] < -EPSILON)\
                            failwith("Imaginary eigenvalues"); }

/* for unboxing abstract types - get pointer or actual value */
#define FM_val(v) (*((struct matrix**)Data_custom_val(v)))
#define FM_ptr(v)   ((struct matrix**)Data_custom_val(v))
#define ML_ptr(v)   ((struct ml**)Data_custom_val(v))
#define ML_val(v) (*((struct ml**)Data_custom_val(v)))
/* extra boxing/unboxing for option types */
#define Val_none Val_int(0)
#define Some_val(v) Field(v,0)

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
void printmatrixf( const double* Z, const int n, const int m)
{
    int i,j;
    for (i=0; i<m; ++i) {
        putchar('\t');
        for (j=0; j<n; ++j)
            printf("[%11.10f] ", Z[i*n+j]);
        putchar('\n'); 
    }
}
void printmatrixi( const int* Z, const int n, const int m)
{
    int i,j;
    for (i=0; i<m; ++i) {
        putchar('\t');
        for (j=0; j<n; ++j)
            printf("[%d] ", Z[i*n+j]);
        putchar('\n'); 
    }
}
/* prints an array horizontally */
void printarrayf( const double* a, const int n )
{
    int i;
    for(i=0;i<n;++i)
        printf("[%6.5f] ", a[i]);
    putchar('\n');
}
/* prints an array horizontally */
void printarrayi( const int* a, const int n )
{
    int i;
    for(i=0;i<n;++i)
        printf("[%d] ", a[i]);
    putchar('\n');
}

/* print out a caml value */
void CAML_debug( value s )
{
    mll* a;
    a = ML_val( s );
    printf("DEBUG:\n\tR: %d\tN: %d\n\tL: %d\n",
                a->rates,a->stride,a->c_len);
    printmatrixf(a->lv_s,a->stride,a->c_len*a->rates);
    printarrayi(a->lv_invar,a->c_len);
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
 * happens */
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
    if ( 1 == val->invar)
        free (val->lv_invar);
    free (val);
}

/* serialization of struct ml */
void
likelihood_CAML_serialize(value v,unsigned long* wsize_32,unsigned long* wsize_64)
{
    mll* s;
    int i,j;
    s = ML_val(v);
    caml_serialize_int_4( s->stride );
    caml_serialize_int_4( s->c_len );
    caml_serialize_int_4( s->rates );
    caml_serialize_int_4( s->invar );

    j = s->stride * s->c_len * s->rates;
    for(i = 0; i < j; ++i)
        caml_serialize_float_8( s->lv_s[i] );
    if( s->invar > 0 ){
        for(i = 0; i < s->c_len; ++i)
            caml_serialize_int_4( s->lv_invar[i] );
    }
}

/* deserialization of struct ml */  
unsigned long likelihood_CAML_deserialize( void* v )
{
    mll* s;
    double* lk_vec;
    int *pinvar,i,j;
    j = s->stride*s->c_len*s->rates;
    lk_vec = (double*) malloc( sizeof(double)*j );
    CHECK_MEM(lk_vec);
    s = (mll*) malloc( sizeof(mll) );
    CHECK_MEM(s);
    s->stride = caml_deserialize_sint_4();
    s->c_len  = caml_deserialize_sint_4();
    s->rates  = caml_deserialize_sint_4();
    s->invar  = caml_deserialize_sint_4();
    for(i = 0; i < j; ++i)
        lk_vec[i] = caml_deserialize_float_8();
    s->lv_s = lk_vec;
    if( s->invar > 0 ){
        pinvar = (int*) malloc( sizeof(int)*s->c_len );
        CHECK_MEM(pinvar);
        for(i = 0; i < s->c_len; ++i)
            pinvar[i] = caml_deserialize_sint_4();
    }
    s->lv_invar = pinvar;
    v = s; //TODO really?
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
    ret = c1->rates - c2->rates;
    if (ret != 0)
        return ret;
    if (c2->invar == c1->invar)
        return ret;
    n = c1->c_len * c1->stride * c1->rates;
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
    custom_hash_default,            //hash
    (&likelihood_CAML_serialize),   //serialize
    (&likelihood_CAML_deserialize)  //deserialize
};

/* registration with ocaml GC */
value likelihood_CAML_register (value u) {
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

value likelihood_CAML_StoBigarray( value s, value p )
{
    CAMLparam1( s );
    CAMLlocal3( res1,res2,tuple );
    mll* work;
    long dims3[3],dims1[1];
    double *newone;
    int *newtwo;

    work = ML_val( s );
    dims3[0] = work->rates;
    dims3[1] = dims1[0] = work->c_len;
    dims3[2] = work->stride;
    newone = (double*) malloc( dims3[0] * dims3[1] * dims3[2] * sizeof(double));
    memcpy(newone,work->lv_s, dims3[0] * dims3[1] * dims3[2] * sizeof(double));
    res1 = alloc_bigarray(BIGARRAY_FLOAT64 | BIGARRAY_C_LAYOUT,3,newone,dims3);
    if (work->invar > 0){
        newtwo = (int*) malloc( dims3[0] * sizeof(int));
        memcpy(newtwo,work->lv_invar, dims3[0] * sizeof(int));
        res2 = Val_some (alloc_bigarray(BIGARRAY_CAML_INT | BIGARRAY_C_LAYOUT,1,newtwo,dims1));
    } else {
        res2 = Val_none;
    }
    tuple = caml_alloc_tuple(2);
    Store_field(tuple,0,res1);
    Store_field(tuple,1,res2);
    CAMLreturn( tuple );
}

/**
 *  convert bigarray.array2 * bigarray.array1 to struct ml
 */
value likelihood_CAML_BigarraytoS( value A, value B )
{
    CAMLparam2( A,B );
    CAMLlocal1( lk ); //the abstract type
    double *lkvec,*l_stuff;
    int *i_stuff,*invar,len;
    mll* ret;

    //basic data
    l_stuff = (double*) Data_bigarray_val( A );
    ret = (mll*) malloc( sizeof(mll) );
    CHECK_MEM(ret);
    ret->rates = Bigarray_val(A)->dim[0];
    ret->c_len = Bigarray_val(A)->dim[1];
    ret->stride = Bigarray_val(A)->dim[2];

    //allocate new and copy
    lkvec = (double*) malloc( ret->rates * ret->c_len * ret->stride * sizeof(double));
    CHECK_MEM(lkvec);
    memcpy( lkvec, l_stuff, ret->rates * ret->c_len * ret->stride * sizeof(double));
    ret->lv_s = lkvec;

    //check if we have invar, then allocate new and copy
    if( Val_none == B ){
        ret->invar = 0;
    } else {
        i_stuff = (int*) Data_bigarray_val( Some_val(B) );
        len = Bigarray_val( Some_val(B))->dim[0];
        assert( len == ret->c_len );
        invar = (int*) malloc( ret->c_len * sizeof(int));
        CHECK_MEM(lkvec);
        memcpy( invar, i_stuff, ret->c_len * sizeof(int));
        ret->lv_invar = invar;
        ret->invar = 1;
    }

    //set up abstract type
    lk = caml_alloc_custom(&likelihood_custom_operations, (sizeof(mll*)),
                            ret->rates*ret->c_len*ret->stride, CAML_ALLOC_2); 
    ML_val( lk ) = ret;

    assert( ML_val(lk) == ret );
    CAMLreturn( lk );
}

/**
 * filters an ml struct with all the indexes in the passed array --ordered
 *
 * TODO with invar
 */
value likelihood_CAML_filter(value as, value ibs)
{  
    CAMLparam2( as, ibs );
    CAMLlocal1( package );
    mll *norm, *ret;
    double* area;
    int m,i,j,n,k,r,size;

    m = Wosize_val( ibs );
    norm = ML_val( as );
    n = norm->c_len;
    size = norm->c_len * norm->stride;
    area = (double*) malloc(sizeof(double)*(norm->c_len-m)*norm->stride*norm->rates);
    CHECK_MEM(area);

    for(r=0;r<norm->rates;r++){
        for(i=0,j=0;i<n;++i){     
            if( i != Int_val(Field(ibs, j)) ){
                for(k=0;k<norm->stride;k++)
                    area[r*size + i*norm->stride + k] = norm->lv_s[ r*size + i*norm->stride + k];
            } else if(j++ == m)
                break;
        }
        for(;i<m;++i)
            area[i] = norm->lv_s[i];
    }
    package = caml_alloc_custom(&likelihood_custom_operations, (sizeof(mll*)),
                                            (n-m)*norm->stride*norm->rates, CAML_ALLOC_2);
    ret = ML_val( package );
    ret->stride = norm->stride;
    ret->c_len = norm->c_len - m;
    ret->rates = norm->rates;

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
 * diagonal elements as scale factors for VL Conversly this can be done by
 * scaling VR, by scaling though each column --more of a fortran style though. 
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
    assert(a->rates == b->rates);

    prop = k = 0;
    //i holds current char, j holds current array location, k holds curr char
    for(i=0;i<a->c_len*a->rates;i++){
        maxa = maxb = 0;
        s = t = 0;
        for(j=0;j<a->stride;j++){
            if (a->lv_s[k+j] > maxa){ maxa = a-> lv_s[k+j]; }
            if (b->lv_s[k+j] > maxb){ maxb = b-> lv_s[k+j]; }
        }
        for(j=0;j<a->stride;j++,k++){
            if(a->lv_s[k] == maxa && b->lv_s[k] == maxb){ s++; t+=2; }
            else if (a->lv_s[k] == maxa || b->lv_s[k] == maxb){ t++; }
        }
        prop += (2*s)/t;
    }
    prop = prop / (a->c_len*a->rates);
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
    //wi = register_section( space, n, 0 );
    //U  = register_section( space, n*n, 0 );
    wi = (double*) malloc( sizeof(double) * n);
    U = (double*) malloc( sizeof(double) * n * n);

    //find the optimal work (call func with -1 in lwork)
    lwork = -1;
    //D holds the real values eigen values through the computation
    dgeev_(&jobv_,&jobv_,&n,A,&n,D,wi,U,&n,Ui,&n,&work_size,&lwork,&info);
    if( info == 0 ) {
        lwork = (int)work_size;
        //work = register_section( space, lwork, 1);
        work = (double*) malloc( lwork * sizeof(double));
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
 * UNCHECKED EXCEPTIONS ::
 * [n]x[n] == dim([D]) == dim([P]) == dim([U]) == dim([Ui]) == dim([TMP])
 * t > 0
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


/* calculates the likelihood for the invariant site class */
double
loglikelihood_site_invar(const mll* l,const double* pi,const int i)
{
    int j,num,set;
    double ret;
    assert( 1 == l->invar );
    j = set = 0;
    num = 1;
    ret = 0.0;
    if( 0 == l->lv_invar[i]){ return ret; }
    while( j < 32 ){ //TODO:: abstract max alphabet size
        if ( num & l->lv_invar[i] ){
            ret += pi[j]; //similar to line 846, but lv_s[..] = 1
            set++;
        }
        j++;
        num *= 2;
    }
    return ret / set;
}
/* [loglikelihood_site ml p prob]
 * calculate the likelihood for a single site */
double 
loglikelihood_site(const mll* l,const double* pi,const double* prob,const double pinvar,const int i)
{
    int j,h;
    double tmp1,tmp2;
    tmp2 = 0.0;
    for(h = 0; h < l->rates; ++h) {
        tmp1 = 0.0;
        for(j=0;j < l->stride; ++j)
            tmp1 += l->lv_s[(l->c_len * l->stride * h) + (l->stride*i) + j] * pi[j];
        tmp2 += tmp1 * prob[h];
    }
    if( 1 == l->invar ){
        assert( pinvar >= 0.0 );
        tmp1 = loglikelihood_site_invar ( l , pi, i );
        tmp2 = MAX( tmp2 * (1 - pinvar), 1e-300 ) + (pinvar * tmp1);
    } else {
        tmp2 = MAX( tmp2 , 1e-300 );
    }
    return ( log( tmp2 ) ) /* times weight */ ;
}
/* calculates the likelihood for a rate class */
double
loglikelihood_rate(const mll* l,const double* pi,const int h)
{
    int i,j;
    double tmp1, ret;
    ret = 0.0;
    for( i = 0; i < l->c_len;++i){
        tmp1 = 0.0;
         for(j=0;j < l->stride; ++j)
            tmp1 += l->lv_s[(l->c_len * l->stride * h) + (l->stride*i) + j] * pi[j];
         ret -= log ( MAX(ret, 1e-300) );
    }
    return ret;
}
double
loglikelihood_invar( const mll* l, const double *pi )
{
    int i;
    double ret,tmp;
    if (l->invar == 0){ return 0.0; }
    for( i = 0; i < l->c_len;++i){
        tmp = loglikelihood_site_invar( l , pi, i );
        ret += (tmp == 0) ? 0 : log(tmp);
    }
    return ret;
}

/** [loglikelihood ml p prob %]
 * returns loglikelihood of a character set. |p| = l->stride
 *                                           |prob| = l->rates
 *                                           if pinvar > 0.0 then l-> rates == 1
 *
 * uses Kahan Summation to reduce floating point addition error. reduces error
 * of a summation to 2e
 */
double loglikelihood( const mll* l, const double* pi, const double* prob,const double pinvar )
{
    int i; 
    double ret_s,ret_c,ret_y,ret_t; /* for kahan summation */

    ret_s = loglikelihood_site(l,pi,prob,pinvar,0);
    ret_c = 0.0;
    for(i = 1;i<l->c_len; i++)
    {
        ret_y = (loglikelihood_site(l,pi,prob,pinvar,i)) - ret_c;
        ret_t = ret_s + ret_y;
        ret_c = (ret_t - ret_s) - ret_y;
        ret_s = ret_t;
    }

    return ( -ret_s );
}

/* [likelihoood_CAML_loglikelihood s p] wrapper for loglikelihood */
value likelihood_CAML_loglikelihood(value s, value pi, value prob, value pinvar)
{
    CAMLparam3( s, pi, prob );
    CAMLlocal1( mle );
    double ll;
    ll = loglikelihood( ML_val(s),
                        (double*) Data_bigarray_val(pi),
                        (double*) Data_bigarray_val(prob),
                        Double_val( pinvar )
                      );
    mle = caml_copy_double( ll );
    CAMLreturn( mle );
}


/** [median_c Pa Pb a b c]
 * Finds the likelihood vector [c] based on vectors of its children [a] and
 * [b] with the probability matrices [Pa] and [Pb], respectively, for a single
 * rate class [rate_idx]. [Pa] and [Pb] already contain rate/time information, and
 * the [rate_idx] serves to indicate the location in the likelihood vector.
 *
 * [tmp1] is size of alphabet.
 * [tmp2] is size of alphabet.
 */
void
median_charset(const double* Pa,const double* Pb, const mll* a,const mll* b,
                mll* c, double* tmp1, double* tmp2,const int rate_idx )
{
    int i,j,len,oth;

    oth = rate_idx*(a->stride*a->c_len);
    for(i=0; i < a->c_len; ++i){
        len = oth + (i*a->stride);
        /* find half medians for each character */
        median_h( Pa, a->lv_s, len, tmp1, a->stride );
        median_h( Pb, b->lv_s, len, tmp2, a->stride );
        /* pairwise multiplication */
        for(j=0;j < a->stride;++j)
            c->lv_s[len+j] = MAX((tmp2[j]*tmp1[j]),0);
    }
}
void 
median_invar(const mll* a, const mll* b, mll* c){
    int i;
    //assert( (a->invar == b->invar) && (a->invar == 1) );
    for(i=0;i < c->c_len;++i)
        c->lv_invar[i] = (a->lv_invar[i]) & (b->lv_invar[i]);
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
    assert( cta >= 0 && ctb >= 0 );
    /* character sets */
    a = ML_val( ml_a );
    b = ML_val( ml_b );
    assert( a->stride == b->stride );
    assert( a->c_len == b->c_len );
    assert( a->rates == b->rates );
    assert( a->invar == b->invar );
    /* create new character set, c */
    ml_c = caml_alloc_custom(&likelihood_custom_operations,(sizeof(mll*)),
                                a->rates*a->stride*a->c_len,CAML_ALLOC_2);
    c = (mll*) malloc( sizeof(mll) );
    CHECK_MEM(c);
    ML_val( ml_c ) = c;
    c->stride = a->stride;
    c->c_len = a->c_len;
    c->rates = a->rates;
    c->invar = a->invar;
    c->lv_s = (double*) malloc( c->rates * c->c_len * c->stride * sizeof(double));
    CHECK_MEM(c->lv_s);
    /* register temp variables */
    PA = (double*) register_section(space, a->stride*a->stride, 1);
    PB = (double*) register_section(space, a->stride*a->stride, 1);
    tmp1 = (double*) register_section(space, a->stride*a->stride, 1);
    /** 
     *  F(x|r_i,P) = f_right(x*r_i|P) * f_left(x*r_i|P)
     *
     *      P    = probability matrix
     *      r_i  = the rate of i
     *      x    = the branch length
     */
    //printf ("rates:%d\tchars:%d\tstride:%d\n",a->rates,a->c_len,a->stride);
    for(i=0;i<num_rates;i++){
        compose_sym( PA, c_U, c_D, cta*g_rs[i], a->stride,tmp1 );
        compose_sym( PB, c_U, c_D, ctb*g_rs[i], b->stride,tmp1 );
        median_charset( PA, PB, a, b, c, tmp1, &(tmp1[a->stride]), i );
    }
    if(a->invar == 1){
        c->lv_invar = (int*) malloc ( c->c_len * sizeof(double));
        median_invar( a,b,c );
    }
    //printarray(c->lv_s, a->stride * a->c_len * a->rates);
    /* free up space, return */
    free_all( space );
    assert( c == ML_val(ml_c));
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
    assert(a->rates == b->rates);
    assert(a->invar == b->invar);
    /* create new character set */
    ml_c = caml_alloc_custom(&likelihood_custom_operations, (sizeof(mll*)), 
                                            a->rates*a->stride*a->c_len, CAML_ALLOC_2);
    c = (mll*) malloc( sizeof(mll) );
    CHECK_MEM(c);
    ML_val( ml_c ) = c;
    c->lv_s = malloc(a->c_len * a->stride * a->rates * sizeof(double)); 
    c->stride = a->stride;
    c->c_len  = a->c_len;
    c->rates  = a->rates;
    c->invar = a->invar;
    /* register temp variables */
    PA = register_section( space, b->stride * b->stride, 1 );
    PB = register_section( space, b->stride * b->stride, 1 );
    tmp1 = register_section( space, b->stride * b->stride, 1 );
    /* main loop, see symmetric for description */
    //printf("R C S:%d\t%d\t%d\n[",a->rates,a->c_len,a->stride);
    for(i=0;i<num_rates;i++){
        compose_gtr( PA, c_U, c_D, c_Ui, cta*g_rs[i], a->stride, tmp1);
        compose_gtr( PB, c_U, c_D, c_Ui, ctb*g_rs[i], b->stride, tmp1);
        median_charset( PA, PB, a,b,c, tmp1,&(tmp1[b->stride]), i);
    }
    if(a->invar == 1){
        c->lv_invar = (int*) malloc ( c->c_len * sizeof(double));
        median_invar( a,b,c );
    }
    //printf("]\n");
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
    const mll *b,const double time_a,const double time_b,const double *gam,const double *prob,
    const double *pi,const int g_n,double percent,double *tmp)
{
    int i;
    simp->time = time_a;
    for(i=0;i<g_n;i++){
        compose_sym( PA, U, D, time_a*gam[i], a->stride,tmp );
        compose_sym( PB, U, D, time_b*gam[i], b->stride,tmp );
        median_charset( PA, PB, a,b,simp->vs,tmp,&(tmp[a->stride]),i );
    }
    if(a->invar == 1) median_invar( a,b,simp->vs );
    simp->ll = loglikelihood( simp->vs,pi,prob, percent );
}
void single_gtr(ptr *simp,double *PA,double *PB,const double *U,const double *D,const double *Ui,
    const mll *a,const mll *b,const double time_a,const double time_b,const double *gam,
    const double *prob,const double *pi,const int g_n,double percent,double *tmp)
{
    int i;
    simp->time = time_a;
    for(i=0;i<g_n;i++){
        compose_gtr( PA, U, D, Ui, time_a*gam[i], a->stride, tmp );
        compose_gtr( PB, U, D, Ui, time_b*gam[i], b->stride, tmp );
        median_charset( PA, PB, a,b,simp->vs,tmp,&(tmp[a->stride]), i );
    }
    if(a->invar == 1) median_invar( a,b,simp->vs );
    simp->ll = loglikelihood( simp->vs, pi, prob, percent );
}

#define golden_exterior(a,b) a->time + ((b->time - a->time)*2.0/(sqrt(5.0)-1))

void
readjust_brents_sym(mat *space,const double* U,const double* D,const mll* a,
        const mll* b,mll* c,double* b_ta,const double b_tb,double* b_mle,
        const double* rates, const double *prob, const int g_n,const double* pi,
        const double pinvar)
{
    /* variables */
    ptr *left,*right,*middle,*best,*temp,*temp_ptr;
    double *PA,*PB,*TMP,ttime1,ttime2,deno,numr;
    int max_iter,size,bracketed;
    max_iter = MAX_ITER;
    size = a->c_len*a->stride*a->rates; //size of the entire array
    /* register temp space and matrices */
    PA = register_section(space, a->stride*a->stride , 0);
    PB = register_section(space, a->stride*a->stride , 0);
    TMP = register_section(space, a->stride*a->stride, 0);
    /* set up the 'points' */ //TODO --convert to stack space
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
    /* set lengths of vectors in new likelihood structures */
    (best->vs)->rates = (middle->vs)->rates = a->rates;
    (left->vs)->rates = (right->vs)->rates = (temp->vs)->rates = a->rates;
    (best->vs)->c_len = (middle->vs)->c_len = a->c_len;
    (left->vs)->c_len = (right->vs)->c_len = (temp->vs)->c_len = a->c_len;
    (best->vs)->stride = (middle->vs)->stride = a->stride;
    (left->vs)->stride = (right->vs)->stride = (temp->vs)->stride = a->stride;
    (best->vs)->invar = (middle->vs)->invar = a->invar;
    (left->vs)->invar = (right->vs)->invar = (temp->vs)->invar = a->invar;
    /* set initial times and constants */
    best->time = *b_ta;
    middle->time = *b_ta;
    left->time = middle->time / 10;
    right->time = middle->time * 2;
    best->ll = *b_mle;
    middle->ll = *b_mle;
    // set up best vectors
    memcpy( (best->vs)->lv_s, c->lv_s, size * sizeof(double) );
    memcpy( (middle->vs)->lv_s, c->lv_s, size * sizeof(double) );
    // fill in initial data of new points
    single_sym(left,PA,PB,U,D,a,b,left->time,b_tb,rates,prob,pi,g_n,pinvar,TMP);
    single_sym(right,PA,PB,U,D,a,b,right->time,b_tb,rates,prob,pi,g_n,pinvar,TMP);
    /* bracket minimum */
    bracketed = 1; // true, since we are going to expect the best
    while( (max_iter-- != 0) && !(left->ll > middle->ll && middle->ll < right->ll)){
        if (left->ll <= middle->ll && middle->ll <= right->ll){
            // increasing,  move all points down
            temp_ptr = right;
            right = middle;
            middle = left;
            left = temp_ptr;
            ttime1 = golden_exterior( middle, right );
            if (ttime1 < EPSILON) { ttime1 = EPSILON; break; }
            single_sym(left,PA,PB,U,D,a,b,ttime1,b_tb,rates,prob,pi,g_n,pinvar,TMP);
        } else if ( left->ll >= middle->ll && middle->ll >= right->ll) {
            // decreasing,  move all points up (middle,right,new)
            temp_ptr = left;
            left = middle;
            middle = right;
            right = temp_ptr;
            ttime1 = golden_exterior( left, middle );
            single_sym(right,PA,PB,U,D,a,b,ttime1,b_tb,rates,prob,pi,g_n,pinvar,TMP);
        } else {
            //printf("Curvature?:%f/%f -- %f/%f -- %f/%f\n",left->time,
            //        left->ll,middle->time,middle->ll,right->time,right->ll);
            bracketed = 0;
            if (left->ll < right->ll){ //new left middle
                temp_ptr = right;
                right = middle;
                middle = left;
                left = temp_ptr;
                ttime1 = golden_exterior (middle, right);
                if (ttime1 < EPSILON) { ttime1 = EPSILON; break; }
                single_sym(left,PA,PB,U,D,a,b,ttime1,b_tb,rates,prob,pi,g_n,pinvar,TMP);
            } else { //middle right new
                temp_ptr = left;
                left = middle;
                middle = right;
                right = temp_ptr;
                ttime1 = golden_exterior( left, middle );
                single_sym(right,PA,PB,U,D,a,b,ttime1,b_tb,rates,prob,pi,g_n,pinvar,TMP);
            }
            break;
        }
    }
    //printf("Bracketed in %d iterations\n",MAX_ITER - max_iter);
    //printf("Some New Likelihood: %f(%f)\t%f(%f)\t%f(%f)\n",
    //        left->ll,left->time, middle->ll,middle->time, right->ll,right->time);
    /* parabolic interpolation */
    while( bracketed && ((max_iter--) > 0 ) &&
                      (fabs (middle->time - left->time) > EPSILON) &&
                      (fabs (right->time - middle->time) > EPSILON) ){
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
        single_sym(temp,PA,PB,U,D,a,b,ttime1,b_tb,rates,prob,pi,g_n,pinvar,TMP);
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
    *b_ta = best->time;
    *b_mle= best->ll;

    free(best->vs); free(middle->vs); free(left->vs); free(right->vs); free(temp->vs);
    free(best); free(middle); free(left); free(right); free(temp);
}

void
readjust_brents_gtr(mat * space,const double* U,const double* D,const double* Ui,
        const mll* a, const mll* b, mll* c,double* b_ta,const double b_tb,double* b_mle,
        const double *rates,const double *prob,const int g_n,const double* pi,const double pinvar)
{
}

//-----------------------------------------------------------------------------
/** [..._readjust_gtr... tmp U D Ui A B C ta tb gammas probs pis ll] **/
value
likelihood_CAML_readjust_gtr_wrapped
    (value tmp,value U,value D,value Ui,value A,value B,value C,value ta,value tb,
        value pinv, value g,value p,value pis,value ll)
{   
    /* ocaml macros */
    CAMLparam5(U,D,Ui,A,B);
    CAMLxparam5(C,ta,tb,g,p);
    CAMLxparam3(ll,pis,tmp);
    CAMLlocal1( res );
    /* declare/unbox variables to update/return */
    double cta,likelihood;
    cta = Double_val( ta );
    likelihood = Double_val( ll );
    /* check rates, probabilities and priors */
    assert( (ML_val(A)->stride) == Bigarray_val(pis)->dim[0] );
    assert( Bigarray_val(g)->dim[0] == Bigarray_val(p)->dim[0] );
    /* readjust */
    //printf("%f\t%f\n",cta,,likelihood);
    readjust_brents_gtr(
            FM_val(tmp),
            Data_bigarray_val( U ),
            Data_bigarray_val( D ),
            Data_bigarray_val( Ui ),
            ML_val( A ), ML_val( B ), ML_val( C ),
            &cta, Double_val( tb ), &likelihood,
            Data_bigarray_val(g),
            Data_bigarray_val(p),
            Bigarray_val( g )->dim[0],
            Data_bigarray_val(pis),
            Double_val( pinv )
    );
    //printf("%f\t\t%f\n",cta,likelihood);
    /* construct tuple / freespace */
    res = caml_alloc_tuple( 2 );
    Store_field(res, 0, caml_copy_double(cta));
    Store_field(res, 1, caml_copy_double(likelihood));
    free_all( FM_val(tmp) );
    CAMLreturn(res);
}
value
likelihood_CAML_readjust_sym_wrapped
    (value tmp,value U,value D,value A,value B,value C,value ta,value tb,value pinv,value g,value p,value pis,value ll)
{
    /* ocaml macros for GC */
    CAMLparam5(U,D,ll,A,B);
    CAMLxparam5(C,ta,tb,g,p);
    CAMLxparam2( pis,tmp );
    CAMLlocal1( res );
    /* declare/unbox variables to return */
    double cta,likelihood;
    cta = Double_val( ta );
    likelihood = Double_val( ll );
    /* check probabilities, rates, priors */
    assert( ML_val(A)->stride == Bigarray_val( pis )->dim[0]);
    assert( Bigarray_val( g )->dim[0] == Bigarray_val( p )->dim[0] );
    /* readjust loop */
    //printf("s: %f\t%f\t%f\n", cta,ctb,likelihood);
    readjust_brents_sym(
            FM_val(tmp),
            Data_bigarray_val( U ), Data_bigarray_val( D ),
            ML_val(A), ML_val(B), ML_val(C),
            &cta, Double_val( tb ), &likelihood,
            Data_bigarray_val( g ),
            Data_bigarray_val( p ),
            Bigarray_val( g )->dim[0],
            Data_bigarray_val( pis ),
            Double_val( pinv )
    );
    //printf("e: %f\t%f\t%f\n", cta,ctb,likelihood);
    /* construct tuple / free space */
    res = caml_alloc_tuple( 2 );
    Store_field(res, 0, caml_copy_double (cta));
    Store_field(res, 1, caml_copy_double (likelihood));
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
