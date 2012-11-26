/* POY 5.0 Alpha. A phylogenetic analysis program using Dynamic Homologies.   *\
(* Copyright (C) 2011 Andrés Varón, Lin Hong, Nicholas Lucaroni, Ward Wheeler,*)
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
#include <stdlib.h>         //malloc, calloc, srand, RAND_MAX, posix_memalign
#include <string.h>         //memcpy, memset
#include <assert.h>         
#include <time.h>

#undef DEBUG                //will print_loading/unloading SSE information

//** these do not work properly yet. 
#undef __SSE4__
#undef __SSE3__

#if defined( __SSE4__ )
    #define LK_DATA_ALIGNMENT
    #include <smmintrin.h>  //header for core2 duo SSE4.1 intrinsics
#elif defined( __SSE3__ )
    #define LK_DATA_ALIGNMENT
    #include <pmmintrin.h>  //SSE3 instrinsics (Pentium 4)
#endif

#include "config.h"         //defines, if likelihood, USE_LIKELIHOOD
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

#ifdef USE_LIKELIHOOD

#include "floatmatrix.h"
#include "likelihood.h"     //includes floatmatrix

//CONSTANTS
#define EPSILON      1e-10   // error for numerical calculations
#define MAX_ITER     200     // number of iterations for brents method
#define BL_MIN       1e-13   // minimum branch length 
#define BL_MAX       100     // maximum branch length

#define NEGINF       log(0)
#define INF          1/0
#define TRUE         1
#define FALSE        0

/* Cost Fn Mode Codes */
#define MPLCOST      1
#define MALCOST      0

#define KAHANSUMMATION       /* reduce error in sum over likelihood vector? */
 
#define MAX(a,b) (a >= b)?a:b
#define MIN(a,b) (a <= b)?a:b
#define SHIFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d)
#define SIGN(a,b) (b>0)?fabs(a):-fabs(a)
#define EQUALS(a,b) fabs(a-b)<EPSILON

/* ~CHECK_ZEROES    -- verify array is all zeroes with EPSILON error */
/*                     used to ensure no imaginary roots when diagonalizing */
#define CHECK_ZEROES(a,t,n); for(t=0;t<n;t++){if(a[t] > EPSILON || a[t] < -EPSILON){failwith("Imaginary eigenvalues");}}

/* for unboxing abstract types - get pointer or actual value */
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

/** memory allocation and verification; posix memalign returns 0 on success */
#define CHECK_MEM(a); if(a==NULL || a==0){printf("LK:failed on %d, ", __LINE__);failwith("I cannot allocate more memory.");}

/** MAC alloc functions automatically align data for SSE and MMX type functions **/
#if defined( __APPLE__ )
	#define lk_malloc(x,y);x=(double*)malloc(y);if(0==x){printf("LK:failed on %d",__LINE__);failwith("I cannot allocate more memory.");}
#elif defined( __SSE4__ )
	#define lk_malloc(x,y);if(0!=posix_memalign((void*)&x,16,y)){printf("LK:failed on %d",__LINE__);failwith("I cannot allocate more memory.");}
#elif defined( __SSE3__ )
	#define lk_malloc(x,y);if(0!=posix_memalign((void*)&x,16,y)){printf("LK:failed on %d",__LINE__);failwith("I cannot allocate more memory.");}
#else
	#define lk_malloc(x,y);x=(double*)malloc(y);if(0==x){printf("LK:failed on %d",__LINE__);failwith("I cannot allocate more memory.");}
#endif
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
 * D is a diagonal matrix (usually of eigenvalues)
 * U are eigenvectors (column major, Ut is rowmajor)
 * pi are the priors as a vector
 * Dpi are the priors as a diagonal matrix
 *
 * variable:
 *      n/m is the number of rows/columns or size of the alphabet
 *      t is the branch length
 */

//------------------------------------------------------------------------------
/* prints a matrix (block format) */
void printmatrixf( const double* Z, const int s, const int n, const int m)
{
    int i,j;
    for (i=0; i<m; ++i) {
        putchar('\t');
        for (j=0; j<n; ++j)
            printf("[%11.10f] ", Z[i*s+j]);
        putchar('\n');
    }
}
void printmatrixi( const int* Z, const int s, const int n, const int m)
{
    int i,j;
    for (i=0; i<m; ++i) {
        putchar('\t');
        for (j=0; j<n; ++j)
            printf("[%d] ", Z[i*s+j]);
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
void printarrayi( const int* a, const int n )
{
    int i;
    for(i=0;i<n;++i)
        printf("[%d] ", a[i]);
    putchar('\n');
}

/* print out a vector */
void debug( const mll* a )
{
    printf("DEBUG:\n\tR: %d\tS: %d\tA: %d\tL: %d\n",
                a->rates,a->stride,a->alph,a->c_len);
    printmatrixf(a->lv_s,a->stride,a->alph,a->c_len*a->rates);
    if(!(0 == a->invar)){
        printarrayi(a->lv_invar,a->c_len);
    }
    putchar('\n');
}

/* print out a caml value */
value likelihood_CAML_debug( value s )
{
    CAMLparam1( s );
    debug( ML_val( s ) );
    CAMLreturn( Val_unit );
}

/** creates a random substituation rate matrix
 * each row sums to 0 on diagonal elements */
void rand_sub_mat_gtr( double* A, const int n, const int s)
{
    int i,j;
    double diag;
    srand( time(NULL) );
    for(i=0;i<n;++i){
        diag = 0;
        for(j=0;j<n;++j){
            if( i == j ){ continue; }
            A[i*n+j] = (double)rand() / (double)RAND_MAX;
            diag = diag + A[i*n+j];
        }
        A[i*n+i] = -diag;
    }
}

/** transpose of matrix [A] of [n]x[n], must be square */
void transpose( double *A, const int n, const int s )
{
    double tmp;
    int i,j;
    for(i=0;i<n;++i){
        for(j=i+1;j<n;++j){
            tmp = A[i*s+j];
            A[i*s+j] = A[j*s+i];
            A[j*s+i] = tmp;
        }
    }
}

/** creates a random symmetric sub-rate matrix
 * each row sums to 0 on diagonal elements */
void rand_sub_mat_sym( double* A, const int n, const int s)
{
    int i,j;
    double temp,diag;
    srand( time(NULL) );
    for(i=0;i<n;++i){
        diag = 0;
        for(j=0;j<n;++j){
            if( i <= j ) { continue; }
            temp = ((double)rand()/(double)RAND_MAX);
            A[i*s+j] = temp;
            A[j*s+i] = temp;
            diag = diag + temp + temp;
        }
        A[i*s+i] = -diag;
    }
}


//------------------------------------------------------------------------------
// O'CAML interface functions for Garbage Collection + callback for a GC

int CAML_ALLOC = 10000; /* # of nodes to alloc before GC kicks in */
/* n - #taxa, r - #rates, c - #characters, a - alphabet, i - invar?, 
 * x - holds the number of nodes to keep before GC occurs */
value likelihood_GC_custom_max( value n )
{
    CAMLparam1( n );
    CAML_ALLOC = Int_val( n );
    CAMLreturn( Val_unit );
}

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
    caml_serialize_int_4( s->stride);
    caml_serialize_int_4( s->alph  );
    caml_serialize_int_4( s->c_len );
    caml_serialize_int_4( s->rates );
    caml_serialize_int_4( s->invar );
    /* alph is not used since we serialize empty stride buffer */
    j = s->stride * s->c_len * s->rates;
    for(i = 0; i < j; ++i)
        caml_serialize_float_8( s->lv_s[i] );
    if( s->invar > 0 ){
        for(i = 0; i < s->c_len; ++i)
            caml_serialize_int_4( s->lv_invar[i] );
    }
    return;
}

/* deserialization of struct ml */  
unsigned long likelihood_CAML_deserialize( void* v )
{
    mll* s;
    double* lk_vec;
    int *pinvar,i,j;
    s = (mll *) v;
    s->stride = caml_deserialize_sint_4();
    s->alph   = caml_deserialize_sint_4();
    s->c_len  = caml_deserialize_sint_4();
    s->rates  = caml_deserialize_sint_4();
    s->invar  = caml_deserialize_sint_4();
    j = s->stride * s->c_len * s->rates;
    lk_malloc( lk_vec, sizeof(double)*j);
    for(i = 0; i < j; ++i)
        lk_vec[i] = caml_deserialize_float_8();
    s->lv_s = lk_vec;
    if( s->invar > 0 ){
        pinvar = (int*) malloc( sizeof(int)*s->c_len );
        CHECK_MEM(pinvar);
        for(i = 0; i < s->c_len; ++i)
            pinvar[i] = caml_deserialize_sint_4();
        s->lv_invar = pinvar;
    }
    return (sizeof( mll ));
}

/* compares two character sets */
int compare_chars( const mll* c1, const mll* c2)
{
    int ret,i,n;
    ret = c1->alph - c2->alph;
    if (ret != 0)
        return ret;
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

int likelihood_CAML_compare( value c1, value c2 ){
    return compare_chars( ML_val(c1), ML_val(c2) );
}

/* custom garbage collection */
static struct custom_operations likelihood_custom_operations  = {
    "http://www.amnh.org/poy/likelihood/likelihood.0.1", //identifier
    (&likelihood_CAML_free),        //finalize
    (&likelihood_CAML_compare),     //compare
    custom_hash_default,            //hash
    (&likelihood_CAML_serialize),   //serialize
    (&likelihood_CAML_deserialize)  //deserialize
};

/* copy a likelihood vector */
value likelihood_CAML_copy( value x )
{
    CAMLparam1( x );
    CAMLlocal1( r );
    mll *cur, *ret;

    cur = ML_val( x );
    ret = (mll*) malloc( sizeof(mll) );
    CHECK_MEM(ret);
    ret->rates  = cur->rates;
    ret->alph   = cur->alph;
    ret->c_len  = cur->c_len;
    ret->stride = cur->stride;
    ret->invar  = cur->invar;
    lk_malloc( ret->lv_s, cur->rates * cur->c_len * cur->stride * sizeof(double) );
    memcpy( ret->lv_s, cur->lv_s, ret->c_len * cur->rates * cur->stride * sizeof(double));
    if (1 == ret->invar){
        ret->lv_invar = (int*) malloc( ret->c_len * sizeof(int));
        CHECK_MEM(ret->lv_s);
        memcpy( ret->lv_invar, cur->lv_invar, ret->c_len * sizeof(int) );
    }
    r = caml_alloc_custom(&likelihood_custom_operations, sizeof(mll*),1,CAML_ALLOC);
    ML_val( r ) = ret;
    CAMLreturn( r );
}

/* registration with ocaml GC */
value likelihood_CAML_register (value u)
{
    CAMLparam1(u);
    register_custom_operations (&likelihood_custom_operations);
    CAMLreturn (Val_unit);
}

/* convert S to ocaml types (bigarray) for printing */
value likelihood_CAML_StoBigarray( value s, value p )
{
    CAMLparam1( s );
    CAMLlocal3( res1,res2,tuple );
    mll* work;
    long dims3[3],dims1[1];
    double *newone;
    int *newtwo, r, c, a, size;

    work = ML_val( s );
    dims3[0] = work->rates;
    dims3[1] = dims1[0] = work->c_len;
    dims3[2] = work->alph;
    size = dims3[0] * dims3[1] * dims3[2];
    newone = (double*) malloc( size * sizeof(double));
    if (work->alph == work->stride){
        memcpy(newone,work->lv_s, size * sizeof(double));
    } else {
        for( r=0; r < work->rates; ++r ){
            for( c=0; c < work->c_len; ++c ){
                for( a=0; a < work->alph; ++a ){
                    newone[(r*(work->alph*work->c_len))+(c*work->alph)+a]
                        = work->lv_s[(r*(work->c_len*work->stride))+(c*work->stride)+a];
                }
            }
        }
    }
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

/*  convert bigarray.array2 * bigarray.array1 to struct ml */
value likelihood_CAML_BigarraytoS( value A, value B, value mpl )
{
    CAMLparam2( A,B );
    CAMLlocal1( lk ); //the abstract type
    double *lkvec, *l_stuff;
    int *i_stuff, *invar, len, r;
    mll* ret;

    #ifdef LK_DATA_ALIGNMENT

        int c,a;
        //basic data
        l_stuff = (double*) Data_bigarray_val( A );
        ret = (mll*) malloc( sizeof(mll) );
        CHECK_MEM(ret);
        ret->rates = Bigarray_val(A)->dim[0];
        ret->c_len = Bigarray_val(A)->dim[1];
        ret->alph  = Bigarray_val(A)->dim[2];
        ret->stride = ret->alph + (4 - (ret->alph % 4)); //4 doubles = 16 bytes
        //allocate new and copy
        lk_malloc(lkvec, ret->rates * ret->c_len * ret->stride * sizeof(double));
        ret->lv_s = lkvec;
        if( MPLCOST == Int_val(mpl) ){
            for( r=0; r<ret->rates; ++r ){
                for( c=0; c<ret->c_len; ++c ){
                    for( a=0; a<ret->alph; ++a ){
                        ret->lv_s[(r*ret->c_len*ret->stride)+(c*ret->stride)+a]
                            = (l_stuff[(r*ret->c_len*ret->alph)+(c*ret->alph)+a] >= 1.0)?0.0:NEGINF;
                    }
                }
            }
        } else {
            for( r=0; r<ret->rates; ++r ){
                for( c=0; c<ret->c_len; ++c ){
                    for( a=0; a<ret->alph; ++a ){
                        ret->lv_s[(r*ret->c_len*ret->stride)+(c*ret->stride)+a]
                            = l_stuff[(r*ret->c_len*ret->alph)+(c*ret->alph)+a];
                    }
                }
            }
        }
    #else

        //basic data
        l_stuff = (double*) Data_bigarray_val( A );
        ret = (mll*) malloc( sizeof(mll) );
        CHECK_MEM(ret);
        ret->rates = Bigarray_val(A)->dim[0];
        ret->c_len = Bigarray_val(A)->dim[1];
        ret->alph  = Bigarray_val(A)->dim[2];

        ret->stride = ret->alph;
        //allocate new and copy
        lk_malloc(lkvec, ret->rates * ret->c_len * ret->stride * sizeof(double));
        ret->lv_s = lkvec;
        if( MPLCOST == Int_val(mpl) ){
            for( r = 0; r < (ret->rates*ret->c_len*ret->stride); ++r){
                ret->lv_s[r] = (l_stuff[r] >= 1.0)?0.0:NEGINF;
            }
        } else {
            memcpy( lkvec, l_stuff, ret->rates * ret->c_len * ret->stride * sizeof(double));
        }
    #endif

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
    lk = caml_alloc_custom(&likelihood_custom_operations, (sizeof(mll*)),1,CAML_ALLOC); 
    ML_val( lk ) = ret;

    assert( ML_val(lk) == ret );
    CAMLreturn( lk );
}

/** filters an ml struct with all the indexes in the passed array --ordered */
value likelihood_CAML_filter(value as, value ibs)
{  
    CAMLparam2( as, ibs );
    CAMLlocal1( package );
    mll *norm, *ret;
    double* area;
    int m,i,j,k,r,n,size;
    int* invr = 0; /* will not be used it norm->invar = 0 */

    m = Wosize_val( ibs );
    norm = ML_val( as );
    size = norm->c_len * norm->stride;
    
    area = (double*) malloc(sizeof(double)*(norm->c_len-m)*norm->stride*norm->rates);
    CHECK_MEM(area);
    if(!(0 == norm->invar)){
        invr = (int*) malloc(sizeof(int)*(norm->c_len-m));
        CHECK_MEM(invr);
    }

    //for each character, do ... 
    for(i=0,j=0,n=0;i<norm->c_len;++i){
        //printf("Testing index: %d...", i);
        if( i != Int_val(Field(ibs, j)) ){
            /* copy */
            //printf("Copying: ");
            for(r=0;r<norm->rates;r++){
                for(k=0;k<norm->alph;k++){
                    //printf("%5f\t", norm->lv_s[ r*size + i*norm->stride + k] );
                    area[r*size + n*norm->stride + k] = norm->lv_s[ r*size + i*norm->stride + k];
                }
            }
            //printf("\n");
            if(!(0 == norm->invar))
                invr[n] = norm->lv_invar[i];
            ++n;
        } else {
            //printf("filtering\n");
            if( j == m ){ break; } else { ++j; }
        }
    }
    //fill in left-over; naivly is fine, function isn't used much
    for(;i<norm->c_len;i++,n++){
        //printf("Testing index: %d...Copying", i);
        for(r=0;r<norm->rates;r++){
            for(k=0;k<norm->alph;k++){
                //printf("%5f\t", norm->lv_s[ r*size + i*norm->stride + k] );
                area[r*size + n*norm->stride + k] = norm->lv_s[ r*size + i*norm->stride + k];
            }
        }
        //printf("\n");
        if(!(0 == norm->invar))
            invr[n] = norm->lv_invar[i];
    }
    assert( n == norm->c_len-m );

    package = caml_alloc_custom(&likelihood_custom_operations, (sizeof(mll*)),1,CAML_ALLOC);
    ret = (mll*) malloc( sizeof(mll) );
    CHECK_MEM(ret);
    ML_val( package ) = ret;

    ret->lv_s  = area;
    ret->stride= norm->stride;
    ret->alph  = norm->alph;
    ret->c_len = norm->c_len - m;
    ret->rates = norm->rates;
    ret->invar = norm->invar;
    if( !(0 == norm->invar) )
        ret->lv_invar = invr;

    assert( ret == ML_val(package));
    CAMLreturn( package );
}

/**  [mk_diag diag mat n m] ~ makes a diagonal square matrix from the first [n]
 * elements of the matrix. Put first [n] elements of [M] along diagonal of [M].
 * [s] becomes the stride factor of the matrix. */
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
void scale_VL( const double *VR, double *VL, const int n, const int s)
{   NOT USED CURRENTLY --calculating the inverse is better.
    double *alphas,a,b;
    int i,j; 
    char ntran, _tran;

    ntran = 'N'; _tran = 'T';
    a = 1; b = 0;
    alphas = (double*) malloc(n*n*sizeof(double));
    CHECK_MEM( alphas );

    cblas_dgemm( CblasRowMajor, CblasTrans, CblasNoTrans,
            n, n, n, a, VL, s, VR, s, b, alphas, n );

    for(i = 0;i < n; i++){
        for(j=0;j<n;j++)
            VL[i*s+j] = VL[i*s+j] / alphas[i*n+i];
    }
    free(alphas);
} */

/** [mk_inverse VL D VR] makes the inverse of VL from the diagonalization
 *
 * Since lapack doesn't guarentee that VL = inv(VR) when we have duplicate
 * eigenvalues (this happens less often in GTR, but symmetric matrices have
 * an easy solution, Ut = Ui). Otherwise we just scale Ui's rows (scale_VL).
 *
 * output is in VL.  */
int
mk_inverse(mat *space,double *VL, const double *VR, int n)
{
    //LU factorization and inverse via DGETRF and DGETRI
    int i,lwork,*pivot;
    double work_size, *work;

    lwork = -1;
    pivot = (int*) malloc(n * sizeof(int));
    memcpy(VL, VR, n*n*sizeof(double) );
    dgetrf_(&n, &n, VL, &n, pivot, &i);
    if( 0 == i ){
        dgetri_(&n, VL, &n, pivot, &work_size, &lwork, &i); //optimal work
        if( 0 == i ){
            lwork = MIN( (int)work_size, free_space(space) );
            work = register_section( space, lwork, 1 );
            dgetri_(&n, VL, &n, pivot, work, &lwork, &i);
            if ( i < 0 ){
                failwith ( "dgetri_ argument failed. 1" );
            } else if (i > 0){
                failwith ( "dgetri_ matrix is singular and inverse cannot be computed. 1" );
            }
        } else if (i < 0) { 
            failwith ( "dgetri_ argument failed. 2" );
        } else {
            failwith ( "dgetri_ matrix is singular and inverse cannot be computed. 2" );
        }
    } else {
        failwith ( "dgetrf_ unknown error." );
    }
    free (pivot);
    return i;
}

/**  [apply_exp diag n m t] Multiplies the diagonal of [diag], an [n]x[m]
 * matrix by [t] and applies exp() */
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

/** [log_matrix m] Take the log of each value of the matrix */
#ifdef _WIN32
__inline void 
#else
inline void 
#endif
log_matrix(double* PA, const int n, const int m){
    int i,s;
    s = m * n;
    for( i=0;i<s;++i )
        PA[i] = log( PA[i] );
}

/** [create_identity P n] Fill a square matrix with I; width && height = n */
#ifdef _WIN32
__inline void 
#else
inline void 
#endif
create_identity(double* P, const int n){
    int i,j;
    for( i=0;i<n;++i ){
        for( j=0;j<n;++j )
            P[i*n+j] = ( i == j )? 1.0 : 0.0;
    }
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
        for(j=0;j<a->alph;j++){
            if (a->lv_s[k+j] > maxa){ maxa = a->lv_s[k+j]; }
            if (b->lv_s[k+j] > maxb){ maxb = b->lv_s[k+j]; }
        }
        for(j=0;j<a->stride;j++,k++){
            if(a->lv_s[k] == maxa && b->lv_s[k] == maxb){ s++; t+=2; }
            else if (a->lv_s[k] >= maxa || b->lv_s[k] >= maxb){ t++; }
        }
        prop += (t > 0) ? ((2*s)/t) : (0);
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

/** Return the minimum branch length used in the optimization routines. */
value likelihood_CAML_minimum_bl( value unit )
{
    CAMLparam1(unit);
    CAMLlocal1(minimum);
    minimum = caml_copy_double( BL_MIN );
    CAMLreturn( minimum );
}

/** set the smoothness for SMPL; this may or may not be used, and could possibly
 * be boot-strapped so that MAX(a,b) + Epsilon < SMAX(a,b); thus, smoothness has
 * occured, and is great enough to be effective **/
double SMPL_K = 100.0;
value likelihood_CAML_set_smoothness( value a ){
    CAMLparam1( a );
    SMPL_K = Double_val( a );
    CAMLreturn( Val_unit );
}

/** Set the TOLERANCE used in the optimization functions for the convergence of
 * the brents method functions used in optimizing the branch lengths. */
double BRENT_TOL = 1e-6;
value likelihood_CAML_set_tol( value a ){
    CAMLparam1( a );
    BRENT_TOL = Double_val( a );
    CAMLreturn( Val_unit );
}

/** Get the TOLERANCE used in the optimization functions for the convergence of
 * the brents method functions used in optimizing the branch lengths. */
value likelihood_CAML_get_tol( value unit ){
    CAMLparam1( unit );
    CAMLlocal1( brent_tol );
    brent_tol = caml_copy_double( BRENT_TOL );
    CAMLreturn( brent_tol );
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
int diagonalize_sym( mat *space,  double* A, double* D, int n)
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
        expand_matrix( space, lwork ); 
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
value likelihood_CAML_diagonalize_sym( value tmp, value Q, value D)
{
    CAMLparam3( tmp, Q, D );
    diagonalize_sym(
            FM_val( tmp ),
            (double*) Data_bigarray_val( Q ),
            (double*) Data_bigarray_val( D ),
            Bigarray_val( Q )->dim[0]);
    CAMLreturn( Val_unit );
}

int diagonalize_gtr(mat *space, double* A, double* D, double* Ui, int n)
{
    char jobv_ = 'V'; //we went to find left and right eigenvectors
    double *wi,*U,*work,work_size;
    int lwork,info;

    //find the optimal work (call func with -1 in lwork)
    lwork = -1;
    //D holds the real values eigen values through the computation
    expand_matrix( space, n + (n * n) );
    wi = register_section( space, n, 1 );
    U  = register_section( space, n*n, 1 );
    dgeev_(&jobv_,&jobv_,&n,A,&n,D,wi,U,&n,Ui,&n,&work_size,&lwork,&info);
    if( info == 0 ) {
        lwork = (int)work_size;
        expand_matrix( space, (2*n) + n*n + (2*lwork) ); //reallocate wi and U
        wi = register_section( space, n, 1 );
        U  = register_section( space, n*n, 1 );
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
        CHECK_ZEROES(wi,lwork,n);
        if( 0 == info ) {
            mk_diag(D,n);
            mk_inverse(space, U, Ui, n);
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
value
likelihood_CAML_diagonalize_gtr( value tmp, value Q, value D, value Qi)
{
    CAMLparam4( tmp, Q, D, Qi );
    diagonalize_gtr(
            FM_val( tmp ),
            (double*) Data_bigarray_val( Q ),
            (double*) Data_bigarray_val( D ),
            (double*) Data_bigarray_val( Qi),
            Bigarray_val( Q )->dim[0]);
    CAMLreturn( Val_unit );
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
    if(t == -1.0){ /** Instantaneous rate matrix **/
        dgemm_(&ntran,&ntran, &n, &n, &n, &alpha, U, &n, P, &n, &beta, TMP, &n );
        dgemm_(&ntran,&_tran,&n,&n,&n,&alpha,TMP,&n,U,&n,&beta,P,&n);
/*    } else if(t <= EPSILON){*/
/*        create_identity( P , n );*/
    } else {
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
}
value
likelihood_CAML_compose_sym(value tmp,value U, value D, value t)
{   /* only called for testing purposes */
    CAMLparam4( tmp,U,D,t );
    CAMLlocal1( res );
    double *c_P, c_t,*c_D,*c_U,*c_T;
    mat *space;
    int n;
    long dims[2];

    n = Bigarray_val( U )->dim[0];
    space = FM_val (tmp);
    expand_matrix (space, n*n);//ensure space is available

    c_t = Double_val( t );
    c_U = (double *) Data_bigarray_val( U );
    c_D = (double *) Data_bigarray_val( D );
    c_P = (double*) malloc( n*n*sizeof(double));
    CHECK_MEM(c_P);
    c_T = register_section( space, n*n, 0 );

    compose_sym(c_P,c_U,c_D,c_t,n,c_T);
    dims[0] = n; dims[1] = n;
    res = alloc_bigarray(BIGARRAY_FLOAT64 | BIGARRAY_C_LAYOUT, 2, c_P, dims);
    CAMLreturn ( res );
}

void
compose_gtr(double* P, const double* U, const double* D, const double* Ui, 
        const double t, const int n,double *tmp)
{
    double alpha,beta; char ntran;
    alpha = 1; beta = 0; ntran = 'N';
    memcpy(P,D,n*n*sizeof(double));
    if(t == -1.0){
        dgemm_(&ntran,&ntran,&n,&n,&n,&alpha,Ui,&n,P,&n,&beta,tmp,&n);
        dgemm_(&ntran,&ntran,&n,&n,&n,&alpha,tmp,&n,U,&n,&beta,P,&n);
    } else if(t >= EPSILON || t == -1.0){
        apply_exp(P,n,n,t);
        //S is U*exp(D)...
        dgemm_(&ntran,&ntran,&n,&n,&n,&alpha,Ui,&n,P,&n,&beta,tmp,&n);
        //P becomes U*expD*Ui... done --note: VL = inv(VR)
        dgemm_(&ntran,&ntran,&n,&n,&n,&alpha,tmp,&n,U,&n,&beta,P,&n);
    } else { //identity matrix
        create_identity( P, n );
    }
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

    n = Bigarray_val( U )->dim[0];
    space = FM_val( tmp );
    expand_matrix( space, n*n );

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

    CAMLreturn( res );
}



/**  [median_h P l c nl a]  
 * Finds half the likelihood vector of child [l] and probability matrix [P] into [nl]
 *
 *  [a] - length of each vector (in practice, the alphabet size) [P] is [a]x[a]
 *  [c] - the start of the character to do work on */
#ifdef _WIN32
__inline void 
#else
inline void 
#endif
median_h( const double* P, const double* l, const int c_start, double* tmp, const int a)
{
    int i,j; double elm;
    //for each row of P
    for(i=0; i<a; ++i){
        elm = 0;
        for(j=0;j<a;++j)
            elm += P[(i*a)+j] * l[c_start+j];
        tmp[i] = elm;
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
    while( j < 32 ){ //TODO:: max alphabet size?
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
#ifdef _WIN32
__inline double
#else
inline double
#endif
loglikelihood_site(const mll* l, const double weight, const double* pi,
                    const double* prob, const double pinvar,const int i )
{
    int j,h,r_start,c_start;
    double tmp1,tmp2;
    tmp2 = 0.0;
    r_start = 0;
    for(h = 0; h < l->rates; ++h,r_start += l->c_len ) {
        tmp1 = 0.0;
        c_start = (r_start + i) * l->stride;
        for(j=0;j < l->alph; ++j,++c_start)
            tmp1 += l->lv_s[c_start] * pi[j];
        tmp2 += tmp1 * prob[h];
    }
    if( 1 == l->invar ){
        assert( pinvar >= 0.0 );
        tmp1 = loglikelihood_site_invar ( l , pi, i );
        tmp2 = MAX( tmp2 * (1 - pinvar), 1e-300 ) + (pinvar * tmp1);
    } else {
        tmp2 = MAX( tmp2 , 1e-300 );
    }
    return ( log( tmp2 ) * weight );
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
        for(j=0;j < l->alph; ++j)
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
    ret = 0.0;
    if (l->invar == 0){ return ret; }
    for( i = 0; i < l->c_len;++i){
        tmp = loglikelihood_site_invar( l , pi, i );
        ret += (tmp == 0) ? 0 : log(tmp);
    }
    return ret;
}

double
logMPL_site( const mll* l, const double weight, const double* pi,
                    const double* prob, const int i )
{
    int r, j, c;
    double max_v;
    max_v = c = NEGINF; /* placeholder for C; to avoid compiler warnings */
    for(r=0; r < l->rates;++r){
        c = (r * (l->stride * l->c_len)) + (l->stride * i);
        for(j=0; j < l->alph; ++j){
            max_v = MAX (l->lv_s[c+j] + log(pi[j]), max_v);
        }
    }
    max_v= max_v * weight;
    return max_v;
}

/** Smooth/soft Max function; does not contain a sharp discontinuity;

                        log ( \summation_{i=0}^{n} \exp^{kx_i} )
    smax(x_0...x_n|k) = ----------------------------------------
                                            k

        where k is some smoothness factor. If we set k to the actual maximum,
    we can avoid some terms and reduce error, and better reflect the actual 
    max function.

    See, http://www.johndcook.com/blog/2010/01/13/soft-maximum/
    and, http://www.johndcook.com/blog/2010/01/20/how-to-compute-the-soft-maximum/ */
double
logSMPL_site( const mll* l, const double weight, const double* pi,
                    const double* prob, const int i )
{
    int r, j, c;
    double maximum, s_max;
    //maximum = logMPL_site( l, 1.0, pi, prob, i );
    maximum = 0.0;
    s_max = 0;
    for(r=0; r < l->rates;++r){
        c = (r * (l->stride * l->c_len)) + (l->stride * i);
        for(j=0; j < l->alph; ++j){
            s_max += exp( ((l->lv_s[c+j] + log(pi[j])) - maximum) * SMPL_K );
        }
    }
    s_max = (maximum + (log(s_max)/SMPL_K)) * weight;
    return s_max;
}

/** [loglikelihood ml p prob %]
 * returns loglikelihood of a character set. |p| = l->stride
 *                                           |prob| = l->rates
 *                                           if pinvar > 0.0 then l-> rates == 1
 *
 * uses Kahan Summation to reduce floating point addition error. reduces error
 * of a summation to 2e
 */
double logMALlikelihood( const mll* l, const double* ws, const double* pi, 
                                const double* prob, const double pinvar )
{
#ifdef KAHANSUMMATION
    int i; 
    double ret_s,ret_c,ret_y,ret_t;

    ret_s = loglikelihood_site(l,ws[0],pi,prob,pinvar,0);
    ret_c = 0.0;
    for(i = 1;i<l->c_len; i++)
    {
        ret_y = (loglikelihood_site(l,ws[i],pi,prob,pinvar,i)) - ret_c;
        ret_t = ret_s + ret_y;
        ret_c = (ret_t - ret_s) - ret_y;
        ret_s = ret_t;
    }
    return ( -ret_s );
#else
    int i; double ret;
    ret = 0.0;
    for(i=0;i<l->c_len;i++)
        ret += loglikelihood_site(l,ws[i],pi,prob,pinvar,i);
    return ( -ret );
#endif
}

double logMPLlikelihood( const mll* l,const double* ws,const double* pi,const double* prob)
{
#ifdef KAHANSUMMATION
    int i; 
    double ret_s,ret_c,ret_y,ret_t;

    ret_s = logMPL_site(l,ws[0],pi,prob,0);
    ret_c = 0.0;
    for(i = 1;i<l->c_len; i++)
    {
        ret_y = (logMPL_site(l,ws[i],pi,prob,i)) - ret_c;
        ret_t = ret_s + ret_y;
        ret_c = (ret_t - ret_s) - ret_y;
        ret_s = ret_t;
    }
    return ( -ret_s );
#else
    int i;
    double ret;
    ret = 0.0;
    for(i=0;i<l->c_len;i++)
        ret += logMPL_site(l,ws[i],pi,prob,i);
    return ( -ret );
#endif
}

double logSMPLlikelihood( const mll* l,const double* ws,const double* pi,const double* prob)
{
#ifdef KAHANSUMMATION
    int i; 
    double ret_s,ret_c,ret_y,ret_t;

    ret_s = logSMPL_site(l,ws[0],pi,prob,0);
    ret_c = 0.0;
    for(i = 1;i<l->c_len; i++)
    {
        ret_y = (logSMPL_site(l,ws[i],pi,prob,i)) - ret_c;
        ret_t = ret_s + ret_y;
        ret_c = (ret_t - ret_s) - ret_y;
        ret_s = ret_t;
    }
    return ( -ret_s );
#else
    int i;
    double ret;
    ret = 0.0;
    for(i=0;i<l->c_len;i++)
        ret += logMPL_site(l,ws[i],pi,prob,i);
    return ( -ret );
#endif
}


double
loglikelihood( const mll* l,const double* ws,const double* pi,const double* prob,
                const double pinvar, const int cost)
{
    double total_cost;
    switch( cost ){
        case MALCOST : 
            total_cost = logMALlikelihood( l, ws, pi, prob, pinvar );
            break;
        case MPLCOST: 
            total_cost = logMPLlikelihood( l, ws, pi, prob );
            break;
        default :
            assert( FALSE );
    }
    return total_cost;
}

/* [likelihoood_CAML_loglikelihood s p] wrapper for loglikelihood */
value
likelihood_CAML_loglikelihood_wrapped(value s,value w,value pi,value prob,value pinvar,value mpl)
{
    CAMLparam5( s, w, pi, prob, pinvar);
    CAMLxparam1( mpl );
    CAMLlocal1( mle );
    double ll;
    ll = loglikelihood( ML_val(s), (double*) Data_bigarray_val(w),
                (double*) Data_bigarray_val(pi), (double*) Data_bigarray_val(prob),
                Double_val( pinvar ), Int_val(mpl));
    mle = caml_copy_double( ll );
    CAMLreturn( mle );
}

value
likelihood_CAML_loglikelihood( value * argv, int argn )
{
    return likelihood_CAML_loglikelihood_wrapped
        (argv[0], argv[1],argv[2],argv[3],argv[4],argv[5]);
}

/** Calculate the variance of the ratio of the sites of two vectors, presumably,
 * the leaf taxa-set is the same in each vector; this is the running variance
 * algorithm so only one pass is done on the vector. */
float variance_ratio(const mll* a, const mll* b, const double* ws,
                     const double* pi, const double* pr, const double pinvar)
{
    int h;
    double avg=0,var=0,a_h,b_h,v_h;

    for(h = 0; h < a->c_len; ++h){
        a_h = loglikelihood_site( a, ws[h], pi, pr, pinvar, h );
        b_h = loglikelihood_site( b, ws[h], pi, pr, pinvar, h );
        avg+= a_h - b_h;
    }
    avg = avg / a->c_len;
    for(h = 0; h < a->c_len; ++h){
        a_h = loglikelihood_site( a, ws[h], pi, pr, pinvar, h );
        b_h = loglikelihood_site( b, ws[h], pi, pr, pinvar, h );
        v_h = (a_h - b_h) - avg;
        var+= (v_h * v_h);
    }
    var = var * (a->c_len / (a->c_len-1));
    return var;
}

/** Calculate the variance of a likelihood vector **/
float variance( const mll* a, const double* ws, const double* pi,
                    const double* pr, const double pinvar)
{
    int h;
    double avg=0,var=0,a_h;

    for(h = 0; h < a->c_len; ++h)
        avg+= loglikelihood_site( a, ws[h], pi, pr, pinvar, h );
    avg = avg / a->c_len;
    for(h = 0; h < a->c_len; ++h){
        a_h = loglikelihood_site( a, ws[h], pi, pr, pinvar, h );
        var+= (a_h - avg) * (a_h - avg);
    }
    var = var * (a->c_len / (a->c_len-1));
    return var;
}

/** SSE functions for medians */
#ifdef __SSE3__

    #ifndef DEBUG
        #define load_next2(x) _mm_load_pd( x )
        #define load_next1(x) _mm_load1_pd( x )
        #define _mm_hmul_pd(x) _mm_mul_pd(x,_mm_shuffle_pd(x,x,_MM_SHUFFLE2(0,1)))

        #define load_next2u(x) _mm_loadu_pd( x )

    #else
        __m128d load_next2(const double* x)
        {
            printf("Loading, %f and %f\n", x[0], x[1]);
            return _mm_load_pd( x );
        }

        __m128d load_next1(const double* x)
        {
            printf("Loading, %f\n", x[0]);
            return _mm_load1_pd(x);
        }

        __m128d load_next2u(const double* x)
        {
            printf("Loading, %f and %f\n", x[0], x[1]);
            return _mm_loadu_pd(x);
        }

        __m128d load_next1p(const double* x)
        {
            printf("Loading, %f\n", x[0]);
            return _mm_load1_pd(x);
        }

        __m128d _mm_hmul_pd( __m128d x )
        {
            double* d;
            __m128d y;
            d = (double*) malloc( sizeof(double) * 4 );
            y = _mm_mul_pd(x,_mm_shuffle_pd(x,x,_MM_SHUFFLE2(0,1)));

            _mm_store_pd( d, x );
            _mm_store_pd( &d[2], y );
            printf("HMUL: %f, %f => [%f|%f]\n", d[0], d[1], d[2] , d[3]);
            return y;
        }
    #endif

    inline __m128d dotproduct2x2(__m128d a, __m128d b, __m128d pa, __m128d pb){
        a = _mm_mul_pd( a, pa );        b = _mm_mul_pd( b, pb );
        return _mm_hadd_pd( a, b );
    }
    inline __m128d dotproduct2x1(__m128d a, __m128d b, __m128d pa, __m128d pb){
         a = _mm_mul_sd( a, pa );       b = _mm_mul_sd( b, pb );
        return _mm_shuffle_pd( a, b, 0x00 );
    }

#endif

/** [median_c Pa Pb a b c]
 * Finds the likelihood vector [c] based on vectors of its children [a] and
 * [b] with the probability matrices [Pa] and [Pb], respectively, for a single
 * rate class [rate_idx]. [Pa] and [Pb] already contain rate/time information, and
 * the [rate_idx] serves to indicate the location in the likelihood vector.
 *
 * [tmp1] is size of alphabet.
 */
void
median_MAL(const double* Pa,const double* Pb, const mll* amll,const mll* bmll,
                mll* cmll, const int rate_idx )
{

#if defined( __SSE4__ )

    printf("Using SSE4 instructions");

#elif defined( __SSE3__ )

    /* GENERATION III -- SSE3; aligned data, unaligned P matrices */
    int i, j, k, c_loc, p_loc;
    __m128d acc;
    
    c_loc = rate_idx * amll->stride * amll->c_len;
    for(i=0; i < amll->c_len; ++i){
        p_loc = 0;
        for(j=0; j < amll->alph; ++j){
            acc = _mm_setzero_pd();
            for( k=0; k < (amll->alph-1); k+=2 ){
                acc = _mm_add_pd( acc, 
                            dotproduct2x2( 
                                load_next2( &amll->lv_s[(c_loc+k)]),
                                load_next2( &bmll->lv_s[(c_loc+k)]),
                                load_next2u(&Pa[(p_loc+k)]),
                                load_next2u(&Pb[(p_loc+k)]) ));
            }
            switch( amll->alph - k){
                case 1 :
                    acc = _mm_add_pd( acc, 
                                dotproduct2x1( 
                                    load_next1( &amll->lv_s[(c_loc+k)]),
                                    load_next1( &bmll->lv_s[(c_loc+k)]),
                                    load_next1(&Pa[(p_loc+k)]),
                                    load_next1(&Pb[(p_loc+k)]) ));
                case 0 :
                    break;
                default:
                    assert( 0 == 1 );
            }
            _mm_store_sd( &cmll->lv_s[c_loc+j], _mm_hmul_pd( acc ) );
            p_loc += amll->alph; //Ps to use stride; not aligned
            //printf("%d -> %f\n",(c_loc+j), cmll->lv_s[c_loc+j] );
        }
        c_loc += amll->stride;
    }

#else

    /* GENERATION II --inlined median_h, avoid recalc array access locations */
    // although, we can assume stride=alph, we don't because it doesn't add
    // complexity to the algorithm.
    int i, j, k, c_loc, p_loc;
    double tmp2, tmp1;
    c_loc = rate_idx * amll->stride * amll->c_len;
    for(i=0; i < amll->c_len; ++i){
        p_loc = 0;
        for(j=0; j < amll->alph; ++j){
            tmp1 = tmp2 = 0.0;
            for(k=0; k < amll->alph; ++k){
                tmp1 += Pa[p_loc+k] * amll->lv_s[c_loc+k];
                tmp2 += Pb[p_loc+k] * bmll->lv_s[c_loc+k];
            }
            cmll->lv_s[c_loc+j] = MAX( tmp1 * tmp2, 0.0 );
            p_loc += amll->alph;
            //printf("%d -> %f\n",(c_loc+j), cmll->lv_s[c_loc+j] );
        }
        c_loc += amll->stride;
    }
#endif

    /* GENERATION I -naive implementation
    int i,j,len,oth;
    double *tmp2, *tmp1;
    tmp1 = (double*) malloc( sizeof(double) * amll->stride );
    tmp2 = (double*) malloc( sizeof(double) * amll->stride );
    oth = rate_idx * amll->stride * amll->c_len;
    for(i=0; i < amll->c_len; ++i){
        len = oth + (i*amll->stride);
        tmp1[i] = tmp2[i] = 0.0;
        median_h( Pa, amll->lv_s, len, tmp1, amll->stride );
        median_h( Pb, bmll->lv_s, len, tmp2, amll->stride );
        for(j=0;j < amll->stride;++j)
            cmll->lv_s[len+j] = MAX((tmp2[j]*tmp1[j]),0);
    } */
}

void
median_MPL(const double* Pa,const double* Pb, const mll* amll,const mll* bmll,
                 mll* cmll, const int rate_idx )
{
#if defined( __SSE4__ )
    printf("Using SSE4 instructions");
#elif __SSE3__
    printf("Using SSE3 instructions");
#else
#endif

    int i, j, k, c_start, p_start;
    double tmp2, tmp1;
    c_start = rate_idx * amll->stride * amll->c_len;
    for(i=0; i < amll->c_len; ++i){
        p_start = 0;
        for(j=0; j < amll->stride; ++j){
            tmp1 = tmp2 = NEGINF;
            for(k=0; k < amll->stride; ++k){
                tmp1 = MAX (tmp1, log(Pa[p_start+k]) + amll->lv_s[c_start+k]);
                tmp2 = MAX (tmp2, log(Pb[p_start+k]) + bmll->lv_s[c_start+k]);
            }
            cmll->lv_s[c_start+j] = tmp1 + tmp2;
            p_start += amll->stride;
        }
        c_start += amll->stride;
    }
}

void 
median_invar(const mll* a, const mll* b, mll* c){
    int i;
    assert( (a->invar == b->invar) && (a->invar == 1) );
    for(i=0;i < c->c_len;++i)
        c->lv_invar[i] = (a->lv_invar[i]) & (b->lv_invar[i]);
}


#ifdef _WIN32
__inline
#else
inline
#endif
void
median(const double* PA, const double* PB, const mll* amll, const mll* bmll,
        mll* cmll, const int cost,const int rate_idx)
{
    switch( cost ){
        case MALCOST : 
            median_MAL( PA, PB, amll, bmll, cmll, rate_idx );
            break;
        case MPLCOST: 
            median_MPL( PA, PB, amll, bmll, cmll, rate_idx );
            break;
        default :
            assert( FALSE );
    }
}

value
likelihood_CAML_median2_wrapped_sym( value tmp, value U, value D, value ta,
        value tb, value ml_a, value ml_b, value rates, value mpl)
{
    /* ocaml macros */
    CAMLparam5( tmp, U, D, ta, tb );
    CAMLxparam4( ml_a, ml_b, rates, mpl );
    CAMLlocal1( ml_c );
    /* all variables */
    double cta,ctb,*c_U,*c_D,*PA,*PB,*tmp1,*g_rs;
    mll *a,*b,*c;
    mat *space;
    int num_rates,i;
    /* rates */
    num_rates = Bigarray_val(rates)->dim[0];
    g_rs = (double*) Data_bigarray_val( rates );
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
    assert( a->alph == b->alph );
    assert( a->c_len == b->c_len );
    assert( a->rates == b->rates );
    assert( a->invar == b->invar );
    /* expand scratch space for all function operations */
    space = FM_val( tmp );
    expand_matrix( space, 3 * (a->alph * a->alph) );
    /* create new character set, c */
    ml_c = caml_alloc_custom(&likelihood_custom_operations,(sizeof(mll*)),1,CAML_ALLOC);
    c = (mll*) malloc( sizeof(mll) );
    CHECK_MEM(c);
    ML_val( ml_c ) = c;
    c->stride = a->stride;
    c->alph  = a->alph;
    c->c_len = a->c_len;
    c->rates = a->rates;
    c->invar = a->invar;
    //printf ("rates:%d\tchars:%d\tstride:%d\n",a->rates,a->c_len,a->stride);
    lk_malloc(c->lv_s, c->rates * c->c_len * c->stride * sizeof(double));
    /* register temp variables */
    PA = register_section(space, a->alph*a->alph, 1);
    PB = register_section(space, a->alph*a->alph, 1);
    tmp1 = register_section( space, b->alph * b->alph, 1 );
    /** 
     *  F(x|r_i,P) = f_right(x*r_i|P) * f_left(x*r_i|P)
     *
     *      P    = probability matrix
     *      r_i  = the rate of i
     *      x    = the branch length
     */
    //printf ("rates:%d\tchars:%d\tstride:%d\n",a->rates,a->c_len,a->stride);
    for(i=0;i<num_rates;++i){
        compose_sym( PA, c_U, c_D, cta*g_rs[i], a->alph, tmp1 );
        compose_sym( PB, c_U, c_D, ctb*g_rs[i], b->alph, tmp1 );
        median( PA, PB, a, b, c, Int_val(mpl), i );
    }
    if(a->invar == 1){
        c->lv_invar = (int*) malloc ( c->c_len * sizeof(int));
        CHECK_MEM(c->lv_invar);
        median_invar( a,b,c );
    }
    //printarray(c->lv_s, a->stride * a->c_len * a->rates);
    /* free up space, return */
    assert( c == ML_val(ml_c));
    CAMLreturn(ml_c);
}

value likelihood_CAML_median2_wrapped_gtr( value tmp, value U, value D, value Ui,
        value ta, value tb, value ml_a, value ml_b, value rates, value mpl)
{
    /* ocaml macros */
    CAMLparam5( tmp,U,D,Ui,ta );
    CAMLxparam5( tb,ml_a,ml_b,rates,mpl );
    CAMLlocal1( ml_c );
    /* declare variables */
    double cta,ctb,*c_U,*c_D,*c_Ui,*PA,*PB,*tmp1,*g_rs;
    mll *a,*b,*c;
    mat *space;
    int num_rates,i;
    /* rates and probability vectors */
    num_rates = Bigarray_val(rates)->dim[0];
    g_rs= (double*) Data_bigarray_val(rates);
    /* diagonalized transition matrix */
    c_U = (double*) Data_bigarray_val( U );
    c_D = (double*) Data_bigarray_val( D );
    c_Ui= (double*) Data_bigarray_val( Ui);
    /* branch lengths */
    cta = Double_val( ta );
    ctb = Double_val( tb );
    assert( cta >= 0 && ctb >= 0 );
    /* character set */
    a = ML_val( ml_a );
    b = ML_val( ml_b );
    assert(a->stride == b->stride);
    assert(a->alph == b->alph);
    assert(b->c_len == a->c_len);
    assert(a->rates == b->rates);
    assert(a->invar == b->invar);
    /* expand scratch space for all function purposes */
    space = FM_val (tmp);
    expand_matrix(space, 3 * (a->alph * a->alph) );
    /* create new character set */
    ml_c = caml_alloc_custom(&likelihood_custom_operations, (sizeof(mll*)),1,CAML_ALLOC);
    c = (mll*) malloc( sizeof(mll) );
    CHECK_MEM(c);
    ML_val( ml_c ) = c;
    c->stride = a->stride;
    c->alph   = a->alph;
    c->c_len  = a->c_len;
    c->rates  = a->rates;
    c->invar  = a->invar;
    lk_malloc( c->lv_s, c->c_len * c->stride * c->rates * sizeof(double));
    /* register temp variables */
    PA = register_section( space, b->alph * b->alph, 1 );
    PB = register_section( space, b->alph * b->alph, 1 );
    tmp1 = register_section( space, b->alph * b->alph, 1 );
    /* main loop, see symmetric for description */
    //printf("R C S:%d\t%d\t%d\n[",a->rates,a->c_len,a->stride);
    for(i=0;i<num_rates;++i){
        compose_gtr( PA, c_U, c_D, c_Ui, cta*g_rs[i], a->alph, tmp1);
        compose_gtr( PB, c_U, c_D, c_Ui, ctb*g_rs[i], b->alph, tmp1);
        median( PA, PB, a,b,c, Int_val(mpl), i);
    }
    if(a->invar == 1){
        c->lv_invar = (int*) malloc( c->c_len * sizeof(int));
        CHECK_MEM(c->lv_invar);
        median_invar( a,b,c );
    }
    //printf("]\n");
    /* free all variables */
    assert( c == ML_val(ml_c) );
    CAMLreturn(ml_c);
}

/* [likelihood_CAML_median_sym ,,,] argument wrapper for median_sym */
value likelihood_CAML_median2_sym(value * argv, int argn)
{
    return likelihood_CAML_median2_wrapped_sym
        ( argv[0],argv[1],argv[2],argv[3],argv[4],argv[5],argv[6],argv[7],argv[8] );
}

/* [likelihood_CAML_median_gtr ...] argument wrapper for median_gtr */
value likelihood_CAML_median2_gtr(value * argv, int argn)
{
    return likelihood_CAML_median2_wrapped_gtr
        ( argv[0],argv[1],argv[2],argv[3],argv[4],argv[5],argv[6],argv[7],argv[8],argv[9] );
}

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

void
median3_MPL( const double* PA, const double* PB, const double* PC, 
             const mll* amll, const mll* bmll, const mll* cmll,
             mll* dmll, int rate_idx )
{
    int i, j, k, c_start, p_start;
    double tmp2, tmp1, tmp3;
    c_start = rate_idx * amll->stride * amll->c_len;
    for(i=0; i < amll->c_len; ++i){
        p_start = 0;
        for(j=0; j < amll->alph; ++j){
            tmp1 = tmp2 = tmp3 = NEGINF;
            for(k=0; k < amll->alph; ++k){
                tmp1 = MAX (tmp1, log(PA[p_start+k]) + amll->lv_s[c_start+k]);
                tmp2 = MAX (tmp2, log(PB[p_start+k]) + bmll->lv_s[c_start+k]);
                tmp3 = MAX (tmp3, log(PC[p_start+k]) + cmll->lv_s[c_start+k]);
            }
            dmll->lv_s[c_start+j] = tmp1 + tmp2 + tmp3;
            p_start += amll->alph;
        }
        c_start += amll->stride;
    }
}

#ifdef _WIN32
__inline
#else
inline
#endif
void
median3( const double* PA, const double* PB, const double* PC, 
             const mll* amll, const mll* bmll, const mll* cmll,
             mll* dmll, int mpl, int rate_idx )
{
    if( MPLCOST == mpl ){
        median3_MPL( PA,PB,PC,amll,bmll,cmll,dmll,rate_idx );
    } else if( MALCOST == mpl ){
        assert( FALSE );
    } else {
        assert( FALSE );
    }
}

void
median_3_sym(mat* space,const double* U,const double* D, const mll* amll,
        const mll* bmll, const mll* cmll, const double bl_a, const double bl_b, 
        const double bl_c, mll* dmll, const double* rates, const int num_rates,
        const int mpl)
{
    int i;
    double *PA, *PB, *PC, *TMP;

    /* allocate D */
    dmll->stride = amll->stride;
    dmll->c_len  = amll->c_len;
    dmll->rates  = amll->rates;
    dmll->invar = amll->invar;
    dmll->alph  = amll->alph;
    lk_malloc( dmll->lv_s, dmll->c_len * dmll->stride * dmll->rates * sizeof(double));

    /* determine space requirements */
    expand_matrix(space, 4 * (amll->alph * amll->alph) );
    PA = register_section( space, bmll->alph * bmll->alph, 1);
    PB = register_section( space, bmll->alph * bmll->alph, 1);
    PC = register_section( space, bmll->alph * bmll->alph, 1);
    TMP= register_section( space, bmll->alph * bmll->alph, 1);
    
    /* center of median three */
    for(i=0;i<num_rates;++i){
        compose_sym( PA, U, D, bl_a * rates[i], amll->alph, TMP);
        compose_sym( PB, U, D, bl_b * rates[i], amll->alph, TMP);
        compose_sym( PC, U, D, bl_c * rates[i], amll->alph, TMP);
        median3( PA, PB, PC, amll,bmll,cmll,dmll, mpl, i );
    }
    if( amll->invar == 1 ){
        dmll->lv_invar = (int*) malloc(dmll->c_len * sizeof(int));
        CHECK_MEM(dmll->lv_invar);
        //median3_invar(amll,bmll,cmll,dmll);
    }
}

void
median_3_gtr(mat* space,const double* U,const double* D, const double *Ui,
        const mll* amll, const mll* bmll, const mll* cmll, const double bl_a,
        const double bl_b, const double bl_c, mll* dmll, const double* rates,
        const int num_rates, const int mpl)
{
    int i;
    double *PA, *PB, *PC, *TMP;

    /* allocate D */
    dmll->stride = amll->stride;
    dmll->alph = amll->alph;
    dmll->c_len  = amll->c_len;
    dmll->rates  = amll->rates;
    dmll->invar = amll->invar;
    lk_malloc( dmll->lv_s, dmll->c_len * dmll->stride * dmll->rates * sizeof(double));

    /* determine space requirements */
    expand_matrix(space, 4 * (amll->stride * amll->stride) );
    PA = register_section( space, bmll->alph * bmll->alph, 1);
    PB = register_section( space, bmll->alph * bmll->alph, 1);
    PC = register_section( space, bmll->alph * bmll->alph, 1);
    TMP= register_section( space, bmll->alph * bmll->alph, 1);
    
    /* center of median three */
    for(i=0;i<num_rates;++i){
        compose_gtr( PA, U, D, Ui, bl_a * rates[i], amll->alph, TMP);
        compose_gtr( PB, U, D, Ui, bl_b * rates[i], amll->alph, TMP);
        compose_gtr( PC, U, D, Ui, bl_c * rates[i], amll->alph, TMP);
        median3( PA, PB, PC, amll,bmll,cmll,dmll, mpl, i );
    }
    if( amll->invar == 1 ){
        dmll->lv_invar = (int*) malloc(dmll->c_len * sizeof(int));
        CHECK_MEM(dmll->lv_invar);
        //median3_invar(amll,bmll,cmll,dmll);
    }
}

value
likelihood_CAML_median3_wrapped_gtr( value tmp, value U, value D, value Ui, value ta,
        value tb, value tc, value mlla, value mllb, value mllc, value rates, value mpl)
{
    CAMLparam5(tmp,U,D,ta,tb);
    CAMLxparam5(tc,mlla,mllb,mllc,rates);
    CAMLxparam1(mpl);
    CAMLlocal2( res, ml_d );
    double cta,ctb,ctc;
    int n_rates;

    n_rates = Bigarray_val(rates)->dim[0];
    mll *d;
    d = (mll*) malloc( sizeof(mll) ); CHECK_MEM(d);
    cta = Double_val( ta );
    ctb = Double_val( tb );
    ctc = Double_val( tc );

    median_3_gtr( FM_val(tmp), Data_bigarray_val(U), Data_bigarray_val(D), 
            Data_bigarray_val(Ui), ML_val(mlla), ML_val(mllb), ML_val(mllc),
            cta, ctb, ctc, d, Data_bigarray_val(rates), n_rates, Int_val(mpl) );

    ml_d = caml_alloc_custom(&likelihood_custom_operations,sizeof(mll*),1,CAML_ALLOC);
    ML_val( ml_d ) = d;
    CAMLreturn(ml_d);
}

value 
likelihood_CAML_median3_gtr(value * argv, int argn)
{
    return likelihood_CAML_median3_wrapped_gtr
        ( argv[0],argv[1],argv[2],argv[3],argv[4],argv[5],argv[6],argv[7],
          argv[8],argv[9],argv[10],argv[11] );
}


value
likelihood_CAML_median3_wrapped_sym( value tmp, value U, value D, value ta,
        value tb, value tc, value mlla, value mllb, value mllc, value rates, value mpl)
{
    CAMLparam5(tmp,U,D,ta,tb);
    CAMLxparam5(tc,mlla,mllb,mllc,rates);
    CAMLxparam1(mpl);
    CAMLlocal2( res, ml_d );
    double cta,ctb,ctc;
    int n_rates;

    n_rates = Bigarray_val(rates)->dim[0];
    mll *d;
    d = (mll*) malloc( sizeof(mll) ); CHECK_MEM(d);
    cta = Double_val( ta );
    ctb = Double_val( tb );
    ctc = Double_val( tc );

    median_3_sym( FM_val(tmp), Data_bigarray_val(U), Data_bigarray_val(D), 
            ML_val(mlla), ML_val(mllb), ML_val(mllc), cta, ctb, ctc,
            d, Data_bigarray_val(rates), n_rates, Int_val(mpl) );

    ml_d = caml_alloc_custom(&likelihood_custom_operations,sizeof(mll*),1,CAML_ALLOC);
    ML_val( ml_d ) = d;
    CAMLreturn(ml_d);
}

value 
likelihood_CAML_median3_sym(value * argv, int argn)
{
    return likelihood_CAML_median3_wrapped_sym
        ( argv[0],argv[1],argv[2],argv[3],argv[4],argv[5],argv[6],argv[7],
          argv[8],argv[9],argv[10] );
}

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
void
median1_MPL(const double* PA,const mll* amll,const mll* bmll,mll* dmll,int rate_idx)
{
    int i, j, k, c_start, p_start;
    double tmp;
    c_start = rate_idx * amll->stride * amll->c_len;
    for(i=0; i < amll->c_len; ++i){
        p_start = 0;
        for(j=0; j < amll->stride; ++j){
            tmp = NEGINF;
            for(k=0; k < amll->stride; ++k){
                tmp = MAX (tmp, log(PA[p_start+k]) + amll->lv_s[c_start+k]);
            }
            dmll->lv_s[c_start+j] = tmp + bmll->lv_s[c_start+j];
            p_start += amll->stride;
        }
        c_start += amll->stride;
    }
}

#ifdef _WIN32
__inline
#else
inline
#endif
void
median1( const double* PA, const mll* amll, const mll* bmll, 
             mll* dmll, int cost, int rate_idx )
{
    switch( cost ){
        case MPLCOST: 
            median1_MPL( PA,amll,bmll,dmll,rate_idx );
            break;
        default :
            assert( FALSE );
    }
}


void
median_1_sym(mat* space, const double* U, const double* D, const mll* amll,
        const mll* bmll, const double bl_a, mll* dmll, const double* rates,
        const int num_rates, const int mpl)
{
    int i;
    double *PA, *TMP;

    /* allocate D */
    dmll->stride = amll->stride;
    dmll->alph = amll->alph;
    dmll->c_len  = amll->c_len;
    dmll->rates  = amll->rates;
    dmll->invar = amll->invar;
    lk_malloc( dmll->lv_s, dmll->c_len * dmll->stride * dmll->rates * sizeof(double));

    /* determine space requirements */
    expand_matrix(space, 4 * (amll->alph * amll->alph) );
    PA = register_section( space, bmll->alph * bmll->alph, 1);
    TMP= register_section( space, bmll->alph * bmll->alph, 1);
    
    /* center of median three */
    for(i=0;i<num_rates;++i){
        compose_sym( PA, U, D, bl_a * rates[i], amll->alph, TMP);
        median1( PA, amll, bmll, dmll, mpl, i );
    }
    if( amll->invar == 1 ){
        dmll->lv_invar = (int*) malloc(dmll->c_len * sizeof(int));
        CHECK_MEM(dmll->lv_invar);
        //median1_invar(amll,bmll,dmll);
    }
}

void
median_1_gtr(mat* space, const double* U, const double* D, const double* Ui,
        const mll* amll, const mll* bmll, const double bl_a, mll* dmll,
        const double* rates, const int num_rates, const int mpl)
{
    int i;
    double *PA, *TMP;

    /* allocate D */
    dmll->stride = amll->stride;
    dmll->alph = amll->alph;
    dmll->c_len  = amll->c_len;
    dmll->rates  = amll->rates;
    dmll->invar = amll->invar;
    lk_malloc( dmll->lv_s, dmll->c_len * dmll->stride * dmll->rates * sizeof(double));

    /* determine space requirements */
    expand_matrix(space, 4 * (amll->alph * amll->alph) );
    PA = register_section( space, bmll->alph * bmll->alph, 1);
    TMP= register_section( space, bmll->alph * bmll->alph, 1);
    
    /* center of median three */
    for(i=0;i<num_rates;++i){
        compose_gtr( PA, U, D, Ui, bl_a * rates[i], amll->alph, TMP);
        median1( PA, amll,bmll,dmll, mpl, i );
    }
    if( amll->invar == 1 ){
        dmll->lv_invar = (int*) malloc(dmll->c_len * sizeof(int));
        CHECK_MEM(dmll->lv_invar);
        //median1_invar(amll,bmll,dmll);
    }
}


value
likelihood_CAML_median1_wrapped_gtr( value tmp, value U, value D, value Ui, value ta,
        value mlla, value mllb, value rates, value mpl)
{
    CAMLparam5(tmp,U,D,Ui,ta);
    CAMLxparam4(mlla,mllb,rates,mpl);
    CAMLlocal2( res, ml_d );
    double cta;
    int n_rates;

    n_rates = Bigarray_val(rates)->dim[0];
    mll *d;
    d = (mll*) malloc( sizeof(mll) ); CHECK_MEM(d);
    cta = Double_val( ta );

    median_1_gtr( FM_val(tmp), Data_bigarray_val(U), Data_bigarray_val(D), 
            Data_bigarray_val(Ui), ML_val(mlla), ML_val(mllb),
            cta, d, Data_bigarray_val(rates), n_rates, Int_val(mpl) );

    ml_d = caml_alloc_custom(&likelihood_custom_operations,sizeof(mll*),1,CAML_ALLOC);
    ML_val( ml_d ) = d;
    CAMLreturn(ml_d);
}

value 
likelihood_CAML_median1_gtr(value * argv, int argn)
{
    return likelihood_CAML_median1_wrapped_gtr
        (argv[0],argv[1],argv[2],argv[3],argv[4],argv[5],argv[6],argv[7],argv[8]);
}


value
likelihood_CAML_median1_wrapped_sym( value tmp, value U, value D, value ta,
        value mlla, value mllb, value rates, value mpl)
{
    CAMLparam5(tmp,U,D,ta,mlla);
    CAMLxparam3(mllb,rates,mpl);
    CAMLlocal2( res, ml_d );
    double cta;
    int n_rates;

    n_rates = Bigarray_val(rates)->dim[0];
    mll *d;
    d = (mll*) malloc( sizeof(mll) ); CHECK_MEM(d);
    cta = Double_val( ta );

    median_1_sym( FM_val(tmp), Data_bigarray_val(U), Data_bigarray_val(D), 
            ML_val(mlla), ML_val(mllb), cta, d, Data_bigarray_val(rates), 
            n_rates, Int_val(mpl) );

    ml_d = caml_alloc_custom(&likelihood_custom_operations,sizeof(mll*),1,CAML_ALLOC);
    ML_val( ml_d ) = d;
    CAMLreturn(ml_d);
}

value 
likelihood_CAML_median1_sym(value * argv, int argn)
{
    return likelihood_CAML_median1_wrapped_sym
        ( argv[0],argv[1],argv[2],argv[3],argv[4],argv[5],argv[6],argv[7] );
}

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
void single_sym(ptr *simp, double *PA, double *PB, const double *U, const double *D,
    const mll *a, const mll *b, const double time_a, const double time_b,
    const double *ws, const double *gam, const double* prob, const double *pi, 
    const int g_n, const double percent, const int mpl, double *tmp)
{ 
    int i;
    simp->time = time_a;
    for(i=0;i<g_n;++i){
        compose_sym( PA, U, D, time_a*gam[i], a->alph, tmp );
        compose_sym( PB, U, D, time_b*gam[i], b->alph, tmp );
        median( PA, PB, a, b, simp->vs, mpl, i );
    }
    /* invar isn't iterated
     * if(a->invar == 1) median_invar( a,b,simp->vs ); */
    simp->ll = loglikelihood( simp->vs, ws, pi, prob, percent, mpl );
}
void single_gtr(ptr *simp, double *PA, double *PB, const double *U, const double *D,
    const double *Ui, const mll *a, const mll *b, const double time_a, 
    const double time_b, const double* ws, const double* gam, const double *prob,
    const double *pi, const int g_n, const double percent, const int mpl, double *tmp)
{
    int i;
    simp->time = time_a;
    for(i=0;i<g_n;++i){
        compose_gtr( PA, U, D, Ui, time_a*gam[i], a->alph, tmp );
        compose_gtr( PB, U, D, Ui, time_b*gam[i], b->alph, tmp );
        median( PA, PB, a, b, simp->vs, mpl, i );
    }
    /* invar isn't iterated
     * if(a->invar == 1) median_invar( a,b,simp->vs ); */
    simp->ll = loglikelihood( simp->vs, ws, pi, prob, percent, mpl );
}
#define golden (2.0/(sqrt(5.0)+1))
#define golden_exterior_l(b,c) ((b - (c * golden)) / (1 - golden))
#define golden_exterior_r(a,b) (a + ((b - a) / golden))

void
readjust_brents_sym(mat *space,const double* Um,const double* D,const mll* data_c1,
        const mll* data_c2,mll* data_p,double* b_tc1,const double b_tc2,double* b_mle,
        const double* ws,const double* rates,const double *prob,const int g_n,
        const double* pi,const double pinvar, const int mpl)
{
    int iter,size,bracketed;
    ptr temp; mll temp_;
    double x,w,v,u,fx,fw,fv,fu,xm,a,b,fa,fb,fi,pu;
    double *PA,*PB,*TMP;
    double r,q,p,d,e,tol,tmp;

    size = data_c1->c_len * data_c1->stride * data_c1->rates;
    /* register temp space and matrices --1 vector, 3 matrices */
    expand_matrix( space, (size) + (3*data_c1->stride*data_c1->stride) );
    PA  = register_section(space, data_c1->stride*data_c1->stride, 0);
    PB  = register_section(space, data_c1->stride*data_c1->stride, 0);
    TMP = register_section(space, data_c1->stride*data_c1->stride, 0);
    (temp.vs)           = &temp_;
    (temp.vs)->lv_s     = register_section( space, size, 0 );
    (temp.vs)->stride   = data_c1->stride;
    (temp.vs)->c_len    = data_c1->c_len;
    (temp.vs)->invar    = data_c1->invar;
    (temp.vs)->alph     = data_c1->alph;
    (temp.vs)->lv_invar = data_c1->lv_invar;
    (temp.vs)->rates    = data_c1->rates;

    e = d = pu = 0.0;
    /* initial bracket */
    x  = *b_tc1;
    fx = fi = *b_mle;
    a = MAX( BL_MIN, x / 10.0 );
    b = MIN( BL_MAX, x * 10.0 );
    single_sym(&temp,PA,PB,Um,D,data_c1,data_c2,a,b_tc2,ws,rates,prob,pi,g_n,pinvar,mpl,TMP);
    fa = temp.ll;
    single_sym(&temp,PA,PB,Um,D,data_c1,data_c2,b,b_tc2,ws,rates,prob,pi,g_n,pinvar,mpl,TMP);
    fb = temp.ll;
    bracketed = 1;
    /* bracket a region; a and b become the bracket with temp */
    iter = 0;
    while( iter <= MAX_ITER ){
        //printf("Bracketing: %f(%f)\t%f(%f)\t%f(%f)\n",a,fa,x,fx,b,fb);
        /* bracketed */
        if( (fa > fx) && (fx < fb)){
            break;
        /* check if they are equal; the branch is a nussance paramater */
        } else if( EQUALS(fa,fx) && EQUALS(fb,fx) ){
            bracketed = 0; break;
        /* increasing; move left: (new,left,middle) */
        } else if (fa < fx && fx < fb){
            b = x; fb = fx;
            x = a; fx = fa;
            a = MAX( BL_MIN, golden_exterior_l( x,b ) );
            //assert ( a < x || EQUALS(a,x) );
            single_sym(&temp,PA,PB,Um,D,data_c1,data_c2,a,b_tc2,ws,rates,prob,pi,g_n,pinvar,mpl,TMP);
            fa = temp.ll;
        /* decreasing; move right: (middle,right,new) */
        } else if( fa > fx && fx > fb ){ 
            a = x; fa = fx;
            x = b; fx = fb;
            b = MAX (BL_MIN, golden_exterior_r( a, x ) );
            //assert ( b > x || EQUALS(b,x) );
            single_sym(&temp,PA,PB,Um,D,data_c1,data_c2,a,b_tc2,ws,rates,prob,pi,g_n,pinvar,mpl,TMP);
            fb = temp.ll;
        /* convergence or curvature issue */
        } else {
            bracketed = 0; break;
        }
        ++iter;
    }

    /* if( iter >= MAX_ITER ){ printf("HIT MAX COUNT IN BRENT w/ BRACKETING!\n"); }*/
    /* printf("Bracketed(%d): %f(%f)\t%f(%f)\t%f(%f)\n",bracketed,a,fa,x,fx,b,fb); */
    /* we have a bracketed region; (a < b < c) && (fa > fb < fc) */
    v=w=x;
    fv=fw=fx;
    while( (iter <= MAX_ITER) && (1 == bracketed)){
        xm = 0.5 * (a + b);
        tol = BRENT_TOL * fabs(x) + EPSILON;
        /* determine convergence; x contains best score */
        if( fabs(x - xm) <= ((2.0*tol) - 0.5 * (b - a)) ){ break; }

        /* construct parabolic fit */
        if(fabs(e) > tol){
            r = (x-w) * (fx-fv);
            q = (x-v) * (fx-fw);
            p = (x-v)*q-(x-w)*r;
            q = 2.0*(q-r);
            if( q > 0.0){ p = -p; }
            q = fabs(q);
            tmp = e;
            e = d; 
            /* determine acceptability of fit; else do golden section search */
            if((fabs(p) >= fabs(0.5*q*tmp)) || (p <= q*(a-x)) || (p >= q*(b-x))){
                /* printf("\tGolden Section\n");*/
                e = (x >= xm) ? (a-x) : (b-x);
                d = golden * e;
            } else {
                /* printf("\tParabolic Fit\n");*/
                d = p / q;
                u = x + d;
                if( (u-a < 2.0*tol) || (b-u < 2.0*tol) ){
                    d = SIGN(tol,xm-x);
                }
            }
            /* golden section search */
        } else {
            /* printf("\tGolden Section\n");*/
            e = (x >= xm) ? (a-x) : (b-x);
            d = golden * e;
        }
        /* function evaluation; replace data in u. */
        u = (fabs(d) >= tol) ? MAX( BL_MIN, (x+d) ) : MAX( BL_MIN, (x+SIGN(tol,d)) );
        single_sym(&temp,PA,PB,Um,D,data_c1,data_c2,u,b_tc2,ws,rates,prob,pi,g_n,pinvar,mpl,TMP);
        fu = temp.ll;
        /** printf("\tIteration(%d): %f(%f)\t[%f(%f)]\t%f(%f)\n",iter,a,fa,u,fu,b,fb); **/

        /* Determine Convergence to a singularity */
        if( pu == u ){
            if( fx < fu ) { break; }
                     else { fx = fu; x = u; break;}
        }
        /* move variables around for next motion */
        if( fu <= fx ){
            if( u >= x ){ a = x; fa = fx; } else { b = x; fb = fx; }
            SHIFT(v,w,x,u);
            SHIFT(fv,fw,fx,fu);
        } else {
            if( u < x ){ a = u; fa = fu; } else { b = u; fb = fu; }
            if( (fu <= fw) || (w == x) ){
                v=w; fv=fw;
                w=u; fw=fu;
            } else if( (fu <= fv) || (v == x) || (v==w) ){
                v=u; fv=fu;
            }
        }
        pu = u;
        ++iter;
    }
    /* if( iter >= MAX_ITER ){ printf("\tHIT MAX COUNT IN BRENT!\n"); }*/
    memcpy( data_p->lv_s, (temp.vs)->lv_s, size * sizeof(double));
    *b_tc1 = x;
    *b_mle = fx;
}

void
readjust_brents_gtr(mat * space,const double* Um,const double* D,const double* Ui,
        const mll* data_c1,const mll* data_c2,mll* data_p,double* b_tc1,const double b_tc2,
        double* b_mle,const double* ws,const double *rates,const double *prob,
        const int g_n,const double* pi,const double pinvar,const int mpl)
{
    int iter,size,bracketed;
    ptr temp; mll temp_;
    double x,w,v,u,fx,fw,fv,fu,xm,a,b,fa,fb,fi,pu;
    double *PA,*PB,*TMP;
    double r,q,p,d,e,tol,tmp;

    size = data_c1->c_len * data_c1->stride * data_c1->rates;
    /* register temp space and matrices --1 vector, 3 matrices */
    expand_matrix( space, (size) + (3*data_c1->stride*data_c1->stride) );
    PA  = register_section(space, data_c1->stride*data_c1->stride, 0);
    PB  = register_section(space, data_c1->stride*data_c1->stride, 0);
    TMP = register_section(space, data_c1->stride*data_c1->stride, 0);
    (temp.vs)           = &temp_;
    (temp.vs)->lv_s     = register_section( space, size, 0 );
    (temp.vs)->stride   = data_c1->stride;
    (temp.vs)->c_len    = data_c1->c_len;
    (temp.vs)->invar    = data_c1->invar;
    (temp.vs)->alph     = data_c1->alph;
    (temp.vs)->lv_invar = data_c1->lv_invar;
    (temp.vs)->rates    = data_c1->rates;

    e = d = pu = 0.0;
    /* initial bracket */
    x = *b_tc1;
    fx = fi = *b_mle;
    a = MAX( BL_MIN, x / 100.0 );
    b = MIN( BL_MAX, x * 10.0 );
    single_gtr(&temp,PA,PB,Um,D,Ui,data_c1,data_c2,a,b_tc2,ws,rates,prob,pi,g_n,pinvar,mpl,TMP);
    fa = temp.ll;
    single_gtr(&temp,PA,PB,Um,D,Ui,data_c1,data_c2,b,b_tc2,ws,rates,prob,pi,g_n,pinvar,mpl,TMP);
    fb = temp.ll;
    bracketed = 1;
    /* bracket a region; a and b become the bracket with temp */
    iter = 0;
    while( iter <= MAX_ITER ){
    /* printf("Bracketing: %f(%f)\t%f(%f)\t%f(%f)\n",a,fa,x,fx,b,fb);*/
        /* bracketed */
        if( (fa > fx) && (fx < fb)){
            break;
        /* increasing; move left: (new,left,middle) */
        } else if (fa < fx && fx < fb){
            b = x; fb = fx;
            x = a; fx = fa;
            a = MAX( BL_MIN, golden_exterior_l( x,b ) );
            assert ( a <= x );
            single_gtr(&temp,PA,PB,Um,D,Ui,data_c1,data_c2,a,b_tc2,ws,rates,prob,pi,g_n,pinvar,mpl,TMP);
            fa = temp.ll;
        /* decreasing; move right: (middle,right,new) */
        } else if( fa > fx && fx > fb ){ 
            a = x; fa = fx;
            x = b; fx = fb;
            b = MAX (BL_MIN, golden_exterior_r( a, x ) );
            assert ( b >= x );
            single_gtr(&temp,PA,PB,Um,D,Ui,data_c1,data_c2,a,b_tc2,ws,rates,prob,pi,g_n,pinvar,mpl,TMP);
            fb = temp.ll;
        /* convergence or curvature issue */
        } else {
            bracketed = 0; break;
        }
        ++iter;
    }

    /* if( iter >= MAX_ITER ){ printf("HIT MAX COUNT IN BRENT w/ BRACKETING!\n"); }*/
    /* printf("Bracketed(%d): %f(%f)\t%f(%f)\t%f(%f)\n",bracketed,a,fa,x,fx,b,fb);*/
    /* we have a bracketed region; (a < b < c) && (fa > fb < fc) */
    v=w=x;
    fv=fw=fx;
    while( (iter <= MAX_ITER) && (1 == bracketed)){
        xm = 0.5 * (a + b);
        tol = BRENT_TOL * fabs(x) + EPSILON;
        /* determine convergence; x contains best score */
        if( fabs(x - xm) <= ((2.0*tol) - 0.5 * (b - a)) ){ break; }

        /* construct parabolic fit */
        if(fabs(e) > tol){
            r = (x-w) * (fx-fv);
            q = (x-v) * (fx-fw);
            p = (x-v)*q-(x-w)*r;
            q = 2.0*(q-r);
            if( q > 0.0){ p = -p; }
            q = fabs(q);
            tmp = e;
            e = d; 
            /* determine acceptability of fit; else do golden section search */
            if((fabs(p) >= fabs(0.5*q*tmp)) || (p <= q*(a-x)) || (p >= q*(b-x))){
                /* printf("\tGolden Section\n");*/
                e = (x >= xm) ? (a-x) : (b-x);
                d = golden * e;
            } else {
                /* printf("\tParabolic Fit\n");*/
                d = p / q;
                u = x + d;
                if( (u-a < 2.0*tol) || (b-u < 2.0*tol) ){
                    d = SIGN(tol,xm-x);
                }
            }
            /* golden section search */
        } else {
            /* printf("\tGolden Section\n");*/
            e = (x >= xm) ? (a-x) : (b-x);
            d = golden * e;
        }
        /* function evaluation; replace data in u. */
        u = (fabs(d) >= tol) ? MAX( BL_MIN, (x+d) ) : MAX( BL_MIN, (x+SIGN(tol,d)) );
        single_gtr(&temp,PA,PB,Um,D,Ui,data_c1,data_c2,u,b_tc2,ws,rates,prob,pi,g_n,pinvar,mpl,TMP);
        fu = temp.ll;
        /* printf("\tIteration(%d): %f(%f)\t[%f(%f)]\t%f(%f)\n",iter,a,fa,u,fu,b,fb);*/
        /* Determine Convergence to a singularity */
        if( pu == u ){
            if( fx < fu ) { break; }
                     else { fx = fu; x = u; break;}
        }
        /* move variables around for next motion */
        if( fu <= fx ){
            if( u >= x ){ a = x; } else { b = x; }
            SHIFT(v,w,x,u);
            SHIFT(fv,fw,fx,fu);
        } else {
            if( u < x ){ a = u; } else { b = u; }
            if( (fu <= fw) || (w == x) ){
                v=w; fv=fw;
                w=u; fw=fu;
            } else if( (fu <= fv) || (v == x) || (v==w) ){
                v=u; fv=fu;
            }
        }
        ++iter;
        pu = u;
    }
    /* if( iter >= MAX_ITER ){ printf("\tHIT MAX COUNT IN BRENT!\n"); }*/
    memcpy( data_p->lv_s, (temp.vs)->lv_s, size * sizeof(double));
    *b_tc1 = x;
    *b_mle = fx;
}

//-----------------------------------------------------------------------------
/** [..._readjust_gtr... tmp U D Ui A B C ta tb gammas probs pis ll] **/
value
likelihood_CAML_readjust_gtr_wrapped
    (value tmp,value U,value D,value Ui,value A,value B,value C,value ta,value tb,
        value pinv, value ws, value g,value p,value pis,value ll,value mpl)
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
    assert( (ML_val(A)->alph) == Bigarray_val(pis)->dim[0] );
    assert( Bigarray_val(g)->dim[0] == Bigarray_val(p)->dim[0] );
    /* readjust */
    //printf("S:%f\t%f\t",cta,,likelihood);
    readjust_brents_gtr(
            FM_val(tmp),
            Data_bigarray_val( U ),
            Data_bigarray_val( D ),
            Data_bigarray_val( Ui ),
            ML_val( A ), ML_val( B ), ML_val( C ),
            &cta, Double_val( tb ), &likelihood,
            Data_bigarray_val(ws),
            Data_bigarray_val(g),
            Data_bigarray_val(p),
            Bigarray_val( g )->dim[0],
            Data_bigarray_val(pis),
            Double_val( pinv ),
            Int_val(mpl)
    );
    //printf("E:%f\t\t%f\n",cta,likelihood);
    /* construct tuple / freespace */
    res = caml_alloc_tuple( 2 );
    Store_field(res, 0, caml_copy_double(cta));
    Store_field(res, 1, caml_copy_double(likelihood));
    CAMLreturn(res);
}
value
likelihood_CAML_readjust_sym_wrapped
    (value tmp,value U,value D,value A,value B,value C,value ta,value tb,
        value pinv,value ws,value g,value p,value pis,value ll,value mpl)
{
    /* ocaml macros for GC */
    CAMLparam5(U,D,ll,A,B);
    CAMLxparam5(C,ta,tb,g,p);
    CAMLxparam3( pis,tmp,ws );
    CAMLlocal1( res );
    /* declare/unbox variables to return */
    double cta,likelihood;
    cta = Double_val( ta );
    likelihood = Double_val( ll );
    /* check probabilities, rates, priors */
    assert( ML_val(A)->alph == Bigarray_val( pis )->dim[0]);
    assert( Bigarray_val( g )->dim[0] == Bigarray_val( p )->dim[0] );
    /* readjust loop */
    //printf("\tS:%f\t%f\t%f\n",cta,Double_val(tb),likelihood);
    readjust_brents_sym(
            FM_val(tmp),
            Data_bigarray_val( U ), Data_bigarray_val( D ),
            ML_val(A), ML_val(B), ML_val(C),
            &cta, Double_val( tb ), &likelihood,
            Data_bigarray_val( ws),
            Data_bigarray_val( g ),
            Data_bigarray_val( p ),
            Bigarray_val( g )->dim[0],
            Data_bigarray_val( pis ),
            Double_val( pinv ),
            Int_val(mpl)
    );
    //printf("\tE:%f\t%f\t%f\n",cta,Double_val(tb),likelihood);
    /* construct tuple / free space */
    res = caml_alloc_tuple( 2 );
    Store_field(res, 0, caml_copy_double (cta));
    Store_field(res, 1, caml_copy_double (likelihood));
    CAMLreturn(res);
}

/** functions for interface **/
value
likelihood_CAML_readjust_gtr(value * argv, int argn)
{
    return likelihood_CAML_readjust_gtr_wrapped
        (argv[0], argv[1],argv[2],argv[3],argv[4],argv[5],argv[6],argv[7],
         argv[8],argv[9],argv[10],argv[11],argv[12],argv[13],argv[14],argv[15]);
}
value
likelihood_CAML_readjust_sym(value * argv, int argn)
{
    return likelihood_CAML_readjust_sym_wrapped
        (argv[0], argv[1],argv[2],argv[3],argv[4],argv[5],argv[6],
         argv[7],argv[8],argv[9],argv[10],argv[11],argv[12],argv[13],argv[14]);
}

#else

double BRENT_TOL = 1e-6;
value likelihood_CAML_set_tol( value a ){
    CAMLparam1( a );
    BRENT_TOL = Double_val( a );
    CAMLreturn( Val_unit );
}

/** Get the TOLERANCE used in the optimization functions for the convergence of
 * the brents method functions used in optimizing the branch lengths. */
value likelihood_CAML_get_tol( value unit ){
    CAMLparam1( unit );
    CAMLlocal1( brent_tol );
    brent_tol = caml_copy_double( BRENT_TOL );
    CAMLreturn( brent_tol );
}

#endif /* USE_LIKELIHOOD */
