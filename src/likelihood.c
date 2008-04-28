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
#include <stdlib.h> //malloc, calloc, srand, RAND_MAX
#include <string.h> //memcpy, memset
#include <assert.h>

#include "config.h"
#ifdef USE_LIKELIHOOD   
#include <math.h>   //log,exp

#include <caml/mlvalues.h>
#include <caml/memory.h>
#include <caml/bigarray.h>
#include <caml/custom.h>
#include <caml/intext.h>
#include <caml/fail.h>

#include "array_pool.h"
#include "likelihood.h"

//decent values (time/accuracy) for gamma calculations
#define STEP 1e-4     //for numerical integration
#define EPSILON 1e-4  //error for numerical calculations
#define MAX_ITER 10000
#define ML_ptr(v)   ((struct ml**)Data_custom_val(a))
#define ML_val(v) (*((struct ml**)Data_custom_val(v)))

#define CAML_ALLOC_2 1000000

//------------------------------------------------------------------------------
/* assert checking 
 *
 * ~CHECK_MEM       -- verify enough memory
 * ~CHECK_POSITIVE  -- verify array is all positive values with EPSILON error
 * ~CHECK_ZEROES    -- verify array is all zeroes with EPSILON error
 * ~CHECK_MEAN      -- verify the mean of the gamma rates == 1
 */

/*  asserts  */
#define CHECK_MEM(a) if(a==NULL) failwith("I can't allocate more memory.")
//#define CHECK_POSITIVE(a,n); int Z;for(Z=0;Z<n;Z++){ if( a[Z] < -EPSILON)\
                            failwith("Negative Likelihood"); }
#define CHECK_ZEROES(a,n); int Z;for(Z=0;Z<n;Z++){ if(a[Z] > EPSILON || a[Z] < -EPSILON)\
                            failwith("Imaginary eigenvalue"); }
//#define CHECK_MEAN(a,n); int Z;double SUM;for(Z=0,SUM=0;Z<n;Z++){ SUM += a[Z]; }\
                            if( SUM/(double)n > 1.0+EPSILON ){ \
                            failwith("Incorrect Mean of Gamma Rates"); }

/*  printline  */
//#define CHECK_MEM(a) if(a==NULL) printf("I can't allocate more memory. %d",__LINE__)
//#define CHECK_POSITIVE(a,n); int Z;for(Z=0;Z<n;Z++){ if( a[Z] < -EPSILON){ printf("Negative Likelihood :: %f\n",a[Z]);} }
//#define CHECK_ZEROES(a,n); int Z;for(Z=0;Z<n;Z++){ if(a[Z] > EPSILON || a[Z] < -EPSILON){ printf("Imaginary eigenvalue :: %f\n", a[Z]); } }
//#define CHECK_MEAN(a,n); int Z;double SUM;for(Z=0,SUM=0;Z<n;Z++){ SUM += a[Z]; } if( SUM/(double)n > 1.0+EPSILON ){ printf("Mean of Rates Error :: %f\n",SUM/(double)n); }

/*  no action  */
//#define CHECK_MEM(a);
#define CHECK_POSITIVE(a,n);
//#define CHECK_ZEROES(a,n);
#define CHECK_MEAN(a,b);

//------------------------------------------------------------------------------
/** 
 * CONVENTIONS:
 *
 * (double*) variables that start with a:
 *        ...uppercase letter --> MATRICES
 *        ...lowercase letter --> VECTORS
 *
 * Xt is the transpose of X
 * Xi is the inverse of X
 *
 * D is a diagonal matrix of eigenvalues
 * U are eigenvectors (column major, Ut is rowmajor)
 * pi are the priors as a vector
 * PI are the priors as a diagonal matrix
 *
 * (int) variable:
 *        n is the number of columns/size of alphabet
 *        m is the number of rows
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
        printf("\n"); 
    }
}
/* prints an array horizontally */
void printarray( const double* a, const int n )
{
    int i;
    for(i=0;i<n;++i)
        printf("[%6.5f] ", a[i]);
    putchar('\n');
}
/* print out a caml value */
void CAML_debug( value s )
{
    mll* a; int i; 
    a = ML_val( s );
    //printf("DEBUG:\n\tN: %d\n\tL: %d", a->stride, a->c_len );
    printf("DEBUG:");
    for (i=0;i< (a->stride*a->c_len);i++){
        if( i % a->stride == 0 )
            printf("\n\t");
        printf(" [ %f ] ", a->lv_s[i]);
    }
    printf("\n\n");
}
/** creates a random substituation rate matrix
 * each row sums to 0 on diagonal elements
 */
void rand_sub_mat_gtr( double* A, const int n)
{
    srand(time(NULL) + getpid());
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
    srand(time(NULL) + getpid());
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
    int i,j=s->stride*s->c_len;
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
//------------------------------------------------------------------------------

/**  [mk_diag diag mat n m] ~ makes a diagonal matrix from an array
 * Places the [m] elements of [diag] along diagonal of [mat] with [n] cols
 */
#ifdef _WIN32
__inline void 
#else
inline void 
#endif
mk_diag(const double* diag, double* M, const int n, const int m)
{
    int i;
    for(i =0; i < m; ++i)
        M[i*n+i] = diag[i];
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
 */
void scale_VL( const double *VR, double *VL, const int n)
{
    double *alphas,a = 1,b = 0;
    int i,j; 
    char ntran = 'N', _tran = 'T';
    alphas = (double*) malloc(n*n*sizeof(double));
    CHECK_MEM( alphas );

    dgemm_( &_tran, &ntran, &n, &n, &n,
            &a, VL, &n, VR, &n, &b, alphas, &n );

    for(i = 0;i < n; i++){
        for(j=0;j<n;j++)
            VL[i*n+j] = VL[i*n+j] / alphas[i*n+i];
    }
    free(alphas);
}

/** [mk_inverse VL D VR] makes the inverse of VL from the diagonalization
 *
 * Since lapack doesn't guarentee that VL = inv(VR) when we have duplicate
 * eigenvalues (this happens less often in GTR, but symmetric matrices have
 * an easy solution, Ut = Ui). Otherwise we just scale Ui's rows (scale_VL).
 *
 * output is in VL.
 */
int
mk_inverse( double *VL, const double *D, const double *VR, const int n)
{
    int i;
    //LU factorization and inverse via DGETRF and DGETRI
    memcpy(VL, VR, n*n*sizeof(double) );
    int *piv,lwork = -1; 
    double work_size, *work;
    piv = malloc(n*sizeof(double));
    dgetrf_(&n, &n, VL, &n, piv, &i);
    if( !i ){
        dgetri_(&n, VL, &n, piv, &work_size, &lwork, &i); //optimal work
        if( !i ){
            lwork = (int)work_size;
            work = (double*) calloc( lwork , sizeof( double) );
            dgetri_(&n, VL, &n, piv, work, &lwork, &i);
            free( work );
        }
    }

    free( piv );
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
    int i,s = m * n;
    for(i=0;i<s;i=i+n+1 )
        D[i] = exp( D[i] * t );
}

//------------------------------------------------------------------------------
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

    ret = (mll*) malloc( sizeof(mll));
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

/**
 * compares two character sets
 */
int compare_chars( const mll* c1, const mll* c2)
{
    int ret,i,n;
    //first check sizes are the same
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
{
    mll *ccs1,*ccs2;
    int i;
    ccs1 = ML_val(c1);
    ccs2 = ML_val(c2);
    return compare_chars( ccs1, ccs2 );
}

//------------------------------------------------------------
/** [gamma z]
 * Computes the gamma function of z using Lanczos Approximation. (e~1e-15)
 */
double gamma( double z )
{
    double x,Lg; int k; int g = 7;
    static double C[] = 
    {  0.99999999999980993, 676.5203681218851, -1259.1392167224028,
        771.32342877765313, -176.61502916214059, 12.507343278686905,
        -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7 };

    if( z < .5 )
        return M_PI / (sin(M_PI*z)*gamma(1-z));
    z--;
    Lg = C[0];
    for(k=1;k<g+2;++k)
        Lg = Lg + (C[k] / (z+k));

    x = z + g + 0.5;
    return sqrt(2*M_PI) * pow(x, z+0.5) * exp(-x) * Lg;
}

/** [gamma_pdf r alpha beta]
 * Return the probability density function of the gamma distribution
 */
double gamma_pdf(const double r, const double alpha, const double beta)
{
    return (pow(beta,alpha)*pow(r, alpha-1))/(exp(beta*r)*gamma(alpha));
}

/** [gamma_pp out k alpha beta step]
 * Finds the cut points, into [out], that create fractions of the percent total
 */
void gamma_pp(double* out,const int k,const double alpha,const double beta)
{
    double x,y1,y2,y,z,f; int i;
    x = STEP; y = 0; z = 0; f = STEP/2; //initial values; gamma_pdf(0) = inf

    for(i=0;i<(k-1);i++){
        z = ((double)i+1) / ((double)k);
        y1 = gamma_pdf(x,alpha,beta);
        while( y <= z ){
            x += STEP;
            y2 = gamma_pdf(x,alpha,beta);
            y += f*(y1+y2);
            y1 = y2;
        }
        out[i] = x;
        z += 1/k;
    }
}

/** confluent hypergeometric (for incomplete gamma ratio)
 *   ___inf
 *   \       ______x^r_______ 
 *   /      (p+1)(p+2)...(p+r)
 *   ---r=1
 */
double gamma_M(const double z, const double a) 
{
    double bottom,sum,z_pow,prev; int iter;
    iter = 0; bottom = 1; z_pow = 1; sum = 0; prev = 0;

    do {
        prev = sum;
        sum += z_pow / bottom;
        iter++;
        bottom = bottom * (a+iter);
        z_pow = z_pow * z;
    } while( fabs(sum-prev) > EPSILON );

    return sum;
}

/** [gamma_i x a ]
 * calculates the incomplete gamma ratio, based on methods described in
 *      ~Bhattacharjee (1970), Algorithm AS32
 *
 *      a - alpha parameter
 *      x - bound of integral
 */
double gamma_i(const double x, const double a)
{
    return (exp(-x)*pow(x,a)/gamma(a+1)) * gamma_M(x,a);
}

/** [gamma_rates rates alpha beta cuts k]
 * calculates the rates into [rates] with percentage points [cuts] of [k]
 * categories and shape parameters [alpha] and [beta]
 */ 
void 
gamma_rates(double* rates,const double a,const double b,const double* cuts,const int k)
{
    double fac, *ingam;
    fac = a*((double)k)/b;
    ingam = (double*) malloc( k * sizeof(double));
    CHECK_MEM(ingam);
    int j;

    //calculate: rj = (A*k/B)*(I(bB,A+1)-I(aB,A+1))
    for(j=0;j<(k-1);j++)
        ingam[j] = gamma_i( cuts[j]*b, a+1 );
    rates[0]   = ingam[0] * fac;            //lower rate
    rates[k-1] = (1 - ingam[k-2]) * fac;    //upper rate
    for(j=1;j<k-1;j++)
        rates[j] = (ingam[j] - ingam[j-1]) * fac;
    CHECK_MEAN( rates, k );
    free( ingam );
}

value gamma_CAML_rates( value a, value b, value c )
{
    CAMLparam3(a,b,c); 
    CAMLlocal1( rates );
    double alpha,beta,*rate_ray,*pcut_ray;
    int cats,j;

    alpha = Double_val( a );
    beta = Double_val( b );
    cats = Int_val( c );
    rate_ray = (double*) malloc( sizeof(double)*c ); 
    CHECK_MEM(rate_ray);
    pcut_ray = (double*) malloc( sizeof(double)*c );
    CHECK_MEM(pcut_ray);

    gamma_pp( pcut_ray, cats, alpha, beta );
    gamma_rates( rate_ray, alpha, beta, pcut_ray, cats );

    free( pcut_ray );
    long dims[1]; dims[0] = cats;
    rates = alloc_bigarray(BIGARRAY_FLOAT64 | BIGARRAY_C_LAYOUT,1,rate_ray,dims);

    CAMLreturn( rates );
}
//------------------------------------------------------------------------------

/**  [diagonalize_sym A D N]
 * Diagonalizes [A] that is [N]x[N] and upper symmetric and places
 * the resultant eigenvalues in the matrix [D] along the diagonal.
 * Overwrites [A] with the eigenvectors (rowmajor).
 *
 * @returns =0 on success
 *          <0 ith argument had an issue
 *          >0 failed in convergence
 */
int diagonalize_sym(double* A, double* D, const int n)
{
    char jobz = 'V';
    char uplo = 'U';
    int info = 0,lwork=-1;
    double *eigen,*work,work_size;
    eigen = (double*) malloc(n*sizeof(double));
    CHECK_MEM(eigen);

    //find the optimal work (call func with -1 in lwork)
    //result ends up in work_size
    dsyev_(&jobz, &uplo, &n, A, &n, eigen, &work_size, &lwork, &info);
    if( info == 0 ){
        lwork = (int)work_size;
        work = (double*) calloc( lwork , sizeof( double) );
        CHECK_MEM(work);

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
        dsyev_(&jobz, &uplo, &n, A, &n, eigen, work, &lwork, &info);
        if( info == 0 ){
            mk_diag(eigen, D, n,n);
        }
        free( work );
    }
    free( eigen );
    return info;
}

int diagonalize_gtr(double* A, double* D, double* Ui, const int n)
{
    char jobv_ = 'V'; //we went to find left and right eigenvectors
    double *wi,*wr,*U,*work,work_size;
    int lwork,info;
    wi = (double*) malloc( n*sizeof(double));
    CHECK_MEM(wi);
    wr = (double*) malloc( n*sizeof(double));
    CHECK_MEM(wr);
    U  = (double*) calloc( n*n,sizeof(double));
    CHECK_MEM(U);

    //find the optimal work (call func with -1 in lwork)
    lwork = -1;
    dgeev_(&jobv_,&jobv_,&n,A,&n,wr,wi,U,&n,Ui,&n,&work_size,&lwork,&info);
    if( info == 0 ) {

        lwork = (int)work_size;
        work = (double*) calloc( lwork , sizeof( double) );
        CHECK_MEM(work);
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
         * WR       - real parts of eigenvalues (**OUTPUT**)
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
        dgeev_(&jobv_,&jobv_,&n,A,&n,wr,wi,U,&n,Ui,&n,work,&lwork,&info);
        CHECK_ZEROES(wi,n);
        if( info == 0 ) {
            //imaginary eigenvals are ignored --since the Matrix is similar to 
            //some symmetric matrix, thus have same eigenvalues (Keilson 1979)
            mk_diag(wr,D, n,n);
            mk_inverse(U, D, Ui, n);
        }
        free( work );
    }
    memcpy(A,U,n*n*sizeof(double));
    free(U); free(wr); free(wi);

    return info;
}

/** [mk_probmat_*** P U D [Ui] t]
 * Finds the probability matrix based on the diagonalized substitution
 * matrix, [U] and [D], with branch length [t]. returns result in [P].
 *
 * UNCHECKED EXCEPTION ::
 * [n]x[n] == dim([D]) == dim([P]) == dim([U]) == dim([Ui])
 */  
void
compose_sym(double* P,const double* U,const double* D,const float t,const int n,double *tmp)
{
    char _tran = 'T'; double alpha = 1;
    char ntran = 'N'; double beta = 0;
    char ctran = 'C';

    memcpy(P, D, n*n*sizeof(double) );
    apply_exp(P,n,n,t); //exp(D*t); along diagonal only
    //calculates: P = op(A)*op(B)*a + C*b
    dgemm_(&ntran,&ntran,       //format, op(A), op(B)
            &n, &n, &n, &alpha,  //rows A, cols B, cols A, multiplier of (A*B)
            U, &n,              //MATRIX A, stride for A
            P, &n,              //MATRIX B, stride for B
            &beta, tmp, &n );     //multiplier for C, MATRIX C, stride for C
    //if scalor mult of C == 0, C isn't used to add
    dgemm_(&ntran,&_tran,&n,&n,&n,&alpha,tmp,&n,U,&n,&beta,P,&n);
}

void
compose_gtr(double* P, const double* U, const double* D, const double* Ui, 
        const double t, const int n,double *tmp)
{
    //constants for computation
    char _tran = 'T'; double alpha = 1;
    char ntran = 'N'; double beta = 0;
    char ctran = 'C';

    memcpy(P,D,n*n*sizeof(double));
    apply_exp(P,n,n,t);
    //S is U*exp(D)...
    dgemm_(&ntran,&ntran,&n,&n,&n,&alpha,Ui,&n,P,&n,&beta,tmp,&n);
    //P becomes U*expD*Ui... done --note: VL = inv(VR)
    dgemm_(&ntran,&ntran,&n,&n,&n,&alpha,tmp,&n,U,&n,&beta,P,&n);
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
        //for col element in Pi and l
        for(j=0;j<a;++j)
            elm = elm + (P[(i*a)+j] * l[c+j]);
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
        //printf("\t%f  * %f ---> %f ---> %f\n",l->lv_s[i],p[j],l->lv_s[i]*p[j],tmp); 
    }
    //putchar('\n');
    ret -= log (tmp);
    return ret;
}

/** [median_c Pa Pb a b c]
 * Finds the likelihood vector [c] based on vectors of its children [a] and
 * [b] with the probability matrices [Pa] and [Pb], respectively, and
 * probability [p].
 */
void
median_charset(const double* Pa,const double* Pb,
        const mll* a,const mll* b, mll* c,const double p)
{
    assert( a->stride == b->stride );
    assert( b->c_len == a->c_len );

    int i,j,len;
    double *tmp1,*tmp2;
    tmp1 = (double*) malloc( a->stride * sizeof(double) );
    CHECK_MEM(tmp1);
    tmp2 = (double*) malloc( a->stride * sizeof(double) );
    CHECK_MEM(tmp2);

    for(i=0; i < a->c_len; ++i){
        len = i*a->stride; //find half medians for each character
        median_h( Pa, a->lv_s, len, tmp1, a->stride );
        median_h( Pb, b->lv_s, len, tmp2, a->stride );
        //pairwise multiplication
        for(j=0;j < a->stride;++j){
            c->lv_s[len+j] += (tmp2[j]*tmp1[j])*p;
            //printf("Character %d, Index %d : %f * %f\n",i,j,tmp2[j],tmp1[j]);
        }
    }
    c->c_len = a->c_len;
    c->stride = a->stride;
    CHECK_POSITIVE( c->lv_s, c->stride*c->c_len );
    free( tmp1 ); free( tmp2 );
}
//------------------------------------------------------------------------------
void SWITCH(ptr *a, ptr *b)
{
    ptr *c;
    //printf("B: %f <--> %f \n",b->ll,a->ll);
    c = (ptr*) malloc( sizeof(ptr));
    c->ta = a->ta;  c->tb = a->tb;  c->ll = a->ll;  c->vs = a->vs;
    a->ta = b->ta;  a->tb = b->tb;  a->ll = b->ll;  b->vs = c->vs;
    b->ta = c->ta;  b->tb = c->tb;  b->ll = c->ll;  b->vs = c->vs;
    //printf("A: %f <--> %f \n",b->ll,a->ll);
    free( c );
}
void simplex_sym(ptr *simp,double *PA,double *PB,const double *U,const double *D,const mll *a,
    const mll *b,const double r0,const double r1,const double *gam,const double *prob,
    const double *pi,const int g_n,double *tmp)
{
    int i =0;
    mll *c; c = simp->vs;
    simp->ta = r0;
    simp->tb = r1;
    memset(c->lv_s,0,a->c_len*a->stride*sizeof(double));
    for(;i<g_n;i++){
        compose_sym( PA, U, D, r0*gam[i], a->stride,tmp );
        //if( cta != ctb ) PB = PA; else
        compose_sym( PB, U, D, r1*gam[i], b->stride,tmp );
        median_charset( PA, PB, a,b,c,prob[i] ); //adds median to c with rate tmp
    }
    simp->ll = loglikelihood( c,pi );
    //printf( "updating %d: score %f with %f and %f\n",c->id, simp->ll, simp->ta, simp->tb);
}
void simplex_gtr(ptr *simp,double *PA,double *PB,const double *U,const double *D,const double *Ui,
    const mll *a,const mll *b,const double r0,const double r1,const double *gam,const double *prob,
    const double *pi,const int g_n,double *tmp)
{
    int i = 0;
    mll *c; c = simp->vs;
    simp->ta = r0;
    simp->tb = r1;
    memset(c->lv_s,0,a->c_len*a->stride*sizeof(double));
    for(;i<g_n;i++){
        compose_sym( PA, U, D, r0*gam[i], a->stride,tmp );
        //if( cta != ctb ) PB = PA; else
        compose_sym( PB, U, D, r1*gam[i], b->stride,tmp );
        median_charset( PA, PB, a,b,c,prob[i] ); //adds median to c with rate tmp
    }
    simp->ll = loglikelihood( c,pi );
    //printf( "updating %d: score %f with %f and %f\n",c->id, simp->ll, simp->ta, simp->tb);
}
void
readjust_simplex_sym(const double* U,const double* D,const mll* a,const mll* b,
        mll* c,double* b_ta,double* b_tb,double* b_mle,const double* rates,
        const double *prob, const int g_n,const double* pi){}

void
readjust_simplex_gtr(const double* u,const double* d,const double* ui,const mll* a,
        const mll* b, mll* c,double* b_ta,double* b_tb,double* b_mle,
        const double *rates,const double *prob,const int g_n,const double* pi){}

//-----------------------------------------------------------------------------
/* Diagonlize matrix interfaces: symmetric */
value
likelihood_CAML_diagonalize_sym(value Q, value D)
{
    CAMLparam2( Q, D );
    int n = Bigarray_val( Q )->dim[0];
    double* c_Q = (double*) Data_bigarray_val( Q );
    double* c_D = (double*) Data_bigarray_val( D );

    diagonalize_sym( c_Q, c_D, n );
    CAMLreturn0;
}

/* Diagonalize matrix interface: general */
value
likelihood_CAML_diagonalize_gtr(value Q, value D, value Qi)
{
    CAMLparam3( Q, D, Qi );
    int n = Bigarray_val( Q )->dim[0];
    double* c_Q = (double*) Data_bigarray_val( Q );
    double* c_D = (double*) Data_bigarray_val( D );
    double* c_Qi= (double*) Data_bigarray_val( Qi);

    diagonalize_gtr( c_Q, c_D, c_Qi, n );
    CAMLreturn0;
}
value likelihood_CAML_median_wrapped_sym
    (value U,value D,value ta,value tb,value ml_a,value ml_b,value rates,value probs)
{
    CAMLparam5(U,D,ta,tb, ml_a);
    CAMLxparam3( ml_b,rates,probs );
    CAMLlocal1( ml_c );
    double cta,ctb,*c_U,*c_D,*c_lv,*PA,*PB,*tmp,*g_rs,*p_rs;
    mll *a,*b,*c;
    int num_rates,num_probs,i=0;

    num_rates = Bigarray_val(rates)->dim[0];
    num_probs = Bigarray_val(probs)->dim[0];
    assert( num_rates == num_probs );

    c_U = (double*) Data_bigarray_val( U );
    c_D = (double*) Data_bigarray_val( D );
    g_rs = (double*) Data_bigarray_val( rates );
    p_rs = (double*) Data_bigarray_val( probs );
    cta = Double_val( ta );
    ctb = Double_val( tb );
//    printf("TA:%f\nTB:%f\n",cta,ctb);
    assert( cta > 0 && ctb > 0 );

    a = ML_val( ml_a );
    b = ML_val( ml_b );
    assert( a->stride == b->stride );
    assert( a->c_len == b->c_len );

    ml_c = caml_alloc_custom(&likelihood_custom_operations, (sizeof(mll*)), 
                                                    a->stride*a->c_len, CAML_ALLOC_2);
    c = (mll*) malloc( sizeof(mll) );
    CHECK_MEM(c);
    ML_val( ml_c ) = c;

    c->stride = a->stride;
    c->c_len = a->c_len;

    assert( a->c_len == b->c_len );
    assert( b->stride == a->stride );

    c->lv_s = (double*) calloc( c->c_len * c->stride, sizeof(double));
    CHECK_MEM(c->lv_s);
    PA = (double*) calloc( c->stride * c->stride, sizeof(double));
    CHECK_MEM(PA);
    PB = (double*) calloc( c->stride * c->stride, sizeof(double));
    CHECK_MEM(PB);

    //below, used for less malloc calls in intermediary multiplications
    tmp = (double*) calloc( a->stride * b->stride, sizeof(double));
    CHECK_MEM(tmp);

    /** 
     *   ___ num_rates
     *   \     
     *   /__ pr_i * F(x | r_i,P )
     *    i
     *        where:  F(x|r,P) = f_right(x*r|P) * f_left(x*r|P)
     *              f_j(x*r|P) = SUM(P_jk,ke{ACTG})
     *
     *              P    = probability matrix
     *              pr_i = the probability of the rate
     *              r_i  = the rate of i
     *              x    = the branch length
     *              P_j  = is the row jth row of P
     */
    for(;i<num_rates;i++){
        compose_sym( PA, c_U, c_D, cta*g_rs[i], a->stride,tmp );
        //if( cta != ctb ) PB = PA; else
        compose_sym( PB, c_U, c_D, ctb*g_rs[i], b->stride,tmp );
        median_charset( PA, PB, a,b,c,p_rs[i] ); //adds median to c with rate tmp
    }

    free( PA ); free( PB ); free( tmp );
    assert( c == ML_val(ml_c));
    CAMLreturn(ml_c);
}

value likelihood_CAML_median_wrapped_gtr
    (value U,value D,value Ui,value ta,value tb,value ml_a,value ml_b,value rates,value probs)
{
    CAMLparam5(U,D,Ui,ta,tb);
    CAMLxparam4( ml_a,ml_b,rates,probs );
    CAMLlocal1( ml_c );
    double cta,ctb,*c_U,*c_D,*c_Ui,*c_lv,*PA,*PB,*tmp,*g_rs,*p_rs;
    mll *a,*b,*c;
    int num_rates,num_probs,i = 0;

    num_rates = Bigarray_val(rates)->dim[0];
    num_probs = Bigarray_val(probs)->dim[0];
    assert( num_rates == num_probs );

    c_U = (double*) Data_bigarray_val( U );
    c_D = (double*) Data_bigarray_val( D );
    c_Ui= (double*) Data_bigarray_val( Ui);
    g_rs= (double*) Data_bigarray_val(rates);
    p_rs= (double*) Data_bigarray_val(probs);
    cta = Double_val( ta );
    ctb = Double_val( tb );
    assert( cta > 0 && ctb > 0 );

    a = ML_val( ml_a );
    b = ML_val( ml_b );
    assert(a->stride == b->stride);
    assert(b->c_len == a->c_len);

    ml_c = caml_alloc_custom(&likelihood_custom_operations, (sizeof(mll*)), 
                                                    a->stride*a->c_len, CAML_ALLOC_2);
    c = (mll*) malloc( sizeof(mll) );
    CHECK_MEM(c);
    ML_val( ml_c ) = c;

    c->lv_s = calloc(a->c_len * a->stride, sizeof(double)); 
    PA = (double*) calloc( a->stride * a->stride, sizeof(double));
    CHECK_MEM(PA);
    PB = (double*) calloc( b->stride * b->stride, sizeof(double));
    CHECK_MEM(PB);
    tmp = (double*) calloc( b->stride * b->stride, sizeof(double));
    CHECK_MEM(tmp);

    for(;i<num_rates;i++){
        compose_gtr( PA, c_U, c_D, c_Ui, cta*g_rs[i], a->stride, tmp);
        //if( cta != ctb ) PB = PA; else
        compose_gtr( PB, c_U, c_D, c_Ui, ctb*g_rs[i], b->stride, tmp);
        median_charset( PA, PB, a,b,c, p_rs[i] );
    }

    free( PA ); free( PB ); free( tmp );
    assert( c == ML_val(ml_c) );
    CAMLreturn(ml_c);
}

/* [likelihood_CAML_median_sym ,,,] argument wrapper for median_sym */
value likelihood_CAML_median_sym(value * argv, int argn)
{
    return likelihood_CAML_median_wrapped_sym
        ( argv[0],argv[1],argv[2],argv[3],argv[4],argv[5],argv[6],argv[7] ); 
}

/* [likelihood_CAML_median_gtr ...] argument wrapper for median_gtr */
value likelihood_CAML_median_gtr(value * argv, int argn)
{
    return likelihood_CAML_median_wrapped_gtr
        ( argv[0],argv[1],argv[2],argv[3],argv[4],argv[5],argv[6],argv[7],argv[8] );
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

/* [likelihood_CAML_compose_sym u d t] testing function. composes a probability 
 * matrix with U D and t */
value likelihood_CAML_compose_sym(value U, value D, value t)
{
    CAMLparam3( U,D,t );
    CAMLlocal1( res );
    double *c_P, c_t,*c_D,*c_U,*c_T;
    int n;

    n = Bigarray_val( U ) -> dim[0];
    c_t = Double_val( t );
    c_U = (double *) Data_bigarray_val( U );
    c_D = (double *) Data_bigarray_val( D );
    c_P = (double*) malloc( n*n*sizeof(double));
    CHECK_MEM(c_P);
    c_T = (double*) malloc( n*n*sizeof(double));

    compose_sym(c_P,c_U,c_D,c_t,n,c_T);
    long dims[2];
    dims[0] = n; dims[1] = n;
    res = alloc_bigarray(BIGARRAY_FLOAT64 | BIGARRAY_C_LAYOUT, 2, c_P, dims);
    CAMLreturn ( res );
}
value 
likelihood_CAML_compose_gtr(value U, value D, value Ui, value t)
{
    CAMLparam4( U,D,Ui,t );
    CAMLlocal1( res );
    double *c_P, c_t,*c_D,*c_U,*c_Ui,*c_T;
    int n;

    n = Bigarray_val( U ) -> dim[0];
    c_t = Double_val( t );
    c_U = (double*) Data_bigarray_val( U );
    c_D = (double*) Data_bigarray_val( D );
    c_Ui= (double*) Data_bigarray_val( Ui);
    c_P = (double*) malloc( n*n*sizeof(double));
    CHECK_MEM(c_P);
    c_T = (double*) malloc( n*n*sizeof(double));
    CHECK_MEM(c_T);

    compose_gtr(c_P,c_U,c_D,c_Ui,c_t,n,c_T);
    long dims[2];
    dims[0] = n; dims[1] = n;
    res = alloc_bigarray(BIGARRAY_FLOAT64 | BIGARRAY_C_LAYOUT, 2, c_P, dims);
    CAMLreturn( res );
}

/** [..._readjust_gtr... U D Ui A B C ta tb gammas probs pis ll] **/
value
likelihood_CAML_readjust_gtr_wrapped
    (value U,value D,value Ui,value A,value B,value C,value ta,value tb,value g,value p,value pis,value ll)
{   

    double cta,ctb,*rates,*pi,likelihood,*c_U,*c_D,*c_Ui,ml,*probs;
    mll *a,*b,*c;
    int g_n, pi_n,p_n;

    //macros for GC
    CAMLparam5(U,D,Ui,A,B);
    CAMLxparam5(C,ta,tb,g,p);
    CAMLxparam2(ll,pis);
    CAMLlocal1( res );

    //macros to 'convert' values
    a = ML_val( A );
    b = ML_val( B );
    c = ML_val( C );
    cta = Double_val( ta );
    ctb = Double_val( tb );

    c_U = (double*) Data_bigarray_val( U );
    c_D = (double*) Data_bigarray_val( D );
    c_Ui= (double*) Data_bigarray_val( Ui);

    pi_n = Bigarray_val(pis) ->dim[0];
    pi = (double*) Data_bigarray_val(pis);

    g_n = Bigarray_val( g ) ->dim[0];
    rates = (double*) Data_bigarray_val(g);
    p_n = Bigarray_val( p ) ->dim[0];
    probs = (double*) Data_bigarray_val(p);
    likelihood = Double_val( ll );

    //printf("%f\t%f\t%f\n", cta,ctb,likelihood);
    readjust_simplex_gtr(c_U,c_D,c_Ui,a,b,c,&cta,&ctb,&likelihood,rates,probs,g_n,pi);
    //printf("%f\t%f\t%f\n", cta,ctb,likelihood);
    
    res = caml_alloc_tuple( 3 );
    Store_field(res, 0, caml_copy_double(cta));
    Store_field(res, 1, caml_copy_double(ctb));
    Store_field(res, 2, caml_copy_double(likelihood));
    CAMLreturn(res);
}
value
likelihood_CAML_readjust_sym_wrapped
    (value U,value D,value A,value B,value C,value ta,value tb,value g,value p,value pis,value ll)
{
    //macros for GC
    CAMLparam5(U,D,ll,A,B);
    CAMLxparam5(C,ta,tb,g,p);
    CAMLxparam1( pis );
    CAMLlocal1( res );
    
    double cta,ctb,*rates,*probs,*pi,likelihood,*c_U,*c_D,ml;
    mll *a,*b,*c;
    int g_n, p_n, pi_n;

    //macros to 'convert' values
    a = ML_val( A );
    b = ML_val( B );
    c = ML_val( C );
    cta = Double_val( ta );
    ctb = Double_val( tb );
    c_U = (double*) Data_bigarray_val( U );
    c_D = (double*) Data_bigarray_val( D );
    pi_n = Bigarray_val( pis )->dim[0];
    pi = (double*)Data_bigarray_val( pis );
    g_n = Bigarray_val( g )->dim[0];
    rates = (double*) Data_bigarray_val(g);
    p_n = Bigarray_val( p )->dim[0];
    probs = (double*) Data_bigarray_val(p);
    likelihood = Double_val( ll );

    assert( p_n == g_n );
    assert( a->stride == b->stride && a->stride == c->stride );
    assert( a->stride == pi_n );

    //printf("%f\t%f\t%f\n", cta,ctb,likelihood);
    readjust_simplex_sym(c_U,c_D,a,b,c,&cta,&ctb,&likelihood,rates,probs,g_n,pi);
    //printf("%f\t%f\t%f\n", cta,ctb,likelihood);
    
    res = caml_alloc_tuple( 3 );
    Store_field(res, 0, caml_copy_double (cta));
    Store_field(res, 1, caml_copy_double (ctb));
    Store_field(res, 2, caml_copy_double (likelihood));
    CAMLreturn(res);
}


/** functions for interface **/
value
likelihood_CAML_readjust_gtr(value * argv, int argn)
{
    return likelihood_CAML_readjust_gtr_wrapped
        (argv[0], argv[1],argv[2],argv[3],argv[4],argv[5],argv[6],argv[7],argv[8],argv[9],argv[10],argv[11]);
}
value
likelihood_CAML_readjust_sym(value * argv, int argn)
{
    return likelihood_CAML_readjust_sym_wrapped
        (argv[0], argv[1],argv[2],argv[3],argv[4],argv[5],argv[6],argv[7],argv[8],argv[9],argv[10]);
}

#endif /* USE_LIKELIHOOD */
