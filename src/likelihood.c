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
#define EPSILON 1e-5 //error for numerical calculations
#define ML_val(v) (*((struct ml**)Data_custom_val(v)))
typedef struct ml mll; 

/* assert checking 
 *
 * ~CHECK_MEM       -- verify enough memory
 * ~CHECK_POSITIVE  -- verify array is all positive values with EPSILON error
 * ~CHECK_ZEROES    -- verify array is all zeroes with EPSILON error
 */

/*  asserts  */
#define CHECK_MEM(a) if(a==NULL) failwith("I can't allocate more memory.")
//#define CHECK_POSITIVE(a,n); int Z;for(Z=0;Z<n;Z++){ if( a[Z] < -EPSILON) failwith("Negative Likelihood"); }
//#define CHECK_ZEROES(a,n); int Z;for(Z=0;Z<n;Z++){ if(a[Z] > EPSILON || a[Z] < -EPSILON) failwith("Imaginary eigenvalue"); }
//#define CHECK_MULTIPLCITY(a); if( a > (1.0+EPSILON) ){ failwith("Eigensystem doesn't satisfy eigenvalues"); }

/*  printline  */
//#define CHECK_MEM(a) if(a==NULL) printf("I can't allocate more memory. %d",__LINE__)
#define CHECK_POSITIVE(a,n); int Z;for(Z=0;Z<n;Z++){ if( a[Z] < -EPSILON){ printf("Negative Likelihood :: %f\n",a[Z]);} }
#define CHECK_ZEROES(a,n); int Z;for(Z=0;Z<n;Z++){ if(a[Z] > EPSILON || a[Z] < -EPSILON){ printf("Imaginary eigenvalue :: %f\n", a[Z]); } }
//#define CHECK_MULTIPLCITY(a); if( a > (1.0+EPSILON) ){ printf("Algebraic Multiplicity > 1"); }

/*  no action  */
//#define CHECK_MEM(a);
//#define CHECK_POSITIVE(a,n);
//#define CHECK_ZEROES(a,n);
#define CHECK_MULTIPLCITY(a);


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

//---------------------------------------------------------------------------------------

/* prints a matrix (block format) */
void
printmatrix( const double* Z, const int n, const int m)
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
void
printarray( const double* a, const int n )
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
void
rand_sub_mat_gtr( double* A, const int n)
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
void
rand_sub_mat_sym( double* A, const int n)
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
/* release struct ml->lv_s array to pool */
void likelihood_CAML_free( value v )
{
    mll* val;
    val = ML_val(v);
    //free (val->lv_s); //causing double free?!?
    return;
}
/* serialization of struct ml */
void likelihood_CAML_serialize(value v, unsigned long* wsize_32, unsigned long* wsize_64)
{
    CAMLparam1(v);
    mll* s;
    s = ML_val(v);
    caml_serialize_int_4( s->stride );
    caml_serialize_int_4( s->c_len );
    int i,j=s->stride*s->c_len;
    for(i=0;i<j;++i)
        caml_serialize_float_8( s->lv_s[i] );
    CAMLreturn0;
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
void
scale_VL( const double *VR, double *VL, const int n)
{
    double *alphas,a = 1,b = 0;
    int i,j; 
    char ntran = 'N', _tran = 'T';
    alphas = (double*) malloc(n*n*sizeof(double));
    CHECK_MEM( alphas );

    dgemm_( &_tran, &ntran, &n, &n, &n,
                &a, VL, &n, VR, &n, &b, alphas, &n );

    for(i = 0;i < n; i++){
        CHECK_MULTIPLCITY( alphas[i*n+i] );
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
   /* 
    //check for algebriac multiplicity
    int i,j,k = 1;
    for(i=0;i<n-1;++i){
        for(j=i+1;j<n;++j){
            if( VL[i*n+i] == VL[j*n+j] ){
                k = 0;
                break;
    }   }   }
    if( k == 1 ){
        scale_VL(VR, VL, n);
        transpose( VR );
        return 0;
    } */ int i;

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
    int s = m * n;
   int i = 0;
    for(; i<s; i=i+n+1 )
        D[i] = exp( D[i] * t );
}

/**
 *  convert struct ml to a bigarray.array2
 */
value StoBigarray( const mll* charss )
{
    CAMLlocal1( res );
    long dims[2];
    dims[0] = charss->c_len;
    dims[1] = charss->stride;
    res = alloc_bigarray( BIGARRAY_FLOAT64 | BIGARRAY_C_LAYOUT, 2, charss->lv_s, dims );
    return (res);
}
value likelihood_CAML_StoBigarray( value s )
{
    CAMLparam1( s );
    CAMLlocal1( res );
    mll* work;
    long dims[2];

    work = ML_val( s );
    dims[0] = work->c_len;
    dims[1] = work->stride;
    res = alloc_bigarray( BIGARRAY_FLOAT64 | BIGARRAY_C_LAYOUT, 2, work->lv_s, dims );
    CAMLreturn( res );
}

/**
 *  convert bigarray.array2 to struct ml
 */
value likelihood_CAML_BigarraytoS( value A )
{
    CAMLparam1( A );
    CAMLlocal1( s );
    double* stuff;
    mll* ret;

    s = caml_alloc_custom(&likelihood_custom_operations, (sizeof(mll*)), 0, 1); 
    ret = (mll*) malloc( sizeof(mll));
    CHECK_MEM(ret);
    ML_val(s) = ret;
    stuff = (double*) Data_bigarray_val( A );

    ret->c_len = Bigarray_val(A)->dim[0];
    ret->stride = Bigarray_val(A)->dim[1];
    ret->lv_s = stuff;

    CAMLreturn( s ); 
}

/**
 * filters an ml struct with all the indexes in the passed array --should be in order
 */
value
likelihood_CAML_filter(value as, value ibs)
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

    package = caml_alloc_custom(&likelihood_custom_operations, (sizeof(mll*)), 0, 1);
    ret = ML_val( package );
    ret->stride = norm->stride;
    ret->c_len = norm->c_len - m;
    CAMLreturn( package );
}

/**
 * compares two character sets
 */
int
compare_chars( const mll* c1, const mll* c2)
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

int  
likelihood_CAML_compare( value c1, value c2 )
{
    CAMLparam2( c1, c2 );
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
 *
 */
void
gamma_pp(double* out,const int k,const double alpha,const double beta)
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
    } while( abs(sum-prev) > EPSILON );

    return sum;
}

/** [gamma_i x a max_iter]
 * calculates the incomplete gamma ratio, based on methods described in
 *      ~Bhattacharjee (1970), Algorithm AS32
 *      ~Kostlan/Gokhman (Nov. 1987)
 *      ~Numerical Recipes in C
 */
double gamma_i(const double x, const double a)
{ 
//    if( (a <= x && x <= 1) || x < a ){ //series expansion
        return (exp(-x)*pow(x,a)/gamma(a+1)) * gamma_M(x,a);
//    } else { //continued fraction
//        return 1 - (exp(-x)*pow(x,a))/gamma(p)*gamma_Cf(x);
//    }
}

/** [gamma_rates rates alpha beta cuts k]
 * calculates the rates into [rates] with percentage points [cuts] of [k]
 * categories and shape parameters [alpha] and [beta]
 */ 
void 
gamma_rates(double* rates,const double a,const double b, const double* cuts, const int k)
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
    free( ingam );
}

value
gamma_CAML_rates( value a, value b, value c )
{
    CAMLparam3(a,b,c); //need this here?
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
    rates = alloc( cats*Double_wosize, Double_array_tag );
    for(j=0;j<cats;j++)
        Store_double_field(rates,j,rate_ray[j]);
    free( rate_ray ); 
    CAMLreturn( rates );
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
int 
diagonalize_sym(double* A, double* D, const int n)
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

int
diagonalize_gtr(double* A, double* D, double* Ui, const int n)
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
         * A        - matrix of LDAxN (**OUTPUT**)
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
compose_sym( double* P, const double* U, const double* D, const float t, const int n)
{
    char _tran = 'T'; double alpha = 1;
    char ntran = 'N'; double beta = 0;
    char ctran = 'C';

    memcpy(P, D, n*n*sizeof(double) );
    apply_exp(P,n,n,t); //exp(D*t); along diagonal only
    double *I;
    I = (double*) malloc( n*n*sizeof(double));
    CHECK_MEM(I);
                        //calculates: P = op(A)*op(B)*a + C*b
    dgemm_(&ntran,&ntran,       //format, op(A), op(B)
           &n, &n, &n, &alpha,  //rows A, cols B, cols A, scalor multiplier of (A*B)
            U, &n,              //MATRIX A, stride for A
            P, &n,              //MATRIX B, stride for B
            &beta, I, &n );     //scalor multiplier for C, MATRIX C, stride for C
                                //if scalor mult of C == 0, C isn't used to add
    dgemm_(&ntran,&_tran,&n,&n,&n,&alpha,I,&n,U,&n,&beta,P,&n);
    free( I );
}

void
compose_gtr(double* P, const double* U, const double* D, const double* Ui, 
        const double t, const int n)
{
    //constants for computation
    char _tran = 'T'; double alpha = 1;
    char ntran = 'N'; double beta = 0;
    char ctran = 'C';

    memcpy(P,D,n*n*sizeof(double));
    apply_exp(P,n,n,t);

    double *S;
    S = (double*) malloc( n*n*sizeof(double));
    CHECK_MEM(S);

    //S is U*exp(D)...
    dgemm_(&ntran,&ntran,&n,&n,&n,&alpha,Ui,&n,P,&n,&beta,S,&n);
    //P becomes U*expD*Ui... done --note: VL = inv(VR)
    dgemm_(&ntran,&ntran,&n,&n,&n,&alpha,S,&n,U,&n,&beta,P,&n);
    
    free(S);
}

/**  [median_h P l n nl nn a]
 * Finds half the likelihood vector of child [l] and probability matrix [P] into [nl]
 *
 *  [n] - starting location for vector in  [l]
 *  [nn]- starting location for vector in [nl]
 *  [a] - length of each vector (in practice, the alphabet size) [P] is [a]x[a]
 */
void 
median_h( const double* P, const double* l, double* nl, const int a)
{
    int i,j; double elm;
    //for each row of P
    for(i=0; i<a; ++i){
        elm = 0;
        //for col element in Pi and l
        for(j=0;j<a;++j)
            elm = elm + (P[i*a+j] * l[j]);
        nl[i] = elm;
    }
}

/** [loglikelihood ml p]
 * returns loglikelihood of a character set. |p| = l->stride
 */
double
loglikelihood( const mll* l, const double* p )
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
 * [b] with the probability matrices [Pa] and [Pb], respectively.
 */
void
median_charset(const double* Pa,const double* Pb,
                    const mll* a,const mll* b, mll* c,const double r)
{
    assert( a->stride == b->stride );
    assert( b->c_len == a->c_len );

    int i,j,len;
    double *tmp1,*tmp2,g_r;
    tmp1 = (double*) malloc( a->stride * sizeof(double) );
    CHECK_MEM(tmp1);
    tmp2 = (double*) malloc( a->stride * sizeof(double) );
    CHECK_MEM(tmp2);

    g_r = gamma( r );
    for(i=0; i < a->c_len ; ++i){
        len = i*a->stride;
        //find half medians for each character
        median_h( Pa, a->lv_s, tmp1, a->stride );
        median_h( Pb, b->lv_s, tmp2, a->stride );
        //pairwise multiplication
        for(j=0;j < a->stride;++j){
            c->lv_s[len+j] += tmp2[j] * tmp1[j];//* g_r;
        }
    }
    c->c_len = a->c_len;
    c->stride = a->stride;
    CHECK_POSITIVE( c->lv_s, c->stride*c->c_len );

    //ensure these are set properly
    free( tmp1 ); free(tmp2);
}

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
value
likelihood_CAML_median_wrapped_sym
    (value U,value D,value ta,value tb,value ml_a,value ml_b,value rates)
{
    CAMLparam5(U,D,ta,tb, ml_a);
    CAMLxparam2( ml_b,rates );
    CAMLlocal1( ml_c );
    double cta,ctb,*c_U,*c_D,*c_lv,*PA,*PB;
    mll *a,*b,*c;
    int num_rates,i=0;

    num_rates = Wosize_val( rates ) / 2;
    c_U = (double*) Data_bigarray_val( U );
    c_D = (double*) Data_bigarray_val( D );
    cta = Double_val( ta );
    ctb = Double_val( tb );
    a = ML_val( ml_a );
    b = ML_val( ml_b );

    ml_c = caml_alloc_custom(&likelihood_custom_operations, (sizeof(mll*)), 0, 1);
    c = (mll*) malloc( sizeof(mll) );
    CHECK_MEM(c);
    ML_val( ml_c ) = c;

    c->stride = a->stride;
    c->c_len = a->c_len;
    c->lv_s =  (double*) calloc( a->stride * a->c_len, sizeof(double));
    CHECK_MEM(c->lv_s);
    PA = (double*) calloc( a->stride * a->stride, sizeof(double));
    CHECK_MEM(PA);
    PB = (double*) calloc( a->stride * b->stride, sizeof(double));
    CHECK_MEM(PB);

    double tmp;
    for(;i<num_rates;i++){
        tmp = Double_field( rates, i );
        compose_sym( PA, c_U, c_D, cta*tmp, a->stride );
      //if( cta != ctb )
      //    PB = PA;
      //else
            compose_sym( PB, c_U, c_D, ctb*tmp, b->stride );
        median_charset( PA, PB, a,b,c,tmp ); //adds median to c with rate tmp
    }


    //divide by # of cat...
    for(i=0;i<c->c_len*c->stride;i++)
        c->lv_s[i] = c->lv_s[i] / ((double) num_rates);

    free( PA ); free( PB );
    CAMLreturn(ml_c);
}

value
likelihood_CAML_median_wrapped_gtr
    (value U,value D,value Ui,value ta,value tb,value ml_a,value ml_b,value rates)
{
    CAMLparam5(U,D,Ui,ta,tb);
    CAMLxparam3( ml_a,ml_b,rates );
    CAMLlocal1( ml_c );
    double cta,ctb,*c_U,*c_D,*c_Ui,*c_lv,*PA,*PB;
    mll *a,*b,*c;
    int num_rates,i = 0;

    num_rates = Wosize_val( rates ) / 2;
    c_U = (double*) Data_bigarray_val( U );
    c_D = (double*) Data_bigarray_val( D );
    c_Ui= (double*) Data_bigarray_val( Ui);
    cta = Double_val( ta );
    ctb = Double_val( tb );
    a = ML_val( ml_a );
    b = ML_val( ml_b );

    ml_c = caml_alloc_custom(&likelihood_custom_operations, (sizeof(mll*)), 0, 1);
    c = (mll*) malloc( sizeof(mll) );
    CHECK_MEM(c);
    ML_val( ml_c ) = c;
    c->lv_s = calloc(a->stride * a->c_len, sizeof(double)); 
    PA = (double*) calloc( a->stride * a->stride, sizeof(double));
    CHECK_MEM(PA);
    PB = (double*) calloc( b->stride * b->stride, sizeof(double));
    CHECK_MEM(PB);

    double tmp;

    for(;i<num_rates;i++){
        tmp = Double_field( rates, i );
        compose_gtr( PA, c_U, c_D, c_Ui, cta*tmp, a->stride);
      //if( cta != ctb )
      //    PB = PA;
      //else
            compose_gtr( PB, c_U, c_D, c_Ui, ctb*tmp, b->stride);
        median_charset( PA, PB, a,b,c, tmp );
    }

    //divide by # of cat...
    for(i=0;i<c->c_len*c->stride;i++)
        c->lv_s[i] = c->lv_s[i] / ((double) num_rates);

    free( PA ); free( PB );
    CAMLreturn(ml_c);
}

/* argument wrapper */
value
likelihood_CAML_median_sym(value * argv, int argn)
{
    return likelihood_CAML_median_wrapped_sym
            ( argv[0],argv[1],argv[2],argv[3],argv[4],argv[5],argv[6] ); 
}

/* argument wrapper */
value 
likelihood_CAML_median_gtr(value * argv, int argn)
{
    return likelihood_CAML_median_wrapped_gtr
            ( argv[0],argv[1],argv[2],argv[3],argv[4],argv[5],argv[6],argv[7] );
}

/* wrapper for loglikelihood */
value
likelihood_CAML_loglikelihood(value s, value p)
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

value 
likelihood_CAML_compose_sym(value U, value D, value t)
{
    CAMLparam3( U,D,t );
    CAMLlocal1( res );
    double *c_P, c_t,*c_D,*c_U;
    int n;

    n = Bigarray_val( U ) -> dim[0];
    c_t = Double_val( t );
    c_U = (double *) Data_bigarray_val( U );
    c_D = (double *) Data_bigarray_val( D );
    c_P = (double*) malloc( n*n*sizeof(double));
    CHECK_MEM(c_P);

    compose_sym(c_P,c_U,c_D,c_t,n);
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
    double *c_P, c_t,*c_D,*c_U,*c_Ui;
    int n;

    n = Bigarray_val( U ) -> dim[0];
    c_t = Double_val( t );
    c_U = (double *) Data_bigarray_val( U );
    c_D = (double *) Data_bigarray_val( D );
    c_Ui= (double *) Data_bigarray_val( Ui);
    c_P = (double*) malloc( n*n*sizeof(double));
    CHECK_MEM(c_P);

    compose_gtr(c_P,c_U,c_D,c_Ui,c_t,n);
    long dims[2];
    dims[0] = n; dims[1] = n;
    res = alloc_bigarray(BIGARRAY_FLOAT64 | BIGARRAY_C_LAYOUT, 2, c_P, dims);
    CAMLreturn( res );
}

#endif /* USE_LIKELIHOOD */
