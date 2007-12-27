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
#include <stdlib.h> //malloc, srand, RAND_MAX
#include <string.h> //memcpy, memset
#include <math.h>   //log,exp
#include <cblas.h>  //dgemm

#include <caml/mlvalues.h>
#include <caml/memory.h>
#include <caml/bigarray.h>


//---------------------------------------------------------------------------------------
/** mlchar
 * struct holding values of a character set
 */
struct c_ml
{
    int cc;     //character code
    CAMLprim value l_v;  //bigarray of likelihood vector
};
//---------------------------------------------------------------------------------------
//prints a matrix (block format)
void
printmatrix( const double* Z, const int n, const int m)
{
    int i,j;
    for (i=0; i<m; ++i) {
        putchar('\t');
        for (j=0; j<n; ++j)
            printf("[%6.5f] ", Z[i*n+j]);
        printf("\n"); 
    }
}
//prints an array horizontally
void
printarray( const double* a, const int n )
{
    int i;
    for(i=0;i<n;++i)
        printf("[%6.5f] ", a[i]);
    putchar('\n');
}
//creates a random substituation rate matrix
// *each row sums to 0 on diagonal elements*
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
//creates a random substituation rate matrix
// *each row sums to 0 on diagonal elements*
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
            tmep = ((double)rand()/(double)RAND_MAX);
            A[i*n+j] = temp; A[j*n+i] = temp;
            diag = diag + A[i*n+j];
        }
        A[i*n+i] = -diag;
    }
}

//wrapper to convert double* to bigarray with proper flags (one dim)
CAMLprim value
doublearray_to_bigarray( double* ray, int m, long* n )
{
    return alloc_bigarray( BIGARRAY_FLOAT64 | BIGARRAY_C_LAYOUT, m, ray, n );
}

//function for caml_alloc_array; used to map over 'c_ml array' -> 'mlchar array'
CAMLprim value
struct_to_value( struct c_ml * ptr )
{
    CAMLlocal1(res);
    res = caml_alloc_tuple(2);
    Store_field(res,0, Val_int( ptr->cc ) );
    Store_field(res,1, ptr->l_v );
    return res;
}

//---------------------------------------------------------------------------------------
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
 * D is a diagonal matrix of eigenvectors
 * U are eigenvectors (column major, Ut is rowmajor)
 * pi are the priors as a vector
 * PI are the priors as a diagonal matrix
 *
 * (int) variable:
 *        n is the number of columns/size of alphabet
 *        m is the number of rows
 *        t is the branch length
 */

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
    int lwork = n*n;
    double* work[lwork];
    char jobz = 'V';
    char uplo = 'U';
    int info = 0;
    double* eigen = (double*) malloc(n*sizeof(double));

    /** dsyev - calculate eigenvalues and vectors
     *          possibly use PDSYGVX or PDSYEVX --less precision)
     *
     * jobz   - V to return eigenvectors (N otherwise)
     * uplo   - Upper/Lower triangular? (U/L)
     * n      - order of the matrix
     * E_vecs - (in/out) dim(lda,n) --note: lda==n
     * n      - lda
     * e_vals - (out) eigenvalues
     * work   - (workspace/out) = dim(n)
     * lwork  - length of work >= (NB+2)*N
     * info   - retval; see comments above
     */
    dsyev_(&jobz, &uplo, &n, A, &n, eigen, &work, &lwork, &info);
    if( info != 0 )
        return info;
    mk_diag(eigen, D, n,n);

    free( eigen );
    return info;
}
int
diagonalize_gtr(double* A, double* D, double* Ui, const int n)
{
    char jobv_ = 'V'; //we went to find left and right eigenvectors
    double* wi = (double*) malloc( n*sizeof(double));
    double* wr = (double*) malloc( n*sizeof(double));

    int info, lwork = n*n;
    double* work = (double*) malloc( lwork*sizeof(double));
    double* U  = (double*) malloc( n*n*sizeof(double));

    /** dgeev   - A * v(j) = lambda(j) * v(j)
     *            u(j)**H * A = lambda(j) * u(j)**H
     *            where:
     *             ` **H is conjugate transpose)
     *             ` u(j) is a left eigenvector
     *             ` v(j) is a right eigenvector
     *             ` lambda(j) is an eigenvalue
     *
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
     *            >0 if QR failed; i+1 is firt eigenvalue
     */
    dgeev_(&jobv_,&jobv_,&n,A,&n,wr,wi,U,&n,Ui,&n,work,&lwork,&info);
    if( info != 0 )
        return info;

    //imaginary eigenvalues ignored since the Matrix is 'similar' to 
    //some symmetric matrix and will have the same eigenvalues 
    //          ~(Keilson 1979)
    mk_diag(wr,D, n,n);
    memcpy(A,U, n*n*sizeof(double));

    free(U); free(wr); free(wi); free(work);
    return 0;
}

/**  [median_h P l nl n]
 * Finds the likelihood vector [ret] for a child with likelihood [l]
 * using the probability matrix [p], and dimension [n]
 *
 * UNCHECKED EXCEPTION ::
 * ||[l]|| == ||[nl]|| == [n] == dim([P]) == col([P])
 *
 */
#ifdef _WIN32
__inline void 
#else
inline void 
#endif
median_h( const double* P, const double* l, double* nl, const int n )
{
    int i,j; double elm;
    //for each row of P
    for(i=0; i<n; ++i){
        elm = 0;
        //for each element in Pi and l
        for(j=0;j<n;++j)
            elm = elm + (P[i*n+j] * l[j]);
        nl[i] = elm;
    }
}

/** [loglikelihood l p]
 * finds the negative log likelihood of [l] with [p] as the priors
 *
 * UNCHECKED EXCEPTION ::
 * [n] == ||[l]|| == ||[p]||
 */
#ifdef _WIN32
__inline double
#else
inline double
#endif
loglikelihood( const double* l, const double* p, const int n )
{
    int j;
    double ret = 0;
    for(j=0;j<n;++j)
        ret = ret + (p[j] * l[j]);
    ret = -log( ret );
    return ret;
}

/** [mk_probmat_*** P U D [Ui] t]
 * Finds the probability matrix based on the diagonalized substitution
 * matrix, [U] and [D], with branch length [t]. returns result in [P].
 *
 * UNCHECKED EXCEPTION ::
 * [n]x[n] == dim([D]) == dim([P]) == dim([U]) == dim([Ui])
 *
 */
#define NTran  111
#define Tran   112
#define RMajor 101
#define CMajor 102
void
mk_probmat_sym( double* P, const double* U, const double* D, const float t, const int n)
{
    memcpy(P, D, n*n*sizeof(double) );
    apply_exp(P,n,n,t); //exp(D*t); along diagonal only

    double* I = (double*) malloc( n*n*sizeof(double));
    memset(I, 0, n*n*sizeof(double) );
    cblas_dgemm(RMajor,Tran,NTran, //format, op(A), op(B)
                n, n, n, 1,        //rows A, cols B, cols A, scalor multiplier of (A*B)
                U, n,              //MATRIX A, stride for A
                P, n,              //MATRIX B, stride for B
                1, I, n );         //scalor multiplier for C, MATRIX C, stride for C

    memset(P, 0, n*n*sizeof(double) );
    cblas_dgemm(RMajor,NTran,NTran, n,n,n,1, I,n, U,n, 1,P,n );
    free( I );
}
void
mk_probmat_gtr(double* P, const double* U, const double* D, const double* Ui, 
        const float t, const int n)
{
    memcpy(P,D,n*n*sizeof(double));
    apply_exp(P,n,n,t);

    double* S = (double*) malloc( n*n*sizeof(double));
    memset(S, 0, n*n*sizeof(double));
    cblas_dgemm(RMajor,NTran,NTran, n,n,n,1, U,n, P,n, 1,S,n);
    //S is U*exp(D)...

    memset(P,0,n*n*sizeof(double));
    cblas_dgemm(RMajor,NTran,Tran, n,n,n,1, S,n, Ui,n, 1,P,n);
    //P becomes U*expD*Ui... done
    free(S);
}

/** [median_c Pa Pb a b at bt c n]
 * Finds the likelihood vector [c] based on characters of its children [a] and
 * [b] with branch lenghts [at] and [bt] respectively, with the probability 
 * matrix, with branch lengths defined in [Pa] and [Pb] respectively.
 * 
 * UNCHECKED EXCEPTION ::
 * dim([Px]) == n == ||a|| == ||c|| == ||b|| == col([Px])
 */
void
median_char( const double* Pa, const double* Pb, 
                const double* a, const double* b, const int n, double* c)
{
    //temporary variable for first op, second op uses ret val
    double* nla= (double*) malloc(n*sizeof(double));
    //find likelihood vectors for each side
    median_h( Pa, a, nla, n);
    median_h(Pb, b, c, n);
    //pairwise multiplication...
    int i;
    for(i=0;i<n;++i)
        c[i] = nla[i] * c[i];
    free( nla );
}

/* ------------------------------------------------------ *
 *               (direct) INTERFACE WITH ML               *
 * ------------------------------------------------------ *
 */

/* Diagonlize matrix interfaces: symmetric */
CAMLprim value
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
CAMLprim value
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

/* general method for calculating likelihood of the current node */
CAMLprim value
likelihood_CAML_innerloop( const double* c_Pa, const double* c_Pb,
                                value a_chars, value b_chars, const int n)
{
    CAMLparam2( a_chars, b_chars );
    //loop over each element in a_chars + b_chars
    int i,k,code;
    k = Wosize_val(a_chars);
    struct c_ml* ptrs[k+1];

    //used for dimension of bigarray
    long dim[1];
    dim[0] = n;

    for(i=0;i<k;++i)
    {
        //convert next set of variables
        code = Int_val (Field( Field(a_chars,i), 0 ) );
        double* lv_b = (double*) Data_bigarray_val( Field(Field(b_chars,i),1));
        double* lv_a = (double*) Data_bigarray_val( Field(Field(a_chars,i),1));

        //create new vector and struct
        double* lv_n = (double*) malloc( n * sizeof(double));
        median_char(c_Pa,c_Pb,lv_a,lv_b,n,lv_n);

        //store it..
        struct c_ml *ptr = malloc( sizeof(struct c_ml) );
        ptr->cc = code;
        ptr->l_v = doublearray_to_bigarray( lv_n, 1, dim );
        ptrs[i] = ptr;
    }
    ptrs[k] = '\0'; //make array null terminated
    CAMLreturn( alloc_array(struct_to_value, ptrs) );
}

/** interface for SYM matrices **/
CAMLprim value
likelihood_CAML_median_wrapped_sym(value Q, value D, value at, value bt, 
                                            value a_chars, value b_chars )
{
    CAMLparam4(Q,D,at,bt); //a/b chars are done in inner loop

    //convert variables
    int n = Bigarray_val(Q)->dim[0];
    double* c_D = (double*) Data_bigarray_val( D );
    double* c_Q = (double*) Data_bigarray_val( Q );
    double c_at = Double_val(at);
    double c_bt = Double_val(bt);

    //make the probability matrices for a and b
    double* c_Pb = (double*) malloc( n*n*sizeof(double));
    double* c_Pa = (double*) malloc( n*n*sizeof(double));
    mk_probmat_sym(c_Pa,c_Q,c_D,c_at,n);
    mk_probmat_sym(c_Pb,c_Q,c_D,c_bt,n);

    return likelihood_CAML_innerloop(c_Pa, c_Pb, a_chars, b_chars, n);
}
/* argument wrapper */
CAMLprim value likelihood_CAML_median_sym(value * argv, int argn){
    return likelihood_CAML_median_wrapped_sym
            (argv[0],argv[1],argv[2],argv[3],argv[4],argv[5]); 
}

/* interface for GTR matrices */
CAMLprim value
likelihood_CAML_median_wrapped_gtr(value Q, value D, value Qi, value at,
                                     value bt, value a_chars, value b_chars)
{
    CAMLparam5(Q,D,Qi,at,bt); //a/b chars are done in inner loop

    //convert variables
    int n = Bigarray_val(Q)->dim[0];
    double* c_D = (double*) Data_bigarray_val( D );
    double* c_Q = (double*) Data_bigarray_val( Q );
    double* c_Qi= (double*) Data_bigarray_val( Qi);
    double c_at = Double_val(at);
    double c_bt = Double_val(bt);

    //make the probability matrices for a and b
    double* c_Pb = (double*) malloc( n*n*sizeof(double));
    double* c_Pa = (double*) malloc( n*n*sizeof(double));
    mk_probmat_gtr(c_Pa,c_Q,c_D,c_Qi,c_at,n);
    mk_probmat_gtr(c_Pb,c_Q,c_D,c_Qi,c_bt,n);

    return likelihood_CAML_innerloop(c_Pa,c_Pb,a_chars,b_chars,n);
}
/* argument wrapper */
CAMLprim value likelihood_CAML_median_gtr(value * argv, int argn){
    return likelihood_CAML_median_wrapped_gtr
            (argv[0],argv[1],argv[2],argv[3],argv[4],argv[5],argv[6]);
}
