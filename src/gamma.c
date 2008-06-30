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
#include <math.h>   //log,exp
#include <caml/mlvalues.h>
#include <caml/memory.h>
#include <caml/bigarray.h>

#define STEP 1e-4     //for numerical integration
#define EPSILON 1e-4  //error for numerical calculations

#define CHECK_MEM(a) if(a==NULL) failwith("I can't allocate more memory.")

/* -- check that the sum of the mean == 1.0 +/- EPSILON */
//#define CHECK_MEAN(a,n); int Z;double SUM;for(Z=0,SUM=0;Z<n;Z++){ SUM += a[Z]; }\
                            if( SUM/(double)n > 1.0+EPSILON ){ \
                                failwith("Incorrect Mean of Gamma Rates"); }
//#define CHECK_MEAN(a,n); int Z;double SUM;for(Z=0,SUM=0;Z<n;Z++){ SUM += a[Z]; }\
                            if( SUM/(double)n > 1.0+EPSILON ){ \
                                printf("Mean of Rates Error :: %f\n",SUM/(double)n); }
#define CHECK_MEAN(a,b);

//------------------------------------------------------------
//  SPECIAL FUNCTIONS -- 
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
