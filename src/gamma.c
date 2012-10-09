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
#include <stdlib.h> /* malloc, calloc, srand, RAND_MAX */
#include <string.h> /* memcpy, memset */
#include <assert.h>
#include <malloc.h>

#include "config.h" //defines if likelihood, use_.., et cetera
#include <math.h>   //log,exp
#include <time.h>   //for random init

//caml specific headers
#include <caml/alloc.h>     //copy_double, et cetera
#include <caml/mlvalues.h>
#include <caml/memory.h>
#include <caml/bigarray.h>
#include <caml/fail.h>

#define EPSILON 3e-7  //error for numerical calculations
#define FABMIN  1e-30 //floating point absolute value min  

#ifndef M_PI
#define M_PI    3.14159265358979323846264338327
#endif

#define CHECK_MEM(a) if(a==NULL) failwith("I can't allocate more memory.")

void CHECK_MEAN(double*a, int n){
    int Z;
    double SUM;
    for( Z=0,SUM=0; Z<n; Z++){ SUM += a[Z]; }
    if( SUM/(double)n > 1.0+EPSILON ){
        failwith("Incorrect Mean of Gamma Rates");
    }
}

//------------------------------------------------------------
//  SPECIAL FUNCTIONS -- 
//------------------------------------------------------------

/** [gamma z]
 * Computes the gamma function of z using Lanczos Approximation. (e~1e-15)
 *      caveate -- max z = 142.0
 */
double gamma( double z )
{
    double x,Lg; int k;
    static double C[] = 
    {  0.99999999999980993, 676.5203681218851, -1259.1392167224028,
        771.32342877765313, -176.61502916214059, 12.507343278686905,
        -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7 };

    if( z < .5 )
        return M_PI / (sin(M_PI*z)*gamma(1-z));
    z--;
    Lg = C[0];
    for(k=1;k<9;++k)
        Lg = Lg + (C[k] / (z+k));

    x = z + 7.5;
    x = sqrt(2*M_PI) * pow(x, z+0.5) * exp(-x) * Lg;
    return (x);
}
value gamma_CAML_gamma( value v_x ) 
{
    CAMLparam1( v_x ); 
    CAMLlocal1( v_g );
    v_g = caml_copy_double( gamma( Double_val( v_x ) ) );
    CAMLreturn( v_g );
}

/** [lngamma z]
 * Computes the ln gamma function of z using Lanczos Approximation.
 *
 * "Numerical Recipes in C", Section 6.1
double lngamma( const double xx )
{
    double x,y,tmp,ser;
    int j;
    static double cof[6] =
    {
      76.18009172947146, -86.50532032941677, 24.01409824083091,
      -1.231739572450155, 0.1208650973866179e-2, -0.5395239384953e-5 };

    y = x = xx;
    tmp = x+5.5;
    tmp -= (x+0.5) * log(tmp);
    ser=1.000000000190015;
    for(j=0;j<=5;j++) ser += cof[j]/++y;
    return -tmp+log(2.5066282746310005*ser/x);
}
*/

double lngamma( const double xx )
{
    int j;
    double x,tmp,y,ser;
    static const double cof[14] = { 
        57.1562356658629235,    -59.5979603554754912,     14.1360979747417471, 
        -0.491913816097620199,    0.339946499848118887e-4, 0.465236289270485756e-4,
        -0.983744753048795646e-4, 0.158088703224912494e-3,-0.210264441724104883e-3,
         0.217439618115212643e-3,-0.164318106536763890e-3, 0.844182239838527433e-4,
        -0.261908384015814087e-4, 0.368991826595316234e-5
    };
    if (xx <= 0) failwith("bad arg in gammln");
    y=x=xx;
    tmp = x+5.24218750000000000; /* Rational 671/128 */
    tmp = (x+0.5)*log(tmp)-tmp;
    ser = 0.999999999999997092;
    for (j=0;j<14;j++) ser += cof[j]/++y;
    return tmp+log(2.5066282746310005*ser/x);
}
value gamma_CAML_lngamma( value v_x ) 
{
    CAMLparam1( v_x ); 
    CAMLlocal1( v_lng );
    v_lng = caml_copy_double( lngamma( Double_val( v_x ) ) );
    CAMLreturn( v_lng );
}

/** [gamma_pdf r alpha beta]
 * Return the probability density function of the gamma distribution */
double gamma_pdf(const double r, const double alpha, const double beta)
{
    return (pow(beta,alpha)*pow(r, alpha-1))/(exp(beta*r)*gamma(alpha));
}

/** [lngamma_pdf r alpha beta]
 * Return the probability density function of the gamma distribution using ln
 * gamma. Since, exp(b*r)*gamma(a) = exp(b*r)*e(lngam(a)) = exp(b*r+lngam(a)) */
double lngamma_pdf(const double r, const double alpha, const double beta)
{
    return (pow(beta,alpha)*pow(r, alpha-1)) / (exp(beta*r + lngamma(alpha)) );
}

/** [rand_normal m s]
 * generate a random number in a given gaussian/normal distribution */
double rand_normal( const double mean, const double stdev )
{
    double u1,u2, r,theta;
    assert( stdev > 0.0 );
    srand( time(NULL) );

    u1 = rand();
    u2 = rand();
    r = sqrt( -2.0 * log(u1) );
    theta = 2.0 * M_PI * u2;

    return (mean + stdev * r * sin(theta));
}
value gamma_CAML_randnormal( value m, value s )
{
    CAMLparam2( m, s );
    CAMLlocal1( r );
    r = caml_copy_double( rand_normal( Double_val(m), Double_val(s) ) );
    return r;
}

/* [rand_exponential mean]
 * generate a random exponential value from a mean */
double rand_exp( const double mean )
{
    assert( mean > 0.0 );
    srand( time(NULL) );
    return (-mean * log(rand()));
}
value gamma_CAML_randexp( value m )
{
    CAMLparam1( m );
    CAMLlocal1( r );
    r = caml_copy_double( rand_exp( Double_val(m) ) );
    return r;
}

/** [rand_gamma a b]
 * Implementation based on "A Simple Method for Generating Gamma Variables"
 * by George Marsaglia and Wai Wan Tsang.  
 * ACM Transactions on Mathematical Software
 * Vol 26, No 3, September 2000, pages 363-372.
 */
double rand_gamma( const double shape, const double scale )
{
    double g,w,x,v,c,d,xsq,u;

    assert( shape > 0.0 );
    assert( scale > 0.0 );
    srand( time(NULL) );

    if( shape >= 1.0 ){
        d = shape - 1.0 / 3.0;
        c = 1.0 / (sqrt( 9.0 * d));
        do {
            x = rand_normal(0, 1);
            v = 1.0 + c *x;
            while( v <= 0.0 ){
                x = rand_normal(0, 1);
                v = 1.0 + c * x;
            }
            v = v*v*v;
            u = rand();
            xsq = x*x;
        } while((u < (1.0 - 0.331 * xsq * xsq))
            || (log(u) < (0.5*xsq + d*(1.0-v+log(v))))); 
    } else {
        g = rand_gamma( shape + 1.0, 1.0 );
        w = rand();
        return (scale*g*pow(w,1.0/shape));
    }
    return (scale * d * v);
}
value gamma_CAML_randgamma( value sh, value sc )
{
    CAMLparam2( sh,sc );
    CAMLlocal1( r );
    r = caml_copy_double( rand_gamma( Double_val(sh), Double_val(sc) ) );
    return r;
}



/** confluent hypergeometric (for incomplete gamma ratio)
 *   ___inf
 *   \       ______x^r_______ 
 *   /      (p+1)(p+2)...(p+r)
 *   ---r=1
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
} */


/** [gamma_i x a ]
 * calculates the incomplete gamma ratio, based on methods described in
 *      ~Bhattacharjee (1970), Algorithm AS32
 *
 *      a - alpha parameter
 *      x - bound of integral
double gamma_i(const double x, const double a)
{
    //if ( x < p || (a <= x && x <= 1.0) ) {
        return (exp(-x)*pow(x,a)/gamma(a+1)) * gamma_M(x,a);
    //else
    //    return 1 - exp(-x)*pow(x,a) / gamma(a) * gamma_C(x,a);
} */


/**
 * Retuns the incomplete gamma function evaluated by its series representation
 * also, ln gamma in gln.
 *
 * "Numerical Recipes in C", Section 6.2
 */
void gser(double *gam,double a,double x,double *gln)
{
    double sum,del,ap;
    *gln = lngamma ( a );
    if (x <= 0.0 ) { *gam = 0.0; return; }
    ap = a;
    del = sum = 1.0 / a;
    for(;;) {
        ++ap;
        del *= x/ap;
        sum += del;
        if (fabs(del) < fabs(sum)*EPSILON){
            *gam = sum * exp(-x+a*log(x)-(*gln));
            break;
        }
    }
}

/**
 * Returns the incomplete gamma function evaluated by its continued fraction
 * also, ln gamma in gln.
 *
 * "Numerical Recipes in C", Section 6.2
 */
void gcf( double *gam,double a,double x,double *gln )
{
    int i;
    double an,b,c,d,del,h;

    *gln = lngamma( a );
    b = x + 1.0 - a;
    c = 1.0 / FABMIN;
    d = 1.0 / b;
    h=d;
    for(i=1;;i++){
        an = -i*(i-a);
        b += 2.0;
        d = an*d+b;
        if ( fabs(d) < FABMIN ){ d=FABMIN; }
        c = b+an/c;
        if ( fabs(c) < FABMIN ){ c=FABMIN; }
        d = 1.0 / d;
        del = d*c;
        h *= del;
        if (fabs(del-1.0) <= EPSILON){ break; }
    }
    *gam = exp(-x+a*log(x)-(*gln))*h;
}

/** 
 * Returns the incomplete gamma function P(a,x)
 *
 * "Numerical Recipes in C", Section 6.2
 */
double gammap( const double x, const double a )
{
    double gam,gln;
    if (x < 0.0 || a <= 0.0)
        failwith("Invalid argument for incomplete gamma");
    if ( x < (a+1.0)) {
        gser(&gam,a,x,&gln);
    } else {
        gcf (&gam,a,x,&gln); //inverse of continued fraction representation
        gam = 1.0 - gam;
    }
    return gam;
}

/** [point_normal p]
 * Finds the percentage point [p] of the normal distribution.
 *
 * Algorithm AS241: The Percentage Points of the Normal Distribution.
 *     Accurate to about 1 part in 10**16
 */
#ifdef _WIN32
__inline
#else
inline
#endif
double point_normal_eq( const double *a, const double *b, const double r, const double q){
    double numr, deno;
    numr = ((((((a[7]*r+a[6])*r+a[5])*r+a[4])*r+a[3])*r+a[2])*r+a[1])*r+a[0];
    deno = ((((((b[7]*r+b[6])*r+b[5])*r+b[4])*r+b[3])*r+b[2])*r+b[1])*r+1.0 ;
    return q * numr / deno;
}

double point_normal( double p ){
    double q, r, ppnd;
    static double a[8] = {
        3.3871328727963666080,
        1.3314166789178437745e+2,
        1.9715909503065514427e+3,
        1.3731693765509461125e+4,
        4.5921953931549871457e+4,
        6.7265770927008700853e+4,
        3.3430575583588128105e+4,
        2.5090809287301226727e+3
    };
    static double b[8] = { 
        0.0,
        4.2313330701600911252e+1,
        6.8718700749205790830e+2,
        5.3941960214247511077e+3,
        2.1213794301586595867e+4,
        3.9307895800092710610e+4,
        2.8729085735721942674e+4,
        5.2264952788528545610e+3
    };
    static double c[8] = {
        1.42343711074968357734,
        4.63033784615654529590,
        5.76949722146069140550,
        3.64784832476320460504,
        1.27045825245236838258,
        2.41780725177450611770e-1,
        2.27238449892691845833e-2,
        7.74545014278341407640e-4
    };
    static double d[8] = { 
        0.0,
        2.05319162663775882187,
        1.67638483018380384940,
        6.89767334985100004550e-1,
        1.48103976427480074590e-1,
        1.51986665636164571966e-2,
        5.47593808499534494600e-4,
        1.05075007164441684324e-9
    };
    static double e[8] = {
        6.65790464350110377720,
        5.46378491116411436990,
        1.78482653991729133580,
        2.96560571828504891230e-1,
        2.65321895265761230930e-2,
        1.24266094738807843860e-3,
        2.71155556874348757815e-5,
        2.01033439929228813265e-7
    };
    static double f[8] = { 
        0.0,
        5.99832206555887937690e-1,
        1.36929880922735805310e-1,
        1.48753612908506148525e-2,
        7.86869131145613259100e-4,
        1.84631831751005468180e-5,
        1.42151175831644588870e-7,
        2.04426310338993978564e-15
    };
  
    q = p - 0.5;
    if (fabs(q) <= 0.425)
        return point_normal_eq( a, b, (0.180625-q*q), q);

    r = p;
    if (q >= 0.0) r = 1.0 - r;
    if (r <= 0.0) return 0.0;

    r = sqrt(-log(r));
    if (r <= 0.5)
        ppnd = point_normal_eq( c, d, (r-1.6), 1.0);
    else
        ppnd = point_normal_eq( e, f, (r-0.5), 1.0);

    if (q < 0.0) ppnd = -ppnd;
    return ppnd;
}

/** [chi_pp p v]
 * Finds the percentage point [p] of the chi squared distribution of [v] degrees
 * of freedom.  The gamma is related to this distribution by the define below.
 *
 * Algorithm AS91: The Percentage Points of the chi^2 Distribution
 *      (translated to C, and removed goto's ~nrl) */
double chi_pp( double p, double v ){
    double ch,s1,s2,s3,s4,s5,s6;
    double e,aa,xx,c,g,x,p1,a,q,p2,t,ig,b;

    assert( v > 0.0 );
    if (p < 0.000002 || p > 0.999998) failwith("Chi^2 Percentage Points incorrect.1");

    e = 0.5e-6;         /** error term **/
    aa= 0.6931471805;
    xx = 0.5 * v;
    c  = xx - 1.0;
    g  = lngamma( xx );

    if( v < -1.24 * log(p) ){
        ch = pow(p * xx * exp ( g + xx * aa), 1.0/xx);
        if( ch - e < 0 ) return ch;
    } else if( v > 0.32) {
        x = point_normal( p );
        p1 = 0.222222 / v;
        ch = v * pow( x * sqrt( p1 ) + 1 - p1, 3.0) ;
        if (ch > 2.2 * v + 6)
            ch = -2.0 * (log(1-p) - c*log(0.5*ch)+g);
    } else {
        ch  = 0.4;
        a = log (1 - p);
        do{
            q  = ch;
            p1 = 1 + ch * (4.67 + ch);
            p2 = ch * (6.73 + ch * (6.66 + ch));
            t  = -0.5 + (4.67 + 2*ch)/p1 - (6.73 + ch*(13.32 + 3*ch))/p2;
            ch = ch - (1- exp( a + g + 0.5*ch+c*aa) * p2/p1)/t;
        } while( fabs( q/ch - 1) - 0.01 > 0.0 );
    }

    do{
        q  = ch;
        p1 = .5*ch;
        ig = gammap( p1, xx );
        if (ig < 0){ failwith("Chi^2 Percentage Points incorrect.2"); }
        p2 = p - ig;
        t  = p2 * exp( xx*aa + g + p1 - c*log(ch));
        b  = t / ch;
        a  = (0.5*t) - (b*c);
        /* Seven terms of the Taylor series */
        s1 = (210 + a*(140 + a*(105 + a*(84 + a*(70 + 60*a))))) / 420.0;
        s2 = (420 + a*(735 + a*(966 + a*(1141 + 1278*a)))) / 2520.0;
        s3 = (210 + a*(462 + a*(707 + 932*a))) / 2520.0;
        s4 = (252 + a*(672 + 1182*a) + c*(294 + a*(889 + 1740*a))) / 5040.0;
        s5 = ( 84 + 264*a + c*(175 + 606*a)) / 2520.0;
        s6 = (120 + c*(346 + 127*c)) / 5040.0;
        ch+= t*(1+0.5*t*s1-b*c*(s1-b*(s2-b*(s3-b*(s4-b*(s5-b*s6))))));
    } while( fabs(q / ch - 1.0) > e);

    return (ch);
}
#define gamma_pp(prob,alpha,beta) chi_pp(prob,2.0*alpha)/(2.0*beta)

/** [gamma_rates rates alpha beta cuts k]
 * calculates the rates into [rates] with percentage points [cuts] of [k]
 * categories and shape parameters [alpha] and [beta] */ 
void 
gamma_rates(double* rates,const double a,const double b,const double* cuts,const int k)
{
    double fac, *ingam;
    int j;
    fac = a*((double)k)/b;
    ingam = (double*) malloc( k * sizeof(double));
    CHECK_MEM(ingam);

    //calculate: rj = (A*k/B)*(I(bB,A+1)-I(aB,A+1))
    for(j=0;j<(k-1);j++){
        ingam[j] = gammap( cuts[j]*b, a+1 );
        //printf("LG: %f\tNL: %f\n", ingam[j], gamma_i( cuts[j]*b, a+1 ) );
    }
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
    long dims[1];

    alpha = Double_val( a );
    beta =  Double_val( b );
    cats = Int_val( c );

    assert( cats > 0 );
    rate_ray = (double*) malloc( sizeof(double)*cats ); 
    CHECK_MEM(rate_ray);
    
    if( 1 == cats ){
        rate_ray[0] = 1.0;
    } else {
        pcut_ray = (double*) malloc( sizeof(double)*cats );
        CHECK_MEM(pcut_ray);
        for(j=1;j<cats;++j)
            pcut_ray[j-1] = gamma_pp( (double)j/(double)cats, alpha, beta );
        gamma_rates( rate_ray, alpha, beta, pcut_ray, cats );
        free( pcut_ray );
    }

    dims[0] = cats;
    rates = alloc_bigarray(BIGARRAY_FLOAT64 | BIGARRAY_C_LAYOUT,1,rate_ray,dims);
    CAMLreturn( rates );
}

