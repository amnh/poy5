#include <stdio.h>
#include <stdlib.h> //malloc, srand, RAND_MAX
#include <string.h> //memcpy, memset
#include <assert.h>
#include <math.h>   //log,exp,sin

#include <caml/mlvalues.h>
#include <caml/memory.h>
#include <caml/bigarray.h>
#include <caml/custom.h>
#include <caml/intext.h>

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

