/* POY 5.1.1. A phylogenetic analysis program using Dynamic Homologies.       *\
(* Copyright (C) 2011  Andrés Varón, Lin Hong, Nicholas Lucaroni, Ward Wheeler*)
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

#include "config.h"         //defines, if likelihood, USE_LIKELIHOOD
#ifdef USE_LIKELIHOOD   
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

#endif /* USE_LIKELIHOOD */
