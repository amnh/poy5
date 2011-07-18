/* POY 4.0 Beta. A phylogenetic analysis program using Dynamic Homologies.    */
/* Copyright (C) 2007  Andrés Varón, Le Sy Vinh, Illya Bomash, Ward Wheeler,  */
/* and the American Museum of Natural History.                                */
/*                                                                            */
/* This program is free software; you can redistribute it and/or modify       */
/* it under the terms of the GNU General Public License as published by       */
/* the Free Software Foundation; either version 2 of the License, or          */
/* (at your option) any later version.                                        */
/*                                                                            */
/* This program is distributed in the hope that it will be useful,            */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of             */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the              */
/* GNU General Public License for more details.                               */
/*                                                                            */
/* You should have received a copy of the GNU General Public License          */
/* along with this program; if not, write to the Free Software                */
/* Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301   */
/* USA                                                                        */

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

#include "seq.h"
#include "fm.h"
#include "matrices.h"

struct fcm {
    fm*        fmat;  /** costs and assignments for alignment           **/
    matricest  direc; /** Matrix of the directions attached to the cost **/
    double*    costs; /** Costs for the substring alignment             **/
    int*       nukk;  /** Number of gaps in the subseq for the median   **/
};
typedef struct fcm fcmt;

double full_falign( const seqt x, const seqt y, fcmt *FA );

void full_backtrace (const seqt x, const seqt y, seqt e1, seqt e2, seqt m, fcmt *FA);

#endif
