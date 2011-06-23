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

#ifdef USE_LIKELIHOOD

#include <caml/alloc.h>     //copy_double, et cetera
#include <caml/mlvalues.h>
#include <caml/memory.h>    //caml_param, et cetera
#include <caml/bigarray.h>
#include <caml/custom.h>    //abstract types
#include <caml/intext.h>    //serialization
#include <caml/fail.h>      //failwith('')

#include "falgn.h"

#define GET(x,i) seq_get(x,i)

double min3(CDIR* d_ptr, CASN* a_ptr, double a,CASN aasn, double i,CASN iasn, double d,CASN dasn )
{
    *a_ptr = aasn;
    *d_ptr = ALIGN;
    return a;
}

double min2_ins(CDIR* d_ptr, CASN* a_ptr, double a,CASN aasn, double i,CASN iasn )
{
    if( a == i ){
        *a_ptr = aasn | iasn;
        *d_ptr = ALIGN | INSERT;
        return a;
    } else if ( a < i ){
        *a_ptr = aasn;
        *d_ptr = ALIGN;
        return a;
    } else {
        *a_ptr = iasn;
        *d_ptr = INSERT;
        return i;
    }
}

double min2_del(CDIR* d_ptr, CASN* a_ptr, double a,CASN aasn, double d,CASN dasn )
{
    if( a == d ){
        *a_ptr = aasn | dasn;
        *d_ptr = ALIGN | DELETE;
        return a;
    } else if ( a < d ){
        *a_ptr = aasn;
        *d_ptr = ALIGN;
        return a;
    } else {
        *a_ptr = dasn;
        *d_ptr = DELETE;
        return d;
    }
}

double cost( fm *FM, CASN* a_ptr, CASN x, CASN y )
{
    if( 1 == FM->comb ){
        *a_ptr = FM->cost_asgn[ x * FM->size + y ];
        return FM->cost[ x * FM->size + y ];
    } else {
        return 0.0;
        //return calculate_cost( a_ptr, x, y );
    }
}

/** Do a full alignment; filling the entire cost matrix; FM is allocated **/
void full_falign( seqt x, seqt y, fcmt *FA )
{
    int i, j, place;
    double cost_algn, cost_ins, cost_del;
    CASN   asgn_algn, asgn_ins, asgn_del;

    CASN temp_asn;
    CDIR temp_dir;

    for( i = 0; i < x->len; ++i ){
        for( j = 0; j < y->len; ++j ){
            place = i*FA->fmat->size + j;
            cost_algn = cost( FA->fmat, &asgn_algn, GET(x,i), GET(y,j));
            cost_ins  = cost( FA->fmat, &asgn_ins,  GET(x,i), FA->fmat->gap );
            cost_del  = cost( FA->fmat, &asgn_del,  FA->fmat->gap , GET(y,j));
            FA->costs[ place ] = 
                min3( &temp_dir, &temp_asn, 
                    /* ALIGNMENT */ cost_algn, asgn_algn,
                    /* INSERTION */ cost_ins,  asgn_ins,
                    /* DELETION  */ cost_del,  asgn_del );
            FA->assgn[ place ] = temp_asn;
            FA->direc[ place ] = temp_dir;
        }
    }
}

void update_exclude_row( fcmt *FM, int i, int j )
{

}

void update_exclude_col( fcmt *FM, int i, int j )
{

}

void update_all( fcmt *FM, int i, int j )
{

}

void ukk_align( seqt x, seqt y, fm *FM )
{

}

#endif
