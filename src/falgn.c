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

#include <assert.h>         

#include <caml/alloc.h>     //copy_double, et cetera
#include <caml/mlvalues.h>
#include <caml/memory.h>    //caml_param, et cetera
#include <caml/bigarray.h>
#include <caml/intext.h>    //serialization
#include <caml/fail.h>      //failwith('')

#define Val_none    Val_int(0)
#define Some_val(v) Field(v,0)

#include "falgn.h"

/** constants **/
#define EPSILON         1e-6
#ifndef INFINITY
    #define INFINITY    1/0
#endif

/** dor this file CDIR == DIRECTION_MATRIX in matiches.h **/

/** defines for ease of use **/
#define CDIR         DIRECTION_MATRIX
#define MAX(a,b)     (a>=b)?a:b
#define MIN(a,b)     (a<=b)?a:b
#define EQUALS(a,b)  fabs(a-b)<EPSILON
#define ISSET(x,i)   0<((1<<i)&x)
#define GET(x,i)     seq_get(x,i)

#define DSET(f,p,d)  f->direc->matrix_d[p]=(CDIR)d
#define DGET(f,p)    f->direc->matrix_d[p]

//used for debugging with electric fense; we cannot use scratch space
/*#define register_section(x,a,c) (double*) malloc( sizeof(double)*a)*/

void
print_alignment_matrixf(const double* costs, matricest mat, const int w, const int h)
{
    int i,j;
    for (i = 0; i < h; i++) {
        for (j = 0; j < w; j++)
            fprintf (stdout, "| %d %8.5f[00] ", mat->matrix_d[(w*i)+j],costs[(w*i)+j]);
        fprintf (stdout, "|\n");
    }
    fprintf (stdout, "\n");
}

double min2_ins(CDIR* d_ptr, const double a, const double i)
{
    if( EQUALS(a,i) ){
        *d_ptr = ALIGN | INSERT;
        return (MIN(a,i));
    } else if ( a < i ){
        *d_ptr = ALIGN;
        return a;
    } else {
        *d_ptr = INSERT;
        return i;
    }
}

double min2_del( CDIR* d_ptr, const double a, const double d)
{
    if( EQUALS(a,d) ){
        *d_ptr = ALIGN | DELETE;
        return (MIN(a,d));
    } else if ( a < d ){
        *d_ptr = ALIGN;
        return a;
    } else {
        *d_ptr = DELETE;
        return d;
    }
}

double min3( CDIR* d_ptr, const double a, const double i, const double d )
{
    double min_ai;
    min_ai = min2_ins( d_ptr, a, i );
    if( EQUALS(min_ai,d) ){
        *d_ptr = *d_ptr | DELETE;
        return (MIN(d,min_ai));
    } else if ( d < a ){
        *d_ptr = DELETE;
        return d;
    } else {
        return min_ai;
    }
}

/** Create a full backtrace with median and edited sequences; polymorphic **/
void
median_backtrace (const seqt x, const seqt y, seqt m, fcmt *FA)
{
    int x_loc, y_loc, linear;
    x_loc = x->len - 1;
    y_loc = y->len - 1;

    seq_clear( m  );
    while( !(x_loc == 0 && y_loc == 0) )
    {
        linear = x_loc * y->len + y_loc;
/*        printf("(%d,%d):%d -- %d !e { A:%d, D:%d, I:%d }\n",
                x_loc,y_loc,linear,DGET(FA,linear),ALIGN,DELETE,INSERT); */
        if( INSERT & DGET(FA,linear) ){
            seq_prepend( m , fasn(FA->fmat,FA->fmat->gap,GET(y,y_loc)) );
            --y_loc;
        } else if( ALIGN & DGET(FA,linear) ){
            seq_prepend( m , fasn(FA->fmat,GET(x,x_loc),GET(y,y_loc)) );
            --x_loc;
            --y_loc;
        } else if( DELETE & DGET(FA,linear) ){
            seq_prepend( m , fasn(FA->fmat,GET(x,x_loc),FA->fmat->gap) );
            --x_loc;
        } else {
            assert( 1 == 0 );
        }
    }
    seq_prepend( m, FA->fmat->gap); 
}


/** Create a full backtrace with median and edited sequences; polymorphic **/
void
full_backtrace (const seqt x, const seqt y, seqt e1, seqt e2, seqt m, fcmt *FA)
{
    int x_loc, y_loc, linear;
    x_loc = x->len - 1;
    y_loc = y->len - 1;

    seq_clear( m  );
    seq_clear( e1 );
    seq_clear( e2 );

    while( !(x_loc == 0 && y_loc == 0) )
    {
        linear = x_loc * y->len + y_loc;
/*        printf("(%d,%d):%d -- %d !e { A:%d, D:%d, I:%d }\n",
                x_loc,y_loc,linear,DGET(FA,linear),ALIGN,DELETE,INSERT); */
        if( INSERT & DGET(FA,linear) ){
            seq_prepend( e1, FA->fmat->gap );
            seq_prepend( e2, GET(y, y_loc) );
            seq_prepend( m , fasn(FA->fmat,FA->fmat->gap,GET(y,y_loc)) );
            --y_loc;
        } else if( ALIGN & DGET(FA,linear) ){
            seq_prepend( e1, GET(x, x_loc) );
            seq_prepend( e2, GET(y, y_loc) );
            seq_prepend( m , fasn(FA->fmat,GET(x,x_loc),GET(y,y_loc)) );
            --x_loc;
            --y_loc;
        } else if( DELETE & DGET(FA,linear) ){
            seq_prepend( e1, GET(x, x_loc) );
            seq_prepend( e2, FA->fmat->gap );
            seq_prepend( m , fasn(FA->fmat,GET(x,x_loc),FA->fmat->gap) );
            --x_loc;
        } else {
            assert( 1 == 0 );
        }
    }
    seq_prepend( m, FA->fmat->gap); 
    seq_prepend(e2, FA->fmat->gap); 
    seq_prepend(e1, FA->fmat->gap); 
}

/** Determine the cost of the alignment of A and B across branch length t **/
#ifdef _WIN32
__inline void 
#else
inline void 
#endif
full_cost(const seqt x, const seqt y, fcmt *FA, const int i, const int j)
{
    double cost_algn, cost_ins, cost_del, temp_min;
    int place;
    CDIR temp_dir;
    place = (y->len * i) + j;
    if( 0 == i && 0 == j){
        FA->costs[ place ] = 0.0;
    } else if( 0 == i ){
        FA->costs[ place ] = FA->costs[ j-1 ] + fcost( FA->fmat, FA->fmat->gap, GET(y,j) );
        DSET( FA, place, INSERT );
    } else if( 0 == j ){
        FA->costs[ place ] = FA->costs[ (i-1)*y->len ] + fcost( FA->fmat, GET(x,i), FA->fmat->gap );
        DSET( FA, place, DELETE );
    } else {
        cost_algn = FA->costs[(i-1)*y->len + j-1] + fcost( FA->fmat, GET(x,i),      GET(y,j)      );
        cost_ins  = FA->costs[ i   *y->len + j-1] + fcost( FA->fmat, FA->fmat->gap, GET(y,j)      );
        cost_del  = FA->costs[(i-1)*y->len + j  ] + fcost( FA->fmat, GET(x,i),      FA->fmat->gap );
        temp_min = min3( &temp_dir, cost_algn, cost_ins, cost_del );
/*        printf("(%d,%d): A:%f\tI:%f\tD:%f --> %f[%d]\n",i,j,cost_algn,cost_ins,cost_del,temp_min,temp_dir); */
        FA->costs[ place ] = temp_min;
        DSET( FA, place, temp_dir );
    }
}

/** Do a full alignment; filling the entire cost matrix; FCM is allocated **/
double full_falign( const seqt x, const seqt y, fcmt *FA )
{
    int i, j;
    for( i = 0; i < x->len; ++i ){
        for( j = 0; j < y->len; ++j ){
            full_cost( x, y, FA, i, j );
        }
    }
    return FA->costs[ x->len*y->len - 1 ];
}

void update_exclude_row( fcmt *FA, seqt x, seqt y, const int i, const int j )
{
    int place;
    double cost_algn, cost_del;
    CDIR temp_dir;

    place = i * FA->fmat->size + j;
    cost_algn = fcost( FA->fmat, GET(x,i), GET(y,j) );
    cost_del  = fcost( FA->fmat, FA->fmat->gap , GET(y,j) );
    FA->costs[ place ] = min2_del( &temp_dir, cost_algn, cost_del );
    DSET( FA, place, temp_dir );
}

void update_exclude_col( fcmt *FA, seqt x, seqt y, const int i, const int j )
{
    int place;
    double cost_algn, cost_ins;
    CDIR temp_dir;

    place = i * FA->fmat->size + j;
    cost_algn = fcost( FA->fmat, GET(x,i), GET(y,j) );
    cost_ins  = fcost( FA->fmat, GET(x,i), FA->fmat->gap );
    FA->costs[ place ] = min2_ins( &temp_dir, cost_algn, cost_ins );
    DSET( FA, place, temp_dir );
}

void update_all( fcmt *FA, seqt x, seqt y, const int i, const int j )
{
    int place;
    double cost_algn, cost_ins, cost_del;
    CDIR temp_dir;

    place = i * FA->fmat->size + j;
    cost_algn = fcost( FA->fmat, GET(x,i), GET(y,j) );
    cost_ins  = fcost( FA->fmat, GET(x,i), FA->fmat->gap );
    cost_del  = fcost( FA->fmat, FA->fmat->gap , GET(y,j) );
    FA->costs[ place ] = min3( &temp_dir, cost_algn, cost_ins, cost_del );
    DSET( FA, place, temp_dir );
}

double nukk_falign( const seqt x, const seqt y, fcmt *FA )
{
    return 0.0;
}

value
falign_CAML_median(value oSpace, value oMat, value oU, value oD, value oUi,
                   value ota, value otb, value oa, value ob, value om, value ogap)
{
    CAMLparam5( oSpace, oMat, oU, oD, oUi );
    CAMLxparam5( ota, otb, oa, ob, om );
    CAMLxparam1( ogap );

    seqt a,b, m;
    double ta, tb, *PA, *PB, *TMP;
    fcmt results;
    mat *space;
    int mat_size, alph, fcm_size, gap;

    Seq_custom_val(a,oa);
    Seq_custom_val(b,ob);
    Seq_custom_val(m, om);
    space = FM_val( oSpace );
    results.direc = Matrices_struct( oMat );
    alph = Bigarray_val(oU)->dim[0];
    gap = Int_val( ogap );

    ta = Double_val( ota );
    tb = Double_val( otb );

    /* we need to ensure that we didn't "clear" the memory earlier; although we
     * never actually clear it, we need to at least ensure that an alignment is
     * stored. To re-register that section, since it must start at 0, 'free' and
     * register; this allows us to re-register space for composing matrices.
     *      TODO: can we be safe and re-use these cost matrices? */
    mat_size = a->len * b->len;
    assert( space->loc > mat_size );
    free_all( space );
    fcm_size = ((1 << alph)-1);

    results.fmat = (fm*) malloc( sizeof(fm) );
    results.fmat->alph = alph;
    results.fmat->size = fcm_size;
    results.fmat->comb = 1;
    results.fmat->gap  = 1 << gap;
    results.costs = register_section( space, mat_size, 0 );
    mat_setup_size (results.direc, a->len, b->len, 0, 0, alph, 0);
    results.fmat->cost = register_section( space, fcm_size*fcm_size, 0 );
    results.fmat->cost_asgn = (int*) malloc( sizeof(int)*fcm_size*fcm_size );

    /** pre-calculate the matrix **/
    PA = register_section( space, alph*alph, 0 );
    PB = register_section( space, alph*alph, 0 );
    TMP= register_section( space, alph*alph, 0 );
    if( Val_none == oUi ){
        compose_sym( PA, Data_bigarray_val(oU), Data_bigarray_val(oD), ta, alph, TMP );
        compose_sym( PB, Data_bigarray_val(oU), Data_bigarray_val(oD), tb, alph, TMP );
    } else {
        compose_gtr( PA, Data_bigarray_val(oU), Data_bigarray_val(oD),
                     Data_bigarray_val(Some_val(oUi)), ta, alph, TMP );
        compose_gtr( PB, Data_bigarray_val(oU), Data_bigarray_val(oD),
                     Data_bigarray_val(Some_val(oUi)), tb, alph, TMP );
    }
    neg_log_comp( PA, alph, alph );
    neg_log_comp( PB, alph, alph );
    precalc( results.fmat, PA, PB );

    median_backtrace( a, b, m, &results );
    free(results.fmat->cost_asgn);
    free(results.fmat);
    CAMLreturn( Val_unit );
}
                        
value falign_CAML_median_wrapper( value* argv, int argn )
{
    return falign_CAML_median
        (argv[0], argv[1], argv[2], argv[3], argv[4], argv[5], argv[6],
            argv[7], argv[8], argv[9], argv[10]);
}


value
falign_CAML_backtrace( value oSpace, value oMat, value oU, value oD, value oUi,
                        value ota, value otb, value oa, value ob, value oea,
                         value oeb, value om, value ogap)
{
    CAMLparam5( oSpace, oMat, oU, oD, oUi );
    CAMLxparam5( ota, otb, oa, ob, oea );
    CAMLxparam3( oeb, om, ogap );

    seqt a,b, ea, eb, m;
    double ta, tb, *PA, *PB, *TMP;
    fcmt results;
    mat *space;
    int mat_size, alph, fcm_size, gap;

    Seq_custom_val(a,oa);
    Seq_custom_val(b,ob);
    Seq_custom_val(ea, oea);
    Seq_custom_val(eb, oeb);
    Seq_custom_val(m, om);
    space = FM_val( oSpace );
    results.direc = Matrices_struct( oMat );
    alph = Bigarray_val(oU)->dim[0];
    gap = Int_val( ogap );

    ta = Double_val( ota );
    tb = Double_val( otb );

    /* we need to ensure that we didn't "clear" the memory earlier; although we
     * never actually clear it, we need to at least ensure that an alignment is
     * stored. To re-register that section, since it must start at 0, 'free' and
     * register; this allows us to re-register space for composing matrices.
     *      TODO: can we be safe and re-use these cost matrices? */
    mat_size = a->len * b->len;
    assert( space->loc > mat_size );
    free_all( space );
    fcm_size = ((1 << alph)-1);

    results.fmat = (fm*) malloc( sizeof(fm) );
    results.fmat->alph = alph;
    results.fmat->size = fcm_size;
    results.fmat->comb = 1;
    results.fmat->gap  = 1 << gap;
    results.costs = register_section( space, mat_size, 0 );
    mat_setup_size (results.direc, a->len, b->len, 0, 0, alph, 0);
    results.fmat->cost = register_section( space, fcm_size*fcm_size, 0 );
    results.fmat->cost_asgn = (int*) malloc( sizeof(int)*fcm_size*fcm_size );

    /** pre-calculate the matrix **/
    PA = register_section( space, alph*alph, 0 );
    PB = register_section( space, alph*alph, 0 );
    TMP= register_section( space, alph*alph, 0 );
    if( Val_none == oUi ){
        compose_sym( PA, Data_bigarray_val(oU), Data_bigarray_val(oD), ta, alph, TMP );
        compose_sym( PB, Data_bigarray_val(oU), Data_bigarray_val(oD), tb, alph, TMP );
    } else {
        compose_gtr( PA, Data_bigarray_val(oU), Data_bigarray_val(oD),
                     Data_bigarray_val(Some_val(oUi)), ta, alph, TMP );
        compose_gtr( PB, Data_bigarray_val(oU), Data_bigarray_val(oD),
                     Data_bigarray_val(Some_val(oUi)), tb, alph, TMP );
    }
    neg_log_comp( PA, alph, alph );
    neg_log_comp( PB, alph, alph );
    precalc( results.fmat, PA, PB );

    full_backtrace( a, b, ea, eb, m, &results );
    free(results.fmat->cost_asgn);
    free(results.fmat);
    CAMLreturn( Val_unit );
}
 
value falign_CAML_backtrace_wrapper( value* argv, int argn )
{
    return falign_CAML_backtrace
        (argv[0], argv[1], argv[2], argv[3], argv[4], argv[5], argv[6],
            argv[7], argv[8], argv[9], argv[10], argv[11], argv[12]);
}



value
falign_CAML_align_2(value oSpace, value oMat, value oU, value oD, value oUi,
                        value ota, value otb, value oa, value ob, value ogap)
{

    CAMLparam5( oSpace, oMat, oU, oD, oUi );
    CAMLxparam5( ota, otb, oa, ob, ogap );
    CAMLlocal1( cost );

    double ta, tb, *PA, *PB, *TMP, min_cost;
    int alph, mat_size, fcm_size, gap;
    seqt a, b;
    fcmt results;
    mat *space;

    ta = Double_val( ota );
    tb = Double_val( otb );
    gap = Int_val( ogap );
    Seq_custom_val( a, oa );
    Seq_custom_val( b, ob );

    space = FM_val( oSpace );
    alph = Bigarray_val(oU)->dim[0];
    mat_size = a->len * b->len;
    fcm_size = ((1 << alph) - 1);

    /** register scratch space **/
    expand_matrix( space, (3 * alph * alph) + mat_size + fcm_size *fcm_size );
    results.direc = Matrices_struct( oMat );

    /** set up the cost matrix */
    results.fmat = (fm*) malloc( sizeof(fm) );
    results.fmat->alph = alph;
    results.fmat->size = fcm_size;
    results.fmat->comb = 1;
    results.fmat->gap = 1 << gap;
    mat_setup_size (results.direc, a->len, b->len, 0, 0, alph, 0);
    results.fmat->cost_asgn = (int*) malloc( sizeof(int)*fcm_size*fcm_size );
    mat_clean_direction_matrix( results.direc );

    results.costs = register_section( space, mat_size, 0 );
    PA = register_section( space, alph*alph, 0 );
    PB = register_section( space, alph*alph, 0 );
    TMP= register_section( space, alph*alph, 0 );
    results.fmat->cost = register_section( space, fcm_size*fcm_size, 0 );

    /** pre-calculate the matrix **/
    if( Val_none == oUi ){
        compose_sym( PA, Data_bigarray_val(oU), Data_bigarray_val(oD), ta, alph, TMP );
        compose_sym( PB, Data_bigarray_val(oU), Data_bigarray_val(oD), tb, alph, TMP );
    } else {
        compose_gtr( PA, Data_bigarray_val(oU), Data_bigarray_val(oD),
                     Data_bigarray_val(Some_val(oUi)), ta, alph, TMP );
        compose_gtr( PB, Data_bigarray_val(oU), Data_bigarray_val(oD),
                     Data_bigarray_val(Some_val(oUi)), tb, alph, TMP );
    }
    neg_log_comp( PA, alph, alph );
    neg_log_comp( PB, alph, alph );
    precalc( results.fmat, PA, PB );
    
    min_cost = full_falign( a, b, &results );
    cost = caml_copy_double( min_cost );

    print_alignment_matrixf( results.costs, results.direc, b->len, a->len );

    free(results.fmat->cost_asgn);
    free(results.fmat);
    CAMLreturn( cost );
}

value falign_CAML_align_2_wrapper( value* argv, int argn )
{
    return falign_CAML_align_2
        (argv[0],argv[1],argv[2],argv[3],argv[4],argv[5],argv[6],argv[7],argv[8],argv[9]);
}

value
falign_CAML_nukk_2(value oSpace, value oMat, value oU, value oD, value oUi,
                        value ota, value otb, value oa, value ob, value ogap)
{
    CAMLparam5( oSpace, oMat, oU, oD, oUi );
    CAMLxparam4( ota, otb, oa, ob );
    CAMLxparam1( ogap );

    CAMLreturn( caml_copy_double( 0.0 ) );
}

value falign_CAML_nukk_2_wrapper( value* argv, int argn )
{
    return falign_CAML_nukk_2
        (argv[0],argv[1],argv[2],argv[3],argv[4],argv[5],argv[6],argv[7],argv[8],argv[9]);
}

value
falign_CAML_cost_2(value oSpace, value oU, value oD, value oUi, value ota, value otb, value oa, value ob, value ogap)
{
    CAMLparam5( oSpace, oU, oD, oUi, ogap );
    CAMLxparam4( ota, otb, oa, ob );

    CAMLreturn( caml_copy_double( 0.0 ) );
}

value falign_CAML_cost_2_wrapper( value* argv, int argn )
{
    return falign_CAML_cost_2
        (argv[0],argv[1],argv[2],argv[3],argv[4],argv[5],argv[6],argv[7],argv[8]);
}

value
falign_CAML_closest(value oSpace, value oU, value oD, value oUi, value ot, value op, value om, value ores, value ogap)
{
    CAMLparam5( oSpace, oU, oD, oUi, ogap );
    CAMLxparam4( ot, op, om, ores );
    CAMLreturn( caml_copy_double(0.0) );
}

value falign_CAML_closest_wrapper( value* argv, int argn )
{
    return falign_CAML_closest
        (argv[0],argv[1],argv[2],argv[3],argv[4],argv[5],argv[6],argv[7],argv[8]);
}

#endif
