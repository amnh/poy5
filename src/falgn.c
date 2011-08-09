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
#ifndef INFINITY
    #define INFINITY    1/0
#endif

#define EPSILON         1e-6
#define MIN_UKK         1
#define MAX_UKK         INT_MAX

#define CLEAR           1

/** dor this file CDIR == DIRECTION_MATRIX in matiches.h **/

/** defines for ease of use **/
#define CDIR         DIRECTION_MATRIX
#define MAX(a,b)     (a>=b)?a:b
#define MIN(a,b)     (a<=b)?a:b
#define EQUALS(a,b)  fabs(a-b)<EPSILON
#define ISSET(x,i)   ((1<<i)&x)
#define GET(x,i)     seq_get(x,i)

#define DSET(f,p,d)  f->direc->matrix_d[p]=(CDIR)d
#define DGET(f,p)    f->direc->matrix_d[p]

/**
 * Conventions:
 *
 *  horizontal axis of the matrix has sequence labeled X, vertically Y. X is
 *  associated with the index, i, and j for Y. Insertions follow a veritical
 *  movement --a gap in sequence X, consuming sequence Y-- and alternatively
 *  Deletions horizontally and added to sequence Y.
 **/


void
print_alignment_matrixf(const double* costs, matricest mat, const int w, const int h)
{
    int i,j;
    for (i = 0; i < h; i++) {
        for (j = 0; j < w; j++)
            fprintf (stdout, "| %d %8.4f[00] ", mat->matrix_d[(w*i)+j],costs[(w*i)+j]);
        fprintf (stdout, "|\n");
    }
    fprintf (stdout, "\n");
    fflush( stdout );
}

void
print_alignment_indel_matrixf(const double* costs, matricest mat, const int* indels, const int w, const int h)
{
    int i,j;
    fprintf( stdout, "C:\n" );
    for (i = 0; i < h; i++) {
        for (j = 0; j < w; j++)
            fprintf (stdout, "| %d %8.4f[%2d] ",mat->matrix_d[(w*i)+j],costs[(w*i)+j],indels[(w*i)+j]);
        fprintf (stdout, "|\n");
    }
    fprintf (stdout, "\n");
    fflush( stdout );
}


void clear_matrixi( int *mat, const int x, const int y ){
    int i,j;

    for( i = 0; i < y; ++i ){
        for (j = 0; j < x; ++j ){
            mat[x*i + j] = 0;
        }
    }
}

void clear_matrixf( double *mat, const int x, const int y ){
    int i,j;

    for( i = 0; i < y; ++i ){
        for (j = 0; j < x; ++j ){
            mat[x*i + j] = 0;
        }
    }
}


/** Create a full backtrace with median and edited sequences; polymorphic **/
void
median_backtrace (const seqt x, const seqt y, seqt m, fcmt *FA)
{
    int x_loc, y_loc, linear;
    x_loc = x->len - 1;
    y_loc = y->len - 1;

    seq_clear( m );
    while( !(x_loc == 0 && y_loc == 0) )
    {
        linear = y_loc * x->len + x_loc;
        if( DELETE & DGET(FA,linear) ){
            /*printf("(%d,%d):%d=%d -I- %d\n", x_loc,y_loc,linear,DGET(FA,linear), 
                                     fasn(FA->fmat,GET(x,x_loc),FA->fmat->gap));*/
            seq_prepend( m , fasn(FA->fmat,FA->fmat->gap,GET(y,y_loc)));
            --y_loc;

        } else if( ALIGN & DGET(FA,linear) ){
            /*printf("(%d,%d):%d=%d -A- %d\n", x_loc,y_loc,linear,DGET(FA,linear), 
                                        fasn(FA->fmat,GET(x,x_loc),GET(y,y_loc)));*/
            seq_prepend( m , fasn(FA->fmat,GET(x,x_loc),GET(y,y_loc)));
            --x_loc;
            --y_loc;

        } else if( INSERT & DGET(FA,linear) ){
            /*printf("(%d,%d):%d=%d -D- %d\n", x_loc,y_loc,linear,DGET(FA,linear), 
                                              fasn(FA->fmat,FA->fmat->gap,GET(y,y_loc)));*/
            seq_prepend( m , fasn(FA->fmat,GET(x,x_loc),FA->fmat->gap));
            --x_loc;

        } else {
            printf("(%d,%d):%d=%d - Error!\n", x_loc, y_loc, linear, DGET(FA,linear) );
            failwith( "Unknown Direction" );
        }
        fflush(stdout);
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
        linear = y_loc * x->len + x_loc;
        if( DELETE & DGET(FA,linear) ){
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

        } else if( INSERT & DGET(FA,linear) ){
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
general_cost(const seqt x, const seqt y, fcmt *FA, const int i, const int j)
{
    double cost_aln, cost_ins, cost_del, temp_min;
    int place;
    CDIR temp_dir;
    place = (x->len * j) + i;
    if( 0 == i && 0 == j){
        FA->costs[ place ] = 0.0;
    } else if( 0 == i ){
        FA->costs[ place ] = FA->costs[ (j-1)*x->len ] + fcost( FA->fmat, FA->fmat->gap, GET(y,j) );
        DSET( FA, place, DELETE );
    } else if( 0 == j ){
        FA->costs[ place ] = FA->costs[ i-1 ] + fcost( FA->fmat, GET(x,i), FA->fmat->gap );
        DSET( FA, place, INSERT );
    } else {
        cost_aln = FA->costs[(j-1)*x->len + i-1] + fcost( FA->fmat, GET(x,i),      GET(y,j)      );
        cost_ins = FA->costs[ j   *x->len + i-1] + fcost( FA->fmat, GET(x,i),      FA->fmat->gap );
        cost_del = FA->costs[(j-1)*x->len + i  ] + fcost( FA->fmat, FA->fmat->gap, GET(y,j)      );

        if( EQUALS(cost_aln,cost_ins) ){
            temp_dir = ALIGN | INSERT;
            temp_min = MIN(cost_aln,cost_ins);
        } else if ( cost_aln < cost_ins ){
            temp_dir = ALIGN;
            temp_min = cost_aln;
        } else {
            temp_dir = INSERT;
            temp_min = cost_ins;
        }
        if( EQUALS( temp_min, cost_del ) ){
            temp_dir |=  DELETE;
            temp_min  = MIN(temp_min,cost_del);
        } else if ( cost_del < temp_min ){
            temp_dir = DELETE;
            temp_min = cost_del;
        }

        FA->costs[ place ] = temp_min;
        DSET( FA, place, temp_dir );
        /*printf("(%d,%d)[%d]: A:%f\tI:%f\tD:%f --> %f[%d]\n",j,i,place,cost_aln,cost_ins,cost_del,temp_min,temp_dir);*/
    }
}


/** Do a full alignment; filling the entire cost matrix; FCM is allocated **/
double full_falign( const seqt x, const seqt y, fcmt *FA )
{
    int i, j;
    for( j = 0; j < y->len; ++j ){
        for( i = 0; i < x->len; ++i ){
            general_cost( x, y, FA, i, j );
        }
    }
    return FA->costs[ x->len*y->len - 1 ];
}


/** exclude i-1 **/
void update_exclude_row( fcmt *FA, const seqt x, const seqt y, const int i, const int j )
{
    int place;
    double cost_aln, cost_del;

    place = (x->len * j) + i;
    cost_aln = FA->costs[(j-1)*x->len + i-1] + fcost( FA->fmat, GET(x,i),      GET(y,j));
    cost_del = FA->costs[(j-1)*x->len + i  ] + fcost( FA->fmat, FA->fmat->gap, GET(y,j));

    if( EQUALS(cost_aln, cost_del) ){
        DSET( FA, place, ALIGN | DELETE );
        FA->costs[place] = MIN( cost_aln, cost_del );
        FA->nukk[place] = MIN( 1 + FA->nukk[(j-1)*x->len+i], 
                               ((GET(x,i)&GET(y,j))>0?0:1) + FA->nukk[(j-1)*x->len+(i-1)]);
    } else if ( cost_aln < cost_del ){
        DSET( FA, place, ALIGN );
        FA->costs[place] = cost_aln;
        FA->nukk[place] = ((GET(x,i)&GET(y,j))>0?0:1) + FA->nukk[(j-1)*x->len+(i-1)];
    } else {
        DSET( FA, place, DELETE );
        FA->costs[place] = cost_del;
        FA->nukk[place] = 1 + FA->nukk[(j-1)*x->len+i];
    }
/*    printf("\tC:(%d,%d): A:%f\tI:%f\tD:%f --> %f[%d]\n",*/
/*            i,j,cost_aln,0.0,cost_del,FA->costs[place],DGET(FA,place) );*/
}


/** exclude j-1 **/
void update_exclude_col( fcmt *FA, const seqt x, const seqt y, const int i, const int j )
{
    int place;
    double cost_aln, cost_ins;

    place = (x->len * j) + i;
    cost_aln = FA->costs[(j-1)*x->len + i-1] + fcost( FA->fmat, GET(x,i), GET(y,j)      );
    cost_ins = FA->costs[ j   *x->len + i-1] + fcost( FA->fmat, GET(x,i), FA->fmat->gap );

    if( EQUALS(cost_aln, cost_ins) ){
        DSET( FA, place, ALIGN | INSERT );
        FA->costs[place] = MIN( cost_aln, cost_ins );
        FA->nukk[place] = MIN( 1 + FA->nukk[j*x->len + (i-1)],
                               ((GET(x,i)&GET(y,j))>0?0:1) + FA->nukk[(j-1)*x->len+(i-1)]);
    } else if ( cost_aln < cost_ins ){
        DSET( FA, place, ALIGN );
        FA->costs[place] = cost_aln;
        FA->nukk[place] = ((GET(x,i)&GET(y,j))>0?0:1) + FA->nukk[(j-1)*x->len+(i-1)];
    } else {
        DSET( FA, place, INSERT );
        FA->costs[place] = cost_ins;
        FA->nukk[place] = 1 + FA->nukk[j*x->len + (i-1)]; 
    }
/*    printf("\tC:(%d,%d): A:%f\tI:%f\tD:%f --> %f[%d]\n",*/
/*            i,j,cost_aln,cost_ins,0.0,FA->costs[place],DGET(FA,place) );*/
}


void update_all( fcmt *FA, const seqt x, const seqt y, const int i, const int j )
{
    int place, temp_gap;
    double cost_aln, cost_ins, cost_del, temp_min;
    CDIR temp_dir;

/*    printf("\tC:A:%f\tI:%f\tD:%f\n", FA->costs[(j-1)*x->len + i-1],*/
/*           FA->costs[j*x->len + i-1], FA->costs[(j-1)*x->len + i]);*/

    place = j * x->len + i;
    cost_aln = FA->costs[(j-1)*x->len + i-1] + fcost( FA->fmat, GET(x,i),      GET(y,j)      );
    cost_ins = FA->costs[ j   *x->len + i-1] + fcost( FA->fmat, GET(x,i),      FA->fmat->gap );
    cost_del = FA->costs[(j-1)*x->len + i  ] + fcost( FA->fmat, FA->fmat->gap, GET(y,j)      );

    if( EQUALS(cost_aln,cost_ins) ){
        temp_dir = ALIGN | INSERT;
        temp_min = MIN(cost_aln,cost_ins);
        temp_gap =  MIN( 1 + FA->nukk[j*x->len + (i-1)],
                         ((GET(x,i)&GET(y,j))>0?0:1) + FA->nukk[(j-1)*x->len+(i-1)]);
    } else if ( cost_aln < cost_ins ){
        temp_dir = ALIGN;
        temp_min = cost_aln;
        temp_gap = ((GET(x,i)&GET(y,j))>0?0:1) + FA->nukk[(j-1)*x->len+(i-1)];
    } else {
        temp_dir = INSERT;
        temp_min = cost_ins;
        temp_gap =  1 + FA->nukk[j*x->len + (i-1)];
    }
    if( EQUALS( temp_min, cost_del ) ){
        temp_dir |=  DELETE;
        temp_min  = MIN(temp_min,cost_del);
        temp_gap  = MIN( 1 + FA->nukk[(j-1)*x->len+i], temp_gap );
    } else if ( cost_del < temp_min ){
        temp_dir = DELETE;
        temp_min = cost_del;
        temp_gap = 1 + FA->nukk[(j-1)*x->len+i];
    }
/*    printf("\tC:(%d,%d): A:%f\tI:%f\tD:%f --> %f[%d]\n",*/
/*            i,j, cost_aln, cost_ins, cost_del, temp_min, temp_dir);*/

    FA->nukk[place]  = temp_gap;
    FA->costs[place] = temp_min;
    DSET( FA, place, temp_dir );
}


void build_strip( const seqt x, const seqt y, fcmt *FA, const int k )
{
    int i, b, j, i_min, i_max, p_max;
    b = (k - (x->len - y->len)) / 2;

    FA->costs[ 0 ] = 0.0;
    FA->nukk[ 0 ]  = 0;
    for( i=1; i < x->len; ++i ){ /** update row 0 **/
        FA->costs[ i ] = FA->costs[ i-1 ] + fcost( FA->fmat, GET(x,i), FA->fmat->gap );
        DSET( FA, i, INSERT );
        FA->nukk[i] = i;
    }
    for( j=1; j < y->len; ++j ){ /** update col 0 **/
        FA->costs[ j*x->len ] = FA->costs[ (j-1)*x->len ] + fcost( FA->fmat, FA->fmat->gap, GET(y,j) );
        DSET( FA, j*x->len, DELETE );
        FA->nukk[j*x->len] = j;
    }

    p_max = 1;
    for( j=1; j < y->len; ++j ) {
        i_min = MAX( 1, j-b-1 );
        i_max = MIN( x->len-1, j + b + (x->len - y->len) );
        /*printf("C:Building Strip: (%d,%d) -> (%d,%d)\n",i_min,j,i_max,j);*/
        if( i_min == 1 ){
            update_all( FA, x, y, i_min, j );
        } else {
            update_exclude_row( FA, x, y, i_min, j );
        }
        for( i = i_min+1; i < i_max; ++i ){
            update_all( FA, x, y, i, j );
        }
        if( p_max == (x->len-1) ){
            update_all( FA, x, y, i_max, j );
        } else {
            update_exclude_col( FA, x, y, i_max, j );
        }
        p_max = i_max;
    }
}


void update_barrier( const seqt x, const seqt y, fcmt *FA, const int k, const int ok )
{
    int i_min, i_max, i, b, ob, j, p_max;

    b  =  (k - (y->len - x->len)) / 2;
    ob = (ok - (y->len - x->len)) / 2;

    p_max = 1;
    for( j=1; j < y->len; ++j ) {
        i_min = MAX( 1, j-b-1 );
        i_max = MIN( x->len-1, j + b + (x->len - y->len) );
        /*printf("C:Updating Strip: (%d,%d) -> (%d,%d)\n",i_min,j,i_max,j);*/
        if( i_min == 1 ){
            update_all( FA, x, y, i_min, j );
        } else {
            update_exclude_row( FA, x, y, i_min, j );
        }
        for( i = i_min+1; i < i_max; ++i ){
            update_all( FA, x, y, i, j );
        }
        if( p_max == (x->len-1) ){
            update_all( FA, x, y, i_max, j );
        } else {
            update_exclude_col( FA, x, y, i_max, j );
        }
        p_max = i_max;
    }
}


double nukk_falign( const seqt x, const seqt y, fcmt *FA )
{
    int k;

    assert( x->len >= y->len );
    k = MAX( MIN_UKK, (x->len-y->len)+1 );
    if( CLEAR ){
        clear_matrixi( FA->nukk, x->len, y->len );
        clear_matrixf( FA->costs, x->len, y->len );
        mat_clean_direction_matrix( FA->direc );
    }
    build_strip( x, y, FA, k );
    /*printf("C:Finished Building Strip k = %d\n", k);*/

    while( k <= FA->nukk[ x->len*y->len - 1 ] && (k < MAX_UKK) ){
        /*printf("C:Increasing K to %d\n", k*2);*/
        update_barrier( x, y, FA, k*2, k );
        k = k*2;
    }
    /*printf("C:Finished Alignment in k:%d = %f\n", k, FA->costs[ x->len*y->len - 1 ]);*/
    return FA->costs[ x->len*y->len - 1 ];
}

value
falign_CAML_median(value oSpace, value oMat, value oU, value oD, value oUi,
                   value ota, value otb, value oa, value ob, value om, value ogap)
{
    CAMLparam5( oSpace, oMat, oU, oD, oUi );
    CAMLxparam5( ota, otb, oa, ob, om );
    CAMLxparam1( ogap );

    seqt a,b, m;
    double ta, tb, *PA, *PB;
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
    fcm_size = ((1 << alph)-1);
    assert( space->loc == ((3 * alph * alph) + mat_size + fcm_size *fcm_size) );

    results.fmat = (fm*) malloc( sizeof(fm) );
    results.fmat->alph = alph;
    results.fmat->size = fcm_size;
    results.fmat->comb = 1;
    results.fmat->gap  = 1 << gap;
    mat_setup_size (results.direc, a->len, b->len, 0, 0, alph, 0);
    results.fmat->cost_asgn = (int*) malloc( sizeof(int)*fcm_size*fcm_size );

    /** re-register the matrix **/
    /** CARE MUST BE TAKEN THAT THESE ARE REGISTERED THE SAME WAY AS DOWNPASS! **/
    free_all( space ); /** reset location; no need to expand **/
    results.costs      = register_section( space, mat_size, 0 );
    results.fmat->cost = register_section( space, fcm_size*fcm_size, 0 );
    PA                 = register_section( space, alph*alph, 0 );
    PB                 = register_section( space, alph*alph, 0 );
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
falign_CAML_print_alignments( value oSpace, value oMat, value oAlen, value oBlen )
{
    CAMLparam4( oSpace, oMat, oAlen, oBlen );
    int aLen, bLen;
    mat      *align;
    matricest direc;

    aLen = Int_val( oAlen );
    bLen = Int_val( oBlen );
    direc = Matrices_struct( oMat );
    align = FM_val( oSpace );

    /* space must start with the alignment matrix; as is usual */
    print_alignment_matrixf( align->mat, direc, aLen, bLen );

    CAMLreturn( Val_unit );
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
    double ta, tb, *PA, *PB;
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
    fcm_size = ((1 << alph)-1);

    results.fmat = (fm*) malloc( sizeof(fm) );
    results.fmat->alph = alph;
    results.fmat->size = fcm_size;
    results.fmat->comb = 1;
    results.fmat->gap  = 1 << gap;
    mat_setup_size (results.direc, a->len, b->len, 0, 0, alph, 0);
    results.fmat->cost_asgn = (int*) malloc( sizeof(int)*fcm_size*fcm_size );

    /** re-register the matrix **/
    /** CARE MUST BE TAKEN THAT THESE ARE REGISTERED THE SAME WAY AS DOWNPASS! **/
    free_all( space ); /** reset location; no need to expand **/
    results.costs      = register_section( space, mat_size, 0 );
    results.fmat->cost = register_section( space, fcm_size*fcm_size, 0 );
    PA                 = register_section( space, alph*alph, 0 );
    PB                 = register_section( space, alph*alph, 0 );
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

    results.costs      = register_section( space, mat_size, 0 );
    results.fmat->cost = register_section( space, fcm_size*fcm_size, 0 );
    PA                 = register_section( space, alph*alph, 0 );
    PB                 = register_section( space, alph*alph, 0 );
    TMP                = register_section( space, alph*alph, 0 );

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

    /*print_alignment_matrixf( results.costs, results.direc, a->len, b->len );*/

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

    /** register scratch space; free first to 'erase' information. **/
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
    results.nukk = (int*) malloc( sizeof(int)*mat_size );
    mat_clean_direction_matrix( results.direc );

    results.costs      = register_section( space, mat_size, 0 );
    results.fmat->cost = register_section( space, fcm_size*fcm_size, 0 );
    PA                 = register_section( space, alph*alph, 0 );
    PB                 = register_section( space, alph*alph, 0 );
    TMP                = register_section( space, alph*alph, 0 );

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
    
    min_cost = nukk_falign( a, b, &results );
    cost = caml_copy_double( min_cost );

    /*print_alignment_indel_matrixf ( results.costs, results.direc, results.nukk, a->len, b->len );*/

    free(results.fmat->cost_asgn);
    free(results.fmat);
    free(results.nukk);

    CAMLreturn( cost );
}

value falign_CAML_nukk_2_wrapper( value* argv, int argn )
{
    return falign_CAML_nukk_2
        (argv[0],argv[1],argv[2],argv[3],argv[4],argv[5],argv[6],argv[7],argv[8],argv[9]);
}


value
falign_CAML_cost_2(value oSpace, value oU, value oD, value oUi, value ot, value oa, value ob, value ogap)
{
    CAMLparam5( oSpace, oU, oD, oUi, ogap );
    CAMLxparam3( ot, oa, ob );
    CAMLlocal1( ocosts );

    double t, *P, *TMP, costs, min_cost;
    int alph, gap, i, aj, bj;
    seqt a, b;
    mat *space;

    t = Double_val( ot );
    gap = Int_val( ogap );
    Seq_custom_val( a, oa );
    Seq_custom_val( b, ob );
    assert( a->len == b->len );

    space = FM_val( oSpace );
    alph = Bigarray_val(oU)->dim[0];

    /** register scratch space **/
    expand_matrix( space, (2 * alph * alph) );
    P   = register_section( space, alph*alph, 0 );
    TMP = register_section( space, alph*alph, 0 );

    /** pre-calculate the matrix **/
    if( Val_none == oUi ){
        compose_sym( P, Data_bigarray_val(oU), Data_bigarray_val(oD), t, alph, TMP );
    } else {
        compose_gtr( P, Data_bigarray_val(oU), Data_bigarray_val(oD),
                     Data_bigarray_val(Some_val(oUi)), t, alph, TMP );
    }
    neg_log_comp( P, alph, alph );

    costs = 0.0;
    for( i = 1; i < a->len; ++i ){
        min_cost = INFINITY;
        for( aj = 0; aj < alph; ++aj ){
            if( ISSET(GET(a,i),aj) ){
                for( bj = 0; bj < alph; ++bj ){
                    if( ISSET(GET(b,i),bj) ){
                        min_cost = MIN( min_cost, P[ aj*alph + bj ]);
                    }
                }
            }
        }
        costs += min_cost;
    }
    ocosts = caml_copy_double( costs );
    CAMLreturn( ocosts );
}

value falign_CAML_cost_2_wrapper( value* argv, int argn )
{
    return falign_CAML_cost_2
        (argv[0],argv[1],argv[2],argv[3],argv[4],argv[5],argv[6],argv[7]);
}


double
closest( seqt n_child, const seqt child, const seqt parent, const double* cm, const int alph, const int gap )
{
    int i, g, tmp;
    double cost, s_cost;

    assert( child->len == parent->len );
    seq_clear( n_child );
    g = (1 << gap);

    /* end to begining to utilize seq prepend */
    cost = 0.0;
    for( i = child->len-1; i > 0; --i ){
        if( g == GET(child,i) || g == GET(parent,i) ){
            tmp = g;
        } else {
            tmp = (~g) & GET(child, i);
        }
        tmp = get_closest(&s_cost, tmp, GET(parent,i), cm, alph );
        cost += s_cost;
        /* printf( "C:%d\t P:%d M:%d --> %d : %f\n",i, GET(parent,i), GET(child,i), tmp, s_cost ); */
        seq_prepend( n_child, tmp );
    }
    seq_prepend( n_child, gap );
    return cost;
}


value
falign_CAML_closest(value oSpace, value oU, value oD, value oUi, value ot, value op, value om, value on, value ogap)
{
    CAMLparam5( oSpace, oU, oD, oUi, ogap );
    CAMLxparam4( ot, op, om, on );
    CAMLlocal1( cost );

    double t, *P, *T;
    seqt p, m, n;
    int alph, gap;
    mat *space;

    space = FM_val( oSpace );
    t = Double_val( ot );
    alph = Bigarray_val(oU)->dim[0];
    gap = Int_val( ogap );

    Seq_custom_val( p, op );
    Seq_custom_val( m, om );
    Seq_custom_val( n, on );

    expand_matrix( space, (2*alph*alph) );
    P = register_section( space, alph*alph, 0 );
    T = register_section( space, alph*alph, 0 );

    if( Val_none == oUi ){
        compose_sym( P, Data_bigarray_val(oU), Data_bigarray_val(oD), t, alph, T );
    } else {
        compose_gtr( P, Data_bigarray_val(oU), Data_bigarray_val(oD),
                     Data_bigarray_val(Some_val(oUi)), t, alph, T );
    }
    neg_log_comp( P, alph, alph );

    cost = caml_copy_double( closest( n, p, m, P, alph, gap ) );
    CAMLreturn( cost );
}

value falign_CAML_closest_wrapper( value* argv, int argn )
{
    return falign_CAML_closest
        (argv[0],argv[1],argv[2],argv[3],argv[4],argv[5],argv[6],argv[7],argv[8]);
}

#endif
