/* POY 5.1.1. A phylogenetic analysis program using Dynamic Homologies.       *\
(* Copyright (C) 2014 Andrés Varón, Lin Hong, Nicholas Lucaroni, Ward Wheeler,*)
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
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <time.h>

#include <caml/alloc.h>
#include <caml/mlvalues.h>
#include <caml/memory.h>
#include <caml/custom.h>
#include <caml/intext.h>
#include <caml/fail.h>

#include "floatmatrix.h"

#define CHECK_MEM(a); if(a==NULL || a==0){ printf("FM:failed on %d, ",__LINE__); failwith("I cannot allocate more memory");}
#define DEBUG 0

/* compare two float matrices */
int floatmatrix_CAML_compare( value one, value two )
{
    mat *a, *b;
    a = FM_val (one);
    b = FM_val (two);
    return (a->size - b->size);
}

/* free matrix */
void floatmatrix_CAML_free( value v )
{
    mat* s;
    s = FM_val (v);
    if( DEBUG )
        printf ("Freeing Data: %d --> 0\n",s->size);
    s = FM_val(v);
    free( s->mat );
    free( s );
}

/* serialize a matrix */
void floatmatrix_CAML_serialize(value v,unsigned long* wsize_32,unsigned long* wsize_64)
{
    mat* s; int i;
    s = FM_val (v);
    caml_serialize_int_4( s->size );
    for(i=0;i<s->size;++i)
        caml_serialize_float_8( s->mat[i] );
}

/* deserialize a matrix. allocate the size read in */
unsigned long floatmatrix_CAML_deserialize( void* v )
{
    mat* s; double* m; int i;
    s = (mat*) malloc( sizeof(mat) );    
    CHECK_MEM (s);
    s->size = caml_deserialize_sint_4();
    m = (double*) malloc( sizeof(double) * s->size );
    CHECK_MEM( m );
    for(i=0;i<s->size;++i)
        m[i] = caml_deserialize_float_8();
    return (sizeof(struct mat **));
}  

/* custom operations for float matrix */
static struct custom_operations floatmatrix_custom_operations  = {
    "http://www.amnh.org/poy/floatmatrix/floatmatrix.0.1", //identifier
    (&floatmatrix_CAML_free),        //finalize
    (&floatmatrix_CAML_compare),     //compare
    custom_hash_default,             //hash
    (&floatmatrix_CAML_serialize),   //serialize
    (&floatmatrix_CAML_deserialize)  //deserialize
};
/* register to the garbage collector */
value floatmatrix_CAML_register (value u)
{
    CAMLparam1( u );
    register_custom_operations (&floatmatrix_custom_operations);
    CAMLreturn (Val_unit);
}

/** --------------------------------- **
 ** start of non-gc related functions **
 ** --------------------------------- **/
 
/* resets the loc counter to 0 --doesn't actually free anything */
void free_all( mat *m )
{
    if (DEBUG)
        printf ("'Freeing' %d of %d data\n",m->loc,m->size);
    m->loc = 0; 
}

/* expand array to size [s] */
void expand_matrix (mat* m, const int s)
{
    double *tmp;
    free_all( m );
    if (m->size < s){
        if( DEBUG )
            printf ("Expanding Data: %d --> %d\n",m->size,s);
        if( m->loc > 0 ){
            printf("ERROR: Reallocation of memory while some is used\n");
            assert( 1 == 0 );
        }
        m->size = s;
        tmp = (double*) realloc(m->mat,sizeof(double)*s);
        CHECK_MEM( tmp );
        m->mat = tmp;
    }
}

/* clear a section of an array */
void clear_section( mat* m, int l, int h )
{
    if( DEBUG )
        printf ("Clearing Data: %d --> %d\n",l,h);
    assert ( (h-l) > 0 );
    memset( &(m->mat[l]), 0, sizeof(double)* (h-l) );
}

/* return how much free space we have */
int free_space( mat* m ) { return (m->size - m->loc); }

/* register a section of the matrix --clear section if asked */
double* register_section( mat* m, int s, int c )
{
    double *ptr;
    //all memory must be expanded explicitly
    if ( (m->loc + s) > m->size ){
        printf ("ERROR: not enough memory in current allocation: has %d@%d, needs %d\n",
                    m->size,m->loc,s);
        assert( 1 == 0 );
    }
    if ( c ) /* clear the section too */
        clear_section(m, m->loc, m->loc+s);
    if (DEBUG)
        printf ("Registering: %d -- %d of %d\n",m->loc,m->loc+s-1,m->size);
    ptr = &(m->mat[m->loc]);
    m->loc = m->loc + s;
    return ptr;
}

/* create a matrix with size of [s] */
value floatmatrix_CAML_create (value s)
{
    CAMLparam1( s );
    CAMLlocal1 ( ml_m );
    mat *m; int size;

    size = Int_val ( s );
    ml_m = caml_alloc_custom( &floatmatrix_custom_operations, sizeof(mat*), 1, 10000 );
    m = (mat*) malloc( sizeof(mat) );
    CHECK_MEM(m);
    FM_val ( ml_m ) = m;
    m->size = size;
    m->loc = 0;
    m->mat = (double*) malloc( s * sizeof(double) );
    CAMLreturn( ml_m );
}

void print_mat_all (mat* m)
{
    int i;
    for(i = 0;i<m->size;++i)
        printf("[%f] ", m->mat[i]);
    printf("\n");
}
void print_mat (mat* m, int off, int w, int h)
{
    int i,j, place;
    for(i = 0; i < w; ++i){
        place = i * w + h;
        for(j = 0; j < h; ++j){
            printf( "[%f] ", m->mat[place+j] );
        }
        printf("\n");
    }
    printf("\n");
}

/**
 *  Some functions used in testing the functionality.
*/
value floatmatrix_CAML_random (value m_val)
{
    CAMLparam1(m_val);
    mat *m; int i;
    srand( time(NULL) );
    m = FM_val (m_val);
    for(i=0;i<m->size;++i)
        m->mat[i] = ((double)rand()/(double)RAND_MAX);
    CAMLreturn (Val_unit);
}
value floatmatrix_CAML_print (value m_val)
{
    CAMLparam1( m_val );
    print_mat_all( FM_val( m_val ) );
    CAMLreturn( Val_unit );
}
value floatmatrix_CAML_clear (value m, value i,value j)
{
    CAMLparam3(m,i,j);
    clear_section( FM_val (m), Int_val (i), Int_val (j) );
    CAMLreturn( Val_unit );
}
value floatmatrix_CAML_expand(value m,value i)
{
    CAMLparam2(m,i);
    expand_matrix( FM_val(m), Int_val(i));
    CAMLreturn( Val_unit );
}
value floatmatrix_CAML_freeall(value m)
{
    CAMLparam1( m );
    free_all( FM_val(m) );
    CAMLreturn( Val_unit );
}

value floatmatrix_CAML_register_test(value m,value i,value v)
{
    CAMLparam3(m,i,v);
    mat *a;
    double *sec;
    int j;

    a = FM_val( m );
    sec = register_section (a, Int_val(i), 0);
    for(j = 0; j < Int_val(i); ++j)
        sec[j] = Double_val( v );
    CAMLreturn( Val_unit );
}

value floatmatrix_CAML_getsize (value m)
    { CAMLparam1(m); CAMLreturn( Int_val((FM_val(m))->size) );}
value floatmatrix_CAML_getused (value m)
    { CAMLparam1(m); CAMLreturn( Int_val((FM_val(m))->loc ) );}
