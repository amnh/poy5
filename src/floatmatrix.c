#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include <caml/alloc.h>
#include <caml/mlvalues.h>
#include <caml/memory.h>
#include <caml/custom.h>
#include <caml/intext.h>
#include <caml/fail.h>

#include "floatmatrix.h"

#define FM_val(v) (*((struct matrix**)Data_custom_val(v)))
#define FM_ptr(v)   ((struct matrix**)Data_custom_val(v))
#define CHECK_MEM(a) if(a==NULL) failwith("I cannot allocate more memory")
#define DEBUG 0

/* compare two float matrices, TODO: look inside array. I doubt this is really
 * necessary since it's just scratch space. */
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
 
/* expand array to size [s] */
void expand_matrix (mat* m, int s)
{
    if (m->size < s){
        if( DEBUG )
            printf ("Expanding Data: %d --> %d\n",m->size,s);
        m->size = s * 2;
        m->mat = (double*) realloc(m->mat,sizeof(double)*s*2);
        CHECK_MEM( m->mat );
    }
}

/* clear a section of an array */
void clear_section( mat* m, int l, int h )
{
    assert ( (h-l) > 0 );
    memset( &(m->mat[l]), 0, sizeof(double)* (h-l) );
}

/* register a section of the matrix --clear section if asked */
double* register_section( mat* m, int s, int c )
{
    double *ptr;
    if ( (m->loc + s) > m->size )
        expand_matrix (m, s + m->loc);
    if ( c ) /* clear the section too */
        clear_section(m, m->loc, m->loc+s);
    if (DEBUG)
        printf ("Registering: %d -- %d\n",m->loc,m->loc+s-1);
    ptr = &(m->mat[m->loc]);
    m->loc = m->loc + s;
    return ptr;
}

/* resets the loc counter to 0 --doesn't actually free anything */
void free_all( mat *m ){ m->loc = 0; }

/* create a matrix with size of [s] */
value floatmatrix_CAML_create (value s)
{
    CAMLparam1( s );
    CAMLlocal1 ( ml_m );
    mat *m; int size;

    size = Int_val ( s );
    ml_m = caml_alloc_custom( &floatmatrix_custom_operations, sizeof(mat*), 1, 10000 );
    m = (mat*) malloc( sizeof(mat) );
    FM_val ( ml_m ) = m;
    m->size = size;
    m->loc = 0;
    m->mat = (double*) malloc( s * sizeof(double) );
    CAMLreturn( ml_m );
}


/**
 *  Some functions used in testing the functionality.

value floatmatrix_CAML_random (value m_val) *
{
    CAMLparam1(m_val);
    mat *m; int i;
    srand(time(NULL) + getpid() );
    m = FM_val (m_val);
    for(i=0;i<m->size;++i)
        m->mat[i] = ((double)rand()/(double)RAND_MAX);
    CAMLreturn (Val_unit);
}
value floatmatrix_CAML_print (value m_val)
{
    CAMLparam1( m_val );
    mat *m; int i;
    m = FM_val( m_val );
    for(i = 0;i<m->size;++i)
        printf("[%f] ", m->mat[i]);
    printf("\n");
    CAMLreturn( Val_unit );
}
value floatmatrix_CAML_clear (value m, value i,value j)
{
    CAMLparam3(m,i,j);
    clear_section( FM_val (m), Int_val (i), Int_val (j) );
    CAMLreturn( Val_unit );
}
value get_size (value m){ CAMLparam1(m); CAMLreturn( Int_val((FM_val(m))->size) );}
value get_used (value m){ CAMLparam1(m); CAMLreturn( Int_val((FM_val(m))->loc ) );}
*/
