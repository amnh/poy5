/* $Id: lists.c 48 2005-04-13 15:56:25Z ron $
   Written by Adam Siepel, Spring 2001
   Copyright 2001, Adam Siepel */

/* Simple list-handling functions, allowing indexing or either FIFO or
   LIFO behavior. */

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "lists.h"
#include <caml/fail.h>
//check if the list is full, full return 1, else return 0
int
is_full ( List * q )
{
     if (( q->ridx >= q->CAPACITY ) && ( q->lidx <= 0 ))
         return 1;
     else
         return 0;
}


int
is_empty (List * q)
{
     if ( q->lidx >= q->ridx ) return 1;
     else return 0;
}

void
push ( List * q,  ElementUnion v )
{
    int i=0;
    if ( q->ridx >= q->CAPACITY )
    {
        if ( q->lidx > 0 )
        {
            for ( i = q->lidx; i < q->ridx; i++ )
                q->array[i - q->lidx] = q->array[i];
            q->ridx -= q->lidx;
            q->lidx = 0;
        }
        else
        {
            fprintf ( stderr, "ERROR: Exceeded list capacity,list's right idx=%d, left idx=%d.\n" );
            assert ( 0 );
        }
    }
    q->array[q->ridx++] = v;

}

int
list_size ( List * l )
{
    return ( l->ridx - l->lidx );
}

//void *
ElementUnion 
list_get ( List * l, int i )
{
    if ( i >= list_size ( l ) )
        //return NULL;
    { 
        fprintf(stderr,"ERROR: Try to get list[%d], but the list only has %d elements\n", i, list_size(l));
        assert(0);
    }
      return (l->array[l->lidx + i]);
    //return res;
}

//void *
ElementUnion
pop_queue ( List * q )
{
    if ( q->lidx >= q->ridx )
    //    return NULL;
     { 
        fprintf(stderr,"ERROR: Try to pop an empty list\n ");
        assert(0);
     }
    return (q->array[q->lidx++]);
 }

//void *
ElementUnion
peek_queue ( List * q )
{
    if ( q->lidx >= q->ridx )
      //  return NULL;
      { 
        fprintf(stderr,"ERROR: Try to peek an empty list\n ");
        assert(0);
     }
     return q->array[q->lidx];
}

//void *
ElementUnion
pop_stack ( List * q )
{
    if ( q->ridx <= q->lidx )
    //    return NULL;
    { 
        fprintf(stderr,"ERROR: Try to pop an empty stack\n ");
        assert(0);
     }
    return q->array[--q->ridx];
}

//void * 
ElementUnion
peek_stack ( List * s )
{
    if ( s->ridx <= s->lidx )
   //     return NULL;
    { 
        fprintf(stderr,"ERROR: Try to peek an empty stack\n ");
        assert(0);
     }
   return s->array[s->ridx-1];
}

int
empty ( List * q )
{
    return ( q->lidx >= q->ridx );
}

/* Must be executed on a new list before it is usable!!! */
void
init_list ( List * q, int nelements,  enum datatype dtype )
{
    q->dtype = dtype;
    q->ridx = q->lidx = 0;
    q->CAPACITY = nelements;
    q->array = ( ElementUnion * ) malloc ( nelements * sizeof(ElementUnion) );
    if ( q->array == NULL)
        fprintf(stderr, "cannot malloc list (funciton init_list in grappa/lists.c) \n");
}

void
copy_list ( List * new, List * old )
{
    int i;
    init_list ( new, old->CAPACITY, old->dtype );
    for ( i = 0; i < list_size ( old ); i++ )
    {
        push ( new, list_get ( old, i ) );
    }
}

void
free_list ( List * q )
{
    free(q->array);
}

void
clear_list ( List * l )
{
    l->ridx = l->lidx = 0;
}

void
list_delete ( List * l, int idx )
{
    int i;
    if ( idx >= list_size ( l ) )
        return;
     (l->array)[l->lidx + i - 1] = (l->array)[l->lidx + i];
    l->ridx--;
}

/* note: pointer-based */
int
//list_contains ( List * l, void *ptr )
list_contains ( List * l, ElementUnion* ptr )
{
    int i;
    ElementUnion * res; List * browse_lst;
    switch ( l->dtype)
    {
        case IntData:
            for ( i = 0; i < list_size ( l ); i++ )
            { 
                  if (list_get ( l, i ).intelement == ptr->intelement)
                      return 1;
            }
            break;
        case ReversalData:
            for ( i = 0; i < list_size ( l ); i++ )
            { 
                  if(
                      (list_get ( l, i ).revelement.start == ptr->revelement.start)
                      &&
                      (list_get ( l, i ).revelement.stop == ptr->revelement.stop)
                    )
                      return 1;
            }
            break;
        case ListData: //on what condition ppl consider two lists are the same?
            fprintf(stderr,"ERROR: we don't compare two list now\n");
            assert(0);
   /*         for ( i = 0; i < list_size ( l ); i++ )
            { 
                browse_lst = list_get ( l, i ).listelement;
              //  input = ptr->listelement;
                int j; int sign = 1;
                if( list_size(browse_lst) == list_size(ptr->listelement) )
                {
                    for (j=0; j<list_size(ptr->listelement); j++)
                    {
                        if ( list_contains(browse_list,list_get(ptr->listelement,j)))
                        {  sign =1; }
                        else { sign = 0; }
                    }
                    if (1==sign)  return 1;
                }
            }*/
            break;
        case VertexData: //on what condition ppl consider two lists are the same?
            fprintf(stderr,"ERROR: we don't compare two vertex now\n");
            assert(0);
            break;
        case IntarrayData:
            fprintf(stderr,"ERROR: we don't compare two int array now\n");
            assert(0);
            break;
        default:
             fprintf(stderr,"ERROR: try to identify an element with undefined data type\n");
             assert(0);
             break;

    }
    return 0;
}
