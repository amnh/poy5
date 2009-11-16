/* $Id: lists.h 48 2005-04-13 15:56:25Z ron $
   Written by Adam Siepel, Spring 2001
   Copyright 2001, Adam Siepel */

/* Simple list-handling functions, allowing indexing or either FIFO or
   LIFO behavior. */

#ifndef LISTS_H
#define LISTS_H

#include "vertex_factory.h"

enum datatype { IntData = 1, ReversalData=2, ListData=3, VertexData=4, IntarrayData = 5} ;

typedef struct reversal_struct Reversal;

typedef union elementunion ElementUnion;


typedef struct list List;
struct list
{
    enum datatype dtype;
    //void **array;
   // ArrayUnion array;
    ElementUnion * array;
    int lidx;
    int ridx;
    int CAPACITY;
   // int elementsz;
};



struct reversal_struct
{
    int start;
    int stop;
};

union elementunion
{
   // enum datatype dtype;
    int intelement;
    Reversal  revelement;
    List listelement;
    Vertex vertexelement;
    int * intarrelement;
    //int * intlistelement;
} ;


int check ( List * q );
int is_empty ( List * q );
//void init_list ( List * q, int nelements, int elementsz );
void init_list ( List * q, int nelements, enum datatype dtype );
void free_list ( List * q );
//void push ( List * q, void *v );
void push ( List * q, ElementUnion v );
//void *
ElementUnion pop_stack ( List * q );
//void *
ElementUnion pop_queue ( List * q );
//void *
ElementUnion peek_queue ( List * q );
//void *
ElementUnion peek_stack ( List * q );
int empty ( List * q );
int list_size ( List * l );
void copy_list ( List * old, List * new );
//void *
ElementUnion list_get ( List * l, int i );
void clear_list ( List * l );
void list_delete ( List * l, int idx );
//int list_contains ( List * l, void *ptr );
int list_contains ( List * l, ElementUnion* ptr );

#endif
