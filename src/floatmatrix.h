/* POY 5.0 Beta. A phylogenetic analysis program using Dynamic Homologies.    *\
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

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include <caml/mlvalues.h>
#include <caml/memory.h>
#include <caml/custom.h>
#include <caml/intext.h>

#define FM_val(v) (*((struct matrix**)Data_custom_val(v)))
#define FM_ptr(v)   ((struct matrix**)Data_custom_val(v))

struct matrix {
    int size;       /* current maximum size of array */
    int loc;        /* current portion used */
    double *mat;    /* data */
};
typedef struct matrix mat;

void floatmatrix_CAML_free_floatmatrix( value v );
int  floatmatrix_CAML_compare( value one,value two );
void floatmatrix_CAML_serialize(value v, unsigned long* wsize_32, unsigned long* wsize_64);
unsigned long floatmatrix_CAML_deserialize( void* dst );
value floatmatrix_CAML_register (value u);

/* registers section and returns ptr */
double* register_section (mat* m,int s,int c);
void expand_matrix( mat* m, int s );
int free_space (mat *m);
void free_all (mat* m);
/* void clear_section( mat* m, int l, int h ); -- register should do this */
