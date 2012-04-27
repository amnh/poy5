/* POY 5.0 Alpha. A phylogenetic analysis program using Dynamic Homologies.   *\
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


/*
 newkkonen.h and newkkonen.c are alignment c functions with our version of Ukkonen.
 Ukkonen algorithm is based on ukkonen's paper "Algorithms for Approximate String Matching".
 instead of updating the whole matrix, we only update a band of 
 diagonals in each iteration. 
 algn.c has functions for ukkonen alignment, too.
 since we only need part of the matrix, we only need to allocate that part of memory. 
 that's difference between this and those in algn.c. 
 so far, this is still ukkonen. but we made a twist here. 
 there is a threshold T in ukkonen, when we get a cost C<T, we stop, 
 otherwise, we double T, continue updating more diagonals. 
 T is set to diff_len * delta. 
 diff_len is the difference of length of two sequences. delta is min(substitution cost, gap cost).
 we are using a new threshold K here. K is the possible number of gaps between two sequence. 
 K = num_gaps + num_diff_char. 
 num_diff_char is the number of different charactors between two sequences during this iteration. 
 since we only update a band of diagonals, say bandsize = B. when B/2 > K , 
 updating more diagonals won't give us better cost.
 NOTE: K will change with iterations, while T remains the same.
 */
#include "seq.h"
#include "cm.h"

#define Newkkmat_struct(a) ((struct newkkmat *) Data_custom_val(a))
#define DIRECTION_MATRIX unsigned short

#ifdef USE_LONG_SEQUENCES
#define MAT_SIZE long int
#else 
#define MAT_SIZE int
#endif
//each cost_dir is a diagonal of ukkmatrix.
struct cost_dir {
    int len; //len of this diagonal,start from 1
    int * costarr; //cost 
    DIRECTION_MATRIX * dirarr;//direction
    DIRECTION_MATRIX * gapnumarr1;//max gap number we could have in alignment
    DIRECTION_MATRIX * gapnumarr2;//max gap number we could have in alignment
    int * affParr; //affine matrix P, vertical
    int * affQarr; //affine matrix Q, horizontal
    int * affDarr; //affine matrix Diagonal
};

typedef struct cost_dir * ukkdiag_p;
struct newkkmat {
    MAT_SIZE baseband;  // this doesn't change. we start from (lenY-lenX+1) diagonals. 
    //we define baseband as MAT_SIZE here just because total_len are the result
    //of baseband * lenX. if baseband is int, even with "use_long_sequence", we still
    //won't get a long int result of total_len.
    MAT_SIZE total_len;          /* how many ukk cells are there */
    MAT_SIZE total_len_in_use;     //how many cells are in use.
    MAT_SIZE k;                  //current k.define it MAT_SIZE like baseband.
    MAT_SIZE diag_size_in_use;       //just how many diagonal are there in use. start from 1
    MAT_SIZE diag_size; //size of array diagonal_size_arr allocated.
    int * diagonal_size_arr;    
    // size array of diagonal, size = baseband+K*2
    // here is the layout of diagonal:
    // 0,1,...,baseband-1, k=1,k=-1, k=2,k=-2, .....,k=K,k=-K.
    ukkdiag_p diagonal; //pointer arrays to cost_dir. 
    int * pool_cost;
    DIRECTION_MATRIX * pool_dir;
    DIRECTION_MATRIX * pool_gapnum;
    // matrix P and Q for affine cost. 
    // according to "An Improved Algorithm for Matching Biological Sequence".
    // w_k = gap_cost * num_of_gaps(=k) + gap_opening  
    // D_m_n = Min( D_m-1_n-1 + d(a_m,b_n), P_m_n , Q_m_n )
    // P_m_n = Min( D_m-k_n + w_k ), k=[1,m] = Min( D_m-1_n + w_1 , P_m-1_n + gap_cost )
    // Q_m_n = Min( D_m_n-k + w_k ), k=[1,m] = Min( D_m_n-1 + w_1 , Q_m_n-1 + gap_cost )
    int total_len_affine;
    int * pool_affP;
    int * pool_affQ;
    int * pool_affD;
};

typedef struct newkkmat * newkkmat_p;

#ifdef _WIN32
__inline int 
#else
inline int 
#endif
newkk_algn (const seqt s1, const seqt s2, MAT_SIZE s1_len, MAT_SIZE s2_len, int go, const cmt c, newkkmat_p m, int swaped);


