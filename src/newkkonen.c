/* POY 5.0. A phylogenetic analysis program using Dynamic Homologies.         *\
(* Copyright (C) 2013 Andrés Varón, Lin Hong, Nicholas Lucaroni, Ward Wheeler,*)
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
#include <malloc.h>
#include <limits.h>
#include <string.h>
#include <assert.h>
#include <caml/mlvalues.h>
#include <caml/memory.h>
#include <caml/custom.h>
#include <caml/fail.h>

#include "seq.h"
#include "queue_with_linkedlist.h"
#include "newkkonen.h"
//trivial case
#define MUCH_LONGER 100

//has gap
#define GAPCODE 16
#define HAS_GAP(code) ( code&GAPCODE )

//use this to kick gap out of code
#define NOGAP 15

//directions
#define START 0
#define ALIGN_TO_ALIGN 1
#define ALIGN_TO_INSERT 2
#define ALIGN_TO_DELETE 4
#define ALIGN_TO_DIAGONAL 8
#define DO_DELETE 16
#define DO_INSERT 32
//just like in algn.c
//END_INSERT, END_DELETE and END_DIAG are where gap opening start
//we need these for affine backtrace, just put a sign END_xxx to the position
//with gap opening, where the gaps 'start'. I know the name END_xxxx is a little
//confusing, but for backtrace function, those positions are the 'end' of aligning gaps.
#define END_INSERT 64
#define END_DELETE 128
#define END_DIAG 256
#define DO_DIAG 512
//when gap extention cost in horizontal and vertical are the same
#define INSERT_EQ_DELETE 1024

#define has_flag(dir,flag) (dir&flag)

inline int my_add(a, b) {
    if(a>=INT_MAX/2) 
        return INT_MAX/2; 
    else if (b>=INT_MAX/2) 
        return INT_MAX/2; 
    else return(a+b);
}

enum MODE { m_todo, m_delete, m_insert, m_diagonal, m_align } backtrace_mode;

void print_mode (enum MODE mode)
{
    if (mode == m_todo) printf("mode = m_todo\n");
    if (mode == m_delete) printf("mode = m_delete\n");
    if (mode == m_insert) printf("mode = m_insert\n");
    if (mode == m_align) printf("mode = m_align\n");
}
//some explanation for directions and affine
//
//ALIGN_TO_ALIGN
//like other ALIGN_TO_xxxx, this direction sign is from algn.c, we can think it as 'DO_ALIGN'.
//
//ALIGN_TO_INSERT and ALIGN_TO_DELETE
//these two are from algn.c, they are not usefull in newkkonen.c
//
//
//
//END_HORIZONTAL and END_VERTICAL
//consider the case below 
//let gap opening cost 10, other operation cost 1
//     0  1  2  3  4  5  6  7
//     -  t  a  g  c  a  t  g
//0 -    11 12 13 14 15 16 17 
//1 t 11  0 11 12 13 14 15 16   
//2 a 12 11  0 11 12 13 14 15
//3 a 13 12 11  1 12 12 14 15
//                 |____| 
//
//pos (3,6) get cost=14 from two directions, one is alignment, the other is
//insertion.For insertion, it's not directly from pos(3,5), but from two gap extensions away:pos (3,4).
//we acctually favor insert/delete over alignment, so when we reach
//pos(3,6), next step will be pos(3,5), but cost=12 in pos(3,5) is only from
//alignment, follow that will take us to pos(2,4), which is wrong. 
//we need a way to continue insertion. that's why
//we have flag END_HORIZONTAL and END_VERTICAL and a list of mode = {m_todo,m_vertical,m_xxxx,blablabal} 
//mode is the one decide which direction we are going next
//we put a 'END_HORIZONTAL' sign in algn.c/'END_INSERT' sign in newkkonen.c at the position where gap opening happens, 
//in the case above, pos(3,4). 
//later in backtrace, we set mode to m_todo when END_HORIZONTAL show up in current position.
//remember we only change mode based on direction_matrix when mode=m_todo, otherwise we continue
//what we are doing -- either insertion or deletion.
//
//
//INSERT_EQ_DELETE
//The logic above has a flaw, if gap extention cost in horizontal and vertical are the same for some
//position. we do need to stop what we are doing(insertion or deletion), pick direction again. 
//So we introduce another sign 'INSERT_EQ_DELETE'. when insert ext and delete ext cost
//are the same, put this down for traceback.
//

/*
 * GO : gap opening cost passed in
 *
 * v  : vertical, from right
 * h  : horizontal, from left
 * 
 *
 * go is the gap opening cost function
 * go(Ai) = [ if i=1 and Ai contains gap then 0
 *            if i>1 and Ai-1 contains NO gap, Ai contains gap then 0
 *            otherwise GO ]
 *
 * ge is the gap extension cost, which in our case is just the cost from align no-gap base to a gap
 * ge(x) = [ if x has gap then 0 else b ]
 *
 * go' is the extra cost funtion when not selecting gap in Ai means splitting an indel block
   go'(Ai,Bj) = subst(Ai no gap,Bj) + [ if Ai contains no Gap then 0, otherwise GO ]


   g[i,j] = min {
      g[i-1,j-1] + subst(Ai,Bj)
      d[i-1,j-1] + subst(Ai,Bj) + go(A,i) + go(B,j)
      v[i-1,j-1] + go'(Bj,Ai)
      h[i-1,j-1] + go'(Ai,Bj)
     }

   d[i,j] = diag(Ai,Bj) + min [d[i-1,j-1], g[i-1,j-1] + go(Ai) + go(Bi)]

   diag(Ai,Bj) = [ if Ai has gap and Bj has gap then 0
                    else Inf ]
   
   h[i.j] = min [ h[i,j-1] + ge(Bj), d[i,j-1] + ge(Bj) + go(Bj)  ]

   v[i,j] = min [ v[i-1,j-1] + ge(Ai) , d[i-1,j] + ge(Ai) + go(Ai) ]

*/


#define my_prepend(a,b) assert (a->cap > a->len); \
    (a)->begin = (((a)->begin) - 1); \
    ((a)->len = 1 + (a)->len); *((a)->begin) = b
#define my_get(a,b) ((a)->begin)[b]

#define MAX(a,b) ( ((a)>(b))? (a):(b) )
#define MIN(a,b) ( ((a)<(b))? (a):(b) )
#define MAX3(a,b,c) ( (MAX(a,b)>MAX(b,c))? (MAX(a,b)):(MAX(b,c)) )
#define MAX4(a,b,c,d) ( MAX3(a,b,c)>d ) ? ( MAX3(a,b,c):d )

//get gap openning cost when start a indel block
#ifdef _WIN32
__inline int
#else
inline int
#endif
get_go(int base,int prebase, int idx, int gap_opening)
{
    if (idx ==1 && HAS_GAP(base) ) return 0;
    else if ( idx > 1 && !HAS_GAP(prebase) && HAS_GAP(base) ) return 0;
    else return gap_opening;
}

//get gap openning cost when extending a indel block
#ifdef _WIN32
__inline int
#else
inline int
#endif
get_extgo(int base,int prebase, int idx, int gap_opening)
{
    if ( HAS_GAP(prebase) && !HAS_GAP(base) ) return gap_opening;
    else return 0;
}


#ifdef _WIN32
__inline int
#else
inline int
#endif
get_diag (int base1, int base2) {
    if (HAS_GAP(base1) && HAS_GAP(base2)) return 0;
    else return INT_MAX/2;
}




void print_direction (DIRECTION_MATRIX dir)
{
    printf("[");
    //if (has_flag(dir,DO_ALIGN)) printf("align,");
    if (has_flag(dir,ALIGN_TO_ALIGN)) printf("align,");
    if (has_flag(dir,DO_DELETE)) printf("delete,");
    if (has_flag(dir,DO_INSERT)) printf("insert,");
    if (has_flag(dir,END_INSERT)) printf("end insert,");
    if (has_flag(dir,END_DELETE)) printf("end delete,");
    if (has_flag(dir,INSERT_EQ_DELETE)) printf("insert ext = delete ext,");
    if (has_flag(dir,ALIGN_TO_INSERT)) printf("align to insert,");
    if (has_flag(dir,ALIGN_TO_DELETE)) printf("align to delete,");
    if (has_flag(dir,ALIGN_TO_DIAGONAL)) printf("align to diagnoal,");
    printf("] ");
}


#ifdef _WIN32
__inline int
#else
inline int
#endif
get_k (newkkmat_p mp)
{
    if(mp == NULL) failwith ("newkkonen.get_k,NULL newkkmat");
    return mp->k;
};

value
newkkonen_CAML_get_k (value m)
{
    CAMLparam1(m);
    newkkmat_p tmp;
    tmp = Newkkmat_struct(m);
    CAMLreturn(Val_int(get_k(tmp)));
}

void print_diagonal (ukkdiag_p diag)
{
    printf("diagonal,len=%d",diag->len); 
    
};

void print_int_arr (int * arr, int size)
{
    int i; 
    printf("[");
    for (i=0;i<size;i++)
        printf("%d,",arr[i]);
    printf("]\n");
    fflush(stdout);
};

//we need to reset cells to 0 if we are doing speeding up version of newkkonen.
//in the speeding up version, we compare the old cost to new cost in newkkmat,
//old cost are supposed to be 0 if we haven't fill it yet, but we reuse memory all the time, that might not be true.
void
reset_mat (newkkmat_p m)
{
    int debug = 0;
    MAT_SIZE total_len = m->total_len;
    if (debug) printf ("reset matrix to 0, total_len=%li\n",total_len);
    //reset cost,we compare old cost to new cost in update_a_cell
    memset(m->pool_cost,(-1),total_len*sizeof(int));
    //these two are not necessary.
    //DIRECTION_MATRIX v = (DIRECTION_MATRIX) (-1);
    //memset(m->pool_dir,v,total_len*sizeof(DIRECTION_MATRIX));
    //memset(m->pool_gapnum,v,total_len*sizeof(DIRECTION_MATRIX));
    return;
}

void
expand_mat (newkkmat_p m,MAT_SIZE newk,MAT_SIZE oldk,int affine)
{
    int debug = 0;
    MAT_SIZE expand_diag_size = 1;
    MAT_SIZE expand_diagonal = 1;
    ukkdiag_p thisdiag;
    MAT_SIZE i=0;
    int diaglen;
    MAT_SIZE oldsize = m->diag_size_in_use;
    //get diagonal array size, and set diag_size_in_use
    MAT_SIZE newsize = oldsize+(newk-oldk)*2;
    m->k = newk;
    m->diag_size_in_use = newsize; 
    //if (debug) { printf ("newkkonen.c, expand_mat,newk=%li(oldk=%li),newsize of diagarr=%li(old=%li),baseband=%li,affine=%d\n", newk,oldk,newsize,oldsize,baseband,affine); fflush(stdout); }
    //set expand sign
    //if (debug) { printf ("newsize=%li <=> m->diag_size=%li\n", newsize, m->diag_size); fflush(stdout);}
    if (newsize <= m->diag_size) expand_diag_size = 0;
    else 
    {
        if (debug) {printf ("we NEED to expand diagonal_size_arr\n"); fflush(stdout); }
        expand_diag_size = 1;
        m->diagonal_size_arr = (int *)realloc(m->diagonal_size_arr,newsize*sizeof(int));
        m->diagonal = (ukkdiag_p ) realloc (m->diagonal,newsize*sizeof(struct cost_dir));
    }
    if ((m->diagonal==NULL)|| (m->diagonal_size_arr==NULL)) 
	failwith ("newkkonen.expand_mat,NULL diagnoal array,or NULL diagonal_size_arr");
    for (i=oldsize;i<newsize;i=i+2)
    {
        (m->diagonal_size_arr)[i] =  (m->diagonal_size_arr)[i-1] - 1;
        (m->diagonal_size_arr)[i+1] =  (m->diagonal_size_arr)[i-1] - 1;
    }
    
    for (i=oldsize;i<newsize;i++) //update total_len_in_use
    {
        diaglen = (m->diagonal_size_arr)[i];
        m->total_len_in_use += diaglen;
        thisdiag = m->diagonal;
        thisdiag += i;
        thisdiag->len = diaglen;
    }
    //set expand sign
/*    if (debug) {
        printf ("total_len_in_use=%li <=> total_len=%li\n", m->total_len_in_use, m->total_len);
        fflush(stdout); }*/
    if (m->total_len_in_use <= m->total_len) expand_diagonal = 0;
    else 
    {
        //if (debug) { printf ("we NEED to expand diagonal,total_len=%li,sizeof int=%li,sizeof MAT_SIZE=%lu,sizeof DIRECTION_MATRIX=%lu\n",m->total_len_in_use,sizeof(int),sizeof(MAT_SIZE),sizeof(DIRECTION_MATRIX)); fflush(stdout);}
        expand_diagonal = 1;
        MAT_SIZE len_we_need = m->total_len_in_use;
        if (len_we_need < 0) { printf("Warning: size to allocation < 0 , might be a integer overflow, are you runing on 32bit machine?\n"); fflush(stdout); }
        m->pool_cost = realloc (m->pool_cost,len_we_need*sizeof(int));
        m->pool_dir = realloc (m->pool_dir,len_we_need*sizeof(DIRECTION_MATRIX));
        //we have two gapnum array now
        m->pool_gapnum = realloc (m->pool_gapnum,2*len_we_need*sizeof(DIRECTION_MATRIX));
        if (affine) {
            m->pool_affP = realloc (m->pool_affP,len_we_need*sizeof(int));
            m->pool_affQ = realloc (m->pool_affQ,len_we_need*sizeof(int));
            m->pool_affED = realloc (m->pool_affED,len_we_need*sizeof(int));
            m->pool_affCD = realloc (m->pool_affCD,len_we_need*sizeof(int));
        }
    }
        ukkdiag_p emptydiag;
        int * thiscost= m->pool_cost;
        //define these two pointer anyway 
        int * thisaffP = m->pool_affP;
        int * thisaffQ = m->pool_affQ;
        int * thisaffED = m->pool_affED;
        int * thisaffCD = m->pool_affCD;
        DIRECTION_MATRIX * thisdir = m->pool_dir;
        DIRECTION_MATRIX * thisgap = m->pool_gapnum;
        //we need to set pointers of each diagonal again. realloc might move the
        //whole block to some new place, those diagonals pointing to the old
        //block of memory will be useless.
        emptydiag = m->diagonal;
        for (i=0;i<oldsize;i++)
        {
            diaglen = emptydiag->len;
            assert(diaglen>0);
            //printf ("set diag#%d(%d);",i,diaglen); fflush(stdout);
            emptydiag -> costarr = thiscost;
            emptydiag -> dirarr = thisdir;
            emptydiag -> gapnumarr1 = thisgap;
            thisgap += diaglen;
            emptydiag -> gapnumarr2 = thisgap;
            if (affine) {
            emptydiag -> affParr = thisaffP;
            emptydiag -> affQarr = thisaffQ;
            emptydiag -> affEDarr = thisaffED;
            emptydiag -> affCDarr = thisaffCD;
            thisaffP += diaglen;
            thisaffQ += diaglen;
            thisaffED += diaglen;
            thisaffCD += diaglen;
            }
            thiscost += diaglen;
            thisdir += diaglen;
            thisgap += diaglen;
            emptydiag++;
        }
        for (i=oldsize;i<newsize;i++)
        {
            diaglen = emptydiag->len;
            if(diaglen<=0) failwith ("ERROR: newkkonen.expand_mat , diaglen<0");
            //printf ("set diag#%d(%d);",i,diaglen); fflush(stdout);
            emptydiag -> costarr = thiscost;
            //init new cost array to (-1),we compare old cost to new cost in
            //update_a_cell, if old cost is not initialized, there might be bus
            //error.
            memset(thiscost,(-1),diaglen*sizeof(int));
            emptydiag -> dirarr = thisdir;
            emptydiag -> gapnumarr1 = thisgap;
            thisgap += diaglen;
            emptydiag -> gapnumarr2 = thisgap;
            if (affine) {
            emptydiag -> affParr = thisaffP;
            emptydiag -> affQarr = thisaffQ;
            emptydiag -> affEDarr = thisaffED;
            emptydiag -> affCDarr = thisaffCD;
            }
            if (i<newsize-1) //don't move out of range
            { 
            emptydiag ++;
            thiscost += diaglen ;
            thisdir += diaglen;
            thisgap += diaglen; 
            if (affine) {
                thisaffP += diaglen;
                thisaffQ += diaglen; 
                thisaffED += diaglen; 
                thisaffCD += diaglen; 
            }
            }
        }
    //}
    //update diag_size and total_len if we expand our memory
    if (expand_diag_size) { m->diag_size = m->diag_size_in_use; }
    if (expand_diagonal) 
    { m->total_len = m->total_len_in_use; 
        if (affine) m->total_len_affine=m->total_len;
    }
    /*if (debug) 
    { printf ("Total size of memory we are using :%lu(%lu) cells, diag_size=%d(%d)\n",
            m->total_len_in_use,m->total_len,m->diag_size_in_use,m->diag_size); 
        fflush(stdout);} */
    /*for (i=0;i<newsize;i++)
    {
        fprintf (stdout,"i=%d,diaglen=%d;",i,(m->diagonal_size_arr)[i]);
    }
    printf("\n"); fflush(stdout);*/
};


int get_delta (const cmt c)
{
    int delta = cm_get_min_non0_cost(c);
    return delta; 
};



//[get_idx] fill in the idx in diag arrays -- (whichdiag,idx_in_my_diag) based on the input idx -- (i,j). 
// also return sign showing if this cell is on the left/right border of the new band
//if oldk>=0 , we also return first_time_update=1 if this cell is in the new band but not in the old band, thus is never been updated before.
void 
get_idx (int oldk, int i,int j,newkkmat_p m, int * whichdiag,int * idx_in_my_diag, int * at_leftborder, int * at_rightborder, int * first_time_update)
{
    int currentK = m->k;
    if(currentK<0) failwith("ERROR:newkkonen.get_idx,k<0");
    //though baseband is defined as MAT_SIZE in newkkonen.h, which might be long int. but actually we don't need that much, int is fine. if baseband itself is bigger than max_int, the memory requirement will blow off our computer anyway. 
    int bb = m->baseband;
    if(bb<=0) failwith("ERROR:newkkonen.get_idx,bb<=0"); 
    int k_ij= i-j, k_ji = j-i;
    //see if this cell is in the new left band, or new right band with newk
    int in_new_band_l=0, in_new_band_r=0;
    if(i>j) { 
        *whichdiag = 2*(i-j)+bb-1 ;
    } else if ( (i<=j)&&( (j-i)<bb ) ) { 
        *whichdiag = j-i; 
    } else {
        *whichdiag = bb + 2*(j-i+1 - bb) -1 -1;
    }
    if(*whichdiag <0) failwith ("ERROR,newkkonen.get_idx,whichdiag<0");
    if(*whichdiag > m->diag_size) failwith("ERROR,newkkonen.get_idx,whichdiag>m->diag_size");
    //check if this cell is the new band area--thus it will get updated for the
    //first time
    int tmpr = k_ji-bb+1;
    //if we care about in_new_band_r/in_new_band_l
    if ( oldk >= 0 ) {
        if (tmpr > oldk ) { in_new_band_r = 1; } else { in_new_band_r = 0; }
        if (k_ij > oldk) { in_new_band_l = 1; } else { in_new_band_l = 0; }
        if (in_new_band_r||in_new_band_l) { *first_time_update = 1; } else { *first_time_update = 0; }
    }
    //check if this cell is on the new left/right border
    if ( tmpr==currentK ) { *at_rightborder=1;}
        else *at_rightborder=0;
    if (k_ij==currentK) {  *at_leftborder = 1; } 
        else *at_leftborder = 0;
    int idx; 
    if(i<=j) idx=i;  
    else idx=j;  
    *idx_in_my_diag = idx;
}


#ifdef _WIN32
__inline void 
#else
inline void 
#endif
get_ukkcost (int whichdiag, int idx_in_my_diag, newkkmat_p m, int * cost, DIRECTION_MATRIX * dir, DIRECTION_MATRIX *  s1_max_gapnum,DIRECTION_MATRIX *  s2_max_gapnum, int * affP, int * affQ, int * affED, int * affCD,int affine )
{
    int debug = 0;
    ukkdiag_p thisdiag = m->diagonal;
    if(whichdiag >= m->diag_size) failwith("ERROR: get_ukkcost,whichdiag >= diag size");
    thisdiag += whichdiag;
    int * thiscostarr = thisdiag->costarr;
    DIRECTION_MATRIX * thisdirarr = thisdiag->dirarr;
    DIRECTION_MATRIX * thisgapnumarr1 = thisdiag->gapnumarr1;
    DIRECTION_MATRIX * thisgapnumarr2 = thisdiag->gapnumarr2;
    if(idx_in_my_diag >= thisdiag->len) failwith("ERROR : get_ukkcost,idx_in_my_diag >= length of this diag");
    *cost = thiscostarr[idx_in_my_diag];
    *dir = thisdirarr[idx_in_my_diag];
    *s1_max_gapnum = thisgapnumarr1[idx_in_my_diag];
    *s2_max_gapnum = thisgapnumarr2[idx_in_my_diag];
    if (affine) {
    int * thisParr = thisdiag->affParr;
    int * thisQarr = thisdiag->affQarr;
    int * thisEDarr = thisdiag->affEDarr;
    int * thisCDarr = thisdiag->affCDarr;
    *affP = thisParr[idx_in_my_diag];
    *affQ = thisQarr[idx_in_my_diag];
    *affCD = thisCDarr[idx_in_my_diag];
    *affED = thisEDarr[idx_in_my_diag];
    }
    if (debug) { printf("get ukkcost,diag#.%d,idx#.%d,cost=%d,dir=%d,mapgapnum=%d/%d,affP=%d,affQ=%d,affED=%d,affCD=%d,affine=%d\n",
            whichdiag,idx_in_my_diag,*cost,*dir,*s1_max_gapnum,*s2_max_gapnum,*affP,*affQ,*affED,*affCD,affine); fflush(stdout);}
};


void set_ukkcost (int whichdiag, int idx_in_my_diag, newkkmat_p m, int cost, DIRECTION_MATRIX dir,DIRECTION_MATRIX  s1_max_gapnum,DIRECTION_MATRIX  s2_max_gapnum, int affP, int affQ, int affED, int affCD, int affine)
{
    int debug = 0;
    ukkdiag_p thisdiag = m->diagonal;
    if ((idx_in_my_diag >= thisdiag->len) || (idx_in_my_diag<0)) { debug=1; }
    if (debug) { 
        printf ("set ukkcost, diag#.%d,idx#.%d,cost=%d,dir=%d,gapnum=%d/%d,affP=%d,affQ=%d,affED=%d,affCD=%d,",
                whichdiag,idx_in_my_diag,cost,dir,s1_max_gapnum,s2_max_gapnum,affP,affQ,affED,affCD); 
        print_direction(dir);
        printf ("\n");
        fflush(stdout);}
    thisdiag += whichdiag;
    int * thiscostarr = thisdiag->costarr;
    DIRECTION_MATRIX * thisdirarr = thisdiag->dirarr;
    DIRECTION_MATRIX * thisgapnumarr1 = thisdiag->gapnumarr1;
    DIRECTION_MATRIX * thisgapnumarr2 = thisdiag->gapnumarr2;
    if(idx_in_my_diag >= thisdiag->len) failwith("ERROR:set_ukkcost,idx in my diag >= this diag's len");
    thiscostarr[idx_in_my_diag] = cost;
    thisdirarr[idx_in_my_diag] = dir;
    thisgapnumarr1[idx_in_my_diag] = s1_max_gapnum;
    thisgapnumarr2[idx_in_my_diag] = s2_max_gapnum;
    if (affine) {
        int * thisParr = thisdiag->affParr;
        int * thisQarr = thisdiag->affQarr;
        int * thisEDarr = thisdiag->affEDarr;
        int * thisCDarr = thisdiag->affCDarr;
    thisParr[idx_in_my_diag] = affP;
    thisQarr[idx_in_my_diag] = affQ;
    thisEDarr[idx_in_my_diag] = affED;
    thisCDarr[idx_in_my_diag] = affCD;
    }

};

void sanity_check2 (int x)
{
    if ((x==1)||(x==2)||(x==4)||(x==8)||(x==16)) {}
    else {
        printf ("WARNING: we have a x=%d, not a dna code\n", x);
        fflush(stdout);
    }
};

//get_cmcost return cost between code a and code b from cost matrix c.
//when exclude_gap_for_a is true, we get rid of the gap from a before getting the cost, exclude_gap_for_b is similar
void get_cmcost (const cmt c, int a, int b, int * res, int exclude_gap_for_a, int exclude_gap_for_b)
{
    if (exclude_gap_for_a) a = NOGAP & a;
    if (exclude_gap_for_b) b = NOGAP & b;
    if(cm_check_level(c) == 1)
        *res = cm_get_cost (c->cost, a, b, c->map_sz+1);
    else
        *res = cm_calc_cost (c->cost, a, b, c->lcm);
        
};

int in_non_change_zone (int i, int j,int oldk, int lenX, int lenY)
{
    int lefttopi=0, lefttopj=0;
    int leftbottomi=oldk;
    int righttopj=lenY-lenX+oldk;
    return ( i>=lefttopi) && (i<= leftbottomi) && (j>= lefttopj) && (j<=righttopj);
};

int outside_diagonal_area (int i, int j, int newk, int lenX, int lenY)
{
    return ( (i>j)&&( (i-j)>newk) ) || ( (i<j)&&((j-i)>(lenY-lenX+newk)) );
};


void assign_best_cost_and_direction 
(int costM, int costDiag, int costL, int costR, 
 int thisP, int thisQ, int thisED, int thisCD,
 int open_costDiag, int ext_costDiag, int open_costL, int ext_costL, int open_costR, int ext_costR, 
 DIRECTION_MATRIX gapnum1_fromM,DIRECTION_MATRIX gapnum1_fromL, DIRECTION_MATRIX gapnum1_fromR,
 DIRECTION_MATRIX gapnum2_fromM,DIRECTION_MATRIX gapnum2_fromL,DIRECTION_MATRIX gapnum2_fromR, 
 int whichdiag,int idx_in_my_diag, newkkmat_p m, int affine,int debug)
{
    if(debug) printf("\nassign_best_cost_and_direction,costM=%d,costDiag=%d,costL=%d,costR=%d,ext_costR=%d,open_costR=%d,ext_costL=%d,open_costL=%d,ext_costDiag=%d,open_costDiag=%d\n",
     costM,costDiag,costL,costR,ext_costR,open_costR,ext_costL,open_costR,ext_costDiag,open_costDiag);
    //do gaps on horizontal
    DIRECTION_MATRIX bestdir = DO_INSERT;
    int bestcost = costL;
    int resgapnum1 = gapnum1_fromL+1; 
    int resgapnum2 = gapnum2_fromL;
    //do gaps on vertical
    if (costR<=bestcost) {
        if (costR<bestcost) {
            bestcost = costR;
            bestdir = DO_DELETE;
            resgapnum1 = gapnum1_fromR;
            resgapnum2 = gapnum2_fromR+1;
        }
        else {
            bestdir = bestdir | DO_DELETE;
            resgapnum1 = MAX(resgapnum1,gapnum1_fromR);
            resgapnum2 = MAX(resgapnum2,gapnum2_fromR+1);
        }
    }
    //do alignment on diagonal
    if (costM<=bestcost) {
        if (costM<bestcost) {
            bestcost = thisCD;
            bestdir = ALIGN_TO_ALIGN;
            resgapnum1 = gapnum1_fromM;
            resgapnum2 = gapnum2_fromM;
        }
        else {
            //this is not complete, for affine, thisCD could mean Delete2Algn,
            //  Insert2Algn or Algn2Algn, but for traceback, all three will lead
            //  us from i,j to i-1,j-1, so why bother to tell them apart?
            bestdir = bestdir | ALIGN_TO_ALIGN; 
            resgapnum1 = MAX(resgapnum1,gapnum1_fromM);
            resgapnum2 = MAX(resgapnum2,gapnum2_fromM);
        }
    }
    //do gaps on diagonal
    if (costDiag<=bestcost) {
        if (costDiag<bestcost) {
            bestcost = costDiag;
            bestdir = DO_DIAG;
            resgapnum1 = gapnum1_fromM;
            resgapnum2 = gapnum2_fromM;
        }
        else {
            bestdir = bestdir | DO_DIAG;
            resgapnum1 = MAX(resgapnum1,gapnum1_fromM);
            resgapnum2 = MAX(resgapnum2,gapnum2_fromM);
        }
    }
    if (ext_costR>=open_costR) bestdir = bestdir|END_DELETE;
    if (ext_costL>=open_costL) bestdir = bestdir|END_INSERT;
    if (ext_costDiag>=open_costDiag) bestdir = bestdir|END_DIAG;
    //for affine, we favor gap extention over all other directions.
    //that's why in the old code, once we start insert/delete, 
    //we won't stop until we reach the gap opening position. 
    //But if ext_costL = ext_costR we do need to stop, 
    //choose the direction again.
    if ( (ext_costR==ext_costL)&&(ext_costR==bestcost) ) 
        bestdir = bestdir|INSERT_EQ_DELETE;
    if (debug) {
        printf("best cost = %d, best dir = ",bestcost);
        print_direction(bestdir); 
        printf("\n");  fflush(stdout);
    }
    //update cost arrays
    set_ukkcost(whichdiag,idx_in_my_diag, m, bestcost, bestdir, resgapnum1,resgapnum2, thisP, thisQ, thisED, thisCD, affine);
}

//this is the function update internal/left border/right border cell.
void update_internal_cell (const seqt s1, const seqt s2,newkkmat_p m, const cmt c, int i, int j, int newk,int oldk, int go,int at_leftborder, int at_rightborder)
{
    int debug = 0; //for this function
    int debug2 = 0;//for assign_best_cost_and_direction
    //if (i==35&&j>=50) { debug = 1; debug2 = 1; }
    int gapcode = cm_get_gap (c);
    //get s1[i],s1[i-1],s2[j],s2[j-1]
    int base1 = seq_get(s1,i); int base2 = seq_get(s2,j);
    int prev_base1 = seq_get(s1,i-1); int prev_base2 = seq_get(s2,j-1);
    //addcost will be reused in all three directions
    int addcost = 0;
    //set affine to 1 if gap opening is set
    int affine=0;
    if (go>=0) affine=1;
    int realgo=0;
    if (go>0) realgo=go;
    //init affP and affQ anyway
    int thisP=INT_MAX/2, thisQ=INT_MAX/2, thisED=INT_MAX/2, thisCD=INT_MAX/2;
    int affP=INT_MAX/2, affQ=INT_MAX/2, affED=INT_MAX/2, affCD=INT_MAX/2;
    //init gap opening cost 
    int go_base1=0, go_base2=0;
    int extgo_base1=0, extgo_base2=0;
    if(affine) {
        go_base1 = get_go (base1,prev_base1,i,realgo);
        go_base2 = get_go (base2,prev_base2,j,realgo);
	extgo_base1 = get_extgo(base1,prev_base1,i,realgo);
	extgo_base2 = get_extgo(base2,prev_base2,i,realgo);
    }
    if (debug) { printf(" update cell at [%d,%d], affine:%d,at left border:%d,at right border:%d\n",
            i,j,affine,at_leftborder,at_rightborder); fflush(stdout); }
    //fill costL with cost from left (insert)
    int costfromL=0; 
    DIRECTION_MATRIX dirL; DIRECTION_MATRIX gapnum1_fromL=0;DIRECTION_MATRIX gapnum2_fromL=0;
    int costL=0, ext_costL=0, open_costL=0;
    if (at_leftborder) {
        //we are at left border, there is no cost from left, set costL to big
        //number, also set ext_costL and open_costL to the same big number, make
        //sure that later function [assign_best_cost_and_direction] add 'END_INSERT' to this cell
        costfromL = costL = thisQ = ext_costL = open_costL = INT_MAX/2; 
        gapnum1_fromL=0; gapnum2_fromL=0;
    }
    else {
        if (j==0) {//leftmost column of matrix
            if (debug) printf ("j==0,set costfromL = costL = thisQ = ext_costL = open_costL <- INT_MAX/2\n");
            costfromL = costL = thisQ = ext_costL = open_costL = INT_MAX/2; 
            gapnum1_fromL=0; gapnum2_fromL=0;
        }
        else {
            int whichdiagL=0, idx_in_my_diagL=0, at_leftborderL=0, at_rightborderL=0;
            //for affine, get h[i,j-1]:affQ and d[i,j-1]:affED
            //for regular, coustfromL is enough
            get_idx(-1,i,j-1,m,&whichdiagL,&idx_in_my_diagL,&at_leftborderL,&at_rightborderL,NULL);
            get_ukkcost(whichdiagL,idx_in_my_diagL, m, &costfromL, &dirL, &gapnum1_fromL,&gapnum2_fromL, &affP, &affQ, &affED, &affCD, affine);
            //get gap exten cost: ge(Bj) 
            get_cmcost (c,base2,gapcode,&addcost,0,0);
            if (affine) {
                // h[i.j] = min [ h[i,j-1] + ge(Bj), d[i,j-1] + ge(Bj) + go(Bj)  ]
                ext_costL = my_add(affQ,addcost+extgo_base2);
                open_costL = my_add(affCD, addcost+go_base2);
                costL = MIN(open_costL,ext_costL);
                thisQ = costL;
                if (debug) {  printf ("costL <- min (%d,%d) = %d\n",ext_costL,open_costL,costL); fflush(stdout); }
            }
            else
                costL = costfromL + addcost;
        }
    }//end of at left border or not
    //fill costR with cost from right (delete)
    int costfromR=0;
    int costR = 0, ext_costR=0, open_costR=0;
    DIRECTION_MATRIX dirR; DIRECTION_MATRIX gapnum1_fromR=0; DIRECTION_MATRIX gapnum2_fromR=0;
    if (at_rightborder) {
        //we are at right border, there is no cost from right, set costR to big
        //number, also set ext_costR and open_costR to the same big number, make
        //sure that later function [assign_best_cost_and_direction] add 'END_DELETE' to this cell
        costfromR = costR = thisP = open_costR = ext_costR= INT_MAX/2; 
        gapnum1_fromR=0; gapnum2_fromR=0;
    }
    else { 
        //first line of matrix
        if (i==0) {
            if (debug) printf ("i==0,costfromR = costR = thisP = open_costR = ext_costR <- INT_MAX/2\n");
            costfromR = costR = thisP = open_costR = ext_costR= INT_MAX/2; 
            gapnum1_fromR=0; gapnum2_fromR=0;
        }
        else {
            int whichdiagR=0, idx_in_my_diagR=0, at_leftborderR=0, at_rightborderR=0;
            //for affine, get v[i-1,j] : affP and d[i-1,j] : affED
            //for regular, costfromR is enough.
            get_idx(-1,i-1,j,m,&whichdiagR,&idx_in_my_diagR,&at_leftborderR,&at_rightborderR,NULL);
            get_ukkcost(whichdiagR,idx_in_my_diagR, m, &costfromR, &dirR, &gapnum1_fromR, &gapnum2_fromR, &affP, &affQ, &affED,&affCD, affine); 
            //get gap exten cost : ge(Ai)
            get_cmcost (c,base1,gapcode,&addcost,0,0);
            if (affine) {
                //v[i,j] = min [ v[i-1,j-1] + ge(Ai) , d[i-1,j] + ge(Ai) + go(Ai) ]
                ext_costR = my_add(affP, addcost+extgo_base1);
                open_costR = my_add(affCD,addcost+go_base1);
                costR = MIN(open_costR,ext_costR);
                thisP = costR;
                if (debug) printf ("costR <- min (%d,%d)=%d\n",ext_costR,open_costR,costR);
            }
            else
                costR = costfromR + addcost;
        }
    }//end of at right border or not
    //fill costM and costDiag from diagonal 
    //for affine, 
    //g[i,j] = min {
    //  case1:g[i-1,j-1] + subst(Ai,Bj)
    //  case2:d[i-1,j-1] + subst(Ai,Bj) + go(A,i) + go(B,j)
    //  case3:v[i-1,j-1] + go'(Bj,Ai)
    //  case4:h[i-1,j-1] + go'(Ai,Bj)
    // }
    //d[i,j] = diag(Ai,Bj) + min [d[i-1,j-1], g[i-1,j-1] + go(Ai) + go(Bi)]
    int costfromM=0; 
    DIRECTION_MATRIX dirM; DIRECTION_MATRIX gapnum1_fromM=0; DIRECTION_MATRIX gapnum2_fromM=0;
    int costDiag=0, costM=0; 
    int ext_costDiag=0, open_costDiag=0;
    int costDiagfromR=0, costDiagfromL=0;
    int costDiagfromCD=0, costDiagfromED=0;
    if (i*j==0) { //first line or leftmost column of matrix
        if (debug) printf ("i==0||j==0,set costfromM = costDiag <- INT_MAX/2\n");
        costfromM = costM = costDiag = ext_costDiag = open_costDiag = costDiagfromR = costDiagfromL = costDiagfromCD = costDiagfromED = INT_MAX/2; 
        if (i+j==0) failwith ("newkkonen.c,update_internal_cell,we don't update (0,0) here\n");
    }
    else {
        //addcost : align (base1, base2)
        get_cmcost (c,base1,base2,&addcost,1,1);
        int whichdiagM=0, idx_in_my_diagM=0, at_leftborderM=0, at_rightborderM=0;
        get_idx(-1,i-1,j-1,m,&whichdiagM,&idx_in_my_diagM,&at_leftborderM,&at_rightborderM,NULL);
        //for affine, get previous d,g,v,h as affED,affCD,affQ,affP
        //for regular, costfromM is enough
        get_ukkcost(whichdiagM,idx_in_my_diagM, m, &costfromM, &dirM, &gapnum1_fromM, &gapnum2_fromM, &affP, &affQ, &affED, &affCD, affine);
        if (affine) {
            //g[i,j] case1: cost at [i-1,j-1] could from align Ai-1 and Bj-1
            //align Ai,Bj base on the previous align of Ai-1,Bj-1 : g[i-1,j-1] + subst(Ai,Bj)
            costDiagfromCD=my_add(affCD,addcost);
            thisCD = costDiagfromCD;
            //g[i,j] case 3&4: cost at [i-1,j-1] could from deletion or insertion
            //align Ai,Bj base on the previous deletion: costDiagfromR: v[i-1,j-1] + go'(Bj,Ai)
            //align Ai,Bj base on the previous insertion: costDiagfromL: h[i-1,j-1] + go'(Ai,Bj)
            //addcost_base1nogap : align (base1 without its gap , base2)
            //addcost_base2nogap : align (base2 without its gap , base1)
            int base1nogap = base1 & NOGAP;
            int base2nogap = base2 & NOGAP;
            int addcost_base1nogap, addcost_base2nogap;
            get_cmcost (c,base1nogap,base2,&addcost_base1nogap,0,0);
            get_cmcost (c,base1,base2nogap,&addcost_base2nogap,0,0);
            if (base1nogap == base1)
                costDiagfromR = affP + addcost_base1nogap;
            else
                costDiagfromR = affP + addcost_base1nogap + realgo;
            if (base2nogap == base2)
                costDiagfromL = affQ + addcost_base2nogap;
            else
                costDiagfromL = affQ + addcost_base2nogap + realgo;
            thisCD = MIN(thisCD,costDiagfromR);
            thisCD = MIN(thisCD,costDiagfromL);
            //g[i,j] case2: cost at [i-1,j-1] could from gap extension or gap opening
            //align Ai,Bj base on the previous diag gap extenion of opening: 
            //d[i-1,j-1] + subst(Ai,Bj) + go(A,i) + go(B,j)
            costDiagfromED = my_add(affED, go_base1+go_base2+addcost);
            thisCD = MIN(thisCD,costDiagfromED);
            //d[i,j] : align base1's gap with base2's gap
            int diag = get_diag (base1, base2);
            open_costDiag = go_base1 + go_base2;
            open_costDiag = my_add(affCD,open_costDiag);
            ext_costDiag = my_add(affED,diag);
            thisED = MIN(ext_costDiag,open_costDiag);
            //remember to udpate costM and costDiag
            costM = thisCD;
            costDiag = thisED;
            if (debug&&affine) printf(" thisCD <- min[%d,%d,%d,%d], thisED <-min[%d,%d] \n",
                    costDiagfromCD,costDiagfromED,costDiagfromR,costDiagfromL,open_costDiag,ext_costDiag );
        }
        else
            costM = addcost + costfromM;
    }
    //if (debug&&affine) printf("open_costL=%d,ext_costL=%d,open_costR=%d,ext_costR=%d,costM=%d,costDiag=%d\n", 
    //open_costL,ext_costL,open_costR,ext_costR,costM,costDiag);
    //update cost&dir for position (i,j)
    //no need to get whichdiag&idx_in_my_diag again, we can pass this in
    int whichdiag=0, idx_in_my_diag=0, this_at_leftborder=0, this_at_rightborder=0;
    get_idx(-1,i,j,m,&whichdiag,&idx_in_my_diag,&this_at_leftborder,&this_at_rightborder,NULL);
    assign_best_cost_and_direction(costM,costDiag,costL,costR,
                                   thisP,thisQ,thisED,thisCD,
                                   open_costDiag,ext_costDiag,open_costL,ext_costL,open_costR,ext_costR,
                                   gapnum1_fromM,gapnum1_fromL,gapnum1_fromR,
                                   gapnum2_fromM,gapnum2_fromL,gapnum2_fromR,
                                   whichdiag,idx_in_my_diag,m,affine,debug2);
};




int update_a_cell (const seqt s1, const seqt s2,newkkmat_p m, const cmt c, int i, int j, int newk, int oldk, int go, int getcost, int swaped)
{
    int debug = 0;
    int oldcost=0,newcost=0;
    int oldaffP=0, oldaffQ=0, oldaffED=0, oldaffCD=0;
    int newaffP=0, newaffQ=0, newaffED=0, newaffCD=0;
    DIRECTION_MATRIX oldgapnum1=0, oldgapnum2=0;
    int whichdiag=0, idx_in_my_diag=0, at_leftborder=0, at_rightborder=0, first_time_update=0;
    get_idx(oldk,i,j,m,&whichdiag,&idx_in_my_diag,&at_leftborder,&at_rightborder,&first_time_update);
    //if it's the first time we update this cell, there is no old cost to get
    if(first_time_update==1) getcost = 0;
    //define affine
    int affine=0;
    if (go>=0) affine=1; else affine=0;
    //sanity check
    if ((whichdiag<0)||(whichdiag>= m->diag_size_in_use)) {debug=1;}
    if (debug) { printf("update_a_cell,(%d,%d),affine=%d,diag#.%d,idx#.%d,leftB=%d,rightB=%d\n",i,j,affine,whichdiag,idx_in_my_diag,at_leftborder,at_rightborder); fflush(stdout);}
    if (getcost) {
        //get old cost first
        DIRECTION_MATRIX olddir;
        get_ukkcost(whichdiag,idx_in_my_diag, m, &oldcost, &olddir, &oldgapnum1,&oldgapnum2, &oldaffP, &oldaffQ, &oldaffED, &oldaffCD, affine); 
    }
    if ( i==0 && j==0 ) return 0;//we don't update (0,0)
    else 
    {
        /*if (at_leftborder&&at_rightborder) 
            update_central_diagonal_cell (s1, s2, m, c,i, j, newk, go);
        else*/ 
            update_internal_cell (s1, s2, m, c,i, j, newk, oldk, go, at_leftborder, at_rightborder);
        /*else if (at_leftborder)
            update_left_border_cell (s1, s2, m, c,i, j, newk, go);
        else if (at_rightborder)
            update_right_border_cell (s1, s2, m, c,i, j, newk, go);
        else
            update_internal_cell (s1, s2, m, c,i, j, newk, oldk, go, swaped); */
        if (getcost) {
            DIRECTION_MATRIX newdir, newgapnum1,newgapnum2;
            get_ukkcost(whichdiag,idx_in_my_diag, m, &newcost, &newdir, &newgapnum1,&newgapnum2, &newaffP, &newaffQ, &newaffED, &newaffCD, affine);
            if (debug) {
                printf("print direction:"); fflush(stdout);
                print_direction(newdir);
                printf("newcost=%d(gapnum=%d/%d,aff=P:%d/Q:%d/ED:%d/CD:%d),oldcost=%d(gapnum=%d/%d,aff=P:%d/Q:%d/ED:%d/CD:%d)\n",newcost,newgapnum1,newgapnum2,newaffP,newaffQ, newaffED, newaffCD, oldcost,oldgapnum1,oldgapnum2,oldaffP,oldaffQ,oldaffED,oldaffCD); 
                fflush(stdout);
                assert(newcost>=0);
            }
            if (affine) {
                if ( (newaffCD==oldaffCD)&&(newaffED==oldaffED)&&(newaffQ==oldaffQ)&&(newaffP==oldaffP)&&(newcost==oldcost)&&(newgapnum1==oldgapnum1)&&(newgapnum2==oldgapnum2) ) 
                { return 0; }
                else 
                {  return 1; }
            }
            else {
                if ((newcost==oldcost)&&(newgapnum1==oldgapnum1)&&(newgapnum2==oldgapnum2)) 
                {   return 0; }
                else {  return 1; }
            }
        }
        else
        {
            if (debug) {
                DIRECTION_MATRIX newdir, newgapnum1,newgapnum2;
                get_ukkcost(whichdiag,idx_in_my_diag, m, &newcost, &newdir, &newgapnum1,&newgapnum2, &newaffP, &newaffQ, &newaffED, &newaffCD, affine);
                printf("print direction:"); fflush(stdout);
                print_direction(newdir);
                printf("newcost=%d(gapnum=%d/%d,aff=%d/%d/%d/%d)\n",newcost,newgapnum1,newgapnum2,newaffP,newaffQ,newaffED, newaffCD); 
                fflush(stdout);}
                if(first_time_update==1) //first time updating this cell, set costchange to 1
                    return 1;
                else //we don't care about costchange, just return 0
                    return 0;
        }
    }
};


//function for speeding up start

/* this can be done with diag# and K
 * return true if (i,j) is in the new left band of newk, including old left border of oldk*/
int in_new_left_band (int i, int j,int newk,int oldk)
{
    if ( (i>=j)&&((i-j)<=newk)&&((i-j)>=oldk) ) return 1;
    else return 0;
};


//update a row in newkkonen matrix. if this is the first run, which means we are
//filling cells in baseband, we need to update each one of them. if not, we use
//two queue, prevq and thisq, to decide which cell we are going to update next.
void update_a_row (const seqt s1, const seqt s2,newkkmat_p m, const cmt c,int i,int startj,int endj,int last_updated_j,q_t prevq,q_t thisq,int newk,int oldk, int go, int swaped)
{
    int debug = 0;
    int next_j;//which cell we should update at this run
    int follow_prevq;//if we update this cell by following previous queue
    QUEUE_DATA msb=0;
    int get_cost_change = 0;
    if (newk == oldk)//first run, update all cells
    {
        if (debug) printf(" update_a_row,first run,update all,i=%d\n",i);
        //iteration is faster than recurtion
        int j=0;
        for (j=startj;j<=endj;j++)
        {
            update_a_cell (s1,s2,m,c,i,j,newk,oldk,go,get_cost_change,swaped);
        }
        /* recurtion way of doing it
        next_j = last_updated_j + 1;
        if (next_j<=endj) {
             update_a_cell (s1,s2,m,c,i,next_j,newk,go,0);
             update_a_row(s1,s2,m,c,i,startj,endj,next_j,prevq,thisq,newk,oldk,go);
        }
        else {}//done with this row
        */
    }
    else
    {
        if (debug) 
        {printf("\nupdate a row,i=%d,(startj,endj)=(%d,%d),last_updated_j=%d\n",i,startj,endj,last_updated_j); fflush(stdout);}
        //let's do this in iteration
        //init next_j with startj, follow_prevq set to false
        next_j = startj; follow_prevq=0;//last_updated_j=startj-1;
        do
        {
            //update cost of a cell
            if (debug) { printf("update this cell.(%d,%d)\n",i,next_j); fflush(stdout); }
            get_cost_change = 1;
            int costchange = update_a_cell (s1,s2,m,c,i,next_j,newk,oldk,go,get_cost_change, swaped);
            if (debug) { printf("costchange=%d\n",costchange); fflush(stdout); }
            //push next_j to thisq when necessary
            if ( costchange && !(in_new_left_band(i,next_j,newk,oldk)) ) {
                /*  we only push j to this queue when the cost changes
                    * in current cell. But if current cell is in the new left
                    * band of newk, including left border of old k, there is no
                    * need to push j.
                    * for example, oldk=0,newk=2,lenY-lenX=3
                    * when we are working on row 2. 
                    * N=newcost,O=oldcost,X= not in diagonal
                    *   0  1  2  3  4  5  6  7
                    * 0 ----------------------- row 0 is done.
                    * 1 ----------------------- row 1 is done.  
                    * 2 N  N  N  O  N  O  N  N  row 2 is being updated
                    * 3 X  ?  ?  ?  ?  ?  ?  ?  row 3 is next
                    *
                    * for this queue in row2, we push <0,1,2,4,6,7>, then pass it
                    * as previous queue to row3, but 0,1,2 are not necessary.
                    * because for row 3, we have to update pos(3,1),(3,2) anyway,
                    * for they are the new left band for newk=2.
                    * also pos(3,3) is the old left border with oldk=0, since
                    * it's left neighbor:pos(3,2) is always updated with a new cost,
                    * it will be updated always. 
                 */
                enqueue(thisq,next_j);
                if (debug) {printf("push %d to thisq\n",next_j); fflush(stdout); }
            }
            else {} //no need to add next_j to queue
            //update last_updated_j
            if ((costchange)) //set last_updated_j to next_j if we have better cost, 
            { last_updated_j = next_j; }
            else if (follow_prevq) {
                /*if this cell is being updated because the one
                * right above it--from preq, even the cost does
                * not change, we still need to theck the cell
                * right next to this one.
                *         j     j+1   
                *   i     N     O
                * i+1     O     ?--> might have better cost from (i).(j)
                */
                last_updated_j = next_j; }
            else 
                { last_updated_j = -1; }
            //update next_j
            if (last_updated_j == (-1) )
            {//no change in previous cell of current row, check preq of the row above us
                if (is_emptyqueue(prevq)) { msb=endj+1; follow_prevq=0; }
                else
                { dequeue(prevq,&msb); follow_prevq=1; }
                if (debug) { printf("no change in pre cell,follow msb=%d\n",msb); fflush(stdout); }
                next_j=msb;
            }
            else
            {//some change in previous cell
                //get head of previous queue, if there is any
                if (is_emptyqueue(prevq)) { msb = endj + 1; follow_prevq=0; }
                else
                { peekqueue(prevq,&msb); follow_prevq=1; }
                if (debug) {printf ("msb = %d,last_j=%d;",msb,last_updated_j); fflush(stdout); }
                if ( (last_updated_j+1)<msb ) {//continue with last_updated_j
                    next_j = last_updated_j + 1; follow_prevq=0; }
                else {//follow previous queue,
                    if(follow_prevq) { dequeue(prevq,&msb); }
                    else {}//just in case preq is empty
                    next_j = msb;
                }
            }
        } while ((next_j!=endj+1)&&(last_updated_j!=endj));
        //end of iteration 
    }//end of not working on baseband
    return;
};

//function for speeding up end

void ukktest (const seqt s1, const seqt s2,newkkmat_p m, const cmt c,int currentT, int p, int lenX, int lenY, int go, int * res_cost, DIRECTION_MATRIX * res_gapnum, int swaped)
{
    int debug = 0;
    //if go=(-1), we are not doing affine
    //set affine to 1 if gap opening is set
    int affine=0;
    if (go>=0) affine=1;
    int newk=p;
    if (p>=lenX) newk=lenX-1;
    MAT_SIZE bb = m->baseband; 
    MAT_SIZE old_size = m->diag_size_in_use;
    MAT_SIZE oldk;
    oldk = (old_size-bb)/2;
    if (debug) { printf("\n ukktest,newk=%d,bb=%li,oldsize=%li,oldk=%li,go=%d\n",newk,bb,old_size,oldk,go);fflush(stdout); }
    expand_mat (m,newk,oldk,affine);
    int i;
    int startj,endj,oldendj;
    //init two queue for speeding up function
    q_t prevq = (q_t) malloc(sizeof(struct Queue));
    init_empty_queue(prevq);
    q_t thisq = (q_t) malloc(sizeof(struct Queue));
    init_empty_queue(thisq);
    for (i=0;i<lenX;i++)
    {
        startj = ( (i-newk)>0 )?(i-newk):0 ;
        endj = ( (i+lenY-lenX+newk)<(lenY-1) )?(i+lenY-lenX+newk):(lenY-1) ;
        oldendj = ( (i+lenY-lenX+oldk)<(lenY-1) )?(i+lenY-lenX+oldk):(lenY-1) ;
        /*
	if (affine) {
            for (j= MAX(i-newk,0);j<=MIN(i+(lenY-lenX+newk),lenY-1);j++)
            {
            if (in_non_change_zone(i,j,oldk,lenX,lenY) ) {}
            else
                update_a_cell (s1,s2,m,c,i,j,newk,go,0);
            }
        }
        else { */
            // * call speed up function
            if ((debug)&&(is_emptyqueue(thisq)==0)) printf("Warning,thisQ is not empty before update_a_row\n");
            // special case: update the first row of matrix. 
            // for example if we only did j=0~5 during last
            // run(ukktest), now we work on j=0~7, 
            // nothing new can happen to 0~5, we need to start updating from (0,6),
            //
            // -1             N
            //    0 1 2 3 4 5 6 7
            //  0 O O O O O O ? ? -- working on row#0
            //  1   
            // remember that we empty prevq at the end of each run with
            // 'transfer_queue', last_updated_j is set to startj-1 to make sure
            // the first cell of each row will get updated -- that won't
            // help us here. hence here the solution, we imagine there is a row#-1,
            // we push cell 6 to prevq(6 = old end j + 1),[update_a_row] will pick
            // that up, start to update the cell right below (-1,6), which is (0,6).
            // since this is the first time we are updating (0,6), costchange will be true,
            // it get pushed to thisq, and last_updated_j is also set to 6.
            // we will continue to update (0,7) after (0,6).
            if (i==0) enqueue(prevq,oldendj+1);
            update_a_row (s1,s2,m,c,i,startj,endj,startj-1,prevq,thisq,newk,oldk,go,swaped);
            transfer_queue (thisq,prevq);
            //end of calling speed up function */
            /* old code start, no speeding up
            for (j= MAX(i-newk,0);j<=MIN(i+(lenY-lenX+newk),lenY-1);j++)
            {
                if (in_non_change_zone(i,j,oldk,lenX,lenY) ) {}
                else
                    update_a_cell (s1,s2,m,c,i,j,newk,go,0);
            }
            old code end */
        //}
    }
    int whichdiag=0, idx_in_my_diag=0, at_leftborder=0, at_rightborder=0;
    get_idx(-1,lenX-1,lenY-1,m,&whichdiag,&idx_in_my_diag,&at_leftborder,&at_rightborder,NULL);
    int cost; DIRECTION_MATRIX dir; DIRECTION_MATRIX gapnum1; DIRECTION_MATRIX gapnum2;
    int affP=INT_MAX/2, affQ=INT_MAX/2, affED=INT_MAX/2,affCD=INT_MAX/2; //init affP and affQ anyway
    get_ukkcost(whichdiag,idx_in_my_diag,m, &cost, &dir, &gapnum1,&gapnum2,&affP, &affQ, &affED, &affCD, affine);
    *res_cost = cost;
    if (gapnum1>gapnum2) *res_gapnum = gapnum1;
    else *res_gapnum = gapnum2;
};

int increaseT (const seqt s1, const seqt s2,newkkmat_p m, const cmt c,int newT, int lenX, int lenY, int go, int swaped)
{
    int debug = 0;
    int p = (newT - (lenY-lenX))/2;
    if (debug) { printf ("increaseT, lenX=%d, lenY=%d, newT=%d,p=%d,go=%d,swaped=%d\n",lenX, lenY, newT,p,go,swaped); fflush(stdout); }
    int res_cost=-1; 
    DIRECTION_MATRIX res_gapnum=-1;
    ukktest(s1,s2,m,c,newT,p,lenX,lenY,go,&res_cost,&res_gapnum, swaped);
    int newp = ( newT*2 - (lenY-lenX) )/2;
    assert(res_cost>=0); 
//    if ((res_cost<=newT)  //we accept the cost -- threshold by ukkonen
    //|| (res_cost<newT*2) //I know we will accept this cost in next round anyway,but this is not a shortcut,because the cost might be better in next round, we cannot stop here. 
//        || (newp>res_gapnum)) //max possible gap number -- threshold by ward's ukkonen -- newkkonen
    if ((p>res_gapnum)||(newp-lenY+1>=0))
    {
        if (debug) 
        {printf ("cost=%d,end by p(%d)>res_gapnum(%d) or newp(%d)>=lenY(%d)-1\n",res_cost,p,res_gapnum,newp,lenY);
        fflush(stdout); }
        return res_cost;
    }
    else //double the T, continue
    {
        if (debug) { printf ("cost=%d,p=%d,res_gapnum=%d,continue \n",res_cost,p,res_gapnum);
        fflush(stdout); }
        return increaseT (s1,s2,m, c,2*newT, lenX, lenY, go, swaped);
    };
};


void 
init_affine_mat (newkkmat_p m)
{
    int debug = 0;
    int oldlen = m->total_len;
    int oldafflen = m->total_len_affine;
    int i=0;
    if (oldlen>oldafflen)//need to alloc affine part of matrix
    {
        m->pool_affP = realloc (m->pool_affP,oldlen*sizeof(int));
        m->pool_affQ = realloc (m->pool_affQ,oldlen*sizeof(int));
        m->pool_affED = realloc (m->pool_affED,oldlen*sizeof(int));
        m->pool_affCD = realloc (m->pool_affCD,oldlen*sizeof(int));
        int * thisaffP = m->pool_affP;
        int * thisaffQ = m->pool_affQ; 
        int * thisaffED = m->pool_affED; 
        int * thisaffCD = m->pool_affCD; 
        ukkdiag_p thisdiag = m->diagonal;
        thisdiag -> affParr = thisaffP;
        thisdiag -> affQarr = thisaffQ;
        thisdiag -> affEDarr = thisaffED;
        thisdiag -> affCDarr = thisaffCD;
        int diaglen;
        int diag_arr_size = m->diag_size;
        if (debug) { printf("init space for affine, oldlen=%d,diag_arr_size=%d\n",oldlen,diag_arr_size); fflush(stdout);}
        for (i=0;i<diag_arr_size;i++) {
            //printf ("init set diag#%d(%d);",i,thisdiag->len); fflush(stdout);
            diaglen = thisdiag->len;
            assert(diaglen>0);
            thisdiag -> affParr = thisaffP;
            thisdiag -> affQarr = thisaffQ;
            thisdiag -> affEDarr = thisaffED;
            thisdiag -> affCDarr = thisaffCD;
            if (i<diag_arr_size-1) //don't move out of range
            {
            thisaffP += diaglen;
            thisaffQ += diaglen;
            thisaffED += diaglen;
            thisaffCD += diaglen;
            thisdiag ++;
            }
        }
        m->total_len_affine = oldlen;
    }
}

void
init_mat (MAT_SIZE lenX, MAT_SIZE lenY,newkkmat_p m,int affine)
{
    int debug = 0;
    int i=0;
    //expand sign
    int expand_diag_size=1;
    int expand_diagonal=1;
    MAT_SIZE len=0, oldlen=0;
    MAT_SIZE baseband=0;
    baseband = lenY-lenX+1;
    oldlen = m->total_len;
    len = lenX * baseband;
    //if (debug) { printf ("newkkonen init_mat,oldlen = %li, lenx=%li,leny=%li,baseband=%li,oldlenaff=%d\n",oldlen,lenX,lenY,baseband,m->total_len_affine); fflush(stdout); }
    assert(m != NULL);
    //k is init to 0
    m->k=0;
    //set baseband , this won't change during this alignment
    m->baseband = baseband;
    //init diag_size_in_use with baseband,may increase later
    m->diag_size_in_use = baseband; //array of diagonal size, start with [lenx,lenx,....]
    //set total_len_in_use to len, may increase later
    m->total_len_in_use = len;
/*    if (debug) { printf ("sizeofint=%lu,sizeoflong=%lu,len=%li <=> oldlen=%li\n", sizeof(int),sizeof(long),len, oldlen); fflush(stdout);}
*/
    //init affine part of matrix if necessray
    if (affine) init_affine_mat(m);
    //set expand sign
    if (len<=oldlen) 
    {
        expand_diagonal = 0;
        if (debug) { printf ("NO need to add new diagonal\n"); fflush(stdout);}
    }
    else 
    {
        expand_diagonal = 1;
        m->total_len = m->total_len_in_use;
        if (affine) m->total_len_affine = m->total_len;
        if (debug) { printf ("we NEED to add new diagonal\n"); fflush(stdout);}
    }
    if (debug) { 
        printf ("baseband=%li <=> m->diag_size=%li\n", baseband, m->diag_size); fflush(stdout); }
    if (baseband <= m->diag_size) 
    { 
        expand_diag_size=0; 
        if (debug) { printf("NO need to increase diag size\n"); fflush(stdout);}
    }
    else 
    {   
        expand_diag_size=1;
        m->diag_size = m->diag_size_in_use;
        if (debug) { printf("we NEED to increase diag size\n"); fflush(stdout);}
        m->diagonal_size_arr = (int *) realloc (m->diagonal_size_arr,baseband*sizeof(int));
        m->diagonal = (ukkdiag_p ) realloc (m->diagonal,baseband*sizeof(struct cost_dir));
    }
    //if (expand_diag_size==1)
    //{ 
    //    m->diagonal_size_arr = (int *) realloc (m->diagonal_size_arr,baseband*sizeof(int));
    //}
    assert(m->diagonal_size_arr != NULL);
    assert(m->diagonal != NULL);
    for (i=0;i<baseband;i++)
    {  
        int * tmp = m->diagonal_size_arr;
        tmp[i] = lenX;  }
    //set len of each diag 
    for (i=0;i<baseband;i++) 
    {
        ukkdiag_p thisdiag = m->diagonal;
        thisdiag += i;
        thisdiag -> len = lenX;
    }
    //set pointer of cost/dir/gap of each diag
    if (expand_diagonal==1) 
    {
        if (len<0) {printf("Warning: len to allocation < 0 , might be integer overflow. are you runing on 32bit machine? \n"); fflush(stdout);}
        m->pool_cost = realloc (m->pool_cost,len*sizeof(int));
        m->pool_dir = realloc (m->pool_dir,len*sizeof(DIRECTION_MATRIX));
        m->pool_gapnum = realloc (m->pool_gapnum,2*len*sizeof(DIRECTION_MATRIX));//2* because we need gapnum for both seq
        if (affine) {
            if (debug) {printf("expand affine P and Q diagnonal\n"); fflush(stdout);}
            m->pool_affP = realloc (m->pool_affP,len*sizeof(int));
            m->pool_affQ = realloc (m->pool_affQ,len*sizeof(int));
            m->pool_affED = realloc (m->pool_affED,len*sizeof(int));
            m->pool_affCD = realloc (m->pool_affCD,len*sizeof(int));
        }
        else
        {
            m->pool_affP = NULL;
            m->pool_affQ = NULL;
            m->pool_affED = NULL;
            m->pool_affCD = NULL;
        }
    }   
        int * thiscost = m->pool_cost;
        int * thisaffP = m->pool_affP;
        int * thisaffQ = m->pool_affQ; 
        int * thisaffED = m->pool_affED; 
        int * thisaffCD = m->pool_affCD; 
        DIRECTION_MATRIX * thisdir = m->pool_dir;
        DIRECTION_MATRIX * thisgap = m->pool_gapnum;
        ukkdiag_p thisdiag = m->diagonal;
        for (i=0;i<baseband;i++) {
            //printf ("init set diag#%d(%d);",i,thisdiag->len); fflush(stdout);
            thisdiag -> costarr = thiscost;
            thisdiag -> dirarr = thisdir;
            thisdiag -> gapnumarr1 = thisgap;
            thisgap += lenX;
            thisdiag -> gapnumarr2 = thisgap;
            if (affine) {
            thisdiag -> affParr = thisaffP;
            thisdiag -> affQarr = thisaffQ;
            thisdiag -> affEDarr = thisaffED;
            thisdiag -> affCDarr = thisaffCD;
            }
            set_ukkcost(i,0,m,0,0,0,0,0,0,0,0,affine);
            if (i<baseband-1) //don't move out of range
            {
            thisdiag ++;
            thiscost += lenX ;
            thisdir += lenX;
            thisgap += lenX; 
            if (affine) {
                thisaffP += lenX;
                thisaffQ += lenX;
                thisaffED += lenX;
                thisaffCD += lenX;
            }
            }
        }
    //}
};


int get_indel_cost (const seqt s, int len, const cmt c)
{
    int indelcost = 0;
    int thiscost, thiscode;
    int gapcode = cm_get_gap (c);
    int i;
    for (i=0;i<len;i++)
    {
        thiscode = seq_get (s,i);
        thiscost=0;
        get_cmcost(c,thiscode,gapcode,&thiscost,1,0);
        indelcost += thiscost;
    }
    return indelcost;

}

int trivial_algn (const seqt s1, const seqt s2, int s1_len, int s2_len, const cmt c)
{
    int debug = 0;
    if (debug)  { printf("trivial algn\n");fflush(stdout);}
    int indelcost1 = get_indel_cost (s1,s1_len,c);
    int indelcost2 = get_indel_cost (s2,s2_len,c);
    if (debug) { printf("end of trivial algn, cost=%d\n",indelcost1+indelcost2); fflush(stdout);}
    return (indelcost1 + indelcost2);
}



int
newkk_algn (const seqt s1, const seqt s2, MAT_SIZE s1_len, MAT_SIZE s2_len, int go, const cmt c, newkkmat_p m, int swaped) {
    if(s1_len>s2_len) failwith("ERROR: newkkonen.newkk_algn, s1 len > s2 len");
    int debug = 0;
    int debug2 = 0;
    //set affine to 1 if gap openning is set
    int affine=0;
    if (go>=0) affine=1; 
    int realgo=0;
    if (go>0) realgo=go;
    //if (debug)  { printf("newkk_algn,lenx=%li,leny=%li,affine=%d\n",s1_len,s2_len,affine);fflush(stdout);}
    if (s1_len * MUCH_LONGER < s2_len) 
    {
        return trivial_algn (s1,s2,s1_len,s2_len,c);
    }
    else 
    {
        //reset newkkmat to 0 if we are reusing the memory block.
        init_mat (s1_len,s2_len,m,affine);
        if (affine &&(m->total_len > 0)) reset_mat (m);
        int gapcode = cm_get_gap (c); 
        int delta = get_delta (c);
        int i; 
        int bb = m->baseband;
        assert(bb>0);
        if (debug)
        {
            printf("newkkonen.newkk_algn, baseband=%d,gapcode=%d,delta=%d,go=%d(realgo=%d),swaped=%d\n",bb,gapcode,delta,go,realgo,swaped);
            fflush(stdout);
        }
        //init first row , algn gaps to longer sequence
        //g[0,0]=0, d[0,0]=INF, v[0,0]=GO, h[0,0]=GO
        set_ukkcost(0,0,m,0,START,0,0,realgo,realgo,INT_MAX/2,0,affine);
        if (debug) {printf("init first row\n"); fflush(stdout);}
        for (i=1;i<bb;i++)
        {
           if (debug2) { printf("init (0,%d) ",i); fflush(stdout); }
            int whichdiag=0, idx_in_my_diag=0, at_leftborder=0, at_rightborder=0;
            get_idx(-1,0,i-1,m,&whichdiag,&idx_in_my_diag,&at_leftborder,&at_rightborder,NULL);
            int costfromL; DIRECTION_MATRIX predir; DIRECTION_MATRIX pre_gapnum1,pre_gapnum2; 
            //we init affP and affQ even we are not doing affine
            //g[0,i] = d[0,i] = v[0,i] = INF, h[0,i] = h[0,i-1] + GO
            int affP=INT_MAX/2, affQ=INT_MAX/2, affED=INT_MAX/2, affCD=INT_MAX/2; 
            get_ukkcost(whichdiag,idx_in_my_diag,m, &costfromL, &predir, &pre_gapnum1,&pre_gapnum2, &affP, &affQ, &affED, &affCD, affine);
            whichdiag=0; idx_in_my_diag=0; at_leftborder=0; at_rightborder=0;
            get_idx(-1,0,i,m,&whichdiag,&idx_in_my_diag,&at_leftborder,&at_rightborder,NULL);
            int costL=0, gocost = 0, addcost = 0, open_costL=0, ext_costL=0;
            int thisP=INT_MAX/2,thisQ=0,thisED=INT_MAX/2, thisCD=INT_MAX/2;
            int thisbase = seq_get (s2,i);
            int prevbase = seq_get (s2,i-1);
            int flag = gapcode & prevbase;
            int flag2 = gapcode & thisbase;
            get_cmcost(c,thisbase,gapcode,&addcost,0,0);
            if (debug2) { printf("prev has gap:%d, this has gap:%d,",flag,flag2); }
            if(affine) {
                if (i==1) {//we append the first gap to each seq, that doesn't count as prev base contains gap
                    if (flag2) gocost = 0; else gocost = realgo;
                }
                else {
                     if (!flag && flag2) gocost = 0; 
                            else gocost = realgo;

                }
                if ( flag && !flag2 ) ext_costL = affQ + addcost + realgo;
                else ext_costL = affQ + addcost;
                open_costL = costfromL + addcost + gocost;
                costL = MIN(open_costL, ext_costL); 
                thisQ = costL; }//end of if (affine)
            else { costL = costfromL + addcost; }//end of not affine
            if (debug2) { printf("with cost:%d [whichdiag=%d,idx=%d], gapnum:%d/%d\n ", costL,whichdiag,idx_in_my_diag,pre_gapnum1+1,pre_gapnum2+1); fflush(stdout);}
            set_ukkcost(whichdiag,idx_in_my_diag,m, costL, DO_INSERT | END_INSERT | END_DELETE | END_DIAG, pre_gapnum1+1, pre_gapnum2,thisP,thisQ,thisED,thisCD,affine);
        }
        int iniT = (m->baseband) * delta;
        if (debug) {printf ("done with first row, call increaseT with iniT=%d\n",iniT); fflush(stdout);}
        int rescost = increaseT (s1,s2,m,c,iniT,s1_len,s2_len,go,swaped); 
        if (debug) {printf("end of newkk_algn,cost=%d\n\n",rescost); fflush(stdout);}
        /* this is for cost compare test with full alignment
        FILE * outf;
        outf = fopen ("poy.out","w");
        fprintf(outf,"%d",rescost);
        fclose(outf);
        */
        return rescost;
    }
};



value 
newkkonen_CAML_algn (value s1, value s2, value c, value a,value swaped)
{
    CAMLparam5 (s1,s2,c,a,swaped);
    cmt tc;
    newkkmat_p mp;
    tc = Cost_matrix_struct(c);
    mp = Newkkmat_struct(a);
    int spd = Int_val(swaped);
    seqt s1p, s2p;
    Seq_custom_val(s1p,s1);
    Seq_custom_val(s2p,s2);
    MAT_SIZE s1_len, s2_len;
    s1_len = seq_get_len (s1p);
    s2_len = seq_get_len (s2p);
    assert (s2_len >= s1_len);
    int res;
    int debug = 0;
    if (debug) {
    printf("newkkonen_CAML_algn,no gap opening\n"); fflush(stdout);}
    res = newkk_algn (s1p, s2p, s1_len,s2_len,(-1),tc, mp,spd);
    CAMLreturn(Val_int(res));
};

value 
newkkonen_CAML_algn_affine (value s1, value s2, value c, value a, value swaped)
{
    CAMLparam5 (s1,s2,c,a, swaped);
    cmt tc;
    newkkmat_p mp;
    tc = Cost_matrix_struct(c);
    mp = Newkkmat_struct(a);
    int spd = Int_val(swaped);
    seqt s1p, s2p;
    Seq_custom_val(s1p,s1);
    Seq_custom_val(s2p,s2);
    MAT_SIZE s1_len, s2_len;
    s1_len = seq_get_len (s1p);
    s2_len = seq_get_len (s2p);
    assert (s2_len >= s1_len);
    int res;
    int gap_open = tc->gap_open;
    int debug = 0;
    if (debug) {
    printf("newkkonen_CAML_algn_affine,gap_open=%d\n",gap_open); fflush(stdout);}
    res = newkk_algn (s1p, s2p, s1_len,s2_len,gap_open,tc, mp, spd);
    CAMLreturn(Val_int(res));
};


void trivial_backtrace(const seqt s1, const seqt s2, seqt alis1, seqt alis2, 
        int len1, int len2, const cmt c) 
{
    int debug = 0;
    if(debug) {
       printf("trivial backtrace,len1=%d,len2=%d\n",len1,len2); fflush(stdout); }
    int i;
    int add1, add2;
    for (i=0;i<len1;i++)
    {
        add1 = seq_get(s1,i);
        add2 = cm_get_gap(c);
        my_prepend(alis1,add1);
        my_prepend(alis2,add2);
    }
    for (i=0;i<len2;i++)
    {
        add1 = cm_get_gap(c);
        add2 = seq_get(s2,i);
        my_prepend(alis1,add1);
        my_prepend(alis2,add2);
    }
    if(debug) {printf ("end of trivial backtrace\n"); fflush(stdout);}
}

void newkk_follow_deletion (const seqt s1, seqt alis1, seqt alis2, const cmt c, int * i)
{
    int debug = 0;
   int add1 = 0,add2 = 0;
    if(debug) { printf ("Delete,"); fflush(stdout);}
        add1 = seq_get(s1,*i);
        add2 = cm_get_gap(c);
        my_prepend(alis1,/*NOGAP &*/ add1);
        my_prepend(alis2,add2);
        *i=*i-1;

}

void newkk_follow_insertion (const seqt s2, seqt alis1, seqt alis2,const cmt c, int * j)
{
    int debug = 0;
   int add1 = 0,add2 = 0;
    if(debug) { printf ("Insert,"); fflush(stdout);}
        add2 = seq_get(s2,*j);
        add1 = cm_get_gap(c);
        my_prepend(alis1,add1);
        my_prepend(alis2,/*NOGAP &*/ add2);
        *j=*j-1;
}

void newkk_follow_deletion_or_insertion (int swaped,int has_insert, int has_delete,const seqt s1, const seqt s2, seqt alis1, seqt alis2,const cmt c, int * i, int * j)
{
    if (!swaped) {
        if (has_delete)//(dir==DO_DELETE)
            newkk_follow_deletion(s1,alis1,alis2,c,i);
        else{
            assert(has_insert);//(dir==DO_INSERT);
            newkk_follow_insertion(s2,alis1,alis2,c,j);
        }
    } else {
        if (has_insert)//(dir==DO_INSERT)
            newkk_follow_insertion(s2,alis1,alis2,c,j);
        else {
            assert(has_delete);//(dir==DO_DELETE);
            newkk_follow_deletion(s1,alis1,alis2,c,i);
        }
    }
}

void newkk_follow_deletion_or_insertion_affine (enum MODE * mode, int insert_delete_are_the_same, int has_end_insert, int has_end_delete, int swaped,int has_insert, int has_delete,const seqt s1, const seqt s2, seqt alis1, seqt alis2,const cmt c, int * i, int * j)
{
    if (!swaped) {
        if (has_delete)//(dir==DO_DELETE)
        {
            newkk_follow_deletion(s1,alis1,alis2,c,i);
            if (has_end_delete||insert_delete_are_the_same) *mode = m_todo;
        }
        else
        {
            assert(has_insert);//(dir==DO_INSERT);
            newkk_follow_insertion(s2,alis1,alis2,c,j);
            if (has_end_insert||insert_delete_are_the_same) *mode = m_todo;
        }
    } else {
        if (has_insert)//(dir==DO_INSERT) 
        {
            newkk_follow_insertion(s2,alis1,alis2,c,j);
            if (has_end_insert||insert_delete_are_the_same) *mode = m_todo;
        }
        else 
        {
            assert(has_delete);//(dir==DO_DELETE);
            newkk_follow_deletion(s1,alis1,alis2,c,i);
            if (has_end_delete||insert_delete_are_the_same) *mode = m_todo;
        }
    }
}


//backtrace function for non-affine
void backtrace (const seqt s1, const seqt s2, seqt alis1, seqt alis2, 
                int len1, int len2, const cmt c, const newkkmat_p m,int swaped)
{
    int add1 = 0,add2 = 0, cost, affine;
    int whichdiag, idx_in_my_diag, at_leftborder, at_rightborder;
    DIRECTION_MATRIX dir, gapnum1, gapnum2; 
    int affP=INT_MAX/2, affQ=INT_MAX/2, affED=INT_MAX/2, affCD=INT_MAX/2;
    int i=len1-1, j=len2-1;
    int has_insert, has_delete;
    assert(len1<=len2);
    if (len2 > MUCH_LONGER * len1) {
        trivial_backtrace(s1,s2,alis1,alis2,len1,len2,c);
    } else {
        while(i>=0&&j>=0) {
            whichdiag = 0, idx_in_my_diag=0, at_leftborder=0, at_rightborder=0;
            cost = 0; dir = 0; gapnum1 = 0; gapnum2 = 0;
            get_idx(-1,i,j,m,&whichdiag,&idx_in_my_diag,&at_leftborder,&at_rightborder,NULL);
            get_ukkcost(whichdiag,idx_in_my_diag,m, &cost, &dir, &gapnum1, &gapnum2,&affP, &affQ, &affED, &affCD, affine);
            has_insert = has_flag(dir,DO_INSERT);
            has_delete = has_flag(dir,DO_DELETE);
            if (dir==START) {
                if ( (i!=0) ||(j!=0) ) 
                    failwith ("ERROR: newkkonen.backtrace,we are expecting i==0 && j==0");
                add1 = add2 = cm_get_gap(c);
                my_prepend(alis1,add1);
                my_prepend(alis2,add2);
                i--; j--;
            } else { //follow what algn.c is doing -- for affine, we favor 
                     //insert/delete over algn; otherwise we pick algn first
                if (has_flag(dir,ALIGN_TO_ALIGN))//(dir==DO_ALIGN)
                {
                    add1 = seq_get(s1,i);
                    add2 = seq_get(s2,j);
                    my_prepend(alis1,add1);
                    my_prepend(alis2,add2);
                    i--; j--;
                } else {
                    newkk_follow_deletion_or_insertion (swaped,has_insert,has_delete,s1,s2,alis1,alis2,c,&i,&j);
                }
            }
       }//end of while(i>=0&&j>=0)
       //at the end of while loop, we should reach (-1,-1) from (0,0)
       assert( (i== -1)&&(j==-1));
   }//end of non-trivial case
}

void copy_dir_to_mode (enum MODE * mode, int has_insert, int has_delete, int has_algn, int has_diag, int swaped)
{
    if (!swaped) {
        if (has_delete)  *mode = m_delete;
        else if (has_insert)  *mode = m_insert;
        else if (has_diag) *mode = m_diagonal;
        else if (has_algn) *mode = m_align;
        else 
            failwith ("newkkonen.c,copy_dir_to_mode,invalid dir");
    }
    else {
        if (has_insert) *mode = m_insert;
        else if (has_delete) *mode = m_delete;
        else if (has_diag) *mode = m_diagonal;
        else if (has_algn) *mode = m_align;
        else
            failwith ("newkkonen.c,copy_dir_to_mode,invalid dir");
    }

}

void backtrace_affine (const seqt s1, const seqt s2, seqt alis1, seqt alis2, 
        int len1, int len2, const cmt c, const newkkmat_p m,int swaped)
{
   int debug = 0;
   int affine = 1;
   if(debug) {
       printf("\nnewkkonen backtrace affine,len1=%d,len2=%d,swaped=%d\n",len1,len2,swaped); fflush(stdout);       }
   assert(len1<=len2);
   enum MODE mode = m_todo;
   if (len2 > MUCH_LONGER * len1) 
   {
       trivial_backtrace(s1,s2,alis1,alis2,len1,len2,c);
   }//end of trivial case
   else {
   int add1 = 0,add2 = 0;
   int whichdiag, idx_in_my_diag, at_leftborder, at_rightborder;
   int cost; DIRECTION_MATRIX dir; DIRECTION_MATRIX gapnum1; DIRECTION_MATRIX gapnum2; 
    //set affP and affQ anyway
   int affP=INT_MAX/2, affQ=INT_MAX/2, affED=INT_MAX/2, affCD=INT_MAX/2;
   int i=len1-1, j=len2-1;
   int mode_delete=0, mode_insert=0;
   int has_diag=0, has_insert=0, has_delete=0, has_algn=0;
   int has_end_delete=0, has_end_insert=0;
   int insert_delete_are_the_same = 0;
   while(i>=0&&j>=0)
   {
       //if (i==35&&j==62) debug = 1;
        whichdiag = 0, idx_in_my_diag=0, at_leftborder=0, at_rightborder=0;
        cost = 0; dir = 0; gapnum1 = 0; gapnum2 = 0;
        get_idx(-1,i,j,m,&whichdiag,&idx_in_my_diag,&at_leftborder,&at_rightborder,NULL);
        get_ukkcost(whichdiag,idx_in_my_diag,m, &cost, &dir, &gapnum1, &gapnum2,&affP, &affQ, &affED, &affCD, affine);
        if(debug) {printf ("i=%d,j=%d :",i,j); print_direction(dir); fflush(stdout); }
        if (dir==START) 
        {
            if(debug) { printf ("Start\n"); fflush(stdout);}
            if ( (i!=0) ||(j!=0) ) 
                failwith ("ERROR: newkkonen.backtrace,we are expecting i==0 && j==0");
            add1 = add2 = cm_get_gap(c);
            my_prepend(alis1,add1);
            my_prepend(alis2,add2);
            i--; j--;
        }
        else {//follow what algn.c is doing -- for affine, we favor insert/delete over algn
            if (debug) print_mode (mode);
            if (mode == m_delete) mode_delete = 1; else mode_delete = 0;
            if (mode == m_insert) mode_insert = 1; else mode_insert = 0;
            if (mode == m_todo)
            {
                 has_insert = has_flag(dir,DO_INSERT);
                 has_delete = has_flag(dir,DO_DELETE);
                 has_algn = has_flag(dir,ALIGN_TO_ALIGN);
                 has_diag = has_flag(dir,DO_DIAG);
                 copy_dir_to_mode (&mode,has_insert,has_delete,has_algn,has_diag,swaped);
            }
            else if (mode_delete||mode_insert)
            {
                if (mode_delete) {
                    has_delete=1;
                } 
                else 
                    has_delete=0;
                if (mode_insert) { 
                    has_insert=1; 
                } 
                else 
                    has_insert=0;
                has_end_delete = has_flag(dir,END_DELETE);
                has_end_insert = has_flag(dir,END_INSERT);
                insert_delete_are_the_same = has_flag(dir,INSERT_EQ_DELETE);
                if(debug) printf("mode_insert=%d,mode_delete=%d,has_end_insert=%d,has_end_delete=%d,has_insert=%d,has_delete=%d\n",mode_insert,mode_delete,has_end_insert,has_end_delete,has_insert,has_delete);
                newkk_follow_deletion_or_insertion_affine (&mode,insert_delete_are_the_same,has_end_insert,has_end_delete,swaped,has_insert,has_delete,s1,s2,alis1,alis2,c,&i,&j);
                 
            }
            else if (mode == m_diagonal)
            {
                if (has_flag(dir,END_DIAG)) mode = m_todo;
                add1 = seq_get(s1,i);
                add2 = seq_get(s2,j);
                my_prepend(alis1,add1);
                my_prepend(alis2,add2);
                i--; j--;
            }
            else // (mode == m_align)
            {
                assert(mode == m_align);
                if(debug) { printf ("Align\n");fflush(stdout);}
                add1 = seq_get(s1,i);
                add2 = seq_get(s2,j);
                my_prepend(alis1,NOGAP & add1);
                my_prepend(alis2,NOGAP & add2);
                i--; j--;
                mode = m_todo;
            }
       }//end of if (dir==START)
   }//end of while(i>=0&&j>=0)
   }//end of non-trivial case
   if(debug) { printf ("end of backtrace\n"); fflush(stdout);}
};



value
newkkonen_CAML_backtrace (value s1, value s2, value s1p, value s2p, value c, value a, value swaped) {
    CAMLparam5(s1, s2, s1p, s2p, c);
    CAMLxparam2(a,swaped);
    seqt ss1, ss2, ss1p, ss2p;
    newkkmat_p mp;
    cmt cc;
    mp = Newkkmat_struct(a);
    Seq_custom_val(ss1,s1);
    Seq_custom_val(ss2,s2);
    Seq_custom_val(ss1p,s1p);
    Seq_custom_val(ss2p,s2p);
    cc = Cost_matrix_struct(c);
    int len1 = seq_get_len(ss1);
    int len2 = seq_get_len(ss2);
    int spd = Int_val(swaped);
    backtrace (ss1, ss2, ss1p, ss2p, len1, len2, cc, mp, spd);
    CAMLreturn(Val_unit);
};

value
newkkonen_CAML_backtrace_affine (value s1, value s2, value s1p, value s2p, value c, value a, value swaped) {
    CAMLparam5(s1, s2, s1p, s2p, c);
    CAMLxparam2(a,swaped);
    seqt ss1, ss2, ss1p, ss2p;
    newkkmat_p mp;
    cmt cc;
    mp = Newkkmat_struct(a);
    Seq_custom_val(ss1,s1);
    Seq_custom_val(ss2,s2);
    Seq_custom_val(ss1p,s1p);
    Seq_custom_val(ss2p,s2p);
    cc = Cost_matrix_struct(c);
    int spd = Int_val(swaped);
    backtrace_affine (ss1, ss2, ss1p, ss2p, seq_get_len(ss1), seq_get_len(ss2), cc, mp,spd);
    CAMLreturn(Val_unit);
};


value 
newkkonen_CAML_backtrace_bc (value *argv, int argn) {
    return (newkkonen_CAML_backtrace (argv[0], argv[1], argv[2], argv[3], argv[4], argv[5], argv[6]) );
};

value 
newkkonen_CAML_backtrace_affine_bc (value *argv, int argn) {
    return (newkkonen_CAML_backtrace (argv[0], argv[1], argv[2], argv[3], argv[4], argv[5], argv[6]) );
};


void
newkkmat_CAML_free (value m) {
    int debug = 0;
    if (debug) { printf("newkkmat_CAML_free\n"); fflush(stdout); }
    newkkmat_p tmp;
    tmp = Newkkmat_struct(m);
    int i = 0; 
    //reset pointers to memory pool to NULL
    int size = tmp->diag_size;
    for (i=0;i<size;i++)
    {
        ukkdiag_p t2 = tmp->diagonal;
        t2 = t2 + i ;
        t2->costarr = NULL;
        t2->dirarr = NULL;
        t2->gapnumarr1 = NULL;
        t2->gapnumarr2 = NULL;
        //free( t2->costarr );
        //free( t2->dirarr );
        //free( t2->gapnumarr );
    }
    tmp->diag_size=0;
    tmp->diag_size_in_use=0;
    tmp->k=0;
    tmp->total_len=0;
    tmp->total_len_in_use=0;
    tmp->total_len_affine=0;
    tmp->baseband=0;
    //free cost_dir array
    free (tmp->diagonal);
    free (tmp->diagonal_size_arr);
    //free memory pool
    free(tmp->pool_cost);
    free(tmp->pool_dir);
    free(tmp->pool_gapnum);
    if (tmp->pool_affP!=NULL) free(tmp->pool_affP);
    if (tmp->pool_affQ!=NULL) free(tmp->pool_affQ);
    if (tmp->pool_affED!=NULL) free(tmp->pool_affED);
    if (tmp->pool_affCD!=NULL) free(tmp->pool_affCD);
    return;
};


static struct custom_operations newkk_alignment_matrix = {
    "http://www.amnh.org/poy/",
    &newkkmat_CAML_free,
    custom_compare_default,
    custom_hash_default,
    custom_serialize_default,
    custom_deserialize_default
};


value 
newkkonen_CAML_create_general (value a) {
    CAMLparam1(a);
    CAMLlocal1(res);
    newkkmat_p m;
    res = 
        alloc_custom (&newkk_alignment_matrix, sizeof(struct newkkmat), 1, 1000000);
    m = Newkkmat_struct(res);
    m->total_len = m->total_len_in_use = m->k = m->diag_size = m->diag_size_in_use = m->baseband = -1;
    m->diagonal = NULL;
    m->diagonal_size_arr = NULL;
    m->pool_cost = NULL;
    m->pool_dir = NULL;
    m->pool_gapnum = NULL;
    m->total_len_affine = -1;
    m->pool_affP = NULL;
    m->pool_affQ = NULL;
    CAMLreturn(res);
}

value 
newkkonen_CAML_initialize (value unit) {
    CAMLparam1(unit);
    caml_register_custom_operations (&newkk_alignment_matrix);
    CAMLreturn(Val_unit);
}



