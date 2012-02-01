/* POY 5.0 Alpha. A phylogenetic analysis program using Dynamic Homologies.   *\
(* Copyright (C) 2011 Andrés Varón, Lin Hong, Nicholas Lucaroni, Ward Wheeler,*)
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
#include <stdio.h>
#include <limits.h>
#include <string.h>
#define NDEBUG 1
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
//alignment
#define START 0
#define DO_ALIGN 1
#define DO_DELETE 2
#define DO_INSERT 4


#define my_prepend(a,b) assert (a->cap > a->len); \
    (a)->begin = (((a)->begin) - 1); \
    ((a)->len = 1 + (a)->len); *((a)->begin) = b
#define my_get(a,b) ((a)->begin)[b]

#define MAX(a,b) ( ((a)>(b))? (a):(b) )
#define MIN(a,b) ( ((a)<(b))? (a):(b) )

#ifdef _WIN32
__inline int
#else
inline int
#endif
get_k (newkkmat_p mp)
{
    assert(mp != NULL);
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
reset_mat_to_0 (newkkmat_p m)
{
    int debug = 0;
    DIRECTION_MATRIX v = (DIRECTION_MATRIX) (-1);
    MAT_SIZE total_len = m->total_len;
    if (debug) printf ("reset matrix to 0, total_len=%d\n",total_len);
    //reset cost
    memset(m->pool_cost,(-1),total_len*sizeof(int));
    //these two are not necessary.
    memset(m->pool_dir,v,total_len*sizeof(DIRECTION_MATRIX));
    memset(m->pool_gapnum,v,total_len*sizeof(DIRECTION_MATRIX));
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
    MAT_SIZE baseband = m->baseband;
    MAT_SIZE oldsize = m->diag_size_in_use;
    //get diagonal array size, and set diag_size_in_use
    MAT_SIZE newsize = oldsize+(newk-oldk)*2;
    m->k = newk;
    m->diag_size_in_use = newsize; 
    if (debug) 
    { printf ("newkkonen.c, expand_mat,newk=%d(oldk=%d),newsize of diagarr=%d(old=%d),baseband=%d,affine=%d\n",
            newk,oldk,newsize,oldsize,baseband,affine); fflush(stdout); }
    //set expand sign
    if (debug) { printf ("newsize=%d <=> m->diag_size=%d\n", newsize, m->diag_size); fflush(stdout);}
    if (newsize <= m->diag_size) expand_diag_size = 0;
    else 
    {
        if (debug) {printf ("we NEED to expand diagonal_size_arr\n"); fflush(stdout); }
        expand_diag_size = 1;
        m->diagonal_size_arr = (int *)realloc(m->diagonal_size_arr,newsize*sizeof(int));
        m->diagonal = (ukkdiag_p ) realloc (m->diagonal,newsize*sizeof(struct cost_dir));
    }
    assert(m->diagonal!=NULL);
    assert(m->diagonal_size_arr!=NULL);
    for (i=oldsize;i<newsize;i=i+2)
    {
        (m->diagonal_size_arr)[i] =  (m->diagonal_size_arr)[i-1] - 1;
        (m->diagonal_size_arr)[i+1] =  (m->diagonal_size_arr)[i-1] - 1;
    }
    
    for (i=oldsize;i<newsize;i++) //update total_len_in_use
    {
        diaglen = (m->diagonal_size_arr)[i];
        m->total_len_in_use += diaglen;
        //if (affine) m->total_len_affine += diaglen;
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
        if (debug) { printf ("we NEED to expand diagonal\n"); fflush(stdout);}
        expand_diagonal = 1;
        MAT_SIZE len_we_need = m->total_len_in_use;
        m->pool_cost = realloc (m->pool_cost,len_we_need*sizeof(int));
        m->pool_dir = realloc (m->pool_dir,len_we_need*sizeof(DIRECTION_MATRIX));
        m->pool_gapnum = realloc (m->pool_gapnum,len_we_need*sizeof(DIRECTION_MATRIX));
        if (affine) {
            m->pool_affP = realloc (m->pool_affP,len_we_need*sizeof(int));
            m->pool_affQ = realloc (m->pool_affQ,len_we_need*sizeof(int));
        }
    }
        ukkdiag_p emptydiag;
        int * thiscost= m->pool_cost;
        //define these two pointer anyway if (affine) {
        int * thisaffP = m->pool_affP;
        int * thisaffQ = m->pool_affQ;//}
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
            emptydiag -> gapnumarr = thisgap;
            if (affine) {
            emptydiag -> affParr = thisaffP;
            emptydiag -> affQarr = thisaffQ;
            thisaffP += diaglen;
            thisaffQ += diaglen;
            }
            thiscost += diaglen;
            thisdir += diaglen;
            thisgap += diaglen;
            
            emptydiag++;
        }
        for (i=oldsize;i<newsize;i++)
        {
            diaglen = emptydiag->len;
            assert(diaglen>0);
            //printf ("set diag#%d(%d);",i,diaglen); fflush(stdout);
            emptydiag -> costarr = thiscost;
            emptydiag -> dirarr = thisdir;
            emptydiag -> gapnumarr = thisgap;
            if (affine) {
            emptydiag -> affParr = thisaffP;
            emptydiag -> affQarr = thisaffQ;}
            if (i<newsize-1) //don't move out of range
            { 
            emptydiag ++;
            thiscost += diaglen ;
            thisdir += diaglen;
            thisgap += diaglen; 
            if (affine) {
                thisaffP += diaglen;
            thisaffQ += diaglen; }
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

void
get_idx (int i,int j,newkkmat_p m, int * whichdiag,int * idx_in_my_diag, int * at_leftborder, int * at_rightborder)
{
    int debug = 0;
    int currentK = m->k;
    assert(currentK>=0);
    //though baseband is defined as MAT_SIZE in newkkonen.h, which might be long int. but actually we don't need that much, int is fine. if baseband itself is bigger than max_int, the memory requirement will blow off our computer anyway. 
    int bb = m->baseband;
    assert(bb>0); 
    if (debug) { printf("get_idx,i=%d,j=%d.bb=%d.currentk=%d,",i,j,bb,currentK); fflush(stdout);}
    int k_ij= i-j, k_ji = j-i;
    if(i>j) 
    { 
        *whichdiag = 2*(i-j)+bb-1 ; 
        if (debug) { printf("k_ij=%d,whichdiag=%d,",k_ij,*whichdiag); fflush(stdout);}
    }
    else if ( (i<=j)&&( (j-i)<bb ) )
    { 
        *whichdiag = j-i; 
        if (debug) {printf("whichdiag=%d,",*whichdiag); fflush(stdout);}
    }
    else {
        *whichdiag = bb + 2*(j-i+1 - bb) -1 -1;
        if (debug) {printf("k_ji=%d,whichdiag=%d,",k_ji,*whichdiag); fflush(stdout); }
    }
    assert(*whichdiag >=0);
    assert(*whichdiag < m->size);
    if ( (k_ji-bb+1)==currentK ) { *at_rightborder=1;}
        else *at_rightborder=0;
    if (k_ij==currentK) {  *at_leftborder = 1; } 
        else *at_leftborder = 0;
    int idx; 
    if(i<=j) idx=i;  
    else idx=j;  
    *idx_in_my_diag = idx;
    if (debug) printf("idx=%d\n",idx);
};


void 
get_ukkcost (int whichdiag, int idx_in_my_diag, newkkmat_p m, int * cost, DIRECTION_MATRIX * dir, DIRECTION_MATRIX *  max_gapnum, int * affP, int * affQ, int affine )
{
    int debug = 0;
    ukkdiag_p thisdiag = m->diagonal;
    assert(whichdiag < m->size);
    thisdiag += whichdiag;
    int * thiscostarr = thisdiag->costarr;
    DIRECTION_MATRIX * thisdirarr = thisdiag->dirarr;
    DIRECTION_MATRIX * thisgapnumarr = thisdiag->gapnumarr;
    assert(idx_in_my_diag < thisdiag->len);
    *cost = thiscostarr[idx_in_my_diag];
    *dir = thisdirarr[idx_in_my_diag];
    *max_gapnum = thisgapnumarr[idx_in_my_diag];
    if (affine) {
    int * thisParr = thisdiag->affParr;
    int * thisQarr = thisdiag->affQarr;
    *affP = thisParr[idx_in_my_diag];
    *affQ = thisQarr[idx_in_my_diag];}
    if (debug) { printf("get ukkcost,diag#.%d,idx#.%d,cost=%d,dir=%d,mapgapnum=%d,affP=%d,affQ=%d,affine=%d\n",whichdiag,idx_in_my_diag,*cost,*dir,*max_gapnum,affine); fflush(stdout);}
};

void set_ukkcost (int whichdiag, int idx_in_my_diag, newkkmat_p m, int cost, DIRECTION_MATRIX dir,DIRECTION_MATRIX  max_gapnum, int affP, int affQ, int affine)
{
    int debug = 0;
    ukkdiag_p thisdiag = m->diagonal;
    if ((idx_in_my_diag >= thisdiag->len) || (idx_in_my_diag<0)) { debug=1; }
    if (debug) { printf ("set ukkcost, diag#.%d,idx#.%d,cost=%d,dir=%d,mapgapnum=%d,affP=%d,affQ=%d\n",whichdiag,idx_in_my_diag,cost,dir,max_gapnum,affP,affQ); fflush(stdout);}
    thisdiag += whichdiag;
    int * thiscostarr = thisdiag->costarr;
    DIRECTION_MATRIX * thisdirarr = thisdiag->dirarr;
    DIRECTION_MATRIX * thisgapnumarr = thisdiag->gapnumarr;
    assert(idx_in_my_diag < thisdiag->len);
    thiscostarr[idx_in_my_diag] = cost;
    thisdirarr[idx_in_my_diag] = dir;
    thisgapnumarr[idx_in_my_diag] = max_gapnum;
    if (affine) {
        int * thisParr = thisdiag->affParr;
        int * thisQarr = thisdiag->affQarr;
    thisParr[idx_in_my_diag] = affP;
    thisQarr[idx_in_my_diag] = affQ;}

};

void sanity_check2 (int x)
{
    if ((x==1)||(x==2)||(x==4)||(x==8)||(x==16)) {}
    else {
        printf ("WARNING: we have a x=%d, not a dna code\n", x);
        fflush(stdout);
    }
};

void get_cmcost (const cmt c, int a, int b, int * res)
{
    //if (a==b) { *res = 0; } else { *res = 2 ; }
    //sanity_check2(a);
    //sanity_check2(b);
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

void update_internal_cell (const seqt s1, const seqt s2,newkkmat_p m, const cmt c, int i, int j, int newk,int go)
{
    int debug = 0;
    if (debug) printf (" update internal cell: ");
    int gapcode = cm_get_gap (c);
    int addcost = 0;
    int costL=0, costR=0, costM=0;
    int thisP=INT_MAX/2, thisQ=INT_MAX/2;
//set affine to 1 is gap opening is set
    int affine=0;
    if (go>=0) affine=1;
    int realgo=0;
    if (go>0) realgo=go;
    //init affP and affQ anyway
    int affP=INT_MAX/2, affQ=INT_MAX/2;
    int costfromL=0; DIRECTION_MATRIX dirL; DIRECTION_MATRIX gapnum_fromL=0;
    if (j==0) {
        if (debug) printf ("j==0,set costfromL = costL = thisQ <- INT_MAX/2\n");
        costfromL = costL = thisQ = INT_MAX/2; }//leftmost column of matrix
    else {
    int whichdiagL=0, idx_in_my_diagL=0, at_leftborderL=0, at_rightborderL=0;
    get_idx(i,j-1,m,&whichdiagL,&idx_in_my_diagL,&at_leftborderL,&at_rightborderL);
    get_ukkcost(whichdiagL,idx_in_my_diagL, m, &costfromL, &dirL, &gapnum_fromL, &affP, &affQ,affine);
    get_cmcost (c,seq_get(s2,j),gapcode,&addcost);
    costL = MIN(costfromL + addcost + realgo, affQ + addcost);
    thisQ = costL;
    if (debug) printf ("costL <- min (costfromL<%d> + addcost<%d> + go <%d>,affQ<%d>+addcost<%d>)=%d\n",
    costfromL,addcost,realgo,affQ,addcost,costL);
    }
    int costfromR=0; DIRECTION_MATRIX dirR; DIRECTION_MATRIX gapnum_fromR=0;
    if (i==0) {
        if (debug) printf ("i==0,costfromR = costR = thisP <- INT_MAX/2\n");
        costfromR = costR = thisP = INT_MAX/2; }//first line of matrix
    else {
    int whichdiagR=0, idx_in_my_diagR=0, at_leftborderR=0, at_rightborderR=0;
    get_idx(i-1,j,m,&whichdiagR,&idx_in_my_diagR,&at_leftborderR,&at_rightborderR);
    get_ukkcost(whichdiagR,idx_in_my_diagR, m, &costfromR, &dirR, &gapnum_fromR, &affP, &affQ,affine); 
    get_cmcost (c,seq_get(s1,i),gapcode,&addcost);
    costR = MIN(costfromR + addcost + realgo, affP + addcost);
    thisP = costR;
    if (debug) printf ("costR <- min (costfromR<%d> + addcost<%d> + realgo <%d>,affP<%d>+addcost<%d>)=%d\n",
    costfromR,addcost,realgo,affP,addcost,costR);
    }
    int costfromM=0; DIRECTION_MATRIX dirM; DIRECTION_MATRIX gapnum_fromM=0;
    if (i*j==0) { //first line or leftmost column of matrix
        if (debug) printf ("i==0||j==0,set costfromM = costM <- INT_MAX/2\n");
        costfromM = costM = INT_MAX/2; 
        if (i+j==0) failwith ("newkkonen.c,update_internal_cell,we don't update (0,0) here\n");
    }
    else {
    int whichdiagM=0, idx_in_my_diagM=0, at_leftborderM=0, at_rightborderM=0;
    get_idx(i-1,j-1,m,&whichdiagM,&idx_in_my_diagM,&at_leftborderM,&at_rightborderM);
    get_ukkcost(whichdiagM,idx_in_my_diagM, m, &costfromM, &dirM, &gapnum_fromM, &affP, &affQ,affine); 
    get_cmcost (c,seq_get(s1,i),seq_get(s2,j),&addcost);
    costM = costfromM + addcost;
    if (debug) printf ("costM <- costfromM<%d> + addcost = %d\n",costfromM,addcost,costM);
    }
    int whichdiag=0, idx_in_my_diag=0, at_leftborder=0, at_rightborder=0;
    get_idx(i,j,m,&whichdiag,&idx_in_my_diag,&at_leftborder,&at_rightborder);
    if ((costL<=costR) && (costL<=costM)) 
    {
        set_ukkcost(whichdiag,idx_in_my_diag, m, costL, DO_INSERT, gapnum_fromL+1, thisP, thisQ, affine);
        if (debug) printf ("cost <---- costL\n");
    }
    else if ((costR<=costL) && (costR<=costM)) 
    {
        set_ukkcost(whichdiag,idx_in_my_diag, m, costR, DO_DELETE, gapnum_fromR+1, thisP, thisQ, affine);
        if (debug) printf ("cost <---- costR\n");
    }
    else
    {
        int newgapnum;
        if (costM==costfromM) newgapnum=gapnum_fromM;
        else newgapnum = gapnum_fromM + 1;
        set_ukkcost(whichdiag,idx_in_my_diag, m, costM, DO_ALIGN, newgapnum, thisP, thisQ, affine);
        if (debug) printf ("cost <---- costM\n");
    }

};

void update_left_border_cell (const seqt s1, const seqt s2,newkkmat_p m, const cmt c, int i, int j, int newk,int go)
{
    int debug = 0;
    if (debug) printf (" update left border cell, ");
    int gapcode = cm_get_gap (c);
    int addcost = 0;
    int thisP=INT_MAX/2, thisQ=INT_MAX/2;
    //set affine to 1 if gap opening is set
    int affine=0;
    if (go>=0) affine=1;
    int realgo=0;
    if (go>0) realgo=go;
    //set affP and affQ anyway
    int affP=INT_MAX/2, affQ=INT_MAX/2;
    int costR=0, costM=0;
    int costfromR=0; DIRECTION_MATRIX dirR; DIRECTION_MATRIX gapnum_fromR=0;
    if (i==0) { failwith ("newkkonen.c,update_left_border_cell,we don't update (0,0)\n"); }
    else {
    int whichdiagR=0, idx_in_my_diagR=0, at_leftborderR=0, at_rightborderR=0;
    get_idx(i-1,j,m,&whichdiagR,&idx_in_my_diagR,&at_leftborderR,&at_rightborderR);
    get_ukkcost(whichdiagR,idx_in_my_diagR, m, &costfromR, &dirR, &gapnum_fromR, &affP, &affQ, affine); 
    get_cmcost (c,seq_get(s1,i),gapcode,&addcost);
    //int gap_opening = 0;if (dirR==DO_DELETE) { gap_opening = 0; }else { gap_opening = go; }
    //costR = costfromR + addcost;// + gap_opening;
    costR = MIN(costfromR + addcost + realgo, affP + addcost);
    thisP = costR;
    }
    int costfromM=0; DIRECTION_MATRIX dirM; DIRECTION_MATRIX gapnum_fromM=0;
    if (j==0) { costfromM = costM = INT_MAX/2; }//leftmost column of matrix
    else {
    int whichdiagM=0, idx_in_my_diagM=0, at_leftborderM=0, at_rightborderM=0;
    get_idx(i-1,j-1,m,&whichdiagM,&idx_in_my_diagM,&at_leftborderM,&at_rightborderM);
    get_ukkcost(whichdiagM,idx_in_my_diagM, m, &costfromM, &dirM, &gapnum_fromM, &affP, &affQ, affine); 
    get_cmcost (c,seq_get(s1,i),seq_get(s2,j),&addcost);
    costM = costfromM + addcost;
    }
    int whichdiag=0, idx_in_my_diag=0, at_leftborder=0, at_rightborder=0;
    get_idx(i,j,m,&whichdiag,&idx_in_my_diag,&at_leftborder,&at_rightborder);
    assert(at_leftborder);
    if (costR<=costM) 
    {
        set_ukkcost(whichdiag,idx_in_my_diag, m, costR, DO_DELETE, gapnum_fromR+1,thisP, thisQ, affine);
    }
    else
    {
        int newgapnum;
        if (costM==costfromM) newgapnum=gapnum_fromM;
        else newgapnum = gapnum_fromM + 1;
        set_ukkcost(whichdiag,idx_in_my_diag, m, costM, DO_ALIGN, newgapnum,thisP, thisQ, affine);
    } 
};

void update_right_border_cell (const seqt s1, const seqt s2,newkkmat_p m, const cmt c, int i, int j, int newk,int go)
{
    int debug = 0;
    if (debug) printf (" update right border cell ,");
    int gapcode = cm_get_gap (c);
    int addcost = 0;
    int thisP=INT_MAX/2, thisQ=INT_MAX/2;
    //set affine to 1 if gap opening is set
    int affine=0;
    if (go>=0) affine=1;
    int realgo=0;
    if (go>0) realgo=go;
    //set affP and affQ anyway
    int affP=INT_MAX/2, affQ=INT_MAX/2;
    int costL=0, costM=0;
    int costfromL=0; DIRECTION_MATRIX dirL; DIRECTION_MATRIX gapnum_fromL=0;
    if (j==0) { failwith ("newkkonen.c update_right_border_cell, cannot be at left border as well\n"); }
    else {
    int whichdiagL=0, idx_in_my_diagL=0, at_leftborderL=0, at_rightborderL=0;
    get_idx(i,j-1,m,&whichdiagL,&idx_in_my_diagL,&at_leftborderL,&at_rightborderL);
    get_ukkcost(whichdiagL,idx_in_my_diagL, m, &costfromL, &dirL, &gapnum_fromL,&affP, &affQ, affine);
    get_cmcost (c,seq_get(s2,j),gapcode,&addcost);
    //int gap_opening = 0;if (dirL==DO_INSERT) { gap_opening = 0; }else { gap_opening = go; }
    //costL = costfromL + addcost;// + gap_opening;
    costL = MIN(costfromL + addcost + realgo, affQ + addcost);
    thisQ = costL;
    }
    int costfromM=0; DIRECTION_MATRIX dirM; DIRECTION_MATRIX gapnum_fromM=0;
    if (i==0) { costfromM = costM = INT_MAX/2; }//top line of matrix
    else {
    int whichdiagM=0, idx_in_my_diagM=0, at_leftborderM=0, at_rightborderM=0;
    get_idx(i-1,j-1,m,&whichdiagM,&idx_in_my_diagM,&at_leftborderM,&at_rightborderM);
    get_ukkcost(whichdiagM,idx_in_my_diagM, m, &costfromM, &dirM, &gapnum_fromM,&affP, &affQ, affine); 
    get_cmcost (c,seq_get(s1,i),seq_get(s2,j),&addcost);
    costM = costfromM + addcost;
    }
    int whichdiag=0, idx_in_my_diag=0, at_leftborder=0, at_rightborder=0;
    get_idx(i,j,m,&whichdiag,&idx_in_my_diag,&at_leftborder,&at_rightborder);
    assert(at_rightborder);
    if (costL<=costM)
    {
        if (debug) {printf("costL=%d <= costM=%d\n",costL,costM); fflush(stdout);}
        set_ukkcost(whichdiag,idx_in_my_diag, m, costL, DO_INSERT, gapnum_fromL+1, thisP, thisQ, affine);
    }
    else
    {
        if (debug) {printf("costL=%d > costM=%d\n",costL,costM); fflush(stdout);}
        int newgapnum;
        if (costM==costfromM) newgapnum=gapnum_fromM;
        else newgapnum = gapnum_fromM + 1;
        set_ukkcost(whichdiag,idx_in_my_diag, m, costM, DO_ALIGN, newgapnum, thisP, thisQ, affine);
    }

};

void update_central_diagonal_cell (const seqt s1, const seqt s2,newkkmat_p m, const cmt c, int i, int j, int newk, int go)
{
    int debug = 0;
    int costM=0;
    int addcost = 0;
    int thisP=INT_MAX/2, thisQ=INT_MAX/2;
     //set affine to 1 if gap opening is set
    int affine=0;
    if (go>=0) affine=1;
    //set affP and affQ anyway
    int affP=INT_MAX/2, affQ=INT_MAX/2;
    int costfromM=0; 
    DIRECTION_MATRIX dirM; DIRECTION_MATRIX gapnum_fromM=0;
    if ((i==0)||(j==0)) { failwith ("newkkonen.c, update_central_diagonal_cell, we don't need update (0,0)\n"); }
    else {
    int whichdiagM=0, idx_in_my_diagM=0, at_leftborderM=0, at_rightborderM=0;
    get_idx(i-1,j-1,m,&whichdiagM,&idx_in_my_diagM,&at_leftborderM,&at_rightborderM);
    get_ukkcost(whichdiagM,idx_in_my_diagM, m, &costfromM, &dirM, &gapnum_fromM, &affP, &affQ, affine); 
    get_cmcost (c,seq_get(s1,i),seq_get(s2,j),&addcost);
    costM = costfromM + addcost;
    }
    int whichdiag=0, idx_in_my_diag=0, at_leftborder=0, at_rightborder=0;
    get_idx(i,j,m,&whichdiag,&idx_in_my_diag,&at_leftborder,&at_rightborder);
    assert(at_rightborder);
    assert(at_leftborder);
    int newgapnum;
    if (costM==costfromM) newgapnum=gapnum_fromM;
    else newgapnum = gapnum_fromM + 1;
    set_ukkcost(whichdiag,idx_in_my_diag, m, costM, DO_ALIGN, newgapnum, thisP, thisQ, affine);
    if (debug) printf("update_central_diagonal,cost<-costM=%d\n",costM);
};

int update_a_cell (const seqt s1, const seqt s2,newkkmat_p m, const cmt c, int i, int j, int newk, int go)
{
    int debug = 0;
    int whichdiag=0, idx_in_my_diag=0, at_leftborder=0, at_rightborder=0;
    get_idx(i,j,m,&whichdiag,&idx_in_my_diag,&at_leftborder,&at_rightborder);
    //define affine
    int affine=0;
    //get old cost first
    int oldcost;
    DIRECTION_MATRIX olddir, oldgapnum;
    int oldaffP, oldaffQ;
    get_ukkcost(whichdiag,idx_in_my_diag, m, &oldcost, &olddir, &oldgapnum, &oldaffP, &oldaffQ, affine);
    if ((whichdiag<0)||(whichdiag>= m->diag_size_in_use)) {debug=1;}
    if (debug) { printf("update_a_cell,(%d,%d),diag#.%d,idx#.%d,leftB=%d,rightB=%d,oldcost=%d\n",i,j,whichdiag,idx_in_my_diag,at_leftborder,at_rightborder,oldcost); fflush(stdout);}
    if ( i==0 && j==0 ) return 0;//we don't update (0,0)
    else 
    {
        if (at_leftborder&&at_rightborder) 
            update_central_diagonal_cell (s1, s2, m, c,i, j, newk, go);
        else if (at_leftborder)
            update_left_border_cell (s1, s2, m, c,i, j, newk, go);
        else if (at_rightborder)
            update_right_border_cell (s1, s2, m, c,i, j, newk, go);
        else
            update_internal_cell (s1, s2, m, c,i, j, newk, go);
        int newcost;
        DIRECTION_MATRIX newdir, newgapnum;
        int newaffP, newaffQ;
        get_ukkcost(whichdiag,idx_in_my_diag, m, &newcost, &newdir, &newgapnum, &newaffP, &newaffQ, affine);
        if (debug) {printf("newcost=%d,oldcost=%d\n",newcost,oldcost); fflush(stdout);}
        if (newcost==oldcost) return 0;
        else return 1;
    }
};


//function for speeding up start

/*
 * return true if (i,j) is in the new left band of newk, including old left border of oldk*/
int in_new_left_band (int i, int j,int newk,int oldk)
{
    if ( (i>=j)&&((i-j)<=newk)&&((i-j)>=oldk) ) return 1;
    else return 0;
};


void update_a_row (const seqt s1, const seqt s2,newkkmat_p m, const cmt c,int i,int startj,int endj,int last_updated_j,q_t prevq,q_t thisq,int newk,int oldk, int go)
{
    int debug = 0;
    int next_j;//which cell we should update at this run
    int follow_prevq;//if we update this cell by previous queue
    QUEUE_DATA msb=0;
    if (newk == oldk)//first run, update all cells
    {
        if (debug) printf("update_a_row,first run,update all,i=%d\n",i);
        next_j = last_updated_j + 1;
        if (next_j<=endj) {
             update_a_cell (s1,s2,m,c,i,next_j,newk,go);
             update_a_row(s1,s2,m,c,i,startj,endj,next_j,prevq,thisq,newk,oldk,go);
        }
        else {}//done with this row
    }
    else
    {
        if (debug) 
        {printf("update a row,i=%d,(startj,endj)=(%d,%d),last_updated_j=%d\n",i,startj,endj,last_updated_j); fflush(stdout);}
        if (last_updated_j == (startj-1) )
        {//this is the first move on current row
            next_j = startj; follow_prevq=0;
        }
        else if (last_updated_j == (-1) ) 
        {//no change in previous cell of current row, check preq of the row above us
            if (is_emptyqueue(prevq)) { msb=endj+1; follow_prevq=0; }
            else
                { dequeue(prevq,&msb); follow_prevq=1; }
            if (debug) { printf("no change in pre cell,follow msb=%d,",msb); fflush(stdout); }
            next_j=msb;
        }
        else
        {
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
        if (next_j==endj+1) {
            if (debug) {printf("next_j=endj+1,done with this row\n"); fflush(stdout);}  
        }
        else {
            int costchange = update_a_cell (s1,s2,m,c,i,next_j,newk,go);
            if (debug) { printf("update this cell.(%d,%d),costchange=%d\n",i,next_j,costchange); fflush(stdout); }
            //add next_j to this queue if necessary
            if ( costchange && !(in_new_left_band(i,next_j,newk,oldk)) ) {
                /*  we only push j to this queue when the cost changes
                    * in current cell. also if current cell is in the new left
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
            if ((costchange)||(startj==0)) { last_updated_j = next_j; }
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
            //see if we are done with this row
            if ((next_j==endj)||(last_updated_j==endj)) {}
            else {//not done, call itself again
                update_a_row(s1,s2,m,c,i,startj,endj,last_updated_j,prevq,thisq,newk,oldk,go);
            }
        }
    }
    return;
};

//function for speeding up end

void ukktest (const seqt s1, const seqt s2,newkkmat_p m, const cmt c,int currentT, int p, int lenX, int lenY, int go, int * res_cost, DIRECTION_MATRIX * res_gapnum)
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
    if (debug) { printf("\n ukktest,newk=%d,bb=%d,oldsize=%d,oldk=%d,go=%d\n",newk,bb,old_size,oldk,go);
    fflush(stdout); }
    expand_mat (m,newk,oldk,affine);
    int i,j;
    int startj,endj;
    //init two queue for speeding up function
    q_t prevq = (q_t) malloc(sizeof(struct Queue));
    init_empty_queue(prevq);
    q_t thisq = (q_t) malloc(sizeof(struct Queue));
    init_empty_queue(thisq);
    for (i=0;i<lenX;i++)
    {
        startj = ( (i-newk)>0 )?(i-newk):0 ;
        endj = ( (i+lenY-lenX+newk)<(lenY-1) )?(i+lenY-lenX+newk):(lenY-1) ;
        if (affine) {
            for (j= MAX(i-newk,0);j<=MIN(i+(lenY-lenX+newk),lenY-1);j++)
            {
            if (in_non_change_zone(i,j,oldk,lenX,lenY) ) {}
            else
                update_a_cell (s1,s2,m,c,i,j,newk,go);
            }
        }
        else {
            if ((debug)&&(is_emptyqueue(thisq)==0)) printf("Warning,thisQ is not empty before update_a_row\n");
            update_a_row (s1,s2,m,c,i,startj,endj,startj-1,prevq,thisq,newk,oldk,go);
            transfer_queue (thisq,prevq);
        }
    }
    int whichdiag=0, idx_in_my_diag=0, at_leftborder=0, at_rightborder=0;
    get_idx(lenX-1,lenY-1,m,&whichdiag,&idx_in_my_diag,&at_leftborder,&at_rightborder);
    int cost; DIRECTION_MATRIX dir; DIRECTION_MATRIX gapnum;
    int affP=INT_MAX/2, affQ=INT_MAX/2; //init affP and affQ anyway
    get_ukkcost(whichdiag,idx_in_my_diag,m, &cost, &dir, &gapnum,&affP, &affQ, affine);
    *res_cost = cost;
    *res_gapnum = gapnum;
};

int increaseT (const seqt s1, const seqt s2,newkkmat_p m, const cmt c,int newT, int lenX, int lenY, int go)
{
    int debug = 0;
    int p = (newT - (lenY-lenX))/2;
    if (debug) { printf ("increaseT, newT=%d,p=%d,",newT,p); fflush(stdout); }
    int res_cost=-1; 
    DIRECTION_MATRIX res_gapnum=-1;
    ukktest(s1,s2,m,c,newT,p,lenX,lenY,go,&res_cost,&res_gapnum);
    int newp = ( newT*2 - (lenY-lenX) )/2;
    assert(res_cost>=0); assert(res_gapnum>=0);
    if ((res_cost<=newT)  //we accept the cost -- threshold by ukkonen
        || (res_cost<newT*2) //we will accept this cost in next round anyway 
        || (newp>res_gapnum)) //max possible gap number -- threshold by ward's ukkonen -- newkkonen
    {
        if (debug) 
        {printf ("cost=%d,end by K(%d>%d) or T(%d<%d)\n",res_cost,newp,res_gapnum,res_cost,newT*2);
        fflush(stdout); }
        return res_cost;
    }
    else //double the T, continue
    {
        if (debug) { printf ("cost=%d,continue \n",res_cost);
        fflush(stdout); }
        return increaseT (s1,s2,m, c,2*newT, lenX, lenY, go);
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
        int * thisaffP = m->pool_affP;
        int * thisaffQ = m->pool_affQ; 
        ukkdiag_p thisdiag = m->diagonal;
        thisdiag -> affParr = thisaffP;
        thisdiag -> affQarr = thisaffQ;
        int diaglen;
        int diag_arr_size = m->diag_size;
        if (debug) { printf("init space for affine, oldlen=%d,diag_arr_size=%d\n",oldlen,diag_arr_size); fflush(stdout);}
        for (i=0;i<diag_arr_size;i++) {
            //printf ("init set diag#%d(%d);",i,thisdiag->len); fflush(stdout);
            diaglen = thisdiag->len;
            assert(diaglen>0);
            thisdiag -> affParr = thisaffP;
            thisdiag -> affQarr = thisaffQ;
            if (i<diag_arr_size-1) //don't move out of range
            {
            thisaffP += diaglen;
            thisaffQ += diaglen;
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
    if (debug) { printf ("newkkonen init_mat,oldlen = %li, lenx=%d,leny=%d,baseband=%d,oldlenaff=%d\n",oldlen,lenX,lenY,baseband,m->total_len_affine); fflush(stdout); }
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
    init_affine_mat(m);
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
        printf ("baseband=%d <=> m->diag_size=%d\n", baseband, m->diag_size); fflush(stdout); }
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
        m->pool_cost = realloc (m->pool_cost,len*sizeof(int));
        m->pool_dir = realloc (m->pool_dir,len*sizeof(DIRECTION_MATRIX));
        m->pool_gapnum = realloc (m->pool_gapnum,len*sizeof(DIRECTION_MATRIX));
        if (affine) {
            if (debug) {printf("expand affine P and Q diagnonal\n"); fflush(stdout);}
            m->pool_affP = realloc (m->pool_affP,len*sizeof(int));
        m->pool_affQ = realloc (m->pool_affQ,len*sizeof(int));}
        else
        {
            m->pool_affP = NULL;
            m->pool_affQ = NULL;
        }
    }   
        int * thiscost = m->pool_cost;
        //def these two anyway if (affine) { 
        int * thisaffP = m->pool_affP;
        int * thisaffQ = m->pool_affQ; //}
        DIRECTION_MATRIX * thisdir = m->pool_dir;
        DIRECTION_MATRIX * thisgap = m->pool_gapnum;
        ukkdiag_p thisdiag = m->diagonal;
        for (i=0;i<baseband;i++) {
            //printf ("init set diag#%d(%d);",i,thisdiag->len); fflush(stdout);
            thisdiag -> costarr = thiscost;
            thisdiag -> dirarr = thisdir;
            thisdiag -> gapnumarr = thisgap;
            if (affine) {
            thisdiag -> affParr = thisaffP;
            thisdiag -> affQarr = thisaffQ;}
            set_ukkcost(i,0,m,0,0,0,0,0,affine);
            if (i<baseband-1) //don't move out of range
            {
            thisdiag ++;
            thiscost += lenX ;
            thisdir += lenX;
            thisgap += lenX; 
            if (affine) {
                thisaffP += lenX;
                thisaffQ += lenX;}
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
        get_cmcost(c,thiscode,gapcode,&thiscost);
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
newkk_algn (const seqt s1, const seqt s2, MAT_SIZE s1_len, MAT_SIZE s2_len, int go, const cmt c, newkkmat_p m) {
    assert(s1_len<=s2_len);
    int debug = 0;
    int debug2 = 0;
    //set affine to 1 if gap openning is set
    int affine=0;
    if (go>=0) affine=1; 
    int realgo=0;
    if (go>0) realgo=go;
    if (debug)  { printf("newkk_algn,lenx=%d,leny=%d,affine=%d\n",s1_len,s2_len,affine);fflush(stdout);}
    if (s1_len * MUCH_LONGER < s2_len) 
    {
        return trivial_algn (s1,s2,s1_len,s2_len,c);
    }
    else 
    {
        //reset newkkmat to 0 if we are reusing the memory block.
        if (m->total_len > 0) reset_mat_to_0 (m);
        init_mat (s1_len,s2_len,m,affine);
        int gapcode = cm_get_gap (c); 
        int delta = get_delta (c);
        int i; 
        int bb = m->baseband;
        assert(bb>0);
        if (debug)
        {
            printf("baseband=%d,gapcode=%d,delta=%d,go=%d\n",bb,gapcode,delta,go);
            fflush(stdout);
        }
        set_ukkcost(0,0,m,0,START,0,0,0,affine);
        if (debug) {printf("init first row\n"); fflush(stdout);}
        for (i=1;i<bb;i++)
        {
           //if ((s1_len==234)&&(s2_len==10578)&&(i==3648)) debug2=1;
           //these 4 Int will be used again and again, reset them before use.
           if (debug2) { printf("init (0,%d) ",i); fflush(stdout); }
            int whichdiag=0, idx_in_my_diag=0, at_leftborder=0, at_rightborder=0;
            get_idx(0,i-1,m,&whichdiag,&idx_in_my_diag,&at_leftborder,&at_rightborder);
            int precost; DIRECTION_MATRIX predir; DIRECTION_MATRIX pre_gapnum; 
            int affP=INT_MAX/2, affQ=INT_MAX/2; //we init affP and affQ even we are not doing affine
            get_ukkcost(whichdiag,idx_in_my_diag,m, &precost, &predir, &pre_gapnum, &affP, &affQ,affine);
            whichdiag=0; idx_in_my_diag=0; at_leftborder=0; at_rightborder=0;
            get_idx(0,i,m,&whichdiag,&idx_in_my_diag,&at_leftborder,&at_rightborder);
            int thiscost=0;
            int thisP=INT_MAX/2,thisQ=0;
            DIRECTION_MATRIX dir=DO_INSERT;
            int thiscode = seq_get (s2,i);
            get_cmcost(c,thiscode,gapcode,&thiscost);
            thiscost += precost;
            if (i==1) { thiscost += realgo; }
            thisQ = thiscost;
            if (debug2) { printf("with %d; whichdiag=%d,idx=%d,dir=%d,gapnum=%d\n ", thiscost,whichdiag,idx_in_my_diag,dir,pre_gapnum+1); fflush(stdout);}
            set_ukkcost(whichdiag,idx_in_my_diag,m, thiscost, dir, pre_gapnum+1,thisP,thisQ,affine);
        }
        int iniT = (m->baseband) * delta;
        if (debug) {printf ("done with first row, call increaseT with iniT=%d\n",iniT); fflush(stdout);}
        int rescost = increaseT (s1,s2,m,c,iniT,s1_len,s2_len,go); 
        if (debug) {printf("end of newkk_algn,cost=%d\n\n",rescost); fflush(stdout);}
        FILE * outf;
        outf = fopen ("poy.out","w");
        fprintf(outf,"%d",rescost);
        fclose(outf);
        return rescost;
    }
};



value 
newkkonen_CAML_algn (value s1, value s2, value c, value a)
{
    CAMLparam4 (s1,s2,c,a);
    cmt tc;
    newkkmat_p mp;
    tc = Cost_matrix_struct(c);
    mp = Newkkmat_struct(a);
    seqt s1p, s2p;
    Seq_custom_val(s1p,s1);
    Seq_custom_val(s2p,s2);
    MAT_SIZE s1_len, s2_len;
    s1_len = seq_get_len (s1p);
    s2_len = seq_get_len (s2p);
    assert (s2_len >= s1_len);
    int res;
    //int gap_open = tc->gap_open; no need to pass gap opening if we are not doing affine
    int debug = 0;
    if (debug) {
    printf("newkkonen_CAML_algn,no gap opening\n"); fflush(stdout);}
    res = newkk_algn (s1p, s2p, s1_len,s2_len,(-1),tc, mp);
    CAMLreturn(Val_int(res));
};

value 
newkkonen_CAML_algn_affine (value s1, value s2, value c, value a)
{
    CAMLparam4 (s1,s2,c,a);
    cmt tc;
    newkkmat_p mp;
    tc = Cost_matrix_struct(c);
    mp = Newkkmat_struct(a);
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
    res = newkk_algn (s1p, s2p, s1_len,s2_len,gap_open,tc, mp);
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

void backtrace (const seqt s1, const seqt s2, seqt alis1, seqt alis2, 
        int len1, int len2, const cmt c, const newkkmat_p m,int affine)
{
   int debug = 0;
   if(debug) {
       printf("backtrace,len1=%d,len2=%d\n",len1,len2); fflush(stdout); }
   assert(len1<=len2);
   if (len2 > MUCH_LONGER * len1) 
   {
       trivial_backtrace(s1,s2,alis1,alis2,len1,len2,c);
   }
   else {
   int add1 = 0;
   int add2 = 0;
   int whichdiag, idx_in_my_diag, at_leftborder, at_rightborder;
   int cost; DIRECTION_MATRIX dir; DIRECTION_MATRIX gapnum; 
    //set affP and affQ anyway
   int affP=INT_MAX/2, affQ=INT_MAX/2;
   int i=len1-1, j=len2-1;
   while(i>=0&&j>=0)
   {
    if(debug) {printf ("i=%d,j=%d :",i,j); fflush(stdout); }
    whichdiag = 0, idx_in_my_diag=0, at_leftborder=0, at_rightborder=0;
    cost = 0; dir = 0; gapnum = 0;
    get_idx(i,j,m,&whichdiag,&idx_in_my_diag,&at_leftborder,&at_rightborder);
    get_ukkcost(whichdiag,idx_in_my_diag,m, &cost, &dir, &gapnum,&affP, &affQ, affine);
    if (dir==START) 
    {
        if(debug) { printf ("Start\n"); fflush(stdout);}
        assert (i==0 && j==0); 
        add1 = add2 = cm_get_gap(c);
        my_prepend(alis1,add1);
        my_prepend(alis2,add2);
        i--; j--;
    }
    else if (dir==DO_ALIGN)
    {
        if(debug) { printf ("Align\n");fflush(stdout);}
        add1 = seq_get(s1,i);
        add2 = seq_get(s2,j);
        my_prepend(alis1,add1);
        my_prepend(alis2,add2);
        i--; j--;
    }
    else if (dir==DO_DELETE)
    {
        if(debug) { printf ("Delete\n"); fflush(stdout);}
        add1 = seq_get(s1,i);
        add2 = cm_get_gap(c);
        my_prepend(alis1,add1);
        my_prepend(alis2,add2);
        i--;
    }
    else 
    {
        if(debug) { printf ("Insert\n"); fflush(stdout);}
        add2 = seq_get(s2,j);
        add1 = cm_get_gap(c);
        my_prepend(alis1,add1);
        my_prepend(alis2,add2);
        j--;
    }
   }
   }//non-trivial case
   if(debug) { printf ("end of backtrace\n"); fflush(stdout);}
};

value
newkkonen_CAML_backtrace (value s1, value s2, value s1p, value s2p, value c, value a) {
    CAMLparam5(s1, s2, s1p, s2p, c);
    CAMLxparam1(a);
    seqt ss1, ss2, ss1p, ss2p;
    newkkmat_p mp;
    cmt cc;
    mp = Newkkmat_struct(a);
    Seq_custom_val(ss1,s1);
    Seq_custom_val(ss2,s2);
    Seq_custom_val(ss1p,s1p);
    Seq_custom_val(ss2p,s2p);
    cc = Cost_matrix_struct(c);
    int affine=0;
    backtrace (ss1, ss2, ss1p, ss2p, seq_get_len(ss1), seq_get_len(ss2), cc, mp, affine);
    CAMLreturn(Val_unit);
};

value
newkkonen_CAML_backtrace_affine (value s1, value s2, value s1p, value s2p, value c, value a) {
    CAMLparam5(s1, s2, s1p, s2p, c);
    CAMLxparam1(a);
    seqt ss1, ss2, ss1p, ss2p;
    newkkmat_p mp;
    cmt cc;
    mp = Newkkmat_struct(a);
    Seq_custom_val(ss1,s1);
    Seq_custom_val(ss2,s2);
    Seq_custom_val(ss1p,s1p);
    Seq_custom_val(ss2p,s2p);
    cc = Cost_matrix_struct(c);
    int affine=1;
    backtrace (ss1, ss2, ss1p, ss2p, seq_get_len(ss1), seq_get_len(ss2), cc, mp, affine);
    CAMLreturn(Val_unit);
};


value 
newkkonen_CAML_backtrace_bc (value *argv, int argn) {
    return (newkkonen_CAML_backtrace (argv[0], argv[1], argv[2], argv[3], argv[4], argv[5]) );
};

value 
newkkonen_CAML_backtrace_affine_bc (value *argv, int argn) {
    return (newkkonen_CAML_backtrace (argv[0], argv[1], argv[2], argv[3], argv[4], argv[5]) );
};


void
newkkmat_CAML_free (value m) {
    int debug=1;
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
        t2->gapnumarr = NULL;
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



