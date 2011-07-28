
#include <stdio.h>
#include <malloc.h>
#include <stdio.h>
#include <limits.h>
#define NDEBUG 1
#include <assert.h>
#include <caml/mlvalues.h>
#include <caml/memory.h>
#include <caml/custom.h>
#include <caml/fail.h>

#include "seq.h"
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

void
expand_mat (newkkmat_p m,MAT_SIZE newk,MAT_SIZE oldk)
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
    { printf ("newkkonen.c, expand_mat,newk=%d(oldk=%d),newsize of diagarr=%d(old=%d),baseband=%d\n",
            newk,oldk,newsize,oldsize,baseband); fflush(stdout); }
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
        thisdiag = m->diagonal;
        thisdiag += i;
        thisdiag->len = diaglen;
    }
    //set expand sign
    if (debug) { 
        printf ("total_len_in_use=%li <=> total_len=%li\n", m->total_len_in_use, m->total_len); 
        fflush(stdout); }
    if (m->total_len_in_use <= m->total_len) expand_diagonal = 0;
    else 
    {
        if (debug) { printf ("we NEED to expand diagonal\n"); fflush(stdout);}
        expand_diagonal = 1;
        MAT_SIZE len_we_need = m->total_len_in_use;
        m->pool_cost = realloc (m->pool_cost,len_we_need*sizeof(int));
        m->pool_dir = realloc (m->pool_dir,len_we_need*sizeof(DIRECTION_MATRIX));
        m->pool_gapnum = realloc (m->pool_gapnum,len_we_need*sizeof(DIRECTION_MATRIX));
    }ukkdiag_p emptydiag;
        int * thiscost= m->pool_cost;
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
            if (i<newsize-1) //don't move out of range
            { 
            emptydiag ++;
            thiscost += diaglen ;
            thisdir += diaglen;
            thisgap += diaglen; }
        }
    //}
    //update diag_size and total_len if we expand our memory
    if (expand_diag_size) { m->diag_size = m->diag_size_in_use; }
    if (expand_diagonal) { m->total_len = m->total_len_in_use; }
    if (debug) 
    { printf ("Total size of memory we are using :%lu(%lu) cells, diag_size=%d(%d)\n",
            m->total_len_in_use,m->total_len,m->diag_size_in_use,m->diag_size); 
        fflush(stdout);}
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
    int debug = 0 ;
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
        if (debug) { printf("k_ij=%d,whichdiag=%d\n",k_ij,*whichdiag); fflush(stdout);}
    }
    else if ( (i<=j)&&( (j-i)<bb ) )
    { 
        *whichdiag = j-i; 
        if (debug) {printf("whichdiag=%d\n",*whichdiag); fflush(stdout);}
    }
    else {
        *whichdiag = bb + 2*(j-i+1 - bb) -1 -1;
        if (debug) {printf("k_ji=%d,whichdiag=%d\n",k_ji,*whichdiag); fflush(stdout); }
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
};


void 
get_ukkcost (int whichdiag, int idx_in_my_diag, newkkmat_p m, int * cost, DIRECTION_MATRIX * dir, DIRECTION_MATRIX *  max_gapnum)
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
    if (debug) { printf("get ukkcost,diag#.%d,idx#.%d,cost=%d,dir=%d,mapgapnum=%d\n",whichdiag,idx_in_my_diag,*cost,*dir,*max_gapnum); fflush(stdout);}
};

void set_ukkcost (int whichdiag, int idx_in_my_diag, newkkmat_p m, int cost, DIRECTION_MATRIX dir,DIRECTION_MATRIX  max_gapnum)
{
    int debug = 0;
    ukkdiag_p thisdiag = m->diagonal;
    if ((idx_in_my_diag >= thisdiag->len) || (idx_in_my_diag<0)) { debug=1; }
    if (debug) { printf ("set ukkcost, diag#.%d,idx#.%d,cost=%d,dir=%d,mapgapnum=%d\n",whichdiag,idx_in_my_diag,cost,dir,max_gapnum); fflush(stdout);}
    thisdiag += whichdiag;
    int * thiscostarr = thisdiag->costarr;
    DIRECTION_MATRIX * thisdirarr = thisdiag->dirarr;
    DIRECTION_MATRIX * thisgapnumarr = thisdiag->gapnumarr;
    assert(idx_in_my_diag < thisdiag->len);
    thiscostarr[idx_in_my_diag] = cost;
    thisdirarr[idx_in_my_diag] = dir;
    thisgapnumarr[idx_in_my_diag] = max_gapnum;

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
    sanity_check2(a);
    sanity_check2(b);
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

void update_internal_cell (const seqt s1, const seqt s2,newkkmat_p m, const cmt c, int i, int j, int newk)
{
    int debug = 0;
    if (debug) printf (" update internal cell ");
    int gapcode = cm_get_gap (c);
    int addcost = 0;
    int costL=0, costR=0, costM=0;
    int costfromL=0; DIRECTION_MATRIX dirL; DIRECTION_MATRIX gapnum_fromL=0;
    if (j==0) { costfromL = costL = INT_MAX/2; }
    else {
    int whichdiagL=0, idx_in_my_diagL=0, at_leftborderL=0, at_rightborderL=0;
    get_idx(i,j-1,m,&whichdiagL,&idx_in_my_diagL,&at_leftborderL,&at_rightborderL);
    get_ukkcost(whichdiagL,idx_in_my_diagL, m, &costfromL, &dirL, &gapnum_fromL);
    get_cmcost (c,seq_get(s2,j),gapcode,&addcost);
    costL = costfromL + addcost;
    }
    int costfromR=0; DIRECTION_MATRIX dirR; DIRECTION_MATRIX gapnum_fromR=0;
    if (i==0) { costfromR = costR = INT_MAX/2; }
    else {
    int whichdiagR=0, idx_in_my_diagR=0, at_leftborderR=0, at_rightborderR=0;
    get_idx(i-1,j,m,&whichdiagR,&idx_in_my_diagR,&at_leftborderR,&at_rightborderR);
    get_ukkcost(whichdiagR,idx_in_my_diagR, m, &costfromR, &dirR, &gapnum_fromR); 
    get_cmcost (c,seq_get(s1,i),gapcode,&addcost);
    costR = costfromR + addcost;
    }
    int costfromM=0; DIRECTION_MATRIX dirM; DIRECTION_MATRIX gapnum_fromM=0;
    if ((i==0)||(j==0)) { costfromM = costM = INT_MAX/2; }
    else {
    int whichdiagM=0, idx_in_my_diagM=0, at_leftborderM=0, at_rightborderM=0;
    get_idx(i-1,j-1,m,&whichdiagM,&idx_in_my_diagM,&at_leftborderM,&at_rightborderM);
    get_ukkcost(whichdiagM,idx_in_my_diagM, m, &costfromM, &dirM, &gapnum_fromM); 
    get_cmcost (c,seq_get(s1,i),seq_get(s2,j),&addcost);
    costM = costfromM + addcost;
    }
    int whichdiag=0, idx_in_my_diag=0, at_leftborder=0, at_rightborder=0;
    get_idx(i,j,m,&whichdiag,&idx_in_my_diag,&at_leftborder,&at_rightborder);
    if ((costL<=costR) && (costL<=costM)) 
    {
        set_ukkcost(whichdiag,idx_in_my_diag, m, costL, DO_INSERT, gapnum_fromL+1);
    }
    else if ((costR<=costL) && (costR<=costM)) 
    {
        set_ukkcost(whichdiag,idx_in_my_diag, m, costR, DO_DELETE, gapnum_fromR+1);
    }
    else
    {
        int newgapnum;
        if (costM==costfromM) newgapnum=gapnum_fromM;
        else newgapnum = gapnum_fromM + 1;
        set_ukkcost(whichdiag,idx_in_my_diag, m, costM, DO_ALIGN, newgapnum);
    }

};

void update_left_border_cell (const seqt s1, const seqt s2,newkkmat_p m, const cmt c, int i, int j, int newk)
{
    int debug = 0;
    if (debug) printf (" update left cell ");
    int gapcode = cm_get_gap (c);
    int addcost = 0;
    int costR=0, costM=0;
    int costfromR=0; DIRECTION_MATRIX dirR; DIRECTION_MATRIX gapnum_fromR=0;
    if (i==0) { costfromR = costR = INT_MAX/2; }
    else {
    int whichdiagR=0, idx_in_my_diagR=0, at_leftborderR=0, at_rightborderR=0;
    get_idx(i-1,j,m,&whichdiagR,&idx_in_my_diagR,&at_leftborderR,&at_rightborderR);
    get_ukkcost(whichdiagR,idx_in_my_diagR, m, &costfromR, &dirR, &gapnum_fromR); 
    get_cmcost (c,seq_get(s1,i),gapcode,&addcost);
    costR = costfromR + addcost;
    }
    int costfromM=0; DIRECTION_MATRIX dirM; DIRECTION_MATRIX gapnum_fromM=0;
    if ((i==0)||(j==0)) { costfromM = costM = INT_MAX/2; }
    else {
    int whichdiagM=0, idx_in_my_diagM=0, at_leftborderM=0, at_rightborderM=0;
    get_idx(i-1,j-1,m,&whichdiagM,&idx_in_my_diagM,&at_leftborderM,&at_rightborderM);
    get_ukkcost(whichdiagM,idx_in_my_diagM, m, &costfromM, &dirM, &gapnum_fromM); 
    get_cmcost (c,seq_get(s1,i),seq_get(s2,j),&addcost);
    costM = costfromM + addcost;
    }
    int whichdiag=0, idx_in_my_diag=0, at_leftborder=0, at_rightborder=0;
    get_idx(i,j,m,&whichdiag,&idx_in_my_diag,&at_leftborder,&at_rightborder);
    assert(at_leftborder);
    if (costR<=costM) 
    {
        set_ukkcost(whichdiag,idx_in_my_diag, m, costR, DO_DELETE, gapnum_fromR+1);
    }
    else
    {
        int newgapnum;
        if (costM==costfromM) newgapnum=gapnum_fromM;
        else newgapnum = gapnum_fromM + 1;
        set_ukkcost(whichdiag,idx_in_my_diag, m, costM, DO_ALIGN, newgapnum);
    } 
};

void update_right_border_cell (const seqt s1, const seqt s2,newkkmat_p m, const cmt c, int i, int j, int newk)
{
    int debug = 0;
    if (debug) printf (" update right cell ");
    int gapcode = cm_get_gap (c);
    int addcost = 0;
    int costL=0, costM=0;
    int costfromL=0; DIRECTION_MATRIX dirL; DIRECTION_MATRIX gapnum_fromL=0;
    if (j==0) { costfromL = costL = INT_MAX/2; }
    else {
    int whichdiagL=0, idx_in_my_diagL=0, at_leftborderL=0, at_rightborderL=0;
    get_idx(i,j-1,m,&whichdiagL,&idx_in_my_diagL,&at_leftborderL,&at_rightborderL);
    get_ukkcost(whichdiagL,idx_in_my_diagL, m, &costfromL, &dirL, &gapnum_fromL);
    get_cmcost (c,seq_get(s2,j),gapcode,&addcost);
    costL = costfromL + addcost;
    }
    int costfromM=0; DIRECTION_MATRIX dirM; DIRECTION_MATRIX gapnum_fromM=0;
    if ((i==0)||(j==0)) { costfromM = costM = INT_MAX/2; }
    else {
    int whichdiagM=0, idx_in_my_diagM=0, at_leftborderM=0, at_rightborderM=0;
    get_idx(i-1,j-1,m,&whichdiagM,&idx_in_my_diagM,&at_leftborderM,&at_rightborderM);
    get_ukkcost(whichdiagM,idx_in_my_diagM, m, &costfromM, &dirM, &gapnum_fromM); 
    get_cmcost (c,seq_get(s1,i),seq_get(s2,j),&addcost);
    costM = costfromM + addcost;
    }
    int whichdiag=0, idx_in_my_diag=0, at_leftborder=0, at_rightborder=0;
    get_idx(i,j,m,&whichdiag,&idx_in_my_diag,&at_leftborder,&at_rightborder);
    assert(at_rightborder);
    if (costL<=costM)
    {
        if (debug) {printf("costL=%d <= costM=%d",costL,costM); fflush(stdout);}
        set_ukkcost(whichdiag,idx_in_my_diag, m, costL, DO_INSERT, gapnum_fromL+1);
    }
    else
    {
        if (debug) {printf("costL=%d > costM=%d",costL,costM); fflush(stdout);}
        int newgapnum;
        if (costM==costfromM) newgapnum=gapnum_fromM;
        else newgapnum = gapnum_fromM + 1;
        set_ukkcost(whichdiag,idx_in_my_diag, m, costM, DO_ALIGN, newgapnum);
    }

};

void update_central_diagonal_cell (const seqt s1, const seqt s2,newkkmat_p m, const cmt c, int i, int j, int newk)
{
    int costM=0;
    int addcost = 0;
    int costfromM=0; DIRECTION_MATRIX dirM; DIRECTION_MATRIX gapnum_fromM=0;
    if ((i==0)||(j==0)) { costfromM = costM = INT_MAX/2; }
    else {
    int whichdiagM=0, idx_in_my_diagM=0, at_leftborderM=0, at_rightborderM=0;
    get_idx(i-1,j-1,m,&whichdiagM,&idx_in_my_diagM,&at_leftborderM,&at_rightborderM);
    get_ukkcost(whichdiagM,idx_in_my_diagM, m, &costfromM, &dirM, &gapnum_fromM); 
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
    set_ukkcost(whichdiag,idx_in_my_diag, m, costM, DO_ALIGN, newgapnum);
};

void update_a_cell (const seqt s1, const seqt s2,newkkmat_p m, const cmt c, int i, int j, int newk)
{
    int debug = 0;
    int whichdiag=0, idx_in_my_diag=0, at_leftborder=0, at_rightborder=0;
    get_idx(i,j,m,&whichdiag,&idx_in_my_diag,&at_leftborder,&at_rightborder);
    if ((whichdiag<0)||(whichdiag>= m->diag_size_in_use)) {debug=1;}
    if (debug) { printf("update_a_cell,(%d,%d),diag#.%d,idx#.%d,leftB=%d,rightB=%d\n",i,j,whichdiag,idx_in_my_diag,at_leftborder,at_rightborder); fflush(stdout);}
    if (at_leftborder&&at_rightborder) 
        update_central_diagonal_cell (s1, s2, m, c,i, j, newk);
    else if (at_leftborder)
        update_left_border_cell (s1, s2, m, c,i, j, newk);
    else if (at_rightborder)
        update_right_border_cell (s1, s2, m, c,i, j, newk);
    else
        update_internal_cell (s1, s2, m, c,i, j, newk);
};

void ukktest (const seqt s1, const seqt s2,newkkmat_p m, const cmt c,int currentT, int p, int lenX, int lenY, int * res_cost, DIRECTION_MATRIX * res_gapnum)
{
    int debug = 0;
    int newk=p;
    if (p>=lenX) newk=lenX-1;
    MAT_SIZE bb = m->baseband; 
    MAT_SIZE old_size = m->diag_size_in_use;
    MAT_SIZE oldk;
    oldk = (old_size-bb)/2;
    if (debug) { printf("\n ukktest,newk=%d,bb=%d,oldsize=%d,oldk=%d\n",newk,bb,old_size,oldk);
    fflush(stdout); }
    expand_mat (m,newk,oldk);
    int i,j;
    for (i=0;i<lenX;i++)
        for (j= MAX(i-newk,0);j<=MIN(i+(lenY-lenX+newk),lenY-1);j++)
        //for(j=0;j<lenY;j++)
        {
            if (/* outside_diagonal_area(i,j,newk,lenX,lenY) ||*/in_non_change_zone(i,j,oldk,lenX,lenY) ) {}
            else
                update_a_cell (s1,s2,m,c,i,j,newk);
        }
    int whichdiag=0, idx_in_my_diag=0, at_leftborder=0, at_rightborder=0;
    get_idx(lenX-1,lenY-1,m,&whichdiag,&idx_in_my_diag,&at_leftborder,&at_rightborder);
    int cost; DIRECTION_MATRIX dir; DIRECTION_MATRIX gapnum;
    get_ukkcost(whichdiag,idx_in_my_diag,m, &cost, &dir, &gapnum);
    *res_cost = cost;
    *res_gapnum = gapnum;
};

int increaseT (const seqt s1, const seqt s2,newkkmat_p m, const cmt c,int newT, int lenX, int lenY)
{
    int debug = 1;
    int p = (newT - (lenY-lenX))/2;
    if (debug) { printf ("increaseT, newT=%d,p=%d,",newT,p); fflush(stdout); }
    int res_cost=-1; 
    DIRECTION_MATRIX res_gapnum=-1;
    ukktest(s1,s2,m,c,newT,p,lenX,lenY,&res_cost,&res_gapnum);
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
        increaseT (s1,s2,m, c,2*newT, lenX, lenY);
    };
    //this does not work, why?
    //assert(0==1);
    //return -1;
};

void
init_mat (MAT_SIZE lenX, MAT_SIZE lenY,newkkmat_p m)
{
    int debug = 1;
    int i=0;
    //expand sign
    int expand_diag_size=1;
    int expand_diagonal=1;
    MAT_SIZE len=0, oldlen=0;
    MAT_SIZE baseband=0;
    baseband = lenY-lenX+1;
    oldlen = m->total_len;
    len = lenX * baseband;
    if (debug) { printf ("newkkonen init_mat,oldlen = %li, lenx=%d,leny=%d,baseband=%d,\n",oldlen,lenX,lenY,baseband); fflush(stdout); }
    assert(m != NULL);
    //k is init to 0
    m->k=0;
    //set baseband , this won't change during this alignment
    m->baseband = baseband;
    //init diag_size_in_use with baseband,may increase later
    m->diag_size_in_use = baseband; //array of diagonal size, start with [lenx,lenx,....]
    //set total_len_in_use to len, may increase later
    m->total_len_in_use = len;
    if (debug) { printf ("sizeofint=%lu,sizeoflong=%lu,len=%li <=> oldlen=%li\n", sizeof(int),sizeof(long),len, oldlen); fflush(stdout);}
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
    }   
    int * thiscost = m->pool_cost;
        DIRECTION_MATRIX * thisdir = m->pool_dir;
        DIRECTION_MATRIX * thisgap = m->pool_gapnum;
        ukkdiag_p thisdiag = m->diagonal;
        for (i=0;i<baseband;i++) {
            //printf ("init set diag#%d(%d);",i,thisdiag->len); fflush(stdout);
            thisdiag -> costarr = thiscost;
            thisdiag -> dirarr = thisdir;
            thisdiag -> gapnumarr = thisgap;
            set_ukkcost(i,0,m,0,0,0);
            if (i<baseband-1) //don't move out of range
            {
            thisdiag ++;
            thiscost += lenX ;
            thisdir += lenX;
            thisgap += lenX; }
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
    int debug = 1;
    if (debug)  { printf("trivial algn\n");fflush(stdout);}
    int indelcost1 = get_indel_cost (s1,s1_len,c);
    int indelcost2 = get_indel_cost (s2,s2_len,c);
    if (debug) { printf("end of trivial algn, cost=%d\n",indelcost1+indelcost2); fflush(stdout);}
    return (indelcost1 + indelcost2);
}



int
newkk_algn (const seqt s1, const seqt s2, MAT_SIZE s1_len, MAT_SIZE s2_len, const cmt c, newkkmat_p m) {
    assert(s1_len<=s2_len);
    int debug = 1;
    int debug2 = 0;
    if (debug)  { printf("newkk_algn,lenx=%d,leny=%d\n",s1_len,s2_len);fflush(stdout);}
    if (s1_len * MUCH_LONGER < s2_len) 
    {
        trivial_algn (s1,s2,s1_len,s2_len,c);
    }
    else {
    init_mat (s1_len,s2_len,m);
    int gapcode = cm_get_gap (c); 
    int delta = get_delta (c);
    int i; 
    int bb = m->baseband;
    assert(bb>0);
    if (debug)
    {
        printf("baseband=%d,gapcode=%d,delta=%d\n",bb,gapcode,delta);
    fflush(stdout);
    }
    set_ukkcost(0,0,m,0,START,0);
    if (debug) {printf("init first row\n"); fflush(stdout);}
    for (i=1;i<bb;i++)
    {
        
       //if ((s1_len==234)&&(s2_len==10578)&&(i==3648)) debug2=1;
       //these 4 Int will be used again and again, reset them before use.
       if (debug2) { printf("init (0,%d) ",i); fflush(stdout); }
        int whichdiag=0, idx_in_my_diag=0, at_leftborder=0, at_rightborder=0;
        get_idx(0,i-1,m,&whichdiag,&idx_in_my_diag,&at_leftborder,&at_rightborder);
        int precost; DIRECTION_MATRIX predir; DIRECTION_MATRIX pre_gapnum;
        get_ukkcost(whichdiag,idx_in_my_diag,m, &precost, &predir, &pre_gapnum);
        whichdiag=0; idx_in_my_diag=0; at_leftborder=0; at_rightborder=0;
        get_idx(0,i,m,&whichdiag,&idx_in_my_diag,&at_leftborder,&at_rightborder);
        int thiscost=0; 
        DIRECTION_MATRIX dir=DO_INSERT;
        int thiscode = seq_get (s2,i);
        get_cmcost(c,thiscode,gapcode,&thiscost);
        thiscost += precost;
        if (debug2) { printf("with %d; whichdiag=%d,idx=%d,dir=%d,gapnum=%d\n ", thiscost,whichdiag,idx_in_my_diag,dir,pre_gapnum+1); fflush(stdout);}
        set_ukkcost(whichdiag,idx_in_my_diag,m, thiscost, dir, pre_gapnum+1);
    }
    int iniT = (m->baseband) * delta;
    if (debug) {printf ("done with first row, call increaseT with iniT=%d\n",iniT); fflush(stdout);}
    int rescost = increaseT (s1,s2,m,c,iniT,s1_len,s2_len); 
    if (debug) {printf("end of newkk_algn,cost=%d\n\n",rescost); fflush(stdout);}
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
    res = newkk_algn (s1p, s2p, s1_len,s2_len,tc, mp);
    CAMLreturn(Val_int(res));
};

void trivial_backtrace(const seqt s1, const seqt s2, seqt alis1, seqt alis2, 
        int len1, int len2, const cmt c) 
{
    int debug = 1;
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
        int len1, int len2, const cmt c, const newkkmat_p m)
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
   int i=len1-1, j=len2-1;
   while(i>=0&&j>=0)
   {
    if(debug) {printf ("i=%d,j=%d :",i,j); fflush(stdout); }
    whichdiag = 0, idx_in_my_diag=0, at_leftborder=0, at_rightborder=0;
    cost = 0; dir = 0; gapnum = 0;
    get_idx(i,j,m,&whichdiag,&idx_in_my_diag,&at_leftborder,&at_rightborder);
    get_ukkcost(whichdiag,idx_in_my_diag,m, &cost, &dir, &gapnum);
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
   printf ("end of backtrace\n"); fflush(stdout);
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
    backtrace (ss1, ss2, ss1p, ss2p, seq_get_len(ss1), seq_get_len(ss2), cc, mp);
    CAMLreturn(Val_unit);
};

value 
newkkonen_CAML_backtrace_bc (value *argv, int argn) {
    return (newkkonen_CAML_backtrace (argv[0], argv[1], argv[2], argv[3], argv[4], argv[5]) );
};


void
newkkmat_CAML_free (value m) {
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
    //free cost_dir array
    free (tmp->diagonal);
    free (tmp->diagonal_size_arr);
    //free memory pool
    free(tmp->pool_cost);
    free(tmp->pool_dir);
    free(tmp->pool_gapnum);
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
    CAMLreturn(res);
}

value 
newkkonen_CAML_initialize (value unit) {
    CAMLparam1(unit);
    caml_register_custom_operations (&newkk_alignment_matrix);
    CAMLreturn(Val_unit);
}



