/* POY 5.1.1. A phylogenetic analysis program using Dynamic Homologies.       */
/* Copyright (C) 2014 Andrés Varón, Lin Hong, Nicholas Lucaroni, Ward Wheeler,*/
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

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <string.h>
#include <assert.h>
#include <caml/mlvalues.h>
#include <caml/memory.h>
#include <caml/custom.h>
#include <caml/fail.h>
#include <caml/bigarray.h>
#include <caml/intext.h>

#include "sankoff.h"


#define median_3_su 0

value sankoff_GC_custom_max( value n )
{
    CAMLparam1( n );
    alloc_custom_max = Int_val( n );
    CAMLreturn( Val_unit );
}


void free_elt (elt_p ep)
{
    if(ep->states!=NULL) free(ep->states);
    if(ep->leftstates!=NULL) free(ep->leftstates);
    if(ep->rightstates!=NULL) free(ep->rightstates);
    //for new median_3
    if (median_3_su) {
        if(ep->left_costdiff_mat!=NULL) free(ep->left_costdiff_mat);
        if(ep->right_costdiff_mat!=NULL) free(ep->right_costdiff_mat);
    }
    if(ep->beta!=NULL) free(ep->beta);
    if(ep->e!=NULL) free(ep->e);
    if(ep->m!=NULL) free(ep->m);
    //we don't free(ep) in this function
    return;
}


void free_eltarr(eltarr_p eap) {
    int i;
    if (eap->tcm!=NULL)
        free(eap->tcm);
    for (i=0;i<eap->num_elts;i++)
        free_elt(&((eap->elts)[i]));
    free(eap->elts);
    return;
}

void sankoff_CAML_free_eltarr (value v) {
    eltarr_p eap;
    eap = Sankoff_return_eltarr(v);
    assert(eap!=NULL);
    if (eap->tcm!=NULL) free(eap->tcm);
    int i;
    for (i=0;i<eap->num_elts;i++)
        free_elt(&((eap->elts)[i]));
    free(eap->elts);
    //if (debug) {printf("free eap\n"); fflush(stdout);}
    free(eap);
    return;
}

//similar to free_elt(elt_p ep), but we free ep itself here
void
sankoff_CAML_free_elt (value v) {
    elt_p ep;
    ep = Sankoff_return_elt(v);
    if(ep->states!=NULL) free(ep->states);
    if(ep->leftstates!=NULL) free(ep->leftstates);
    if(ep->rightstates!=NULL) free(ep->rightstates);
    // free for new median_3
    if (median_3_su) { 
        if(ep->left_costdiff_mat!=NULL) free(ep->left_costdiff_mat);
        if(ep->right_costdiff_mat!=NULL) free(ep->right_costdiff_mat);
    }
    if(ep->beta!=NULL) free(ep->beta);
    if(ep->e!=NULL) free(ep->e);
    if(ep->m!=NULL) free(ep->m);
    free(ep);
    return;
}


static struct custom_operations sankoff_custom_operations_elt = {
    "http://www.amnh.org/poy/",
    &sankoff_CAML_free_elt,
    custom_compare_default,
    custom_hash_default,
    custom_serialize_default, 
    custom_deserialize_default
};

static struct custom_operations sankoff_custom_operations_eltarr = {
    "http://www.amnh.org/poy/",
    &sankoff_CAML_free_eltarr,
    custom_compare_default,
    custom_hash_default,
    custom_serialize_default, 
    custom_deserialize_default
};

value sankoff_CAML_register_elt (value u) {
    CAMLparam1(u);
    register_custom_operations (&sankoff_custom_operations_elt);
    CAMLreturn (Val_unit);
}

value sankoff_CAML_register_eltarr (value u) {
    CAMLparam1(u);
    register_custom_operations (&sankoff_custom_operations_eltarr);
    CAMLreturn (Val_unit);
}


//return 0 if two int array are the same, 1 otherwise
#ifdef _WIN32
__inline int
#elif __clang__
int
#else
inline int
#endif 
sankoff_compare_two_int_array (int * arr1, int * arr2, int size) {
    int res=0;//init res to 0
    int i;
    for (i=0;i<size;i++) {
        if(arr1[i]!=arr2[i]) res=1;
        else {}
    }
    return res;
}

void 
sankoff_print_int_array (char * str, int * arr, int size)
{
    printf("%s,",str);
    int i, tmp;
    for (i=0;i<size;i++)
    {
        tmp = arr[i];
        if(is_infinity(tmp)) printf ("inf,");
        else
        printf("%d,",tmp);
    }
    printf("\n");
    fflush(stdout);
}

//return 1 if x is a member of array. 0 otherwise
#ifdef _WIN32
__inline int
#elif __clang__
int
#else
inline int
#endif
int_array_is_mem (int * arr, int size, int x) {
    int i; int res=0;
    for(i=0;i<size;i++) {
        if(x==arr[i]) res=1;
    }
    return res;
}


//given a matrix of sizex * sizey, return pointer to the start of line#.i,
//i start from 0
#ifdef _WIN32
__inline int *
#elif __clang__
int*
#else
inline int *
#endif
sankoff_move_to_line_i (int * arrarr, int sizex, int sizey, int i) {
    assert(i<sizex);//i start from 0.
    return ( arrarr + sizey * i);
}

//given a int mat of sizex*sizey, return value on mat.(i).(j).
//i,j start from 0
#ifdef _WIN32
__inline int
#elif __clang__
int
#else
inline int
#endif
sankoff_return_value (int * arrarr, int sizex, int sizey, int i, int j) {
    assert(i<sizex);
    int * start_pos = sankoff_move_to_line_i(arrarr, sizex,sizey,i);
    assert(j<sizey);//j start form 0
    return (*(start_pos + j));
}

//return 1 if epN is left child of epA, 0 if it's right child
#ifdef _WIN32
__inline int
#elif __clang__
int
#else
inline int
#endif
sankoff_is_left_or_right_child (eltarr_p eapN, eltarr_p eapA) {
    int debug = 0;
    if(debug) printf("sankoff_is_left_or_right_child, A's child left/right = %d/%d, N's taxoncode=%d\n",
            eapA->left_taxon_code,eapA->right_taxon_code,eapN->taxon_code);
    if (eapN->taxon_code == eapA->left_taxon_code)
        return 1;
    else if(eapN->taxon_code == eapA->right_taxon_code) 
        return 0;
    else 
        failwith ("sankoff, node N is neither left nor right child of node A");
}

//return 1 if this node is a leaf, for a leaf node , its left and right child
//taxon code should be both 0.
#ifdef _WIN32
__inline int
#elif __clang__
int
#else
inline int
#endif
sankoff_is_leaf_node (eltarr_p eap) {
    if (eap->left_taxon_code == eap->right_taxon_code)
    {
        assert(eap->left_taxon_code==0);
        return 1;
    }
    else 
        return 0;
}


void
sankoff_print_elt (elt_p ep, int printe, int printbeta, int print_costdiff, int print_best_child_states)
{
    printf("{ ");
    //printf("[ ecode=%d, number of states=%d ] \n",ep->ecode,ep->num_states);
    int num_states = ep->num_states;
    sankoff_print_int_array("states:",ep->states,num_states);
    if (printe) sankoff_print_int_array("e:",ep->e,num_states);
    fflush(stdout);
    if (printbeta) sankoff_print_int_array("beta:",ep->beta,num_states);
    fflush(stdout);
    if (ep->m_already_set) sankoff_print_int_array("m:",ep->m,num_states);
    fflush(stdout);
    if (print_best_child_states) {
        sankoff_print_int_array("left child states :",ep->leftstates,num_states);
        sankoff_print_int_array("right child states :",ep->rightstates, num_states);
    }
    fflush(stdout);
    // for new median_3
    if (median_3_su) { 
        if (print_costdiff) {
            int * tmp; int i;
            printf("cost diff for left child\n");
            for (i=0;i<num_states;i++) {
            tmp = sankoff_move_to_line_i(ep->left_costdiff_mat,num_states,num_states,i);
            sankoff_print_int_array("", tmp,num_states);
            }
            printf("cost diff for right child\n");
            for (i=0;i<num_states;i++) {
            tmp = sankoff_move_to_line_i(ep->right_costdiff_mat,num_states,num_states,i);
            sankoff_print_int_array("",
                   tmp,num_states);
            }
        } 
    }
    //if (printbs) sankoff_print_int_array("best_states:",ep->best_states,num_states);
    printf("}\n");
    fflush(stdout);
    return;
}

void 
sankoff_print_eltarr (eltarr_p eap, int printe, int printbeta, int print_costdiff, int print_best_child_states)
{
    printf("taxon code = %d, code=%d,num_elts=%d,sum cost = %li, left child = %d, right child = %d\n",
            eap->taxon_code,eap->code,eap->num_elts,eap->sum_cost,eap->left_taxon_code,eap->right_taxon_code);
    int i;
    int num_elts = eap->num_elts;
    int printtcm = 1;
    if(printtcm) {
        int num_states = eap->num_states;
        sankoff_print_int_array("tcm = ", eap->tcm, num_states*num_states);
    }
    for (i=0;i<num_elts;i++)
        sankoff_print_elt(&((eap->elts)[i]),printe,printbeta,print_costdiff,print_best_child_states);
    fflush(stdout);
    return;
}


//create new empty elt to pointer newelt.
//allocate memory for what its pointers pointing to -- states
//array, beta array, e array, m array and best_states array.
//init num_states and ecode. ecode could be (-1).
#ifdef _WIN32
__inline void
#elif __clang__
void
#else
inline void
#endif
sankoff_create_empty_elt (elt_p newelt, int num_states,int ecode) {
    newelt->ecode = ecode;
    newelt->num_states = num_states;
    newelt->states = (int*)calloc(num_states,sizeof(int));
    newelt->leftstates = (int*)calloc(num_states,sizeof(int));
    newelt->rightstates = (int*)calloc(num_states,sizeof(int));
    // for new median_3
    if (median_3_su) { 
        newelt->left_costdiff_mat = (int*)calloc(num_states*num_states,sizeof(int));
    newelt->right_costdiff_mat = (int*)calloc(num_states*num_states,sizeof(int));
    }
    newelt->beta = (int*)calloc(num_states,sizeof(int));
    newelt->e = (int*)calloc(num_states,sizeof(int));
    newelt->m = (int*)calloc(num_states,sizeof(int));
    return;
}

//init a new eltarr with necessary memory and code&tcm
//neweltarr itself must be already allocated, we only allocate what its pointers
//pointing to.
//allocate memory for cost matrix and elts.
//fill in  code,copy cost matrix.
#ifdef _WIN32
__inline void
#elif __clang__
void
#else
inline void
#endif
sankoff_init_eltarr (eltarr_p neweltarr, int num_states, int num_elts, int code, int taxon_code, int left_taxon_code, int right_taxon_code, int * tcm, int is_id) {
    int debug = 0;
    if(debug){
        printf("sankoff_init_eltarr, ltcode=%d, rtcode=%d, tcode=%d, code = %d\n",
                left_taxon_code, right_taxon_code,taxon_code,code);
    }
    neweltarr->code = code;
    neweltarr->taxon_code = taxon_code;
    neweltarr->left_taxon_code = left_taxon_code;
    neweltarr->right_taxon_code = right_taxon_code;
    neweltarr->sum_cost = 0;
    neweltarr->num_states=num_states;
    neweltarr->num_elts=num_elts;
    neweltarr->tcm = (int*)calloc(num_states*num_states,sizeof(int));
    neweltarr->is_identity = is_id;
    memcpy(neweltarr->tcm,tcm,sizeof(int)*num_states*num_states);
    neweltarr->elts = (elt_p)calloc(num_elts,sizeof(struct elt));
    return;
}

value
sankoff_CAML_get_code (value this_eltarr)
{
    CAMLparam1(this_eltarr);
    CAMLlocal1(res);
    eltarr_p eap;
    Sankoff_eltarr_custom_val(eap,this_eltarr);
    CAMLreturn(Val_int(eap->code));
}

value
sankoff_CAML_get_num_elts (value this_eltarr)
{
    CAMLparam1(this_eltarr);
    CAMLlocal1(res);
    eltarr_p eap;
    Sankoff_eltarr_custom_val(eap,this_eltarr);
    CAMLreturn(Val_int(eap->num_elts));
}

value
sankoff_CAML_get_taxon_code (value this_eltarr)
{
    CAMLparam1(this_eltarr);
    CAMLlocal1(res);
    eltarr_p eap;
    Sankoff_eltarr_custom_val(eap,this_eltarr);
    CAMLreturn(Val_int(eap->taxon_code));
}

value
sankoff_CAML_get_ecode(value this_elt) {
    CAMLparam1(this_elt);
    CAMLlocal1(res);
    elt_p ep;
    Sankoff_elt_custom_val(ep,this_elt);
    CAMLreturn(Val_int(ep->ecode));
}

value
sankoff_CAML_get_states(value this_elt, value this_or_left_or_right) {
    CAMLparam2(this_elt,this_or_left_or_right);
    CAMLlocal1(res);
    elt_p ep;
    Sankoff_elt_custom_val(ep,this_elt);
    int whatstate = Int_val(this_or_left_or_right);
    int num_states = ep->num_states;
    long dims[1];
    dims[0] = num_states;
    if (whatstate==1) 
        CAMLreturn( alloc_bigarray(BIGARRAY_INT32 | BIGARRAY_C_LAYOUT, 1, ep->states, dims));
    else if (whatstate==2) 
        CAMLreturn( alloc_bigarray(BIGARRAY_INT32 | BIGARRAY_C_LAYOUT, 1, ep->leftstates, dims));
    else if (whatstate==3) 
        CAMLreturn( alloc_bigarray(BIGARRAY_INT32 | BIGARRAY_C_LAYOUT, 1, ep->rightstates, dims));
    else
        failwith ("sankoff_CAML_get_states,must pick which state array you need (1=states,2=leftstates,3=rightstates)");
}




//copy ep1 to ep2, ep2 must be init with enough memory already for what its
//pointers pointing to
#ifdef _WIN32
__inline void
#elif __clang__
void
#else
inline void
#endif
sankoff_clone_elt (elt_p ep2, elt_p ep1) {
    ep2->ecode = ep1->ecode;
    ep2->num_states = ep1->num_states;
    ep2->m_already_set = ep1->m_already_set;
    int num_states = ep1->num_states;
    memcpy(ep2->states,ep1->states,sizeof(int)*num_states);
    memcpy(ep2->leftstates,ep1->leftstates,sizeof(int)*num_states);
    memcpy(ep2->rightstates,ep1->rightstates,sizeof(int)*num_states);
    // for new median_3
    if (median_3_su) { 
        memcpy(ep2->left_costdiff_mat,ep1->left_costdiff_mat,sizeof(int)*num_states*num_states);
    memcpy(ep2->right_costdiff_mat,ep1->right_costdiff_mat,sizeof(int)*num_states*num_states);
    }
    memcpy(ep2->e,ep1->e,sizeof(int)*num_states);
    memcpy(ep2->beta,ep1->beta,sizeof(int)*num_states);
    memcpy(ep2->m,ep1->m,sizeof(int)*num_states);
}

//return 0 if two eltarr are the same, 1 otherwise
//we compare array states, e, beta. not m and best_states
int
sankoff_compare_elt(elt_p ep1, elt_p ep2) {
    int num_states = ep1->num_states;
    if (ep1->num_states != ep2->num_states) 
        return 1;
    else if (ep1->ecode != ep2->ecode)
        return 1;
    else if (ep1->m_already_set != ep2->m_already_set)
        return 1;
    else {
        int res = 0;
        res = sankoff_compare_two_int_array(ep1->states,ep2->states,num_states);
        if (res==1) return 1;
        else {
            res = sankoff_compare_two_int_array(ep1->e,ep2->e,num_states);
            if (res==1) return 1;
            else {
                res = sankoff_compare_two_int_array(ep1->beta,ep2->beta,num_states);
                if (res==1) return 1;
                else {
                    int res1 = sankoff_compare_two_int_array(ep1->leftstates,ep2->leftstates,num_states);
                    int res2 = sankoff_compare_two_int_array(ep1->rightstates,ep2->rightstates,num_states);
                    res = res1 + res2;
                    if (res>=1) {
                        return 1;
                    } else {
                        if(ep1->m_already_set){//m is set 
                            res = sankoff_compare_two_int_array(ep1->m,ep2->m,num_states);
                            return (res==1)?1:0;
                        } else { //m is not set,all other arrays are the same, return 0
                            return 0;
                        }
                    }
                }//left&right states are the same
            }//beta are the same
        }//states are the same
    }//num_states,ecode,m_already_set are the same
}


//return 0 if two eltarr are the same, 1 otherwise
value
sankoff_CAML_compare_eltarr(value eltarr1, value eltarr2) {
    CAMLparam2(eltarr1,eltarr2);
    eltarr_p eap1; 
    Sankoff_eltarr_custom_val(eap1,eltarr1);
    eltarr_p eap2; 
    Sankoff_eltarr_custom_val(eap2,eltarr2);
    int res=0;//init res to 0
    if (eap1->code != eap2->code ) res=1;
    else if (eap1->num_states != eap2->num_states) res=1;
    else if (eap1->num_elts != eap2->num_elts) res=1;
    else if (eap1->taxon_code != eap2->taxon_code) res=1;
    else if (eap1->left_taxon_code != eap2->left_taxon_code) res=1;
    else if (eap1->right_taxon_code != eap2->right_taxon_code) res=1;
    else {//if any of the elt is different, set res to 1
        int i, tmp;
        for (i=0;i<eap1->num_elts;i++) {
            tmp = sankoff_compare_elt(&((eap1->elts)[i]),&((eap2->elts)[i]));
            if (tmp==1) res=1;
        }
    }
    CAMLreturn(Val_int(res));
}

//please note that we use alloc_custom_max here like we did for get_eltarr. I
//know that each eltarr might have more than 1 elt, so this number should be
//bigger. But
//we only get elt directly when dealing with fixed_state, or diagnosis output to
//screen. for fixed_state, each eltarr only has one elt. for diagnosis, we only
//print result best tree(s) with Preliminary&Final&Single states. 
//either way, we don't need that many. 
value
sankoff_CAML_get_elt (value this_eltarr,value idx)
{
    CAMLparam2(this_eltarr,idx);
    CAMLlocal1(res);
    eltarr_p eap; 
    Sankoff_eltarr_custom_val(eap,this_eltarr);
    int debug = 0;
    int i = Int_val(idx);
    int num_states = eap->num_states;
    elt_p ep;
    //alloc ep before create_emtpy_elt
    ep = (elt_p)calloc(1,sizeof(struct elt));
    sankoff_create_empty_elt (ep,num_states,-1);
    sankoff_clone_elt( ep, &((eap->elts)[i]) );
    if (debug) { 
        printf("sankoff_CAML_get_elt NO.%d from eltarr, res = ",i);
        sankoff_print_elt(ep,1,1,1,1);
    }
    res = caml_alloc_custom(&sankoff_custom_operations_elt,sizeof(elt_p),1,3*alloc_custom_max);
    Sankoff_return_elt(res) = ep;
    CAMLreturn(res);
}

value
sankoff_CAML_get_tcm (value this_eltarr)
{
    CAMLparam1(this_eltarr);
    CAMLlocal1(res);
    eltarr_p eap;
    Sankoff_eltarr_custom_val(eap,this_eltarr);
    int num_states=eap->num_states;
    long dims[2];
    dims[0] = num_states; dims[1] = num_states;
    CAMLreturn(alloc_bigarray(BIGARRAY_INT32 | BIGARRAY_C_LAYOUT, 2, eap->tcm, dims));
}

value
sankoff_CAML_get_sumcost (value this_eltarr)
{
    CAMLparam1(this_eltarr);
    CAMLlocal1(res);
    eltarr_p eap;
    Sankoff_eltarr_custom_val(eap,this_eltarr);
    long int sum_cost=eap->sum_cost;
    CAMLreturn(Val_long(sum_cost));
    
}

value 
sankoff_CAML_filter_character(value this_eltarr, value ecode_bigarr, value get_comp) {
    CAMLparam3(this_eltarr,ecode_bigarr,get_comp);
    CAMLlocal1(res);
    int get_complementary = Int_val(get_comp);
    int * ecode_arr = (int*) Data_bigarray_val(ecode_bigarr);
    eltarr_p eap;
    Sankoff_eltarr_custom_val(eap,this_eltarr);
    int num_elts=eap->num_elts; 
    int res_num_elts=0;//must init to 0
    int i;
    elt_p ep;
    int * sign_arr = (int*)calloc(num_elts,sizeof(int));
    for (i=0;i<num_elts;i++) {
        ep = &((eap->elts)[i]);
        if(get_complementary) {
            if( int_array_is_mem(ecode_arr,num_elts,ep->ecode) )
            {  sign_arr[i]=0;  }
            else { sign_arr[i]=1; res_num_elts++; }
        }
        else {
            if( int_array_is_mem(ecode_arr,num_elts,ep->ecode) )
            {  sign_arr[i]=1; res_num_elts++; }
            else sign_arr[i]=0;
        }
        
    }
    eltarr_p res_eap;
    res_eap = (eltarr_p)calloc(1,sizeof(struct elt_arr));
    int num_states = eap->num_states;
    res_eap->code = eap->code;
    res_eap->num_states = eap->num_states;
    res_eap->num_elts = res_num_elts;
    res_eap->tcm = (int*)calloc(num_states*num_states,sizeof(int));
    res_eap->is_identity = eap->is_identity;
    memcpy (res_eap->tcm,eap->tcm,sizeof(int)*num_states*num_states);
    res_eap->elts = (elt_p)calloc(res_num_elts,sizeof(struct elt));
    elt_p res_elts = res_eap->elts;
    int j=0;
    for (i=0;i<num_elts;i++) {
        if(sign_arr[i]==1) {
         sankoff_create_empty_elt(&(res_elts[j]),num_states,-1);
         sankoff_clone_elt(&(res_elts[j]),&((eap->elts)[i]));
         j++;
        }
    }
    free(sign_arr);
    assert(j==res_num_elts);
    res = caml_alloc_custom(&sankoff_custom_operations_eltarr,sizeof(eltarr_p),1,alloc_custom_max);
    Sankoff_return_eltarr(res) = res_eap;
    CAMLreturn(res);
}


#ifdef _WIN32
__inline void
#elif __clang__
void
#else
inline void
#endif
sankoff_set_state (elt_p ep, int pos, int c)
{
    assert(pos<ep->num_states);
    (ep->states)[pos] = c;
    return;
}

#ifdef _WIN32
__inline void
#elif __clang__
void
#else
inline void
#endif
sankoff_set_leftstate (elt_p ep, int pos, int c)
{
    assert(pos<ep->num_states);
    (ep->leftstates)[pos] = c;
    return;
}

#ifdef _WIN32
__inline void
#elif __clang__
void
#else
inline void
#endif
sankoff_set_rightstate (elt_p ep, int pos, int c)
{
    assert(pos<ep->num_states);
    (ep->rightstates)[pos] = c;
    return;
}



#ifdef _WIN32
__inline int
#elif __clang__
int
#else
inline int
#endif 
sankoff_get_state (elt_p ep, int pos)
{
    assert(pos<ep->num_states);
    return (ep->states)[pos];
}

#ifdef _WIN32
__inline int
#elif __clang__
int
#else
inline int
#endif 
sankoff_get_leftstate (elt_p ep, int pos)
{
    assert(pos<ep->num_states);
    return (ep->leftstates)[pos];
}

#ifdef _WIN32
__inline int
#elif __clang__
int
#else
inline int
#endif 
sankoff_get_rightstate (elt_p ep, int pos)
{
    assert(pos<ep->num_states);
    return (ep->rightstates)[pos];
}

value
sankoff_CAML_init_state (value this_elt, value position, value cost) {
    CAMLparam3(this_elt,position,cost);
    elt_p ep;
    Sankoff_elt_custom_val(ep,this_elt);
    int pos, c;
    pos = Int_val(position);
    c = Int_val(cost);
    sankoff_set_state (ep,pos,c);
    CAMLreturn (Val_unit);
}

#ifdef _WIN32
__inline void
#elif __clang__
void
#else
inline void
#endif 
sankoff_set_e (elt_p ep, int pos, int c)
{
    (ep->e)[pos] = c;
    return;
}

#ifdef _WIN32
__inline int
#else
inline int
#endif
sankoff_get_e (elt_p ep, int pos)
{
    return (ep->e)[pos];
}

value
sankoff_CAML_init_e (value this_elt, value position, value cost)
{
    CAMLparam3(this_elt,position,cost);
    elt_p ep;
    Sankoff_elt_custom_val(ep,this_elt);
    int pos, c;
    pos = Int_val(position);
    c = Int_val(cost);
    sankoff_set_e(ep,pos,c);
    CAMLreturn (Val_unit);
}

#ifdef _WIN32
__inline void
#elif __clang__
void
#else
inline void
#endif 
sankoff_set_beta (elt_p ep, int pos, int c)
{
    (ep->beta)[pos] = c;
    return;
}

#ifdef _WIN32
__inline int
#elif __clang__
int
#else
inline int
#endif
sankoff_get_beta (elt_p ep, int pos)
{
    return (ep->beta)[pos];
}

value
sankoff_CAML_init_beta (value this_elt, value position, value cost)
{
    CAMLparam3(this_elt,position,cost);
    elt_p ep;
    Sankoff_elt_custom_val(ep,this_elt);
    int pos, c;
    pos = Int_val(position);
    c = Int_val(cost);
    sankoff_set_beta(ep,pos,c);
    CAMLreturn (Val_unit);
}

#ifdef _WIN32
__inline void
#elif __clang__
void
#else
inline void
#endif 
sankoff_set_m (elt_p ep, int pos, int c)
{
    (ep->m)[pos] = c;
    return;
}

#ifdef _WIN32
__inline int
#elif __clang__
int
#else
inline int
#endif
sankoff_get_m (elt_p ep, int pos)
{
    return (ep->m)[pos];
}

//return 1 if a<b or b is infinity but a is not,
//return 0 if a>=b or if a is infinity 
#ifdef _WIN32
__inline int
#elif __clang__
int
#else
inline int
#endif
cost_less (int a, int b)
{
    if (is_infinity(a)) return 0;
    else if (is_infinity(b)) return 1;
    else return (a<b);
}

//return 1 if a<b or b is infinity but a is not,
//return 0 if a>=b or if a is infinity 
#ifdef _WIN32
__inline int
#elif __clang__
int
#else
inline int
#endif
cost_less_or_equal (int a, int b)
{
    if (is_infinity(a)) return 0;
    else if (is_infinity(b)) return 1;
    else return (a<=b);
}

//return min(a,b)
#ifdef _WIN32
__inline int
#elif __clang__
int
#else
inline int
#endif
cost_min (int a, int b)
{
    if (cost_less(a,b)) return a;
    else return b;
}

#ifdef _WIN32
__inline int
#elif __clang__
int
#else
inline int
#endif
cost_plus (int a, int b)
{
    if ((is_infinity(a))||(is_infinity(b))) return infinity;
    else return(a+b);
}


#ifdef _WIN32
__inline int
#elif __clang__
int
#else
inline int
#endif
cost_minus (int a, int b)
{
    if (is_infinity(a)) {
        assert(!(is_infinity(b)));
        return INT_MAX;
    }
    else if (is_infinity(b)) {
        return 0;
    }
    else return(a-b);
}

//store min(a,*b) to b. Note that b is a pointer. return 1 if b is changed.
#ifdef _WIN32
__inline int
#elif __clang__
int
#else
inline int
#endif
store_min (int a, int * b) {
    if(cost_less(a,*b)) 
    { 
        *b=a;
        return 1;
    }
    else { return 0; };
}


//fill in cost with min cost, idx with which state give us min cost
#ifdef _WIN32
__inline void
#elif __clang__
void
#else
inline void
#endif
sankoff_get_min_state (elt_p ep, int * cost, int * idx)
{
    int i;
    int best=infinity; 
    int thisstate;
    for (i=0;i<ep->num_states;i++)
    {
        thisstate = sankoff_get_state(ep,i);
        if( cost_less(thisstate,best) ) {
            *cost=thisstate; *idx = i; best=thisstate;
        }
    }
    return;
}

value
sankoff_CAML_get_e_array (value a) {
    CAMLparam1(a);
    elt_p ep;
    ep = Sankoff_return_elt(a);
    int num_states = ep->num_states;
    long dims[1];
    dims[0] = num_states;
    CAMLreturn(alloc_bigarray(BIGARRAY_INT32 | BIGARRAY_C_LAYOUT, 1, ep->e, dims));
}

//return min cost between same states, if num_samestates=0, return inf
#ifdef _WIN32
__inline int
#elif __clang__
int
#else
inline int
#endif
sankoff_get_min_cost_between_same_states (int * samestates, int num_samestates, int * tcm, int num_states)
{
    int i;
    int res=infinity;
    int tmp;
    for (i=0;i<num_samestates;i++) {
       tmp = sankoff_return_value(tcm,num_states,num_states,samestates[i],samestates[i]);
       store_min(tmp,&res);
    }
    return res;
}

void
sankoff_canonize (elt_p ep, int * cm)
{
    int i; int j;
    int state; int newe; 
    int newbeta;
    int mincost=infinity,state_w_mincost=0;
    sankoff_get_min_state(ep,&mincost,&state_w_mincost);
    int num_states = ep->num_states;
    for (i=0;i<num_states;i++) {
        state = sankoff_get_state(ep,i);
        newe = (is_infinity(state)) ? infinity : state-mincost;
        sankoff_set_e(ep,i,newe);
    }
    int best; int thise;
    for (i=0;i<num_states;i++) {
        newbeta = infinity;
        for (j=num_states-1;j>=0;j--) {
           thise = sankoff_get_e(ep,j);
           best = (is_infinity(thise)) ? thise : thise + sankoff_return_value(cm,num_states,num_states,i,j);
           store_min(best,&newbeta);
        }
        sankoff_set_beta(ep,i,newbeta);
    }
    return;
}

value
sankoff_CAML_create_eltarr (value is_identity, value taxon_code, value code, value number_of_states, value ecode_bigarr, value states_bigarr, value tcm_bigarr) {
    CAMLparam5(is_identity,taxon_code,code,number_of_states,ecode_bigarr);
    CAMLxparam2(states_bigarr,tcm_bigarr);
    CAMLlocal1(res);
    int num_states;
    num_states = Int_val(number_of_states);
    int tcode = Int_val(taxon_code);
    int iside = Int_val(is_identity);
    int mycode = Int_val(code); 
    int * cost_mat; int dimcm1, dimcm2;
    int * states_arrarr; int dims1, dims2;
    int * ecode_arr; int dim;
    ecode_arr = (int*) Data_bigarray_val(ecode_bigarr);
    dim = Bigarray_val(ecode_bigarr)->dim[0];//number of elts
    states_arrarr = (int*) Data_bigarray_val(states_bigarr);
    dims1 = Bigarray_val(states_bigarr)->dim[0]; //number of elts
    dims2 = Bigarray_val(states_bigarr)->dim[1]; //number of states in each elt
    if (dim!=dims1) failwith ("sankoff.c, size of ecode array != number of charactors");
    if (dims2!= num_states) failwith ("sankoff.c, size of states array != number of states");
    cost_mat = (int*) Data_bigarray_val(tcm_bigarr);
    dimcm1 = Bigarray_val(tcm_bigarr)->dim[0];//number of states
    dimcm2 = Bigarray_val(tcm_bigarr)->dim[1];//number of states
    if ((dimcm1!=dimcm2)||(dimcm1!=dims2)) 
        failwith ("sankoff.c, wrong size of costmat between states");
    eltarr_p neweltarr;
    //alloc struct elt_arr 
    neweltarr = (eltarr_p)calloc(1,sizeof(struct elt_arr));
    neweltarr->code = mycode;
    neweltarr->taxon_code = tcode;
    neweltarr->left_taxon_code = tcode;
    neweltarr->right_taxon_code = tcode;
    neweltarr->sum_cost = 0;
    neweltarr->num_states = dimcm1;
    neweltarr->num_elts = dim;
    neweltarr->is_identity = iside;
    //alloc its pointers
    neweltarr->tcm = (int*)calloc(dimcm1*dimcm2,sizeof(int));
    memcpy(neweltarr->tcm,cost_mat,sizeof(int) * dimcm1 * dimcm2);
    neweltarr->elts = (elt_p)calloc(dim,sizeof(struct elt));
    int i; int j;
    int * states_arr;
    elt_p newelt;
    for (i=0;i<dim;i++)
    {
        newelt = &((neweltarr->elts)[i]);
        assert(newelt!=NULL);
        newelt->ecode = ecode_arr[i];
        newelt->num_states = num_states;
        newelt->states = (int*)calloc( num_states, sizeof(int) );
        newelt->leftstates = (int*)calloc( num_states, sizeof(int) );
        newelt->rightstates = (int*)calloc( num_states, sizeof(int) );
        //for new median_3
        if (median_3_su) { 
            newelt->left_costdiff_mat = (int*)calloc(num_states*num_states,sizeof(int));
            newelt->right_costdiff_mat = (int*)calloc(num_states*num_states,sizeof(int));
        }
        states_arr = sankoff_move_to_line_i(states_arrarr,dims1,dims2,i);
        //the infinity on ocaml side is diff from here, so we pass -1 instead
        //memcpy(newelt->states,states_arr,sizeof(int)*num_states);
        for (j=0;j<num_states;j++) {
            (newelt->states)[j] = ( states_arr[j]==(-1) ) ? infinity : states_arr[j];
        }   
        newelt->beta = (int*)calloc(num_states,sizeof(int));
        newelt->e = (int*)calloc(num_states,sizeof(int));
        newelt->m = (int*)calloc(num_states,sizeof(int));
        sankoff_canonize(newelt,cost_mat);
    }
    res = caml_alloc_custom (&sankoff_custom_operations_eltarr,sizeof (eltarr_p), 1,alloc_custom_max);
    Sankoff_return_eltarr(res) = neweltarr;
    CAMLreturn(res);
}


value
sankoff_CAML_create_eltarr_bytecode (value * argv, int argn){
    return (sankoff_CAML_create_eltarr (argv[0],argv[1], argv[2], argv[3], argv[4], argv[5], argv[6]));
}

//store the shared states between two eltarr. return the number of shared states
#ifdef _WIN32
__inline int
#elif __clang__
int
#else
inline int
#endif
elt_return_shared_states(elt_p ep1, elt_p ep2, int * samestates) {
    int num_states = ep1->num_states;
    int idx = 0;//idx must init to 0
    int i; 
    int e1,e2;
    for (i=0;i<num_states;i++) {
       e1 = sankoff_get_e(ep1,i);
       if (e1==0) {
           e2 = sankoff_get_e(ep2,i);
           if(e2==0) {
               samestates[idx]=i;
               idx++;
           }
       }
    }
    //idx is the true size of array for samestates,
    //when there is no same state, idx should be 0.
    return idx;
}

//return the extra cost by adding a root if the root share states with its
//left/right child, but cost between same states is non-zero.
int
sankoff_elt_get_extra_cost_for_root(elt_p eproot,int * tcm) {
    int * leftchildstates = eproot->leftstates;
    int * rightchildstates = eproot->rightstates;
    int cost=0, beststate=0;
    int num_states = eproot->num_states;
    sankoff_get_min_state(eproot,&cost,&beststate);
    assert(beststate<num_states);

    int beststateleft = leftchildstates[beststate];
    int beststateright = rightchildstates[beststate];

    int extra_cost=0;
    if (beststateleft==beststate) {
        extra_cost = sankoff_return_value(tcm,num_states,num_states,beststate,beststateleft);
    } else if (beststateright==beststate) {
        extra_cost = sankoff_return_value(tcm,num_states,num_states,beststate,beststateright);
    }
    return (extra_cost);
}

int
sankoff_get_extra_cost_for_root(eltarr_p eapRoot) {
    int acc=0;
    int i;
    int num_elts = eapRoot->num_elts;
    int iside = eapRoot->is_identity;
    if (!(iside==1)){
        for (i=0;i<num_elts;i++) {
            acc = acc + sankoff_elt_get_extra_cost_for_root(&((eapRoot->elts)[i]),eapRoot->tcm);
        }
    }
    return acc;
}

//Calculates array m -- the extra cost when joining a new branch/node R to existing
//subtree on edge(nodeA,nodeD), creating a new node M between node A and D. 
int
sankoff_elt_dist_2 (elt_p epD, elt_p epA, elt_p epR, int * tcm) {
    int num_states = epD->num_states;
    int * shared_states_rd = (int*)calloc(num_states,sizeof(int));
    int * shared_states_ra = (int*)calloc(num_states,sizeof(int));
    int num_srd = elt_return_shared_states(epD,epR,shared_states_rd);
    int num_sra = elt_return_shared_states(epA,epR,shared_states_ra);
    if ((num_srd>0)||(num_sra>0)) {
        int costrd = sankoff_get_min_cost_between_same_states(shared_states_rd,num_srd,tcm,num_states);
        int costra = sankoff_get_min_cost_between_same_states(shared_states_ra,num_sra,tcm,num_states);
        free(shared_states_rd);
        free(shared_states_ra);
        return (costrd<costra) ? costrd : costra;
    } else {
        //free these two first
        free(shared_states_rd); free(shared_states_ra);
        //if we don't have array m yet, this part will cost us O(n^3)+O(n^2) time, n is the number of states. if array m is already there, we just need O(n^2) time.
        int y,s,i,x;//idx for loops
        //fill in array m if it hasn't been taken cared of.
        if (epD->m_already_set) {}//do nothing
        else {//this part only be called once for each node.
            int tcm_ix, tcm_is, tcm_ss;
            int d_si;
            //fill in array m for nodeD, arrays and integers for each s loop
            int * tbeta = (int*)calloc(num_states,sizeof(int));
            //remember what is state i when D(s,i) reach mininum.
            int * besti_arr = (int*)calloc(num_states,sizeof(int));
            int min_tbeta;  int min_dsi;  int besti_arr_size;
            for (s=0;s<num_states;s++) {
                min_dsi = infinity;//D.s = min[D(s,i)], i~[]
                besti_arr_size=0; //reset idx of besti_arr to 0
                for (i=0;i<num_states;i++) {
                    min_tbeta=infinity;
                    for (x=0;x<num_states;x++) {
                        tcm_ix = sankoff_return_value(tcm,num_states,num_states,i,x);
                        tbeta[x] = cost_plus(tcm_ix, sankoff_get_beta(epD,x));
                        store_min(tbeta[x],&min_tbeta);
                    }
                    tcm_is = sankoff_return_value(tcm,num_states,num_states,i,s);
                    d_si = sankoff_get_e(epA,i) + tbeta[s] - min_tbeta;
                    if (cost_minus(d_si,min_dsi)) {
                        besti_arr[besti_arr_size]=i;
                        besti_arr_size ++;
                    }
                    store_min(d_si,&min_dsi);
                }
                if(int_array_is_mem(besti_arr,besti_arr_size,s)) {
                    tcm_ss = sankoff_return_value(tcm,num_states,num_states,s,s);
                    min_dsi = min_dsi + tcm_ss;
                }
                sankoff_set_m(epD,s,min_dsi);
            }
            free(tbeta);
            free(besti_arr);
            epD->m_already_set = 1;
        }
        int best=infinity;
        int tcm_sy; int costRM;
        int e_R_y; int m_D_s;
        for (y=0;y<num_states;y++) {
            e_R_y = sankoff_get_e(epR,y);
            for (s=0;s<num_states;s++) {
                m_D_s = sankoff_get_m(epD,s);
                tcm_sy = sankoff_return_value(tcm,num_states,num_states,s,y);
                costRM = cost_plus(tcm_sy,e_R_y);
                costRM = cost_plus(m_D_s,costRM);
                store_min(costRM,&best);
            }
        }
        return best;
    }
}

#ifdef _WIN32
__inline int
#elif __clang__
int
#else
inline int
#endif
sankoff_dist_2 (eltarr_p eapD,eltarr_p eapA, eltarr_p eapR) {
    int acc=0;
    int i;
    int num_elts;
    num_elts = eapD->num_elts;
    for (i=0;i<num_elts;i++) {
        acc = acc + sankoff_elt_dist_2(&((eapD->elts)[i]),&((eapA->elts)[i]),&((eapR->elts)[i]),eapD->tcm);
    }
    return acc;
}

/* This algorithm is taken from Goloboff 1998.  Calculate D', then the final E
 * value from that.*/
void
sankoff_elt_median_3(elt_p epA, elt_p epN, elt_p epL, elt_p epR, elt_p newepN, int * tcm,int is_left_child) {
    int num_states = epN->num_states;
    //copy epN to newepN, we are going to update array e later
    sankoff_clone_elt (newepN,epN);
    // for new  median3
    if (median_3_su) { 
        int i, s;
        int eAi;
        int di, min_di;
        int * cost_diff_mat; int cost_diff;
        if(is_left_child) cost_diff_mat = epA->left_costdiff_mat;
        else cost_diff_mat = epA->right_costdiff_mat;
        for (s=0;s<num_states;s++) {
            min_di = infinity;//reset min to inf
            for (i=0;i<num_states;i++) {
                eAi = sankoff_get_e(epA,i);
                cost_diff = sankoff_return_value(cost_diff_mat,num_states,num_states,i,s);
                di = cost_plus(cost_diff,eAi);
                store_min(di,&min_di);
            }
            sankoff_set_e(newepN,s,min_di);
        }
        //end new median3 
    } else {
        //old median3 start
        int i,s;
        //store L.beta.s + R.beta.s to betasum
        int * betasum = (int*)calloc(num_states,sizeof(int));
        int x;
        int betaL,betaR;
        for (x=0;x<num_states;x++) {
            betaL = sankoff_get_beta(epL,x);
            betaR = sankoff_get_beta(epR,x);
            betasum[x] = cost_plus(betaL,betaR);
        }
        //min_tbeta_i : fill in min(t.i.x + L.beta.x + R.beta.x) for each i
        int * min_tbeta_i = (int*)calloc(num_states,sizeof(int));
        int tix;
        int tbeta;
        int min_betasum;
        for (i=0;i<num_states;i++) {
            min_betasum = infinity; //reset min to inf
            for (x=0;x<num_states;x++) {
                tix = sankoff_return_value(tcm,num_states,num_states,i,x);
                tbeta = cost_plus(tix,betasum[x]);
                store_min(tbeta,&min_betasum);
            }
            min_tbeta_i[i] = min_betasum;
        }
        //calc D_s_i for each pair of (s,i)
        int eAi;
        int min_di;
        int d_i;
        int tis;
        for (s=0;s<num_states;s++) {
            min_di = infinity;//reset min to inf
            for (i=0;i<num_states;i++) {
                eAi = sankoff_get_e(epA,i);
                tis = sankoff_return_value(tcm,num_states,num_states,i,s);
                d_i = cost_plus(eAi,tis);
                d_i = cost_plus(d_i,betasum[s]);
                d_i = cost_minus(d_i,min_tbeta_i[i]);
                store_min(d_i,&min_di);
            }
            //fill in best D_s_i as E for each s
            sankoff_set_e(newepN,s,min_di);
        }
        free(betasum);
        free(min_tbeta_i);
        //end of old median3 
    }
    return;
}


void
sankoff_median_3(eltarr_p eapA,eltarr_p eapN, eltarr_p eapL, eltarr_p eapR, eltarr_p neweapN) {
    int i;
    int num_states, num_elts;
    num_states = eapN->num_states;
    num_elts = eapN->num_elts;
    if (sankoff_is_leaf_node(eapN)) {
        //alloc neweapN's pointers 
        sankoff_init_eltarr (neweapN, num_states, num_elts, eapN->code, eapN->taxon_code,eapN->left_taxon_code,eapN->right_taxon_code, eapN->tcm, eapN->is_identity);
        for (i=0;i<num_elts;i++) {
            sankoff_create_empty_elt(&((neweapN->elts)[i]),num_states,-1);
            sankoff_clone_elt(&((neweapN->elts)[i]),&((eapN->elts)[i]));
        }
        return;
    } else {
        assert(eapA->code == eapN->code);
        assert(eapA->num_elts == eapN->num_elts);
        assert(eapA->num_elts == eapL->num_elts);
        assert(eapA->num_elts == eapR->num_elts);
        int is_left_child = sankoff_is_left_or_right_child(eapN,eapA);
        sankoff_init_eltarr (neweapN, num_states, num_elts, eapN->code, eapN->taxon_code,eapL->taxon_code,eapR->taxon_code, eapN->tcm, eapN->is_identity);
        for (i=0;i<num_elts;i++) {
            sankoff_create_empty_elt(&((neweapN->elts)[i]),num_states,-1);
            sankoff_elt_median_3 (&((eapA->elts)[i]),&((eapN->elts)[i]),&((eapL->elts)[i]),&((eapR->elts)[i]),&((neweapN->elts)[i]),eapA->tcm,is_left_child);
        }
        return;
    }
}



/* Calculates the median between to characters a and b. a and b have to be
* homologous (same code), and must also have the same transformation cost matrix
* associated. If they are homologous, it must also be the case that they hold
* the same number of valid states. */
int
sankoff_elt_median(elt_p ep1, elt_p ep2, int * tcm, elt_p newep) {
    int debug = 0;
    //ecode was init with (-1) before this function
    newep->ecode = ep1->ecode;
    int num_states = ep1->num_states;
    int i,j,k;
    int costij,costik;
    int costj,costk;
    int median_cost;
    int best_cost=infinity;//best cost for all possible states
    if(debug) printf("sankoff_elt_median\n");
    //find best states
    //this way we need only O(n^2) 
    int bestj_i=0, bestk_i=0;
    int best_cost_i_j, best_cost_i_k;
    int anychange;
    int * lcm; int * rcm; 
    int tmp;
    for (i=0;i<num_states;i++) {
        best_cost_i_j = infinity;
        //for new median3 speedup
        if (median_3_su) {
            lcm = sankoff_move_to_line_i(newep->left_costdiff_mat,num_states,num_states,i);
        }
        for(j=0;j<num_states;j++) {
            costij = sankoff_return_value(tcm,num_states,num_states,i,j);
            costj = sankoff_get_state(ep1,j);
            tmp = cost_plus(costij,costj);
            //for new median3 speedup
            if (median_3_su) { lcm[j]= tmp; }
            anychange=store_min(tmp,&best_cost_i_j);
            if(anychange) bestj_i=j;
        }
        //*for new median3 speedup
        if (median_3_su) {
            for(j=0;j<num_states;j++) {
                lcm[j] = cost_minus(lcm[j],best_cost_i_j);
            } 
        }
        best_cost_i_k = infinity;
        //*for new median3 speedup
        if (median_3_su) {
            rcm = sankoff_move_to_line_i(newep->right_costdiff_mat,num_states,num_states,i); 
        }
        for(k=0;k<num_states;k++) {
            costik = sankoff_return_value(tcm,num_states,num_states,i,k);
            costk = sankoff_get_state(ep2,k);
            tmp = cost_plus(costik,costk);
            //*for new median3 speedup
            if (median_3_su) {rcm[k] = tmp; }
            anychange=store_min(tmp,&best_cost_i_k);
            if(anychange) bestk_i=k;
        }
        //for new median3 speedup
        if (median_3_su) { 
        for(k=0;k<num_states;k++) {
            rcm[k] = cost_minus(rcm[k],best_cost_i_k);
        }
        }
        median_cost = cost_plus(best_cost_i_j,best_cost_i_k);
        sankoff_set_state(newep,i,median_cost);
        sankoff_set_leftstate(newep,i,bestj_i);
        sankoff_set_rightstate(newep,i,bestk_i);
        store_min(median_cost,&best_cost);
    }
    //fill in array e
    int coste;
    for (i=0;i<num_states;i++) {
        coste = cost_minus(sankoff_get_state(newep,i),best_cost); 
        sankoff_set_e(newep,i,coste);
    }
    //fill in array beta
    int costbetai; int costej; 
    for (i=0;i<num_states;i++) {
        costbetai = infinity;
        for (j=0;j<num_states;j++) {
            costej = sankoff_get_e (newep,j);
            tmp = sankoff_return_value(tcm,num_states,num_states,i,j);
            tmp = cost_plus(tmp,costej);
            store_min (tmp,&costbetai);
        }
        sankoff_set_beta(newep,i,costbetai);
    }
    int cost=infinity,cost_idx=0;
    sankoff_get_min_state(newep,&cost,&cost_idx);
    if(debug) { printf("end of elt median\n"); sankoff_print_elt(newep,1,1,0,1); }
    return cost;
}

/* Calculates the distance between two characters a and b 
* Note that we only calculate the _added_ distance */
int
sankoff_elt_distance(elt_p ep1, elt_p ep2, int * tcm, elt_p newep) {
    int debug = 0;
    if (debug) {
        printf("sankoff_elt_distance,call elt_median first\n");
        sankoff_print_elt(ep1,1,1,1,1);
        sankoff_print_elt(ep2,1,1,1,1);
    }
    if (debug) {
        printf("end of elt median\n"); 
        sankoff_print_elt(newep,1,1,1,1);
    }
    int state_w_mincost=0, mincost=infinity;
    sankoff_get_min_state(newep,&mincost,&state_w_mincost);
    if (debug) { printf("mincost = %d, state_w_mincost = %d,",mincost,state_w_mincost); fflush(stdout); }
    int leftstate = sankoff_get_leftstate(newep,state_w_mincost);
    int rightstate = sankoff_get_rightstate(newep,state_w_mincost);
    if (debug) printf("leftstate = %d, rightstate = %d,",leftstate,rightstate);
    int leftcost = sankoff_get_state(ep1,leftstate);
    int rightcost = sankoff_get_state(ep2,rightstate);
    if (debug) printf("cost = mincost =%d - leftcost = %d - rightcost = %d\n", mincost,leftcost,rightcost);
    return (mincost - leftcost - rightcost);
}

void
sankoff_median(int median_node_tcode,eltarr_p eap1,eltarr_p eap2, eltarr_p neweltarr) {
    int debug = 0;
    if (debug) { printf("===== sankoff_median =====\n on eap1:\n");
        sankoff_print_eltarr(eap1,1,0,0,0);
        printf(" and eap2:\n");
        sankoff_print_eltarr(eap2,1,0,0,0);
    }
    if(debug && eap1->code != eap2->code) printf("code of eap1 = %d/%d, of eap2 = %d/%d\n",
            eap1->code,eap1->taxon_code,eap2->code,eap2->taxon_code);
    assert(eap1->code == eap2->code);
    assert(eap1->num_elts == eap2->num_elts);
    long int sum_cost=0, elt_cost;
    int i; int num_states, num_elts;
    num_states = eap1->num_states;
    num_elts = eap1->num_elts;
    //remember which one is our left child , which one is right.
    sankoff_init_eltarr (neweltarr, num_states, num_elts, eap1->code, median_node_tcode, eap1->taxon_code,eap2->taxon_code, eap1->tcm, eap1->is_identity);
    for (i=0;i<num_elts;i++)
    {
        sankoff_create_empty_elt(&((neweltarr->elts)[i]),num_states,-1);
        elt_cost = sankoff_elt_median(&((eap1->elts)[i]),&((eap2->elts)[i]),eap1->tcm, &((neweltarr->elts)[i]));
        sum_cost = sum_cost + (long int)elt_cost;
    }
    neweltarr->sum_cost = sum_cost;
    if (debug) { printf("====== return median ====== \n");
        sankoff_print_eltarr(neweltarr,1,1,1,1);
    }
    return;
}

#ifdef _WIN32
__inline int
#elif __clang__
int
#else
inline int
#endif
sankoff_get_best_child_state (eltarr_p eap, int getleftchild) {
    assert(eap->num_elts == 1);//this function is for fixed states only
    elt_p ep = &(eap->elts)[0];
    int cost=0, idx=0;
    sankoff_get_min_state(ep,&cost,&idx);
    if (getleftchild) 
        return (ep->leftstates)[idx];
    else
        return (ep->rightstates)[idx];
}


//this function is for fixed states only
value
sankoff_CAML_get_best_child_state (value a, value ch_tcode) {
    CAMLparam2(a,ch_tcode);
    eltarr_p eap;
    eap = Sankoff_return_eltarr(a);
    //int position = Int_val(p);
    int child_taxon_code = Int_val(ch_tcode);
    int leftORright=0;
    if (eap->left_taxon_code == child_taxon_code) leftORright=1;
    else if(eap->right_taxon_code == child_taxon_code) leftORright = 0;
    else failwith("sankoff,get_best_child_state, neither left nor right");
    int beststate = sankoff_get_best_child_state(eap,leftORright);
    CAMLreturn(Val_int(beststate));
}


int
sankoff_distance(eltarr_p eap1, eltarr_p eap2, eltarr_p neweltarr) {
    int acc=0;
    int i;
    int num_states, num_elts;
    num_states = eap1->num_states;
    num_elts = eap1->num_elts;
    sankoff_init_eltarr (neweltarr, num_states, num_elts, eap1->code, 0, eap1->taxon_code,eap2->taxon_code,eap1->tcm,eap1->is_identity);
    for (i=0;i<num_elts;i++) {
        sankoff_create_empty_elt(&((neweltarr->elts)[i]),num_states,-1);
        sankoff_elt_median(&((eap1->elts)[i]),&((eap2->elts)[i]),eap1->tcm, &((neweltarr->elts)[i]));
        acc = acc + sankoff_elt_distance (&((eap1->elts)[i]),&((eap2->elts)[i]),eap1->tcm, &((neweltarr->elts)[i]));
    }
    return acc;
}

        
value
sankoff_CAML_median(value code, value a, value b) {
   CAMLparam3(code,a,b);
   CAMLlocal1(res);
   eltarr_p eap1;
   eltarr_p eap2;
   eap1 = Sankoff_return_eltarr(a);
   eap2 = Sankoff_return_eltarr(b);
   eltarr_p neweltarr;
   neweltarr = (eltarr_p)calloc(1,sizeof(struct elt_arr));
   sankoff_median(Int_val(code),eap1,eap2,neweltarr);
   res = caml_alloc_custom (&sankoff_custom_operations_eltarr,sizeof (eltarr_p), 1,alloc_custom_max);
   Sankoff_return_eltarr(res) = neweltarr; 
   CAMLreturn(res);
}

value
sankoff_CAML_median_3(value a, value n, value l, value r) {
   CAMLparam4(a,n,l,r);
   CAMLlocal1(res);
   eltarr_p eapA;
   eltarr_p eapN;
   eltarr_p eapL;
   eltarr_p eapR;
   eapA = Sankoff_return_eltarr(a);
   eapN = Sankoff_return_eltarr(n);
   eapL = Sankoff_return_eltarr(l);
   eapR = Sankoff_return_eltarr(r);
   eltarr_p neweltarr;
   neweltarr = (eltarr_p)calloc(1,sizeof(struct elt_arr));
   sankoff_median_3(eapA,eapN,eapL,eapR,neweltarr);
   res = caml_alloc_custom (&sankoff_custom_operations_eltarr,sizeof (eltarr_p), 1,alloc_custom_max);
   Sankoff_return_eltarr(res) = neweltarr;
   CAMLreturn(res);
}



/* [distance a b] return the sankoff distance. Note that it calls
* [elt_distance], which will call [elt_median]. if you need median and distance,
* don't call two functions seperately, use [distance_and_median] instead*/
value
sankoff_CAML_distance(value a, value b) {
   CAMLparam2(a,b);
   eltarr_p eap1;
   eltarr_p eap2;
   eap1 = Sankoff_return_eltarr(a);
   eap2 = Sankoff_return_eltarr(b);
   eltarr_p neweltarr;
   neweltarr = (eltarr_p)calloc(1,sizeof(struct elt_arr)); 
   int res = sankoff_distance(eap1,eap2,neweltarr);
   free_eltarr(neweltarr);
   free(neweltarr);
   CAMLreturn(Val_int(res));
}


value
sankoff_CAML_dist_2(value a, value b,value c) {
   CAMLparam3(a,b,c);
   eltarr_p eapD;
   eltarr_p eapA;
   eltarr_p eapR;
   eapD = Sankoff_return_eltarr(a);
   eapA = Sankoff_return_eltarr(b);
   eapR = Sankoff_return_eltarr(c);
   int res = sankoff_dist_2(eapD,eapA,eapR);
   CAMLreturn(Val_int(res));
}

value
sankoff_CAML_get_extra_cost_for_root(value a) {
    CAMLparam1(a);
    CAMLlocal1(ores);
    eltarr_p eapRoot;
    eapRoot = Sankoff_return_eltarr(a);
    int res = sankoff_get_extra_cost_for_root(eapRoot);
    ores = Val_int(res);
    CAMLreturn(ores);
}
