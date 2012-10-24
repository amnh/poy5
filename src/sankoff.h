#include <caml/custom.h>

#define Sankoff_return_elt(x) (*((elt_p*)Data_custom_val(x)))
#define Sankoff_elt_custom_val(to_asgn,a)  to_asgn = Sankoff_return_elt(a)

#define Sankoff_return_eltarr(x) (*((eltarr_p*)Data_custom_val(x)))
#define Sankoff_eltarr_custom_val(to_asgn,a) to_asgn = Sankoff_return_eltarr(a)

#define infinity (INT_MAX/4)
#define is_infinity(x) (x>=infinity)

struct elt {
    int ecode;
    int num_states;
    /* The array of states, each index in the array corresponds to a single
     * state and the integer in the array corresponds to the cost of the that
     * state given the neighbors of the character in the node that contains it in
     * the tree.
     */
    int * states;
    /* The following values come from Goloboff 1998: Tree Searches under Sankoff
     * Parsimony. */
    /* This is a downpass-calculated intermediate value */
    int * beta;
    /* e.(s) is the preliminary added cost for using state s
     * (uppass-calculated) */
    int * e;
    /* M.(s) is a cached value; it depends on E and beta, so we reset it during
     * the downpass. */
    int m_already_set;
    int * m;
    /* Scratch area for the diagnosis, it is safe to do funcky things with it
    * because we only use it during the diagnosis */
    //int * best_states;
    //we need to remember when we pick state i(i might be different from best state for current node), 
    //which state of left child give us best cost
    int * leftstates;
    //we need to remember which state of right child give us best cost
    int * rightstates;
    //n^2 matrix
    //(i,j) when this elt choose state i, left child choose state j~[0,n-1], what's the
    //extra cost comparing to choosing best j for left child.
    int * left_costdiff_mat;
    //similar matrix for right child
    int * right_costdiff_mat;
};

typedef struct elt * elt_p;

struct elt_arr {
    int code;//mycode from sankCS.of_parser
    //taxon_code of leaf node is from sankCS.of_parser, positive numbers from 1
    //taxon_code of median node is passed from Ocaml side, it's a negative
    //number,increased as number of internal node increase
    int taxon_code; 
    int left_taxon_code; //taxon code of left child
    int right_taxon_code; //taxon code of right child
    int num_states;//each elt should have the same state
    int * tcm; //cost matrix between states;
    int is_identity; //1 if there is no cost between same states.
    int num_elts;//number of elts
    elt_p elts;//array of elts
    long int sum_cost;//sum of cost from each elt
};

typedef struct elt_arr * eltarr_p;

int alloc_custom_max = 10000;
