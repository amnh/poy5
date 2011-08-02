
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
    DIRECTION_MATRIX * gapnumarr;//max gap number we could have in alignment
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
};

typedef struct newkkmat * newkkmat_p;

#ifdef _WIN32
__inline int 
#else
inline int 
#endif
newkk_algn (const seqt s1, const seqt s2, MAT_SIZE s1_len, MAT_SIZE s2_len, int go, const cmt c, newkkmat_p m);


