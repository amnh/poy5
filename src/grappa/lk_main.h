#ifndef LK_MAIN_H
#define LK_MAIN_H

#ifdef USE_CONCORDE

#include "machdefs.h"
#include "linkern.h"
#include "util.h"
#include "kdtree.h"
#include "edgegen.h"
#include "macrorus.h"
#include "tsp.h"


int chlinkern ( int ncount, int **adj, int *tour, int *incycle,
                int *outcycle );
int greedylk ( int ncount, int **adj, int *tour, int *incycle,
               int *outcycle );
#endif

#endif
