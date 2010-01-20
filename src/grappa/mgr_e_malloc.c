/* e_malloc.c
*    Memory allocation with error abort if necessary.
*
* Copyright (C) 2001-2006 The Regents of the University of California
* by Glenn Tesler
*
* See file COPYRIGHT for details.
*****************************************************************************
* This file is part of GRIMM.
*
* GRIMM is free software; you can redistribute it and/or modify
* it under the terms of the GNU General Public License, Version 2,
* dated June 1991, as published by the Free Software Foundation.
*
* GRIMM is distributed in the hope that it will be useful, but
* WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
* General Public License for more details.
*
* You should have received a copy of the GNU General Public License along
* with this program; if not, write to the Free Software Foundation, Inc.,
* 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
*/

/* Last modified on Tue Aug 1, 2006, by Glenn Tesler
*/

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include "mgrstructs.h"
#include "mgr_e_malloc.h"
//mcstructs.h under MGR share some structures with structs.h under GRAPPA


void *e_malloc(int s, char *msg) {
  void *p;

  p = malloc(s);
  if (p == (void *)NULL) {
    fflush(stdout);
    fclose(stdout);

    if (errno == ENOMEM)
      fprintf(stderr,"ERROR: Insufficient memory on this system.");
    else if (errno == EAGAIN)
      fprintf(stderr,"ERROR: Insufficient memory at this time.");
    else {
      fprintf(stderr,"ERROR: Unexpected error during malloc:");
      perror("mcdist");
    }
    fprintf(stderr," %s NULL\n", msg);

    exit(-1);
  }
  return p;
}


