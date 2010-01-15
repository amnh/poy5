/* check.c
*
* Copyright (C) 2004-2006 Genome Insitute of Singapore
* Copyright (C) 2000-2002 University of Southern California
*
* See file COPYRIGHT for details. 
*****************************************************************************
* This file is part of MGR.
*
* MGR is free software; you can redistribute it and/or modify
* it under the terms of the GNU General Public License, Version 2,
* dated June 1991, as published by the Free Software Foundation.
*
* MGR is distributed in the hope that it will be useful, but
* WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
* General Public License for more details.
*
* You should have received a copy of the GNU General Public License along
* with this program; if not, write to the Free Software Foundation, Inc.,
* 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
*/

/* Last modified on Wed Aug 2, 2006, by Guillaume Bourque
*/

#include <stdio.h>
#include <stdlib.h>
#include "mgrstructs.h"



void *ckalloc(int amount)
{
	void *p;

	if(amount == 0)	{
		amount = (unsigned) 100;
	}
	if ((p = (void *) calloc( (unsigned) amount, 1)) == NULL)	{
		printf("Ran out of memory.\n");
                printf("There may be errors as follows:\n");
                printf("1) Not enough memory.\n");
                printf("2) The ARRAY may be overrode.\n");
                printf("3) The wild pointers.\n");
                exit(-1);
	}
	return(p);
}
