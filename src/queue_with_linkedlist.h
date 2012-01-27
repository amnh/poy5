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

#define QUEUE_DATA int

struct node
{
   QUEUE_DATA data;
   struct node *link;
};
typedef struct node * n_t;

struct Queue 
{
    int size;
    n_t front;
    n_t rear;
};
typedef struct Queue * q_t;


void init_empty_queue (q_t emptyq);

void enqueue (q_t thisq, QUEUE_DATA value);

void dequeue (q_t thisq, QUEUE_DATA *value);

void peekqueue(q_t thisq, QUEUE_DATA *value);


#ifdef _WIN32
__inline int 
#else
inline int 
#endif
is_empty (q_t thisq);


void transfer_queue(q_t q1, q_t q2);
