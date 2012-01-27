/*************************************************************************************/
/*  C program to implement a Queue Implementation using Linked List                  */
/*  modify from www.jave2s.com 
 *  by lin lhong@amnh.org */
/*************************************************************************************/

# include <stdio.h>
# include <stdlib.h>
#include "queue_with_linkedlist.h"

# define DEBUG 0



//just make sure front and rear are set to NULL, and size is 0 for empty queue.
void init_empty_queue (q_t emptyq)
{
    emptyq -> front = NULL;
    emptyq -> rear = NULL;
    emptyq -> size = 0;
}



//maybe we can reuse nodes
void enqueue (q_t thisq, QUEUE_DATA value)
{
   n_t* front = &(thisq->front); 
   n_t* rear = &(thisq->rear);
   n_t temp;
   temp=(struct node *)malloc(sizeof(struct node));
   if(temp==NULL)
   {
      printf("ERROR in queue.c , No Memory available \n");
      exit(0);
   }
   thisq->size ++;
   temp->data = value;
   temp->link=NULL;
   if(*rear == NULL)
   {
      *rear = temp;
      *front = *rear;
   }
   else
   {
      (*rear)->link = temp;
      *rear = temp;
   }
   if (DEBUG) 
   { printf ("enqueue %d to queue,size=%d\n",value,thisq->size); fflush(stdout); }
}

//void dequeue(struct node **front, struct node **rear, QUEUE_DATA *value)
void dequeue(q_t thisq, QUEUE_DATA *value)
{
   n_t* front = &(thisq->front); 
   n_t* rear = &(thisq->rear);
   n_t temp;
   thisq->size --;
   if((*front == *rear) && (*rear == NULL))
   {
      printf("ERROR cannot dequeue element from empty queue\n");
      exit(0);
   }
   *value = (*front)->data;
   temp = *front;
   *front = (*front)->link;
   if(*rear == temp)
   *rear = (*rear)->link;
   free(temp);
   if (DEBUG) 
   { printf ("dequeue front(=%d) from queue,size=%d\n",*value,thisq->size); fflush(stdout); }

}

void peekqueue(q_t thisq, QUEUE_DATA *value)
{
   n_t* front = &(thisq->front); 
   n_t* rear = &(thisq->rear);
   if((*front == *rear) && (*rear == NULL))
   {
      printf("ERROR cannot peek into empty queue\n");
      exit(0);
   }
   *value = (*front)->data;
}

#ifdef _WIN32
__inline int 
#else
inline int 
#endif
is_emptyqueue (q_t thisq)
{//maybe we should check front==rear && rear==NULL as well
    if ( thisq->size == 0 ) 
    {
        if ((thisq->front==NULL)&&(thisq->rear==NULL))
            return 1;
        else 
        {
            printf ("ERROR:front and rear !=NULL in empty queue\n");
            exit(0);
        }
    }
    else return 0;
}

void clear_queue (q_t thisq)
{
   QUEUE_DATA tmp;
   while (thisq->size > 0)
   {
      dequeue(thisq,&tmp);
   }
   thisq->front = NULL;
   thisq->rear = NULL;
   if (thisq->size!=0) 
   {
       printf ("ERROR: queue size should be 0 at the end of clear_queue\n");
       exit(0);
   }
}

//transfer every thing in q1 to q2 
void transfer_queue(q_t q1, q_t q2)
{
    //empty q2 first
    clear_queue(q2);
    //copy everything
    q2->size = q1->size;
    q2->front = q1->front;
    q2->rear = q1->rear;
    //reset q1's pointer and size
    q1->size=0;
    q1->front=NULL;
    q1->rear=NULL;
}

/*
int main()
{
   //n_t front=NULL; n_t rear = NULL;
   QUEUE_DATA value;

   q_t myqueue = (q_t) malloc(sizeof(struct Queue));
   init_empty_queue(myqueue);

   enqueue(myqueue,1);
   enqueue(myqueue,10);
   enqueue(myqueue,100);
   enqueue(myqueue,1000);

   //dequeue(&front,&rear,&value);
   dequeue(myqueue,&value);
   printf("The value dequeued is %d\n",value);

   dequeue(myqueue,&value);//dequeue(&front,&rear,&value);
   printf("The value dequeued is %d\n",value);

   dequeue(myqueue,&value);//dequeue(&front,&rear,&value);
   printf("The value dequeued is %d\n",value);

   dequeue(myqueue,&value);//dequeue(&front,&rear,&value);
   printf("The value dequeued is %d\n",value);

   dequeue(myqueue,&value);//dequeue(&front,&rear,&value);
   printf("The value dequeued is %d\n",value);

   return 0;

}*/

