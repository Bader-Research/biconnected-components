/**********************************************************************************
 NOTES:  This listranking code is modified to return the prefix value of the last element in the list after the List Ranking operation is done. This prefix value is the total no of leaves for the rake operation. It does not return the list head pointer.

Last changed July 30 2002

**********************************************************************************/

#ifndef _LISTRANK_H
#define _LISTRANK_H

#include "simple.h"

typedef long LDATA;

#define LISTRANK_TYPE     LDATA
#define LISTRANK_OPERATOR +
#define LISTRANK_IDENTITY 0

typedef struct list_d
{
  LISTRANK_TYPE prefix;
  LDATA succ;
} list_t;

/* n = number of elements in the list;
   k = multplier of THREADS, where number of sublists (s) == k*THREADS
   List = the List
   //return value = the index of the list head ptr
   return value = total no of leaves given by prefix value of last element in the list (modified for rake operation)
*/

LDATA list_ranking(LDATA n, int k, list_t *List, THREADED);

#endif

