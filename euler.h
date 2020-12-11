
#ifndef _EULER_H_
#define _EULER_H_

#include <sys/types.h>
#include "simple.h"
#include <thread.h>
#include "graph.h"

typedef struct euler_struct {
	int v1,v2;
	int next;
	int twin;
	int value;
	int tmp_buff;
	} Euler_E;
	

#endif /*_EULER_H_*/
