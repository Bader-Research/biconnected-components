
#include "simple.h"
#include "graph.h"
#include "stack.h"
#define NANO 1000000000
#define THRESHOLD 50

/* A balanced breadth-first search*/
int span_gw(V* graph,int nVertices,THREADED)
{
#define S_POINTS (THREADS*THREADS*THREADS*2)
#define MYCOLOR (MYTHREAD+1)

  hrtime_t start,end;
  double interval;

  int * color, first_time, work_to_steal,myroot,b,start_sr,start_sl, *stack,top,count,bottom=-1;
  int i,j,n,root,walks,r,counter=0,visited,neighbor,ret=0,n_neighbors;
  double power;
  int ** stack_M, **top_M, **bottom_M;
  int * finished;
  unsigned int seed=MYCOLOR;
  int * count_M;

  stack_M = node_malloc(THREADS*sizeof(int *),TH);
  top_M = node_malloc(THREADS*sizeof(int *),TH);
  bottom_M=node_malloc(THREADS*sizeof(int *),TH);

  stack_M[MYTHREAD]=malloc(nVertices*sizeof(int));
  stack=stack_M[MYTHREAD];
  top_M[MYTHREAD]=&top;
  bottom_M[MYTHREAD]=&bottom;
  count_M=node_malloc(THREADS*sizeof(int),TH);

  color=node_malloc(sizeof(int)*nVertices,TH); 

  start=gethrtime();
  pardo(i,0,nVertices,1){
    color[i]=0;
  }
  end=gethrtime();
  interval=end-start;
  on_one printf("The time used for setting up is %f \n",interval/NANO);
 
  bottom=-1;
  top=-1;
  count=0;

  node_Barrier();

  seed=gethrtime()/(MYTHREAD+1)+MYTHREAD;
  if(ret==0) root=(rand_r(&seed)%nVertices);
  
  /*lets select a point to start in the graph*/
  on_one_thread {
    color[root]=MYCOLOR;
    graph[root].parent=root;
    myroot=root;
    j=0;
    push(myroot,stack_M[j],top_M[j]);
    j=(j+1)%THREADS;
    i=0;
    for(i=0;i<S_POINTS;i++)
      {
	n_neighbors=graph[myroot].n_neighbors;
	r=rand_r(&seed);
	if(r%2==0){
	  for(r=0;r<n_neighbors;r++)
	    {
	      visited++;
	      n=graph[myroot].my_neighbors[r];
	      if(color[n]==0){
		graph[n].parent=myroot;
		color[n]=MYCOLOR;
		myroot=n;
		count++;
		push(myroot,stack_M[j],top_M[j]);
		/*printf("sub tree : %d \n",myroot);*/
		j=(j+1)%THREADS;
		break;
	      } 
	    }
	  if(r==n_neighbors){
	    r=(rand_r(&seed)%n_neighbors);
	    myroot=graph[myroot].my_neighbors[r];
	  }
	}
	else {
	  for(r=n_neighbors-1;r>=0;r--)
	    {
	      visited++;
	      n=graph[myroot].my_neighbors[r];
	      if(color[n]==0){
		graph[n].parent=myroot;
		color[n]=MYCOLOR;
		myroot=n;
		count++;
		push(myroot,stack_M[j],top_M[j]);
		/*printf("sub tree : %d \n",myroot);*/
		j=(j+1)%THREADS;
		break;
	      }     
	    }
	  if(r<0){
	    r=(rand_r(&seed)%n_neighbors);
	    myroot=graph[myroot].my_neighbors[r];
	  }
	}		
      }
    end=gethrtime();
    interval=end-start;
  }
	       
  node_Barrier();
  on_one printf("Now walking...\n");
  
  start=gethrtime();
  first_time=1;
  work_to_steal=0;
  start_sr=MYTHREAD; /*when I am out of work, where do i start to search. r towards right, l towards left*/
  start_sl=MYTHREAD;
  
  while(first_time || work_to_steal)
  {
      /*if(work_to_steal) printf("stealing work\n");*/
      while(!is_empty( stack, &top, &bottom))
	  {
	    n=pop(stack,&top,bottom);
	    if(n==-1){
	      printf("THREAD %d:stack overflow\n",MYTHREAD);
	      bottom=-1;
	      top=-1;
	      break;
	    }
	    visited+=graph[n].n_neighbors;
	    for(i=0;i<graph[n].n_neighbors;i++)
	    {
	      neighbor=graph[n].my_neighbors[i];
	      if(color[neighbor]==0) {/*found new frontier*/
		    color[neighbor]=1;
		    graph[neighbor].parent=n;
		    push(neighbor,stack,&top);
		    count++;
	      }
	    } 
	  }
      /*printf("Thread %d:done with this stack, current count is %d \n", MYTHREAD, count);*/

      if(first_time) first_time=0;
      work_to_steal=0;
      if(MYTHREAD%2==0)
	  {
	    for(j=0;j<THREADS;j++)
	    {
	      i=(j+start_sr)%THREADS; /*start searching in circular from my neighbor*/
	      if(i==MYTHREAD) continue;

	      n=*(top_M[i]);
	      b=*(bottom_M[i]);	
	      if(n-b<THRESHOLD) continue;	      	       
	      if(count>nVertices/THREADS) /*I did my share*/
		     r=b+(n-b)/THREADS;
	      else r=b+max((n-b)/THREADS,min(nVertices/THREADS-count,(n-b)*3/4));
	      /*printf("Thread%d: thread %d's stack is %d tall\n",MYTHREAD,i,n-b);*/
	      if(r<(*top_M[i])) {
		   (*bottom_M[i])=max(r-1,-1);
		   /*printf("THREAD %d: I am taking %d elements\n",MYTHREAD, r-b);*/
		   work_to_steal=1;		 
		   while((r--)>max(0,b))
		     push(stack_M[i][r],stack,&top);
		   break;		
	      }
	    }
	} else{
	  for(j=THREADS-1;j>0;j--)
	    {
	      i=(j+start_sl)%THREADS;
	      if(i==MYTHREAD) continue;
	      n=*(top_M[i]);
	      b=(*bottom_M[i]);
	      if(n-b<THRESHOLD) continue;	     
	      /*printf("Thread%d: thread %d's stack is %d tall\n",MYTHREAD,i,n-b);*/
	      if(count>nVertices/THREADS) /*I did my share*/
		r=b+(n-b)/THREADS;
	      else r=b+max((n-b)/THREADS,min(nVertices/THREADS-count,(n-b)*3/4));
	      if(r<(*top_M[i])) {
		(*bottom_M[i])=max(r-1,-1);
		/*printf("THREAD %d: I am taking %d elements\n",MYTHREAD, r-b);*/
		work_to_steal=1;
		while((r--)>max(b,0))
		  push(stack_M[i][r],stack,&top);
		break;
	      } 
	    }
	}
      start_sr=(start_sr+1)%THREADS;
      start_sl=(start_sl-1+THREADS)%THREADS;

    }
  count_M[MYTHREAD]=count;
  end=gethrtime();
  interval=end-start;
  printf("THREAD %d: I AM DONE..., time used: %f\n",MYTHREAD,interval/1000000000);
  node_Barrier();
  printf("Thread %d count is %d, visited is %d\n",MYTHREAD, count,visited);  

  on_one_thread{
    int max=0, min=nVertices;
    for(i=0;i<THREADS;i++)
      {
	if(count_M[i]>max) max=count_M[i];
	if(count_M[i]<min) min=count_M[i];
      }
    printf("METRICS===span_gw:The difference between counts is %d\n",max-min);
  }
  node_Barrier();

  node_free(count_M,TH);
  node_free(color, TH);
  node_free(stack_M,TH);
  node_free(top_M,TH);
  free(stack);
}




