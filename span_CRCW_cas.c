#include "simple.h"
#include "graph.h"

#define DEBUG_CORRECTNESS 1 

extern inline int CASW(int * a, int e, int n);

int spanning_tree_CRCW_cas(V* graph,E* El, int nVertices,int n_edges,THREADED)
{
  int *p_changed,*buffer,*D,*assign,*done,changed,i,j,n,iter_count=0;
  V * tree;
  E *El_tmp;

#if DEBUG_TIMING
  hrtime_t start,end;
  double interval1=0,interval2=0,interval3=0;
#endif

#if DEBUG_CORRECTNESS
  int edge_count=0;
#endif
  
  p_changed=node_malloc(sizeof(int),TH);
  on_one_thread (*p_changed)=1;

  buffer=node_malloc(sizeof(int)*nVertices,TH);
  D=node_malloc(sizeof(int)*nVertices,TH); /*used for Di*/
  assign=node_malloc(sizeof(int)*nVertices,TH); /*used to check which processor should do this*/ 
  done=node_malloc(sizeof(int)*nVertices,TH);

  pardo(i,0,nVertices,1) /*initialize the D values, Jaja p217*/
    {     
      done[i]=0;
      graph[i].v_attribute=i;
      D[i]=i;
    }

  node_Barrier();
  while((*p_changed)==1)
    {
      iter_count++;
      node_Barrier();

      /*1. phase one*/
      on_one_thread (*p_changed)=0;
      changed=0;

      node_Barrier();

    pardo(n,0,n_edges,1)
	{
	  if(El[n].workspace==1) continue;
	  i=El[n].v1;
	  j=El[n].v2;
      
	  if(D[j]<D[i] && D[i]==D[D[i]] && CASW(&done[D[i]],0,MYTHREAD+1)==0)
	    {
	      changed=1;
		  El[n].workspace=1;
	      D[D[i]]=D[j];
#if DEBUG_CORRECTNESS
	      edge_count++;
#endif
	    }
	}

    if(changed) (*p_changed)=1;
      node_Barrier();
      
     /*3.phase 3*/	
    pardo(i,0,nVertices,1)
	{
	  while(D[i]!=D[D[i]]) D[i]=D[D[i]];
    }

      node_Barrier();

      /*phase 4: shrink the edge list to reduce the size*/
#if 0
      if(iter_count%2==0){
	El_tmp=shrink_edgeL(El, &n_edges,D,TH);
	node_free(El,TH);
	El=El_tmp;
      }
#endif

    } /*while*/

#if DEBUG_CORRECTNESS
    node_Barrier();
    pardo(i,0,nVertices,1)
    {
#if 0
	 if(graph[i].v_attribute!=graph[0].v_attribute) {
	   printf("ERROR in SPANNING TREE\n");
	   break;
	 }
#endif
	 if(D[i]!=D[0] ){
	   printf("ERROR in SPANNING TREE\n");
	   break;
	 }
    }
#endif

    printf("number of iterations:%d\n",iter_count);


#if DEBUG_CORRECTNESS
	edge_count = node_Reduce_i(edge_count,SUM,TH);
    on_one printf("Total edges got is %d\n",edge_count);
#endif

    node_free(p_changed,TH);
    node_free(buffer, TH);
    node_free(D, TH);
    node_free(assign, TH);
    node_free(done, TH);
}

