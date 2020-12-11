#include "simple.h"
#include "graph.h"


/*===============================================================================*/
/*the spanning tree without lock, without phase 2, and El.workspace will be modified*/
/*modified to handle edgelist that only have one copy for each edge, that is, no two anti-parallel edge*/
int spanning_tree_CRCW(V* graph,E* El, int nVertices,int n_edges,THREADED)
{
  int *p_changed,*buffer,*D,*assign,*done,changed,i,j,n,iter_count=0;
  V * tree;
  E *El_tmp;

#if DEBUG_TIMING
  hrtime_t start,end;
  double interval1=0,interval2=0,interval3=0;
#endif

#if DEBUG_WAIT
  hrtime_t start1,end1;
  double interval_wait=0;
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
#if DEBUG_WAIT
      start1= gethrtime();
#endif

      node_Barrier();

#if  DEBUG_WAIT
      end1=gethrtime();
      interval_wait+=end1-start1;
#endif

      /*1. phase one*/

#if DEBUG_TIMING
      start= gethrtime();
#endif

      on_one_thread (*p_changed)=0;
      changed=0;

#if DEBUG_WAIT
      start1= gethrtime();
#endif
      node_Barrier();

#if  DEBUG_WAIT
      end1=gethrtime();
      interval_wait+=end1-start1;
#endif

      // pardo(i,0,nVertices,1){
      //	assign[i]=0;
      //	D[i]=graph[i].v_attribute;
      //}

#if DEBUG_WAIT
      start1= gethrtime();
#endif

      node_Barrier();

#if  DEBUG_WAIT
      end1=gethrtime();
      interval_wait+=end1-start1;
#endif
    pardo(n,0,n_edges,1)
	{
	  if(El[n].workspace==1) continue; /*this edge has been added to tree*/
	  i=El[n].v1;
	  j=El[n].v2;

	  if(D[j]<D[i] && D[i]==D[D[i]])
	    {
	      done[D[i]]=0;
	      assign[D[i]]=MYTHREAD;/*here all threads that want to change the DDi value should contend for it, one of them will win*/
	    }	  
	    
	}

    node_Barrier();

    pardo(n,0,n_edges,1)
	{
	  if(El[n].workspace==1) continue;
	  i=El[n].v1;
	  j=El[n].v2;
      
	  if(D[j]<D[i] && D[i]==D[D[i]] && assign[D[i]]==MYTHREAD && done[D[i]]==0)
	    {      
	      changed=1;
	      //graph[D[i]].v_attribute=D[j];
	      D[D[i]]=D[j];
		  done[D[i]]=1;
	      El[n].workspace=1;
#if DEBUG_CORRECTNESS
	      edge_count++;
#endif
	    }
	  j=El[n].v1;
	  i=El[n].v2;
      
	}

    if(changed) (*p_changed)=1;
#if DEBUG_WAIT
    start1= gethrtime();
#endif
      node_Barrier();

#if  DEBUG_WAIT
    end1=gethrtime();
    interval_wait+=end1-start1;
#endif

#if DEBUG_TIMING
      end = gethrtime();
      interval1+=end-start;
      printf("this iteration for phase 1 on thread %d is %f\n",MYTHREAD,interval1/1000000000);
#endif
      
     /*3.phase 3*/

#if DEBUG_GRAPH
      printf("phase 3:====\n");
#endif

#if DEBUG_TIMING
      start = gethrtime();
#endif


#if DEBUG_WAIT
      start1=gethrtime();
#endif

      node_Barrier();

#if DEBUG_WAIT    
      end1 = gethrtime();
      interval_wait+=end1-start1;
#endif
	
    pardo(i,0,nVertices,1)
	{
#if 0
	  while(graph[i].v_attribute!=graph[graph[i].v_attribute].v_attribute)
	  graph[i].v_attribute=graph[graph[i].v_attribute].v_attribute;
#endif
	  while(D[i]!=D[D[i]]) D[i]=D[D[i]];
    }

#if DEBUG_WAIT
      start1=gethrtime();
#endif
      node_Barrier();

#if DEBUG_WAIT    
      end1 = gethrtime();
      interval_wait+=end1-start1;
#endif

#if DEBUG_TIMING
      end = gethrtime();
      interval3+=end-start;
#endif

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

    on_one printf("number of iterations:%d\n",iter_count);
#if DEBUG_TIMING
    on_one_thread{
      printf("total time for phase 1 :%f s\n", interval1/1000000000);
      printf("total time for phase 2 :%f s\n", interval2/1000000000);
      printf("total time for phase 3 :%f s\n", interval3/1000000000);
    }
#endif


#if DEBUG_WAIT
    printf("Barrier waiting time is %f \n", interval_wait/1000000000);
#endif
    

#if DEBUG_CORRECTNESS
    printf("Total edges got is %d\n",edge_count);
#endif

    node_free(p_changed,TH);
    node_free(buffer, TH);
    node_free(D, TH);
    node_free(assign, TH);
    node_free(done, TH);
}

