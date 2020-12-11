#include "simple.h"
#include "graph.h"

#define DEBUG_CORRECTNESS 0

/*Nov 7,03: ordering of modifying ddi and done is changed in the second run*/
/*===============================================================================*/
/*AS: CRCW rooted spanning tree*/

int rooted_spanning_tree(int* Parent,E* El, int nVertices,int n_edges,THREADED)
{
  int *p_changed,*buffer,*D,*assign,*done,changed,i,j,n,iter_count=0;
  int p,c,t;
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
  D=node_malloc(sizeof(int)*nVertices,TH); 
  assign=node_malloc(sizeof(int)*nVertices,TH); /*used to check which processor should graft*/ 
  done=node_malloc(sizeof(int)*nVertices,TH);

  pardo(i,0,nVertices,1) /*initialize the D values*/
    {     
      done[i]=0;
      D[i]=i;
    }

  node_Barrier();
  while((*p_changed)==1)
    {
      iter_count++;

      node_Barrier();

      /*1. phase one*/

#if DEBUG_TIMING
      start= gethrtime();
#endif

      on_one_thread (*p_changed)=0;
      changed=0;
      node_Barrier();

    pardo(n,0,n_edges,1)
	{
	  if(El[n].workspace==1) continue; /*this edge has been added to tree*/
	  i=El[n].v1;
	  j=El[n].v2;

	  if(D[j]<D[i] && D[i]==D[D[i]])
	    {
	      done[D[i]]=0;
	      assign[D[i]]=MYTHREAD;/*here all threads that want to change the D[D[i]] value should contend for it, one of them will win*/
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
	      D[D[i]]=D[j];
		  done[D[i]]=1;
		  if(Parent[i]==i ) Parent[i]=j;
		  else {
		  	/* reverse the parent pointers for all nodes on the path to root */
			p=Parent[i];
			c=i;
		  	Parent[i]=j;
			while(Parent[p]!=p)
			{
				t = Parent[p];
				Parent[p]=c;
				c = p;
				p = t;
			}
			Parent[p]=c;
				
		  }
	      El[n].workspace=1;
#if DEBUG_CORRECTNESS
		  printf("THREAD %d: %d -->%d\n",MYTHREAD,i,j);
	      edge_count++;
#endif
	    }
	}

    if(changed) (*p_changed)=1;
      node_Barrier();
	
#if DEBUG_TIMING
      end = gethrtime();
      interval1+=end-start;
      printf("this iteration for phase 1 on thread %d is %f\n",MYTHREAD,interval1/1000000000);
#endif
      
     /*3.phase 3*/

#if DEBUG_TIMING
      start = gethrtime();
#endif

      node_Barrier();
	
    pardo(i,0,nVertices,1)
	{
	  while(D[i]!=D[D[i]]) D[i]=D[D[i]];
    }

      node_Barrier();


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


#if DEBUG_CORRECTNESS
    printf("Total edges got is %d\n",edge_count);
#endif

    node_free(p_changed,TH);
    node_free(buffer, TH);
    node_free(D, TH);
    node_free(assign, TH);
    node_free(done, TH);
}

