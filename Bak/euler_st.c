#include "simple.h"
#include "graph.h"
#include "types.h"
#include "listrank.h"

#define debug 0
#define DEBUG_CORRECTNESS 0
int prefix_sum(int *buff, int n, THREADED);

/*===============================================================================*/
/*Run the spanning tree routine twice to get 
 (1) parent relationship
 (2) euler-tour setup : in the first run, each vertex will remember which vertex that its supervertex is grafted to.(i,j)
 cause D[i] grafted to j.
  in the second run, we will check the reverse direction (the twin of the edges being
 selected), if it causes grafting, then we now it is the twin of the previously selected edges. 
 graft_to is a n*logn array that records in each iteration where a vertex is grafted to.
 */

int st_euler(V* graph,E* El, int nVertices,int n_edges,int * graft_to,int * Buff,THREADED)
{
  int *p_changed,*D,*assign,*done,changed,i,j,n,iter_count=0,label;
  V * tree;
  E *El_tmp;
  
  hrtime_t start,end;
  double interval1=0,interval2=0,interval3=0,interval;
  int edge_count=0;

  p_changed=node_malloc(sizeof(int),TH);
  on_one_thread (*p_changed)=1;
	 
  D=node_malloc(sizeof(int)*nVertices,TH); /*used for Di*/
  assign=node_malloc(sizeof(int)*nVertices,TH); /*used to check which processor should do this*/ 
  done=node_malloc(sizeof(int)*nVertices,TH);

  pardo(i,0,nVertices,1) /*initialize the D values, Jaja p217*/
    {     
      done[i]=0;
      D[i]=i;
    }

  node_Barrier(); 
  while((*p_changed)==1)
  {
      node_Barrier();

      /*1. phase one*/

      on_one_thread (*p_changed)=0;
      changed=0;
	  label = (iter_count+1)*nVertices;
      node_Barrier();
      
    pardo(n,0,n_edges,1)
	{
	  if(El[n].workspace==1) continue; /*this edge has been added to tree*/
	  i=El[n].v1;
	  j=El[n].v2;
	  if(D[j]<D[i] && D[i]==D[D[i]])
	    {
		  if(!changed) changed =1;
		  done[D[i]]=0;
	      assign[D[i]]=MYTHREAD;/*here all threads that want to change the DDi value should contend for it, one of them will win*/
	    }	  
	}
    node_Barrier();

    if(changed) (*p_changed)=1;
	node_Barrier();
	if(*p_changed==0) break;
	
    pardo(n,0,n_edges,1)
	{
	  if(El[n].workspace==1) continue;
	  i=El[n].v1;
	  j=El[n].v2;
      
	  if(D[j]<D[i] && D[i]==D[D[i]] && assign[D[i]]==MYTHREAD && done[D[i]]==0)
	    {
		  /*printf("edge (%d,%d)\n",j,i);*/  	  	      
	      D[D[i]]=D[j];
		  done[D[i]]=1;
	      El[n].workspace=1;
		  graft_to[i]=label+j; /*record where a vertex is grafted,and make it unique for each round*/
	      edge_count++;
		  
	    }
	}
	node_Barrier();

	n=node_Reduce_i(edge_count,SUM,TH);
	on_one printf("Thread %d:edge_count is %d\n",MYTHREAD,n);
	/*pardo(i,0,nVertices,1)
	  printf("gt[%d]=%d \n", i,graft_to[i]);*/
	pardo(n,0,n_edges,1)
	{
	  if(El[n].workspace==1) continue;
	  j=El[n].v1;
	  i=El[n].v2;    
	  /*printf("gt[%d]=%d ", i,graft_to[i]);*/
	  if( graft_to[i]==(label+j) )
	    {			
		  /*printf("edge (%d,%d)\n",i,j);*/	  
	      El[n].workspace=1;
	      edge_count++;
		 
	    }
	}
	
    node_Barrier();
	n=node_Reduce_i(edge_count,SUM,TH);
	on_one printf("Thread %d:edge_count is %d\n",MYTHREAD,n);   
	pardo(i,0,nVertices,1)
	{
	  while(D[i]!=D[D[i]]) D[i]=D[D[i]];
    }
    node_Barrier();

#if 1	
	if(iter_count>0) {
	pardo(i,0,n_edges,1)
	{
		if(El[i].workspace!=1 && D[El[i].v1]!=D[El[i].v2]) Buff[i]=1;
		else Buff[i]=0;
	}
	node_Barrier();
	prefix_sum(Buff,n_edges,TH);
	on_one printf("after prefix sum\n");
	n= Buff[n_edges-1];
	if(n<=n_edges/2){
		on_one printf(" n is %d\n",n);
		El_tmp = node_malloc(sizeof(E)*n,TH);
		pardo(i,0,n_edges,1)
		{
			if((i==0 && Buff[i]==1) || (i!=0 && Buff[i]!=Buff[i-1]))
			El_tmp[Buff[i]-1]=El[i];
		}
		node_Barrier();
		n_edges=n;
		El=El_tmp;
	}
	}
#endif
	
	iter_count++;
	node_Barrier();
	
  } /*while*/

	
#if DEBUG_CORRECTNESS
    node_Barrier();
    pardo(i,0,nVertices,1)
    {
	 if(D[i]!=D[0] ){
	   printf("ERROR in SPANNING TREE\n");
	   break;
	 }
    }
#endif

    on_one printf("number of iterations:%d\n",iter_count);    

#if DEBUG_CORRECTNESS
	edge_count = node_Reduce_i(edge_count,SUM,TH);
    printf("Total edges got is %d\n",edge_count);
#endif

    node_free(p_changed,TH);
    node_free(D, TH);
    node_free(assign, TH);
    node_free(done, TH);	
}



#if 0

/* Set up the structure used for the euler technique. Assuming edges adjacent to one vertex are 
	in consequtive locations*/

void construct_Euler_path ( ET* treeEL, int n_edges,THREADED)
{

#define M 8192
#define P 8
#define SAMPLES 8

	ET* treeEL_tmp,*treeEL_tmp2;
	ET* treeEL_sorted;
	int * Buff,i,next,twin,orig_index1,orig_index2;
	ET max_dummy, max_dummy_m1;
	
	hrtime_t start,end;
	double interval;
	
	Buff = node_malloc(sizeof(int)*n_edges,TH);
	
	treeEL_tmp = node_malloc(sizeof(ET)*(n_edges+4*(n_edges/M+THREADS*THREADS)),TH);
	treeEL_tmp2 = node_malloc(sizeof(ET)*(n_edges+4*(n_edges/M+THREADS*THREADS)),TH);
	
	pardo(i,0,n_edges,1)
	{	
		if(i==0 || treeEL[i].v1!=treeEL[i-1].v1) Buff[i]=i;
		else Buff[i]=0;	
		treeEL_tmp[i]=treeEL[i];
		treeEL_tmp[i].workspace=i; /*workspace is now my index in the orig edge list*/
	}
	node_Barrier();
	
	prefix_sum_max(Buff,n_edges,TH);
	node_Barrier();

#if DEBUG
	
	on_one {
		for(i=0;i<n_edges;i++)
			printf(" %d   ", Buff[i]);
	printf("\n");
	}
#endif
	
	pardo(i,0,n_edges,1)
	{
		if(i<n_edges-1 && treeEL[i].v1==treeEL[i+1].v1) 
			treeEL[i].next=i+1;
		else treeEL[i].next=Buff[i]; /* makes "next" circular */
	}

#if DEBUG	
	on_one {
		printf("next:\n");
		for(i=0;i<n_edges;i++)
			printf(" %d   ", treeEL[i].next);
	printf("\n");
	}
#endif
	
	node_Barrier();
	
	max_dummy.v1=2*n_edges+2;
	max_dummy.v2=2*n_edges+2;
	
	max_dummy_m1.v1=2*n_edges+1;
	max_dummy_m1.v2=2*n_edges+1;
	
	on_one init_sort(max_dummy,max_dummy_m1);
	node_Barrier();
	start = gethrtime();
	sample_sort(n_edges,M,P,SAMPLES,treeEL_tmp,treeEL_tmp2,&treeEL_sorted,TH); /*according to (min(v1,v2),max(v1,v2))*/
	end = gethrtime();
	interval = end -start;
	on_one printf("METRICS: time used for sorting is %f s\n", interval/1000000000);
	
	node_Barrier();
	
	if(treeEL_tmp !=treeEL_sorted) {
		treeEL_tmp2=treeEL_tmp;
		treeEL_tmp=treeEL_sorted;
	}
	
#if DEBUG	
	on_one {
		printf("Sorted results:\n");
		for(i=0;i<n_edges;i++)
			printf("(%d,%d,%d) \n",treeEL_tmp[i].v1,treeEL_tmp[i].v2,treeEL_tmp[i].workspace);
		}
#endif

# if 0	
		printf("check sort,treeEL_tmp is %d\n",treeEL_tmp);
		pardo(i,0,n_edges-1,1)
		{
			int a1,a2, b1,b2;
			a1 = min(treeEL_tmp[i].v1,treeEL_tmp[i].v2);
			a2 = min(treeEL_tmp[i+1].v1,treeEL_tmp[i+1].v2);
			b1 = max(treeEL_tmp[i].v1,treeEL_tmp[i].v2);
			b2 = max(treeEL_tmp[i+1].v1,treeEL_tmp[i+1].v2);
			if(	!(a1 <a2 || (a1==a2 && b1<=b2)) ) printf("sorting error %d \n",i);
		}

		pardo(i,0,n_edges,2)
		{
			if(treeEL_tmp[i].v1!=treeEL_tmp[i+1].v2) printf("not in pairs\n");
		}
		printf("check done\n");
#endif
	
	node_Barrier();
	
	pardo(i,0,n_edges,2)
	{
		orig_index1 = treeEL_tmp[i].workspace;
		orig_index2 = treeEL_tmp[i+1].workspace;
		treeEL[orig_index1].twin=orig_index2;
		treeEL[orig_index2].twin=orig_index1;
	}
	
	node_Barrier();
	
#if DEBUG	
	printf("set twin done\n");
	on_one {
		for(i=0;i<n_edges;i++)
			printf(" %d   \n", treeEL[i].twin);
	printf("\n");
	}
#endif
	
	pardo(i,0,n_edges,1)
	{
		twin = treeEL[i].twin;
		next = treeEL[twin].next;
		treeEL[next].pred=i;
	}
	node_Barrier();
	node_free(treeEL_tmp,TH);
	node_free(treeEL_tmp2,TH);
	node_free(Buff,TH);
		
}

#endif
