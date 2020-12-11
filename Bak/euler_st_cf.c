#include "simple.h"
#include "graph.h"
#include "types.h"
#include "listrank.h"

#define debug 0


int prefix_sum(int *buff, int n, THREADED);

/*===============================================================================*/
/* Try to make it more cache friendly by making a array of structs instead of
many seperate arrays. */

int st_euler_cf(V* graph,E* El, int nVertices,int n_edges,THREADED)
{
  typedef struct aux_r {
   int D, assign, done, graft_to;} aux_t;
   
  int *p_changed,*D,*assign,*done,changed,i,j,n,iter_count=0,logn;
  V * tree;
  E *El_tmp;
  aux_t * Aux;
  
  hrtime_t start,end;
  double interval1=0,interval2=0,interval3=0,interval;

#if DEBUG_CORRECTNESS
  int edge_count=0;
#endif

  p_changed=node_malloc(sizeof(int),TH);
  on_one_thread (*p_changed)=1;

  Aux= node_malloc(sizeof(aux_t)*nVertices,TH);
  
 #if 0 
  D=node_malloc(sizeof(int)*nVertices,TH); /*used for Di*/
  assign=node_malloc(sizeof(int)*nVertices,TH); /*used to check which processor should do this*/ 
  done=node_malloc(sizeof(int)*nVertices,TH);
 #endif
  
  
  pardo(i,0,nVertices,1) /*initialize the D values, Jaja p217*/
    {     
      Aux[i].done=0;
      Aux[i].D=i;
    }

  node_Barrier(); 
  while((*p_changed)==1)
  {
      node_Barrier();

      /*1. phase one*/

      on_one_thread (*p_changed)=0;
      changed=0;
      node_Barrier();

    pardo(n,0,n_edges,1)
	{
	  if(El[n].workspace==1) continue; /*this edge has been added to tree*/
	  i=El[n].v1;
	  j=El[n].v2;
	  if(Aux[j].D<Aux[i].D && Aux[i].D==Aux[Aux[i].D].D)
	    {
	      Aux[Aux[i].D].done=0;
	      Aux[Aux[i].D].assign=MYTHREAD;/*here all threads that want to change the DDi value should contend for it, one of them will win*/
	    }	  
	}

    node_Barrier();

    pardo(n,0,n_edges,1)
	{
	  if(El[n].workspace==1) continue;
	  i=El[n].v1;
	  j=El[n].v2;
      
	  if(Aux[j].D<Aux[i].D && Aux[i].D==Aux[Aux[i].D].D && Aux[Aux[i].D].assign==MYTHREAD && Aux[Aux[i].D].done==0)
	    {
		  
	      Aux[Aux[i].D].done=1;
	      changed=1;
	      Aux[Aux[i].D].D=Aux[j].D;
	      El[n].workspace=1;
#if DEBUG_CORRECTNESS
	      edge_count++;
#endif
		  Aux[Aux[i].D].graft_to=j; /*record where a super vertex is grafted*/
	    }
	}

    if(changed) (*p_changed)=1;
	iter_count++;
    node_Barrier();
    pardo(i,0,nVertices,1)
	{
	  while(Aux[i].D!=Aux[Aux[i].D].D) Aux[i].D=Aux[Aux[i].D].D;
    }

      node_Barrier();
  } /*while*/

  /* The second run*/
  node_Barrier();
  pardo(i,0,nVertices,1) 
  {     
      Aux[i].done=0;
      Aux[i].D=i;
  }
  iter_count=0;
  on_one *p_changed=1;
  node_Barrier();
  while((*p_changed)==1)
  {
      node_Barrier();

      /*1. phase one*/

      on_one_thread (*p_changed)=0;
      changed=0;
      node_Barrier();

      pardo(n,0,n_edges,1)
	  {
	    if(El[n].workspace==1) continue; /*this edge has been added to tree*/
	    j=El[n].v1;
	    i=El[n].v2;
	    if(Aux[j].D<Aux[i].D && Aux[i].D==Aux[Aux[i].D].D && Aux[Aux[i].D].done==0 && Aux[Aux[i].D].graft_to==j)
	    {
	      Aux[Aux[i].D].done=1;
	      changed=1;
	      Aux[Aux[i].D].D=Aux[j].D;
	      El[n].workspace=1;
#if DEBUG_CORRECTNESS
	      edge_count++;
#endif
	    }	  
	  }

      node_Barrier();

      if(changed) (*p_changed)=1;
	  iter_count++;
      node_Barrier();
      pardo(i,0,nVertices,1)
	  {
	    while(Aux[i].D!=Aux[Aux[i].D].D) Aux[i].D=Aux[Aux[i].D].D;
      }

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
    printf("Total edges got is %d\n",edge_count);
#endif

    node_free(p_changed,TH);
    node_free(Aux,TH);
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
