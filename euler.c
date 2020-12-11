#include "simple.h"
#include "graph.h"
#include "types.h"
#include "listrank.h"

#define debug 0

int prefix_sum(int *buff, int n, THREADED);
void sample_sort(int size,int M,int p,int samples,TYPE *Buffer1,
                           TYPE *Buffer2,TYPE **out_array_ptr,THREADED);

/*spanning tree algorithm only finds each edge once without 
  the twin edges, here label the twin edges also */

void label_twin_edges(E* EL, V* G, int n_edges, THREADED)
{
	int i,j,v1,v2,start, *Found;
	Found = node_malloc(sizeof(int)*n_edges,TH);
	
	pardo(i,0,n_edges,1)
		if(EL[i].workspace==1) Found[i]=1;
		else Found[i]=0;
		
	/* Just copy the tree edge list in the reverse direction and sort 
	might be a more theoretically sound solution*/
		
	node_Barrier();
	pardo(i,0,n_edges,1)
	{
		if(Found[i]==1) {
			v2 = EL[i].v2;
			v1 = EL[i].v1;
			start = G[v2].edge_start;
			for(j=start;EL[j].v1==v2;j++)
				if(EL[j].v2==v1) {
					EL[j].workspace=1;
					break;
				}
		
		}
	}
	node_Barrier();
	node_free(Found,TH);
}


/* Pick out all the spanning tree edges in the graph edge list rep */
ET*  pick_tree_edges(E* EL,int n_edges,THREADED)
{
	
	int * Buff,i,j,new_edges;
	ET* newEL;
	
	Buff = node_malloc(sizeof(int)*n_edges,TH);
	
	pardo(i,0,n_edges,1)
		Buff[i]=EL[i].workspace;
	
	node_Barrier();
	
	prefix_sum(Buff,n_edges,TH);
	node_Barrier();
	
	new_edges=Buff[n_edges-1];
	/*on_one printf("new_edges is %d\n",new_edges);*/
	
	newEL = node_malloc(sizeof(ET)*new_edges,TH);
	
	pardo(i,0,n_edges,1)
		if(EL[i].workspace==1){
			newEL[Buff[i]-1].v1=EL[i].v1;
			newEL[Buff[i]-1].v2=EL[i].v2;
			newEL[Buff[i]-1].workspace=0;
	}
	node_Barrier();
	node_free(Buff,TH);
	return(newEL);
}


/* Pick out all the spanning tree edges in the graph edge list rep . Here El only has one copy for each edge*/
ET*  pick_tree_edges_s(E* EL,int n_edges,THREADED)
{
	
	int * Buff,i,j,new_edges;
	ET* newEL;
	
	Buff = node_malloc(sizeof(int)*n_edges,TH);
	
	pardo(i,0,n_edges,1)
		Buff[i]=EL[i].workspace;
	
	node_Barrier();
	
	prefix_sum(Buff,n_edges,TH);
	node_Barrier();
	
	new_edges=Buff[n_edges-1];
	on_one printf("new_edges is %d\n",new_edges);
	
	newEL = node_malloc(sizeof(ET)*2*new_edges,TH);
	
	pardo(i,0,n_edges,1)
		if(EL[i].workspace==1){
			newEL[2*(Buff[i]-1)].v1=EL[i].v1;
			newEL[2*(Buff[i]-1)].v2=EL[i].v2;
			newEL[2*(Buff[i]-1)].workspace=0;
			newEL[2*(Buff[i]-1)+1].v2=EL[i].v1;
			newEL[2*(Buff[i]-1)+1].v1=EL[i].v2;
			newEL[2*(Buff[i]-1)+1].workspace=0;
	}
	node_Barrier();
	node_free(Buff,TH);
	return(newEL);
}
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

/* the list-ranking code we use actually start from the beginnig of the list*/
/* a->b->c->d is a list defined by successor. Yet the value got is v(a)=a, v(b)=
a+b,etc*/

void Tree_to_list(ET* treeEL, list_t * List,int n_edges,THREADED)
{
        int i;
        pardo(i,0,n_edges,1){
                List[treeEL[i].pred].succ=i;
        }
}

/* root the tree at the first vertex */
void Euler_root_tree(int *Parent,ET* treeEL,list_t * List,int n_edges,THREADED)
{
	int i,j;

	pardo(i,0,n_edges,1){
		List[i].prefix=1;	
	}
	on_one List[treeEL[0].pred].succ=-1;
	
	node_Barrier();	
#if debug	
	on_one {
		for(i=0;i<n_edges;i++)
			printf(" (%d,%d,%d ) ",i,List[i].succ, List[i].prefix);
	}
	
	node_Barrier();
#endif

	list_ranking(n_edges,8,List,TH);	
	node_Barrier();

    pardo(i,0,n_edges,1)
		treeEL[i].value = List[i].prefix;
	
	node_Barrier();		
	pardo(i,0,n_edges,1)
	{
		if(treeEL[i].value < treeEL[treeEL[i].twin].value) 
			Parent[treeEL[i].v2]=treeEL[i].v1;
	}

#if debug	
	on_one {
		printf("The values:\n");
		for(i=0;i<n_edges;i++)
			printf(" %d  ", treeEL[i].value);
		printf("\n");
	}
#endif	
	node_Barrier();
}

/* get the size of each subtree*/
void Euler_size_tree(int *Parent,int * Size,ET* treeEL,list_t * List,int n_edges,THREADED)
{
	int i,j;
	
	List[treeEL[0].pred].succ=-1;	
	node_Barrier();
	pardo(i,0,n_edges,1){
		if(treeEL[i].v1==Parent[treeEL[i].v2])
			List[i].prefix=0;
		else List[i].prefix=1;
	}
	node_Barrier();
		
	list_ranking(n_edges,8,List,TH);
	node_Barrier();

    pardo(i,0,n_edges,1)
		treeEL[i].value = List[i].prefix;
	
	node_Barrier();		
	pardo(i,0,n_edges,1)
	{
		if(Parent[treeEL[i].v1]==treeEL[i].v2) {
			Size[treeEL[i].v1]= treeEL[i].value-treeEL[treeEL[i].twin].value; 
		}
	}

#if debug	
	on_one {
		printf("The values:\n");
		for(i=0;i<n_edges;i++)
			printf(" %d  ", treeEL[i].value);
		printf("\n");
	}
#endif	
	node_Barrier();
	on_one Size[0]=n_edges;
}

/* get the size of each subtree*/
void Euler_preorder(int *Parent,int * Preorder,ET* treeEL,list_t * List,int n_edges,THREADED)
{
	int i,j;
	
	List[treeEL[0].pred].succ=-1;	
	node_Barrier();
	pardo(i,0,n_edges,1){
		if(treeEL[i].v1==Parent[treeEL[i].v2])
			List[i].prefix=1;
		else List[i].prefix=0;
	}
	node_Barrier();
		
	list_ranking(n_edges,8,List,TH);
	node_Barrier();

    pardo(i,0,n_edges,1)
		treeEL[i].value = List[i].prefix;
	
	node_Barrier();		
	pardo(i,0,n_edges,1)
	{
		if(Parent[treeEL[i].v2]==treeEL[i].v1) {
			Preorder[treeEL[i].v2]= treeEL[i].value; 
		}
	}

#if debug	
	on_one {
		printf("The values:\n");
		for(i=0;i<n_edges;i++)
			printf(" %d  ", treeEL[i].value);
		printf("\n");
	}
#endif	
	node_Barrier();
	on_one Preorder[0] = 0;
}

/*get the low of each node, jaja 238.*/
int Euler_get_lowhigh(E* El, int *Parent,int *order,int n,int n_edges, int root,int *low,int *high, THREADED)
{
  int i,j,pos,ret,done=0;
  int *finished,*low_,*high_,*lock_array;
  hrtime_t start,end;
  float interval;
  int * D, *leaf;
  
  
  low_=node_malloc(sizeof(int)*n,TH); 
  high_=node_malloc(sizeof(int)*n,TH); 
  lock_array = node_malloc(sizeof(int)*n,TH);
  
  pardo(i,0,n,1)
    {
      low[i]=order[i];
      high[i]=order[i];
	  low_[i]=order[i];
	  high_[i]=order[i];
	  lock_array[i]=0;
    }
  node_Barrier();

  pardo(i,0,n_edges,1)
  {
  	if(El[i].workspace!=1) {
		spin_lock(&lock_array[El[i].v1], MYTHREAD+1);
		if(low_[El[i].v1] > low[El[i].v2]) low_[El[i].v1] = low[El[i].v2];
		if(high_[El[i].v1]< high[El[i].v2]) high_[El[i].v1]=high[El[i].v2];
		spin_unlock(&lock_array[El[i].v1]);
	}
  }
  
  node_Barrier();

  pardo(i,0,n,1)
    {
      low[i]=low_[i];
      high[i]=high_[i];
#if DEBUG_EULER 
      printf("after the first step low[%d]=%d, high[%d]=%d\n",i,low[i],i,high[i]);
#endif
    }
  node_Barrier();
  D=Parent;
 
  pardo(i,0,n,1)
  {
  	j=i;
	while(D[j]!=D[D[j]]){
  		spin_lock(&lock_array[D[j]], MYTHREAD+1);
		if(low[D[j]]> low[j] || high[D[j]]<high[j]){
			if(low[D[j]]> low[j]) low[D[j]]=low[j];
			if(high[D[j]]<high[j]) high[D[j]]=high[j];
			spin_unlock(&lock_array[D[j]]);
			j = D[j];
		}
		else {
			spin_unlock(&lock_array[D[j]]);
			break;
		}
	}
  }

#if 0
  pardo(i,0,n,1)
  {
  	j=i;
	while(D[j]!=D[D[j]]){
  		spin_lock(&lock_array[D[j]], MYTHREAD+1);
		if(high[D[j]]< high[j]){
			high[D[j]]=high[j];
			spin_unlock(&lock_array[D[j]]);
			j = D[j];
		}
		else {
			spin_unlock(&lock_array[D[j]]);
			break;
		}
	}
  }
 #endif
 
  
  node_Barrier();
  node_free(lock_array,TH);
  node_free(low_,TH);
  node_free(high_,TH);

}

/*this is for the version where El has only one copy of eah edge*/
int Euler_get_lowhigh_s(E* El, int *Parent,int *order,int n,int n_edges, int root,int *low,int *high, THREADED)
{ 
	Euler_get_lowhigh_filter(El, Parent, order,n, n_edges,root,low,high,TH);
}


/*get the low of each node, jaja 238. Here is the version for bicc_filter. El is just the critical edges*/
int Euler_get_lowhigh_filter(E* El, int *Parent,int *order,int n,int n_edges, int root,int *low,int *high, THREADED)
{
  int i,j,u,v,pos,ret,done=0;
  int *finished,*low_,*high_,*lock_array;
  hrtime_t start,end;
  float interval;
  int * D, *leaf;
  
  
  low_=node_malloc(sizeof(int)*n,TH); 
  high_=node_malloc(sizeof(int)*n,TH); 
  lock_array = node_malloc(sizeof(int)*n,TH);
  
  pardo(i,0,n,1)
    {
      low[i]=order[i];
      high[i]=order[i];
	  low_[i]=order[i];
	  high_[i]=order[i];
	  lock_array[i]=0;
    }
  node_Barrier();

  pardo(i,0,n_edges,1)
  {
  	u = El[i].v1;
	v = El[i].v2;	
	spin_lock(&lock_array[u], MYTHREAD+1);
	if(low_[u] > low[v]) low_[u] = low[v];
	if(high_[u]< high_[v]) high_[u]=high[v];
	spin_unlock(&lock_array[u]);
	
	spin_lock(&lock_array[v],MYTHREAD+1);
	if(low_[v] > low[u]) low_[v] = low[u];
	if(high_[v]< high[u]) high_[v]=high[u];
	spin_unlock(&lock_array[v]);
  }
  
  node_Barrier();

  pardo(i,0,n,1)
    {
      low[i]=low_[i];
      high[i]=high_[i];
#if DEBUG_EULER 
      printf("after the first step low[%d]=%d, high[%d]=%d\n",i,low[i],i,high[i]);
#endif
    }
  node_Barrier();
  D=Parent;
 
  pardo(i,0,n,1)
  {
  	j=i;
	while(D[j]!=D[D[j]]){
  		spin_lock(&lock_array[D[j]], MYTHREAD+1);
		if(low[D[j]]> low[j] || high[D[j]]<high[j]){
			if(low[D[j]]> low[j]) low[D[j]]=low[j];
			if(high[D[j]]<high[j]) high[D[j]]=high[j];
			spin_unlock(&lock_array[D[j]]);
			j = D[j];
		}
		else {
			spin_unlock(&lock_array[D[j]]);
			break;
		}
	}
  }

  node_Barrier();
  node_free(lock_array,TH);
  node_free(low_,TH);
  node_free(high_,TH);

}

