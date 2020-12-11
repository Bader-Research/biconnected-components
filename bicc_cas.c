#include <sys/types.h>
#include "simple.h"
#include "graph.h"
#include "listrank.h"

#define NANO 1000000000

int initialize_graph_edgelist(V* graph, int n_vertices,E **pEL,int *pn_edge, THREADED);
ET*  pick_tree_edges(E* EL,int n_edges,THREADED);
void construct_Euler_path ( ET* treeEL, int n_edges,THREADED);
void Euler_root_tree(int *Parent,ET* treeEL,list_t *List,int n_edges,THREADED);
void label_twin_edges(E* EL, V* G, int n_edges, THREADED);
V* r_graph(int n,int m);
V* k_graph(int n, int k);
V* torus(int k);


/* still the tarjan viskin algorithm. Yet use lock free spanning tree
generation*/

int bicc_cas(E* El, V* G, int n_vertices, int n_edges,THREADED)
{
  int i,t,j,u,v,N,k;
  int *Low, *High,*Parent, *Size, * Preorder,logn,opt,l,s;
  hrtime_t start,end, s1,t1;
  double interval,total=0;
  int * D,tree_edges;
  int * graft_to,*Buff;
  long seed;
  ET * treeEL;
  E* El_tmp,*El_sel;
  list_t * List;
  int * K;

  tree_edges = 2*(n_vertices-1);
  logn = (int) log(n_vertices);
  
  K=(int*) node_malloc(sizeof(int)*THREADS,TH);
  Parent = (int *) node_malloc(sizeof(int)*n_vertices,TH);
  Low = (int *) node_malloc(sizeof(int)*n_vertices,TH);
  High = (int *) node_malloc(sizeof(int)*n_vertices,TH);
  Preorder = (int *) node_malloc(sizeof(int)*n_vertices,TH);
  Size = (int *) node_malloc(sizeof(int)*n_vertices,TH);
  List = (list_t *) node_malloc(sizeof(list_t)*2*n_vertices,TH);
  El_tmp = malloc(sizeof(E)*n_edges*2);
  
  pardo(i,0,n_vertices,1){
  	Parent[i]=i;
	Preorder[i]=0;
	Size[i]=0;
  }
    
  graft_to = node_malloc(sizeof(int)*n_vertices,TH);
  pardo(i,0,n_vertices,1) graft_to[i]=0;
  Buff = node_malloc(sizeof(int)*n_edges,TH);
  pardo(i,0,n_edges,1) Buff[i]=0;
  node_Barrier();
  
  start = gethrtime();
  s1 = gethrtime();
  spanning_tree_CRCW_cas(G,El,n_vertices,n_edges,TH);
  node_Barrier();
  t1 = gethrtime();
  interval = t1-s1;
  on_one printf("METRICS:Time used for spanning tree is %f s\n", interval/NANO);
  
  s1 = gethrtime();
  label_twin_edges(El,G,n_edges,TH);	
  node_Barrier();
  t1 = gethrtime();
  interval = t1-s1;
  on_one printf("METRICS:Time used for label twin edges is %f s\n", interval/NANO);


  s1 = gethrtime();
  treeEL= pick_tree_edges(El,n_edges,TH); 
  node_Barrier();
  t1 = gethrtime();
  interval = t1-s1;
  on_one printf("METRICS:Time used for pick_tree edges is %f s\n", interval/NANO);
  
  s1 = gethrtime();
  construct_Euler_path (treeEL,tree_edges,TH);
  node_Barrier();
  t1 = gethrtime();
  interval = t1-s1;
  on_one printf("METRICS:Time used for construct euler path is %f s\n", interval/NANO);
  
  Tree_to_list(treeEL,List,tree_edges,TH);
  node_Barrier();
  
  s1 = gethrtime();
  on_one printf("tree edges is %d \n", tree_edges);
  Euler_root_tree(Parent,treeEL,List,tree_edges,TH);
  t1 = gethrtime();
  interval = t1-s1;
  on_one printf("METRICS:Time used for rooting tree is %f s\n", interval/NANO);
  
  end = gethrtime();
  interval=end-start;
  on_one_thread
  	printf("METRICS:Time for spanning_tree+euler_tour is %f\n",interval/NANO);
	
  node_Barrier();
  s1 = gethrtime();
  Euler_preorder(Parent,Preorder,treeEL,List,tree_edges,TH);
  t1 = gethrtime();
  interval = t1-s1;
  on_one printf("METRICS:Time used for preorder tree is %f s\n", interval/NANO);
  
  
  node_Barrier();
  s1 = gethrtime();
  Euler_size_tree(Parent,Size,treeEL,List,tree_edges,TH);
  t1 = gethrtime();
  interval = t1-s1;
  on_one printf("METRICS:Time used for size tree is %f s\n", interval/NANO);
  on_one printf("size[0]=%d\n", Size[0]);
  
  s1 = gethrtime();
  Euler_get_lowhigh(El,Parent,Preorder,n_vertices,n_edges,0,Low,High, TH);
  t1 = gethrtime();
  interval = t1-s1;
  on_one printf("METRICS:Time used for Euler_get_lowhigh is %f s\n", interval/NANO);
  
  k=0;
  pardo(i,0,n_edges,1)
  {
  	if(El[i].workspace!=1 && Preorder[El[i].v2]<Preorder[El[i].v1]) {
		El_tmp[k++]=El[i];
		El_tmp[k].v2=El[i].v1;
		El_tmp[k].v1=Parent[El[i].v1];
  		k++;
	}
	if(El[i].workspace!=1 && Preorder[El[i].v2]+Size[El[i].v2] <=Preorder[El[i].v1]) {
		El_tmp[k].v1=El[i].v1;
		El_tmp[k].v2=Parent[El[i].v1];
		k++;
		El_tmp[k].v1=El[i].v2;
		El_tmp[k].v2=Parent[El[i].v2];
		k++;
	}
	if(El[i].workspace==1 && El[i].v2!=0 && El[i].v2==Parent[El[i].v1]  ) {
  		u = El[i].v1;
		v= El[i].v2;
		if(Low[u]<Preorder[v] || High[u]>=Preorder[v]+Size[v]){
			El_tmp[k++]=El[i];
			El_tmp[k].v1=v;
			El_tmp[k].v2=Parent[v];
			k++;
		}
	}
  }
  K[MYTHREAD]=k;
  
  printf(" k =%d\n", k);
  node_Barrier();
  
#if 0 
  El_sel=(E*)node_Bcast_ip((int*)El_tmp,TH);
  on_one{
  	for(i=1;i<THREADS;i++)
		K[i]+=K[i-1];
  }
  node_Barrier();
  
  if(MYTHREAD!=0)
  {   
	for(i=0;i<k;i++)
		El_sel[K[MYTHREAD-1]+i]=El_tmp[i];
  } 
  node_Barrier();
#endif
   
  connected_comp(El_tmp,n_vertices,k,TH);
  end = gethrtime();
  interval = end - start;
  interval = interval /NANO;
  on_one printf("METRICS: time used for biconn is %f s\n", interval);
  
#if CHECK
	pardo(i,0,n_vertices,1)
		while(Parent[i]!=Parent[Parent[i]]) Parent[i]=Parent[Parent[i]];
	node_Barrier();
	pardo(i,1,n_vertices,1)
		if(Parent[i]!=Parent[i-1]) printf("%d,%d,%d,wrong results\n",i, Parent[i],Parent[i-1]);
	node_Barrier();
	on_one printf("check done for euler tour rooting\n");
#endif

  node_free(graft_to,TH);
  node_free(Parent,TH);
  node_free(Size,TH);
  node_free(Preorder,TH);
  node_free(List,TH);
  node_Barrier();
  node_free(Parent,TH);
  node_free(Low,TH);
  node_free(High,TH);
  node_free(K,TH);
  free(El_tmp);
  
}

