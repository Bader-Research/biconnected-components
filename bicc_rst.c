#include <sys/types.h>
#include "simple.h"
#include "graph.h"
#include "listrank.h"

#define NANO 1000000000

V* r_graph(int n,int m);
V* k_graph(int n, int k);
V* torus(int k);
E* span_gw_euler(V* graph,int nVertices,THREADED);
void Euler_scan_size_tree(int *Parent,int * Size,E* Tour,int n_edges,THREADED);
void Euler_scan_preorder(int *Parent,int * Preorder,E* Tour,int n_edges,THREADED);

/* still the tarjan viskin algorithm. Use the rooted spanning tree to euler tour approach*/

int bicc_rst(E* El, V* G, int n_vertices, int n_edges,THREADED)
{
  int i,t,j,u,v,N,k;
  int *Low, *High,*Parent,*Size,* Preorder,* D,*K,tree_edges,logn,opt,l,s;
  hrtime_t start,end, s1,t1;
  double interval,d1,d2,d3,total=0;
  long seed;
  E* El_tmp, *final_tour;
  
  pardo(i,0,n_vertices,1)
  	G[i].v_attribute=-1;
  node_Barrier();
  	
  s1 = gethrtime();
  final_tour = span_gw_euler(G,n_vertices,TH);
  node_Barrier();
  t1 = gethrtime();
  interval=t1-s1;
  on_one printf("METRICS1: Time for span_gw_euler is %f\n",interval/NANO);
  
  
  tree_edges = 2*(n_vertices-1);
  logn = (int) log(n_vertices);
  
  K=(int*) node_malloc(sizeof(int)*THREADS,TH);
  Parent = (int *) node_malloc(sizeof(int)*n_vertices,TH);
  Low = (int *) node_malloc(sizeof(int)*n_vertices,TH);
  High = (int *) node_malloc(sizeof(int)*n_vertices,TH);
  Preorder = (int *) node_malloc(sizeof(int)*n_vertices,TH);
  Size = (int *) node_malloc(sizeof(int)*n_vertices,TH);
  El_tmp = malloc(sizeof(E)*n_edges);
  
  pardo(i,0,n_vertices,1){
  	Parent[i]=G[i].parent;
	Preorder[i]=0;
	Size[i]=0;
  }
	
  node_Barrier();
	
  
  pardo(i,0,n_edges,1)
  {
  	if(Parent[El[i].v1]==El[i].v2 || Parent[El[i].v2]==El[i].v1)
		El[i].workspace=1;
  }
  
  node_Barrier();
  
  start = gethrtime();
  s1 = gethrtime();
  Euler_scan_preorder(Parent,Preorder,final_tour,tree_edges,TH);
  t1 = gethrtime();
  interval = t1-s1;
  on_one printf("METRICS:Time used for preorder tree is %f s\n", interval/NANO);
  
  node_Barrier();
  s1 = gethrtime();
  Euler_scan_size_tree(Parent,Size,final_tour,tree_edges,TH);
  t1 = gethrtime();
  interval = t1-s1;
  on_one printf("METRICS:Time used for size tree is %f s\n", interval/NANO);
  interval = t1 - start;
  on_one printf("METRICS1:Time used for tree computation is %f s\n", interval/NANO);
  
  s1 = gethrtime();
  Euler_get_lowhigh(El,Parent,Preorder,n_vertices,n_edges,0,Low,High, TH);
  t1 = gethrtime();
  interval = t1-s1;
  on_one printf("METRICS1:Time used for Euler_get_lowhigh is %f s\n", interval/NANO);
  
  
  s1 = gethrtime();
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
  j = node_Reduce_i(k, SUM,TH);
  on_one printf("METRICS: number of comp edges is %d\n", j);
  t1 = gethrtime();
  interval = t1-s1;
  on_one printf("METRICS1:Time used for labeling comp edges is %f s\n", interval/NANO);
  node_Barrier();
  
  s1 = gethrtime(); 
  connected_comp(El_tmp,n_vertices,k,TH);
  t1 = gethrtime();
  interval = t1 - s1;
  interval = interval /NANO;
  on_one printf("METRICS1: time used for conn_comps is %f s\n", interval);

  node_free(Parent,TH);
  node_free(Size,TH);
  node_free(Preorder,TH);
  node_Barrier();
  node_free(Parent,TH);
  node_free(Low,TH);
  node_free(High,TH);
  node_free(K,TH);
  node_free(final_tour,TH);
  free(El_tmp);
  
}


/* get the size of each subtree. Only prefix sum is enough because of the consecutive layout of the tour*/
void Euler_scan_size_tree(int *Parent,int * Size,E* Tour,int n_edges,THREADED)
{
#define twin workspace

	int i,j;
	int * buff;
	
	buff = node_malloc(sizeof(int)*n_edges,TH);	
	pardo(i,0,n_edges,1){
		if(Tour[i].v1==Parent[Tour[i].v2])
			buff[i]=0;
		else buff[i]=1;
		
	}
	node_Barrier();	
		
	prefix_sum(buff,n_edges,TH);
	node_Barrier();		
	
	pardo(i,0,n_edges,1)
	{
		if(Parent[Tour[i].v1]==Tour[i].v2) {
			Size[Tour[i].v1]= buff[i]-buff[Tour[i].twin]; 
		}
	}

#if debug	
	on_one {
		printf("The values:\n");
		for(i=0;i<n_edges;i++)
			printf(" %d  ", buff[i]);
		printf("\n");
	}
#endif	

	node_Barrier();
	node_free(buff,TH);
	on_one Size[0]=n_edges;
#undef twin
}

/* get the size of each subtree*/
void Euler_scan_preorder(int *Parent,int * Preorder,E* Tour,int n_edges,THREADED)
{
	int i,j,*buff;
	
	buff = node_malloc(sizeof(int)*n_edges,TH);
	
	pardo(i,0,n_edges,1){
		if(Tour[i].v1==Parent[Tour[i].v2])
			buff[i]=1;
		else buff[i]=0;
	}
	node_Barrier();
		
	prefix_sum(buff, n_edges,TH);
	node_Barrier();

	node_Barrier();		
	pardo(i,0,n_edges,1)
	{
		if(Parent[Tour[i].v2]==Tour[i].v1) {
			Preorder[Tour[i].v2]= buff[i]; 
		}
	}

#if debug	
	on_one {
		printf("The values:\n");
		for(i=0;i<n_edges;i++)
			printf(" %d  ", buff[i]);
		printf("\n");
	}
#endif	
	node_Barrier();
	on_one Preorder[0] = 0;
	node_free(buff,TH);
}
