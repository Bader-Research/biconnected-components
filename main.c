#include <sys/types.h>
#include "simple.h"
#include "graph.h"
#include "listrank.h"

#define NANO 1000000000
#define CHECK 0

V* G;
E* El;
int n_edges;
int initialize_graph_edgelist(V* graph, int n_vertices,E **pEL,int *pn_edge, THREADED);
V* r_graph(int n,int m);
V* k_graph(int n, int k);
V* torus(int k);


void *SIMPLE_main(THREADED)
{
  int i,t,j,u,v, n_vertices,N,k,opt;
  hrtime_t start,end, s1,t1;
  double interval,total=0;
  char * input_file;
  long seed;
  

# if 0 
  /*initialize graph from input file */  
  input_file = THARGV[0];  
  on_one_thread{
	t = initialize_graph(input_file,&G,&n_vertices);
	if(t!=0) exit(0);
  }
  node_Barrier();
#endif
  
#if 1 
  on_one{
  	opt = atoi(THARGV[0]);
	seed = gethrtime();
	seed = -712456314 ;
	/*seed =660;*/
	seed=568;
	srand(seed);
	printf("METRICS: seed is %d \n", seed);
  	switch(opt){
  	case 0: n_vertices=atoi(THARGV[1]); n_edges=atoi(THARGV[2]);
				  G=r_graph(n_vertices,n_edges); break; 
  	case 1: k=atoi(THARGV[1]); n_vertices=k*k; 
				  G = torus(k); break;
  	case 2: n_vertices=atoi(THARGV[1]); k = atoi(THARGV[2]);
				  G = k_graph(n_vertices,k); break;
	default: printf("unknown graph type, exit\n"); exit(1);
  	}
	printf("n_edges=%d\n",n_edges);
  }     
  node_Barrier();
#endif
  
  n_vertices=node_Bcast_i(n_vertices,TH);
 
  start = gethrtime();
  initialize_graph_edgelist(G, n_vertices,&El,&n_edges, TH);
  node_Barrier();
  end = gethrtime();
  interval=end-start;
  on_one printf("Time for initialization(graph edge list) is %f, n_vertices=%d,n_edges = %d\n",interval/NANO,n_vertices,n_edges);
 
  node_Barrier();

  start = gethrtime();
  bicc_tv(El,G, n_vertices, n_edges,TH);
  end = gethrtime();
  interval=end-start;
  node_Barrier();
  on_one printf("METRICS:bicc_tv uses %f s\n", interval/NANO);
  
  pardo(i,0,n_edges,1)
  	El[i].workspace=0;

  node_Barrier();
  start = gethrtime();
  bicc_rst(El,G, n_vertices, n_edges,TH);
  end = gethrtime();
  interval=end-start;
  node_Barrier();
  on_one printf("METRICS:bicc_rst uses %f s\n", interval/NANO);

  
  pardo(i,0,n_edges,1)
  	El[i].workspace=0;

  node_Barrier();
  start = gethrtime();
  bicc_filter(El,G, n_vertices, n_edges,TH);
  end = gethrtime();
  interval=end-start;
  node_Barrier();
  on_one printf("METRICS:bicc_filter uses %f s\n", interval/NANO);
  
  on_one_thread{
	delete_graph(G,n_vertices);
	printf("delete_graph done\n");
	if(El) free(El);
	printf("delete EL done\n");
  } 
 
  SIMPLE_done(TH);
}


