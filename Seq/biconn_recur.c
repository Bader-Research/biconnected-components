#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include "../graph.h"
struct s_ent{
	int a, b;};

V* G;
int *dfsnumber, *highwater,lastdfsnumber=0,top=-1,count=0,*cc;
struct s_ent * stack;
void push(int a, int b);
void pop(int *a, int *b);
void BiCC(int v, int p_v);
V* r_graph(int n,int m);
V* torus(int k);
V* k_graph(int n, int k);

int main(int argc, char** argv)
{
	int opt, n_vertices, n_edges, i, j,k;
	long seed;
	hrtime_t start,end;
	double interval;
	
	opt = atoi(argv[1]);
	seed = gethrtime();
	seed = -712456314 ;
	/*seed =660;*/
	seed=568;
	srand(seed);
	printf("METRICS: seed is %d \n", seed);
  	switch(opt){
  	case 0: n_vertices=atoi(argv[2]); n_edges=atoi(argv[3]);
				  G = r_graph(n_vertices,n_edges); break; 
  	case 1: k=atoi(argv[2]); n_vertices=k*k; 
				  G = torus(k); n_edges=n_vertices*4;break;
  	case 2: n_vertices=atoi(argv[2]); k = atoi(argv[3]);
				  G = k_graph(n_vertices,k); n_edges=n_vertices*k;break;
	default: printf("unknown graph type, exit\n"); exit(1);
  	}
	
	dfsnumber = malloc(sizeof(int)*n_vertices);
	highwater = malloc(sizeof(int)*n_vertices);
	cc = malloc(sizeof(int)*n_vertices);
	stack = malloc(sizeof(struct s_ent)*(n_edges+1000));

	for(i=0;i<n_vertices;i++){
		dfsnumber[i]=0;
		cc[i]=0;
	}

	start = gethrtime();
	BiCC(0,-1);
	end = gethrtime();
	interval = end - start;
	printf("METRICS: bicc_s done in %f s, found %d comps\n", interval/1000000000,count);
	
	free(dfsnumber);
	free(highwater);
	free(stack);
	free(cc);
	
	return(0);
}

void BiCC(int v, int u)
{
	int i,j,w,a,b;
	
	lastdfsnumber++;
	dfsnumber[v]=lastdfsnumber;
	highwater[v]=lastdfsnumber;
	
	for(i=0; i<G[v].n_neighbors;i++)
	{
		w=G[v].my_neighbors[i];
		if(dfsnumber[w]==0)
		{
			push(v,w);
			BiCC(w,v);
			if(highwater[w]<highwater[v]) highwater[v]=highwater[w];
			if(dfsnumber[v]<=highwater[w])
			{
				count++;
				pop(&a,&b);
				cc[a]=count;
				cc[b]=count;
				while(!(a==v && b==w)){
					pop(&a,&b);
					if(cc[a]==0) cc[a]=count;
					if(cc[b]==0) cc[b]=count;
				}
			} 
		}else if(dfsnumber[w]<dfsnumber[v] && w!=u)
					{
						push(v,w);
						if(dfsnumber[w]<highwater[v]) highwater[v]=dfsnumber[w];
					}
	}
		
}


void push(int a, int b)
{
   top++;
   stack[top].a=a;
   stack[top].b=b;
}

void pop(int * a, int *b)
{
	*a = stack[top].a;
	*b = stack[top].b;
	top--;
}
