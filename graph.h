#ifndef _GRAPH_H_
#define _GRAPH_H_

#define DIRECTED 1
#define UNDIRECTED 0
#define DEBUG_BALANCE 1

typedef struct vertex{
   int self;
   int parent;    /*used only in spanning tree*/
   int v_attribute;
   int n_neighbors;
   int *my_neighbors;
   int *is_tree_edge; /* if<i,graph[i].my_neighbors[k]> is a tree edge in spanning tree*/
   int *pal_index;          /*for edge <i,graph[i].my_neighbors[k]> , what is the index for i in i,graph[i].my_neighbors[k]>'s neighborlist*/ 
   int edge_start; /* the position of the first edge going out from the vertex in edgelist*/
} V; /*tree structure for graph and spanning tree*/

typedef struct edge{
	int v1,v2;
	int workspace; 
	} E;

typedef struct tree_edge{
	int v1,v2;
	int workspace; 

	int next;
	int twin; /* next and twin are for euler tour technique of tree*/
	int pred; /*predecessor in the euler circuit*/
	int value,value_b; /*value_b and pred_b are for sync p-prefix*/
	int pred_b;
	} ET;
 
int initialize_graph(const char * file,V** graph, int *nVertices);

#endif /*_GRAPH_H_*/
