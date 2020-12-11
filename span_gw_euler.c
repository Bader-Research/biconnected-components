
#include "simple.h"
#include "graph.h"
#include "stack.h"
#define NANO 1000000000
#define THRESHOLD 50

#define DEBUG 0
#define CHECK 0
	
typedef struct rep_node{
  int v;
  int n_children;
  int n_st_size; /* size of the piece of the small sub tree that it represents*/
  int n_gt_size; /* size of whole sub tree of the tree rooted at v, could be larger than n_st_size if has other rep as
  child*/
  int parent;
  int local_dfs_number;
  int global_start;
  int * children; 
  E *tour;
} rep_t;
typedef E* twin_t[2];

extern inline int spin_lock(int * lock, int id);
extern inline int spin_unlock(int * lock);


/* In this version, we use spinlocks to further help load-balancing. that is, work is not explicitly partitioned on the
processors. Instead, while there is work to do, they will each fetch a piece. Almost perfect load-balancing, as in
span_gw_euler.c.7. Here we try to eliminate one round of copy. For each rep vertex, previously its tour is first copied into its
field:tour, then all these tours are copied into the final tour. Here we do not copy to tour in to its field:tour, instead,
just mark where it begins in the auxilary tour, and have just one final round of copy. This should have better scalability
on SUN*/

/*And we intend this to be used for biconnected comps, so we need to do more
book keeping for it than in just rooted spanning tree to euler tour. The major
thing is we need to label the tree edges in edgelist that corresponds to tree
edges scattered in the adjacency list because later on biconn use edge lists*/

/* Also twin edge was not set in the ICPP paper, Here we need to set twins.
because computing size of tree need that. The basic idea is that when traverse
the tree block, keep track the twin information.*/

/*twin_relation is an array of size nVertices, with each element be a 1-d array
of 2 elements, each being a pointer to E*/

/* there are loopholes in selecting the rep vertices which may cause seg fault.
can be fixed easily though*/

E* span_gw_euler(V* graph, int nVertices,THREADED)
{
#define position workspace
#define S_POINTS (THREADS*THREADS*THREADS*2)
#define MYCOLOR (MYTHREAD+1)
  hrtime_t start,end;
  double interval,power;
  
  int i,j,k,p,l,v,u,s,n,r,ret=0,t=0;
  int root,walks,counter=0,visited=0,neighbor,n_neighbors,n_rep,add_n,n_local_tours=0;
  int first_time, work_to_steal,myroot,b,start_sr,start_sl;
  int top,count,bottom=-1,path_top=-1,path_bottom=-1;
  int ** stack_M, **top_M, **bottom_M,*stack,* finished,* path_stack,* count_M,*color;
  int * DFS_order, *Post_order,*Post_visited, *lock_array,*done;
  E * Tour,*final_Tour;
  rep_t * Rep_tree,rep;
  unsigned int seed=MYCOLOR;
  twin_t * twin_relation; 
  
  stack_M = node_malloc(THREADS*sizeof(int *),TH);
  top_M = node_malloc(THREADS*sizeof(int *),TH);
  bottom_M=node_malloc(THREADS*sizeof(int *),TH);
  stack_M[MYTHREAD]=malloc(nVertices*sizeof(int));
  stack=stack_M[MYTHREAD];
  top_M[MYTHREAD]=&top;
  bottom_M[MYTHREAD]=&bottom;
  count_M=node_malloc(THREADS*sizeof(int),TH);
  color=node_malloc(sizeof(int)*nVertices,TH); 
  twin_relation = node_malloc(sizeof(twin_t)*nVertices,TH);
    
  pardo(i,0,nVertices,1){
    color[i]=0;
  }
  
  start = gethrtime();  
  bottom=-1;
  top=-1;
  count=0;
  node_Barrier();

  seed=gethrtime()/(MYTHREAD+1)+MYTHREAD;
  if(ret==0) root=(rand_r(&seed)%nVertices);
  root = 0;
  
  /*lets select a point to start in the graph*/
  on_one_thread {
    color[root]=MYCOLOR;
    graph[root].parent=root;
    myroot=root;
    j=0;
    push(myroot,stack_M[j],top_M[j]);
    j=(j+1)%THREADS;
    i=0;
    for(i=0;i<S_POINTS;i++)
    {
		n_neighbors=graph[myroot].n_neighbors;
		r=rand_r(&seed);
		if(r%2==0){
	  	for(r=0;r<n_neighbors;r++)
	    {
	      visited++;
	      n=graph[myroot].my_neighbors[r];
		  		 
	      if(color[n]==0){
		  	graph[myroot].is_tree_edge[r]=1;
		  	graph[n].parent=myroot;
		  	color[n]=MYCOLOR;
		  	myroot=n;
		  	count++;
		  	push(myroot,stack_M[j],top_M[j]);
		  	/*printf("sub tree : %d \n",myroot);*/
		  	j=(j+1)%THREADS;
		  	break;
	    	} 
	 	}
	    if(r==n_neighbors){
	    	r=(rand_r(&seed)%n_neighbors);
	    	myroot=graph[myroot].my_neighbors[r];
	  	}
		}
	 	else {
	  		for(r=n_neighbors-1;r>=0;r--)
	    	{
	      		visited++;
	      		n=graph[myroot].my_neighbors[r];
	      		if(color[n]==0){
				graph[myroot].is_tree_edge[r]=1;
				graph[n].parent=myroot;
				color[n]=MYCOLOR;
				myroot=n;
				count++;
				push(myroot,stack_M[j],top_M[j]);
				/*printf("sub tree : %d \n",myroot);*/
				j=(j+1)%THREADS;
				break;
	     	 }     
	    	}
	  		if(r<0){
	    		r=(rand_r(&seed)%n_neighbors);
	    		myroot=graph[myroot].my_neighbors[r];
	  		}
		}		
      }
    end=gethrtime();
    interval=end-start;
  }
	       
  node_Barrier();
  on_one printf("Now walking...\n");
  
  start=gethrtime();
  first_time=1;
  work_to_steal=0;
  start_sr=MYTHREAD; /*when I am out of work, where do i start to search. r towards right, l towards left*/
  start_sl=MYTHREAD;
  
  while(first_time || work_to_steal)
  {
      while(!is_empty( stack, &top, &bottom))
	  {
	    n=pop(stack,&top,bottom);
	    if(n==-1){
	      printf("THREAD %d:stack overflow\n",MYTHREAD);
	      bottom=-1;
	      top=-1;
	      break;
	    }
	    visited+=graph[n].n_neighbors;
	    for(i=0;i<graph[n].n_neighbors;i++)
	    {
	      neighbor=graph[n].my_neighbors[i];
	      if(color[neighbor]==0) {/*found new frontier*/
			graph[n].is_tree_edge[i]=1;  /* in fact this records the number of children I have*/
		    color[neighbor]=1;
		    graph[neighbor].parent=n;
		    push(neighbor,stack,&top);
		    count++;
	      }
	    } 
	  }
      if(first_time) first_time=0;
      work_to_steal=0;
      if(MYTHREAD%2==0)
	  {
	    for(j=0;j<THREADS;j++)
	    {
	      i=(j+start_sr)%THREADS; /*start searching in circular from my neighbor*/
	      if(i==MYTHREAD) continue;

	      n=*(top_M[i]);
	      b=*(bottom_M[i]);	
	      if(n-b<THRESHOLD) continue;	      	       
	      if(count>nVertices/THREADS) /*I did my share*/
		     r=b+(n-b)/THREADS;
	      else r=b+max((n-b)/THREADS,min(nVertices/THREADS-count,(n-b)*3/4));
	      if(r<(*top_M[i])) {
		   (*bottom_M[i])=max(r-1,-1);
		   work_to_steal=1;		 
		   while((r--)>max(0,b))
		     push(stack_M[i][r],stack,&top);
		   break;		
	      }
	    }
	} else{
	  for(j=THREADS-1;j>0;j--)
	    {
	      i=(j+start_sl)%THREADS;
	      if(i==MYTHREAD) continue;
	      n=*(top_M[i]);
	      b=(*bottom_M[i]);
	      if(n-b<THRESHOLD) continue;	     
	      if(count>nVertices/THREADS) /*I did my share*/
		r=b+(n-b)/THREADS;
	      else r=b+max((n-b)/THREADS,min(nVertices/THREADS-count,(n-b)*3/4));
	      if(r<(*top_M[i])) {
		(*bottom_M[i])=max(r-1,-1);
		work_to_steal=1;
		while((r--)>max(b,0))
		  push(stack_M[i][r],stack,&top);
		break;
	      } 
	    }
	}
      start_sr=(start_sr+1)%THREADS;
      start_sl=(start_sl-1+THREADS)%THREADS;

    }
  count_M[MYTHREAD]=count;
  end=gethrtime();
  interval=end-start;
  node_Barrier();
  /*printf("Thread %d count is %d, visited is %d\n",MYTHREAD, count,visited);*/  
  on_one printf("METRICS1:Time used for span-gw: %f\n",interval/1000000000);

#if 0  
  on_one_thread{
    int max=0, min=nVertices;
    for(i=0;i<THREADS;i++)
      {
	if(count_M[i]>max) max=count_M[i];
	if(count_M[i]<min) min=count_M[i];
      }
    printf("METRICS===span_gw:The difference between counts is %d\n",max-min);
  }
#endif  
  node_Barrier();

  start = gethrtime();
  /* to clean up the mislabeled is_tree_edge because there might be race conditions*/ 
  pardo(i,0,nVertices,1)
  {
  	for(j=0;j<graph[i].n_neighbors;j++)
		if( graph[i].is_tree_edge[j]){
			if(graph[graph[i].my_neighbors[j]].parent!=i) 
				graph[i].is_tree_edge[j]=0;
		}
  
  }
  node_Barrier();
  end = gethrtime();
  interval = end-start;
  on_one printf("METRICS: time to set up correct children is %f \n", interval/NANO);
  
  
  /* Now lets build a euler-tour for the found spanning tree*/ 
  p = 40;
  s = p*THREADS;
  Rep_tree = node_malloc(sizeof(struct rep_node)*(s+1),TH);
  
  pardo(i,0,s+1,1)
  {
  	Rep_tree[i].n_children = 0;
	Rep_tree[i].children = malloc(sizeof(int)*(s+1));
	Rep_tree[i].global_start=0;
  	Rep_tree[i].n_st_size=0;
  	Rep_tree[i].parent=i;
  	Rep_tree[i].local_dfs_number=0;
  }
  
  path_stack = malloc(sizeof(int)*nVertices/THREADS);
  
  /* each processor chooses p non-leaf vertices, pathological case: not enough non-leaves
  vertices to choose. Wouldn't be likely if there are many vertices
  */ 	
  pardo(n,0,s,1) 
  {
    i= nVertices*n/(THREADS*p);	
    k=0;
    for(j=0;j<graph[i].n_neighbors;j++)
	  if(graph[i].is_tree_edge[j]==1) k++;
	while(k<2 && i<nVertices){
		k=0;
		i++;
		for(j=0;j<graph[i].n_neighbors;j++)
	  		if(graph[i].is_tree_edge[j]==1) k++;
	}	
	graph[i].v_attribute=n; /* label this as a rep vertex, and also at the same time show which rep it is*/
#if CHECK
	printf("Rep %d\n",i);
#endif
	Rep_tree[n].v=i;
  }
  
  node_Barrier();
  
  n_rep=s;
  on_one{
  	if(graph[root].v_attribute==-1) {
  		graph[root].v_attribute=s;
  		Rep_tree[s].v=root;
		n_rep++;
 	 }
  }
  n_rep = node_Bcast_i(n_rep,TH);
  node_Barrier();
/*  on_one printf("n_rep is %d\n", n_rep);*/
  lock_array = node_malloc(sizeof(int)*n_rep,TH);
  done = node_malloc(sizeof(int)*n_rep,TH);
  pardo(i,0,n_rep,1){
  	lock_array[i]=0;
	done[i]=0;
  }
  node_Barrier();
  
  start = gethrtime(); 
  /* do local depth-first search */
  visited=0;
  t = 0;
  Tour = malloc(sizeof(E)* 2*nVertices);
  for(n=0; n<n_rep; n++){
  	spin_lock(&(lock_array[n]),MYTHREAD+1);
	if(done[n]) {
		spin_unlock(&(lock_array[n]));
		continue;
	}else{
		done[n]=1;
		spin_unlock(&lock_array[n]);
	}
	
  	v=Rep_tree[n].v;
	r=v;
	/*printf("Thread %d: my sub tree starts at %d\n", MYTHREAD, v);*/

	top=-1;
	bottom=-1;
	push(v,stack,&top);
	l=0; /* for local length */
	Rep_tree[n].n_st_size=0;
	while(!is_empty( stack, &top, &bottom))
	{
		v=pop(stack,&top,bottom);
		visited++;
		if(v==-1) printf("stack overflow\n");
		/*printf(" pop out %d \n",v);*/
		if(l!=0){			
			Tour[t+l-1].v1=graph[v].parent;
			Tour[t+l-1].v2=v;
			Tour[t+l-1].position=l-1;
			twin_relation[v][0]=&(Tour[t+l-1]);
		}
		Rep_tree[n].n_st_size++;
			
		if(graph[v].v_attribute!=-1 && v!=r){ /* found a rep*/
			Rep_tree[graph[v].v_attribute].local_dfs_number= l;
			Rep_tree[graph[v].v_attribute].parent=n;
			Rep_tree[n].children[Rep_tree[n].n_children++]= graph[v].v_attribute;/* setup the rep tree relation*/
		} else{
			for(i=0;i<graph[v].n_neighbors;i++)
			{
				if(!graph[v].is_tree_edge[i]) continue;
				neighbor=graph[v].my_neighbors[i];
	      		if(neighbor!=graph[v].parent) {
		    		push(neighbor,stack,&top);
	      		}
			}
		}
		l++;
		if(is_empty(stack,&top,&bottom) || graph[stack[top]].parent!=v){ /*back track*/
			if(l!=0){		
				Tour[t+l-1].v1=v;
				Tour[t+l-1].v2=graph[v].parent;
				Tour[t+l-1].position=l-1;	
				twin_relation[v][1]=&(Tour[t+l-1]);		
			}
			l++;
			while(!is_empty(path_stack,&path_top,&path_bottom) &&
			(is_empty(stack,&top,&bottom) ||
			path_stack[path_top]!=graph[stack[top]].parent)) {
				u=pop(path_stack,&path_top,path_bottom);
				if(l!=0 && !is_empty(path_stack,&path_top,&path_bottom)){				
					Tour[t+l-1].v1=u;
					Tour[t+l-1].v2=path_stack[path_top];
					Tour[t+l-1].position=l-1;	
					twin_relation[u][1]=&(Tour[t+l-1]);			
					
				}
				l++;
			}
		} else push(v,path_stack,&path_top); /* keep track of the current path */
	}

#if 0
	Rep_tree[n].tour=malloc(sizeof(E)*2*(Rep_tree[n].n_st_size-1));
	for(i=0;i<2*(Rep_tree[n].n_st_size-1);i++){
		Rep_tree[n].tour[i]=Tour[i];
		/*printf(" THREAD %d: <%d,%d > \n", MYTHREAD,Tour[i].v1,Tour[i].v2);*/
    }
#endif

	Rep_tree[n].tour = &(Tour[t]);		
	t+=l;

#if CHECK	
	n_local_tours+=2*(Rep_tree[n].n_st_size-1);
	for(i=0;i<2*(Rep_tree[n].n_st_size-1)-1;i++)
		if(Rep_tree[n].tour[i].v2!=Rep_tree[n].tour[i+1].v1) 
			printf("ERROR in LOCAL DFS:[%d].v1=%d,[%d].v2=%d \n", i,Rep_tree[n].tour[i].v2,i+1,Rep_tree[n].tour[i+1].v1);
#endif
		
  }
  end = gethrtime();
  interval = end -start;
  
 /* printf("METRICS: Thread %d: Time used to do local DFS is %f s, visited is
 %d\n", MYTHREAD,interval/NANO,visited);*/
  node_Barrier();  
  
  interval = node_Reduce_i(interval,MAX,TH);
  on_one printf("METRICS: time used on local-dfs search is %f s\n", interval/1000000000);
  
  /* process the Rep_tree structure. First we do a depth-first search of the Rep_tree, and 
     list the rep vertices in that order. For each Rep vertex we want to find out where is its starting location for the 
	 euler tour of the subtree.
	 If in the list my pred is my parent, then my starting position is pred's start+ local_length
	 Else my start is my pred's start+ pred's size(to be exact, 2(n_st_size-1))+local_length */
  
  on_one{

	start = gethrtime();
	/*printf("\n post-order:\n");*/
	
	/* find the post ordering of the rep vertices, reuse the piece of memory of DFS_order.
	We use the post ordering of the rep_tree to find n_gt_size. Probably we don't need a very efficient impl for this*/
	Post_order = malloc(sizeof(int)*n_rep); /* shows a sequence of vertices in post ordering*/
	Post_visited=malloc(sizeof(int)*n_rep);
	for(i=0;i<n_rep;i++)
		Post_visited[i]=0;
		
	n = graph[root].v_attribute; /*root*/
	top = -1;
	bottom =-1;
	i=0;
	j=0;
	push(n,stack,&top);
	while(!is_empty( stack, &top, &bottom))
	{
		v = stack[top];
		k=0;
		for(i=Rep_tree[v].n_children-1;i>=0;i--)
			if(!Post_visited[Rep_tree[v].children[i]]){
				push(Rep_tree[v].children[i],stack,&top);
				k++;
		}
		if(k==0) {
			v=pop(stack,&top,bottom);
#if DEBUG
			printf("%d - ", Rep_tree[v].v);
#endif
			Post_visited[v]=1;
			Post_order[j++]=v;
			
		}	
		
	}
	
#if DEBUG 
	printf("\n");
#endif

/* set n_gt_size */	
	for(i=0;i<n_rep;i++)
	{
		Rep_tree[Post_order[i]].n_gt_size=Rep_tree[Post_order[i]].n_st_size;
		for(k=0;k<Rep_tree[Post_order[i]].n_children;k++)
			Rep_tree[Post_order[i]].n_gt_size+=(Rep_tree[Rep_tree[Post_order[i]].children[k]].n_gt_size-1);
#if DEBUG
		v = Rep_tree[Post_order[i]].v;
		printf("[%d].n_gt_size = %d,[%d].n_st_size =%d",v,Rep_tree[Post_order[i]].n_gt_size,v,Rep_tree[Post_order[i]].n_st_size);		
#endif
	}
  
    /* find the DFS ordering of the rep vertices, e.g., DFS_order[0] points to the
	root */
#if DEBUG
	printf("\nDFS_order:\n");
#endif

  	DFS_order = malloc(sizeof(struct rep_node)*n_rep);
    n = graph[root].v_attribute; /*root*/
	top = -1;
	bottom =-1;
	i=0;
	push(n,stack,&top);
	while(!is_empty( stack, &top, &bottom))
	{
		v = pop(stack,&top,bottom);
		DFS_order[i++]=v;
#if DEBUG
		printf("%d - ", Rep_tree[v].v);
#endif
		
		for(k=Rep_tree[v].n_children-1;k>=0;k--)
			if(Rep_tree[v].children[k]!=v) push(Rep_tree[v].children[k],stack,&top);
	}
	
#if DEBUG
	printf("\n");
#endif
	
	/* compute the starting position for each rep_vertex */
	
	Rep_tree[DFS_order[0]].global_start =0;
	for(i=1;i< n_rep;i++){
		if( Rep_tree[DFS_order[i]].parent==DFS_order[i-1] ) {
			Rep_tree[DFS_order[i]].global_start +=
			(Rep_tree[DFS_order[i-1]].global_start+Rep_tree[DFS_order[i]].local_dfs_number);
		} else {
			v = Rep_tree[DFS_order[i]].parent;
			Rep_tree[DFS_order[i]].global_start=Rep_tree[v].global_start;
			Rep_tree[DFS_order[i]].global_start+=Rep_tree[DFS_order[i]].local_dfs_number;
			for(j=0;j<Rep_tree[v].n_children;j++) {		
				if(Rep_tree[v].children[j]==DFS_order[i]) break;
				Rep_tree[DFS_order[i]].global_start+=2*(Rep_tree[Rep_tree[v].children[j]].n_gt_size-1);			
			}
		}
	}
	
#if DEBUG
	for(i=0;i<n_rep;i++)
	{
		v = Rep_tree[DFS_order[i]].v;
		printf("[%d].global_start=%d,[%d].local_dfs_number=%d\n",v,Rep_tree[DFS_order[i]].global_start,v,Rep_tree[DFS_order[i]].local_dfs_number);		
	}
#endif
	free(DFS_order);
	free(Post_visited);
	free(Post_order);
 	
	end = gethrtime();
	interval = end - start;
  }
  node_Barrier();
  
  /* Now copy the edges to the appropriate place */
  node_Barrier();
  
  start = gethrtime();
  final_Tour = node_malloc(sizeof(E)*2*(nVertices-1),TH);

 #if CHECK
   pardo(i,0,2*(nVertices-1),1)
   {
   	Tour[i].position=-1;
   }
 #endif
 
  pardo(n,0,n_rep,1)
  {
    
  	add_n = Rep_tree[n].global_start;
  	for(i=0;i<(Rep_tree[n].n_st_size-1)*2;i++)
	{
		Rep_tree[n].tour[i].position+=add_n;
	}
	
	for(k=0;k<Rep_tree[n].n_children;k++)
	{
		add_n = 2*(Rep_tree[Rep_tree[n].children[k]].n_gt_size-1);
		for(i=Rep_tree[Rep_tree[n].children[k]].local_dfs_number;i<(Rep_tree[n].n_st_size-1)*2;i++)
			Rep_tree[n].tour[i].position+=add_n;
	}
	
	for(i=0;i<(Rep_tree[n].n_st_size-1)*2;i++)
	  final_Tour[Rep_tree[n].tour[i].position]= Rep_tree[n].tour[i];
	
  }
  node_Barrier();
  
  /*set the twin relation: in final_Tour, we use the 'position' field to note the twin*/
  /* and starts at 1 because no tour for root 0*/
  pardo(i,1,nVertices,1)
  {
	final_Tour[twin_relation[i][0]->position].position=twin_relation[i][1]->position;
	final_Tour[twin_relation[i][1]->position].position=twin_relation[i][0]->position;
  }
  node_Barrier();
  
  end = gethrtime();
  interval = end - start;
  on_one printf("METRICS: Time used on copying to appropriate location is %f s\n", interval/NANO);
  
#if CHECK  
  n_local_tours=node_Reduce_i(n_local_tours,SUM,TH);
  on_one printf(" Total number of elements in the tour is %d \n", n_local_tours);
  pardo(i,0,2*(nVertices-1)-1,1)
  {
  	if(final_Tour[i].v2!=final_Tour[i+1].v1) printf("error in the Tour: %d:<%d,%d,%d>->%d<%d,%d,%d>\n",
	i,final_Tour[i].v1,final_Tour[i].v2, final_Tour[i].position,i+1,final_Tour[i+1].v1,final_Tour[i+1].v2,final_Tour[i+1].position);
  }
#endif


#if DEBUG
on_one{
	for(i=0;i<2*(nVertices-1);i++)
		printf(" <%d,%d> \n", final_Tour[i].v1,final_Tour[i].v2);
}
node_Barrier();
#endif
		
  node_free(twin_relation,TH);  
  node_free(count_M,TH);
  node_free(color, TH);
  node_free(stack_M,TH);
  node_free(top_M,TH);
  on_one{
  	for(i=0;i<n_rep;i++){
		free(Rep_tree[i].children);
	}
  }
  node_Barrier();
  node_free(Rep_tree,TH);
  node_free(lock_array,TH);
  node_free(done,TH);
  free(Tour);
  free(stack);
  free(path_stack);
  
  return(final_Tour);
#undef position
}
