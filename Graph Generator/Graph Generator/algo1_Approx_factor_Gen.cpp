#include<algorithm>
#include<fstream>
#include<iostream>
#include<iterator>
#include<limits.h>
#include<list>
#include<map>
#include<math.h>
#include<stdio.h>
#include<string>
#include<stdlib.h>
#include<time.h>
#include<vector>
#include<set>
#include<stack>
#include<unordered_map>

using namespace std;
 
// Number of vertices in the graph
unsigned __int64 V = 5 ;
unsigned __int64 MAX_NODES = 200;


struct sort_pair_by_first_element
{
    bool operator()(const pair<unsigned __int64,unsigned __int64> &left, const pair<unsigned __int64,unsigned __int64> &right) 
	{
		return left.first < right.first;
    }
};
struct sort_pair_by_second_element
{
    bool operator()(const pair<unsigned __int64,unsigned __int64> &left, const pair<unsigned __int64,unsigned __int64> &right) 
	{
        return left.second < right.second;
    }
};

//struct sort_edges_by_first_vertex
//{
//    bool operator()(const pair<pair<unsigned __int64,unsigned __int64>,unsigned __int64> &left, const pair<pair<unsigned __int64,unsigned __int64>,unsigned __int64> &right) 
//	{
//		return left.first.first < right.first.first;
//    }
//};
//struct sort_edges_by_second_vertex
//{
//    bool operator()(const pair<pair<unsigned __int64,unsigned __int64>,unsigned __int64> &left, const pair<pair<unsigned __int64,unsigned __int64>,unsigned __int64> &right) 
//	{
//		return left.first.second < right.first.second;
//    }
//};

struct sort_edges_by_like_set
{
    bool operator()(const pair<pair<unsigned __int64,unsigned __int64>,unsigned __int64> &left, const pair<pair<unsigned __int64,unsigned __int64>,unsigned __int64> &right) 
	{
		if(left.first.first == right.first.first)
			return left.first.second < right.first.second;
		else
			return left.first.first < right.first.first;
    }
};

// A utility function to find the vertex with minimum key value, from
// the set of vertices not yet included in MST
unsigned __int64 minKey(unsigned __int64 key[], bool mstSet[])
{
	// Initialize min value
	unsigned __int64 min = INT_MAX, min_index;
    
	for (unsigned __int64 v = 0; v < V; v++)
		if (mstSet[v] == false && key[v] < min)
			min = key[v], min_index = v;

	return min_index;
}

//int printMST(int parent[], int n, int **graph)
////int printMST(int parent[], int n, int graph[V][V])
//{
//	printf("Edge   Weight\n");
//		for (int i = 1; i < V; i++)
//			printf("%d - %d    %d \n", parent[i], i, graph[i][parent[i]]);
//	return 0;
//}
unsigned __int64 printMST(unsigned __int64 parent[], unsigned __int64 n, vector<vector<pair<unsigned __int64, unsigned __int64> > > &graph)
//unsigned __int64 printMST(unsigned __int64 parent[], unsigned __int64 n, vector<map<unsigned __int64, unsigned __int64> > &graph)
//unsigned __int64 printMST(unsigned __int64 parent[], unsigned __int64 n, vector<vector<unsigned __int64> > &graph)
//int printMST(int parent[], int n, int graph[V][V])
{
	//printf("Edge   Weight\n");
	//	for (unsigned __int64 i = 1; i < V; i++)
	//		printf("%d - %d    %d \n", parent[i], i, graph[i][parent[i]]);

	printf("Edge   Weight\n");
		for (unsigned __int64 i = 1; i < V; i++)
		{	
			__int64 find_v_in_u = 0;
			for(find_v_in_u = 0; find_v_in_u< graph[i].size(); find_v_in_u++)
				if(graph[i][find_v_in_u].first == parent[i])
					break;
			printf("%d - %d    %d \n", parent[i], i, graph[i][find_v_in_u].second);
		}
	return 0;
}
//int findWeight(int parent[], int n, int **graph)
////int printMST(int parent[], int n, int graph[V][V])
//{
//	int weight = 0;
//	//printf("Edge   Weight\n");
//		for (int i = 1; i < V; i++)
//			weight += graph[i][parent[i]] ;
//	return weight ;
//}

unsigned __int64 findWeight(unsigned __int64 parent[], unsigned __int64 n, vector<vector<pair<unsigned __int64, unsigned __int64> > > &graph)
//unsigned __int64 findWeight(unsigned __int64 parent[], unsigned __int64 n, vector<map<unsigned __int64, unsigned __int64> > &graph)
//unsigned __int64 findWeight(unsigned __int64 parent[], unsigned __int64 n, vector<vector<unsigned __int64> > &graph)
//int printMST(int parent[], int n, int graph[V][V])
{
	//unsigned __int64 weight = 0;
	////printf("Edge   Weight\n");
	//	for (unsigned __int64 i = 1; i < V; i++)
	//		weight += graph[i][parent[i]] ;
	//return weight ;

	unsigned __int64 weight = 0;
	//printf("Edge   Weight\n");
		for (unsigned __int64 i = 1; i < V; i++)
		{			
			__int64 find_v_in_u = 0;
			for(find_v_in_u = 0; find_v_in_u< graph[i].size(); find_v_in_u++)
				if(graph[i][find_v_in_u].first == parent[i])
					break;
			weight += graph[i][find_v_in_u].second ;
		}
			
	return weight ;
}
 
// Function to construct and print MST for a graph represented using adjacency
// matrix representation
//template<int R, int C>
//int primMST(int **graph)
////void primMST(int **graph)
////void primMST(int graph[V][V])
//{
//	int *parent = new int[V]; // Array to store constructed MST
//	int *key = new int[V];   // Key values used to pick minimum weight edge in cut
//	bool *mstSet = new bool[V];  // To represent set of vertices not yet included in MST
//
//	// Initialize all keys as INFINITE
//	for (int i = 0; i < V; i++)
//		key[i] = INT_MAX, mstSet[i] = false;
// 
//	// Always include first 1st vertex in MST.
//	key[0] = 0;     // Make key 0 so that this vertex is picked as first vertex
//	parent[0] = -1; // First node is always root of MST
//
//	// The MST will have V vertices
//	for (int count = 0; count < V-1; count++)
//	{
//		// Pick thd minimum key vertex from the set of vertices
//		// not yet included in MST
//		int u = minKey(key, mstSet);
//
//		// Add the picked vertex to the MST Set
//		mstSet[u] = true;
//
//		// Update key value and parent index of the adjacent vertices of
//		// the picked vertex. Consider only those vertices which are not yet
//		// included in MST
//		for (int v = 0; v < V; v++)
//
//		// graph[u][v] is non zero only for adjacent vertices of m
//		// mstSet[v] is false for vertices not yet included in MST
//		// Update the key only if graph[u][v] is smaller than key[v]
//			if (graph[u][v] && mstSet[v] == false && graph[u][v] <  key[v])
//				parent[v]  = u, key[v] = graph[u][v];
//	}
// 
//	// print the constructed MST
//	//printMST(parent, V, graph);
//	
//	//std::cout<<findWeight(parent, V, graph);
//	//std::cout<<endl;
//
//	return findWeight(parent, V, graph);
//}
unsigned __int64 primMST(vector<vector<pair<unsigned __int64, unsigned __int64> > > &graph)
//unsigned __int64 primMST(vector<map<unsigned __int64, unsigned __int64> >  &graph)
//unsigned __int64 primMST(vector<vector<unsigned __int64> > &graph)
//void primMST(int **graph)
//void primMST(int graph[V][V])
{
	unsigned __int64 *parent = new unsigned __int64[V]; // Array to store constructed MST
	unsigned __int64 *key = new unsigned __int64[V];   // Key values used to pick minimum weight edge in cut
	bool *mstSet = new bool[V];  // To represent set of vertices not yet included in MST

	// Initialize all keys as INFINITE
	for (unsigned __int64 i = 0; i < V; i++)
		key[i] = INT_MAX, mstSet[i] = false;
 
	// Always include first 1st vertex in MST.
	key[0] = 0;     // Make key 0 so that this vertex is picked as first vertex
	parent[0] = -1; // First node is always root of MST

	// The MST will have V vertices
	for (unsigned __int64 count = 0; count < V-1; count++)
	{
		// Pick thd minimum key vertex from the set of vertices
		// not yet included in MST
		unsigned __int64 u = minKey(key, mstSet);

		// Add the picked vertex to the MST Set
		mstSet[u] = true;

		// Update key value and parent index of the adjacent vertices of
		// the picked vertex. Consider only those vertices which are not yet
		// included in MST

		//for (int v = 0; v < V; v++)
		//{
		//// graph[u][v] is non zero only for adjacent vertices of m
		//// mstSet[v] is false for vertices not yet included in MST
		//// Update the key only if graph[u][v] is smaller than key[v]
		//	if (graph[u][v] && mstSet[v] == false && graph[u][v] <  key[v])
		//		parent[v]  = u, key[v] = graph[u][v];
		//}

		for (unsigned __int64 v = 0; v < V; v++)
		{
		// graph[u][v] is non zero only for adjacent vertices of m
		// mstSet[v] is false for vertices not yet included in MST
		// Update the key only if graph[u][v] is smaller than key[v]
			__int64 find_v_in_u = 0;
			for(find_v_in_u = 0; find_v_in_u< graph[u].size(); find_v_in_u++)
				if(graph[u][find_v_in_u].first == v)
				{
					if (graph[u][find_v_in_u].second && mstSet[v] == false && graph[u][find_v_in_u].second <  key[v])
						parent[v]  = u, key[v] = graph[u][find_v_in_u].second;
					break;
				}
		}
	}
 
	// print the constructed MST
	//printMST(parent, V, graph);
	
	//std::cout<<findWeight(parent, V, graph);
	//std::cout<<endl;

	return findWeight(parent, V, graph);
}
 
//void generate_complete_graph(unsigned __int64 &v_count, vector<unsigned __int64> &vertices, vector<pair<unsigned __int64,unsigned __int64> > &all_edges, unsigned __int64 &max_weight)
//{
//	//clearing previous data
//	vertices.clear();
//	all_edges.clear();
//
//	for(unsigned __int64 i=0; i<v_count; i++)
//	//for(unsigned __int64 i=1; i<=v_count; i++)
//		vertices.push_back(i);
//	for(unsigned __int64 i=0; i<v_count; i++)
//		for(unsigned __int64 j=i+1; j<v_count; j++)
//			all_edges.push_back(pair<unsigned __int64,unsigned __int64> (i,j));
//}

void generate_random_edges_making_one_connected_component(unsigned __int64 &v_count, vector<pair<pair<unsigned __int64,unsigned __int64>, unsigned __int64 > > &edges, unsigned __int64 &total_edges, unsigned __int64 &max_weight)
{
	//clearing previous data
	edges.clear();

	set<pair<unsigned __int64,unsigned __int64> > random_edges;
	
	//making sure that it is 1 connected component
	for(int i=0; i<v_count; i++)
	{
		random_edges.insert(pair<unsigned __int64,unsigned __int64>(i,(i+1)%v_count)) ;
		random_edges.insert(pair<unsigned __int64,unsigned __int64>((i+1)%v_count,i)) ;
	}

	while(random_edges.size() <= total_edges)
	{
		__int64 u, v ;
		u = rand()%v_count;
		do
		{
			v = rand()%v_count;
		}while(u==v);
		random_edges.insert(pair<unsigned __int64,unsigned __int64>(u,v));
		random_edges.insert(pair<unsigned __int64,unsigned __int64>(v,u));
	}
	//for(set<pair<unsigned __int64,unsigned __int64> >::iterator it = random_edges.begin(); it != random_edges.end(); it++)
	//{
	//	edges.push_back(pair<pair<unsigned __int64,unsigned __int64>, unsigned __int64 > (*it,rand()%max_weight + 1));
	//}

	vector<pair<unsigned __int64,unsigned __int64> > random_edges_vector(random_edges.begin(),random_edges.end());
	do
	{
		__int64 u = random_edges_vector[0].first, v = random_edges_vector[0].second ;		
		__int64 index = find(random_edges_vector.begin(), random_edges_vector.end(), pair<unsigned __int64,unsigned __int64>(v,u)) - random_edges_vector.begin() ;
		__int64 wt = rand()%max_weight + 1 ;
		pair<unsigned __int64,unsigned __int64> pair_uv(u,v), pair_vu(v,u);
		edges.push_back(pair<pair<unsigned __int64,unsigned __int64>, unsigned __int64 > (pair_uv,wt));
		edges.push_back(pair<pair<unsigned __int64,unsigned __int64>, unsigned __int64 > (pair_vu,wt));
		random_edges_vector.erase(random_edges_vector.begin() + index);
		random_edges_vector.erase(random_edges_vector.begin() + 0);
	}while(random_edges_vector.size() != 0);
	stable_sort(edges.begin(),edges.end(),sort_edges_by_like_set());


	//for(unsigned __int64 i=1; i<=total_edges; i++)
	//{
	//	unsigned __int64 random_edge = rand()%all_edges.size();
	//	unsigned __int64 wt = rand()%max_weight + 1 ;

	//	pair<unsigned __int64,unsigned __int64> p1(all_edges[random_edge].first , all_edges[random_edge].second);
	//	//fout<<all_edges[random_edge].first<<" "<<all_edges[random_edge].second<<" "<<wt<<endl;
	//	edges.push_back(pair<pair<unsigned __int64,unsigned __int64>,unsigned __int64>(p1,wt));

	//	pair<unsigned __int64,unsigned __int64> p2(all_edges[random_edge].second , all_edges[random_edge].first);
	//	//fout<<all_edges[random_edge].second<<" "<<all_edges[random_edge].first<<" "<<wt<<endl;
	//	edges.push_back(pair<pair<unsigned __int64,unsigned __int64>,unsigned __int64>(p2,wt));

	//	//all_edges.erase(all_edges.begin() + random_edge);
	//}
	//stable_sort(edges.begin(),edges.end(),sort_edges_by_like_set());
}

__int64 Prims_algo_MST(vector<pair<pair<unsigned __int64,unsigned __int64>, unsigned __int64 > > &edges, __int64 v_count, __int64 total_edges)
{
		V = v_count ;
	
	unsigned __int64 E;
	E = total_edges ;

	//int **graph = new int*[V] ; 
	//for(int i=0; i<V; i++)
	//{
	//	int *g = new int[V];
	//	for(int j=0; j<V; j++)
	//	{
	//		//std::cin>>g[j];
	//		g[j] = 0;
	//	}
	//	graph[i] = g;
	//}
	
	////adj-matrix
	//vector<vector<unsigned __int64> > graph;
	//vector<unsigned __int64> vet;
	//for(unsigned __int64 i=0; i<V; i++)
	//{
	//	vet.push_back(0);
	//}
	//for(unsigned __int64 i=0; i<V; i++)
	//{
	//	graph.push_back(vet);
	//}
	////

	////adj-List
	//vector<map<unsigned __int64, unsigned __int64> > graph;
	//for(unsigned __int64 i=0; i<V; i++)
	//{	
	//	map<unsigned __int64, unsigned __int64> vet;
	//	vet[(i+1)%V] = 0 ; 
	//	graph.push_back(vet);
	//}
	////

	//for(unsigned __int64 i=0; i<edges.size(); i++)
	//{
	//	unsigned __int64 aa,ab,ac;
	//	aa=edges[i].first.first;	ab=edges[i].first.second;	ac=edges[i].second;
	//	graph[edges[i].first.first][edges[i].first.second] = edges[i].second ;
	//	graph[edges[i].first.second][edges[i].first.first] = edges[i].second ;
	//}


	//adj-List as pair
	vector<vector<pair<unsigned __int64, unsigned __int64> > > graph;
	__int64 current_vertex = 0 ;
	vector<pair<unsigned __int64, unsigned __int64> > vet ;
	for(unsigned __int64 i=0; i<edges.size(); i++)
	{
		bool new_vertex_entered = edges[i].first.first == current_vertex+1 ?  true : false ;
		
		if(new_vertex_entered)
		{
			current_vertex++;
			graph.push_back(vet) ; 
			vet.clear();
		}
		
		vet.push_back(pair<unsigned __int64, unsigned __int64>(edges[i].first.second,edges[i].second) );		
	}
	graph.push_back(vet) ; 
	//


	
	//E = atoi(argv[argv_to_read++]);
	//std::cin>>E;
	//fin>>E ;
	
	//int e_count = 0;
	//while(e_count<E)
	//{
	//	int u,v,wt;
	//	
	//	//u = atoi(argv[argv_to_read++]);
	//	//v = atoi(argv[argv_to_read++]);
	//	//wt = atoi(argv[argv_to_read++]);
	//	//std::cin>>u>>v>>wt;
	//	fin>>u>>v>>wt;

	//	graph[u][v] = wt; 
	//	graph[v][u] = wt; 
	//	e_count++ ;
	//}

	//making sure that graph is connected
	//for(int i=0;i<V; i++)
	//	graph[i][(i+1)%V] == graph[i][(i+1)%V] == 0 ? 1 : graph[i][(i+1)%V] ;

	unsigned __int64 weight_MST_Prims = primMST(graph);
	return weight_MST_Prims;
}

class connected_components
{
private:
	set<unsigned __int64> vertices;
	//unordered_map<unsigned __int64, &unsigned __int64> vertex_address;
	list<unsigned __int64> connection;
};

class Graph
{
public:

	//unordered_map<unsigned __int64, bool> mark;
	//unordered_map<unsigned __int64, unsigned __int64> edges;
	//vector<list<unsigned __int64> > connected_components;
	//set<unordered_map<unsigned __int64, unsigned __int64> > mst;


	set<unsigned __int64> vertices;
	unordered_map<unsigned __int64, unordered_map<unsigned __int64, unsigned __int64> > edgeList;
	vector<pair<pair<unsigned __int64,unsigned __int64>, unsigned __int64 > > sorted_edgeList;
	set<pair<pair<unsigned __int64,unsigned __int64>, unsigned __int64 > > sorted_edgeList_set;

	bool hasVertex(unsigned __int64 v)
	{
		if(vertices.find(v) != vertices.end())
			return true;
		else return false;
	}

	//returns truee if directed edge added successfully
	bool addEdge(unsigned __int64 u, unsigned __int64 v, unsigned __int64 wt)
	{

		//pair<unsigned __int64,unsigned __int64> pair_uv_edge_set(u,v),pair_vu_edge_set(v,u);
		//sorted_edgeList_set.insert(pair<pair<unsigned __int64,unsigned __int64>, unsigned __int64 >(pair_uv_edge_set, wt));
		//sorted_edgeList_set.insert(pair<pair<unsigned __int64,unsigned __int64>, unsigned __int64 >(pair_vu_edge_set, wt));

		if(edgeList.find(u) != edgeList.end())
		{
			if(edgeList[u].find(v) == edgeList[u].end())
			{
				edgeList[u].insert(pair<unsigned __int64,unsigned __int64>(v,wt));
				return true;
			}
			else
			{
				//cout<<v<<" is already present while inserting "<<u<<","<<v<<","<<wt<<endl;
				return false;
			}
		}
		else
		{
			//cout<<u<<" is not present in vertices while inserting "<<u<<","<<v<<","<<wt<<endl;
			return false;
		}

		//edgeList[u].insert(pair<unsigned __int64,unsigned __int64>(v,wt));
		//return true;


		//else
		//{
		//	unordered_map<unsigned __int64, unsigned __int64> UV_Edge;
		//	UV_Edge.insert(pair<unsigned __int64, unsigned __int64>(v,wt));
		//	edgeList.insert(pair<unsigned __int64, unordered_map<unsigned __int64, unsigned __int64> >(u, UV_Edge)) ;
		//	unordered_map<unsigned __int64, unsigned __int64> VU_Edge;
		//	VU_Edge.insert(pair<unsigned __int64, unsigned __int64>(u,wt));
		//	edgeList.insert(pair<unsigned __int64, unordered_map<unsigned __int64, unsigned __int64> >(v, VU_Edge)) ;
		//	return true;
		//}

		//for(unsigned __int64 i=0; i<edgeList.size(); i++ )
		//{
		//	if(edgeList[i].find(u) != edgeList[i][0].end())
		//	{
		//		if(edgeList[i][0])
		//	}
		//}
	}
	void addEdge_pair_in_set(unsigned __int64 u, unsigned __int64 v, unsigned __int64 wt)
	{

		pair<unsigned __int64,unsigned __int64> pair_uv_edge_set(u,v),pair_vu_edge_set(v,u);
		sorted_edgeList_set.insert(pair<pair<unsigned __int64,unsigned __int64>, unsigned __int64 >(pair_uv_edge_set, wt));
		sorted_edgeList_set.insert(pair<pair<unsigned __int64,unsigned __int64>, unsigned __int64 >(pair_vu_edge_set, wt));
	}

	//returns true if there is a path from u to v
	bool findPath(unsigned __int64 u, unsigned __int64 v, stack<unsigned __int64> &path, set<unsigned __int64> &mark)
	{
		if(u==v)
		{
			//path.push(v);
			return true ;
		}
		else 
		{
			for(unordered_map<unsigned __int64, unsigned __int64>::iterator it = edgeList[u].begin(); it!=edgeList[u].end(); it++)
			{
				if(mark.find(it->first) == mark.end())
				{
					mark.insert(it->first);
					path.push(it->first);
					if(findPath(it->first,v,path,mark))
						return true;
					else
						path.pop();
				}
			}
			return false;
		}
	}
	bool findPath_using_set(unsigned __int64 u, unsigned __int64 v, stack<unsigned __int64> &path, set<unsigned __int64> &mark)
	{
		//cout<<"findPath_using_set\n";
		if(u==v)
		{
			//path.push(v);
			return true ;
		}
		else 
		{
			set<pair<pair<unsigned __int64,unsigned __int64>, unsigned __int64 > >::iterator it ;
			for(it = sorted_edgeList_set.begin(); it != sorted_edgeList_set.end(); it++)
			{
				if(it->first.first == u )
				{
					break;
				}
			}
			do
			{					
				if(mark.find(it->first.second) == mark.end())
				{
					mark.insert(it->first.second);
					path.push(it->first.second);
					if(findPath(it->first.second,v,path,mark))
						return true;
					else
						path.pop();
				}
				it++;
			}while(it->first.first == u);
			return false;
		}
	}
	
	pair<unsigned __int64, unsigned __int64> findHeavyEdge(unsigned __int64 u, unsigned __int64 v)
	{
		stack<unsigned __int64> s;
		set<unsigned __int64> mark;
		mark.insert(u);
		s.push(u);
		findPath(u,v,s,mark);
		unsigned __int64 wt = 0;
		unsigned __int64 e1,e2, h1,h2;
		e1 = s.top();
		s.pop();
		while(!s.empty())
		{
			e2 = s.top();
			if(edgeList[e1][e2] > wt )
			{
				h1=e1;
				h2=e2;
				wt=edgeList[e1][e2];
			}
			e1 = s.top();
			s.pop();
		}
		return pair<unsigned __int64, unsigned __int64>(h1,h2);
	}
	pair<pair<unsigned __int64, unsigned __int64>,unsigned __int64> findHeavyEdge_for_set(unsigned __int64 u, unsigned __int64 v)
	{
		//cout<<"findHeavyEdge_for_set\n";
		stack<unsigned __int64> s;
		set<unsigned __int64> mark;
		mark.insert(u);
		s.push(u);
		//findPath(u,v,s,mark);
		findPath_using_set(u,v,s,mark);
		unsigned __int64 wt = 0;
		unsigned __int64 e1,e2, h1,h2;
		e1 = s.top();
		s.pop();
		while(!s.empty())
		{
			e2 = s.top();
			if(edgeList[e1][e2] > wt )
			{
				h1=e1;
				h2=e2;
				wt=edgeList[e1][e2];
			}
			e1 = s.top();
			s.pop();
		}
		pair<unsigned __int64, unsigned __int64> edg(h1,h2);
		return pair<pair<unsigned __int64, unsigned __int64>, unsigned __int64>(edg,wt);
	}
	
	void removeHeavyEdge(unsigned __int64 u, unsigned __int64 v)
	{
		pair<unsigned __int64, unsigned __int64> p = findHeavyEdge(u, v);
		edgeList[p.first].erase(p.second);
		edgeList[p.second].erase(p.first);
	}
	void removeHeavyEdge_for_set(unsigned __int64 u, unsigned __int64 v)
	{

	}

	void make_sorted_edgeList()
	{
		for(unordered_map<unsigned __int64, unordered_map<unsigned __int64, unsigned __int64> >::iterator it_u = edgeList.begin();
			it_u != edgeList.end();
			it_u++)
		{
			for(unordered_map<unsigned __int64, unsigned __int64>::iterator it_v = it_u->second.begin();
				it_v != it_u->second.end();
				it_v++)
			{
				pair<unsigned __int64,unsigned __int64> p(it_u->first,it_v->first);
				sorted_edgeList.push_back(pair<pair<unsigned __int64,unsigned __int64>, unsigned __int64 >(p,it_v->second) ) ;
			}
		}
		stable_sort(sorted_edgeList.begin(),sorted_edgeList.end(),sort_edges_by_like_set());
	}

	__int64 weigt_of_graph_using_set()
	{
		__int64 wt = 0;
		for(set<pair<pair<unsigned __int64,unsigned __int64>, unsigned __int64 > >::iterator it = sorted_edgeList_set.begin();
			it != sorted_edgeList_set.end();
			it++)
		{
			wt += it->second;
		}
		wt /= 2;
		return wt;
	}

};

__int64 approx_MST_algo1(vector<pair<pair<unsigned __int64,unsigned __int64>, unsigned __int64 > > &edges)
{
	
	vector<Graph *> connected_components;
	//Graph *g = new Graph;
	//connected_components.push_back(g);
	//unsigned __int64 V_total = v_count, E_total = edges.size();
	
	//cin>>V_total>>E_total;
	//fin>>V_total>>E_total;
	//unsigned __int64 u, v, wt ;
	
	
	//std::cin>>u>>v>>wt;
	//fin>>u>>v>>wt;
	int current_edge_no = 0; 
	do
	{
		unsigned __int64 u, v, wt ;
		u=edges[current_edge_no].first.first, v=edges[current_edge_no].first.second, wt=edges[current_edge_no].second ;			
		current_edge_no++;

		unsigned __int64 g1, g2 ;
		g1 = g2 = connected_components.size();
		for(unsigned __int64 i=0ULL; i<connected_components.size(); i++)
		{
			if(connected_components[i]->hasVertex(u))
				g1 = i ;
			if(connected_components[i]->hasVertex(v))
				g2 = i ;
			if(g1!=connected_components.size() && g2 != connected_components.size())
				break;
		}


		if((connected_components.size()==0) || (g1==g2 && g1==connected_components.size()) )	//edge belongs to no existing components	//adding a new component	
		{																					
			Graph *g = new Graph;
			g->vertices.insert(u);
			g->vertices.insert(v);
			//unordered_map<unsigned __int64, unsigned __int64> UV_Edge;
			//UV_Edge.insert(pair<unsigned __int64, unsigned __int64>(v,wt));
			//g->edgeList.insert(pair<unsigned __int64, unordered_map<unsigned __int64, unsigned __int64> >(u, UV_Edge)) ;
			//unordered_map<unsigned __int64, unsigned __int64> VU_Edge;
			//VU_Edge.insert(pair<unsigned __int64, unsigned __int64>(u,wt));
			//g->edgeList.insert(pair<unsigned __int64, unordered_map<unsigned __int64, unsigned __int64> >(v, VU_Edge)) ;
			
			g->addEdge_pair_in_set(u,v,wt);

			connected_components.push_back(g);
		}
		else if(g1 == g2)																		//edge is in the same component	//add edge and check for loop
		{																						
			//pair<unsigned __int64, unsigned __int64> heavyEdge = connected_components[g1]->findHeavyEdge(u,v);
			//if(connected_components[g1]->edgeList[heavyEdge.first][heavyEdge.second] > wt)							//removing the heavy edge and adding the fresh one
			//{
			//	connected_components[g1]->edgeList[heavyEdge.first].erase(heavyEdge.second);
			//	connected_components[g1]->edgeList[heavyEdge.second].erase(heavyEdge.first);
			//	connected_components[g1]->addEdge(u,v,wt);
			//	connected_components[g1]->addEdge(v,u,wt);
			//}

			pair<pair<unsigned __int64, unsigned __int64>, unsigned __int64> heavyEdge_with_wt = connected_components[g1]->findHeavyEdge_for_set(u,v) ;
			if(heavyEdge_with_wt.second > wt)							//removing the heavy edge and adding the fresh one
			{
				//removing the heavy edge 
				connected_components[g1]->sorted_edgeList_set.erase
					(connected_components[g1]->sorted_edgeList_set.find(heavyEdge_with_wt));			
				pair<unsigned __int64, unsigned __int64>  reverse_of_heavyEdge_without_wt(heavyEdge_with_wt.first.second,heavyEdge_with_wt.first.first) ;
				connected_components[g1]->sorted_edgeList_set.erase
					(pair<pair<unsigned __int64, unsigned __int64>, unsigned __int64>(reverse_of_heavyEdge_without_wt,heavyEdge_with_wt.second));

				//and adding the fresh one
				connected_components[g1]->addEdge_pair_in_set(u,v,wt);
			}

		}
		else if(g1==connected_components.size() || g2==connected_components.size())				//one vertex of the edge belong to a component and the other one is new
		{
			if(g1 == connected_components.size() )
			{
				connected_components[g2]->vertices.insert(u);
				//connected_components[g2]->addEdge(v,u,wt);
				//unordered_map<unsigned __int64, unsigned __int64> p;
				//p.insert(pair<unsigned __int64, unsigned __int64> (v,wt));
				//connected_components[g2]->edgeList.insert( pair<unsigned __int64, unordered_map<unsigned __int64, unsigned __int64> > (u,p) );
				//connected_components[g2]->addEdge(u,v,wt);

				connected_components[g2]->addEdge_pair_in_set(u,v,wt);
				
			}
			else
			{
				connected_components[g1]->vertices.insert(v);
				//connected_components[g1]->addEdge(u,v,wt);
				//connected_components[g1]->addEdge(v,u,wt);

				connected_components[g1]->addEdge_pair_in_set(u,v,wt);

			}
		}
		else																					//edge belong to different components	//merge these components
		{																						
			connected_components[g1]->vertices.insert(connected_components[g2]->vertices.begin(), connected_components[g2]->vertices.end()) ;		
			//connected_components[g1]->addEdge(u,v,wt);	
			//connected_components[g1]->edgeList.insert(connected_components[g2]->edgeList.begin(), connected_components[g2]->edgeList.end());
			//connected_components[g1]->addEdge(v,u,wt);


			connected_components[g1]->addEdge_pair_in_set(u,v,wt);
			connected_components[g1]->sorted_edgeList_set.insert(connected_components[g2]->sorted_edgeList_set.begin(), connected_components[g2]->sorted_edgeList_set.end());

			Graph *g = connected_components[g2];
			connected_components.erase(connected_components.begin()+g2);
			delete(g);
		}
	}while(current_edge_no < edges.size());
	if(connected_components.size() > 1)
	{
		//cout<<"data insufficient"<<endl;
		return 0;
	}
	else
	{
		//cout<<"MST fount"<<endl;
		unsigned __int64 total_weight = 0 ;
		unsigned __int64 total_weight_using_set = 0 ;
		//entering one-by-one in each component
		for(int i=0; i<connected_components.size(); i++)
		{
			total_weight_using_set += connected_components[i]->weigt_of_graph_using_set();
			continue;
			//in this current component
			//entering into it's edge-list(entring into edge-list's component one-by-one)
			for(unordered_map<unsigned __int64, unordered_map<unsigned __int64, unsigned __int64> >::iterator it_j = connected_components[i]->edgeList.begin();
				it_j != connected_components[i]->edgeList.end();
				it_j++)
			{

				//in this instance(vertex u) of current edge-list
				//entering into edges connected to this vertex ( which are stored as (v,wt) )
				for(unordered_map<unsigned __int64, unsigned __int64>::iterator it_k = it_j->second.begin() ; 
					it_k != it_j->second.end();
					it_k++)
				{

					//adding one edge(uv) of this vertex u
					total_weight += it_k->second;
				}				
			}			
		}
		total_weight /= 2 ;
		//cout<<total_weight<<endl;
		return total_weight_using_set ;
	}
}

pair<unsigned __int64, unsigned __int64> Prim_vs_Approx(unsigned __int64 &max_weight, unsigned __int64 &v_count, unsigned __int64 &total_edges)
{

	V = v_count ;
	MAX_NODES = v_count ;
	vector<unsigned __int64> vertices;
	vector<pair<pair<unsigned __int64,unsigned __int64>, unsigned __int64 > > edges;
	generate_random_edges_making_one_connected_component(v_count, edges, total_edges, max_weight);
	total_edges = edges.size();	
	__int64 MST_Prims = Prims_algo_MST(edges, v_count, total_edges) ;
	__int64 MST_Algo1 = approx_MST_algo1(edges);
	return pair<unsigned __int64, unsigned __int64>(MST_Prims, MST_Algo1);
}

pair<pair<unsigned __int64, unsigned __int64>, long double> avg_vector_pair(vector<pair<pair<unsigned __int64, unsigned __int64>, long double> > p_vec)
{
	pair<pair<unsigned __int64, unsigned __int64>, long double> avg;
	for(int i=0; i<p_vec.size(); i++)
	{
		avg.first.first += p_vec[i].first.first;
		avg.first.second += p_vec[i].first.second;
		avg.second += p_vec[i].second;
	}
	avg.first.first /= p_vec.size();
	avg.first.second /= p_vec.size();
	avg.second /= p_vec.size();

	return avg;
}

int main(int argn,char **argv)
{
	////constants
	////------------------------------------------------------
	srand(100);
	unsigned __int64 max_weight = 100 ; //100 tha results nikalte time now 10 hua hai
	unsigned __int64 v_count = 10 ;
	unsigned __int64 total_edges = 25 ; //around 
	V = v_count ;
	MAX_NODES = v_count ;
	////--------------------------------------------------------------------------------------------------------------------------------------

	////generating a random graph
	////------------------------------------------------------

	////unsigned __int64 total_edges = v_count ;
	////int file_conter = 0;
	////added to avoid genreting/reading vertices and edges to file   // can be commented in case reading of vertices and edges are from a file
	//vector<unsigned __int64> vertices;
	////vector<pair<unsigned __int64,unsigned __int64> > all_edges ;
	//vector<pair<pair<unsigned __int64,unsigned __int64>, unsigned __int64 > > edges;
	////complete_graph generated(previous data cleared)
	////generate_complete_graph(v_count,vertices,all_edges,max_weight);	
	////required edges generated(previous data cleared)
	//generate_random_edges_making_one_connected_component(v_count, edges, total_edges, max_weight);
	//total_edges = edges.size();	
	////making sure that graph is connected
	////for(int i=0;i<v_count; i++)
	////	if(edges)
	////		edges[i][(i+1)%V] == graph[i][(i+1)%V] == 0 ? 1 : graph[i][(i+1)%V] ;

	////--------------------------------------------------------------------------------------------------------------------------------------


	////Prim's Algorithm
	////------------------------------------------------------

	////ifstream fin;
	////fin.open("extra.txt");	
	////V = atoi(argv[argv_to_read++]);
	////std::cin>>V;
	////fin>>V;
	//V = v_count ;	
	//unsigned __int64 E;
	//E = total_edges ;
	////int **graph = new int*[V] ; 
	////for(int i=0; i<V; i++)
	////{
	////	int *g = new int[V];
	////	for(int j=0; j<V; j++)
	////	{
	////		//std::cin>>g[j];
	////		g[j] = 0;
	////	}
	////	graph[i] = g;
	////}	
	//////adj-matrix
	////vector<vector<unsigned __int64> > graph;
	////vector<unsigned __int64> vet;
	////for(unsigned __int64 i=0; i<V; i++)
	////{
	////	vet.push_back(0);
	////}
	////for(unsigned __int64 i=0; i<V; i++)
	////{
	////	graph.push_back(vet);
	////}
	//////
	//////adj-List
	////vector<map<unsigned __int64, unsigned __int64> > graph;
	////for(unsigned __int64 i=0; i<V; i++)
	////{	
	////	map<unsigned __int64, unsigned __int64> vet;
	////	vet[(i+1)%V] = 0 ; 
	////	graph.push_back(vet);
	////}
	//////
	////for(unsigned __int64 i=0; i<edges.size(); i++)
	////{
	////	unsigned __int64 aa,ab,ac;
	////	aa=edges[i].first.first;	ab=edges[i].first.second;	ac=edges[i].second;
	////	graph[edges[i].first.first][edges[i].first.second] = edges[i].second ;
	////	graph[edges[i].first.second][edges[i].first.first] = edges[i].second ;
	////}
	////adj-List as pair
	//vector<vector<pair<unsigned __int64, unsigned __int64> > > graph;
	//__int64 current_vertex = 0 ;
	//vector<pair<unsigned __int64, unsigned __int64> > vet ;
	//for(unsigned __int64 i=0; i<edges.size(); i++)
	//{
	//	bool new_vertex_entered = edges[i].first.first == current_vertex+1 ?  true : false ;
	//	
	//	if(new_vertex_entered)
	//	{
	//		current_vertex++;
	//		graph.push_back(vet) ; 
	//		vet.clear();
	//	}
	//	
	//	vet.push_back(pair<unsigned __int64, unsigned __int64>(edges[i].first.second,edges[i].second) );		
	//}
	//graph.push_back(vet);
	////
	////E = atoi(argv[argv_to_read++]);
	////std::cin>>E;
	////fin>>E ;	
	////int e_count = 0;
	////while(e_count<E)
	////{
	////	int u,v,wt;
	////	
	////	//u = atoi(argv[argv_to_read++]);
	////	//v = atoi(argv[argv_to_read++]);
	////	//wt = atoi(argv[argv_to_read++]);
	////	//std::cin>>u>>v>>wt;
	////	fin>>u>>v>>wt;
	////	graph[u][v] = wt; 
	////	graph[v][u] = wt; 
	////	e_count++ ;
	////}
	////making sure that graph is connected
	////for(int i=0;i<V; i++)
	////	graph[i][(i+1)%V] == graph[i][(i+1)%V] == 0 ? 1 : graph[i][(i+1)%V] ;
	//unsigned __int64 weight_MST_Prims = primMST(graph);

	////--------------------------------------------------------------------------------------------------------------------------------------


	////Prim'Approx-Algorithm 1
	////------------------------------------------------------

	//vector<Graph *> connected_components;
	////Graph *g = new Graph;
	////connected_components.push_back(g);
	////unsigned __int64 V_total = v_count, E_total = edges.size();
	//
	////cin>>V_total>>E_total;
	////fin>>V_total>>E_total;
	////unsigned __int64 u, v, wt ;
	//
	//
	////std::cin>>u>>v>>wt;
	////fin>>u>>v>>wt;
	//int current_edge_no = 0; 
	//do
	//{
	//	unsigned __int64 u, v, wt ;
	//	u=edges[current_edge_no].first.first, v=edges[current_edge_no].first.second, wt=edges[current_edge_no].second ;			
	//	current_edge_no++;
	//	unsigned __int64 g1, g2 ;
	//	g1 = g2 = connected_components.size();
	//	for(unsigned __int64 i=0ULL; i<connected_components.size(); i++)
	//	{
	//		if(connected_components[i]->hasVertex(u))
	//			g1 = i ;
	//		if(connected_components[i]->hasVertex(v))
	//			g2 = i ;
	//		if(g1!=connected_components.size() && g2 != connected_components.size())
	//			break;
	//	}
	//	if((connected_components.size()==0) || (g1==g2 && g1==connected_components.size()) )	//edge belongs to no existing components	//adding a new component	
	//	{																					
	//		Graph *g = new Graph;
	//		g->vertices.insert(u);
	//		g->vertices.insert(v);
	//		unordered_map<unsigned __int64, unsigned __int64> UV_Edge;
	//		UV_Edge.insert(pair<unsigned __int64, unsigned __int64>(v,wt));
	//		g->edgeList.insert(pair<unsigned __int64, unordered_map<unsigned __int64, unsigned __int64> >(u, UV_Edge)) ;
	//		unordered_map<unsigned __int64, unsigned __int64> VU_Edge;
	//		VU_Edge.insert(pair<unsigned __int64, unsigned __int64>(u,wt));
	//		g->edgeList.insert(pair<unsigned __int64, unordered_map<unsigned __int64, unsigned __int64> >(v, VU_Edge)) ;
	//		g->addEdge_pair_in_set(u,v,wt);
	//		connected_components.push_back(g);
	//	}
	//	else if(g1 == g2)																		//edge is in the same component	//add edge and check for loop
	//	{																						
	//		pair<unsigned __int64, unsigned __int64> heavyEdge = connected_components[g1]->findHeavyEdge(u,v);
	//		if(connected_components[g1]->edgeList[heavyEdge.first][heavyEdge.second] > wt)							//removing the heavy edge and adding the fresh one
	//		{
	//			connected_components[g1]->edgeList[heavyEdge.first].erase(heavyEdge.second);
	//			connected_components[g1]->edgeList[heavyEdge.second].erase(heavyEdge.first);
	//			connected_components[g1]->addEdge(u,v,wt);
	//			connected_components[g1]->addEdge(v,u,wt);
	//		}
	//	}
	//	else if(g1==connected_components.size() || g2==connected_components.size())				//one vertex of the edge belong to a component and the other one is new
	//	{
	//		if(g1 == connected_components.size() )
	//		{
	//			connected_components[g2]->vertices.insert(u);
	//			connected_components[g2]->addEdge(v,u,wt);
	//			unordered_map<unsigned __int64, unsigned __int64> p;
	//			p.insert(pair<unsigned __int64, unsigned __int64> (v,wt));
	//			connected_components[g2]->edgeList.insert( pair<unsigned __int64, unordered_map<unsigned __int64, unsigned __int64> > (u,p) );
	//			connected_components[g2]->addEdge(u,v,wt);
	//		}
	//		else
	//		{
	//			connected_components[g1]->vertices.insert(v);
	//			connected_components[g1]->addEdge(u,v,wt);
	//			connected_components[g1]->addEdge(v,u,wt);
	//		}
	//	}
	//	else																					//edge belong to different components	//merge these components
	//	{																						
	//		connected_components[g1]->vertices.insert(connected_components[g2]->vertices.begin(), connected_components[g2]->vertices.end()) ;		
	//		connected_components[g1]->addEdge(u,v,wt);	
	//		connected_components[g1]->edgeList.insert(connected_components[g2]->edgeList.begin(), connected_components[g2]->edgeList.end());
	//		connected_components[g1]->addEdge(v,u,wt);
	//		Graph *g = connected_components[g2];
	//		connected_components.erase(connected_components.begin()+g2);
	//		delete(g);
	//	}
	//}while(current_edge_no < edges.size());
	//
	//for(int i=0; i<connected_components.size(); i++)
	//{
	//	connected_components[i]->make_sorted_edgeList();
	//}
	//
	//if(connected_components.size() > 1)
	//{
	//	std::cout<<"data insufficient"<<endl;
	//}
	//else
	//{
	//	std::cout<<"MST fount"<<endl;
	//	unsigned __int64 total_weight = 0 ;
	//	//entering one-by-one in each component
	//	for(int i=0; i<connected_components.size(); i++)
	//	{
	//		//in this current component
	//		//entering into it's edge-list(entring into edge-list's component one-by-one)
	//		for(unordered_map<unsigned __int64, unordered_map<unsigned __int64, unsigned __int64> >::iterator it_j = connected_components[i]->edgeList.begin();
	//			it_j != connected_components[i]->edgeList.end();
	//			it_j++)
	//		{
	//			//in this instance(vertex u) of current edge-list
	//			//entering into edges connected to this vertex ( which are stored as (v,wt) )
	//			for(unordered_map<unsigned __int64, unsigned __int64>::iterator it_k = it_j->second.begin() ; 
	//				it_k != it_j->second.end();
	//				it_k++)
	//			{
	//				//adding one edge(uv) of this vertex u
	//				total_weight += it_k->second;
	//			}				
	//		}			
	//	}
	//	total_weight /= 2 ;
	//	std::cout<<total_weight<<endl;
	//	E = approx_MST_algo1(edges);
	//	cout<<E;
	//	std::cout<<total_weight<<endl;
	//}

 	////--------------------------------------------------------------------------------------------------------------------------------------


	//__int64 MST_P = Prims_algo_MST(edges, v_count, total_edges) ;
	//__int64 MST_Algo1 = approx_MST_algo1(edges);
	
	//vector<ofstream> fout_vector;

	//ofstream fout_0;
	////fout_vector.push_back(fout_0);
	for(v_count = 100; v_count <=100; v_count += 50)
	{
		//string file_name = "output_" + to_string(v_count) + ".txt";
		////fout_vector[fout_vector.size()-1].open(file_name.c_str());
		//fout_0.open(file_name.c_str());
		vector<pair<pair<unsigned __int64, unsigned __int64>, long double> > p_vec;
		unsigned __int64 avg_total_edges = 0;
		for(total_edges = (v_count-1)*(v_count-2)/2; total_edges > 2*v_count-1 ; total_edges -= 100)
		{
			pair<unsigned __int64, unsigned __int64> p = Prim_vs_Approx(max_weight, v_count, total_edges);
			////fout_vector[fout_vector.size()-1] << v_count << " , " << total_edges << " , " << p.first << " , " << p.second << " , " << p.second/p.first << endl ;
			//fout_0 << v_count << " , " << total_edges << " , " << p.first << " , " << p.second << " , " << ((long double)p.second)/((long double)p.first) << endl ;
			
			ofstream fout_local;
			string file_name = "output_max_wt_100_" + to_string(v_count) + "_" + to_string(total_edges) + ".txt";
			fout_local.open(file_name.c_str());
			fout_local<<v_count << " , " << total_edges << " , " << p.first << " , " << p.second << " , " << ((long double)p.second)/((long double)p.first) << endl ;
			fout_local.close();

			p_vec.push_back(pair<pair<unsigned __int64, unsigned __int64>, unsigned __int64>(p,((long double)p.second)/((long double)p.first)) );
			avg_total_edges += total_edges ;
		}
		pair<pair<unsigned __int64, unsigned __int64>, long double> avg = avg_vector_pair(p_vec);
		////fout_vector[fout_vector.size()-1] << v_count << " , " << total_edges << " , " << avg.first.first << " , " << avg.first.second << " , " << avg.second << endl ;
		//fout_0 << v_count << " , " << avg_total_edges/p_vec.size() << " , " << avg.first.first << " , " << avg.first.second << " , " << avg.second << endl ;
		
		////fout_vector[fout_vector.size()-1].close();
		////ofstream fout_next;
		////fout_vector.push_back(fout_next);
		//fout_0.close();
	}

	//for(v_count = 500; v_count<=1000; v_count+=500)
	//{
	//	for(total_edges = v_count-1; total_edges<(v_count*(v_count-1))/2; total_edges+=v_count )
	//	{
	//		V = v_count ;
	//		MAX_NODES = v_count ;
	//		vector<unsigned __int64> vertices;
	//		vector<pair<pair<unsigned __int64,unsigned __int64>, unsigned __int64 > > edges;
	//		generate_random_edges_making_one_connected_component(v_count, edges, total_edges, max_weight);
	//		total_edges = edges.size();	
	//	}
	//}
	//cin>>V;
    return 0;
}