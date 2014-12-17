#include<iostream>
#include<vector>
#include<iterator>
#include<algorithm>
#include<time.h>
#include<fstream>
#include<unordered_map>
#include<string>
#include<string.h>
#include<stdlib.h>
#include<set>

using namespace std;

unsigned __int64 V = 5 ;
unsigned __int64 MAX_NODES = 200;

struct sort_pair_by_first_element
{
    bool operator()(const pair<__int64,__int64> &left, const pair<__int64,__int64> &right) 
	{
		return left.first < right.first;
    }
};
struct sort_pair_by_second_element
{
    bool operator()(const pair<__int64,__int64> &left, const pair<__int64,__int64> &right) 
	{
        return left.second < right.second;
    }
};

struct sort_edges_by_first_vertex
{
    bool operator()(const pair<pair<__int64,__int64>,__int64> &left, const pair<pair<__int64,__int64>,__int64> &right) 
	{
		return left.first.first < right.first.first;
    }
};
struct sort_edges_by_second_vertex
{
    bool operator()(const pair<pair<__int64,__int64>,__int64> &left, const pair<pair<__int64,__int64>,__int64> &right) 
	{
		return left.first.second < right.first.second;
    }
};
struct sort_edges_by_weights
{
    bool operator()(const pair<pair<__int64,__int64>,__int64> &left, const pair<pair<__int64,__int64>,__int64> &right) 
	{
		return left.second < right.second;
    }
};

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

void increse_edges_vertices_by_one(vector<unsigned __int64> *vertices, vector<pair<pair<unsigned __int64,unsigned __int64>, unsigned __int64 > > *edges)
{
	for(int i=0; i<vertices->size(); i++)
		vertices->at(i)++;
	for(int i=0; i<edges->size(); i++)
	{
		edges->at(i).first.first++;
		edges->at(i).first.second++;
	}
}

void reduce_edges_vertices_by_one(vector<unsigned __int64> *vertices, vector<pair<pair<unsigned __int64,unsigned __int64>, unsigned __int64 > > *edges)
{
	for(int i=0; i<vertices->size(); i++)
		vertices->at(i)--;
	for(int i=0; i<edges->size(); i++)
	{
		edges->at(i).first.first--;
		edges->at(i).first.second--;
	}
}

//void erase_zero_from_edges(vector<pair<pair<unsigned __int64,unsigned __int64>, unsigned __int64 > > &edges)
//{
//	vector<pair<pair<unsigned __int64,unsigned __int64>, unsigned __int64 > > new_edges;
//	for(unsigned __int64 i=0; i<edges.size(); i++)
//	{
//		if(edges.at(i).first.first ==0 || edges.at(i).first.second ==0)
//		{
//		}
//		else
//		{
//			new_edges.push_back(edges[i]);
//		}
//	}
//	edges.clear();
//	for(int i=0; i<new_edges.size(); i++)
//	{
//		edges.push_back(new_edges[i]);
//	}
//}

void map_vertices_edges_from_1_to_n(vector<unsigned __int64> *vertices, vector<pair<pair<unsigned __int64,unsigned __int64>, unsigned __int64 > > *edges)
{
	unordered_map<unsigned __int64, unsigned __int64> vet;
	for(int i=0; i<vertices->size(); i++)
	{
		vet[vertices->at(i)] = i+1;
		vertices->at(i) = i+1;
	}

	for(unsigned __int64 i=0; i<edges->size(); i++)
	{
		edges->at(i).first.first = vet[edges->at(i).first.first];
		edges->at(i).first.second= vet[edges->at(i).first.second];
	}

	//erase_zero_from_edges(*edges);
}

void making_edges_even(vector<pair<pair<unsigned __int64,unsigned __int64>, unsigned __int64 > > &edges)
{
	set<pair<pair<unsigned __int64,unsigned __int64>, unsigned __int64 > > even_edges;
	for(int i=0; i<edges.size(); i++)
	{
		even_edges.insert(edges[i]);
		pair<unsigned __int64,unsigned __int64> p(edges[i].first.second, edges[i].first.first);
		even_edges.insert(pair<pair<unsigned __int64,unsigned __int64>, unsigned __int64 >(p,edges[i].second));
	}
	edges.clear();
	for(set<pair<pair<unsigned __int64,unsigned __int64>, unsigned __int64 > >::iterator it=even_edges.begin(); it!=even_edges.end(); it++ )
		edges.push_back(*it);
}

__int64 myrandom (__int64 i) 
{ 
	return rand()%i;
}

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
	std::stable_sort(edges.begin(),edges.end(),sort_edges_by_like_set());


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

//pair<vector<vector<unsigned __int64> >, unsigned __int64> find_connected_components(vector<unsigned __int64> &vertices, vector<pair<pair<unsigned __int64,unsigned __int64>, unsigned __int64 > > &edges, unsigned __int64 &max_weight)
pair<vector<set<unsigned __int64> >, unsigned __int64> find_connected_components(vector<unsigned __int64> &vertices, vector<pair<pair<unsigned __int64,unsigned __int64>, unsigned __int64 > > &edges)
{
	unsigned __int64 total_edges = edges.size();

	//vector<unsigned __int64> vertices;
	//for(int i=0; i<v_count; i++)
	//	vertices.push_back(i);

    //vertices should be entered as continuous natural number strating from 1 	
			
	//vector<__int64> vertices;
	//ifstream fin;
	//fin.open("V.txt");
	//while(fin.good())
	//{
	//	__int64 no;
	//	fin>>no;
	//	vertices.push_back(no);
	//}
	//fin.close();

	//stable_sort(vertices.begin(),vertices.end());


	//edges are given in the form of ((u,v),uv_weight)
	//((u,v),uv_weight) and ((v,u)uvu_weight) both shuld be given as graph is undirected	
			
	//vector<pair<pair<__int64,__int64>, __int64 > > edges;
	//fin.open("E.txt");
	//while(fin.good())
	//{
	//	__int64 u,v,wt;
	//	fin>>u>>v>>wt;
	//	pair<__int64,__int64> p(u,v);
	//	edges.push_back(pair<pair<__int64,__int64>,__int64> (p,wt));
	//}
	//fin.close();

	//total_vertices = vertices.size();
	unsigned __int64 v_count = vertices.size();
	unsigned __int64 total_vertices = v_count ;
	
	//genrating random numbers
	vector<unsigned __int64> random_numbers;
	for(unsigned __int64 i=1; i<=vertices.size(); i++)
	{
		random_numbers.push_back(i);
	}
	//srand(unsigned(time(0)));
	random_shuffle ( random_numbers.begin(), random_numbers.end(), myrandom );

	//V-R Stream/Streams
	vector<pair<__int64,__int64> > V_R;
	for(__int64 i=0; i<vertices.size(); i++)
	{
		V_R.push_back(pair<__int64,__int64> (vertices[i],random_numbers[i]) );
	}

	vector<pair<__int64,__int64> > V_R_sorted_by_vertices(V_R);
	stable_sort(V_R_sorted_by_vertices.begin(),V_R_sorted_by_vertices.end(),sort_pair_by_first_element());

	//vector<pair<__int64,__int64> > V_R_sorted_by_random_numbers(V_R);
	//stable_sort(V_R_sorted_by_random_numbers.begin(),V_R_sorted_by_random_numbers.end(),sort_pair_by_second_element());

	//sorting edges as first vertex
	vector<pair<pair<unsigned __int64,unsigned __int64>, unsigned __int64 > > sorted_edges(edges);
	stable_sort(sorted_edges.begin(),sorted_edges.end(),sort_edges_by_first_vertex());


	bool change;
	__int64 count = -1;			
	do{

		//vector<pair<pair<__int64,__int64>,__int64> > R_E;
		//for(__int64 i=0; i<sorted_edges.size(); i++)
		//{
		//	pair<__int64,__int64> p;
		//	p.first = sorted_edges[i].first.first ;
		//	//__int64 secod_vertex_index = vertices sorted_edges[i].first.second
		//	p.second = random_numbers[sorted_edges[i].first.second - 1] ;
		//	R_E.push_back( pair<pair<__int64,__int64>,__int64> (p,sorted_edges[i].second) ) ;
		//}


		//Creating E' where first component of each edge, vi, is replaced by corresponding random number, ri
		vector<pair<pair<unsigned __int64,unsigned __int64>, unsigned __int64 > > edges_dash(sorted_edges);
		for(__int64 i=0; i<sorted_edges.size(); i++)
		{
			//pair<__int64,__int64> p;
			//p.first = random_numbers[sorted_edges[i].first.first-1];
			//p.first = V_R_sorted_by_vertices[sorted_edges[i].first.first-1].second;
			edges_dash[i].first.first = V_R_sorted_by_vertices[sorted_edges[i].first.first-1].second;
			//p.second = sorted_edges[i].first.second ;
			//edges_dash.push_back(pair<pair<__int64,__int64>,__int64> (p,sorted_edges[i].second));
		}

		stable_sort(edges_dash.begin(),edges_dash.end(),sort_edges_by_second_vertex());

	
		//Creating new(smallest) lable for each vertex
		change = false;
		//vector<pair<__int64,__int64> > new_V_R(V_R_sorted_by_vertices);
		__int64 min = V_R_sorted_by_vertices[0].second;
		min = edges_dash[0].first.first < min ? edges_dash[0].first.first : min ;
		for(__int64 i=1; i<edges_dash.size(); i++)
		{
			//new vertex
				if(edges_dash[i-1].first.second != edges_dash[i].first.second) 
				{
					if(V_R_sorted_by_vertices[edges_dash[i-1].first.second - 1].second != min)
					{
						change = true;
						V_R_sorted_by_vertices[edges_dash[i-1].first.second - 1].second = min ;
					}
					min = V_R_sorted_by_vertices[edges_dash[i].first.second - 1].second;
				}
			min = edges_dash[i].first.first < min ? edges_dash[i].first.first : min ;
		}
		if(V_R_sorted_by_vertices[edges_dash[edges_dash.size()-1].first.second - 1].second != min)
		{
			change = true;
			V_R_sorted_by_vertices[edges_dash[edges_dash.size()-1].first.second - 1].second = min ;
		}

		count++;
	}while (change);

	//pair<vector<vector<unsigned __int64> >, unsigned __int64> answer;
	pair<vector<set<unsigned __int64> >, unsigned __int64> answer;

	answer.second = count;
	unsigned __int64 current_connected_component_number = V_R_sorted_by_vertices[0].second ;

	for(unsigned __int64 i=0; i<V_R_sorted_by_vertices.size(); i++)
	{
		if(answer.first.size()==0)
		{	
			//vector<unsigned __int64> vec;
			set<unsigned __int64> vec;
			
			//vec.push_back(V_R_sorted_by_vertices[0].first);
			vec.insert(V_R_sorted_by_vertices[0].first);

			answer.first.push_back(vec);
			current_connected_component_number = V_R_sorted_by_vertices[0].second ;
		}
		else if (current_connected_component_number != V_R_sorted_by_vertices[i].second)
		{
			//vector<unsigned __int64> vec;
			set<unsigned __int64> vec;
			
			//vec.push_back(V_R_sorted_by_vertices[i].first);
			vec.insert(V_R_sorted_by_vertices[i].first);

			answer.first.push_back(vec);
			current_connected_component_number = V_R_sorted_by_vertices[i].second ;
		}
		else
		{
			//answer.first[answer.first.size()-1].push_back(V_R_sorted_by_vertices[i].first);
			answer.first[answer.first.size()-1].insert(V_R_sorted_by_vertices[i].first);
		}
	}

	return answer ;
}


void removing_heavy_duplicate_edges(vector<pair<pair<unsigned __int64,unsigned __int64>, unsigned __int64 > > &curr_edges)
{
	if(curr_edges.size() <= 1)
		return;

	vector<pair<pair<unsigned __int64,unsigned __int64>, unsigned __int64 > > new_edges;
	pair<pair<unsigned __int64,unsigned __int64>, unsigned __int64 > prev_pair(curr_edges[0]);
	pair<pair<unsigned __int64,unsigned __int64>, unsigned __int64 > min(curr_edges[0]);
	for(unsigned __int64 i=1; i<curr_edges.size(); i++)
	{
		if(prev_pair.first.first==curr_edges[i].first.first && prev_pair.first.second==curr_edges[i].first.second)
		{
			min.second = min.second > curr_edges[i].second ? curr_edges[i].second : min.second ;
		}
		else
		{
			new_edges.push_back(min);

			prev_pair.first.first = min.first.first = curr_edges[i].first.first ;
			prev_pair.first.second= min.first.second= curr_edges[i].first.second;
			prev_pair.second = min.second = curr_edges[i].second ;
		}
	}
	curr_edges.clear();
	for(unsigned __int64 i=1; i<new_edges.size(); i++)
		curr_edges.push_back(new_edges[i]);
}

//pair<unsigned __int64,set<pair<pair<unsigned __int64,unsigned __int64>, unsigned __int64 > > > approx_MST_algo2(vector<unsigned __int64> &vertices, vector<pair<pair<unsigned __int64,unsigned __int64>, unsigned __int64 > > &edges,  unsigned __int64 &max_weight)
unsigned __int64 approx_MST_algo2(vector<unsigned __int64> &vertices, vector<pair<pair<unsigned __int64,unsigned __int64>, unsigned __int64 > > &edges, double &edges_factor)
{
	//if(edges.size()!=0 && find_connected_components(vertices,edges).first.size() == 1 && edges.size()<=100)
	//{
	//	//use prims algorithm
	//	reduce_edges_vertices_by_one(&vertices, &edges);
	//	unsigned __int64 total_edges = edges.size();
	//	unsigned __int64 v_count = vertices.size();
	//	//set<pair<pair<unsigned __int64,unsigned __int64>, unsigned __int64 > > edges_temp(edges.begin(),edges.end());
	//	unsigned __int64 wt = Prims_algo_MST(edges,v_count,total_edges);
	//	increse_edges_vertices_by_one(&vertices, &edges);
	//	return wt;
	//}
	////if(vertices.size() >= (edges.size()/2)+1)
	//else if(vertices.size() <= 3 || vertices.size() >= (edges.size()/2)+1)
	//{
	//	unsigned __int64 wt = 0;

	//	for(int i=0; i<edges.size(); i++)
	//		wt += edges[i].second;
	//	return wt/2;
	//}

	//if(edges.size() <= 4)
	//{
	//	unsigned __int64 wt = 0;

	//	for(int i=0; i<edges.size(); i++)
	//		wt += edges[i].second;
	//	return wt/2;
	//}

	//if(vertices.size() <= 1)
	//{
	//	return 0;
	//}
	//else if(vertices.size() == 2)
	//{
	//	unsigned __int64 wt = 0;
	//	for(int i=0; i<edges.size(); i++)
	//	{
	//		if(edges[i].first.first == vertices[0] && edges[i].first.second == vertices[1])
	//			wt += edges[i].second;
	//	}
	//	return wt/2;
	//}

	if(edges.size()<=4 || vertices.size() <= 3 || vertices.size() >= (edges.size()/2)+1)
	{
		unsigned __int64 wt = 0;

		for(int i=0; i<edges.size(); i++)
			wt += edges[i].second;
		return wt/2;
	}
	else
	{
		unsigned __int64 total_edges = edges.size();

		std::stable_sort(edges.begin(),edges.end(),sort_edges_by_weights());
		//3*5*64 = 960
		vector<pair<pair<unsigned __int64,unsigned __int64>, unsigned __int64 > > half_edges(edges.begin(),edges.begin() + (((int)(960*edges_factor))*edges.size())/960);

		unsigned __int64 half_edges_total = half_edges.size();
		//vector<unsigned __int64> vertices;
		//for(unsigned __int64 i=0; i<v_count; i++)
		//	vertices.push_back(i);

		//if(vertices.size() == 37)
		//{
		//	vertices.size();
		//}

		//std::stable_sort(half_edges.begin(),half_edges.end(),sort_edges_by_like_set());
		//erase_zero_from_edges(edges);
		making_edges_even(half_edges);
		pair<vector<set<unsigned __int64> >, unsigned __int64> connected_component = find_connected_components(vertices,half_edges);

		//if(connected_component.first.size() == 1)
		//{
		//	//pair<unsigned __int64,set<pair<pair<unsigned __int64,unsigned __int64>, unsigned __int64 > > > weight_mst = approx_MST_algo2(vertices,half_edges,max_weight);
		//	map_vertices_edges_from_1_to_n(&vertices,&half_edges);
		//	unsigned __int64 weight_mst = approx_MST_algo2(vertices,half_edges,edges_factor);
		//	return weight_mst ;
		//}
		//else
		//{
			unsigned __int64 total_wt_MST = 0;
			
			//MST foreach connected components
			for(unsigned __int64 cc=0; cc<connected_component.first.size(); cc++)
			{
				set<pair<pair<unsigned __int64,unsigned __int64>, unsigned __int64 > > edges_for_current_component;
				set<unsigned __int64>::iterator current_vertex = connected_component.first[cc].begin();
				for(unsigned __int64 ei=0; ei<half_edges.size()-1; ei++)
				//for(unsigned __int64 ei=0; ei<edges.size()-1; ei++)
				//for(unsigned __int64 i=0; i<edges.size(); i++)
				{
					if(half_edges[ei].first.first == *current_vertex)
					//if(edges[ei].first.first == *current_vertex)
					{
						//edges_for_current_component.insert(edges[ei]);
						//pair<unsigned __int64, unsigned __int64> p(edges[ei].first.second, edges[ei].first.first);
						//edges_for_current_component.insert(pair<pair<unsigned __int64, unsigned __int64>, unsigned __int64>(p,edges[ei].second));
						//if(edges[ei+1].first.first != *current_vertex)
						//	current_vertex++;
						
						edges_for_current_component.insert(half_edges[ei]);
						pair<unsigned __int64, unsigned __int64> p(half_edges[ei].first.second, half_edges[ei].first.first);
						edges_for_current_component.insert(pair<pair<unsigned __int64, unsigned __int64>, unsigned __int64>(p,half_edges[ei].second));
						if(half_edges[ei+1].first.first != *current_vertex)
							current_vertex++;
					}
					if(current_vertex == connected_component.first[cc].end())
						break;
				}
				//vector<unsigned __int64> vt(connected_component.first.at(i).begin() , connected_component.first.at(i).end());
				vector<unsigned __int64> curr_vertices(connected_component.first[cc].begin(), connected_component.first[cc].end());
				vector<pair<pair<unsigned __int64,unsigned __int64>, unsigned __int64 > > curr_edges_temp(edges_for_current_component.begin(), edges_for_current_component.end());
				vector<pair<pair<unsigned __int64,unsigned __int64>, unsigned __int64 > > curr_edges;
				for(unsigned __int64 i=0; i<curr_edges_temp.size(); i++)
				{
					if(connected_component.first[cc].find(curr_edges_temp[i].first.first) != connected_component.first[cc].end()
						&&
						connected_component.first[cc].find(curr_edges_temp[i].first.second) != connected_component.first[cc].end())
						curr_edges.push_back(curr_edges_temp[i]);
				}
				map_vertices_edges_from_1_to_n(&curr_vertices,&curr_edges);

				unsigned __int64 weight_mst = approx_MST_algo2(curr_vertices, curr_edges, edges_factor);
				total_wt_MST += weight_mst ;
			}
			//

			//MST of remaining set
			set<set<unsigned __int64> > conn_comp_vet_set(connected_component.first.begin(), connected_component.first.end());
			set<pair<pair<unsigned __int64,unsigned __int64>, unsigned __int64 > > edges_for_remaining;
			//first vertex will be representative of the set
			set<unsigned __int64> set_new_vertices;
			//for(set<set<unsigned __int64> >::iterator it=conn_comp_vet_set.begin(); it!=conn_comp_vet_set.end(); it++)
			//{
			//	//first vertex is marker of the set
			//	set_new_vertices.insert( *(it->begin()) );
			//}
			vector<pair<pair<unsigned __int64, unsigned __int64>, unsigned __int64> > remaining_edges(edges.begin()+(((int)(960*edges_factor))*edges.size())/960, edges.end());
			making_edges_even(remaining_edges);
			for(unsigned __int64 i=0; i<remaining_edges.size()-1; i++)
			{
				unsigned __int64 u=remaining_edges[i].first.first, v=remaining_edges[i].first.first ;
				set<unsigned __int64>::iterator it_u = conn_comp_vet_set.begin()->find(u);
				set<unsigned __int64>::iterator it_v = conn_comp_vet_set.begin()->find(v);
				unsigned __int64 u_rep, v_rep;
				bool u_found=false, v_found=false ;
				for(set<set<unsigned __int64> >::iterator it=conn_comp_vet_set.begin(); it!=conn_comp_vet_set.end(); it++)
				{
					if(u_found==false)
						it_u = it->find(u);
					
					if(v_found==false)
						it_v = it->find(v);

					if(u_found==false && it_u != it->end())
					{
						u_rep = *(it->begin());
						u_found=true;
					}
					if(v_found==false && it_v != it->end())
					{
						v_rep = *(it->begin());
						v_found=true;
					}

					if(u_found && v_found)
						break;
				}
				if(u_rep != v_rep)
				{
					set_new_vertices.insert(u_rep);
					set_new_vertices.insert(v_rep);
					pair<unsigned __int64, unsigned __int64> uv_p(u_rep,v_rep), vu_p(v_rep,u_rep) ;
					edges_for_remaining.insert(pair<pair<unsigned __int64, unsigned __int64>, unsigned __int64>(uv_p,remaining_edges[i].second));
					edges_for_remaining.insert(pair<pair<unsigned __int64, unsigned __int64>, unsigned __int64>(vu_p,remaining_edges[i].second));
				}
			}

			vector<unsigned __int64> curr_ver(set_new_vertices.begin(), set_new_vertices.end());
			vector<pair<pair<unsigned __int64,unsigned __int64>, unsigned __int64 > > curr_edges(edges_for_remaining.begin(),edges_for_remaining.end());

			removing_heavy_duplicate_edges(curr_edges);

			making_edges_even(curr_edges);
			map_vertices_edges_from_1_to_n(&curr_ver,&curr_edges);
			unsigned __int64 weight_mst = approx_MST_algo2(curr_ver, curr_edges, edges_factor);
			total_wt_MST += weight_mst;
			
			return total_wt_MST ;
		//}	
	}

}


pair<unsigned __int64, unsigned __int64> Prim_vs_Approx2(unsigned __int64 &max_weight, unsigned __int64 &v_count, unsigned __int64 &total_edges, double &edges_factor)
{

	V = v_count ;
	MAX_NODES = v_count ;
	vector<unsigned __int64> vertices;
	for(int i=0; i<v_count; i++)
		vertices.push_back(i);
	vector<pair<pair<unsigned __int64,unsigned __int64>, unsigned __int64 > > edges;
	generate_random_edges_making_one_connected_component(v_count, edges, total_edges, max_weight);
	total_edges = edges.size();	
	__int64 MST_Prims = Prims_algo_MST(edges, v_count, total_edges) ;

	increse_edges_vertices_by_one(&vertices,&edges);
	map_vertices_edges_from_1_to_n(&vertices,&edges);
	__int64 MST_Algo2 = approx_MST_algo2(vertices, edges, edges_factor);
	return pair<unsigned __int64, unsigned __int64>(MST_Prims, MST_Algo2);
}

pair<unsigned __int64, unsigned __int64> avg_vec(vector<pair<unsigned __int64,unsigned __int64> > &avg)
{
	pair<unsigned __int64, unsigned __int64> ans(0,0);
	for(unsigned __int64 i=0; i<avg.size(); i++)
	{
		ans.first += avg[i].first;
		ans.second += avg[i].second;
	}
	ans.first /= avg.size();
	ans.second/= avg.size();
	return ans;
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

void write_into_file(unsigned __int64 &v_count, unsigned __int64 &total_edges, unsigned __int64 &max_weight, double &edges_factor)
{
	srand(100);
	ofstream fout;
	string file_name = "output_10_100.txt";
	fout.open(file_name.c_str());

	//fout_vector.push_back(fout_0);
	for(v_count = 100; v_count <=100; v_count += 100)
	{
		//string file_name = "output_" + to_string(v_count) + ".txt";
		//fout.open(file_name.c_str());

		vector<pair<pair<unsigned __int64, unsigned __int64>, long double> > p_vec;
		unsigned __int64 avg_total_edges = 0;
		for(total_edges = v_count-1; total_edges < (v_count-1)*(v_count-2)/2 ; total_edges += 100)
		{
			pair<unsigned __int64, unsigned __int64> p = Prim_vs_Approx2(max_weight, v_count, total_edges, edges_factor);
			//fout_vector[fout_vector.size()-1] << v_count << " , " << total_edges << " , " << p.first << " , " << p.second << " , " << p.second/p.first << endl ;
			//fout << v_count << " , " << total_edges << " , " << p.first << " , " << p.second << " , " << ((long double)p.second)/((long double)p.first) << endl ;
			fout << v_count << " , " << total_edges << " , " << p.first << " , " << p.second << " , jugaad , " << ((long double)p.second)/((long double)p.first) + 1<< endl ;
			cout << v_count << " , " << total_edges << " , " << p.first << " , " << p.second << " , " << ((long double)p.second)/((long double)p.first) << endl ;
			
			p_vec.push_back(pair<pair<unsigned __int64, unsigned __int64>, unsigned __int64>(p,((long double)p.second)/((long double)p.first)) );
			avg_total_edges += total_edges ;
		}
		pair<pair<unsigned __int64, unsigned __int64>, long double> avg = avg_vector_pair(p_vec);
		////fout_vector[fout_vector.size()-1] << v_count << " , " << total_edges << " , " << avg.first.first << " , " << avg.first.second << " , " << avg.second << endl ;
	//fout << v_count << " , " << avg_total_edges/p_vec.size() << " , " << avg.first.first << " , " << avg.first.second << " , " << avg.second << endl ;
		//
		////fout_vector[fout_vector.size()-1].close();
		////ofstream fout_next;
		////fout_vector.push_back(fout_next);
		//
		////fout.close();

	}
	fout.close();

}

void write_into_files()
{
}

int main()
{
	//srand (100);
	unsigned __int64 v_count = 100 ;
	unsigned __int64 max_weight = 10 ;
	unsigned __int64 total_edges = 2000;
	double edges_factor = 0.5;
	
	
	//ofstream fout;
	//fout.open("results_count_100.txt");
	//for(v_count; v_count<=100; v_count +=10)
	//{
	//	vector<pair<unsigned __int64,unsigned __int64> > avg;
	//	vector<unsigned __int64> vertices;
	//	for(unsigned __int64 i=1; i<=v_count; i++)
	//		vertices.push_back(i);
	//	for(total_edges=v_count-1; total_edges<(v_count-1)*(v_count-2)/2; total_edges += 50)
	//	{
	//		vector<pair<pair<unsigned __int64,unsigned __int64>, unsigned __int64 > > edges;
	//		generate_random_edges_making_one_connected_component(v_count, edges, total_edges,max_weight);
	//		for(unsigned __int64 i=0; i<edges.size(); i++)
	//		{
	//			edges.at(i).first.first++;
	//			edges.at(i).first.second++;
	//		}
	//		pair<vector<set<unsigned __int64> >, unsigned __int64> ans = find_connected_components(vertices,edges);
	//		//pair<vector<vector<unsigned __int64> >, unsigned __int64> ans = find_connected_components(vertices,edges);
	//		fout << v_count << " , " << edges.size() << " , " << ans.second << endl;
	//		avg.push_back(pair<unsigned __int64,unsigned __int64>(edges.size(),ans.second));
	//	}
	//	fout << v_count << " , " << avg_vec(avg).first << " , " << avg_vec(avg).second << endl;
	//}
	//fout.close();

	//write_into_file(v_count,total_edges, max_weight, edges_factor);

	srand(100);
	for(v_count = 1000; v_count<=10000; v_count+=1000)
	{
		for(total_edges = v_count*unsigned __int64(log(v_count)/log(2)); total_edges<=(v_count*v_count)/(unsigned __int64)sqrt(v_count); total_edges+=5*v_count)
		{
			V = v_count ;
			MAX_NODES = v_count ;
			vector<unsigned __int64> vertices;
			for(int i=0; i<v_count; i++)
				vertices.push_back(i);
			vector<pair<pair<unsigned __int64,unsigned __int64>, unsigned __int64 > > edges;
			generate_random_edges_making_one_connected_component(v_count, edges, total_edges, max_weight);
			total_edges = edges.size();	
			increse_edges_vertices_by_one(&vertices,&edges);
			map_vertices_edges_from_1_to_n(&vertices,&edges);

			pair<vector<set<unsigned __int64> >, unsigned __int64> conn_comp( find_connected_components(vertices, edges) );

			ofstream fout;
			string file_name = "iterations_" + to_string(v_count) + "_" + to_string(total_edges) + ".txt";
			fout.open(file_name.c_str());
			fout << v_count << " , " << total_edges << " , " << conn_comp.second << endl;
			fout.close();
		}
	}

	//cin>>V;
	return 0;
}