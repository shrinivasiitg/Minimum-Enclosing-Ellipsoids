//#include<algorithm>
//#include<fstream>
//#include<iostream>
//#include<iterator>
//#include<limits.h>
//#include<list>
//#include<map>
//#include<math.h>
//#include<stdio.h>
//#include<string>
//#include<stdlib.h>
//#include<time.h>
//#include<vector>
//#include<set>
//#include<stack>
//#include<unordered_map>
//
//using namespace std;
// 
//// Number of vertices in the graph
//unsigned __int64 V = 5 ;
//unsigned __int64 MAX_NODES = 200;
//
//
//struct sort_pair_by_first_element
//{
//    bool operator()(const pair<unsigned __int64,unsigned __int64> &left, const pair<unsigned __int64,unsigned __int64> &right) 
//	{
//		return left.first < right.first;
//    }
//};
//struct sort_pair_by_second_element
//{
//    bool operator()(const pair<unsigned __int64,unsigned __int64> &left, const pair<unsigned __int64,unsigned __int64> &right) 
//	{
//        return left.second < right.second;
//    }
//};
//
//struct sort_edges_by_like_set
//{
//    bool operator()(const pair<pair<unsigned __int64,unsigned __int64>,unsigned __int64> &left, const pair<pair<unsigned __int64,unsigned __int64>,unsigned __int64> &right) 
//	{
//		if(left.first.first == right.first.first)
//			return left.first.second < right.first.second;
//		else
//			return left.first.first < right.first.first;
//    }
//};
//
//// A utility function to find the vertex with minimum key value, from
//// the set of vertices not yet included in MST
//unsigned __int64 minKey(unsigned __int64 key[], bool mstSet[])
//{
//	// Initialize min value
//	unsigned __int64 min = INT_MAX, min_index;
//    
//	for (unsigned __int64 v = 0; v < V; v++)
//		if (mstSet[v] == false && key[v] < min)
//			min = key[v], min_index = v;
//
//	return min_index;
//}
//
//unsigned __int64 printMST(unsigned __int64 parent[], unsigned __int64 n, vector<vector<pair<unsigned __int64, unsigned __int64> > > &graph)
//{
//	//printf("Edge   Weight\n");
//	//	for (unsigned __int64 i = 1; i < V; i++)
//	//		printf("%d - %d    %d \n", parent[i], i, graph[i][parent[i]]);
//
//	printf("Edge   Weight\n");
//		for (unsigned __int64 i = 1; i < V; i++)
//		{	
//			__int64 find_v_in_u = 0;
//			for(find_v_in_u = 0; find_v_in_u< graph[i].size(); find_v_in_u++)
//				if(graph[i][find_v_in_u].first == parent[i])
//					break;
//			printf("%d - %d    %d \n", parent[i], i, graph[i][find_v_in_u].second);
//		}
//	return 0;
//}
//
//unsigned __int64 findWeight(unsigned __int64 parent[], unsigned __int64 n, vector<vector<pair<unsigned __int64, unsigned __int64> > > &graph)
//{
//	unsigned __int64 weight = 0;
//		for (unsigned __int64 i = 1; i < V; i++)
//		{			
//			__int64 find_v_in_u = 0;
//			for(find_v_in_u = 0; find_v_in_u< graph[i].size(); find_v_in_u++)
//				if(graph[i][find_v_in_u].first == parent[i])
//					break;
//			weight += graph[i][find_v_in_u].second ;
//		}
//			
//	return weight ;
//}
// 
//unsigned __int64 primMST(vector<vector<pair<unsigned __int64, unsigned __int64> > > &graph)
//{
//	unsigned __int64 *parent = new unsigned __int64[V]; // Array to store constructed MST
//	unsigned __int64 *key = new unsigned __int64[V];   // Key values used to pick minimum weight edge in cut
//	bool *mstSet = new bool[V];  // To represent set of vertices not yet included in MST
//
//	// Initialize all keys as INFINITE
//	for (unsigned __int64 i = 0; i < V; i++)
//		key[i] = INT_MAX, mstSet[i] = false;
// 
//	// Always include first 1st vertex in MST.
//	key[0] = 0;     // Make key 0 so that this vertex is picked as first vertex
//	parent[0] = -1; // First node is always root of MST
//
//	// The MST will have V vertices
//	for (unsigned __int64 count = 0; count < V-1; count++)
//	{
//		// Pick thd minimum key vertex from the set of vertices
//		// not yet included in MST
//		unsigned __int64 u = minKey(key, mstSet);
//
//		// Add the picked vertex to the MST Set
//		mstSet[u] = true;
//
//		// Update key value and parent index of the adjacent vertices of
//		// the picked vertex. Consider only those vertices which are not yet
//		// included in MST
//
//		for (unsigned __int64 v = 0; v < V; v++)
//		{
//		// graph[u][v] is non zero only for adjacent vertices of m
//		// mstSet[v] is false for vertices not yet included in MST
//		// Update the key only if graph[u][v] is smaller than key[v]
//			__int64 find_v_in_u = 0;
//			for(find_v_in_u = 0; find_v_in_u< graph[u].size(); find_v_in_u++)
//				if(graph[u][find_v_in_u].first == v)
//				{
//					if (graph[u][find_v_in_u].second && mstSet[v] == false && graph[u][find_v_in_u].second <  key[v])
//						parent[v]  = u, key[v] = graph[u][find_v_in_u].second;
//					break;
//				}
//		}
//	}
//	return findWeight(parent, V, graph);
//}
//
//void generate_random_edges_making_one_connected_component(unsigned __int64 &v_count, vector<pair<pair<unsigned __int64,unsigned __int64>, unsigned __int64 > > &edges, unsigned __int64 &total_edges, unsigned __int64 &max_weight)
//{
//	edges.clear();
//
//	set<pair<unsigned __int64,unsigned __int64> > random_edges;
//	
//	//making sure that it is 1 connected component
//	for(int i=0; i<v_count; i++)
//	{
//		random_edges.insert(pair<unsigned __int64,unsigned __int64>(i,(i+1)%v_count)) ;
//		random_edges.insert(pair<unsigned __int64,unsigned __int64>((i+1)%v_count,i)) ;
//	}
//
//	while(random_edges.size() <= total_edges)
//	{
//		__int64 u, v ;
//		u = rand()%v_count;
//		do
//		{
//			v = rand()%v_count;
//		}while(u==v);
//		random_edges.insert(pair<unsigned __int64,unsigned __int64>(u,v));
//		random_edges.insert(pair<unsigned __int64,unsigned __int64>(v,u));
//	}
//	vector<pair<unsigned __int64,unsigned __int64> > random_edges_vector(random_edges.begin(),random_edges.end());
//	do
//	{
//		__int64 u = random_edges_vector[0].first, v = random_edges_vector[0].second ;		
//		__int64 index = find(random_edges_vector.begin(), random_edges_vector.end(), pair<unsigned __int64,unsigned __int64>(v,u)) - random_edges_vector.begin() ;
//		__int64 wt = rand()%max_weight + 1 ;
//		pair<unsigned __int64,unsigned __int64> pair_uv(u,v), pair_vu(v,u);
//		edges.push_back(pair<pair<unsigned __int64,unsigned __int64>, unsigned __int64 > (pair_uv,wt));
//		edges.push_back(pair<pair<unsigned __int64,unsigned __int64>, unsigned __int64 > (pair_vu,wt));
//		random_edges_vector.erase(random_edges_vector.begin() + index);
//		random_edges_vector.erase(random_edges_vector.begin() + 0);
//	}while(random_edges_vector.size() != 0);
//	stable_sort(edges.begin(),edges.end(),sort_edges_by_like_set());
//}
//
//__int64 Prims_algo_MST(vector<pair<pair<unsigned __int64,unsigned __int64>, unsigned __int64 > > &edges, __int64 v_count, __int64 total_edges)
//{
//	V = v_count ;
//	unsigned __int64 E;
//	E = total_edges ;
//	vector<vector<pair<unsigned __int64, unsigned __int64> > > graph;
//	__int64 current_vertex = 0 ;
//	vector<pair<unsigned __int64, unsigned __int64> > vet ;
//	for(unsigned __int64 i=0; i<edges.size(); i++)
//	{
//		bool new_vertex_entered = edges[i].first.first == current_vertex+1 ?  true : false ;
//		
//		if(new_vertex_entered)
//		{
//			current_vertex++;
//			graph.push_back(vet) ; 
//			vet.clear();
//		}
//		
//		vet.push_back(pair<unsigned __int64, unsigned __int64>(edges[i].first.second,edges[i].second) );		
//	}
//	graph.push_back(vet) ; 
//
//	unsigned __int64 weight_MST_Prims = primMST(graph);
//	return weight_MST_Prims;
//}
//
//class Graph
//{
//public:
//
//	set<unsigned __int64> vertices;
//	unordered_map<unsigned __int64, unordered_map<unsigned __int64, unsigned __int64> > edgeList;
//	vector<pair<pair<unsigned __int64,unsigned __int64>, unsigned __int64 > > sorted_edgeList;
//	set<pair<pair<unsigned __int64,unsigned __int64>, unsigned __int64 > > sorted_edgeList_set;
//
//	bool hasVertex(unsigned __int64 v)
//	{
//		if(vertices.find(v) != vertices.end())
//			return true;
//		else return false;
//	}
//
//	//returns truee if directed edge added successfully
//	bool addEdge(unsigned __int64 u, unsigned __int64 v, unsigned __int64 wt)
//	{
//		if(edgeList.find(u) != edgeList.end())
//		{
//			if(edgeList[u].find(v) == edgeList[u].end())
//			{
//				edgeList[u].insert(pair<unsigned __int64,unsigned __int64>(v,wt));
//				return true;
//			}
//			else
//			{
//				return false;
//			}
//		}
//		else
//		{
//			//cout<<u<<" is not present in vertices while inserting "<<u<<","<<v<<","<<wt<<endl;
//			return false;
//		}
//	}
//	void addEdge_pair_in_set(unsigned __int64 u, unsigned __int64 v, unsigned __int64 wt)
//	{
//
//		pair<unsigned __int64,unsigned __int64> pair_uv_edge_set(u,v),pair_vu_edge_set(v,u);
//		sorted_edgeList_set.insert(pair<pair<unsigned __int64,unsigned __int64>, unsigned __int64 >(pair_uv_edge_set, wt));
//		sorted_edgeList_set.insert(pair<pair<unsigned __int64,unsigned __int64>, unsigned __int64 >(pair_vu_edge_set, wt));
//	}
//
//	//returns true if there is a path from u to v
//	bool findPath(unsigned __int64 u, unsigned __int64 v, stack<unsigned __int64> &path, set<unsigned __int64> &mark)
//	{
//		if(u==v)
//		{
//			return true ;
//		}
//		else 
//		{
//			for(unordered_map<unsigned __int64, unsigned __int64>::iterator it = edgeList[u].begin(); it!=edgeList[u].end(); it++)
//			{
//				if(mark.find(it->first) == mark.end())
//				{
//					mark.insert(it->first);
//					path.push(it->first);
//					if(findPath(it->first,v,path,mark))
//						return true;
//					else
//						path.pop();
//				}
//			}
//			return false;
//		}
//	}
//	bool findPath_using_set(unsigned __int64 u, unsigned __int64 v, stack<unsigned __int64> &path, set<unsigned __int64> &mark)
//	{
//		cout<<"findPath_using_set\n";
//		if(u==v)
//		{
//			return true ;
//		}
//		else 
//		{
//			set<pair<pair<unsigned __int64,unsigned __int64>, unsigned __int64 > >::iterator it ;
//			for(it = sorted_edgeList_set.begin(); it != sorted_edgeList_set.end(); it++)
//			{
//				if(it->first.first == u )
//				{
//					break;
//				}
//			}
//			do
//			{					
//				if(mark.find(it->first.second) == mark.end())
//				{
//					mark.insert(it->first.second);
//					path.push(it->first.second);
//					if(findPath(it->first.second,v,path,mark))
//						return true;
//					else
//						path.pop();
//				}
//				it++;
//			}while(it->first.first == u);
//			return false;
//		}
//	}
//	
//	pair<unsigned __int64, unsigned __int64> findHeavyEdge(unsigned __int64 u, unsigned __int64 v)
//	{
//		stack<unsigned __int64> s;
//		set<unsigned __int64> mark;
//		mark.insert(u);
//		s.push(u);
//		findPath(u,v,s,mark);
//		unsigned __int64 wt = 0;
//		unsigned __int64 e1,e2, h1,h2;
//		e1 = s.top();
//		s.pop();
//		while(!s.empty())
//		{
//			e2 = s.top();
//			if(edgeList[e1][e2] > wt )
//			{
//				h1=e1;
//				h2=e2;
//				wt=edgeList[e1][e2];
//			}
//			e1 = s.top();
//			s.pop();
//		}
//		return pair<unsigned __int64, unsigned __int64>(h1,h2);
//	}
//	pair<pair<unsigned __int64, unsigned __int64>,unsigned __int64> findHeavyEdge_for_set(unsigned __int64 u, unsigned __int64 v)
//	{
//		cout<<"findHeavyEdge_for_set\n";
//		stack<unsigned __int64> s;
//		set<unsigned __int64> mark;
//		mark.insert(u);
//		s.push(u);
//		//findPath(u,v,s,mark);
//		findPath_using_set(u,v,s,mark);
//		unsigned __int64 wt = 0;
//		unsigned __int64 e1,e2, h1,h2;
//		e1 = s.top();
//		s.pop();
//		while(!s.empty())
//		{
//			e2 = s.top();
//			if(edgeList[e1][e2] > wt )
//			{
//				h1=e1;
//				h2=e2;
//				wt=edgeList[e1][e2];
//			}
//			e1 = s.top();
//			s.pop();
//		}
//		pair<unsigned __int64, unsigned __int64> edg(h1,h2);
//		return pair<pair<unsigned __int64, unsigned __int64>, unsigned __int64>(edg,wt);
//	}
//	
//	void removeHeavyEdge(unsigned __int64 u, unsigned __int64 v)
//	{
//		pair<unsigned __int64, unsigned __int64> p = findHeavyEdge(u, v);
//		edgeList[p.first].erase(p.second);
//		edgeList[p.second].erase(p.first);
//	}
//	void removeHeavyEdge_for_set(unsigned __int64 u, unsigned __int64 v)
//	{
//
//	}
//
//	void make_sorted_edgeList()
//	{
//		for(unordered_map<unsigned __int64, unordered_map<unsigned __int64, unsigned __int64> >::iterator it_u = edgeList.begin();
//			it_u != edgeList.end();
//			it_u++)
//		{
//			for(unordered_map<unsigned __int64, unsigned __int64>::iterator it_v = it_u->second.begin();
//				it_v != it_u->second.end();
//				it_v++)
//			{
//				pair<unsigned __int64,unsigned __int64> p(it_u->first,it_v->first);
//				sorted_edgeList.push_back(pair<pair<unsigned __int64,unsigned __int64>, unsigned __int64 >(p,it_v->second) ) ;
//			}
//		}
//		stable_sort(sorted_edgeList.begin(),sorted_edgeList.end(),sort_edges_by_like_set());
//	}
//
//	__int64 weigt_of_graph_using_set()
//	{
//		__int64 wt = 0;
//		for(set<pair<pair<unsigned __int64,unsigned __int64>, unsigned __int64 > >::iterator it = sorted_edgeList_set.begin();
//			it != sorted_edgeList_set.end();
//			it++)
//		{
//			wt += it->second;
//		}
//		wt /= 2;
//		return wt;
//	}
//
//};
//
//__int64 approx_MST_algo1(vector<pair<pair<unsigned __int64,unsigned __int64>, unsigned __int64 > > &edges)
//{
//	
//	vector<Graph *> connected_components;
//	int current_edge_no = 0; 
//	do
//	{
//		unsigned __int64 u, v, wt ;
//		u=edges[current_edge_no].first.first, v=edges[current_edge_no].first.second, wt=edges[current_edge_no].second ;			
//		current_edge_no++;
//
//		unsigned __int64 g1, g2 ;
//		g1 = g2 = connected_components.size();
//		for(unsigned __int64 i=0ULL; i<connected_components.size(); i++)
//		{
//			if(connected_components[i]->hasVertex(u))
//				g1 = i ;
//			if(connected_components[i]->hasVertex(v))
//				g2 = i ;
//			if(g1!=connected_components.size() && g2 != connected_components.size())
//				break;
//		}
//
//
//		if((connected_components.size()==0) || (g1==g2 && g1==connected_components.size()) )	//edge belongs to no existing components	//adding a new component	
//		{																					
//			Graph *g = new Graph;
//			g->vertices.insert(u);
//			g->vertices.insert(v);
//		
//			g->addEdge_pair_in_set(u,v,wt);
//
//			connected_components.push_back(g);
//		}
//		else if(g1 == g2)																		//edge is in the same component	//add edge and check for loop
//		{																						
//			pair<pair<unsigned __int64, unsigned __int64>, unsigned __int64> heavyEdge_with_wt = connected_components[g1]->findHeavyEdge_for_set(u,v) ;
//			if(heavyEdge_with_wt.second > wt)							//removing the heavy edge and adding the fresh one
//			{
//				//removing the heavy edge 
//				connected_components[g1]->sorted_edgeList_set.erase
//					(connected_components[g1]->sorted_edgeList_set.find(heavyEdge_with_wt));			
//				pair<unsigned __int64, unsigned __int64>  reverse_of_heavyEdge_without_wt(heavyEdge_with_wt.first.second,heavyEdge_with_wt.first.first) ;
//				connected_components[g1]->sorted_edgeList_set.erase
//					(pair<pair<unsigned __int64, unsigned __int64>, unsigned __int64>(reverse_of_heavyEdge_without_wt,heavyEdge_with_wt.second));
//
//				//and adding the fresh one
//				connected_components[g1]->addEdge_pair_in_set(u,v,wt);
//			}
//		}
//		else if(g1==connected_components.size() || g2==connected_components.size())				//one vertex of the edge belong to a component and the other one is new
//		{
//			if(g1 == connected_components.size() )
//			{
//				connected_components[g2]->vertices.insert(u);
//				connected_components[g2]->addEdge_pair_in_set(u,v,wt);
//			}
//			else
//			{
//				connected_components[g1]->vertices.insert(v);
//				connected_components[g1]->addEdge_pair_in_set(u,v,wt);
//			}
//		}
//		else																					//edge belong to different components	//merge these components
//		{																						
//			connected_components[g1]->vertices.insert(connected_components[g2]->vertices.begin(), connected_components[g2]->vertices.end()) ;		
//			connected_components[g1]->addEdge_pair_in_set(u,v,wt);
//			connected_components[g1]->sorted_edgeList_set.insert(connected_components[g2]->sorted_edgeList_set.begin(), connected_components[g2]->sorted_edgeList_set.end());
//			Graph *g = connected_components[g2];
//			connected_components.erase(connected_components.begin()+g2);
//			delete(g);
//		}
//	}while(current_edge_no < edges.size());
//	if(connected_components.size() > 1)
//	{
//		cout<<"data insufficient"<<endl;
//		return 0;
//	}
//	else
//	{
//		cout<<"MST fount"<<endl;
//		unsigned __int64 total_weight = 0 ;
//		unsigned __int64 total_weight_using_set = 0 ;
//		for(int i=0; i<connected_components.size(); i++)
//		{
//			total_weight_using_set += connected_components[i]->weigt_of_graph_using_set();
//			continue;
//			for(unordered_map<unsigned __int64, unordered_map<unsigned __int64, unsigned __int64> >::iterator it_j = connected_components[i]->edgeList.begin();
//				it_j != connected_components[i]->edgeList.end();
//				it_j++)
//			{
//				for(unordered_map<unsigned __int64, unsigned __int64>::iterator it_k = it_j->second.begin() ; 
//					it_k != it_j->second.end();
//					it_k++)
//				{
//					total_weight += it_k->second;
//				}				
//			}			
//		}
//		total_weight /= 2 ;
//		cout<<total_weight<<endl;
//		return total_weight_using_set ;
//	}
//}
//
//pair<unsigned __int64, unsigned __int64> Prim_vs_Approx(unsigned __int64 &max_weight, unsigned __int64 &v_count, unsigned __int64 &total_edges)
//{
//
//	V = v_count ;
//	MAX_NODES = v_count ;
//	vector<unsigned __int64> vertices;
//	vector<pair<pair<unsigned __int64,unsigned __int64>, unsigned __int64 > > edges;
//	generate_random_edges_making_one_connected_component(v_count, edges, total_edges, max_weight);
//	total_edges = edges.size();	
//	__int64 MST_Prims = Prims_algo_MST(edges, v_count, total_edges) ;
//	__int64 MST_Algo1 = approx_MST_algo1(edges);
//	return pair<unsigned __int64, unsigned __int64>(MST_Prims, MST_Algo1);
//}
//
//int main(int argn,char **argv)
//{
//	////constants
//	////------------------------------------------------------
//	srand(100);
//	unsigned __int64 max_weight = 100 ;
//	unsigned __int64 v_count = 10 ;
//	unsigned __int64 total_edges = 25 ; //around 
//	V = v_count ;
//	MAX_NODES = v_count ;
//	////--------------------------------------------------------------------------------------------------------------------------------------
//	cin>>V;
//    return 0;
//}