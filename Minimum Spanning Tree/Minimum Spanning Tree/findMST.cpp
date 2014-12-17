#include<iostream>
#include<vector>
#include<unordered_map>
#include<time.h>
#include<math.h>
#include<fstream>
#include<list>
#include<set>
#include<iterator>
#include<unordered_map>
#include<stack>

unsigned __int64 MAX_NODES = 200;
using namespace std;

//struct linked_list
//{
//	unsigned __int64 value;
//	struct linked_list *next;
//};

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

	bool hasVertex(unsigned __int64 v)
	{
		if(vertices.find(v) != vertices.end())
			return true;
		else return false;
	}

	//returns truee if directed edge added successfully
	bool addEdge(unsigned __int64 u, unsigned __int64 v, unsigned __int64 wt)
	{
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
	
	void removeHeavyEdge(unsigned __int64 u, unsigned __int64 v)
	{
		pair<unsigned __int64, unsigned __int64> p = findHeavyEdge(u, v);
		edgeList[p.first].erase(p.second);
		edgeList[p.second].erase(p.first);
	}

};

int main()
{
	//ofstream fout;
	//fout.open("output1.csv");
	//double Log2 = log(2.0);
	//fout<<1<<","<<2<<","<<3<<","<<4<<","<<5<<endl; 
	//for(unsigned long long i=1ULL; i<=1000ULL; i++)
	//{
	//	double Log = (double)i*log((double)i) ;
	//	fout<<i<<","<<Log<<","<<(double)(100ULL/i)*Log<<","<<i*i<<","<<i*i*i<<endl; 
	//}
	//fout.close();


	ifstream fin;
	fin.open("extra.txt");

	ofstream fout;
	fout.open("changed.txt");
	while(fin.good())
	{
		int a,b,c;
		fin>>a>>b>>c;
		fout<<a<<" "<<b<<" "<<c<<endl<<b<<" "<<a<<" "<<c<<endl;
	}
	fin.close();
	fout.close();
	fin.open("changed.txt");

	vector<Graph *> connected_components;
	//Graph *g = new Graph;
	//connected_components.push_back(g);
	unsigned __int64 V_total, E_total;
	
	//cin>>V_total>>E_total;
	fin>>V_total>>E_total;
	unsigned __int64 u, v, wt ;
	
	
	//std::cin>>u>>v>>wt;
	fin>>u>>v>>wt;
	
	srand(time(NULL));
	//ends when u=v
	while(fin.good())
	//while(u!=v)
	{
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
			unordered_map<unsigned __int64, unsigned __int64> UV_Edge;
			UV_Edge.insert(pair<unsigned __int64, unsigned __int64>(v,wt));
			g->edgeList.insert(pair<unsigned __int64, unordered_map<unsigned __int64, unsigned __int64> >(u, UV_Edge)) ;
			unordered_map<unsigned __int64, unsigned __int64> VU_Edge;
			VU_Edge.insert(pair<unsigned __int64, unsigned __int64>(u,wt));
			g->edgeList.insert(pair<unsigned __int64, unordered_map<unsigned __int64, unsigned __int64> >(v, VU_Edge)) ;
			connected_components.push_back(g);
		}
		else if(g1 == g2)																		//edge is in the same component	//add edge and check for loop
		{																						
			pair<unsigned __int64, unsigned __int64> heavyEdge = connected_components[g1]->findHeavyEdge(u,v);
			if(connected_components[g1]->edgeList[heavyEdge.first][heavyEdge.second] > wt)							//removing the heavy edge and adding the fresh one
			{
				connected_components[g1]->edgeList[heavyEdge.first].erase(heavyEdge.second);
				connected_components[g1]->edgeList[heavyEdge.second].erase(heavyEdge.first);
				connected_components[g1]->addEdge(u,v,wt);
				connected_components[g1]->addEdge(v,u,wt);
			}
		}
		else if(g1==connected_components.size() || g2==connected_components.size())				//one vertex of the edge belong to a component and the other one is new
		{
			if(g1 == connected_components.size() )
			{
				connected_components[g2]->vertices.insert(u);
				connected_components[g2]->addEdge(v,u,wt);
				unordered_map<unsigned __int64, unsigned __int64> p;
				p.insert(pair<unsigned __int64, unsigned __int64> (v,wt));
				connected_components[g2]->edgeList.insert( pair<unsigned __int64, unordered_map<unsigned __int64, unsigned __int64> > (u,p) );
				connected_components[g2]->addEdge(u,v,wt);
			}
			else
			{
				connected_components[g1]->vertices.insert(v);
				connected_components[g1]->addEdge(u,v,wt);
				connected_components[g1]->addEdge(v,u,wt);
			}
		}
		else																					//edge belong to different components	//merge these components
		{																						
			connected_components[g1]->vertices.insert(connected_components[g2]->vertices.begin(), connected_components[g2]->vertices.end()) ;		
			connected_components[g1]->addEdge(u,v,wt);	
			connected_components[g1]->edgeList.insert(connected_components[g2]->edgeList.begin(), connected_components[g2]->edgeList.end());
			connected_components[g1]->addEdge(v,u,wt);
			Graph *g = connected_components[g2];
			connected_components.erase(connected_components.begin()+g2);
			delete(g);
		}

		//std::cin>>u>>v>>wt;
		fin>>u>>v>>wt;

	}
	fin.close();
	if(connected_components.size() > 1)
	{
		cout<<"data insufficient"<<endl;
	}
	else
	{
		cout<<"MST fount"<<endl;
		unsigned __int64 total_weight = 0 ;

		//entering one-by-one in each component
		for(int i=0; i<connected_components.size(); i++)
		{

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
		cout<<total_weight<<endl;
	}

	//list<int> comp ;
	//comp.push_back(2);
	//list<int>::iterator it = comp.begin();
	//it++;
	//comp.insert(it,3);
	//for(it = comp.begin(); it!=comp.end(); it++)
	//	cout<<*it<<endl;

	//Graph g;
	//unsigned __int64 vert;
	//cin>>vert;
	//while(vert!=0)
	//{
	//	g.vertices.insert(vert);
	//	cin>>vert;
	//}
	//cin>>u>>v>>wt;
	//while (u!=0)
	//{
	//	g.addEdge(u,v,wt);
	//	g.addEdge(v,u,wt);
	//	cin>>u>>v>>wt;
	//}
	//stack<unsigned __int64> s;
	//set<unsigned __int64> mark;
	//cin>>u>>v;
	//mark.insert(u);
	//s.push(u);
	//if(g.findPath(u,v,s,mark))
	//{
	//	while(!s.empty())
	//	{
	//		cout<<s.top()<<"<--";
	//		s.pop();
	//	}
	//}
	//else
	//{
	//	cout<<"No Path";
	//}

	int i;
	cin>>i;
	return 1;
}