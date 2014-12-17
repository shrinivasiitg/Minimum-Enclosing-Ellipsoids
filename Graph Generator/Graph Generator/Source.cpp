////#include<iostream>
////#include<time.h>
////#include<vector>
////#include<fstream>
////
////using namespace std;
////
////int main()
////{
////	srand(time(0));
////	int vertices = 100;
////	int max_weight = 100;
////
////	int edges = rand()%(vertices*((vertices-1)/2));
////	edges = edges < vertices ? edges+vertices : edges ;
////
////	ofstream fout;
////	
////	fout.open("V.txt");
////	for(int i=1; i<=vertices; i++)
////		fout<<i<<endl;
////	fout.close();
////
////	vector<pair<int,int> > all_edges ;
////	for(int i=1; i<=vertices; i++)
////		for(int j=i+1; j<=vertices; j++)
////			all_edges.push_back(pair<int,int> (i,j));
////
////	fout.open("E.txt");
////	for(int i=1; i<=edges; i++)
////	{
////		__int64 random_edge = rand()%all_edges.size();
////		int wt = rand()%max_weight + 1 ;
////		fout<<all_edges[random_edge].first<<" "<<all_edges[random_edge].second<<" "<<wt<<endl;
////		fout<<all_edges[random_edge].second<<" "<<all_edges[random_edge].first<<" "<<wt<<endl;
////		all_edges.erase(all_edges.begin() + random_edge);
////	}
////	fout.close();
////	
////	return 0;
////}