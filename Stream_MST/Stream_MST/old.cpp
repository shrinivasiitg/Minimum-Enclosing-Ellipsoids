//#include<iostream>
//#include<vector>
//#include<iterator>
//#include<algorithm>
//#include<time.h>
//#include<fstream>
//#include<string>
//#include<string.h>
//#include<stdlib.h>
//
//using namespace std;
//
//struct sort_pair_by_first_element
//{
//    bool operator()(const pair<__int64,__int64> &left, const pair<__int64,__int64> &right) 
//	{
//		return left.first < right.first;
//    }
//};
//struct sort_pair_by_second_element
//{
//    bool operator()(const pair<__int64,__int64> &left, const pair<__int64,__int64> &right) 
//	{
//        return left.second < right.second;
//    }
//};
//
//struct sort_edges_by_first_vertex
//{
//    bool operator()(const pair<pair<__int64,__int64>,__int64> &left, const pair<pair<__int64,__int64>,__int64> &right) 
//	{
//		return left.first.first < right.first.first;
//    }
//};
//struct sort_edges_by_second_vertex
//{
//    bool operator()(const pair<pair<__int64,__int64>,__int64> &left, const pair<pair<__int64,__int64>,__int64> &right) 
//	{
//		return left.first.second < right.first.second;
//    }
//};
//
//__int64 myrandom (__int64 i) 
//{ 
//	return rand()%i;
//}
//
//
//int main()
//{
//	__int64 v_count = 10 ;
//	fstream fout_results[13];
//	for(int i=0; i<=12; i++)
//	{
//		string file_name = "Results_";
//		file_name +=  to_string(i);
//		file_name +=  ".txt" ;
//		fout_results[i].open(file_name.c_str());
//	}
//	
//	do
//	{
//		__int64 total_edges = v_count ;
//		int file_conter = 0;
//		do
//		{
//			//added to avoid genreting/reading vertices and edges to file   // can be commented in case reading of vertices and edges are from a file
//			vector<__int64> vertices;
//			vector<pair<pair<__int64,__int64>, __int64 > > edges;
//
//
//			//random graph generator code
//
//			srand(unsigned(time(0)));
//			__int64 total_vertices = v_count;
//			__int64 max_weight = 100;
//
//			//__int64 total_edges ;
//			//total_edges = rand()%(total_vertices*((total_vertices-1)/2));
//			//total_edges = total_vertices - 1 ;
//			//total_edges = total_edges < total_vertices ? total_edges+total_vertices : total_edges ;
//
//			//ofstream fout;
//			//fout.open("V.txt");
//			for(__int64 i=1; i<=total_vertices; i++)
//			{
//				//fout<<i<<endl;
//				vertices.push_back(i);
//			}
//			//fout.close();
//
//			vector<pair<__int64,__int64> > all_edges ;
//			for(__int64 i=1; i<=total_vertices; i++)
//				for(__int64 j=i+1; j<=total_vertices; j++)
//					all_edges.push_back(pair<__int64,__int64> (i,j));
//
//			//fout.open("E.txt");
//			for(__int64 i=1; i<=total_edges; i++)
//			{
//				__int64 random_edge = rand()%all_edges.size();
//				__int64 wt = rand()%max_weight + 1 ;
//
//				pair<__int64,__int64> p1(all_edges[random_edge].first , all_edges[random_edge].second);
//				//fout<<all_edges[random_edge].first<<" "<<all_edges[random_edge].second<<" "<<wt<<endl;
//				edges.push_back(pair<pair<__int64,__int64>,__int64>(p1,wt));
//
//				pair<__int64,__int64> p2(all_edges[random_edge].second , all_edges[random_edge].first);
//				//fout<<all_edges[random_edge].second<<" "<<all_edges[random_edge].first<<" "<<wt<<endl;
//				edges.push_back(pair<pair<__int64,__int64>,__int64>(p2,wt));
//
//				all_edges.erase(all_edges.begin() + random_edge);
//			}
//			//fout.close();
//
//
//
//
//
//
//
//			//vertices should be entered as continuous natural number strating from 1 	
//			
//			//vector<__int64> vertices;
//			//ifstream fin;
//			//fin.open("V.txt");
//			//while(fin.good())
//			//{
//			//	__int64 no;
//			//	fin>>no;
//			//	vertices.push_back(no);
//			//}
//			//fin.close();
//
//			stable_sort(vertices.begin(),vertices.end());
//
//
//			//edges are given in the form of ((u,v),uv_weight)
//			//((u,v),uv_weight) and ((v,u)uvu_weight) both shuld be given as graph is undirected	
//			
//			//vector<pair<pair<__int64,__int64>, __int64 > > edges;
//			//fin.open("E.txt");
//			//while(fin.good())
//			//{
//			//	__int64 u,v,wt;
//			//	fin>>u>>v>>wt;
//			//	pair<__int64,__int64> p(u,v);
//			//	edges.push_back(pair<pair<__int64,__int64>,__int64> (p,wt));
//			//}
//			//fin.close();
//
//			total_vertices = vertices.size();
//
//	
//			//genrating random numbers
//			vector<__int64> random_numbers;
//			for(__int64 i=1; i<=vertices.size(); i++)
//			{
//				random_numbers.push_back(i);
//			}
//			srand(unsigned(time(0)));
//			random_shuffle ( random_numbers.begin(), random_numbers.end(), myrandom );
//
//			//V-R Stream/Streams
//			vector<pair<__int64,__int64> > V_R;
//			for(__int64 i=0; i<vertices.size(); i++)
//			{
//				V_R.push_back(pair<__int64,__int64> (vertices[i],random_numbers[i]) );
//			}
//
//			vector<pair<__int64,__int64> > V_R_sorted_by_vertices(V_R);
//			stable_sort(V_R_sorted_by_vertices.begin(),V_R_sorted_by_vertices.end(),sort_pair_by_first_element());
//
//			//vector<pair<__int64,__int64> > V_R_sorted_by_random_numbers(V_R);
//			//stable_sort(V_R_sorted_by_random_numbers.begin(),V_R_sorted_by_random_numbers.end(),sort_pair_by_second_element());
//
//			//sorting edges as first vertex
//			vector<pair<pair<__int64,__int64>,__int64 > > sorted_edges(edges);
//			stable_sort(sorted_edges.begin(),sorted_edges.end(),sort_edges_by_first_vertex());
//
//
//			bool change;
//			__int64 count = -1;
//			do{
//
//				//vector<pair<pair<__int64,__int64>,__int64> > R_E;
//				//for(__int64 i=0; i<sorted_edges.size(); i++)
//				//{
//				//	pair<__int64,__int64> p;
//				//	p.first = sorted_edges[i].first.first ;
//				//	//__int64 secod_vertex_index = vertices sorted_edges[i].first.second
//				//	p.second = random_numbers[sorted_edges[i].first.second - 1] ;
//				//	R_E.push_back( pair<pair<__int64,__int64>,__int64> (p,sorted_edges[i].second) ) ;
//				//}
//
//
//				//Creating E' where first component of each edge, vi, is replaced by corresponding random number, ri
//				vector<pair<pair<__int64,__int64>, __int64> > edges_dash(sorted_edges);
//				for(__int64 i=0; i<sorted_edges.size(); i++)
//				{
//					//pair<__int64,__int64> p;
//					//p.first = random_numbers[sorted_edges[i].first.first-1];
//					//p.first = V_R_sorted_by_vertices[sorted_edges[i].first.first-1].second;
//					edges_dash[i].first.first = V_R_sorted_by_vertices[sorted_edges[i].first.first-1].second;
//					//p.second = sorted_edges[i].first.second ;
//					//edges_dash.push_back(pair<pair<__int64,__int64>,__int64> (p,sorted_edges[i].second));
//				}
//
//				stable_sort(edges_dash.begin(),edges_dash.end(),sort_edges_by_second_vertex());
//
//	
//				//Creating new(smallest) lable for each vertex
//				change = false;
//				//vector<pair<__int64,__int64> > new_V_R(V_R_sorted_by_vertices);
//				__int64 min = V_R_sorted_by_vertices[0].second;
//				min = edges_dash[0].first.first < min ? edges_dash[0].first.first : min ;
//				for(__int64 i=1; i<edges_dash.size(); i++)
//				{
//					//new vertex
//						if(edges_dash[i-1].first.second != edges_dash[i].first.second) 
//						{
//							if(V_R_sorted_by_vertices[edges_dash[i-1].first.second - 1].second != min)
//							{
//								change = true;
//								V_R_sorted_by_vertices[edges_dash[i-1].first.second - 1].second = min ;
//							}
//							min = V_R_sorted_by_vertices[edges_dash[i].first.second - 1].second;
//						}
//					min = edges_dash[i].first.first < min ? edges_dash[i].first.first : min ;
//				}
//				if(V_R_sorted_by_vertices[edges_dash[edges_dash.size()-1].first.second - 1].second != min)
//				{
//					change = true;
//					V_R_sorted_by_vertices[edges_dash[edges_dash.size()-1].first.second - 1].second = min ;
//				}
//
//				count++;
//			}while (change);
//
//
//			fout_results[file_conter]<<total_vertices<<","<<total_edges<<","<<count<<endl;
//			
//			//total_edges = 2*total_edges < v_count*((v_count-1)/2) ?  2*total_edges : v_count*((v_count-1)/2) ;
//			total_edges = v_count*((v_count-1)/2) ;
//			
//			file_conter++;
//		}while(total_edges!=v_count*((v_count-1)/2));
//		v_count = v_count + 1000 ;
//	}while(v_count<=1000);
//
//	for(int i=0; i<=12; i++)
//		fout_results[i].close();
//
//	return 0;
//}