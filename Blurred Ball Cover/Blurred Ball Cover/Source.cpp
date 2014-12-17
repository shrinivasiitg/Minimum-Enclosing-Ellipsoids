//#include <iostream>
//#include<fstream>
//#include <windows.h>
//#include<time.h>
//#include<vector>
//#include<math.h>
////#include<unordered_set>
//#include<set>
//#include<unordered_map>
//#include<stdio.h>
//
//using namespace std;
//
////int DIMENSIONS ;
//
////typedef struct _FILETIME {
////  DWORD dwLowDateTime;
////  DWORD dwHighDateTime;
////} FILETIME, *PFILETIME;
////VOID WINAPI GetSystemTimePreciseAsFileTime(  _Out_  LPFILETIME lpSystemTimeAsFileTime);
////BOOL WINAPI FileTimeToSystemTime(  _In_   const FILETIME *lpFileTime,  _Out_  LPSYSTEMTIME lpSystemTime);
////LARGE_INTEGER KeQueryPerformanceCounter(  _Out_opt_  PLARGE_INTEGER PerformanceFrequency);
//
////void print(vector<double>* v)
////{
////	cout<<v->at(0)<<endl;
////}
//
//double distance(vector<double>* point1, vector<double>* point2)
//{
//	double distance_square = 0.0;
//	for(int i=0; i<point1->size(); i++)
//	{
//		distance_square += (point1->at(i)-point2->at(i) )*(point1->at(i)-point2->at(i) ) ;
//	}
//	//cout << distance_square << "," <<sqrt(distance_square) <<endl ;
//	double sqroot =  sqrt(distance_square);
//	return sqroot ;
//}
//double distance_from_origin(vector<double>* point)
//{
//	vector<double> origin ;
//	for(int i=0; i<point->size(); i++)
//	{
//		origin.push_back(0.0) ;
//	}
//	return distance(point, &origin) ;
//}
//
//struct ball
//{
//	//int d;
//	vector<double> center;
//	double radius;
//	//ball( const vector<double>& vet_center = vector<double>() ,  const double& double_radius = NULL) : center(vet_center), radius(double_radius) {}
//	//friend bool operator< (const ball & l, const ball & r) 
//	//{
//	//	vector<double> left_center(l.center);
//	//	vector<double> right_center(r.center);
//	//	return (distance_from_origin(&left_center) < distance_from_origin(&right_center)) && l.radius < r.radius;
//	//}
//	//friend bool operator== (const ball & l, const ball & r) 
//	//{
//	//	return (l.center == r.center) && (l.radius == r.radius) ;
//	//}
//};
//
//
//struct ball* approx_MEB(vector<vector<double> >* Total_Points, double epsilon)
//{
//	//srand(time(NULL));
//	vector<double> current_center(Total_Points->at(rand()%Total_Points->size() ) );
//	//vector<double> current_center(Total_Points->at(rand()%(Total_Points->size()-1) ) );
//	vector<double> next_center(current_center);
//	
//	for(double i=1.0; i<=( 1/(epsilon*epsilon) ); i++)
//	//for(int i=1; i<=(int)( 1/(epsilon*epsilon) ); i++)
//	{
//		current_center = Total_Points->at(rand()%Total_Points->size()) ;
//		//current_center = Total_Points->at(rand()%(Total_Points->size()-1) ) ;
//		vector<double> max_distant_point(current_center);
//		double max_distance = 0;
//		for(int k=0; k<Total_Points->size(); k++)
//		{
//			double curr_distance = distance(&current_center,&Total_Points->at(k)); 
//			max_distant_point = max_distance > curr_distance ? max_distant_point : Total_Points->at(k) ;
//			max_distance = max_distance > curr_distance ? max_distance : curr_distance ;
//		}
//		for(int j=0; j<current_center.size(); j++)
//		{
//			next_center[j] = current_center[j] + (1/((double)i+1))*(max_distant_point[j]-current_center[j]) ;
//		}
//		current_center = next_center ;
//	}
//	double max_radius = 0 ;
//	for(int k=0; k<Total_Points->size(); k++)
//	{
//		double curr_distance = distance(&current_center,&Total_Points->at(k) ); 
//		max_radius = max_radius > curr_distance ? max_radius : curr_distance ;
//	}
//	struct ball *B;
//	B = new ball ; //(current_center,max_radius) ;
//	//B.d = DIMENSIONS;
//	B->center = current_center ;
//	//B->radius = max_radius ;
//	B->radius = max_radius/(1+epsilon) ;
//
//	return B;
//}
//
//bool is_point_in_Ball(struct ball* B, vector<double>* point)
//{
//	return B->radius >= distance(&(B->center),point) ;
//}
//
//bool is_point_in_epsilon_expansion_of_Ball(struct ball* B, vector<double>* point, double epsilon)
////bool is_point_in_epsilon_expansion_of_Ball(struct ball* B, vector<double>* center, double radius, vector<double>* point, double epsilon)
//{
//	//vector<double> center1(B.center);
//	//return true ;
//	//if((1+epsilon)*(B.radius) >
//	//	distance(center1, point))
//	//{
//	//	return true;
//	//}
//
//	//cout << B->radius << "," << B->center[0] << "," << B->center[1] << "," << point->at(0) << "," << point->at(1) << endl; 
//	//cout << (1+epsilon)*(B->radius) << "," << distance(&(B->center),point) << endl ;
//	return (1+epsilon)*(B->radius) >= distance(&(B->center),point) ;
//}
//
//void sort_points(vector<vector<double> > &points, unsigned __int64 d)
//{
//
//}
//void generate_random_points(vector<vector<double> > &points, unsigned __int64 d, unsigned __int64 &n, unsigned __int64 &max_cord)
//{
//	points.clear();
//	vector<double> coords(d);
//	srand(100);
//	for (int i=0; i<n; ++i) 
//	{
//		for (int j=0; j<d; ++j) 
//		{
//			coords[j] = rand()%max_cord;
//			coords[j] *= 1.0;
//		}
//		points.push_back(coords);		
//	}
//}
//
//
//
////vector<vector<double> > approx_MEB_point_set()
////{
////}
//
//
//vector<pair<set<vector<double> >,struct ball > > approx_MEB(int &DIMENSIONS, int &BUFFER_SIZE, int &STREAM_SIZE, int &MAX_X, double &epsilon, vector<vector<double> > &sorted_points_vector)
//{
//
//
//	//unordered_set<unorder
//	//vector<vector<double> > input_stream;
//	//set<>
//	vector<pair<set<vector<double> >,struct ball > > set_Ks_and_Ball_Bs ;
//	//set<set<vector<double> > > set_Ks ;
//	//vector<struct ball> MEBs_set_Ks;
//	vector<vector<double> > buffer_A ;
//
//	for(int i=0; i<STREAM_SIZE; i++)
//	{
//
//		//getting new point from stream
//		//srand(time(NULL));
//		//vector<double> new_point;
//		//for(int i=0; i<DIMENSIONS; i++)
//		//{
//		//	double p;
//		//	char c;
//		//	//p = rand()%MAX_X;
//		//	fin>>p;
//		//	new_point.push_back(p);
//		//	fout<<p<<" , ";
//		//}
//		//fout<<endl;
//		//--------
//		vector<double> new_point(sorted_points_vector[i]);
//
//		buffer_A.push_back(new_point);
//		if(buffer_A.size() >= BUFFER_SIZE)
//		{
//			if(set_Ks_and_Ball_Bs.size() == 0)
//			{
//				struct ball *b_address = approx_MEB(&buffer_A, epsilon/3);
//				struct ball b_value;
//				b_value.center = b_address->center;
//				b_value.radius = b_address->radius;
//				//MEBs_set_Ks.push_back(b_value);
//				set<vector<double> > current_K_i;
//				//for(int i=0; (i<buffer_A.size() ) && (is_point_in_Ball(b_address,&buffer_A[i]) == false ); i++)
//				for(int i=0; i<buffer_A.size(); i++)
//				//for(int i=0; (i<buffer_A.size() ) && (is_point_in_Ball(b_address,&buffer_A[i]) == true ); i++)
//				{
//					if(is_point_in_Ball(b_address,&buffer_A[i]) == true )
//					{
//						current_K_i.insert(buffer_A[i]);
//						//buffer_A.erase(buffer_A.begin()+ i);
//					}
//				}
//				//set_Ks.insert(current_K_i);
//				set_Ks_and_Ball_Bs.push_back(pair<set<vector<double> >,struct ball >(current_K_i,b_value));
//				//set_Ks_and_Ball_Bs.push_back(pair<set<vector<double> >,struct ball >(buffer_A,b_value));
//			}
//			else
//			{
//				bool any_point_of_buffer_A_is_outside = false ; 
//				for(int i=0; i<buffer_A.size(); i++)
//				{
//					//any_point_of_buffer_A_is_outside = is_point_in_Ball(MEB_K[MEB_K.size()-1], buffer_A[i]) == true ? false : true ;
//				
//					//if(1)
//					if(is_point_in_epsilon_expansion_of_Ball(&(set_Ks_and_Ball_Bs[set_Ks_and_Ball_Bs.size()-1].second), &buffer_A[i], epsilon) == false )
//					{
//						any_point_of_buffer_A_is_outside = true ;
//						break;
//					}
//				}
//				if(any_point_of_buffer_A_is_outside == true)
//				{				
//					//UPDATE is invoked
//					//approx_MEB()
//					//set<vector<double> > new_K;
//					//for(set<set<vector<double> > >::iterator it=set_Ks.begin(); it!=set_Ks.end(); it++)
//					//{
//					//	new_K.insert(it->begin(),it->end());
//					//}
//					//new_K.insert(buffer_A.begin(),buffer_A.end());
//
//					//struct ball *B_for_K_unioun_buffer_A = approx_MEB(&vector<vector<double> >(new_K.begin(),new_K.end()), epsilon/3);
// 
//					//set<vector<double> > new_K;
//					set<vector<double> > K_unioun_buffer_A ;
//					for(int i=0; i<set_Ks_and_Ball_Bs.size(); i++)
//					{
//						//new_K.insert(set_Ks_and_Ball_Bs[i].first.begin(), set_Ks_and_Ball_Bs[i].first.end());
//						K_unioun_buffer_A.insert(set_Ks_and_Ball_Bs[i].first.begin(), set_Ks_and_Ball_Bs[i].first.end());
//					}
//					K_unioun_buffer_A.insert(buffer_A.begin(),buffer_A.end());
//					//new_K.insert(buffer_A.begin(),buffer_A.end());
//
//					//calculating B*
//					struct ball *B_address_for_K_unioun_buffer_A = approx_MEB(&vector<vector<double> >(K_unioun_buffer_A.begin(),K_unioun_buffer_A.end()), epsilon/3);
//					//struct ball *B_address_for_K_unioun_buffer_A = approx_MEB(&vector<vector<double> >(new_K.begin(),new_K.end()), epsilon/3);
//					struct ball B_value_for_K_unioun_buffer_A;
//					B_value_for_K_unioun_buffer_A.center = B_address_for_K_unioun_buffer_A->center;
//					B_value_for_K_unioun_buffer_A.radius = B_address_for_K_unioun_buffer_A->radius;
//
//					//calculating K*
//					set<vector<double> > new_K_star;
//					//for(set<vector<double> >::iterator it=new_K.begin();
//					//	(it!=new_K.end() ) && (is_point_in_Ball(B_address_for_K_unioun_buffer_A, &vector<double>(*it) ) ); 
//					//	it++)
//					//{					
//					//	new_K_star.insert(*it);
//					//}
//
//					for(set<vector<double> >::iterator it=K_unioun_buffer_A.begin(); it!=K_unioun_buffer_A.end() ;it++)
//					//for(set<vector<double> >::iterator it=new_K.begin(); it!=new_K.end() ;it++)
//					{					
//						if(is_point_in_Ball(B_address_for_K_unioun_buffer_A, &vector<double>(*it) ) )
//							new_K_star.insert(*it);
//					}
//
//					//adding K* to KK
//					set_Ks_and_Ball_Bs.push_back(pair<set<vector<double> >,struct ball >(new_K_star,B_value_for_K_unioun_buffer_A));
//
//					//deleting smaller blurr-balls
//					for(int i=0; i<set_Ks_and_Ball_Bs.size(); i++)
//					{
//						if(set_Ks_and_Ball_Bs[i].second.radius <= (epsilon/4)*B_value_for_K_unioun_buffer_A.radius)
//						{
//							set_Ks_and_Ball_Bs.erase(set_Ks_and_Ball_Bs.begin() + i) ;
//						}
//					}
//
//				}
//				////UPDATE is invoked
//				////approx_MEB()
//				////set<vector<double> > new_K;
//				////for(set<set<vector<double> > >::iterator it=set_Ks.begin(); it!=set_Ks.end(); it++)
//				////{
//
//			}
//
//				////	new_K.insert(it->begin(),it->end());
//				////}
//				////new_K.insert(buffer_A.begin(),buffer_A.end());
//
//			//	//struct ball *B_for_K_unioun_buffer_A = approx_MEB(&vector<vector<double> >(new_K.begin(),new_K.end()), epsilon/3);
// 
//			//	set<vector<double> > new_K;
//			//	for(int i=0; i<set_Ks_and_Ball_Bs.size(); i++)
//			//	{
//			//		new_K.insert(set_Ks_and_Ball_Bs[i].first.begin(), set_Ks_and_Ball_Bs[i].first.end());
//			//	}
//			//	new_K.insert(buffer_A.begin(),buffer_A.end());
//
//			//	//calculating B*
//			//	struct ball *B_address_for_K_unioun_buffer_A = approx_MEB(&vector<vector<double> >(new_K.begin(),new_K.end()), epsilon/3);
//			//	struct ball B_value_for_K_unioun_buffer_A;
//			//	B_value_for_K_unioun_buffer_A.center = B_address_for_K_unioun_buffer_A->center;
//			//	B_value_for_K_unioun_buffer_A.radius = B_address_for_K_unioun_buffer_A->radius;
//
//			//	//calculating K*
//			//	set<vector<double> > new_K_star;
//			//	for(set<vector<double> >::iterator it=new_K.begin();
//			//		(it!=new_K.end() ) && (is_point_in_Ball(B_address_for_K_unioun_buffer_A, &vector<double>(*it) ) ); 
//			//		it++)
//			//	{					
//			//		new_K_star.insert(*it);
//			//	}
//
//			//	//adding K* to KK
//			//	set_Ks_and_Ball_Bs.push_back(pair<set<vector<double> >,struct ball >(new_K_star,B_value_for_K_unioun_buffer_A));
//
//			//	deleting smaller blurr-balls
//			//	for(int i=0; 
//			//		i<set_Ks_and_Ball_Bs.size() && (set_Ks_and_Ball_Bs[i].second.radius <= (epsilon/4)*B_value_for_K_unioun_buffer_A.radius); 
//			//		i++)
//			//	{
//			//		set_Ks_and_Ball_Bs.erase(set_Ks_and_Ball_Bs.begin() + i) ;
//			//	}
//			//}
//			buffer_A.clear();
//		}
//	}
//	
//	return set_Ks_and_Ball_Bs;
//}
//
//int main()
//{
//	//srand(time(NULL)) ;
//	srand(164);
//
//	////testing code
//	////------------------------------------------------------
//	////set<vector<double> > A;
//	//set<set<vector<double> > > AB ;
//	//for(int i=0; i<30; i++)
//	//{	
//	//	set<vector<double> > A;
//	//	vector<double> vet1, vet2 ;
//	//	double x1, y1, z1, x2, y2, z2 ;
//	//	cin>>x1>>y1>>z1 ;
//	//	vet1.push_back(x1) ;
//	//	vet1.push_back(y1) ;
//	//	vet1.push_back(z1) ;
//	//	A.insert(vet1);
//	//	cin>>x2>>y2>>z2 ;
//	//	vet2.push_back(x2) ;
//	//	vet2.push_back(y2) ;
//	//	vet2.push_back(z2) ;
//	//	A.insert(vet2);
//	//	if(AB.find(A)!=AB.end())
//	//		cout << "Repet" << endl ;
//	//	AB.insert(A);
//	//	cout << AB.size() << endl ;
//	//}
//	//cin>>DIMENSIONS;
//	////identifies repetions in set
//
//	//DIMENSIONS = 2;
//	//vector<double> v1,v2;
//	//v1.push_back(1.0);
//	//v1.push_back(1.0);
//	//v2.push_back(2.0);
//	//v2.push_back(2.0);
//	//vector<vector<double> > vv;
//	//vv.push_back(v1);
//	//vv.push_back(v2);
//	////print(&v1);
//	//struct ball *BB = approx_MEB(&vv,0.9);
//	//std::cout << BB->center[0] << "\t" << BB->center[1] << endl << BB->radius << "\t" << 1.9*BB->radius<<endl;
//	////cout << is_point_in_epsilon_expansion_of_Ball(BB, &BB->center, BB->radius ,&vector<double>(3,3), 0.9) << endl ;
//	//vector<double> vp;
//	//vp.push_back(2.0) ;
//	//vp.push_back(2.0) ;
//	//if(is_point_in_epsilon_expansion_of_Ball(BB, &vp, 0.9) == true )
//	////if(is_point_in_epsilon_expansion_of_Ball(BB, &vector<double>(30,30), 0.9) == true )
//	//{	//cout << is_point_in_epsilon_expansion_of_Ball(BB, &vp, 0.9) << endl ;
//	//	std::cout << "YUP\n" ; 
//	//}
//	////cout << is_point_in_epsilon_expansion_of_Ball(BB, &BB->center, BB->radius ,&vector<double>(1.5,1.6), 0.9);
//	////cout << is_point_in_epsilon_expansion_of_Ball(BB, &vector<double>(1.5,1.6), 0.9);
//	//std::cin>>DIMENSIONS;
//	////all functions checked
//
//	////max memory usage permission checking
//	//set<set<vector<double> > > all_ks;
//	//set<struct ball > set_bs;
//	//for(int i=0; i<1000; i++)
//	//{
//	//	double x, y, r ;
//	//	cin >> x >> y >> r ;
//	//	//struct ball *b_add = new ball(vector<double>(x,y), r);
//	//	//struct ball b(vector<double>(x,y), r);
//	//	struct ball *b_add = new ball(vector<double>(), NULL);
//	//	struct ball b(vector<double>(), NULL);
//	//	b.center.push_back(x);
//	//	b.center.push_back(y);
//	//	b.radius = r ;
//	//	if(set_bs.find(b) != set_bs.end())
//	//		cout << "repet" <<endl ;
//	//	set_bs.insert(b);
//	//}
//	//cin >> DIMENSIONS;
//	////--------------------------------------------------------------------------------------------------------------------------------------
//	
//	////constants
//	////------------------------------------------------------
//	int DIMENSIONS = 3;
//	int BUFFER_SIZE = 1;		//in terms of no of points
//	int STREAM_SIZE = 1000;	//total input points
//	int MAX_X = 10;			//max coordinate in any dimension
//	double epsilon = 0.1 ;
//	////--------------------------------------------------------------------------------------------------------------------------------------
//
//	////time counters
//	////------------------------------------------------------
//	////method 1
//	//int time_span = 1000000;
//	//LARGE_INTEGER lFreq, lStart;
//	//LARGE_INTEGER lEnd;
//	//double d;
//	//QueryPerformanceFrequency(&lFreq);
//	//QueryPerformanceCounter(&lStart);
//	///* do something ... */
//	//for(int i=0;i<=time_span;i++)
//	//{
//	//}
//	//QueryPerformanceCounter(&lEnd);
//	//d = ((double)lEnd.QuadPart - (double)lStart.QuadPart) / (double)lFreq.QuadPart;
//	//////in micro-seconds (1/1000000)
//	//cout<<d*1000*1000<<endl;
//	//
//	//////method 2
//	//long int before = GetTickCount();
//	//cout<<before<<endl;
//	///* do something ... */
//	//for(int i=0;i<=time_span;i++)
//	//{
//	//}
//	//long int after = GetTickCount();
//	//cout<<after-before<<endl;
//	////--------------------------------------------------------------------------------------------------------------------------------------
//
//
//
//	////Reading from stream
//	////------------------------------------------------------
//	
//		ofstream fout;
//	ifstream fin;
//	fin.open("input_points.txt");
//
//	//sorting the file
//	set<vector<double> > sorted_points;
//	vector<vector<double> > sorted_points_vector;
//	for(int i=0; i<STREAM_SIZE; i++)
//	{
//		vector<double> new_point;
//		for(int i=0; i<DIMENSIONS; i++)
//		{
//			int p;
//			fin>>p;
//			new_point.push_back(p);
//		}
//		sorted_points.insert(new_point);
//		sorted_points_vector.push_back(new_point);
//	}
//	fin.close();
//	fout.open("sorted_points.txt");
//	for(set<vector<double> >::iterator it = sorted_points.begin(); it!= sorted_points.end(); it++)
//	{
//		for(int i=0; i<DIMENSIONS; i++)
//		{
//			fout<<it->at(i)<<" ";		
//		}
//		fout<<endl;
//	}
//	fin.open("sorted_points.txt");
//	//
//
//	fout.open("points.txt");
//	
//	//unordered_set<unorder
//	//vector<vector<double> > input_stream;
//	//set<>
//	vector<pair<set<vector<double> >,struct ball > > set_Ks_and_Ball_Bs ;
//	//set<set<vector<double> > > set_Ks ;
//	//vector<struct ball> MEBs_set_Ks;
//	vector<vector<double> > buffer_A ;
//
//	for(int i=0; i<STREAM_SIZE; i++)
//	{
//
//		//getting new point from stream
//		//srand(time(NULL));
//		vector<double> new_point;
//		for(int i=0; i<DIMENSIONS; i++)
//		{
//			double p;
//			char c;
//			//p = rand()%MAX_X;
//			fin>>p;
//			new_point.push_back(p);
//			fout<<p<<" , ";
//		}
//		fout<<endl;
//		//--------
//
//		buffer_A.push_back(new_point);
//		if(buffer_A.size() >= BUFFER_SIZE)
//		{
//			if(set_Ks_and_Ball_Bs.size() == 0)
//			{
//				struct ball *b_address = approx_MEB(&buffer_A, epsilon/3);
//				struct ball b_value;
//				b_value.center = b_address->center;
//				b_value.radius = b_address->radius;
//				//MEBs_set_Ks.push_back(b_value);
//				set<vector<double> > current_K_i;
//				//for(int i=0; (i<buffer_A.size() ) && (is_point_in_Ball(b_address,&buffer_A[i]) == false ); i++)
//				for(int i=0; i<buffer_A.size(); i++)
//				//for(int i=0; (i<buffer_A.size() ) && (is_point_in_Ball(b_address,&buffer_A[i]) == true ); i++)
//				{
//					if(is_point_in_Ball(b_address,&buffer_A[i]) == true )
//					{
//						current_K_i.insert(buffer_A[i]);
//						//buffer_A.erase(buffer_A.begin()+ i);
//					}
//				}
//				//set_Ks.insert(current_K_i);
//				set_Ks_and_Ball_Bs.push_back(pair<set<vector<double> >,struct ball >(current_K_i,b_value));
//				//set_Ks_and_Ball_Bs.push_back(pair<set<vector<double> >,struct ball >(buffer_A,b_value));
//			}
//			else
//			{
//				bool any_point_of_buffer_A_is_outside = false ; 
//				for(int i=0; i<buffer_A.size(); i++)
//				{
//					//any_point_of_buffer_A_is_outside = is_point_in_Ball(MEB_K[MEB_K.size()-1], buffer_A[i]) == true ? false : true ;
//				
//					//if(1)
//					if(is_point_in_epsilon_expansion_of_Ball(&(set_Ks_and_Ball_Bs[set_Ks_and_Ball_Bs.size()-1].second), &buffer_A[i], epsilon) == false )
//					{
//						any_point_of_buffer_A_is_outside = true ;
//						break;
//					}
//				}
//				if(any_point_of_buffer_A_is_outside == true)
//				{				
//					//UPDATE is invoked
//					//approx_MEB()
//					//set<vector<double> > new_K;
//					//for(set<set<vector<double> > >::iterator it=set_Ks.begin(); it!=set_Ks.end(); it++)
//					//{
//					//	new_K.insert(it->begin(),it->end());
//					//}
//					//new_K.insert(buffer_A.begin(),buffer_A.end());
//
//					//struct ball *B_for_K_unioun_buffer_A = approx_MEB(&vector<vector<double> >(new_K.begin(),new_K.end()), epsilon/3);
// 
//					//set<vector<double> > new_K;
//					set<vector<double> > K_unioun_buffer_A ;
//					for(int i=0; i<set_Ks_and_Ball_Bs.size(); i++)
//					{
//						//new_K.insert(set_Ks_and_Ball_Bs[i].first.begin(), set_Ks_and_Ball_Bs[i].first.end());
//						K_unioun_buffer_A.insert(set_Ks_and_Ball_Bs[i].first.begin(), set_Ks_and_Ball_Bs[i].first.end());
//					}
//					K_unioun_buffer_A.insert(buffer_A.begin(),buffer_A.end());
//					//new_K.insert(buffer_A.begin(),buffer_A.end());
//
//					//calculating B*
//					struct ball *B_address_for_K_unioun_buffer_A = approx_MEB(&vector<vector<double> >(K_unioun_buffer_A.begin(),K_unioun_buffer_A.end()), epsilon/3);
//					//struct ball *B_address_for_K_unioun_buffer_A = approx_MEB(&vector<vector<double> >(new_K.begin(),new_K.end()), epsilon/3);
//					struct ball B_value_for_K_unioun_buffer_A;
//					B_value_for_K_unioun_buffer_A.center = B_address_for_K_unioun_buffer_A->center;
//					B_value_for_K_unioun_buffer_A.radius = B_address_for_K_unioun_buffer_A->radius;
//
//					//calculating K*
//					set<vector<double> > new_K_star;
//					//for(set<vector<double> >::iterator it=new_K.begin();
//					//	(it!=new_K.end() ) && (is_point_in_Ball(B_address_for_K_unioun_buffer_A, &vector<double>(*it) ) ); 
//					//	it++)
//					//{					
//					//	new_K_star.insert(*it);
//					//}
//
//					for(set<vector<double> >::iterator it=K_unioun_buffer_A.begin(); it!=K_unioun_buffer_A.end() ;it++)
//					//for(set<vector<double> >::iterator it=new_K.begin(); it!=new_K.end() ;it++)
//					{					
//						if(is_point_in_Ball(B_address_for_K_unioun_buffer_A, &vector<double>(*it) ) )
//							new_K_star.insert(*it);
//					}
//
//					//adding K* to KK
//					set_Ks_and_Ball_Bs.push_back(pair<set<vector<double> >,struct ball >(new_K_star,B_value_for_K_unioun_buffer_A));
//
//					//deleting smaller blurr-balls
//					for(int i=0; i<set_Ks_and_Ball_Bs.size(); i++)
//					{
//						if(set_Ks_and_Ball_Bs[i].second.radius <= (epsilon/4)*B_value_for_K_unioun_buffer_A.radius)
//						{
//							set_Ks_and_Ball_Bs.erase(set_Ks_and_Ball_Bs.begin() + i) ;
//						}
//					}
//
//				}
//				////UPDATE is invoked
//				////approx_MEB()
//				////set<vector<double> > new_K;
//				////for(set<set<vector<double> > >::iterator it=set_Ks.begin(); it!=set_Ks.end(); it++)
//				////{
//
//			}
//
//				////	new_K.insert(it->begin(),it->end());
//				////}
//				////new_K.insert(buffer_A.begin(),buffer_A.end());
//
//			//	//struct ball *B_for_K_unioun_buffer_A = approx_MEB(&vector<vector<double> >(new_K.begin(),new_K.end()), epsilon/3);
// 
//			//	set<vector<double> > new_K;
//			//	for(int i=0; i<set_Ks_and_Ball_Bs.size(); i++)
//			//	{
//			//		new_K.insert(set_Ks_and_Ball_Bs[i].first.begin(), set_Ks_and_Ball_Bs[i].first.end());
//			//	}
//			//	new_K.insert(buffer_A.begin(),buffer_A.end());
//
//			//	//calculating B*
//			//	struct ball *B_address_for_K_unioun_buffer_A = approx_MEB(&vector<vector<double> >(new_K.begin(),new_K.end()), epsilon/3);
//			//	struct ball B_value_for_K_unioun_buffer_A;
//			//	B_value_for_K_unioun_buffer_A.center = B_address_for_K_unioun_buffer_A->center;
//			//	B_value_for_K_unioun_buffer_A.radius = B_address_for_K_unioun_buffer_A->radius;
//
//			//	//calculating K*
//			//	set<vector<double> > new_K_star;
//			//	for(set<vector<double> >::iterator it=new_K.begin();
//			//		(it!=new_K.end() ) && (is_point_in_Ball(B_address_for_K_unioun_buffer_A, &vector<double>(*it) ) ); 
//			//		it++)
//			//	{					
//			//		new_K_star.insert(*it);
//			//	}
//
//			//	//adding K* to KK
//			//	set_Ks_and_Ball_Bs.push_back(pair<set<vector<double> >,struct ball >(new_K_star,B_value_for_K_unioun_buffer_A));
//
//			//	deleting smaller blurr-balls
//			//	for(int i=0; 
//			//		i<set_Ks_and_Ball_Bs.size() && (set_Ks_and_Ball_Bs[i].second.radius <= (epsilon/4)*B_value_for_K_unioun_buffer_A.radius); 
//			//		i++)
//			//	{
//			//		set_Ks_and_Ball_Bs.erase(set_Ks_and_Ball_Bs.begin() + i) ;
//			//	}
//			//}
//			buffer_A.clear();
//		}
//	}
//	////--------------------------------------------------------------------------------------------------------------------------------------
//
//
//
//	int wait_check;
//	cin>>wait_check;
//	return 0;
//}