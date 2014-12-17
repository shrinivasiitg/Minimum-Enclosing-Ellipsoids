#include <iostream>     // std::cout
#include <algorithm>    // std::random_shuffle
#include <vector>       // std::vector
#include <ctime>        // std::time
#include <cstdlib> 
#include<fstream>

using namespace std;

int myrandom (int i) 
{ 
	return rand()%i;
}

int main()
{
	srand(time(NULL));
	vector<double> data;
	ifstream fin;
	fin.open("input.txt");
	int counter=0;
	while(fin.good())
	{
		double d;
		fin>>d;
		
		d += counter;
		
		data.push_back(d);
		counter++;
	}
	fin.close();
	//random_shuffle ( data.begin(), data.end(), myrandom);
    ofstream fout;
	fout.open("output.txt");
	for(long long i=0; i<data.size(); i++)
	{
			fout<<data[i]<<endl;
	}
	fout<<endl;
	fout.close();
	return 0;
}