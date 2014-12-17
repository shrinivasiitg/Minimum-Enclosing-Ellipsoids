#include <iostream>
#include<fstream>

using namespace std;

int main()
{
	ifstream fin;
	ofstream fout;
	fin.open("float.txt");
	fout.open("int.txt");
	while(fin.good())
	{
		double d;
		fin>>d;
		fout<<(long long)(d-2)<<endl;
	}
	fin.close();
	fout.close();
	return 0;
}