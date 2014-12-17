#include<iostream>
#include<fstream>
#include<stdlib.h>

using namespace std;

int main()
{
	ifstream fin;
	fin.open("extra.txt");
	string arg;
	fin>>arg;
	arg = " " ;
	while(fin.good())
	{
		string s;
		fin>>s;
		arg += s;
		arg += " ";
	}
	string cmd = "D:\\BTP\\Prims Algo\\algo.exe " + arg ;
	system(cmd.c_str());
	return 0;
}