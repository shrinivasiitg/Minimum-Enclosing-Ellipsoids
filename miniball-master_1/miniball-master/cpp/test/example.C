// Synopsis: Example program illustrating how to use the Seb library
//
// Authors: Martin Kutz <kutz@math.fu-berlin.de>,
//          Kaspar Fischer <kf@iaeth.ch>

#include <iostream>
#include <cstdio>
#include<fstream>

#include "Seb.h"
#include "Seb_debug.C" // ... only needed because we use Seb::Timer below

int main(int argn,char **argv) {
  srand(time(NULL));
  typedef double FT;
  typedef Seb::Point<FT> Point;
  typedef std::vector<Point> PointVector;
  typedef Seb::Smallest_enclosing_ball<FT> Miniball;

  using std::cout;
  using std::endl;
  using std::vector;
  using std::string;
  using std::stringstream;
  ////cout << "====================================================" << endl
  ////     << "Seb example" << endl;

  // Check for right number of arguments ...
  ////if (argn < 3) {
  ////  cout << "Usage: " << argv[0] << " number-of-points dimension [boundary]" << endl
  ////       << "If 'boundary' is given, all points will be on the boundary of a sphere." << endl
	 ////<< "====================================================" << endl;
  ////  return 1;
  ////}
  ////cout << "====================================================" << endl;
  ////const int n = std::atoi(argv[1]), d = std::atoi(argv[2]);
  ////const bool on_boundary = argn > 3 && std::string(argv[3]) == "boundary";

  bool on_boundary ;
  int max_cord;
  if(argn > 3)
  {
	max_cord = argv[3] ;
  }

  // Construct n random points in dimension d
  
  //n = 100;
  //d = 100;
  max_cord = 100;
  for(int n=500; n<=10000; n+=500)
  {  
	  for(int d=500; d<=10000; d+=500)
	  {
		  PointVector S;
		  vector<double> coords(d);
		  //srand(clock());
		  //srand(100);
		  //printing points in file
  
		  std::ofstream fout;
		  string file_name = "results_MSB_";
		  stringstream ss;
		  ss<<"results_MSB_"<<n<<"_"<<d<<".txt";
		  //file_name += to_string(n) + "_" to_string(d) + ".txt" ;
		  fout.open(ss.str().c_str());
		  fout<< n << " , " << d << " , " ;

		  srand(100);
		  for (int i=0; i<n; ++i) {

			// Generate coordindates in [-1,1]
			double len = 0;
			for (int j=0; j<d; ++j) {
			  coords[j] = static_cast<FT>(rand()%max_cord);
			  coords[j] *= 1.0;
			  //fout<<coords[j]<<" ";
			  //coords[j] = static_cast<FT>(2.0*rand()/RAND_MAX - 1.0);
			  len += coords[j]*coords[j];
			}
			//fout<<endl;
			// Normalize length to "almost" 1 (makes it harder for the algorithm)
			/*
			if (on_boundary) {
			  const double Wiggle = 1e-2;
			  len = 1/(std::sqrt(len)+Wiggle*rand()/RAND_MAX);
			  for (int j=0; j<d; ++j)
				coords[j] *= len;
			}
			*/
			S.push_back(Point(d,coords.begin()));
		  }
		  //fout.close();
  
		  ////cout << "Starting computation..." << endl
		  ////     << "====================================================" << endl;
		  ////Seb::Timer::instance().start("all");

		  // Compute the miniball by inserting each value
		  Miniball mb(d, S);

		  // Output
		  FT rad = mb.radius();
		  FT rad_squared = mb.squared_radius();
		  fout<<rad<<endl;
		  ////cout << "Running time: " << Seb::Timer::instance().lapse("all") << "s" << endl
		  ////     << "Radius = " << rad << " (squared: " << rad_squared << ")" << endl
		  ////     << "Center:" << endl;
		  ////Miniball::Coordinate_iterator center_it = mb.center_begin();
		  ////for (int j=0; j<d; ++j) 
		  ////  cout << "  " << center_it[j] << endl;
		  ////cout << "=====================================================" << endl;

		  ////mb.verify();
		  ////cout << "=====================================================" << endl;
		  fout.close();
	  }
  }
}



