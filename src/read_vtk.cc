#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <mpi.h>
#include <sys/stat.h>
#include <sys/types.h>
using namespace std;

#include "grid.h"
#include <cgnslib.h>

extern Grid grid;
extern int np, rank;
extern string int2str(int number) ;

void read_vtk(int timeStep) {

// Read partitionMap.dat
vector<vector<int>> partitionMap;
double proc,ncells;
ifstream file;
file.open("partitionMap.dat",ios::in);
file >> proc;
file >> ncells;
cout << proc << "\t" << ncells << endl;
	

}