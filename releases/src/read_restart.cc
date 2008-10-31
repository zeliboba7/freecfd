#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
using namespace std;

#include "grid.h"
#include <cgnslib.h>

extern Grid grid;
extern int np, Rank;
extern IndexMaps maps;
extern string int2str(int number) ;

void read_restart(int restart, double &time) {
	
	fstream file;

	// Read partitionMap
	string fileName="./restart/"+int2str(restart)+"/partitionMap.dat";
	int nprocs,cellGlobalId;
	int ncells[np],nnodes[np];
	file.open(fileName.c_str(),ios::in);
	file >> nprocs;
	vector<int> partitionMap[nprocs];
	for (int p=0;p<nprocs;++p) {
		file >> nnodes[p];
		file >> ncells[p];
		for (unsigned int c=0;c<ncells[p];++c) {
			file >> cellGlobalId;
			partitionMap[p].push_back(cellGlobalId);
		}
	}
	file.close();
	
	fileName="./restart/"+int2str(restart)+"/field.dat";
	string data="";
	double dummy;
	file.open(fileName.c_str(),ios::in);
	// Skip first two header lines
	getline(file,data,'\n'); getline(file,data,'\n');
	// Skip 3 "=" signs and read solution time
	for (unsigned int i=0; i<3; ++i) getline(file,data,'=');
	file >> time;
	// Rewind the file
	file.seekg(0,ios::beg);
	// Skip the variable header line
	getline(file,data,'\n');
	int id;
	for (unsigned int p=0;p<nprocs;++p) {
		// Skip two header lines
		getline(file,data,'\n'); getline(file,data,'\n');
		// Skip node data
		for (unsigned int n=0;n<3*nnodes[p];++n) file >> dummy;
		// Read density
		for (unsigned int c=0;c<ncells[p];++c) {
			// Get local cell id
			id=-1;
			if (maps.cellGlobal2Local.find(partitionMap[p][c])!=maps.cellGlobal2Local.end()) {
				id=maps.cellGlobal2Local[partitionMap[p][c]];
			}
			// If id is negative, that means the cell currently lies on another partition
			if (id>=0) { file >> grid.cell[id].rho; } else { file >> dummy; }
		}
		// Read u-velocity
		for (unsigned int c=0;c<ncells[p];++c) {
			id=-1;
			if (maps.cellGlobal2Local.find(partitionMap[p][c])!=maps.cellGlobal2Local.end()) {
				id=maps.cellGlobal2Local[partitionMap[p][c]];
			}
			if (id>=0) { file >> grid.cell[id].v.comp[0]; } else { file >> dummy; }
		}
		// Read v-velocity
		for (unsigned int c=0;c<ncells[p];++c) {
			id=-1;
			if (maps.cellGlobal2Local.find(partitionMap[p][c])!=maps.cellGlobal2Local.end()) {
				id=maps.cellGlobal2Local[partitionMap[p][c]];
			}
			if (id>=0) { file >> grid.cell[id].v.comp[1]; } else { file >> dummy; }
		}
		// Read w-velocity
		for (unsigned int c=0;c<ncells[p];++c) {
			id=-1;
			if (maps.cellGlobal2Local.find(partitionMap[p][c])!=maps.cellGlobal2Local.end()) {
				id=maps.cellGlobal2Local[partitionMap[p][c]];
			}
			if (id>=0) { file >> grid.cell[id].v.comp[2]; } else { file >> dummy; }
		}		
		// Read pressure
		for (unsigned int c=0;c<ncells[p];++c) {
			id=-1;
			if (maps.cellGlobal2Local.find(partitionMap[p][c])!=maps.cellGlobal2Local.end()) {
				id=maps.cellGlobal2Local[partitionMap[p][c]];
			}
			if (id>=0) { file >> grid.cell[id].p; } else { file >> dummy; }
		}
		// Skip connectivity list
		getline(file,data,'\n');
		for (unsigned int c=0;c<ncells[p];++c) getline(file,data,'\n');
	}

	file.close();

}
