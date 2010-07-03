/************************************************************************
	
	Copyright 2007-2010 Emre Sozer

	Contact: emresozer@freecfd.com

	This file is a part of Free CFD

	Free CFD is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    any later version.

    Free CFD is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    For a copy of the GNU General Public License,
    see <http://www.gnu.org/licenses/>.

*************************************************************************/
#include <iostream>
#include <fstream>
using namespace std;

#include "utilities.h"
#include "vec3d.h"
#include "grid.h"
#include "inputs.h"
#include "ns.h"
#include "hc.h"
#include "commons.h"

extern vector<Grid> grid;
extern InputFile input;
extern vector<NavierStokes> ns;
extern vector<HeatConduction> hc;
extern vector<Variable<double> > dt;
extern vector<int> equations;

void read_restart(int gid,int restart_step,double time) {

	fstream file;
		
	// Read partitionMap
	string fileName="./restart/partitionMap_"+int2str(gid+1)+".dat";
	int nprocs,cellGlobalId;
	int ncells[np],nnodes[np];
	file.open(fileName.c_str(),ios::in);	
	if (!file.is_open()) {
		if (Rank==0) cerr << "[E] Restart file " << fileName  << " couldn't be opened" << endl;
		exit(1);
	}
	file >> nprocs;
	vector<vector<int> > partitionMap;
	partitionMap.resize(nprocs);
	for (int p=0;p<nprocs;++p) {
		file >> nnodes[p];
		file >> ncells[p];
		partitionMap[p].resize(ncells[p]);
		for (int c=0;c<ncells[p];++c) {
			file >> cellGlobalId;
			partitionMap[p][c]=cellGlobalId;
		}
	}
	file.close();
	
	for (int g=0;g<grid.size();++g) {
		if (equations[gid]==NS) ns[gid].read_restart(restart_step,partitionMap);
	}
	
	return;
}


