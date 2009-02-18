/************************************************************************
	
	Copyright 2007-2009 Emre Sozer & Patrick Clark Trizila

	Contact: emresozer@freecfd.com , ptrizila@freecfd.com

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
#include "commons.h"
#include <iostream>
#include <fstream>
#include <iomanip>
using namespace std;

#include <cgnslib.h>

extern IndexMaps maps;
extern string int2str(int number) ;

void read_restart(double &time) {
	
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
		// Read pressure
		for (unsigned int c=0;c<ncells[p];++c) {
			// Get local cell id
			id=-1;
			if (maps.cellGlobal2Local.find(partitionMap[p][c])!=maps.cellGlobal2Local.end()) {
				id=maps.cellGlobal2Local[partitionMap[p][c]];
			}
			// If id is negative, that means the cell currently lies on another partition
			if (id>=0) { file >> grid.cell[id].p; } else { file >> dummy; }
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
		// Read temperature
		for (unsigned int c=0;c<ncells[p];++c) {
			id=-1;
			if (maps.cellGlobal2Local.find(partitionMap[p][c])!=maps.cellGlobal2Local.end()) {
				id=maps.cellGlobal2Local[partitionMap[p][c]];
			}
			if (id>=0) { 
				file >> grid.cell[id].T; 
				grid.cell[id].rho=eos.rho(grid.cell[id].p,grid.cell[id].T);
			
			} else { file >> dummy; }
		}
		
		if (TURBULENCE_MODEL!=NONE) {
			// Read k
			for (unsigned int c=0;c<ncells[p];++c) {
				id=-1;
				if (maps.cellGlobal2Local.find(partitionMap[p][c])!=maps.cellGlobal2Local.end()) {
					id=maps.cellGlobal2Local[partitionMap[p][c]];
				}
				if (id>=0) { 
					file >> grid.cell[id].k;
				
				} else { file >> dummy; }
			}
			
			// Read omega
			for (unsigned int c=0;c<ncells[p];++c) {
				id=-1;
				if (maps.cellGlobal2Local.find(partitionMap[p][c])!=maps.cellGlobal2Local.end()) {
					id=maps.cellGlobal2Local[partitionMap[p][c]];
				}
				if (id>=0) { 
					file >> grid.cell[id].omega; 
				
				} else { file >> dummy; }
			}
		}
		
		// Skip connectivity list
		getline(file,data,'\n');
		for (unsigned int c=0;c<ncells[p];++c) getline(file,data,'\n');
	}

	file.close();

}
