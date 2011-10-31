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
#include "rans.h"
#include "hc.h"
#include "commons.h"

extern vector<Grid> grid;
extern InputFile input;
extern vector<NavierStokes> ns;
extern vector<HeatConduction> hc;
extern vector<RANS> rans;
extern vector<bool> turbulent;
extern vector<Variable<double> > dt;
extern vector<int> equations;

void write_restart(int gid,int timeStep,double time) {

	string fileName;
	
	// Create the restart folder
	ofstream file;
	mkdir("./restart",S_IRWXU);
	string dirname="./restart/"+int2str(timeStep);
	mkdir(dirname.c_str(),S_IRWXU);

	// write partition map
	for (int p=0;p<np;++p) {
		if (Rank==p) {
			fileName="./restart/partitionMap_"+int2str(gid+1)+".dat";
			if (Rank==0) { 
				file.open(fileName.c_str());
				file << np << endl;
			}
			else { file.open(fileName.c_str(), ios::app); }
			file << grid[gid].nodeCount << "\t" << grid[gid].cellCount << endl;
			for (int c=0;c<grid[gid].cellCount;++c) file << grid[gid].cell[c].globalId << endl;
			file.close();
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}
	
	// Write time file
	fileName=dirname+"/time.dat";
	if (Rank==0) { 
		file.open(fileName.c_str());
		// Write physical time
		file << time << endl;
		// Write residual normalization information
		if (equations[gid]==NS) {
			file << ns[gid].first_residuals[0] << "\t" << ns[gid].first_residuals[1] << "\t" << ns[gid].first_residuals[2] << endl;
			if (turbulent[gid]) file << rans[gid].first_residuals[0] << "\t" << rans[gid].first_residuals[1] << endl;
			else file << -1. << "\t" << -1. << endl;
		} else {
			file << -1. << "\t" << -1. << "\t" << -1. << endl;	
			file << -1. << "\t" << -1. << endl;	
		}
		if (equations[gid]==HEAT) file << hc[gid].first_residual << endl;
		else file << -1 << endl;
		file.close();
	}		
	
	string gs="."+int2str(gid+1);
	dt[gid].dump_cell_data(dirname+"/dt"+gs);
	if (equations[gid]==NS) {
		ns[gid].write_restart(timeStep);
		if (turbulent[gid]) rans[gid].write_restart(timeStep);
	}
	if (equations[gid]==HEAT) hc[gid].write_restart(timeStep);	
}


