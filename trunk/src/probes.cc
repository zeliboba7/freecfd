/************************************************************************
	
	Copyright 2007-2008 Emre Sozer & Patrick Clark Trizila

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
#include "inputs.h"
#include "probe.h"
using namespace std;

extern InputFile input;
extern vector<Probe> probes;
extern vector<Load> loads;

extern string int2str(int number) ;

void set_probes(void) {

	Probe temp;
	double distance=1.e20;
	double distanceTest;
	int probeCount;
	
	struct { 
		double value; 
		int   rank; 
	} in,out;

	probeCount=input.section("probes").subsection("probe",0).count;
	
	// Loop through all the probes
	for (int p=0;p<probeCount;++p) {
		temp.coord=input.section("probes").subsection("probe",p).get_Vec3D("coord");
		// Find the nearest cell
		for (unsigned int c=0;c<grid.cellCount;++c) {
			distanceTest=fabs(temp.coord-grid.cell[c].centroid);
			if (distanceTest<distance) {
				temp.nearestCell=c;
				distance=distanceTest;
			}
		}
		in.value=distance; in.rank=Rank;
		MPI_Allreduce(&in,&out,1, MPI_DOUBLE_INT,MPI_MINLOC,MPI_COMM_WORLD);
		if (out.rank==Rank) { 
			temp.Rank=Rank;
			temp.id=p;
			probes.push_back(temp);
		}
	}

	// Create probe files
	for (int p=0;p<probeCount;++p) {
		probes[p].fileName="probe"+int2str(probes[p].id+1)+".dat";
		ofstream file;
		file.open((probes[p].fileName).c_str(),ios::out);
		file.close();
	}
		
	return;
}

void set_loads(void) {
	
	if (Rank==0) {
		int loadCount=input.section("loads").subsection("load",0).count;
		// Create load files
		for (int n=0;n<loadCount;++n) {
			Load temp;
			temp.bc=input.section("loads").subsection("load",n).get_int("bc");
			temp.fileName="loads_bc_"+int2str(temp.bc)+".dat";
			loads.push_back(temp);
			ofstream file;
			file.open((loads[n].fileName).c_str(),ios::out);
			file.close();
		}
	}
		
	return;
}