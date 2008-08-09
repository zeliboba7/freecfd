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
#include <vector>
#include <iostream>
#include <fstream>
#include "grid.h"
#include "inputs.h"
#include "probe.h"
using namespace std;

extern InputFile input;
extern Grid grid;
extern vector<Probe> probes;
extern int Rank;

extern string int2str(int number) ;

void set_probes(void) {

	Probe temp;
	double distance=1.e20;
	double distanceTest;
	int probeCount;

	probeCount=input.section["probes"].numberedSubsections["probe"].count;
	
	// Loop through all the probes
	for (int p=0;p<probeCount;++p) {
		temp.coord=input.section["probes"].numberedSubsections["probe"].Vec3Ds[p]["coord"];
		// Find the nearest cell
		for (unsigned int c=0;c<grid.cellCount;++c) {
			distanceTest=fabs(temp.coord-grid.cell[c].centroid);
			if (distanceTest<distance) {
				temp.nearestCell=c;
				distance=distanceTest;
			}
		}
		MPI_Allreduce(&distance,&distanceTest,1, MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
		if (distance<=distanceTest) {
			temp.Rank=Rank;
			temp.id=p;
			probes.push_back(temp);
		}
	}

	// Create probe files
	for (int p=0;p<probeCount;++p) {
		if (Rank==probes[p].Rank) {
		string fileName;
		fileName="probe"+int2str(p)+".dat";
		ofstream file;
		//probeFiles.push_back((new ofstream));
		file.open((fileName).c_str(),ios::out);
		file.close();
		}
	}
		
	return;
}
