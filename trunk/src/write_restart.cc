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
#include <string>
#include <iomanip>
#include <mpi.h>
#include <sys/stat.h>
#include <sys/types.h>
using namespace std;

#include <cgnslib.h>
#include "rans.h"
#include "flamelet.h"

extern RANS rans;
extern Flamelet flamelet;

extern string int2str(int number) ;

void write_restart(double time) {

	ofstream file;
	mkdir("./restart",S_IRWXU);
	string dirname="./restart/"+int2str(timeStep);
	mkdir(dirname.c_str(),S_IRWXU);
	string fileName="./restart/"+int2str(timeStep)+"/field.dat";
	// Proc 0 creates the output file and writes variable list
	int nVar=8;
	if (Rank==0) {
		file.open((fileName).c_str(),ios::out); 
		file << "VARIABLES = \"x\", \"y\", \"z\",\"p\",\"u\",\"v\",\"w\",\"T\" " ;
		if (TURBULENCE_MODEL!=NONE) {
			file << "\"k\", \"omega\" " ; nVar+=2;
		}
		if (FLAMELET) {
			file << "\"Mixture Fraction\", \"Mixture Fraction Variance\" " ; nVar+=2;
		}
		file << endl; // DEBUG
	} else {
		file.open((fileName).c_str(),ios::app);
	}
	
	// Write header (each proc has its own zone)
	file << "ZONE, T=\"Partition " << Rank << "\"" ;
	file << ", N=" << grid.nodeCount << ", E=" << grid.cellCount << ", VARLOCATION=([4-" << nVar << "]=CellCentered)" << endl;
	file << "DATAPACKING=BLOCK, ZONETYPE=FEBRICK, SOLUTIONTIME=" << time << endl;

	// Write coordinates
	for (unsigned int i=0;i<3;++i) {
		for (unsigned int n=0;n<grid.nodeCount;++n) {
			file << setw(16) << setprecision(8) << scientific << grid.node[n].comp[i] << endl;
		}
	}

	// Write variables

	for (unsigned int c=0;c<grid.cellCount;++c) file << setw(16) << setprecision(8) << scientific << grid.cell[c].p << endl;
	for (unsigned int c=0;c<grid.cellCount;++c) file << setw(16) << setprecision(8) << scientific << grid.cell[c].v.comp[0] << endl;
	for (unsigned int c=0;c<grid.cellCount;++c) file << setw(16) << setprecision(8) << scientific << grid.cell[c].v.comp[1] << endl;
	for (unsigned int c=0;c<grid.cellCount;++c) file << setw(16) << setprecision(8) << scientific << grid.cell[c].v.comp[2] << endl;
	for (unsigned int c=0;c<grid.cellCount;++c) file << setw(16) << setprecision(8) << scientific << grid.cell[c].T << endl;
	if (TURBULENCE_MODEL!=NONE) {
		for (unsigned int c=0;c<grid.cellCount;++c) file << setw(16) << setprecision(8) << scientific << rans.cell[c].k << endl;
		for (unsigned int c=0;c<grid.cellCount;++c) file << setw(16) << setprecision(8) << scientific << rans.cell[c].omega << endl;
	}
	if (FLAMELET) {
		for (unsigned int c=0;c<grid.cellCount;++c) file << setw(16) << setprecision(8) << scientific << flamelet.cell[c].Z << endl;
		for (unsigned int c=0;c<grid.cellCount;++c) file << setw(16) << setprecision(8) << scientific << flamelet.cell[c].Zvar << endl;
	}
	// Write coonnectivity
	for (unsigned int c=0;c<grid.cellCount;++c) {
		if (grid.cell[c].nodeCount==4) {
			file << grid.cell[c].nodes[0]+1 << "\t" ;
			file << grid.cell[c].nodes[2]+1 << "\t" ;
			file << grid.cell[c].nodes[1]+1 << "\t" ;
			file << grid.cell[c].nodes[1]+1 << "\t" ;
			file << grid.cell[c].nodes[3]+1 << "\t" ;
			file << grid.cell[c].nodes[3]+1 << "\t" ;
			file << grid.cell[c].nodes[3]+1 << "\t" ;
			file << grid.cell[c].nodes[3]+1 << "\t" ;
		}
		else if (grid.cell[c].nodeCount==5) {
			file << grid.cell[c].nodes[0]+1 << "\t" ;
			file << grid.cell[c].nodes[1]+1 << "\t" ;
			file << grid.cell[c].nodes[2]+1 << "\t" ;
			file << grid.cell[c].nodes[3]+1 << "\t" ;
			file << grid.cell[c].nodes[4]+1 << "\t" ;
			file << grid.cell[c].nodes[4]+1 << "\t" ;
			file << grid.cell[c].nodes[4]+1 << "\t" ;
			file << grid.cell[c].nodes[4]+1 << "\t" ;
		}
		else if (grid.cell[c].nodeCount==6) {
			file << grid.cell[c].nodes[0]+1 << "\t" ;
			file << grid.cell[c].nodes[1]+1 << "\t" ;
			file << grid.cell[c].nodes[2]+1 << "\t" ;
			file << grid.cell[c].nodes[2]+1 << "\t" ;
			file << grid.cell[c].nodes[3]+1 << "\t" ;
			file << grid.cell[c].nodes[4]+1 << "\t" ;
			file << grid.cell[c].nodes[5]+1 << "\t" ;
			file << grid.cell[c].nodes[5]+1 << "\t" ;
		} else if (grid.cell[c].nodeCount==8) {
			for (unsigned int i=0;i<8;++i) {
				file << grid.cell[c].nodes[i]+1 << "\t";
			}
		}
		file << endl;
	}
	
	file.close();

	// write partition map
	for (int p=0;p<np;++p) {
		if (Rank==p) {
			fileName="./restart/"+int2str(timeStep)+"/partitionMap.dat";
			if (Rank==0) { file.open(fileName.c_str(), ios::out); file << np << endl;}
			else { file.open(fileName.c_str(), ios::app); }
			file << grid.nodeCount << "\t" << grid.cellCount << endl;
			for (unsigned int c=0;c<grid.cellCount;++c) file << grid.cell[c].globalId << endl;
			file.close();
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}


}
