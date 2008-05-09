#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <mpi.h>
using namespace std;

#include "grid.h"
#include <cgnslib.h>

extern Grid grid;
extern int np, rank;
extern string int2str(int number) ;

void write_tec(int timeStep, double time) {

	ofstream file;
	string fileName="./output/out"+int2str(timeStep)+".tec";
	// Proc 0 creates the output file and writes variable list
	if (rank==0) {
		file.open((fileName).c_str(),ios::out); 
		file << "VARIABLES = \"x\", \"y\", \"z\",\"rho\",\"u\",\"v\",\"w\",\"p\" " << endl;
		//file << "VARIABLES = \"x\", \"y\", \"z\" " << endl;
	} else {
		file.open((fileName).c_str(),ios::app);
	}

	// Write header (each proc has its own zone)
	file << "ZONE, T=\"Partition " << rank << "\"" ;
	file << ", N=" << grid.nodeCount << ", E=" << grid.cellCount << ", VARLOCATION=([4-8]=CellCentered)" << endl;
	file << "DATAPACKING=BLOCK, ZONETYPE=FEBRICK, SOLUTIONTIME=" << time << endl;

	// Write coordinates
	for (unsigned int i=0;i<3;++i) {
		for (unsigned int n=0;n<grid.nodeCount;++n) {
			file << setw(16) << setprecision(8) << scientific << grid.node[n].comp[i] << endl;
		}
	}

	// Write variables
	for (unsigned int c=0;c<grid.cellCount;++c) file << setw(16) << setprecision(8) << scientific << grid.cell[c].rho << endl;
	for (unsigned int c=0;c<grid.cellCount;++c) file << setw(16) << setprecision(8) << scientific << grid.cell[c].v.comp[0] << endl;
	for (unsigned int c=0;c<grid.cellCount;++c) file << setw(16) << setprecision(8) << scientific << grid.cell[c].v.comp[1] << endl;
	for (unsigned int c=0;c<grid.cellCount;++c) file << setw(16) << setprecision(8) << scientific << grid.cell[c].v.comp[2] << endl;
	for (unsigned int c=0;c<grid.cellCount;++c) file << setw(16) << setprecision(8) << scientific << grid.cell[c].p << endl;

	// Write coonnectivity
	for (unsigned int c=0;c<grid.cellCount;++c) {
		if (grid.cell[c].type==PENTA_6) {
			file << grid.cell[c].nodes[0]+1 << "\t" ;
			file << grid.cell[c].nodes[1]+1 << "\t" ;
			file << grid.cell[c].nodes[2]+1 << "\t" ;
			file << grid.cell[c].nodes[2]+1 << "\t" ;
			file << grid.cell[c].nodes[3]+1 << "\t" ;
			file << grid.cell[c].nodes[4]+1 << "\t" ;
			file << grid.cell[c].nodes[5]+1 << "\t" ;
			file << grid.cell[c].nodes[5]+1 << "\t" ;
		} else if (grid.cell[c].type==HEXA_8) {
			for (unsigned int i=0;i<8;++i) {
				file << grid.cell[c].nodes[i]+1 << "\t";
			}
		}
		file << endl;
	}
	
	file.close();

}
