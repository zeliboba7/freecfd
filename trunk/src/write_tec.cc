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

void write_tec(string fileName, double time) {
	fstream file;

	if (rank==0) file.open((fileName).c_str(), ios::out); else file.open((fileName).c_str(), ios::out | ios::in );

	file << "VARIABLES = \"x\", \"y\", \"z\",\"rho\",\"u\",\"v\",\"w\",\"p\" " << endl;
	//file << "VARIABLES = \"x\", \"y\", \"z\" " << endl;
	file << "ZONE n="<<grid.globalNodeCount<< "  e=" <<grid.globalCellCount<< "  VARLOCATION=([4-8]=CellCentered)" << endl;
	if (grid.cell[0].type==TETRA_4) {
		file << "DATAPACKING=BLOCK, ZONETYPE=FETETRAHEDRON, SOLUTIONTIME=" << time << endl;
	} else if (grid.cell[0].type==HEXA_8) {
		file << "DATAPACKING=BLOCK, ZONETYPE=FEBRICK, SOLUTIONTIME=" << time << endl;
	} else if (grid.cell[0].type==PENTA_6) {
		file << "DATAPACKING=BLOCK, ZONETYPE=FEBRICK, SOLUTIONTIME=" << time << endl;
	}

	int header_offset=file.tellg();
	
	if (rank==0) {
		int count=1;
		for (unsigned int i=0;i<3;++i) {
			for (unsigned int n=0;n<grid.globalNodeCount;++n) {
				file << setw(16) << setprecision(8) << scientific << 0.123 ; // just a filler
				if (count==100) {
					file << "\n";
					count=1;
				} else {
					file << "\t";
					++count;
				}
			}
		}
		file << "\n";
		
		count=1;
		for (unsigned int i=0;i<5;++i) {
			for (unsigned int c=0;c<grid.globalCellCount;++c) {
				file << setw(16) << setprecision(8) << scientific << 0.123 ; // just a filler
				if (count==100) {
					file << "\n";
					count=1;
				} else {
					file << "\t";
					++count;
				}
			}
		file << "\n";
		}

	}

	MPI_Barrier(MPI_COMM_WORLD);
	
	int pos=0;
	for (unsigned int i=0;i<3;++i) {
		for (unsigned int n=0;n<grid.nodeCount;++n) {
			pos=header_offset+(i*grid.globalNodeCount*17)+(grid.node[n].globalId*17);
			file.seekg(pos);
			file << setw(16) << setprecision(8) << scientific << grid.node[n].comp[i] ;
		}
	}
	
	for (unsigned int c=0;c<grid.cellCount;++c) {
		pos=header_offset+(3*grid.globalNodeCount*17)+1+(grid.cell[c].globalId*17);
		file.seekg(pos);
		file << setw(16) << setprecision(8) << scientific << grid.cell[c].rho ;
	}
	
	for (unsigned int i=0;i<3;++i) {
		for (unsigned int c=0;c<grid.cellCount;++c) {
			pos=header_offset+(3*grid.globalNodeCount*17)+1+((i+1)*grid.globalCellCount*17)+(grid.cell[c].globalId*17)+i+1;
			file.seekg(pos);
			file << setw(16) << setprecision(8) << scientific << grid.cell[c].v.comp[i] ;
		}
	}
	
	for (unsigned int c=0;c<grid.cellCount;++c) {
		pos=header_offset+(3*grid.globalNodeCount*17)+1+(4*grid.globalCellCount*17)+(grid.cell[c].globalId*17)+4;
		file.seekg(pos);
		file << setw(16) << setprecision(8) << scientific << grid.cell[c].p ;
	}
	
	MPI_Barrier(MPI_COMM_WORLD);
	
	
	if (rank==np-1) {
		file.seekg(0,ios_base::end);
		fstream connFile;
		connFile.open("connectivity.dat", ios::in);
		char line[100];
		for (unsigned int c=0;c<grid.globalCellCount;++c) {
			connFile.getline(line,100);
			file << line << "\n" ;
		}
		connFile.close();
	}
	
/*
	count=1;
	for (unsigned int c=0;c<grid.cellCount;++c) {
		file << setprecision(8) << scientific << grid.cell[c].rho;
		if (count==100) {
			file << "\n";
			count=1;
		} else {
			file << "\t";
			++count;
		}
	}
	file << "\n";
	count =1;
	for (unsigned int i=0;i<3;++i) {
		for (unsigned int c=0;c<grid.cellCount;++c) {
			file << setprecision(8) << scientific << grid.cell[c].v.comp[i] ;
			if (count==100) {
				file << "\n";
				count=1;
			} else {
				file << "\t";
				++count;
			}
		}
		file << "\n";
	}
	count=1;
	for (unsigned int c=0;c<grid.cellCount;++c) {
		file << setprecision(8) << scientific << grid.cell[c].p ;
		if (count==100) {
			file << "\n";
			count=1;
		} else {
			file << "\t";
			++count;
		}
	}
	*/
	/*
	file << "\n";
	for (unsigned int c=0;c<grid.cellCount;++c) {
		if (grid.cell[0].type==PENTA_6) {
			file << grid.cell[c].nodes[0]+1 << "\t" ;
			file << grid.cell[c].nodes[1]+1 << "\t" ;
			file << grid.cell[c].nodes[2]+1 << "\t" ;
			file << grid.cell[c].nodes[2]+1 << "\t" ;
			file << grid.cell[c].nodes[3]+1 << "\t" ;
			file << grid.cell[c].nodes[4]+1 << "\t" ;
			file << grid.cell[c].nodes[5]+1 << "\t" ;
			file << grid.cell[c].nodes[5]+1 << "\t" ;
		} else {
			for (unsigned int n=0;n<grid.cell[c].nodes.size();++n) {
				file << grid.cell[c].nodes[n]+1 << "\t";
			}
		}
		file << endl;
	}
*/
	file.close();

}
