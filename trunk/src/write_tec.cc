#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
using namespace std;

#include "grid.h"
#include <cgnslib.h>

extern Grid grid;

void write_tec(string fileName, double time) {
fstream file;

file.open((fileName).c_str(), ios::out);

file << "VARIABLES = \"x\", \"y\", \"z\",\"rho\",\"u\",\"v\",\"w\",\"p\" " << endl;
file << "ZONE n="<<grid.nodeCount<< "  e=" <<grid.cellCount<< "  VARLOCATION=([4-8]=CellCentered)" << endl;
if (grid.cell[0].type==TETRA_4) {
    file << "DATAPACKING=BLOCK, ZONETYPE=FETETRAHEDRON, SOLUTIONTIME=" << time << endl;
} else if (grid.cell[0].type==HEXA_8) {
    file << "DATAPACKING=BLOCK, ZONETYPE=FEBRICK, SOLUTIONTIME=" << time << endl;
} else if (grid.cell[0].type==PENTA_6) {
    file << "DATAPACKING=BLOCK, ZONETYPE=FEBRICK, SOLUTIONTIME=" << time << endl;
}

int count=1;
for (unsigned int i=0;i<3;++i) {
	for (unsigned int n=0;n<grid.nodeCount;++n) {
		file << setprecision(8) << scientific << grid.node[n].comp[i] ;
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
file << "\n";
for (unsigned int c=0;c<grid.cellCount;++c) {
  if  (grid.cell[0].type==PENTA_6) {
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

file.close();

}
