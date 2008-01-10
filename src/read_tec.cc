#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
using namespace std;

#include "grid.h"
#include <cgnslib.h>

extern Grid grid;

void read_tec(string fileName, double &time) {
	fstream file;
	string data="";

	file.open((fileName).c_str(), ios::in);

	for (unsigned int i=0; i<8; ++i) getline(file,data,'=');

	file >> time; getline(file,data,'\n');

	int lineCount;

// Skip the node data
	lineCount=grid.nodeCount*3/int (100);
	if ((grid.nodeCount*3) %int (100) !=0) lineCount++;
	for (unsigned int n=0;n<lineCount;++n) getline(file,data,'\n');

	for (unsigned int c=0;c<grid.cellCount;++c) file >> grid.cell[c].rho;
	for (unsigned int i=0;i<3;++i) for (unsigned int c=0;c<grid.cellCount;++c) file >> grid.cell[c].v.comp[i];
	for (unsigned int c=0;c<grid.cellCount;++c) file >> grid.cell[c].p;

	file.close();

}
