#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
using namespace std;

#include "grid.h"
#include <cgnslib.h>

extern Grid grid;

void read_tec(string fileName, int global2local[], double &time) {
	fstream file;
	string data="";

	file.open((fileName).c_str(), ios::in);

	for (unsigned int i=0; i<8; ++i) getline(file,data,'=');

	file >> time; getline(file,data,'\n');

	int lineCount;

	// Skip the node data
	lineCount=grid.globalNodeCount*3/int (100);
	if ((grid.nodeCount*3) %int (100) !=0) lineCount++;
	for (unsigned int n=0;n<lineCount;++n) getline(file,data,'\n');

	double temp;
	for (unsigned int c=0;c<grid.globalCellCount;++c) {
		file >> temp;
		if (global2local[c]>=0) grid.cell[global2local[c]].rho=temp;
	}
	
	for (unsigned int i=0;i<3;++i) {
		for (unsigned int c=0;c<grid.globalCellCount;++c) {
			file >> temp;
			if (global2local[c]>=0) grid.cell[global2local[c]].v.comp[i]=temp;
		}
	}
	
	for (unsigned int c=0;c<grid.globalCellCount;++c) {
		file >> temp;
		if (global2local[c]>=0) grid.cell[global2local[c]].p=temp;
	}
	
	file.close();

}
