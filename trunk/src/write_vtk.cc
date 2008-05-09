#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <mpi.h>
#include <sys/stat.h>
#include <sys/types.h>
using namespace std;

#include "grid.h"
#include <cgnslib.h>

extern Grid grid;
extern int np, rank;
extern string int2str(int number) ;

void write_vtk(int timeStep) {

	string filePath="./output/"+int2str(timeStep);
	string fileName=filePath+"/proc"+int2str(rank)+".vtu";
		
	mkdir(filePath.c_str(),S_IRWXU);
	
 	ofstream file;
	file.open((fileName).c_str(),ios::out);
	file << "<?xml version=\"1.0\"?>" << endl;
	file << "<VTKFile type=\"UnstructuredGrid\">" << endl;
	file << "<UnstructuredGrid>" << endl;
	file << "<Piece NumberOfPoints=\"" << grid.nodeCount << "\" NumberOfCells=\"" << grid.cellCount << "\">" << endl;
	file << "<Points>" << endl;
	file << "<DataArray NumberOfComponents=\"3\" type=\"Float32\" format=\"ascii\" >" << endl;
 	for (unsigned int n=0;n<grid.nodeCount;++n) {
		for (unsigned int i=0; i<3; ++i) file<< setw(16) << setprecision(8) << scientific << grid.node[n].comp[i] << endl;
	}
	file << "</DataArray>" << endl;
	file << "</Points>" << endl;
	file << "<Cells>" << endl;
	
	file << "<DataArray Name=\"connectivity\" type=\"Int32\" format=\"ascii\" >" << endl;
	for (unsigned int c=0;c<grid.cellCount;++c) {
		for (unsigned int n=0;n<grid.cell[c].nodeCount;++n) {
			file << grid.cell[c].nodes[n] << "\t";
		}
		file << endl;
	}
	
	file << "</DataArray>" << endl;
	file << "<DataArray Name=\"offsets\" type=\"Int32\" format=\"ascii\" >" << endl;
	int offset=0;
	for (unsigned int c=0;c<grid.cellCount;++c) {
		offset+=grid.cell[c].nodeCount;
		file << offset << "\t";
	}
	file << endl;
	file << "</DataArray>" << endl;
			
	file << "<DataArray Name=\"types\" type=\"UInt8\" format=\"ascii\" >" << endl;
	for (unsigned int c=0;c<grid.cellCount;++c) {
		if (grid.cell[c].type==TETRA_4) file << "10\t";
		if (grid.cell[c].type==HEXA_8) file << "12\t";
		if (grid.cell[c].type==PENTA_6) file << "13\t";
	}
	file << endl;
	file << "</DataArray>" << endl;;
	
	file << "</Cells>" << endl;

	file << "<CellData Scalars=\"Density\" Vectors=\"Velocity\" format=\"ascii\">" << endl;
	file << "<DataArray Name=\"Density\" type=\"Float32\" format=\"ascii\" >" << endl;
	for (unsigned int c=0;c<grid.cellCount;++c) file << setw(16) << setprecision(8) << scientific << grid.cell[c].rho << endl;
	file << "</DataArray>" << endl;
	file << "<DataArray Name=\"Velocity\" NumberOfComponents=\"3\" type=\"Float32\" format=\"ascii\" >" << endl;
	for (unsigned int c=0;c<grid.cellCount;++c) for (unsigned int i=0;i<3;++i) file << setw(16) << setprecision(8) << scientific << grid.cell[c].v.comp[i] << endl;
	file << "</DataArray>" << endl;
	file << "<DataArray Name=\"Pressure\" type=\"Float32\" format=\"ascii\" >" << endl;
	for (unsigned int c=0;c<grid.cellCount;++c) file << setw(16) << setprecision(8) << scientific << grid.cell[c].p << endl;
	file << "</DataArray>" << endl;
	file << "</CellData>" << endl;
	
	file << "</Piece>" << endl;
	file << "</UnstructuredGrid>" << endl;
	file << "</VTKFile>" << endl;
	file.close();
	
}

void write_vtk_parallel(int timeStep) {

	string filePath="./output/"+int2str(timeStep);
	string fileName=filePath+"/out"+int2str(timeStep)+".pvtu";
		
	mkdir("./output",S_IRWXU);
	mkdir(filePath.c_str(),S_IRWXU);
	
	ofstream file;
	file.open((fileName).c_str(),ios::out);
	file << "<?xml version=\"1.0\"?>" << endl;
	file << "<VTKFile type=\"PUnstructuredGrid\">" << endl;
	file << "<PUnstructuredGrid GhostLevel=\"0\">" << endl;
	file << "<PPoints>" << endl;
	file << "<DataArray NumberOfComponents=\"3\" type=\"Float32\" format=\"ascii\" />" << endl;
	file << "</PPoints>" << endl;

	file << "<PCellData Scalars=\"Density\" Vectors=\"Velocity\" format=\"ascii\">" << endl;
	file << "<DataArray Name=\"Density\" type=\"Float32\" format=\"ascii\" />" << endl;
	file << "<DataArray Name=\"Velocity\" NumberOfComponents=\"3\" type=\"Float32\" format=\"ascii\" />" << endl;
	file << "<DataArray Name=\"Pressure\" type=\"Float32\" format=\"ascii\" />" << endl;
	file << "</PCellData>" << endl;
	for (int p=0;p<np;++p) file << "<Piece Source=\"proc" << int2str(p) << ".vtu\" />" << endl;
	file << "</PUnstructuredGrid>" << endl;
	file << "</VTKFile>" << endl;
	file.close();
	
}
