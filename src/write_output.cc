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
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <cmath>
#include <mpi.h>
#include <sys/stat.h>
#include <sys/types.h>
using namespace std;

#include "grid.h"
#include "bc.h"
#include "inputs.h"
#include <cgnslib.h>

extern Grid grid;
extern BC bc;
extern int np, rank;
extern double Gamma;
extern double Pref;
extern string int2str(int number) ;
void write_tec(int timeStep, double time);
void write_vtk(int timeStep);
void write_vtk_parallel(int timeStep);

void write_output(int timeStep, double time, InputFile input) {
	mkdir("./output",S_IRWXU);
	if (input.section["output"].strings["format"]=="vtk") {
		// Write vtk output file
		if (rank==0) write_vtk_parallel(timeStep);
		write_vtk(timeStep);
	} else if(input.section["output"].strings["format"]=="tecplot") {
		// Write tecplot output file
		for (int p=0;p<np;++p) {
			if(rank==p) write_tec(timeStep,time);
			MPI_Barrier(MPI_COMM_WORLD);
		}
	}
	return;
}

void write_tec(int timeStep, double time) {

	ofstream file;
	string fileName="./output/out"+int2str(timeStep)+".dat";
	// Proc 0 creates the output file and writes variable list
	if (rank==0) {
		file.open((fileName).c_str(),ios::out); 
		file << "VARIABLES = \"x\", \"y\", \"z\",\"rho\",\"u\",\"v\",\"w\",\"p\",\"Ma\" " << endl;
	} else {
		file.open((fileName).c_str(),ios::app);
	}

	// Write header (each proc has its own zone)
	file << "ZONE, T=\"Partition " << rank << "\"" ;
	file << ", N=" << grid.nodeCount << ", E=" << grid.cellCount << endl;
	file << "DATAPACKING=POINT, ZONETYPE=FEBRICK, SOLUTIONTIME=" << time << endl;

	// Write variables
	map<int,double>::iterator it;
	set<int>::iterator sit;
	double rho_node,p_node;
    Vec3D v_node;
	int count_rho,count_v,count_p;
	double Ma;
	for (unsigned int n=0;n<grid.nodeCount;++n) {
		rho_node=0.;v_node=0.;p_node=0.;
		for ( it=grid.node[n].average.begin() ; it != grid.node[n].average.end(); it++ ) {
			if ((*it).first>=0) { // if contribution is coming from a real cell
				rho_node+=(*it).second*grid.cell[(*it).first].rho;
				v_node+=(*it).second*grid.cell[(*it).first].v;
				p_node+=(*it).second*grid.cell[(*it).first].p;
			} else { // if contribution is coming from a ghost cell
				rho_node+=(*it).second*grid.ghost[-1*((*it).first+1)].rho;
				v_node+=(*it).second*grid.ghost[-1*((*it).first+1)].v;
				p_node+=(*it).second*grid.ghost[-1*((*it).first+1)].p;
			}
		}
		
		count_rho=0; count_v=0; count_p=0;
		for (sit=grid.node[n].bcs.begin();sit!=grid.node[n].bcs.end();sit++) {
			if (bc.region[(*sit)].type=="inlet") {
				rho_node=bc.region[(*sit)].rho;
				v_node=bc.region[(*sit)].v;
				count_rho++; count_v++;
			}
			if (bc.region[(*sit)].type=="noslip") {
				if (count_v==0) v_node=0.;
				count_v++;
			}

			if (count_rho>0) rho_node/=double(count_rho);
			if (count_v>0) v_node/=double(count_v);
			if (count_p>0) p_node/=double(count_p);
		}
		
		Ma=sqrt((v_node.dot(v_node))/(Gamma*(p_node+Pref)/rho_node));
					
		file << setw(16) << setprecision(8) << scientific;
		file << grid.node[n].comp[0] << "\t";
		file << grid.node[n].comp[1] << "\t";
		file << grid.node[n].comp[2] << "\t";
		file << rho_node << "\t" ;
		file << v_node.comp[0] << "\t";
		file << v_node.comp[1] << "\t";
		file << v_node.comp[2] << "\t";
		file << p_node << "\t";
		file << Ma << "\t";
		file << endl;
	}

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

