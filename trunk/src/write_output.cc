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

#include <iomanip>
#include <cmath>
#include <mpi.h>
#include <sys/stat.h>
#include <sys/types.h>
using namespace std;

#include "bc.h"
#include "inputs.h"
#include "rans.h"
#include "flamelet.h"

extern BC bc;
extern RANS rans;
extern Flamelet flamelet;

extern string int2str(int number) ;
void write_tec(double time);
void write_vtk(void);
void write_vtk_parallel(void);

void write_output(double time, InputFile input) {
	mkdir("./output",S_IRWXU);
	if (OUTPUT_FORMAT==VTK) {
		// Write vtk output file
		if (Rank==0) write_vtk_parallel();
		write_vtk();
	} else if(OUTPUT_FORMAT==TECPLOT) {
		// Write tecplot output file
		for (int p=0;p<np;++p) {
			if(Rank==p) write_tec(time);
			MPI_Barrier(MPI_COMM_WORLD);
		}
	}
	return;
}

void write_tec(double time) {
	
	ofstream file;
	string fileName="./output/out"+int2str(timeStep)+".dat";
	
	// Shorthand for turbulence model test
	bool turb=(TURBULENCE_MODEL!=NONE);
	
	// Proc 0 creates the output file and writes variable list
	
	if (Rank==0) {
		file.open((fileName).c_str(),ios::out); 
		file << "VARIABLES = \"x\", \"y\", \"z\",\"p\",\"u\",\"v\",\"w\",\"T\",\"rho\",\"Ma\"";
		if (turb) file << ",\"k\",\"omega\" ";
		if (FLAMELET) file << ",\"Mixture Fraction\",\"Mixture Fraction Variance\" ";
		file << endl;
	} else {
		file.open((fileName).c_str(),ios::app);
	}

	// Write header (each proc has its own zone)
	file << "ZONE, T=\"Partition " << Rank << "\"" ;
	file << ", N=" << grid.nodeCount << ", E=" << grid.cellCount << endl;
	file << "DATAPACKING=POINT, ZONETYPE=FEBRICK, SOLUTIONTIME=" << time << endl;

	int zoneNeighborCount=0;
	 for (unsigned int f=0;f<grid.faceCount; ++f) {                    
                if (grid.face[f].bc==GHOST) zoneNeighborCount++; 
        }

	//file << "FACENEIGHBORMODE=GLOBALONETOONE, FACENEIGHBORCONNECTIONS=" << zoneNeighborCount << endl;

	// Write variables
	map<int,double>::iterator it;
	set<int>::iterator sit;
	double p_node,T_node,rho_node,k_node,omega_node,Z_node,Zvar_node;
    	Vec3D v_node;
	int count_p,count_v,count_T,count_k,count_omega,count_rho,count_Z,count_Zvar;
	double Ma;
	for (unsigned int n=0;n<grid.nodeCount;++n) {
		p_node=0.;v_node=0.;T_node=0.;k_node=0.;omega_node=0.;Z_node=0.;Zvar_node=0.;
		for ( it=grid.node[n].average.begin() ; it != grid.node[n].average.end(); it++ ) {
			if ((*it).first>=0) { // if contribution is coming from a real cell
				p_node+=(*it).second*grid.cell[(*it).first].p;
				v_node+=(*it).second*grid.cell[(*it).first].v;
				T_node+=(*it).second*grid.cell[(*it).first].T;
				if (turb) k_node+=(*it).second*rans.cell[(*it).first].k;
				if (turb) omega_node+=(*it).second*rans.cell[(*it).first].omega;
				if (FLAMELET) Z_node+=(*it).second*flamelet.cell[(*it).first].Z;
				if (FLAMELET) Zvar_node+=(*it).second*flamelet.cell[(*it).first].Zvar;
			} else { // if contribution is coming from a ghost cell
				p_node+=(*it).second*grid.ghost[-1*((*it).first+1)].p;
				v_node+=(*it).second*grid.ghost[-1*((*it).first+1)].v;
				T_node+=(*it).second*grid.ghost[-1*((*it).first+1)].T;
				if (turb) k_node+=(*it).second*rans.ghost[-1*((*it).first+1)].k;
				if (turb) omega_node+=(*it).second*rans.ghost[-1*((*it).first+1)].omega;
				if (FLAMELET) Z_node+=(*it).second*flamelet.ghost[-1*((*it).first+1)].Z;
				if (FLAMELET) Zvar_node+=(*it).second*flamelet.ghost[-1*((*it).first+1)].Zvar;
			}
		}
		rho_node=eos.rho(p_node,T_node);
		k_node=max(k_node,kLowLimit);
		omega_node=max(omega_node,omegaLowLimit);
		Z_node=max(Z_node,0.);
		Z_node=min(Z_node,1.);
		Zvar_node=max(Zvar_node,0.);
		
		count_p=0; count_v=0; count_T=0; count_rho=0.; count_k=0; count_omega=0; count_Z=0; count_Zvar=0;
		for (sit=grid.node[n].bcs.begin();sit!=grid.node[n].bcs.end();sit++) {
			if (bc.region[(*sit)].specified==BC_RHO) {
				rho_node=bc.region[(*sit)].rho; count_rho++;
				T_node=eos.T(p_node,rho_node); count_T++;
			} else if (bc.region[(*sit)].specified==BC_T){
				T_node=bc.region[(*sit)].T; count_T++;
				rho_node=eos.rho(p_node,T_node); count_rho++;
			}
			if (bc.region[(*sit)].type==INLET) {
				v_node=bc.region[(*sit)].v; count_v++;
				k_node=bc.region[(*sit)].k; count_k++;
				omega_node=bc.region[(*sit)].omega; count_omega++;
			}
			if (bc.region[(*sit)].type==NOSLIP) {
				if (count_v==0) v_node=0.; count_v++;
				if (count_k==0) k_node=0.; count_k++;
			}

			if (count_rho>0) rho_node/=double(count_rho);
			if (count_T>0) T_node/=double(count_T);
			if (count_v>0) v_node/=double(count_v);
			if (count_k>0) k_node/=double(count_k);
			if (count_omega>0) omega_node/=double(count_omega);
			if (count_Z>0) Z_node/=double(count_Z);
			if (count_Zvar>0) Zvar_node/=double(count_Zvar);
		}
		
		Ma=sqrt((v_node.dot(v_node))/(Gamma*(p_node+Pref)/rho_node));
					
		file << setw(16) << setprecision(8) << scientific;
		file << grid.node[n].comp[0] << "\t";
		file << grid.node[n].comp[1] << "\t";
		file << grid.node[n].comp[2] << "\t";
		file << p_node << "\t" ;
		file << v_node.comp[0] << "\t";
		file << v_node.comp[1] << "\t";
		file << v_node.comp[2] << "\t";
		file << T_node << "\t";
		file << rho_node << "\t";
		file << Ma << "\t";
		if (turb) file << k_node << "\t";
		if (turb) file << omega_node << "\t";
		if (FLAMELET) file << Z_node << "\t";
		if (FLAMELET) file << Zvar_node << "\t";
		file << endl;
	}

	// Write connectivity
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

	// Write zone connectivity
	for (unsigned int f=0;f<grid.faceCount; ++f) {
		if (grid.face[f].bc==GHOST) {
			int g=-1*grid.face[f].neighbor-1;
	//		file << grid.face[f].parent+1 << "\t" << grid.face[f].parentIndex+1 << "\t" << grid.ghost[g].partition+1 << "\t" << grid.ghost[g].matrix_id-grid.partitionOffset[grid.ghost[g].partition]+1 << endl;
		}
	}
	file.close();

}

void write_vtk(void) {

	string filePath="./output/"+int2str(timeStep);
	string fileName=filePath+"/proc"+int2str(Rank)+".vtu";
		
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
		file << offset << endl;
	}
	file << "</DataArray>" << endl;
			
	file << "<DataArray Name=\"types\" type=\"UInt8\" format=\"ascii\" >" << endl;
	for (unsigned int c=0;c<grid.cellCount;++c) {
		if (grid.cell[c].nodeCount==4) file << "10" << endl; // Tetra
		if (grid.cell[c].nodeCount==8) file << "12" << endl; // Hexa
		if (grid.cell[c].nodeCount==6) file << "13" << endl; // Prism (Wedge)
		if (grid.cell[c].nodeCount==5) file << "14" << endl; // Pyramid (Wedge)
	}
	file << endl;
	file << "</DataArray>" << endl;;
	
	file << "</Cells>" << endl;

	file << "<CellData Scalars=\"Pressure\" Vectors=\"Velocity\" format=\"ascii\">" << endl;
	
	file << "<DataArray Name=\"Pressure\" type=\"Float32\" format=\"ascii\" >" << endl;
	for (unsigned int c=0;c<grid.cellCount;++c) file << setw(16) << setprecision(8) << scientific << grid.cell[c].p << endl;
	file << "</DataArray>" << endl;
	
	file << "<DataArray Name=\"Velocity\" NumberOfComponents=\"3\" type=\"Float32\" format=\"ascii\" >" << endl;
	for (unsigned int c=0;c<grid.cellCount;++c) for (unsigned int i=0;i<3;++i) file << setw(16) << setprecision(8) << scientific << grid.cell[c].v[i] << endl;
	file << "</DataArray>" << endl;	
	
	file << "<DataArray Name=\"Temperature\" type=\"Float32\" format=\"ascii\" >" << endl;
	for (unsigned int c=0;c<grid.cellCount;++c) for (unsigned int i=0;i<3;++i) file << setw(16) << setprecision(8) << scientific << grid.cell[c].T << endl;
	file << "</DataArray>" << endl;
	
	file << "<DataArray Name=\"Density\" type=\"Float32\" format=\"ascii\" >" << endl;
	for (unsigned int c=0;c<grid.cellCount;++c) file << setw(16) << setprecision(8) << scientific << grid.cell[c].rho << endl;
	file << "</DataArray>" << endl;


	file << "</CellData>" << endl;
	
	file << "</Piece>" << endl;
	file << "</UnstructuredGrid>" << endl;
	file << "</VTKFile>" << endl;
	file.close();
	
}

void write_vtk_parallel(void) {

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

	file << "<PCellData Scalars=\"Pressure\" Vectors=\"Velocity\" format=\"ascii\">" << endl;
	file << "<DataArray Name=\"Pressure\" type=\"Float32\" format=\"ascii\" />" << endl;
	file << "<DataArray Name=\"Velocity\" NumberOfComponents=\"3\" type=\"Float32\" format=\"ascii\" />" << endl;
	file << "<DataArray Name=\"Temperature\" type=\"Float32\" format=\"ascii\" />" << endl;
	file << "<DataArray Name=\"Density\" type=\"Float32\" format=\"ascii\" />" << endl;
	file << "</PCellData>" << endl;
	for (int p=0;p<np;++p) file << "<Piece Source=\"proc" << int2str(p) << ".vtu\" />" << endl;
	file << "</PUnstructuredGrid>" << endl;
	file << "</VTKFile>" << endl;
	file.close();
	
}

