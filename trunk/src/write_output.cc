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

#include <cmath>
#include <mpi.h>

using namespace std;

#include "bc.h"
#include "inputs.h"
#include "rans.h"
#include "flamelet.h"

extern BC bc;
extern RANS rans;
extern Flamelet flamelet;

extern string int2str(int number) ;
void write_tec_vars(int &nVar);
void write_tec_cells(void);
void write_tec_bcs(int bcNo,int nVar);
void write_vtk(void);
void write_vtk_parallel(void);

void write_output(double time, InputFile input) {
	mkdir("./output",S_IRWXU);
	if (OUTPUT_FORMAT==VTK) {
		// Write vtk output file
		if (Rank==0) write_vtk_parallel();
		write_vtk();
	} else if(OUTPUT_FORMAT==TECPLOT) {
		int nVar=0;
		// Write tecplot output file
		for (int p=0;p<np;++p) {
			if(Rank==p) write_tec_vars(nVar);
			MPI_Barrier(MPI_COMM_WORLD);
		}
		for (int p=0;p<np;++p) {
			if(Rank==p) write_tec_cells();
			MPI_Barrier(MPI_COMM_WORLD);
		}
		for (int nbc=0;nbc<grid.globalBoundaryFaceCount.size();++nbc) {
			for (int p=0;p<np;++p) {
				if(Rank==p) write_tec_bcs(nbc,nVar);
				MPI_Barrier(MPI_COMM_WORLD);
			}
		}
	}
	return;
}

void write_tec_vars(int &nVar) {
	
	ofstream file;
	string fileName="./output/out"+int2str(timeStep)+".dat";
	
	// Shorthand for turbulence model test
	bool turb=(TURBULENCE_MODEL!=NONE);
	
	// Proc 0 creates the output file and writes variable list
	
	if (Rank==0) {
		file.open((fileName).c_str(),ios::out); 
		nVar=11;
		file << "VARIABLES = \"x\", \"y\", \"z\",\"p\",\"u\",\"v\",\"w\",\"T\",\"rho\",\"Ma\",\"dt\"";
		if (turb) {
			file << ",\"k\",\"omega\",\"mu_t\"";
			nVar+=3;
			if (TURBULENCE_FILTER!=NONE) {
				file << ",\"filterFunction\"";
				nVar++;
			}
		}
		if (FLAMELET) {
			file << ",\"viscosity\", \"Mixture Fraction\",\"Mixture Fraction Variance\" ";
			nVar+=3;
		}
		file << endl;
		file << "ZONE, T=\"Volume " << Rank << "\" ZONETYPE=FEBRICK DATAPACKING=POINT" << endl;
		file << "NODES=" << grid.globalNodeCount << " ELEMENTS=" << grid.globalCellCount << endl;
	} else {
		file.open((fileName).c_str(),ios::app);
	}
	
		// Write variables
		map<int,double>::iterator it;
		set<int>::iterator sit;
		double p_node,T_node,rho_node,k_node,omega_node,Z_node,Zvar_node,visc_node,dt_node;
		double mu_t_node, filter_node,R_node;
		Vec3D v_node;
		int count_p,count_v,count_T,count_k,count_omega,count_mu_t,count_rho,count_Z,count_Zvar;
		double Ma;
		
		for (unsigned int n=0;n<grid.nodeCount;++n) {
			if (maps.nodeGlobal2Output[grid.node[n].globalId]>=grid.nodeCountOffset) {
				
				p_node=0.;v_node=0.;T_node=0.;k_node=0.;omega_node=0.;Z_node=0.;Zvar_node=0.;
				mu_t_node=0.; filter_node=0.; dt_node=0.; R_node=0.;
				double real_weight_sum=0.;
				for ( it=grid.node[n].average.begin() ; it != grid.node[n].average.end(); it++ ) {
					if ((*it).first>=0) { // if contribution is coming from a real cell
						p_node+=(*it).second*grid.cell[(*it).first].p;
						v_node+=(*it).second*grid.cell[(*it).first].v;
						T_node+=(*it).second*grid.cell[(*it).first].T;
						dt_node+=(*it).second*grid.cell[(*it).first].dt;
						if (turb) {
							k_node+=(*it).second*rans.cell[(*it).first].k;
							omega_node+=(*it).second*rans.cell[(*it).first].omega;
							mu_t_node+=(*it).second*rans.cell[(*it).first].mu_t;
							filter_node+=(*it).second*rans.cell[(*it).first].filterFunction;
						}
						if (FLAMELET) {
							Z_node+=(*it).second*flamelet.cell[(*it).first].Z;
							Zvar_node+=(*it).second*flamelet.cell[(*it).first].Zvar;
							R_node+=(*it).second*flamelet.cell[(*it).first].R;
						}
						real_weight_sum+=(*it).second;
					} else { // if contribution is coming from a ghost cell
						p_node+=(*it).second*grid.ghost[-1*((*it).first+1)].p;
						v_node+=(*it).second*grid.ghost[-1*((*it).first+1)].v;
						T_node+=(*it).second*grid.ghost[-1*((*it).first+1)].T;
						if (turb) {
							k_node+=(*it).second*rans.ghost[-1*((*it).first+1)].k;
							omega_node+=(*it).second*rans.ghost[-1*((*it).first+1)].omega;
							mu_t_node+=(*it).second*rans.ghost[-1*((*it).first+1)].mu_t;
						}
						if (FLAMELET) {
							Z_node+=(*it).second*flamelet.ghost[-1*((*it).first+1)].Z;
							Zvar_node+=(*it).second*flamelet.ghost[-1*((*it).first+1)].Zvar;
							R_node+=(*it).second*flamelet.ghost[-1*((*it).first+1)].R;
						}
					}
				}
				
				rho_node=eos.rho(p_node,T_node);
				
				if (real_weight_sum>=1.e-10) {
					filter_node/=real_weight_sum;
					dt_node/=real_weight_sum;
				} else { filter_node=1.; dt_node=dt_current; }
				
				filter_node=min(1.,filter_node);
				filter_node=max(0.,filter_node);
				
				k_node=max(k_node,kLowLimit);
				omega_node=max(omega_node,omegaLowLimit);
				Z_node=max(Z_node,0.);
				Z_node=min(Z_node,1.);
				Zvar_node=max(Zvar_node,0.);
				
				if (FLAMELET) {
					double Chi=2.0*rans.kepsilon.beta_star*omega_node*Zvar_node;
					flamelet.table.get_rho_T_comp(p_node,Z_node,Zvar_node,Chi,rho_node,T_node);
					visc_node=flamelet.table.get_mu(Z_node,Zvar_node,Chi,false);
				}
				
				count_p=0; count_v=0; count_T=0; count_rho=0; count_k=0; count_omega=0; count_Z=0; count_Zvar=0; count_mu_t=0;
				for (sit=grid.node[n].bcs.begin();sit!=grid.node[n].bcs.end();sit++) {
					if (bc.region[(*sit)].specified==BC_STATE) {
						if (count_p>0) p_node+=bc.region[(*sit)].p; else p_node+=bc.region[(*sit)].p;
						if (count_rho>0) rho_node+=bc.region[(*sit)].rho; else rho_node=bc.region[(*sit)].rho;
						if (count_T>0) T_node+=bc.region[(*sit)].T; else T_node=bc.region[(*sit)].T;
						count_p++; count_rho++; count_T++;
					} else if (bc.region[(*sit)].specified==BC_RHO) {
						if (count_rho>0) rho_node+=bc.region[(*sit)].rho; else rho_node=bc.region[(*sit)].rho;
						count_rho++;
					} else if (bc.region[(*sit)].specified==BC_T) {
						if (count_T>0) T_node+=bc.region[(*sit)].T; else T_node=bc.region[(*sit)].T;
						count_T++;
					} else if (bc.region[(*sit)].specified==BC_P) {
						if (count_p>0) p_node+=bc.region[(*sit)].p; else p_node+=bc.region[(*sit)].p;
						count_p++;
					}

					if (bc.region[(*sit)].type==INLET) {
						if (bc.region[(*sit)].kind==VELOCITY) {
							if (count_v>0) v_node+=bc.region[(*sit)].v; else v_node=bc.region[(*sit)].v;
							count_v++;
						}
						if (turb) {
							if (count_k>0) k_node+=bc.region[(*sit)].k; else k_node+=bc.region[(*sit)].k;
							if (count_omega>0) omega_node+=bc.region[(*sit)].omega; else omega_node+=bc.region[(*sit)].omega;
							count_k++; count_omega++;
						}
					}
					
					if (bc.region[(*sit)].type==NOSLIP) {
						v_node=0.; count_v++;
						k_node=0.; count_k++;
						mu_t_node=0.; count_mu_t++;
					}
		
					if (count_p>0) p_node/=double(count_p);
					if (count_rho>0) rho_node/=double(count_rho);
					if (count_T>0) T_node/=double(count_T);
					if (count_v>0) v_node/=double(count_v);
					if (count_k>0) k_node/=double(count_k);
					if (count_omega>0) omega_node/=double(count_omega);
					if (count_mu_t>0) mu_t_node/=double(count_mu_t);
					if (count_Z>0) Z_node/=double(count_Z);
					if (count_Zvar>0) Zvar_node/=double(count_Zvar);
				}
				
				Ma=sqrt((v_node.dot(v_node))/(Gamma*(p_node+Pref)/rho_node));
		
				file << setw(16) << setprecision(8) << scientific;
				file << grid.node[n][0] << "\t";
				file << grid.node[n][1] << "\t";
				file << grid.node[n][2] << "\t";
				file << p_node << "\t" ;
				file << v_node.comp[0] << "\t";
				file << v_node.comp[1] << "\t";
				file << v_node.comp[2] << "\t";
				file << T_node << "\t";
				file << rho_node << "\t";
				file << Ma << "\t";
				file << dt_node << "\t";
				if (turb) {
					file << k_node << "\t";
					file << omega_node << "\t";
					file << mu_t_node << "\t";
					if (TURBULENCE_FILTER!=NONE) file << filter_node << "\t";
				}
				if (FLAMELET) {
					file << visc_node << "\t";
					file << Z_node << "\t";
					file << Zvar_node << "\t";
				}
				file << endl;
			}
		}

	file.close();

}

void write_tec_cells() {
	
	ofstream file;
	string fileName="./output/out"+int2str(timeStep)+".dat";
	
	file.open((fileName).c_str(),ios::app);

		// Write connectivity
	for (unsigned int c=0;c<grid.cellCount;++c) {
		if (grid.cell[c].nodeCount==4) {
			file << maps.nodeGlobal2Output[grid.cell[c].node(0).globalId]+1 << "\t" ;
			file << maps.nodeGlobal2Output[grid.cell[c].node(2).globalId]+1 << "\t" ;
			file << maps.nodeGlobal2Output[grid.cell[c].node(1).globalId]+1 << "\t" ;
			file << maps.nodeGlobal2Output[grid.cell[c].node(1).globalId]+1 << "\t" ;
			file << maps.nodeGlobal2Output[grid.cell[c].node(3).globalId]+1 << "\t" ;
			file << maps.nodeGlobal2Output[grid.cell[c].node(3).globalId]+1 << "\t" ;
			file << maps.nodeGlobal2Output[grid.cell[c].node(3).globalId]+1 << "\t" ;
			file << maps.nodeGlobal2Output[grid.cell[c].node(3).globalId]+1 << "\t" ;
		}
		else if (grid.cell[c].nodeCount==5) {
			file << maps.nodeGlobal2Output[grid.cell[c].node(0).globalId]+1 << "\t" ;
			file << maps.nodeGlobal2Output[grid.cell[c].node(1).globalId]+1 << "\t" ;
			file << maps.nodeGlobal2Output[grid.cell[c].node(2).globalId]+1 << "\t" ;
			file << maps.nodeGlobal2Output[grid.cell[c].node(3).globalId]+1 << "\t" ;
			file << maps.nodeGlobal2Output[grid.cell[c].node(4).globalId]+1 << "\t" ;
			file << maps.nodeGlobal2Output[grid.cell[c].node(4).globalId]+1 << "\t" ;
			file << maps.nodeGlobal2Output[grid.cell[c].node(4).globalId]+1 << "\t" ;
			file << maps.nodeGlobal2Output[grid.cell[c].node(4).globalId]+1 << "\t" ;
		}
		else if (grid.cell[c].nodeCount==6) {
			file << maps.nodeGlobal2Output[grid.cell[c].node(0).globalId]+1 << "\t" ;
			file << maps.nodeGlobal2Output[grid.cell[c].node(1).globalId]+1 << "\t" ;
			file << maps.nodeGlobal2Output[grid.cell[c].node(2).globalId]+1 << "\t" ;
			file << maps.nodeGlobal2Output[grid.cell[c].node(2).globalId]+1 << "\t" ;
			file << maps.nodeGlobal2Output[grid.cell[c].node(3).globalId]+1 << "\t" ;
			file << maps.nodeGlobal2Output[grid.cell[c].node(4).globalId]+1 << "\t" ;
			file << maps.nodeGlobal2Output[grid.cell[c].node(5).globalId]+1 << "\t" ;
			file << maps.nodeGlobal2Output[grid.cell[c].node(5).globalId]+1 << "\t" ;
		} else if (grid.cell[c].nodeCount==8) {
			for (unsigned int i=0;i<8;++i) {
				file << maps.nodeGlobal2Output[grid.cell[c].node(i).globalId]+1 << "\t" ;
			}
		}
		file << endl;
	}
	
	file.close();

}

void write_tec_bcs(int bcNo,int nVar) {
	
	ofstream file;
	string fileName="./output/out"+int2str(timeStep)+".dat";
	
	file.open((fileName).c_str(),ios::app);
	if (Rank==0) {
		file << "ZONE T=\"BC_" << bcNo+1 << "\" " << "NODES=" << grid.globalNodeCount << " ELEMENTS=" << grid.globalBoundaryFaceCount[bcNo] << " ZONETYPE=FEQUADRILATERAL" << endl;
		file << "VARSHARELIST = ([1-" << nVar << "]=1)" << endl;
	}

	// Write connectivity
	for (unsigned int f=0;f<grid.faceCount;++f) {
		if (grid.face[f].bc==bcNo) {
			if (grid.face[f].nodeCount==3) {
				file << maps.nodeGlobal2Output[grid.face[f].node(0).globalId]+1 << "\t" ;
				file << maps.nodeGlobal2Output[grid.face[f].node(1).globalId]+1 << "\t" ;
				file << maps.nodeGlobal2Output[grid.face[f].node(2).globalId]+1 << "\t" ;
				file << maps.nodeGlobal2Output[grid.face[f].node(2).globalId]+1 << "\t" ;
			} else if (grid.face[f].nodeCount==4) {
				for (unsigned int i=0;i<4;++i) {
					file << maps.nodeGlobal2Output[grid.face[f].node(i).globalId]+1 << "\t" ;
				}
			}
			file << endl;
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

