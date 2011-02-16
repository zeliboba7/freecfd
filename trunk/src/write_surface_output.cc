/************************************************************************
	
	Copyright 2007-2010 Emre Sozer 
 
	Contact: emresozer@freecfd.com

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
using namespace std;

#include "utilities.h"
#include "vec3d.h"
#include "grid.h"
#include "inputs.h"
#include "ns.h"
#include "rans.h"
#include "hc.h"
#include "commons.h"

extern vector<Grid> grid;
extern InputFile input;
extern vector<NavierStokes> ns;
extern vector<HeatConduction> hc;
extern vector<RANS> rans;
extern vector<Variable<double> > dt;
extern vector<int> equations;
extern vector<Loads> loads;

namespace surface_output {
	int timeStep,gid;
	vector<string> varList;
	vector<bool> var_is_vec3d;
}
using namespace surface_output;

void write_surface_tec_header(int b);
void write_surface_tec_nodes(int i);
void write_surface_tec_var(int ov,int i,int b);
void write_surface_tec_cells(int b);

void write_surface_output(int gridid, int step) {
	mkdir("./output",S_IRWXU);
	gid=gridid;
	timeStep=step;
	varList=input.section("grid",gid).subsection("writeoutput").get_stringList("surfacevariables");
	var_is_vec3d.resize(varList.size());
	int nVar=0;
	for (int var=0; var<varList.size(); ++var) {
		if (varList[var]!="null") nVar++;
		var_is_vec3d[var]=false;
		if (varList[var]=="V" || varList[var]=="tau") {
			var_is_vec3d[var]=true;
			nVar+=2;
		}
	}

	string format=input.section("grid",gid).subsection("writeoutput").get_string("format");

	if (format=="tecplot") {
	
		for (int b=0;b<grid[gid].bcCount;++b) {
			// Write tecplot output file
			if (Rank==0) write_surface_tec_header(b);
			
			if (b==0) {
				for (int i=0;i<3;++i) {
					for (int p=0;p<np;++p) {
						if(Rank==p) write_surface_tec_nodes(i);
						MPI_Barrier(MPI_COMM_WORLD);
					}
				}
			}
			
			for (int ov=0;ov<varList.size();++ov) {
				int nn=1;
				if (var_is_vec3d[ov]) nn=3;
				for (int i=0;i<nn;++i) {
					for (int p=0;p<np;++p) {
						if(Rank==p) write_surface_tec_var(ov,i,b);
						MPI_Barrier(MPI_COMM_WORLD);
					}
					
				}
			}
			
			for (int p=0;p<np;++p) {
				if(Rank==p) write_surface_tec_cells(b);
				MPI_Barrier(MPI_COMM_WORLD);
			}
		}
	} 	
	return;
}

void write_surface_tec_header(int b) {
	
	ofstream file;
	string fileName="./output/surface_"+int2str(timeStep)+"_"+int2str(gid+1)+".dat";
	
	// Proc 0 creates the output file and writes variable list
	int nvars=3;

	if (b==0) {
		file.open((fileName).c_str(),ios::out);
		file << "VARIABLES = \"x\", \"y\", \"z\" ";
	}
	
	for (int var=0;var<varList.size();++var) {
		if (var_is_vec3d[var]) {
			if (b==0) {
				file << ",\"" << varList[var] << "_x\" "; 
				file << ",\"" << varList[var] << "_y\" "; 
				file << ",\"" << varList[var] << "_z\" "; 
			}
			nvars+=3;
		} else if (varList[var]!="null") {
			if (b==0) file << ",\"" << varList[var] << "\" "; 
			nvars++;
		}
	}
	file << endl;
	

	if (b!=0) file.open((fileName).c_str(),ios::app);
	
	file << "ZONE, T=\"BC_" << b+1 << "\", ZONETYPE=FEQUADRILATERAL, DATAPACKING=BLOCK" << endl;
	if (b==0) file << "NODES=" << grid[gid].global_bc_nodeCount << ", ";
	file << "ELEMENTS=" << grid[gid].globalBoundaryFaceCount[b] << endl;
	if (nvars==4) {
		file << "VARLOCATION=([4]=CELLCENTERED)" << endl;	
	} else {
		file << "VARLOCATION=([4-" << nvars << "]=CELLCENTERED)" << endl;
	}
	if (b!=0) file << "VARSHARELIST=([1-3]=1)" << endl;
	return;
}
	
void write_surface_tec_nodes(int i) {
	
	ofstream file;
	string fileName="./output/surface_"+int2str(timeStep)+"_"+int2str(gid+1)+".dat";
	
	file.open((fileName).c_str(),ios::app);
	file << scientific << setprecision(8);
	
	int count=0;
	for (int n=0;n<grid[gid].nodeCount;++n) {
		// Note that some nodes are repeated in different partitions
		if (grid[gid].node[n].bc_output_id>=grid[gid].node_bc_output_offset) {
			// Write node coordinates
			file << grid[gid].node[n][i];
			count++;
			if (count%10==0) file << "\n";
			else file << "\t";
		}
	}
	file << endl;
	file.close();
	
	return;
}
			
			
void write_surface_tec_var(int ov,int i,int b) {

	ofstream file;
	string fileName="./output/surface_"+int2str(timeStep)+"_"+int2str(gid+1)+".dat";
	
	file.open((fileName).c_str(),ios::app);
	file << scientific << setprecision(8);
	
	if (varList[ov]=="p") {
		for (int bf=0;bf<grid[gid].boundaryFaces[b].size();++bf) {
			//int parent=grid[gid].face[grid[gid].boundaryFaces[b][bf]].parent;
			//file << ns[gid].p.cell(parent); // DEBUG
			file << ns[gid].p.face(grid[gid].boundaryFaces[b][bf]);
			if ((bf+1)%10==0) file << "\n";
			else file << "\t";
		}
	} else if (varList[ov]=="T") {
		if (equations[gid]==NS)	{
			for (int bf=0;bf<grid[gid].boundaryFaces[b].size();++bf) {
				file << ns[gid].T.face(grid[gid].boundaryFaces[b][bf]);
				if ((bf+1)%10==0) file << "\n";
				else file << "\t";
			}
		} else if (equations[gid]==HEAT) {
			for (int bf=0;bf<grid[gid].boundaryFaces[b].size();++bf) {
				file << hc[gid].T.face(grid[gid].boundaryFaces[b][bf]);
				if ((bf+1)%10==0) file << "\n";
				else file << "\t";
			}
		}
	} else if (varList[ov]=="rho") {
		for (int bf=0;bf<grid[gid].boundaryFaces[b].size();++bf) {
			file << ns[gid].rho.face(grid[gid].boundaryFaces[b][bf]);
			if ((bf+1)%10==0) file << "\n";
			else file << "\t";
		}
	} else if (varList[ov]=="dt") {
		for (int bf=0;bf<grid[gid].boundaryFaces[b].size();++bf) {
			file << dt[gid].face(grid[gid].boundaryFaces[b][bf]);
			if ((bf+1)%10==0) file << "\n";
			else file << "\t";
		}
	} else if (varList[ov]=="mu") {
		for (int bf=0;bf<grid[gid].boundaryFaces[b].size();++bf) {
			file << ns[gid].material.viscosity(ns[gid].T.face(grid[gid].boundaryFaces[b][bf]));
			if ((bf+1)%10==0) file << "\n";
			else file << "\t";
		}
	} else if (varList[ov]=="lambda") {
		for (int bf=0;bf<grid[gid].boundaryFaces[b].size();++bf) {
			file << ns[gid].material.therm_cond(ns[gid].T.face(grid[gid].boundaryFaces[b][bf]));
			if ((bf+1)%10==0) file << "\n";
			else file << "\t";
		}
	} else if (varList[ov]=="Cp") {
		for (int bf=0;bf<grid[gid].boundaryFaces[b].size();++bf) {
			file << ns[gid].material.Cp(ns[gid].T.face(grid[gid].boundaryFaces[b][bf]));
			if ((bf+1)%10==0) file << "\n";
			else file << "\t";
		}
	} else if (varList[ov]=="resp") {
		for (int bf=0;bf<grid[gid].boundaryFaces[b].size();++bf) {
			file << ns[gid].update[0].face(grid[gid].boundaryFaces[b][bf]);
			if ((bf+1)%10==0) file << "\n";
			else file << "\t";
		}
	} else if (varList[ov]=="resT") {
		for (int bf=0;bf<grid[gid].boundaryFaces[b].size();++bf) {
			file << ns[gid].update[4].face(grid[gid].boundaryFaces[b][bf]);
			if ((bf+1)%10==0) file << "\n";
			else file << "\t";
		}
	} else if (varList[ov]=="limiterp") {
		for (int bf=0;bf<grid[gid].boundaryFaces[b].size();++bf) {
			file << ns[gid].limiter[0].face(grid[gid].boundaryFaces[b][bf]);
			if ((bf+1)%10==0) file << "\n";
			else file << "\t";
		}
	} else if (varList[ov]=="limiterT") {
		for (int bf=0;bf<grid[gid].boundaryFaces[b].size();++bf) {
			file << ns[gid].limiter[4].face(grid[gid].boundaryFaces[b][bf]);
			if ((bf+1)%10==0) file << "\n";
			else file << "\t";
		}
	} else if (varList[ov]=="mdot") {
		for (int bf=0;bf<grid[gid].boundaryFaces[b].size();++bf) {
			file << ns[gid].mdot.face(grid[gid].boundaryFaces[b][bf]);
			if ((bf+1)%10==0) file << "\n";
			else file << "\t";
		}
	} else if (varList[ov]=="qdot") {
		if (equations[gid]==NS)	{
			for (int bf=0;bf<grid[gid].boundaryFaces[b].size();++bf) {
				file << ns[gid].qdot.face(grid[gid].boundaryFaces[b][bf]);
				if ((bf+1)%10==0) file << "\n";
				else file << "\t";
			}
		} else if (equations[gid]==HEAT) {
			for (int bf=0;bf<grid[gid].boundaryFaces[b].size();++bf) {
				file << hc[gid].qdot.face(grid[gid].boundaryFaces[b][bf]);
				if ((bf+1)%10==0) file << "\n";
				else file << "\t";
			}
		}
	} else if (varList[ov]=="yplus") {
		for (int bf=0;bf<grid[gid].boundaryFaces[b].size();++bf) {
			file << rans[gid].yplus.face(grid[gid].boundaryFaces[b][bf]);
			if ((bf+1)%10==0) file << "\n";
			else file << "\t";
		}
	} else if (varList[ov]=="k") {
		for (int bf=0;bf<grid[gid].boundaryFaces[b].size();++bf) {
			file << rans[gid].k.face(grid[gid].boundaryFaces[b][bf]);
			if ((bf+1)%10==0) file << "\n";
			else file << "\t";
		}
	} else if (varList[ov]=="omega") {
		for (int bf=0;bf<grid[gid].boundaryFaces[b].size();++bf) {
			file << rans[gid].omega.face(grid[gid].boundaryFaces[b][bf]);
			if ((bf+1)%10==0) file << "\n";
			else file << "\t";
		}
	} else if (varList[ov]=="mu_t") {
		for (int bf=0;bf<grid[gid].boundaryFaces[b].size();++bf) {
			file << rans[gid].mu_t.face(grid[gid].boundaryFaces[b][bf]);
			if ((bf+1)%10==0) file << "\n";
			else file << "\t";
		}
	} else if (varList[ov]=="Mach") {
		for (int bf=0;bf<grid[gid].boundaryFaces[b].size();++bf) {
			file << fabs(ns[gid].V.face(grid[gid].boundaryFaces[b][bf]))/ns[gid].material.a(ns[gid].p.face(grid[gid].boundaryFaces[b][bf]),ns[gid].T.face(grid[gid].boundaryFaces[b][bf]));
			if ((bf+1)%10==0) file << "\n";
			else file << "\t";
		}
	} else if (varList[ov]=="rank") {
		for (int bf=0;bf<grid[gid].boundaryFaces[b].size();++bf) {
			file << Rank;
			if ((bf+1)%10==0) file << "\n";
			else file << "\t";
		}
	} else if (varList[ov]=="V") {
		for (int bf=0;bf<grid[gid].boundaryFaces[b].size();++bf) {
			file << ns[gid].V.face(grid[gid].boundaryFaces[b][bf])[i];
			if ((bf+1)%10==0) file << "\n";
			else file << "\t";
		}
	} else if (varList[ov]=="gradp") {
		for (int bf=0;bf<grid[gid].boundaryFaces[b].size();++bf) {
			file << ns[gid].gradp.face(grid[gid].boundaryFaces[b][bf])[i];
			if ((bf+1)%10==0) file << "\n";
			else file << "\t";
		}
	} else if (varList[ov]=="gradu") {
		for (int bf=0;bf<grid[gid].boundaryFaces[b].size();++bf) {
			file << ns[gid].gradu.face(grid[gid].boundaryFaces[b][bf])[i];
			if ((bf+1)%10==0) file << "\n";
			else file << "\t";
		}
	} else if (varList[ov]=="gradv") {
		for (int bf=0;bf<grid[gid].boundaryFaces[b].size();++bf) {
			file << ns[gid].gradv.face(grid[gid].boundaryFaces[b][bf])[i];
			if ((bf+1)%10==0) file << "\n";
			else file << "\t";
		}
	} else if (varList[ov]=="gradw") {
		for (int bf=0;bf<grid[gid].boundaryFaces[b].size();++bf) {
			file << ns[gid].gradw.face(grid[gid].boundaryFaces[b][bf])[i];
			if ((bf+1)%10==0) file << "\n";
			else file << "\t";
		}
	} else if (varList[ov]=="gradT") {
		if (equations[gid]==NS) {
			for (int bf=0;bf<grid[gid].boundaryFaces[b].size();++bf) {
				file << ns[gid].gradT.face(grid[gid].boundaryFaces[b][bf])[i];
				if ((bf+1)%10==0) file << "\n";
				else file << "\t";
			}
		} else if (equations[gid]==HEAT) {
			for (int bf=0;bf<grid[gid].boundaryFaces[b].size();++bf) {
				file << hc[gid].gradT.face(grid[gid].boundaryFaces[b][bf])[i];
				if ((bf+1)%10==0) file << "\n";
				else file << "\t";
			}
		}
	} else if (varList[ov]=="gradrho") {
		for (int bf=0;bf<grid[gid].boundaryFaces[b].size();++bf) {
			file << ns[gid].gradrho.face(grid[gid].boundaryFaces[b][bf])[i];
			if ((bf+1)%10==0) file << "\n";
			else file << "\t";
		}
	} else if (varList[ov]=="resV") {
		for (int bf=0;bf<grid[gid].boundaryFaces[b].size();++bf) {
			file << ns[gid].update[i+1].face(grid[gid].boundaryFaces[b][bf]);
			if ((bf+1)%10==0) file << "\n";
			else file << "\t";
		}
	} else if (varList[ov]=="limiterV") {
		for (int bf=0;bf<grid[gid].boundaryFaces[b].size();++bf) {
			file << ns[gid].limiter[i+1].face(grid[gid].boundaryFaces[b][bf]);
			if ((bf+1)%10==0) file << "\n";
			else file << "\t";
		}
	} else if (varList[ov]=="gradk") {
		for (int bf=0;bf<grid[gid].boundaryFaces[b].size();++bf) {
			file << rans[gid].gradk.face(grid[gid].boundaryFaces[b][bf])[i];
			if ((bf+1)%10==0) file << "\n";
			else file << "\t";
		}
	} else if (varList[ov]=="gradomega") {
		for (int bf=0;bf<grid[gid].boundaryFaces[b].size();++bf) {
			file << rans[gid].gradomega.face(grid[gid].boundaryFaces[b][bf])[i];
			if ((bf+1)%10==0) file << "\n";
			else file << "\t";
		}
	} else if (varList[ov]=="tau") {
		for (int bf=0;bf<grid[gid].boundaryFaces[b].size();++bf) {
			file << ns[gid].tau.face(grid[gid].boundaryFaces[b][bf])[i];
			if ((bf+1)%10==0) file << "\n";
			else file << "\t";
		}
	}
	
	file.close();
		
	return;
} 
			
void write_surface_tec_cells(int b) {
	
	ofstream file;
	string fileName="./output/surface_"+int2str(timeStep)+"_"+int2str(gid+1)+".dat";
	
	file.open((fileName).c_str(),ios::app);
	
	// Transform this logic: map node indices to this output node list
	
	for (int bf=0;bf<grid[gid].boundaryFaces[b].size();++bf) {
		// Write connectivity
		int f=grid[gid].boundaryFaces[b][bf];
		
		file << grid[gid].faceNode(f,0).bc_output_id+1 << "\t" ;
		file << grid[gid].faceNode(f,1).bc_output_id+1 << "\t" ;
		file << grid[gid].faceNode(f,2).bc_output_id+1 << "\t" ;
		
		if (grid[gid].face[f].nodeCount==4) {
			file << grid[gid].faceNode(f,3).bc_output_id+1 << "\t" ;			
		} else if (grid[gid].face[f].nodeCount==3) {
			file << grid[gid].faceNode(f,2).bc_output_id+1 << "\t" ;
		}

		file << endl;
	}
	
	file.close();
	
	return;
	
}
