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

int timeStep,gid;
vector<string> varList;
vector<bool> var_is_vec3d;

void write_tec_vars(void);
void write_tec_cells(void);
void write_tec_bcs(int nbc, int nVar);

void write_output(int gridid, int step) {
	mkdir("./output",S_IRWXU);
	gid=gridid;
	timeStep=step;
	varList=input.section("grid",gid).subsection("writeoutput").get_stringList("variables");
	var_is_vec3d.resize(varList.size());
	int nVar=0;
	for (int var=0; var<varList.size(); ++var) {
		nVar++;
		var_is_vec3d[var]=false;
		if (varList[var]=="V" || varList[var].substr(0,4)=="grad" || varList[var]=="resV" || varList[var]=="limiterV") {
			var_is_vec3d[var]=true;
			nVar+=2;
		}
		// The following ghost values only need to be updated when writing node centered output
		if  (varList[var]=="dt") dt[gid].mpi_update(); 
		if  (varList[var]=="resp") ns[gid].update[0].mpi_update(); 
		if  (varList[var]=="resT") ns[gid].update[4].mpi_update();
		if  (varList[var]=="resV") for (int i=1;i<4;++i) ns[gid].update[i].mpi_update();
		if  (varList[var]=="limiterp") ns[gid].limiter[0].mpi_update(); 
		if  (varList[var]=="limiterT") ns[gid].limiter[4].mpi_update();
		if  (varList[var]=="limiterV") for (int i=1;i<4;++i) ns[gid].limiter[i].mpi_update();

	}
	// Write tecplot output file
	for (int p=0;p<np;++p) {
		if(Rank==p) write_tec_vars();
		MPI_Barrier(MPI_COMM_WORLD);
	}
	for (int p=0;p<np;++p) {
		if(Rank==p) write_tec_cells();
		MPI_Barrier(MPI_COMM_WORLD);
	}
	// Write boundary condition regions
	for (int nbc=0;nbc<input.section("grid",gid).subsection("BC",0).count;++nbc) {
		for (int p=0;p<np;++p) {
			if(Rank==p) write_tec_bcs(nbc,nVar);
			MPI_Barrier(MPI_COMM_WORLD);
		}
	}

	return;
}

void write_tec_vars(void) {
	
	ofstream file;
	string fileName="./output/out_"+int2str(timeStep)+"_"+int2str(gid+1)+".dat";
	
	// Proc 0 creates the output file and writes variable list
	
	if (Rank==0) {
		file.open((fileName).c_str(),ios::out);
		file << "VARIABLES = \"x\", \"y\", \"z\" ";
		for (int var=0;var<varList.size();++var) {
			if (var_is_vec3d[var]) {
				file << ",\"" << varList[var] << "_x\" "; 
				file << ",\"" << varList[var] << "_y\" "; 
				file << ",\"" << varList[var] << "_z\" "; 
			} else {
				file << ",\"" << varList[var] << "\" "; 
			}
		}
		file << endl;
		file << "ZONE, T=\"Grid " << gid+1 << "\", ZONETYPE=FEBRICK, DATAPACKING=POINT" << endl;
		file << "NODES=" << grid[gid].globalNodeCount << ", ELEMENTS=" << grid[gid].globalCellCount << endl;
	} else {
		file.open((fileName).c_str(),ios::app);
	}
	
	file << scientific << setprecision(8);
	
	for (int n=0;n<grid[gid].nodeCount;++n) {
		// Note that some nodes are repeated in different partitions
		if (grid[gid].maps.nodeGlobal2Output[grid[gid].node[n].globalId]>=grid[gid].nodeCountOffset) {
			
			// Write node coordinates
			file << grid[gid].node[n][0] << "\t";
			file << grid[gid].node[n][1] << "\t";
			file << grid[gid].node[n][2] << "\t";
			// Write Navier-Stokes variables
			
			// Node temperature and pressure are repeatedly used. Store them first
			double p_node,T_node;
			
			if (equations[gid]==NS)	{p_node=ns[gid].p.node(n); T_node=ns[gid].T.node(n);}
			else if (equations[gid]==HEAT) T_node=hc[gid].T.node(n);
			
			for (int ov=0;ov<varList.size();++ov) {
				// Scalars
				if (varList[ov]=="p") file << p_node << "\t";
				else if (varList[ov]=="T") file << T_node << "\t";
				else if (varList[ov]=="rho") file << ns[gid].material.rho(p_node,T_node) << "\t";
				else if (varList[ov]=="dt") file << dt[gid].node(n) << "\t";
				else if (varList[ov]=="mu") file << ns[gid].material.viscosity(T_node) << "\t";
				else if (varList[ov]=="lambda") file << ns[gid].material.therm_cond(T_node) << "\t";	
				else if (varList[ov]=="Cp") file << ns[gid].material.Cp(T_node) << "\t";
				else if (varList[ov]=="resp") file << ns[gid].update[0].node(n) << "\t";
				else if (varList[ov]=="resT") file << ns[gid].update[4].node(n) << "\t";
				else if (varList[ov]=="limiterp") file << ns[gid].limiter[0].node(n) << "\t"; 
				else if (varList[ov]=="limiterT") file << ns[gid].limiter[4].node(n) << "\t";
				else if (varList[ov]=="k") file << rans[gid].k.node(n) << "\t";
				else if (varList[ov]=="omega") file << rans[gid].omega.node(n) << "\t";
				else if (varList[ov]=="mu_t") file << rans[gid].mu_t.node(n) << "\t";
				else if (varList[ov]=="Mach") {
					file << fabs(ns[gid].V.node(n))/ns[gid].material.a(p_node,T_node) << "\t";
				}
				// Vectors
				else if (varList[ov]=="V") {
					file << ns[gid].V.node(n)[0] << "\t";
					file << ns[gid].V.node(n)[1] << "\t";
					file << ns[gid].V.node(n)[2] << "\t";
				} 
				else if (varList[ov]=="gradp") {
					file << ns[gid].gradp.node(n)[0] << "\t";
					file << ns[gid].gradp.node(n)[1] << "\t";
					file << ns[gid].gradp.node(n)[2] << "\t";
				}
				else if (varList[ov]=="gradu") {
					file << ns[gid].gradu.node(n)[0] << "\t";
					file << ns[gid].gradu.node(n)[1] << "\t";
					file << ns[gid].gradu.node(n)[2] << "\t";
				}
				else if (varList[ov]=="gradv") {
					file << ns[gid].gradv.node(n)[0] << "\t";
					file << ns[gid].gradv.node(n)[1] << "\t";
					file << ns[gid].gradv.node(n)[2] << "\t";
				}
				else if (varList[ov]=="gradw") {
					file << ns[gid].gradw.node(n)[0] << "\t";
					file << ns[gid].gradw.node(n)[1] << "\t";
					file << ns[gid].gradw.node(n)[2] << "\t";
				}
				else if (varList[ov]=="gradT") {
					if (equations[gid]==NS) {
						file << ns[gid].gradT.node(n)[0] << "\t";
						file << ns[gid].gradT.node(n)[1] << "\t";
						file << ns[gid].gradT.node(n)[2] << "\t";
					} else if (equations[gid]==HEAT) {
						file << hc[gid].gradT.node(n)[0] << "\t";
						file << hc[gid].gradT.node(n)[1] << "\t";
						file << hc[gid].gradT.node(n)[2] << "\t";
					}
				}				
				else if (varList[ov]=="gradrho") {
					file << ns[gid].gradrho.node(n)[0] << "\t";
					file << ns[gid].gradrho.node(n)[1] << "\t";
					file << ns[gid].gradrho.node(n)[2] << "\t";
				}
				else if (varList[ov]=="resV") {
					file << ns[gid].update[1].node(n) << "\t";
					file << ns[gid].update[2].node(n) << "\t";
					file << ns[gid].update[3].node(n) << "\t";
				}
				else if (varList[ov]=="limiterV") {
					file << ns[gid].limiter[1].node(n) << "\t";
					file << ns[gid].limiter[2].node(n) << "\t";
					file << ns[gid].limiter[3].node(n) << "\t";
				}
				else if (varList[ov]=="gradk") {
					file << rans[gid].gradk.node(n)[0] << "\t";
					file << rans[gid].gradk.node(n)[1] << "\t";
					file << rans[gid].gradk.node(n)[2] << "\t";
				}
				else if (varList[ov]=="gradomega") {
					file << rans[gid].gradomega.node(n)[0] << "\t";
					file << rans[gid].gradomega.node(n)[1] << "\t";
					file << rans[gid].gradomega.node(n)[2] << "\t";
				}				
			}
		}
		file << endl;
	}
		 
	file.close();
}

void write_tec_cells() {
	
	ofstream file;
	string fileName="./output/out_"+int2str(timeStep)+"_"+int2str(gid+1)+".dat";
	
	file.open((fileName).c_str(),ios::app);
	
	// Write connectivity
	for (int c=0;c<grid[gid].cellCount;++c) {
		if (grid[gid].cell[c].nodeCount==4) {
			file << grid[gid].maps.nodeGlobal2Output[grid[gid].cellNode(c,0).globalId]+1 << "\t" ;
			file << grid[gid].maps.nodeGlobal2Output[grid[gid].cellNode(c,2).globalId]+1 << "\t" ;
			file << grid[gid].maps.nodeGlobal2Output[grid[gid].cellNode(c,1).globalId]+1 << "\t" ;
			file << grid[gid].maps.nodeGlobal2Output[grid[gid].cellNode(c,1).globalId]+1 << "\t" ;
			file << grid[gid].maps.nodeGlobal2Output[grid[gid].cellNode(c,3).globalId]+1 << "\t" ;
			file << grid[gid].maps.nodeGlobal2Output[grid[gid].cellNode(c,3).globalId]+1 << "\t" ;
			file << grid[gid].maps.nodeGlobal2Output[grid[gid].cellNode(c,3).globalId]+1 << "\t" ;
			file << grid[gid].maps.nodeGlobal2Output[grid[gid].cellNode(c,3).globalId]+1 << "\t" ;
		}
		else if (grid[gid].cell[c].nodeCount==5) {
			file << grid[gid].maps.nodeGlobal2Output[grid[gid].cellNode(c,0).globalId]+1 << "\t" ;
			file << grid[gid].maps.nodeGlobal2Output[grid[gid].cellNode(c,1).globalId]+1 << "\t" ;
			file << grid[gid].maps.nodeGlobal2Output[grid[gid].cellNode(c,2).globalId]+1 << "\t" ;
			file << grid[gid].maps.nodeGlobal2Output[grid[gid].cellNode(c,3).globalId]+1 << "\t" ;
			file << grid[gid].maps.nodeGlobal2Output[grid[gid].cellNode(c,4).globalId]+1 << "\t" ;
			file << grid[gid].maps.nodeGlobal2Output[grid[gid].cellNode(c,4).globalId]+1 << "\t" ;
			file << grid[gid].maps.nodeGlobal2Output[grid[gid].cellNode(c,4).globalId]+1 << "\t" ;
			file << grid[gid].maps.nodeGlobal2Output[grid[gid].cellNode(c,4).globalId]+1 << "\t" ;
		}
		else if (grid[gid].cell[c].nodeCount==6) {
			file << grid[gid].maps.nodeGlobal2Output[grid[gid].cellNode(c,0).globalId]+1 << "\t" ;
			file << grid[gid].maps.nodeGlobal2Output[grid[gid].cellNode(c,1).globalId]+1 << "\t" ;
			file << grid[gid].maps.nodeGlobal2Output[grid[gid].cellNode(c,2).globalId]+1 << "\t" ;
			file << grid[gid].maps.nodeGlobal2Output[grid[gid].cellNode(c,2).globalId]+1 << "\t" ;
			file << grid[gid].maps.nodeGlobal2Output[grid[gid].cellNode(c,3).globalId]+1 << "\t" ;
			file << grid[gid].maps.nodeGlobal2Output[grid[gid].cellNode(c,4).globalId]+1 << "\t" ;
			file << grid[gid].maps.nodeGlobal2Output[grid[gid].cellNode(c,5).globalId]+1 << "\t" ;
			file << grid[gid].maps.nodeGlobal2Output[grid[gid].cellNode(c,5).globalId]+1 << "\t" ;
		} else if (grid[gid].cell[c].nodeCount==8) {
			for (int i=0;i<8;++i) {
				file << grid[gid].maps.nodeGlobal2Output[grid[gid].cellNode(c,i).globalId]+1 << "\t" ;
			}
		}
		file << endl;
	}
	
	file.close();
	
	return;
	
}

void write_tec_bcs(int bcNo,int nVar) {
	
	ofstream file;
	string fileName="./output/out_"+int2str(timeStep)+"_"+int2str(gid+1)+".dat";

	file.open((fileName).c_str(),ios::app);
	if (Rank==0) {
		file << "ZONE T=\"Grid " << gid+1 << " BC_" << bcNo+1 << "\" " << ", ZONETYPE=FEQUADRILATERAL, VARSHARELIST=([1-" << nVar+3 << "]=1)" << endl;
		file << "NODES=" << grid[gid].globalNodeCount << ", ELEMENTS=" << grid[gid].globalBoundaryFaceCount[bcNo] << endl;
	}
	
	// Write connectivity
	for (int f=0;f<grid[gid].faceCount;++f) {
	   if (grid[gid].face[f].bc==bcNo) {
		   if (grid[gid].face[f].nodeCount==3) {
			   file << grid[gid].maps.nodeGlobal2Output[grid[gid].faceNode(f,0).globalId]+1 << "\t" ;
			   file << grid[gid].maps.nodeGlobal2Output[grid[gid].faceNode(f,1).globalId]+1 << "\t" ;
			   file << grid[gid].maps.nodeGlobal2Output[grid[gid].faceNode(f,2).globalId]+1 << "\t" ;
			   file << grid[gid].maps.nodeGlobal2Output[grid[gid].faceNode(f,2).globalId]+1 << "\t" ;
		   } else if (grid[gid].face[f].nodeCount==4) {
			   for (int i=0;i<4;++i) {
				   file << grid[gid].maps.nodeGlobal2Output[grid[gid].faceNode(f,i).globalId]+1 << "\t" ;
			   }
		   }
		   file << endl;
	   }
	}

	file.close();

	return;

}
	
