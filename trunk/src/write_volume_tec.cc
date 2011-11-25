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
#include <iomanip>
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

namespace volume_output {
	int timeStep,gid;
	vector<string> varList;
	vector<bool> var_is_vec3d;
}
using namespace volume_output;
		
void write_tec_header(void);
void write_tec_nodes(int i);
void write_tec_var(int ov, int i);
void write_tec_face_node_counts (void);
void write_tec_face_nodes(void);
void write_tec_left(void);
void write_tec_right(void);

void write_loads(int gid,int step,double time) {
	ofstream file;
	for (int b=0;b<loads[gid].include_bcs.size();++b) {
		double force_x,force_y,force_z,moment_x,moment_y,moment_z;
		MPI_Allreduce (&loads[gid].force[b][0],&force_x,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
		MPI_Allreduce (&loads[gid].force[b][1],&force_y,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
		MPI_Allreduce (&loads[gid].force[b][2],&force_z,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
		MPI_Allreduce (&loads[gid].moment[b][0],&moment_x,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
		MPI_Allreduce (&loads[gid].moment[b][1],&moment_y,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
		MPI_Allreduce (&loads[gid].moment[b][2],&moment_z,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
		if (Rank==0) {
			string fileName="./output/loads_grid_"+int2str(gid+1)+"_BC_"+int2str(loads[gid].include_bcs[b]+1)+".dat";
			file.open((fileName).c_str(),ios::app);
			file << scientific << setprecision(8);
			file << step << "\t" << time << "\t";
			file << force_x << "\t" << force_y << "\t" << force_z << "\t" << moment_x << "\t" << moment_y << "\t" << moment_z << endl;
			file.close();
		}
	}
	return;
}

void write_volume_output(int gridid, int step) {
	mkdir("./output",S_IRWXU);
	gid=gridid;
	timeStep=step;
	varList=input.section("grid",gid).subsection("writeoutput").get_stringList("volumevariables");
	var_is_vec3d.resize(varList.size());
	int nVar=0;
	for (int var=0; var<varList.size(); ++var) {
		if (varList[var]!="null") nVar++;
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
		if  (varList[var]=="gradrho") ns[gid].gradrho.mpi_update();
	}
	
	string format=input.section("grid",gid).subsection("writeoutput").get_string("format");
	if (format=="tecplot") {
		// Write tecplot output file
		if (Rank==0) write_tec_header();
		for (int i=0;i<3;++i) {
			for (int p=0;p<np;++p) {
				if(Rank==p) write_tec_nodes(i);
				MPI_Barrier(MPI_COMM_WORLD);
			}
		}
		for (int ov=0;ov<varList.size();++ov) {
			int nn=1;
			if (var_is_vec3d[ov]) nn=3;
			for (int i=0;i<nn;++i) {
				for (int p=0;p<np;++p) {
					if(Rank==p) write_tec_var(ov,i);
					MPI_Barrier(MPI_COMM_WORLD);
				}
				
			}
		}
		
		for (int p=0;p<np;++p) {
			if(Rank==p) write_tec_face_node_counts();
			MPI_Barrier(MPI_COMM_WORLD);
		}
		 
		for (int p=0;p<np;++p) {
			if(Rank==p) write_tec_face_nodes();
			MPI_Barrier(MPI_COMM_WORLD);
		}

		for (int p=0;p<np;++p) {
			if(Rank==p) write_tec_left();
			MPI_Barrier(MPI_COMM_WORLD);
		}

		for (int p=0;p<np;++p) {
			if(Rank==p) write_tec_right();
			MPI_Barrier(MPI_COMM_WORLD);
		}
	} 

	return;
}

void write_tec_header(void) {
	
	ofstream file;
	string fileName="./output/volume_"+int2str(timeStep)+"_"+int2str(gid+1)+".dat";
	string link_comm="ln -sf "+fileName+" ./volume_latest_"+int2str(gid+1)+".dat";
	system(link_comm.c_str());
	
	// Proc 0 creates the output file and writes variable list
	
	file.open((fileName).c_str(),ios::out);
	file << "VARIABLES = \"x\", \"y\", \"z\" ";
	int nvars=3;
	for (int var=0;var<varList.size();++var) {
		if (var_is_vec3d[var]) {
			file << ",\"" << varList[var] << "_x\" "; 
			file << ",\"" << varList[var] << "_y\" "; 
			file << ",\"" << varList[var] << "_z\" "; 
			nvars+=3;
		} else if (varList[var]!="null") {
			file << ",\"" << varList[var] << "\" "; 
			nvars++;
		}
	}
	file << endl;
	file << "ZONE T=\"Grid_" << gid+1 << "\"" <<  endl;
	file << "Nodes=" << grid[gid].globalNodeCount << ", Faces=" << grid[gid].globalFaceCount << ", Elements=" << grid[gid].globalCellCount << ", ZONETYPE=FEPolyhedron" << endl;
	file << "DATAPACKING=BLOCK" << endl;
	file << "TotalNumFaceNodes=" << grid[gid].globalNumFaceNodes << ", NumConnectedBoundaryFaces=0, TotalNumBoundaryConnections=0 " << endl;
	file << "DT=(";
	for (int var=0;var<nvars;++var) file << "DOUBLE ";
	file << ")" << endl;
	if (nvars==4) {
		file << "VARLOCATION=([4]=CELLCENTERED)" << endl;	
	} else {
		file << "VARLOCATION=([4-" << nvars << "]=CELLCENTERED)" << endl;
	}
	
	return;
}
	
void write_tec_nodes(int i) {
	
	ofstream file;
	string fileName="./output/volume_"+int2str(timeStep)+"_"+int2str(gid+1)+".dat";
	
	file.open((fileName).c_str(),ios::app);
	file << scientific << setprecision(8);

	int count=0;
	for (int n=0;n<grid[gid].nodeCount;++n) {
		// Note that some nodes are repeated in different partitions
		if (grid[gid].node[n].output_id>=grid[gid].node_output_offset) {
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
			
			
void write_tec_var(int ov, int i) {

	ofstream file;
	string fileName="./output/volume_"+int2str(timeStep)+"_"+int2str(gid+1)+".dat";
	
	file.open((fileName).c_str(),ios::app);
	file << scientific << setprecision(8);
	
	if (varList[ov]=="p") {
		for (int c=0;c<grid[gid].cellCount;++c) {
			file << ns[gid].p.cell(c);
			if ((c+1)%10==0) file << "\n";
			else file << "\t";
		}
	} else if (varList[ov]=="T") {
		if (equations[gid]==NS)	{
			for (int c=0;c<grid[gid].cellCount;++c) {
				file << ns[gid].T.cell(c);
				if ((c+1)%10==0) file << "\n";
				else file << "\t";
			}
		} else if (equations[gid]==HEAT) {
			for (int c=0;c<grid[gid].cellCount;++c) {
				file << hc[gid].T.cell(c);
				if ((c+1)%10==0) file << "\n";
				else file << "\t";
			}
		}
	} else if (varList[ov]=="rho") {
		for (int c=0;c<grid[gid].cellCount;++c) {
			file << ns[gid].rho.cell(c);
			if ((c+1)%10==0) file << "\n";
			else file << "\t";
		}
	} else if (varList[ov]=="dt") {
		for (int c=0;c<grid[gid].cellCount;++c) {
			file << dt[gid].cell(c);
			if ((c+1)%10==0) file << "\n";
			else file << "\t";
		}
	} else if (varList[ov]=="mu") {
		for (int c=0;c<grid[gid].cellCount;++c) {
			file << ns[gid].material.viscosity(ns[gid].T.cell(c));
			if ((c+1)%10==0) file << "\n";
			else file << "\t";
		}
	} else if (varList[ov]=="lambda") {
		for (int c=0;c<grid[gid].cellCount;++c) {
			file << ns[gid].material.therm_cond(ns[gid].T.cell(c));
			if ((c+1)%10==0) file << "\n";
			else file << "\t";
		}
	} else if (varList[ov]=="Cp") {
		for (int c=0;c<grid[gid].cellCount;++c) {
			file << ns[gid].material.Cp(ns[gid].T.cell(c));
			if ((c+1)%10==0) file << "\n";
			else file << "\t";
		}
	} else if (varList[ov]=="resp") {
		for (int c=0;c<grid[gid].cellCount;++c) {
			file << ns[gid].update[0].cell(c);
			if ((c+1)%10==0) file << "\n";
			else file << "\t";
		}
	} else if (varList[ov]=="resT") {
		for (int c=0;c<grid[gid].cellCount;++c) {
			file << ns[gid].update[4].cell(c);
			if ((c+1)%10==0) file << "\n";
			else file << "\t";
		}
	} else if (varList[ov]=="limiterp") {
		for (int c=0;c<grid[gid].cellCount;++c) {
			file << ns[gid].limiter[0].cell(c);
			if ((c+1)%10==0) file << "\n";
			else file << "\t";
		}
	} else if (varList[ov]=="limiterT") {
		for (int c=0;c<grid[gid].cellCount;++c) {
			file << ns[gid].limiter[4].cell(c);
			if ((c+1)%10==0) file << "\n";
			else file << "\t";
		}
	} else if (varList[ov]=="k") {
		for (int c=0;c<grid[gid].cellCount;++c) {
			file << rans[gid].k.cell(c);
			if ((c+1)%10==0) file << "\n";
			else file << "\t";
		}
	} else if (varList[ov]=="omega") {
		for (int c=0;c<grid[gid].cellCount;++c) {
			file << rans[gid].omega.cell(c);
			if ((c+1)%10==0) file << "\n";
			else file << "\t";
		}
	} else if (varList[ov]=="mu_t") {
		for (int c=0;c<grid[gid].cellCount;++c) {
			file << rans[gid].mu_t.cell(c);
			if ((c+1)%10==0) file << "\n";
			else file << "\t";
		}
	} else if (varList[ov]=="Mach") {
		for (int c=0;c<grid[gid].cellCount;++c) {
			file << fabs(ns[gid].V.cell(c))/ns[gid].material.a(ns[gid].p.cell(c),ns[gid].T.cell(c));
			if ((c+1)%10==0) file << "\n";
			else file << "\t";
		}
	} else if (varList[ov]=="rank") {
		for (int c=0;c<grid[gid].cellCount;++c) {
			file << Rank;
			if ((c+1)%10==0) file << "\n";
			else file << "\t";
		}
	} else if (varList[ov]=="volume") {
		for (int c=0;c<grid[gid].cellCount;++c) {
			file << grid[gid].cell[c].volume;
			if ((c+1)%10==0) file << "\n";
			else file << "\t";
		}
	} else if (varList[ov]=="percent_grad_error") {
		if (gradient_test==LINEAR) {
			for (int c=0;c<grid[gid].cellCount;++c) {
				file << (ns[gid].gradp.cell(c)[0]-1.)*100.;
				if ((c+1)%10==0) file << "\n";
				else file << "\t";
			}
		} else if (gradient_test==QUADRATIC) {
 			for (int c=0;c<grid[gid].cellCount;++c) {
				file << (ns[gid].gradp.cell(c)[0]-(2.*grid[gid].cell[c].centroid[0]+3.*(max_x-min_x)))*100.;
				if ((c+1)%10==0) file << "\n";
				else file << "\t";
			}
		}
	} else if (varList[ov]=="V") {
		for (int c=0;c<grid[gid].cellCount;++c) {
			file << ns[gid].V.cell(c)[i];
			if ((c+1)%10==0) file << "\n";
			else file << "\t";
		}
	} else if (varList[ov]=="gradp") {
		for (int c=0;c<grid[gid].cellCount;++c) {
			file << ns[gid].gradp.cell(c)[i];
			if ((c+1)%10==0) file << "\n";
			else file << "\t";
		}
	} else if (varList[ov]=="gradu") {
		for (int c=0;c<grid[gid].cellCount;++c) {
			file << ns[gid].gradu.cell(c)[i];
			if ((c+1)%10==0) file << "\n";
			else file << "\t";
		}
	} else if (varList[ov]=="gradv") {
		for (int c=0;c<grid[gid].cellCount;++c) {
			file << ns[gid].gradv.cell(c)[i];
			if ((c+1)%10==0) file << "\n";
			else file << "\t";
		}
	} else if (varList[ov]=="gradw") {
		for (int c=0;c<grid[gid].cellCount;++c) {
			file << ns[gid].gradw.cell(c)[i];
			if ((c+1)%10==0) file << "\n";
			else file << "\t";
		}
	} else if (varList[ov]=="gradT") {
		if (equations[gid]==NS) {
			for (int c=0;c<grid[gid].cellCount;++c) {
				file << ns[gid].gradT.cell(c)[i];
				if ((c+1)%10==0) file << "\n";
				else file << "\t";
			}
		} else if (equations[gid]==HEAT) {
			for (int c=0;c<grid[gid].cellCount;++c) {
				file << hc[gid].gradT.cell(c)[i];
				if ((c+1)%10==0) file << "\n";
				else file << "\t";
			}
		}
	} else if (varList[ov]=="gradrho") {
		for (int c=0;c<grid[gid].cellCount;++c) {
			file << ns[gid].gradrho.cell(c)[i];
			if ((c+1)%10==0) file << "\n";
			else file << "\t";
		}
	} else if (varList[ov]=="resV") {
		for (int c=0;c<grid[gid].cellCount;++c) {
			file << ns[gid].update[i+1].cell(c);
			if ((c+1)%10==0) file << "\n";
			else file << "\t";
		}
	} else if (varList[ov]=="limiterV") {
		for (int c=0;c<grid[gid].cellCount;++c) {
			file << ns[gid].limiter[i+1].cell(c);
			if ((c+1)%10==0) file << "\n";
			else file << "\t";
		}
	} else if (varList[ov]=="gradk") {
		for (int c=0;c<grid[gid].cellCount;++c) {
			file << rans[gid].gradk.cell(c)[i];
			if ((c+1)%10==0) file << "\n";
			else file << "\t";
		}
	} else if (varList[ov]=="gradomega") {
		for (int c=0;c<grid[gid].cellCount;++c) {
			file << rans[gid].gradomega.cell(c)[i];
			if ((c+1)%10==0) file << "\n";
			else file << "\t";
		}
	} else if (varList[ov]=="grad") {
		for (int c=0;c<grid[gid].cellCount;++c) {
			file << ns[gid].gradp.cell(c)[i];
			if ((c+1)%10==0) file << "\n";
			else file << "\t";
		}
	}
	
	file.close();
		
	return;
} 

void write_tec_face_node_counts () {

	ofstream file;
	string fileName="./output/volume_"+int2str(timeStep)+"_"+int2str(gid+1)+".dat";
	file.open((fileName).c_str(),ios::app);

	if (Rank==0) file << "\n# node count per face" << endl;
	int count=0;
	int g;
	bool write;
	for (int f=0;f<grid[gid].faceCount;++f) {
		write=true;
		if (grid[gid].face[f].bc==GHOST_FACE) {
			g=-1*grid[gid].face[f].neighbor-1;
			if (grid[gid].ghost[g].partition<Rank) write=false;
		}
		if (write) {
			file << grid[gid].face[f].nodeCount << endl;
		}
	}
	

	file.close();
	
	return;
}
			
void write_tec_face_nodes() {
	
	ofstream file;
	string fileName="./output/volume_"+int2str(timeStep)+"_"+int2str(gid+1)+".dat";
	file.open((fileName).c_str(),ios::app);

	if (Rank==0) file << "# face nodes" << endl;
	int g;
	bool write;	
	for (int f=0;f<grid[gid].faceCount;++f) { 
		write=true;
		if (grid[gid].face[f].bc==GHOST_FACE) {
			g=-1*grid[gid].face[f].neighbor-1;
			if (grid[gid].ghost[g].partition<Rank) write=false;
		}

		if (write) {
			for (int fn=0;fn<grid[gid].face[f].nodeCount;++fn) { 
				file << grid[gid].faceNode(f,fn).output_id+1 << " ";
			}
			file << endl;
		}
	}
				
	file.close();
	
	return;
	
}

void write_tec_left() {
	
	ofstream file;
	string fileName="./output/volume_"+int2str(timeStep)+"_"+int2str(gid+1)+".dat";
	file.open((fileName).c_str(),ios::app);

	if (Rank==0) file << "# left elements" << endl;
	int g;
	bool write;
	for (int f=0;f<grid[gid].faceCount;++f) { 
		write=true;
		if (grid[gid].face[f].bc==GHOST_FACE) {
			g=-1*grid[gid].face[f].neighbor-1;
			if (grid[gid].ghost[g].partition<Rank) write=false;
		}

		if (write) file << grid[gid].face[f].parent+grid[gid].partitionOffset[Rank]+1 << endl;
	}

	file.close();
	
	return;
}

void write_tec_right() {
	
	ofstream file;
	string fileName="./output/volume_"+int2str(timeStep)+"_"+int2str(gid+1)+".dat";
	file.open((fileName).c_str(),ios::app);

	if (Rank==0) file << "# right elements" << endl;
	int g;
	bool write;	
	for (int f=0;f<grid[gid].faceCount;++f) { 
		write=true;
		if (grid[gid].face[f].bc==GHOST_FACE) {
			g=-1*grid[gid].face[f].neighbor-1;
			if (grid[gid].ghost[g].partition<Rank) write=false;
		}

		if (write) {
			if (grid[gid].face[f].bc==GHOST_FACE) file << grid[gid].ghost[g].id_in_owner+grid[gid].partitionOffset[grid[gid].ghost[g].partition]+1 << endl;
			else if (grid[gid].face[f].bc>=0) file << 0 << endl;
			else file << grid[gid].face[f].neighbor+grid[gid].partitionOffset[Rank]+1 << endl;
		}
	}

	file.close();
	
	return;
}
