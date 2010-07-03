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
#include <cmath>
#include <iomanip>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <limits>
#include <mpi.h>
#include <vector>
#include "inputs.h"
using namespace std;
#include "grid.h"
#include "ns.h"
#include "hc.h"
#include "bc.h"
#include "commons.h"
#include "bc_interface.h"

// Function prototypes
void read_inputs(void);
void set_bcs(int gid);
void write_output(int gid, int step);
void write_restart(int gid,int timeStep,int restart_step,double time);
void read_restart(int gid,int restart_step,double time);
void set_lengthScales(int gid);
void set_time_step_options(void);
void update_time_step(int timeStep,int gid);
void bc_interface_sync(void);

// Global declerations
InputFile input;
vector<InputFile> material_input;
vector<Grid> grid;
vector<vector<BCregion> > bc;
vector<NavierStokes> ns;
vector<HeatConduction> hc;
// Time step for each grid
vector<Variable<double> > dt;
vector<int> equations;
vector<vector<BC_Interface> > interface; // for each grid
int Rank,np;

// Equation options
#define NONE -1
#define NS 1
#define HEAT 2

static char help[] = "Free CFD\n - A free general purpose computational fluid dynamics code";

int main(int argc, char *argv[]) {

    // Initialize mpi
	MPI_Init(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &Rank);
	MPI_Comm_size(MPI_COMM_WORLD, &np);
	
	// Initialize petsc
	PetscInitialize(&argc,&argv,(char *)0,help);
	
	string inputFileName;
	inputFileName.assign(argv[1]);
	inputFileName+=".in";

	int restart_step=0;
	if (argc>2) restart_step=atoi(argv[2]);
	
	// Read the input file
	input.setFile(inputFileName);
	read_inputs();
	
	equations.resize(input.section("grid",0).count);
	for (int gid=0;gid<input.section("grid",0).count;++gid) {
		if (input.section("grid",gid).get_string("equations")=="navierstokes") {
			equations[gid]=NS;
		} else if (input.section("grid",gid).get_string("equations")=="heatconduction") {
			equations[gid]=HEAT;
		}
	}

	// Allocate space for each grid and equations and BC's on each grid
	grid.resize(input.section("grid",0).count);
	ns.resize(grid.size());
	hc.resize(grid.size());
	bc.resize(grid.size());
	interface.resize(grid.size());
	
	// Read the grid and initialize
	for (int gid=0;gid<grid.size();++gid) {
		grid[gid].dimension=input.section("grid",gid).get_int("dimension");
		grid[gid].gid=gid;
		// Read the grid raw data from file
		grid[gid].read(input.section("grid",gid).get_string("file"));
		// Do the transformations
		int tcount=input.section("grid",gid).subsection("transform",0).count;
		
		for (int t=0;t<tcount;++t) {
			if (input.section("grid",gid).subsection("transform",t).get_string("function")=="translate") {
				Vec3D begin=input.section("grid",gid).subsection("transform",t).get_Vec3D("begin");
				Vec3D end=input.section("grid",gid).subsection("transform",t).get_Vec3D("end");
				grid[gid].translate(begin,end);
			} else if (input.section("grid",gid).subsection("transform",t).get_string("function")=="scale") {
				Vec3D anchor=input.section("grid",gid).subsection("transform",t).get_Vec3D("anchor");
				double factor=input.section("grid",gid).subsection("transform",t).get_double("factor");
				grid[gid].scale(anchor,factor);
			} else if (input.section("grid",gid).subsection("transform",t).get_string("function")=="rotate") {
				Vec3D anchor=input.section("grid",gid).subsection("transform",t).get_Vec3D("anchor");
				Vec3D axis=input.section("grid",gid).subsection("transform",t).get_Vec3D("axis");
				double angle=input.section("grid",gid).subsection("transform",t).get_double("angle");
				grid[gid].rotate(anchor,axis,angle);
			} 
		}
		// Establish connectivity, area, volume etc... all the needed information
		grid[gid].setup();
		set_bcs(gid);
		set_lengthScales(gid);
		if (Rank==0) cout << "[I grid=" << gid+1 << " ] Calculating node averaging metrics, this might take a while..." << endl;
		grid[gid].nodeAverages(); // Linear triangular (tetrahedral) + idw blended mode
        grid[gid].faceAverages();
	}

	for (int gid=0;gid<grid.size();++gid) {
		for (int i=0;i<interface[gid].size();++i) {
			interface[gid][i].setup();
		}
	}

	double time=0.;
	for (int gid=0;gid<grid.size();++gid) {
		if (equations[gid]==NS) {
			if (Rank==0) cout << "[I grid=" << gid+1 << " ] Initializing Navier Stokes solver" << endl; 
			ns[gid].gid=gid;
			ns[gid].initialize();
		}
		if (equations[gid]==HEAT) {
			if (Rank==0) cout << "[I grid=" << gid+1 << " ] Initializing Heat Conduction solver" << endl;
			hc[gid].gid=gid;
			hc[gid].initialize();
		}
		if (restart_step>0) read_restart(gid,restart_step,time);
		set_time_step_options();
	}
	
	if (Rank==0) cout << "[I] Beginning time loop\n" << endl; 
	int timeStepMax=input.section("timemarching").get_int("numberofsteps");
	// Get the output frequency for each grid
	vector<int> plotFreq,restartFreq;
	for (int gid=0;gid<grid.size();++gid) {
		plotFreq.push_back(input.section("grid",gid).subsection("writeoutput").get_int("plotfrequency"));
		restartFreq.push_back(input.section("grid",gid).subsection("writeoutput").get_int("restartfrequency"));	
	}
	/*****************************************************************************************/
	// Begin time loop
	/*****************************************************************************************/
	
    bool lastTimeStep=false;
	cout << setprecision(3);
	for (int timeStep=restart_step+1;timeStep<=timeStepMax+restart_step;++timeStep) {
		if (timeStep==(timeStepMax+restart_step)) lastTimeStep=true;
		for (int gid=0;gid<grid.size();++gid) {
			update_time_step(timeStep,gid);
			if (equations[gid]==NS) ns[gid].solve(timeStep);
			if (equations[gid]==HEAT) hc[gid].solve(timeStep);
			bc_interface_sync();
			if (timeStep%plotFreq[gid]==0 || lastTimeStep) {
				if (Rank==0) cout << "\nwriting output for grid=" << gid+1;
				  write_output(gid,timeStep);
			} // end if
			if (timeStep%restartFreq[gid]==0 || lastTimeStep) {
				if (Rank==0) cout << "\nwriting restart for grid=" << gid+1;
				write_restart(gid,timeStep,restart_step,time);
			} // end if
		} // end grid loop
		if (Rank==0) cout << endl;
	}
	/*****************************************************************************************/
	// End time loop
	/*****************************************************************************************/	
		
	PetscFinalize();
	
	return 0;
}                	
                 	
void set_lengthScales(int gid) {
	// Loop through the cells and calculate length scales for each cell
	for (int c=0;c<grid[gid].cellCount;++c) {
		grid[gid].cell[c].lengthScale=1.e20;
		double height;
		int f;
		for (int cf=0;cf<grid[gid].cell[c].faceCount;++cf) {
			f=grid[gid].cell[c].faces[cf];
			height=fabs(grid[gid].face[f].normal.dot(grid[gid].face[f].centroid-grid[gid].cell[c].centroid));
			bool skipScale=false;
			if (grid[gid].face[f].bc>=0) {
				if (bc[gid][grid[gid].face[f].bc].type==SYMMETRY) {
					skipScale=true;
				}
			}
			if (!skipScale) grid[gid].cell[c].lengthScale=min(grid[gid].cell[c].lengthScale,height);
		}
	}
	
	// Find out the global, grid cell length scale
	if (grid[gid].dimension==3) {
		grid[gid].lengthScale=pow(grid[gid].globalTotalVolume,1./3.);
	} else {
		// If the problem is 2D, finding the grid length scale is a bit more challenging
		// Find the symmetry BC region with the largest area, and take the square root
		double maxTotalArea=0.;
		double totalArea;
		for (int b=0;b<bc[gid].size();++b) {
			totalArea=0.;
			if (bc[gid][b].kind==SYMMETRY) {
				MPI_Allreduce (&bc[gid][b].area,&totalArea,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
				maxTotalArea=max(maxTotalArea,totalArea);
			}
		}
		grid[gid].lengthScale=sqrt(maxTotalArea);
	}
	
	return;
} 
