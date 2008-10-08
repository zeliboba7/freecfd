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
#include <map>
#include <cmath>
#include <iomanip>
using namespace std;

#include "grid.h"
#include "inputs.h"
#include "bc.h"
#include "mpi_functions.h"
#include "petsc_functions.h"
#include "probe.h"

// Function prototypes
void read_inputs(InputFile &input);
void initialize(Grid &grid, InputFile &input);
void write_output(int timeStep, double time, InputFile input);
void write_restart(int timeStep, double time);
void read_restart(int restart, int global2local[], double &time);
void set_bcs(Grid& grid, InputFile &input, BC &bc);
bool within_box(Vec3D centroid, Vec3D box_1, Vec3D box_2);
double minmod(double a, double b);
double maxmod(double a, double b);
double superbee(double a, double b);
void update(double dt);
void updatePrimitive(double dt);
string int2str(int number) ;
void assemble_linear_system(void);
void initialize_linear_system(void);
void get_dt(string type);
void get_CFL(void);
void set_probes(void);
void set_loads(void);

// Initiate grid class
Grid grid;
// Allocate container for boundary conditions
BC bc;
InputFile input;

int np, Rank;
double Gamma,dt,CFL,CFLtarget;
double Pref;
string fluxFunction;
vector<Probe> probes;
int timeStep,restart;

bool grad_test=false; // DEBUG

int main(int argc, char *argv[]) {

	// Initialize mpi
	mpi_init(argc,argv);

	string inputFileName, gridFileName;
	inputFileName.assign(argv[1]);
	gridFileName.assign(argv[1]);
	inputFileName+=".in";
	gridFileName+=".cgns";

	restart=0;
	if (argc>2) restart=atoi(argv[2]);

	input.setFile(inputFileName);
	read_inputs(input);
	
	double time = 0.;
	Gamma = input.section["fluidProperties"].doubles["gamma"];
	dt = input.section["timeMarching"].doubles["step"];
	CFLtarget= input.section["timeMarching"].doubles["CFL"];
	if (restart!=0 && input.section["timeMarching"].strings["type"]=="CFLramp") {
		input.section["timeMarching"].strings["type"]="CFL";
	}
	int timeStepMax = input.section["timeMarching"].ints["numberOfSteps"];
	double mu;
	if (input.section["fluidProperties"].subsections["viscosity"].strings["type"]=="fixed") {
		mu=input.section["fluidProperties"].subsections["viscosity"].doubles["value"];
	}
	Pref=input.section["fluidProperties"].doubles["Pref"];
	string limiter=input.section["numericalOptions"].strings["limiter"];
	string order=input.section["numericalOptions"].strings["order"];
	double sharpeningFactor=input.section["numericalOptions"].doubles["sharpeningFactor"];
	int outFreq = input.section["output"].ints["outFreq"];
	int restartFreq = input.section["output"].ints["restartFreq"];
	fluxFunction=input.section["numericalOptions"].strings["flux"];
	int probeFreq=input.section["probes"].ints["frequency"];
	int probeCount=input.section["probes"].numberedSubsections["probe"].count;
	int loadFreq=input.section["loads"].ints["frequency"];
	int loadCount=input.section["loads"].numberedSubsections["load"].count;
	
	grid.read(gridFileName);
	
	// Hand shake with other processors
	// (agree on which cells' data to be shared)
	mpi_handshake();
	// Map global cell id's to local id's
	// Fills in array global2local[globalID]=localId
	// If the queried globalId is not on this partition, returns -1 as localId
	// If the queried globalId belongs to a ghost cell, return local ghost id
	mpi_map_global2local();
	
	initialize(grid,input);
	if (Rank==0) cout << "[I] Applied initial conditions" << endl;

	set_bcs(grid,input,bc);
	if (Rank==0) cout << "[I] Set boundary conditions" << endl;

	set_probes();

	grid.lengthScales();
	grid.nodeAverages();
	grid.faceAverages();
	grid.gradMaps();
	
	if (restart!=0) {
		for (int p=0;p<np;++p) {
			if (Rank==p) read_restart(restart,global2local,time);
			MPI_Barrier(MPI_COMM_WORLD);
		}
	}

	// Initialize petsc
	petsc_init(argc,argv,
				input.section["linearSolver"].doubles["relTolerance"],
				input.section["linearSolver"].doubles["absTolerance"],
				input.section["linearSolver"].ints["maxIterations"]);

	MPI_Barrier(MPI_COMM_WORLD);
	double timeRef, timeEnd;
	timeRef=MPI_Wtime();

	unsigned int parent,neighbor;
	
	double CFLramp;

	
	if (input.section["timeMarching"].strings["type"]=="CFLramp") CFL=1.;
	if (input.section["timeMarching"].strings["type"]=="CFL") CFL=CFLtarget;
	if (input.section["timeMarching"].strings["type"]=="CFLlocal") CFL=CFLtarget;

	if (Rank==0) cout << "[I] Beginning time loop" << endl;
	
	/*****************************************************************************************/
	// Begin time loop
	/*****************************************************************************************/
	
	for (timeStep=restart+1;timeStep<=input.section["timeMarching"].ints["numberOfSteps"]+restart;++timeStep) {

		int nIter; // Number of linear solver iterations
		double rNorm; // Residual norm of the linear solver
if (grad_test) { // DEBUG
	grid.gradients(); // DEBUG
	write_restart(timeStep,time); // DEBUG
	break;
} else {
		// Update the primitive variables of the ghost cells
		mpi_update_ghost_primitives();

		// Get time step (will handle fixed, CFL and CFLramp)
		get_dt(input.section["timeMarching"].strings["type"]);	

		// Gradient are calculated for NS and/or second order schemes
		if (input.section["equations"].strings["set"]=="NS" |
			input.section["numericalOptions"].strings["order"]=="second" ) {
			// Calculate all the cell gradients for each variable
			grid.gradients();
			// Update gradients of the ghost cells
			mpi_update_ghost_gradients();
		}
		
		// Limit gradients (limited gradients are stored separately)
		if (input.section["numericalOptions"].strings["order"]=="second" ) {
			grid.limit_gradients(limiter,sharpeningFactor);
			// Update limited gradients of the ghost cells
			mpi_update_ghost_limited_gradients();
		}

		initialize_linear_system();
		assemble_linear_system();
		petsc_solve(nIter,rNorm);
		updatePrimitive(dt);

		// Advance physical time
		time += dt;
		if (input.section["timeMarching"].strings["type"]=="fixed") get_CFL();
		if (timeStep==(restart+1)) if (Rank==0) cout << "step" << "\t" << "time" << "\t\t" << "dt" << "\t\t" << "CFL" << "\t\t" << "nIter" << "\t" << "rNorm" << endl;
		if (Rank==0) cout << timeStep << "\t" << setprecision(4) << scientific << time << "\t" << dt << "\t" << CFL << "\t" << nIter << "\t" << rNorm << endl;
		
		if ((timeStep) % restartFreq == 0) write_restart(timeStep,time);

		
		if ((timeStep) % outFreq == 0) write_output(timeStep,time,input);
		
		if ((timeStep) % probeFreq == 0 && Rank==0) {
			string fileName;
			for (int p=0;p<probes.size();++p) {
				unsigned int c=probes[p].nearestCell;
				fileName="probe"+int2str(probes[p].id+1)+".dat";
				ofstream file;
				file.open((fileName).c_str(),ios::app);
				file << timeStep << setw(16) << setprecision(8)  << "\t" << time << "\t" << grid.cell[c].rho << "\t" << grid.cell[c].v.comp[0] << "\t" << grid.cell[c].v.comp[1] << "\t" << grid.cell[c].v.comp[2] << "\t" << grid.cell[c].p << endl;
				file.close();
			}
		}

		if ((timeStep) % loadFreq == 0 && Rank==0) {
			string fileName;
			for (int n=0;n<loadCount;++n) {
				int load_bc=input.section["loads"].numberedSubsections["load"].ints[n]["bc"];
				string fileName;
				fileName="loads_bc_"+int2str(load_bc)+".dat";
				ofstream file;
				file.open((fileName).c_str(),ios::app);
				file << timeStep << setw(16) << setprecision(8)  << "\t" << time << "\t" << bc.region[load_bc-1].momentum.comp[0];
				file << "\t" << bc.region[load_bc-1].momentum.comp[1];
				file << "\t" << bc.region[load_bc-1].momentum.comp[2] << endl;
				file.close();
			}
		}

		// Flush boundary forces
		for (int i=0;i<input.section["boundaryConditions"].numberedSubsections["BC"].count;++i) {
			//cout << "BC_" << i+1 << " momentum flux:  " << bc.region[i].momentum << endl;
			//cout << i+1 << "\t" << bc.region[i].area << endl;
			bc.region[i].momentum=0.;
		}
}// DEBUG


	}

	/*****************************************************************************************/
	// End time loop
	/*****************************************************************************************/
	
	//writeTecplotMacro(restart,timeStepMax, outFreq);
	// Syncronize the processors
	MPI_Barrier(MPI_COMM_WORLD);

	// Report the wall time
	if (Rank==0) {
		timeEnd=MPI_Wtime();
		cout << "* Wall time: " << timeEnd-timeRef << " seconds" << endl;
	}
		
	return 0;
}

double minmod(double a, double b) {
	if ((a*b)<0) {
		return 0.;
	} else if (fabs(a) < fabs(b)) {
		return a;
	} else {
		return b;
	}
}

double maxmod(double a, double b) {
	if ((a*b)<0) {
		return 0.;
	} else if (fabs(a) < fabs(b)) {
		return b;
	} else {
		return a;
	}
}

double superbee(double a, double b) {
	return minmod(maxmod(a,b),minmod(a,b));
}

bool within_box(Vec3D centroid, Vec3D box_1, Vec3D box_2) {
	double tocorners_x=fabs(centroid.comp[0]-box_1.comp[0]) +fabs(centroid.comp[0]-box_2.comp[0]);
	double tocorners_y=fabs(centroid.comp[1]-box_1.comp[1]) +fabs(centroid.comp[1]-box_2.comp[1]);
	double tocorners_z=fabs(centroid.comp[2]-box_1.comp[2]) +fabs(centroid.comp[2]-box_2.comp[2]);
	if (tocorners_x<=fabs(box_1.comp[0]-box_2.comp[0])) {
		if (tocorners_y<=fabs(box_1.comp[1]-box_2.comp[1])) {
			if (tocorners_z<=fabs(box_1.comp[2]-box_2.comp[2])) {
				// The cell centroid is inside the box region
				return true;
			}
		}
	}
	return false;
}

void update(double dt) {

	double conservative[5];
	for (unsigned int c = 0;c < grid.cellCount;++c) {
		conservative[0] = grid.cell[c].rho;
		conservative[1] = grid.cell[c].rho * grid.cell[c].v.comp[0];
		conservative[2] = grid.cell[c].rho * grid.cell[c].v.comp[1];
		conservative[3] = grid.cell[c].rho * grid.cell[c].v.comp[2];
		conservative[4] = 0.5 * grid.cell[c].rho * grid.cell[c].v.dot(grid.cell[c].v) + (grid.cell[c].p+Pref) / (Gamma - 1.);
		for (int i = 0;i < 5;++i) {
			conservative[i] -= dt / grid.cell[c].volume * grid.cell[c].flux[i];
			grid.cell[c].flux[i] = 0.;
		}
		grid.cell[c].rho = conservative[0];
		grid.cell[c].v.comp[0] = conservative[1] / conservative[0];
		grid.cell[c].v.comp[1] = conservative[2] / conservative[0];
		grid.cell[c].v.comp[2] = conservative[3] / conservative[0];
		grid.cell[c].p = (conservative[4] - 0.5 * conservative[0] * grid.cell[c].v.dot(grid.cell[c].v)) * (Gamma - 1.)-Pref;
	} // cell loop
	return;
}


void updatePrimitive(double dt) {

	for (unsigned int c = 0;c < grid.cellCount;++c) {
		grid.cell[c].rho +=grid.cell[c].flux[0];
		grid.cell[c].v.comp[0] +=grid.cell[c].flux[1];
		grid.cell[c].v.comp[1] +=grid.cell[c].flux[2];
		grid.cell[c].v.comp[2] +=grid.cell[c].flux[3];
		grid.cell[c].p += grid.cell[c].flux[4];
 		grid.cell[c].k += grid.cell[c].flux[5];
 		grid.cell[c].omega += grid.cell[c].flux[6];
		for (int i = 0;i < 7;++i) {
			grid.cell[c].flux[i] = 0.;
		}
	} // cell loop
	return;
}

string int2str(int number) {
	char dummy[12];
	// Print integer to character
	sprintf(dummy, "%12d", number);
	// Convert character to string and erase leading whitespaces
	string name = dummy;
	name.erase(0, name.rfind(" ", name.length()) + 1);
	return name;
}

void get_dt(string type) {

	if (type=="fixed") return;
	
	if (type=="CFL" | type=="CFLramp") {
		// Determine time step with CFL condition
		double lengthScale;
		dt=1.E20;
		for (unsigned int c=0;c<grid.cellCount;++c) {
			double a=sqrt(Gamma*(grid.cell[c].p+Pref)/grid.cell[c].rho);
			lengthScale=grid.cell[c].lengthScale;
			dt=min(dt,CFL*lengthScale/(fabs(grid.cell[c].v.comp[0])+a));
			dt=min(dt,CFL*lengthScale/(fabs(grid.cell[c].v.comp[1])+a));
			dt=min(dt,CFL*lengthScale/(fabs(grid.cell[c].v.comp[2])+a));
		}
		MPI_Allreduce(&dt,&dt,1, MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
	}

	if (type=="CFLramp") {
		CFL=min(CFL*1.1,CFLtarget);
	}
	
	return;
} // end get_dt

void get_CFL(void) {

	// Determine time step with CFL condition
	double lengthScale;
	CFL=0.;
	for (unsigned int c=0;c<grid.cellCount;++c) {
		double a=sqrt(Gamma*(grid.cell[c].p+Pref)/grid.cell[c].rho);
		lengthScale=grid.cell[c].lengthScale;
		CFL=max(CFL,(fabs(grid.cell[c].v.comp[0])+a)*dt/lengthScale);
		CFL=max(CFL,(fabs(grid.cell[c].v.comp[1])+a)*dt/lengthScale);
		CFL=max(CFL,(fabs(grid.cell[c].v.comp[2])+a)*dt/lengthScale);;
	}
	MPI_Allreduce(&CFL,&CFL,1, MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);

	return;
}

