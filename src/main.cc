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
#define DEBUG 0

#include <iostream>
#include <fstream>
#include <map>
#include <cmath>
#include <iomanip>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <limits>

#include "commons.h"
#include "inputs.h"
#include "bc.h"
#include "mpi_functions.h"
#include "petsc_functions.h"
#include "probe.h"
using namespace std;

// Function prototypes
void read_inputs(InputFile &input);
void initialize(InputFile &input);
void write_output(double time, InputFile input);
void write_restart(double time);
void read_restart(double &time);
void setBCs(InputFile &input, BC &bc);
bool withinBox(Vec3D point,Vec3D corner_1,Vec3D corner_2);
bool withinCylinder(Vec3D point,Vec3D center,double radius,Vec3D axisDirection,double height);
bool withinSphere(Vec3D point,Vec3D center,double radius);
double minmod(double a, double b);
double maxmod(double a, double b);
double superbee(double a, double b);
void update(void);
void updatePrimitive(void);
string int2str(int number) ;
void assemble_linear_system(void);
void initialize_linear_system(void);
void get_dt(void);
void get_CFLmax(void);
void set_probes(void);
void set_loads(void);

// Allocate container for boundary conditions
BC bc;
InputFile input;

vector<Probe> probes;
vector<Load> loads;

bool grad_test=false; // DEBUG

int main(int argc, char *argv[]) {

	// Set global parameters
	sqrt_machine_error=sqrt(std::numeric_limits<double>::epsilon());
	
	// Initialize mpi
	mpi_init(argc,argv);

	string inputFileName, gridFileName;
	inputFileName.assign(argv[1]);
	gridFileName.assign(argv[1]);
	inputFileName+=".in";
	gridFileName+=".cgns";
	//gridFileName+=".dat";

	restart=0;
	if (argc>2) restart=atoi(argv[2]);

	input.setFile(inputFileName);
	read_inputs(input);
	// Physical time
	double time=0.;

	if (restart!=0 && ramp) {
		ramp=false;
	}
	int timeStepMax=input.section("timeMarching").get_int("numberOfSteps");
	
	grid.read(gridFileName);
	
	// Hand shake with other processors
	// (agree on which cells' data to be shared)
	mpi_handshake();
	
	initialize(input);
	if (Rank==0) cout << "[I] Applied initial conditions" << endl;

	setBCs(input,bc);
	if (Rank==0) cout << "[I] Set boundary conditions" << endl;

	set_probes();
	set_loads();

	grid.lengthScales();
	grid.nodeAverages();       // Linear triangular (tetrahedral) + idw blended mode
	//grid.nodeAverages_idw(); // Inverse distance based mode
	grid.faceAverages();
	grid.gradMaps();

	if (restart!=0) {
		for (int p=0;p<np;++p) {
			if (Rank==p) read_restart(time);
			MPI_Barrier(MPI_COMM_WORLD);
		}
	}

	// Initialize petsc
	petsc_init(argc,argv,
				input.section("linearSolver").get_double("relTolerance"),
				input.section("linearSolver").get_double("absTolerance"),
				input.section("linearSolver").get_int("maxIterations"));
	
	// Initialize time step or CFL size if ramping
	if (ramp) {
		if (TIME_STEP_TYPE==FIXED) dt=ramp_initial;
		if (TIME_STEP_TYPE==CFL_MAX) CFLmax=ramp_initial;
		if (TIME_STEP_TYPE==CFL_LOCAL) CFLlocal=ramp_initial;
	}

	if (Rank==0) cout << "[I] Beginning time loop" << endl;
	// Begin timing
	MPI_Barrier(MPI_COMM_WORLD);
	double timeRef, timeEnd;
	timeRef=MPI_Wtime();
	/*****************************************************************************************/
	// Begin time loop
	/*****************************************************************************************/
	
	for (timeStep=restart+1;timeStep<=timeStepMax+restart;++timeStep) {

		int nIter; // Number of linear solver iterations
		double rNorm; // Residual norm of the linear solver
if (grad_test) { // DEBUG
	grid.gradients(); // DEBUG
	write_restart(time); // DEBUG
	break;
} else {
		// Update the primitive variables of the ghost cells
		if (DEBUG) cout << "before mpi_update_ghost_primitives" << endl;
		mpi_update_ghost_primitives();
		
		// Get time step (will handle fixed, CFL and CFLramp)
		if (DEBUG) cout << "before get_dt" << endl;
		get_dt();
		if (TIME_STEP_TYPE==FIXED) { get_CFLmax(); } 
		else if (TIME_STEP_TYPE==CFL_LOCAL) { CFLmax=CFLlocal;}
		
		// Gradient are calculated for NS and/or second order schemes
		if (EQUATIONS==NS |	order==2 ) {
			// Calculate all the cell gradients for each variable
			if (DEBUG) cout << "before grid.gradients" << endl;
			grid.gradients();
			// Update gradients of the ghost cells
			if (DEBUG) cout << "before mpi_update_ghost_gradients" << endl;
			mpi_update_ghost_gradients();
		}
		
		// Limit gradients (limited gradients are stored separately)
		if (order==2) {
			if (DEBUG) cout << "before grid.limit_gradients" << endl;
			grid.limit_gradients();
			// Update limited gradients of the ghost cells
			if (DEBUG) cout << "before mpi_update_ghost_limited_gradients" << endl;
			mpi_update_ghost_limited_gradients();
		}

		if (DEBUG) cout << "before initialize_linear_system" << endl;
		initialize_linear_system();
		if (DEBUG) cout << "before assemble_linear_system" << endl;
		assemble_linear_system();
		if (DEBUG) cout << "before petsc_solve" << endl;
		petsc_solve(nIter,rNorm);
		if (DEBUG) cout << "before update_primitive" << endl;
		updatePrimitive();

		// Advance physical time
		time += dt;
		if (Rank==0) {
			if (timeStep==(restart+1))  cout << "step" << "\t" << "time" << "\t\t" << "dt" << "\t\t" << "CFLmax" << "\t\t" << "nIter" << "\t" << "rNorm" << endl;
			cout << timeStep << "\t" << setprecision(4) << scientific << time << "\t" << dt << "\t" << CFLmax << "\t" << nIter << "\t" << rNorm << endl;
		}
		
		// Ramp-up if needed
		if (ramp) {
				if (TIME_STEP_TYPE==FIXED) dt=min(dt*ramp_growth,dtTarget);
				if (TIME_STEP_TYPE==CFL_MAX) CFLmax=min(CFLmax*ramp_growth,CFLmaxTarget);
				if (TIME_STEP_TYPE==CFL_LOCAL) CFLlocal=min(CFLlocal*ramp_growth,CFLlocalTarget);
		}
			
		if ((timeStep) % restartFreq == 0) {
			for (int p=0;p<np;p++) {
				if (p==Rank) write_restart(time);
				// Syncronize the processors
				MPI_Barrier(MPI_COMM_WORLD);
			}
		}

		if ((timeStep) % outFreq == 0) write_output(time,input);
		
		if ((timeStep) % probeFreq == 0 && Rank==0) {
			for (int p=0;p<probes.size();++p) {
				Cell &c=grid.cell[probes[p].nearestCell];
				ofstream file;
				file.open((probes[p].fileName).c_str(),ios::app);
				file << timeStep << setw(16) << setprecision(8)  << "\t" << time << "\t" << c.rho << "\t" << c.v[0] << "\t" << c.v[1] << "\t" << c.v[2] << "\t" << c.p << endl;
				file.close();
			}
		}

		if ((timeStep) % loadFreq == 0 && Rank==0) {
			for (int n=0;n<loads.size();++n) {
				ofstream file;
				file.open((loads[n].fileName).c_str(),ios::app);
				file << timeStep << setw(16) << setprecision(8)  << "\t" << time << "\t" << bc.region[loads[n].bc-1].momentum[0];
				file << "\t" << bc.region[loads[n].bc-1].momentum[1];
				file << "\t" << bc.region[loads[n].bc-1].momentum[2] << endl;
				file.close();
			}
		}
		
		// Flush boundary forces
		for (int i=0;i<bcCount;++i) bc.region[i].momentum=0.;
		
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


void update(void) {

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


void updatePrimitive(void) {

	for (unsigned int c = 0;c < grid.cellCount;++c) {
		grid.cell[c].rho +=grid.cell[c].update[0];
		grid.cell[c].v.comp[0] +=grid.cell[c].update[1];
		grid.cell[c].v.comp[1] +=grid.cell[c].update[2];
		grid.cell[c].v.comp[2] +=grid.cell[c].update[3];
		grid.cell[c].p += grid.cell[c].update[4];
		grid.cell[c].k += grid.cell[c].update[5];
		grid.cell[c].omega += grid.cell[c].update[6];
	} // cell loop
	return;
}

string int2str(int number) {
	std::stringstream ss;
	ss << number;
	return ss.str();
}

void get_dt(void) {
	
	if (TIME_STEP_TYPE==CFL_MAX) {
		// Determine time step with CFL condition
		dt=1.E20;
		for (cit=grid.cell.begin();cit<grid.cell.end();cit++) {
			double a=sqrt(Gamma*((*cit).p+Pref)/(*cit).rho);
			for (int i=0;i<3;++i) dt=min(dt,CFLmax*(*cit).lengthScale/(fabs((*cit).v[i])+a));
		}
		MPI_Allreduce(&dt,&dt,1, MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
	} 

	return;
} // end get_dt

void get_CFLmax(void) {

	CFLmax=0.;
	for (cit=grid.cell.begin();cit<grid.cell.end();cit++) {
		double a=sqrt(Gamma*((*cit).p+Pref)/(*cit).rho);
		//cout << Gamma <<  "\t" << (*cit).p << "\t" << Pref << "\t" << (*cit).rho << endl;
		for (int i=0;i<3;++i) CFLmax=max(CFLmax,(fabs((*cit).v[0])+a)*dt/(*cit).lengthScale);
	}
	MPI_Allreduce(&CFLmax,&CFLmax,1, MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);

	return;
}

bool withinBox(Vec3D point,Vec3D corner_1,Vec3D corner_2) {
	Vec3D toCorners;
	bool check=true;
	for (int i=0;i<3;++i) {
		toCorners[i]=fabs(point[i]-corner_1[i])+fabs(point[i]-corner_2[i]);
		if (toCorners[i]>fabs(corner_1[i]-corner_2[i])) check=false;
	}
	return check;
}

bool withinCylinder(Vec3D point,Vec3D center,double radius,Vec3D axisDirection,double height) {
	bool check=true;
	Vec3D radialPoint=point-center;
	radialPoint=radialPoint.dot(axisDirection);
	if (fabs(axisDirection.dot(point-center))>0.5*height) check=false;
	if (fabs(radialPoint)>radius) check=false;
	return check;
}

bool withinSphere(Vec3D point,Vec3D center,double radius) {
	return (fabs(point-center)<=radius);
}
